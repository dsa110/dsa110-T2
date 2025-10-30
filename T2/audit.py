#!/usr/bin/env python3
"""
Audit FRB injections against the T2 pipeline.

We track each injected FRB through the pipeline gates and keep:
- an append-only timeline CSV (audit_log_<date>.csv)
- a rolling snapshot CSV (injections.csv)
- optional per-injection JSON state on disk (state/<inj_id>.json) if persist_json=True

Gates (-1 unknown, 0 failed, 1 passed):

  G0_injected        we injected this FRB into the data
  G1_parsed          T1 parse produced any rows for this gulp
  G2_T1_detected     we found a T1 candidate matching this injection's time/dm/beam windows
  G3_beam_kept       that same candidate survives beam-flagging (stamped in finalize)
  G4_clustered       clustering produced at least one peak consistent with this injection
  G5_filters_passed  a peak consistent with this injection survives filter_clustered
  G6_cooldown_ok     cooldown gate says we're allowed to trigger
  G7_triggered       trigger actually fired
"""

from __future__ import annotations

import os
import io
import json
import fcntl
import datetime
import logging
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, List, Tuple

import numpy as np
import pandas as pd
from astropy.time import Time
from astropy.table import Table
from astropy.io import ascii
from astropy.io.ascii.core import InconsistentTableError

from T2 import cluster_heimdall

# default output root for production runs (socket.py can override with absolute path)
DEFAULT_AUDIT_DIR = Path("audit_results")

# Append-only timeline columns, written to audit_log_<day>.csv
_AUDIT_LOG_HEADER = [
    "audit_ts_iso",
    "inj_id",
    "inj_mjd",
    "inj_beam",
    "inj_dm",
    "inj_snr_req",
    "inj_width",
    "inj_spec_ind",
    "inj_frbno",
    "gulp",
    "host",
    "candname",
    "G0_injected",
    "G1_parsed",
    "G2_T1_detected",
    "G3_beam_kept",
    "G4_clustered",
    "G5_filters_passed",
    "G6_cooldown_ok",
    "G7_triggered",
    "t1_best_snr",
    "t1_best_dm",
    "t1_best_ibox",
    "t1_best_ibeam",
    "t1_dt_sec",
    "t1_dbeam",
    "t1_ddm",
    "npoints_tab",
    "npoints_after_flag",
    "nclusters",
    "cluster_size",
    "cntb",
    "cntc",
    "peak_snr",
    "peak_dm",
    "peak_ibox",
    "peak_ibeam",
    "nbeams_this_gulp",
    "nbeams_queue_snapshot",
    "nbeams_queue_sum",
    "min_snr",
    "min_snr_wide",
    "min_snr_1arm",
    "wide_ibox",
    "max_ibox",
    "max_cntb",
    "max_cntb0",
    "max_ncl",
    "min_dm",
    "max_nbeams",
    "cooldown_wait_s",
    "t2_output_file",
    "drop_reason",
    "exception",
]

# Rolling snapshot columns, rewritten every update
_INJECTIONS_HEADER = [
    "inj_id",
    "candname",
    "inj_mjd",
    "inj_beam",
    "inj_dm",
    "inj_snr_req",
    "inj_width",
    "inj_spec_ind",
    "inj_frbno",
    "G0_injected",
    "G1_parsed",
    "G2_T1_detected",
    "G3_beam_kept",
    "G4_clustered",
    "G5_filters_passed",
    "G6_cooldown_ok",
    "G7_triggered",
    "drop_reason",
    "t1_best_snr",
    "t1_best_dm",
    "t1_best_ibox",
    "t1_best_ibeam",
    "t1_dt_sec",
    "t1_dbeam",
    "t1_ddm",
    "peak_snr",
    "peak_dm",
    "peak_ibox",
    "peak_ibeam",
    "cooldown_wait_s",
    "t2_output_file",
    "exception",
]


def _merge_gate(st: Dict[str, Any], gate_name: str, new_val: int):
    """
    Gate update: -1/0/1, keep the max.
    """
    gates = st.setdefault("gates", {})
    old_val = int(gates.get(gate_name, -1))
    gates[gate_name] = int(max(old_val, int(new_val)))


def _match_rows_for_injection(
    tab: Table,
    inj_mjd: float,
    inj_dm: float,
    inj_beam: int,
    time_window_s: float,
    dm_window: float,
    beam_window: int,
) -> Tuple[Optional[int], Dict[str, Any]]:
    """
    Find the best T1 row consistent with (time, DM, beam) windows.
    Ranking: higher SNR, then smaller |dt|, |ddm|, |dbeam|.
    Returns (row_index_or_None, summary_dict).
    """
    if tab is None or len(tab) == 0:
        return None, {}

    n = len(tab)

    # time offset
    if "mjds" in tab.colnames:
        dt_sec = np.abs(tab["mjds"] - inj_mjd) * 86400.0
    else:
        dt_sec = np.full(n, np.inf)

    # DM offset
    if "dm" in tab.colnames:
        ddm = np.abs(tab["dm"] - inj_dm)
    else:
        ddm = np.full(n, np.inf)

    if "ibeam" in tab.colnames:
        dbeam0 = np.abs(tab["ibeam"] - inj_beam)
        dbeam1 = np.abs(tab["ibeam"] - (inj_beam + 256))
        dbeam_min = np.minimum(dbeam0, dbeam1)
        which_arm = np.where(dbeam0 <= dbeam1, 0, 1)
    else:
        dbeam_min = np.full(n, np.inf)
        which_arm = np.zeros(n, dtype=int)

    time_mask = dt_sec <= float(time_window_s)
    dm_mask = ddm <= float(dm_window)
    beam_mask = dbeam_min <= int(beam_window)

    win_mask = time_mask & dm_mask & beam_mask
    if not np.any(win_mask):
        return None, {}

    cand_ix = np.where(win_mask)[0]

    if "snr" in tab.colnames:
        snr_arr = np.array([float(tab["snr"][i]) for i in cand_ix])
    else:
        snr_arr = np.full(len(cand_ix), -np.inf)

    order = np.lexsort(
        (
            dbeam_min[cand_ix],
            ddm[cand_ix],
            dt_sec[cand_ix],
            -snr_arr,
        )
    )
    best_idx = int(cand_ix[order[0]])

    # report arm-corrected beam distance
    arm = which_arm[best_idx]
    if arm == 0:
        dbeam_report = int(np.abs(tab["ibeam"][best_idx] - inj_beam))
    else:
        dbeam_report = int(np.abs(tab["ibeam"][best_idx] - (inj_beam + 256)))

    out = dict(
        t1_best_snr=float(tab["snr"][best_idx]) if "snr" in tab.colnames else "",
        t1_best_dm=float(tab["dm"][best_idx]) if "dm" in tab.colnames else "",
        t1_best_ibox=int(tab["ibox"][best_idx]) if "ibox" in tab.colnames else "",
        t1_best_ibeam=int(tab["ibeam"][best_idx]) if "ibeam" in tab.colnames else "",
        t1_dt_sec=float(dt_sec[best_idx]) if np.isfinite(dt_sec[best_idx]) else "",
        t1_dbeam=dbeam_report,
        t1_ddm=float(ddm[best_idx]) if np.isfinite(ddm[best_idx]) else "",
    )

    return best_idx, out


def _extract_peak_info_for_injection(
    tab_peak: Optional[Table],
    inj_mjd: float,
    inj_dm: float,
    inj_beam: int,
    time_window_s: float,
    dm_window: float,
    beam_window: int,
) -> Dict[str, Any]:
    """Pick the cluster peak matching this injection; same windows as T1."""

    best_idx, _ = _match_rows_for_injection(
        tab_peak,
        inj_mjd,
        inj_dm,
        inj_beam,
        time_window_s,
        dm_window,
        beam_window,
    )

    if best_idx is None or tab_peak is None or len(tab_peak) == 0:
        return {}

    r = tab_peak[best_idx]
    out = {}
    out["snr"] = float(r["snr"]) if "snr" in r.colnames else ""
    out["dm"] = float(r["dm"]) if "dm" in r.colnames else ""
    out["ibox"] = int(r["ibox"]) if "ibox" in r.colnames else ""
    out["ibeam"] = int(r["ibeam"]) if "ibeam" in r.colnames else ""
    out["cntb"] = int(r["cntb"]) if "cntb" in r.colnames else ""
    out["cntc"] = int(r["cntc"]) if "cntc" in r.colnames else ""
    return out


def _earliest_drop_reason(gates: Dict[str, int]) -> Optional[str]:
    """
    Get reason for early gate failure, if any.
    """
    g1 = gates.get("G1_parsed", -1)
    g2 = gates.get("G2_T1_detected", -1)

    if g1 == 0:
        return "no_parse"
    if g2 == 0:
        return "no_T1_in_window"
    return None


def _stable_inj_id(mjd: float, beam: int, dm: float) -> str:
    """Create a stable, unique id for an injection. Had to do this because FRBno was not unique."""
    return f"{float(mjd):.9f}_B{int(beam)}_DM{float(dm):.2f}"


def _read_legacy_injection_table(path: str) -> Table:
    """Read legacy injection file."""
    names = "MJD   Beam   DM    SNR   Width_fwhm   spec_ind  FRBno".split()
    try:
        tab = ascii.read(path, guess=True, fast_reader=False)
        # If the file has no column names, attach the expected ones:
        if any(str(c).startswith("col") for c in tab.colnames) or len(tab.colnames) != 7:
            tab = ascii.read(path, names=names, guess=False, fast_reader=False)
    except InconsistentTableError:
        tab = ascii.read(path, names=names, guess=False, fast_reader=False)
    return tab


class Auditor:
    """
    Make one of these in tests or pipeline.

    self.injections[inj_id] holds:
      {
        "inj": {
          "MJD": ...,
          "Beam": ...,
          "DM": ...,
          "SNR": ...,
          "Width_fwhm": ...,
          "spec_ind": ...,
          "FRBno": ...
        },
        "gates": {...},
        "candname": "...",
        "drop_reason": "...",
        "t1_best": {...},
        "stats": {...}
      }
    """

    def __init__(
        self,
        audit_dir: Optional[str] = None,
        time_window_s: int = 300, 
        dm_window: float = 20.0,
        beam_window: int = 2,
        persist_json: bool = False,
    ):
        self.audit_dir = Path(audit_dir) if audit_dir else DEFAULT_AUDIT_DIR
        self.audit_dir.mkdir(parents=True, exist_ok=True)

        self.persist_json = bool(persist_json)

        # only create state dir if we're actually persisting JSON
        self.state_dir = None
        if self.persist_json:
            self.state_dir = self.audit_dir / "state"
            self.state_dir.mkdir(parents=True, exist_ok=True)

        self.time_window_s = int(time_window_s)
        self.dm_window = float(dm_window)
        self.beam_window = int(beam_window)

        # inj_id => dict state
        self.injections: Dict[str, Dict[str, Any]] = {}

        # logger setup
        self.log = logging.getLogger("AUDIT")
        if not self.log.handlers:
            h = logging.StreamHandler()
            fmt = "[%(asctime)s][AUDIT][%(levelname)s] %(message)s"
            dfmt = "%Y-%m-%d %H:%M:%S"
            h.setFormatter(logging.Formatter(fmt=fmt, datefmt=dfmt))
            self.log.addHandler(h)
            self.log.setLevel(logging.INFO)

        self.log.info(
            f"[AUDIT] init audit_dir={self.audit_dir} "
            f"persist_json={self.persist_json} "
            f"time_window_s={self.time_window_s} "
            f"dm_window={self.dm_window} "
            f"beam_window={self.beam_window}"
        )
        print(
            f"[AUDIT] init audit_dir={self.audit_dir} "
            f"persist_json={self.persist_json} "
            f"time_window_s={self.time_window_s} "
            f"dm_window={self.dm_window} "
            f"beam_window={self.beam_window}"
        )

    def ingest_legacy_injections(self, legacy_path: str) -> int:
        """
        Read the old injection_list.txt and assign a unique id for each of them. FRBNo was not unique.
        This is used to create a copy downstream, so we can track the new injections from audits.
        """
        if not legacy_path or not os.path.exists(legacy_path):
            self.log.info(f"[AUDIT][ingest] legacy file missing: {legacy_path}")
            return 0

        try:
            tab = _read_legacy_injection_table(legacy_path)
        except Exception as e:
            self.log.warning(f"[AUDIT][ingest] failed reading legacy file {legacy_path}: {e}")
            return 0

        added = 0
        for row in tab:
            try:
                mjd = float(row["MJD"])
                beam = int(row["Beam"])
                dm = float(row["DM"])
                inj_id = _stable_inj_id(mjd, beam, dm)

                if inj_id in self.injections:
                    continue

                inj_meta = dict(
                    MJD=mjd,
                    Beam=beam,
                    DM=dm,
                    SNR=float(row["SNR"]) if "SNR" in row.colnames else "",
                    Width_fwhm=float(row["Width_fwhm"]) if "Width_fwhm" in row.colnames else "",
                    spec_ind=float(row["spec_ind"]) if "spec_ind" in row.colnames else "",
                    FRBno=str(row["FRBno"]) if "FRBno" in row.colnames else "",
                )
                self.seed_injection(inj_id, inj_meta)
                added += 1
            except Exception as e:
                self.log.warning(f"[AUDIT][ingest] skipping row due to parse error: {e}")
                continue

        self.log.info(f"[AUDIT][ingest] added {added} new injections from legacy file")
        print(f"[AUDIT][ingest] added {added} new injections from legacy file")
        return added

    def _now_iso(self) -> str:
        return datetime.datetime.utcnow().isoformat()

    def _audit_log_path(self) -> Path:
        day = datetime.datetime.utcnow().strftime("%Y%m%d")
        return self.audit_dir / f"audit_log_{day}.csv"

    def _ensure_audit_header(self):
        """
        Make sure audit_log_<date>.csv exists and has header.
        We keep the leading two '#' comment lines for documentation,
        then write the CSV header row.
        """
        p = self._audit_log_path()
        if not p.exists():
            with open(p, "w") as fh:
                fh.write("# Audit timeline for injected FRBs vs pipeline gates.\n")
                fh.write("# update_from_tab stamps G1,G2; finalize stamps G3..G7.\n")
                fh.write(",".join(_AUDIT_LOG_HEADER) + "\n")

    def _injections_csv_path(self) -> Path:
        return self.audit_dir / "injections.csv"

    def _state_json_path(self, inj_id: str) -> Path:
        if self.state_dir is None:
            raise RuntimeError("state dir requested but persist_json=False")
        return self.state_dir / f"{inj_id}.json"

    def _flock_append_csvrow(self, path: Path, row_dict: Dict[str, Any]):
        """
        Append one CSV row (dict with keys matching _AUDIT_LOG_HEADER)
        using flock to avoid writer collisions.
        """
        df = pd.DataFrame([row_dict], columns=_AUDIT_LOG_HEADER)
        csv_buf = io.StringIO()
        df.to_csv(csv_buf, index=False, header=False)
        line = csv_buf.getvalue()

        with open(path, "a") as fh:
            fcntl.flock(fh, fcntl.LOCK_EX)
            fh.write(line if line.endswith("\n") else line + "\n")
            fcntl.flock(fh, fcntl.LOCK_UN)

    def _median_mjd_from_tab(self, tab: Table) -> float:
        if tab is not None and len(tab) and ("mjds" in tab.colnames):
            try:
                return float(np.median(tab["mjds"]))
            except Exception:
                pass
        # fallback: current time
        return Time.now().mjd

    def _iter_injections_in_window(
        self,
        mjd_center: float,
    ) -> Iterable[Tuple[str, Dict[str, Any]]]:
        """
        Yield all injections whose MJD is within last 2 days
        """
        tw_days = self.time_window_s / 86400.0
        mjd_min = mjd_center - tw_days
        mjd_max = mjd_center + tw_days
        for inj_id, st in self.injections.items():
            inj_mjd = float(st["inj"]["MJD"])
            if mjd_min <= inj_mjd <= mjd_max:
                yield inj_id, st

    def _persist_state_json(self, inj_id: str):
        """
        Write out current state for an injection, only if persist_json=True.
        """
        if not self.persist_json:
            return
        try:
            pj = self._state_json_path(inj_id)
            with open(pj, "w") as fh:
                json.dump(self.injections[inj_id], fh, indent=2)
            msg = f"[AUDIT] persisted state JSON for inj_id={inj_id} -> {pj}"
            self.log.debug(msg)
            print(msg)
        except Exception as e:
            self.log.warning(
                f"[AUDIT] failed to persist state JSON for inj_id={inj_id}: {e}"
            )
            print(f"[AUDIT] failed to persist state JSON for inj_id={inj_id}: {e}")

    # public API

    def seed_injection(self, inj_id: str, inj_meta: Dict[str, Any]):
        """
        Called by the injection code or by legacy ingestion.
        Create initial record with gates all -1 except G0 (injected=1).
        """
        st = dict(
            inj=dict(
                MJD=float(inj_meta["MJD"]),
                Beam=int(inj_meta["Beam"]),
                DM=float(inj_meta["DM"]),
                SNR=float(inj_meta.get("SNR", "")) if str(inj_meta.get("SNR", "")) != "" else "",
                Width_fwhm=float(inj_meta.get("Width_fwhm", "")) if str(inj_meta.get("Width_fwhm", "")) != "" else "",
                spec_ind=float(inj_meta.get("spec_ind", "")) if str(inj_meta.get("spec_ind", "")) != "" else "",
                FRBno=str(inj_meta.get("FRBno", "")),
            ),
            gates={
                "G0_injected": 1,
                "G1_parsed": -1,
                "G2_T1_detected": -1,
                "G3_beam_kept": -1,
                "G4_clustered": -1,
                "G5_filters_passed": -1,
                "G6_cooldown_ok": -1,
                "G7_triggered": -1,
            },
            candname="",
            drop_reason="",
            t1_best={},
            stats={},
        )

        self.injections[inj_id] = st

        if self.persist_json:
            self._persist_state_json(inj_id)

        # log snapshot
        self._append_audit_log(inj_id, st, meta={}, exception="")

        # refresh rolling CSV
        self._write_injections_csv()

        self.log.info(f"[AUDIT][seed] new inj_id={inj_id}")
        print(f"[AUDIT][seed] new inj_id={inj_id}")

    def update_from_tab(
        self,
        host: str,
        gulp: int,
        tab: Table,
        nbeams_queue_snapshot: Iterable[int],
        nbeams_queue_sum: int,
        runtime_thresholds: Dict[str, Any],
        prev_trig_time: Optional[Time],
    ):
        """
        Called after parse_candsfile. We only update early gates:
        G1_parsed        parser produced any rows for this gulp
        G2_T1_detected   our injection shows up in that T1 output
        """
        try:
            mjd_now = self._median_mjd_from_tab(tab)
            self.log.info(
                f"[AUDIT][update] start gulp={gulp} host={host} len(tab)={len(tab)} mjd_now={mjd_now}"
            )
            print(
                f"[AUDIT][update] start gulp={gulp} host={host} len(tab)={len(tab)} mjd_now={mjd_now}"
            )

            for inj_id, st in self._iter_injections_in_window(mjd_now):
                inj = st["inj"]
                inj_mjd = float(inj["MJD"])
                inj_dm = float(inj["DM"])
                inj_beam = int(inj["Beam"])

                npoints_tab = int(len(tab))
                gate_G1_parsed = 1 if npoints_tab > 0 else 0

                # match injection to tab rows
                match_idx, t1_best = _match_rows_for_injection(
                    tab,
                    inj_mjd=inj_mjd,
                    inj_dm=inj_dm,
                    inj_beam=inj_beam,
                    time_window_s=self.time_window_s,
                    dm_window=self.dm_window,
                    beam_window=self.beam_window,
                )
                gate_G2_T1_detected = 1 if match_idx is not None else 0

                if gate_G2_T1_detected == 1 and t1_best:
                    st["t1_best"] = dict(t1_best)

                _merge_gate(st, "G1_parsed", gate_G1_parsed)
                _merge_gate(st, "G2_T1_detected", gate_G2_T1_detected)

                # earliest failure reason, if any, after this stage
                er = _earliest_drop_reason(st["gates"])
                if er:
                    st["drop_reason"] = er

                st["stats"]["npoints_tab"] = npoints_tab
                st["stats"]["npoints_after_flag"] = ""  # filled in finalize
                st["stats"]["nbeams_queue_snapshot"] = list(nbeams_queue_snapshot)
                st["stats"]["nbeams_queue_sum"] = int(nbeams_queue_sum)
                st["stats"]["runtime_thresholds"] = dict(runtime_thresholds)

                if self.persist_json:
                    self._persist_state_json(inj_id)

                meta = {
                    "gulp": gulp,
                    "host": host,
                    "npoints_tab": npoints_tab,
                    "npoints_after_flag": "",
                    "nbeams_queue_snapshot": json.dumps(list(nbeams_queue_snapshot)),
                    "nbeams_queue_sum": int(nbeams_queue_sum),
                    "drop_reason": st.get("drop_reason", ""),
                }
                meta.update(t1_best)
                self._append_audit_log(inj_id, st, meta, exception="")

            self._write_injections_csv()

            self.log.info(
                f"[AUDIT][update] done: considered_injections={len(list(self.injections.keys()))}"
            )
            print(
                f"[AUDIT][update] done: considered_injections={len(list(self.injections.keys()))}"
            )

        except Exception as e:
            self.log.warning(f"[AUDIT][update] fatal error: {e}")
            print(f"[AUDIT][update] fatal error: {e}")

    def finalize_from_cluster_result(
        self,
        host: str,
        gulp: int,
        tab_pre_filter: Table,
        tab_after_beam_flag: Optional[Table],  
        tab_peak: Optional[Table],
        tab_after_filters: Optional[Table],
        lastname: Optional[str],
        trigtime: Optional[Time],
        triggered: bool,
        nbeams_this_gulp: Optional[int],
        nbeams_queue_snapshot: Iterable[int],
        nbeams_queue_sum: int,
        thresholds: Dict[str, Any],
        cooldown_wait_s: Optional[float],
        candname_for_injection: Optional[str],
    ):
        """
        After clustering: decide G3..G7 in order.
        """
        try:
            mjd_now = self._median_mjd_from_tab(tab_pre_filter)
            self.log.info(
                f"[AUDIT][finalize] start gulp={gulp} host={host} triggered={triggered} "
                f"mjd_now={mjd_now} len(pre)={(len(tab_pre_filter) if tab_pre_filter is not None else 'NA')} "
                f"len(flag)={(len(tab_after_beam_flag) if tab_after_beam_flag is not None else 'NA')} "
                f"len(peak)={(len(tab_peak) if tab_peak is not None else '0')} "
                f"len(after)={(len(tab_after_filters) if tab_after_filters is not None else '0')}"
            )
            print(
                f"[AUDIT][finalize] start gulp={gulp} host={host} triggered={triggered} "
                f"mjd_now={mjd_now} len(pre)={(len(tab_pre_filter) if tab_pre_filter is not None else 'NA')} "
                f"len(flag)={(len(tab_after_beam_flag) if tab_after_beam_flag is not None else 'NA')} "
                f"len(peak)={(len(tab_peak) if tab_peak is not None else '0')} "
                f"len(after)={(len(tab_after_filters) if tab_after_filters is not None else '0')}"
            )

            # record npoints_after_flag stats for all considered injections
            npoints_after_flag_val = (int(len(tab_after_beam_flag)) if tab_after_beam_flag is not None else "")

            for inj_id, st in self._iter_injections_in_window(mjd_now):
                inj = st["inj"]
                inj_mjd = float(inj["MJD"])
                inj_dm = float(inj["DM"])
                inj_beam = int(inj["Beam"])

                gates = st["gates"]

                if gates.get("G2_T1_detected", -1) == 1 and tab_after_beam_flag is not None:
                    match_idx_flag, _ = _match_rows_for_injection(
                        tab_after_beam_flag,
                        inj_mjd=inj_mjd,
                        inj_dm=inj_dm,
                        inj_beam=inj_beam,
                        time_window_s=self.time_window_s,
                        dm_window=self.dm_window,
                        beam_window=self.beam_window,
                    )
                    gate_G3_beam_kept = 1 if match_idx_flag is not None else 0
                else:
                    # unknown at this point (e.g., not detected at G2 or missing flag table)
                    gate_G3_beam_kept = -1

                _merge_gate(st, "G3_beam_kept", gate_G3_beam_kept)
                st["stats"]["npoints_after_flag"] = npoints_after_flag_val

                # G1, G2 (Get Drop Reason)
                if not (gates.get("G1_parsed", -1) == 1 and gates.get("G2_T1_detected", -1) == 1):
                    er = _earliest_drop_reason(st["gates"])
                    if er:
                        st["drop_reason"] = er
                #G3 (Get Drop Reason)
                elif st["gates"].get("G3_beam_kept", -1) != 1:
                    if st["gates"].get("G3_beam_kept", -1) == 0:
                        st["drop_reason"] = "beam_flagged"

                else:
                    # summarize peak rows that match THIS injection
                    peak_info = _extract_peak_info_for_injection(
                        tab_peak,
                        inj_mjd=inj_mjd,
                        inj_dm=inj_dm,
                        inj_beam=inj_beam,
                        time_window_s=self.time_window_s,
                        dm_window=self.dm_window,
                        beam_window=self.beam_window,
                    )

                    # G4_clustered:
                    gate_G4_clustered = 1 if peak_info else 0

                    # G5_filters_passed:
                    match_idx_after, _tmp_after = _match_rows_for_injection(
                        tab_after_filters,
                        inj_mjd=inj_mjd,
                        inj_dm=inj_dm,
                        inj_beam=inj_beam,
                        time_window_s=self.time_window_s,
                        dm_window=self.dm_window,
                        beam_window=self.beam_window,
                    )
                    gate_G5_filters_passed = 1 if match_idx_after is not None else 0

                    max_nbeams_allowed = int(thresholds.get("max_nbeams", 40))

                    # G6_cooldown_ok:
                    if cooldown_wait_s is None:
                        gate_G6_cooldown_ok = 0
                    else:
                        min_timedelt = float(thresholds.get("min_timedelt", 60.0))
                        gate_G6_cooldown_ok = 1 if cooldown_wait_s >= min_timedelt else 0

                    # G7_triggered:
                    gate_G7_triggered = 1 if bool(triggered) else 0

                    # update gate states
                    _merge_gate(st, "G4_clustered", gate_G4_clustered)
                    _merge_gate(st, "G5_filters_passed", gate_G5_filters_passed)
                    _merge_gate(st, "G6_cooldown_ok", gate_G6_cooldown_ok)
                    _merge_gate(st, "G7_triggered", gate_G7_triggered)

                    # name the candidate we're associating
                    if candname_for_injection:
                        st["candname"] = candname_for_injection

                    # compute late-stage drop_reason only if early path did not fail
                    drop_reason = ""
                    if gate_G4_clustered == 0:
                        drop_reason = "no_cluster_peak"
                    elif gate_G5_filters_passed == 0:
                        drop_reason = "filtered_out"
                    elif nbeams_queue_sum > max_nbeams_allowed:
                        drop_reason = (
                            f"nbeams_gate_exceeded({nbeams_queue_sum}>{max_nbeams_allowed})"
                        )
                    elif gate_G6_cooldown_ok == 0:
                        if cooldown_wait_s is None:
                            drop_reason = "cooldown_unknown"
                        else:
                            drop_reason = (
                                f"cooldown({cooldown_wait_s:.2f}s<"
                                f"{float(thresholds.get('min_timedelt',60.0))}s)"
                            )
                    elif gate_G7_triggered == 0:
                        drop_reason = "not_triggered"
                    
                    er = _earliest_drop_reason(st["gates"])
                    st["drop_reason"] = er if er else drop_reason

                    # stats snapshot
                    st["stats"]["nclusters"] = (
                        int(len(tab_peak)) if tab_peak is not None else 0
                    )
                    st["stats"]["cluster_size"] = peak_info.get("cntc", "")
                    st["stats"]["cntb"] = peak_info.get("cntb", "")
                    st["stats"]["cntc"] = peak_info.get("cntc", "")
                    st["stats"]["peak_snr"] = peak_info.get("snr", "")
                    st["stats"]["peak_dm"] = peak_info.get("dm", "")
                    st["stats"]["peak_ibox"] = peak_info.get("ibox", "")
                    st["stats"]["peak_ibeam"] = peak_info.get("ibeam", "")
                    st["stats"]["nbeams_this_gulp"] = (
                        int(nbeams_this_gulp) if nbeams_this_gulp is not None else ""
                    )
                    st["stats"]["nbeams_queue_snapshot"] = list(nbeams_queue_snapshot)
                    st["stats"]["nbeams_queue_sum"] = int(nbeams_queue_sum)
                    st["stats"]["cooldown_wait_s"] = (
                        cooldown_wait_s if cooldown_wait_s else ""
                    )
                    st["stats"]["t2_output_file"] = thresholds.get(
                        "t2_output_file", ""
                    )
                    st["stats"]["exception"] = thresholds.get("exception", "")

                # force-sync JSON after updating gates or stats
                if self.persist_json:
                    self._persist_state_json(inj_id)

                # write timeline row for this finalize snapshot
                meta = dict(
                    gulp=gulp,
                    host=host,
                    nclusters=st["stats"].get("nclusters", ""),
                    cluster_size=st["stats"].get("cluster_size", ""),
                    cntb=st["stats"].get("cntb", ""),
                    cntc=st["stats"].get("cntc", ""),
                    peak_snr=st["stats"].get("peak_snr", ""),
                    peak_dm=st["stats"].get("peak_dm", ""),
                    peak_ibox=st["stats"].get("peak_ibox", ""),
                    peak_ibeam=st["stats"].get("peak_ibeam", ""),
                    nbeams_this_gulp=st["stats"].get("nbeams_this_gulp", ""),
                    nbeams_queue_snapshot=json.dumps(list(nbeams_queue_snapshot)),
                    nbeams_queue_sum=int(nbeams_queue_sum),
                    cooldown_wait_s=st["stats"].get("cooldown_wait_s", ""),
                    t2_output_file=st["stats"].get("t2_output_file", ""),
                    drop_reason=st.get("drop_reason", ""),
                    exception=st["stats"].get("exception", ""),
                )

                if "t1_best" in st and st["t1_best"]:
                    meta.update(st["t1_best"])

                self._append_audit_log(
                    inj_id=inj_id,
                    st=st,
                    meta=meta,
                    exception=st["stats"].get("exception", ""),
                )

            self._write_injections_csv()

            self.log.info(
                f"[AUDIT][finalize] done: considered_injections={len(list(self.injections.keys()))}"
            )
            print(
                f"[AUDIT][finalize] done: considered_injections={len(list(self.injections.keys()))}"
            )

        except Exception as e:
            self.log.warning(f"[AUDIT][finalize] fatal error: {e}")
            print(f"[AUDIT][finalize] fatal error: {e}")

    def _append_audit_log(
        self,
        inj_id: str,
        st: Dict[str, Any],
        meta: Dict[str, Any],
        exception: str,
    ):
        """
        Append a new row to audit_log_<day>.csv capturing current state.
        Uses pandas for CSV escaping and fcntl for locking.
        """

        inj = st["inj"]
        gates = st["gates"]
        telem = st.get("stats", {})
        t1_best = st.get("t1_best", {})

        row = {
            "audit_ts_iso": self._now_iso(),
            "inj_id": inj_id,
            "inj_mjd": f"{float(inj['MJD']):.9f}",
            "inj_beam": int(inj["Beam"]),
            "inj_dm": float(inj["DM"]),
            "inj_snr_req": float(inj.get("SNR", "")) if str(inj.get("SNR", "")) != "" else "",
            "inj_width": float(inj.get("Width_fwhm", "")) if str(inj.get("Width_fwhm", "")) != "" else "",
            "inj_spec_ind": float(inj.get("spec_ind", "")) if str(inj.get("spec_ind", "")) != "" else "",
            "inj_frbno": str(inj.get("FRBno", "")),
            "gulp": meta.get("gulp", ""),
            "host": meta.get("host", ""),
            "candname": st.get("candname", ""),
            "G0_injected": gates.get("G0_injected", -1),
            "G1_parsed": gates.get("G1_parsed", -1),
            "G2_T1_detected": gates.get("G2_T1_detected", -1),
            "G3_beam_kept": gates.get("G3_beam_kept", -1),
            "G4_clustered": gates.get("G4_clustered", -1),
            "G5_filters_passed": gates.get("G5_filters_passed", -1),
            "G6_cooldown_ok": gates.get("G6_cooldown_ok", -1),
            "G7_triggered": gates.get("G7_triggered", -1),
            "t1_best_snr": t1_best.get("t1_best_snr", ""),
            "t1_best_dm": t1_best.get("t1_best_dm", ""),
            "t1_best_ibox": t1_best.get("t1_best_ibox", ""),
            "t1_best_ibeam": t1_best.get("t1_best_ibeam", ""),
            "t1_dt_sec": t1_best.get("t1_dt_sec", ""),
            "t1_dbeam": t1_best.get("t1_dbeam", ""),
            "t1_ddm": t1_best.get("t1_ddm", ""),
            "npoints_tab": telem.get("npoints_tab", meta.get("npoints_tab", "")),
            "npoints_after_flag": telem.get(
                "npoints_after_flag", meta.get("npoints_after_flag", "")
            ),
            "nclusters": telem.get("nclusters", meta.get("nclusters", "")),
            "cluster_size": telem.get("cluster_size", meta.get("cluster_size", "")),
            "cntb": telem.get("cntb", meta.get("cntb", "")),
            "cntc": telem.get("cntc", meta.get("cntc", "")),
            "peak_snr": telem.get("peak_snr", meta.get("peak_snr", "")),
            "peak_dm": telem.get("peak_dm", meta.get("peak_dm", "")),
            "peak_ibox": telem.get("peak_ibox", meta.get("peak_ibox", "")),
            "peak_ibeam": telem.get("peak_ibeam", meta.get("peak_ibeam", "")),
            "nbeams_this_gulp": telem.get(
                "nbeams_this_gulp", meta.get("nbeams_this_gulp", "")
            ),
            "nbeams_queue_snapshot": meta.get(
                "nbeams_queue_snapshot",
                json.dumps(telem.get("nbeams_queue_snapshot", [])),
            ),
            "nbeams_queue_sum": meta.get(
                "nbeams_queue_sum", telem.get("nbeams_queue_sum", "")
            ),
            "min_snr": telem.get("runtime_thresholds", {}).get(
                "min_snr", meta.get("min_snr", "")
            ),
            "min_snr_wide": telem.get("runtime_thresholds", {}).get(
                "min_snr_wide", meta.get("min_snr_wide", "")
            ),
            "min_snr_1arm": telem.get("runtime_thresholds", {}).get(
                "min_snr_1arm", meta.get("min_snr_1arm", "")
            ),
            "wide_ibox": telem.get("runtime_thresholds", {}).get(
                "wide_ibox", meta.get("wide_ibox", "")
            ),
            "max_ibox": telem.get("runtime_thresholds", {}).get(
                "max_ibox", meta.get("max_ibox", "")
            ),
            "max_cntb": telem.get("runtime_thresholds", {}).get(
                "max_cntb", meta.get("max_cntb", "")
            ),
            "max_cntb0": telem.get("runtime_thresholds", {}).get(
                "max_cntb0", meta.get("max_cntb0", "")
            ),
            "max_ncl": telem.get("runtime_thresholds", {}).get(
                "max_ncl", meta.get("max_ncl", "")
            ),
            "min_dm": telem.get("runtime_thresholds", {}).get(
                "min_dm", meta.get("min_dm", "")
            ),
            "max_nbeams": telem.get("runtime_thresholds", {}).get(
                "max_nbeams", meta.get("max_nbeams", "")
            ),
            "cooldown_wait_s": telem.get(
                "cooldown_wait_s", meta.get("cooldown_wait_s", "")
            ),
            "t2_output_file": telem.get(
                "t2_output_file", meta.get("t2_output_file", "")
            ),
            "drop_reason": st.get("drop_reason", meta.get("drop_reason", "")),
            "exception": exception or telem.get("exception", meta.get("exception", "")),
        }

        self._ensure_audit_header()
        self._flock_append_csvrow(self._audit_log_path(), row)

        self.log.info(f"[AUDIT][append_audit_log] wrote audit row for inj_id={inj_id}")
        print(f"[AUDIT][append_audit_log] wrote audit row for inj_id={inj_id}")

    def _write_injections_csv(self):
        """
        Rewrite injections.csv with one line per injection reflecting latest state.

        Also, if persist_json=True, force-write every per-injection JSON snapshot.
        """
        rows: List[Dict[str, Any]] = []

        for inj_id, st in self.injections.items():
            inj = st["inj"]
            gates = st["gates"]
            telem = st.get("stats", {})
            t1_best = st.get("t1_best", {})

            row = {
                "inj_id": inj_id,
                "candname": st.get("candname", ""),
                "inj_mjd": f"{float(inj['MJD']):.9f}",
                "inj_beam": int(inj["Beam"]),
                "inj_dm": float(inj["DM"]),
                "inj_snr_req": float(inj.get("SNR", "")) if str(inj.get("SNR", "")) != "" else "",
                "inj_width": float(inj.get("Width_fwhm", "")) if str(inj.get("Width_fwhm", "")) != "" else "",
                "inj_spec_ind": float(inj.get("spec_ind", "")) if str(inj.get("spec_ind", "")) != "" else "",
                "inj_frbno": str(inj.get("FRBno", "")),
                "G0_injected": gates.get("G0_injected", -1),
                "G1_parsed": gates.get("G1_parsed", -1),
                "G2_T1_detected": gates.get("G2_T1_detected", -1),
                "G3_beam_kept": gates.get("G3_beam_kept", -1),
                "G4_clustered": gates.get("G4_clustered", -1),
                "G5_filters_passed": gates.get("G5_filters_passed", -1),
                "G6_cooldown_ok": gates.get("G6_cooldown_ok", -1),
                "G7_triggered": gates.get("G7_triggered", -1),
                "drop_reason": st.get("drop_reason", ""),
                "t1_best_snr": t1_best.get("t1_best_snr", ""),
                "t1_best_dm": t1_best.get("t1_best_dm", ""),
                "t1_best_ibox": t1_best.get("t1_best_ibox", ""),
                "t1_best_ibeam": t1_best.get("t1_best_ibeam", ""),
                "t1_dt_sec": t1_best.get("t1_dt_sec", ""),
                "t1_dbeam": t1_best.get("t1_dbeam", ""),
                "t1_ddm": t1_best.get("t1_ddm", ""),
                "peak_snr": telem.get("peak_snr", ""),
                "peak_dm": telem.get("peak_dm", ""),
                "peak_ibox": telem.get("peak_ibox", ""),
                "peak_ibeam": telem.get("peak_ibeam", ""),
                "cooldown_wait_s": telem.get("cooldown_wait_s", ""),
                "t2_output_file": telem.get("t2_output_file", ""),
                "exception": telem.get("exception", ""),
            }

            rows.append(row)

        p = self._injections_csv_path()
        df = pd.DataFrame(rows, columns=_INJECTIONS_HEADER)
        df.to_csv(p, index=False)

        self.log.info(f"[AUDIT][write_injections_csv] wrote {len(rows)} rows to {p}")
        print(f"[AUDIT][write_injections_csv] wrote {len(rows)} rows to {p}")

        if self.persist_json:
            for inj_id in self.injections.keys():
                self._persist_state_json(inj_id)
