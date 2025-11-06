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
import csv
import threading

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

def _roundish(x, nd=3):
    try:
        return round(float(x), nd)
    except Exception:
        return x

def _as_int_seconds(x):
    try:
        return int(float(x))
    except Exception:
        return ""
        
def _state_fingerprint(st: dict) -> tuple:
    gates = st.get("gates", {})
    t1b = st.get("t1_best", {})
    telem = st.get("stats", {})

    gates_tup = tuple(gates.get(k, -1) for k in [
        "G0_injected","G1_parsed","G2_T1_detected","G3_beam_kept",
        "G4_clustered","G5_filters_passed","G6_cooldown_ok","G7_triggered"
    ])

    t1_tup = (
        _roundish(t1b.get("t1_best_snr","")),
        _roundish(t1b.get("t1_best_dm","")),
        t1b.get("t1_best_ibox",""),
        t1b.get("t1_best_ibeam",""),
        _roundish(t1b.get("t1_dt_sec","")),
        _roundish(t1b.get("t1_ddm","")),
        t1b.get("t1_dbeam",""),
    )

    telem_tup = (
        telem.get("npoints_tab",""),
        telem.get("npoints_after_flag",""),
        telem.get("nclusters",""),
        telem.get("cluster_size",""),
        telem.get("cntb",""),
        telem.get("cntc",""),
        _roundish(telem.get("peak_snr","")),
        _roundish(telem.get("peak_dm","")),
        telem.get("peak_ibox",""),
        telem.get("peak_ibeam",""),
        _as_int_seconds(telem.get("cooldown_wait_s","")),
        telem.get("t2_output_file",""),
    )

    return (gates_tup, t1_tup, telem_tup, st.get("drop_reason",""), st.get("candname",""))


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
        base_local = int(inj_beam) % 256      # normalize to 0..255
        ew_global  = base_local               # 0..255
        ns_global  = base_local + 256         # 256..511
        dbeam_ew = np.abs(tab["ibeam"] - ew_global)
        dbeam_ns = np.abs(tab["ibeam"] - ns_global)
        dbeam_min = np.minimum(dbeam_ew, dbeam_ns)
        which_arm = np.where(dbeam_ew <= dbeam_ns, 0, 1)  # 0=EW match, 1=NS match
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
    arm = int(which_arm[best_idx])  # 0=EW, 1=NS
    matched_global = ew_global if arm == 0 else ns_global
    dbeam_report = int(abs(int(tab["ibeam"][best_idx]) - matched_global))

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



def _earliest_drop_reason(gates: Dict[str, int], *, empty_gulp: bool|None=None) -> Optional[str]:
    """
    If empty_gulp is True, T1 returned zero rows this gulp.
    If False, T1 had rows but (potentially) no match.
    If None, we don't know (keep legacy behavior).
    """
    g1 = gates.get("G1_parsed", -1)
    g2 = gates.get("G2_T1_detected", -1)

    if g1 == 0:
        # distinguish the two early cases if we know it
        return "t1_empty" if empty_gulp else "no_parse"
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
        self._source_file_path: Optional[str] = None
        self._source_mtime: Optional[float] = None

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
        self._lock = threading.RLock()


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

    @staticmethod
    def _parse_injection_file(path: str):
        """
        Yield dicts from
        `injections_for_audit.txt`. Skips lines starting with '#'.
        """
        if not os.path.exists(path):
            return
        with open(path, "r") as fh:
            for ln in fh:
                ln = ln.strip()
                if not ln or ln.startswith("#"):
                    continue
                # Expect 7 tokens: MJD Beam DM SNR Width_fwhm spec_ind FRBno
                toks = ln.split()
                if len(toks) < 7:
                    continue
                try:
                    mjd        = float(toks[0])
                    beam       = int(toks[1])
                    dm         = float(toks[2])
                    snr        = float(toks[3])
                    width_fwhm = float(toks[4])
                    spec_ind   = float(toks[5])
                    frbno      = str(toks[6])
                except Exception:
                    continue
                yield dict(
                    MJD=mjd, Beam=beam, DM=dm, SNR=snr,
                    Width_fwhm=width_fwhm, spec_ind=spec_ind, FRBno=frbno
                )

    def _injections_csv_path(self) -> Path:
        return self.audit_dir / "injections.csv"

    def _coerce_int_or(self, v, default=""):
        try:
            if v == "" or pd.isna(v):
                return default
            return int(v)
        except Exception:
            return default

    def _coerce_float_or(self, v, default=""):
        try:
            if v == "" or pd.isna(v):
                return default
            return float(v)
        except Exception:
            return default

    def _coerce_str_or(self, v, default=""):
        try:
            if v == "" or pd.isna(v):
                return default
            return str(v)
        except Exception:
            return default

    def _reconstruct_state_from_snapshot_row(self, row: pd.Series) -> dict:
        """Turn a row from injections.csv back into our in-memory state dict."""
        st = dict(
            inj=dict(
                MJD=self._coerce_float_or(row.get("inj_mjd", "")),
                Beam=self._coerce_int_or(row.get("inj_beam", "")),
                DM=self._coerce_float_or(row.get("inj_dm", "")),
                SNR=self._coerce_float_or(row.get("inj_snr_req", "")),
                Width_fwhm=self._coerce_float_or(row.get("inj_width", "")),
                spec_ind=self._coerce_float_or(row.get("inj_spec_ind", "")),
                FRBno=self._coerce_str_or(row.get("inj_frbno", "")),
            ),
            gates={
                "G0_injected": self._coerce_int_or(row.get("G0_injected", -1), -1),
                "G1_parsed":   self._coerce_int_or(row.get("G1_parsed",   -1), -1),
                "G2_T1_detected": self._coerce_int_or(row.get("G2_T1_detected", -1), -1),
                "G3_beam_kept":   self._coerce_int_or(row.get("G3_beam_kept",   -1), -1),
                "G4_clustered":   self._coerce_int_or(row.get("G4_clustered",   -1), -1),
                "G5_filters_passed": self._coerce_int_or(row.get("G5_filters_passed", -1), -1),
                "G6_cooldown_ok": self._coerce_int_or(row.get("G6_cooldown_ok", -1), -1),
                "G7_triggered":   self._coerce_int_or(row.get("G7_triggered",   -1), -1),
            },
            candname=self._coerce_str_or(row.get("candname", "")),
            drop_reason=self._coerce_str_or(row.get("drop_reason", "")),
            t1_best=dict(
                t1_best_snr=self._coerce_float_or(row.get("t1_best_snr", "")),
                t1_best_dm=self._coerce_float_or(row.get("t1_best_dm", "")),
                t1_best_ibox=self._coerce_int_or(row.get("t1_best_ibox", "")),
                t1_best_ibeam=self._coerce_int_or(row.get("t1_best_ibeam", "")),
                t1_dt_sec=self._coerce_float_or(row.get("t1_dt_sec", "")),
                t1_dbeam=self._coerce_int_or(row.get("t1_dbeam", "")),
                t1_ddm=self._coerce_float_or(row.get("t1_ddm", "")),
            ),
            stats=dict(
                peak_snr=self._coerce_float_or(row.get("peak_snr", "")),
                peak_dm=self._coerce_float_or(row.get("peak_dm", "")),
                peak_ibox=self._coerce_int_or(row.get("peak_ibox", "")),
                peak_ibeam=self._coerce_int_or(row.get("peak_ibeam", "")),
                cooldown_wait_s=self._coerce_float_or(row.get("cooldown_wait_s", "")),
                t2_output_file=self._coerce_str_or(row.get("t2_output_file", "")),
                exception=self._coerce_str_or(row.get("exception", "")),
            ),
        )
        # If G3 already decided, freeze further G1/G2-only noise
        if st["gates"].get("G3_beam_kept", -1) in (0, 1):
            st["_frozen"] = True
        return st

    def bootstrap_from_snapshot(self) -> int:
        """
        Load self.injections from injections.csv if present.
        Idempotent: does not write audit_log and does not touch disk.
        Returns number of rows restored.
        """
        p = self._injections_csv_path()
        if not p.exists():
            return 0

        try:
            df = pd.read_csv(p)
        except Exception as e:
            self.log.warning(f"[AUDIT] bootstrap: failed to read {p}: {e}")
            return 0

        restored = 0
        with self._lock:
            for _, row in df.iterrows():
                inj_id = self._coerce_str_or(row.get("inj_id", ""))
                if not inj_id:
                    continue
                st = self._reconstruct_state_from_snapshot_row(row)
                self.injections[inj_id] = st
                restored += 1

        self.log.info(f"[AUDIT] bootstrap: restored {restored} injections from snapshot")
        print(f"[AUDIT] bootstrap: restored {restored} injections from snapshot")
        return restored

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

    def attach_injection_source(self, path: str):
        """
        Register the on-disk file we should watch and refresh from.
        Does not force an immediate load; call refresh_from_source() to load.
        """
        with self._lock:
            self._source_file_path = str(path)
            try:
                self._source_mtime = os.path.getmtime(self._source_file_path)
            except Exception:
                self._source_mtime = None

    def refresh_from_source(self, force: bool = False) -> int:
        """
        If the source file changed (mtime bumped) or force=True,
        parse it and seed/update injections. Returns number of
        new/updated rows applied.
        - New rows: seed as new injection (G0=1, others -1).
        - Existing rows (same inj_id): update meta fields (SNR, width, spec_ind, FRBno)
        without touching gate states.
        """
        with self._lock:
            if not self._source_file_path:
                return 0

            try:
                mtime = os.path.getmtime(self._source_file_path)
            except Exception:
                return 0

            if (not force) and (self._source_mtime is not None) and (mtime <= self._source_mtime):
                return 0  # no change

        applied = 0
        metas = list(self._parse_injection_file(self._source_file_path))
        with self._lock:
            for meta in metas:
                inj_id = _stable_inj_id(meta["MJD"], meta["Beam"], meta["DM"])
                if inj_id not in self.injections:
                    # brand-new
                    self.seed_injection(inj_id, meta)
                    applied += 1
                else:
                    # update meta (non-destructive to gates/candname/stats)
                    st = self.injections[inj_id]
                    inj = st.get("inj", {})
                    inj["SNR"]        = float(meta.get("SNR", inj.get("SNR", ""))) if str(meta.get("SNR", "")) != "" else inj.get("SNR", "")
                    inj["Width_fwhm"] = float(meta.get("Width_fwhm", inj.get("Width_fwhm", ""))) if str(meta.get("Width_fwhm", "")) != "" else inj.get("Width_fwhm", "")
                    inj["spec_ind"]   = float(meta.get("spec_ind", inj.get("spec_ind", ""))) if str(meta.get("spec_ind", "")) != "" else inj.get("spec_ind", "")
                    inj["FRBno"]      = str(meta.get("FRBno", inj.get("FRBno", "")))
                    st["inj"] = inj
                    applied += 1

            # bump mtime and rewrite rolling CSV snapshot
            self._source_mtime = mtime
            self._write_injections_csv()
        return applied

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
        fh = open(p, "a+")
        try:
            fcntl.flock(fh, fcntl.LockFlags.LOCK_EX if hasattr(fcntl, "LockFlags") else fcntl.LOCK_EX)
            fh.seek(0, os.SEEK_END)
            empty = fh.tell() == 0
            if empty:
                fh.write("# Audit timeline for injected FRBs vs pipeline gates.\n")
                fh.write("# update_from_tab stamps G1,G2; finalize stamps G3..G7.\n")
                fh.write(",".join(_AUDIT_LOG_HEADER) + "\n")
        finally:
            try:
                fcntl.flock(fh, fcntl.LockFlags.LOCK_UN if hasattr(fcntl, "LockFlags") else fcntl.LOCK_UN)
            finally:
                fh.close()

    def _state_json_path(self, inj_id: str) -> Path:
        if self.state_dir is None:
            raise RuntimeError("state dir requested but persist_json=False")
        return self.state_dir / f"{inj_id}.json"

    def _flock_append_csvrow(self, path: Path, row_dict: Dict[str, Any]):
        with self._lock:
            df = pd.DataFrame([row_dict], columns=_AUDIT_LOG_HEADER)
            csv_buf = io.StringIO()
            # Force quotes around string fields so lists like "[1,2,3]" don't explode columns
            df.to_csv(csv_buf, index=False, header=False, quoting=csv.QUOTE_MINIMAL)
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
        # empty tab or no mjds column
        return None

    def _iter_injections_in_window(
        self,
        mjd_center: float,
    ) -> Iterable[Tuple[str, Dict[str, Any]]]:
        """
            Yield all injections within +/- time_window_s seconds of mjd_center.
        """
        window_days = self.time_window_s / 86400.0
        mjd_min = mjd_center - window_days
        mjd_max = mjd_center + window_days
        with self._lock:
            items = list(self.injections.items())
        for inj_id, st in items:
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
            with self._lock:
                pj = self._state_json_path(inj_id)
                tmp = pj.with_suffix(".json.tmp")
                # snapshot state under lock
                st_snapshot = json.loads(json.dumps(self.injections[inj_id]))
            # write outside the lock
            with open(tmp, "w") as fh:
                json.dump(st_snapshot, fh, indent=2)
            os.replace(tmp, pj)
            msg = f"[AUDIT] persisted state JSON for inj_id={inj_id} -> {pj}"
            self.log.debug(msg)
            print(msg)
        except Exception as e:
            self.log.warning(f"[AUDIT] failed to persist state JSON for inj_id={inj_id}: {e}")
            print(f"[AUDIT] failed to persist state JSON for inj_id={inj_id}: {e}")

    # public API
    def seed_injection(self, inj_id: str, inj_meta: Dict[str, Any]):
        """
        Called by the injection code or by legacy ingestion.
        Create initial record with gates all -1 except G0 (injected=1).
        """
        with self._lock:
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
            if mjd_now is None:
                self.log.warning("[AUDIT][update] cannot determine mjd_now from tab likely because it is empty.")
                mjd_now = Time.now().mjd
                self.log.info("[AUDIT][update] empty tab; using now() for mjd_center")
                print("[AUDIT][update] empty tab; using now() for mjd_center")

            self.log.info(f"[AUDIT][update] start gulp={gulp} host={host} len(tab)={len(tab)} mjd_now={mjd_now}")
            print(f"[AUDIT][update] start gulp={gulp} host={host} len(tab)={len(tab)} mjd_now={mjd_now}")

            # compute once
            npoints_tab = int(len(tab))
            empty_gulp = (npoints_tab == 0)
            gate_G1_parsed = 1 if npoints_tab > 0 else 0

            for inj_id, st in self._iter_injections_in_window(mjd_now):
                inj = st["inj"]
                match_idx, t1_best = _match_rows_for_injection(
                    tab,
                    inj_mjd=float(inj["MJD"]),
                    inj_dm=float(inj["DM"]),
                    inj_beam=int(inj["Beam"]),
                    time_window_s=self.time_window_s,
                    dm_window=self.dm_window,
                    beam_window=self.beam_window,
                )
                gate_G2_T1_detected = 1 if match_idx is not None else 0

                with self._lock:
                    st2 = self.injections.get(inj_id)
                    if not st2 or st2.get("_frozen", False):
                        continue

                    if gate_G2_T1_detected == 1 and t1_best:
                        st2["t1_best"] = dict(t1_best)

                    _merge_gate(st2, "G1_parsed", gate_G1_parsed)
                    _merge_gate(st2, "G2_T1_detected", gate_G2_T1_detected)

                    er = _earliest_drop_reason(st2["gates"], empty_gulp=empty_gulp)
                    if er:
                        st2["drop_reason"] = er

                    st2["stats"]["npoints_tab"] = npoints_tab
                    st2["stats"]["npoints_after_flag"] = ""  # filled in finalize
                    st2["stats"]["nbeams_queue_snapshot"] = list(nbeams_queue_snapshot)
                    st2["stats"]["nbeams_queue_sum"] = int(nbeams_queue_sum)
                    st2["stats"]["runtime_thresholds"] = dict(runtime_thresholds)

                    if self.persist_json:
                        self._persist_state_json(inj_id)

                    fp = _state_fingerprint(st2)
                    already_finalized = st2["gates"].get("G3_beam_kept", -1) in (0, 1)
                    if not (npoints_tab == 0 and already_finalized) and fp != st2.get("_last_fp"):
                        meta = {
                            "gulp": gulp,
                            "host": host,
                            "npoints_tab": npoints_tab,
                            "npoints_after_flag": "",
                            "nbeams_queue_snapshot": json.dumps(list(nbeams_queue_snapshot)),
                            "nbeams_queue_sum": int(nbeams_queue_sum),
                            "drop_reason": st2.get("drop_reason", ""),
                        }
                        if t1_best:
                            meta.update(t1_best)
                        self._append_audit_log(inj_id, st2, meta, exception="")
                        st2["_last_fp"] = fp

            # write snapshot once
            self._write_injections_csv()

            self.log.info(f"[AUDIT][update] done: considered_injections={len(list(self.injections.keys()))}")
            print(f"[AUDIT][update] done: considered_injections={len(list(self.injections.keys()))}")

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
            if mjd_now is None:
                mjd_now = Time.now().mjd

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

            npoints_after_flag_val = (int(len(tab_after_beam_flag)) if tab_after_beam_flag is not None else "")

            for inj_id, st in self._iter_injections_in_window(mjd_now):
                inj = st["inj"]
                inj_mjd = float(inj["MJD"])
                inj_dm = float(inj["DM"])
                inj_beam = int(inj["Beam"])

                # compute decisions without mutating
                if st["gates"].get("G2_T1_detected", -1) == 1 and tab_after_beam_flag is not None:
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
                    gate_G3_beam_kept = -1

                peak_info = _extract_peak_info_for_injection(
                    tab_peak,
                    inj_mjd=inj_mjd,
                    inj_dm=inj_dm,
                    inj_beam=inj_beam,
                    time_window_s=self.time_window_s,
                    dm_window=self.dm_window,
                    beam_window=self.beam_window,
                ) if gate_G3_beam_kept == 1 else {}

                if gate_G3_beam_kept == 1:
                    match_idx_after, _tmp_after = _match_rows_for_injection(
                        tab_after_filters,
                        inj_mjd=inj_mjd,
                        inj_dm=inj_dm,
                        inj_beam=inj_beam,
                        time_window_s=self.time_window_s,
                        dm_window=self.dm_window,
                        beam_window=self.beam_window,
                    )
                else:
                    match_idx_after = None

                gate_G4_clustered = 1 if peak_info else 0
                gate_G5_filters_passed = 1 if match_idx_after is not None else 0

                max_nbeams_allowed = int(thresholds.get("max_nbeams", 40))
                min_timedelt = float(thresholds.get("min_timedelt", 60.0))
                if cooldown_wait_s is None:
                    gate_G6_cooldown_ok = 1
                else:
                    try:
                        cw = float(cooldown_wait_s)
                        gate_G6_cooldown_ok = 1 if (not np.isfinite(cw)) or (cw >= min_timedelt) else 0
                    except Exception:
                        gate_G6_cooldown_ok = 1

                gate_G7_triggered = 1 if bool(triggered) else 0

                # now mutate under lock
                with self._lock:
                    st2 = self.injections.get(inj_id)
                    if not st2 or st2.get("_frozen", False):
                        continue

                    _merge_gate(st2, "G3_beam_kept", gate_G3_beam_kept)
                    st2["stats"]["npoints_after_flag"] = npoints_after_flag_val

                    # early drop reasons
                    if not (st2["gates"].get("G1_parsed", -1) == 1 and st2["gates"].get("G2_T1_detected", -1) == 1):
                        er = _earliest_drop_reason(st2["gates"])
                        if er:
                            st2["drop_reason"] = er
                    elif st2["gates"].get("G3_beam_kept", -1) != 1:
                        if st2["gates"].get("G3_beam_kept", -1) == 0:
                            st2["drop_reason"] = "beam_flagged"
                    else:
                        _merge_gate(st2, "G4_clustered", gate_G4_clustered)
                        _merge_gate(st2, "G5_filters_passed", gate_G5_filters_passed)
                        _merge_gate(st2, "G6_cooldown_ok", gate_G6_cooldown_ok)
                        _merge_gate(st2, "G7_triggered", gate_G7_triggered)

                        # decide drop_reason hierarchy
                        er = _earliest_drop_reason(st2["gates"])
                        if er:
                            st2["drop_reason"] = er
                        elif st2["gates"]["G3_beam_kept"] == 0:
                            st2["drop_reason"] = "beam_flagged"
                        elif st2["gates"]["G4_clustered"] == 0:
                            st2["drop_reason"] = "no_cluster_peak"
                        elif st2["gates"]["G5_filters_passed"] == 0:
                            st2["drop_reason"] = "filtered_out"
                        elif nbeams_queue_sum > max_nbeams_allowed:
                            st2["drop_reason"] = f"nbeams_gate_exceeded({nbeams_queue_sum}>{max_nbeams_allowed})"
                        elif st2["gates"]["G6_cooldown_ok"] == 0:
                            st2["drop_reason"] = f"cooldown({cooldown_wait_s:.2f}s<{min_timedelt}s)"
                        elif st2["gates"]["G7_triggered"] == 0:
                            st2["drop_reason"] = "not_triggered"
                        else:
                            st2["drop_reason"] = ""

                        if candname_for_injection:
                            st2["candname"] = candname_for_injection

                        st2["stats"]["nclusters"] = int(len(tab_peak)) if tab_peak is not None else 0
                        st2["stats"]["cluster_size"] = peak_info.get("cntc", "")
                        st2["stats"]["cntb"] = peak_info.get("cntb", "")
                        st2["stats"]["cntc"] = peak_info.get("cntc", "")
                        st2["stats"]["peak_snr"] = peak_info.get("snr", "")
                        st2["stats"]["peak_dm"] = peak_info.get("dm", "")
                        st2["stats"]["peak_ibox"] = peak_info.get("ibox", "")
                        st2["stats"]["peak_ibeam"] = peak_info.get("ibeam", "")
                        st2["stats"]["nbeams_this_gulp"] = int(nbeams_this_gulp) if nbeams_this_gulp is not None else ""
                        st2["stats"]["nbeams_queue_snapshot"] = list(nbeams_queue_snapshot)
                        st2["stats"]["nbeams_queue_sum"] = int(nbeams_queue_sum)
                        st2["stats"]["cooldown_wait_s"] = cooldown_wait_s if cooldown_wait_s else ""
                        st2["stats"]["t2_output_file"] = thresholds.get("t2_output_file", "")
                        st2["stats"]["exception"] = thresholds.get("exception", "")

                    if self.persist_json:
                        self._persist_state_json(inj_id)

                    meta = {
                        "gulp": gulp,
                        "host": host,
                        "nclusters": st2["stats"].get("nclusters", ""),
                        "cluster_size": st2["stats"].get("cluster_size", ""),
                        "cntb": st2["stats"].get("cntb", ""),
                        "cntc": st2["stats"].get("cntc", ""),
                        "peak_snr": st2["stats"].get("peak_snr", ""),
                        "peak_dm": st2["stats"].get("peak_dm", ""),
                        "peak_ibox": st2["stats"].get("peak_ibox", ""),
                        "peak_ibeam": st2["stats"].get("peak_ibeam", ""),
                        "nbeams_this_gulp": st2["stats"].get("nbeams_this_gulp", ""),
                        "nbeams_queue_snapshot": json.dumps(list(nbeams_queue_snapshot)),
                        "nbeams_queue_sum": int(nbeams_queue_sum),
                        "cooldown_wait_s": st2["stats"].get("cooldown_wait_s", ""),
                        "t2_output_file": st2["stats"].get("t2_output_file", ""),
                        "drop_reason": st2.get("drop_reason", ""),
                        "exception": st2["stats"].get("exception", ""),
                    }
                    if st2.get("t1_best"):
                        meta.update(st2["t1_best"])

                    fp = _state_fingerprint(st2)
                    if fp != st2.get("_last_fp"):
                        self._append_audit_log(inj_id, st2, meta, exception=st2["stats"].get("exception", ""))
                        st2["_last_fp"] = fp

                    if st2["gates"].get("G3_beam_kept", -1) in (0, 1):
                        st2["_frozen"] = True

            self._write_injections_csv()

            self.log.info(f"[AUDIT][finalize] done: considered_injections={len(list(self.injections.keys()))}")
            print(f"[AUDIT][finalize] done: considered_injections={len(list(self.injections.keys()))}")

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
        with self._lock:
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
        """
        # snapshot rows under lock
        with self._lock:
            rows: List[Dict[str, Any]] = []
            for inj_id, st in self.injections.items():
                inj = st["inj"]
                gates = st["gates"]
                telem = st.get("stats", {})
                t1_best = st.get("t1_best", {})

                rows.append({
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
                })

        p = self._injections_csv_path()
        tmp = p.with_suffix(".tmp")
        df = pd.DataFrame(rows, columns=_INJECTIONS_HEADER)
        df.to_csv(tmp, index=False)
        os.replace(tmp, p)

        self.log.info(f"[AUDIT][write_injections_csv] wrote {len(rows)} rows to {p}")
        print(f"[AUDIT][write_injections_csv] wrote {len(rows)} rows to {p}")

        if self.persist_json:
            # persist each JSON outside the CSV lock to avoid deadlocks
            for inj_id in list(self.injections.keys()):
                self._persist_state_json(inj_id)
