#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Offline T2 replay script.

Replays T2 logic on T1 candidates stored in a CSV, **without** sending real
voltage triggers, etcd writes, Slack alerts, or touching live T2 directories.

Everything this script outputs is sandboxed into an offline directory.

Assumes:
    - The original T2 'socket.py' has been renamed to 't2_socket.py'
      and is imported as `from T2 import t2_socket`.
"""

import argparse
import os
import re
import sys
import logging
import socket  # stdlib; needed so pyarrow/pandas don't accidentally import T2.socket

import pandas as pd
from astropy.table import Table
from astropy.time import Time

from T2 import t2_socket
from T2 import cluster_heimdall


# ============================================================
# Logging Setup
# ============================================================

logger = logging.getLogger("offline_t2")
logger.setLevel(logging.INFO)

handler = logging.StreamHandler(sys.stdout)
handler.setFormatter(logging.Formatter("[%(asctime)s][%(levelname)s] %(message)s"))
logger.addHandler(handler)


# ============================================================
# Global State
# ============================================================

OFFLINE_TRIGGERS = []
_ORIG_SEND_TRIGGER = cluster_heimdall.send_trigger


# ============================================================
# Monkeypatch: Record Voltage Triggers
# ============================================================

def offline_send_trigger(output_dict=None, outputfile=None):
    """
    Override for cluster_heimdall.send_trigger.

    Instead of talking to etcd/Slack/corr, this records the would-be
    voltage trigger into OFFLINE_TRIGGERS, including the current state
    of the T2 nbeams_queue.
    """
    if not output_dict:
        return

    candname = list(output_dict.keys())[0]
    info = dict(output_dict[candname])

    info["candname"] = candname
    info["event_type"] = "voltage"

    # Capture actual T2 nbeams_queue at trigger time
    try:
        q = list(t2_socket.nbeams_queue)
        info["nbeams_sum10"] = int(sum(q))
        info["nbeams_queue"] = q
    except Exception:
        info["nbeams_sum10"] = None
        info["nbeams_queue"] = None

    OFFLINE_TRIGGERS.append(info)
    logger.debug(
        f"Recorded offline trigger: {candname}, "
        f"nbeams_queue={info['nbeams_queue']}"
    )


def enable_offline_mode():
    cluster_heimdall.send_trigger = offline_send_trigger


def disable_offline_mode():
    cluster_heimdall.send_trigger = _ORIG_SEND_TRIGGER


# ============================================================
# Helpers
# ============================================================

def extract_date_from_filename(path):
    """
    Extract YYYYMMDD from a filename like t1_candidates_20251111.csv.
    """
    m = re.search(r"(\d{8})", os.path.basename(path))
    return m.group(1) if m else "unknown"


def df_to_astropy_tab(df_g):
    """
    Convert one gulp group's rows to the Astropy Table format expected by T2.

    Required columns:
        snr, if, itime, mjds, ibox, idm, dm, ibeam
    """
    cols = ["snr", "if", "itime", "mjds", "ibox", "idm", "dm", "ibeam"]
    missing = [c for c in cols if c not in df_g]
    if missing:
        raise ValueError(f"Missing required T1 CSV columns: {missing}")

    df2 = df_g[cols].copy()
    for c in ["ibeam", "idm", "ibox", "itime"]:
        df2[c] = df2[c].astype(int)

    return Table.from_pandas(df2)


# ============================================================
# Main Replay Logic
# ============================================================

def run_offline_t2(
    t1_csv,
    output_csv,
    outroot,
    source_catalog=None,
    dry_run=False,
):
    logger.info(f"Reading T1 CSV: {t1_csv}")

    df = pd.read_csv(t1_csv)
    if "gulp" not in df:
        raise ValueError("Input CSV must contain a 'gulp' column.")
    if "mjds" not in df:
        raise ValueError("Input CSV must contain an 'mjds' column.")

    # Sort by time first, then gulp, then ibeam to respect temporal order
    sort_cols = ["mjds", "gulp"]
    if "ibeam" in df:
        sort_cols.append("ibeam")
    df = df.sort_values(sort_cols)

    OFFLINE_TRIGGERS.clear()

    # Reset the actual T2 nbeams_queue so replay starts from a clean state
    try:
        t2_socket.nbeams_queue.clear()
    except Exception:
        logger.warning("Could not clear t2_socket.nbeams_queue; using existing state.")

    try:
        lastname = cluster_heimdall.names.get_lastname()
    except Exception:
        lastname = None

    enable_offline_mode()
    prev_trig_time = None  # time gating explicitly ignored in offline replay

    # Group gulps the way the online system would see them:
    # all rows with the same (time, gulp) together, processed in time order.
    if "recv_ts_iso" in df.columns:
        group_cols = ["recv_ts_iso", "gulp"]
    else:
        # fallback: (mjds, gulp); assumes mjds per gulp is effectively constant
        group_cols = ["mjds", "gulp"]

    logger.info(f"Grouping by {group_cols} in chronological order.")

    # Because df is already sorted by mjds, gulp, ibeam, groupby(sort=False)
    # will yield groups in that order of first appearance.
    grouped = df.groupby(group_cols, sort=False)

    n_groups = len(grouped)
    logger.info(f"Found {n_groups} gulp groups to process.")

    for idx, (keys, df_g) in enumerate(grouped, start=1):
        if df_g.empty:
            continue

        # keys is either (recv_ts_iso, gulp) or (mjds, gulp)
        if isinstance(keys, tuple):
            time_key, gulp = keys
        else:
            time_key, gulp = None, keys

        if idx % 200 == 0:
            logger.info(
                f"Processing group {idx}/{n_groups}: "
                f"gulp={gulp}, time_key={time_key}"
            )

        tab = df_to_astropy_tab(df_g)

        try:
            lastname, trigtime, triggered = t2_socket.cluster_and_plot(
                tab,
                gulp=int(gulp),                    # keep the real gulp number
                selectcols=["itime", "idm", "ibox"],
                outroot=outroot,
                plot_dir=None,
                trigger=True,
                lastname=lastname,
                max_ncl=None,
                cat=source_catalog,
                beam_model=None,
                coords=None,
                snrs=None,
                prev_trig_time=prev_trig_time,
            )
        except Exception as e:
            logger.error(
                f"Error in cluster_and_plot for gulp={gulp}, time_key={time_key}: {e}"
            )
            continue

    disable_offline_mode()

    n_triggers = len(OFFLINE_TRIGGERS)

    if dry_run:
        if n_triggers:
            logger.info(
                f"[dry-run] Would write {n_triggers} triggers to {output_csv} "
                f"(not writing due to --dry-run)."
            )
        else:
            logger.warning("[dry-run] No triggers recorded. Nothing would be written.")
        return

    # Normal mode: write CSV
    if n_triggers:
        dfout = pd.DataFrame(OFFLINE_TRIGGERS)

        # Add ISO timestamp from mjds as a final step
        if "mjds" in dfout.columns:
            dfout["recv_ts_iso"] = dfout["mjds"].apply(
                lambda x: Time(x, format="mjd").isot
            )
        else:
            dfout["recv_ts_iso"] = None

        # Ensure useful columns appear early
        cols = list(dfout.columns)
        for k in ["event_type", "candname", "recv_ts_iso", "nbeams_sum10", "nbeams_queue"]:
            if k in cols:
                cols.insert(0, cols.pop(cols.index(k)))
        dfout = dfout[cols]

        os.makedirs(os.path.dirname(output_csv), exist_ok=True)
        dfout.to_csv(output_csv, index=False)
        logger.info(f"Wrote {n_triggers} triggers to {output_csv}")
    else:
        pd.DataFrame(
            {"event_type": [], "candname": [], "recv_ts_iso": []}
        ).to_csv(output_csv, index=False)
        logger.warning(f"No triggers found. Wrote empty file: {output_csv}")


# ============================================================
# CLI
# ============================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Offline T2 replay — re-run T2 trigger decisions on T1 candidates."
    )

    parser.add_argument("--t1-csv", required=True, help="Input T1 candidates CSV.")
    parser.add_argument("--output-dir", default="offline_T2",
                        help="Directory for offline outputs.")

    parser.add_argument("--output-csv",
                        default=None,
                        help="Optional explicit output CSV path.")

    parser.add_argument("--tmp-outroot",
                        default=None,
                        help="Directory for intermediate T2 files (cand/json). "
                             "Default: <output-dir>/t2_tmp")

    parser.add_argument("--source-catalog",
                        default=None,
                        help="Optional source catalog for T2 source rejection.")

    parser.add_argument("--quiet", action="store_true",
                        help="Suppress INFO logs; only warnings/errors.")

    parser.add_argument("--dry-run", action="store_true",
                        help="Run replay and log what would be written, "
                             "but do not write the offline trigger CSV.")

    return parser.parse_args()


def main():
    args = parse_args()

    if args.quiet:
        logger.setLevel(logging.WARNING)
    else:
        logger.setLevel(logging.INFO)

    t1_csv = args.t1_csv
    outdir = os.path.abspath(args.output_dir)
    os.makedirs(outdir, exist_ok=True)

    date_str = extract_date_from_filename(t1_csv)

    if args.output_csv is None:
        output_csv = os.path.join(outdir, f"t2_replay_offline_triggers_{date_str}.csv")
    else:
        output_csv = os.path.abspath(args.output_csv)

    if args.tmp_outroot is None:
        outroot = os.path.join(outdir, "t2_tmp") + os.sep
    else:
        outroot = os.path.abspath(args.tmp_outroot) + os.sep

    logger.info("========== Offline T2 Replay ==========")
    logger.info(f"T1 CSV:             {t1_csv}")
    logger.info(f"Output CSV:         {output_csv}")
    logger.info(f"T2 temp outroot:    {outroot}")
    logger.info(f"Dry run mode:       {args.dry_run}")
    if args.source_catalog:
        logger.info(f"Using source catalog: {args.source_catalog}")

    run_offline_t2(
        t1_csv=t1_csv,
        output_csv=output_csv,
        outroot=outroot,
        source_catalog=args.source_catalog,
        dry_run=args.dry_run,
    )


if __name__ == "__main__":
    main()
