import T2
import pytest
import os
import os.path
import sys
import glob
import shutil
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd

_install_dir = os.path.abspath(os.path.dirname(__file__))

def test_T2():
    """
    Test T2 functions. Read in heimdall giants, cluster them, send results to outputfile and save plots. 

    Parameters
    ----------
    mode : String. "file" or "socket". 
    outputfile : default is "T2_output.txt".
    candsfile : default is "giants.cand.".
    host, ports : use the same host, port with heimdall -coincidencer HOST:PORT 
    plot : default is True.
    plot_dir : default is "./".
    """

    outputfile = os.path.join(_install_dir, 'data/T2_output.txt')
    candsfile = os.path.join(_install_dir, 'data/giants_1.cand')

    # read in giants 
    tab = T2.cluster_heimdall.parse_candsfile(candsfile)
    flagged_tab = T2.cluster_heimdall.flag_beams(tab)
    tab = flagged_tab
    
    # T2 clustering
    T2.cluster_heimdall.cluster_data(tab, allow_single_cluster=True)
    tab2 = T2.cluster_heimdall.get_peak(tab)
    nbeams_gulp = T2.cluster_heimdall.get_nbeams(tab)
    #print number of rows of astropy table
    print(f"Before Clustering: We have {len(tab)} candidates from {nbeams_gulp} beams")
    tab3 = T2.cluster_heimdall.filter_clustered(tab2)
    print(f"After Clustering: We have {len(tab3)} candidates from {nbeams_gulp} beams")
    
    # send T2 cluster results to outputfile
    row, candname, trigtime = T2.cluster_heimdall.dump_cluster_results_json(tab3, outputfile=outputfile)

    assert candname is not None

def test_cluster_and_plot():
    """
    Run socket.cluster_and_plot
    """

    candsfile = os.path.join(_install_dir, 'data/giants_1.cand')

    # read in giants 
    tab = T2.cluster_heimdall.parse_candsfile(candsfile)

    lastname, trigtime, triggered = T2.socket.cluster_and_plot(tab, 0) 

    assert lastname is None   # if too many, as for giants_1.cand


def test_lastname():
    """
    Run socket.cluster_and_plot and use lastname
    """

    candsfile = os.path.join(_install_dir, 'data/T1_output1744907347.csv')

    # read in giants 
    tab = T2.cluster_heimdall.parse_candsfile(candsfile)

    lastname, trigtime, triggered = T2.socket.cluster_and_plot(tab, 0, outroot='tests/test_', max_ncl=100000)
    lastname2, trigtime, triggered = T2.socket.cluster_and_plot(tab, 0, outroot='tests/test_', max_ncl=100000, lastname=lastname)
    
    assert lastname is not None
    assert lastname != lastname2


def test_atomic_csv_write_concurrent():
    """
    Verify _atomic_csv_write produces a valid CSV even when many threads
    write to the same path simultaneously.
    """
    tmpdir = tempfile.mkdtemp(prefix="t2_atomic_test_")
    target = os.path.join(tmpdir, "concurrent.csv")
    n_threads = 8
    n_writes_per_thread = 20

    def writer(thread_id):
        for i in range(n_writes_per_thread):
            df = pd.DataFrame({
                "thread": [thread_id] * 5,
                "iter": [i] * 5,
                "value": list(range(5)),
            })
            T2.socket._atomic_csv_write(df, target)

    with ThreadPoolExecutor(max_workers=n_threads) as pool:
        futs = [pool.submit(writer, t) for t in range(n_threads)]
        for f in as_completed(futs):
            f.result()

    # The file must exist and be a well-formed CSV (last writer wins)
    assert os.path.exists(target)
    df = pd.read_csv(target)
    assert list(df.columns) == ["thread", "iter", "value"]
    assert len(df) == 5

    shutil.rmtree(tmpdir)


def test_aggregate_locking():
    """
    Submit several cluster_and_plot calls concurrently and verify the
    aggregate CSVs are well-formed (no interleaved or garbled rows).
    """
    candsfile = os.path.join(_install_dir, "data/T1_output1744907347.csv")
    tab = T2.cluster_heimdall.parse_candsfile(candsfile)

    tmpdir = tempfile.mkdtemp(prefix="t2_lock_test_")
    outroot = os.path.join(tmpdir, "lock_")
    n_threads = 4

    def run_one(idx):
        return T2.socket.cluster_and_plot(
            tab, gulp=idx, outroot=outroot, max_ncl=100000
        )

    with ThreadPoolExecutor(max_workers=n_threads) as pool:
        futs = [pool.submit(run_one, i) for i in range(n_threads)]
        for f in as_completed(futs):
            f.result()

    agg = os.path.join(tmpdir, "lock_cluster_output.csv")
    assert os.path.exists(agg), "aggregate CSV was not created"

    df = pd.read_csv(agg)
    assert len(df) > 0, "aggregate CSV is empty"

    # Every row must have the expected column count (no interleaved lines)
    with open(agg) as fh:
        header = fh.readline()
        n_cols = len(header.strip().split(","))
        for lineno, line in enumerate(fh, start=2):
            parts = line.strip().split(",")
            assert len(parts) == n_cols, (
                f"line {lineno}: expected {n_cols} columns, got {len(parts)}"
            )

    # Daily CSV must also be valid
    daily_csvs = glob.glob(os.path.join(tmpdir, "lock_[0-9]*.csv"))
    for csv_path in daily_csvs:
        df_daily = pd.read_csv(csv_path)
        assert len(df_daily) > 0

    shutil.rmtree(tmpdir)

