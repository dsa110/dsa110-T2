import T2
import pytest
import os.path
import sys
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

