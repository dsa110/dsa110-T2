import T2
import pytest
import os.path

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

    # T2 cluster 
    clusterer = T2.cluster_heimdall.cluster_data(tab, min_cluster_size=10, min_samples=10,
                                                 selectcols=['itime', 'dm', 'ibox', 'ibeam'],
                                                 metric='euclidean', allow_single_cluster=True,
                                                 return_clusterer=True)
    tab2 = T2.cluster_heimdall.get_peak(tab)
    
    # send T2 cluster results to outputfile
    row, candname, prev_trig_time = T2.cluster_heimdall.dump_cluster_results_json(tab, outputfile=outputfile, output_cols=['mjds', 'snr', 'ibox', 'ibeam', 'dm'])

    assert candname is not None
    
#    if plot: 
#        T2.cluster_heimdall.plot_giants(tab, plot_dir=plot_dir) # plot giants      
#        T2.cluster_heimdall.plot_clustered(clusterer, clsnr, snrs, data, tab, cols=['itime', 'idm', 'ibox'], plot_dir=plot_dir) # plot cluster results  


def test_cluster_and_plot():
    """
    Run socket.cluster_and_plot
    """

    candsfile = os.path.join(_install_dir, 'data/giants_1.cand')

    # read in giants 
    tab = T2.cluster_heimdall.parse_candsfile(candsfile)

    inputname = '200101aaaa'
    lastname, prev_trig_time = T2.socket.cluster_and_plot(tab, 0, outroot='tests/test_', lastname=inputname)

    assert lastname is not inputname   # if too many, as for giants_1.cand


def test_lastname():
    """
    Run socket.cluster_and_plot and use lastname
    """

    candsfile = os.path.join(_install_dir, 'data/giants_1.cand')

    # read in giants 
    tab = T2.cluster_heimdall.parse_candsfile(candsfile)

    lastname, prev_trig_time = T2.socket.cluster_and_plot(tab, 0, outroot='tests/test_')
    lastname2, prev_trig_time = T2.socket.cluster_and_plot(tab, 0, outroot='tests/test_', lastname=lastname)

    assert lastname != lastname2

def test_import():
    """
    Test import
    """

    import T2.cluster_heimdall
    import T2.socket
    import T2.T2

    assert True
    