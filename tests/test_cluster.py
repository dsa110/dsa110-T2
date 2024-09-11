import pytest
import os.path
from astropy import table
from T2 import cluster_heimdall, plotting

_install_dir = os.path.abspath(os.path.dirname(__file__))

@pytest.fixture(scope="module")
def tab():
    candsfile = os.path.join(_install_dir, 'data/giants.cand')
    tab = cluster_heimdall.parse_candsfile(candsfile)
    return tab


@pytest.fixture(scope="module")
def tabs():
    candsfile1 = os.path.join(_install_dir, 'data/giants_1.cand')
    candsfile2 = os.path.join(_install_dir, 'data/giants_2.cand')
    tab1 = cluster_heimdall.parse_candsfile(candsfile1)
    tab2 = cluster_heimdall.parse_candsfile(candsfile2)
    return tab1, tab2


def test_parse(tab):
    assert len(tab) == 3411
    assert len(tab[0]) == 8


def test_cluster1(tab):
    cluster_heimdall.cluster_data(tab, return_clusterer=False)
    assert len(tab[0]) == 11


def test_cluster_2arm(tabs):
    tab = table.hstack(tabs)
    # find cluster ignoring ibeam
    cluster_heimdall.cluster_data(tab, return_clusterer=False, metric='euclidean', min_cluster_size=2, min_samples=5,
                                  allow_single_cluster=True, selectcols=["itime", "idm", "ibox"])
    tabp = cluster_heimdall.get_peak(tab)  # returns peak snr of each cluster for two beam sets (0-255, 256-511)
    tabpf = cluster_heimdall.filter_clustered(tabp)


def test_cluster_2arm_alt(tabs):
    tab1, tab2 = tabs
    # find cluster per beam set (ns, ew)
    cluster_heimdall.cluster_data(tab1, return_clusterer=False, metric='euclidean', min_cluster_size=2, min_samples=5,
                                  allow_single_cluster=True, selectcols=["itime", "idm", "ibox", "ibeam"])
    cluster_heimdall.cluster_data(tab2, return_clusterer=False, metric='euclidean', min_cluster_size=2, min_samples=5,
                                  allow_single_cluster=True, selectcols=["itime", "idm", "ibox", "ibeam"])

    # reduce to snr peak per cluster per arm
    tab1p = cluster_heimdall.get_peak(tab1)
    tab2p = cluster_heimdall.get_peak(tab2)

    # TODO: decide whether to coincidence again, then filter or vice versa
    tab1pf = cluster_heimdall.filter_clustered(tab1p)
    tab2pf = cluster_heimdall.filter_clustered(tab2p)


def test_peak(tab):
    cluster_heimdall.cluster_data(tab, return_clusterer=False)

    tab2 = cluster_heimdall.get_peak(tab)
    assert len(tab2) == 1
    assert len(tab2[0]) == 11
    assert tab2 == tab[1380]
    assert tab2['snr'] == 117.613
    # TODO
    # assert cb ==
    # assert cc ==


def test_filter(tab):
    cluster_heimdall.cluster_data(tab, return_clusterer=False)
    tab2 = cluster_heimdall.get_peak(tab)
    assert len(tab2) == 1
    tab2 = cluster_heimdall.filter_clustered(tab, min_snr=1000)  # remove all
    assert len(tab2) == 0


def test_json(tab):
    outfile = os.path.join(_install_dir, 'trigger.json')
    cluster_heimdall.cluster_data(tab, return_clusterer=False)

    tab2 = cluster_heimdall.get_peak(tab)
    cluster_heimdall.dump_cluster_results_json(tab2, outfile)

    assert os.path.exists(outfile)

    
def test_json(tab):
    outfile = os.path.join(_install_dir, 'output.cand')
    cluster_heimdall.cluster_data(tab, return_clusterer=False)

    tab2 = cluster_heimdall.get_peak(tab)
    cluster_heimdall.dump_cluster_results_heimdall(tab2, outfile)
    assert os.path.exists(outfile)


def test_plot_dmhist(tab):
    plotting.plot_dm_hist(tab, plot_dir=os.path.join(_install_dir, 'plot_'))


def test_plot_bt(tab):
    plotting.plot_beam_time(tab, plot_dir=os.path.join(_install_dir, 'plot_'))


def test_giantst(tab):
    plotting.plot_giants(tab, plot_dir=os.path.join(_install_dir, 'plot_'))
