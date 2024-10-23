import pytest
import os.path
from astropy import table
from T2 import cluster_heimdall, plotting

_install_dir = os.path.abspath(os.path.dirname(__file__))

@pytest.fixture(scope="module")
def tab():
    candsfile = os.path.join(_install_dir, 'data/T1_output1729695471.csv')
    tab = cluster_heimdall.parse_candsfile(candsfile)
    return tab


@pytest.fixture(scope="module")
def tabs():
    candsfile1 = os.path.join(_install_dir, 'data/T1_output1729695471.csv')
    candsfile2 = os.path.join(_install_dir, 'data/T1_output1729719874.csv')
    tab1 = cluster_heimdall.parse_candsfile(candsfile1)
    tab2 = cluster_heimdall.parse_candsfile(candsfile2)
    return tab1, tab2


def test_parse(tab):
    assert len(tab) == 1000
    assert len(tab[0]) == 11


def test_cluster1(tab):
    cluster_heimdall.cluster_data(tab, return_clusterer=False)
    assert len(tab[0]) == 11


def test_peak(tab):
    cluster_heimdall.cluster_data(tab, return_clusterer=False)
    tab2 = cluster_heimdall.get_peak(tab)
    assert len(tab2) == 1
    assert len(tab2[0]) == 21
    assert tab2['snr'].max() == 17.6562
    # TODO
    # assert cb ==
    # assert cc ==


def test_filter(tab):
    cluster_heimdall.cluster_data(tab, return_clusterer=False)
    tab2 = cluster_heimdall.get_peak(tab)
    tab2 = cluster_heimdall.filter_clustered(tab2, min_snr=1000)  # remove all
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
