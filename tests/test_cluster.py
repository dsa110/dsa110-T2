import pytest
import os.path
from T2 import cluster_heimdall, plotting

_install_dir = os.path.abspath(os.path.dirname(__file__))

@pytest.fixture(scope="module")
def datasnrs():
    candsfile = os.path.join(_install_dir, 'data/giants.cand')
    tab, data, snrs = cluster_heimdall.parse_candsfile(candsfile)
    return tab, data, snrs


def test_parse(datasnrs):
    tab, data, snrs = datasnrs
    assert len(data)
    assert len(data) == len(snrs)


def test_cluster1(datasnrs):
    tab, data, snrs = datasnrs
    datal = cluster_heimdall.cluster_data(data, return_clusterer=False)

    assert len(datal) == len(data)
    assert len(datal[0]) == 7


def test_peak(datasnrs):
    tab, data, snrs = datasnrs
    datal = cluster_heimdall.cluster_data(data, return_clusterer=False)

    clsnr = cluster_heimdall.get_peak(datal, snrs)
    assert len(clsnr) == 1

    i, snr, cb, cc = clsnr[0]

    assert i == 1380
    assert snr == 117.613
    # TODO
    # assert cb ==
    # assert cc ==


def test_filter(datasnrs):
    tab, data, snrs = datasnrs
    datal = cluster_heimdall.cluster_data(data, return_clusterer=False)

    clsnr = cluster_heimdall.get_peak(datal, snrs)
    assert len(clsnr) == 1
    clsnr = cluster_heimdall.filter_clustered(clsnr, min_snr=1000)
    assert len(clsnr) == 0


def test_json(datasnrs):
    outfile = os.path.join(_install_dir, 'trigger.json')
    tab, data, snrs = datasnrs
    datal = cluster_heimdall.cluster_data(data, return_clusterer=False)

    clsnr = cluster_heimdall.get_peak(datal, snrs)
    cluster_heimdall.dump_cluster_results_json(tab, clsnr, outfile)

    assert os.path.exists(outfile)

    
def test_json(datasnrs):
    outfile = os.path.join(_install_dir, 'outout.cand')
    tab, data, snrs = datasnrs
    datal = cluster_heimdall.cluster_data(data, return_clusterer=False)

    clsnr = cluster_heimdall.get_peak(datal, snrs)
    cluster_heimdall.dump_cluster_results_heimdall(tab, clsnr, outputfile)
    assert os.path.exists(outfile)


def test_plot_dmhist(datasnrs):
    tab, data, snrs = datasnrs
    plotting.plot_dm_hist(tab)


def test_plot_bt(datasnrs):
    tab, data, snrs = datasnrs
    plotting.plot_beam_time(tab, plot_dir=_install_dir + 'plot_')


def test_giantst(datasnrs):
    tab, data, snrs = datasnrs
    plotting.plot_giants(tab, plot_dir=_install_dir + 'plot_')
