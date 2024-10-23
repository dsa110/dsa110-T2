import pytest
from T2 import cluster_heimdall, triggering
import os.path

_install_dir = os.path.abspath(os.path.dirname(__file__))
    
def test_on_clustered(source_check=False, catalog=None, beam_model=None, max_ncl=10, outputfile=os.path.join(_install_dir, 'tmp.dat')):
    candsfile = os.path.join(_install_dir, 'data/T1_output1729719680.csv')

    tab = cluster_heimdall.parse_candsfile(candsfile)
    nbeams_gulp = cluster_heimdall.get_nbeams(tab)
    tab2 = cluster_heimdall.get_peak(tab)
    tab3 = cluster_heimdall.filter_clustered(tab2, min_snr=8.0, min_dm=20.0, max_ibox=33.0, max_ncl=30)
    tab3['gulp_ncl'] = [len(tab)] * len(tab3)

    if len(tab2) and len(tab2)<max_ncl:
        maxsnr = tab2['snr'].max()
        imaxsnr = np.where(tab2['snr'] == maxsnr)[0][0]
        tab3 = tab2[imaxsnr:imaxsnr+1]

        if source_check is True and catalog is not None:
            coords,snrs = triggering.parse_catalog(catalog)
            tab4 = triggering.check_clustered_sources(tab3, coords, snrs, beam_model=beam_model)

            if len(tab4):
                tab4.write(outputfile, format='ascii.no_header')

        else:
            tab3.write(outputfile, format='ascii.no_header')

            
