import pytest
from T2 import cluster_heimdall, triggering
import os.path
import numpy as np

_install_dir = os.path.abspath(os.path.dirname(__file__))
    
def test_on_clustered(source_check=False, catalog=None, beam_model=None, outputfile=os.path.join(_install_dir, 'tmp.dat')):
    candsfile = os.path.join(_install_dir, 'data/giants.cand')

    tab = cluster_heimdall.parse_candsfile(candsfile)
    nbeams_gulp = cluster_heimdall.get_nbeams(tab)
    tab2 = cluster_heimdall.filter_clustered(tab, min_snr=8.0, min_dm=20.0, max_ibox=33.0)
    tab2['gulp_ncl'] = [len(tab)] * len(tab2)

    if len(tab2):
        maxsnr = tab2['snr'].max()
        imaxsnr = np.where(tab2['snr'] == maxsnr)[0][0]
        tab3 = tab2[imaxsnr:imaxsnr+1]

        if source_check is True and catalog is not None:
            coords,snrs = triggering.parse_catalog(catalog)
            tab4 = triggering.check_clustered_sources(tab3, coords, snrs, beam_model=beam_model)

            if len(tab4):
                tab4.write(outputfile, format='ascii.no_header', overwrite=True)

        else:
            tab3.write(outputfile, format='ascii.no_header', overwrite=True)

            
