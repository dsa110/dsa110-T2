import os

import numpy as np
import glob 

import dsautils.dsa_store as ds
import time
from astropy.time import Time

d = ds.DsaStore()
fmt_out = '%5.9f  %d  %0.2f %0.1f %0.3f %0.2f %s\n'
fnout = '/home/ubuntu/injection_list.txt'

if not os.path.exists(fnout):
    f = open(fnout,'w+')
    f.write('# MJD   Beam   DM    SNR   Width_fwhm   spec_ind  FRBno\n')
    f.close()

params = np.genfromtxt('/home/ubuntu/simulated_frb_params.txt')
flist = glob.glob('/home/ubuntu/data/test_inj*.dat')
flist = ['/home/ubuntu/data/test_inj_0010.dat','/home/ubuntu/data/test_inj_0008.dat']
nfrb = len(flist)
nfrb=200

print(f"Starting loop over {nfrb} injections with names {flist[0]} and {flist[1]} at MJD: {Time.now().mjd}")

for zz in range(100):
    for kk in range(17,21):
        for ii in range(nfrb):
            f = open(fnout,'a')
            subbeam = (2*ii+1) % 64
            beam = 64*(kk-17)+subbeam
            print("Injecting into beam %d"%beam)
            fn = flist[int(ii%2)]
            frbno = fn.split('_')[-1][:4]
            ind = np.where(params[:,-1]==float(frbno))[0]
            DM, SNR, Width_fwhm, spec_ind = params[ind][0][0],params[ind][0][1],params[ind][0][2],params[ind][0][3]
            print("pushing injection to command to etcd")
            d.put_dict('/cmd/corr/%d'%kk,{'cmd':'inject','val':'%d-%s-'%(subbeam,fn)})
            imjd = Time.now().mjd
#            inj_dict = {'mjd':imjd,
#                        'ibeam':beam,
#                        'dm':DM,
#                        'snr':SNR,
#                        'width':Width_fwhm,
#                        'spec_ind':spec_ind,
#                        'frbno':frbno}
#            d.put_dict('/mon/corr/injection',inj_dict)
            print("writing parameters to disk")
            f.write(fmt_out % (imjd, beam, DM, SNR, Width_fwhm, spec_ind, frbno))
            f.close()
            time.sleep(2700)

f.close()
