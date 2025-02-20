import os

import numpy as np
import glob 

import dsautils.dsa_store as ds
import time
from astropy.time import Time
import slack_sdk as slack

# set up slack client
slack_file = '{0}/.config/slack_api'.format(os.path.expanduser("~"))
if not os.path.exists(slack_file):
    raise RuntimeError("Could not find file with slack api token at {0}".format(slack_file))

with open(slack_file) as sf_handler:
    slack_token = sf_handler.read()
    slack_client = slack.WebClient(token=slack_token)

d = ds.DsaStore()
fmt_out = '%5.9f  %d  %0.2f %0.1f %0.3f %0.2f %s\n'

# This file has parameters for FRBs that have been injected
fnout = '/home/ubuntu/data/injections/injection_list.txt'

if not os.path.exists(fnout):
    f = open(fnout,'w+')
    f.write('# MJD   Beam   DM    SNR   Width_fwhm   spec_ind  FRBno\n')
    f.close()

# This file has the parameters for simulated FRBs
scfac = 0.2
params = np.genfromtxt('/home/ubuntu/simulated_frb_params.txt')
print(params)
flist = ['/home/ubuntu/data/burst_0.inject','/home/ubuntu/data/burst_1.inject','/home/ubuntu/data/burst_2.inject','/home/ubuntu/data/burst_3.inject','/home/ubuntu/data/burst_4.inject']
nfrb = len(flist)

for kk in [17,18]:
    for ii in np.arange(128):
            
        f = open(fnout,'a')
        subbeam = ii
        beam = 128*(kk-17)+subbeam
        print("Injecting into beam %d"%beam)
        fn = flist[int(ii%5)]
        frbno = fn.split('_')[-1].split('.')[0]
        ind = int(ii%5)#np.where(params[:,-1]==float(frbno))[0]
        DM, SNR, Width_fwhm, spec_ind = params[ind][0],params[ind][1],params[ind][2],params[ind][3]
        print("pushing injection to command to etcd")
        slack_client.chat_postMessage(channel='candidates', text=f'Sending injection to beam {subbeam} with DM={DM} and SNR={SNR}')
        d.put_dict('/cmd/corr/%d'%kk,{'cmd':'inject','val':'%d-%s-%f-'%(subbeam,fn,scfac)})
        d.put_dict('/cmd/corr/%d'%(kk+2),{'cmd':'inject','val':'%d-%s-%f-'%(subbeam,fn,scfac*35./47.)})
        imjd = Time.now().mjd
        print("writing parameters to disk")
        f.write(fmt_out % (imjd, beam, DM, SNR, Width_fwhm, spec_ind, frbno))
        f.close()
        print("Waiting to inject...")            
        time.sleep(1600)

f.close()
