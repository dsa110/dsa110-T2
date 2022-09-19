#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 11:20:32 2020

@author: liamconnor
"""

import socket 
from astropy.io import ascii 
import numpy as np
import  pandas as pd

from T2 import cluster_heimdall
from astropy.io import ascii

HOST = "127.0.0.1"
PORT = 12346

gulpsize=16384

s=socket.socket(socket.AF_INET, socket.SOCK_STREAM)
#s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
# Create a UDP socket
s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
# Bind the socket to the port
server_address = (HOST, PORT)
#s.close()
s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
s.bind(server_address)
print("Do Ctrl+c to exit the program !!")

n = 0
x=''

col_heimdall = ['snr', 'if', 'itime', 'mjds', 'ibox', 'idm', 'dm', 'ibeam']
candsfile = ''
gulp = 0

min_dm = 50
max_ibox = 64
min_snr = 8.0
min_snr_t2out = 10.0
max_ncl = np.inf
max_cntb = np.inf
max_cntb0 = np.inf
target_params = (50., 100., 20.)  # Galactic bursts

while True:
#    print("####### Server is listening #######")
    data, address = s.recvfrom(4096)
    tab = data.decode('utf-8')
    itime = int(tab.split('\t')[2])
    iff = int(tab.split('\t')[1])
    candsfile += tab
    if itime//gulpsize != gulp:
        print("GULP", gulp)
        TAB = ascii.read(candsfile, names=col_heimdall,
                     guess=True, fast_reader=False,
                     format='no_header')
        print(len(TAB))
        cluster_heimdall.cluster_data(TAB, metric='euclidean', allow_single_cluster=True, return_clusterer=False)
        print(len(TAB))
        tab2 = cluster_heimdall.get_peak(TAB)
        print(len(tab2))
        tab3 = cluster_heimdall.filter_clustered(tab2, min_snr=min_snr, min_dm=min_dm, max_ibox=max_ibox, max_cntb=max_cntb,
                                             max_cntb0=max_cntb0, max_ncl=max_ncl, target_params=target_params)

        print(len(tab3))
        itimes = tab3['itime']
        maxsnr = tab3['snr'].max()
        imaxsnr = np.where(tab3['snr'] == maxsnr)[0][0]        
        itime_imax = str(itimes[imaxsnr])
        mjd = tab3['mjds'][imaxsnr]
        
        gulp = itime//gulpsize
        print("")
        print(tab3[imaxsnr])
        print("")
        trigger = False
        lastname = 'blegh'
        cat = None
        beam_model = None
        coords = None
        snrs = None
        outroot = './'
        nbeams_queue = 0
        prev_trig_time = None
        min_timedelt = 60.
        print("Starting dump step")
        X = cluster_heimdall.dump_cluster_results_json(
                                                tab3,
                                                trigger=trigger,
                                                lastname=lastname,
                                                cat=cat,
                                                beam_model=beam_model,
                                                coords=coords,
                                                snrs=snrs,
                                                outroot=outroot,
                                                nbeams=1,
                                                frac_wide=0.0,
                                            )
        print("beep", X)
        
exit()
