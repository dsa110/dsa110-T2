import numpy as np
import socket 
from T2 import cluster_heimdall, plotting

import dsautils.dsa_syslog as dsl
logger = dsl.DsaSyslogger()
logger.subsystem('software')
logger.app('T2')


def parse_socket(host, ports, selectcols=['itime', 'idm', 'ibox', 'ibeam'], outputfile=None, plot_dir=None, trigger=False):
    """ 
    Takes standard MBHeimdall giants socket output and returns full table, classifier inputs and snr tables.
    selectcols will take a subset of the standard MBHeimdall output for cluster. 
    
    host, ports: same with heimdall -coincidencer host:port
    ports can be list of integers.
    selectcol: list of str.  Select columns for clustering. 
    """

    if isinstance(ports, int):
        ports = [ports]

    assert isinstance(ports, list)

    ss = []
    for port in ports:
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM) 
        s.bind((host, port))    # assigns the socket with an address
        s.listen(5)             # accept no. of incoming connections
        ss.append(s)

    while True:
        cls = []
        try:
            for s in ss:
                clientsocket, address = s.accept() # stores the socket details in 2 variables
                logger.info(f"Connection from {address} has been established")
                cls.append(clientsocket)
        except KeyboardInterrupt:
            logger.info("Escaping socket connection")
            break

        # read in heimdall socket output  
        logger.info(f"Reading candsfile from {len(cls)} sockets...")
        candsfile = ''
        gulps = []
        for i, cl in enumerate(cls):
            cf = recvall(cl, 100000000).decode('utf-8')
            cl.close()

            gulp, *lines = cf.split('\n')
            try:
                gulp = int(gulp)
                gulps.append(gulp)
#                print(f"received gulp {gulp} with {len(lines)-1} lines")
            except ValueError:
                print(f'Could not get int from gulp: {gulp}.')
                continue

            if len(lines) > 1:
                if len(lines[0]) > 0:
                    candsfile += '\n'.join(lines)

        print(f'Received gulp_i {gulps}')
        if len(gulps) != len(cls):
            logger.info(f"not all clients are gulping gulp {gulp_i}.")
            print(f"not all clients are gulping gulp {gulp_i}. Skipping...")
            continue
        else:
            #print(f"receiving from all ports")
            pass

        if len(set(gulps)) > 1:
            logger.info(f"not all clients received from same gulp: {set(gulps)}")
            print(f"not all clients received from same gulp: {set(gulps)}.")
            continue

        if candsfile == '\n' or candsfile == '':  # skip empty candsfile
            continue

        try:
            tab = cluster_heimdall.parse_candsfile(candsfile)
            logger.info(f"Table has {len(tab)} rows")
            if len(tab) == 0:
                continue
            cluster_and_plot(tab, set(gulps).pop(), selectcols=selectcols, outputfile=outputfile, plot_dir=plot_dir,
                             trigger=trigger)
        except KeyboardInterrupt:
            logger.info("Escaping parsing and plotting")
            break
        except OverflowError:
            logger.warning("overflowing value. Skipping this gulp...")
            continue


def cluster_and_plot(tab, gulp_i, selectcols=['itime', 'idm', 'ibox', 'ibeam'], outputfile=None, plot_dir=None,
                     trigger=False):
    """ 
    Run clustering and plotting on read data.
    Can optionally save clusters as heimdall candidate table before filtering and json version of buffer trigger.
    """

    # TODO: put these in json config file
    min_dm = 100.0  # smallest dm in filtering
    max_ibox = 20  # largest ibox in filtering
    min_snr = 8.  # smallest snr in filtering
    max_ncl = 10  # largest number of clusters allowed in triggering

    # cluster
    cluster_heimdall.cluster_data(tab, metric='euclidean', allow_single_cluster=True, return_clusterer=False)
    tab2 = cluster_heimdall.get_peak(tab)
    tab3 = cluster_heimdall.filter_clustered(tab2, min_snr=min_snr, min_dm=min_dm, max_ibox=max_ibox)

    col_trigger = np.zeros(len(tab2), dtype=int)
    if outputfile is not None and len(tab3):
        tab4 = cluster_heimdall.dump_cluster_results_json(tab3, outputfile+str(gulp_i)+".json", trigger=trigger, max_ncl=max_ncl)
        if tab4 is not None:
            col_trigger = np.where(tab4 == tab2, 1, 0)  # if trigger, then overload

    # send T2 cluster results to outputfile
    if outputfile is not None and len(tab2):
        tab2['trigger'] = col_trigger
        cluster_heimdall.dump_cluster_results_heimdall(tab2, outputfile+str(gulp_i)+".cand")

#    if plot_dir is not None: 
#         plotting.plot_giants(tab, plot_dir=plot_dir+str(gulp_i)+"_") # plot giants      
#         plotting.plot_clustered(clusterer, clsnr, snrs, data, tab, cols=['itime', 'idm', 'ibox'], plot_dir=plot_dir+str(gulp_i)+"_") # plot cluster results  


def recvall(sock, n):
    """
    helper function to receive all bytes from a socket
    sock: open socket
    n: maximum number of bytes to expect. you can make this ginormous!
    """

    data = bytearray()
    while len(data) < n:
        packet = sock.recv(n - len(data))
        if not packet:
            return data
        data.extend(packet)

    return data
