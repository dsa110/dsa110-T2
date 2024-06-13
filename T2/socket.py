import socket

import numpy as np

from T2 import cluster_heimdall

try:
    from T2 import triggering
except ModuleNotFoundError:
    print("not importing triggering")
import datetime
import time

from astropy.time import Time
from astropy import units
from dsautils import cnf, dsa_store, dsa_syslog
from etcd3.exceptions import ConnectionFailedError
from event import names
import os
import pandas

ds = dsa_store.DsaStore()
logger = dsa_syslog.DsaSyslogger()
logger.subsystem("software")
logger.app("T2")
my_cnf = cnf.Conf(use_etcd=True)
try:
    t2_cnf = my_cnf.get("t2")
except (KeyError, ConnectionFailedError):
    print("Cannot find t2 cnf using etcd. Falling back to hard coded values.")
    logger.warning(
        "Cannot find t2 cnf using etcd. Falling back to hard coded values."
    )
    my_cnf = cnf.Conf(use_etcd=False)
    t2_cnf = my_cnf.get("t2")

from collections import deque

nbeams_queue = deque(maxlen=10)


def parse_socket(
    host,
    ports,
    selectcols=["itime", "idm", "ibox", "ibeam"],
    outroot=None,
    plot_dir=None,
    trigger=False,
    source_catalog=None,
):
    """
    Takes standard MBHeimdall giants socket output and returns full table, classifier inputs and snr tables.
    selectcols will take a subset of the standard MBHeimdall output for cluster.

    host, ports: same with heimdall -coincidencer host:port
    ports can be list of integers.
    selectcol: list of str.  Select columns for clustering.
    source_catalog: path to file containing source catalog for source rejection. default None
    """

    # startup time
    min_timedelt = 10. ## TODO put this in etcd
    prev_trig_time = Time.now()
    
    if isinstance(ports, int):
        ports = [ports]

    assert isinstance(ports, list)

    lastname = names.get_lastname()
    lastname_cleared = lastname

    ss = []

    # pre-calculate beam model and get source catalog
    if source_catalog is not None:
        # model = triggering.get_2Dbeam_model()
        # model = triggering.read_beam_model()
        model = None
        coords, snrs = triggering.parse_catalog(source_catalog)
    else:
        print("No source catalog found. No model generated.")
        model = None
        coords = None
        snrs = None

    logger.info(f"Reading from {len(ports)} sockets...")
    print(f"Reading from {len(ports)} sockets...")
    while True:
        if len(ss) != len(ports):
            for port in ports:
                s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                s.bind((host, port))  # assigns the socket with an address
                s.listen(1)  # accept no. of incoming connections
                ss.append(s)

        ds.put_dict(
            "/mon/service/T2service",
            {"cadence": 60, "time": Time(datetime.datetime.utcnow()).mjd},
        )

        cls = []
        try:
            for s in ss:
                (
                    clientsocket,
                    address,
                ) = s.accept()  # stores the socket details in 2 variables
                cls.append(clientsocket)
        except KeyboardInterrupt:
            print("Escaping socket connection")
            logger.info("Escaping socket connection")
            break

        # read in heimdall socket output
        candsfile = ""
        gulps = []
        for cl in cls:
            cf = recvall(cl, 100000000).decode("utf-8")

            gulp, *lines = cf.split("\n")
            try:
                gulp = int(gulp)
            #                print(f"received gulp {gulp} with {len(lines)-1} lines")
            except ValueError:
                print(
                    f"Could not get int from this read ({gulp}). Skipping this client."
                )
                continue

            gulps.append(gulp)
            cl.close()

            if len(lines) > 1:
                if len(lines[0]) > 0:
                    candsfile += "\n".join(lines)

        print(f"Received gulp_i {gulps}; prev_trig_time {prev_trig_time}")
        if len(gulps) != len(cls):
            print(f"not all clients are gulping gulp {gulps}. Skipping...")
            gulp_status(1)
            continue

        if len(set(gulps)) > 1:
            logger.info(
                f"not all clients received from same gulp: {set(gulps)}. Restarting socket connections."
            )
            print(
                f"not all clients received from same gulp: {set(gulps)}. Restarting socket connections."
            )

            for s in ss:
                s.close()
            ss = []
            for port in ports:
                s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                try:
                    s.bind((host, port))  # assigns the socket with an address
                except OSError:
                    print("socket bind failed.")
                    continue
                s.listen(1)  # accept no. of incoming connections
                ss.append(s)
            gulp_status(2)
            continue
        else:
            gulp = set(gulps).pop()  # get gulp number
            ds.put_dict(
                "/mon/service/T2gulp",
                {"cadence": 60, "time": Time(datetime.datetime.utcnow()).mjd},
            )

        if candsfile == "\n" or candsfile == "":  # skip empty candsfile
            print(f"candsfile is empty. Skipping.")
            logger.info(f"candsfile is empty. Skipping.")

            gulp_status(0)
            continue

        # send flush trigger after min_timedelt (once per candidate)
        if Time.now() - prev_trig_time > min_timedelt*units.s and lastname_cleared != lastname:
            ds.put_dict('/cmd/corr/0', {'cmd': 'trigger', 'val': '0-flush-'})
            lastname_cleared = lastname   # reset to avoid continuous calls
            prev_trig_time = Time.now()  # pass this on to log extra triggers in second latency window

        try:
            tab = cluster_heimdall.parse_candsfile(candsfile)
            lastname, trigtime = cluster_and_plot(
                tab,
                gulp=gulp,
                selectcols=selectcols,
                outroot=outroot,
                plot_dir=plot_dir,
                trigger=trigger,
                lastname=lastname,
                cat=source_catalog,
                beam_model=model,
                coords=coords,
                snrs=snrs,
                prev_trig_time=prev_trig_time
            )
            if trigtime is not None:
                prev_trig_time = trigtime
        except KeyboardInterrupt:
            print("Escaping parsing and plotting")
            logger.info("Escaping parsing and plotting")
            break
        except OverflowError:
            print("overflowing value. Skipping this gulp...")
            logger.warning("overflowing value. Skipping this gulp...")

            print(candsfile)
            gulp_status(3)
            continue
        gulp_status(0)  # success!


def cluster_and_plot(
        tab,
        gulp=None,
        selectcols=["itime", "idm", "ibox", "ibeam"],
        outroot=None,
        plot_dir=None,
        trigger=False,
        lastname=None,
        cat=None,
        beam_model=None,
        coords=None,
        snrs=None,
        prev_trig_time=None
    ):
    """
    Run clustering and plotting on read data.
    Can optionally save clusters as heimdall candidate table before filtering and json version of buffer trigger.
    lastname is name of previously triggered/named candidate.
    cat: path to source catalog (default None)
    beam_model: pre-calculated beam model (default None)
    coords and snrs: from source catalog (default None)
    """

    # TODO: put these in json config file
    min_timedelt = 60. ## TODO put this in etcd
    trigtime = None
    columns = ['snr','if','specnum','mjds','ibox','idm','dm','ibeam','cl','cntc','cntb','trigger']
    
    use_gal_dm = t2_cnf.get("use_gal_dm", None)  # not used currently
    if use_gal_dm:
        dm_mw = ds.get_dict('/mon/array/gal_dm')['gal_dm']
        logger.debug(f'Using DM of {dm_mw}')
    else:
        dm_mw = 0
    min_dm = max(50., dm_mw*0.75)
    max_ibox = t2_cnf.get("max_ibox", None)
    min_snr = t2_cnf.get("min_snr", None)
    min_snr_t2out = t2_cnf.get("min_snr_t2out", None)  # write T2 output cand file above this snr
    max_cntb = t2_cnf.get("max_ctb", None)
    writeT1 = t2_cnf.get("writeT1", False)
    nbeams_max = t2_cnf.get("nbeams_max", None)

    # cluster
    cluster_heimdall.cluster_data(
        tab,
        metric="euclidean",
        allow_single_cluster=True,
        return_clusterer=False,
    )
    if writeT1:
        output_file = outroot + "T1_output" + str(np.floor(time.time()).astype("int")) + ".cand"
        outputted = cluster_heimdall.dump_cluster_results_heimdall(tab,
                                                                   output_file,
                                                                   min_snr_t2out=min_snr_t2out)

    # filter the peaks
    tab2 = cluster_heimdall.get_peak(tab)
    nbeams_gulp = cluster_heimdall.get_nbeams(tab2)
    nbeams_queue.append(nbeams_gulp)
    print(f"nbeams_queue: {nbeams_queue}")
    tab3 = cluster_heimdall.filter_clustered(
        tab2,
        min_dm=min_dm,
        min_snr=min_snr,
        max_ibox=max_ibox,
        max_cntb=max_cntb,
    )

    # trigger decision
    col_trigger = np.zeros(len(tab2), dtype=int)
    if outroot is not None and len(tab3):
        tab4, lastname, trigtime = cluster_heimdall.dump_cluster_results_json(
            tab3,
            trigger=trigger,
            lastname=lastname,
            gulp=gulp,
            cat=cat,
            beam_model=beam_model,
            coords=coords,
            snrs=snrs,
            outroot=outroot,
            nbeams=sum(nbeams_queue),
            nbeams_max=nbeams_max,
            prev_trig_time=prev_trig_time,
            min_timedelt=min_timedelt
        )
        if tab4 is not None and trigger:
            col_trigger = np.where(
                tab4 == tab2, lastname, 0
            )  # if trigger, then overload

    # write T2 clustered/filtered results
    if outroot is not None and len(tab2):
        tab2["trigger"] = col_trigger
        output_file = outroot + "cluster_output" + str(np.floor(time.time()).astype("int")) + ".cand"
        outputted = cluster_heimdall.dump_cluster_results_heimdall(tab2,
                                                                   output_file,
                                                                   min_snr_t2out=min_snr_t2out)

        # aggregate files
        if outputted:
            a = Time.now().mjd
            output_mjd = str(int(a))
            old_mjd = str(int(a)-1)
            fl1 = outroot+old_mjd+".csv"
            fl2 = outroot+output_mjd+".csv"
            ofl = outroot+"cluster_output.csv"

            df0 = pandas.read_csv(output_file, delimiter=' ', names=columns)

            dfs = [df0]
            if os.path.exists(fl1):  # accumulate to yesterday's for rolling 2-day file
                df1 = pandas.read_csv(fl1)
                dfs.append(df1)

            if os.path.exists(fl2):  # accumulate to today's for 1-day file
                df2 = pandas.read_csv(fl2)
                dfs.append(df2)
                dfc2 = pandas.concat( (df0, df2) )
                dfc2.to_csv(fl2, index=False)
            else:
                df0.to_csv(fl2, index=False)

            dfc = pandas.concat(dfs)
            dfc.to_csv(ofl, index=False)

    return lastname, trigtime

def cross_match_peaks(tab_ew, tab_ns, 
                      bin_width_sec=1,
                      dmmax=1000):

    # Convert MJD to seconds starting at zero
    mjd_min = min(tab_ew['mjds'].min(), tab_ns['mjds'].min())
    time_ew_sec = np.array(tab_ew['mjds'][:] - mjd_min) * 86400
    time_ns_sec = np.array(tab_ns['mjds'][:] - mjd_min) * 86400

    max_time_sec = max(time_ew_sec.max(), time_ns_sec.max())
    nbin_time = max_time_sec / bin_width_sec
    nbin_time = int(nbin_time)

    nbin_dm = 256

    # Grid the EW arm candidates
    arr_ew, time_bins, dm_bins = np.histogram2d(time_ew_sec, tab_ew['dm'][:], 
                            bins=(nbin_time, nbin_dm),
                            range=((0, max_time_sec),
                                   (0, dmmax)))
    
    # Grid the NS arm candidates
    arr_ns, time_bins, dm_bins = np.histogram2d(time_ns_sec, tab_ns['dm'][:], 
                            bins=(nbin_time, nbin_dm),
                            range=((0, max_time_sec),
                                   (0, dmmax)))

    # cross match
    cross_match_arr = arr_ew * arr_ns

    return cross_match_arr, time_bins, dm_bins

def cluster_twoarms(
        tab_ew,
        tab_ns,
        gulp=None,
        outroot=None,
        trigger=False,
        lastname=None,
        cat=None,
        beam_model=None,
        coords=None,
        snrs=None,
        prev_trig_time=None
        cluster_selection_epsilon=10,
    ):
    """
    Run clustering and plotting on read data.
    Can optionally save clusters as heimdall candidate table before filtering and json version of buffer trigger.
    lastname is name of previously triggered/named candidate.
    cat: path to source catalog (default None)
    beam_model: pre-calculated beam model (default None)
    coords and snrs: from source catalog (default None)
    """

    # TODO: put these in json config file
    min_timedelt = 60. ## TODO put this in etcd
    trigtime = None
    columns = ['snr','if','specnum','mjds','ibox','idm','dm','ibeam','cl','cntc','cntb','trigger']
    
    use_gal_dm = t2_cnf.get("use_gal_dm", None)  # not used currently
    if use_gal_dm:
        dm_mw = ds.get_dict('/mon/array/gal_dm')['gal_dm']
        logger.debug(f'Using DM of {dm_mw}')
    else:
        dm_mw = 0
    min_dm = max(50., dm_mw*0.75)
    max_ibox = t2_cnf.get("max_ibox", None)
    min_snr = t2_cnf.get("min_snr", None)
    min_snr_t2out = t2_cnf.get("min_snr_t2out", None)  # write T2 output cand file above this snr
    max_cntb = t2_cnf.get("max_ctb", None)
    writeT1 = t2_cnf.get("writeT1", False)
    nbeams_max = t2_cnf.get("nbeams_max", None)

    # cluster east-west arm, ignore beam number
    cluster_data(
                tab_ew,
                selectcols=["itime", "idm", "ibox"],
                min_cluster_size=2,
                min_samples=5,
                metric="euclidean",
                return_clusterer=False,
                allow_single_cluster=True,
                cluster_selection_epsilon=cluster_selection_epsilon,
                )
    
    # cluster north-south arm, ignore beam number
    cluster_data(
                tab_ns,
                selectcols=["itime", "idm", "ibox"],
                min_cluster_size=2,
                min_samples=5,
                metric="euclidean",
                return_clusterer=False,
                allow_single_cluster=True,
                cluster_selection_epsilon=cluster_selection_epsilon,
                )

    if writeT1:
        output_file_ew = outroot + "T1_output_ew" + str(np.floor(time.time()).astype("int")) + ".cand"
        outputted = cluster_heimdall.dump_cluster_results_heimdall(tab_ew,
                                                                   output_file_ew,
                                                                   min_snr_t2out=min_snr_t2out)
        
        output_file_ns = outroot + "T1_output_ns" + str(np.floor(time.time()).astype("int")) + ".cand"
        outputted = cluster_heimdall.dump_cluster_results_heimdall(tab_ns,
                                                                   output_file_ns,
                                                                   min_snr_t2out=min_snr_t2out)

    # filter the peaks
    tab2_ew = cluster_heimdall.get_peak(tab_ew)
    tab2_ns = cluster_heimdall.get_peak(tab_ns)
    
    cross_match_peaks(tab2_ew, tab2_ns)


#    nbeams_gulp = cluster_heimdall.get_nbeams(tab2_ew)
#    nbeams_gulp = cluster_heimdall.get_nbeams(tab2_ew)
    nbeams_queue.append(nbeams_gulp)
    print(f"nbeams_queue: {nbeams_queue}")
    tab3 = cluster_heimdall.filter_clustered(
        tab2,
        min_dm=min_dm,
        min_snr=min_snr,
        max_ibox=max_ibox,
        max_cntb=max_cntb,
    )

    # trigger decision
    col_trigger = np.zeros(len(tab2), dtype=int)
    if outroot is not None and len(tab3):
        tab4, lastname, trigtime = cluster_heimdall.dump_cluster_results_json(
            tab3,
            trigger=trigger,
            lastname=lastname,
            gulp=gulp,
            cat=cat,
            beam_model=beam_model,
            coords=coords,
            snrs=snrs,
            outroot=outroot,
            nbeams=sum(nbeams_queue),
            nbeams_max=nbeams_max,
            prev_trig_time=prev_trig_time,
            min_timedelt=min_timedelt
        )
        if tab4 is not None and trigger:
            col_trigger = np.where(
                tab4 == tab2, lastname, 0
            )  # if trigger, then overload

    # write T2 clustered/filtered results
    if outroot is not None and len(tab2):
        tab2["trigger"] = col_trigger
        output_file = outroot + "cluster_output" + str(np.floor(time.time()).astype("int")) + ".cand"
        outputted = cluster_heimdall.dump_cluster_results_heimdall(tab2,
                                                                   output_file,
                                                                   min_snr_t2out=min_snr_t2out)

        # aggregate files
        if outputted:
            a = Time.now().mjd
            output_mjd = str(int(a))
            old_mjd = str(int(a)-1)
            fl1 = outroot+old_mjd+".csv"
            fl2 = outroot+output_mjd+".csv"
            ofl = outroot+"cluster_output.csv"

            df0 = pandas.read_csv(output_file, delimiter=' ', names=columns)

            dfs = [df0]
            if os.path.exists(fl1):  # accumulate to yesterday's for rolling 2-day file
                df1 = pandas.read_csv(fl1)
                dfs.append(df1)

            if os.path.exists(fl2):  # accumulate to today's for 1-day file
                df2 = pandas.read_csv(fl2)
                dfs.append(df2)
                dfc2 = pandas.concat( (df0, df2) )
                dfc2.to_csv(fl2, index=False)
            else:
                df0.to_csv(fl2, index=False)

            dfc = pandas.concat(dfs)
            dfc.to_csv(ofl, index=False)

    return lastname, trigtime

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


def gulp_status(status):
    """Set etcd key to track gulp status.
    0 means good, non-zero means some kind of failure for a gulp.
    1 means not all clients are gulping
    2 means different gulps received, so restarting clients
    3 means overflow error during parsing of table.
    t2_num is the process number running T2. Only one for now.
    """

    t2_num = 1

    ds.put_dict(
        f"/mon/T2/{t2_num}",
        {
            "gulp_status": int(status),
            "t2_num": t2_num,
            "time": Time(datetime.datetime.utcnow()).mjd,
        },
    )
