import socket
import numpy as np

from T2 import cluster_heimdall
import func_timeout as ft
try:
    from T2 import triggering
except ModuleNotFoundError:
    print("not importing triggering")
import datetime
import time

from astropy.time import Time
from astropy import units
try:
    from dsautils import cnf, dsa_store, dsa_syslog
    ds = dsa_store.DsaStore()
    logger = dsa_syslog.DsaSyslogger()
    logger.subsystem("software")
    logger.app("T2")
    my_cnf = cnf.Conf(use_etcd=True)
except:
    ds = None
    import logging
    logger = logging.getLogger()
    
from etcd3.exceptions import ConnectionFailedError
from event import names
import os
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import pandas
from astropy.table import Table

try:
    t2_cnf = my_cnf.get("t2")
except (KeyError, ConnectionFailedError):
    print("Cannot find t2 cnf using etcd. Falling back to hard coded values.")
    logger.warning(
        "Cannot find t2 cnf using etcd. Falling back to hard coded values."
    )
    my_cnf = cnf.Conf(use_etcd=False)
    t2_cnf = my_cnf.get("t2")

#Auditor import
from T2.audit import Auditor  
INJECTION_FILE = "/operations/T2/injection_audit_results/injections_for_audit.txt"  

from collections import deque

nbeams_queue = deque(maxlen=10)


def parse_socket(
    host,
    ports,
    selectcols=["itime", "idm", "ibox"],
    outroot=None,
    plot_dir=None,
    trigger=False,
    source_catalog=None,
    #audit injection parameters
    audit_injections=False,
    audit_dump_json=False,
    audit_dir="/operations/T2/injection_audit_results/",
    dump_t1=False,
    t1_dump_dir="/operations/T2/T1_DUMPS_VISHNU_TEST/",
):
    """
    Takes standard MBHeimdall giants socket output and returns full table, classifier inputs and snr tables.
    selectcols will take a subset of the standard MBHeimdall output for cluster.

    host, ports: same with heimdall -coincidencer host:port
    ports can be list of integers.
    selectcol: list of str.  Select columns for clustering.
    source_catalog: path to file containing source catalog for source rejection. default None
    audit_injections: bool, whether to enable injection auditing (default False)
    audit_dump_json: bool, whether to dump per-injection JSON files (default False)
    audit_dir: str, directory to store audit results (default "/operations/T2/injection_audit_results/")

    """

    # startup time
    min_timedelt = 60. ## TODO put this in etcd
    prev_trig_time = Time.now()
    last_real_trigger_time = None
    last_real_trigger_name = None
    last_flushed_for_name = None
    
    # count of output - separate from gulps
    globct = 0

    if isinstance(ports, int):
        ports = [ports]

    assert isinstance(ports, list)

    lastname = names.get_lastname()
    lastname_cleared = lastname

    ss = []
    cls = []
    
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

    # Audit Injection Setup
    audit_enabled = bool(audit_injections)
    auditor = None
    injection_list_csv = None  # canonical injection list file path

    if audit_enabled:
        # Window parameters: try etcd -> t2_cnf -> hard defaults
        try:
            time_window_s = int(ds.get_dict('/cnf/t2').get('audit_time_window_s', 300))
        except Exception:
            time_window_s = int(t2_cnf.get('audit_time_window_s', 300))

        try:
            dm_window = float(ds.get_dict('/cnf/t2').get('audit_dm_window', 20.0))
        except Exception:
            dm_window = float(t2_cnf.get('audit_dm_window', 20.0))

        try:
            beam_window = int(ds.get_dict('/cnf/t2').get('audit_beam_window', 2))
        except Exception:
            beam_window = int(t2_cnf.get('audit_beam_window', 2))

        # Ensure audit_dir exists
        try:
            os.makedirs(audit_dir, exist_ok=True)
        except Exception as _e:
            logger.warning(f"Could not create audit_dir {audit_dir}: {_e}")

        # Canonical injection list lives inside audit_dir
        injection_list_csv = os.path.join(audit_dir, "injections.csv")

        # Initialize the file if missing, with a commented header
        if not os.path.exists(injection_list_csv):
            try:
                with open(injection_list_csv, "w") as f:
                    f.write("# MJD   Beam   DM    SNR   Width_fwhm   spec_ind   FRBno\n")
            except Exception as _e:
                logger.warning(f"Could not initialize injection list {injection_list_csv}: {_e}")

        # Instantiate Auditor with chosen persistence mode
        from T2.audit import Auditor
        auditor = Auditor(
            audit_dir=audit_dir,
            time_window_s=time_window_s,
            dm_window=dm_window,
            beam_window=beam_window,
            persist_json=bool(audit_dump_json),
        )
        
        auditor.attach_injection_source(INJECTION_FILE)
        added0 = auditor.refresh_from_source(force=True)
        if added0:
            print(f"[AUDIT] initial load: applied {added0} injections from {INJECTION_FILE}")
    # end of audit setup

    logger.info(f"Reading from {len(ports)} sockets...")
    print(f"Reading from {len(ports)} sockets...")

    pool = ThreadPoolExecutor(max_workers=10)
    futures = {}
    trigtime = None
    triggered_ret = False
    trigtime_ret = None
    lastname_ret = None
    while True:

        if len(futures):
            lastname_ret, trigtime_ret, triggered_ret, futures = manage_futures(lastname, trigtime, futures)

        if triggered_ret:
            last_real_trigger_time = trigtime_ret
            last_real_trigger_name = lastname_ret

        if trigtime is not None:
            prev_trig_time = trigtime
            

        # set up socket connections if needed
        if len(ss) != len(ports):
            for port in ports:
                s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                s.bind((host, port))  # assigns the socket with an address
                s.listen(1)  # accept no. of incoming connections
                ss.append(s)
                #print(f"Appended socket for port {port} on host {host}")

        # accept new socket connections
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


                
        ds.put_dict(
            "/mon/service/T2service",
            {"cadence": 60, "time": Time(datetime.datetime.utcnow()).mjd},
        )

        # read in heimdall socket output
        candsfile = ""
        gulps = []
        for cl in cls:

            cf = recvall(cl, 1000000000).decode("utf-8")                
            gulp, *lines = cf.split("\n")
            #print(cl,gulp)

            try:
                gulp = int(gulp)                    
            #                print(f"received gulp {gulp} with {len(lines)-1} lines")
            except ValueError:
                print(
                    f"Could not get int from this read ({gulp}). Skipping this client."
                )
                continue
                                        
            cl.close()
            gulps.append(gulp)                
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

            for cl in cls:
                cl.close()
            time.sleep(0.1)
            cls = []
            for s in ss:
                s.close()
            time.sleep(0.1)
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
        
        # Flush exactly once after a REAL trigger, separated by min_timedelt
        if (last_real_trigger_time is not None
            and (Time.now() - last_real_trigger_time).to_value('sec') >= min_timedelt
            and last_flushed_for_name != last_real_trigger_name):
            print("Sending flush trigger")
            logger.info("Sending flush trigger")
            ds.put_dict('/cmd/corr/0', {'cmd': 'trigger', 'val': '0-flush-'})
            last_flushed_for_name = last_real_trigger_name
        else:
            dt = (Time.now() - last_real_trigger_time).to_value('sec') if last_real_trigger_time is not None else None
            print(f"Cannot send flush: dt_since_real_trigger={dt} min_timedelt={min_timedelt}")
            logger.info(f"Cannot send flush: dt_since_real_trigger={dt} min_timedelt={min_timedelt}")
            print(f"last_flushed_for_name={last_flushed_for_name} last_real_trigger_name={last_real_trigger_name}")
            logger.info(f"last_flushed_for_name={last_flushed_for_name} last_real_trigger_name={last_real_trigger_name}")


        
        if candsfile == "\n" or candsfile == "":  # skip empty candsfile
            print(f"candsfile is empty. Skipping.")
            logger.info(f"candsfile is empty. Skipping.")
            # Advance the rolling window for this gulp
            nbeams_queue.append(0)
            if audit_enabled and auditor is not None:
                auditor.update_from_tab(
                    host=host,
                    gulp=gulp,
                    tab=Table(), #This is an empty table passed for audit update
                    nbeams_queue_snapshot=list(nbeams_queue),
                    nbeams_queue_sum=sum(nbeams_queue) if len(nbeams_queue) else 0,
                    runtime_thresholds={},
                    prev_trig_time=prev_trig_time,
                )

            gulp_status(0)
            continue


        tab = cluster_heimdall.parse_candsfile(candsfile)
        try:
            added = auditor.refresh_from_source(force=False)
            if added:
                print(f"[AUDIT] reload: applied {added} injection rows from {INJECTION_FILE}")
        except Exception as e:
            logger.warning(f"[AUDIT] refresh_from_source failed: {e}")

        #injection auditor changes...
        if audit_enabled and auditor is not None:
            try:
                # Early-stage audit update (G1/G2)
                auditor.update_from_tab(
                    host=host,
                    gulp=gulp,
                    tab=tab,
                    nbeams_queue_snapshot=list(nbeams_queue),
                    nbeams_queue_sum=sum(nbeams_queue) if len(nbeams_queue) else 0,
                    runtime_thresholds={},  # thresholds filled at finalize stage
                    prev_trig_time=prev_trig_time,
                )
            except Exception as e:
                logger.warning(f"[AUDIT] update_from_tab failed: {e}")

        # to handle too many futures
        if len(futures)>1:
            print(f"Waiting for >=2 futures to finish -- skipping {gulps}")
        else:
            now = Time.now()
            key = f'{gulp}-{globct}-{now}'
            future = pool.submit(cluster_and_plot, tab, gulp=gulp, selectcols=selectcols,
                                 outroot=outroot, plot_dir=plot_dir, trigger=trigger, lastname=lastname,
                                 cat=source_catalog, beam_model=model, coords=coords, snrs=snrs,
                                 prev_trig_time=last_real_trigger_time, auditor=auditor, audit_enabled=audit_enabled, dump_t1=dump_t1, t1_dump_dir=t1_dump_dir)
            globct += 1
            futures[key] = future
            print(f'Processing {len(futures)} gulps')

            # try:
            #     lastname, trigtime, futures = manage_futures(lastname, trigtime, futures)  # returns latest result from iteration over futures
            # except:
            #     print('Caught error in manage_futures. Closing sockets.')
            #     for cl in cls:
            #         cl.close()

            if trigtime is not None:
                prev_trig_time = trigtime


def manage_futures(old_lastname, old_trigtime, futures):
    """ Take list of cluster_and_plot futures and handle the output.
    Currently returns only one (lastname, trigtime) tuple for all futures that are done.
    Small chance that lastname or prev_trig_time will not be updated correctly.
    """

    done = []
    triggered = False
    ret_lastname = old_lastname
    ret_trigtime = None
    for k, future in futures.items():
        if future.done():
            done.append(k)
            try:
                lastname,trigtime,did_trigger = future.result()
                if trigtime is not None:
                    gulp_status(0)  # success!
                if lastname is not None:
                    ret_lastname = lastname
                if trigtime is not None:
                    ret_trigtime = trigtime
                triggered = triggered or bool(did_trigger)
            except KeyboardInterrupt:
                print("Escaping parsing and plotting")
                logger.info("Escaping parsing and plotting")
            except OverflowError:
                print("overflowing value. Skipping this gulp...")
                logger.warning("overflowing value. Skipping this gulp...")
                gulp_status(3)

    if len(done):
        for k in done:
            _ = futures.pop(k)
        print(f'{len(done)} gulp future(s) completed')

    return ret_lastname, ret_trigtime, triggered, futures


def cluster_and_plot(tab, gulp=None, selectcols=["itime", "idm", "ibox"],
                     outroot=None, plot_dir=None, trigger=False, lastname=None,
                     max_ncl=None, cat=None, beam_model=None, coords=None,
                     snrs=None, prev_trig_time=None, auditor=None, audit_enabled=False, dump_t1=False, t1_dump_dir=None):
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
    triggered = False
    #injection auditing related
    output_file = ""
    candname_for_audit = None
    columns = ['snr','if','specnum','mjds','ibox','idm','dm','ibeam','cl','cntc','cntb','snrs0','beams0','snrs1','beams1','snrs2','beams2','snrs3','beams3','snrs4','beams4','snrs5','beams5','snrs6','beams6','snrs7','beams7','snrs8','beams8','snrs9','beams9','trigger']
    
    # obtain this from etcd
    # TODO: try a timeout exception
    try:
        max_ibox = ds.get_dict('/cnf/t2')["max_ibox"]  # largest ibox in filtering
    except:
        max_ibox = t2_cnf["max_ibox"]
    
    try:
        min_snr = ds.get_dict('/cnf/t2')["min_snr"]  # smallest snr in filtering
    except:
        min_snr = t2_cnf["min_snr"]

    try:
        min_snr_wide = ds.get_dict('/cnf/t2')["min_snr_wide"]  # smallest snr in filtering
    except:
        min_snr_wide = t2_cnf["min_snr_wide"]

    # adjust min dm according to t2 cnf in etcd
    try:
        use_gal_dm = ds.get_dict('/cnf/t2')["use_gal_dm"]
    except:
        use_gal_dm = 1

    if use_gal_dm == 0:
        min_dm = 20.
    else:
        # Take min DM to be either 0.75 times MW DM or 50., whatever
        # is higher.
        try:
            dm_mw = ds.get_dict('/mon/array/gal_dm')['gal_dm']
            min_dm = max(50., dm_mw*0.75)
        except:
            min_dm = 50.
    
    wide_ibox = t2_cnf["wide_ibox"]  # min ibox for wide snr thresholding
    min_snr_t2out = t2_cnf["min_snr_t2out"]  # write T2 output cand file above this snr
    if max_ncl is None:
        max_ncl = t2_cnf["max_ncl"]  # largest number of clusters allowed in triggering
    max_cntb0 = t2_cnf["max_ctb0"]
    max_cntb = t2_cnf["max_ctb"]
    #target_params = (50.0, 100.0, 20.0)  # Galactic bursts
    target_params = None

    # --- Required for injection auditing ---
    #Vishnu: These two variables (min_snr_1arm, max_nbeams) are currently hardcoded but I think they eventually belong to t2_cnf or etcd. Adding this assuming it will be done later.
    try:
        min_snr_1arm = ds.get_dict('/cnf/t2').get("min_snr_1arm", 10)
    except Exception:
        min_snr_1arm = t2_cnf.get("min_snr_1arm", 10)

    try:
        max_nbeams_allowed = int(ds.get_dict('/cnf/t2').get("max_nbeams", 40))
    except Exception:
        max_nbeams_allowed = int(t2_cnf.get("max_nbeams", 40))

    
    thresholds_snapshot = {
        "min_snr": float(min_snr),
        "min_snr_wide": float(min_snr_wide),
        "min_snr_1arm": float(min_snr_1arm),
        "wide_ibox": int(wide_ibox),
        "max_ibox": int(max_ibox),
        "max_cntb": int(max_cntb),
        "max_cntb0": int(max_cntb0),
        "max_ncl": int(max_ncl) if max_ncl is not None else "",
        "min_dm": float(min_dm),
        "max_nbeams": int(max_nbeams_allowed),
        "min_timedelt": float(60.0),  # matches local default below; but should really belong to etcd/t2_cnf
    }

    #ind = np.where(tab["ibox"]<32)[0]
    #tab = tab[ind]

    # how many points
    print(f"cluster_and_plot: have {len(tab)} inputs")
    logger.info(f"cluster_and_plot: have {len(tab)} inputs")

    #Optional Dumping of T1 candidates to disk.
    if dump_t1 and len(tab):
        try:
            print("Dumping T1 candidates to disk...")
            cluster_heimdall.dump_T1_csv(tab, gulp=gulp, dump_dir=t1_dump_dir)
        except Exception as _e:
            logger.warning(f"T1 dump failed: {_e}")

    
    # raise SNR threshold in case of bright events
    #max_snr = tab["snr"].max()
    #snrthresh = max_snr-10.
    #good = tab["snr"] > snrthresh
    #tab = tab[good]
    
    # how many points
    #print(f"cluster_and_plot: have {len(tab)} inputs ABOVE {snrthresh}")
    #logger.info(f"cluster_and_plot: have {len(tab)} inputs ABOVE {snrthresh}")

    tab_pre_flag = tab.copy()
    # flag beams
    tab_after_flag = cluster_heimdall.flag_beams(tab)
    
    # cluster
    cluster_heimdall.cluster_data(
        tab_after_flag,
        #metric="euclidean",
        allow_single_cluster=True,
        return_clusterer=False,
    )
    tab2 = cluster_heimdall.get_peak(tab_after_flag)
    nbeams_gulp = cluster_heimdall.get_nbeams(tab2, threshold=min_snr)
    nbeams_queue.append(nbeams_gulp)

    # Liam edit to preserve real FRBs during RFI storm:
    # if nbeam > 100 and frac_wide < 0.8: do not discard
    #maxsnr = tab["snr"].max()
    #imaxsnr = np.where(tab["snr"] == maxsnr)[0][0]
    #cl_max = tab["cl"][imaxsnr]
    #frac_wide = np.sum(tab["ibox"][tab["cl"] == cl_max] >= 32) / float(
    #    len(tab["ibox"][tab["cl"] == cl_max])
    #)

    #if len(tab["ibox"][tab["cl"] == cl_max]) == 1:
#    frac_wide = 0.0

    # Width filter for false positives
    #ibox64_filter = False
    #if len(tab2):
    #    ibox64_cnt = np.sum(tab2["ibox"] == 64) / float(len(tab2["ibox"]))
#        print("here", ibox64_cnt, tab2["ibox"])
    #    if ibox64_cnt > 0.85 and len(tab2["ibox"]) > 15:
    #        ibox64_filter = True

    # Done
    tab3 = cluster_heimdall.filter_clustered(
        tab2,
        min_dm=min_dm,
        min_snr=min_snr,
        min_snr_wide=min_snr_wide,
        wide_ibox=wide_ibox,
        max_ibox=max_ibox,
        max_cntb=max_cntb,
        max_cntb0=max_cntb0,
        max_ncl=max_ncl,
        target_params=target_params,
    )  # max_ncl rows returned

    col_trigger = np.zeros(len(tab2), dtype=int)
    if outroot is not None and len(tab3):# and not ibox64_filter:
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
            prev_trig_time=prev_trig_time,
            min_timedelt=min_timedelt
        )
        #Updating candname for injection audits.
        candname_for_audit = lastname
        # "triggered" means a real voltage trigger (trigtime set by dump_cluster_results_json only for real triggers)
        triggered = bool(trigtime is not None)

        if tab4 is not None and trigger:
            col_trigger = np.where(
                tab4 == tab2, lastname, 0
            )  # if trigger, then overload

            # write all T1 cands
            outputted = cluster_heimdall.dump_cluster_results_heimdall(tab, outroot + f"T1_output{str(np.floor(time.time()).astype('int'))}.csv")

           

    # write T2 clustered/filtered results
    if outroot is not None and len(tab2):
        tab2["trigger"] = col_trigger
        output_file = outroot + "cluster_output" + str(np.floor(time.time()).astype("int")) + ".cand"
        outputted = cluster_heimdall.dump_cluster_results_heimdall(tab2,
                                                                   output_file,
                                                                   min_snr_t2out=min_snr_t2out,
                                                                   max_ncl=max_ncl)

        # aggregate files
        if outputted:
            a = Time.now().mjd
            output_mjd = str(int(a))
            old_mjd = str(int(a)-1)
            fl1 = outroot+old_mjd+".csv"
            fl2 = outroot+output_mjd+".csv"
            ofl = outroot+"cluster_output.csv"

#            os.system("cat "+output_file+" >> "+outroot+output_mjd+".csv")
#            os.system("if ! grep -Fxq 'snr,if,specnum,mjds,ibox,idm,dm,ibeam,cl,cntc,cntb,trigger' "+outroot+output_mjd+".csv; then sed -i '1s/^/snr\,if\,specnum\,mjds\,ibox\,idm\,dm\,ibeam\,cl\,cntc\,cntb\,trigger\\n/' "+outroot+output_mjd+".csv; fi")

            df0 = pandas.read_csv(output_file, delimiter=' ', names=columns, on_bad_lines='warn')

            dfs = [df0]
            if os.path.exists(fl1):  # accumulate to yesterday's for rolling 2-day file
                df1 = pandas.read_csv(fl1, on_bad_lines='warn')
                dfs.append(df1)

            if os.path.exists(fl2):  # accumulate to today's for 1-day file
                df2 = pandas.read_csv(fl2, on_bad_lines='warn')
                dfs.append(df2)
                dfc2 = pandas.concat( (df0, df2) )
                dfc2.to_csv(fl2, index=False)
            else:
                df0.to_csv(fl2, index=False)

            dfc = pandas.concat(dfs)
            dfc.to_csv(ofl, index=False)
    

    # --- Finalize audit (stamps G3..G7) ---
    if audit_enabled and auditor is not None:
        try:
            if prev_trig_time is None:
                # No prior real trigger this session: cooldown is effectively satisfied
                cooldown_wait_s = float('inf')
            else:
                try:
                    cooldown_wait_s = float((Time.now() - prev_trig_time).to_value('sec'))
                except Exception:
                    cooldown_wait_s = float('inf')

            nbeams_this_gulp = int(nbeams_gulp)
            thresholds_for_finalize = dict(thresholds_snapshot)
            thresholds_for_finalize["t2_output_file"] = output_file if outroot is not None else ""
            thresholds_for_finalize["exception"] = ""

            auditor.finalize_from_cluster_result(
                host="T2",                              
                gulp=int(gulp) if gulp is not None else -1,
                tab_pre_filter=tab_pre_flag,
                tab_after_beam_flag=tab_after_flag,
                tab_peak=tab2,
                tab_after_filters=tab3,
                lastname=lastname,                     
                trigtime=trigtime,
                triggered=bool(triggered),
                nbeams_this_gulp=nbeams_this_gulp,
                nbeams_queue_snapshot=list(nbeams_queue),
                nbeams_queue_sum=sum(nbeams_queue) if len(nbeams_queue) else 0,
                thresholds=thresholds_for_finalize,
                cooldown_wait_s=cooldown_wait_s,
                candname_for_injection=(candname_for_audit if 'candname_for_audit' in locals() else None),
            )
        except Exception as e:
            logger.warning(f"[AUDIT] finalize_from_cluster_result failed: {e}")


    return lastname, trigtime, triggered


def recvall(sock, n):
    """
    helper function to receive all bytes from a socket
    sock: open socket
    n: maximum number of bytes to expect. you can make this ginormous!
    """
    
    data = bytearray()
    while len(data) < n-1:
        try:
            packet = sock.recv(n - len(data))
            if not packet:
                return data
        except socket.timeout:
            print(f"Only received {n} of 1e6")
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
