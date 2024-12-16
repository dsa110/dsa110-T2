#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# dsahead python 3.7
import json
import os.path

# from sklearn import cluster  # for dbscan
#import hdbscan
from sklearn.cluster import DBSCAN
import numpy as np
from astropy import time, units, table
from astropy.io import ascii
from astropy.io.ascii.core import InconsistentTableError

try:
    from T2 import triggering
except ModuleNotFoundError:
    print("not importing triggering")
try:
    from dsautils import coordinates, dsa_store
    ds = dsa_store.DsaStore()
    import dsautils.dsa_syslog as dsl
    logger = dsl.DsaSyslogger()
    logger.subsystem("software")
    logger.app("T2")
except:
    print("not importing dsautils")
    ds = None
    import logging
    logger = logging.getLogger()

from event import names  # TODO: add event to get DSAEvent class

# half second at heimdall time resolution (after march 18)
offset = 1907
downsample = 4
NSNR = 10


def parse_candsfile(candsfile):
    """Takes standard MBHeimdall giants output and returns full table, classifier inputs and snr tables.
    (Can add cleaning here, eventually)
    """

    if os.path.exists(candsfile):
        logger.debug(f"Candsfile {candsfile} is path, so opening it")
        candsfile = open(candsfile, "r").read()
    else:
        ncands = len(candsfile.split("\n")) - 1
        logger.debug(f"Received {ncands} candidates")
    #    candsfile = '\n'.join([line for line in candsfile.split('\n') if line.count(' ') == 7])
    #    print(f'Received {ncands0} candidates, removed {ncands0-ncands} lines.')
    col_heimdall = ["snr", "if", "itime", "mjds", "ibox", "idm", "dm", "ibeam"]
    col_T2old = ["snr", "if", "itime", "mjds", "ibox", "idm", "dm", "ibeam",
                 "cl", "cntc", "cntb"]
    col_T2nobeams = ["snr", "if", "itime", "mjds", "ibox",
                     "idm", "dm", "ibeam", "cl", "cntc", "cntb", "trigger"]
    col_T2 = ['snr','if','itime','mjds','ibox','idm','dm',
              'ibeam','cl','cntc','cntb','snrs0','beams0','snrs1',
              'beams1','snrs2','beams2','snrs3','beams3','snrs4',
              'beams4','snrs5','beams5','snrs6','beams6','snrs7',
              'beams7','snrs8','beams8','snrs9','beams9','trigger']

    # flag for heimdall file
    hdfile = False

    try:
        tab = ascii.read(
            candsfile,
            names=col_heimdall,
            guess=True,
            fast_reader=False,
            format="no_header",
        )
        hdfile = True
        logger.debug("Read with heimdall columns")
    except InconsistentTableError:
        try:
            tab = ascii.read(
                candsfile,
                names=col_T2,
                guess=True,
                fast_reader=False,
                format="no_header",
            )
            hdfile = False
            logger.debug("Read with T2 columns")
        except InconsistentTableError:
            try:
                tab = ascii.read(
                    candsfile,
                    names=col_T2nobeams,
                    guess=True,
                    fast_reader=False,
                    format="no_header",
                )
                hdfile = False
                logger.debug("Read with old style T2 columns")
            except InconsistentTableError:
                try:
                    tab = ascii.read(
                        candsfile,
                        names=col_T2old,
                        guess=True,
                        fast_reader=False,
                        format="no_header",
                    )
                    hdfile = False
                    logger.debug("Read with old style T2 columns")
                except InconsistentTableError:
                    logger.warning("Inconsistent table. Skipping...")
                    return []

    tab["ibeam"] = tab["ibeam"].astype(int)
    tab["idm"] = tab["idm"].astype(int)
    tab["ibox"] = tab["ibox"].astype(int)
    tab["itime"] = tab["itime"].astype(int)
    if hdfile is True:
        try:
            ret_time = (
                ds.get_dict("/mon/snap/1/armed_mjd")["armed_mjd"]
                + float(ds.get_dict("/mon/snap/1/utc_start")["utc_start"])
                * 4.0
                * 8.192e-6
                / 86400.0
            )
        except:
            ret_time = 55000.0
#        print(ret_time)
        tab["mjds"] = tab["mjds"] + ret_time

    #
    #    snrs = tab['snr']
    # how to use ibeam?

    #    return tab, data, snrs
    return tab

def flag_beams(tab,stat=5e5):

    tab2 = tab
    bms = np.unique(np.asarray(tab["ibeam"]))
    stds = np.zeros(len(bms))
    for i in np.arange(len(bms)):
        tt = tab[tab["ibeam"]==bms[i]]["dm"]
        stds[i] = np.std(tt)
        stds[i] *= len(tt)
        if stds[i]>stat:
            print(f"Flagging beam {bms[i]}")
            tab2 = tab2[tab2["ibeam"]!=bms[i]]
            
    return tab2
                                                                                            

def cluster_data(
    tab,
    selectcols=["itime", "idm", "ibox"],  # no ibeam for 2-arm clustering
    min_cluster_size=2,
    min_samples=2,
    metric="cityblock",
    return_clusterer=False,
    allow_single_cluster=True,
):
    """Take data from parse_candsfile and identify clusters via hamming metric.
    selectcols will take a subset of the standard MBHeimdall output
    """

    tt = tab[selectcols]
#    print(tt[:10])
    data = np.lib.recfunctions.structured_to_unstructured(
        tab[selectcols].as_array(), dtype=np.int
    )  # ok for single dtype (int)
#    np.savez("test.npz",data=data)
    
    try:
        #clusterer = hdbscan.HDBSCAN(
        #    metric=metric,
        #    min_cluster_size=min_cluster_size,
        #    min_samples=min_samples,
        #    cluster_selection_method="eom",
        #    allow_single_cluster=allow_single_cluster,
        #).fit(data)
        clusterer = DBSCAN(metric=metric, min_samples=min_samples,
                           eps=10, algorithm='auto', leaf_size=23).fit(data)

        nclustered = np.max(clusterer.labels_ + 1)
        nunclustered = len(np.where(clusterer.labels_ == -1)[0])
        cl = clusterer.labels_
    except ValueError:
        print("Clustering did not run. Each point assigned to unique cluster.")
        logger.info(
            "Clustering did not run. Each point assigned to unique cluster."
        )
        cl = np.arange(len(data))
        nclustered = 0
        nunclustered = len(cl)

    #    logger.info(f'Found {nclustered} clustered and {nunclustered} unclustered rows')

    bl = tab['ibeam'].astype(int).data  # use tab to get ibeam values
    cntb = np.zeros((len(data), 1), dtype=int)
    cntc = np.zeros((len(data), 1), dtype=int)
    ucl = np.unique(cl)

    for i in ucl:
        ww = np.where(i == cl)
        cntc[ww] = len(ww[0])   # TODO: figure out how to count for ns and ew separately
        ubl = np.unique(bl[ww])
        cntb[ww] = len(ubl)

    # append useful metastats to original data
    #    data_labeled = np.hstack((data, cl[:,None], cntb, cntc))
    # modifies tab in place
    tab["cl"] = cl.tolist()
    tab["cntc"] = cntc.flatten().tolist()
    tab["cntb"] = cntb.flatten().tolist()

    if return_clusterer:
        #        return clusterer, data_labeled
        return clusterer


def get_peak(tab, nsnr=NSNR):
    """Given labeled data of search output from two arms, find snr distribution per cluster
    Adds in count of candidates in same beam and same cluster.
    Add SNR from top nsnr beams as new columns (snr1, snr2, ..., snr<nsnr>) per arm.
    """

    beams_ns = []  # nsnr-tuple per cluster
    snrs_ns = []
    beams_ew = []
    snrs_ew = []
    inds_peak = []

    cl = tab["cl"].astype(int)
    ncl = len(np.unique(cl))
    snrs = tab["snr"]
    ipeak = []

    print("finding unique clusters")
    
    for c in np.unique(cl):
        if c == -1:
            continue

        clusterinds = np.where(c == cl)[0]
        maxsnr = snrs[clusterinds].max()
        imaxsnr = np.where(snrs == maxsnr)[0][0]
        ipeak.append(imaxsnr)

#    ipeak += [i for i in range(len(tab)) if cl[i] == -1]  # TODO figure out how to handle unclustered events
    tab2 = tab[ipeak]
    ncl = len(ipeak)

    # get top beam snrs and attach as new columns
    beams = np.zeros((nsnr, ncl), dtype=int)
    snrs = np.zeros((nsnr, ncl), dtype=float)
    iii = 0

    print("getting max snr per beam")
    
    for c in np.unique(cl):
        if c == -1:
            continue
        
        # iterate to get max snr per beam
        ss = np.zeros(512, dtype=float)
        tcl = tab[tab['cl'] == c]   # TODO this is probably very slow        
        for i in np.unique(tcl['ibeam']):
            maxsnr = tcl[tcl['ibeam'] == i]['snr'].max()
            ss[i] = maxsnr
        beams[:, iii] = ss.argsort()[::-1][:nsnr]
        snrs[:, iii] = ss[ss.argsort()[::-1]][:nsnr]
        iii += 1
        
    try:
        for i in range(nsnr):
            tab2[f'snrs{i}'] = list(snrs[i])
            tab2[f'beams{i}'] = list(beams[i])
    except:
        print("Error in adding beam SNRs")
        for i in range(nsnr):
            print(list(snrs[i]))
        beams = np.zeros((nsnr, ncl), dtype=int)
        snrs = np.zeros((nsnr, ncl), dtype=float)
        for i in range(nsnr):
            tab2[f'snrs{i}'] = list(snrs[i])
            tab2[f'beams{i}'] = list(beams[i])

    logger.info(f"Found {len(ipeak)} cluster peaks")
    print(f"Found {len(ipeak)} cluster peaks")
    logger.info(f"Got top {nsnr} beams from {ncl} clusters")
    print(f"Got top {nsnr} beams from {ncl} clusters")

    return tab2


def filter_clustered(
        tab,
        min_dm=50,
        min_snr=8,
        min_snr_wide=9,
        min_snr_1arm=10,
        wide_ibox=17,
        max_ibox=33,
        min_cntb=None,
        max_cntb=None,
        max_cntb0=None,
        min_cntc=None,
        max_cntc=None,
        max_ncl=None,
        target_params=None,
        frac_wide=0.0,
        nsnr=NSNR
):
    """Function to select a subset of clustered output.
    Can set minimum SNR, min/max number of beams in cluster, min/max total count in cluster.
    target_params is a tuple (min_dmt, max_dmt, min_snrt) for custom snr threshold for target.
    max_ncl is maximum number of clusters returned (sorted by SNR).
    """

    if target_params is not None:
        min_dmt, max_dmt, min_snrt = target_params
    else:
        min_dmt, max_dmt, min_snrt = None, None, None

    good = [True] * len(tab)

    if min_snr is not None:
        if min_snrt is None:
            # snr limit for narrow and wide, with requirement of at least two beams
            df = tab.to_pandas()
            nsarr = ((df[[f'beams{i}' for i in range(nsnr)]].values > 255)) & (df[[f'snrs{i}' for i in range(nsnr)]].values > 0)
            ewarr = ((df[[f'beams{i}' for i in range(nsnr)]].values <= 255)) & (df[[f'snrs{i}' for i in range(nsnr)]].values > 0)
            twoarm = ewarr.any(axis=1) & nsarr.any(axis=1)
            twoarm = (ewarr.any(axis=1) & nsarr.any(axis=1)) | (df['snr'].values > min_snr_1arm).any()
            print(f'nsarr: {nsarr}, ewarr: {ewarr}, twoarm: {twoarm}')

            good0 = (tab["snr"] > min_snr) * (tab["ibox"] < wide_ibox)
            good1 = (tab["snr"] > min_snr_wide) * (tab["ibox"] >= wide_ibox)
            print(f'good0: {good0}; good1: {good1}')
            good0 *= twoarm
            good1 *= twoarm
            good *= good0 + good1
        else:
            # print(f'min_snr={min_snr}, min_snrt={min_snrt}, min_dmt={min_dmt}, max_dmt={max_dmt}, tab={tab[["snr", "dm"]]}')
            good0 = (tab["snr"] > min_snr) * (tab["dm"] > max_dmt)
            good1 = (tab["snr"] > min_snr) * (tab["dm"] < min_dmt)
            good2 = (
                (tab["snr"] > min_snrt)
                * (tab["dm"] > min_dmt)
                * (tab["dm"] < max_dmt)
            )
            good *= good0 + good1 + good2
            # print('good0, good1, good2, good:')
            # print(good0, good1, good2, good)
    if min_dm is not None:
        good *= tab["dm"] > min_dm
    if max_ibox is not None:
#        print(f"using max_ibox {max_ibox}")
        good *= tab["ibox"] < max_ibox
    if min_cntb is not None:
        good *= tab["cntb"] > min_cntb
    if max_cntb is not None:
        good *= tab["cntb"] < max_cntb
    if min_cntc is not None:
        good *= tab["cntc"] > min_cntc
    if max_cntc is not None:
        good *= tab["cntc"] < max_cntc

    tab_out = tab[good]

    if max_ncl is not None:
        if len(tab_out) > max_ncl:
            min_snr_cl = sorted(tab_out["snr"])[-max_ncl]
            good = tab_out["snr"] >= min_snr_cl
            tab_out = tab_out[good]
            print(
                f"Limiting output to {max_ncl} clusters with snr>{min_snr_cl}."
            )

    logger.info(
        f"Filtering clusters from {len(tab)} to {len(tab_out)} candidates."
    )
    print(f"Filtering clusters from {len(tab)} to {len(tab_out)} candidates.")

    return tab_out


def get_nbeams(tab, threshold=7.5):
    """Calculate number of beams in table above SNR threshold."""

    goody = [True] * len(tab)
    goody *= tab["snr"] > threshold
    tab_out2 = tab[goody]
    if len(tab_out2) > 0:
        ibeams = np.asarray(tab_out2["ibeam"])
        nbeams = len(np.unique(ibeams))
    else:
        nbeams = 0

    return nbeams


def dump_cluster_results_json(
        tab,
        outputfile=None,
        output_cols=["mjds", "snr", "ibox", "dm", "ibeam", "cntb", "cntc"] + [f'snrs{i}' for i in range(NSNR)] + [f'beams{i}' for i in range(NSNR)],
        trigger=False,
        lastname=None,
        gulp=None,
        cat=None,
        beam_model=None,
        coords=None,
        snrs=None,
        outroot="",
        nbeams=0,
        max_nbeams=40,
        frac_wide=0.0,
        injectionfile='/home/ubuntu/data/injections/injection_list.txt',
        prev_trig_time=None,
        min_timedelt=1.    
    ):
    """
    Takes tab from parse_candsfile and clsnr from get_peak,
    json file will be named with generated name, unless outputfile is set
    TODO: make cleaner, as it currently assumes NSNR at compile time.
    candidate name and specnum is calculated. name is unique.
    trigger is bool to update DsaStore to trigger data dump.
    cat is path to source catalog (default None)
    beam_model is pre-calculated beam model (default None)
    coords and snrs are parsed source file input
    injectionfile is path to info on injects and controls whether trigger is compared to that
    returns row of table that triggered, along with name generated for candidate.
    """

    if coords is None or snrs is None:
        if cat is not None:
            coords, snrs = triggering.parse_catalog(cat)

    itimes = tab["itime"]
    maxsnr = tab["snr"].max()
    imaxsnr = np.where(tab["snr"] == maxsnr)[0][0]
    itime = str(itimes[imaxsnr])
    specnum = (int(itimes[imaxsnr]) - offset) * downsample
    mjd = tab["mjds"][imaxsnr]
    snr = tab["snr"][imaxsnr]
    dm = tab["dm"][imaxsnr]
    ibeam = tab["ibeam"][imaxsnr]
    
    isinjection = False
    if injectionfile is not None:
        # check candidate against injectionfile
        try:
            tab_inj = ascii.read(injectionfile)
        except:
            tab_inj = ascii.read(injectionfile, names="MJD   Beam   DM    SNR   Width_fwhm   spec_ind  FRBno".split())
        finally:
            assert all([col in tab_inj.columns for col in ["MJD", "Beam", "DM", "SNR", "FRBno"]])

        # is candidate proximal to any in tab_inj?
        t_close = 300 # seconds  TODO: why not 1 sec?
        dm_close = 20 # pc/cm3
        beam_close = 2 # number
        sel_t = np.abs(tab_inj["MJD"] - mjd) < t_close/(3600*24)
        sel_dm = np.abs(tab_inj["DM"] - dm) < dm_close
        sel_beam = np.abs(tab_inj["Beam"] - ibeam) < beam_close
        sel_beam_2 = np.abs(tab_inj["Beam"]+256 - ibeam) < beam_close
        print(f"INJECTION TEST: min abs time diff {np.abs((tab_inj['MJD']-mjd)*24*3600).min()} seconds. Sel any? t {sel_t.any()}, dm {sel_dm.any()}, beam {sel_beam.any()}, beam2 {sel_beam_2.any()}")
        sel = sel_t*sel_dm*sel_beam
        sel2 = sel_t*sel_dm*sel_beam_2
        if len(np.where(sel)[0]):
            isinjection = True
            selt = sel
        if len(np.where(sel2)[0]):
            isinjection = True
            selt = sel2

        if time.Time.now().mjd - mjd > 13:
            logger.warning("Event MJD is {mjd}, which is more than 13 days in the past. SNAP counter overflow?")

    if isinjection:
        basename = names.increment_name(mjd, lastname=lastname)
        print(tab_inj[selt])
        candname = f"{basename}_inj{tab_inj[selt][-1]['FRBno']}"
        print(f"Candidate identified as injection. Naming it {candname}")
        if len(selt) > 1:
            print(f"Found {len(selt)} injections coincident with this event. Using first.")
        # if injection is found, skip the voltage trigger via etcd
    else:
        # if no injection file or no coincident injection
        candname = names.increment_name(mjd, lastname=lastname)

    output_dict = {candname: {}}
    if outputfile is None:
        outputfile = f"{outroot}cluster_output{candname}.json"

    row = tab[imaxsnr]
    red_tab = tab[imaxsnr : imaxsnr + 1]
    for col in output_cols:
        if type(row[col]) == np.int64:
            output_dict[candname][col] = int(row[col])
        else:
            output_dict[candname][col] = row[col]

    output_dict[candname]["specnum"] = specnum
    (
        output_dict[candname]["ra"],
        output_dict[candname]["dec"],
    ) = get_radec(beamnum=127)  # quick and dirty

    if gulp is not None:
        output_dict[candname]["gulp"] = gulp

    if isinjection:  # add in any case?
        output_dict[candname]['injected'] = isinjection

    nbeams_condition = False
    if nbeams > max_nbeams:
        nbeams_condition = True
        print(f"nbeams condition is {nbeams_condition} ({nbeams}>{max_nbeams}). Will not trigger.")
        # Liam edit to preserve real FRBs during RFI storm:
        # if nbeam > max_nbeams and frac_wide < 0.8: do not discard because
        # most FPs are wide
#        if frac_wide < 0.8:
#            nbeams_condition = False

    if len(tab) and nbeams_condition is False:
        print(red_tab)

        # TODO: create DSAEvent here and use it instead of output_dict

        if cat is not None and red_tab is not None:
            # beam_model = triggering.read_beam_model(beam_model)
            tab_checked = triggering.check_clustered_sources(
                red_tab, coords, snrs, beam_model=beam_model, do_check=False
            )
            if len(tab_checked):

                if prev_trig_time is not None:
                    if time.Time.now()-prev_trig_time < min_timedelt*units.s:
                        print(f"Not triggering because of short wait time")
                        logger.info(f"Not triggering because of short wait time")
                        return None, candname, None
                    else:
                        trigtime = None
                        
                with open(outputfile, "w") as f:  # encoding='utf-8'
                    print(
                        f"Writing trigger file for index {imaxsnr} with SNR={maxsnr}"
                    )
                    logger.info(
                        f"Writing trigger file for index {imaxsnr} with SNR={maxsnr}"
                    )
                    json.dump(output_dict, f, ensure_ascii=False, indent=4)   # could replace this with DSAEvent method

                if trigger and time.Time.now().mjd - mjd < 13:  #  and not isinjection ?
                    send_trigger(output_dict=output_dict)
                    trigtime = time.Time.now()
                else:
                    trigtime = None

                return row, candname, trigtime

            else:
                print(f"Not triggering on source in beam")
                logger.info(f"Not triggering on source in beam")
                return None, candname, None

        else:
            if prev_trig_time is not None:
                if time.Time.now()-prev_trig_time < min_timedelt*units.s:
                    print(f"Not triggering because of short wait time")
                    logger.info(f"Not triggering because of short wait time")
                    return None, candname, None
                else:
                    trigtime = None
                
            with open(outputfile, "w") as f:  # encoding='utf-8'
                print(
                    f"Writing trigger file for index {imaxsnr} with SNR={maxsnr}"
                )
                logger.info(
                    f"Writing trigger file for index {imaxsnr} with SNR={maxsnr}"
                )
                json.dump(output_dict, f, ensure_ascii=False, indent=4)

            if trigger and time.Time.now().mjd - mjd < 13:  #  and not isinjection ?
                send_trigger(output_dict=output_dict)
                trigtime = time.Time.now()
            else:
                trigtime = None

            return row, candname, trigtime

    else:
        print(
            f"Not triggering on block with {len(tab)} candidates and nbeams {nbeams}>{max_nbeams} beam count sum"
        )
        logger.info(
            f"Not triggering on block with {len(tab)} candidates and nbeams {nbeams}>{max_nbeams} beam count sum"
        )
        return None, lastname, None

    print("Not triggering on nbeams condition")
    return None, lastname, None


def get_radec(mjd=None, beamnum=None):
    """Use time, beam number, and and antenna elevation to get RA, Dec of beam."""

    if mjd is not None:
        print("Using time to get ra,dec")
        tt = time.Time(mjd, format="mjd")
    else:
        tt = None

    ra, dec = coordinates.get_pointing(ibeam=beamnum, obstime=tt)    
    
    return ra.value, dec.value


def send_trigger(output_dict=None, outputfile=None):
    """Use either json file or dict to send trigger for voltage dumps via etcd."""

    if outputfile is not None:
        print("Overloading output_dict trigger info with that from outputfile")
        logger.info(
            "Overloading output_dict trigger info with that from outputfile"
        )
        with open(outputfile, "w") as f:
            output_dict = json.load(f)

    candname = list(output_dict)[0]
    val = output_dict.get(candname)
    print(candname, val)
    print(
        f"Sending trigger for candidate {candname} with specnum {val['specnum']}"
    )
    logger.info(
        f"Sending trigger for candidate {candname} with specnum {val['specnum']}"
    )

    with open(f"/home/ubuntu/data/T2test/{candname}.json", "w") as f:  # encoding='utf-8'
        print(
            f"Writing dump dict"
        )
        json.dump({"cmd": "trigger", "val": f'{val["specnum"]}-{candname}-'}, f, ensure_ascii=False, indent=4)

    ds.put_dict(
        "/cmd/corr/0",
        {"cmd": "trigger", "val": f'{val["specnum"]}-{candname}-'},
    )  # triggers voltage dump in corr.py
    ds.put_dict(
        "/mon/corr/1/trigger", output_dict
    )  # tells look_after_dumps.py to manage data


def dump_cluster_results_heimdall(
    tab, outputfile, min_snr_t2out=None, max_ncl=None
):
    """
    Takes tab from parse_candsfile and clsnr from get_peak,
    output T2-clustered results.
    Both T1 (heimdall) and T2 outputs are named "*.cand" and can be parsed by parse_candsfile().
    min_snr_t2out is a min snr on candidates to write.
    max_ncl is number of rows to write.
    """

    tab["itime"] = (tab["itime"] - offset) * downsample  # transform to specnum

    if min_snr_t2out is not None:
        good = [True] * len(tab)
        good *= tab["snr"] > min_snr_t2out
        tab = tab[good]
        if not all(good) and len(tab):
            print(
                f"Limiting output to SNR>{min_snr_t2out} with {len(tab)} clusters."
            )

    if max_ncl is not None:
        if len(tab) > max_ncl:
            min_snr_cl = sorted(tab["snr"])[-max_ncl]
            good = (tab["snr"] >= min_snr_cl) + [
                str(tt) != "0" for tt in tab["trigger"]
            ]  # keep trigger
            tab = tab[good]
            print(
                f"Limiting output to {max_ncl} clusters with snr>{min_snr_cl}."
            )
    else:
        print("max_ncl not set. Not filtering heimdall output file.")

    if len(tab) > 0:
        tab.write(outputfile, format="ascii.no_header", overwrite=True)
        return True

    return False
        
    
