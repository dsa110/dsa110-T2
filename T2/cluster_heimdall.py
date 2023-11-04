#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# dsahead python 3.7
import json
import os.path

# from sklearn import cluster  # for dbscan
import hdbscan
import numpy as np
from astropy import time, units
from astropy.io import ascii
from astropy.io.ascii.core import InconsistentTableError

try:
    from T2 import triggering
except ModuleNotFoundError:
    print("not importing triggering")
from dsautils import coordinates, dsa_store
from event import names  # TODO: add event to get DSAEvent class

ds = dsa_store.DsaStore()
import dsautils.dsa_syslog as dsl

logger = dsl.DsaSyslogger()
logger.subsystem("software")
logger.app("T2")

# half second at heimdall time resolution (after march 18)
offset = 1907
downsample = 4


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
    col_T2old = [
        "snr",
        "if",
        "itime",
        "mjds",
        "ibox",
        "idm",
        "dm",
        "ibeam",
        "cl",
        "cntc",
        "cntb",
    ]
    col_T2 = [
        "snr",
        "if",
        "itime",
        "mjds",
        "ibox",
        "idm",
        "dm",
        "ibeam",
        "cl",
        "cntc",
        "cntb",
        "trigger",
    ]

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
                    names=col_T2old,
                    guess=True,
                    fast_reader=False,
                    format="no_header",
                )
                hdfile = False
                logger.debug("Read with old style T2 columns")
            except InconsistentTableError:
                logger.warning("Inconsistent table. Skipping...")
                return ([], [], [])

    tab["ibeam"] = tab["ibeam"].astype(int)
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
        tab["mjds"] = tab["mjds"] / 86400.0 + ret_time

    #
    #    snrs = tab['snr']
    # how to use ibeam?

    #    return tab, data, snrs
    return tab


def cluster_data(
    tab,
    selectcols=["itime", "idm", "ibox", "ibeam"],
    min_cluster_size=2,
    min_samples=5,
    metric="hamming",
    return_clusterer=False,
    allow_single_cluster=True,
):
    """Take data from parse_candsfile and identify clusters via hamming metric.
    selectcols will take a subset of the standard MBHeimdall output
    """

    data = np.lib.recfunctions.structured_to_unstructured(
        tab[selectcols].as_array()
    )  # ok for single dtype (int)
    try:
        clusterer = hdbscan.HDBSCAN(
            metric=metric,
            min_cluster_size=min_cluster_size,
            min_samples=min_samples,
            cluster_selection_method="eom",
            allow_single_cluster=allow_single_cluster,
        ).fit(data)
        #        clusterer = cluster.DBSCAN(metric='chebyshev', min_samples=min_samples,
        #                                   eps=14, algorithm='auto', leaf_size=23).fit(data)

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

    # hack assumes fixed columns
    bl = data[:, 3]
    cntb, cntc = np.zeros((len(data), 1), dtype=int), np.zeros(
        (len(data), 1), dtype=int
    )
    ucl = np.unique(cl)

    for i in ucl:
        ww = np.where(i == cl)
        cntc[ww] = len(ww[0])
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


#    else:
#        return data_labeled


def get_peak(tab):
    """Given labeled data, find max snr row per cluster
    Adds in count of candidates in same beam and same cluster.
    Puts unclustered candidates in as individual events.
    """

    #    clsnr = []
    #    cl = datal[:, 4].astype(int)   # hack. should really use table.
    #    cnt_beam = datal[:, 5].astype(int)
    #    cnt_cl = datal[:, 6].astype(int)
    cl = tab["cl"].astype(int)
    #    cnt_beam = tab['cntb'].astype(int)
    #    cnt_cl = tab['cntc'].astype(int)
    snrs = tab["snr"]
    ipeak = []
    for i in np.unique(cl):
        if i == -1:
            continue
        clusterinds = np.where(i == cl)[0]
        maxsnr = snrs[clusterinds].max()
        imaxsnr = np.where(snrs == maxsnr)[0][0]
        ipeak.append(imaxsnr)
    #        clsnr.append((imaxsnr, maxsnr, cnt_beam[imaxsnr], cnt_cl[imaxsnr]))
    ipeak += [i for i in range(len(tab)) if cl[i] == -1]  # append unclustered
    logger.info(f"Found {len(ipeak)} cluster peaks")
    print(f"Found {len(ipeak)} cluster peaks")

    return tab[ipeak]


def filter_clustered(
        tab,
        min_dm=50,
        min_snr=8,
        min_snr_wide=9,
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
#            good *= tab["snr"] > min_snr
            good0 = (tab["snr"] > min_snr) * (tab["ibox"] < wide_ibox)
            good1 = (tab["snr"] > min_snr_wide) * (tab["ibox"] >= wide_ibox)
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
        print(f"using max_ibox {max_ibox}")
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
        output_cols=["mjds", "snr", "ibox", "dm", "ibeam", "cntb", "cntc"],
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
        injectionfile='/home/ubuntu/injection_list.txt',
        prev_trig_time=None,
        min_timedelt=1.    
    ):
    """
    Takes tab from parse_candsfile and clsnr from get_peak,
    json file will be named with generated name, unless outputfile is set
    candidate name and specnum is calculated. name is unique.
    trigger is bool to update DsaStore to trigger data dump.
    cat is path to source catalog (default None)
    beam_model is pre-calculated beam model (default None)
    coords and snrs are parsed source file input
    injectionfile is path to info on injects and controls whether trigger is compared to that
    returns row of table that triggered, along with name generated for candidate.
    """

    if coords is None or snrs is None:
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
        t_close = 20 # seconds  TODO: why not 1 sec?
        dm_close = 20 # pc/cm3
        beam_close = 2 # number
        sel_t = np.abs(tab_inj["MJD"] - mjd) < t_close/(3600*24)
        sel_dm = np.abs(tab_inj["DM"] - dm) < dm_close
        sel_beam = np.abs(tab_inj["Beam"] - ibeam) < beam_close
        sel = sel_t*sel_dm*sel_beam
        if len(np.where(sel)[0]):
            isinjection = True

        if time.Time.now().mjd - mjd > 13:
            logger.warning("Event MJD is {mjd}, which is more than 13 days in the past. SNAP counter overflow?")

    if isinjection:
        basename = names.increment_name(mjd, lastname=lastname)
        candname = f"{basename}_inj{tab_inj[sel]['FRBno'][0]}"
        print(f"Candidate identified as injection. Naming it {candname}")
        if len(sel) > 1:
            print(f"Found {len(sel)} injections coincident with this event. Using first.")
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
    print(f"Checking nbeams condition: {nbeams}>{max_nbeams}")
    if nbeams > max_nbeams:
        nbeams_condition = True
        # Liam edit to preserve real FRBs during RFI storm:
        # if nbeam > max_nbeams and frac_wide < 0.8: do not discard because
        # most FPs are wide
        if frac_wide < 0.8:
            nbeams_condition = False

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
    output T2-clustered results with the same columns as heimdall.cand into a file outputfile.
    The output is in pandas format with column names in the 1st row.
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
        
    
