#!/usr/bin/env python3
import os, time, random
import numpy as np
from astropy.time import Time
import dsautils.dsa_store as ds
import slack_sdk as slack

# ----------------- config -----------------
SLACK_CHANNEL = "candidates"
SLACK_TOKEN_FILE = f"{os.path.expanduser('~')}/.config/slack_api"
INJ_LIST = "/home/ubuntu/data/injections/injection_list.txt"
AUDIT_INJECTION_FILENAME = "/operations/T2/injection_audit_results/injections_for_audit.txt" 
PARAMS_TXT = "/home/ubuntu/simulated_frb_params.txt"  # columns: DM SNR width_fwhm spec_ind
TEMPLATES = [
    "/home/ubuntu/data/burst_0.inject",
    "/home/ubuntu/data/burst_1.inject",
    "/home/ubuntu/data/burst_2.inject",
    "/home/ubuntu/data/burst_3.inject",
    "/home/ubuntu/data/burst_4.inject",
]

SLEEP_SEC = 600  # 10 minutes

NODES = [17, 18, 19, 20]
PAIR  = {17: 19, 18: 20, 19: 17, 20: 18}
LOCAL_BEAMS_PER_NODE = 128

# Legacy output formatting (unchanged)
FMT_OUT = "%5.9f  %d  %0.2f %0.1f %0.3f %0.2f %s\n"

# --- SNR→scale mapping (Recovered_SNR ≈ k * scale) ---
K_DEFAULT = 135.0   # 27/.2 from old script
SCALE_MIN = 0.075
SCALE_MAX = 0.3
niterations = 288
scale = np.random.uniform(SCALE_MIN, SCALE_MAX, niterations)
# ------------------------------------------------------

def ensure_slack():
    if not os.path.exists(SLACK_TOKEN_FILE):
        raise RuntimeError(f"Could not find file with slack api token at {SLACK_TOKEN_FILE}")
    with open(SLACK_TOKEN_FILE, "r") as sf_handler:
        slack_token = sf_handler.read()
    return slack.WebClient(token=slack_token)

def slack_msg(cli, text):
    try:
        cli.chat_postMessage(channel=SLACK_CHANNEL, text=text)
    except Exception as e:
        print(f"[slack] {e}")

def ensure_injection_list(filename):
    if not os.path.exists(filename):
        with open(filename, "w") as f:
            f.write("# MJD   Beam   DM    SNR   Width_fwhm   spec_ind  FRBno\n")

def global_to_node_local(g):
    if not (0 <= g <= 511):
        raise ValueError("global beam must be in [0, 511]")
    group = g // LOCAL_BEAMS_PER_NODE
    node  = NODES[group]
    local = g %  LOCAL_BEAMS_PER_NODE
    return node, local

def load_params_and_templates():
    params = np.genfromtxt(PARAMS_TXT)  # (N,4): DM, SNR, width_fwhm, spec_ind
    n = min(len(params), len(TEMPLATES))
    if n == 0:
        raise RuntimeError("No params/templates found")
    if n < len(params):
        print(f"[warn] Truncating params to {n} to match templates")
    return params[:n], TEMPLATES[:n]

def snr_to_scale(target_snr: float) -> float:
    sc = float(target_snr) / float(K_DEFAULT)
    if sc < SCALE_MIN: sc = SCALE_MIN
    if sc > SCALE_MAX: sc = SCALE_MAX
    return sc

def scale_to_snr(scale: float) -> float:
    return float(scale * K_DEFAULT)

def run():
    random.seed()
    np.random.seed()
    ensure_injection_list(INJ_LIST)
    ensure_injection_list(AUDIT_INJECTION_FILENAME)

    slack_cli = ensure_slack()
    store = ds.DsaStore()

    params, templates = load_params_and_templates()

    for i in range(niterations):
        current_scale = scale[i]
        print(f"[info] Starting iteration {i+1}/{niterations}")
        for idx, row in enumerate(params):
            DM, SNR, width_fwhm, spec_ind = map(float, row)   # SNR is fixed from file (e.g., 15.0)
            template = templates[idx]
            frbno = os.path.splitext(os.path.basename(template))[0].split("_")[-1]

            # choose a random global beam, then derive EW/NS pair
            gbeam = random.randint(0, 511)
            print(f"[info] Injecting FRB {frbno} with DM={DM}, SNR={SNR}, into global beam {gbeam} using template {template}")

            if gbeam <= 255:
                ew_beam, ns_beam = gbeam, gbeam + 256
            else:
                ns_beam, ew_beam = gbeam, gbeam - 256

            node, local = global_to_node_local(gbeam)
            partner = PAIR[node]

            
            effective_snr = scale_to_snr(current_scale)

             # event bookkeeping (legacy append)
            mjd = Time.now().mjd
            with open(INJ_LIST, "a") as f:
                # Beam column uses the global beam you selected
                f.write(FMT_OUT % (mjd, gbeam, DM, effective_snr, width_fwhm, spec_ind, frbno))
            
            with open(AUDIT_INJECTION_FILENAME, "a") as f:
                f.write(FMT_OUT % (mjd, gbeam, DM, effective_snr, width_fwhm, spec_ind, frbno))


            injection_cmd = {"cmd": "inject", "val": f"{local}-{template}-{current_scale:.3f}-"}
            store.put_dict(f"/cmd/corr/{node}", injection_cmd)
            store.put_dict(f"/cmd/corr/{partner}", injection_cmd)
            

            slack_msg(
                slack_cli,
                (f"Sending injection to beam {ew_beam} (EW) and {ns_beam} (NS) "
                  f"with DM={DM:.1f} and SNR={effective_snr:.1f}"
                )
            )

            time.sleep(SLEEP_SEC)

if __name__ == "__main__":
    run()
