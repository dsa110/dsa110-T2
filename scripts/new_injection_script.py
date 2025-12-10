#!/usr/bin/env python3
import os, time, random
import numpy as np
from astropy.time import Time
import dsautils.dsa_store as ds
import slack_sdk as slack

# ================== user knobs ==================
# tune only these two to change the population
USER_SNR_MIN = 10.0
USER_SNR_MAX = 25.0

# dry run mode:
#   - no writes to INJ_LIST / AUDIT_INJECTION_FILENAME
#   - no Slack posts
#   - no put_dict() to correlator
#   - no sleep between injections
DRY_RUN = False
# =================================================

# ----------------- beam / node layout -----------------
# Global beams:
#   0–255   : East-West (EW) arm
#       0–127   -> node 17 (EW)
#       128–255 -> node 18 (EW)
#   256–511: North-South (NS) arm
#       256–383 -> node 19 (NS)
#       384–511 -> node 20 (NS)
#
# EW/NS beam pairing:
#   For each EW beam b in [0, 255], NS partner is b + 256.
#
# Node pairing:
#   17 <-> 19
#   18 <-> 20
#
# The EW and NS nodes for a given EW/NS beam pair are separated by +2 in node-id.
EW_TO_NS_NODE_SEP_FACTOR = 2
# -------------------------------------------------------

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
   # "/home/ubuntu/data/burst_4.inject", DM 2500 is not used in search
]

SLEEP_SEC = 1800  # 30 minutes

NODES = [17, 18, 19, 20]
PAIR  = {17: 19, 18: 20, 19: 17, 20: 18}
LOCAL_BEAMS_PER_NODE = 128

# Legacy output formatting
FMT_OUT = "%5.9f  %d  %0.2f %0.1f %0.3f %0.2f %s\n"

# --- amplitude/SNR calibration ---
# recovered_snr ≈ K_DEFAULT * scale
K_DEFAULT = 45.0  # 135/3, adjusted for new beam rms.

# derive scale bounds from user SNR range
# scale = snr / K
_DERIVED_SCALE_MIN = USER_SNR_MIN / K_DEFAULT
_DERIVED_SCALE_MAX = USER_SNR_MAX / K_DEFAULT

# global safety rails for scale
GLOBAL_SCALE_MIN = 0.01
GLOBAL_SCALE_MAX = 0.5

SCALE_MIN = max(_DERIVED_SCALE_MIN, GLOBAL_SCALE_MIN)
SCALE_MAX = min(_DERIVED_SCALE_MAX, GLOBAL_SCALE_MAX)

niterations = 5000
# -----------------------------------


def draw_truncated_euclidean_snr(smin=USER_SNR_MIN, smax=USER_SNR_MAX):
    u = np.random.uniform(0.0, 1.0)
    a = smin**(-1.5)
    b = smax**(-1.5)
    s = (a - u * (a - b)) ** (-2.0 / 3.0)
    return s


def draw_euclidean_snr(smin=USER_SNR_MIN, smax=USER_SNR_MAX):
    # Decided not to use this because creates a huge bump at smax
    u = np.random.uniform(0.0, 1.0)
    s = smin * u**(-2.0/3.0)
    if s > smax:
        s = smax
    return s


def snr_to_scale_clipped(snr):
    sc = snr / K_DEFAULT
    if sc < SCALE_MIN:
        sc = SCALE_MIN
    if sc > SCALE_MAX:
        sc = SCALE_MAX
    return sc


def ensure_slack():
    if not os.path.exists(SLACK_TOKEN_FILE):
        raise RuntimeError(f"Could not find file with slack api token at {SLACK_TOKEN_FILE}")
    with open(SLACK_TOKEN_FILE, "r") as sf_handler:
        slack_token = sf_handler.read().strip()
    return slack.WebClient(token=slack_token)


def slack_msg(cli, text):
    if DRY_RUN or cli is None:
        print(f"[dry-run][slack] {text}")
        return
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
    params = np.genfromtxt(PARAMS_TXT)
    n = min(len(params), len(TEMPLATES))
    if n == 0:
        raise RuntimeError("No params/templates found")
    if n < len(params):
        print(f"[warn] Truncating params to {n} to match templates")
    return params[:n], TEMPLATES[:n]


def run():
    random.seed()
    np.random.seed()

    if not DRY_RUN:
        ensure_injection_list(INJ_LIST)
        ensure_injection_list(AUDIT_INJECTION_FILENAME)

    slack_cli = None
    if not DRY_RUN:
        slack_cli = ensure_slack()

    store = ds.DsaStore()

    

    params, templates = load_params_and_templates()

    print(f"[cfg] USER_SNR_MIN={USER_SNR_MIN}, USER_SNR_MAX={USER_SNR_MAX}")
    print(f"[cfg] SCALE_MIN={SCALE_MIN:.4f}, SCALE_MAX={SCALE_MAX:.4f}, K_DEFAULT={K_DEFAULT}")
    print(f"[cfg] DRY_RUN={DRY_RUN}")

    for i in range(niterations):
        print(f"[info] Starting iteration {i+1}/{niterations}")
        # Get current EW/NS stddevs
        beam_stats = store.get_dict("/mon/corr/101")
        stddev_ew = beam_stats["stddev_ew"]
        stddev_ns = beam_stats["stddev_ns"]

        # NS should get a smaller scale so that SNR is equalized
        ns_scale_suppression_factor = stddev_ns / stddev_ew
        # east-west beams have a higher stddev. Just the way the beamformer was setup.
        for idx, row in enumerate(params):
            DM, SNR_file, width_fwhm, spec_ind = map(float, row)

            # 1) draw target SNR from Euclidean in the user range
            effective_snr = draw_truncated_euclidean_snr()

            # 2) map to scale with clipping (EW reference)
            ew_scale = snr_to_scale_clipped(effective_snr)
            ns_scale = ew_scale * ns_scale_suppression_factor

            template = templates[idx]
            frbno = os.path.splitext(os.path.basename(template))[0].split("_")[-1]

            # draw an EW beam (0–255) and derive NS partner
            ew_beam = random.randint(0, 255)
            ns_beam = ew_beam + 256

            # node/local for each arm from global beam indices
            ew_node, ew_local = global_to_node_local(ew_beam)
            ns_node, ns_local = global_to_node_local(ns_beam)

            # assertions for mapping consistency
            expected_ns_node = ew_node + EW_TO_NS_NODE_SEP_FACTOR
            assert expected_ns_node == ns_node, (
                f"ns_node mismatch: ew_node={ew_node}, "
                f"expected ns_node={expected_ns_node}, got ns_node={ns_node}"
            )
            assert ew_local == ns_local, (
                f"local index mismatch: ew_local={ew_local}, ns_local={ns_local}"
            )

            mjd = Time.now().mjd

            # Log meta info to console always
            print(
                f"[inj] MJD={mjd:.5f} FRB={frbno} DM={DM:.1f} "
                f"target_SNR={effective_snr:.1f} "
                f"EW: beam={ew_beam}, ew_scale={ew_scale:.3f}, "
                f"NS: beam={ns_beam}, ns_scale={ns_scale:.3f}"
            )

            # Append to injection lists only if not dry-run
            if not DRY_RUN:
                with open(INJ_LIST, "a") as f:
                    f.write(FMT_OUT % (mjd, ew_beam, DM, effective_snr,
                                       width_fwhm, spec_ind, frbno))
                with open(AUDIT_INJECTION_FILENAME, "a") as f:
                    f.write(FMT_OUT % (mjd, ew_beam, DM, effective_snr,
                                       width_fwhm, spec_ind, frbno))

            # Build injection commands
            ew_injection_cmd = {
                "cmd": "inject",
                "val": f"{ew_local}-{template}-{ew_scale:.3f}-"
            }
            ns_injection_cmd = {
                "cmd": "inject",
                "val": f"{ns_local}-{template}-{ns_scale:.3f}-"
            }

            # In dry-run mode, just print the commands
            if DRY_RUN:
                print(f"[dry-run][cmd] EW node {ew_node}: {ew_injection_cmd}")
                print(f"[dry-run][cmd] NS node {ns_node}: {ns_injection_cmd}")
            else:
                # EW node (17/18) gets EW scale, NS node (19/20) gets NS scale
                store.put_dict(f"/cmd/corr/{ew_node}", ew_injection_cmd)
                store.put_dict(f"/cmd/corr/{ns_node}", ns_injection_cmd)

            # Slack notification (or dry-run echo)
            slack_msg(
                slack_cli,
                (
                    f"Sending injection to EW-beam {ew_beam} (scale={ew_scale:.2f}; stddev={stddev_ew:.1f}) and "
                    f"NS-beam {ns_beam} (scale={ns_scale:.2f}; stddev={stddev_ns:.1f}) with DM={DM:.1f} and SNR={effective_snr:.1f}"
                )
            )

            if not DRY_RUN:
                time.sleep(SLEEP_SEC)


if __name__ == "__main__":
    run()
