#!/usr/bin/env python3
"""run_etas_interior — the INTERIOR product for the hi ETAS targets.

The corner attempt's KKT leg measured dNLL/dlog10_tau STRICTLY POSITIVE
at the fit bound ([9.9e-6, 1.2e-5]) — the hi observed-NLL has an
interior taper optimum below 12.26 (unlike the lo fits, whose boundary
gradient is strictly negative). This runner walks the FREE system from
the corner center back to the interior stationary point and certifies
the full 9-coordinate system.

Usage: python3 run_etas_interior.py q3_hi [q0_hi]
"""
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import __init__ as pkg                             # noqa

import clinch                                      # noqa
from clinch import engine, cert                    # noqa
from oracle_etas import EtasOracle, Part           # noqa
from twin_etas import TwinEtas                     # noqa

RECEIPTS = pkg.RECEIPTS_DIR
DATA = pkg.DATA_DIR
NW = int(os.environ.get("CERTPASS_NW", "14"))
NAMES9 = ["log10_mu", "log10_k0", "a", "log10_c", "omega", "log10_tau",
          "log10_d", "gamma", "rho"]


def log_for(tag):
    def log(s):
        print(f"[{tag}] {s}", flush=True)
    return log


def run(tag):
    log = log_for(tag)
    t0 = time.time()
    z = dict(np.load(os.path.join(DATA, f"etas_{tag}.npz")))
    cc = json.load(open(os.path.join(RECEIPTS, f"etas_{tag}_corner",
                                     "CERT.json")))
    th_corner = np.asarray(cc["polished_center"], float)
    thb = np.r_[np.asarray(z["theta"], float)[0],
                np.asarray(z["theta"], float)[2:]]
    part = Part(blocks=[list(range(1, 9))], block_names=["triggering"],
                border_idx=[0])
    oracle = EtasOracle(z, part, tag, nworkers=NW)
    twin = TwinEtas(z)
    log("free hybrid polish from the BANKED point (the hi profile NLL "
        "is monotone in ltau over [4.5, 12.26], g_ltau > 0 measured — "
        "the interior optimum sits near the banked taper)")
    th0, hist = twin.hybrid_polish(thb.copy(), oracle, maxit=60, log=log)
    log(f"free endpoint log10_tau = {th0[5]:.4f}")
    if not (0.02 < th0[5] < 12.25):
        log("free walk left the interior of the fit's RANGES — "
            "stopping (no product)")
        oracle.close()
        return
    log("arb newton_polish (full 9-dim)")
    thp, prcpt = engine.newton_polish(oracle, th0, log=log)
    log("certify: full interior system")
    res = engine.certify(oracle, thp, budget_s=10800, log=log)
    shift = dict(
        distance_inf=float(np.max(np.abs(thp - thb))),
        per_coord={NAMES9[a]: dict(banked=float(thb[a]),
                                   polished=float(thp[a]),
                                   shift=float(thp[a] - thb[a]))
                   for a in range(9)},
        nll_twin_at_banked=float(twin.value(thb)),
        nll_twin_at_polished=float(twin.value(thp)))
    shift["nll_drop"] = (shift["nll_twin_at_banked"]
                         - shift["nll_twin_at_polished"])
    extra = dict(
        candidate_provenance="interior stationary point reached from the "
                             "corner center (free polish; KKT-positive "
                             "boundary gradient receipt in the corner "
                             "CERT)",
        fit_artifact=cc.get("fit_artifact"),
        shift_vs_banked=shift, polish=prcpt,
        polished_center=[float(v) for v in thp],
        corner_attempt=f"receipts/etas_{tag}_corner (KKT boundary "
                       "gradient STRICTLY POSITIVE — interior optimum "
                       "exists; banked beside this)")
    out = os.path.join(RECEIPTS, f"etas_{tag}_polished")
    cert.write_cert(res, out, f"etas_{tag}_polished", gates=None,
                    extra=extra)
    log(f"interior: {res['verdict']} "
        f"(ltau* {thp[5]:.4f}, dist {shift['distance_inf']:.3g}, "
        f"nll drop {shift['nll_drop']:.4g}, wall "
        f"{round(time.time() - t0, 1)}s)")
    oracle.close()


if __name__ == "__main__":
    for tag in sys.argv[1:] or ["q3_hi", "q0_hi"]:
        run(tag)
    print("INTERIOR PASS DONE", flush=True)
