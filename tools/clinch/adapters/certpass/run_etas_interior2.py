#!/usr/bin/env python3
"""run_etas_interior2 — profile bisection on log10_tau for the hi ETAS
interior optimum. The profile NLL along the taper coordinate falls from
the banked point (boundary gradient negative there) and rises at the fit
bound (+1.1e-5 measured) — a single interior minimum whose quasi-flat
curvature (~2.5e-8) defeats plain Newton over the ridge. Bisect the SIGN
of the full-system dNLL/dlog10_tau at per-trial conditional optima of
the other 8 coordinates, then polish free and certify the full 9-dim
system.

Usage: python3 run_etas_interior2.py q3_hi [q0_hi]
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


def g_ltau_at(twin, oracle, th):
    sv, sg, _ = oracle.sources_f64(th, order=1)
    pv, pg, _ = twin.pairs_val_grad_hess(th)
    return float((pg + sg)[5]), pg + sg


def run(tag):
    log = log_for(tag)
    t0 = time.time()
    z = dict(np.load(os.path.join(DATA, f"etas_{tag}.npz")))
    cc = json.load(open(os.path.join(RECEIPTS, f"etas_{tag}_corner",
                                     "CERT.json")))
    th = np.asarray(cc["polished_center"], float)
    thb = np.r_[np.asarray(z["theta"], float)[0],
                np.asarray(z["theta"], float)[2:]]
    part = Part(blocks=[list(range(1, 9))], block_names=["triggering"],
                border_idx=[0])
    oracle = EtasOracle(z, part, tag, nworkers=NW)
    twin = TwinEtas(z)
    lo, hi = 4.5, 12.26          # g(hi) > 0 measured; g(lo) < 0 expected
    bis = []
    for it in range(22):
        mid = 0.5 * (lo + hi)
        th_try, _ = twin.hybrid_polish(th, oracle, maxit=12, tol_g=3e-5,
                                       pinned={5: mid},
                                       log=lambda s: None)
        g5, _ = g_ltau_at(twin, oracle, th_try)
        bis.append(dict(ltau=mid, g_ltau=g5))
        log(f"bisect it{it}: ltau {mid:.5f} g_ltau {g5:+.3e}")
        if abs(g5) < 2e-7 or (hi - lo) < 2e-4:
            th = th_try
            break
        if g5 > 0:
            hi = mid
        else:
            lo = mid
        th = th_try
    th[5] = 0.5 * (lo + hi) if abs(bis[-1]["g_ltau"]) >= 2e-7 else th[5]
    log("free hybrid polish (full 9-dim) from the bisection point")
    th0, _ = twin.hybrid_polish(th, oracle, maxit=30, log=log)
    log(f"free endpoint log10_tau = {th0[5]:.5f}")
    if not (4.0 < th0[5] < 12.26):
        log("endpoint left the interior window — no product")
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
        candidate_provenance="interior stationary point located by "
                             "profile bisection on log10_tau (the "
                             "quasi-flat ridge defeats plain Newton; "
                             "receipted trace below)",
        bisection_trace=bis,
        fit_artifact=cc.get("fit_artifact"),
        shift_vs_banked=shift, polish=prcpt,
        polished_center=[float(v) for v in thp],
        corner_attempt=f"receipts/etas_{tag}_corner (KKT boundary "
                       "gradient STRICTLY POSITIVE, the interior "
                       "evidence; banked beside this)")
    out = os.path.join(RECEIPTS, f"etas_{tag}_polished")
    cert.write_cert(res, out, f"etas_{tag}_polished", gates=None,
                    extra=extra)
    log(f"interior: {res['verdict']} (ltau* {thp[5]:.5f}, dist "
        f"{shift['distance_inf']:.3g}, nll drop "
        f"{shift['nll_drop']:.4g}, wall {round(time.time()-t0,1)}s)")
    oracle.close()


if __name__ == "__main__":
    for tag in sys.argv[1:] or ["q3_hi", "q0_hi"]:
        run(tag)
    print("INTERIOR2 PASS DONE", flush=True)
