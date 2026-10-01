#!/usr/bin/env python3
"""recertify_taupinned — the reduced pinned-taper product for the ETAS
targets whose corner KKT leg is DEGENERATE (the taper multiplier is
exact-strictly-negative at the center but smaller than the gradient's
variation over any admissible certificate box — the measured form of
ETAS taper non-identifiability).

Product: engine.certify on the MASKED oracle (log10_tau pinned at the
fit's own bound 12.26): a certificate that the remaining 8 coordinates'
optimum is a strict local minimum of the pinned-taper system. The
KKT-degeneracy numbers ride the receipt; the corner attempt's honest
refusal stays banked beside it.

Usage: python3 recertify_taupinned.py q3_lo [q0_lo ...]
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
from clinch.mask import MaskedOracle               # noqa
from oracle_etas import EtasOracle, Part           # noqa
from baller.hygiene import ctx_guard               # noqa
import baller.certify as bc                        # noqa

RECEIPTS = pkg.RECEIPTS_DIR
DATA = pkg.DATA_DIR
NW = int(os.environ.get("CERTPASS_NW", "10"))
LTAU_BOUND = 12.26


def log_for(tag):
    def log(s):
        print(f"[{tag}] {s}", flush=True)
    return log


def run(tag):
    log = log_for(tag)
    z = dict(np.load(os.path.join(DATA, f"etas_{tag}.npz")))
    cpath = os.path.join(RECEIPTS, f"etas_{tag}_corner", "CERT.json")
    cc = json.load(open(cpath))
    th = np.asarray(cc["polished_center"], float)
    assert th[5] == LTAU_BOUND
    part = Part(blocks=[list(range(1, 9))], block_names=["triggering"],
                border_idx=[0])
    oracle = EtasOracle(z, part, tag, nworkers=NW)
    masked = MaskedOracle(oracle, {5: LTAU_BOUND})
    # exact-center full-gradient receipt at the pin (the degeneracy row)
    bk = bc.block_krawczyk
    mp_ = masked.part
    with ctx_guard(prec=192):
        box0 = bk.box_around(
            ([[float(th[t]) for t in mp_.blocks[0]]], [float(th[0])]),
            0.0)
        Fz, Fg = masked.F_full(box0)
        g5 = Fz[0][4]          # base part: block=[1..8], ltau at pos 4
    kkt_deg = dict(
        grad_ltau_center_mid=float(g5.mid()),
        grad_ltau_center_rad=float(g5.rad()),
        strictly_negative_at_center=bool(g5 < 0),
        corner_kkt_enclosure=cc["reason"],
        note="the multiplier is EXACT-register strictly negative at the "
             "center but smaller than the gradient's variation over the "
             "certificate box — strict complementarity unprovable at any "
             "admissible box: the measured taper non-identifiability")
    log(f"g_ltau at pin: {kkt_deg['grad_ltau_center_mid']:.3e} "
        f"(rad {kkt_deg['grad_ltau_center_rad']:.1e}, "
        f"strictly<0: {kkt_deg['strictly_negative_at_center']})")
    log("certify: reduced pinned-taper system")
    res = engine.certify(masked, th, budget_s=10800, log=log)
    extra = dict(cc.get("shift_vs_banked") and
                 {"shift_vs_banked": cc["shift_vs_banked"]} or {},
                 pin={"log10_tau": LTAU_BOUND,
                      "register": "the fit's own L-BFGS-B upper bound "
                                  "(etas.inversion.RANGES)"},
                 kkt_degeneracy=kkt_deg,
                 corner_attempt="receipts/etas_%s_corner (REFUSED "
                                "kkt:coord_5, banked beside this)" % tag,
                 polished_center=[float(v) for v in th],
                 free_walk_diagnostic=cc.get("free_walk_diagnostic"))
    out = os.path.join(RECEIPTS, f"etas_{tag}_taupinned")
    cert.write_cert(res, out, f"etas_{tag}_taupinned (8-coord reduced "
                    "system, log10_tau pinned at the fit bound "
                    f"{LTAU_BOUND})", gates=None, extra=extra)
    log(f"taupinned: {res['verdict']} margin="
        f"{(res.get('krawczyk') or {}).get('max_margin')}")
    oracle.close()
    return res["verdict"]


if __name__ == "__main__":
    t0 = time.time()
    for tag in sys.argv[1:] or ["q3_lo", "q0_lo"]:
        run(tag)
    print("taupinned pass done in", round(time.time() - t0, 1), "s",
          flush=True)
