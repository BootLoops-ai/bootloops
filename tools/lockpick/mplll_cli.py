#!/usr/bin/env python3
"""mplll_cli.py — CLI for mplll v2 (qrbkz/rawlll/rawbkz/graded/cvp, batch, dmin-sweep)."""
import argparse, json, os, sys
import mpmath as mp
from . import mplll, mplll_batch, mplll_graded  # g3 package-relative


def main():
    ap = argparse.ArgumentParser(description="multi-point LLL analytic-regression fitter")
    ap.add_argument("--basis"); ap.add_argument("--target")
    ap.add_argument("--dps", type=int, default=100)
    ap.add_argument("--height", type=int, default=10 ** 8)
    ap.add_argument("--held-out", default=""); ap.add_argument("--fit-points", default="")
    ap.add_argument("--method", default="qrbkz",
                    choices=["qrbkz", "rawlll", "rawbkz", "graded", "cvp"])
    ap.add_argument("--backend", default=None); ap.add_argument("--no-controls", action="store_true")
    ap.add_argument("--weights", default=None, help="FILE.json: list[int] or {name:w}")
    ap.add_argument("--ls-coeffs", default=None, help="cvp: FILE.json list of K reals")
    ap.add_argument("--denom", type=int, default=1, help="cvp: integer denominator D")
    ap.add_argument("--dmin-sweep", action="store_true")
    a = ap.parse_args()
    d = json.load(open(a.basis)); names = d["names"]
    F = [d["values"][n] for n in names]
    tgt = json.load(open(a.target))
    multi = isinstance(tgt, dict) and "targets" in tgt
    I = None if multi else (tgt.get("values", tgt) if isinstance(tgt, dict) else tgt)
    N = len(F[0]); ho = [int(x) for x in a.held_out.split(",") if x] or list(range(N))[-2:]
    fit = [int(x) for x in a.fit_points.split(",") if x] or [j for j in range(N) if j not in ho]
    mp.mp.dps = a.dps + 30
    if not a.no_controls and not multi:
        c = mplll.controls(I, F, names, fit, ho, a.dps, a.height, a.backend,
                           "qrbkz" if a.method in ("graded", "cvp") else a.method)
        if not c["ok"]:
            print(json.dumps({"status": "CONTROLS_FAILED", "detail": c}, default=str)); sys.exit(2)
    if multi:
        r = mplll_batch.qrbkz_batch(tgt["targets"], F, names, fit, ho, a.dps,
                                    a.height, backend=a.backend)
    elif a.dmin_sweep:
        r = mplll_batch.dmin_sweep(I, F, names, fit, ho, a.dps, method=a.method,
                                   backend=a.backend)
    elif a.method == "cvp":
        cls = json.load(open(a.ls_coeffs)) if a.ls_coeffs else None
        r = mplll_batch.cvp_fit(I, F, names, fit, ho, a.dps, a.denom,
                                c_LS=cls, backend=a.backend)
    elif a.method == "graded":
        w = None
        if a.weights:
            wd = json.load(open(a.weights))
            w = wd if isinstance(wd, list) else [wd[n] for n in names]
        r = mplll_graded.graded_fit(I, F, names, w, fit, ho, a.dps,
                                    hmax=a.height, backend=a.backend)
    else:
        r = mplll.mplll_fit(I, F, names, fit, ho, a.dps, a.height,
                            backend=a.backend, method=a.method)
    print(json.dumps(r, indent=1, default=str))
    st = r.get("status") if isinstance(r, dict) else None
    sys.exit(0 if st == "HIT" or multi or a.dmin_sweep else 1)


if __name__ == "__main__":
    main()
