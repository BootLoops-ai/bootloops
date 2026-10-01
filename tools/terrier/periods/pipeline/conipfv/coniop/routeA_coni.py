#!/usr/bin/env python3
"""routeA_coni.py — restricted PF operator L_s from the graded coni series.

Engine: tools/annihilator (annihilator.py) via the pipeline restrict_op helpers
(mod-p ode-order-first (r,s) pin at TWO primes x TWO truncations, multi-prime
CRT exact reconstruction, exact all-window annihilation proof with true
held-out tail >= 30%, theta-form + indicial factorization).
Usage: routeA_coni.py <card>.series.json [rmax] [smax]   (gseries.py output)
"""
import json, os, sys, time
from fractions import Fraction as Fr
_TOOLS = os.environ.get("TERRIER_TOOLS_ROOT") or os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), *[os.pardir] * 5))
sys.path.insert(0, os.path.join(_TOOLS, "annihilator"))
sys.path.insert(0, _TOOLS)
_PIPE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, _PIPE)
from restrict_op import (pin_rs, reconstruct, prove_annihilation,
                         indicial_factor)
from annihilator import rec_to_theta

HERE = os.path.dirname(os.path.abspath(__file__))


def run(series_file, rmax, smax, holdout_frac=0.30):
    d = json.load(open(series_file))
    A = [Fr(x) for x in d["a"]]
    M = len(A) - 1
    walls = {}
    t0 = time.time()
    r1, s1 = pin_rs(A, rmax, smax, p2=2147483629, trunc2=int(0.72 * len(A)))
    walls["pin_s"] = round(time.time() - t0, 1)
    print(f"[pin] (r,s)=({r1},{s1}) @p31 N={M+1} AND @p2147483629 "
          f"N={int(0.72*len(A))}   {walls['pin_s']}s", flush=True)
    nunk = (r1 + 1) * (s1 + 1)
    nfit = min(r1 + nunk + 40, int((M - r1) * (1 - holdout_frac)))
    assert nfit >= r1 + nunk + 25, "series too short for reconstruction"
    t0 = time.time()
    coeffs, bits = reconstruct(A, r1, s1, nfit, n_primes=max(
        6, (nunk * 110) // 62 // 30 + 2))
    walls["reconstruct_s"] = round(time.time() - t0, 1)
    print(f"[exact] {len(coeffs)} nonzero c_jk, CRT {bits} bits   "
          f"{walls['reconstruct_s']}s", flush=True)
    t0 = time.time()
    nwin, nheld = prove_annihilation(A, coeffs, r1, nfit)
    walls["prove_s"] = round(time.time() - t0, 1)
    frac = nheld / nwin
    print(f"[p1 GATE] PASS: all {nwin} windows exact; {nheld} held-out "
          f"({100*frac:.0f}%)   {walls['prove_s']}s", flush=True)
    assert frac >= holdout_frac - 0.02, "held-out fraction below gate"
    ode = rec_to_theta(r1, s1, coeffs)
    assert all(c.denominator == 1 for c in ode.values())
    lead, fac, mult0 = indicial_factor(ode)
    print(f"[indicial@s=0] lead {lead} factors {fac} mult0 {mult0}")
    out = {"card": d["card"], "series_M": M, "r_shift": r1, "s_theta": s1,
           "recurrence": {f"{j},{k}": str(c) for (j, k), c in coeffs.items()},
           "theta_form": {f"{i},{k}": str(c) for (i, k), c in ode.items()},
           "fit_window": nfit, "held_out_windows": nheld, "windows": nwin,
           "indicial_s0": {"lead": lead, "factors": fac, "mult0": mult0},
           "pins": {"p31_full": [r1, s1], "p2147483629_trunc": [r1, s1]},
           "walls_s": walls}
    fn = os.path.join(HERE, f"{d['card']}.coni_op.json")
    json.dump(out, open(fn, "w"), indent=1)
    print("wrote", fn, flush=True)
    return out


if __name__ == "__main__":
    run(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 40,
        int(sys.argv[3]) if len(sys.argv) > 3 else 10)
