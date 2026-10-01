"""MIXALOT worked examples — large g, large k, large N, and exact-vs-MC.

Runnable as shipped, from ANY clean scratch cwd, printing live exact
computations at pilot-sized points, plus one live exact-vs-MC comparison.
Print-only: writes nothing anywhere.

Run: python3 examples/worked_examples.py   (from the tools/mixalot directory)
"""
import os
import random
import sys
import time
from fractions import Fraction as F
from math import log

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import mixalot  # noqa: E402
from mixalot.engines import load  # noqa: E402


def lnfrac(z):
    n, d = z.numerator, z.denominator
    shn = max(0, n.bit_length() - 500)
    shd = max(0, d.bit_length() - 500)
    return (log(n >> shn) + shn * log(2)) - (log(d >> shd) + shd * log(2))


def main():
    mixalot.verify(quiet=True)
    print("=" * 72)
    print("MIXALOT worked examples — exact computations live (sha-pinned")
    print("vendored engines), plus one live exact-vs-MC comparison")
    print("=" * 72)

    # ---------------- 1. LARGE N (1-var closed form / series route)
    print("\n[1] LARGE N — exact evidence in sub-second time")
    cf = load("closed_form_1var")
    for N in (1000, 100000):
        u = N // 2
        t0 = time.time()
        z = cf.Z_closed(u, u)
        print(f"    N={N:>7,}: exact logZ = {lnfrac(z):.6f}   "
              f"({time.time()-t0:.3f}s, exact rational held)")
    print("    -> the exact column is the referee: a sampler that drifts at "
          "large N is CAUGHT by it, never by itself (section 5 runs one "
          "such comparison live).")

    # ---------------- 2. LARGE K (bigk2 route; safe output handling)
    print("\n[2] LARGE K — thousands of states (no sampler reaches this)")
    bigk2 = load("bigk2")
    for k, N in ((100, 10000), (300, 30000)):
        random.seed(42)
        cuts = sorted(random.sample(range(1, N + k), k - 1))
        U = [b - a - 1 for a, b in zip([0] + cuts, cuts + [N + k])]
        t0 = time.time()
        Zv = bigk2.Z_fast2(U)
        # NOTE (avoids the 4300-digit repr crash): NEVER print the raw
        # Fraction at this scale — report lnZ; the exact rational stays
        # in memory for exact downstream arithmetic.
        print(f"    k={k:>5} N={N:>7,}: lnZ = {lnfrac(Zv):.6f}   "
              f"({time.time()-t0:.2f}s; numerator ~"
              f"{Zv.numerator.bit_length():,} bits — print lnZ, not repr())")

    # ---------------- 3. LARGE G + the exact DP limit
    print("\n[3] LARGE G — any component count, and the exact DP limit")
    bg = load("bigg")
    for g in (2, 3, 4):
        print(f"    Z((2,2,2), g={g}) = {bg.Z_bigg([2, 2, 2], g)}")
    zdpm = mixalot.core.Zdpm([2, 1], 1)     # front door: exact coercion
    gf = load("w4_blind_gf")
    print(f"    Z_DPM((2,1), alpha=1) = {zdpm}  (exact Fraction)")
    for g in (64, 256):
        gap = g * (gf.Zg_gf([2, 1], F(1), g) - zdpm)
        print(f"    g*(Z_g - Z_DPM) at g={g}: {gap}  (the exact -1/48 law)")
    print("    float route to g=1e6 components: engines.swap_route (k=2, "
          "certified vs exact at N=100 — a route demonstration; see MANUAL)")

    # ---------------- 4. The LSX flagship route (runnable)
    print("\n[4] LSX flagship (coin10) through the runnable front door:")
    ok = mixalot.core.lsx_example("coin10")
    print(f"    mixalot.core.lsx_example('coin10') -> PASS={ok}")

    # ---------------- 5. Exact-vs-MC, live at a small, quick point
    print("\n[5] EXACT vs MC, live at one small, quick point (m1, k=2, "
          "U=(50,50)):")
    truth = lnfrac(cf.Z_closed(50, 50))
    est = load("estimators")
    r = est.run("m1", 2, [50, 50], "nested", seed=1)
    z = (r["logZ_hat"] - truth) / r["err_est"]
    print(f"    exact truth  : {truth:.11f}   (closed form, instant)")
    print(f"    nested seed=1: {r['logZ_hat']:.11f} +- {r['err_est']:.4f} "
          f"(z={z:+.2f}, wall {r['wall_s']}s)")
    print("    at THIS N the sampler is calibrated — at larger N that can "
          "stop being true, and the exact number never moves.")
    print("\nDone. (No files were written.)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
