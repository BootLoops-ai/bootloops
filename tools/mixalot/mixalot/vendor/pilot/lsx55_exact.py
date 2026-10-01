#!/usr/bin/env python3
"""lsx55_exact.py — reproduce LSX Example 5.5 (the 3x3 patients-with-
schizophrenia table of Evans, Gilula & Guttman (1989), N=132; their cost:
16 DAYS in 2009 Maple) exactly, via the mod-p + CRT +
balanced-rational-reconstruction route (route 2 of the Phase-1 gate design).

Gates: (a) ratrec STABLE across 3 consecutive prime-count checkpoints;
(b) verification at 2 HELD-OUT primes never used in the CRT;
(c) fraction equality vs the paper's printed 143-digit / 262-digit pair
    (transcribed; equality as FRACTIONS — LSX do not claim coprimality);
(d) float sanity ~ 2.2625e-119.
Discipline: phylo ratrec (balanced lift, B = isqrt(M//2)); workers niced.
"""

import json
import math
import multiprocessing as mp
import os
import sys
import time
from fractions import Fraction

# sibling imports resolve to the vendored copies beside this file
# (MIXALOT_PILOT_DIR overrides the default sibling path)
sys.path.insert(0, os.environ.get("MIXALOT_PILOT_DIR",
                                  os.path.dirname(os.path.abspath(__file__))))

LSX55_TABLE = [[43, 16, 3], [6, 11, 10], [9, 18, 16]]

PAPER_NUM = int(
    "278019488531063389120643600324989329103876140805"
    "285242839582092569357265886675322845874097528033"
    "99493069713103633199906939405711180837568853737")
PAPER_DEN = int(
    "12288402873591935400678094796599848745442833177572204"
    "50448819979286456995185542195946815073112429169997801"
    "33503900169921912167352239204153786645029153951176422"
    "43298328046163472261962028461650432024356339706541132"
    "34375318471880274818667657423749120000000000000000")


def gen_primes_25bit(count):
    ps, n = [], (1 << 25) - 1
    while len(ps) < count:
        if all(n % q for q in range(2, int(n ** 0.5) + 1)):
            ps.append(n)
        n -= 2
    return ps


def worker(p):
    os.nice(5)
    from zseries import Z_table_modp
    t0 = time.time()
    z = Z_table_modp(LSX55_TABLE, p)
    return p, z, round(time.time() - t0, 1)


def crt_list(rems, mods):
    R, M = rems[0], mods[0]
    for r, m in zip(rems[1:], mods[1:]):
        g, x, _ = ext_gcd(M, m)
        assert g == 1
        R = (R + (r - R) * x % m * M) % (M * m)
        M *= m
    return R, M


def ext_gcd(a, b):
    if b == 0:
        return a, 1, 0
    g, x, y = ext_gcd(b, a % b)
    return g, y, x - (a // b) * y


def ratrec(x, M, B_num=None):
    """Rational reconstruction. Evidence fractions are ASYMMETRIC (num << den
    since Z ~ 10^-119): symmetric isqrt bound wastes half the modulus — a
    farm 'instability' is exactly this symptom. Pass B_num ~ 2^(expected
    numerator bits + margin); den bound is implied M/(2*B_num). Default
    (B_num=None) = balanced isqrt(M//2)."""
    if B_num is None:
        B_num = math.isqrt(M // 2)
    den_bound = M // (2 * B_num)
    a0, a1 = M, x % M
    p0, p1 = 0, 1
    while a1 > B_num:
        q = a0 // a1
        a0, a1 = a1, a0 - q * a1
        p0, p1 = p1, p0 - q * p1
    if abs(p1) > den_bound or math.gcd(a1, abs(p1)) != 1 or p1 == 0:
        return None
    num, den = (a1, p1) if p1 > 0 else (-a1, -p1)
    return Fraction(num, den)


def main():
    nprimes, nheld, nworkers = 64, 2, 12
    primes = gen_primes_25bit(nprimes + nheld)
    held, use = primes[:nheld], primes[nheld:]
    t00 = time.time()
    with mp.Pool(nworkers) as pool:
        results = dict()
        for p, z, dt in pool.imap_unordered(worker, use + held):
            results[p] = z
            print(f"  p={p} z={z} ({dt}s) [{len(results)}/{len(primes)}]",
                  flush=True)
    wall = time.time() - t00

    rems = [results[p] for p in use]
    # stability checkpoints
    fracs = {}
    for cut in (nprimes - 8, nprimes - 4, nprimes):
        R, M = crt_list(rems[:cut], use[:cut])
        fracs[cut] = ratrec(R, M, B_num=1 << 520)
    stable = (fracs[nprimes - 8] is not None
              and fracs[nprimes - 8] == fracs[nprimes - 4] == fracs[nprimes])
    Z = fracs[nprimes]
    print(f"ratrec stable across {nprimes-8}/{nprimes-4}/{nprimes} primes: {stable}")
    assert stable, "ratrec not stable — need more primes"

    heldout_ok = all(
        Z.numerator % p * pow(Z.denominator % p, p - 2, p) % p == results[p]
        for p in held)
    print(f"held-out prime verification ({nheld} primes): "
          f"{'PASS' if heldout_ok else 'FAIL'}")

    paper = Fraction(PAPER_NUM, PAPER_DEN)
    paper_eq = (Z == paper)
    print(f"== paper 143/262-digit fraction (Ex 5.5): "
          f"{'PASS' if paper_eq else 'FAIL'}")
    import mpmath as mpm
    mpm.mp.dps = 30
    fz = mpm.mpf(Z.numerator) / mpm.mpf(Z.denominator)
    print(f"float: {mpm.nstr(fz, 10)} (paper ~2.2625e-119)")
    print(f"digits: num {len(str(Z.numerator))} den {len(str(Z.denominator))} "
          f"(paper printed pair 143/262, coprimality not claimed)")
    print(f"total wall {wall:.0f}s with {nworkers} niced workers "
          f"(LSX 2009: 16 days)")

    json.dump({"table": LSX55_TABLE, "num": str(Z.numerator),
               "den": str(Z.denominator), "stable": stable,
               "heldout_ok": heldout_ok, "paper_eq": paper_eq,
               "wall_s": round(wall, 1), "nprimes": nprimes,
               "primes_held": held},
              open("LSX55_EXACT_RESULT.json", "w"), indent=1)
    ok = stable and heldout_ok and paper_eq
    print("LSX55:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
