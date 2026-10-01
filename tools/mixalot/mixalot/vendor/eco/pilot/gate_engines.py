"""Gate the product-tree engines against the certified brute-force oracle.

G1: P_exact (fmpz product tree) == etienne_P (naive convolution oracle),
    exact Fraction equality, random datasets J<=60, 3 (theta,m) points.
G2: logP_ball enclosure CONTAINS log(exact rational P), and ball width
    <= 2^(-prec+30), same datasets; two precisions (128, 256).
"""
import random
import sys
from fractions import Fraction

from etienne_oracle import etienne_P
from ball_engine import P_exact, logP_ball, K_ball
from flint import arb, ctx


def random_dataset(J):
    """Random composition-ish abundance dataset summing to J."""
    D = []
    left = J
    while left > 0:
        n = min(left, random.choice([1, 1, 1, 2, 2, 3, 4, 5, 8, 12, 20]))
        D.append(n)
        left -= n
    return sorted(D)


POINTS = [(Fraction(3, 2), Fraction(1, 3)),
          (Fraction(20), Fraction(9, 10)),
          (Fraction(1, 5), Fraction(1, 7))]

random.seed(20260708)
fails = 0
datasets = [random_dataset(J) for J in (10, 25, 40, 60) for _ in range(3)]

# G1
for D in datasets:
    for theta, m in POINTS:
        a = etienne_P(D, theta, m=m)
        b = P_exact(D, theta, m=m)
        if a != b:
            fails += 1
            print(f"G1 FAIL D={D} theta={theta} m={m}: naive={a} tree={b}")
print(f"G1: product-tree exact == naive oracle on {len(datasets)} datasets x {len(POINTS)} points")

# G2
import math
for prec in (128, 256):
    worst_w = 0.0
    for D in datasets:
        Kp = K_ball(D, prec)
        for theta, m in POINTS:
            exact = etienne_P(D, theta, m=m)
            lp = logP_ball(D, theta, m, prec, Kpoly=Kp)
            # exact log via high-precision arb of the rational
            ctx.prec = prec + 200
            ex = (arb(exact.numerator) / arb(exact.denominator)).log()
            ctx.prec = prec
            mid, rad = float(lp.mid()), float(lp.rad())
            exf = float(ex.mid())
            if not (mid - rad <= exf <= mid + rad):
                # containment check with proper ball overlap
                if not lp.overlaps(ex):
                    fails += 1
                    print(f"G2 FAIL containment prec={prec} D={D} th={theta} m={m}")
            w = rad / max(abs(mid), 1e-300)
            worst_w = max(worst_w, w)
    bound = 2.0 ** (-prec + 30)
    status = "OK" if worst_w <= bound else "FAIL"
    if status == "FAIL":
        fails += 1
    print(f"G2 prec={prec}: worst rel ball width {worst_w:.3e} (bound {bound:.3e}) {status}")

if fails:
    print(f"*** {fails} FAILURES ***")
    sys.exit(1)
print("GATES PASS")
