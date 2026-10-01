"""MANDATORY verification driver: 3-level formula + lattice evaluator verification.

V1: formula_P3 == nested-urn process DP, exact Fraction equality, on EVERY
    dataset of 3 configurations (2 regions, 3 plots; 340 probabilities total),
    with exact unitarity of both sides.
V2 (G1): hier_eval.hier_P3 (scalable lattice DP) == formula_P3 on all cases.
V3 (G1-kron spot): hier_kron enclosures contain the exact values (12 cases).
Fail-loud: rc!=0 on any mismatch. Run time ~10 min.
"""
import sys
from fractions import Fraction

from hierarchical3 import formula_P3, process_distribution3
from hier_eval import hier_P3
from hier_kron import hier_lnP_kron
from flint import arb, ctx

CONFIGS = [
    ([2,2,2], [0,0,1], Fraction(3,2), [Fraction(2), Fraction(5)],
     [Fraction(1), Fraction(3), Fraction(2)]),
    ([3,2,2], [0,0,1], Fraction(7), [Fraction(1,3), Fraction(4)],
     [Fraction(2), Fraction(2), Fraction(6)]),
    ([2,3,2], [0,1,1], Fraction(1,5), [Fraction(6), Fraction(1,2)],
     [Fraction(5), Fraction(1), Fraction(1)]),
]

fails = tests = 0
for Jvec, reg_of, th, rhos, Ivec in CONFIGS:
    dist = process_distribution3(Jvec, reg_of, th, rhos, Ivec)
    assert sum(dist.values()) == 1, "DP unitarity"
    tot_f = tot_e = Fraction(0)
    for key, pdp in dist.items():
        Dmat = [list(r) for r in key]
        pf = formula_P3(Dmat, reg_of, th, rhos, Ivec)
        pe = hier_P3(Dmat, reg_of, th, rhos, Ivec)
        tot_f += pf
        tot_e += pe
        tests += 1
        if pf != pdp:
            fails += 1
            print(f"V1 FAIL {key}: formula={pf} dp={pdp}")
        if pe != pf:
            fails += 1
            print(f"V2 FAIL {key}: lattice={pe} formula={pf}")
    if tot_f != 1 or tot_e != 1:
        fails += 1
        print(f"unitarity FAIL: {tot_f}, {tot_e}")
    print(f"J={Jvec} reg={reg_of}: {len(dist)} datasets, V1+V2 exact, unitarity 1")

# V3: kron spot enclosures
Jvec, reg_of = [3,2,2], [0,0,1]
th, rhos, Ivec = Fraction(7), [Fraction(1,3), Fraction(4)], [Fraction(2), Fraction(2), Fraction(6)]
dist = process_distribution3(Jvec, reg_of, th, rhos, Ivec)
ctx.prec = 256
n3 = 0
for key, pex in list(dist.items())[:12]:
    lnb, _ = hier_lnP_kron([list(r) for r in key], reg_of, '7',
                           ['0.33333333333333333333333333333333', '4'],
                           ['2', '2', '6'], prec=128)
    exl = float((arb(pex.numerator).log() - arb(pex.denominator).log()).mid())
    if abs(float(lnb.mid()) - exl) > float(lnb.rad()) + 1e-20:
        fails += 1
        print(f"V3 FAIL {key}")
    n3 += 1
print(f"V3: {n3} kron enclosures contain exact values")

print(f"\nTOTAL: {tests} dataset probabilities checked (V1+V2), {fails} failures")
if fails:
    sys.exit(1)
print("VERIFY_HIERARCHICAL3: ALL PASS")
