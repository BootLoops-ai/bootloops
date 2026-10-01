#!/usr/bin/env python3
"""Decisive integrand-identity test for the cone decomposition.

For each simplicial subcone (rays R, kappa, Q-tables built at the GENERIC
point) and random z in (0,inf)^4:
    y = exp(-z), w = R z (columns=rays... here rows->transpose), x = exp(w)
    claim:  prod_i y_i^{kappa_i} * prod_t Q_t(y)^{s_t}
         ==  prod_{all 14 nonconst minors} P_t(x)^{s_t}
with P_t the exact sympy minors. 30-digit check, all subcones, 2 points each.
This validates every Q-table entry, kappa, and the ray/orientation bookkeeping
(fan coverage and the measure factor |det R| are validated separately by the
20k-direction coverage test and the Gamma-form control gate).
"""
import mpmath as mp
import numpy as np
import json, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from run_pilot import generic_point
from cone_engine import subcone_tasks, FAN
from tropical_cones import POLY, PHAT_MONOS
import sympy as sp
import itertools

mp.mp.dps = 40

# exact minors from sympy (independent rebuild, same as setup_x36.py)
x11, x12, x21, x22 = sp.symbols('x11 x12 x21 x22', positive=True)
m33 = 1 + x11*(1+x12) + x21*(1+x12+x22)
M = sp.Matrix([[0,0,1,1,1,1],
               [0,-1,0,1,1+x11,1+x11*(1+x12)],
               [1,0,0,1,1+x11+x21,m33]])
NAMES = {}
for t in itertools.combinations(range(1,7), 3):
    d = sp.expand(M[:, [t[0]-1, t[1]-1, t[2]-1]].det())
    NAMES[''.join(map(str, t))] = d

s = generic_point()
tasks = subcone_tasks(s)
# rebuild ray lists in the same order used by subcone_tasks
ray_sets = []
for c in FAN:
    for simp in c['simplices']:
        ray_sets.append([c['rays'][i] for i in simp])
assert len(ray_sets) == len(tasks)

rng = np.random.default_rng(99)
worst = 0.0
for (task, rays) in zip(tasks, ray_sets):
    for rep in range(2):
        z = [mp.mpf(str(v)) for v in rng.uniform(0.2, 2.5, 4)]
        y = [mp.e**(-zi) for zi in z]
        w = [sum(mp.mpf(rays[i][j])*z[i] for i in range(4)) for j in range(4)]
        xv = {x11: mp.e**w[0], x12: mp.e**w[1],
              x21: mp.e**w[2], x22: mp.e**w[3]}
        lhs = mp.mpf(1)
        for i in range(4):
            k = Fr(task['kappa'][i])
            lhs *= y[i]**(mp.mpf(k.numerator)/k.denominator)
        for gs, monos in task['facs']:
            gam = Fr(gs)
            Q = mp.mpf(0)
            for e in monos:
                Q += y[0]**e[0]*y[1]**e[1]*y[2]**e[2]*y[3]**e[3]
            lhs *= Q**(mp.mpf(gam.numerator)/gam.denominator)
        rhs = mp.mpf(1)
        for t, stv in s.items():
            P = NAMES[t]
            Pv = mp.mpf(sp.srepr(P) and str(sp.N(P.subs(
                {x11: sp.Float(str(xv[x11]), 45), x12: sp.Float(str(xv[x12]), 45),
                 x21: sp.Float(str(xv[x21]), 45), x22: sp.Float(str(xv[x22]), 45)}), 45)))
            rhs *= Pv**(mp.mpf(stv.numerator)/stv.denominator)
        rel = abs(lhs-rhs)/abs(rhs)
        worst = max(worst, float(rel))
print(f"subcones tested: {len(tasks)} x2 points; worst rel deviation: {worst:.3e}")
print("IDENTITY:", "PASS" if worst < 1e-25 else "FAIL")
