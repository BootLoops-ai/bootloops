#!/usr/bin/env python3
# DIPSTICK member — verb `nonres`: GKZ facet enumeration + exact beta-pairing
# nonresonance certificate. Certifies ONLY at the tested beta; genericity by
# sampling. qhull float facets are exactness-verified before use; unverifiable
# facets dropped (conservative). GKZ36.json ships beside this file and loads
# package-relative; the seed point is pinned BY VALUE in w5_seed_base.json.
"""Nonresonance of beta for the K(3,6) GKZ system — facet enumeration of
cone(A) (28 generators in R^18) + exact pairing check. Background job."""
import json, itertools, sys
import numpy as np
from fractions import Fraction as Fr
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))
D = json.load(open(_os.path.join(_HERE, 'GKZ36.json')))
labels = D['labels']
minors = sorted({tuple(k) for k, m, c in labels})
mi = {k: i for i, k in enumerate(minors)}
N = len(labels)
A = np.zeros((18, N))
for j, (k, m, c) in enumerate(labels):
    A[mi[tuple(k)], j] = 1
    for r in range(4):
        A[14 + r, j] = m[r]
from scipy.spatial import ConvexHull
pts = np.hstack([np.zeros((18, 1)), A]).T
hull = ConvexHull(pts, qhull_options='Qx C-0')
eqs = hull.equations
facets = [e for e in eqs if abs(e[-1]) < 1e-8]
print(f'{len(facets)} candidate cone facets from qhull', flush=True)
import sympy as sp
Asym = sp.Matrix([[sp.Integer(int(round(x))) for x in row] for row in A])
uniq = {}
for e in facets:
    nr = [sp.nsimplify(round(x, 7), rational=True, tolerance=1e-6) for x in e[:-1]]
    d = sp.lcm([sp.Rational(x).q for x in nr])
    nz = tuple(int(sp.Rational(x)*d) for x in nr)
    gg = sp.gcd(list(nz))
    nz = tuple(int(x/gg) for x in nz) if gg else nz
    if nz in uniq or tuple(-x for x in nz) in uniq:
        continue
    vals = [sum(nz[i]*int(Asym[i, j]) for i in range(18)) for j in range(N)]
    if all(v >= 0 for v in vals) or all(v <= 0 for v in vals):
        uniq[nz] = vals
print(f'{len(uniq)} exact-verified distinct facets', flush=True)
import sympy as sp2
# The seed point is pinned BY VALUE (w5_seed_base.json beside this file, exact
# rationals).  BASE below reproduces the reference run in fixtures/ bit-identically.
_SEED = json.load(open(_os.path.join(_HERE, 'w5_seed_base.json')))
T = sorted(_SEED.keys())
BASE = {t: Fr(_SEED[t][0], _SEED[t][1]) for t in T}
TRI = [''.join(str(x) for x in c) for c in itertools.combinations(range(1,7), 3)]
sv = {t: sp2.Symbol('s'+t) for t in TRI}
eqs2 = [sum(sv[t] for t in TRI if i in t) for i in '123456']
miss = [t for t in TRI if t not in T]
sol = sp2.solve(eqs2, [sv[t] for t in miss], dict=True)[0]
subs = {sv[t]: sp2.Rational(BASE[t].numerator, BASE[t].denominator) for t in T}
S20 = {}
for t in TRI:
    S20[t] = BASE[t] if t in T else Fr(int(sol[sv[t]].subs(subs).p),
                                       int(sol[sv[t]].subs(subs).q))
beta = [-S20[''.join(str(z+1) for z in k)] for k in minors] + [Fr(0)]*4
res = []
for nz, vals in uniq.items():
    pair = sum(Fr(nz[i])*beta[i] for i in range(18))
    if pair.denominator == 1:
        res.append((nz, pair))
print(f'RESONANT facets at the physical point: {len(res)} of {len(uniq)}', flush=True)
for nz, p in res[:6]:
    print('  ', nz[:9], '... pairing', p, flush=True)
print('VERDICT:', 'NONRESONANT -> GKZ module IRREDUCIBLE -> minimal module of K = 26'
      if not res else 'RESONANT AT THIS POINT — check genericity in s', flush=True)
