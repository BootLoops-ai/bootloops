# exact_polys.py -- P_x(u1,u2,v;p) exact trivariate polys, p-Taylor shift at
# pbar, Bernstein-basis conversion. Degrees (4,3,1) per char, p-deg 4.
import json
from math import comb
from fractions import Fraction as Fr

import os as _os
BASE = _os.environ.get("POSQ_STAGE0", _os.path.dirname(_os.path.abspath(__file__)))
# The production pattern-count data is NOT shipped; the polynomial machinery
# below needs no data file. Point POSQ_STAGE0 at a directory containing
# QUARTET_COLLAPSE.json to enable the production PATS table (else None).
_QC_PATH = f"{BASE}/QUARTET_COLLAPSE.json"
if _os.path.exists(_QC_PATH):
    QC = json.load(open(_QC_PATH))
    PATS = [(tuple(int(ch) for ch in p), int(n))
            for p, n in QC["pattern_counts"].items()]
    assert sum(n for _, n in PATS) == 318 and len(PATS) == 14
else:
    QC, PATS = None, None
DU1, DU2, DV = 4, 3, 1                   # per-char degrees in (u1,u2,v)

def padd(a, b):
    d = {m: dict(c) for m, c in a.items()}
    for m, c in b.items():
        t = d.setdefault(m, {})
        for k, v in c.items():
            t[k] = t.get(k, Fr(0)) + v
    return d

def pmul(a, b):
    d = {}
    for m1, c1 in a.items():
        for m2, c2 in b.items():
            m = (m1[0] + m2[0], m1[1] + m2[1], m1[2] + m2[2])
            t = d.setdefault(m, {})
            for k1, v1 in c1.items():
                for k2, v2 in c2.items():
                    t[k1 + k2] = t.get(k1 + k2, Fr(0)) + v1 * v2
    return d

def pconst(cp):
    return {(0, 0, 0): dict(cp)}

def scalec(poly, cp):
    d = {}
    for m, c in poly.items():
        t = {}
        for k1, v1 in c.items():
            for k2, v2 in cp.items():
                t[k1 + k2] = t.get(k1 + k2, Fr(0)) + v1 * v2
        d[m] = t
    return d

U1 = {(1, 0, 0): {0: Fr(1)}}; U2 = {(0, 1, 0): {0: Fr(1)}}; U3 = {(0, 0, 1): {0: Fr(1)}}
PI1 = {1: Fr(1)}; PI0 = {0: Fr(1), 1: Fr(-1)}
PI = (PI0, PI1)

def trans(um):
    out = [[None, None], [None, None]]
    for i in (0, 1):
        for j in (0, 1):
            delta = {0: Fr(1)} if i == j else {0: Fr(0)}
            dmj = {k: delta.get(k, Fr(0)) - PI[j].get(k, Fr(0))
                   for k in set(delta) | set(PI[j])}
            out[i][j] = padd(pconst(PI[j]), scalec(um, dmj))
    return out

U12 = pmul(U1, U2); U123 = pmul(U12, U3)
P1 = trans(U1); P2 = trans(U2); P3 = trans(U3)
P12 = trans(U12); P123 = trans(U123)

def pat_poly_uv(x):
    """poly in (u1,u2,v) [v = u3^2, parity-checked], coeffs = p-poly dicts."""
    a, b, c, d = x
    tot = {}
    for r in (0, 1):
        for w in (0, 1):
            for v in (0, 1):
                t = pmul(pmul(pmul(P3[r][w], P2[w][v]), pmul(P1[v][a], P1[v][b])),
                         pmul(P12[w][c], P123[r][d]))
                tot = padd(tot, scalec(t, PI[r]))
    pv = {}
    for (i, j, k), cp in tot.items():
        cp = {kk: v for kk, v in cp.items() if v != 0}
        if not cp:
            continue
        assert k % 2 == 0, ("u3-parity fail", x, (i, j, k))
        t = pv.setdefault((i, j, k // 2), {})
        for kk, v in cp.items():
            t[kk] = t.get(kk, Fr(0)) + v
    return pv

def taylor_shift(pv, pbar):
    """coeff p-polys -> t-polys, p = pbar + t; returns m -> {k: Fr}."""
    out = {}
    for m, cp in pv.items():
        t = {}
        for k, v in cp.items():
            for j in range(k + 1):
                t[j] = t.get(j, Fr(0)) + v * comb(k, j) * pbar ** (k - j)
        out[m] = {k: v for k, v in t.items() if v != 0}
    return out

def bk_power_grid(shifted, kmax=4):
    """b_k(u) power-basis coeff grid: list over k of dict m->Fr."""
    out = [dict() for _ in range(kmax + 1)]
    for m, tc in shifted.items():
        for k, v in tc.items():
            out[k][m] = v
    return out

def to_bernstein(grid):
    """power->Bernstein (deg 4,3,1 per var): c^B[i,j,l]."""
    cB = {}
    for (i, j, l) in [(a, b, c) for a in range(DU1 + 1)
                      for b in range(DU2 + 1) for c in range(DV + 1)]:
        s = Fr(0)
        for (a, b, c), v in grid.items():
            if a <= i and b <= j and c <= l:
                s += v * (Fr(comb(i, a), comb(DU1, a)) * Fr(comb(j, b), comb(DU2, b))
                          * Fr(comb(l, c), comb(DV, c)))
        if s != 0:
            cB[(i, j, l)] = s
    return cB

def eval_ppoly(pv, p, u1, u2, v):
    s = Fr(0)
    for (i, j, l), cp in pv.items():
        cval = sum(vv * p ** k for k, vv in cp.items())
        s += cval * u1 ** i * u2 ** j * v ** l
    return s
