# exact_polys_bal.py -- BALANCED-shape exact pattern polynomials (the GAUSS-
# COLLAPSE min/diff decomposition, PINNED_SPEC pin 3 product-domain convention).
#
# Wedge case "first cherry OLDER": cherry(L1,L2) at height m+d, cherry(L3,L4)
# at height m, root at h = m+d+s3; m = min(u1,u2) ~ Exp(2*lam), d = |u1-u2| ~
# Exp(lam), s3 ~ Exp(lam) (iid-Exp min/diff factorization, each ordering w.p.
# 1/2). Variables w1 = e^-beta*m, w2 = e^-beta*d, w3 = e^-beta*s3.
# Edge decays: leaves L1,L2: w1*w2; leaves L3,L4: w1; root->old: w3;
# root->young: w2*w3. Root sum under stationarity kills odd w3 powers
# (asserted), v = w3^2. Per-char degrees (w1,w2,v) = (4,3,1) (asserted);
# p-degree measured and exported as KMAX_BAL.
#
# Z_bal(p) = c^3 * (tot_case1 + tot_case2)/2,
#   tot_case = sum w1w2w3-rule-weighted F_case over the tensor rule with
#   u1-weight w1^(2c-1) [alpha = 2c-1], u2-weight w2^(c-1), v-weight v^(c/2-1)
#   -- SAME assembly form Z = c^3 * tot/2 as the caterpillar engine, since
#   prefactor 2c * c * (c/2) = c^3 and Fbar = (F1+F2)/2.
from fractions import Fraction as Fr
from exact_polys import (padd, pmul, pconst, scalec, trans, taylor_shift,
                         bk_power_grid, U1, U2, U3, PI, DU1, DU2, DV)
from math import comb

E1 = pmul(U1, U2)      # old-cherry leaf decay  w1*w2
E2 = U1                # young-cherry leaf decay w1
E3 = U3                # root -> old-cherry edge w3
E4 = pmul(U2, U3)      # root -> young-cherry edge w2*w3
PE1 = trans(E1); PE2 = trans(E2); PE3 = trans(E3); PE4 = trans(E4)

def bal_pat_poly(x):
    """poly in (w1,w2,v) [v = w3^2, parity-asserted]; coeffs = p-poly dicts.
    x = (x1,x2,x3,x4) in wedge slot order (old1, old2, young1, young2)."""
    x1, x2, x3, x4 = x
    tot = {}
    for r in (0, 1):
        for a in (0, 1):
            for b in (0, 1):
                t = pmul(pmul(pmul(PE3[r][a], PE1[a][x1]), pmul(PE1[a][x2],
                         PE4[r][b])), pmul(PE2[b][x3], PE2[b][x4]))
                tot = padd(tot, scalec(t, PI[r]))
    pv = {}
    for (i, j, k), cp in tot.items():
        cp = {kk: v for kk, v in cp.items() if v != 0}
        if not cp:
            continue
        assert k % 2 == 0, ("w3-parity fail", x, (i, j, k))
        t = pv.setdefault((i, j, k // 2), {})
        for kk, v in cp.items():
            t[kk] = t.get(kk, Fr(0)) + v
    return {m: {k: v for k, v in cp.items() if v != 0} for m, cp in pv.items()
            if any(v != 0 for v in cp.values())}

# degree audit on ALL 16 patterns + sum-to-one exact identity
_md = [0, 0, 0]
_pd = 0
_tot = {}
_ALL16 = [tuple(int(ch) for ch in f"{i:04b}") for i in range(16)]
for _x in _ALL16:
    _pv = bal_pat_poly(_x)
    for (_i, _j, _l), _cp in _pv.items():
        _md = [max(_md[0], _i), max(_md[1], _j), max(_md[2], _l)]
        _pd = max(_pd, max(_cp))
    for _m, _cp in _pv.items():
        _t = _tot.setdefault(_m, {})
        for _k, _v in _cp.items():
            _t[_k] = _t.get(_k, Fr(0)) + _v
assert _md == [DU1, DU2, DV], ("balanced per-char degrees", _md)
KMAX_BAL = _pd
_tot = {m: {k: v for k, v in cp.items() if v != 0} for m, cp in _tot.items()}
_tot = {m: cp for m, cp in _tot.items() if cp}
assert _tot == {(0, 0, 0): {0: Fr(1)}}, "balanced sum-to-one identity FAILED"

def to_bernstein_bal(grid):
    """power->Bernstein (deg 4,3,1 per var) -- same law as exact_polys but
    local so the balanced tables never depend on the caterpillar module's
    globals beyond the degree constants."""
    cB = {}
    for (i, j, l) in [(a, b, c) for a in range(DU1 + 1)
                      for b in range(DU2 + 1) for c in range(DV + 1)]:
        s = Fr(0)
        for (a, b, c), v in grid.items():
            if a <= i and b <= j and c <= l:
                s += v * (Fr(comb(i, a), comb(DU1, a))
                          * Fr(comb(j, b), comb(DU2, b))
                          * Fr(comb(l, c), comb(DV, c)))
        if s != 0:
            cB[(i, j, l)] = s
    return cB

def eval_bal_ppoly(pv, p, w1, w2, v):
    s = Fr(0)
    for (i, j, l), cp in pv.items():
        cval = sum(vv * p ** k for k, vv in cp.items())
        s += cval * w1 ** i * w2 ** j * v ** l
    return s
