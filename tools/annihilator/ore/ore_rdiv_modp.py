#!/usr/bin/env python3
"""
ore_rdiv_modp.py — Ore-algebra pseudo right-division over GF(p)[x].

Given L, R with ord(L) ≥ ord(R), compute (Q, S, e) with
    lc(R)^e · L  =  Q · R  +  S,     ord(S) < ord(R),
all with polynomial coefficients in GF(p)[x].
If S ≡ 0, R right-divides L over GF(p)(x) and Q is (a poly-multiple of) L//R.
"""
import flint
from math import comb


def _to_nm(op_il, p):
    return [flint.nmod_poly([int(c) for c in ci], p) for ci in op_il]


def _ord(L):
    for i in range(len(L) - 1, -1, -1):
        if L[i] != 0:
            return i
    return -1


def _Dk_times(R, k, p):
    """D^k · R  as an operator (Leibniz):  (D^k·R)_m = Σ_{s} C(k,s) D^s R_{m-k+s}."""
    r = len(R) - 1
    # precompute D^s R_j
    dR = [[flint.nmod_poly(list(R[j]), p) for j in range(r + 1)]]
    for _ in range(k):
        dR.append([c.derivative() for c in dR[-1]])
    out = [flint.nmod_poly([0], p) for _ in range(r + k + 1)]
    for m in range(r + k + 1):
        for s in range(k + 1):
            j = m - k + s
            if 0 <= j <= r:
                out[m] += dR[s][j] * comb(k, s)
    return out


def ore_rightdiv_modp(L_il, R_il, p, verbose=False):
    """Pseudo right-division.  Returns (Q_il, S_il, e, content_removed)."""
    L = _to_nm(L_il, p)
    R = _to_nm(R_il, p)
    n = _ord(L); r = _ord(R)
    if n < r:
        return [[0]], L_il, 0, None
    lcR = R[r]
    Q = [flint.nmod_poly([0], p) for _ in range(n - r + 1)]
    e = 0
    cur = [flint.nmod_poly(list(c), p) for c in L]
    for k in range(n - r, -1, -1):
        top = cur[r + k] if r + k < len(cur) else flint.nmod_poly([0], p)
        if top == 0:
            continue
        # multiply cur and Q by lcR, then subtract top · (D^k · R)
        cur = [c * lcR for c in cur]
        Q = [q * lcR for q in Q]
        e += 1
        DkR = _Dk_times(R, k, p)
        for m in range(len(DkR)):
            cur[m] -= top * DkR[m]
        Q[k] += top
        if verbose:
            print(f"  k={k}: ord(cur)={_ord(cur)} deg(top)={top.degree()}")
        assert cur[r + k] == 0, "top not killed"
    # remainder = cur[0..r-1]
    S = cur[:r]
    # content removal on Q
    g = Q[0]
    for q in Q[1:]:
        g = g.gcd(q)
    content = g
    if g.degree() >= 0 and g != 0:
        Q = [q // g for q in Q]
    return ([_il(q) for q in Q], [_il(s) for s in S], e,
            _il(content) if content != 0 else [0])


def _il(nm):
    d = nm.degree()
    return [int(nm[j]) for j in range(d + 1)] if d >= 0 else [0]


if __name__ == '__main__':
    # Synthetic planted-factor control, per the package law: validate
    # right-division on a known-good factor FIRST (S=0) before trusting
    # S!=0 negatives.
    import sys, os
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from lclm_modp import op_to_nmod, op_from_nmod, op_mul_nm
    p = 1048573
    # Plant L = Qf·R (Ore product), so R right-divides L exactly.
    R_il = [[1, 2, 3], [0, 1], [5, 0, 1]]           # order 2
    Qf_il = [[2, 1], [1, 1, 1], [3], [0, 0, 1]]     # order 3
    L_il = op_from_nmod(op_mul_nm(op_to_nmod(Qf_il, p), op_to_nmod(R_il, p), p))
    Q, S, e, cont = ore_rightdiv_modp(L_il, R_il, p)
    s_zero = not any(any(row) for row in S)
    print(f"S=0 control (planted factor): S vanishes = {s_zero} "
          f"(Q order={len(Q)-1}, e={e})")
    assert s_zero, "planted right factor must give S = 0"
    # Control: a perturbed operator must NOT right-divide L.
    R2_il = [[7, 2, 3], [0, 1], [5, 0, 1]]
    Q2, S2, e2, _ = ore_rightdiv_modp(L_il, R2_il, p)
    s2_nonzero = any(any(row) for row in S2)
    print(f"S!=0 control (non-factor): S nonzero = {s2_nonzero}")
    assert s2_nonzero, "non-factor must give S != 0"
    print("SELFTEST PASS")
