#!/usr/bin/env python3
"""
lclm_modp.py — Ore-algebra pairwise LCLM over GF(p)[x], via Sylvester nullspace
in nmod_poly.  Each stage verified by right-divisibility.
"""
import numpy as np
import flint


# ════════════════════════════════════════════════════════════════════════════
#  Operator representation:  list of nmod_poly, low→high derivative order.
# ════════════════════════════════════════════════════════════════════════════
def op_to_nmod(op_intlists, p):
    return [flint.nmod_poly([int(c) for c in ci], p) for ci in op_intlists]


def op_from_nmod(op_nm):
    return [[int(c[j]) for j in range(c.degree() + 1)] if c.degree() >= 0 else [0]
            for c in op_nm]


def op_order(A):
    for i in range(len(A) - 1, -1, -1):
        if A[i].degree() >= 0 or A[i] != 0:
            return i
    return 0


def op_primitive_nm(L, p):
    """Remove polynomial content (gcd of coefficient polys) from operator L."""
    g = L[0]
    for c in L[1:]:
        g = g.gcd(c)
    if g.degree() <= 0:
        return L
    return [c // g for c in L]


def op_mul_nm(L, R, p):
    """Ore product L·R:  (L·R)[f] = L[R[f]].  Leibniz on the right."""
    a, b = len(L) - 1, len(R) - 1
    # Precompute D^k(R_j) for k up to a.
    dR = [[flint.nmod_poly(list(R[j]), p) for j in range(b + 1)]]
    for k in range(1, a + 1):
        dR.append([c.derivative() for c in dR[-1]])
    from math import comb
    out = [flint.nmod_poly([0], p) for _ in range(a + b + 1)]
    # L·R = Σ_i L_i D^i · Σ_j R_j D^j = Σ_i L_i Σ_j Σ_{s≤i} C(i,s) D^s(R_j) D^{i-s+j}
    for i in range(a + 1):
        for j in range(b + 1):
            for s in range(i + 1):
                k = (i - s) + j
                out[k] += L[i] * dR[s][j] * comb(i, s)
    return out


def right_divides_nm(L, A, p, extra_deg=6):
    """Check if A right-divides L over GF(p)(x).  Exact: Ore pseudo-division
    (ore_rdiv_modp), True iff the remainder S vanishes identically.

    Divisibility over GF(p)(x) permits a RATIONAL cofactor (denominator
    lc(A)^e after content removal), so solving a linear system for a
    polynomial cofactor of bounded degree can false-negative on valid LCLMs;
    the Ore pseudo-division test is exact.
    `extra_deg` is retained for API compatibility and ignored.
    """
    try:
        from ore_rdiv_modp import ore_rightdiv_modp
    except ImportError:  # package-style import
        from .ore_rdiv_modp import ore_rightdiv_modp
    L_il = op_from_nmod(L); A_il = op_from_nmod(A)
    _, S, _, _ = ore_rightdiv_modp(L_il, A_il, p)
    return not any(any(row) for row in S)


def lclm2_nm(A, B, p, ord_hint=None, deg_cap=None):
    """Pairwise LCLM(A,B) over GF(p)(x) via Sylvester nullspace on nmod_poly.

    L = P·A = Q·B  with ord(P)=n-a, ord(Q)=n-b, n=ord_hint (default a+b).
    Solve the (a+b+1)-block system for P,Q's rational-function coeffs at each
    x-evaluation point → NO: solve directly over nmod_poly by flattening to a
    coefficient linear system.

    We flatten: P_i(x), Q_j(x) polynomials of degree ≤ deg_cap.
    Equations:  Σ_k coeff( (P·A - Q·B)_k, x^m ) = 0  for all k=0..n, m=0..(deg_cap+dA).
    """
    a, b = len(A) - 1, len(B) - 1
    n = ord_hint if ord_hint is not None else a + b
    pa, pb = n - a, n - b
    dA = max(c.degree() for c in A)
    dB = max(c.degree() for c in B)
    if deg_cap is None:
        # LCLM degree ≤ dA + dB (rough); cofactor deg ≤ that.
        deg_cap = dA + dB + n
    dP, dQ = deg_cap, deg_cap
    from math import comb
    # Precompute D^s A_j, D^s B_j
    def dks(op, kmax):
        out = [[flint.nmod_poly(list(c), p) for c in op]]
        for _ in range(kmax):
            out.append([c.derivative() for c in out[-1]])
        return out
    dAop = dks(A, pa); dBop = dks(B, pb)
    # (P·A)_k = Σ_{i,j: (i-s)+j=k, s≤i} C(i,s)·P_i·D^s A_j
    #         = Σ_{i≥?} P_i · [ Σ_{s≤i, j=k-i+s} C(i,s) D^s A_j ]
    # Let  MA[k][i] = Σ_{s=0}^{i} C(i,s) D^s A_{k-i+s}  (poly in x)
    def prodmat(dop, po_, ao_):
        MA = [[flint.nmod_poly([0], p) for _ in range(po_ + 1)] for _ in range(n + 1)]
        for k in range(n + 1):
            for i in range(po_ + 1):
                acc = flint.nmod_poly([0], p)
                for s in range(i + 1):
                    j = k - i + s
                    if 0 <= j <= ao_:
                        acc += dop[s][j] * comb(i, s)
                MA[k][i] = acc
        return MA
    MA = prodmat(dAop, pa, a)
    MB = prodmat(dBop, pb, b)
    # Unknowns: P_0..P_{pa} (each deg ≤ dP), Q_0..Q_{pb} (each deg ≤ dQ)
    ncolP = (pa + 1) * (dP + 1)
    ncolQ = (pb + 1) * (dQ + 1)
    ncol = ncolP + ncolQ
    # Equations: for each k=0..n, each x^m coeff of  Σ_i P_i·MA[k][i] − Σ_j Q_j·MB[k][j] = 0
    dRow = deg_cap + max(dA, dB)
    nrow = (n + 1) * (dRow + 1)
    ent = np.zeros((nrow, ncol), dtype=np.int64)
    for k in range(n + 1):
        for i in range(pa + 1):
            poly = MA[k][i]
            dp_ = poly.degree()
            for r in range(dp_ + 1):
                cr = int(poly[r])
                if cr == 0: continue
                # P_i(x)·x^r contributes to x^{m}=x^{r+t} with P_i coeff t
                for t in range(dP + 1):
                    m = r + t
                    if m > dRow: break
                    ent[k*(dRow+1) + m, i*(dP+1) + t] = (ent[k*(dRow+1)+m, i*(dP+1)+t] + cr) % p
        for j in range(pb + 1):
            poly = MB[k][j]
            dp_ = poly.degree()
            for r in range(dp_ + 1):
                cr = int(poly[r])
                if cr == 0: continue
                for t in range(dQ + 1):
                    m = r + t
                    if m > dRow: break
                    ent[k*(dRow+1) + m, ncolP + j*(dQ+1) + t] = \
                        (ent[k*(dRow+1)+m, ncolP + j*(dQ+1)+t] - cr) % p
    M = flint.nmod_mat(nrow, ncol, ent.flatten().tolist(), p)
    V, nul = M.nullspace()
    if nul == 0:
        return None, 0
    # Take the nullvector giving nonzero leading L_n = (P·A)_n = P_{pa}·A_a.
    # Choose the column whose P_{pa} block is nonzero.
    for kk in range(nul):
        v = np.array([int(V[i, kk]) % p for i in range(ncol)], dtype=np.int64)
        Pcoeffs = [v[i*(dP+1):(i+1)*(dP+1)] for i in range(pa + 1)]
        if np.any(Pcoeffs[pa] % p):
            break
    else:
        return None, nul
    # Assemble L = P·A, then remove polynomial content.
    Pnm = [flint.nmod_poly([int(c) for c in Pcoeffs[i]], p) for i in range(pa + 1)]
    L = op_mul_nm(Pnm, A, p)
    L = op_primitive_nm(L, p)
    return L, nul


def lclm_modp(summands_il, p, ord_hint_total=None, verbose=True):
    """Iterated pairwise LCLM.  Returns (L_nm, order, degs)."""
    ops = [op_to_nmod(A, p) for A in summands_il]
    L = ops[0]
    for k in range(1, len(ops)):
        A = ops[k]
        # Search orders ASCENDING from the minimum possible LCLM order, so
        # the first verified success is the MINIMAL-order LCLM (a descending
        # loop can accept a non-minimal common multiple at the top order once
        # caps are adequate).  At each order, escalate deg_cap — the default
        # (dA+dB+n) can undersize the true LCLM coefficient degree (apparent
        # singularities).
        base = (len(L) - 1) + (len(A) - 1)
        lo = max(len(L) - 1, len(A) - 1)
        Lnew = None
        for n_try in range(lo, base + 1):
            dL_ = max(c.degree() for c in L)
            dA_ = max(c.degree() for c in A)
            cap0 = dL_ + dA_ + n_try
            for cap in (cap0, 2*cap0, 4*cap0):
                Lnew, nul = lclm2_nm(L, A, p, ord_hint=n_try, deg_cap=cap)
                if Lnew is not None:
                    break
                if verbose:
                    print(f"  lclm stage {k}: ord={n_try} cap={cap} → nul={nul}"
                          f" (no valid nullvec)")
            if verbose and Lnew is None:
                print(f"  lclm stage {k}: ord={n_try} → nul={nul} (no valid nullvec)")
            if Lnew is not None:
                # verify right-divisibility by BOTH inputs
                okL = right_divides_nm(Lnew, L, p)
                okA = right_divides_nm(Lnew, A, p)
                if verbose:
                    dmax = max(c.degree() for c in Lnew)
                    print(f"  lclm stage {k}: ord={n_try} nul={nul} "
                          f"deg={dmax}  rdiv(prev)={okL} rdiv(new)={okA}")
                if okL and okA:
                    L = Lnew
                    break
        else:
            raise RuntimeError(f"lclm stage {k} failed")
    return L


if __name__ == '__main__':
    # Selftest of record: build the order-11 LCLM of the four arXiv:1110.1705
    # direct-summand operators (known_factors fixture) and verify each one
    # right-divides the result.  This fixture exercises both the deg_cap
    # escalation (stage 3 needs deg 247 vs the base cap) and the
    # ascending-order minimality search (a non-minimal order-12 common
    # multiple must be rejected; the true LCLM has order 11).
    import sys, os, time
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from known_factors import L3_1, L3_3, WU_2, LVU_6
    def opil(op): return [[int(c) for c in poly.all_coeffs()[::-1]] for poly in op]
    p = 1048573
    t0 = time.time()
    print("Building L11 = LCLM(L3_1, L3_3, WU_2, LVU_6) mod p...")
    L11 = lclm_modp([opil(L3_1), opil(L3_3), opil(WU_2), opil(LVU_6)], p)
    print(f"L11 order={len(L11)-1}, degs={[c.degree() for c in L11]}, "
          f"wall={time.time()-t0:.1f}s")
    assert len(L11) - 1 == 11, f"expected minimal LCLM order 11, got {len(L11)-1}"
    # Verify against each summand
    all_ok = True
    for name, A in [('L3_1',L3_1),('L3_3',L3_3),('WU_2',WU_2),('LVU_6',LVU_6)]:
        ok = right_divides_nm(L11, op_to_nmod(opil(A), p), p)
        all_ok = all_ok and ok
        print(f"  L11 = P·{name}: {ok}")
    assert all_ok, "a summand fails to right-divide the LCLM"
    print("SELFTEST PASS")
