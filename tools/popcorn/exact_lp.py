#!/usr/bin/env python3
"""Exact rational LP feasibility for hull membership (phase-1 simplex, Bland).

membership(points, q):
  decide exists lambda >= 0, sum lambda = 1, sum lambda_i P_i = q, ALL in gmpy2.mpq.
  Returns ('FEASIBLE', lam) with lam independently re-verified, or
          ('INFEASIBLE', (w, w0)) with the affine Farkas functional
          F(x) = <w, x> + w0 satisfying F(P_i) >= 0 for all i, F(q) < 0,
          independently re-verified in plain rational arithmetic.
No floats anywhere; scipy deliberately not imported — this module is the
slow, fully rational validation oracle for the float-assisted siblings.
"""
from gmpy2 import mpq

def _phase1(cols, b):
    """min 1'a s.t. Sum_j x_j cols[j] + I a = b, x,a >= 0 (b >= 0 required).
    Full-tableau, Bland's rule. Returns (z, x_basic dict, y_dual list)."""
    m = len(b)
    N = len(cols)
    # tableau rows: [cols | I | b]
    T = [[cols[j][i] for j in range(N)] + [mpq(1) if k == i else mpq(0) for k in range(m)]
         + [b[i]] for i in range(m)]
    basis = [N + i for i in range(m)]
    ZERO = mpq(0)
    cost = lambda j: mpq(1) if j >= N else ZERO
    it = 0
    while True:
        it += 1
        # duals y_i = sum_r cB_r * Binv[r][i] = cost-weighted artificial columns
        y = [sum((T[r][N + i] for r in range(m) if basis[r] >= N), ZERO) for i in range(m)]
        # reduced costs r_j = c_j - y'A_j ; Bland: first j with r_j < 0
        enter = -1
        for j in range(N + m):
            if j in basis:
                continue
            yA = sum((y[i] * (cols[j][i] if j < N else (mpq(1) if j - N == i else ZERO))
                      for i in range(m)), ZERO)
            if cost(j) - yA < 0:
                enter = j
                break
        if enter < 0:
            break
        # ratio test (Bland tie-break: smallest basis var)
        leave, best = -1, None
        for r in range(m):
            a = T[r][enter]
            if a > 0:
                ratio = T[r][-1] / a
                if best is None or ratio < best or (ratio == best and basis[r] < basis[leave]):
                    best, leave = ratio, r
        if leave < 0:
            raise RuntimeError("phase-1 unbounded (impossible)")
        piv = T[leave][enter]
        T[leave] = [v / piv for v in T[leave]]
        for r in range(m):
            if r != leave and T[r][enter] != 0:
                f = T[r][enter]
                T[r] = [T[r][c] - f * T[leave][c] for c in range(N + m + 1)]
        basis[leave] = enter
        if it > 20000:
            raise RuntimeError("iteration cap")
    z = sum((T[r][-1] for r in range(m) if basis[r] >= N), ZERO)
    x = {basis[r]: T[r][-1] for r in range(m) if basis[r] < N}
    y = [sum((T[r][N + i] for r in range(m) if basis[r] >= N), ZERO) for i in range(m)]
    return z, x, y

def membership(points, q):
    d = len(q)
    N = len(points)
    for p in points:
        assert len(p) == d
    # rows: d coordinate equations + 1 sum row
    cols = []
    for p in points:
        cols.append([mpq(v) for v in p] + [mpq(1)])
    b = [mpq(v) for v in q] + [mpq(1)]
    # flip rows to b >= 0
    sgn = [mpq(1)] * (d + 1)
    for i in range(d + 1):
        if b[i] < 0:
            sgn[i] = mpq(-1)
            b[i] = -b[i]
            for j in range(N):
                cols[j][i] = -cols[j][i]
    z, x, y = _phase1(cols, b)
    if z == 0:
        lam = [x.get(j, mpq(0)) for j in range(N)]
        # ---- independent verification ----
        assert all(l >= 0 for l in lam)
        assert sum(lam) == 1
        for i in range(d):
            s = sum((lam[j] * points[j][i] for j in range(N)), mpq(0))
            assert s == q[i], f"lambda re-verify failed at coord {i}"
        return 'FEASIBLE', lam
    # infeasible: y'A_j <= 0 for all j, y'b > 0 (min problem optimality)
    # unflip: yu_i = sgn_i * y_i acts on original rows
    yu = [sgn[i] * y[i] for i in range(d + 1)]
    w = [-yu[i] for i in range(d)]
    w0 = -yu[d]
    # ---- independent verification ----
    Fq = sum((w[i] * q[i] for i in range(d)), mpq(0)) + w0
    assert Fq < 0, f"Farkas re-verify: F(q) = {Fq} not < 0"
    for j in range(N):
        Fp = sum((w[i] * points[j][i] for i in range(d)), mpq(0)) + w0
        assert Fp >= 0, f"Farkas re-verify failed at point {j}"
    return 'INFEASIBLE', (w, w0)

def _F(w, w0, p):
    return sum((w[i] * p[i] for i in range(len(p))), mpq(0)) + w0

def membership_colgen(points, q, seed=None, batch=6, log=None):
    """Exact column generation: solve small exact LPs on an active subset,
    verify every certificate against the FULL point set in exact arithmetic.
    FEASIBLE lambda is globally valid (subset of columns); INFEASIBLE Farkas
    is accepted only when it separates ALL points. Returns same API as
    membership(), with lam indexed over the full point list."""
    N = len(points)
    d = len(q)
    if seed is None:
        step = max(1, N // 24)
        seed = list(range(0, N, step))
    active = sorted(set(seed))
    rounds = 0
    while True:
        rounds += 1
        sub = [points[i] for i in active]
        st, res = membership(sub, q)
        if st == 'FEASIBLE':
            lam_full = [mpq(0)] * N
            for loc, glob in enumerate(active):
                lam_full[glob] = res[loc]
            return 'FEASIBLE', lam_full
        w, w0 = res
        viol = []
        for j in range(N):
            if j in active:
                continue
            v = _F(w, w0, points[j])
            if v < 0:
                viol.append((v, j))
        if not viol:
            # Farkas already exactly verified on active inside membership();
            # loop above verified the rest -> certificate valid on all points
            return 'INFEASIBLE', (w, w0)
        viol.sort()
        newcols = [j for _, j in viol[:batch]]
        active = sorted(set(active) | set(newcols))
        if log:
            log(f"    colgen round {rounds}: +{len(newcols)} cols (|active|={len(active)})")
        if rounds > 200:
            raise RuntimeError("colgen cap")
