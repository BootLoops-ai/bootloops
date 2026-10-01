#!/usr/bin/env python3
"""Exact hull membership via integer cone form + exact dual simplex (flint).

BASIS COMPLETION (the load-bearing design choice): _basis_complete does NOT
do greedy per-column exact rank checks (measured ~51 s PER rank check at
n=50 -> tens of minutes just to seed a basis). It extracts the pivot columns
of ONE exact fraction-free RREF over the candidate columns, ordered warm
columns first and remaining columns by ascending bit size (measured at n=50:
one bit-sorted rref = 63 s vs 264 s unsorted vs ~40 min greedy). Bit-ordering
also keeps the working basis small-bit, which speeds every subsequent exact
B.solve.

q in conv{P_i}  <=>  exists mu >= 0 : G mu = h, where
  g_i = D_i [xi(t_i); s(t_i)] in Z^m (positive column scaling, D_i clears
        denominators), h = D_q [q; 1] in Z^m,  m = n.
Dual simplex with cost 0: any basis is dual-feasible; iterate
  leaving  = Bland smallest basic index with negative basic value,
  entering = Bland smallest nonbasic j with (row of B^{-1}) . g_j < 0;
no entering column -> the row IS a Farkas certificate. Every terminal answer
is re-verified independently in gmpy2.mpq. float64 enters ONLY through
the optional warm-start basis hint; every decision comparison is exact.
"""
import flint
from gmpy2 import mpq, lcm

class LPStall(Exception):
    pass

def col_from_point(P):
    """P: exact simplex point (mpq tuple). Integer column P*D, plus D.
    NOTE: the sum row is intentionally dropped — every point lies on the
    simplex, so the normalization row is linearly dependent; convexity is
    implied (sum both sides coordinatewise: sum mu'_i = 1)."""
    D = 1
    for x in P:
        D = lcm(D, x.denominator)
    g = [int(x * D) for x in P]
    return g, D

def _dot(wq, g):
    return sum(wq[i] * g[i] for i in range(len(g)))

def _rank_of(cols, idxs, m):
    M = flint.fmpz_mat(m, len(idxs))
    for r in range(m):
        for k, j in enumerate(idxs):
            M[r, k] = cols[j][r]
    return M.rank()

def _basis_complete(cols, order, m):
    """Nonsingular basis from candidate order (indices into cols).
    Phase-0b change: warm head fast path kept; otherwise ONE exact RREF
    (warm-first, then bit-ascending columns) replaces the greedy per-column
    rank loop of phase 0."""
    head = order[:m]
    if len(head) == m and _rank_of(cols, head, m) == m:
        return list(head)
    # completion from scratch: pure bit-ascending order (fat deep-ladder
    # columns cost the most in every later exact solve), rref on a growing
    # prefix so the fat columns only enter the elimination if the cheap ones
    # cannot carry rank m.
    allc = sorted(set(order),
                  key=lambda j: max(abs(v).bit_length() for v in cols[j]))
    K = min(len(allc), 2 * m)
    while True:
        full = allc[:K]
        M = flint.fmpz_mat(m, len(full))
        for r in range(m):
            for k, j in enumerate(full):
                M[r, k] = cols[j][r]
        R, _den, rank = M.rref()
        if rank == m:
            break
        if K >= len(allc):
            raise LPStall(f"cannot complete basis (rank {rank} < {m})")
        K = min(len(allc), 2 * K)
    chosen = []
    k0 = 0
    for r in range(rank):
        k = k0
        while k < len(full) and R[r, k] == 0:
            k += 1
        chosen.append(full[k])
        k0 = k + 1
    if len(chosen) < m:
        raise LPStall(f"cannot complete basis (pivots {len(chosen)} < {m})")
    return chosen

def dual_simplex(cols, h, warm=None, cap=1500, log=None):
    """cols: list of integer columns (len-m int lists); h: int list.
    Returns ('FEASIBLE', {j: mu_j (mpq)}) or ('INFEASIBLE', w (list mpq, len m))
    with w.g_j >= 0 for all j and w.h < 0 — caller re-verifies."""
    m = len(h)
    N = len(cols)
    order = list(warm or []) + list(range(0, N, max(1, N // (2 * m)))) + list(range(N))
    seen = set(); order = [j for j in order if not (j in seen or seen.add(j))]
    basis = _basis_complete(cols, order, m)
    hmat = flint.fmpz_mat(m, 1)
    for r in range(m):
        hmat[r, 0] = h[r]
    visited = set()
    bland = False
    for it in range(cap):
        key = frozenset(basis)
        if key in visited:
            bland = True  # cycling guard: greedy revisited a basis -> Bland
        visited.add(key)
        B = flint.fmpz_mat(m, m)
        for r in range(m):
            for k in range(m):
                B[r, k] = cols[basis[k]][r]
        x = B.solve(hmat)  # fmpq_mat
        xv = [mpq(int(x[r, 0].p), int(x[r, 0].q)) for r in range(m)]
        neg = [k for k in range(m) if xv[k] < 0]
        if not neg:
            return 'FEASIBLE', {basis[k]: xv[k] for k in range(m)}
        if bland:
            k_leave = min(neg, key=lambda k: basis[k])
        else:
            k_leave = min(neg, key=lambda k: xv[k])  # most negative
        e = flint.fmpz_mat(m, 1)
        for r in range(m):
            e[r, 0] = 1 if r == k_leave else 0
        z = B.transpose().solve(e)
        wq = [mpq(int(z[r, 0].p), int(z[r, 0].q)) for r in range(m)]
        # w.g_basis[k] = delta_{k,k_leave}; w.h = xv[k_leave] < 0
        enter, best = None, None
        for j in range(N):
            if j in basis:
                continue
            dj = _dot(wq, cols[j])
            if dj < 0:
                if bland:
                    enter = j
                    break
                sc = dj / max(abs(g) for g in cols[j])  # scale-invariant violation
                if best is None or sc < best:
                    best, enter = sc, j
        if enter is None:
            return 'INFEASIBLE', wq
        basis[k_leave] = enter
        if log and it % 50 == 49:
            log(f"      dual-simplex iter {it+1}{' (bland)' if bland else ''}")
    raise LPStall("dual simplex cap hit")

def membership_exact(points, q, warm=None, log=None, _cols=None, _scales=None):
    """points: list of exact simplex tuples (mpq); q: exact simplex point.
    Fully exact; returns ('FEASIBLE', lam list) | ('INFEASIBLE', (w, w0)) and
    the winning basis via third element for warm starts."""
    if _cols is None:
        _cols, _scales = [], []
        for P in points:
            g, D = col_from_point(P)
            _cols.append(g)
            _scales.append(D)
    cols, scales = _cols, _scales
    hg, Dq = col_from_point(q)
    st, res = dual_simplex(cols, hg, warm=warm, log=log)
    d = len(q)
    if st == 'FEASIBLE':
        lam = [mpq(0)] * len(points)
        for j, mu in res.items():
            lam[j] = mu * scales[j] / Dq
        # ---- independent exact verification ----
        assert all(l >= 0 for l in lam)
        assert sum(lam) == 1
        for i in range(d):
            assert sum((lam[j] * points[j][i] for j in range(len(points))), mpq(0)) == q[i]
        return 'FEASIBLE', lam, sorted(res.keys())
    w_cone = res
    w = [w_cone[i] for i in range(d)]
    w0 = mpq(0)
    # ---- independent exact verification ----
    Fq = sum((w[i] * q[i] for i in range(d)), mpq(0)) + w0
    assert Fq < 0, "Farkas re-verify failed at q"
    for j, P in enumerate(points):
        assert sum((w[i] * P[i] for i in range(d)), mpq(0)) + w0 >= 0, \
            f"Farkas re-verify failed at point {j}"
    return 'INFEASIBLE', (w, w0), None
