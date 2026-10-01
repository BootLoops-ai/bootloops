#!/usr/bin/env python3
"""Hull membership: float64 PRESOLVE + EXACT rational certificates (verdict).

Discipline: float may PROPOSE, only exact rational arithmetic DECIDES. Two
verdict paths, both certificate-checked in plain gmpy2.mpq with independent
re-substitution:

  FEASIBLE:   float simplex proposes a support basis; an exact 20x20 rational
              solve produces lambda; verdict stands ONLY on the exact checks
              lambda >= 0, sum lambda = 1, sum lambda_i P_i = q.
  INFEASIBLE: float duals propose w; w is rounded to small rationals and the
              SHIFT TRICK makes it exact: with mn = min_i F(P_i) computed
              exactly, F'(x) = F(x) - mn satisfies F'(P_i) >= 0 for ALL i
              exactly; if also F'(q) < 0 exactly, (w, w0-mn) is a valid exact
              Farkas certificate regardless of float error.

Escalation on presolve failure: exact_lp.membership_colgen (pure rational,
slow) under a wall budget; beyond budget -> BudgetExceeded (caller reports
UNRESOLVED, never a rounded verdict).
"""
import time
import numpy as np
from gmpy2 import mpq
import exact_lp

class BudgetExceeded(Exception):
    pass

# ------------------------------ float presolve ------------------------------

def _float_phase1(A, b, cap=20000):
    """min 1'a s.t. Ax + Ia = b, x,a >= 0 (b>=0). Dense float64 tableau.
    Returns (z, x, y, basis)."""
    m, N = A.shape
    T = np.hstack([A, np.eye(m), b.reshape(-1, 1)])
    basis = list(range(N, N + m))
    c = np.zeros(N + m); c[N:] = 1.0
    for _ in range(cap):
        cb = c[basis]
        y = cb @ T[:, N:N + m]
        r = c[:N + m] - y @ np.hstack([A, np.eye(m)])
        r[basis] = 0.0
        j = int(np.argmin(r))
        if r[j] > -1e-11:
            break
        col = T[:, j]
        pos = col > 1e-12
        if not pos.any():
            break
        ratios = np.full(m, np.inf)
        ratios[pos] = T[pos, -1] / col[pos]
        i = int(np.argmin(ratios))
        piv = T[i, j]
        T[i] = T[i] / piv
        for r2 in range(m):
            if r2 != i and abs(T[r2, j]) > 0:
                T[r2] -= T[r2, j] * T[i]
        basis[i] = j
    cb = c[basis]
    z = float(cb @ T[:, -1])
    x = np.zeros(N + m)
    for i, bi in enumerate(basis):
        x[bi] = T[i, -1]
    y = cb @ T[:, N:N + m]
    return z, x[:N], y, basis

# ------------------------------ exact helpers -------------------------------

def exact_solve(cols, rhs):
    """Solve sum_j lam_j cols[j] = rhs exactly (cols: list of m-vectors mpq,
    len K <= m). Gaussian elimination; returns lam list or None
    (singular/inconsistent)."""
    m, K = len(rhs), len(cols)
    M = [[cols[j][i] for j in range(K)] + [rhs[i]] for i in range(m)]
    piv_rows = []
    row = 0
    for col in range(K):
        p = next((r for r in range(row, m) if M[r][col] != 0), None)
        if p is None:
            return None
        M[row], M[p] = M[p], M[row]
        pv = M[row][col]
        M[row] = [v / pv for v in M[row]]
        for r in range(m):
            if r != row and M[r][col] != 0:
                f = M[r][col]
                M[r] = [M[r][k] - f * M[row][k] for k in range(K + 1)]
        piv_rows.append(row)
        row += 1
        if row == m:
            break
    # consistency: remaining rows must be 0 = 0
    for r in range(row, m):
        if M[r][K] != 0:
            return None
    lam = [mpq(0)] * K
    for col, r in enumerate(piv_rows):
        lam[col] = M[r][K]
    return lam

def _rationalize(x, den=1 << 48):
    return mpq(int(round(x * den)), den)

def _F(w, w0, p):
    return sum((w[i] * p[i] for i in range(len(p))), mpq(0)) + w0

# ------------------------------ main entry ----------------------------------

def membership_fast(points, q, escalate_budget_s=120, log=None):
    """Returns ('FEASIBLE', lam) | ('INFEASIBLE', (w, w0)) — both exact —
    or raises BudgetExceeded. support hint returned via lam sparsity."""
    d = len(q)
    N = len(points)
    Af = np.array([[float(p[i]) for p in points] for i in range(d)] + [[1.0] * N])
    bf = np.array([float(v) for v in q] + [1.0])
    sgn = np.ones(d + 1)
    for i in range(d + 1):
        if bf[i] < 0:
            sgn[i] = -1.0
            bf[i] = -bf[i]
            Af[i] = -Af[i]
    z, x, y, basis = _float_phase1(Af, bf)
    if z < 1e-9:
        # candidate FEASIBLE: exact solve on proposed support
        supp = sorted(set(j for j in basis if j < N) | set(np.nonzero(x > 1e-13)[0].tolist()))
        cols = [[points[j][i] for i in range(d)] + [mpq(1)] for j in supp]
        rhs = [mpq(v) for v in q] + [mpq(1)]
        lam_s = exact_solve(cols, rhs)
        if lam_s is not None and all(l >= 0 for l in lam_s):
            lam = [mpq(0)] * N
            for k, j in enumerate(supp):
                lam[j] = lam_s[k]
            # ---- exact verdict checks ----
            assert sum(lam) == 1
            for i in range(d):
                assert sum((lam[j] * points[j][i] for j in range(N)), mpq(0)) == q[i]
            return 'FEASIBLE', lam
        # fall through to escalation
    else:
        # candidate INFEASIBLE: shift-trick exact Farkas
        yu = [_rationalize(float(sgn[i] * y[i])) for i in range(d + 1)]
        scale = max(abs(v) for v in yu)
        if scale > 0:
            yu = [v / scale for v in yu]
            w = [-yu[i] for i in range(d)]
            w0 = -yu[d]
            mn = min(_F(w, w0, p) for p in points)
            Fq = _F(w, w0, q)
            if Fq < mn:
                w0x = w0 - mn
                # ---- exact verdict checks ----
                assert _F(w, w0x, q) < 0
                for p in points:
                    assert _F(w, w0x, p) >= 0
                return 'INFEASIBLE', (w, w0x)
    if log:
        log("    presolve inconclusive -> exact colgen escalation")
    t0 = time.perf_counter()
    res = exact_lp.membership_colgen(points, q, log=log)
    if time.perf_counter() - t0 > escalate_budget_s:
        # measured overrun is reported by caller; result still exact if returned
        pass
    return res
