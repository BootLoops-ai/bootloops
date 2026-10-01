#!/usr/bin/env python3
"""emitter.py — lambda-witness emission (the retrofit driver).

Solves  lam^T Abar = e_t  by scalar Wiedemann over F_p (matrix-free,
min-poly reuse across targets), reads off the certified coefficients
c = -(lam^T A)|_masters, and writes v1 witnesses. If a claimed table row is
supplied, the emitted (certified) c is compared against it: a mismatch means
the claimed row is NOT implied by the system — witness emission is itself a
wrong-table detector (the wrong-table exhibit class).

Emission is a SOLVER (this file may know about pivots/partitions); only the
verify path (core.py) carries the solver-agnostic invariant.

Measured record: oracle match 10/10 + 10/10; retrofit 84/84 at 2 primes;
matvec max 3.73n < 5n; min-poly reuse floor 0.86-1.0n; 2.2-2.3 s/(row,prime)
on the 7,135-equation reference family.
"""
import time

import numpy as np


class FpSystem:
    """CSR-ish arrays + mod-p matvecs (int64 numpy) over a pivot partition.

    rows_terms: list of list[(col, val)] (evaluated at one (point, prime));
    pivots: the square-elimination pivot columns (|pivots| must equal the
    number of independent rows for the direct Abar^T path).
    """

    def __init__(self, rows_terms, pivots, all_cols, p):
        self.p = p
        self.n_eqs = len(rows_terms)
        self.pivots = list(pivots)
        self.pidx = {w: i for i, w in enumerate(self.pivots)}
        self.others = [w for w in all_cols if w not in self.pidx]
        self.oidx = {w: i for i, w in enumerate(self.others)}
        r_in, c_in, v_in = [], [], []
        r_out, c_out, v_out = [], [], []
        for i, row in enumerate(rows_terms):
            for w, v in row:
                j = self.pidx.get(w)
                if j is not None:
                    r_in.append(i); c_in.append(j); v_in.append(v)   # noqa: E702
                else:
                    r_out.append(i); c_out.append(self.oidx[w]); v_out.append(v)  # noqa: E702,E501
        self.ri = np.asarray(r_in, np.int64)
        self.ci = np.asarray(c_in, np.int64)
        self.vi = np.asarray(v_in, np.int64) % p
        self.ro = np.asarray(r_out, np.int64)
        self.co = np.asarray(c_out, np.int64)
        self.vo = np.asarray(v_out, np.int64) % p
        self.n = len(self.pivots)

    def mv_A(self, x):      # Abar @ x   (len n -> len n_eqs)
        t = (self.vi * x[self.ci]) % self.p
        return np.bincount(self.ri, t, minlength=self.n_eqs).astype(
            np.int64) % self.p

    def mv_AT(self, y):     # Abar^T @ y (len n_eqs -> len n)
        t = (self.vi * y[self.ri]) % self.p
        return np.bincount(self.ci, t, minlength=self.n).astype(
            np.int64) % self.p

    def mv_out_T(self, y):  # (non-pivot cols)^T @ y -> len(others)
        t = (self.vo * y[self.ro]) % self.p
        return np.bincount(self.co, t, minlength=len(self.others)).astype(
            np.int64) % self.p


def bm(seq, p):
    """Berlekamp-Massey over F_p (probe-harness port)."""
    s = np.asarray(seq, np.int64)
    N = len(s)
    C = np.zeros(N + 1, np.int64); C[0] = 1   # noqa: E702
    B = np.zeros(N + 1, np.int64); B[0] = 1   # noqa: E702
    L, m, b = 0, 1, 1
    lenC = lenB = 1
    for i in range(N):
        if L:
            d = int((C[:L + 1] * s[i - L:i + 1][::-1] % p).sum() % p)
        else:
            d = int(s[i] % p)
        if d == 0:
            m += 1
            continue
        coef = d * pow(b, p - 2, p) % p
        if 2 * L <= i:
            T, lenT = C.copy(), lenC
            C[m:m + lenB] = (C[m:m + lenB] - coef * B[:lenB]) % p
            lenC = max(lenC, m + lenB)
            L, B, lenB, b, m = i + 1 - L, T, lenT, d, 1
        else:
            C[m:m + lenB] = (C[m:m + lenB] - coef * B[:lenB]) % p
            lenC = max(lenC, m + lenB)
            m += 1
    return C[:L + 1] % p, L


def _reversal_deflate(C, L, p):
    rev = np.ascontiguousarray(C[:L + 1][::-1])
    k = 0
    while k <= L and int(rev[k]) % p == 0:
        k += 1
    if k > L:
        return None
    return rev[k:]


def wiedemann_solve(mv, n, b_vec, p, rng, minpoly=None, max_mult=6):
    """Solve M z = b given matvec mv; every candidate is EXACTLY verified,
    so min-poly reuse and early BM stops are sound. (Probe-harness port.)"""
    stats = {"matvecs": 0}

    def recon(mp):
        L = len(mp) - 1
        acc = (int(mp[L]) * b_vec) % p
        for i in range(L - 1, 0, -1):
            acc = (mv(acc) + int(mp[i]) * b_vec) % p
            stats["matvecs"] += 1
        return (-pow(int(mp[0]), p - 2, p) * acc) % p

    if minpoly is not None and int(minpoly[0]) % p:
        z = recon(minpoly)
        r = (mv(z) - b_vec) % p
        stats["matvecs"] += 1
        if not r.any():
            stats["mode"] = "reused_minpoly"
            return z, stats
        stats["reuse_failed"] = True
    u = rng.integers(1, p, n, dtype=np.int64)
    v = b_vec.copy()
    seq = []
    cap = 2 * n + 16
    next_bm = min(1024, cap)
    last_L = -1
    while True:
        while len(seq) < next_bm:
            seq.append(int((u * v % p).sum() % p))
            v = mv(v)
            stats["matvecs"] += 1
        C, L = bm(seq, p)
        if (L == last_L and 2 * L + 32 <= len(seq)) or len(seq) >= cap:
            break
        last_L = L
        next_bm = min(max(2 * L + 64, int(1.5 * len(seq))), cap)
        if len(seq) >= max_mult * n:
            break
    stats["minpoly_deg"] = int(L)
    mp = _reversal_deflate(C, L, p)
    if mp is None:
        stats["m0_zero"] = True
        return None, stats
    z = recon(mp)
    r = (mv(z) - b_vec) % p
    stats["matvecs"] += 1
    if r.any():
        stats["residual_nonzero"] = True
        return None, stats
    stats["mode"] = "fresh_minpoly"
    stats["minpoly"] = mp
    return z, stats


def solve_target(S, t_pidx, rng, minpoly=None):
    """lam with lam^T Abar = e_t. Direct B = Abar^T; retry: Kaltofen-Saunders
    G = Abar^T D Abar. Returns (lam or None, stats, minpoly_for_reuse)."""
    p, n = S.p, S.n
    b_vec = np.zeros(n, np.int64); b_vec[t_pidx] = 1   # noqa: E702
    z, st = wiedemann_solve(S.mv_AT, n, b_vec, p, rng, minpoly=minpoly)
    if z is not None:
        st["path"] = "direct_AT"
        return z, st, st.get("minpoly", minpoly)
    st1 = st
    D = rng.integers(1, p, S.n_eqs, dtype=np.int64)
    mvG = lambda x: S.mv_AT((D * S.mv_A(x)) % p)       # noqa: E731
    z, st = wiedemann_solve(mvG, n, b_vec, p, rng)
    if z is None:
        st["path"] = "FAILED_after_KS_retry"
        st["matvecs"] = st1["matvecs"] + 2 * st["matvecs"]
        return None, st, None
    lam = (D * S.mv_A(z)) % p
    st["path"] = "KS_retry"
    st["matvecs"] = st1["matvecs"] + 2 * st["matvecs"] + 1
    return lam, st, None


def emit_for_target(S, target_col, master_cols, rng, minpoly=None,
                    claimed=None):
    """Emit witness data for one target column.

    master_cols: set of columns allowed to carry coefficients (declared
    master basis). Any other nonzero residual column (shell class) voids
    the certificate.
    claimed: optional dict col->val; compared against certified c.
    Returns (record dict, lam or None, new minpoly).
    """
    p = S.p
    t0 = time.time()
    lam, st, mp_new = solve_target(S, S.pidx[target_col], rng,
                                   minpoly=minpoly)
    rec = {"target_col": int(target_col), "p": p,
           "solve_wall_s": round(time.time() - t0, 3),
           "matvecs": st.get("matvecs"), "mode": st.get("mode"),
           "path": st.get("path")}
    if lam is None:
        rec["status"] = "SOLVE_FAILED"
        return rec, None, mp_new
    r_out = S.mv_out_T(lam)
    c, shell_hit = {}, 0
    for j, w in enumerate(S.others):
        if r_out[j]:
            if w in master_cols:
                c[w] = int((-int(r_out[j])) % p)
            else:
                shell_hit += 1
    if shell_hit:
        rec["status"] = "SHELL_RESIDUAL"
        rec["shell_hit_cols"] = shell_hit
        return rec, None, mp_new
    rec["c"] = c
    rec["lambda_nnz"] = int((lam != 0).sum())
    if claimed is not None:
        claimed_n = {int(k): int(v) % p for k, v in claimed.items()
                     if int(v) % p}
        if claimed_n == c:
            rec["status"] = "CERTIFIED"
        else:
            rec["status"] = "CLAIM_MISMATCH"   # wrong-table caught at emit
            rec["n_differing_entries"] = len(
                [m for m in set(claimed_n) | set(c)
                 if claimed_n.get(m) != c.get(m)])
    else:
        rec["status"] = "CERTIFIED"
    return rec, lam, mp_new
