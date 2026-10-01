#!/usr/bin/env python3
"""rf_lib.py — order-4 right-factor extraction machinery for theta-form
operators.

Operators in THETA form: rows[k] = ascending z-coeff list of P_k(z),
L = sum_k P_k(z) theta^k.  Semantics note (verified on a quintic control):
the annihilator.py convention is (rec shift order = ODE z-degree,
n-degree = ODE theta-order) — beware tables that swap the two slots.

Route B (MUM-tower): formal log-solutions of L at z=0 mod p capped at
log^3; y = sum_m f_m(z) log^m z, coeff 4-vectors a_n, theta acts as
(nI + Lam), (Lam a)_m = (m+1) a_{m+1}.  Seed pure-log^3 needs indicial
root 0 mult >= 4.  Fit order-4 theta-form R (deg scan) killing the top
tower solution; verify by Ore right-division mod p (ore_rdiv_modp,
validated on known-good factors first — its law); multi-prime CRT +
rational reconstruction -> exact R; exact gates downstream.

Imports resolve through the annihilator package (the lifted ore_rdiv_modp
primitive); no paths outside the package are required.
"""
import sys, os
from fractions import Fraction as Fr
from math import gcd, comb

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))))  # tools/ — so `annihilator` resolves
import numpy as np
from annihilator import find_recurrence_fast, nullspace_mod_fast, P31
from annihilator.ore.ore_rdiv_modp import ore_rightdiv_modp

PRIMES = [2147483647, 2147483629, 2147483587, 2147483579]  # all < 2^31


def rec_to_theta_modp(r, d, v, p):
    """(r,d,v) from find_recurrence_fast -> theta rows mod p.
    ODE: sum_j z^{r-j} c_j(theta - j); c_j(n) = sum_i v[j(d+1)+i] n^i."""
    rows = [[0] * (r + 1) for _ in range(d + 1)]
    for j in range(r + 1):
        cj = [int(v[j * (d + 1) + i]) % p for i in range(d + 1)]
        # c_j(n - j): expand via binomials
        sh = [0] * (d + 1)
        for i, ci in enumerate(cj):
            if ci:
                for k in range(i + 1):
                    sh[k] = (sh[k] + ci * comb(i, k) * pow(-j % p, i - k, p)) % p
        for k in range(d + 1):
            rows[k][r - j] = sh[k]
    return [row for row in rows]


def theta_to_D_modp(rows, p):
    """theta rows -> D-form op_il (list over D^j of ascending z-polys) mod p.
    theta^k = sum_j S2(k,j) z^j D^j."""
    K = len(rows) - 1
    S2 = [[0] * (K + 1) for _ in range(K + 1)]
    S2[0][0] = 1
    for n in range(1, K + 1):
        for j in range(1, n + 1):
            S2[n][j] = (S2[n - 1][j - 1] + j * S2[n - 1][j]) % p
    out = []
    for j in range(K + 1):
        # coeff of D^j: z^j * sum_k S2[k][j] P_k(z)
        acc = [0] * (max(len(r) for r in rows) + j)
        for k in range(j, K + 1):
            if S2[k][j]:
                for m, c in enumerate(rows[k]):
                    acc[m + j] = (acc[m + j] + S2[k][j] * (c % p)) % p
        out.append(acc)
    while len(out) > 1 and not any(out[-1]):
        out.pop()
    return out


def theta_mul(A, B):
    """Ore product A*B in theta form over Fraction: theta z^b = z^b(theta+b).
    A,B = theta rows (lists of ascending z-coeff lists)."""
    Ka, Kb = len(A) - 1, len(B) - 1
    out = [[Fr(0)] for _ in range(Ka + Kb + 1)]
    for ka in range(Ka + 1):
        for ja, ca in enumerate(A[ka]):
            if ca == 0:
                continue
            for kb in range(Kb + 1):
                for jb, cb in enumerate(B[kb]):
                    if cb == 0:
                        continue
                    # z^ja th^ka z^jb th^kb = z^{ja+jb} (th+jb)^ka th^kb
                    for t in range(ka + 1):
                        co = Fr(ca) * cb * comb(ka, t) * jb ** (ka - t)
                        k, j = t + kb, ja + jb
                        while len(out[k]) <= j:
                            out[k].append(Fr(0))
                        out[k][j] += co
    return out


# ---------------- mod-p MUM tower (log-degree <= 3) ------------------------
def _Qj_mats(rows, p):
    """rows mod p -> list over z-power j of 4x4 matrix-poly evaluator input:
    Qj(theta) coeffs qj[k]."""
    smax = max(len(r) for r in rows) - 1
    return [[rows[k][j] % p if j < len(rows[k]) else 0
             for k in range(len(rows))] for j in range(smax + 1)]


def _mat_theta(n, p):
    """4x4 matrix of theta on log-vectors at z^n: nI + Lam, (Lam a)_m=(m+1)a_{m+1}."""
    M = [[0] * 4 for _ in range(4)]
    for m in range(4):
        M[m][m] = n % p
        if m < 3:
            M[m][m + 1] = m + 1
    return M


def _mmul(A, B, p):
    return [[sum(A[i][k] * B[k][j] for k in range(4)) % p for j in range(4)]
            for i in range(4)]


def _mpolyval(q, n, p):
    """Q(nI+Lam) for Q coeffs q[k], as 4x4 mod p (Horner)."""
    X = _mat_theta(n, p)
    R = [[0] * 4 for _ in range(4)]
    for c in reversed(q):
        R = _mmul(R, X, p)
        for i in range(4):
            R[i][i] = (R[i][i] + c) % p
    return R


def _solve4(M, b, p):
    """solve 4x4 mod p; returns (sol, freedom_dim) or None if inconsistent.
    Free vars set to 0."""
    A = [[M[i][j] % p for j in range(4)] + [b[i] % p] for i in range(4)]
    piv, prow = [], 0
    for col in range(4):
        r = next((i for i in range(prow, 4) if A[i][col]), None)
        if r is None:
            continue
        A[prow], A[r] = A[r], A[prow]
        inv = pow(A[prow][col], p - 2, p)
        A[prow] = [x * inv % p for x in A[prow]]
        for i in range(4):
            if i != prow and A[i][col]:
                f = A[i][col]
                A[i] = [(A[i][j] - f * A[prow][j]) % p for j in range(5)]
        piv.append(col)
        prow += 1
    for i in range(prow, 4):
        if A[i][4]:
            return None
    x = [0] * 4
    for i, col in enumerate(piv):
        x[col] = A[i][4]
    return x, 4 - len(piv)


def indicial_modp(rows, p):
    """indicial poly at 0 = Q_0(theta) coeffs mod p (ascending)."""
    return [rows[k][0] % p if len(rows[k]) > 0 else 0 for k in range(len(rows))]


def _poly_eval(q, n, p):
    v = 0
    for c in reversed(q):
        v = (v * n + c) % p
    return v


def mum_tower_modp(rows, p, N):
    """Top log^3 formal solution of L at 0 mod p, N terms.
    Returns (avec list of 4-vectors, resonances list, None) or (None, res, why)."""
    Qj = _Qj_mats(rows, p)
    smax = len(Qj) - 1
    ind = indicial_modp(rows, p)
    # mult of root 0
    m0 = 0
    while m0 < len(ind) and ind[m0] % p == 0:
        m0 += 1
    if m0 < 4:
        return None, [], f"indicial 0-mult {m0} < 4 (no log^3 MUM tower)"
    a = [[0, 0, 0, 1]]                      # seed e_3 (pure log^3)
    M00 = _mpolyval(Qj[0], 0, p)
    if any((sum(M00[i][j] * a[0][j] for j in range(4))) % p for i in range(4)):
        return None, [], "seed e_3 not killed at n=0"
    resonances = []
    for n in range(1, N):
        rhs = [0, 0, 0, 0]
        for j in range(1, min(n, smax) + 1):
            Mj = _mpolyval(Qj[j], n - j, p)
            for i in range(4):
                rhs[i] = (rhs[i] + sum(Mj[i][t] * a[n - j][t]
                                       for t in range(4))) % p
        rhs = [(-x) % p for x in rhs]
        if _poly_eval(ind, n, p) == 0:
            sol = _solve4(_mpolyval(Qj[0], n, p), rhs, p)
            if sol is None:
                return None, resonances, f"inconsistent resonance at n={n}"
            a.append(sol[0])
            resonances.append((n, sol[1]))
        else:
            # back-substitute upper-triangular Q0(nI+Lam)
            M = _mpolyval(Qj[0], n, p)
            x = [0] * 4
            for i in range(3, -1, -1):
                s = sum(M[i][j] * x[j] for j in range(i + 1, 4)) % p
                x[i] = (rhs[i] - s) * pow(M[i][i], p - 2, p) % p
            a.append(x)
    return a, resonances, None


def fit_order4_modp(avec, p, dmax, order=4, verify_rows=30, dlist=None):
    """Fit R = sum_{i<=order, j<=d} r_ij z^j th^i killing the tower top.
    Row (n, m-component): sum_j sum_i r_ij [ (n-j I+Lam)^i a_{n-j} ]_m = 0.
    Scan d ascending; return (d, theta_rows_modp, nulldim) or None."""
    N = len(avec)
    # precompute powers (nI+Lam)^i a_n for i<=order
    pw = [[avec[n]] for n in range(N)]
    for n in range(N):
        X = _mat_theta(n, p)
        for i in range(order):
            v = pw[n][-1]
            pw[n].append([sum(X[r][c] * v[c] for c in range(4)) % p
                          for r in range(4)])
    for d in (dlist if dlist is not None else range(dmax + 1)):
        ncols = (order + 1) * (d + 1)
        if 4 * N < ncols + verify_rows:
            continue
        A = np.zeros((4 * N, ncols), dtype=np.int64)
        for n in range(N):
            for j in range(min(d, n) + 1):
                for i in range(order + 1):
                    col = i * (d + 1) + j
                    for m in range(4):
                        A[4 * n + m, col] = pw[n - j][i][m]
        v = nullspace_mod_fast(A % p, p)
        if v is None:
            continue
        res = (A @ np.array([int(x) for x in v], dtype=object)) % p
        if any(int(x) % p for x in res):
            continue
        rows = [[int(v[i * (d + 1) + j]) % p for j in range(d + 1)]
                for i in range(order + 1)]
        return d, rows, 1
    return None


# ------------- multi-prime CRT + rational reconstruction -------------------
def _ratrec(c, m):
    """Wang rational reconstruction of c mod m; |num|,|den| <= sqrt(m/2)."""
    if c % m == 0:
        return Fr(0)
    a, b, x, y = m, c % m, 0, 1
    bnd = int((m // 2) ** 0.5)
    while b > bnd:
        a, (q, b) = b, divmod(a, b)
        x, y = y, x - q * y
    if b == 0 or abs(y) > bnd or gcd(b, abs(y)) != 1:
        return None
    return Fr(b, y) if y > 0 else Fr(-b, -y)


def rows_ratrec(rows_by_p, primes):
    """rows_by_p: list (per prime) of theta rows mod p, SAME shape.
    Normalize each by a common pivot (first nonzero of highest theta row),
    CRT + ratrec -> Fraction rows, or None."""
    shp = [(len(r), [len(x) for x in r]) for r in rows_by_p]
    if len(set(map(str, shp))) != 1:
        return None
    K = len(rows_by_p[0]) - 1
    piv = next((j for j, c in enumerate(rows_by_p[0][K]) if c), None)
    if piv is None or any(rw[K][piv] == 0 for rw in rows_by_p):
        return None
    norm = []
    for rw, p in zip(rows_by_p, primes):
        inv = pow(rw[K][piv], p - 2, p)
        norm.append([[c * inv % p for c in row] for row in rw])
    M = 1
    for p in primes:
        M *= p
    out = []
    for k in range(K + 1):
        row = []
        for j in range(len(norm[0][k])):
            c, mm = 0, 1
            for rw, p in zip(norm, primes):
                # CRT combine
                t = (rw[k][j] - c) * pow(mm % p, p - 2, p) % p
                c, mm = c + mm * t, mm * p
            q = _ratrec(c % M, M)
            if q is None:
                return None
            row.append(q)
        out.append(row)
    return int_rows(out)


def int_rows(rows):
    """Fraction theta rows -> primitive integer rows, sign: lead of top row >0."""
    den = 1
    for r in rows:
        for c in r:
            den = den * Fr(c).denominator // gcd(den, Fr(c).denominator)
    g = 0
    ir = [[int(Fr(c) * den) for c in r] for r in rows]
    for r in ir:
        for c in r:
            g = gcd(g, abs(c))
    g = g or 1
    ir = [[c // g for c in r] for r in ir]
    lead = next((c for c in reversed(ir[-1]) if c), 1)
    if lead < 0:
        ir = [[-c for c in r] for r in ir]
    return [[Fr(c) for c in r] for r in ir]


def apply_theta_exact(rows, a):
    """coefficient list of L(sum a_n z^n): out[n] = sum_j Qj(n-j) a_{n-j}."""
    K = len(rows) - 1
    smax = max(len(r) for r in rows) - 1
    out = []
    for n in range(len(a)):
        s = Fr(0)
        for j in range(min(n, smax) + 1):
            q = sum((rows[k][j] if j < len(rows[k]) else 0) * Fr(n - j) ** k
                    for k in range(K + 1))
            s += q * a[n - j]
        out.append(s)
    return out


# ------------- exact Ore right-division over Q (D-form) --------------------
def _pmulQ(a, b):
    out = [Fr(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return out


def _pdiff(a):
    return [Fr(i) * a[i] for i in range(1, len(a))] or [Fr(0)]


def _ptrim(a):
    while len(a) > 1 and a[-1] == 0:
        a = a[:-1]
    return a


def _Dk_times_Q(R, k):
    """D^k * R over Q, D-form lists."""
    r = len(R) - 1
    dR = [[list(c) for c in R]]
    for _ in range(k):
        dR.append([_pdiff(c) for c in dR[-1]])
    out = [[Fr(0)] for _ in range(r + k + 1)]
    for m in range(r + k + 1):
        for s in range(k + 1):
            j = m - k + s
            if 0 <= j <= r:
                t = [Fr(comb(k, s)) * c for c in dR[s][j]]
                n = max(len(out[m]), len(t))
                out[m] = [(out[m][i] if i < len(out[m]) else Fr(0)) +
                          (t[i] if i < len(t) else Fr(0)) for i in range(n)]
    return out


def ore_rightdiv_Q(L, R):
    """Exact pseudo right-division over Q, D-form (port of ore_rdiv_modp):
    lc(R)^e L = Q R + S. Returns (Q, S, e); S all-zero <=> R | L on the right."""
    L = [[Fr(c) for c in row] for row in L]
    R = [[Fr(c) for c in row] for row in R]
    n, r = len(L) - 1, len(R) - 1
    lcR = _ptrim(R[r])
    Q = [[Fr(0)] for _ in range(n - r + 1)]
    cur = [list(c) for c in L]
    e = 0
    for k in range(n - r, -1, -1):
        top = _ptrim(cur[r + k]) if r + k < len(cur) else [Fr(0)]
        if top == [Fr(0)] or not any(top):
            continue
        cur = [_pmulQ(c, lcR) for c in cur]
        Q = [_pmulQ(q, lcR) for q in Q]
        e += 1
        DkR = _Dk_times_Q(R, k)
        for m in range(len(DkR)):
            t = _pmulQ(top, DkR[m])
            nn = max(len(cur[m]), len(t))
            cur[m] = [(cur[m][i] if i < len(cur[m]) else Fr(0)) -
                      (t[i] if i < len(t) else Fr(0)) for i in range(nn)]
        nq = max(len(Q[k]), len(top))
        Q[k] = [(Q[k][i] if i < len(Q[k]) else Fr(0)) +
                (top[i] if i < len(top) else Fr(0)) for i in range(nq)]
        assert not any(_ptrim(cur[r + k])), "top not killed"
    S = [_ptrim(c) for c in cur[:r]]
    return Q, S, e


def theta_to_D_Q(rows):
    """exact theta rows -> D-form over Q."""
    K = len(rows) - 1
    S2 = [[0] * (K + 1) for _ in range(K + 1)]
    S2[0][0] = 1
    for nn in range(1, K + 1):
        for j in range(1, nn + 1):
            S2[nn][j] = S2[nn - 1][j - 1] + j * S2[nn - 1][j]
    out = []
    for j in range(K + 1):
        acc = [Fr(0)] * (max(len(r) for r in rows) + j)
        for k in range(j, K + 1):
            if S2[k][j]:
                for m, c in enumerate(rows[k]):
                    acc[m + j] += S2[k][j] * Fr(c)
        out.append(_ptrim(acc))
    while len(out) > 1 and not any(out[-1]):
        out.pop()
    return out


# --------- affine MUM tower (resonance + seed freedom carried) --------------
def _rref_constraints(M, p):
    """RREF of 4x4 M mod p -> (row-ops matrix E (4x4), pivot cols, zero-rows).
    E@M = R with R in RREF; zero-rows of R give constraint functionals E[i]."""
    A = [[M[i][j] % p for j in range(4)] for i in range(4)]
    E = [[1 if i == j else 0 for j in range(4)] for i in range(4)]
    prow = 0
    piv = []
    for col in range(4):
        r = next((i for i in range(prow, 4) if A[i][col]), None)
        if r is None:
            continue
        A[prow], A[r] = A[r], A[prow]
        E[prow], E[r] = E[r], E[prow]
        inv = pow(A[prow][col], p - 2, p)
        A[prow] = [x * inv % p for x in A[prow]]
        E[prow] = [x * inv % p for x in E[prow]]
        for i in range(4):
            if i != prow and A[i][col]:
                f = A[i][col]
                A[i] = [(A[i][j] - f * A[prow][j]) % p for j in range(4)]
                E[i] = [(E[i][j] - f * E[prow][j]) % p for j in range(4)]
        piv.append(col)
        prow += 1
    return A, E, piv, list(range(prow, 4))


def mum_tower_affine_modp(rows, p, N):
    """Affine family of log^3 formal solutions at 0 mod p.
    Returns (cols, why): cols = list over n of list over l of 4-vectors;
    l=0 particular (seed e_3), l>=1 free parameters (seed e_0..e_2 +
    resonance freedom).  Constraints eliminate parameters in place."""
    Qj = _Qj_mats(rows, p)
    smax = len(Qj) - 1
    ind = indicial_modp(rows, p)
    m0 = 0
    while m0 < len(ind) and ind[m0] % p == 0:
        m0 += 1
    if m0 < 4:
        return None, f"indicial 0-mult {m0} < 4"
    A = [[[0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]]]
    for n in range(1, N):
        k1 = len(A[0])
        rhs = [[0] * 4 for _ in range(k1)]
        for j in range(1, min(n, smax) + 1):
            Mj = _mpolyval(Qj[j], n - j, p)
            for l in range(k1):
                v = A[n - j][l]
                for i in range(4):
                    rhs[l][i] = (rhs[l][i] - sum(Mj[i][t] * v[t]
                                                 for t in range(4))) % p
        M = _mpolyval(Qj[0], n, p)
        if _poly_eval(ind, n, p) != 0:
            new = []
            for l in range(k1):
                x = [0] * 4
                for i in range(3, -1, -1):
                    s = sum(M[i][j2] * x[j2] for j2 in range(i + 1, 4)) % p
                    x[i] = (rhs[l][i] - s) * pow(M[i][i], p - 2, p) % p
                new.append(x)
            A.append(new)
            continue
        R, E, piv, zrows = _rref_constraints(M, p)
        for zr in zrows:                     # constraints g.(1,c)=0
            g = [sum(E[zr][i] * rhs[l][i] for i in range(4)) % p
                 for l in range(k1)]
            if not any(g):
                continue
            ls = next((l for l in range(k1 - 1, 0, -1) if g[l]), None)
            if ls is None:
                return None, f"inconsistent resonance at n={n}"
            inv = pow(g[ls], p - 2, p)
            for m in range(len(A)):          # substitute c_ls, drop column
                for l in range(k1):
                    if l != ls and g[l]:
                        f = g[l] * inv % p
                        A[m][l] = [(A[m][l][i] - f * A[m][ls][i]) % p
                                   for i in range(4)]
                del A[m][ls]
            for l in range(k1):
                if l != ls and g[l]:
                    f = g[l] * inv % p
                    rhs[l] = [(rhs[l][i] - f * rhs[ls][i]) % p
                              for i in range(4)]
            del rhs[ls]
            k1 -= 1
        # particular per column via RREF rows (consistent now)
        new = []
        for l in range(k1):
            b = [sum(E[i][t] * rhs[l][t] for t in range(4)) % p
                 for i in range(4)]
            x = [0] * 4
            for i, col in enumerate(piv):
                x[col] = (b[i] - sum(R[i][j2] * x[j2]
                                     for j2 in range(col + 1, 4))) % p
            new.append(x)
        # new freedom: nullspace of M -> fresh columns (zero history)
        free_cols = [c for c in range(4) if c not in piv]
        for fc in free_cols:
            v = [0] * 4
            v[fc] = 1
            for i, col in enumerate(piv):
                v[col] = (-R[i][fc]) % p
            for m in range(len(A)):
                A[m].append([0, 0, 0, 0])
            new.append(v)
        A.append(new)
    return A, None


def nullspace_basis_modp(A, p, cap=8):
    """Full nullspace of A mod p: (basis, dim).  dim = exact nullity; basis
    = one RREF vector per free column, TRAILING free cols first (basis[0]
    == nullspace_mod_fast's vector), at most cap kept.  Same elimination
    as annihilator.nullspace_mod_fast (single-product int64, p<2^31 safe)."""
    A = np.array(A, dtype=np.int64) % p
    m, n = A.shape
    pivs = []
    row = 0
    for col in range(n):
        sel = next((r for r in range(row, m) if A[r, col]), -1)
        if sel < 0:
            continue
        if sel != row:
            A[[row, sel]] = A[[sel, row]]
        inv = pow(int(A[row, col]), p - 2, p)
        A[row] = (A[row] * inv) % p
        mask = np.ones(m, dtype=bool)
        mask[row] = False
        fac = A[mask, col].copy()
        nz = fac != 0
        idx = np.where(mask)[0][nz]
        if len(idx):
            A[idx] = (A[idx] - fac[nz, None] * A[row]) % p
        pivs.append(col)
        row += 1
        if row == m:
            break
    pset = set(pivs)
    free = [c for c in range(n) if c not in pset]
    basis = []
    for f in reversed(free):
        if len(basis) >= cap:
            break
        v = np.zeros(n, dtype=np.int64)
        v[f] = 1
        for i, pc in enumerate(pivs):
            v[pc] = (-int(A[i, f])) % p
        basis.append(v.tolist())
    return basis, len(free)


def _affine_fit_matrix(cols, p, dR, order=4, verify_rows=30):
    """Lambda-closure linearized fit matrix (16N x (order+1)(dR+1)k1),
    unknowns u^l_ij = r_ij*(1,c)_l (see fit_order4_affine).  (A,k1)|None."""
    N = len(cols)
    k1 = len(cols[0])
    ncols = (order + 1) * (dR + 1) * k1
    if 16 * N < ncols + verify_rows:
        return None
    fac = [[1, 1, 1, 1], [1, 2, 3, 0], [2, 6, 0, 0], [6, 0, 0, 0]]
    A = np.zeros((16 * N, ncols), dtype=np.int64)
    for s in range(4):
        # Lam^s y: comp m <- ((m+s)!/m!) a_{m+s}
        pw = []
        for n in range(N):
            X = _mat_theta(n, p)
            pn = []
            for l in range(k1):
                v0 = [fac[s][m] * cols[n][l][m + s] % p if m + s < 4 else 0
                      for m in range(4)]
                vs = [v0]
                for i in range(order):
                    v = vs[-1]
                    vs.append([sum(X[r][c] * v[c] for c in range(4)) % p
                               for r in range(4)])
                pn.append(vs)
            pw.append(pn)
        for n in range(N):
            for j in range(min(dR, n) + 1):
                for i in range(order + 1):
                    for l in range(k1):
                        col = (i * (dR + 1) + j) * k1 + l
                        for m in range(4):
                            A[4 * (s * N + n) + m, col] = pw[n - j][l][i][m]
    return A, k1


def _U_of(v, order, dR, k1, p):
    """u-vector -> U matrix, rows a=(i,j), cols l."""
    return [[int(v[(i * (dR + 1) + j) * k1 + l]) % p for l in range(k1)]
            for i in range(order + 1) for j in range(dR + 1)]


def _rank1_col(U, p):
    """U rank<=1 and nonzero -> a nonzero column index l0, else None."""
    i0, l0 = next(((a, b) for a in range(len(U)) for b in range(len(U[0]))
                   if U[a][b]), (None, None))
    if i0 is None:
        return None
    for l in range(len(U[0])):
        for a in range(len(U)):
            if (U[a][l] * U[i0][l0] - U[i0][l] * U[a][l0]) % p:
                return None
    return l0


def fit_order4_affine(cols, p, dR, order=4, verify_rows=30, info=None):
    """Linearized fit with LAMBDA-CLOSURE rows: find R (order<=order,
    zdeg<=dR) and params c with R.(Lam^s (y_0 + sum c_l y_l)) = 0 for
    s = 0..3; genuine right factor <=> rank-1 U = R x (1,c) (column-wise
    syzygy tuples, a measured rank-1 killer on control families, are not).
    Uses the full nullspace basis (cap 192) + rank-1 search IN THE
    SPAN — an arbitrary nullvector is generically NOT rank-1 when dim>1.
    Primary: structured u0-split (_rank1_u0_split; measured control-family
    nullspaces are 20+-dim parameter-column syzygy modules, k_u0=0); fallback: dim<=3
    exact minor variety.  info dict gets nulldim[dR], k_u0, rank1_via,
    n_rank1, basis_truncated.
    Returns (Rrows_modp, l0) | ('exists', None) | None."""
    built = _affine_fit_matrix(cols, p, dR, order, verify_rows)
    if built is None:
        return None
    A, k1 = built
    basis, dim = nullspace_basis_modp(A, p, cap=192)
    if info is not None:
        info.setdefault("nulldim", {})[dR] = dim
    if not basis:
        return None
    Us = [_U_of(v, order, dR, k1, p) for v in basis[:8]]
    cands = [U for U in Us if _rank1_col(U, p) is not None]
    how = "basis"
    if not cands and len(basis) == dim:
        st = _rank1_u0_split(basis, order, dR, k1, p, info=info)
        if st:
            cands, how = st, "u0-split"
        if not cands and 2 <= dim <= 3:
            cands, how = _rank1_span_search(Us[:dim], p), "span"
    elif not cands and info is not None:
        info["basis_truncated"] = dim
    if not cands:
        res = (A @ np.array(basis[0], dtype=object)) % p
        if any(int(x) % p for x in res):
            return None                    # integrity guard on RREF
        return ("exists", None)
    nR = (order + 1) * (dR + 1)
    U = min(cands, key=lambda u: tuple(
        sorted((a, l) for a in range(nR) for l in range(k1) if u[a][l])))
    res = (A @ np.array([U[a][l] for a in range(nR) for l in range(k1)],
                        dtype=object)) % p
    if any(int(x) % p for x in res):
        return None                        # selected element must be exact
    l0 = _rank1_col(U, p)
    if info is not None:
        info["rank1_via"] = how
        info["n_rank1"] = len(cands)
    rows = [[U[i * (dR + 1) + j][l0] for j in range(dR + 1)]
            for i in range(order + 1)]
    while len(rows) > 1 and not any(rows[-1]):
        rows.pop()
    return rows, l0


# --------- rank-1-in-span minor-variety solver (stage-B fix) ----------------
def _pt_(f):
    while len(f) > 1 and f[-1] == 0:
        f.pop()
    return f


def _pmulp(f, g, p):
    out = [0] * (len(f) + len(g) - 1)
    for i, x in enumerate(f):
        if x:
            for j, y in enumerate(g):
                out[i + j] = (out[i + j] + x * y) % p
    return _pt_(out)


def _pdivmodp(f, g, p):
    f = list(f)
    g = _pt_(list(g))
    dg = len(g) - 1
    inv = pow(g[-1], p - 2, p)
    q = [0] * max(1, len(f) - dg)
    while len(_pt_(f)) - 1 >= dg and any(f):
        c = f[-1] * inv % p
        sh = len(f) - 1 - dg
        if c:
            q[sh] = c
            for i, gc in enumerate(g):
                f[sh + i] = (f[sh + i] - c * gc) % p
        f.pop()
    return _pt_(q), _pt_(f or [0])


def _pgcdp(f, g, p):
    f, g = _pt_(list(f)), _pt_(list(g))
    while any(g):
        f, g = g, _pdivmodp(f, g, p)[1]
    if any(f) and f[-1] != 1:
        inv = pow(f[-1], p - 2, p)
        f = [c * inv % p for c in f]
    return f


def _ppowmodp(b, e, m, p):
    r = [1]
    b = _pdivmodp(b, m, p)[1]
    while e:
        if e & 1:
            r = _pdivmodp(_pmulp(r, b, p), m, p)[1]
        b = _pdivmodp(_pmulp(b, b, p), m, p)[1]
        e >>= 1
    return r


def _proots(f, p, rng):
    """All roots of f in F_p (Cantor-Zassenhaus on the linear-factor part)."""
    f = _pt_(list(f))
    if len(f) < 2 or not any(f):
        return []
    xp = _ppowmodp([0, 1], p, f, p)          # x^p mod f
    xpx = list(xp) + [0] * (2 - len(xp))
    xpx[1] = (xpx[1] - 1) % p                # x^p - x mod f
    g = _pgcdp(f, xpx, p) if any(_pt_(list(xpx))) else \
        [c * pow(f[-1], p - 2, p) % p for c in f]
    roots, stack, guard = [], [g], 0
    while stack and guard < 512:
        h = stack.pop()
        guard += 1
        if len(h) < 2:
            continue
        if len(h) == 2:
            roots.append((-h[0]) * pow(h[1], p - 2, p) % p)
            continue
        for _ in range(64):
            a = rng.randrange(p)
            t = _ppowmodp([a, 1], (p - 1) // 2, h, p)
            t = list(t)
            t[0] = (t[0] - 1) % p
            d = _pgcdp(h, _pt_(t), p)
            if 0 < len(d) - 1 < len(h) - 1:
                stack.append(d)
                stack.append(_pdivmodp(h, d, p)[0])
                break
    return sorted(set(roots))


def _paddp(f, g, p, sub=False):
    n = max(len(f), len(g))
    s = -1 if sub else 1
    return _pt_([((f[i] if i < len(f) else 0) +
                  s * (g[i] if i < len(g) else 0)) % p for i in range(n)])


def _pencil_minor(M0, M1, e1, e2, p):
    """2x2 minor of M0 + t*M1 at entries e1=(a,l), e2=(a2,l2), as t-poly."""
    (a, l), (a2, l2) = e1, e2
    x0, x1 = M0[a][l], M1[a][l]
    y0, y1 = M0[a2][l2], M1[a2][l2]
    u0, u1 = M0[a][l2], M1[a][l2]
    w0, w1 = M0[a2][l], M1[a2][l]
    return _pt_([(x0 * y0 - u0 * w0) % p,
                 (x0 * y1 + x1 * y0 - u0 * w1 - u1 * w0) % p,
                 (x1 * y1 - u1 * w1) % p])


def _pencil_rank1(M0, M1, p, rng):
    """Verified rank-1 elements of {M0 + t M1 : t in F_p} U {M1}.
    Complete over P^1: gcd of a subset of minor polys can only ENLARGE the
    root set (true common roots survive any subset); every candidate is
    fully verified by _rank1_col."""
    nR, k1 = len(M0), len(M0[0])
    ents = [(a, l) for a in range(nR) for l in range(k1)
            if M0[a][l] or M1[a][l]]
    g, cnt, it = None, 0, 0
    for i in range(len(ents)):
        for j in range(i + 1, len(ents)):
            it += 1
            if it > 60000 or cnt > 400 or g == [1]:
                break
            if ents[i][0] == ents[j][0] or ents[i][1] == ents[j][1]:
                continue
            q = _pencil_minor(M0, M1, ents[i], ents[j], p)
            if not any(q):
                continue
            g = q if g is None else _pgcdp(g, q, p)
            cnt += 1
        if it > 60000 or cnt > 400 or g == [1]:
            break
    if g is None:
        ts = [0]                    # every minor identically 0: pencil rank<=1
    elif len(g) == 1:
        ts = []                     # constant gcd: no common root
    else:
        ts = _proots(g, p, rng)
    out = []
    for t in ts:
        U = [[(M0[a][l] + t * M1[a][l]) % p for l in range(k1)]
             for a in range(nR)]
        if any(any(r) for r in U) and _rank1_col(U, p) is not None:
            out.append(U)
    if _rank1_col(M1, p) is not None:
        out.append([row[:] for row in M1])
    return out


def _tri_minor(Ms, e1, e2, p):
    """2x2 minor of M0 + s M1 + t M2 at (e1,e2): [C0(s),C1(s),C2(s)] asc-t."""
    (a, l), (a2, l2) = e1, e2
    X = [Ms[k][a][l] for k in range(3)]
    Y = [Ms[k][a2][l2] for k in range(3)]
    U = [Ms[k][a][l2] for k in range(3)]
    W = [Ms[k][a2][l] for k in range(3)]

    def prod(P, Q):
        return ([P[0] * Q[0] % p, (P[0] * Q[1] + P[1] * Q[0]) % p,
                 P[1] * Q[1] % p],
                [(P[0] * Q[2] + P[2] * Q[0]) % p,
                 (P[1] * Q[2] + P[2] * Q[1]) % p],
                [P[2] * Q[2] % p])
    A1, A2 = prod(X, Y), prod(U, W)
    return [_paddp(A1[k], A2[k], p, sub=True) for k in range(3)]


def _det_poly(M, p):
    if len(M) == 1:
        return M[0][0]
    out = [0]
    for j in range(len(M)):
        if not any(M[0][j]):
            continue
        t = _pmulp(M[0][j], _det_poly([r[:j] + r[j + 1:] for r in M[1:]], p), p)
        if j % 2:
            t = [(-c) % p for c in t]
        out = _paddp(out, t, p)
    return out


def _res_t(m1, m2, p):
    """Resultant in t of two t-polys with s-poly coeffs; None if either deg 0."""
    c1 = [list(c) for c in m1]
    c2 = [list(c) for c in m2]
    while len(c1) > 1 and not any(c1[-1]):
        c1.pop()
    while len(c2) > 1 and not any(c2[-1]):
        c2.pop()
    d1, d2 = len(c1) - 1, len(c2) - 1
    if d1 == 0 or d2 == 0:
        return None
    n = d1 + d2
    S = [[[0] for _ in range(n)] for _ in range(n)]
    for r in range(d2):
        for k in range(d1 + 1):
            S[r][r + k] = c1[d1 - k]
    for r in range(d1):
        for k in range(d2 + 1):
            S[d2 + r][r + k] = c2[d2 - k]
    return _det_poly(S, p)


def _rank1_dim3(Ms, p, rng):
    """Verified rank-1 in span{M0,M1,M2}: chart (1,s,t) via resultants in t
    -> s-candidates -> pencil solve; s=inf charts via _pencil_rank1(M1,M2).
    Positive-dim rank-1 locus fallback: pencil solve at random s."""
    nR, k1 = len(Ms[0]), len(Ms[0][0])
    ents = [(a, l) for a in range(nR) for l in range(k1)
            if any(M[a][l] for M in Ms)]
    minors, sgcd, tries = [], None, 0
    while len(minors) < 10 and tries < 4000 and len(ents) >= 2:
        tries += 1
        e1, e2 = rng.sample(ents, 2)
        if e1[0] == e2[0] or e1[1] == e2[1]:
            continue
        m = _tri_minor(Ms, e1, e2, p)
        if not any(any(c) for c in m):
            continue
        if not any(m[1]) and not any(m[2]):      # t-free: constrains s only
            sgcd = m[0] if sgcd is None else _pgcdp(sgcd, m[0], p)
            continue
        minors.append(m)
    for i in range(len(minors) - 1):
        if sgcd == [1]:
            break
        r = _res_t(minors[i], minors[i + 1], p)
        if r is None or not any(r):
            continue                              # degenerate/identically 0
        sgcd = r if sgcd is None else _pgcdp(sgcd, r, p)
    out = []
    ss = []
    if sgcd is None:
        ss = [rng.randrange(p) for _ in range(3)]  # curve fallback
    elif len(sgcd) > 1:
        ss = _proots(sgcd, p, rng)
    for s in ss:
        M0s = [[(Ms[0][a][l] + s * Ms[1][a][l]) % p for l in range(k1)]
               for a in range(nR)]
        out += _pencil_rank1(M0s, Ms[2], p, rng)
    out += _pencil_rank1(Ms[1], Ms[2], p, rng)     # s = inf charts
    return out


def _rank1_span_search(Us, p, seed=20260712):
    """Verified rank-1 elements of span(Us), deduped by support.
    dim 2: complete (P^1 charts).  dim 3: complete up to random-choice
    degeneracy (2 GL(3) mixes).  dim>3: PARTIAL — basis triples only."""
    import random
    rng = random.Random(seed)
    D = len(Us)
    nR, k1 = len(Us[0]), len(Us[0][0])
    found = []
    if D == 2:
        found = _pencil_rank1(Us[0], Us[1], p, rng)
    elif D == 3:
        found = _rank1_dim3(Us, p, rng)
        if not found:                              # retry under random mix
            G = [[rng.randrange(1, p) for _ in range(3)] for _ in range(3)]
            Mx = [[[sum(G[r][k] * Us[k][a][l] for k in range(3)) % p
                    for l in range(k1)] for a in range(nR)] for r in range(3)]
            found = _rank1_dim3(Mx, p, rng)
    else:
        from itertools import combinations
        for tri in combinations(range(D), 3):
            found += _rank1_dim3([Us[t] for t in tri], p, rng)
            if found:
                break
    dedup, seen = [], set()
    for U in found:
        piv = next((U[a][l] for a in range(nR) for l in range(k1)
                    if U[a][l]), None)
        if piv is None:
            continue
        inv = pow(piv, p - 2, p)
        key = tuple(U[a][l] * inv % p for a in range(nR) for l in range(k1))
        if key not in seen:                  # ray identity, scale-invariant
            seen.add(key)
            dedup.append(U)
    return dedup


# --------- structured u0-split rank-1 search (large-dim cure) ---------------
def _rref_rows_modp(M, p):
    """RREF of list-of-int-rows mod p -> (nonzero rows, pivot cols)."""
    Mx = [row[:] for row in M]
    rr, piv = 0, []
    for c in range(len(Mx[0]) if Mx else 0):
        pr = next((i for i in range(rr, len(Mx)) if Mx[i][c]), None)
        if pr is None:
            continue
        Mx[rr], Mx[pr] = Mx[pr], Mx[rr]
        iv = pow(Mx[rr][c], p - 2, p)
        Mx[rr] = [x * iv % p for x in Mx[rr]]
        for i in range(len(Mx)):
            if i != rr and Mx[i][c]:
                f = Mx[i][c]
                Mx[i] = [(x - f * y) % p for x, y in zip(Mx[i], Mx[rr])]
        piv.append(c)
        rr += 1
    return Mx[:rr], piv


def _solve_affine_modp(Arows, rhs, p):
    """One solution x of A x = rhs mod p (free vars 0) or None.
    Via nullspace of [A|rhs]: trailing basis vector has last coord 1 iff
    the rhs column is free (= consistent); then x = -v[:-1]."""
    A = np.array(Arows, dtype=np.int64) % p
    b = (np.array(rhs, dtype=np.int64) % p).reshape(-1, 1)
    basis, _ = nullspace_basis_modp(np.hstack([A, b]), p, cap=1)
    n = A.shape[1]
    if basis and basis[0][n] % p == 1:
        return [(-c) % p for c in basis[0][:n]]
    return None


def _u0_try_r(B, D, r, order, dR, k1, p):
    """Completion for FIXED r: find x (nullspace coords) + c with
    u^0 = r, u^l = c_l r.  Returns verified rank-1 U or None."""
    nR = (order + 1) * (dR + 1)
    rows, rhs = [], []
    for a in range(nR):
        for l in range(k1):
            eq = [B[d][a * k1 + l] for d in range(D)] + [0] * (k1 - 1)
            if l:
                eq[D + l - 1] = (-r[a]) % p
            rows.append(eq)
            rhs.append(r[a] if l == 0 else 0)
    z = _solve_affine_modp(rows, rhs, p)
    if z is None:
        return None
    u = [sum(z[d] * B[d][col] for d in range(D)) % p
         for col in range(nR * k1)]
    U = _U_of(u, order, dR, k1, p)
    return U if _rank1_col(U, p) is not None else None


def _u0_pencil_ts(B, W, D, order, dR, k1, p, rng):
    """k=2 case: r(t) = W0 + t W1.  Numerically eliminate x once (x-columns
    are t-independent), leaving (G0 + t G1) c = h0 + t h1; consistent t are
    common roots of the (generic-rank)-minors of the augmented pencil
    [G|h](t).  Returns candidate t list (each later verified by full solve)."""
    nR = (order + 1) * (dR + 1)
    nc = D + 2 * (k1 - 1) + 2      # [x | C0 | C1 | h0 | h1]
    rows = []
    for a in range(nR):
        for l in range(k1):
            r = [0] * nc
            for d in range(D):
                r[d] = B[d][a * k1 + l]
            if l:
                r[D + (l - 1)] = (-W[0][a]) % p
                r[D + (k1 - 1) + (l - 1)] = (-W[1][a]) % p
            else:
                r[D + 2 * (k1 - 1)] = W[0][a]
                r[D + 2 * (k1 - 1) + 1] = W[1][a]
            rows.append(r)
    M = np.array(rows, dtype=np.int64)
    m = M.shape[0]
    rr = 0
    for col in range(D):               # pivots restricted to x-columns
        sel = next((i for i in range(rr, m) if M[i, col]), -1)
        if sel < 0:
            continue
        if sel != rr:
            M[[rr, sel]] = M[[sel, rr]]
        inv = pow(int(M[rr, col]), p - 2, p)
        M[rr] = (M[rr] * inv) % p
        mask = np.ones(m, dtype=bool)
        mask[rr] = False
        fac = M[mask, col].copy()
        nz = fac != 0
        idx = np.where(mask)[0][nz]
        if len(idx):
            M[idx] = (M[idx] - fac[nz, None] * M[rr]) % p
        rr += 1
    T = M[rr:]
    T = T[~np.any(T[:, :D] != 0, axis=1)]   # keep fully x-eliminated rows
    kc = k1 - 1
    P0 = np.hstack([T[:, D:D + kc], T[:, D + 2 * kc:D + 2 * kc + 1]])
    P1 = np.hstack([T[:, D + kc:D + 2 * kc], T[:, D + 2 * kc + 1:]])
    live = np.any(P0 != 0, axis=1) | np.any(P1 != 0, axis=1)
    P0, P1 = P0[live], P1[live]
    if P0.shape[0] == 0:
        return [rng.randrange(p)]           # unconstrained: generic t works
    def rank_at(t):
        Z = (P0 + t * P1) % p
        _, dnull = nullspace_basis_modp(Z, p, cap=0)
        return Z.shape[1] - dnull           # rank = ncols - nullity
    rho = max(rank_at(rng.randrange(p)) for _ in range(2))
    if rho == 0:
        return [rng.randrange(p)]
    g, folds, att = None, 0, 0
    nrowsT = P0.shape[0]
    while folds < 16 and att < 200 and g != [1]:
        att += 1
        ridx = rng.sample(range(nrowsT), min(rho, nrowsT))
        cidx = rng.sample(range(kc + 1), min(rho, kc + 1))
        if len(ridx) < rho or len(cidx) < rho:
            break
        sub = [[[int(P0[i][j]) % p, int(P1[i][j]) % p]
                for j in cidx] for i in ridx]
        dpoly = _det_poly([[_pt_(e) for e in row] for row in sub], p)
        if not any(dpoly) or len(dpoly) == 1:
            continue
        g = dpoly if g is None else _pgcdp(g, dpoly, p)
        folds += 1
    if g is None or not any(g):
        return [rng.randrange(p) for _ in range(3)]   # degenerate: sample
    if len(g) == 1:
        return []
    return _proots(g, p, rng)


def _rank1_u0_split(basis, order, dR, k1, p, info=None, seed=20260712):
    """MUM-structured rank-1 search in the nullspace span.  A genuine
    right factor kills y_0 + sum c_l y_l, so its tuple has u^0 = r != 0:
    quotient the (possibly huge) parameter-column syzygy module N_0 =
    {u : u^0 = 0} out by working with k = dim of the u^0-image W.
      k=0 : NO y_0-participating element -> no MUM factor at this dR
            (certificate; returns []).
      k=1 : r pinned up to scale -> completion is one LINEAR solve.
      k=2 : r(t) on P^1 -> consistency-minor pencil -> roots -> linear.
      k>=3: not implemented -> returns None (recorded in info).
    Every candidate verified rank-1; caller re-verifies residual."""
    import random
    rng = random.Random(seed)
    D = len(basis)
    nR = (order + 1) * (dR + 1)
    B = [[int(c) % p for c in w] for w in basis]
    W, _ = _rref_rows_modp([[B[d][a * k1] for a in range(nR)]
                            for d in range(D)], p)
    k = len(W)
    if info is not None:
        info["k_u0"] = k
    if k == 0:
        return []
    if k == 1:
        U = _u0_try_r(B, D, W[0], order, dR, k1, p)
        return [U] if U is not None else []
    if k == 2:
        out = []
        for t in _u0_pencil_ts(B, W, D, order, dR, k1, p, rng):
            r = [(W[0][a] + t * W[1][a]) % p for a in range(nR)]
            U = _u0_try_r(B, D, r, order, dR, k1, p)
            if U is not None:
                out.append(U)
        U = _u0_try_r(B, D, W[1], order, dR, k1, p)   # t = infinity chart
        if U is not None:
            out.append(U)
        return out
    return None
