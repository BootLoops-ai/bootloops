# baller transport engine — exact (no-float) certified-transport layer.
#!/usr/bin/env python3
"""
pipe_lib.py — EXACT layer of the certified transport, parametrized over
(operator, curve, flux, conventions) via a family card instead of DKMM
constants; generalized from a fixed-instance pilot.  No floats anywhere.

Provides (after configure(card, op_json, tower_json, gv_json, curve_json)):
  * PK / R_ORD / DEG          — theta-form operator, general order & degree
  * R_d, R_id, ind_m          — recurrences; indicial via exact integer-root
                                factorization (lead * prod (m - e_i))
  * PHI                       — restricted rho-jet tower Phi_alpha, |alpha|<=3
  * model_periods             — pinned-prepotential model side, general h:
                                Pi = (F0, F_a, X0, X^a) * w0 * v^3
  * solve_frame               — exact integral-frame dictionary
                                v^3 Pi_i = sum c_{alpha,p,q} v^p z3^q Phi_alpha
                                (fit prefix; held-out gate f1; L-annihilation
                                gate f2; xi-mutation gate f3 via run_pipe)
  * extend_tower              — exact recurrence extension with regression
All algorithms are the pilot's, with 2 -> h, 6 -> 2h+2, order/degree
parametric; for the DKMM card every output must equal the pilot's bank.
Run single-thread, exact fractions; ulimit -v 32505856, nice 5.
"""
from fractions import Fraction as Fr
from math import comb, factorial
import json, os

# ------------------------------------------------ configured globals
H = None          # h_eff
NU = None         # curve exponents
PK = None         # theta-form: {k: dense s-coeff list}
R_ORD = None      # operator order
DEG = None        # max s-degree
IND_LEAD = None   # indicial leading integer
IND_FAC = None    # [(root, mult)]
KAPPA = None; A_MAT = None; B_VEC = None; XI_COEF = None
GV = None         # {d tuple: int}
FFLUX = None; HFLUX = None
M_TOW = None; C_JET = None; BETAS = None; PHI = None
NPER = None       # 2h+2


def configure(card, op, tower, gv, curve):
    """op/tower/gv/curve = parsed JSON banks from the exact stages."""
    global H, NU, PK, R_ORD, DEG, IND_LEAD, IND_FAC, KAPPA, A_MAT, B_VEC
    global XI_COEF, GV, FFLUX, HFLUX, M_TOW, C_JET, BETAS, PHI, NPER
    H = card.h
    NU = tuple(curve["nu"])
    NPER = 2 * H + 2
    R_ORD = max(int(k.split(",")[1]) for k in op["theta_form"])
    DEG = max(int(k.split(",")[0]) for k in op["theta_form"])
    PK = {k: [0] * (DEG + 1) for k in range(R_ORD + 1)}
    for key, c in op["theta_form"].items():
        d, k = map(int, key.split(","))
        PK[k][d] = int(c)
    IND_LEAD = int(op["indicial"]["lead"])
    IND_FAC = [(int(e), int(m)) for e, m in op["indicial"]["factors"]]
    assert sum(m for _, m in IND_FAC) == R_ORD
    KAPPA = card.kappa
    A_MAT = card.a_mat
    B_VEC = card.b_vec
    XI_COEF = card.xi_coef
    GV = {tuple(map(int, k.split(","))): int(v) for k, v in gv["n_d"].items()}
    for d, v in card.gv_pinned.items():
        assert GV.get(d) == v, f"GV mismatch at {d}"
    FFLUX = list(curve["F"])
    HFLUX = list(curve["H"])
    M_TOW = tower["M"]
    C_JET = {}
    for key, lst in tower["C"].items():
        b = tuple(map(int, key.split(",")))
        C_JET[b] = [{tuple(map(int, kk.split(","))): Fr(v)
                     for kk, v in dd.items()} for dd in lst]
    BETAS = sorted(C_JET)
    PHI = {al: build_phi(al) for al in BETAS}


def kap(i, j, k):
    return KAPPA.get(tuple(sorted((i, j, k))), Fr(0))

# ------------------------------------------------ operator recurrences
def R_d(m, d):
    return sum(PK[k][d] * (m - d) ** k for k in range(R_ORD + 1) if PK[k][d])


def R_id(m, d, i):
    return sum(comb(k, i) * PK[k][d] * (m - d) ** (k - i)
               for k in range(i, R_ORD + 1) if PK[k][d])


def ind_m(m):
    v = R_d(m, 0)
    chk = IND_LEAD
    for e, mu in IND_FAC:
        chk *= (m - e) ** mu
    assert v == chk, f"indicial factorization mismatch at m={m}"
    return v


def falling(a, i):
    r = 1
    for t in range(i):
        r *= (a - t)
    return r

# ------------------------------------------------ jet tower -> Phi_alpha
def build_phi(al):
    """Phi_alpha as tower dict (lpow, e2, e3) -> coeff list [0..M_TOW]:
    Phi_a = a! sum_{b<=a} prod_c nu_c^{a_c-b_c}/(a_c-b_c)! l^{|a-b|} C_b
    (pilot convention: C_b already carries the b! factor)."""
    out = {}
    def rec(i, b):
        if i == H:
            bb = tuple(b)
            if bb not in C_JET:
                return
            pref = Fr(1)
            lp = 0
            for c in range(H):
                pref *= Fr(factorial(al[c]),
                           factorial(b[c]) * factorial(al[c] - b[c]))
                pref *= Fr(NU[c]) ** (al[c] - b[c])
                lp += al[c] - b[c]
            for m, dd in enumerate(C_JET[bb]):
                for (e2, e3), v in dd.items():
                    key = (lp, e2, e3)
                    if key not in out:
                        out[key] = [Fr(0)] * (M_TOW + 1)
                    out[key][m] += pref * v
            return
        for x in range(al[i] + 1):
            b.append(x); rec(i + 1, b); b.pop()
    rec(0, [])
    return out

# ------------------------------------------------ tower ops (l = log s)
def theta_tower(T, M):
    out = {}
    for (lp, P, Q), a in T.items():
        dst = out.setdefault((lp, P, Q), [Fr(0)] * (M + 1))
        for m in range(M + 1):
            if a[m]:
                dst[m] += m * a[m]
        if lp > 0:
            dst2 = out.setdefault((lp - 1, P, Q), [Fr(0)] * (M + 1))
            for m in range(M + 1):
                if a[m]:
                    dst2[m] += lp * a[m]
    return out


def apply_L(T, M):
    acc = {}
    cur = T
    for k in range(R_ORD + 1):
        if k > 0:
            cur = theta_tower(cur, M)
        for key, a in cur.items():
            for d in range(DEG + 1):
                c = PK[k][d]
                if c == 0:
                    continue
                dst = acc.setdefault(key, [Fr(0)] * (M + 1))
                for m in range(M + 1 - d):
                    if a[m]:
                        dst[m + d] += c * a[m]
    return acc


def tower_is_zero(T, Mchk):
    return all(a[m] == 0 for a in T.values() for m in range(Mchk + 1))

# ------------------------------------------- plain rational series helpers
def smul(x, y, M):
    out = [Fr(0)] * (M + 1)
    for i, xi in enumerate(x[:M + 1]):
        if xi:
            for j in range(0, M + 1 - i):
                if y[j]:
                    out[i + j] += xi * y[j]
    return out


def sinv(x, M):
    assert x[0] == 1
    out = [Fr(0)] * (M + 1); out[0] = Fr(1)
    for m in range(1, M + 1):
        out[m] = -sum(x[k] * out[m - k] for k in range(1, m + 1) if x[k])
    return out

def sexp(x, M):
    assert x[0] == 0
    out = [Fr(0)] * (M + 1); out[0] = Fr(1)
    for m in range(1, M + 1):
        out[m] = sum((Fr(k) * x[k] * out[m - k] for k in range(1, m + 1)
                      if x[k]), Fr(0)) / m
    return out


def spow(x, k, M):
    out = [Fr(0)] * (M + 1); out[0] = Fr(1)
    for _ in range(k):
        out = smul(out, x, M)
    return out

# --------------------------------- ring series: key (lpow, vpow, z3pow)
def rs_add(A, B, M, sgn=1):
    out = {k: v[:] for k, v in A.items()}
    for k, v in B.items():
        dst = out.setdefault(k, [Fr(0)] * (M + 1))
        for m in range(M + 1):
            if v[m]:
                dst[m] += sgn * v[m]
    return out


def rs_scale(A, M, c=Fr(1), dl=0, dP=0, dQ=0):
    out = {}
    for (lp, P, Q), v in A.items():
        out[(lp + dl, P + dP, Q + dQ)] = [c * x for x in v]
    return out


def rs_mul(A, B, M):
    out = {}
    for (l1, P1, Q1), va in A.items():
        for (l2, P2, Q2), vb in B.items():
            key = (l1 + l2, P1 + P2, Q1 + Q2)
            dst = out.setdefault(key, [Fr(0)] * (M + 1))
            for i, xi in enumerate(va):
                if xi:
                    for j in range(0, M + 1 - i):
                        if vb[j]:
                            dst[i + j] += xi * vb[j]
    return out


def rs_from_rat(x, M, lp=0, P=0, Q=0):
    return {(lp, P, Q): [x[m] if m < len(x) else Fr(0) for m in range(M + 1)]}


def rs_trim(A, M):
    return {k: v for k, v in A.items() if any(v[m] for m in range(M + 1))}

# --------------------------------------------- model periods on the curve
def model_periods(MF):
    """Exact ring-series (to order MF) of v^3 * (F0, F_a, X0, X^a) * w0 on the
    monomial curve, pinned frame; needs GV classes with nu.d <= MF (window
    sufficiency asserted by the caller)."""
    h = H
    z = (0,) * h
    def jet_rat(b):
        return [C_JET[b][m].get((0, 0), Fr(0)) for m in range(MF + 1)]
    def ebas(a):
        return tuple(1 if i == a else 0 for i in range(h))
    w0 = jet_rat(z)
    w0i = sinv(w0, MF)
    Sg = [smul(jet_rat(ebas(a)), w0i, MF) for a in range(h)]
    T = [rs_add(rs_from_rat([Fr(NU[a])], MF, lp=1, P=-1),
                rs_from_rat(Sg[a], MF, P=-1), MF) for a in range(h)]
    W0 = rs_from_rat(w0, MF)
    E = [sexp(Sg[a], MF) for a in range(h)]
    li3 = [Fr(0)] * (MF + 1)
    li2d = [[Fr(0)] * (MF + 1) for _ in range(h)]
    for d, nd in sorted(GV.items()):
        lead = sum(NU[a] * d[a] for a in range(h))
        if lead == 0 or lead > MF:
            continue
        qd = [Fr(1) if m == 0 else Fr(0) for m in range(MF + 1)]
        for a in range(h):
            qd = smul(qd, spow(E[a], d[a], MF), MF)
        qd = [qd[m - lead] if m >= lead else Fr(0) for m in range(MF + 1)]
        xk = qd[:]
        k = 1
        while k * lead <= MF:
            for m in range(MF + 1):
                if xk[m]:
                    li3[m] += Fr(nd, k ** 3) * xk[m]
                    for a in range(h):
                        li2d[a][m] += Fr(nd * d[a], k ** 2) * xk[m]
            xk = smul(xk, qd, MF)
            k += 1
    Finst = rs_from_rat([-x for x in li3], MF, P=-3)
    FinstA = [rs_from_rat([-x for x in li2d[a]], MF, P=-2) for a in range(h)]
    # F_a = -1/2 kap_abc t_b t_c + a_ab t_b + b_a + Finst_a
    TT = [[rs_mul(T[b], T[c], MF) for c in range(h)] for b in range(h)]
    FA = []
    for a in range(h):
        acc = FinstA[a]
        for b in range(h):
            for c in range(h):
                acc = rs_add(acc, rs_scale(TT[b][c], MF, c=-kap(a, b, c)/2), MF)
            acc = rs_add(acc, rs_scale(T[b], MF, c=A_MAT[a][b]), MF)
        acc = rs_add(acc, rs_from_rat([B_VEC[a]], MF), MF)
        FA.append(acc)
    # F = -1/6 kap ttt + 1/2 a tt + b t + xi + Finst
    Fc = rs_add(Finst, rs_from_rat([XI_COEF], MF, P=-3, Q=1), MF)
    for b in range(h):
        for c in range(h):
            for e in range(h):
                Fc = rs_add(Fc, rs_scale(rs_mul(TT[b][c], T[e], MF), MF,
                                         c=-kap(b, c, e)/6), MF)
            Fc = rs_add(Fc, rs_scale(TT[b][c], MF, c=A_MAT[b][c]/2), MF)
        Fc = rs_add(Fc, rs_scale(T[b], MF, c=B_VEC[b]), MF)
    F0 = rs_scale(Fc, MF, c=Fr(2))
    for a in range(h):
        F0 = rs_add(F0, rs_scale(rs_mul(T[a], FA[a], MF), MF, c=Fr(-1)), MF)
    P_list = [rs_mul(W0, F0, MF)] + [rs_mul(W0, FA[a], MF) for a in range(h)] \
        + [W0] + [rs_mul(W0, T[a], MF) for a in range(h)]
    return [rs_trim(rs_scale(Pi, MF, dP=3), MF) for Pi in P_list]

# ------------------------------------------------------------ frame solve
def solve_frame(M_FIT_LO=14, M_FIT_HI=20, L_CHK=177):
    """v^3 Pi_i = sum_{alpha,p,q} c^i_{alpha,p,q} v^p z3^q Phi_alpha, exact.
    Fit m <= M_FIT_LO; held-out to M_FIT_HI (gate f1); L_s-annihilation of the
    fitted towers to L_CHK (gate f2).  Returns (coeffs, towers)."""
    MP = model_periods(M_FIT_HI)
    unknowns = [(al, p, q) for al in BETAS for p in range(0, 4) for q in (0, 1)]

    def col_tower(u, M):
        al, p, q = u
        out = {}
        for (lp, e2, e3), v in PHI[al].items():
            key = (lp, p + 2 * e2, q + e3)
            dst = out.setdefault(key, [Fr(0)] * (M + 1))
            f = Fr(-1, 24) ** e2
            for m in range(M + 1):
                if v[m]:
                    dst[m] += f * v[m]
        return out

    cols = [col_tower(u, M_FIT_HI) for u in unknowns]
    coeffs, towers = [], []
    for i in range(NPER):
        rhs = MP[i]
        keys = sorted(set(k for c in cols for k in c) | set(rhs))
        rows, bvec = [], []
        for key in keys:
            for m in range(M_FIT_LO + 1):
                row = [c.get(key, None) for c in cols]
                row = [(r[m] if r else Fr(0)) for r in row]
                b = rhs.get(key, None)
                b = b[m] if b else Fr(0)
                if any(row) or b:
                    rows.append(row); bvec.append(b)
        sol = gauss_solve(rows, bvec, len(unknowns))
        for key in keys:
            for m in range(M_FIT_HI + 1):
                lhs = sum(sol[j] * cols[j].get(key, [Fr(0)]*(M_FIT_HI+1))[m]
                          for j in range(len(unknowns)) if sol[j])
                b = rhs.get(key, [Fr(0)]*(M_FIT_HI+1))[m]
                if lhs != b:
                    raise AssertionError(
                        f"frame gate f1 FAIL period {i} key {key} m {m}")
        Tfull = {}
        for j, u in enumerate(unknowns):
            if sol[j] == 0:
                continue
            ct = col_tower(u, M_TOW)
            for key, v in ct.items():
                dst = Tfull.setdefault(key, [Fr(0)] * (M_TOW + 1))
                for m in range(M_TOW + 1):
                    if v[m]:
                        dst[m] += sol[j] * v[m]
        Tfull = rs_trim(Tfull, M_TOW)
        res = apply_L(Tfull, M_TOW)
        assert tower_is_zero(res, L_CHK), f"frame gate f2 FAIL period {i}"
        coeffs.append({unknowns[j]: sol[j]
                       for j in range(len(unknowns)) if sol[j]})
        towers.append(Tfull)
    return coeffs, towers

def gauss_solve(rows, bvec, ncol):
    """exact Gaussian elimination; free vars -> 0; raises on inconsistency."""
    A = [list(r) + [bvec[i]] for i, r in enumerate(rows)]
    n = len(A); piv = []; rl = 0
    for c in range(ncol):
        pr = next((i for i in range(rl, n) if A[i][c] != 0), None)
        if pr is None:
            continue
        A[rl], A[pr] = A[pr], A[rl]
        pv = A[rl][c]
        A[rl] = [x / pv for x in A[rl]]
        for i in range(n):
            if i != rl and A[i][c] != 0:
                f = A[i][c]
                A[i] = [A[i][t] - f * A[rl][t] for t in range(ncol + 1)]
        piv.append(c); rl += 1
    for i in range(rl, n):
        if A[i][ncol] != 0:
            raise AssertionError("frame solve: inconsistent system")
    sol = [Fr(0)] * ncol
    for i, c in enumerate(piv):
        sol[c] = A[i][ncol]
    return sol

# ------------------------------------------------------------- extension
def extend_tower(T, N, reg_lo=None):
    """Extend every (P,Q)-graded slice of tower T from M_TOW to N terms via
    the exact recurrence; top-down in log-power within each grading.
    Regression: recurrence must reproduce banked terms on [reg_lo, M_TOW]."""
    if reg_lo is None:
        reg_lo = M_TOW - 40
    emax = max(e for e, _ in IND_FAC)
    assert reg_lo > emax + DEG, "regression window touches resonances"
    grades = sorted(set((P, Q) for (lp, P, Q) in T))
    out = {}
    for g in grades:
        layers = sorted([lp for (lp, P, Q) in T if (P, Q) == g], reverse=True)
        F = {lp: T[(lp,) + g] + [Fr(0)] * (N - M_TOW) for lp in layers}
        for lp in layers:
            higher = [l2 for l2 in layers if l2 > lp]
            def rhs_at(m):
                num = Fr(0)
                for d in range(1, DEG + 1):
                    if m - d >= 0:
                        num += R_d(m, d) * F[lp][m - d]
                for l2 in higher:
                    i = l2 - lp
                    ff = falling(l2, i)
                    for d in range(0, DEG + 1):
                        if m - d >= 0:
                            r = R_id(m, d, i)
                            if r:
                                num += ff * r * F[l2][m - d]
                return num
            for m in range(reg_lo, M_TOW + 1):
                assert ind_m(m) * F[lp][m] + rhs_at(m) == 0, \
                    f"recurrence regression FAIL lp={lp} g={g} m={m}"
            for m in range(M_TOW + 1, N + 1):
                F[lp][m] = -rhs_at(m) / ind_m(m)
        for lp in layers:
            out[(lp,) + g] = F[lp]
    return out


def bank_json(path, obj):
    with open(path, "w") as f:
        json.dump(obj, f)


def tower_to_json(T):
    return {f"{k[0]},{k[1]},{k[2]}": [str(x) for x in v] for k, v in T.items()}


def tower_from_json(d, N):
    out = {}
    for k, v in d.items():
        lp, P, Q = map(int, k.split(","))
        out[(lp, P, Q)] = [Fr(x) for x in v]
        assert len(out[(lp, P, Q)]) == N + 1
    return out
