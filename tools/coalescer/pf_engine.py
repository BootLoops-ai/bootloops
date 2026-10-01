#!/usr/bin/env python3
r"""
pf_engine.py — shared, ORDER-generic PF/Frobenius/transport machinery for the
banana coalescence ladder.  Self-contained (no imports from the K3 dir, so the
positive-control re-derivation is honest: SAME code path for K3 and CY3).

Conventions (matched to k3-119-funclevel/):
  z = 1/p^2 (MUM coord),   s = -p^2  =>  z = -1/s,   theta_z = z d/dz = -theta_t = -s d/ds.
  BFKNS holomorphic period:  varpi_0(z) = sum_n a_n z^n,  a_n = (-1)^n c_n,
    c_n = sum_{k1+..+kL=n} (n!/prod k_i!)^2 prod m_i^{2 k_i}.
  Threshold p^2=0 sits at z=inf (s=0); MUM at z=0 (s=inf).
"""
import numpy as np
import sympy as sp
import mpmath as mp
from fractions import Fraction as Fr
from math import factorial, comb

z = sp.symbols('z')
PRIMES = [2147483647, 2147483629, 2147483587, 2147483579, 2147483563,
          2147483549, 2147483543, 2147483497, 2147483489, 2147483477]


# ---------------------------------------------------------------------------
# EXACT PF operator via mod-p nullspace + CRT/ratrecon (build_cy4_collapse engine).
# ---------------------------------------------------------------------------
def wms_series_int(msq, N):
    """a_n = (-1)^n c_n as exact Python ints, c_n = sum_{k_1+..+k_L=n} (n!/prod k_i!)^2 prod m_i^{2k_i}.
    All-integer recursion, one mass at a time: (n!/(k_1!..k_L!))^2 = C(n,k_L)^2 ((n-k_L)!/(k_1!..k_{L-1}!))^2,
    so c^(L)_n = sum_j C(n,j)^2 m_L^j c^(L-1)_{n-j} with c^(1)_n = m_1^n (identical output to the former
    Fraction convolution, ~12x faster at N = 500)."""
    msq = [int(m) for m in msq]
    binsq = [[comb(n, j) ** 2 for j in range(n + 1)] for n in range(N + 1)]
    c = [msq[0] ** n for n in range(N + 1)]
    for m in msq[1:]:
        mpow = [m ** j for j in range(N + 1)]
        c = [sum(binsq[n][j] * mpow[j] * c[n - j] for j in range(n + 1)) for n in range(N + 1)]
    return [((-1) ** n) * c[n] for n in range(N + 1)]


def _np_nullspace(a_int, order, degz, p, extra=120):
    cols = [(j, k) for k in range(degz + 1) for j in range(order + 1)]
    ncol = len(cols)
    N = ncol + order + extra
    M = np.zeros((N + 1, ncol), dtype=np.int64)
    for ci, (j, k) in enumerate(cols):
        for m in range(N + 1):
            mm = m - k
            if 0 <= mm < len(a_int):
                if mm == 0 and j > 0:
                    val = 0
                else:
                    val = (a_int[mm] % p) * pow(mm % p, j, p) % p
                M[m, ci] = val % p
    nrow = M.shape[0]; r = 0; pivcols = []
    for c in range(ncol):
        piv = -1
        for i in range(r, nrow):
            if M[i, c] % p != 0:
                piv = i; break
        if piv < 0:
            continue
        M[[r, piv]] = M[[piv, r]]
        inv = pow(int(M[r, c]), p - 2, p)
        M[r] = (M[r] * inv) % p
        col = M[:, c].copy(); col[r] = 0
        nz = np.nonzero(col)[0]
        for i in nz:
            M[i] = (M[i] - M[i, c] * M[r]) % p
        pivcols.append(c); r += 1
        if r == nrow:
            break
    return ncol - r, cols, M, pivcols, r


def _ratrecon(a, m):
    if a % m == 0:
        return Fr(0)
    from math import isqrt as _isqrt  # integer sqrt: a float bound overflows here
    a %= m; N = _isqrt(m // 2)
    r0, r1 = m, a; s0, s1 = 0, 1
    while r1 > N:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if s1 != 0 and abs(s1) <= N:
        return Fr(r1, s1)
    return Fr(a - m) if a > m // 2 else Fr(a)


def find_minimal_op(msq, ord_max=12, degz_max=20, nprime=8):
    """Minimal-(order,degz) search; returns (order, degz, Pj, residual)."""
    Nser = (ord_max + 1) * (degz_max + 1) + ord_max + 150
    a_int = wms_series_int(msq, Nser)
    found = None
    for order in range(2, ord_max + 1):
        for degz in range(1, degz_max + 1):
            d1, *_ = _np_nullspace(a_int, order, degz, PRIMES[0])
            if d1 > 0:
                d2, *_ = _np_nullspace(a_int, order, degz, PRIMES[1])
                if d2 > 0:
                    found = (order, degz); break
        if found:
            break
    if not found:
        return None
    order, degz = found
    Pj = build_exact_op(a_int, order, degz, PRIMES[:nprime])
    res = verify_annihilation(a_int, Pj, min(len(a_int) - 2, (order + 1) * (degz + 1) + order + 80))
    return order, degz, Pj, res, a_int


def build_exact_op(a_int, order, degz, primes):
    cols = [(j, k) for k in range(degz + 1) for j in range(order + 1)]
    ncol = len(cols)
    sols = []
    for p in primes:
        dim, _, M, pivcols, rank = _np_nullspace(a_int, order, degz, p)
        if dim <= 0:
            return None
        free = [c for c in range(ncol) if c not in pivcols][0]
        v = np.zeros(ncol, dtype=object); v[free] = 1
        for ri, pc in enumerate(pivcols):
            v[pc] = int((-M[ri, free]) % p)
        sols.append((p, [int(x) % p for x in v]))
    P = 1
    for p, _ in sols:
        P *= p
    rec = []
    for ci in range(ncol):
        x = 0
        for p, v in sols:
            Mp = P // p
            x += v[ci] * Mp * pow(Mp, p - 2, p)
        x %= P
        rec.append(_ratrecon(x, P))
    L = 1
    for f in rec:
        L = L * f.denominator // sp.igcd(L, f.denominator)
    iv = [int(f * L) for f in rec]
    g = 0
    for x in iv:
        g = sp.igcd(g, x)
    if g:
        iv = [x // g for x in iv]
    Pj = {j: sp.Integer(0) for j in range(order + 1)}
    for ci, (j, k) in enumerate(cols):
        Pj[j] += iv[ci] * z ** k
    return Pj


def verify_annihilation(a_int, Pj, Ncheck):
    Pcoef = {}
    for j, pp in Pj.items():
        if pp == 0:
            continue
        cs = sp.Poly(pp, z).all_coeffs()[::-1]
        for k, c in enumerate(cs):
            Pcoef.setdefault(k, {})[j] = int(c)
    maxr = 0
    for m in range(Ncheck):
        s = 0
        for k, jc in Pcoef.items():
            mm = m - k
            if 0 <= mm < len(a_int):
                for j, cj in jc.items():
                    s += cj * a_int[mm] * (mm ** j if not (mm == 0 and j > 0) else 0)
        maxr = max(maxr, abs(s))
    return maxr


def indicial_threshold(Pj):
    """Indicial at z=inf (= threshold t=p^2=0) in t-exponent rho (theta_z=-theta_t)."""
    degz = max(sp.Poly(p, z).degree() for p in Pj.values() if p != 0)
    rho = sp.symbols('rho')
    ind = sum(sp.Poly(Pj[j], z).coeff_monomial(z ** degz) * (-rho) ** j
              for j in Pj if Pj[j] != 0)
    return sp.factor(ind), sp.roots(sp.Poly(ind, rho)), degz


def indicial_mum(Pj):
    def lowdeg(p):
        po = sp.Poly(p, z)
        return min(t for t in range(po.degree() + 1) if po.coeff_monomial(z ** t) != 0)
    md = min(lowdeg(p) for p in Pj.values() if p != 0)
    r = sp.symbols('r')
    ind = sum(sp.Poly(Pj[j], z).coeff_monomial(z ** md) * r ** j for j in Pj if Pj[j] != 0)
    return sp.factor(ind), sp.roots(sp.Poly(ind, r)), md


# ---------------------------------------------------------------------------
# Threshold Frobenius: clean power branch t^rho * sum a_n t^n  (exact-Q recursion).
# ---------------------------------------------------------------------------
def theta_t_terms(Pj, order):
    """Operator L = sum_j P_j(z) theta_z^j, z=1/t, theta_z=-theta_t.
    Returns list of (Fraction coeff, t_power=-z_power, theta_t_degree=j)."""
    terms = []
    for j in range(order + 1):
        if Pj.get(j, 0) == 0:
            continue
        for (k,), c in sp.Poly(Pj[j], z).terms():
            terms.append((Fr(int(c)), -int(k), j))
    return terms


def frob_power_branch(terms, rho, N):
    """f = sum_{n>=0} a_n t^{rho+n}, a_0=1. Exact-Q recursion. Returns (a, obstr)."""
    rho = sp.Rational(rho)
    a = {0: sp.Integer(1)}
    tp_min = min(tp for (_, tp, _) in terms)
    for target_n in range(1, N + 1):
        m = target_n + tp_min
        coeff_unknown = sp.Integer(0)
        known = sp.Integer(0)
        for (c, tp, jdeg) in terms:
            n = m - tp
            if n == target_n:
                coeff_unknown += sp.Rational(c) * (-(rho + n)) ** jdeg
            elif 0 <= n < target_n and n in a:
                known += sp.Rational(c) * (-(rho + n)) ** jdeg * a[n]
        if coeff_unknown == 0:
            return a, target_n
        a[target_n] = sp.cancel(-known / coeff_unknown)
    return a, None


# ---------------------------------------------------------------------------
# theta_z-state of a power-branch Phi_rho at Euclidean s_dec on a chosen branch
# of t = -s (arg t = +/- pi).  theta_z = -theta_t.
# ---------------------------------------------------------------------------
def phi_theta_state(a_coeffs, rho, s_dec, order, branch='lower', N=None):
    """Return col vector [theta_z^k Phi_rho]_{k=0..order-1} at s=s_dec, with
    Phi_rho(t) = t^rho sum_n a_n t^n,  t = s_dec * e^{i*arg},  arg = -pi (lower) / +pi (upper).
    All computed as theta_t^k then * (-1)^k.  Exact-Q a_n -> mp at CURRENT mp.mp.dps."""
    rho_q = sp.Rational(rho)
    arg = -mp.pi if branch == 'lower' else mp.pi
    t = mp.mpc(s_dec) * mp.expj(arg)
    keys = sorted(a_coeffs.keys())
    if N is not None:
        keys = [k for k in keys if k <= N]
    Y = mp.matrix(order, 1)
    for n in keys:
        v = sp.Rational(a_coeffs[n])
        cn = mp.mpf(v.p) / mp.mpf(v.q)
        ex = mp.mpf(rho_q.p) / mp.mpf(rho_q.q) + n
        term = cn * t ** ex
        exk = mp.mpf(1)
        for k in range(order):
            Y[k, 0] += exk * term
            exk *= ex
    # theta_z^k = (-1)^k theta_t^k
    for k in range(order):
        Y[k, 0] *= (-1) ** k
    return Y


# ---------------------------------------------------------------------------
# BFKNS seed state in mpmath (high-precision floats).
# ---------------------------------------------------------------------------
def wms_series_mp(msq, N):
    prod = [mp.mpf(0)] * (N + 1)
    prod[0] = mp.mpf(1)
    for mu in msq:
        mu = mp.mpf(mu)
        sfac = [mu ** k / (mp.factorial(k) ** 2) for k in range(N + 1)]
        new = [mp.mpf(0)] * (N + 1)
        for i in range(N + 1):
            if prod[i] == 0:
                continue
            pi = prod[i]
            for j in range(0, N + 1 - i):
                new[i + j] += pi * sfac[j]
        prod = new
    return [((-1) ** n) * prod[n] * (mp.factorial(n) ** 2) for n in range(N + 1)]


def make_period_state(msq, order, Nser):
    a_mp = wms_series_mp(list(msq), Nser)

    def state(sv):
        sv = mp.mpc(sv)
        zz = mp.mpf(-1) / sv
        Y = mp.matrix(order, 1)
        zn = mp.mpf(1)
        for n in range(len(a_mp)):
            term = a_mp[n] * zn
            nk = mp.mpf(1)
            for k in range(order):
                Y[k, 0] += nk * term
                nk *= n
            zn *= zz
        return Y
    return state


# ---------------------------------------------------------------------------
# Companion-system transport in s-coordinate, theta_z-state.
# ---------------------------------------------------------------------------
def build_ode_s(Pj, order):
    """dY/ds = A(s) Y,  Y=(f, theta_z f, ..., theta_z^{order-1} f),  z=-1/s.
    theta_z f_{k} = f_{k+1} for k<order-1; last row from L=0.
    A = -B/s where B_{k,k+1}=1, B_{order-1,j} = -P_j(-1/s)/P_order(-1/s)."""
    s = sp.symbols('s')
    Pj_s = {j: sp.cancel(Pj[j].subs(z, -sp.Integer(1) / s)) for j in Pj}
    Pord = Pj_s[order]
    B = sp.zeros(order, order)
    for k in range(order - 1):
        B[k, k + 1] = 1
    for j in range(order):
        B[order - 1, j] = sp.cancel(-Pj_s.get(j, 0) / Pord)
    A = sp.Matrix(order, order, lambda i, j: sp.cancel(-B[i, j] / s))
    n = order
    num = [[None] * n for _ in range(n)]
    den = [[None] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            e = sp.cancel(sp.together(A[i, j]))
            nn, dd = sp.fraction(e)
            num[i][j] = sp.Poly(sp.expand(nn), s)
            den[i][j] = sp.Poly(sp.expand(dd), s)
    # singularities in s: s=0 (threshold) + roots of P_order(z) at z=-1/s
    lead = sp.Poly(Pj[order], z)
    sings = [mp.mpc(0)]
    for r in sp.nroots(lead, n=50):
        rr = complex(r)
        if abs(rr) > 1e-40:
            sings.append(mp.mpc(-1) / mp.mpc(rr.real, rr.imag))
    return num, den, s, sings


def _shift_poly_coeffs(poly, s0, ordr):
    coeffs = poly.all_coeffs()[::-1]
    deg = len(coeffs) - 1
    out = [mp.mpc(0)] * (ordr + 1)
    for d in range(deg + 1):
        cd = mp.mpf(int(coeffs[d]))
        for k in range(min(d, ordr) + 1):
            out[k] += cd * mp.binomial(d, k) * s0 ** (d - k)
    return out


def _series_ratio(numc, denc, ordr):
    inv0 = 1 / denc[0]
    g = [mp.mpc(0)] * (ordr + 1)
    g[0] = inv0
    for k in range(1, ordr + 1):
        s = mp.mpc(0)
        for m in range(1, k + 1):
            dm = denc[m] if m < len(denc) else mp.mpc(0)
            s += dm * g[k - m]
        g[k] = -s * inv0
    r = [mp.mpc(0)] * (ordr + 1)
    for k in range(ordr + 1):
        acc = mp.mpc(0)
        for m in range(k + 1):
            nm = numc[m] if m < len(numc) else mp.mpc(0)
            acc += nm * g[k - m]
        r[k] = acc
    return r


def _A_taylor_at(num, den, s0, ordr):
    n = len(num)
    Ms = [mp.zeros(n, n) for _ in range(ordr + 1)]
    for i in range(n):
        for j in range(n):
            numc = _shift_poly_coeffs(num[i][j], s0, ordr)
            denc = _shift_poly_coeffs(den[i][j], s0, ordr)
            r = _series_ratio(numc, denc, ordr)
            for k in range(ordr + 1):
                Ms[k][i, j] = r[k]
    return Ms


def step_taylor(num, den, Y0, s0, h, ordr, ncols=1):
    n = len(num)
    Ms = _A_taylor_at(num, den, s0, ordr)
    a = [None] * (ordr + 2)
    a[0] = mp.matrix(Y0)
    for nn in range(ordr + 1):
        acc = mp.zeros(n, ncols)
        for k in range(nn + 1):
            acc += Ms[k] * a[nn - k]
        a[nn + 1] = acc / (nn + 1)
    Y = mp.zeros(n, ncols)
    hp = mp.mpc(1)
    for nn in range(ordr + 2):
        Y += a[nn] * hp
        hp *= h
    return Y


def transport_adaptive(num, den, Y0, waypoints, ordr, sings, frac=0.3, ncols=1):
    sings = [mp.mpc(s_) for s_ in sings]
    Y = mp.matrix(Y0)
    zc = mp.mpc(waypoints[0])
    for tgt in waypoints[1:]:
        tgt = mp.mpc(tgt)
        nstep = 0
        while True:
            d = abs(tgt - zc)
            dist = min(abs(zc - s_) for s_ in sings)
            nstep += 1
            if nstep > 200000:
                raise RuntimeError("step limit")
            if d <= frac * dist:
                Y = step_taylor(num, den, Y, zc, tgt - zc, ordr, ncols)
                zc = tgt; break
            h = (tgt - zc) * (frac * dist / d)
            Y = step_taylor(num, den, Y, zc, h, ordr, ncols)
            zc += h
    return Y


def digits(x, y):
    x = mp.mpc(x); y = mp.mpc(y)
    if x == y:
        return float(mp.mp.dps)
    return float(-mp.log10(abs(x - y) / (abs(y) + mp.mpf(10) ** (-mp.mp.dps))))
