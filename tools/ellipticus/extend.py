"""ellipticus.extend — exact extended inhomogeneous-PF systems for window
moments of a family of (hyper)elliptic charts, + exact local series seeds.

The construction is family-agnostic; window endpoints and the parameter
symbol are arguments, not module globals.  The moment-system machinery has
been gated at 133.5d (kernel) / 172.7d (two-precision) against independent
quadrature.

Objects: for a squarefree-in-x family F(x;z) (deg 3 or 4 in x, coefficients
in Q(z)) and an exact rational window [x1, x2],

    I_k(z) = int_{x1}^{x2} x^k dx / sqrt(F(x;z)),  k = 0..kmax,
    B_i(z) = F(x_i; z)^(-1/2)                       (boundary states),

satisfy a CLOSED first-order (kmax+3)x(kmax+3) system Y' = M(z) Y exactly over
Q(z) (Bezout split of -x^k F_z/2 = al F + be F' + integration by parts + the
exact-form moment reduction d[x^n sqrt F] with window boundary terms).  The
same construction with x1,x2 = branch points gives the homogeneous
Picard-Fuchs block for complete cycles (B_i columns then vanish).
"""
import json
from fractions import Fraction
import sympy as sp
import flint
import mpmath as mp

from .curve import exact


# ------------------------------------------------------------ system builder
def moment_system(Fpoly, kmax, x, z, X1, X2):
    """Exact first-order system for the window moments of y^2 = F(x;z).
    Fpoly: sympy Poly in x (coefficients rational in z); X1, X2: exact
    sympy Rationals (window endpoints).  Returns sympy Matrix M(z) with
    state (I_0..I_kmax, B_1, B_2).  Verbatim logic of the gated
    s1_extend.moment_system with (x, z, X1, X2) parameterized."""
    DOM = sp.QQ.frac_field(z)
    Fx = Fpoly.as_expr()
    dF = sp.diff(Fx, x)
    Fz = sp.diff(Fx, z)
    d = sp.degree(Fx, x)
    n = kmax + 3
    M = sp.zeros(n, n)
    Px = sp.Poly(Fx, x, domain=DOM)
    dPx = sp.Poly(dF, x, domain=DOM)
    sB, tB, g = Px.gcdex(dPx)
    assert g.degree() == 0, 'F not squarefree over Q(z)'
    sB = sp.Poly(sB.as_expr() / g.as_expr(), x, domain=DOM)
    tB = sp.Poly(tB.as_expr() / g.as_expr(), x, domain=DOM)

    topneed = d - 1 + kmax
    red = {}
    for m in range(kmax + 1):
        v = [sp.Integer(0)] * (kmax + 1)
        v[m] = sp.Integer(1)
        red[m] = (v, sp.Integer(0), sp.Integer(0))

    def reduce_m(m):
        if m in red:
            return red[m]
        nn_ = m - d + 1
        assert nn_ >= 0
        rel = sp.expand(nn_ * x ** (nn_ - 1) * Fx + x ** nn_ * dF / 2) \
            if nn_ > 0 else sp.expand(dF / 2)
        relp = sp.Poly(rel, x)
        lead = relp.coeff_monomial(x ** m)
        assert sp.simplify(lead) != 0
        v = [sp.Integer(0)] * (kmax + 1)
        b1 = -X1 ** nn_ * Fpoly.as_expr().subs(x, X1) / lead
        b2 = X2 ** nn_ * Fpoly.as_expr().subs(x, X2) / lead
        for mm in range(m - 1, -1, -1):
            cj = relp.coeff_monomial(x ** mm) if mm <= m - 1 else 0
            if cj == 0:
                continue
            vv, c1, c2 = reduce_m(mm)
            for t_ in range(kmax + 1):
                v[t_] -= cj * vv[t_] / lead
            b1 -= cj * c1 / lead
            b2 -= cj * c2 / lead
        red[m] = ([sp.cancel(e) for e in v], sp.cancel(b1), sp.cancel(b2))
        return red[m]

    for m in range(topneed + 1):
        reduce_m(m)

    for k in range(kmax + 1):
        Nkp = sp.Poly(sp.expand(-x ** k * Fz / 2), x, domain=DOM)
        be = (Nkp * tB).rem(Px)
        al = (Nkp - be * dPx).quo(Px)
        assert (al * Px + be * dPx - Nkp).is_zero, 'Bezout split failed'
        ie = sp.expand(al.as_expr() + 2 * sp.diff(be.as_expr(), x))
        dg = sp.degree(ie, x) if ie != 0 else -1
        for m in range(dg + 1):
            cm = sp.cancel(ie.coeff(x, m))
            if cm == 0:
                continue
            vv, c1, c2 = red[m]
            for t_ in range(kmax + 1):
                M[k, t_] += cm * vv[t_]
            M[k, kmax + 1] += cm * c1
            M[k, kmax + 2] += cm * c2
        M[k, kmax + 1] += 2 * be.as_expr().subs(x, X1)
        M[k, kmax + 2] += -2 * be.as_expr().subs(x, X2)
    for i, XI in enumerate((X1, X2)):
        M[kmax + 1 + i, kmax + 1 + i] = sp.cancel(
            -Fz.subs(x, XI) / (2 * Fx.subs(x, XI)))
    return sp.Matrix([[sp.cancel(sp.together(M[i, j])) for j in range(n)]
                      for i in range(n)])


def polyform(M, z):
    """Common-denominator polynomial form Den(z) Y' = NumM(z) Y as exact
    integer-string coefficient lists (descending), + singular inventory.
    In-memory; no cache file."""
    n = M.rows
    dens = []
    for i in range(n):
        for j in range(n):
            if M[i, j] != 0:
                dens.append(sp.denom(sp.together(M[i, j])))
    Den = sp.Integer(1)
    for dd in dens:
        Den = sp.lcm(Den, dd)
    Den = sp.Poly(Den, z)
    NumM = [[sp.Poly(sp.cancel(sp.together(M[i, j] * Den.as_expr())), z)
             if M[i, j] != 0 else None for j in range(n)] for i in range(n)]
    sing = set()
    for dd in dens:
        if sp.degree(dd, z) > 0:
            for f_, e_ in sp.factor_list(dd, z)[1]:
                if sp.degree(f_, z) > 0:
                    sing.add(sp.Poly(f_, z).monic().as_expr())
    return dict(
        n=n,
        Den=[str(c) for c in Den.all_coeffs()],
        NumM=[[([str(c) for c in NumM[i][j].all_coeffs()]
                if NumM[i][j] is not None else None)
               for j in range(n)] for i in range(n)],
        den_factors=[str(s) for s in sorted(sing, key=str)])


# ------------------------------------------------------ exact series seeds
def _fq(c):
    f = Fraction(c)
    return flint.fmpq(f.numerator, f.denominator)


def _fqp(coeffs_asc):
    return flint.fmpq_poly([_fq(c) for c in coeffs_asc])


def h_series(Fz, N):
    """g_n recurrence for h = (F/F_0)^(-1/2):
       (n+1) g_{n+1} = - sum_{i>=1} ((n+1) - i/2) F_i g_{n+1-i} F_0^{i-1},
    h_n = g_n / F_0^n, g_0 = 1.  Fz = [F_0, F_1, ...] flint fmpq_polys in x
    (F expanded in the LOCAL parameter).  Returns [g_0..g_N]."""
    F0 = Fz[0]
    F0pow = [flint.fmpq_poly([1])]
    for i in range(1, len(Fz)):
        F0pow.append(F0pow[-1] * F0)
    g = [flint.fmpq_poly([1])]
    for n in range(N):
        acc = flint.fmpq_poly([])
        for i in range(1, min(n + 1, len(Fz) - 1) + 1):
            if Fz[i].degree() < 0:
                continue
            cf = _fq(Fraction(2 * (n + 1) - i, 2))
            acc = acc + Fz[i] * g[n + 1 - i] * F0pow[i - 1] \
                * flint.fmpq_poly([cf])
        g.append(acc * flint.fmpq_poly([_fq(Fraction(-1, n + 1))]))
    return g


def _taylor_shift(p, a):
    cs = list(p.coeffs())[::-1] if p.degree() >= 0 else [flint.fmpq(0)]
    out = [cs[0]]
    for c in cs[1:]:
        new = [None] * (len(out) + 1)
        new[0] = out[0] * a + c
        for k in range(1, len(out)):
            new[k] = out[k] * a + out[k - 1]
        new[len(out)] = out[-1]
        out = new
    return flint.fmpq_poly(out)


def _series_inv(p, N):
    cs = list(p.coeffs())
    assert cs and cs[0] != 0
    inv = [1 / cs[0]]
    for n in range(1, N):
        s = flint.fmpq(0)
        for k in range(1, min(n, len(cs) - 1) + 1):
            s += cs[k] * inv[n - k]
        inv.append(-s / cs[0])
    return flint.fmpq_poly(inv)


def _poly_coeffs_frac(p, N):
    cs = list(p.coeffs())
    out = [Fraction(int(c.p), int(c.q)) for c in cs]
    while len(out) < N:
        out.append(Fraction(0))
    return out


def int_rational(numer, dfacs, lo, hi):
    """EXACT integral of numer(x)/prod (x-r)^m over [lo,hi]:
    numer flint fmpq_poly, dfacs [(Fraction r, int m)], lo/hi mp values at
    ambient precision; partial fractions EXACT in fmpq.  No pole in [lo,hi]
    (caller guarantees)."""
    den = flint.fmpq_poly([1])
    for r, m in dfacs:
        lin = flint.fmpq_poly([_fq(-r), _fq(1)])
        for _ in range(m):
            den = den * lin
    quo = numer // den
    rem = numer - quo * den
    val = mp.mpc(0)
    qc = _poly_coeffs_frac(quo, quo.degree() + 1) if quo.degree() >= 0 else []
    for k, c in enumerate(qc):
        cm = mp.mpf(c.numerator) / mp.mpf(c.denominator)
        val += cm * (hi ** (k + 1) - lo ** (k + 1)) / (k + 1)
    for (r, m) in dfacs:
        other = flint.fmpq_poly([1])
        for (r2, m2) in dfacs:
            if r2 == r:
                continue
            lin = flint.fmpq_poly([_fq(-r2), _fq(1)])
            for _ in range(m2):
                other = other * lin
        rs = _fq(r)
        remS = _taylor_shift(rem, rs)
        othS = _taylor_shift(other, rs)
        inv = _series_inv(othS, m)
        ser = remS * inv
        cs = _poly_coeffs_frac(ser, m)
        rm = mp.mpf(r.numerator) / mp.mpf(r.denominator)
        for j0 in range(m):
            aj = cs[j0]
            p = m - j0
            if aj == 0:
                continue
            cm = mp.mpf(aj.numerator) / mp.mpf(aj.denominator)
            if p == 1:
                val += cm * (mp.log(hi - rm) - mp.log(lo - rm))
            else:
                val += cm * ((hi - rm) ** (1 - p) - (lo - rm) ** (1 - p)) \
                    / (1 - p)
    return val


def family_zpolys(Fpoly, x, z, zbase, N):
    """F(x; zbase + t) expanded in t to order N: list of flint fmpq_polys in
    x, [F_0(x), F_1(x), ...].  Coefficients may be RATIONAL functions of z:
    F = num(x,z)/den(z); den(zb+t) is series-inverted
    exactly in fmpq.  zbase exact rational; den(zbase) != 0 required."""
    zb = sp.Rational(exact(zbase).numerator, exact(zbase).denominator)
    Fx = sp.together(sp.cancel(Fpoly.as_expr()))
    num, den = sp.fraction(Fx)
    if den.has(x):
        raise ValueError("family_zpolys: denominator depends on x")
    t = sp.Symbol('_t_local')
    dx = sp.degree(num, x)
    # numerator: per x-power, exact t-polynomial
    numt = sp.expand(num.subs(z, zb + t))
    dent = sp.expand(den.subs(z, zb + t))
    dp = sp.Poly(dent, t)
    dco = [Fraction(sp.Rational(c).p, sp.Rational(c).q)
           for c in dp.all_coeffs()[::-1]]
    assert dco[0] != 0, "family_zpolys: den(zbase) == 0"
    dinv = _series_inv(_fqp(dco), N + 1)          # 1/den as t-series
    dinv_c = _poly_coeffs_frac(dinv, N + 1)
    out_cols = []                                  # per x-power t-series
    for k in range(dx + 1):
        ck = sp.expand(sp.expand(numt).coeff(x, k))
        pk = sp.Poly(ck, t) if ck != 0 else None
        nco = ([Fraction(sp.Rational(c).p, sp.Rational(c).q)
                for c in pk.all_coeffs()[::-1]] if pk is not None else [])
        # multiply by dinv, truncate to N+1
        col = []
        for j in range(N + 1):
            s = Fraction(0)
            for i in range(min(j, len(nco) - 1) + 1 if nco else 0):
                s += nco[i] * dinv_c[j - i]
            col.append(s)
        out_cols.append(col)
    out = []
    for j in range(N + 1):
        out.append(_fqp([out_cols[k][j] for k in range(dx + 1)]))
    return out


def analytic_seed_series(Fpoly, x, z, zbase, X1, X2, kmax, N):
    """Exact z-Taylor of (I_0..I_kmax)(zbase + t) at an ANALYTIC base fibre.
    Two exact charts (both gated constructions):

    chart S (perfect square):
        F(x; zbase) = const * (poly(x))^2, poly root-free on [X1, X2] —
        the integrand's z-coefficients are RATIONAL in x; exact partial
        fractions (int_rational).
    chart SL (square x linear):
        F(x; zbase) = const * S(x)^2 * L(x), L = alpha x + beta linear,
        L > 0 on the window; substitution v = sqrt(L(x)) rationalizes:
        x = (v^2 - beta)/alpha; requires each root r_i of S to give a
        RATIONAL rho_i = sqrt(beta + alpha r_i) (else refuse).

    Returns dict k -> [c_0..c_N] (mp values at ambient precision).
    Raises ValueError (loud scope certificate) on any other base fibre."""
    zb = sp.Rational(exact(zbase).numerator, exact(zbase).denominator)
    F0 = sp.factor(sp.cancel(Fpoly.as_expr().subs(z, zb)))
    fl = sp.factor_list(F0, x)
    const = sp.Rational(fl[0])
    odd = [(b_, e_) for b_, e_ in fl[1] if e_ % 2 == 1]
    X1f, X2f = exact(X1), exact(X2)
    lo_q = sp.Rational(X1f.numerator, X1f.denominator)
    hi_q = sp.Rational(X2f.numerator, X2f.denominator)
    xm = (lo_q + hi_q) / 2
    Fz = family_zpolys(Fpoly, x, z, zbase, N)
    g = h_series(Fz, N)
    out = {k: [] for k in range(kmax + 1)}

    def sq_factors():
        """(lead, dfacs_base) of the square-root polynomial S(x) with
        F0 = const * S^2 * (odd part); linear factors only."""
        lead = sp.Integer(1)
        dfb = []
        for b_, e_ in fl[1]:
            if e_ < 2:
                continue
            pb = sp.Poly(b_, x)
            if pb.degree() == 1:
                cc = pb.all_coeffs()
                r_ = -sp.Rational(cc[1]) / sp.Rational(cc[0])
                lead = lead * sp.Rational(cc[0]) ** (e_ // 2)
                dfb.append((Fraction(r_.p, r_.q), e_ // 2))
            elif pb.degree() == 0:
                lead = lead * sp.Rational(b_) ** (e_ // 2)
            else:
                raise ValueError("non-linear repeated factor in base fibre: "
                                 "unsupported seed chart")
        return lead, dfb

    if not odd:
        # ---- chart S ----
        sq = sp.sqrt(const)
        if not sq.is_rational:
            raise ValueError("perfect-square constant not a rational square")
        lead, dfacs_base = sq_factors()
        lead = lead * sq
        # sign: sqrt(F0) positive on window
        sgn_probe = lead * sp.prod((xm - sp.Rational(r.numerator,
                                                     r.denominator)) ** m
                                   for r, m in dfacs_base)
        if sgn_probe < 0:
            lead = -lead
        leadF = Fraction(sp.Rational(lead).p, sp.Rational(lead).q)
        lo = mp.mpf(X1f.numerator) / mp.mpf(X1f.denominator)
        hi = mp.mpf(X2f.numerator) / mp.mpf(X2f.denominator)
        for n in range(N + 1):
            m = 2 * n + 1
            dfacs = [(r_, mi * m) for (r_, mi) in dfacs_base]
            sc = mp.mpf(1) / (mp.mpf(leadF.numerator)
                              / mp.mpf(leadF.denominator)) ** m
            for k in range(kmax + 1):
                numer = g[n] * flint.fmpq_poly([0] * k + [1])
                v = int_rational(numer, dfacs, lo, hi)
                out[k].append(mp.re(v) * sc)
        return out

    if len(odd) == 1 and sp.Poly(odd[0][0], x).degree() == 1 \
            and odd[0][1] == 1:
        # ---- chart SL ----
        Lb = sp.Poly(odd[0][0], x)
        alpha = sp.Rational(Lb.all_coeffs()[0])
        beta = sp.Rational(Lb.all_coeffs()[1])
        if const < 0:
            # fold the sign into L (e.g. -x^2(x-1) = x^2 (1-x))
            const, alpha, beta = -const, -alpha, -beta
        if alpha * xm + beta < 0:
            raise ValueError("linear factor negative on window: base fibre "
                             "not positive there")
        sq = sp.sqrt(const)
        if not sq.is_rational:
            raise ValueError("SL chart: const not a rational square")
        lead, dfacs_base = sq_factors()
        lead = lead * sq
        sgn_probe = lead * sp.prod((xm - sp.Rational(r.numerator,
                                                     r.denominator)) ** m
                                   for r, m in dfacs_base)
        if sgn_probe < 0:
            lead = -lead
        leadF = Fraction(sp.Rational(lead).p, sp.Rational(lead).q)
        # v-roots of the square factors: rho_i = sqrt(L(r_i)) rational
        vfacs_base = []          # [(rho or -rho Fractions, mult per (x-r))]
        for r_, mi in dfacs_base:
            Lr = alpha * sp.Rational(r_.numerator, r_.denominator) + beta
            rho = sp.sqrt(Lr)
            if not rho.is_rational:
                raise ValueError("SL chart: sqrt(L(r_i)) irrational — "
                                 "unsupported seed chart")
            rho = Fraction(sp.Rational(rho).p, sp.Rational(rho).q)
            # (x - r) = (v^2 - L(r))/alpha = (v-rho)(v+rho)/alpha
            vfacs_base.append((rho, mi))
        av = Fraction(alpha.p, alpha.q)
        bv = Fraction(beta.p, beta.q)
        with mp.workdps(mp.mp.dps):
            lov = mp.sqrt(av * mp.mpf(X1f.numerator) / X1f.denominator + bv)
            hiv = mp.sqrt(av * mp.mpf(X2f.numerator) / X2f.denominator + bv)
        vlo, vhi = (lov, hiv) if lov < hiv else (hiv, lov)
        orient = 1 if lov < hiv else -1
        # x(v) = (v^2 - beta)/alpha as fmpq_poly in v
        xv = _fqp([-bv / av, Fraction(0), 1 / av])
        for n in range(N + 1):
            m = 2 * n + 1
            # denominator: lead^m * prod_i ((v-rho)(v+rho)/alpha)^(mi*m)
            #              * v^(2n)  [from L^(n+1/2) = v^(2n+1), one v into dx]
            dfacs = []
            alpha_pow = 0
            for rho, mi in vfacs_base:
                dfacs.append((rho, mi * m))
                dfacs.append((-rho, mi * m))
                alpha_pow += mi * m
            if 2 * n > 0:
                dfacs.append((Fraction(0), 2 * n))
            scF = Fraction(2) * av ** alpha_pow / (av * leadF ** m)
            sc = mp.mpf(scF.numerator) / mp.mpf(scF.denominator)
            # g_n(x(v)), x(v)^k composed exactly in fmpq
            gn = g[n]
            gnc = _poly_coeffs_frac(gn, gn.degree() + 1) \
                if gn.degree() >= 0 else [Fraction(0)]
            gnv = flint.fmpq_poly([0])
            pw = flint.fmpq_poly([1])
            for c in gnc:
                gnv = gnv + pw * flint.fmpq_poly([_fq(c)])
                pw = pw * xv
            for k in range(kmax + 1):
                numer = gnv
                for _ in range(k):
                    numer = numer * xv
                v = int_rational(numer, dfacs, vlo, vhi)
                out[k].append(mp.re(v) * sc * orient)
        return out

    raise ValueError(
        "analytic_seed_series: base fibre is neither a perfect square nor "
        "square x linear — use a Frobenius/regular seed or another base")
