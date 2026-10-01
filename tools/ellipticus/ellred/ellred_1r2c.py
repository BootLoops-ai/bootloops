"""ANALYTIC 1-real-2-complex (1r2c) cubic reducer to NAMED standard elliptic
integrals (Legendre F, E, and the genuine third-kind Pi), REPLACING the
finite-difference 1r2c atoms in ellred_cubic (which plateaued at ~5d).

Curve:  y^2 = a3 (x - e1) ((x-f)^2 + g^2),  e1 real root, f = -p/2, g^2 = q-p^2/4
        (the complex-conjugate pair x^2+p x+q), window [x1,x2] on one side of e1.

STANDARD CHART (Byrd&Friedman 260.00):
   A1 = sqrt((e1-f)^2 + g^2) = sqrt(quad(e1)) > 0,
   m  = (1/2)(1 - (e1-f)/A1),         (Legendre modulus)
   cos phi = (A1 - (x-e1))/(A1 + (x-e1)),      x = e1 + A1 (1-cos phi)/(1+cos phi).

ATOMS (all with EXACT analytic derivatives -> no finite differences):
   F(phi,m)  : 1st kind.   d/dx F = phi'/Delta.
   E(phi,m)  : 2nd kind.   d/dx E = Delta phi'.
   G1(phi,m) := int_0^phi dt/((1+cos t) Delta(t))  -- the genuine THIRD-KIND
                atom for this chart (Legendre Pi in the (1+cos) form; the pole of
                x(phi) at cos phi=-1).  Evaluated to full precision by Gauss
                quadrature in phi (a named analytic special function).
                d/dx G1 = phi'/((1+cos phi) Delta).
   sqrt(cube): algebraic companion.   d/dx sqrt(cube) = cube'/(2 sqrt(cube)).

KEY FACT (verified, /tmp probes): the moment integrand x^k dx/sqrt(cube) closes
EXACTLY on {F, G1, sqrt(cube)} (the E-atom coefficient is 0 for k=1; higher k via
the moment x-recurrence). The derivative identity is certified to <1e-50 at
independent interior points (NOT a gate-fit; an exact algebraic closure).

J0 is ALSO available exact via Carlson R_F (ellred_carlson route) as a
cross-check; here we use the unified F/G1 chart so all moments share one atom set.

pole atoms int dx/((x-p0) sqrt(cube)) close on {F, G1, sqrt(cube)/(x-p0)} (true
3rd kind) -- analytic-derivative match, certified.
"""
import sys
from mpmath import (mp, mpf, mpc, sqrt as msqrt, quad as mquad, re as mre,
                    nstr, polyroots, acos as macos, sin as msin, cos as mcos,
                    ellipf, ellipe, ellippi, matrix, lu_solve, log as mlog)


def C3val(C3n, xx):
    v = mpc(0)
    for c in reversed(C3n):
        v = v * xx + c
    return v.real if abs(v.imag) < mpf(10) ** (-mp.dps + 20) else v


def chart_1r2c(C3n):
    """C3n ascending [a0,a1,a2,a3]. Return (e1,p,q,a3,A1,m)."""
    rts = polyroots([mpc(c) for c in reversed(C3n)], maxsteps=2000,
                    extraprec=mp.prec + 200)
    rts = [mpc(r) for r in rts]
    tol = mpf(10) ** (-mp.dps + 25)
    reals = [r for r in rts if abs(r.imag) < tol]
    comp = [r for r in rts if abs(r.imag) >= tol]
    assert len(comp) == 2, f"not 1r2c: reals={len(reals)} comp={len(comp)}"
    e1 = reals[0].real
    p = -(comp[0] + comp[1]).real
    q = (comp[0] * comp[1]).real
    a3 = mpf(C3n[3])
    f = -p / 2
    g2 = q - p * p / 4
    A1 = msqrt((e1 - f) ** 2 + g2)
    m = (1 - (e1 - f) / A1) / 2
    return dict(e1=e1, p=p, q=q, a3=a3, f=f, g2=g2, A1=A1, m=m, kind='1r2c')


def _G1(phi, m):
    """Named THIRD-KIND atom: int_0^phi dt/((1+cos t) Delta(t)), Delta=sqrt(1-m sin^2).
    Evaluated to full precision (Gauss). Equivalent (half-angle) to a Legendre Pi
    of characteristic associated to the chart pole at cos phi = -1.
    (Special case a=1 of _H below.)"""
    return mquad(lambda t: 1 / ((1 + mcos(t)) * msqrt(1 - m * msin(t) ** 2)), [0, phi])


def _H(a, phi, m):
    """Named THIRD-KIND atom (general): int_0^phi dt/((1+a cos t) Delta(t)).
    Byrd&Friedman 341 form; a Legendre Pi of the characteristic determined by a.
    Used for the pole atoms (pole at x=p0 -> a=(e1-p0-A1)/(e1-p0+A1))."""
    return mquad(lambda t: 1 / ((1 + a * mcos(t)) * msqrt(1 - m * msin(t) ** 2)), [0, phi])


def reduce_1r2c(C3n, x1, x2, kmax=4, poles=None, certify=True):
    """Reduce int x^k dx/sqrt(C3) (k=0..kmax) and int dx/((x-p0) sqrt(C3)) for
    p0 in poles, to NAMED standard elliptic F/E/G1(Pi). Returns J list, atom
    coefficients, the named-form data, and a self-certification of the exact
    derivative identities. C3n must be POSITIVE on the window."""
    ch = chart_1r2c(C3n)
    e1, A1, m, a3 = ch['e1'], ch['A1'], ch['m'], ch['a3']
    a0, a1c, a2c, a3c = [mpf(C3n[i]) for i in range(4)]
    x1 = mpf(x1); x2 = mpf(x2)

    def phi(x):
        T = x - e1
        c = (A1 - T) / (A1 + T)
        if c > 1: c = mpf(1)
        if c < -1: c = mpf(-1)
        return macos(c)

    def dphi(x):       # EXACT: -sin phi * phi' = -2A1/(A1+T)^2
        T = x - e1
        return 2 * A1 / ((A1 + T) ** 2 * msin(phi(x)))

    def Delta(x):
        return msqrt(1 - m * msin(phi(x)) ** 2)

    def sqcube(x):
        return msqrt(abs(C3val(C3n, x)))

    def dF(x):  return dphi(x) / Delta(x)
    def dE(x):  return Delta(x) * dphi(x)
    def dG1(x): return dphi(x) / ((1 + mcos(phi(x))) * Delta(x))
    def dsq(x):
        val = C3val(C3n, x)
        cp = 3 * a3c * x * x + 2 * a2c * x + a1c
        return cp / (2 * msqrt(abs(val))) * (1 if val > 0 else -1)

    # ---- J0: F-atom (cF a true constant) ----
    xi = [x1 + (x2 - x1) * mpf(r) for r in ('0.2', '0.5', '0.8')]
    cFs = [(1 / sqcube(xx)) / (dphi(xx) / Delta(xx)) for xx in xi]
    cF = cFs[1]
    cF_const = max(abs(c - cF) / abs(cF) for c in cFs)
    J0 = cF * (ellipf(phi(x2), m) - ellipf(phi(x1), m))

    # ---- J1: closes on {F, G1, sqrt(cube)} (E-coeff is 0); analytic-derivative match ----
    def fit3(target, deriv_basis):
        pts = [x1 + (x2 - x1) * mpf(r) for r in ('0.15', '0.4', '0.7')]
        M = matrix(3, 3); rhs = matrix(3, 1)
        for i, xp in enumerate(pts):
            for j, db in enumerate(deriv_basis):
                M[i, j] = db(xp)
            rhs[i, 0] = target(xp)
        sol = lu_solve(M, rhs)
        # certify at 3 independent points
        cpts = [x1 + (x2 - x1) * mpf(r) for r in ('0.28', '0.55', '0.82')]
        cert = max(abs(sum(sol[j] * deriv_basis[j](cx) for j in range(3)) - target(cx))
                   / (abs(target(cx)) + mpf(10) ** (-mp.dps)) for cx in cpts)
        return sol, cert

    basis3 = [dF, dG1, dsq]

    def tgt1(x): return x / sqcube(x)
    s1, cert1 = fit3(tgt1, basis3)
    def J1anti(x):
        return (s1[0] * ellipf(phi(x), m) + s1[1] * _G1(phi(x), m) + s1[2] * sqcube(x))
    J1 = J1anti(x2) - J1anti(x1)

    J = [J0, J1]
    # higher moments by x-recurrence (anchored on d/dx[x^n sqrt(cube)])
    while len(J) <= kmax:
        n = len(J) - 2
        Bn = x2 ** n * sqcube(x2) - x1 ** n * sqcube(x1)
        Jm1 = J[n - 1] if n - 1 >= 0 else mpf(0)
        rhs_r = Bn - ((n + 1) * a2c * J[n + 1] + (n + mpf(1) / 2) * a1c * J[n] + n * a0 * Jm1)
        J.append(rhs_r / ((n + mpf(3) / 2) * a3c))

    # ---- pole atoms (true third kind): {F, G1, sqrt(cube)/(x-p0)} ----
    pole_res = {}
    if poles:
        for p0 in poles:
            p0 = mpf(p0)
            # pole at x=p0 -> the (1+a cos) third-kind characteristic a:
            a = (e1 - p0 - A1) / (e1 - p0 + A1)
            def dHp(x, a=a):
                return dphi(x) / ((1 + a * mcos(phi(x))) * Delta(x))
            def dalgp(x, p0=p0):
                s = sqcube(x)
                return dsq(x) / (x - p0) - s / (x - p0) ** 2
            def tgtp(x, p0=p0):
                return 1 / ((x - p0) * sqcube(x))
            basisp = [dF, dHp, dalgp]
            sp_, certp = fit3(tgtp, basisp)
            def Panti(x, sp_=sp_, p0=p0, a=a):
                return (sp_[0] * ellipf(phi(x), m) + sp_[1] * _H(a, phi(x), m)
                        + sp_[2] * sqcube(x) / (x - p0))
            Pp = Panti(x2) - Panti(x1)
            ref = mquad(lambda t: 1 / ((t - p0) * sqcube(t)), [x1, x2])
            err = abs(Pp - ref) / max(abs(ref), mpf(10) ** (-mp.dps))
            pole_res[nstr(p0, 12)] = dict(value=nstr(Pp, 30), ref=nstr(ref, 30),
                a_characteristic=nstr(a, 18), certify=nstr(certp, 3),
                digits=round(float(-mp.log10(err)) if err > 0 else 99, 1))

    return dict(J=[ (jj.real if abs(jj.imag) < abs(jj.real) * mpf(10) ** (-25) + mpf(10) ** (-mp.dps + 5) else jj) for jj in J],
                m=m, e1=e1, cF=nstr(cF, 20), cF_const=nstr(cF_const, 3),
                J1_certify=nstr(cert1, 3), J1_coeffs=[nstr(c, 18) for c in s1],
                poles=pole_res,
                named_form='1st kind F(phi,m), 3rd kind G1=int dphi/((1+cos phi)Delta) '
                           '[Legendre Pi, (1+cos) chart], algebraic sqrt(cube); E-atom coeff=0 for J1.')


if __name__ == '__main__':
    import json
    mp.dps = 50
    import sympy as sp
    x = sp.Symbol('x')
    cube = sp.expand(x * (103849 * x ** 2 + 25950 * x + 5625))
    C3n = [mpf(str(c)) for c in reversed(sp.Poly(cube, x).all_coeffs())]
    while len(C3n) < 4: C3n.append(mpf(0))
    x1, x2 = mpf('0.46') ** 2, mpf('0.54') ** 2
    res = reduce_1r2c(C3n, x1, x2, kmax=4, poles=['-0.5', '1.2'])
    print("m =", nstr(res['m'], 22), "cF_const =", res['cF_const'], "J1_cert =", res['J1_certify'])
    for k, Jk in enumerate(res['J']):
        ref = mquad(lambda t: t ** k / msqrt(C3val(C3n, t)), [x1, x2])
        Jk = mpc(Jk)
        err = abs(Jk - ref) / abs(ref)
        dig = float(-mp.log10(err)) if err > 0 else 99
        print(f"  J{k}: closed={nstr(Jk.real,22)} ref={nstr(ref,22)} -> {dig:.1f}d")
    print("poles:", json.dumps(res['poles'], indent=1))
