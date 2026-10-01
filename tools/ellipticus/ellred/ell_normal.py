"""ell_normal — vendored Legendre-chart layer for the ellred member:
P4_exact (demonstration quartic family), quartic_pairs, legendre_chart,
I0_legendre, G_pval.

Note: this module never mutates ambient mp.dps at import time — the
ellipticus package contract is "ambient mp.dps never left modified";
callers set dps explicitly.
"""
import sympy as sp
from mpmath import (mp, mpf, mpc, sqrt as msqrt, atan as matan, ellipf,
                    ellipe, ellippi, quad as mquad, nstr, polyroots,
                    matrix, qr_solve, norm, sin as msin, cos as mcos)

tau, z23, x = sp.symbols('tau z23 x', positive=True)

# ---------------------------------------------------------------- (1) exact P4
def P4_exact():
    """x-core (numerator) of Q1(dL) with z34 = (1-z23) x."""
    Z = 1 - z23
    z34 = Z*x
    z24 = Z*(1 - x)
    c = 1 - 2*z23/((1 - z34)*(1 - z24))
    s2 = 4*z23*z24*z34/((1 - z34)**2*(1 - z24)**2)
    b = (1 - tau)/(1 + tau)
    l1 = (1 - tau)*(1 - z34)
    l0 = (1 + z34) - tau*(1 - z34)
    dL = -l0/l1
    Q1dL = sp.together(dL**2 - 2*b*c*dL + b**2 - s2)
    nn, dd = sp.fraction(sp.cancel(Q1dL))
    nn = sp.expand(nn)
    # strip square factors: keep the odd-multiplicity core
    fl = sp.factor_list(nn, x)
    core = fl[0]
    sqpart = sp.Integer(1)
    for bb, ee in fl[1]:
        if ee % 2:
            core = core*bb**1
        sqpart = sqpart*bb**(ee//2)
    # core currently includes the constant content fl[0]
    if isinstance(core, sp.Integer) or not core.has(x):
        raise RuntimeError('core lost x')
    return sp.expand(core), sqpart, sp.factor(dd)

# ------------------------------------------------- (2) Legendre reduction
def quartic_pairs(P4n):
    """P4 mpc coeffs ascending -> (A, [(p1,q1),(p2,q2)]) with
    P4 = A (x^2+p1x+q1)(x^2+p2x+q2), conjugate pairs combined (real)."""
    # full-precision coeffs (mpc), NOT complex(c) -- the float64 cast was the
    # dps-independent ~3e-17 floor in the I0_legendre reference (roots good only
    # to ~16d -> a1,b1,a2,b2 good to ~16d).  polyroots respects mp.prec on mpc.
    rts = polyroots([mpc(c) for c in reversed(P4n)], maxsteps=600,
                    extraprec=mp.prec + 100)
    rts = [mpc(r) for r in rts]
    used = [False]*4
    quads = []
    for i in range(4):
        if used[i]:
            continue
        for j in range(i + 1, 4):
            if used[j] and False:
                continue
            if not used[j] and abs(rts[j] - rts[i].conjugate()) < \
                    mpf(10)**(-mp.dps + 25)*max(1, abs(rts[i])):
                used[i] = used[j] = True
                p = -(rts[i] + rts[j]).real
                q = (rts[i]*rts[j]).real
                quads.append((mpf(p), mpf(q)))
                break
        else:
            raise RuntimeError('roots not in conjugate pairs')
    A = P4n[-1]
    return A, quads

def legendre_chart(P4n):
    """returns dict with Mobius (t1,t2), a_i,b_i, modulus m, and the
    atom evaluators F,E,Pi as functions of x on a panel."""
    A, ((p1, q1), (p2, q2)) = quartic_pairs(P4n)
    # lambda-resolvent: S1 - lam S2 perfect square:
    # disc: (p1 - lam p2)^2 - 4(1-lam)(q1 - lam q2) = 0
    lam = sp.Symbol('lam')
    dd = sp.expand((p1 - lam*p2)**2 - 4*(1 - lam)*(q1 - lam*q2))
    cs = sp.Poly(dd, lam).all_coeffs()   # desc: a lam^2 + b lam + c
    # the resolvent dd is built from p1,q1,p2,q2 (mpf from quartic_pairs); its
    # coeffs are sympy Floats at mp.dps precision -- mpf(str(.)) keeps them
    # (str of an mp-dps Float is full width).  Use mpmathify to avoid any
    # repr rounding.
    a_, b_, c_ = [mpf(mp.mpf(str(t))) for t in cs]
    disc = b_*b_ - 4*a_*c_
    assert disc > 0, 'lambda resolvent: complex (extend)'
    l1_ = (-b_ + msqrt(disc))/(2*a_)
    l2_ = (-b_ - msqrt(disc))/(2*a_)
    ts = []
    for lv in (l1_, l2_):
        ts.append((p1 - lv*p2)/(2*(lv - 1)))   # S1-lam S2 = (1-lam)(x-t)^2
    t1, t2 = ts
    # S_i = a_i (x-t1)^2 + b_i (x-t2)^2: solve 2x2 linear (exact relations)
    # at x = t2: S_i(t2) = a_i (t2-t1)^2; at x = t1: S_i(t1) = b_i (t1-t2)^2
    d12 = (t1 - t2)**2
    S1 = lambda xx: (xx + p1)*xx + q1
    S2 = lambda xx: (xx + p2)*xx + q2
    a1 = S1(t2)/d12; b1 = S1(t1)/d12
    a2 = S2(t2)/d12; b2 = S2(t1)/d12
    return dict(A=A, t1=t1, t2=t2, a1=a1, b1=b1, a2=a2, b2=b2,
                p1=p1, q1=q1, p2=p2, q2=q2)

def I0_legendre(P4n, x1, x2):
    """int_{x1}^{x2} dx/sqrt(P4) as c*F-difference; returns (value, data)."""
    ch = legendre_chart(P4n)
    A, t1, t2 = ch['A'], ch['t1'], ch['t2']
    a1, b1, a2, b2 = ch['a1'], ch['b1'], ch['a2'], ch['b2']
    # x = (t2 zeta - t1)/(zeta - 1), zeta = (x - t1)/(x - t2)
    # dx/sqrt(P4) = sgn * dzeta /(sqrt(A) sqrt((a1 z^2+b1)(a2 z^2+b2)))
    # (the (x-t2)^2 factors cancel exactly)
    zof = lambda xx: (xx - t1)/(xx - t2)
    # S_i(x) = (x - t2)^2 (a_i zeta^2 + b_i)  (verified A*S1*S2 == P4 and the
    # per-factor (x-t2)^2 identity).  So sqrt(P4) = sqrt(A)*(x-t2)^2*sqrt(prod)
    # and dx = (x-t2)^2/(t1-t2) dzeta  (zeta'(x) = (t1-t2)/(x-t2)^2), hence
    #   dx/sqrt(P4) = dzeta / ( (t1-t2) sqrt(A) sqrt(prod) ).
    # S14 carried only 1/sqrt(A) and DROPPED the 1/(t1-t2) Mobius Jacobian -- that
    # was the [verify] bug (derivative-cert ratio == exactly (t1-t2): 2.27 at A,
    # 17.8 at C, ...).  jac restores it; derivative ratio is now +-1 to dps.
    jac = 1/(t1 - t2)
    # int dz/sqrt((a1 z^2 + b1)(a2 z^2 + b2)): assume a_i z^2 + b_i > 0 forms
    # = 1/sqrt(a1 b2) F(atan(sqrt(a1/b1) z), 1 - (a2 b1)/(a1 b2))
    al2 = a1/b1
    m_ = 1 - (a2*b1)/(b2*a1)
    cF = 1/msqrt(a1*b2)
    def Fat(xx):
        zz = zof(xx)
        return jac*cF*ellipf(matan(msqrt(al2)*zz), m_)
    val = (Fat(x2) - Fat(x1))/msqrt(A)
    # derivative certification + sign fix
    h = mpf(10)**(-mp.dps//3)
    xm = (x1 + x2)/2
    der = (Fat(xm + h) - Fat(xm - h))/(2*h)/msqrt(A)
    tru = 1/msqrt(G_pval(P4n, xm))
    sgn = der/tru
    assert abs(abs(sgn) - 1) < mpf(10)**(-mp.dps//3 + 8), \
        f'I0 chart derivative mismatch {nstr(der,10)} vs {nstr(tru,10)}'
    sgn = mpf(1) if abs(sgn - 1) < 1 else mpf(-1)
    return sgn*val, dict(ch=ch, m=m_, alpha2=al2, cF=cF, jac=jac, sgn=sgn)

def G_pval(p, xx):
    v = mpc(0)
    for c in reversed(p):
        v = v*xx + c
    return v if abs(v.imag) > 0 else v.real
