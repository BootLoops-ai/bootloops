# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision.
r"""
bench_box.py -- benchmark the fixed Gauss-product engine on the one-loop
massive box.

After the closed-form loop-parameter integration the box is a smooth 2-D
Feynman integral; the fixed Gauss-Legendre product rule (hiprec.
integrate_gauss_product) then converges spectrally to 40+ digits.  Two
independent reductions are shipped: `box` (Cheng-Wu x3=1, loop parameter x0
integrated out) and `box_crosschart` (Cheng-Wu x0=1, mass parameter x2
integrated out).  They share no algebra beyond U and F, so their agreement
validates both the reduction and the quadrature; the pinned references below
are cross-validated this way to 46+ digits.

The box (matched to the pySecDec normalization, eps^0):
    Box(s,t,m2,M5) = \int_simplex 1/F^2 ,
    F = (M5 x0 + m2 x1 + m2 x2 + m2 x3) U - s x0 x2 - u x1 x3 ,  u=-s-t .

PRECISION CONTRACT: pass kinematics as EXACT inputs (int, Fraction, decimal
string, or a binary-exact float like -1.0).  An mpf built at a lower ambient
precision (e.g. mpf(-1)/3 at the mpmath default dps) is a *different rational
number* at the 1e-16 level; every engine then converges spectrally to the
integral of THAT number, silently capping agreement with the true value near
the input's own precision.  All conversions here happen at working precision.
"""
import sys, time
from fractions import Fraction
import mpmath as mp


def _I2(P, Q, R):
    r"""\int_0^inf dx/(P x^2 + Q x + R)^2 , P>0, positive quadratic on [0,inf).

    Closed form on both discriminant branches; near D = 4PR - Q^2 = 0 the two
    closed-form terms cancel catastrophically (each ~1/D), so for
    |D| <= Q^2/8 the kernel switches to the series in e = D/(4P^2) around the
    double-root point, which is uniformly accurate there and converges
    geometrically at ratio <= 1/8 (~0.9 digits/term):

        I2 = (1/P^2) sum_{k>=0} (-1)^k (k+1)/(2k+3) e^k / c^{2k+3},
        c = Q/(2P),  radius |e| < c^2 (used for |e| <= c^2/8).
    """
    D = 4*P*R - Q*Q
    c = Q/(2*P)
    c2 = c*c
    e = D/(4*P*P)
    if abs(e) <= c2/8:
        # series branch: geometric convergence, no 1/D cancellation
        eps = mp.mpf(10)**(-(mp.mp.dps + 5))
        acc = mp.mpf(0)
        w = c**-3                      # (-e)^k / c^{2k+3} at k=0
        k = 0
        while True:
            term = (mp.mpf(k+1)/(2*k+3)) * w
            acc += term
            if abs(term) < eps*abs(acc):
                break
            w *= -e/c2
            k += 1
            if k > 100*mp.mp.dps:
                raise ArithmeticError("_I2 series failed to converge")
        return acc/(P*P)
    if D > 0:
        sD = mp.sqrt(D)
        J1 = (2/sD)*(mp.pi/2 - mp.atan(Q/sD))
    else:
        sN = mp.sqrt(-D)
        J1 = -(1/sN)*mp.log((Q - sN)/(Q + sN))
    return (-Q/(D*R)) + (2*P/D)*J1


def _as_mpf(v):
    """Convert an input to mpf at the CURRENT working precision.
    Fraction, int and decimal strings convert exactly; float and mpf are taken
    at face value (their binary value is used bit-exactly -- a value that was
    already rounded, like mpf(-1)/3 at low ambient dps, stays rounded)."""
    if isinstance(v, Fraction):
        return mp.mpf(v.numerator)/v.denominator
    return mp.mpf(v)


def box(s, t, m2, M5, dps=40):
    """eps^0 box via 1 analytic loop integration + 2-D FIXED Gauss-Legendre product.

    Deep-Euclidean, u = -s-t < 4 m2 (smooth integrand, no threshold crossing).
    Pass EXACT kinematics (int, Fraction, decimal string) -- see the module
    docstring's precision contract.

    The deterministic Gauss product rule (hiprec.integrate_gauss_product) is used
    rather than adaptive nested tanh-sinh: the latter has a ~17-digit nested
    error floor (the inner quad's adaptive tolerance injects noise the outer quad
    cannot integrate through), whereas the fixed product rule converges
    SPECTRALLY and self-stabilises to the full requested dps -- agreeing with the
    independent cross-chart reduction `box_crosschart` to 46+ d (measured at
    dps=45/60 on the M5=1,2,5 references).
    """
    import os, sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from hiprec import integrate_gauss_product
    work = dps + 25
    mp.mp.dps = work
    s = _as_mpf(s); t = _as_mpf(t); m2 = _as_mpf(m2); M5 = _as_mpf(M5)
    u = -s - t
    if not (s < 0 and t < 0 and m2 > 0 and M5 > 0):
        raise ValueError("box: require deep-Euclidean s<0,t<0,m2>0,M5>0")
    if u >= 4*m2:
        raise ValueError("box: u>=4 m2 pseudo-threshold regime not supported")

    # Integrate over (x1, X) with X = x1 + x2 + 1 >= x1 + 1.  Matching the inner
    # lower limit to x1+1 and the tail scale to 1+x1 is what makes the fixed
    # Gauss rule spectral (40+ d); a flat [0,inf)/L=2 box plateaus at ~5 d.
    one = mp.mpf(1)
    a = M5 + m2 - s

    def integrand(x1, X):
        return _I2(M5, a*X + s*(x1+1), m2*X*X - u*x1)

    def limits(fixed):
        if not fixed:               # outer x1 axis
            return (mp.mpf(0), mp.mpf(2))
        x1 = fixed[0]               # inner X axis: lives on [x1+1, inf)
        return (x1 + one, one + x1)
    return integrate_gauss_product(integrand, 2, dps=dps, limits=limits)


def box_crosschart(s, t, m2, M5, dps=40):
    """Same box integral through an INDEPENDENT reduction: Cheng-Wu chart x0=1
    (instead of x3=1) and closed-form integration of the mass parameter x2
    (kernel P=m2) instead of the loop parameter x0 (kernel P=M5).  The
    surviving 2-D integral runs over (x1, X) with X = 1 + x1 + x3, using the
    same matched-limits trick as `box`.

    Beyond U and F the two routes share no algebra, so digit-for-digit
    agreement of `box` and `box_crosschart` validates reduction + quadrature
    together (measured: 46-47 d agreement at dps=45 on the pinned references).
    """
    import os, sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from hiprec import integrate_gauss_product
    work = dps + 25
    mp.mp.dps = work
    s = _as_mpf(s); t = _as_mpf(t); m2 = _as_mpf(m2); M5 = _as_mpf(M5)
    u = -s - t
    if not (s < 0 and t < 0 and m2 > 0 and M5 > 0):
        raise ValueError("box_crosschart: require deep-Euclidean s<0,t<0,m2>0,M5>0")
    if u >= 4*m2:
        raise ValueError("box_crosschart: u>=4 m2 pseudo-threshold regime not supported")
    one = mp.mpf(1)

    # F(x0=1) as a quadratic in x2: P = m2, Q = m2 X + A - s, R = A X - u x1 x3,
    # with A = M5 + m2 (x1 + x3) and X = 1 + x1 + x3 (so x3 = X - 1 - x1).
    def integrand(x1, X):
        A = M5 + m2*(X - one)
        return _I2(m2, m2*X + A - s, A*X - u*x1*(X - one - x1))

    def limits(fixed):
        if not fixed:               # outer x1 axis
            return (mp.mpf(0), mp.mpf(2))
        x1 = fixed[0]               # inner X axis: lives on [x1+1, inf)
        return (x1 + one, one + x1)
    return integrate_gauss_product(integrand, 2, dps=dps, limits=limits)


if __name__ == '__main__':
    # Pinned references at s=-1, t=-1/3, m2=1: 40 significant digits each,
    # cross-validated to 46+ d by (a) dps=45 vs dps=60 runs of `box` and
    # (b) the independent `box_crosschart` reduction at dps=45.  The digits
    # are quotable in full; the quadrature is spectral, so requested dps
    # beyond 40 reproduces and extends them (two-dps check recommended
    # before quoting new depth).
    REF = {
        '1.0': '0.1780502267932306438845714594088460679896',
        '2.0': '0.1245570572327092697104024441545542480922',
        '5.0': '0.06920789507050544783290083200087922071809',
    }
    s, t, m2 = -1, Fraction(-1, 3), 1
    dps = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    print(f"=== one-loop box benchmark (hiprec Gauss product, dps={dps}) ===", flush=True)
    for k in ['1.0', '2.0', '5.0']:
        t0 = time.time()
        v = box(s, t, m2, int(float(k)), dps=dps)
        dt = time.time()-t0
        mp.mp.dps = dps + 25
        ref = mp.mpf(REF[k])
        nref = len(REF[k].split('.')[1])
        agree = (-int(mp.log10(abs(v-ref))) if v != ref else 99)
        print(f"M5={k:>4}  I={mp.nstr(v, dps)}", flush=True)
        print(f"        ref={REF[k]}  (agree>= {min(agree, nref)} d, ref has {nref} d)  [{dt:.1f}s]", flush=True)
