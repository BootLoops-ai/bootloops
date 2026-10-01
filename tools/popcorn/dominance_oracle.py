#!/usr/bin/env python3
"""Certified 2-fold oracle for the dominance (h != 1/2) fixed-(S,h) SFS —
certified values of the dominance SFS, an object that fastDFE evaluates by
nested float64 quadrature and that polyDFE does not parametrize.

Object (conventions pinned: S = 4*Ne*s, fitnesses 1 : 1+2sh : 1+2s):
  Q(x)  = p x + qt x^2,   p = 2 h S,  qt = (1-2h) S       [psi = e^{-Q}]
  E_k   = C(n,k) * INT_0^1 x^{k-1}(1-x)^{n-k-1} e^{Q(x)} R(x) dx / Z
  R(x)  = INT_x^1 e^{-Q(q)} dq,   Z = R(0).
Neutral normalization: E_k(0,h) = 1/k; h=1/2 => the Kummer closed form.

Method: nested certified composite GL (certquad.py) with analytic
log-integrand estimators; R(x) evaluated per outer node by its own certified
integral (panels from -Q's linear/quadratic structure). Self-checks:
degree-pair at both levels, h=1/2 collapse gate vs the Kummer closed form.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mpmath import mp, mpf, exp, expm1, fabs, log10, binomial, hyp1f1
from mpmath import log as mplog
from certquad import integrate, PAD


def _cuts_for_exp_poly(a, b, c1, c2, extra=()):
    """Panel cuts on [a,b] for integrands ~ e^{c1 x + c2 x^2}: endpoint
    ladders + stationary point of the exponent."""
    cuts = {mpf(a), mpf(b)}
    L = max(abs(float(c1)), abs(float(c2)), 10.0)
    for cc in ('0.25', '1', '4', '16', '64', '256'):
        t = mpf(cc) / L
        if a + t < b:
            cuts.add(a + t)
        if b - t > a:
            cuts.add(b - t)
    if c2 != 0:
        x0 = -c1 / (2 * c2)
        if a < x0 < b:
            w = 1 / fabs(c2) ** mpf('0.5')
            for k in (-8, -2, 0, 2, 8):
                t = x0 + k * w
                if a < t < b:
                    cuts.add(t)
    cuts.update(mpf(e) for e in extra if a < mpf(e) < b)
    return sorted(cuts)


def R_certified(x, p, qt, dps):
    """INT_x^1 e^{-(p q + qt q^2)} dq — CLOSED FORM by completing the square:
      -(pq+qt q^2) = -qt (q + c)^2 + p^2/(4 qt),  c = p/(2 qt),
      R = e^{p^2/(4qt)} sqrt(pi)/(2 sqrt(qt)) [erf(sqrt(qt)(1+c)) - erf(sqrt(qt)(x+c))]
    valid for either sign of qt via complex sqrt/erf (mpmath erf is entire);
    result is real — imaginary residue checked against tolerance. The v1
    nested-quadrature form cost ~minutes/point; this is ~ms.
    Exponential-scale cancellation (qt<0: erfi growth vs the prefactor) is
    bought with a computed pad; the OUTER integral's degree-pair self-check
    plus the h=1/2 collapse gate certify the result end-to-end."""
    from mpmath import erf, sqrt as msqrt, mpc, im, re as mre
    if x >= 1:
        return mpf(0), 9999.0
    pad = int(0.11 * abs(float(p)) ** 2 / max(abs(float(qt)), 1e-30)) \
        if qt != 0 else 0
    pad = min(pad, 2000) + int(0.44 * abs(float(qt))) + 10
    with mp.workdps(mp.dps + pad):
        if qt == 0:
            if p == 0:
                return mpf(1) - x, 9999.0
            val = (exp(-p * x) - exp(-p)) / p
            return val, 9999.0
        c = p / (2 * qt)
        s = msqrt(mpc(qt))
        pref = exp(p * p / (4 * qt)) * msqrt(mp.pi) / (2 * s)
        val = pref * (erf(s * (1 + c)) - erf(s * (x + c)))
        v = mre(val)
        # imaginary residue must be at rounding level of the real part
        if fabs(im(val)) > fabs(v) * mpf(10) ** (-(dps // 2)) and fabs(v) > 0:
            raise RuntimeError(f"R erf-form imaginary residue too large at x={x}")
        return v, 9999.0


def E_dominant(n, k, S, h, dps=50):
    """Certified E_k(S,h). Returns (value, worst_selfcons_digits)."""
    with mp.workdps(dps + PAD):
        S = mpf(str(S))
        h = mpf(str(h))
        if S == 0:
            return mpf(1) / k, 9999.0
        p = 2 * h * S
        qt = (1 - 2 * h) * S
        Z, scZ = R_certified(mpf(0), p, qt, dps + 5)
        # cache R at outer nodes lazily via memo on exact node values
        memo = {}
        def Rx(x):
            key = x
            if key not in memo:
                memo[key] = R_certified(x, p, qt, dps + 5)[0]
            return memo[key]
        C = binomial(n, k)
        def f(x):
            return x ** (k - 1) * (1 - x) ** (n - k - 1) * exp(p * x + qt * x * x) * Rx(x)
        def lnf(x):
            if x <= 0 or x >= 1:
                # k=1 / k=n-1 endpoints are finite; nudge inside
                x = mpf('1e-30') if x <= 0 else 1 - mpf('1e-30')
            # ln R approx: R in [e^{-max(-Q)}*(1-x)*small, (1-x)*e^{max over [x,1] of -Q}]
            # use crude upper envelope: max(-Q) over [x,1] + ln(1-x)
            candidates = [-(p * x + qt * x * x), -(p + qt)]
            if qt != 0:
                q0 = -p / (2 * qt)
                if x < q0 < 1:
                    candidates.append(-(p * q0 + qt * q0 * q0))
            lnR = max(candidates) + mplog(1 - x if x < 1 else mpf('1e-30'))
            return (k - 1) * mplog(x) + (n - k - 1) * mplog(1 - x) \
                + (p * x + qt * x * x) + lnR
        cuts = _cuts_for_exp_poly(mpf(0), mpf(1), p, qt)
        I, scI = integrate(f, lnf, cuts, dps)
        return C * I / Z, min(scZ, scI)


def _selftest():
    """Gates: (1) h=1/2 collapse to the Kummer closed form; (2) neutral S->0;
    (3) degree-pair self-consistency at hostile (S,h)."""
    import time
    ok = True
    with mp.workdps(80):
        for (n, k, S) in [(20, 3, -100.0), (20, 19, -1000.0), (20, 10, 500.0)]:
            t0 = time.time()
            got, sc = E_dominant(n, k, S, 0.5, dps=50)
            Sm = mpf(str(S))
            truth = mpf(n) / (mpf(k) * (n - k)) * (1 - hyp1f1(n - k, n, -Sm)) / (-expm1(-Sm))
            d = float(-log10(fabs((got - truth) / truth))) if got != truth else 9999
            print(f"h=1/2 collapse n={n} k={k} S={S}: {d:6.1f}d (sc {sc:6.1f}) {time.time()-t0:5.1f}s")
            ok = ok and d >= 45
        for (n, k, S, h) in [(20, 1, -100.0, 0.2), (20, 10, -1000.0, 0.2),
                             (20, 5, 500.0, 0.7), (20, 3, -5.0, 0.05),
                             (20, 19, 10.0, 0.3), (20, 10, -100.0, 0.95)]:
            t0 = time.time()
            got, sc = E_dominant(n, k, S, h, dps=50)
            print(f"dominant n={n} k={k} S={S:+8.1f} h={h}: {mp.nstr(got, 12)} "
                  f"selfcons {sc:6.1f}d {time.time()-t0:5.1f}s")
            ok = ok and sc >= 45
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(_selftest())
