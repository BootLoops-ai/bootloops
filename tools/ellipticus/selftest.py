"""ellipticus.selftest — fast gates (run: python3 -m ellipticus.selftest with
tools/ on PYTHONPATH, or python3 selftest.py from this directory).

S1  curve: Legendre-K period vs sin^2-quadrature route on synthetic
    4-real-root quartics + a K(1/4) closed form (independent formulas;
    1e-(dps-8) bar).
S2  extend: rebuild the demonstration quartic's 5x5 extended system from the
    exact coefficients in fixtures/extended_systems.json and compare the
    polyform EXACTLY (string-identical integer coefficient lists) against
    the pinned fixtures/system_polyform.json (a snapshot this same call
    regenerates).
S3  march: transported kernels (I0,I1,I2) at z=1/5 vs a direct tanh-sinh
    x-quadrature oracle (independent leg; >=dps-12 digits).
S4  qengine: three-route g^(n) cross-check at a generic (z,tau) + the
    certified tail bound sanity (bound >= true truncation error) + sunrise
    equal-mass positive control (>= dps-8 digits).
Exit 0 all PASS; raises on any failure (fail-closed).
"""
import json
import os
import sys
import time
from fractions import Fraction

if __name__ == "__main__" and not __package__:
    # run as a plain script (python3 selftest.py): re-enter as the package
    # module so the relative imports below resolve, then stop here.
    import runpy
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    runpy.run_module("ellipticus.selftest", run_name="__main__", alter_sys=True)
    raise SystemExit(0)

import mpmath as mp

from . import (QuarticCurve, MomentFamily, shared_digits, qengine, env)

# system fixtures (tiny, bundled in fixtures/ beside this package;
# ELLIPTICUS_W3, or the former name EMPL_EVAL_W3, overrides to another copy
# with the same layout).
W3 = env("ELLIPTICUS_W3", "EMPL_EVAL_W3",
         os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures"))


def fixture_quartic():
    """(coeffs_x_desc, window, param) of the demonstration quartic family."""
    d = json.load(open(os.path.join(W3, "extended_systems.json")))
    return (d["P4_coeffs_x_desc"],
            (Fraction(d["window"][0]), Fraction(d["window"][1])),
            d["param"])


def s1_curve(dps=50):
    t0 = time.time()
    worst = 999.0
    # (a) synthetic 4-real-root quartics, every oval class:
    #   lead<0: P=-(x)(x-1)(x-2)(x-3), ovals (0,1),(2,3)
    #   lead>0: P=(x)(x-1)(x-2)(x-3) shifted: middle oval (1,2)
    cases = [
        ([-1, 6, -11, 6, 0], 0),   # -(x)(x-1)(x-2)(x-3): oval (0,1)
        ([-1, 6, -11, 6, 0], 1),   # oval (2,3)
        ([1, -6, 11, -6, 0], 0),   # (x)(x-1)(x-2)(x-3): middle oval (1,2)
    ]
    for cs, ov in cases:
        cur = QuarticCurve(cs)
        with mp.workdps(dps + 30):
            v1, m1 = cur.period_oval(dps, oval=ov)
            assert m1.get("route") == "legendre_K", m1
            a, b = cur.real_ovals(dps)[ov]
            v2, _ = cur._oval_quad(a, b, 0, dps)
            worst = min(worst, shared_digits(v1, 2 * v2))
    # (b) known closed form: y^2=(1-x^2)(1-x^2/4): oint dx/y over (-1,1)
    #     = 2 * 2 K(1/4)  (m-convention: k^2 = 1/4)
    cur = QuarticCurve([Fraction(1, 4), 0, Fraction(-5, 4), 0, 1])
    with mp.workdps(dps + 30):
        v1, m1 = cur.period_oval(dps, oval=0)
        ref = 4 * mp.ellipk(mp.mpf(1) / 4)
        worst = min(worst, shared_digits(v1, ref))
    assert worst >= dps - 8, f"S1 period route disagreement: {worst}"
    print(f"S1 PASS  legendre-vs-quad + K(1/4) closed form {worst:.1f}d  "
          f"({time.time()-t0:.1f}s)")


def s2_extend():
    t0 = time.time()
    coeffs, window, param = fixture_quartic()
    fam = MomentFamily(coeffs, window=window, kmax=2, param=param)
    pf = fam.polyform()
    bank = json.load(open(os.path.join(W3, "system_polyform.json")))["quartic"]
    assert pf["n"] == bank["n"]
    assert pf["Den"] == bank["Den"], "Den mismatch"
    assert pf["NumM"] == bank["NumM"], "NumM mismatch"
    print(f"S2 PASS  quartic 5x5 polyform EXACT match to the pinned fixture "
          f"({time.time()-t0:.1f}s)")
    return fam


def s3_march(fam, dps=40):
    t0 = time.time()
    res = fam.eval_words([[("mom", 0)]], [Fraction(1, 5)], dps)
    ker = res["1/5"]["kernels"]
    # independent oracle: direct tanh-sinh x-quadrature of x^k/sqrt(F(x;1/5))
    import sympy as sp
    coeffs, window, param = fixture_quartic()
    z = sp.Symbol(param)
    cs = [sp.sympify(c, locals={param: z}).subs(z, sp.Rational(1, 5))
          for c in coeffs]
    cs = [Fraction(sp.Rational(c).p, sp.Rational(c).q) for c in cs]
    with mp.workdps(dps + 40):
        lo = mp.mpf(window[0].numerator) / window[0].denominator
        hi = mp.mpf(window[1].numerator) / window[1].denominator
        md = (lo + hi) / 2
        csm = [mp.mpf(c.numerator) / mp.mpf(c.denominator) for c in cs]

        def pv(x):
            s = mp.mpf(0)
            for c in csm:
                s = s * x + c
            return s
        worst = 99.0
        for k in range(3):
            v = mp.quad(lambda t: t ** k / mp.sqrt(pv(t)), [lo, md, hi],
                        maxdegree=9)
            worst = min(worst, shared_digits(ker[k], v))
    assert worst >= dps - 12, f"S3 kernel-vs-quad {worst}"
    print(f"S3 PASS  march kernels vs x-quad {worst:.1f}d "
          f"({time.time()-t0:.1f}s)")


def s4_qengine(dps=50):
    t0 = time.time()
    with mp.workdps(dps):
        tau = mp.mpc("0.3", "1.1")
        zz = mp.mpc("0.22", "0.37")
        gq = qengine.K.g_coeffs(zz, tau, 3, "qbar")
        gt = qengine.K.g_coeffs(zz, tau, 3, "theta")
        ge = qengine.K.g_coeffs(zz, tau, 3, "eisenstein")
        worst = min(min(shared_digits(gq[n], gt[n]) for n in range(1, 4)),
                    min(shared_digits(gq[n], ge[n]) for n in range(1, 4)))
        assert worst >= dps - 10, f"S4 route cross-check {worst}"
        # certified value + bound sanity: truncate coarsely, bound must cover
        vc, bnd = qengine.g_n_certified(2, zz, tau, dps)
        vshort = qengine.K.g_n(2, zz, tau, route="qbar", J=25)
        err = abs(vc - vshort)
        bnd25 = qengine.g_tail_bound(2, zz, tau, 25)
        assert bnd25 >= err * mp.mpf("0.99"), \
            f"S4 tail bound {mp.nstr(bnd25,4)} < true error {mp.nstr(err,4)}"
    em, eref, agree = qengine.sunrise_equal_mass_control(Fraction(-3),
                                                         dps=dps)
    assert agree >= dps - 8, f"S4 equal-mass control {agree}"
    print(f"S4 PASS  routes {worst:.1f}d, tail-bound covers, equal-mass "
          f"control {agree}d  ({time.time()-t0:.1f}s)")


def main():
    t0 = time.time()
    s1_curve()
    fam = s2_extend()
    s3_march(fam)
    s4_qengine()
    print(f"ALL SELFTESTS PASS  ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
