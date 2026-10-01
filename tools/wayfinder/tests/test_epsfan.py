#!/usr/bin/env python3
"""
Tests for tools/wayfinder/epsfan.py.

E1  known Laurent fan     — f = pi^2/e^2 - gammaE/e + zeta(3) + sqrt2*e
                            + log2*e^2 + Catalan*e^3 + zeta(5)*e^4;
                            recover kmin=-2..kmax=4 at dps=100 to >=90 d,
                            diag reports honest (>=90 d) agreement.
E2  deliberately aliased  — same f + BIG*e^10 with nodes=8 (window -2..4):
                            10 = 2 (mod 8) corrupts c_2 through R^8 alias
                            leverage.  The (r,delta) diag MUST catch it:
                            agree_digits[2] collapses while clean orders
                            stay honest, and the reported floor matches the
                            REAL error vs the known c_2.
E3  vector f + half_nodes — vector fan returns per-component lists; the
                            check='half_nodes' second pass works; zero
                            coefficients report dps-level agreement (no
                            false alarms on 0/0).
E4  validation            — nodes < window, bad R, bad check mode.
E5  vandermonde realgrid  — second validated mode (the production
                            Laurent-gate node
                            design: geometric decade grid, eps_max ~1e-3,
                            interleaved verify spares): recovers the known
                            fan >= 90d, self_digits honest; explicit 1/q
                            prime nodes (the production grid) also work;
                            validation (n_verify=0, span<=1, complex node).
E6  anchor round-trip     — boundary cancels: clean coefficients round-trip
                            at ~dps; a 1e-20-relative corruption of c_{-2}
                            is CAUGHT at ~20d (must-fail control);
                            precomputed-value (non-callable) path works.
E7  schwarz_real mirror   — opt-in conjugate mirror: full digits on the
                            real fan with 17-of-32 evaluations per pass;
                            on a COMPLEX-Laurent fan the mirror corrupts
                            the complex order to <5 true digits while the
                            per-k diag stays deceptively clean there (the
                            documented blind spot) — and the full-circle
                            default extracts it correctly.

All references computed from mpmath constants in-test (no fabricated
digits).  Global mp.dps asserted untouched.
"""
import os
import sys

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
sys.path.insert(0, PKG)

import epsfan  # noqa: E402

FAILS = []


def check(name, cond, msg=""):
    tag = "PASS" if cond else "FAIL"
    print("  [%s] %s%s" % (tag, name, (": " + msg) if msg else ""))
    if not cond:
        FAILS.append(name)


def agree_d(a, b, wp):
    with mp.workdps(wp):
        a, b = mp.mpc(a), mp.mpc(b)
        diff = abs(a - b)
        if diff == 0:
            return float(wp)
        s = max(abs(a), abs(b))
        return float(-mp.log10(diff / s)) if s > 0 else float(wp)


DPS0 = mp.mp.dps


def make_f():
    """Laurent fan with irrational mpmath coefficients, evaluated at the
    caller's ambient precision (cauchy_laurent calls f inside workdps)."""
    def coeffs():
        return {-2: mp.pi ** 2,
                -1: -mp.euler,
                0: mp.zeta(3),
                1: mp.sqrt(2),
                2: mp.log(2),
                3: mp.catalan,
                4: mp.zeta(5)}

    def f(e):
        c = coeffs()
        return (c[-2] / e ** 2 + c[-1] / e + c[0] + c[1] * e
                + c[2] * e ** 2 + c[3] * e ** 3 + c[4] * e ** 4)
    return f, coeffs


# ---------------------------------------------------------------------------
def test_known_fan():
    dps = 100
    f, coeffs = make_f()
    out = epsfan.cauchy_laurent(f, -2, 4, dps)  # defaults R='1e-3', nodes=32
    with mp.workdps(dps + 15):
        ref = coeffs()
    worst = 1e9
    for k in range(-2, 5):
        ad = agree_d(out["coeffs"][k], ref[k], dps + 10)
        worst = min(worst, ad)
    check("E1-recovered", worst >= 90, "worst coeff agree=%.1f d" % worst)
    check("E1-diag-honest", out["diag"]["min_agree_digits"] >= 90,
          "min_agree=%.1f" % out["diag"]["min_agree_digits"])
    # imaginary parts must be at the noise floor (real fan)
    with mp.workdps(dps + 10):
        imax = max(abs(out["coeffs"][k].imag) for k in range(-2, 5))
        rmax = max(abs(out["coeffs"][k]) for k in range(-2, 5))
        im_ok = imax < rmax * mp.mpf(10) ** (-90)
    check("E1-real", bool(im_ok))


def test_aliased_fan():
    dps = 100
    f, coeffs = make_f()
    BIG = 10 ** 6

    def g(e):
        return f(e) + BIG * mp.pi * e ** 10

    # nodes=8, window -2..4 (width 7 <= 8): e^10 aliases onto k=2 (10-8)
    # with weight R^8 = 1e-24 -> corrupts c_2 at the ~1e-18 ABS level,
    # i.e. ~17 correct digits instead of 100.  The diag must measure this.
    out = epsfan.cauchy_laurent(g, -2, 4, dps, R='1e-3', nodes=8)
    ag = out["diag"]["agree_digits"]
    with mp.workdps(dps + 15):
        ref = coeffs()
    err2 = agree_d(out["coeffs"][2], ref[2], dps + 10)
    check("E2-c2-really-corrupted", err2 < 30,
          "c_2 true agreement %.1f d (alias is real)" % err2)
    check("E2-diag-catches", ag[2] < 30,
          "diag agree_digits[2]=%.1f" % ag[2])
    check("E2-diag-localizes", ag[-2] >= 90 and ag[-1] >= 90 and ag[0] >= 90,
          "clean orders: %s" % {k: round(ag[k], 1) for k in (-2, -1, 0)})
    check("E2-min-floor", out["diag"]["min_agree_digits"] < 30)
    # the diag floor must not be wildly optimistic vs the real error
    check("E2-floor-honest", ag[2] <= err2 + 8,
          "diag=%.1f vs real=%.1f" % (ag[2], err2))


def test_vector_and_half_nodes():
    dps = 60

    def fv(e):
        return [1 / e + mp.pi, e * mp.euler]

    out = epsfan.cauchy_laurent(fv, -1, 1, dps, R='1e-2', nodes=16,
                                check='half_nodes')
    c = out["coeffs"]
    with mp.workdps(dps + 10):
        ref = {(-1, 0): mp.mpf(1), (0, 0): mp.pi, (1, 0): mp.mpf(0),
               (-1, 1): mp.mpf(0), (0, 1): mp.mpf(0), (1, 1): mp.euler}
        ok = True
        gscale = max(abs(v) for k in c for v in c[k])
        for k in (-1, 0, 1):
            for comp in (0, 1):
                r = ref[(k, comp)]
                if r == 0:
                    ok = ok and abs(c[k][comp]) < gscale * mp.mpf(10) ** (-50)
                else:
                    ok = ok and agree_d(c[k][comp], r, dps + 5) >= 50
    check("E3-vector-values", bool(ok))
    check("E3-diag-no-false-alarm", out["diag"]["min_agree_digits"] >= 50,
          "min_agree=%.1f (zero coeffs must not alarm)"
          % out["diag"]["min_agree_digits"])
    check("E3-nodes2", out["diag"]["nodes2"] == 8 and
          out["diag"]["check"] == "half_nodes")


def test_validation():
    def expect(name, fn):
        try:
            fn()
            check(name, False, "no exception")
        except ValueError:
            check(name, True)

    f = lambda e: 1 / e  # noqa: E731
    expect("E4-window-too-wide",
           lambda: epsfan.cauchy_laurent(f, -3, 5, 50, nodes=8))
    expect("E4-bad-R",
           lambda: epsfan.cauchy_laurent(f, -1, 1, 50, R='0'))
    expect("E4-bad-check",
           lambda: epsfan.cauchy_laurent(f, -1, 1, 50, check='thrice'))
    expect("E4-half-nodes-too-few",
           lambda: epsfan.cauchy_laurent(f, -2, 4, 50, nodes=8,
                                         check='half_nodes'))


def test_vandermonde_realgrid():
    dps = 100
    f, coeffs = make_f()
    out = epsfan.vandermonde_laurent(f, -2, 4, dps)  # default decade grid
    with mp.workdps(dps + 15):
        ref = coeffs()
    worst = min(agree_d(out["coeffs"][k], ref[k], dps + 10)
                for k in range(-2, 5))
    check("E5-recovered", worst >= 90, "worst coeff agree=%.1f d" % worst)
    check("E5-self-honest", out["diag"]["self_digits"] >= 90,
          "self_digits=%.1f" % out["diag"]["self_digits"])
    check("E5-diag-shape", out["diag"]["mode"] == "vandermonde_realgrid"
          and out["diag"]["n_fit"] == 7 and out["diag"]["n_verify"] == 7)

    # explicit 1/q prime nodes — the production grid (first 11 primes)
    qs = [1013, 1117, 1237, 1367, 1511, 1663, 1847, 2039, 2243, 2503, 2741]
    out2 = epsfan.vandermonde_laurent(f, -2, 4, dps,
                                      eps_nodes=["1/%d" % q for q in qs])
    worst2 = min(agree_d(out2["coeffs"][k], ref[k], dps + 10)
                 for k in range(-2, 5))
    check("E5-explicit-nodes", worst2 >= 90, "worst=%.1f d" % worst2)
    check("E5-explicit-verify", out2["diag"]["n_verify"] == 4
          and out2["diag"]["n_nodes"] == 11)

    # vector fan
    def fv(e):
        return [1 / e + mp.pi, e * mp.euler]
    outv = epsfan.vandermonde_laurent(fv, -1, 1, 60, n_verify=3)
    with mp.workdps(70):
        okv = (agree_d(outv["coeffs"][-1][0], mp.mpf(1), 65) >= 50
               and agree_d(outv["coeffs"][0][0], mp.pi, 65) >= 50
               and agree_d(outv["coeffs"][1][1], mp.euler, 65) >= 50)
    check("E5-vector", bool(okv))

    def expect(name, fn):
        try:
            fn()
            check(name, False, "no exception")
        except ValueError:
            check(name, True)
    expect("E5-no-verify-refused",
           lambda: epsfan.vandermonde_laurent(f, -2, 4, 50, n_verify=0))
    expect("E5-bad-span",
           lambda: epsfan.vandermonde_laurent(f, -2, 4, 50, span='1'))
    expect("E5-bad-window",
           lambda: epsfan.vandermonde_laurent(f, 4, -2, 50))
    expect("E5-too-few-explicit",
           lambda: epsfan.vandermonde_laurent(
               f, -2, 4, 50, eps_nodes=['1e-3'] * 7))
    expect("E5-complex-node",
           lambda: epsfan.vandermonde_laurent(
               f, -1, 1, 50, eps_nodes=['1e-3', '2e-3', '3e-3', 1e-3j]))


def test_anchor_roundtrip():
    dps = 100
    f, coeffs = make_f()
    out = epsfan.cauchy_laurent(f, -2, 4, dps)
    # clean coefficients: window-limited fan -> truncation zero, boundary
    # cancels, residual = pure extraction error ~ dps
    rt = epsfan.anchor_roundtrip(f, out["coeffs"], dps)
    check("E6-clean", rt["agree_digits"] >= 90,
          "roundtrip %.1f d" % rt["agree_digits"])
    # must-FAIL control: corrupt c_{-2} by 1e-20 relative — the probe must
    # measure ~20 digits, not stay silent
    with mp.workdps(dps + 10):
        bad = dict(out["coeffs"])
        bad[-2] = bad[-2] * (1 + mp.mpf(10) ** (-20))
    rt2 = epsfan.anchor_roundtrip(f, bad, dps)
    check("E6-corrupt-caught", 10 <= rt2["agree_digits"] < 30,
          "corrupted roundtrip %.1f d (expect ~20)" % rt2["agree_digits"])
    # precomputed-value path (fixed-eps transported value, f not callable)
    with mp.workdps(dps + 20):
        direct = mp.nstr(f(mp.mpmathify('6.1e-4')), dps + 15)
    rt3 = epsfan.anchor_roundtrip(direct, out["coeffs"], dps,
                                  eps_anchor='6.1e-4')
    check("E6-precomputed", rt3["agree_digits"] >= 90,
          "precomputed-value roundtrip %.1f d" % rt3["agree_digits"])
    # vector coeffs path
    def fv(e):
        return [1 / e + mp.pi, e * mp.euler]
    outv = epsfan.cauchy_laurent(fv, -1, 1, 60, R='1e-2', nodes=16)
    rtv = epsfan.anchor_roundtrip(fv, outv["coeffs"], 60, eps_anchor='3.3e-3')
    check("E6-vector", rtv["agree_digits"] >= 50
          and len(rtv["per_component"]) == 2,
          "vector roundtrip %.1f d" % rtv["agree_digits"])


def test_schwarz_mirror():
    dps = 100
    f, coeffs = make_f()
    with mp.workdps(dps + 15):
        ref = coeffs()

    calls = [0]

    def f_count(e):
        calls[0] += 1
        return f(e)

    # real (Schwarz) fan: mirror gives full digits at 17-of-32 evals/pass
    out_m = epsfan.cauchy_laurent(f_count, -2, 4, dps, schwarz_real=True)
    worst = min(agree_d(out_m["coeffs"][k], ref[k], dps + 10)
                for k in range(-2, 5))
    check("E7-real-correct", worst >= 90, "worst=%.1f d" % worst)
    check("E7-evals-halved", calls[0] == 34
          and out_m["diag"]["n_evals"] == 34,
          "n_evals=%d (expect 2 passes x 17)" % calls[0])
    check("E7-flag", out_m["diag"]["schwarz_real"] is True)

    out_full = epsfan.cauchy_laurent(f, -2, 4, dps)
    check("E7-default-full-circle",
          out_full["diag"]["schwarz_real"] is False
          and out_full["diag"]["n_evals"] == 64)

    # COMPLEX-Laurent fan (c_1 = sqrt2 + i*Catalan): the mirror silently
    # destroys the imaginary content — and the per-k diag at the corrupted
    # order is DECEPTIVELY clean (R-independent corruption).  This is the
    # documented blind spot behind the opt-in warning.
    def g(e):
        return f(e) + mp.mpc(0, 1) * mp.catalan * e

    with mp.workdps(dps + 15):
        c1_true = mp.sqrt(2) + mp.mpc(0, 1) * mp.catalan
    out_bad = epsfan.cauchy_laurent(g, -2, 4, dps, schwarz_real=True)
    err1 = agree_d(out_bad["coeffs"][1], c1_true, dps + 10)
    check("E7-complex-mirror-corrupts", err1 < 5,
          "mirror c_1 true agreement %.2f d (catastrophic, not degraded)"
          % err1)
    check("E7-blind-spot-documented",
          out_bad["diag"]["agree_digits"][1] >= 30,
          "diag[1]=%.1f looks clean while c_1 is WRONG — why the mirror "
          "is opt-in" % out_bad["diag"]["agree_digits"][1])
    # full-circle default extracts the complex row correctly
    out_good = epsfan.cauchy_laurent(g, -2, 4, dps)
    err1g = agree_d(out_good["coeffs"][1], c1_true, dps + 10)
    check("E7-full-circle-correct", err1g >= 90, "c_1 %.1f d" % err1g)


def test_zz_all_checks_enforced():
    """pytest enforcement (runs last): the check() accumulator must be
    empty — without this, check() failures in the functions above would be
    invisible to pytest (they only set the __main__ exit code)."""
    assert not FAILS, "epsfan self-checks failed: %s" % ",".join(FAILS)
    assert mp.mp.dps == DPS0, "global mp.dps mutated (import-dps footgun)"


if __name__ == "__main__":
    print("== test_epsfan ==")
    test_known_fan()
    test_aliased_fan()
    test_vector_and_half_nodes()
    test_validation()
    test_vandermonde_realgrid()
    test_anchor_roundtrip()
    test_schwarz_mirror()
    check("global-dps-untouched", mp.mp.dps == DPS0,
          "dps %d -> %d" % (DPS0, mp.mp.dps))
    print("== %s ==" % ("ALL PASS" if not FAILS else "FAILURES: " + ",".join(FAILS)))
    sys.exit(0 if not FAILS else 1)
