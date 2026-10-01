#!/usr/bin/env python3
"""
Tests for the E4 certified layer of tools/wayfinder/epsfan.py
.

Synthetic Laurent function with KNOWN EXACT coefficients (Fractions — no
fabricated digits, no mpmath constants needed for the trap arithmetic):

    f(e) = 3/7 e^-2 - 9/5 e^-1 + 22/7 - 1/3 e + 7/9 e^2 + 5/11 e^3
           - 13/8 e^4  +  (1/1000) e^66            <- PLANTED ALIAS TRAP

The trap reproduces the measured frozen-grid disease in miniature: with the FROZEN
legacy configuration r=1/32, M=64 nodes, the trap order 66 = 2 + 64 aliases
into c_2 with weight (1/32)^64 * 10^-3 = 10^-99.33 — a dps-INDEPENDENT
~99-digit error floor (two-dps frozen reruns agree with each other and are
both wrong past it).  The DE-side singularity data places the nearest
eps-singularity at R0 = 1/4 (denominator 1 - 4*eps).

C1  planted-alias floor + break — frozen legacy (r=1/32, M=64) at dps=150
    reproduces the ~99d floor on c_2 (measured 90..110 d vs the exact
    Fraction); cauchy_laurent_certified at the SAME seeds escalates M
    geometrically (>= 256), breaks the floor, and recovers EVERY window
    coefficient to >= 145 d.
C2  R0-refusal — no R0 and no eps_polys raises ValueError naming R0
    (never silently assume); same for the Vandermonde variant and for
    eps_analyticity_radius on empty data.
C3  cap-raise — max_nodes=128 with r-shrink disabled cannot certify
    10^-(150+10): EpsFanCertifyError kind='alias-cap' with (bound, tol,
    r, M, R0) NAMED in .info and present in the message.
C4  cross-radius belt catch — non-analytic contamination
    + 1e-60*conj(e) passes the Cauchy alias bound (control-circle
    magnitudes unchanged) but corrupts c_{-1} radius-dependently
    (delta*r^2): kind='belt' raise at dps=80; the clean f at identical
    settings passes.
C5  R0 machinery — eps_analyticity_radius on exact data: [1,-4] -> 1/4;
    eps^m factors (roots AT 0) excluded exactly; constant polys -> inf;
    float coefficients refused; desystem_eps_polys duck-typed on de_load
    monomial-style entries (den poly in eps at fixed rational x) and
    NotImplementedError on unknown entry classes.
C6  shrink leg  — max_nodes=256 forces the escalation to shrink r below
    the seed (q = sqrt(r/R0) drops) and still certify: coefficients
    >= 140 d, diag r < r_seed.
C7  Vandermonde certified — clean exact-window fan at dps=60, R0=1/4:
    n_grid1=24 vs n_grid2=36 shifted (node count f(dps) seeded at 24),
    measured cond(V) drives wp above the seed (item 5), all window
    coefficients >= 55 d; contamination 50*e^4 beyond the window raises
    (vand-verify/vand-grid); max_wp cap raises kind='vand-wp-cap'.
C8  global mp.dps untouched (import-dps footgun).
"""
import os
import sys
from fractions import Fraction

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
sys.path.insert(0, PKG)

import epsfan  # noqa: E402

FAILS = []
DPS0 = mp.mp.dps

# KNOWN EXACT coefficients (window -2..4) + planted trap at order 66
COEFFS_EXACT = {
    -2: Fraction(3, 7),
    -1: Fraction(-9, 5),
    0: Fraction(22, 7),
    1: Fraction(-1, 3),
    2: Fraction(7, 9),
    3: Fraction(5, 11),
    4: Fraction(-13, 8),
}
TRAP_ORDER = 66            # 66 = 2 + 64: aliases into c_2 under frozen M=64
TRAP_COEFF = Fraction(1, 1000)
R0_EXACT = '1/4'           # denominator 1 - 4*eps -> nearest eps-singularity
EPS_POLYS = [[1, -4]]      # the same data in exact-poly form


def check(name, cond, msg=""):
    tag = "PASS" if cond else "FAIL"
    print("  [%s] %s%s" % (tag, name, (": " + msg) if msg else ""))
    if not cond:
        FAILS.append(name)


def frac_mpc(fr):
    return mp.mpf(fr.numerator) / mp.mpf(fr.denominator)


def f_trap(e):
    """Synthetic fan with the planted alias trap, exact rational
    coefficients, evaluated at ambient precision."""
    tot = mp.mpc(0)
    for k, c in COEFFS_EXACT.items():
        tot += frac_mpc(c) * mp.mpc(e) ** k
    tot += frac_mpc(TRAP_COEFF) * mp.mpc(e) ** TRAP_ORDER
    return tot


def err_digits(got, k, wp):
    """Measured agreement digits of got vs the EXACT coefficient c_k."""
    with mp.workdps(wp):
        exact = frac_mpc(COEFFS_EXACT[k])
        diff = abs(mp.mpc(got) - exact)
        if diff == 0:
            return float(wp)
        return float(-mp.log10(diff / abs(exact)))


def test_c1_planted_alias_floor_and_break():
    print("C1: planted alias trap — frozen (r=1/32, M=64) floor vs "
          "certified break (dps=150)")
    dps = 150
    # FROZEN legacy config: the ~99d-style floor (dps-independent alias)
    frozen = epsfan.cauchy_laurent(f_trap, -2, 4, dps, R='1/32', nodes=64)
    d2 = err_digits(frozen["coeffs"][2], 2, dps + 30)
    check("C1-frozen-floor", 90.0 <= d2 <= 110.0,
          "frozen c_2 true error %.2f d (predicted alias floor 99.33 d)"
          % d2)
    # clean orders are NOT floor-limited at the same frozen config
    d0 = err_digits(frozen["coeffs"][0], 0, dps + 30)
    check("C1-frozen-clean-order", d0 >= 130.0,
          "frozen c_0 %.1f d (floor is order-selective: 66 = 2 mod 64)" % d0)

    # certified engine, SAME seeds — must escalate M and break the floor
    out = epsfan.cauchy_laurent_certified(
        f_trap, -2, 4, dps, R0=R0_EXACT, r_seed='1/32', nodes_seed=64)
    dg = out["diag"]
    check("C1-M-escalated", dg["nodes"] >= 256,
          "M=%d from seed 64 (geometric escalation)" % dg["nodes"])
    worst = min(err_digits(out["coeffs"][k], k, dps + 30)
                for k in COEFFS_EXACT)
    check("C1-floor-broken", worst >= 145.0,
          "certified worst window coefficient %.2f d (>= 145 required, "
          "floor was ~99)" % worst)
    check("C1-diag-certified", dg.get("certified") is True
          and dg["R0"].startswith("0.25"))
    with mp.workdps(30):
        check("C1-alias-bound-below-tol",
              mp.mpmathify(dg["alias_bound"]) <=
              mp.mpmathify(dg["alias_tol"]),
              "bound %s vs tol %s" % (dg["alias_bound"], dg["alias_tol"]))
    check("C1-belt-honest", dg["belt_min_agree_digits"] >= dps,
          "belt min agree %.1f d" % dg["belt_min_agree_digits"])
    # eps_polys route (runtime R0 from exact singularity data) — same result
    out2 = epsfan.cauchy_laurent_certified(
        f_trap, -2, 4, dps, eps_polys=EPS_POLYS, r_seed='1/32',
        nodes_seed=64)
    check("C1-eps-polys-route",
          out2["diag"]["R0_source"] == "eps_polys"
          and out2["diag"]["R0"].startswith("0.25")
          and min(err_digits(out2["coeffs"][k], k, dps + 30)
                  for k in COEFFS_EXACT) >= 145.0)


def test_c2_r0_refusal():
    print("C2: R0 refusal — no singularity data and no caller R0")
    for fn, tag in ((epsfan.cauchy_laurent_certified, "cauchy"),
                    (epsfan.vandermonde_laurent_certified, "vand")):
        try:
            fn(f_trap, -2, 4, 40)
            check("C2-%s-refuses" % tag, False, "did NOT raise")
        except ValueError as e:
            check("C2-%s-refuses" % tag,
                  "R0" in str(e) and "REFUS" in str(e).upper(),
                  str(e)[:80])
    try:
        epsfan.eps_analyticity_radius([])
        check("C2-empty-polys", False, "did NOT raise")
    except ValueError as e:
        check("C2-empty-polys", "EMPTY" in str(e).upper())
    try:
        epsfan.eps_analyticity_radius(None)
        check("C2-none-polys", False, "did NOT raise")
    except ValueError:
        check("C2-none-polys", True)


def test_c3_cap_raise():
    print("C3: cap raise — max_nodes=128, shrink disabled, dps=150")
    try:
        epsfan.cauchy_laurent_certified(
            f_trap, -2, 4, 150, R0=R0_EXACT, r_seed='1/32', nodes_seed=64,
            max_nodes=128, max_r_halvings=0)
        check("C3-raises", False, "did NOT raise")
        return
    except epsfan.EpsFanCertifyError as e:
        check("C3-raises", True)
        check("C3-kind", e.kind == 'alias-cap', "kind=%r" % e.kind)
        info_ok = all(key in e.info for key in
                      ("bound", "tol", "r", "M", "R0"))
        check("C3-info-named", info_ok, "info keys: %s" % sorted(e.info))
        msg = str(e)
        check("C3-message-names-numbers",
              all(s in msg for s in ("bound=", "tol=", "r=", "M=", "R0=")),
              msg[:120])
        check("C3-M-at-cap", e.info["M"] == "128", "M=%s" % e.info["M"])


def test_c4_belt_catch():
    print("C4: cross-radius belt catches non-analytic contamination "
          "(dps=80)")
    delta = mp.mpf('1e-60')

    def f_contaminated(e):
        return f_trap(e) + delta * mp.conj(mp.mpc(e))

    # clean control at identical settings must PASS
    out = epsfan.cauchy_laurent_certified(
        f_trap, -2, 4, 80, R0=R0_EXACT, r_seed='1/32', nodes_seed=64)
    check("C4-clean-passes", out["diag"]["belt_min_agree_digits"] >= 80.0,
          "clean belt %.1f d" % out["diag"]["belt_min_agree_digits"])
    try:
        epsfan.cauchy_laurent_certified(
            f_contaminated, -2, 4, 80, R0=R0_EXACT, r_seed='1/32',
            nodes_seed=64)
        check("C4-belt-raises", False, "did NOT raise — the alias bound "
              "alone cannot see conj(e) contamination")
    except epsfan.EpsFanCertifyError as e:
        check("C4-belt-raises", True)
        check("C4-kind", e.kind == 'belt', "kind=%r" % e.kind)
        check("C4-names-k", e.info.get("k") == "-1",
              "corrupted order k=%s (conj(e) -> delta*r^2/e on the circle)"
              % e.info.get("k"))


def test_c5_r0_machinery():
    print("C5: eps_analyticity_radius + desystem_eps_polys")
    with mp.workdps(40):
        r = epsfan.eps_analyticity_radius([[1, -4]])
        check("C5-simple-root", abs(r - mp.mpf(1) / 4) < mp.mpf('1e-25'),
              "R0=%s (exact 1/4)" % mp.nstr(r, 12))
        # eps^2 * (1 - 4 eps): roots AT 0 excluded exactly
        r2 = epsfan.eps_analyticity_radius([[0, 0, 1, -4]])
        check("C5-zero-roots-excluded",
              abs(r2 - mp.mpf(1) / 4) < mp.mpf('1e-25'))
        # min over several polys, 'p/q' strings accepted
        r3 = epsfan.eps_analyticity_radius([[1, -4], ['1', '-8'], [7]])
        check("C5-min-over-polys",
              abs(r3 - mp.mpf(1) / 8) < mp.mpf('1e-25'),
              "R0=%s (exact 1/8)" % mp.nstr(r3, 12))
        # entire in eps -> inf
        check("C5-entire-inf",
              epsfan.eps_analyticity_radius([[5], [3, 0]]) == mp.inf)
    try:
        epsfan.eps_analyticity_radius([[1.0, -4.0]])
        check("C5-floats-refused", False, "did NOT raise")
    except ValueError as e:
        check("C5-floats-refused", "float" in str(e).lower())

    # duck-typed de_load monomial entries: den = (1 - 4 eps) + x*(...)
    class MonoEnt:
        def __init__(self, num, den):
            self.num, self.den = num, den

    class FakeSys:
        def __init__(self, entries):
            self._entries = entries

    # den(eps, x) = 1 - 4 eps - 2 x eps  ->  at x = 1/2:  1 - 5 eps
    ent = MonoEnt(num=[(0, 0, Fraction(1))],
                  den=[(0, 0, Fraction(1)), (1, 0, Fraction(-4)),
                       (1, 1, Fraction(-2))])
    polys = epsfan.desystem_eps_polys(FakeSys({(0, 0): ent}), '1/2')
    check("C5-desys-mono", polys == [[Fraction(1), Fraction(-5)]],
          "polys=%r" % polys)
    with mp.workdps(40):
        rr = epsfan.eps_analyticity_radius(polys)
        check("C5-desys-R0", abs(rr - mp.mpf(1) / 5) < mp.mpf('1e-25'),
              "R0=%s (exact 1/5)" % mp.nstr(rr, 12))

    class WeirdEnt:
        pass

    try:
        epsfan.desystem_eps_polys(FakeSys({(0, 0): WeirdEnt()}), 0)
        check("C5-unknown-entry-raises", False, "did NOT raise")
    except NotImplementedError as e:
        check("C5-unknown-entry-raises", "WeirdEnt" in str(e))


def test_c6_shrink_lane():
    print("C6: r-shrink leg — max_nodes=256 forces r below the seed "
          "(dps=150)")
    out = epsfan.cauchy_laurent_certified(
        f_trap, -2, 4, 150, R0=R0_EXACT, r_seed='1/32', nodes_seed=64,
        max_nodes=256, max_r_halvings=8)
    with mp.workdps(30):
        r_final = mp.mpmathify(out["diag"]["r"])
        check("C6-r-shrunk", r_final < mp.mpf(1) / 32,
              "r=%s < seed 1/32, M=%d" % (out["diag"]["r"],
                                          out["diag"]["nodes"]))
    worst = min(err_digits(out["coeffs"][k], k, 180) for k in COEFFS_EXACT)
    check("C6-still-certified", worst >= 140.0,
          "worst window coefficient %.2f d" % worst)


def f_window(e):
    """Exact-window fan (orders -2..3 only) for the Vandermonde variant."""
    tot = mp.mpc(0)
    for k, c in COEFFS_EXACT.items():
        if k <= 3:
            tot += frac_mpc(c) * mp.mpc(e) ** k
    return tot


def test_c7_vandermonde_certified():
    print("C7: Vandermonde certified — two-grid gate + cond monitor "
          "(dps=60)")
    dps = 60
    out = epsfan.vandermonde_laurent_certified(
        f_window, -2, 3, dps, R0=R0_EXACT)
    dg = out["diag"]
    check("C7-grids-24-36", dg["n_grid1"] == 24 and dg["n_grid2"] == 36,
          "n1=%d n2=%d (f(dps) seeded at 24)" % (dg["n_grid1"],
                                                 dg["n_grid2"]))
    worst = min(err_digits(out["coeffs"][k], k, dps + 30)
                for k in COEFFS_EXACT if k <= 3)
    check("C7-recovers", worst >= 55.0, "worst coefficient %.2f d" % worst)
    check("C7-wp-from-cond", dg["wp"] >= dps + 4 * 6 + 20 + 10,
          "wp=%d, measured cond=%s (item 5: conditioning-driven margin)"
          % (dg["wp"], dg["cond"]))
    check("C7-gate-digits", dg["grid_min_agree_digits"] >= dps,
          "two-grid min agree %.1f d" % dg["grid_min_agree_digits"])

    # beyond-window contamination: the two grids (and held-out spares)
    # disagree -> fail-closed raise, never a silent wrong window
    def f_contam(e):
        return f_window(e) + 50 * mp.mpc(e) ** 4

    try:
        epsfan.vandermonde_laurent_certified(f_contam, -2, 3, dps,
                                             R0=R0_EXACT)
        check("C7-contam-raises", False, "did NOT raise")
    except epsfan.EpsFanCertifyError as e:
        check("C7-contam-raises", e.kind in ('vand-grid', 'vand-verify'),
              "kind=%r" % e.kind)

    # conditioning-driven wp above a crippled cap -> named cap raise
    try:
        epsfan.vandermonde_laurent_certified(f_window, -2, 3, dps,
                                             R0=R0_EXACT, max_wp=dps + 30)
        check("C7-wp-cap-raises", False, "did NOT raise")
    except epsfan.EpsFanCertifyError as e:
        check("C7-wp-cap-raises", e.kind == 'vand-wp-cap',
              "kind=%r" % e.kind)
        check("C7-wp-cap-named",
              all(k in e.info for k in ("cond", "wp_req", "max_wp", "R0")),
              "info keys: %s" % sorted(e.info))


def test_zz_all_checks_enforced():
    """pytest enforcement (runs last): the check() accumulator must be
    empty, and global mp.dps untouched (import-dps footgun)."""
    assert not FAILS, "epsfan certified self-checks failed: %s" % \
        ",".join(FAILS)
    assert mp.mp.dps == DPS0, "global mp.dps mutated (import-dps footgun)"


if __name__ == "__main__":
    print("== test_epsfan_certified (E4 certified layer) ==")
    test_c1_planted_alias_floor_and_break()
    test_c2_r0_refusal()
    test_c3_cap_raise()
    test_c4_belt_catch()
    test_c5_r0_machinery()
    test_c6_shrink_lane()
    test_c7_vandermonde_certified()
    check("C8-global-dps-untouched", mp.mp.dps == DPS0,
          "dps %d -> %d" % (DPS0, mp.mp.dps))
    print("== %s ==" % ("ALL PASS" if not FAILS
                        else "FAILURES: " + ",".join(FAILS)))
    sys.exit(0 if not FAILS else 1)
