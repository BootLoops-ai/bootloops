#!/usr/bin/env python3
"""
Tests for tools/wayfinder/de_load.py.

Reference values are NEVER pasted-in digits: every expected number is derived
in-test from exact Fraction arithmetic on formulas read by eye from the raw
source files (cited per test), or from independent raw-JSON re-evaluation
through a second code path. Fast (dps <= 120, seconds) — heavy controls are
the Verify phase's job.

Run:  python3 tests/test_de_load.py        (or pytest tests/test_de_load.py)
"""
import json
import os
import sys
from fractions import Fraction as F

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import de_load
from de_load import (load_monomial_json, load_amatrix_json, load_kira_targets)

# Optional reference fixtures (large production DE files, not shipped): point
# WAYFINDER_REFERENCE_FIXTURES at a directory holding the named files to run
# them; without it these legs SKIP (the format loaders are still exercised by
# the synthetic legs below).
_FIX = os.environ.get("WAYFINDER_REFERENCE_FIXTURES", "")
CONE = os.path.join(_FIX, "cone_As_t13.json")
CLOSURE = os.path.join(_FIX, "cone_closure.json")
AOFY = os.path.join(_FIX, "A_of_y.json")
KIRA6 = os.path.join(_FIX, "kira_target.m")


def _agree_digits(got_mpc, want_frac_re, want_frac_im, dps):
    """Digits of agreement between an mpc and an exact Gaussian rational."""
    with mp.workdps(dps + 20):
        want = mp.mpc(mp.mpf(want_frac_re.numerator) / want_frac_re.denominator,
                      mp.mpf(want_frac_im.numerator) / want_frac_im.denominator)
        diff = abs(got_mpc - want)
        scale = max(abs(want), mp.mpf(1))
        if diff == 0:
            return dps + 20
        return int(-mp.log10(diff / scale))


def _require(path):
    if not os.path.isfile(path):
        print(f"  SKIP (missing {path})")
        return False
    return True


# ---------------------------------------------------------------------------
# monomial-json (cone_As_t13.json)
# ---------------------------------------------------------------------------

def test_monomial_load_and_spot_values():
    """Spot formulas read by eye from the raw JSON:
        A_s['0,0']: num [[1,0,'-1']], den [[0,1,'1']]      -> -eps / s
        A_s['4,4']: num [[1,1,-2],[1,0,-4],[0,1,2],[0,0,4]],
                    den [[0,2,1],[0,1,-4]]
                    -> (-2*eps*s - 4*eps + 2*s + 4) / (s^2 - 4*s)
    computed independently here with Fractions."""
    if not _require(CONE):
        return
    desys = load_monomial_json(CONE, dps_check=50)
    assert desys.n == 22, desys.n
    assert desys.var == "s", desys.var
    assert desys.meta["format"] == "monomial-json"
    assert len(desys.meta["sha256"]) == 64
    dps = 60
    for (s, eps) in ((F(1, 3), F(1, 7)), (F(-2, 5), F(3, 11))):
        rows = desys.A(s, eps, dps)
        assert len(rows) == 22 and all(len(r) == 22 for r in rows)
        exp00 = -eps / s
        exp44 = (F(-2) * eps * s - 4 * eps + 2 * s + 4) / (s * s - 4 * s)
        assert _agree_digits(rows[0][0], exp00, F(0), dps) >= dps - 2, \
            (s, eps, rows[0][0], exp00)
        assert _agree_digits(rows[4][4], exp44, F(0), dps) >= dps - 2, \
            (s, eps, rows[4][4], exp44)
        # exact-path cross-check (Fraction, no rounding at all)
        q00 = desys.entry_exact(0, 0, s, eps)
        assert q00.re == exp00 and q00.im == 0, (q00, exp00)
        q44 = desys.entry_exact(4, 4, s, eps)
        assert q44.re == exp44 and q44.im == 0, (q44, exp44)
    print("  PASS monomial spot values (2 rational points, exact + rounded)")


def test_monomial_denominators_nonzero():
    """Full-matrix eval at 2 rational points: any vanishing denominator
    raises ZeroDivisionError, so a clean pass IS the nonzero check."""
    if not _require(CONE):
        return
    desys = load_monomial_json(CONE, dps_check=0)
    for (s, eps) in ((F(1, 3), F(1, 7)), (F(-2, 5), F(3, 11))):
        rows = desys.A(s, eps, 40)
        with mp.workdps(40):
            for r in rows:
                for v in r:
                    assert mp.isfinite(v.real) and mp.isfinite(v.imag)
    print("  PASS monomial denominators nonzero at 2 rational points")


def test_monomial_complex_point_exact():
    """Gaussian-rational path: at s = 1/3 + i/5 (exact), A[0][0] = -eps/s,
    computed here with exact complex-Fraction arithmetic."""
    if not _require(CONE):
        return
    desys = load_monomial_json(CONE, dps_check=0)
    eps = F(1, 7)
    sre, sim = F(1, 3), F(1, 5)
    d2 = sre * sre + sim * sim
    exp_re = -eps * sre / d2          # -eps * conj(s) / |s|^2
    exp_im = eps * sim / d2
    dps = 60
    # exact Gaussian input via the _QQi hook (complex() literals are dyadic-
    # only, so 1/3 + i/5 must go in as exact Fractions)
    q = desys.entry_exact(0, 0, de_load._QQi(sre, sim), eps)
    assert q.re == exp_re and q.im == exp_im, (q, exp_re, exp_im)
    # dyadic complex input through the public API: s = 1/4 + i/2 is exact
    rows = desys.A(complex(0.25, 0.5), eps, dps)
    sre2, sim2 = F(1, 4), F(1, 2)
    d22 = sre2 * sre2 + sim2 * sim2
    got = rows[0][0]
    assert _agree_digits(got, -eps * sre2 / d22, eps * sim2 / d22, dps) >= dps - 2
    print("  PASS monomial complex-point (exact Gaussian + dyadic mpc path)")


def test_monomial_mpf_input_matches_fraction():
    """mpf('0.375') is dyadic == 3/8 exactly: both inputs must agree exactly."""
    if not _require(CONE):
        return
    desys = load_monomial_json(CONE, dps_check=0)
    dps = 50
    with mp.workdps(30):
        x_mpf = mp.mpf("0.375")
    a = desys.A(x_mpf, F(1, 7), dps)
    b = desys.A(F(3, 8), F(1, 7), dps)
    with mp.workdps(dps):
        for i in (0, 4):
            assert a[i][i] == b[i][i], (i, a[i][i], b[i][i])
    print("  PASS mpf dyadic input == Fraction input")


def test_monomial_closure_aux_and_no_global_dps():
    if not _require(CONE) or not _require(CLOSURE):
        return
    before = mp.mp.dps
    desys = load_monomial_json(CONE, dps_check=40, closure_json=CLOSURE)
    assert "closure" in desys.meta and "tcone_to_scone" in desys.meta["closure"]
    desys.A(F(1, 3), F(1, 7), 90)
    assert mp.mp.dps == before, f"global dps mutated: {before} -> {mp.mp.dps}"
    print("  PASS closure aux loaded; global mp.dps untouched")


# ---------------------------------------------------------------------------
# amatrix-json (A_of_y.json)
# ---------------------------------------------------------------------------

def _raw_amatrix_entry(rec, y, d):
    """INDEPENDENT re-evaluation of one A_of_y entry straight from raw JSON
    (schema: fit_A_of_y.py lines 110-112; low->high Horner = its ev())."""
    def ev(coeffs, x):
        r = F(0)
        for c in reversed(coeffs):
            r = r * x + c
        return r

    def side(items):
        return sum((ev([F(c) for c in it["num"]], d)
                    / ev([F(c) for c in it["den"]], d)) * y ** k
                   for k, it in enumerate(items))

    return side(rec["N"]) / side(rec["Q"])


def test_amatrix_load_and_crosscheck():
    if not _require(AOFY):
        return
    desys = load_amatrix_json(AOFY)
    assert desys.n == 79, desys.n     # the reference 79-master system
    assert desys.var == "y"
    assert len(desys.meta["masters_hash"]) == 64
    assert desys.meta["format"] == "amatrix-json"
    with open(AOFY) as f:
        raw = json.load(f)
    y, eps = F(1, 2), F(1, 1009)      # d0 = 4 - 2/1009: the reference q=1009 point
    d = F(4) - 2 * eps
    dps = 80
    rows = desys.A(y, eps, dps)
    for key in ("10,0", "0,0" if "0,0" in raw["entries"] else "10,3"):
        if key not in raw["entries"]:
            continue
        i, j = (int(t) for t in key.split(","))
        want = _raw_amatrix_entry(raw["entries"][key], y, d)
        assert _agree_digits(rows[i][j], want, F(0), dps) >= dps - 2, \
            (key, rows[i][j], want)
        qq = desys.entry_exact(i, j, y, eps)
        assert qq.re == want and qq.im == 0, (key, qq.re, want)
    # unstored key must be exact zero
    allkeys = set(raw["entries"])
    for i in range(79):
        for j in range(79):
            if f"{i},{j}" not in allkeys:
                with mp.workdps(dps):
                    assert rows[i][j] == mp.mpc(0)
                break
        else:
            continue
        break
    print("  PASS amatrix n=79, raw-JSON cross-check exact at d0=4-2/1009")


# ---------------------------------------------------------------------------
# kira-target-m (eta-DE)
# ---------------------------------------------------------------------------

def test_kira_system6_closed_form():
    """system_6 is the single-mass vacuum bubble subsystem:
        propagators D1 = k2^2 - eta, D2 = (k2-k1)^2  (config yaml)
        masters (file order): [1,1], [1,0]
        kira_target.m:
          [2,1] -> [1,0]*((-d+2)/(2*eta^2+2*eta)) + [1,1]*((d-3)/(eta+1))
          [2,0] -> [1,0]*((d-2)/(2*eta))
    With d/deta I[a] = +a_1 I[a+e_1] (c_1 = -1) the eta-DE is EXACTLY
        A = [[ (d-3)/(eta+1),  (-d+2)/(2*eta^2+2*eta) ],
             [ 0,               (d-2)/(2*eta)          ]]
    (basis order [1,1], [1,0]). Independent physics check: the [1,0] row is
    the massive-tadpole scaling d/deta log(eta^(d/2-1)) = (d-2)/(2*eta)."""
    if not _require(KIRA6):
        return
    desys = load_kira_targets(KIRA6)
    assert desys.n == 2, desys.n
    assert desys.var == "eta"
    assert desys.meta["format"] == "kira-target-m"
    assert desys.meta["eta_coefficients"] == [-1, 0]
    assert desys.meta["masters"][0].endswith("[1,1]")
    assert desys.meta["masters"][1].endswith("[1,0]")
    dps = 60
    for (eta, eps) in ((F(3), F(1, 2)), (F(5, 7), F(1, 3))):
        d = F(4) - 2 * eps
        exp = [[(d - 3) / (eta + 1), (-d + 2) / (2 * eta * eta + 2 * eta)],
               [F(0), (d - 2) / (2 * eta)]]
        rows = desys.A(eta, eps, dps)
        for i in range(2):
            for j in range(2):
                assert _agree_digits(rows[i][j], exp[i][j], F(0), dps) >= dps - 2, \
                    (eta, eps, i, j, rows[i][j], exp[i][j])
        # exact equality through the exact hook
        q = desys.entry_exact(1, 1, eta, eps)
        assert q.re == (d - 2) / (2 * eta) == (d / 2 - 1) / eta  # tadpole scaling
        assert desys.entry_exact(1, 0, eta, eps).is_zero()
    print("  PASS kira system_6 eta-DE == hand closed form (2 points, exact)")


def test_kira_complex_eta_point():
    """Waypoint-style complex eta (dyadic 1/2 + i/4) through the public API."""
    if not _require(KIRA6):
        return
    desys = load_kira_targets(KIRA6)
    eps = F(1, 5)
    d = F(4) - 2 * eps
    dps = 50
    rows = desys.A(complex(0.5, 0.25), eps, dps)
    # A[1][1] = (d-2)/(2*eta), eta = 1/2 + i/4, exact Gaussian arithmetic:
    er, ei = F(1, 2), F(1, 4)
    n2 = 4 * (er * er + ei * ei)      # |2*eta|^2
    exp_re = (d - 2) * 2 * er / n2
    exp_im = -(d - 2) * 2 * ei / n2
    assert _agree_digits(rows[1][1], exp_re, exp_im, dps) >= dps - 2
    print("  PASS kira complex-eta waypoint point")


def test_error_paths():
    # unknown scalar type
    if _require(CONE):
        desys = load_monomial_json(CONE, dps_check=0)
        try:
            desys.A(object(), F(1, 7), 30)
            raise AssertionError("expected TypeError")
        except TypeError:
            pass
        try:
            desys.A(F(1, 3), F(1, 7), 1)
            raise AssertionError("expected ValueError on dps=1")
        except ValueError:
            pass
    # kira grammar guard: unparseable construct must raise NotImplementedError
    try:
        de_load._parse_kira_rules("{fam[1,0] -> Sqrt[fam[2,0]]}", "<test>")
        raise AssertionError("expected NotImplementedError")
    except NotImplementedError:
        pass
    # eta-placement guard
    try:
        de_load._eta_coefficient("k1^2 - 2*eta")
        raise AssertionError("expected NotImplementedError")
    except NotImplementedError:
        pass
    assert de_load._eta_coefficient("(k2-k1)^2") == 0
    assert de_load._eta_coefficient("k2^2 - eta") == -1
    assert de_load._eta_coefficient("k2^2-eta") == -1
    assert de_load._eta_coefficient("eta") == 1
    # production shapes: eta additive but NOT trailing (SingleMass
    # insertion on an already-massive line: mass -> -(eta+1))
    assert de_load._eta_coefficient(
        "(l1 - l2 + p1 + p2 + p3)^2 - eta-1") == -1
    assert de_load._eta_coefficient("(l2 - p2 - p3)^2 - eta - 1") == -1
    assert de_load._eta_coefficient("l1^2 + eta - 1") == 1
    for bad in ("k1^2 - eta*s", "k1^2 - eta^2", "(a - eta)*b"):
        try:
            de_load._eta_coefficient(bad)
            raise AssertionError(f"expected NotImplementedError on {bad!r}")
        except NotImplementedError:
            pass
    print("  PASS error paths (bad scalar, bad dps, grammar + eta guards)")


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        print(t.__name__)
        t()
    print(f"ALL {len(tests)} de_load test groups done")


if __name__ == "__main__":
    main()
