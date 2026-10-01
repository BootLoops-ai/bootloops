#!/usr/bin/env python3
"""
test_fastpath.py — DESystem.enable_fast_path (the A_series/singular_points
protocol).

Discipline: references are exact Fractions computed in-test from the same
closed forms already pinned by test_de_load.py (the kira vacuum-bubble
hand form, monomial-json spot entries) — no pasted digits. dps <= 120 (Verify phase owns
heavy controls). Every check would fail if the fast path aliased, zero-padded,
ran at the wrong eps, or broke the plain .A() path.
"""
import os
import sys
from fractions import Fraction as F

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from de_load import load_monomial_json, load_kira_targets, load_amatrix_json  # noqa: E402
from transport import transport_fixed_eps  # noqa: E402

# Optional reference fixtures, not shipped (legs SKIP when absent).
_FIX = os.environ.get("WAYFINDER_REFERENCE_FIXTURES", "")
CONE = os.path.join(_FIX, "cone_As_t13.json")
AOFY = os.path.join(_FIX, "A_of_y.json")
KIRA6 = os.path.join(_FIX, "kira_target.m")


def _require(path):
    if not os.path.isfile(path):
        print(f"  SKIP (missing {path})")
        return False
    return True


def _digits(got, want_re, want_im, dps):
    with mp.workdps(dps + 20):
        want = mp.mpc(mp.mpf(want_re.numerator) / want_re.denominator,
                      mp.mpf(want_im.numerator) / want_im.denominator)
        diff = abs(got - want)
        if diff == 0:
            return dps + 20
        return int(-mp.log10(diff / max(abs(want), mp.mpf(1))))


def test_kira_fastpath_series_and_sings():
    """system_6 closed form (pinned in test_de_load):
        A = [[(d-3)/(eta+1), (2-d)/(2*eta^2+2*eta)], [0, (d-2)/(2*eta)]].
    Taylor about eta0: A11 coeffs (d-3)(-1)^m/(1+eta0)^{m+1};
    A22 coeffs (d-2)/2 * (-1)^m/eta0^{m+1} — exact Fractions here.
    Singular points must include eta = 0 and eta = -1."""
    if not _require(KIRA6):
        return
    desys = load_kira_targets(KIRA6)
    assert not hasattr(desys, "A_series"), \
        "plain DESystem must NOT expose A_series (transport duck-typing)"
    assert not hasattr(desys, "singular_points")
    eps = F(1, 3)
    d = F(4) - 2 * eps
    desys.enable_fast_path(eps, wp=120)
    assert hasattr(desys, "A_series") and hasattr(desys, "singular_points")
    dps = 60
    M = 5
    for eta0 in (F(3), F(5, 7)):
        ser = desys.A_series(eta0, eps, dps, M)
        assert len(ser) == M + 1
        rows = desys.A(eta0, eps, dps)  # plain path unchanged
        for m in range(M + 1):
            c11 = (d - 3) * F(-1) ** m / (1 + eta0) ** (m + 1)
            c22 = (d - 2) / 2 * F(-1) ** m / eta0 ** (m + 1)
            assert _digits(ser[m][0][0], c11, F(0), dps) >= dps - 3, \
                (eta0, m, ser[m][0][0], c11)
            assert _digits(ser[m][1][1], c22, F(0), dps) >= dps - 3, \
                (eta0, m, ser[m][1][1], c22)
            assert ser[m][1][0] == 0
        # order 0 == pointwise A
        with mp.workdps(dps):
            for i in range(2):
                for j in range(2):
                    dd = abs(mp.mpc(ser[0][i][j]) - mp.mpc(rows[i][j]))
                    assert dd <= mp.mpf(10) ** (-(dps - 3)) * \
                        max(1, abs(mp.mpc(rows[i][j]))), (i, j)
    with mp.workdps(30):
        pts = [mp.mpc(s) for s in desys.singular_points]
        for want in (mp.mpc(0), mp.mpc(-1)):
            assert any(abs(p - want) < mp.mpf("1e-15") for p in pts), \
                (want, desys.singular_points)
    print("  PASS kira fast path: Taylor coeffs == closed form (2 pts, M=5), "
          "sings {0,-1} found")


def test_kira_fastpath_transport_equivalence():
    """transport_fixed_eps through the fast path must agree with the plain
    (Cauchy-circle, two-radius cross-checked) path — the alias-safe claim."""
    if not _require(KIRA6):
        return
    eps = F(1, 3)
    dps = 50
    plain = load_kira_targets(KIRA6)
    fast = load_kira_targets(KIRA6).enable_fast_path(eps, wp=120)
    y0 = ["0.25", "1"]
    a = transport_fixed_eps(plain, eps, "3", "1", y0, dps)
    b = transport_fixed_eps(fast, eps, "3", "1", y0, dps)
    with mp.workdps(dps + 10):
        for i in range(2):
            diff = abs(mp.mpc(a[i]) - mp.mpc(b[i]))
            scale = max(abs(mp.mpc(a[i])), mp.mpf(1))
            assert diff <= scale * mp.mpf(10) ** (-(dps - 5)), \
                (i, a[i], b[i])
    print("  PASS kira fast path == plain path through transport_fixed_eps")


def test_fastpath_guards():
    """Wrong-eps and over-dps calls must raise (no silent wrong data)."""
    if not _require(KIRA6):
        return
    desys = load_kira_targets(KIRA6).enable_fast_path(F(1, 3), wp=100)
    try:
        desys.A_series(F(2), F(1, 5), 40, 2)
        raise AssertionError("wrong-eps A_series call did not raise")
    except ValueError as e:
        assert "eps" in str(e)
    try:
        desys.A_series(F(2), F(1, 3), 300, 2)
        raise AssertionError("dps > wp A_series call did not raise")
    except ValueError as e:
        assert "wp" in str(e)
    # eps differing only at rounding level (dyadic parse of 1/3) must PASS
    with mp.workdps(80):
        eps_dyadic = mp.mpf(1) / 3
    ser = desys.A_series(F(2), eps_dyadic, 60, 1)
    assert len(ser) == 2
    print("  PASS fast-path guards (wrong eps / over-dps raise; dyadic eps ok)")


def test_monomial_fastpath():
    """cone_As_t13 spot entries (pinned in test_de_load):
        A[0][0] = -eps/s        -> Taylor about s0: -eps (-1)^m / s0^{m+1}
        A[4][4] den = s^2 - 4 s -> singular points include 0 and 4."""
    if not _require(CONE):
        return
    desys = load_monomial_json(CONE, dps_check=0)
    eps = F(1, 7)
    desys.enable_fast_path(eps, wp=100)
    dps = 50
    s0 = F(1, 3)
    ser = desys.A_series(s0, eps, dps, 3)
    for m in range(4):
        want = -eps * F(-1) ** m / s0 ** (m + 1)
        assert _digits(ser[m][0][0], want, F(0), dps) >= dps - 3, \
            (m, ser[m][0][0], want)
    with mp.workdps(30):
        pts = [mp.mpc(s) for s in desys.singular_points]
        for want in (mp.mpc(0), mp.mpc(4)):
            assert any(abs(p - want) < mp.mpf("1e-12") for p in pts), want
    print("  PASS monomial fast path: -eps/s Taylor + sings {0,4}")


def test_amatrix_fastpath_spot():
    """A_of_y fast path: order-0 A_series == plain A at y=1/2, d0=4-2/1009
    (the reference q=1009 point, same as test_de_load)."""
    if not _require(AOFY):
        return
    desys = load_amatrix_json(AOFY)
    eps = F(1, 1009)
    desys.enable_fast_path(eps, wp=100)
    dps = 50
    y0 = F(1, 2)
    ser = desys.A_series(y0, eps, dps, 1)
    rows = desys.A(y0, eps, dps)
    with mp.workdps(dps):
        worst = 0
        for i in range(desys.n):
            for j in range(desys.n):
                dd = abs(mp.mpc(ser[0][i][j]) - mp.mpc(rows[i][j]))
                sc = max(abs(mp.mpc(rows[i][j])), mp.mpf(1))
                assert dd <= sc * mp.mpf(10) ** (-(dps - 5)), (i, j)
    assert len(desys.singular_points) > 0
    print(f"  PASS amatrix fast path: order-0 == A on all {desys.n}^2 "
          f"entries; {len(desys.singular_points)} singular points")


def test_no_global_dps_mutation():
    if not _require(KIRA6):
        return
    before = mp.mp.dps
    desys = load_kira_targets(KIRA6).enable_fast_path(F(1, 4), wp=80)
    desys.A_series(F(3), F(1, 4), 60, 4)
    assert mp.mp.dps == before, "global mp.dps mutated by fast path"
    print("  PASS fast path leaves global mp.dps untouched")


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        print(t.__name__)
        t()
    print(f"ALL {len(tests)} fast-path test groups done")


if __name__ == "__main__":
    main()
