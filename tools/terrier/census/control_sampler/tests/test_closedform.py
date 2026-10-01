"""Closed-form acceptance-rate battery (synthetic, fixed seeds).

Synthetic AD-form density (2 - x)(1 + y) = det(-R - omega 1) with
omega = I, R(z) = [[x - 3, 0], [0, -2 - y]] (all exact rationals).
Box [0,1]^2:  integral = 9/4 exact, envelope M = 4,
              overall acceptance = (9/4)/(4 * 1) = 9/16 EXACT.
Triangle {0 <= y <= x < 1}: integral = 7/8 exact,
              overall acceptance = 7/32, in-domain fraction = 1/2.
"""
import os
import sys
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import seeds as S
from ad_density import ADDensity, det_fraction
from domain import BoxDomain, PredicateDomain
from rejection import sample

SEEDS_JSON = os.environ.get("TERRIER_SEEDS_JSON", "")
if not SEEDS_JSON:
    raise SystemExit("REFUSE: env TERRIER_SEEDS_JSON unset — must point at "
                     "the seeds JSON file (reference data not included in "
                     "the package)")


def _R(z):
    x, y = z
    return [[x - 3, 0], [0, -2 - y]]


def _W(z):
    return [[1, 0], [0, 1]]


DENS = ADDensity(2, _R, _W)


def _gen(bin_label):
    return S.bin_generator(S.load_master(SEEDS_JSON), "SYNTH", bin_label)


def test_density_exact_values():
    assert DENS([Fraction(1, 2), Fraction(1, 2)]) == Fraction(9, 4)
    assert DENS([Fraction(0), Fraction(0)]) == 2
    assert DENS([Fraction(1), Fraction(1)]) == 2
    assert det_fraction([[2, 1], [1, 1]]) == 1
    assert det_fraction([[0, 1], [1, 0]]) == -1
    assert det_fraction([[1, 2], [2, 4]]) == 0


def test_exact_closed_form_rates():
    # the closed-form arithmetic itself, exact (no sampling):
    assert Fraction(9, 4) / (4 * BoxDomain([0, 0], [1, 1]).volume()) \
        == Fraction(9, 16)
    assert Fraction(7, 8) / (4 * 1) == Fraction(7, 32)


def test_box_acceptance_matches_closed_form():
    pts, c = sample(BoxDomain([0, 0], [1, 1]), DENS, 4, _gen(1), 2000)
    assert c["out_of_domain"] == 0
    assert c["accepted"] == 2000 and len(pts) == 2000
    p, n = 9.0 / 16.0, c["proposed"]
    p_hat = c["accepted"] / n
    assert abs(p_hat - p) <= 4 * (p * (1 - p) / n) ** 0.5, (p_hat, n)


def test_triangle_acceptance_and_domain_rejection():
    tri = PredicateDomain([0, 0], [1, 1], lambda z: z[1] <= z[0])
    pts, c = sample(tri, DENS, 4, _gen(2), 1500)
    n = c["proposed"]
    p = 7.0 / 32.0
    assert abs(c["accepted"] / n - p) <= 4 * (p * (1 - p) / n) ** 0.5
    q = 0.5  # in-domain fraction of the bounding box
    assert abs(c["out_of_domain"] / n - q) <= 4 * (q * (1 - q) / n) ** 0.5
    assert all(tri.contains(x) for x in pts)
    assert all(isinstance(v, Fraction) for x in pts for v in x)


def test_envelope_violation_halts():
    try:
        sample(BoxDomain([0, 0], [1, 1]), DENS, 2, _gen(3), 500)
    except ValueError as e:
        assert "envelope violated" in str(e)
    else:
        raise AssertionError("bad envelope M=2 < sup f = 4 did not halt")


def test_negative_density_halts():
    bad = ADDensity(2, lambda z: [[3, 0], [0, -2]], _W)  # det(-R-I) = -4
    try:
        bad([Fraction(1, 2), Fraction(1, 2)])
    except ValueError as e:
        assert "sign law" in str(e)
    else:
        raise AssertionError("negative index density did not halt")
