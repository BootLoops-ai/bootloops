"""CP machinery vs published tables + closed forms.

Closed forms (Clopper & Pearson 1934, k = 0 case): one-sided upper
p_U = 1 - alpha^(1/n); two-sided upper = 1 - (alpha/2)^(1/n).
Published two-sided 95% CP values (standard tables, 4 dp; also the
beta-quantile identity p_U = BetaInv(1 - alpha/2, k+1, n-k)):
  (k=5, n=10): (0.1871, 0.8129)   (k=10, n=20): (0.2720, 0.7280)
  (k=1, n=20): (0.0013, 0.2487)
"""
import os
import sys
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from cp_bound import (ALPHA0, binom_cdf, cp_interval, cp_lower, cp_on_fctl,
                      cp_upper, exceedance, min_n_for_significance)


def test_binom_cdf_exact():
    assert binom_cdf(1, 2, Fraction(1, 2)) == Fraction(3, 4)
    assert binom_cdf(0, 10, Fraction(1, 4)) == Fraction(3, 4) ** 10
    assert binom_cdf(10, 10, Fraction(1, 3)) == 1
    assert binom_cdf(-1, 5, Fraction(1, 2)) == 0


def test_k0_closed_forms():
    lo, hi = cp_upper(0, 10)                     # one-sided 95%
    v = 1 - 0.05 ** 0.1                          # 0.258866
    assert float(lo) <= v <= float(hi) + 1e-12
    assert abs(float(hi) - v) < 1e-9
    lo2, hi2 = cp_upper(0, 100)                  # 0.029513
    assert abs(float(hi2) - (1 - 0.05 ** 0.01)) < 1e-9
    itv = cp_interval(0, 10)                     # (0, 0.308497)
    assert itv[0] == 0
    assert abs(float(itv[1]) - (1 - 0.025 ** 0.1)) < 1e-9


def test_published_two_sided_tables():
    refs = [(5, 10, 0.1871, 0.8129), (10, 20, 0.2720, 0.7280),
            (1, 20, 0.0013, 0.2487)]
    for k, n, lo_ref, hi_ref in refs:
        lo, hi = cp_interval(k, n)
        assert abs(float(lo) - lo_ref) < 6e-5, (k, n, float(lo))
        assert abs(float(hi) - hi_ref) < 6e-5, (k, n, float(hi))


def test_edge_cases_and_exceedance():
    assert cp_upper(10, 10) == (1, 1)
    assert cp_lower(0, 10) == (0, 0)
    assert exceedance(Fraction(1, 2), 10) == Fraction(1, 1024)
    assert min_n_for_significance(Fraction(1, 2)) == 10
    n66 = min_n_for_significance(Fraction(9, 10))
    assert n66 == 66
    assert Fraction(9, 10) ** 66 < ALPHA0 <= Fraction(9, 10) ** 65
    assert min_n_for_significance(Fraction(1)) is None
    assert cp_on_fctl(0, 10) == cp_upper(0, 10)[1]


def test_beta_identity_crosscheck():
    try:
        from scipy.stats import beta
    except ImportError:
        return  # optional cross-check only; exact tests above are binding
    for k, n in [(3, 50), (0, 100), (7, 12)]:
        ref = float(beta.ppf(0.95, k + 1, n - k))
        assert abs(float(cp_upper(k, n)[1]) - ref) < 1e-9, (k, n)
