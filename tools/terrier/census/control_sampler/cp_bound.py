"""Clopper-Pearson exact binomial bounds.

CP-on-F_ctl: the H0-a exceedance is evaluated at the Clopper-Pearson upper
95% bound on F_ctl, the on-locus fraction of the N_ctl-point control
rejection sample (drawn at fixed seeds). This module is that machinery:
exact binomial
arithmetic (Fraction CDF via term-ratio recurrence, dyadic bisection)
returning CERTIFIED BRACKETS [lo, hi] that provably contain the CP root;
quote hi for an upper bound, lo for a lower bound (conservative both ways).
"""
from fractions import Fraction
from math import log

ALPHA0 = Fraction(1, 1000)  # pinned significance floor alpha0 (power gate)


def binom_cdf(k, n, p):
    """P(Bin(n, p) <= k), exact Fraction; term-ratio recurrence."""
    p = Fraction(p)
    if k < 0:
        return Fraction(0)
    if k >= n:
        return Fraction(1)
    if p == 0:
        return Fraction(1)
    if p == 1:
        return Fraction(0)
    q = 1 - p
    term = q ** n
    s = term
    for i in range(1, k + 1):
        term = term * (n - i + 1) * p / (Fraction(i) * q)
        s += term
    return s


def cp_upper(k, n, conf=Fraction(95, 100), tol=Fraction(1, 2 ** 40)):
    """One-sided CP upper bound: smallest p with P(X <= k | n, p) <= 1-conf.
    Returns (lo, hi), root in [lo, hi], hi - lo <= tol. k >= n -> (1, 1).
    binom_cdf(k, n, .) is strictly decreasing in p on (0, 1)."""
    alpha = 1 - Fraction(conf)
    if k >= n:
        return Fraction(1), Fraction(1)
    lo, hi = Fraction(0), Fraction(1)
    while hi - lo > tol:
        mid = (lo + hi) / 2
        if binom_cdf(k, n, mid) > alpha:
            lo = mid
        else:
            hi = mid
    return lo, hi


def cp_lower(k, n, conf=Fraction(95, 100), tol=Fraction(1, 2 ** 40)):
    """One-sided CP lower bound: largest p with P(X >= k | n, p) <= 1-conf.
    Returns (lo, hi) bracket; quote lo (conservative). k <= 0 -> (0, 0)."""
    alpha = 1 - Fraction(conf)
    if k <= 0:
        return Fraction(0), Fraction(0)
    lo, hi = Fraction(0), Fraction(1)
    while hi - lo > tol:
        mid = (lo + hi) / 2
        if 1 - binom_cdf(k - 1, n, mid) <= alpha:
            lo = mid
        else:
            hi = mid
    return lo, hi


def cp_interval(k, n, conf=Fraction(95, 100), tol=Fraction(1, 2 ** 40)):
    """Two-sided equal-tail CP interval at confidence conf.
    Returns (lower_lo, upper_hi): the conservative OUTER interval."""
    side = 1 - (1 - Fraction(conf)) / 2
    return cp_lower(k, n, side, tol)[0], cp_upper(k, n, side, tol)[1]


def cp_on_fctl(k_on_locus, n_ctl, tol=Fraction(1, 2 ** 40)):
    """Pinned quantity: upper 95% CP bound on F_ctl from the fixed-seed
    N_ctl control sample. Quotes the conservative bracket end."""
    return cp_upper(k_on_locus, n_ctl, Fraction(95, 100), tol)[1]


def exceedance(p, n_bin):
    """Per-bin exceedance P(all N on-locus | p) = p^N, exact."""
    return Fraction(p) ** n_bin


def min_n_for_significance(p_u, alpha0=ALPHA0):
    """Smallest N with p_u^N < alpha0 (power gate: a bin with
    N < this value is INSUFFICIENT-N). p_u >= 1 -> None (never powered).
    Float log gives the starting guess; the answer is verified EXACTLY."""
    p_u, alpha0 = Fraction(p_u), Fraction(alpha0)
    if p_u >= 1:
        return None
    if p_u == 0:
        return 1
    n = max(1, int(log(float(alpha0)) / log(float(p_u))) - 2)
    while p_u ** n >= alpha0:
        n += 1
    while n > 1 and p_u ** (n - 1) < alpha0:
        n -= 1
    return n
