"""dist_cert.ball -- exact rational ball/interval arithmetic.

Primitives for the certified-distance module dist_cert.
NO census data is touched here. All arithmetic is exact (fractions.Fraction);
sqrt/ln enclosures are certified rational bounds (directed, never floats).
"""
from fractions import Fraction
from math import isqrt

F = Fraction


def sqrt_enclosure(x, bits=64):
    """Certified [lo, hi] Fractions with lo^2 <= x <= hi^2; exact when possible.

    x: nonnegative Fraction. Uses sqrt(p/q) = sqrt(p*q)/q and integer isqrt:
    s = floor(sqrt(p*q*M^2)), M = 2^bits => s/(q*M) <= sqrt(x) <= (s+1)/(q*M).
    """
    x = F(x)
    if x < 0:
        raise ValueError("sqrt_enclosure: negative argument")
    if x == 0:
        return F(0), F(0)
    p, q = x.numerator, x.denominator
    M = 1 << bits
    n = p * q * M * M
    s = isqrt(n)
    if s * s == n:
        v = F(s, q * M)
        return v, v
    return F(s, q * M), F(s + 1, q * M)


def _atanh_enclosure(u, terms=64):
    """Certified [lo, hi] for atanh(u), 0 <= u < 1 Fraction (series, all terms > 0).

    Tail after N terms bounded by u^(2N+1)/((2N+1)(1-u^2)).
    """
    u = F(u)
    assert 0 <= u < 1
    if u == 0:
        return F(0), F(0)
    s = F(0)
    upow = u
    u2 = u * u
    for k in range(terms):
        s += upow / (2 * k + 1)
        upow *= u2
    tail = upow / ((2 * terms + 1) * (1 - u2))
    return s, s + tail


# ln 2 = 2 atanh(1/3); enclosure computed once at import (exact rational bounds)
_LN2_LO, _LN2_HI = (2 * b for b in _atanh_enclosure(F(1, 3), 64))


def ln_enclosure(x, terms=64):
    """Certified [lo, hi] Fractions for ln(x), x > 0 Fraction.

    Argument-reduce by powers of 2 into [2/3, 4/3], then ln(y) = 2 atanh(u),
    u = (y-1)/(y+1) with |u| <= 1/5 -- fast certified series.
    """
    x = F(x)
    if x <= 0:
        raise ValueError("ln_enclosure: nonpositive argument")
    k = 0
    y = x
    while y > F(4, 3):
        y /= 2
        k += 1
    while y < F(2, 3):
        y *= 2
        k -= 1
    u = (y - 1) / (y + 1)
    if u >= 0:
        alo, ahi = _atanh_enclosure(u, terms)
        lo, hi = 2 * alo, 2 * ahi
    else:
        alo, ahi = _atanh_enclosure(-u, terms)
        lo, hi = -2 * ahi, -2 * alo
    if k >= 0:
        return lo + k * _LN2_LO, hi + k * _LN2_HI
    return lo + k * _LN2_HI, hi + k * _LN2_LO


class RIv:
    """Closed rational interval [lo, hi] (exact Fractions)."""

    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None):
        self.lo = F(lo)
        self.hi = F(lo if hi is None else hi)
        if self.lo > self.hi:
            raise ValueError("RIv: lo > hi")

    def __add__(self, o):
        o = o if isinstance(o, RIv) else RIv(o)
        return RIv(self.lo + o.lo, self.hi + o.hi)

    def __sub__(self, o):
        o = o if isinstance(o, RIv) else RIv(o)
        return RIv(self.lo - o.hi, self.hi - o.lo)

    def __mul__(self, o):
        o = o if isinstance(o, RIv) else RIv(o)
        c = (self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi)
        return RIv(min(c), max(c))

    def contains(self, o):
        o = o if isinstance(o, RIv) else RIv(o)
        return self.lo <= o.lo and o.hi <= self.hi

    def overlaps(self, o):
        return self.lo <= o.hi and o.lo <= self.hi

    def __repr__(self):
        return f"RIv[{self.lo}, {self.hi}]"


def sqrt_iv(iv, bits=64):
    """Certified sqrt of a nonnegative RIv."""
    lo, _ = sqrt_enclosure(max(F(0), iv.lo), bits)
    _, hi = sqrt_enclosure(iv.hi, bits)
    return RIv(lo, hi)


class CBall:
    """Complex ball: exact rational center (re, im) + exact rational radius."""

    __slots__ = ("re", "im", "rad")

    def __init__(self, re, im=0, rad=0):
        self.re, self.im, self.rad = F(re), F(im), F(rad)
        if self.rad < 0:
            raise ValueError("CBall: negative radius")

    def __sub__(self, o):
        return CBall(self.re - o.re, self.im - o.im, self.rad + o.rad)

    def abs_iv(self, bits=64):
        """Certified RIv enclosing |z| for every z in the ball. Exact at 0."""
        m2 = self.re * self.re + self.im * self.im
        clo, chi = sqrt_enclosure(m2, bits)
        return RIv(max(F(0), clo - self.rad), chi + self.rad)

    def __repr__(self):
        return f"CBall({self.re}, {self.im}; r={self.rad})"


def dist_iv(a, b, bits=64):
    """Certified RIv enclosing |z_a - z_b| over the two balls. Exact 0 when
    both are exact points (rad 0) with identical centers."""
    return (a - b).abs_iv(bits)
