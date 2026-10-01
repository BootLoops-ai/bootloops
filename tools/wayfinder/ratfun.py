#!/usr/bin/env python3
"""
ratfun.py — fixed-eps rational-function layer (mpc-coefficient polynomials).

THE LAYER
=========
The "polynomial / rational-function layer": functions ptrim/padd/pmul/pscale/
peval/pshift/pseries_div and class RF (with the audit fixes kept).
The originating production driver attached a driver-local RatFun A_series
fast path + denominator-root singular points to wayfinder DESystems
(both gates PASS at full stored ~70 digits); this module is the
package promotion of that layer so loaders can offer the same alias-safe
fast path without driver-local monkey-patching (see
de_load.DESystem.enable_fast_path). The driver originals are production
artifacts and are NOT modified or imported at runtime.

DISCIPLINE
==========
* All coefficients are mpc AT FIXED eps — build once per eps value, at an
  explicit working precision; the caller owns workdps (import-dps footgun;
  this module never touches
  global mp.dps — every function must be called inside mp.workdps(...)).
* No sympy. Pure mpmath; sizes here are small (per-entry num/den polys).
"""

from mpmath import mp, mpf, mpc

__all__ = ["ptrim", "padd", "pmul", "pscale", "peval", "pshift",
           "pseries_div", "RF"]


def ptrim(c):
    while len(c) > 1 and c[-1] == 0:
        c.pop()
    return c


def padd(a, b):
    n = max(len(a), len(b))
    out = [mpc(0)] * n
    for i, v in enumerate(a):
        out[i] = out[i] + v
    for i, v in enumerate(b):
        out[i] = out[i] + v
    return ptrim(out)


def pmul(a, b):
    if (len(a) == 1 and a[0] == 0) or (len(b) == 1 and b[0] == 0):
        return [mpc(0)]
    out = [mpc(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x == 0:
            continue
        for j, y in enumerate(b):
            if y != 0:
                out[i + j] += x * y
    return ptrim(out)


def pscale(a, s):
    return ptrim([v * s for v in a])


def peval(a, x):
    r = mpc(0)
    for c in reversed(a):
        r = r * x + c
    return r


def pshift(a, z0):
    """Coefficients of a(z0 + u) in u (repeated synthetic division)."""
    n = len(a)
    out = []
    c = list(a)
    for k in range(n):
        # divide c by (x - z0): remainder = c(z0)
        rem = mpc(0)
        for i in range(len(c) - 1, -1, -1):
            rem = rem * z0 + c[i]
        out.append(rem)
        # deflate: q coefficients
        q = [mpc(0)] * (len(c) - 1)
        acc = mpc(0)
        for i in range(len(c) - 1, 0, -1):
            acc = c[i] + acc * z0
            q[i - 1] = acc
        c = q if q else [mpc(0)]
        if len(out) == n:
            break
    return ptrim(out)


def pseries_div(num, den, M):
    """Taylor series of num(u)/den(u) about u=0 to order M (den[0] != 0)."""
    d0 = den[0]
    if d0 == 0:
        raise ZeroDivisionError("series division: den(0) == 0")
    inv0 = 1 / d0
    out = []
    for k in range(M + 1):
        s = num[k] if k < len(num) else mpc(0)
        top = min(k, len(den) - 1)
        for j in range(1, top + 1):
            dj = den[j]
            if dj != 0:
                s -= dj * out[k - j]
        out.append(s * inv0)
    return out


class RF:
    """Rational function num(x)/den(x), mpc coefficients (fixed eps).
    num/den are NOT
    reduced to lowest terms — common factors survive arithmetic, so
    denominator roots may include APPARENT points (conservative for step
    control; see de_load.enable_fast_path docstring)."""
    __slots__ = ("num", "den")

    def __init__(self, num, den=None):
        self.num = ptrim(list(num))
        self.den = ptrim(list(den)) if den is not None else [mpc(1)]
        if len(self.den) == 1 and self.den[0] == 0:
            raise ZeroDivisionError("RF with zero denominator")

    def is_zero(self):
        return len(self.num) == 1 and self.num[0] == 0

    def __add__(self, o):
        if self.den == o.den:
            return RF(padd(self.num, o.num), self.den)
        return RF(padd(pmul(self.num, o.den), pmul(o.num, self.den)),
                  pmul(self.den, o.den))

    def __neg__(self):
        return RF(pscale(self.num, mpc(-1)), self.den)

    def __sub__(self, o):
        return self + (-o)

    def __mul__(self, o):
        return RF(pmul(self.num, o.num), pmul(self.den, o.den))

    def __truediv__(self, o):
        if o.is_zero():
            raise ZeroDivisionError("RF division by zero RF")
        return RF(pmul(self.num, o.den), pmul(self.den, o.num))

    def __pow__(self, k):
        if k == 0:
            return RF([mpc(1)])
        b = self if k > 0 else RF(self.den, self.num)
        k = abs(k)
        out = RF([mpc(1)])
        while k:
            if k & 1:
                out = out * b
            b = b * b
            k >>= 1
        return out

    def eval(self, x):
        dv = peval(self.den, x)
        if dv == 0:
            raise ZeroDivisionError("RF eval at denominator zero")
        return peval(self.num, x) / dv

    def taylor(self, z0, M):
        """Taylor coefficients about z0 to order M."""
        return pseries_div(pshift(self.num, z0), pshift(self.den, z0), M)

    def series_inf(self, M):
        """Coefficients m_p of sum_{p>=1} m_p x^{-p}; raises if the
        expansion at infinity has any x^{>=0} term (not Fuchsian-ready)."""
        dn, dd = len(self.num) - 1, len(self.den) - 1
        v = dd - dn                       # leading power ~ x^{-v}
        nr = list(reversed(self.num))     # num(1/w) * w^dn
        dr = list(reversed(self.den))
        ser = pseries_div(nr, dr, M + max(0, -v) + 2)
        # entry(x) = w^{v} * ser(w),  w = 1/x
        out = [mpc(0)] * (M + 1)          # out[p] = coeff of x^{-p}, p=0..M
        for k, cv in enumerate(ser):
            p = v + k
            if p < 0 and cv != 0:
                raise ValueError("entry grows at x=inf (power %d)" % (-p))
            if 0 <= p <= M:
                out[p] = cv
        if out[0] != 0:
            raise ValueError("entry has a constant term at x=inf")
        return out

    def laurent0(self, M):
        """(min_order, coeffs) about x=0: entry = sum c_k x^{min_order+k}."""
        vn = 0
        while vn < len(self.num) - 1 and self.num[vn] == 0:
            vn += 1
        if self.num[vn] == 0:
            return (0, [mpc(0)] * (M + 1))
        vd = 0
        while vd < len(self.den) - 1 and self.den[vd] == 0:
            vd += 1
        mo = vn - vd
        ser = pseries_div(self.num[vn:], self.den[vd:], M)
        return (mo, ser)
