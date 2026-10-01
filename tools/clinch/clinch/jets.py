"""clinch.jets — N-variable second-order jets over python-flint arb balls.

A Jet carries (value, gradient, upper-triangular Hessian) with respect to a
small set of ACTIVE variables (N <= 3 in the v3.1 oracle: survival
(eta, log_nu), advance (eta,), fecundity (eta, a_fec1, log_phi)). Every
component is an arb ball, so a Jet evaluated over a box yields rigorous
enclosures of the function and its first two derivatives over that box.

Order-1 mode: h=None skips all second-order work (the gradient-only oracle
pass); mixing orders raises (typed) — a silent order drop would fake rigor.

LAW (banked in baller MANUAL): never x**n on a ball that may contain zero
(python-flint pow routes through exp*log and NaNs) — powers here are
explicit products.
"""
from flint import arb

__all__ = ["Jet", "JetError", "jconst", "jvar"]

_HPAIRS = {1: ((0, 0),),
           2: ((0, 0), (0, 1), (1, 1)),
           3: ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))}
_HIDX = {n: {p: i for i, p in enumerate(_HPAIRS[n])} for n in _HPAIRS}


class JetError(TypeError):
    """Typed jet misuse (order/arity mix, unsupported operand)."""


class Jet:
    __slots__ = ("v", "g", "h", "n")

    def __init__(self, v, g, h, n):
        self.v, self.g, self.h, self.n = v, g, h, n

    # -- construction -----------------------------------------------------
    @staticmethod
    def const(v, n, order2=True):
        z = arb(0)
        return Jet(v if isinstance(v, arb) else arb(v),
                   [z] * n, ([z] * len(_HPAIRS[n]) if order2 else None), n)

    @staticmethod
    def var(v, i, n, order2=True):
        g = [arb(0)] * n
        g[i] = arb(1)
        return Jet(v if isinstance(v, arb) else arb(v), g,
                   ([arb(0)] * len(_HPAIRS[n]) if order2 else None), n)

    def _like(self, other):
        if not isinstance(other, Jet):
            raise JetError(f"jet op with {type(other).__name__}")
        if other.n != self.n or ((self.h is None) != (other.h is None)):
            raise JetError("jet arity/order mismatch")
        return other

    # -- ring ops ---------------------------------------------------------
    def __neg__(self):
        return Jet(-self.v, [-x for x in self.g],
                   None if self.h is None else [-x for x in self.h], self.n)

    def __add__(self, other):
        if isinstance(other, Jet):
            self._like(other)
            return Jet(self.v + other.v,
                       [a + b for a, b in zip(self.g, other.g)],
                       None if self.h is None else
                       [a + b for a, b in zip(self.h, other.h)], self.n)
        return Jet(self.v + other, list(self.g),
                   None if self.h is None else list(self.h), self.n)

    __radd__ = __add__

    def __sub__(self, other):
        return self.__add__(-other if isinstance(other, Jet) else -arb(other))

    def __rsub__(self, other):
        return (-self).__add__(arb(other))

    def __mul__(self, other):
        if isinstance(other, Jet):
            self._like(other)
            v = self.v * other.v
            g = [self.g[i] * other.v + other.g[i] * self.v
                 for i in range(self.n)]
            if self.h is None:
                return Jet(v, g, None, self.n)
            h = []
            for k, (i, j) in enumerate(_HPAIRS[self.n]):
                h.append(self.h[k] * other.v + other.h[k] * self.v
                         + self.g[i] * other.g[j] + self.g[j] * other.g[i])
            return Jet(v, g, h, self.n)
        c = other if isinstance(other, arb) else arb(other)
        return Jet(self.v * c, [x * c for x in self.g],
                   None if self.h is None else [x * c for x in self.h],
                   self.n)

    __rmul__ = __mul__

    def __truediv__(self, other):
        if isinstance(other, Jet):
            return self.__mul__(other._recip())
        c = other if isinstance(other, arb) else arb(other)
        return self.__mul__(1 / c)

    def _recip(self):
        iv = 1 / self.v
        iv2 = iv * iv
        g = [-x * iv2 for x in self.g]
        if self.h is None:
            return Jet(iv, g, None, self.n)
        iv3_2 = 2 * iv2 * iv
        h = []
        for k, (i, j) in enumerate(_HPAIRS[self.n]):
            h.append(self.g[i] * self.g[j] * iv3_2 - self.h[k] * iv2)
        return Jet(iv, g, h, self.n)

    # -- chain rule for f(u): f', f'' supplied ----------------------------
    def _chain(self, fv, fp, fpp):
        g = [x * fp for x in self.g]
        if self.h is None:
            return Jet(fv, g, None, self.n)
        h = []
        for k, (i, j) in enumerate(_HPAIRS[self.n]):
            h.append(self.h[k] * fp + self.g[i] * self.g[j] * fpp)
        return Jet(fv, g, h, self.n)

    def exp(self):
        e = self.v.exp()
        return self._chain(e, e, e)

    def log(self):
        iv = 1 / self.v
        return self._chain(self.v.log(), iv,
                           None if self.h is None else -iv * iv)

    def log1p(self):
        w = 1 + self.v
        iw = 1 / w
        return self._chain(self.v.log1p(), iw,
                           None if self.h is None else -iw * iw)

    def expm1(self):
        e = self.v.exp()
        return self._chain(self.v.expm1(), e, e)

    def expit(self):
        """1/(1+exp(-x)): p' = p(1-p), p'' = p(1-p)(1-2p)."""
        p = 1 / (1 + (-self.v).exp())
        q = p * (1 - p)
        return self._chain(p, q, None if self.h is None else q * (1 - 2 * p))

    def lgamma(self):
        """log Gamma on a POSITIVE ball; digamma + trigamma via certified
        MONOTONE ENDPOINT HULLS (dig_ball/trig_ball/lgam_ball) — arb's own
        ball digamma degrades fast and NaNs at rad ~ 0.85 (measured,
        python-flint 0.8.0; [10 +/- 0.85].digamma() = nan while the true
        range is [2.158, 2.337]); psi/psi' are monotone on (0, inf) so the
        endpoint hull is rigorous AND tight."""
        v = self.v
        fp = dig_ball(v)
        fpp = None if self.h is None else trig_ball(v)
        return self._chain(lgam_ball(v), fp, fpp)

    def pow_const(self, c):
        """self**c via exp(c*log(self)) — self must be provably > 0."""
        if not (self.v > 0):
            raise JetError("pow_const on a non-positive-provable ball")
        return (self.log() * c).exp()


_NAN = None


def _nan():
    global _NAN
    if _NAN is None:
        _NAN = arb(0) / arb(0)
    return _NAN


def dig_ball(x):
    """Certified digamma of a positive arb ball via the monotone endpoint
    hull (psi strictly increasing on (0, inf)); non-positive -> NaN
    (fail-closed: downstream finiteness checks refuse)."""
    if not (x > 0):
        return _nan()
    return x.lower().digamma().union(x.upper().digamma())


def trig_ball(x):
    """Certified trigamma (Hurwitz zeta(2, .)) of a positive ball via the
    monotone endpoint hull (psi' strictly decreasing on (0, inf))."""
    if not (x > 0):
        return _nan()
    two = arb(2)
    return two.zeta(x.upper()).union(two.zeta(x.lower()))


_XSTAR = None


def lgam_ball(x):
    """Certified log-Gamma of a positive ball: endpoint hull, plus the
    interior minimum at x* = 1.46163... when the ball may straddle it
    (lgamma is convex on (0, inf))."""
    global _XSTAR
    if not (x > 0):
        return _nan()
    lo, hi = x.lower(), x.upper()
    out = lo.lgamma().union(hi.lgamma())
    if _XSTAR is None:
        _XSTAR = arb("1.4616321449683623413", "1e-18")
    if not (hi < _XSTAR or lo > _XSTAR):
        out = out.union(_XSTAR.lgamma())
    return out


def jconst(v, n, order2=True):
    return Jet.const(v, n, order2)


def jvar(v, i, n, order2=True):
    return Jet.var(v, i, n, order2)


def hidx(n, i, j):
    """Index of (i,j) (i<=j) in the packed upper-triangular Hessian."""
    return _HIDX[n][(i, j) if i <= j else (j, i)]
