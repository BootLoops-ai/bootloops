"""Truncated Laurent series in eps over mpmath (complex-capable).

Generalized from the L class of the single-constant c_M prototype.
Differences: complex coefficients allowed, NMAX is an instance-carried cap
(class default kept for drop-in compatibility), and taylor() precision
headroom scales with the ambient dps (no fixed precision guards; see
the dps_lint member of tools/baller, `baller.hygiene.dps_lint`).
"""
import mpmath as mp


class L:
    """sum_{n>=nmin} c[n-nmin] eps^n, truncated above NMAX."""
    NMAX = 6

    def __init__(self, nmin, coeffs, nmax=None):
        self.nmin = nmin
        self.c = [mp.mpmathify(x) for x in coeffs]
        self.nmax = L.NMAX if nmax is None else nmax

    @staticmethod
    def taylor(f, N):
        # headroom scales with working dps (was fixed 4*DPS in the prototype;
        # kept the same law, derived from ambient dps at call time)
        with mp.workdps(4 * mp.mp.dps):
            c = mp.taylor(f, 0, N)
        return L(0, c)

    def __mul__(self, other):
        if not isinstance(other, L):
            return L(self.nmin, [other * x for x in self.c], self.nmax)
        nmin = self.nmin + other.nmin
        out = [mp.mpf(0)] * (len(self.c) + len(other.c) - 1)
        for i, ci in enumerate(self.c):
            for j, cj in enumerate(other.c):
                out[i + j] += ci * cj
        nmax = min(self.nmax, other.nmax)
        keep = nmax - nmin + 1
        return L(nmin, out[:max(keep, 0)], nmax)

    __rmul__ = __mul__

    def __add__(self, other):
        nmin = min(self.nmin, other.nmin)
        nmax_len = max(self.nmin + len(self.c), other.nmin + len(other.c))
        out = [mp.mpf(0)] * (nmax_len - nmin)
        for i, ci in enumerate(self.c):
            out[self.nmin - nmin + i] += ci
        for i, ci in enumerate(other.c):
            out[other.nmin - nmin + i] += ci
        return L(nmin, out, min(self.nmax, other.nmax))

    def coeff(self, n):
        i = n - self.nmin
        return self.c[i] if 0 <= i < len(self.c) else mp.mpf(0)

    def __repr__(self):
        return " + ".join(
            f"({mp.nstr(c, 18)})*eps^{self.nmin + i}" for i, c in enumerate(self.c)
        )
