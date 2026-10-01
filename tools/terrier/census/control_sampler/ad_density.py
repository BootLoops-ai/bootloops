"""Ashok-Douglas index density factor (the H0-a index-density control null).

Pinned formula [AD hep-th/0307049 eq (4.10), harvmac label dmuIresult]:
  d mu_I(z) = pi^-n det(-R + kappa omega 1),  <W W*> = e^(kappa K),
with KAPPA = -1 for the flux ensemble  =>  integrand det(-R - omega 1).

This module evaluates det(-R(z) - omega(z)) EXACTLY (Fraction arithmetic)
given exact-rational matrix callables R(z), omega(z) in the pinned chart.
The pi^-n prefactor and the coordinate measure drop out of rejection
sampling and of F_ctl (both consume only density RATIOS / count ratios),
so no float constant is ever introduced.

SIGN LAW: eq (4.10) is an INDEX density; the sampler density must be
nonnegative on the sampling domain. A negative determinant is
therefore a DATA ERROR here -- ADDensity raises, never silently |.|'s.
"""
from fractions import Fraction

KAPPA = -1  # pinned: <W W*> = e^(-K) for the flux ensemble


def det_fraction(M):
    """Exact determinant, fraction-preserving Gaussian elimination."""
    n = len(M)
    A = [[Fraction(M[i][j]) for j in range(n)] for i in range(n)]
    det = Fraction(1)
    for c in range(n):
        piv = next((r for r in range(c, n) if A[r][c] != 0), None)
        if piv is None:
            return Fraction(0)
        if piv != c:
            A[c], A[piv] = A[piv], A[c]
            det = -det
        det *= A[c][c]
        inv = 1 / A[c][c]
        for r in range(c + 1, n):
            f = A[r][c] * inv
            if f:
                for j in range(c, n):
                    A[r][j] -= f * A[c][j]
    return det


class ADDensity:
    """det(-R(z) - omega(z)) with R, omega exact n x n matrix callables."""

    def __init__(self, n, R_fn, omega_fn):
        self.n, self.R_fn, self.omega_fn = n, R_fn, omega_fn

    def __call__(self, z):
        R, w = self.R_fn(z), self.omega_fn(z)
        A = [[-Fraction(R[i][j]) - Fraction(w[i][j])
              for j in range(self.n)] for i in range(self.n)]
        d = det_fraction(A)
        if d < 0:
            raise ValueError("AD index density negative at %r (sign law)" % (z,))
        return d
