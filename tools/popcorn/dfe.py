"""popcorn.dfe — certified DFE-mixing kernel + certified quadrature.

Aliases, in place and by identity (see popcorn._loader), loaded lazily:
  dfe_layer — certified DFE-mixing kernel cache + exact param gradients
              (incl. exact d_s gamma(s,x)); tails bounded analytically.
              Needs python-flint (it drives the ball route per node).
  certquad  — certified composite Gauss-Legendre quadrature with log-range-
              bounded panels + degree-pair self-check (mp.quad silently
              underconverges on exponential-ramp panels — its error estimate
              is absolute; that failure is why this module exists).
              mpmath only.

The kernel's evidentiary controls are the built-in degree-pair self-check,
the analytic tail bounds, and the finite-difference gradient gate in the
battery. For production claims, also gate against an oracle that shares no
code with this layer.
"""
from ._loader import import_in_place

__all__ = ["dfe_layer", "certquad", "CertifiedKernel"]


def __getattr__(name):
    if name in ("dfe_layer", "certquad"):
        mod = import_in_place(name)
        globals()[name] = mod
        return mod
    if name == "CertifiedKernel":
        return __getattr__("dfe_layer").CertifiedKernel
    raise AttributeError(f"module 'popcorn.dfe' has no attribute {name!r}")
