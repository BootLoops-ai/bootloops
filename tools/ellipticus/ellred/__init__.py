"""ellipticus.ellred — ellred member: standard-elliptic (Legendre/Carlson)
reduction of one-fold elliptic integrals  int rat(x) dx / sqrt(cubic|quartic).

Every reduction is gated against direct high-precision quadrature of the
defining integrals (measured classes on this copy: 45-51d quartic atoms,
30-38d cubic moments, ~30d 1r2c modulus, ~50-60d route agreements).

Modules:
  ell_normal    — vendored Legendre chart layer (quartic_pairs/legendre_chart/
                  I0_legendre/G_pval + the P4_exact demonstration family)
  ellred_engine — quartic E: chart, F/E/Pi atoms (closed + calibrated),
                  J-moment recurrence, pole atoms
  ellred_cubic  — cubic Legendre moments: 3-real-root (both regions) and
                  1-real-2-complex charts, Jmoments_cubic
  ellred_reduce — x-space reducer: int rat(x) dx/sqrt(P4) -> F/E/Pi + algebraic
  ellred_1r2c   — analytic 1r2c reducer: F/E/H moments + 3rd-kind pole atoms,
                  exact analytic-derivative closures (certified 1e-50..1e-103)
  ellred_carlson— Carlson R_F route cross-check of the 1r2c reduction
  ellred_kernels— consolidated re-exports + P4_at family fixture
  derive_atoms  — closed-form E/Pi atom derivation gate (script, run via -m)

Battery: leg L5 in ellipticus.battery exercises these code paths against
direct quadrature + the pinned reference strings in
fixtures/battery/ellred_controls.json.
"""
from . import ell_normal
from . import ellred_engine
from . import ellred_cubic
from . import ellred_1r2c
from . import ellred_carlson
from . import ellred_reduce
from . import ellred_kernels

__all__ = ["ell_normal", "ellred_engine", "ellred_cubic", "ellred_1r2c",
           "ellred_carlson", "ellred_reduce", "ellred_kernels"]
