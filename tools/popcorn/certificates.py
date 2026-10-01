"""popcorn.certificates — the exact certificate-instrument family.

Aliased in place and by identity (see popcorn._loader), loaded lazily so
`import popcorn.certificates` stays dependency-light (scipy is required only
when cert_lp is touched; flint only for positivity route R and cone_lp_b):

  positivity  — exact nonnegativity of sparse rational polynomials on [0,1],
              two routes: B = rational Bernstein/de Casteljau subdivision
              with touch-root deflation; R = certified flint root isolation +
              exact sign sampling. Pure rational.
  exact_lp    — pure-rational LP membership + Farkas witness (phase-1
              simplex, Bland); the slow validation oracle.
  fast_lp     — float may PROPOSE, only exact rational arithmetic DECIDES;
              shift-trick Farkas certificates.
  cone_lp_b   — exact integer-cone dual simplex (one exact bit-sorted RREF
              basis completion — the fast host variant; a greedy per-column
              rank-check variant of the same simplex measures ~40x slower at
              n=50, so do not "simplify" the basis completion back).
  cert_lp     — CertHull certify-after-float membership head with
              certificate reuse; scipy-HiGHS is PROPOSER-ONLY.
  region      — exact-rational REGION certificates
              (region_certificates.py): Bernstein box positivity + dyadic
              subdivision, certified_sup with LP-dual rationalization,
              Cauchy-Schwarz certified chi^2 lower bound, rational-function
              reconstruction with held-out degree certification + Sturm root
              counting. Dependency-light (stdlib; sympy for root counting).

The decide() law for all membership instruments: verdicts are
FEASIBLE/INFEASIBLE (or a named refusal) with independently re-substituted
exact certificates — never a rounded verdict, never an averaged one.
"""
from ._loader import import_in_place

_MODULES = ("positivity", "exact_lp", "fast_lp", "cone_lp_b", "cert_lp")

__all__ = list(_MODULES) + ["region", "CertHull", "membership_exact",
                            "membership_fast"]


def __getattr__(name):
    if name == "region":
        from . import region_certificates as mod
        globals()["region"] = mod
        return mod
    if name in _MODULES:
        mod = import_in_place(name)
        globals()[name] = mod
        return mod
    if name == "CertHull":                      # the production instrument
        return __getattr__("cert_lp").CertHull
    if name == "membership_exact":              # exact cone decider (host)
        return __getattr__("cone_lp_b").membership_exact
    if name == "membership_fast":               # float-propose/exact-decide
        return __getattr__("fast_lp").membership_fast
    raise AttributeError(f"module 'popcorn.certificates' has no attribute {name!r}")
