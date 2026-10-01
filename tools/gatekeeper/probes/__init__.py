"""
gatekeeper.probes — Gatekeeper's probe kit: pre-farm quadrature honesty.

  ell_subst  — weight-matched 1D substitution selector for sqrt-product
               chains (Gauss-Chebyshev / Jacobi-sn via cross-ratio /
               one-sided Jacobi-sn^2; honest 'unsupported' instead of the
               silent tanh-sinh ~3 d/level trap).
  quad_probe — per-chain quadrature convergence prober (vary ONE chain's
               level, others cancel exactly in consecutive diffs; <5
               digits/level = METHOD-MISMATCHED).  Probe BEFORE farming.

The flat tools/ell_subst.py + tools/quad_probe.py paths are
attribute-complete shims that forward here.
NOTE: ell_subst needs sympy on top of the harness's stdlib+mpmath
baseline; quad_probe is stdlib+mpmath only.
"""
__all__ = ["select", "cubic_under_sqrt_K", "SubstPlan",
           "probe", "EnvEvaluator"]


def __getattr__(name):
    # lazy imports, same pattern as gatekeeper/__init__.py — keeps
    # `python -m gatekeeper.probes.quad_probe` from double-importing
    if name in ("select", "cubic_under_sqrt_K", "SubstPlan"):
        from . import ell_subst as _e
        return getattr(_e, name)
    if name in ("probe", "EnvEvaluator"):
        from . import quad_probe as _q
        return getattr(_q, name)
    raise AttributeError(name)
