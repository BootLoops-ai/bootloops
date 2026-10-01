"""formglue — glue around the external FORM engine + independent MZV sums.

Glue code that CALLS external FORM (Vermaseren et al. / form-dev; binary at
FORM5_BIN env var, default: `form` on PATH — needs FORM 5 with float support).
The FORM binary, its invocation paths, and its name are untouched — external
tools are called, never absorbed.

Submodules:

  formglue.form_io     — shared line-unwrap fix for FORM output (#write /
                         Print continuation-line footgun).
  formglue.form_oracle — strict wrapper around FORM 5's arbitrary-precision
                         MZV engine: zeta / mzv / euler, rc-coded, GUARD=10.
  formglue.form_route  — momentum-routing glue for FORM 5 diagrams_ output
                         (spanning-tree loop assignment + 5 rc-coded checks).
  formglue.tornheim    — fast Mordell-Tornheim-Witten sums + depth-2 MZVs
                         (mpmath engine, independent of FORM; an
                         independent cross-check oracle).

Submodule access is lazy: `import formglue` stays light; `formglue.tornheim`
imports mpmath only on first touch.

Tests: tools/formglue/tests/test_form_{io,oracle,route}.py (strict, rc-coded).
"""

__all__ = ["form_io", "form_oracle", "form_route", "tornheim"]


def __getattr__(name):
    if name in __all__:
        import importlib
        return importlib.import_module("." + name, __name__)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(list(globals()) + list(__all__))
