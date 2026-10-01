"""BALLER — arbitrary precision ball arithmetic for science (front door).

Ball-arithmetic instruments ONLY — never a header for Arb;
eras (Taylor forms) is a SEPARATE consumer package.

FRONT DOOR:
BALLER never aliases/shadows Arb's namespace, but it OWNS the
user-facing ball-arithmetic layer that builds with it:
  baller.run(f, x0=None, n=None, dps=50)   expression/iteration in ball
                    arithmetic at working precision (thin over arb/acb)
  baller.solve(f, target_digits, ...)      adaptive certified driver —
                    radius-watching precision escalation; SolveRefused typed
  baller.render(x, digits=None, ...)       fail-closed printing — certified
                    digits only; uncertified marked/refused, never bare
(baller.frontdoor module docstring carries the full contract INCLUDING the
documented zero-divisor DEVIATION from Arb's nan convention: on a divisor
ball containing 0 the radius goes to +inf and the central value continues
as the fixed-precision stream — a valid enclosure.)

Modules:
  baller.quad       certified nested quadrature — the gated POSQ core,
                    ALIASED in place (its gates are passed; tools/posq stays
                    the tool of record)
  baller.transport  certified period/ODE transport — vendored kklt engines
                    (pipe_transport, pipe_lib, march_lib with the gated
                    complex-step fix)
  baller.certify    Krawczyk / interval-Newton unique-zero machinery —
                    vendored kklt pipe_vac engine half
  baller.contract   dual-path tripwires that halt rather than average —
                    wayfinder gate+manifest aliased in place, geo mc
                    (batched MC + Clopper-Pearson certified quantiles)
                    vendored; the tripwire + mutation harness are native
  baller.hygiene    ambient-dps lints + ball footgun checks — the
                    dps_lint member (mpmath import-time-dps AST linter,
                    `python3 -m baller.hygiene lint ...`), the purity_scan
                    member, plus native enforced checks

Integrity: baller.verify() — fail-closed on the vendored set, advisory-loud
on aliased-source drift; vendored engines load via exec(compile(source)),
never the bytecode cache.
"""
from ._core import verify, VendorTamperError, load_vendored, alias_registered

__version__ = "1.0.0"
_SUBMODULES = ("quad", "transport", "certify", "contract", "hygiene",
               "frontdoor")
_FRONTDOOR = ("run", "solve", "render", "certified_digits", "RunResult",
              "SolveResult", "SolveRefused", "RenderRefused")


def __getattr__(name):
    if name in _SUBMODULES:
        import importlib
        mod = importlib.import_module(f".{name}", __name__)
        globals()[name] = mod
        return mod
    if name in _FRONTDOOR:
        from . import frontdoor
        obj = getattr(frontdoor, name)
        globals()[name] = obj
        return obj
    raise AttributeError(f"module 'baller' has no attribute {name!r}")


__all__ = (list(_SUBMODULES) + list(_FRONTDOOR)
           + ["verify", "VendorTamperError", "__version__"])
