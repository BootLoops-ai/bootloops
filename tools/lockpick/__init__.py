"""Lockpick — the house integer-relation / lattice value-fit family.

NOTHING deleted: every old flat name keeps working through an
attribute-complete shim at tools/<old_name>.py (WINNOW-rename pattern).

Submodules (original basenames kept):

  lockpick.mplll         — multi-point LLL/BKZ analytic-regression fitter
                           (v2 core: qrbkz/rawlll/rawbkz; held-out >=30d gate,
                           two-precision, capacity, +/- controls).
  lockpick.mplll_batch   — multi-target LQ cache, CVP/Babai snap, d_min sweep.
  lockpick.mplll_cli     — CLI front end (qrbkz/rawlll/rawbkz/graded/cvp,
                           batch, --dmin-sweep).
  lockpick.mplll_graded  — weight-graded block reduction (--method graded);
                           localises span-deficiency per weight.
  lockpick.mplll_lattice — integer-lattice LLL/BKZ backends (fpylll ->
                           fplll-cli -> Nemo.lll -> pure fallback).
  lockpick.mplll_points  — D-optimal/Leja sample-point selector (basis-side,
                           BEFORE target eval).
  lockpick.pslq_gate     — shared PSLQ closure harness: canonicalize,
                           two_prec_stable, capacity, members_audit,
                           DegeneratePoolError, controls, reverify; --selftest
                           refinds planted synthetic relations it builds
                           itself (positive controls) + member-integrity legs.
  lockpick.ring15        — ring-basis MEMBER (subpackage):
                           conductor-15 / Q(sqrt(-15)) Chowla-Selberg CM
                           constant ring (core/hecke/lfun/orchestrator).
                           Eichler.jl dirichlet_L stays Eichler's by
                           identity — linked, never copied.

Script-member directories (run as scripts, not importable subpackages):

  sealrun/               — MEMBER: preregistered sealed-closure runner on the
                           pslq_gate engine (t6_close.py + t6_basket.py +
                           frozen BASKET_T6.json + committed selftest
                           targets); CLOSED-with-name vs
                           NULL-with-exclusion-floor. See GUIDE.md.
  bilmine/               — MEMBER: prereg-sealed bilinear-relation lattice
                           miner on certified period frames (bilmine_lib +
                           drivers; f1/ = the Z[sqrt15] ring-window stage);
                           synthetic selftest.py rides the member. Data
                           contract via BILMINE_* env vars, no defaults
                           shipped. See GUIDE.md.

Ambient-dps law (tests/test_mplll_ambient_dps.py):
all lattice-build paths own their workdps(d+35) — API calls are safe at any
ambient mp.dps.

Submodule access is lazy: `import lockpick` stays light; `lockpick.mplll`
imports mpmath only on first touch.
"""

__all__ = ["mplll", "mplll_batch", "mplll_cli", "mplll_graded",
           "mplll_lattice", "mplll_points", "pslq_gate", "ring15"]


def __getattr__(name):
    if name in __all__:
        import importlib
        return importlib.import_module("." + name, __name__)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(list(globals()) + list(__all__))
