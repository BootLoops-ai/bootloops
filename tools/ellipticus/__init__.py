"""Ellipticus — certified evaluator for elliptic multiple polylogarithms and
iterated integrals on an elliptic curve (BootLoops toolkit package; import
name `ellipticus`, formerly `empl_eval`).

Given an exact curve/family spec (quartic or cubic y^2 = F(x;z) with
coefficients rational in z), letters, a word, an evaluation point and a
precision, every verb returns values with a certificate.

Three certified engines behind one front door (see SUMMARY.md + GUIDE.md
for word-class coverage, measured rates, caveats):

  engine i  (march, the workhorse): exact inhomogeneous-PF extended systems
            (extend.moment_system) + certified adaptive-Taylor
            transport with dlog-letter states, tangential-base regularization
            and complex-waypoint detours (march.System).
  engine q  (fast tau-words): vendored Kronecker-Eisenstein kernels + AGW
            shuffle-regularized tau-iterated integrals + BMSW sunrise closed
            forms with PROVEN q-tail bounds (qengine).
  engine agm (weight-1/period layer): exact curve layer, AGM/Legendre-K
            periods, sin^2-substitution certified cycle moments (curve).

PURITY BAR: exact inputs only (Fraction/int/str); floats are refused at the
API boundary.  Every verb takes dps explicitly; the two-precision self-gate
is built into evaluate().

ENVIRONMENT: every ELLIPTICUS_* variable (ELLIPTICUS_CACHE, ELLIPTICUS_W3,
ELLIPTICUS_ARTIFACTS, ELLIPTICUS_KRON_DIR, ELLIPTICUS_ROWS_DIR,
ELLIPTICUS_BATTERY_OUT) is also read under its former EMPL_EVAL_* /
EMPL_BATTERY_OUT name; the new name wins when both are set (see env()).
"""
import os as _os


def env(new, old, default=None):
    """Read an environment variable under its current name `new`, falling
    back to the former name `old` (kept working for existing callers), then
    to `default`.  Empty strings count as unset."""
    return _os.environ.get(new) or _os.environ.get(old) or default


from .curve import QuarticCurve, sunrise_frame, exact
from .extend import moment_system, polyform, analytic_seed_series, h_series
from .march import System, letter_base_log, letter_base_r
from .api import MomentFamily, evaluate, shared_digits
from . import qengine

__version__ = "1.0.0"

__all__ = [
    "env",
    "QuarticCurve", "sunrise_frame", "exact",
    "moment_system", "polyform", "analytic_seed_series", "h_series",
    "System", "letter_base_log", "letter_base_r",
    "MomentFamily", "evaluate", "shared_digits",
    "qengine",
]
