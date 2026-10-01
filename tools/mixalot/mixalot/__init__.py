"""MIXALOT — exact Bayesian evidence for mixture models, in a box.

A set of validated mixture-evidence engines, vendored byte-exact and
sha-pinned, behind one front door. A user with a
mixture question — how many components? at large k, large g, large N —
computes with the original validated instruments instead of rebuilding them
from a paper's description.

LAW: mixalot.verify() re-hashes every vendored engine against mixalot._pins
FAIL-CLOSED (tampered vendor = tool refuses to compute) and checks the
upstream sources ADVISORY-LOUD (source drift is reported, never bricks the
tool). Vendor bytes are never edited; wrappers parameterize cwd/output paths.

Front door:
    mixalot.verify()                    -> integrity report (raises on tamper)
    mixalot.Z(U, g)                     -> exact evidence, small/medium regime
    mixalot.gstar(U, gmax, prior)       -> exact component-count posterior
    mixalot.bayes_factor(U, g1, g2)     -> exact BF
    mixalot.engines.*                   -> the vendored engines, importable
    mixalot.advise(N=, g=)              -> exact-frontier advisories (quoted)
    mixalot.plant / recover_blind       -> synthetic-truth discipline
Scale routes (large k / large g / large N) are fronted via mixalot.engines
with worked recipes in MANUAL.md; capability battery in battery/.
"""
__version__ = "1.0.0"

from ._verify import verify                    # noqa: F401
from .core import Z, gstar, bayes_factor       # noqa: F401
from .advisory import advise                   # noqa: F401
from .plants import plant, recover_blind       # noqa: F401
from . import engines                          # noqa: F401
