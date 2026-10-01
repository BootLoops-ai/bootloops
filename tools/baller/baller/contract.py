"""baller.contract — dual-path tripwires that HALT rather than average.

Aliased in place (registered tool wayfinder stays of record):
  gate      two-precision digit gates against vendored reference strings —
            the reusable halt-not-average discipline
  manifest  hash-pinned input manifests (the reference-table pin machinery)

Vendored (src-sha header):
  mc        generic batched MC + Clopper-Pearson certified quantiles
            (vendor/geo)

NATIVE:
  tripwire  — dual-path halt-not-average API (dual(), Tripwire sampled
              oracle, DualPathDisagreement) — baller.tripwire
  mutation  — mutation-tested-gate harness (scratch-copies-only law) —
              baller.mutation
"""
from .tripwire import dual, Tripwire, DualPathDisagreement, agree_digits
from .mutation import run_mutation_gate, MutationHarnessError
import importlib
import os
import sys

from ._core import load_vendored, TOOLS

_DT = os.path.join(TOOLS, "wayfinder")

# wayfinder is a package: import via its parent (tools/) so submodule
# resolution uses the package machinery, then identity-check the file home.
# the insert is REMOVED after the import — a permanent sys.path[0] entry
# would shadow user modules like their own eta.py process-wide; identity
# survives via sys.modules.
_added = TOOLS not in sys.path
if _added:
    sys.path.insert(0, TOOLS)
try:
    _dt = importlib.import_module("wayfinder")
    gate = importlib.import_module("wayfinder.gate")
    manifest = importlib.import_module("wayfinder.manifest")
finally:
    if _added and TOOLS in sys.path:
        sys.path.remove(TOOLS)
for _m in (gate, manifest):
    _home = os.path.dirname(os.path.abspath(_m.__file__))
    if _home != _DT:
        raise ImportError(f"baller identity violation: {_m.__name__} resolved "
                          f"to {_home}, not {_DT}")


def __getattr__(name):
    if name == "mc":
        mod = load_vendored("mc")
        globals()[name] = mod
        return mod
    raise AttributeError(f"module 'baller.contract' has no attribute {name!r}")


__all__ = ["gate", "manifest", "mc", "dual", "Tripwire",
           "DualPathDisagreement", "agree_digits", "run_mutation_gate",
           "MutationHarnessError"]
