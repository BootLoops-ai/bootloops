"""CLINCH — Certified Local Interval-Newton Convergence for Hierarchies.

Given a fitted MAP optimum theta* and model oracles evaluable in ball
arithmetic, CLINCH either issues a CERTIFICATE — inside B(theta*, r) a
true stationary point of the penalized objective EXISTS, is UNIQUE, and
the Hessian is positive definite there (a certified local minimum) — or
REFUSES, honestly, with the failed contraction named. No third outcome.

The generic
certification kernel is baller.certify.block_krawczyk (consumed by
IDENTITY import, never copied). Modules:
  clinch.jets        N-var second-order jets over arb balls
  clinch.oracle_v31  the v3.1 model oracle (gradient + block-arrow Hessian)
  clinch.adapter_v31 reference-data bridge + oracle identity gates
  clinch.engine      radius policy + certificate driver (code-pinned, no knobs)
  clinch.cert        CERT.json + human statement

Integrity: clinch.verify() — fail-closed source-sha check of every module
against _pins.PINS; bytecode caches are purged at import and never
written (sys.dont_write_bytecode). Boundary, stated plainly: an
in-process actor who can rewrite _pins.py or __init__.py can also rewrite
verify itself — no in-process check beats that class (the same
boundary BALLER documents).
"""
import hashlib
import os
import shutil
import sys

sys.dont_write_bytecode = True
_HERE = os.path.dirname(os.path.abspath(__file__))
_pyc = os.path.join(_HERE, "__pycache__")
if os.path.isdir(_pyc):
    shutil.rmtree(_pyc, ignore_errors=True)

__version__ = "1.0.0"


class ClinchTamperError(RuntimeError):
    """A pinned clinch module's source bytes do not match _pins.PINS."""


def verify(quiet=False):
    """Fail-closed source integrity: sha256 of every pinned module."""
    from . import _pins
    bad = []
    for name, want in _pins.PINS.items():
        p = os.path.join(_HERE, name)
        try:
            got = hashlib.sha256(open(p, "rb").read()).hexdigest()[:16]
        except OSError as e:
            raise ClinchTamperError(f"{name}: unreadable ({e})")
        if got != want:
            bad.append((name, want, got))
    if bad:
        raise ClinchTamperError(
            "TAMPER: " + "; ".join(f"{n} pinned {w} got {g}"
                                   for n, w, g in bad))
    if not quiet:
        print(f"clinch.verify: {len(_pins.PINS)} modules pinned OK")
    return {"pinned_ok": sorted(_pins.PINS)}


_SUBMODULES = ("jets", "oracle_v31", "adapter_v31", "engine", "cert")


def __getattr__(name):
    if name in _SUBMODULES:
        import importlib
        verify(quiet=True)
        mod = importlib.import_module(f".{name}", __name__)
        globals()[name] = mod
        return mod
    raise AttributeError(f"module 'clinch' has no attribute {name!r}")


__all__ = list(_SUBMODULES) + ["verify", "ClinchTamperError", "__version__"]
