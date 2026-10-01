"""TRUST — streamed-IBP certification triad, assembly package.

The three checks on one reduction table (CORE lineages disjoint — censused
by the battery; see the check-2 as-executed caveat below):

  1. kira crosscheck      — strata's full-table + dictionary held-out check
                            (lineage: Kira 3.1 C++ + Fermat, adopted 3rd party)
  2. lambda-witness audit — receipt core vs the CALLER'S OWN re-parse
                            (CORE lineage: stdlib-only pure-python ints;
                            CAVEAT: the as-executed row re-parse rides
                            strata's loader/CoeffEvaluator — common substrate
                            with check 1's kira-side eval; named caveat +
                            battery-side independent re-parse spot-check in
                            trust/receipt.py docstring)
  3. lp_syz syzygy oracle — parametric exact IBP at rational (d0,eta0)
                            (lineage: Singular syz + flint fmpq RREF;
                            FLINT-substrate caveat -> fraction_oracle leg)

tools/trust/ IS lp_syz's home: the certifier set lives in vendor/lpsyz;
strata/ and receipt/ beside this package are its two member trees.
strata / winnow(ibplapper) / receipt are linked BY IDENTITY, never
copied — see trust/strata.py (incl. the ONE-LINEAGE strata/winnow caveat and
the arbitration table) and trust/receipt.py.

Certifier scope: production-scale lp_syz claims are OUT — the 8-var
top-node syzygy wall stands (>1200 s, measured).

Usage:
    import trust
    trust.verify()                       # fail-closed pins gate
    m = trust.load_vendored("lp_syz_431")  # pinned engine, in-process
    trust.run_vendored("lp_syz_431", ["--stage","control452","--D0","6"])
                                          # child mode; needs TRUST_OUT_ROOT
    trust.strata.load(); trust.receipt.core()   # identity fronts
"""
import os as _os
import sys as _sys


class EnvPoisonError(ImportError):
    """Typed refusal: the importing interpreter's environment can shadow the
    stdlib, so every sha pin this package checks is meaningless."""


class StdlibShadowError(ImportError):
    """Typed refusal: a security-critical stdlib module resolved from OUTSIDE
    the STDLIB roots — the sys.path[0] shadow class (script dir in script
    mode; '' = cwd under -c, the REPL and stdin) needs NO poisoned env var,
    so the env gate cannot see it; module identity can. The roots are
    sysconfig stdlib/platstdlib ONLY — a critical stdlib module resolving
    from ANY site-packages (system or user) is never legitimate, so a
    .pth-reordered site plant refuses here too."""


# --- env perimeter at IMPORT time ------------------------------------------
# A PYTHONPATH/PYTHONSTARTUP/PYTHONHOME world can shadow the stdlib before
# any gate runs: a pins-aware hashlib shadow on PYTHONPATH greens a
# byte-tampered vendor through trust.verify() and executes the tampered
# bytes via load_vendored. This guard covers the package's own in-process
# code paths (verify / load_vendored / oracle_k2disp drivers), not just child
# spawns. By trust-import time the host process has long since imported
# hashlib, so scrubbing cannot help: the only sound response is refusal,
# applied here.
# Exemption: under -E/-I (sys.flags.ignore_environment) the interpreter
# ignored these vars — the environment did NOT poison it (battery children).
# MEMBERSHIP test, not truthiness — a set-but-EMPTY PYTHONPATH happens to be
# path-inert on CPython 3.12, but that inertness is undocumented interpreter
# behavior; the conservative form refuses it too.
_BAD_ENV = [k for k in ("PYTHONPATH", "PYTHONSTARTUP", "PYTHONHOME")
            if k in _os.environ]
if _BAD_ENV and not _sys.flags.ignore_environment:
    raise EnvPoisonError(
        f"trust REFUSES to import: {_BAD_ENV} set in the environment — the "
        f"stdlib-shadow class defeats sha pinning (set-but-empty also "
        f"refuses). Unset and rerun from a clean environment, or launch "
        f"the interpreter -E.")
del _BAD_ENV

# --- post-import stdlib IDENTITY at IMPORT time -----------------------------
# The env gate above closes only the env-var vector. CPython puts sys.path[0]
# AHEAD of the stdlib in every non--P/-I invocation, so a pins-aware
# hashlib.py planted in the directory the operator runs from (or the
# directory a driver script lives in) shadows the stdlib with ZERO poisoned
# env vars, greens a byte-tampered vendor, and executes its bytes through
# exactly this hole. Poison VECTORS are not enumerable; module IDENTITY is:
# every security-critical stdlib module this package rides must have
# resolved from the STDLIB roots.
# A whitelist of ALL interpreter-owned roots is not enough — site.py
# executes 'import' lines in .pth files at interpreter startup, BEFORE any
# user code, so a 2-line .pth in any site dir reorders sys.path and a
# site-packages-planted hashlib WINS import with __file__ INSIDE such a
# whitelist. site-packages is NOT a stdlib home: authenticate against
# sysconfig stdlib/platstdlib (+ the stdlib zip fallback) ONLY — a critical
# stdlib module resolving from ANY site dir, system or user, is never
# legitimate.
# Builtin/frozen modules (no __file__) are interpreter-owned by construction.
# KNOWN LIMITS (named, not silently claimed away): (a) __file__ is spoofable
# by an already-hostile interpreter — this gate authenticates the HONEST-
# import case; (b) .pth code in the HOST interpreter's site dirs runs before
# this gate and is arbitrary — on a typical host the system site dirs are
# root-owned (root-write residual) and the user site is same-user-write
# (an attacker already writing as the user is outside the pin perimeter).
# Child mode (run_vendored: -E -P -s -B) never site-processes the
# user-writable site at all; its residual is the root-owned system .pth.
_STDLIB_CRITICAL = (
    "hashlib", "json", "os", "os.path", "shutil", "subprocess", "threading",
    "types", "importlib", "importlib.util", "ast", "re", "time", "itertools",
    "fractions", "site",
)


def _stdlib_identity_gate():
    import importlib as _il
    import sysconfig as _sc
    # BASE-installation scheme vars: a venv's own lib dir is NOT a stdlib
    # home (venvs share the base installation's stdlib; the venv lib dir
    # legitimately contains ONLY site-packages — yet the default 'venv'
    # scheme reports platstdlib INSIDE the venv, which would prefix-cover
    # the venv's site-packages: exactly the hole the BASE-scheme roots
    # close).
    _scheme = "posix_prefix" if _os.name == "posix" else "nt"
    _bx = getattr(_sys, "base_exec_prefix", _sys.base_prefix)
    _bv = {"base": _sys.base_prefix, "platbase": _bx,
           "installed_base": _sys.base_prefix, "installed_platbase": _bx}
    roots = set()
    for _key in ("stdlib", "platstdlib"):
        try:
            p = _sc.get_path(_key, _scheme, vars=_bv)   # GUARD STDLIB-ROOTS (battery mutation anchor)
        except (KeyError, ValueError):
            p = None
        if p:
            roots.add(_os.path.realpath(p))
    if not roots:
        raise StdlibShadowError(
            "trust REFUSES to import: sysconfig reports no stdlib root — "
            "module identity is unverifiable (fail closed)")
    # lib-dynload lives under platstdlib (prefix match covers it); the stdlib
    # zip fallback (pythonXY.zip beside the stdlib dir) is named explicitly —
    # a zipimported module's __file__ keeps the .zip path segment
    for r in list(roots):
        roots.add(_os.path.join(_os.path.dirname(r),
                                "python%d%d.zip" % _sys.version_info[:2]))
    for name in _STDLIB_CRITICAL:
        mod = _sys.modules.get(name)
        if mod is None:
            mod = _il.import_module(name)
        f = getattr(mod, "__file__", None)
        if f is None:   # builtin/frozen: interpreter-owned by construction
            continue
        rp = _os.path.realpath(f)
        # site-segment belt: site-packages/dist-packages are NEVER a stdlib
        # home even when nested UNDER a stdlib root (Debian ships
        # /usr/lib/python3.12/dist-packages inside the stdlib dir) or under
        # a venv lib dir — a critical stdlib module resolving from any such
        # dir is a plant, whoever owns the dir.
        _parts = rp.split(_os.sep)
        if any(_seg in _parts for _seg in ("site-packages", "dist-packages")):   # GUARD SITE-SEG (battery mutation anchor)
            raise StdlibShadowError(
                f"trust REFUSES to import: stdlib module {name!r} resolved "
                f"to {rp} — a site-packages/dist-packages dir is NEVER a "
                f"stdlib home (the site-.pth reorder class: site.py "
                f"executes .pth 'import' lines at startup, before any "
                f"user code). Audit that dir for planted "
                f"modules and .pth files.")
        if not any(rp == r or rp.startswith(r + _os.sep) for r in roots):
            raise StdlibShadowError(
                f"trust REFUSES to import: stdlib module {name!r} resolved "
                f"to {rp} — OUTSIDE every STDLIB root (sysconfig stdlib/"
                f"platstdlib under the BASE installation; site-packages is "
                f"NOT a stdlib home). The sys.path[0]/cwd shadow class "
                f"and the site-.pth reorder class both land here; no env "
                f"var needed. Rerun "
                f"from a clean directory, or launch the interpreter -P/-I, "
                f"and audit site-packages for planted modules/.pth files.")


_stdlib_identity_gate()   # GATE-CALL (battery mutation anchor)

from ._core import (VendorTamperError, OutputRootError, verify,  # noqa: E402
                    load_vendored, run_vendored, vendored_path, out_root,
                    alias_identity)
from ._pins import PINS, LINKED, LINKED_SHAS  # noqa: E402
from . import strata, receipt, fraction_oracle  # noqa: E402

__all__ = ["verify", "load_vendored", "run_vendored", "vendored_path",
           "out_root", "alias_identity", "PINS", "LINKED", "LINKED_SHAS",
           "strata", "receipt", "fraction_oracle",
           "VendorTamperError", "OutputRootError", "EnvPoisonError",
           "StdlibShadowError"]
