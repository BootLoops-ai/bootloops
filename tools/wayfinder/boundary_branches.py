#!/usr/bin/env python3
r"""
boundary_branches.py — thin, import-guarded delegator to
tools/frobenius_boundary (var = INFINITY Frobenius-branch machinery).

This module vendors NOTHING: every call goes straight to
frobenius_boundary.core (which has its own selftest —
`cd tools && python3 -m frobenius_boundary.selftest`, with pinned
acceptance numbers). It exists so wayfinder users find
the infinity-branch route next to the finite-point routes, with the routing
and the standing caveat in ONE place.

WHICH FROBENIUS TOOL WHEN
=========================
  * wayfinder.frobenius   — MATRIX first-order systems LANDING at a FINITE
    regular-singular point: fundamental basis Y(u) = P(u) u^M (+ shearing
    for nonzero-integer indicial gaps), least-squares landing of a
    transported vector. Real-eps0 matching policy.
  * wayfinder.frob_scalar — SCALAR operators / P-recurrences with full
    resonant integer towers: transport of a series value TO the singular
    point (per-root eps-jet basis; validated at 260d two-precision).
  * frobenius_boundary (THIS wrapper) — branch CLASSIFICATION at
    var = INFINITY for cut-block DEs: exact poly-DE
    from a Kira table, exponent spectrum via D1 = R_Delta - diag(L*d/2+ai),
    multi-eps rational-linear classification, strata exclusion, branch
    series + Phi matrix for the anchor solve. Use it when the boundary
    anchor sits at infinity of a mass-like variable and the physical
    solution is NOT analytic there.

THE RESOLVED RULE: "exponents = eig(D1)"
is exact only when the N-frame is Fuchsian. With a GENUINE nilpotent D_0
block the TRUE exponent multiset = eig of the shear-reduced
residue, and the killed complex/irrational eig(D1) values are exponents of
NO solution (verified independently on a genuine-D_0 case:
counting is 22 survivors + 5 flagged + 3 killed = 27 allowed, NOT a
naive 12). The cli's classify/exclusion stages now run on the reduced
residue automatically when a genuine D_0 is detected (see
spec_json['classification_frame']); the CORE API delegated here still
classifies whatever eigenvalue lists YOU pass it — for genuine-D_0
systems pass the reduced-residue spectrum, not eig(D1). Also do NOT set
kill_integer_eps_indep=True unless the vacuum theorem holds for YOUR
family (verify with a kira_vacuum IBP run).

MINOR UPSTREAM WART (flagged, not fixed here — this wrapper never edits the
delegate): classify()'s GROUP-level 'verify_maxres_alleps' maps an EXACT-ZERO
family residual to 1 (`x['verify_maxres_alleps'] or 1`, classify_mod.py:96 —
0.0 is falsy). Read the per-FAMILY residuals for the honest all-eps
verification number; the group field is only trustworthy when nonzero.

IMPORT GUARD / dps HYGIENE
==========================
frobenius_boundary mutates ambient state at import (raises global mp.dps
to >= 50, sets OPENBLAS_NUM_THREADS via setdefault) and pulls sympy. This
wrapper therefore imports it LAZILY inside mp.workdps(>=50): the package's
`if mp.mp.dps < 50` guard sees a compliant ambient value and global mp.dps
is NEVER mutated (wayfinder design rule 4). Callers should still pass
explicit prec=/dps arguments or wrap calls in mp.workdps — several
frobenius_boundary entry points default their working precision to the
AMBIENT mp.mp.dps at call time.

Delegated surface (signatures documented in frobenius_boundary/core.py):
    build_poly_DE, spectrum, classify, exclude_strata, branch_series,
    phi_matrix, seed_vectors, eps_to_d
plus available() / get_module() from this wrapper.
"""

import os
import sys

import mpmath as mp

__all__ = ["available", "get_module", "build_poly_DE", "spectrum",
           "classify", "exclude_strata", "branch_series", "phi_matrix",
           "seed_vectors", "eps_to_d"]

_MOD = None
_IMPORT_ERR = None


def _fb():
    """Lazy, guarded import of frobenius_boundary.core (see module doc)."""
    global _MOD, _IMPORT_ERR
    if _MOD is not None:
        return _MOD
    if _IMPORT_ERR is not None:
        raise _IMPORT_ERR
    tools_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if tools_dir not in sys.path:
        sys.path.insert(0, tools_dir)
    try:
        # inside workdps(>=50) the package's `if mp.mp.dps < 50` import
        # guard is a no-op => ambient global mp.dps stays untouched
        with mp.workdps(max(mp.mp.dps, 50)):
            from frobenius_boundary import core as _core  # noqa: PLC0415
        _MOD = _core
        return _MOD
    except ImportError as exc:
        _IMPORT_ERR = ImportError(
            f"boundary_branches: tools/frobenius_boundary is not importable "
            f"({exc}). It lives at {tools_dir}/frobenius_boundary (own "
            f"selftest: cd {tools_dir} && python3 -m "
            f"frobenius_boundary.selftest) and needs sympy. This wrapper "
            f"delegates only — nothing is vendored here.")
        raise _IMPORT_ERR from exc


def available():
    """True iff tools/frobenius_boundary imports cleanly (lazy, cached)."""
    try:
        _fb()
        return True
    except ImportError:
        return False


def get_module():
    """The frobenius_boundary.core module itself (for anything not wrapped)."""
    return _fb()


# ---- thin delegates (docs: frobenius_boundary/core.py + its README) ----

def build_poly_DE(*args, **kwargs):
    """Exact polynomial DE D(var)*dM/dvar = G(var)*M from a per-point Kira
    table. Delegates to frobenius_boundary.core.build_poly_DE; extras of the
    most recent build live in get_module().LAST_BUILD."""
    return _fb().build_poly_DE(*args, **kwargs)


def spectrum(*args, **kwargs):
    """Exponent spectrum at var=inf via the D1 eigenproblem (SEE the eig(D1)
    FLAG in the module docstring). Delegates to core.spectrum."""
    return _fb().spectrum(*args, **kwargs)


def classify(*args, **kwargs):
    """Multi-eps rational-linear classification lambda = a + b*eps, verified
    at ALL provided eps (4-eps pattern). kill_integer_eps_indep=True ONLY
    when the vacuum theorem holds for your family. Delegates to
    core.classify."""
    return _fb().classify(*args, **kwargs)


def exclude_strata(*args, **kwargs):
    """Survivor strata + exclusion hardening rows. Delegates to
    core.exclude_strata (re-export of classify_mod.exclude_strata)."""
    return _fb().exclude_strata(*args, **kwargs)


def branch_series(*args, **kwargs):
    """Branch solution series at var=inf (resonance-tolerant D_p-block
    recursion; genuine nilpotent D_0 handled by ITS OWN exact shear — see
    the frobenius_boundary README). Delegates to core.branch_series."""
    return _fb().branch_series(*args, **kwargs)


def phi_matrix(*args, **kwargs):
    """N x J branch-solution matrix Phi for the anchor solve
    R_rows*Phi*c = O_rows. Delegates to branches_mod.phi_matrix."""
    return _fb().phi_matrix(*args, **kwargs)


def seed_vectors(*args, **kwargs):
    """Nullspace of (D1 - lam I): Frobenius seeds of branch lam. Delegates
    to core.seed_vectors."""
    return _fb().seed_vectors(*args, **kwargs)


def eps_to_d(*args, **kwargs):
    """'1/101' | Fraction | float -> (dd_rat sympy, dd mp, eps mp).
    Delegates to core.eps_to_d."""
    return _fb().eps_to_d(*args, **kwargs)
