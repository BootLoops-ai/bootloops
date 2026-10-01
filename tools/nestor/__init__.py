#!/usr/bin/env python3
r"""NESTOR -- NESted Tanh-sinh ORacle.

The pipeline's high-precision integration oracle, consolidated under one
name: one entry point over the node-farmed
tanh-sinh engine (oracle_par.py), the per-chain convergence prober
(quad_probe), and the wp-scaled precision ladder of the originating
certification evaluators.  It reports only the digits successive refinements agree on,
and refuses -- named, fail-closed, receipted -- when they do not agree.

FRONT DOOR
==========
    import nestor
    v   = nestor.integrate(f, 0, 1, dps=30)          # value or refusal
    rec = nestor.oracle(f, 0, 1, dps=30, spec=...)   # full JSON receipt

    nestor.probe(...)   -> per-chain convergence prober (quad_probe),
                           delegated to its canonical home
                           tools/gatekeeper/probes/quad_probe.py --
                           ONE source of truth, no duplicated math.

MODULES
=======
    ladder.py    per-level precision ladder, refine-until-agreement,
                 fail-closed (src: the originating evaluators, shas in header)
    farm.py      node farm: explicit tanh-sinh grid over a multiprocessing
                 Pool, grid-nesting self-consistency, orphan-proof pool
                 discipline (src: the originating oracle_par.py, sha in header)
    singcheck.py load-time singularity verification of the integrand spec
    receipts.py  JSON receipts: per-level walls, agreement digits, refusals
    dispersion/  singularity-subtracted dispersion quadrature: (1/pi) int
                 rho K dw' with known threshold turn-ons in rho subtracted on
                 a sub-panel and added back in closed form (exact power-log
                 moments x kernel Taylor coefficients); mpmath only
                 (from nestor.dispersion import disp_subtracted, ...)
    examples/    box1_dilog.py -- a worked analytic dispersion kernel (the
                 one-loop light-by-light box in 1D dilogarithmic form)
    selftest.py  positive + mutation controls (python3 -m nestor.selftest),
                 including the dispersion legs in tests/

PROVENANCE
==========
    nestor.farm ports the node-farm engine of an originating parallel
    oracle (not part of this repo) verbatim, generalizing only the plumbing;
    nestor.ladder extracts the precision-ladder pattern of the originating
    certified-quadrature evaluators, with the source sha recorded in its
    header; nestor.probe delegates to the canonical
    tools/gatekeeper/probes/quad_probe.py at import time (computed sibling
    path).

No absolute paths anywhere in this package: the only path
arithmetic is computed relative to __file__ at runtime.
"""
from __future__ import annotations

import os as _os
import sys as _sys

from .ladder import (AbortRefusal, CertFail, EXIT_GATE, EXIT_MISSING,
                     EXIT_OK, EXIT_PIN, EvalRefusal, GateRefusal,
                     MissingRefusal, NestorRefusal, PinRefusal, SpecRefusal,
                     cert_quad, crank_check, digits, escalation_schedule,
                     inner_tol_exp, ladder, prec_for, wp_for)
from .farm import farm_campaign, farm_integrate, ts_nodes
from .oracle import integrate, oracle
from .receipts import NESTOR_VERSION, SRC_SHAS
from .singcheck import verify as singcheck
from . import dispersion
from .dispersion import disp_subtracted

__version__ = NESTOR_VERSION

# --- nestor.probe: delegate to the canonical quad_probe (one source of
# truth).  tools/nestor/ and tools/gatekeeper/ are siblings; the entry
# is COMPUTED from __file__, never hard-coded.
_TOOLS_DIR = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))


def _load_probe():
    if _TOOLS_DIR not in _sys.path:
        _sys.path.insert(0, _TOOLS_DIR)
    try:
        from gatekeeper.probes import quad_probe as _qp
        return _qp
    except Exception as e:                     # pragma: no cover
        raise MissingRefusal(
            "nestor.probe: canonical quad_probe unavailable "
            "(gatekeeper.probes.quad_probe next to the nestor package): "
            "%s" % e)


try:
    quad_probe = _load_probe()
    probe = quad_probe.probe
    EnvEvaluator = quad_probe.EnvEvaluator
except MissingRefusal:                          # deferred: refuse at call
    quad_probe = None

    def probe(*a, **k):
        _load_probe()                          # raises the named refusal

    EnvEvaluator = None

__all__ = [
    "integrate", "oracle", "probe", "quad_probe", "EnvEvaluator",
    "cert_quad", "ladder", "crank_check", "digits", "escalation_schedule",
    "wp_for", "prec_for", "inner_tol_exp",
    "farm_integrate", "farm_campaign", "ts_nodes", "singcheck",
    "dispersion", "disp_subtracted",
    "NestorRefusal", "PinRefusal", "GateRefusal", "MissingRefusal",
    "SpecRefusal", "CertFail", "AbortRefusal", "EvalRefusal",
    "EXIT_OK", "EXIT_PIN", "EXIT_GATE", "EXIT_MISSING",
    "NESTOR_VERSION", "SRC_SHAS",
]
