#!/usr/bin/env python3
"""
transport.py (periods wing) — certified-transport chassis facade.

Consumes (operator, curve, flux, frame) GENERICALLY via the card + stored
JSON layer; delegates to the pipeline modules:
  * pipeline/pipe_transport.py — stepped majorant legs: MUM envelopes from
    the FACTORED indicial (lead prod(1-e_i/m)^mult), D-form local legs with
    PROVED descriptor tails (abs-majorant, FC1 semantics), ball-trim,
    schedule(s0, hints, s_star), run_route, exact symplectic contract.
    Ball arithmetic (python-flint arb/acb) END-TO-END — never reduce a ball
    through float64 (arb-to-float64 reduction is a known pitfall).
  * pipeline/pipe_vac.py — jet layer (generic normal-form closure: singles
    cascade + coupled Gaussian over Q(z), leftover-relations gate; r3
    D-module-route == annihilator-route identity = ring-identity frame
    fixing), connection matrices, curve M(s), certified V0 + off-curve leg,
    ball-AD, KRAWCZYK existence+uniqueness in 2h+2 real coords, FD-vs-AD
    Jacobian gate.
Operator derivation DELEGATES to tools/annihilator/ through
pipeline/restrict_op.py (wrapper only; 2-prime x 2-truncation pin + CRT/Wang
+ exact annihilation proof).  Julia counterpart:
upgrades/Eichler.jl/src/cy_transport.jl, whose last-window x4 tail is
empirical; the proven tail envelope is envelope_certified.py
(PROVEN_MAJORANT.md).
Battery: selftest_transport.py — reproduce the reference DKMM certified
balls EXACTLY (midpoint string-identical, radius ratio 1.00) through this
chassis: W0 on-curve + certified vacuum.
"""
import os, sys
_PIPE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pipeline")
sys.path.insert(0, _PIPE)
import pipe_transport as PT                              # noqa: F401
import pipe_lib as PL                                    # noqa: F401
from pipe_transport import (schedule, run_route, contract,   # noqa: F401
                            mum_eval, local_leg, sympl, fr2arb, fr2acb)


def configure(card, op, tower, gv, curve):
    PL.configure(card, op, tower, gv, curve)
    PT.configure()
    return PL, PT


def vac_layer():
    """Krawczyk vacuum-certification layer — returns pipe_vac (lazy import:
    it pulls sympy + flint).  Certified entry points are orchestrated per
    card as in pipeline/run_pipe.py jets()/vac()/fdgate()."""
    import pipe_vac as PV
    return PV
