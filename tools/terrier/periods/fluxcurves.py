#!/usr/bin/env python3
"""
fluxcurves.py (periods wing) — PFFV curve + exact flux contraction
+ frame dictionary, as code.  Delegates to the pipeline modules:

  * pffv_curve(card)  [pipeline/curve_from_flux.py]: N = kappa.M,
    p = N^-1 K, nu = d*p primitive -> monomial curve z~_a = s^{nu_a};
    integer flux vectors F/H (freeze C1); exact tadpole check vs Q_D3.
    EXACT rational linear algebra only (no floats).
  * frame dictionary [pipeline/pipe_lib.py solve_frame after configure]:
    integral-frame dictionary v^3*Pi_i = sum c * v^p z3^q * Phi_alpha —
    the dressed(zeta3/xi counterterms)/undressed(pure Phi) translation.
    XI-TOGGLE as code: xi_flipped(PL, card) context = mutation control f3
    (flipping xi MUST demand zeta3 counterterms or fail the solve; a clean
    flipped fit is a FAIL of the battery: the xi sign is a fixed
    convention, not a fit parameter).
  * exact symplectic contraction [pipeline/pipe_transport.py sympl/contract]
    (2h+2 frame; W0/absW/emK gauge outputs).
Battery: selftest_fluxcurves.py = G2 exact values + stored curve/flux
regression + stage1 G4 frame/f3/towers gates (run_pipe stage1).
coni-PFV cards (K.N^-1.K != 0): OUT OF SCOPE here — they route through
conifold.py / pipeline/conipfv (gates in pipeline/conipfv/VALIDATION.md).
"""
import os, sys
from contextlib import contextmanager
_PIPE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pipeline")
sys.path.insert(0, _PIPE)
from curve_from_flux import pffv_curve, solve_exact      # noqa: F401
import pipe_lib as PL                                    # noqa: F401
import pipe_transport as PT                              # noqa: F401


def configure(card, op, tower, gv, curve):
    """One-stop chassis configure: exact layer + certified layer."""
    PL.configure(card, op, tower, gv, curve)
    PT.configure()
    return PL, PT


def solve_frame(card, **kw):
    """Frame dictionary solve (card must be configured first)."""
    return PL.solve_frame(**kw)


@contextmanager
def xi_flipped(card):
    """f3 mutation control: temporarily flip the xi pin. Callers must verify
    the solve either FAILS or demands zeta3 counterterms (see battery)."""
    PL.XI_COEF = -card.xi_coef
    try:
        yield
    finally:
        PL.XI_COEF = card.xi_coef
