#!/usr/bin/env python3
"""
conifold.py (coni_pfv, periods wing) — coni-PFV chart: curve/frame
gates CP0-CP10, conifold Frobenius block, fail-closed transport router,
resonant log-tower machinery.  Delegates to the modules under
pipeline/conipfv/:

  * coni_pfv_curve(card)  [conipfv/coni_curve.py]: CP0 routing (PFFV cards
    with K.N^-1.K = 0 return NOT-CONI and stay with fluxcurves/g2), CP1-CP6
    exact lemma gates (unique conifold class, Diophantine K.p = 0, flux
    integrality, tadpole, C1 flux vectors, crosswalk gate).  EXACT.
  * g3_series(card, cur)  [conipfv/coni_frobenius.py]: CP7-CP10 — exact
    racetrack ladder, Frobenius block (c_log = -M_cf n_cf, c_tau = lam
    exact identities), joint-F-term vev + z_cf vev + W0 racetrack.  The
    numeric legs are LABELED ESTIMATES (pinned GV window); certified value
    = transport legs, never this.
  * transport_ready / landing_ball  [conipfv/coni_transport.py]: FAIL-CLOSED
    router (leg 1 refuses without coni_op/coni_towers/coni_routes/s_star)
    + leg-3 certified landing ball (acb end-to-end; never reduced through float64).
  * frobenius_tower / verify_tower / theta_from_json / routes_and_star /
    c2_bound  [conipfv/coniop/coni_pack.py]: exact eps-jet resonant
    Frobenius towers from a theta-form operator, with the L_s-annihilation
    verifier as the fail-closed gate (log-branch mutations MUST be caught —
    see battery).  DKMM-regression-proven.
  * run_card  [conipfv/run_conipfv.py]: per-card orchestrator + verdicts
    (NOT-CONI / CONI-PENDING / CONI-READY / CONI-BLOCKED / CONI-FAIL).

HONESTY LINES (verdict scope):
  - float pin scoping: vacuum values (tau, z) taken from a float64 source
    are POINT PINS; any W0 verdict is CONDITIONAL on them, floor ~1e-13
    rel — never quote tighter.
  - m2 amendment: on coni charts the 56-tower rank-saturation fit is
    STRUCTURALLY UNAVAILABLE (level-infinite jet towers); uniqueness is by
    zero-kernel EXACT identities (on-ray tau-sector cancellation, c_tau =
    lam, crosswalk-unique scan) — genuinely weaker than a rank proof; say so.
  - s_star / coni_C2-tail / route clearances are labeled estimates.
Battery: selftest_conifold.py (DKMM NOT-CONI control, DKMM tower +
log-branch-flip mutation, fail-closed transport).
"""
import os, sys

_PIPE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pipeline")
_CONI = os.path.join(_PIPE, "conipfv")
for _p in (_PIPE, _CONI, os.path.join(_CONI, "coniop")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from coni_curve import coni_pfv_curve                       # noqa: F401
from coni_frobenius import (g3_series, ladder,              # noqa: F401
                            frobenius_block)
from coni_transport import (transport_ready, run_leg1,      # noqa: F401
                            landing_ball, REQUIRED)
from coni_pack import (theta_from_json, frobenius_tower,    # noqa: F401
                       verify_tower, routes_and_star, c2_bound, Jet)
from run_conipfv import run_card                            # noqa: F401
