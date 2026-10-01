#!/usr/bin/env python3
"""
verdicts.py (periods wing) — wing-facing verdict layer.

DELEGATES every receipt/gate primitive to common/verdict.py: parse_ball,
matched_digits, ball_gate (like-for-like; radius bar), route_gate
(cross-route; NO radius bar — independent routes carry independent radii,
so only midpoint digits are compared), two_dps_gate, code_sha, receipt
(atomic, sha-stamped, schema terrier-receipt-v1).

Adds the periods-wing MUTATION-CONTROL HOOKS (a certification run must fail
under each mutation to count):
  * xi-flip (f3): fluxcurves.xi_flipped — solve must fail or demand zeta3
    counterterms;
  * abs-majorant control (FC1): naive tail summation must FAIL where the
    proved abs-majorant certifies (DKMM example: s_vac ratio 0.976, radius
    0.0138);
  * two-dps: lo-dps vs hi-dps midpoints beyond bar;
  * dual-route: independent R/C transport routes agree (route_gate).
Battery: selftest_verdicts.py = primitive unit tests + the full 13/13 DKMM
scoreboard (run_pipe gates) regenerated in place.
"""
import os, sys
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "common"))
sys.path.insert(0, os.path.join(_HERE, "pipeline"))
from verdict import (verdict, parse_ball, matched_digits,     # noqa: F401
                     ball_gate, route_gate, two_dps_gate, code_sha, receipt)

MUTATION_CONTROLS = {
    "xi_flip_f3": "fluxcurves.xi_flipped context; clean flipped fit = FAIL",
    "abs_majorant_FC1": "naive tail sum must fail inside proved envelope",
    "two_dps": "two_dps_gate(lo, hi, bar)",
    "dual_route": "route_gate(R, C, bar) — independent transport routes",
}
