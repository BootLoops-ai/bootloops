#!/usr/bin/env python3
"""
geffseries.py (periods wing) — effective Gamma-series facade.

Delegates to pipeline/geff_series.py:
h-parameter Gamma-series as exact rho-jets over Q[gamma, zeta2, zeta3]
(|alpha| <= 3: CY3 towers stop at log^3 for any h_eff), curve restriction,
h-dim HKTY GV extraction with sign-frame pin, restricted mirror map + N_m
peel.  Includes the NON-SIMPLICIAL handling (g1 route):
  * eff_route(card): SIMPLICIALITY CHECK FIRST — exact extreme-ray reduction
    + facet H-rep (no floats) from optional card field `eff_rays`; cards
    without it keep the legacy simplicial contract byte-identically.
  * gv_extract on non-simplicial cards FAIL-CLOSES unless the provenance-
    gated `gv_table` card field is present (cone-containment gated);
    otherwise the card stays PENDING-CARD (fail-closed semantics).
Mod-p holonomic-rank probes DELEGATE to tools/pf_rank.jl; operator
derivation DELEGATES to tools/annihilator/ via pipeline/restrict_op.py.
Battery: selftest_geffseries.py = run_pipe g1 (DKMM Gamma-series/GV/w0/tower/
racetrack/mirror byte-exact vs the stored reference series).
"""
import os, sys
_PIPE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pipeline")
sys.path.insert(0, _PIPE)
import geff_series as _gs
from geff_series import *          # noqa: F401,F403 — full public surface

# stable aliases for the suite-facing names
pin_sign_frame = _gs.pin_sign_frame
mod = _gs                          # escape hatch: geffseries.mod.<anything>
