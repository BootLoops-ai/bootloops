#!/usr/bin/env python3
"""
cards.py (periods wing) — generalized model-card contract.

A card is one JSON object describing an h_eff-parameter effective CY family +
flux vacuum, with EXACT rational data only (strings parsed as Fraction; no
floats).  Field names are family-generic (nothing KKLT-specific).

Loading/parsing DELEGATES to pipeline/family.py; this facade adds the schema
contract + a validator usable on any directory of cards.
Fail-closed semantics: non-simplicial cards without provenance-gated
`eff_rays`+`gv_table` stay PENDING-CARD (geff_series.gv_extract fail-closes).
Battery: selftest_cards.py (schema over every card in pipeline/cards/ + G2
exact PFFV values on the DKMM card).
"""
import os, sys, json
_PIPE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pipeline")
sys.path.insert(0, _PIPE)
from family import Card, load_card                       # noqa: F401

REQUIRED = {
    "name": "card id string",
    "h": "h_eff (# Mori charge vectors)",
    "charge_vectors": "h Mori vectors, entries int, each sums to 0 (CY)",
    "kappa": "triple intersections {'i,j,k': rational-string}",
    "a_mat": "prepotential quadratic-term matrix (rationals)",
    "c2D": "second-Chern pairings (ints; b = c2D/24)",
    "chi_A": "Euler characteristic (xi = chi_A/2 * zeta3 / v^3)",
    "gv_pinned": "pinned GV leaders {'d1,..,dh': int}",
    "gv_window": "GV extraction window (ints)",
    "M_flux": "PFFV M vector (rational-strings)",
    "K_flux": "PFFV K vector (rational-strings)",
    "Q_D3": "tadpole bound (int)",
}
# required additionally once a card is TRANSPORT-READY (has routes):
TRANSPORT_REQUIRED = {
    "series_depth": "depth params {N_EXT, M_FIT_LO, M_FIT_HI, ...}",
    "s_star": "target point", "tau_pin": "tau pin",
}
OPTIONAL = {
    "n_num": "# numerator Gamma columns (default 1)",
    "sigma": "sign-frame pin", "eps": "sign-frame parity pin",
    "s_star": "target point (rational-string)", "tau_pin": "tau pin",
    "racetrack_Nm": "racetrack integers {m: N_m}",
    "routes": "transport routes {R|C: {s0, rho0, hints}}",
    "ideal": "effective PF ideal (sympy strings; vacuum layer)",
    "module_basis": "D-module basis exponent tuples",
    "op_search": "annihilator search params (rmax/smax/primes/fit_window)",
    "eff_rays": "extreme rays of CY-effective cone (non-simplicial route)",
    "gv_table": "provenance-gated GV classes (required iff non-simplicial)",
}


# Slice-operator cards: a second card CLASS in pipeline/cards/ (schema tag
# "slice-operator-v1"; e.g. hv4-diag-L5.json, the HV4 diagonal-slice order-5
# operator). These are
# derived-operator cards (exact ODE data + receipts), NOT flux-family cards:
# no Card() parse, no transport routes; validated against their own key set
# and returned as a plain dict. Battery: hv4-diag-L5_bank/battery_l5.py.
SLICE_OPERATOR_REQUIRED = {
    "name": "card id string", "schema": "slice-operator-v1",
    "scope": "slice locus statement (what the operator does NOT cover)",
    "n_params": "slice dimension", "order": "ODE order",
    "operator_normal_form": "exact operator data (rational strings)",
    "symbol": "leading symbol", "singular_points": "singular locus",
    "local_exponents": "exponents at each singular point",
    "mum_series_head_c0_c40": "series head (exact ints)",
    "provenance": "derivation receipts", "battery": "battery receipt pointer",
    "status": "honesty line (slice-only scope, transport state)",
}


def validate_card(path):
    """Schema + exactness gate. Flux-family cards return the parsed Card
    (family.py asserts CY condition + gamma-cancellation); slice-operator
    cards validate their own key set and return the raw dict (fail-closed:
    they never enter flux/transport code paths). Raises on any violation."""
    with open(path) as f:
        d = json.load(f)
    if d.get("schema") == "slice-operator-v1":
        miss = [k for k in SLICE_OPERATOR_REQUIRED if k not in d]
        assert not miss, f"{os.path.basename(path)}: missing required {miss}"
        return d
    miss = [k for k in REQUIRED if k not in d]
    if d.get("routes"):
        miss += [k for k in TRANSPORT_REQUIRED if k not in d]
    assert not miss, f"{os.path.basename(path)}: missing required {miss}"
    return Card(d)
