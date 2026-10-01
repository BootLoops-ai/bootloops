#!/usr/bin/env python3
"""Battery (cards.py): schema over every card in pipeline/cards/ + G2 exact."""
import glob, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cards
from fractions import Fraction as Fr

n_ok, n_pend, n_slice = 0, 0, 0
for fn in sorted(glob.glob(os.path.join(HERE, "pipeline", "cards", "*.json"))):
    if fn.endswith(".pre-nonsimp-route"):
        continue
    c = cards.validate_card(fn)
    n_ok += 1
    if isinstance(c, dict):  # slice-operator class (returns raw dict)
        n_slice += 1
        continue
    if c.eff_rays is not None and c.gv_table is None:
        n_pend += 1          # non-simplicial without provenance table:
        # PENDING-CARD (gv_extract fail-closes) — schema-valid, not runnable
print(f"[cards schema] PASS  {n_ok} cards validate "
      f"({n_pend} fail-closed PENDING-CARD, {n_slice} slice-operator)")
# shipped fixture set: dkmm, ads-5-81-3213, hv4-diag-L5
assert n_ok >= 3, "shipped fixture card set shrank below 3"
assert n_slice >= 1, "slice-operator card (hv4-diag-L5) missing"

# G2 exact PFFV values on the DKMM card
sys.path.insert(0, os.path.join(HERE, "pipeline"))
from curve_from_flux import pffv_curve
dk = cards.load_card(os.path.join(HERE, "pipeline", "cards", "dkmm.json"))
pc = pffv_curve(dk)
assert pc["p"] == [Fr(2, 5), Fr(3, 10)] and pc["nu"] == (4, 3)
assert pc["F"] == [7, 3, -24, 0, -16, 50] and pc["H"] == [0, 3, -4, 0, 0, 0]
assert pc["tadpole"] == 124 <= 138
print("[cards G2] PASS  DKMM PFFV curve/flux/tadpole exact")
print("SELFTEST cards.py: ALL PASS")
