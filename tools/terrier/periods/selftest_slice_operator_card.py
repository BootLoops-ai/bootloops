#!/usr/bin/env python3
"""Battery PER-SC (slice-operator card): replay the hv4-diag-L5 card battery
(battery_l5.py, in the card's data directory) — fit == printed JKK/GvdH L5,
minimality, symbol, exponents, MUM regeneration, jet gates, K3/L3 subsector,
pfaffian payload — and require the regenerated BATTERY_RECEIPT.json to equal
the reference receipt field-for-field.
Data: periods/pipeline/cards/hv4-diag-L5_bank/ (sha-pinned)."""
import json, os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
BANK = os.path.join(HERE, "pipeline", "cards", "hv4-diag-L5_bank")
r = subprocess.run([sys.executable, "battery_l5.py"], cwd=BANK,
                   capture_output=True, text=True, timeout=300)
out = r.stdout + r.stderr
print(out[-1500:])
assert r.returncode == 0, "battery_l5 rc != 0"
assert "BATTERY PASS" in out, "BATTERY PASS line missing"
got = json.load(open(os.path.join(BANK, "BATTERY_RECEIPT.json")))
want = json.load(open(os.path.join(BANK, "BATTERY_RECEIPT_reference.json")))
VOLATILE = ("elapsed_s",)                    # tolerated if present; the receipt no longer carries it
for k in VOLATILE:
    got.pop(k, None); want.pop(k, None)
assert set(got) == set(want), f"receipt key sets differ: {set(got) ^ set(want)}"
bad = [k for k in want if got[k] != want[k]]
assert not bad, f"receipt fields != reference: {bad}"
print("SELFTEST slice_operator_card: ALL PASS "
      f"(battery PASS + {len(want)} receipt fields == reference)")
