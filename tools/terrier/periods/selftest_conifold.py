#!/usr/bin/env python3
"""selftest_conifold.py — coni_pfv battery (see conifold.py).
C1 DKMM NOT-CONI control routes correctly (== stored gates JSON)
C4 DKMM tower regression: exact resonant Frobenius to s^80, L_s-annihilation
   verified, branch-0 == stored curve series
C5 MUTATION (must-fail): log-branch flip AND branch swap both caught
C6 transport stays FAIL-CLOSED on a card without coni_op (no shipped card
   carries it, by design)
(Legs C2 and C3, the conifold-PFV card replays, are not part of this release:
their example card is not included in the package.)"""
import json, os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import conifold as C

PIPE = os.path.join(HERE, "pipeline")
CONI = os.path.join(PIPE, "conipfv")
sys.path.insert(0, PIPE)
from family import load_card
NP = 0


def ok(name, cond, detail=""):
    global NP
    assert cond, f"{name} FAIL {detail}"
    NP += 1
    print(f"[PASS] {name}" + (f" -- {detail}" if detail else ""))


# C1 — DKMM NOT-CONI control
rec = C.run_card(os.path.join(PIPE, "cards", "dkmm.json"), None)
bank = json.load(open(os.path.join(CONI, "out", "dkmm_P11169.gates.json")))
ok("C1 dkmm NOT-CONI control", rec["verdict"] == "NOT-CONI" and rec == bank,
   "routing + note == stored")

# C4 — DKMM tower regression (coni_pack machinery on the stored operator)
op = json.load(open(os.path.join(PIPE, "dkmm", "operator_LS.json")))
theta = C.theta_from_json(op)
mult0 = op["indicial"]["mult0"]
extra = sum(m for e, m in op["indicial"]["factors"] if e != 0)
T = C.frobenius_tower(theta, mult0, 80, mult0 + 2 * extra, extra)
C.verify_tower(theta, T, 80)
ser = json.load(open(os.path.join(PIPE, "dkmm", "series_curve.json")))
a = [Fr(x) for x in ser["a"]]
n = min(81, len(a))
ok("C4 DKMM tower + L_s-annihilation + branch-0 == curve series",
   all(T[0][m] == a[m] for m in range(n)),
   f"{mult0} log branches to s^80; {n} branch-0 terms exact")

# C5 — mutation controls MUST be caught
for name, mut in (("log-branch flip", lambda t: [t[0]] +
                   [[-x for x in t[1]]] + t[2:]),
                  ("branch swap", lambda t: [t[0], t[2], t[1], t[3]])):
    Tm = mut([row[:] for row in T])
    caught = False
    try:
        C.verify_tower(theta, Tm, 80)
    except AssertionError:
        caught = True
    ok(f"C5 mutation caught: {name}", caught)

# C6 — fail-closed transport (no shipped card carries coni_op)
card = load_card(os.path.join(PIPE, "cards", "dkmm.json"))
ok("C6a transport_ready fail-closed", C.transport_ready(card) is False)
closed = False
try:
    C.run_leg1(card, None, 200)
except AssertionError:
    closed = True
ok("C6b run_leg1 refuses without coni_op", closed)

print(f"selftest_conifold: {NP}/{NP} PASS")
