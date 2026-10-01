#!/usr/bin/env python3
"""selftest_pfaffian.py — battery for periods/pfaffian.py (see PFAFFIAN.md).
P0 route-choice front door (scalar-elimination refusal law)
P2 fresh-prime gate replay: ads-5-81 (512 nodes), 0 diffs
P3 MUTATION control: one perturbed connection coeff MUST be caught by P2 gate
P6 ads-5-81 certified W0 ball dps 60+90 == stored verdict (mid AND radius)
(The route-B legs P1/P4/P5 — lift replay, annihilation replay and landing
balls of the rank-16 conifold-frame example — are not part of this release:
that example is not included in the package.)
Needs the reference receipts, which are not included in the package: set
TERRIER_PFAFFIAN_BANK to their root (read-only). Resource caps: ulimit -v
32505856, nice >= 5."""
import ast, copy, json, os, sys, time
from fractions import Fraction as Fr
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "pipeline"))
import pfaffian as PF

BANK = os.environ.get("TERRIER_PFAFFIAN_BANK")
if not BANK:
    sys.exit("REFUSE: reference data not included in the package "
             "(Pfaffian transport receipts: per-prime samples, exact "
             "connections, gate records, landing balls) -- set "
             "TERRIER_PFAFFIAN_BANK to run")
ADS = os.path.join(BANK, "pfaffian", "ads581")
FAIL = []


def gate(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name} {detail}", flush=True)
    if not ok:
        FAIL.append(name)


# ------------------------------------------------------------------- P0
def run_P0():
    r = PF.route_choice("point-value")
    gate("P0a point-value -> matrix transport", r["route"] == "matrix-transport")
    try:
        PF.route_choice("scalar-operator", probe=(13, 4052))
        gate("P0b ads-5-81 probe refused", False, "no refusal raised")
    except PF.ScalarEliminationRefusal as e:
        gate("P0b ads-5-81 probe refused", "REFUSED" in e.advice
             and "matrix" in e.advice.lower(), f"slots 13x4053")
    try:
        PF.route_choice("scalar-operator")
        gate("P0c scalar without probe refused", False, "no refusal raised")
    except PF.ScalarEliminationRefusal:
        gate("P0c scalar without probe refused", True)
    r = PF.route_choice("scalar-operator", probe=(4, 60))
    gate("P0d small scalar accepted", r["route"] == "scalar-elimination-ok"
         and r["slots"] == 244)


# ------------------------------------------------------------------- P2/P3
def run_P2_P3():
    AXa = json.load(open(os.path.join(ADS, "A_exact.json")))
    Da = np.load(os.path.join(
        BANK, "cards", "ads-5-81-3213_bank", "dmodule",
        "conn_samples_2147483489.npz"))
    t0 = time.time()
    ra = PF.fresh_prime_gate(AXa, 2147483489, Da["sig"], Da["A"])
    ga = json.load(open(os.path.join(ADS, "gate_PA1.json")))
    gate("P2b ads-5-81 fresh-prime gate 0 diffs",
         ra["pass"] and ra["diffs"] == ga["diffs"] == 0
         and ra["nodes"] == ga["nodes"],
         f"{ra['entry_values_checked']} values ({time.time()-t0:.0f}s)")
    MU = copy.deepcopy(AXa)
    k = sorted(MU["entries"])[0]
    MU["entries"][k][0] = str(Fr(MU["entries"][k][0]) + 1)
    rm = PF.fresh_prime_gate(MU, 2147483489, Da["sig"], Da["A"])
    gate("P3 mutation control: perturbed coeff CAUGHT", rm["diffs"] > 0,
         f"entry {k} +1 -> {rm['diffs']} diffs (must-fail law)")


# ------------------------------------------------------------------- P6
def run_P6():
    import family
    import pipe_lib as PL
    card = family.Card(json.load(open(os.path.join(
        HERE, "pipeline", "cards", "ads-5-81-3213.json"))))
    g2 = json.load(open(os.path.join(
        BANK, "cards", "ads-5-81-3213_bank", "gate_g2.json")))
    MTOW = 120
    tw = json.load(open(os.path.join(ADS, "tower_120.json")))
    C = {tuple(map(int, k.split(","))):
         [{tuple(map(int, kk.split(","))): Fr(v) for kk, v in d.items()}
          for d in lst] for k, lst in tw["C"].items()}
    PL.H = card.h; PL.NU = (9, 10, 10, 10, 10)
    PL.M_TOW = MTOW; PL.C_JET = C
    fr = json.load(open(os.path.join(ADS, "f3_frame.json")))
    coeffs = [{ast.literal_eval(k): Fr(v) for k, v in cd.items()}
              for cd in fr["coeffs"]]
    als = []
    for cd in coeffs:
        for (al, p, q) in cd:
            if al not in als:
                als.append(al)
    phi = {al: PL.build_phi(al) for al in als}
    t0 = time.time()
    towers = PF.frame_towers(coeffs, phi, MTOW)
    print(f"  [P6] towers rebuilt from stored f3_frame "
          f"({time.time()-t0:.0f}s)", flush=True)
    out = {}
    for dps in (60, 90):
        t0 = time.time()
        Pi, v = PF.direct_sum_periods(towers, coeffs, card.s_star, MTOW,
                                      dps, PF.bound_level_ads581)
        w = PF.contract_w0(Pi, v, g2["F"], g2["H"], 12, card.tau_pin)
        out[dps] = w["W0"]
        print(f"  [P6] dps {dps}: W0 = {w['W0']} ({time.time()-t0:.0f}s)",
              flush=True)
    bank = json.load(open(os.path.join(ADS, "verdict_raw.json")))
    gate("P6 ads-5-81 W0 ball dps 60+90 byte-exact vs stored verdict",
         str(out[60]) == bank["W0_60"] and str(out[90]) == bank["W0_90"],
         f"W0_90[:40] {str(out[90])[:40]}")


if __name__ == "__main__":
    T0 = time.time()
    run_P0()
    run_P2_P3()
    run_P6()
    n = 7
    print(f"[selftest_pfaffian] {n - len(FAIL)}/{n} passed "
          f"({time.time()-T0:.0f}s wall)"
          + ("" if not FAIL else f"  FAILED: {FAIL}"))
    sys.exit(1 if FAIL else 0)
