#!/usr/bin/env python3
"""Selftest for tools/posq/closure_interface.py — consumer-runnable, no
phylo context needed (no data, no trees, no network; numpy+scipy only).

Gates
  G1 PIN       upstream sha pins verify (sanctioned build present, no drift)
  G2 SPEC      check_closure passes all fixtures; plants P1 (shape) and
               P2 (NaN) fire
  G3 EVIDENCE  |logZ - exact truth| <= recorded tol, per fixture x pinned seeds
  G4 REGISTER  every result carries register='R3-certified-calibrated',
               a calibration_version, and the calibration-transfer caveat
  G5 MISSING   plant P3: infer_hazard_field / posterior_draws raise
               ClosureObjectMissing with non-empty reason + unlock

Exit 0 iff every gate passes. Run:
  python3 <your-checkout>/tools/posq/selftest_closure_interface.py
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import numpy as np  # noqa: E402

import closure_interface as ci  # noqa: E402

FIXTURES = os.path.join(HERE, "closure_fixtures", "fixtures.json")


def build_closure(kind, params):
    """Rebuild a fixture closure from its JSON spec (no pickles)."""
    if kind == "separable_power":
        a = np.asarray(params["a"], dtype=float)
        return lambda X: np.sum(a * np.log(np.clip(X, 1e-300, None)), axis=1)
    if kind == "constant":
        c = float(params["c"])
        return lambda X: np.full(X.shape[0], c)
    if kind == "hazard7_toy":
        d = np.asarray(params["events"], dtype=float)
        E = np.asarray(params["exposure"], dtype=float)
        return lambda X: np.sum(
            d * np.log(np.clip(X, 1e-300, None)) - np.clip(X, 0.0, None) * E,
            axis=1)
    raise ValueError(f"unknown fixture kind {kind!r}")


def main():
    with open(FIXTURES) as fh:
        spec = json.load(fh)
    fails = []
    receipt = []

    def gate(name, ok, detail):
        receipt.append(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}")
        if not ok:
            fails.append(name)

    # ---- G1 PIN ----
    try:
        pins = ci.verify_pins()
        gate("G1-PIN", True, f"{len(pins)} upstream files sha-verified "
                             f"against the pinned R3 build")
    except Exception as e:  # noqa: BLE001 — report, then hard stop
        gate("G1-PIN", False, repr(e))
        print("SELFTEST closure_interface")
        print("\n".join(receipt))
        print("RESULT: FAIL (pins) — no further gates run")
        return 1

    # ---- G2 SPEC (fixtures pass, plants fire) ----
    for fx in spec["fixtures"]:
        try:
            r = ci.check_closure(build_closure(fx["kind"], fx["params"]),
                                 fx["dim"])
            gate(f"G2-SPEC-{fx['id']}", bool(r["ok"]),
                 f"spec-compliant, checks={len(r['checks'])}")
        except Exception as e:  # noqa: BLE001
            gate(f"G2-SPEC-{fx['id']}", False, repr(e))
    # P1 wrong shape
    try:
        ci.check_closure(lambda X: np.zeros((X.shape[0], 2)), 3)
        gate("G2-PLANT-P1-shape", False, "planted shape violation NOT caught")
    except ci.ClosureInterfaceError:
        gate("G2-PLANT-P1-shape", True, "planted (n,2) output fired")
    # P2 NaN
    try:
        ci.check_closure(lambda X: np.full(X.shape[0], np.nan), 3)
        gate("G2-PLANT-P2-nan", False, "planted NaN NOT caught")
    except ci.ClosureInterfaceError:
        gate("G2-PLANT-P2-nan", True, "planted NaN output fired")

    # ---- G3 EVIDENCE + G4 REGISTER ----
    for fx in spec["fixtures"]:
        f = build_closure(fx["kind"], fx["params"])
        truth = float(fx["truth_lnZ"])
        g = fx["gate"]
        worst = 0.0
        reg_ok = True
        for seed in g["seeds"]:
            res = ci.evaluate_evidence(f, fx["dim"], budget=g["budget"],
                                       seed=seed, name=fx["id"])
            worst = max(worst, abs(res["logZ"] - truth))
            reg_ok &= (res["register"] == ci.MC_REGISTER
                       and bool(res["calibration_version"])
                       and bool(res["calibration_transfer"])
                       and res["certificate"]["register"] == ci.MC_REGISTER)
        gate(f"G3-EVIDENCE-{fx['id']}", worst <= g["tol"],
             f"max|logZ-truth| {worst:.3e} vs tol {g['tol']:.0e} "
             f"({len(g['seeds'])} seeds @ budget {g['budget']})")
        gate(f"G4-REGISTER-{fx['id']}", reg_ok,
             "register label + calibration_version + transfer caveat on "
             "every result" if reg_ok else "register labeling incomplete")

    # ---- G5 MISSING (plant P3) ----
    for fn_name in ("infer_hazard_field", "posterior_draws"):
        try:
            getattr(ci, fn_name)()
            gate(f"G5-MISSING-{fn_name}", False, "did NOT raise")
        except ci.ClosureObjectMissing as e:
            ok = bool(e.reason) and bool(e.unlock)
            gate(f"G5-MISSING-{fn_name}", ok,
                 "ClosureObjectMissing with reason+unlock" if ok
                 else "raised but reason/unlock empty")
        except Exception as e:  # noqa: BLE001
            gate(f"G5-MISSING-{fn_name}", False, f"wrong exception {e!r}")

    print("SELFTEST closure_interface — read-only "
          "evaluation tier")
    print("\n".join(receipt))
    n_gates = len(receipt)
    if fails:
        print(f"RESULT: FAIL ({len(fails)}/{n_gates} gates): "
              + ", ".join(fails))
        return 1
    print(f"RESULT: PASS ({n_gates}/{n_gates} gates)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
