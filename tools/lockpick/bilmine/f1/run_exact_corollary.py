#!/usr/bin/env python3
# lockpick bilmine member (F1 ring stage) — exact corollary of the parent exact stage over Q(sqrt15) (rank of an integer system is field-independent).
"""run_exact_corollary.py — BILMINE-F1 exact stage (prereg sec 7): the
Q(sqrt15) corollary of the parent's stored exact invariant computation.
No new nullspace is computed: the invariance equations form a fixed INTEGER
linear system, and the rank of an integer matrix is field-independent
(largest nonzero minor), so dims over Q(sqrt15) equal the stored dims over
Q. Re-verifies E1 (G^T P G = P, all 10, exact) as a deterministic COPAIR."""
import json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bilmine_f1_lib import (leg_dir_f1, parent_dir, require_env, emit,
                            sha_file, exact_invariance_check)


def main():
    PARENT_EX = f"{parent_dir()}/work/EXACT_INVARIANTS.json"
    ex = json.load(open(PARENT_EX))
    pf = json.load(open(require_env(
        "BILMINE_PFRAME", "the exact generator-frame JSON")))
    gens = [pf["G_integer"][k] for k in pf["frame_specification"]["GENS_order"]]
    P = pf["invariant_form"]["P"]

    e1 = exact_invariance_check(P, gens, "GtSG")
    ref = {
        "dim_sym_Q": ex["W_invariant_spaces"]["GtSG"]["dim_sym"],
        "dim_anti_Q": ex["W_invariant_spaces"]["GtSG"]["dim_anti"],
        "End_Q_dim": ex["commutant_dim_EndQ"],
        "P_in_W_sym": ex["P_in_W_sym"]["GtSG"],
        "detP": ex["detP"],
        "opposite_orientation": {
            "dim_sym": ex["W_invariant_spaces"]["GSGt"]["dim_sym"],
            "dim_anti": ex["W_invariant_spaces"]["GSGt"]["dim_anti"],
        },
    }
    ok = (e1 and ref["dim_sym_Q"] == 1 and ref["dim_anti_Q"] == 0
          and ref["End_Q_dim"] == 1 and ref["P_in_W_sym"])
    out = {
        "stage": "BILMINE-F1 exact corollary (prereg sec 7)",
        "parent_exact_invariants": {"path": PARENT_EX,
                                    "sha256": sha_file(PARENT_EX)},
        "E1_GtPG_eq_P_all10_exact_recheck": bool(e1),
        "reference_dims_over_Q": ref,
        "corollary_over_Qsqrt15": {
            "W_sym_dim": 1, "W_sym_generator": "P (reference canonical)",
            "W_anti_dim": 0, "End_dim": 1,
            "argument": ("the g-invariance equations are a fixed integer "
                         "linear system; rank of an integer matrix is "
                         "field-independent (largest nonzero minor), so the "
                         "Q(sqrt15) solution spaces are the Q(sqrt15)-spans "
                         "of the reference Q-bases: W_sym tensor Q(sqrt15) = "
                         "Q(sqrt15)*P, W_anti tensor Q(sqrt15) = 0, End = "
                         "Q(sqrt15)"),
        },
        "abstract_prediction": ("any global monodromy-invariant pairing with "
                                "Z[sqrt15] coefficients is (a + b*sqrt15)*P; "
                                "a mined rank exceeding 1 in SYM or 0 in ANTI "
                                "at the global gate is an automatic failure "
                                "for the excess"),
        "consistency": bool(ok),
    }
    emit(out, f"{leg_dir_f1()}/work/EXACT_COROLLARY.json",
         os.path.abspath(__file__))
    print("E1 recheck:", e1, "| reference dims (sym,anti,End):",
          ref["dim_sym_Q"], ref["dim_anti_Q"], ref["End_Q_dim"],
          "| corollary consistent:", ok)
    return 0 if ok else 2


if __name__ == "__main__":
    sys.exit(main())
