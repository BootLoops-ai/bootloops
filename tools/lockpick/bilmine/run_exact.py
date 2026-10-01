#!/usr/bin/env python3
# lockpick bilmine member — exact stage: deterministic integer algebra on the certified generators (invariant bilinear spaces, commutant dimension).
"""run_exact.py — BILMINE exact stage (prereg sec 7): deterministic integer
algebra on the stored EXACT generators. E1: re-verify G^T P G = P for all 10
(COPAIR to the stored EXACT_PFRAME re-check). E2: exact invariant bilinear
spaces W_sym / W_anti (both orientations) + commutant dimension. No period
number is touched here."""
import json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bilmine_lib import (leg_dir, require_env, emit, sha_file,
                         invariant_space, commutant_dim,
                         exact_invariance_check, hnf_basis, canonicalize)


def main():
    d = json.load(open(require_env(
        "BILMINE_PFRAME", "the exact generator-frame JSON")))
    order = d["frame_specification"]["GENS_order"]
    gens = [d["G_integer"][k] for k in order]
    P = d["invariant_form"]["P"]

    # E1 — exact re-check, all 10 generators
    e1 = exact_invariance_check(P, gens, "GtSG")
    e1_alt = exact_invariance_check(P, gens, "GSGt")
    import sympy
    Pm = sympy.Matrix(P)
    detP = int(Pm.det())
    eig_sig = sum(1 for x in Pm.eigenvals(multiple=True) if x > 0)

    # E2 — exact invariant spaces (both orientations) + commutant
    W = {}
    for orient in ("GtSG", "GSGt"):
        W[orient] = {
            "sym_basis_pairsorder": invariant_space(gens, "sym", orient),
            "anti_basis_pairsorder": invariant_space(gens, "anti", orient),
        }
        W[orient]["dim_sym"] = len(W[orient]["sym_basis_pairsorder"])
        W[orient]["dim_anti"] = len(W[orient]["anti_basis_pairsorder"])
    cdim = commutant_dim(gens)

    # P membership in W (its own orientation)
    from bilmine_lib import sym_pairs
    Pvec = list(canonicalize([P[k][l] for (k, l) in sym_pairs()]))
    in_w = {}
    for orient in ("GtSG", "GSGt"):
        basis = W[orient]["sym_basis_pairsorder"]
        aug = hnf_basis(basis + [Pvec]) if basis else ([Pvec] if any(Pvec) else [])
        in_w[orient] = (len(aug) == len(basis))

    out = {
        "stage": "BILMINE exact invariants (prereg sec 7)",
        "gens_order": order,
        "E1_GtPG_eq_P_all10_exact": bool(e1),
        "E1_GPGt_eq_P_all10_exact": bool(e1_alt),
        "E1_copair_banked": bool(d["exact_checks_fresh"]["GtPG_eq_P_all10_exact"]),
        "detP": detP,
        "P_signature_plus": eig_sig,
        "W_invariant_spaces": W,
        "commutant_dim_EndQ": cdim,
        "P_in_W_sym": in_w,
        "note": ("bases are primitive HNF rows in the sym_pairs()/anti_pairs() "
                 "coordinate order of bilmine_lib"),
    }
    emit(out, f"{leg_dir()}/work/EXACT_INVARIANTS.json", os.path.abspath(__file__))
    print("E1 GtPG==P all10:", e1, "| alt orientation:", e1_alt,
          "| detP", detP, "| sig+", eig_sig)
    for orient in W:
        print(orient, "dim_sym", W[orient]["dim_sym"],
              "dim_anti", W[orient]["dim_anti"], "| P in W:", in_w[orient])
    print("commutant dim:", cdim)


if __name__ == "__main__":
    main()
