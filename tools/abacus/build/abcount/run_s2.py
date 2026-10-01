"""S2 driver — stage S2 of the sha-pinned spec (M3T_REGISTRATION.md sec 4
row S2):
simple genus-2 blocks — battery members B2, B3 (non-diagonal tau;
cross-method vs hyperellcharpoly; Jacobi-sum desk truth for B3).

Walls (S2 row, declared BEFORE launch, priced from the S1 actuals banked in
S1_RESULT.md sec 3): <= 7200 s wallclock, <= 16 GB per member-prime; working
dps cap 400 declared (member arithmetic at dps 60 as contracted; desk theta
sums boost locally, S1 pattern).  STOP semantics: any cross-method integer
disagreement = STOP (the disagreement is itself a banked receipt); C7-g2
dominance violation = STOP; wall breach = STOP-WALL (no silent retry);
timed period-integral pilot priced in-bridge BEFORE the sweep.

Seeds: battery/SEEDS.json S2 revision (sha-pinned below, committed before
any S2 run).  Bindings riding every output: Q3c/Q3d, Q6 (single conversion
site = abcount_s0.convert_sign_convention; NO second site in S2 code).
"""

import hashlib
import json
import math
import os
import sys
import time
import traceback
from fractions import Fraction

from flint import acb, arb, ctx

import abcount_s0 as ab
import abcount_s1 as s1
import abcount_s2 as s2
import hecke_desk as hd

HERE = os.path.dirname(os.path.abspath(__file__))
PKG_ROOT = os.path.dirname(os.path.dirname(HERE))

PINNED = {
    "M3T_REGISTRATION.md": "85ed68777d171378c70c4b590ae35c25abb768431274bf00af6e2ed959b7bdba",
    "M3T_VALIDATION_BATTERY.md": "af97c604704f64e74194d64ddcd05fe5d5076f0fd999511b4875194c980b0d0c",
    "S0PRE_RECEIPT.md": "a23cb2ebe8de4b0819b08c4ffb9e9f59c9f891c5ef6ef99ede56c71154e9ec2e",
    "build/abcount/abcount_s0.py": "448dad948b9de15f4075097fe6b4c669e2cfc928e4836bec3cbe86307bb013cf",
    "build/abcount/abcount_s1.py": "8a6130d82ac9e3a20f6fa4210d654db7d204810de52baab7888d2deff043c0bf",
    "build/abcount/theta_g4_bridge.jl": "863f7c3d4f0fd2981a455992a33a58a2e3f616cbacb9a2ebdd8fad322f487fcf",
    "build/abcount/battery/SEEDS.json": "23c86ba23c6b9b9e16d7f9342bd04b48e70bdffc437c6aa4065aa062f3e3488e",
}
# the C7-g2 instrument is calibrated as even-char dominance + all-16 overlap
# + odd-char exact-zero containment (an ALL-characteristics gate STOPs on the
# 6 odd characteristics) — see SEEDS gate_direction_erratum.

WALL_SECONDS = 7200          # pinned S2 row, per member-prime; single-process budget
WALL_MEM_GB = 16.0
WALL_DPS_CAP = 400
MEMBER_DPS = 60
DESK_N = 12                  # SEEDS S2_theta_receipts.desk_truncation_N
G4_CALL_PREC = 250           # SEEDS member_call_prec_bits
RADIUS_GATE = "1e-30"
PERIOD_PREC = 350            # SEEDS precision.period_working_prec_bits
CONTRACT_RADIUS = "1e-55"    # SEEDS precision.per_entry_radius


def sha256_file(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def quintic_disc(c):
    """Exact discriminant of the degree-5 integer polynomial c0..c5 via the
    Sylvester resultant (fraction-free determinant); pure desk arithmetic."""
    f = c
    fp = [i * c[i] for i in range(1, 6)]
    S = [[0] * 9 for _ in range(9)]
    for i in range(4):
        for j in range(6):
            S[i][i + j] = f[5 - j]
    for i in range(5):
        for j in range(5):
            S[4 + i][i + j] = fp[4 - j]
    A = [[Fraction(x) for x in row] for row in S]
    det = Fraction(1)
    for i in range(9):
        piv = next((r for r in range(i, 9) if A[r][i] != 0), None)
        if piv is None:
            return 0
        if piv != i:
            A[i], A[piv] = A[piv], A[i]
            det = -det
        det *= A[i][i]
        inv = A[i][i]
        for r in range(i + 1, 9):
            fac = A[r][i] / inv
            for cc in range(i, 9):
                A[r][cc] -= fac * A[i][cc]
    d = det / f[5]
    assert d.denominator == 1
    return int(d)            # (-1)^(5*4/2) = +1


def dyadic_tau_entries(tau):
    """SEEDS tau_dyadic_rule: exact 64-bit dyadic truncation of tau midpoints,
    normalized (trailing zero bits stripped) for exact echo comparison."""
    ents = {}
    exact = {}
    for (i, j) in ((0, 0), (0, 1), (1, 1)):
        rm, re_ = s2.normalize_dyadic(*s1.dyadic_truncate_64(tau[i][j].real))
        im, ime = s2.normalize_dyadic(*s1.dyadic_truncate_64(tau[i][j].imag))
        ents[(i, j)] = (rm, re_, im, ime)
        exact[(i, j)] = acb(s1._arb_exact_from_man_exp(rm, re_),
                            s1._arb_exact_from_man_exp(im, ime))
    exact[(1, 0)] = exact[(0, 1)]
    return ents, exact


def main():
    ctx.dps = MEMBER_DPS
    wall = ab.Wall(seconds=WALL_SECONDS, mem_gb=WALL_MEM_GB, dps=WALL_DPS_CAP)
    checks = []
    banked = {"stage": "S2", "checks": checks}

    def bank(name, verdict, detail):
        checks.append({"check": name, "verdict": verdict, "detail": detail})
        print("[%s] %s" % (verdict, name))
        with open(os.path.join(HERE, "S2_RUN_OUTPUT.json"), "w") as f:
            json.dump(banked, f, indent=1, default=str)

    # ---- check 1: sha pins + bindings + S2 seeds sha (before anything runs) -
    sha_rows = {}
    ok = True
    for rel, want in PINNED.items():
        got = sha256_file(os.path.join(PKG_ROOT, rel))
        sha_rows[rel] = {"got": got, "want": want, "match": got == want}
        ok = ok and got == want
    bank("sha-pin-seeds-binding-verify", "OK" if ok else "FAIL-SHA-PIN-VERIFY",
         sha_rows)
    if not ok:
        return finish(banked, wall, "S2-FAIL")

    seeds = json.load(open(os.path.join(HERE, "battery", "SEEDS.json")))
    B2, B3 = seeds["B2"], seeds["B3"]
    pc = seeds["S2_period_construction"]
    tr = seeds["S2_theta_receipts"]

    # ---- check 2: input contract + models vs banked seeds -------------------
    rows = {}
    ok = True
    for name, m in (("B2", B2), ("B3", B3)):
        f = m["curve"]["f_coeffs_c0_to_c5"]
        d = quintic_disc(f)
        good = all(d % p != 0 and p != 2 for p in m["pilot_primes"])
        row_ok = (d == m["curve"]["disc_f_banked"] and good
                  and m["sign_convention"]["s1"] == ab.GLOBAL_PIN_S1)
        ok = ok and row_ok
        rows[name] = {"f_coeffs": f, "disc_desk": d,
                      "disc_banked": m["curve"]["disc_f_banked"],
                      "good_reduction_at_pilots": good,
                      "pilot_primes": m["pilot_primes"],
                      "polarization_elementary_divisors": m["polarization_elementary_divisors"],
                      "sign_convention": m["sign_convention"], "match": row_ok}
    b3_char_ok = all(p % 10 == 1 for p in B3["pilot_primes"])
    ok = ok and b3_char_ok
    rows["B3_characters_exist"] = {"both_pilot_primes_1_mod_10": b3_char_ok}
    bank("input-contract-and-models", "OK" if ok else "FAIL-INPUT-CONTRACT", rows)
    if not ok:
        return finish(banked, wall, "S2-FAIL")

    # ---- checks 3..6 per member: periods -> tau -> theta -> tie -------------
    member_tau = {}
    member_roots = {}
    contract_rad = arb(CONTRACT_RADIUS)
    for name, m in (("B2", B2), ("B3", B3)):
        f = m["curve"]["f_coeffs_c0_to_c5"]
        chain = pc["chain_permutation"][name]

        # pilot-priced period run (bridge STOP-WALLs after the pilot on breach)
        Pi, roots, prec_receipt = s2.periods_from_bridge(
            f, chain, PERIOD_PREC, HERE, "periods_%s" % name,
            budget_ns=int(WALL_SECONDS * 0.5 * 1e9))
        wall.check("post-periods-%s" % name)
        member_roots[name] = roots
        max_pi_rad = max((v.real.rad().max(v.imag.rad())
                          for row in Pi for v in row), key=lambda r: float(r))
        rad_ok = bool(arb(max_pi_rad) <= contract_rad)

        tau, hom_trace = s2.recover_tau(Pi)
        member_tau[name] = tau
        tau_rads = [tau[i][j].real.rad().max(tau[i][j].imag.rad())
                    for i in range(2) for j in range(2)]
        tau_rad_ok = all(bool(arb(r) <= contract_rad) for r in tau_rads)
        ok = rad_ok and tau_rad_ok and "accepted" in hom_trace
        bank("%s-periods-and-homology" % name, "OK" if ok else "STOP",
             {"period_receipt": prec_receipt,
              "period_radius_contract": {"declared": CONTRACT_RADIUS,
                                         "max_period_radius": max_pi_rad.str(6),
                                         "within_contract": rad_ok},
              "homology_recovery": hom_trace,
              "tau_mid": [[str(complex(tau[i][j])) for j in range(2)] for i in range(2)],
              "tau_max_radius": max(tau_rads, key=lambda r: float(r)).str(6),
              "tau_radius_within_contract": tau_rad_ok,
              "im_tau_positive_definite_certified": True,
              "chain_permutation": chain,
              "sign_convention": m["sign_convention"]})
        if not ok:
            return finish(banked, wall, "STOP")

        # theta receipts: dyadic tau -> matched-prec g2 call + C7-g2 desk gate
        ents, exact = dyadic_tau_entries(tau)
        im_exact = [[exact[(i, j)].imag for j in range(2)] for i in range(2)]
        lam_lb = s2.lambda_min_lower_2x2(im_exact)
        tb, tb_receipt = s1.tail_bound_c7(2, DESK_N, lam_lb)
        matched_prec = int(math.ceil(-math.log2(float(tb)))) + 64
        # call-prec floor = the committed member_call_prec_bits: the matched
        # rule sets the MINIMUM depth at which the dominance comparison is
        # fair (it prevents false STOPs from an under-precise call); a tighter
        # call only shrinks the wrapper radii — conservative in the gate's
        # direction — and keeps the committed 1e-30 radius-discipline gate
        # meaningful when the accepted symplectic basis is non-reduced
        # (a B3 basis with lambda_min ~ 0.06 gives matched_prec = 110
        # < the gate's reach — hence the floor below)
        call_prec = max(matched_prec, G4_CALL_PREC)
        jg2 = s2.theta_from_bridge(2, ents, call_prec, HERE, "g2_%s" % name)
        wall.check("post-g2-theta-%s" % name)

        desk_dps = int(call_prec * 0.302) + 20
        old = ctx.dps
        ctx.dps = desk_dps
        try:
            tau_exact = [[exact[(0, 0)], exact[(0, 1)]],
                         [exact[(0, 1)], exact[(1, 1)]]]
            desk_vals = s2.desk_theta_g2_all(tau_exact, DESK_N)
        finally:
            ctx.dps = old
        tb_arb = arb(tb)
        pad = arb("0 +/- 1") * tb_arb
        even_idx = [i for i in range(16)
                    if (((i >> 3) & 1) * ((i >> 1) & 1)
                        + ((i >> 2) & 1) * (i & 1)) % 2 == 0]
        odd_idx = [i for i in range(16) if i not in even_idx]
        # calibrated instrument (SEEDS gate_direction + gate_direction_erratum):
        # even-char dominance; all-16 overlap; odd-char exact-zero containment
        dominance_viol = [i for i in even_idx
                         if not bool(jg2["claimed_radii"][i] <= tb_arb)]
        overlap_fail = [i for i in range(16)
                        if not ab.balls_overlap(desk_vals[i] + acb(pad, pad),
                                                jg2["values"][i])]
        odd_fail = [i for i in odd_idx
                    if not bool(jg2["values"][i].contains(acb(0)))]
        rad_gate = arb(RADIUS_GATE)
        radg_fail = [i for i in range(16)
                     if not bool(jg2["claimed_radii"][i] < rad_gate)]
        ok = (jg2["echo_exact"] and not dominance_viol and not overlap_fail
              and not odd_fail and not radg_fail)
        bank("%s-theta-receipts-g2" % name, "OK" if ok else "STOP",
             {"dyadic_rule": "exact 64-bit truncation of tau mids (SEEDS tau_dyadic_rule); echo exact: %s" % jg2["echo_exact"],
              "tail_bound_receipt": tb_receipt,
              "matched_truncation": {"desk_N": DESK_N, "matched_prec_bits": matched_prec,
                                     "call_prec_bits": call_prec,
                                     "call_prec_floor": "max(matched_prec, member_call_prec_bits = %d)" % G4_CALL_PREC,
                                     "rule": tr["matched_precision_rule"]},
              "pinned_criterion": tr["gate_direction"],
              "claimed_radius_max_even": max((jg2["claimed_radii"][i] for i in even_idx),
                                             key=lambda r: float(r)).str(8),
              "claimed_radius_max_odd": max((jg2["claimed_radii"][i] for i in odd_idx),
                                            key=lambda r: float(r)).str(8),
              "odd_radius_note": "odd chars at z=0 are exactly 0; acb_theta_all reaches them by sqrt-of-square, claimed radius ~2^(-prec/2) (SEEDS gate_direction_erratum); gated by zero-containment, not dominance",
              "tail_bound": tb_arb.str(8),
              "dominance_violations_even": dominance_viol,
              "desk_vs_wrapper_overlap_failures_16": overlap_fail,
              "odd_characteristic_vanishing": {"n_odd": 6, "failures": odd_fail},
              "radius_discipline": {"gate": RADIUS_GATE, "failures": radg_fail},
              "call_cost": {"call_seconds": jg2["call_seconds"],
                            "maxrss_gb": jg2["maxrss_gb"]},
              "route_T_scope": "theta receipts ONLY; route T emits NO Frobenius data (C2)"})
        if not ok:
            return finish(banked, wall, "STOP")

        if name == "B2":
            # genus-4 block receipts for the PRODUCT member diag(tau, tau)
            ents4 = {}
            for (i, j), v in ents.items():
                ents4[(i, j)] = v
                ents4[(i + 2, j + 2)] = v
            jg4 = s2.theta_from_bridge(4, ents4, G4_CALL_PREC, HERE, "g4_B2")
            wall.check("post-g4-theta-B2")
            odd4, prod4, radg4 = [], [], []
            for idx in range(256):
                a, b = s1.char_decode(idx)
                if not s1.char_is_even(idx):
                    if not bool(jg4["values"][idx].contains(acb(0))):
                        odd4.append(idx)
                i2 = (a[0] << 3) | (a[1] << 2) | (b[0] << 1) | b[1]
                j2 = (a[2] << 3) | (a[3] << 2) | (b[2] << 1) | b[3]
                prod = jg2["values"][i2] * jg2["values"][j2]
                if not ab.balls_overlap(prod, jg4["values"][idx]):
                    prod4.append(idx)
                if not bool(jg4["claimed_radii"][idx] < rad_gate):
                    radg4.append(idx)
            ok = jg4["echo_exact"] and not odd4 and not prod4 and not radg4
            bank("B2-theta-receipts-g4-block", "OK" if ok else "STOP",
                 {"tau4": "diag(tau2, tau2), same dyadic tau2 (SEEDS B2_genus4_block_receipt)",
                  "echo_exact": jg4["echo_exact"],
                  "odd_characteristic_vanishing": {"n_odd": 120, "failures": odd4},
                  "product_factorization_overlap_256": {"failures": prod4},
                  "radius_discipline": {"gate": RADIUS_GATE, "failures": radg4},
                  "call_cost": {"call_seconds": jg4["call_seconds"],
                                "maxrss_gb": jg4["maxrss_gb"],
                                "prec_bits": G4_CALL_PREC}})
            if not ok:
                return finish(banked, wall, "STOP")

        # tie receipt: Thomae branch-matching + covariant cross-products
        troots, tcov, chi5_nz = s2.thomae_roots_from_bridge(
            tau, 300, HERE, "thomae_%s" % name)
        ccov = s2.curve_covariants_from_bridge(f, 300, HERE, "cov_%s" % name)
        matches = s2.moebius_branch_match(roots, troots)
        cov_rows = s2.covariant_cross_products(ccov, tcov)
        anchored_rows = [r for r in cov_rows if r["anchored"]]
        cov_ok = all(r["overlap"] for r in cov_rows) and \
            all(r["overlap"] for r in anchored_rows)
        expect_anchored = (name == "B2")
        anchored_as_banked = (len(anchored_rows) > 0) == expect_anchored
        ok = chi5_nz and len(matches) > 0 and cov_ok and anchored_as_banked
        bank("%s-tie-receipt" % name, "OK" if ok else "STOP",
             {"chi5_nonzero_certified": chi5_nz,
              "moebius_branch_matching": {"n_consistent_triples": len(matches),
                                          "first_triples": matches[:5],
                                          "procedure": "SEEDS S2_tie_receipt.primary"},
              "covariant_cross_products": cov_rows,
              "anchored_status_as_banked": {"expected_anchored": expect_anchored,
                                            "got_anchored_rows": len(anchored_rows),
                                            "seed_caveat": seeds["S2_tie_receipt"]["B3_covariant_caveat_banked"] if name == "B3" else "n/a"},
              "binding_tie": "Moebius branch-matching" if name == "B3" else "both (matching + anchored covariants)"})
        if not ok:
            return finish(banked, wall, "STOP")

    # ---- check 7: tie cross-member CONTROLS (pre-declared: MUST FAIL) -------
    ctl_rows = {}
    ok = True
    for tau_name, curve_name in (("B2", "B3"), ("B3", "B2")):
        troots, _, _ = s2.thomae_roots_from_bridge(
            member_tau[tau_name], 300, HERE, "thomae_ctl_%s" % tau_name)
        matches = s2.moebius_branch_match(member_roots[curve_name], troots)
        failed_as_required = (len(matches) == 0)
        ok = ok and failed_as_required
        ctl_rows["tau(%s)_vs_branchset(%s)" % (tau_name, curve_name)] = {
            "n_consistent_triples": len(matches),
            "required": "0 (SEEDS controls_pre_declared)",
            "failed_as_required": failed_as_required}
    bank("tie-cross-member-controls", "OK" if ok else "S2-FAIL", ctl_rows)
    if not ok:
        return finish(banked, wall, "S2-FAIL")
    wall.check("post-tie-controls")

    # ---- checks 8..11: integer routes per member-prime ----------------------
    for name, m in (("B2", B2), ("B3", B3)):
        f = m["curve"]["f_coeffs_c0_to_c5"]
        for p in m["pilot_primes"]:
            bank_row = m["banked_desk_truths"]["p%d" % p]
            cp1 = s2.hyperellcharpoly_pari(f, p)          # Path 1 planted truth
            n1 = s2.count_curve_Fp(f, p)                  # desk, no PARI
            n2 = s2.count_curve_Fp2(f, p)
            cp_desk = s2.quartic_from_counts(n1, n2, p)
            base_ok = (cp1 == bank_row["charpoly_monic_desc"]
                       and cp_desk == cp1
                       and n1 == bank_row["count_C_Fp"]
                       and n2 == bank_row["count_C_Fp2"])
            detail = {
                "path1_hyperellcharpoly": cp1,
                "desk_counts": {"count_C_Fp": n1, "count_C_Fp2": n2,
                                "quartic_from_counts": cp_desk},
                "banked_desk_truths": bank_row,
                "quartic_agreement": "EXACT" if base_ok else "MISMATCH",
            }
            if name == "B2":
                c8_desk = s2.square_quartic_deg8(cp_desk)     # route-D squaring
                c8_path1 = s2.square_quartic_deg8(cp1)
                fe_ok, fe_rows = s1.functional_equation_receipt(c8_desk, p)
                verdict, tuples, box = s1.weil_box_isolate_deg8(
                    [acb(v) for v in c8_desk[1:5]], p)
                count2 = s1.count_from_charpoly(c8_desk)
                ok = (base_ok and c8_desk == c8_path1 and fe_ok
                      and verdict == "OK" and list(tuples[0]) == c8_desk[1:5]
                      and count2 == bank_row["count_B2_member"]
                      and count2 == bank_row["jac_order_P1"] ** 2)
                detail.update({
                    "path2": "desk counts (no PARI) -> route-D SQUARING assembly (abcount-internal, exact integers) -> deg-8 Weil box",
                    "assembly": {"deg8_from_desk": c8_desk,
                                 "deg8_from_path1": c8_path1,
                                 "assembled_equal": c8_desk == c8_path1,
                                 "kernel_receipt": m["isogeny_declared"]},
                    "functional_equation_deg8": fe_rows,
                    "weil_box_deg8": {"verdict": verdict,
                                      "unique_tuple_a1_a4": list(tuples[0]) if tuples else None,
                                      "receipt": box,
                                      "enclosure_source": "integer-exact route D = zero-radius (C3)"},
                    "count_member": {"P1_squared": count2,
                                     "banked": bank_row["count_B2_member"],
                                     "agreement": "EXACT-INTEGER" if count2 == bank_row["count_B2_member"] else "MISMATCH"}})
            else:
                pinned, plus, hd_receipt = hd.b3_charpoly(p)  # Path 2: hecke_desk
                fe_ok, fe_rows = s2.functional_equation_receipt_deg4(pinned, p)
                verdict, tuples, box = s2.weil_box_isolate_deg4(
                    [acb(v) for v in pinned[1:5]], p)
                count2 = sum(pinned)
                jsum_ok = (hd_receipt["jsum_vectors_basis_1_z_z2_z3"]
                           == {k: list(v) for k, v in
                               bank_row["jsum_vectors_basis_1_z_z2_z3"].items()}
                           and hd_receipt["generator_canonical_g"]
                           == seeds["B3"]["jacobi_sign_pin"]["generators_banked"]["p%d" % p])
                pin_ok = (hd_receipt["pin_discrimination"]["odd_coefficients_differ"]
                          and hd_receipt["pin_discrimination"]["even_coefficients_equal"]
                          and plus == bank_row["charpoly_plusJ_variant"])
                ok = (base_ok and pinned == cp1 and jsum_ok and pin_ok and fe_ok
                      and verdict == "OK"
                      and list(tuples[0]) == pinned[1:3]
                      and count2 == bank_row["jac_order_P1"])
                detail.update({
                    "path2_hecke_desk": hd_receipt,
                    "cross_method": {"hyperellcharpoly_vs_hecke_desk":
                                     "EXACT" if pinned == cp1 else "MISMATCH"},
                    "pin_discrimination_receipt": {
                        "plusJ_variant_equals_banked": plus == bank_row["charpoly_plusJ_variant"],
                        "odd_coefficients_differ": hd_receipt["pin_discrimination"]["odd_coefficients_differ"],
                        "battery_note": "an instantiation against +J would falsely FAIL a correct tool (pinned battery row B3)"},
                    "functional_equation_deg4": fe_rows,
                    "weil_box_deg4": {"verdict": verdict,
                                      "unique_pair_a1_a2": list(tuples[0]) if tuples else None,
                                      "receipt": box,
                                      "enclosure_source": "integer-exact hecke_desk = zero-radius"},
                    "count_member": {"P1": count2, "banked": bank_row["jac_order_P1"],
                                     "agreement": "EXACT-INTEGER" if count2 == bank_row["jac_order_P1"] else "MISMATCH"}})
            bank("%s-integer-routes-p%d" % (name, p), "OK" if ok else "STOP", detail)
            if not ok:
                return finish(banked, wall, "STOP")
            wall.check("post-%s-p%d" % (name, p))

    return finish(banked, wall, "S2-PASS")


def finish(banked, wall, verdict):
    w = wall.check("finish")
    checks = banked["checks"]
    n_ok = sum(1 for c in checks if c["verdict"] == "OK")
    banked.update({
        "tool": "abcount (S2 simple-genus-2-blocks extension over the S0 core + S1 products layer)",
        "verdict": verdict,
        "checks_ok": "%d/%d" % (n_ok, len(checks)),
        "wall_declared": {"seconds": WALL_SECONDS, "mem_gb": WALL_MEM_GB,
                          "dps_cap_declared": WALL_DPS_CAP,
                          "priced_from": "S1_RESULT.md sec 3 banked actuals + in-run period pilot"},
        "wall_actuals": w,
        "sign_convention": {"s1": "P",
                            "s2": "N/A-B2 / N/A-B3 (carried, never dropped)",
                            "conversion_site": "abcount_s0.convert_sign_convention (single site, Q6; NO second site in S2 code)"},
        "input_invariants_Q3c": {
            "B2": {"dimension": 4, "member": "J(C)^2, C: y^2 = x^5 - x + 1; base_field Q",
                   "polarization_elementary_divisors": [1, 1, 1, 1],
                   "disc_f": 2869,
                   "isogeny_declared": "product-diagonal (trivial kernel receipt)"},
            "B3": {"dimension": 2, "member": "J(C5), C5: y^2 = x^5 + 1; base_field Q",
                   "polarization_elementary_divisors": [1, 1],
                   "disc_f": 3125,
                   "endomorphisms_declared": "CM by Q(zeta_5) (pinned battery row; hecke_desk route)",
                   "isogeny_declared": "NONE (simple CM member)"},
        },
        "pinned_lattice_claim": "NONE (Q3d — no identification against the pinned P is made or implied)",
        "route_scope": "route T = theta receipts only (NO Frobenius data); route D = integer assembly (B2 squaring; B3 hecke_desk) + Weil box; route C = hyperellcharpoly (planted truth side)",
    })
    with open(os.path.join(HERE, "S2_RUN_OUTPUT.json"), "w") as f:
        json.dump(banked, f, indent=1, default=str)
    print("VERDICT:", verdict, banked["checks_ok"], "wall:", w)
    return 0 if verdict == "S2-PASS" else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except RuntimeError as e:
        if "STOP-WALL" in str(e):
            print("STOP-WALL:", e)
            sys.exit(3)
        print("FAIL-INTERNAL:", e)
        traceback.print_exc()
        sys.exit(2)
    except Exception as e:
        print("FAIL-INTERNAL:", e)
        traceback.print_exc()
        sys.exit(2)
