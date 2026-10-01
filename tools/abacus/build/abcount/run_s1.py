"""S1 driver — stage S1 of the sha-pinned spec (M3T_REGISTRATION.md sec 4
row S1):
products — battery member B1 (block-diagonal genus-4 tau);
ONE acb_theta genus-4 call cost measured FIRST (timed pilot, one member, one
prime = 13); independent C7 tail-bound cross-check on that one call;
then the member sweep (B1 at the pre-committed pilot primes {13, 31}).

Walls (S1 row, declared BEFORE launch): <= 3600 s wallclock, <= 8 GB per
member-prime; working dps declared 400 cap (the S1 row carries no dps wall;
member arithmetic runs at dps 60 as contracted, pilot desk sums boost
locally).  STOP semantics: wall breach = STOP-WALL receipt (no silent
retry); tail-bound dominance violation = STOP; any cross-route integer
disagreement = STOP.

Seeds: battery/SEEDS.json (sha-pinned, committed before any run).
Bindings riding every output: Q3c/Q3d, Q6 (single conversion site =
abcount_s0.convert_sign_convention; NO second site exists in S1 code).
"""

import hashlib
import json
import math
import os
import sys
import time
import traceback

from flint import acb, arb, ctx, fmpz

import abcount_s0 as ab
import abcount_s1 as s1

HERE = os.path.dirname(os.path.abspath(__file__))
PKG_ROOT = os.path.dirname(os.path.dirname(HERE))

PINNED = {
    "M3T_REGISTRATION.md": "85ed68777d171378c70c4b590ae35c25abb768431274bf00af6e2ed959b7bdba",
    "M3T_VALIDATION_BATTERY.md": "af97c604704f64e74194d64ddcd05fe5d5076f0fd999511b4875194c980b0d0c",
    "S0PRE_RECEIPT.md": "a23cb2ebe8de4b0819b08c4ffb9e9f59c9f891c5ef6ef99ede56c71154e9ec2e",
    "build/abcount/abcount_s0.py": "448dad948b9de15f4075097fe6b4c669e2cfc928e4836bec3cbe86307bb013cf",
    "build/abcount/battery/SEEDS.json": "6f823872744f36a7be49d03257f3dffa31e2f34a780be5241edb09ca24c96751",
}

WALL_SECONDS = 3600
WALL_MEM_GB = 8.0
WALL_DPS_CAP = 400
MEMBER_DPS = 60
PILOT_N = 12           # SEEDS S1_C7_pilot.desk_truncation_N
MEMBER_PREC_BITS = 250
RADIUS_DISCIPLINE_GATE = "1e-30"


def sha256_file(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def main():
    ctx.dps = MEMBER_DPS
    wall = ab.Wall(seconds=WALL_SECONDS, mem_gb=WALL_MEM_GB, dps=WALL_DPS_CAP)
    checks = []
    banked = {"stage": "S1", "checks": checks}

    def bank(name, verdict, detail):
        checks.append({"check": name, "verdict": verdict, "detail": detail})
        print("[%s] %s" % (verdict, name))
        with open(os.path.join(HERE, "S1_RUN_OUTPUT.json"), "w") as f:
            json.dump(banked, f, indent=1, default=str)

    # ---- check 1: sha pins + bindings + seeds sha (before anything runs) ----
    sha_rows = {}
    ok = True
    for rel, want in PINNED.items():
        got = sha256_file(os.path.join(PKG_ROOT, rel))
        sha_rows[rel] = {"got": got, "want": want, "match": got == want}
        ok = ok and got == want
    bank("sha-pin-seeds-binding-verify", "OK" if ok else "FAIL-SHA-PIN-VERIFY",
         sha_rows)
    if not ok:
        return finish(banked, wall, "S1-FAIL")

    seeds = json.load(open(os.path.join(HERE, "battery", "SEEDS.json")))
    b1 = seeds["B1"]
    curves = b1["curves"]
    primes = b1["pilot_primes"]
    assert primes == [13, 31]

    # ---- check 2: input contract + model invariants vs banked seeds ---------
    invs = []
    rows = []
    ok = True
    for c in curves:
        inv = ab.weierstrass_invariants(*c["weierstrass"])
        invs.append(inv)
        match = int(inv["disc"]) == c["disc_banked"]
        ok = ok and match
        rows.append({"label": c["label"], "model": c["weierstrass"],
                     "disc": int(inv["disc"]), "disc_banked": c["disc_banked"],
                     "match": match,
                     "j": "%d/%d" % (inv["j_num"], inv["j_den"])})
    js = set((r["j"] for r in rows))
    ok = ok and len(js) == 4
    contract_ok = (b1["sign_convention"]["s1"] == ab.GLOBAL_PIN_S1
                   and b1["polarization_elementary_divisors"] == [1, 1, 1, 1])
    bank("input-contract-and-models",
         "OK" if (ok and contract_ok) else "FAIL-INPUT-CONTRACT",
         {"sign_convention": b1["sign_convention"],
          "base_field": "Q",
          "polarization_elementary_divisors": [1, 1, 1, 1],
          "product_polarization": "principal (product of four (1)s); E = J_8 on the product basis",
          "curves": rows, "four_distinct_j": len(js) == 4,
          "pairwise_non_isogenous": "pinned member declaration (battery row B1)"})
    if not (ok and contract_ok):
        return finish(banked, wall, "S1-FAIL")

    # ---- check 3: per-block tau in H_1 + genus-1 theta receipts -------------
    taus = []
    g1_thetas = []
    rows = []
    ok = True
    for c in curves:
        w1, w2, prec_receipt = ab.pari_periods(tuple(c["weierstrass"]), MEMBER_DPS)
        tau, tau_receipt = ab.tau_from_periods(w1, w2)
        taus.append(tau)
        t2a, t3a, t4a = ab.theta_nulls_arb(tau)
        t2d, t3d, t4d, tail_receipt = ab.theta_nulls_desk(tau, N=PILOT_N)
        g1_thetas.append((t2a, t3a, t4a))
        overlaps = [ab.balls_overlap(x, y) for x, y in
                    ((t2a, t2d), (t3a, t3d), (t4a, t4d))]
        jac_ok = bool((t3a ** 4 - t2a ** 4 - t4a ** 4).contains(acb(0)))
        blk_ok = all(overlaps) and jac_ok and bool(tau.imag > 0)
        ok = ok and blk_ok
        rows.append({"label": c["label"], "tau_mid": str(complex(tau)),
                     "im_tau_positive": bool(tau.imag > 0),
                     "arb_vs_desk_overlap": overlaps,
                     "jacobi_identity_contains_0": jac_ok,
                     "per_entry_radius": prec_receipt["per_entry_radius"],
                     "sl2": tau_receipt["sl2_matrix_abcd"]})
    bank("blocks-tau-H1-genus1-theta", "OK" if ok else "STOP", rows)
    if not ok:
        return finish(banked, wall, "STOP")
    wall.check("post-blocks")

    # ---- check 4: TIMED PILOT — the ONE genus-4 acb_theta call, C7 gate ----
    # pilot tau = exact 64-bit dyadic truncations (SEEDS pilot_input_rule)
    pilot_entries = []
    pilot_taus = []
    for tau in taus:
        rm, re_e = s1.dyadic_truncate_64(tau.real)
        im, im_e = s1.dyadic_truncate_64(tau.imag)
        pilot_entries.append({"re_man": rm, "re_exp": re_e,
                              "im_man": im, "im_exp": im_e,
                              "rad_man": 0, "rad_exp": 0})
        pilot_taus.append(acb(s1._arb_exact_from_man_exp(rm, re_e),
                              s1._arb_exact_from_man_exp(im, im_e)))
    lam_lb = min((t.imag for t in pilot_taus), key=lambda x: float(x)).lower()
    tb, tb_receipt = s1.tail_bound_c7(4, PILOT_N, lam_lb)
    matched_prec = int(math.ceil(-math.log2(float(tb)))) + 64
    jres = s1.julia_theta_all_g4(pilot_entries, matched_prec, HERE, "pilot")
    wall.check("post-pilot-call")

    # echo receipt: the wrapper worked on EXACTLY the pilot tau (no rounding)
    echo_ok = True
    for k, e in enumerate(jres["tau_echo"]):
        sent = pilot_entries[k]
        same = (s1._arb_exact_from_man_exp(e[1], e[2])
                == s1._arb_exact_from_man_exp(sent["re_man"], sent["re_exp"])) \
            and (s1._arb_exact_from_man_exp(e[3], e[4])
                 == s1._arb_exact_from_man_exp(sent["im_man"], sent["im_exp"])) \
            and e[5] == 0 and e[7] == 0
        echo_ok = echo_ok and bool(same)

    # desk side at matched depth
    desk_dps = int(matched_prec * 0.302) + 20
    old_dps = ctx.dps
    ctx.dps = desk_dps
    try:
        desk_vals = s1.desk_theta_g4_all(pilot_taus, PILOT_N)
    finally:
        ctx.dps = old_dps
    tb_arb = arb(tb)
    pad = arb("0 +/- 1") * tb_arb
    dominance_viol = [i for i in range(256)
                      if not bool(jres["claimed_radii"][i] <= tb_arb)]
    overlap_fail = []
    for i in range(256):
        desk_encl = desk_vals[i] + acb(pad, pad)
        if not ab.balls_overlap(desk_encl, jres["values"][i]):
            overlap_fail.append(i)
    max_claimed = max(jres["claimed_radii"], key=lambda r: float(r))
    c7_ok = echo_ok and not dominance_viol and not overlap_fail
    bank("C7-pilot-tail-bound-cross-check",
         "OK" if c7_ok else "STOP",
         {"pilot_input": "exact 64-bit dyadic truncation of member tau mids (SEEDS S1_C7_pilot)",
          "tail_bound_receipt": tb_receipt,
          "matched_truncation": {"desk_N": PILOT_N,
                                 "matched_prec_bits": matched_prec,
                                 "rule": "prec = ceil(-log2 TailBound(4,N,lambda)) + 64 guard"},
          "pinned_criterion": "acb_theta claimed radius <= desk tail bound at matched truncation (review item C7; ANY violation = STOP)",
          "claimed_radius_max": max_claimed.str(8),
          "tail_bound": tb_arb.str(8),
          "dominance_violations": dominance_viol,
          "desk_vs_wrapper_overlap_failures": overlap_fail,
          "tau_echo_exact": echo_ok,
          "pilot_cost": {"call_seconds": jres["call_seconds"],
                         "load_seconds": jres["load_seconds"],
                         "subprocess_wall_seconds": jres["subprocess_wall_seconds"],
                         "maxrss_gb": jres["maxrss_gb"],
                         "prec_bits": matched_prec}})
    if not c7_ok:
        return finish(banked, wall, "STOP")

    # ---- check 5: pilot extrapolation (price the sweep BEFORE running it) --
    remaining_calls = 1     # one member-run ball call (theta is p-independent)
    est = jres["subprocess_wall_seconds"] * remaining_calls + 60.0
    est_ok = est < WALL_SECONDS and jres["maxrss_gb"] < WALL_MEM_GB
    bank("pilot-sweep-pricing", "OK" if est_ok else "STOP-WALL",
         {"pilot_actuals_bank_for_S2": {
             "genus4_theta_call_seconds_at_%d_bits" % matched_prec: jres["call_seconds"],
             "julia_load_seconds": jres["load_seconds"],
             "maxrss_gb": jres["maxrss_gb"]},
          "member_sweep_estimate_seconds": est,
          "wall_seconds": WALL_SECONDS,
          "note": "member-run call at %d bits is cheaper than the pilot; estimate uses pilot cost as upper bound" % MEMBER_PREC_BITS})
    if not est_ok:
        return finish(banked, wall, "STOP-WALL")

    # ---- check 6: member-run genus-4 theta receipts (route T; no Frobenius) -
    member_entries = []
    for tau in taus:
        rm, re_e = s1._man_exp_of_mid(tau.real)
        im, im_e = s1._man_exp_of_mid(tau.imag)
        rad = tau.real.rad().abs_upper().max(tau.imag.rad().abs_upper())
        rman, rexp = (rad.man_exp() if not rad.is_zero() else (0, 0))
        member_entries.append({"re_man": rm, "re_exp": re_e,
                               "im_man": im, "im_exp": im_e,
                               "rad_man": int(rman), "rad_exp": int(rexp)})
    jmem = s1.julia_theta_all_g4(member_entries, MEMBER_PREC_BITS, HERE, "member")
    wall.check("post-member-call")

    odd_fail, prod_fail = [], []
    rad_gate = arb(RADIUS_DISCIPLINE_GATE)
    rad_fail = [i for i in range(256)
                if not bool(jmem["claimed_radii"][i] < rad_gate)]
    for i in range(256):
        a, b = s1.char_decode(i)
        prod = acb(1)
        for k in range(4):
            t2a, t3a, t4a = g1_thetas[k]
            prod = prod * s1.genus1_theta_char(a[k], b[k], t2a, t3a, t4a)
        if not s1.char_is_even(i):
            if not bool(jmem["values"][i].contains(acb(0))):
                odd_fail.append(i)
        if not ab.balls_overlap(prod, jmem["values"][i]):
            prod_fail.append(i)
    pos_ok = all(bool(t.imag > 0) for t in taus)
    ok = not odd_fail and not prod_fail and not rad_fail and pos_ok
    bank("member-theta-receipts-route-T", "OK" if ok else "STOP",
         {"odd_characteristic_vanishing": {"n_odd": 120, "failures": odd_fail},
          "product_factorization_overlap_256": {"failures": prod_fail},
          "radius_discipline": {"gate": RADIUS_DISCIPLINE_GATE,
                                "max_claimed_radius": max(jmem["claimed_radii"], key=lambda r: float(r)).str(8),
                                "failures": rad_fail},
          "positivity_Im_tau_block_diagonal": pos_ok,
          "riemann_relation_receipts": "odd-char vanishing + product factorization + per-block Jacobi identity (check 3)",
          "route_T_scope": "theta receipts ONLY; route T emits NO Frobenius data (review item C2)",
          "member_call_cost": {"call_seconds": jmem["call_seconds"],
                               "maxrss_gb": jmem["maxrss_gb"],
                               "prec_bits": MEMBER_PREC_BITS}})
    if not ok:
        return finish(banked, wall, "STOP")

    # ---- check 7: pslq_gate recognition positive control (held-out gate) ----
    sys.path.insert(0, os.path.abspath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..", "..", "..", "lockpick")))
    import pslq_gate
    import mpmath as mp
    old_dps = ctx.dps
    ctx.dps = 235
    try:
        w1, w2, _ = ab.pari_periods(tuple(curves[0]["weierstrass"]), 235)
        tau_hi, _ = ab.tau_from_periods(w1, w2)
        t2h, t3h, _ = ab.theta_nulls_arb(tau_hi)
        j_ball = ab.j_from_theta(t2h, t3h)
        j_str = j_ball.real.mid().str(228, radius=False)
    finally:
        ctx.dps = old_dps
    with mp.workdps(240):
        pool = {"1": mp.nstr(mp.mpf(1), 228, strip_zeros=False),
                "pi": mp.nstr(mp.pi, 228, strip_zeros=False)}
    audit = pslq_gate.members_audit(pool, ["1", "pi"])
    proto = dict(pslq_gate.DEFAULT, maxcoeff=10 ** 9)
    ctl = pslq_gate.controls(pool, proto)
    vec = pslq_gate.two_prec_stable(j_str, ["1"], pool,
                                    dps_pair=(120, 200), maxcoeff=10 ** 9)
    inv0 = invs[0]
    jn, jd = inv0["j_num"], inv0["j_den"]
    if jd < 0:
        jn, jd = -jn, -jd
    g = math.gcd(abs(jn), jd)
    expected = (jd // g, -jn // g) if jd // g > 0 else (-jd // g, jn // g)
    rev = pslq_gate.reverify({"target": j_str, "members": ["1"],
                              "vector": list(vec)} if vec else
                             {"target": j_str, "members": ["1"], "vector": [0, 0]},
                             pool)
    ok = (audit == [] and ctl["ok"] and vec == expected and rev["ok"])
    bank("pslq-gate-j-recognition-control", "OK" if ok else "FAIL-RECOGNITION-CONTROL",
         {"block": curves[0]["label"],
          "held_out": "expected vector computed from the ACTUAL input model's exact j AFTER the scan",
          "members_audit": audit, "controls": ctl["detail"],
          "vector_got": vec, "vector_expected": expected,
          "exact_j": "%d/%d" % (inv0["j_num"], inv0["j_den"]),
          "reverify": rev,
          "discipline": "recognition is a held-out gate, never an L-factor emitter (route T re-scope)"})
    if not ok:
        return finish(banked, wall, "S1-FAIL")
    wall.check("post-pslq")

    # ---- checks 8..9: B1 integer routes at the pre-committed primes ---------
    truth13 = b1["banked_desk_truths"]["a_13"]
    truth31 = b1["banked_desk_truths"]["a_31"]
    pinned_counts = {13: b1["banked_desk_truths"]["count_F13"],
                     31: b1["banked_desk_truths"]["count_F31_derived"]}
    count_provenance = {13: "pinned ground truth (SEEDS.json)",
                        31: "derived at seed time from the pinned a_31 list"}
    for p in primes:
        truth = truth13 if p == 13 else truth31
        path1, path2, factor_rows = [], [], []
        agree = True
        for c in curves:
            a_ellap = ab.ap_pari(tuple(c["weierstrass"]), p)
            a_desk, npts = ab.ap_naive(*c["weierstrass"], p)
            banked_a = truth[c["label"]]
            row_ok = (a_ellap == a_desk == banked_a)
            agree = agree and row_ok
            path1.append(a_ellap)
            path2.append(a_desk)
            factor_rows.append({"label": c["label"], "ellap": a_ellap,
                                "desk_naive": a_desk, "banked": banked_a,
                                "n_points": npts, "exact_match": row_ok})
        # route-D assembly (abcount-internal, NO PARI in this layer), from the
        # DESK per-factor inputs (path 2)
        c2 = s1.assemble_lp_deg8(path2, p)
        c1 = s1.assemble_lp_deg8(path1, p)
        fe_ok, fe_rows = s1.functional_equation_receipt(c2, p)
        verdict, tuples, box_receipt = s1.weil_box_isolate_deg8(
            [acb(v) for v in c2[1:5]], p)
        count2 = s1.count_from_charpoly(c2)
        count1 = s1.count_from_charpoly(c1)
        prod_receipt = math.prod(p + 1 - a for a in path2)
        ok = (agree and c1 == c2 and fe_ok and verdict == "OK"
              and list(tuples[0]) == c2[1:5]
              and count1 == count2 == prod_receipt == pinned_counts[p])
        bank("B1-integer-routes-p%d" % p, "OK" if ok else "STOP",
             {"per_factor": factor_rows,
              "path1": "PARI ellap per factor (planted truth; padic_counters row)",
              "path2": "in-house naive desk count per factor (no PARI) -> abcount route-D assembly",
              "assembly": {"L_coeffs_c0_c8_from_desk_inputs": c2,
                           "L_coeffs_from_ellap_inputs": c1,
                           "assembled_equal": c1 == c2,
                           "assembly_layer": "exact integer polynomial product, abcount-internal, NO PARI call",
                           "kernel_receipt": "trivial-product (identity isogeny; route-D trivial split)"},
              "functional_equation_a8k_eq_p4k_ak": fe_rows,
              "weil_box": {"verdict": verdict,
                           "unique_tuple_a1_a4": list(tuples[0]) if tuples else None,
                           "receipt": box_receipt,
                           "enclosure_source": "integer-exact route D = zero-radius (registration C3)"},
              "count_B1_Fp": {"from_desk_assembly_P(1)": count2,
                              "from_ellap_assembly_P(1)": count1,
                              "product_receipt_prod(p+1-a)": prod_receipt,
                              "pinned_or_derived_truth": pinned_counts[p],
                              "truth_provenance": count_provenance[p],
                              "agreement": "EXACT-INTEGER"}})
        if not ok:
            return finish(banked, wall, "STOP")
        wall.check("post-p%d" % p)

    return finish(banked, wall, "S1-PASS")


def finish(banked, wall, verdict):
    w = wall.check("finish")
    checks = banked["checks"]
    n_ok = sum(1 for c in checks if c["verdict"] == "OK")
    banked.update({
        "tool": "abcount (S1 products extension over the S0 core)",
        "verdict": verdict,
        "checks_ok": "%d/%d" % (n_ok, len(checks)),
        "wall_declared": {"seconds": WALL_SECONDS, "mem_gb": WALL_MEM_GB,
                          "dps_cap_declared": WALL_DPS_CAP},
        "wall_actuals": w,
        "sign_convention": {"s1": "P", "s2": "N/A-B1 (carried, never dropped)",
                            "conversion_site": "abcount_s0.convert_sign_convention (single site, Q6; NO second site in S1 code)"},
        "input_invariants_Q3c": {
            "dimension": 4,
            "member": "B1 = 11a1 x 14a1 x 15a1 x 17a1 (product; base_field Q)",
            "polarization_elementary_divisors": [1, 1, 1, 1],
            "discs": [-161051, -21952, 50625, -83521],
            "isogeny_declared": "NONE (product member; route-D split is the trivial identity)",
        },
        "pinned_lattice_claim": "NONE (Q3d — no identification against the pinned P is made or implied)",
        "route_scope": "route T = theta receipts only (NO Frobenius data); route D = integer assembly + Weil box; route C = per-factor ellap (planted truth side)",
    })
    with open(os.path.join(HERE, "S1_RUN_OUTPUT.json"), "w") as f:
        json.dump(banked, f, indent=1, default=str)
    print("VERDICT:", verdict, banked["checks_ok"], "wall:", w)
    return 0 if verdict == "S1-PASS" else 1


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
