"""S0 driver — stage S0 of the sha-pinned spec (M3T_REGISTRATION.md sec 4
row S0):
dimension-1 sanity, ONE elliptic curve: period lattice -> tau in H_1 ->
theta nulls -> j recognition -> ellap comparison at ONE prime (p=11-class).

Instance (declared BEFORE run, from the pinned spec texts):
  curve = Cremona 11a1, model [0,-1,1,-10,-20] — a B1 factor of the pinned
  battery; prime p = 13 — the battery Scope pre-commits {13,31} with p = 11
  EXCLUDED (11a1 bad at 11) and banks the desk truth a_13(11a1) = 4,
  #E(F_13) = 10; "p=11-class" = a prime of that size class with good
  reduction, hence 13.  ONE curve, ONE prime.
Wall (S0 row): minutes wallclock (300 s declared), <= 2 GB, dps <= 60.
STOP semantics: any theta-vs-ellap mismatch or radius blow-up = STOP.
"""

import json
import sys
import time
import traceback

from flint import acb, arb, ctx

import abcount_s0 as ab

ctx.dps = 60  # wall: dps <= 60

PINNED_DESK_TRUTH_A13 = 4      # battery Scope, banked desk truth (independent)
PINNED_DESK_TRUTH_NPTS = 10    # #B1(F_13) factor for 11a1 in "10*18*16*16"

MODEL = (0, -1, 1, -10, -20)   # 11a1
PRIME = 13
DPS = 60
WALL_SECONDS = 300


def main():
    wall = ab.Wall(seconds=WALL_SECONDS, mem_gb=2.0, dps=60)
    checks = []
    def bank(name, verdict, detail):
        checks.append({"check": name, "verdict": verdict, "detail": detail})
        print("[%s] %s" % (verdict, name))

    # ---- input instance (format-B-shaped, S0 subset of the sec-6 contract) --
    w1, w2, period_receipt = ab.pari_periods(MODEL, DPS)
    inv = ab.weierstrass_invariants(*MODEL)
    instance = {
        "sign_convention": {"s1": "P", "s2": ab.S2_DIM1},
        "base_field": "Q",
        "prime": PRIME,
        "precision": {"dps": DPS, "per_entry_radius": period_receipt["per_entry_radius"]},
        "polarization_elementary_divisors": [1],
        "model_provenance": {"weierstrass": list(MODEL), "label_declared": "11a1"},
        "isogeny_receipt": None,
        "periods": (w1, w2),
        "tau": None,
    }

    # ---- check 1: input contract fields (Q6 mandatory field present) --------
    ok = (instance["sign_convention"]["s1"] in ab.SIGN_S1_VALUES
          and instance["sign_convention"]["s1"] == ab.GLOBAL_PIN_S1
          and instance["polarization_elementary_divisors"] == [1])
    bank("input-contract-fields", "OK" if ok else "FAIL-INPUT-CONTRACT",
         {"sign_convention": instance["sign_convention"],
          "base_field": "Q", "prime": PRIME,
          "polarization_elementary_divisors": [1],
          "period_import": period_receipt})
    if not ok:
        return finish(checks, wall, "S0-FAIL")

    # ---- check 2: tau in H_1 (positivity receipt) ---------------------------
    tau, tau_receipt = ab.tau_from_periods(w1, w2)
    instance["tau"] = tau
    covol = (w1.conjugate() * w2).imag
    bank("tau-in-H1-positivity", "OK", dict(tau_receipt,
         lattice_covolume_abs=str(abs(covol))[:40]))
    wall.check("post-tau")

    # ---- check 3: theta nulls, ARB route vs DESK route (proved tail) --------
    t2a, t3a, t4a = ab.theta_nulls_arb(tau)
    t2d, t3d, t4d, tail_receipt = ab.theta_nulls_desk(tau, N=12)
    overlaps = [ab.balls_overlap(x, y) for x, y in
                ((t2a, t2d), (t3a, t3d), (t4a, t4d))]
    jac = t3a ** 4 - t2a ** 4 - t4a ** 4      # Jacobi identity (dim-1 Riemann rel.)
    jac_ok = bool(jac.contains(acb(0)))
    rad_ok = all(bool(arb(t.real.rad()) < arb("1e-10")) for t in (t2a, t3a, t4a))
    ok = all(overlaps) and jac_ok and rad_ok
    bank("theta-arb-vs-desk-enclosures",
         "OK" if ok else ("STOP" if not all(overlaps) else "FAIL-THETA-RADIUS"),
         {"overlap_t2_t3_t4": overlaps, "jacobi_identity_contains_0": jac_ok,
          "arb_radii_below_1e-10": rad_ok, "desk_tail": tail_receipt,
          "theta3_arb": str(t3a)[:60]})
    if not ok:
        return finish(checks, wall, "STOP")
    wall.check("post-theta")

    # ---- check 4: j recognition against the ACTUAL input's exact j (Q3c) ----
    j_ball = ab.j_from_theta(t2a, t3a)
    jr_ok, jr_receipt = ab.j_recognition(j_ball, inv["j_num"], inv["j_den"])
    j_arb = acb(tau).modular_j()               # independent arb j receipt
    j_cross = ab.balls_overlap(j_ball / 1, j_arb)
    bank("j-recognition-actual-input",
         "OK" if (jr_ok and j_cross) else "FAIL-J-RECOGNITION",
         dict(jr_receipt, arb_modular_j_cross_overlap=j_cross,
              exact_j="%d/%d" % (inv["j_num"], inv["j_den"])))
    if not (jr_ok and j_cross):
        return finish(checks, wall, "STOP")
    wall.check("post-j")

    # ---- check 5: point count, ellap vs independent desk count --------------
    ap_c = ab.ap_pari(MODEL, PRIME)
    ap_d, npts = ab.ap_naive(*MODEL, PRIME)
    weil_ok = ap_c * ap_c <= 4 * PRIME
    ok = (ap_c == ap_d == PINNED_DESK_TRUTH_A13
          and npts == PINNED_DESK_TRUTH_NPTS and weil_ok)
    bank("count-ellap-vs-desk-p13",
         "OK" if ok else "STOP",
         {"ellap": ap_c, "desk_naive": ap_d, "n_points": npts,
          "pinned_battery_desk_truth": PINNED_DESK_TRUTH_A13,
          "pinned_battery_npts": PINNED_DESK_TRUTH_NPTS,
          "weil_bound_a2_le_4p": weil_ok, "agreement": "EXACT-INTEGER"})
    if not ok:
        return finish(checks, wall, "STOP")

    # ---- check 6: Weil-box isolation layer (enumeration semantics) ----------
    v_exact, c_exact = ab.weil_box_isolate_g1(acb(ap_c), PRIME)   # zero-radius
    wide = acb(arb("0 +/- 5"))                                    # diagnostic
    v_wide, c_wide = ab.weil_box_isolate_g1(wide, PRIME)
    ok = (v_exact == "OK" and c_exact == [PINNED_DESK_TRUTH_A13]
          and v_wide == "UNDECIDED-PRECISION")
    bank("weil-box-isolation-layer",
         "OK" if ok else "FAIL-WEIL-BOX-SEMANTICS",
         {"zero_radius_enclosure": {"verdict": v_exact, "candidates": c_exact,
                                    "note": "exact route = zero-radius, isolation trivially unique (reg. sec 1)"},
          "inflated_diagnostic": {"verdict": v_wide, "n_candidates": len(c_wide),
                                  "expected": "UNDECIDED-PRECISION (refuse, never guess)"}})
    if not ok:
        return finish(checks, wall, "S0-FAIL")

    # ---- check 7: Q6 single-site sign conversion (receipted, involutive) ----
    conv, conv_receipt = ab.convert_sign_convention(instance)
    tau_c = conv["tau"]
    j_conv = ab.j_from_theta(*ab.theta_nulls_arb(tau_c)[:2])
    jn, jd = inv["j_num"], inv["j_den"]
    if jd < 0:
        jn, jd = -jn, -jd
    jc_ok = ab.ball_contains_rational(j_conv, jn, jd)
    back, back_receipt = ab.convert_sign_convention(conv)
    ok = (conv["sign_convention"]["s1"] == "P(-1)"
          and bool(tau_c.imag > 0) and jc_ok
          and back["sign_convention"]["s1"] == "P")
    bank("sign-convention-single-site-Q6",
         "OK" if ok else "FAIL-SIGN-PIN",
         {"forward": conv_receipt, "back": back_receipt,
          "tau_conj_in_H1": bool(tau_c.imag > 0),
          "j_invariant_under_conversion_contains_exact_j": jc_ok,
          "involution_s1_restored": back["sign_convention"]["s1"] == "P"})
    if not ok:
        return finish(checks, wall, "S0-FAIL")

    return finish(checks, wall, "S0-PASS")


def finish(checks, wall, verdict):
    w = wall.check("finish")
    n_ok = sum(1 for c in checks if c["verdict"] == "OK")
    inv = ab.weierstrass_invariants(*MODEL)
    out = {
        "stage": "S0",
        "tool": "abcount (S0 dimension-1 elliptic core)",
        "verdict": verdict,
        "checks_ok": "%d/%d" % (n_ok, len(checks)),
        "wall_actuals": w,
        "sign_convention": {"s1": "P", "s2": ab.S2_DIM1,
                            "conversion_site": "abcount_s0.convert_sign_convention (single site, Q6)"},
        "input_invariants_Q3c": {
            "dimension": 1,
            "polarization_elementary_divisors": [1],
            "model": list(MODEL),
            "c4": int(inv["c4"]), "c6": int(inv["c6"]),
            "disc": int(inv["disc"]),
            "j": "%d/%d" % (inv["j_num"], inv["j_den"]),
            "base_field": "Q", "prime": PRIME,
            "isogeny_declared": "NONE",
        },
        "pinned_lattice_claim": "NONE (Q3d — no identification against the pinned P is made or implied)",
        "checks": checks,
    }
    with open("S0_RUN_OUTPUT.json", "w") as f:
        json.dump(out, f, indent=1, default=str)
    print("VERDICT:", verdict, out["checks_ok"], "wall:", w)
    return 0 if verdict == "S0-PASS" else 1


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
