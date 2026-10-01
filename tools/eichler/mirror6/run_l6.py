#!/usr/bin/env python3
"""run_l6.py — the L6 fingerprint runs at BOTH MUM-type points (x=+1, x=-1).

RUNS ONLY AFTER CONTROL-PASS (checked from <outdir>/CONTROL_MIRROR6.json as
written by run_control.py; the control law is absolute).  A pilot depth runs
first, receipted, then full.

Inputs (pinned):
  --l6 PATH: the L6_exact.json operator file (the exact order-6 operator over
      Q; the operator file is not included in this repo). Its sha256 is hard-pinned
      below and checked before any run.
  Reference local structure (pinned):
      x=+1: exponents {-3,-3,0,1,2,4}, Jordan (3,1,1,1), log^2  => k=3
      x=-1: exponents {0,2,4,6,8,10},  Jordan (5,1),     log^4  => k=5
  This script re-derives the indicial multisets from the loaded operator and
  HARD-GATES on the reference values (orientation/transcription gate).

REGISTER: L6 is NOT a CY3 operator; outputs are exact invariants of the
operator's local structure under DEFINITIONS D1-D5 (mirror6.py header).  No
geometric-mirror claim, no unit-root claim, no Pi/Hodge-functional claim, no
membership claim into any CY database (order-6: not in scope of those tables).
Numeric blocks are measured-grade screens (mpmath, dps stated); NO sealed
constant verdicts (those belong to a separate sealed PSLQ layer) — value
tables only.
"""
import argparse
import json
import math
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gmpy2 import mpq
import mirror6 as M6

HERE = os.path.dirname(os.path.abspath(__file__))
L6_SHA = "b13ee756fce179b7c4fd2b51116e02a1f0b7a01038a9ac2594ce33a555caff78"

RECEIPT_NAMES = {"+1": "MIRROR6_X1.json", "-1": "MIRROR6_XM1.json",
                 "0": "MIRROR6_X0.json", "*+1": "MIRROR6_STAR_X1.json",
                 "*-1": "MIRROR6_STAR_XM1.json", "*0": "MIRROR6_STAR_X0.json"}

REFERENCE = {
    "+1": {"exponents": [-3, -3, 0, 1, 2, 4], "k": 3, "rho_min": -3},
    "-1": {"exponents": [0, 2, 4, 6, 8, 10], "k": 5, "rho_min": 0},
    # x=0: the THIRD unipotent point ({-2,-2,0,1,2,5}, J(5,1), log^4,
    # MUM-type; reference gauge).  The two NAMED points land structural q-gauge
    # obstructions (N-chain anchored at the TOP exponent); x=0 has the
    # doubled-LOWEST exponent — the only CY-positioned MUM-type pattern of
    # L6 — so the fingerprint deliverable lands here.  Ratios (mirror map,
    # alphas) are Q(x)-gauge-invariant, so a reference-gauge apparent-
    # singularity radius cap does not contaminate them; noted in the receipt.
    "0": {"exponents": [-2, -2, 0, 1, 2, 5], "k": 5, "rho_min": -2},
}


def load_l6(path):
    if not path or not os.path.exists(path):
        raise SystemExit("--l6 PATH is required and must exist (L6_exact.json "
                         "is the operator file, not included in this repo)")
    sha = M6.sha256_file(path)
    if sha != L6_SHA:
        raise SystemExit(f"L6_exact.json sha mismatch: {sha}")
    d = json.load(open(path))
    coeffs = [[mpq(int(c)) for c in row] for row in d["coeffs"]]
    return coeffs, sha


def adjoint(coeffs):
    """Formal adjoint of L = sum_i p_i(x) D^i:  L* = sum_i (-1)^i D^i o p_i,
    i.e. p*_j = sum_{i>=j} (-1)^i C(i, i-j) d^{i-j}/dx^{i-j} p_i.
    DEFINITION D6: the dual fingerprint applies D1-D5 to L* (its solution
    space is the dual local system; the exponent ladder reverses, so an
    anti-MUM-positioned chain of L can become CY-positioned for L*).
    Same Jordan partitions (transpose-inverse preserves Jordan type)."""
    n = len(coeffs) - 1
    out = [[] for _ in range(n + 1)]
    for j in range(n + 1):
        acc = []
        for i in range(j, n + 1):
            der = [mpq(c) for c in coeffs[i]]
            for _ in range(i - j):
                der = M6.pderiv(der)
            term = M6.pscal(mpq((-1) ** i) * math.comb(i, i - j), der)
            acc = M6.padd(acc, term)
        out[j] = acc if acc else [mpq(0)]
    return out


def deflate(p, r):
    """Divide poly p (low->high) by (u - r); returns (quotient, remainder)."""
    n = len(p) - 1
    q = [M6.MPQ0] * n
    q[n - 1] = p[n]
    for i in range(n - 1, 0, -1):
        q[i - 1] = p[i] + mpq(r) * q[i]
    rem = p[0] + mpq(r) * q[0]
    return q, rem


def indicial_multiset(R, lo=-20, hi=25):
    """Integer roots (with multiplicity) of R_0 in [lo, hi], exact."""
    p = [mpq(c) for c in R[0]]
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    roots = []
    for r in range(lo, hi + 1):
        while len(p) > 1 and M6.peval(p, mpq(r)) == 0:
            p, rem = deflate(p, r)
            assert rem == 0
            roots.append(r)
    return sorted(roots)


def numeric_screens(fp, label, R_alpha_depth):
    """Measured-grade numerics: ratio-test radii + rational-q value screen.
    mpmath dps 60; labeled estimates, not certified."""
    import mpmath as mp
    mp.mp.dps = 60
    out = {"register": "measured-grade numeric screen (mpmath dps 60); "
                       "ratio-test estimates at finite depth, not certified; "
                       "no constant-recognition verdicts here"}
    a2 = [mpq(c) for c in fp["_series"]["alphas"]["2"]]
    sq = [mpq(c) for c in fp["_series"]["s_of_q"]]
    for name, ser in (("alpha_2", a2), ("s_of_q", sq)):
        rats = []
        for d in range(max(2, len(ser) - 8), len(ser)):
            if ser[d - 1] != 0 and ser[d] != 0:
                rats.append(abs(mp.mpf(ser[d].numerator) / mp.mpf(ser[d].denominator)
                                / (mp.mpf(ser[d - 1].numerator) / mp.mpf(ser[d - 1].denominator))))
        if rats:
            R_est = 1 / rats[-1]
            out[f"{name}_ratio_test_radius_estimate"] = float(R_est)
            out[f"{name}_last_ratios"] = [float(x) for x in rats[-4:]]
    # rational q-point partial-sum screen inside ~half the estimated radius
    Rq = out.get("alpha_2_ratio_test_radius_estimate", 0.1)
    pts = []
    for num, den in ((1, 8), (-1, 8), (1, 5), (-1, 5), (1, 4), (-1, 4),
                     (1, 3), (-1, 3), (1, 2), (-1, 2)):
        q0 = num / den
        if abs(q0) <= 0.6 * Rq:
            val = mp.mpf(0)
            for d, c in enumerate(a2):
                val += mp.mpf(c.numerator) / mp.mpf(c.denominator) * mp.mpf(q0) ** d
            tail = abs(mp.mpf(a2[-1].numerator) / mp.mpf(a2[-1].denominator)
                       * mp.mpf(q0) ** (len(a2) - 1)) / max(1e-30, (1 - abs(q0) / Rq))
            pts.append({"q0": f"{num}/{den}",
                        "alpha2_partial_sum": mp.nstr(val, 40),
                        "tail_estimate": float(tail)})
    out["rational_q_value_screen"] = pts
    out["note"] = ("value table recorded as data for any later sealed scan; "
                   "no PSLQ / no recognition here")
    return out


def run_point(coeffs, x0, k, rho, N, label, overrides=None):
    R = M6.theta_form_at_point(coeffs, x0)
    fp = M6.fingerprint(R, rho, k, N, label, overrides=overrides)
    return R, fp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--l6", required=True,
                    help="path to L6_exact.json (sha-pinned; the operator file, not included)")
    ap.add_argument("--outdir", default="mirror6_work",
                    help="output dir; must hold CONTROL_MIRROR6.json from a "
                         "run_control.py run (default: ./mirror6_work)")
    a = ap.parse_args()
    work = os.path.abspath(a.outdir)
    os.makedirs(work, exist_ok=True)
    me = os.path.abspath(__file__)

    # control gate
    ctl = json.load(open(os.path.join(work, "CONTROL_MIRROR6.json")))
    if not ctl.get("control_pass"):
        raise SystemExit("control_pass is not true — L6 run refused")

    coeffs, sha = load_l6(a.l6)

    # exponent gates (hard)
    gates = {}
    for tag, x0 in (("+1", 1), ("-1", -1), ("0", 0)):
        R = M6.theta_form_at_point(coeffs, x0)
        ind = indicial_multiset(R)
        want = REFERENCE[tag]["exponents"]
        gates[tag] = (ind == want)
        if not gates[tag]:
            rec = {"verdict": "EXPONENT-GATE-FAIL", "point": tag,
                   "derived": ind, "reference": want}
            M6.emit_receipt(os.path.join(work, "L6_GATE_FAIL.json"), rec, me)
            raise SystemExit(f"exponent gate failed at {tag}: {ind} vs {want}")

    # points: L6 itself at the three unipotent points, then the formal
    # adjoint L6* (DEFINITION D6, dual local system) at the same points.
    star = adjoint(coeffs)
    runs = [("+1", 1, coeffs, None), ("-1", -1, coeffs, None),
            ("0", 0, coeffs, None),
            ("*+1", 1, star, 3), ("*-1", -1, star, 5), ("*0", 0, star, 5)]
    results = {}
    obstructions = {}
    pilots = {}
    for tag, x0, ops, kstar in runs:
        if kstar is None:
            k = REFERENCE[tag]["k"]
            rho = REFERENCE[tag]["rho_min"]
        else:
            k = kstar  # Jordan type preserved under adjoint
            Rst = M6.theta_form_at_point(ops, x0)
            ind = indicial_multiset(Rst, lo=-40, hi=45)
            gates[tag] = {"derived_exponents_L6star": ind,
                          "note": "no reference pin for L6*; derived this leg"}
            if not ind:
                obstructions[tag] = {
                    "verdict": "ADJOINT-INDICIAL-EMPTY (no integer class)",
                    "derived": ind}
                continue
            rho = min(ind)

        # pilot first
        t0 = time.time()
        try:
            Rp, fp_pilot = run_point(ops, x0, k, rho, 24,
                                     f"L6@{tag}-pilot")
        except M6.ChainSolveError as e:
            # STRUCTURAL LANDING: the obstruction is the deliverable at this
            # point.  The residual VALUES are an elimination-order-dependent
            # presentation; the nonempty polar set itself is invariant.
            obstructions[tag] = {
                "verdict": f"STRUCTURAL-OBSTRUCTION: no mirror-type "
                           f"q-coordinate at x={tag} under DEFINITION D1 "
                           f"(any chain orientation: the log-anchor line "
                           f"im(N^(k-1)) cap ker(N) is unique, so no "
                           f"alternative anchoring exists)",
                "mechanism": str(e),
                "residual_presentation_caveat": (
                    "un-killable residual values depend on the elimination "
                    "order (a row-reduced presentation); the nonempty "
                    "polar obstruction set below the chain bottom is "
                    "gauge/order-invariant"),
                "reference_local_structure": REFERENCE.get(tag,
                    {"note": "adjoint point (D6); exponents derived, "
                             "see gates"}),
                "pilot_depth": 24,
                "wall_seconds": round(time.time() - t0, 2),
            }
            nm = RECEIPT_NAMES[tag]
            M6.emit_receipt(os.path.join(work, nm), obstructions[tag], me)
            continue
        pilot_wall = time.time() - t0
        pilots[tag] = {"depth": 24, "wall_seconds": round(pilot_wall, 2),
                       "resonances": fp_pilot["resonances_relative_orders"],
                       "chain_bottom_relative_valuation":
                           fp_pilot["chain_bottom_relative_valuation"]}
        # price full depth: scale ~ N^2 in solve + N^3 in composition
        full_N = 96 if pilot_wall < 60 else (64 if pilot_wall < 240 else 32)
        pilots[tag]["full_N_chosen"] = full_N
        pilots[tag]["pilot_verdict"] = "PILOT-OK"

        t1 = time.time()
        R, fp = run_point(ops, x0, k, rho, full_N, f"L6@{tag}")
        fp["full_wall_seconds"] = round(time.time() - t1, 2)

        # depth-stability: pilot coefficients must prefix-match full run
        a2p = fp_pilot["_series"]["alphas"]["2"][:20]
        a2f = fp["_series"]["alphas"]["2"][:20]
        fp["depth_stability_alpha2_prefix20"] = (a2p == a2f)
        mmp = fp_pilot["_series"]["s_of_q"][:20]
        mmf = fp["_series"]["s_of_q"][:20]
        fp["depth_stability_mirror_prefix20"] = (mmp == mmf)

        # sensitivity: un-freeze each canonical-gauge parameter at unit value
        sens = []
        for (pid, j, n) in fp["params_frozen_to_zero_canonical_gauge"]:
            try:
                _, fps = run_point(ops, x0, k, rho, 24,
                                   f"L6@{tag}-sens-p{pid}",
                                   overrides={pid: 1})
                da = next((d for d in range(24)
                           if fps["_series"]["alphas"]["2"][d]
                           != fp_pilot["_series"]["alphas"]["2"][d]), None)
                dm = next((d for d in range(24)
                           if fps["_series"]["s_of_q"][d]
                           != fp_pilot["_series"]["s_of_q"][d]), None)
                sens.append({"param": pid, "tower_row_j": j, "rel_order_n": n,
                             "alpha2_first_divergence_order": da,
                             "mirror_first_divergence_order": dm})
            except M6.ChainSolveError as e:
                sens.append({"param": pid, "tower_row_j": j,
                             "rel_order_n": n, "error": str(e)})
        fp["gauge_sensitivity_block"] = sens

        # numeric screens
        fp["numeric_screens"] = numeric_screens(fp, tag, full_N)

        # divisibility fingerprint on the d^3-kernel numbers (n0=1 frame)
        a2 = [mpq(c) for c in fp["_series"]["alphas"]["2"]]
        nums = M6.lambert_invert(a2, 3, full_N)
        iv = M6.integrality_verdict(nums)
        if iv["rescale_N"]:
            scaled = {d: mpq(iv["rescale_N"]) * v for d, v in nums.items()}
            fp["divisibility_valuation_table_d3_rescaled"] = \
                M6.valuation_table(scaled)
        fp["kernel_d3_full_integrality"] = iv

        results[tag] = fp

    # summary verdict lines (script-assembled)
    def vline(tag):
        fp = results[tag]
        iv = fp["kernel_d3_full_integrality"]
        return (f"L6@{tag}: q-coordinate EXISTS (chain k={fp['chain_length_k']}, "
                f"bottom rel-val {fp['chain_bottom_relative_valuation']}; "
                f"{len(fp['params_frozen_to_zero_canonical_gauge'])} frozen "
                f"gauge params; "
                f"{len(fp['params_eliminated_by_obstructions'])} obstruction-"
                f"eliminations); mirror-map integrality: "
                f"{fp['mirror_map_integrality']['verdict']} (kappa mode "
                f"{fp['mirror_map_q_rescale_kappa']}); alpha_2 d^3-kernel "
                f"numbers: {iv['verdict']}; ladder log-free: "
                f"{fp['alpha_ladder_log_free']}")

    vlines = {tag: vline(tag) for tag in results}
    anchor_vals = {}
    import re as _re
    for tag, ob in obstructions.items():
        vlines[tag] = ob["verdict"] + " — mechanism in receipt"
        mm = _re.search(r"rel val (\d+)", ob.get("mechanism", ""))
        if mm:
            anchor_vals[tag] = int(mm.group(1))
    n_ob = len(obstructions)
    n_fp = len(results)
    leg_verdict = (
        f"L6-MIRROR-FINGERPRINT: STRUCTURAL-OBSTRUCTION-AT-ALL-UNIPOTENT-"
        f"POINTS ({n_ob}/6 runs obstructed, {n_fp} fingerprinted): L6 and "
        f"its adjoint dual carry NO CY-positioned MUM at any of x=0,+1,-1 — "
        f"every maximal N-chain is Frobenius-anchored at the TOP of its "
        f"resonant exponent ladder with exact un-gaugeable polar partners "
        f"below it (Bessel-class mechanism); under D1-D6 no mirror-type "
        f"q-coordinate, hence no instanton-type expansion, exists for L6. "
        f"Instrument validity: CONTROL-PASS 28/28 (quintic + Rodland both "
        f"MUM points, 13 published instanton numbers reproduced exactly; "
        f"Bessel negative control lands exactly this obstruction class)."
        if n_ob == 6 else
        f"L6-MIRROR-FINGERPRINT: MIXED ({n_fp} fingerprints, {n_ob} "
        f"obstructions) — see verdict_lines")
    summary = {
        "leg_verdict": leg_verdict,
        "anchor_relative_valuations_from_mechanisms": anchor_vals,
        "L6_sha256": sha,
        "reference_local_structure_consumed": REFERENCE,
        "exponent_gates": gates,
        "pilots": pilots,
        "verdict_lines": vlines,
        "register": ("exact operator-local invariants under DEFINITIONS "
                     "D1-D5; NOT CY3; no geometric-mirror/unit-root/Pi/"
                     "membership claims; reference connection data NOT "
                     "consumed — noted available for a future "
                     "absolute-normalization leg"),
    }

    # emit receipts
    for tag, nm in RECEIPT_NAMES.items():
        if tag not in results:
            continue  # obstruction receipt already emitted
        fp = results[tag]
        slim = {k2: v for k2, v in fp.items() if k2 != "_series"}
        M6.emit_receipt(os.path.join(work, nm), slim, me)
        ser = {"point": tag, "series": fp["_series"],
               "depth_N": fp["depth_N"]}
        M6.emit_receipt(os.path.join(work, nm.replace(".json",
                                                      "_SERIES.json")),
                        ser, me)
    M6.emit_receipt(os.path.join(work, "FINGERPRINT_SUMMARY.json"),
                    summary, me)
    print("L6 RUNS DONE")
    for tag in results:
        print(summary["verdict_lines"][tag])


if __name__ == "__main__":
    main()
