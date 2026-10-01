#!/usr/bin/env python3
"""run_l6_mine.py — stages B+C: pilot (p=53) then the split-prime fan on the
L6 candidate series.  Consumes work/L6_GAUGE.json (stage A).
Emits work/ASD_L6_PILOT.json then work/ASD_L6_TABLE.json.
Points/gauges: (x0=0, rho=5), (x0=-1, rho=10) mined on the modular route at
N(p)=p^2+p, S_target=4; (x0=+1, rho=0) is p-denominator-walled at every menu
prime (gauge receipt) — the pilot demonstrates the mechanical wall verdict
from the exact series; fan rows cite the measured wall.  Weight hypothesis
k=4 (rank-2 class); mined gamma mod p^2 is k-blind for k>=3.
Set FROB_L6_OPERATOR (see run_l6_gauge.py)."""
import json, os, subprocess, sys, time, datetime
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import asdlib
from asdlib import (load_operator, shift_operator, frobenius_exact,
                    frobenius_modular, mine, emp_depth_profile, sha256_file)

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.environ.get("FROB_ASDMINER_WORK") or os.path.join(HERE, "work")
L6PATH = os.environ.get("FROB_L6_OPERATOR")
K = 4
S_TARGET = 4
BUDGET_S = 3600   # wall-clock cap for the fan; trim-from-top + state if over

def utc():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

def run_point_prime(Cs, rho, reszero, p, note, N=None):
    """Modular-route mine of one (point, prime). Returns table row."""
    if N is None:
        N = p * p + p
    t = time.time()
    try:
        bm, F, W, V = frobenius_modular(Cs, rho, N, p, S_TARGET, reszero)
    except RuntimeError as e:
        return {"p": p, "N": N, "verdict": "P-DENOMINATOR-WALL",
                "wall_step": str(e), "wall_s": round(time.time() - t, 2),
                "note": note}
    aval = lambda n: bm[n - 1]
    r = mine(aval, p, N, K, want_profile=False)
    gt = r["guaranteed"]; xt = r["exploratory"]
    row = {"p": p, "N": N, "floor_F": F, "working_W": W, "V_total": V,
           "verdict": gt["verdict"], "D": gt["depth_D"],
           "rows_used": gt["rows_used"], "rows_degenerate": gt["rows_degenerate"],
           "n_witnesses": gt["n_witnesses"],
           "gamma_mod_pD": gt.get("gamma_mod_pD"),
           "gamma_class": gt.get("gamma_class"),
           "ordinary": gt.get("ordinary"), "pins": gt.get("pins"),
           "weil_box_k4": gt.get("weil_box_k"),
           "weil_unique": gt.get("weil_unique"),
           "witnesses": gt.get("witnesses", [])[:4],
           "exploratory": {"verdict": xt["verdict"], "D": xt["depth_D"],
                           "gamma_mod_pD": xt.get("gamma_mod_pD"),
                           "gamma_class": xt.get("gamma_class"),
                           "n_witnesses": xt["n_witnesses"],
                           "drift": xt.get("drift_detector")},
           "note": note, "wall_s": round(time.time() - t, 2)}
    c = asdlib.profile_candidate(gt, p, K) if str(gt["verdict"]).startswith("CONSISTENT") else None
    if c is not None:
        row["profile_candidate"] = c
        row["empirical_depth_capped_at_F"] = emp_depth_profile(
            aval, p, N, K, c, cap=F)
    return row

def main():
    t0 = time.time()
    if not L6PATH or not os.path.exists(L6PATH):
        raise SystemExit("run_l6_mine: the L6-class operator is not distributed "
                         "with the repo; set FROB_L6_OPERATOR to a JSON operator file")
    gauge = json.load(open(os.path.join(WORK, "L6_GAUGE.json")))
    menu = gauge["menu_primes_split_Qsqrt15"]
    d, C = load_operator(L6PATH)
    pts = {}
    for x0 in (0, -1, 1):
        gp = gauge["points"][str(x0)]
        rho = gp["candidate"]["rho"]
        reszero = [r["resonances_zeroed"] for r in gp["rho_scan"]
                   if r["status"] == "LOG-FREE"][0]
        pts[x0] = {"Cs": shift_operator(C, x0), "rho": rho,
                   "reszero": reszero,
                   "pden_profile": gp["denominator_law"]["menu_prime_profile"]}

    # ---------------- stage B: measured timing PILOT at p=53
    pilot = {"stamp_utc_start": utc(), "pilot_prime": 53, "k": K,
             "S_target": S_TARGET, "points": {}}
    for x0 in (0, -1):
        pilot["points"][str(x0)] = run_point_prime(
            pts[x0]["Cs"], pts[x0]["rho"], pts[x0]["reszero"], 53,
            f"x0={x0}, rho={pts[x0]['rho']}, modular route")
        pr = pilot["points"][str(x0)]
        print(f"PILOT x0={x0} p=53: {pr['verdict']} D={pr.get('D')} "
              f"gamma={pr.get('gamma_mod_pD')} class={pr.get('gamma_class')} "
              f"rows={pr.get('rows_used')} wall={pr['wall_s']}s")
    # index-gauge probe rider: the one residual mod-p
    # gauge freedom is the index convention; probe fg vs upow vs +-1 shifts
    pilot["index_gauge_probe_p53"] = {}
    for x0 in (0, -1):
        rho = pts[x0]["rho"]
        bm, F, W, V = frobenius_modular(pts[x0]["Cs"], rho, 53 * 54, 53,
                                        S_TARGET, pts[x0]["reszero"])
        probe = {}
        for name, sh in (("fg_b_nm1", 0), (f"upow_b_nm1m{rho}", rho),
                         ("shift_plus1", -1), ("shift_minus1", 1)):
            def aval(n, sh=sh):
                i = n - 1 - sh
                return bm[i] if 0 <= i < len(bm) else 0
            r = mine(aval, 53, 53 * 54, K, want_profile=False)
            gt = r["guaranteed"]
            probe[name] = {"verdict": gt["verdict"], "rows": gt["rows_used"],
                           "n_witnesses": gt["n_witnesses"]}
        pilot["index_gauge_probe_p53"][str(x0)] = probe
        print(f"index probe x0={x0}:",
              {k: v["verdict"] for k, v in probe.items()})

    # x=+1 wall demonstration from the exact series (mechanical verdict)
    b1, res1, _ = frobenius_exact(pts[1]["Cs"], pts[1]["rho"], 600)
    a1 = lambda n: b1[n - 1]
    r1 = mine(a1, 53, 600, K, want_profile=False)
    pilot["points"]["1"] = {
        "p": 53, "N": 600, "route": "exact (wall demonstration)",
        "verdict": r1["guaranteed"]["verdict"],
        "rows_used": r1["guaranteed"]["rows_used"],
        "pden_rows_head": r1["guaranteed"]["pden_rows"],
        "note": "x0=+1 rho=0: every row's a(np) hits the measured "
                "m=p-3 double-root p^2 denominator (gauge receipt); "
                "wall verdict is the mechanical outcome"}
    print("PILOT x0=+1 p=53 (exact route):", r1["guaranteed"]["verdict"])
    pilot["producer"] = {"component": "frobenius-boundary/asdminer",
                         "script": "asdminer/run_l6_mine.py",
                         "script_sha256": sha256_file(os.path.abspath(__file__)),
                         "asdlib_sha256": sha256_file(asdlib.__file__),
                         "stamp_utc": utc()}
    outp = os.path.join(WORK, "ASD_L6_PILOT.json")
    with open(outp, "w") as f:
        json.dump(pilot, f, indent=1)
    print(asdlib.lint_check(outp))

    # ---------------- stage C: FAN (budget-projected from the pilot)
    w53 = sum(pilot["points"][str(x0)]["wall_s"] for x0 in (0, -1))
    proj = sum(w53 * ((64 * p) / (53 * 53 + 53)) ** 1.7 for p in menu)
    fan_menu = list(menu)
    trim_note = None
    while proj > BUDGET_S and len(fan_menu) > 1:
        drop = fan_menu.pop()
        proj -= w53 * ((64 * drop) / (53 * 53 + 53)) ** 1.7
        trim_note = (f"fan trimmed at the wall-clock budget {BUDGET_S}s: primes above "
                     f"{fan_menu[-1]} not run (raise BUDGET_S to run wider)")
    table = {"stamp_utc_start": utc(), "k_hypothesis": K,
             "k_blindness_note": "mined gamma mod p^2 is identical for every "
                                 "k>=3; k=4 labels the rank-2 class",
             "fan_depth_law": "adaptive: stage-1 "
                              "N=64p, escalate to p^2+p on CONSISTENT*; "
                              "stage-1 D capped at 1 for p>64",
             "menu_primes": menu, "fan_menu_run": fan_menu,
             "trim_note": trim_note, "S_target": S_TARGET,
             "projection_s_from_pilot": round(proj, 1),
             "points": {"0": {"gauge": "x0=0, rho=5, u=x, formal-group index",
                              "rows": []},
                        "-1": {"gauge": "x0=-1, rho=10, u=x+1, formal-group index",
                               "rows": []},
                        "1": {"gauge": "x0=+1, rho=0, u=x-1, formal-group index",
                              "rows": []}},
             "copair_rider": "COPAIR-READY: an independent series-side "
                             "partner for any operator-side route; this run "
                             "reads nothing from such a route."}
    for p in fan_menu:
        for x0 in (0, -1):
            prof = pts[x0]["pden_profile"][str(p)]
            if prof["first_m_with_p_denominator"] is not None:
                table["points"][str(x0)]["rows"].append({
                    "p": p, "verdict": "P-DENOMINATOR-WALL",
                    "mechanism": f"measured first p-denominator at m="
                                 f"{prof['first_m_with_p_denominator']} "
                                 f"(gauge receipt), v_max_600="
                                 f"{prof['v_max_to_600']}"})
                continue
            if p == 53:
                row = dict(pilot["points"][str(x0)])
                row["note"] = (row.get("note") or "") + " [pilot row reused]"
                table["points"][str(x0)]["rows"].append(row)
                continue
            row = run_point_prime(pts[x0]["Cs"], pts[x0]["rho"],
                                  pts[x0]["reszero"], p,
                                  f"x0={x0}, rho={pts[x0]['rho']}, "
                                  f"stage-1 N=64p (adaptive law)", N=64 * p)
            if str(row.get("verdict", "")).startswith("CONSISTENT"):
                row_full = run_point_prime(pts[x0]["Cs"], pts[x0]["rho"],
                                           pts[x0]["reszero"], p,
                                           f"x0={x0}, ESCALATED to N=p^2+p "
                                           f"on stage-1 consistency")
                row_full["stage1_row"] = row
                row = row_full
            table["points"][str(x0)]["rows"].append(row)
            print(f"fan x0={x0} p={p}: {row['verdict']} D={row.get('D')} "
                  f"gamma={row.get('gamma_mod_pD')} class={row.get('gamma_class')} "
                  f"({row['wall_s']}s)")
        # x=+1 wall rows from the gauge receipt (measured, per prime)
        prof1 = pts[1]["pden_profile"][str(p)]
        table["points"]["1"]["rows"].append({
            "p": p, "verdict": "P-DENOMINATOR-WALL",
            "mechanism": f"(n+3)^2 double-root division at m=p-3="
                         f"{prof1['first_m_with_p_denominator']} uncancelled "
                         f"(measured); every row's a(np) is p-denominated"})

    # ---------------- per-point summary verdict lines (script-computed)
    def summarize(rows):
        cons = [r for r in rows if str(r.get("verdict", "")).startswith("CONSISTENT")]
        incons = [r for r in rows if r.get("verdict") == "INCONSISTENT"]
        wall = [r for r in rows if r.get("verdict") == "P-DENOMINATOR-WALL"]
        classes = sorted({str(r.get("gamma_class")) for r in cons})
        return {"n_primes": len(rows), "n_consistent": len(cons),
                "n_inconsistent": len(incons), "n_walled": len(wall),
                "consistent_classes": classes,
                "verdict_line": (
                    f"CONSISTENT at {len(cons)}/{len(rows)} fan primes "
                    f"(classes {classes}), INCONSISTENT at {len(incons)}, "
                    f"walled {len(wall)}")}
    for x0 in ("0", "-1", "1"):
        table["points"][x0]["summary"] = summarize(table["points"][x0]["rows"])
    table["walls_s_total"] = round(time.time() - t0, 2)
    table["producer"] = {"component": "frobenius-boundary/asdminer",
                         "script": "asdminer/run_l6_mine.py",
                         "script_sha256": sha256_file(os.path.abspath(__file__)),
                         "asdlib_sha256": sha256_file(asdlib.__file__),
                         "stamp_utc": utc()}
    outt = os.path.join(WORK, "ASD_L6_TABLE.json")
    with open(outt, "w") as f:
        json.dump(table, f, indent=1)
    print(asdlib.lint_check(outt))
    for x0 in ("0", "-1", "1"):
        print(f"SUMMARY x0={x0}:", table["points"][x0]["summary"]["verdict_line"])
    return 0

if __name__ == "__main__":
    sys.exit(main())
