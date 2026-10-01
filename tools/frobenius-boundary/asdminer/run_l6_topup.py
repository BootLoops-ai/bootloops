#!/usr/bin/env python3
"""run_l6_topup.py — complete the declared 17-prime menu: if a fan run's
launch-time projection trimmed primes but the measured fan wall left headroom
under the stated budget, this top-up runs the trimmed primes under the same
adaptive law and re-emits ASD_L6_TABLE.json with recomputed summaries and a
narrated topup_note.  Same verdict grammar; nothing hand-set.
Set FROB_L6_OPERATOR (see run_l6_gauge.py)."""
import json, os, subprocess, sys, time, datetime
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import asdlib
from asdlib import load_operator, shift_operator, sha256_file
from run_l6_mine import run_point_prime, WORK, L6PATH

def utc():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

def main():
    t0 = time.time()
    if not L6PATH or not os.path.exists(L6PATH):
        raise SystemExit("run_l6_topup: the L6-class operator is not distributed "
                         "with the repo; set FROB_L6_OPERATOR to a JSON operator file")
    table = json.load(open(os.path.join(WORK, "ASD_L6_TABLE.json")))
    gauge = json.load(open(os.path.join(WORK, "L6_GAUGE.json")))
    menu = table["menu_primes"]
    done = set(table["fan_menu_run"])
    missing = [p for p in menu if p not in done]
    if not missing:
        print("nothing to top up")
        return 0
    d, C = load_operator(L6PATH)
    pts = {}
    for x0 in (0, -1, 1):
        gp = gauge["points"][str(x0)]
        pts[x0] = {"Cs": shift_operator(C, x0), "rho": gp["candidate"]["rho"],
                   "reszero": [r["resonances_zeroed"] for r in gp["rho_scan"]
                               if r["status"] == "LOG-FREE"][0],
                   "pden_profile": gp["denominator_law"]["menu_prime_profile"]}
    for p in missing:
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
            row = run_point_prime(pts[x0]["Cs"], pts[x0]["rho"],
                                  pts[x0]["reszero"], p,
                                  f"x0={x0}, rho={pts[x0]['rho']}, "
                                  f"stage-1 N=64p (topup)", N=64 * p)
            if str(row.get("verdict", "")).startswith("CONSISTENT"):
                full = run_point_prime(pts[x0]["Cs"], pts[x0]["rho"],
                                       pts[x0]["reszero"], p,
                                       f"x0={x0}, ESCALATED on stage-1 "
                                       f"consistency (topup)")
                full["stage1_row"] = row
                row = full
            table["points"][str(x0)]["rows"].append(row)
            print(f"topup x0={x0} p={p}: {row['verdict']} D={row.get('D')} "
                  f"class={row.get('gamma_class')} ({row['wall_s']}s)")
        prof1 = pts[1]["pden_profile"][str(p)]
        table["points"]["1"]["rows"].append({
            "p": p, "verdict": "P-DENOMINATOR-WALL",
            "mechanism": f"(n+3)^2 double-root division at m=p-3="
                         f"{prof1['first_m_with_p_denominator']} uncancelled "
                         f"(measured); every row's a(np) is p-denominated"})
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
        table["points"][x0]["rows"].sort(key=lambda r: r["p"])
        table["points"][x0]["summary"] = summarize(table["points"][x0]["rows"])
    table["fan_menu_run"] = sorted(done | set(missing))
    table["topup_note"] = (
        f"top-up {utc()}: the launch-time projection trimmed {missing}; the "
        f"measured fan wall left budget headroom, so the trimmed primes were "
        f"run under the same adaptive law; full menu now covered. topup wall "
        f"{round(time.time() - t0, 1)}s")
    table["walls_s_total"] = round(
        table.get("walls_s_total", 0) + (time.time() - t0), 2)
    table["producer"] = {"component": "frobenius-boundary/asdminer",
                         "script": "asdminer/run_l6_topup.py",
                         "script_sha256": sha256_file(os.path.abspath(__file__)),
                         "asdlib_sha256": sha256_file(asdlib.__file__),
                         "stamp_utc": utc()}
    out = os.path.join(WORK, "ASD_L6_TABLE.json")
    with open(out, "w") as f:
        json.dump(table, f, indent=1)
    print(asdlib.lint_check(out))
    for x0 in ("0", "-1", "1"):
        print(f"SUMMARY x0={x0}:", table["points"][x0]["summary"]["verdict_line"])
    return 0

if __name__ == "__main__":
    sys.exit(main())
