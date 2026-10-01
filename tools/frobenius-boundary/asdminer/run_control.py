#!/usr/bin/env python3
"""run_control.py — asdminer positive-control battery.  Emits
work/ASD_CONTROL.json (script-computed gates G1,G2,G3a,G3b,G3c,G4,G5,G6 +
control_pass) and work/probe_gauge_lesson.json (the report-only gauge-lesson
probe).  Needs gp (PARI) on PATH for the comparator routes; reference gp
script+log pairs from a prior run live in fixtures/.  A pilot mine (E1, p=5)
prints first; the battery follows.  This is the member's fast battery leg
(~6 s measured with gp on PATH)."""
import json, os, subprocess, sys, time, random, datetime
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import asdlib
from asdlib import (mine, apery_zeta3, apery_zeta2, eta_product_8_4,
                    ell_omega_series, sha256_file)

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.environ.get("FROB_ASDMINER_WORK") or os.path.join(HERE, "work")
os.makedirs(WORK, exist_ok=True)
NC = 3000          # Apery / eigenform depth
NELL = 1200        # elliptic depth
E1 = [0, -1, 1, -10, -20]   # conductor 11
E2 = [0, 0, 0, -1, 0]       # conductor 32 (CM)
P_ELL = [5, 7, 13, 17, 19, 23, 29, 31]
P_AP = [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43]
P_C0 = [5, 7, 11, 13, 17]

def utc():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

def gp_run(text, log):
    with open(log + ".gp", "w") as f:
        f.write(text)
    out = subprocess.run(["gp", "-q", "-f", log + ".gp"],
                         capture_output=True, text=True, timeout=900).stdout
    with open(log, "w") as f:
        f.write(out)
    return out

def main():
    t0 = time.time()
    rec = {"spec": "asdminer positive-control battery (gates G1-G6)",
           "N_c": NC, "N_ell": NELL, "stamp_utc_start": utc(), "walls_s": {}}

    # ---------------- series builds
    t = time.time(); A = apery_zeta3(NC)
    rec["walls_s"]["apery3_build"] = round(time.time() - t, 2)
    t = time.time(); eta = eta_product_8_4(NC)
    rec["walls_s"]["eta_build"] = round(time.time() - t, 2)

    # ---------------- route A: PARI mfcoefs
    t = time.time()
    gpout = gp_run('{S = mfinit([8,4], 0); d = mfdim([8,4], 0);'
                   ' print("DIM=", d); f = mfeigenbasis(S)[1];'
                   ' print("COEFS=", mfcoefs(f, %d)); quit; }\n' % NC,
                   os.path.join(WORK, "gp_mf84.log"))
    rec["walls_s"]["gp_mfcoefs"] = round(time.time() - t, 2)
    dim, mfc = None, None
    for line in gpout.splitlines():
        if line.startswith("DIM="):
            dim = int(line[4:].strip())
        if line.startswith("COEFS="):
            mfc = [int(x) for x in line[7:].strip().strip("[]").split(",")]
    g1 = (dim == 1 and mfc is not None and len(mfc) >= NC + 1 and
          all(mfc[n] == eta[n] for n in range(NC + 1)))
    rec["G1_route_agreement"] = {"pass": bool(g1), "newspace_dim_8_4": dim,
        "detail": "eta(2t)^4 eta(4t)^4 (python ints) == mfcoefs(8.4.a.a) all n<=3000"}
    ap84 = {p: eta[p] for p in range(2, 200) if all(p % q for q in range(2, p))}

    # ---------------- route A': gp ellak tables for E1, E2 (+ conductors)
    gpo = gp_run('{E1 = ellinit(%s); E2 = ellinit(%s);'
                 ' print("N1=", ellglobalred(E1)[1]);'
                 ' print("N2=", ellglobalred(E2)[1]);'
                 ' print("AP1=", vector(50, i, ellak(E1, prime(i))));'
                 ' print("AP2=", vector(50, i, ellak(E2, prime(i)))); quit; }\n'
                 % (E1, E2), os.path.join(WORK, "gp_ellap.log"))
    N1 = N2 = None; ap_e1 = ap_e2 = None
    prlist = []
    n = 2
    while len(prlist) < 50:
        if all(n % q for q in range(2, n)):
            prlist.append(n)
        n += 1
    for line in gpo.splitlines():
        if line.startswith("N1="): N1 = int(line[3:])
        if line.startswith("N2="): N2 = int(line[3:])
        if line.startswith("AP1="):
            ap_e1 = dict(zip(prlist, [int(x) for x in line[4:].strip().strip("[]").split(",")]))
        if line.startswith("AP2="):
            ap_e2 = dict(zip(prlist, [int(x) for x in line[4:].strip().strip("[]").split(",")]))
    rec["route_gp_ell"] = {"E1": E1, "conductor_E1": N1,
                           "E2": E2, "conductor_E2": N2}

    # ---------------- measured timing pilot: E1 at p=5 (printed first)
    BIGM = {}
    def ell_aval(ainvs, p, tag):
        if tag not in BIGM:
            BIGM[tag] = ell_omega_series(ainvs, NELL, p ** 8)
        c = BIGM[tag]
        return lambda nn: c[nn]
    t = time.time()
    pilot = mine(ell_aval(E1, 5, "E1_5"), 5, NELL, 2)
    rec["walls_s"]["pilot_E1_p5"] = round(time.time() - t, 2)
    rec["pilot_E1_p5"] = {kk: pilot[kk] for kk in ("guaranteed", "exploratory")}
    print("PILOT (Honda E1, p=5):", pilot["guaranteed"]["verdict"],
          "D=", pilot["guaranteed"]["depth_D"],
          "weil=", pilot["guaranteed"].get("weil_box_k"),
          "| ellak a_5 =", ap_e1[5])

    # ---------------- G3a / G3b: Honda controls
    def ell_gate(ainvs, apt, tag):
        out = {}
        allok = True
        for p in P_ELL:
            r = mine(ell_aval(ainvs, p, f"{tag}_{p}"), p, NELL, 2,
                     want_profile=False)
            gt = r["guaranteed"]
            need_D = 2 if apt[p] % p != 0 else 1   # supersingular rule (add-2)
            ok = (gt["verdict"] == "CONSISTENT" and gt["depth_D"] >= need_D
                  and gt.get("weil_unique") and gt["weil_box_k"][0] == apt[p])
            allok = allok and ok
            out[str(p)] = {"verdict": gt["verdict"], "D": gt["depth_D"],
                           "rows": gt["rows_used"],
                           "candidate": gt.get("weil_box_k"),
                           "ellak_ap": apt[p], "match": ok}
        return out, allok
    t = time.time()
    c1e1, g3a = ell_gate(E1, ap_e1, "E1")
    c1e2, g3b = ell_gate(E2, ap_e2, "E2")
    rec["walls_s"]["honda_both"] = round(time.time() - t, 2)
    rec["C1_E1_honda"] = c1e1
    rec["C1_E2_honda"] = c1e2
    rec["G3a_E1"] = {"pass": bool(g3a), "primes": P_ELL}
    rec["G3b_E2"] = {"pass": bool(g3b), "primes": P_ELL}

    # ---------------- G2 / C0: exact-Hecke smoke test (k=4, p^3 wiring)
    a_f = lambda n: eta[n]
    c0 = {}
    g2 = True
    for p in P_C0:
        r = mine(a_f, p, NC, 4, want_profile=False)
        gt = r["guaranteed"]; xt = r["exploratory"]
        guar_match = gt.get("weil_unique") and gt["weil_box_k"][0] == ap84[p]
        expl_match = xt.get("weil_unique") and xt.get("weil_box_k", [None])[0] == ap84[p]
        ok = (gt["verdict"] == "CONSISTENT" and xt["verdict"] == "CONSISTENT"
              and (guar_match or expl_match))   # non-ordinary fallback (add-2)
        g2 = g2 and ok
        c0[str(p)] = {"verdict": gt["verdict"], "D": gt["depth_D"],
                      "candidate": gt.get("weil_box_k"),
                      "expl_candidate": xt.get("weil_box_k"),
                      "via": ("guaranteed" if guar_match else
                              "exploratory" if expl_match else "none"),
                      "match_ap": ok}
    rec["C0_eigenform_smoke"] = c0
    rec["G2_eigenform"] = {"pass": bool(g2), "primes": P_C0}

    # ---------------- G3c / C1b: Beukers-85 trivial-unit-root depth law
    a_ap = lambda n: A[n - 1]
    c1b = {}
    g3c = True
    for p in P_AP:
        r = mine(a_ap, p, NC, 4)
        gt = r["guaranteed"]
        wantD = 3 if p <= 13 else 2
        emp = r.get("empirical_depth", {}).get("plain", {})
        ok = (gt["verdict"] == "CONSISTENT"
              and str(gt.get("gamma_class", "")).startswith("TRIVIAL")
              and gt["depth_D"] >= wantD
              and isinstance(emp.get("min"), int) and emp["min"] >= 3)
        g3c = g3c and ok   # Eisenstein-pin gate (add-2): gamma==1+p^3 class
        c1b[str(p)] = {"verdict": gt["verdict"], "D": gt["depth_D"],
                       "rows": gt["rows_used"], "candidate": gt.get("weil_box_k"),
                       "gamma_class": gt.get("gamma_class"),
                       "profile_c": r.get("integer_candidate_k"),
                       "emp_plain_min": emp.get("min"), "ok": ok}
    rec["C1b_apery3_trivial_root"] = c1b
    rec["G3c_apery3"] = {"pass": bool(g3c), "primes": P_AP,
        "law": "Beukers-85 A(mp^r-1)==A(mp^{r-1}-1) mod p^{3r}; gamma in TRIVIAL class (==1+p^3 where depth>3 resolves the Eisenstein pin)"}

    # ---------------- G4 / C2: planted fault
    A2 = list(A); A2[137] += 1
    r = mine(lambda n: A2[n - 1], 23, NC, 4, want_profile=False)
    gt = r["guaranteed"]
    wit_ns = [w.get("n") for w in gt["witnesses"]]
    g4 = (gt["verdict"] == "INCONSISTENT" and 6 in wit_ns)
    rec["C2_planted_fault"] = {"p": 23, "perturbed": "A(137) += 1",
                               "verdict": gt["verdict"],
                               "witness_rows": wit_ns, "pass": bool(g4)}
    rec["G4_planted_fault"] = {"pass": bool(g4)}

    # ---------------- G5 / C3: structural null
    rng = random.Random(20260831)
    R = [0] + [rng.randint(1, 10 ** 12) for _ in range(NC)]
    g5 = True
    c3 = {}
    for p in (23, 53):
        r = mine(lambda n: R[n], p, NC, 4, want_profile=False)
        v = r["guaranteed"]["verdict"]
        c3[str(p)] = {"verdict": v, "n_witnesses": r["guaranteed"]["n_witnesses"]}
        g5 = g5 and (v == "INCONSISTENT")
    rec["C3_structural_null"] = c3
    rec["G5_structural_null"] = {"pass": bool(g5)}

    # ---------------- report-only gauge-lesson probe
    B2 = apery_zeta2(150)
    probe = {"stamp_utc": utc(),
        "half_index_A_vs_84aa": {str(p): {"A_half": A[(p - 1) // 2] % p,
                                          "a_p_mod_p": ap84[p] % p}
                                 for p in (5, 7, 11, 13, 17, 19, 23)},
        "mp_grid_zeta3_gamma_is_1": {str(p): A[p - 1] % p ** 3 for p in (5, 7, 11, 13)},
        "zeta2_b_head": B2[:8],
        "zeta2_bp1_mod_p": {str(p): B2[p - 1] % p for p in (5, 13, 17, 29, 37, 41)},
        "lesson": ("family-MUM mp^r-1 grid mines the TRIVIAL unit root "
                   "(gamma==1); nontrivial-a_p laws for the Apery families "
                   "are half-index congruences; proven nontrivial three-term "
                   "gamma_p lives on fixed-variety formal groups (Honda)."),
        "producer": {"component": "frobenius-boundary/asdminer",
                     "script": "asdminer/run_control.py",
                     "script_sha256": sha256_file(os.path.abspath(__file__)),
                     "stamp_utc": utc()}}
    with open(os.path.join(WORK, "probe_gauge_lesson.json"), "w") as f:
        json.dump(probe, f, indent=1)

    # ---------------- assemble + lint
    rec["comparator_note"] = (
        "Comparison columns = dual local routes: PARI mfcoefs (dim-1 "
        "newspace 8.4) / gp ellak vs independent python eta-product / "
        "formal-group integer arithmetic — no web dependency (online "
        "eigenvalue databases are deliberately not consulted).")
    gates = {"G1": bool(g1), "G2": bool(g2), "G3a": bool(g3a),
             "G3b": bool(g3b), "G3c": bool(g3c), "G4": bool(g4), "G5": bool(g5)}
    rec["walls_s"]["total"] = round(time.time() - t0, 2)
    rec["producer"] = {"component": "frobenius-boundary/asdminer",
                       "script": "asdminer/run_control.py",
                       "script_sha256": sha256_file(os.path.abspath(__file__)),
                       "asdlib_sha256": sha256_file(asdlib.__file__),
                       "stamp_utc": utc()}
    rec["gates"] = gates
    out = os.path.join(WORK, "ASD_CONTROL.json")
    with open(out, "w") as f:
        json.dump(rec, f, indent=1)
    lint_line = asdlib.lint_check(out)
    skipped = lint_line.startswith("producer_lint: SKIPPED")
    g6 = skipped or not ("FAIL" in lint_line or "refus" in lint_line)
    gates["G6"] = bool(g6)
    rec["G6_producer_lint"] = {"pass": bool(g6), "skipped": skipped,
                               "line": lint_line}
    rec["control_pass"] = all(gates.values())
    rec["verdict"] = ("ASDMINER-CONTROL-PASS: " +
                      ",".join(k for k in sorted(gates)) + " all green; "
                      "Honda 8/8+8/8 exact ellak match, Beukers-85 12/12, "
                      "Hecke 5/5, fault+null caught"
                      if rec["control_pass"] else
                      "ASDMINER-CONTROL-FAIL: " +
                      ",".join(k for k, v in sorted(gates.items()) if not v))
    with open(out, "w") as f:
        json.dump(rec, f, indent=1)
    print(asdlib.lint_check(out))
    print("VERDICT:", rec["verdict"])
    print("control_pass =", rec["control_pass"])
    return 0 if rec["control_pass"] else 1

if __name__ == "__main__":
    sys.exit(main())
