#!/usr/bin/env python3
"""run_l6_gauge.py — stage A: the L6 gauge receipt.
Derives, per point x0 in {+1, 0, -1}: local exponent data (exact indicial
polynomial + all rational roots), the log-free-exponent scan with exact
obstruction certificates, the chosen candidate series (exact, N=600),
denominator laws (small-prime slopes + per-menu-prime integrality profile),
and the exact-vs-modular co-check at the pilot prime 53.
The operator is NOT distributed with the repo: set FROB_L6_OPERATOR to a
JSON operator file ({"coeffs": [[...], ...]}); refuses loudly otherwise.
Emits work/L6_GAUGE.json."""
import json, os, subprocess, sys, time, datetime
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import asdlib
from asdlib import (load_operator, shift_operator, op_local_data,
                    frobenius_exact, frobenius_modular, vp_int, sha256_file)

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.environ.get("FROB_ASDMINER_WORK") or os.path.join(HERE, "work")
L6PATH = os.environ.get("FROB_L6_OPERATOR")
NEX = 600
PILOT_P = 53

def utc():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

def kron60(p):
    """kronecker(60,p) for odd prime p via Euler criterion."""
    return pow(60 % p, (p - 1) // 2, p)

def menu_primes():
    out = []
    for p in range(53, 200):
        if all(p % q for q in range(2, p)) and kron60(p) == 1:
            out.append(p)
    return out

def rational_roots(Ipoly):
    """All rational roots of the integer polynomial (exact), with mult."""
    P = list(Ipoly)
    roots = []
    def val(P, s):  # s Fraction
        v = Fraction(0)
        for cf in reversed(P):
            v = v * s + cf
        return v
    # integer + small-denominator scan (denominators dividing lead coeff)
    lead = P[-1]
    dens = sorted({d for d in range(1, 13) if lead % d == 0})
    cands = set()
    for b in dens:
        for a in range(-40 * b, 40 * b + 1):
            cands.add(Fraction(a, b))
    for s in sorted(cands):
        while val(P, s) == 0:
            roots.append(str(s))
            # deflate: synthetic division by (x - s) exactly
            Q = [Fraction(0)] * (len(P) - 1)
            acc = Fraction(0)
            for i in range(len(P) - 1, 0, -1):
                acc = acc * s + P[i]
                Q[i - 1] = acc
            P = Q
            if len(P) <= 1:
                break
        if len(P) <= 1:
            break
    return roots, len(P) - 1  # residual degree (unresolved part)

def main():
    t0 = time.time()
    if not L6PATH or not os.path.exists(L6PATH):
        raise SystemExit("run_l6_gauge: the L6-class operator is not distributed "
                         "with the repo; set FROB_L6_OPERATOR to a JSON operator file")
    os.makedirs(WORK, exist_ok=True)
    d, C = load_operator(L6PATH)
    rec = {"stamp_utc_start": utc(),
           "operator": {"path": L6PATH, "sha256": sha256_file(L6PATH),
                        "order": d["order"], "deg": d["deg"],
                        "gauge_note_carried": d["note_indicial"]},
           "menu_primes_split_Qsqrt15": menu_primes(),
           "N_exact": NEX, "pilot_prime": PILOT_P, "points": {}, "walls_s": {}}
    series_store = {}
    for x0 in (1, 0, -1):
        t = time.time()
        Cs = shift_operator(C, x0)
        loc = op_local_data(Cs)
        rr, resid = rational_roots(loc["Ipoly"])
        pt = {"local_coordinate": f"u = x - ({x0})",
              "d0": loc["d0"], "dmax": loc["dmax"],
              "indicial_rational_roots": rr,
              "indicial_residual_degree": resid,
              "rho_scan": []}
        chosen = None
        for rho in sorted(set(int(Fraction(r)) for r in rr
                              if Fraction(r).denominator == 1
                              and Fraction(r) >= 0)):
            try:
                b, reslist, _ = frobenius_exact(Cs, rho, NEX)
                pt["rho_scan"].append({"rho": rho, "status": "LOG-FREE",
                                       "resonances_zeroed": reslist})
                chosen = (rho, b)
                break
            except RuntimeError as e:
                pt["rho_scan"].append({"rho": rho, "status": str(e)})
        if chosen is None:
            pt["candidate"] = "NONE-LOG-FREE (wall)"
            rec["points"][str(x0)] = pt
            continue
        rho, b = chosen
        series_store[x0] = (Cs, rho, b)
        pt["candidate"] = {"rho": rho,
            "normalization": "b_0 = 1; resonant coefficients zeroed "
                             "(Frobenius gauge); mining index "
                             "a(n) = b_(n-1) (formal-group law)"}
        # denominator laws
        dens = [x.denominator for x in b]
        def vq(x, q):
            v = 0
            while x % q == 0:
                x //= q; v += 1
            return v
        pt["denominator_law"] = {
            "small_prime_slopes_v_at_300_600": {
                str(q): [vq(dens[300], q), vq(dens[600], q)] for q in (2, 3, 5)},
            "menu_prime_profile": {}}
        for p in rec["menu_primes_split_Qsqrt15"]:
            first = next((m for m in range(len(dens)) if dens[m] % p == 0), None)
            vmax = max(vq(dens[m], p) for m in range(len(dens)))
            pt["denominator_law"]["menu_prime_profile"][str(p)] = {
                "first_m_with_p_denominator": first, "v_max_to_600": vmax}
        pt["p_integral_at_all_menu_primes_to_600"] = all(
            v["first_m_with_p_denominator"] is None
            for v in pt["denominator_law"]["menu_prime_profile"].values())
        rec["points"][str(x0)] = pt
        rec["walls_s"][f"x0_{x0}"] = round(time.time() - t, 2)
    # COPAIR at pilot prime for each minable point (exact vs modular route)
    rec["copair_p53"] = {}
    for x0, (Cs, rho, b) in series_store.items():
        prof = rec["points"][str(x0)]["denominator_law"]["menu_prime_profile"][str(PILOT_P)]
        if prof["first_m_with_p_denominator"] is not None:
            rec["copair_p53"][str(x0)] = {
                "mode": "wall-step-equality",
                "exact_first_pden_m": prof["first_m_with_p_denominator"]}
            try:
                frobenius_modular(Cs, rho, NEX, PILOT_P, 4,
                                  [rz["resonances_zeroed"] for rz in
                                   rec["points"][str(x0)]["rho_scan"]
                                   if rz["status"] == "LOG-FREE"][0])
                rec["copair_p53"][str(x0)]["modular_wall_step"] = None
                rec["copair_p53"][str(x0)]["agree"] = False
            except RuntimeError as e:
                step = int(str(e).split("n=")[1]) if "n=" in str(e) else None
                rec["copair_p53"][str(x0)]["modular_wall_step"] = step
                rec["copair_p53"][str(x0)]["agree"] = (
                    step == prof["first_m_with_p_denominator"])
        else:
            t = time.time()
            reszero = [rz["resonances_zeroed"] for rz in
                       rec["points"][str(x0)]["rho_scan"]
                       if rz["status"] == "LOG-FREE"][0]
            bm, F, W, V = frobenius_modular(Cs, rho, NEX, PILOT_P, 4, reszero)
            pf = PILOT_P ** F
            ok = all(bm[m] % pf == (b[m].numerator *
                     pow(b[m].denominator, -1, pf)) % pf
                     for m in range(NEX + 1))
            rec["copair_p53"][str(x0)] = {
                "mode": "integral-equality", "floor_F": F, "working_W": W,
                "V_total": V, "n_overlap": NEX + 1, "agree": bool(ok),
                "wall_s": round(time.time() - t, 2)}
    agree_all = all(v.get("agree") for v in rec["copair_p53"].values())
    rec["copair_verdict"] = ("COPAIR-PASS: exact and modular routes agree at "
                             "p=53 on every candidate point"
                             if agree_all and rec["copair_p53"] else
                             "COPAIR-FAIL")
    rec["walls_s"]["total"] = round(time.time() - t0, 2)
    rec["producer"] = {"component": "frobenius-boundary/asdminer",
                       "script": "asdminer/run_l6_gauge.py",
                       "script_sha256": sha256_file(os.path.abspath(__file__)),
                       "asdlib_sha256": sha256_file(asdlib.__file__),
                       "stamp_utc": utc()}
    out = os.path.join(WORK, "L6_GAUGE.json")
    with open(out, "w") as f:
        json.dump(rec, f, indent=1)
    print(asdlib.lint_check(out))
    print("COPAIR:", rec["copair_verdict"])
    for x0 in (1, 0, -1):
        pt = rec["points"][str(x0)]
        c = pt.get("candidate")
        print(f"x0={x0}: cand={'rho='+str(c['rho']) if isinstance(c,dict) else c}"
              f" p-integral(menu,<=600)={pt.get('p_integral_at_all_menu_primes_to_600')}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
