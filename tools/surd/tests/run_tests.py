#!/usr/bin/env python3
"""Surd test driver.  One command; exit code 0 = every selected test that ran passed (skipped tests are named).
   python3 tests/run_tests.py [--only h,p,g,n,a,b,c,d2,e,s12] [--data DIR] [--n4slice PATH] [--quick] [--dps 60]

Hermetic tests (default; read only the package data, write nothing, ~5 s total):
 (h)   hyperlog evaluator self-check (gate/scripts/hlog_eval.selfgate): G(w;1) and ZIP_reg words with real/complex algebraic letters against
       closed forms (log, Li2, Li3) and nested quadrature at 50 digits, every difference < 1e-45;
 (p)   the on-path 'letter - i0' prescription of slice/scripts/hpath: G(1/2 - i0; 1) = -i pi and G(0, 1/2 - i0; 1) = -Li2(2 + i0) against an
       independent contour quadrature (path deformed into Im > 0), 45 digits;
 (g)   included cell data: the exact Cheng-Wu cells of channel q_qbpqpgq in the three charts x_1 = 1, x_2 = 1, x_4 = 1 (gate/cells/*.pkl) give the
       SAME degree -4 homogeneous integrand at random exact rational points and slice kinematics (exact Fraction identity), and the first-pair
       group (34, chart 2) has the frozen census (6 cells, 1857 numerator monomials);
 (n)   data/N4_SLICE.json: sha256 prefix and 34 canonical irreducible slice letters.
Engine test (opt-in, --only s12, ~50 s, writes under $SURD_WORK or a temporary directory):
 (s12) stages 1-2 of the group q_qbpqpgq_s1234_p34_g2 and the stage 1-2 check at t = 1/2, x = 7/3: symbolic one-fold object vs the independent
       residue-route fiber function, >= 30 digits (tests/smoke_test.sh runs the whole chain including stage 3).
Data tests (need the assembled files out/*.json.gz and SLICE_S5_LETTERS.json, ~0.5 GB, NOT part of the package; give --data DIR or $SURD_DATA, else
they SKIP by name; all five ~56 min single thread, 2.1 GiB; b + d2 ~8 min):
 (a)   n4 closed form vs frozen CMSYZ (arXiv:2401.06463) values at t = 1/2, 2/3, 4/5, 3/7: |24 F/(5 CMSYZ_Re) - 1| <= 2e-26, and F(1/2) to >= 40 digits
       of the frozen file value;
 (b)   N=4 prime support == the 34 CMSYZ slice letters (SLICE_S5_LETTERS.json vs data/N4_SLICE.json);
 (c)   symtest (fixed labeling) on n4 and quark at (e/7, pi/9) and (1/sqrt7, pi/9): frozen verdict pattern (every Table-1 class CANCELS <= 1e-58,
       the two genuine classes SURVIVE), LLL verified == found; --quick: n4 only;
 (d2)  quark jet at t = 1/2 vs the frozen S0 reference q_jet value (slice_oracle.py) (>= 36 digits) AND identical value (<= 1e-50) with every root list REVERSED;
 (e)   evaluate.py at representation-singular t: quark@3/5 -> 'symmetric limit', n4@3/5 -> 'direct', q_qbqgq@5/8 -> 'symmetric limit', each
       >= 29 digits of the frozen values.
Single thread (OMP_NUM_THREADS=1 is set here)."""
import os, sys, json, gzip, time, argparse, traceback, hashlib, tempfile
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); SCRIPTS = os.path.join(ROOT, "slice", "scripts")
if not os.path.exists(os.path.join(SCRIPTS, "slice_prov.py")): sys.exit("run_tests.py: cannot find slice/scripts next to %s" % HERE)
os.environ.setdefault("OMP_NUM_THREADS", "1")
sys.path.insert(0, SCRIPTS)
from fractions import Fraction as Fr
import mpmath as mp

FROZEN = {
 "cmsyz_Re": {"1/2": "1468.1493450846724778184508646", "2/3": "1016.80770027829552600981296509", "4/5": "734.300496586676388389623415671", "3/7": "1734.992090049981175025454775084718202052"},
 "n4_file_1_2": "305.86444689264009954551059734352338942858597093132",
 "S0_qjet_1_2": "349.346865156319860592390498268165579897",           # S0 oracle (slice_oracle.py) q_jet at t = 1/2, 36 certified digits
 "n4_3_5": "245.87527586150966006840229747961171011",                 # evaluate.py direct
 "q_qbqgq_5_8": "8.184320525040340232233572530595188909",              # evaluate.py symmetric limit, delta = 1e-20
 "quark_3_5": "282.5818179378116526279331062518124624777",            # evaluate.py symmetric limit
 "symtest_w3": {"quark": {"('25*t**2 + 7', 't + 1', 't - 1', 't**2 + 1')": "SURVIVES", "('-1', '7*t**2 + 25', 't + 1', 't - 1', 't**2 + 1')": "SURVIVES", "('-1', '8*t + 5', '8*t - 5')": "CANCELS", "('-1', '35*t**2 + 6*t - 45', '5*t**2 - 6*t + 5')": "CANCELS", "('-1', '35*t**2 - 6*t - 45', '5*t**2 + 6*t + 5')": "CANCELS", "('45*t**2 + 6*t - 35', '5*t**2 + 6*t + 5')": "CANCELS", "('5*t + 8', '5*t - 8')": "CANCELS", "('45*t**2 - 6*t - 35', '5*t**2 - 6*t + 5')": "CANCELS"},
                "n4": {"('-1', '7*t**2 + 25', 't + 1', 't - 1', 't**2 + 1')": "SURVIVES", "('25*t**2 + 7', 't + 1', 't - 1', 't**2 + 1')": "SURVIVES"}},
 "n4slice_sha16": "cc3bd6b61fec4a1d",
 "group_q_qbpqpgq_p34_g2": {"n_cells": 6, "N_monomials": 1857, "n_cells_chart": 30, "n_monomials_chart": 45746},
 "s12_gate_t_x": ("1/2", "7/3"),
}
FILES = {"quark": "E4C_LO_QCD_dipole_slice_quark.json.gz", "gluon": "E4C_LO_QCD_dipole_slice_gluon.json.gz", "n4": "E4C_LO_N4_dipole_slice_n4.json.gz", "q_qbqgq": "E4C_LO_QCD_dipole_slice_q_qbqgq.json.gz"}
HERMETIC = ("h", "p", "g", "n"); DATA_TESTS = ("a", "b", "c", "d2", "e"); OPTIN = ("s12",)
def digits(a, b):
    a, b = mp.mpf(a), mp.mpf(b); return 99.0 if a == b else float(-mp.log10(abs(a - b) / abs(b)))
class Suite:
    def __init__(self, data, n4slice, dps, quick):
        self.data = data; self.n4slice = n4slice; self.dps = dps; self.quick = quick; self.docs = {}; self.results = []
    def path(self, lab): return os.path.join(self.data, "out", FILES[lab])
    def doc(self, lab):
        if lab not in self.docs: t0 = time.time(); self.docs[lab] = json.load(gzip.open(self.path(lab), "rt")); print("   loaded %s (%d terms, %.0fs)" % (lab, len(self.docs[lab]["terms"]), time.time() - t0), flush=True)
        return self.docs[lab]
    def rec(self, name, ok, **kw): self.results.append({"test": name, "PASS": bool(ok), **kw}); print("[%s] %s %s" % ("PASS" if ok else "FAIL", name, json.dumps(kw, default=str)[:600]), flush=True)
    def skip(self, name, why): self.results.append({"test": name, "SKIP": why}); print("[SKIP] %s: %s" % (name, why), flush=True)
    # ---------------------------------------------------------------- hermetic
    def test_h(self):
        import slice_prov, hlog_eval          # slice_prov puts gate/scripts on sys.path
        t0 = time.time(); r = hlog_eval.selfgate(50); mx = max(r.values()); self.rec("h_hlog_selfgate", mx < 1e-45, max_abs_diff=mx, n_checks=len(r), wall_s=round(time.time() - t0, 1))
    def test_p(self):
        import hpath
        t0 = time.time(); mp.mp.dps = 50; H = hpath.HPath(50); a = Fr(1, 2)
        g1 = H.G1([(a, "below")]); d1 = float(abs(g1 - mp.mpc(0, -1) * mp.pi))
        g2 = H.G1([(Fr(0), None), (a, "below")]); am = mp.mpf(1) / 2
        ref = mp.quad(lambda s: mp.log(1 - s / am) / s, [0, mp.mpc("0.5", "0.25"), 1]); d2 = float(abs(g2 - ref))
        li = -mp.polylog(2, 2); d3 = float(abs(g2 - mp.mpc(li.real, -abs(li.imag))))          # -Li2(2 + i0): Im Li2(x + i0) = +pi log x for x > 1
        ok = d1 < 1e-45 and d2 < 1e-45 and d3 < 1e-45
        self.rec("p_hpath_minus_i0_prescription", ok, **{"|G(1/2-i0;1) + i pi|": d1, "|G(0,1/2-i0;1) - contour quad|": d2, "|G(0,1/2-i0;1) + Li2(2+i0)|": d3, "wall_s": round(time.time() - t0, 1)})
    def test_g(self):
        import random, slice_common as sc, oracle_cells as OC
        t0 = time.time(); fz = FROZEN["group_q_qbpqpgq_p34_g2"]
        cells = {g: OC.build_cells("q_qbpqpgq", g, verbose=False) for g in (1, 2, 4)}
        cen_ok = all(c["n_cells"] == fz["n_cells_chart"] and c["n_monomials"] == fz["n_monomials_chart"] for c in cells.values())
        def phi(cellres, y, d):
            g = cellres["gauge"]; xs = {j: Fr(y[j - 1]) / Fr(y[g - 1]) for j in (1, 2, 3, 4)}
            return OC.eval_cells_exact(cellres, xs, d) / Fr(y[g - 1]) ** 4
        rng = random.Random(5); eq = []
        for _ in range(3):
            t = Fr(rng.randint(1, 19), rng.randint(2, 23)); sig = tuple(rng.sample([1, 2, 3, 4], 4)); ds = sc.sigma_d(sc.d_of_t(t), sig)
            y = [Fr(rng.randint(1, 9), rng.randint(1, 9)) for _ in range(4)]; v = {g: phi(cells[g], y, ds) for g in (1, 2, 4)}
            eq.append({"t": str(t), "sigma": "".join(map(str, sig)), "charts_1_2_4_identical": v[1] == v[2] == v[4], "value_float": float(v[1])})
        G = sc.build_group("q_qbpqpgq", 2, "34", (1, 2, 3, 4))
        grp_ok = G["n_cells"] == fz["n_cells"] and len(G["N"]) == fz["N_monomials"]
        ok = cen_ok and grp_ok and all(e["charts_1_2_4_identical"] for e in eq)
        self.rec("g_cells_cross_chart_exact", ok, chart_census_ok=cen_ok, group_census_ok=grp_ok, group={"n_cells": G["n_cells"], "N_monomials": len(G["N"])}, points=eq, wall_s=round(time.time() - t0, 1))
    def test_n(self):
        import sympy as sp; T = sp.Symbol("t")
        t0 = time.time(); cz = json.load(open(self.n4slice)); sh = hashlib.sha256(open(self.n4slice, "rb").read()).hexdigest()[:16]
        def canon(expr):
            P = sp.Poly(sp.sympify(str(expr).replace("^", "**")), T)
            if P.degree() <= 0: return None
            P = P.primitive()[1]; P = -P if P.LC() < 0 else P; return str(P.as_expr())
        theirs = set(filter(None, (canon(k) for k in cz["letters"])))
        ok = sh == FROZEN["n4slice_sha16"] and len(theirs) == 34
        self.rec("n_N4_SLICE_data", ok, sha16=sh, sha_expected=FROZEN["n4slice_sha16"], n_letters_in_file=len(cz["letters"]), n_canonical_irreducible=len(theirs), wall_s=round(time.time() - t0, 1))
    # ---------------------------------------------------------------- opt-in engine test
    def test_s12(self):
        import subprocess, slice_common as sc, slice_eval, slice_gate12
        work = os.environ["SURD_WORK"]; tag = "q_qbpqpgq_s1234_p34_g2"; t0 = time.time()
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "slice_run.py"), "--channel", "q_qbpqpgq", "--sigma", "1234", "--pair", "34"], cwd=SCRIPTS, capture_output=True, text=True)
        if p.returncode != 0: self.rec("s12_stages12_and_S1_gate", False, rc=p.returncode, stderr=p.stderr[-400:]); return
        w12 = time.time() - t0
        C, letters, ts, extra = sc.ts_load(os.path.join(work, "ckpt", tag, "stage2.pkl")); G = sc.build_group("q_qbpqpgq", extra["gauge"], "34", (1, 2, 3, 4))
        tv, xv = (Fr(s) for s in FROZEN["s12_gate_t_x"]); S = slice_eval.Specialized(C, letters, ts, tv, "x%d" % extra["order"][2]); ef = sc.ext_factor_value(G, tv)
        val, mx = S.value(xv, 60, return_terms=True); mine = val / (mp.mpf(ef.numerator) / ef.denominator)
        d = sc.sigma_d(sc.d_of_t(tv), (1, 2, 3, 4)); best = -1.0
        for roles, fv, st in slice_gate12.fibre_real("q_qbpqpgq", extra["gauge"], "34", d, extra["order"][2], xv, Fr(1, 32), 4.7, 34):
            ref = mp.mpf(fv.mid().str(70, radius=False)); rel = abs(mine - ref) / abs(ref); best = max(best, 99.0 if rel == 0 else float(-mp.log10(rel)))
        self.rec("s12_stages12_and_S1_gate", best >= 30, n_stage2_keys=len(ts), symbolic_over_ext=mp.nstr(mine.real, 35), agree_digits_vs_residue_route=round(best, 1), stages12_wall_s=round(w12), wall_s=round(time.time() - t0), work=work)
    # ---------------------------------------------------------------- data tests
    def test_a(self):
        import assemble; mp.mp.dps = self.dps; ok = True; rows = {}
        for t in ("1/2", "2/3", "4/5", "3/7"):
            t0 = time.time(); F = assemble.evaluate_file(self.path("n4"), Fr(t), self.dps, doc=self.doc("n4")); C = mp.mpf(FROZEN["cmsyz_Re"][t])
            r = abs(24 * F.real / (5 * C) - 1); good = r <= mp.mpf("2e-26") and abs(F.imag) < mp.mpf(10) ** -40; ok &= bool(good)
            rows[t] = {"F": mp.nstr(F.real, 42), "|24F/5C-1|": mp.nstr(r, 3), "abs_Im": mp.nstr(abs(F.imag), 3), "ok": bool(good), "wall_s": round(time.time() - t0)}
            if t == "1/2": d = digits(F.real, FROZEN["n4_file_1_2"]); rows[t]["digits_vs_frozen_file"] = round(d, 1); ok &= d >= 40
        self.rec("a_n4_vs_frozen_CMSYZ", ok, rows=rows)
    def test_b(self):
        import sympy as sp; T = sp.Symbol("t")
        def canon(expr):
            P = sp.Poly(sp.sympify(str(expr).replace("^", "**")), T)
            if P.degree() <= 0: return None
            P = P.primitive()[1]; P = -P if P.LC() < 0 else P; return str(P.as_expr())
        LT = json.load(open(os.path.join(self.data, "SLICE_S5_LETTERS.json")))["functions"]["n4"]; cz = json.load(open(self.n4slice))
        sh = hashlib.sha256(open(self.n4slice, "rb").read()).hexdigest()[:16]
        theirs = set(filter(None, (canon(k) for k in cz["letters"])))
        ours = set()
        for f in LT["coefficient_denominator_t_letters"]["irreducible_factors(count of denominators containing it)"]: ours.add(canon(f))
        for f in LT["rational_point_t_letters"]["factors_of_norms_of_points_and_1+point"]: ours.add(canon(f))
        for f in LT["rational_point_t_letters"]["factors_of_norms_of_differences"]: ours.add(canon(f))
        for v in LT["algebraic_letters"].values():
            for f in (v.get("sqrt_field_kernel(odd factors)") or []) + (v.get("disc_square_class(odd factors)") or []): ours.add(canon(f))
        ours.discard(None)
        ok = (len(ours & theirs) == 34) and not (ours - theirs) and not (theirs - ours)
        self.rec("b_n4_support_34_of_34", ok, both=len(ours & theirs), only_ours=sorted(ours - theirs), only_cmsyz=sorted(theirs - ours), N4_SLICE_sha16=sh, sha_expected=FROZEN["n4slice_sha16"], sha_matches=(sh == FROZEN["n4slice_sha16"]))
    def test_c(self):
        import symtest; ok = True; rows = {}
        for lab in (("n4",) if self.quick else ("n4", "quark")):
            for pair in (("e/7", "pi/9"), ("1/sqrt7", "pi/9")):
                t0 = time.time(); r = symtest.run(lab, self.path(lab), dps=self.dps, tstar_name=pair[0], tstar2_name=pair[1], continuation=True); mp.mp.dps = self.dps
                verd = {cls: rec["weights"]["3"]["verdict"].split(" (")[0] for cls, rec in r["classes"].items() if "weights" in rec}
                lll_ok = all(rec["LLL"]["n_relations_found_at_tstar"] == rec["LLL"]["n_verified_at_tstar2"] for rec in r["classes"].values() if rec.get("LLL"))
                maxcancel = max([float(v["ratio"]) for cls, rec in r["classes"].items() if verd.get(cls) == "CANCELS" for v in rec["weights"]["3"]["functionals"].values()] or [0.0])
                low_ok = all(rec["weights"][w]["verdict"].startswith("no terms") for rec in r["classes"].values() if "weights" in rec for w in ("1", "2"))
                good = (verd == FROZEN["symtest_w3"][lab]) and lll_ok and maxcancel <= 1e-58 and low_ok; ok &= good
                rows["%s@%s,%s" % (lab, pair[0], pair[1])] = {"verdicts_match_frozen": verd == FROZEN["symtest_w3"][lab], "LLL_verified==found": lll_ok, "max_cancel_ratio": maxcancel, "weights_1_2_empty": low_ok, "wall_s": round(time.time() - t0)}
        self.rec("c_symtest_two_pairs", ok, rows=rows)
    def test_d2(self):
        import assemble; mp.mp.dps = self.dps
        t0 = time.time(); F = assemble.evaluate_file(self.path("quark"), Fr(1, 2), self.dps, doc=self.doc("quark")); w1 = time.time() - t0
        d_s0 = digits(F.real, FROZEN["S0_qjet_1_2"])
        orig = mp.polyroots
        def rev(*a, **k): return list(reversed(orig(*a, **k)))
        mp.polyroots = rev
        try: t0 = time.time(); Fr_ = assemble.evaluate_file(self.path("quark"), Fr(1, 2), self.dps, doc=self.doc("quark")); w2 = time.time() - t0
        finally: mp.polyroots = orig
        d_perm = digits(Fr_.real, F.real)
        ok = d_s0 >= 36 and d_perm >= 50
        self.rec("d2_quark_half_vs_S0_and_slot_permutation", ok, F=mp.nstr(F.real, 45), digits_vs_S0=round(d_s0, 1), F_reversed_roots=mp.nstr(Fr_.real, 45), digits_reversed_vs_direct=round(d_perm, 1), abs_Im=mp.nstr(abs(F.imag), 3), wall_s=[round(w1), round(w2)])
    def test_e(self):
        import evaluate as EV; ok = True; rows = {}
        for lab, t, mode, frozen in (("n4", "3/5", "direct", FROZEN["n4_3_5"]), ("q_qbqgq", "5/8", "symmetric limit", FROZEN["q_qbqgq_5_8"]), ("quark", "3/5", "symmetric limit", FROZEN["quark_3_5"])):
            t0 = time.time(); r = EV.evaluate(self.path(lab), t, self.dps, 30, doc=self.doc(lab)); mp.mp.dps = self.dps
            val = r.get("value") or ""; good = r["mode"].startswith(mode) and r["digits_est"] >= 30 and bool(val) and digits(val, frozen) >= 29
            ok &= bool(good); rows["%s@%s" % (lab, t)] = {"mode": r["mode"], "value": val[:42], "frozen": frozen[:42], "digits_est": r["digits_est"], "ok": bool(good), "wall_s": round(time.time() - t0)}
        self.rec("e_singular_t_trio", ok, rows=rows)
def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--data", default=os.environ.get("SURD_DATA", ""), help="directory holding out/*.json.gz and SLICE_S5_LETTERS.json for the data tests")
    ap.add_argument("--n4slice", default=os.environ.get("SURD_N4SLICE", "") or os.path.join(ROOT, "data", "N4_SLICE.json"))
    ap.add_argument("--only", default=",".join(HERMETIC + DATA_TESTS), help="comma list of tests (default: hermetic + data tests; s12 is opt-in)"); ap.add_argument("--quick", action="store_true"); ap.add_argument("--dps", type=int, default=60); a = ap.parse_args()
    names = [x.strip() for x in a.only.split(",") if x.strip()]
    if "s12" in names and not os.environ.get("SURD_WORK"): os.environ["SURD_WORK"] = tempfile.mkdtemp(prefix="surd_tests.")   # the only test that writes; must be set before the slice modules are imported
    data = a.data if (a.data and os.path.exists(os.path.join(a.data, "out", FILES["n4"]))) else ""
    print("Surd tests: scripts %s | data %s | N4_SLICE %s | dps %d | only %s%s" % (SCRIPTS, data or "(none: data tests skip)", a.n4slice, a.dps, a.only, " (quick)" if a.quick else ""), flush=True)
    S = Suite(data, a.n4slice, a.dps, a.quick); T0 = time.time()
    for name in names:
        fn = getattr(S, "test_" + name, None)
        if fn is None: print("unknown test", name); continue
        if name in DATA_TESTS and not data: S.skip(name, "needs --data DIR or SURD_DATA (assembled files out/*.json.gz; not part of the package)"); continue
        try: fn()
        except Exception as ex: S.rec(name, False, error=repr(ex)[:300], traceback=traceback.format_exc()[-600:])
    ran = [r for r in S.results if "PASS" in r]; skipped = [r["test"] for r in S.results if "SKIP" in r]
    allok = bool(ran) and all(r["PASS"] for r in ran)
    print(json.dumps({"results": S.results, "ALL_PASS": allok, "n_ran": len(ran), "skipped": skipped, "wall_s": round(time.time() - T0)}, default=str)[:4000])
    print("RESULT:", "PASS" if allok else "FAIL", "(%d ran, %d skipped: %s)" % (len(ran), len(skipped), ",".join(skipped) or "-")); sys.exit(0 if allok else 1)
if __name__ == "__main__":
    main()
