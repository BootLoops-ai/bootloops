"""The letter table on the slice.  From the assembled closed forms (out/*.json.gz): (i) the algebraic X-letters (degree >= 2 polynomials over Q(i)(t)
whose roots are hyperlog letters) per function, their discriminants factored over Q[t] and matched against the twelve rational tangency letters P(t)
predicted for the dipole slice by the Landau analysis (TABLE1 below; called the 'Table-1 letters' throughout this package);
(ii) cubic letters: discriminant square class and field isomorphism vs the CMSYZ cubic P restricted to the slice; (iii) the rational t-letters:
irreducible Q[t]-factors of coefficient denominators and of norms of rational points / their pairwise differences inside words; (iv) per quadratic
letter: the net contribution S_L^(w)(t) of all terms containing its roots, per transcendental weight w (numeric at gate t's): S_L == 0 <=> the
sqrt(disc L) letters cancel at function level in this representation.  Writes SLICE_S5_LETTERS.json."""
import os, sys, json, gzip, time, argparse, itertools, collections
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import slice_prov, assemble, k3field
from fractions import Fraction as Fr
import sympy as sp, mpmath as mp
from flint import fmpq_poly, fmpq
from k3field import GT
WORK = slice_prov.WORK; DATA = slice_prov.DATA; OUT = os.path.join(DATA, "out")
T = sp.Symbol("t")
TABLE1 = {"t+5/8": 8*T+5, "t-5/8": 8*T-5, "t+8/5": 5*T+8, "t-8/5": 5*T-8, "t^2+2t/15-7/9": 45*T**2+6*T-35, "t^2-2t/15-7/9": 45*T**2-6*T-35,
          "t^2+6t/35-9/7": 35*T**2+6*T-45, "t^2-6t/35-9/7": 35*T**2-6*T-45, "t^2-1/2": 2*T**2-1, "t^2-2": T**2-2,
          "quartic+": 225*T**4+60*T**3-826*T**2+60*T+225, "quartic-": 225*T**4-60*T**3-826*T**2-60*T+225}
def fp2sp(p): return sum(sp.Rational(int(c.p), int(c.q)) * T ** k for k, c in enumerate(p.coeffs()))
def gt2sp(g): return (fp2sp(g.A) + sp.I * fp2sp(g.B)) / fp2sp(g.D)
def _flint_factor_sp(pol):
    """irreducible primitive factors of a sympy Poly in t via FLINT (fast for the degree-400 denominators): list of (sympy expr, multiplicity)"""
    P = sp.Poly(pol, T); cs = [sp.Rational(c) for c in reversed(P.all_coeffs())]
    fp = fmpq_poly([fmpq(int(c.p), int(c.q)) for c in cs]); out = []
    cont, facs = fp.factor()
    for f, e in facs:
        if f.degree() < 1: continue
        co = [sp.Rational(int(f[i].p), int(f[i].q)) for i in range(f.degree() + 1)]
        fe = sp.Poly(sum(c * T ** i for i, c in enumerate(co)), T).primitive()[1]
        if fe.LC() < 0: fe = -fe
        out.append((fe.as_expr(), int(e)))
    return out
def qfactors(expr):
    """irreducible Q[t] factors (primitive) with multiplicity of a rational function in t (numerator +, denominator -); FLINT factoring"""
    num, den = sp.fraction(sp.cancel(sp.together(expr))); out = collections.Counter()
    for pol, sgn in ((num, 1), (den, -1)):
        if sp.Poly(pol, T).degree() <= 0: continue
        for f, e in _flint_factor_sp(pol): out[str(f)] += sgn * e
    return out
def table1_hits(factors):
    hits = {}
    for nm, pol in TABLE1.items():
        key = str(sp.Poly(sp.Poly(pol, T).primitive()[1], T).as_expr())
        for f, e in factors.items():
            if sp.simplify(sp.Poly(sp.sympify(f), T).primitive()[1].as_expr() - sp.Poly(pol, T).primitive()[1].as_expr()) == 0 or sp.simplify(sp.Poly(sp.sympify(f), T).primitive()[1].as_expr() + sp.Poly(pol, T).primitive()[1].as_expr()) == 0: hits[nm] = e
    return hits
def analyze_file(path, label, tvals=("1/2", "4/5"), dps=60, do_S=True):
    doc = json.load(gzip.open(path, "rt")); res = {"file": path, "label": label, "complete": doc.get("complete"), "census": doc["census"]}
    letters = doc["letters"]; alg = {}
    # (i) algebraic X-letters
    for L, info in letters.items():
        if info["deg_x"] < 2: continue
        cs = [GT.deser(c) for c in info["x_coeffs(GT ser, increasing powers of X)"]]; real = all(c.is_real() for c in cs)
        X = sp.Symbol("X"); pol = sum(gt2sp(c) * X ** k for k, c in enumerate(cs))
        rec = {"poly": info["poly"], "deg_x": info["deg_x"], "coefficients_real(Q(t))": real}
        if info["deg_x"] == 2:
            disc = sp.cancel(gt2sp(cs[1]) ** 2 - 4 * gt2sp(cs[2]) * gt2sp(cs[0]))
            if real:
                fac = qfactors(disc); rec["disc_Q[t]_factors(exponent)"] = dict(fac); odd = {f: e for f, e in fac.items() if e % 2}; rec["sqrt_field_kernel(odd factors)"] = sorted(odd)
                rec["table1_letters_in_kernel"] = table1_hits(collections.Counter(odd))
            else:
                nrm = sp.cancel(sp.expand(disc * disc.subs(sp.I, -sp.I))); fac = qfactors(nrm); rec["disc_in_Q(i)(t)"] = str(disc)[:300]; rec["norm(disc)_Q[t]_factors"] = dict(fac); rec["table1_letters_in_norm(odd)"] = table1_hits(collections.Counter({f: e for f, e in fac.items() if e % 2}))
        elif info["deg_x"] == 3:
            Pl = sp.Poly(pol, X); disc = sp.cancel(sp.discriminant(Pl, X)); rec["coefficients_real(Q(t))"] = real
            if real:
                fac = qfactors(disc); rec["disc_Q[t]_factors"] = dict(fac); rec["disc_square_class(odd factors)"] = sorted(f for f, e in fac.items() if e % 2); rec["table1_letters_in_disc(odd)"] = table1_hits(collections.Counter({f: e for f, e in fac.items() if e % 2}))
            else: rec["disc"] = str(disc)[:300]
        alg[L] = rec
    res["algebraic_letters"] = alg
    # (iii) rational t-letters: coefficient denominators; rational points and adjacent differences in words
    dens = collections.Counter(); ptset = {}; pairset = set(); denpolys = {}
    for tm in doc["terms"]:
        for cs_ in tm["coef"]:
            D = fmpq_poly([fmpq(p, q) for p, q in cs_[2]])
            if D.degree() >= 1: k_ = str(D); dens[k_] += 1; denpolys.setdefault(k_, D)
        words = [w for w in tm["Z"]] + [a[1] for a in tm["const"] if a[0] == "G0"]
        for w in words:
            prev = None
            for p in w:
                key = json.dumps(p); ptset[key] = p
                if prev is not None and prev != key: pairset.add((prev, key) if prev < key else (key, prev))
                prev = key
    tfac = collections.Counter()
    # factor the distinct denominators with FLINT directly
    for Ds, cnt in dens.items():
        fp = denpolys[Ds]; cont, facs = fp.factor()
        for f, e in facs:
            if f.degree() < 1: continue
            co = [sp.Rational(int(f[i].p), int(f[i].q)) for i in range(f.degree() + 1)]
            fe = sp.Poly(sum(c * T ** i for i, c in enumerate(co)), T).primitive()[1]
            if fe.LC() < 0: fe = -fe
            tfac[str(fe.as_expr())] += 1
    res["coefficient_denominator_t_letters"] = {"n_distinct_denominators": len(dens), "irreducible_factors(count of denominators containing it)": dict(tfac), "table1_hits": table1_hits(tfac)}
    # rational points: norms of numerator and denominator; differences of adjacent rational points
    pfac = collections.Counter(); npt = 0
    def gt_of(p): return GT.deser(p[1])
    for key, p in ptset.items():
        if p[0] != "Q": continue
        g = gt_of(p)
        if g.is_zero(): continue
        npt += 1; ex = gt2sp(g); nrm = sp.cancel(sp.expand(ex * ex.subs(sp.I, -sp.I))) if not g.is_real() else ex
        for f, m in qfactors(nrm).items(): pfac[f] += 1
        onep = sp.cancel(1 + ex); nrm1 = sp.cancel(sp.expand(onep * onep.subs(sp.I, -sp.I))) if not g.is_real() else onep
        for f, m in qfactors(nrm1).items(): pfac[f] += 1
    dfac = collections.Counter(); nd = 0
    for k1, k2 in pairset:
        p1, p2 = ptset[k1], ptset[k2]
        if p1[0] != "Q" or p2[0] != "Q": continue
        d = gt2sp(gt_of(p1)) - gt2sp(gt_of(p2)); d = sp.cancel(d)
        if d == 0: continue
        nd += 1; nrm = sp.cancel(sp.expand(d * d.subs(sp.I, -sp.I)))
        for f, m in qfactors(nrm).items(): dfac[f] += 1
    res["rational_point_t_letters"] = {"n_rational_points": npt, "factors_of_norms_of_points_and_1+point": dict(pfac), "table1_hits_points": table1_hits(pfac), "n_adjacent_rational_pairs": nd, "factors_of_norms_of_differences": dict(dfac), "table1_hits_differences": table1_hits(dfac)}
    # (iv) per algebraic letter: net contribution per weight at the gate t's
    if do_S:
        S = {}
        for L in alg:
            def filt(tm, L=L):
                if any(g[0] == L for g in tm["gens"]): return True
                for w in tm["Z"]:
                    if any(p[0] == "A" and p[1] == L for p in w): return True
                for a in tm["const"]:
                    if a[0] == "G0" and any(p[0] == "A" and p[1] == L for p in a[1]): return True
                return False
            nL = sum(1 for tm in doc["terms"] if filt(tm)); row = {"n_terms_with_letter": nL, "per_t": {}}
            for tv in tvals:
                t1 = time.time(); tot, byw = assemble.evaluate_file(path, Fr(tv), dps, term_filter=filt, by_weight=True, doc=doc)
                row["per_t"][tv] = {"S_L_total": mp.nstr(tot, 20), "abs": mp.nstr(abs(tot), 3), "by_weight": {str(w): mp.nstr(abs(v), 3) for w, v in sorted(byw.items())}, "wall_s": round(time.time() - t1, 1)}
            S[L] = row; print("   letter deg %d %s: terms %d, |S_L|(t) = %s" % (alg[L]["deg_x"], alg[L]["poly"][:50], nL, {tv: row["per_t"][tv]["abs"] for tv in tvals}), flush=True)
        res["net_contribution_of_terms_with_algebraic_letter"] = S
    return res
def cmsyz_cubic_compare(functions):
    """our non-linear X-letters vs the CMSYZ cubic P restricted to the slice (ckpt/cmsyz_cubics_slice.json: 6 distinct polynomials over the 24
    labelings, 4 irreducible cubics and 2 that factor as (eta+1) x quadratic): exact polynomial identity up to eta in {X, 1/X, -X, -1/X} and an
    overall constant."""
    p = slice_prov.rpath("ckpt", "cmsyz_cubics_slice.json")          # optional input (format: {"cubics": [{"sigmas": [...], "factor": "poly in X, t"}, ...]}); WORK copy else DATA
    if not os.path.exists(p): return {"note": "cmsyz_cubics_slice.json not present (optional comparison skipped)"}
    cz = json.load(open(p))["cubics"]; X = sp.Symbol("X"); ETA = sp.Symbol("eta")
    cfacs = []
    for c in cz:
        expr = sp.sympify(c["P_slice(eta,t)"]); cc, facs = sp.factor_list(expr)
        for f, e in facs:
            Pf = sp.Poly(f, ETA, T)
            if Pf.degree(ETA) >= 2: cfacs.append({"sigmas": c["sigmas"], "factor": str(f), "deg_eta": Pf.degree(ETA), "full": c["P_slice(eta,t)"]})
    def canon(poly2):   # primitive integer poly in (X,t), sign-normalized
        Pp = sp.Poly(poly2, X, T); Pp = Pp.primitive()[1]
        if Pp.LC() < 0: Pp = -Pp
        return Pp
    out = {"cmsyz_factors_deg>=2": cfacs, "matches": {}}
    for lab, r in functions.items():
        for Lk, v in r["algebraic_letters"].items():
            items = json.loads(Lk); ours = sum(sp.Rational(pq[0], pq[1]) * X ** ex * T ** et * sp.I ** ej for (ex, et, ej), pq in items)
            if ours.has(sp.I): continue
            Po = canon(ours); hits = []
            for cf in cfacs:
                f = sp.sympify(cf["factor"]); d = cf["deg_eta"]
                for name, sub in (("eta=X", X), ("eta=1/X", 1 / X), ("eta=-X", -X), ("eta=-1/X", -1 / X)):
                    g = sp.expand(sp.cancel(f.subs(ETA, sub) * (X ** d if "1/X" in name else 1)))
                    try: Pg = canon(g)
                    except Exception: continue
                    if Pg == Po: hits.append({"cmsyz_sigmas": cf["sigmas"], "map": name, "cmsyz_factor": cf["factor"][:120]})
            if hits: out["matches"].setdefault(v["poly"], {"deg_x": v["deg_x"], "in_functions": [], "cmsyz": hits})["in_functions"].append(lab)
            else: out["matches"].setdefault(v["poly"], {"deg_x": v["deg_x"], "in_functions": [], "cmsyz": []})["in_functions"].append(lab)
    return out
def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--files", default="quark,gluon,n4"); ap.add_argument("--t", default="1/2,4/5"); ap.add_argument("--dps", type=int, default=60); ap.add_argument("--noS", action="store_true")
    a = ap.parse_args(); slice_prov.log_pid("letters"); t0 = time.time(); out = {"table1": {k: str(v) for k, v in TABLE1.items()}, "functions": {}}
    files = {"quark": "E4C_LO_QCD_dipole_slice_quark.json.gz", "gluon": "E4C_LO_QCD_dipole_slice_gluon.json.gz", "n4": "E4C_LO_N4_dipole_slice_n4.json.gz"}
    for ch_ in ("q_qbpqpgq", "q_qbqgq", "q_gggq", "g_qbpqpqbq", "g_qbqqbq", "g_qbggq", "g_gggg"): files[ch_] = "E4C_LO_QCD_dipole_slice_%s.json.gz" % ch_
    for lab in a.files.split(","):
        p = os.path.join(OUT, files[lab])
        if not os.path.exists(p): print("missing", p); continue
        print("==", lab, flush=True); out["functions"][lab] = analyze_file(p, lab, tuple(a.t.split(",")), a.dps, do_S=not a.noS)
    # shared / only tables
    sets = {lab: set((v["poly"]) for v in r["algebraic_letters"].values()) for lab, r in out["functions"].items()}
    if sets:
        allp = set().union(*sets.values())
        out["algebraic_letter_membership"] = {p_: sorted(lab for lab in sets if p_ in sets[lab]) for p_ in sorted(allp)}
    tl = {lab: set(r["coefficient_denominator_t_letters"]["irreducible_factors(count of denominators containing it)"]) | set(r["rational_point_t_letters"]["factors_of_norms_of_points_and_1+point"]) | set(r["rational_point_t_letters"]["factors_of_norms_of_differences"]) for lab, r in out["functions"].items()}
    if tl:
        allt = set().union(*tl.values()); out["rational_t_letter_membership"] = {p_: sorted(lab for lab in tl if p_ in tl[lab]) for p_ in sorted(allt)}
    out["cmsyz_compare"] = cmsyz_cubic_compare(out["functions"]); out["wall_s"] = round(time.time() - t0, 1)
    slice_prov.write_receipt(os.path.join(WORK, "SLICE_S5_LETTERS.json"), out, inputs=[os.path.join(OUT, f) for f in files.values() if os.path.exists(os.path.join(OUT, f))])
    print("wrote SLICE_S5_LETTERS.json")
if __name__ == "__main__":
    main()
