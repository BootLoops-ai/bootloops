"""N=4 slice closed form vs the CMSYZ slice alphabet.  CMSYZ side: data/N4_SLICE.json (the 63 GPL/q letters of the arXiv:2401.06463 ancillary files
restricted to the dipole slice over the 12 labeling classes: 34 irreducible Q[t]-norm factors, root letters via the cubic resultant) and the cubic
P|slice per labeling class (ckpt/cmsyz_cubics_slice.json).  Our side (SLICE_S5_LETTERS.json, function label n4 and optionally quark/gluon):
the Q[t]-prime support of the representation = irreducible factors of (i) coefficient denominators, (ii) norms of rational fibration points p and
1+p, (iii) norms of adjacent point differences inside hyperlog words, (iv) discriminant square classes of the quadratic X-letters and
discriminants of the cubic X-letters.  NOTE: (i)-(iii) are the prime support of ONE representation (fibration basis, last variable of each
group, base point x0 = 11/2 excluded by construction here: Z0/G_0(.;x0) constants are not scanned), a superset of the function's symbol
alphabet; the comparison is therefore 'support vs support', stated as such.  Writes SLICE_S5_N4ALPHABET.json."""
import os, sys, json, argparse
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import slice_prov
import sympy as sp
WORK = slice_prov.WORK; T = sp.Symbol("t")
# N4_SLICE.json: $SURD_N4SLICE if set, else the copy shipped in data/
N4S = os.environ.get("SURD_N4SLICE") or os.path.join(slice_prov.ROOT, "data", "N4_SLICE.json")
def canon(expr):
    P = sp.Poly(sp.sympify(str(expr).replace("^", "**")), T)
    if P.degree() <= 0: return None
    P = P.primitive()[1]
    if P.LC() < 0: P = -P
    return str(P.as_expr())
def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--labels", default="n4"); a = ap.parse_args(); slice_prov.log_pid("n4_alphabet")
    LTP = slice_prov.rpath("SLICE_S5_LETTERS.json")
    if not os.path.exists(LTP): raise SystemExit(f"n4_alphabet: {LTP} not found; run letters.py first (it writes SLICE_S5_LETTERS.json under the work directory)")
    if not os.path.exists(N4S): raise SystemExit(f"n4_alphabet: {N4S} not found; set SURD_N4SLICE or keep data/N4_SLICE.json in the package")
    LT = json.load(open(LTP)); cz = json.load(open(N4S))          # WORK copy if letters.py ran there, else DATA
    cm = {}
    for k, v in cz["letters"].items():
        c = canon(k)
        if c: cm[c] = {"cmsyz_letter": k, "degree": v["degree"], "disc": v.get("disc"), "from_q": v.get("from_q"), "n_sources": v.get("n_sources")}
    out = {"method": __doc__, "cmsyz_slice_letters(n)": len(cm), "cmsyz_source": N4S, "cmsyz_source_sha256": slice_prov.sha(N4S), "functions": {}}
    for lab in a.labels.split(","):
        r = LT["functions"].get(lab)
        if not r: continue
        sup = {}
        def add(f, origin):
            c = canon(f)
            if c: sup.setdefault(c, set()).add(origin)
        for f in r["coefficient_denominator_t_letters"]["irreducible_factors(count of denominators containing it)"]: add(f, "coef-denominator")
        for f in r["rational_point_t_letters"]["factors_of_norms_of_points_and_1+point"]: add(f, "norm(point) or norm(1+point)")
        for f in r["rational_point_t_letters"]["factors_of_norms_of_differences"]: add(f, "norm(adjacent point difference)")
        alg = []
        for L, v in r["algebraic_letters"].items():
            if v["deg_x"] == 2:
                for f in (v.get("sqrt_field_kernel(odd factors)") or []): add(f, "disc class of a quadratic X-letter")
                alg.append({"deg": 2, "poly": v["poly"][:120], "sqrt_class": v.get("sqrt_field_kernel(odd factors)") or v.get("disc_in_Q(i)(t)")})
            elif v["deg_x"] == 3:
                for f in (v.get("disc_square_class(odd factors)") or []): add(f, "disc class of a cubic X-letter")
                alg.append({"deg": 3, "poly": v["poly"][:120], "disc_square_class": v.get("disc_square_class(odd factors)")})
        ours = set(sup); theirs = set(cm)
        both = sorted(ours & theirs, key=lambda s: (len(s), s)); only_ours = sorted(ours - theirs, key=lambda s: (len(s), s)); only_cmsyz = sorted(theirs - ours, key=lambda s: (len(s), s))
        out["functions"][lab] = {"complete": r.get("complete"), "n_prime_support(ours)": len(ours), "in_both": [{"letter": x, "origins(ours)": sorted(sup[x]), "cmsyz": cm[x]} for x in both],
                                 "only_in_our_representation": [{"letter": x, "origins": sorted(sup[x])} for x in only_ours], "only_in_cmsyz_slice_alphabet": [{"letter": x, **cm[x]} for x in only_cmsyz],
                                 "algebraic_X_letters(ours)": alg, "counts": {"both": len(both), "only_ours": len(only_ours), "only_cmsyz": len(only_cmsyz)}}
        print(lab, out["functions"][lab]["counts"]); print("  both:", both); print("  only ours:", only_ours[:40]); print("  only CMSYZ:", only_cmsyz)
    slice_prov.write_receipt(os.path.join(WORK, "SLICE_S5_N4ALPHABET.json"), out, inputs=[LTP, N4S])
if __name__ == "__main__":
    main()
