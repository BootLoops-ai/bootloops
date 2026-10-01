"""Which Landau family produces each Table-1 slice letter in the closed forms.  For every group (SLICE_S12 receipt) and every quadratic
last-variable letter L(X;t) it used whose X-discriminant square class contains a Table-1 factor, search the group's denominator atoms
(slice_common.build_group: Om, x_i, ell_S, s_S composed with d_sigma(t), chart x_gauge = 1, elimination order of the receipt) for the chain
that generates L:  (a) '(1,1,2) tangency': L == s_S restricted to the line {A = B = 0} of two atoms linear in the eliminated variables
(A, B in {x_i, ell_S, Om, and the coordinate hyperplanes x_a = 0 / x_a = oo of the integration domain}); (b) generic two-step resultant chains
Res_{x_b}(F1, F2) with F1, F2 in {atoms free of x_a, leading/trailing x_a-coefficients, x_a-discriminants, Res_{x_a}(P, Q)} (irreducible factors).
Output per Table-1 class: channels, number of groups, the producing families (atom names) with counts.  Writes SLICE_S5_LANDAU.json."""
import os, sys, json, glob, time, argparse, collections, itertools
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import slice_prov, slice_common as sc, letters as LT
import sympy as sp
from flint import fmpq, fmpq_mpoly_ctx
WORK = slice_prov.WORK; T = sp.Symbol("t")
T1 = {nm: sp.Poly(pol, T).primitive()[1] for nm, pol in LT.TABLE1.items()}
def canon(p):
    """primitive, sign-normalized fmpq_mpoly (leading coefficient positive in the ctx term order)"""
    if p.is_zero(): return p
    cont, facs = p.factor()
    q = None
    for f, e in facs:
        q = f ** int(e) if q is None else q * f ** int(e)
    if q is None: return p.context().from_dict({(0,) * len(p.context().names()): fmpq(1)})
    lc = q.leading_coefficient() if hasattr(q, "leading_coefficient") else list(q.to_dict().values())[0]
    if lc < 0: q = -q
    return q
def irr_factors(p):
    if p.is_zero() or p.total_degree() <= 0: return []
    cont, facs = p.factor(); return [canon(f) for f, e in facs if f.total_degree() > 0]
def deg(p, v): return int(p.degrees()[p.context().names().index(v)]) if v in p.context().names() else 0
def coeffs_in(p, v):
    """coefficients of p as a polynomial in variable v: list index k -> fmpq_mpoly (same ctx, v-free)"""
    C = p.context(); i = C.names().index(v); out = {}
    for e, c in p.to_dict().items():
        k = e[i]; e2 = list(e); e2[i] = 0; out.setdefault(k, {}); out[k][tuple(e2)] = out[k].get(tuple(e2), fmpq(0)) + c
    return {k: C.from_dict(d) for k, d in out.items()}
def table1_class_of(Lpoly, X):
    """square class of disc_X(L) over Q[t]: list of odd primitive factors; Table-1 names hit"""
    cs = coeffs_in(Lpoly, X)
    if max(cs) != 2: return None, []
    C = Lpoly.context(); zero = C.from_dict({(0,) * len(C.names()): fmpq(0)})
    a, b, c = cs.get(2, zero), cs.get(1, zero), cs.get(0, zero); disc = b * b - 4 * a * c
    if disc.is_zero(): return None, []
    cont, facs = disc.factor(); odd = []
    if cont < 0: odd.append("-1")
    for f, e in facs:
        if int(e) % 2:
            ex = sp.Poly(sp.sympify(str(f).replace("^", "**")), T).primitive()[1]
            if ex.LC() < 0: ex = -ex
            odd.append(str(ex.as_expr()))
    hits = [nm for nm, P in T1.items() if any(o != "-1" and (sp.Poly(sp.sympify(o), T) - P).is_zero for o in odd)]
    return sorted(odd), hits
def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--channels", default="q_qbpqpgq,q_qbqgq,q_gggq,g_qbpqpqbq,g_qbqqbq,g_qbggq,g_gggg,n4")
    ap.add_argument("--all-sigma", dest="one_sigma", action="store_false", help="run the chain search for every sigma (default: one sigma per (channel, pair))"); ap.add_argument("--all-letters", dest="all_letters", action="store_true")
    a = ap.parse_args()
    slice_prov.log_pid("landau_family"); t0 = time.time()
    out = {"method": __doc__, "per_class": {}, "per_group": {}, "channels": a.channels.split(",")}
    recs = []
    for f in sorted(glob.glob(os.path.join(WORK, "SLICE_S12_*.json"))):
        r = json.load(open(f))
        if r["channel"] in out["channels"]: recs.append(r)
    percls = collections.defaultdict(lambda: {"table1": None, "channels": collections.Counter(), "groups": 0, "families": collections.Counter(), "unexplained_groups": []})
    prov_done = {}     # (channel, pair) -> True: the atom-chain search is run for ONE sigma per (channel, first-pair group); other sigma differ by the d-relabeling only
    for r in recs:
        nl = {k: v for k, v in r["nonlinear_last_variable_letters_used"].items() if v["deg_x"] == 2 and not v["has_j"]}
        if not nl: continue
        last = r["order"][2]; xa, xb = r["order"][0], r["order"][1]
        G = None; hits_here = []
        for Ls, info in nl.items():
            # parse letter into the group ctx lazily (need ctx names): build group once
            if G is None:
                G = sc.build_group(r["channel"], r["gauge"], r["pair"], tuple(int(c) for c in r["sigma"])); C = G["ctx"]; names = C.names()
                gens = dict(zip(names, C.gens()))
                def parse(s):
                    e = sp.Poly(sp.sympify(s.replace("^", "**")), *[sp.Symbol(n) for n in names]); d = {}
                    for mon, cf in e.terms(): d[tuple(int(m) for m in mon)] = fmpq(int(sp.Rational(cf).p), int(sp.Rational(cf).q))
                    return C.from_dict(d)
                atoms = dict(G["atoms"]); one = C.from_dict({(0,) * len(names): fmpq(1)})
                # coordinate hyperplanes of the domain as 'atoms' for the line search
                lin = {nm: p for nm, p in atoms.items() if all(int(dg) <= 1 for dg in p.degrees()[:3]) and p.total_degree() >= 1}
                for v in (xa, xb, last): lin["{%s=0}" % v] = gens[v]
                quads = {nm: p for nm, p in atoms.items() if nm.startswith("s")}
            Lp = canon(parse(Ls)); cls, hits = table1_class_of(Lp, last)
            if not hits: continue
            key = str(tuple(cls)); pc = percls[key]; pc["table1"] = hits; pc["channels"][r["channel"]] += 1; pc["groups"] += 1
            fams = []
            rep = prov_done.setdefault(("rep", r["channel"], r["pair"]), r["tag"])
            if a.one_sigma and rep != r["tag"]:
                hits_here.append({"letter": Ls[:100], "class": cls, "table1": hits, "families": "see representative sigma %s of this (channel, pair)" % rep}); continue
            # (a) (1,1,2): s_S restricted to {A=B=0}, A,B linear: solve A for xa (if present) then B for xb
            SA, SB_, SX = sp.Symbol(xa), sp.Symbol(xb), sp.Symbol(last)
            def tosp(p): return sp.sympify(str(p).replace("^", "**"))
            for (na, A), (nb, B) in itertools.combinations(sorted(lin.items()), 2):
                sol = sp.solve([tosp(A), tosp(B)], [SA, SB_], dict=True)
                if not sol or len(sol) != 1: continue
                sol = sol[0]
                if SA not in sol or SB_ not in sol: continue
                for ns, Sq in quads.items():
                    e = sp.together(tosp(Sq).subs(sol)); num, den = sp.fraction(e)
                    if num == 0: continue
                    Pn = sp.Poly(num, SX, T)
                    if Pn.degree(SX) != 2: continue
                    try: cand = canon(parse(str(Pn.as_expr())))
                    except Exception: continue
                    for fct in irr_factors(cand):
                        if fct == Lp: fams.append("(1,1,2): %s | line {%s = %s = 0}" % (ns, na, nb))
            # (b) generic resultant chains (only if (a) found nothing)
            if not fams:
                st1 = {}
                for nm, p in atoms.items():
                    if deg(p, xa) == 0: st1[nm] = [canon(p)]; continue
                    cs = coeffs_in(p, xa); parts = []
                    for k_ in (max(cs), min(cs)): parts += irr_factors(cs[k_])
                    if max(cs) == 2: parts += irr_factors(p.discriminant(xa))
                    st1["coef/disc_%s(%s)" % (xa, nm)] = parts
                for (n1, p1), (n2, p2) in itertools.combinations(sorted(atoms.items()), 2):
                    if deg(p1, xa) and deg(p2, xa): st1["Res_%s(%s,%s)" % (xa, n1, n2)] = irr_factors(p1.resultant(p2, xa))
                flat = [(k_, q) for k_, lst in st1.items() for q in lst if not q.is_zero() and q.total_degree() > 0]
                for (k1, q1), (k2, q2) in itertools.combinations(flat, 2):
                    if deg(q1, xb) and deg(q2, xb):
                        for fct in irr_factors(q1.resultant(q2, xb)):
                            if fct == Lp: fams.append("Res_%s[ %s , %s ]" % (xb, k1, k2))
                    for q in (q1,):
                        pass
                for k1, q1 in flat:
                    if deg(q1, xb) == 2:
                        for fct in irr_factors(q1.discriminant(xb)):
                            if fct == Lp: fams.append("disc_%s[ %s ]" % (xb, k1))
                    if deg(q1, xb) == 0 and q1 == Lp: fams.append("direct[ %s ]" % k1)
            fams = sorted(set(fams)); prov_done[(r["channel"], r["pair"], Ls)] = True
            for fm in fams: pc["families"][fm] += 1
            if not fams: pc["unexplained_groups"].append(r["tag"])
            hits_here.append({"letter": Ls[:100], "class": cls, "table1": hits, "families": fams})
        if hits_here: out["per_group"][r["tag"]] = {"channel": r["channel"], "order": r["order"], "letters": hits_here}; print(r["tag"], [(h["table1"], h["families"][:3]) for h in hits_here], flush=True)
    out["per_class"] = {k: {"table1": v["table1"], "channels(groups)": dict(v["channels"]), "n_groups": v["groups"], "families(count over groups x letters)": dict(v["families"]), "unexplained_groups": v["unexplained_groups"]} for k, v in percls.items()}
    out["wall_s"] = round(time.time() - t0, 1)
    print(json.dumps(out["per_class"], indent=1))
    slice_prov.write_receipt(os.path.join(WORK, "SLICE_S5_LANDAU.json"), out, inputs=sorted(glob.glob(os.path.join(WORK, "SLICE_S12_*.json"))))
if __name__ == "__main__":
    main()
