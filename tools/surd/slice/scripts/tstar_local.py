"""FUNCTION-level behavior of an assembled slice closed form (out/*.json.gz) at the square-root letter classes, on the PHYSICAL sheet.
For each square class C = odd Q[t]-factors of disc_X of the quadratic X-letters (as symtest.py), tau = swap of the two root slots of every
letter in C (a field automorphism over Q(i)(t, other roots)).  With S_C(t) = sum of the terms that carry a root of a C-letter (all other terms
are tau-fixed), the decomposition F = E + O, E = (S_C + S_C^tau)/2 + rest, O = (S_C - S_C^tau)/2 = sqrt(Delta_C) x B(t) is exact, and
   F analytic at a simple root t* of Delta_C on this sheet  <=>  O == 0 near t* (hence on the whole real segment of this sheet) and E bounded.
Numerically (hpath evaluator of assemble.py with permuted slots; -i0 rule for on-axis letters as everywhere in this package):
  (i)  O(t)/scale at the gate points t = 1/2, 4/5 (scale = max |term|): digits of vanishing;
  (ii) at t = t* -+ h, h = 1e-2..1e-5 (rational t), both sides where inside (0,1): O/scale, |Im S_C|, boundedness of E_C = (S_C+S_C^tau)/2,
       the termwise pole order (slope of log max|term| vs log h; individual coefficients carry inverse powers of the root difference
       delta = sqrt(Delta)/a, i.e. the capacity for the Landau exponent -3/2 = delta^-3), and the slope of log|O| if O != 0.
Table-1 letters (letters.TABLE1, the tangency letters of the dipole slice) are flagged per class with their roots; classes whose real roots all lie outside (0,1] are related to
their (13)(24)-images (t -> 1/t) and reported as such.  Writes SLICE_S5_TSTAR_<label>.json."""
import os, sys, json, gzip, time, argparse, collections, math, resource
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import slice_prov, hpath, letters as LT
from fractions import Fraction as Fr
import mpmath as mp, sympy as sp
from k3field import GT
WORK = slice_prov.WORK; DATA = slice_prov.DATA; OUT = os.path.join(DATA, "out"); T = sp.Symbol("t"); X0 = Fr(11, 2)
FILES = {"quark": "E4C_LO_QCD_dipole_slice_quark.json.gz", "gluon": "E4C_LO_QCD_dipole_slice_gluon.json.gz", "n4": "E4C_LO_N4_dipole_slice_n4.json.gz"}
for ch_ in ("q_qbpqpgq", "q_qbqgq", "q_gggq", "g_qbpqpqbq", "g_qbqqbq", "g_qbggq", "g_gggg"): FILES[ch_] = "E4C_LO_QCD_dipole_slice_%s.json.gz" % ch_

def gt2sp(ser):
    g = GT.deser(ser); f = lambda P: sum(sp.Rational(int(c.p), int(c.q)) * T ** k for k, c in enumerate(P.coeffs()))
    return (f(g.A) + sp.I * f(g.B)) / f(g.D)
def sqclass(expr):
    num, den = sp.fraction(sp.cancel(sp.together(expr))); ker = []
    for pol in (num, den):
        P = sp.Poly(pol, T)
        if P.degree() <= 0:
            if sp.Rational(P.as_expr()) < 0: ker.append("-1")
            continue
        cst, facs = sp.factor_list(P)
        if sp.Rational(cst) < 0: ker.append("-1")
        for f_, e in facs:
            if int(e) % 2: ker.append(str(sp.Poly(f_, T).primitive()[1].as_expr()))
    return tuple(sorted(ker))
def classes_of(doc):
    quad = {}
    for L, info in doc["letters"].items():
        if info["deg_x"] != 2: continue
        cs = [gt2sp(c) for c in info["x_coeffs(GT ser, increasing powers of X)"]]
        real = all(GT.deser(c).is_real() for c in info["x_coeffs(GT ser, increasing powers of X)"])
        disc = sp.cancel(cs[1] ** 2 - 4 * cs[2] * cs[0]); quad[L] = {"poly": info["poly"], "real": real, "disc": disc, "class": sqclass(disc) if real else ("GAUSSIAN", str(disc)[:80])}
    cl = collections.defaultdict(list)
    for L, q in quad.items(): cl[q["class"]].append(L)
    return quad, cl

def evaluate(doc, terms, tval, dps, swap=frozenset()):
    """sum of `terms` of doc at rational t; roots of letters in `swap` (degree 2) with slots exchanged. Returns (total, max|term|, sum |Im| not needed)"""
    tval = Fr(tval); mp.mp.dps = dps; H = hpath.HPath(dps); tol = mp.mpf(10) ** (-(dps - 8)); roots = {}
    def rts(L):
        if L not in roots:
            cs = [GT.deser(c).mp(tval, mp) for c in doc["letters"][L]["x_coeffs(GT ser, increasing powers of X)"]]
            with mp.workdps(dps + 40): rr = mp.polyroots(list(reversed(cs)), maxsteps=500, extraprec=4 * dps)
            rr = [mp.mpc(r) for r in rr]
            if L in swap: assert len(rr) == 2; rr = [rr[1], rr[0]]
            roots[L] = rr
        return roots[L]
    def ptval(p):
        if p[0] == "Q":
            g = GT.deser(p[1])
            if g.is_real():
                re, im = g(tval); return (re, "below" if re > 0 else None)
            v = g.mp(tval, mp)
            if abs(v.imag) < tol * (1 + abs(v.real)) and v.real > 0: return (mp.mpc(v.real, 0), "below")
            return (v, None)
        v = rts(p[1])[p[2]]
        if abs(v.imag) < tol * (1 + abs(v.real)): return (mp.mpc(v.real, 0), "below" if v.real > 0 else None)
        return (v, None)
    x0m = mp.mpf(X0.numerator) / X0.denominator
    def hl_at_x0(h):
        (p, q), E = h; v = mp.mpc(mp.mpf(p) / q)
        if not E: return Fr(p, q)
        for L, e in E:
            items = json.loads(L); val = mp.mpc(0)
            for (ex, et, ej), (pp, qq) in items: val += (mp.mpf(pp) / qq) * x0m ** ex * (mp.mpf(tval.numerator) / tval.denominator) ** et * mp.mpc(0, 1) ** ej
            v = v * val ** e
        return v
    cZ = {}; cC = {}
    def Zval(w):
        k = json.dumps(w)
        if k not in cZ:
            lets = []
            for p in w:
                if p[0] == "Q" and GT.deser(p[1]).is_zero(): lets.append((Fr(0), None)); continue
                lets.append(ptval(p))
            cZ[k] = H.zip_value(lets)
        return cZ[k]
    def Cval(a):
        k = json.dumps(a)
        if k in cC: return cC[k]
        if a[0] == "Z0":
            lets = []
            for h in a[1]:
                lv = hl_at_x0(((h[0][0], h[0][1]), tuple((x[0], x[1]) for x in h[1])))
                if isinstance(lv, Fr): lets.append((lv, "below" if lv > 0 else None))
                elif abs(lv.imag) < tol * (1 + abs(lv.real)): lets.append((mp.mpc(lv.real, 0), "below" if lv.real > 0 else None))
                else: lets.append((lv, None))
            v = H.zip_value(lets)
        else:
            ZER = ("E", 0); toks = tuple(ZER if (p[0] == "Q" and GT.deser(p[1]).is_zero()) else ("P", json.dumps(p)) for p in a[1])
            tot = mp.mpc(0); lx = mp.log(x0m)
            for (k0, u), c in hpath.decomp_trailing(toks, ZER).items():
                lets = []
                for tk in u:
                    if tk == ZER: lets.append((Fr(0), None)); continue
                    val, sd = ptval(json.loads(tk[1])); lets.append(((val / X0) if isinstance(val, Fr) else val / x0m, sd))
                g = H.G1(lets) if lets else 1
                tot += (mp.mpf(c.numerator) / c.denominator) * lx ** k0 / math.factorial(k0) * g
            v = tot
        cC[k] = v; return v
    tot = mp.mpc(0); mx = mp.mpf(0); byw = collections.defaultdict(lambda: mp.mpc(0))
    for tm in terms:
        c = mp.mpc(0)
        for cs_ in tm["coef"]: c += GT.deser(cs_).mp(tval, mp)
        for L, s, e in tm["gens"]: c *= rts(L)[s] ** e
        for a_ in tm["const"]: c *= Cval(a_)
        for w in tm["Z"]: c *= Zval(w)
        tot += c; mx = max(mx, abs(c)); byw[sum(len(w) for w in tm["Z"]) + sum(len(a_[1]) for a_ in tm["const"])] += c
    return tot, mx, dict(byw)

def frac_near(x, den=10 ** 14): return Fr(int(round(float(x) * den)), den) if not isinstance(x, Fr) else x
def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--label", required=True); ap.add_argument("--dps", type=int, default=60); ap.add_argument("--gate-t", dest="gt", default="1/2,4/5")
    ap.add_argument("--hs", default="1e-2,1e-3,1e-4,1e-5"); ap.add_argument("--classes", default="all", help="'table1' = only classes containing a Table-1 factor")
    a = ap.parse_args(); slice_prov.log_pid("tstar_local:" + a.label); t00 = time.time()
    path = os.path.join(OUT, FILES.get(a.label, a.label)); doc = json.load(gzip.open(path, "rt")); quad, cl = classes_of(doc)
    T1 = {nm: sp.Poly(pol, T).primitive()[1].as_expr() for nm, pol in LT.TABLE1.items()}
    out = {"label": a.label, "file": path, "file_sha256": slice_prov.sha(path), "complete": doc.get("complete"), "method": __doc__, "classes": {}}
    for cls, Ls in cl.items():
        if cls and cls[0] == "GAUSSIAN": out["classes"][str(cls)] = {"letters": [quad[L]["poly"] for L in Ls], "note": "non-real quadratic letter: not tested"}; continue
        facs = [f for f in cls if f != "-1"]
        t1hits = [nm for nm, pe in T1.items() if any(sp.simplify(sp.sympify(f) - pe) == 0 or sp.simplify(sp.sympify(f) + pe) == 0 for f in facs)]
        if a.classes == "table1" and not t1hits: continue
        Lset = frozenset(Ls)
        def involves(tm):
            if any(g[0] in Lset for g in tm["gens"]): return True
            for w in tm["Z"]:
                if any(p[0] == "A" and p[1] in Lset for p in w): return True
            for c_ in tm["const"]:
                if c_[0] == "G0" and any(p[0] == "A" and p[1] in Lset for p in c_[1]): return True
            return False
        # class factor polynomials (flint) for the exact denominator test: a term whose coefficient denominator D(t) is divisible by a class
        # factor carries a pole at its roots and belongs to the local partial sum S'_C even if it carries no class root
        from flint import fmpq_poly as _fp, fmpq as _fq
        fpolys = []      # only the Table-1 factors of the class enter the denominator criterion (t, t+-1, t^2+1 divide most denominators and are not letters under test)
        for f in facs:
            P_ = sp.Poly(sp.sympify(f), T)
            if not any(sp.simplify(P_.primitive()[1].as_expr() - pe) == 0 or sp.simplify(P_.primitive()[1].as_expr() + pe) == 0 for pe in T1.values()): continue
            fpolys.append((f, _fp([_fq(int(sp.Rational(c_).p), int(sp.Rational(c_).q)) for c_ in reversed(P_.all_coeffs())])))
        def den_order(tm):
            """max multiplicity of any class factor in any coefficient denominator of the term"""
            m = 0
            for cs_ in tm["coef"]:
                D = _fp([_fq(p_, q_) for p_, q_ in cs_[2]])
                if D.degree() < 1: continue
                for f, fp_ in fpolys:
                    k = 0; Q_ = D
                    while Q_.degree() >= fp_.degree():
                        qq, rr = divmod(Q_, fp_)
                        if rr != 0: break
                        k += 1; Q_ = qq
                    m = max(m, k)
            return m
        terms = []; n_root = 0; n_den_only = 0; max_den = 0; den_hist = collections.Counter()
        for tm in doc["terms"]:
            iv = involves(tm); m = den_order(tm)
            if iv or m > 0:
                terms.append(tm); n_root += int(iv); n_den_only += int((not iv) and m > 0); max_den = max(max_den, m); den_hist[m] += 1
        rec = {"letters": [quad[L]["poly"] for L in Ls], "sqrt_class(odd Q[t] factors of disc_X)": list(cls), "table1_letters": t1hits, "n_terms_with_class_roots": n_root, "n_terms_denominator_only": n_den_only, "n_terms_local_partial_sum": len(terms),
               "max_class_factor_multiplicity_in_coefficient_denominators": max_den, "denominator_multiplicity_histogram(terms of S'_C)": dict(den_hist),
               "weight_histogram": dict(collections.Counter(sum(len(w) for w in tm["Z"]) + sum(len(c_[1]) for c_ in tm["const"]) for tm in terms))}
        # real roots of the class factors in (0, 1]
        roots = []
        for f in facs:
            P = sp.Poly(sp.sympify(f), T)
            for r in P.nroots(n=40):
                if abs(sp.im(r)) < 1e-30 and 0 < sp.re(r) <= 1: roots.append((f, sp.re(r)))
        rec["real_roots_in_(0,1]"] = [{"factor": f, "t*": str(r)[:22]} for f, r in roots]
        if not roots:
            imgs = []
            for f in facs:
                P = sp.Poly(sp.sympify(f), T); d = P.degree(); img = sp.Poly(sp.expand(T ** d * P.as_expr().subs(T, 1 / T)), T).primitive()[1]
                imgs.append(str(img.as_expr()))
            rec["note"] = "no real root in (0,1]; (13)(24)-image factors (t -> 1/t): %s" % imgs
        print("== class %s  Table-1 %s  letters %d  terms %d (root %d, den-only %d, max den mult %d)  roots %s" % (list(cls), t1hits, len(Ls), len(terms), n_root, n_den_only, max_den, [str(r)[:8] for f, r in roots]), flush=True)
        # (i) gate points
        rows = []
        def probe(tv, tag):
            t1_ = time.time(); F, mx, bw = evaluate(doc, terms, tv, a.dps); Ft, mxt, bwt = evaluate(doc, terms, tv, a.dps, swap=Lset)
            O = (F - Ft) / 2; E = (F + Ft) / 2; sc_ = max(mx, mxt)
            row = {"tag": tag, "t": str(tv), "t_float": "%.16g" % float(tv), "S_C": mp.nstr(F, 25), "S_C^tau": mp.nstr(Ft, 25), "abs_O": mp.nstr(abs(O), 5), "abs_E": mp.nstr(abs(E), 8), "scale(max|term|)": mp.nstr(sc_, 5),
                   "O/scale": mp.nstr(abs(O) / sc_, 5) if sc_ > 0 else None, "digits_O_vanishes(-log10 |O|/scale)": (round(float(-mp.log10(abs(O) / sc_)), 1) if abs(O) > 0 else 99.0) if sc_ > 0 else None,
                   "abs_Im_S_C": mp.nstr(abs(F.imag), 5), "O_by_weight(|.|)": {str(w): mp.nstr(abs((bw.get(w, 0) - bwt.get(w, 0)) / 2), 5) for w in sorted(set(bw) | set(bwt))}, "wall_s": round(time.time() - t1_, 1)}
            rows.append(row); print("   %s t=%s: |O|/scale %s (%s d) |E| %s scale %s |Im S| %s (%.0fs)" % (tag, row["t_float"], row["O/scale"], row["digits_O_vanishes(-log10 |O|/scale)"], row["abs_E"], row["scale(max|term|)"], row["abs_Im_S_C"], time.time() - t1_), flush=True)
            return row
        for tv in a.gt.split(","): probe(Fr(tv), "gate")
        # (ii) near each root
        local = []
        for f, r in roots:
            rstar = float(r); side_rows = {"+": [], "-": []}
            for hs in a.hs.split(","):
                h = float(hs)
                for sgn in ("+", "-"):
                    tv = rstar + h if sgn == "+" else rstar - h
                    if not (0 < tv < 1): continue
                    row = probe(frac_near(tv), "t*%s%s [%s]" % (sgn, hs, f)); row["h"] = h; row["side"] = sgn; row["root_of"] = f; side_rows[sgn].append(row)
            fit = {}
            for sgn, lst in side_rows.items():
                if len(lst) >= 2:
                    hs_ = [math.log10(x["h"]) for x in lst]; sc_ = [math.log10(float(mp.mpf(x["scale(max|term|)"]))) for x in lst]; Es = [float(mp.mpf(x["abs_E"])) for x in lst]
                    slope_sc = (sc_[-1] - sc_[0]) / (hs_[-1] - hs_[0])
                    Os = [float(mp.mpf(x["abs_O"])) for x in lst]
                    slope_O = ((math.log10(Os[-1]) - math.log10(Os[0])) / (hs_[-1] - hs_[0])) if all(o > 0 for o in Os) else None
                    erel = (abs(Es[-1] - Es[-2]) / abs(Es[-2])) if Es[-2] != 0 else None
                    fit[sgn] = {"termwise_slope dlog(max|term|)/dlog h": round(slope_sc, 3), "termwise_pole_order_in_delta=sqrt(Delta) (= -2 slope)": round(-2 * slope_sc, 2), "E_values(h decreasing)": ["%.12g" % e for e in Es],
                                "E_rel_change_last_two_h": erel, "E_bounded(<0.1 rel change)": (erel is not None and erel < 0.1), "slope dlog|O|/dlog h": (round(slope_O, 3) if slope_O is not None else None), "min digits_O_vanishes": min(x["digits_O_vanishes(-log10 |O|/scale)"] for x in lst)}
            local.append({"factor": f, "t*": str(r)[:22], "fits": fit})
        rec["gate_rows"] = [r_ for r_ in rows if r_["tag"] == "gate"]; rec["local_rows"] = [r_ for r_ in rows if r_["tag"] != "gate"]; rec["local_fits"] = local
        dig = [r_["digits_O_vanishes(-log10 |O|/scale)"] for r_ in rows if r_.get("digits_O_vanishes(-log10 |O|/scale)") is not None]
        eb = all(all(ft.get("E_bounded(<0.1 rel change)") for ft in x["fits"].values()) for x in local) if local else None
        rec["E_bounded_all_roots"] = eb
        rec["verdict_physical_sheet"] = ("sqrt(Delta_C)-odd part O == 0 at all %d probe points (min %.1f digits below the termwise scale): closed form tau-EVEN (=> F meromorphic in t at the class roots on this sheet); local partial sum S'_C (class-root terms + terms with the class factor in a coefficient denominator, max multiplicity %d) bounded as h -> 0 at every probed root: %s%s" % (len(dig), min(dig), max_den, eb, " => analytic at t*" if eb else (" (S'_C grows: either a pole of F at t* or other letters degenerate at the same t*; see rows)" if eb is False else ""))) if dig and min(dig) >= a.dps // 2 - 5 else "O != 0 at some probe point (min digits %s): representation NOT tau-even (inconsistency)" % (min(dig) if dig else None)
        out["classes"][str(cls)] = rec; print("   verdict:", rec["verdict_physical_sheet"], flush=True)
    out["wall_s"] = round(time.time() - t00, 1); out["peak_rss_MB"] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1)
    slice_prov.write_receipt(os.path.join(WORK, "SLICE_S5_TSTAR_%s.json" % a.label), out, inputs=[path])
    print("wrote SLICE_S5_TSTAR_%s.json %.0fs" % (a.label, time.time() - t00))
if __name__ == "__main__":
    main()
