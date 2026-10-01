"""Assembly: F_ch(t) = w_flavor(nf=5) * S_ch * sum_{sigma classes} mult * sum_{first-pair groups} I_group(t),  I_group = (stage-3 closed form)/ext(t),
and the jets F_q = sum quark channels, F_g = sum gluon channels (normalization: the xL^3 = (max d_ab)^3 prefactor of oracle.py replaced by d34^3 = 1, i.e. the xL^3-normalized value = max(1,t^2)^3 F).
All group term stores are mapped to CANONICAL keys (letter polynomials as strings in (X = last variable, t, j), points = (poly, slot) or exact
Q(i)(t) values, hyperlog words over points, constants ZIP[w(x0)] / G_0(v;x0)) and identical hyperlogs are collected with exact Q(i)(t) coefficients.
Outputs out/E4C_LO_QCD_dipole_slice_{quark,gluon}.json (+ per channel, + n4) in the documented format, and a numeric evaluator of the files."""
import os, sys, json, time, argparse, pickle, collections, hashlib, gzip
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import slice_prov, slice_common as sc, k3field, hpath
from fractions import Fraction as Fr
import mpmath as mp
from flint import fmpq, fmpq_poly, fmpq_mpoly_ctx
from k3field import GT
import oracle_cells as OC
WORK = slice_prov.WORK; OUT = os.path.join(WORK, "out")
X0 = Fr(11, 2)

def poly_key(dct, names, last):
    """canonical string of a letter polynomial: monomials in (X, t, j) sorted"""
    jx, jt, jj = names.index(last), names.index("t"), names.index("j")
    items = sorted(((int(e[jx]), int(e[jt]), int(e[jj])), (int(p), int(q))) for e, p, q in dct)
    return json.dumps(items, separators=(",", ":"))
def poly_str(key):
    items = json.loads(key); terms = []
    for (ex, et, ej), (p, q) in items:
        c = Fr(p, q); mono = "*".join([s for s, e in (("X^%d" % ex, ex), ("t^%d" % et, et), ("j^%d" % ej, ej)) if e]).replace("^1*", "*")
        if mono.endswith("^1"): mono = mono[:-2]
        terms.append(("%s" % c) + ("*" + mono if mono else ""))
    return " + ".join(terms).replace("+ -", "- ")

class Assembler:
    def __init__(self):
        self.terms = collections.defaultdict(list)               # (const, Z, gens) -> list of GT summands (one per contributing group component)
        self.letters = {}                                        # polykey -> {"deg_x", "x_coeffs": [GT ser]}
        self.groups = []
    def add_group(self, tag, scale_gt, meta):
        d = pickle.load(open(os.path.join(WORK, "ckpt", tag, "stage3.pkl"), "rb"))
        names = d["ctx_names"]; last = d["last_var"]; x0 = Fr(*d["x0"]); assert x0 == X0, "base point differs"
        pk = {}
        for i, lst in enumerate(d["alphabet_dicts"]):
            key = poly_key(lst, names, last); pk[i] = key
            if i in map(int, d["letters_x_coeffs"].keys()) or str(i) in d["letters_x_coeffs"]:
                cs = d["letters_x_coeffs"].get(i) or d["letters_x_coeffs"].get(str(i))
                if key not in self.letters: self.letters[key] = {"deg_x": len(cs) - 1, "x_coeffs": cs, "poly": poly_str(key)}
        ext = fmpq_poly([fmpq(0)]); 
        # external factor polynomial (string of fmpq_poly in x) -> GT
        import sympy as sp
        e = sp.Poly(sp.sympify(d["ext_factor_poly_t"].replace("^", "**")), sp.Symbol("x")); cs = [sp.Rational(c) for c in reversed(e.all_coeffs())]
        extgt = GT(fmpq_poly([fmpq(int(c.p), int(c.q)) for c in cs]))
        scale = scale_gt / extgt
        def cpt(p):
            return ("Q", tuple(tuple(tuple(x) for x in comp) for comp in p[1])) if p[0] == "Q" else ("A", pk[p[1]], p[2])
        def chl(h):
            return ((h[0][0], h[0][1]), tuple(sorted(((pk[i], e) for i, e in h[1]))))
        n = 0
        for ck, zk, te in d["terms"]:
            cat = tuple(sorted(((("Z0", tuple(chl(h) for h in a[1])) if a[0] == "Z0" else ("G0", tuple(cpt(p) for p in a[1]))) for a in ck), key=repr))
            zz = tuple(sorted((tuple(cpt(p) for p in w) for w in zk), key=repr))
            for mono, gser in te:
                gens = tuple(sorted(((pk[g0], g1, ex) for g0, g1, ex in mono)))
                c = GT.deser(gser) * scale
                key = (cat, zz, gens); self.terms[key].append(c); n += 1
        self.groups.append({"tag": tag, "n_components": n, "scale": repr(scale_gt), **meta})
        return n
    def finalize(self):
        """collect identical hyperlog monomials: coefficients stay lists of Q(i)(t) summands (compact); a key is dropped only if its exact sum vanishes"""
        out = {}; self.n_cancelled = 0
        for k, lst in self.terms.items():
            if len(lst) > 1:
                tot = GT()
                for c in lst: tot = tot + c
                if tot.is_zero(): self.n_cancelled += 1; continue
                a_, b_ = tot.degs()
                if a_ + b_ <= sum(sum(c.degs()) for c in lst): lst = [tot]        # keep the summed form only when it is not larger
            out[k] = lst
        self.terms = out
    def census(self):
        Z = set(); C = set(); gens = set(); md = [0, 0]; wt = collections.Counter()
        ncomp = 0
        for (cat, zz, g), lst in self.terms.items():
            for w in zz: Z.add(w)
            for a in cat: C.add(a)
            for pk_, s, e in g: gens.add((pk_, s))
            for c in lst:
                ncomp += 1; a_, b_ = c.degs(); md[0] = max(md[0], a_); md[1] = max(md[1], b_)
            wt[sum(len(w) for w in zz) + sum(len(a[1]) for a in cat)] += 1
        pts = set(p for w in Z for p in w) | set(p for a in C if a[0] == "G0" for p in a[1])
        alg = sorted(set(p[1] for p in pts if p[0] == "A"))
        return {"n_terms(distinct hyperlog monomials)": len(self.terms), "n_coefficient_summands": ncomp, "n_monomials_cancelled_exactly_between_groups": getattr(self, "n_cancelled", None), "n_distinct_hyperlogs_Z": len(Z), "n_distinct_constants": len(C), "n_points": len(pts), "n_algebraic_point_letters": len(alg), "algebraic_point_letters": [{"poly": poly_str(k), "deg_x": self.letters.get(k, {}).get("deg_x")} for k in alg],
                "n_root_generators": len(gens), "max_t_degree_coefficients(num,den)": md, "weight_histogram(terms)": dict(sorted(wt.items())), "n_groups": len(self.groups)}
    def dump(self, path, header):
        def sp_(p): return ["Q", [list(map(list, comp)) for comp in p[1]]] if p[0] == "Q" else ["A", p[1], p[2]]
        terms = []
        for (cat, zz, g), lst in self.terms.items():
            terms.append({"coef": [c.ser() for c in lst], "gens": [list(x) for x in g],
                          "const": [["Z0", [[list(h[0]), [list(x) for x in h[1]]] for h in a[1]]] if a[0] == "Z0" else ["G0", [sp_(p) for p in a[1]]] for a in cat],
                          "Z": [[sp_(p) for p in w] for w in zz]})
        doc = {"format": {
            "value": "F(t) = sum_terms coef(t) * prod_{(L,s,e) in gens} rho_{L,s}(t)^e * prod const * prod_{w in Z} Z(w)",
            "coef": "LIST of summands, each [A, B, D] = coefficient lists (increasing powers of t) of (A(t) + i B(t))/D(t) in Q(i)(t); the coefficient is their sum",
            "gens": "[letter polynomial key L, slot s, exponent e]: rho_{L,s} = the s-th root in X of the letter polynomial L(X; t, j=i) (any fixed bijection slots<->roots, used consistently everywhere)",
            "points": "['Q', [A,B,D]] = the exact point (A+iB)/D in Q(i)(t);  ['A', L, s] = rho_{L,s}",
            "Z": "Z(w) = ZIP_reg[w] = shuffle-regularized G(w_1,...,w_n; oo) from 0 (HyperFLINT ZeroInfPeriod semantics: log-divergences at 0 and oo dropped); = sum over letters a of {1 (coef -1)} + {a/(1+a) (coef +1, absent if a=-1)} of G(...;1)",
            "const": "['Z0', word] = ZIP_reg of the word of hyperlog letters c*prod L^e (rational functions of X,t,j) evaluated at X = x0 = %s;  ['G0', v] = G_0(v; x0) = iterated integral from 0 to x0 with letters = points (G(0;x0) = log x0)" % X0,
            "branches": "every letter that lies exactly on the positive real integration path is 'letter - i0' (path passes above); all other hyperlogs on their principal branch",
            "letters": "letter polynomial keys are JSON lists of [[deg_X, deg_t, deg_j], [p, q]] monomials; 'letters' gives their X-coefficients as Q(i)(t) elements"},
            **header, "letters": {k: {"poly": v["poly"], "deg_x": v["deg_x"], "x_coeffs(GT ser, increasing powers of X)": v["x_coeffs"]} for k, v in self.letters.items()},
            "census": self.census(), "groups": self.groups, "terms": terms}
        path = path + ".gz" if not path.endswith(".gz") else path
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True); tmp = path + ".tmp%d" % os.getpid()
        with gzip.open(tmp, "wt", compresslevel=6) as fh: json.dump(doc, fh)
        os.replace(tmp, path)
        return path

# ---------------------------------------------------------------- numeric evaluation of an assembled file (independent of the group pickles)
def evaluate_file(path, tval, dps=50, term_filter=None, by_weight=False, doc=None):
    doc = doc or json.load(gzip.open(path, "rt") if path.endswith(".gz") else open(path)); tval = Fr(tval); mp.mp.dps = dps; H = hpath.HPath(dps); tol = mp.mpf(10) ** (-(dps - 8))
    roots = {}
    def rts(L):
        if L not in roots:
            cs = [GT.deser(c).mp(tval, mp) for c in doc["letters"][L]["x_coeffs(GT ser, increasing powers of X)"]]
            with mp.workdps(dps + 40): rr = mp.polyroots(list(reversed(cs)), maxsteps=400, extraprec=4 * dps)
            roots[L] = [mp.mpc(r) for r in rr]
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
    x0 = X0; x0m = mp.mpf(x0.numerator) / x0.denominator
    def hl_at_x0(h):
        (p, q), E = h; v = mp.mpc(mp.mpf(p) / q)
        if not E: return Fr(p, q)
        for L, e in E:
            items = json.loads(L); val = mp.mpc(0)
            for (ex, et, ej), (pp, qq) in items:
                val += (mp.mpf(pp) / qq) * x0m ** ex * (mp.mpf(tval.numerator) / tval.denominator) ** et * mp.mpc(0, 1) ** ej
            v = v * val ** e
        return v
    cZ = {}; cC = {}
    import math
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
            ZER = ("E", 0); toks = tuple(ZER if (p[0] == "Q" and GT.deser(p[1]).is_zero()) else ("P", tuple(map(lambda z: z if not isinstance(z, list) else json.dumps(z), p)), p) for p in a[1])
            toks = tuple(ZER if (p[0] == "Q" and GT.deser(p[1]).is_zero()) else ("P", json.dumps(p)) for p in a[1])
            tot = mp.mpc(0); lx = mp.log(x0m)
            for (k0, u), c in hpath.decomp_trailing(toks, ZER).items():
                lets = []
                for tk in u:
                    if tk == ZER: lets.append((Fr(0), None)); continue
                    val, sd = ptval(json.loads(tk[1])); lets.append(((val / x0) if isinstance(val, Fr) else val / x0m, sd))
                g = H.G1(lets) if lets else 1
                tot += (mp.mpf(c.numerator) / c.denominator) * lx ** k0 / math.factorial(k0) * g
            v = tot
        cC[k] = v; return v
    tot = mp.mpc(0); byw = {}
    for tm in doc["terms"]:
        if term_filter is not None and not term_filter(tm): continue
        c = mp.mpc(0)
        for cs_ in tm["coef"]: c += GT.deser(cs_).mp(tval, mp)
        for L, s, e in tm["gens"]: c *= rts(L)[s] ** e
        for a in tm["const"]: c *= Cval(a)
        for w in tm["Z"]: c *= Zval(w)
        tot += c
        if by_weight:
            wgt = sum(len(w) for w in tm["Z"]) + sum(len(a[1]) for a in tm["const"]); byw[wgt] = byw.get(wgt, mp.mpc(0)) + c
    return (tot, byw) if by_weight else tot

CHW = {ch: OC.flavor_weight(ch, 5) * OC.SFAC[ch] for ch in OC.QUARK_JET + OC.GLUON_JET + ["n4"]}
def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--channels", default="all"); ap.add_argument("--eval-t", dest="evalt", default=""); ap.add_argument("--dps", type=int, default=50)
    a = ap.parse_args(); slice_prov.log_pid("assemble"); t0 = time.time()
    CL = json.load(open(os.path.join(WORK, "SLICE_CLASSES.json")))["channels"]
    chans = (OC.QUARK_JET + OC.GLUON_JET + ["n4"]) if a.channels == "all" else a.channels.split(",")
    report = {"channels": {}, "jets": {}}
    jets = {"quark": Assembler(), "gluon": Assembler()}
    for ch in chans:
        A = Assembler(); missing = []; used = []
        for cl in CL[ch]["classes"]:
            for pair, gi in CL[ch]["first_pair_groups"].items():
                tag = "%s_s%s_p%s_g%d" % (ch, cl["rep"], pair, gi["gauge"])
                if not os.path.exists(os.path.join(WORK, "ckpt", tag, "stage3.pkl")) or not os.path.exists(os.path.join(WORK, "SLICE_S3_%s.json" % tag)): missing.append(tag); continue
                sc_ = GT(Fr(cl["mult"]) * CHW[ch])
                A.add_group(tag, sc_, {"channel": ch, "sigma": cl["rep"], "mult": cl["mult"], "pair": pair}); used.append(tag)
                if ch in OC.QUARK_JET: jets["quark"].add_group(tag, sc_, {"channel": ch, "sigma": cl["rep"], "mult": cl["mult"], "pair": pair})
                if ch in OC.GLUON_JET: jets["gluon"].add_group(tag, sc_, {"channel": ch, "sigma": cl["rep"], "mult": cl["mult"], "pair": pair})
        A.finalize()
        hdr = {"object": "LO collinear E4C on the dipole slice, channel %s, F_ch(t) = w_flavor(nf=5)*S_ch*sum_sigma I_sigma(d(t)); xL^3-normalized value = max(1,t^2)^3 * F" % ch, "channel": ch, "flavor_weight_x_sfac": str(CHW[ch]), "complete": not missing, "missing_groups": missing, "x0": str(X0)}
        path = A.dump(os.path.join(OUT, "E4C_LO_%s_dipole_slice_%s.json" % ("N4" if ch == "n4" else "QCD", ch)), hdr)
        cen = A.census(); report["channels"][ch] = {"file": path, "size_MB": round(os.path.getsize(path) / 1e6, 2), "complete": not missing, "n_missing": len(missing), "n_groups": len(used), "census": cen}
        print(ch, "groups", len(used), "missing", len(missing), {k: v for k, v in cen.items() if k != "algebraic_point_letters"}, flush=True)
    for jn, A in jets.items():
        chs = OC.QUARK_JET if jn == "quark" else OC.GLUON_JET
        if not any(c in chans for c in chs): continue
        A.finalize(); miss = [m for c in chs for m in report["channels"].get(c, {}).get("n_missing", 0) * [c]]
        complete = all(report["channels"].get(c, {}).get("complete") for c in chs)
        hdr = {"object": "LO collinear E4C on the dipole slice, %s jet (n_f=5): F(t) = sum_channels w_flavor*S_ch*sum_sigma I_sigma(d(t)); d12=t^2, d34=1, d13=d24=t^2/4+3t/10+1/4, d14=d23=t^2/4-3t/10+1/4; xL^3-normalized (oracle.py) value = max(1,t^2)^3 F(t); F(1/t) = t^6 F(t)" % jn, "jet": jn, "channels": chs, "complete": complete, "x0": str(X0)}
        path = A.dump(os.path.join(OUT, "E4C_LO_QCD_dipole_slice_%s.json" % jn), hdr); cen = A.census()
        report["jets"][jn] = {"file": path, "size_MB": round(os.path.getsize(path) / 1e6, 2), "complete": complete, "census": cen}
        with open(os.path.join(OUT, "E4C_LO_QCD_dipole_slice_%s.txt" % jn), "w") as f:
            f.write("LO collinear E4C, %s jet, dipole slice (cos phi = 3/5), exact closed form: see the .json (format block inside).\ncomplete: %s\n" % (jn, complete))
            for k, v in cen.items(): f.write("%s: %s\n" % (k, v if k != "algebraic_point_letters" else ""))
            f.write("algebraic (non-linear in X) point letters:\n")
            for L in cen["algebraic_point_letters"]: f.write("  deg %s : %s\n" % (L["deg_x"], L["poly"]))
        print(jn, "jet:", {k: v for k, v in cen.items() if k != "algebraic_point_letters"}, "file", path, report["jets"][jn]["size_MB"], "MB", flush=True)
    if a.evalt:
        for ts_ in a.evalt.split(","):
            for jn in jets:
                p = report["jets"].get(jn, {}).get("file")
                if p: t1 = time.time(); v = evaluate_file(p, Fr(ts_), a.dps); report["jets"][jn].setdefault("values", {})[ts_] = mp.nstr(v, 40); print(jn, "t=%s" % ts_, mp.nstr(v, 35), "(%.0fs)" % (time.time() - t1), flush=True)
    report["wall_s"] = round(time.time() - t0, 1)
    import resource; report["peak_rss_MB"] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1)
    # merge with an existing receipt (separate invocations per jet keep memory bounded); entries of this run replace older ones
    rp = os.path.join(WORK, "SLICE_ASSEMBLY.json")
    if os.path.exists(rp) and a.channels != "all":
        try:
            old = json.load(open(rp)); old.pop("producer", None)
            for k in ("channels", "jets"):
                for kk, vv in old.get(k, {}).items():
                    if kk not in report[k] and os.path.exists(vv.get("file", "")): vv["from_earlier_invocation"] = old.get("runs", [{}])[-1].get("launch") if old.get("runs") else True; report[k][kk] = vv
            report["runs"] = old.get("runs", []) + [{"launch": " ".join(sys.argv), "wall_s": report["wall_s"], "peak_rss_MB": report["peak_rss_MB"], "stamp_utc": slice_prov.stamp()}]
        except Exception as ex: report["merge_note"] = repr(ex)
    for k in ("channels", "jets"):
        for kk, vv in report[k].items():
            if "file" in vv and os.path.exists(vv["file"]) and "sha256" not in vv: vv["sha256"] = slice_prov.sha(vv["file"])
    slice_prov.write_receipt(rp, report, inputs=[os.path.join(WORK, "SLICE_CLASSES.json")] + [v["file"] for k in ("channels", "jets") for v in report[k].values() if "file" in v and os.path.exists(v["file"])])
if __name__ == "__main__":
    main()
