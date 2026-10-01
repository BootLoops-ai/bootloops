"""Stage 3 (function level) for one slice group: I(t) = int_0^oo dx F(x;t), F = stage-2 one-fold object (ckpt/<tag>/stage2.pkl).
  A. fibration of every ZIP-word product into G_0(points; x) with constants at the base point x0 (fib3.Fibrator);
  B. exact integration of every (rational function) x (basis word) with algebraic partial fractions (int3.Integrator);
Output = sum over terms  coef(t, roots) * CONST * prod Z(words), stored exactly in ckpt/<tag>/stage3.pkl:
  coef  : element of the tensor ring TRing over Q(i)(t) in the root generators (letter, slot) (k3field), serialized;
  CONST : product of ZIP[w(x0)] (t-only hyperlog letters at x = x0) and G_0(v; x0) (points scaled by 1/x0);
  Z(u)  : reg_{oo} G_0(u; oo) over points (ZIP semantics), letters = 0, rational points in Q(i)(t), roots of the nonlinear x-letters.
Then numeric evaluation at the gate t's (roots by mpmath polyroots, hyperlogs by hpath with the universal -i0 rule for on-axis letters)
divided by the external factor prod d^ext(t), against the per-group piece oracle (slice_po) -> writes SLICE_S3_<tag>.json."""
import os, sys, json, time, argparse, pickle, resource, collections
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import slice_prov, slice_common as sc, analyze_s2, slice_eval, fib3, int3, k3field, hpath, slice_po
from fractions import Fraction as Fr
import mpmath as mp
from k3field import GT, TRing, Points
from fibr import ZERO_HL
WORK = slice_prov.WORK
def rss_mb(): return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0

def ser_point(p): return ["Q", p[1].ser()] if p[0] == "Q" else ["A", int(p[1]), int(p[2])]
def deser_point(s): return ("Q", GT.deser(s[1])) if s[0] == "Q" else ("A", s[1], s[2])
def ser_hl(h): return [[h[0][0], h[0][1]], [[int(i), int(e)] for i, e in h[1]]]
def deser_hl(s): return ((s[0][0], s[0][1]), tuple((i, e) for i, e in s[1]))
def ser_atom(a): return ["Z0", [ser_hl(h) for h in a[1]]] if a[0] == "Z0" else ["G0", [ser_point(p) for p in a[1]]]
def deser_atom(s): return ("Z0", tuple(deser_hl(h) for h in s[1])) if s[0] == "Z0" else ("G0", tuple(deser_point(p) for p in s[1]))

class Stage3:
    def __init__(self, tag, x0=Fr(11, 2), log=print):
        self.tag = tag; self.x0 = Fr(x0); self.log = log
        self.C, self.AL, self.E, self.ts, self.extra = analyze_s2.load_engine(tag)
        self.names = list(self.C.names()); self.last = "x%d" % self.extra["order"][2]
        self.jx, self.jt, self.jj = self.names.index(self.last), self.names.index("t"), self.names.index("j")
        self.F = fib3.Fibrator(self.E, self.jx, self.jt, self.jj, self.x0); self.P = self.F.P; self.R = TRing(self.P); self.I = int3.Integrator(self.R)
        self.terms = {}; self.stats = collections.Counter()
    def xdeg(self, i): return int(self.AL.letters[i].degrees()[self.jx])
    def run(self):
        R = self.R; t0 = time.time(); nk = 0
        keys = sorted(self.ts.keys(), key=repr)
        for (Dk, Wk) in keys:
            N = self.ts[(Dk, Wk)]; nk += 1; t1 = time.time()
            fp = self.F.fib(tuple(Wk))                                  # {v: {constkey: Fr}}
            cx = GT(1); xf = []; pts_of = {}
            for i, e in Dk:
                if self.xdeg(i) == 0: cx = cx * GT.from_tj_poly(self.AL.letters[i], self.jt, self.jj) ** (-e)
                else:
                    pts = self.F.letter_points(i); pts_of[len(xf)] = (i, pts); xf.append((self.P.L[i], e))
            Ncoef = k3field.xcoeffs_GT(N, self.jx, self.jt, self.jj)
            polypart, poles = self.I.rational_parts(Ncoef, xf, pts_of)
            self.stats["keys"] += 1; self.stats["poles"] += len(poles); self.stats["polypart_max_deg"] = max(self.stats["polypart_max_deg"], len(polypart) - 1)
            for v, cv in fp.items():
                B = self.I.integrate_word(polypart, poles, v)
                if not B: continue
                for ck, q in cv.items():
                    sc_ = cx * q
                    for zk, te in B.items():
                        key = (ck, zk); w = self.terms.get(key); te2 = R.scal(te, sc_)
                        if w is None: self.terms[key] = te2
                        else:
                            w = R.add(w, te2)
                            if w: self.terms[key] = w
                            else: del self.terms[key]
            if self.log and (nk % 10 == 0 or time.time() - t1 > 30): self.log("  key %d/%d (weight %d, %d poles, %d fib words): %.1fs; terms %d; RSS %.0f MB" % (nk, len(keys), sum(len(w) for w in Wk), len(poles), len(fp), time.time() - t1, len(self.terms), rss_mb()))
        self.stats["wall_AB_s"] = round(time.time() - t0, 1); self.stats["terms"] = len(self.terms)
        return self.terms
    # ---- census of the exact result
    def census(self):
        zw = set(); consts = set(); gens = set(); maxdeg = [0, 0]; ncomp = 0; wts = collections.Counter(); onlett = collections.Counter()
        for (ck, zk), te in self.terms.items():
            for a in ck: consts.add(a)
            for w in zk: zw.add(w)
            wt = sum(len(w) for w in zk) + sum((len(a[1]) if a[0] == "G0" else len(a[1])) for a in ck); wts[wt] += 1
            for k, v in te.items():
                ncomp += 1
                for g, e in k: gens.add(g)
                a_, b_ = v.degs(); maxdeg[0] = max(maxdeg[0], a_); maxdeg[1] = max(maxdeg[1], b_)
        pts = set(p for w in zw for p in w) | set(p for a in consts if a[0] == "G0" for p in a[1])
        letters_in_points = sorted(set(p[1] for p in pts if p[0] == "A"))
        return {"n_terms": len(self.terms), "n_coefficient_components": ncomp, "max_t_degree_coefficients(num,den)": maxdeg, "n_distinct_Zwords": len(zw), "n_distinct_constants": len(consts),
                "Zword_weights": dict(collections.Counter(len(w) for w in zw)), "total_transcendental_weight_histogram(terms)": dict(wts),
                "generators(letter,slot)": sorted([list(g) for g in gens]), "nonlinear_letters_in_points": {str(self.AL.letters[i])[:90]: {"idx": i, "deg_x": self.P.deg[i], "has_j": int(self.AL.letters[i].degrees()[self.jj]) > 0} for i in letters_in_points},
                "n_rational_points": len([p for p in pts if p[0] == "Q"]), "n_algebraic_points": len([p for p in pts if p[0] == "A"]), "x0": str(self.x0)}
    def load(self, path):
        d = pickle.load(open(path, "rb"))
        assert d["ctx_names"] == self.names and len(d["alphabet_dicts"]) >= len(self.AL.letters)
        from flint import fmpq
        for i, lst in enumerate(d["alphabet_dicts"]):
            Lp = self.C.from_dict({tuple(e): fmpq(p, q) for e, p, q in lst})
            if i < len(self.AL.letters): assert Lp == self.AL.letters[i], ("alphabet drift at", i)
            else:
                j = self.AL.add(Lp, "stage3pickle"); assert j == i
        for i, cs in d["letters_x_coeffs"].items():
            self.P.L[int(i)] = [GT.deser(c) for c in cs]; self.P.deg[int(i)] = len(self.P.L[int(i)]) - 1
        self.terms = {}
        for ck, zk, te in d["terms"]:
            key = (tuple(deser_atom(x) for x in ck), tuple(tuple(deser_point(p) for p in w) for w in zk)); self.terms[key] = TRing.deser(te)
        self.stats.update(d.get("stats", {})); return self.terms
    def save(self, path):
        out = {"tag": self.tag, "x0": [self.x0.numerator, self.x0.denominator], "last_var": self.last, "ctx_names": self.names, "alphabet": [str(L) for L in self.AL.letters],
               "alphabet_dicts": [[(list(int(z) for z in e), int(c.p), int(c.q)) for e, c in L.to_dict().items()] for L in self.AL.letters],
               "letters_x_coeffs": {int(i): [c.ser() for c in cs] for i, cs in self.P.L.items()},
               "terms": [[[ser_atom(a) for a in ck], [[ser_point(p) for p in w] for w in zk], TRing.ser(te)] for (ck, zk), te in self.terms.items()],
               "ext_factor_poly_t": self.extra["ext_factor_poly_t"], "extra": self.extra, "stats": dict(self.stats)}
        tmp = path + ".tmp%d" % os.getpid()
        with open(tmp, "wb") as fh: pickle.dump(out, fh, protocol=4)
        os.replace(tmp, path); return path
    # ---- numeric evaluation
    def evaluate(self, tval, dps=50):
        tval = Fr(tval); mp.mp.dps = dps
        N = fib3.NumCtx(self.P, None, tval, dps); gv = N.genval()
        S = slice_eval.Specialized(self.C, self.AL.letters, {}, tval, self.last)
        x0m = mp.mpf(self.x0.numerator) / self.x0.denominator
        cache_c = {}; cache_z = {}
        def atom_val(a):
            if a in cache_c: return cache_c[a]
            if a[0] == "Z0":
                lets = []
                for h in a[1]:
                    lv = S.hl_value(h, x0m)
                    if isinstance(lv, Fr): lets.append((lv, bool(lv > 0)))
                    elif abs(lv.imag) < N.tol * (1 + abs(lv.real)): lets.append((mp.mpc(lv.real, 0), bool(lv.real > 0)))
                    else: lets.append((lv, False))
                v = N.H.zip_value(lets)
            else: v = N.G0(a[1], self.x0)
            cache_c[a] = v; return v
        def zval(w):
            if w in cache_z: return cache_z[w]
            v = N.Z(w); cache_z[w] = v; return v
        tot = mp.mpc(0); mx = mp.mpf(0)
        for (ck, zk), te in self.terms.items():
            c = TRing.evalf(te, gv, tval, mp)
            for a in ck: c *= atom_val(a)
            for w in zk: c *= zval(w)
            tot += c; mx = max(mx, abs(c))
        ef = self.extra["ext_factor_poly_t"]
        G = sc.build_group(self.extra["channel"], self.extra["gauge"], self.extra["pair"], tuple(int(c) for c in self.extra["sigma"]))
        efv = sc.ext_factor_value(G, tval)
        return tot / (mp.mpf(efv.numerator) / efv.denominator), mx, {"hpath": dict(N.H.stats), "n_oncontour_points": sum(1 for w in cache_z for p in w if N.ptval(p)[1] if not (p[0] == "Q" and p[1].is_zero()))}

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--tag", required=True); ap.add_argument("--x0", default="11/2"); ap.add_argument("--t", default="1/2,4/5"); ap.add_argument("--dps", type=int, default=70)
    ap.add_argument("--no-po", dest="po", action="store_false"); ap.add_argument("--reuse", action="store_true")
    a = ap.parse_args(); slice_prov.log_pid("stage3:" + a.tag); t00 = time.time()
    logf = open(os.path.join(slice_prov.wdir("logs"), "stage3_%s.log" % a.tag), "a")
    def log(s):
        l = "%s %s" % (time.strftime("%H:%M:%S"), s); print(l, flush=True); logf.write(l + "\n"); logf.flush()
    S3 = Stage3(a.tag, Fr(a.x0), log=log)
    log("%s: %d stage-2 keys, last var %s, x0=%s" % (a.tag, len(S3.ts), S3.last, a.x0))
    pkl = os.path.join(WORK, "ckpt", a.tag, "stage3.pkl")
    if a.reuse and os.path.exists(pkl): S3.load(pkl); log("loaded %s (%d terms)" % (pkl, len(S3.terms)))
    else: S3.run(); S3.save(pkl)
    cen = S3.census()
    log("A+B done: %s ; %.1fs RSS %.0f MB" % ({k: v for k, v in cen.items() if k not in ("generators(letter,slot)", "nonlinear_letters_in_points")}, S3.stats.get("wall_AB_s"), rss_mb()))
    rows = []
    for ts_ in a.t.split(","):
        tv = Fr(ts_); t1 = time.time(); val, mx, est = S3.evaluate(tv, a.dps); te = time.time() - t1
        row = {"t": ts_, "closed_form_value/ext": mp.nstr(val, 42), "abs_Im": mp.nstr(abs(val.imag), 3), "max_term": mp.nstr(mx, 5), "eval_wall_s": round(te, 1), "eval_stats": est}
        if a.po:
            po = slice_po.piece_oracle(S3.extra["channel"], S3.extra["sigma"], S3.extra["pair"], tv); ref = mp.mpf(po["routes"]["A"]["value_mid"])
            rel = abs(val - ref) / abs(ref) if ref != 0 else abs(val - ref); dg = 99.0 if rel == 0 else float(-mp.log10(rel))
            row.update({"piece_oracle_A": po["routes"]["A"]["value"], "piece_oracle_A_vs_B": po["A_vs_B_rel"], "agree_digits": round(dg, 1), "po_wall_s": po["wall_s"]})
        rows.append(row); log("t=%s: closed form %s | oracle %s | digits %s (eval %.0fs)" % (ts_, mp.nstr(val, 25), row.get("piece_oracle_A", "-")[:30], row.get("agree_digits"), te))
    rec = {"tag": a.tag, "channel": S3.extra["channel"], "sigma": S3.extra["sigma"], "pair": S3.extra["pair"], "gauge": S3.extra["gauge"], "last_var": S3.last, "x0": a.x0, "census": cen, "stats": dict(S3.stats),
           "gate_rows": rows, "min_agree_digits": min((r["agree_digits"] for r in rows if "agree_digits" in r), default=None), "PASS(>=30)": bool(rows and all(r.get("agree_digits", 0) >= 30 for r in rows)),
           "stage3_pickle": pkl, "wall_s": round(time.time() - t00, 1), "peak_rss_MB": round(rss_mb(), 1)}
    out = os.path.join(WORK, "SLICE_S3_%s.json" % a.tag)
    slice_prov.write_receipt(out, rec, inputs=[os.path.join(WORK, "ckpt", a.tag, "stage2.pkl"), pkl] + [os.path.join(HERE, m) for m in ("fib3.py", "int3.py", "k3field.py", "hpath.py", "fibr_gi.py", "alphabet_gi.py")])
    log("wrote %s min digits %s PASS %s wall %.0fs RSS %.0f MB" % (out, rec["min_agree_digits"], rec["PASS(>=30)"], time.time() - t00, rss_mb()))
if __name__ == "__main__":
    main()
