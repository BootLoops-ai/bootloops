"""Structure census of a slice stage-2 one-fold object F(x;t) (input of stage 3): pole letters in x (degree, multiplicity, root fields,
real-positive roots on 0<t<1?), the alphabet letters that enter the hyperlog word letters a, 1+a and differences a-b (their roots are the
fibration points), polynomial-part degrees, and on-contour crossings of word letters for x in (0,oo) at sample t."""
import os, sys, json, argparse, itertools
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import slice_prov, slice_common as sc, alphabet_gi, fibr_gi
from fractions import Fraction as Fr
import mpmath as mp, sympy as sp
from flint import fmpq, fmpq_mpoly_ctx
from fibr import ZERO_HL
WORK = slice_prov.WORK
def load_engine(tag):
    C, letters, ts, extra = sc.ts_load(os.path.join(WORK, "ckpt", tag, "stage2.pkl"))
    AL = alphabet_gi.GIAlphabet(C, cache_path=os.path.join(WORK, "ckpt", "gaussian_split_cache_%s.json" % "".join(C.names()[:3])))
    assert len(AL.letters) == len(C.names())
    for i, Lp in enumerate(letters):
        if i < len(C.names()): assert str(Lp) == str(AL.letters[i]); continue
        j = AL.add(Lp, "pickle"); assert j == i, ("alphabet index drift", i, j)
    E = fibr_gi.Engine(AL)
    return C, AL, E, ts, extra
def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--tag", required=True); ap.add_argument("--tsample", default="1/2,4/5,1/5"); a = ap.parse_args()
    C, AL, E, ts, extra = load_engine(a.tag); names = list(C.names()); last = "x%d" % extra["order"][2]; jx = names.index(last); jt = names.index("t"); jj = names.index("j")
    def xdeg(i): return int(AL.letters[i].degrees()[jx])
    # pole letters
    poles = {}
    polydeg = -10 ** 9
    for (Dk, Wk), N in ts.items():
        dx = int(N.degrees()[jx]); dd = 0
        for i, e in Dk:
            if xdeg(i) > 0:
                p = poles.setdefault(i, {"deg_x": xdeg(i), "deg_t": int(AL.letters[i].degrees()[jt]), "has_j": int(AL.letters[i].degrees()[jj]) > 0, "max_mult": 0, "poly": str(AL.letters[i])[:120]}); p["max_mult"] = max(p["max_mult"], e); dd += e * xdeg(i)
        polydeg = max(polydeg, dx - dd)
    # word letters -> alphabet letters entering a, 1+a, a-b
    hls = set(h for (Dk, Wk) in ts for wd in Wk for h in wd)
    entering = {}
    def note(i, why):
        if xdeg(i) > 0: entering.setdefault(i, set()).add(why)
    one = E.hl(Fr(1), {})
    for h in hls:
        if h == ZERO_HL: continue
        for i, e in E.hl_E(h).items(): note(i, "a")
        c1, E1 = E.hl_diff_mono(h, E.hl(Fr(-1), {}))      # a - (-1) = 1 + a
        for i, e in E1.items(): note(i, "1+a")
    words = set(wd for (Dk, Wk) in ts for wd in Wk)
    for wd in words:
        for h1, h2 in itertools.combinations([h for h in wd if h != ZERO_HL], 2):
            if h1 == h2: continue
            cd, Ed = E.hl_diff_mono(h1, h2)
            for i, e in Ed.items(): note(i, "a-b(same word)")
    # differences across words that get shuffled together (products of two weight-1 words)
    for (Dk, Wk), N in ts.items():
        if len(Wk) == 2:
            for h1 in Wk[0]:
                for h2 in Wk[1]:
                    if ZERO_HL in (h1, h2) or h1 == h2: continue
                    cd, Ed = E.hl_diff_mono(h1, h2)
                    for i, e in Ed.items(): note(i, "a-b(shuffle)")
    pts = {}
    for i, why in entering.items():
        pts[i] = {"deg_x": xdeg(i), "deg_t": int(AL.letters[i].degrees()[jt]), "has_j": int(AL.letters[i].degrees()[jj]) > 0, "enters_via": sorted(why), "poly": str(AL.letters[i])[:120], "also_pole": i in poles}
    # numeric root census at sample t: real positive roots?
    xs = sp.Symbol("x")
    def roots_at(i, tv):
        P = AL.letters[i]; d = P.to_dict(); poly = 0
        for e, cc in d.items():
            term = sp.Rational(int(cc.p), int(cc.q)) * xs ** int(e[jx]) * sp.Rational(tv) ** int(e[jt]) * sp.I ** int(e[jj])
            poly += term
        poly = sp.Poly(sp.expand(poly), xs)
        return [complex(r) for r in sp.polys.polytools.nroots(poly, n=30)]
    rootinfo = {}
    for i in sorted(set(poles) | set(pts)):
        ri = {}
        for tv_ in a.tsample.split(","):
            rr = roots_at(i, Fr(tv_)); ri[tv_] = {"roots": ["%.6g%+.6gi" % (r.real, r.imag) for r in rr], "real_positive_root": any(abs(r.imag) < 1e-20 and r.real > 0 for r in rr)}
        rootinfo[str(AL.letters[i])[:60]] = ri
    # on-contour word letters: sign of a(x) on a grid (real part positive & imaginary ~ 0)
    import slice_eval
    oncont = {}
    for tv_ in a.tsample.split(","):
        S = slice_eval.Specialized(C, AL.letters, ts, Fr(tv_), last); cnt = 0; ex = []
        for h in hls:
            if h == ZERO_HL: continue
            hit = []
            for xv in [Fr(k, 8) for k in range(1, 17)] + [Fr(3), Fr(5), Fr(9), Fr(20), Fr(60)]:
                v = S.hl_value(h, mp.mpf(xv.numerator) / xv.denominator)
                if abs(v.imag) < 1e-25 and v.real > 0: hit.append(str(xv))
            if hit: cnt += 1; ex.append((E.hl_str(h)[:70], hit[:6]))
        oncont[tv_] = {"n_word_letters_on_contour_somewhere_on_grid": cnt, "n_word_letters": len(hls), "examples": ex[:8]}
    rec = {"tag": a.tag, "last_var": last, "n_keys": len(ts), "pole_letters_x": {str(AL.letters[i])[:60]: v for i, v in poles.items()}, "max_polynomial_part_degree(deg_x N - deg_x D)": polydeg,
           "fibration_point_letters(x-letters entering word letters a, 1+a, a-b)": {str(AL.letters[i])[:60]: v for i, v in pts.items()},
           "summary": {"pole_degrees": sorted(set(v["deg_x"] for v in poles.values())), "point_degrees": sorted(set(v["deg_x"] for v in pts.values())), "max_pole_mult_by_degree": {d: max(v["max_mult"] for v in poles.values() if v["deg_x"] == d) for d in set(v["deg_x"] for v in poles.values())},
                       "nonlinear_pole_letters": [k for k, v in poles.items() if v["deg_x"] >= 2], "nonlinear_point_letters": [k for k, v in pts.items() if v["deg_x"] >= 2]},
           "roots_at_sample_t": rootinfo, "on_contour_word_letters": oncont}
    rec["summary"]["nonlinear_pole_letters"] = [str(AL.letters[i])[:80] for i in rec["summary"]["nonlinear_pole_letters"]]; rec["summary"]["nonlinear_point_letters"] = [str(AL.letters[i])[:80] for i in rec["summary"]["nonlinear_point_letters"]]
    rec["summary"]["any_real_positive_pole_root_at_samples"] = any(ri[tv_]["real_positive_root"] for k, ri in rootinfo.items() for tv_ in ri)
    out = os.path.join(WORK, "SLICE_A2_%s.json" % a.tag); slice_prov.write_receipt(out, rec, inputs=[os.path.join(WORK, "ckpt", a.tag, "stage2.pkl")])
    print(json.dumps(rec["summary"], indent=1)); print("polypart deg", polydeg); print(json.dumps(oncont, indent=0)[:1500]); print(json.dumps({k: {t_: v_["real_positive_root"] for t_, v_ in v.items()} for k, v in rootinfo.items()}, indent=0)[:3000])
if __name__ == "__main__":
    main()
