"""Stage 3, step A: fibration of the stage-2 ZIP-word products with x-dependent letters into the basis G_0(v; x) (iterated integrals from 0,
letters = fibration POINTS = roots in x of alphabet letters, incl. 0) with x-independent constants fixed at a real base point x0:
    P(x) := prod_k ZIP[w_k(x)] = sum_v c_v G_0(v; x),   c_v in the Z-module spanned by constant monomials
    {ZIP[w(x0)] (t-only letters after x -> x0)} x {G_0(v; x0)}.
Recursion on the weight via the total differential of ZIP words (Goncharov; verified numerically in selftest):
    d ZIP[a_1..a_n] = sum_i ZIP[..^a_i..] ( dlog(a_{i-1} - a_i) [i>=2]  -  dlog(a_{i+1} - a_i) [a_{n+1} := 0] ),
    dlog(c prod_l L_l^{e_l}) = sum_l e_l sum_{roots r of L_l} dx/(x - r)   (x-free letters drop),
so dP/dx = sum_r 1/(x-r) Q_r(x) and P = C + sum_r sum_u fib(Q_r)[u] G_0(r,u; x), C = P(x0) - sum (...) G_0(r,u; x0).
No arithmetic on algebraic points is needed here (integer coefficients); differences of letters are factored over the GI alphabet (auto-extended)."""
import os, sys, itertools
import slice_prov
from fractions import Fraction as Fr
from fibr import ZERO_HL
import k3field
from k3field import GT, Points

class Fibrator:
    def __init__(self, E, jx, jt, jj, x0):
        """E: fibr_gi.Engine over the group's GI alphabet; jx/jt/jj: indices of the last variable, t, j in the ctx; x0: Fraction base point"""
        self.E = E; self.A = E.A; self.jx, self.jt, self.jj = jx, jt, jj; self.x0 = Fr(x0)
        self.P = Points({}); self._memo = {(): {(): {(): Fr(1)}}}; self._dz = {}; self._pts = {}
        self.MINUS1 = E.hl(Fr(-1), {})
    # ---- points of a letter (lazy registration)
    def xdeg(self, i): return int(self.A.letters[i].degrees()[self.jx])
    def letter_points(self, i):
        if i not in self._pts:
            if i not in self.P.L:
                cs = k3field.xcoeffs_GT(self.A.letters[i], self.jx, self.jt, self.jj); self.P.L[i] = k3field.ptrim(cs); self.P.deg[i] = len(self.P.L[i]) - 1
            self._pts[i] = self.P.points_of_letter(i)
        return self._pts[i]
    def dlog_terms(self, c, Emono):
        """d/dx log(c prod L^e) = sum e_l/(x - r) over x-dependent letters and their roots: list of (Fraction e, point)"""
        out = []
        for l, e in Emono.items():
            if e == 0 or self.xdeg(l) == 0: continue
            for pt in self.letter_points(l): out.append((Fr(e), pt))
        return out
    # ---- total differential of one ZIP word: list of (subword, coefficient Fr, point)
    def dzip(self, word):
        if word in self._dz: return self._dz[word]
        n = len(word); out = []
        for i in range(n):
            sub = word[:i] + word[i + 1:]
            a = word[i]
            if i >= 1:
                prev = word[i - 1]
                if prev != a:
                    c, Em = self.E.hl_diff_mono(prev, a)
                    if c != 0:
                        for e, pt in self.dlog_terms(c, Em): out.append((sub, e, pt))
            nxt = word[i + 1] if i + 1 < n else ZERO_HL
            if nxt != a:
                if nxt == ZERO_HL:
                    c, Em = self.E.hl_c(a), self.E.hl_E(a)
                else:
                    c, Em = self.E.hl_diff_mono(nxt, a)
                if c != 0:
                    for e, pt in self.dlog_terms(c, Em): out.append((sub, -e, pt))
        self._dz[word] = out
        return out
    # ---- fibration of a product of ZIP words (sorted tuple of words)
    @staticmethod
    def cmul(c1, c2, scale=Fr(1)):
        """product of two constant-coefficient dicts {constkey: Fr}"""
        out = {}
        for k1, v1 in c1.items():
            for k2, v2 in c2.items():
                k = tuple(sorted(k1 + k2, key=repr)); out[k] = out.get(k, Fr(0)) + v1 * v2 * scale
        return {k: v for k, v in out.items() if v != 0}
    @staticmethod
    def cadd(acc, c, scale=Fr(1)):
        for k, v in c.items():
            w = acc.get(k, Fr(0)) + v * scale
            if w == 0: acc.pop(k, None)
            else: acc[k] = w
        return acc
    def fib(self, prod):
        prod = tuple(sorted(w for w in prod if w))
        if prod in self._memo: return self._memo[prod]
        if any(w == (self.MINUS1,) for w in prod):        # ZIP[-1] = 0
            self._memo[prod] = {}; return {}
        res = {}
        for wi, w in enumerate(prod):
            rest = prod[:wi] + prod[wi + 1:]
            for sub, e, pt in self.dzip(w):
                Q = rest + ((sub,) if sub else ())
                fq = self.fib(Q)
                for u, cu in fq.items():
                    v = (pt,) + u; acc = res.setdefault(v, {}); self.cadd(acc, cu, e)
                    if not acc: del res[v]
        # constant: P(x0) - sum_v c_v G_0(v; x0)
        const = {tuple(sorted((("Z0", w) for w in prod), key=repr)): Fr(1)}
        for v, cv in res.items():
            g = (("G0", v),)
            for k, c in cv.items():
                kk = tuple(sorted(k + g, key=repr)); const[kk] = const.get(kk, Fr(0)) - c
        const = {k: c for k, c in const.items() if c != 0}
        if const: res[()] = const
        self._memo[prod] = res
        return res

# ---------------------------------------------------------------- numeric side (for verification and final evaluation)
class NumCtx:
    """numeric context at rational t: generator values (roots of letters, slot -> k-th root in mp.polyroots order), point values, hyperlogs via hpath"""
    def __init__(self, fibr_or_points, letters_xcoeffs_fn, tval, dps, x0=None):
        import mpmath as mp, hpath
        self.mp = mp; mp.mp.dps = dps; self.dps = dps; self.tv = Fr(tval); self.H = hpath.HPath(dps); self.x0 = x0
        self.P = fibr_or_points; self._roots = {}; self.tol = mp.mpf(10) ** (-(dps - 8))
    def roots(self, i):
        if i not in self._roots:
            mp = self.mp; cs = [c.mp(self.tv, mp) for c in self.P.L[i]]
            with mp.workdps(self.dps + 40):
                rr = mp.polyroots(list(reversed(cs)), maxsteps=400, extraprec=4 * self.dps)
            self._roots[i] = [mp.mpc(r) for r in rr]
        return self._roots[i]
    def genval(self):
        class GV(dict):
            def __init__(s2, ctx): super().__init__(); s2.ctx = ctx
            def __missing__(s2, g):
                v = s2.ctx.roots(g[0])[g[1]]; s2[g] = v; return v
        return GV(self)
    def ptval(self, p):
        """(value, below flag): exact positive rationals and numerically-real positive values get the -i0 prescription"""
        mp = self.mp
        if p[0] == "Q":
            g = p[1]
            if g.is_real():
                re, im = g(self.tv)
                return (re, bool(re > 0))                       # exact Fraction letter (keeps ZIP[-1] = 0 and endpoint letters exact)
            re, im = g(self.tv); v = mp.mpc(mp.mpf(re.numerator) / re.denominator, mp.mpf(im.numerator) / im.denominator)
            if abs(v.imag) < self.tol * (1 + abs(v.real)) and v.real > 0: return (mp.mpc(v.real, 0), True)
            return (v, False)
        v = self.roots(p[1])[p[2]]
        if abs(v.imag) < self.tol * (1 + abs(v.real)):
            return (mp.mpc(v.real, 0), bool(v.real > 0))
        return (v, False)
    def G0(self, v, x):
        """G_0(v; x) for real x > 0 (iterated integral from 0, G(0;x) = log x): trailing zeros -> powers of log x (shuffle), then G(u/x; 1)"""
        import hpath, math
        mp = self.mp; xx = mp.mpf(Fr(x).numerator) / Fr(x).denominator
        ZER = ("E", 0)
        toks = tuple(ZER if (p[0] == "Q" and p[1].is_zero()) else ("P", p) for p in v)
        tot = mp.mpc(0); lx = mp.log(xx)
        for (k, u), c in hpath.decomp_trailing(toks, ZER).items():
            lets = []
            for tk in u:
                if tk == ZER: lets.append((Fr(0), False)); continue
                val, below = self.ptval(tk[1]); lets.append(((val / Fr(x)) if isinstance(val, Fr) else val / xx, below))
            g = self.H.G1(lets) if lets else 1
            tot += (mp.mpf(c.numerator) / c.denominator) * lx ** k / math.factorial(k) * g
        return tot
    def Z(self, v):
        """ZIP_reg over points: reg G_0(v; oo)"""
        lets = []
        for p in v:
            if p[0] == "Q" and p[1].is_zero(): lets.append((Fr(0), False)); continue
            val, below = self.ptval(p); lets.append((val, below))
        return self.H.zip_value(lets)

def selftest(tag="q_qbqgq_s1234_p34_g2", tval="1/2", x0="11/2", nprod=12, xs=("3", "3/5", "9")):
    """fibrate word products of a real stage-2 object and check  prod ZIP[w(x)] == sum_v c_v G_0(v;x)  numerically at several x"""
    import time, random, mpmath as mp
    import analyze_s2, slice_eval
    C, AL, E, ts, extra = analyze_s2.load_engine(tag)
    names = list(C.names()); last = "x%d" % extra["order"][2]; jx, jt, jj = names.index(last), names.index("t"), names.index("j")
    F = Fibrator(E, jx, jt, jj, Fr(x0)); dps = 50; mp.mp.dps = dps
    prods = sorted(set(tuple(sorted(Wk)) for (Dk, Wk) in ts if Wk), key=repr)
    random.Random(3).shuffle(prods); prods = prods[:nprod]
    S = slice_eval.Specialized(C, AL.letters, ts, Fr(tval), last)      # for HL letter values at x
    N = NumCtx(F.P, None, Fr(tval), dps)
    import hpath; H = hpath.HPath(dps)
    def zip_words_at(prod, x):
        tot = mp.mpc(1); xx = mp.mpf(Fr(x).numerator) / Fr(x).denominator
        for w in prod:
            lets = []
            for h in w:
                lv = S.hl_value(h, xx)
                if isinstance(lv, Fr): lets.append((lv, bool(lv > 0)))
                elif abs(lv.imag) < N.tol * (1 + abs(lv.real)): lets.append((mp.mpc(lv.real, 0), bool(lv.real > 0)))
                else: lets.append((lv, False))
            tot *= H.zip_value(lets)
        return tot
    def const_val(ck):
        v = mp.mpc(1)
        for atom in ck:
            if atom[0] == "Z0": v *= zip_words_at((atom[1],), F.x0)
            else: v *= N.G0(atom[1], F.x0)
        return v
    out = {}; t0 = time.time()
    for pi_, prod in enumerate(prods):
        S.letters = [S.letters[i] if i < len(S.letters) else None for i in range(len(S.letters))]
        fp = F.fib(prod)
        # the alphabet may have grown (difference letters): refresh the specialized letters
        S2 = slice_eval.Specialized(C, AL.letters, {}, Fr(tval), last); S.letters = S2.letters
        for x in xs:
            direct = zip_words_at(prod, Fr(x)); tot = mp.mpc(0)
            for v, cv in fp.items():
                cnum = sum((mp.mpf(c.numerator) / c.denominator) * const_val(ck) for ck, c in cv.items())
                tot += cnum * (N.G0(v, Fr(x)) if v else 1)
            out["prod%d[w%s]_x=%s" % (pi_, "+".join(str(len(w)) for w in prod), x)] = float(abs(tot - direct) / (1 + abs(direct)))
    out["_n_points_registered"] = sum(len(v) for v in F._pts.values()); out["_wall_s"] = round(time.time() - t0, 1); out["_memo"] = len(F._memo)
    return out
if __name__ == "__main__":
    import json, argparse
    ap = argparse.ArgumentParser(); ap.add_argument("--tag", default="q_qbqgq_s1234_p34_g2"); ap.add_argument("--t", default="1/2"); ap.add_argument("--x0", default="11/2"); ap.add_argument("--n", type=int, default=12)
    a = ap.parse_args(); r = selftest(a.tag, a.t, a.x0, a.n); print(json.dumps(r, indent=1)); print("PASS", all(v < 1e-40 for k, v in r.items() if not k.startswith("_")))
