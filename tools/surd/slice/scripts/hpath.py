"""hpath — numeric hyperlogarithms G(b_1..b_n; 1) and ZIP_reg words whose letters may lie ON the integration path with a side prescription.
Letters are (value, below) pairs: value = exact Fraction or mpmath number; below=True means the letter is 'value - i0', i.e. the integration
path passes ABOVE it (the universal rule of this package: every spurious on-contour pole is passed above, = contour deformed into Im > 0).
Method: if no letter is on the open path (0,1), delegate to gate/scripts/hlog_eval.GEvalAlg (Hoelder bisection + series, shuffle regularization of
leading 1 / trailing 0 via zip_num).  Otherwise deform the path 0 -> m -> 1, m = 1/2 + i c with c > 0 below every upper-half-plane letter near the
path, and use path concatenation  I_{0->m->1}(a_1..a_n) = sum_k I_{m->1}(a_1..a_k) I_{0->m}(a_{k+1}..a_n)  with the segment integrals as ordinary
G(.;1):  I_{0->m}(w) = G(w/m; 1),  I_{m->1}(w) = G((w-m)/(1-m); 1)  (regularization at the true endpoints 0, 1 is inherited).
ZIP_reg[a_1..a_n] = sum over letters of {1 (coef -1), a/(1+a) (coef +1, omitted if a == -1)} of G(...;1)  (zip_num semantics); a - i0 -> a/(1+a) - i0."""
import sys
import slice_prov
from fractions import Fraction as Fr
import mpmath as mp
import hlog_eval
from fibrate_dep import zip_num

def _n(a): return mp.mpf(a.numerator) / a.denominator if isinstance(a, Fr) else a
def _isreal_onpath(v, lo=0, hi=1):
    if isinstance(v, Fr): return lo < v and (hi is None or v < hi)
    z = mp.mpc(v); return abs(z.imag) < mp.mpf(10) ** (-(mp.mp.dps - 8)) and lo < z.real and (hi is None or z.real < hi)

def _ins_after_first(u, x):
    """all words obtained from u by inserting letter x at some position >= 1 (after u's first letter), with multiplicity"""
    out = {}
    for i in range(1, len(u) + 1):
        w = u[:i] + (x,) + u[i:]; out[w] = out.get(w, 0) + 1
    return out
_DL = {}
def decomp_leading(w, one):
    """w -> {(k, v): c} with  w = sum c * (one^k shuffle v) as iterated-integral words, v not starting with `one`  (so G(w) = sum c G(one)^k/k! G(v))."""
    key = (w, one)
    if key in _DL: return _DL[key]
    r = 0
    while r < len(w) and w[r] == one: r += 1
    if r == 0: res = {(0, w): Fr(1)}
    elif r == len(w): res = {(r, ()): Fr(1)}          # one^r = (one)^{shuffle r}/r!  -> G = G(one)^r/r!
    else:
        u = w[r:]; res = {}
        # (one) sh (one^{r-1} u) = r * one^r u + sum_{ins after first of u} one^{r-1} u'
        for (k, v), c in decomp_leading(w[1:], one).items():       # one^{r-1} u, then shuffled with one more `one`
            res[(k + 1, v)] = res.get((k + 1, v), Fr(0)) + c * (k + 1) / r    # (one) sh (one^k sh v) = (k+1) one^{k+1} sh v  (as normalized words)
        for u2, mlt in _ins_after_first(u, one).items():
            for (k, v), c in decomp_leading(w[1:r] + u2, one).items():
                res[(k, v)] = res.get((k, v), Fr(0)) - c * mlt / r
        res = {kv: c for kv, c in res.items() if c != 0}
    _DL[key] = res
    return res
def decomp_trailing(w, zero):
    """mirror: w = sum c * (v shuffle zero^k), v not ending with `zero`"""
    rev = tuple(reversed(w)); d = decomp_leading(rev, zero)
    return {(k, tuple(reversed(v))): c for (k, v), c in d.items()}

class HPath:
    def __init__(self, dps):
        self.dps = dps; self.ge = hlog_eval.GEvalAlg(dps); self.Z = zip_num.ZipEval(); self.cache = {}; self.stats = {"deformed": 0, "direct": 0}
    def key(self, letters):
        return tuple((("F", a.numerator, a.denominator) if isinstance(a, Fr) else ("N", mp.nstr(mp.mpc(a), self.dps + 5)), self._side(b)) for a, b in letters)
    # ---- G(w;1) with regularization (leading exact 1, trailing exact 0) for letters given as plain values (no on-path letters)
    def G_reg_plain(self, w):
        tot = mp.mpc(0)
        for v, c0 in self.Z.strip_trailing_zeros(tuple(w)).items():
            for u, c in self.Z.reg_leading_ones(v).items():
                q = c0 * c
                tot += (mp.mpf(q.numerator) / q.denominator) * self.ge.G1(tuple(u))
        return tot
    @staticmethod
    def _side(b):
        if b is True or b == "below": return "below"
        if b == "above": return "above"
        return None
    def G1(self, letters):
        """G(b_1..b_n; 1), letters = [(value, side)], side in {False/None, True/'below', 'above'}: 'below' = letter is value - i0 (path passes above it),
        'above' = letter is value + i0 (path passes below).  Exact Fractions 0 and 1 allowed at the ends (shuffle-regularized in the original variable)."""
        letters = [(a, self._side(b)) for a, b in letters]; k = self.key(letters)
        if k in self.cache: return self.cache[k]
        onpath = [(i, a, sd) for i, (a, sd) in enumerate(letters) if _isreal_onpath(a) and not (isinstance(a, Fr) and a in (0, 1))]
        if not onpath:
            self.stats["direct"] += 1
            val = self.G_reg_plain([a for a, b in letters])
        else:
            self.stats["deformed"] += 1
            for i, a, sd in onpath: assert sd in ("below", "above"), ("letter on the path (0,1) without a side prescription", a)
            # detour height: below every off-axis letter near the path
            ims = []
            for a, b in letters:
                if isinstance(a, Fr): continue
                z = mp.mpc(a)
                if abs(z.imag) > mp.mpf(10) ** (-(mp.mp.dps - 8)) and -0.25 < z.real < 1.25: ims.append(abs(z.imag))
            h = mp.mpf(1) / 4
            if ims: h = min(h, min(ims) / 3)
            assert h > mp.mpf(10) ** (-12), ("off-axis letter too close to the path for the polygonal detour", h)
            # waypoints: above (side below) or below (side above) each on-path letter, ordered by real part; merge equal positions
            pts = sorted(set((mp.mpf(mp.re(_n(a))), sd) for i, a, sd in onpath), key=lambda x: x[0])
            way = [mp.mpc(0)] + [mp.mpc(re, h if sd == "below" else -h) for re, sd in pts] + [mp.mpc(1)]
            nseg = len(way) - 1; n = len(letters); val = mp.mpc(0)
            ONE = ("E", 1); ZER = ("E", 0)
            def tok(a):
                if isinstance(a, Fr) and a == 1: return ONE
                if isinstance(a, Fr) and a == 0: return ZER
                return ("V", mp.mpc(_n(a)))
            toks = tuple(tok(a) for a, b in letters)
            import math, itertools as _it
            m_last = way[-2]; m_first = way[1]
            Lone = -mp.log(1 - m_last)          # G(1; segment m_last->1) regularized in the ORIGINAL variable at the endpoint 1
            Lzero = mp.log(m_first)             # G(0; segment 0->m_first) regularized in the original variable at 0
            def seg_val(piece, j):
                """iterated integral of the sub-word `piece` along segment j (way[j] -> way[j+1])"""
                p, q = way[j], way[j + 1]; d = q - p
                if not piece: return mp.mpc(1)
                if j == nseg - 1:            # last segment: ends at 1, leading exact 1's regularized
                    tot = mp.mpc(0)
                    for (k1, v), c1 in decomp_leading(piece, ONE).items():
                        lv = [(Fr(1) if tk == ONE else (((mp.mpc(0) if tk == ZER else tk[1]) - p) / d)) for tk in v]
                        tot += (mp.mpf(c1.numerator) / c1.denominator) * Lone ** k1 / math.factorial(k1) * (self.ge.G1(tuple(lv)) if lv else 1)
                    return tot
                if j == 0:                   # first segment: starts at 0, trailing exact 0's regularized
                    tot = mp.mpc(0)
                    for (k0, v), c0 in decomp_trailing(piece, ZER).items():
                        rv = [(Fr(0) if tk == ZER else ((mp.mpc(1) if tk == ONE else tk[1]) / d)) for tk in v]
                        tot += (mp.mpf(c0.numerator) / c0.denominator) * Lzero ** k0 / math.factorial(k0) * (self.ge.G1(tuple(rv)) if rv else 1)
                    return tot
                lv = [(((mp.mpc(1) if tk == ONE else (mp.mpc(0) if tk == ZER else tk[1])) - p) / d) for tk in piece]
                return self.ge.G1(tuple(lv))
            if nseg == 1:
                val = self.G_reg_plain([a for a, b in letters])
            else:
                # Chen: I_{gamma_1...gamma_nseg}(a_1..a_n) = sum over cut positions 0<=c_1<=...<=c_{nseg-1}<=n of prod_j I_{gamma_{nseg-j}}(piece_j),
                # piece_0 = a_1..a_{c_1} on the LAST segment, ..., last piece on the FIRST segment
                for cuts in _it.combinations_with_replacement(range(n + 1), nseg - 1):
                    bounds = (0,) + cuts + (n,); prod = mp.mpc(1)
                    for jpiece in range(nseg):
                        piece = toks[bounds[jpiece]:bounds[jpiece + 1]]; seg = nseg - 1 - jpiece
                        prod *= seg_val(piece, seg)
                        if prod == 0: break
                    val += prod
        self.cache[k] = val
        return val
    def G1_forced(self, letters, c=0.2):
        """debug: always use the deformed path (for words with no on-path letter the result must equal the direct value)"""
        saved = _isreal_onpath
        letters = list(letters)
        # emulate by temporarily appending nothing: replicate the deformed branch
        globals()["_FORCE"] = True
        try:
            k = ("forced",) + self.key(letters)
            onp = [i for i, (a, b) in enumerate(letters)]
            # copy of the deformed branch with fixed c
            m = mp.mpc(mp.mpf(1) / 2, c); one_m = 1 - m; n = len(letters); val = mp.mpc(0)
            ONE = ("E", 1); ZER = ("E", 0)
            def tok(a):
                if isinstance(a, Fr) and a == 1: return ONE
                if isinstance(a, Fr) and a == 0: return ZER
                return ("V", mp.mpc(_n(a)))
            toks = tuple(tok(a) for a, b in letters); Lone = -mp.log(one_m); Lzero = mp.log(m)
            import math
            for kk in range(n + 1):
                left = toks[:kk]; right = toks[kk:]; vl = mp.mpc(0)
                for (k1, v), c1 in decomp_leading(left, ONE).items():
                    lv = [(Fr(1) if tk == ONE else (((mp.mpc(0) if tk == ZER else tk[1]) - m) / one_m)) for tk in v]
                    vl += (mp.mpf(c1.numerator) / c1.denominator) * Lone ** k1 / math.factorial(k1) * (self.ge.G1(tuple(lv)) if lv else 1)
                vr = mp.mpc(0)
                for (k0, v), c0 in decomp_trailing(right, ZER).items():
                    rv = [(Fr(0) if tk == ZER else ((mp.mpc(1) if tk == ONE else tk[1]) / m)) for tk in v]
                    vr += (mp.mpf(c0.numerator) / c0.denominator) * Lzero ** k0 / math.factorial(k0) * (self.ge.G1(tuple(rv)) if rv else 1)
                val += vl * vr
            return val
        finally:
            globals()["_FORCE"] = False
    def zip_value(self, letters):
        """ZIP_reg[a_1..a_n], letters = [(value, below)]; values must not be exactly on (0,oo) unless below=True."""
        letters = list(letters)
        exps = [((), Fr(1))]
        for a, b in letters:
            opts = [((Fr(1), None), Fr(-1))]
            if not (isinstance(a, Fr) and a == -1):
                if isinstance(a, Fr): mob = (a / (1 + a), b)
                else: mob = (mp.mpc(a) / (1 + mp.mpc(a)), b)
                if _isreal_onpath(a, 0, None) and not self._side(b):
                    raise AssertionError(("on-contour ZIP letter without prescription", a))
                opts.append((mob, Fr(1)))
            exps = [(ww + (l,), c * s) for (ww, c) in exps for (l, s) in opts]
        tot = mp.mpc(0)
        for ww, c in exps: tot += (mp.mpf(c.numerator) / c.denominator) * self.G1(list(ww))
        return tot

def selftest(dps=40):
    mp.mp.dps = dps; H = HPath(dps); Zr = hlog_eval.ZipEvalAlg(dps); out = {}
    # 1. off-path words: forced deformation vs direct
    import random; rng = random.Random(5)
    for trial in range(10):
        n = rng.choice([1, 2, 3]); w = []
        for pos in range(n):
            kind = rng.choice(["neg", "cplx", "zero", "one"])
            if kind == "neg": w.append(Fr(-rng.randint(1, 9), rng.randint(1, 7)))
            elif kind == "zero": w.append(Fr(0))
            elif kind == "one": w.append(Fr(1))
            else: w.append(mp.mpc(rng.uniform(-2, 2), rng.choice([-1, 1]) * rng.uniform(0.3, 2)))
        if all(isinstance(a, Fr) and a == 0 for a in w) or all(isinstance(a, Fr) and a == 1 for a in w): continue
        direct = H.G_reg_plain(w)
        val = H.G1_forced([(a, False) for a in w], c=rng.choice([0.2, 0.11, 0.35]))
        out["concat_vs_direct_%d_%s" % (trial, "".join("1" if (isinstance(a, Fr) and a == 1) else ("0" if (isinstance(a, Fr) and a == 0) else "x") for a in w))] = float(abs(val - direct) / (abs(direct) + 1))
    # 2. weight 1 on path: G(b - i0; 1) = log(1 - 1/(b - i0)) = log((1-b)/b) + i*pi  (b in (0,1))
    b = Fr(2, 5); g = H.G1([(b, True)]); ref = mp.log((1 - _n(b)) / _n(b)) - mp.mpc(0, 1) * mp.pi     # int_0^1 dt/(t - b + i0) = PV - i pi
    out["w1_onpath"] = float(abs(g - ref))
    # 3. ZIP[s - i0] = -log(-(s - i0)) = -(log s + i pi)
    s = Fr(3, 2); z = H.zip_value([(s, True)]); out["zip1_onpath"] = float(abs(z - (-(mp.log(_n(s)) + mp.mpc(0, 1) * mp.pi))))
    # 4. weight 2 on path vs nested quadrature with letters displaced by -i*eta (eta = 1e-12, quad at dps 24)
    eta = mp.mpf("1e-12")
    for w in ([(Fr(2, 5), True), (Fr(-1, 3), False)], [(Fr(-2), False), (Fr(3, 7), True)], [(Fr(1, 3), True), (Fr(2, 3), True)], [(mp.mpc(0.3, 0.5), False), (Fr(1, 2), True)]):
        g = H.G1(w)
        wl = [(_n(a) - mp.mpc(0, 1) * eta) if b else mp.mpc(_n(a)) for a, b in w]
        with mp.workdps(24):
            a1, a2 = wl
            inner = lambda t1: mp.quad(lambda t2: 1 / (t2 - a2), [0, mp.re(a2), t1] if 0 < mp.re(a2) < 1 else [0, t1])
            pts = [0] + sorted(set([mp.re(x) for x in wl if 0 < mp.re(x) < 1])) + [1]
            q = mp.quad(lambda t1: inner(t1) / (t1 - a1), pts)
        out["w2_onpath_vs_eta_quad_%d" % len(out)] = float(abs(g - q))
    # 5. weight 3 on path: detour-height independence (c = 1/4 default vs forced smaller height via a dummy far letter is not available) ->
    #    compare against letters displaced by -i*eta evaluated by the plain evaluator (eta = 1e-6; agreement O(eta log^2 eta))
    for w in ([(Fr(2, 5), True), (Fr(-1, 3), False), (Fr(1, 4), True)], [(mp.mpc(0.3, 0.5), False), (Fr(1, 2), True), (Fr(0), False)], [(Fr(1), False), (Fr(1, 2), True), (Fr(-3), False)]):
        g = H.G1(w); eta6 = mp.mpf("1e-6")
        wl = [((_n(a) - mp.mpc(0, 1) * eta6) if b else (a if isinstance(a, Fr) else mp.mpc(a))) for a, b in w]
        q = H.G_reg_plain(wl)
        out["w3_onpath_vs_eta1e-6_plain_%d" % len(out)] = float(abs(g - q))
    return out
if __name__ == "__main__":
    import json; r = selftest(); print(json.dumps(r, indent=1)); print("PASS", all((v < 1e-30 if ("eta" not in k) else (v < 1e-9 if "w2" in k else v < 1e-4)) for k, v in r.items()))
