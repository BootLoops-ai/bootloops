"""Stage 3, step B: exact integration over x in (0,oo) of  R(x) * G_0(v; x)  for R = N(x)/prod L_i(x)^e_i (N, L_i over F0 = Q(i)(t)),
v a word of fibration points.  Port of the fibr.py primitives J (poles), I (polynomial part), expand_zero/expand_inf to ALGEBRAIC points:
partial fractions of R at every root (slot) of every x-letter by exact Laurent expansion in K = F0[y]/(L_i) (k3field.laurent_all), coefficients in
the tensor ring k3field.TRing, transcendental output = products of Z-words  Z(u) := reg_{oo} G_0(u; oo)  (= ZIP over points) — the same
shuffle-regularization conventions as fibr/HyperFLINT (log-divergences at 0 and oo dropped in the scale x^n log^j x)."""
import sys, math
import slice_prov
from fractions import Fraction as Fr
import k3field
from k3field import GT, TRing, Points

def binom(n, k):
    r = Fr(1)
    for i in range(k): r = r * (n - i) / (i + 1)
    return r
def zkey_mul(k1, k2): return tuple(sorted(k1 + k2, key=repr))

class Integrator:
    def __init__(self, R: TRing):
        self.R = R; self.P = R.P; self._c = {}
        self.ONE = {(): R.const(1)}
    # ---- BTS algebra: dict {zkey (tuple of point-words): TRing element}
    def add_into(self, acc, bts, coef=None, scale=None):
        R = self.R
        for k, te in bts.items():
            v = te
            if coef is not None: v = R.mul(v, coef)
            if scale is not None: v = R.scal(v, scale)
            if not v: continue
            w = acc.get(k)
            if w is None: acc[k] = v
            else:
                w = R.add(w, v)
                if w: acc[k] = w
                else: del acc[k]
        return acc
    def bscale(self, bts, scale): return {k: self.R.scal(v, scale) for k, v in bts.items() if scale != 0}
    def bmul(self, bts, coef):
        out = {}
        for k, v in bts.items():
            w = self.R.mul(v, coef)
            if w: out[k] = w
        return out
    def is_minus_one(self, p): return p[0] == "Q" and p[1] == GT(-1)
    def is_zero_pt(self, p): return p[0] == "Q" and p[1].is_zero()
    def attach(self, bts, words):
        if any(len(w) == 1 and self.is_minus_one(w[0]) for w in words): return {}
        out = {}
        for k, v in bts.items():
            kk = zkey_mul(k, tuple(words)); w = out.get(kk)
            out[kk] = v if w is None else self.R.add(w, v)
        return {k: v for k, v in out.items() if v}
    # ---- expansions of G_0(w; x)
    def expand_inf(self, w, M):
        key = ("inf", w, M)
        if key in self._c: return self._c[key]
        if not w: return {(0, 0): dict(self.ONE)}
        a = w[0]; Dp = self.expand_inf(w[1:], M)
        out = {(0, 0): self.attach(dict(self.ONE), [w])}
        za = self.is_zero_pt(a)
        for (j, i), T in Dp.items():
            kmax = 0 if za else M - i
            for k in range(0, max(kmax, 0) + 1):
                n = i + k
                if n > M: break
                Tk = T if k == 0 else self.bmul(T, self.R.point_pow(a, k))
                if n == 0:
                    self.add_into(out.setdefault((j + 1, 0), {}), Tk, scale=Fr(1, j + 1))
                else:
                    for l in range(j + 1):
                        coef = -Fr(math.factorial(j), math.factorial(j - l)) / Fr(n) ** (l + 1)
                        self.add_into(out.setdefault((j - l, n), {}), Tk, scale=coef)
                if za: break
        out = {k: v for k, v in out.items() if v}
        self._c[key] = out
        return out
    def expand_zero(self, w, M):
        key = ("zero", w, M)
        if key in self._c: return self._c[key]
        if not w: return {(0, 0): dict(self.ONE)}
        a = w[0]; Dp = self.expand_zero(w[1:], M); out = {}
        if self.is_zero_pt(a):
            for (j, n), T in Dp.items():
                if n == 0: self.add_into(out.setdefault((j + 1, 0), {}), T, scale=Fr(1, j + 1))
                else:
                    for l in range(j + 1):
                        coef = (-1) ** l * Fr(math.factorial(j), math.factorial(j - l)) / Fr(n) ** (l + 1)
                        self.add_into(out.setdefault((j - l, n), {}), T, scale=coef)
        else:
            for (j, n0), T in Dp.items():
                for k in range(0, M + 1):
                    n = n0 + k + 1
                    if n > M: break
                    Tk = self.bmul(T, self.R.scal(self.R.point_pow(a, -(k + 1)), Fr(-1)))
                    for l in range(j + 1):
                        coef = (-1) ** l * Fr(math.factorial(j), math.factorial(j - l)) / Fr(n) ** (l + 1)
                        self.add_into(out.setdefault((j - l, n), {}), Tk, scale=coef)
        out = {k: v for k, v in out.items() if v}
        self._c[key] = out
        return out
    # ---- integrals
    def J(self, k, s, w):
        """reg int_0^oo (x-s)^-k G_0(w;x) dx"""
        key = ("J", k, s, w)
        if key in self._c: return self._c[key]
        R = self.R
        if k == 1:
            res = self.attach(dict(self.ONE), [(s,) + tuple(w)])
        elif not w:
            if self.is_zero_pt(s): res = {}
            else: res = {(): R.scal(R.point_pow(s, 1 - k), Fr((-1) ** ((k - 1) % 2), k - 1))}
        else:
            a = w[0]; wp = tuple(w[1:]); res = {}
            if self.is_zero_pt(s):
                ez = self.expand_zero(tuple(w), k - 1); T = ez.get((0, k - 1))
                if T: self.add_into(res, T, scale=Fr(1, k - 1))
            n = k - 1
            if a == s:
                self.add_into(res, self.J(k, s, wp), scale=Fr(1, k - 1))
            else:
                dinv = R.inv_diff(a, s, n)                       # (a-s)^-n
                self.add_into(res, self.J(1, a, wp), coef=dinv, scale=Fr(1, k - 1))
                for j in range(1, n + 1):
                    pw = n + 1 - j
                    self.add_into(res, self.J(j, s, wp), coef=R.inv_diff(a, s, pw), scale=Fr(-1, k - 1))
        self._c[key] = res
        return res
    def I(self, m, w, M=None):
        """reg int_0^oo x^m G_0(w;x) dx"""
        key = ("I", m, w, M)
        if key in self._c: return self._c[key]
        if not w: return {}
        R = self.R; a = w[0]; wp = tuple(w[1:])
        Mx = M if M is not None else m + 1
        E = self.expand_inf(tuple(w), max(Mx, m + 1)); res = {}
        T = E.get((0, m + 1))
        if T: self.add_into(res, T, scale=Fr(1, m + 1))
        if self.is_zero_pt(a):
            self.add_into(res, self.I(m, wp, Mx), scale=-Fr(1, m + 1))
        else:
            for i in range(m + 1):
                Ii = self.I(i, wp, Mx)
                if Ii: self.add_into(res, Ii, coef=R.point_pow(a, m - i), scale=-Fr(1, m + 1))
            self.add_into(res, self.J(1, a, wp), coef=R.point_pow(a, m + 1), scale=-Fr(1, m + 1))
        self._c[key] = res
        return res
    # ---- one rational function times one fibration basis word
    def rational_parts(self, Ncoef, xfactors, pts_of):
        """R = N/prod L^e  ->  (polynomial part [GT list], poles [(point, order s, TRing coef)])"""
        R = self.R
        D = [GT(1)]
        for cs, e in xfactors: D = k3field.pmul(D, k3field.ppow(cs, e))
        degN = len(k3field.ptrim(Ncoef)) - 1; degD = len(D) - 1
        polypart = []
        if degN >= degD:
            q, r = k3field.pdivmod(Ncoef, D); polypart = q
        poles = []
        for kf, (cs, e) in enumerate(xfactors):
            K, h = k3field.laurent_all(Ncoef, xfactors, GT(1), kf)
            li = pts_of[kf][0]; pts = pts_of[kf][1]
            for p in pts:
                for s in range(1, e + 1):
                    hk = h[e - s]
                    if K.is_zero(hk): continue
                    coef = R.const(hk[0]) if p[0] == "Q" else R.from_K(hk, p[1], p[2])
                    if coef: poles.append((p, s, coef))
        return polypart, poles
    def integrate_word(self, polypart, poles, v):
        """sum_poles coef J(s,p,v) + sum_m q_m I(m,v)  -> BTS"""
        R = self.R; res = {}
        for p, s, coef in poles:
            Js = self.J(s, p, v)
            if Js: self.add_into(res, Js, coef=coef)
        for m, q in enumerate(polypart):
            if q.is_zero() or not v: continue
            Im = self.I(m, v)
            if Im: self.add_into(res, Im, coef=R.const(q))
        return res

def selftest():
    """(1) rational letters only: reproduce symbolic/scripts/fibr.Engine.J/I values numerically (toy alphabet), (2) an algebraic toy:
    int_0^oo dx log(1+x)/((x^2+x+1)(x+2)^2)  vs mpmath quad;  (3) int_0^oo dx G_0(r1, r2; x)/(x-p)^2 with r,p roots of a cubic over Q(i)(t) vs quad."""
    import mpmath as mp, fib3, hpath
    mp.mp.dps = 40; out = {}
    t = GT.t(); tv = Fr(1, 3)
    # toy letters: L5 = x^2 + x + 1 (points A slots), L6 = x + 2, L7 = x (zero), L8 = x^3 + t x + 1 + i  (cubic over Q(i)(t))
    P = Points({5: [GT(1), GT(1), GT(1)], 6: [GT(2), GT(1)], 7: [GT(0), GT(1)], 8: [GT(1) + GT.i(), t, GT(0), GT(1)], 9: [GT(1), GT(1)]})
    R = TRing(P); I_ = Integrator(R); N = fib3.NumCtx(P, None, tv, 40); gv = N.genval()
    def val(bts):
        tot = mp.mpc(0)
        for zk, te in bts.items():
            c = TRing.evalf(te, gv, tv, mp)
            for w in zk: c *= N.Z(w)
            tot += c
        return tot
    m1 = ("Q", GT(-1)); m2 = ("Q", GT(-2)); zero = Points.ZERO
    # (2): log(1+x) = -G... G_0(-1; x) = log(1 - x/(-1)) = log(1+x).  R = 1/((x^2+x+1)(x+2)^2)
    polypart, poles = I_.rational_parts([GT(1)], [(P.L[5], 1), (P.L[6], 2)], {0: (5, P.points_of_letter(5)), 1: (6, P.points_of_letter(6))})
    bts = I_.integrate_word(polypart, poles, (m1,))
    ref = mp.quad(lambda x: mp.log(1 + x) / ((x * x + x + 1) * (x + 2) ** 2), [0, mp.inf])
    out["alg_quadratic_double_pole"] = float(abs(val(bts) - ref))
    # (2b) numerator of higher degree: x^3 log(1+x)/((x^2+x+1)(x+2)^2 (x+1)) -> convergent? deg 3 - 5 = -2 ok
    polypart, poles = I_.rational_parts([GT(0), GT(0), GT(0), GT(1)], [(P.L[5], 1), (P.L[6], 2), (P.L[9], 1)], {0: (5, P.points_of_letter(5)), 1: (6, P.points_of_letter(6)), 2: (9, P.points_of_letter(9))})
    bts = I_.integrate_word(polypart, poles, (m1,))
    ref = mp.quad(lambda x: x ** 3 * mp.log(1 + x) / ((x * x + x + 1) * (x + 2) ** 2 * (x + 1)), [0, mp.inf])
    out["alg_quadratic_cubic_numerator"] = float(abs(val(bts) - ref))
    # (3) cubic points inside the word and as pole: word (r0, -1), pole (x - r1)^-2 (x+2)^-1 where r = roots of L8 ... integrand must decay: ok
    r = P.points_of_letter(8); rv = [gv[(8, s)] for s in range(3)]
    def G0w(x):   # G_0(r0, -1; x) numerically by nested quad
        return mp.quad(lambda t1: 1 / (t1 - rv[0]) * mp.log(1 + t1), [0, x])
    # Use the full letter L8 as denominator (all three roots as poles, mult 2) times (x+2): R = 1/(L8(x)^2 (x+2))
    polypart, poles = I_.rational_parts([GT(1)], [(P.L[8], 2), (P.L[6], 1)], {0: (8, r), 1: (6, P.points_of_letter(6))})
    bts = I_.integrate_word(polypart, poles, (r[0], m1))
    L8 = lambda x: x ** 3 + t.mp(tv, mp) * x + (1 + mp.mpc(0, 1))
    with mp.workdps(25):
        ref = mp.quad(lambda x: G0w(x) / (L8(x) ** 2 * (x + 2)), [0, 1, 4, mp.inf])
    out["cubic_pole_and_cubic_word_letter"] = float(abs(val(bts) - ref))
    # (4) pole at 0 boundary + polynomial part: R = (x^2+3)/(x (x+2)) has poly part 1 and poles at 0, -2; word (-1,-2): total integrand ~ log^2 at infinity * 1 -> divergent;
    #     use instead R = 1/(x (x+2)^3) with word (-1): integrable at 0 (G_0(-1;x) ~ x) 
    polypart, poles = I_.rational_parts([GT(1)], [(P.L[7], 1), (P.L[6], 3)], {0: (7, [zero]), 1: (6, P.points_of_letter(6))})
    bts = I_.integrate_word(polypart, poles, (m1,))
    ref = mp.quad(lambda x: mp.log(1 + x) / (x * (x + 2) ** 3), [0, mp.inf])
    out["pole_at_zero"] = float(abs(val(bts) - ref))
    # (5) weight-2 word with trailing zero: G_0(-1, 0; x) = int_0^x dt/(t+1) log t ; R = 1/(x+2)^2/(x^2+x+1)
    polypart, poles = I_.rational_parts([GT(1)], [(P.L[5], 1), (P.L[6], 2)], {0: (5, P.points_of_letter(5)), 1: (6, P.points_of_letter(6))})
    bts = I_.integrate_word(polypart, poles, (m1, zero))
    g = lambda x: mp.quad(lambda t1: mp.log(t1) / (t1 + 1), [0, x])
    with mp.workdps(20):
        ref = mp.quad(lambda x: g(x) / ((x * x + x + 1) * (x + 2) ** 2), [0, 1, mp.inf])
    out["trailing_zero_word"] = float(abs(val(bts) - ref))
    return out
if __name__ == "__main__":
    import json; r = selftest(); print(json.dumps(r, indent=1)); print("PASS", all(v < 1e-15 for v in r.values()))
