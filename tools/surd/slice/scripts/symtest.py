"""Symbol-level survival test of a square-root letter class in an assembled slice function (out/*.json.gz), per transcendental weight.
Field class C = square class of Delta(t) in Q(t)^*/squares; tau = the automorphism sqrt(Delta) -> -sqrt(Delta): swaps the two root slots of EVERY
quadratic X-letter whose discriminant is in C, fixes Q(i)(t), and fixes the roots of all other letters (checked: their discriminants / the cubic
discriminants are not in C).  Every symbol entry is a product of atoms; an atom gamma moved by tau contributes the ODD atom u_gamma = gamma/tau(gamma)
(log gamma = 1/2 log N(gamma) + 1/2 log u_gamma).  Relations among the odd atoms (only torsion can relate them to tau-fixed elements) are found by
LLL on 300-digit logs at a transcendental t* = e/7 and verified at t** = pi/9; U_free = odd atoms
modulo relations.  The weight-w symbol S_w = sum_terms c_T (e_1 (x) ... (x) e_w) is probed by the homomorphism-valued functionals
   Phi_p = sum_T c_T(t*) * oddvec(e_p) * prod_{q != p} log|e_q(t*)|   in  C (x) Q^{U_free},  p = 1..w,
which vanish identically on every tensor without odd content; Phi_p != 0 for some p  =>  sqrt(Delta) letters SURVIVE in the weight-w symbol;
all Phi_p = 0 at two transcendental points (and with log|.| replaced by a second homomorphism, the valuation-free log at t**) => absent
(up to an accidental zero of the functionals).  Coefficients c_T are the exact tensor-ring coefficients evaluated numerically at t*."""
import os, sys, json, gzip, time, argparse, itertools, collections, math
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import slice_prov, hpath
from fractions import Fraction as Fr
import mpmath as mp, sympy as sp
from flint import fmpz_mat
from k3field import GT
WORK = slice_prov.WORK; DATA = slice_prov.DATA; OUT = os.path.join(DATA, "out")
T = sp.Symbol("t")

# ---------------------------------------------------------------- numeric context at a real (transcendental) t
class Num:
    def __init__(self, doc, tstar, dps):
        self.doc = doc; self.t = tstar; self.dps = dps; self._roots = {}; self.x0 = mp.mpf(11) / 2
    def gt(self, ser):
        A, B, D = [[Fr(p, q) for p, q in lst] for lst in ser]
        def ev(cs):
            acc = mp.mpf(0)
            for c in reversed(cs): acc = acc * self.t + mp.mpf(c.numerator) / c.denominator
            return acc
        return mp.mpc(ev(A), ev(B)) / ev(D) if D else mp.mpc(ev(A), ev(B))
    ROOT_DPS = 340      # fixed working precision for X-coefficients and polyroots, independent of the caller's dps (keeps the root labeling identical across calls)
    def roots(self, L):
        """CANONICAL slot order: coefficients and roots at the fixed precision ROOT_DPS whatever mp.mp.dps is, roots sorted by
        (Re quantised at 1e-25, Im): a complex-conjugate pair is always (Im<0, Im>0), real roots ascend.  A subclass/hook may then
        re-order by continuation from a reference point (see ContNum)."""
        if L not in self._roots:
            with mp.workdps(max(self.ROOT_DPS, self.dps + 40)):
                cs = [self.gt(c) for c in self.doc["letters"][L]["x_coeffs(GT ser, increasing powers of X)"]]
                rr = mp.polyroots(list(reversed(cs)), maxsteps=600, extraprec=2 * max(self.ROOT_DPS, self.dps + 40))
                rr = [mp.mpc(r) for r in rr]
                rr.sort(key=lambda r: (int(mp.nint(mp.re(r) * mp.mpf(10) ** 25)), mp.im(r)))
            self._roots[L] = rr
        return self._roots[L]
    def point(self, p):
        return self.gt(p[1]) if p[0] == "Q" else self.roots(p[1])[p[2]]

    def hl(self, h, x0=None):
        """hyperlog letter of a Z0 constant: c * prod L(x0,t,j)^e with L given by polynomial keys"""
        x0 = self.x0 if x0 is None else mp.mpf(x0)
        (pn, pd), E = h; v = mp.mpc(mp.mpf(pn) / pd)
        for L, e in E:
            items = json.loads(L); val = mp.mpc(0)
            for (ex, et, ej), (pp, qq) in items: val += (mp.mpf(pp) / qq) * x0 ** ex * self.t ** et * mp.mpc(0, 1) ** ej
            v = v * val ** e
        return v

class ContNum(Num):
    """Num at t whose slot bijection is CARRIED from a reference Num (at t_ref) by root continuation along the real segment t_ref -> t
    (nstep sub-steps, nearest-neighbor matching in the chordal metric), so 'slot s of letter L' denotes the same analytic root branch at both
    points.  Records which letters were re-ordered relative to the canonical order."""
    def __init__(self, doc, tval, dps, ref, nstep=96):
        super().__init__(doc, tval, dps); self.ref = ref; self.nstep = nstep; self.reordered = []
    @staticmethod
    def _chord(a, b): return abs(a - b) / mp.sqrt((1 + abs(a) ** 2) * (1 + abs(b) ** 2))
    def roots(self, L):
        if L not in self._roots:
            canon = Num.roots(self, L)            # caches canonical order in self._roots[L]
            ref = self.ref.roots(L)
            if len(ref) > 1:
                prev = [mp.mpc(r) for r in ref]
                with mp.workdps(60):
                    for k in range(1, self.nstep + 1):
                        tk = self.ref.t + (self.t - self.ref.t) * k / self.nstep
                        cur = Num(self.doc, tk, 40).roots(L); left = list(range(len(cur))); order = []
                        for r0 in prev:
                            j = min(left, key=lambda jj: ContNum._chord(cur[jj], r0)); left.remove(j); order.append(j)
                        prev = [cur[j] for j in order]
                left = list(range(len(canon))); order = []
                for r0 in prev:
                    j = min(left, key=lambda jj: ContNum._chord(canon[jj], r0)); left.remove(j); order.append(j)
                if order != list(range(len(canon))) and L not in self.reordered: self.reordered.append(L)
                self._roots[L] = [canon[j] for j in order]
        return self._roots[L]

POINTS = {"e/7": lambda: mp.e / 7, "pi/9": lambda: mp.pi / 9, "1/sqrt7": lambda: 1 / mp.sqrt(7), "sqrt2/3": lambda: mp.sqrt(2) / 3, "e/8": lambda: mp.e / 8, "ln2/2": lambda: mp.log(2) / 2}
def point_value(expr):
    """FREE-FORM generic point: a POINTS name, or any expression in integers, decimals, + - * / ( ), sqrt(.), 'sqrtN'
    shorthand, pi, e, E, ln(.)/log(.), e.g. 'sqrt(3)/5', '2/(3+sqrt5)', '0.3771', '3/7'; evaluated at the CURRENT mp precision (call inside workdps)."""
    if expr in POINTS: return POINTS[expr]()
    import re as _re
    src = _re.sub(r"sqrt(\d+)", r"sqrt(\1)", expr.strip())
    if not _re.fullmatch(r"[0-9eEpilnogsqrt+\-*/(). ]+", src): raise ValueError("unsupported characters in point expression %r" % expr)
    loc = {"sqrt": sp.sqrt, "pi": sp.pi, "e": sp.E, "E": sp.E, "ln": sp.log, "log": sp.log}
    val = sp.sympify(src, locals=loc, rational=True)      # decimals become exact rationals
    if val.free_symbols: raise ValueError("point expression %r has free symbols %s" % (expr, val.free_symbols))
    v = mp.mpf(str(sp.N(val, mp.mp.dps + 10)))
    if not (0 < v < 1): print("   note: point %s = %s lies outside (0,1)" % (expr, mp.nstr(v, 12)), flush=True)
    return v
def point_tag(expr): return __import__("re").sub(r"[^0-9A-Za-z]+", "", expr.replace("/", "d").replace("sqrt", "r").replace(".", "p"))[:24]

# ---------------------------------------------------------------- atoms and entries
# atom = ("D", pkey, qkey) p - q | ("P", pkey) p | ("1P", pkey) 1 + p | ("XP", pkey) x0 - p | ("X0",) x0 | ("H", h) | ("1H", h) | ("DH", h, h')   (h = Z0 letters)
# entry = dict atom -> int exponent (a formal monomial); () empty = 1 (dropped)
def ent(*pairs):
    d = {}
    for a, e in pairs: d[a] = d.get(a, 0) + e
    return frozenset((a, e) for a, e in d.items() if e)
def sym_G(letters, diff, is_zero, is_one=None):
    """Goncharov symbol of G(a_1..a_n; 1): dict {tuple of entries: int}; letters are opaque keys; diff(u,v) -> entry of (u - v) or None if zero;
    a_0 = 1 ('ONE'), a_{n+1} = 0 ('ZERO')."""
    n = len(letters)
    if n == 0: return {(): 1}
    a = ["ONE"] + list(letters) + ["ZERO"]; out = collections.Counter()
    for i in range(1, n + 1):
        sub = sym_G(tuple(letters[:i - 1] + letters[i:]), diff, is_zero)
        for (uu, vv), sg in (((a[i - 1], a[i]), 1), ((a[i + 1], a[i]), -1)):
            e = diff(uu, vv)
            if e is None or len(e) == 0: continue
            for k, c in sub.items(): out[k + (e,)] += sg * c
    return {k: c for k, c in out.items() if c}
def shuffle(k1, k2):
    out = collections.Counter(); n, m = len(k1), len(k2)
    for pos in itertools.combinations(range(n + m), n):
        wl = [None] * (n + m); i1 = iter(k1); i2 = iter(k2); ps = set(pos)
        for i in range(n + m): wl[i] = next(i1) if i in ps else next(i2)
        out[tuple(wl)] += 1
    return out
def sym_mul(S1, S2):
    out = collections.Counter()
    for k1, c1 in S1.items():
        for k2, c2 in S2.items():
            for k, m in shuffle(k1, k2).items(): out[k] += c1 * c2 * m
    return {k: c for k, c in out.items() if c}

class SymBuilder:
    def __init__(self):
        self.cache = {}
    # --- points
    @staticmethod
    def pk(p): return json.dumps(p)
    def diff_points(self, u, v):
        """entry for u - v where u,v in {'ONE'-> handled by caller... } here u,v are ('pt', key) or 'ZERO'"""
        raise NotImplementedError
    def sym_Z(self, w):
        """symbol of ZIP_reg over points w: ZIP[a_1..a_n] = sum over letter choices {1 (coef -1), a/(1+a) (coef +1; absent if a = -1)} of G(b;1)"""
        key = ("Z", json.dumps(w))
        if key in self.cache: return self.cache[key]
        opts = []
        for p in w:
            o = [("ONE", -1)]
            if not (p[0] == "Q" and GT.deser(p[1]) == GT(-1)): o.append((("M", self.pk(p), tuple(map(tuple, [p])) and json.dumps(p)), 1))
            opts.append(o)
        S = collections.Counter()
        for combo in itertools.product(*opts):
            letters = tuple(l for l, s in combo); sign = 1
            for l, s in combo: sign *= s
            for k, c in sym_G(letters, self.diff_mob, None).items(): S[k] += sign * c
        S = {k: c for k, c in S.items() if c}; self.cache[key] = S; return S
    def diff_mob(self, u, v):
        """entries for differences of Moebius letters b = a/(1+a) (u, v in {'ONE', 'ZERO', ('M', pkey, pjson)})"""
        def isz(x): return x == "ZERO" or (x != "ONE" and json.loads(x[2])[0] == "Q" and GT.deser(json.loads(x[2])[1]).is_zero())
        if u == v: return None
        if u == "ONE" and v == "ONE": return None
        if isz(u) and isz(v): return None
        # b_u - b_v
        if u == "ONE":
            if isz(v): return ent()                                    # 1 - 0 = 1
            p = v[2]; return ent((("1P", p), -1))                      # 1 - a/(1+a) = 1/(1+a)
        if v == "ONE":
            if isz(u): return ent()                                    # 0 - 1 = -1 (torsion dropped)
            p = u[2]; return ent((("1P", p), -1))
        if isz(v):
            p = u[2]; return ent((("P", p), 1), (("1P", p), -1))       # a/(1+a)
        if isz(u):
            p = v[2]; return ent((("P", p), 1), (("1P", p), -1))
        p, q = u[2], v[2]
        if p == q: return None
        a, b = sorted([p, q]); return ent((("D", a, b), 1), (("1P", p), -1), (("1P", q), -1))   # (a_p - a_q)/((1+a_p)(1+a_q)), sign dropped
    def sym_G0(self, v, x0="11/2"):
        """symbol of G_0(v; x0) = G(v/x0; 1)"""
        key = ("G0", json.dumps(v), x0)
        if key in self.cache: return self.cache[key]
        letters = tuple(("S", json.dumps(p)) for p in v)
        def dz(u, w_):
            def isz(x): return x == "ZERO" or x == ("E", 0) or (x != "ONE" and json.loads(x[1])[0] == "Q" and GT.deser(json.loads(x[1])[1]).is_zero())
            if u == w_: return None
            if isz(u) and isz(w_): return None
            if u == "ONE":
                if isz(w_): return ent()
                return ent((("XP", w_[1], x0), 1), (("X0", x0), -1))            # 1 - p/x0 = (x0 - p)/x0
            if w_ == "ONE":
                if isz(u): return ent()
                return ent((("XP", u[1], x0), 1), (("X0", x0), -1))
            if isz(w_): return ent((("P", u[1]), 1), (("X0", x0), -1))
            if isz(u): return ent((("P", w_[1]), 1), (("X0", x0), -1))
            if u[1] == w_[1]: return None
            a, b = sorted([u[1], w_[1]]); return ent((("D", a, b), 1), (("X0", x0), -1))
        # trailing zero points: G(u, 0^k; x0) = sum_j c_j log(x0)^j/j! G(u_j; x0) (shuffle regularization, G(0;x) = log x), G(u_j; x0) = G(u_j/x0; 1)
        ZER = ("E", 0)
        toks = tuple(ZER if (json.loads(l[1])[0] == "Q" and GT.deser(json.loads(l[1])[1]).is_zero()) else l for l in letters)
        S = collections.Counter()
        for (j, u), c in hpath.decomp_trailing(toks, ZER).items():
            Su = sym_G(tuple(u), dz, None) if u else {(): 1}
            Lx = {((("LOGX0",),) * 0): 1} if j == 0 else {tuple(ent((("X0", x0), 1)) for _ in range(j)): 1}
            for k_, m_ in sym_mul(Su, Lx).items(): S[k_] += c * m_
        S = {k_: (int(v_) if getattr(v_, "denominator", 1) == 1 else v_) for k_, v_ in S.items() if v_ != 0}
        self.cache[key] = S; return S
    def sym_Z0(self, hw, x0="11/2"):
        """symbol of ZIP[h_1..h_n] with explicit Q(i)(t) letters (Z0 constants): atoms H, 1H, DH"""
        key = ("Z0", json.dumps(hw), x0)
        if key in self.cache: return self.cache[key]
        opts = []
        for h in hw:
            hj = json.dumps([h, x0]); o = [("ONE", -1)]
            if not (not h[1] and Fr(h[0][0], h[0][1]) == -1): o.append((("MH", hj), 1))
            opts.append(o)
        def dz(u, v):
            def isz(x): return x == "ZERO" or (x != "ONE" and (not json.loads(x[1])[0][1]) and json.loads(x[1])[0][0][0] == 0)
            if u == v: return None
            if isz(u) and isz(v): return None
            if u == "ONE":
                if isz(v): return ent()
                return ent((("1H", v[1]), -1))
            if v == "ONE":
                if isz(u): return ent()
                return ent((("1H", u[1]), -1))
            if isz(v): return ent((("H", u[1]), 1), (("1H", u[1]), -1))
            if isz(u): return ent((("H", v[1]), 1), (("1H", v[1]), -1))
            if u[1] == v[1]: return None
            a, b = sorted([u[1], v[1]]); return ent((("DH", a, b), 1), (("1H", u[1]), -1), (("1H", v[1]), -1))
        S = collections.Counter()
        for combo in itertools.product(*opts):
            letters = tuple(l for l, s in combo); sign = 1
            for l, s in combo: sign *= s
            for k, c in sym_G(letters, dz, None).items(): S[k] += sign * c
        S = {k: c for k, c in S.items() if c}; self.cache[key] = S; return S

# ---------------------------------------------------------------- the test
def run(label, path, weights=(1, 2, 3), dps=60, lll_dps=300, table_only=True, tstar_name="e/7", tstar2_name="pi/9", continuation=True):
    doc = json.load(gzip.open(path, "rt")); letters = doc["letters"]
    X = sp.Symbol("X")
    # quadratic letters and their discriminant square classes over Q[t] (real coefficients) — group into field classes
    def gt2sp(ser):
        g = GT.deser(ser); f = lambda P: sum(sp.Rational(int(c.p), int(c.q)) * T ** k for k, c in enumerate(P.coeffs()))
        return (f(g.A) + sp.I * f(g.B)) / f(g.D)
    def sqclass(expr):
        num, den = sp.fraction(sp.cancel(sp.together(expr))); ker = []
        for pol in (num, den):
            P = sp.Poly(pol, T)
            if P.degree() <= 0:
                c = sp.Rational(P.as_expr())   # rational constant: square class over Q matters only for sign/-1; keep sign
                if c < 0: ker.append("-1")
                continue
            cst, facs = sp.factor_list(P)
            if sp.Rational(cst) < 0: ker.append("-1")
            for f_, e in facs:
                if int(e) % 2: ker.append(str(sp.Poly(f_, T).primitive()[1].as_expr()))
        return tuple(sorted(ker))
    quad = {}; cub = {}
    for L, info in letters.items():
        if info["deg_x"] == 2:
            cs = [gt2sp(c) for c in info["x_coeffs(GT ser, increasing powers of X)"]]
            real = all(GT.deser(c).is_real() for c in info["x_coeffs(GT ser, increasing powers of X)"])
            disc = sp.cancel(cs[1] ** 2 - 4 * cs[2] * cs[0])
            quad[L] = {"poly": info["poly"], "real": real, "class": sqclass(disc) if real else ("GAUSSIAN", str(disc)[:80])}
        elif info["deg_x"] == 3:
            cs = [gt2sp(c) for c in info["x_coeffs(GT ser, increasing powers of X)"]]; real = all(GT.deser(c).is_real() for c in info["x_coeffs(GT ser, increasing powers of X)"])
            disc = sp.cancel(sp.discriminant(sp.Poly(sum(c * X ** k for k, c in enumerate(cs)), X), X)) if real else None
            cub[L] = {"poly": info["poly"], "real": real, "class": sqclass(disc) if real else None}
    classes = collections.defaultdict(list)
    for L, q in quad.items(): classes[q["class"]].append(L)
    results = {"label": label, "file": path, "quadratic_letters": {L: {"poly": q["poly"], "class": q["class"]} for L, q in quad.items()}, "cubic_letters": {L: {"poly": c["poly"], "disc_class": c["class"]} for L, c in cub.items()}, "classes": {}}
    with mp.workdps(lll_dps): tstar = point_value(tstar_name); tstar2 = point_value(tstar2_name)      # generic points at lll_dps precision (E6/P2: any expression; recorded below)
    results["points"] = {"tstar": tstar_name, "tstar2": tstar2_name, "tstar_value": mp.nstr(tstar, 40), "tstar2_value": mp.nstr(tstar2, 40), "precision_digits": lll_dps, "slot_bijection_tstar_to_tstar2": "root continuation (ContNum, 96 steps, chordal matching)" if continuation else "canonical order at each point independently", "root_order": "canonical: X-coefficients and polyroots at %d digits, sorted by (Re@1e-25, Im)" % Num.ROOT_DPS}
    SB = SymBuilder()
    for cls, Ls in classes.items():
        if cls and cls[0] == "GAUSSIAN": results["classes"][str(cls)] = {"letters": Ls, "note": "quadratic letter with non-real coefficients: not tested"}; continue
        Lset = set(Ls)
        cub_conflict = [L for L, c in cub.items() if c["class"] == cls]
        rec = {"letters": [quad[L]["poly"] for L in Ls], "sqrt_class(odd factors)": list(cls), "cubic_letters_with_same_disc_class": cub_conflict, "weights": {}}
        # terms that involve a root of a class letter
        def involves(tm):
            if any(g[0] in Lset for g in tm["gens"]): return True
            for w in tm["Z"]:
                if any(p[0] == "A" and p[1] in Lset for p in w): return True
            for a in tm["const"]:
                if a[0] == "G0" and any(p[0] == "A" and p[1] in Lset for p in a[1]): return True
            return False
        terms = [tm for tm in doc["terms"] if involves(tm)]
        rec["n_terms_involving_class"] = len(terms)
        mp.mp.dps = lll_dps
        N1 = Num(doc, tstar, lll_dps); N2 = ContNum(doc, tstar2, lll_dps, ref=N1) if continuation else Num(doc, tstar2, lll_dps)
        # tau image of a point: swap slot for class letters
        def tau_pt(pj):
            p = json.loads(pj)
            if p[0] == "A" and p[1] in Lset: return json.dumps(["A", p[1], 1 - p[2]])
            return pj
        # odd atoms: for an atom gamma with tau(gamma) != gamma: u = gamma/tau(gamma); canonical representative: the lexicographically smaller of (gamma, tau gamma) with sign
        def tau_atom(a):
            if a[0] == "D": x, y = tau_pt(a[1]), tau_pt(a[2]); x2, y2 = sorted([x, y]); return ("D", x2, y2)
            if a[0] in ("P", "1P"): return (a[0], tau_pt(a[1]))
            if a[0] == "XP": return ("XP", tau_pt(a[1]), a[2])
            return a
        def atom_val(a, N):
            if a[0] == "D": return N.point(json.loads(a[1])) - N.point(json.loads(a[2]))
            if a[0] == "P": return N.point(json.loads(a[1]))
            if a[0] == "1P": return 1 + N.point(json.loads(a[1]))
            if a[0] == "XP": return mp.mpf(Fr(a[2]).numerator) / Fr(a[2]).denominator - N.point(json.loads(a[1]))
            if a[0] == "X0": return mp.mpc(mp.mpf(Fr(a[1]).numerator) / Fr(a[1]).denominator)
            def hval(hs): h_, x0_ = json.loads(hs); return N.hl(h_, mp.mpf(Fr(x0_).numerator) / Fr(x0_).denominator)
            if a[0] == "H": return hval(a[1])
            if a[0] == "1H": return 1 + hval(a[1])
            if a[0] == "DH": return hval(a[1]) - hval(a[2])
            raise ValueError(a)
        # NOTE: ("D", x, y) after tau may flip order -> gamma/tau(gamma) picks up a sign (torsion): irrelevant for oddvec (torsion) and for log|.|
        odd_atoms = {}      # canonical odd generator key -> (gamma, taugamma)
        def odd_key(a):
            ta = tau_atom(a)
            if ta == a: return None, 0
            k1, k2 = repr(a), repr(ta)
            if k1 < k2: odd_atoms.setdefault(k1, (a, ta)); return k1, 1
            odd_atoms.setdefault(k2, (ta, a)); return k2, -1
        # pass 1: collect symbols of all involved terms, per weight
        per_w = {w_: [] for w_ in weights}
        t0 = time.time()
        for tm in terms:
            S = {(): 1}
            for wz in tm["Z"]: S = sym_mul(S, SB.sym_Z(wz))
            x0s = "%d/%d" % tuple(tm.get("x0", [11, 2]))
            for a in tm["const"]:
                S = sym_mul(S, SB.sym_G0(a[1], x0s) if a[0] == "G0" else SB.sym_Z0(a[1], x0s))
            wgt = sum(len(wz) for wz in tm["Z"]) + sum(len(a[1]) for a in tm["const"])
            if wgt in per_w and S: per_w[wgt].append((tm, S))
        rec["symbol_build_s"] = round(time.time() - t0, 1)
        # collect odd atoms
        for w_, lst in per_w.items():
            for tm, S in lst:
                for key in S:
                    for e in key:
                        for a, n in e: odd_key(a)
        oa = sorted(odd_atoms)
        rec["n_odd_atoms"] = len(oa)
        if not oa: rec["verdict"] = "no tau-moved atom occurs: sqrt letters ABSENT trivially"; results["classes"][str(cls)] = rec; continue
        # LLL relations among u_k = gamma/tau gamma (complex logs) + torsion
        def uval(k, N): g, tg = odd_atoms[k]; return atom_val(g, N) / atom_val(tg, N)
        with mp.workdps(lll_dps):
            logs = [mp.log(uval(k, N1)) for k in oa]; LAM = mp.mpf(10) ** (lll_dps - 60); tors = 2 * mp.pi / 24
            n = len(oa); rows = []
            for jn in range(n):
                row = [0] * (n + 1); row[jn] = 1; rows.append(row + [int(mp.nint(LAM * mp.re(logs[jn]))), int(mp.nint(LAM * mp.im(logs[jn])))])
            rows.append([0] * n + [1] + [0, int(mp.nint(LAM * tors))])
            M = fmpz_mat(rows).lll(); rels = []; n_pass1 = 0; failed2 = []
            logs2 = [mp.log(uval(k, N2)) for k in oa]
            for i in range(M.nrows()):
                r = [int(M[i, jn]) for jn in range(n + 3)]; cf = r[:n]; kk = r[n]
                if all(c == 0 for c in cf) or max(abs(c) for c in cf) > 10 ** 4: continue
                res1 = abs(sum(c * l for c, l in zip(cf, logs)) + kk * mp.mpc(0, 1) * tors)
                if res1 > mp.mpf(10) ** (-(lll_dps - 80)): continue
                n_pass1 += 1
                # verify at t** (imaginary part mod 2 pi i /24)
                v2 = sum(c * l for c, l in zip(cf, logs2)); re2 = abs(mp.re(v2)); im2 = mp.im(v2) / tors; im2 = abs(im2 - mp.nint(im2))
                if re2 < mp.mpf(10) ** (-(lll_dps - 80)) and im2 < mp.mpf(10) ** (-(lll_dps - 90)): rels.append(cf)
                else: failed2.append({"cf": cf, "re2": mp.nstr(re2, 3), "im2_frac": mp.nstr(im2, 3)})
        rec["LLL"] = {"n_relations_found_at_tstar": n_pass1, "n_verified_at_tstar2": len(rels), "failed_verification_at_tstar2": failed2[:40], "letters_reordered_by_continuation_at_tstar2": (list(map(lambda L_: quad.get(L_, cub.get(L_, {})).get("poly", L_)[:60], N2.reordered)) if continuation else None)}
        if n_pass1 != len(rels): print("   WARNING: %d of %d LLL relations failed verification at t** (labeling or accidental relation)" % (n_pass1 - len(rels), n_pass1), flush=True)
        rec["n_relations_among_odd_atoms(LLL at t*, verified at t**)"] = len(rels); rec["odd_atoms"] = [repr(odd_atoms[k][0])[:200] for k in oa]; rec["relations"] = rels
        # free quotient: rref of relation matrix over Q -> pivot atoms eliminated
        import sympy as _sp
        if rels:
            Mr = _sp.Matrix(rels); rr, piv = Mr.rref(); piv = list(piv)
            elim = {}
            for irow, pc in enumerate(piv):
                rowv = rr.row(irow); elim[pc] = {jn: -rowv[jn] for jn in range(n) if jn != pc and rowv[jn] != 0}
        else: elim = {}
        free = [jn for jn in range(n) if jn not in elim]; rec["n_free_odd_atoms"] = len(free)
        def oddvec(e):
            """entry (frozenset of (atom, exp)) -> dict free-index -> Fraction (1/2 * sum n * sign * [u])"""
            v = collections.defaultdict(lambda: Fr(0))
            for a, nexp in e:
                k, sg = odd_key(a)
                if k is None: continue
                jn = oa.index(k); coef = Fr(nexp * sg, 2)
                if jn in elim:
                    for j2, c2 in elim[jn].items(): v[j2] += coef * Fr(int(c2.p), int(c2.q)) if hasattr(c2, "p") else coef * Fr(c2)
                else: v[jn] += coef
            return {jn: c for jn, c in v.items() if c != 0}
        def logabs(e, N):
            s = mp.mpf(0)
            for a, nexp in e: s += nexp * mp.log(abs(atom_val(a, N)))
            return s
        mp.mp.dps = dps
        for w_, lst in per_w.items():
            best = {}; norms = {}
            for N, nm in ((N1, "t*=" + tstar_name), (N2, "t**=" + tstar2_name)):
                # no cache reset: slot order is canonical/continued and precision-independent, so the Phi pass uses the SAME bijection as the LLL pass
                Phi = [collections.defaultdict(lambda: mp.mpc(0)) for _ in range(w_)]; scale = mp.mpf(0)
                for tm, S in lst:
                    c = mp.mpc(0)
                    for cs_ in tm["coef"]: c += N.gt(cs_)
                    for Lg, s, e in tm["gens"]: c *= N.roots(Lg)[s] ** e
                    for key, mult in S.items():
                        la = [logabs(e, N) for e in key]
                        for p in range(w_):
                            ov = oddvec(key[p])
                            if not ov: continue
                            rest = mp.mpf(1)
                            for q in range(w_):
                                if q != p: rest *= la[q]
                            for jn, fc in ov.items():
                                val = c * mult * (mp.mpf(fc.numerator) / fc.denominator) * rest; Phi[p][jn] += val; scale = max(scale, abs(val))
                mx = max((abs(v) for P in Phi for v in P.values()), default=mp.mpf(0))
                norms[nm] = {"Phi_per_slot_atom": {"%d:%d" % (p, jn): mp.nstr(v, 8) for p, P in enumerate(Phi) for jn, v in P.items()}, "max|Phi| over slots and free odd atoms": mp.nstr(mx, 5), "scale(max single contribution)": mp.nstr(scale, 5), "ratio": mp.nstr(mx / scale, 5) if scale > 0 else None,
                             "n_nonzero(>1e-%d*scale)" % (dps // 2): sum(1 for P in Phi for v in P.values() if abs(v) > scale * mp.mpf(10) ** (-(dps // 2)))}
            verdict = "SURVIVES" if any(float(v_["ratio"] or 0) > 10.0 ** (-(dps // 2)) for v_ in norms.values()) else ("CANCELS (all functionals zero to %d digits at both points)" % (dps // 2) if lst else "no terms of this weight")
            rec["weights"][str(w_)] = {"n_terms": len(lst), "functionals": norms, "verdict": verdict}
            print("   class %s weight %d: %s  %s" % (list(cls), w_, verdict, {k: v["ratio"] for k, v in norms.items()}), flush=True)
        results["classes"][str(cls)] = rec
    return results
if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--files", default="quark"); ap.add_argument("--dps", type=int, default=60)
    ap.add_argument("--tstar", default="e/7", help="generic point t*: a name %s or any expression, e.g. 'sqrt(3)/5', '2/(3+sqrt5)', '0.377', '3/7'" % sorted(POINTS)); ap.add_argument("--tstar2", default="pi/9", help="second generic point t** (same syntax)"); ap.add_argument("--no-continuation", dest="cont", action="store_false")
    ap.add_argument("--out", default="", help="receipt name (default SLICE_S5_SYMTEST.json for (e/7, pi/9), else SLICE_S5_SYMTEST_<t*>_<t**>.json)")
    ap.add_argument("--data", default="", help="directory holding out/*.json.gz (same meaning as tests/run_tests.py --data and $SURD_DATA); the receipt goes to $SURD_WORK (default ./surd_work)")
    a = ap.parse_args()
    if a.data: OUT = os.path.join(os.path.abspath(a.data), "out")
    slice_prov.log_pid("symtest"); files = {"quark": "E4C_LO_QCD_dipole_slice_quark.json.gz", "gluon": "E4C_LO_QCD_dipole_slice_gluon.json.gz", "n4": "E4C_LO_N4_dipole_slice_n4.json.gz"}
    for ch_ in ("q_qbpqpgq", "q_qbqgq", "q_gggq", "g_qbpqpqbq", "g_qbqqbq", "g_qbggq", "g_gggg"): files[ch_] = "E4C_LO_QCD_dipole_slice_%s.json.gz" % ch_
    out = {}
    tag = lambda x: x.replace("/", "").replace("sqrt", "r") if x in POINTS else point_tag(x)
    rname = a.out or ("SLICE_S5_SYMTEST.json" if (a.tstar, a.tstar2) == ("e/7", "pi/9") else "SLICE_S5_SYMTEST_%s_%s.json" % (tag(a.tstar), tag(a.tstar2)))
    path = os.path.join(WORK, rname)
    if os.path.exists(path):
        try: out = json.load(open(path)); out.pop("producer", None)
        except Exception: out = {}
    for lab in a.files.split(","):
        p = files[lab] if lab in files else lab; p = p if os.path.isabs(p) else os.path.join(OUT, p); print("==", lab, "(t*=%s, t**=%s, continuation=%s)" % (a.tstar, a.tstar2, a.cont), flush=True); t0 = time.time(); r = run(lab, p, dps=a.dps, tstar_name=a.tstar, tstar2_name=a.tstar2, continuation=a.cont); r["wall_s"] = round(time.time() - t0, 1); out[lab] = r
    slice_prov.write_receipt(path, out, inputs=[os.path.join(OUT, f) for f in files.values() if os.path.exists(os.path.join(OUT, f))])
    print("wrote", path)
