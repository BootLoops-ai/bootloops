"""fibr_gi — COPY of symbolic/scripts/fibr.py with all polynomial arithmetic done in R = Q[x...,t,j]/(j^2+1)
(j = sqrt(-1)): every product is reduced mod j^2+1 (Engine.red), letter cancellation uses exact division in R (GIAlphabet.gi_divides), and the
alphabet is alphabet_gi.GIAlphabet (irreducible letters over Q(i), auto-extending).  Everything else is verbatim.

ORIGINAL DOCSTRING: fibr — fibration-basis hyperlogarithm integrator over (0,oo) for integrands
      R(x; other vars) x prod ZIP-words,   R = N / prod letters^e  (letters from a fixed ALPHABET, linear in x),
with sparse fmpq_mpoly numerators and denominators kept as EXPONENT VECTORS over the alphabet (never expanded, never GCD'd
across supports).  Semantics of words: ZIP[a1,...,an] = reg. int_0^oo dx1/(x1-a1) int_0^x1 dx2/(x2-a2) ... (HyperFLINT
ZeroInfPeriod, shuffle-regularized: log-divergences at oo and 0 dropped in the scale {x^n log^j x}); a term's value is
N/D * prod_words ZIP[word].  Hyperlog/pole letters are ratios of alphabet monomials  c * prod L^e  (e in Z), so equality
is decidable by key.

Term store TS: dict {(Dkey, Wkey): N}   Dkey = tuple(sorted((idx, e>0))), Wkey = tuple(sorted(words)), word = tuple(HL),
HL = ((p, q), tuple(sorted((idx, e))))  [c = p/q; zero letter = ((0,1),())].

Primitives:
  pf            partial fractions in x by Laurent expansion at each root (denominators: resultants a b'-a' b and leading
                coefficients, factored over the alphabet -> refuse if absent) + polynomial part by expansion at oo;
  activate      ZIP[s(x)] (weight-1 word whose letter depends on x) -> ZIP[s(0)] - sum_i e_i G(rho_i; x);
  J             reg. int_0^oo (x-s)^-k G(w;x) dx  (k=1: ZIP[s,w]; k>=2: integration by parts, boundary expansions at 0);
  I             reg. int_0^oo x^m G(w;x) dx  (expansion of G at oo in x^-i log^j x with ZIP constants);
  integrate     one Brown-linear elimination of a whole TS; consolidate = common denominator ONLY within an identical
                letter support (no new letters), then exact cancellation of letter factors by trial division."""
import sys, os, time, math
from fractions import Fraction as Fr
from flint import fmpq, fmpq_mpoly
from alphabet import Alphabet, OffAlphabet
from alphabet_gi import GIAlphabet

ZERO_HL = ((0, 1), ())

def binom(n, k):
    """generalized binomial for integer n (may be negative), k>=0"""
    r = Fr(1)
    for i in range(k):
        r = r * (n - i) / (i + 1)
    return r

class Engine:
    def __init__(self, A: Alphabet):
        self.A = A; self.ctx = A.ctx; self.nv = len(A.names)
        self.one = self.ctx.from_dict({(0,) * self.nv: fmpq(1)})
        self._powcache = {}; self._hlpoly = {}; self._diffcache = {}; self._res = {}
        self.stats = {"factor_calls": 0, "new_denominators": {}}
        self.sign_checker = None     # optional callable(hl or poly) for domain-sign assertions
        self.use_pf2 = True          # partial fractions over one common denominator per root (pf = reference implementation)
        self.red = A.red if hasattr(A, 'red') else (lambda P: P)
        self.gi = hasattr(A, 'gi_divides')

    # ------------------------------------------------------------ basic helpers
    def cancel_letters(self, P, D):
        """exact cancellation of denominator letters from numerator P by trial division (with the evaluation pre-test)"""
        for l in list(D):
            L = self.A.letters[l]
            if self.gi:
                k, P = self.A.gi_divide_maxpow(P, L, D.get(l, 0))
                if k:
                    D[l] -= k
                    if D[l] == 0: del D[l]
                continue
            while D.get(l, 0) > 0:          # flint divmod exits fast when L does not divide P (measured 0.01 s on 152k terms)
                q, r = divmod(P, L)
                if r.is_zero():
                    P = q; D[l] -= 1
                    if D[l] == 0: del D[l]
                else: break
        return P, D
    def qq(self, f: Fr):
        return fmpq(f.numerator, f.denominator)
    def lpow(self, idx, e):
        k = (idx, e)
        if k not in self._powcache:
            self._powcache[k] = self._rpow(self.A.letters[idx], e)
        return self._powcache[k]
    def _rpow(self, p, e):
        r = self.one; b = self.red(p)
        while e:
            if e & 1: r = self.red(r * b)
            e >>= 1
            if e: b = self.red(b * b)
        return r
    def ppow(self, p, e, tag):
        k = ("P", tag, e)
        if k not in self._powcache:
            self._powcache[k] = self._rpow(p, e)
        return self._powcache[k]
    @staticmethod
    def dkey(D):
        return tuple(sorted((i, e) for i, e in D.items() if e))
    def factor(self, p, context):
        self.stats["factor_calls"] += 1
        c, E = self.A.factor(p, context)
        for i in E:
            self.stats["new_denominators"].setdefault(context.split(":")[0], set()).add(i)
        return c, E

    # ------------------------------------------------------------ TS algebra
    @staticmethod
    def ts_add_into(acc, ts, scale=None):
        for k, N in ts.items():
            if scale is not None:
                N = N * scale
            v = acc.get(k)
            if v is None:
                acc[k] = N
            else:
                v = v + N
                if v.is_zero():
                    del acc[k]
                else:
                    acc[k] = v
        return acc
    def ts_scale(self, ts, f: Fr):
        if f == 0: return {}
        q = self.qq(f)
        return {k: N * q for k, N in ts.items()}
    def ts_mul_mono(self, ts, c: Fr, E: dict):
        """ts * (c * prod letters^E), E signed"""
        if c == 0: return {}
        out = {}
        q = self.qq(c)
        for (Dk, Wk), N in ts.items():
            D = dict(Dk); P = None
            for i, e in E.items():
                if e < 0:
                    D[i] = D.get(i, 0) - e
                elif e > 0:
                    have = D.get(i, 0)
                    cnl = min(have, e)
                    if cnl:
                        if have - cnl: D[i] = have - cnl
                        else: del D[i]
                    r = e - cnl
                    if r:
                        P = self.lpow(i, r) if P is None else self.red(P * self.lpow(i, r))
            N2 = N * q
            if P is not None: N2 = self.red(N2 * P)
            k = (self.dkey(D), Wk)
            v = out.get(k)
            if v is None: out[k] = N2
            else:
                v = v + N2
                if v.is_zero(): del out[k]
                else: out[k] = v
        return out
    def ts_mul_poly(self, ts, P):
        if P.is_zero(): return {}
        out = {}
        for k, N in ts.items():
            v = self.red(N * P)
            if not v.is_zero(): out[k] = v
        return out
    MINUS_ONE_WORD = (((-1, 1), ()),)
    @staticmethod
    def ts_attach(ts, words):
        """multiply by prod ZIP[word] for word in words (ZIP[-1] = -log 1 = 0 kills the term)"""
        if any(tuple(w) == Engine.MINUS_ONE_WORD for w in words): return {}
        out = {}
        for (Dk, Wk), N in ts.items():
            k = (Dk, tuple(sorted(Wk + tuple(words))))
            v = out.get(k)
            out[k] = N if v is None else v + N
        return {k: v for k, v in out.items() if not v.is_zero()}
    def ts_one(self):
        return {((), ()): self.one}
    def ts_mono(self, c, E):
        return self.ts_mul_mono(self.ts_one(), c, E)
    @staticmethod
    def ts_nmon(ts):
        return sum(len(N) for N in ts.values())

    # ------------------------------------------------------------ hyperlog letters (ratios of alphabet monomials)
    def hl(self, c: Fr, E: dict):
        if c == 0: return ZERO_HL
        return ((c.numerator, c.denominator), tuple(sorted((i, e) for i, e in E.items() if e)))
    @staticmethod
    def hl_c(h): return Fr(h[0][0], h[0][1])
    @staticmethod
    def hl_E(h): return dict(h[1])
    def hl_from_alpha_beta(self, alpha, beta, context):
        """sigma = -beta/alpha"""
        if beta.is_zero(): return ZERO_HL
        cb, Eb = self.factor(beta, "beta:" + context); ca, Ea = self.factor(alpha, "alpha:" + context)
        E = dict(Eb)
        for i, e in Ea.items(): E[i] = E.get(i, 0) - e
        return self.hl(-cb / ca, E)
    def hl_polys(self, h):
        """(num poly incl. constant, den poly) of the letter"""
        if h in self._hlpoly: return self._hlpoly[h]
        c = self.hl_c(h); E = self.hl_E(h)
        num = self.one * self.qq(c); den = self.one
        for i, e in E.items():
            if e > 0: num = self.red(num * self.lpow(i, e))
            else: den = self.red(den * self.lpow(i, -e))
        self._hlpoly[h] = (num, den)
        return num, den
    def hl_diff_mono(self, a, s):
        """a - s as (c, E) monomial over the alphabet (factor the numerator); raises OffAlphabet"""
        k = (a, s)
        if k in self._diffcache: return self._diffcache[k]
        na, da = self.hl_polys(a); ns, ds = self.hl_polys(s)
        Ea = {i: -e for i, e in self.hl_E(a).items() if e < 0}; Es = {i: -e for i, e in self.hl_E(s).items() if e < 0}
        L = dict(Ea)
        for i, e in Es.items(): L[i] = max(L.get(i, 0), e)
        ca = self.one; cs = self.one
        for i, e in L.items():
            if e - Ea.get(i, 0): ca = self.red(ca * self.lpow(i, e - Ea.get(i, 0)))
            if e - Es.get(i, 0): cs = self.red(cs * self.lpow(i, e - Es.get(i, 0)))
        num = self.red(na * ca - ns * cs)
        if num.is_zero():
            res = (Fr(0), {})
        else:
            c, E = self.factor(num, "diff:%s|%s" % (self.hl_str(a), self.hl_str(s)))
            E = dict(E)
            for i, e in L.items(): E[i] = E.get(i, 0) - e
            res = (c, E)
        self._diffcache[k] = res
        return res
    def hl_str(self, h):
        if h == ZERO_HL: return "0"
        c = self.hl_c(h); E = self.hl_E(h)
        num = [str(c)] + ["(%s)^%d" % (self.A.describe(i), e) for i, e in E.items() if e > 0]
        den = ["(%s)^%d" % (self.A.describe(i), -e) for i, e in E.items() if e < 0]
        s = "*".join(num)
        return s + ("/(" + "*".join(den) + ")" if den else "")
    def hl_depends(self, h, var):
        return any(self.A.depends(i, var) for i, e in h[1])

    # ------------------------------------------------------------ partial fractions in var
    def resultant_mono(self, i, j, var):
        """R_ij = alpha_i beta_j - alpha_j beta_i, factored: (c, E)"""
        k = (i, j, var)
        if k in self._res: return self._res[k]
        ai, bi = self.A.linear_parts(i, var); aj, bj = self.A.linear_parts(j, var)
        R = self.red(ai * bj - aj * bi)
        if R.is_zero(): raise ValueError("zero resultant between letters %d %d" % (i, j))
        res = self.factor(R, "res:%s|%s|%s" % (var, self.A.describe(i)[:40], self.A.describe(j)[:40]))
        self._res[k] = res
        return res
    def pf(self, Ncoef, dens, var, poly_only=False):
        """N(x)/prod_i L_i^k_i, N = sum_n Ncoef[n] x^n (Ncoef[n]: TS, x-free), dens = {idx: k}.
        Returns (polypart {p: TS}, poles {(idx, s): TS}) meaning  sum_p polypart[p] x^p + sum poles[(i,s)] / L_i(x)^s."""
        poles = {}; polypart = {}
        idxs = sorted(dens)
        lin = {i: self.A.linear_parts(i, var) for i in idxs}
        afac = {}; 
        for i in idxs:
            afac[i] = self.factor(lin[i][0], "alpha:%s|%s" % (var, self.A.describe(i)[:40]))
        degN = max(Ncoef) if Ncoef else -1
        K = sum(dens.values())
        for i in (idxs if not poly_only else []):
            ki = dens[i]; alpha, beta = lin[i]; ca, Ea = afac[i]
            # S_N(t): coefficient of t^m, m < ki:  sum_{n>=m} N_n binom(n,m) (-beta)^(n-m) alpha^(-n)
            SN = [dict() for _ in range(ki)]
            for n, Nn in Ncoef.items():
                inv_an = (Fr(1) / ca ** n, {l: -e * n for l, e in Ea.items()}) if n else (Fr(1), {})
                base = self.ts_mul_mono(Nn, inv_an[0], inv_an[1]) if n else Nn
                for m in range(min(n, ki - 1) + 1):
                    coef = binom(n, m) * (-1) ** (n - m)
                    if n - m == 0:
                        t = self.ts_scale(base, coef)
                    elif beta.is_zero():
                        continue
                    else:
                        t = self.ts_mul_poly(self.ts_scale(base, coef), self.ppow(beta, n - m, ("beta", i, var)))
                    self.ts_add_into(SN[m], t)
            # multiply by prod_{j != i} alpha_i^{k_j} (R_ij + alpha_j t)^(-k_j)
            S = SN
            Kp = sum(dens[j] for j in idxs if j != i)
            if Kp:
                S = [self.ts_mul_mono(s, ca ** Kp, {l: e * Kp for l, e in Ea.items()}) for s in S]
            for j in idxs:
                if j == i: continue
                kj = dens[j]; cR, ER = self.resultant_mono(i, j, var); aj = lin[j][0]
                # series coefficients: binom(-kj, m) alpha_j^m R^(-kj-m)
                ser = []
                for m in range(ki):
                    cm = binom(-kj, m) / cR ** (kj + m)
                    Em = {l: -e * (kj + m) for l, e in ER.items()}
                    ser.append((cm, Em, m))
                newS = [dict() for _ in range(ki)]
                for m1, s1 in enumerate(S):
                    if not s1: continue
                    for (cm, Em, m) in ser:
                        if m1 + m >= ki: break
                        t = self.ts_mul_mono(s1, cm, Em)
                        if m: t = self.ts_mul_poly(t, self.ppow(aj, m, ("alpha", j, var)))
                        self.ts_add_into(newS[m1 + m], t)
                S = newS
            for s in range(1, ki + 1):
                if S[ki - s]:
                    poles[(i, s)] = S[ki - s]
        if degN >= K:
            # expansion at infinity: N_n x^(n-K) prod_i alpha_i^-k_i (1 + beta_i/(alpha_i x))^-k_i
            base_c = Fr(1); base_E = {}
            for i in idxs:
                ca, Ea = afac[i]; base_c /= ca ** dens[i]
                for l, e in Ea.items(): base_E[l] = base_E.get(l, 0) - e * dens[i]
            # series in y = 1/x: prod_i sum_m binom(-k_i, m) beta_i^m alpha_i^-m y^m, truncated at order degN - K
            M = degN - K
            ser = {0: self.ts_mono(base_c, base_E)}
            for i in idxs:
                ca, Ea = afac[i]; alpha, beta = lin[i]
                if beta.is_zero(): continue
                fac = []
                for m in range(M + 1):
                    fac.append((binom(-dens[i], m) / ca ** m, {l: -e * m for l, e in Ea.items()}, m))
                new = {}
                for o, t0 in ser.items():
                    for (cm, Em, m) in fac:
                        if o + m > M: break
                        t = self.ts_mul_mono(t0, cm, Em)
                        if m: t = self.ts_mul_poly(t, self.ppow(beta, m, ("beta", i, var)))
                        self.ts_add_into(new.setdefault(o + m, {}), t)
                ser = new
            for n, Nn in Ncoef.items():
                for o, t0 in ser.items():
                    p = n - K - o
                    if p < 0: continue
                    # t0 * Nn : TS times TS(pure numerators, W=()) -> multiply numerators pairwise
                    prod = {}
                    for (Dk, Wk), N0 in t0.items():
                        for (Dk2, Wk2), N1 in Nn.items():
                            D = dict(Dk)
                            for l, e in Dk2: D[l] = D.get(l, 0) + e
                            self.ts_add_into(prod, {(self.dkey(D), tuple(sorted(Wk + Wk2))): self.red(N0 * N1)})
                    self.ts_add_into(polypart.setdefault(p, {}), prod)
        return polypart, poles

    # ------------------------------------------------------------ partial fractions: one common denominator per root (no intermediate key explosion)
    def pf2(self, Npoly, dens, var):
        """Same contract as pf but the input numerator coefficients are plain polynomials Npoly = {n: fmpq_mpoly (x-free)} and every
        root i is processed over ONE common denominator  alpha_i^a * prod_j R_ij^(k_j+k_i-1)  (R, alpha factored over the alphabet),
        so the sums over n and over the series orders merge monomials immediately.  Output: (polypart {p: TS}, poles {(i,s): TS})."""
        poles = {}; polypart = {}
        idxs = sorted(dens)
        lin = {i: self.A.linear_parts(i, var) for i in idxs}
        afac = {i: self.factor(lin[i][0], "alpha:%s|%s" % (var, self.A.describe(i)[:40])) for i in idxs}
        degN = max(Npoly) if Npoly else -1
        K = sum(dens.values())
        zero = self.ctx.from_dict({(0,) * self.nv: fmpq(0)})
        for i in idxs:
            ki = dens[i]; alpha, beta = lin[i]; ca, Ea = afac[i]
            nmax = degN
            # S_N(t) = sum_n N_n alpha^(nmax-n) (t - beta)^n  truncated at t^(ki-1)   [ / alpha^nmax ],  by Horner in n:
            #   S <- S*(t - beta) + N_n alpha^(nmax-n)   (n = nmax .. 0); each step multiplies by the SMALL beta only.
            import time as _t; _t0 = _t.time()
            SN = [zero for _ in range(ki)]
            for n in range(nmax, -1, -1):
                # S *= (t - beta): new[m] = S[m-1] - beta*S[m]
                newS = [zero for _ in range(ki)]
                for m in range(ki):
                    acc = zero
                    if m >= 1 and not SN[m - 1].is_zero(): acc = SN[m - 1]
                    if not SN[m].is_zero() and not beta.is_zero(): acc = acc - self.red(SN[m] * beta)
                    newS[m] = acc
                SN = newS
                Nn = Npoly.get(n)
                if Nn is not None:
                    SN[0] = SN[0] + (self.red(Nn * self.ppow(alpha, nmax - n, ("alpha", i, var))) if nmax - n else Nn)
            self.stats.setdefault("t_SN", 0.0); self.stats["t_SN"] += _t.time() - _t0; _t0 = _t.time()
            # denominators so far: alpha_i^nmax ; numerator factor alpha_i^Kp from the other roots
            Kp = sum(dens[j] for j in idxs if j != i)
            Dc = Fr(1); DE = {}
            net = nmax - Kp                      # alpha_i^(-net)
            if net > 0:
                Dc *= ca ** net
                for l, e in Ea.items(): DE[l] = DE.get(l, 0) + e * net
            elif net < 0:
                P = self.ppow(alpha, -net, ("alpha", i, var)); SN = [self.red(c * P) for c in SN]
            S = SN
            for j in idxs:
                if j == i: continue
                kj = dens[j]; aj = lin[j][0]
                R = self.red(lin[i][0] * lin[j][1] - lin[j][0] * lin[i][1])
                cR, ER = self.resultant_mono(i, j, var)          # validates R over the alphabet (refuse otherwise)
                top = kj + ki - 1                                  # common power R^top
                ser = []
                for m in range(ki):
                    t = self.ppow(aj, m, ("alpha", j, var)) if m else self.one
                    if top - kj - m: t = self.red(t * self.ppow(R, top - kj - m, ("R", i, j, var)))
                    ser.append(t * self.qq(binom(-kj, m)))
                newS = [zero for _ in range(ki)]
                for m1 in range(ki):
                    if S[m1].is_zero(): continue
                    for m in range(ki - m1):
                        newS[m1 + m] = newS[m1 + m] + self.red(S[m1] * ser[m])
                S = newS
                self.stats.setdefault("t_Rser", 0.0); self.stats["t_Rser"] += _t.time() - _t0; _t0 = _t.time()
                # denominator R^top = (cR prod L^ER)^top
                Dc *= cR ** top
                for l, e in ER.items(): DE[l] = DE.get(l, 0) + e * top
            inv = Fr(1) / Dc
            for sidx in range(1, ki + 1):
                P = S[ki - sidx]
                if P.is_zero(): continue
                P, D = self.cancel_letters(P, dict(DE))       # exact cancellation of denominator letters (trial division)
                poles[(i, sidx)] = {(self.dkey(D), ()): P * self.qq(inv)}
            self.stats.setdefault("t_cancel", 0.0); self.stats["t_cancel"] += _t.time() - _t0
        if degN >= K:
            # polynomial part: reuse pf's expansion at infinity on TS-wrapped coefficients (small: only degN-K+1 orders)
            Ncoef = {n: {((), ()): P} for n, P in Npoly.items()}
            pp, _unused = self.pf(Ncoef, dens, var, poly_only=True)
            polypart = pp
        return polypart, poles

    # ------------------------------------------------------------ expansions of G(w;x)
    def expand_inf(self, w, M):
        """{(j, i): TS} : G(w;x) = sum TS_{j,i} log^j x  x^-i + O(x^-(M+1)), constants ZIP[subwords] as words."""
        key = ("inf", w, M)
        if key in self._powcache: return self._powcache[key]
        if not w:
            return {(0, 0): self.ts_one()}
        a = w[0]; Dp = self.expand_inf(w[1:], M)
        out = {(0, 0): self.ts_attach(self.ts_one(), [w])}
        ca = self.hl_c(a); Ea = self.hl_E(a)
        for (j, i), T in Dp.items():
            kmax = 0 if a == ZERO_HL else M - i   # need n = i + k <= M  (n>=1 terms) ; n = 0 only if i = k = 0
            for k in range(0, max(kmax, 0) + 1):
                n = i + k
                if n > M: break
                Tk = T if k == 0 else self.ts_mul_mono(T, ca ** k, {l: e * k for l, e in Ea.items()})
                if n == 0:
                    self.ts_add_into(out.setdefault((j + 1, 0), {}), Tk, self.qq(Fr(1, j + 1)))
                else:
                    for l in range(j + 1):
                        coef = -Fr(math.factorial(j), math.factorial(j - l)) / Fr(n) ** (l + 1)
                        self.ts_add_into(out.setdefault((j - l, n), {}), Tk, self.qq(coef))
                if a == ZERO_HL: break
        out = {k: v for k, v in out.items() if v}
        self._powcache[key] = out
        return out
    def expand_zero(self, w, M):
        """{(j, n): TS}: G(w;x) = sum TS log^j x x^n + O(x^(M+1)) at x -> 0 (G(0;x) = log x)."""
        key = ("zero", w, M)
        if key in self._powcache: return self._powcache[key]
        if not w:
            return {(0, 0): self.ts_one()}
        a = w[0]; Dp = self.expand_zero(w[1:], M)
        out = {}
        if a == ZERO_HL:
            # int_0^x dt/t  t^n log^j t
            for (j, n), T in Dp.items():
                if n == 0:
                    self.ts_add_into(out.setdefault((j + 1, 0), {}), T, self.qq(Fr(1, j + 1)))
                else:
                    for l in range(j + 1):
                        coef = (-1) ** l * Fr(math.factorial(j), math.factorial(j - l)) / Fr(n) ** (l + 1)
                        self.ts_add_into(out.setdefault((j - l, n), {}), T, self.qq(coef))
        else:
            # 1/(t-a) = -(1/a) sum_k (t/a)^k
            ca = self.hl_c(a); Ea = self.hl_E(a)
            for (j, n0), T in Dp.items():
                for k in range(0, M + 1):
                    n = n0 + k + 1          # after integration t^(n0+k) -> x^(n0+k+1)
                    if n > M: break
                    Tk = self.ts_mul_mono(T, -Fr(1) / ca ** (k + 1), {l: -e * (k + 1) for l, e in Ea.items()})
                    for l in range(j + 1):
                        coef = (-1) ** l * Fr(math.factorial(j), math.factorial(j - l)) / Fr(n) ** (l + 1)
                        self.ts_add_into(out.setdefault((j - l, n), {}), Tk, self.qq(coef))
        out = {k: v for k, v in out.items() if v}
        self._powcache[key] = out
        return out

    # ------------------------------------------------------------ integrals
    def J(self, k, s, w):
        """reg int_0^oo (x-s)^-k G(w;x) dx -> TS (coefficient 1 context)"""
        key = ("J", k, s, w)
        if key in self._powcache: return self._powcache[key]
        if k == 1:
            res = self.ts_attach(self.ts_one(), [(s,) + tuple(w)])
        elif not w:
            if s == ZERO_HL: res = {}
            else:
                cs = self.hl_c(s); Es = self.hl_E(s)
                res = self.ts_mono((-cs) ** (1 - k) / (k - 1), {l: e * (1 - k) for l, e in Es.items()})
        else:
            a = w[0]; wp = tuple(w[1:])
            res = {}
            # boundary at 0 (only when s == 0): (1/(k-1)) [x^(k-1) log^0] of G(w;x)
            if s == ZERO_HL:
                ez = self.expand_zero(tuple(w), k - 1)
                T = ez.get((0, k - 1))
                if T: self.ts_add_into(res, T, self.qq(Fr(1, k - 1)))
            n = k - 1
            if a == s:
                self.ts_add_into(res, self.J(k, s, wp), self.qq(Fr(1, k - 1)))
            else:
                cd, Ed = self.hl_diff_mono(a, s)          # d = a - s
                # 1/((x-s)^n (x-a)) = d^-n/(x-a) - sum_{j=1}^{n} d^-(n+1-j) / (x-s)^j
                t = self.ts_mul_mono(self.J(1, a, wp), Fr(1) / cd ** n, {l: -e * n for l, e in Ed.items()})
                self.ts_add_into(res, t, self.qq(Fr(1, k - 1)))
                for j in range(1, n + 1):
                    pw = n + 1 - j
                    t = self.ts_mul_mono(self.J(j, s, wp), -Fr(1) / cd ** pw, {l: -e * pw for l, e in Ed.items()})
                    self.ts_add_into(res, t, self.qq(Fr(1, k - 1)))
        self._powcache[key] = res
        return res
    def I(self, m, w, M=None):
        """reg int_0^oo x^m G(w;x) dx -> TS"""
        key = ("I", m, w)
        if key in self._powcache: return self._powcache[key]
        if not w:
            return {}
        a = w[0]; wp = tuple(w[1:])
        Mx = M if M is not None else m + 1
        E = self.expand_inf(tuple(w), max(Mx, m + 1))
        res = {}
        T = E.get((0, m + 1))
        if T: self.ts_add_into(res, T, self.qq(Fr(1, m + 1)))
        if a == ZERO_HL:
            self.ts_add_into(res, self.I(m, wp, Mx), self.qq(-Fr(1, m + 1)))
        else:
            ca = self.hl_c(a); Ea = self.hl_E(a)
            for i in range(m + 1):
                Ii = self.I(i, wp, Mx)
                if Ii:
                    t = self.ts_mul_mono(Ii, ca ** (m - i), {l: e * (m - i) for l, e in Ea.items()})
                    self.ts_add_into(res, t, self.qq(-Fr(1, m + 1)))
            t = self.ts_mul_mono(self.J(1, a, wp), ca ** (m + 1), {l: e * (m + 1) for l, e in Ea.items()})
            self.ts_add_into(res, t, self.qq(-Fr(1, m + 1)))
        self._powcache[key] = res
        return res

    # ------------------------------------------------------------ activation of x-dependent weight-1 words
    def activate_word(self, word, var):
        """ZIP[word] with var-dependent letters -> list of (coef Fr, const_words tuple, active G-word tuple).
        weight 1 only: ZIP[s(x)] = ZIP[s(0)] - sum_i e_i G(rho_i; x)."""
        if not any(self.hl_depends(h, var) for h in word):
            return [(Fr(1), (word,), ())]
        if len(word) != 1:
            raise NotImplementedError("activation of x-dependent words of weight >= 2 (needs the symbol-level fibration change): %s" % (word,))
        h = word[0]; c = self.hl_c(h); E = self.hl_E(h)
        out = []
        c0 = c; E0 = {}
        for i, e in E.items():
            if self.A.depends(i, var):
                alpha, beta = self.A.linear_parts(i, var)
                rho = self.hl_from_alpha_beta(alpha, beta, "act:%s|%s" % (var, self.A.describe(i)[:40]))
                out.append((Fr(-e), (), (rho,)))
                # L_i = alpha (x - rho): L_i(0) = beta  -> constant part gets beta^e
                if beta.is_zero():
                    raise NotImplementedError("activation with a letter vanishing at x=0 (log x piece)")
                cb, Eb = self.factor(beta, "beta:act")
                c0 *= cb ** e
                for l, ee in Eb.items(): E0[l] = E0.get(l, 0) + ee * e
            else:
                E0[i] = E0.get(i, 0) + e
        s0 = self.hl(c0, E0)
        out.insert(0, (Fr(1), ((s0,),), ()))
        return out

    # ------------------------------------------------------------ one elimination
    def integrate(self, ts, var, log=None, consolidate=True, tag=""):
        """integrate every term of ts over var in (0,oo); returns TS free of var."""
        vj = self.A.names.index(var)
        out = {}
        nterm = 0; t0 = time.time()
        for (Dk, Wk), N in ts.items():
            nterm += 1
            # split denominators
            dens = {}; D0 = {}
            for i, e in Dk:
                (dens if self.A.depends(i, var) else D0)[i] = e
            # activate words
            combos = [(Fr(1), (), ())]
            for word in Wk:
                opts = self.activate_word(word, var)
                combos = [(c1 * c2, cw1 + cw2, self._shuffle_pair(aw1, aw2)) for (c1, cw1, aw1) in combos for (c2, cw2, aw2) in opts]
            # flatten shuffles: aw may be dict word->mult after shuffle
            flat = []
            for c, cw, aw in combos:
                if isinstance(aw, dict):
                    for wd, m in aw.items(): flat.append((c * m, cw, wd))
                else:
                    flat.append((c, cw, aw))
            # numerator as polynomial in var with TS coefficients (pure numerators)
            # cheaper: build each N_n in one go
            byn = {}
            for e, cf in N.to_dict().items():
                n = int(e[vj]); e2 = list(int(t) for t in e); e2[vj] = 0
                byn.setdefault(n, {})[tuple(e2)] = cf
            Npoly = {n: self.ctx.from_dict(d) for n, d in byn.items()}
            Ncoef = {n: {((), ()): P} for n, P in Npoly.items()}
            if not dens:
                if any(aw for _, _, aw in flat):
                    polypart, poles = {n: T for n, T in Ncoef.items()}, {}
                else:
                    continue   # pure power divergence int x^n dx = 0  (no x-dependence at all in denominators and words)
            elif self.use_pf2:
                polypart, poles = self.pf2(Npoly, dens, var)
            else:
                polypart, poles = self.pf(Ncoef, dens, var)
            for c, cw, aw in flat:
                acc = {}
                for p, T in polypart.items():
                    if not aw: continue
                    Ip = self.I(p, aw)
                    if Ip: self.ts_add_into(acc, self._ts_mul_ts(T, Ip))
                for (i, s), T in poles.items():
                    alpha, beta = self.A.linear_parts(i, var)
                    sig = self.hl_from_alpha_beta(alpha, beta, "pole:%s|%s" % (var, self.A.describe(i)[:40]))
                    ca, Ea = self.factor(alpha, "alpha:int")
                    Js = self.J(s, sig, aw)
                    if not Js: continue
                    Tp = self.ts_mul_mono(T, Fr(1) / ca ** s, {l: -e * s for l, e in Ea.items()})
                    self.ts_add_into(acc, self._ts_mul_ts(Tp, Js))
                if not acc: continue
                acc = self.ts_mul_mono(acc, c, {l: -e for l, e in D0.items()})
                if cw: acc = self.ts_attach(acc, list(cw))
                self.ts_add_into(out, acc)
            if log and nterm % 1 == 0:
                log("  [%s] term %d/%d  out keys %d  monomials %d  %.1fs" % (tag, nterm, len(ts), len(out), self.ts_nmon(out), time.time() - t0))
        if consolidate:
            out = self.consolidate(out, log=log)
        return out
    @staticmethod
    def _shuffle_pair(u, v):
        from itertools import product
        if isinstance(u, dict) or isinstance(v, dict):
            du = u if isinstance(u, dict) else {u: 1}; dv = v if isinstance(v, dict) else {v: 1}
            out = {}
            for a, ma in du.items():
                for b, mb in dv.items():
                    for wd, m in Engine._sh(a, b).items(): out[wd] = out.get(wd, 0) + m * ma * mb
            return out
        if not u: return v
        if not v: return u
        return Engine._sh(u, v)
    @staticmethod
    def _sh(u, v):
        if not u: return {v: 1}
        if not v: return {u: 1}
        out = {}
        for w, c in Engine._sh(u[1:], v).items(): out[(u[0],) + w] = out.get((u[0],) + w, 0) + c
        for w, c in Engine._sh(u, v[1:]).items(): out[(v[0],) + w] = out.get((v[0],) + w, 0) + c
        return out
    def _ts_mul_ts(self, T, U):
        out = {}
        for (Dk, Wk), N0 in T.items():
            for (Dk2, Wk2), N1 in U.items():
                D = dict(Dk)
                for l, e in Dk2: D[l] = D.get(l, 0) + e
                k = (self.dkey(D), tuple(sorted(Wk + Wk2)))
                v = out.get(k); p = self.red(N0 * N1)
                out[k] = p if v is None else v + p
        return {k: v for k, v in out.items() if not v.is_zero()}

    # ------------------------------------------------------------ consolidation (same letter support only) + factor cancellation
    def consolidate(self, ts, log=None, cancel=True):
        groups = {}
        for (Dk, Wk), N in ts.items():
            sup = tuple(sorted(i for i, e in Dk))
            groups.setdefault((sup, Wk), []).append((dict(Dk), N))
        out = {}
        for (sup, Wk), lst in groups.items():
            Dmax = {}
            for D, N in lst:
                for i, e in D.items(): Dmax[i] = max(Dmax.get(i, 0), e)
            acc = None
            for D, N in lst:
                P = N
                for i, e in Dmax.items():
                    r = e - D.get(i, 0)
                    if r: P = self.red(P * self.lpow(i, r))
                acc = P if acc is None else acc + P
            if acc.is_zero(): continue
            D = dict(Dmax)
            if cancel:
                acc, D = self.cancel_letters(acc, D)
            k = (self.dkey(D), Wk)
            v = out.get(k)
            out[k] = acc if v is None else v + acc
        out = {k: v for k, v in out.items() if not v.is_zero()}
        if log: log("  consolidate: %d -> %d keys, monomials %d -> %d" % (len(ts), len(out), self.ts_nmon(ts), self.ts_nmon(out)))
        return out

    # ------------------------------------------------------------ output canonical form: one fraction per word product
    def per_word_form(self, ts, log=None):
        """merge all keys with the same word product over the lcm of their (letter-monomial) denominators + exact letter cancellation.
        Used ONLY as the final output normal form (after the last elimination of a run), never between eliminations."""
        byw = {}
        for (Dk, Wk), N in ts.items(): byw.setdefault(Wk, []).append((dict(Dk), N))
        out = {}; dropped = 0
        for Wk, lst in byw.items():
            Dmax = {}
            for D, N in lst:
                for i, e in D.items(): Dmax[i] = max(Dmax.get(i, 0), e)
            acc = None
            for D, N in lst:
                P = N
                for i, e in Dmax.items():
                    r = e - D.get(i, 0)
                    if r: P = self.red(P * self.lpow(i, r))
                acc = P if acc is None else acc + P
            if acc.is_zero(): dropped += 1; continue
            P, D = self.cancel_letters(acc, dict(Dmax))
            out[(self.dkey(D), Wk)] = P
        if log: log("  per-word form: %d keys/%d monomials -> %d words/%d monomials (%d word products cancel exactly)" % (len(ts), self.ts_nmon(ts), len(out), self.ts_nmon(out), dropped))
        return out, dropped

    # ------------------------------------------------------------ census
    def census(self, ts):
        words = set(); letters = set(); dens = set(); maxw = 0
        for (Dk, Wk), N in ts.items():
            for i, e in Dk: dens.add(i)
            for w in Wk:
                words.add(w); maxw = max(maxw, len(w))
                for h in w: letters.add(h)
        return {"n_keys": len(ts), "n_monomials": self.ts_nmon(ts), "n_distinct_words": len(words), "n_distinct_hlog_letters": len(letters),
                "max_word_weight": maxw, "denominator_letters": sorted(dens), "max_numerator_terms": max((len(N) for N in ts.values()), default=0)}
