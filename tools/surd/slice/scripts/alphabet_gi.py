"""GI alphabet: auto-extending letter table over the ring R = Q[x..., t, j]/(j^2+1) (j = sqrt(-1)); letters are irreducible over Q(i).
A polynomial of j-degree <= 1 is factored as  c * prod letters^e  (c in Q; an overall power of the unit j is absorbed by multiplying one
odd-exponent letter by j, or, if all exponents are even, by returning c with a 'unit' flag -> we then multiply the FIRST letter's exponent... see factor()).
Q-irreducible j-free polynomials are split over Q(i) with sympy.factor_list(gaussian=True) (cached, small polynomials only: alphabet letters).
Products of letters are only meaningful modulo j^2+1: use red()."""
import os, sys, json
from fractions import Fraction as Fr
from flint import fmpq, fmpq_mpoly_ctx
import sympy as sp
import alphabet
from alphabet import normalize, OffAlphabet

class GIAlphabet(alphabet.Alphabet):
    def __init__(self, ctx, jname="j", cache_path=None):
        self.jname = jname
        super().__init__(ctx)
        self.nv = len(self.names); self.jj = self.names.index(jname)
        self.J = ctx.gens()[self.jj]; self.J2 = self.J * self.J + ctx.from_dict({(0,) * self.nv: fmpq(1)})
        self.one = ctx.from_dict({(0,) * self.nv: fmpq(1)})
        self.auto_added = []; self._gsplit = {}; self._conj = {}; self.cache_path = cache_path; self._gcache = {}
        if cache_path and os.path.exists(cache_path):
            try: self._gcache = json.load(open(cache_path))
            except Exception: self._gcache = {}
        self.stats = {"gaussian_factor_calls": 0, "gaussian_split": 0, "gaussian_inert": 0}
    # ---------------------------------------------------------------- ring helpers
    def red(self, P):
        if int(P.degrees()[self.jj]) < 2: return P
        return divmod(P, self.J2)[1]
    def jdeg(self, P): return int(P.degrees()[self.jj])
    def conj(self, P):
        """j -> -j (P reduced, deg_j <= 1)"""
        B = P.derivative(self.jj); return P - (self.J + self.J) * B
    def parts(self, P):
        """P = A + j B (deg_j <= 1) -> (A, B) j-free (C-level: B = dP/dj, A = P - j B)"""
        assert self.jdeg(P) <= 1
        B = P.derivative(self.jj); A = P - self.J * B
        return A, B
    def norm(self, P):
        A, B = self.parts(P); return A * A + B * B
    def _spec_images(self):
        """random integer specialization of every variable except j and the 'main' variable (chosen per letter): cached images per main var"""
        if not hasattr(self, "_spec"):
            import random
            rng = random.Random(20260921); self._spec_vals = [fmpq(rng.randint(7, 97)) for _ in self.names]; self._spec = {}
        return self._spec_vals
    def _specialize(self, P, main):
        """P -> (A, B) fmpq_poly in the main variable after substituting random integers for the other non-j variables (P = A + j B)"""
        from flint import fmpq_poly
        vals = self._spec_images()
        C2 = fmpq_mpoly_ctx.get((self.names[main], self.jname)); gm, gj = C2.gens()
        imgs = []
        for k, n in enumerate(self.names):
            if k == main: imgs.append(gm)
            elif k == self.jj: imgs.append(gj)
            else: imgs.append(C2.from_dict({(0, 0): vals[k]}))
        P2 = P.compose(*imgs, ctx=C2); d = P2.to_dict(); A = {}; B = {}
        for e, c in d.items(): (A if int(e[1]) == 0 else B)[int(e[0])] = c
        def up(dd):
            if not dd: return fmpq_poly([])
            m = max(dd); return fmpq_poly([dd.get(i, fmpq(0)) for i in range(m + 1)])
        return up(A), up(B)
    def _pretest(self, P, L, li=None):
        """necessary condition for L | P in R via a univariate specialization (main variable = first variable L depends on); True = may divide"""
        degs = [int(z) for z in L.degrees()]
        mains = [k for k in range(self.nv) if k != self.jj and degs[k] > 0]
        if not mains: return True
        main = mains[0]
        key = ("specL", str(L))
        if key not in self._fact_cache: self._fact_cache[key] = self._specialize(L, main)
        a, b = self._fact_cache[key]; A, B = self._specialize(P, main)
        n = a * a + b * b
        if n.degree() <= 0: return True
        ra = divmod(A * a + B * b, n)[1]; rb = divmod(B * a - A * b, n)[1]
        return ra.degree() < 0 and rb.degree() < 0
    def _lbn(self, L, k):
        """(conj(L)^k reduced, norm(L)^k) cached by letter string"""
        key = ("lbn", str(L), k)
        if key not in self._fact_cache:
            if k == 1:
                Lb = self.conj(L); self._fact_cache[key] = (Lb, self.red(L * Lb))
            else:
                h = k // 2; Lb1, n1 = self._lbn(L, h); Lb2, n2 = self._lbn(L, k - h)
                self._fact_cache[key] = (self.red(Lb1 * Lb2), n1 * n2)
        return self._fact_cache[key]
    def gi_divides(self, P, L, k=1):
        """exact division in R: returns Q with P == Q*L^k mod j^2+1, or None.  P, L reduced (deg_j <= 1)."""
        self.stats["gi_div_calls"] = self.stats.get("gi_div_calls", 0) + 1
        if self.jdeg(L) == 0:
            q, r = divmod(P, L if k == 1 else L ** k)
            return q if r.is_zero() else None
        if k == 1 and not self._pretest(P, L):
            self.stats["gi_div_pretest_reject"] = self.stats.get("gi_div_pretest_reject", 0) + 1
            return None
        Lb, n = self._lbn(L, k)
        num = self.red(P * Lb); A, B = self.parts(num)
        qa, ra = divmod(A, n)
        if not ra.is_zero(): return None
        qb, rb = divmod(B, n)
        if not rb.is_zero(): return None
        return qa + self.J * qb
    def gi_divide_maxpow(self, P, L, emax):
        """largest k <= emax with L^k | P in R, galloping; returns (k, P / L^k)"""
        if emax <= 0: return 0, P
        q = self.gi_divides(P, L, 1)
        if q is None: return 0, P
        k = 1; P = self.red(q); step = 2
        while k < emax:
            s = min(step, emax - k)
            q = self.gi_divides(P, L, s)
            if q is None:
                if s == 1: break
                step = max(1, s // 2)
                # retry smaller chunk
                while step >= 1:
                    s = min(step, emax - k); q = self.gi_divides(P, L, s)
                    if q is not None: k += s; P = self.red(q); break
                    if s == 1: return k, P
                    step //= 2
                else: break
                continue
            k += s; P = self.red(q); step *= 2
        return k, P
    # ---------------------------------------------------------------- gaussian split of a Q-irreducible j-free polynomial
    def _to_sympy(self, P):
        syms = [sp.Symbol(n) for n in self.names]; expr = 0
        for e, c in P.to_dict().items():
            term = sp.Rational(int(c.p), int(c.q))
            for k, ee in enumerate(e):
                if int(ee): term = term * syms[k] ** int(ee)
            expr = expr + term
        return expr, syms
    def _from_sympy(self, expr, syms):
        expr = sp.expand(expr).subs(sp.I, syms[self.jj])
        P = sp.Poly(sp.expand(expr), *syms)
        return self.ctx.from_dict({tuple(int(z) for z in m): fmpq(int(sp.Rational(c).p), int(sp.Rational(c).q)) for m, c in P.as_dict().items()})
    def gaussian_split(self, f):
        """f: Q-irreducible, j-free, normalized.  -> None if irreducible over Q(i), else g (normalized, deg_j = 1) with f = c g conj(g) mod j^2+1."""
        k = str(f)
        if k in self._gsplit: return self._gsplit[k]
        # quick necessary condition: even degree in every variable
        if any(int(d) % 2 for d in f.degrees()) or f.total_degree() == 0:
            self._gsplit[k] = None; return None
        if k in self._gcache:
            g = None if self._gcache[k] is None else self.ctx.from_dict({tuple(e): fmpq(p, q) for e, p, q in self._gcache[k]})
        else:
            self.stats["gaussian_factor_calls"] += 1
            expr, syms = self._to_sympy(f)
            c, facs = sp.factor_list(expr, gaussian=True)
            cplx = [(ff, m) for ff, m in facs if ff.has(sp.I)]
            if not cplx: g = None
            else:
                assert len(cplx) == 2 and all(m == 1 for ff, m in cplx), ("unexpected gaussian factorization", facs)
                g = self._from_sympy(cplx[0][0], syms); g = normalize(g)[1]
                chk = self.red(g * self.conj(g)); cc, E = super().factor(chk, "gsplit-check") if False else (None, None)
                q, r = divmod(chk, f)
                assert r.is_zero() and q.total_degree() == 0, "gaussian split check failed"
            self._gcache[k] = None if g is None else [(list(int(z) for z in e), int(c.p), int(c.q)) for e, c in g.to_dict().items()]
            if self.cache_path:
                tmp = self.cache_path + ".tmp%d" % os.getpid(); json.dump(self._gcache, open(tmp, "w")); os.replace(tmp, self.cache_path)
        self._gsplit[k] = g
        self.stats["gaussian_split" if g is not None else "gaussian_inert"] += 1
        return g
    # ---------------------------------------------------------------- letters
    def add(self, p, source):
        idx = super().add(p, source)
        return idx
    def _addq(self, q, context):
        kk = str(q)
        if kk not in self.key2idx:
            idx = self.add(q, "auto:" + context.split(":")[0]); self.auto_added.append((idx, kk, context[:60]))
        return self.key2idx[kk]
    def factor(self, p, context=""):
        """p (deg_j <= 1 after red) -> (c, {idx: e}) with p == c * j^u * prod letters^e mod (j^2+1); the unit j^u is folded in: letter 'j' itself is
        letter index self.jj (a variable letter), so u is returned as an exponent of that letter (j^2 = -1 handled by c)."""
        p = self.red(p); k = str(p)
        if k in self._fact_cache: return self._fact_cache[k]
        if p.is_zero(): raise ValueError("factor(0) " + context)
        c = Fr(1); E = {}
        def addE(i, e):
            E[i] = E.get(i, 0) + e
        if self.jdeg(p) == 0:
            cont, facs = p.factor(); c *= Fr(int(cont.p), int(cont.q))
            for f, e in facs:
                cf, q = normalize(f); c *= cf ** int(e)
                if f.total_degree() == 0: continue
                g = self.gaussian_split(q)
                if g is None: addE(self._addq(q, context), int(e))
                else:
                    gb = normalize(self.conj(g))[1]
                    # q = cq * g * gb mod j^2+1 : determine cq
                    prod = self.red(g * gb); qq, r = divmod(prod, q); assert r.is_zero() and qq.total_degree() == 0
                    cq = Fr(1) / Fr(int(qq.to_dict()[(0,) * self.nv].p), int(qq.to_dict()[(0,) * self.nv].q))
                    c *= cq ** int(e); addE(self._addq(g, context), int(e)); addE(self._addq(gb, context), int(e))
        else:
            A, B = self.parts(p)
            # j-free content
            gAB = A.gcd(B)
            if gAB.total_degree() > 0:
                cg, Eg = self.factor(gAB, context)      # recursion on a j-free polynomial
                c *= cg
                for i, e in Eg.items(): addE(i, e)
                A = divmod(A, gAB)[0]; B = divmod(B, gAB)[0]; p1 = A + self.J * B
            else:
                p1 = p
            if p1.total_degree() == 1 and self.parts(p1)[0].is_zero():
                # p1 = b*j, b constant
                bb = self.parts(p1)[1].to_dict()[(0,) * self.nv]; c *= Fr(int(bb.p), int(bb.q)); addE(self.jj, 1)
            elif self.jdeg(p1) == 0:
                cf, q = normalize(p1); c *= cf
                if q.total_degree() > 0:
                    c2, E2 = self.factor(q, context); c *= c2
                    for i, e in E2.items(): addE(i, e)
            else:
                n = self.red(p1 * self.conj(p1))          # j-free norm; all its Q-irreducible factors split (gcd(A,B)=1) except possibly constants
                cont, facs = n.factor(); rem = p1
                for f, e in facs:
                    if f.total_degree() == 0: continue
                    q = normalize(f)[1]; g = self.gaussian_split(q)
                    if g is None:
                        # inert prime dividing the norm of a primitive element: must divide p1 to an even power... cannot happen with gcd(A,B)=1 unless q | both
                        raise OffAlphabet("inert factor %s in norm of primitive %s (%s)" % (str(q)[:80], str(p1)[:80], context))
                    gb = normalize(self.conj(g))[1]
                    etot = int(e)      # multiplicity of q in the norm = mult of g in p1 + mult of gb in p1 (q = norm of g up to const, g != gb up to units when q is Q-irreducible of even degree... )
                    for cand in (g, gb):
                        while etot > 0:
                            qt = self.gi_divides(rem, cand)
                            if qt is None: break
                            rem = self.red(qt); addE(self._addq(cand, context), 1); etot -= 1
                    assert etot == 0, ("norm multiplicity not exhausted", str(q)[:60], context)
                # rem is now a unit of R: a + b j constant
                A_, B_ = self.parts(rem); assert A_.total_degree() <= 0 and B_.total_degree() <= 0, ("non-constant remainder", str(rem)[:80])
                A0, B0 = self.parts(rem); a0 = A0.to_dict().get((0,) * self.nv, fmpq(0)); b0 = B0.to_dict().get((0,) * self.nv, fmpq(0))
                if b0 == 0: c *= Fr(int(a0.p), int(a0.q))
                elif a0 == 0: c *= Fr(int(b0.p), int(b0.q)); addE(self.jj, 1)
                else:
                    cq, u, primes = self.gauss_const_factor(Fr(int(a0.p), int(a0.q)), Fr(int(b0.p), int(b0.q)))
                    c *= cq
                    if u: addE(self.jj, u)
                    for (pa, pb), e in primes.items():
                        cf, qn = normalize(self.one * fmpq(pa) + self.J * fmpq(pb)); c *= cf ** e
                        addE(self._addq(qn, "gaussprime"), e)
        # fold j^2 = -1
        if E.get(self.jj, 0):
            u = E[self.jj]; c *= Fr(-1) ** (u // 2); u = u % 2
            if u: E[self.jj] = u
            else: del E[self.jj]
        res = (c, {i: e for i, e in E.items() if e})
        self._fact_cache[k] = res
        return res
    def gauss_const_factor(self, a, b):
        """a + b j (rationals, both nonzero) = cq * j^u * prod pi^e with canonical Gaussian primes pi = (pa, pb): for p = 2: (1,1); for p = 1 mod 4:
        (x, y) and (x, -y) with x > y > 0, x^2 + y^2 = p.  Returns (cq Fraction, u in 0..3, {(pa,pb): e})."""
        from math import gcd
        import sympy as _sp
        den = a.denominator * b.denominator // gcd(a.denominator, b.denominator)
        A = int(a * den); B = int(b * den); g = gcd(abs(A), abs(B)); A //= g; B //= g
        cq = Fr(g, den); primes = {}
        N = A * A + B * B
        za, zb = A, B          # current Gaussian integer
        def gdiv(xa, xb, pa, pb):
            # (xa + xb i)/(pa + pb i) exact?  = (xa+xb i)(pa - pb i)/n
            n = pa * pa + pb * pb; ra = xa * pa + xb * pb; rb = xb * pa - xa * pb
            if ra % n == 0 and rb % n == 0: return ra // n, rb // n
            return None
        for p, e in _sp.factorint(N).items():
            if p == 2: cands = [(1, 1)]
            elif p % 4 == 1:
                from sympy.solvers.diophantine.diophantine import cornacchia
                sols = cornacchia(1, 1, p); (x, y) = max((max(s_), min(s_)) for s_ in sols)
                cands = [(x, y), (x, -y)]
            else: raise ValueError("prime 3 mod 4 divides a primitive Gaussian integer norm?")
            left = e
            for (pa, pb) in cands:
                while left > 0:
                    r = gdiv(za, zb, pa, pb)
                    if r is None: break
                    za, zb = r; primes[(pa, pb)] = primes.get((pa, pb), 0) + 1; left -= (2 if p == 2 else 1) if False else 1
                    if p == 2: left -= 0
            # for p == 2 the norm exponent e counts (1+i) twice per... N((1+i)) = 2 so each division removes one factor of 2 from N: fine
            assert left == 0, ("gauss prime exhaustion", p, e, left)
        # remaining unit
        units = {(1, 0): 0, (0, 1): 1, (-1, 0): 2, (0, -1): 3}
        assert (za, zb) in units, ("not a unit", za, zb)
        return cq, units[(za, zb)], primes
    def check_factor(self, p, c, E):
        """verify p == c * prod letters^E mod j^2+1 (debug)"""
        q = self.one * fmpq(c.numerator, c.denominator)
        for i, e in E.items():
            for _ in range(e): q = self.red(q * self.letters[i])
        return self.red(p - q).is_zero()
