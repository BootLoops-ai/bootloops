#!/usr/bin/env python3
# galois member — exact MZV relation engine (depth<=3, weight<=21; NO PSLQ).
r"""mzvring.py — exact MZV relation engine, depth<=3, weight<=21. NO PSLQ.

Symbols: ('pi2',) [= pi^2]; ('z',n) n odd>=3; ('m2',a,b); ('m3',a,b,c)
(outer-first convention: m2(a,b)=sum_{m>n} m^-a n^-b, a>=2).
Monomial = sorted tuple of symbols, () = 1.  Vector = {monomial: Fraction}.
zeta(2k) is normalized to rational*pi2^k at construction (Bernoulli).

Relation families, generated per weight w and Gauss-eliminated over Q with
single m2/m3 symbols as the only eliminable columns (products/pi2-towers/
z(odd) are permanent basis monomials):
  R2 stuffle  z(a)z(b) = m2(a,b)+m2(b,a)+z(w)
  R3 shuffle  z(a)z(b) = sum_shuffle m2(comps)           [x^{a-1}y sh x^{b-1}y]
  R4 stuffle  z(a)m2(b,c) = m3(a,b,c)+m3(b,a,c)+m3(b,c,a)+m2(a+b,c)+m2(b,a+c)
  R5 shuffle  z(a)m2(b,c) = sum_shuffle m3(comps)
  R6 Hoffman (T^0 regularized double shuffle, one y): for convergent u of
     weight w-1, depth 1|2:  [y sh word(u)] - [1 * u] = 0 after the single
     divergent word y.word(u) cancels between the two expansions.
Elimination order prefers KEEPING: products, z(w), all-odd no-1 m2/m3 with
lexicographically largest args (so e.g. m2(5,3) survives at w=8).
Lower-weight symbols inside product monomials are canonicalized recursively
(tables built weight-ascending by build_all)."""
from fractions import Fraction as Fr
from functools import lru_cache
import sympy as sp

PI2 = ('pi2',)

def mmul(m1, m2):
    return tuple(sorted(m1+m2))

def dadd(d, m, c):
    if c:
        d[m] = d.get(m, Fr(0)) + c
        if d[m] == 0:
            del d[m]

def dsum(*ds):
    o = {}
    for d in ds:
        for m, c in d.items():
            dadd(o, m, c)
    return o

def dscale(d, f):
    return {} if f == 0 else {m: c*f for m, c in d.items()}

def dmul(d1, d2):
    o = {}
    for m1, c1 in d1.items():
        for m2, c2 in d2.items():
            dadd(o, mmul(m1, m2), c1*c2)
    return o

@lru_cache(maxsize=None)
def z2rat(k):
    """zeta(2k)/pi^{2k} as Fraction."""
    B = sp.Rational(sp.bernoulli(2*k))
    return Fr(int(B.p), int(B.q))*Fr((-1)**(k+1)*2**(2*k-1), int(sp.factorial(2*k)))

def znorm(n):
    """zeta(n) as a vector (even n -> rational * pi2^k)."""
    if n % 2 == 0:
        return {(PI2,)*(n//2): z2rat(n//2)}
    return {(('z', n),): Fr(1)}

def sweight(s):
    return 2 if s[0] == 'pi2' else sum(s[1:])

def mweight(m):
    return sum(sweight(s) for s in m)

# ---- words: tuple of 0(=x)/1(=y); comp (s1..sk) -> x^{s1-1}y ... ----
def comp_word(c):
    w = ()
    for s in c:
        w += (0,)*(s-1)+(1,)
    return w

def word_comp(w):
    out, run = [], 0
    for L in w:
        if L == 0:
            run += 1
        else:
            out.append(run+1)
            run = 0
    assert run == 0
    return tuple(out)

@lru_cache(maxsize=None)
def shuffle(w1, w2):
    if not w1:
        return {w2: 1}
    if not w2:
        return {w1: 1}
    o = {}
    for w, c in shuffle(w1[1:], w2).items():
        o[(w1[0],)+w] = o.get((w1[0],)+w, 0)+c
    for w, c in shuffle(w1, w2[1:]).items():
        o[(w2[0],)+w] = o.get((w2[0],)+w, 0)+c
    return o

def comp_sym(c):
    assert len(c) in (2, 3), f"comp depth {len(c)} not a symbol: {c}"
    return ('m2',)+c if len(c) == 2 else ('m3',)+c

def stuffle1(u):
    """index-terms of 1 * (s1..sk), divergent leading-1 term EXCLUDED."""
    k = len(u)
    terms = []
    for i in range(1, k+1):                      # insert 1 after position i
        terms.append(u[:i]+(1,)+u[i:])
    for i in range(k):                           # merge into s_i
        terms.append(u[:i]+(u[i]+1,)+u[i+1:])
    return terms

def elim_key(s):
    """Eliminate smaller keys first; survivors = largest keys."""
    args = s[1:]
    depth_rank = 0 if s[0] == 'm3' else 1        # kill m3 before m2
    no_even = 0 if any(a % 2 == 0 for a in args) else 1
    no_one = 0 if any(a == 1 for a in args) else 1
    return (depth_rank, no_even, no_one, args)


class Ring:
    """Canonicalization tables, built weight-ascending."""

    def __init__(self, wmax=21):
        self.canon = {}          # sym -> vector over canonical monomials
        self.survivors = {}      # w -> [syms kept]
        self.wmax = wmax
        for w in range(3, wmax+1):
            self._build(w)

    # -- canonical form of symbols / monomials / vectors (lower weights ready)
    def csym(self, s):
        if s in self.canon:
            return self.canon[s]
        if s[0] == 'z' and s[1] % 2 == 0:      # zeta(2k) -> rational*pi2^k
            return znorm(s[1])
        return {(s,): Fr(1)}

    def cmono(self, m):
        v = {(): Fr(1)}
        for s in m:
            v = dmul(v, self.csym(s))
        return v

    def cvec(self, d):
        return dsum(*(dscale(self.cmono(m), c) for m, c in d.items())) if d else {}

    def _rows(self, w):
        """[(U:{sym:Fr}, R:vector)] meaning sum U[s]*s = R."""
        rows = []
        for a in range(2, w//2+1):               # R2,R3: z(a)z(b), a<=b
            b = w-a
            prod = self.cvec(dmul(znorm(a), znorm(b)))
            U = {}
            dadd(U, ('m2', a, b), Fr(1)); dadd(U, ('m2', b, a), Fr(1))
            rows.append((U, dsum(prod, dscale(znorm(w), Fr(-1)))))
            U = {}
            for wd, c in shuffle(comp_word((a,)), comp_word((b,))).items():
                dadd(U, comp_sym(word_comp(wd)), Fr(c))
            rows.append((U, dict(prod)))
        for a in range(2, w-2):                  # R4,R5: z(a)*m2(b,c)
            for b in range(2, w-a):
                c = w-a-b
                if c < 1:
                    continue
                prod = self.cvec(dmul(znorm(a), self.csym(('m2', b, c))))
                U = {}
                for t in [('m3', a, b, c), ('m3', b, a, c), ('m3', b, c, a),
                          ('m2', a+b, c), ('m2', b, a+c)]:
                    dadd(U, t, Fr(1))
                rows.append((U, prod))
                U = {}
                for wd, cnt in shuffle(comp_word((a,)), comp_word((b, c))).items():
                    dadd(U, comp_sym(word_comp(wd)), Fr(cnt))
                rows.append((U, dict(prod)))
        us = [(w-1,)] if w >= 3 else []          # R6 Hoffman: u of weight w-1
        us += [(b, w-1-b) for b in range(2, w-1) if w-1-b >= 1]
        for u in us:
            U = {}
            R = {}
            yw = comp_word(u)
            for wd, cnt in shuffle((1,), yw).items():
                if wd == (1,)+yw:
                    assert cnt == 1
                    continue
                dadd(U, comp_sym(word_comp(wd)), Fr(cnt))
            for t in stuffle1(u):
                if len(t) == 1:                  # depth-1 merge -> RHS z(w)
                    R = dsum(R, znorm(t[0]))
                else:
                    dadd(U, comp_sym(t), Fr(-1))
            rows.append((U, R))
        return rows

    def _build(self, w):
        idx = {}

        def oi(s):
            if s not in idx:
                idx[s] = elim_key(s)
            return idx[s]
        piv = {}                                 # sym -> (U with coeff 1, R)

        def substitute(U, R):
            while True:
                hit = [s for s in U if s in piv]
                if not hit:
                    return U, R
                for s in hit:
                    c = U.pop(s, None)
                    if c is None:
                        continue
                    pU, pR = piv[s]
                    for s2, c2 in pU.items():
                        if s2 != s:
                            dadd(U, s2, -c*c2)
                    R = dsum(R, dscale(pR, -c))
        for U, R in self._rows(w):
            U = dict(U)
            U, R = substitute(U, R)
            if not U:
                assert not R, f"inconsistent relation at w={w}: residual {R}"
                continue
            lead = min(U, key=oi)
            f = Fr(1)/U[lead]
            U = {s: c*f for s, c in U.items()}
            R = dscale(R, f)
            piv[lead] = (U, R)
        for _ in range(len(piv)):                # back-substitute to closure
            changed = False
            for s, (U, R) in list(piv.items()):
                hit = [x for x in U if x != s and x in piv]
                if hit:
                    U2 = {k: v for k, v in U.items() if k != s}
                    U2, R2 = substitute(U2, R)
                    piv[s] = (dsum(U2, {s: Fr(1)}), R2)
                    changed = True
            if not changed:
                break
        surv = set()
        for s, (U, R) in piv.items():
            v = dict(R)
            for s2, c in U.items():
                if s2 != s:
                    dadd(v, (s2,), -c)
                    surv.add(s2)
            self.canon[s] = v
        self.survivors[w] = sorted(surv | {s for s in idx if s not in piv
                                           and s not in self.canon},
                                   key=elim_key, reverse=True)
