"""Stage ALPHABET table for the fibration-basis integrator.

Letters = irreducible polynomials (fmpq_mpoly, primitive, sign-normalized) in the cell variables + shape symbols
(x_a, x_b, x_c, u, ub, v, vb).  Sources: (1) the variables themselves; (2) the per-support-group compatibility-graph
reduction of a linear-reducibility run (gate/lr/out_groups_symbolic.json: stages[k].letters for the chosen channel/group/gauge/perm);
(3) the input atom polynomials (factored).  Every denominator the integrator ever creates (resultants a b' - a' b,
leading coefficients [L,oo], constant terms [L,0]) is FACTORED (flint) and each irreducible factor is LOOKED UP here;
absent => OffAlphabet is raised (refuse) by design.  Nothing is ever expanded over a common denominator."""
import json, os
from fractions import Fraction as Fr
from flint import fmpq, fmpq_mpoly_ctx, fmpq_mpoly
import sympy as sp

import os as _os
# <root>/gate/lr relative to this file (<root>/symbolic/scripts); override with SURD_LR
LR_SYMBOLIC = _os.environ.get("SURD_LR", _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))), "gate", "lr", "out_groups_symbolic.json"))

class OffAlphabet(Exception):
    pass

def normalize(p):
    """primitive, sign-normalized representative of the class {c*p}; returns (unit c such that p = c * q, q)."""
    if p.is_zero():
        raise ValueError("zero polynomial")
    cont = p.content() if hasattr(p, "content") else None
    # fmpq_mpoly: make integer-primitive with positive leading coefficient
    d = p.to_dict()
    dens = [Fr(int(c.p), int(c.q)) for c in d.values()]
    from math import gcd
    from functools import reduce
    L = reduce(lambda a, b: a * b // gcd(a, b), [f.denominator for f in dens], 1)
    G = reduce(gcd, [abs(f.numerator * (L // f.denominator)) for f in dens])
    c = Fr(G, L)
    q = p * fmpq(L, G)
    if q.leading_coefficient() < 0:
        q = -q; c = -c
    return c, q

class Alphabet:
    def __init__(self, ctx):
        self.ctx = ctx
        self.names = list(ctx.names())
        self.letters = []          # fmpq_mpoly (normalized)
        self.key2idx = {}          # str(poly) -> idx
        self.meta = []             # dict(source=set, vars=set of var names)
        self._lin = {}             # (idx, var) -> (alpha poly, beta poly)
        self._fact_cache = {}
        for i, n in enumerate(self.names):
            self.add(ctx.gen(i) if hasattr(ctx, "gen") else ctx.gens()[i], "variable")

    def add(self, p, source):
        c, q = normalize(p)
        k = str(q)
        if k in self.key2idx:
            self.meta[self.key2idx[k]]["source"].add(source)
            return self.key2idx[k]
        idx = len(self.letters)
        self.letters.append(q); self.key2idx[k] = idx
        d = q.to_dict(); nv = len(self.names)
        vars_ = {self.names[j] for j in range(nv) if any(e[j] for e in d)}
        degs = {self.names[j]: max(int(e[j]) for e in d) for j in range(nv)}
        self.meta.append({"source": {source}, "vars": vars_, "deg": degs, "nterms": len(d)})
        return idx

    def add_from_string(self, s, source):
        expr = sp.sympify(s.replace("//", "/").replace("^", "**"), locals={n: sp.Symbol(n) for n in self.names})
        return self.add_sympy(expr, source)

    def add_sympy(self, expr, source):
        gens = [sp.Symbol(n) for n in self.names]
        P = sp.Poly(sp.together(expr).as_numer_denom()[0], *gens)
        d = {tuple(int(e) for e in m): fmpq(int(sp.Rational(c).p), int(sp.Rational(c).q)) for m, c in P.as_dict().items()}
        p = self.ctx.from_dict(d)
        # the LR letters are irreducible by construction; factor anyway and add each factor (robust to products in the input list)
        cont, facs = p.factor()
        ids = []
        for f, e in facs:
            if f.total_degree() == 0: continue
            ids.append(self.add(f, source))
        return ids

    def load_lr(self, channel, group, gauge, perm, path=LR_SYMBOLIC):
        d = json.load(open(path))
        hit = None
        for r in d["runs"]:
            if r["channel"] == channel and r["group"] == group and int(r["gauge"]) == int(gauge) and list(r["perm"]) == list(perm) and r["shape"] == "symbolic":
                hit = r; break
        if hit is None:
            raise KeyError("no LR run for %s %s g%s perm %s" % (channel, group, gauge, perm))
        self.lr_run = {k: hit[k] for k in ("channel", "group", "gauge", "perm", "order", "input_polys_sha256", "verdict")}
        for nm, s in hit["input_polys"].items():
            self.add_from_string(s, "lr_input:%s" % nm)
        self.stage_letter_ids = {}
        for st in hit["stages"]:
            if not st.get("computed"): continue
            k = int(st["stage"]); ids = []
            for s in st["letters"]:
                ids += self.add_from_string(s, "lr_stage%d" % k)
            self.stage_letter_ids[k] = sorted(set(ids))
        return self.lr_run

    # ---------------------------------------------------------------- factor over the alphabet
    def factor(self, p, context=""):
        """p (fmpq_mpoly) -> (Fraction c, {idx: e}) with p == c * prod letters^e ; raises OffAlphabet."""
        k = str(p)
        if k in self._fact_cache:
            return self._fact_cache[k]
        if p.is_zero():
            raise ValueError("factor(0) " + context)
        cont, facs = p.factor()
        c = Fr(int(cont.p), int(cont.q)); E = {}
        for f, e in facs:
            if f.total_degree() == 0:
                cf, _ = normalize(f); c *= cf ** int(e); continue
            cf, q = normalize(f)
            c *= cf ** int(e)
            kk = str(q)
            if kk not in self.key2idx:
                raise OffAlphabet("polynomial factor not in the alphabet (%s): %s" % (context, kk))
            i = self.key2idx[kk]; E[i] = E.get(i, 0) + int(e)
        res = (c, E)
        self._fact_cache[k] = res
        return res

    def linear_parts(self, idx, var):
        """letter idx = alpha*var + beta (must be of degree <= 1 in var) -> (alpha, beta) as fmpq_mpoly; raises if degree > 1."""
        key = (idx, var)
        if key in self._lin:
            return self._lin[key]
        p = self.letters[idx]; j = self.names.index(var)
        d = p.to_dict()
        if max(int(e[j]) for e in d) > 1:
            raise OffAlphabet("letter %s is not linear in %s" % (str(p), var))
        da, db = {}, {}
        for e, c in d.items():
            e2 = list(int(t) for t in e)
            if e2[j] == 1:
                e2[j] = 0; da[tuple(e2)] = c
            else:
                db[tuple(e2)] = c
        alpha = self.ctx.from_dict(da) if da else self.ctx.from_dict({(0,) * len(self.names): fmpq(0)})
        beta = self.ctx.from_dict(db) if db else self.ctx.from_dict({(0,) * len(self.names): fmpq(0)})
        self._lin[key] = (alpha, beta)
        return alpha, beta

    def depends(self, idx, var):
        return var in self.meta[idx]["vars"]

    def poly_of(self, E):
        """product of letters^e for e>0 entries (expanded fmpq_mpoly)"""
        one = self.ctx.from_dict({(0,) * len(self.names): fmpq(1)})
        r = one
        for i, e in E.items():
            if e > 0: r = r * self.letters[i] ** e
        return r

    def describe(self, idx):
        return str(self.letters[idx])
