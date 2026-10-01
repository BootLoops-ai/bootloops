#!/usr/bin/env python3
"""wayfinder.flintexport (formerly the standalone flint_export) — FLINT fmpq_mpoly
fraction walker for bulk exact rational-function algebra (connection assembly,
reduction-row contraction, big substitute-compare).  Member of the Wayfinder
package: it is the exporter that writes the exact monomial-json connections
wayfinder.de_load.load_monomial_json reads.  Requires python-flint (and sympy
for the input trees); the rest of Wayfinder does not, so the package exposes
this module LAZILY (from wayfinder import flintexport).

Measured: 285.6 s for a full 22x22 3-var connection assembly (~1000x sympy).

Method: walk each sympy expression tree ONCE with values in the fraction field
Q(x1..xk) represented as (num, den) fmpq_mpoly pairs, gcd-reducing as we go
(FLINT gcd). Substitutions (e.g. d -> 4-2*eps, kinematic pins) are applied at
the leaves via the env. sympy is used only for the input tree; ALL arithmetic
is FLINT. Julia consumers: tools/counterweight/src/MpolyFeed.jl feed_matrix_from_chunks.

API:
  ctx = mk_ctx(('eps','s2','s12'))
  env = make_env(ctx, {sym_s2: gen(ctx,1), sym_d: poly(ctx,{(0,)*k:4,(1,)+(0,)*(k-1):-2}), ...})
  fr  = walk(sympy_expr, env, cache={}, ctx)   # Frac (num,den) pair, lazily reduced
  num_terms, den_terms = canonical_terms(fr)   # exact integer-coprime export form

Import:  from wayfinder.flintexport import mk_ctx, make_env, walk, canonical_terms
     or  sys.path.insert(0, <tools/wayfinder>); import flintexport   (flat)
Self-test: `python3 tools/wayfinder/flintexport.py` (or `--selftest`; exact
identities incl. a 3-var (d,s2,s12) cancel; exits nonzero on any mismatch).
Under the package battery it runs as tests/test_flintexport.py (skips by name
without python-flint or sympy).
"""
import sys
import math
from fractions import Fraction

from flint import fmpq, fmpq_mpoly_ctx, Ordering


def mk_ctx(names, ordering=Ordering.lex):
    """fmpq_mpoly context over the given variable names."""
    return fmpq_mpoly_ctx.get(tuple(names), ordering)


def _nvars(ctx):
    return len(ctx.gens())


def gen(ctx, i):
    """i-th generator as an fmpq_mpoly."""
    k = _nvars(ctx)
    return ctx.from_dict({tuple(1 if j == i else 0 for j in range(k)): fmpq(1)})


def const(ctx, q):
    """Constant fmpq as an fmpq_mpoly."""
    k = _nvars(ctx)
    return ctx.from_dict({(0,) * k: q if isinstance(q, fmpq) else fmpq(q)})


def poly(ctx, monom_dict):
    """{exponent_tuple: coeff} -> fmpq_mpoly (coeffs int/Fraction/fmpq)."""
    d = {}
    for m, c in monom_dict.items():
        if isinstance(c, fmpq):
            d[tuple(m)] = c
        elif isinstance(c, Fraction):
            d[tuple(m)] = fmpq(c.numerator, c.denominator)
        else:
            d[tuple(m)] = fmpq(c)
    return ctx.from_dict(d)


class Frac:
    """(num, den) fmpq_mpoly pair; gcd-reduce on add, cross-cancel on mul."""
    __slots__ = ("n", "d")

    def __init__(self, n, dd):
        self.n = n
        self.d = dd

    def reduce(self):
        if self.n.is_zero():
            return Frac(self.n, self.d * 0 + 1)
        g = self.n.gcd(self.d)
        return Frac(self.n / g, self.d / g)

    def add(self, o):
        g = self.d.gcd(o.d)
        a, b = self.d / g, o.d / g
        return Frac(self.n * b + o.n * a, g * a * b).reduce()

    def mul(self, o):
        g1 = self.n.gcd(o.d)
        g2 = o.n.gcd(self.d)
        return Frac((self.n / g1) * (o.n / g2), (self.d / g2) * (o.d / g1))

    def inv(self):
        return Frac(self.d, self.n)

    def pow(self, k):
        if k < 0:
            return self.inv().pow(-k)
        return Frac(self.n ** k, self.d ** k)


def make_env(ctx, subs_syms):
    """{sympy_symbol: fmpq_mpoly value} -> walk env (values wrapped as Frac)."""
    one = const(ctx, fmpq(1))
    return {sym: Frac(val, one) for sym, val in subs_syms.items()}


def walk(e, env, cache, ctx):
    """Evaluate sympy tree e in Q(ctx vars) via env. cache: id(e)->Frac (per-tree)."""
    key = id(e)
    got = cache.get(key)
    if got is not None:
        return got
    one = const(ctx, fmpq(1))
    if e.is_Symbol:
        r = env[e]
    elif e.is_Integer:
        r = Frac(const(ctx, fmpq(int(e))), one)
    elif e.is_Rational:
        r = Frac(const(ctx, fmpq(int(e.p), int(e.q))), one)
    elif e.is_Add:
        r = Frac(const(ctx, fmpq(0)), one)
        for a in e.args:
            r = r.add(walk(a, env, cache, ctx))
    elif e.is_Mul:
        r = Frac(one, one)
        for a in e.args:
            r = r.mul(walk(a, env, cache, ctx))
    elif e.is_Pow:
        if not e.exp.is_Integer:
            raise ValueError("non-integer exponent %s in Pow %s "
                             "(only integer powers are exact here)" % (e.exp, e))
        r = walk(e.base, env, cache, ctx).pow(int(e.exp))
    else:
        raise ValueError("unsupported sympy node %s" % type(e))
    cache[key] = r
    return r


def canonical_terms(fr):
    """Reduced Frac -> (num_terms, den_terms), integer coprime coefficients,
    den leading (ordering-max) coefficient positive. None if num == 0.
    terms = [[exponent_list, 'coeff_str'], ...] sorted monomial-descending."""
    fr = fr.reduce()
    nd, dd = fr.n.to_dict(), fr.d.to_dict()
    if not nd:
        return None
    L = 1
    for c in list(nd.values()) + list(dd.values()):
        c = Fraction(str(c))
        L = L * c.denominator // math.gcd(L, c.denominator)
    ni = {m: Fraction(str(c)) * L for m, c in nd.items()}
    di = {m: Fraction(str(c)) * L for m, c in dd.items()}
    G = 0
    for c in list(ni.values()) + list(di.values()):
        G = math.gcd(G, int(c))
    if G == 0:
        G = 1
    lead = max(di)
    sign = -1 if di[lead] < 0 else 1
    ni = {m: int(c) // G * sign for m, c in ni.items()}
    di = {m: int(c) // G * sign for m, c in di.items()}
    return ([[list(map(int, m)), str(c)] for m, c in sorted(ni.items(), reverse=True)],
            [[list(map(int, m)), str(c)] for m, c in sorted(di.items(), reverse=True)])


def _selftest():
    import sympy as sp
    d, s2, s12, eps = sp.symbols('d s2 s12 eps')
    ctx = mk_ctx(('eps', 's2', 's12'))
    d4 = poly(ctx, {(0, 0, 0): 4, (1, 0, 0): -2})       # 4 - 2 eps
    env = make_env(ctx, {s2: gen(ctx, 1), s12: gen(ctx, 2), d: d4})
    fails = 0

    # 1: exact zero after expansion
    e1 = (s2 + 1)**2 - s2**2 - 2*s2 - 1
    if canonical_terms(walk(e1, env, {}, ctx)) is not None:
        fails += 1
        print("FAIL 1: nonzero")

    # 2: gcd cancel (s2^2 - s12^2)/(s2 - s12) == s2 + s12
    e2 = (s2**2 - s12**2) / (s2 - s12) - (s2 + s12)
    r2 = walk(e2, env, {}, ctx).reduce()
    if not r2.n.is_zero():
        fails += 1
        print("FAIL 2: cancel")

    # 3: 3-var rational identity with the d -> 4-2eps pin
    #    (d-2)/(s2*s12) - (2-2*eps)/(s2*s12) == 0
    e3 = (d - 2) / (s2 * s12) - 2 * (1 - eps) / (s2 * s12)
    env3 = dict(env)
    env3[eps] = Frac(gen(ctx, 0), const(ctx, 1))
    r3 = walk(e3, env3, {}, ctx).reduce()
    if not r3.n.is_zero():
        fails += 1
        print("FAIL 3: 3-var pin")

    # 4: canonical export shape + positive-lead law
    e4 = (2*s2 + 2) / (-4*s2 + 4*s12)
    nt, dt = canonical_terms(walk(e4, env, {}, ctx))
    if not (nt and dt and all(len(t) == 2 for t in nt + dt)):
        fails += 1
        print("FAIL 4: export shape")

    # 5: non-integer exponents refuse by name (never silently floor)
    for bad in (sp.sqrt(s2) * s12, s2 ** sp.Rational(3, 2)):
        try:
            walk(bad, env, {}, ctx)
            fails += 1
            print("FAIL 5: non-integer exponent accepted:", bad)
        except ValueError:
            pass

    print("wayfinder.flintexport selftest:", "PASS" if fails == 0 else f"{fails} FAILS")
    return fails


if __name__ == '__main__':
    sys.exit(_selftest())
