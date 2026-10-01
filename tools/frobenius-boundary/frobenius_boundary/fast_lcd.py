#!/usr/bin/env python3
"""fast_lcd.py — fast denominator-LCD: eval Kira denominator strings in the
ring Q[var] via operator overloading (Fraction coefficients), then sympy-lcm
the deduped coefficient dicts.  Replaces sympify+expand per string (the
~600s-CPU/eps bottleneck measured on the reference family: 1690 strings, 10-23KB each,
Horner-form) with ~seconds of exact Fraction arithmetic.  Charset verified
on P1: ()*+-/0-9^ d var only; no polynomial denominators ('/(' absent).
"""
import re
from fractions import Fraction

_INT = re.compile(r'(?<![\w])(\d+)')


class P:
    """Sparse polynomial in one variable over Fraction: {deg: coeff}."""
    __slots__ = ('c',)

    def __init__(self, c):
        self.c = {k: v for k, v in c.items() if v != 0}

    @staticmethod
    def _of(o):
        return o if isinstance(o, P) else P({0: Fraction(o)})

    def __add__(self, o):
        o = P._of(o)
        c = dict(self.c)
        for k, v in o.c.items():
            c[k] = c.get(k, Fraction(0)) + v
        return P(c)
    __radd__ = __add__

    def __neg__(self):
        return P({k: -v for k, v in self.c.items()})

    def __sub__(self, o):
        return self + (-P._of(o))

    def __rsub__(self, o):
        return P._of(o) + (-self)

    def __mul__(self, o):
        o = P._of(o)
        c = {}
        for k1, v1 in self.c.items():
            for k2, v2 in o.c.items():
                k = k1 + k2
                c[k] = c.get(k, Fraction(0)) + v1 * v2
        return P(c)
    __rmul__ = __mul__

    def __truediv__(self, o):
        if isinstance(o, P):
            if set(o.c) - {0}:
                raise ValueError('division by non-constant polynomial')
            o = o.c.get(0, Fraction(0))
        return P({k: v / o for k, v in self.c.items()})

    def __rtruediv__(self, o):
        raise ValueError('1/polynomial not supported')

    def __pow__(self, e):
        e = int(e)
        r = P({0: Fraction(1)}); b = self
        while e:
            if e & 1:
                r = r * b
            b = b * b
            e >>= 1
        return r


def parse_poly(ds, d_val, var='m2'):
    """Kira denominator string -> {deg: Fraction} at d = d_val (Fraction)."""
    s = ds.replace('^', '**')
    s = _INT.sub(r'F(\1)', s)
    env = {'__builtins__': {}, 'F': Fraction,
           var: P({1: Fraction(1)}), 'd': P({0: Fraction(d_val)})}
    v = eval(compile(s, '<den>', 'eval'), env)
    return P._of(v).c


def lcd_at_d(denom_strings, dd_rat, var='m2', fallback_sympy=True):
    """LCD over QQ of the denominator polynomials at rational d.
    Returns (coeffs ascending as sympy.Rational, degree)."""
    import sympy as sp
    xs, d = sp.symbols(var + ' d')
    d_val = Fraction(int(dd_rat.p), int(dd_rat.q))
    seen = {}
    for ds in denom_strings:
        try:
            c = parse_poly(ds, d_val, var)
        except Exception:
            if not fallback_sympy:
                raise
            e = sp.sympify(ds.replace('^', '**')).subs(d, dd_rat)
        else:
            if not c:                      # zero/constant denominator
                continue
            seen[frozenset(c.items())] = c
            continue
        p = sp.Poly(e, xs, domain='QQ')    # sympy fallback path
        c = {int(k[0]): Fraction(int(v.p), int(v.q))
             for k, v in p.as_dict().items()}
        seen[frozenset(c.items())] = c
    D = sp.Poly(1, xs, domain='QQ')
    for c in seen.values():
        if set(c) == {0}:
            continue                       # nonzero constants don't matter
        pc = sp.Poly({k: sp.Rational(v.numerator, v.denominator)
                      for k, v in c.items()}, xs, domain='QQ')
        D = sp.lcm(D, pc)
    D = D.monic()
    return [sp.Rational(c) for c in D.all_coeffs()[::-1]], D.degree()
