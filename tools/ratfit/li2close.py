#!/usr/bin/env python3
r"""li2close.py -- close finite sums of log/linear ATOMS in ReLi2 + log
products, with an endpoint-divergence cancellation ledger and a pure-mpmath
evaluator emitter.

Generalized from an internal production closure script (wall constants
K0(u), K1).

INPUT: a list of atoms, each meaning

    c * ln|P0 + P1*t| / (W0 + W1*t)   integrated over t in (lo, hi),

with c, P0, P1, W0, W1 rational numbers or sympy expressions rational in one
parameter (e.g. u); lo/hi rational, lo may be 0 (endpoint eps-regularized) and
hi may be oo (Lambda-regularized).  Such atoms are what 1-d rational x log
integrals reduce to after partial fractions -- e.g. every corner-web wall-term
kernel of the face-flux anomaly recipe.

CLOSED FORM per atom (t0 = -W0/W1 the denominator root):
  generic (P(t0) != 0, P1 != 0):
    F(t) = (c/W1) [ ln|P(t0)| ln|W0+W1 t| - ReLi2( -P1 (t-t0)/P(t0) ) ]
  const-log (P1 == 0):     F(t) = (c/W1) ln|P0| ln|W0+W1 t|
  proportional (P(t0)==0): F(t) = (c/(2 W1)) ln|P0+P1 t|^2
with ReLi2(x) = Re Li2(x + i0), d/dx ReLi2 = -ln|1-x|/x.

ENDPOINT LEDGER: divergences at t->0+ and t->oo are carried symbolically as
LE = ln(eps), LL = ln(Lambda) (with their squares, from Li2(+-oo) asymptotics
including the branch constants -pi^2/6 / +pi^2/3).  close_atoms() sums all
atoms, collects the Poly in (LE, LL), and REQUIRES every non-finite monomial
coefficient to vanish -- symbolically if simplify can see it, else numerically
at rational sample points of the parameter (< tol, default 1e-40 at dps 60).
If any coefficient survives, DivergenceError is raised with the full ledger:
a non-cancelling web means the atom list is incomplete, not "regularizable".

SPURIOUS-POLE CAVEAT: a
W-root inside an atom's support is only admissible if a GROUP of atoms sharing
that W cancels it (sum_i c_i ln|P_i(t0)| = 0); this module does not re-check
interior regularity -- the originating script checked groups symbolically
and pointwise vs quadrature.  Keep doing both.

API
  Atom(c, W0, W1, P0, P1, lo, hi, tag='')     (sympy-fied on construction)
  ilin_atoms(...)   -- mechanical atom generation for _ilin-style
                       ln-difference/linear-denominator pieces (see docstring)
  close_atoms(atoms, param=None, domain=None, ...) -> Closure
      .finite (sympy expr in RLi2/log/Abs), .ledger, .pieces
  emit_mpmath(expr, name='K_closed', param=None) -> str   (standalone source)
  emitted_function(src, name) -> callable    (exec the emitted source)

POSITIVE CONTROLS (tools/tests/test_li2close.py): known integrals, e.g.
int_0^oo ln(1+t)/(t(1+t)) dt = pi^2/6 with the LL^2 cancellation across the
partial-fraction pair; and re-derivation of the reference production wall
constants K0(1/2), K1 from archived atom lists, matched to the archived
evaluators at >= 45 digits.
"""
import sympy as sp
from sympy import Abs, log, oo, pi

LE, LL = sp.symbols('LE LL')          # ln(eps), ln(Lambda)


class DivergenceError(Exception):
    """Endpoint divergences do not cancel across the atom sum."""

    def __init__(self, msg, ledger=None):
        super().__init__(msg)
        self.ledger = ledger


class RLi2(sp.Function):
    """Real dilogarithm ReLi2(x) = Re Li2(x + i0)."""
    nargs = 1

    def fdiff(self, argindex=1):
        z = self.args[0]
        return -log(Abs(1 - z)) / z


class Atom:
    """c * ln|P0 + P1*t| / (W0 + W1*t) on the open support (lo, hi)."""

    def __init__(self, c, W0, W1, P0, P1, lo, hi, tag=''):
        self.c, self.W0, self.W1 = map(sp.nsimplify, (c, W0, W1))
        self.P0, self.P1 = map(sp.nsimplify, (P0, P1))
        self.lo, self.hi = lo, hi
        self.tag = tag

    @classmethod
    def from_obj(cls, o):
        """Duck-typed adapter for foreign atom objects (e.g. assemble.Atom)."""
        return cls(o.c, o.W0, o.W1, o.P0, o.P1, o.lo, o.hi,
                   getattr(o, 'tag', ''))

    def integrand(self, tv):
        return self.c * log(Abs(self.P0 + self.P1 * tv)) / (self.W0
                                                            + self.W1 * tv)

    def __repr__(self):
        return (f"Atom[{self.tag}] {self.c}*ln|{self.P0}+{self.P1} t|"
                f"/({self.W0}+{self.W1} t) on [{self.lo},{self.hi}]")


def linco(expr, tv):
    """expr must be linear in tv: return (P0, P1) with expr = P0 + P1*tv."""
    e = sp.expand(expr)
    P1 = sp.simplify(sp.diff(e, tv))
    P0 = sp.simplify(e - P1 * tv)
    if sp.simplify(sp.diff(P1, tv)) != 0:
        raise ValueError(f'not linear in {tv}: {expr}')
    return P0, P1


def ilin_atoms(c_outer, A0, A1, B0, B1, lo_i, hi_i, tv, sup, tag):
    """Atoms of c_outer * _ilin(A0,A1,B0,B1,lo_i,hi_i) as a function of tv.

    _ilin = [ln|B0+B1*hi| - ln|B0+B1*lo| - ln|A0+A1*hi| + ln|A0+A1*lo|]/Wr,
    Wr = A0*B1 - A1*B0 (must be linear in tv); A1, B1 constants in tv;
    A0, B0, lo_i, hi_i at most linear in tv; sup = (lo, hi) outer support.
    This is the closed form of int dc/((A0+A1 c)(B0+B1 c)) over [lo_i, hi_i].
    """
    Wr0, Wr1 = linco(A0 * B1 - A1 * B0, tv)
    out = []
    for sgn, form in [(1, B0 + B1 * hi_i), (-1, B0 + B1 * lo_i),
                      (-1, A0 + A1 * hi_i), (1, A0 + A1 * lo_i)]:
        P0, P1 = linco(form, tv)
        out.append(Atom(sgn * c_outer, Wr0, Wr1, P0, P1, sup[0], sup[1], tag))
    return out


# ------------------------------------------------------------ antiderivative
def _pt0(a):
    return sp.simplify(a.P0 - a.P1 * a.W0 / a.W1)


def _fval(a, x):
    """Exact antiderivative at a finite regular point x (sympy expr)."""
    if a.P1 == 0:
        return a.c / a.W1 * log(Abs(a.P0)) * log(Abs(a.W0 + a.W1 * x))
    pt0 = _pt0(a)
    if pt0 == 0:
        return a.c / (2 * a.W1) * log(Abs(a.P0 + a.P1 * x)) ** 2
    t0 = -a.W0 / a.W1
    z = sp.simplify(-a.P1 * (x - t0) / pt0)
    return a.c / a.W1 * (log(Abs(pt0)) * log(Abs(a.W0 + a.W1 * x)) - RLi2(z))


def _fval_eps(a):
    """F at t = eps -> 0+, with LE = ln(eps).  Only for atoms with lo == 0."""
    if a.W0 != 0:                          # regular at 0
        return _fval(a, sp.Integer(0))
    if a.P1 == 0:
        return a.c / a.W1 * log(Abs(a.P0)) * (log(Abs(a.W1)) + LE)
    if _pt0(a) == 0:                       # P = (P1/W1) W; W0=0 forces P0=0
        return a.c / (2 * a.W1) * (log(Abs(a.P1)) + LE) ** 2
    # generic, W0 = 0: ln|P(0)| (ln|W1| + LE) - Li2(0)
    return a.c / a.W1 * log(Abs(a.P0)) * (log(Abs(a.W1)) + LE)


def _sign_on_domain(zs, param, samples):
    zsign = sp.sign(zs)
    if zsign in (sp.Integer(1), sp.Integer(-1)):
        return zsign
    if param is None or not samples:
        raise ValueError(f'cannot fix sign of {zs}: need param+domain samples')
    signs = {sp.sign(zs.subs(param, q)) for q in samples}
    if len(signs) != 1:
        raise ValueError(f'z-sign not constant on the domain: {zs}')
    return signs.pop()


def _fval_inf(a, param, samples):
    """F at t = Lambda -> +oo, with LL = ln(Lambda)."""
    if a.P1 == 0:
        return a.c / a.W1 * log(Abs(a.P0)) * (LL + log(Abs(a.W1)))
    pt0 = _pt0(a)
    if pt0 == 0:
        return a.c / (2 * a.W1) * (log(Abs(a.P1)) + LL) ** 2
    # z = -P1 (t - t0)/P(t0) -> sign(-P1/pt0) * oo; ReLi2 asymptotics:
    #   z -> -oo: -pi^2/6 - (1/2) ln^2|z|;   z -> +oo: +pi^2/3 - (1/2) ln^2|z|
    zsign = _sign_on_domain(sp.simplify(-a.P1 / pt0), param, samples)
    cpi = -pi ** 2 / 6 if zsign == -1 else pi ** 2 / 3
    return a.c / a.W1 * (log(Abs(pt0)) * (LL + log(Abs(a.W1))) - cpi
                         + sp.Rational(1, 2)
                         * (LL + log(Abs(a.P1 / pt0))) ** 2)


def atom_integral(a, param=None, samples=()):
    hi = _fval_inf(a, param, samples) if a.hi == oo else _fval(a, a.hi)
    lo = _fval_eps(a) if a.lo == 0 else _fval(a, a.lo)
    return hi - lo


# ------------------------------------------------------------------- closure
class Closure:
    def __init__(self, finite, ledger, pieces, param):
        self.finite = finite      # sympy expr (RLi2 / log / Abs / rationals)
        self.ledger = ledger      # list of dicts, one per divergent monomial
        self.pieces = pieces      # {(pow_LE, pow_LL): coeff expr}
        self.param = param


def _domain_samples(domain, nsamples):
    lo, hi = map(sp.nsimplify, domain)
    return [lo + sp.Rational(k, nsamples + 1) * (hi - lo)
            for k in range(1, nsamples + 1)]


def _numeric_zero(expr, param, samples, dps, tol):
    """max |expr| over samples (or at the point, if param-free) < tol?"""
    import mpmath
    f = sp.lambdify((param,) if param is not None else (), expr,
                    modules=['mpmath'])
    old = mpmath.mp.dps
    try:
        mpmath.mp.dps = dps
        pts = samples if param is not None else [None]
        worst = 0
        for q in pts:
            v = f(mpmath.mpf(sp.Rational(q).p) / sp.Rational(q).q) \
                if q is not None else f()
            worst = max(worst, abs(v))
        return worst < tol, float(worst)
    finally:
        mpmath.mp.dps = old


def close_atoms(atoms, param=None, domain=None, nsamples=6, dps=60,
                tol=None):
    """Assemble the atom sum; verify endpoint-divergence cancellation.

    atoms:  list of Atom (or duck-typed objects -> Atom.from_obj applied).
    param:  the sympy Symbol the atoms depend on (None = pure numbers).
    domain: (lo, hi) open interval of validity for param; REQUIRED whenever a
            sign or numeric-zero decision needs sampling.
    Returns Closure; raises DivergenceError if any LE/LL monomial coefficient
    fails both the symbolic-zero and the numeric-zero (< tol at dps) test.
    """
    atoms = [a if isinstance(a, Atom) else Atom.from_obj(a) for a in atoms]
    if tol is None:
        tol = 1e-40
    samples = _domain_samples(domain, nsamples) if domain else ()
    if param is not None and not samples:
        raise ValueError('param given but no domain to sample')
    totals = {}
    natoms = {}                     # divergent monom -> # contributing atoms
    for a in atoms:
        p = sp.Poly(sp.expand(atom_integral(a, param, samples)), LE, LL)
        for monom, coef in p.terms():
            totals[monom] = totals.get(monom, sp.Integer(0)) + coef
            if monom != (0, 0):
                natoms[monom] = natoms.get(monom, 0) + 1
    pieces, ledger = {}, []
    for monom in sorted(totals, reverse=True):
        coef = sp.expand(totals[monom])
        pieces[monom] = coef
        if monom == (0, 0):
            continue
        entry = {'monom': f'LE^{monom[0]} LL^{monom[1]}',
                 'atoms': natoms.get(monom, 0)}
        if coef == 0:
            entry['verdict'] = 'cancelled-at-expand'
            ledger.append(entry)
            continue
        c2 = sp.simplify(sp.expand_log(coef, force=True))
        if c2 == 0:
            entry['verdict'] = 'symbolic-zero'
        else:
            ok, worst = _numeric_zero(c2, param, samples, dps, tol)
            entry['verdict'] = (f'numeric-zero(max|.|={worst:.1e}'
                                f'<{tol:.0e} at {len(samples) or 1} pts)'
                                if ok else
                                f'FAIL(max|.|={worst:.1e}>={tol:.0e})')
            entry['residual_expr'] = sp.srepr(c2)
            if not ok:
                ledger.append(entry)
                raise DivergenceError(
                    f"divergence does not cancel: {entry['monom']} "
                    f'coefficient max|.| = {worst:.3e}', ledger)
        ledger.append(entry)
    finite = sp.expand(pieces.get((0, 0), sp.Integer(0)))
    return Closure(finite, ledger, pieces, param)


# ------------------------------------------------------------------- emitter
class _MpPrinter(sp.printing.str.StrPrinter):
    """sympy -> pure-mpmath source (exact rationals as mpf(p)/mpf(q))."""

    def _print_Rational(self, e):
        return f'(mpf({e.p})/mpf({e.q}))'

    def _print_Half(self, e):
        return '(mpf(1)/mpf(2))'

    def _print_Integer(self, e):
        return f'mpf({e.p})'

    def _print_Pi(self, e):
        return 'pi'

    def _print_Abs(self, e):
        return f'abs({self._print(e.args[0])})'

    def _print_Function(self, e):
        if e.func.__name__ == 'RLi2':
            return f'_rli2({self._print(e.args[0])})'
        if e.func.__name__ == 'log':
            return f'log({self._print(e.args[0])})'
        raise NotImplementedError(f'cannot emit function {e.func.__name__}')


_HEADER = '''\
# auto-generated by tools/li2close.py -- pure-mpmath closed-form evaluator
# (no fitted constants: exact rationals + log/ReLi2 of explicit arguments)
from fractions import Fraction as _F

from mpmath import mp, mpc, mpf, log, pi, polylog


def _rli2(x):
    """Real dilogarithm ReLi2(x) = Re Li2(x + i0)."""
    x = mpf(x)
    if x <= 1:
        return polylog(2, x)
    return polylog(2, mpc(x)).real
'''


def emit_mpmath(expr, name='K_closed', param=None):
    """Emit standalone pure-mpmath source for expr (a Closure.finite)."""
    src = _MpPrinter().doprint(sp.expand(expr))
    if param is not None:
        p = str(param)
        body = (f'def {name}({p}, dps=50):\n'
                f'    """closed form; {p} rational/float; mpf at dps."""\n'
                f'    old = mp.dps\n'
                f'    try:\n'
                f'        mp.dps = dps + 15\n'
                f'        {p} = _F({p})\n'
                f'        {p} = mpf({p}.numerator)/{p}.denominator\n'
                f'        v = ({src})\n'
                f'        mp.dps = dps\n'
                f'        return +v\n'
                f'    finally:\n'
                f'        mp.dps = old\n')
    else:
        body = (f'def {name}(dps=50):\n'
                f'    old = mp.dps\n'
                f'    try:\n'
                f'        mp.dps = dps + 15\n'
                f'        v = ({src})\n'
                f'        mp.dps = dps\n'
                f'        return +v\n'
                f'    finally:\n'
                f'        mp.dps = old\n')
    return _HEADER + '\n\n' + body


def emitted_function(src, name):
    """exec emitted source, return the named evaluator."""
    ns = {}
    exec(src, ns)
    return ns[name]
