#!/usr/bin/env python3
"""de_core.py — DE assembly + fast gmpy2 evaluation (family-agnostic).

The DE variable name and the set of massive propagator slots are parameters
(reference family: var='m2', mass_slots=0..7; dM/d(m2) = sum over massive lines of
a_l * M(line l dotted), reduced back to the basis by the Kira table).
"""
import re
import mpmath as mp
mp.mp.dps = max(mp.mp.dps, 50)      # set FIRST (import-dps footgun)
import gmpy2

from .kira_parse import parse_kira2math

_RAT = re.compile(r'(?<![\w.])(\d+)(?=\s*/)')


def compile_coef(cs):
    """Compile a Kira coefficient string; integer numerators before '/' are
    wrapped in gmpy2.mpfr (as _m) so division happens at full precision."""
    s2 = cs.replace('^', '**')
    s2 = _RAT.sub(r'_m(\1)', s2)
    return compile(s2, '<coef>', 'eval')


def load_masters(path_or_list):
    """Masters file: one FAM[i1,...,in] per line -> list of index tuples."""
    if isinstance(path_or_list, (list, tuple)):
        return [tuple(m) for m in path_or_list]
    out = []
    for L in open(path_or_list):
        m = re.search(r'\[([^\]]+)\]', L)
        if m:
            out.append(tuple(int(x) for x in m.group(1).split(',')))
    return out


def alpha_int(basis):
    """Integer part of the homogeneity weight: ai = (#ISP powers) - (sum dots)."""
    return [sum(-x for x in m if x < 0) - sum(x for x in m if x > 0)
            for m in basis]


def assemble_A(kira_path, basis, mass_slots, family=None):
    """A[i][j] = list[(int_prefactor, code|None)] such that
    dM_i/d(var) = sum_j (sum of prefactor*coef) M_j  on the cut block.
    Raising each massive line: d/d(var) (q^2-var)^-a = a*(q^2-var)^-(a+1)."""
    red = parse_kira2math(kira_path, family)
    IDX = {m: i for i, m in enumerate(basis)}
    N = len(basis)
    A = [[[] for _ in range(N)] for _ in range(N)]
    miss = []
    for i, m in enumerate(basis):
        for l in mass_slots:
            a = m[l]
            if a == 0:
                continue
            sh = list(m); sh[l] += 1; sh = tuple(sh)
            if sh in IDX:
                A[i][IDX[sh]].append((a, None))      # None => literal 1
            elif sh in red:
                for rhs, cs in red[sh]:
                    if rhs in IDX:
                        A[i][IDX[rhs]].append((a, compile_coef(cs)))
                    else:
                        miss.append((sh, rhs))
            else:
                miss.append((sh, None))
    return A, N, miss


def g2mp(v):
    """gmpy2.mpfr/mpc -> mp.mpf/mpc preserving full precision."""
    if isinstance(v, (int, float)):
        return mp.mpf(v)
    if isinstance(v, gmpy2.mpc):
        return mp.mpc(g2mp(v.real), g2mp(v.imag))
    if v == 0 or gmpy2.is_nan(v):
        return mp.mpf(0)
    m, e = gmpy2.mpfr(v).as_mantissa_exp()
    return mp.ldexp(mp.mpf(int(m)), int(e))


def mp2g(v):
    if isinstance(v, mp.mpc):
        return gmpy2.mpc(gmpy2.mpfr(mp.nstr(v.real, mp.mp.dps + 5)),
                         gmpy2.mpfr(mp.nstr(v.imag, mp.mp.dps + 5)))
    return gmpy2.mpfr(mp.nstr(mp.mpf(v), mp.mp.dps + 5))


def evalA(A, N, x, dd, var='m2'):
    """Evaluate the connection at var=x, dim d=dd.  Coefs evaluated in raw
    gmpy2 (~7x faster than mpmath), lifted to mp.matrix.  Precision follows
    the CURRENT mp.mp.dps."""
    gmpy2.get_context().precision = int(mp.mp.dps * 3.33) + 30
    if isinstance(x, mp.mpc) and x.imag != 0:
        gx = gmpy2.mpc(gmpy2.mpfr(mp.nstr(x.real, mp.mp.dps + 5)),
                       gmpy2.mpfr(mp.nstr(x.imag, mp.mp.dps + 5)))
    else:
        rx = x.real if isinstance(x, mp.mpc) else x
        gx = gmpy2.mpfr(mp.nstr(mp.mpf(rx), mp.mp.dps + 5))
    gd = gmpy2.mpfr(mp.nstr(mp.mpf(dd), mp.mp.dps + 5))
    M = mp.matrix(N, N)
    G = {'__builtins__': {}, '_m': gmpy2.mpfr}
    L = {var: gx, 'd': gd}
    for i in range(N):
        for j in range(N):
            if not A[i][j]:
                continue
            v = 0
            for pre, co in A[i][j]:
                v += pre * (eval(co, G, L) if co else 1)
            M[i, j] = g2mp(v)
    return M
