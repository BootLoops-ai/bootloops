# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision.
"""
parametric.py -- GENERAL finite-parametric pathway for hiprec (Longhand's
route A: the UF-spec front end).

Extends hiprec beyond the a=0,b=2 scalar case of feynman.py: evaluates ANY
finite projective parametric integrand given as the UF-spec schema:

    spec = {
      'name': str, 'n_den': E,
      'polys': { tag: [ [monomial_exponents(list of E ints), [num, den]], ... ] },
      'numerator': { 'gamma': int,
                     'terms': [ [ [cn, cd], [ [tag, power], ... ] ], ... ] }
    }

    I_param = gamma * sum_t (cn/cd) * T[ prod_k P_tag^power ]

where T[.] is the projective simplex integral (each term must have total
projective degree -E; asserted).  Evaluation is done in the Cheng-Wu chart
x_chart = 1 over [0,inf)^(E-1), which equals the projective integral exactly,
and handed to hiprec.integrate_qmc (rank-1 CBC lattice, Korobov periodization,
random shifts -> statistical error estimate) in arbitrary precision (gmpy2).

This is what lets hiprec run 12-dim irreducible 4-loop integrands (13
denominators + an ISP numerator via the Mellin/derivative trick) as an
INDEPENDENT quadrature cross-check of a tropical-MC route on the identical
parametric object.  Measured reach: QMC at 12-D delivers a few digits, NOT a
30-d certificate -- the same honest limit as the 6-D self-energy benchmark.

Positive control: `box_spec()` re-derives the one-loop box (E=4 -> 3-dim) as a
spec; QMC on it must reproduce bench_box 0.12455705723270926971...
(selftest_parametric.py checks this at 5 sigma + >=4 d).
"""
import json
from fractions import Fraction

import gmpy2
from gmpy2 import mpfr


def load_spec(path, family=None):
    with open(path) as fh:
        d = json.load(fh)
    if family is not None:
        d = d[family]
    return d


def check_spec(spec):
    """Homogeneity/positivity checks: every poly homogeneous, every numerator
    term of total projective degree -E.  Returns (E, degrees dict)."""
    E = spec['n_den']
    degs = {}
    for tag, terms in spec['polys'].items():
        dset = {sum(m) for m, _ in terms}
        assert len(dset) == 1, f"{tag} not homogeneous: {dset}"
        assert all(len(m) == E for m, _ in terms), f"{tag} monomial length != {E}"
        degs[tag] = dset.pop()
    for (cn, cd), facs in spec['numerator']['terms']:
        pdeg = sum(p * degs[tag] for tag, p in facs)
        assert pdeg == -E, f"term {facs}: projective degree {pdeg} != -{E}"
    return E, degs


def make_eval(spec, chart=None):
    """Return (f, dim): gmpy2 callable f(*x) over [0,inf)^dim, dim = E-1,
    the FULL integrand (including gamma) in the Cheng-Wu chart x_chart = 1.
    Built at the CURRENT gmpy2 context precision (set it before calling).
    Sparse power-table evaluator: no sympy at eval time, fork-safe."""
    E, _ = check_spec(spec)
    chart = E - 1 if chart is None else chart
    live = [i for i in range(E) if i != chart]
    pos = {v: k for k, v in enumerate(live)}   # var index -> position in x
    maxdeg = [0] * len(live)
    compiled = {}
    for tag, terms in spec['polys'].items():
        ent = []
        for m, (num, den) in terms:
            sparse = tuple((pos[i], e) for i, e in enumerate(m) if i != chart and e)
            for p, e in sparse:
                if e > maxdeg[p]:
                    maxdeg[p] = e
            ent.append((mpfr(num) / mpfr(den), sparse))
        compiled[tag] = ent
    gamma = mpfr(spec['numerator']['gamma'])
    nterms = [(mpfr(cn) / mpfr(cd), tuple((tag, p) for tag, p in facs))
              for (cn, cd), facs in spec['numerator']['terms']]
    one = mpfr(1)

    def f(*x):
        pw = []
        for p, xi in enumerate(x):
            row = [one]
            acc = one
            for _ in range(maxdeg[p]):
                acc = acc * xi
                row.append(acc)
            pw.append(row)
        vals = {}
        for tag, ent in compiled.items():
            s = mpfr(0)
            for c, sparse in ent:
                t = c
                for p, e in sparse:
                    t = t * pw[p][e]
                s += t
            vals[tag] = s
        acc = mpfr(0)
        for c, facs in nterms:
            t = c
            for tag, pwr in facs:
                t = t * vals[tag] ** pwr
            acc += t
        return gamma * acc
    return f, len(live)


def spec_factory(cfg):
    """integrand_factory for hiprec.integrate_qmc.  cfg (picklable dict):
    {'json': path, 'family': name-or-None, 'chart': int-or-None}."""
    spec = load_spec(cfg['json'], cfg.get('family'))
    f, _ = make_eval(spec, cfg.get('chart'))
    return f


def integrate_spec(json_path, family=None, chart=None, N=8009, n_shifts=8,
                   korobov_p=3, dps=25, nproc=None, seed=1):
    """Convenience driver: spec file -> hiprec QMC result dict."""
    try:                       # imported as the package member longhand.parametric
        from . import hiprec
    except ImportError:        # imported flat, with this directory on sys.path
        import hiprec
    spec = load_spec(json_path, family)
    E, _ = check_spec(spec)
    return hiprec.integrate_qmc(spec_factory,
                                {'json': json_path, 'family': family, 'chart': chart},
                                E - 1, N=N, n_shifts=n_shifts, korobov_p=korobov_p,
                                dps=dps, nproc=nproc, seed=seed)


# --------------------------------------------------------------------------- #
#  positive control: the one-loop box as a spec (E=4 -> 3-dim)
# --------------------------------------------------------------------------- #
def box_spec(s=-1, t=Fraction(-1, 3), m2=1, M5=2):
    """One-loop massive box, I = Gamma(2) T[F0^-2] (a=0, b=2, E=4).
    F = (M5 x0 + m2 x1 + m2 x2 + m2 x3) U - s x0 x2 - u x1 x3,  U = sum x."""
    s, t, m2, M5 = Fraction(s), Fraction(t), Fraction(m2), Fraction(M5)
    u = -s - t
    mass = [M5, m2, m2, m2]
    F = {}
    for i in range(4):
        for jj in range(4):
            m = [0, 0, 0, 0]
            m[i] += 1
            m[jj] += 1
            F[tuple(m)] = F.get(tuple(m), Fraction(0)) + mass[i]
    for (i, jj, c) in [(0, 2, -s), (1, 3, -u)]:
        m = [0, 0, 0, 0]
        m[i] += 1
        m[jj] += 1
        F[tuple(m)] = F.get(tuple(m), Fraction(0)) + c
    terms = [[list(m), [c.numerator, c.denominator]]
             for m, c in sorted(F.items()) if c != 0]
    return {'name': 'box1L', 'n_den': 4, 'polys': {'F0': terms},
            'numerator': {'gamma': 1, 'terms': [[[1, 1], [['F0', -2]]]]}}
