# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision.
r"""
bench_se3l.py -- arbitrary-precision QMC evaluation of a 3-loop self-energy
benchmark at eps^0.

The benchmark is the full 3-loop, 8-propagator scalar light-by-light
self-energy integral (propagator list below).  It is FINITE at
eps^0:  I = Gamma(2) * \int_simplex 1/F^2 .  After Cheng-Wu (x7=1) and the
closed-form integration of x6 (F is linear in x6 -> 1/(A B)), the surviving
integral is 6-dimensional over (x0..x5) in [0,inf)^6 of  1/(A B).

Target point:  s=-1, t=-1/3, m2=1   ->  I = 0.52495777354... (disteval, ~8 d).

This module builds the reduced 6-D integrand once (JSON-cached beside this
file; the build needs pySecDec + sympy), then runs the median-shift
CBC-lattice QMC engine from hiprec.py.
"""
import os, sys, time, json

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import sympy as sp
import gmpy2
from gmpy2 import mpfr

import hiprec

CACHE = os.path.join(HERE, 'se3l_reduced.json')

PROPAGATORS = [
    '(ell+p1)**2 - m2', '(ell+p1+p2)**2 - m2', '(ell-p4)**2 - m2',
    'k1**2 - m2', 'k2**2 - m2', '(ell-k1-k2)**2 - m2',
    '(ell-k1)**2', '(ell-k2)**2',
]
LOOPS = ['ell', 'k1', 'k2']
EXT = ['p1', 'p2', 'p3', 'p4']
RR = [
    ('p1*p1', 0), ('p2*p2', 0), ('p3*p3', 0), ('p4*p4', 0),
    ('p1*p2', 's/2'), ('p1*p3', 't/2'), ('p2*p3', '(-s-t)/2'),
    ('p1*p4', '(-s-t)/2'), ('p2*p4', 't/2'), ('p3*p4', 's/2'),
]


def build_reduced(s_val, t_val, m2_val):
    """Build & cache the reduced 6-D integrand source (1/(A B)) for given kinematics."""
    key = [str(s_val), str(t_val), str(m2_val)]
    if os.path.exists(CACHE):
        d = json.load(open(CACHE))
        if d.get('key') == key:
            return d['src'], d['vars']
    from feynman import build_UF, reduce_linear
    U, F, xs, L, N = build_UF(PROPAGATORS, LOOPS, EXT, RR)
    s, t, m2 = sp.symbols('s t m2')
    sub = {s: sp.nsimplify(s_val), t: sp.nsimplify(t_val), m2: sp.nsimplify(m2_val)}
    Fn = F.subs(sub)
    den_factors, remaining = reduce_linear(Fn, xs, b=2)
    assert len(den_factors) == 2, "expected exactly one linear elimination"
    A, B = den_factors
    expr = sp.expand(1/(A*B))
    src = str(expr)
    varnames = [str(v) for v in remaining]
    json.dump({'key': key, 'src': src, 'vars': varnames}, open(CACHE, 'w'))
    return src, varnames


# ----- integrand factory (module-level so fork workers can rebuild) ----------
_SRC = None
_VARS = None


def _factory(_):
    e = sp.sympify(_SRC)
    syms = sp.symbols(' '.join(_VARS))
    f = sp.lambdify(syms, e, modules=[{'sqrt': gmpy2.sqrt}, 'math'], cse=True)
    return f


def run(s_val=-1, t_val=sp.Rational(-1, 3), m2_val=1,
        N=64007, n_shifts=16, korobov_p=3, dps=40, nproc=None, seed=1):
    global _SRC, _VARS
    _SRC, _VARS = build_reduced(s_val, t_val, m2_val)
    dim = len(_VARS)
    res = hiprec.integrate_qmc(_factory, None, dim, N=N, n_shifts=n_shifts,
                               korobov_p=korobov_p, dps=dps, nproc=nproc, seed=seed)
    return res


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--N', type=int, default=64007)
    ap.add_argument('--shifts', type=int, default=16)
    ap.add_argument('--p', type=int, default=3)
    ap.add_argument('--dps', type=int, default=40)
    ap.add_argument('--nproc', type=int, default=0)
    ap.add_argument('--out', type=str, default='')
    args = ap.parse_args()
    nproc = args.nproc if args.nproc > 0 else None
    t0 = time.time()
    res = run(N=args.N, n_shifts=args.shifts, korobov_p=args.p,
              dps=args.dps, nproc=nproc)
    dt = time.time()-t0
    TARGET = '0.52495777354'
    print(f"=== 3-loop self-energy eps^0  s=-1,t=-1/3,m2=1  (CBC-QMC, 6-D) ===")
    print(f"N={args.N}  shifts={res['n_shifts']}  korobov_p={args.p}  dps={args.dps}")
    print(f"I     = {str(res['value'])[:48]}")
    print(f"error = {gmpy2.digits(res['error'], 10, 4)[0] if res['error']>0 else 'n/a'}  (~{res['digits']} digits)")
    print(f"target= {TARGET} (disteval ~8 d)")
    agree = -int(gmpy2.log10(abs(res['value']-mpfr(TARGET))+mpfr('1e-50')))
    print(f"agree vs disteval target: ~{agree} d")
    print(f"time  = {dt:.0f}s ({res['n_shifts']} shifts x {args.N} pts)")
    if args.out:
        json.dump({'value': str(res['value']), 'error': str(res['error']),
                   'digits': res['digits'], 'N': args.N, 'shifts': res['n_shifts'],
                   'korobov_p': args.p, 'dps': args.dps, 'time_s': dt},
                  open(args.out, 'w'), indent=1)
        print(f"wrote {args.out}")
