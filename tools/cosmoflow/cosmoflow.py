#!/usr/bin/env python3
"""cosmoflow.py -- thin CLI over polytope (general) and the triangle worked
example (examples/triangle: oracle, maxcut).  Letter alphabets for an
arbitrary site graph: `python3 -m cosmoflow.alphabet graph --nv N --edges ...`.

Subcommands:
    build    n_s               -- print B(y;X) for the n_s-gon (sympy, any n_s)
    oracle   key a lam c [dps] -- eval one triangle master (worked example)
    maxcut   a lam c [dps]     -- residue period at q_{G12}=0 (worked example)
    periods  a lam             -- (varpi0, varpi1, K^2) at (a,lam) (worked example)
"""
import sys, json
import mpmath as mp
import sympy as sp

from . import polytope
from .examples.triangle import oracle, maxcut     # the worked example


def cmd_build(argv):
    n_s = int(argv[0])
    ys = sp.symbols(f'y1:{n_s+1}')
    Xs = sp.symbols(f'X1:{n_s*(n_s-1)//2 + 1}')
    B = polytope.baikov_B(n_s, ys, Xs)
    print(f"B_{n_s}(y;X) =")
    print(sp.expand(B))


def cmd_oracle(argv):
    key, a, lam, c = argv[0], argv[1], argv[2], int(argv[3])
    dps = int(argv[4]) if len(argv) > 4 else 30
    mp.mp.dps = dps + 20
    v, d, w = oracle.eval_master(3, key, mp.mpf(a), mp.mpf(lam), c, dps=dps)
    print(json.dumps({'key': key, 'a': a, 'lam': lam, 'c': c,
                      'value': mp.nstr(v, dps), 'digits_est': round(d, 2),
                      'wall_s': round(w, 2)}))


def cmd_maxcut(argv):
    a, lam, c = argv[0], argv[1], int(argv[2])
    dps = int(argv[3]) if len(argv) > 3 else 80
    mp.mp.dps = dps + 20
    I = maxcut.period_residue(mp.mpf(a), mp.mpf(lam), c, dps=dps)
    print(json.dumps({'a': a, 'lam': lam, 'c': c,
                      'period_residue': mp.nstr(I, dps)}))


def cmd_periods(argv):
    a, lam = argv[0], argv[1]
    mp.mp.dps = 60
    aM, lM = mp.mpf(a), mp.mpf(lam)
    print(json.dumps({'a': a, 'lam': lam,
                      'K2': mp.nstr(maxcut.K2_modulus(aM, lM), 40),
                      'varpi0': mp.nstr(maxcut.varpi0(aM, lM), 40),
                      'varpi1': mp.nstr(maxcut.varpi1(aM, lM), 40)}))


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if not argv or argv[0] in ('-h', '--help'):
        print(__doc__); return 0
    cmd, rest = argv[0], argv[1:]
    {'build': cmd_build, 'oracle': cmd_oracle,
     'maxcut': cmd_maxcut, 'periods': cmd_periods}[cmd](rest)
    return 0


if __name__ == '__main__':
    sys.exit(main())
