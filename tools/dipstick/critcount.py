#!/usr/bin/env python3
# DIPSTICK member — verb `count`: Lagrange/ML-degree critical-system RANK
# probe, n>=5 regime; emits a msolve .ms — the msolve run + degree parse are
# the external steps. RANK PROBE, NOT A SOLVER.
"""critcount — critical-point / holonomic-rank counter for many-factor
Euler integrals prod p_i^{s_i}, WITHOUT symbolic expansion.
Lagrange (ML-degree) formulation: variables (x, lam), equations
  lam_i * p_i = u_i          (one per factor; excludes p_i = 0)
  sum_i lam_i * dp_i/dx_j = 0  (one per x)
All equations sparse, degree <= 1 + deg p. msolve counts solutions.
Control: X(3,6) must give 26 (verified independently by the product method)."""
import sympy as sp, itertools, random, sys, subprocess
def build(npts, seed=11):
    nm = npts - 4
    xs = sp.symbols(f'x1:{2*nm+1}')
    cols = [sp.Matrix([1,0,0]), sp.Matrix([0,1,0]), sp.Matrix([0,0,1]),
            sp.Matrix([1,1,1])]
    for i in range(nm):
        cols.append(sp.Matrix([1, xs[2*i], xs[2*i+1]]))
    Q = []
    for c in itertools.combinations(range(npts), 3):
        d = sp.expand(sp.Matrix.hstack(*[cols[i] for i in c]).det())
        if d.free_symbols:
            Q.append(d)
    random.seed(seed)
    u = [random.randint(2, 97) for _ in Q]
    lam = sp.symbols(f'l1:{len(Q)+1}')
    E = []
    for i, (q, ui) in enumerate(zip(Q, u)):
        E.append(sp.expand(lam[i]*q - ui))
    for xv in xs:
        E.append(sp.expand(sum(lam[i]*sp.diff(Q[i], xv)
                               for i in range(len(Q)))))
    allv = list(xs) + list(lam)
    return allv, E, len(Q)
if __name__ == '__main__':
    npts = int(sys.argv[1])
    allv, E, nq = build(npts)
    print(f'X(3,{npts}): {nq} factors, {len(allv)} vars, {len(E)} eqs, '
          f'max deg {max(sp.Poly(e, *allv).total_degree() for e in E)}',
          flush=True)
    sysstr = ',\n'.join(str(e).replace('**', '^') for e in E)
    fn = f'crit3{npts}.ms'
    open(fn, 'w').write(','.join(str(v) for v in allv)
                        + '\n1073741827\n' + sysstr + '\n')
    print(f'{fn} written', flush=True)
