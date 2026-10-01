#!/usr/bin/env python3
"""Cone-decomposed evaluation of the X(3,6) Grassmannian string integral.

I = sum_{simplicial subcones} |det R| int_{(0,1)^4} prod y_i^{kappa_i-1}
        prod_{t active} Q_t(y)^{s_t} dy
with Q_t(0)=1, all coefficients +1 (analytic, >=1 on closed cube):
tensor Gauss-Jacobi per axis converges geometrically.

Stages: control / generic (kinematics from run_pilot.py).
"""
import gmpy2 as g
from gmpy2 import mpfr
from fractions import Fraction as Fr
import numpy as np
import json, time, os, sys, argparse
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from x36_engine import gauss_jacobi_01, axis_exponents, TRIPLES, VARS
import x36_engine
from x36_engine import log
from tropical_cones import PHAT_MONOS, POLY  # POLY order MUST match fan keys

FAN = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                  'fan_x36.json')))['cones']

def mpfr_str(v, nd=50):
    m, e, _ = v.digits(10, nd)
    neg = m.startswith('-')
    if neg:
        m = m[1:]
    return ('-' if neg else '') + f"0.{m}e{e}"

def subcone_tasks(s):
    """Build per-subcone data: rays, det, kappa (Fr), factor exponent tables."""
    A, B = axis_exponents(s)
    Avec = [A[v] for v in VARS]
    tasks = []
    for c in FAN:
        key = c['key']
        vsel = {t: PHAT_MONOS[t][key[ti]] for ti, t in enumerate(POLY)}
        # combined exponent vector A + sum_t s_t v_t
        comb = [Avec[i] + sum(s[t]*Fr(vsel[t][i]) for t in POLY)
                for i in range(4)]
        for simp in c['simplices']:
            rays = [c['rays'][i] for i in simp]
            det = abs(int(round(np.linalg.det(np.array(rays, dtype=float)))))
            assert det != 0
            kappa = []
            for r in rays:
                k = -sum(comb[i]*Fr(r[i]) for i in range(4))
                assert k > 0, f"kappa<=0: cone {key} ray {r}: {k}"
                kappa.append(k)
            facs = []
            for t in POLY:
                if s[t] == 0:
                    continue
                monos = []
                for m in PHAT_MONOS[t]:
                    e = []
                    for r in rays:
                        ei = sum((vsel[t][i]-m[i])*r[i] for i in range(4))
                        assert ei >= 0, (key, t, m, r, ei)
                        e.append(ei)
                    monos.append(tuple(e))
                facs.append((str(s[t]), monos))
            tasks.append(dict(det=det, kappa=[str(k) for k in kappa],
                              facs=facs))
    return tasks

def eval_subcone(args):
    task, n, prec = args
    g.set_context(g.context(precision=prec))
    kappa = [Fr(k) for k in task['kappa']]
    nodes, wts = [], []
    for i in range(4):
        nd, wt = gauss_jacobi_01(n, kappa[i], Fr(1), prec)
        nodes.append(nd); wts.append(wt)
    facs = [(mpfr(Fr(gs).numerator)/mpfr(Fr(gs).denominator), monos)
            for gs, monos in task['facs']]
    # per-axis power tables
    maxe = [0]*4
    for _, monos in facs:
        for e in monos:
            for i in range(4):
                maxe[i] = max(maxe[i], e[i])
    powt = []
    for i in range(4):
        tab = []
        for y in nodes[i]:
            row = [mpfr(1)]
            for p in range(maxe[i]):
                row.append(row[-1]*y)
            tab.append(row)
        powt.append(tab)
    # group each factor's monomials by inner-axis (axis 3) exponent
    gfacs = []
    for gam, monos in facs:
        by3 = {}
        for e in monos:
            by3.setdefault(e[3], []).append(e[:3])
        gfacs.append((gam, sorted(by3.items())))
    one = mpfr(1)
    tot = mpfr(0)
    for i0 in range(n):
        p0 = powt[0][i0]; w0 = wts[0][i0]
        for i1 in range(n):
            p1 = powt[1][i1]; w01 = w0*wts[1][i1]
            for i2 in range(n):
                p2 = powt[2][i2]
                w012 = w01*wts[2][i2]
                # coefficient lists per factor
                coefs = []
                for gam, by3 in gfacs:
                    cl = []
                    for e3, elist in by3:
                        sacc = mpfr(0)
                        for e in elist:
                            sacc += p0[e[0]]*p1[e[1]]*p2[e[2]]
                        cl.append((e3, sacc))
                    coefs.append((gam, cl))
                acc = mpfr(0)
                for i3 in range(n):
                    p3 = powt[3][i3]
                    E = mpfr(0)
                    for gam, cl in coefs:
                        Q = mpfr(0)
                        for e3, cv in cl:
                            Q += cv*p3[e3]
                        E += gam*g.log(Q)
                    acc += wts[3][i3]*g.exp(E)
                tot += acc*w012
    return g.to_binary(tot*task['det'])

def run_cone_integral(s, n, prec, procs, tag):
    t0 = time.time()
    tasks = subcone_tasks(s)
    log(f"{tag}: {len(tasks)} subcones; n={n} prec={prec} procs={procs} "
        f"pts={len(tasks)*n**4}")
    args = [(t, n, prec) for t in tasks]
    if procs > 1:
        with Pool(procs) as pool:
            parts = pool.map(eval_subcone, args)
    else:
        parts = [eval_subcone(a) for a in args]
    g.set_context(g.context(precision=prec))
    tot = mpfr(0)
    for p in parts:
        tot += g.from_binary(p)
    wall = time.time()-t0
    log(f"{tag}: value={mpfr_str(tot, 45)} wall={wall:.2f}s "
        f"({len(tasks)*n**4/wall:.3g} pts/s)")
    return tot, wall

if __name__ == '__main__':
    from run_pilot import split_point, split_oracle, generic_point, \
        solve_const_exponents, conservation_check, polytope_check
    ap = argparse.ArgumentParser()
    ap.add_argument('stage', choices=['control', 'generic'])
    ap.add_argument('--n', type=int, default=16)
    ap.add_argument('--prec', type=int, default=160)
    ap.add_argument('--procs', type=int, default=16)
    args = ap.parse_args()
    x36_engine.LOG = open(os.path.join(
        os.environ.get('CEGM_GJ_OUT', os.getcwd()), 'pilot.log'), 'a')
    log(f"=== CONE stage={args.stage} n={args.n} prec={args.prec} "
        f"procs={args.procs} pid={os.getpid()} ===")
    if args.stage == 'control':
        s, us, meta = split_point()
        import mpmath as mp
        oracle = split_oracle(us, dps=80)
        val, wall = run_cone_integral(s, args.n, args.prec, args.procs,
                                      f"CONE-CONTROL n={args.n} prec={args.prec}")
        v = mp.mpf(mpfr_str(val, 55))
        rel = abs(v-oracle)/oracle
        digits = float(-mp.log10(rel)) if rel > 0 else 99.0
        log(f"CONE-CONTROL n={args.n} prec={args.prec}: "
            f"agree_digits={digits:.2f} wall={wall:.2f}s")
        print(json.dumps(dict(n=args.n, prec=args.prec, digits=digits,
                              wall=wall, value=mp.nstr(v, 40),
                              oracle=mp.nstr(oracle, 45))))
    else:
        s = generic_point()
        sfull = dict(s); sfull.update(solve_const_exponents(s))
        ok, a, tot = conservation_check(sfull)
        assert ok
        assert polytope_check(s)
        val, wall = run_cone_integral(s, args.n, args.prec, args.procs,
                                      f"CONE-GENERIC n={args.n} prec={args.prec}")
        import mpmath as mp
        mp.mp.dps = 60
        v = mp.mpf(mpfr_str(val, 55))
        log(f"CONE-GENERIC n={args.n} prec={args.prec}: "
            f"value={mp.nstr(v, 45)} wall={wall:.2f}s")
        print(json.dumps(dict(n=args.n, prec=args.prec, wall=wall,
                              value=mp.nstr(v, 45))))
