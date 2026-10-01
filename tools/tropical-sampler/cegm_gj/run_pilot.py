#!/usr/bin/env python3
"""X(3,6) pilot driver: control gate + generic seed benchmark.

Stages (CLI):
  control  --n N --prec BITS [--procs P]   engine vs Gamma-form oracle
  generic  --n N --prec BITS [--procs P]   seed benchmark point
  qmc      --logn K                        float64 Sobol cross-check (generic)
Kinematics fixed in this file; all exponents exact Fractions; alpha'=1.
"""
import sys, os, json, time, argparse
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import x36_engine as E
from x36_engine import TRIPLES, POLY, VARS, DEG, axis_exponents, log
import x36_engine

ALL20 = ['123','124','125','126','134','234'] + TRIPLES  # 6 const + 14

def conservation_check(s20):
    """s20: dict all 20 triples -> Fr. Verify sum_{t ni a} s_t = 0 for a=1..6."""
    for a in '123456':
        tot = sum(v for t, v in s20.items() if a in t)
        if tot != 0:
            return False, a, tot
    return True, None, None

def solve_const_exponents(s14):
    """Solve the 6 constant-minor exponents from conservation (sympy)."""
    import sympy as sp
    cs = {t: sp.Symbol('c'+t) for t in ['123','124','125','126','134','234']}
    eqs = []
    for a in '123456':
        tot = sum(Fr(s14[t]) for t in TRIPLES if a in t)
        expr = sum(v for t, v in cs.items() if a in t) + sp.Rational(tot)
        eqs.append(expr)
    sol = sp.solve(eqs, list(cs.values()), dict=True)[0]
    return {t: Fr(str(sp.Rational(sol[cs[t]]))) for t in cs}

# ---------------- kinematic points ------------------------------------------
def split_point():
    """GS split kinematics: u_{i,j} from chosen positive (u0,u1,u2) per triangle.
    (u01,u11,u21)=(3/4,1/3,1/2), (u02,u12,u22)=(2/3,2/5,1/2)."""
    u01, u11, u21 = Fr(3,4), Fr(1,3), Fr(1,2)
    u02, u12, u22 = Fr(2,3), Fr(2,5), Fr(1,2)
    u31 = -(u01+u11+u21)
    u32 = -(u02+u12+u22)
    s = {t: Fr(0) for t in TRIPLES}
    s['156'] = u12
    s['456'] = u22
    s['145'] = u11 - u12 - u22
    s['346'] = u32
    s['345'] = u21 - u32 - u22
    s['235'] = u31
    meta = dict(u1=(str(u01), str(u11), str(u21)),
                u2=(str(u02), str(u12), str(u22)))
    return s, (u01, u11, u21, u02, u12, u22), meta

def split_oracle(us, dps=60):
    import mpmath as mp
    mp.mp.dps = dps
    u01, u11, u21, u02, u12, u22 = [mp.mpf(u.numerator)/u.denominator
                                    for u in us]
    J1 = mp.gamma(u01)*mp.gamma(u11)*mp.gamma(u21)/mp.gamma(u01+u11+u21)
    J2 = mp.gamma(u02)*mp.gamma(u12)*mp.gamma(u22)/mp.gamma(u02+u12+u22)
    return J1*J2

def generic_point():
    """Generic non-split rational point; all 10 poly-minor exponents negative
    (AHL convergence chamber), monomial exponents positive."""
    s = {'135': Fr(-2,5), '136': Fr(-3,7), '146': Fr(-4,9), '235': Fr(-5,11),
         '236': Fr(-6,13), '245': Fr(-1,3), '246': Fr(-2,7), '256': Fr(-3,8),
         '346': Fr(-4,11), '356': Fr(-5,13),
         '145': Fr(1,4), '156': Fr(2,5), '345': Fr(1,1), '456': Fr(1,2)}
    return s

# Newton polytope (content-extracted Phat) monomial exponent vectors (a,b,c,d)
PHAT_MONOS = {
 '135': [(0,0,0,0),(1,0,0,0)],
 '136': [(0,0,0,0),(1,0,0,0),(1,1,0,0)],
 '146': [(0,0,0,0),(0,1,0,0)],
 '235': [(0,0,0,0),(1,0,0,0),(0,0,1,0)],
 '245': [(1,0,0,0),(0,0,1,0)],
 '236': [(0,0,0,0),(1,0,0,0),(0,0,1,0),(1,1,0,0),(0,1,1,0),(0,0,1,1)],
 '246': [(1,0,0,0),(0,0,1,0),(1,1,0,0),(0,1,1,0),(0,0,1,1)],
 '256': [(1,1,0,0),(0,1,1,0),(0,0,1,1)],
 '346': [(0,0,0,0),(0,1,0,0),(0,0,0,1)],
 '356': [(0,1,0,0),(0,0,0,1),(1,0,0,1)],
}

def polytope_check(s, tries=9):
    """AHL Claim-1 convergence: A in int( sum_t (-s_t) N(Phat_t) ), s_t<0."""
    import numpy as np
    from scipy.optimize import linprog
    assert all(s[t] < 0 for t in POLY), "need all poly-minor exponents < 0"
    A, B = axis_exponents(s)
    Avec = np.array([float(A[v]) for v in VARS])
    cols, bounds_c = [], []
    tcount = []
    for t in POLY:
        for m in PHAT_MONOS[t]:
            cols.append(np.array(m, dtype=float))
        tcount.append(len(PHAT_MONOS[t]))
    ncols = len(cols)
    Amat = np.zeros((4+len(POLY), ncols))
    idx = 0
    for ti, t in enumerate(POLY):
        for m in PHAT_MONOS[t]:
            Amat[:4, idx] = m
            Amat[4+ti, idx] = 1.0
            idx += 1
    rng = np.random.default_rng(20260708)
    ok_all = True
    for k in range(tries):
        target = Avec if k == 0 else Avec*(1 + 0.01*rng.uniform(-1, 1, 4))
        bvec = np.concatenate([target,
                               np.array([float(-s[t]) for t in POLY])])
        r = linprog(np.zeros(ncols), A_eq=Amat, b_eq=bvec,
                    bounds=[(0, None)]*ncols, method='highs')
        if not r.success:
            ok_all = False
    return ok_all

def qmc_check(s, logn=21):
    """float64 Sobol estimate of the full integral (independent path check)."""
    import numpy as np
    from scipy.stats import qmc as sq
    from scipy.special import betaincinv, betaln
    A, B = axis_exponents(s)
    Af = {v: float(A[v]) for v in VARS}
    Bf = {v: float(B[v]) for v in VARS}
    D = {v: sum(float(s[t])*(DEG[t].get(v,(0,0))[1]-DEG[t].get(v,(0,0))[0])
                for t in POLY) for v in VARS}
    sob = sq.Sobol(4, scramble=True, seed=7)
    U = sob.random(2**logn)
    U = np.clip(U, 1e-14, 1-1e-14)
    T = {}
    for i, v in enumerate(VARS):
        T[v] = betaincinv(Af[v], Bf[v], U[:, i])
    X = {v: T[v]/(1-T[v]) for v in VARS}
    a, b, c, d = X['a'], X['b'], X['c'], X['d']
    lnG = sum(D[v]*np.log(1-T[v]) for v in VARS)
    P = {'135': 1+a, '136': 1+a+a*b, '146': 1+b, '235': 1+a+c, '245': a+c,
         '346': 1+b+d, '236': 1+a+c+a*b+b*c+c*d, '246': a+c+a*b+b*c+c*d,
         '256': a*b+b*c+c*d, '356': b+d+a*d}
    for t in POLY:
        st = float(s[t])
        if st != 0:
            lnG = lnG + st*np.log(P[t])
    Gv = np.exp(lnG)
    norm = np.exp(sum(betaln(Af[v], Bf[v]) for v in VARS))
    est = Gv.mean()*norm
    # crude error: half-sample spread
    h1 = Gv[:len(Gv)//2].mean()*norm
    h2 = Gv[len(Gv)//2:].mean()*norm
    return est, abs(h1-h2)

# ---------------- stages ------------------------------------------------------
def stage_control(n, prec, procs):
    s, us, meta = split_point()
    sfull = dict(s); sfull.update(solve_const_exponents(s))
    ok, a, tot = conservation_check(sfull)
    assert ok, f"conservation fails at particle {a}: {tot}"
    log(f"CONTROL point: {meta}; conservation OK; "
        f"s14={{{', '.join(t+':'+str(s[t]) for t in TRIPLES if s[t]!=0)}}}")
    import mpmath as mp
    oracle = split_oracle(us, dps=80)
    log(f"CONTROL oracle (Gamma-form) = {mp.nstr(oracle, 45)}")
    val, wall = E.run_integral(s, n, prec, procs, f"CONTROL n={n} prec={prec}")
    try:
        v = mp.mpf(f"{val:.55e}")
    except (ValueError, TypeError):
        v = mp.mpf(str(val))
    rel = abs(v-oracle)/oracle
    digits = float(-mp.log10(rel)) if rel > 0 else 99.0
    log(f"CONTROL n={n} prec={prec}: value={mp.nstr(v,40)} "
        f"agree_digits={digits:.2f} wall={wall:.2f}s")
    return dict(n=n, prec=prec, value=mp.nstr(v, 40),
                oracle=mp.nstr(oracle, 45), digits=digits, wall=wall)

def stage_generic(n, prec, procs):
    s = generic_point()
    sfull = dict(s); sfull.update(solve_const_exponents(s))
    ok, a, tot = conservation_check(sfull)
    assert ok, f"conservation fails at particle {a}: {tot}"
    consts = {t: str(sfull[t]) for t in ['123','124','125','126','134','234']}
    assert all(Fr(v) != 0 for v in consts.values()), "degenerate const exponent"
    pc = polytope_check(s)
    log(f"GENERIC point: s14={{{', '.join(t+':'+str(s[t]) for t in TRIPLES)}}}")
    log(f"GENERIC const-minor exponents (conservation): {consts}")
    log(f"GENERIC polytope interior check (AHL Claim 1 + 8 perturbations): "
        f"{'PASS' if pc else 'FAIL'}")
    assert pc
    val, wall = E.run_integral(s, n, prec, procs, f"GENERIC n={n} prec={prec}")
    import mpmath as mp
    mp.mp.dps = 60
    try:
        v = mp.mpf(f"{val:.55e}")
    except (ValueError, TypeError):
        v = mp.mpf(str(val))
    log(f"GENERIC n={n} prec={prec}: value={mp.nstr(v,45)} wall={wall:.2f}s")
    return dict(n=n, prec=prec, value=mp.nstr(v, 45), wall=wall)

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('stage', choices=['control', 'generic', 'qmc'])
    ap.add_argument('--n', type=int, default=16)
    ap.add_argument('--prec', type=int, default=160)
    ap.add_argument('--procs', type=int, default=16)
    ap.add_argument('--logn', type=int, default=21)
    args = ap.parse_args()
    x36_engine.LOG = open(os.path.join(
        os.environ.get('CEGM_GJ_OUT', os.getcwd()), 'pilot.log'), 'a')
    log(f"=== stage={args.stage} n={args.n} prec={args.prec} "
        f"procs={args.procs} pid={os.getpid()} ===")
    if args.stage == 'control':
        r = stage_control(args.n, args.prec, args.procs)
    elif args.stage == 'generic':
        r = stage_generic(args.n, args.prec, args.procs)
    else:
        s = generic_point()
        est, err = qmc_check(s, args.logn)
        log(f"QMC generic 2^{args.logn}: est={est!r} half-spread={err:.2e}")
        r = dict(est=est, spread=err)
    print(json.dumps(r, default=str))
