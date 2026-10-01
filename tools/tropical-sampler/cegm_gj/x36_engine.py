#!/usr/bin/env python3
"""X(3,6) Grassmannian string integral -- 4-dim tensor Gauss-Jacobi engine.

Definition (AHL 1912.08707, "Grassmannian string integrals", alpha'=1):
  I_{3,6}(s) = int_{R+^4} prod_{i} dx_i/x_i  prod_{a<b<c} X_abc(x)^{s_abc}
in the GS 2501.10805 positive network parametrization (setup_x36.py).
Variables (a,b,c,d) := (x11,x12,x21,x22).

All 20 minors are multilinear with nonnegative integer coefficients:
  const(=1): 123,124,125,126,134,234   (exponents fixed by momentum cons.)
  mono: 145=a, 156=ab, 345=c, 456=acd
  poly: 135=1+a; 136=1+a+ab; 146=a(1+b); 235=1+a+c; 245=a+c;
        346=c(1+b+d); 236=1+a+c+ab+bc+cd; 246=a+c+ab+bc+cd;
        256=ab+bc+cd; 356=c(b+d+ad)
Axis map x=t/(1-t).  Per-axis weight t^{A-1}(1-t)^{B-1} with
  A_i = sum_t s_t * mindeg_i(minor_t)   (monomial content exponent)
  B_i = -sum_t s_t * maxdeg_i(minor_t)
Residual G(t) = prod_polyfactors Phat^s * prod_i (1-t_i)^{D_i},
  D_i = sum_poly s_t * maxdeg_i(Phat_t),  computed in ln-space per point.
"""
import gmpy2 as g
from gmpy2 import mpfr, mpq
from fractions import Fraction as Fr
import numpy as np
import time, json, sys, os, argparse
from multiprocessing import Pool

LOG = None
def log(msg):
    line = "[%.3f] %s" % (time.time(), msg)
    print(line, flush=True)
    if LOG:
        LOG.write(line + "\n"); LOG.flush()

# ---------------- minors: multilinear structure tables ----------------------
# per-variable (mindeg, maxdeg) for each nonconstant minor, vars order (a,b,c,d)
TRIPLES = ['145','156','345','456','135','136','146','235','245','346',
           '236','246','256','356']
DEG = {  # var -> (min,max) per minor
 '145': {'a':(1,1)}, '156': {'a':(1,1),'b':(1,1)}, '345': {'c':(1,1)},
 '456': {'a':(1,1),'c':(1,1),'d':(1,1)},
 '135': {'a':(0,1)}, '136': {'a':(0,1),'b':(0,1)}, '146': {'a':(1,1),'b':(0,1)},
 '235': {'a':(0,1),'c':(0,1)}, '245': {'a':(0,1),'c':(0,1)},
 '346': {'b':(0,1),'c':(1,1),'d':(0,1)},
 '236': {'a':(0,1),'b':(0,1),'c':(0,1),'d':(0,1)},
 '246': {'a':(0,1),'b':(0,1),'c':(0,1),'d':(0,1)},
 '256': {'a':(0,1),'b':(0,1),'c':(0,1),'d':(0,1)},
 '356': {'a':(0,1),'b':(0,1),'c':(1,1),'d':(0,1)},
}
POLY = ['135','136','146','235','245','346','236','246','256','356']
VARS = ['a','b','c','d']

def axis_exponents(s):
    """s: dict triple->Fraction. Return A_i, B_i (Fractions)."""
    A, B = {}, {}
    for v in VARS:
        A[v] = sum(s[t]*Fr(DEG[t].get(v,(0,0))[0]) for t in TRIPLES)
        B[v] = -sum(s[t]*Fr(DEG[t].get(v,(0,0))[1]) for t in TRIPLES)
    return A, B

# ---------------- Gauss-Jacobi nodes/weights in gmpy2 ------------------------
def jacobi_pn_dpn(n, al, be, x):
    """P_n^{(al,be)}(x) and derivative, gmpy2 recurrence."""
    P0, dP0 = mpfr(1), mpfr(0)
    if n == 0: return P0, dP0
    P1 = (al+be+2)/2*x + (al-be)/2
    dP1 = (al+be+2)/2
    if n == 1: return P1, dP1
    for k in range(2, n+1):
        c0 = 2*k*(k+al+be)*(2*k+al+be-2)
        c1 = (2*k+al+be-1)*(2*k+al+be)*(2*k+al+be-2)
        c2 = (2*k+al+be-1)*(al*al-be*be)
        c3 = 2*(k+al-1)*(k+be-1)*(2*k+al+be)
        P2 = ((c1*x+c2)*P1 - c3*P0)/c0
        dP2 = ((c1*x+c2)*dP1 + c1*P1 - c3*dP0)/c0
        P0, P1, dP0, dP1 = P1, P2, dP1, dP2
    return P1, dP1

def gauss_jacobi_01(n, Afr, Bfr, prec):
    """Nodes/weights for int_0^1 t^{A-1}(1-t)^{B-1} f(t) dt, gmpy2 at prec bits.
    Jacobi (al,be)=(B-1,A-1) on [-1,1], t=(1+x)/2."""
    from scipy.special import roots_jacobi
    alf, bef = float(Bfr-1), float(Afr-1)
    xs, _ = roots_jacobi(n, alf, bef)
    with g.context(precision=prec+40):
        al, be = mpfr(Fr(Bfr-1).numerator)/mpfr(Fr(Bfr-1).denominator), \
                 mpfr(Fr(Afr-1).numerator)/mpfr(Fr(Afr-1).denominator)
        # C = 2^{al+be+1} G(n+al+1)G(n+be+1) / (n! G(n+al+be+1))
        lC = (al+be+1)*g.log(mpfr(2)) + g.lgamma(n+al+1)[0] + g.lgamma(n+be+1)[0] \
             - g.lgamma(mpfr(n+1))[0] - g.lgamma(n+al+be+1)[0]
        C = g.exp(lC)
        nodes, wts = [], []
        for x0 in xs:
            x = mpfr(x0)
            for _ in range(8):
                P, dP = jacobi_pn_dpn(n, al, be, x)
                dx = P/dP
                x = x - dx
                if abs(dx) < mpfr(2)**(-(prec+20)):
                    break
            P, dP = jacobi_pn_dpn(n, al, be, x)
            w = C/((1-x*x)*dP*dP)
            t = (1+x)/2
            nodes.append(t); wts.append(w)
        # map: int_0^1 t^{A-1}(1-t)^{B-1} f = 2^{-(A+B-1)} * sum w_i f(t_i)
        sc = g.exp(-(al+be+1)*g.log(mpfr(2)))
        wts = [w*sc for w in wts]
    return nodes, wts

def gj_selftest(n, Afr, Bfr, prec):
    """check moments m=0..min(2n-1,12) against Beta(A+m,B)."""
    nodes, wts = gauss_jacobi_01(n, Afr, Bfr, prec)
    worst = 0.0
    with g.context(precision=prec+40):
        for m in range(0, min(2*n, 13)):
            q = sum(w*t**m for t, w in zip(nodes, wts))
            A = mpfr(Fr(Afr).numerator)/mpfr(Fr(Afr).denominator)
            B = mpfr(Fr(Bfr).numerator)/mpfr(Fr(Bfr).denominator)
            ex = g.exp(g.lgamma(A+m)[0]+g.lgamma(B)[0]-g.lgamma(A+B+m)[0])
            rel = abs(q-ex)/ex
            worst = max(worst, float(rel))
    return worst

# ---------------- driver ----------------------------------------------------
def run_integral(s, n, prec, procs, tag):
    """s: dict triple->Fraction (14 nonconstant exponents). Returns mpfr."""
    t0 = time.time()
    A, B = axis_exponents(s)
    for v in VARS:
        assert A[v] > 0 and B[v] > 0, f"axis {v}: A={A[v]} B={B[v]} not >0"
    log(f"{tag}: axis exponents A={ {v:str(A[v]) for v in VARS} } "
        f"B={ {v:str(B[v]) for v in VARS} }")
    nodes, wts = {}, {}
    for v in VARS:
        nodes[v], wts[v] = gauss_jacobi_01(n, A[v], B[v], prec)
        st = gj_selftest(min(n,24), A[v], B[v], prec)
        log(f"{tag}: GJ axis {v} n={n} A={A[v]} B={B[v]} "
            f"moment-selftest(n=24) worst rel err={st:.2e}")
        assert st < 10**(-(prec*0.28)), f"GJ selftest failed axis {v}: {st}"
    nodesS = {v: [g.to_binary(x) for x in nodes[v]] for v in VARS}
    wtsS = {v: [g.to_binary(x) for x in wts[v]] for v in VARS}
    sv = {t: str(s[t]) for t in TRIPLES}
    t1 = time.time()
    log(f"{tag}: GJ setup done in {t1-t0:.2f}s; starting tensor sum "
        f"n^4={n**4} procs={procs}")
    chunks = []
    step = max(1, n//(procs*3))
    lo = 0
    while lo < n:
        hi = min(n, lo+step)
        chunks.append((lo, hi, prec, sv, nodesS, wtsS))
        lo = hi
    if procs > 1:
        with Pool(procs) as pool:
            parts = pool.map(worker_entry, chunks)
    else:
        parts = [worker_entry(ch) for ch in chunks]
    with g.context(precision=prec):
        tot = mpfr(0)
        for p in parts:
            tot += g.from_binary(p)
    t2 = time.time()
    log(f"{tag}: tensor sum done in {t2-t1:.2f}s  ({n**4/(t2-t1):.3g} pts/s) "
        f"value={g.format(tot,'.35e') if hasattr(g,'format') else repr(tot)}")
    return tot, t2-t1

def worker_entry(args):
    (ia_lo, ia_hi, prec, sv, nodesS, wtsS) = args
    ctx = g.context(precision=prec)
    g.set_context(ctx)
    s = {t: Fr(v) for t, v in sv.items()}
    gam = {t: mpfr(v.numerator)/mpfr(v.denominator)
           for t, v in ((tt, Fr(vv)) for tt, vv in sv.items())}
    N = {v: [g.from_binary(x) for x in nodesS[v]] for v in VARS}
    W = {v: [g.from_binary(x) for x in wtsS[v]] for v in VARS}
    na, nb, nc, nd = (len(N[v]) for v in VARS)
    Dv = {}
    for v in VARS:
        Dv[v] = sum(Fr(sv[t])*(DEG[t].get(v,(0,0))[1]-DEG[t].get(v,(0,0))[0])
                    for t in POLY)
    X = {v: [t/(1-t) for t in N[v]] for v in VARS}
    WD = {}
    for v in VARS:
        dv = mpfr(Fr(Dv[v]).numerator)/mpfr(Fr(Dv[v]).denominator)
        WD[v] = [W[v][i]*g.exp(dv*g.log(1-N[v][i])) for i in range(len(N[v]))]
    a_, b_, c_, d_ = X['a'], X['b'], X['c'], X['d']
    ga = gam
    z = mpfr(0)
    tot = mpfr(0)
    active = {t: (Fr(sv[t]) != 0) for t in POLY}
    BC = [[b_[j]*c_[k] for k in range(nc)] for j in range(nb)]
    CD = [[c_[k]*d_[l] for l in range(nd)] for k in range(nc)]
    if active['346']:
        L346 = [[ga['346']*g.log(1+b_[j]+d_[l]) for l in range(nd)]
                for j in range(nb)]
    W_bd = [[WD['b'][j]*WD['d'][l] for l in range(nd)] for j in range(nb)]
    use4 = active['236'] or active['246'] or active['256'] or active['356']
    for ia in range(ia_lo, ia_hi):
        a = a_[ia]; wa = WD['a'][ia]
        E135 = ga['135']*g.log(1+a) if active['135'] else z
        AB = [a*b_[j] for j in range(nb)]
        E136 = [ga['136']*g.log(1+a+AB[j]) for j in range(nb)] \
               if active['136'] else [z]*nb
        E146 = [ga['146']*g.log(1+b_[j]) for j in range(nb)] \
               if active['146'] else [z]*nb
        E235 = [ga['235']*g.log(1+a+c_[k]) for k in range(nc)] \
               if active['235'] else [z]*nc
        E245 = [ga['245']*g.log(a+c_[k]) for k in range(nc)] \
               if active['245'] else [z]*nc
        AD1 = [d_[l]*(1+a) for l in range(nd)]
        for jb in range(nb):
            Eab = E135 + E136[jb] + E146[jb]
            bj = b_[jb]; BCj = BC[jb]; Wj = W_bd[jb]
            L346j = L346[jb] if active['346'] else None
            for kc in range(nc):
                Eabc = Eab + E235[kc] + E245[kc]
                wac = wa*WD['c'][kc]
                S2 = a + c_[kc]
                uv = AB[jb] + BCj[kc]
                CDk = CD[kc]
                acc = z
                if use4:
                    for ld in range(nd):
                        P256 = uv + CDk[ld]
                        P246 = P256 + S2
                        E4 = ga['236']*g.log(1+P246) + ga['246']*g.log(P246) \
                             + ga['256']*g.log(P256) \
                             + ga['356']*g.log(bj+AD1[ld])
                        if L346j is not None: E4 += L346j[ld]
                        acc += Wj[ld]*g.exp(Eabc+E4)
                elif L346j is not None:
                    for ld in range(nd):
                        acc += Wj[ld]*g.exp(Eabc+L346j[ld])
                else:
                    sWj = sum(Wj)
                    acc = sWj*g.exp(Eabc)
                tot += acc*wac
    return g.to_binary(tot)

if __name__ == '__main__':
    print("module; use via run_pilot.py")
