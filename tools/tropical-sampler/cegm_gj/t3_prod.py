#!/usr/bin/env python3
"""Production driver — X(3,6) Grassmannian string integral values.

REUSES the certified pilot integrator (co-located modules): fan_x36.json
(52 unimodular simplicial subcones, certified), cone_engine.subcone_tasks,
x36_engine.gauss_jacobi_01, run_pilot kinematics + AHL Claim-1 certificate.
Adds ONLY: (a) outer-axis (i0) chunking so >52 procs load-balance and the
run emits timestamped progress lines for rate-fit ETAs; (b) per-run
JSON artifacts with provenance.  Nothing about the decomposition, kinematics,
or per-point arithmetic is re-derived; the inner loop is verbatim
cone_engine.eval_subcone with an i0 range parameter.

Stages:
  sanity                      chunked engine must reproduce pilot n=10 value
  cert   --point 2            certificate for the second kinematic point
  run    --point {1,2} --n N --prec BITS [--procs P] [--cpp C]
"""
import sys, os, json, time, argparse, hashlib
from fractions import Fraction as Fr
from multiprocessing import Pool

# The pilot modules + fan_x36.json are co-located in this directory (fan
# sha256 04807e002de2...). Run artifacts (prod.log, certificates, run JSONs)
# go to $CEGM_GJ_OUT or the caller's cwd, never into the tools tree.
HERE = os.path.dirname(os.path.abspath(__file__))
PILOT = HERE
OUT = os.environ.get('CEGM_GJ_OUT', os.getcwd())
sys.path.insert(0, PILOT)

import gmpy2 as g
from gmpy2 import mpfr
import x36_engine
from x36_engine import gauss_jacobi_01, log, TRIPLES
from cone_engine import subcone_tasks, mpfr_str
import run_pilot as RP

x36_engine.LOG = open(os.path.join(OUT, 'prod.log'), 'a')

# pilot n=10/160b generic value — sanity target
PILOT_N10 = '127.399414826922828643207977258924035788138193'
# pilot n=32/192b generic value — cert partner
PILOT_N32 = '127.399414828184508306057221011681649157001134'


def fan_hash():
    return hashlib.sha256(
        open(os.path.join(PILOT, 'fan_x36.json'), 'rb').read()).hexdigest()


def point1():
    """Pilot seed point s* — imported, not retyped."""
    return RP.generic_point()


def point2():
    """Second generic rational point. All 10 polynomial
    minor exponents < 0 (AHL chamber); 14 distinct rationals, denominators
    largely disjoint from point 1 to keep eventual PSLQ closure safe from
    point-specific accidents."""
    return {'135': Fr(-5, 12), '136': Fr(-3, 10), '146': Fr(-7, 16),
            '235': Fr(-4, 9),  '236': Fr(-8, 17), '245': Fr(-5, 14),
            '246': Fr(-2, 9),  '256': Fr(-7, 18), '346': Fr(-6, 17),
            '356': Fr(-3, 11),
            '145': Fr(5, 16),  '156': Fr(4, 11),  '345': Fr(7, 8),
            '456': Fr(5, 9)}


POINTS = {1: point1, 2: point2}


def eval_chunk(args):
    """cone_engine.eval_subcone verbatim, with outer-axis range [lo,hi)."""
    task, n, prec, lo, hi = args
    g.set_context(g.context(precision=prec))
    kappa = [Fr(k) for k in task['kappa']]
    nodes, wts = [], []
    for i in range(4):
        nd, wt = gauss_jacobi_01(n, kappa[i], Fr(1), prec)
        nodes.append(nd); wts.append(wt)
    facs = [(mpfr(Fr(gs).numerator)/mpfr(Fr(gs).denominator), monos)
            for gs, monos in task['facs']]
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
    gfacs = []
    for gam, monos in facs:
        by3 = {}
        for e in monos:
            by3.setdefault(e[3], []).append(e[:3])
        gfacs.append((gam, sorted(by3.items())))
    tot = mpfr(0)
    for i0 in range(lo, hi):
        p0 = powt[0][i0]; w0 = wts[0][i0]
        for i1 in range(n):
            p1 = powt[1][i1]; w01 = w0*wts[1][i1]
            for i2 in range(n):
                p2 = powt[2][i2]
                w012 = w01*wts[2][i2]
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


def eval_chunk_hoist(args):
    """Exact reordering of eval_chunk:
    factors depending on fewer than all 4 axes get gam*log(Q) precomputed on
    their own subgrid (n^|S| logs instead of n^4); only full-4-axis factors
    pay a live log per point.  Associativity-only change — same nodes, same
    arithmetic content; gated vs the plain path before production use."""
    task, n, prec, lo, hi = args
    g.set_context(g.context(precision=prec))
    kappa = [Fr(k) for k in task['kappa']]
    nodes, wts = [], []
    for i in range(4):
        nd, wt = gauss_jacobi_01(n, kappa[i], Fr(1), prec)
        nodes.append(nd); wts.append(wt)
    facs = [(mpfr(Fr(gs).numerator)/mpfr(Fr(gs).denominator), monos)
            for gs, monos in task['facs']]
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

    def qlog(gam, monos, i0, i1, i2, i3):
        idx = (i0, i1, i2, i3)
        Q = mpfr(0)
        for e in monos:
            t = mpfr(1)
            for ax in range(4):
                if e[ax]:
                    t *= powt[ax][idx[ax]][e[ax]]
            Q += t
        return gam*g.log(Q)

    live4 = []                     # full-support factors: live log per point
    tabf = {0: [], 1: [], 2: []}   # (S,T) complete at axis level max(S), 3∉S
    tab3 = []                      # (Srest,T): 3∈S, |S|<4; T[..][i3] lists
    r0 = range(lo, hi)
    rn = range(n)
    for gam, monos in facs:
        S = tuple(sorted({ax for m in monos for ax in range(4) if m[ax]}))
        if len(S) == 4:
            by3 = {}
            for e in monos:
                by3.setdefault(e[3], []).append(e[:3])
            live4.append((gam, sorted(by3.items())))
        elif 3 not in S:
            if len(S) == 0:
                T = qlog(gam, monos, 0, 0, 0, 0)
            elif len(S) == 1:
                a0, = S
                rr = r0 if a0 == 0 else rn
                T = {i: qlog(gam, monos, *[i if a == a0 else 0
                                           for a in range(3)], 0) for i in rr}
            elif len(S) == 2:
                a0, a1 = S
                T = {}
                for i in (r0 if a0 == 0 else rn):
                    T[i] = {j: qlog(gam, monos,
                                    *[i if a == a0 else j if a == a1 else 0
                                      for a in range(3)], 0) for j in rn}
            else:
                a0, a1, a2 = S
                T = {}
                for i in (r0 if a0 == 0 else rn):
                    T[i] = {j: {k: qlog(gam, monos,
                                *[i if a == a0 else j if a == a1 else
                                  k if a == a2 else 0 for a in range(3)], 0)
                                for k in rn}
                            for j in rn}
            tabf[max(S) if S else 0].append((S, T))
        else:
            Srest = tuple(a for a in S if a != 3)
            if len(Srest) == 0:
                T = [qlog(gam, monos, 0, 0, 0, k) for k in rn]
            elif len(Srest) == 1:
                a0, = Srest
                T = {i: [qlog(gam, monos, *[i if a == a0 else 0
                                            for a in range(3)], k)
                         for k in rn] for i in (r0 if a0 == 0 else rn)}
            else:
                a0, a1 = Srest
                T = {}
                for i in (r0 if a0 == 0 else rn):
                    T[i] = {j: [qlog(gam, monos,
                                     *[i if a == a0 else j if a == a1 else 0
                                       for a in range(3)], k) for k in rn]
                            for j in rn}
            tab3.append((Srest, T))

    def look(S, T, iv):
        v = T
        for a in S:
            v = v[iv[a]]
        return v

    tot = mpfr(0)
    z = mpfr(0)
    for i0 in r0:
        w0 = wts[0][i0]; p0 = powt[0][i0]
        iv = (i0, 0, 0)
        E0 = z
        for S, T in tabf[0]:
            E0 += look(S, T, iv) if S else T
        for i1 in rn:
            w01 = w0*wts[1][i1]; p1 = powt[1][i1]
            iv = (i0, i1, 0)
            E1 = E0
            for S, T in tabf[1]:
                E1 += look(S, T, iv)
            for i2 in rn:
                w012 = w01*wts[2][i2]; p2 = powt[2][i2]
                iv = (i0, i1, i2)
                E2 = E1
                for S, T in tabf[2]:
                    E2 += look(S, T, iv)
                slices = [look(S, T, iv) if S else T for S, T in tab3]
                coefs = []
                for gam, by3 in live4:
                    cl = []
                    for e3, elist in by3:
                        sacc = mpfr(0)
                        for e in elist:
                            sacc += p0[e[0]]*p1[e[1]]*p2[e[2]]
                        cl.append((e3, sacc))
                    coefs.append((gam, cl))
                acc = z
                for i3 in rn:
                    p3 = powt[3][i3]
                    E = E2
                    for sl in slices:
                        E += sl[i3]
                    for gam, cl in coefs:
                        Q = mpfr(0)
                        for e3, cv in cl:
                            Q += cv*p3[e3]
                        E += gam*g.log(Q)
                    acc += wts[3][i3]*g.exp(E)
                tot += acc*w012
    return g.to_binary(tot*task['det'])


_CLIB = None


def _clib():
    global _CLIB
    if _CLIB is None:
        import ctypes
        so = os.path.join(HERE, 'cchunk.so')
        lib = ctypes.CDLL(so)
        lib.eval_chunk_c.restype = ctypes.c_void_p
        lib.eval_chunk_c.argtypes = [
            ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int,
            ctypes.c_long,
            ctypes.POINTER(ctypes.c_char_p), ctypes.POINTER(ctypes.c_char_p),
            ctypes.c_int, ctypes.POINTER(ctypes.c_char_p),
            ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]
        lib.free_str.argtypes = [ctypes.c_void_p]
        _CLIB = (lib, ctypes)
    return _CLIB


def eval_chunk_cext(args):
    """C-MPFR implementation of the hoisted algorithm (cchunk.c).  Nodes and
    weights come from the SAME certified gauss_jacobi_01; C only replays the
    summation.  Gated vs the Python paths at the roundoff floor."""
    task, n, prec, lo, hi = args
    lib, ctypes = _clib()
    g.set_context(g.context(precision=prec + 40))
    nd = int((prec + 40)*0.30103) + 8
    kappa = [Fr(k) for k in task['kappa']]
    nodes_l, wts_l = [], []
    for i in range(4):
        nds, wt = gauss_jacobi_01(n, kappa[i], Fr(1), prec)
        nodes_l += [mpfr_str(x, nd).encode() for x in nds]
        wts_l += [mpfr_str(x, nd).encode() for x in wt]
    NodesArr = (ctypes.c_char_p * (4*n))(*nodes_l)
    WtsArr = (ctypes.c_char_p * (4*n))(*wts_l)
    gams, nms, monos_flat = [], [], []
    for gs, monos in task['facs']:
        gams.append(str(Fr(gs)).encode())
        nms.append(len(monos))
        for e in monos:
            monos_flat += list(e)
    GamArr = (ctypes.c_char_p * len(gams))(*gams)
    NmArr = (ctypes.c_int * len(nms))(*nms)
    MonoArr = (ctypes.c_int * len(monos_flat))(*monos_flat)
    p = lib.eval_chunk_c(n, prec, lo, hi, int(task['det']),
                         NodesArr, WtsArr, len(gams), GamArr, NmArr, MonoArr)
    sval = ctypes.string_at(p).decode()
    lib.free_str(p)
    g.set_context(g.context(precision=prec))
    return g.to_binary(mpfr(sval))


EVAL = {'plain': eval_chunk, 'hoist': eval_chunk_hoist, 'c': eval_chunk_cext}


def run_value(s, n, prec, procs, cpp, tag, engine='plain'):
    t0 = time.time()
    tasks = subcone_tasks(s)
    if cpp <= 0:
        cpp = max(1, (procs*4 + len(tasks) - 1)//len(tasks))
    bounds = [(i*n)//cpp for i in range(cpp + 1)]
    chunks = []
    for t in tasks:
        for i in range(cpp):
            lo, hi = bounds[i], bounds[i+1]
            if hi > lo:
                chunks.append((t, n, prec, lo, hi))
    K = len(chunks)
    pts = len(tasks)*n**4
    log(f"{tag}: {len(tasks)} subcones x cpp={cpp} = {K} tasks; n={n} "
        f"prec={prec} procs={procs} pts={pts} engine={engine} "
        f"load={os.getloadavg()[0]:.0f}")
    parts, done = [], 0
    step = max(1, K//30)
    with Pool(procs) as pool:
        for p in pool.imap_unordered(EVAL[engine], chunks):
            parts.append(p); done += 1
            if done % step == 0 or done == K:
                el = time.time() - t0
                log(f"{tag}: {done}/{K} chunks {el:.1f}s "
                    f"proj_total={el*K/done:.0f}s rate={pts*done/K/el:.3g}pts/s")
    g.set_context(g.context(precision=prec))
    tot = mpfr(0)
    for p in parts:
        tot += g.from_binary(p)
    wall = time.time() - t0
    log(f"{tag}: value={mpfr_str(tot, 75)} wall={wall:.2f}s "
        f"({pts/wall:.3g} pts/s)")
    return tot, wall, pts


def certificate(pid):
    s = POINTS[pid]()
    sfull = dict(s); sfull.update(RP.solve_const_exponents(s))
    ok, a, totc = RP.conservation_check(sfull)
    assert ok, f"conservation fails at particle {a}: {totc}"
    consts = {t: str(sfull[t]) for t in
              ['123', '124', '125', '126', '134', '234']}
    assert all(Fr(v) != 0 for v in consts.values()), "degenerate const exp"
    pc = RP.polytope_check(s)
    assert pc, "AHL Claim-1 polytope interior check FAILED"
    tasks = subcone_tasks(s)   # asserts kappa>0 in every subcone ray
    kmin = min(Fr(k) for t in tasks for k in t['kappa'])
    log(f"CERT point{pid}: s14={{{', '.join(t+':'+str(s[t]) for t in TRIPLES)}}}")
    log(f"CERT point{pid}: const-minor exponents (conservation): {consts}")
    log(f"CERT point{pid}: AHL Claim-1 LP + 8 perturbations: PASS; "
        f"kappa>0 all {len(tasks)} subcones (min {kmin}); conservation exact")
    return dict(point=pid, s14={t: str(s[t]) for t in TRIPLES},
                const_minors=consts, conservation='PASS-exact',
                ahl_claim1_lp='PASS(+8 perturbations)',
                kappa_positive_all_subcones=True, kappa_min=str(kmin),
                n_subcones=len(tasks), fan_sha256=fan_hash())


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('stage', choices=['sanity', 'cert', 'run'])
    ap.add_argument('--point', type=int, default=1)
    ap.add_argument('--n', type=int, default=10)
    ap.add_argument('--prec', type=int, default=160)
    ap.add_argument('--procs', type=int, default=16)
    ap.add_argument('--cpp', type=int, default=0)
    ap.add_argument('--engine', choices=['plain', 'hoist', 'c'],
                    default='plain')
    args = ap.parse_args()
    log(f"=== stage={args.stage} point={args.point} n={args.n} "
        f"prec={args.prec} procs={args.procs} pid={os.getpid()} "
        f"nice={os.nice(0)} ===")
    if args.stage == 'sanity':
        s = point1()
        v, wall, pts = run_value(s, 10, 160, args.procs, 2,
                                 f"SANITY engine={args.engine}",
                                 engine=args.engine)
        import mpmath as mp
        mp.mp.dps = 60
        rel = abs(mp.mpf(mpfr_str(v, 55)) - mp.mpf(PILOT_N10))/mp.mpf(PILOT_N10)
        ok = rel < mp.mpf('1e-38')
        log(f"SANITY: rel-vs-pilot-n10={mp.nstr(rel, 3)} "
            f"gate<1e-38: {'PASS' if ok else 'FAIL'}")
        print(json.dumps(dict(rel=float(rel), gate='PASS' if ok else 'FAIL')))
        sys.exit(0 if ok else 1)
    elif args.stage == 'cert':
        c = certificate(args.point)
        fn = os.path.join(OUT, f'point{args.point}_certificate.json')
        json.dump(c, open(fn, 'w'), indent=1)
        print(json.dumps(c))
    else:
        c = certificate(args.point)
        v, wall, pts = run_value(POINTS[args.point](), args.n, args.prec,
                                 args.procs, args.cpp,
                                 f"P{args.point} n={args.n} prec={args.prec}"
                                 f" [{args.engine}]", engine=args.engine)
        out = dict(point=args.point, n=args.n, prec_bits=args.prec,
                   procs=args.procs, wall_s=round(wall, 2), pts=pts,
                   pts_per_s=round(pts/wall, 1), engine=args.engine,
                   value=mpfr_str(v, 75), alpha_prime='1',
                   s14=c['s14'], const_minors=c['const_minors'],
                   fan_sha256=c['fan_sha256'],
                   loadavg=os.getloadavg()[0], nice=os.nice(0),
                   timestamp=time.strftime('%Y-%m-%d %H:%M:%S'))
        fn = os.path.join(OUT, f'run_p{args.point}_n{args.n}_b{args.prec}'
                          + ('' if args.engine == 'plain'
                             else '_' + args.engine) + '.json')
        json.dump(out, open(fn, 'w'), indent=1)
        print(json.dumps(dict(value=out['value'][:40], wall=out['wall_s'],
                              file=fn)))
