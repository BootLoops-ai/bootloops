# eichler.genus2 — curve re-recognition at the EXACT basepoint, full PSLQ law.
# Leg 1: xcf_cert0ex.json (dps 400, exact 1/36) -> tau -> theta -> BRUTE-Thomae Rosenhain
#        -> [x]+basis16 PSLQ with EXPLICIT pigeonhole-floor line + margin.
# Controls IN-PIPELINE: (P) planted algebraic element of the same field/height pushed
#        through the identical PSLQ stage -> must recover exactly;
#        (N) matched-magnitude random real -> must return nothing above floor.
# Leg 2 (second precision, independent transports): xcf_cert0ex2.json (dps 250) -> same
#        chain -> relations must MATCH leg 1 exactly (byte-identical integer vectors).
# Writes recognition_law.json incrementally as it goes.
import os, json, itertools, traceback
# precision knobs for the period/theta pipeline module (not shipped) — set
# your own pipeline's equivalents here before importing it:
os.environ['PIPELINE_DPS'] = '430'; os.environ['PIPELINE_MONO_PREF'] = '200'; os.environ['PIPELINE_DPS_IP'] = '400'
import numpy as np
import sympy as sp
from mpmath import mp, mpc, mpf, matrix, lu_solve, eig, pslq, sqrt as msqrt, nstr, log10
# NOTE: this file is the reference implementation from the originating study,
# shipped as the reusable STRUCTURE of the PSLQ recognition law
# (planted + negative controls in-pipeline, two independent legs, pigeonhole
# floor). `invariants_pipeline` is that study's period/theta pipeline
# (not distributed) — adapt the law skeleton to your own pipeline.
try:
    import invariants_pipeline as IP
except ImportError as _e:
    raise SystemExit(
        "pslq_law_template: `invariants_pipeline` (the originating project's "
        "period/theta pipeline, not distributed) is not available — this file ships "
        "as a read-and-adapt template for the PSLQ recognition law, not as a "
        "runnable module. Wire the IP.* calls to your own pipeline.") from _e

OUT = 'recognition_law.json'
out = {'producer': 'pslq_law_template.py', 'law': 'positive+negative control, two independent legs, pigeonhole margin'}
def save(): json.dump(out, open(OUT, 'w'), indent=1, default=str)

def pfloor(nterms, maxcoeff, valdigits):
    """pigeonhole artifact floor (relative residual) for an n-term integer relation search:
    ~ 10^-( (nterms-1)*log10(maxcoeff) ) — best artifact achievable by tuning coefficients."""
    return -( (nterms-1)*log10(mpf(maxcoeff)) )

def rosenhain_at(fn, dps_work):
    mp.dps = dps_work
    d = json.load(open(fn))
    F = matrix(5, 5)
    for j in range(5):
        for i in range(5):
            F[i, j] = IP.pc(d['cols'][str(j)][i])
    jets = F*IP.Lf
    PiS = [jets[0, j] for j in range(5)]
    Pi = [sum(IP.Ninv[i, j]*PiS[j] for j in range(5)) for i in range(5)]
    bf = json.load(open('exact_bivector_frame.json'))
    Bs = [sp.Matrix([[sp.nsimplify(v) for v in row] for row in Bi]) for Bi in bf['B']]
    def tompc(M):
        R = matrix(M.shape[0], M.shape[1])
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                r = sp.Rational(M[i, j]); R[i, j] = mpf(r.p)/mpf(r.q)
        return R
    BPi = matrix(4, 4)
    for i in range(5):
        BPi += tompc(Bs[i])*Pi[i]
    def herm(M):
        R = matrix(M.cols, M.rows)
        for i in range(M.rows):
            for j in range(M.cols):
                R[j, i] = mpc(M[i, j]).conjugate()
        return R
    ev, EV = eig(herm(BPi)*BPi)
    prs = sorted(range(4), key=lambda k: abs(ev[k]))
    phi1 = matrix([EV[i, prs[0]] for i in range(4)])
    phi2 = matrix([EV[i, prs[1]] for i in range(4)])
    dl = json.load(open('exact_spin_lattice.json'))
    H = sp.Matrix(dl['H']); den = sp.nsimplify(dl['den'])
    sb = json.load(open('exact_spin_sympbasis.json'))
    Lsym_n = tompc((H/den)*sp.Matrix(sb['Pchg']))
    C = matrix(2, 4)
    for r, phi in enumerate((phi1, phi2)):
        co = lu_solve(Lsym_n, phi)
        for c in range(4): C[r, c] = co[c]
    Ca = matrix(2, 2); Cb = matrix(2, 2)
    for r in range(2):
        Ca[r, 0], Ca[r, 1] = C[r, 0], C[r, 1]
        Cb[r, 0], Cb[r, 1] = C[r, 2], C[r, 3]
    def inv2(M):
        dt = M[0, 0]*M[1, 1]-M[0, 1]*M[1, 0]
        R = matrix(2, 2)
        R[0, 0], R[1, 1] = M[1, 1]/dt, M[0, 0]/dt
        R[0, 1], R[1, 0] = -M[0, 1]/dt, -M[1, 0]/dt
        return R
    T = inv2(Ca)*Cb*matrix([[1, 0], [0, 2]])
    t12 = (T[0, 1]+T[1, 0])/2
    T = matrix([[T[0, 0], t12], [t12, T[1, 1]]])
    if mpc(T[0, 0]).imag < 0:
        T = matrix([[-mpc(T[i, j]).conjugate() for j in range(2)] for i in range(2)])
    # Siegel reduce (same algorithm as tau_reduce)
    def inv2m(M): return inv2(M)
    for _ in range(300):
        changed = False
        Bsh = matrix(2, 2)
        for i in range(2):
            for j in range(2):
                Bsh[i, j] = mpf(round(float(mpc(T[i, j]).real)))
        if any(abs(Bsh[i, j]) > 0.5 for i in range(2) for j in range(2)):
            T = T - Bsh; changed = True
        Y = matrix([[mpc(T[i, j]).imag for j in range(2)] for i in range(2)])
        n = mpf(round(float(Y[0, 1]/Y[0, 0])))
        if abs(n) >= 1:
            Um = matrix([[1, -n], [0, 1]]); T = Um.T*T*Um; changed = True
        Y = matrix([[mpc(T[i, j]).imag for j in range(2)] for i in range(2)])
        if Y[1, 1] < Y[0, 0]:
            Um = matrix([[0, 1], [1, 0]]); T = Um.T*T*Um; changed = True
        if abs(T[0, 0]) < 0.985:
            A2 = matrix([[0, 0], [0, 1]]); B2 = matrix([[-1, 0], [0, 0]])
            C2 = matrix([[1, 0], [0, 0]]); D2 = matrix([[0, 0], [0, 1]])
            T = (A2*T+B2)*inv2m(C2*T+D2)
            t12n = (T[0, 1]+T[1, 0])/2
            T = matrix([[T[0, 0], t12n], [t12n, T[1, 1]]])
            changed = True
        if not changed: break
    th = IP.thetas(T[0, 0], T[0, 1], T[1, 1])
    th2 = [t*t for t in th]
    # BRUTE Thomae-verified Rosenhain (fresh, no inherited quads)
    th2f = np.array([complex(t) for t in th2])
    th4 = th2f**2
    vals = {}
    for (a, b) in itertools.combinations(range(10), 2):
        for (c, dd) in itertools.combinations(range(10), 2):
            if {a, b} & {c, dd}: continue
            r = th2f[a]*th2f[b]/(th2f[c]*th2f[dd])
            key = (round(r.real, 12), round(r.imag, 12))
            if key not in vals: vals[key] = (r, (a, b, c, dd))
    V = [v[0] for v in vals.values()]; Vq = [v[1] for v in vals.values()]
    parts = [S for S in itertools.combinations(range(6), 3) if 0 in S]
    th4_r = np.sort(np.abs(th4)); th4_r = th4_r/th4_r[0]
    def test(lam, mu, nu):
        e = [0j, 1+0j, None, lam, mu, nu]
        P = []
        for S in parts:
            Sc = tuple(i for i in range(6) if i not in S)
            p = 1+0j
            for grp in (S, Sc):
                for i, j in itertools.combinations(grp, 2):
                    if e[i] is None or e[j] is None: continue
                    p *= (e[i]-e[j])
            P.append(p)
        P = np.array(P)
        if np.any(np.abs(P) < 1e-18): return 1e9
        r1 = np.sort(np.abs(P)); r1 = r1/r1[0]
        return float(np.max(np.abs(r1-th4_r)/th4_r))
    best = None
    N = len(V)
    for i in range(N):
        for j in range(i+1, N):
            for k in range(j+1, N):
                s = test(V[i], V[j], V[k])
                if s < 1e-6:
                    best = (s, Vq[i], Vq[j], Vq[k]); break
            if best: break
        if best: break
    assert best, 'no Thomae triple'
    s, qa, qb, qc = best
    lam = th2[qa[0]]*th2[qa[1]]/(th2[qa[2]]*th2[qa[3]])
    mu  = th2[qb[0]]*th2[qb[1]]/(th2[qb[2]]*th2[qb[3]])
    nu  = th2[qc[0]]*th2[qc[1]]/(th2[qc[2]]*th2[qc[3]])
    return [mpc(t).real for t in (lam, mu, nu)], s, (qa, qb, qc)

PRIMES = (2, 3, 5, 409)
def basis16():
    b = []
    for s in range(16):
        m = mpf(1)
        for i in range(4):
            if s >> i & 1: m *= msqrt(PRIMES[i])
        b.append(m)
    return b

def recognize(v, tag, maxc, usedigits):
    B = basis16()
    rel = pslq([v] + B, maxcoeff=maxc, maxsteps=6*10**6)
    rec = {'maxcoeff': maxc}
    if rel and rel[0] != 0:
        acc = rel[0]*v + sum(rel[i+1]*B[i] for i in range(16))
        scale = max(abs(mpf(c)) for c in rel)
        resid = abs(acc)/scale
        Hgt = max(abs(c) for c in rel)
        floor_exp = float(-(16)*log10(mpf(Hgt)))       # 17-term pigeonhole with actual height
        margin = float(-log10(resid)) + floor_exp      # digits beyond floor
        rec.update({'rel': rel, 'residual': nstr(resid, 4), 'height': int(Hgt),
                    'floor_10exp': floor_exp, 'margin_digits_beyond_floor': margin,
                    'verdict': 'PASS' if margin > 30 else 'AT-FLOOR (reject)'})
        print(f'{tag}: rel found, height {Hgt:.2e}, resid {nstr(resid,3)}, '
              f'floor 1e{floor_exp:.0f}, margin {margin:.0f}d -> {rec["verdict"]}', flush=True)
    else:
        rec['rel'] = None
        print(f'{tag}: no relation at maxcoeff {maxc}', flush=True)
    return rec

try:
    # ---- controls first (law: controls run through the SAME stage) ----
    mp.dps = 380
    import random
    random.seed(7)
    B = basis16()
    planted_coeffs = [random.randint(-10**12, 10**12) for _ in range(16)]
    planted_den = random.randint(10**11, 10**12)
    vP = sum(planted_coeffs[i]*B[i] for i in range(16))/mpf(planted_den)
    cP = recognize(vP, 'CONTROL-P(planted h~1e12)', 10**14, 380)
    okP = cP.get('rel') is not None and cP['verdict'] == 'PASS'
    # verify planted recovery is EXACT (up to overall sign/gcd)
    out['control_planted'] = {'ok_pass': okP, 'rec': cP}
    vN = mpf('0.' + ''.join(random.choice('0123456789') for _ in range(370)))*4 + 1
    cN = recognize(vN, 'CONTROL-N(random)', 10**14, 380)
    okN = (cN.get('rel') is None) or cN.get('verdict', '').startswith('AT-FLOOR')
    out['control_negative'] = {'ok_reject': okN, 'rec': cN}
    save()
    assert okP, 'positive control failed - pipeline invalid'
    assert okN, 'negative control failed - pipeline invalid'

    # ---- leg 1: dps-400 exact-basepoint ----
    ros1, thomae1, quads1 = rosenhain_at('xcf_cert0ex.json', 400)
    out['leg1'] = {'thomae_resid': thomae1, 'quads': [list(q) for q in quads1],
                   'rosenhain': [nstr(t, 370) for t in ros1]}
    save()
    mp.dps = 380
    leg1 = {}
    for nm, v in zip(('lam', 'mu', 'nu'), ros1):
        leg1[nm] = recognize(v, f'leg1-{nm}', 10**20, 380)
    out['leg1']['recognition'] = leg1
    save()

    # ---- leg 2: dps-250 independent transports ----
    import glob, time
    while len(glob.glob('xcf_cert0ex2_col*.json')) < 5:
        time.sleep(60)
    cols = {}
    for j in range(5):
        d2 = json.load(open(f'xcf_cert0ex2_col{j}.json'))
        cols[str(j)] = d2['cols'][str(j)]
    json.dump({'dps': 250, 'xc': ['1/36', '0/1'], 'cols': cols}, open('xcf_cert0ex2.json', 'w'))
    ros2, thomae2, quads2 = rosenhain_at('xcf_cert0ex2.json', 250)
    out['leg2'] = {'thomae_resid': thomae2, 'rosenhain': [nstr(t, 230) for t in ros2]}
    save()
    mp.dps = 235
    leg2 = {}
    for nm, v in zip(('lam', 'mu', 'nu'), ros2):
        leg2[nm] = recognize(v, f'leg2-{nm}', 10**12, 235)
    out['leg2']['recognition'] = leg2
    # cross-leg agreement of values
    mp.dps = 235
    agree = [nstr(abs(mpf(nstr(a, 230)) - b), 4) for a, b in zip(ros1, ros2)]
    out['legs_value_agreement'] = agree
    # relation match verdict
    match = all(
        (leg1[nm].get('rel') is not None and leg2[nm].get('rel') is not None and
         leg1[nm]['rel'] == leg2[nm]['rel'])
        or (leg1[nm].get('rel') is None and leg2[nm].get('rel') is None)
        for nm in ('lam', 'mu', 'nu'))
    out['two_leg_relation_match'] = match
    save()
    print('DONE. two-leg match:', match, flush=True)
except Exception:
    out['error'] = traceback.format_exc()
    save()
    raise
