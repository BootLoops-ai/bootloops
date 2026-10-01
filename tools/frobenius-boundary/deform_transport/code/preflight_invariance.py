#!/usr/bin/env python3
"""preflight_invariance.py — STAGE 0 (mandatory pre-flight before the
rate-pilot stage of run_cal_d2.py).

Re-receipts, as fresh script-emitted measurements on the EXISTING
controls (their receipted configurations, the pinned kit):
 (A) charpoly-invariance of the residual commutant directions — rank-2
     Legendre p in {23,29} and Sym^5-Legendre p in {23,29}: kernel ladder
     census + STRUCTURE-LEVEL classification of every depth-0 residual
     direction (INVARIANT / SCALE-ONLY / MOVING) at modulus p^s12 grade,
     t in {1,2}; fiber-level corroboration column on rank-2.
 (B) the Sym^5 mult-1 stall run (same NF/s as the mult-5 endpoint —
     generous-window contrast); d_star contrast recorded.
 (C) (in-receipt part D) L6 J(5,1) commutant structure census
     at p=61 — exact charpoly-affecting dimension (the K3 adjudication
     frame the pilot's gauge-aware GE consumes).

VERDICT (script-emitted): PREFLIGHT-PASS licenses stage 1; any MOVING
residual direction = PREFLIGHT-DESIGN-WALL (K3 class) and the run STOPS.
The kit is imported as-is and its sha recorded — NOT edited."""
import json, os, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
MEMBER_DIR = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
import deform_kit as dk
from emit_receipt import emit as lk_emit

COMPONENT = 'frobenius-boundary/deform_transport (preflight)'
t0 = time.time()


def berkowitz_exact(M):
    """Division-free charpoly det(xI - M) over Z (kit algorithm, no modulus);
    ascending coefficients, leading 1."""
    n_ = len(M)
    C = [1, -M[0][0]]
    for k in range(1, n_):
        R = [M[k][j] for j in range(k)]
        Ccol = [M[i][k] for i in range(k)]
        A_ = [row[:k] for row in M[:k]]
        s = [1, -M[k][k]]
        v = Ccol[:]
        for t in range(1, k + 1):
            dot = sum(R[i] * v[i] for i in range(k))
            s.append(-dot)
            if t < k:
                v = [sum(A_[i][j] * v[j] for j in range(k)) for i in range(k)]
        newC = [0] * (len(C) + len(s) - 1)
        for i, si in enumerate(s):
            if si:
                for j, cj in enumerate(C):
                    newC[i + j] += si * cj
        C = newC[:k + 2]
    return list(reversed(C))


def E_of(evec, Ebasis, n_):
    return [[sum(evec[b] * Ebasis[b][i][j] for b in range(len(Ebasis)))
             for j in range(n_)] for i in range(n_)]


def cp_weights(cp_asc):
    """b_k = coeff vector with scaling weight k under E -> uE:
    cp ascending a_0..a_n (a_n=1); b_k = a_{n-k}."""
    return list(reversed(cp_asc))


def classify_pair(b_ref, b_got, p, K):
    """INVARIANT / SCALE-ONLY / MOVING at modulus p^K on weight-indexed
    charpoly vectors (b_0=1).  SCALE-ONLY = the exact multiplicative
    signature of E -> uE: b'_j^k * b_k^j == b_j^k * b'_k^j for all j,k."""
    mod = p ** K
    if all((x - y) % mod == 0 for x, y in zip(b_ref, b_got)):
        return 'INVARIANT'
    n = len(b_ref) - 1
    informative = 0
    for j in range(1, n + 1):
        for k in range(j + 1, n + 1):
            lhs = (pow(b_got[j] % mod, k, mod) * pow(b_ref[k] % mod, j, mod)) % mod
            rhs = (pow(b_ref[j] % mod, k, mod) * pow(b_got[k] % mod, j, mod)) % mod
            if (lhs - rhs) % mod:
                return 'MOVING'
            if b_ref[j] % mod or b_got[j] % mod or b_ref[k] % mod or b_got[k] % mod:
                informative += 1
    return 'SCALE-ONLY' if informative else 'DEGENERATE-BOTH'


def build_pin_rows(tr):
    """Rebuild the transport's OWN pinning rows (kit row law verbatim):
    base rows (raw, below the Laurent offset) + twisted tail rows."""
    S, T, n_, p = tr['S'], tr['T'], tr['n'], tr['p']
    offset = tr['offset']
    Ebasis, Ptw = tr['Ebasis'], tr['Ptw']
    nb = len(Ebasis)
    KE = tr['KE']
    s12v = tr['sigma'] + tr['sigma2']
    pinrows_idx = list(range(min(n_, 3)))
    Frows = [tr['assemble_rows'](E, pinrows_idx) for E in Ebasis]
    FrowsT = [[[Frows[b][ri][j] * Ptw for j in range(n_)]
               for ri in range(len(pinrows_idx))] for b in range(nb)]
    minv = None
    for b in range(nb):
        for ri in range(len(pinrows_idx)):
            for j in range(n_):
                poly = FrowsT[b][ri][j]
                for idx in range(T - 220, T, 4):
                    cv = int(poly[idx])
                    if cv:
                        v = dk.valp(cv, p) or 0
                        minv = v if minv is None else min(minv, v)
    dmax_try = KE - 4
    modcap = p ** (s12v + dmax_try)
    base_rows = []
    stride_neg = max(1, offset // 48) if offset else 1
    for ri in range(len(pinrows_idx)):
        for j in range(n_):
            for idx in range(0, offset, stride_neg):
                row = [int(Frows[b][ri][j][idx]) % modcap for b in range(nb)]
                if any(row):
                    base_rows.append(row)
    wlen = max(24, min(220, (T - offset) // 6))
    tail_rows = []
    for ri in range(len(pinrows_idx)):
        for j in range(n_):
            for idx in range(T - wlen, T - 4, 2):
                row = [int(FrowsT[b][ri][j][idx]) % modcap for b in range(nb)]
                if any(row):
                    tail_rows.append(row)
    return base_rows + tail_rows, s12v, dmax_try, minv, len(base_rows), len(tail_rows)


def invariance_measurement(tag, op, p, K, expect, fiber_lams=None, twist=None,
                           NF=None, s_twist=None):
    """One control-prime: transport re-run + reproducibility columns + kernel
    ladder census + structure-level invariance classification (+ optional
    fiber corroboration).  Returns the receipt row."""
    tw0 = time.time()
    row = {'control': tag, 'p': p, 'K_target': K}
    tr = dk.transport(op, p, K, twist=twist, NF=NF, s_twist=s_twist)
    r = tr['res']
    row['transport'] = {k: r[k] for k in ('NF', 'Ntot', 'jordan', 'nilpotent',
                                          'E_space_dim')}
    row['pin_reproduced'] = r['pin']
    row['reproducibility_vs_expected'] = {
        'expected': expect,
        'measured': {'d_star': r['pin']['depth_ladder_dstar'],
                     'E_space_dim': r['E_space_dim'],
                     'free_cols': r['pin']['free_cols'],
                     'n_gens': r['pin']['n_gens']},
    }
    row['reproducibility_vs_expected']['MATCH'] = bool(
        expect is None or
        (expect.get('d_star') == r['pin']['depth_ladder_dstar'] and
         expect.get('E_space_dim') == r['E_space_dim'] and
         expect.get('free_cols') == r['pin']['free_cols']))
    allrows, s12v, dmax_try, minv, nbase, ntail = build_pin_rows(tr)
    row['pin_rows_rebuilt'] = {'n_base': nbase, 'n_tail': ntail,
                               'min_stored_val_measured': minv, 's12': s12v,
                               'bounded_integrality_all_directions':
                               bool(minv is not None and minv >= s12v - 8)}
    d_star = tr['d_star']
    ladder = []
    dset = sorted(set([0, 1, 2, 3] + [max(0, d_star - 1), d_star,
                                      min(dmax_try, d_star + 1),
                                      min(dmax_try, d_star + 2)]))
    for d in dset:
        kk = dk.padic_kernel(allrows, p, s12v + d)
        ladder.append({'depth': d, 'n_gens': len(kk['gens']),
                       'free_cols': kk['free_cols'],
                       'n_pivots': len(kk['pivots'])})
        if d == 0:
            k0 = kk
        if d == d_star:
            kstar = kk
    row['kernel_ladder_census'] = ladder
    ray = kstar['gens'][0]
    row['ray_reproduced_from_rebuilt_rows'] = bool(
        len(kstar['gens']) == len(tr['gens']) and
        all((a - b) % (p ** (s12v + d_star)) == 0
            for ga, gb in zip(kstar['gens'], tr['gens']) for a, b in zip(ga, gb)))
    ray_cols = set(kstar['free_cols'])
    resid_cols = [c for c in k0['free_cols'] if c not in ray_cols]
    row['residual_directions_at_depth0'] = resid_cols
    n_ = tr['n']
    Eb = tr['Ebasis']
    Kcmp = s12v
    modc = p ** Kcmp
    Eref = E_of(ray, Eb, n_)
    cp_ref = berkowitz_exact(Eref)
    b_ref = cp_weights(cp_ref)
    detref = b_ref[n_]
    row['ray_structure_charpoly'] = {
        'det_val': dk.valp(detref % modc, p) if detref % modc else 'ZERO-mod-p^s12',
        'trace_val': dk.valp(b_ref[1] % modc, p) if b_ref[1] % modc else 'ZERO-mod-p^s12'}
    gmap = {c: g for c, g in zip(k0['free_cols'], k0['gens'])}
    dirs = {}
    any_moving = False
    for c in resid_cols:
        gc = gmap[c]
        drow = {'t_results': {}}
        Ealone = E_of(gc, Eb, n_)
        b_alone = cp_weights(berkowitz_exact(Ealone))
        drow['alone_det_val'] = (dk.valp(b_alone[n_] % modc, p)
                                 if b_alone[n_] % modc else 'ZERO-mod-p^s12')
        cls_set = set()
        for t in (1, 2):
            evec = [(a + t * b) for a, b in zip(ray, gc)]
            Et = E_of(evec, Eb, n_)
            b_t = cp_weights(berkowitz_exact(Et))
            cls = classify_pair(b_ref, b_t, p, Kcmp)
            drow['t_results'][str(t)] = cls
            cls_set.add(cls)
        if len(cls_set) == 1:
            drow['classification'] = cls_set.pop()
        else:
            drow['classification'] = 'MOVING(t-inconsistent: %s)' % sorted(cls_set)
        if drow['classification'].startswith('MOVING'):
            any_moving = True
        dirs[str(c)] = drow
    row['direction_classification'] = dirs
    row['certified_modulus'] = 'p^%d (generators canonical mod p^s12; comparisons exact at that depth)' % Kcmp
    row['any_MOVING'] = any_moving
    if fiber_lams:
        fib = {}
        F_ref = dk.assemble_full_F(tr, Eref)
        for lam in fiber_lams:
            frow = {}
            try:
                Fm_r, cert_r, tail_r = dk.eval_F_at(tr, F_ref, lam)
                cp_r = dk.berkowitz_charpoly(Fm_r, p ** cert_r)
                frow['ray'] = {'cert': cert_r, 'tail_val': tail_r}
                for c in resid_cols:
                    evec = [(a + b) for a, b in zip(ray, gmap[c])]
                    Fc = dk.assemble_full_F(tr, E_of(evec, Eb, n_))
                    try:
                        Fm_c, cert_c, tail_c = dk.eval_F_at(tr, Fc, lam)
                        cc = min(cert_r, cert_c)
                        cp_c = dk.berkowitz_charpoly(Fm_c, p ** cert_c)
                        agree = all((a - b) % (p ** cc) == 0
                                    for a, b in zip(cp_r, cp_c))
                        deep = None
                        for a, b in zip(cp_r, cp_c):
                            dv = dk.valp((a - b) % (p ** min(cert_r, cert_c, 12)), p)
                            if dv is not None:
                                deep = dv if deep is None else min(deep, dv)
                        frow['dir_%d' % c] = {
                            'combo_cert': cert_c, 'combo_tail_val': tail_c,
                            'equal_at_min_cert': bool(agree),
                            'raw_agreement_depth_uncertified': deep if deep is not None else '>=cap',
                            'register': 'corroboration-grade only (non-decaying direction caps the certificate)'}
                    except AssertionError as ex:
                        frow['dir_%d' % c] = {'integrality_break': str(ex.args)[:120]}
            except AssertionError as ex:
                frow['ray_eval_break'] = str(ex.args)[:120]
            fib[str(lam)] = frow
        row['fiber_corroboration'] = fib
    row['wall_s'] = round(time.time() - tw0, 1)
    print('[%s p=%d] d_star=%d resid_dirs=%s classes=%s wall=%.1fs' %
          (tag, p, d_star, resid_cols,
           {c: d['classification'] for c, d in dirs.items()}, row['wall_s']),
          flush=True)
    return row, tr


rec = {'spec': 'stage-0 pre-flight for the rate pilot',
       'register': 'stage-0 pre-flight; every verdict script-emitted; kit sha recorded in engine_shas',
       'engine_shas': {
           'deform_kit.py': dk.sha256_file(os.path.join(HERE, 'deform_kit.py')),
           'lfact_kit.py': dk.sha256_file(os.path.join(HERE, 'lfact_kit.py')),
           'SYM5_OP.json': dk.sha256_file(os.path.join(MEMBER_DIR, 'work', 'SYM5_OP.json'))},
       'controls': {}}

# ---- (A) rank-2 Legendre, CAL-A configuration (K=10, default twist mult 1, s=18)
opL = dk.op_legendre()
expA = {23: {'d_star': 20, 'E_space_dim': 2, 'free_cols': [1]},
        29: {'d_star': 20, 'E_space_dim': 2, 'free_cols': [1]}}
for p in (23, 29):
    row, _tr = invariance_measurement('rank2-Legendre(CAL-A cfg)', opL, p, 10,
                                      expA[p], fiber_lams=[4])
    rec['controls']['rank2_p%d' % p] = row

# ---- (A) Sym^5, CAL-B configuration (K=8, default twist mult 5, s=16)
opS = dk.op_sym5()
expB = {23: {'d_star': 6, 'E_space_dim': 6, 'free_cols': [5]},
        29: {'d_star': 3, 'E_space_dim': 6, 'free_cols': [5]}}
tr_s23 = None
for p in (23, 29):
    row, tr = invariance_measurement('Sym5-Legendre(CAL-B cfg)', opS, p, 8,
                                     expB[p], fiber_lams=None)
    rec['controls']['sym5_p%d' % p] = row
    if p == 23:
        tr_s23 = tr

# ---- (B) the Sym^5 mult-1 stall run (same NF/s as the mult-5 endpoint)
NF5 = tr_s23['res']['NF']
tw0 = time.time()
tr_m1 = dk.transport(opS, 23, 8, twist=[([1, -1], '1-x', 1)], NF=NF5, s_twist=16)
rec['stall_contrast_sym5_p23'] = {
    'config': 'same NF=%d, same s_twist=16; twist mult 1 vs 5 (generous-window contrast: '
              'mult-1 gets MORE window per twist-degree, so a stall is conservative)' % NF5,
    'd_star_mult5': tr_s23['d_star'],
    'd_star_mult1': tr_m1['d_star'],
    'pin_mult1': {'free_cols': tr_m1['res']['pin']['free_cols'],
                  'n_gens': tr_m1['res']['pin']['n_gens'],
                  'pivots': tr_m1['res']['pin']['pivots']},
    'stall_measured': bool(tr_m1['d_star'] < tr_s23['d_star']),
    'prose_expectation': 'd_star stalls at ~s/5 with mult 1 (s=16 -> ~3-class), tracks with mult 5',
    'register': 'measured-on-controls rate-law contrast; NOT a stop condition',
    'wall_s': round(time.time() - tw0, 1)}
print('[stall] mult5 d_star=%d mult1 d_star=%d' %
      (tr_s23['d_star'], tr_m1['d_star']), flush=True)

# ---- (C/D) L6 J(5,1) commutant structure census at p=61 (exact, no march)
tw0 = time.time()
op6 = dk.op_L6()
sb6 = dk.seed_basis(op6, op6.order, max(op6.resonances) + 10 if op6.resonances else 10)
sb6 = dk.echelonize_basis(sb6, op6)
from fractions import Fraction
import math as _math
Nfrac = sb6['N_struct']
den = 1
for r_ in Nfrac:
    for x in r_:
        den = den * x.denominator // _math.gcd(den, x.denominator)
Eb6 = dk.commutant_basis([[Fraction(x) for x in r_] for r_ in Nfrac], 61, op6.order)
E0 = [[sum((b + 1) * Eb6[b][i][j] for b in range(len(Eb6))) for j in range(6)]
      for i in range(6)]
b0 = cp_weights(berkowitz_exact(E0))
census = {'jordan': sb6['jordan'], 'commutant_dim': len(Eb6),
          'base_point': 'sum_(b+1)*E_b', 'base_det_nonzero': bool(b0[6] != 0),
          'directions': {}}
cnt = {'INVARIANT': 0, 'SCALE-ONLY': 0, 'MOVING': 0}
for b in range(len(Eb6)):
    cls_set = set()
    for t in (1, 2):
        Et = [[E0[i][j] + t * Eb6[b][i][j] for j in range(6)] for i in range(6)]
        bt = cp_weights(berkowitz_exact(Et))
        # exact integer classification (no modulus): INVARIANT / cross-ratio / MOVING
        if bt == b0:
            cls_set.add('INVARIANT')
            continue
        ok = True
        for j in range(1, 7):
            for k in range(j + 1, 7):
                if bt[j] ** k * b0[k] ** j != b0[j] ** k * bt[k] ** j:
                    ok = False
                    break
            if not ok:
                break
        cls_set.add('SCALE-ONLY' if ok else 'MOVING')
    cls = cls_set.pop() if len(cls_set) == 1 else 'MOVING(t-inconsistent)'
    census['directions'][str(b)] = cls
    cnt[cls.split('(')[0]] = cnt.get(cls.split('(')[0], 0) + 1
census['counts'] = cnt
census['charpoly_affecting_dim_measured'] = cnt['SCALE-ONLY'] + cnt['MOVING']
census['k3_frame'] = ('the pilot gauge-aware GE fires K3 iff the PINNED kernel retains '
                      'MOVING freedom (more than the one gamma-normalizable scale); this census '
                      'is the exact structure-level bound on what the pin can leave open')
census['wall_s'] = round(time.time() - tw0, 1)
rec['l6_commutant_census_p61'] = census
print('[L6 census] dim=%d counts=%s' % (len(Eb6), cnt), flush=True)

# ---- verdict (script-emitted)
movers = []
for key, row in rec['controls'].items():
    for c, d in row['direction_classification'].items():
        if d['classification'].startswith('MOVING'):
            movers.append('%s dir %s: %s' % (key, c, d['classification']))
repro_all = all(row['reproducibility_vs_expected']['MATCH'] for row in rec['controls'].values())
rec['reproducibility_all_match'] = repro_all
if movers:
    rec['verdict'] = ('PREFLIGHT-DESIGN-WALL: residual commutant direction(s) MOVE the charpoly '
                      'beyond the one normalizable scale on the controls [%s] — K3-class '
                      'design wall; the (3,5,1,1) pilot premise dissolves; the run STOPS'
                      % '; '.join(movers))
else:
    rec['verdict'] = ('PREFLIGHT-PASS: every residual commutant direction on the controls '
                      '(rank-2 p23/p29 CAL-A cfg; Sym^5 p23/p29 CAL-B cfg) classifies '
                      'charpoly-INVARIANT or SCALE-ONLY at p^s12 certified grade, t in {1,2} — '
                      'the K2 recipe-gap premise stands; stage-1 (3,5,1,1) pilot licensed to fire')
rec['wall_s_total'] = round(time.time() - t0, 1)
out = os.path.join(MEMBER_DIR, 'work', 'PREFLIGHT_RECEIPT.json')
lk_emit(rec, out, component=COMPONENT, script_path=os.path.abspath(__file__),
        member_dir=MEMBER_DIR)
dk.lint(out)
print('VERDICT:', rec['verdict'][:200])
