#!/usr/bin/env python3
"""run_cal_d2.py — STAGE 1: the (3,5,1,1)-RATE PILOT (license gate).
Fires ONLY on the stage-0 PREFLIGHT-PASS receipt (preflight_invariance.py).

The L6-class operator at the declared pilot law (p=61, x0=2 + riders 3,
1/2), K_target=7, s law unchanged (s_twist = K+8 = 15), twist mults
(3,5,1,1) — the theta-lead-multiplicity-law rate.  Two truncation passes
(NF scale 1.0, 0.8).
GAUGE-AWARE GE: if the depth ladder leaves n_gens > 1, the surviving
generators are classified at structure level (charpoly(E) mod p^s12:
INVARIANT / SCALE-ONLY / MOVING vs a det-nondegenerate representative);
gauge/scale-only freedom + a nondegenerate representative = the transported
charpoly is well-defined up to the gamma-normalizable scale (GE-gauge PASS);
any MOVING freedom = K3' design wall.  All other bars verbatim CAL-D.
The kit is imported as-is and its sha recorded — NOT edited.
K2' band: measured-wall class bands with an in-runner stop if pass A exceeds
2x its band ceiling before pass B fires.  Comparator inputs (env-gated, loud
refusal): DEFORM_R2_GRID, DEFORM_F1_RECEIPTS, FROB_L6_OPERATOR — see
run_cal_d.py.  Emits work/CAL_D2_RECEIPT.json (snapshot emitter) with the
pricing row incl. measured max RSS."""
import json, os, sys, time, resource

HERE = os.path.dirname(os.path.abspath(__file__))
MEMBER_DIR = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
import deform_kit as dk
import lfact_kit as lk
from emit_receipt import emit as lk_emit

COMPONENT = 'frobenius-boundary/deform_transport (rate pilot)'

# ---- STAGE-0 GATE (mechanical): the pilot fires on the pre-flight receipt
pf_path = os.path.join(MEMBER_DIR, 'work', 'PREFLIGHT_RECEIPT.json')
pf = json.load(open(pf_path))
assert pf['verdict'].startswith('PREFLIGHT-PASS'), (
    'stage-0 gate refused: ' + pf['verdict'][:160])

t0 = time.time()
p = 61
K = int(sys.argv[1]) if len(sys.argv) > 1 else 7
mults = [int(x) for x in sys.argv[2].split(',')] if len(sys.argv) > 2 else [3, 5, 1, 1]
polys = [([1, -1], '1-x'), ([1, 1], '1+x'), ([1, 1, 1], 'Phi3'), ([1, 1, 1, 1, 1], 'Phi5')]
twist = [(c, nm, m) for (c, nm), m in zip(polys, mults) if m > 0]
BAND = {'pass_A_s': [20320, 29029], 'pass_B_s': [10297, 14710],
        'basis': '1.4-2x band (class estimate) on measured prior-recipe walls 14514.5/7354.8 s',
        'stop_trigger_s': {'pass_A': 58058, 'pass_B': 29420}}
rec = {'spec': 'stage 1 rate pilot (declared pilot law: p=61, x0=2, riders 3, 1/2)',
       'register': 'license gate; every verdict script-emitted; kit sha recorded',
       'stage0_gate': {'receipt': 'work/PREFLIGHT_RECEIPT.json',
                       'verdict_consumed': pf['verdict'][:200],
                       'producer_stamp': pf['producer']['stamp_utc']},
       'engine_shas': {
           'deform_kit.py': dk.sha256_file(os.path.join(HERE, 'deform_kit.py')),
           'lfact_kit.py': dk.sha256_file(os.path.join(HERE, 'lfact_kit.py')),
           'l6_operator': lk.l6_sha()},
       'pilot': {'p': p, 'x0': 2, 'riders': ['3', '1/2'], 'K_target': K,
                 'twist_mults': mults, 's_twist_law': 'K+8 (unchanged; rate carried by mults)'},
       'pricing_band': BAND}
op = dk.op_L6()

def _env_json(var, what):
    path = os.environ.get(var)
    if not path or not os.path.exists(path):
        raise SystemExit('run_cal_d2: comparator input missing — set %s to %s '
                         '(not distributed with the repo)' % (var, what))
    return json.load(open(path))

grid61 = _env_json('DEFORM_R2_GRID', 'a level-1 unit-root grid JSON')['grid']['61']


def berkowitz_exact(M):
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


def classify_pair(b_ref, b_got, p, Kc):
    mod = p ** Kc
    if all((x - y) % mod == 0 for x, y in zip(b_ref, b_got)):
        return 'INVARIANT'
    n = len(b_ref) - 1
    for j in range(1, n + 1):
        for k in range(j + 1, n + 1):
            lhs = (pow(b_got[j] % mod, k, mod) * pow(b_ref[k] % mod, j, mod)) % mod
            rhs = (pow(b_ref[j] % mod, k, mod) * pow(b_got[k] % mod, j, mod)) % mod
            if (lhs - rhs) % mod:
                return 'MOVING'
    return 'SCALE-ONLY'


def gauge_aware_E(tr):
    """Gauge-aware E selection.  Returns (E0, evec, ge_block).
    ge_block['mode'] in: classic-ray-unique | gauge-well-defined | K3-MOVING |
    no-nondegenerate-representative."""
    gens, Eb, n_, p_ = tr['gens'], tr['Ebasis'], tr['n'], tr['p']
    s12 = tr['sigma'] + tr['sigma2']
    d_star = tr['d_star']
    ge = {'n_gens': len(gens), 'd_star': d_star, 's12': s12}
    if len(gens) == 1 and (d_star or 0) >= 1:
        ge['mode'] = 'classic-ray-unique'
        return E_of(gens[0], Eb, n_), gens[0], ge
    # candidate det-nondegenerate representative: gens, pairwise sums, total
    cands = [(('g%d' % i,), g) for i, g in enumerate(gens)]
    for i in range(len(gens)):
        for j in range(i + 1, len(gens)):
            cands.append((('g%d+g%d' % (i, j),),
                          [a + b for a, b in zip(gens[i], gens[j])]))
    if len(gens) > 2:
        tot = [sum(g[b] for g in gens) for b in range(len(gens[0]))]
        cands.append((('sum-all',), tot))
    modc = p_ ** s12
    best = None
    for tag, vec in cands:
        Et = E_of(vec, Eb, n_)
        det = berkowitz_exact(Et)[0] * (1 if n_ % 2 == 0 else -1)
        dv = dk.valp(det % modc, p_) if det % modc else None
        row = {'candidate': tag[0], 'det_val_mod_ps12': dv if dv is not None else '>=s12'}
        if dv is not None and (best is None or dv < best[2]):
            best = (tag[0], vec, dv)
        ge.setdefault('candidate_scan', []).append(row)
    if best is None:
        ge['mode'] = 'no-nondegenerate-representative'
        return None, None, ge
    ge['reference'] = {'tag': best[0], 'det_val': best[2]}
    ref = best[1]
    b_ref = list(reversed(berkowitz_exact(E_of(ref, Eb, n_))))
    classes = {}
    any_moving = False
    for i, g in enumerate(gens):
        cls_set = set()
        for t in (1, 2):
            vec = [a + t * b for a, b in zip(ref, g)]
            b_t = list(reversed(berkowitz_exact(E_of(vec, Eb, n_))))
            cls_set.add(classify_pair(b_ref, b_t, p_, s12))
        classes['g%d' % i] = (cls_set.pop() if len(cls_set) == 1
                              else 'MOVING(t-inconsistent)')
        if classes['g%d' % i].startswith('MOVING'):
            any_moving = True
    ge['direction_classes'] = classes
    ge['certified_modulus'] = 'p^%d' % s12
    if any_moving:
        ge['mode'] = 'K3-MOVING'
        return None, None, ge
    ge['mode'] = 'gauge-well-defined'
    return E_of(ref, Eb, n_), ref, ge


def one_pass(NF_scale, tag):
    out = {'tag': tag}
    tw0 = time.time()
    NF = None
    if NF_scale != 1.0:
        dtarget = K + 6
        s_twist = dtarget + 2
        degD = sum((len(c) - 1) * m for (c, nm, m) in twist)
        NF = int(p * (s_twist * degD + 3.0 * dtarget + 8) * NF_scale)
    tr = dk.transport(op, p, K, twist=twist, NF=NF, progress=2000)
    r = tr['res']
    out['transport'] = {k: r[k] for k in ('NF', 'Ntot', 'jordan', 'nilpotent',
                                          'E_space_dim', 'inv_residual_certified_depth')}
    out['march'] = r['march']
    out['pin'] = r['pin']
    out['twist_law'] = r['twist']
    E0, evec, ge = gauge_aware_E(tr)
    out['GE_gauge'] = ge
    if E0 is None:
        out['wall_s'] = round(time.time() - tw0, 1)
        out['E_OBJECT'] = 'REFUSED(%s)' % ge['mode']
        return out, tr, None, None
    F = dk.assemble_full_F(tr, E0)
    E0, F, baseline = dk.normalize_ray(tr, E0, F)
    out['ray_baseline'] = baseline
    gd = dk.gd_residual(tr, F)
    out['GD_residual_worst_deficit'] = gd
    out['GD_pass'] = bool((gd or 0) <= tr['keep'] - tr['keep_target'] + 64)
    gb = dk.gamma_branches(E0, p, wt=4, n=6, Kg=max(6, min((tr['d_star'] or 6) + 2, 15)))
    out['gamma'] = {'val_detE0': gb['val_detE0'], 'target_val': gb['target_val'],
                    'gamma_valuation': gb.get('gamma_valuation'),
                    'n_branches': len(gb['branches'])}
    fibs = {}
    for x0 in (2, 3, (1, 2)):
        key = '2' if x0 == 2 else ('3' if x0 == 3 else '1/2')
        Fm, cert, tailv = dk.eval_F_at(tr, F, x0)
        mod = p ** cert
        cp = dk.berkowitz_charpoly(Fm, mod)
        c = dk.rp_coeffs(cp, p, cert)
        row = {'cert_depth': cert, 'tail_val': tailv}
        branches = []
        for br in gb['branches']:
            cn = dk.apply_gamma(c, br, p, cert)
            if cn is None:
                continue
            iso = dk.isolate_factor(cn, p, min(cert, K), wt=4, n=6)
            branches.append({'sign': br['sign'], 'root_mod_p': br['root_mod_p'],
                             'gamma_val': br['gamma_val'], 'iso': iso})
        row['branches'] = branches
        row['unit_root_mod_p3'] = None
        try:
            mod3 = p ** min(3, cert)
            cps = [int(x) % mod3 for x in dk.berkowitz_charpoly(Fm, mod3)]
            for r0 in range(1, p):
                v = sum(cps[i] * pow(r0, i, p) for i in range(7)) % p
                if v == 0:
                    x_ = r0
                    for _ in range(4):
                        fx = sum(cps[i] * pow(x_, i, mod3) for i in range(7)) % mod3
                        fpx = sum(i * cps[i] * pow(x_, i - 1, mod3) for i in range(1, 7)) % mod3
                        if fpx % p == 0:
                            break
                        x_ = (x_ - fx * pow(fpx, -1, mod3)) % mod3
                    row['unit_root_mod_p3'] = int(x_)
                    break
        except Exception as ex:
            row['unit_root_err'] = str(ex)[:80]
        fibs[key] = row
    out['fibers'] = fibs
    out['wall_s'] = round(time.time() - tw0, 1)
    return out, tr, F, gb


tA0 = time.time()
try:
    passA, trA, FA, gbA = one_pass(1.0, 'full')
except AssertionError as ex:
    passA, trA = {'tag': 'full', 'WALL': str(ex.args)[:400],
                  'wall_s': round(time.time() - tA0, 1)}, None
passA['rss_max_kb'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
rec['pass_A'] = passA
print('pass A done %.1fs rss=%dMB GE=%s' % (passA['wall_s'],
      passA['rss_max_kb'] // 1024,
      passA.get('GE_gauge', {}).get('mode', passA.get('WALL', '?')[:40])), flush=True)

k2_stop = passA['wall_s'] > BAND['stop_trigger_s']['pass_A']
k3_stop = passA.get('E_OBJECT', '').startswith('REFUSED(K3')
if k2_stop:
    rec['K2_stop'] = {'fired': True, 'pass_A_wall_s': passA['wall_s'],
                      'trigger_s': BAND['stop_trigger_s']['pass_A'],
                      'action': 'pass B not fired; re-price receipt owed'}
elif k3_stop:
    rec['pass_B'] = {'skipped_by': 'K3-MOVING in pass A (design wall; premise dead, '
                                   'no further spend priced for that branch)'}
elif trA is None:
    rec['pass_B'] = {'skipped_by': 'pass A WALL (assert) — see pass_A.WALL'}
else:
    tB0 = time.time()
    try:
        passB, trB, FB, gbB = one_pass(0.8, 'shorter-truncation')
    except AssertionError as ex:
        passB = {'tag': 'shorter-truncation', 'WALL': str(ex.args)[:400],
                 'wall_s': round(time.time() - tB0, 1)}
    rec['pass_B'] = passB
    print('pass B done %.1fs GE=%s' % (passB['wall_s'],
          passB.get('GE_gauge', {}).get('mode', passB.get('WALL', '?')[:40])), flush=True)

# ---- gates (verbatim CAL-D where reachable)
gates = {}
if trA is not None:
    sols = trA['mr']['sols']
    deff = trA['sb']['deff']
    j7 = deff.index(7)
    sigma = trA['sigma']
    mod2 = p ** (sigma + 6)
    c1rows = {}
    for x0, key in ((2, '2'), (3, '3'), ((1, 2), '1/2')):
        xh = dk.teich(x0, p, sigma + 6)
        B = 0
        po = 1
        for m in range(0, 7 + p):
            v = sols[j7][m][0] if sols[j7][m] else 0
            B = (B + v * po) % mod2
            po = po * xh % mod2
        Tt = 0
        po = 1
        for m in range(0, 7 + p * p):
            v = sols[j7][m][0] if sols[j7][m] else 0
            Tt = (Tt + v * po) % mod2
            po = po * xh % mod2
        vB = dk.valp(B, p) or 0
        u = (Tt // p ** vB) * pow(B // p ** vB, -1, p * p) % (p * p)
        bu = grid61.get(key, {}).get('u')
        c1rows[key] = {'my_level1_ratio': int(u), 'comparator_u': bu,
                       'MATCH': bool(bu is not None and u == bu)}
else:
    c1rows = {}
gates['C1prime_comparator_reproduction'] = c1rows
gates['C1prime_pass'] = bool(c1rows) and all(v['MATCH'] for v in c1rows.values()
                                             if v['comparator_u'] is not None)
pc = lk.pcurv_jordan(p, 2)
gates['C3_pcurv_pilot_fiber'] = {'jordan_type': pc.get('jordan_type'),
                                 'profile': str(pc.get('profile'))}
_f1dir = os.environ.get('DEFORM_F1_RECEIPTS')
if not _f1dir or not os.path.exists(os.path.join(_f1dir, 'FANF1_RECEIPT.json')):
    raise SystemExit('run_cal_d2: comparator input missing — set DEFORM_F1_RECEIPTS to a '
                     'directory holding FANF1_RECEIPT.json (not distributed with the repo)')
fanf1 = json.load(open(os.path.join(_f1dir, 'FANF1_RECEIPT.json')))
gates['C3_fanf1_majority'] = fanf1['scan']['61']['majority_jordan_type']
gates['C3_pass'] = bool(pc.get('jordan_type') == [5, 1] == fanf1['scan']['61']['majority_jordan_type']
                        and trA is not None and trA['sb']['jordan'] == [5, 1])
chi61 = dk.f2.chi_table(p)
ap_e = dk.f2.legendre_ap_fast(4, p, chi61)
ap2_e = ap_e * ap_e - 2 * p
gates['C2_elliptic_exact'] = {'ap': ap_e, 'ap2_identity': ap2_e,
                              'ap2_charsum': dk.f2.ap2_charsum(4, p, chi61),
                              'COPAIR': bool(ap2_e == dk.f2.ap2_charsum(4, p, chi61))}
c2rows = []
if 'fibers' in rec['pass_A']:
    for br in rec['pass_A']['fibers']['2']['branches']:
        iso = br['iso']
        if iso.get('Rp_integer'):
            cont = dk.elliptic_containment(iso['Rp_integer'], p, ap_e)
            c2rows.append({'sign': br['sign'], 'root_mod_p': br['root_mod_p'],
                           'containment': cont})
gates['C2_containment_rows'] = c2rows
gates['C2_some_branch_contains'] = any(r['containment']['S1_div_1_papt2_p4'] or
                                       r['containment']['S2_div_1_pap_p3'] or
                                       r['containment']['S2b_div_1_p2ap_p5']
                                       for r in c2rows)
haveB = 'pass_B' in rec and 'fibers' in rec.get('pass_B', {})
stab = {'unit_root_A': rec['pass_A'].get('fibers', {}).get('2', {}).get('unit_root_mod_p3'),
        'unit_root_B': rec['pass_B']['fibers']['2'].get('unit_root_mod_p3') if haveB else None}
stab['unit_root_stable_mod_p3'] = bool(stab['unit_root_A'] is not None and
                                       stab['unit_root_A'] == stab['unit_root_B'])
pairsA = {(b['sign'], b['root_mod_p']): b['iso'].get('Rp_integer')
          for b in rec['pass_A'].get('fibers', {}).get('2', {}).get('branches', [])}
pairsB = ({(b['sign'], b['root_mod_p']): b['iso'].get('Rp_integer')
           for b in rec['pass_B']['fibers']['2']['branches']} if haveB else {})
stab['isolated_integers_agree'] = bool(pairsA and pairsA == pairsB)
stab['some_branch_isolates'] = any(v for v in pairsA.values())
gates['stability'] = stab

# E-object bar (gauge-aware GE)
def e_ok(prow):
    ge = prow.get('GE_gauge', {})
    if ge.get('mode') == 'classic-ray-unique':
        return True
    return ge.get('mode') == 'gauge-well-defined' and (ge.get('d_star') or 0) >= 1
gates['E_object_ok_A'] = e_ok(rec['pass_A'])
gates['E_object_ok_B'] = e_ok(rec['pass_B']) if haveB else False
rec['gates'] = gates

clean = bool(gates['C1prime_pass'] and gates['C3_pass'] and
             rec['pass_A'].get('GD_pass') and (haveB and rec['pass_B'].get('GD_pass')) and
             gates['E_object_ok_A'] and gates['E_object_ok_B'] and
             stab['unit_root_stable_mod_p3'] and stab['isolated_integers_agree'] and
             stab['some_branch_isolates'] and gates['C2_some_branch_contains'] and
             not k2_stop)
rec['verdict'] = ('CAL-D2-CLEAN: L6 pilot (p=61, x0=2 + riders) at the (3,5,1,1) rate is '
                  'two-truncation STABLE with all structural gates green (GE %s/%s), C1\' exact, '
                  'C2 containment, C3 (5,1)' % (rec['pass_A'].get('GE_gauge', {}).get('mode'),
                                                rec.get('pass_B', {}).get('GE_gauge', {}).get('mode'))) \
                 if clean else 'CAL-D2-NOT-CLEAN: see gates (K1\' class: one debug cycle governs; license rides on this)'
rec['named_branch_p29_class'] = {
    'd_star_A': rec['pass_A'].get('pin', {}).get('depth_ladder_dstar'),
    'd_star_B': rec.get('pass_B', {}).get('pin', {}).get('depth_ladder_dstar') if haveB else None,
    'prior_recipe_d_star_at_2311': 0,
    'recurrence_semantics': 'd_star still 0 at the full rate after the one honest debug cycle -> '
                            'LICENSE-REFUSED-FINAL, deeper-NF-window named the live branch'}
rec['pricing_row'] = {'pass_A_wall_s': rec['pass_A'].get('wall_s'),
                      'pass_B_wall_s': rec.get('pass_B', {}).get('wall_s') if haveB else None,
                      'band': BAND,
                      'rss_max_kb_end': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                      'note': 'per-prime production cost class = pass walls at that prime\'s T-scaling '
                              '(projection law wall_p ~ pilot * (T_p/T_61)^1.7); RSS is the memory-width input'}
rec['wall_s_total'] = round(time.time() - t0, 1)
out = os.path.join(MEMBER_DIR, 'work', 'CAL_D2_RECEIPT.json')
lk_emit(rec, out, component=COMPONENT, script_path=os.path.abspath(__file__),
        member_dir=MEMBER_DIR)
dk.lint(out)
print('VERDICT:', rec['verdict'][:160])
