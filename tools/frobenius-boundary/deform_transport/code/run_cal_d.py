#!/usr/bin/env python3
"""run_cal_d.py — CAL-D (license gate 2): the L6-class operator at the
declared pilot (p=61, x0=2 + riders 3, 1/2).  Matrix-transport STABILITY
CLEAN = two-truncation agreement of the isolated integers + all structural
gates green + C1' exact + C2 containment + C3 (5,1)/Newton.
usage: run_cal_d.py [K] [mults e.g. 3,5,1,1]
Comparator inputs (NOT distributed with the repo; refuses loudly without
them): DEFORM_R2_GRID = a level-1 unit-root grid JSON ({'grid': {'61':
{fiber: {'u': int}}}}); DEFORM_F1_RECEIPTS = a directory with
FANF1_RECEIPT.json (majority Jordan type per prime).  FROB_L6_OPERATOR
selects the operator.  Emits work/CAL_D_RECEIPT.json; prints the measured
per-prime cost row (the production pricing input)."""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import deform_kit as dk
import lfact_kit as lk

t0 = time.time()
p = 61
K = int(sys.argv[1]) if len(sys.argv) > 1 else 9
mults = [int(x) for x in sys.argv[2].split(',')] if len(sys.argv) > 2 else [3, 5, 1, 1]
polys = [([1, -1], '1-x'), ([1, 1], '1+x'), ([1, 1, 1], 'Phi3'), ([1, 1, 1, 1, 1], 'Phi5')]
twist = [(c, nm, m) for (c, nm), m in zip(polys, mults) if m > 0]
rec = {'producer': dk.producer('code/run_cal_d.py'),
       'spec': 'CAL-D (declared pilot law: p=61, x0=2, riders 3, 1/2)',
       'register': 'license gate 2; every verdict script-emitted',
       'pilot': {'p': p, 'x0': 2, 'riders': ['3', '1/2'], 'K_target': K,
                 'twist_mults': mults}}
op = dk.op_L6()

def _env_json(var, what):
    path = os.environ.get(var)
    if not path or not os.path.exists(path):
        raise SystemExit('run_cal_d: comparator input missing — set %s to %s '
                         '(not distributed with the repo)' % (var, what))
    return json.load(open(path))

grid61 = _env_json('DEFORM_R2_GRID', 'a level-1 unit-root grid JSON')['grid']['61']

def one_pass(NF_scale, tag):
    out = {'tag': tag}
    tw0 = time.time()
    NF = None
    if NF_scale != 1.0:
        dtarget = K + 6
        s_twist = dtarget + 2
        degD = sum((len(c) - 1) * m for (c, nm, m) in twist)
        NF = int(p * (s_twist * degD + 3.0 * dtarget + 8) * NF_scale)
    tr = dk.transport(op, p, K, twist=twist, NF=NF)
    r = tr['res']
    out['transport'] = {k: r[k] for k in ('NF', 'Ntot', 'jordan', 'nilpotent',
                                          'E_space_dim', 'inv_residual_certified_depth')}
    out['march'] = r['march']
    out['pin'] = r['pin']
    out['twist_law'] = r['twist']
    E0, evec = dk.pinned_E(tr)
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
        # unit root mod p^3 (slope-0 eigenvalue of Fm): from charpoly via Newton
        # iteration on the unit root of cp: x satisfies cp(x)=0 with x unit
        row['unit_root_mod_p3'] = None
        try:
            mod3 = p ** min(3, cert)
            # unit root = root of charpoly that is a unit: Hensel from the
            # mod-p root of cp with unit value
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

passA, trA, FA, gbA = one_pass(1.0, 'full')
rec['pass_A'] = passA
print('pass A done %.1fs' % passA['wall_s'], flush=True)
passB, trB, FB, gbB = one_pass(0.8, 'shorter-truncation')
rec['pass_B'] = passB
print('pass B done %.1fs' % passB['wall_s'], flush=True)

# ---- gates
gates = {}
# C1': comparator level-1 ratio reproduction from pass-A's own series
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
gates['C1prime_comparator_reproduction'] = c1rows
gates['C1prime_pass'] = all(v['MATCH'] for v in c1rows.values() if v['comparator_u'] is not None)
# C3: independent p-curvature + exhaustive-scan majority row
pc = lk.pcurv_jordan(p, 2)
gates['C3_pcurv_pilot_fiber'] = {'jordan_type': pc.get('jordan_type'),
                                 'profile': pc.get('profile')}
_f1dir = os.environ.get('DEFORM_F1_RECEIPTS')
if not _f1dir or not os.path.exists(os.path.join(_f1dir, 'FANF1_RECEIPT.json')):
    raise SystemExit('run_cal_d: comparator input missing — set DEFORM_F1_RECEIPTS to a '
                     'directory holding FANF1_RECEIPT.json (not distributed with the repo)')
fanf1 = json.load(open(os.path.join(_f1dir, 'FANF1_RECEIPT.json')))
gates['C3_fanf1_majority'] = fanf1['scan']['61']['majority_jordan_type']
gates['C3_pass'] = bool(pc.get('jordan_type') == [5, 1] == fanf1['scan']['61']['majority_jordan_type']
                        and trA['sb']['jordan'] == [5, 1])
# C2: elliptic containment + declared sign adjudication (on any branch with an isolated factor)
chi61 = dk.f2.chi_table(p)
ap_e = dk.f2.legendre_ap_fast(4, p, chi61)
ap2_e = ap_e * ap_e - 2 * p
gates['C2_elliptic_exact'] = {'ap': ap_e, 'ap2_identity': ap2_e,
                              'ap2_charsum': dk.f2.ap2_charsum(4, p, chi61),
                              'COPAIR': bool(ap2_e == dk.f2.ap2_charsum(4, p, chi61))}
c2rows = []
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
# STABILITY: two-truncation agreement of isolated integers (per branch tag) + unit root mod p^3
stab = {'unit_root_A': rec['pass_A']['fibers']['2'].get('unit_root_mod_p3'),
        'unit_root_B': rec['pass_B']['fibers']['2'].get('unit_root_mod_p3')}
stab['unit_root_stable_mod_p3'] = bool(stab['unit_root_A'] is not None and
                                       stab['unit_root_A'] == stab['unit_root_B'])
pairsA = {(b['sign'], b['root_mod_p']): b['iso'].get('Rp_integer')
          for b in rec['pass_A']['fibers']['2']['branches']}
pairsB = {(b['sign'], b['root_mod_p']): b['iso'].get('Rp_integer')
          for b in rec['pass_B']['fibers']['2']['branches']}
stab['isolated_integers_agree'] = bool(pairsA and pairsA == pairsB)
stab['some_branch_isolates'] = any(v for v in pairsA.values())
gates['stability'] = stab
rec['gates'] = gates
clean = bool(gates['C1prime_pass'] and gates['C3_pass'] and
             rec['pass_A']['GD_pass'] and rec['pass_B']['GD_pass'] and
             rec['pass_A']['pin']['ray_unique_up_to_scalar'] and
             stab['unit_root_stable_mod_p3'] and stab['isolated_integers_agree'] and
             stab['some_branch_isolates'] and gates['C2_some_branch_contains'])
rec['verdict'] = ('CAL-D-CLEAN: L6 pilot (p=61, x0=2 + riders) matrix transport is two-truncation '
                  'STABLE with all structural gates green, C1\' exact, C2 containment, C3 (5,1)') \
                 if clean else 'CAL-D-NOT-CLEAN: see gates (K1 class: one debug cycle governs; license rides on this)'
rec['pricing_row'] = {'pass_A_wall_s': rec['pass_A']['wall_s'],
                      'pass_B_wall_s': rec['pass_B']['wall_s'],
                      'note': 'per-prime production cost class = pass_A + pass_B at that prime\'s T-scaling'}
rec['wall_s_total'] = round(time.time() - t0, 1)
dk.emit(rec, 'work/CAL_D_RECEIPT.json')
print('VERDICT:', rec['verdict'][:130])
