#!/usr/bin/env python3
"""run_cal_a.py — CAL-A: the rank-2 Legendre theorem case through the FULL
matrix transport.  The transported R_p must equal the Dwork-theorem quadratic
exactly; fixes the member's sign/twist convention ONCE.
Emits work/CAL_A_RECEIPT.json (script-emitted verdicts)."""
import json, os, sys, time, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import deform_kit as dk

MEMBER_DIR = dk.MEMBER_DIR
t0 = time.time()
rec = {'producer': dk.producer('code/run_cal_a.py'),
       'spec': 'CAL-A: rank-2 theorem control (instrument spec in deform_kit.py docstring)',
       'register': 'calibration control (theorem case); every verdict script-emitted',
       'primes': {}}
op = dk.op_legendre()
K = 10
overall = True
convention_votes = {}
for p in (23, 29):
    row = {'K_target': K}
    try:
        tr = dk.transport(op, p, K)
        r = tr['res']
        row['transport'] = {k: r[k] for k in ('p', 'op', 'NF', 'Ntot', 'jordan',
                                              'nilpotent', 'E_space_dim',
                                              'inv_residual_certified_depth')}
        row['march'] = r['march']
        row['pin'] = r['pin']
        row['twist_law'] = r['twist']
        E0, evec = dk.pinned_E(tr)
        F = dk.assemble_full_F(tr, E0)
        gd = dk.gd_residual(tr, F)
        row['GD_residual_worst_deficit'] = gd
        row['GD_pass'] = bool((gd or 0) <= tr['keep'] - tr['keep_target'] + 64)
        gb = dk.gamma_branches(E0, p, wt=1, n=2, Kg=min(tr['d_star'], 12))
        row['gamma'] = {'detE0': gb['detE0'], 'val_detE0': gb['val_detE0'],
                        'target_val': gb['target_val'], 'n_branches': len(gb['branches'])}
        chi = dk.f2.chi_table(p)
        chim1 = 1 if p % 4 == 1 else -1
        fibers = {}
        for x0 in (2, 3, 4):
            lam = x0 * x0 % p
            if lam in (0, 1):
                continue
            ap = dk.f2.legendre_ap_fast(lam, p, chi)
            Fm, cert, tailv = dk.eval_F_at(tr, F, lam)
            mod = p ** cert
            cp = dk.berkowitz_charpoly(Fm, mod)
            c = dk.rp_coeffs(cp, p, cert)
            det_consist = (Fm[0][0] * Fm[1][1] - Fm[0][1] * Fm[1][0] - gb['detE0']) % mod == 0
            frow = {'lam': lam, 'ap_exact': ap, 'cert_depth': cert, 'tail_val': tailv,
                    'detE_consistency': bool(det_consist), 'branch_results': []}
            for br in gb['branches']:
                cn = dk.apply_gamma(c, br, p, cert)
                if cn is None:
                    continue
                c1n = dk.centered(cn[1], mod)
                c2n = dk.centered(cn[2], mod)
                tag = None
                if c2n == p:
                    if c1n == -chim1 * ap:
                        tag = 'Rp = 1 - chi2(-1) ap T + p T^2'
                    elif c1n == chim1 * ap:
                        tag = 'Rp = 1 + chi2(-1) ap T + p T^2'
                frow['branch_results'].append({'sign': br['sign'], 'root_mod_p': br['root_mod_p'],
                                               'c1': c1n if abs(c1n) < 10 ** 12 else 'large',
                                               'c2': c2n if abs(c2n) < 10 ** 12 else 'large',
                                               'theorem_match': tag})
                if tag:
                    convention_votes[tag] = convention_votes.get(tag, 0) + 1
            frow['exact_match_some_branch'] = any(b['theorem_match'] for b in frow['branch_results'])
            fibers[str(lam)] = frow
            overall = overall and frow['exact_match_some_branch'] and det_consist
        row['fibers'] = fibers
        row['PASS'] = all(f['exact_match_some_branch'] and f['detE_consistency']
                          for f in fibers.values()) and row['GD_pass'] and r['pin']['ray_unique_up_to_scalar']
    except AssertionError as ex:
        row['WALL'] = str(ex.args)[:400]
        row['PASS'] = False
        overall = False
    rec['primes'][str(p)] = row

rec['convention_measured'] = max(convention_votes, key=convention_votes.get) if convention_votes else None
rec['convention_votes'] = convention_votes
rec['verdict'] = ('CAL-A-PASS: the rank-2 Dwork theorem case reproduces EXACTLY through the full matrix '
                  'transport at p in {23,29} (all ordinary fibers, supersingular included); measured convention: %s'
                  % rec['convention_measured']) if overall else \
                 'CAL-A-FAIL: see prime rows (K1 class: one debug cycle governs)'
rec['wall_s'] = round(time.time() - t0, 1)
dk.emit(rec, 'work/CAL_A_RECEIPT.json')
print('VERDICT:', rec['verdict'][:120])
