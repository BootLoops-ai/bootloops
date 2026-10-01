#!/usr/bin/env python3
"""selftest.py — ACCEPTANCE TEST: reproduce the reference-family Frobenius-
boundary numbers from its on-disk kira_targets_all.m.

Regression targets are the CONFIRMED reduced-spectrum counts
(independently verified by disjoint methods):
  degD=38, degG=40, poly-DE fit relerr <= 1e-70; 30 exponents at every eps;
  genuine nilpotent D_0 (rank 5), shear reduction in 3 steps;
  TRUE exponent multiset (eig of the reduced residue Dt1), classified:
    22 SURVIVORS  {1-eps x3, 1-2eps x6, 2-eps x2, 7/4+eps/2 x3,
                   9/4+eps/2 x5, 11/4+eps/2 x2, 7/2+2eps x1}
     5 FLAGGED    {3/2 x3, 5/2 x2}  (half-integer eps-indep — no recorded kill)
     3 KILLED     {2 x3}            (integer eps-indep, vacuum theorem)
  all semisimple; trace rules tr D1 = 69-10eps, tr Dt1 = 54-10eps;
  branch stage: the 12 reference branch series (requests labeled by eig(D1))
  still reproduce, no log-branch flag, Phi(196) 30x12 with FULL COLUMN RANK 12
  (column-scaled SVD; an argmax seed->column fetch collides here and drops
  the rank to 11 — seed_index is required) and an injective per-family
  seed->orbit-column map.
The naive expectation "12 survivors of 30" is WRONG as a counting claim;
the 12 branch series themselves are genuine.
Run:  python3 -m frobenius_boundary.selftest   (from tools/)   ~15-25 min.
"""
import json, os, sys, time
os.environ.setdefault('OPENBLAS_NUM_THREADS', '2')
import mpmath as mp
mp.mp.dps = 50      # set FIRST

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from frobenius_boundary import cli, io_util

# TRUE exponent multiset (reduced-residue frame, gauge as _reduce_fuchs emits)
EXP_SURV = {('1', '-1'): 3, ('1', '-2'): 6, ('2', '-1'): 2,
            ('7/4', '1/2'): 3, ('9/4', '1/2'): 5, ('11/4', '1/2'): 2,
            ('7/2', '2'): 1}
EXP_FLAG = {('3/2', '0'): 3, ('5/2', '0'): 2}
EXP_KILL = {('2', '0'): 3}
# the 12 reference branch requests (eig(D1) labels; series confirmed genuine)
EXP_BRANCH_FAMS = {'1 + (-1)*eps': 1, '1 + (-2)*eps': 6, '2 + (-1)*eps': 1,
                   '11/4 + (1/2)*eps': 2, '13/4 + (1/2)*eps': 1,
                   '9/2 + (2)*eps': 1}

def main():
    # FROB_CONFIG env override: the regression targets above belong to a
    # reference reduction table that does not ship with this repo; build your
    # own per-point Kira table, fill in a copy of config_example.json
    # (kira_targets_m + masters), and point FROB_CONFIG at it.
    cfg_path = os.environ.get('FROB_CONFIG',
                              os.path.join(HERE, 'config_example.json'))
    cfg = json.load(open(cfg_path))
    # The regression targets above include 3 KILLED rows of class
    # INT_EPS_INDEP__KILLED_vacuum_theorem, which classify_mod emits only when
    # kill_integer_eps_indep is true; config_example.json ships it false (the
    # kill needs a per-family vacuum check — never blanket).  Force it here so
    # a filled-in copy of config_example.json reproduces the targets.
    cfg['kill_integer_eps_indep'] = True
    if not os.path.exists(cfg.get('kira_targets_m', '')):
        print('SELFTEST DATA-GATED REFUSAL: the configured Kira reduction '
              'table does not exist:\n  kira_targets_m = '
              f"{cfg.get('kira_targets_m')}\nThis acceptance test needs a "
              'per-point Kira reduction table (kira2math output), which is '
              'not distributed with the repo. Fill in a copy of '
              'config_example.json and set FROB_CONFIG to it.')
        sys.exit(3)
    t0 = time.time()
    spec, br = cli.run(cfg)
    fails = []
    def chk(name, cond, detail=''):
        print(f"  [{'PASS' if cond else 'FAIL'}] {name} {detail}", flush=True)
        if not cond:
            fails.append(name)
    print('\n=== frobenius_boundary selftest vs the reference family ===')
    for e in cfg['eps_list']:
        pe = spec['per_eps'][e]
        chk(f'degD@{e}', pe['degD'] == 38, f"got {pe['degD']}")
        chk(f'degG@{e}', pe['degG'] == 40, f"got {pe['degG']}")
        chk(f'fit_err@{e}', pe['fit_err_log10'] <= -70,
            f"log10={pe['fit_err_log10']:.1f}")
        chk(f'n_eig@{e}', len(pe['eigvals']) == 30)
        chk(f'genuine D_0 detected@{e}', pe.get('fwd_blocks') == [0],
            f"got {pe.get('fwd_blocks')}")
        chk(f'D_0 rank 5@{e}', pe.get('D0_rank') == 5,
            f"got {pe.get('D0_rank')}")
        chk(f'shear reduction 3 steps@{e}', pe.get('reduce_jmax') == 3,
            f"got {pe.get('reduce_jmax')}")
        chk(f'n_eig_reduced@{e}',
            len(pe.get('eigvals_reduced') or []) == 30,
            f"got {len(pe.get('eigvals_reduced') or [])}")
    chk('classification frame = reduced residue',
        'reduced residue' in spec.get('classification_frame', ''),
        spec.get('classification_frame', 'MISSING'))
    surv = [g for g in spec['groups'] if 'SURVIVES' in g['class']]
    flag = [g for g in spec['groups'] if 'FLAG' in g['class']]
    kill = [g for g in spec['groups'] if 'KILLED' in g['class']]
    chk('survivor families == reduced-audit list',
        {(g['a'], g['b']): g['mult'] for g in surv} == EXP_SURV,
        f"got {sorted((g['a'], g['b'], g['mult']) for g in surv)}")
    chk('flagged families == half-int eps-indep {3/2 x3, 5/2 x2}',
        {(g['a'], g['b']): g['mult'] for g in flag} == EXP_FLAG,
        f"got {sorted((g['a'], g['b'], g['mult']) for g in flag)}")
    chk('killed == integer eps-indep {2 x3} (vacuum theorem)',
        {(g['a'], g['b']): g['mult'] for g in kill} == EXP_KILL
        and all('vacuum' in g['class'] for g in kill),
        f"got {sorted((g['a'], g['b'], g['mult'], g['class']) for g in kill)}")
    chk('n_free_constants == 22 (NOT the refuted 12)',
        spec['n_free_constants'] == 22, f"got {spec['n_free_constants']}")
    chk('n_flagged == 5', spec.get('n_flagged') == 5,
        f"got {spec.get('n_flagged')}")
    chk('n_killed == 3', spec['n_killed'] == 3, f"got {spec['n_killed']}")
    chk('22 + 5 + 3 == 30 (full local space)',
        spec['n_free_constants'] + spec.get('n_flagged', 0)
        + spec['n_killed'] == 30)
    for g in spec['groups']:
        chk(f"verify {g['a']}+({g['b']})eps <= 1e-24",
            g['verify_maxres_alleps'] <= 1e-24,
            f"res={g['verify_maxres_alleps']:.1e}")
        chk(f"semisimple {g['a']}+({g['b']})eps",
            g['n_jordan_blocks_geo'] == g['mult'],
            f"geo={g['n_jordan_blocks_geo']} mult={g['mult']}")
    chk('trace rule = 69 + (-10)*eps', spec['trace_rule']['tr_D1'] == '69 + (-10)*eps'
        and spec['trace_rule']['pass'], str(spec['trace_rule']))
    chk('reduced trace rule = 54 + (-10)*eps',
        spec.get('trace_rule_reduced', {}).get('tr_D1') == '54 + (-10)*eps'
        and spec.get('trace_rule_reduced', {}).get('pass', False),
        str(spec.get('trace_rule_reduced')))
    chk('no unclassified exponents (exclusion rows empty: all 30 rational)',
        len(spec['exclusion_rows']) == 0, f"got {len(spec['exclusion_rows'])}")
    chk('no exclusion flags', len(spec['exclusion_flags']) == 0)
    # branch stage (requests labeled by eig(D1); the 12 reference genuine series)
    chk('branch stage ran', br is not None)
    if br:
        chk('12 branch series', br['n_branches'] == 12, f"got {br['n_branches']}")
        fams_got = {}
        for b in br['branches']:
            fams_got[b['family']] = fams_got.get(b['family'], 0) + 1
        chk('branch families == reference 12-request list',
            fams_got == EXP_BRANCH_FAMS, f"got {sorted(fams_got.items())}")
        chk('no log-branch flag',
            not any(b['log_branch_required'] for b in br['branches']))
        chk('branch refine residuals small',
            all(b['refine_rfinal'] is None or b['refine_rfinal'] < 1e-25
                for b in br['branches']),
            f"max={max((b['refine_rfinal'] or 0) for b in br['branches']):.1e}")
        chk('seed residuals small',
            all(b['seed_resid'] < 1e-12 for b in br['branches']),
            f"max={max(b['seed_resid'] for b in br['branches']):.1e}")
        # guard: seed->orbit-column map must be injective per family
        used = {}
        for b in br['branches']:
            used.setdefault(b['family'], []).append(b.get('seed_col_used'))
        inj = all(None not in cols and len(set(cols)) == len(cols)
                  for cols in used.values())
        chk('seed->column map injective per family (seed_col_used)', inj,
            f"got {sorted(used.items())}")
        phi_ok = all(mp.isfinite(mp.mpf(re)) and mp.isfinite(mp.mpf(im))
                     for row in br['Phi'] for re, im in row)
        chk('Phi(196) finite 30x12', phi_ok and len(br['Phi']) == 30
            and len(br['Phi'][0]) == 12)
        # guard: rank(Phi) must equal the column count — counts INDEPENDENT
        # columns, not requested series (a duplicate column would pass a
        # naive count: e.g. a duplicate 1-2eps column with column-scaled
        # sv 9.4e-63 vs 5.9e-5).
        if phi_ok:
            save = mp.mp.dps; mp.mp.dps = 60
            Phi = io_util.s2mat(br['Phi'])
            for j in range(Phi.cols):
                m = max(abs(Phi[i, j]) for i in range(Phi.rows))
                if m > 0:
                    for i in range(Phi.rows):
                        Phi[i, j] /= m
            S = mp.svd_c(Phi, compute_uv=False)
            rank = sum(1 for i in range(len(S))
                       if S[i] > S[0] * mp.mpf(10) ** -30)
            mp.mp.dps = save
            chk('rank(Phi) == n_branches (column-scaled SVD, rel 1e-30)',
                rank == br['n_branches'],
                f"rank={rank} smin/smax={mp.nstr(S[len(S)-1] / S[0], 3)}")
    verdict = 'PASS' if not fails else f'FAIL ({len(fails)}: {fails})'
    print(f'\n=== SELFTEST {verdict}  ({time.time()-t0:.0f}s) ===')
    json.dump({'verdict': verdict, 'fails': fails,
               'elapsed_s': round(time.time() - t0, 1)},
              open(os.path.join(cfg['out_dir'], 'SELFTEST_VERDICT.json'), 'w'),
              indent=1)
    sys.exit(0 if not fails else 1)

if __name__ == '__main__':
    main()
