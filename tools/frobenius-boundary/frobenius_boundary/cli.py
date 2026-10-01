#!/usr/bin/env python3
"""cli.py — end-to-end Frobenius-boundary run from a config json.

Usage: python3 -m frobenius_boundary.cli config.json
Config keys (see config_example.json):
  kira_targets_m, masters, family, var, mass_slots, n_loops, s_t_numeric,
  eps_list, dps, out_dir, tag, branch_n_terms, Mstart, kill_integer_eps_indep,
  ncpu, dps_eig
Writes <out_dir>/SPECTRUM_<tag>.json and BRANCHES_<tag>.json.
Classification frame: when a genuine nilpotent D_0 forward block
is detected, the counting classification (groups / n_free_constants / n_killed
/ n_flagged / exclusions) runs on eig of the SHEAR-REDUCED residue Dt1 (the
true exponent multiset; see spec_json['classification_frame']); the eig(D1)
classification is kept as 'groups_eig_D1' and still labels the branch stage's
requests.  Fuchsian systems (no genuine forward blocks): unchanged, eig(D1).
"""
import json, os, sys, time
os.environ.setdefault('OPENBLAS_NUM_THREADS', '2')
import mpmath as mp
mp.mp.dps = max(mp.mp.dps, 50)      # set FIRST
import multiprocessing as mpc

from . import core, de_core, poly_de, infinity, classify_mod, branches_mod, io_util

_CFG = {}


def _worker(eps):
    cfg = _CFG
    mp.mp.dps = cfg['dps'] + 30
    t0 = time.time()
    dd_rat, dd, epsv = core.eps_to_d(eps)
    basis = de_core.load_masters(cfg['masters'])
    ai = de_core.alpha_int(basis)
    r = poly_de.build_poly_DE_full(cfg['kira_targets_m'], basis, dd_rat,
                                   cfg['dps'] + 30, cfg['mass_slots'],
                                   var=cfg.get('var', 'm2'),
                                   family=cfg.get('family'),
                                   verbose=(eps == cfg['eps_list'][0]))
    N = r['N']
    Rj = infinity.laurent_xA(r['Dc'], r['Gc'], N, 10, cfg['dps'] + 70)
    Delta = core.delta_of(ai)
    D1 = infinity.D1_matrix(Rj, ai, dd, Delta, N, cfg.get('n_loops', 3))
    ev = infinity.eig_spectrum(D1, cfg.get('dps_eig', 60))
    out = {'eps': eps, 'degD': r['degD'], 'degG': r['degG'],
           'fit_err_log10': r['fit_err_log10'],
           'D1': io_util.mat2s(D1, 45),
           'eigvals': [(mp.nstr(z.real, 30), mp.nstr(mp.im(z), 30)) for z in ev]}
    # Reduced-residue spectrum:
    # with a genuine nilpotent D_0 the N-frame is NOT Fuchsian at u=0 and
    # eig(D1) is NOT the exponent multiset — the true exponents are eig of the
    # residue of the shear-reduced frame (the same _reduce_fuchs the branch
    # stage uses).  Computed here per eps so classification can consume it.
    Dp = infinity.build_Dp(Rj, ai, dd, Delta, N, 10, cfg.get('n_loops', 3))
    fwd = branches_mod.genuine_fwd_blocks(Dp)
    out['fwd_blocks'] = fwd
    if fwd == [0]:
        scale = max(branches_mod._matnorm(Dp[p]) for p in (1, 2) if p in Dp)
        S0 = mp.svd_c(mp.matrix(Dp[0]), compute_uv=False)
        out['D0_rel'] = float(branches_mod._matnorm(Dp[0]) / scale)
        out['D0_rank'] = sum(1 for i in range(len(S0))
                             if S0[i] > S0[0] * mp.mpf(10) ** -30)
        red = branches_mod._reduce_fuchs(Dp, N, cfg['dps'] + 40)
        if red is not None:
            Dt, Tj = red
            evr = infinity.eig_spectrum(Dt[1], cfg.get('dps_eig', 60))
            out['eigvals_reduced'] = [(mp.nstr(z.real, 30), mp.nstr(mp.im(z), 30))
                                      for z in evr]
            out['Dt1'] = io_util.mat2s(Dt[1], 45)
            out['reduce_jmax'] = len(Tj) - 1
        else:
            out['reduce_failed'] = True
    out['t_s'] = round(time.time() - t0, 1)
    if eps == cfg['eps_list'][0]:        # keep D,G for the branch stage
        out['Dc'] = io_util.coeffs2s(r['Dc'], cfg['dps'] + 30)
        out['Gc'] = [io_util.mat2s(g, cfg['dps'] + 30) for g in r['Gc']]
    return out


def run(cfg):
    global _CFG
    _CFG = cfg
    t0 = time.time()
    out_dir = cfg.get('out_dir', '.'); tag = cfg.get('tag', 'run')
    os.makedirs(out_dir, exist_ok=True)
    basis = de_core.load_masters(cfg['masters'])
    ai = de_core.alpha_int(basis)
    # pre-warm caches in parent so fork()ed workers inherit (168MB parse once)
    poly_de.get_A(cfg['kira_targets_m'], basis, cfg['mass_slots'], cfg.get('family'))
    poly_de.get_denom_strings(cfg['kira_targets_m'])
    eps_list = cfg['eps_list']
    ncpu = min(cfg.get('ncpu', 4), len(eps_list))
    if ncpu > 1:
        with mpc.get_context('fork').Pool(ncpu) as pool:
            specs = pool.map(_worker, eps_list)
    else:
        specs = [_worker(e) for e in eps_list]
    specs = {s['eps']: s for s in specs}
    mp.mp.dps = cfg['dps'] + 30
    eps_vals = [core.eps_to_d(e)[2] for e in eps_list]
    EVs = [[mp.mpc(mp.mpf(re), mp.mpf(im)) for re, im in specs[e]['eigvals']]
           for e in eps_list]
    D1s = [io_util.s2mat(specs[e]['D1']) for e in eps_list]
    fams, groups = classify_mod.classify(
        EVs, eps_vals, kill_integer_eps_indep=cfg.get('kill_integer_eps_indep', True))
    groups = classify_mod.jordan_info(D1s[0], groups, eps_vals[0])
    tr = classify_mod.trace_rule(D1s, eps_vals)
    strata = classify_mod.exclude_strata(groups, EVs, eps_vals)
    # ---- classification frame: with a genuine nilpotent D_0 the
    # counting truth is the shear-reduced residue spectrum; eig(D1) keeps only
    # the role of branch-request labels (the constructed series are unchanged).
    reduced = all(specs[e].get('eigvals_reduced') for e in eps_list)
    warn = None
    if reduced:
        EVr = [[mp.mpc(mp.mpf(re), mp.mpf(im))
                for re, im in specs[e]['eigvals_reduced']] for e in eps_list]
        Dt1s = [io_util.s2mat(specs[e]['Dt1']) for e in eps_list]
        fams_r, groups_r = classify_mod.classify(
            EVr, eps_vals,
            kill_integer_eps_indep=cfg.get('kill_integer_eps_indep', True))
        groups_r = classify_mod.jordan_info(Dt1s[0], groups_r, eps_vals[0])
        tr_r = classify_mod.trace_rule(Dt1s, eps_vals)
        strata_r = classify_mod.exclude_strata(groups_r, EVr, eps_vals)
        groups_active, strata_active = groups_r, strata_r
        frame = ('shear-reduced residue Dt1 (genuine nilpotent D_0: '
                 'true exponent multiset; eig(D1) is NOT a counting basis)')
    else:
        groups_active, strata_active = groups, strata
        frame = 'D1 (Fuchsian N-frame: exponents = eig(D1))'
        bad = [e for e in eps_list if specs[e].get('fwd_blocks')]
        if bad:
            warn = (f'genuine forward blocks at eps {bad} not shear-reduced '
                    f'(reduce_failed or fdepth>1): eig(D1) counts are '
                    f'UNRELIABLE for this system — do not reuse as a '
                    f'counting claim')
    spec_json = {
        'config': {k: v for k, v in cfg.items()}, 'date': time.strftime('%F %T'),
        'convention': ('M_k ~ var^(L*d/2 + ai_k + lambda_j), lambda_j = eig(D1)'
                       if not reduced else
                       'M_k ~ var^(L*d/2 + ai_k + lambda_j); counting spectrum '
                       '= eig(shear-reduced residue Dt1); branch requests are '
                       'still labeled by eig(D1)'),
        'classification_frame': frame,
        'per_eps': {e: {k: specs[e][k] for k in
                        ('degD', 'degG', 'fit_err_log10', 'eigvals', 't_s',
                         'fwd_blocks', 'D0_rel', 'D0_rank', 'reduce_jmax',
                         'eigvals_reduced', 'reduce_failed')
                        if k in specs[e]}
                    for e in eps_list},
        'D1_first_eps': specs[eps_list[0]]['D1'],
        'trace_rule': tr,
        'groups': [{k: v for k, v in g.items() if k not in ('a_num', 'b_num')}
                   for g in groups_active],
        'n_free_constants': strata_active['n_free_constants'],
        'n_killed': strata_active['n_killed'],
        'n_flagged': sum(g['mult'] for g in strata_active['flagged']),
        'exclusion_rows': strata_active['exclusion_rows'],
        'exclusion_flags': strata_active['exclusion_flags'],
        'elapsed_s': round(time.time() - t0, 1)}
    if reduced:
        spec_json['Dt1_first_eps'] = specs[eps_list[0]]['Dt1']
        spec_json['trace_rule_reduced'] = tr_r
        spec_json['groups_eig_D1'] = [
            {k: v for k, v in g.items() if k not in ('a_num', 'b_num')}
            for g in groups]
    if warn:
        spec_json['classification_warning'] = warn
    fn = os.path.join(out_dir, f'SPECTRUM_{tag}.json')
    json.dump(spec_json, open(fn, 'w'), indent=1, default=str)
    print(f'[cli] wrote {fn} ({spec_json["elapsed_s"]}s)', flush=True)
    # ---- branch stage at eps_list[0]
    nb = cfg.get('branch_n_terms', 0)
    if not nb:
        return spec_json, None
    e0 = eps_list[0]
    _, dd, epsv = core.eps_to_d(e0)
    Dc = io_util.s2coeffs(specs[e0]['Dc'])
    Gc = [io_util.s2mat(g) for g in specs[e0]['Gc']]
    N = len(ai); Delta = core.delta_of(ai)
    prec = cfg['dps'] + 40
    Rj = infinity.laurent_xA(Dc, Gc, N, nb, prec)
    Dp = infinity.build_Dp(Rj, ai, dd, Delta, N, nb, cfg.get('n_loops', 3))
    D1 = D1s[0]
    branches, blist = [], []
    # Branch requests deliberately stay in the eig(D1) frame (strata, not
    # strata_active): the request lam must be an eig(D1) value for the seed /
    # orbit machinery, and the constructed series are confirmed genuine
    # (independent audit).  The COUNTING outputs above use strata_active.
    # seed_index=c makes the seed->orbit-column fetch injective.
    for g in strata['survivors']:
        lam = classify_mod._lam_of_group(g, eps_vals[0])
        V, ndim, sres = infinity.seed_vectors(D1, lam)
        for c in range(ndim):
            seed = mp.matrix([V[r, c] for r in range(N)])
            Nn, diag = branches_mod.branch_series(Dp, N, lam, seed, nb, prec,
                                                  eigvals=EVs[0], seed_index=c)
            blist.append((lam, Nn))
            branches.append({'family': f"{g['a']} + ({g['b']})*eps",
                             'lambda': mp.nstr(lam, 30), 'seed_col': c,
                             'seed_col_used': diag.get('seed_col_used'),
                             'seed_dim': ndim, 'seed_resid': sres,
                             'resonances': diag['resonances'],
                             'refine_rfinal': diag.get('refine_rfinal'),
                             'log_branch_required': diag.get('log_branch_required', False),
                             'Nn_norm_profile':
                                 [float(max(abs(Nn[n][i, 0]) for i in range(N)))
                                  for n in range(nb + 1)],
                             'Nn': {str(n): io_util.vec2s(
                                 [Nn[n][i, 0] for i in range(N)], cfg['dps'])
                                 for n in range(nb + 1)}})
    br_json = {'eps': e0, 'n_terms': nb, 'n_branches': len(blist),
               'branches': branches}
    if cfg.get('Mstart'):
        Phi = branches_mod.phi_matrix(mp.mpf(str(cfg['Mstart'])), blist, ai, dd,
                                      cfg.get('n_loops', 3))
        br_json['Mstart'] = cfg['Mstart']
        br_json['Phi'] = io_util.mat2s(Phi, cfg['dps'])
    fn2 = os.path.join(out_dir, f'BRANCHES_{tag}.json')
    json.dump(br_json, open(fn2, 'w'), indent=1, default=str)
    print(f'[cli] wrote {fn2} ({len(blist)} branch solutions)', flush=True)
    return spec_json, br_json


if __name__ == '__main__':
    cfg = json.load(open(sys.argv[1]))
    run(cfg)
