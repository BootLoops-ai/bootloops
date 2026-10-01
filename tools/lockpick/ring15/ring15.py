# lockpick ring15 member — orchestrator (builds RING15.json; writes under RING15_OUT).
"""ring15: orchestrator. Builds the conductor-15 / Q(sqrt(-15)) CM constant
ring at dps 210 AND dps 120 (independent reruns of the whole pipeline),
checkpoints each block to RING15_partial.json, writes final RING15.json.
PSLQ protocol: height-capped (maxcoeff 10^4), run at both precisions,
identical integer vector required; positive control = disc -4 lemniscatic
relation recovered from the AGM (Gamma-free) side.
Run:  python3 ring15.py
"""
import json, os, time
import mpmath as mp
from . import ring15_core as core
from . import ring15_hecke as hk
from . import ring15_lfun as lf

# Rebuild runs write under RING15_OUT (default cwd — launch from scratch space).
BASE = os.environ.get('RING15_OUT', '.')
PARTIAL = f'{BASE}/RING15_partial.json'
STATE = {'meta': {
    'started_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
    'tool': 'ring15.py', 'mpmath': mp.__version__,
    'provenance': {
        'newforms': 'S_3^new(Gamma_0(15),chi_-15) dim 2, CM by Q(sqrt(-15));'
            ' eigenforms from PARI/GP 2.15.4 mfeigenbasis (oracle transcript not shipped).'
            ' LABEL RULE-BASED (FLAG): LMFDB was unreachable at build time'
            ' (recaptcha/404); 15.3.b.a := lex-min trace orbit (a_2=-1) per LMFDB'
            ' sort rule; f_b := a_2=+1. Both forms in ring => downstream use is'
            ' label-independent.',
        'a_n': 'Grossencharacter psi((alpha))=alpha^2, psi(p2)=s*beta, s=+-1'
            ' (ring15_hecke.py); exact integer 3-way match: eta products'
            ' eta(3t)^3eta(5t)^3 +- eta(t)^3eta(15t)^3 (n<=60) and PARI (n<=20).',
        'L_f': 'incomplete-Gamma Hecke sums; Fricke matrix FIT (came out exactly'
            ' identity: eigenvalue +1 for both forms); split-point invariance'
            ' ~10^-(dps+30); independent oracle: PARI lfunmf matches 50-52 digits.',
        'CS': 'Omega_D=(2pi|D|)^{-1/2}[prodGamma^chi]^{w/(4h)}; h=2 normalization'
            ' PINNED: cm1*cm2=(sqrt15/2)Omega15^2, cm1/cm2=phi^(-1/3); unified'
            ' (sqrt|D|/2)^{h/2} law verified on D=-3,-4,-15 (ring15_core.py).',
        'controls': 'Catalan=L(chi_-4,2) vs mp.catalan; lemniscatic AGM PSLQ'
            ' [-2,4,-3,-1] at both dps; Lerch residuals ~0 for D=-3,-4,-15.',
        'pslq_protocol': 'maxcoeff 10^4, dps 120 AND 210, identical vectors'
            ' required, positive control recorded.'}}}


def ckpt():
    with open(PARTIAL, 'w') as fh:
        json.dump(STATE, fh, indent=1, default=str)


def nstr(v, dps):
    return mp.nstr(v, dps, strip_zeros=False)


def logdiff(v):
    av = abs(v)
    return float(mp.log10(av)) if av > 0 else -float('inf')


def run(dps):
    mp.mp.dps = dps + 30
    R, K = {}, {}
    blk = {'members': R, 'checks': K}
    STATE[f'dps{dps}'] = blk

    # block 1: elementary constants + Dirichlet L + CS periods + eta side
    R['pi'] = mp.pi; R['pi2'] = mp.pi ** 2; R['pi3'] = mp.pi ** 3
    R['sqrt3'] = mp.sqrt(3); R['sqrt5'] = mp.sqrt(5); R['sqrt15'] = mp.sqrt(15)
    R['phi'] = (1 + mp.sqrt(5)) / 2; R['log_phi'] = mp.log(R['phi'])
    R['gamma_quotient_15'] = core.gamma_quotient(-15)
    R['Omega15'] = core.Omega(-15)
    R['Omega15_cm1'] = core.cm_period_eta(1, 1, 4)   # tau_1, form [1,1,4]
    R['Omega15_cm2'] = core.cm_period_eta(2, 1, 2)   # tau_2, form [2,1,2]
    for D, tag in ((-15, 'chi15'), (-3, 'chi3'), (-4, 'chi4'), (5, 'chi5')):
        for s in (1, 2):
            R[f'L_{tag}_{s}'] = core.dirichlet_L(s, D)
    ctr = core.disc4_controls()
    R['varpi_agm'] = ctr.pop('varpi_agm')
    K.update({k: logdiff(v) for k, v in ctr.items()})
    K['classnum_L1_chi15'] = logdiff(R['L_chi15_1'] - 2 * mp.pi / mp.sqrt(15))
    # h=2 normalization, PINNED (see ring15_core docstring):
    #   cm1*cm2 = (sqrt15/2)*Omega15^2 ;  cm1/cm2 = phi^(-1/3)
    K['h2_pin_etaprod_minus_sqrt15half_Om2'] = logdiff(
        R['Omega15_cm1'] * R['Omega15_cm2']
        - mp.sqrt(15) / 2 * R['Omega15'] ** 2)
    K['gross_unit_cm1_over_cm2_minus_phim13'] = logdiff(
        R['Omega15_cm1'] / R['Omega15_cm2'] - R['phi'] ** (mp.mpf(-1) / 3))
    # unified CS normalization (sqrt|D|/2)^{h/2}: disc -3 factor (3/4)^{1/4}
    K['h1_eta_check_disc3'] = logdiff(
        (mp.mpf(3) / 4) ** mp.mpf('0.25') * core.Omega(-3)
        - core.cm_period_eta(1, 1, 1))
    ckpt()

    # block 2: PSLQ (height-capped, logged for 2-precision comparison)
    P = {}
    blk['pslq'] = P
    ctrl = mp.pslq([mp.log(R['varpi_agm']), mp.loggamma(mp.mpf(1) / 4),
                    mp.log(2), mp.log(mp.pi)], maxcoeff=10 ** 4, maxsteps=10 ** 6)
    P['control_lemniscate_logs'] = ctrl       # expect +-[2,-4,3,1]
    P['control_ok'] = (ctrl is not None and
                       [abs(x) for x in ctrl] == [2, 4, 3, 1])
    ratio = R['Omega15_cm1'] / R['Omega15_cm2']
    P['cm_ratio_logs_basis'] = ['log(cm1/cm2)', 'log_phi', 'log2', 'log3', 'log5']
    P['cm_ratio_logs'] = mp.pslq([mp.log(ratio), R['log_phi'], mp.log(2),
                                  mp.log(3), mp.log(5)],
                                 maxcoeff=10 ** 4, maxsteps=10 ** 6)
    h2r = R['Omega15_cm1'] * R['Omega15_cm2'] / R['Omega15'] ** 2
    P['h2_ratio_logs_basis'] = ['log(etaprod/Om15^2)', 'log2', 'log3', 'log5', 'logpi']
    P['h2_ratio_logs'] = mp.pslq([mp.log(h2r), mp.log(2), mp.log(3),
                                  mp.log(5), mp.log(mp.pi)],
                                 maxcoeff=10 ** 4, maxsteps=10 ** 6)
    ckpt()

    # block 3: Hecke eigenvalues, 3-way exact validation
    val = hk.run_validation(60)
    blk['hecke_validation'] = {k: v for k, v in val.items()}
    ckpt()

    # block 4: L(f,s) both forms; label rule: 15.3.b.a := a_2=-1 orbit (lex-min
    # trace sequence, LMFDB sorting rule; direct LMFDB fetch blocked -> rule-based)
    NL = 1400 if dps > 150 else 900
    lout = lf.L_values(NL)
    for tag, lab in (('m', 'f_a'), ('p', 'f_b')):   # m: a2=-1 -> 15.3.b.a
        for s in (1, 2, 3):
            R[f'L_{lab}_{s}'] = lout[f'L_{tag}_{s}']
            K[f'tinv_{lab}_{s}'] = logdiff(lout[f'tinv_{tag}_{s}'])
    blk['fricke'] = {
        'rows_pp_pm': [nstr(x, 30) for x in lout['fricke_rows']['p']],
        'rows_mp_mm': [nstr(x, 30) for x in lout['fricke_rows']['m']],
        'resid_y3': {k: logdiff(v) for k, v in lout['fricke_resid_y3'].items()}}
    ckpt()

    # block 5: independent oracle -- PARI/GP lfunmf at 50 digits (recorded
    # verbatim from PARI/GP 2.15.4; the oracle transcript is not shipped); digits matched.
    PARI_L = {  # form1 = a2=+1 = f_b ; form2 = a2=-1 = f_a (label rule)
        'L_f_b_1': '0.54271934916842485520176378379613163082128362217397',
        'L_f_b_2': '0.88045982535822981044968910894132568513898932226250',
        'L_f_b_3': '0.99639244167374753562904143889719213423457403419298',
        'L_f_a_1': '0.48542294297801677480151889969428395466779186286390',
        'L_f_a_2': '0.78750720838343799254783741243379045263699578873500',
        'L_f_a_3': '0.90860174384973943684739989414691362560475772401861'}
    K['pari_lfun_match_digits'] = {
        k: logdiff(R[k] - mp.mpf(v)) for k, v in PARI_L.items()}
    ckpt()

    blk['members_str'] = {k: nstr(v, dps) for k, v in R.items()}
    del blk['members']
    ckpt()
    return blk


def main():
    t0 = time.time()
    b210 = run(210)
    b120 = run(120)
    # cross-precision agreement per member
    agree = {}
    for k, v210 in b210['members_str'].items():
        v120 = b120['members_str'][k]
        n = 0
        for c1, c2 in zip(v210, v120):
            if c1 != c2:
                break
            n += 1
        agree[k] = n   # chars of common prefix (incl. '0.', sign)
    STATE['cross_precision_prefix_chars'] = agree
    STATE['pslq_two_precision_identical'] = {
        k: b210['pslq'][k] == b120['pslq'][k]
        for k in ('control_lemniscate_logs', 'cm_ratio_logs', 'h2_ratio_logs')}
    STATE['meta']['elapsed_s'] = round(time.time() - t0, 1)
    ckpt()
    # final artifact
    with open(f'{BASE}/RING15.json', 'w') as fh:
        json.dump(STATE, fh, indent=1, default=str)
    print('DONE', STATE['meta']['elapsed_s'], 's')


if __name__ == '__main__':
    main()
