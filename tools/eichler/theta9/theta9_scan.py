# theta9_scan.py — frame-covariant theta-ratio PSLQ scan layer + the Watson
# 15:10:6 positive-control verdict.
#
# Every verdict field in the emitted receipt is computed and written by THIS
# script.
#
# TWO LEGS, by design:
#   REPRO   — byte-faithful re-run of the reference deep_theta protocol (dps 46,
#             raw mp.pslq, maxcoeff 1e5 ratio2 / 1e6 pairprod, imag and
#             denominator filters as in the reference deep_theta protocol) with ONE
#             substitution: the theta constants come from the THETA9 certified
#             engine (FLINT acb_theta balls) instead of the mpmath lattice sum.
#             Output compared canonically against reference deep_theta.json.
#             (Historical note, recorded not repaired: the July protocol's
#             maxcoeff exceeds the pslq_gate capacity law at dps 46; it is
#             reproduced as-is because reproduction IS the control. The
#             HYGIENE leg below is the capacity-lawful version.)
#   HYGIENE — the lockpick.pslq_gate discipline the engine ships with:
#             members_audit at high dps, two-precision canonicalized legs
#             (two_prec_stable), capacity-refused height rungs, planted
#             positive control INSIDE the scan (coeffs >=10x below maxcoeff),
#             sha512-derived negative control (must NULL).
#
# SIGNATURE: from the surviving pairprod relations the script derives the
# radicand triple of the pair-product triangle A=(4,5), B=(6,7), C=(8,9)
# (deep_theta positions; FLINT chars (4,6),(8,9),(12,15)) and checks the
# reference Watson rate signature 15:10:6 = BC:AC:AB.
#
# USAGE: python3 theta9_scan.py --engine-mid out_p320.json --engine-ball out_u320.json \
#          [--reference fixtures/deep_theta.json] --receipt RECEIPT.json [--pilot N]
import argparse, hashlib, json, os, sys, time
from mpmath import mp, mpf, mpc, pslq, sqrt as msqrt

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(_HERE, '..', '..')))
from lockpick import pslq_gate as pg

BNAMES = ['1', 'sqrt2', 'sqrt3', 'sqrt5', 'sqrt6', 'sqrt10', 'sqrt15', 'sqrt30']
BRADS = [1, 2, 3, 5, 6, 10, 15, 30]


def basis_at(dps):
    with mp.workdps(dps):
        return [mpf(1), msqrt(2), msqrt(3), msqrt(5), msqrt(6), msqrt(10),
                msqrt(15), msqrt(30)]


def load_theta(engine_json):
    # parse at high working precision — ambient-dps parsing truncates the
    # engine's ~94-digit midpoints (the lockpick ambient-dps law)
    d = json.load(open(engine_json))
    th2 = {}
    with mp.workdps(120):
        for e in d['theta_sq_even']:
            th2[e['pos']] = mpc(mpf(e['v']['m_re']), mpf(e['v']['m_im']))
    return d, th2


def canon(rel):
    return tuple(pg.canonicalize(list(rel)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--engine-mid', required=True)
    ap.add_argument('--engine-ball', required=True)
    ap.add_argument('--reference',
                    default=os.path.join(_HERE, 'fixtures', 'deep_theta.json'))
    ap.add_argument('--receipt', required=True)
    ap.add_argument('--pilot', type=int, default=0)
    a = ap.parse_args()

    t0 = time.time()
    rec = {'tool': 'theta9', 'module': 'scan',
           'inputs': {'engine_mid': a.engine_mid, 'engine_ball': a.engine_ball,
                      'reference': a.reference}}
    mid_doc, th2m = load_theta(a.engine_mid)
    ball_doc = json.load(open(a.engine_ball))
    reference = json.load(open(a.reference))
    rec['engine_gates'] = {'mid': mid_doc['gates'], 'ball': ball_doc['gates']}

    # certified digit floor per pairprod observable (mode-U balls)
    ball_pp = {tuple(e['quad']): max(float(e['v']['rad_re']), float(e['v']['rad_im']))
               for e in ball_doc.get('pairprod', [])}
    worst_rad = max(ball_pp.values()) if ball_pp else None
    rec['certified_support'] = {
        'mode_U_tau_rad': ball_doc['tau_rad'],
        'pairprod_worst_rad': worst_rad,
        'pairprod_median_rad': sorted(ball_pp.values())[len(ball_pp) // 2] if ball_pp else None,
        'note': 'REPRO leg runs at dps 46 (historic protocol); worst-case certified '
                'support is -log10(worst_rad) digits; the 40-digit HYGIENE leg is '
                'fully inside certified support everywhere.'}

    # ================= REPRO leg (deep_theta.py protocol, dps 46) ============
    mp.dps = 46
    basis = basis_at(46)
    th2 = {p: mpc(v) for p, v in th2m.items()}   # midpoints, truncated at dps 46
    import itertools
    ratio2_hits = []
    n_r2 = 0
    for i in range(10):
        for j in range(10):
            if i == j:
                continue
            if abs(th2[j]) < mpf(10) ** (-15):
                continue
            r = th2[i] / th2[j]
            if abs(r.imag) > mpf(10) ** (-25) * max(1, abs(r.real)):
                continue
            n_r2 += 1
            if a.pilot and n_r2 > a.pilot:
                break
            rel = pslq([r.real] + basis, maxcoeff=10 ** 5, maxsteps=80000)
            if rel and rel[0] != 0:
                ratio2_hits.append({'ij': (i, j), 'rel': list(rel)})
    pairprod_hits = []
    n_pp = 0
    for (ia, ib) in itertools.combinations(range(10), 2):
        for (ic, id_) in itertools.combinations(range(10), 2):
            if {ia, ib} & {ic, id_}:
                continue
            den = th2[ic] * th2[id_]
            if abs(den) < mpf(10) ** (-15):
                continue
            r = (th2[ia] * th2[ib]) / den
            if abs(r.imag) > mpf(10) ** (-25) * max(1, abs(r.real)):
                continue
            n_pp += 1
            if a.pilot and n_pp > a.pilot:
                break
            rel = pslq([r.real] + basis, maxcoeff=10 ** 6, maxsteps=100000)
            if rel and rel[0] != 0:
                pairprod_hits.append({'quad': (ia, ib, ic, id_), 'rel': list(rel),
                                      'val': mp.nstr(r.real, 30)})
    rec['repro'] = {'dps': 46, 'n_ratio2_scanned': n_r2, 'n_pairprod_scanned': n_pp,
                    'ratio2_hits': ratio2_hits, 'pairprod_hits': pairprod_hits,
                    'protocol': 'deep_theta reference params, engine-substituted'}

    # ---- canonical comparison vs reference
    bank_r2 = {tuple(h['ij']): canon(h['rel']) for h in reference['ratio2_hits']}
    bank_pp = {tuple(h['quad']): canon(h['rel']) for h in reference['pairprod_hits']}
    mine_r2 = {tuple(h['ij']): canon(h['rel']) for h in ratio2_hits}
    mine_pp = {tuple(h['quad']): canon(h['rel']) for h in pairprod_hits}
    repro_match = (mine_r2 == bank_r2) and (mine_pp == bank_pp)
    rec['repro_vs_reference'] = {
        'reference_n': {'ratio2': len(bank_r2), 'pairprod': len(bank_pp)},
        'mine_n': {'ratio2': len(mine_r2), 'pairprod': len(mine_pp)},
        'identical_canonical': bool(repro_match),
        'missing': [list(k) for k in bank_pp.keys() - mine_pp.keys()],
        'extra': [list(k) for k in mine_pp.keys() - bank_pp.keys()],
        'rel_mismatch': [list(k) for k in bank_pp if k in mine_pp and bank_pp[k] != mine_pp[k]]}

    # ---- 15:10:6 signature derivation (script-derived, then compared)
    PAIRS = {'A': (4, 5), 'B': (6, 7), 'C': (8, 9)}
    sig = {}
    for x, y in (('A', 'B'), ('A', 'C'), ('B', 'C')):
        q = PAIRS[x] + PAIRS[y]
        h = mine_pp.get(tuple(q))
        if h is None:
            sig[x + y] = None
            continue
        # relation c0*r + sum ci*b_i = 0 with exactly one nonzero ci expected
        support = [(i, c) for i, c in enumerate(h[1:]) if c != 0]
        sig[x + y] = BRADS[support[0][0]] if len(support) == 1 else \
            [(BNAMES[i], c) for i, c in support]
    sig_expected = {'AB': 6, 'AC': 10, 'BC': 15}
    sig_pass = (sig == sig_expected)
    rec['signature'] = {'derived': sig, 'expected_reference': sig_expected,
                        'reference_source': 'reference Watson rate signature '
                                         '(15:10:6 = BC:AC:AB)',
                        'pass': bool(sig_pass)}

    # ================= HYGIENE leg (pslq_gate discipline) ====================
    hy = {'dps_pair': (40, 46), 'maxcoeff': 100, 'maxsteps': 200000, 'resid_slack': 2}
    audit_dps = 160
    with mp.workdps(audit_dps + 40):
        members = {n: mp.nstr(v, audit_dps + 20, strip_zeros=False)
                   for n, v in zip(BNAMES, basis_at(audit_dps + 40))}
    audit = pg.members_audit(members, BNAMES, dps=audit_dps)
    audit_clean = (audit == [])   # [] iff Q-independent to height 1e4 at audit_dps
    hy['members_audit'] = {'clean': bool(audit_clean), 'raw': audit}
    cap = pg.capacity(min(hy['dps_pair']), len(BNAMES), hy['maxcoeff'])
    hy['capacity'] = cap
    cap1000 = pg.capacity(min(hy['dps_pair']), len(BNAMES), 1000)
    hy['capacity_next_rung_1e3'] = {'ok': cap1000['ok'],
                                    'note': 'refused' if not cap1000['ok'] else 'allowed'}

    # scan members at 66 digits (>= max dps_pair + margin)
    with mp.workdps(90):
        mem66 = {n: mp.nstr(v, 66, strip_zeros=False) for n, v in zip(BNAMES, basis_at(90))}

    def hygiene_scan(target_mpf):
        with mp.workdps(90):
            ts = mp.nstr(target_mpf, 60, strip_zeros=False)
        return pg.two_prec_stable(ts, BNAMES, mem66, dps_pair=hy['dps_pair'],
                                  maxcoeff=hy['maxcoeff'], maxsteps=hy['maxsteps'],
                                  resid_slack=hy['resid_slack'])

    # planted positive control (coeffs <=10 = 10x below maxcoeff 100)
    with mp.workdps(90):
        b90 = basis_at(90)
        planted_val = (3 * b90[4] - 2 * b90[6] + 5) / 7     # 7t - 3*sqrt6 + 2*sqrt15 - 5 = 0
    planted_rel = hygiene_scan(planted_val)
    planted_expect = canon([7, -5, 0, 0, 0, -3, 0, 2, 0])
    planted_pass = (planted_rel is not None and canon(planted_rel) == planted_expect)
    # negative control: sha512-derived digits, must NULL. Seed re-derived
    # deliberately at repo adoption (renaming the seed changes the control
    # value by design; the new value must still NULL and is checked here).
    hsh = hashlib.sha512(b'theta9 negative control (repo re-derivation)').hexdigest()
    with mp.workdps(90):
        neg_val = mpf('0.' + str(int(hsh[:64], 16))[:70]) + mpf('0.5')
    neg_rel = hygiene_scan(neg_val)
    neg_pass = (neg_rel is None)
    hy['planted'] = {'target': '(3*sqrt6 - 2*sqrt15 + 5)/7',
                     'found': list(planted_rel) if planted_rel else None,
                     'expected_canonical': list(planted_expect), 'pass': bool(planted_pass)}
    hy['negative'] = {'target': 'sha512-derived, 0.5+frac', 'found':
                      list(neg_rel) if neg_rel else None, 'pass': bool(neg_pass)}

    # hygiene re-scan of the 6 reference-hit observables + null spot-checks
    with mp.workdps(90):
        b46full = None
    hy_hits = {}
    for q in bank_pp:
        r = (th2[q[0]] * th2[q[1]]) / (th2[q[2]] * th2[q[3]])
        rel = hygiene_scan(r.real)
        hy_hits[str(list(q))] = {'rel': list(rel) if rel else None,
                                 'match_reference': bool(rel is not None and
                                                      canon(rel) == bank_pp[q])}
    hy['reference_hit_rescan'] = hy_hits
    hy_hits_pass = all(v['match_reference'] for v in hy_hits.values())
    # null spot-checks: 12 disjoint quads recorded as null must stay null
    null_checks = {}
    n_done = 0
    for (ia, ib) in itertools.combinations(range(10), 2):
        if n_done >= 12:
            break
        for (ic, id_) in itertools.combinations(range(10), 2):
            if {ia, ib} & {ic, id_} or (ia, ib, ic, id_) in bank_pp:
                continue
            r = (th2[ia] * th2[ib]) / (th2[ic] * th2[id_])
            if abs(r.imag) > mpf(10) ** (-25) * max(1, abs(r.real)):
                continue
            rel = hygiene_scan(r.real)
            null_checks[str([ia, ib, ic, id_])] = (list(rel) if rel else None)
            n_done += 1
            break
    hy['null_spot_checks'] = null_checks
    nulls_pass = all(v is None for v in null_checks.values())
    hy['nulls_pass'] = bool(nulls_pass)
    rec['hygiene'] = hy

    # ================= VERDICT (script-emitted) ==============================
    rec['verdict'] = {
        'repro_identical_to_reference': bool(repro_match),
        'signature_15_10_6_pass': bool(sig_pass),
        'hygiene_members_audit_clean': bool(audit_clean),
        'hygiene_planted_pass': bool(planted_pass),
        'hygiene_negative_pass': bool(neg_pass),
        'hygiene_hit_rescan_pass': bool(hy_hits_pass),
        'hygiene_null_spots_pass': bool(nulls_pass),
        'pilot_mode': bool(a.pilot),
        'CONTROL_PASS': bool(repro_match and sig_pass and audit_clean and
                             planted_pass and neg_pass and hy_hits_pass and
                             nulls_pass and not a.pilot),
    }
    rec['wall_s'] = round(time.time() - t0, 2)
    json.dump(rec, open(a.receipt, 'w'), indent=1, default=str)
    print('scan receipt ->', a.receipt)
    for k, v in rec['verdict'].items():
        print(' ', k, '=', v)
    print('wall_s:', rec['wall_s'])
    sys.exit(0 if (rec['verdict']['CONTROL_PASS'] or
                   (a.pilot and repro_match is not None)) else 1)


if __name__ == '__main__':
    main()
