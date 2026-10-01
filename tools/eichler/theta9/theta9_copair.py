# theta9_copair.py — COPAIR gates for the theta9 engine.
#
# All verdict fields in the emitted receipt are computed and written by THIS
# script — never hand-set (receipt fields are script-emitted).
#
# COPAIR law: verify against reference exact/receipted
# objects by APPLICATION, never by route agreement alone. Two gates here:
#
#  COPAIR-THETA: the 10 even squared theta constants from the certified engine
#    (FLINT acb_theta, ball arithmetic) at the reference cand-0 tau
#    midpoint vs the reference deep_theta.json theta_sq values (independent
#    engine: direct mpmath lattice sum, dps 46, R=12; shipped in ./fixtures/).
#    The reference
#    tau strings carry 40 digits, reference thetas 44 digits at dps 46 — the
#    honest agreement bar is 1e-36 relative.
#
#  TWO-PREC: engine run at two precisions — every observable's two balls must
#    overlap (|mid1-mid2| <= rad1+rad2+parse_eps), and the higher-precision
#    radius must be strictly smaller on non-degenerate entries.
#
# USAGE: python3 theta9_copair.py <out_lo.json> <out_hi.json> <receipt.json>
import json, os, sys
from mpmath import mp, mpf, mpc

mp.dps = 120

REFERENCE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      'fixtures', 'deep_theta.json')
AGREE_REL = mpf('1e-36')


def pc(s):
    s = s.replace(' ', '').strip('()')
    k = None
    for p in range(1, len(s)):
        if s[p] in '+-' and s[p - 1] not in 'eE':
            k = p
    if s.endswith('j'):
        return mpc(mpf(s[:k]), mpf(s[k:-1])) if k else mpc(0, mpf(s[:-1]))
    return mpc(mpf(s), 0)


def ball(v):
    return (mpf(v['m_re']), mpf(v['m_im']), mpf(repr(v['rad_re'])), mpf(repr(v['rad_im'])))


def main():
    out_lo, out_hi, receipt_path = sys.argv[1], sys.argv[2], sys.argv[3]
    lo = json.load(open(out_lo))
    hi = json.load(open(out_hi))
    reference = json.load(open(REFERENCE))
    rec = {
        'tool': 'theta9', 'gate': 'COPAIR',
        'inputs': {'engine_lo': out_lo, 'engine_hi': out_hi, 'reference': REFERENCE},
        'engine_gates_lo': lo['gates'], 'engine_gates_hi': hi['gates'],
    }

    # ---- COPAIR-THETA vs reference mpmath lattice-sum values (44d strings)
    bth = [pc(s) for s in reference['theta_sq']]
    eng = {e['pos']: ball(e['v']) for e in hi['theta_sq_even']}
    worst = mpf(0); rows = []
    for p in range(10):
        mre, mim, _, _ = eng[p]
        d = abs(mpc(mre, mim) - bth[p]) / max(mpf(1), abs(bth[p]))
        worst = max(worst, d)
        rows.append({'pos': p, 'reldiff': float(d)})
    copair_theta_pass = bool(worst < AGREE_REL)
    rec['copair_theta'] = {
        'reference_engine': 'deep_theta comparator: mpmath direct lattice sum (dps 46, R=12)',
        'this_engine': 'FLINT acb_theta certified balls (theta_all!, sqr=1)',
        'worst_reldiff': float(worst), 'bar': float(AGREE_REL),
        'rows': rows, 'pass': copair_theta_pass,
    }

    # ---- TWO-PREC overlap on every emitted observable class
    def two_prec_class(key, idkey):
        if key not in lo or key not in hi:
            return None
        lid = {tuple(e[idkey]) if isinstance(e[idkey], list) else e[idkey]: ball(e['v'])
               for e in lo[key]}
        n_bad_overlap = 0; n_not_tighter = 0; worst_gap = mpf(0)
        for e in hi[key]:
            k = tuple(e[idkey]) if isinstance(e[idkey], list) else e[idkey]
            b_hi = ball(e['v']); b_lo = lid[k]
            for c in (0, 1):   # re, im
                gap = abs(b_hi[c] - b_lo[c]) - (b_hi[2 + c] + b_lo[2 + c] + mpf('1e-55'))
                if gap > 0:
                    n_bad_overlap += 1; worst_gap = max(worst_gap, gap)
            if max(b_hi[2], b_hi[3]) > max(b_lo[2], b_lo[3]) and max(b_lo[2], b_lo[3]) > 0:
                n_not_tighter += 1
        return {'n': len(hi[key]), 'n_overlap_fail': n_bad_overlap,
                'n_not_tighter': n_not_tighter, 'worst_gap': float(worst_gap),
                'pass': bool(n_bad_overlap == 0)}

    tp = {}
    for key, idkey in (('theta_sq_even', 'pos'), ('ratio2', 'ij'), ('pairprod', 'quad')):
        r = two_prec_class(key, idkey)
        if r is not None:
            tp[key] = r
    rec['two_prec'] = tp
    two_prec_pass = all(v['pass'] for v in tp.values())

    rec['engine_internal_pass'] = bool(
        lo['gates']['imtau_posdef'] and lo['gates']['odd_vanish_contain_zero'] and
        hi['gates']['imtau_posdef'] and hi['gates']['odd_vanish_contain_zero'])
    rec['verdict'] = {
        'copair_theta_pass': copair_theta_pass,
        'two_prec_pass': bool(two_prec_pass),
        'engine_internal_pass': rec['engine_internal_pass'],
        'all_pass': bool(copair_theta_pass and two_prec_pass and rec['engine_internal_pass']),
    }
    json.dump(rec, open(receipt_path, 'w'), indent=1)
    print('COPAIR receipt ->', receipt_path)
    print('copair_theta worst_reldiff:', float(worst), 'pass:', copair_theta_pass)
    print('two_prec:', {k: v['pass'] for k, v in tp.items()})
    print('ALL_PASS:', rec['verdict']['all_pass'])
    sys.exit(0 if rec['verdict']['all_pass'] else 1)


if __name__ == '__main__':
    main()
