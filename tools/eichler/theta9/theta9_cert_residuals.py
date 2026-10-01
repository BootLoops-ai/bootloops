# theta9_cert_residuals.py — certified enclosures of the Watson 15:10:6
# relations from the mode-U engine balls.
#
# All receipt fields are script-emitted.
#
# For each reference hit (c0, member) the certified residual enclosure is
#   |c0*r + c_m*sqrt(m)| <= |c0*mid + c_m*sqrt(m)| + |c0|*rad
# (sqrt(m) at 90 dps, exact to 1e-88 — negligible), r the mode-U ball
# (midpoint mid, radius rad) of the pair-product observable. This upgrades
# the PSLQ hit to a CERTIFIED near-identity at the frame's tau ball: the
# relation holds to within the printed enclosure for EVERY tau in the input
# ball (tau_rad = transport uncertainty). Not a proof of exact equality —
# the exactness claim stays with the reference Watson receipts; this is the
# certified-transport version of it.
import json, sys
from mpmath import mp, mpf, sqrt as msqrt

mp.dps = 90
BRADS = {4: 6, 5: 10, 6: 15}   # rel[1:] position -> radicand (BNAMES order)


def main():
    eng, reference_p, out = sys.argv[1], sys.argv[2], sys.argv[3]
    d = json.load(open(eng))
    reference = json.load(open(reference_p))
    pp = {tuple(e['quad']): e['v'] for e in d['pairprod']}
    rows = []
    worst = mpf(0)
    for h in reference['pairprod_hits']:
        q = tuple(h['quad'])
        rel = h['rel']
        c0 = rel[0]
        support = [(i, c) for i, c in enumerate(rel[1:]) if c != 0]
        (bi, cm), = support
        v = pp[q]
        mid_re, mid_im = mpf(v['m_re']), mpf(v['m_im'])
        rad = max(mpf(repr(v['rad_re'])), mpf(repr(v['rad_im'])))
        resid = abs(c0 * mid_re + cm * msqrt(BRADS[bi])) + abs(c0) * rad
        imag_bound = abs(mid_im) + rad
        worst = max(worst, resid, imag_bound)
        rows.append({'quad': list(q), 'relation': rel,
                     'member': f'sqrt{BRADS[bi]}',
                     'resid_encl': mp.nstr(resid, 4),
                     'imag_encl': mp.nstr(imag_bound, 4)})
    rec = {'tool': 'theta9', 'module': 'cert_residuals',
           'engine_ball': eng, 'tau_rad': d['tau_rad'],
           'rows': rows, 'worst_enclosure': mp.nstr(worst, 4),
           'statement': 'each reference 15:10:6 relation holds to within '
                        'resid_encl, and each observable is real to within '
                        'imag_encl, CERTIFIED for every tau in the input ball',
           'pass_1e-40': bool(worst < mpf('1e-40'))}
    json.dump(rec, open(out, 'w'), indent=1)
    print('cert_residuals ->', out, '| worst enclosure:', mp.nstr(worst, 4),
          '| pass_1e-40:', rec['pass_1e-40'])
    sys.exit(0 if rec['pass_1e-40'] else 1)


if __name__ == '__main__':
    main()
