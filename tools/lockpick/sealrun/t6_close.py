#!/usr/bin/env python3
# lockpick sealrun member — preregistered sealed-closure runner (CLOSED-with-name vs NULL-with-exclusion-floor against a frozen basket) on the pslq_gate scan engine.
"""t6_close.py — sealed closure runner: one high-precision target vs the frozen 24-member basket.

Usage:
    python3 t6_close.py --target-string-file FILE [--label NAME] [--out DIR]
    python3 t6_close.py --target-json FILE:key.path [--label NAME] [--out DIR]

FILE holds (or key.path addresses) a single high-precision decimal string
(>=150 digits required; >=300 typical).  If the value is complex, pass re and
im separately (one run each).

Protocol (pslq_gate discipline, preregistered BEFORE any target exists — seal
discipline):
  1. sha-pin the basket (BASKET_T6.json from t6_basket.py, gates green inside);
  2. fresh members_audit at scan dps (degenerate-basket NULL is VOID);
  3. fresh planted positive control INSIDE the basket at >=10x headroom, and a
     matched-magnitude negative control (must NULL) — controls fail => exit 2,
     no verdict of any kind is reported;
  4. two-precision canonicalized scans over a HEIGHT LADDER, each rung
     capacity-checked; fit legs use at most (target_digits - 40) digits so the
     final reverify holds out >= 30 genuine digits never seen by any fit leg;
  5. HIT => reverify at (target_digits - 5) dps, pigeonhole-floor line printed,
     canonical vector + both-leg receipt recorded;
  6. NULL through the ladder => honest exclusion floor: "no integer relation
     with |coeff| <= H_max over these 24 members at this precision", with the
     capacity arithmetic printed.
Exit: 0 HIT (closed), 1 NULL (exclusion floor recorded), 2 controls failed,
      3 degenerate pool, 4 input/digits error.
"""
import argparse, hashlib, json, os, sys, time
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))  # package parent
from lockpick import pslq_gate  # noqa: E402

LADDER = [10 ** 4, 10 ** 5, 10 ** 6, 10 ** 8, 10 ** 10, 10 ** 12]


def get_target(args):
    if args.target_string_file:
        s = open(args.target_string_file).read().strip().split()[0]
        return s
    fn, key = args.target_json.split(':', 1)
    d = json.load(open(fn))
    for part in key.split('.'):
        d = d[int(part)] if isinstance(d, list) else d[part]
    return str(d).strip()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--target-string-file')
    ap.add_argument('--target-json')
    ap.add_argument('--label', default='target')
    ap.add_argument('--out', default='.',
                    help='directory for the CLOSE_*.json receipt (default: cwd)')
    a = ap.parse_args()
    if not (a.target_string_file or a.target_json):
        ap.error('need --target-string-file or --target-json')
    t0 = time.time()
    bpath = os.path.join(HERE, 'BASKET_T6.json')
    bsha = hashlib.sha256(open(bpath, 'rb').read()).hexdigest()
    B = json.load(open(bpath))
    assert not B['members_audit']['found'] and B['planted_control']['ok'] \
        and B['negative_control']['ok'], 'basket gates not green'
    M = B['members']
    names = sorted(M)
    n = len(names)
    assert n == 24

    ts = get_target(a)
    if any(c in ts for c in 'jJ(') :
        print('complex target: split re/im and run separately'); sys.exit(4)
    nd = pslq_gate.ndigits(ts)
    min_member_d = min(B['member_digits'].values())
    out = {'producer': 't6_close.py (lockpick sealrun member)',
           'stamp_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
           'label': a.label, 'target_digits': nd, 'basket_sha256': bsha,
           'target_head40': ts[:42], 'rungs': [], 'n_members': n}
    if nd < 150:
        print(f'target carries only {nd} digits; refusing (need >=150)'); sys.exit(4)

    rev_dps = min(nd - 5, min_member_d - 5)
    hi = min(nd - 40, rev_dps - 30)
    lo = hi - 30
    out['dps_plan'] = {'fit_pair': [lo, hi], 'reverify_dps': rev_dps,
                       'held_out_digits': rev_dps - hi}
    assert rev_dps - hi >= 30, 'held-out law violated'

    # fresh audit + controls at scan precision, same pool/protocol
    aud = pslq_gate.members_audit(M, names, dps=min(300, hi), h=1e5)
    out['members_audit'] = aud
    if aud:
        print('DEGENERATE POOL', json.dumps(aud, default=str))
        json.dump(out, open(os.path.join(a.out, f'CLOSE_{stamp_tag(a)}.json'), 'w'),
                  indent=1, default=str)
        sys.exit(3)
    with mp.workdps(rev_dps + 40):
        synth = (3 * mp.mpf(M['pi^2']) - 7 * mp.mpf(M['catalan'])
                 + 2 * mp.mpf(M['log2^3'])) / 5
        sstr = mp.nstr(synth, rev_dps + 30, strip_zeros=False)
        neg = mp.mpf('0.' + hashlib.sha512(ts.encode()).hexdigest()
                     .translate(str.maketrans('abcdef', '123456')) * 8)
        negs = mp.nstr(neg * mp.mpf(M['pi^2']), rev_dps + 30, strip_zeros=False)
    want = {nm: 0 for nm in names}
    want.update({'pi^2': -3, 'catalan': 7, 'log2^3': -2})
    wantvec = pslq_gate.canonicalize([5] + [want[nm] for nm in names])
    pv = pslq_gate.two_prec_stable(sstr, names, M, dps_pair=(lo, hi), maxcoeff=10 ** 6)
    nv = pslq_gate.two_prec_stable(negs, names, M, dps_pair=(lo, hi), maxcoeff=10 ** 6)
    out['controls'] = {'planted_ok': pv == wantvec, 'planted_got': pv,
                       'negative_ok': nv is None, 'negative_got': nv,
                       'headroom_x': 10 ** 6 // 7}
    if not (out['controls']['planted_ok'] and out['controls']['negative_ok']):
        print('CONTROLS FAILED — refusing any verdict', json.dumps(out['controls'],
              default=str))
        json.dump(out, open(os.path.join(a.out, f'CLOSE_{stamp_tag(a)}.json'), 'w'),
                  indent=1, default=str)
        sys.exit(2)

    hit = None
    for h in LADDER:
        cap = pslq_gate.capacity(lo, n, h)
        rung = {'height': h, 'capacity': cap}
        if not cap['ok']:
            rung['skipped'] = 'over capacity at fit legs'
            out['rungs'].append(rung)
            break
        try:
            v = pslq_gate.two_prec_stable(ts, names, M, dps_pair=(lo, hi),
                                          maxcoeff=h, maxsteps=400000)
        except pslq_gate.DegeneratePoolError as e:
            rung['degenerate'] = {'relation': e.relation, 'members': e.members}
            out['rungs'].append(rung)
            print('DEGENERATE at height', h, e)
            json.dump(out, open(os.path.join(a.out, f'CLOSE_{stamp_tag(a)}.json'), 'w'),
                      indent=1, default=str)
            sys.exit(3)
        rung['result'] = list(v) if v else None
        out['rungs'].append(rung)
        if v:
            hit = (h, v)
            break

    if hit:
        h, v = hit
        rv = pslq_gate.reverify({'target': ts, 'members': names, 'vector': list(v)},
                                M, dps=rev_dps, tol=f'1e-{rev_dps - 60}')
        maxc = max(abs(int(c)) for c in v)
        with mp.workdps(50):
            pfloor = -(n * mp.log10(mp.mpf(maxc)))
        out['HIT'] = {'height_rung': h, 'vector': dict(zip(['target'] + names, v)),
                      'reverify': rv, 'max_coeff': maxc,
                      'pigeonhole_floor_log10': float(pfloor),
                      'held_out_digits': rev_dps - hi,
                      'margin_note': (f'relation residual at {rev_dps}d vs pigeonhole '
                                      f'artifact floor 1e{float(pfloor):.0f} — margin is '
                                      'the gap between them')}
        verdict = 'CLOSED' if rv['ok'] else 'HIT-UNVERIFIED (reverify failed — treat as NULL-class, escalate)'
    else:
        hmax = max(r['height'] for r in out['rungs'] if 'result' in r) if any(
            'result' in r for r in out['rungs']) else 0
        out['NULL'] = {'exclusion_height': hmax,
                       'statement': (f'no integer relation with |coeff| <= {hmax:g} '
                                     f'links {a.label} to the 24-member basket at '
                                     f'fit legs ({lo},{hi})d (controls green, audit '
                                     'clean); capacity-checked per rung'),
                       'height_growth_note': ('~x100/weight-step law (lockpick manual): '
                                              'a NULL here does not exclude higher-height '
                                              'relations; escalation needs more digits')}
        verdict = 'NULL-with-floor'
    out['verdict'] = verdict
    out['wall_s'] = round(time.time() - t0, 2)
    fn = os.path.join(a.out, f'CLOSE_{stamp_tag(a)}.json')
    json.dump(out, open(fn, 'w'), indent=1, default=str)
    print('T6 CLOSE', verdict, '->', fn, f"({out['wall_s']}s)")
    sys.exit(0 if hit and verdict == 'CLOSED' else 1)


def stamp_tag(a):
    return (a.label.replace(' ', '_').replace('=', '').replace('/', '_')
            + time.strftime('_%H%M%SZ', time.gmtime()))


if __name__ == '__main__':
    sys.exit(main())
