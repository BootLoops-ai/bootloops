# theta9_stub.py — application entry of the theta9 engine — emits
# frame-covariant theta-ratio observables for downstream sealed scans.
# EMISSION ONLY: no PSLQ verdict of any kind is produced here (constant
# verdicts belong to a separate sealed PSLQ layer, e.g. lockpick.pslq_gate).
#
# Verdict-free by design; the emitted gates (posdef, odd-vanish, two-prec) are
# engine consistency gates, script-emitted by the chain below.
#
# ENTRY A (tau direct, g in {2,4}):
#   python3 theta9_stub.py --label L --tau "t11re,t11im,t12re,t12im,t22re,t22im" \
#       [--g 2] [--prec-bits 320] [--tau-rad 1e-48] --outdir DIR
#   (g=4: 20 comma-separated re,im pairs, row-major upper triangle
#    (1,1),(1,2),(1,3),(1,4),(2,2),(2,3),(2,4),(3,3),(3,4),(4,4))
#
# ENTRY B (period point in a marked frame — the (Pi, N) interface):
#   python3 theta9_stub.py --label L --pi-json FILE:key --gram-json FILE:key \
#       [--normalize-slot 1] [--tube-slots 2,3,4] [--n-cands 4] \
#       [--prec-bits 320] [--tau-rad 1e-48] --outdir DIR
#   Runs the Watson frame law ((Pi,N) -> tube coords -> paramodular candidate
#   fan, theta9_frame.py), then the certified engine on the first --n-cands
#   fan members (complexity order). The right member self-identifies
#   downstream by producing low-height algebraic ratios — this stub does not
#   pick it.
#
# OUTPUT per candidate: DIR/OBSERVABLES_<label>_cand<k>.json =
#   engine JSON (theta_sq_even + ratio2 + pairprod certified balls, gates)
#   + 'real_candidates': the observables whose imag ball contains 0 or has
#     |mid_im| <= 1e-25*max(1,|mid_re|) (the deep_theta realness filter) —
#     THESE are the PSLQ targets a downstream sealed scan consumes
#   + provenance block (inputs, frame metadata, stamps).
import argparse, json, os, subprocess, sys, time

CODE = os.path.dirname(os.path.abspath(__file__))
# Eichler.jl project for the certified engine: EICHLER_PROJECT env overrides;
# the default is the repo's upgrades/Eichler.jl relative to this family node.
EICHLER_PROJECT = os.environ.get('EICHLER_PROJECT', os.path.abspath(
    os.path.join(CODE, '..', '..', '..', 'upgrades', 'Eichler.jl')))
JULIA = ['julia', '--project=' + EICHLER_PROJECT,
         os.path.join(CODE, 'theta9_theta.jl')]


def run_engine(label, g, pairs, prec_bits, tau_rad, outdir, emit_pp):
    req = os.path.join(outdir, f'req_{label}.txt')
    out = os.path.join(outdir, f'engine_{label}.json')
    names = [(i, j) for i in range(1, g + 1) for j in range(i, g + 1)]
    lines = [f'label {label}', f'g {g}', f'prec_bits {prec_bits}',
             f'tau_rad {tau_rad}', 'emit_ratio2 1',
             f'emit_pairprod {1 if emit_pp else 0}']
    for (i, j), (re, im) in zip(names, pairs):
        lines.append(f'tau_{i}_{j}_re {re}')
        lines.append(f'tau_{i}_{j}_im {im}')
    open(req, 'w').write('\n'.join(lines) + '\n')
    r = subprocess.run(JULIA + [req, out], capture_output=True, text=True,
                       timeout=1800)
    if r.returncode != 0:
        sys.stderr.write(r.stderr[-2000:])
        raise RuntimeError(f'engine failed for {label}')
    return out


def real_candidates(doc):
    from mpmath import mp, mpf
    mp.dps = 60
    out = []
    for key, idkey in (('ratio2', 'ij'), ('pairprod', 'quad')):
        for e in doc.get(key, []):
            v = e['v']
            mre, mim = mpf(v['m_re']), mpf(v['m_im'])
            rim = mpf(repr(v['rad_im']))
            contains0 = abs(mim) <= rim
            heur = abs(mim) <= mpf('1e-25') * max(1, abs(mre))
            if contains0 or heur:
                out.append({'class': key, 'id': e[idkey],
                            'certified_real': bool(contains0),
                            'm_re': v['m_re'], 'rad': max(float(v['rad_re']),
                                                          float(v['rad_im']))})
    return out


def sp(s):
    s = s.replace(' ', '').strip('()')
    k = None
    for p in range(1, len(s)):
        if s[p] in '+-' and s[p - 1] not in 'eE':
            k = p
    return (s[:k], (s[k:-1] if s[k] == '-' else s[k + 1:-1])) if s.endswith('j') \
        else (s, '0')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--label', required=True)
    ap.add_argument('--tau', default=None)
    ap.add_argument('--g', type=int, default=2)
    ap.add_argument('--pi-json', default=None)
    ap.add_argument('--gram-json', default=None)
    ap.add_argument('--normalize-slot', type=int, default=1)
    ap.add_argument('--tube-slots', default='2,3,4')
    ap.add_argument('--n-cands', type=int, default=4)
    ap.add_argument('--prec-bits', type=int, default=320)
    ap.add_argument('--tau-rad', default='1e-48')
    ap.add_argument('--outdir', required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    stamp = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
    cands = []
    frame_meta = None
    if a.tau:
        parts = a.tau.split(',')
        need = a.g * (a.g + 1)   # re,im per upper-triangle entry
        assert len(parts) == need, f'--tau needs {need} numbers for g={a.g}'
        pairs = [(parts[2 * k], parts[2 * k + 1]) for k in range(need // 2)]
        cands = [('cand0', pairs)]
    else:
        assert a.pi_json and a.gram_json, 'entry B needs --pi-json and --gram-json'
        fr_out = os.path.join(a.outdir, f'frame_{a.label}.json')
        r = subprocess.run([sys.executable, os.path.join(CODE, 'theta9_frame.py'),
                            '--pi-json', a.pi_json, '--gram-json', a.gram_json,
                            '--normalize-slot', str(a.normalize_slot),
                            '--tube-slots', a.tube_slots, '--nstr', '60',
                            '--out', fr_out], capture_output=True, text=True)
        print(r.stdout.strip())
        if r.returncode != 0:
            sys.stderr.write(r.stderr[-2000:])
            raise RuntimeError('frame module failed (gates?)')
        fr = json.load(open(fr_out))
        frame_meta = {'frame_json': fr_out, 'gates': fr['gates'],
                      'n_candidates': fr['n_candidates']}
        for k, cd in enumerate(fr['tau_candidates'][:a.n_cands]):
            pairs = [sp(cd['tau11']), sp(cd['tau12']), sp(cd['tau22'])]
            cands.append((f'cand{k}', pairs, {kk: cd[kk] for kk in
                                              ('zset', 'order', 'scales')}))
    emitted = []
    for tup in cands:
        ck, pairs = tup[0], tup[1]
        eng = run_engine(f'{a.label}_{ck}', a.g, pairs, a.prec_bits, a.tau_rad,
                         a.outdir, emit_pp=(a.g == 2))
        doc = json.load(open(eng))
        doc['provenance'] = {
            'tool': 'eichler.theta9', 'stamp_utc': stamp,
            'entry': 'A' if a.tau else 'B', 'inputs': {
                'tau_arg': bool(a.tau), 'pi_json': a.pi_json,
                'gram_json': a.gram_json, 'tau_rad': a.tau_rad,
                'prec_bits': a.prec_bits},
            'frame': frame_meta, 'cand_meta': (tup[2] if len(tup) > 2 else None),
            'consumer_note': 'PSLQ targets = real_candidates[*].m_re truncated '
                             'to the certified digit count (-log10(rad)); '
                             'sealed verdicts belong to a separate PSLQ '
                             'gate layer, never this stub.'}
        doc['real_candidates'] = real_candidates(doc)
        outp = os.path.join(a.outdir, f'OBSERVABLES_{a.label}_{ck}.json')
        json.dump(doc, open(outp, 'w'), indent=1)
        emitted.append({'cand': ck, 'path': outp,
                        'gates': doc['gates'],
                        'n_real_candidates': len(doc['real_candidates'])})
        print(f'{ck}: -> {outp} | gates {doc["gates"]} | '
              f'real candidates {len(doc["real_candidates"])}')
    idx = {'tool': 'theta9', 'label': a.label, 'stamp_utc': stamp,
           'entry': 'A' if a.tau else 'B', 'emitted': emitted,
           'frame': frame_meta}
    json.dump(idx, open(os.path.join(a.outdir, f'STUB_INDEX_{a.label}.json'),
                        'w'), indent=1)
    print('stub index ->', os.path.join(a.outdir, f'STUB_INDEX_{a.label}.json'))


if __name__ == '__main__':
    main()
