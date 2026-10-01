# theta9_frame.py — marked-frame period point -> tube coords -> tau candidate
# fan — the (Pi, N) entry of the engine (stub entry B).
#
# Fresh implementation of the reference Watson frame law — no reference tau consumed
# here; input is the marked-frame period vector Pi and the exact marked Gram N.
# All emitted verdict fields computed by this script.
#
# THE FRAME LAW (Watson, receipted): with a marked isometry N = U + U(2) + <-4>
# in basis order (e1, f1, e2, f2, g),
#   1. gates: q_{N^-1}(Pi, Pi) ~ 0 (isotropy) and q_{N^-1}(Pi, conj Pi) > 0;
#   2. omega coords v = N^{-1} Pi, normalized at the f1 slot: vn = v / v[slot];
#   3. tube coords z = (vn[z1_slot], vn[z2_slot], vn[z3_slot]);
#   4. paramodular dictionary ambiguity -> candidate fan
#        tau = [[a*w1, c*w3], [c*w3, b*w2]],  (w = s(z), s in {+, -, conj, -conj}),
#        scales (a,b,c) in {1,2,1/2,4,1/4}^3, both (tau11,tau22) orderings,
#      filtered to Im tau > 0 (pos-def), deduped, sorted by dyadic complexity.
# The correct member of the fan SELF-IDENTIFIES downstream (low-height
# algebraic theta ratios); this module only produces the fan.
#
# USAGE:
#   python3 theta9_frame.py --pi-json FILE:key --gram-json FILE:key \
#       [--normalize-slot 1] [--tube-slots 2,3,4] [--dps 120] \
#       [--match-tau "t11re,t11im,t12re,t12im,t22re,t22im" --match-bar 1e-36] \
#       --out OUT.json
# Pi is a list of complex-number strings (python repr '(re+imj)' or 're+imj').
import argparse, json, sys
from mpmath import mp, mpc, mpf, matrix, lu_solve, sqrt as msqrt


def pc(s):
    s = str(s).replace(' ', '').strip('()')
    k = None
    for p in range(1, len(s)):
        if s[p] in '+-' and s[p - 1] not in 'eE':
            k = p
    if s.endswith('j'):
        return mpc(mpf(s[:k]), mpf(s[k:-1])) if k else mpc(0, mpf(s[:-1]))
    return mpc(mpf(s), 0)


def walk(spec):
    path, _, key = spec.partition(':')
    d = json.load(open(path))
    if key:
        for part in key.split('.'):
            d = d[int(part)] if isinstance(d, list) else d[part]
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pi-json', required=True)
    ap.add_argument('--gram-json', required=True)
    ap.add_argument('--normalize-slot', type=int, default=1)
    ap.add_argument('--tube-slots', default='2,3,4')
    ap.add_argument('--dps', type=int, default=120)
    ap.add_argument('--match-tau', default=None)
    ap.add_argument('--match-bar', default='1e-36')
    ap.add_argument('--nstr', type=int, default=50)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    mp.dps = a.dps

    Pi = [pc(s) for s in walk(a.pi_json)]
    n = len(Pi)
    Ng = walk(a.gram_json)
    N = matrix(n, n)
    for i in range(n):
        for j in range(n):
            N[i, j] = mpf(int(Ng[i][j]))

    rec = {'tool': 'theta9', 'module': 'frame', 'n': n,
           'inputs': {'pi_json': a.pi_json, 'gram_json': a.gram_json,
                      'normalize_slot': a.normalize_slot, 'tube_slots': a.tube_slots,
                      'dps': a.dps}}

    # N^{-1}
    Ninv = matrix(n, n)
    for j in range(n):
        e = matrix(n, 1); e[j] = 1
        col = lu_solve(N, e)
        for i in range(n):
            Ninv[i, j] = col[i]

    def quad(M, u, v):
        return sum(u[i] * M[i, j] * v[j] for i in range(n) for j in range(n))

    Pib = [mpc(t).conjugate() for t in Pi]
    scale = sum(abs(t) ** 2 for t in Pi)
    iso = abs(quad(Ninv, Pi, Pi)) / scale
    pos = quad(Ninv, Pi, Pib)
    gate_iso = bool(iso < mpf('1e-40'))
    gate_pos = bool(pos.real > 0 and abs(pos.imag) < mpf('1e-40') * abs(pos.real))
    rec['gates'] = {'isotropy_resid': mp.nstr(iso, 6), 'isotropy_pass': gate_iso,
                    'positivity': mp.nstr(pos, 8), 'positivity_pass': gate_pos}

    v = [sum(Ninv[i, j] * Pi[j] for j in range(n)) for i in range(n)]
    ns = a.normalize_slot
    if abs(v[ns]) < mpf(10) ** (-40):
        rec['warn'] = f'v[{ns}] ~ 0; normalized by v[0] instead'
        vn = [t / v[0] for t in v]
    else:
        vn = [t / v[ns] for t in v]
    ts = [int(x) for x in a.tube_slots.split(',')]
    z = tuple(vn[t] for t in ts)
    rec['tube_z'] = [mp.nstr(t, a.nstr) for t in z]

    # candidate fan
    from itertools import product
    zsets = [z, tuple(-t for t in z), tuple(mpc(t).conjugate() for t in z),
             tuple(-mpc(t).conjugate() for t in z)]
    znames = ['+z', '-z', 'conj', '-conj']
    scales = [mpf(1), mpf(2), mpf(1) / 2, mpf(4), mpf(1) / 4]
    cands = []
    for zi, zz in enumerate(zsets):
        for al, be, ga in product(scales, repeat=3):
            t1, t3, t2 = al * zz[0], be * zz[1], ga * zz[2]
            for oi, (T1, T2, T3) in enumerate([(t1, t2, t3), (t3, t2, t1)]):
                y11, y12, y22 = mpc(T1).imag, mpc(T2).imag, mpc(T3).imag
                if y11 > 0 and (y11 * y22 - y12 ** 2) > 0:
                    cands.append({'tau11': mp.nstr(T1, a.nstr), 'tau12': mp.nstr(T2, a.nstr),
                                  'tau22': mp.nstr(T3, a.nstr),
                                  'zset': znames[zi], 'order': oi,
                                  'scales': [mp.nstr(al, 6), mp.nstr(be, 6), mp.nstr(ga, 6)],
                                  '_t': (T1, T2, T3)})
    seen = set(); uniq = []
    for cd in cands:
        key = (cd['tau11'][:30], cd['tau12'][:30], cd['tau22'][:30])
        if key not in seen:
            seen.add(key); uniq.append(cd)
    import math
    def cx(cd):
        return sum(abs(math.log2(float(mpf(s)))) for s in cd['scales'])
    uniq.sort(key=cx)
    rec['n_candidates'] = len(uniq)

    # optional COPAIR match against a reference tau
    if a.match_tau:
        parts = [mpf(x) for x in a.match_tau.split(',')]
        ref = (mpc(parts[0], parts[1]), mpc(parts[2], parts[3]), mpc(parts[4], parts[5]))
        bar = mpf(a.match_bar)
        best = None
        for ci, cd in enumerate(uniq):
            d = max(abs(cd['_t'][k] - ref[k]) for k in range(3))
            if best is None or d < best[1]:
                best = (ci, d)
        rec['match'] = {'best_cand': best[0], 'maxdiff': mp.nstr(best[1], 6),
                        'bar': a.match_bar, 'pass': bool(best[1] < bar),
                        'best_cand_meta': {k: uniq[best[0]][k] for k in
                                           ('zset', 'order', 'scales')}}
    for cd in uniq:
        del cd['_t']
    rec['tau_candidates'] = uniq
    rec['verdict'] = {'gates_pass': bool(gate_iso and gate_pos),
                      'match_pass': (rec['match']['pass'] if a.match_tau else None)}
    json.dump(rec, open(a.out, 'w'), indent=1, default=str)
    print('frame ->', a.out, '| n_candidates:', len(uniq),
          '| gates_pass:', rec['verdict']['gates_pass'],
          '| match:', rec.get('match', {}).get('pass'), rec.get('match', {}).get('maxdiff', ''))
    ok = rec['verdict']['gates_pass'] and (rec['verdict']['match_pass'] is not False)
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
