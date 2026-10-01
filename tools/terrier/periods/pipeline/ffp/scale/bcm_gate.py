# BCM_GATE: positive-control + negative-control receipts for the BCM pass-sets.
# Positive: AESZ4 z*=-1/5832, AESZ11 z*=-1/432 must fire at all 12 primes; alpha must
# match the identified wt-2 curve; AESZ11 beta must match 180.4.a.e where archived.
# Negative: pilot-style rational controls must fail the AND.
import json, os
from fractions import Fraction
from check_candidate import load_passsets, check

HERE = os.path.dirname(os.path.abspath(__file__))

def ap_curve(a1, a2, a3, a4, a6, p):
    n = 1
    for x in range(p):
        rhs = (x * x * x + a2 * x * x + a4 * x + a6) % p
        B = (a1 * x + a3) % p
        d = (B * B + 4 * rhs) % p
        n += 1 if d == 0 else (2 if pow(d, (p - 1) // 2, p) == 1 else 0)
    return p + 1 - n

CURVES = {'AESZ4': (1, -1, 0, -3, 3),    # conductor 54, a5=3: ADEK g (54.2, matches qexp)
          'AESZ11': (0, 0, 0, 0, 1)}     # conductor 36, CM sqrt(-3): 36.2.a class
NF180 = {int(k): v for k, v in json.load(open(os.path.join(HERE, 'nf_180_4_a_e.json')))['ap'].items()}
NEG = ['-2', '3', '2/3', '-5/11', '1/7', '-5']

def main():
    rows, ok = [], True
    for op, ctrl in (('AESZ4', Fraction(-1, 5832)), ('AESZ11', Fraction(-1, 432))):
        ps = load_passsets(op)
        r = check(dict(kind='rational', z=ctrl), ps)
        arow = {}
        for p in sorted(ps):
            d = json.load(open(os.path.join(HERE, f'passsets/{op}.p{p}.json')))
            sp = d['control'].get('splits', [])
            am = any(al == ap_curve(*CURVES[op], p) for al, be in sp)
            bm = (any(be == NF180[p] for al, be in sp) if op == 'AESZ11' and p in NF180
                  else None)
            arow[str(p)] = dict(splits=sp, alpha_match_curve=am, beta_match_180_4_a_e=bm)
            ok &= am and (bm is not False)
        g = r['and_pass'] and r['n_checkable'] == len(ps)
        ok &= g
        rows.append(dict(gate=f'{op}-POSITIVE-CONTROL', z=str(ctrl), verdict='PASS' if g else 'FAIL',
                         wt2_curve=CURVES[op], control_receipts=arow,
                         **{k: v for k, v in r.items() if k != 'alphas_at_in'}))
        for z in NEG:
            rn = check(dict(kind='rational', z=Fraction(z)), ps)
            gate = not rn['and_pass']
            ok &= gate
            rows.append(dict(gate=f'{op}-NEG-{z}', verdict='PASS' if gate else 'FAIL',
                             n_checkable=rn['n_checkable'], n_in=rn['n_in']))
    hdr = dict(gate='BCM-GATE', ops=['AESZ4', 'AESZ11'],
               label='CONTROLLED (known points fire; alpha/beta anchored)',
               note='ADEK 2203.09426 prints g in S2new(Gamma0(32)) for the AESZ11 point; '
                    'exact 12-prime alpha identification gives conductor 36 (CM -3) instead; '
                    'beta matches their f=180.4.a.e at all archived primes 61,73,97.',
               verdict='PASS' if ok else 'FAIL')
    with open(os.path.join(HERE, 'BCM_GATE.jsonl'), 'w') as fh:
        fh.write(json.dumps(hdr) + '\n')
        for row in rows:
            fh.write(json.dumps(row) + '\n')
    print('BCM_GATE:', hdr['verdict'])

if __name__ == '__main__':
    main()
