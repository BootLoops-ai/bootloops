# Regression gate: pilot points through the production checker.
# Expect: T1, T2 -> IN at all checkable primes; 9 pilot controls -> AND fails.
import json, os
from fractions import Fraction
from check_candidate import load_passsets, check, uvd_minpoly

HERE = os.path.dirname(os.path.abspath(__file__))

CASES = [
    ('T1(-1/7)', dict(kind='rational', z=Fraction(-1, 7)), 'ALL-IN'),
    ('T2(33+-8s17)', dict(kind='quadratic', minpoly=uvd_minpoly(33, 8, 17)), 'ALL-IN'),
    ('C-rat-1/7', dict(kind='rational', z=Fraction(1, 7)), 'AND-FAIL'),
    ('C-rat--2', dict(kind='rational', z=Fraction(-2)), 'AND-FAIL'),
    ('C-rat-3', dict(kind='rational', z=Fraction(3)), 'AND-FAIL'),
    ('C-rat--5', dict(kind='rational', z=Fraction(-5)), 'AND-FAIL'),
    ('C-rat-2/3', dict(kind='rational', z=Fraction(2, 3)), 'AND-FAIL'),
    ('C-rat--5/11', dict(kind='rational', z=Fraction(-5, 11)), 'AND-FAIL'),
    ('C-alg-30+8s17', dict(kind='quadratic', minpoly=uvd_minpoly(30, 8, 17)), 'AND-FAIL'),
    ('C-alg-5+2s17', dict(kind='quadratic', minpoly=uvd_minpoly(5, 2, 17)), 'AND-FAIL'),
    ('C-alg-33+8s13', dict(kind='quadratic', minpoly=uvd_minpoly(33, 8, 13)), 'AND-FAIL'),
]

def main():
    ps = load_passsets('AESZ34')
    rows, ok = [], True
    for name, cand, expect in CASES:
        r = check(cand, ps)
        if expect == 'ALL-IN':
            gate = r['and_pass'] and r['n_checkable'] >= 5
        else:
            gate = not r['and_pass']
        ok &= gate
        rows.append(dict(name=name, expected=expect, gate='PASS' if gate else 'FAIL',
                         spec=dict(kind=cand['kind'],
                                   z=str(cand.get('z', '')),
                                   minpoly=list(cand.get('minpoly', []))), **r))
        print(name, expect, '->', 'PASS' if gate else 'FAIL',
              f"({r['n_in']}/{r['n_checkable']} IN)")
    with open(os.path.join(HERE, 'CHECKER_GATE.jsonl'), 'w') as fh:
        fh.write(json.dumps(dict(gate='CHECKER-REGRESSION', primes=sorted(ps),
                                 n_primes=len(ps), verdict='PASS' if ok else 'FAIL')) + '\n')
        for row in rows:
            fh.write(json.dumps(row) + '\n')
    print('CHECKER_GATE verdict:', 'PASS' if ok else 'FAIL')
    return ok

if __name__ == '__main__':
    main()
