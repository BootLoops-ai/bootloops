# FFP candidate checker: reduce z at each pass-set prime, per-prime IN/OUT + AND verdict.
# Quadratic points: split primes only (both Galois reductions must fire); inert/ramified/
# bad-denominator/singular reductions are SKIPPED (not checkable), matching pilot semantics.
import json, os, sys, glob
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from modsqrt import modsqrt

MIN_CHECKABLE = 5

def load_passsets(op='AESZ34'):
    ps = {}
    for fn in sorted(glob.glob(os.path.join(HERE, 'passsets', f'{op}.p*.json'))):
        d = json.load(open(fn))
        if d.get('status') == 'OK':
            ps[d['p']] = dict(fires=set(d['fires']), sing=set(d['sing']),
                              alphas=d.get('alphas', {}))
    return ps

def residues_mod_p(cand, p):
    """-> (status, [residues]). cand: {'kind':'rational','z':Fraction} or
    {'kind':'quadratic','minpoly':(A,B,C)} with A z^2 + B z + C = 0."""
    if cand['kind'] == 'rational':
        z = cand['z']
        if z.denominator % p == 0:
            return 'skip-den-bad', []
        return 'ok', [z.numerator * pow(z.denominator, p - 2, p) % p]
    A, B, C = cand['minpoly']
    if A % p == 0:
        return 'skip-lead-bad', []
    disc = (B * B - 4 * A * C)
    if disc % p == 0:
        return 'skip-ramified', []
    s = modsqrt(disc % p, p)
    if s is None:
        return 'skip-inert', []
    i2a = pow(2 * A % p, p - 2, p)
    return 'ok', sorted({(-B + s) * i2a % p, (-B - s) * i2a % p})

def check(cand, passsets):
    per, alphas = {}, {}
    for p, ps in sorted(passsets.items()):
        st, res = residues_mod_p(cand, p)
        if st != 'ok':
            per[p] = st; continue
        if any(r == 0 or r in ps['sing'] for r in res):
            per[p] = 'skip-singular'; continue
        inn = all(r in ps['fires'] for r in res)
        per[p] = 'IN' if inn else 'OUT'
        if inn:
            alphas[p] = sorted(set(sum((ps['alphas'].get(str(r), []) for r in res), [])))
    checkable = [p for p, s in per.items() if s in ('IN', 'OUT')]
    n_in = sum(per[p] == 'IN' for p in checkable)
    and_pass = bool(checkable) and n_in == len(checkable)
    return dict(per_prime={str(p): s for p, s in per.items()},
                n_checkable=len(checkable), n_in=n_in, and_pass=and_pass,
                survivor=and_pass and len(checkable) >= MIN_CHECKABLE,
                alphas_at_in={str(p): a for p, a in alphas.items()})

def parse_spec(kind, z=None, minpoly=None):
    """kind 'rational': z='num/den'. kind 'quadratic': minpoly=(A,B,C) ints, or
    (u,v,d) via uvd= for z=u+v*sqrt(d)."""
    if kind == 'rational':
        return dict(kind='rational', z=Fraction(z))
    return dict(kind='quadratic', minpoly=tuple(int(x) for x in minpoly))

def uvd_minpoly(u, v, d):
    u, v = Fraction(u), Fraction(v)
    mp = (Fraction(1), -2 * u, u * u - d * v * v)
    den = 1
    for f in mp:
        den = den * f.denominator // __import__('math').gcd(den, f.denominator)
    return tuple(int(f * den) for f in mp)

if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--op', default='AESZ34')
    ap.add_argument('--rational', help='num/den')
    ap.add_argument('--minpoly', help='A,B,C')
    ap.add_argument('--uvd', help='u,v,d for z=u+v*sqrt(d)')
    a = ap.parse_args()
    ps = load_passsets(a.op)
    if a.rational:
        c = parse_spec('rational', z=a.rational)
    elif a.minpoly:
        c = parse_spec('quadratic', minpoly=a.minpoly.split(','))
    else:
        c = dict(kind='quadratic', minpoly=uvd_minpoly(*a.uvd.split(',')))
    print(json.dumps(dict(op=a.op, spec=str(c), primes=sorted(ps), **check(c, ps))))
