# Independent arithmetic anchors.
# - Elliptic curve 14.a (14a1: y^2 + xy + y = x^3 + 4x - 6): a_p = wt-2 form 14.2.a.a coeffs.
# - nf_14_4_a_a.json: wt-4 14.4.a.a a_p (LMFDB qexp, 168 primes).
# - nf_34_2_b_a.json: wt-2 34.2.b.a a_p at split primes (17|p)=1 (LMFDB trace form, p<100).
import json, os

HERE = os.path.dirname(os.path.abspath(__file__))

def ap_14a(p):
    """a_p of elliptic curve 14a1 by direct count."""
    n = 0
    for x in range(p):
        rhs = (x * x * x + 4 * x - 6) % p
        # y^2 + (x+1) y - rhs = 0  -> disc = (x+1)^2 + 4 rhs
        d = ((x + 1) * (x + 1) + 4 * rhs) % p
        if d == 0:
            n += 1
        elif pow(d, (p - 1) // 2, p) == 1:
            n += 2
    n += 1  # point at infinity
    return p + 1 - n

def load_nf(name):
    with open(os.path.join(HERE, name)) as f:
        return {int(k): v for k, v in json.load(f).items()}

if __name__ == '__main__':
    nf144 = load_nf('nf_14_4_a_a.json')
    for p in [5, 11, 13, 19, 53, 67, 89, 101]:
        print(p, 'alpha(14a) =', ap_14a(p), ' beta(14.4.a.a) =', nf144.get(p))
