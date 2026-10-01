# BCM validation battery:
# V1 Legendre F_{p^2}: H_{p^2}(lam) == a_p(lam)^2 - 2p for all lam (f=2 machinery).
# V2 d=4 unit-root congruence: trace4(z) == sum_{n<p} A_n z^n mod p, ALL z, per op/prime;
#    variant scan (sign, twist rule in {1, chi(-1)}, dict t=Mz vs 1/(Mz)) -> unique+uniform.
# V3 joint F_p/F_p2: b=(a^2-s2)/(2p) integral + weil_circle at ALL smooth z.
import json, numpy as np
from fractions import Fraction
from bcm import GaussField, h_coeffs, h_all_fp
from bcm_calibrate import brute_ap_legendre
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from weil import weil_circle

OPS = {'h33': (['1/3', '1/3', '2/3', '2/3'], 729),
       'h34': (['1/3', '2/3', '1/4', '3/4'], 1728),
       'h5': (['1/5', '2/5', '3/5', '4/5'], 3125)}

def series_trunc_mod_p(alpha, M, p):
    A = np.zeros(p, dtype=np.int64); A[0] = 1
    fracs = [Fraction(a) for a in alpha]
    cur = 1
    for n in range(1, p):
        c = M % p
        for a in fracs:
            fr = a + (n - 1)
            c = c * (fr.numerator % p) % p * pow(fr.denominator, p - 2, p) % p
        c = c * pow(pow(n, 4, p), p - 2, p) % p
        cur = cur * c % p
        A[n] = cur
    return A

def horner_all(A, p):
    z = np.arange(p, dtype=np.int64)
    S = np.full(p, A[p - 1], dtype=np.int64)
    for n in range(p - 2, -1, -1):
        S = (S * z + A[n]) % p
    return S

def trace_all(op, p, deg):
    alpha, M = OPS[op]
    gf = GaussField(p, deg)
    H = h_all_fp(gf, h_coeffs(gf, alpha))
    resid = float(np.max(np.abs(H[1:] - np.round(H.real[1:]))))
    return np.round(H.real).astype(np.int64), resid

def v1():
    ok = True
    for p in (13, 19):
        ap = brute_ap_legendre(p)
        H2 = h_all_fp(GaussField(p, 2), h_coeffs(GaussField(p, 2), ['1/2', '1/2']))
        r = float(np.max(np.abs(H2[2:] - np.round(H2.real[2:]))))
        good = all(int(round(H2[l].real)) == ap[l] ** 2 - 2 * p for l in ap if l != 1)
        print(f'V1 p={p}: match={good} resid={r:.2e}'); ok &= good
    return ok

def v2(op, primes):
    alpha, M = OPS[op]
    sets = []
    for p in primes:
        H, resid = trace_all(op, p, 1)
        S = horner_all(series_trunc_mod_p(alpha, M, p), p)
        chim1 = 1 if p % 4 == 1 else -1
        cur = set()
        for s in (1, -1):
            for rule in ('1', 'phi'):
                tw = 1 if rule == '1' else chim1
                for dic in ('Mz', '1/(Mz)'):
                    good = True
                    for z in range(1, p):
                        t = M * z % p
                        if dic != 'Mz':
                            t = pow(t, p - 2, p) if t else 0
                        if t in (0, 1):
                            continue
                        if (s * tw * H[t] - S[z]) % p:
                            good = False; break
                    if good:
                        cur.add((s, rule, dic))
        print(f'V2 {op} p={p}: resid={resid:.2e} variants={sorted(cur)}')
        sets.append(cur)
    inter = set.intersection(*sets)
    print(f'V2 {op} intersection: {sorted(inter)}')
    return inter

def v3(op, p, s_sign=1):
    """Full joint test at p=1 mod 12: a=-trace_p, s2=trace_p2, all smooth z."""
    alpha, M = OPS[op]
    H1, r1 = trace_all(op, p, 1)
    H2, r2 = trace_all(op, p, 2)
    Minv = pow(M % p, p - 2, p)
    nb = ni = nw = 0
    for z in range(1, p):
        t = M * z % p
        if t in (0, 1):
            continue
        a, s2 = -s_sign * int(H1[t]), s_sign * int(H2[t])
        nb += 1
        if (a * a - s2) % (2 * p):
            continue
        ni += 1
        if weil_circle(a, (a * a - s2) // (2 * p), p):
            nw += 1
    print(f'V3 {op} p={p}: resid=({r1:.1e},{r2:.1e}) smooth={nb} b-integral={ni} weil={nw}')
    return nb == ni == nw

if __name__ == '__main__':
    ok = v1()
    i33 = v2('h33', [61, 67, 73, 79])
    i34 = v2('h34', [61, 73, 97])
    i5 = v2('h5', [11, 31, 41, 61, 71])
    ok3 = v3('h33', 61) & v3('h34', 61) & v3('h33', 73) & v3('h34', 73)
    print('V1', 'PASS' if ok else 'FAIL', '| V3', 'PASS' if ok3 else 'FAIL')
