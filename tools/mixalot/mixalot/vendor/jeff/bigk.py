"""Large-k demo: exact evidence of the 2-mixture of one k-state variable, uniform priors.
Identity: since C(u,x)*x!*(u-x)! = u!,
  Z(U) = prod_v u_v! * (k-1)!^2 / (N+1)! * sum_m R_U(m) * m!*(N-m)!/[(m+k-1)!*(N-m+k-1)!]
where R_U(m) = #{0<=x<=U, |x|=m} = [t^m] prod_v (1+...+t^{u_v}).
Gate: against check_evidence-style brute x-sum on random small instances."""
import sys, time, random
from fractions import Fraction
from math import comb, factorial, lgamma, log

def R_poly(U):
    R = [1]
    for u in U:
        newdeg = len(R) + u
        out = [0]*newdeg
        # multiply R by (1 + t + ... + t^u) via prefix sums
        pref = [0]*(len(R)+1)
        for i, c in enumerate(R): pref[i+1] = pref[i] + c
        for j in range(newdeg):
            lo, hi = max(0, j-u), min(len(R)-1, j)
            out[j] = pref[hi+1] - pref[lo]
        R = out
    return R

def Z_fast(U):
    k, N = len(U), sum(U)
    R = R_poly(U)
    s = Fraction(0)
    # rising factorial denominators built incrementally
    for m in range(N+1):
        num = R[m]
        den = 1
        for t in range(1, k): den *= (m+t)
        for t in range(1, k): den *= (N-m+t)
        s += Fraction(num, den)
    pref = Fraction(factorial(k-1)**2, factorial(N+1))
    for u in U: pref *= factorial(u)
    # fold the m!(N-m)!/(m+k-1)!(N-m+k-1)! = 1/(rising products) — already in den
    return pref * s

def Z_brute(U):
    # independent brute x-sum in THEIR convention (check_evidence.py: D_unif, E_unif_block)
    from itertools import product as iproduct
    k, N = len(U), sum(U)
    def D(a,b): return Fraction(factorial(a)*factorial(b), factorial(a+b+1))
    def E(b):
        num = factorial(k-1)
        for x in b: num *= factorial(x)
        return Fraction(num, factorial(sum(b)+k-1))
    tot = Fraction(0)
    for x in iproduct(*[range(u+1) for u in U]):
        c = 1
        for u, xi in zip(U, x): c *= comb(u, xi)
        tot += c * D(sum(x), N-sum(x)) * E(x) * E(tuple(u-xi for u,xi in zip(U,x)))
    return tot

if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'gate':
    random.seed(7)
    for trial in range(12):
        k = random.randint(2, 5)
        U = tuple(random.randint(0, 5) for _ in range(k))
        a, b = Z_fast(U), Z_brute(U)
        assert a == b, (U, a, b)
        print('ok', U, a)
    print('GATE PASS: Z_fast == brute x-sum on 12 random instances')
elif __name__ == '__main__':
    k, N = int(sys.argv[1]), int(sys.argv[2])
    random.seed(42)
    # random composition of N into k parts
    cuts = sorted(random.sample(range(1, N+k), k-1))
    U = [b-a-1 for a,b in zip([0]+cuts, cuts+[N+k])]
    assert sum(U) == N and len(U) == k
    t0 = time.time()
    Z = Z_fast(U)
    t1 = time.time()
    def lnint(x):
        sh = max(0, x.bit_length() - 500)
        return log(x >> sh) + sh * 0.6931471805599453
    lnZv = lnint(Z.numerator) - lnint(Z.denominator)
    print(f"k={k:5d} N={N:8d}  lnZ={lnZv:.6f}  eval={t1-t0:8.3f}s  (exact rational, {Z.numerator.bit_length()+Z.denominator.bit_length()} bits)")
