import sys, time, random
from fractions import Fraction
from math import factorial, log
from bigk import R_poly, Z_fast

def Z_fast2(U):
    k, N = len(U), sum(U)
    R = R_poly(U)
    # r[m] = prod_{t=1}^{k-1} (m+t), incremental: r[m+1] = r[m]*(m+k)//(m+1)
    r = [0]*(N+1)
    r0 = 1
    for t in range(1, k): r0 *= t
    r[0] = r0
    for m in range(N):
        r[m] = r[m] if m == 0 else r[m]
        r[m+1] = r[m]*(m+k)//(m+1)
    terms = [Fraction(R[m], r[m]*r[N-m]) for m in range(N+1)]
    while len(terms) > 1:
        nxt = [terms[i]+terms[i+1] for i in range(0, len(terms)-1, 2)]
        if len(terms) % 2: nxt.append(terms[-1])
        terms = nxt
    pref = Fraction(factorial(k-1)**2, factorial(N+1))
    for u in U: pref *= factorial(u)
    return pref * terms[0]

if __name__ == '__main__':
    k, N = int(sys.argv[1]), int(sys.argv[2])
    random.seed(42)
    cuts = sorted(random.sample(range(1, N+k), k-1))
    U = [b-a-1 for a,b in zip([0]+cuts, cuts+[N+k])]
    if len(sys.argv) > 3 and sys.argv[3] == 'check':
        assert Z_fast2(U) == Z_fast(U); print('agree with v1')
    t0 = time.time(); Z = Z_fast2(U); t1 = time.time()
    def lnint(x):
        sh = max(0, x.bit_length()-500)
        return log(x >> sh) + sh*0.6931471805599453
    print(f"k={k:5d} N={N:8d}  lnZ={lnint(Z.numerator)-lnint(Z.denominator):.6f}  eval={t1-t0:8.3f}s")
