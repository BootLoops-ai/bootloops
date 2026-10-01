"""General-g mixture evidence: brute force vs collapsed contingency-table formula.

Model (bigk.py / check_evidence.py convention: NORMALIZED Dir(1,..,1) probability
measures on all simplices; Z = marginal prob of one fixed ordered sequence):
  Z(U) = int ds prod_a dtheta^(a)  prod_v [ sum_a s_a theta^(a)_v ]^{u_v},
  ds       = (g-1)! * Lebesgue on Delta_{g-1},
  dtheta^a = (k-1)! * Lebesgue on Delta_{k-1}.

CLAIMED COLLAPSE (to be gated):
  Z(U) = [prod_v u_v!] * (g-1)! * ((k-1)!)^g / (N+g-1)!
         * SUM_{m_1+..+m_g=N} CT_U(m) * prod_a [ m_a! / (m_a+k-1)! ],
  CT_U(m) = #{k x g nonneg integer tables, row sums u_v, col sums m_a}
          = [z^m] prod_v h_{u_v}(z_1..z_g).
"""
import sys, time, random
from fractions import Fraction
from math import comb, factorial, lcm, log
from itertools import product as iproduct
from collections import defaultdict

try:
    from gmpy2 import mpz
except ImportError:
    mpz = int

# ---------------- brute force: direct allocation-table enumeration ----------------

def compositions(n, parts):
    """All tuples of `parts` nonneg ints summing to n."""
    if parts == 1:
        yield (n,)
        return
    for first in range(n + 1):
        for rest in compositions(n - first, parts - 1):
            yield (first,) + rest

def Z_brute_g(U, g):
    """Exact evidence by summing over all k x g allocation tables.
    Same measure convention as bigk.py Z_brute, extended to g components:
      mixing integral D_g(m) = (g-1)! prod_a m_a! / (N+g-1)!   (prob. measure on Delta_{g-1})
      component integral E(x_col) = (k-1)! prod_v x_v! / (|x|+k-1)!  (prob. measure on Delta_{k-1})
      per-category multinomial u_v!/prod_a x_va!."""
    k, N = len(U), sum(U)
    def D_g(m):
        num = factorial(g - 1)
        for ma in m:
            num *= factorial(ma)
        return Fraction(num, factorial(N + g - 1))
    def E(col):
        num = factorial(k - 1)
        for x in col:
            num *= factorial(x)
        return Fraction(num, factorial(sum(col) + k - 1))
    tot = Fraction(0)
    percat = [list(compositions(u, g)) for u in U]
    for table in iproduct(*percat):          # table[v] = (x_v1..x_vg)
        c = 1
        for u, row in zip(U, table):         # multinomial u!/prod_a x_va!
            mult = factorial(u)
            for x in row:
                mult //= factorial(x)
            c *= mult
        cols = tuple(sum(row[a] for row in table) for a in range(g))
        w = D_g(cols)
        for a in range(g):
            w *= E(tuple(row[a] for row in table))
        tot += c * w
    return tot

# ---------------- collapsed formula, general g: dict DP over g-1 margins ----------------

def CT_table_dict(U, g):
    """dict {(m_1..m_{g-1}): CT_U(m)} via product of dehomogenized simplex polys
    T_u(z_1..z_{g-1}) = sum_{|x|<=u} z^x  (= h_u(z_1..z_g) with z_g -> 1 tracked implicitly)."""
    P = {(0,) * (g - 1): 1}
    for u in U:
        shifts = [x for s in range(u + 1) for x in compositions(s, g - 1)] if g > 1 else [()]
        newP = defaultdict(int)
        for m, c in P.items():
            for x in shifts:
                key = tuple(mi + xi for mi, xi in zip(m, x))
                newP[key] += c
        P = dict(newP)
    return P

def Z_margin_g(U, g):
    """Collapsed formula via contingency-table DP (general g, dict-based)."""
    k, N = len(U), sum(U)
    P = CT_table_dict(U, g)
    S = Fraction(0)
    for m, ct in P.items():
        mg = N - sum(m)
        if mg < 0:
            continue
        w = Fraction(1)
        for ma in list(m) + [mg]:
            w *= Fraction(factorial(ma), factorial(ma + k - 1))
        S += ct * w
    pref = Fraction(factorial(g - 1) * factorial(k - 1) ** g, factorial(N + g - 1))
    for u in U:
        pref *= factorial(u)
    return pref * S

# ---------------- fast g=3 path: 2D triangular DP with anti-diagonal prefix trick ---------

def CT_rows_g3(U):
    """CT_U(m1,m2) as triangular rows P[a][b], a+b<=N (m3 = N-a-b implicit).
    Multiply by T_u(z1,z2) = sum_{i+j<=u} z1^i z2^j per category, using
      (1-z1)(1-z2) T_u = 1 - s_{u+1} + z1 z2 s_u   (s_d = sum_{i+j=d} z1^i z2^j),
    so P*T_u = double-prefix-sum(P - P*s_{u+1} + z1 z2 (P*s_u)), and the s_d
    convolutions are segment sums along anti-diagonals (O(1) each after
    per-anti-diagonal prefix sums).  O(size) big-int ops per category."""
    Z0, n = mpz(0), 0
    P = [[mpz(1)]]
    for u in U:
        n2 = n + u
        # anti-diagonal prefix sums of P: AD[d][r] = sum_{r'<=r} P[r'][d-r']
        AD = []
        for d in range(n + 1):
            acc, row = Z0, []
            for r in range(d + 1):
                acc += P[r][d - r]
                row.append(acc)
            AD.append(row)
        def seg(d, rlo, rhi):
            # sum_{r=rlo..rhi} P[r][d-r], clipped to valid cells
            if d < 0 or d > n:
                return Z0
            rlo, rhi = max(rlo, 0, d - n), min(rhi, d)  # d-r <= n-r always ok; r<=d ensures b>=0
            if rlo > rhi:
                return Z0
            s = AD[d][rhi]
            if rlo > 0:
                s -= AD[d][rlo - 1]
            return s
        # Q = P * (1 - s_{u+1} + z1 z2 s_u) on triangle a+b <= n2
        Q = []
        for a in range(n2 + 1):
            Pa = P[a] if a <= n else None
            row = []
            for b in range(n2 - a + 1):
                v = Pa[b] if (Pa is not None and b <= n - a) else Z0
                v = v - seg(a + b - u - 1, a - u - 1, a)          # - P*s_{u+1}
                if a > 0 and b > 0:
                    v = v + seg(a + b - u - 2, a - 1 - u, a - 1)  # + z1 z2 (P*s_u)
                row.append(v)
            Q.append(row)
        # double prefix sums: along b in each row, then along a for each column
        for a in range(n2 + 1):
            row, acc = Q[a], Z0
            for b in range(len(row)):
                acc += row[b]
                row[b] = acc
        for a in range(1, n2 + 1):
            prev, row = Q[a - 1], Q[a]
            for b in range(len(row)):
                row[b] += prev[b]
        P, n = Q, n2
    return P

def Z_margin_g3(U):
    """Collapsed formula for g=3 via fast triangular DP + integer-common-denominator sum."""
    g, k, N = 3, len(U), sum(U)
    P = CT_rows_g3(U)
    # weights w(m) = m!/(m+k-1)! = 1/prod_{j=1..k-1}(m+j); common denominator L
    rising = [1] * (N + 1)
    for m in range(N + 1):
        r = 1
        for j in range(1, k):
            r *= m + j
        rising[m] = r
    L = 1
    for r in rising:
        L = lcm(L, r)
    B = [mpz(L // r) for r in rising]
    S_num = mpz(0)
    for a in range(N + 1):
        row, inner = P[a], mpz(0)
        for b in range(N - a + 1):
            ct = row[b]
            if ct:
                inner += ct * B[b] * B[N - a - b]
        S_num += B[a] * inner
    S = Fraction(int(S_num), L ** 3)
    pref = Fraction(factorial(g - 1) * factorial(k - 1) ** g, factorial(N + g - 1))
    for u in U:
        pref *= factorial(u)
    return pref * S

# ---------------- fast g=4 path: z3-slices of 2D triangles ----------------

def CT_slices_g4(U):
    """CT_U(m1,m2,m3) as P[c][a][b] (c = z3 degree, a+b <= n-c), m4 = N-a-b-c implicit.
    Uses T3_u(z1,z2,z3) = sum_{x3=0}^{u} z3^{x3} T2_{u-x3}(z1,z2), so
      (P*T3_u) slice c = prefix2( sum_{x3=0}^{min(u,c)} P_{c-x3} * C2_{u-x3} ),
    with C2_i = 1 - s_{i+1} + z1 z2 s_i handled by anti-diagonal prefix sums per
    source slice.  O(u * size) big-int ops per category."""
    Z0, n = mpz(0), 0
    P = [[[mpz(1)]]]                     # P[c][a][b]
    for u in U:
        n2 = n + u
        # anti-diagonal prefix sums per source slice j: ADS[j][d][r]
        ADS = []
        for j in range(n + 1):
            Pj, nj = P[j], n - j
            AD = []
            for d in range(nj + 1):
                acc, row = Z0, []
                for r in range(d + 1):
                    acc += Pj[r][d - r]
                    row.append(acc)
                AD.append(row)
            ADS.append(AD)
        def seg(j, d, rlo, rhi):
            nj = n - j
            if d < 0 or d > nj:
                return Z0
            rlo, rhi = max(rlo, 0), min(rhi, d)
            if rlo > rhi:
                return Z0
            AD = ADS[j][d]
            s = AD[rhi]
            if rlo > 0:
                s -= AD[rlo - 1]
            return s
        out = []
        for c in range(n2 + 1):
            nc = n2 - c                      # triangle bound for output slice c
            acc = [[Z0] * (nc - a + 1) for a in range(nc + 1)]
            for x3 in range(min(u, c) + 1):
                j, i = c - x3, u - x3        # source slice j, conv with C2_i
                if j > n:
                    continue
                Pj, nj = P[j], n - j
                for a in range(nc + 1):
                    Pa = Pj[a] if a <= nj else None
                    row = acc[a]
                    for b in range(nc - a + 1):
                        v = Pa[b] if (Pa is not None and b <= nj - a) else Z0
                        v = v - seg(j, a + b - i - 1, a - i - 1, a)
                        if a > 0 and b > 0:
                            v = v + seg(j, a + b - i - 2, a - 1 - i, a - 1)
                        row[b] += v
            # double prefix sums on acc
            for a in range(nc + 1):
                row, s = acc[a], Z0
                for b in range(len(row)):
                    s += row[b]
                    row[b] = s
            for a in range(1, nc + 1):
                prev, row = acc[a - 1], acc[a]
                for b in range(len(row)):
                    row[b] += prev[b]
            out.append(acc)
        P, n = out, n2
    return P

def Z_margin_g4(U):
    """Collapsed formula for g=4 via sliced triangular DP + common-denominator sum."""
    g, k, N = 4, len(U), sum(U)
    P = CT_slices_g4(U)
    rising = [1] * (N + 1)
    for m in range(N + 1):
        r = 1
        for j in range(1, k):
            r *= m + j
        rising[m] = r
    L = 1
    for r in rising:
        L = lcm(L, r)
    B = [mpz(L // r) for r in rising]
    S_num = mpz(0)
    for c in range(N + 1):
        Pc, sc = P[c], mpz(0)
        for a in range(N - c + 1):
            row, inner = Pc[a], mpz(0)
            for b in range(N - c - a + 1):
                ct = row[b]
                if ct:
                    inner += ct * B[b] * B[N - a - b - c]
            sc += B[a] * inner
        S_num += B[c] * sc
    S = Fraction(int(S_num), L ** 4)
    pref = Fraction(factorial(g - 1) * factorial(k - 1) ** g, factorial(N + g - 1))
    for u in U:
        pref *= factorial(u)
    return pref * S

# ---------------- helpers ----------------

def random_U(k, umax, rng):
    return tuple(rng.randint(0, umax) for _ in range(k))

def lnfrac(q):
    def lnint(x):
        sh = max(0, x.bit_length() - 500)
        return log(x >> sh) + sh * 0.6931471805599453
    return lnint(q.numerator) - lnint(q.denominator)

def balanced_U(k, N, seed=42):
    rng = random.Random(seed)
    cuts = sorted(rng.sample(range(1, N + k), k - 1))
    U = [b - a - 1 for a, b in zip([0] + cuts, cuts + [N + k])]
    assert sum(U) == N and len(U) == k
    return tuple(U)

# ---------------- gates ----------------

def gate():
    rng = random.Random(20260718)
    print("== GATE g=3: Z_brute_g == Z_margin_g == Z_margin_g3 (exact Fractions) ==")
    n_ok = 0
    while n_ok < 10:
        k = rng.randint(2, 4)
        U = random_U(k, 3, rng)
        if sum(U) == 0:
            continue
        zb = Z_brute_g(U, 3)
        zm = Z_margin_g(U, 3)
        zf = Z_margin_g3(U)
        status = "OK" if (zb == zm == zf) else "MISMATCH"
        print(f"  g=3 k={k} U={U}: brute={zb}  margin={zm}  fast={zf}  {status}")
        if status != "OK":
            print(f"    ratio margin/brute = {zm / zb},  fast/brute = {zf / zb}")
            return False
        n_ok += 1
    print("== GATE g=4: Z_brute_g == Z_margin_g ==")
    n_ok = 0
    while n_ok < 5:
        k = rng.randint(2, 3)
        U = random_U(k, 2, rng)
        if sum(U) == 0:
            continue
        zb = Z_brute_g(U, 4)
        zm = Z_margin_g(U, 4)
        zf = Z_margin_g4(U)
        status = "OK" if zb == zm == zf else "MISMATCH"
        print(f"  g=4 k={k} U={U}: brute={zb}  margin={zm}  fast={zf}  {status}")
        if status != "OK":
            print(f"    ratio margin/brute = {zm / zb}")
            return False
        n_ok += 1
    print("== GATE g=2: Z_margin_g == bigk.Z_fast ==")
    import bigk
    for _ in range(5):
        k = rng.randint(2, 5)
        U = random_U(k, 5, rng)
        za, zm = bigk.Z_fast(U), Z_margin_g(U, 2)
        status = "OK" if za == zm else "MISMATCH"
        print(f"  g=2 k={k} U={U}: Z_fast={za}  margin={zm}  {status}")
        if status != "OK":
            print(f"    ratio margin/Z_fast = {zm / za}")
            return False
    # mid-size cross-check of the fast g=3 DP against the dict DP
    print("== CROSS-CHECK g=3 fast vs dict DP, mid-size ==")
    for U in [(5, 4, 3), (7, 0, 2, 6), (3, 3, 3, 3, 3)]:
        zm, zf = Z_margin_g(U, 3), Z_margin_g3(U)
        status = "OK" if zm == zf else "MISMATCH"
        print(f"  U={U}: dict={zm}  fast={zf}  {status}")
        if status != "OK":
            return False
    print("== CROSS-CHECK g=4 fast vs dict DP, mid-size ==")
    for U in [(4, 3, 2), (2, 0, 3, 2), (5, 5)]:
        zm, zf = Z_margin_g(U, 4), Z_margin_g4(U)
        status = "OK" if zm == zf else "MISMATCH"
        print(f"  U={U}: dict={zm}  fast={zf}  {status}")
        if status != "OK":
            return False
    print("ALL GATES PASS")
    return True

# ---------------- timings ----------------

def timing_g3(k, N):
    U = balanced_U(k, N)
    t0 = time.time()
    Z = Z_margin_g3(U)
    dt = time.time() - t0
    bits = Z.numerator.bit_length() + Z.denominator.bit_length()
    print(f"g=3 k={k:4d} N={N:6d}  lnZ={lnfrac(Z):.6f}  eval={dt:8.2f}s  (exact rational, {bits} bits)")
    return dt

def timing_g4(k, N):
    U = balanced_U(k, N)
    t0 = time.time()
    Z = Z_margin_g4(U)
    dt = time.time() - t0
    bits = Z.numerator.bit_length() + Z.denominator.bit_length()
    print(f"g=4 k={k:4d} N={N:6d}  lnZ={lnfrac(Z):.6f}  eval={dt:8.2f}s  (exact rational, {bits} bits)")
    return dt

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "gate":
        ok = gate()
        sys.exit(0 if ok else 1)
    elif len(sys.argv) > 1 and sys.argv[1] == "time3":
        timing_g3(int(sys.argv[2]), int(sys.argv[3]))
    elif len(sys.argv) > 1 and sys.argv[1] == "time4":
        timing_g4(int(sys.argv[2]), int(sys.argv[3]))
