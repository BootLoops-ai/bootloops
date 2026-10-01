#!/usr/bin/env python3
"""bigg.py — W2: the general-g evaluator at scale, from the gated identity.

Notation of record: k = states, g = components, N = sum U, m = allocation
counts, gf dummy z.  Gated identity (GATE_W1_IDENTITY.json all_pass,
normalized Dir(1) convention):

  Z(U) = [prod_v u_v!] (g-1)! [(k-1)!]^g / (N+g-1)!
         * sum_{m composition of N into g parts} CT_U(m) prod_a m_a!/(m_a+k-1)!
  CT_U(m) = # k x g tables, row sums U, col sums m
          = [z^m] prod_v h_{u_v}(z_1..z_g).

ALGORITHM — CT-generating-function DP over states.  Track the g-1 exponents
(e_1..e_{g-1}) of z_1..z_{g-1}; z_g's exponent is determined because the
product is homogeneous of degree S = partial sum of the u_v processed so far
(truncation to total degree N is automatic).  In tracked coordinates
w_a = z_a/z_g, multiplying by h_u(z_1..z_g) is multiplying by
T_u(w) = sum_{|x| <= u} w^x.  All coefficients exact big ints.

Per-state multiply, three paths:
  g = 2 (1D):  (1-w) T_u = 1 - w^{u+1}       -> subtract shift, one prefix sum.
  g = 3 (2D):  (1-w1)(1-w2) T_u = 1 - s_{u+1} + w1 w2 s_u
               (s_d = sum_{i+j=d} w1^i w2^j), the s_d convolutions are
               antidiagonal segment sums (O(1)/cell after per-antidiagonal
               prefix sums), then two prefix passes.  O(1) big-int ops per
               cell per state => total ~ sum_v C(S_v+2,2) ~ k N^2/6 ops,
               INDEPENDENT of per-state u_v granularity.  Same trick as
               biggen.py (vendored beside this file);
               re-implemented here independently and cross-gated against the
               w1 references and the generic path below.
  g >= 4 (generic, any g): Newton recurrence on complete homogeneous polys,
               h_j = sum_{i=1..g} (-1)^{i-1} e_i h_{j-i}: walk the degree up
               one unit at a time, each step a (2^g - 1)-point stencil in
               tracked coordinates (e_i = sum of squarefree monomials; z_g in
               a monomial shifts nothing tracked).  Dict-keyed slices.
               O((2^g-1) * N * C(N+g-1,g-1)) big-int ops — fine for gate/probe
               sizes; a production g=4 path would extend the boundary-
               correction identity (prod_a (1-w_a) T_u = sum_j (-1)^j
               [coeff C(sigma-1,j-1) at |x|=u+j] ...) to O(k * array) as well.

ASSEMBLY — big-int numerator over ONE common denominator (chosen over
per-term fractions.Fraction, justification):
  the weight is prod_a 1/R[m_a] with R[m] = (m+1)(m+2)...(m+k-1) a (k-1)-term
  rising-factorial window (~(k-1) log2 N bits).  Let L = lcm_m R[m]
  (~1.44 (N+k) + O(k) bits, by the psi-function bound plus window prime-power
  excess) and B[m] = L // R[m].  Then
      sum_m CT(m) prod 1/R[m_a] = [sum_m CT(m) prod_a B[m_a]] / L^g,
  a pure big-int accumulation — no per-term gcd (Fraction would gcd ~L-sized
  operands every add), one Fraction at the very end.  For g = 3 the sum runs
  over sorted triples m_1 <= m_2 <= m_3 only (CT_U(m) is symmetric under
  column permutations — relabeling table columns is a bijection), a ~6x cut;
  symmetry is spot-asserted at runtime.

MEMORY / PRACTICAL FRONTIER (31 GB house ulimit): live arrays are ~3
triangles (P, its antidiagonal prefix table, Q) of C(N+g-1,g-1) Python ints
of b ~ log2 prod_v C(u_v+g-1,g-1) bits, i.e. ~3 * C(N+g-1,g-1) *
(44 + 4*ceil(b/30)) bytes.  g=3: N=3000 (~1 GB) easy, N=6000 ~ 5 GB,
N=10^4 ~ 25-30 GB — at the ulimit edge, do not attempt on a <32 GB host.  g=4 via the
generic path: 5 dict slices of C(N+3,3) entries -> N ~ 500-700 is the memory
edge; time binds first.

Checks built into every evaluation: CT checksum sum_m CT_U(m) ==
prod_v C(u_v+g-1,g-1), and (g=3) column-symmetry spot check.
"""

import random
from fractions import Fraction
from math import comb, factorial, lcm

__all__ = ["ct_slice_g2", "ct_slice_g3", "ct_slice_dict", "Z_bigg"]


# ---------------------------------------------------------------- g = 2 path
def ct_slice_g2(U):
    """CT_U(m) as a list P, P[a] = CT at m = (a, N-a).  1D prefix DP."""
    P = [1]
    n = 0
    for u in U:
        if u == 0:
            continue
        n2 = n + u
        Q = [0] * (n2 + 1)
        for b in range(n2 + 1):
            v = P[b] if b <= n else 0
            if b - u - 1 >= 0:
                v -= P[b - u - 1]
            Q[b] = v
        acc = 0
        for b in range(n2 + 1):
            acc += Q[b]
            Q[b] = acc
        P, n = Q, n2
    return P


# ---------------------------------------------------------------- g = 3 path
def ct_slice_g3(U):
    """CT_U(m) as triangular rows P[a][b] (a+b <= N; m = (a, b, N-a-b)).

    Per state u: multiply by T_u via
        P' = prefix_a( prefix_b( P - P*s_{u+1} + w1 w2 (P*s_u) ) ),
    P*s_d read off antidiagonal prefix sums of P:
        (P*s_d)[a,b] = sum_{r=a-d..a} P[r][a+b-d-r]   (clipped).
    """
    P = [[1]]
    n = 0
    for u in U:
        if u == 0:
            continue
        n2 = n + u
        # antidiagonal prefix sums: AD[d][r] = sum_{r'<=r} P[r'][d-r']
        AD = []
        for d in range(n + 1):
            acc = 0
            row = []
            for r in range(d + 1):
                acc += P[r][d - r]
                row.append(acc)
            AD.append(row)

        def seg(d, rlo, rhi):
            if d < 0 or d > n:
                return 0
            if rlo < 0:
                rlo = 0
            if rhi > d:
                rhi = d
            if rlo > rhi:
                return 0
            s = AD[d][rhi]
            if rlo:
                s -= AD[d][rlo - 1]
            return s

        Q = []
        for a in range(n2 + 1):
            Pa = P[a] if a <= n else None
            la = n - a  # max valid b in old row a
            row = []
            for b in range(n2 - a + 1):
                v = Pa[b] if (Pa is not None and b <= la) else 0
                v -= seg(a + b - u - 1, a - u - 1, a)           # - P*s_{u+1}
                if a and b:
                    v += seg(a + b - u - 2, a - u - 1, a - 1)   # + w1w2 P*s_u
                row.append(v)
            Q.append(row)
        # prefix along b, then along a
        for row in Q:
            acc = 0
            for b in range(len(row)):
                acc += row[b]
                row[b] = acc
        for a in range(1, n2 + 1):
            prev, row = Q[a - 1], Q[a]
            for b in range(len(row)):
                row[b] += prev[b]
        P, n = Q, n2
    return P


# ------------------------------------------------------------- generic g path
def ct_slice_dict(U, g):
    """CT_U(m) as dict {(m_1..m_{g-1}): CT} (m_g = N - sum, determined).

    Newton recurrence h_j = sum_{i=1..min(j,g)} (-1)^{i-1} e_i h_{j-i}; per
    state u, walk D_j = (partial product) * h_j for j = 1..u keeping a window
    of g slices.  e_i shifts = squarefree degree-i monomials in z_1..z_g,
    tracked coordinates only (z_g contributes no tracked shift).
    """
    from itertools import combinations

    gm = g - 1
    shifts = []  # shifts[i-1] = tracked shift tuples of e_i
    for i in range(1, g + 1):
        shifts.append([tuple(1 if a in T else 0 for a in range(gm))
                       for T in combinations(range(g), i)])
    D = {(0,) * gm: 1}
    for u in U:
        if u == 0:
            continue
        window = [D] + [None] * (g - 1)  # window[i-1] = D_{j-i}
        for _ in range(u):
            new = {}
            for i in range(1, g + 1):
                src = window[i - 1]
                if not src:
                    continue
                sgn = 1 if (i & 1) else -1
                for sh in shifts[i - 1]:
                    for e, c in src.items():
                        f = tuple(x + y for x, y in zip(e, sh))
                        new[f] = new.get(f, 0) + sgn * c
            window = [new] + window[: g - 1]
        D = window[0]
    return D


# ------------------------------------------------------------------- assembly
def _weight_ints(k, N):
    """R[m] = (m+1)...(m+k-1) = (m+k-1)!/m!; L = lcm; B[m] = L // R[m]."""
    R = [0] * (N + 1)
    r = factorial(k - 1)
    R[0] = r
    for m in range(N):
        r = r * (m + k) // (m + 1)
        R[m + 1] = r
    L = 1
    for r in R:
        L = lcm(L, r)
    B = [L // r for r in R]
    return B, L


def Z_bigg(U, g, measure="dirichlet"):
    """Exact Z(U) for a g-component mixture via the CT-gf DP + assembly."""
    U = list(U)
    k = len(U)
    N = sum(U)
    B, L = _weight_ints(k, N)
    tables = 1
    for u in U:
        tables *= comb(u + g - 1, g - 1)

    if g == 2:
        P = ct_slice_g2(U)
        assert sum(P) == tables, "CT checksum failed (g=2)"
        S_num = sum(P[a] * B[a] * B[N - a] for a in range(N + 1))
    elif g == 3:
        P = ct_slice_g3(U)
        assert sum(map(sum, P)) == tables, "CT checksum failed (g=3)"
        # column-symmetry spot check (CT_U(m) symmetric under permuting m)
        rng = random.Random(1234)
        for _ in range(min(20, (N + 1) ** 2)):
            a = rng.randint(0, N)
            b = rng.randint(0, N - a)
            assert P[a][b] == P[b][a], "CT column-symmetry violated"
        # sum over sorted triples a <= b <= c only, with multiplicity
        S_num = 0
        for a in range(N // 3 + 1):
            row, Ba = P[a], B[a]
            inner = 0
            for b in range(a, (N - a) // 2 + 1):
                ct = row[b]
                if not ct:
                    continue
                c = N - a - b
                if a == b == c:
                    mult = 1
                elif a == b or b == c or a == c:
                    mult = 3
                else:
                    mult = 6
                inner += mult * ct * B[b] * B[c]
            S_num += Ba * inner
    else:
        D = ct_slice_dict(U, g)
        assert sum(D.values()) == tables, "CT checksum failed (dict path)"
        S_num = 0
        for e, ct in D.items():
            mg = N - sum(e)
            w = B[mg]
            for x in e:
                w *= B[x]
            S_num += ct * w

    S = Fraction(S_num, L ** g)
    pref = Fraction(factorial(g - 1) * factorial(k - 1) ** g,
                    factorial(N + g - 1))
    for u in U:
        pref *= factorial(u)
    Z = pref * S
    if measure == "lebesgue":
        Z /= factorial(g - 1) * factorial(k - 1) ** g
    elif measure != "dirichlet":
        raise ValueError(measure)
    return Z


# ------------------------------------------------------------------- selftest
if __name__ == "__main__":
    import os
    import sys
    # sibling engines live beside this file (env-overridable)
    sys.path.insert(0, os.environ.get("MIXALOT_JEFF_DIR")
                    or os.path.dirname(os.path.abspath(__file__)))
    from w1_collapsed import Z_g2, Z_general

    for U, g in [([2, 1], 2), ([3, 0, 2], 2), ([2, 1], 3), ([1, 2, 1], 3),
                 ([2, 2], 4), ([1, 1, 2], 4), ([3], 5)]:
        zb = Z_bigg(U, g)
        zr = Z_general(U, g)
        assert zb == zr, (U, g, zb, zr)
        if g == 2:
            assert zb == Z_g2(U)
    # cross-path: fast g3 slice == dict slice
    for U in [[2, 1, 3], [4, 2], [1, 1, 1, 1]]:
        P = ct_slice_g3(U)
        D = ct_slice_dict(U, 3)
        Nn = sum(U)
        for a in range(Nn + 1):
            for b in range(Nn - a + 1):
                assert P[a][b] == D.get((a, b), 0), (U, a, b)
    # lebesgue convention
    assert Z_bigg([2, 1], 3, "lebesgue") * factorial(2) * factorial(1) ** 3 \
        == Z_bigg([2, 1], 3)
    print("bigg selftest PASS")
