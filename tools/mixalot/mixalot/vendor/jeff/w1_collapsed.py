#!/usr/bin/env python3
"""w1_collapsed.py — the two TARGET IDENTITIES, implemented exactly as stated.

Notation of record: k = states, g = components, N = sum U, dummy variable z.
Prefactors are for NORMALIZED Dir(1,..,1) measures on every simplex
(measure='dirichlet'); measure='lebesgue' is the bare-Lebesgue variant, which
is smaller by exactly (g-1)! * [(k-1)!]^g.

(g = 2):
  Z(U) = [prod_v u_v!] * (k-1)!^2 / (N+1)!
         * sum_{m=0}^N R_U(m) * m!(N-m)! / [(m+k-1)!(N-m+k-1)!],
  R_U(m) = #{x : 0 <= x_v <= u_v, sum_v x_v = m}
         = [z^m] prod_v (1 + z + ... + z^{u_v}).

(general g):
  Z(U) = [prod_v u_v!] * (g-1)! [(k-1)!]^g / (N+g-1)!
         * sum_{m = (m_1..m_g) |= N} CT_U(m) * prod_a m_a!/(m_a+k-1)!,
  CT_U(m) = # of k x g nonnegative-integer contingency tables with row sums U
  and column sums m; computed by dynamic programming over the k states,
  tracking the g-1 running column sums (m_g is determined by N).
"""

from collections import defaultdict
from fractions import Fraction
from math import factorial

from w1_brute import compositions


def R_poly(U):
    """Coefficient list of prod_v (1 + z + ... + z^{u_v}); R_poly(U)[m] = R_U(m)."""
    coeffs = [1]
    for u in U:
        new = [0] * (len(coeffs) + u)
        for i, c in enumerate(coeffs):
            for j in range(u + 1):
                new[i + j] += c
        coeffs = new
    return coeffs


def Z_g2(U, measure="dirichlet"):
    """The g=2 identity, verbatim."""
    U = list(U)
    k = len(U)
    N = sum(U)
    R = R_poly(U)
    s = Fraction(0)
    for m in range(N + 1):
        s += Fraction(
            R[m] * factorial(m) * factorial(N - m),
            factorial(m + k - 1) * factorial(N - m + k - 1),
        )
    pref = Fraction(1)
    for u in U:
        pref *= factorial(u)
    Z = pref * Fraction(factorial(k - 1) ** 2, factorial(N + 1)) * s
    if measure == "lebesgue":
        Z /= factorial(2 - 1) * factorial(k - 1) ** 2
    elif measure != "dirichlet":
        raise ValueError(measure)
    return Z


def ct_table(U, g):
    """dict: (m_1..m_{g-1}) -> CT_U(m), tables with row sums U, col sums m
    (m_g = N - sum of the tracked g-1). DP over states."""
    dp = {(0,) * (g - 1): 1}
    for u in U:
        ndp = defaultdict(int)
        for state, cnt in dp.items():
            for split in compositions(u, g):
                key = tuple(state[a] + split[a] for a in range(g - 1))
                ndp[key] += cnt
        dp = dict(ndp)
    return dp


def Z_general(U, g, measure="dirichlet"):
    """The general-g identity, verbatim (CT by DP)."""
    U = list(U)
    k = len(U)
    N = sum(U)
    dp = ct_table(U, g)
    s = Fraction(0)
    for state, ct in dp.items():
        m = list(state) + [N - sum(state)]
        term = Fraction(ct)
        for ma in m:
            term *= Fraction(factorial(ma), factorial(ma + k - 1))
        s += term
    pref = Fraction(1)
    for u in U:
        pref *= factorial(u)
    Z = pref * Fraction(factorial(g - 1) * factorial(k - 1) ** g, factorial(N + g - 1)) * s
    if measure == "lebesgue":
        Z /= factorial(g - 1) * factorial(k - 1) ** g
    elif measure != "dirichlet":
        raise ValueError(measure)
    return Z


if __name__ == "__main__":
    print("Z_g2([2,1]) =", Z_g2([2, 1]))
    print("Z_general([2,1], g=2) =", Z_general([2, 1], 2))
    assert Z_g2([2, 1]) == Z_general([2, 1], 2)
    print("g=2 forms agree on [2,1]  PASS")
