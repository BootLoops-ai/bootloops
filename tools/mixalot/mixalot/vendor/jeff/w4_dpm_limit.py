#!/usr/bin/env python3
"""w4_dpm_limit.py — W4 THEORY probe: large-g limit of the mixture evidence
under symmetric Dir(alpha/g,..,alpha/g) mixing weights (states stay Dir(1)).

Notation of record: k = states, g = components, t = occupied clusters,
m = allocation counts, N = sum U.

CLAIM CHAIN (each leg checked exactly, fractions.Fraction throughout):

  [A] Brute definition with Dir(alpha/g) weights
        Z_g(U;a) = E[ prod_v (sum_a w_a p_{a,v})^{u_v} ],
        w ~ Dir(a/g,..,a/g), p_a ~ Dir(1,..,1)
      equals the collapsed occupied-cluster regrouping, EXACT AT EVERY g:
        Z_g(U;a) = [prod_v u_v!]/(a)_N * sum_{t=0}^{N} (g)_t_falling / t!
                   * sum_{m=(m_1..m_t), m_i>=1, sum=N} CT_U(m)
                     * prod_i (a/g)_{m_i} (k-1)!/(m_i+k-1)!
      with CT_U(m) = # k x t contingency tables, row sums U, col sums m
      (same CT object as the gated identity, g -> t columns).

  [B] alpha = g reproduces the gated Dir(1,..,1) identity: Z_g(U;g) == Z_general(U,g).

  [C] The g -> infinity limit
        Z_DPM(U;a) = [prod_v u_v!]/(a)_N * sum_{t=1}^{N} (a^t/t!)
                     * sum_{m positive comp of N into t} CT_U(m)
                       * prod_i (m_i-1)! (k-1)!/(m_i+k-1)!
      equals the exact finite-N DPM evidence computed INDEPENDENTLY by
      enumerating set partitions of the N labeled observations with EPPF
      a^t prod_B (|B|-1)!/(a)_N and per-block Dirichlet-multinomial
      (k-1)! prod_v n_{B,v}! / (|B|+k-1)!.

  [D] For every integer g >= 1, Z_g(U;a) is a POLYNOMIAL in 1/g of degree
      <= N-1 with constant term Z_DPM(U;a); the 1/g coefficient is
        C_1 = [prod u_v!]/(a)_N * sum_{t,m} (a^t/t!) CT_U(m)
              * prod_i (m_i-1)!(k-1)!/(m_i+k-1)!
              * [ a*sum_i H_{m_i-1} - t(t-1)/2 ],
      H_n the n-th harmonic number. Verified by exact polynomial
      interpolation of Z_g at N distinct g values + prediction of a held-out g.
"""

from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import factorial

import os
import sys
# sibling engines live beside this file (env-overridable)
sys.path.insert(0, os.environ.get("MIXALOT_JEFF_DIR")
                or os.path.dirname(os.path.abspath(__file__)))
from w1_brute import compositions
from w1_collapsed import ct_table, Z_general


def rising(x, n):
    """(x)_n = x (x+1) ... (x+n-1), Fraction-safe."""
    out = Fraction(1)
    for j in range(n):
        out *= x + j
    return out


def positive_compositions(n, parts):
    """Compositions of n into exactly `parts` positive parts."""
    for c in compositions(n - parts, parts):
        yield tuple(ci + 1 for ci in c)


# ---------------------------------------------------------------- leg A: brute
def Z_brute_alpha(U, g, alpha):
    """Definition, full multinomial expansion; weights ~ Dir(alpha/g,..)."""
    U = list(U)
    k = len(U)
    N = sum(U)
    a = Fraction(alpha)
    total = Fraction(0)
    per_state = [list(compositions(u, g)) for u in U]
    for table in product(*per_state):
        coeff = 1
        for v in range(k):
            c = factorial(U[v])
            for x in table[v]:
                c //= factorial(x)
            coeff *= c
        m = [sum(table[v][comp] for v in range(k)) for comp in range(g)]
        # weight moment under Dir(a/g,..,a/g)
        term = Fraction(coeff)
        for ma in m:
            term *= rising(a / g, ma)
        term /= rising(a, N)
        # state moments under Dir(1,..,1), normalized
        for comp in range(g):
            exps = tuple(table[v][comp] for v in range(k))
            num = factorial(k - 1)
            for e in exps:
                num *= factorial(e)
            term *= Fraction(num, factorial(sum(exps) + k - 1))
        total += term
    return total


# ------------------------------------------- collapsed CT sums (shared kernel)
def ct_positive_sums(U):
    """For each t=1..N and each positive composition m of N into t parts,
    yield (t, m, CT_U(m)). CT via the gated DP with t columns."""
    U = list(U)
    N = sum(U)
    for t in range(1, N + 1):
        dp = ct_table(U, t) if t > 1 else {(): 1}
        if t == 1:
            # single column: CT = 1 (the column is U itself), m = (N,)
            yield 1, (N,), 1
            continue
        for m in positive_compositions(N, t):
            ct = dp.get(tuple(m[:-1]), 0)
            if ct:
                yield t, m, ct


def Z_g_regroup(U, g, alpha):
    """Leg A right-hand side: exact occupied-cluster regrouping, any g >= 1."""
    U = list(U)
    k = len(U)
    N = sum(U)
    a = Fraction(alpha)
    if N == 0:
        return Fraction(1)
    s = Fraction(0)
    for t, m, ct in ct_positive_sums(U):
        fall = Fraction(1)
        for i in range(t):
            fall *= g - i
        if fall == 0:
            continue
        term = fall / factorial(t) * ct
        for mi in m:
            term *= rising(a / Fraction(g), mi) * Fraction(
                factorial(k - 1), factorial(mi + k - 1)
            )
        s += term
    pref = Fraction(1)
    for u in U:
        pref *= factorial(u)
    return pref / rising(a, N) * s


# ------------------------------------------------------------- leg C: the limit
def Z_dpm_ct(U, alpha):
    """The conjectured limit: collapsed occupied-cluster table sum."""
    U = list(U)
    k = len(U)
    N = sum(U)
    a = Fraction(alpha)
    if N == 0:
        return Fraction(1)
    s = Fraction(0)
    for t, m, ct in ct_positive_sums(U):
        term = a**t / factorial(t) * ct
        for mi in m:
            term *= Fraction(
                factorial(mi - 1) * factorial(k - 1), factorial(mi + k - 1)
            )
        s += term
    pref = Fraction(1)
    for u in U:
        pref *= factorial(u)
    return pref / rising(a, N) * s


def set_partitions(items):
    """All set partitions of a list of labeled items."""
    if not items:
        yield []
        return
    first, rest = items[0], items[1:]
    for part in set_partitions(rest):
        for i in range(len(part)):
            yield part[:i] + [part[i] + [first]] + part[i + 1 :]
        yield part + [[first]]


def Z_dpm_crp(U, alpha):
    """Independent leg: exact DPM evidence by set-partition enumeration.
    EPPF a^t prod_B (|B|-1)!/(a)_N, block factor Dirichlet-multinomial."""
    U = list(U)
    k = len(U)
    N = sum(U)
    a = Fraction(alpha)
    if N == 0:
        return Fraction(1)
    obs = []
    for v, u in enumerate(U):
        obs += [v] * u
    total = Fraction(0)
    for part in set_partitions(obs):
        t = len(part)
        term = a**t / rising(a, N)
        for block in part:
            b = len(block)
            term *= factorial(b - 1)
            counts = [0] * k
            for v in block:
                counts[v] += 1
            num = factorial(k - 1)
            for c in counts:
                num *= factorial(c)
            term *= Fraction(num, factorial(b + k - 1))
        total += term
    return total


# ------------------------------------------------------ leg D: 1/g coefficient
def harmonic(n):
    return sum((Fraction(1, j) for j in range(1, n + 1)), Fraction(0))


def C1_formula(U, alpha):
    """Predicted coefficient of 1/g: DPM partition sum with the insertion
    a*sum_i H_{m_i-1} - t(t-1)/2."""
    U = list(U)
    k = len(U)
    N = sum(U)
    a = Fraction(alpha)
    if N == 0:
        return Fraction(0)
    s = Fraction(0)
    for t, m, ct in ct_positive_sums(U):
        term = a**t / factorial(t) * ct
        for mi in m:
            term *= Fraction(
                factorial(mi - 1) * factorial(k - 1), factorial(mi + k - 1)
            )
        insertion = a * sum(harmonic(mi - 1) for mi in m) - Fraction(
            t * (t - 1), 2
        )
        s += term * insertion
    pref = Fraction(1)
    for u in U:
        pref *= factorial(u)
    return pref / rising(a, N) * s


def poly_fit_in_invg(U, alpha, gs):
    """Exact Lagrange interpolation of Z_g_regroup(U,g,alpha) as a polynomial
    in x = 1/g through the nodes gs; returns coefficient list c[0..deg]."""
    xs = [Fraction(1, g) for g in gs]
    ys = [Z_g_regroup(U, g, alpha) for g in gs]
    n = len(xs)
    # Newton -> monomial coefficients, exact
    coef = [Fraction(0)] * n
    # divided differences
    dd = list(ys)
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            dd[i] = (dd[i] - dd[i - 1]) / (xs[i] - xs[i - j])
    # expand Newton form
    poly = [Fraction(0)] * n
    for i in reversed(range(n)):
        # poly = poly*(x - xs[i]) + dd[i]
        newp = [Fraction(0)] * n
        for d in range(n - 1):
            newp[d + 1] += poly[d]
            newp[d] -= poly[d] * xs[i]
        newp[0] += dd[i]
        poly = newp
    return poly


def main():
    alphas = [Fraction(1, 2), Fraction(1), Fraction(3, 2), Fraction(2)]
    cases = [[2, 1], [3, 1], [2, 2], [0, 3, 2], [1, 1, 1, 2]]

    print("== Leg A: brute Dir(alpha/g) == occupied-cluster regrouping (every g) ==")
    ok = True
    for U in cases:
        for g in (1, 2, 3):
            if sum(U) * g > 14 and g > 2:
                continue
            for a in alphas:
                zb = Z_brute_alpha(U, g, a)
                zr = Z_g_regroup(U, g, a)
                match = zb == zr
                ok &= match
                print(f"  U={U} g={g} alpha={a}: {zr}  match={match}")
    print("  LEG A", "PASS" if ok else "FAIL")

    print("== Leg B: alpha=g reproduces gated Dir(1) identity ==")
    okB = True
    for U in cases:
        for g in (2, 3, 4):
            zr = Z_g_regroup(U, g, Fraction(g))
            zg = Z_general(U, g)
            okB &= zr == zg
            print(f"  U={U} g={g}: regroup(alpha=g)={zr} gated={zg} match={zr==zg}")
    print("  LEG B", "PASS" if okB else "FAIL")

    print("== Leg C: collapsed CT limit == CRP set-partition enumeration ==")
    okC = True
    for U in cases + [[4, 3]]:
        for a in alphas:
            zc = Z_dpm_ct(U, a)
            zp = Z_dpm_crp(U, a)
            okC &= zc == zp
            print(f"  U={U} alpha={a}: Z_DPM={zc}  match={zc==zp}")
    # sanity: k=1 must be exactly 1; N=1 must be 1/k
    for a in alphas:
        okC &= Z_dpm_ct([5], a) == 1
        okC &= Z_dpm_ct([1, 0, 0], a) == Fraction(1, 3)
    print("  sanity k=1 -> 1, N=1 -> 1/k:", okC)
    print("  LEG C", "PASS" if okC else "FAIL")

    print("== Leg D: Z_g polynomial in 1/g, deg<=N-1, c0=Z_DPM, c1=C1 formula ==")
    okD = True
    for U in [[2, 1], [2, 2], [0, 3, 2]]:
        N = sum(U)
        for a in [Fraction(1, 2), Fraction(3, 2)]:
            gs = list(range(N, 2 * N))  # N nodes, all >= N
            poly = poly_fit_in_invg(U, a, gs)
            c0, c1 = poly[0], (poly[1] if len(poly) > 1 else Fraction(0))
            zd = Z_dpm_ct(U, a)
            c1f = C1_formula(U, a)
            # held-out prediction at g = 2N+3
            gtest = 2 * N + 3
            pred = sum(c * Fraction(1, gtest) ** i for i, c in enumerate(poly))
            actual = Z_g_regroup(U, gtest, a)
            # small-g check: polynomial must ALSO reproduce g < N (falling factorial zeros)
            small_ok = all(
                sum(c * Fraction(1, gg) ** i for i, c in enumerate(poly))
                == Z_g_regroup(U, gg, a)
                for gg in range(1, N)
            )
            legs = (c0 == zd, c1 == c1f, pred == actual, small_ok)
            okD &= all(legs)
            print(
                f"  U={U} alpha={a}: c0==Z_DPM {legs[0]}, c1==C1 {legs[1]}, "
                f"held-out g={gtest} {legs[2]}, small-g {legs[3]}; c1={c1}"
            )
    print("  LEG D", "PASS" if okD else "FAIL")

    allp = ok and okB and okC and okD
    print("ALL", "PASS" if allp else "FAIL")
    return 0 if allp else 1


if __name__ == "__main__":
    raise SystemExit(main())
