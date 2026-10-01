"""frozen_comp_v1.py — exact evidence for frozen components + Dirichlet weights.

Implementation 1 (category route) of the independently double-derived
frozen-component evidence (the derivation note is not distributed here).

Z(u; p, beta) = E_w[ prod_v (sum_a w_a p_{a,v})^{u_v} ],  w ~ Dirichlet(beta),
computed EXACTLY: dynamic programming over categories expands the integrand into
monomials w^m (m a weak composition of N into g parts) with rational coefficients
c_m, then Z = sum_m c_m * prod_a (beta_a)_{m_a} / (B)_N.

Definition pin: sequence evidence, NO multinomial coefficient N!/prod u_v!.
All arithmetic in fractions.Fraction; p rows and betas must be rational.

Route: category-at-a-time polynomial multiplication with explicit multinomial
coefficients per category, then a single closed-form moment contraction at the end.
(The blind implementation uses the observation-at-a-time Polya-urn route instead.)
"""

from fractions import Fraction
from math import comb


def rising(x, n):
    """Rising factorial (Pochhammer) x(x+1)...(x+n-1), exact for Fraction x."""
    out = Fraction(1)
    for j in range(n):
        out *= x + j
    return out


def _weak_compositions(total, parts):
    """Yield all weak compositions of `total` into `parts` nonnegative ints."""
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for rest in _weak_compositions(total - first, parts - 1):
            yield (first,) + rest


def _category_poly(pcol, uv, g):
    """Multinomial expansion of (sum_a w_a * pcol[a])**uv as {n: coeff}.

    n runs over weak compositions of uv into g parts;
    coeff = multinomial(uv; n) * prod_a pcol[a]**n[a], exact Fractions.
    """
    poly = {}
    for n in _weak_compositions(uv, g):
        coeff = Fraction(1)
        rem = uv
        for na in n:
            coeff *= comb(rem, na)
            rem -= na
        for a in range(g):
            if n[a]:
                coeff *= Fraction(pcol[a]) ** n[a]
        if coeff:
            poly[n] = coeff
    return poly


def expand_coefficients(p, u):
    """DP over categories: coefficient table {m: c_m} of prod_v (w . p_v)^{u_v}.

    p: list of g rows, each a list of k rationals (component distributions).
    u: list of k nonnegative ints (category counts).
    Returns dict mapping weak compositions m of N=sum(u) into g parts -> Fraction c_m.
    """
    g = len(p)
    k = len(p[0])
    assert all(len(row) == k for row in p), "ragged p"
    assert len(u) == k, "u/k mismatch"
    table = {(0,) * g: Fraction(1)}
    for v in range(k):
        if u[v] == 0:
            continue
        pcol = [p[a][v] for a in range(g)]
        cpoly = _category_poly(pcol, u[v], g)
        new = {}
        for m, cm in table.items():
            for n, cn in cpoly.items():
                key = tuple(m[a] + n[a] for a in range(g))
                w = cm * cn
                if key in new:
                    new[key] += w
                else:
                    new[key] = w
        table = new
    return table


def dirichlet_moment(m, beta):
    """E[prod_a w_a^{m_a}] for w ~ Dirichlet(beta), general rational beta. Exact."""
    beta = [Fraction(b) for b in beta]
    B = sum(beta)
    N = sum(m)
    num = Fraction(1)
    for a, ma in enumerate(m):
        num *= rising(beta[a], ma)
    return num / rising(B, N)


def evidence(p, beta, u):
    """Exact Z(u; p, beta) as a Fraction. See module docstring."""
    assert len(beta) == len(p), "beta/g mismatch"
    beta = [Fraction(b) for b in beta]
    B = sum(beta)
    N = sum(u)
    table = expand_coefficients(p, u)
    denom = rising(B, N)
    total = Fraction(0)
    for m, cm in table.items():
        num = Fraction(1)
        for a, ma in enumerate(m):
            num *= rising(beta[a], ma)
        total += cm * num
    return total / denom


def presence_bayes_factor(p_with, beta_with, p_without, beta_without, u):
    """BF = Z(u; p_with, beta_with) / Z(u; p_without, beta_without), exact Fraction.

    p_with: component set S including the queried signature; p_without: S minus it.
    The sequence-vs-counts multinomial coefficient cancels in this ratio.
    """
    return evidence(p_with, beta_with, u) / evidence(p_without, beta_without, u)
