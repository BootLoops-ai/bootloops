"""frozen_comp_blind.py — BLIND implementation of the frozen-component evidence.

Independent second derivation, written without re-opening v1's source.

DISCIPLINE NOTE: written from DERIVATION_FROZEN_COMP.md alone, AFTER frozen_comp_v1.py
was finished, WITHOUT re-opening or referencing v1's source. Same author
(that limit on independence is acknowledged), but a deliberately different route with
no shared code: v1 expands category-at-a-time with explicit multinomial coefficients
and contracts against a closed-form moment table at the end; THIS file processes the N
observations one at a time (Polya-urn / sequential predictive factorization), never
forms a multinomial coefficient, and accumulates the Pochhammer factors incrementally
inside the recursion. A coefficient bug in one route has no analogue in the other.

Object (definition pin — sequence evidence, no multinomial coefficient):
    Z(u; p, beta) = E_w[ prod_v (sum_a w_a p_{a,v})^{u_v} ],  w ~ Dirichlet(beta).

Route: expand u into its N individual observations v_1..v_N (order irrelevant). Carry
    W_t[m] = sum over component-assignments z of the first t observations with
             count vector m of  prod_{i<=t} p_{z_i, v_i} * prod_a (beta_a)_{m_a},
via  W_{t+1}[m + e_a] += W_t[m] * p_{a, v_{t+1}} * (beta_a + m_a),
then  Z = sum_m W_N[m] / (B)_N,  where (x)_n is the rising factorial.
Exactness: every step is a finite sum of products of rationals; fractions.Fraction
throughout.
"""

from fractions import Fraction


def _rising_from(x, n):
    """x(x+1)...(x+n-1) as an exact Fraction (used only for the final (B)_N)."""
    acc = Fraction(1)
    for i in range(n):
        acc *= x + i
    return acc


def evidence_blind(p, beta, u):
    """Exact Z(u; p, beta) as a Fraction via the observation/urn recursion."""
    g = len(p)
    if g == 0:
        raise ValueError("need at least one component")
    k = len(p[0])
    if any(len(row) != k for row in p):
        raise ValueError("ragged p")
    if len(u) != k:
        raise ValueError("u length != k")
    if len(beta) != g:
        raise ValueError("beta length != g")
    betas = [Fraction(b) for b in beta]
    B = sum(betas)
    # observation list: category v repeated u_v times (order does not matter)
    obs = []
    for v in range(k):
        obs.extend([v] * u[v])
    N = len(obs)

    zero = tuple(0 for _ in range(g))
    layer = {zero: Fraction(1)}
    for v in obs:
        nxt = {}
        for m, wgt in layer.items():
            for a in range(g):
                pav = Fraction(p[a][v])
                if pav == 0:
                    continue
                contrib = wgt * pav * (betas[a] + m[a])
                key = m[:a] + (m[a] + 1,) + m[a + 1:]
                if key in nxt:
                    nxt[key] += contrib
                else:
                    nxt[key] = contrib
        layer = nxt

    return sum(layer.values(), Fraction(0)) / _rising_from(B, N)


def presence_bayes_factor_blind(p_with, beta_with, p_without, beta_without, u):
    """BF = Z(u; S) / Z(u; S minus queried signature), both by the urn route."""
    return evidence_blind(p_with, beta_with, u) / evidence_blind(
        p_without, beta_without, u
    )
