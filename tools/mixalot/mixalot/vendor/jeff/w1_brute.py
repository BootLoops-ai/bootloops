#!/usr/bin/env python3
"""w1_brute.py — Z from the DEFINITION, by full expansion. No shortcuts.

MODEL (notation of record: k states, g components, N = sum U, dummy z):
  one k-state categorical variable; data = counts U = (u_1..u_k);
  component a has p_a ~ Dir(1,..,1) on k states (independent across a);
  mixing weights w ~ Dir(1,..,1) on g components;
  Z(U) = E[ prod_v (sum_a w_a p_{a,v})^{u_v} ].

METHOD: for each state v expand (sum_a w_a p_{a,v})^{u_v} multinomially over
the split x_{v,*} = (x_{v,1},..,x_{v,g}) of the u_v observations of state v
among the g components:
  (sum_a w_a p_{a,v})^{u_v}
      = sum_{x_{v,*}} u_v!/prod_a x_{v,a}! * prod_a w_a^{x_{v,a}} p_{a,v}^{x_{v,a}}.
Each monomial is integrated EXACTLY with the Dirichlet-multinomial formula:
  normalized Dir(1,..,1) on the (d-1)-simplex:
      E[prod_i y_i^{n_i}] = (d-1)! * prod_i n_i! / (sum_i n_i + d - 1)!
  bare Lebesgue on the same simplex (measure='lebesgue'):
      int prod_i y_i^{n_i} dy = prod_i n_i! / (sum_i n_i + d - 1)!
(the two differ by the Dir(1) density constant (d-1)!). All arithmetic in
fractions.Fraction. Cost = prod_v C(u_v+g-1, g-1) terms: fine for g<=4,
sum U <= ~12 (probe-class).
"""

from fractions import Fraction
from itertools import product
from math import factorial


def compositions(n, parts):
    """All tuples of `parts` nonnegative ints summing to n."""
    if parts == 1:
        yield (n,)
        return
    for first in range(n + 1):
        for rest in compositions(n - first, parts - 1):
            yield (first,) + rest


def _simplex_moment(exps, measure):
    """int prod_i y_i^{exps[i]} over the (d-1)-simplex, d = len(exps).
    measure='dirichlet': normalized Dir(1,..,1) expectation.
    measure='lebesgue' : bare Lebesgue integral."""
    d = len(exps)
    num = 1
    for n in exps:
        num *= factorial(n)
    val = Fraction(num, factorial(sum(exps) + d - 1))
    if measure == "dirichlet":
        val *= factorial(d - 1)
    elif measure != "lebesgue":
        raise ValueError(measure)
    return val


def Z_definition(U, g, measure="dirichlet"):
    """Exact Z(U) for a g-component mixture, straight from the definition."""
    U = list(U)
    k = len(U)
    total = Fraction(0)
    per_state = [list(compositions(u, g)) for u in U]
    for table in product(*per_state):
        # table[v][a] = x_{v,a}; multinomial coefficient per state
        coeff = 1
        for v in range(k):
            c = factorial(U[v])
            for x in table[v]:
                c //= factorial(x)
            coeff *= c
        # column sums m_a = total allocations to component a
        m = [sum(table[v][a] for v in range(k)) for a in range(g)]
        term = Fraction(coeff) * _simplex_moment(tuple(m), measure)
        for a in range(g):
            term *= _simplex_moment(tuple(table[v][a] for v in range(k)), measure)
        total += term
    return total


if __name__ == "__main__":
    # smoke: N=0 must give 1 under the normalized measure, any g
    for g in (2, 3, 4):
        assert Z_definition([0, 0], g) == 1
    print("w1_brute smoke: Z([0,0],g)=1 for g=2,3,4  PASS")
    print("Z([2,1], g=2) =", Z_definition([2, 1], 2))
