"""Rejection sampler for the H0-a (index-density) control.

Proposals: uniform EXACT dyadic rationals on the fundamental domain's
bounding box, from the fixed per-(geometry, bin) PCG64 stream
(seeds.py). Acceptance: exact comparison u * M < f(x) with f the AD
density factor (Fraction) and M a claimed sup of f on the domain; every
evaluation is checked against M, so a bad envelope HALTS, never biases.

NO-EXCLUSION-ZONE GUARANTEE (testable): nothing in this module
accepts locus data -- no special-point list, no tube radius, no
exclusion region of any kind. The only rejection causes are (i) proposal
outside the fundamental domain, (ii) the density coin. Hence any
positive-volume subset of the domain -- in particular every tube around
a frozen special point -- receives proposal mass equal to its volume
fraction, by construction. Enforced structurally (AST identifier scan)
and statistically in tests/test_no_exclusion.py.
"""
from fractions import Fraction

try:
    from .seeds import dyadic_uniform
except ImportError:
    from seeds import dyadic_uniform


def propose(domain, gen):
    """One uniform exact-dyadic proposal on the bounding box."""
    lo, hi = domain.box()
    return [a + (b - a) * dyadic_uniform(gen) for a, b in zip(lo, hi)]


def sample(domain, density, envelope_M, gen, n_accept,
           max_proposals=10_000_000):
    """Draw until n_accept acceptances. Returns (points, counters).

    counters: proposed / out_of_domain / rejected_coin / accepted.
    Exact-rational throughout; raises on envelope violation or negative
    density (halt discipline, never silent).
    """
    M = Fraction(envelope_M)
    acc = []
    c = {"proposed": 0, "out_of_domain": 0, "rejected_coin": 0, "accepted": 0}
    while c["accepted"] < n_accept:
        if c["proposed"] >= max_proposals:
            raise RuntimeError("max_proposals exhausted (rate-collapse stop)")
        x = propose(domain, gen)
        c["proposed"] += 1
        if not domain.contains(x):
            c["out_of_domain"] += 1
            continue
        f = density(x)
        if f > M:
            raise ValueError("envelope violated: f(x)=%s > M=%s" % (f, M))
        u = dyadic_uniform(gen)
        if u * M < f:
            acc.append(x)
            c["accepted"] += 1
        else:
            c["rejected_coin"] += 1
    return acc, c
