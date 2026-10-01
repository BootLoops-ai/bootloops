"""nullcal — null calibration of mixed-vs-pure Bayes factors.

Any claim "object X shows BF(mixed : best pure) = B" is read against the
distribution of the SAME statistic over KNOWN-PURE objects under the same
model and conventions (leave-one-out wherever profiles are trained on the
known objects). A mixed-vs-pure factor inside the known-pure tail is not
evidence of mixture: the mixture's extra freedom also absorbs
object-to-object rate variation that the iid idealisation leaves out.

Reference instance (Federalist, exact): No. 55's blend factors
10^2.0-10^4.0 are reached by 3, 6, and 9 of the 65 known single-author papers
under the three evaluation word lists (fixtures_nullcal.json).

Two layers:
  calibrate(claim_bf, known_bfs)   pure arithmetic on precomputed log10 BFs
  blend_vs_pure_bf(U, gmix=2, Z=?) exact log10 BF(g_mix : 1) for exchangeable
      count vectors. Inside the mixalot package the engine resolves through
      the integrity-checked loader (engines.load("bigg")); standalone callers
      MUST pass Z explicitly — there is no path fallback.
"""
from fractions import Fraction
from math import log10


class EngineUnavailable(RuntimeError):
    """Raised when no exact engine is resolvable: outside the mixalot
    package, pass Z=callable(U, g) -> Fraction explicitly."""


def calibrate(claim_bf, known_bfs):
    """Locate claim_bf in the known-pure distribution of the same statistic.

    Returns counts, quantiles and a verdict:
      INSIDE-NULL  claim <= q95 of the known-pure distribution
      TAIL         q95 < claim <= max
      ABOVE-NULL   claim > every known-pure value
    """
    v = sorted(float(x) for x in known_bfs)
    n = len(v)
    if n == 0:
        raise ValueError('empty known_bfs')
    claim = float(claim_bf)
    n_at_or_above = sum(1 for x in v if x >= claim)

    def q(p):
        return v[min(n - 1, int(p * n))]

    if claim <= q(0.95):
        verdict = 'INSIDE-NULL'
    elif claim <= v[-1]:
        verdict = 'TAIL'
    else:
        verdict = 'ABOVE-NULL'
    return {
        'n_known': n,
        'claim': claim,
        'n_at_or_above': n_at_or_above,
        'percentile': 100.0 * (n - n_at_or_above) / n,
        'median': q(0.5),
        'q95': q(0.95),
        'max': v[-1],
        'verdict': verdict,
    }


def _log10_fraction(fr):
    assert fr > 0
    return log10(fr.numerator) - log10(fr.denominator)


def _resolve_engine():
    try:
        from .engines import load          # package context only
    except ImportError as e:
        raise EngineUnavailable(
            'REFUSED [nullcal.engine-unresolved]: not running inside the '
            'mixalot package and no Z callable was passed') from e
    return load('bigg').Z_bigg


def _z1_closed(U):
    from math import factorial
    k, N = len(U), sum(U)
    z = Fraction(factorial(k - 1))
    for u in U:
        z *= factorial(u)
    return z / factorial(N + k - 1)


def blend_vs_pure_bf(U, gmix=2, Z=None):
    """Exact log10 BF(gmix : 1) for one exchangeable count vector."""
    if Z is None:
        Z = _resolve_engine()
    zg = Z(list(U), gmix)
    try:
        z1 = Z(list(U), 1)
    except Exception:
        z1 = _z1_closed(U)
    return _log10_fraction(Fraction(zg)) - _log10_fraction(Fraction(z1))


def calibrate_objects(U_claim, known_Us, gmix=2, Z=None):
    """End-to-end: exact BFs for the claim object and every known-pure
    object, then calibrate. For frozen-profile blends (Federalist-style),
    compute the leave-one-out BFs with the corpus machinery and call
    calibrate() directly instead."""
    known = [blend_vs_pure_bf(u, gmix, Z) for u in known_Us]
    return calibrate(blend_vs_pure_bf(U_claim, gmix, Z), known)
