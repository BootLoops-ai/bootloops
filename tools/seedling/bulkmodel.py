"""Seed-box bulk combinatorics + two-phase memory preflight model.

Closed form (validated against an 11-cell measured staging matrix on a
3-loop 15-propagator family; completed cells match RSS ratios to 5-15%):

  N_sigma(r,s,d) = C(min(d, r-t)+t, t) * C(s+m, m),   m = n_sp - t
  bulk B = n_ops * sum_sigma N_sigma,                 n_ops = L*(L+E)

Two-phase memory model (calibrated on the same family):
  select-phase RSS ~ C_SEL * B      (equation skeletons; the "load floor")
  generate-phase RSS ~ A_GEN + C_GEN * B   (coefficients resident)

Calibration provenance: C_SEL from one cell (35 GB @ 1207M eqns) predicted the
independently measured 156/151 GB monolith floor to 4-7% on two machines.
Constants are family-portable HYPOTHESES until re-validated on a second family —
the preflight prediction should be treated as an estimate with ~2x safety factor.
"""
from math import comb

C_SEL = 29.0e-9   # GB per equation (skeletons)
C_GEN = 398.0e-9  # GB per equation (coefficients resident)
A_GEN = 1.7       # GB offset


def seeds_per_sector(t, n_sp, r, s, d):
    if s < 0 or r < 0 or d < 0 or t < 0:
        # degenerate-input guard: negative
        # support is a caller bug; comb() would raise an opaque error or,
        # worse, a silent 0 could poison a preflight. Fail loudly.
        raise ValueError(f"seeds_per_sector: negative support "
                         f"(t={t}, r={r}, s={s}, d={d})")
    if t > n_sp:
        raise ValueError(f"seeds_per_sector: t={t} exceeds n_sp={n_sp}")
    dhat = min(d, r - t)
    if dhat < 0:
        return 0
    m = n_sp - t
    return comb(dhat + t, t) * comb(s + m, m)


def bulk(census, n_sp, n_ops, r, s, d):
    """census: iterable of t values (one per nontrivial sector)."""
    return n_ops * sum(seeds_per_sector(t, n_sp, r, s, d) for t in census)


def preflight(census, n_sp, n_ops, r, s, d):
    """Predict staging memory (GB). Returns dict; treat as estimate (~2x band)."""
    B = bulk(census, n_sp, n_ops, r, s, d)
    return {
        "bulk_eqns": B,
        "select_floor_GB": C_SEL * B,
        "generate_peak_GB": A_GEN + C_GEN * B,
    }
