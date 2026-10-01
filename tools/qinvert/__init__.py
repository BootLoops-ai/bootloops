"""qinvert — exact consistency certificates for published summary statistics
over partially released entity-level tables.

Given statistics an agency published EXACTLY (quantiles, Tukey fences,
min/max/range, means) and the partial dataset it released, solve for the
minimal delta families (additions with value windows, removals of released
rows) that make the recomputed statistics match exactly — with parity
constraints per quantile definition, cross-statistic stacking across
datasets, and a global no-regression gate.

Promoted and generalized from the originating analysis's inverse_membership2.py
+ membership_stacker.py (validated against published CMS Star Ratings fence
tables; receipts in the origin tree, internal reference).  See README.md for the problem
statement, the algebra, and the honest limits.

Not for person-level data or any table suppressed to protect the
confidentiality of individuals (small-cell / statistical-disclosure-control
suppression); the tool operates on organization-level releases and cannot
identify any true suppressed record.

Exact arithmetic (fractions.Fraction) for every certified statement; no
float participates in any accept/reject decision anywhere in this package.
Multiprocessing uses the spawn context only.
"""

from .targets import (ALL_DEFINITIONS, Interval, LinearCombo, MaxStat, Mean,
                      MinStat, OrderStat, Quantile, RangeStat, R_DEFINITIONS,
                      SAS_DEFINITIONS, Statistic, TukeyFence, as_fraction,
                      evaluate_targets, invert_fence_pair, targets_exact)
from .solver import DEFAULT_CAPS, solve, solve_many
from .stacker import Dataset, StackerError, stack

__version__ = "1.0.0"

REGISTER = {
    "targets": "documented-rule (SAS QNTLDEF 1-5 doc formulas; Hyndman-Fan "
               "R types 1-9); sas5 additionally validated end-to-end "
               "vs CMS K-5/K-6 fence tables; out-of-range "
               "OrderStat evaluates not-exact, never a crash",
    "solver": "planted-truth validated (all 14 definitions, removals, "
              "parity/averaging, equal-pinned-value span merge, k=30 "
              "family) + negative controls; completeness bounded by LOUD "
              "caps incl. k_cap-exhausted on every witness-free search "
              "(README 'Honest limits'); the equal-pinned-value span merge "
              "covers equal-quartile scenarios",
    "stacker": "planted-truth e2e (shared entities, coupling, removal "
               "gate) + hard no-regression gate receipts; greedy assembly "
               "certifies feasibility, NOT entity-count minimality",
}

__all__ = [
    "ALL_DEFINITIONS", "Interval", "LinearCombo", "MaxStat", "Mean",
    "MinStat", "OrderStat", "Quantile", "RangeStat", "R_DEFINITIONS",
    "SAS_DEFINITIONS", "Statistic", "TukeyFence", "as_fraction",
    "evaluate_targets", "invert_fence_pair", "targets_exact",
    "DEFAULT_CAPS", "solve", "solve_many",
    "Dataset", "StackerError", "stack",
    "REGISTER",
]
