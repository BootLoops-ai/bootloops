"""Winnow — standalone upgraded-Laporta F_p eliminator library with
native receipts (package `ibplapper`; never bare `lapper` — registry squat).

    import ibplapper as lap
    res = lap.eliminate(lap.System(rows, p, order, forbid),
                        lap.Schedule(policy="B2FT", bank=lap.Bank(path)),
                        witnesses="retrofit", witness_targets=[...])

Sparse F_p rows + order/stratum/forbid hooks in; closed substitutions +
master-relations + per-row lambda-witnesses + audit ledger out. Everything
upstream of the rows (seeding, symmetries, dictionaries, kinematics) and
downstream of the per-(point,prime) tables (CRT, rational reconstruction)
stays caller-side — see README.md for the measured reasons.

Regression bar: byte-equal interface-table (T) bank replays + the 84/84
retrofit witness battery + the recorded wrong-table MUST-FAIL. Verifier =
the vendored RECEIPT component (imported, not forked).
"""
from .api import (Bank, CensoredError, Result, Schedule, System,  # noqa: F401
                  eliminate)
from .ledger import LIB_VERSION as __version__  # noqa: F401
from . import adapters, certify, coverage, engine, receipt, witness  # noqa: F401,E501
