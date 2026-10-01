"""trust.strata — identity front for the LIVE streamed-IBP certifier.

Home: the strata/ member tree beside this package — the front resolves there
by default; TRUST_STRATA_ROOT overrides for an external strata tree (a
sys.path ROOT: the directory CONTAINING the strata package).
Fronted BY IMPORT-IDENTITY: `trust.strata.load()` returns the very module
objects every other consumer imports (sha advisory only — strata is a live
tool and advances on its own; trust.verify() reports drift, never
freezes it).

ARBITRATION TABLE (fixed design decision —
TRUST packages this decision, it does not re-litigate it):

    strata   = END-TO-END certification of a streamed stratified Route-B
               F_p solve vs kira (4-slice 2-prime + weight<->integral
               dictionary crosscheck + planted faults).
    receipt  = per-(kinematic point, prime) ROW audit (lambda-witness
               against the caller's OWN re-parse of the system).
    lp_syz   = INDEPENDENT exact oracle when kira itself is in question
               (parametric syzygy route; catches basis rank-deficiency /
               non-uniqueness — the sec39/sec45/LP(431) class — which the
               other two structurally cannot see).

ONE-LINEAGE CAVEAT (certify explicitly, never present as independence):
strata and winnow (ibplapper) are ONE lineage — strata's engine.py/bank.py/
weights are VERBATIM LIFTS into winnow (winnow manual RELATED; byte-pinned
lifts re-verified in winnow's own battery). Agreement between a strata solve
and a winnow solve is a REGRESSION check of shared bytes, NOT an independent
confirmation. Independence in the TRUST triad comes only from: kira (adopted
third-party C++/Fermat), receipt core (stdlib-only, reimplementable from
WITNESS_FORMAT.md alone), and lp_syz (Singular syz + flint fmpq RREF, with
the FLINT common-substrate caveat mitigated by the stdlib fraction_oracle
leg, see trust.fraction_oracle).
"""
from . import _core

HOME = _core.LINKED["strata"]


def load():
    """Import the live `strata` package in place by identity."""
    return _core.alias_identity("strata")


def certify():
    """strata.certify (lambda-certificate check module) by identity."""
    load()
    return _core.alias_identity("strata.certify")


def fp_eliminate():
    """strata.fp_eliminate (kira-free stratified eliminator) by identity."""
    load()
    return _core.alias_identity("strata.fp_eliminate")
