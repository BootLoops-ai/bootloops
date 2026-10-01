"""
gatekeeper — reusable data-integrity + held-out-CV honesty layer for
amflow/oracle per-master-per-point JSON dumps.

Two cheap honesty wins:
  1. scan_integrity(dirs)   — hash-based dedup audit (catches byte-identical
                              dumps mis-collected under two master labels).
  2. heldout_certify(...)   — leave-one-out CV + PSLQ rational fit; a master
                              is CERTIFIED only if every held-out point
                              matches the N-1-point fit to ≥ gate digits.

Dependency-light: stdlib + mpmath only.
"""
__all__ = ["scan_integrity", "IntegrityReport",
           "heldout_certify", "load_oracle_values", "CVResult",
           "lyndon_basis", "auto_condition", "pslq_with_lll",
           "amf_laurent", "amf_arb", "amf_ball"]


def __getattr__(name):
    # lazy import so `python -m gatekeeper.oracle_integrity` doesn't
    # double-import the submodule
    if name in ("scan_integrity", "IntegrityReport"):
        from .oracle_integrity import scan as scan_integrity, IntegrityReport
        return {"scan_integrity": scan_integrity,
                "IntegrityReport": IntegrityReport}[name]
    if name in ("heldout_certify", "load_oracle_values", "CVResult",
                "lyndon_basis", "auto_condition", "pslq_with_lll"):
        from . import heldout_cv as _h
        return getattr(_h, name)
    if name in ("amf_laurent", "amf_arb", "amf_ball"):
        from . import amf_result as _a
        return getattr(_a, name)
    raise AttributeError(name)
