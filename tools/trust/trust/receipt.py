"""trust.receipt — identity front for the winnow.receipt verification layer.

TWO registered live homes carry THE SAME receipt core BY DESIGN
(byte-identical vendor, test-enforced both ways — winnow manual GATES):

    tools/trust/receipt  the receipt MEMBER tree beside this package
                         (solver-agnostic, also runs standalone: core.py +
                          receipt.py CLI + adapters/strata.py;
                          WITNESS_FORMAT.md v1.0; TRUST_RECEIPT_ROOT
                          points the front at an external receipt tree)
    ibplapper.receipt    env-pointed via TRUST_IBPLAPPER_ROOT
                         (winnow's vendored copy — the site winnow row names
                          `winnow.receipt`; the TRUST page links it by identity)

Both are fronted BY IMPORT-IDENTITY, never copied (R1). The witness schema is
v1.0 (R5): syzygy-side provenance rides the EXTENSION slot —
`system.source` (free-text provenance, optional per WITNESS_FORMAT.md) plus
the optional `family`/`point`/label fields.

byte_identity() re-checks the two cores are byte-identical at call time
(their own test batteries enforce it; this is the cheap cross-front check).

STRATA-LOADER SUBSTRATE CAVEAT (check 2 AS EXECUTED —
named parallel to the FLINT caveat): the receipt CORE is stdlib-only and
lineage-clean, but the verify-path ROWS in the pilots and the battery come
from receipt/adapters/strata.py::load_rows (the receipt member beside
this package), which is strata's OWN
loader.load_system + loader.CoeffEvaluator — the same strata-lineage
parse/eval that stages the eliminator's rows — and that adapter imports
fp_eliminate (the eliminator) at module level. So the independence claim
stops at the CORE: a systematic loader/CoeffEvaluator bug is common-mode
across check 1's kira-side evaluation AND check 2's row re-parse (live
precedent: a k2disp float-division exactness bug recorded in
CoeffEvaluator's own docstring). Only check 3 (Singular/sympy/flint, own
parser) escapes that substrate. Never present check 2 as executed as
"no shared code" with the eliminator. Mitigations wired battery-side:
an INDEPENDENT stdlib re-parse of a sha-seeded FOREIGN SYSTEM-row subsample
(T4, no loader import, ast-based modular arithmetic) and the T7C
import-graph census that names the adapter->loader/fp_eliminate edge.
"""
import hashlib
import importlib.util
import os
import re as _re
import sys

from . import _core

HOME_TOOLS = _core.LINKED["receipt"]
# Unset TRUST_IBPLAPPER_ROOT leaves the winnow home unregistered: the
# winnow fronts then refuse by name (alias_identity fail-closed;
# byte_identity guards below) instead of failing at import.
HOME_WINNOW = _core.LINKED.get("ibplapper", "")


def core():
    """The receipt member's core.py by identity (module name
    'trust_receipt_core' to avoid squatting the generic 'core' name — the
    receipt dir is not a package)."""
    path = os.path.join(HOME_TOOLS, "core.py")
    name = "trust_receipt_core"
    if name in sys.modules:
        mod = sys.modules[name]
        if os.path.abspath(getattr(mod, "__file__", "")) != os.path.abspath(path):
            raise ImportError(f"trust identity violation: {name} shadowed")
        return mod
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    old = sys.dont_write_bytecode
    sys.dont_write_bytecode = True
    try:
        sys.modules[name] = mod
        spec.loader.exec_module(mod)
    finally:
        sys.dont_write_bytecode = old
    return mod


def winnow_core():
    """ibplapper's vendored receipt core by identity (the `winnow.receipt`
    the site names). Import path: ibplapper.receipt.core."""
    _core.alias_identity("ibplapper")
    return _core.alias_identity("ibplapper.receipt.core")


def winnow():
    """The ibplapper (winnow) package itself by identity."""
    return _core.alias_identity("ibplapper")


def byte_identity():
    """True iff the two receipt cores are byte-identical (test-enforced by
    their own batteries; re-checked here cheaply) — up to the TWO recorded
    docstring deltas (a lineage line and a design-invariant attribution;
    code bytes identical, any other difference still fails)."""
    if not HOME_WINNOW:
        raise ImportError(
            "TRUST_IBPLAPPER_ROOT unset — byte_identity needs the winnow "
            "(ibplapper) home to compare cores; set TRUST_IBPLAPPER_ROOT")

    def _norm(p):
        raw = open(p, "rb").read()
        # Collapse the two recorded delta lines to fixed tags by PATTERN
        # (either copy's wording matches; no other byte may differ).
        raw = _re.sub(rb"Lineage: checker\.py of [^\n]*",
                      b"Lineage: checker.py of <lineage-run>", raw)
        raw = _re.sub(rb"DESIGN INVARIANT \([^)\n]*2026-07-07\):",
                      b"DESIGN INVARIANT <2026-07-07>:", raw)
        return raw
    a = os.path.join(HOME_TOOLS, "core.py")
    b = os.path.join(HOME_WINNOW, "ibplapper", "receipt", "core.py")
    ha = hashlib.sha256(_norm(a)).hexdigest()
    hb = hashlib.sha256(_norm(b)).hexdigest()
    return ha == hb, {"tools_receipt_core": ha, "winnow_receipt_core": hb}
