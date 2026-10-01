#!/usr/bin/env python3
r"""nestor.receipts -- JSON receipts: per-level walls, agreement digits,
refusals; atomic checkpoint books.

Every nestor run can leave a machine-readable receipt: what was integrated,
at what levels/dps, how long each rung took, what agreed with what, and --
when a run refuses -- the NAMED refusal (class + message + exit code)
instead of a value.  A refusal receipt is a first-class artifact: the
absence of a number is a result.

Atomic-write pattern:
json -> <path>.tmp -> os.replace(path); checkpoint book = {key: record},
resume = skip.

"""
from __future__ import annotations

import json
import os
import time

NESTOR_VERSION = "1.0.0"

# integrity manifest of the engine modules, stamped into every receipt:
# which bytes produced this number (computed at import, cannot go stale)
import hashlib as _hashlib

_HERE = os.path.dirname(os.path.abspath(__file__))
SRC_SHAS = {
    _f: _hashlib.sha256(
        open(os.path.join(_HERE, _f), "rb").read()).hexdigest()
    for _f in ("oracle.py", "farm.py", "ladder.py", "singcheck.py")
}


def new_receipt(kind, config=None):
    """Fresh receipt skeleton for one nestor evaluation."""
    return {
        "tool": "nestor",
        "version": NESTOR_VERSION,
        "kind": kind,                    # 'integrate' | 'farm' | 'selftest'...
        "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "config": config or {},
        "src_shas": dict(SRC_SHAS),
        "singcheck": None,
        "rungs": [],                     # per-level: level/dps/wall_s/
                                         #   value_str/agreement_digits
        "value": None,                   # set ONLY on success
        "agreed_digits": None,
        "wall_s": None,
        "refusal": None,                 # set ONLY on refusal
    }


def add_rung(receipt, **rung):
    """Append one per-level rung record (level, dps, wall_s, value_str,
    agreement_digits -- whatever the stage measured)."""
    receipt["rungs"].append(rung)
    return receipt


def add_refusal(receipt, exc):
    """Record a NAMED refusal; the receipt then carries NO value."""
    receipt["refusal"] = {
        "class": type(exc).__name__,
        "exit_code": getattr(exc, "exit_code", 3),
        "message": str(exc),
    }
    receipt["value"] = None
    receipt["agreed_digits"] = None
    return receipt


def write(receipt, path):
    """Atomic write (ported tmp+os.replace pattern)."""
    d = os.path.dirname(os.path.abspath(path))
    if d and not os.path.isdir(d):
        os.makedirs(d, exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(receipt, fh, indent=1, default=str)
    os.replace(tmp, path)
    return path


def checkpoint_load(path):
    """Load a checkpoint book ({key: record}); missing file = empty book
    (ported: oracle_par.main resume logic)."""
    if path and os.path.exists(path):
        with open(path) as fh:
            return json.load(fh)
    return {}


def checkpoint_save(path, book):
    """Atomic book save (ported tmp+os.replace)."""
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(book, fh, indent=1, default=str)
    os.replace(tmp, path)
    return path
