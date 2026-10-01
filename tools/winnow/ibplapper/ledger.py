"""ledger.py — the Winnow audit ledger.

One JSON-serializable dict per eliminate() call:
  - library version, policy + config echo, prime, caps;
  - input row-set sha256 (receipt.core.system_fingerprint — the same
    canonical hash the witness sidecars carry);
  - per-stratum Stats (fill, peak_live_nnz/rows, ops, pivots, walls),
    pass-down counts, interface-row counts (all straight from
    engine.stratified_solve's stats_per);
  - banked-block sha256s (from the Bank index entries appended this run);
  - leftover master-relations census (count + support sample —
    NEVER silently absorbed);
  - witness summary (retrofit mode);
  - censor record when a cap fired (partial ledger, loud — the
    print-loud-rc0 rule: a censored run must never look complete).

Append-only file form (JSONL, one ledger per line) = the RECEIPT sidecar of
a run directory.
"""
import json
import time

from .receipt.core import system_fingerprint

LIB_VERSION = "1.0.0"


def new_ledger(*, p, n_rows, n_cols, policy, caps, forbid_n, strata_n,
               rows_fingerprint, meta_echo=None):
    return {
        "tool": "ibplapper",
        "version": LIB_VERSION,
        "ts": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "p": int(p),
        "n_rows_in": int(n_rows),
        "n_cols": int(n_cols),
        "policy": policy,
        "caps": caps or None,
        "n_forbid": int(forbid_n),
        "n_strata": int(strata_n),
        "input_rows_sha256": rows_fingerprint,
        "meta_echo": meta_echo,
        "per_stratum": {},          # stratum -> Stats dict (+pass-down etc.)
        "banked_blocks": [],        # [{key, sha256, n_rows}] appended this run
        "leftover": None,           # filled at close
        "witnesses": None,          # filled by witness stage
        "censored": None,           # filled if a cap fires
        "wall_s": None,
    }


def fingerprint_rows(rows, p):
    """Canonical sha256 of the input row set (row order as given, columns
    sorted within a row) — identical convention to the witness sidecars."""
    return system_fingerprint(rows, p)


def record_leftover(ledger, leftover, sample=5):
    ledger["leftover"] = {
        "n": len(leftover),
        "support_sample": [sorted(r) for r in leftover[:sample]],
        "note": ("relations among masters / rank-deficiency evidence "
                 "(sec39 class) — reported, never absorbed") if leftover
                else None,
    }


def record_censor(ledger, stratum, kind, measured, cap):
    ledger["censored"] = {"stratum": stratum, "cap": kind,
                          "measured": measured, "limit": cap,
                          "note": "loud CensoredError raised; ledger PARTIAL"}


def save(ledger, path):
    """Append the ledger as one JSONL line (append-only sidecar form)."""
    with open(path, "a") as fh:
        fh.write(json.dumps(ledger) + "\n")
    return path
