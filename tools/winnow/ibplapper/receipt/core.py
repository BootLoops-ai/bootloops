#!/usr/bin/env python3
"""core.py — solver-agnostic lambda-witness verification core.

DESIGN INVARIANT: this module makes ZERO
assumptions about how the reduction table was produced. It never imports
solver or adapter code (kira/B2FT/strata). Its entire world is:

    system rows   R_i : list of sparse F_p rows (dict col_id -> value mod p)
    witness       (target col t, coefficient dict c, multiplier vector lam)
    the identity  sum_i lam_i R_i  ==  e_t - sum_m c_m e_m   (mod p)

If the identity holds, the row  I_t = sum_m c_m I_m  is an exact linear
consequence of the system at this (point, prime) — whatever solver claimed it.
This is what lets ANY future reduction engine be a drop-in replacement:
correctness is checked per row against the system itself, not against trust
in the producer (see README: the confident-wrong-table exhibit).

Scope: certificates are per-(kinematic point, prime). No cross-point claim.
Symbolic tables are out of scope.

Measured record: 84/84 archived rows verified at 2 primes (0.96 ms/row mean);
the 6-class sabotage battery is replayed in tests/.
"""
import hashlib
import json
import os

WITNESS_VERSION = "1.0"

# rc discipline (tests/test_core.py pins this):
RC_OK = 0          # every row verified PASS
RC_FAIL = 1        # at least one row FAILED verification
RC_MALFORMED = 2   # input could not be interpreted (bad JSON/fields/files)


class MalformedInput(Exception):
    """Raised when a witness/table/system file cannot be interpreted."""


# --------------------------------------------------------------- witness IO
_REQUIRED = ("receipt_version", "p", "target", "c", "lam")


def load_witness(path):
    """Load + structurally validate a v1 witness JSON. Raises MalformedInput.

    Returns dict with normalized fields:
      p:int, target_col:int, c: dict int->int, lam_idx: list int,
      lam_val: list int, n_rows:int, meta: raw dict.
    """
    try:
        with open(path) as fh:
            w = json.load(fh)
    except (OSError, json.JSONDecodeError) as e:
        raise MalformedInput(f"{path}: unreadable witness ({e})")
    if not isinstance(w, dict):
        raise MalformedInput(f"{path}: witness is not a JSON object")
    missing = [k for k in _REQUIRED if k not in w]
    if missing:
        raise MalformedInput(f"{path}: missing fields {missing}")
    if str(w["receipt_version"]).split(".")[0] != "1":
        raise MalformedInput(
            f"{path}: unsupported receipt_version {w['receipt_version']!r}")
    try:
        p = int(w["p"])
        tcol = int(w["target"]["col"])
        c = {int(k): int(v) % p for k, v in w["c"].items()}
        lam = w["lam"]
        n_rows = int(lam["n_rows"])
        idx = [int(i) for i in lam["idx"]]
        val = [int(v) % p for v in lam["val"]]
    except (KeyError, TypeError, ValueError) as e:
        raise MalformedInput(f"{path}: bad field content ({e})")
    if p < 2:
        raise MalformedInput(f"{path}: p={p} is not a modulus")
    if len(idx) != len(val):
        raise MalformedInput(f"{path}: lam idx/val length mismatch")
    if any(i < 0 or i >= n_rows for i in idx):
        raise MalformedInput(f"{path}: lam index out of range [0,{n_rows})")
    if len(set(idx)) != len(idx):
        raise MalformedInput(f"{path}: duplicate lam indices")
    return {"p": p, "target_col": tcol, "c": c, "lam_idx": idx,
            "lam_val": val, "n_rows": n_rows, "meta": w, "path": path}


def write_witness(path, *, p, target_col, c, lam_idx, lam_val, n_rows,
                  family=None, point=None, target_label=None, labels=None,
                  system=None):
    """Serialize a v1 witness (see WITNESS_FORMAT.md)."""
    w = {
        "receipt_version": WITNESS_VERSION,
        "kind": "lambda-witness",
        "identity":
            "sum_i lam[i]*R_i == e_target - sum_m c[m]*e_m (mod p)",
        "p": int(p),
        "target": {"col": int(target_col)},
        "c": {str(int(k)): int(v) % p for k, v in c.items() if int(v) % p},
        "lam": {"n_rows": int(n_rows),
                "idx": [int(i) for i in lam_idx],
                "val": [int(v) % p for v in lam_val]},
    }
    if family is not None:
        w["family"] = family
    if point is not None:
        w["point"] = {k: int(v) for k, v in point.items()}
    if target_label is not None:
        w["target"]["label"] = list(target_label)
    if labels:
        w["labels"] = {str(int(k)): list(v) for k, v in labels.items()}
    if system:
        w["system"] = system
    with open(path, "w") as fh:
        json.dump(w, fh)
    return w


# ------------------------------------------------------------ system access
def load_system_jsonl(path):
    """Generic system snapshot: one JSON object per line, {"<col>": val, ...}
    (row order = line order = the order lam indexes). Returns list of dict."""
    rows = []
    try:
        with open(path) as fh:
            for ln, line in enumerate(fh):
                line = line.strip()
                if not line:
                    continue
                obj = json.loads(line)
                if not isinstance(obj, dict):
                    raise MalformedInput(
                        f"{path}:{ln + 1}: row is not an object")
                rows.append({int(k): int(v) for k, v in obj.items()})
    except (OSError, json.JSONDecodeError, TypeError, ValueError) as e:
        raise MalformedInput(f"{path}: unreadable system snapshot ({e})")
    return rows


def system_fingerprint(rows, p):
    """Canonical sha256 of the evaluated system at one (point, prime):
    row order as given, columns sorted within a row."""
    h = hashlib.sha256()
    h.update(f"p={p}\n".encode())
    for i, row in enumerate(rows):
        h.update(f"{i}:".encode())
        for c in sorted(row):
            h.update(f"{c}={row[c] % p};".encode())
        h.update(b"\n")
    return h.hexdigest()


# ------------------------------------------------------------------- verify
def verify_row(rows, wit, table_row=None, check_fingerprint=None):
    """Verify ONE witness against the system rows (and optionally against a
    claimed table row). Pure function of its inputs; no solver knowledge.

    rows: list of dict col->val (mod p), in the row order lam indexes.
    wit:  dict from load_witness().
    table_row: optional dict col->val = the CLAIMED c for this target.
               If given, the witness must both hold against the system AND
               carry exactly this c — a certified-but-different c means the
               claimed table row is WRONG (the wrong-table exhibit mode).
    check_fingerprint: optional fingerprint of `rows`; compared to the
               witness's recorded system fingerprint when both exist.

    Returns (passed: bool, detail: dict).
    """
    p = wit["p"]
    detail = {"target_col": wit["target_col"], "p": p,
              "lam_nnz": len(wit["lam_idx"])}

    if wit["n_rows"] != len(rows):
        detail["reason"] = (f"lam n_rows {wit['n_rows']} != system rows "
                            f"{len(rows)}")
        return False, detail

    recorded = (wit["meta"].get("system") or {}).get("fingerprint")
    if check_fingerprint and recorded and check_fingerprint != recorded:
        detail["reason"] = "system fingerprint mismatch"
        return False, detail

    # residual r = sum_i lam_i R_i (sparse dict accumulation)
    r = {}
    for i, li in zip(wit["lam_idx"], wit["lam_val"]):
        li %= p
        if not li:
            continue
        for ccol, v in rows[i].items():
            r[ccol] = (r.get(ccol, 0) + li * v) % p
    r = {ccol: v for ccol, v in r.items() if v}

    expect = {wit["target_col"]: 1}
    for m, cm in wit["c"].items():
        cm %= p
        if cm:
            expect[m] = (-cm) % p
    ok = (r == expect)
    if not ok:
        detail["reason"] = "identity failed: residual != e_t - sum c_m e_m"
        detail["residual_nnz"] = len(r)
        return False, detail

    if table_row is not None:
        claimed = {int(k): int(v) % p for k, v in table_row.items()
                   if int(v) % p}
        if claimed != wit["c"]:
            detail["reason"] = ("claimed table row != certified row "
                                "(witness holds; the TABLE is wrong)")
            diff = sorted(m for m in set(claimed) | set(wit["c"])
                          if claimed.get(m) != wit["c"].get(m))
            detail["n_differing_entries"] = len(diff)
            detail["first_diff_col"] = diff[0] if diff else None
            return False, detail

    return True, detail


def verify_many(rows, witnesses, table=None, fingerprint=None):
    """Verify a list of loaded witnesses. table: optional dict
    target_col -> {col: val}. Returns (rc, report dict)."""
    results, n_pass = [], 0
    for wit in witnesses:
        trow = table.get(wit["target_col"]) if table else None
        ok, det = verify_row(rows, wit, table_row=trow,
                             check_fingerprint=fingerprint)
        det["witness"] = os.path.basename(wit.get("path", "?"))
        det["pass"] = bool(ok)
        results.append(det)
        n_pass += bool(ok)
    rc = RC_OK if n_pass == len(results) and results else RC_FAIL
    if not results:
        rc = RC_MALFORMED
    return rc, {"n_rows": len(results), "n_pass": n_pass,
                "all_pass": n_pass == len(results) and bool(results),
                "results": results}
