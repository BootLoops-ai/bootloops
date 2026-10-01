#!/usr/bin/env python3
"""adapters/strata.py — kira/STRATA bridge for the receipt tool.

Everything kira-specific lives HERE, never in core.py (design invariant:
the verify core is solver-agnostic; this adapter only supplies (a) the
system rows evaluated at a (point, prime) and (b) label<->column-id
bookkeeping for human-readable index tuples).

Depends on the strata house tools (loader.py, weights.py,
fp_eliminate.py, coverage.py, corpus_bench.py) — by default the strata
member of the trust package, the sibling of this receipt member
(tools/trust/strata, i.e. ../strata from the receipt directory);
RECEIPT_STRATA_PATH overrides for an external strata tree. Semantics are
the reference ones: trivial-sector drop, ordinal-master decode, shell =
decode-class outside the solve box, survivor masters via masters-file
(sector,class) buckets, closed-vector anchoring for tuple->weight.
"""
import glob
import gzip
import json
import os
import sys

# Default: the strata member of the trust package, the sibling of this
# receipt member (tools/trust/receipt/../strata = tools/trust/strata);
# RECEIPT_STRATA_PATH overrides. Without either, the adapter refuses by
# name (generic systems go through --system-jsonl instead).
_DEFAULT_STRATA = os.path.abspath(os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "..", "strata"))
STRATA_PATH = os.environ.get("RECEIPT_STRATA_PATH", "")
if not STRATA_PATH and os.path.isdir(_DEFAULT_STRATA):
    STRATA_PATH = _DEFAULT_STRATA
if not STRATA_PATH:
    raise ImportError(
        "RECEIPT_STRATA_PATH unset and no strata member found beside this "
        "tool — the strata house tools are required for this adapter "
        "(generic systems: use --system-jsonl)")
if STRATA_PATH not in sys.path:
    sys.path.insert(0, STRATA_PATH)

from corpus_bench import eval_vec, parse_oracle_table   # noqa: E402
from coverage import assemble_closed                    # noqa: E402
from fp_eliminate import eliminate_fast                 # noqa: E402
from loader import (CoeffEvaluator, load_system,        # noqa: E402
                    load_trivial_sectors, parse_masters_file)
from weights import MASTER_WEIGHT_CEILING, decode, popcount  # noqa: E402


# ------------------------------------------------------------ system access
def load_rows(art_dir, family, p, d0, eta0):
    """Independent-parser rows for the VERIFY path (house loader semantics).
    Aborts loudly if the loader dropped empty rows (row-index integrity)."""
    sysd = load_system(art_dir, family, p, d0, eta0)
    assert sysd["counters"]["cancelled_to_empty"] == 0, \
        "row indexing mismatch risk: loader dropped cancelled-to-empty rows"
    return sysd


def parse_structural(art_dir, family):
    """SYSTEM_*.gz -> (rows_terms, exprs, sector_of) without evaluation.
    rows_terms: list of list[(expr_id, col)]; trivial-sector terms dropped.
    (Probe-harness port; used by the EMIT path only.)"""
    trivial = load_trivial_sectors(art_dir, family)
    files = sorted(glob.glob(os.path.join(art_dir, "tmp", family,
                                          "SYSTEM_*.gz")))
    assert files, f"no SYSTEM dumps under {art_dir}/tmp/{family}"
    exprs, eid = [], {}
    rows, sector_of = [], {}
    for f in files:
        with gzip.open(f, "rt") as fh:
            lines = fh.read().split("\n")
        i, n = 0, len(lines)
        while i < n:
            if lines[i] != "Eq":
                i += 1
                continue
            nt = int(lines[i + 2])
            row = []
            for j in range(nt):
                parts = lines[i + 3 + j].split()
                assert len(parts) == 6, f"bad term line {lines[i+3+j]!r}"
                coeff, w, sec = parts[1], int(parts[2]), int(parts[3])
                if sec in trivial:
                    continue
                k = eid.get(coeff)
                if k is None:
                    k = eid[coeff] = len(exprs)
                    exprs.append(coeff)
                row.append((k, w))
                sector_of[w] = sec
            i += 3 + nt
            rows.append(row)
    return rows, exprs, sector_of


def classify(sector_of, staging, art_dir, pivot_set, family):
    """Column partition per the reference corpus_bench semantics (probe port):
    pivots (given), shell (decode class outside solve box), masters_map
    (weight -> index tuple: preferred by ordinal order + survivors by
    masters-file (sector,class) buckets)."""
    box = json.load(open(os.path.join(staging, "box.json")))["solve"]
    weights = sorted(sector_of)
    cls = {w: decode(w) for w in weights}
    ordinal = [w for w in weights if w < MASTER_WEIGHT_CEILING]
    lead_set = set(pivot_set)
    shell = set()
    for w in weights:
        if w < MASTER_WEIGHT_CEILING or w in lead_set:
            continue
        dots, s = cls[w]
        t = popcount(sector_of[w])
        if dots > box["d"] or s > box["s"] or t + dots > box["r"]:
            shell.add(w)
    preferred = parse_masters_file(os.path.join(staging, "preferred"))
    fm_path = os.path.join(art_dir, "results", family, "masters")
    file_masters = (parse_masters_file(fm_path)
                    if os.path.exists(fm_path) else [])
    pref_idx = {tuple(i) for i, _ in preferred}
    masters_map = {}
    assert len(ordinal) == len(preferred), \
        f"ordinal count {len(ordinal)} != preferred {len(preferred)}"
    for w, (idx, sec) in zip(sorted(ordinal), preferred):
        assert sector_of[w] == sec, f"pref-order map broken w={w}"
        masters_map[w] = tuple(idx)
    idx_class = lambda idx: (sum(a - 1 for a in idx if a > 0),   # noqa: E731
                             -sum(a for a in idx if a < 0))
    buckets = {}
    for idx, sec in file_masters:
        if tuple(idx) in pref_idx:
            continue
        buckets.setdefault((sec, idx_class(idx)), []).append(tuple(idx))
    survivors = {}
    for w in weights:
        if w in shell or w < MASTER_WEIGHT_CEILING or w in lead_set:
            continue
        b = buckets.get((sector_of[w], cls[w]), [])
        if b:
            survivors[w] = b.pop(0)
    unclassified = [w for w in weights
                    if w >= MASTER_WEIGHT_CEILING and w not in lead_set
                    and w not in shell and w not in survivors]
    masters_map.update(survivors)
    return {"pivots": sorted(lead_set), "shell": shell,
            "masters_map": masters_map, "box": box,
            "n_unclassified": len(unclassified)}


def stage_partition(art_dir, staging, family, p, d0, eta0, sector_of):
    """One eliminator pass at the anchor slice -> pivot set +
    classification + closed rows (for tuple anchoring). A production
    partition would come from a full elimination census; this stand-in is
    measured seconds-class at the 7,135-equation reference scale."""
    sysd = load_rows(art_dir, family, p, d0, eta0)
    subs, leftover, _ = eliminate_fast(sysd["rows"], sysd["order"], p,
                                       forbid=set(sysd["masters"]))
    assert not leftover, f"leftover master relations: {len(leftover)}"
    info = classify(sector_of, staging, art_dir, set(subs), family)
    closed = assemble_closed(subs, sysd["order"], p)
    return sysd, info, closed


def anchor_table(orc_rows, info, closed, sector_of, ev, p):
    """kira_target.m rows -> [(idx_tuple, terms, target_col)] via
    closed-vector anchoring (probe/retrofit semantics, incl. duplicate-class
    disambiguation). Returns (targets, skip_census)."""
    masters_map = info["masters_map"]
    midx2w = {idx: w for w, idx in masters_map.items()}
    mm_inv = {v: k for k, v in masters_map.items()}
    piv = set(info["pivots"])
    sig_index = {}
    for w, r in closed.items():
        if w in piv:
            sig_index.setdefault(frozenset(r.items()), []).append(w)
    targets = []
    skip = {"master_or_empty": 0, "foreign_master": 0, "zero_or_denhit": 0,
            "unanchored": 0}
    for tgt, terms in orc_rows:
        tgt = tuple(tgt)
        if tgt in mm_inv or not terms:
            skip["master_or_empty"] += 1
            continue
        if any(tuple(mi) not in set(mm_inv) for _s, mi, _e in terms):
            skip["foreign_master"] += 1
            continue
        vec = eval_vec(terms, ev, midx2w, p)
        if isinstance(vec, str) or not vec:
            skip["zero_or_denhit"] += 1
            continue
        ws = sig_index.get(frozenset(vec.items()), [])
        if not ws:
            skip["unanchored"] += 1
            continue
        if len(ws) > 1:
            t_sec = sum(1 << i for i, a in enumerate(tgt) if a > 0)
            t_cls = (sum(a - 1 for a in tgt if a > 0),
                     -sum(a for a in tgt if a < 0))
            good = [w for w in ws
                    if sector_of[w] == t_sec and decode(w) == t_cls]
            ws = [min(good)] if good else [min(ws)]
        targets.append((tgt, terms, ws[0]))
    return targets, skip


def oracle_c(terms, ev, midx2w, p):
    """kira table row -> claimed c dict {master_col: value} in the receipt
    convention I_t = sum c_m I_m (eval_vec returns the NEGATED subs-convention
    vector, so c = -vec mod p). Returns dict or a status string."""
    vec = eval_vec(terms, ev, midx2w, p)
    if isinstance(vec, str):
        return vec
    return {w: (-v) % p for w, v in vec.items()}


# --------------------------------------------------- legacy witness support
def read_legacy_pair(c_json_path, masters_map, n_rows):
    """Probe-era witness pair (c_*.json + sibling lam_*.npy) -> the loaded-
    witness dict shape used by core.verify_row. Tuple-keyed c is resolved to
    column ids via masters_map (index tuple -> weight)."""
    import numpy as np
    meta = json.load(open(c_json_path))
    lam = np.load(c_json_path.replace("c_", "lam_").replace(".json", ".npy"))
    p = int(meta["p"])
    inv_map = {tuple(idx): w for w, idx in masters_map.items()}
    c = {}
    for k, v in meta["c"].items():
        idx = tuple(json.loads(k))
        w = inv_map.get(idx)
        assert w is not None, f"unresolvable master tuple {idx} in {c_json_path}"
        c[w] = int(v) % p
    nz = [int(i) for i in lam.nonzero()[0]]
    return {"p": p, "target_col": int(meta["weight"]), "c": c,
            "lam_idx": nz, "lam_val": [int(lam[i]) % p for i in nz],
            "n_rows": int(n_rows) if n_rows else len(lam),
            "meta": {"legacy": True, "point": {"d": meta.get("d0"),
                                               "eta": meta.get("eta0")}},
            "path": c_json_path}
