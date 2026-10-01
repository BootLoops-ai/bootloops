"""Test-harness oracle plumbing — parse_oracle_table, eval_vec, idx_class,
plus column-partition/anchoring helpers over ibplapper.adapters.kira. This
is HARNESS code: it never enters the library (the weight<->integral
dictionary is caller-side by design)."""
import json
import os
import re

from ibplapper.adapters.kira import (CoeffEvaluator, parse_masters_file)
from ibplapper.adapters.weights import (MASTER_WEIGHT_CEILING, decode,
                                        popcount)

_ROW_HEAD = re.compile(r"^\s*(\w+)\[([-\d, ]+)\]\s*->\s*(.*)$")
_TERM = re.compile(r"([+-])?\s*(\w+)\[([-\d, ]+)\](?:\*\((.*)\))?\s*,?\s*$")


def parse_oracle_table(path):
    """kira_target.m -> list of (target_idx, [(sign, master_idx, expr), ...])."""
    rows, cur = [], None
    for raw in open(path):
        line = raw.strip()
        if line in ("{", "}", ""):
            continue
        if line == ",":
            cur = None
            continue
        m = _ROW_HEAD.match(line)
        if m:
            tgt = tuple(int(x) for x in m.group(2).split(","))
            cur = (tgt, [])
            rows.append(cur)
            rest = m.group(3).strip().rstrip(",")
            if rest in ("", "+", "0"):
                continue
            line = rest
        if cur is None:
            raise ValueError(f"orphan line in {path}: {raw!r}")
        body = line.rstrip(",").strip()
        if body in ("0", "+ 0", "+0", "- 0", "-0") or not body:
            continue
        t = _TERM.match(body)
        if not t:
            raise ValueError(f"unparsed term in {path}: {raw!r}")
        sign = -1 if t.group(1) == "-" else 1
        midx = tuple(int(x) for x in t.group(3).split(","))
        cur[1].append((sign, midx, t.group(4) or "1"))
    return rows


def idx_class(idx):
    return sum(a - 1 for a in idx if a > 0), -sum(a for a in idx if a < 0)


def eval_vec(terms, ev, midx2w, p):
    """(sign, midx, expr) list -> NEGATED F_p vector over master weights
    (matching the subs convention pivot + rest = 0), or ('denhit'|'nomaster')."""
    vec = {}
    for sign, midx, expr in terms:
        try:
            v = ev(expr)
        except ZeroDivisionError:
            return "denhit"
        mw = midx2w.get(midx)
        if mw is None:
            return "nomaster"
        vec[mw] = (vec.get(mw, 0) - sign * v) % p
    return {k: v for k, v in vec.items() if v}


def classify(sector_of, staging, art_dir, pivot_set, family):
    """Column partition: pivots (given), shell (decode class
    outside the solve box), masters_map (preferred ordinals + survivor
    buckets)."""
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


def anchor_table(orc_rows, info, closed, sector_of, ev, p):
    """Oracle rows -> [(idx_tuple, terms, target_col)] via closed-vector
    anchoring, incl. duplicate-class disambiguation (receipt strata lift)."""
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
            t_cls = idx_class(tgt)
            good = [w for w in ws
                    if sector_of[w] == t_sec and decode(w) == t_cls]
            ws = [min(good)] if good else [min(ws)]
        targets.append((tgt, terms, ws[0]))
    return targets, skip


def stage_partition(art_dir, staging, family, p, d0, e0):
    """Load + eliminate (B2F) + classify + close — the anchoring stand-in
    (declared probe semantics, seconds-class at k2disp scale)."""
    import ibplapper as lap
    from ibplapper.adapters import kira as kad
    from ibplapper.coverage import assemble_closed
    system, sysd = kad.load(art_dir, family, p, d0, e0)
    res = lap.eliminate(system, lap.Schedule(policy="B2F"))
    assert not res.leftover, f"leftover master relations: {len(res.leftover)}"
    info = classify(sysd["sector_of"], staging, art_dir, set(res.subs),
                    family)
    return system, sysd, info, res.subs
