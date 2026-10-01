#!/usr/bin/env python3
"""Two-engine dedup reconciler — the disagreement path.

Laws implemented:
- RAW and DEDUPED counts both reported -> N_raw + per-engine N_dedup are
  first-class in every report; must-fail plants are inserted PRE-dedup and
  must survive to the distance gate as their OWN orbit class; the
  orbit-size distribution is reported split on/off locus (diagnostic).
- TWO-STACK: every census count, RAW and DEDUPED, reproduced by an
  independent exact stack (Engine-A/Engine-B discipline).
- Any RED control kills -> every RED status here carries halt=True; a RED
  reconciliation is reported, never silently consumed.

DISAGREEMENT PATH:
  partitions equal                      -> GREEN-AGREE (two-stack GREEN)
  unequal, every extra merge's generator-word witness VERIFIES exactly
    (witness.py referee)                -> RED-RECONCILED-MISS: the engine
    lacking a witnessed merge UNDER-merged (bounded search miss); the
    reconciled partition (closure of witnessed merges) is reported for the
    audit, but a two-stack count mismatch is a RED control => halt.
  any witness FAILS                     -> RED-WITNESS-FAIL (integrity
    halt; an unwitnessable merge is never accepted).
  a plant shares a class with an honest candidate -> RED-PLANT-MERGED
    (must-fail law: control failure, halt).
"""
import witness

_INV = {"S": "s", "s": "S", "T": "t", "t": "T"}

def inv_word(word):
    out = []
    for tok in reversed(word):
        if tok in _INV:
            out.append(_INV[tok])
        elif tok[0] == "M":
            out.append("m" + tok[1:])
        else:
            out.append("M" + tok[1:])
    return tuple(out)

def _canon(part):
    return sorted(sorted(c) for c in part)

def _cls_of(part, n):
    m = [-1] * n
    for ci, c in enumerate(part):
        for i in c:
            m[i] = ci
    return m

def _extra_pairs(psrc, pother, n):
    """Index pairs merged by psrc but split by pother."""
    co = _cls_of(pother, n)
    out = []
    for c in psrc:
        for j in c[1:]:
            if co[j] != co[c[0]]:
                out.append((c[0], j))
    return out

def _closure(parts, n):
    uf = list(range(n))
    def find(x):
        while uf[x] != x:
            uf[x] = uf[uf[x]]
            x = uf[x]
        return x
    for part in parts:
        for c in part:
            for j in c[1:]:
                uf[find(c[0])] = find(j)
    cls = {}
    for i in range(n):
        cls.setdefault(find(i), []).append(i)
    return _canon(cls.values())

def _plant_audit(part, tags):
    ok, rows = True, []
    for ci, c in enumerate(part):
        if any(tags[i] == "plant" for i in c):
            pure = all(tags[i] == "plant" for i in c)
            ok = ok and pure
            for i in c:
                if tags[i] == "plant":
                    rows.append({"plant_idx": i, "class": ci,
                                 "own_class": pure})
    return ok, rows

def _hist(part, labels):
    out = {}
    for c in part:
        labs = sorted(set(labels[i] for i in c))
        lab = labs[0] if len(labs) == 1 else "mixed"
        d = out.setdefault(lab, {})
        d[str(len(c))] = d.get(str(len(c)), 0) + 1
    return out

def reconcile(states, gens, resA, resB, tags=None, labels=None):
    """Adjudicate the two independent dedup stacks. Returns the report."""
    n = len(states)
    pa, pb = _canon(resA["partition"]), _canon(resB["partition"])
    tab = witness.make_table(gens, len(states[0][0]))
    rep = {"N_raw": n, "N_dedup_A": len(pa), "N_dedup_B": len(pb),
           "truncated_A": resA.get("truncated", False),
           "truncated_B": resB.get("truncated", False),
           "engine_miss_pairs": {"A": [], "B": []}, "witness_failures": []}
    if pa == pb:
        status, part = "GREEN-AGREE", pa
    else:
        for name, psrc, pother, res, other in (
                ("A", pa, pb, resA, "B"), ("B", pb, pa, resB, "A")):
            for (i, j) in _extra_pairs(psrc, pother, n):
                w = res["words"][i] + inv_word(res["words"][j])
                if witness.verify(states[i], w, states[j], gens, tab):
                    rep["engine_miss_pairs"][other].append([i, j])
                else:
                    rep["witness_failures"].append([name, i, j])
        if rep["witness_failures"]:
            status, part = "RED-WITNESS-FAIL", None
        else:
            status, part = "RED-RECONCILED-MISS", _closure([pa, pb], n)
    if part is not None and tags is not None:
        ok, rows = _plant_audit(part, tags)
        rep["plant_own_class"], rep["plant_rows"] = ok, rows
        if not ok:
            status = "RED-PLANT-MERGED"
    if part is not None and labels is not None:
        rep["orbit_size_hist_by_label"] = _hist(part, labels)
    rep["partition"] = part
    rep["status"] = status
    rep["halt"] = status != "GREEN-AGREE"
    return rep
