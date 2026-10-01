"""qinvert.stacker — cross-dataset stacking of delta families.

An agency publishes statistics for MANY datasets (measures, years, cohorts)
computed from overlapping entity sets it only partially releases.  A single
suppressed entity contributes one value to each dataset it participates in.
This module finds a small shared set of added entities (per-dataset values,
class-gated participation) plus removals of released entities that makes
every missed published target exact while provably leaving every
already-exact dataset exactly unchanged (hard no-regression gate).

The org-class / participation / value-coupling model is pluggable via
predicates:

  participation(cls, dataset_name) -> bool   may an entity of class `cls`
                                             carry a value in that dataset?
  couplings = [(super_name, sub_name), ...]  an entity carrying a value in
                                             `sub` must carry the SAME value
                                             in `super` (verify a coupling
                                             in-data before enforcing it)

Exactness: all accept/reject decisions in fractions.Fraction; the final
battery re-evaluates every dataset from scratch through targets.py.
"""

from collections import Counter
from fractions import Fraction as F

from .targets import as_fraction, targets_exact
from .solver import solve, solve_many

__all__ = ["Dataset", "stack", "StackerError"]


class StackerError(RuntimeError):
    pass


class Dataset:
    """One released dataset: rows (entity_id, value) + published targets."""

    def __init__(self, name, rows, targets, register=None,
                 value_lo=None, value_hi=None):
        self.name = name
        self.rows = [(eid, as_fraction(v)) for eid, v in rows]
        self.targets = targets
        self.register = as_fraction(register) if register is not None else None
        self.value_lo = as_fraction(value_lo) if value_lo is not None else None
        self.value_hi = as_fraction(value_hi) if value_hi is not None else None

    def values(self, removed_eids=frozenset(), adds=()):
        return sorted([v for e, v in self.rows if e not in removed_eids] +
                      [as_fraction(a) for a in adds])

    def exact(self, removed_eids=frozenset(), adds=()):
        return targets_exact(self.values(removed_eids, adds), self.targets)

    def entity_ids(self):
        return {e for e, _ in self.rows}


def _family_values(fam):
    return sorted(spec["value"] for spec in fam["adds"])


def _sub_multiset(small, big):
    cs, cb = Counter(small), Counter(big)
    return all(cb[v] >= c for v, c in cs.items())


def _pick_removal_entities(ds, rem_values, datasets, already_removed):
    """Entities of `ds` realizing rem_values such that every EXACT dataset
    stays exact when they vanish everywhere (the removal-is-globally-safe
    gate).  Returns list of eids or None."""
    pool = {}
    for e, v in ds.rows:
        pool.setdefault(v, []).append(e)
    chosen = []
    exact_now = [d for d in datasets
                 if d.exact(removed_eids=set(already_removed))]
    for v, cnt in sorted(Counter(rem_values).items()):
        got = 0
        for e in pool.get(v, []):
            if e in already_removed or e in chosen:
                continue
            trial = set(already_removed) | set(chosen) | {e}
            if all(d.exact(removed_eids=trial) for d in exact_now
                   if e in d.entity_ids()):
                chosen.append(e)
                got += 1
                if got == cnt:
                    break
        if got < cnt:
            return None
    return chosen


def stack(datasets, classes=("ENT",), participation=None, couplings=(),
          k_cap=8, rem_cap=0, max_witnesses=8, workers=0, solver_caps=None):
    """Assemble a shared delta set across datasets.

    Returns {"entities": [{"id","cls","values":{ds: str}}...],
             "removed_entities": [eid...],
             "battery": [{"name","before_exact","after_exact"}...],
             "exact_before", "exact_after", "evaluable",
             "unsolved": [name...], "notes": [...], "regressions": []}

    Raises StackerError if the assembled delta would regress an
    already-exact dataset (this is a hard gate; by construction it should
    never fire — the exception is the receipt that the gate is real).
    """
    if participation is None:
        participation = lambda cls, ds_name: True
    notes = []
    status = {d.name: ("exact" if d.exact() else "miss") for d in datasets}
    misses = [d for d in datasets if status[d.name] == "miss"]
    unsolved = []

    # ---- phase 1: per-dataset solve --------------------------------------
    jobs = []
    for d in misses:
        jobs.append({"job_id": d.name, "values": [v for _, v in d.rows],
                     "targets": d.targets, "register": d.register,
                     "value_lo": d.value_lo, "value_hi": d.value_hi,
                     "k_cap": k_cap, "rem_cap": rem_cap,
                     "max_witnesses": max_witnesses, "caps": solver_caps})
    results = solve_many(jobs, workers=workers) if workers > 1 else \
        [None] * len(jobs)
    if workers <= 1:
        results = [None] * len(jobs)
        for i, j in enumerate(jobs):
            r = solve(j["values"], j["targets"], register=j["register"],
                      value_lo=j["value_lo"], value_hi=j["value_hi"],
                      k_cap=k_cap, rem_cap=rem_cap,
                      max_witnesses=max_witnesses, caps=solver_caps)
            r["job_id"] = j["job_id"]
            results[i] = r
    res_by_name = {r["job_id"]: r for r in results}

    # ---- phase 2: removal entity selection (globally gated) --------------
    removed = []
    dsmap = {d.name: d for d in datasets}
    for d in misses:
        r = res_by_name[d.name]
        if r["k"] is None:
            unsolved.append(d.name)
            continue
        add_only = [f for f in r["families"] if not f["removed"]]
        if add_only:
            continue  # no removals needed for this dataset
        picked = False
        for fam in r["families"]:
            eids = _pick_removal_entities(d, fam["removed"], datasets,
                                          removed)
            if eids is not None:
                removed.extend(eids)
                notes.append(f"removals for {d.name}: {eids} "
                             f"(values {[str(v) for v in fam['removed']]})")
                picked = True
                break
        if not picked:
            notes.append(f"{d.name}: no globally safe removal entities")
            unsolved.append(d.name)
    removed = list(dict.fromkeys(removed))

    # ---- phase 3: re-solve misses on the removal-reduced bases -----------
    final = {}   # name -> list of families (add-only on reduced base)
    for d in misses:
        if d.name in unsolved:
            continue
        r = solve(d.values(removed_eids=set(removed)), d.targets,
                  register=d.register, value_lo=d.value_lo,
                  value_hi=d.value_hi, k_cap=k_cap, rem_cap=0,
                  max_witnesses=max_witnesses, caps=solver_caps)
        if r["k"] is None:
            unsolved.append(d.name)
            notes.append(f"{d.name}: unsolved after removals")
            continue
        final[d.name] = r["families"]

    # ---- phase 4: coupling (sub multiset rides super's shared entities) --
    chosen = {}          # name -> chosen multiset (list of Fractions)
    coupled_shared = {}  # (sup, sub) -> shared multiset S
    for sup_name, sub_name in couplings:
        sup_miss, sub_miss = sup_name in final, sub_name in final
        if not (sup_miss or sub_miss):
            continue
        if sup_miss and sub_miss:
            match = None
            for fam_m in final[sup_name]:
                M = _family_values(fam_m)
                for fam_s in final[sub_name]:
                    S = _family_values(fam_s)
                    if _sub_multiset(S, M):
                        match = (M, S)
                        break
                if match:
                    break
            if match is None:
                # force-include: pre-add each sub witness to super's base
                # and solve the remainder (fallback)
                sup_ds = dsmap[sup_name]
                for fam_s in final[sub_name]:
                    S = _family_values(fam_s)
                    r = solve(sup_ds.values(removed_eids=set(removed),
                                            adds=S),
                              sup_ds.targets, register=sup_ds.register,
                              value_lo=sup_ds.value_lo,
                              value_hi=sup_ds.value_hi,
                              k_cap=max(0, k_cap - len(S)), rem_cap=0,
                              max_witnesses=2, caps=solver_caps)
                    if r["k"] is not None:
                        M = sorted(S + _family_values(r["families"][0]))
                        match = (M, S)
                        break
            if match is None:
                notes.append(f"coupling {sup_name}<-{sub_name}: FAILED")
                unsolved.extend([sup_name, sub_name])
                final.pop(sup_name, None)
                final.pop(sub_name, None)
                continue
            M, S = match
            chosen[sup_name], chosen[sub_name] = sorted(M), sorted(S)
            coupled_shared[(sup_name, sub_name)] = sorted(S)
        elif sub_miss:
            # super exact: forced coupling adds the sub values to super too;
            # allowed only if super stays exact (neutrality check)
            sup_ds = dsmap[sup_name]
            ok = None
            for fam_s in final[sub_name]:
                S = _family_values(fam_s)
                if sup_ds.exact(removed_eids=set(removed), adds=S):
                    ok = S
                    break
            if ok is None:
                notes.append(f"coupling {sup_name}<-{sub_name}: sub witness "
                             f"not neutral in exact super; unsolved")
                unsolved.append(sub_name)
                final.pop(sub_name, None)
                continue
            chosen[sub_name] = sorted(ok)
            chosen[sup_name] = sorted(ok)   # forced neutral participation
            coupled_shared[(sup_name, sub_name)] = sorted(ok)
            notes.append(f"coupling {sup_name}<-{sub_name}: neutral forced "
                         f"participation in exact {sup_name}")
        else:
            # super miss, sub exact: put super's adds on entities that do not
            # participate in sub (participation freedom); handled at
            # assignment via class choice
            pass

    for name, fams in final.items():
        if name not in chosen:
            chosen[name] = _family_values(fams[0])

    # ---- phase 5: slot assignment (shared entities, class-gated) ---------
    entities = []   # {"id","cls","values":{name: Fraction}}

    def new_entity(cls):
        e = {"id": f"ADD-{cls}-{len(entities)+1:02d}", "cls": cls,
             "values": {}}
        entities.append(e)
        return e

    def coupling_ok(e, name, v):
        for sup_name, sub_name in couplings:
            if name == sub_name and sup_name in e["values"] and \
                    e["values"][sup_name] != v:
                return False
            if name == sup_name and sub_name in e["values"] and \
                    e["values"][sub_name] != v:
                return False
        return True

    def assign(name, vals, allowed_cls, paired=None):
        """Assign sorted vals of dataset `name` to entities of allowed
        classes; paired = (other_name) to co-assign equal values (coupled)."""
        for v in vals:
            placed = False
            for e in entities:
                if e["cls"] not in allowed_cls:
                    continue
                if name in e["values"]:
                    continue
                if not coupling_ok(e, name, v):
                    continue
                if paired and (paired in e["values"] or
                               not coupling_ok(e, paired, v)):
                    continue
                e["values"][name] = v
                if paired:
                    e["values"][paired] = v
                placed = True
                break
            if not placed:
                e = new_entity(allowed_cls[0])
                e["values"][name] = v
                if paired:
                    e["values"][paired] = v

    handled = set()
    for (sup_name, sub_name), S in coupled_shared.items():
        both = [c for c in classes if participation(c, sup_name)
                and participation(c, sub_name)]
        if not both:
            raise StackerError(f"no class may participate in both "
                               f"{sup_name} and {sub_name}")
        assign(sub_name, S, both, paired=sup_name)
        handled |= {sub_name}
        M = chosen.get(sup_name, S)
        surplus = sorted((Counter(M) - Counter(S)).elements())
        if surplus:
            solo = [c for c in classes if participation(c, sup_name)
                    and not participation(c, sub_name)]
            if not solo:
                solo = both
                notes.append(f"coupling {sup_name}<-{sub_name}: surplus "
                             f"rides a shared class (participation freedom "
                             f"invoked: those entities skip {sub_name})")
            assign(sup_name, surplus, solo)
        handled |= {sup_name}

    for name in sorted(chosen):
        if name in handled:
            continue
        allowed = [c for c in classes if participation(c, name)]
        if not allowed:
            raise StackerError(f"no class may participate in {name}")
        # super-miss/sub-exact coupling: prefer classes outside the sub
        for sup_name, sub_name in couplings:
            if name == sup_name and status.get(sub_name) == "exact":
                pref = [c for c in allowed if not participation(c, sub_name)]
                if pref:
                    allowed = pref
                else:
                    notes.append(f"{name}: adds may force participation in "
                                 f"exact {sub_name} (freedom invoked)")
        assign(name, chosen[name], allowed)

    # ---- phase 6: final battery + hard no-regression gate ----------------
    battery = []
    n_exact_after = 0
    regressions = []
    removed_set = set(removed)
    for d in datasets:
        adds = [e["values"][d.name] for e in entities if d.name in e["values"]]
        ok = d.exact(removed_eids=removed_set, adds=adds)
        n_exact_after += ok
        row = {"name": d.name, "before_exact": status[d.name] == "exact",
               "after_exact": ok,
               "adds": [str(v) for v in sorted(adds)],
               "n_removed_here": sum(1 for e, _ in d.rows
                                     if e in removed_set)}
        battery.append(row)
        if row["before_exact"] and not row["after_exact"]:
            regressions.append(d.name)
    if regressions:
        raise StackerError(f"no-regression gate violated: {regressions}")

    return {"entities": [{"id": e["id"], "cls": e["cls"],
                          "values": {k: str(v)
                                     for k, v in sorted(e["values"].items())}}
                         for e in entities],
            "removed_entities": removed,
            "battery": battery,
            "exact_before": sum(1 for d in datasets
                                if status[d.name] == "exact"),
            "exact_after": n_exact_after,
            "evaluable": len(datasets),
            "unsolved": sorted(set(unsolved)),
            "notes": notes,
            "regressions": regressions}
