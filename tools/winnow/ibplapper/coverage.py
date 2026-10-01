#!/usr/bin/env python3
"""coverage.py — coverage-counted fresh-identity gate + planted-fault trials.

The standing objection: a Freivalds-style residual gate only certifies what the
sampled identities actually TOUCH. Required fix: count, for every banked
interface row (pivot), how many independent fresh identities reach it, report
the coverage profile per run, and demand >= c_min on every banked row before
the word "certified" is used; back the argument with planted-fault mutation
trials on scratch copies (memory rule: mutation tests NEVER in place).

Touch semantics: identity I touches pivot u iff u is reachable from I's
directly-contained pivots through the banked (per-stratum, pre-closure) row
dependency graph  w -> {pivots referenced by w's banked row}.  This is the
correct detectability relation for faults planted in BANKED rows when fresh
identities are reduced against the ASSEMBLED (closed) substitutions: a fault
in banked row u contaminates exactly the closed rows of pivots w with
u in closure(w).
"""
from .engine import row_axpy   # sole lift edit: package-relative import


def reduce_against(row, subs, p):
    """Reduce dict-row against CLOSED subs; returns (residual, direct_hits).
    subs rows must reference only non-pivot (master) columns."""
    row = dict(row)
    hits = set()
    for c in [c for c in row if c in subs]:
        f = row.pop(c)
        hits.add(c)
        row_axpy(row, subs[c], (p - f) % p, p)
    return {c: v for c, v in row.items() if v}, hits


def banked_dependency_closure(banked_rows):
    """banked_rows: dict pivot -> dict(col->val) (per-stratum, pre-closure).
    Returns dict pivot -> frozenset of banked pivots in its closure (self incl)."""
    dep = {u: [c for c in r if c in banked_rows] for u, r in banked_rows.items()}
    closure = {}

    def visit(u, stack):
        if u in closure:
            return closure[u]
        if u in stack:          # cycle would be an eliminator bug
            raise AssertionError(f"cyclic banked dependency at {u}")
        stack.add(u)
        acc = {u}
        for v in dep[u]:
            acc |= visit(v, stack)
        stack.discard(u)
        closure[u] = frozenset(acc)
        return closure[u]

    for u in banked_rows:
        visit(u, set())
    return closure


def coverage_gate(fresh_rows, closed_subs, banked_rows, masters, p, c_min=2,
                  excluded_cols=frozenset()):
    """Run the fresh-identity residual gate with coverage counting.

    fresh_rows : list of dict-rows (fresh identities NOT used in the solve),
                 already filtered to rows whose every weight is in
                 closed_subs | masters (usable identities).
    closed_subs: assembled substitutions (pivot -> master-only row).
    banked_rows: pivot -> banked (pre-closure) row, i.e. what T stores.
    excluded_cols: unresolved SHELL columns (kira dump-closure gap): an
                 identity whose residual is supported ONLY there is
                 shell-blocked — valid but unverifiable; it earns NO touch
                 credit and is counted separately, never as a pass.
    Returns report dict; report["pass"] requires zero residuals AND full
    coverage at c_min.
    """
    closure = banked_dependency_closure(banked_rows)
    touch = {u: 0 for u in banked_rows}
    n_zero = n_shell_blocked = 0
    bad = []
    for k, row in enumerate(fresh_rows):
        residual, hits = reduce_against(row, closed_subs, p)
        real = {c: v for c, v in residual.items() if c not in excluded_cols}
        if real:
            bad.append({"identity": k, "residual_nnz": len(real)})
            continue
        if residual:                    # nonzero only on shell columns
            n_shell_blocked += 1
            continue
        n_zero += 1
        reached = set()
        for h in hits:
            reached |= closure.get(h, frozenset({h}))
        for u in reached:
            if u in touch:
                touch[u] += 1
    counts = sorted(touch.values())
    below = [u for u, c in touch.items() if c < c_min]
    report = {
        "n_fresh_used": len(fresh_rows),
        "n_zero_residual": n_zero,
        "n_shell_blocked": n_shell_blocked,
        "n_nonzero_residual": len(bad),
        "nonzero_residuals": bad[:10],
        "n_banked_rows": len(banked_rows),
        "coverage_min": counts[0] if counts else None,
        "coverage_median": counts[len(counts) // 2] if counts else None,
        "coverage_c_min": c_min,
        "n_rows_below_c_min": len(below),
        "rows_below_c_min_sample": below[:10],
        "coverage_ok": not below,
        "pass": (not bad) and (not below),
    }
    return report, touch


def assemble_closed(banked_rows, order, p):
    """Close banked per-stratum rows to master-only rows (assembly step).

    DEPENDENCY-DRIVEN substitution: the pivot-reference graph is
    a DAG (cross-block references point at pivots solved in earlier blocks
    or lower in-block order; pre-substitution removes the reverse edges),
    but it is NOT weight-monotone once initiate-census pseudo-masters are
    demoted (loader): their pivot rows carry tiny kira ordinals yet
    reference lower-stratum pivots that a plain ascending sweep would close
    later (measured, E5 validity probe #1). Iterative DFS closes each
    pivot's dependencies first; a cycle is an eliminator bug and asserts
    loudly. Ascending start order keeps the monotone common case
    deterministic and identical to the old sweep."""
    closed = {}
    GRAY, BLACK = 1, 2
    state = {}
    for u0 in sorted(banked_rows, key=lambda w: order.get(w, -1)):
        if state.get(u0) == BLACK:
            continue
        stack = [u0]
        while stack:
            u = stack[-1]
            st = state.get(u, 0)
            if st == BLACK:
                stack.pop()
                continue
            if st == 0:                 # first visit: expand deps
                state[u] = GRAY
                deps = [c for c in banked_rows[u]
                        if c in banked_rows and state.get(c) != BLACK]
                for c in deps:
                    assert state.get(c) != GRAY, \
                        f"cyclic pivot dependency: {u} <-> {c}"
                stack.extend(deps)
                continue
            # GRAY second visit: all deps BLACK -> close u
            row = dict(banked_rows[u])
            for c in [c for c in row if c in closed]:
                f = row.pop(c)
                row_axpy(row, closed[c], (p - f) % p, p)
            for c in row:
                assert c not in banked_rows, \
                    f"assembly not closed: pivot {u} still references " \
                    f"pivot {c}"
            closed[u] = {c: v for c, v in row.items() if v}
            state[u] = BLACK
            stack.pop()
    return closed


def planted_fault_trials(fresh_rows, banked_rows, order, masters, p,
                         n_trials=3, seed=20260706, c_min=2,
                         excluded_cols=frozenset()):
    """Mutation trials on SCRATCH COPIES of the banked rows: flip one term in
    one covered banked row, re-assemble, re-run the gate; the gate MUST detect
    (nonzero residual). Detection < 100 pct on covered rows = STRATA-K2 fires."""
    import random
    rng = random.Random(seed)
    baseline_closed = assemble_closed(banked_rows, order, p)
    base, touch = coverage_gate(fresh_rows, baseline_closed, banked_rows,
                                masters, p, c_min=c_min,
                                excluded_cols=excluded_cols)
    covered = [u for u in banked_rows if touch.get(u, 0) >= 1]
    if not covered:
        return {"baseline_gate": base, "trials": [], "detection_rate": None,
                "all_detected": False, "note": "no covered rows to mutate"}
    trials = []
    for t in range(n_trials):
        u = rng.choice(covered)
        mutated = {k: dict(v) for k, v in banked_rows.items()}  # scratch copy
        row = mutated[u]
        if row:
            c = rng.choice(list(row))
            row[c] = (row[c] + 1 + rng.randrange(p - 2)) % p or 1
        else:               # empty banked row: plant a spurious master term
            mc = rng.choice(list(masters))
            row[mc] = 1 + rng.randrange(p - 1)
        closed_m = assemble_closed(mutated, order, p)
        rep, _ = coverage_gate(fresh_rows, closed_m, mutated, masters, p,
                               c_min=c_min, excluded_cols=excluded_cols)
        trials.append({"trial": t, "mutated_pivot": u,
                       "detected": rep["n_nonzero_residual"] > 0})
    return {"baseline_gate": base, "trials": trials,
            "detection_rate": (sum(tr["detected"] for tr in trials)
                               / max(1, len(trials))),
            "all_detected": all(tr["detected"] for tr in trials)}
