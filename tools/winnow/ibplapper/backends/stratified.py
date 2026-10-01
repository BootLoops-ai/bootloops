#!/usr/bin/env python3
"""stratified.py — dense-backed mirror of engine.stratified_solve.

stratified_solve_dense mirrors engine.stratified_solve EXACTLY for the
B2F/B2FT policies (the byte-equal T-bank pair): block assignment by lead
stratum, pre-substitution of higher-stratum pivots into each block BEFORE
solving (fixpoint loop, engine.row_axpy), pass-down of leftover rows to
lower strata with the same below-current-stratum assert, per-stratum
table.append_rows with the identical key format (f"{node_key}/t{t}") and
row serialization ([(c, sorted(r.items())) ...]), B2FT = triangular global
store / B2F = closed global store with the same closure loop. The ONLY
difference: each stratum block is solved by backends.dense.eliminate_dense
when its dense buffer (n_rows * n_union_cols * 8 bytes) fits
dense_budget_bytes, else by engine.eliminate_fast — recorded per stratum
in stats as "path": "dense" | "fallback-cpu". Both block solvers are
contract-identical INCLUDING pivot-row choice and leftover emission order
(eliminate_dense mirrors eliminate_fast's (nnz, bucket-arrival) pivot rule
— required for byte-equality on rank-deficient blocks, see the MEASURED
corner in backends/dense.py), so subs / leftover / banked bytes / pass-down
row order are the same on every path mix.

transcript is NOT a parameter — dense mode has no witness recording
(in-elimination lambda tracking needs per-axpy provenance, which the
batched outer-product kernel does not carry). Certificates remain
available via RETROFIT mode on the results: witness.retrofit is post-hoc —
it lambda-solves against the caller's own rows and verifies through the
vendored receipt core, independent of which backend produced the rows.
"""
from ..engine import eliminate_fast, row_axpy
from .dense import eliminate_dense


def stratified_solve_dense(rows, order, p, stratum_of, policy="B2F",
                           forbid=(), table=None, cell_key=None,
                           node_key=None, device="cpu",
                           dense_budget_bytes=8 * 2**30):
    """Mirror of engine.stratified_solve (B2F/B2FT only; see module
    docstring). Returns (subs_all, leftover, per_stratum_stats); the stats
    dicts carry the eliminate_dense/eliminate_fast fields plus "path" and
    "dense_bytes" per stratum."""
    if policy not in ("B2F", "B2FT"):
        raise ValueError(
            f"stratified_solve_dense: policy {policy!r} not in "
            "('B2F', 'B2FT') — the dense mirror covers the fast-path pair "
            "only (B2/B3 are cost knobs of the same answer; use "
            "engine.stratified_solve for them)")
    keep_closed = policy != "B2FT"
    forbid_set = set(forbid)

    by_stratum = {}
    for r in rows:
        if not r:
            continue
        lead = max((c for c in r if c in order), key=lambda c: order[c])
        by_stratum.setdefault(stratum_of[lead], []).append(dict(r))

    subs_all, leftovers, stats_per = {}, [], {}
    all_strata = sorted(set(stratum_of.values()) | set(by_stratum),
                        reverse=True)
    for t in all_strata:
        block = by_stratum.get(t, [])
        if not block:
            continue
        # substitute previously solved (higher-stratum) pivots INTO this
        # block — identical fixpoint loop to engine.stratified_solve
        for r in block:
            hits = [c for c in r if c in subs_all]
            while hits:
                for c in hits:
                    if c in r:
                        f = r.pop(c)
                        row_axpy(r, subs_all[c], (p - f) % p, p)
                hits = [c for c in r if c in subs_all]
        block = [r for r in block if r]
        order_t = {c: order[c] for c in order if stratum_of.get(c) == t}
        union = set()
        for r in block:
            union |= set(r)
        dense_bytes = len(block) * len(union) * 8
        if dense_bytes <= dense_budget_bytes:
            subs, leftover, stt = eliminate_dense(block, order_t, p,
                                                  forbid=forbid_set,
                                                  device=device)
            stt["path"] = "dense"
        else:
            subs, leftover, stt = eliminate_fast(block, order_t, p,
                                                 forbid=forbid_set)
            stt["path"] = "fallback-cpu"
        stt["dense_bytes"] = dense_bytes
        stats_per[t] = stt
        # pass-down: leftover rows live entirely below stratum t
        n_down = 0
        for r in leftover:
            elig = [c for c in r if c in order and c not in forbid_set]
            if not elig:
                leftovers.append(r)
                continue
            lead = max(elig, key=lambda c: order[c])
            td = stratum_of[lead]
            assert td < t, "pass-down row not below current stratum"
            by_stratum.setdefault(td, []).append(r)
            n_down += 1
        stats_per[t]["rows_passed_down"] = n_down
        leftover = []
        stats_per[t]["interface_rows_out"] = sum(
            1 for c, r in subs.items()
            if any(stratum_of.get(cc, -1) < t for cc in r))
        if table is not None:
            table.append_rows(cell_key, f"{node_key}/t{t}",
                              [(c, sorted(r.items())) for c, r in subs.items()])
        # merge into global subs (close under substitution) — skipped for
        # B2FT (triangular store; single global closure via assemble_closed)
        if keep_closed:
            for pc, pr in subs_all.items():
                hits = [c for c in pr if c in subs]
                for c in hits:
                    f = pr.pop(c)
                    row_axpy(pr, subs[c], (p - f) % p, p)
        subs_all.update(subs)
        leftovers.extend(leftover)
    return subs_all, leftovers, stats_per
