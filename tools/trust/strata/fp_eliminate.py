#!/usr/bin/env python3
"""fp_eliminate.py — STRATA Route-B F_p eliminator.

Own mod-p Gaussian eliminator for kira-free stratified Laporta: kira resume
cannot load only B_k, so STRATA solves its own strata.

Two pivot policies, shipped as an A/B pair:
  B2 "laporta"   : eliminate columns strictly in canonical Laporta order
                   (highest weight first); among rows exposing the column,
                   pick the shortest row (classic complexity tie-break).
  B3 "markowitz" : fill-in-aware: among eliminable columns, pick (row, col)
                   minimizing  (nnz_row - 1)*(colcount - 1) + alpha*shadow(col),
                   shadow(col) = dots + s of the integral (degree-shadow: a
                   proxy for the coefficient-degree the pivot drags into the
                   symbolic lift; alpha=0 gives pure Markowitz).
Both policies FINISH with a canonicalization pass to full RREF in the
canonical column order, so their outputs are IDENTICAL by construction if both
are correct — the A/B comparison is about cost (fill, peak live nnz, ops),
never about the answer. Comparison hooks emit JSONL rows.

Masters are never pivoted on (forbidden set); surviving relations among
masters (rank deficiency) are reported, not silently absorbed (sec39 class).

Scale target: unit gates + small systems (<= 10^5 nnz); pure-python
dict-rows.
"""
import json
import time


# ---------------------------------------------------------------- core algebra
def row_axpy(dst, src, factor, p):
    """dst += factor*src (mod p) in-place on dict rows; returns ops count."""
    ops = 0
    for c, v in src.items():
        nv = (dst.get(c, 0) + factor * v) % p
        ops += 1
        if nv:
            dst[c] = nv
        else:
            dst.pop(c, None)
    return ops


def normalize(row, col, p):
    inv = pow(row[col], p - 2, p)
    for c in list(row):
        row[c] = row[c] * inv % p
    return row


class Stats:
    def __init__(self):
        self.fill = 0            # nnz of stored substitution rows (total)
        self.peak_live_nnz = 0
        self.peak_live_rows = 0
        self.ops = 0
        self.pivots = 0
        self.t0 = time.time()

    def live(self, rows_iter):
        nnz = sum(len(r) for r in rows_iter)
        if nnz > self.peak_live_nnz:
            self.peak_live_nnz = nnz

    def done(self):
        return {"fill": self.fill, "peak_live_nnz": self.peak_live_nnz,
                "peak_live_rows": self.peak_live_rows, "ops": self.ops,
                "pivots": self.pivots, "wall_s": round(time.time() - self.t0, 3)}


# ---------------------------------------------------------------- elimination
def eliminate(rows, order, p, policy="B2", alpha=0, shadow=None, forbid=(),
              stats=None):
    """Forward elimination of `rows` (list of dict col->val mod p).

    order: dict col -> rank; HIGHER rank = eliminated FIRST (Laporta: rank =
           kira weight). Columns not in `order` or in `forbid` are never pivots.
    Returns (subs, leftover, stats_dict):
      subs: dict pivot_col -> normalized row (pivot col removed) meaning
            pivot_col = -sum(coeff*col)  [row stored as the RHS combination]
      leftover: rows with no eliminable column left (e.g. master-only
            relations — rank-deficiency evidence, NEVER dropped silently).
    """
    st = stats or Stats()
    shadow = shadow or {}
    forbid = set(forbid)
    live = [dict(r) for r in rows if r]
    subs = {}

    def reduce_by_subs(row):
        # substitute all known pivots out of row (repeat until fixpoint).
        # subs[c] stores rest with  c + rest = 0, so f*c -> -f*rest.
        changed = True
        while changed:
            changed = False
            hits = [c for c in row if c in subs]
            for c in sorted(hits, key=lambda c: -order.get(c, -1)):
                if c not in row:
                    continue
                f = row.pop(c)
                st.ops += row_axpy(row, subs[c], (p - f) % p, p)
                changed = True
        return row

    while True:
        st.live(live)
        st.peak_live_rows = max(st.peak_live_rows, len(live))
        # candidate pivots: (row_idx, col)
        cands = []
        colcount = {}
        for i, r in enumerate(live):
            for c in r:
                if c not in forbid and c in order:
                    colcount[c] = colcount.get(c, 0) + 1
        if not colcount:
            break
        if policy == "B2":
            # strict Laporta: next column = highest order among eliminable
            col = max(colcount, key=lambda c: order[c])
            rows_with = [(len(live[i]), i) for i, r in enumerate(live) if col in r]
            _, i = min(rows_with)
            pick = (i, col)
        elif policy == "B3":
            best = None
            for i, r in enumerate(live):
                nr = len(r)
                for c in r:
                    if c in forbid or c not in order:
                        continue
                    score = (nr - 1) * (colcount[c] - 1) + alpha * shadow.get(c, 0)
                    key = (score, -order[c], nr, i)
                    if best is None or key < best[0]:
                        best = (key, (i, c))
            pick = best[1]
        else:
            raise ValueError(policy)
        i, col = pick
        row = live.pop(i)
        row = reduce_by_subs(row)
        if col not in row:      # pivot vanished under substitution; requeue
            if row:
                live.append(row)
            continue
        normalize(row, col, p)
        row.pop(col)                # row now = rest, with  col + rest = 0
        # substitute new pivot into existing subs (keep subs closed)
        for pc, pr in subs.items():
            if col in pr:
                f = pr.pop(col)
                st.ops += row_axpy(pr, row, (p - f) % p, p)
        subs[col] = row
        st.pivots += 1
        st.fill += len(row)
        # eliminate col from remaining live rows (f*col -> -f*rest)
        for r in live:
            if col in r:
                f = r.pop(col)
                st.ops += row_axpy(r, row, (p - f) % p, p)
        live = [r for r in live if r]

    leftover = []
    for r in live:
        r = reduce_by_subs(r)
        if r:
            leftover.append(r)
    return subs, leftover, st.done()


def eliminate_fast(rows, order, p, forbid=(), stats=None, transcript=None):
    """B2 FAST PATH: strict-Laporta-order elimination via lead
    bucketing — contract-identical to eliminate(policy='B2') (CLOSED subs +
    reduced leftover; unique RREF given the rank-determined pivot set), but
    O(ops) instead of a full candidate rescan per pivot.

    transcript: optional certify.Transcript (default None = OFF, zero
    behavior change). When given, every linear step (entry copy, normalize
    scale, axpy) is recorded so closed rows can be re-expressed as F_p
    combinations of the ORIGINAL rows (certify.extract_lambdas). Input rows
    must already be registered with the transcript (reg_new/reg_copy).

    Key invariant making lazy == eager: a row containing column c has lead
    >= c, so when c becomes the current (descending) column, the rows
    containing c are EXACTLY bucket[c] and all higher pivots have already
    been substituted out of them. v0's eager per-pivot sweep touches the
    same rows with the same axpys.

    Fill semantics match v0: fill is counted at pivot creation, when the row
    is (identically to v0) reduced against all existing higher pivots.
    peak_live_rows / peak_live_nnz count bucketed (not-yet-pivot) rows only,
    like v0's `live`.
    """
    import heapq

    st = stats or Stats()
    forbid = set(forbid)

    def lead_of(r):
        best, bo = None, -1
        for c in r:
            if c in forbid:
                continue
            o = order.get(c)
            if o is not None and o > bo:
                bo, best = o, c
        return best

    buckets, leftover = {}, []
    n_live = total_nnz = 0
    for r in rows:
        if not r:
            continue
        rc = dict(r)
        if transcript is not None:
            transcript.reg_copy(rc, r)
        r = rc
        lead = lead_of(r)
        if lead is None:
            leftover.append(r)
        else:
            buckets.setdefault(lead, []).append(r)
            n_live += 1
            total_nnz += len(r)
    heap, inheap = [], set()

    def push(c):
        if c not in inheap:
            heapq.heappush(heap, (-order[c], c))
            inheap.add(c)

    for c in buckets:
        push(c)
    subs = {}
    st.peak_live_rows = max(st.peak_live_rows, n_live)
    st.peak_live_nnz = max(st.peak_live_nnz, total_nnz)
    while heap:
        _, col = heapq.heappop(heap)
        inheap.discard(col)
        rows_c = buckets.pop(col, None)
        if not rows_c:
            continue
        assert col not in subs, "column re-processed — heap logic broken"
        rows_c.sort(key=len)
        piv, rest = rows_c[0], rows_c[1:]
        if transcript is not None:
            transcript.scale(piv, pow(piv[col], p - 2, p))
        normalize(piv, col, p)
        piv.pop(col)
        subs[col] = piv
        st.pivots += 1
        st.fill += len(piv)
        n_live -= 1
        total_nnz -= len(piv) + 1
        for r in rest:
            f = r.pop(col)
            before = len(r) + 1
            st.ops += row_axpy(r, piv, (p - f) % p, p)
            if transcript is not None:
                transcript.axpy(r, piv, (p - f) % p)
            total_nnz += len(r) - before
            if not r:
                n_live -= 1
                continue
            lead = lead_of(r)
            if lead is None:
                leftover.append(r)
                n_live -= 1
                total_nnz -= len(r)
            else:
                buckets.setdefault(lead, []).append(r)
                push(lead)
        st.peak_live_rows = max(st.peak_live_rows, n_live)
        st.peak_live_nnz = max(st.peak_live_nnz, total_nnz)
    # closing pass (ascending): referenced pivot rows are already closed, and
    # closed rows contain no subs columns, so one sweep per row suffices.
    for col in sorted(subs, key=lambda c: order[c]):
        row = subs[col]
        hits = sorted((c for c in row if c in subs),
                      key=lambda c: -order[c])
        for c in hits:
            f = row.pop(c)
            st.ops += row_axpy(row, subs[c], (p - f) % p, p)
            if transcript is not None:
                transcript.axpy(row, subs[c], (p - f) % p)
    # leftover rows contain no eliminable columns by construction (v0's
    # final reduce_by_subs is a no-op on them: subs keys are eliminable).
    return subs, [r for r in leftover if r], st.done()


def rref_canonical(rows, order, p):
    """Full RREF in canonical order (policy-independent canonical form):
    returns sorted list of normalized rows [(lead_col, ((col,val),...)), ...].
    Used for gate comparisons; masters allowed as leads (no forbid)."""
    subs, leftover, _ = eliminate(rows, order, p, policy="B2")
    assert not leftover or all(not r for r in leftover), \
        "rref_canonical: nonzero rows without eliminable columns"
    out = []
    for c, r in subs.items():
        terms = tuple(sorted(((cc, vv) for cc, vv in r.items())))
        out.append((c, terms))
    return sorted(out)


# ------------------------------------------------------- stratified driver v0
def stratified_solve(rows, order, p, stratum_of, policy="B2", alpha=0,
                     shadow=None, forbid=(), table=None, cell_key=None,
                     node_key=None, transcript=None):
    """Solve stratum by stratum, TOP-DOWN in stratum index (deepest sectors
    first — their pivots are the highest-weight columns), banking each
    stratum's interface rows (pivot rows referencing lower strata) into the
    interface table T if given. Returns (subs_all, leftover, per_stratum_stats).

    stratum_of: dict col -> stratum id (int; HIGHER = eliminated earlier,
    e.g. t of the home sector). Equations are assigned to the stratum of
    their highest-order column. Within a stratum block, ONLY in-stratum
    columns are pivots; rows that reduce to lower-stratum-only content are
    passed DOWN to that stratum's block (never eliminated across strata —
    that is the whole point of the interface).

    transcript: optional certify.Transcript (default None = OFF, zero
    behavior change). Records row origins (index into `rows`) and every
    linear step so per-row multiplier certificates over the ORIGINAL rows
    can be extracted (certify.extract_lambdas). B2F/B2FT paths only.
    """
    # B2FT: like B2F but the GLOBAL subs
    # store stays TRIANGULAR — the per-stratum "merge into global subs"
    # closure is skipped (block pre-substitution already chases to fixpoint,
    # so banked per-stratum rows are IDENTICAL; the single global closure
    # happens once in assemble_closed, which downstream runs anyway).
    keep_closed = policy != "B2FT"
    if policy == "B2FT":
        policy = "B2F"
    if transcript is not None and policy != "B2F":
        raise ValueError("transcript recording supported only for B2F/B2FT")
    by_stratum = {}
    for i, r in enumerate(rows):
        if not r:
            continue
        lead = max((c for c in r if c in order), key=lambda c: order[c])
        rc = dict(r)
        if transcript is not None:
            transcript.reg_new(rc, i)
        by_stratum.setdefault(stratum_of[lead], []).append(rc)

    subs_all, leftovers, stats_per = {}, [], {}
    all_strata = sorted(set(stratum_of.values()) | set(by_stratum), reverse=True)
    for t in all_strata:
        block = by_stratum.get(t, [])
        if not block:
            continue
        # substitute previously solved (higher-stratum) pivots INTO this block
        for r in block:
            hits = [c for c in r if c in subs_all]
            while hits:
                for c in hits:
                    if c in r:
                        f = r.pop(c)
                        row_axpy(r, subs_all[c], (p - f) % p, p)
                        if transcript is not None:
                            transcript.axpy(r, subs_all[c], (p - f) % p)
                hits = [c for c in r if c in subs_all]
        block = [r for r in block if r]
        order_t = {c: order[c] for c in order if stratum_of.get(c) == t}
        if policy == "B2F":
            subs, leftover, stt = eliminate_fast(block, order_t, p,
                                                 forbid=forbid,
                                                 transcript=transcript)
        else:
            subs, leftover, stt = eliminate(block, order_t, p, policy=policy,
                                            alpha=alpha, shadow=shadow,
                                            forbid=forbid)
        stats_per[t] = stt
        # pass-down: leftover rows live entirely below stratum t
        n_down = 0
        for r in leftover:
            elig = [c for c in r if c in order and c not in forbid]
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
                    if transcript is not None:
                        transcript.axpy(pr, subs[c], (p - f) % p)
        subs_all.update(subs)
        leftovers.extend(leftover)
    return subs_all, leftovers, stats_per


# ------------------------------------------------------------- A/B table hook
def b2b3_compare(rows, order, p, shadow, forbid, tag, out_jsonl=None,
                 alphas=(0, 1, 4)):
    """Run B2 and B3 (at each alpha) on copies of the same system; emit the
    comparison table rows. Correctness cross-check: equal rank (pivot count)
    and IDENTICAL canonical RREF of the reconstructed substitution systems
    (pivot SETS may legitimately differ between policies when free columns
    exist; the canonical form may not)."""
    def canon(subs, leftover):
        recon = []
        for pc, pr in subs.items():
            row = {pc: 1}
            row.update(pr)          # NOT dict(**pr): int keys break keyword form
            recon.append(row)
        recon.extend(dict(r) for r in leftover)
        return rref_canonical(recon, order, p)

    out = []
    subs2, left2, s2 = eliminate([dict(r) for r in rows], order, p, "B2",
                                 forbid=forbid)
    canon2 = canon(subs2, left2)
    out.append({"tag": tag, "policy": "B2", **s2, "n_subs": len(subs2),
                "n_leftover": len(left2)})
    for a in alphas:
        subs3, left3, s3 = eliminate([dict(r) for r in rows], order, p, "B3",
                                     alpha=a, shadow=shadow, forbid=forbid)
        assert len(subs3) == len(subs2), \
            f"rank mismatch B2 vs B3(alpha={a}) — eliminator bug"
        assert canon(subs3, left3) == canon2, \
            f"canonical-form mismatch B2 vs B3(alpha={a}) — eliminator bug"
        out.append({"tag": tag, "policy": f"B3_a{a}", **s3,
                    "n_subs": len(subs3), "n_leftover": len(left3)})
    if out_jsonl:
        with open(out_jsonl, "a") as fh:
            for r in out:
                fh.write(json.dumps(r) + "\n")
    return out
