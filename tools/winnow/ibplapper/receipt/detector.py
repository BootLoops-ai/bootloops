#!/usr/bin/env python3
"""detector.py — structural false-master / wrong-basis detector.

The core lesson (measured on an archived audit): "incompleteness is not
self-signaling" — an insufficient system produces COMPLETE-LOOKING reduction
rows onto a too-large 'master' set (+20 false masters in one audited system),
with every internal consistency check passing. The fix is structural, not
statistical: eliminate the system with the DECLARED masters forbidden and
audit what survives.

Alarms:
  A1 master_relations   leftover rows supported only on declared masters =
                        the declared basis is rank-deficient / dependent
                        (also fires when a reducible column is
                        declared master — its defining row descends to a
                        master-only relation).
  A2 pseudo_masters     columns surviving in closed rows that are NOT
                        declared masters and NOT declared shell = false
                        masters (the solver would have
                        silently called them masters).
  A3 uncovered_targets  requested target columns with no closing pivot
                        (no covering identity
                        chain in the staged system).
  A4 foreign_support    claimed table rows carrying coefficients on columns
                        outside the declared master basis.

detect() is a DIAGNOSTIC: it uses elimination internally but assumes nothing
about how the audited TABLE was produced. The trust anchor remains the
per-row lambda witness (core.py).

Eliminator: self-contained port of the strata B2F bucketed strict-order
eliminator (fp_eliminate.eliminate_fast), held to a byte-equal output contract against
the reference eliminator.
"""
import heapq


def _row_axpy(dst, src, factor, p):
    for c, v in src.items():
        nv = (dst.get(c, 0) + factor * v) % p
        if nv:
            dst[c] = nv
        else:
            dst.pop(c, None)


def eliminate(rows, order, p, forbid=()):
    """Strict-order forward elimination (B2F port). rows: list dict col->val.
    order: dict col->rank (higher rank eliminated first); cols not in order
    or in forbid are never pivots. Returns (subs, leftover):
      subs: pivot_col -> CLOSED row (pivot removed; col + row = 0)
      leftover: fully-reduced rows with no eliminable column (never dropped).
    """
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
    for r in rows:
        if not r:
            continue
        r = dict(r)
        lead = lead_of(r)
        if lead is None:
            leftover.append(r)
        else:
            buckets.setdefault(lead, []).append(r)
    heap, inheap = [], set()

    def push(c):
        if c not in inheap:
            heapq.heappush(heap, (-order[c], c))
            inheap.add(c)

    for c in buckets:
        push(c)
    subs = {}
    while heap:
        _, col = heapq.heappop(heap)
        inheap.discard(col)
        rows_c = buckets.pop(col, None)
        if not rows_c:
            continue
        rows_c.sort(key=len)
        piv, rest = rows_c[0], rows_c[1:]
        inv = pow(piv[col], p - 2, p)
        for c in list(piv):
            piv[c] = piv[c] * inv % p
        piv.pop(col)
        subs[col] = piv
        for r in rest:
            f = r.pop(col)
            _row_axpy(r, piv, (p - f) % p, p)
            if not r:
                continue
            lead = lead_of(r)
            if lead is None:
                leftover.append(r)
            else:
                buckets.setdefault(lead, []).append(r)
                push(lead)
    # closing pass (ascending): make every subs row master/survivor-only
    for col in sorted(subs, key=lambda c: order[c]):
        row = subs[col]
        hits = sorted((c for c in row if c in subs),
                      key=lambda c: -order[c])
        for c in hits:
            f = row.pop(c)
            _row_axpy(row, subs[c], (p - f) % p, p)
    return subs, [r for r in leftover if r]


def detect(rows, p, masters, order=None, shell=(), targets=(), table=None):
    """Run the structural audit. Returns (alarm: bool, report dict).

    rows:    list of dict col->val (mod p).
    masters: declared master basis (iterable of cols).
    order:   dict col->rank; default = col id itself for all non-master cols.
    shell:   declared out-of-closure columns (dump-truncation class), never
             alarmed on.
    targets: columns whose closure is demanded (A3).
    table:   optional claimed table {target_col: {col: val}} for A4.
    """
    masters = set(masters)
    shell = set(shell)
    if order is None:
        allcols = {c for r in rows for c in r}
        order = {c: c for c in allcols if c not in masters}
    subs, leftover = eliminate(rows, order, p, forbid=masters | shell)

    # A1: relations among declared masters (rank deficiency / fake master)
    a1 = [sorted(r) for r in leftover]

    # A2: survivors in closed rows outside declared masters+shell
    survivors = set()
    for row in subs.values():
        for c in row:
            if c not in masters and c not in shell:
                survivors.add(c)
    a2 = sorted(survivors)

    # A3: demanded targets that did not close
    a3 = sorted(t for t in targets if t not in subs)

    # A4: claimed rows with support outside the declared basis
    a4 = {}
    if table:
        for t, crow in table.items():
            bad = sorted(int(m) for m in crow if int(m) not in masters)
            if bad:
                a4[int(t)] = bad

    alarm = bool(a1 or a2 or a3 or a4)
    report = {
        "n_rows": len(rows), "n_pivots": len(subs),
        "n_declared_masters": len(masters), "n_shell": len(shell),
        "alarm": alarm,
        "A1_master_relations": {"fired": bool(a1), "n": len(a1),
                                "sample": a1[:3]},
        "A2_pseudo_masters": {"fired": bool(a2), "n": len(a2),
                              "sample": a2[:10]},
        "A3_uncovered_targets": {"fired": bool(a3), "n": len(a3),
                                 "sample": a3[:10]},
        "A4_foreign_support": {"fired": bool(a4), "n": len(a4),
                               "sample": dict(list(a4.items())[:3])},
    }
    return alarm, report
