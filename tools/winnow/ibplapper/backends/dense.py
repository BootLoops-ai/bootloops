#!/usr/bin/env python3
"""dense.py — dense int64 torch eliminator (the GPU/dense backend kernel).

eliminate_dense(rows, order, p, forbid=(), device='cpu', stats=None) is
CONTRACT-IDENTICAL to engine.eliminate_fast on a stratum block: same
`subs` dict (keys, values, INSERTION ORDER) and the same `leftover` rows
(values AND emission order). Only subs/leftover are contract; the stats
fields ops/peak_live_* are BACKEND-DEFINED cost measures (a dense
outer-product axpy touches whole rows, so its op count is not comparable
to the dict path's per-nonzero count). fill/pivots are honest: pivots =
number of pivot rows (rank of the eligible-column projection,
backend-independent), fill = nnz of the CLOSED substitution rows (the dict
path counts fill at pivot creation, pre-closure; documented, not contract).

Pivot choice: the license and its measured limit
------------------------------------------------
engine.py's b2b3_compare docstring/asserts license answer-identity across
pivot POLICIES: "Both policies FINISH with a canonicalization pass to full
RREF in the canonical column order, so their outputs are IDENTICAL by
construction" — the pivot-column SET is rank-determined by the column
ordering (descending `order` rank, ties ascending column id — the same
(-order[c], c) key as eliminate_fast's heap), and the canonical closed
form is choice-invariant. That license makes the pivot ROW choice free for
the CANONICAL form, and byte-sufficient on FULL-ROW-RANK blocks (no
leftover): there the closed subs rows are the unique RREF basis.

MEASURED corner (synthetic 3-strata fixture): when a block is
rank-deficient (leftover rows exist), the closed pivot rows are unique
only MODULO span(leftover) — adding a leftover row (zero on every eligible
column) to a pivot row preserves 1-at-pivot/0-at-other-pivots. Different
pivot ROW choices then land on different (all individually correct)
representatives, and the banked bytes differ. So this backend does NOT
exercise the free choice: it MIRRORS eliminate_fast's pivot rule exactly —
among the active rows holding the column, take the minimal (nnz,
bucket-arrival) row, where bucket-arrival replays the engine's bucket
append order (input order first, then displaced rows in the sorted order
the engine re-buckets them) — and emits leftover rows in the engine's
order (eligible-free-at-entry rows in input order, then rows in death
order). With the same pivot rows, forward values, closure fixpoint, and
leftover values are identical row-for-row, so subs AND leftover are
byte-equal on every input, rank-deficient or not.

Semantics (mirrors eliminate_fast step for step)
------------------------------------------------
  (a) eligible columns = in `order` and not in `forbid`;
  (b) the block's column union is sorted eligible-first by DESCENDING
      order rank (ties: ascending column id), non-eligible columns after
      (they are never pivots; their internal order is cosmetic);
  (c) forward pass, columns left->right (descending rank): a column with
      no nonzero among active rows is skipped (not a pivot); otherwise the
      minimal-(nnz, arrival) active row becomes the pivot (the engine
      mirror above), is normalized to pivot coefficient 1 (mod-p inverse
      via pow(v, p-2, p)), and the column is eliminated from all other
      ACTIVE rows by one batched outer-product axpy mod p;
  (d) closure/backward pass to the full CLOSED form: pivots in ascending
      rank order, each pivot column eliminated from all other pivot rows
      (full RREF; the row used is already closed — back-substitution);
  (e) subs insertion order = descending pivot-column rank, exactly
      eliminate_fast's heap pop order (the T-bank append iterates
      subs.items(), so this order is part of the byte-equal contract);
  (f) leftover = rows with zero coefficients on ALL eligible columns after
      reduction; returned zero-dropped, supported on master/non-order
      columns only (an active row is zeroed at each eligible column the
      moment that column is processed and never repopulated there, because
      pivot rows carry no support on already-processed eligible columns).

int64 safety
------------
Requires p < 2^31 (raises ValueError otherwise — the typed refusal): all
residues are in [0, p), so a*b <= (2^31-1)^2 < 2^62 and |A - f*pivrow| <
2^62 + 2^31 < 2^63; torch.remainder(..., p) is applied after every
product-add, so no intermediate ever leaves the int64 range.

Input assumption (same as the dict engine's canonical rows): row values
are nonzero mod p. (The dict path has an unreachable-in-practice quirk for
an explicit 0-valued lead entry; dense treats a 0 residue as absent.)

device='cpu' is fully supported and tested; device='cuda' runs the same
code path (GPU walls are collected on a GPU host, not here).
"""
import torch

from ..engine import Stats

_P_MAX = 2**31


def _axpy_block(A, sel, jpos, pivrow, p):
    """A[sel] -= A[sel, jpos] (outer) pivrow, mod p — the single mod-
    arithmetic kernel of the backend (and the mutation-test hook point).
    Returns the number of cell updates (backend-defined ops measure)."""
    f = A[sel, jpos]
    A[sel] = torch.remainder(A[sel] - f.unsqueeze(1) * pivrow.unsqueeze(0), p)
    return int(sel.numel()) * int(pivrow.numel())


def eliminate_dense(rows, order, p, forbid=(), device="cpu", stats=None):
    """Dense torch eliminator; contract-identical to engine.eliminate_fast
    on a stratum block. Returns (subs, leftover, stats_dict); stats_dict is
    the engine Stats().done() shape plus "backend": "dense-torch/<device>".
    See module docstring for the full semantics + the pivot-mirror rule."""
    p = int(p)
    if not (2 <= p < _P_MAX):
        raise ValueError(
            f"eliminate_dense: p={p} outside [2, 2^31) — int64 safety needs "
            "a*b < 2^62 and a*b + c < 2^63 (module docstring); use the "
            "dict engine for larger moduli")
    st = stats or Stats()
    forbid = set(forbid)
    live = [r for r in rows if r]

    # (b) column union, eligible-first by descending rank
    elig_set, other_set = set(), set()
    for r in live:
        for c in r:
            (elig_set if (c in order and c not in forbid)
             else other_set).add(c)
    elig = sorted(elig_set, key=lambda c: (-order[c], c))
    others = sorted(other_set)
    col_list = elig + others
    pos = {c: k for k, c in enumerate(col_list)}
    n, m = len(live), len(col_list)
    n_elig = len(elig)
    if n == 0:
        sd = st.done()
        sd["backend"] = f"dense-torch/{device}"
        return {}, [], sd

    ii, jj, vv = [], [], []
    for i, r in enumerate(live):
        for c, v in r.items():
            ii.append(i)
            jj.append(pos[c])
            vv.append(v % p)
    A = torch.zeros((n, m), dtype=torch.int64, device=device)
    A[torch.tensor(ii, device=device), torch.tensor(jj, device=device)] = \
        torch.tensor(vv, dtype=torch.int64, device=device)
    st.peak_live_rows = max(st.peak_live_rows, n)
    st.peak_live_nnz = max(st.peak_live_nnz, len(vv))  # input nnz (backend-defined)

    # bucket-arrival mirror: input order first, displaced rows renumbered
    # in the engine's re-bucket (sorted) order as columns are processed
    arrival = list(range(n))
    counter = n
    leftover_seq = []                 # engine leftover emission order
    active = torch.ones(n, dtype=torch.bool, device=device)
    nz_tot0 = (A != 0).sum(1)
    nz_el0 = (A[:, :n_elig] != 0).sum(1) if n_elig else nz_tot0 * 0
    for i in range(n):
        if int(nz_tot0[i]) == 0:
            active[i] = False         # all-zero residue row: dropped
        elif int(nz_el0[i]) == 0:
            leftover_seq.append(i)    # eligible-free at entry (input order)

    # (c) forward pass
    pivots = []                       # (col, row_idx, col_pos), descending rank
    for jpos, c in enumerate(elig):
        cand = torch.nonzero((A[:, jpos] != 0) & active).flatten()
        if cand.numel() == 0:
            continue                  # skipped column: not a pivot
        # engine pivot rule: minimal (nnz, bucket-arrival) row
        nnz_c = (A[cand] != 0).sum(1).tolist()
        cand_l = cand.tolist()
        sorted_k = sorted(range(len(cand_l)),
                          key=lambda k: (nnz_c[k], arrival[cand_l[k]]))
        r = cand_l[sorted_k[0]]
        inv = pow(int(A[r, jpos]), p - 2, p)
        A[r] = torch.remainder(A[r] * inv, p)
        rest_l = [cand_l[k] for k in sorted_k[1:]]
        if rest_l:
            rest = torch.tensor(rest_l, device=device)
            st.ops += _axpy_block(A, rest, jpos, A[r], p)
            # re-bucket mirror: new arrivals + death detection, in the
            # engine's (sorted) iteration order
            nz_el = (A[rest, :n_elig] != 0).sum(1).tolist()
            nz_tot = (A[rest] != 0).sum(1).tolist()
            for k, ri in enumerate(rest_l):
                arrival[ri] = counter
                counter += 1
                if nz_tot[k] == 0:
                    active[ri] = False          # vanished: dropped
                elif nz_el[k] == 0:
                    leftover_seq.append(ri)     # death order
        active[r] = False
        pivots.append((c, r, jpos))

    # (d) closure/backward pass (ascending rank; used rows already closed)
    if pivots:
        piv_rows_t = torch.tensor([r for _, r, _ in pivots], device=device)
        for k in range(len(pivots) - 1, 0, -1):
            _c, r, jpos = pivots[k]
            earlier = piv_rows_t[:k]  # only higher-rank pivot rows can hold jpos
            nz = torch.nonzero(A[earlier, jpos] != 0).flatten()
            if nz.numel() == 0:
                continue
            st.ops += _axpy_block(A, earlier[nz], jpos, A[r], p)

    # (e) subs, insertion order = descending pivot-column rank
    subs, fill = {}, 0
    for c, r, jpos in pivots:
        rowv = A[r]
        assert int(rowv[jpos]) == 1, "pivot not normalized — backend bug"
        d = {}
        for kk in torch.nonzero(rowv).flatten().tolist():
            if kk != jpos:
                d[col_list[kk]] = int(rowv[kk])
        subs[c] = d
        fill += len(d)
    st.pivots = len(pivots)
    st.fill = fill

    # (f) leftover: engine emission order, zero-dropped
    leftover = []
    for i in leftover_seq:
        nzp = torch.nonzero(A[i]).flatten().tolist()
        if not nzp:
            continue
        assert all(kk >= n_elig for kk in nzp), \
            "leftover row retains eligible support — backend bug"
        leftover.append({col_list[kk]: int(A[i][kk]) for kk in nzp})
    # accounting: at the end active rows are exactly the leftover-emitted
    # ones (pivot rows and vanished rows were deactivated)
    assert int(active.sum()) == len(leftover_seq), \
        "row accounting broken — some live row neither pivot nor leftover"

    sd = st.done()
    sd["backend"] = f"dense-torch/{device}"
    return subs, leftover, sd
