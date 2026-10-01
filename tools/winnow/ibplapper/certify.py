"""certify.py — native in-elimination lambda tracking (witnesses="native").

The engine's transcript hooks (engine.eliminate_fast / stratified_solve,
`transcript=` parameter) record every linear step of a B2F/B2FT solve;
this module supplies the Transcript those hooks call and the extraction
that turns the recording into per-row multiplier certificates over the
ORIGINAL input rows — no post-hoc lambda solve (contrast witness.retrofit,
which re-derives each lambda by a Wiedemann solve after the fact).

Invariant maintained by construction: every tracked row dict r satisfies

    sum_i lam[i] * R_i  ==  r            (a live / leftover row)
    sum_i lam[i] * R_i  ==  e_u + r      (a pivot row; u = its popped pivot)

because each engine step is exactly one of reg_new / reg_copy / scale /
axpy, and each updates the lambda by the same linear operation it applies
to the row. extract_lambdas() then closes pivot lambdas over the
pivot-reference DAG (coverage.assemble_closed mirrored in lambda space),
giving per target t the WITNESS_FORMAT v1.0 identity

    sum_i lam[i] * R_i  ==  e_t - sum_m c_m e_m   (mod p),

which is verified through the vendored receipt core exactly like a
retrofit certificate — trust attaches to the receipt, never to this
recording.

Cost accounting: native tracking multiplies fill — every tracked row
carries a lambda vector over original-row indices — so the cost is
measured, never guessed. stats() reports
registered rows, lambda ops, and total/peak lambda nnz; api.eliminate
copies the numbers into ledger["witnesses"]["transcript_fill"], so every
native run ships its own measurement.

Scope: B2F/B2FT on the cpu backend only — those are the paths with
transcript hooks. B2/B3 and the dense backend refuse loudly; retrofit
remains the policy/backend-independent mode.
"""
import os

from .receipt import core


class Transcript:
    """Recorder for the engine's transcript hooks.

    Keys rows by object identity and holds a strong reference to every
    registered row dict, so an id is never reused while the transcript is
    alive. index_map (optional) translates the engine's row position into
    the caller's row index — api.eliminate uses it so lambda indices match
    the FULL System.rows list even when some rows bypass the solve.
    """

    def __init__(self, p, index_map=None):
        if p < 2:
            raise ValueError(f"p={p} is not a modulus")
        self.p = int(p)
        self.index_map = None if index_map is None else list(index_map)
        self._lam = {}      # id(row) -> {orig_row_index: coeff mod p}
        self._keep = {}     # id(row) -> row (strong ref: id stays unique)
        self.n_registered = 0
        self.lam_ops = 0
        self.lam_nnz = 0        # current total nnz over tracked lambdas
        self.peak_lam_nnz = 0

    # ------------------------------------------------------------ internals
    def _set(self, row, lam):
        k = id(row)
        old = self._lam.get(k)
        if old is not None:             # re-registration of a tracked dict
            self.lam_nnz -= len(old)
        self._lam[k] = lam
        self._keep[k] = row
        self.lam_nnz += len(lam)
        if self.lam_nnz > self.peak_lam_nnz:
            self.peak_lam_nnz = self.lam_nnz

    def lam_of(self, row):
        """Lambda vector of a tracked row (KeyError = row never registered:
        an engine path without hooks touched it — a bug, never silence)."""
        return self._lam[id(row)]

    # ------------------------------------------------------- engine hooks
    def reg_new(self, row, i):
        """Row copied from input row i (engine position; index_map applies)."""
        if self.index_map is not None:
            i = self.index_map[i]
        self._set(row, {int(i): 1})
        self.n_registered += 1

    def reg_copy(self, row, src):
        """Row copied from an already-tracked row."""
        self._set(row, dict(self._lam[id(src)]))
        self.n_registered += 1

    def scale(self, row, factor):
        """row *= factor was applied (normalize)."""
        p = self.p
        f = factor % p
        lam = self._lam[id(row)]
        for i in list(lam):
            v = lam[i] * f % p
            self.lam_ops += 1
            if v:
                lam[i] = v
            else:                       # f == 0 never happens on normalize,
                del lam[i]              # but keep the algebra total
                self.lam_nnz -= 1

    def axpy(self, dst, src, factor):
        """dst += factor*src was applied (row_axpy)."""
        p = self.p
        f = factor % p
        ld = self._lam[id(dst)]
        before = len(ld)
        for i, v in self._lam[id(src)].items():
            nv = (ld.get(i, 0) + f * v) % p
            self.lam_ops += 1
            if nv:
                ld[i] = nv
            else:
                ld.pop(i, None)
        self.lam_nnz += len(ld) - before
        if self.lam_nnz > self.peak_lam_nnz:
            self.peak_lam_nnz = self.lam_nnz

    # ---------------------------------------------------------- reporting
    def stats(self):
        """The fill measurement (per-run, ledger-bound — never a guess)."""
        return {"n_registered": self.n_registered,
                "lam_ops": self.lam_ops,
                "lam_nnz_total": self.lam_nnz,
                "peak_lam_nnz": self.peak_lam_nnz}


def extract_lambdas(transcript, pivot_rows, targets=None):
    """Close pivot lambdas over the pivot-reference DAG.

    pivot_rows: dict pivot_col -> tracked row dict (Result.subs_raw:
    triangular for B2FT, closed for B2F — both work; a closed row simply
    has no pivot references left).

    Returns dict target -> lam ({orig_row_index: coeff mod p}) satisfying
        sum_i lam[i] * R_i == e_target + closed_row(target).
    Mirrors coverage.assemble_closed in lambda space (same DFS, same
    single-pass substitution — closed rows are master-only, so one pass per
    dependency suffices); a cyclic pivot reference asserts loudly.
    """
    p = transcript.p
    want = list(pivot_rows if targets is None else targets)
    for t in want:
        if t not in pivot_rows:
            raise ValueError(f"extract_lambdas: {t} is not a solved pivot")
    lam_closed = {}
    GRAY, BLACK = 1, 2
    state = {}
    for t0 in sorted(want):
        if state.get(t0) == BLACK:
            continue
        stack = [t0]
        while stack:
            u = stack[-1]
            st = state.get(u, 0)
            if st == BLACK:
                stack.pop()
                continue
            if st == 0:                 # first visit: expand dependencies
                state[u] = GRAY
                deps = [c for c in pivot_rows[u]
                        if c in pivot_rows and state.get(c) != BLACK]
                for c in deps:
                    assert state.get(c) != GRAY, \
                        f"cyclic pivot dependency: {u} <-> {c}"
                stack.extend(deps)
                continue
            # GRAY second visit: all deps BLACK -> close u's lambda
            lam = dict(transcript.lam_of(pivot_rows[u]))
            for c, f in pivot_rows[u].items():
                if c not in pivot_rows:
                    continue
                for i, v in lam_closed[c].items():
                    nv = (lam.get(i, 0) - f * v) % p
                    if nv:
                        lam[i] = nv
                    else:
                        lam.pop(i, None)
            lam_closed[u] = lam
            state[u] = BLACK
            stack.pop()
    return {t: lam_closed[t] for t in want}


def native(system, transcript, subs_raw, closed_subs, targets, masters=None,
           write_dir=None, family=None, point=None):
    """Emit + independently verify native witnesses for `targets`.

    Same record shape and refusal semantics as witness.retrofit: default
    masters = every non-pivot column; a closed row touching a non-master
    column is refused (SHELL_RESIDUAL), never certified against it; every
    certificate is verified through the vendored receipt core against the
    system rows before it is written. Returns dict target -> record.
    """
    p = system.p
    rows = system.rows
    pivot_set = set(closed_subs)
    for t in targets:
        if t not in pivot_set:
            raise ValueError(f"witness target {t} is not a solved pivot")
    if masters is None:
        masters = set(system.all_cols) - pivot_set
    else:
        masters = set(masters)
    fpr = core.system_fingerprint(rows, p)
    if write_dir:
        os.makedirs(write_dir, exist_ok=True)

    lam_by_t = extract_lambdas(transcript, subs_raw, list(targets))
    out = {}
    for t in targets:
        rec = {"target_col": int(t), "p": p, "path": "native_transcript"}
        crow = {m: v % p for m, v in closed_subs[t].items() if v % p}
        shell_hit = sum(1 for m in crow if m not in masters)
        if shell_hit:
            rec["status"] = "SHELL_RESIDUAL"
            rec["shell_hit_cols"] = shell_hit
            out[t] = rec
            continue
        c = {int(m): int((-v) % p) for m, v in crow.items()}
        lam = {int(i): int(v % p) for i, v in lam_by_t[t].items() if v % p}
        lam_idx = sorted(lam)
        lam_val = [lam[i] for i in lam_idx]
        rec["c"] = c
        rec["lambda_nnz"] = len(lam_idx)
        rec["lam"] = lam
        rec["status"] = "CERTIFIED"
        wit = {"p": p, "target_col": int(t), "c": c, "lam_idx": lam_idx,
               "lam_val": lam_val, "n_rows": len(rows), "meta": {},
               "path": "<in-memory>"}
        ok, det = core.verify_row(rows, wit)
        rec["verify_pass"] = bool(ok)
        if not ok:                              # loud, never silent
            rec["status"] = "VERIFY_FAILED"
            rec["verify_detail"] = det
        if write_dir and ok:
            wp = os.path.join(write_dir, f"w_{p}_{t}.json")
            core.write_witness(
                wp, p=p, target_col=t, c=c, lam_idx=lam_idx,
                lam_val=lam_val, n_rows=len(rows), family=family,
                point=point,
                system={"n_rows": len(rows), "fingerprint": fpr})
            rec["witness_path"] = wp
        out[t] = rec
    return out
