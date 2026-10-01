"""api.py — the Winnow v1 library surface.

    import ibplapper as lap

    sys = lap.System(rows, p, order, forbid, stratum_of=None, meta=None)
    sched = lap.Schedule(policy="B2FT", bank=lap.Bank(path), caps=...)
    res = lap.eliminate(sys, sched, witnesses="retrofit",
                        witness_targets=[cols...])

    res.subs       # pivot -> CLOSED row (master/forbidden cols only)
    res.leftover   # relations among masters — NEVER silently absorbed
    res.witnesses  # per requested target: lambda vector + c + verify status
    res.ledger     # audit ledger (JSON-serializable; JSONL sidecar form)

Column ids are opaque ints; order / stratum_of / forbid / meta are the
entire physics surface. The engine modules (engine.py, bank.py, coverage.py)
carry the elimination semantics — behavior is pinned by the byte-equal
interface-table (T) bank regressions; THIS module is shell code only.

Contract invariants:
  - B2/B3/B2F/B2FT are answer-identical (canonical RREF equality);
    policy is a COST knob only.
  - Masters are never pivots; rank deficiency surfaces as `leftover`.
  - Any cap/censor is a loud CensoredError with a partial ledger, never
    rc=0 (print-loud-rc0 rule).
  - Witnesses verify via the vendored RECEIPT core against the CALLER's
    independent parse of the rows (the confident-wrong-table lesson).

Scope pins: single (prime, point) per call. witnesses = "retrofit"
(post-hoc Wiedemann lambda solve; policy- and backend-independent) or
"native" (in-elimination lambda tracking via the engine transcript hooks —
B2F/B2FT on the cpu backend only; its fill cost is MEASURED per run and
recorded in ledger["witnesses"]["transcript_fill"], see certify.py).
"""
import time

from . import engine, ledger as _ledger
from .bank import InterfaceTable
from .coverage import assemble_closed

VALID_POLICIES = ("B2", "B3", "B2F", "B2FT")


class CensoredError(RuntimeError):
    """A per-stratum cap fired. Carries the PARTIAL audit ledger — the run
    is loudly incomplete, never silently truncated."""

    def __init__(self, msg, partial_ledger=None):
        super().__init__(msg)
        self.ledger = partial_ledger


class Bank(InterfaceTable):
    """The Bank protocol, formalized (design §2: fp_eliminate duck-types
    `table.append_rows`). Any object with this surface is a valid
    Schedule.bank:

        append_rows(cell_key, node_key, rows, family=None) -> key
            rows = [(lead_col, [(col, val), ...]), ...]; append-only —
            re-appending an existing key MUST be refused.
        get(cell_key, node_key, family=None, verify=True) -> rows
            read-back MUST be content-hash verified.
        has(cell_key, node_key, family=None) -> bool

    This class IS the strata interface_table (verbatim lift): mmap
    append-only binary + sha256-verified read-back. v0 scaling note: the
    idx json is rewritten per append (fine <= 1e4 keys)."""


class System:
    """Sparse F_p identity rows + the physics hooks. Everything upstream of
    the rows (seeding, symmetries, dictionaries, kinematics evaluation)
    stays caller-side by design — see README §caller-side."""

    def __init__(self, rows, p, order, forbid=(), stratum_of=None, meta=None):
        if p < 2:
            raise ValueError(f"p={p} is not a modulus")
        self.p = int(p)
        self.rows = [dict(r) for r in rows]
        self.order = dict(order)
        self.forbid = set(forbid)
        bad = self.forbid & set(self.order)
        if bad:
            raise ValueError(
                f"order and forbid overlap on {len(bad)} columns "
                f"(sample {sorted(bad)[:5]}): a master must not carry an "
                "elimination rank")
        if stratum_of is None:
            stratum_of = {c: 0 for c in self.order}     # single stratum
        self.stratum_of = dict(stratum_of)
        missing = [c for c in self.order if c not in self.stratum_of]
        if missing:
            raise ValueError(f"stratum_of missing {len(missing)} order "
                             f"columns (sample {missing[:5]})")
        self.meta = meta
        self.all_cols = sorted({c for r in self.rows for c in r})

    def fingerprint(self):
        return _ledger.fingerprint_rows(self.rows, self.p)


class Schedule:
    """Elimination schedule: policy + optional Bank + optional caps.

    caps (all optional; enforced at STRATUM boundaries in v1 — granularity
    documented, loud CensoredError, never silent truncation):
        stratum_wall_s   max wall for any single stratum block
        total_wall_s     max wall for the whole solve
        stratum_fill     max banked nnz for any single stratum block
    bank_cell / bank_node: interface-table key parts (default "cell"/"node";
    a common convention is cell="tiny", node="p{p}_d{d}_e{e}").
    """

    def __init__(self, policy="B2FT", alpha=0, shadow=None, bank=None,
                 caps=None, bank_cell="cell", bank_node="node"):
        if policy not in VALID_POLICIES:
            raise ValueError(f"policy {policy!r} not in {VALID_POLICIES}")
        self.policy = policy
        self.alpha = alpha
        self.shadow = shadow
        self.bank = bank
        self.caps = dict(caps) if caps else None
        self.bank_cell = bank_cell
        self.bank_node = bank_node


class Result:
    def __init__(self, subs, subs_raw, leftover, stats_per, ledger,
                 witnesses=None):
        self.subs = subs            # CLOSED substitutions
        self.subs_raw = subs_raw    # as returned by the engine (triangular
        #                             for B2FT; closed otherwise)
        self.leftover = leftover    # master-relation rows (sec39 class)
        self.stats_per = stats_per  # per-stratum engine Stats
        self.ledger = ledger
        self.witnesses = witnesses or {}

    def save_ledger(self, path):
        return _ledger.save(self.ledger, path)


class _CapBank:
    """Bank-protocol proxy that enforces Schedule.caps at stratum boundaries
    (append_rows is called exactly once per stratum block by the engine).
    Passing NO caps means the real bank is handed to the engine directly —
    this proxy is never on the byte-equal path."""

    def __init__(self, inner, caps, led):
        self.inner = inner
        self.caps = caps
        self.led = led
        self.t0 = time.time()
        self.t_last = self.t0

    def append_rows(self, cell_key, node_key, rows, family=None):
        now = time.time()
        stratum_wall = now - self.t_last
        self.t_last = now
        nnz = sum(len(terms) for _lead, terms in rows)
        cap = self.caps.get("stratum_wall_s")
        if cap is not None and stratum_wall > cap:
            _ledger.record_censor(self.led, node_key, "stratum_wall_s",
                                  round(stratum_wall, 3), cap)
            raise CensoredError(
                f"stratum {node_key}: wall {stratum_wall:.3f}s > cap {cap}s",
                self.led)
        cap = self.caps.get("stratum_fill")
        if cap is not None and nnz > cap:
            _ledger.record_censor(self.led, node_key, "stratum_fill",
                                  nnz, cap)
            raise CensoredError(
                f"stratum {node_key}: banked nnz {nnz} > cap {cap}", self.led)
        cap = self.caps.get("total_wall_s")
        if cap is not None and now - self.t0 > cap:
            _ledger.record_censor(self.led, node_key, "total_wall_s",
                                  round(now - self.t0, 3), cap)
            raise CensoredError(
                f"total wall {now - self.t0:.3f}s > cap {cap}s", self.led)
        if self.inner is not None:
            return self.inner.append_rows(cell_key, node_key, rows,
                                          family=family)
        return None


def eliminate(system, schedule=None, witnesses="none", witness_targets=(),
              witness_masters=None, witness_dir=None, witness_seed=20260706,
              backend="cpu", device="cpu", dense_budget=8 * 2**30):
    """Run the stratified Laporta elimination per the Schedule; optionally
    emit per-row lambda witnesses. Returns Result.

    witnesses: "none" | "retrofit" | "native".
        "retrofit": post-hoc Wiedemann lambda solve per target (witness.py)
            — any policy, any backend.
        "native": in-elimination lambda tracking via the engine transcript
            hooks (certify.py) — B2F/B2FT on the cpu backend only (the only
            paths with hooks); every certificate is still verified through
            the vendored receipt core, and the transcript's fill cost is
            measured per run into ledger["witnesses"]["transcript_fill"].
    witness_targets: pivot columns to certify.
    witness_masters: columns allowed to carry witness coefficients; default
        = every non-pivot column (pass an explicit set to exclude shell
        columns — see README honesty notes).
    backend: "cpu" (default — the verbatim engine path, byte-identical to
        pre-backend builds) | "dense" (backends.stratified dense-torch
        mirror; B2F/B2FT only; requires p < 2^31 and torch). device /
        dense_budget apply to backend="dense" only ("cpu" or "cuda"; the
        per-stratum dense-buffer budget in bytes, else engine fallback).
        Caps are honored identically on both backends (the _CapBank proxy
        wraps the bank the same way — append_rows fires once per stratum
        on either path). witnesses="retrofit" works with backend="dense"
        (retrofit is post-hoc; witness.py never sees the solve path);
        witnesses="native" would need an in-elimination transcript, which
        the dense backend does not record -> NotImplementedError.
    """
    sched = schedule or Schedule()
    if backend not in ("cpu", "dense"):
        raise ValueError(f"backend={backend!r} not in ('cpu', 'dense')")
    if backend == "dense":
        if witnesses == "native":
            raise NotImplementedError(
                "witnesses='native' requires in-elimination transcript "
                "recording; the dense backend has none "
                "(stratified_solve_dense takes no transcript parameter). "
                "Use witnesses='retrofit' — retrofit certificates are "
                "post-hoc and backend-independent.")
        if sched.policy not in ("B2F", "B2FT"):
            raise NotImplementedError(
                f"backend='dense' mirrors policies B2F/B2FT only (got "
                f"{sched.policy!r}); B2/B3 are answer-identical cost knobs "
                "— run them on backend='cpu'.")
    if witnesses == "native" and sched.policy not in ("B2F", "B2FT"):
        raise NotImplementedError(
            "witnesses='native' replays the elimination transcript, which "
            f"exists on the B2F/B2FT fast path only (policy="
            f"{sched.policy!r}); run B2F/B2FT — answer-identical, policy is "
            "a cost knob — or use witnesses='retrofit' "
            "(policy-independent).")
    sys_ = system
    p, order, stratum_of = sys_.p, sys_.order, sys_.stratum_of

    led = _ledger.new_ledger(
        p=p, n_rows=len(sys_.rows), n_cols=len(sys_.all_cols),
        policy=sched.policy, caps=sched.caps, forbid_n=len(sys_.forbid),
        strata_n=len(set(stratum_of.values())),
        rows_fingerprint=sys_.fingerprint(),
        meta_echo=(sys_.meta if isinstance(sys_.meta, (str, int)) else None))

    # rows with no eliminable-order column at all go straight to leftover
    # (relations among masters as GIVEN — the engine's lead scan requires a
    # lead; behavior of rows with leads is untouched).
    solve_rows, solve_idx, pre_leftover = [], [], []
    for i, r in enumerate(sys_.rows):
        if any(c in order for c in r):
            solve_rows.append(r)
            solve_idx.append(i)
        elif r:
            pre_leftover.append(dict(r))

    transcript = None
    if witnesses == "native":
        from .certify import Transcript
        # index_map: lambda indices always name positions in the FULL
        # System.rows list (the caller's parse), never the filtered solve
        # list — the receipt-core length check depends on it.
        transcript = Transcript(p, index_map=solve_idx)

    table_arg = sched.bank
    banked_before = set(sched.bank.idx["keys"]) if sched.bank else set()
    if sched.caps:
        table_arg = _CapBank(sched.bank, sched.caps, led)

    t0 = time.time()
    if backend == "dense":
        from .backends.stratified import stratified_solve_dense
        led["backend"] = f"dense-torch/{device}"
        subs_all, leftover, stats_per = stratified_solve_dense(
            solve_rows, order, p, stratum_of, policy=sched.policy,
            forbid=sys_.forbid, table=table_arg, cell_key=sched.bank_cell,
            node_key=sched.bank_node, device=device,
            dense_budget_bytes=dense_budget)
    else:
        subs_all, leftover, stats_per = engine.stratified_solve(
            solve_rows, order, p, stratum_of, policy=sched.policy,
            alpha=sched.alpha, shadow=sched.shadow, forbid=sys_.forbid,
            table=table_arg, cell_key=sched.bank_cell,
            node_key=sched.bank_node, transcript=transcript)
    wall_solve = time.time() - t0

    leftover = pre_leftover + [dict(r) for r in leftover]
    led["per_stratum"] = {str(k): v for k, v in stats_per.items()}
    led["wall_s"] = round(wall_solve, 4)
    _ledger.record_leftover(led, leftover)
    if sched.bank:
        for k in sorted(set(sched.bank.idx["keys"]) - banked_before):
            ent = sched.bank.idx["keys"][k]
            led["banked_blocks"].append({"key": k, "sha256": ent["sha256"],
                                         "n_rows": ent["n_rows"]})

    t0 = time.time()
    closed = assemble_closed(subs_all, order, p)
    led["wall_close_s"] = round(time.time() - t0, 4)

    wit = {}
    if witnesses == "retrofit":
        from .witness import retrofit
        wit = retrofit(sys_, closed, list(witness_targets),
                       masters=witness_masters, write_dir=witness_dir,
                       seed=witness_seed)
    elif witnesses == "native":
        from .certify import native
        wit = native(sys_, transcript, subs_all, closed,
                     list(witness_targets), masters=witness_masters,
                     write_dir=witness_dir)
    elif witnesses != "none":
        raise ValueError(f"witnesses={witnesses!r}")
    if witnesses in ("retrofit", "native"):
        led["witnesses"] = {
            "mode": witnesses,
            "n_targets": len(wit),
            "n_certified": sum(1 for w in wit.values()
                               if w["status"] == "CERTIFIED"),
            "n_verified": sum(1 for w in wit.values()
                              if w.get("verify_pass")),
        }
        if transcript is not None:
            # native fill cost: recorded per run, never estimated
            led["witnesses"]["transcript_fill"] = transcript.stats()

    return Result(closed, subs_all, leftover, stats_per, led, wit)
