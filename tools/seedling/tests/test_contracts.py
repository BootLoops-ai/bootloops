"""Regression tests for the hardened contracts:

- receipt-routing: backend chosen by stratum-size-vs-budget; dense alloc
  failure on a receipt pass falls back to the engine path — NEVER OOM-dies.
- gate bookkeeping: every gate emission names its receipt/evidence file;
  SKIPPED gates name a reason.
- degenerate-input guards on fitting/projection code (margins prediction,
  bulkmodel combinatorics).
Synthetic fixtures except the two frozen-model legs, which need the
registered model + labels (SEEDLING_ML_DIR + SEEDLING_LABELS) and skip
without them."""
import json
import math
import os
import sys
import tempfile

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", ".."))

from seedling import bulkmodel, cli, receipts  # noqa: E402


# ---------------------------------------------------------------- fixtures
class FakeSystem:
    """Duck-typed ibplapper System: .rows (list of dict col->coeff),
    .stratum_of (dict col->stratum)."""

    def __init__(self, rows, stratum_of):
        self.rows = rows
        self.stratum_of = stratum_of


def _fake_system(n_rows, n_cols, stratum=0):
    cols = [f"c{stratum}_{j}" for j in range(n_cols)]
    stratum_of = {c: stratum for c in cols}
    rows = [{cols[i % n_cols]: 1, cols[(i + 1) % n_cols]: 2}
            for i in range(n_rows)]
    return FakeSystem(rows, stratum_of)


class FakeResult:
    def __init__(self, tg=("x", "y"), certified=True):
        self.subs = {t: {} for t in tg}
        self.witnesses = {t: {"status": "CERTIFIED" if certified else "FAILED",
                              "verify_pass": certified} for t in tg}


# ------------------------------------------------------------- D1: routing
def test_backend_choice_small_system_dense():
    sysm = _fake_system(10, 5)
    backend, info = receipts.choose_backend(sysm, budget_bytes=8 * 2**30)
    assert backend == "dense"
    assert info["max_stratum_bytes_est"] > 0


def test_backend_choice_big_stratum_routes_cpu():
    # one stratum whose dense estimate exceeds a tiny budget -> cpu
    sysm = _fake_system(1000, 100)
    backend, info = receipts.choose_backend(sysm, budget_bytes=1000)
    assert backend == "cpu"
    assert info["max_stratum_bytes_est"] > info["effective_cap_bytes"]


def test_backend_choice_respects_va_headroom(monkeypatch):
    # nominal budget huge, but measured VA headroom tiny -> cpu
    # (the L3 failure mode: 8G budget "fine", allocator refused 6.98G)
    sysm = _fake_system(1000, 100)
    monkeypatch.setattr(receipts, "va_headroom_bytes", lambda: 2000)
    backend, info = receipts.choose_backend(sysm, budget_bytes=8 * 2**30)
    assert backend == "cpu"
    assert info["effective_cap_bytes"] == 1000  # safety 0.5 * 2000


def test_emit_receipts_dense_oom_falls_back_to_cpu(tmp_path):
    sysm = _fake_system(10, 5)
    calls = []

    def fake_eliminate(system, sched, backend, device, dense_budget, **kw):
        calls.append(backend)
        if backend == "dense":
            raise RuntimeError(
                "DefaultCPUAllocator: can't allocate memory: you tried "
                "to allocate 6981091584 bytes")
        return FakeResult()

    rec = receipts.emit_receipts(sysm, str(tmp_path / "w"),
                                 eliminate=fake_eliminate,
                                 schedule=lambda policy: policy)
    assert rec["ok"] is True
    assert rec["backend_used"] == "cpu-fallback"
    assert rec["fallbacks"] and rec["fallbacks"][0]["from"] == "dense"
    # first call dense (routed), second cpu (fallback), pass2 stays cpu
    assert calls[0] == "dense" and set(calls[1:]) == {"cpu"}
    # D2/D4 bookkeeping: record names its evidence + measured pass walls
    assert rec["witness_dir"] == str(tmp_path / "w")
    assert "pass1_wall_s" in rec and "pass2_wall_s" in rec


def test_emit_receipts_cpu_failure_is_fail_closed(tmp_path):
    sysm = _fake_system(10, 5)

    def fake_eliminate(system, sched, backend, device, dense_budget, **kw):
        raise RuntimeError("DefaultCPUAllocator: can't allocate memory")

    monkey_headroom = receipts.va_headroom_bytes
    try:
        receipts.va_headroom_bytes = lambda: 0     # force cpu routing
        with pytest.raises(RuntimeError):
            receipts.emit_receipts(sysm, str(tmp_path / "w"),
                                   eliminate=fake_eliminate,
                                   schedule=lambda policy: policy)
    finally:
        receipts.va_headroom_bytes = monkey_headroom


def test_emit_receipts_non_alloc_error_raises(tmp_path):
    sysm = _fake_system(10, 5)

    def fake_eliminate(system, sched, backend, device, dense_budget, **kw):
        raise ValueError("some genuine bug")        # not an alloc failure

    with pytest.raises(ValueError):
        receipts.emit_receipts(sysm, str(tmp_path / "w"),
                               eliminate=fake_eliminate,
                               schedule=lambda policy: policy)


# -------------------------------------------------- D4: gate bookkeeping
def test_gate_record_without_evidence_raises(tmp_path):
    led = str(tmp_path / "ledger.jsonl")
    with pytest.raises(ValueError, match="names no evidence"):
        cli._ledger_write(led, {"phase": "two_slice_gate",
                                "verdict": "PASS", "fails1": []})
    assert not os.path.exists(led)                  # nothing half-written


def test_gate_record_with_evidence_passes(tmp_path):
    led = str(tmp_path / "ledger.jsonl")
    cli._ledger_write(led, {"phase": "two_slice_gate", "verdict": "PASS",
                            "evidence": {"fresh": "f.json", "ref": "r.json"}})
    cli._ledger_write(led, {"phase": "receipts", "ok": True,
                            "witness_dir": "/x/receipts"})
    rows = [json.loads(l) for l in open(led)]
    assert len(rows) == 2


def test_skipped_gate_requires_reason(tmp_path):
    led = str(tmp_path / "ledger.jsonl")
    with pytest.raises(ValueError, match="reason"):
        cli._ledger_write(led, {"phase": "receipts", "verdict": "SKIPPED"})
    cli._ledger_write(led, {"phase": "receipts", "verdict": "SKIPPED",
                            "reason": "no --receipt-hook supplied"})


def test_non_gate_phases_unaffected(tmp_path):
    led = str(tmp_path / "ledger.jsonl")
    cli._ledger_write(led, {"phase": "preflight", "bulk_eqns": 42})
    cli._ledger_write(led, {"phase": "failed", "stage": "receipts"})
    assert len(open(led).read().splitlines()) == 2


# ------------------------------------------- D5: degenerate-input guards
def test_bulkmodel_negative_support_raises():
    with pytest.raises(ValueError, match="negative support"):
        bulkmodel.seeds_per_sector(3, 11, 8, -1, 0)
    with pytest.raises(ValueError, match="negative support"):
        bulkmodel.seeds_per_sector(3, 11, -8, 1, 0)


def test_bulkmodel_t_above_nsp_raises():
    with pytest.raises(ValueError, match="exceeds n_sp"):
        bulkmodel.seeds_per_sector(12, 11, 12, 0, 0)


def test_bulkmodel_benign_degenerates():
    assert bulkmodel.seeds_per_sector(9, 11, 8, 5, 0) == 0   # r < t sector
    assert bulkmodel.bulk([], 11, 18, 8, 5, 0) == 0          # empty census


def test_margins_const_arm_selfcontained(monkeypatch):
    """The const arm is fully self-contained: with NO model dir and NO
    labels anywhere, predict_margins(mode='const') emits the registered
    constant margin, and m2 mode still refuses loudly naming the env."""
    from seedling import margins as m, support
    monkeypatch.setattr(m, "_ML_DIR", "")
    monkeypatch.setattr(m, "DEFAULT_LABELS", "")
    monkeypatch.delenv("SEEDLING_ML_DIR", raising=False)
    monkeypatch.delenv("SEEDLING_LABELS", raising=False)
    census = {3: 2, 5: 2, 7: 3}
    targets = [support.Target("toy", (1, 1, 0))]
    tab = support.support_table(targets)
    rep = m.predict_margins(census, tab, targets, [7], mode="const", n_sp=3,
                            labels_path="")
    assert rep["b_const"] == m.B_CONST == 6
    assert rep["per_class"] and all(
        v["margin"] == 6 and v["raw"] is None and not v["floor_applied"]
        for v in rep["per_class"].values())
    assert rep["provenance"]["mode"] == "const"
    assert rep["s_floor"] == 1              # sub-top target: chain floor
    with pytest.raises(RuntimeError, match="SEEDLING_"):
        m.predict_margins(census, tab, targets, [7], mode="m2", n_sp=3,
                          labels_path="")


_NEEDS_MODEL = pytest.mark.skipif(
    not (os.environ.get("SEEDLING_ML_DIR") and os.environ.get("SEEDLING_LABELS")),
    reason="registered frozen model + labels not provided "
           "(set SEEDLING_ML_DIR and SEEDLING_LABELS)")


@_NEEDS_MODEL
def test_margins_nonfinite_raw_raises_loudly(monkeypatch):
    """A degenerate feature row reaching the frozen model must fail loudly,
    never emit a NaN-derived margin (sxx class). The model itself is frozen —
    we patch its predict at the module boundary only."""
    from seedling import margins as m
    frozen = m.fit_frozen()                     # real frozen model (cached)
    assert frozen["b_const"] == 6               # registered constant arm
    sys.path.insert(0, os.environ["SEEDLING_ML_DIR"])
    import train_eval
    monkeypatch.setattr(train_eval, "m2_predict",
                        lambda model, x: float("nan"))
    census = {3: 2, 7: 3}                       # toy census; top sector 7
    from seedling import support
    with tempfile.TemporaryDirectory() as td:
        tf = os.path.join(td, "targets")
        with open(tf, "w") as fh:
            fh.write("toy[1,1,0]\n")
        targets = support.parse_targets(tf)
    tab = support.support_table(targets)
    with pytest.raises(ValueError, match="non-finite"):
        m.predict_margins(census, tab, targets, [7], mode="m2", n_sp=3)


@_NEEDS_MODEL
def test_margins_const_mode_unaffected_by_predict(monkeypatch):
    from seedling import margins as m, support
    sys.path.insert(0, os.environ["SEEDLING_ML_DIR"])
    import train_eval
    monkeypatch.setattr(train_eval, "m2_predict",
                        lambda model, x: float("nan"))
    census = {3: 2, 7: 3}
    with tempfile.TemporaryDirectory() as td:
        tf = os.path.join(td, "targets")
        with open(tf, "w") as fh:
            fh.write("toy[1,1,0]\n")
        targets = support.parse_targets(tf)
    tab = support.support_table(targets)
    rep = m.predict_margins(census, tab, targets, [7], mode="const", n_sp=3)
    assert all(v["margin"] == 6 for v in rep["per_class"].values())
