"""Fail-closed contract of `seedling run`:
ANY stage exception => rc != 0 AND ledger phase=failed with stage name.
Synthetic fixtures only."""
import json
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", ".."))

from seedling import cli, margins, pinner  # noqa: E402


def _fixture(td):
    cfg = os.path.join(td, "config")
    os.makedirs(cfg)
    with open(os.path.join(cfg, "integralfamilies.yaml"), "w") as fh:
        fh.write("integralfamilies:\n - name: \"toy\"\n"
                 "   loop_momenta: [l]\n   top_level_sectors: [3]\n"
                 "   propagators:\n    - [\"l^2\", 0]\n"
                 "    - [\"(l+p1)^2\", 0]\n")
    with open(os.path.join(td, "targets"), "w") as fh:
        fh.write("toy[1,1]\n")
    with open(os.path.join(td, "nts"), "w") as fh:
        fh.write("3 2\n1 1\n")
    return cfg


def _run(td, out):
    return cli.main(["run", "--config", _fixture(td),
                     "--targets", os.path.join(td, "targets"),
                     "--out", out, "--nts", os.path.join(td, "nts")])


def _ledger(out):
    return [json.loads(l) for l in open(os.path.join(out, "ledger.jsonl"))]


def test_failclosed_margins_stage():
    orig = margins.predict_margins
    margins.predict_margins = lambda *a, **k: (_ for _ in ()).throw(
        RuntimeError("boom-margins"))
    try:
        with tempfile.TemporaryDirectory() as td:
            out = os.path.join(td, "run")
            rc = _run(td, out)
            led = _ledger(out)
            assert rc != 0, "swallowed rc on margins-stage exception"
            fail = [r for r in led if r["phase"] == "failed"]
            assert fail and fail[-1]["stage"] == "margins" \
                and "boom-margins" in fail[-1]["error"]
    finally:
        margins.predict_margins = orig


def test_failclosed_schedule_stage():
    orig = pinner.build_schedule
    orig_m = margins.predict_margins
    margins.predict_margins = lambda *a, **k: {
        "mode": "m2", "s_floor": 1, "b_const": 6,
        "per_class": {1: {"raw": 2.0, "margin": 2, "floor_applied": False},
                      2: {"raw": 2.5, "margin": 3, "floor_applied": False}},
        "provenance": {"model_hash": "stub"}}     # synthetic; no external IO
    pinner.build_schedule = lambda *a, **k: (_ for _ in ()).throw(
        RuntimeError("boom-schedule"))
    try:
        with tempfile.TemporaryDirectory() as td:
            out = os.path.join(td, "run")
            rc = _run(td, out)
            led = _ledger(out)
            assert rc != 0
            fail = [r for r in led if r["phase"] == "failed"]
            assert fail and fail[-1]["stage"] == "schedule" \
                and "boom-schedule" in fail[-1]["error"]
    finally:
        pinner.build_schedule = orig
        margins.predict_margins = orig_m


def test_missing_dep_error_names_the_dep():
    try:
        cli._ensure_dep("definitely_not_a_module_xyz", "nowhere")
        raised = False
    except RuntimeError as e:
        raised = "definitely_not_a_module_xyz" in str(e)
    assert raised


if __name__ == "__main__":
    test_failclosed_margins_stage()
    test_failclosed_schedule_stage()
    test_missing_dep_error_names_the_dep()
    print("test_run_failclosed: ALL PASS")
