#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision.
"""
Battery of the loomcheck member (loomcheck/): the Yangian/loom/fishnet
applicability screen for position-space conformal integrals.  One
hand-checkable graph per screened condition, asserted from the published
definitions and needing no data files:

  (2) planarity of the FULL graph, numerator edges counted: K5 with the edge
      (1,2) written as a numerator (power -1) is still K5, non-planar
      (Kuratowski), the numerator-edge list is [[1, 2]], its one internal
      vertex (5) weighs 4 = D, and the graph is excluded;
  (1) internal-vertex conformal weight: the 4-point star's vertex 5 weighs
      4 = D with unit powers and 5 != D with one propagator squared;
  (3) the level-one-momentum face condition sum_face a_e = (k-2) D/2: a unit
      4-cycle has two quadrilateral faces, both satisfying it (0 violations);
      a unit 5-cycle's two pentagon faces sum to 5 != 6 (2 violations) and
      the graph is excluded;
  the summary invariant yangian_class == (weight and planar and 0 face
  violations) on each; the basis-file parser on a synthetic two-entry file in
  the arXiv:2607.11645 ancillary format (entry = numerator/(denominators),
  single-digit vertex labels); the CLI surface (`-m loomcheck --selftest`,
  the script form, `--basis/--targets/--json`); the package self-test's
  [loomcheck] leg registered and green.  networkx absent -> the graph cases
  SKIP by name (the package's declared dependency).
"""
import json
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
TOOL_DIR = os.path.dirname(HERE)
sys.path.insert(0, TOOL_DIR)
import topology_audit as ta  # noqa: E402

try:
    import pytest
except ImportError:  # the __main__ battery runner does not need pytest
    pytest = None

NX_SKIP = ta.NX_SKIP_REASON


class NxMissing(Exception):
    pass


def _lc():
    """The member module, imported as the `loomcheck` package from the tool dir."""
    if not ta._HAVE_NX:
        raise NxMissing(NX_SKIP)
    import loomcheck  # noqa: E402  (tools/dogtag on sys.path)
    return loomcheck


def nx_required(fn):
    if pytest is None:
        return fn
    return pytest.mark.skipif(not ta._HAVE_NX, reason=NX_SKIP)(fn)


K5_NUM = {**{(i, j): 1 for i in range(1, 6) for j in range(i + 1, 6)}, (1, 2): -1}
STAR = {(1, 5): 1, (2, 5): 1, (3, 5): 1, (4, 5): 1}
SQUARE = {(1, 2): 1, (2, 3): 1, (3, 4): 1, (1, 4): 1}
PENTAGON = {(1, 2): 1, (2, 3): 1, (3, 4): 1, (4, 5): 1, (1, 5): 1}
KEYS = {"conformal_internal_weight", "planar_full_graph", "numerator_edges",
        "faces", "phat_face_violations", "yangian_class"}


@nx_required
def test_k5_with_a_numerator_edge_is_non_planar_and_excluded():
    lc = _lc()
    r = lc.screen_graph(dict(K5_NUM), D=4)
    assert r["planar_full_graph"] is False
    assert r["numerator_edges"] == [[1, 2]]
    assert r["conformal_internal_weight"] is True      # vertex 5: four unit lines
    assert r["faces"] is None and r["phat_face_violations"] is None
    assert r["yangian_class"] is False
    assert set(r) == KEYS


@nx_required
def test_star_internal_vertex_weight_unit_vs_squared():
    lc = _lc()
    assert lc.screen_graph(dict(STAR), D=4)["conformal_internal_weight"] is True
    assert lc.screen_graph({**STAR, (1, 5): 2}, D=4)["conformal_internal_weight"] is False
    # explicit internal list overrides the ">= 5" convention: vertex 1 has weight 1 != 4
    assert lc.screen_graph(dict(STAR), D=4, internal=[1])["conformal_internal_weight"] is False
    # D is a parameter: unit star weighs 4, so D=3 fails and D=4 passes
    assert lc.screen_graph(dict(STAR), D=3)["conformal_internal_weight"] is False


@nx_required
def test_face_condition_quadrilateral_passes_pentagon_fails():
    lc = _lc()
    sq = lc.screen_graph(dict(SQUARE), D=4)
    assert sq["planar_full_graph"] is True and sq["faces"] == 2
    assert sq["phat_face_violations"] == 0
    assert sq["yangian_class"] is True          # no internal vertices, planar, faces fine
    pent = lc.screen_graph(dict(PENTAGON), D=4)
    assert pent["faces"] == 2 and pent["phat_face_violations"] == 2
    assert pent["yangian_class"] is False


@nx_required
def test_yangian_class_is_the_conjunction_of_the_three_conditions():
    lc = _lc()
    for g in (K5_NUM, STAR, SQUARE, PENTAGON, {**STAR, (1, 5): 2}):
        r = lc.screen_graph(dict(g), D=4)
        want = bool(r["conformal_internal_weight"] and r["planar_full_graph"]
                    and r["phat_face_violations"] == 0)
        assert r["yangian_class"] == want, (g, r)


@nx_required
def test_selftest_cases_all_pass_and_selftest_rc_zero():
    lc = _lc()
    cases = lc.selftest_cases(D=4)
    assert len(cases) >= 10
    bad = [(n, g, w) for n, g, w in cases if g != w]
    assert not bad, bad
    assert lc.selftest() == 0


def _synthetic_basis(path):
    """Two entries in the ancillary format: a header comment, then
    {{expr/(dens), tag}, ...}.  Entry 1 = the unit 4-cycle on the external
    vertices 1..4 (dens x12 x23 x34 x14: planar, two quadrilateral faces, no
    internal vertex -- in class by the three conditions); entry 2 = the
    4-point star with x15 squared and a numerator x12 (weight at the internal
    vertex 5: 2+1+1+1 = 5 != 4 -- excluded by condition (1))."""
    txt = ("(* synthetic basis in the arXiv:2607.11645 ancillary format *)\n"
           "{{1/(x[1,2] x[2,3] x[3,4] x[1,4]), I[1]},\n"
           " {x[1,2]/(x[1,5]^2 x[2,5] x[3,5] x[4,5]), I[2]}}\n")
    with open(path, "w") as fh:
        fh.write(txt)


@nx_required
def test_parser_on_a_synthetic_basis_file():
    lc = _lc()
    with tempfile.TemporaryDirectory() as td:
        p = os.path.join(td, "basis.m")
        _synthetic_basis(p)
        entries = lc.load_entries(p)
        assert len(entries) == 2
        assert lc.parse(entries, 1) == dict(SQUARE)
        p2 = lc.parse(entries, 2)
        assert p2 == {(1, 5): 2, (2, 5): 1, (3, 5): 1, (4, 5): 1, (1, 2): -1}
        r1 = lc.screen_graph(lc.parse(entries, 1), D=4)
        r2 = lc.screen_graph(p2, D=4)
        assert r1 == lc.screen_graph(dict(SQUARE), D=4) and r1["yangian_class"] is True
        assert r2["conformal_internal_weight"] is False and r2["numerator_edges"] == [[1, 2]]
        assert r2["yangian_class"] is False


def _run(args, cwd=TOOL_DIR, env=None):
    e = dict(os.environ)
    if env:
        e.update(env)
    return subprocess.run([sys.executable] + args, capture_output=True, text=True, cwd=cwd, env=e)


@nx_required
def test_cli_selftest_module_and_script_forms():
    lc = _lc()
    m = _run(["-m", "loomcheck", "--selftest"])
    assert m.returncode == 0, m.stdout[-600:] + m.stderr[-600:]
    assert "SELFTEST PASS: " in m.stdout and "FAIL " not in m.stdout
    n = len(lc.selftest_cases(D=4))
    assert f"SELFTEST PASS: {n}/{n} checks passed" in m.stdout
    s = _run([os.path.join(TOOL_DIR, "loomcheck", "loomcheck.py"), "--selftest"], cwd=HERE)
    assert s.returncode == 0, s.stdout[-600:] + s.stderr[-600:]
    assert s.stdout == m.stdout


@nx_required
def test_cli_basis_targets_json_receipt_and_env_default():
    lc = _lc()
    with tempfile.TemporaryDirectory() as td:
        p = os.path.join(td, "basis.m")
        _synthetic_basis(p)
        out = os.path.join(td, "receipt.json")
        r = _run(["-m", "loomcheck", "--basis", p, "--targets", "1,2", "--json", out])
        assert r.returncode == 0, r.stderr[-600:]
        assert "I1:" in r.stdout and "I2:" in r.stdout
        j = json.load(open(out))
        assert set(j) == {"I1", "I2"}
        # the receipt is the API's screen of the parsed entries, key by key
        entries = lc.load_entries(p)
        for t in (1, 2):
            api = lc.screen_graph(lc.parse(entries, t), D=4)
            assert j[f"I{t}"] == json.loads(json.dumps(api)), (t, j[f"I{t}"], api)
        n_in = sum(1 for v in j.values() if v["yangian_class"])
        assert f"SUMMARY: {n_in}/2 in a proven Yangian/SoV class" in r.stdout
        assert n_in == 1 and j["I1"]["yangian_class"] is True and j["I2"]["yangian_class"] is False
        # LOOMCHECK_BASIS is the --basis default
        r2 = _run(["-m", "loomcheck", "--targets", "1"], env={"LOOMCHECK_BASIS": p})
        assert r2.returncode == 0 and "SUMMARY: 1/1 in a proven" in r2.stdout, r2.stderr[-600:]
        # the script form from another cwd gives the same lines
        r3 = _run([os.path.join(TOOL_DIR, "loomcheck", "loomcheck.py"), "--basis", p,
                   "--targets", "1,2"], cwd=td)
        assert r3.returncode == 0 and r3.stdout == r.stdout, r3.stderr[-600:]


@nx_required
def test_self_test_leg_loomcheck_is_registered_and_green():
    _lc()
    legs = dict((lid, fn) for lid, _t, fn in ta.EXTRA_SELF_TEST_LEGS)
    assert "loomcheck" in legs and legs["loomcheck"] is ta._selftest_loomcheck
    out = legs["loomcheck"]()
    assert out["ok"] is True and out["skipped"] == []
    assert out["executed"] >= 14
    rows = out["results"]["loomcheck_battery"]
    assert not [r for r in rows if r["status"] != "PASS"], rows
    wrapped = ta._run_extra_leg(ta._selftest_loomcheck)
    assert wrapped["status"] == "PASS", wrapped["detail"]


def test_member_layout_is_in_the_package():
    """The member is a subfolder of this package (no top-level tools/loomcheck):
    loomcheck/{__init__,__main__,loomcheck}.py + README.md beside topology_audit.py."""
    d = os.path.join(TOOL_DIR, "loomcheck")
    for f in ("__init__.py", "__main__.py", "loomcheck.py", "README.md"):
        assert os.path.isfile(os.path.join(d, f)), f
    assert not os.path.exists(os.path.join(os.path.dirname(TOOL_DIR), "loomcheck"))


# ---------------------------------------------------------------------------
#  battery runner (python3 tests/test_loomcheck.py): PASS / SKIP (named) / FAIL
# ---------------------------------------------------------------------------
BATTERY = [
    test_k5_with_a_numerator_edge_is_non_planar_and_excluded,
    test_star_internal_vertex_weight_unit_vs_squared,
    test_face_condition_quadrilateral_passes_pentagon_fails,
    test_yangian_class_is_the_conjunction_of_the_three_conditions,
    test_selftest_cases_all_pass_and_selftest_rc_zero,
    test_parser_on_a_synthetic_basis_file,
    test_cli_selftest_module_and_script_forms,
    test_cli_basis_targets_json_receipt_and_env_default,
    test_self_test_leg_loomcheck_is_registered_and_green,
    test_member_layout_is_in_the_package,
]


def run_battery(verbose=True):
    rows = []
    for fn in BATTERY:
        name = fn.__name__
        try:
            fn()
            rows.append({"name": name, "status": "PASS", "detail": ""})
        except NxMissing as exc:
            rows.append({"name": name, "status": "SKIP", "detail": str(exc)})
        except Exception as exc:
            rows.append({"name": name, "status": "FAIL", "detail": f"{type(exc).__name__}: {exc}"})
        if verbose:
            r = rows[-1]
            print(f"   {r['status']:4s} {name}" + (f"  ({r['detail']})" if r["detail"] else ""), flush=True)
    return rows


if __name__ == "__main__":
    out = run_battery(verbose=True)
    n_fail = sum(1 for r in out if r["status"] == "FAIL")
    n_skip = sum(1 for r in out if r["status"] == "SKIP")
    n_pass = sum(1 for r in out if r["status"] == "PASS")
    print(f"loomcheck battery: {len(out)} cases, {n_pass} PASS, {n_fail} FAIL, {n_skip} SKIP")
    sys.exit(1 if n_fail else (2 if n_skip else 0))
