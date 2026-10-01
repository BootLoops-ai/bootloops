#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision.
"""
pytest wrapper of the compare battery (tests/census_battery.py): one test per
family-level entry of the identity-census corpus (the verdict STRING equal to
the census's -- every pair the manifest carries, the dual cases inside their
entries), one per row-level rollup, the planted controls and the unit values
(the manifest's counts re-derived, the rows not compared absent from it),
plus the CLI surface (--drawn --verdict / --json 'verdict' block).
Cases that need networkx are marked nx_required (SKIPPED by name without
it); a missing fixture is SKIP by name, a drifted one a refusal.
"""
import json
import os
import re
import subprocess
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import census_battery as cb  # noqa: E402
import topology_audit as ta  # noqa: E402

PKG = os.path.dirname(HERE)
TOOL = os.path.join(PKG, "topology_audit.py")


def _man():
    try:
        return cb.load_manifest()
    except cb.FixtureMissing as exc:
        pytest.skip(str(exc))


MAN = _man()
ENTRY_IDS = [e["key"] for e in MAN["entries"]]
ROLLUP_IDS = [r["key"] for r in MAN["rollups"]]


def _run(fn):
    try:
        return fn()
    except cb.FixtureMissing as exc:
        pytest.skip(str(exc))


@pytest.mark.parametrize("key", ENTRY_IDS)
@cb.nx_required
def test_census_entry_verdict_string_equals_the_census(key):
    entry = [e for e in MAN["entries"] if e["key"] == key][0]
    r = _run(lambda: cb.check_entry(entry, MAN))
    assert r["got"] == r["expected"], (key, r["got"], r["expected"], r["line"])
    for x in r["extras"]:
        assert x["got"] == x["expected"], (key, x["drawn"], x["got"], x["expected"], x["line"])
    for a in r["also"]:
        assert a["got"] == a["expected"], (key, a["drawn"], a["got"], a["expected"], a["line"])
    assert r["result"]["verdict"] in ta.VERDICT_STRINGS
    assert r["result"]["line"].startswith("VERDICT: " + r["got"])


@pytest.mark.parametrize("key", ROLLUP_IDS)
@cb.nx_required
def test_census_rollup_string_joined_from_the_family_verdicts(key):
    roll = [r for r in MAN["rollups"] if r["key"] == key][0]
    case = cb._rollup_case(roll)
    r = _run(case)
    assert r["ok"], r


def test_manifest_is_the_census_by_object():
    _run(cb.unit_manifest_is_the_census_by_object)


def test_rows_not_compared_have_no_entry_rollup_or_fixture():
    _run(cb.unit_rows_not_compared_have_no_entry_rollup_or_fixture)


@cb.nx_required
def test_control_named_family_is_load_bearing():
    _run(cb.control_named_family_is_load_bearing)


def test_drawn_header_forms_read_by_object():
    _run(cb.unit_drawn_header_forms_read_by_object)


@cb.nx_required
def test_legs_per_vertex_reads_the_row27_photon_counts():
    _run(cb.unit_legs_per_vertex_reads_the_row27_photon_counts)


@cb.nx_required
def test_control_mass_moved():
    _run(cb.control_mass_moved_flips_to_family_mismatch_at_the_mass_level)


@cb.nx_required
def test_control_offshell_leg_moved():
    _run(cb.control_offshell_leg_moved_flips_to_family_mismatch_at_the_leg_level)


@cb.nx_required
def test_control_seven_line_top_sector_is_non_graph():
    _run(cb.control_seven_line_top_sector_is_non_graph)


@cb.nx_required
def test_control_undeclared_leg_classes_not_checkable():
    _run(cb.control_undeclared_leg_classes_are_not_checkable_with_four_legs)


@cb.nx_required
def test_control_two_point_undeclared_class_vacuous():
    _run(cb.control_two_point_undeclared_class_is_vacuous_and_said_by_name)


@cb.nx_required
def test_control_cyclic_order_mismatch_named_and_verdict_kept():
    _run(cb.control_cyclic_order_mismatch_is_a_named_sub_finding_that_keeps_the_verdict)


@cb.nx_required
def test_control_multi_integral_evidence_names_the_index_read():
    _run(cb.control_multi_integral_evidence_names_the_index_read)


@cb.nx_required
def test_control_multi_integral_tampered_expectation_fails_by_name():
    _run(cb.control_multi_integral_tampered_expectation_fails_by_name)


def _cli(args, cwd=PKG):
    return subprocess.run([sys.executable, TOOL] + args, capture_output=True, text=True, cwd=cwd)


@cb.nx_required
def test_cli_multi_integral_evidence_line_default_and_integral_flag():
    """Row 21 through the CLI: the default reading prints 'integral 0 of 3'
    beneath a FAMILY-MISMATCH verdict line, `--integral 2` prints 'integral 2
    of 3' beneath the census's IDENTITY-PASS; --json carries the four keys on
    both sides; a one-integral record (row 28) prints no such line."""
    e = [x for x in MAN["entries"] if x["key"] == "row21"][0]
    fam, drawn = cb.fixture(e["pair"]["family"]), cb.fixture(e["pair"]["drawn"])
    ints = [list(it["indices"]) for it in json.load(open(fam))["integrals"]]
    vs = cb.verdict_strings(MAN)["census"]
    p0 = _run(lambda: _cli([fam, "--drawn", drawn, "--verdict"]))
    assert p0.returncode == 0, p0.stderr[-800:]
    l0 = p0.stdout.splitlines()
    assert l0[0].startswith("VERDICT: " + vs["FAMILY-MISMATCH"]) and " of 3" not in l0[0]
    hit0 = [ln for ln in l0 if ln.strip().startswith("integral 0 of 3: %s; others: " % ints[0])]
    assert len(hit0) == 1 and "1 %s" % ints[1] in hit0[0] and "2 %s" % ints[2] in hit0[0], l0
    ii = e["pair"]["integral_index"]
    p2 = _run(lambda: _cli([fam, "--drawn", drawn, "--integral", str(ii), "--verdict"]))
    l2 = p2.stdout.splitlines()
    assert l2[0].startswith("VERDICT: " + e["pair"]["expected"]) and " of 3" not in l2[0]
    hit2 = [ln for ln in l2 if ln.strip().startswith("integral %d of 3: %s; others: " % (ii, ints[ii]))]
    assert len(hit2) == 1 and "0 %s" % ints[0] in hit2[0], l2
    j = json.loads(_run(lambda: _cli([fam, "--drawn", drawn, "--verdict", "--json"])).stdout)["verdict"]
    for tag in ("family", "drawn"):
        assert {"integral_index", "n_integrals", "integral_indices", "other_integrals"} <= set(j["sides"][tag])
    assert j["sides"]["family"]["integral_index"] == 0 and j["sides"]["family"]["n_integrals"] == 3
    assert j["sides"]["family"]["other_integrals"] == [[1, ints[1]], [2, ints[2]]]
    assert j["sides"]["drawn"]["n_integrals"] == 0
    e28 = [x for x in MAN["entries"] if x["key"] == "row28"][0]
    p28 = _run(lambda: _cli([cb.fixture(e28["pair"]["family"]), "--drawn", cb.fixture(e28["pair"]["drawn"]),
                             "--verdict"]))
    assert p28.stdout.splitlines()[0].startswith("VERDICT: " + e28["pair"]["expected"])
    assert not [ln for ln in p28.stdout.splitlines() if ln.strip().startswith("integral ") and " of " in ln]


@cb.nx_required
def test_cli_verdict_line_first_and_json_verdict_block():
    """`FAMILY --drawn DRAWN [--kinematics K] [--kinematics-drawn KD] --verdict`
    prints the verdict line first; --json carries the compare under 'verdict'
    (and the audit beside it without --verdict); the strings are the census's
    for the three verdict classes (one entry of each)."""
    def entry(key):
        return [e for e in MAN["entries"] if e["key"] == key][0]

    def args_for(e, extra_drawn=None):
        p = e["pair"]
        a = [cb.fixture(p["family"]), "--drawn", cb.fixture(extra_drawn or p["drawn"])]
        if p["kinematics"]:
            a += ["--kinematics", cb.fixture(p["kinematics"])]
        if p["drawn_kinematics"] and not extra_drawn:
            a += ["--kinematics-drawn", cb.fixture(p["drawn_kinematics"])]
        if p.get("integral_index") is not None:
            a += ["--integral", str(p["integral_index"])]
        return a

    counts = MAN["summary_verdict_counts"]
    picks = {}
    for e in MAN["entries"]:
        picks.setdefault(e["pair"]["expected"], e)
    assert set(picks) == set(counts)
    # + the entries whose census invocation named an integral index (--integral)
    # or a kinematics file on either side (--kinematics / --kinematics-drawn)
    with_index = [e for e in MAN["entries"] if e["pair"].get("integral_index") is not None]
    with_kin = [e for e in MAN["entries"] if e["pair"]["kinematics"] and e["pair"]["drawn_kinematics"]]
    assert with_index and with_kin
    picks[with_index[0]["pair"]["expected"] + " (with --integral)"] = with_index[0]
    picks[with_kin[0]["pair"]["expected"] + " (with kinematics both sides)"] = with_kin[0]
    for exp, e in picks.items():
        exp = e["pair"]["expected"]
        p = _run(lambda: _cli(args_for(e) + ["--verdict"]))
        assert p.returncode == 0, p.stderr[-800:]
        first = p.stdout.splitlines()[0]
        assert first.startswith("VERDICT: " + exp), (e["key"], first)
        j = _run(lambda: _cli(args_for(e) + ["--json"]))
        assert j.returncode == 0, j.stderr[-800:]
        d = json.loads(j.stdout)
        # the verdict string and the line's head are the contract; the leg map inside
        # the line is VF2's pick among a symmetric graph's valid maps, which varies
        # between processes (the double box has four), so it is not compared here
        assert d["verdict"]["verdict"] == exp and d["verdict"]["line"].startswith("VERDICT: " + exp)
        assert d["verdict"]["line"].split(" (")[0] == first.split(" (")[0]
        assert isinstance(d.get("audit"), list) and d["audit"] and "canonical_hash" in d["audit"][0]
        jv = _run(lambda: _cli(args_for(e) + ["--json", "--verdict"]))
        dv = json.loads(jv.stdout)
        assert "audit" not in dv and dv["verdict"]["verdict"] == exp
    # a pair that names its family (a multi-family file): the CLI carries the name
    # through --name (the battery's route, load_compare_side name=), the JSON's family
    # side names the family compared and says it was selected by name; the verdict is
    # the census's (the planted --name controls are test_cli_name_selects_... below)
    named = [e for e in MAN["entries"] if e["pair"].get("family_name")]
    assert named
    for e in named:
        fn = e["pair"]["family_name"]
        assert fn in e["pair"]["families_in_file"], (e["key"], e["pair"]["families_in_file"])
        d = json.loads(_run(lambda: _cli(args_for(e) + ["--name", fn, "--json", "--verdict"])).stdout)
        side = d["verdict"]["sides"]["family"]
        assert side["name"] == fn, (e["key"], side["name"])
        assert side["notes"]["family_name"] == fn and side["notes"]["family_selected_by"] == "name", side["notes"]
        assert side["notes"]["n_families_in_file"] == len(e["pair"]["families_in_file"]), side["notes"]
        assert d["verdict"]["verdict"] == e["pair"]["expected"], (e["key"], d["verdict"]["line"])
    # the dual case of row 24's corrected family: the same family, two drawn readings, two strings
    e24 = entry("row24/cured_adjacent_p3p4")
    x = e24["extras"][0]
    p1 = _run(lambda: _cli(args_for(e24) + ["--verdict"]))
    p2 = _run(lambda: _cli(args_for(e24, extra_drawn=x["drawn"]) + ["--verdict"]))
    assert p1.stdout.splitlines()[0].startswith("VERDICT: " + e24["pair"]["expected"])
    assert p2.stdout.splitlines()[0].startswith("VERDICT: " + x["expected"])
    assert e24["pair"]["expected"] != x["expected"]
    # --verdict / --kinematics-drawn without --drawn are refused by argparse
    bad = _cli([cb.fixture(e24["pair"]["family"]), "--verdict"])
    assert bad.returncode != 0 and "--drawn" in bad.stderr


@cb.nx_required
def test_cli_name_selects_the_family_compared_and_refuses_unknown():
    """The planted --name controls on row 18's two-family record (ggH_fam2 +
    ggH_fam5 in one integralfamilies.yaml): `--name <the census's family>`
    prints the census's IDENTITY-PASS; `--name <the other family>` prints
    FAMILY-MISMATCH at the bare level -- the verdict the battery's import
    route (run_pair / load_compare_side name=) gives for that name; a name
    the file does not carry is refused by name, rc != 0, the stderr listing
    the families the file carries and no verdict printed; the JSON verdict
    block names the family compared (sides.family.name, notes.family_name,
    family_selected_by); without --name the first family is compared and
    the notes say so."""
    e = [x for x in MAN["entries"] if x["key"] == "row18/record_family"][0]
    p = e["pair"]
    vs = cb.verdict_strings(MAN)["census"]
    names = p["families_in_file"]
    fn = p["family_name"]
    assert fn in names and len(names) >= 2, (fn, names)
    other = [n for n in names if n != fn][0]
    base = [cb.fixture(p["family"]), "--kinematics", cb.fixture(p["kinematics"]),
            "--drawn", cb.fixture(p["drawn"])]
    # (1) the named family: the census's verdict
    r1 = _run(lambda: _cli(base + ["--name", fn, "--verdict"]))
    assert r1.returncode == 0, r1.stderr[-800:]
    assert r1.stdout.splitlines()[0].startswith("VERDICT: " + p["expected"]), r1.stdout.splitlines()[0]
    assert "family %s (" % fn in r1.stdout, r1.stdout[:400]
    # (2) the other family of the same file: the import route's verdict for it
    _g2, res2, _w = _run(lambda: cb.run_pair(p["family"], p["kinematics"], p["drawn"],
                                             p["drawn_kinematics"], None, other))
    r2 = _run(lambda: _cli(base + ["--name", other, "--json", "--verdict"]))
    assert r2.returncode == 0, r2.stderr[-800:]
    d2 = json.loads(r2.stdout)["verdict"]
    assert d2["verdict"] == res2["verdict"] == vs["FAMILY-MISMATCH"] != p["expected"], (d2["line"], res2["line"])
    assert d2["first_failing_level"] == res2["first_failing_level"] == "bare", d2["line"]
    assert d2["E"] == res2["E"] and d2["V"] == res2["V"], (d2["E"], res2["E"])
    side2 = d2["sides"]["family"]
    assert side2["name"] == other and side2["notes"]["family_name"] == other, side2["notes"]
    assert side2["notes"]["family_selected_by"] == "name" and side2["notes"]["n_families_in_file"] == len(names)
    # (3) a name the file does not carry: refused by name, rc != 0, the names listed
    r3 = _cli(base + ["--name", "no_such_family", "--verdict"])
    assert r3.returncode != 0 and r3.stdout == "", (r3.returncode, r3.stdout[:200])
    assert "no_such_family" in r3.stderr and "no family found" in r3.stderr, r3.stderr[-400:]
    for n in names:
        assert n in r3.stderr, (n, r3.stderr[-400:])
    r3j = _cli(base + ["--name", "no_such_family", "--json", "--verdict"])
    assert r3j.returncode == r3.returncode and r3j.stdout == "", r3j.stdout[:200]
    # the audit path refuses the same way (rc != 0, names listed, no JSON)
    r3a = _cli([cb.fixture(p["family"]), "--name", "no_such_family", "--json"])
    assert r3a.returncode != 0 and r3a.stdout == "" and all(n in r3a.stderr for n in names), r3a.stderr[-400:]
    # (4) without --name: the file's first family, and the notes say so
    d4 = json.loads(_run(lambda: _cli(base + ["--json", "--verdict"])).stdout)["verdict"]
    assert d4["sides"]["family"]["name"] == names[0], d4["sides"]["family"]["name"]
    assert d4["sides"]["family"]["notes"]["family_selected_by"] == "first family in the file"
    assert d4["sides"]["family"]["notes"]["n_families_in_file"] == len(names)


@cb.nx_required
def test_self_test_leg_compare_is_registered_and_counts_the_battery():
    ids = [leg_id for leg_id, _t, _f in ta.EXTRA_SELF_TEST_LEGS]
    assert "compare" in ids
    fn = [f for leg_id, _t, f in ta.EXTRA_SELF_TEST_LEGS if leg_id == "compare"][0]
    assert fn is ta._selftest_compare
    out = ta._run_extra_leg(fn)
    assert out["status"] in ("PASS", "SKIP"), out["detail"]
    rows = out["results"]["census_battery"]
    assert len(rows) == len(cb.BATTERY)
    assert not [r for r in rows if r["status"] == "FAIL"], [r for r in rows if r["status"] == "FAIL"]


def test_no_verdict_string_is_typed_in_the_battery_expectations():
    """The battery's expected strings come from the manifest: the census's
    three strings appear in the battery source only as the vocabulary bridge
    (verdict_strings), never as an expected value; the manifest carries every
    expectation, each one of the census's strings."""
    src = open(os.path.join(HERE, "census_battery.py")).read()
    for s in MAN["summary_verdict_counts"]:
        hits = [m.start() for m in re.finditer(re.escape('"%s"' % s), src)]
        # allowed: the dict-key lookups vs["..."] (the vocabulary bridge)
        for h in hits:
            ctx = src[max(0, h - 4):h]
            assert ctx.endswith("vs[") or ctx.endswith("s["), (s, src[h - 40:h + 40])
    for e in MAN["entries"]:
        assert e["pair"]["expected"] in MAN["summary_verdict_counts"]
        for x in e["extras"] + e["also"]:
            assert x["expected"] in MAN["summary_verdict_counts"]
