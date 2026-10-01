#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision.
"""
Planarity battery for topology_audit.py: the closure-cycle test names the
cyclic order it closes the legs in and where that order came from (the
caller's leg_order, the record's cyclic_leg_order, else the canonical index
order as a named convention); with fewer than four attached legs it says
"not discriminating (N < 4 legs)" by name (the method string, the
PLANARITY-NOT-DISCRIMINATING line, planarity_discriminating False) instead of
printing a bare planar=True; the identity census's legs-joined-at-infinity
criterion (planarity_at_infinity, the census's planarity_closure form) is a
second reported field; the cyclic leg orders are enumerated; the record's
Mandelstam declaration is reported as a candidate order with its own closure
verdict, never enforced.  Every case reproduces a hand-checked object of the
identity census (IDENTITY_CENSUS.json .rows[].fingerprint,
rows/<row>/audit_family_import.json, catalog_check.json, ROW.md sentences),
copied, never paraphrased; every planted control flips or is refused by name.

  row 24  retired hexabox + kinematics_of_record.yaml, cured JSON (p5
          implicit): the census's planarity_on_realized_graph
          {"planar_canonical_leg_order": false, "abstract_planar": true,
           "planar_legs_joined_at_infinity": false,
           "closure_order": ["p1", "p2", "p3", "p4", "p5"], "discriminating": true}
          reproduced by planarity_at_infinity; audit planar False AND
          planar_legs_joined_at_infinity False ("planar with all legs on the
          outer face: False (the internal graph alone is planar: True)")
  rows 15-17  icc/iccg + the record kinematics: planar True (census
          "planar": true, method "... closure cycle in canonical order
          1,2,3") now flagged "not discriminating (N < 4 legs)" by name
  row 23  kite_thr + kinematics: two legs, planar True (census method
          "networkx.check_planarity (abstract; <3 external legs)") flagged
  row 18  fam2 (three legs, inferred set): verdict withheld, the report
          carries discriminating False and n_external_legs 3
  rows 8, 13  two-point JSON legs as data (flagged) / yaml inferred (withheld)
  row 20  A1 + family_kinematics.yaml: planar False, {s: 2, t: 3, u: 3},
          size 7 UNCHANGED (the kinematics contract); legs at infinity False;
          the declaration s = (p1+p2)^2, t = (p1+p3)^2 implies p1,p2,p4,p3,
          where the graph is also non-planar (no leg-order line)
  row 2   the family of record (t = (p1+p3)^2 convention): planar False in
          the canonical order (census "planar": false, method "... canonical
          order 1,2,3,4"), planar True in its own order p1,p2,p4,p3 -- the
          census: "the graph is planar with its legs in the drawn order under
          the leg map above, which sends the family's t=(p1+p3)^2 to the
          paper's t=(p2+p3)^2"; reported under PLANARITY-LEG-ORDER
  theta   the (2,1,1) crossed box: catalog_check.json
          "planar_with_legs_for_some_leg_order": false, "abstract_planar":
          true -- now by enumeration (0 of 3 orders) and at infinity False
  row 28  lbl3se + its Kira kinematics form / the JSON: planar True in the
          drawn order p1,p2,p3,p4 (census "planar_with_legs_closed_in_order":
          true, "closure_order": ["p1", "p2", "p3", "p4"]) while its
          t = (p1+p3)^2 declaration implies p1,p2,p4,p3 (planar False there):
          reported, not enforced -- the evidence that a record's channel
          names do not fix the drawn boundary order
  row 5   the pySecDec graph's leg list is the closure order
          (family.cyclic_leg_order); planar True, 1e7dddf35ab0
  row 35  wpairT4 + kinematics: s = (p1+p2)^2, t = (p2+p3)^2 implies the
          canonical order (agrees); planar False {s: 2, t: 3, u: 3}
  controls  a drawn one-loop box with cyclic_leg_order p1,p2,p4,p3 -> planar
          False (p1,p2,p3,p4 -> True); a leg_order that is not a permutation
          -> REFUSED by name; K3,3 and prism three-point graphs where the
          closure test decides False under the flag; the retired
          legs-3-4-interchanged entry planar in the order p1,p2,p4,p3 and the
          planar box non-planar in it; every 2-/3-leg catalog entry flagged,
          every 4-leg entry's method string unchanged

Runs under pytest (test_* functions) and as a battery from the tool's
--self-test leg [planarity] (run_battery()).  A missing fixture is a SKIP by name;
a drifted fixture is a refusal (FixtureDrift), never a pass.
"""
import functools
import hashlib
import json
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
TOOL_DIR = os.path.dirname(HERE)
sys.path.insert(0, TOOL_DIR)
import topology_audit as ta  # noqa: E402

FIX = os.path.join(HERE, "fixtures")

# sha16 of every vendored fixture this battery reads (PROVENANCE.json in the
# row directory carries the full sha256 of the copy on disk and the record's
# own sha256; both are checked on every access).
PINS = {
    "row02/family_of_record.jl": "ddc75ae24954b081",
    "row02/genuine_crossed_box_theta211.yaml": "05d1300140c4b9ca",
    "row02/catalog_check.json": "892422b5b8ba0fac",
    "row05/family_of_record.yaml": "0ec539a5fec61bdf",
    "row05/psd_build_pentagon.py": "70d6e3fc61ffb04a",
    "row06/retired_sevenline_vertex2L/family_of_record_transcribed.yaml": "e759606f59374bd6",
    "row08/family_of_record.json": "4109f261715a4cfe",
    "row13/family_of_record.yaml": "00fa399594fc002a",
    "row15/family_of_record_kira_icc_eqmass.yaml": "9d0e526a98317304",
    "row15/kinematics_of_record.yaml": "961466e57b445b71",
    "row16/family_of_record_kira_iccg_generic.yaml": "5be2fe7689e8dee6",
    "row17/family_of_record_kira_icc_eqmass.yaml": "9d0e526a98317304",
    "row18/family_of_record_fam2.yaml": "c1c9c4e3ea6eb708",
    "row20/family_of_record.yaml": "d4aef7e8609a1f65",
    "row20/family_kinematics.yaml": "9fb7631fb6a8455a",
    "row23/family_of_record_kira_kite_thr_detransport.yaml": "53657fedefba6137",
    "row23/kinematics_of_record.yaml": "5ee9f2717eca1c02",
    "row24/retired_nonadjacent_p1p3/family_of_record.yaml": "9cb3b9b0d49d4ceb",
    "row24/retired_nonadjacent_p1p3/kinematics_of_record.yaml": "64d62e86c0b095e8",
    "row24/cured_adjacent_p3p4/family_of_record.json": "1864217f1e908a08",
    "row28/family_of_record.yaml": "f07570cc5983b64d",
    "row28/family_of_record_kira_kinematics.yaml": "512fc1c966e86e04",
    "row28/family_of_record.json": "61c4e46fb523550b",
    "row35/family_of_record.yaml": "3d3f3825b7b83326",
    "row35/kinematics_of_record.yaml": "63aee0ed75326706",
}

# The census objects, copied (rows/<row>/... of the identity census).
CENSUS = {
    # rows/row24/retired_nonadjacent_p1p3/audit_family_import.json and
    # rows/row24/cured_adjacent_p3p4/audit_family_import.json,
    # "planarity_on_realized_graph"
    "row24_planarity_on_realized_graph": {
        "planar_canonical_leg_order": False, "abstract_planar": True,
        "planar_legs_joined_at_infinity": False,
        "closure_order": ["p1", "p2", "p3", "p4", "p5"], "discriminating": True},
    # rows/row24/retired_nonadjacent_p1p3/ROW.md
    "row24_sentence": "planar with all legs on the outer face: False (the internal graph "
                      "alone is planar: True), so the NP label holds.",
    # IDENTITY_CENSUS.json .rows["15"|"16"|"17"].fingerprint
    "rows15_17": {"planar": True, "planarity_method": "networkx.check_planarity with "
                  "external-leg closure cycle in canonical order 1,2,3",
                  "canonical_hash_normalized": "8a73188ccb5a"},
    # .rows["23"].fingerprint
    "row23": {"planar": True, "planarity_method": "networkx.check_planarity (abstract; "
              "<3 external legs)", "canonical_hash_normalized": "bb2ae2aad8f4"},
    # .rows["18"].fingerprint (the census import route: p3 added)
    "row18": {"planar": True, "canonical_hash_normalized": "996d57175c92"},
    # .rows["20"].fingerprint.family_of_record (kinematics.yaml applied)
    "row20": {"planar": False, "cut_signature": {"s": 2, "t": 3, "u": 3},
              "canonical_hash": "3c695d99a1e1"},
    # .rows["2"].fingerprint + rows/row02/ROW.md
    "row02": {"planar": False, "planarity_method": "networkx.check_planarity with "
              "external-leg closure cycle in canonical order 1,2,3,4",
              "canonical_hash": "a542a386080e", "cut_signature": {"s": 2, "t": 4, "u": 3}},
    "row02_sentence": ("the family's legs follow Smirnov's convention t=(p1+p3)^2 (boundary "
                       "order p1,p2,p4,p3), so the checker's canonical-order closure test "
                       "reports planar=False"),
    "row02_drawn_order_sentence": ("the graph is planar with its legs in the drawn order under "
                                   "the leg map above, which sends the family's t=(p1+p3)^2 to "
                                   "the paper's t=(p2+p3)^2"),
    # rows/row02/catalog_check.json "genuine_crossed_box_theta211"
    "theta": {"planar_with_legs_for_some_leg_order": False, "abstract_planar": True,
              "cut_signature": {"s": 2, "t": 3, "u": 3}, "canonical_hash": "7790fb92864c"},
    # .rows["28"].fingerprint.connected_realization
    "row28": {"planar_with_legs_closed_in_order": True,
              "closure_order": ["p1", "p2", "p3", "p4"],
              "import_route_canonical_hash": "34c5ba1eaf82"},
    # .rows["5"].fingerprint
    "row05": {"planar": True, "planarity_method": "networkx.check_planarity with "
              "external-leg closure cycle in canonical order 1,2,3,4,5",
              "canonical_hash": "1e7dddf35ab0"},
    # .rows["35"].fingerprint.connected_realization
    "row35": {"planar_with_legs_closed_in_order": False,
              "closure_order": ["p1", "p2", "p3", "p4"],
              "cut_signature": {"s": 2, "t": 3, "u": 3}},
}

FLAG = "not discriminating (N < 4 legs)"
NOT_DISC = "PLANARITY-NOT-DISCRIMINATING"
LEG_ORDER = "PLANARITY-LEG-ORDER"
INFERRED = "LEG-SET-INFERRED"
CANONICAL_SRC = ta._PLANARITY_CANONICAL_ORDER_SOURCE
CYCLIC_SRC = "family.cyclic_leg_order (the record's drawn / declared cyclic order)"
ARG_SRC = "leg_order argument (the caller's declared cyclic order)"


class FixtureDrift(Exception):
    """A vendored fixture no longer matches its pin: refused."""


class FixtureMissing(Exception):
    """A vendored fixture is absent: the case is a SKIP by name."""


def nx_required(fn):
    """Mark a case that needs networkx (a declared dependency of the tool):
    pytest SKIPs it by name when the dependency is absent (skipif), and
    run_battery -- the tool's self-test leg -- reports it SKIP with the same
    reason through the missing-fixture path (FixtureMissing): counted, named,
    never PASS, never FAIL."""
    @functools.wraps(fn)
    def gated():
        why = ta._require_nx()
        if why:
            raise FixtureMissing(why)
        return fn()
    try:
        import pytest
    except ImportError:                                            # pragma: no cover
        return gated
    return pytest.mark.skipif(not ta._HAVE_NX, reason=ta.NX_SKIP_REASON)(gated)


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        h.update(f.read())
    return h.hexdigest()


def fixture(rel):
    """Path of a vendored fixture after checking its sha against PINS and
    the row's PROVENANCE.json; raises FixtureMissing / FixtureDrift."""
    path = os.path.join(FIX, rel)
    if not os.path.exists(path):
        raise FixtureMissing(f"fixture absent: tests/fixtures/{rel}")
    if rel not in PINS:
        raise FixtureDrift(f"fixture not pinned: {rel}")
    sha = _sha256(path)
    if sha[:16] != PINS[rel]:
        raise FixtureDrift(f"fixture drift: {rel} sha16 {sha[:16]} != pin {PINS[rel]}")
    row, within = rel.split("/", 1)
    prov_path = os.path.join(FIX, row, "PROVENANCE.json")
    if not os.path.exists(prov_path):
        raise FixtureDrift(f"PROVENANCE.json absent for {row}")
    prov = json.load(open(prov_path))
    entry = None
    for e in prov.get("files", []) or []:
        if e.get("path") == within:
            entry = e
    if entry is None or entry.get("sha256") != sha:
        raise FixtureDrift(f"fixture drift: {rel} not recorded with sha256 {sha} "
                           f"in {row}/PROVENANCE.json")
    return path


def load(rel, kin=None, **kw):
    return ta.load_family(fixture(rel), kinematics=(fixture(kin) if kin else None), **kw)[0]


def audit(rel, kin=None, leg_perms=False):
    return ta.audit(load(rel, kin), try_leg_perms=leg_perms)


def has(fp, tag):
    return any(w.startswith(tag) for w in fp["warnings"])


def fam_from(name, spec):
    return ta.Family(name, spec["loops"], spec["exts"], spec["ext_subs"],
                     spec["propagators"], "memory:" + name)


def _tmpdir():
    return tempfile.mkdtemp(prefix="topology_audit_planarity_")


def _spec_family(spec, name):
    """A Family through the drawn-graph (edge-list) reader from an in-memory
    spec: the route a drawn graph's cyclic_leg_order takes."""
    d = _tmpdir()
    p = os.path.join(d, name + ".json")
    with open(p, "w") as f:
        json.dump(spec, f)
    try:
        return ta.load_family(p)[0]
    finally:
        shutil.rmtree(d, ignore_errors=True)


BOX_SPEC = {"name": "box_drawn", "vertices": ["a", "b", "c", "d"],
            "edges": [["a", "b", 0], ["b", "c", 0], ["c", "d", 0], ["d", "a", 0]],
            "legs": [["p1", "a"], ["p2", "b"], ["p3", "c"], ["p4", "d"]],
            "loop_momenta": ["k1"]}
K33_SPEC = {"name": "k33_3pt", "vertices": ["a1", "a2", "a3", "b1", "b2", "b3"],
            "edges": [[a, b, 0] for a in ("a1", "a2", "a3") for b in ("b1", "b2", "b3")],
            "legs": [["p1", "a1"], ["p2", "b1"], ["p3", "a2"]],
            "loop_momenta": ["k1", "k2", "k3", "k4"]}
PRISM_SPEC = {"name": "prism_3pt", "vertices": ["t1", "t2", "t3", "b1", "b2", "b3"],
              "edges": [["t1", "t2", 0], ["t2", "t3", 0], ["t3", "t1", 0],
                        ["b1", "b2", 0], ["b2", "b3", 0], ["b3", "b1", 0],
                        ["t1", "b1", 0], ["t2", "b2", 0], ["t3", "b3", 0]],
              "legs": [["p1", "t1"], ["p2", "t2"], ["p3", "b3"]],
              "loop_momenta": ["k1", "k2", "k3", "k4"]}


# ---------------------------------------------------------------- pins --
def test_fixture_provenance_pins():
    """Every file this battery pins: sha256 in code == PROVENANCE.json ==
    bytes on disk (refused on drift)."""
    seen = 0
    for rel in sorted(PINS):
        fixture(rel)
        seen += 1
    assert seen == len(PINS)


# ------------------------------------------------- row 24: legs at infinity --
@nx_required
def test_row24_retired_reproduces_the_census_planarity_closure_dict():
    fam = load("row24/retired_nonadjacent_p1p3/family_of_record.yaml",
               "row24/retired_nonadjacent_p1p3/kinematics_of_record.yaml")
    assert fam.exts == ["p1", "p2", "p3", "p4", "p5"]
    got = ta.planarity_at_infinity(fam)
    assert got == CENSUS["row24_planarity_on_realized_graph"], got
    fp = ta.audit(fam)
    # the sentence: "planar with all legs on the outer face: False (the internal
    # graph alone is planar: True)" -- both fields False on the connected cover
    assert fp["planar"] is False and fp["planar_legs_joined_at_infinity"] is False
    assert fp["planarity_report"]["abstract_planar"] is True
    assert fp["planarity_discriminating"] is True
    assert fp["planarity_closure_order"] == ["p1", "p2", "p3", "p4", "p5"]
    assert fp["planarity_closure_order_source"] == CANONICAL_SRC
    assert fp["planarity_report"]["n_leg_orders"] == 12         # (5-1)!/2 cyclic orders
    assert fp["planarity_report"]["n_planar_leg_orders"] == 0
    assert fp["planarity_report"]["planar_in_some_leg_order"] is False
    assert not has(fp, NOT_DISC) and not has(fp, INFERRED)
    # the declaration s12, s23, s34 implies the canonical order (agrees)
    md = fp["planarity_report"]["mandelstam_declaration"]
    assert md["cyclic_order"] == ["p1", "p2", "p3", "p4", "p5"] and md["agrees_with_closure_order"] is True
    assert not has(fp, LEG_ORDER)


@nx_required
def test_row24_cured_json_reproduces_the_census_planarity_closure_dict():
    fam = load("row24/cured_adjacent_p3p4/family_of_record.json")
    assert fam.exts == ["p1", "p2", "p3", "p4", "p5"]          # p5 completed from conservation
    assert ta.planarity_at_infinity(fam) == CENSUS["row24_planarity_on_realized_graph"]
    fp = ta.audit(fam)
    assert fp["planar"] is False and fp["planar_legs_joined_at_infinity"] is False
    assert fp["planarity_discriminating"] is True and not has(fp, INFERRED)


# ------------------------------------------ rows 15-17, 23: N < 4 by name --
@nx_required
def test_rows15_16_17_with_kinematics_are_flagged_not_discriminating():
    e = CENSUS["rows15_17"]
    for rel in ("row15/family_of_record_kira_icc_eqmass.yaml",
                "row16/family_of_record_kira_iccg_generic.yaml",
                "row17/family_of_record_kira_icc_eqmass.yaml"):
        fp = audit(rel, "row15/kinematics_of_record.yaml")
        assert fp["canonical_hash"] == e["canonical_hash_normalized"]
        assert fp["planar"] is e["planar"]                       # the census value, kept
        assert fp["planarity_method"].startswith(FLAG + ": "), fp["planarity_method"]
        assert "external-leg closure cycle in canonical order" in fp["planarity_method"]
        assert fp["planarity_discriminating"] is False
        assert fp["planarity_report"]["n_external_legs"] == 3
        assert fp["planar_legs_joined_at_infinity"] is True
        assert has(fp, NOT_DISC), fp["warnings"]
        assert not has(fp, INFERRED)
        line = [w for w in fp["warnings"] if w.startswith(NOT_DISC)][0]
        assert "3 external leg(s) attached (N < 4)" in line and "planar=True" in line


@nx_required
def test_row23_with_kinematics_two_legs_flagged():
    e = CENSUS["row23"]
    fp = audit("row23/family_of_record_kira_kite_thr_detransport.yaml",
               "row23/kinematics_of_record.yaml")
    assert fp["canonical_hash"] == e["canonical_hash_normalized"]
    assert fp["planar"] is e["planar"]
    assert fp["planarity_method"] == FLAG + ": " + e["planarity_method"]
    assert fp["planarity_discriminating"] is False
    assert fp["planarity_report"]["n_external_legs"] == 2
    assert fp["planar_legs_joined_at_infinity"] is True
    assert has(fp, NOT_DISC) and not has(fp, INFERRED)


@nx_required
def test_row18_fam2_three_legs_inferred_is_withheld_with_the_report():
    fp = audit("row18/family_of_record_fam2.yaml")
    assert has(fp, INFERRED) and fp["planar"] is None
    assert fp["planarity_method"].startswith("withheld: LEG-SET-INFERRED")
    assert fp["planar_legs_joined_at_infinity"] is None            # withheld with the verdict
    rep = fp["planarity_report"]
    assert rep["withheld"].startswith(INFERRED)
    assert rep["discriminating"] is False and rep["n_external_legs"] == 3
    assert not has(fp, NOT_DISC)                                    # no verdict, no flag line
    # the value the census's import route printed (p3 added): planar True, flagged
    assert fp["planarity_on_inferred_leg_set"]["planar"] is CENSUS["row18"]["planar"]
    assert fp["planarity_on_inferred_leg_set"]["method"].startswith(FLAG + ": ")
    pl = ta.planarity(load("row18/family_of_record_fam2.yaml"))
    assert pl["planar"] is True and pl["discriminating"] is False
    assert pl["method"] == (FLAG + ": networkx.check_planarity with external-leg "
                            "closure cycle in canonical order 1,2,3")


@nx_required
def test_rows08_13_two_point_json_legs_flagged_yaml_inferred_withheld():
    fj = audit("row08/family_of_record.json")                       # legs as data
    assert fj["planar"] is True and not has(fj, INFERRED)
    assert fj["planarity_method"] == (FLAG + ": networkx.check_planarity (abstract; "
                                      "<3 external legs)")
    assert fj["planarity_report"]["n_external_legs"] == 2 and has(fj, NOT_DISC)
    fy = audit("row13/family_of_record.yaml")                       # inferred set
    assert has(fy, INFERRED) and fy["planar"] is None
    assert fy["planarity_report"]["withheld"].startswith(INFERRED)
    assert fy["planarity_report"]["n_external_legs"] == 2 and not has(fy, NOT_DISC)


# --------------------------------------- row 20: the kinematics contract --
@nx_required
def test_row20_with_kinematics_is_unchanged_and_non_planar_at_infinity():
    e = CENSUS["row20"]
    fp = audit("row20/family_of_record.yaml", "row20/family_kinematics.yaml")
    assert fp["planar"] is e["planar"] and dict(fp["cut_signature"]) == e["cut_signature"]
    assert fp["canonical_hash"] == e["canonical_hash"] and fp["canonical_set_size"] == 7
    assert fp["planarity_method"] == ("networkx.check_planarity with external-leg closure "
                                      "cycle in canonical order 1,2,3,4")
    assert fp["planar_legs_joined_at_infinity"] is False
    rep = fp["planarity_report"]
    assert rep["n_leg_orders"] == 3 and rep["n_planar_leg_orders"] == 0
    assert rep["planar_in_some_leg_order"] is False
    md = rep["mandelstam_declaration"]
    assert md["channels"] == {"s": [["p1", "p2"]], "t": [["p1", "p3"]]}
    assert md["cyclic_order"] == ["p1", "p2", "p4", "p3"]
    assert md["planar_in_declared_order"] is False and md["agrees_with_closure_order"] is False
    assert not has(fp, LEG_ORDER)          # both orders give the same verdict: no line


# ------------------------------------------ row 2: its own leg order --
@nx_required
def test_row02_family_planar_in_its_own_leg_order_reported_not_enforced():
    e = CENSUS["row02"]
    fam = load("row02/family_of_record.jl")
    fp = ta.audit(fam, try_leg_perms=True)
    assert fp["canonical_hash"] == e["canonical_hash"]
    assert dict(fp["cut_signature"]) == e["cut_signature"]
    assert fp["planar"] is e["planar"]                              # canonical order, as the record says
    assert fp["planarity_method"] == e["planarity_method"]          # the census string, unchanged
    assert fp["planarity_closure_order"] == ["p1", "p2", "p3", "p4"]
    assert fp["planarity_closure_order_source"] == CANONICAL_SRC
    assert fp["planar_legs_joined_at_infinity"] is True
    rep = fp["planarity_report"]
    assert rep["planar_leg_orders"] == [["p1", "p2", "p4", "p3"]]  # the Smirnov boundary order
    md = rep["mandelstam_declaration"]
    assert md["source"].startswith("AmflowFamily kinematics Dict")
    assert md["channels"] == {"s": [["p1", "p2"]], "t": [["p1", "p3"]]}
    assert md["cyclic_order"] == ["p1", "p2", "p4", "p3"]
    assert md["planar_in_declared_order"] is True and md["agrees_with_closure_order"] is False
    assert has(fp, LEG_ORDER), fp["warnings"]
    line = [w for w in fp["warnings"] if w.startswith(LEG_ORDER)][0]
    assert "['p1', 'p2', 'p4', 'p3']" in line and "planar=True" in line and "planar=False" in line
    # "the graph is planar with its legs in the drawn order": the caller's order
    pl = ta.planarity(fam, leg_order=["p1", "p2", "p4", "p3"])
    assert pl["planar"] is True and pl["closure_order_source"] == ARG_SRC
    assert pl["method"] == ("networkx.check_planarity with external-leg closure cycle in "
                            "the declared cyclic order p1,p2,p4,p3 (" + ARG_SRC + ")")
    assert ta.planarity(fam)["planar"] is False                     # and the default stays


# --------------------------------------------- theta: every leg order --
@nx_required
def test_theta_crossed_box_non_planar_in_every_enumerated_order_and_at_infinity():
    e = CENSUS["theta"]
    ref = json.load(open(fixture("row02/catalog_check.json")))["genuine_crossed_box_theta211"]
    assert ref["planar_with_legs_for_some_leg_order"] is e["planar_with_legs_for_some_leg_order"]
    assert ref["census_fingerprint"]["abstract_planar"] is e["abstract_planar"]
    for fam in (load("row02/genuine_crossed_box_theta211.yaml"),
                fam_from("theta", ta.CATALOG_SOURCES["crossed_tausk_dbox"])):
        fp = ta.audit(fam)
        assert fp["planar"] is False and dict(fp["cut_signature"]) == e["cut_signature"]
        rep = fp["planarity_report"]
        assert rep["abstract_planar"] is True
        assert rep["leg_orders_enumerated"] is True and rep["n_leg_orders"] == 3
        assert rep["n_planar_leg_orders"] == 0 and rep["planar_leg_orders"] == []
        assert rep["planar_in_some_leg_order"] is False            # the abstract shortcut said True
        assert fp["planar_legs_joined_at_infinity"] is False
    assert ta.audit(load("row02/genuine_crossed_box_theta211.yaml"))["canonical_hash"] == e["canonical_hash"]


# ------------------------- row 28: the declaration is not the drawn order --
@nx_required
def test_row28_declaration_reported_not_enforced():
    e = CENSUS["row28"]
    for fam, src in ((load("row28/family_of_record.yaml", "row28/family_of_record_kira_kinematics.yaml"),
                      "kinematics.yaml scalarproduct_rules"),
                     (load("row28/family_of_record.json"), "AMFlow-port JSON replacement block")):
        fp = ta.audit(fam)
        assert fp["canonical_hash"] == e["import_route_canonical_hash"]
        assert fp["planar"] is e["planar_with_legs_closed_in_order"]
        assert fp["planarity_closure_order"] == e["closure_order"]
        assert fp["planar_legs_joined_at_infinity"] is True
        md = fp["planarity_report"]["mandelstam_declaration"]
        assert md["source"].startswith(src), md["source"]
        assert md["channels"] == {"s": [["p1", "p2"]], "t": [["p1", "p3"]]}
        assert md["cyclic_order"] == ["p1", "p2", "p4", "p3"]
        assert md["planar_in_declared_order"] is False               # the box drawn p1,p2,p3,p4
        assert has(fp, LEG_ORDER)                                   # reported, the verdict kept
        assert fp["planarity_report"]["planar_leg_orders"] == [["p1", "p2", "p3", "p4"]]


# ------------------------------------------- row 5: a declared cyclic order --
@nx_required
def test_row05_pysecdec_leg_list_is_the_closure_order():
    e = CENSUS["row05"]
    fpsd = ta.audit(load("row05/psd_build_pentagon.py"))
    assert fpsd["canonical_hash"] == e["canonical_hash"] and fpsd["planar"] is e["planar"]
    assert fpsd["planarity_closure_order_source"] == CYCLIC_SRC
    assert fpsd["planarity_closure_order"] == ["p1", "p2", "p3", "p4", "p5"]
    assert fpsd["planarity_method"] == ("networkx.check_planarity with external-leg closure "
                                        "cycle in the declared cyclic order p1,p2,p3,p4,p5 ("
                                        + CYCLIC_SRC + ")")
    fy = audit("row05/family_of_record.yaml")
    assert fy["planarity_method"] == e["planarity_method"] and fy["planar"] is True
    assert fy["planarity_closure_order_source"] == CANONICAL_SRC
    for fp in (fpsd, fy):
        assert fp["planar_legs_joined_at_infinity"] is True
        assert fp["planarity_report"]["n_leg_orders"] == 12
        assert fp["planarity_report"]["planar_leg_orders"] == [["p1", "p2", "p3", "p4", "p5"]]


@nx_required
def test_row35_declared_channels_agree_with_the_canonical_order():
    e = CENSUS["row35"]
    fp = audit("row35/family_of_record.yaml", "row35/kinematics_of_record.yaml")
    assert fp["planar"] is e["planar_with_legs_closed_in_order"]
    assert dict(fp["cut_signature"]) == e["cut_signature"]
    assert fp["planarity_closure_order"] == e["closure_order"]
    assert fp["planar_legs_joined_at_infinity"] is False
    md = fp["planarity_report"]["mandelstam_declaration"]
    assert md["channels"] == {"s": [["p1", "p2"]], "t": [["p2", "p3"]]}
    assert md["cyclic_order"] == ["p1", "p2", "p3", "p4"] and md["agrees_with_closure_order"] is True
    assert not has(fp, LEG_ORDER)


# ------------------------------------------------------------- controls --
@nx_required
def test_control_declared_cyclic_order_flips_the_one_loop_box():
    """A drawn one-loop box: with cyclic_leg_order p1,p2,p4,p3 the closure
    reads planar False (source family.cyclic_leg_order); with p1,p2,p3,p4 or
    no declared order, planar True."""
    crossed = dict(BOX_SPEC, cyclic_leg_order=["p1", "p2", "p4", "p3"])
    fam = _spec_family(crossed, "box_crossed_order")
    assert fam.cyclic_leg_order == ["p1", "p2", "p4", "p3"]
    fp = ta.audit(fam)
    assert fp["planar"] is False and fp["planarity_closure_order_source"] == CYCLIC_SRC
    assert fp["planarity_closure_order"] == ["p1", "p2", "p4", "p3"]
    assert "declared cyclic order p1,p2,p4,p3" in fp["planarity_method"]
    assert fp["planar_legs_joined_at_infinity"] is True           # the graph itself is a box
    assert fp["planarity_report"]["planar_leg_orders"] == [["p1", "p2", "p3", "p4"]]
    fam2 = _spec_family(dict(BOX_SPEC, cyclic_leg_order=["p1", "p2", "p3", "p4"]), "box_drawn_order")
    assert ta.audit(fam2)["planar"] is True
    fam3 = _spec_family(BOX_SPEC, "box_no_order")
    fp3 = ta.audit(fam3)
    assert fp3["planar"] is True and fp3["planarity_closure_order_source"] == CANONICAL_SRC
    # a declared order that does not name the attached legs is ignored by name
    fam4 = _spec_family(dict(BOX_SPEC, cyclic_leg_order=["q1", "q2", "q3", "q4"]), "box_bad_order")
    pl4 = ta.planarity(fam4)
    assert pl4["closure_order_source"] == CANONICAL_SRC and "does not name the attached legs" in pl4["closure_order_note"]


@nx_required
def test_control_leg_order_argument_not_a_permutation_is_refused():
    fam = fam_from("box", ta.CATALOG_SOURCES["planar_smirnov_dbox"])
    for bad in (["p1", "p2", "p3"], ["p1", "p2", "p3", "p5"], ["p1", "p1", "p2", "p3"]):
        try:
            ta.planarity(fam, leg_order=bad)
        except ValueError as exc:
            assert str(exc).startswith("REFUSED: PLANARITY-LEG-ORDER"), str(exc)
        else:
            raise AssertionError("leg order %s was not refused" % bad)


@nx_required
def test_control_retired_entry_planar_in_the_crossed_order_and_the_box_not():
    """The retired legs-3-4-interchanged entry is planar in p1,p2,p4,p3 and
    the planar box is non-planar there: the leg-order argument flips both."""
    ret = fam_from("retired", ta.CATALOG_RETIRED["planar_smirnov_dbox_legs34_interchanged"])
    box = fam_from("box", ta.CATALOG_SOURCES["planar_smirnov_dbox"])
    assert ta.planarity(ret)["planar"] is False and ta.planarity(box)["planar"] is True
    assert ta.planarity(ret, leg_order=["p1", "p2", "p4", "p3"])["planar"] is True
    assert ta.planarity(box, leg_order=["p1", "p2", "p4", "p3"])["planar"] is False
    assert ta.planarity(ret)["planar_leg_orders"] == [["p1", "p2", "p4", "p3"]]
    assert ta.planarity(box)["planar_leg_orders"] == [["p1", "p2", "p3", "p4"]]
    assert ta.planarity(ret)["planar_legs_joined_at_infinity"] is True


@nx_required
def test_control_three_point_graphs_where_the_closure_test_decides():
    """Under the N < 4 flag the closure value can still be False: K3,3 (not
    planar at all) and the prism with legs on no common face (abstract planar,
    closure False, legs at infinity False)."""
    k33 = ta.audit(_spec_family(K33_SPEC, "k33"))
    assert k33["planar"] is False and k33["planarity_report"]["abstract_planar"] is False
    assert k33["planarity_discriminating"] is False and has(k33, NOT_DISC)
    assert k33["planarity_method"].startswith(FLAG + ": ")
    assert k33["planar_legs_joined_at_infinity"] is False
    pr = ta.audit(_spec_family(PRISM_SPEC, "prism"))
    assert pr["planarity_report"]["abstract_planar"] is True
    assert pr["planar"] is False and pr["planar_legs_joined_at_infinity"] is False
    assert pr["planarity_discriminating"] is False and has(pr, NOT_DISC)
    line = [w for w in pr["warnings"] if w.startswith(NOT_DISC)][0]
    assert "planar=False" in line and "abstract planar True" in line


@nx_required
def test_catalog_entries_flagged_by_leg_count_and_four_leg_strings_unchanged():
    old_form = "networkx.check_planarity with external-leg closure cycle in canonical order "
    for cname, spec in ta.CATALOG_SOURCES.items():
        fp = ta.audit(fam_from(cname, spec))
        assert fp["planar"] is spec["planar"], cname
        rep = fp["planarity_report"]
        n = len(spec["exts"])
        assert rep["n_external_legs"] == n, cname
        if n < 4:
            assert fp["planarity_method"].startswith(FLAG + ": "), cname
            assert fp["planarity_discriminating"] is False and has(fp, NOT_DISC), cname
        else:
            assert fp["planarity_method"] == old_form + ",".join(str(i + 1) for i in range(n)), cname
            assert fp["planarity_discriminating"] is True and not has(fp, NOT_DISC), cname
            assert rep["leg_orders_enumerated"] is True
            # the CAT observation, now asserted: the enumeration equals the declaration
            assert rep["planar_in_some_leg_order"] is spec["planar_in_some_leg_order"], cname
        assert fp["planar_legs_joined_at_infinity"] in (True, False), cname


def test_non_graph_family_reports_none_fields():
    fam = load("row06/retired_sevenline_vertex2L/family_of_record_transcribed.yaml")
    pl = ta.planarity(fam)
    assert pl["planar"] is None and pl["method"] == "realization-failed-or-unavailable"
    assert pl["discriminating"] is None and pl["closure_order"] is None
    assert pl["planar_legs_joined_at_infinity"] is None
    assert ta.planarity_at_infinity(fam) is None
    fp = ta.audit(fam)
    assert fp["planar"] is None and fp["planarity_discriminating"] is None
    assert not has(fp, NOT_DISC) and not has(fp, LEG_ORDER)


@nx_required
def test_planarity_at_infinity_forms_and_the_report_cache():
    fam = load("row24/retired_nonadjacent_p1p3/family_of_record.yaml",
               "row24/retired_nonadjacent_p1p3/kinematics_of_record.yaml")
    G, ok = ta.build_graph(fam)
    assert ok
    ref = CENSUS["row24_planarity_on_realized_graph"]
    assert ta.planarity_at_infinity(G, fam.exts) == ref
    assert ta.planarity_at_infinity(fam, ["p1", "p2", "p3", "p4", "p5"]) == ref
    # the report is cached on the family for the same data and recomputed when it changes
    a = ta.planarity(fam)
    b = ta.planarity(fam)
    assert a == b and a is not b
    key0 = fam._planarity_cache[0]
    fam.propagators = list(fam.propagators)
    fam.propagators[0] = (fam.propagators[0][0], fam.propagators[0][1])
    assert ta.planarity(fam) == a and fam._planarity_cache[0] == key0
    assert ta.planarity(fam, leg_order=["p1", "p2", "p3", "p5", "p4"])["closure_order"] == \
        ["p1", "p2", "p3", "p5", "p4"]
    assert ta.planarity(fam)["closure_order"] == ["p1", "p2", "p3", "p4", "p5"]


def test_cyclic_order_helpers():
    assert ta._cyclic_leg_orders(["p1", "p2", "p3"]) == [["p1", "p2", "p3"]]
    o4 = ta._cyclic_leg_orders(["p1", "p2", "p3", "p4"])
    assert len(o4) == 3 and ["p1", "p2", "p3", "p4"] in o4 and ["p1", "p2", "p4", "p3"] in o4
    assert len(ta._cyclic_leg_orders(list("abcde"))) == 12
    assert ta._same_cyclic_order(["p1", "p2", "p3", "p4"], ["p2", "p3", "p4", "p1"])
    assert ta._same_cyclic_order(["p1", "p2", "p3", "p4"], ["p4", "p3", "p2", "p1"])
    assert not ta._same_cyclic_order(["p1", "p2", "p3", "p4"], ["p1", "p2", "p4", "p3"])


# ---------------------------------------------------------------------------
#  battery runner (used by the tool's --self-test leg [planarity] and by __main__)
# ---------------------------------------------------------------------------
BATTERY = [
    test_fixture_provenance_pins,
    test_row24_retired_reproduces_the_census_planarity_closure_dict,
    test_row24_cured_json_reproduces_the_census_planarity_closure_dict,
    test_rows15_16_17_with_kinematics_are_flagged_not_discriminating,
    test_row23_with_kinematics_two_legs_flagged,
    test_row18_fam2_three_legs_inferred_is_withheld_with_the_report,
    test_rows08_13_two_point_json_legs_flagged_yaml_inferred_withheld,
    test_row20_with_kinematics_is_unchanged_and_non_planar_at_infinity,
    test_row02_family_planar_in_its_own_leg_order_reported_not_enforced,
    test_theta_crossed_box_non_planar_in_every_enumerated_order_and_at_infinity,
    test_row28_declaration_reported_not_enforced,
    test_row05_pysecdec_leg_list_is_the_closure_order,
    test_row35_declared_channels_agree_with_the_canonical_order,
    test_control_declared_cyclic_order_flips_the_one_loop_box,
    test_control_leg_order_argument_not_a_permutation_is_refused,
    test_control_retired_entry_planar_in_the_crossed_order_and_the_box_not,
    test_control_three_point_graphs_where_the_closure_test_decides,
    test_catalog_entries_flagged_by_leg_count_and_four_leg_strings_unchanged,
    test_non_graph_family_reports_none_fields,
    test_planarity_at_infinity_forms_and_the_report_cache,
    test_cyclic_order_helpers,
]


def run_battery(verbose=True):
    rows = []
    for fn in BATTERY:
        name = fn.__name__
        try:
            fn()
            status, detail = "PASS", ""
        except FixtureMissing as exc:
            status, detail = "SKIP", str(exc)
        except FixtureDrift as exc:
            status, detail = "FAIL", "refused: " + str(exc)
        except AssertionError as exc:
            status, detail = "FAIL", (str(exc) or "assertion failed")[:300]
        except Exception as exc:
            status, detail = "FAIL", f"{type(exc).__name__}: {exc}"[:300]
        rows.append({"name": name, "status": status, "detail": detail})
        if verbose:
            print(f"   [planarity] {name}: {status}" + (f"  ({detail})" if detail else ""), flush=True)
    return rows


if __name__ == "__main__":
    out = run_battery(verbose=True)
    n_fail = sum(1 for r in out if r["status"] == "FAIL")
    n_skip = sum(1 for r in out if r["status"] == "SKIP")
    print(f"planarity battery: {len(out)} cases, {n_fail} FAIL, {n_skip} SKIP")
    sys.exit(1 if n_fail else (2 if n_skip else 0))
