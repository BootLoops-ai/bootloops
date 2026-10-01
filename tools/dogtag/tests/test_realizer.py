# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision.
"""
Battery of the vertex realizer (build_graph / _realize_graph /
_exact_vertex_realizer / connected_covers / non_graph_test / cut_signature).

FIXTURES = TESTS.  Every fixture under tests/fixtures/<rowNN>/ is a copy of an
identity-census corpus file (the original's sha256 kept beside the copy's),
pinned by sha256 in that directory's PROVENANCE.json and REFUSED on drift; every
expected value below is
copied from the census's own row record (ROW.md / iso.json / lint_family.json /
lint.json), never paraphrased.  Each positive row is paired with its planted
negative control.  Where a fixture's leg set is not spelled by the propagators
(the yaml transcriptions of JSON configs, the three-point and five-point
families), the legs and the momentum-conservation rule are set on the loaded
Family exactly as the census set them on its import route ("four legs set
explicitly"); the readers that carry them from the file belong to other waves.

Run:  python3 -m pytest tests/test_realizer.py -q     (from tools/dogtag)
"""
import hashlib
import json
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import topology_audit as TA  # noqa: E402

# a case that needs networkx (a declared dependency of the tool) SKIPs by name
# when the dependency is absent -- never PASS, never FAIL
nx_required = pytest.mark.skipif(not TA._HAVE_NX, reason=TA.NX_SKIP_REASON)

FIX = os.path.join(HERE, "fixtures")

FOUR = (["p1", "p2", "p3", "p4"], {"p4": "-p1 - p2 - p3"})
THREE = (["p1", "p2", "p3"], {"p3": "-p1 - p2"})
FIVE = (["p1", "p2", "p3", "p4", "p5"], {"p5": "-p1 - p2 - p3 - p4"})


def pinned(row, rel, root=None):
    """Path of a vendored fixture; refused (test failure by name) on sha256 drift
    against the row's PROVENANCE.json, and on an absent file or pin."""
    root = FIX if root is None else root
    path = os.path.join(root, row, rel)
    prov = os.path.join(root, row, "PROVENANCE.json")
    assert os.path.exists(prov), "PROVENANCE.json absent for %s" % row
    assert os.path.exists(path), "fixture absent: %s/%s" % (row, rel)
    want = None
    for ent in json.load(open(prov))["files"]:
        if ent["path"] == rel:
            want = ent["sha256"]
    have = hashlib.sha256(open(path, "rb").read()).hexdigest()
    assert want is not None, "no sha256 pin for %s/%s" % (row, rel)
    assert have == want, "fixture drift %s/%s: PROVENANCE %s != file %s -- refused" % (
        row, rel, want, have)
    return path


def family(row, rel, legs=None):
    fam = TA.load_family(pinned(row, rel))[0]
    if legs is not None:
        fam.exts, fam.ext_subs = list(legs[0]), dict(legs[1])
    return fam


def realized(fam):
    G, ok = TA.build_graph(fam)
    return ok, G.graph["realization"]


# ----------------------------------------------------------------------------
# fixture pins: every vendored file matches its PROVENANCE sha256 and size
# ----------------------------------------------------------------------------
def test_every_fixture_pinned_and_unchanged():
    rows = sorted(d for d in os.listdir(FIX) if os.path.isdir(os.path.join(FIX, d)))
    assert rows, "no fixtures vendored"
    n = 0
    for row in rows:
        prov = json.load(open(os.path.join(FIX, row, "PROVENANCE.json")))
        for ent in prov["files"]:
            p = os.path.join(FIX, row, ent["path"])
            b = open(p, "rb").read()
            assert hashlib.sha256(b).hexdigest() == ent["sha256"], "drift: %s/%s" % (row, ent["path"])
            assert len(b) == ent["size"], "size drift: %s/%s" % (row, ent["path"])
            n += 1
    assert n >= 18


def test_fixture_drift_is_refused(tmp_path):
    """Planted control of the pin itself: a tampered copy must be refused."""
    src = pinned("row01", "drawn_graph.yaml")
    row = tmp_path / "row01"
    row.mkdir()
    (row / "drawn_graph.yaml").write_bytes(open(src, "rb").read() + b"\n# tampered\n")
    prov = json.load(open(os.path.join(FIX, "row01", "PROVENANCE.json")))
    (row / "PROVENANCE.json").write_text(json.dumps(prov))
    with pytest.raises(AssertionError, match="fixture drift"):
        pinned("row01", "drawn_graph.yaml", root=str(tmp_path))


# ----------------------------------------------------------------------------
# item 5: two-particle-reducible insertions (rows 28, 29, 33) realize CONNECTED
# census iso.json record_connected_realization_fingerprint:
#   row 28 {"connected": true, "cut_signature": {"s": 3, "t": 2, "u": 5}, "V": 6, "E_distinct_lines": 8, "N": 4}
#   row 29 {"cut_signature": {"s": 2, "t": 2, "u": 4}, "V": 7, "E_distinct_lines": 9, "N": 4}
#   row 33 {"cut_signature": {"s": 2, "t": 2, "u": 4}, "V": 6, "E_distinct_lines": 8, "N": 4}
# and "three_V": 18, "two_E_plus_N": 20, "nongraph_test_3V_le_2E_plus_N": true (row 28),
# 21 / 22 (row 29); the tool's cut signature before: row 28 {"s": null, "t": 2, "u": null}.
# ----------------------------------------------------------------------------
@pytest.mark.parametrize("row,rel,V,E,three_V,two_E_plus_N,cuts", [
    ("row28", "family_of_record.yaml", 6, 8, 18, 20, {"s": 3, "t": 2, "u": 5}),
    ("row29", "family_of_record.yaml", 7, 9, 21, 22, {"s": 2, "t": 2, "u": 4}),
    ("row33", "family_of_record.yaml", 6, 8, 18, 20, {"s": 2, "t": 2, "u": 4}),
])
@nx_required
def test_reducible_insertion_realizes_connected(row, rel, V, E, three_V, two_E_plus_N, cuts):
    fam = family(row, rel, FOUR)
    ok, rep = realized(fam)
    assert ok and rep["verdict"] == "GRAPH"
    assert rep["connected"] is True
    assert (rep["V"], rep["E"], rep["N"]) == (V, E, 4)
    assert (rep["three_V"], rep["two_E_plus_N"]) == (three_V, two_E_plus_N)
    assert rep["count_test_3V_le_2E_plus_N"] is True
    assert rep["conservation_at_every_vertex"] is True
    assert rep["legs_per_vertex"][:4] == [1, 1, 1, 1]
    csig = TA.cut_signature(fam)
    assert dict(csig) == cuts
    assert csig.cause is None
    assert None not in csig.values()
    # every distinct connected realization gives the same cut signature
    assert csig.ambiguity is None
    # census: "planar with legs closed in order p1,p2,p3,p4 = True (connected realization)"
    assert TA.planarity(fam)["planar"] is True


@nx_required
def test_row28_tool_realization_was_disconnected_before():
    """The class of the flip: the record's own cover search (first zero-sum
    partition, no connectivity) closed the kite on itself.  With the connected
    requirement the rim line and the box share vertices: the realized graph is
    connected and every leg attaches to a distinct vertex."""
    fam = family("row28", "family_of_record.yaml", FOUR)
    G, ok = TA.build_graph(fam)
    assert ok
    import networkx as nx
    assert nx.is_connected(G)
    attach = [next(iter(G.neighbors("ext_%s" % e))) for e in fam.exts]
    assert len(set(attach)) == 4


@nx_required
def test_row28_negative_rim_line_dropped_is_the_contracted_subsector():
    """The control as the plan words it (item 5): the kite rim line
    "(l - k1)^2" dropped from the top sector "(7 lines, E-V+1 no longer 3 with
    4 legs) must report NON-GRAPH".  RESULT, not a flip: dropping one line is
    the contraction of that line, and the contraction of a connected graph is
    a connected graph with V - 1 vertices and the same loop count (E = 7,
    V = 5, L = 7 - 5 + 1 = 3; 3V = 15 <= 2E + N = 18).  The object is the
    kite with one rim line contracted -- a triangle with a doubled edge
    inserted on the box rung -- and the realizer correctly finds it.  No
    single-line drop can produce a non-graph; the momentum-breaking control
    below is the one that does."""
    fam = family("row28", "family_of_record.yaml", FOUR)
    i = [k for k, (expr, _) in enumerate(fam.propagators) if expr.strip() == "(l - k1)^2"]
    assert len(i) == 1
    fam.physical[i[0]] = False
    ok, rep = realized(fam)
    assert ok and rep["verdict"] == "GRAPH"
    assert (rep["E"], rep["L"], rep["V"]) == (7, 3, 5)
    assert (rep["three_V"], rep["two_E_plus_N"]) == (15, 18)
    assert rep["connected"] is True and rep["conservation_at_every_vertex"] is True
    assert sorted(rep["vertex_degrees"], reverse=True) == [5, 4, 3, 3, 3]


@nx_required
def test_row28_negative_rim_line_broken_is_refused():
    """Planted control that does exercise "no connected cover with
    V = E - L + 1" on row 28: the kite rim line "(l - k1)^2" re-written to
    "(l - k1 + p2)^2" (8 lines, 3 loops, 4 legs: the count test passes, the
    loop rank is 3, but no zero-sum vertex can absorb the shifted line) ->
    NON-GRAPH by name, never a planar verdict."""
    fam = family("row28", "family_of_record.yaml", FOUR)
    i = [k for k, (expr, _) in enumerate(fam.propagators) if expr.strip() == "(l - k1)^2"]
    assert len(i) == 1
    fam.propagators = list(fam.propagators)
    fam.propagators[i[0]] = ("(l - k1 + p2)^2", fam.propagators[i[0]][1])
    ok, rep = realized(fam)
    assert not ok
    assert rep["verdict"] == "NON-GRAPH"
    assert (rep["E"], rep["V"]) == (8, 6)
    assert rep["count_test_3V_le_2E_plus_N"] is True and rep["loop_rank_test"] is True
    assert "no exact cover" in rep["cause"] or "none connected" in rep["cause"]
    assert TA.planarity(fam)["planar"] is None
    csig = TA.cut_signature(fam)
    assert all(v is None for v in csig.values())
    assert csig.cause.startswith("NON-GRAPH")


# ----------------------------------------------------------------------------
# item 12: the non-graph verdict.  Row 6 retired, ROW.md: "E=7 lines, L=2, V=6,
# 3V=18 > 2E+N=17; exact covers of the half-edges found: 20"; census verdict
# "NON-GRAPH".  Row 6 cured, ROW.md: "E=6 lines, L=2, V=5, 3V=15 <= 2E+N=15;
# momentum conservation at every realized vertex: True", fingerprint (import)
# "canonical hash 249f94470cc4, planar True, cut signature {'p1p2': 2, 'p1p3': 2, 'p2p3': 2}".
# ----------------------------------------------------------------------------
@nx_required
def test_row06_retired_sevenline_is_non_graph():
    fam = family("row06", "retired_sevenline_vertex2L/family_of_record_transcribed.yaml", THREE)
    ok, rep = realized(fam)
    assert not ok
    assert rep["verdict"] == "NON-GRAPH"
    assert (rep["E"], rep["L"], rep["V"]) == (7, 2, 6)
    assert (rep["three_V"], rep["two_E_plus_N"]) == (18, 17)
    assert rep["count_test_3V_le_2E_plus_N"] is False
    assert rep["n_exact_covers_found"] == 20
    assert rep["n_covers_with_V_eq_E_minus_L_plus_1"] == 0
    assert "3V = 18 > 2E + N = 17" in rep["cause"]
    csig = TA.cut_signature(fam)
    assert dict(csig) == {"p1p2": None, "p1p3": None, "p2p3": None}
    assert "3V = 18 > 2E + N = 17" in csig.cause
    assert TA.planarity(fam)["planar"] is None
    # the tool before: "planar": true, "cut_signature": {"p1p2": null, "p1p3": 2, "p2p3": null}


@nx_required
def test_row06_cured_ladder_is_a_graph():
    fam = family("row06", "cured_ladder_sudakovPR0/family_of_record.yaml", THREE)
    ok, rep = realized(fam)
    assert ok and rep["verdict"] == "GRAPH"
    assert (rep["E"], rep["L"], rep["V"]) == (6, 2, 5)
    assert (rep["three_V"], rep["two_E_plus_N"]) == (15, 15)
    assert rep["conservation_at_every_vertex"] is True
    assert rep["realization_unique"] is True
    assert rep["n_distinct_connected_realizations"] == 1
    assert rep["legs_per_vertex"] == [1, 1, 1, 0, 0]
    assert dict(TA.cut_signature(fam)) == {"p1p2": 2, "p1p3": 2, "p2p3": 2}
    assert TA.planarity(fam)["planar"] is True
    assert TA.fingerprint(fam)["canonical_hash"] == "249f94470cc4"


@nx_required
def test_row06_cured_negative_top_sector_127_is_non_graph():
    """Planted control (plan item 12): the same yaml with top_level_sectors
    127 (the ISP promoted to a seventh line) -> NON-GRAPH."""
    fam = family("row06", "cured_ladder_sudakovPR0/family_of_record.yaml", THREE)
    fam.physical = [True] * len(fam.propagators)
    ok, rep = realized(fam)
    assert not ok and rep["verdict"] == "NON-GRAPH"
    assert (rep["three_V"], rep["two_E_plus_N"]) == (18, 17)


# ----------------------------------------------------------------------------
# item 25: the drawn double boxes (rows 1, 2, 7) realize with V = 6, not 5.
# census lint.json drawn_graph: "checker_realizer": {"ok": true, "V": 5} (the
# defect), "V": 6, "realizations": {"n_partitions_with_V_eq_E_minus_L_plus_1": 1,
# "n_distinct_connected_graphs": 1, "unique": true}, "census_fingerprint":
# {"abstract_planar": true, "planar_with_legs_in_canonical_order": true,
# "cut_signature": {"s": 2, "t": 3, "u": 4}}; ROW.md "realization unique=True".
# ----------------------------------------------------------------------------
@nx_required
@pytest.mark.parametrize("row", ["row01", "row02", "row07"])
def test_drawn_double_box_realizes_with_six_vertices(row):
    fam = family(row, "drawn_graph.yaml")
    assert fam.exts == ["p1", "p2", "p3", "p4"]
    ok, rep = realized(fam)
    assert ok
    assert rep["V"] == 6
    assert (rep["E"], rep["L"], rep["N"]) == (7, 2, 4)
    assert (rep["three_V"], rep["two_E_plus_N"]) == (18, 18)
    assert rep["n_covers_with_V_eq_E_minus_L_plus_1"] == 1
    assert rep["n_distinct_connected_realizations"] == 1
    assert rep["realization_unique"] is True
    assert rep["vertex_degrees"] == [3, 3, 3, 3, 3, 3]
    assert rep["conservation_at_every_vertex"] is True
    csig = TA.cut_signature(fam)
    assert dict(csig) == {"s": 2, "t": 3, "u": 4}
    plan = TA.planarity(fam)
    assert plan["abstract_planar"] is True
    assert plan["planar"] is True
    # the tool before (pilot): "cut_signature": {"s": null, "t": 4, "u": null}


@nx_required
def test_drawn_double_box_negative_duplicated_ring_line_is_refused():
    """Planted control (plan item 25): one ring line duplicated in the top
    sector (E = 8 with two loop symbols) -> refused by name."""
    fam = family("row01", "drawn_graph.yaml")
    fam.propagators = list(fam.propagators) + [fam.propagators[1]]
    fam.physical = [True] * len(fam.propagators)
    ok, rep = realized(fam)
    assert not ok and rep["verdict"] == "NON-GRAPH"
    assert rep["E"] == 8
    assert ("3V = 21 > 2E + N = 20" in rep["cause"]) or ("repeated line" in rep["cause"])
    assert TA.planarity(fam)["planar"] is None


@nx_required
def test_repeated_line_named_when_the_count_test_passes():
    """The distinct-lines test by itself: a bubble whose second line repeats
    the first (E = 2, L = 1, N = 2: 3V = 6 <= 2E + N = 6 passes) is refused
    as a repeated line, not silently realized."""
    fam = TA.Family("rep", ["k1"], ["p1", "p2"], {"p2": "-p1"},
                    [("k1^2", 0), ("(-k1)^2", 0)], "memory:control")
    ok, rep = realized(fam)
    assert not ok and rep["verdict"] == "NON-GRAPH"
    assert rep["count_test_3V_le_2E_plus_N"] is True
    assert rep["distinct_lines_test"] is False
    assert "repeated line" in rep["cause"]


# ----------------------------------------------------------------------------
# item 25: row 27.  Retired LBL3E, lint_family.json: "n_exact_covers_with_V_eq_E-L+1": 3,
# "n_connected_covers": 2, "n_distinct_connected_realizations": 1, "realization_unique": true,
# "legs_per_vertex_by_realization": [[2, 1, 1, 0], [2, 1, 1, 0]]; ROW.md "photons per
# vertex [2, 1, 1, 0]".  Cured LBL3Ebox: "E=6 massive lines, L=3, V=4, 3V=12 <= 2E+N=16;
# conservation at every vertex True; photons per vertex [1, 1, 1, 1]" and "Fingerprint
# (import): canonical hash a334bce38e65, cut signature {'s': 2, 't': 4, 'u': 6}, planar True".
# ----------------------------------------------------------------------------
@nx_required
def test_row27_retired_reports_three_covers_and_two_photons_on_one_vertex():
    fam = family("row27", "retired_threepoint_LBL3E/family_of_record.yaml", FOUR)
    ok, rep = realized(fam)
    assert ok
    assert (rep["E"], rep["L"], rep["V"]) == (6, 3, 4)
    assert (rep["three_V"], rep["two_E_plus_N"]) == (12, 16)
    assert rep["n_covers_with_V_eq_E_minus_L_plus_1"] == 3
    assert rep["n_connected_covers"] == 2
    assert rep["n_distinct_connected_realizations"] == 1
    assert rep["realization_unique"] is True
    assert rep["legs_per_vertex_by_realization"] == [[2, 1, 1, 0], [2, 1, 1, 0]]
    assert rep["legs_per_vertex"] == [2, 1, 1, 0]
    assert rep["conservation_at_every_vertex"] is True
    csig = TA.cut_signature(fam)
    # two photons on one vertex: the s and u channels have no internal cut;
    # the None entries carry the cause by name (the tool before: no cause)
    assert csig["t"] == 2 and csig["s"] is None and csig["u"] is None
    assert "attach to one vertex" in csig.cause
    assert "legs per vertex [2, 1, 1, 0]" in csig.cause
    assert csig.ambiguity is None


@nx_required
def test_row27_cured_box_realizes_with_one_photon_per_vertex():
    fam = family("row27", "cured_drawn_box_LBL3Ebox/family_of_record.yaml", FOUR)
    ok, rep = realized(fam)
    assert ok
    assert (rep["E"], rep["L"], rep["V"]) == (6, 3, 4)
    assert (rep["three_V"], rep["two_E_plus_N"]) == (12, 16)
    assert rep["conservation_at_every_vertex"] is True
    assert rep["legs_per_vertex"] == [1, 1, 1, 1]
    assert rep["realization_unique"] is True
    csig = TA.cut_signature(fam)
    assert dict(csig) == {"s": 2, "t": 4, "u": 6}
    assert csig.cause is None
    assert TA.planarity(fam)["planar"] is True
    assert TA.fingerprint(fam)["canonical_hash"] == "a334bce38e65"


# ----------------------------------------------------------------------------
# item 25: row 24 (five legs).  lint_family.json (both families): "V_expected_E_minus_L_plus_1": 7,
# "count_test_3V_le_2E_plus_N": {"3V": 21, "2E+N": 21, "passes": true}, "n_exact_covers_with_V_eq_E-L+1": 1,
# "n_connected_covers": 1, "realization_unique": true, "legs_per_vertex_by_realization": [[1, 1, 1, 1, 1, 0, 0]],
# "momentum_conservation_at_every_vertex": true; audit_family_import.json (retired):
# "planarity_on_realized_graph": {"planar_canonical_leg_order": false, "abstract_planar": true};
# the tool before (import route): cut_signature all null, "planar": true (the wrong cover).
# ----------------------------------------------------------------------------
@pytest.mark.parametrize("rel", ["retired_nonadjacent_p1p3/family_of_record.yaml",
                                 "cured_adjacent_p3p4/family_of_record_transcribed.yaml"])
@nx_required
def test_row24_five_point_realizes_unique_with_no_null_cuts(rel):
    fam = family("row24", rel, FIVE)
    ok, rep = realized(fam)
    assert ok
    assert (rep["E"], rep["L"], rep["N"], rep["V"]) == (8, 2, 5, 7)
    assert (rep["three_V"], rep["two_E_plus_N"]) == (21, 21)
    assert rep["n_covers_with_V_eq_E_minus_L_plus_1"] == 1
    assert rep["n_connected_covers"] == 1
    assert rep["realization_unique"] is True
    assert rep["legs_per_vertex"] == [1, 1, 1, 1, 1, 0, 0]
    assert rep["conservation_at_every_vertex"] is True
    csig = TA.cut_signature(fam)
    assert len(csig) == 10 and None not in csig.values()
    assert csig.cause is None
    plan = TA.planarity(fam)
    assert plan["abstract_planar"] is True
    assert plan["planar"] is False


# ----------------------------------------------------------------------------
# item 25: row 21 drawn K4 (four four-valent corners).  iso.json
# record_connected_realization_fingerprint: {"connected": true, "abstract_planar": true,
# "planar_with_legs_closed_in_order": false, "cut_signature": {"s": 4, "t": 4, "u": 4},
# "V": 4, "E_distinct_lines": 6, "N": 4}; the tool before on the drawing: "planar": null.
# ----------------------------------------------------------------------------
@nx_required
def test_row21_drawn_k4_realizes():
    fam = family("row21", "drawn_graph.yaml", FOUR)
    ok, rep = realized(fam)
    assert ok
    assert (rep["V"], rep["E"], rep["N"]) == (4, 6, 4)
    assert rep["vertex_degrees"] == [4, 4, 4, 4]
    assert rep["realization_unique"] is True
    assert dict(TA.cut_signature(fam)) == {"s": 4, "t": 4, "u": 4}
    plan = TA.planarity(fam)
    assert plan["abstract_planar"] is True
    assert plan["planar"] is False


@nx_required
def test_row21_drawn_without_its_fourth_leg_is_not_checkable_not_non_graph():
    """The drawing's propagators spell p1, p2, p4 only; read as-is the leg set
    cannot close.  That is a leg-set defect of the input, reported by name as
    NOT-CHECKABLE (legs do not sum to zero) -- never NON-GRAPH."""
    # the reader's own object (the fallback leg set, before the kinematics
    # wave's completion): three legs that cannot close -> NOT-CHECKABLE
    fam = TA.load_kira_yaml(pinned("row21", "drawn_graph.yaml"))[0]
    assert fam.exts == ["p1", "p2", "p4"] and fam.ext_subs == {}
    ok, rep = realized(fam)
    assert not ok
    assert rep["verdict"] == "NOT-CHECKABLE"
    assert rep["legs_sum_to_zero_test"] is False
    assert "do not sum to zero" in rep["cause"]
    assert rep["steps"] == 0
    assert TA.cut_signature(fam).cause.startswith("NOT-CHECKABLE")
    # load_family completes the unspelled leg (p3 = -p1 - p2 - p4, the gap in
    # the numbering) and marks the set LEG-SET-INFERRED: the graph realizes,
    # the planarity verdict is withheld by audit()
    fam = family("row21", "drawn_graph.yaml")
    assert fam.exts == ["p1", "p2", "p4", "p3"] and fam.ext_subs == {"p3": "-p1 - p2 - p4"}
    assert fam.kinematics_source == "inferred"
    ok, rep = realized(fam)
    assert ok and rep["connected"] is True


# ----------------------------------------------------------------------------
# budget, shape and unchanged-verdict controls
# ----------------------------------------------------------------------------
@nx_required
def test_step_budget_exhaustion_is_named_and_gives_no_verdict():
    fam = family("row01", "drawn_graph.yaml")
    G, ok = TA._realize_graph(fam, step_budget=1)
    rep = G.graph["realization"]
    assert not ok
    assert rep["budget_exhausted"] is True
    assert rep["verdict"] == "UNDECIDED"
    assert "step budget" in rep["cause"]
    assert TA.REALIZER_STEP_BUDGET == 400000


@nx_required
def test_cut_signature_is_a_plain_mapping_for_json():
    fam = family("row01", "drawn_graph.yaml")
    csig = TA.cut_signature(fam)
    assert isinstance(csig, dict)
    assert json.loads(json.dumps(csig)) == {"s": 2, "t": 3, "u": 4}
    assert csig.realization["verdict"] == "GRAPH"
    fp = TA.fingerprint(fam)
    assert fp["cut_signature_multiset"] == [2, 3, 4]


@pytest.mark.parametrize("name,V,planar,cuts", [
    # V, planar and cuts as the tool BEFORE this realizer (the base copy)
    # emitted them for its own catalog entries; they must not move.
    ("planar_smirnov_dbox", 6, True, {"s": 2, "t": 3, "u": 4}),
    # crossed_tausk_dbox is the (2,1,1) theta graph in the catalog unit
    # (CATALOG_SOURCES; the _selftest_catalog leg asserts it): cut signature
    # {s: 2, t: 3, u: 3} (census catalog_check.json); the former list lives
    # in CATALOG_RETIRED and keeps the base values.
    ("crossed_tausk_dbox", 6, False, {"s": 2, "t": 3, "u": 3}),
    ("RETIRED:planar_smirnov_dbox_legs34_interchanged", 6, False, {"s": 2, "t": 4, "u": 3}),
    ("massless_box_1l", 4, True, {"s": 2, "t": 2, "u": 4}),
    ("massless_triangle_1l", 3, True, {"p1p2": 2, "p1p3": 2, "p2p3": 2}),
    ("massless_bubble_1l", 2, True, {"p1p2": 0}),
    ("sunrise_2l", 2, True, {"p1p2": 0}),
    ("ladder_vertex_2l", 5, True, {"p1p2": 2, "p1p3": 2, "p2p3": 2}),
    ("planar_triple_box_3l", 8, True, {"s": 2, "t": 4, "u": 4}),
])
@nx_required
def test_catalog_entries_unchanged(name, V, planar, cuts):
    """No row where the tool agreed with the hand check may change: every
    catalog entry keeps its planarity verdict and cut signature and now
    realizes with V = E - L + 1, connected, unique."""
    if name.startswith("RETIRED:"):
        spec = TA.CATALOG_RETIRED[name[len("RETIRED:"):]]
    else:
        spec = TA.CATALOG_SOURCES[name]
    fam = TA.Family(name, spec["loops"], spec["exts"], spec["ext_subs"],
                    spec["propagators"], "memory:catalog")
    ok, rep = realized(fam)
    assert ok and rep["V"] == V and rep["connected"] is True
    assert rep["realization_unique"] is True
    assert TA.planarity(fam)["planar"] is planar
    assert dict(TA.cut_signature(fam)) == cuts


@nx_required
def test_selftest_leg_registered_and_green():
    legs = dict((lid, fn) for lid, _, fn in TA.EXTRA_SELF_TEST_LEGS)
    assert "realizer" in legs
    out = legs["realizer"]()
    assert out["ok"] is True
    assert out["skipped"] == []
    assert out["executed"] >= 8
