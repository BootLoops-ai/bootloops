#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision.
"""
Leg-permutation battery for topology_audit.py (the --leg-perms search):
momentum conservation re-imposed after every leg permutation
(_dependent_leg_rows / _reimpose_conservation), the leg-count guard in
isomorphism_to and audit (catalog entries with a different number of external
legs are skipped by name, never indexed past the vector), and the
label-preserving restriction of the leg permutations.

Every expected value is the identity census's hand-checked object for that
row, quoted from the census record beside the vendored fixture
(tests/fixtures/<row>/..., PROVENANCE.json pinning the vendored and the
record sha256; refused on drift, skipped by name when absent):

  retired catalog entry  (rows/row02/catalog_check.json, CATALOG_RETIRED)
      == planar_smirnov_dbox under p3 <-> p4 by the tool's own search
  row 21  lbl3x  (rows/row21/iso.json "labeled_isomorphism_drawn_vs_record")
      relabeling "k1 -> -k1, k3 -> -k3; p2 -> p4, p4 -> p2",
      leg_perm {p1: p1, p2: p4, p3: p3, p4: p2}, n_relabelings_tried 126,
      n_leg_perms_allowed_by_labels 24, leg_perm_cyclic_order_preserved true
      controls: 'k - p4 - b' -> 'k - p3 - b' (a K4 with the legs re-seated:
               every match breaks the cyclic order); 'k - p4 - b' -> 'k - p4 + b'
               (not a graph: no match); p1 <-> p2 exchanged in every string ->
               the composed map
  row 5   pentagon  (rows/row05/iso.json "family_of_record_vs_drawn")
      momentum_isomorphic true, n_momentum_relabelings_found 2,
      first_momentum_relabeling "identity"; the second is the p5-moving
      reflection; without re-imposition only identity (the former search)
      control: 'k1 + p1 + p2' -> 'k1 + p1 + p3' -> no match
  rows 20 / 6 retired / Smirnov 2box as printed  (families whose leg count
      differs from the four-leg catalog entries: a naive catalog search would
      index past the leg vector)  -> rc 0, catalog_matches [],
      LEG-COUNT SKIP warning naming the 4-leg entries
  rows 15 / 16  (rows/row15,16/iso.json "pairs"[0]): tool_style "identity" /
      "k1 -> -k1, k2 -> -k2", n_relabelings_tested 48
"""
import collections
import hashlib
import itertools
import json
import os
import subprocess
import sys

import pytest
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
sys.path.insert(0, PKG)

import topology_audit as T  # noqa: E402

# a case that needs networkx (a declared dependency of the tool) SKIPs by name
# when the dependency is absent -- never PASS, never FAIL
nx_required = pytest.mark.skipif(not T._HAVE_NX, reason=T.NX_SKIP_REASON)

FIXTURES = os.path.join(HERE, "fixtures")
PINS = {
    "row21/family_of_record.json":
        "209d2d1741147557a6b7529415606084e3a57660beb6b57b89b1c8bec1dfe42b",
    "row21/drawn_graph.yaml":
        "c307f8008add78e64099b7fb9e7948f91a81c5646bb75d393d9a47a7efea6acf",
    "row05/family_of_record.yaml":
        "0ec539a5fec61bdf5f3c0e5c59f9cda7f8044469200559ec21aadd32d8a8bccc",
    "row05/drawn_graph.yaml":
        "33984b814f95f701aef84f7965acbdf029247cda5e06b73a6ec2343f92e6599d",
    "row20/family_of_record.yaml":
        "d4aef7e8609a1f653af4fa6ba04f1d870c50be86d60ce87a20f1690d4bb8a993",
    "row20/drawn_graph.yaml":
        "c756f4b316e8953f8736a4c3b81cb88b4d12a4d6603fed6da2ecb65a2cebfae9",
    "row06/retired_sevenline_vertex2L/family_of_record_transcribed.yaml":
        "e759606f59374bd60bf8593867ac650a44990107466b46802f550fcb5c1649d5",
    "row02/smirnov_2box_as_printed.yaml":
        "93717fc0180fd92bd72e8956fcbfe745d527afec60073d171da8d347e200d522",
    "row02/family_of_record.jl":
        "ddc75ae24954b08103d01e996f614cae6a6b74a8bbc485c80711cc7e17434917",
    "row15/family_of_record_kira_icc_eqmass.yaml":
        "9d0e526a983173045d5571fd4dd05148c306ca176b414d79353ead17054b82e2",
    "row15/drawn_graph.yaml":
        "2673aab9b205d47eaf2378bb50385b8b6b17cce690f4f4a1925f7a8c580aac34",
    "row16/family_of_record_kira_iccg_generic.yaml":
        "5be2fe7689e8dee6aab40bd7ab9d21f1efa2e7655d80333f0d87c86dfca1e581",
    "row16/drawn_graph.yaml":
        "0d09c246ad19868783db9bb3ff8eeae721aea2f48abcf79da3fb500fed806139",
}

# census objects, quoted
R21_RELABEL = "k1 -> -k1, k3 -> -k3; p2 -> p4, p4 -> p2"
R21_LEG_PERM = {"p1": "p1", "p2": "p4", "p3": "p3", "p4": "p2"}
R21_TRIED = 126
R21_LEG_PERMS = 24
R21_LEGS = (["p1", "p2", "p3", "p4"], {"p4": "-p1-p2-p3"})     # rows/row21/drawn_graph.yaml header
R05_N_FOUND = 2
R05_FIRST = "identity"
R05_P5_MAP = "k1 -> -k1; p1 -> p5, p2 -> p4, p4 -> p2, p5 -> p1"
R15_RELABEL = "identity"
R16_RELABEL = "k1 -> -k1, k2 -> -k2"
R1516_TESTED = 48
FOUR_LEG_ENTRIES = ["planar_smirnov_dbox", "crossed_tausk_dbox"]
# (fixture, legs the tool reads after leg completion, LEG-COUNT SKIP expected):
#   row 20 family + drawn: p1..p3 spelled, the four-leg entries are skipped by name;
#   row 6 retired: p1, p2 spelled and the implicit p3 completed from momentum
#   conservation (LEG-SET-COMPLETED), so the family reads 3 legs and the four-leg
#   entries are skipped by name;
#   Smirnov 2box as printed: the yaml's own external_momenta / momentum_conservation
#   keys give four legs as DATA, so the four-leg entries are SEARCHED (no match: the
#   object is not a graph) and no entry of the same line count is skipped -- the
#   pre-fix IndexError is still the flip (rc 1 -> 0).
RC1_SPECIMENS = [("row20/family_of_record.yaml", 3, True), ("row20/drawn_graph.yaml", 3, True),
                 ("row06/retired_sevenline_vertex2L/family_of_record_transcribed.yaml", 3, True),
                 ("row02/smirnov_2box_as_printed.yaml", 4, False)]
TOOL = os.path.join(PKG, "topology_audit.py")


def fixture(rel):
    """Path of a pinned fixture; refuses on drift, skips by name when absent."""
    path = os.path.join(FIXTURES, rel)
    if not os.path.exists(path):
        pytest.skip(f"fixture tests/fixtures/{rel} absent")
    got = hashlib.sha256(open(path, "rb").read()).hexdigest()
    if got != PINS[rel]:
        pytest.fail(f"fixture tests/fixtures/{rel} drifted: sha256 {got[:16]} != "
                    f"pinned {PINS[rel][:16]} (refused)")
    row, within = rel.split("/", 1)
    prov = json.load(open(os.path.join(FIXTURES, row, "PROVENANCE.json")))
    ent = [e for e in prov["files"] if e["path"] == within]
    assert ent and ent[0]["sha256"] == got, f"{rel} not pinned in PROVENANCE.json with this sha256"
    assert ent[0]["source"].startswith("identity_census/rows/"), ent[0]["source"]
    return path


def fam_from(name, spec):
    return T.Family(name, spec["loops"], spec["exts"], spec["ext_subs"],
                    spec["propagators"], "memory:" + name)


def with_legs(fam, legs):
    return T.Family(fam.name, fam.loops, legs[0], legs[1], fam.propagators, fam.source,
                    fam.physical, nu=fam.nu, leg_virt=fam.leg_virt, mass_values=fam.mass_values)


def edited(fam, old, new):
    props = [(e.replace(old, new), m) for e, m in fam.propagators]
    assert props != list(fam.propagators), f"edit {old!r} -> {new!r} changed nothing"
    return T.Family(fam.name + "_edited", fam.loops, fam.exts, fam.ext_subs, props,
                    "memory:edited", fam.physical)


def search(A, B, find_all=False, dep=True):
    det = {}
    r = T.isomorphism_to(T.momentum_vectors(A)[0], T.momentum_vectors(B)[0], len(A.loops),
                         len(A.exts), try_leg_perms=True,
                         dep_rows=T._dependent_leg_rows(B) if dep else None,
                         find_all=find_all, details=det)
    return r, det


def cyclic_ok(legp):
    n = len(legp)
    seq = list(legp)
    for r in range(n):
        rot = seq[r:] + seq[:r]
        if rot == list(range(n)) or rot == list(range(n))[::-1]:
            return True
    return False


def labeled_ok(A, B, match):
    """Test-side check that a found relabeling also carries the line masses
    (the census's (momentum, mass) multiset comparison, on the tool's primitives)."""
    vA, mA, _ = T.momentum_vectors(A)
    vB, mB, _ = T.momentum_vectors(B)
    dep = T._dependent_leg_rows(B)
    nloop, n_ext = len(A.loops), len(A.exts)
    tB = collections.Counter((T._canon_vec(v), str(sp.nsimplify(m))) for v, m in zip(vB, mB))
    tA = collections.Counter(
        (T._canon_vec(T._apply_signed_perm_to_vec(v, match["loop_perm"], match["signs"],
                                                   match["leg_perm"], nloop, n_ext, dep_rows=dep)),
         str(sp.nsimplify(m))) for v, m in zip(vA, mA))
    return tA == tB


# --------------------------------------------------------------------- pins --
def test_fixture_pins_match_provenance():
    for rel in PINS:
        fixture(rel)


def test_fixture_drift_is_refused(tmp_path):
    rel = "row20/family_of_record.yaml"
    src = fixture(rel)
    row = tmp_path / "row20"
    row.mkdir()
    (row / "family_of_record.yaml").write_bytes(open(src, "rb").read() + b"\n# tampered\n")
    prov = json.load(open(os.path.join(FIXTURES, "row20", "PROVENANCE.json")))
    (row / "PROVENANCE.json").write_text(json.dumps(prov))
    got = hashlib.sha256((row / "family_of_record.yaml").read_bytes()).hexdigest()
    assert got != PINS[rel]


# ------------------------------------------------------------ new units --
def test_dependent_leg_rows_by_value():
    four = T.Family("four", ["k1", "k2"], ["p1", "p2", "p3", "p4"], {"p4": "-p1 - p2 - p3"}, [], "m")
    five = T.Family("five", ["k1"], ["p1", "p2", "p3", "p4", "p5"], {"p5": "-p1 - p2 - p3 - p4"}, [], "m")
    comp = T.Family("comp", ["l1", "l2"], ["P12", "P56", "P34"], {"P34": "P12 - P56"}, [], "m")
    none = T.Family("none", ["k1"], ["p1", "p2"], {}, [], "m")
    assert T._dependent_leg_rows(four) == {3: [-1, -1, -1, 0]}
    assert T._dependent_leg_rows(five) == {4: [-1, -1, -1, -1, 0]}
    assert T._dependent_leg_rows(comp) == {2: [1, -1, 0]}
    assert T._dependent_leg_rows(none) == {}


def test_reimpose_conservation_by_value():
    rows = {3: [-1, -1, -1, 0]}
    one = sp.Integer(1)
    zero = sp.Integer(0)
    # k1 + p4  ->  k1 - p1 - p2 - p3
    assert T._reimpose_conservation((one, zero, zero, zero, zero, one), 2, rows) == (1, 0, -1, -1, -1, 0)
    # nothing in the dependent slot: unchanged
    v = (one, zero, one, zero, zero, zero)
    assert T._reimpose_conservation(v, 2, rows) == v
    # two dependent legs (a composite pair) both re-expressed
    rows2 = {1: [1, 0, 0], 2: [-1, 0, 0]}
    assert T._reimpose_conservation((one, zero, one, one), 1, rows2) == (1, 0, 0, 0)


def test_apply_signed_perm_reimposes_and_refuses_by_name():
    one, zero = sp.Integer(1), sp.Integer(0)
    rows = {3: [-1, -1, -1, 0]}
    # k1 + p3 under p3 -> p4 with conservation re-imposed
    v = (one, zero, zero, zero, one, zero)
    assert T._apply_signed_perm_to_vec(v, (0, 1), (1, 1), (0, 1, 3, 2), 2, 4, dep_rows=rows) == \
        (1, 0, -1, -1, -1, 0)
    # the same without dep_rows: the coefficient sits in the dependent slot (the former behaviour)
    assert T._apply_signed_perm_to_vec(v, (0, 1), (1, 1), (0, 1, 3, 2), 2, 4) == (1, 0, 0, 0, 0, 1)
    # a 2-leg vector under a 4-leg permutation that sends p1 -> p4: refused by name
    with pytest.raises(ValueError, match="external legs"):
        T._apply_signed_perm_to_vec((one, zero, one, zero), (0, 1), (1, 1), (3, 2, 1, 0), 2, 4)


def test_isomorphism_to_leg_count_guard_and_labels():
    one, zero = sp.Integer(1), sp.Integer(0)
    det = {}
    r = T.isomorphism_to([(one, one, zero)], [(one, one, zero, zero)], 1, 2,
                         try_leg_perms=True, details=det)
    assert r is None and det["leg_count_skip"] and "differ" in det["leg_count_skip"]
    det = {}
    assert T.isomorphism_to([(one, one, zero)], [(one, one, zero, zero)], 1, 2,
                            try_leg_perms=True, find_all=True, details=det) == []
    with pytest.raises(ValueError):
        T.isomorphism_to([(one, one, zero)], [(one, one, zero)], 1, 3, try_leg_perms=True)
    smi = fam_from("s", T.CATALOG_SOURCES["planar_smirnov_dbox"])
    det = {}
    T.isomorphism_to(T.momentum_vectors(smi)[0], T.momentum_vectors(smi)[0], 2, 4,
                     try_leg_perms=True, dep_rows=T._dependent_leg_rows(smi),
                     leg_labels=(["0", "0", "offshell", "offshell"], ["0", "0", "offshell", "offshell"]),
                     find_all=True, details=det)
    assert det["n_leg_perms"] == 4
    # labels (0,0,0,off) vs (0,0,off,0) admit the 6 permutations sending p4 -> p3;
    # the planar box is not its own image under p3 <-> p4 (that image is the
    # retired entry), so no match survives the restriction
    det = {}
    r = T.isomorphism_to(T.momentum_vectors(smi)[0], T.momentum_vectors(smi)[0], 2, 4,
                         try_leg_perms=True, dep_rows=T._dependent_leg_rows(smi),
                         leg_labels=(["0", "0", "0", "offshell"], ["0", "0", "offshell", "0"]),
                         details=det)
    assert det["n_leg_perms"] == 6 and r is None and det["n_relabelings_tried"] == 48
    # matching labels: identity found among the 6 admissible permutations
    det = {}
    r = T.isomorphism_to(T.momentum_vectors(smi)[0], T.momentum_vectors(smi)[0], 2, 4,
                         try_leg_perms=True, dep_rows=T._dependent_leg_rows(smi),
                         leg_labels=(["0", "0", "0", "offshell"], ["0", "0", "0", "offshell"]),
                         details=det)
    assert r == "identity" and det["n_leg_perms"] == 6 and det["leg_perm"] == (0, 1, 2, 3)
    import inspect
    assert "try_leg_signs" not in inspect.signature(T.isomorphism_to).parameters


# ---------------------------------------------------- retired catalog entry --
def test_retired_entry_found_under_p3p4_by_the_tool_search():
    ret = T.CATALOG_RETIRED["planar_smirnov_dbox_legs34_interchanged"]
    smi = fam_from("s", T.CATALOG_SOURCES["planar_smirnov_dbox"])
    rel, det = search(fam_from("r", ret), smi)
    assert rel is not None and "p3 -> p4" in rel and "p4 -> p3" in rel, rel
    assert det["leg_perm"] == (0, 1, 3, 2) and det["n_leg_perms"] == 24
    rel_old, _ = search(fam_from("r", ret), smi, dep=False)
    assert rel_old is None            # the former search: the moved leg compared unreduced


@nx_required
def test_row02_family_matches_the_planar_box_under_p3p4_with_leg_perms():
    """Row 2 (Smirnov convention t = (p1+p3)^2): with --leg-perms the family
    now matches planar_smirnov_dbox under the p3 <-> p4 leg map (the census:
    the graph is planar with its legs in the drawn order under a leg map that
    sends the family's t = (p1+p3)^2 to the paper's t = (p2+p3)^2)."""
    fam = T.load_family(fixture("row02/family_of_record.jl"))[0]
    fp = T.audit(fam, try_leg_perms=True)
    mts = [(m["family"], m["relabeling"]) for m in fp["catalog_matches"]]
    assert len(mts) == 1 and mts[0][0] == "planar_smirnov_dbox", mts
    assert "p3 -> p4" in mts[0][1] and "p4 -> p3" in mts[0][1], mts
    assert fp["canonical_hash"] == "a542a386080e" and fp["cut_signature"] == {"s": 2, "t": 4, "u": 3}
    assert not any("crossed" in f for f, _ in mts)
    fp0 = T.audit(fam, try_leg_perms=False)
    assert fp0["catalog_matches"] == []        # identity legs: no match (unchanged)


# ---------------------------------------------------------------- row 21 --
def test_row21_drawn_to_record_relabeling_by_object():
    rec = T.load_family(fixture("row21/family_of_record.json"))[0]
    drawn = with_legs(T.load_family(fixture("row21/drawn_graph.yaml"))[0], R21_LEGS)
    assert rec.exts == ["p1", "p2", "p3", "p4"] and rec.ext_subs == {"p4": "-p1-p2-p3"}
    rel, det = search(drawn, rec)
    assert rel == R21_RELABEL
    legmap = {drawn.exts[j]: rec.exts[det["leg_perm"][j]] for j in range(4)}
    assert legmap == R21_LEG_PERM
    assert det["n_relabelings_tried"] == R21_TRIED
    assert det["n_leg_perms"] == R21_LEG_PERMS
    assert cyclic_ok(det["leg_perm"])
    assert labeled_ok(drawn, rec, det["matches"][0])      # masses carried too
    rel_old, _ = search(drawn, rec, dep=False)
    assert rel_old is None            # the former search misses the p2 <-> p4 class


@nx_required
def test_row21_edited_line_controls():
    rec = T.load_family(fixture("row21/family_of_record.json"))[0]
    drawn = with_legs(T.load_family(fixture("row21/drawn_graph.yaml"))[0], R21_LEGS)
    # the plan's edit: a K4 with the legs re-seated -- every relabeling found
    # breaks the drawn cyclic leg order (census field true -> false)
    neg = edited(drawn, "k - p4 - b", "k - p3 - b")
    allm, det = search(neg, rec, find_all=True)
    assert len(allm) >= 1
    assert not any(cyclic_ok(m["leg_perm"]) for m in det["matches"]), allm
    assert T.build_graph(neg)[1] is True
    # the graph-breaking edit: no relabeling at all
    neg2 = edited(drawn, "k - p4 - b", "k - p4 + b")
    allm2, _ = search(neg2, rec, find_all=True)
    assert allm2 == []
    assert T.build_graph(neg2)[1] is False
    # positive control: p1 <-> p2 exchanged in every drawn string matches
    # under (p1 <-> p2) composed with the census map
    swapped = T.Family(drawn.name + "_p1p2", drawn.loops, drawn.exts, drawn.ext_subs,
                       [(e.replace("p1", "@").replace("p2", "p1").replace("@", "p2"), m)
                        for e, m in drawn.propagators], "memory:p1p2", drawn.physical)
    allm3, det3 = search(swapped, rec, find_all=True)
    maps = [{swapped.exts[j]: rec.exts[m["leg_perm"][j]] for j in range(4)} for m in det3["matches"]]
    assert {"p1": "p4", "p2": "p1", "p3": "p3", "p4": "p2"} in maps, maps
    assert any(m["signs"] == (-1, 1, -1) for m in det3["matches"])


# ----------------------------------------------------------------- row 5 --
def test_row05_two_relabelings_including_the_p5_moving_map():
    f5 = T.load_family(fixture("row05/family_of_record.yaml"))[0]
    d5 = T.load_family(fixture("row05/drawn_graph.yaml"))[0]
    assert f5.exts == ["p1", "p2", "p3", "p4", "p5"] and "p5" in f5.ext_subs
    allm, det = search(f5, d5, find_all=True)
    assert len(allm) == R05_N_FOUND
    assert allm[0] == R05_FIRST
    assert R05_P5_MAP in allm
    p5map = [m for m in det["matches"] if m["relabeling"] == R05_P5_MAP][0]
    assert p5map["leg_perm"] == (4, 3, 2, 1, 0) and p5map["signs"] == (-1,)
    assert cyclic_ok(p5map["leg_perm"])
    assert labeled_ok(f5, d5, p5map)
    allm_old, _ = search(f5, d5, find_all=True, dep=False)
    assert allm_old == ["identity"]            # the former search
    neg = edited(d5, "k1 + p1 + p2", "k1 + p1 + p3")
    allm_neg, _ = search(f5, neg, find_all=True)
    assert allm_neg == []


# ---------------------------------------------------- the rc-1 specimens --
@pytest.mark.parametrize("rel,n_legs,skip_expected", RC1_SPECIMENS)
def test_rc1_specimen_runs_rc0_with_leg_perms(rel, n_legs, skip_expected):
    p = fixture(rel)
    r1 = subprocess.run([sys.executable, TOOL, p, "--leg-perms", "--json"],
                        capture_output=True, text=True)
    r0 = subprocess.run([sys.executable, TOOL, p, "--json"], capture_output=True, text=True)
    assert r1.returncode == 0, r1.stderr[-1500:]
    assert r0.returncode == 0, r0.stderr[-1500:]
    assert "IndexError" not in r1.stderr
    j1, j0 = json.loads(r1.stdout)[0], json.loads(r0.stdout)[0]
    assert j1["catalog_matches"] == []
    assert len(j1["leg_set"]["legs"]) == n_legs, j1["leg_set"]
    warn = [w for w in j1["warnings"] if w.startswith("LEG-COUNT SKIP")]
    skipped = [s["family"] for s in j1["catalog_skipped_leg_count"]]
    if skip_expected:
        assert len(warn) == 1, j1["warnings"]
        assert all(e in warn[0] for e in FOUR_LEG_ENTRIES) and f"reads {n_legs} legs" in warn[0]
        assert all(e in skipped for e in FOUR_LEG_ENTRIES)
    else:
        # four legs read as data: the four-leg entries were searched, not skipped
        assert warn == [], j1["warnings"]
        assert not any(e in skipped for e in FOUR_LEG_ENTRIES), skipped
    for k in ("planar", "cut_signature", "canonical_hash", "canonical_set_size",
              "n_genuine_propagators", "n_massive_lines", "n_loops"):
        assert j1[k] == j0[k], k
    assert not [w for w in j0["warnings"] if w.startswith("LEG-COUNT SKIP")]


def test_three_leg_family_vs_four_leg_catalog_control():
    """A three-leg family (p4 never spelled) against the four-leg entries:
    rc 0, catalog_matches [], the skip named -- the pre-fix rc 1 is the flip."""
    fam = T.Family("three_leg_box", ["k1", "k2"], ["p1", "p2", "p3"], {"p3": "-p1 - p2"},
                   [("k1^2", 0), ("(k1 + p1)^2", 0), ("(k1 + p1 + p2)^2", 0), ("k2^2", 0),
                    ("(k2 + p1 + p2)^2", 0), ("(k2 - p3)^2", 0), ("(k1 - k2)^2", 0)], "memory")
    fp = T.audit(fam, try_leg_perms=True)
    assert fp["catalog_matches"] == []
    names = [s["family"] for s in fp["catalog_skipped_leg_count"]]
    assert all(e in names for e in FOUR_LEG_ENTRIES)
    assert any(w.startswith("LEG-COUNT SKIP") and "planar_smirnov_dbox (4 legs)" in w
               for w in fp["warnings"])


# ----------------------------------------------------------- rows 15 / 16 --
@pytest.mark.parametrize("rel_f,rel_d,expect", [
    ("row15/family_of_record_kira_icc_eqmass.yaml", "row15/drawn_graph.yaml", R15_RELABEL),
    ("row16/family_of_record_kira_iccg_generic.yaml", "row16/drawn_graph.yaml", R16_RELABEL),
])
@nx_required
def test_ice_cream_cone_relabeling_with_composite_legs_normalized(rel_f, rel_d, expect):
    """The drawn header's leg dictionary: p1 = P56, p2 = P34, p3 = -P12, all
    incoming, p3 = -p1 - p2  =>  P56 = p1, P12 = p1 + p2 (P34 eliminated).
    The Kira reader alone gives the two declared legs and no rule (asserted on
    the reader); the test-side normalisation to the drawn dictionary reproduces
    the census relabeling, and so does the family as load_family completes it
    (composite-leg completion: P34 added, P12 resolved outgoing, basis P56 = p1,
    P34 = p2, -P12 = p3) with no test-side rewriting."""
    Fr = T.load_kira_yaml(fixture(rel_f))[0]                   # the reader alone
    D = T.load_family(fixture(rel_d))[0]
    assert Fr.exts == ["P12", "P56"] and Fr.ext_subs == {}
    Fn = T.Family(Fr.name + "_norm", Fr.loops, ["p1", "p2", "p3"], {"p3": "-p1 - p2"},
                  [(e.replace("P56", "(p1)").replace("P12", "(p1 + p2)"), m)
                   for e, m in Fr.propagators], "memory:norm", Fr.physical)
    Dn = with_legs(D, (["p1", "p2", "p3"], {"p3": "-p1 - p2"}))
    rel, det = search(Dn, Fn)
    assert rel == expect
    allm, deta = search(Dn, Fn, find_all=True)
    assert allm == [expect] and deta["n_relabelings_tried"] == R1516_TESTED
    # the completed family (load_family) against the drawing as read: the same
    # relabeling, the same 48 relabelings tried, no rewriting in the test
    F = T.load_family(fixture(rel_f))[0]
    assert F.exts == ["P56", "P34", "P12"] and F.ext_subs == {"P12": "-P56 - P34"}, (F.exts, F.ext_subs)
    rel2, _det2 = search(D, F)
    allm2, deta2 = search(D, F, find_all=True)
    assert rel2 == expect and allm2 == [expect] and deta2["n_relabelings_tried"] == R1516_TESTED


# ---------------------------------------------------------- self-test leg --
@nx_required
def test_selftest_legperms_leg_runs_clean():
    r = T._selftest_legperms(verbose=False)
    assert r["skipped"] == 0, r["skipped_names"]
    assert r["ok"] is True, {k: v for k, v in r["results"].items() if v.get("ok") is False}
    assert r["executed"] >= 15
    assert r["expected_fail"] == 0
    assert any(fn is T._selftest_legperms for _lid, _title, fn in T.EXTRA_SELF_TEST_LEGS)


@nx_required
def test_catalog_leg_3b_is_a_plain_leg_now():
    r = T._selftest_catalog(verbose=False)
    k = "[3b] tool leg-permutation search finds p3<->p4"
    assert k in r["results"] and r["results"][k]["ok"] is True, r["results"].get(k)
    assert "expected_fail" not in r["results"][k]


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
