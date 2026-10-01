#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision.
"""
Labels battery for topology_audit.py: masses on lines, multiplicities
(propagator powers, "dots") and external-leg virtuality classes are LABELS
of the canonical propagator set and of the isomorphism search --
momentum_vectors(fam).nu, canonical_set(vecs, masses, nu) (a labeled
Counter), labeled_set / mass_multiset, isomorphism_to(..., labels=...),
labeled_isomorphism (momentum / mass / labeled levels), the fingerprint's
mass_multiset / nu_multiset / leg_labels / canonical_hash_labeled fields
(n_massive_lines derived from the mass multiset) and the catalog match's
"mass_pattern" (a mass-blind match is stated as such, never a bare identity).

Every case reproduces a hand-checked object of the identity census
(rows/<row>/iso.json, ROW.md, ROW_RECEIPT.json, IDENTITY_CENSUS.json .rows),
copied, never paraphrased; every planted negative control flips or is refused
by name.  The isomorphism class is the tool's signed-permutation one (the
leg-perms unit's conservation re-imposition included): the census identities
that need a loop-momentum basis change (its "stage B" / unimodular maps) are
NOT reached here and are marked EXPECTED-FAIL by name for the wave that adds
the affine search and the realized-graph (VF2) test -- the row-12 diagonal,
rows 20, 31, 32 and the row-14 half-turn control.  Every such mark is strict:
a pass lifts it.

  row 12  (rows/row12/iso.json "pairs" / "controls"): the three families
          sun112 / sun114 / sun123 (masses 1,1,2 / 1,1,4 / 1,2,3 from
          numeric_values) against the three panels: mass multisets equal on
          the diagonal only; every off-diagonal cell fails the labeled compare
          on the mass multiset and passes the mass-blind one under "identity"
          (the census: labelled_found false, unlabelled_found true); the
          diagonal's labeled identity is the census's stage-B map
          (loop_matrix_M [[-1,-1],[1,0]]) -> EXPECTED-FAIL
  row 14  (rows/row14/iso.json, quoted): record kite112 vs the drawing
          labelled_iso "identity"; control_kite_tm3 (equal-mass kite) fails
          the labeled test ("must fail the labelled test against the (1,1,2)
          drawing"); the census's three altered drawings by their
          (edge, mass, momentum) tables
  row 18  (rows/row18/iso.json): fam2 labeled_match true under "identity"
          and "k1 -> -k1, k2 -> -k2; p1 -> p2, p2 -> p1" (16 relabelings
          tested, leg labels on/on/off); fam5 labeled_match false,
          mass_blind_match true, mass multisets quoted
  row 4   (rows/row04/iso.json "extra"): reading (b) vs the drawing
          momentum_isomorphic true (2 relabelings, first "identity"),
          labelled_isomorphic false, the census reason sentence
  row 35  (rows/row35/iso.json): identity with n_leg_perms_allowed_by_labels
          4, n_relabelings_tried 1, nu [1]*7 both sides; control: the
          kinematics with the off-shell leg moved (p2^2 = mW2, p3^2 = 0) fails
  rows 21/28/29/33 (rows/rowNN/iso.json): the record's multiplicities carried
          (row 29 nu3 = 2, row 33 [1,1,1,2,1,2,1,1]); row 21 under the census's
          integral (ROW_RECEIPT "m18 = integrals[2]") gives the census map;
          controls: the nu-mutated JSON copies fail on the multiplicity multiset
  rows 19/31/32/34 (iso.json quoted): 19 / 34 labeled identity; 31 / 32 mass
          multisets equal, momentum identity, the labeled map is the census's
          stage-B reflection -> EXPECTED-FAIL; control: masses [1,1,9,9] fails
  row 20  (rows/row20/iso.json, quoted): leg classes (0,0,0,off) both sides,
          6 of 24 leg permutations admissible; the map is stage B -> EXPECTED-FAIL
  rows 1/7/2/5 (iso.json): (momentum, mass) pairs equal as multisets, the
          signed-relabeling notion False (as the census says), and the catalog
          match states the mass pattern ("massless catalog entry; family
          carries 6 massive lines" / "1 massive lines")
  rows 8/13 controls, row 24 cured drawn pair (leg virtualities alone),
          audit-level pattern strings, the unit values.

Runs under pytest (test_* functions) and as a battery from the tool's
--self-test leg [labels] (run_battery()).  A missing fixture is a SKIP by name;
a drifted fixture is a refusal (FixtureDrift), never a pass.

The four cases marked EXPECTED-FAIL above (the row-12 diagonal, the row-14
half-turn control, rows 31 / 32 and row 20) are plain cases since the
isomorphism unit landed: each still asserts that the signed-permutation class
finds nothing (unchanged) and then asserts the census's map from the affine
search (affine_iso_search) and the realized-graph isomorphism; the strict
xfail mechanism stays in place for any later mark (none is open).
"""
import collections
import functools
import hashlib
import json
import os
import shutil
import sys
import tempfile

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
TOOL_DIR = os.path.dirname(HERE)
sys.path.insert(0, TOOL_DIR)
import topology_audit as ta  # noqa: E402

FIX = os.path.join(HERE, "fixtures")

# sha16 of every vendored fixture this battery reads (PROVENANCE.json in the
# row directory carries the full sha256 of the copy on disk and the record's own
# sha256; both are checked on every access).
PINS = {
    "row12/family_of_record.json": "2f4f811368032eb8",
    "row12/family_of_record.yaml": "21d13748eea6a067",
    "row12/drawn_graph_112.yaml": "3afedc59b3a4ad57",
    "row12/drawn_graph_123.yaml": "e92330d766c142b0",
    "row12/drawn_graph_114.yaml": "0d902f1b2fb881b5",
    "row12/family_also/sun114_d2_tm2_g280_o11b.json": "88bfb98b87004573",
    "row12/family_also/sun123_d2_tm5_g180_o11b.json": "194ec11f95ed6d75",
    "row12/family_also/control_sun112_d2_tm5_g150_o11b.yaml": "21d13748eea6a067",
    "row12/family_also/control_sun114_d2_tm2_g280_o11b.yaml": "d480aab93ae75a57",
    "row12/iso.json": "3da3eba191d43c09",
    "row14/family_of_record.yaml": "a6d50193eabd44f3",
    "row14/drawn_graph.yaml": "e9c6eb631a18d650",
    "row14/family_also/control_kite_tm3_g150_o10.yaml": "00fa399594fc002a",
    "row14/family_also/kira14_sm2.yaml": "07a68b4b83ee5175",
    "row14/family_also/kira14_sm2.json": "fe134467a6731fe5",
    "row18/family_of_record_fam2.yaml": "c1c9c4e3ea6eb708",
    "row18/family_of_record_fam5.yaml": "fde54b3979be755b",
    "row18/drawn_graph.yaml": "d1863cd7d2580b6f",
    "row18/iso.json": "df83c519480ad21b",
    "row04/family_of_record_reading_b_4pt.yaml": "e9c7e5d7aecbaba1",
    "row04/drawn_graph.yaml": "87bcbde125b06c81",
    "row04/iso.json": "63f3370243336941",
    "row35/family_of_record.yaml": "3d3f3825b7b83326",
    "row35/kinematics_of_record.yaml": "63aee0ed75326706",
    "row35/drawn_graph.yaml": "1d180206dd0640e3",
    "row35/iso.json": "b5e635c1398b94f8",
    "row21/family_of_record.json": "209d2d1741147557",
    "row21/drawn_graph.yaml": "c307f8008add78e6",
    "row21/iso.json": "f827e23414c88a93",
    "row28/family_of_record.json": "61c4e46fb523550b",
    "row28/drawn_graph.yaml": "d69e033918b66544",
    "row28/iso.json": "1a9d634686489ae4",
    "row29/family_of_record.json": "768fd40d15f65236",
    "row29/drawn_graph.yaml": "440576817cc5bd06",
    "row29/iso.json": "4ca1cbaa9fe5f2d6",
    "row33/family_of_record.json": "460475b5991282e8",
    "row33/drawn_graph.yaml": "484ccb7f9657f612",
    "row33/iso.json": "287ba65aad4a61ac",
    "row19/family_of_record.yaml": "e50056c83e5bbb9b",
    "row19/family_kinematics.yaml": "2b383b50c72a400c",
    "row19/drawn_graph.yaml": "ab3399318550a5a6",
    "row19/drawn_kinematics.yaml": "08ccb1bef756189c",
    "row31/family_of_record.yaml": "4c748aa2bdb4c57d",
    "row31/family_kinematics.yaml": "d39dcfe38892d94a",
    "row31/drawn_graph.yaml": "4bcc152d1d63c98f",
    "row31/drawn_kinematics.yaml": "08ccb1bef756189c",
    "row32/family_of_record.yaml": "457e82a9a53b2b5f",
    "row32/family_kinematics.yaml": "d39dcfe38892d94a",
    "row32/drawn_graph.yaml": "bc0157b2f5eb579c",
    "row32/drawn_kinematics.yaml": "08ccb1bef756189c",
    "row32/iso.json": "42dd366d6d8dc21b",
    "row34/family_of_record.yaml": "14c51713c29a4f59",
    "row34/family_kinematics.yaml": "2b383b50c72a400c",
    "row34/drawn_graph.yaml": "7170eb1c9bdafff6",
    "row34/drawn_kinematics.yaml": "08ccb1bef756189c",
    "row20/family_of_record.yaml": "d4aef7e8609a1f65",
    "row20/family_kinematics.yaml": "9fb7631fb6a8455a",
    "row20/drawn_graph.yaml": "c756f4b316e8953f",
    "row20/drawn_kinematics.yaml": "b749dd5cf53fee3f",
    "row01/family_of_record.jl": "f252b561257b7839",
    "row01/drawn_graph.yaml": "160f4f486799b5a1",
    "row01/iso.json": "9f467c4aa0d65982",
    "row07/family_of_record.jl": "59dccd6ed629b5d1",
    "row07/drawn_graph.yaml": "047a94a5a2be7935",
    "row07/iso.json": "9f467c4aa0d65982",
    "row02/family_of_record.jl": "ddc75ae24954b081",
    "row02/drawn_graph.yaml": "61dd5328d6ffdcfb",
    "row02/iso.json": "5484bd046ea6e7dd",
    "row05/family_of_record.yaml": "0ec539a5fec61bdf",
    "row05/drawn_graph.yaml": "33984b814f95f701",
    "row05/iso.json": "5b0fcea4d5857d62",
    "row08/family_of_record.json": "4109f261715a4cfe",
    "row08/drawn_graph.yaml": "fa24dcb7064d9cd3",
    "row08/family_also/control_sun112_d2_tm1_g280_o11b.yaml": "a3ff093478731e70",
    "row13/family_of_record.json": "3acdcb030faec81f",
    "row13/drawn_graph.yaml": "f2e0a018ce4e788f",
    "row13/family_also/control_integralfamilies.yaml": "a6d50193eabd44f3",
    "row24/cured_adjacent_p3p4/family_of_record.json": "1864217f1e908a08",
    "row24/cured_adjacent_p3p4/drawn_graph.yaml": "deaf450f7c9f436d",
    "row24/cured_adjacent_p3p4/drawn_graph_adjacent_reading.yaml": "0ecab16ea27b1c9c",
    "row24/cured_adjacent_p3p4/iso.json": "9c5c8bd3dac8d36f",
    "row24/cured_adjacent_p3p4/iso_adjacent_reading.json": "88007a4990b3a70a",
}

# ---- census objects quoted from files that are not vendored (the iso.json
# of rows 8, 13, 14, 19, 20, 31, 34 carry a process phrase inside a record
# string); the census file is named beside each value.
# rows/row14/iso.json "pairs"[0] "labelled_iso" and "controls"
R14_LABELLED_RELABELING = "identity"
R14_MASS_MULTISET = ["0", "0", "1", "1", "2"]           # "mass_multiset_family" == "mass_multiset_drawn"
R14_CONTROLS = [
    # (control text, altered drawing (edge, mass_sq, momentum), expected_labelled, unlabelled_found)
    ("altered drawing with the sqrt2 mass moved from the rim line B-R onto the rung T-B: must "
     "fail the labelled test (the rung-vs-rim question)",
     [("L-T", "1", "k1"), ("T-R", "0", "k2"), ("L-B", "0", "-k1 + p"), ("B-R", "1", "-k2 + p"),
      ("T-B", "2", "k1 - k2")], False, True),
    ("altered drawing with the two rim masses exchanged (L-T <-> B-R): the half-turn symmetry "
     "of the kite maps it back, must pass",
     [("L-T", "2", "k1"), ("T-R", "0", "k2"), ("L-B", "0", "-k1 + p"), ("B-R", "1", "-k2 + p"),
      ("T-B", "1", "k1 - k2")], True, True),
    ("altered drawing with the massless lines adjacent: must fail",
     [("L-T", "1", "k1"), ("T-R", "2", "k2"), ("L-B", "0", "-k1 + p"), ("B-R", "0", "-k2 + p"),
      ("T-B", "1", "k1 - k2")], False, True),
]
R14_CONTROL_KITE = ("the equal-mass kite family must fail the labelled test against the (1,1,2) "
                    "drawing", False, True)     # expected_labelled, unlabelled_found
# rows/row19/iso.json, rows/row31/iso.json, rows/row32/iso.json, rows/row34/iso.json
# "drawn_vs_family" "labeled_masses_and_legs" (leg_map {p: p, pout: pB},
# respect_virtuality true, n_leg_perms_admissible 2) and the ROW.md mass multisets
R19_31_32_34 = {
    "row19": {"masses": ["1", "1", "1", "1"], "n_leg_perms_admissible": 2,
              "loop_relabeling": ["k1 -> -k1 -k2 -k3 +p", "k2 -> k3", "k3 -> k2"],
              "class": "all integer matrices with entries in {-1,0,1} and det = +-1"},
    "row31": {"masses": ["1", "1", "1", "9"], "n_leg_perms_admissible": 2,
              "loop_relabeling": ["k1 -> k3", "k2 -> k2", "k3 -> -k1 -k2 -k3 +p"],
              "class": "all integer matrices with entries in {-1,0,1} and det = +-1"},
    "row32": {"masses": ["1", "1", "1", "1", "16"], "n_leg_perms_admissible": 2,
              "loop_relabeling": ["k1 -> k1", "k2 -> k2", "k3 -> k3", "k4 -> -k1 -k2 -k3 -k4 +p"],
              "class": "signed permutations composed with the two-vertex line reflections R_i (L >= 4)"},
    "row34": {"masses": ["1", "1", "1", "1", "1"], "n_leg_perms_admissible": 2,
              "loop_relabeling": ["k1 -> k1", "k2 -> k2", "k3 -> k3", "k4 -> k4"],
              "class": "signed permutations composed with the two-vertex line reflections R_i (L >= 4)"},
}
# rows/row20/iso.json "drawn_vs_family": labeled_masses_and_legs
# n_leg_perms_admissible 6 (respect_virtuality true) / masses_only_any_leg_perm 24;
# loop_relabeling ["l1 -> k1 -k2 +p3", "l2 -> -k2 +p1 +p2 +p3"], loop_map_matrix [[1,-1],[0,-1]];
# ROW.md: "mass^2 multiset [0,0,0,1,1,1,1] ... legs p1,p2,p3 massless and p4^2 = 7/25"
R20_ADMISSIBLE, R20_ANY = 6, 24
R20_MASSES = ["0", "0", "0", "1", "1", "1", "1"]
R20_LOOP_RELABELING = ["l1 -> k1 -k2 +p3", "l2 -> -k2 +p1 +p2 +p3"]
# rows/row08/ROW.md "1 control(s) behaved as expected"; rows/row13/ROW.md "2 control(s)"
R08_DRAWN_MASSES = ["1", "1", "1"]
R08_CONTROL_MASSES = ["1", "1", "2"]
R13_DRAWN_MASSES = ["0", "0", "1", "1", "1"]
R13_CONTROL_MASSES = ["0", "0", "1", "1", "2"]
# rows/row21/ROW_RECEIPT.json why_of_record: "m18 = integrals[2] = [1,1,1,1,1,1,0x9]"
R21_CENSUS_INTEGRAL = 2
# rows/row35/drawn_graph.yaml header: "legs: p1@BL (p^2=0), p2@TL (p^2=0), p3@TR (p^2=mW2), p4@BR (p^2=mW2)"
R35_DRAWN_LEGS = {"p1": "0", "p2": "0", "p3": "mW2", "p4": "mW2"}
# rows/row18/drawn_graph.yaml header: "leg virtualities (kinematics): {p1: 0, p2: 0, p3: s}"
R18_LEGS = {"p1": "0", "p2": "0", "p3": "s"}
# rows/row24/cured_adjacent_p3p4/drawn_graph*.yaml headers: "leg virtualities: {...}",
# "external legs (all incoming): p1, p2, p3, p4, p5; p5 = -(p1+p2+p3+p4)"
R24_LEGS = (["p1", "p2", "p3", "p4", "p5"], {"p5": "-p1-p2-p3-p4"})
R24_VIRT_TILE = {"p1": "mm", "p2": "0", "p3": "mm", "p4": "0", "p5": "0"}
R24_VIRT_ADJ = {"p1": "0", "p2": "0", "p3": "mm", "p4": "mm", "p5": "0"}
# rows/row21,28,29,33/drawn_graph.yaml headers: "{'p4': '-p1-p2-p3'}"
FOUR_LEGS = (["p1", "p2", "p3", "p4"], {"p4": "-p1-p2-p3"})

CENSUS_REASON_LEGS = ("same lines and masses under some relabeling, but no such relabeling "
                      "carries the off-shell-leg labels onto each other")


class FixtureDrift(Exception):
    pass


class FixtureMissing(Exception):
    pass


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


class ExpectedFail(Exception):
    """Raised by a battery case whose object the signed-permutation class does
    not reach (named for the wave that adds the affine / realized-graph test)."""


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
    if not str(entry.get("source", "")).startswith("identity_census/rows/"):
        raise FixtureDrift(f"{rel}: PROVENANCE source is not census-relative")
    return path


def load(rel, kin=None, **kw):
    return ta.load_family(fixture(rel), kinematics=(fixture(kin) if kin else None), **kw)[0]


def census(rel):
    return json.load(open(fixture(rel)))


def with_legs(fam, legs, nu=None, leg_virt=None):
    """A copy of `fam` in the leg basis `legs` = (exts, ext_subs) (the census
    put every drawing in the record's basis before comparing), optionally with
    the drawing's edge multiplicities / leg virtualities as the census carried
    them."""
    return ta.Family(fam.name, fam.loops, legs[0], legs[1], fam.propagators, fam.source,
                     fam.physical, nu=(nu if nu is not None else fam.nu),
                     leg_virt=(leg_virt if leg_virt is not None else fam.leg_virt),
                     mass_values=fam.mass_values)


def iso(A, B, **kw):
    return ta.labeled_isomorphism(A, B, **kw)


def _tmpdir():
    return tempfile.mkdtemp(prefix="topology_audit_labels_")


def _write(path, text):
    with open(path, "w") as f:
        f.write(text)
    return path


def _found(r, level):
    v = r[level]
    return bool(v) if isinstance(v, list) else (v is not None)


# ---------------------------------------------------------------------------
#  unit values
# ---------------------------------------------------------------------------
def test_momentum_vectors_carries_nu_index_and_unpacks_as_before():
    """A JSON record's index vector rides on .nu (row 29: nu3 = 2, the doubled
    electron stub); a Kira bitmask record gets 1 per physical line, named as
    assumed; the 3-tuple unpack every caller uses is unchanged."""
    f29 = load("row29/family_of_record.json")
    mv = ta.momentum_vectors(f29)
    vecs, masses, isps = mv
    assert vecs is mv[0] and masses is mv[1] and isps == mv[2] == 6
    assert mv.nu == [1, 1, 1, 2, 1, 1, 1, 1, 1] and mv.index == list(range(9))
    assert mv.nu_source == "record"
    assert len(vecs) == 9 and len(masses) == 9
    f14 = load("row14/family_of_record.yaml")
    mv14 = ta.momentum_vectors(f14)
    assert mv14.nu == [1] * 5 and mv14.nu_source.startswith("assumed 1")
    assert isinstance(mv14, tuple) and len(mv14) == 3


def test_canonical_set_is_a_labeled_counter_by_value():
    one, zero, two = sp.Integer(1), sp.Integer(0), sp.Integer(2)
    v1, v2, v3 = (one, zero, zero), (-one, zero, one), (zero, -two, one)
    # _canon_vec unchanged: the first nonzero coefficient made positive
    assert ta._canon_vec(v2) == (1, 0, -1) and ta._canon_vec(v3) == (0, 2, -1)
    bare = ta.canonical_set([v1, v2, v3])
    assert isinstance(bare, collections.Counter)
    assert set(bare) == {(1, 0, 0), (1, 0, -1), (0, 2, -1)} and len(bare) == 3
    lab = ta.canonical_set([v1, v2, v3], ["1", "msq", "0"], [1, 2, 1])
    assert lab == collections.Counter({((1, 0, 0), "1", 1): 1, ((1, 0, -1), "msq", 2): 1,
                                       ((0, 2, -1), "0", 1): 1})
    # the same momenta with one mass moved are a different labeled set
    assert lab != ta.canonical_set([v1, v2, v3], ["msq", "1", "0"], [1, 2, 1])
    # ... and with one multiplicity moved
    assert lab != ta.canonical_set([v1, v2, v3], ["1", "msq", "0"], [2, 1, 1])
    try:
        ta.canonical_set([v1, v2], ["1"], [1, 1])
        raise AssertionError("misaligned labels accepted")
    except ValueError as exc:
        assert "aligned" in str(exc)


def test_leg_class_and_mass_label_by_value():
    assert ta._leg_class("0") == "0" and ta._leg_class("0.0") == "0"
    assert ta._leg_class("mm") == "offshell" and ta._leg_class("7/25") == "offshell"
    assert ta._leg_class("offshell") == "offshell" and ta._leg_class("mW2") == "offshell"
    assert ta._leg_class(None) is None and ta._leg_class("") is None
    f12 = load("row12/family_of_record.json")          # msq1 = 1, msq2 = 1, msq3 = 2
    assert [ta._mass_label(m, f12) for m in ta.momentum_vectors(f12)[1]] == ["1", "1", "2"]
    assert ta._mass_label(sp.Symbol("msq"), None, {"msq": "1"}) == "1"
    assert ta._mass_label("mt2") == "mt2" and ta._mass_label(0) == "0"


def test_labeled_set_and_mass_multiset_by_value():
    f14 = load("row14/family_of_record.yaml")
    L = ta.labeled_set(f14)
    assert L["mass_multiset"] == R14_MASS_MULTISET == ta.mass_multiset(f14)
    assert L["nu_multiset"] == [1] * 5 and L["n_massive_lines"] == 3 and L["n_dotted_lines"] == 0
    assert L["leg_labels"] == {"p": None, "pB": None} and L["leg_classes"] is None
    assert sum(L["lines"].values()) == 5 and all(k[2] == 1 for k in L["lines"])
    f29 = load("row29/family_of_record.json")
    L29 = ta.labeled_set(f29)
    assert L29["mass_multiset"] == ["0", "0"] + ["1"] * 7 and L29["nu_multiset"] == [1] * 8 + [2]
    assert L29["n_dotted_lines"] == 1 and L29["leg_classes"] == ["0", "0", "0", "0"]
    d28 = load("row28/drawn_graph.yaml")
    assert ta.mass_multiset(d28) == ["0", "0"] + ["msq"] * 6
    assert ta.mass_multiset(d28, {"msq": "1"}) == ["0", "0"] + ["1"] * 6


def test_fingerprint_carries_the_labels_and_n_massive_lines_is_derived():
    fp = ta.audit(load("row21/family_of_record.json"))
    assert fp["mass_multiset"] == ["0", "0", "1", "1", "1", "1"]
    assert fp["nu_multiset"] == [1, 1, 1, 1, 1, 2] and fp["n_dotted_lines"] == 1
    assert fp["nu_source"] == "record"
    assert fp["leg_labels"] == {"p1": "0", "p2": "0", "p3": "0", "p4": "0"}
    assert fp["n_massive_lines"] == 4 == sum(1 for m in fp["mass_multiset"] if m != "0")
    assert fp["canonical_hash"] == "ccb63e84071d"          # label-free, unchanged
    assert len(fp["canonical_hash_labeled"]) == 12 and fp["canonical_hash_labeled"] != fp["canonical_hash"]
    # the same momenta with the dot moved hash differently under the labeled hash only
    f2 = load("row21/family_of_record.json", integral_index=R21_CENSUS_INTEGRAL)
    fp2 = ta.audit(f2)
    assert fp2["canonical_hash"] == fp["canonical_hash"]
    assert fp2["canonical_hash_labeled"] != fp["canonical_hash_labeled"]
    assert fp2["nu_multiset"] == [1] * 6


# ---------------------------------------------------------------------------
#  row 12: the 3x3 panel matrix (item 26)
# ---------------------------------------------------------------------------
def _row12():
    fams = {"112": load("row12/family_of_record.json"),
            "114": load("row12/family_also/sun114_d2_tm2_g280_o11b.json"),
            "123": load("row12/family_also/sun123_d2_tm5_g180_o11b.json")}
    panels = {t: load(f"row12/drawn_graph_{t}.yaml") for t in ("112", "123", "114")}
    return fams, panels


def test_row12_panel_matrix_mass_multisets_are_diagonal_and_the_off_diagonal_fails_labeled():
    ref = census("row12/iso.json")
    fams, panels = _row12()
    quoted = {p["drawn_tag"]: p for p in ref["pairs"]}
    for ft, F in fams.items():
        for dt, D in panels.items():
            r = iso(F, D)
            assert r["momentum_match"] == "identity", (ft, dt, r)      # census unlabelled_iso stage A identity
            assert r["mass_multiset_A"] == quoted[ft]["mass_multiset_family"], (ft, r["mass_multiset_A"])
            assert r["mass_multiset_B"] == quoted[dt]["mass_multiset_drawn"], (dt, r["mass_multiset_B"])
            assert r["mass_multisets_equal"] is (ft == dt), (ft, dt)
            if ft != dt:
                assert r["mass_match"] is None and r["labeled_match"] is None, (ft, dt, r)
                assert "mass multiset differs" in r["reason"], r["reason"]
    # the census's two planted controls, by object
    for c in ref["controls"]:
        assert c["expected_labelled"] is False and c["unlabelled_found"] is True and c["control_ok"] is True
    r = iso(fams["112"], panels["123"])
    assert r["labeled_match"] is None and r["momentum_match"] == "identity"
    r = iso(fams["114"], panels["112"])
    assert r["labeled_match"] is None and r["momentum_match"] == "identity"


def test_row12_census_control_transcriptions_behave_as_the_records():
    """The census's control_*.yaml are the same families transcribed with the
    masses substituted; they give the same cells as the JSON records."""
    fams, panels = _row12()
    c112 = load("row12/family_also/control_sun112_d2_tm5_g150_o11b.yaml")
    c114 = load("row12/family_also/control_sun114_d2_tm2_g280_o11b.yaml")
    assert ta.mass_multiset(c112) == ta.mass_multiset(fams["112"]) == ["1", "1", "2"]
    assert ta.mass_multiset(c114) == ta.mass_multiset(fams["114"]) == ["1", "1", "4"]
    for C, own in ((c112, "112"), (c114, "114")):
        for dt, D in panels.items():
            r = iso(C, D)
            assert r["mass_multisets_equal"] is (dt == own)
            assert r["momentum_match"] == "identity"
            if dt != own:
                assert r["labeled_match"] is None


@nx_required
def test_row12_diagonal_labeled_identity_needs_the_affine_class():
    """The census found every diagonal cell at stage B (loop_matrix_M
    [[-1,-1],[1,0]], shift [[1],[0]]), a loop-momentum basis change the
    signed-permutation class does not reach: the signed search still finds no
    labeled map (unchanged), the AFFINE search (the isomorphism unit,
    affine_iso_search) finds the census's map on every diagonal cell."""
    ref = census("row12/iso.json")
    fams, panels = _row12()
    for p in ref["pairs"]:
        assert p["labelled_iso"]["found"] is True
        assert p["labelled_iso"]["loop_matrix_M"] == [[-1, -1], [1, 0]]
        assert p["labelled_iso"]["shift_matrix_S"] == [[1], [0]]
    for t in fams:
        r = iso(fams[t], panels[t])
        assert r["mass_multisets_equal"] is True and r["momentum_match"] == "identity"
        assert r["labeled_match"] is None                 # the signed class: as before
        a = ta.affine_iso_search(fams[t], panels[t])
        assert a["found"] is True, (t, a.get("reason"))
        assert a["loop_map_matrix"] == [[-1, -1], [1, 0]] and a["shift_matrix"] == [[1], [0]]
        assert a["loop_relabeling"] == ["l1 -> -k1 -k2 +p", "l2 -> k1"]
        assert ta.realized_graph_iso(fams[t], panels[t])["with_masses"] is True


# ---------------------------------------------------------------------------
#  row 14: the unequal-mass kite and the census's controls (items 26, 30)
# ---------------------------------------------------------------------------
def _kite_drawing(table):
    props = [(f"({mom})^2", m) for _e, m, mom in table]
    f = ta.Family("altered_drawing", ["k1", "k2"], ["p"], {}, props, "memory:row14-control")
    return ta.complete_legs(f)


def test_row14_record_labeled_identity_and_the_equal_mass_control_fails():
    D = load("row14/drawn_graph.yaml")
    F = load("row14/family_of_record.yaml")
    r = iso(F, D)
    assert r["labeled_match"] == R14_LABELLED_RELABELING == r["mass_match"] == r["momentum_match"]
    assert r["mass_multiset_A"] == r["mass_multiset_B"] == R14_MASS_MULTISET
    for rel in ("row14/family_also/kira14_sm2.yaml", "row14/family_also/kira14_sm2.json"):
        assert iso(load(rel), D)["labeled_match"] == "identity", rel
    text, expected_labelled, unlabelled = R14_CONTROL_KITE
    K = load("row14/family_also/control_kite_tm3_g150_o10.yaml")
    r = iso(K, D)
    assert (r["labeled_match"] is not None) is expected_labelled, (text, r)
    assert (r["momentum_match"] is not None) is unlabelled
    assert r["mass_multiset_A"] == ["0", "0", "1", "1", "1"] and "mass multiset differs" in r["reason"]


def test_row14_altered_drawings_by_the_census_tables():
    F = load("row14/family_of_record.yaml")
    for text, table, expected_labelled, unlabelled in R14_CONTROLS:
        if expected_labelled:
            continue                                   # the half-turn control: next case
        r = iso(F, _kite_drawing(table))
        assert (r["momentum_match"] is not None) is unlabelled, (text, r)
        assert r["mass_multisets_equal"] is True
        assert r["labeled_match"] is None, (text, r)


@nx_required
def test_row14_half_turn_control_needs_the_affine_class():
    """The census's 'rim masses exchanged' drawing must PASS the labelled test
    (the kite's half-turn symmetry maps it back); that symmetry is a
    loop-momentum basis change, not a signed permutation: the signed search
    finds none (unchanged), the affine search and the realized-graph
    isomorphism find it."""
    F = load("row14/family_of_record.yaml")
    text, table, expected_labelled, _u = [c for c in R14_CONTROLS if c[2]][0]
    C = _kite_drawing(table)
    r = iso(F, C)
    assert r["mass_multisets_equal"] is True and r["momentum_match"] is not None
    assert r["labeled_match"] is None                     # the signed class: as before
    a = ta.affine_iso_search(F, C)
    assert a["found"] is expected_labelled is True, (text, a.get("reason"))
    assert a["signed_permutation_match"] is None
    assert ta.realized_graph_iso(F, C)["with_masses"] is True, text


# ---------------------------------------------------------------------------
#  row 18: family 2 identity with the leg labels, family 5 mass-blind (items 2, 17)
# ---------------------------------------------------------------------------
def test_row18_fam2_labeled_identity_and_fam5_mass_blind_only():
    ref = {p["family"]: p for p in census("row18/iso.json")["pairs"] if p["drawn"] == "drawn_row18_fig_ggH"}
    D = load("row18/drawn_graph.yaml")
    F2 = load("row18/family_of_record_fam2.yaml")
    r = iso(F2, D, leg_labelsA=R18_LEGS, leg_labelsB=R18_LEGS, find_all=True)
    assert r["labeled_match"] == [m["tool_style"] for m in ref["fam2"]["labeled_matches"]]
    assert r["labeled_match"] == ["identity", "k1 -> -k1, k2 -> -k2; p1 -> p2, p2 -> p1"]
    assert r["n_relabelings_tried"]["labeled"] == ref["fam2"]["n_relabelings_tested"] == 16
    assert r["leg_classes_A"] == r["leg_classes_B"] == ["0", "0", "offshell"]
    assert r["mass_multiset_A"] == ref["fam2"]["mass_multiset_family"]
    assert r["mass_multiset_B"] == ref["fam2"]["mass_multiset_drawn"]
    F5 = load("row18/family_of_record_fam5.yaml")
    r5 = iso(F5, D, leg_labelsA=R18_LEGS, leg_labelsB=R18_LEGS)
    # iso.json: "labeled_match": false, "mass_and_label_blind_match": true,
    # "mass_labeled_leg_labels_ignored_match": false (IDENTITY_CENSUS.json .rows["18"]
    # iso_result: "labeled_match": false, "mass_blind_match": true)
    assert ref["fam5"]["labeled_match"] is False and ref["fam5"]["mass_and_label_blind_match"] is True
    assert ref["fam5"]["mass_labeled_leg_labels_ignored_match"] is False
    assert r5["labeled_match"] is None and r5["mass_match"] is None
    assert r5["momentum_match"] == ref["fam5"]["tool_isomorphism_to_mass_blind"] == "identity"
    assert r5["mass_multiset_A"] == ref["fam5"]["mass_multiset_family"] == ["0", "M2", "mt2", "mt2", "mt2", "mt2"]
    assert r5["mass_multiset_B"] == ref["fam5"]["mass_multiset_drawn"]
    assert "mass multiset differs" in r5["reason"]


# ---------------------------------------------------------------------------
#  row 4 reading (b): momentum-isomorphic, not labelled (item 23)
# ---------------------------------------------------------------------------
def test_row04_reading_b_momentum_isomorphic_true_labelled_false():
    ref = census("row04/iso.json")["extra"]["family_of_record_reading_b_4pt.yaml"]
    Fb = load("row04/family_of_record_reading_b_4pt.yaml")
    D = load("row04/drawn_graph.yaml")
    r = iso(Fb, D, find_all=True)
    assert ref["momentum_isomorphic"] is True and ref["labelled_isomorphic"] is False
    assert len(r["mass_match"]) == ref["n_momentum_relabelings_found"] == 2
    assert r["mass_match"][0] == ref["first_momentum_relabeling"] == "identity"
    assert r["labeled_match"] == [] and r["labeled_matches"] == []
    assert r["leg_classes_A"] == ["offshell", "0", "offshell", "0"]
    assert r["leg_classes_B"] == ["0", "0", "0", "0"]
    assert r["n_leg_perms_admissible"] == 0
    assert r["reason"] == ref["reason"] == CENSUS_REASON_LEGS


# ---------------------------------------------------------------------------
#  row 35: identity with the W legs off shell; the moved-leg control (item 2)
# ---------------------------------------------------------------------------
def test_row35_identity_with_leg_virtualities_and_the_moved_leg_control_fails():
    ref = census("row35/iso.json")
    F = load("row35/family_of_record.yaml", kin="row35/kinematics_of_record.yaml")
    D = load("row35/drawn_graph.yaml")
    assert F.leg_virt == {"p1": "0", "p2": "0", "p3": "mW2", "p4": "mW2"}
    r = iso(F, D, leg_labelsB=R35_DRAWN_LEGS)
    lab = ref["labeled_isomorphism_drawn_vs_record"]
    assert r["labeled_match"] == lab["relabeling"] == "identity"
    assert r["n_leg_perms_admissible"] == lab["n_leg_perms_allowed_by_labels"] == 4
    assert r["n_relabelings_tried"]["labeled"] == lab["n_relabelings_tried"] == 1
    nu_rec = [n for n in ref["record_nu_of_top_sector_integral"] if n > 0]
    assert r["nu_multiset_A"] == sorted(nu_rec) == [1] * 7
    assert r["nu_multiset_B"] == sorted(ref["drawn_edge_multiplicity"]) == [1] * 7
    assert r["mass_multiset_A"] == r["mass_multiset_B"] == ["0"] * 5 + ["1", "1"]
    # control: the off-shell leg moved (p2^2 = mW2, p3^2 = 0) in a kinematics copy
    d = _tmpdir()
    try:
        txt = open(fixture("row35/kinematics_of_record.yaml")).read()
        assert "- [[p2, p2], 0]" in txt and "- [[p3, p3], mW2]" in txt
        kin = _write(os.path.join(d, "kinematics_moved.yaml"),
                     txt.replace("- [[p2, p2], 0]", "- [[p2, p2], mW2]")
                        .replace("- [[p3, p3], mW2]", "- [[p3, p3], 0]"))
        Fm = ta.load_family(fixture("row35/family_of_record.yaml"), kinematics=kin)[0]
        assert Fm.leg_virt == {"p1": "0", "p2": "mW2", "p3": "0", "p4": "mW2"}
        rm = iso(Fm, D, leg_labelsB=R35_DRAWN_LEGS)
        assert rm["momentum_match"] == "identity" and rm["mass_match"] == "identity"
        assert rm["labeled_match"] is None and rm["reason"] == CENSUS_REASON_LEGS
        assert rm["leg_classes_A"] == ["0", "offshell", "0", "offshell"]
    finally:
        shutil.rmtree(d, ignore_errors=True)


# ---------------------------------------------------------------------------
#  rows 21 / 28 / 29 / 33: multiplicities carried (item 2)
# ---------------------------------------------------------------------------
def _drawn_in_record_basis(rel, ref):
    D = load(rel)
    return with_legs(D, FOUR_LEGS, nu=ref["drawn_edge_multiplicity"],
                     leg_virt={e: "0" for e in FOUR_LEGS[0]})


def test_row21_census_integral_gives_the_census_map_with_multiplicities_carried():
    ref = census("row21/iso.json")
    lab = ref["labeled_isomorphism_drawn_vs_record"]
    D = _drawn_in_record_basis("row21/drawn_graph.yaml", ref)
    F = load("row21/family_of_record.json", integral_index=R21_CENSUS_INTEGRAL)
    assert [n for n in ref["record_nu_of_top_sector_integral"] if n > 0] == [1] * 6
    r = iso(D, F)
    assert r["labeled_match"] == lab["relabeling"] == "k1 -> -k1, k3 -> -k3; p2 -> p4, p4 -> p2"
    assert r["n_relabelings_tried"]["labeled"] == lab["n_relabelings_tried"] == 126
    assert r["n_leg_perms_admissible"] == lab["n_leg_perms_allowed_by_labels"] == 24
    assert r["nu_multisets_equal"] and r["mass_multisets_equal"]
    # the reader's default integral (most positive indices, first among ties)
    # is integrals[0] = [1,1,2,1,1,1,...]: a dotted line the drawing has not
    Fd = load("row21/family_of_record.json")
    assert Fd.nu[:6] == [1, 1, 2, 1, 1, 1]
    rd = iso(D, Fd)
    assert rd["momentum_match"] == lab["relabeling"]
    assert rd["labeled_match"] is None and rd["nu_multisets_equal"] is False
    assert "multiplicity" in rd["reason"]


def test_rows28_29_33_identity_with_the_record_multiplicities():
    for row, nu_rec in (("row28", [1] * 8), ("row29", [1, 1, 1, 2, 1, 1, 1, 1, 1]),
                        ("row33", [1, 1, 1, 2, 1, 2, 1, 1])):
        ref = census(f"{row}/iso.json")
        assert [n for n in ref["record_nu_of_top_sector_integral"] if n > 0] == nu_rec
        assert ref["drawn_edge_multiplicity"] == nu_rec
        F = load(f"{row}/family_of_record.json")
        assert F.nu[:len(nu_rec)] == nu_rec
        D = _drawn_in_record_basis(f"{row}/drawn_graph.yaml", ref)
        r = iso(D, F, mass_mapA={"msq": "1"})
        lab = ref["labeled_isomorphism_drawn_vs_record"]
        assert r["labeled_match"] == lab["relabeling"] == "identity", (row, r)
        assert r["n_relabelings_tried"]["labeled"] == lab["n_relabelings_tried"] == 1
        assert r["n_leg_perms_admissible"] == lab["n_leg_perms_allowed_by_labels"] == 24
        assert r["nu_multiset_A"] == r["nu_multiset_B"] == sorted(nu_rec)


def _mutate_indices(rel, pos, new):
    d = _tmpdir()
    src = fixture(rel)
    js = json.load(open(src))
    assert js["integrals"][0]["indices"][pos] == 2
    js["integrals"][0]["indices"][pos] = new
    path = _write(os.path.join(d, os.path.basename(src)), json.dumps(js, indent=1))
    return d, path


def test_rows29_33_nu_mutated_copies_fail_on_the_multiplicity_multiset():
    for row, pos in (("row29", 3), ("row33", 3)):
        ref = census(f"{row}/iso.json")
        D = _drawn_in_record_basis(f"{row}/drawn_graph.yaml", ref)
        d, path = _mutate_indices(f"{row}/family_of_record.json", pos, 1)
        try:
            Fm = ta.load_family(path)[0]
            r = iso(D, Fm, mass_mapA={"msq": "1"})
            assert r["momentum_match"] == "identity"           # the momenta are untouched
            assert r["labeled_match"] is None and r["nu_multisets_equal"] is False, (row, r)
            assert "multiplicity" in r["reason"], r["reason"]
        finally:
            shutil.rmtree(d, ignore_errors=True)


# ---------------------------------------------------------------------------
#  rows 19 / 31 / 32 / 34: the bananas (item 29)
# ---------------------------------------------------------------------------
def _banana(row):
    F = load(f"{row}/family_of_record.yaml", kin=f"{row}/family_kinematics.yaml")
    D = load(f"{row}/drawn_graph.yaml", kin=f"{row}/drawn_kinematics.yaml")
    return F, D


def test_rows19_34_equal_mass_bananas_labeled_identity():
    for row in ("row19", "row34"):
        F, D = _banana(row)
        r = iso(D, F)
        assert r["mass_multiset_A"] == r["mass_multiset_B"] == R19_31_32_34[row]["masses"]
        assert r["labeled_match"] == "identity" and r["leg_labels_compared"] is True
        assert r["n_leg_perms_admissible"] == R19_31_32_34[row]["n_leg_perms_admissible"] == 2
        assert r["leg_classes_A"] == r["leg_classes_B"] == ["offshell", "offshell"]


def test_rows31_32_mass_multisets_equal_momentum_identity_and_the_control_fails():
    for row in ("row31", "row32"):
        F, D = _banana(row)
        r = iso(D, F)
        assert r["mass_multiset_A"] == r["mass_multiset_B"] == R19_31_32_34[row]["masses"]
        assert r["momentum_match"] == "identity" and r["n_leg_perms_admissible"] == 2
    # control: the row-31 family with masses [1,1,9,9] vs the drawing [1,1,1,9]
    F31, D31 = _banana("row31")
    props = list(F31.propagators)
    assert props[2] == ("k3", "1")
    props[2] = ("k3", "9")
    Fc = ta.Family("BAN_1199_control", F31.loops, F31.exts, F31.ext_subs, props, "memory:control",
                   F31.physical, leg_virt=F31.leg_virt)
    rc = iso(D31, Fc)
    assert ta.mass_multiset(Fc) == ["1", "1", "9", "9"]
    assert rc["momentum_match"] == "identity" and rc["mass_multisets_equal"] is False
    assert rc["labeled_match"] is None and "mass multiset differs" in rc["reason"]


@nx_required
def test_rows31_32_labeled_map_needs_the_affine_class():
    """The census's labeled maps carry the heavy line through a line
    reflection (row 31 'k3 -> -k1 -k2 -k3 +p', row 32 'k4 -> -k1 -k2 -k3 -k4
    +p'), not a signed permutation: the signed search finds none (unchanged),
    the affine search reports the census's map verbatim in the census's
    class, and the realized-graph isomorphism with the leg classes agrees."""
    for row in ("row31", "row32"):
        F, D = _banana(row)
        r = iso(D, F)
        assert r["labeled_match"] is None                 # the signed class: as before
        a = ta.affine_iso_search(D, F)
        assert a["found"] is True, (row, a.get("reason"))
        assert a["loop_relabeling"] == R19_31_32_34[row]["loop_relabeling"], (row, a["loop_relabeling"])
        assert a["loop_map_class"] == R19_31_32_34[row]["class"]
        assert a["leg_map"] == {"p": "p", "pout": "pB"}
        assert a["n_leg_perms_admissible"] == R19_31_32_34[row]["n_leg_perms_admissible"] == 2
        assert ta.realized_graph_iso(D, F)["with_masses_and_legs"] is True, row


# ---------------------------------------------------------------------------
#  row 20: leg classes and the admissible permutations (item 29)
# ---------------------------------------------------------------------------
def test_row20_leg_classes_admit_6_of_24_and_the_moved_leg_changes_the_classes():
    F = load("row20/family_of_record.yaml", kin="row20/family_kinematics.yaml")
    D = load("row20/drawn_graph.yaml", kin="row20/drawn_kinematics.yaml")
    r = iso(D, F)
    assert r["mass_multiset_A"] == r["mass_multiset_B"] == R20_MASSES
    assert r["leg_classes_A"] == r["leg_classes_B"] == ["0", "0", "0", "offshell"]
    assert r["n_leg_perms_admissible"] == R20_ADMISSIBLE and r["n_leg_perms_any"] == R20_ANY
    # control: the off-shell leg moved to p3 in a kinematics copy (the plan's
    # item-29 control, the kinematics unit's mutation: p3.p3 = M and the p2.p3
    # rule re-solved so that p4^2 = 0) -> virtualities {p3: M, p4: 0}
    d = _tmpdir()
    try:
        txt = open(fixture("row20/family_kinematics.yaml")).read()
        assert '- [[p3, p3], 0]' in txt and '- [[p2, p3], "(7/25-s-t)/2"]' in txt
        kin = _write(os.path.join(d, "kinematics_moved.yaml"),
                     txt.replace('- [[p3, p3], 0]', '- [[p3, p3], "M"]').replace(
                         '- [[p2, p3], "(7/25-s-t)/2"]', '- [[p2, p3], "(-M-s-t)/2"]'))
        Fm = ta.load_family(fixture("row20/family_of_record.yaml"), kinematics=kin)[0]
        assert Fm.leg_virt == {"p1": "0", "p2": "0", "p3": "M", "p4": "0"}
        rm = iso(D, Fm)
        assert rm["leg_classes_B"] == ["0", "0", "offshell", "0"] != rm["leg_classes_A"]
        assert rm["n_leg_perms_admissible"] == 6          # the classes still admit 6 maps
        assert rm["labeled_match"] is None
    finally:
        shutil.rmtree(d, ignore_errors=True)


@nx_required
def test_row20_labeled_map_needs_the_affine_class():
    """The census map is stage B (loop_map_matrix [[1,-1],[0,-1]]): the signed
    search finds none (unchanged, 48 tried over the 6 admissible leg
    permutations); the affine search reports the census's map verbatim."""
    F = load("row20/family_of_record.yaml", kin="row20/family_kinematics.yaml")
    D = load("row20/drawn_graph.yaml", kin="row20/drawn_kinematics.yaml")
    r = iso(D, F)
    assert r["labeled_match"] is None and r["n_leg_perms_admissible"] == R20_ADMISSIBLE
    a = ta.affine_iso_search(D, F)
    assert a["found"] is True and a["loop_relabeling"] == R20_LOOP_RELABELING
    assert a["loop_map_matrix"] == [[1, -1], [0, -1]]
    assert a["n_leg_perms_admissible"] == R20_ADMISSIBLE and a["leg_map"] == {p: p for p in F.exts}
    assert ta.realized_graph_iso(D, F)["with_masses_and_legs"] is True


# ---------------------------------------------------------------------------
#  rows 1 / 7 / 2 / 5: (momentum, mass) pairs and the catalog mass pattern (item 23)
# ---------------------------------------------------------------------------
def test_rows01_07_mass_pairs_carried_and_the_signed_notion_false_as_the_census_says():
    for row, masses in (("row01", ["0"] + ["msq"] * 6), ("row07", ["0"] * 6 + ["msq"])):
        ref = census(f"{row}/iso.json")["family_of_record_vs_drawn"]
        F = load(f"{row}/family_of_record.jl")
        D = load(f"{row}/drawn_graph.yaml")
        r = iso(F, D)
        assert r["mass_multiset_A"] == r["mass_multiset_B"] == masses
        assert ref["momentum_isomorphic"] is False and ref["graph"]["with_masses_and_legs"] is True
        assert r["momentum_match"] is None and r["mass_match"] is None      # routing: the affine class
        assert r["leg_classes_A"] is None and r["leg_classes_B"] == ["0", "0", "0", "0"]


def test_rows01_07_02_catalog_match_states_the_mass_pattern():
    fp1 = ta.audit(load("row01/family_of_record.jl"), try_leg_perms=True)
    m = [x for x in fp1["catalog_matches"] if x["family"] == "planar_smirnov_dbox"]
    assert len(m) == 1 and m[0]["relabeling"] == "identity" and m[0]["labels_match"] is False
    assert "massless catalog entry; family carries 6 massive lines" in m[0]["mass_pattern"]
    assert m[0]["mass_pattern"].startswith("MASS-BLIND match")
    assert any("family carries 6 massive lines" in w for w in fp1["warnings"])
    assert fp1["n_massive_lines"] == 6 and fp1["mass_multiset"] == ["0"] + ["msq"] * 6
    fp7 = ta.audit(load("row07/family_of_record.jl"), try_leg_perms=True)
    m7 = [x for x in fp7["catalog_matches"] if x["family"] == "planar_smirnov_dbox"]
    assert len(m7) == 1 and "family carries 1 massive lines" in m7[0]["mass_pattern"]
    assert fp7["n_massive_lines"] == 1
    fp2 = ta.audit(load("row02/family_of_record.jl"), try_leg_perms=True)
    m2 = [x for x in fp2["catalog_matches"] if x["family"] == "planar_smirnov_dbox"]
    assert len(m2) == 1 and fp2["n_massive_lines"] == 0
    assert m2[0]["labeled_relabeling"] is not None                 # masses and multiplicities carried
    assert "leg classes not declared by the record" in m2[0]["mass_pattern"]
    fp5 = ta.audit(load("row05/family_of_record.yaml"), try_leg_perms=True)
    assert fp5["catalog_matches"] == [] and fp5["mass_multiset"] == ["0"] * 5
    ref5 = census("row05/iso.json")["family_of_record_vs_drawn"]
    r5 = iso(load("row05/family_of_record.yaml"), load("row05/drawn_graph.yaml"))
    assert r5["labeled_match"] == ref5["labelled_relabeling"] == "identity"


def test_audit_catalog_match_labels_match_true_when_everything_is_carried():
    smi = ta.CATALOG_SOURCES["planar_smirnov_dbox"]
    fam = ta.Family("planar_dbox_control", smi["loops"], smi["exts"], smi["ext_subs"],
                    smi["propagators"], "memory:control", leg_virt=dict(smi["legs"]))
    fp = ta.audit(fam)
    m = [x for x in fp["catalog_matches"] if x["family"] == "planar_smirnov_dbox"]
    assert len(m) == 1 and m[0]["labels_match"] is True and m[0]["labeled_relabeling"] == "identity"
    assert m[0]["mass_pattern"].startswith("labels carried")
    fp8 = ta.audit(load("row08/family_of_record.json"))
    m8 = [x for x in fp8["catalog_matches"] if x["family"] == "sunrise_2l"]
    assert len(m8) == 1 and m8[0]["labels_match"] is False
    assert "massless catalog entry; family carries 3 massive lines" in m8[0]["mass_pattern"]
    fp4 = ta.audit(load("row04/family_of_record_reading_b_4pt.yaml"), try_leg_perms=True)
    assert fp4["leg_labels"] == {"p1": "offshell", "p2": "0", "p3": "offshell", "p4": "0"}


# ---------------------------------------------------------------------------
#  rows 8 / 13 controls; the row-24 cured drawn pair
# ---------------------------------------------------------------------------
def test_rows08_13_census_controls_fail_on_the_masses():
    D8 = load("row08/drawn_graph.yaml")
    r = iso(load("row08/family_of_record.json"), D8)
    assert r["labeled_match"] == "identity" and r["mass_multiset_B"] == R08_DRAWN_MASSES
    rc = iso(load("row08/family_also/control_sun112_d2_tm1_g280_o11b.yaml"), D8)
    assert rc["mass_multiset_A"] == R08_CONTROL_MASSES and rc["labeled_match"] is None
    assert rc["momentum_match"] == "identity"
    D13 = load("row13/drawn_graph.yaml")
    r13 = iso(load("row13/family_of_record.json"), D13)
    assert r13["labeled_match"] == "identity" and r13["mass_multiset_B"] == R13_DRAWN_MASSES
    rc13 = iso(load("row13/family_also/control_integralfamilies.yaml"), D13)
    assert rc13["mass_multiset_A"] == R13_CONTROL_MASSES and rc13["labeled_match"] is None
    assert rc13["momentum_match"] == "identity"


def test_row24_cured_drawn_pair_separates_on_leg_virtualities_alone():
    ref_tile = census("row24/cured_adjacent_p3p4/iso.json")
    ref_adj = census("row24/cured_adjacent_p3p4/iso_adjacent_reading.json")
    Da = with_legs(load("row24/cured_adjacent_p3p4/drawn_graph.yaml"), R24_LEGS, leg_virt=R24_VIRT_TILE)
    Db = with_legs(load("row24/cured_adjacent_p3p4/drawn_graph_adjacent_reading.yaml"), R24_LEGS,
                   leg_virt=R24_VIRT_ADJ)
    assert [e for e, _m in Da.propagators] == [e for e, _m in Db.propagators]   # the same lines
    r = iso(Da, Db)
    assert r["mass_match"] == "identity" and r["labeled_match"] is None
    assert r["leg_classes_A"] == ["offshell", "0", "offshell", "0", "0"]
    assert r["leg_classes_B"] == ["0", "0", "offshell", "offshell", "0"]
    assert r["reason"] == CENSUS_REASON_LEGS
    assert iso(Db, Db)["labeled_match"] == "identity"
    # the record family vs the two drawings: the tool's signed-permutation search
    # finds nothing (the census: "tool_isomorphism_to_signed_loop_perms_and_leg_perms_no_shifts": null);
    # the leg classes agree with the adjacent reading only
    F = load("row24/cured_adjacent_p3p4/family_of_record.json")
    assert F.leg_virt == {"p1": "0", "p2": "0", "p3": "mm", "p4": "mm", "p5": "0"}
    assert ref_tile["tool_isomorphism_to_signed_loop_perms_and_leg_perms_no_shifts"] is None
    assert ref_adj["tool_isomorphism_to_signed_loop_perms_and_leg_perms_no_shifts"] is None
    ra, rb = iso(F, Da), iso(F, Db)
    assert ra["momentum_match"] is None and rb["momentum_match"] is None
    assert ra["leg_classes_A"] == rb["leg_classes_B"] != ra["leg_classes_B"]
    assert ref_tile["labelled_graph_isomorphism"]["found"] is False
    assert ref_adj["labelled_graph_isomorphism"]["found"] is True


# ---------------------------------------------------------------------------
#  the labeled search on isomorphism_to directly (the primitive)
# ---------------------------------------------------------------------------
def test_isomorphism_to_labels_argument_by_value():
    smi = ta.CATALOG_SOURCES["planar_smirnov_dbox"]
    fam = ta.Family("s", smi["loops"], smi["exts"], smi["ext_subs"], smi["propagators"], "memory")
    vecs, _m, _i = ta.momentum_vectors(fam)
    dep = ta._dependent_leg_rows(fam)
    lab0 = [("0", 1)] * 7
    assert ta.isomorphism_to(vecs, vecs, 2, 4, dep_rows=dep, labels=(lab0, lab0)) == "identity"
    # one line massive on side A only: the mass multiset differs, nothing searched
    det = {}
    labA = [("0", 1)] * 6 + [("msq", 1)]
    assert ta.isomorphism_to(vecs, vecs, 2, 4, dep_rows=dep, labels=(labA, lab0), details=det) is None
    assert det["n_relabelings_tried"] == 0 and "mass multiset differs" in det["reason"]
    # the same multiset on a different line: searched, no relabeling carries it
    labB = [("msq", 1)] + [("0", 1)] * 6
    det = {}
    r = ta.isomorphism_to(vecs, vecs, 2, 4, try_leg_perms=True, dep_rows=dep,
                          labels=(labA, labB), details=det)
    assert r is None and det["n_relabelings_tried"] > 0
    # a dot moved: the multiplicity multiset differs by name
    det = {}
    nuA = [("0", 2)] + [("0", 1)] * 6
    assert ta.isomorphism_to(vecs, vecs, 2, 4, dep_rows=dep, labels=(nuA, lab0), details=det) is None
    assert "multiplicity" in det["reason"]
    # misaligned labels refused
    try:
        ta.isomorphism_to(vecs, vecs, 2, 4, labels=(lab0[:3], lab0))
        raise AssertionError("misaligned labels accepted")
    except ValueError:
        pass


def test_selftest_labels_leg_is_registered():
    assert any(fn is ta._selftest_labels for _lid, _title, fn in ta.EXTRA_SELF_TEST_LEGS)


# ---------------------------------------------------------------------------
BATTERY = [
    (test_momentum_vectors_carries_nu_index_and_unpacks_as_before, None),
    (test_canonical_set_is_a_labeled_counter_by_value, None),
    (test_leg_class_and_mass_label_by_value, None),
    (test_labeled_set_and_mass_multiset_by_value, None),
    (test_fingerprint_carries_the_labels_and_n_massive_lines_is_derived, None),
    (test_row12_panel_matrix_mass_multisets_are_diagonal_and_the_off_diagonal_fails_labeled, None),
    (test_row12_census_control_transcriptions_behave_as_the_records, None),
    (test_row12_diagonal_labeled_identity_needs_the_affine_class, None),
    (test_row14_record_labeled_identity_and_the_equal_mass_control_fails, None),
    (test_row14_altered_drawings_by_the_census_tables, None),
    (test_row14_half_turn_control_needs_the_affine_class, None),
    (test_row18_fam2_labeled_identity_and_fam5_mass_blind_only, None),
    (test_row04_reading_b_momentum_isomorphic_true_labelled_false, None),
    (test_row35_identity_with_leg_virtualities_and_the_moved_leg_control_fails, None),
    (test_row21_census_integral_gives_the_census_map_with_multiplicities_carried, None),
    (test_rows28_29_33_identity_with_the_record_multiplicities, None),
    (test_rows29_33_nu_mutated_copies_fail_on_the_multiplicity_multiset, None),
    (test_rows19_34_equal_mass_bananas_labeled_identity, None),
    (test_rows31_32_mass_multisets_equal_momentum_identity_and_the_control_fails, None),
    (test_rows31_32_labeled_map_needs_the_affine_class, None),
    (test_row20_leg_classes_admit_6_of_24_and_the_moved_leg_changes_the_classes, None),
    (test_row20_labeled_map_needs_the_affine_class, None),
    (test_rows01_07_mass_pairs_carried_and_the_signed_notion_false_as_the_census_says, None),
    (test_rows01_07_02_catalog_match_states_the_mass_pattern, None),
    (test_audit_catalog_match_labels_match_true_when_everything_is_carried, None),
    (test_rows08_13_census_controls_fail_on_the_masses, None),
    (test_row24_cured_drawn_pair_separates_on_leg_virtualities_alone, None),
    (test_isomorphism_to_labels_argument_by_value, None),
    (test_selftest_labels_leg_is_registered, None),
]

# pytest: the expected-fail cases are strict xfail marks (a pass lifts the mark)
try:
    import pytest as _pytest
    for _fn, _why in BATTERY:
        if _why:
            globals()[_fn.__name__] = _pytest.mark.xfail(
                strict=True, raises=ExpectedFail,
                reason=_why + " (the [isomorphism] affine / realized-graph isomorphism unit: "
                              "affine_iso_search / realized_graph_iso)")(_fn)
except ImportError:                                            # pragma: no cover
    pass


def run_battery(verbose=True):
    """Run every case; return [{"name", "status": PASS|FAIL|SKIP|EXPECTED-FAIL,
    "detail"}].  A missing fixture is SKIP (named); drift or any other
    exception is FAIL; a case declared expected-fail that raises ExpectedFail
    is EXPECTED-FAIL by name, and one that PASSES is a FAIL (the mark must be
    lifted)."""
    rows = []
    for fn, why in BATTERY:
        name = fn.__name__
        try:
            fn()
            if why:
                rows.append({"name": name, "status": "FAIL",
                             "detail": "expected to fail but PASSED: lift the mark (" + why + ")"})
            else:
                rows.append({"name": name, "status": "PASS", "detail": ""})
        except FixtureMissing as exc:
            rows.append({"name": name, "status": "SKIP", "detail": str(exc)})
        except ExpectedFail as exc:
            if why:
                rows.append({"name": name, "status": "EXPECTED-FAIL", "detail": str(exc)})
            else:
                rows.append({"name": name, "status": "FAIL", "detail": f"ExpectedFail: {exc}"})
        except Exception as exc:
            rows.append({"name": name, "status": "FAIL",
                         "detail": f"{type(exc).__name__}: {exc}"})
        if verbose:
            r = rows[-1]
            print(f"   {r['status']:13s} {name}" + (f"  ({r['detail']})" if r["detail"] else ""),
                  flush=True)
    return rows


if __name__ == "__main__":
    out = run_battery(verbose=True)
    n_fail = sum(1 for r in out if r["status"] == "FAIL")
    n_skip = sum(1 for r in out if r["status"] == "SKIP")
    n_pass = sum(1 for r in out if r["status"] == "PASS")
    n_x = sum(1 for r in out if r["status"] == "EXPECTED-FAIL")
    print(f"labels battery: {len(out)} cases, {n_pass} PASS, {n_fail} FAIL, {n_skip} SKIP, "
          f"{n_x} EXPECTED-FAIL")
    sys.exit(1 if n_fail else 0)
