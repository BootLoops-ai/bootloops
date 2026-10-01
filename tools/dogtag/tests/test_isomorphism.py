#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision.
"""
Isomorphism battery for topology_audit.py beyond signed loop permutations:
the AFFINE class k_i -> sum_m U_im k_m + sum_j c_ij p_j (unimodular_maps,
_affine_search_vectors, affine_iso_search, isomorphism_to(affine=True)) as
the relabeling DESCRIBER, the REALIZED-GRAPH isomorphism (realized_graph_iso:
VF2 on the realized labeled multigraphs, bare / with masses and
multiplicities / with the leg virtuality classes / with the leg names) as the
identity test, and the ROUTING-INDEPENDENT canonical hash of the realized
labeled graph (_canonical_form, realized_canonical_hash, the fingerprint's
canonical_hash_realized beside the routed canonical_hash).

Every case reproduces a hand-checked object of the identity census
(rows/<row>/iso.json, iso_other_variant.json, ROW.md), copied, never
paraphrased; every planted negative control flips or is refused by name.

  rows 1 / 2 / 7  (iso.json "family_of_record_vs_drawn"): graph-isomorphic
          to their drawings "bare True, with masses True, with masses and leg
          virtualities True" with the census's leg maps proven valid (and the
          tool's own map proven valid the same way); V,E 10,11 both sides;
          the REALIZED hash equal on both sides while the ROUTED hashes differ
          ("propagator sets routing-equivalent under the checker's
          signed-relabeling notion: False")
  row 20  (rows/row20/iso.json, quoted): the affine map ["l1 -> k1 -k2 +p3",
          "l2 -> -k2 +p1 +p2 +p3"], loop_map_matrix [[1, -1], [0, -1]], class
          "all integer matrices with entries in {-1,0,1} and det = +-1",
          n_leg_perms_admissible 6, leg_map identity
  rows 19 / 31 / 32 / 34 (iso.json quoted / vendored for 32): the census's
          loop_relabeling lists and leg_map {"p": "p", "pout": "pB"}, the
          L >= 4 class "signed permutations composed with the two-vertex line
          reflections R_i (L >= 4)"
  rows 9-12 (rows/row12/iso.json vendored; rows 9-11 quoted): the stage-B map
          loop_matrix_M [[-1, -1], [1, 0]], shift_matrix_S [[1], [0]]; the
          row-12 3x3 panel matrix diagonal under the affine class too
  row 14  (quoted in test_labels.py from rows/row14/iso.json): the half-turn
          control PASSES under the affine class, the other two controls fail
  row 30  (iso.json / iso_other_variant.json, both variants vendored): each
          variant graph-isomorphic to its own graphspec (leg map identity,
          mass map {"1": "m2"}), adj vs opp NOT isomorphic at every level
          (labelled found false, unlabelled found false) and under every U
  row 4   (rows/row04/iso.json): family vs drawn NOT isomorphic, "V_A": 10,
          "E_A": 12, "V_B": 12, "E_B": 14; the drawn four-point triple box
          matches planar_triple_box_3l under the affine class (the catalog
          sub-leg [8], formerly EXPECTED-FAIL)
  controls: a moved rung, a mass moved onto another line (same multiset), a
          planted shift whose target does not preserve the labeled set, a
          relabeled copy of a graph (same hash), a non-realizing family (hash
          None with the reason)

Runs under pytest (test_* functions) and as a battery from the tool's
--self-test leg [isomorphism] (run_battery()).  A missing fixture is a SKIP by name;
a drifted fixture is a refusal (FixtureDrift), never a pass.
"""
import collections
import functools
import hashlib
import itertools
import json
import os
import random
import re
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
    "row01/family_of_record.jl": "f252b561257b7839",
    "row01/drawn_graph.yaml": "160f4f486799b5a1",
    "row01/iso.json": "9f467c4aa0d65982",
    "row02/family_of_record.jl": "ddc75ae24954b081",
    "row02/drawn_graph.yaml": "61dd5328d6ffdcfb",
    "row02/iso.json": "5484bd046ea6e7dd",
    "row07/family_of_record.jl": "59dccd6ed629b5d1",
    "row07/drawn_graph.yaml": "047a94a5a2be7935",
    "row07/iso.json": "9f467c4aa0d65982",
    "row04/family_of_record.yaml": "c61b80c50011b900",
    "row04/family_of_record_reading_b_4pt.yaml": "e9c7e5d7aecbaba1",
    "row04/drawn_graph.yaml": "87bcbde125b06c81",
    "row04/iso.json": "63f3370243336941",
    "row20/family_of_record.yaml": "d4aef7e8609a1f65",
    "row20/family_kinematics.yaml": "9fb7631fb6a8455a",
    "row20/drawn_graph.yaml": "c756f4b316e8953f",
    "row20/drawn_kinematics.yaml": "b749dd5cf53fee3f",
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
    "row09/family_of_record.json": "8321e31459aa5e2b",
    "row09/drawn_graph.yaml": "2841e2ae08759216",
    "row10/family_of_record.json": "adb680a25b1132ac",
    "row10/drawn_graph.yaml": "2ad28d93667e43ec",
    "row11/family_of_record.json": "a7acd63a397d3e0f",
    "row11/drawn_graph.yaml": "3d676688bb7f44c7",
    "row12/family_of_record.json": "2f4f811368032eb8",
    "row12/family_also/sun114_d2_tm2_g280_o11b.json": "88bfb98b87004573",
    "row12/family_also/sun123_d2_tm5_g180_o11b.json": "194ec11f95ed6d75",
    "row12/drawn_graph_112.yaml": "3afedc59b3a4ad57",
    "row12/drawn_graph_114.yaml": "0d902f1b2fb881b5",
    "row12/drawn_graph_123.yaml": "e92330d766c142b0",
    "row12/iso.json": "3da3eba191d43c09",
    "row14/family_of_record.yaml": "a6d50193eabd44f3",
    "row14/drawn_graph.yaml": "e9c6eb631a18d650",
    "row14/family_also/control_kite_tm3_g150_o10.yaml": "00fa399594fc002a",
    "row30/lbl3m_adj/family_of_record.yaml": "4f1978a049287764",
    "row30/lbl3m_adj/kinematics_of_record.yaml": "9b5ce587ba9fe3ad",
    "row30/lbl3m_adj/drawn_graph.yaml": "1d1e8738e19d6de9",
    "row30/lbl3m_adj/graphspec_adj_parked_copy.json": "eb021170bf42347c",
    "row30/lbl3m_adj/iso.json": "138ad8582e71121f",
    "row30/lbl3m_adj/iso_other_variant.json": "d31e6972eb5832f0",
    "row30/lbl3m_opp/family_of_record.yaml": "9c62f7a94f1095fa",
    "row30/lbl3m_opp/kinematics_of_record.yaml": "790156fdebd95d1e",
    "row30/lbl3m_opp/drawn_graph.yaml": "19451703ecb416cf",
    "row30/lbl3m_opp/graphspec_opp_parked_copy.json": "52277a61df3959c2",
    "row30/lbl3m_opp/iso.json": "4060f1ad0870da05",
    "row30/lbl3m_opp/iso_other_variant.json": "d31e6972eb5832f0",
    "row06/retired_sevenline_vertex2L/family_of_record.json": "bdfc6ccc26612fff",
    "row21/family_of_record.json": "209d2d1741147557",
    "row21/drawn_graph.yaml": "c307f8008add78e6",
    "row21/iso.json": "f827e23414c88a93",
    "row28/family_of_record.json": "61c4e46fb523550b",
    "row28/drawn_graph.yaml": "d69e033918b66544",
    "row28/iso.json": "1a9d634686489ae4",
}

# ---- census objects quoted from files that are not vendored (their iso.json
# carry a machine path or a process phrase inside a record string; the census
# file is named beside each value; rows 19 / 31 / 34 / 20 / 14 as in
# test_labels.py).
# rows/row20/iso.json "drawn_vs_family" "labeled_masses_and_legs"
R20 = {"loop_relabeling": ["l1 -> k1 -k2 +p3", "l2 -> -k2 +p1 +p2 +p3"],
       "loop_map_matrix": [[1, -1], [0, -1]],
       "leg_map": {"p1": "p1", "p2": "p2", "p3": "p3", "p4": "p4"},
       "n_leg_perms_admissible": 6,
       "loop_map_class": "all integer matrices with entries in {-1,0,1} and det = +-1",
       "census_tried": 20088}
# rows/row19/iso.json, rows/row31/iso.json, rows/row32/iso.json (vendored),
# rows/row34/iso.json "drawn_vs_family" "labeled_masses_and_legs"
R19_31_32_34 = {
    "row19": {"loop_relabeling": ["k1 -> -k1 -k2 -k3 +p", "k2 -> k3", "k3 -> k2"],
              "loop_map_matrix": [[-1, -1, -1], [0, 0, 1], [0, 1, 0]],
              "class": "all integer matrices with entries in {-1,0,1} and det = +-1", "tried": 2804},
    "row31": {"loop_relabeling": ["k1 -> k3", "k2 -> k2", "k3 -> -k1 -k2 -k3 +p"],
              "loop_map_matrix": [[0, 0, 1], [0, 1, 0], [-1, -1, -1]],
              "class": "all integer matrices with entries in {-1,0,1} and det = +-1", "tried": 99321},
    "row32": {"loop_relabeling": ["k1 -> k1", "k2 -> k2", "k3 -> k3", "k4 -> -k1 -k2 -k3 -k4 +p"],
              "loop_map_matrix": [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [-1, -1, -1, -1]],
              "class": "signed permutations composed with the two-vertex line reflections R_i (L >= 4)",
              "tried": 366},
    "row34": {"loop_relabeling": ["k1 -> k1", "k2 -> k2", "k3 -> k3", "k4 -> k4"],
              "loop_map_matrix": [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]],
              "class": "signed permutations composed with the two-vertex line reflections R_i (L >= 4)",
              "tried": 41},
}
BANANA_LEG_MAP = {"p": "p", "pout": "pB"}
# rows/row09/iso.json, rows/row10/iso.json, rows/row11/iso.json "pairs"[0]
# "labelled_iso" (stage B); row 12's is read from the vendored iso.json
SUNRISE_STAGE_B = {"loop_matrix_M": [[-1, -1], [1, 0]], "shift_matrix_S": [[1], [0]],
                   "current_sign": 1}
# rows/row14/iso.json "controls" (as quoted in test_labels.py R14_CONTROLS)
R14_CONTROLS = [
    ("altered drawing with the sqrt2 mass moved from the rim line B-R onto the rung T-B: must "
     "fail the labelled test (the rung-vs-rim question)",
     [("L-T", "1", "k1"), ("T-R", "0", "k2"), ("L-B", "0", "-k1 + p"), ("B-R", "1", "-k2 + p"),
      ("T-B", "2", "k1 - k2")], False),
    ("altered drawing with the two rim masses exchanged (L-T <-> B-R): the half-turn symmetry "
     "of the kite maps it back, must pass",
     [("L-T", "2", "k1"), ("T-R", "0", "k2"), ("L-B", "0", "-k1 + p"), ("B-R", "1", "-k2 + p"),
      ("T-B", "1", "k1 - k2")], True),
    ("altered drawing with the massless lines adjacent: must fail",
     [("L-T", "1", "k1"), ("T-R", "2", "k2"), ("L-B", "0", "-k1 + p"), ("B-R", "0", "-k2 + p"),
      ("T-B", "1", "k1 - k2")], False),
]
FOUR_LEGS = (["p1", "p2", "p3", "p4"], {"p4": "-p1-p2-p3"})
ROW04_VE = {"V_A": 10, "E_A": 12, "V_B": 12, "E_B": 14}
UNIMODULAR_COUNTS = {1: 2, 2: 40, 3: 6960, 4: 3264}


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


def with_legs(fam, legs, leg_virt=None):
    return ta.Family(fam.name, fam.loops, legs[0], legs[1], fam.propagators, fam.source,
                     fam.physical, nu=fam.nu,
                     leg_virt=(leg_virt if leg_virt is not None else fam.leg_virt),
                     mass_values=fam.mass_values)


def relabel_legs(fam, legmap):
    """A copy of `fam` with its external legs RENAMED by legmap (a bijection
    of its leg names onto another family's): the propagator strings, the leg
    list (in the mapped order), the conservation rule and the virtualities.
    Used to PROVE a leg map: the relabeled copy is isomorphic to the other
    family with the legs matched by NAME iff the map is a valid one."""
    def ren(s):
        return re.sub(r'\b(' + "|".join(map(re.escape, legmap)) + r')\b',
                      lambda m: legmap[m.group(1)], str(s))
    props = [(ren(e), m) for e, m in fam.propagators]
    return ta.Family(fam.name + "_relabeled", fam.loops, [legmap[e] for e in fam.exts],
                     {legmap[k]: ren(v) for k, v in fam.ext_subs.items()}, props,
                     "memory:relabeled", fam.physical, nu=fam.nu,
                     leg_virt={legmap[k]: v for k, v in fam.leg_virt.items()},
                     mass_values=fam.mass_values)


def _tmpdir():
    return tempfile.mkdtemp(prefix="topology_audit_iso_")


def _write(path, text):
    with open(path, "w") as f:
        f.write(text)
    return path


def hashes(fam, **kw):
    r = ta.realized_canonical_hash(fam, **kw)
    return r["hash"], r


def family_from_lines(fam, masses=None, name="rebuilt"):
    """A copy of `fam` written from its PHYSICAL lines' momentum vectors
    (the mass in the separate field, `masses` overriding per line): the
    way to plant a mass on another line of a record whose masses sit inside
    the propagator strings (the AmflowFamily .jl form)."""
    vecs, ms, _i = ta.momentum_vectors(fam)
    syms = ta._symbols(fam.loops, fam.exts)
    names = list(fam.loops) + list(fam.exts)
    props = []
    for i, v in enumerate(vecs):
        lin = sum((sp.nsimplify(c) * syms[n] for c, n in zip(v, names)), sp.Integer(0))
        m = masses[i] if masses is not None else ms[i]
        props.append((f"({sp.sstr(lin)})^2", m))
    return ta.Family(name, fam.loops, fam.exts, fam.ext_subs, props, "memory:" + name,
                     leg_virt=fam.leg_virt, mass_values=fam.mass_values)


def valid_leg_maps(A, B, **kw):
    """Every bijection A legs -> B legs under which A, relabeled, is isomorphic
    to B with the legs matched by NAME (deterministic; VF2's own choice among
    them can vary between runs)."""
    out = []
    for perm in itertools.permutations(B.exts):
        lm = dict(zip(A.exts, perm))
        if _leg_map_is_valid(A, B, lm, **kw):
            out.append(lm)
    return out


def _leg_map_is_valid(A, B, legmap, **kw):
    """The map A leg -> B leg is valid iff A with its legs renamed by it is
    isomorphic to B with the external nodes matched by NAME."""
    r = ta.realized_graph_iso(relabel_legs(A, legmap), B, **kw)
    return r["with_masses_and_leg_names"] is True


# ---------------------------------------------------------------------------
#  unit values
# ---------------------------------------------------------------------------
def test_unimodular_maps_counts_and_class_by_value():
    for L, n in UNIMODULAR_COUNTS.items():
        mats, kind = ta.unimodular_maps(L)
        assert len(mats) == n, (L, len(mats))
        assert all(ta._int_det(M) in (1, -1) for M in mats)
        assert len({tuple(x for row in M for x in row) for M in mats}) == n     # no duplicates
        if L <= 3:
            assert kind == "all integer matrices with entries in {-1,0,1} and det = +-1"
            assert all(x in (-1, 0, 1) for M in mats for row in M for x in row)
        else:
            assert kind == ("signed permutations composed with the two-vertex line "
                            "reflections R_i (L >= 4)")
    # every signed permutation is in the class for every L (the former notion is contained)
    for L in (2, 3, 4):
        mats, _k = ta.unimodular_maps(L)
        keys = {tuple(x for row in M for x in row) for M in mats}
        for sigma in itertools.permutations(range(L)):
            for signs in itertools.product([1, -1], repeat=L):
                M = [[0] * L for _ in range(L)]
                for i in range(L):
                    M[i][sigma[i]] = signs[i]
                assert tuple(x for row in M for x in row) in keys, (L, M)
    # the census's own numbers: |U_2| = 40 matrices x 729 shifts = 29160 per leg
    # permutation bounds its row-20 count 20088 (found inside the first leg permutation)
    assert R20["census_tried"] < 40 * 3 ** 6


def test_solve_rational_and_describe_affine_by_value():
    from fractions import Fraction
    c = ta._solve_rational([[2, 0], [0, 1]], [[1, 3], [4, 5]])
    assert c == [[Fraction(1, 2), Fraction(3, 2)], [Fraction(4), Fraction(5)]]
    assert ta._solve_rational([[1, 1], [2, 2]], [[1], [2]]) is None       # singular refused
    assert ta._int_det([[1, -1], [0, -1]]) == -1 and ta._int_det([[-1, -1], [1, 0]]) == 1
    d = ta._describe_affine([[1, -1], [0, -1]], [[0, 0, 1, 0], [1, 1, 1, 0]],
                            ["l1", "l2"], ["k1", "k2"], ["p1", "p2", "p3", "p4"])
    assert d == R20["loop_relabeling"]
    d2 = ta._describe_affine([[1, 0], [0, 1]], [[2, 0], [0, -1]], ["k1", "k2"], ["k1", "k2"],
                             ["p", "pB"])
    assert d2 == ["k1 -> k1 +2*p", "k2 -> k2 -pB"]
    assert ta._describe_affine([[0]], [[0]], ["k1"], ["k1"], ["p"]) == ["k1 -> 0"]


def test_canonical_form_is_invariant_under_relabeling_and_separates_a_moved_edge():
    # the planar double box as a labeled multigraph (six vertices, four legs)
    verts = ["a1", "a2", "a3", "b1", "b2", "b3"]
    edges = [("a1", "a2", "line:m=0;nu=1"), ("a2", "a3", "line:m=0;nu=1"),
             ("b1", "b2", "line:m=0;nu=1"), ("b2", "b3", "line:m=0;nu=1"),
             ("a1", "b1", "line:m=0;nu=1"), ("a3", "b3", "line:m=0;nu=1"),
             ("a2", "b2", "line:m=msq;nu=1"),
             ("ext_p1", "a1", "leg"), ("ext_p2", "b1", "leg"), ("ext_p3", "a3", "leg"),
             ("ext_p4", "b3", "leg")]
    labels = {v: "int" for v in verts}
    labels.update({f"ext_p{j}": "ext" for j in range(1, 5)})
    enc, n_leaves, capped = ta._canonical_form(labels, edges)
    assert not capped and n_leaves >= 1
    rng = random.Random(20260904)
    for _ in range(20):
        names = list(labels)
        perm = names[:]
        rng.shuffle(perm)
        mp = dict(zip(names, perm))
        rng.shuffle(edges)
        e2 = [(mp[u], mp[v], lab) for u, v, lab in edges]
        l2 = {mp[k]: v for k, v in labels.items()}
        enc2, _n, cap2 = ta._canonical_form(l2, e2)
        assert enc2 == enc and not cap2
    # the rung moved (a2-b2 -> a1-b3): a different graph, a different form
    moved = [e for e in edges if not (set(e[:2]) == {"a2", "b2"})] + [("a1", "b3", "line:m=msq;nu=1")]
    assert ta._canonical_form(labels, moved)[0] != enc
    # the rung mass moved onto a ring line: same bare graph, different labeled form
    swapped = [(u, v, ("line:m=msq;nu=1" if {u, v} == {"a1", "a2"} else
                       "line:m=0;nu=1" if {u, v} == {"a2", "b2"} else lab))
               for u, v, lab in edges]
    assert ta._canonical_form(labels, swapped)[0] != enc
    # a parallel edge (multigraph) counts: the banana with 4 vs 5 lines differ
    ban4 = {"v0": "int", "v1": "int", "ext_p": "ext", "ext_pB": "ext"}
    e4 = [("v0", "v1", "line:m=1;nu=1")] * 4 + [("ext_p", "v0", "leg"), ("ext_pB", "v1", "leg")]
    e5 = e4[:4] + [("v0", "v1", "line:m=1;nu=1")] + e4[4:]
    assert ta._canonical_form(ban4, e4)[0] != ta._canonical_form(ban4, e5)[0]
    # the leaf cap is reported by name, never a silent hash
    enc_c, n_c, cap_c = ta._canonical_form(labels, edges, max_leaves=0)
    assert enc_c is None and cap_c is True


# ---------------------------------------------------------------------------
#  rows 1 / 2 / 7: routing-independent identity, the hashes (items 4, 22, 32)
# ---------------------------------------------------------------------------
@nx_required
def test_rows01_02_07_graph_isomorphic_with_equal_realized_and_different_routed_hashes():
    on_shell = {"p1": "0", "p2": "0", "p3": "0", "p4": "0"}
    for row in ("row01", "row02", "row07"):
        ref = census(f"{row}/iso.json")["family_of_record_vs_drawn"]
        F = load(f"{row}/family_of_record.jl")
        D = load(f"{row}/drawn_graph.yaml")
        g = ref["graph"]
        assert g["bare"] is True and g["with_masses"] is True and g["with_masses_and_legs"] is True
        assert (g["V_A"], g["E_A"], g["V_B"], g["E_B"]) == (10, 11, 10, 11)
        assert ref["momentum_isomorphic"] is False               # the signed notion, as the census says
        r = ta.realized_graph_iso(F, D)
        assert r["bare"] is True and r["with_masses"] is True, (row, r)
        assert (r["V_A"], r["E_A"], r["V_B"], r["E_B"]) == (10, 11, 10, 11)
        # the .jl record declares no leg classes: the leg-class level is None by name,
        # never inferred; with the census's leg dictionary it is True
        assert r["with_masses_and_legs"] is None and "not declared" in r["with_masses_and_legs_reason"]
        r2 = ta.realized_graph_iso(F, D, leg_labelsA=on_shell)
        assert r2["with_masses_and_legs"] is True
        # the census's leg map and the tool's own map are both VALID maps: the
        # double box has symmetries (four valid maps of the 24), VF2's choice
        # among them may differ between runs, the SET of valid maps is the object
        valid = valid_leg_maps(F, D)
        assert len(valid) == 4, (row, valid)
        assert g["with_masses_and_legs_leg_map"] in valid, (row, g)
        assert r["with_masses_leg_map"] in valid, (row, r)
        if row in ("row01", "row07"):
            assert g["with_masses_and_legs_leg_map"] == {"p3": "p1", "p4": "p2", "p1": "p3", "p2": "p4"}
        else:
            assert g["with_masses_and_legs_leg_map"] == {"p2": "p1", "p1": "p2", "p3": "p3", "p4": "p4"}
        # a wrong map is refused the same way (positive control of the prover)
        assert {"p1": "p1", "p2": "p3", "p3": "p2", "p4": "p4"} not in valid
        # the hashes: ROUTED differ, REALIZED equal
        fF, fD = ta.fingerprint(F), ta.fingerprint(D)
        assert fF["canonical_hash"] != fD["canonical_hash"], row
        assert fF["canonical_hash_realized"] == fD["canonical_hash_realized"] is not None, row
        assert fF["realized_V"] == fD["realized_V"] == 10 and fF["realized_E"] == fD["realized_E"] == 11
        assert fF["canonical_hash_realized_legs"] is None and "not declared" in fF["canonical_hash_realized_reason"]
        assert fD["canonical_hash_realized_legs"] is not None
        assert fF["realized_hash_alternatives"] == [] and fD["realized_hash_alternatives"] == []
        # the affine class also reaches the record from the drawing (a second routing)
        a = ta.affine_iso_search(F, D, leg_labelsA=on_shell)
        assert a["found"] is True and a["signed_permutation_match"] is None
        assert a["n_unimodular"] == 40 and a["loop_map_class"] == R20["loop_map_class"]
    # the massless (row 2) and outer-mass (row 1) double boxes: same bare graph, different labels
    F1, F2 = load("row01/family_of_record.jl"), load("row02/family_of_record.jl")
    r12 = ta.realized_graph_iso(F1, F2)
    assert r12["bare"] is True and r12["with_masses"] is False
    assert hashes(F1)[0] != hashes(F2)[0]


# ---------------------------------------------------------------------------
#  row 20: the affine map is the census's (item 32)
# ---------------------------------------------------------------------------
@nx_required
def test_row20_affine_map_equals_the_census_map():
    F = load("row20/family_of_record.yaml", kin="row20/family_kinematics.yaml")
    D = load("row20/drawn_graph.yaml", kin="row20/drawn_kinematics.yaml")
    a = ta.affine_iso_search(D, F)
    assert a["found"] is True
    assert a["loop_relabeling"] == R20["loop_relabeling"]
    assert a["loop_map_matrix"] == R20["loop_map_matrix"]
    assert a["leg_map"] == R20["leg_map"]
    assert a["n_leg_perms_admissible"] == R20["n_leg_perms_admissible"]
    assert a["loop_map_class"] == R20["loop_map_class"] and a["n_unimodular"] == 40
    assert a["shift_legs"] == ["p1", "p2", "p3"] and a["shift_matrix"] == [[0, 0, 1], [1, 1, 1]]
    assert a["shift_in_census_range"] is True
    assert a["signed_permutation_match"] is None           # the former notion finds nothing
    assert a["with_mass"] is True and a["respect_virtuality"] is True
    assert a["tried"] >= 1 and a["n_shift_systems_solved"] >= a["tried"]
    assert "40 unimodular matrices" in a["bound"] and "6 leg permutation" in a["bound"]
    # the map is a relabeling of the SET, verified: applying it to the drawing's
    # lines reproduces the family's labeled set (the primitive with find_all
    # lists it too)
    LA, LB = ta.labeled_set(D), ta.labeled_set(F)
    det = {}
    rel = ta.isomorphism_to(LA["vecs"], LB["vecs"], 2, 4, try_leg_perms=True,
                            dep_rows=ta._dependent_leg_rows(F),
                            leg_labels=(LA["leg_classes"], LB["leg_classes"]),
                            labels=(list(zip(LA["masses"], LA["nu"])), list(zip(LB["masses"], LB["nu"]))),
                            affine=True, names=(D.loops, F.loops, F.exts, D.exts), details=det)
    assert rel == ", ".join(R20["loop_relabeling"]) and det["relabeling_class"] == "affine unimodular loop redefinition"
    assert det["affine"]["loop_map_matrix"] == R20["loop_map_matrix"]
    # the realized-graph test agrees and the realized hashes are equal (routed differ)
    r = ta.realized_graph_iso(D, F)
    assert r["with_masses"] is True and r["with_masses_and_legs"] is True
    assert _leg_map_is_valid(D, F, R20["leg_map"]) and _leg_map_is_valid(D, F, r["with_masses_and_legs_leg_map"])
    fF, fD = ta.fingerprint(F), ta.fingerprint(D)
    assert fF["canonical_hash"] == "3c695d99a1e1" != fD["canonical_hash"]
    assert fF["canonical_hash_realized"] == fD["canonical_hash_realized"] is not None
    assert fF["canonical_hash_realized_legs"] == fD["canonical_hash_realized_legs"] is not None
    # the moved off-shell leg (the labels unit's control: the family's kinematics
    # with p3 off shell instead of p4).  A RESULT of the full class: the zgamma
    # graph has a half-turn automorphism (TL,TM,TR) <-> (BL,BM,BR) that carries
    # p3 <-> p4 and p1 <-> p2, so the leg moved from p4 to p3 is the SAME
    # labeled graph under the leg map {p3 -> p4, p4 -> p3} -- the signed class
    # (the labels unit's assertion, unchanged) finds nothing, the affine class
    # and VF2 with the leg classes find the half-turn; the leg moved to p1 or
    # p2 (a vertex with massless lines only) IS refused at the class level
    # while the mass-only level still passes.
    d = _tmpdir()
    try:
        txt = open(fixture("row20/family_kinematics.yaml")).read()
        for leg, expect in (("p3", True), ("p1", False), ("p2", False)):
            kin = _write(os.path.join(d, f"kinematics_moved_{leg}.yaml"),
                         txt.replace(f'- [[{leg}, {leg}], 0]', f'- [[{leg}, {leg}], "M"]').replace(
                             '- [[p2, p3], "(7/25-s-t)/2"]', '- [[p2, p3], "(-M-s-t)/2"]'))
            Fm = ta.load_family(fixture("row20/family_of_record.yaml"), kinematics=kin)[0]
            assert Fm.leg_virt[leg] == "M" and Fm.leg_virt["p4"] == "0"
            am = ta.affine_iso_search(D, Fm)
            assert am["found"] is expect and am["n_leg_perms_admissible"] == 6, (leg, am.get("reason"))
            rm = ta.realized_graph_iso(D, Fm)
            assert rm["with_masses"] is True and rm["with_masses_and_legs"] is expect, leg
            if expect:
                assert am["leg_map"] == {"p1": "p1", "p2": "p2", "p3": "p4", "p4": "p3"}
                assert _leg_map_is_valid(D, Fm, {"p1": "p2", "p2": "p1", "p3": "p4", "p4": "p3"})
                assert ta.labeled_isomorphism(D, Fm)["labeled_match"] is None    # the signed class
            am2 = ta.affine_iso_search(D, Fm, respect_leg_classes=False)
            assert am2["found"] is True and am2["n_leg_perms_admissible"] == 24
    finally:
        shutil.rmtree(d, ignore_errors=True)


# ---------------------------------------------------------------------------
#  rows 19 / 31 / 32 / 34: the bananas' maps (item 32; rows 31 / 32 formerly EXPECTED-FAIL)
# ---------------------------------------------------------------------------
def _banana(row):
    F = load(f"{row}/family_of_record.yaml", kin=f"{row}/family_kinematics.yaml")
    D = load(f"{row}/drawn_graph.yaml", kin=f"{row}/drawn_kinematics.yaml")
    return F, D


@nx_required
def test_rows19_31_32_34_affine_maps_equal_the_census_maps():
    for row, ref in R19_31_32_34.items():
        F, D = _banana(row)
        a = ta.affine_iso_search(D, F)
        assert a["found"] is True, (row, a.get("reason"))
        assert a["loop_relabeling"] == ref["loop_relabeling"], (row, a["loop_relabeling"])
        assert a["loop_map_matrix"] == ref["loop_map_matrix"]
        assert a["leg_map"] == BANANA_LEG_MAP and a["n_leg_perms_admissible"] == 2
        assert a["loop_map_class"] == ref["class"]
        assert a["n_unimodular"] == UNIMODULAR_COUNTS[len(F.loops)]
        r = ta.realized_graph_iso(D, F)
        assert r["with_masses"] is True and r["with_masses_and_legs"] is True
        assert _leg_map_is_valid(D, F, BANANA_LEG_MAP)          # the census's map, proven
        assert _leg_map_is_valid(D, F, r["with_masses_and_legs_leg_map"])
        assert hashes(D)[0] == hashes(F)[0] is not None
    ref32 = census("row32/iso.json")["drawn_vs_family"]["labeled_masses_and_legs"]
    assert ref32["loop_relabeling"] == R19_31_32_34["row32"]["loop_relabeling"]
    assert ref32["loop_map_matrix"] == R19_31_32_34["row32"]["loop_map_matrix"]
    assert ref32["leg_map"] == BANANA_LEG_MAP and ref32["loop_map_class"] == R19_31_32_34["row32"]["class"]
    # the row-31 mass-permuted control (masses [1,1,9,9] vs [1,1,1,9]): refused
    # at the label multiset, nothing searched; the hashes differ
    F31, D31 = _banana("row31")
    props = list(F31.propagators)
    assert props[2] == ("k3", "1")
    props[2] = ("k3", "9")
    Fc = ta.Family("BAN_1199_control", F31.loops, F31.exts, F31.ext_subs, props, "memory:control",
                   F31.physical, leg_virt=F31.leg_virt)
    ac = ta.affine_iso_search(D31, Fc)
    assert ac["found"] is False and "label multisets differ" in ac["reason"] and ac["tried"] == 0
    assert ta.realized_graph_iso(D31, Fc)["with_masses"] is False
    assert hashes(D31)[0] != hashes(Fc)[0]


# ---------------------------------------------------------------------------
#  rows 9-12: the sunrise stage-B map; the row-12 diagonal (items 26, 32)
# ---------------------------------------------------------------------------
@nx_required
def test_rows09_12_sunrise_maps_equal_the_census_stage_b():
    ref12 = census("row12/iso.json")
    for p in ref12["pairs"]:
        assert p["labelled_iso"]["found"] is True
        assert p["labelled_iso"]["loop_matrix_M"] == SUNRISE_STAGE_B["loop_matrix_M"]
        assert p["labelled_iso"]["shift_matrix_S"] == SUNRISE_STAGE_B["shift_matrix_S"]
        assert p["labelled_iso"]["current_sign"] == SUNRISE_STAGE_B["current_sign"]
    pairs = [("row09", "row09/family_of_record.json", "row09/drawn_graph.yaml"),
             ("row10", "row10/family_of_record.json", "row10/drawn_graph.yaml"),
             ("row11", "row11/family_of_record.json", "row11/drawn_graph.yaml"),
             ("row12", "row12/family_of_record.json", "row12/drawn_graph_112.yaml"),
             ("row12/114", "row12/family_also/sun114_d2_tm2_g280_o11b.json", "row12/drawn_graph_114.yaml"),
             ("row12/123", "row12/family_also/sun123_d2_tm5_g180_o11b.json", "row12/drawn_graph_123.yaml")]
    for tag, frel, drel in pairs:
        F, D = load(frel), load(drel)
        a = ta.affine_iso_search(F, D)
        assert a["found"] is True, (tag, a.get("reason"))
        assert a["loop_map_matrix"] == SUNRISE_STAGE_B["loop_matrix_M"], (tag, a["loop_map_matrix"])
        assert a["shift_matrix"] == SUNRISE_STAGE_B["shift_matrix_S"] and a["shift_legs"] == ["p"], tag
        assert a["loop_relabeling"] == ["l1 -> -k1 -k2 +p", "l2 -> k1"]
        assert a["leg_map"] == {"p": "p", "pB": "pB"}
        assert a["signed_permutation_match"] is None               # stage A finds nothing labeled
        r = ta.realized_graph_iso(F, D)
        assert r["with_masses"] is True
        assert hashes(F)[0] == hashes(D)[0] is not None, tag


@nx_required
def test_row12_panel_matrix_is_diagonal_under_the_affine_class_too():
    fams = {"112": load("row12/family_of_record.json"),
            "114": load("row12/family_also/sun114_d2_tm2_g280_o11b.json"),
            "123": load("row12/family_also/sun123_d2_tm5_g180_o11b.json")}
    panels = {t: load(f"row12/drawn_graph_{t}.yaml") for t in ("112", "114", "123")}
    for ft, F in fams.items():
        for dt, D in panels.items():
            a = ta.affine_iso_search(F, D)
            r = ta.realized_graph_iso(F, D)
            assert a["found"] is (ft == dt), (ft, dt, a.get("reason"))
            assert r["with_masses"] is (ft == dt) and r["bare"] is True, (ft, dt)
            assert (hashes(F)[0] == hashes(D)[0]) is (ft == dt)
            if ft != dt:
                assert "label multisets differ" in a["reason"]
                # mass-blind the panels are one graph (the census's unlabelled identity)
                assert ta.affine_iso_search(F, D, with_mass=False)["found"] is True


# ---------------------------------------------------------------------------
#  row 14: the half-turn control passes, the other two fail (item 26)
# ---------------------------------------------------------------------------
def _kite_drawing(table):
    props = [(f"({mom})^2", m) for _e, m, mom in table]
    f = ta.Family("altered_drawing", ["k1", "k2"], ["p"], {}, props, "memory:row14-control")
    return ta.complete_legs(f)


@nx_required
def test_row14_half_turn_control_passes_and_the_others_fail_under_the_affine_class():
    F = load("row14/family_of_record.yaml")
    D = load("row14/drawn_graph.yaml")
    a0 = ta.affine_iso_search(F, D)
    assert a0["found"] is True and a0["signed_permutation_match"] == "identity"
    for text, table, expected in R14_CONTROLS:
        C = _kite_drawing(table)
        a = ta.affine_iso_search(F, C)
        r = ta.realized_graph_iso(F, C)
        assert a["found"] is expected, (text, a.get("reason"))
        assert r["with_masses"] is expected and r["bare"] is True, (text, r)
        assert (hashes(F)[0] == hashes(C)[0]) is expected, text
        if not expected:
            # the search exhausts and says so (the mass multisets agree; no U carries the labels)
            assert a["mass_multiset_A"] == a["mass_multiset_B"]
            assert "no map in the class" in a["reason"] and "%d matrices survive" % a["n_U_surviving_loop_block_filter"] in a["reason"]
    K = load("row14/family_also/control_kite_tm3_g150_o10.yaml")
    ak = ta.affine_iso_search(K, D)
    assert ak["found"] is False and "label multisets differ" in ak["reason"]
    assert ta.realized_graph_iso(K, D)["with_masses"] is False


# ---------------------------------------------------------------------------
#  row 30: adj vs opp NOT isomorphic; each variant is its own graphspec (item 4)
# ---------------------------------------------------------------------------
def _row30():
    FA = load("row30/lbl3m_adj/family_of_record.yaml", kin="row30/lbl3m_adj/kinematics_of_record.yaml")
    FO = load("row30/lbl3m_opp/family_of_record.yaml", kin="row30/lbl3m_opp/kinematics_of_record.yaml")
    GA = load("row30/lbl3m_adj/graphspec_adj_parked_copy.json")
    GO = load("row30/lbl3m_opp/graphspec_opp_parked_copy.json")
    DA = with_legs(load("row30/lbl3m_adj/drawn_graph.yaml"), FOUR_LEGS,
                   leg_virt={"p1": "0", "p2": "0", "p3": "0", "p4": "0"})
    DO = with_legs(load("row30/lbl3m_opp/drawn_graph.yaml"), FOUR_LEGS,
                   leg_virt={"p1": "0", "p2": "0", "p3": "0", "p4": "0"})
    return FA, FO, GA, GO, DA, DO


@nx_required
def test_row30_adj_vs_opp_not_isomorphic_and_each_variant_is_its_own_graphspec():
    FA, FO, GA, GO, DA, DO = _row30()
    m2 = {"m2": "1"}
    for var, F, G, D in (("adj", FA, GA, DA), ("opp", FO, GO, DO)):
        ref = census(f"row30/lbl3m_{var}/iso.json")
        assert ref["labelled_graph_isomorphism"]["found"] is True
        assert ref["labelled_graph_isomorphism"]["chosen"]["leg_map"] == {"p1": "p1", "p2": "p2", "p3": "p3", "p4": "p4"}
        assert ref["labelled_graph_isomorphism"]["chosen"]["mass_map"] == {"1": "m2", "0": "0"}
        assert ref["tool_isomorphism_to_signed_loop_perms_and_leg_perms_no_shifts"] is None
        for other in (G, D):
            r = ta.realized_graph_iso(F, other, mass_mapB=m2)
            assert r["with_masses_and_legs"] is True and r["with_masses_and_leg_names"] is True, (var, r)
            assert (r["V_A"], r["E_A"]) == (r["V_B"], r["E_B"]) == (12, 14)
            a = ta.affine_iso_search(F, other, mass_mapB=m2)
            assert a["found"] is True and a["leg_map"] == {"p1": "p1", "p2": "p2", "p3": "p3", "p4": "p4"}, var
            # (the census's "tool_isomorphism_to ... no_shifts": null was the former
            # tool on the three-leg inferred set; with the kinematics' four legs the
            # graphspec routing is the record's own up to a loop renaming, a signed
            # permutation finds it first; the census's spanning-tree routing of the
            # drawn yaml needs the affine class)
            assert (a["signed_permutation_match"] is not None) is (other is G), var
            assert a["n_unimodular"] == 6960
            assert hashes(F)[0] == hashes(other, mass_map=m2)[0] is not None, var
    # the cross control: the census's iso_other_variant.json
    for var in ("adj", "opp"):
        ref = census(f"row30/lbl3m_{var}/iso_other_variant.json")
        assert ref["labelled_graph_isomorphism"]["found"] is False
        assert ref["unlabelled_graph_isomorphism"]["found"] is False
    for A, B, kw, tag in ((FA, FO, {}, "families"), (GA, GO, {}, "graphspecs"),
                          (FA, GO, {"mass_mapB": m2}, "adj family vs opp graphspec")):
        r = ta.realized_graph_iso(A, B, **kw)
        assert r["bare"] is False and r["with_masses"] is False, (tag, r)
        assert r["with_masses_and_legs"] is False and r["with_masses_and_leg_names"] is False
        a = ta.affine_iso_search(A, B, **kw)
        assert a["found"] is False and a["tried"] == 0, (tag, a)
        assert "no map in the class" in a["reason"]
        assert hashes(A)[0] != hashes(B, mass_map=kw.get("mass_mapB"))[0]
    # the chord hops on the massive ring, as the census read them: [2, 3] vs [2, 4]
    a_adj = ta.affine_iso_search(GA, GO)
    assert a_adj["found"] is False and a_adj["n_unimodular"] == 6960


# ---------------------------------------------------------------------------
#  row 4: NOT isomorphic; the drawn triple box matches the catalog (items 20, 22)
# ---------------------------------------------------------------------------
@nx_required
def test_row04_family_vs_drawn_not_isomorphic_with_the_census_counts():
    ref = census("row04/iso.json")
    g = ref["family_of_record_vs_drawn"]["graph"]
    assert g["bare"] is False and (g["V_A"], g["E_A"], g["V_B"], g["E_B"]) == (10, 12, 12, 14)
    F = load("row04/family_of_record.yaml")
    D = load("row04/drawn_graph.yaml")
    r = ta.realized_graph_iso(F, D)
    assert r["bare"] is False and r["with_masses"] is False and r["with_masses_and_legs"] is False
    assert {k: r[k] for k in ROW04_VE} == ROW04_VE
    assert (r["V_internal_A"], r["E_lines_A"], r["V_internal_B"], r["E_lines_B"]) == (7, 9, 8, 10)
    a = ta.affine_iso_search(F, D)
    assert a["found"] is False and "counts differ" in a["reason"]
    assert hashes(F)[0] != hashes(D)[0]
    # reading (b): the same ten lines and masses (bare / with masses True), the
    # leg classes differ (False), as the census's "extra" says
    ex = ref["extra"]["family_of_record_reading_b_4pt.yaml"]["graph"]
    assert ex["bare"] is True and ex["with_masses"] is True and ex["with_masses_and_legs"] is False
    Fb = load("row04/family_of_record_reading_b_4pt.yaml")
    rb = ta.realized_graph_iso(Fb, D)
    assert rb["bare"] is True and rb["with_masses"] is True and rb["with_masses_and_legs"] is False
    assert rb["with_masses_leg_map"] is not None
    assert hashes(Fb)[0] == hashes(D)[0]                        # the bare-legs hash agrees
    assert hashes(Fb, with_leg_classes=True)[0] != hashes(D, with_leg_classes=True)[0]


@nx_required
def test_row04_drawn_matches_planar_triple_box_in_the_audit():
    D = load("row04/drawn_graph.yaml")
    fp = ta.audit(D, try_leg_perms=True)
    m = [x for x in fp["catalog_matches"] if x["family"] == "planar_triple_box_3l"]
    assert len(m) == 1, fp["catalog_matches"]
    m = m[0]
    assert m["relabeling_class"] == "affine unimodular loop redefinition"
    assert m["relabeling"] is not None and m["labels_match"] is True
    assert m["affine"]["loop_map_class"] == R20["loop_map_class"]
    assert m["affine"]["leg_map"] == {"p1": "p1", "p2": "p2", "p3": "p3", "p4": "p4"}
    rg = m["realized_graph_iso"]
    assert rg["bare"] is True and rg["with_masses"] is True and rg["with_masses_and_legs"] is True
    assert (rg["V_A"], rg["E_A"], rg["V_B"], rg["E_B"]) == (12, 14, 12, 14)
    assert any("affine loop redefinition" in w and "planar_triple_box_3l" in w for w in fp["warnings"])
    assert all(x["family"] != "ud_ladder_3pt_offshell_3l" for x in fp["catalog_matches"])
    pl = ta.planarity(D)
    assert (pl["n_vertices"], pl["n_edges"]) == (12, 14)
    # the map is the one the search reports: applying it to the drawing's lines
    # reproduces the catalog entry's set (find_all lists it among the maps)
    tb = ta.CATALOG_SOURCES["planar_triple_box_3l"]
    TB = ta.Family("tb", tb["loops"], tb["exts"], tb["ext_subs"], tb["propagators"], "memory:catalog",
                   leg_virt=dict(tb["legs"]))
    a = ta.affine_iso_search(D, TB, find_all=True)
    assert a["found"] is True and m["relabeling"] in [x["relabeling"] for x in a["matches"]]
    # the family of record (three legs, nine lines) matches the ud entry under the
    # identity and nothing under the affine class beyond it
    F = load("row04/family_of_record.yaml")
    fpF = ta.audit(F, try_leg_perms=True)
    assert [(x["family"], x["relabeling"], x["relabeling_class"]) for x in fpF["catalog_matches"]] == \
        [("ud_ladder_3pt_offshell_3l", "identity", "signed loop permutation")]


# ---------------------------------------------------------------------------
#  planted controls: the search exhausts / refuses by name
# ---------------------------------------------------------------------------
DBOX_SPEC = {"name": "dbox_edge_list",
             "edges": [["a1", "a2", 0], ["a2", "a3", 0], ["b1", "b2", 0], ["b2", "b3", 0],
                       ["a1", "b1", 0], ["a3", "b3", 0], ["a2", "b2", "msq"]],
             "legs": [["p2", "a1", "0"], ["p1", "b1", "0"], ["p3", "a3", "0"], ["p4", "b3", "0"]],
             "loop_momenta": ["k1", "k2"]}


def _fam_from_spec(d, spec, name):
    p = os.path.join(d, name + ".json")
    with open(p, "w") as f:
        json.dump(spec, f)
    return ta.load_family(p)[0]


@nx_required
def test_moved_rung_is_not_isomorphic_under_every_U_and_the_search_says_so():
    """The row-7 drawing (the C3 double box, rung massive) as an edge list,
    routed by the tool; (i) the rung and a rail re-attached crosswise (a2-b2,
    a3-b3 -> a2-b3, a3-b2): a different graph with the same line count, masses
    and leg count -- not isomorphic under every U (the search exhausts and
    says so), VF2 False, a different hash; (ii) the plan's literal move a2-b2
    -> a1-b3 leaves a2 and b2 two-valent: the realizer refuses it by name
    (NON-GRAPH, repeated line) and the set search finds nothing."""
    D7 = load("row07/drawn_graph.yaml")
    d = _tmpdir()
    try:
        D0 = _fam_from_spec(d, DBOX_SPEC, "orig")
        assert ta.realized_graph_iso(D0, D7)["with_masses"] is True         # the same drawing
        assert ta.affine_iso_search(D0, D7)["found"] is True
        crossed = dict(DBOX_SPEC, edges=[["a1", "a2", 0], ["a2", "a3", 0], ["b1", "b2", 0], ["b2", "b3", 0],
                                         ["a1", "b1", 0], ["a3", "b2", 0], ["a2", "b3", "msq"]])
        Dc = _fam_from_spec(d, crossed, "crossed")
        a = ta.affine_iso_search(D7, Dc)
        assert a["found"] is False and a["mass_multiset_A"] == a["mass_multiset_B"]
        assert "no map in the class" in a["reason"] and a["n_unimodular"] == 40
        assert a["n_leg_perms_admissible"] == 24
        assert ("%d matrices survive the loop-block filter" % a["n_U_surviving_loop_block_filter"]) in a["reason"]
        r = ta.realized_graph_iso(D7, Dc)
        assert r["bare"] is False and r["with_masses"] is False
        assert hashes(D7)[0] != hashes(Dc)[0]
        F = load("row07/family_of_record.jl")
        assert ta.affine_iso_search(F, Dc)["found"] is False
        assert ta.realized_graph_iso(F, Dc)["bare"] is False
        literal = dict(DBOX_SPEC, edges=[["a1", "a2", 0], ["a2", "a3", 0], ["b1", "b2", 0], ["b2", "b3", 0],
                                         ["a1", "b1", 0], ["a3", "b3", 0], ["a1", "b3", "msq"]])
        Dl = _fam_from_spec(d, literal, "literal")
        G, ok = ta.build_graph(Dl)
        assert ok is False and G.graph["realization"]["verdict"] == "NON-GRAPH"
        assert "repeated line momentum" in G.graph["realization"]["cause"]
        rl = ta.realized_graph_iso(D7, Dl)
        assert rl["bare"] is None and "vertex realization failed for B" in rl["error"]
        assert ta.affine_iso_search(D7, Dl)["found"] is False
        assert hashes(Dl)[0] is None and "realization failed" in hashes(Dl)[1]["reason"]
    finally:
        shutil.rmtree(d, ignore_errors=True)


@nx_required
def test_a_different_mass_on_one_line_is_refused_under_every_U():
    F = load("row07/family_of_record.jl")
    vecs, ms, _i = ta.momentum_vectors(F)
    labels = ta.labeled_set(F)["masses"]
    assert labels.count("msq") == 1 and len(labels) == 7
    i0 = labels.index("msq")
    j0 = (i0 + 1) % 7
    F0 = family_from_lines(F, name="row07_rebuilt")
    assert ta.affine_iso_search(F, F0)["signed_permutation_match"] == "identity"
    # the rung's mass moved onto a ring line: the same mass multiset, no U carries it
    moved = list(ms)
    moved[i0], moved[j0] = ms[j0], ms[i0]
    Fm = family_from_lines(F, masses=moved, name="mass_moved_control")
    assert ta.mass_multiset(Fm) == ta.mass_multiset(F)
    a = ta.affine_iso_search(F, Fm)
    assert a["found"] is False and "no map in the class" in a["reason"]
    r = ta.realized_graph_iso(F, Fm)
    assert r["with_masses"] is False and r["bare"] is True
    assert hashes(F)[0] != hashes(Fm)[0]
    # a different mass VALUE on one line: the multiset differs, refused by name, nothing searched
    value = list(ms)
    value[i0] = "M"
    Fv = family_from_lines(F, masses=value, name="mass_value_control")
    av = ta.affine_iso_search(F, Fv)
    assert av["found"] is False and "label multisets differ" in av["reason"] and av["tried"] == 0
    assert ta.realized_graph_iso(F, Fv)["with_masses"] is False


@nx_required
def test_planted_shift_is_found_and_a_tampered_target_is_refused():
    """A copy of the row-20 family with k1 -> k1 + p1 in every line is a second
    routing of the same graph: the search finds the shift (solved, not
    enumerated).  The same copy with one line's momentum tampered no longer
    preserves the labeled set: the shift systems are solved, every candidate
    is verified against the full set and refused."""
    F = load("row20/family_of_record.yaml", kin="row20/family_kinematics.yaml")
    def shifted(expr):
        return re.sub(r'\bk1\b', "(k1+p1)", expr)
    props = [(shifted(e), m) for e, m in F.propagators]
    Fs = ta.Family("shifted", F.loops, F.exts, F.ext_subs, props, "memory:control", F.physical,
                   leg_virt=F.leg_virt)
    a = ta.affine_iso_search(F, Fs)
    assert a["found"] is True and a["signed_permutation_match"] is None
    assert a["loop_map_matrix"] == [[1, 0], [0, 1]]
    # A's k1 -> B's k1 + p1 reproduces B's lines q_A(k1 + p1): the solved shift
    assert a["shift_matrix"] == [[1, 0, 0], [0, 0, 0]] and a["shift_legs"] == ["p1", "p2", "p3"]
    assert a["loop_relabeling"] == ["k1 -> k1 +p1", "k2 -> k2"]
    assert a["n_shift_systems_solved"] >= 1 and a["shift_in_census_range"] is True
    assert hashes(F)[0] == hashes(Fs)[0]
    # the tampered target: one line's momentum changed (its loop block kept, so the
    # loop-block filter passes and a shift is solved) -- the verification refuses
    props_t = list(props)
    k = [i for i, (e, _m) in enumerate(props_t) if e == shifted("k1-p1-p2")][0]
    props_t[k] = ("k1+p1-p1-p3", props_t[k][1])
    Ft = ta.Family("tampered", F.loops, F.exts, F.ext_subs, props_t, "memory:control", F.physical,
                   leg_virt=F.leg_virt)
    at = ta.affine_iso_search(F, Ft)
    assert at["found"] is False
    assert at["n_U_surviving_loop_block_filter"] >= 1 and at["n_shift_systems_solved"] >= 1
    assert at["tried"] >= 1 and "verifications" in at["reason"]
    assert ta.affine_iso_search(F, Ft, respect_leg_classes=False)["found"] is False
    # the tampered set is no graph at all (the realizer says so by name) and the
    # realized-graph test names the failed side
    rt = ta.realized_graph_iso(F, Ft)
    assert rt["with_masses"] is None and "vertex realization failed for B" in rt["error"]
    assert "NON-GRAPH" in rt["error"] and hashes(Ft)[0] is None
    # the primitive on a target shifted by HALF an external momentum (a set no
    # integer loop redefinition reaches): the shift is solved exactly as the
    # rational 1/2, the set equality holds and the map is reported with the
    # fraction spelled out and flagged outside the census's {-1,0,1} range --
    # never silently rounded, never claimed integral
    vecs, _m, _i = ta.momentum_vectors(F)
    half = [tuple(list(v[:2]) + [v[2] + sp.Rational(1, 2) * v[0]] + list(v[3:])) for v in vecs]
    res = ta._affine_search_vectors(vecs, half, 2, 4, [tuple(range(4))],
                                    dep_rows=ta._dependent_leg_rows(F))
    assert res["found"] is True and res["shift_in_census_range"] is False
    assert res["loop_relabeling"] == ["k1 -> k1 +1/2*p1", "k2 -> k2"]
    assert res["shift_matrix"] == [["1/2", 0, 0], [0, 0, 0]]


def test_isomorphism_to_affine_flag_and_details_by_value():
    smi = ta.CATALOG_SOURCES["planar_smirnov_dbox"]
    fam = ta.Family("s", smi["loops"], smi["exts"], smi["ext_subs"], smi["propagators"], "memory")
    vecs, _m, _i = ta.momentum_vectors(fam)
    dep = ta._dependent_leg_rows(fam)
    # a re-routed copy: k2 -> k2 - k1 (unimodular, no shift)
    props = [(re.sub(r'\bk2\b', "(k2-k1)", e), m) for e, m in smi["propagators"]]
    fam2 = ta.Family("s2", smi["loops"], smi["exts"], smi["ext_subs"], props, "memory")
    vecs2, _m, _i = ta.momentum_vectors(fam2)
    det = {}
    assert ta.isomorphism_to(vecs, vecs2, 2, 4, dep_rows=dep, details=det) is None
    assert det["n_relabelings_tried"] == 8 and "affine" not in det
    det = {}
    rel = ta.isomorphism_to(vecs, vecs2, 2, 4, dep_rows=dep, details=det, affine=True)
    assert rel == "k1 -> k1, k2 -> -k1 +k2" or rel == "k1 -> k1 +k2, k2 -> k2", rel
    assert det["relabeling_class"] == "affine unimodular loop redefinition"
    assert det["affine"]["found"] is True and det["affine"]["n_unimodular"] == 40
    assert "40 unimodular matrices" in det["affine"]["bound"]
    # the signed class stays first: an identity is reported as a signed permutation
    det = {}
    assert ta.isomorphism_to(vecs, vecs, 2, 4, dep_rows=dep, details=det, affine=True) == "identity"
    assert det["relabeling_class"] == "signed loop permutation"
    # find_all: the signed matches first, then the affine ones, classes named
    det = {}
    allm = ta.isomorphism_to(vecs, vecs, 2, 4, dep_rows=dep, details=det, affine=True, find_all=True)
    classes = [m["relabeling_class"] for m in det["matches"]]
    assert classes[0] == "signed loop permutation" and "affine unimodular loop redefinition" in classes
    assert len(allm) == len(det["matches"]) and allm[0] == "identity"
    # the leg-count guard still governs (no affine search across leg counts)
    box = ta.CATALOG_SOURCES["massless_box_1l"]
    fb = ta.Family("b", box["loops"], box["exts"], box["ext_subs"], box["propagators"], "memory")
    vb, _m, _i = ta.momentum_vectors(fb)
    det = {}
    assert ta.isomorphism_to(vecs, vb, 2, 4, details=det, affine=True) is None
    assert "propagator counts differ" in det["reason"]


# ---------------------------------------------------------------------------
#  the fingerprint fields; the hash agrees with VF2 over the fixture pairs
# ---------------------------------------------------------------------------
@nx_required
def test_fingerprint_realized_hash_fields_and_the_non_realizing_family():
    # row 21 read at the census's integral (integrals[2], undotted; the default
    # read takes a dotted integral and the multiplicity is a label of the hash)
    fp = ta.fingerprint(load("row21/family_of_record.json", integral_index=2))
    for k in ("canonical_hash_realized", "canonical_hash_realized_legs", "canonical_hash_realized_method",
              "canonical_hash_realized_reason", "realized_hash_alternatives", "realized_V", "realized_E"):
        assert k in fp, k
    assert fp["canonical_hash"] == "ccb63e84071d" and len(fp["canonical_hash_realized"]) == 12
    assert fp["canonical_hash_realized"] != fp["canonical_hash"]
    assert fp["canonical_hash_realized_legs"] is not None and fp["canonical_hash_realized_reason"] is None
    assert fp["canonical_hash_realized_method"].startswith("individualization-refinement canonical form")
    assert (fp["realized_V"], fp["realized_E"]) == (8, 10)
    fp_dot = ta.fingerprint(load("row21/family_of_record.json"))
    assert fp_dot["nu_multiset"] == [1, 1, 1, 1, 1, 2] and fp_dot["canonical_hash"] == fp["canonical_hash"]
    assert fp_dot["canonical_hash_realized"] != fp["canonical_hash_realized"]      # a dot is a label
    # the row-21 drawn K4 (four legs given): the same realized hash, a different routed one
    D = with_legs(load("row21/drawn_graph.yaml"), (["p1", "p2", "p3", "p4"], {"p4": "-p1-p2-p3"}))
    fd = ta.fingerprint(D)
    assert fd["canonical_hash_realized"] == fp["canonical_hash_realized"]
    assert fd["canonical_hash"] != fp["canonical_hash"]
    # a family that does not realize (row 6 retired, NON-GRAPH): hash None, the reason named
    fn = ta.fingerprint(load("row06/retired_sevenline_vertex2L/family_of_record.json"))
    assert fn["canonical_hash_realized"] is None and fn["canonical_hash_realized_legs"] is None
    assert "realization failed" in fn["canonical_hash_realized_reason"] and "NON-GRAPH" in fn["canonical_hash_realized_reason"]
    assert fn["realized_V"] is None
    # the kite (two distinct leg-labelled realizations): one hash, no alternatives
    fk = ta.fingerprint(load("row14/family_of_record.yaml"))
    assert fk["canonical_hash_realized"] is not None and fk["realized_hash_alternatives"] == []
    # the CLI prints both hashes by name
    import subprocess
    p = subprocess.run([sys.executable, os.path.join(TOOL_DIR, "topology_audit.py"),
                        fixture("row21/family_of_record.json"), "--integral", "2"],
                       capture_output=True, text=True)
    assert p.returncode == 0
    assert "canonical_hash=ccb63e84071d  (routed)" in p.stdout
    assert "canonical_hash_realized=%s" % fp["canonical_hash_realized"] in p.stdout
    js = subprocess.run([sys.executable, os.path.join(TOOL_DIR, "topology_audit.py"),
                         fixture("row21/family_of_record.json"), "--integral", "2", "--json"],
                        capture_output=True, text=True)
    out = json.loads(js.stdout)[0]
    assert out["canonical_hash_realized"] == fp["canonical_hash_realized"]


@nx_required
def test_realized_hash_equal_iff_vf2_with_masses_over_the_fixture_pairs():
    """The property gate of the canonical form: over every pair of fixture
    families of equal loop and leg count, the realized hashes are equal
    exactly when VF2 with (mass, multiplicity) on the lines says isomorphic."""
    fams = {
        "r01F": load("row01/family_of_record.jl"), "r01D": load("row01/drawn_graph.yaml"),
        "r02F": load("row02/family_of_record.jl"), "r02D": load("row02/drawn_graph.yaml"),
        "r07F": load("row07/family_of_record.jl"), "r07D": load("row07/drawn_graph.yaml"),
        "r20F": load("row20/family_of_record.yaml", kin="row20/family_kinematics.yaml"),
        "r20D": load("row20/drawn_graph.yaml", kin="row20/drawn_kinematics.yaml"),
        "r21F": load("row21/family_of_record.json"),
        "r28F": load("row28/family_of_record.json"),
        "r04D": load("row04/drawn_graph.yaml"), "r04B": load("row04/family_of_record_reading_b_4pt.yaml"),
        "r19F": _banana("row19")[0], "r19D": _banana("row19")[1],
        "r31F": _banana("row31")[0], "r31D": _banana("row31")[1],
        "r09F": load("row09/family_of_record.json"), "r09D": load("row09/drawn_graph.yaml"),
        "r12F": load("row12/family_of_record.json"), "r12D114": load("row12/drawn_graph_114.yaml"),
        "r14F": load("row14/family_of_record.yaml"), "r14D": load("row14/drawn_graph.yaml"),
        "r14K": load("row14/family_also/control_kite_tm3_g150_o10.yaml"),
    }
    for cname, spec in ta.CATALOG_SOURCES.items():
        fams["cat:" + cname] = ta.Family(cname, spec["loops"], spec["exts"], spec["ext_subs"],
                                         spec["propagators"], "memory:catalog")
    H = {k: hashes(f)[0] for k, f in fams.items()}
    assert all(h is not None for h in H.values()), [k for k, h in H.items() if h is None]
    n_pairs = n_equal = 0
    names = sorted(fams)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            A, B = fams[a], fams[b]
            if len(A.loops) != len(B.loops) or len(A.exts) != len(B.exts):
                continue
            r = ta.realized_graph_iso(A, B)
            assert r["with_masses"] is not None, (a, b, r["error"])
            assert (H[a] == H[b]) is r["with_masses"], (a, b, H[a], H[b], r["with_masses"])
            n_pairs += 1
            n_equal += 1 if r["with_masses"] else 0
    assert n_pairs >= 60 and n_equal >= 8, (n_pairs, n_equal)


@nx_required
def test_selftest_isomorphism_leg_is_registered_and_catalog_leg8_is_plain():
    assert any(fn is ta._selftest_isomorphism for _lid, _title, fn in ta.EXTRA_SELF_TEST_LEGS)
    r = ta._selftest_catalog(verbose=False)
    assert r["expected_fail"] == 0 and r["ok"] is True
    key = "[8] drawn four-point triple box matches planar_triple_box_3l"
    assert key in r["results"] and r["results"][key]["ok"] is True


# ---------------------------------------------------------------------------
BATTERY = [
    test_unimodular_maps_counts_and_class_by_value,
    test_solve_rational_and_describe_affine_by_value,
    test_canonical_form_is_invariant_under_relabeling_and_separates_a_moved_edge,
    test_rows01_02_07_graph_isomorphic_with_equal_realized_and_different_routed_hashes,
    test_row20_affine_map_equals_the_census_map,
    test_rows19_31_32_34_affine_maps_equal_the_census_maps,
    test_rows09_12_sunrise_maps_equal_the_census_stage_b,
    test_row12_panel_matrix_is_diagonal_under_the_affine_class_too,
    test_row14_half_turn_control_passes_and_the_others_fail_under_the_affine_class,
    test_row30_adj_vs_opp_not_isomorphic_and_each_variant_is_its_own_graphspec,
    test_row04_family_vs_drawn_not_isomorphic_with_the_census_counts,
    test_row04_drawn_matches_planar_triple_box_in_the_audit,
    test_moved_rung_is_not_isomorphic_under_every_U_and_the_search_says_so,
    test_a_different_mass_on_one_line_is_refused_under_every_U,
    test_planted_shift_is_found_and_a_tampered_target_is_refused,
    test_isomorphism_to_affine_flag_and_details_by_value,
    test_fingerprint_realized_hash_fields_and_the_non_realizing_family,
    test_realized_hash_equal_iff_vf2_with_masses_over_the_fixture_pairs,
    test_selftest_isomorphism_leg_is_registered_and_catalog_leg8_is_plain,
]


def run_battery(verbose=True):
    """Run every case; return [{"name", "status": PASS|FAIL|SKIP, "detail"}].
    A missing fixture is SKIP (named); drift or any other exception is FAIL."""
    rows = []
    for fn in BATTERY:
        name = fn.__name__
        try:
            fn()
            rows.append({"name": name, "status": "PASS", "detail": ""})
        except FixtureMissing as exc:
            rows.append({"name": name, "status": "SKIP", "detail": str(exc)})
        except Exception as exc:
            rows.append({"name": name, "status": "FAIL",
                         "detail": f"{type(exc).__name__}: {exc}"})
        if verbose:
            r = rows[-1]
            print(f"   {r['status']:5s} {name}" + (f"  ({r['detail']})" if r["detail"] else ""),
                  flush=True)
    return rows


if __name__ == "__main__":
    out = run_battery(verbose=True)
    n_fail = sum(1 for r in out if r["status"] == "FAIL")
    n_skip = sum(1 for r in out if r["status"] == "SKIP")
    n_pass = sum(1 for r in out if r["status"] == "PASS")
    print(f"isomorphism battery: {len(out)} cases, {n_pass} PASS, {n_fail} FAIL, {n_skip} SKIP")
    sys.exit(1 if n_fail else 0)
