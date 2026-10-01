#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision.
"""
Kinematics battery for topology_audit.py: the external-leg set is read as
DATA (a Kira kinematics.yaml, the yaml's own external_momenta /
momentum_conservation / leg_virtualities keys, the JSON legs + conservation)
and completed by complete_legs (the single-current second leg pB = -p, the
implicit p_N, composite legs P12 / P56 with the orientation resolved from the
propagators); the propagator-symbol fallback is named LEG-SET-INFERRED and
gives no planarity verdict.  Every case reproduces a hand-checked object of
the identity census (IDENTITY_CENSUS.json .rows[].fingerprint,
rows/<row>/audit_family*.json import routes, ROW.md sentences), copied, never
paraphrased; every planted negative control flips or is refused by name.

  row 20  A1 (four-point, p4 unspelled) + family_kinematics.yaml
          -> planar False, cut {s: 2, t: 3, u: 3}, canonical_set_size 7,
             hash 3c695d99a1e1, p4^2 = 7/25; without the kinematics the trap
             by object (size 6, hash 0a0dd74496c3) with LEG-SET-INFERRED and
             NO planarity verdict; the drawn graph + drawn_kinematics.yaml
             runs --leg-perms (rc 0) to planar False, {s: 2, t: 3, u: 3},
             hash 2efe9d563281
  row 19  BAN single current -> legs [p, pB], cut {ppB: 0}, hash
          173439ca7eee, size 4, "vertex degrees [5, 5]"; rows 31 / 32 / 34
          -> 173439ca7eee / 778036c88aae / 778036c88aae
  rows 8-14  sunrise / kite two-point families -> the outgoing leg made
          explicit, two-point cut 3 / 2; the record's own-basis hash
          bad374f9a861 / 6f8bbe224443 kept as canonical_hash_declared_legs
  rows 15-17  icc / iccg composite legs P12, P56 -> P34 = P12 - P56 with
          P12 resolved OUTGOING, hash 8a73188ccb5a, cuts (3, 3, 2), planar
          True with the record's kinematics.yaml (virtualities -2, -1, -3);
          control: a kinematics declaring P34 = P12 + P56 -> NON-GRAPH by name
  row 23  kite_thr single current -> hash bb2ae2aad8f4 (the census's
          normalized hash), two-point cut 2, the AuxLeg kinematics form read
  row 24  retired hexabox + kinematics_of_record.yaml -> five legs, p5
          = -p1-p2-p3-p4, virtualities {p1: mm, p2: 0, p3: mm, p4: 0, p5: 0},
          hash 4417de3002c1, size 8; cured JSON (four legs listed) -> p5
          implicit, {p1: 0, p2: 0, p3: mm, p4: mm, p5: 0}
  rows 1, 2, 4, 7  routed drawings: the census label keys READ (no
          LEG-SET-INFERRED); a copy without them -> LEG-SET-INFERRED
  row 21  the yaml transcription without kinematics -> LEG-SET-INFERRED,
          planar withheld (the bare "planar": true is under
          planarity_on_inferred_leg_set); the JSON's four legs unchanged
  row 28  yaml + the Kira kinematics form -> four legs, hash 34c5ba1eaf82
  row 35  four legs spelled: the verdict is withheld without kinematics
          (uniform rule) and restored by kinematics_of_record.yaml (planar
          False, {s: 2, t: 3, u: 3}, p3^2 = p4^2 = mW2)
  row 18  fam2 (two legs spelled) -> p3 implicit, hash 996d57175c92
  controls  kinematics naming [p, q] for a single-current family -> REFUSED
          (LEG-NOT-SPELLED); row-24 kinematics with the p5 rule deleted ->
          LEG-SET-COMPLETED by name; row-20 kinematics with the off-shell leg
          moved -> virtualities {p3: M, p4: 0}
  NP-dbox  the vendored leg-[1] fixture reproduces the record print
          (b6577ea5d3bf, {s: 2, t: 3, u: 4}, planar, planar_smirnov under
          k2 -> -k2)

Runs under pytest (test_* functions) and as a battery from the tool's
--self-test leg [kinematics] (run_battery()).  A missing fixture is a SKIP by name;
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
    "row15/family_of_record_kira_icc_eqmass.yaml": "9d0e526a98317304",
    "row15/kinematics_of_record.yaml": "961466e57b445b71",
    "row15/family_of_record_amflow_icc_eqmass__e1hp_xm1_g320.json": "0c51e5734673787d",
    "row16/family_of_record_kira_iccg_generic.yaml": "5be2fe7689e8dee6",
    "row17/family_of_record_kira_icc_eqmass.yaml": "9d0e526a98317304",
    "row17/family_of_record_kira_iccg_generic.yaml": "5be2fe7689e8dee6",
    "row18/family_of_record_fam2.yaml": "c1c9c4e3ea6eb708",
    "row19/family_of_record.yaml": "e50056c83e5bbb9b",
    "row19/family_kinematics.yaml": "2b383b50c72a400c",
    "row20/family_of_record.yaml": "d4aef7e8609a1f65",
    "row20/family_kinematics.yaml": "9fb7631fb6a8455a",
    "row20/drawn_graph.yaml": "c756f4b316e8953f",
    "row20/drawn_kinematics.yaml": "b749dd5cf53fee3f",
    "row23/family_of_record_kira_kite_thr_detransport.yaml": "53657fedefba6137",
    "row23/kinematics_of_record.yaml": "5ee9f2717eca1c02",
    "row31/family_of_record.yaml": "4c748aa2bdb4c57d",
    "row31/family_kinematics.yaml": "d39dcfe38892d94a",
    "row32/family_of_record.yaml": "457e82a9a53b2b5f",
    "row32/family_kinematics.yaml": "d39dcfe38892d94a",
    "row34/family_of_record.yaml": "14c51713c29a4f59",
    "row34/family_kinematics.yaml": "2b383b50c72a400c",
    "row35/family_of_record.yaml": "3d3f3825b7b83326",
    "row35/kinematics_of_record.yaml": "63aee0ed75326706",
    "npdbox/family.jl": "d9a7d85ff1248e63",
    "row21/family_of_record.yaml": "1ef7a518c2f513a5",
    "row21/family_of_record.json": "209d2d1741147557",
    "row28/family_of_record.yaml": "f07570cc5983b64d",
    "row28/family_of_record_kira_kinematics.yaml": "512fc1c966e86e04",
    "row24/retired_nonadjacent_p1p3/family_of_record.yaml": "9cb3b9b0d49d4ceb",
    "row24/retired_nonadjacent_p1p3/kinematics_of_record.yaml": "64d62e86c0b095e8",
    "row24/cured_adjacent_p3p4/family_of_record.json": "1864217f1e908a08",
    "row08/family_of_record.yaml": "6081fa14a3906a7d",
    "row08/family_of_record.json": "4109f261715a4cfe",
    "row09/family_of_record.yaml": "3b8cc1d4cc9cffc0",
    "row10/family_of_record.yaml": "4d29f572ffa7a4bb",
    "row11/family_of_record.yaml": "4b90f746ee23b672",
    "row12/family_of_record.yaml": "21d13748eea6a067",
    "row13/family_of_record.yaml": "00fa399594fc002a",
    "row14/family_of_record.yaml": "a6d50193eabd44f3",
    "row01/drawn_graph.yaml": "160f4f486799b5a1",
    "row02/drawn_graph.yaml": "61dd5328d6ffdcfb",
    "row04/drawn_graph.yaml": "87bcbde125b06c81",
    "row07/drawn_graph.yaml": "047a94a5a2be7935",
}

# Expected objects, copied from the census record (IDENTITY_CENSUS.json
# .rows[].fingerprint, rows/<row>/audit_family.json import routes,
# rows/<row>/ROW.md), never paraphrased.
EXPECT = {
    "row20": {"hash": "3c695d99a1e1", "cut": {"s": 2, "t": 3, "u": 3}, "planar": False,
              "size": 7, "virt": {"p1": "0", "p2": "0", "p3": "0", "p4": "7/25"},
              "trap_hash": "0a0dd74496c3", "trap_size": 6,
              "drawn_hash": "2efe9d563281"},
    "row19": {"hash": "173439ca7eee", "cut": {"ppB": 0}, "size": 4, "cli_hash": "7e31a3d20145",
              "degrees": [5, 5], "virt": {"p": "pp", "pB": "pp"}},
    "row31": {"hash": "173439ca7eee", "cut": {"ppB": 0}, "masses": ["1", "1", "1", "9"]},
    "row32": {"hash": "778036c88aae", "cut": {"ppB": 0}, "masses": ["1", "1", "1", "1", "16"]},
    "row34": {"hash": "778036c88aae", "cut": {"ppB": 0}, "masses": ["1", "1", "1", "1", "1"]},
    "sunrise": {"cli_hash": "bad374f9a861", "two_point": 3, "n_lines": 3},
    "kite": {"cli_hash": "6f8bbe224443", "two_point": 2, "n_lines": 5},
    "icc": {"hash": "8a73188ccb5a", "cli_hash": "365550ff3e84", "cuts_in_order": [3, 3, 2],
            "planar": True, "legs": ["P56", "P34", "P12"],
            "virt_record": {"P12": "-2", "P56": "-1", "P34": "-3"}},
    "row23": {"hash": "bb2ae2aad8f4", "cli_hash": "6f8bbe224443", "two_point": 2,
              "cut_value": 0, "masses": ["0", "1", "2", "4", "9"]},
    "row24r": {"legs": ["p1", "p2", "p3", "p4", "p5"], "rule": {"p5": "-p1-p2-p3-p4"},
               "virt": {"p1": "mm", "p2": "0", "p3": "mm", "p4": "0", "p5": "0"},
               "hash": "4417de3002c1", "size": 8, "n_lines": 8},
    "row24c": {"virt": {"p1": "0", "p2": "0", "p3": "mm", "p4": "mm", "p5": "0"},
               "hash": "4417de3002c1", "size": 8},
    "row21": {"trap_hash": "77c2b0ec8a8f", "trap_cut": {"p1p2": 3, "p1p3": 3, "p2p3": 3},
              "trap_size": 6, "json_hash": "ccb63e84071d", "json_cut": {"s": 4, "t": 4, "u": 4}},
    "row28": {"hash": "34c5ba1eaf82", "n_lines": 8, "n_massive": 6},
    "row35": {"hash": "8c11e810a1f1", "cut": {"s": 2, "t": 3, "u": 3}, "planar": False,
              # p4^2 = (p1+p2+p3)^2 through the rules = mW2 (the record's own
              # comment: s + t + u = 2 mW2, two massive W legs)
              "virt": {"p1": "0", "p2": "0", "p3": "mW2", "p4": "mW2"}},
    "row18": {"hash": "996d57175c92", "cut": {"p1p2": 2, "p1p3": 2, "p2p3": 2}, "cli_hash": "fb0d97b7be09"},
    "dbox_drawn": {"cut": {"s": 2, "t": 3, "u": 4}, "planar": True},
    "npdbox": {"hash": "b6577ea5d3bf", "cut": {"s": 2, "t": 3, "u": 4}, "planar": True,
               "n_lines": 7, "n_isp": 2, "n_massive": 0, "relabeling": "k2 -> -k2"},
}

INFERRED = "LEG-SET-INFERRED"
COMPLETED = "LEG-SET-COMPLETED"
RESOLVED = "LEG-ORIENTATION-RESOLVED"


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


def _tmpdir():
    return tempfile.mkdtemp(prefix="topology_audit_kinematics_")


def _write(path, text):
    with open(path, "w") as f:
        f.write(text)
    return path


# ---------------------------------------------------------------------------
#  row 20: the four-point family whose top sector never spells p4 (items 28, 18)
# ---------------------------------------------------------------------------
@nx_required
def test_row20_kinematics_read_as_data_gives_the_import_route():
    e = EXPECT["row20"]
    fam = load("row20/family_of_record.yaml", "row20/family_kinematics.yaml")
    assert fam.exts == ["p1", "p2", "p3", "p4"] and fam.ext_subs == {"p4": "-p1-p2-p3"}
    assert fam.kinematics_source.startswith("kinematics:")
    assert fam.leg_virt == e["virt"], fam.leg_virt
    fp = ta.audit(fam)
    assert fp["planar"] is e["planar"], fp["planarity_method"]
    assert dict(fp["cut_signature"]) == e["cut"], dict(fp["cut_signature"])
    assert fp["canonical_set_size"] == e["size"]
    assert fp["canonical_hash"] == e["hash"], fp["canonical_hash"]
    assert not has(fp, INFERRED)
    assert fp["leg_set"]["source"] == fam.kinematics_source
    assert fp["leg_set"]["virtualities"] == e["virt"]


def test_row20_without_kinematics_is_the_named_trap_with_no_planarity_verdict():
    e = EXPECT["row20"]
    fp = audit("row20/family_of_record.yaml")
    assert fp["leg_set"]["source"] == "inferred"
    assert fp["leg_set"]["legs"] == ["p1", "p2", "p3"]
    assert has(fp, INFERRED), fp["warnings"]
    assert fp["planar"] is None and fp["planarity_method"].startswith("withheld: LEG-SET-INFERRED")
    assert fp["canonical_set_size"] == e["trap_size"]          # D5 collapses onto D4
    assert fp["canonical_hash"] == e["trap_hash"]
    assert "planarity_on_inferred_leg_set" in fp


@nx_required
def test_row20_drawn_graph_with_its_kinematics_runs_leg_perms_to_the_census_numbers():
    e = EXPECT["row20"]
    fam = load("row20/drawn_graph.yaml", "row20/drawn_kinematics.yaml")
    assert fam.exts == ["p1", "p2", "p3", "p4"]
    fp = ta.audit(fam, try_leg_perms=True)          # rc 1 IndexError before: four legs now
    assert fp["planar"] is False
    assert dict(fp["cut_signature"]) == e["cut"]
    assert fp["canonical_hash"] == e["drawn_hash"]
    assert fam.leg_virt["p4"] == "M"


def test_row20_control_moved_offshell_leg_changes_the_virtualities():
    src = open(fixture("row20/family_kinematics.yaml")).read()
    assert '- [[p3, p3], 0]' in src and '- [[p2, p3], "(7/25-s-t)/2"]' in src
    mut = src.replace('- [[p3, p3], 0]', '- [[p3, p3], "M"]').replace(
        '- [[p2, p3], "(7/25-s-t)/2"]', '- [[p2, p3], "(-M-s-t)/2"]')
    tmp = _tmpdir()
    try:
        fam = ta.load_family(fixture("row20/family_of_record.yaml"),
                             kinematics=_write(os.path.join(tmp, "k20.yaml"), mut))[0]
        assert fam.leg_virt == {"p1": "0", "p2": "0", "p3": "M", "p4": "0"}, fam.leg_virt
        rec = load("row20/family_of_record.yaml", "row20/family_kinematics.yaml")
        assert rec.leg_virt != fam.leg_virt
        assert ta.fingerprint(rec)["canonical_hash"] == ta.fingerprint(fam)["canonical_hash"]
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# ---------------------------------------------------------------------------
#  single-current two-point families (items 27, 31)
# ---------------------------------------------------------------------------
@nx_required
def test_row19_single_current_gets_its_outgoing_leg():
    e = EXPECT["row19"]
    for kin in (None, "row19/family_kinematics.yaml"):
        fam = load("row19/family_of_record.yaml", kin)
        assert fam.exts == ["p", "pB"] and fam.ext_subs == {"pB": "-p"}, (fam.exts, fam.ext_subs)
        fp = ta.audit(fam)
        assert dict(fp["cut_signature"]) == e["cut"]
        assert fp["canonical_hash"] == e["hash"]
        assert fp["canonical_set_size"] == e["size"]
        assert fp["canonical_hash_declared_legs"] == e["cli_hash"]
        assert fp["two_point_cut_lines_between_the_legs"] == 4
        assert has(fp, COMPLETED)
        rep = fp["cut_signature"].realization
        assert rep["V"] == 2 and rep["connected"] is True
        assert sorted(rep["vertex_degrees"]) == e["degrees"], rep     # "vertex degrees [5, 5]"
        if kin is None:
            assert has(fp, INFERRED) and fp["planar"] is None
        else:
            assert not has(fp, INFERRED)
            assert fp["planar"] is True and "abstract" in fp["planarity_method"]
            assert fam.leg_virt == e["virt"]


@nx_required
def test_rows31_32_34_bananas_reproduce_the_import_route():
    for row in ("row31", "row32", "row34"):
        e = EXPECT[row]
        fam = load(f"{row}/family_of_record.yaml", f"{row}/family_kinematics.yaml")
        assert fam.exts == ["p", "pB"]
        fp = ta.audit(fam)
        assert fp["canonical_hash"] == e["hash"], (row, fp["canonical_hash"])
        assert dict(fp["cut_signature"]) == e["cut"]
        assert fp["planar"] is True
        assert fp["two_point_cut_lines_between_the_legs"] == len(fam.loops) + 1
        _vecs, masses, _i = ta.momentum_vectors(fam)
        assert sorted(str(m) for m in masses) == sorted(e["masses"]), (row, masses)


@nx_required
def test_rows08_14_two_point_cut_counts_and_the_record_basis_hash():
    for row, key in (("row08", "sunrise"), ("row09", "sunrise"), ("row10", "sunrise"),
                     ("row11", "sunrise"), ("row12", "sunrise"), ("row13", "kite"),
                     ("row14", "kite")):
        e = EXPECT[key]
        fp = audit(f"{row}/family_of_record.yaml")
        assert fp["leg_set"]["legs"] == ["p", "pB"], row
        assert fp["two_point_cut_lines_between_the_legs"] == e["two_point"], (row, fp)
        assert fp["canonical_hash_declared_legs"] == e["cli_hash"], (row, fp["canonical_hash_declared_legs"])
        assert fp["n_genuine_propagators"] == e["n_lines"]
        assert dict(fp["cut_signature"]) == {"ppB": 0}
        assert has(fp, INFERRED) and fp["planar"] is None
    fj = audit("row08/family_of_record.json")               # the JSON's legs are data
    assert fj["planar"] is True and not has(fj, INFERRED)
    assert fj["two_point_cut_lines_between_the_legs"] == 3
    assert fj["leg_set"]["virtualities"] == {"p": "s", "pB": "s"}


# ---------------------------------------------------------------------------
#  composite legs P12 / P56 (item 14) and the threshold kite (AuxLeg form)
# ---------------------------------------------------------------------------
def _check_icc(rel, kin=None):
    e = EXPECT["icc"]
    fam = load(rel, kin)
    assert fam.exts == e["legs"], (rel, fam.exts)
    assert fam.ext_subs == {"P12": "-P56 - P34"}, fam.ext_subs
    ls = fam.record_notes["leg_set"]
    assert ls["orientation"]["outgoing"] == ["P12"]
    assert ls["implicit_leg"] == {"P34": "-P12 - P56"}
    assert "P12" in ls["propagator_rewrite"]
    fp = ta.audit(fam)
    assert fp["canonical_hash"] == e["hash"], (rel, fp["canonical_hash"])
    assert fp["canonical_hash_declared_legs"] == e["cli_hash"]
    assert list(fp["cut_signature"].values()) == e["cuts_in_order"], dict(fp["cut_signature"])
    assert has(fp, RESOLVED) and has(fp, COMPLETED)
    return fam, fp


@nx_required
def test_rows15_16_17_composite_legs_resolved_to_the_normalized_hash():
    for rel in ("row15/family_of_record_kira_icc_eqmass.yaml",
                "row16/family_of_record_kira_iccg_generic.yaml",
                "row17/family_of_record_kira_icc_eqmass.yaml",
                "row17/family_of_record_kira_iccg_generic.yaml"):
        fam, fp = _check_icc(rel)
        assert has(fp, INFERRED) and fp["planar"] is None
    # the JSON route (legs [P12, P56], conservation {}) lands on the same object
    fam, fp = _check_icc("row15/family_of_record_amflow_icc_eqmass__e1hp_xm1_g320.json")
    assert not has(fp, INFERRED) and fp["planar"] is True


@nx_required
def test_row15_record_kinematics_gives_planar_and_the_virtualities():
    e = EXPECT["icc"]
    fam, fp = _check_icc("row15/family_of_record_kira_icc_eqmass.yaml", "row15/kinematics_of_record.yaml")
    assert not has(fp, INFERRED)
    assert fp["planar"] is e["planar"]
    assert fam.leg_virt == e["virt_record"], fam.leg_virt
    assert fam.record_notes["leg_set"]["declared_legs"] == ["P12", "P56"]


@nx_required
def test_row15_control_wrong_sign_rule_is_not_realized_right_sign_is():
    tmp = _tmpdir()
    try:
        out = {}
        for tag, rule in (("right", "P12 - P56"), ("wrong", "P12 + P56")):
            k = _write(os.path.join(tmp, f"k_{tag}.yaml"),
                       "kinematics:\n  incoming_momenta: [P12, P56, P34]\n  outgoing_momenta: []\n"
                       f"  momentum_conservation: [P34, {rule}]\n  scalarproduct_rules:\n"
                       '    - [[P12,P12], "-2"]\n    - [[P12,P56], "0"]\n    - [[P56,P56], "-1"]\n')
            fam = ta.load_family(fixture("row15/family_of_record_kira_icc_eqmass.yaml"), kinematics=k)[0]
            out[tag] = (fam, ta.audit(fam))
        fam, fp = out["right"]
        assert fp["canonical_hash"] == EXPECT["icc"]["hash"] and fp["planar"] is True
        assert fam.leg_virt == EXPECT["icc"]["virt_record"]
        fam, fp = out["wrong"]
        rep = fp["cut_signature"].realization
        assert rep["verdict"] == "NON-GRAPH", rep
        assert "no exact cover" in rep["cause"]
        assert fp["planar"] is None
        assert all(v is None for v in fp["cut_signature"].values())
        assert fp["canonical_hash"] != EXPECT["icc"]["hash"]
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


@nx_required
def test_row23_single_current_and_the_auxleg_kinematics_form():
    e = EXPECT["row23"]
    fp = audit("row23/family_of_record_kira_kite_thr_detransport.yaml")
    assert fp["leg_set"]["legs"] == ["p", "pB"]
    assert fp["canonical_hash"] == e["hash"]
    assert fp["canonical_hash_declared_legs"] == e["cli_hash"]
    assert list(fp["cut_signature"].values()) == [e["cut_value"]]
    assert fp["two_point_cut_lines_between_the_legs"] == e["two_point"]
    assert has(fp, INFERRED)
    fam = load("row23/family_of_record_kira_kite_thr_detransport.yaml", "row23/kinematics_of_record.yaml")
    kin = ta.read_kinematics_yaml(fixture("row23/kinematics_of_record.yaml"))
    assert kin["incoming"] == ["p", "kite_thrAuxLeg"] and kin["conservation"] == {"kite_thrAuxLeg": "-(p)"}
    assert fam.exts == ["p", "pB"] and fam.ext_subs == {"pB": "-p"}, (fam.exts, fam.ext_subs)
    assert "kite_thrAuxLeg" in fam.record_notes["leg_set"]["completion"]
    assert fam.leg_virt == {"p": "s", "pB": "s"}
    fp2 = ta.audit(fam)
    assert fp2["canonical_hash"] == e["hash"] and fp2["planar"] is True and not has(fp2, INFERRED)
    _v, masses, _i = ta.momentum_vectors(fam)
    assert sorted(str(m) for m in masses) == e["masses"]


# ---------------------------------------------------------------------------
#  five-point families, p5 implicit (items 8, 11)
# ---------------------------------------------------------------------------
@nx_required
def test_row24_retired_reads_five_legs_from_its_kinematics():
    e = EXPECT["row24r"]
    fam = load("row24/retired_nonadjacent_p1p3/family_of_record.yaml",
               "row24/retired_nonadjacent_p1p3/kinematics_of_record.yaml")
    assert fam.exts == e["legs"] and fam.ext_subs == e["rule"]
    assert fam.leg_virt == e["virt"], fam.leg_virt
    fp = ta.audit(fam)
    assert fp["canonical_hash"] == e["hash"] and fp["canonical_set_size"] == e["size"]
    assert fp["n_genuine_propagators"] == e["n_lines"]
    assert fp["planar"] is False           # legs closed in order p1..p5 (the infinity test is CAT-21's)
    assert not has(fp, INFERRED)
    assert all(v is not None for v in fp["cut_signature"].values())


def test_row24_control_p5_rule_deleted_is_completed_by_name():
    src = open(fixture("row24/retired_nonadjacent_p1p3/kinematics_of_record.yaml")).read()
    assert "momentum_conservation: [p5, -p1-p2-p3-p4]" in src
    tmp = _tmpdir()
    try:
        k = _write(os.path.join(tmp, "k24.yaml"),
                   src.replace("momentum_conservation: [p5, -p1-p2-p3-p4]", "momentum_conservation: []"))
        fam = ta.load_family(fixture("row24/retired_nonadjacent_p1p3/family_of_record.yaml"), kinematics=k)[0]
        assert fam.exts == EXPECT["row24r"]["legs"]
        assert fam.ext_subs == {"p5": "-p1 - p2 - p3 - p4"}
        fp = ta.audit(fam)
        assert has(fp, COMPLETED), fp["warnings"]
        assert fp["canonical_hash"] == EXPECT["row24r"]["hash"]
        assert fam.leg_virt == EXPECT["row24r"]["virt"]
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


@nx_required
def test_row24_cured_json_completes_p5_with_the_replacement_virtualities():
    e = EXPECT["row24c"]
    fam = load("row24/cured_adjacent_p3p4/family_of_record.json")
    assert fam.exts == ["p1", "p2", "p3", "p4", "p5"]
    assert fam.ext_subs == {"p5": "-p1 - p2 - p3 - p4"}
    assert fam.leg_virt == e["virt"], fam.leg_virt
    fp = ta.audit(fam)
    assert fp["canonical_hash"] == e["hash"] and fp["canonical_set_size"] == e["size"]
    assert has(fp, COMPLETED) and not has(fp, INFERRED)
    assert fp["planar"] is False


# ---------------------------------------------------------------------------
#  routed drawings: the census label keys are data (item 18); the trap rows
# ---------------------------------------------------------------------------
@nx_required
def test_rows01_02_04_07_label_keys_are_read_and_their_absence_is_named():
    e = EXPECT["dbox_drawn"]
    for row in ("row01", "row02", "row07"):
        fam = load(f"{row}/drawn_graph.yaml")
        assert fam.kinematics_source == "yaml:external_momenta+momentum_conservation"
        assert fam.exts == ["p1", "p2", "p3", "p4"] and list(fam.ext_subs) == ["p4"]
        assert fam.leg_virt == {"p1": "0", "p2": "0", "p3": "0", "p4": "0"}
        fp = ta.audit(fam)
        assert not has(fp, INFERRED)
        assert fp["planar"] is e["planar"] and dict(fp["cut_signature"]) == e["cut"], (row, dict(fp["cut_signature"]))
    fam = load("row04/drawn_graph.yaml")
    assert fam.kinematics_source == "yaml:external_momenta+momentum_conservation"
    assert fam.exts == ["p1", "p2", "p3", "p4"]
    # control: the same drawing without its label keys
    src = open(fixture("row01/drawn_graph.yaml")).read()
    stripped = "\n".join(l for l in src.splitlines()
                         if not l.strip().startswith(("external_momenta", "momentum_conservation",
                                                      "leg_virtualities"))) + "\n"
    tmp = _tmpdir()
    try:
        fam = ta.load_family(_write(os.path.join(tmp, "drawn01_nolabels.yaml"), stripped))[0]
        assert fam.kinematics_source == "inferred"
        fp = ta.audit(fam)
        assert has(fp, INFERRED) and fp["planar"] is None
        assert fp["planarity_on_inferred_leg_set"]["planar"] is True     # p4 is spelled here
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


@nx_required
def test_row21_yaml_without_kinematics_is_the_trap_and_the_json_is_not():
    e = EXPECT["row21"]
    fp = audit("row21/family_of_record.yaml")
    assert fp["leg_set"]["legs"] == ["p1", "p2", "p3"]
    assert has(fp, INFERRED) and fp["planar"] is None
    assert fp["planarity_on_inferred_leg_set"]["planar"] is True      # the bare true, now withheld
    assert dict(fp["cut_signature"]) == e["trap_cut"]
    assert fp["canonical_hash"] == e["trap_hash"] and fp["canonical_set_size"] == e["trap_size"]
    fj = audit("row21/family_of_record.json")
    assert fj["leg_set"]["legs"] == ["p1", "p2", "p3", "p4"] and not has(fj, INFERRED)
    assert fj["canonical_hash"] == e["json_hash"] and dict(fj["cut_signature"]) == e["json_cut"]
    assert fj["planar"] is False


def test_row28_yaml_with_the_kira_kinematics_form():
    e = EXPECT["row28"]
    fam = load("row28/family_of_record.yaml", "row28/family_of_record_kira_kinematics.yaml")
    assert fam.exts == ["p1", "p2", "p3", "p4"] and fam.ext_subs == {"p4": "-p1-p2-p3"}
    assert fam.mass_values == {"msq": "1"}            # symbol_to_replace_by_one
    assert fam.leg_virt == {"p1": "0", "p2": "0", "p3": "0", "p4": "0"}
    fp = ta.audit(fam)
    assert fp["canonical_hash"] == e["hash"]
    assert fp["n_genuine_propagators"] == e["n_lines"] and fp["n_massive_lines"] == e["n_massive"]
    assert not has(fp, INFERRED)


@nx_required
def test_row35_four_spelled_legs_verdict_withheld_without_kinematics_restored_with_it():
    # the fallback happens to infer the right set here (all four legs are spelled);
    # the rule is uniform: no kinematics data, no planarity verdict -- the record's
    # kinematics_of_record.yaml gives it back, with the off-shell leg p3^2 = mW2
    e = EXPECT["row35"]
    fp = audit("row35/family_of_record.yaml")
    assert has(fp, INFERRED) and fp["planar"] is None
    assert fp["planarity_on_inferred_leg_set"]["planar"] is e["planar"]
    assert fp["canonical_hash"] == e["hash"] and dict(fp["cut_signature"]) == e["cut"]
    fam = load("row35/family_of_record.yaml", "row35/kinematics_of_record.yaml")
    assert fam.exts == ["p1", "p2", "p3", "p4"] and fam.ext_subs == {"p4": "-p1 - p2 - p3"}
    assert fam.leg_virt == e["virt"], fam.leg_virt
    fp = ta.audit(fam)
    assert not has(fp, INFERRED)
    assert fp["planar"] is e["planar"] and dict(fp["cut_signature"]) == e["cut"]
    assert fp["canonical_hash"] == e["hash"]


@nx_required
def test_row18_two_spelled_legs_complete_to_the_normalized_three_point_object():
    e = EXPECT["row18"]
    fp = audit("row18/family_of_record_fam2.yaml")
    assert fp["leg_set"]["legs"] == ["p1", "p2", "p3"] and fp["leg_set"]["ext_subs"] == {"p3": "-p1 - p2"}
    assert fp["canonical_hash"] == e["hash"] and dict(fp["cut_signature"]) == e["cut"]
    assert fp["canonical_hash_declared_legs"] == e["cli_hash"]
    assert has(fp, INFERRED) and has(fp, COMPLETED) and fp["planar"] is None


# ---------------------------------------------------------------------------
#  refusals, the beside-file convention, the reader forms, the NP-dbox leg
# ---------------------------------------------------------------------------
def test_control_kinematics_naming_two_independent_legs_for_a_single_current_is_refused():
    tmp = _tmpdir()
    try:
        k = _write(os.path.join(tmp, "k_pq.yaml"),
                   "kinematics:\n  incoming_momenta: [p, q]\n  outgoing_momenta: []\n"
                   "  momentum_conservation: []\n  scalarproduct_rules:\n"
                   '    - [[p, p], "pp"]\n    - [[q, q], "qq"]\n    - [[p, q], "pq"]\n')
        try:
            ta.load_family(fixture("row19/family_of_record.yaml"), kinematics=k)
        except ValueError as exc:
            assert "REFUSED: LEG-NOT-SPELLED" in str(exc) and "'q'" in str(exc), str(exc)
        else:
            raise AssertionError("a kinematics naming a leg the family never carries was accepted")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


@nx_required
def test_kinematics_yaml_beside_the_family_file_is_read():
    tmp = _tmpdir()
    try:
        d = os.path.join(tmp, "config")
        os.makedirs(d)
        shutil.copy(fixture("row20/family_of_record.yaml"), os.path.join(d, "integralfamilies.yaml"))
        shutil.copy(fixture("row20/family_kinematics.yaml"), os.path.join(d, "kinematics.yaml"))
        fam = ta.load_family(os.path.join(d, "integralfamilies.yaml"))[0]
        assert fam.kinematics_source == "kinematics.yaml beside the family file"
        assert fam.kinematics_path == os.path.join(d, "kinematics.yaml")
        fp = ta.audit(fam)
        assert fp["canonical_hash"] == EXPECT["row20"]["hash"] and fp["planar"] is False
        assert not has(fp, INFERRED)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_read_kinematics_yaml_forms():
    k = ta.read_kinematics_yaml(fixture("row20/family_kinematics.yaml"))
    assert k["incoming"] == ["p1", "p2", "p3", "p4"] and k["legs"] == ["p1", "p2", "p3", "p4"]
    assert k["conservation"] == {"p4": "-p1-p2-p3"}
    assert k["rules"][("p2", "p3")] == "(7/25-s-t)/2" and k["rules"][("p1", "p1")] == "0"
    assert k["invariants"] == ["s", "t"]
    assert k["virtualities"] == {"p1": "0", "p2": "0", "p3": "0", "p4": "7/25"}
    k = ta.read_kinematics_yaml(fixture("row28/family_of_record_kira_kinematics.yaml"))   # "p1*p1" form
    assert k["rules"][("p1", "p2")] == "s/2" and k["symbol_to_replace_by_one"] == "msq"
    k = ta.read_kinematics_yaml(fixture("row24/retired_nonadjacent_p1p3/kinematics_of_record.yaml"))
    assert k["legs"] == ["p1", "p2", "p3", "p4", "p5"] and k["virtualities"]["p5"] == "0"
    k = ta.read_kinematics_yaml(fixture("row15/kinematics_of_record.yaml"))
    assert k["incoming"] == ["P12", "P56"] and k["conservation"] == {} and k["virtualities"] == {"P12": "-2", "P56": "-1"}


@nx_required
def test_npdbox_vendored_fixture_reproduces_the_leg1_record_print():
    e = EXPECT["npdbox"]
    fam = load("npdbox/family.jl")
    assert fam.exts == ["p1", "p2", "p3", "p4"] and fam.kinematics_source == "jl:exts+ext_subs"
    fp = ta.audit(fam)
    assert fp["canonical_hash"] == e["hash"] and dict(fp["cut_signature"]) == e["cut"]
    assert fp["planar"] is e["planar"]
    assert (fp["n_genuine_propagators"], fp["n_isp_or_dotproduct"], fp["n_massive_lines"]) == \
        (e["n_lines"], e["n_isp"], e["n_massive"])
    assert any(m["family"] == "planar_smirnov_dbox" and m["relabeling"] == e["relabeling"]
               for m in fp["catalog_matches"]), fp["catalog_matches"]
    assert not has(fp, INFERRED)


# ---------------------------------------------------------------------------
#  battery runner (used by the tool's --self-test leg [kinematics] and by __main__)
# ---------------------------------------------------------------------------
BATTERY = [
    test_row20_kinematics_read_as_data_gives_the_import_route,
    test_row20_without_kinematics_is_the_named_trap_with_no_planarity_verdict,
    test_row20_drawn_graph_with_its_kinematics_runs_leg_perms_to_the_census_numbers,
    test_row20_control_moved_offshell_leg_changes_the_virtualities,
    test_row19_single_current_gets_its_outgoing_leg,
    test_rows31_32_34_bananas_reproduce_the_import_route,
    test_rows08_14_two_point_cut_counts_and_the_record_basis_hash,
    test_rows15_16_17_composite_legs_resolved_to_the_normalized_hash,
    test_row15_record_kinematics_gives_planar_and_the_virtualities,
    test_row15_control_wrong_sign_rule_is_not_realized_right_sign_is,
    test_row23_single_current_and_the_auxleg_kinematics_form,
    test_row24_retired_reads_five_legs_from_its_kinematics,
    test_row24_control_p5_rule_deleted_is_completed_by_name,
    test_row24_cured_json_completes_p5_with_the_replacement_virtualities,
    test_rows01_02_04_07_label_keys_are_read_and_their_absence_is_named,
    test_row21_yaml_without_kinematics_is_the_trap_and_the_json_is_not,
    test_row28_yaml_with_the_kira_kinematics_form,
    test_row35_four_spelled_legs_verdict_withheld_without_kinematics_restored_with_it,
    test_row18_two_spelled_legs_complete_to_the_normalized_three_point_object,
    test_control_kinematics_naming_two_independent_legs_for_a_single_current_is_refused,
    test_kinematics_yaml_beside_the_family_file_is_read,
    test_read_kinematics_yaml_forms,
    test_npdbox_vendored_fixture_reproduces_the_leg1_record_print,
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
            print(f"   {r['status']:4s} {name}" + (f"  ({r['detail']})" if r["detail"] else ""),
                  flush=True)
    return rows


if __name__ == "__main__":
    out = run_battery(verbose=True)
    n_fail = sum(1 for r in out if r["status"] == "FAIL")
    n_skip = sum(1 for r in out if r["status"] == "SKIP")
    n_pass = sum(1 for r in out if r["status"] == "PASS")
    print(f"kinematics battery: {len(out)} cases, {n_pass} PASS, {n_fail} FAIL, {n_skip} SKIP")
    sys.exit(1 if n_fail else 0)
