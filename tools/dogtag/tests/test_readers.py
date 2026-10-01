#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision.
"""
Reader battery for topology_audit.py: every record form the readers accept is
exercised on a vendored record fixture (tests/fixtures/<row>/..., the
sha-pinned copy of the identity-census corpus file, PROVENANCE.json beside it
pinning both the original's and the copy's sha256, REFUSED on drift) and must
reproduce the hand-checked object of that census
row -- fingerprint hash, line counts, leg set, top-sector index vector, and the
(momentum, mass) list against the census's own transcription of record.  Each
planted negative control (a named mutation of the record) must flip.

  AMFlow-port JSON (load_amflow_json)
    row 21  lbl3x    hash ccb63e84071d, 6 lines / 4 massive, cuts s=t=u=4,
                     planar False with the file's four legs; == transcription
            controls: one top-sector propagator edited -> hash and labeled set
                      change; "conservation" deleted -> refused by name
    row 28  lbl3se   hash 34c5ba1eaf82, 8 lines / 6 massive
    row 29  lbl3kp   hash 5cf5bd67e373, 9 lines / 7 massive, nu = integrals[0]
                     = [1,1,1,2,1,1,1,1,1,0,0,0,0,0,0]
    row 33  lbl3vp   hash 364b988a1df5, 8 lines / 7 massive, nu[:8] =
                     [1,1,1,2,1,2,1,1], cuts {s:2,t:2,u:4}, planar True
    rows 8-12 sunrise  hash bad374f9a861, 3 lines / 3 massive, mask 7,
                     masses from numeric_values == transcription
    rows 13, 14 kite  hash 6f8bbe224443, 5 lines / 3 massive, mask 31
            controls: the census's control families (other mass sets) differ
                      from the drawn graph in the labeled set while the
                      mass-blind hash agrees
    rows 15-17 icc/iccg  legs [P12, P56] declared, mask 15, masses == transcription
            control: row-16 masses permuted in numeric_values -> labeled set
                     changes, hash does not
    row 18  sr2d2    one leg, no conservation key -> loads; == transcription
    row 23  thrkite_h  hash 6f8bbe224443 on the record's leg set; == transcription
    row 24 cured hexaboxP3P4  four legs listed, conservation {} (p5 implicit,
                     left to the kinematics reader), 8 lines, leg virtualities
    row 6 retired vertex2L   three legs with p3 declared, 7 lines, hash
                     284530eb7b89 (the census's import route)
    bare family block (a hand-written propagator list) == the full config
  pySecDec graph (load_pysecdec_graph)
    row 5   pentagon1L  hash 1e7dddf35ab0, the ten pair cuts, planar True,
                     == family_of_record.yaml (the census's routing)
            control: internal line [3,4] -> [3,5] -> hash changes
  drawn-graph edge lists (load_edge_list + route_edges)
    rows 1, 2, 7 (double boxes, V=6 E=7 L=2) and row 4 (triple box, V=8
                     E=10 L=3) from the '# edge list:' headers == the census's
                     routed propagator lists; every vertex conserves momentum
    row 30 graphspec JSON  V=8 E=10 L=3, 8 massive + 2 massless lines, four
                     legs on the nodes named
            control: a rung moved (a2-b2 -> a1-b3) -> hash changes; a
                     disconnected edge list -> refused by name
  dispatch (load_family)  every fixture reaches its reader by extension /
                     content; the routed yamls stay Kira reads

Runs under pytest (test_* functions) and as a battery from the tool's
--self-test leg [8] (run_battery()).  A missing fixture is a SKIP by name;
a drifted fixture is a refusal (FixtureDrift), never a pass.
"""
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

# sha16 (first 16 hex of sha256) of every vendored fixture's bytes on disk;
# PROVENANCE.json in each row directory carries the full sha256 and size of
# the vendored copy, the census record's own sha256 / size, the census-relative
# source name and the list of rewritten lines.  Both are checked on every access.
PINS = {
    "row01/drawn_graph.yaml": "160f4f486799b5a1",
    "row02/drawn_graph.yaml": "61dd5328d6ffdcfb",
    "row04/drawn_graph.yaml": "87bcbde125b06c81",
    "row05/psd_build_pentagon.py": "70d6e3fc61ffb04a",
    "row05/family_of_record.yaml": "0ec539a5fec61bdf",
    "row05/drawn_graph.yaml": "33984b814f95f701",
    "row06/retired_sevenline_vertex2L/family_of_record.json": "bdfc6ccc26612fff",
    "row06/retired_sevenline_vertex2L/family_of_record_transcribed.yaml": "e759606f59374bd6",
    "row07/drawn_graph.yaml": "047a94a5a2be7935",
    "row07/family_of_record.jl": "59dccd6ed629b5d1",
    "row08/family_of_record.json": "4109f261715a4cfe",
    "row08/family_of_record.yaml": "6081fa14a3906a7d",
    "row08/drawn_graph.yaml": "fa24dcb7064d9cd3",
    "row08/family_also/control_sun112_d2_tm1_g280_o11b.yaml": "a3ff093478731e70",
    "row09/family_of_record.json": "8321e31459aa5e2b",
    "row09/family_of_record.yaml": "3b8cc1d4cc9cffc0",
    "row09/drawn_graph.yaml": "2841e2ae08759216",
    "row10/family_of_record.json": "adb680a25b1132ac",
    "row10/family_of_record.yaml": "4d29f572ffa7a4bb",
    "row10/drawn_graph.yaml": "2ad28d93667e43ec",
    "row11/family_of_record.json": "a7acd63a397d3e0f",
    "row11/family_of_record.yaml": "4b90f746ee23b672",
    "row11/drawn_graph.yaml": "3d676688bb7f44c7",
    "row12/family_of_record.json": "2f4f811368032eb8",
    "row12/family_of_record.yaml": "21d13748eea6a067",
    "row13/family_of_record.json": "3acdcb030faec81f",
    "row13/family_of_record.yaml": "00fa399594fc002a",
    "row13/drawn_graph.yaml": "f2e0a018ce4e788f",
    "row13/family_also/control_integralfamilies.yaml": "a6d50193eabd44f3",
    "row14/family_of_record.yaml": "a6d50193eabd44f3",
    "row14/drawn_graph.yaml": "e9c6eb631a18d650",
    "row14/family_also/kira14_sm2.json": "fe134467a6731fe5",
    "row14/family_also/kira14_sm2.yaml": "07a68b4b83ee5175",
    "row14/family_also/control_kite_tm3_g150_o10.yaml": "00fa399594fc002a",
    "row15/family_of_record_amflow_icc_eqmass__e1hp_xm1_g320.json": "0c51e5734673787d",
    "row15/family_of_record_amflow_icc_eqmass.yaml": "3f2ac442c8426abf",
    "row16/family_of_record_amflow_iccg_generic__ge1_xm2.json": "87d4becb98d7d813",
    "row16/family_of_record_amflow_iccg_generic.yaml": "1a3f0921ea6ef2f6",
    "row16/drawn_graph.yaml": "0d09c246ad198687",
    "row17/family_of_record_amflow_icc_eqmass__e1hp_xm1_g320.json": "0c51e5734673787d",
    "row17/family_of_record_amflow_icc_eqmass.yaml": "3f2ac442c8426abf",
    "row17/family_of_record_amflow_iccg_generic__ge1_xm2.json": "87d4becb98d7d813",
    "row17/family_of_record_amflow_iccg_generic.yaml": "1a3f0921ea6ef2f6",
    "row18/family_of_record_amflow_sr2d2_sunrise_subcurve__uncut_d2_in_0.json": "6176cf1243f0b3c0",
    "row18/family_of_record_amflow_sr2d2_sunrise_subcurve.yaml": "8b3206ec548a41fb",
    "row19/family_of_record.yaml": "e50056c83e5bbb9b",
    "row21/family_of_record.json": "209d2d1741147557",
    "row21/family_of_record.yaml": "1ef7a518c2f513a5",
    "row21/drawn_graph.yaml": "c307f8008add78e6",
    "row23/family_of_record_amflow_thrkite_h_oracle__thrkite_h_tm2.json": "d7f51e4e54aa6c74",
    "row23/family_of_record_amflow_thrkite_h_oracle.yaml": "779159f8f8705c16",
    "row24/cured_adjacent_p3p4/family_of_record.json": "1864217f1e908a08",
    "row24/cured_adjacent_p3p4/family_of_record_transcribed.yaml": "ade0416796449ea6",
    "row28/family_of_record.json": "61c4e46fb523550b",
    "row28/family_of_record.yaml": "f07570cc5983b64d",
    "row28/drawn_graph.yaml": "d69e033918b66544",
    "row29/family_of_record.json": "768fd40d15f65236",
    "row29/family_of_record.yaml": "12c82770f0ea5471",
    "row29/drawn_graph.yaml": "440576817cc5bd06",
    "row30/lbl3m_adj/graphspec_adj_parked_copy.json": "eb021170bf42347c",
    "row30/lbl3m_adj/drawn_graph.yaml": "1d1e8738e19d6de9",
    "row30/lbl3m_adj/family_of_record.yaml": "4f1978a049287764",
    "row33/family_of_record.json": "460475b5991282e8",
    "row33/family_of_record.yaml": "7a9c22ebf36bd900",
    "row33/drawn_graph.yaml": "484ccb7f9657f612",
}

# Expected objects, copied from the census record (IDENTITY_CENSUS.json .rows,
# rows/<row>/ROW.md, ROW_RECEIPT.json, iso.json), never paraphrased.
EXPECT = {
    "row21": {"hash": "ccb63e84071d", "n_lines": 6, "n_massive": 4,
              "cut": {"s": 4, "t": 4, "u": 4}, "planar": False,
              "legs": ["p1", "p2", "p3", "p4"], "conservation": {"p4": "-p1-p2-p3"}},
    "row28": {"hash": "34c5ba1eaf82", "n_lines": 8, "n_massive": 6},
    "row29": {"hash": "5cf5bd67e373", "n_lines": 9, "n_massive": 7,
              "nu": [1, 1, 1, 2, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0],
              "cut": {"s": 2, "t": 2, "u": 4}, "planar": True},
    "row33": {"hash": "364b988a1df5", "n_lines": 8, "n_massive": 7,
              "nu8": [1, 1, 1, 2, 1, 2, 1, 1],
              "cut": {"s": 2, "t": 2, "u": 4}, "planar": True},
    "sunrise": {"hash": "bad374f9a861", "n_lines": 3, "n_massive": 3, "mask": 7,
                "nu": [1, 1, 1, 0, 0]},
    "kite": {"hash": "6f8bbe224443", "n_lines": 5, "n_massive": 3, "mask": 31,
             "nu": [1, 1, 1, 1, 1]},
    "icc": {"legs": ["P12", "P56"], "mask": 15, "nu": [1, 1, 1, 1, 0, 0, 0],
            "n_lines": 4, "n_massive": 4, "hash_record_legs": "365550ff3e84"},
    "row23": {"hash_record_legs": "6f8bbe224443", "n_lines": 5, "n_massive": 4,
              "masses": ["0", "1", "2", "4", "9"]},
    "row05": {"hash": "1e7dddf35ab0", "n_lines": 5, "n_massive": 0, "planar": True,
              "cut": {"p1p2": 2, "p1p3": 4, "p1p4": 4, "p1p5": 2, "p2p3": 2,
                      "p2p4": 4, "p2p5": 4, "p3p4": 2, "p3p5": 4, "p4p5": 2}},
    "row24c": {"legs": ["p1", "p2", "p3", "p4"], "n_lines": 8,
               "leg_virt": {"p1": "0", "p2": "0", "p3": "mm", "p4": "mm"}},
    "row06r": {"hash": "284530eb7b89", "n_lines": 7, "legs": ["p1", "p2", "p3"]},
    "row30": {"V": 8, "E": 10, "L": 3, "n_massive": 8, "legs": ["p1", "p2", "p3", "p4"]},
    "dbox_drawn": {"V": 6, "E": 7, "L": 2},
    "ladder_drawn": {"V": 8, "E": 10, "L": 3},
}

try:
    import pytest
except ImportError:                      # direct run without pytest
    pytest = None


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
    # PROVENANCE.json is one canonical shape for the whole package (the shape
    # the self-test's fixture-pin check, _selftest_realizer's drift refusal and
    # tests/census_battery.py read): "files" is a list of {"path", "source",
    # "sha256", "size"}.
    entry = None
    for e in prov.get("files", []) or []:
        if e.get("path") == within:
            entry = e
    if entry is None or entry.get("sha256") != sha:
        raise FixtureDrift(f"fixture drift: {rel} not recorded with sha256 {sha} "
                           f"in {row}/PROVENANCE.json")
    return path


def labeled_set(fam, exts=None, ext_subs=None, mass_values=None):
    """Sorted (canonical momentum vector, mass^2) list of the top-sector
    propagators; exts/ext_subs override the family's leg set (the census
    normalises both sides to one leg set before comparing); mass_values are
    applied on top of the family's own."""
    if exts is not None:
        fam = ta.Family(fam.name, fam.loops, exts, ext_subs or {}, fam.propagators,
                        fam.source, fam.physical, mass_values=fam.mass_values)
    syms = ta._symbols(fam.loops, fam.exts)
    res = ta._resolve_ext_subs(fam.exts, fam.ext_subs, syms)
    subs = {sp.Symbol(k): sp.nsimplify(v) for k, v in (mass_values or {}).items()}
    out = []
    for i, (expr, mass) in enumerate(fam.propagators):
        if i < len(fam.physical) and not fam.physical[i]:
            continue
        parsed = ta._parse_momentum(expr, fam.loops, fam.exts, syms, res, mass)
        if parsed is None:
            out.append((("ISP", expr), ""))
            continue
        vec, _isprop, m = parsed
        if all(c == 0 for c in vec[:len(fam.loops)]):
            out.append((("NO_LOOP_MOMENTUM", expr), ""))
            continue
        m = sp.nsimplify(fam.resolve_mass(m).subs(subs)) if subs else fam.resolve_mass(m)
        out.append((tuple(str(c) for c in ta._canon_vec(vec)), str(m)))
    return sorted(out)


def fp(fam):
    return ta.fingerprint(fam)


def _tmpdir():
    return tempfile.mkdtemp(prefix="topology_audit_readers_")


def _write_json(d, path):
    with open(path, "w") as f:
        json.dump(d, f, indent=1)
    return path


# ---------------------------------------------------------------------------
#  AMFlow-port JSON configs
# ---------------------------------------------------------------------------
@nx_required
def test_row21_json_reads_four_legs_and_import_route_fingerprint():
    fam = ta.load_family(fixture("row21/family_of_record.json"))[0]
    e = EXPECT["row21"]
    assert fam.source_kind == "amflow_json"
    assert fam.kinematics_source == "json:legs+conservation"
    assert fam.exts == e["legs"], fam.exts
    assert {k: v.replace(" ", "") for k, v in fam.ext_subs.items()} == e["conservation"]
    d = json.load(open(fixture("row21/family_of_record.json")))
    assert fam.nu == d["integrals"][0]["indices"]
    assert fam.record_notes["integral_index"] == 0
    assert fam.leg_virt == {"p1": "0", "p2": "0", "p3": "0", "p4": "0"}
    f = fp(fam)
    assert f["canonical_hash"] == e["hash"], f["canonical_hash"]
    assert f["n_genuine_propagators"] == e["n_lines"]
    assert f["n_massive_lines"] == e["n_massive"]
    assert f["cut_signature"] == e["cut"], f["cut_signature"]
    assert f["planar"] is e["planar"]
    # (momentum, mass) list == the census transcription, on the four-leg set
    y = ta.load_family(fixture("row21/family_of_record.yaml"))[0]
    assert labeled_set(fam) == labeled_set(y, e["legs"], {"p4": "-p1-p2-p3"})


def test_row21_control_edited_propagator_changes_hash_and_labeled_set():
    src = fixture("row21/family_of_record.json")
    d = json.load(open(src))
    i = d["family"]["propagators"].index("(k+a-p1)^2 - 1")
    d["family"]["propagators"][i] = "(k+a-p2)^2 - 1"
    tmp = _tmpdir()
    try:
        mut = ta.load_family(_write_json(d, os.path.join(tmp, "row21_edited.json")))[0]
        rec = ta.load_family(src)[0]
        assert fp(mut)["canonical_hash"] != EXPECT["row21"]["hash"]
        assert fp(mut)["n_genuine_propagators"] == 6
        assert labeled_set(mut) != labeled_set(rec)
        assert fp(rec)["canonical_hash"] == EXPECT["row21"]["hash"]   # positive control
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_row21_control_deleted_conservation_is_refused_by_name():
    src = fixture("row21/family_of_record.json")
    d = json.load(open(src))
    del d["family"]["conservation"]
    tmp = _tmpdir()
    try:
        path = _write_json(d, os.path.join(tmp, "row21_noconservation.json"))
        try:
            ta.load_family(path)
        except ValueError as exc:
            assert "conservation" in str(exc) and "REFUSED" in str(exc), str(exc)
        else:
            raise AssertionError("a four-leg family without 'conservation' loaded silently")
        ta.load_family(src)      # positive control: the record loads
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def _check_lbl3_row(row, e):
    fam = ta.load_family(fixture(f"{row}/family_of_record.json"))[0]
    assert fam.exts == ["p1", "p2", "p3", "p4"]
    assert list(fam.ext_subs) == ["p4"]
    f = fp(fam)
    assert f["canonical_hash"] == e["hash"], (row, f["canonical_hash"])
    assert f["n_genuine_propagators"] == e["n_lines"]
    assert f["n_massive_lines"] == e["n_massive"]
    if "cut" in e:
        assert f["cut_signature"] == e["cut"], (row, f["cut_signature"])
        assert f["planar"] is e["planar"]
    y = ta.load_family(fixture(f"{row}/family_of_record.yaml"))[0]
    assert labeled_set(fam) == labeled_set(y, fam.exts, fam.ext_subs, fam.mass_values)
    return fam


def test_row28_json_import_route_fingerprint():
    fam = _check_lbl3_row("row28", EXPECT["row28"])
    assert fam.mass_values == {"msq": "1"}


@nx_required
def test_row29_json_import_route_fingerprint_and_nu():
    fam = _check_lbl3_row("row29", EXPECT["row29"])
    assert fam.nu == EXPECT["row29"]["nu"], fam.nu
    assert fam.record_notes["integral_index"] == 0


@nx_required
def test_row33_json_import_route_fingerprint_and_nu():
    fam = _check_lbl3_row("row33", EXPECT["row33"])
    assert fam.nu[:8] == EXPECT["row33"]["nu8"], fam.nu


def _check_two_point_json(json_rel, yaml_rel, e, extra_yaml=()):
    # the reader's object: the record's own leg set (one current p); the hash
    # quoted by the census is on that basis.  load_family then completes the
    # outgoing leg (pB = -p, the kinematics wave) and keeps the record-basis
    # hash as canonical_hash_declared_legs.
    raw = ta.load_amflow_json(fixture(json_rel))[0]
    assert raw.exts == ["p"] and raw.ext_subs == {}
    assert raw.nu == e["nu"], raw.nu
    assert sum(1 << i for i, ph in enumerate(raw.physical) if ph) == e["mask"]
    f = fp(raw)
    assert f["canonical_hash"] == e["hash"], (json_rel, f["canonical_hash"])
    assert f["n_genuine_propagators"] == e["n_lines"]
    assert f["n_massive_lines"] == e["n_massive"]
    fam = ta.load_family(fixture(json_rel))[0]
    assert fam.exts == ["p", "pB"] and fam.ext_subs == {"pB": "-p"}
    assert ta.audit(fam)["canonical_hash_declared_legs"] == e["hash"]
    y = ta.load_family(fixture(yaml_rel))[0]
    assert labeled_set(fam) == labeled_set(y), (labeled_set(fam), labeled_set(y))
    for rel in extra_yaml:
        z = ta.load_family(fixture(rel))[0]
        assert labeled_set(fam) == labeled_set(z), (rel, labeled_set(fam), labeled_set(z))
    return fam


def test_rows08_to_12_sunrise_json_reproduce_transcription():
    seen = {}
    for row in ("row08", "row09", "row10", "row11", "row12"):
        fam = _check_two_point_json(f"{row}/family_of_record.json",
                                    f"{row}/family_of_record.yaml", EXPECT["sunrise"])
        seen[row] = sorted(m for _v, m in labeled_set(fam))
    # the mass sets the census read off the four families (ROW.md rows 8-12)
    assert seen["row08"] == ["1", "1", "1"]
    assert seen["row09"] == ["1", "1", "2"]
    assert seen["row10"] == ["1", "2", "3"]
    assert seen["row11"] == ["1", "1", "4"]
    assert seen["row12"] == ["1", "1", "2"]


def test_row13_kite_json_reproduces_transcription():
    fam = _check_two_point_json("row13/family_of_record.json",
                                "row13/family_of_record.yaml", EXPECT["kite"])
    assert fam.mass_values == {"msq": "1"}


def test_row14_kite112_json_reproduces_transcription_and_record():
    fam = _check_two_point_json("row14/family_also/kira14_sm2.json",
                                "row14/family_also/kira14_sm2.yaml", EXPECT["kite"],
                                extra_yaml=("row14/family_of_record.yaml",
                                            "row14/drawn_graph.yaml"))
    assert fam.mass_values == {"msq1": "1", "msq3": "1", "msq5": "2"}


def test_rows08_13_14_planted_mass_controls_flip_with_labels_only():
    # the census's own controls: another mass set on the same lines must NOT
    # match the drawn graph once masses are carried, while the mass-blind
    # canonical hash cannot tell them apart
    cases = [("row08/family_of_record.json", "row08/drawn_graph.yaml",
              "row08/family_also/control_sun112_d2_tm1_g280_o11b.yaml"),
             ("row13/family_of_record.json", "row13/drawn_graph.yaml",
              "row13/family_also/control_integralfamilies.yaml"),
             ("row14/family_also/kira14_sm2.json", "row14/drawn_graph.yaml",
              "row14/family_also/control_kite_tm3_g150_o10.yaml")]
    for rec_rel, drawn_rel, ctl_rel in cases:
        rec = ta.load_family(fixture(rec_rel))[0]
        drawn = ta.load_family(fixture(drawn_rel))[0]
        ctl = ta.load_family(fixture(ctl_rel))[0]
        assert labeled_set(rec) == labeled_set(drawn), rec_rel
        assert labeled_set(ctl) != labeled_set(drawn), ctl_rel
        assert fp(ctl)["canonical_hash"] == fp(rec)["canonical_hash"], ctl_rel


def _check_icc(json_rel, yaml_rel, masses):
    raw = ta.load_amflow_json(fixture(json_rel))[0]
    e = EXPECT["icc"]
    assert raw.exts == e["legs"] and raw.ext_subs == {}
    assert raw.nu == e["nu"], raw.nu
    assert sum(1 << i for i, ph in enumerate(raw.physical) if ph) == e["mask"]
    assert raw.leg_virt == {"P12": "p12sq", "P56": "p56sq"}
    f = fp(raw)
    assert f["n_genuine_propagators"] == e["n_lines"]
    assert f["n_massive_lines"] == e["n_massive"]
    assert f["canonical_hash"] == e["hash_record_legs"]
    fam = ta.load_family(fixture(json_rel))[0]        # completed: P34 implicit, P12 outgoing
    assert fam.exts == ["P56", "P34", "P12"]
    assert ta.audit(fam)["canonical_hash_declared_legs"] == e["hash_record_legs"]
    y = ta.load_family(fixture(yaml_rel))[0]
    assert labeled_set(fam) == labeled_set(y, mass_values=fam.mass_values)
    assert sorted(m for _v, m in labeled_set(fam)) == masses
    return fam


@nx_required
def test_rows15_16_17_icc_json_reproduce_transcriptions():
    for row in ("row15", "row17"):
        _check_icc(f"{row}/family_of_record_amflow_icc_eqmass__e1hp_xm1_g320.json",
                   f"{row}/family_of_record_amflow_icc_eqmass.yaml", ["1", "1", "1", "1"])
    for row in ("row16", "row17"):
        fam = _check_icc(f"{row}/family_of_record_amflow_iccg_generic__ge1_xm2.json",
                         f"{row}/family_of_record_amflow_iccg_generic.yaml",
                         ["1", "2", "3", "5"])
        assert fam.mass_values == {"m1sq": "1", "m2sq": "2", "m3sq": "3", "m4sq": "5"}


def test_row16_control_permuted_masses_flip_labeled_set_not_hash():
    src = fixture("row16/family_of_record_amflow_iccg_generic__ge1_xm2.json")
    d = json.load(open(src))
    nv = d["amf_options"]["blackbox"]["numeric_values"]
    nv["m1sq"], nv["m2sq"] = "2", "1"
    tmp = _tmpdir()
    try:
        mut = ta.load_family(_write_json(d, os.path.join(tmp, "row16_permuted.json")))[0]
        rec = ta.load_family(src)[0]
        assert sorted(m for _v, m in labeled_set(mut)) == sorted(m for _v, m in labeled_set(rec))
        assert labeled_set(mut) != labeled_set(rec)
        assert fp(mut)["canonical_hash"] == fp(rec)["canonical_hash"]
        y = ta.load_family(fixture("row16/family_of_record_amflow_iccg_generic.yaml"))[0]
        assert labeled_set(mut) != labeled_set(y, mass_values=mut.mass_values)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_row18_single_leg_without_conservation_key_loads():
    raw = ta.load_amflow_json(fixture(
        "row18/family_of_record_amflow_sr2d2_sunrise_subcurve__uncut_d2_in_0.json"))[0]
    assert raw.exts == ["p"] and raw.ext_subs == {}
    fam = ta.load_family(fixture(
        "row18/family_of_record_amflow_sr2d2_sunrise_subcurve__uncut_d2_in_0.json"))[0]
    assert fam.exts == ["p", "pB"] and fam.ext_subs == {"pB": "-p"}
    assert fam.nu == [1, 1, 1, 0, 0]
    assert fam.mass_values == {"Msq": "9/10"}
    y = ta.load_family(fixture("row18/family_of_record_amflow_sr2d2_sunrise_subcurve.yaml"))[0]
    assert labeled_set(fam) == labeled_set(y, mass_values=fam.mass_values)
    assert sorted(m for _v, m in labeled_set(fam)) == ["1", "1", "9/10"]


def test_row23_threshold_kite_json_reproduces_transcription():
    raw = ta.load_amflow_json(fixture(
        "row23/family_of_record_amflow_thrkite_h_oracle__thrkite_h_tm2.json"))[0]
    assert raw.exts == ["p"] and raw.nu == [1, 1, 1, 1, 1]
    fam = ta.load_family(fixture(
        "row23/family_of_record_amflow_thrkite_h_oracle__thrkite_h_tm2.json"))[0]
    e = EXPECT["row23"]
    assert fam.exts == ["p", "pB"] and fam.nu == [1, 1, 1, 1, 1]
    f = fp(raw)
    assert f["canonical_hash"] == e["hash_record_legs"]
    assert ta.audit(fam)["canonical_hash_declared_legs"] == e["hash_record_legs"]
    assert f["n_genuine_propagators"] == e["n_lines"]
    assert f["n_massive_lines"] == e["n_massive"]
    assert sorted(m for _v, m in labeled_set(fam)) == e["masses"]
    y = ta.load_family(fixture("row23/family_of_record_amflow_thrkite_h_oracle.yaml"))[0]
    assert labeled_set(fam) == labeled_set(y, mass_values=fam.mass_values)


def test_row24_cured_json_reads_four_listed_legs_with_p5_implicit():
    raw = ta.load_amflow_json(fixture("row24/cured_adjacent_p3p4/family_of_record.json"))[0]
    e = EXPECT["row24c"]
    assert raw.exts == e["legs"]
    assert raw.ext_subs == {}                 # conservation declared empty
    assert raw.leg_virt == e["leg_virt"]
    fam = ta.load_family(fixture("row24/cured_adjacent_p3p4/family_of_record.json"))[0]
    assert fam.exts == e["legs"] + ["p5"]     # p5 implicit, completed by the kinematics wave
    assert list(fam.ext_subs) == ["p5"]
    assert fam.mass_values == {}              # every numeric value is kinematic
    assert sum(fam.physical) == e["n_lines"]
    assert fp(fam)["n_genuine_propagators"] == e["n_lines"]
    y = ta.load_family(fixture("row24/cured_adjacent_p3p4/family_of_record_transcribed.yaml"))[0]
    assert labeled_set(fam) == labeled_set(y, fam.exts, fam.ext_subs)


def test_row06_retired_json_reads_three_declared_legs():
    fam = ta.load_family(fixture("row06/retired_sevenline_vertex2L/family_of_record.json"))[0]
    e = EXPECT["row06r"]
    assert fam.exts == e["legs"]
    assert list(fam.ext_subs) == ["p3"]
    assert fam.leg_virt == {"p1": "0", "p2": "0", "p3": "s"}
    f = fp(fam)
    assert f["n_genuine_propagators"] == e["n_lines"]
    assert f["canonical_hash"] == e["hash"]
    y = ta.load_family(fixture(
        "row06/retired_sevenline_vertex2L/family_of_record_transcribed.yaml"))[0]
    assert labeled_set(fam) == labeled_set(y, fam.exts, fam.ext_subs)


def test_bare_family_block_reads_as_a_propagator_list():
    src = fixture("row21/family_of_record.json")
    d = json.load(open(src))
    block = dict(d["family"])
    block["indices"] = d["integrals"][0]["indices"]
    tmp = _tmpdir()
    try:
        bare = ta.load_family(_write_json(block, os.path.join(tmp, "row21_block.json")))[0]
        full = ta.load_family(src)[0]
        assert bare.exts == full.exts and bare.nu == full.nu
        assert labeled_set(bare) == labeled_set(full)
        assert fp(bare)["canonical_hash"] == fp(full)["canonical_hash"]
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_json_integral_index_override():
    src = fixture("row21/family_of_record.json")
    fam = ta.load_family(src, integral_index=2)[0]
    assert fam.nu == [1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    assert fp(fam)["canonical_hash"] == EXPECT["row21"]["hash"]
    try:
        ta.load_family(src, integral_index=7)
    except ValueError as exc:
        assert "REFUSED" in str(exc)
    else:
        raise AssertionError("an out-of-range integral index was accepted")


# ---------------------------------------------------------------------------
#  pySecDec graph
# ---------------------------------------------------------------------------
@nx_required
def test_row05_pysecdec_graph_routes_to_the_record_fingerprint():
    fam = ta.load_family(fixture("row05/psd_build_pentagon.py"))[0]
    e = EXPECT["row05"]
    assert fam.source_kind == "pysecdec_graph"
    assert fam.name == "pentagon1L"
    assert fam.loops == ["k1"]
    assert fam.exts == ["p1", "p2", "p3", "p4", "p5"]
    assert fam.ext_subs == {"p5": "-p1 - p2 - p3 - p4"}
    assert fam.leg_virt == {p: "0" for p in fam.exts}
    assert fam.nu == [1, 1, 1, 1, 1]
    assert all(v == "0" for v in fam.record_notes["vertex_conservation"].values())
    f = fp(fam)
    assert f["canonical_hash"] == e["hash"], f["canonical_hash"]
    assert f["n_genuine_propagators"] == e["n_lines"]
    assert f["n_massive_lines"] == e["n_massive"]
    assert f["planar"] is e["planar"]
    assert f["cut_signature"] == e["cut"], f["cut_signature"]
    y = ta.load_family(fixture("row05/family_of_record.yaml"))[0]
    assert labeled_set(fam) == labeled_set(y)
    d = ta.load_family(fixture("row05/drawn_graph.yaml"))[0]
    assert labeled_set(fam) == labeled_set(d)


def test_row05_control_moved_internal_line_changes_hash():
    src = fixture("row05/psd_build_pentagon.py")
    txt = open(src).read()
    assert "['0',[3,4]]" in txt
    tmp = _tmpdir()
    try:
        path = os.path.join(tmp, "psd_build_pentagon_moved.py")
        open(path, "w").write(txt.replace("['0',[3,4]]", "['0',[3,5]]"))
        mut = ta.load_family(path)[0]
        assert fp(mut)["canonical_hash"] != EXPECT["row05"]["hash"]
        assert labeled_set(mut) != labeled_set(ta.load_family(src)[0])
        assert 2 in mut.record_notes["degrees"].values()   # vertex 4 lost a line
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# ---------------------------------------------------------------------------
#  drawn-graph edge lists
# ---------------------------------------------------------------------------
def _check_drawn_header(rel, e):
    routed = ta.load_edge_list(fixture(rel))[0]
    kira = ta.load_family(fixture(rel))[0]          # the routed yaml stays a Kira read
    assert kira.source_kind == "kira_yaml"
    assert routed.source_kind == "edge_list"
    assert routed.exts == ["p1", "p2", "p3", "p4"]
    assert list(routed.ext_subs) == ["p4"]
    assert routed.leg_virt == {"p1": "0", "p2": "0", "p3": "0", "p4": "0"}
    n = routed.record_notes
    assert (n["V"], n["E"], n["L"]) == (e["V"], e["E"], e["L"]), (rel, n)
    assert all(v == "0" for v in n["vertex_conservation"].values())
    assert labeled_set(routed) == labeled_set(kira), (rel, labeled_set(routed), labeled_set(kira))
    assert fp(routed)["canonical_hash"] == fp(kira)["canonical_hash"]
    return routed


def test_rows01_02_07_drawn_double_boxes_route_as_the_census_did():
    masses = {}
    for row in ("row01", "row02", "row07"):
        fam = _check_drawn_header(f"{row}/drawn_graph.yaml", EXPECT["dbox_drawn"])
        masses[row] = sorted(m for _v, m in labeled_set(fam))
    assert masses["row01"] == ["0"] + ["msq"] * 6      # six massive perimeter lines
    assert masses["row02"] == ["0"] * 7
    assert masses["row07"] == ["0"] * 6 + ["msq"]      # the massive rung


def test_row04_drawn_triple_box_routes_as_the_census_did():
    _check_drawn_header("row04/drawn_graph.yaml", EXPECT["ladder_drawn"])


def test_row30_graphspec_json_routes():
    fam = ta.load_family(fixture("row30/lbl3m_adj/graphspec_adj_parked_copy.json"))[0]
    e = EXPECT["row30"]
    assert fam.source_kind == "edge_list"
    assert fam.name == "lbl3m_adj"
    n = fam.record_notes
    assert (n["V"], n["E"], n["L"]) == (e["V"], e["E"], e["L"])
    assert fam.exts == e["legs"]
    assert n["legs_at"] == {"p1": "1", "p2": "5", "p3": "7", "p4": "8"}
    assert fam.leg_virt == {p: "0" for p in e["legs"]}
    assert all(v == "0" for v in n["vertex_conservation"].values())
    f = fp(fam)
    assert f["n_genuine_propagators"] == e["E"]
    assert f["n_massive_lines"] == e["n_massive"]


def test_drawn_controls_moved_rung_and_disconnected_list():
    src = fixture("row02/drawn_graph.yaml")
    txt = open(src).read()
    assert "a2-b2[0]" in txt
    tmp = _tmpdir()
    try:
        moved = os.path.join(tmp, "drawn_row02_rung_moved.yaml")
        open(moved, "w").write(txt.replace("a2-b2[0]", "a1-b3[0]"))
        mut = ta.load_edge_list(moved)[0]
        rec = ta.load_edge_list(src)[0]
        assert fp(mut)["canonical_hash"] != fp(rec)["canonical_hash"]
        assert labeled_set(mut) != labeled_set(rec)
        disc = os.path.join(tmp, "drawn_row02_disconnected.yaml")
        open(disc, "w").write(txt.replace("a1-b1[0]; ", "").replace("a3-b3[0]; ", "")
                              .replace("a2-b2[0]", "b1-b3[0]"))
        try:
            ta.load_edge_list(disc)
        except ValueError as exc:
            assert "disconnected" in str(exc) and "REFUSED" in str(exc)
        else:
            raise AssertionError("a disconnected edge list was routed")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_structured_edge_list_forms_agree_with_the_header_form():
    ref = ta.load_edge_list(fixture("row07/drawn_graph.yaml"))[0]
    tmp = _tmpdir()
    try:
        spec = {"name": "drawn_row07", "loop_momenta": ["k1", "k2"],
                "edges": [["a1", "a2", 0], ["a2", "a3", 0], ["b1", "b2", 0],
                          ["b2", "b3", 0], ["a1", "b1", 0], ["a3", "b3", 0],
                          {"u": "a2", "v": "b2", "mass": "msq"}],
                "legs": [["p2", "a1", "0"], ["p1", "b1", "0"],
                         {"name": "p3", "vertex": "a3", "virtuality": "0"},
                         ["p4", "b3", "0"]],
                "cyclic_leg_order": ["p1", "p2", "p3", "p4"]}
        js = ta.load_family(_write_json(spec, os.path.join(tmp, "row07_drawn.json")))[0]
        assert js.source_kind == "edge_list"
        assert js.cyclic_leg_order == ["p1", "p2", "p3", "p4"]
        assert labeled_set(js) == labeled_set(ref)
        yml = os.path.join(tmp, "row07_drawn.yaml")
        open(yml, "w").write(
            "drawn_graph:\n  name: drawn_row07\n  loop_momenta: [k1, k2]\n  edges:\n"
            "    - [a1, a2, 0]\n    - [a2, a3, 0]\n    - [b1, b2, 0]\n    - [b2, b3, 0]\n"
            "    - [a1, b1, 0]\n    - [a3, b3, 0]\n    - [a2, b2, msq]\n  legs:\n"
            "    - [p2, a1, 0]\n    - [p1, b1, 0]\n    - [p3, a3, 0]\n    - [p4, b3, 0]\n"
            "  cyclic_leg_order: [p1, p2, p3, p4]\n")
        yf = ta.load_family(yml)[0]
        assert yf.source_kind == "edge_list"
        assert labeled_set(yf) == labeled_set(ref)
        assert fp(yf)["canonical_hash"] == fp(ref)["canonical_hash"]
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# ---------------------------------------------------------------------------
#  dispatch
# ---------------------------------------------------------------------------
def test_load_family_dispatch_by_form():
    assert ta.load_family(fixture("row21/family_of_record.json"))[0].source_kind == "amflow_json"
    assert ta.load_family(fixture("row05/psd_build_pentagon.py"))[0].source_kind == "pysecdec_graph"
    assert ta.load_family(fixture("row30/lbl3m_adj/graphspec_adj_parked_copy.json"))[0].source_kind == "edge_list"
    assert ta.load_family(fixture("row01/drawn_graph.yaml"))[0].source_kind == "kira_yaml"
    jl = ta.load_family(fixture("row07/family_of_record.jl"))[0]
    assert jl.source_kind == "amflow_jl" and jl.kinematics_source == "jl:exts+ext_subs"
    kin = ta.load_family(fixture("row21/family_of_record.json"), kinematics="k.yaml")[0]
    assert kin.kinematics_path == "k.yaml"


# ---------------------------------------------------------------------------
#  battery runner (used by the tool's --self-test leg [8] and by __main__)
# ---------------------------------------------------------------------------
BATTERY = [
    test_row21_json_reads_four_legs_and_import_route_fingerprint,
    test_row21_control_edited_propagator_changes_hash_and_labeled_set,
    test_row21_control_deleted_conservation_is_refused_by_name,
    test_row28_json_import_route_fingerprint,
    test_row29_json_import_route_fingerprint_and_nu,
    test_row33_json_import_route_fingerprint_and_nu,
    test_rows08_to_12_sunrise_json_reproduce_transcription,
    test_row13_kite_json_reproduces_transcription,
    test_row14_kite112_json_reproduces_transcription_and_record,
    test_rows08_13_14_planted_mass_controls_flip_with_labels_only,
    test_rows15_16_17_icc_json_reproduce_transcriptions,
    test_row16_control_permuted_masses_flip_labeled_set_not_hash,
    test_row18_single_leg_without_conservation_key_loads,
    test_row23_threshold_kite_json_reproduces_transcription,
    test_row24_cured_json_reads_four_listed_legs_with_p5_implicit,
    test_row06_retired_json_reads_three_declared_legs,
    test_bare_family_block_reads_as_a_propagator_list,
    test_json_integral_index_override,
    test_row05_pysecdec_graph_routes_to_the_record_fingerprint,
    test_row05_control_moved_internal_line_changes_hash,
    test_rows01_02_07_drawn_double_boxes_route_as_the_census_did,
    test_row04_drawn_triple_box_routes_as_the_census_did,
    test_row30_graphspec_json_routes,
    test_drawn_controls_moved_rung_and_disconnected_list,
    test_structured_edge_list_forms_agree_with_the_header_form,
    test_load_family_dispatch_by_form,
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
    print(f"reader battery: {len(out)} cases, {n_pass} PASS, {n_fail} FAIL, {n_skip} SKIP")
    sys.exit(1 if n_fail else 0)
