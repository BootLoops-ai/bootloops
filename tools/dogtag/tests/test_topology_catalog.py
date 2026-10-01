#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision.
"""
Catalog battery for topology_audit.py: the graph-identity census rows whose
hand-checked values exercise the catalog are reproduced BY OBJECT from the
vendored record files under tests/fixtures/<row>/ (sha256-pinned here and in
each PROVENANCE.json, which also pins the census record's own sha256 and names
any comment or unread-config line rewritten in the vendored copy; any drift is
refused).

Row 2 (massless planar double box, family in the t = (p1+p3)^2 convention)
  - the catalog's crossed-box entry is the (2,1,1) theta graph: cut signature
    {s: 2, t: 3, u: 3}, non-planar with the legs on the outer face in every
    leg order, not isomorphic as a bare graph to the planar double box, no
    planar entry reachable under any leg map (catalog_check.json values)
  - the RETIRED former entry is the planar double box with legs 3 and 4
    interchanged (the named negative: it fails the crossed-box criteria and
    equals planar_smirnov_dbox under p3 <-> p4; bare-graph leg map onto the
    drawn planar box p4->p1, p3->p2, p1->p3, p2->p4)
  - the family of record (hash a542a386080e, cuts {s: 2, t: 4, u: 3}) no
    longer matches a crossed entry
Row 4 (three-loop three-point off-shell ladder)
  - the family of record matches ud_ladder_3pt_offshell_3l under the identity
    (hash fd256fcd39ad, cuts {p1p2: 2, p1p3: 2, p2p3: 2}, V,E 10,12) and not
    the four-point triple box; its leg classes {p1: offshell, p2: 0, p3:
    offshell} are the entry's "legs"
  - the drawn four-point on-shell triple box is not the ud entry; it matches
    planar_triple_box_3l under the affine class (a second routing of the same
    graph, described by the isomorphism unit and certified by the realized-
    graph isomorphism) and the reading-(b) mass-blind match is the labels
    unit's (both former xfails are plain cases)

Every entry of CATALOG_SOURCES declares "legs", "planar" and
"planar_in_some_leg_order"; the last is checked by enumerating the leg orders.
"""
import hashlib
import itertools
import json
import os
import re
import subprocess
import sys

try:
    import networkx as nx
except ImportError:        # the declared dependency; the cases that need it SKIP by name
    nx = None
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
    "row02/genuine_crossed_box_theta211.yaml":
        "05d1300140c4b9cad4929783b25df5225b1e35623491a0b96903097d5c18690e",
    "row02/catalog_check.json":
        "892422b5b8ba0fac09fe148564989962746b13cfc6237fc3cbd7516c6f58ca5f",
    "row02/family_of_record.jl":
        "ddc75ae24954b08103d01e996f614cae6a6b74a8bbc485c80711cc7e17434917",
    "row04/family_of_record.yaml":
        "c61b80c50011b900a1144d396d50d9c1b41e6c77dc40566da6201abefa0d633c",
    "row04/family_of_record_reading_b_4pt.yaml":
        "e9c7e5d7aecbaba1617d0722563448f34612a5dc1b5de38d6823f42f25d5c9a3",
    "row04/drawn_graph.yaml":
        "87bcbde125b06c812b48408a7b0cc2b5a54ca4700aa1730b8341f337d096915b",
    "row04/iso.json":
        "63f33702433369410533e2ae09a47071f332dcd49a3f7761f18bfc0597da625a",
}

# values quoted from the census record (rows/row02/catalog_check.json,
# rows/row02/ROW.md, rows/row04/ROW.md, rows/row04/iso.json)
THETA_CUT = {"s": 2, "t": 3, "u": 3}
THETA_HASH = "7790fb92864c"
RETIRED_CUT = {"s": 2, "t": 4, "u": 3}
RETIRED_HASH = "a542a386080e"
RETIRED_BARE_LEG_MAP = {"p4": "p1", "p3": "p2", "p1": "p3", "p2": "p4"}
ROW02_HASH = "a542a386080e"
ROW02_CUT = {"s": 2, "t": 4, "u": 3}
ROW04_HASH = "fd256fcd39ad"
ROW04_CUT = {"p1p2": 2, "p1p3": 2, "p2p3": 2}
ROW04_LEGS = {"p1": "offshell", "p2": "0", "p3": "offshell"}
ROW04_VE_FAMILY = (10, 12)
ROW04_VE_DRAWN = (12, 14)

E4 = ["p1", "p2", "p3", "p4"]


def fixture(rel):
    """Path of a pinned fixture; refuses on drift, skips by name when absent."""
    path = os.path.join(FIXTURES, rel)
    if not os.path.exists(path):
        pytest.skip(f"fixture tests/fixtures/{rel} absent")
    got = hashlib.sha256(open(path, "rb").read()).hexdigest()
    if got != PINS[rel]:
        pytest.fail(f"fixture tests/fixtures/{rel} drifted: sha256 {got[:16]} != "
                    f"pinned {PINS[rel][:16]} (refused)")
    return path


def fam_from(name, spec):
    return T.Family(name, spec["loops"], spec["exts"], spec["ext_subs"],
                    spec["propagators"], "memory:" + name)


def bare_graph(fam):
    G, ok = T.build_graph(fam)
    assert ok and G is not None, f"vertex realization failed for {fam.name}"
    base = nx.Graph()
    base.add_edges_from((u, v) for u, v, k in G.edges(keys=True))
    return base


def planar_leg_orders(fam):
    """(abstract_planar, n_planar_orders, n_orders): close the legs in every
    cyclic order on the realized graph."""
    base = bare_graph(fam)
    ext = [n for n in base.nodes() if str(n).startswith("ext_")]
    n_ok = n_all = 0
    for perm in itertools.permutations(ext):
        H = base.copy()
        for a in range(len(perm)):
            H.add_edge(perm[a], perm[(a + 1) % len(perm)])
        n_all += 1
        n_ok += 1 if nx.check_planarity(H)[0] else 0
    return bool(nx.check_planarity(base)[0]), n_ok, n_all


def leg_subst(props, legmap, loops, exts):
    syms = {s: sp.Symbol(s) for s in list(loops) + list(exts)}
    res = []
    for e, m in props:
        lin = (T._extract_linear_from_square(e, syms)
               if re.search(r'\^\s*2|\*\*\s*2', e) else sp.sympify(e, locals=syms))
        lin2 = lin.subs({syms[a]: syms[b] for a, b in legmap.items()}, simultaneous=True)
        res.append((f"({sp.sstr(lin2)})^2", m))
    return res


def complete_leg_search(specA, specB):
    """All leg permutations applied to A's strings (momentum conservation
    re-imposed by the Family's ext_subs), then the tool's signed loop
    relabeling search with legs fixed.  [(legmap, relabeling)]."""
    vb, _, _ = T.momentum_vectors(fam_from("b", specB))
    found = []
    exts = list(specA["exts"])
    for perm in itertools.permutations(exts):
        lm = dict(zip(exts, perm))
        fa = T.Family("a", specA["loops"], exts, specA["ext_subs"],
                      leg_subst(specA["propagators"], lm, specA["loops"], exts), "memory:a")
        va, _, _ = T.momentum_vectors(fa)
        rel = T.isomorphism_to(va, vb, len(specA["loops"]), len(exts))
        if rel is not None:
            found.append((lm, rel))
    return found


def crossed_box_criteria(spec):
    fam = fam_from("candidate", spec)
    fp = T.audit(fam)
    ab, n_ok, n_all = planar_leg_orders(fam)
    reach = {n: complete_leg_search(spec, s) for n, s in T.CATALOG_SOURCES.items()
             if s.get("planar") is True and len(s["loops"]) == 2 and len(s["exts"]) == 4}
    reach = {n: r for n, r in reach.items() if r}
    ok = (fp["cut_signature"] == THETA_CUT and fp["planar"] is False and ab
          and n_all > 0 and n_ok == 0 and not reach)
    return ok, dict(cut=fp["cut_signature"], planar=fp["planar"], hash=fp["canonical_hash"],
                    abstract=ab, n_planar_orders=n_ok, n_orders=n_all, reachable=reach)


def matches(fp):
    return [(m["family"], m["relabeling"]) for m in fp["catalog_matches"]]


def leg_virtualities(path):
    m = re.search(r'leg_virtualities:\s*\{([^}]*)\}', open(path).read())
    return dict(re.findall(r'(p\d+):\s*"([^"]*)"', m.group(1))) if m else None


# --------------------------------------------------------------------- pins --
def test_fixture_provenance_pins():
    """Every file this battery pins: sha256 in code == PROVENANCE.json == bytes
    on disk.  PROVENANCE.json is the package's one canonical shape ("files" is
    a list of {"path", "source", "sha256", "size"}; a row directory may also
    carry files vendored by the other batteries, which pin them themselves)."""
    seen = 0
    provs = {}
    for rel in sorted(PINS):
        row, within = rel.split("/", 1)
        if row not in provs:
            provs[row] = json.load(open(os.path.join(FIXTURES, row, "PROVENANCE.json")))
        entries = [e for e in provs[row]["files"] if e.get("path") == within]
        assert len(entries) == 1, f"{rel}: not recorded once in {row}/PROVENANCE.json (refused)"
        e = entries[0]
        path = os.path.join(FIXTURES, rel)
        got = hashlib.sha256(open(path, "rb").read()).hexdigest()
        assert got == e["sha256"] == PINS[rel], f"{rel}: pin mismatch (refused)"
        assert os.path.getsize(path) == e["size"], f"{rel}: size mismatch"
        seen += 1
    assert seen == len(PINS)


# ------------------------------------------------------------- row 2: theta --
@nx_required
def test_theta_entry_is_the_crossed_box():
    spec = T.CATALOG_SOURCES["crossed_tausk_dbox"]
    ok, d = crossed_box_criteria(spec)
    assert d["cut"] == THETA_CUT
    assert d["planar"] is False and d["abstract"] is True
    assert d["n_orders"] == 24 and d["n_planar_orders"] == 0, d
    assert d["reachable"] == {}, d
    assert d["hash"] == THETA_HASH
    assert ok
    assert spec["planar"] is False and spec["planar_in_some_leg_order"] is False
    assert spec["legs"] == {"p1": "0", "p2": "0", "p3": "0", "p4": "0"}


@nx_required
def test_theta_fixture_reproduces_catalog_check():
    ref = json.load(open(fixture("row02/catalog_check.json")))["genuine_crossed_box_theta211"]
    fam = T.load_family(fixture("row02/genuine_crossed_box_theta211.yaml"))[0]
    fp = T.audit(fam, try_leg_perms=True)
    assert fp["cut_signature"] == ref["census_fingerprint"]["cut_signature"] == THETA_CUT
    assert fp["planar"] is ref["census_fingerprint"]["planar_with_legs_in_canonical_order"] is False
    assert fp["canonical_hash"] == THETA_HASH
    assert matches(fp) == [("crossed_tausk_dbox", "identity")]
    ab, n_ok, n_all = planar_leg_orders(fam)
    assert ab is ref["census_fingerprint"]["abstract_planar"] is True
    assert (n_ok > 0) is ref["planar_with_legs_for_some_leg_order"] is False
    # the entry's propagator strings are the census's, line for line
    assert [e for e, _ in fam.propagators] == ref["propagators"]
    assert [re.sub(r'^\((.*)\)\^2$|^(k\d)\^2$', lambda m: m.group(1) or m.group(2), e)
            for e, _ in T.CATALOG_SOURCES["crossed_tausk_dbox"]["propagators"]] == ref["propagators"]
    # bare-graph statements of the record
    gt = bare_graph(fam)
    gr = bare_graph(fam_from("retired", T.CATALOG_RETIRED["planar_smirnov_dbox_legs34_interchanged"]))
    gs = bare_graph(fam_from("planar", T.CATALOG_SOURCES["planar_smirnov_dbox"]))
    assert nx.is_isomorphic(gt, gr) is ref["graph_iso_to_catalog_crossed_tausk_dbox_bare"] is False
    # the record compared against the DRAWN planar box; the catalog planar box
    # is that graph (the retired entry maps onto it bare, see the next test)
    assert nx.is_isomorphic(gt, gs) is ref["graph_iso_to_drawn_planar_box_bare"] is False


@nx_required
def test_retired_entry_is_the_planar_box_relabeled():
    """Named negative: the former 'crossed_tausk_dbox' fails the crossed-box
    criteria and IS planar_smirnov_dbox under p3 <-> p4."""
    ret = T.CATALOG_RETIRED["planar_smirnov_dbox_legs34_interchanged"]
    assert ret["former_name"] == "crossed_tausk_dbox"
    assert "crossed_tausk_dbox" in T.CATALOG_SOURCES
    assert T.CATALOG_SOURCES["crossed_tausk_dbox"]["propagators"] != ret["propagators"]
    assert not set(T.CATALOG_RETIRED) & set(T.CATALOG_SOURCES)
    ok, d = crossed_box_criteria(ret)
    assert ok is False
    assert d["cut"] == RETIRED_CUT and d["hash"] == RETIRED_HASH
    assert d["planar"] is False and d["abstract"] is True
    assert d["n_planar_orders"] == 8 and d["n_orders"] == 24, d
    assert "planar_smirnov_dbox" in d["reachable"]
    legmaps = [lm for lm, _ in d["reachable"]["planar_smirnov_dbox"]]
    assert {"p1": "p1", "p2": "p2", "p3": "p4", "p4": "p3"} in legmaps
    # the record's bare_leg_map is written family->drawn; the search applies
    # the map to the retired strings, so the record's map appears inverted
    inv = {v: k for k, v in RETIRED_BARE_LEG_MAP.items()}
    assert inv in legmaps, legmaps
    smi = T.CATALOG_SOURCES[ret["is"]]
    fsw = T.Family("retired_p3p4", ret["loops"], ret["exts"], ret["ext_subs"],
                   leg_subst(ret["propagators"], ret["is_under_leg_map"], ret["loops"], ret["exts"]),
                   "memory:retired_p3p4")
    assert T.isomorphism_to(T.momentum_vectors(fsw)[0],
                            T.momentum_vectors(fam_from("smi", smi))[0], 2, 4) == "identity"
    assert nx.is_isomorphic(bare_graph(fam_from("r", ret)), bare_graph(fam_from("s", smi)))
    # the retired object is not searched: auditing it matches nothing crossed
    fp = T.audit(fam_from("retired_as_family", ret), try_leg_perms=True)
    assert all("crossed" not in f for f, _ in matches(fp))
    assert all(f in T.CATALOG_SOURCES for f, _ in matches(fp))


def test_tool_leg_perms_find_p3p4_on_the_retired_entry():
    """The tool's own --leg-perms search finds the p3 <-> p4 relabeling once
    the target's momentum conservation is re-imposed after the permutation
    (the leg was an expected-fail until the leg-perms unit re-imposed it)."""
    ret = T.CATALOG_RETIRED["planar_smirnov_dbox_legs34_interchanged"]
    smi = fam_from("s", T.CATALOG_SOURCES["planar_smirnov_dbox"])
    rel = T.isomorphism_to(T.momentum_vectors(fam_from("r", ret))[0],
                           T.momentum_vectors(smi)[0], 2, 4, try_leg_perms=True,
                           dep_rows=T._dependent_leg_rows(smi))
    assert rel is not None and "p3 -> p4" in rel and "p4 -> p3" in rel, rel


@nx_required
def test_row02_family_no_longer_reads_crossed():
    fam = T.load_family(fixture("row02/family_of_record.jl"))[0]
    fp = T.audit(fam, try_leg_perms=True)
    assert fp["canonical_hash"] == ROW02_HASH
    assert fp["cut_signature"] == ROW02_CUT
    assert fp["planar"] is False          # canonical-order closure, as the record says
    assert fp["n_genuine_propagators"] == 7 and fp["n_isp_or_dotproduct"] == 2
    mts = matches(fp)
    assert all("crossed" not in f for f, _ in mts), mts
    assert all(f == "planar_smirnov_dbox" for f, _ in mts), mts
    assert not any("crossed" in w for w in fp["warnings"]), fp["warnings"]


# ------------------------------------------------------------ row 4: ladder --
@nx_required
def test_ud_ladder_entry_fingerprint():
    ud = T.CATALOG_SOURCES["ud_ladder_3pt_offshell_3l"]
    assert ud["legs"] == ROW04_LEGS and ud["planar"] is True
    fam = fam_from("ud", ud)
    fp = T.audit(fam)
    assert fp["canonical_hash"] == ROW04_HASH
    assert fp["cut_signature"] == ROW04_CUT and fp["planar"] is True
    assert fp["n_loops"] == 3 and fp["n_genuine_propagators"] == 9
    pl = T.planarity(fam)
    assert (pl["n_vertices"], pl["n_edges"]) == ROW04_VE_FAMILY
    assert matches(fp) == [("ud_ladder_3pt_offshell_3l", "identity")]


@nx_required
def test_row04_family_matches_ud_not_triple_box():
    path = fixture("row04/family_of_record.yaml")
    fam = T.load_family(path)[0]
    fp = T.audit(fam, try_leg_perms=True)
    assert fp["canonical_hash"] == ROW04_HASH and fp["cut_signature"] == ROW04_CUT
    assert fp["planar"] is True
    assert matches(fp) == [("ud_ladder_3pt_offshell_3l", "identity")]
    assert leg_virtualities(path) == ROW04_LEGS == T.CATALOG_SOURCES["ud_ladder_3pt_offshell_3l"]["legs"]
    tb = fam_from("tb", T.CATALOG_SOURCES["planar_triple_box_3l"])
    assert T.isomorphism_to(T.momentum_vectors(fam)[0], T.momentum_vectors(tb)[0], 3, 4) is None
    # the record: "not isomorphic as graphs (V,E family 10,12 vs drawn 12,14)"
    pl = T.planarity(fam)
    assert (pl["n_vertices"], pl["n_edges"]) == ROW04_VE_FAMILY


def test_row04_drawn_is_not_the_ud_entry():
    """The census FAMILY-MISMATCH direction at the propagator level: the drawn
    four-point on-shell triple box (10 lines, 4 legs) is not the three-point
    off-shell ladder (9 lines, 3 legs), and its legs are all on shell."""
    path = fixture("row04/drawn_graph.yaml")
    fam = T.load_family(path)[0]
    ud = fam_from("ud", T.CATALOG_SOURCES["ud_ladder_3pt_offshell_3l"])
    assert T.isomorphism_to(T.momentum_vectors(fam)[0], T.momentum_vectors(ud)[0], 3, 4) is None
    fp = T.audit(fam, try_leg_perms=True)
    assert all(f != "ud_ladder_3pt_offshell_3l" for f, _ in matches(fp))
    assert fp["n_genuine_propagators"] == 10 and len(fam.exts) == 4
    lv = leg_virtualities(path)
    assert lv == T.CATALOG_SOURCES["planar_triple_box_3l"]["legs"] == {"p1": "0", "p2": "0", "p3": "0", "p4": "0"}
    assert lv != ROW04_LEGS
    ref = json.load(open(fixture("row04/iso.json")))["family_of_record_vs_drawn"]
    assert ref["graph"]["bare"] is False and (ref["graph"]["V_A"], ref["graph"]["E_A"]) == ROW04_VE_FAMILY


@nx_required
def test_row04_drawn_matches_planar_triple_box():
    """The drawn routing and the catalog routing differ by a loop-momentum
    basis change: the isomorphism unit's affine search describes it and the
    realized-graph isomorphism certifies it (the match names its class); the
    realizer returns V,E = 12,14 for this drawing."""
    fam = T.load_family(fixture("row04/drawn_graph.yaml"))[0]
    fp = T.audit(fam, try_leg_perms=True)
    pl = T.planarity(fam)
    assert (pl["n_vertices"], pl["n_edges"]) == ROW04_VE_DRAWN, (pl["n_vertices"], pl["n_edges"])
    assert any(f == "planar_triple_box_3l" for f, _ in matches(fp)), matches(fp)
    m = [x for x in fp["catalog_matches"] if x["family"] == "planar_triple_box_3l"][0]
    assert m["relabeling_class"] == "affine unimodular loop redefinition"
    assert m["realized_graph_iso"]["with_masses_and_legs"] is True
    assert (m["realized_graph_iso"]["V_A"], m["realized_graph_iso"]["E_A"]) == ROW04_VE_DRAWN


def test_row04_reading_b_mass_blind_match_and_leg_label_fail():
    """Reading (b) against the DRAWN four-point triple box (the census object,
    iso.json "extra"): the same ten lines and masses under some relabeling
    ("momentum_isomorphic": true, 2 relabelings, first "identity") but no
    relabeling carries the off-shell-leg labels ("labelled_isomorphic": false)
    -- the labels unit's labeled_isomorphism reports both.  (The catalog
    entry planar_triple_box_3l is routed differently from the drawing; its
    match is the routing-independent isomorphism's, test_row04_drawn_matches_
    planar_triple_box above.)"""
    ref = json.load(open(fixture("row04/iso.json")))["extra"]["family_of_record_reading_b_4pt.yaml"]
    assert ref["momentum_isomorphic"] is True and ref["labelled_isomorphic"] is False
    fb = T.load_family(fixture("row04/family_of_record_reading_b_4pt.yaml"))[0]
    fd = T.load_family(fixture("row04/drawn_graph.yaml"))[0]
    r = T.labeled_isomorphism(fb, fd, find_all=True)
    assert (len(r["mass_match"]) > 0) is ref["momentum_isomorphic"] is True
    assert len(r["mass_match"]) == ref["n_momentum_relabelings_found"] == 2
    assert r["mass_match"][0] == ref["first_momentum_relabeling"] == "identity"
    assert bool(r["labeled_match"]) is ref["labelled_isomorphic"] is False
    assert r["leg_classes_A"] == ["offshell", "0", "offshell", "0"]
    assert r["leg_classes_B"] == ["0", "0", "0", "0"]
    assert r["reason"] == ref["reason"]
    fp = T.audit(fb, try_leg_perms=True)
    assert fp["leg_labels"] == {"p1": "offshell", "p2": "0", "p3": "offshell", "p4": "0"}
    assert fp["mass_multiset"] == ["0"] * 10 and fp["n_massive_lines"] == 0


def test_row04_reading_b_leg_classes_by_object():
    """Provable today: reading (b) carries the census leg classes p1, p3 off
    shell (its record label) while the drawing's are all on shell, so the two
    differ on the legs alone."""
    lv_b = leg_virtualities(fixture("row04/family_of_record_reading_b_4pt.yaml"))
    lv_d = leg_virtualities(fixture("row04/drawn_graph.yaml"))
    assert lv_b == {"p1": "offshell", "p2": "0", "p3": "offshell", "p4": "0"}
    assert lv_d == {"p1": "0", "p2": "0", "p3": "0", "p4": "0"}
    fb = T.load_family(fixture("row04/family_of_record_reading_b_4pt.yaml"))[0]
    fd = T.load_family(fixture("row04/drawn_graph.yaml"))[0]
    assert T.isomorphism_to(T.momentum_vectors(fb)[0], T.momentum_vectors(fd)[0], 3, 4) == "identity"


# ------------------------------------------------------- catalog integrity --
@nx_required
def test_every_entry_declares_legs_and_leg_order_planarity():
    for cname, spec in list(T.CATALOG_SOURCES.items()) + list(T.CATALOG_RETIRED.items()):
        assert list(spec["legs"]) == list(spec["exts"]), cname
        assert all(v in ("0", "offshell") for v in spec["legs"].values()), cname
        assert spec["planar"] in (True, False), cname
        assert spec["planar_in_some_leg_order"] in (True, False), cname
        ab, n_ok, n_all = planar_leg_orders(fam_from(cname, spec))
        assert n_all > 0, cname
        assert (n_ok > 0) is spec["planar_in_some_leg_order"], (cname, n_ok, n_all)
        fp = T.audit(fam_from(cname, spec))
        assert fp["planar"] is spec["planar"], cname


def test_catalog_has_nine_entries_and_one_retired():
    assert len(T.CATALOG_SOURCES) == 9
    assert list(T.CATALOG_RETIRED) == ["planar_smirnov_dbox_legs34_interchanged"]
    names = [f.name for f, _ in T._catalog_families()]
    assert names == list(T.CATALOG_SOURCES)


def test_planted_mutation_of_the_ud_ladder_does_not_match():
    ud = dict(T.CATALOG_SOURCES["ud_ladder_3pt_offshell_3l"])
    props = list(ud["propagators"])
    assert props[3] == ("(k3 - p1)^2", 0)
    props[3] = ("(k3 - p1 - p2)^2", 0)      # the far rung re-attached
    mut = dict(ud, propagators=props)
    fp = T.audit(fam_from("ud_mutated", mut), try_leg_perms=True)
    assert all(f != "ud_ladder_3pt_offshell_3l" for f, _ in matches(fp)), matches(fp)
    assert fp["canonical_hash"] != ROW04_HASH


# ------------------------------------------------------------- self-test --
@nx_required
def test_selftest_catalog_leg_runs_clean():
    r = T._selftest_catalog(verbose=False)
    assert r["skipped"] == 0, r["skipped_names"]
    assert r["ok"] is True, {k: v for k, v in r["results"].items() if v.get("ok") is False}
    assert r["executed"] >= 10
    # registered as (leg id, description, callable) in the package's one
    # EXTRA_SELF_TEST_LEGS list (the _selftest_catalog unit)
    assert any(fn is T._selftest_catalog for _lid, _title, fn in T.EXTRA_SELF_TEST_LEGS)


def test_fixture_env_override_names(monkeypatch, tmp_path):
    """The fixture overrides are read under the package's name first
    (DOGTAG_<X>) and under its former name second (TOPOLOGY_AUDIT_<X>): the
    rename of topology-audit to Dogtag keeps the old spelling working.  Either
    spelling passed to the resolver names the same pair; a set primary wins
    over a set legacy; a set override that does not exist is a named SKIP,
    never a silent fall-back to the vendored copy; with neither set the
    vendored copy is used and an absent one names both spellings."""
    f_new = tmp_path / "new.jl"
    f_old = tmp_path / "old.jl"
    f_new.write_text("x")
    f_old.write_text("y")
    for var in ("DOGTAG_C3_FAMILY", "TOPOLOGY_AUDIT_C3_FAMILY"):
        monkeypatch.delenv(var, raising=False)
    assert T._env_names("DOGTAG_C3_FAMILY") == ("DOGTAG_C3_FAMILY", "TOPOLOGY_AUDIT_C3_FAMILY")
    assert T._env_names("TOPOLOGY_AUDIT_C3_FAMILY") == T._env_names("DOGTAG_C3_FAMILY")
    assert T._env_names("C3_FAMILY") == T._env_names("DOGTAG_C3_FAMILY")
    # neither set: the vendored row-07 copy
    path, why = T._fixture_or_skip("DOGTAG_C3_FAMILY", os.path.join("row07", "family_of_record.jl"), "C3")
    assert why is None and path.endswith(os.path.join("row07", "family_of_record.jl"))
    path, why = T._fixture_or_skip("DOGTAG_C3_FAMILY", os.path.join("nosuchrow", "x.jl"), "C3")
    assert path is None and "set DOGTAG_C3_FAMILY (or TOPOLOGY_AUDIT_C3_FAMILY) to run" in why
    # legacy alone
    monkeypatch.setenv("TOPOLOGY_AUDIT_C3_FAMILY", str(f_old))
    assert T._fixture_or_skip("DOGTAG_C3_FAMILY", None, "C3") == (str(f_old), None)
    assert T._fixture_or_skip("TOPOLOGY_AUDIT_C3_FAMILY", None, "C3") == (str(f_old), None)
    # primary wins over legacy
    monkeypatch.setenv("DOGTAG_C3_FAMILY", str(f_new))
    assert T._fixture_or_skip("DOGTAG_C3_FAMILY", None, "C3") == (str(f_new), None)
    # a set override that does not exist is named, not replaced by the vendored copy
    monkeypatch.setenv("DOGTAG_C3_FAMILY", str(tmp_path / "absent.jl"))
    path, why = T._fixture_or_skip("DOGTAG_C3_FAMILY", os.path.join("row07", "family_of_record.jl"), "C3")
    assert path is None and why.startswith("C3: DOGTAG_C3_FAMILY=") and why.endswith("does not exist")


@nx_required
def test_cli_self_test_green():
    """The CLI self-test runs every registered leg with no FAIL, under the
    contract: "SELF-TEST: ALL PASS" is printed only when 0 legs were skipped
    AND 0 expected-fail sub-legs are open; the exit code is 2 whenever any leg
    was SKIPPED (a skip is an unknown), 1 on any FAIL, else 0 (declared
    expected-fail sub-legs keep rc 0 but hold back ALL PASS).  On a scratch
    install leg [1] (the NP-dbox family, an internal fixture behind
    DOGTAG_NPDBOX_FAMILY, legacy name TOPOLOGY_AUDIT_NPDBOX_FAMILY) is the one
    named SKIP, so the run exits 2
    with "NOT ALL PASS (1 skipped ...)"; with that fixture set nothing is
    skipped and the run exits 0.  An open expected-fail sub-leg would be named
    on the summary line and hold back ALL PASS; none is open ([8]'s drawn
    triple-box matches were lifted by the affine / realized-graph isomorphism
    unit [isomorphism], [3b]'s by the leg-perms unit), so this test reads the skipped
    and expected-fail counts off the run itself and asserts the
    "SELF-TEST: ALL PASS (rc 0)" line when both are 0, the NOT-ALL-PASS
    summary wording otherwise (a scratch install)."""
    p = subprocess.run([sys.executable, os.path.join(PKG, "topology_audit.py"), "--self-test"],
                       capture_output=True, text=True, cwd=HERE)
    assert p.returncode in (0, 2), p.stdout[-2000:] + p.stderr[-2000:]
    assert "SELF-TEST: SOME FAILED" not in p.stdout
    assert "\nFAILED  [" not in p.stdout
    m = re.search(r"^SELF-TEST: (\d+) executed / (\d+) skipped \((\d+) PASS, (\d+) FAIL\) of (\d+) legs$",
                  p.stdout, flags=re.M)
    assert m, p.stdout[-2000:]
    n_skipped, n_fail = int(m.group(2)), int(m.group(4))
    assert n_fail == 0
    ms = re.search(r"^SELF-TEST: registered sub-legs (\d+) executed / (\d+) skipped / (\d+) expected-fail",
                   p.stdout, flags=re.M)
    assert ms, p.stdout[-2000:]
    n_xfail = int(ms.group(3))
    skips = re.findall(r"^SKIPPED \[(\S+)\]", p.stdout, flags=re.M)
    assert len(skips) == n_skipped
    npd = (os.environ.get("DOGTAG_NPDBOX_FAMILY", "")
           or os.environ.get("TOPOLOGY_AUDIT_NPDBOX_FAMILY", ""))
    vendored_npd = os.path.join(PKG, "tests", "fixtures", "npdbox", "family.jl")
    if (npd and os.path.exists(npd)) or os.path.exists(vendored_npd):
        assert skips == [], skips        # leg [1] runs on the vendored NP-dbox fixture
    else:
        assert skips == ["1"], skips
    # the contract, read off the run itself
    assert p.returncode == (2 if n_skipped else 0), (p.returncode, n_skipped)
    if n_skipped == 0 and n_xfail == 0:
        assert "SELF-TEST: ALL PASS (rc 0)" in p.stdout
    else:
        assert "SELF-TEST: ALL PASS" not in p.stdout
        named = ", named above" if n_skipped else ""
        assert (f"SELF-TEST: PASS on the executed legs; NOT ALL PASS ({n_skipped} skipped{named}; "
                f"{n_xfail} open EXPECTED-FAIL sub-legs by name) "
                f"(rc {p.returncode})") in p.stdout
    assert "crossed_tausk_dbox: props=7/7  planar=False (expect False)" in p.stdout
    assert "ud_ladder_3pt_offshell_3l: props=9/9  planar=True (expect True)" in p.stdout
    assert "[catalog] " in p.stdout and "[realizer] " in p.stdout


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
