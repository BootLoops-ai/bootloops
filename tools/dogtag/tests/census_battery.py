#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision.
"""
The COMPARE BATTERY of topology_audit.py's compare mode (compare_drawn).

The identity-census corpus is 30 published Feynman-integral families, each
paired with the graph the paper draws for it (several rows carry more than one
family: a corrected family beside the one it replaced, or two variants).  The
census assigned every family-level pair one verdict string; this battery
reproduces every one of them with the tool and asserts STRING equality.

fixtures/CENSUS_BATTERY.json (the manifest) drives the runs: per family-level
entry it names the sha-pinned fixtures that ARE the corpus's family file,
kinematics and drawn graph (each pinned again in the row's PROVENANCE.json,
which also records the digest of the original corpus file) and the expected
verdict string; the dual cases (one family against a second reading of the
drawing) carry their own expected strings; the row-level "see families: X / Y
[/ Z]" strings are rollups the battery reproduces by joining its own family
verdicts in the manifest's order.  A family file that lists several families
names the one the census judged ("family_name"): the family side is loaded by
that name and the side read is asserted to be it.

The manifest also carries the census's verdict-count summary and its counting
rule; both count summaries (per row: a row's family-level entries replace its
row-level entry; per entry: every entry, the rollups in their own bucket) are
re-derived here from the entries and rollups.  A corpus row the census did not
compare (rows_not_compared: no drawn figure in the paper) has no entry, no
rollup and no fixture in the manifest, and the battery asserts that absence.

No verdict string is typed here: the census's three strings come from the
manifest's summary_verdict_counts (its keys), the fourth (NOT-CHECKABLE, which
the census never emitted) from the tool's own vocabulary VERDICT_STRINGS, said
so by name.

A verdict the tool cannot reproduce is a FAIL of the named entry, by row and
class, never narrowed.

Runs under pytest through test_census_battery.py and as the tool's self-test
leg [compare] (run_battery()).  A missing fixture is SKIP by name; a drifted
fixture is a refusal (FixtureDrift), never a pass.  After an intentional
fixture edit, tests/repin_fixtures.py rewrites the pins.
"""
import functools
import hashlib
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
TOOL_DIR = os.path.dirname(HERE)
sys.path.insert(0, TOOL_DIR)
import topology_audit as ta  # noqa: E402

FIX = os.path.join(HERE, "fixtures")
MANIFEST = os.path.join(FIX, "CENSUS_BATTERY.json")
CENSUS_DIRNAME = "identity_census"
# the census's photons-per-vertex reading of the two row-27 families:
ROW27_PHOTONS_PER_VERTEX = {"retired_threepoint_LBL3E": [2, 1, 1, 0],
                            "cured_drawn_box_LBL3Ebox": [1, 1, 1, 1]}


class FixtureDrift(Exception):
    pass


class FixtureMissing(Exception):
    pass


def nx_required(fn):
    """A case that needs networkx (a declared dependency): pytest SKIPs it by
    name, run_battery reports it SKIP through the missing-fixture path."""
    @functools.wraps(fn)
    def gated(*a, **k):
        why = ta._require_nx()
        if why:
            raise FixtureMissing(why)
        return fn(*a, **k)
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


def load_manifest():
    if not os.path.exists(MANIFEST):
        raise FixtureMissing("fixture absent: tests/fixtures/CENSUS_BATTERY.json")
    return json.load(open(MANIFEST))


def fixture(info):
    """Path of a vendored fixture named by a manifest file record after
    checking its sha256 against the row's PROVENANCE.json and the manifest's
    own pin; raises FixtureMissing / FixtureDrift."""
    rel = info["fixture"]
    path = os.path.join(FIX, rel)
    if not os.path.exists(path):
        raise FixtureMissing(f"fixture absent: tests/fixtures/{rel}")
    sha = _sha256(path)
    if sha[:16] != info.get("vendored_sha16"):
        raise FixtureDrift(f"fixture drift: {rel} sha16 {sha[:16]} != manifest pin {info.get('vendored_sha16')}")
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
        raise FixtureDrift(f"fixture drift: {rel} not recorded with sha256 {sha} in {row}/PROVENANCE.json")
    if entry.get("record_sha256", "")[:16] != info.get("record_sha16"):
        raise FixtureDrift(f"{rel}: PROVENANCE record_sha256 != the manifest's record_sha16")
    if not str(entry.get("source", "")).startswith(CENSUS_DIRNAME + "/rows/"):
        raise FixtureDrift(f"{rel}: PROVENANCE source is not census-relative")
    return path


def verdict_strings(man):
    """The verdict vocabulary BY OBJECT: the census's three strings are the
    keys of its summary.verdict_counts (copied into the manifest); the fourth
    is the tool's own VERDICT_STRINGS entry the census never emitted."""
    counts = man["summary_verdict_counts"]
    out = {}
    for s in ta.VERDICT_STRINGS:
        if s in counts:
            out[s] = s
    extra = [s for s in ta.VERDICT_STRINGS if s not in counts]
    return {"census": out, "tool_only": extra}


def _kin(info):
    return fixture(info) if info else None


def run_pair(fam_info, kin_info, drawn_info, drawn_kin_info, integral_index=None, family_name=None):
    """compare_drawn on a manifest pair; returns (verdict, result, wall_s).
    With `family_name` (a file listing several families) the family side is
    loaded by that name through load_compare_side (a name the file does not
    carry raises by name) and the side the compare read is asserted to be it."""
    t0 = time.time()
    if family_name:
        fam_path = fixture(fam_info)
        A, _ = ta.load_compare_side(fam_path, kinematics=_kin(kin_info), name=family_name,
                                    integral_index=integral_index)
        res = ta.compare_drawn(A, fixture(drawn_info), kinB=_kin(drawn_kin_info),
                               integral_indexA=integral_index)
        side = (res.get("sides") or {}).get("family") or {}
        if side.get("name") != family_name:
            raise FixtureDrift(f"{os.path.basename(fam_path)}: the compare read family {side.get('name')!r}, "
                               f"not the named {family_name!r}")
    else:
        res = ta.compare_drawn(fixture(fam_info), fixture(drawn_info), kinA=_kin(kin_info),
                               kinB=_kin(drawn_kin_info), integral_indexA=integral_index)
    return res["verdict"], res, round(time.time() - t0, 3)


def expected_for(entry_key, man_expected, man):
    """The expected string of an entry or rollup: the manifest's copy of the
    census's verdict (asserted to be one of the census's strings, or a rollup
    joined from them)."""
    census = set(man["summary_verdict_counts"])
    parts = man_expected[len("see families: "):].split(" / ") if man_expected.startswith("see families: ") \
        else [man_expected]
    if not all(p in census for p in parts):
        raise FixtureDrift(f"{entry_key}: manifest expected {man_expected!r} is not made of the census's strings")
    return man_expected


def check_entry(entry, man):
    """One family-level entry: the pair, then the extras and the alsos; returns
    a dict with got / expected per pair and the verdict lines."""
    key = entry["key"]
    exp = expected_for(key, entry["pair"]["expected"], man)
    p = entry["pair"]
    ii = p.get("integral_index")
    fn = p.get("family_name")
    got, res, wall = run_pair(p["family"], p["kinematics"], p["drawn"], p["drawn_kinematics"], ii, fn)
    out = {"key": key, "expected": exp, "got": got, "ok": got == exp, "wall_s": wall,
           "line": res["line"], "result": res, "extras": [], "also": [], "family_name": fn}
    for x in entry.get("extras", []):
        g2, r2, w2 = run_pair(p["family"], p["kinematics"], x["drawn"], x.get("drawn_kinematics"), ii, fn)
        out["extras"].append({"drawn": x["drawn"]["fixture"], "expected": x["expected"],
                              "reading": x["reading"], "got": g2,
                              "ok": g2 == x["expected"], "wall_s": w2, "line": r2["line"], "result": r2})
    for a in entry.get("also", []):
        g3, r3, w3 = run_pair(a["family"], a.get("kinematics"), a["drawn"], None)
        out["also"].append({"family": a["family"]["fixture"], "drawn": a["drawn"]["fixture"],
                            "expected": a["expected"], "got": g3, "ok": g3 == a["expected"],
                            "wall_s": w3, "why": a["why"], "line": r3["line"], "result": r3})
    return out


# ---------------------------------------------------------------------------
#  the cases
# ---------------------------------------------------------------------------
_RESULTS = {}


def _entry_case(entry):
    man = load_manifest()

    @nx_required
    def case():
        r = check_entry(entry, man)
        _RESULTS[entry["key"]] = r
        assert r["ok"], (f"{entry['key']}: tool {r['got']!r} != census {r['expected']!r} "
                         f"({r['line']})")
        for x in r["extras"]:
            assert x["ok"], (f"{entry['key']} vs {x['drawn']}: tool {x['got']!r} != census "
                             f"{x['expected']!r} ({x['reading']}; {x['line']})")
        for a in r["also"]:
            assert a["ok"], (f"{entry['key']} also {a['family']} vs {a['drawn']}: tool {a['got']!r} != "
                             f"census {a['expected']!r} ({a['line']})")
        return r
    case.__name__ = "census_" + entry["key"].replace("/", "__")
    case.__doc__ = f"{entry['key']}: the census's {entry['pair']['expected']!r}"
    return case


def _rollup_case(roll):
    man = load_manifest()

    @nx_required
    def case():
        exp = expected_for(roll["key"], roll["verdict"], man)
        parts = []
        for fk in roll["families_in_order"]:
            r = _RESULTS.get(fk)
            if r is None:
                ent = [e for e in man["entries"] if e["key"] == fk]
                assert ent, f"{roll['key']}: family entry {fk} absent from the manifest"
                r = check_entry(ent[0], man)
                _RESULTS[fk] = r
            parts.append(r["got"])
        got = "see families: " + " / ".join(parts)
        assert got == exp, f"{roll['key']}: rollup {got!r} != census {exp!r}"
        return {"key": roll["key"], "expected": exp, "got": got, "ok": got == exp}
    case.__name__ = "census_rollup_" + roll["key"]
    case.__doc__ = f"{roll['key']}: the census's rollup {roll['verdict']!r} joined from the family verdicts"
    return case


def _entry(man, key):
    for e in man["entries"]:
        if e["key"] == key:
            return e
    raise FixtureMissing(f"manifest entry absent: {key}")


def _side(info, kin=None, integral_index=None):
    fam, notes = ta.load_compare_side(fixture(info), kinematics=_kin(kin),
                                      integral_index=integral_index)
    return fam, notes


def _family_side(e, name=None):
    """The family side of a manifest entry as the census ran it (its
    kinematics, for an AMFlow-port JSON the census's integral index, for a
    multi-family file the census's family by name)."""
    fam, notes = ta.load_compare_side(fixture(e["pair"]["family"]), kinematics=_kin(e["pair"]["kinematics"]),
                                      name=name or e["pair"].get("family_name"),
                                      integral_index=e["pair"].get("integral_index"))
    return fam, notes


def _copy(fam, **over):
    kw = dict(name=fam.name, loops=fam.loops, exts=fam.exts, ext_subs=fam.ext_subs,
              propagators=fam.propagators, source=fam.source, physical=fam.physical, nu=fam.nu,
              leg_virt=fam.leg_virt, mass_values=fam.mass_values,
              kinematics_source=fam.kinematics_source, cyclic_leg_order=fam.cyclic_leg_order,
              source_kind=fam.source_kind, record_notes=fam.record_notes)
    kw.update(over)
    return ta.Family(**kw)


@nx_required
def control_mass_moved_flips_to_family_mismatch_at_the_mass_level():
    """Row 28 (box with the kite self-energy): the drawn graph with one massless
    kite-rim line made massive -- same bare graph, a different mass pattern --
    must read the census's FAMILY-MISMATCH string with the first failing level
    'masses', never the census's IDENTITY-PASS string."""
    man = load_manifest()
    vs = verdict_strings(man)["census"]
    e = _entry(man, "row28")
    A, _ = _family_side(e)
    B, _ = _side(e["pair"]["drawn"])
    base = ta.compare_drawn(A, B)
    assert base["verdict"] == e["pair"]["expected"]
    props = list(B.propagators)
    i0 = [i for i, (_x, m) in enumerate(props) if str(m) in ("0", "0.0")][0]
    props[i0] = (props[i0][0], "msq")
    Bm = _copy(B, propagators=props, name=B.name + "_massmoved")
    r = ta.compare_drawn(A, Bm)
    assert r["verdict"] != e["pair"]["expected"]
    assert r["verdict"] == vs["FAMILY-MISMATCH"], r["line"]
    assert r["first_failing_level"] == "masses", r["line"]
    assert r["legs_per_vertex"]["A"] and r["legs_per_vertex"]["B"], r["line"]


@nx_required
def control_offshell_leg_moved_flips_to_family_mismatch_at_the_leg_level():
    """Row 35 (W pair T4): the family's off-shell W legs p3, p4 read on shell and
    a massless quark leg p1 read off shell instead -- same graph, same masses,
    other leg classes -- must read the census's FAMILY-MISMATCH string at the
    level 'legs'."""
    man = load_manifest()
    vs = verdict_strings(man)["census"]
    e = _entry(man, "row35")
    A, _ = _family_side(e)
    B, _ = _side(e["pair"]["drawn"])
    assert ta.compare_drawn(A, B)["verdict"] == e["pair"]["expected"]
    virt = dict(A.leg_virt)
    off = [k for k, v in virt.items() if ta._leg_class(v) == "offshell"]
    on = [k for k, v in virt.items() if ta._leg_class(v) == "0"]
    assert len(off) == 2 and len(on) == 2, virt
    moved = {k: "0" for k in virt}
    moved[on[0]] = virt[off[0]]
    Am = _copy(A, leg_virt=moved, name=A.name + "_legmoved")
    r = ta.compare_drawn(Am, B)
    assert r["verdict"] != e["pair"]["expected"]
    assert r["verdict"] == vs["FAMILY-MISMATCH"], r["line"]
    assert r["first_failing_level"] == "legs", r["line"]


@nx_required
def control_seven_line_top_sector_is_non_graph():
    """Row 6, corrected family (the six-line Sudakov ladder): the same family
    with all seven propagators in the denominator (top sector 127) must read
    the census's NON-GRAPH string (3V > 2E + N) -- a planted control."""
    man = load_manifest()
    vs = verdict_strings(man)["census"]
    e = _entry(man, "row06/cured_ladder_sudakovPR0")
    A, _ = _family_side(e)
    B, _ = _side(e["pair"]["drawn"])
    assert ta.compare_drawn(A, B)["verdict"] == e["pair"]["expected"]
    A7 = _copy(A, physical=[True] * len(A.propagators), name=A.name + "_sector127")
    r = ta.compare_drawn(A7, B)
    assert r["verdict"] == vs["NON-GRAPH"], r["line"]
    assert "family side" in r["reason"], r["reason"]


@nx_required
def control_undeclared_leg_classes_are_not_checkable_with_four_legs():
    """Row 1 (the outer-mass double box): the family with its leg virtualities
    withheld (four legs, no class declared) against the drawn graph that
    declares them must read NOT-CHECKABLE naming the family side and the legs
    -- nothing is inferred; the two-leg exemption does not apply."""
    man = load_manifest()
    vs = verdict_strings(man)
    e = _entry(man, "row01")
    A, notes = _family_side(e)
    B, _ = _side(e["pair"]["drawn"])
    assert A.leg_virt and notes["leg_virtualities_source"], notes
    assert ta.compare_drawn(A, B)["verdict"] == e["pair"]["expected"]
    A0 = _copy(A, leg_virt={}, name=A.name + "_noclasses")
    r = ta.compare_drawn(A0, B)
    assert "NOT-CHECKABLE" in vs["tool_only"]           # the census emitted none: the tool's word
    assert r["verdict"] == "NOT-CHECKABLE", r["line"]
    assert r["not_checkable"]["side"] == "family", r["not_checkable"]
    assert "undeclared" in r["not_checkable"]["reason"]


@nx_required
def control_two_point_undeclared_class_is_vacuous_and_said_by_name():
    """Row 13 (the equal-mass kite, two legs): the drawn yaml declares no p^2;
    the verdict is the census's IDENTITY-PASS with the leg-class note naming
    the drawn side and N = 2 (both legs carry one p^2 by conservation)."""
    man = load_manifest()
    e = _entry(man, "row13")
    A, _ = _family_side(e)
    B, _ = _side(e["pair"]["drawn"])
    assert len(A.exts) == 2 and not B.leg_virt
    r = ta.compare_drawn(A, B)
    assert r["verdict"] == e["pair"]["expected"], r["line"]
    assert r["leg_class_note"] and "drawn side" in r["leg_class_note"] and "N = 2" in r["leg_class_note"]


@nx_required
def control_cyclic_order_mismatch_is_a_named_sub_finding_that_keeps_the_verdict():
    """Row 21 (LBL3X): the drawn yaml declares the cyclic order p1,p2,p3,p4 and
    the census's map exchanges p2 and p4 (a reflection: preserved).  A family
    copy DECLARING the order p1,p3,p2,p4 maps to a non-cyclic image: the
    sub-finding CYCLIC-ORDER MISMATCH is named and the verdict stays the
    census's IDENTITY-PASS (the census's rule).  A family copy declaring the
    drawn order itself reads preserved True."""
    man = load_manifest()
    e = _entry(man, "row21")
    A, _ = _family_side(e)
    B, _ = _side(e["pair"]["drawn"])
    assert B.cyclic_leg_order, "the drawn yaml declares its cyclic order in the header"
    r = ta.compare_drawn(A, B)
    assert r["verdict"] == e["pair"]["expected"], r["line"]
    assert r["cyclic_order"]["A_declared"] is None and r["cyclic_order"]["drawn_order_in_A_names"]
    ok_order = list(B.cyclic_leg_order)
    Ap = _copy(A, cyclic_leg_order=ok_order)
    rp = ta.compare_drawn(Ap, B)
    assert rp["verdict"] == e["pair"]["expected"] and rp["cyclic_order"]["preserved"] is True, rp["line"]
    bad = [ok_order[0], ok_order[2], ok_order[1], ok_order[3]]
    An = _copy(A, cyclic_leg_order=bad)
    rn = ta.compare_drawn(An, B)
    assert rn["verdict"] == e["pair"]["expected"], rn["line"]
    assert rn["cyclic_order"]["preserved"] is False
    assert rn["cyclic_order"]["finding"] and rn["cyclic_order"]["finding"].startswith(ta.CYCLIC_ORDER_MISMATCH)
    assert ta.CYCLIC_ORDER_MISMATCH in rn["line"]


def _integral_evidence(res, tag):
    """The integral-choice evidence of one side: the four JSON keys of the
    side report and the printed line naming 'integral I of N' (None when the
    side prints none)."""
    s = res["sides"][tag]
    keys = {k: s.get(k) for k in ("integral_index", "n_integrals", "integral_indices", "other_integrals")}
    lines = [ln for ln in ta.verdict_lines(res) if ln.strip().startswith("integral ") and " of " in ln]
    return keys, lines


def _assert_integral_evidence(res, tag, want_index, want_n, ints_by_object):
    """Assert, BY NAME, that the side's evidence names index want_index of
    want_n with the record's own index vector and lists the others."""
    keys, lines = _integral_evidence(res, tag)
    assert keys["integral_index"] == want_index, (
        "%s side: evidence names integral index %r, expected %r" % (tag, keys["integral_index"], want_index))
    assert keys["n_integrals"] == want_n, (
        "%s side: evidence counts %r integrals, expected %r" % (tag, keys["n_integrals"], want_n))
    assert keys["integral_indices"] == ints_by_object[want_index], (
        "%s side: integral_indices %r != the record's integrals[%d].indices %r"
        % (tag, keys["integral_indices"], want_index, ints_by_object[want_index]))
    assert keys["other_integrals"] == [[i, v] for i, v in enumerate(ints_by_object) if i != want_index], (
        "%s side: other_integrals %r do not list the other %d index vectors by index"
        % (tag, keys["other_integrals"], want_n - 1))
    head = "integral %d of %d: %s; others: " % (want_index, want_n, ints_by_object[want_index])
    hit = [ln for ln in lines if ln.strip().startswith(head)]
    assert len(hit) == 1, ("printed evidence lacks the line '%s...' (lines: %r)" % (head, lines))
    for i, v in enumerate(ints_by_object):
        if i != want_index:
            assert "%d %s" % (i, v) in hit[0], ("the line does not list integral %d %s: %s" % (i, v, hit[0]))
    return hit[0]


@nx_required
def control_multi_integral_evidence_names_the_index_read():
    """Row 21 (LBL3X, an AMFlow-port JSON listing three integrals): the
    compare's evidence names WHICH integral the reader took and lists the
    others, so the multiplicity-level FAMILY-MISMATCH of the default reading
    is readable as 'the dotted integral' without the source.  Default (no
    --integral): index 0 of 3 (the integral with the most positive indices,
    first among ties -- the dotted line) with the record's integrals[0].indices
    and the others 1, 2 listed; the verdict is the tool's FAMILY-MISMATCH at
    the multiplicity level (the census's string for this row is reached only
    with its integral index).  --integral 2 (the census's call, the manifest's
    integral_index): index 2 of 3 and the census's IDENTITY-PASS.  The drawn
    side lists no integral and prints no such line; a one-integral record
    (row 28) prints no such line either.  The verdict STRINGS carry none of
    it (the line is evidence beneath the verdict, never part of it)."""
    man = load_manifest()
    vs = verdict_strings(man)["census"]
    e = _entry(man, "row21")
    fam_path = fixture(e["pair"]["family"])
    ints = [list(it["indices"]) for it in json.load(open(fam_path))["integrals"]]   # by object
    assert len(ints) == 3, ints
    ii = e["pair"]["integral_index"]
    assert ii == 2, e["pair"]["integral_index"]
    # default reading: index 0 of 3 (the dotted line), FAMILY-MISMATCH at the multiplicity level
    r0 = ta.compare_drawn(fam_path, fixture(e["pair"]["drawn"]))
    assert r0["verdict"] == vs["FAMILY-MISMATCH"] and r0["first_failing_level"] == "masses", r0["line"]
    line0 = _assert_integral_evidence(r0, "family", 0, 3, ints)
    assert "reader default" in line0, line0
    assert "integral 0 of 3" not in r0["line"] and " of 3" not in r0["line"], r0["line"]
    # the census's call: index 2 of 3, the census's string
    r2 = ta.compare_drawn(fam_path, fixture(e["pair"]["drawn"]), integral_indexA=ii)
    assert r2["verdict"] == e["pair"]["expected"], r2["line"]
    line2 = _assert_integral_evidence(r2, "family", 2, 3, ints)
    assert "--integral" in line2, line2
    assert " of 3" not in r2["line"], r2["line"]
    # the same through the import path the battery uses (a Family object)
    A, _ = _family_side(e)
    rA = ta.compare_drawn(A, fixture(e["pair"]["drawn"]))
    assert rA["verdict"] == e["pair"]["expected"]
    lineA = _assert_integral_evidence(rA, "family", 2, 3, ints)
    assert "--integral" in lineA, lineA
    # the drawn side (a routed yaml) lists no integral: keys present, nothing printed for it
    for r in (r0, r2):
        keys, lines = _integral_evidence(r, "drawn")
        assert keys["n_integrals"] == 0 and keys["integral_index"] is None and keys["other_integrals"] == []
        assert len(lines) == 1, lines
    # a one-integral record prints no such line (row 28: integrals lists one)
    e28 = _entry(man, "row28")
    p28 = fixture(e28["pair"]["family"])
    assert len(json.load(open(p28))["integrals"]) == 1
    r28 = ta.compare_drawn(p28, fixture(e28["pair"]["drawn"]))
    assert r28["verdict"] == e28["pair"]["expected"], r28["line"]
    keys28, lines28 = _integral_evidence(r28, "family")
    assert keys28["n_integrals"] == 1 and keys28["integral_index"] == 0 and keys28["other_integrals"] == []
    assert lines28 == [], lines28
    return {"expected": e["pair"]["expected"], "got": r2["verdict"],
            "line": "default: %s | --integral %d: %s" % (line0.strip()[:60], ii, line2.strip()[:60])}


@nx_required
def control_multi_integral_tampered_expectation_fails_by_name():
    """The planted negative of the evidence check: asserting the WRONG index
    (1 of 3 for row 21's default reading) must fail, and the failure names the
    index it found (0) and the one asserted (1); likewise a tampered count
    (2 integrals) and a tampered index vector."""
    man = load_manifest()
    e = _entry(man, "row21")
    fam_path = fixture(e["pair"]["family"])
    ints = [list(it["indices"]) for it in json.load(open(fam_path))["integrals"]]
    r0 = ta.compare_drawn(fam_path, fixture(e["pair"]["drawn"]))
    _assert_integral_evidence(r0, "family", 0, 3, ints)          # the true reading passes
    failed = []
    for label, args in (("index", (1, 3, ints)), ("count", (0, 2, ints[:2])),
                        ("vector", (0, 3, [[9] * len(ints[0])] + ints[1:]))):
        try:
            _assert_integral_evidence(r0, "family", *args)
        except AssertionError as exc:
            failed.append((label, str(exc)))
    assert [f[0] for f in failed] == ["index", "count", "vector"], failed
    assert "names integral index 0, expected 1" in failed[0][1], failed[0]
    assert "counts 3 integrals, expected 2" in failed[1][1], failed[1]
    assert "integrals[0].indices" in failed[2][1], failed[2]
    return {"expected": "3 tampered expectations fail by name", "got": "%d failed by name" % len(failed)}


@nx_required
def unit_legs_per_vertex_reads_the_row27_photon_counts():
    """legs_per_vertex on the two row-27 families gives the census's photons
    per vertex ([2, 1, 1, 0] for the replaced three-point family, [1, 1, 1, 1]
    for the corrected box family)."""
    man = load_manifest()
    for fam_key, want in ROW27_PHOTONS_PER_VERTEX.items():
        e = _entry(man, "row27/" + fam_key)
        A, _ = _family_side(e)
        lpv = ta.legs_per_vertex(A)
        assert lpv["ok"] and lpv["legs_per_vertex"] == want, (fam_key, lpv)
        assert lpv["realizer_legs_per_vertex"] == want, (fam_key, lpv)
        assert lpv["V"] == 4 and lpv["E"] == 6 and lpv["N"] == 4, lpv


def unit_drawn_header_forms_read_by_object():
    """read_drawn_header on the census's four header forms gives the legs the
    body's family reads after completion, the virtualities' classes matching
    the family side's, and the multiplicity list equal to the record's nu
    multiset where the header carries one (rows 29 / 33)."""
    man = load_manifest()
    seen = set()
    for key in ("row24/cured_adjacent_p3p4", "row21", "row15", "row19", "row29", "row33"):
        e = _entry(man, key)
        hdr = ta.read_drawn_header(fixture(e["pair"]["drawn"]))
        assert hdr and hdr["form"], key
        seen.add(hdr["form"])
        B, notes = _side(e["pair"]["drawn"], e["pair"]["drawn_kinematics"])
        if e["pair"]["drawn_kinematics"] is None:
            assert "legs" in notes["header_applied"], (key, notes)
            assert sorted(hdr["legs"]) == sorted(B.exts), (key, hdr["legs"], B.exts)
            assert all(B.leg_virt.get(lg) == v for lg, v in hdr["virtualities"].items()), (key, B.leg_virt)
        if hdr.get("nu") is not None:
            A, _ = _family_side(e)
            assert sorted(hdr["nu"]) == ta.labeled_set(A)["nu_multiset"], (key, hdr["nu"])
            assert "nu" in notes["header_applied"]
    assert {"external legs", "legs@vertex", "leg dictionary", "External legs"} <= seen, seen


def _per_row_counts(man):
    """The census's own count rule re-derived from the manifest's entries: a
    row's family-level entries replace its row-level entry (the manifest's
    verdict_count_convention.summary_verdict_counts_rule)."""
    by_row = {}
    for e in man["entries"]:
        by_row.setdefault(e["row"], {"row": [], "family": []})[e["level"]].append(e["pair"]["expected"])
    counts = {}
    for row, d in by_row.items():
        for v in (d["family"] or d["row"]):
            counts[v] = counts.get(v, 0) + 1
    return counts


def _pair_counts(man):
    """The battery's pair counts: the expected string of every entry (the
    manifest's verdict_count_convention.pair_verdict_counts rule)."""
    got = {}
    for e in man["entries"]:
        got[e["pair"]["expected"]] = got.get(e["pair"]["expected"], 0) + 1
    return got


def _rollup_bucket(man):
    """The rollup bucket of the census's summary.entry_counts, by object: the
    one key of the manifest's copy that is not a verdict string of
    summary.verdict_counts."""
    conv = man["verdict_count_convention"]
    extra = sorted(set(conv["summary_entry_counts"]) - set(man["summary_verdict_counts"]))
    assert len(extra) == 1, (extra, conv["summary_entry_counts"], man["summary_verdict_counts"])
    return extra[0]


def unit_manifest_is_the_census_by_object():
    """The manifest carries every entry of the census (the family-level pairs
    + the rollups = n_entries_in_rows, over n_rows rows), the census's
    counting rule, its verdict-count summary under the per-row rule
    (re-derived from the entries), its entry counts (re-derived: the pairs +
    the rollups in their own bucket), and every rollup string joined from
    the family verdicts in the manifest's order."""
    man = load_manifest()
    assert man["n_entries_in_rows"] == len(man["entries"]) + len(man["rollups"])
    assert man["n_rows"] == len({e["row"] for e in man["entries"]} | {r["row"] for r in man["rollups"]})
    assert man["n_family_level"] == len(man["entries"]) and man["n_rollups"] == len(man["rollups"])
    conv = man["verdict_count_convention"]
    counts = man["summary_verdict_counts"]
    assert isinstance(conv["counting_rule"], str) and conv["counting_rule"].strip(), conv
    pair = _pair_counts(man)
    assert pair == conv["pair_verdict_counts"], (pair, conv["pair_verdict_counts"])
    per_row = _per_row_counts(man)
    assert per_row == counts == conv["per_row_verdict_counts_rederived"], (per_row, counts, conv)
    bucket = _rollup_bucket(man)
    assert bucket == conv["rollup_bucket"], (bucket, conv["rollup_bucket"])
    ent = dict(pair)
    ent[bucket] = len(man["rollups"])
    assert ent == conv["summary_entry_counts"] == conv["entry_counts_rederived"], (ent, conv)
    for r in man["rollups"]:
        fams = [_entry(man, fk)["pair"]["expected"] for fk in r["families_in_order"]]
        assert r["verdict"] == "see families: " + " / ".join(fams), (r["key"], r["verdict"], fams)
    for e in man["entries"]:
        assert e["pair"]["expected"] in counts, (e["key"], e["pair"]["expected"])
    return {"expected": "manifest counts re-derived", "got": "%d entries + %d rollups over %d rows; %s" % (
        len(man["entries"]), len(man["rollups"]), man["n_rows"], dict(per_row))}


def _all_fixture_infos(man):
    """Every fixture record the manifest's entries name (pairs, extras, alsos)."""
    out = []
    for e in man["entries"]:
        p = e["pair"]
        out += [p[k] for k in ("family", "kinematics", "drawn", "drawn_kinematics") if p.get(k)]
        for x in e.get("extras", []):
            out += [x[k] for k in ("drawn", "drawn_kinematics") if x.get(k)]
        for a in e.get("also", []):
            out += [a[k] for k in ("family", "kinematics", "drawn") if a.get(k)]
    return out


def unit_rows_not_compared_have_no_entry_rollup_or_fixture():
    """A corpus row listed under rows_not_compared (no drawn figure in the
    paper, hence no compare verdict) has NO entry, NO rollup and NO fixture
    in the manifest: no fixture path any entry names starts with its row
    directory.  Its files may still ship under tests/fixtures/rowNN/ for the
    other batteries."""
    man = load_manifest()
    blk = man["rows_not_compared"]
    rows = blk["rows"]
    assert isinstance(rows, list) and rows and blk["reason"], blk
    fixtures = [info["fixture"] for info in _all_fixture_infos(man)]
    for r in rows:
        assert not [e for e in man["entries"] if e["row"] == r], f"row {r} has an entry"
        assert not [x for x in man["rollups"] if x["row"] == r], f"row {r} has a rollup"
        rowdir = "row%02d" % int(r)
        hits = [f for f in fixtures if f.startswith(rowdir + "/")]
        assert not hits, f"row {r}: fixtures named under {rowdir}/: {hits}"
    all_rows = {e["row"] for e in man["entries"]} | {x["row"] for x in man["rollups"]}
    assert not (set(rows) & all_rows), (rows, sorted(all_rows))
    return {"expected": "rows not compared are absent from the manifest",
            "got": "%d row(s) %s: no entry / rollup / fixture" % (len(rows), rows)}


@nx_required
def control_named_family_is_load_bearing():
    """Row 18's record family (a file listing families 2 and 5 of the
    reference; the census judged family 2): named as the census's family the
    pair reads the census's string; named as the OTHER family of the file
    (ggH_fam5, whose first top-level sector 31 has five lines: a shape
    mismatch in loops / legs / line count) it reads the census's
    FAMILY-MISMATCH string at the bare level, so the name is load-bearing,
    never a default; a name the file does not carry is refused by name."""
    man = load_manifest()
    vs = verdict_strings(man)["census"]
    e = [x for x in man["entries"] if x["pair"].get("family_name")]
    assert e, "no manifest entry names its family"
    e = [x for x in e if x["key"] == "row18/record_family"][0]
    p = e["pair"]
    names = p["families_in_file"]
    assert p["family_name"] in names and len(names) >= 2, (p["family_name"], names)
    got, res, _ = run_pair(p["family"], p["kinematics"], p["drawn"], p["drawn_kinematics"], None, p["family_name"])
    assert got == p["expected"], res["line"]
    other = [n for n in names if n != p["family_name"]][0]
    got2, res2, _ = run_pair(p["family"], p["kinematics"], p["drawn"], p["drawn_kinematics"], None, other)
    assert res2["sides"]["family"]["name"] == other
    assert got2 != p["expected"] and got2 == vs["FAMILY-MISMATCH"], res2["line"]
    assert res2["first_failing_level"] == "bare", res2["line"]
    assert res2["E"]["A"] != res2["E"]["B"], res2["line"]           # the line counts differ (5 vs 6)
    try:
        run_pair(p["family"], p["kinematics"], p["drawn"], p["drawn_kinematics"], None, "no_such_family")
    except ValueError as exc:
        assert "no family found" in str(exc), exc
    else:
        raise AssertionError("a family name the file does not carry was not refused")
    return {"expected": p["expected"], "got": got,
            "line": "%s -> %s; %s -> %s (%s)" % (p["family_name"], got, other, got2, res2["first_failing_level"])}


def build_cases():
    man = load_manifest()
    cases = [_entry_case(e) for e in man["entries"]]
    cases += [_rollup_case(r) for r in man["rollups"]]
    cases += [control_mass_moved_flips_to_family_mismatch_at_the_mass_level,
              control_offshell_leg_moved_flips_to_family_mismatch_at_the_leg_level,
              control_seven_line_top_sector_is_non_graph,
              control_undeclared_leg_classes_are_not_checkable_with_four_legs,
              control_two_point_undeclared_class_is_vacuous_and_said_by_name,
              control_cyclic_order_mismatch_is_a_named_sub_finding_that_keeps_the_verdict,
              control_multi_integral_evidence_names_the_index_read,
              control_multi_integral_tampered_expectation_fails_by_name,
              control_named_family_is_load_bearing,
              unit_legs_per_vertex_reads_the_row27_photon_counts,
              unit_drawn_header_forms_read_by_object,
              unit_manifest_is_the_census_by_object,
              unit_rows_not_compared_have_no_entry_rollup_or_fixture]
    return cases


try:
    BATTERY = build_cases()
except FixtureMissing as _exc:                                      # pragma: no cover
    BATTERY = []
    _BATTERY_MISSING = str(_exc)


def run_battery(verbose=True):
    """Run every case; return [{"name", "status": PASS|FAIL|SKIP, "detail",
    "wall_s", "expected", "got"}].  A missing fixture is SKIP (named); drift or
    any other exception is FAIL."""
    rows = []
    if not BATTERY:
        rows.append({"name": "census_battery", "status": "SKIP",
                     "detail": globals().get("_BATTERY_MISSING", "no cases"), "wall_s": 0.0})
        if verbose:
            print("   SKIP  census_battery  (%s)" % rows[-1]["detail"], flush=True)
        return rows
    t_all = time.time()
    for fn in BATTERY:
        name = fn.__name__
        t0 = time.time()
        row = {"name": name, "status": "PASS", "detail": "", "expected": None, "got": None}
        try:
            r = fn()
            if isinstance(r, dict):
                row["expected"], row["got"] = r.get("expected"), r.get("got")
                row["detail"] = (r.get("line") or "")[:160]
                for x in r.get("extras", []) or []:
                    row["detail"] += " | extra %s -> %s" % (os.path.basename(x["drawn"]), x["got"])
                for a in r.get("also", []) or []:
                    row["detail"] += " | also %s -> %s" % (os.path.basename(a["drawn"]), a["got"])
        except FixtureMissing as exc:
            row.update({"status": "SKIP", "detail": str(exc)})
        except AssertionError as exc:
            row.update({"status": "FAIL", "detail": "AssertionError: %s" % exc})
        except Exception as exc:
            row.update({"status": "FAIL", "detail": f"{type(exc).__name__}: {exc}"})
        row["wall_s"] = round(time.time() - t0, 3)
        rows.append(row)
        if verbose:
            print(f"   {row['status']:5s} {name}  [{row['wall_s']} s]"
                  + (f"  ({row['detail']})" if row["detail"] else ""), flush=True)
    if verbose:
        print("   census battery wall %.2f s" % (time.time() - t_all), flush=True)
    return rows


if __name__ == "__main__":
    out = run_battery(verbose=True)
    n_fail = sum(1 for r in out if r["status"] == "FAIL")
    n_skip = sum(1 for r in out if r["status"] == "SKIP")
    n_pass = sum(1 for r in out if r["status"] == "PASS")
    print(f"census battery: {len(out)} cases, {n_pass} PASS, {n_fail} FAIL, {n_skip} SKIP")
    sys.exit(1 if n_fail else (2 if n_skip else 0))
