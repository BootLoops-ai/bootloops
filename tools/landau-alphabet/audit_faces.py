#!/usr/bin/env python3
r"""
audit_faces.py — GUARD against SILENTLY-DROPPED Landau / PLD faces.

Motivation: a hand-written COMPLETE_LANDAU_VARIETY.json / COMPLETE_ALPHABET_FACES.json
can declare itself complete while only a fraction of the kinematic faces were
recorded and the hard faces (the ones that timed out) were never carried as a
NAMED gap — e.g. on the 3-loop K4 graph, 11 of 54 faces recorded with the
leading K4 face and three 5-edge sub-boxes silently missing.  The hard faces
are exactly where the missing letters live.  A truncation that reads as
"covered everything" is an honesty failure mode.

This guard re-runs the face audit on any Landau Alphabet-style output JSON and
HARD-FAILS (exit code 2) on any silent drop.  It handles THREE schemas:

  (1) ENGINE output  (landau_alphabet.py / pld_bridge.jl): has a "faces" list and
      "honesty".  We recompute the COMPLETE flag from the per-face statuses and the
      embedded graph spec (expected_faces), and flag any mismatch with what the file
      asserts (e.g. an honesty block claiming 0 timeouts while a face record is
      timed_out, or attempted < expected => enumeration truncated).

  (2) LEGACY / hand-written analysis JSON (COMPLETE_LANDAU_VARIETY.json,
      COMPLETE_ALPHABET_FACES.json): the "graph" is prose, "honesty" is prose, and
      there is NO machine-checkable per-face list — or a "total_faces" that counts
      only the hand-recorded strata, not the true face lattice.  If a graph spec can
      be located/supplied, we compute the TRUE expected face count and flag the gap.
      Even with no spec, a file whose NAME or keys claim "COMPLETE" but which carries
      no auditable per-face accounting is flagged as UNVERIFIABLE (soft fail unless
      --strict).

  (3) Anything else: reported as "not a face file", skipped.

A face is a SILENT DROP iff it is missing/incomplete AND the file does not surface
it (complete flag still true, or no unresolved list).  The whole point: after this
guard, "complete": true is a CHECKED claim, not a hope.

USAGE
    python audit_faces.py FILE.json [FILE2.json ...] [--spec SPEC.json] [--strict]
    python audit_faces.py --glob 'results/**/COMPLETE_*.json' --strict

EXIT CODES
    0  all audited files are complete (or cleanly flagged-incomplete with no silent drop)
    2  at least one SILENT DROP found (file claims complete but isn't / hidden gap)
    3  at least one file is UNVERIFIABLE under --strict (claims complete, no audit data)
"""
from __future__ import annotations
import json, sys, os, glob, argparse, re
from itertools import combinations

# Reuse the engine's Symanzik face census so the expected count is computed by the
# SAME code the engine enumerates with (no second, drifting implementation).
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    from landau_alphabet import GraphSpec, symanzik_subgraph
except Exception:                                       # pragma: no cover
    GraphSpec = None
    symanzik_subgraph = None


# ---------------------------------------------------------------------------
def expected_face_count(spec: dict, min_face_edges: int = 2) -> int | None:
    """True number of kinematic-or-scaleless faces (connected sub-multigraphs,
    >=min_face_edges edges) of the graph in `spec`.  None if engine unavailable
    or spec is not a real edge-list graph spec."""
    if GraphSpec is None or not isinstance(spec, dict) or "edges" not in spec:
        return None
    try:
        gs = GraphSpec(spec)
    except Exception:
        return None
    ALLE = gs.ALLE
    n = 0
    for k in range(len(ALLE), min_face_edges - 1, -1):
        for face in combinations(ALLE, k):
            if symanzik_subgraph(gs, set(face)) is not None:
                n += 1
    return n


def _looks_like_spec(g) -> bool:
    return isinstance(g, dict) and "edges" in g and "nodes" in g


# ---------------------------------------------------------------------------
def audit_engine_output(d: dict, fname: str, spec_override: dict | None):
    """Audit a landau_alphabet.py / pld_bridge.jl style output. Returns
    (ok: bool, silent_drop: bool, lines: list[str])."""
    lines = []
    faces = d.get("faces", [])
    honesty = d.get("honesty", {}) or {}

    # 1. recompute per-face status tally from the ACTUAL face records (ground truth)
    n_done = n_timeout = n_overbudget = n_scaleless = n_unknown = 0
    incomplete_faces = []
    for rec in faces:
        st = rec.get("face_status")
        if st is None:                       # older engine output: derive from "status"
            s = str(rec.get("status", "")).lower()
            if s.startswith("resolved") or s in ("ok", "done"):
                st = "done"                  # incl. 'resolved-no-letter'
            elif "timeout" in s or "timed_out" in s:
                st = "timed_out"
            elif "scaleless" in s or "no-kinematic" in s:
                st = "scaleless"
            elif "partial" in s:
                st = "over_budget"
            elif s in ("queued", "pending", "running") or "pending" in s or "running" in s:
                st = "pending"               # not yet attempted/finished == an OPEN gap
            elif s in ("candidate",) or "candidate" in s:
                st = "candidate"             # unverified guess == an OPEN gap
            else:
                st = "unknown"
        if st == "done":
            n_done += 1
        elif st == "timed_out":
            n_timeout += 1; incomplete_faces.append((rec.get("face_edges"), "timed_out"))
        elif st == "over_budget":
            n_overbudget += 1; incomplete_faces.append((rec.get("face_edges"), "over_budget"))
        elif st == "scaleless":
            n_scaleless += 1
        else:                                 # pending / candidate / unknown == open gap
            n_unknown += 1
            label = st if st in ("pending", "candidate") else "unknown-status"
            ident = rec.get("face_edges") or rec.get("face") or rec.get("face_id")
            incomplete_faces.append((ident, label))

    n_attempted = len(faces)
    n_incomplete = n_timeout + n_overbudget + n_unknown

    # 2. expected face census from the embedded (or supplied) graph spec
    spec = spec_override or (d.get("graph") if _looks_like_spec(d.get("graph")) else None)
    expected = expected_face_count(spec) if spec is not None else None
    n_missing = (expected - n_attempted) if expected is not None else None

    lines.append(f"  schema: ENGINE output ({n_attempted} face records)")
    lines.append(f"  faces:  done={n_done} timed_out={n_timeout} over_budget="
                 f"{n_overbudget} scaleless={n_scaleless} unknown={n_unknown}")
    if expected is not None:
        lines.append(f"  expected faces (recomputed from graph spec): {expected}"
                     + (f"  -> MISSING {n_missing} NEVER-ATTEMPTED" if n_missing else ""))
    else:
        lines.append("  expected faces: UNKNOWN (no embeddable graph spec; pass --spec to census)")

    # 3. truly complete?
    deadline = bool(honesty.get("global_deadline_hit", False))
    truly_complete = (n_incomplete == 0
                      and (n_missing in (0, None))
                      and not deadline)
    # but if we couldn't census, a "complete" with no missing-info is only PROVISIONAL
    census_known = expected is not None

    if "complete" in d:
        claimed_complete = bool(d.get("complete"))
    elif "complete" in honesty:
        claimed_complete = bool(honesty.get("complete"))
    else:
        # No explicit flag (older engine / hand-rolled output). The file still
        # makes an IMPLICIT completeness claim unless its honesty block names a
        # gap. A clean/absent honesty block over incomplete faces == silent drop.
        named_gap = (honesty.get("bounded_timeout", 0) or honesty.get("bounded_partial", 0)
                     or honesty.get("unresolved_faces") or honesty.get("missing_faces"))
        claimed_complete = not named_gap

    silent_drop = False
    # (a) claims complete but has incomplete/missing faces
    if claimed_complete and not truly_complete:
        silent_drop = True
        lines.append("  *** SILENT DROP: file claims complete but faces are incomplete/missing:")
        for fe, why in incomplete_faces:
            lines.append(f"        - {why}: {fe}")
        if n_missing:
            lines.append(f"        - {n_missing} face(s) NEVER ATTEMPTED (enumeration truncated)")
    # (b) honesty tally disagrees with the actual face records
    h_to = honesty.get("bounded_timeout")
    if h_to is not None and h_to != n_timeout:
        silent_drop = True
        lines.append(f"  *** TALLY MISMATCH: honesty.bounded_timeout={h_to} but "
                     f"{n_timeout} face record(s) are timed_out")
    h_nf = honesty.get("n_faces")
    if h_nf is not None and census_known and h_nf != expected and claimed_complete:
        # an n_faces that equals attempted but not expected, while claiming complete,
        # is the classic "total_faces counts only what I recorded" truncation
        lines.append(f"  *** NOTE: honesty.n_faces={h_nf} != expected={expected} "
                     f"(file counts recorded faces as the total)")

    ok = truly_complete and not silent_drop and census_known
    if truly_complete and not census_known:
        lines.append("  status: incomplete-info (looks complete but face census unverified;"
                     " pass --spec to confirm)")
    elif ok:
        lines.append("  status: COMPLETE (verified: every expected face attempted & resolved)")
    elif not silent_drop:
        lines.append(f"  status: INCOMPLETE but HONESTLY FLAGGED "
                     f"({n_incomplete} unresolved, {n_missing or 0} missing) -- no silent drop")
    return ok, silent_drop, census_known, lines


# ---------------------------------------------------------------------------
def audit_legacy(d: dict, fname: str, spec_override: dict | None, strict: bool):
    """Audit a hand-written 'COMPLETE_*' analysis JSON that has no machine-checkable
    per-face list.  Returns (ok, silent_drop, unverifiable, lines)."""
    lines = ["  schema: LEGACY / hand-written analysis (no auditable per-face list)"]
    base = os.path.basename(fname).lower()
    claims_complete = ("complete" in base) or ("complete" in {k.lower() for k in d})

    # how many faces/strata did it actually record?
    recorded = None
    for key in ("FACES", "faces", "strata"):
        v = _dig(d, key)
        if isinstance(v, (list, dict)):
            recorded = len(v); break
    # did its (prose) honesty leak a "total_faces"?
    leaked_total = _find_int_near(d, "total_faces")
    leaked_to    = _find_int_near(d, "bounded_out") or _find_int_near(d, "bounded_timeout")

    # try to census the true face count
    spec = spec_override
    if spec is None:
        spec = _find_spec_for(fname)        # look for a sibling graph spec
    expected = expected_face_count(spec) if spec is not None else None

    if recorded is not None:
        lines.append(f"  recorded faces/strata: {recorded}")
    if leaked_total is not None:
        lines.append(f"  file's asserted total_faces: {leaked_total}")
    if expected is not None:
        lines.append(f"  TRUE expected faces (recomputed from {os.path.basename(spec.get('__src__','spec'))}): {expected}")

    silent_drop = False
    # the silent-drop failure mode: claims complete, true expected >> recorded/leaked
    cmp_recorded = leaked_total if leaked_total is not None else recorded
    if claims_complete and expected is not None and cmp_recorded is not None \
       and cmp_recorded < expected:
        silent_drop = True
        lines.append(f"  *** SILENT DROP: claims COMPLETE but only {cmp_recorded} of "
                     f"{expected} faces accounted for ({expected - cmp_recorded} dropped).")
        # surface the prose admission if present
        prose = json.dumps(d.get("honesty", ""), default=str)
        for m in re.findall(r"timed out on[^\"]+", prose):
            lines.append(f"        prose admission: {m[:200]}")

    unverifiable = False
    if claims_complete and expected is None:
        unverifiable = True
        lines.append("  *** UNVERIFIABLE: claims COMPLETE but has no auditable per-face list "
                     "and no graph spec found to census against. Provide --spec.")

    ok = (not silent_drop) and (not unverifiable)
    if ok and claims_complete and expected is not None:
        lines.append("  status: COMPLETE (recorded count consistent with face census)")
    elif ok:
        lines.append("  status: not a completeness claim / nothing to falsify")
    return ok, silent_drop, unverifiable, lines


# ---------------------------------------------------------------------------
def _dig(d, key):
    """Find first value for `key` anywhere in nested dict."""
    if isinstance(d, dict):
        if key in d:
            return d[key]
        for v in d.values():
            r = _dig(v, key)
            if r is not None:
                return r
    return None


def _find_int_near(d, name):
    """Find an int field whose key contains `name`, anywhere nested."""
    if isinstance(d, dict):
        for k, v in d.items():
            if name in str(k).lower() and isinstance(v, int):
                return v
            r = _find_int_near(v, name)
            if r is not None:
                return r
    elif isinstance(d, list):
        for v in d:
            r = _find_int_near(v, name)
            if r is not None:
                return r
    return None


def _find_spec_for(fname):
    """Look for a graph spec JSON next to / under the Landau Alphabet graph_specs
    dir whose name matches the analysis file's apparent graph (best-effort)."""
    here = os.path.dirname(os.path.abspath(__file__))
    specdir = os.path.join(here, "graph_specs")
    base = os.path.basename(fname).lower()
    if not os.path.isdir(specdir):
        return None
    # heuristic: match the graph name token (e.g. 'topbox') in the analysis filename
    # or its path; fall back to scanning all specs for a name substring hit.
    path = fname.lower()
    for sp in sorted(glob.glob(os.path.join(specdir, "*.json"))):
        nm = os.path.splitext(os.path.basename(sp))[0].lower()
        if nm in base or nm in path:
            try:
                s = json.load(open(sp)); s["__src__"] = sp
                return s
            except Exception:
                pass
    return None


# ---------------------------------------------------------------------------
def audit_file(fname, spec_override, strict):
    try:
        d = json.load(open(fname))
    except Exception as e:
        print(f"[{fname}]\n  ERROR: cannot parse JSON ({e})")
        return (False, False, True)        # treat as unverifiable
    print(f"[{fname}]")
    if not isinstance(d, dict):
        print("  schema: not a face file (top-level JSON is not an object); skipped")
        return (True, False, False)
    has_face_list = isinstance(d.get("faces"), list) and d.get("faces")
    if has_face_list and isinstance(d["faces"][0], dict) and (
            "face_status" in d["faces"][0] or "status" in d["faces"][0]
            or "face_edges" in d["faces"][0] or "face_id" in d["faces"][0]):
        ok, silent, census_known, lines = audit_engine_output(d, fname, spec_override)
        unverifiable = (not census_known) and strict
    else:
        ok, silent, unverifiable, lines = audit_legacy(d, fname, spec_override, strict)
    for L in lines:
        print(L)
    return (ok, silent, unverifiable)


def main(argv=None):
    ap = argparse.ArgumentParser(description="Guard against silently-dropped Landau faces.")
    ap.add_argument("files", nargs="*", help="output/analysis JSON(s) to audit")
    ap.add_argument("--glob", action="append", default=[], help="glob pattern(s) to audit")
    ap.add_argument("--spec", default=None, help="graph spec JSON to census faces against")
    ap.add_argument("--strict", action="store_true",
                    help="treat 'claims complete but unverifiable' as a failure")
    args = ap.parse_args(argv)

    spec_override = None
    if args.spec:
        spec_override = json.load(open(args.spec)); spec_override["__src__"] = args.spec

    targets = list(args.files)
    for g in args.glob:
        targets += sorted(glob.glob(g, recursive=True))
    if not targets:
        ap.error("no files given (positional or --glob)")

    n_silent = n_unverif = n_ok = 0
    for f in targets:
        ok, silent, unverif = audit_file(f, spec_override, args.strict)
        print()
        if silent:
            n_silent += 1
        elif unverif:
            n_unverif += 1
        elif ok:
            n_ok += 1

    print("=" * 60)
    print(f"audited {len(targets)} file(s): {n_ok} complete, "
          f"{n_silent} SILENT DROP, {n_unverif} unverifiable")
    if n_silent:
        return 2
    if n_unverif and args.strict:
        return 3
    return 0


if __name__ == "__main__":
    sys.exit(main())
