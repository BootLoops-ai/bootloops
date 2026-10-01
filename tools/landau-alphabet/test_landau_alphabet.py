#!/usr/bin/env python3
"""
Regression tests for landau_alphabet.py on SMALL graphs with known alphabets.
Run:  python test_landau_alphabet.py
"""
import os, sys, json
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from landau_alphabet import load_spec, run, canonical_alphabet  # noqa: E402


def _letters_match(got, expected_exprs):
    """True iff every expected expression appears (up to a unit) in `got`."""
    got_e = [sp.sympify(g) for g in got]
    for e in expected_exprs:
        e = sp.sympify(e)
        ok = any(sp.cancel(e / g).is_number for g in got_e if g != 0)
        if not ok:
            return False, str(e)
    return True, None


def test_box1l():
    """1-loop massless on-shell box: PLD alphabet = {s, t} (leading sing 1/(st);
       no u-letter — matches the PLD.jl reference output)."""
    gs = load_spec(os.path.join(HERE, "graph_specs", "box1l.json"))
    log = run(gs, out_path=None, face_timeout=20, global_minutes=5)
    ok, miss = _letters_match(log["alphabet"], ["s", "t"])
    assert ok, f"box1l: missing letter {miss}; got {log['alphabet']}"
    assert set(log["alphabet"]) == {"s", "t"}, f"box1l: extra letters {log['alphabet']}"
    # honesty: every face must be fully resolved on this tiny graph
    assert log["honesty"]["bounded_timeout"] == 0
    assert log["honesty"]["bounded_partial"] == 0
    # top-level complete flag is TRUE only when every expected face is resolved
    assert log["complete"] is True
    assert log["honesty"]["missing_faces"] == 0
    assert log["honesty"]["attempted_faces"] == log["honesty"]["expected_faces"]
    # first-entry should include the bare Mandelstams
    assert "s" in log["first_entry"] and "t" in log["first_entry"]
    print(f"[PASS] box1l  alphabet={log['alphabet']}  "
          f"first_entry={log['first_entry']}  honesty={log['honesty']['resolved']}/"
          f"{log['honesty']['n_faces']} resolved")
    return log


def test_sunrise2l():
    """2-loop equal-mass sunrise: Landau loci = {s, s-m2, s-9m2} (normal+pseudo thresholds)."""
    gs = load_spec(os.path.join(HERE, "graph_specs", "sunrise2l.json"))
    log = run(gs, out_path=None, face_timeout=20, global_minutes=5)
    ok, miss = _letters_match(log["alphabet"], ["s - m2", "s - 9*m2"])
    assert ok, f"sunrise2l: missing letter {miss}; got {log['alphabet']}"
    assert log["honesty"]["bounded_timeout"] == 0
    print(f"[PASS] sunrise2l  alphabet={log['alphabet']}  "
          f"honesty={log['honesty']['resolved']}/{log['honesty']['n_faces']} resolved")
    return log


def test_incremental_flush(tmp=None):
    """Output JSON is written incrementally and matches the in-memory result."""
    if tmp is None:
        import tempfile
        tmp = os.path.join(tempfile.mkdtemp(prefix="lb_test_"), "box1l.json")
    gs = load_spec(os.path.join(HERE, "graph_specs", "box1l.json"))
    log = run(gs, out_path=tmp, face_timeout=20, global_minutes=5)
    with open(tmp) as fh:
        ondisk = json.load(fh)
    assert ondisk["alphabet"] == log["alphabet"]
    assert ondisk["honesty"]["n_faces"] == log["honesty"]["n_faces"]
    print(f"[PASS] incremental flush -> {tmp}")


def test_k4_incomplete_is_flagged():
    """Honesty guard for the bounded sympy chain WITHOUT the Groebner backend:
    with PLD_GROEBNER=off and a SHORT face timeout the 4 hard K4-graph faces
    (leading K4 + three 5-edge sub-boxes) time out and MUST surface as named
    gaps with complete==False. Guards the silent-drop failure mode.

    NOTE: when the msolve/Singular Groebner backend IS present (default), those
    faces now RESOLVE via escalation -- see test_k4_groebner_completes. This
    test pins the *fallback* behavior, so it forces the backend off."""
    old = os.environ.get("PLD_GROEBNER")
    os.environ["PLD_GROEBNER"] = "off"
    try:
        import importlib, landau_alphabet as _L
        importlib.reload(_L)
        gs = _L.load_spec(os.path.join(HERE, "graph_specs", "k4box3l.json"))
        log = _L.run(gs, out_path=None, face_timeout=6, global_minutes=8)
        assert log["complete"] is False, "K4 short-timeout must NOT claim complete"
        assert log["honesty"]["expected_faces"] == 54, log["honesty"]["expected_faces"]
        assert log["honesty"]["missing_faces"] == 0
        assert log["honesty"]["bounded_timeout"] >= 1
        uf = [tuple(sorted(f)) for f in log["honesty"]["unresolved_faces"]]
        assert tuple(sorted(gs.ALLE)) in uf, "leading K4 face must be a named gap"
        print(f"[PASS] k4box3l incomplete flagged (backend OFF): complete={log['complete']} "
              f"timeout={log['honesty']['bounded_timeout']}/{log['honesty']['expected_faces']} faces")
    finally:
        if old is None:
            os.environ.pop("PLD_GROEBNER", None)
        else:
            os.environ["PLD_GROEBNER"] = old
        import importlib, landau_alphabet as _L
        importlib.reload(_L)


def test_k4_groebner_completes():
    """With the Groebner backend present, the leading K4 face that the
    bounded sympy chain TIMES OUT on must now RESOLVE (msolve f4 saturated
    elimination, ~32 s). Skips if neither msolve nor Singular is installed."""
    import landau_alphabet as _L
    if _L._GB is None:
        print("[SKIP] no msolve/Singular backend available")
        return
    gs = _L.load_spec(os.path.join(HERE, "graph_specs", "k4box3l.json"))
    sg = _L.symanzik_subgraph(gs, set(gs.ALLE))     # leading K4 face
    fv = [gs.ev[e] for e in gs.ALLE]
    lf, st = _L.discriminant_letters(gs, sg["F"], fv, 6)  # 6s sympy budget -> escalate
    assert "groebner" in st, f"expected Groebner escalation, got {st}"
    assert lf, "K4 leading face must yield letters"
    print(f"[PASS] k4box3l K4 leading completes via {st}: {len(lf)} letters")


def test_auditor_catches_silent_drop():
    """The guard (audit_faces.py) must exit-2 on a file that claims complete but
    has a timed_out face record, and on a legacy 'COMPLETE_*' that records far
    fewer faces than the true face lattice."""
    import subprocess, tempfile
    audit = os.path.join(HERE, "audit_faces.py")
    spec = os.path.join(HERE, "graph_specs", "k4box3l.json")
    # (1) engine file with a LIE: complete=true but a face is timed_out
    gs = load_spec(spec)
    log = run(gs, out_path=None, face_timeout=6, global_minutes=8)
    # Construct the LIE: summary says complete but a face record is timed_out.
    # (Before the resolved-groebner status fix, the engine itself produced this
    # inconsistency; now we inject it explicitly to test the auditor.)
    log["faces"][-1]["face_status"] = "timed_out"
    log["complete"] = True
    log["honesty"]["complete"] = True
    log["honesty"]["bounded_timeout"] = 0
    log["honesty"]["unresolved_faces"] = []
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
        json.dump(log, fh, default=str); lie = fh.name
    r = subprocess.run([sys.executable, audit, lie, "--spec", spec],
                       capture_output=True, text=True)
    assert r.returncode == 2, f"auditor must exit 2 on tally-lie; got {r.returncode}\n{r.stdout}"
    assert "SILENT DROP" in r.stdout
    # (2) legacy hand-written COMPLETE_* that drops faces
    legacy = {"task": "x", "graph": "K4 prose", "FACES": [1, 2, 3],
              "COMPLETE_LANDAU_VARIETY": {"strata": {"a": 1, "b": 2}},
              "honesty": ["all covered"]}
    with tempfile.NamedTemporaryFile("w", suffix="COMPLETE_k4box3l.json", delete=False) as fh:
        json.dump(legacy, fh); leg = fh.name
    r2 = subprocess.run([sys.executable, audit, leg, "--spec", spec],
                        capture_output=True, text=True)
    assert r2.returncode == 2, f"auditor must exit 2 on legacy drop; got {r2.returncode}\n{r2.stdout}"
    assert "SILENT DROP" in r2.stdout and "of 54" in r2.stdout
    print("[PASS] auditor catches both tally-lie and legacy silent drops (exit 2)")


if __name__ == "__main__":
    test_box1l()
    test_sunrise2l()
    test_incremental_flush()
    test_k4_incomplete_is_flagged()
    test_auditor_catches_silent_drop()
    print("\nALL TESTS PASSED")
