#!/usr/bin/env python3
"""Unit gates for certify.py (lambda certificates for the B2FT eliminator).

Tiny synthetic stratified system (seconds-class, no external paths):
  C1 certificates: B2FT stratified solve WITH transcript + close_all;
     extract_lambdas for every pivot; INDEPENDENT check_certificate passes
     for all; transcripted in-place closure == coverage.assemble_closed.
  C2 default-off: stratified_solve/eliminate_fast WITHOUT the transcript
     kwarg produce byte-identical subs to the transcripted run (the hook
     provably changes nothing by default).
  C3 planted faults on scratch copies: lam flip -> FAIL, c corrupt -> FAIL,
     corrupted transcript AXPY factor -> re-extracted lam FAILs on affected
     rows, PASSes elsewhere.
  C4 serialization: save/load round-trip preserves lam + c exactly.
  C5 standalone eliminate_fast transcript (caller-registered inputs).

Usage: pytest (or python3 -m pytest tests/test_certify.py from tool/strata).
"""
import os
import random
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))

import certify                                          # noqa: E402
from coverage import assemble_closed                    # noqa: E402
from fp_eliminate import eliminate_fast, stratified_solve  # noqa: E402

P = 1048573
MASTERS = {0, 1, 2}
N_COLS = 60


def synth_system(seed=5):
    """Stratified synthetic: 3 strata (col//20), masters {0,1,2} forbidden,
    underdetermined so closed rows carry nonzero master tails."""
    rng = random.Random(seed)
    order = {c: c for c in range(N_COLS) if c not in MASTERS}
    stratum_of = {c: c // 20 for c in range(N_COLS)}
    rows = []
    for t in (2, 1, 0):
        for _ in range(12):
            lead = rng.randrange(max(20 * t, 3), 20 * t + 20)
            r = {lead: rng.randrange(1, P)}
            for _ in range(4):
                c = rng.randrange(0, 20 * t + 20)
                r[c] = rng.randrange(1, P)
            rows.append(r)
    return rows, order, stratum_of


def solve_t(rows, order, stratum_of):
    tr = certify.Transcript(P)
    subs, left, _ = stratified_solve([dict(r) for r in rows], order, P,
                                     stratum_of, policy="B2FT",
                                     forbid=MASTERS, transcript=tr)
    assert not left
    certify.close_all(subs, order, P, transcript=tr)
    return subs, tr


def test_certificates_all_pass():
    rows, order, stratum_of = synth_system()
    subs, tr = solve_t(rows, order, stratum_of)
    # transcripted in-place closure == assemble_closed of a clean run
    subs0, left0, _ = stratified_solve([dict(r) for r in rows], order, P,
                                       stratum_of, policy="B2FT",
                                       forbid=MASTERS)
    assert not left0
    assert subs == assemble_closed(subs0, order, P)
    n_tail = sum(1 for r in subs.values() if any(c in MASTERS for c in r))
    assert n_tail >= 5, "degenerate synthetic (no master tails)"
    lams, _ = certify.extract_lambdas(tr, {w: subs[w] for w in subs},
                                      len(rows), P, chunk=7)
    for w in subs:
        c_dict = {c: (-v) % P for c, v in subs[w].items()}
        assert certify.check_certificate(rows, lams[w], w, c_dict, P), \
            f"certificate FAIL at pivot {w}"


def test_default_off_unchanged():
    rows, order, stratum_of = synth_system()
    subs_t, _tr = solve_t(rows, order, stratum_of)
    # no-kwarg call path (proves default-off identity)
    subs_a, left_a, _ = stratified_solve([dict(r) for r in rows], order, P,
                                         stratum_of, policy="B2FT",
                                         forbid=MASTERS)
    assert not left_a
    assert assemble_closed(subs_a, order, P) == subs_t
    # eliminate_fast standalone, both ways
    s1, l1, _ = eliminate_fast([dict(r) for r in rows], order, P,
                               forbid=MASTERS)
    tr = certify.Transcript(P)
    rows2 = [dict(r) for r in rows]
    for i, r in enumerate(rows2):
        tr.reg_new(r, i)
    s2, l2, _ = eliminate_fast(rows2, order, P, forbid=MASTERS,
                               transcript=tr)
    assert s1 == s2 and l1 == l2


def test_planted_faults():
    rows, order, stratum_of = synth_system()
    subs, tr = solve_t(rows, order, stratum_of)
    lams, _ = certify.extract_lambdas(tr, {w: subs[w] for w in subs},
                                      len(rows), P)
    rng = random.Random(0)
    tails = [w for w in subs if any(c in MASTERS for c in subs[w])]
    for w in tails[:5]:
        c_dict = {c: (-v) % P for c, v in subs[w].items()}
        ok, _ = certify.planted_fault(rows, lams[w], w, c_dict, P,
                                      "lam_flip", rng)
        assert not ok, "lam flip NOT caught"
        ok, _ = certify.planted_fault(rows, lams[w], w, c_dict, P,
                                      "c_corrupt", rng)
        assert not ok, "c corruption NOT caught"
    # transcript factor corruption: brute-force AXPY ops until one changes
    # some lam (tiny system); affected must FAIL, unaffected must PASS
    kinds = tr.finalize()[0]
    hit = None
    for op in [i for i, k in enumerate(kinds) if k == certify.OP_AXPY]:
        tr2 = certify.corrupt_transcript_factor(tr, op)
        lams2, _ = certify.extract_lambdas(tr2, {w: subs[w] for w in subs},
                                           len(rows), P)
        affected = [w for w in subs if lams2[w] != lams[w]]
        if affected:
            hit = (op, lams2, affected)
            break
    assert hit, "no AXPY op affects any lambda"
    _op, lams2, affected = hit
    for w in subs:
        c_dict = {c: (-v) % P for c, v in subs[w].items()}
        ok = certify.check_certificate(rows, lams2[w], w, c_dict, P)
        if w in set(affected):
            assert not ok, f"corrupted-transcript cert PASSED at {w}"
        else:
            assert ok, f"unaffected cert FAILED at {w}"


def test_serialization_roundtrip():
    rows, order, stratum_of = synth_system()
    subs, tr = solve_t(rows, order, stratum_of)
    ws = sorted(subs)[:6]
    lams, _ = certify.extract_lambdas(tr, {w: subs[w] for w in ws},
                                      len(rows), P)
    certs = [{"weight": w,
              "c": {c: (-v) % P for c, v in subs[w].items()},
              "lam": lams[w], "lam_nnz": len(lams[w]),
              "check_wall_s": 0.0} for w in ws]
    with tempfile.TemporaryDirectory() as td:
        certify.save_certificates(td, "toy", [P, 3, 7], certs)
        back = certify.load_certificates(td, "toy")
    assert len(back) == len(certs)
    for a, b in zip(certs, back):
        assert a["weight"] == b["weight"]
        assert a["lam"] == b["lam"]
        assert {int(k): v for k, v in a["c"].items()} == b["c"]
        assert certify.check_certificate(rows, b["lam"], b["weight"],
                                         b["c"], P)


def test_standalone_eliminate_fast_transcript():
    rows, order, stratum_of = synth_system(seed=9)
    tr = certify.Transcript(P)
    rows2 = [dict(r) for r in rows]
    for i, r in enumerate(rows2):
        tr.reg_new(r, i)
    subs, left, _ = eliminate_fast(rows2, order, P, forbid=MASTERS,
                                   transcript=tr)
    assert not left
    lams, _ = certify.extract_lambdas(tr, {w: subs[w] for w in subs},
                                      len(rows), P)
    for w in subs:
        c_dict = {c: (-v) % P for c, v in subs[w].items()}
        assert certify.check_certificate(rows, lams[w], w, c_dict, P)
