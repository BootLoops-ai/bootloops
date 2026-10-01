#!/usr/bin/env python3
"""fold_hnf validation battery (see ../README.md).

Reference data (not included in the package; directory supplied with
TERRIER_G4_DIR): the recorded G4 corner-edge list (256 edges, jsonl) and
the reference per-B raw/sl2 counts (json).

T1 HEXAGONAL-TIE CORNERS: all 256 recorded corner edges merge under
   fold_hnf; the 416 recorded forms fold to 160 classes (fibers 64x2 +
   96x3); recorded HNF labels, pivots, det-signs AND recorded SL2 words
   re-verified on the referee witness.py.
T2 ACCEPTANCE (anchor-set): G4 raw box (Linf<=2, budgets 1..16) folded by
   fold_hnf must give B=1 bin = 6,612 and full-window total = 49,544
   EXACTLY. Raw per-B histogram must reproduce the reference raw counts;
   per-B class counts must be <= the reference N_sl2_per_B (refinement
   direction). Any mismatch is RECORDED (measured AND expected printed)
   and FAILS — never forced.
T3 NO OVER-MERGE: 30 fixed-seed merged raw pairs get explicit SL2
   witness words verified on witness.py (budget asserted per step), plus
   60 rep-membership witnesses and 5 tampered-word must-fail controls.
T4 DETERMINISM + INPUT-ORDER INVARIANCE (fixed shuffle); idempotence
   + key stability; reduce() call-shape words verified; refinement
   spot-check vs engine_a.sl2_reduce (public API, read-only).

Envelope: single core; the battery asserts < 25 CPU-min and finishes far
inside. Fixed SEED = 20260716. Exact integers only.
"""
import json
import os
import random
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
MOD = os.path.dirname(HERE)                       # census/fold_hnf
PKG = os.path.dirname(os.path.dirname(MOD))       # package root
sys.path.insert(0, MOD)
sys.path.insert(0, os.path.join(PKG, "census", "dedup"))

import fold_hnf            # noqa: E402  module under test
import witness             # noqa: E402  referee (read-only API)
import engine_a            # noqa: E402  engine A (read-only API)

SEED = 20260716
# G4 anchor-set reference directory comes from the environment (no default);
# unset => loud refusal (fail-closed).
G4D = os.environ.get("TERRIER_G4_DIR", "")
if not G4D:
    raise SystemExit("REFUSE: env TERRIER_G4_DIR unset — must point at "
                     "the G4 anchor-set reference directory (not included "
                     "in the package)")
EXPECTED_B1, EXPECTED_TOTAL = 6612, 49544
T0 = time.process_time()


def cpu():
    return time.process_time() - T0


def load_pairs():
    with open(os.path.join(G4D, "FOLD_GAP_PAIRS.jsonl")) as fh:
        return [json.loads(line) for line in fh if line.strip()]


def t1():
    recs = load_pairs()
    assert len(recs) == 256, "expected 256 edges, got %d" % len(recs)
    forms = set()
    for r in recs:
        a = (tuple(r["v1_rep_a"][0]), tuple(r["v1_rep_a"][1]))
        b = (tuple(r["v1_rep_b"][0]), tuple(r["v1_rep_b"][1]))
        forms.add(a)
        forms.add(b)
        ka, kb = fold_hnf.key(a), fold_hnf.key(b)
        assert ka == kb, "edge %d endpoints do NOT merge" % r["pair_id"]
        (c1, c2, i, j), det = ka
        assert list(c1) == r["hnf_H"][0] and list(c2) == r["hnf_H"][1], \
            "HNF mismatch vs recorded label, edge %d" % r["pair_id"]
        assert [i, j] == r["hnf_pivots"] and det == r["det_sign"], \
            "pivot/det mismatch vs recorded label, edge %d" % r["pair_id"]
        assert witness.verify(a, tuple(r["sl2_word_a_to_b"]), b, []), \
            "recorded SL2 word failed re-verification, edge %d" % r["pair_id"]
    assert len(forms) == 416, "expected 416 forms, got %d" % len(forms)
    out = fold_hnf.fold(sorted(forms))
    assert out["N_classes"] == 160, \
        "expected 160 classes, got %d" % out["N_classes"]
    hist = {}
    for m in out["members"]:
        hist[len(m)] = hist.get(len(m), 0) + 1
    assert hist == {2: 64, 3: 96}, "fiber histogram %r != {2:64, 3:96}" % hist
    print("T1 PASS: 256/256 edges merge; 416 forms -> 160 classes; fibers "
          "{2:64, 3:96}; recorded labels+words re-verified (%.1f cpu-s)"
          % cpu(), flush=True)


def t2():
    with open(os.path.join(G4D, "COUNTS_G4.json")) as fh:
        c4 = json.load(fh)
    C, BMAX = 2, 16
    rng = range(-C, C + 1)
    vecs = [(a, b, c, d) for a in rng for b in rng for c in rng for d in rng]
    per_b = {b: set() for b in range(1, BMAX + 1)}
    bank1 = {}                    # B=1: key -> first (up to) two raw states
    n_raw = {}
    for h in vecs:
        w0, w1, w2, w3 = h[2], h[3], -h[0], -h[1]
        for f in vecs:
            b = f[0] * w0 + f[1] * w1 + f[2] * w2 + f[3] * w3
            if b < 0:
                b = -b
            if b == 0 or b > BMAX:
                continue
            n_raw[b] = n_raw.get(b, 0) + 1
            k = fold_hnf.key((f, h))
            per_b[b].add(k)
            if b == 1:
                e = bank1.get(k)
                if e is None:
                    bank1[k] = [(f, h)]
                elif len(e) == 1:
                    e.append((f, h))
    got_raw = {str(b): n_raw.get(b, 0) for b in range(1, BMAX + 1)}
    assert got_raw == c4["N_raw_box_per_B"], \
        "raw scan does not reproduce the reference raw counts: %r" % got_raw
    meas = {b: len(per_b[b]) for b in range(1, BMAX + 1)}
    total = sum(meas.values())
    v1 = {int(k): v for k, v in c4["N_sl2_per_B"].items()}
    print("T2 N_hnf per B 1..16:", [meas[b] for b in range(1, 17)], flush=True)
    print("T2 ref N_sl2 per B  :", [v1[b] for b in range(1, 17)], flush=True)
    print("T2 MEASURED B=1 bin = %d (anchor-set %d); full-window total = %d "
          "(anchor-set %d)" % (meas[1], EXPECTED_B1, total, EXPECTED_TOTAL),
          flush=True)
    for b in range(1, BMAX + 1):
        assert meas[b] <= v1[b], \
            "refinement direction broken at B=%d: %d > %d" % (b, meas[b], v1[b])
    assert meas[1] == EXPECTED_B1 and total == EXPECTED_TOTAL, (
        "ACCEPTANCE RECORDED-NOT-FORCED: measured (B1=%d, total=%d) vs "
        "anchor-set (B1=%d, total=%d) — both recorded above, NOT forced"
        % (meas[1], total, EXPECTED_B1, EXPECTED_TOTAL))
    print("T2 PASS: ACCEPTANCE B=1 bin = %d and full-window total = %d "
          "EXACTLY; raw histogram == reference; N_hnf <= ref N_sl2 "
          "per B (%.1f cpu-s)" % (meas[1], total, cpu()), flush=True)
    return bank1


def t3(bank1):
    rnd = random.Random(SEED)
    elig = sorted(k for k, v in bank1.items() if len(v) == 2)
    assert len(elig) >= 30, "only %d multi-member classes" % len(elig)
    picks = rnd.sample(range(len(elig)), 30)
    for ix in picks:
        a, b = bank1[elig[ix]]
        ra, rb = fold_hnf.rep_of(a), fold_hnf.rep_of(b)
        assert ra == rb, "merged pair has different reps"
        wa, wb = fold_hnf.witness_word(a), fold_hnf.witness_word(b)
        assert witness.verify(a, wa, ra, []), "rep-membership witness (a)"
        assert witness.verify(b, wb, rb, []), "rep-membership witness (b)"
        w = wa + fold_hnf.inv_word(wb)
        assert witness.verify(a, w, b, []), "explicit pair witness failed"
    nbad = 0
    for ix in picks[:5]:
        a, b = bank1[elig[ix]]
        w = (fold_hnf.witness_word(a)
             + fold_hnf.inv_word(fold_hnf.witness_word(b)))
        assert not witness.verify(a, w + ("T",), b, []), \
            "tampered word ACCEPTED — referee broken"
        nbad += 1
    print("T3 PASS: 30/30 explicit SL2 pair witnesses + 60 rep-membership "
          "witnesses verified on witness.py; %d/5 tampered words rejected "
          "(%.1f cpu-s)" % (nbad, cpu()), flush=True)


def t4():
    # corpus: Linf<=1 box (all in-window budgets) + the 416 corner forms
    rng = range(-1, 2)
    vecs = [(a, b, c, d) for a in rng for b in rng for c in rng for d in rng]
    corpus = []
    for h in vecs:
        w0, w1, w2, w3 = h[2], h[3], -h[0], -h[1]
        for f in vecs:
            b = f[0] * w0 + f[1] * w1 + f[2] * w2 + f[3] * w3
            if b:
                corpus.append((f, h))
    for r in load_pairs():
        corpus.append((tuple(r["v1_rep_a"][0]), tuple(r["v1_rep_a"][1])))
        corpus.append((tuple(r["v1_rep_b"][0]), tuple(r["v1_rep_b"][1])))
    o1 = fold_hnf.fold(corpus)
    o2 = fold_hnf.fold(corpus)
    assert o1 == o2, "determinism broken (same input, different output)"
    perm = list(range(len(corpus)))
    random.Random(SEED).shuffle(perm)
    o3 = fold_hnf.fold([corpus[p] for p in perm])
    assert o3["class_keys"] == o1["class_keys"], "keys depend on input order"
    assert o3["reps"] == o1["reps"], "reps depend on input order"
    inv = [0] * len(perm)
    for newpos, old in enumerate(perm):
        inv[old] = newpos
    m1 = [sorted(inv[i] for i in m) for m in o1["members"]]
    m3 = [sorted(m) for m in o3["members"]]
    assert m1 == m3, "member partition not permutation-consistent"
    rnd = random.Random(SEED + 1)
    nwit = 0
    for it in range(200):
        st = corpus[rnd.randrange(len(corpus))]
        r = fold_hnf.rep_of(st)
        assert fold_hnf.key(r) == fold_hnf.key(st), "key not rep-stable"
        assert fold_hnf.rep_of(r) == r, "rep_of not idempotent"
        if it < 40:
            rep, toks = fold_hnf.reduce(st)
            assert rep == r and witness.verify(st, toks, rep, []), \
                "reduce() word failed on referee"
            nwit += 1
    chk = 0
    for st in corpus[::53]:       # fixed spot-check stride
        red, _toks = engine_a.sl2_reduce(st)
        assert fold_hnf.key(red) == fold_hnf.key(st), \
            "engine_a.sl2_reduce left the HNF class (refinement broken)"
        chk += 1
    print("T4 PASS: determinism + fixed-shuffle order invariance on %d "
          "states; 200 idempotence/key-stability; %d reduce() words verified; "
          "%d refinement checks vs engine_a.sl2_reduce (%.1f cpu-s)"
          % (len(corpus), nwit, chk, cpu()), flush=True)


def main():
    print("fold_hnf battery start (fixed SEED %d)" % SEED, flush=True)
    t1()
    bank1 = t2()
    t3(bank1)
    t4()
    assert cpu() < 1500, "approaching the 25-CPU-min envelope"
    print("BATTERY GREEN: T1 T2 T3 T4 all PASS in %.1f cpu-s" % cpu(),
          flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
