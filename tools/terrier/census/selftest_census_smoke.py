#!/usr/bin/env python3
"""Battery CEN-X (cross-wing smoke): tiny END-TO-END census on a SYNTHETIC
2x2 example through the census wing, receipted via common/receipts.py.
Pipeline: raw 2x2 window -> fold_hnf (complete SL2 orbit key; words
witness-verified; key invariance under S/T on every state = completeness
at window scale) -> dedup two-stack A+B + plant -> reconcile GREEN-AGREE
-> dist_cert exact ON/OFF -> common receipt (sha-stamped, validated).
Writes census/SMOKE_RECEIPT.json (volatile, unpinned)."""
import os, random, sys
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(HERE, "fold_hnf"))
sys.path.insert(0, os.path.join(HERE, "dedup"))
sys.path.insert(0, os.path.join(HERE, "dist_cert"))
sys.path.insert(0, os.path.join(ROOT, "common"))
import fold_hnf, witness, engine_a, engine_b, reconcile, gens_synth
import lattice_exact as lx
import receipts as rc

# 1. synthetic raw window: every 2x2 integer doublet with entries in
#    {-1,0,1} and nonzero det (columns c1, c2)
from itertools import product
W = [((a, c), (b, d)) for a, b, c, d in product((-1, 0, 1), repeat=4)
     if a * d - b * c != 0]
F = fold_hnf.fold(W)
assert F["N_in"] == len(W) and F["N_classes"] < len(W)
for st in W:                       # soundness: witness word st -> rep
    rep, word = fold_hnf.reduce(st)
    assert witness.verify(st, word, rep, []), f"word fail {st}"
for st in W:                       # completeness cert: key is S/T-invariant
    for tok in ("S", "T", "s", "t"):
        img = witness.step(st, tok, witness.make_table([], 2))
        assert fold_hnf.key(img) == fold_hnf.key(st), f"key moves {st} {tok}"
F2 = fold_hnf.fold(F["reps"])      # idempotence on canonical reps
assert F2["N_classes"] == F["N_classes"] == len(set(F["keys"]))
rng = random.Random(20260717)      # committed-shuffle order invariance
Wsh = list(W)
rng.shuffle(Wsh)
assert fold_hnf.fold(Wsh)["class_keys"] == F["class_keys"]
print(f"[smoke fold] {F['N_in']} raw -> {F['N_classes']} classes; words + "
      f"S/T key-invariance + idempotence + order-invariance PASS")

# 2. dedup two-stack: folded reps + 2 monodromy-transported copies (must
#    merge back) + one far plant (unique budget, must stay own-class)
gens = gens_synth.make_gens(20260717, n=2, count=2)
tab = witness.make_table(gens, 2)
moved = [witness.apply(F["reps"][k], w, gens, tab)
         for k, w in ((0, ("M0", "M1")), (1, ("m1", "M0", "M0")))]
states = list(F["reps"]) + moved + [((97, 2), (3, 89))]
tags = ["honest"] * (len(states) - 1) + ["plant"]
labels = ["rep"] * len(F["reps"]) + ["moved"] * 2 + ["plant"]
resA = engine_a.dedup(states, gens)
resB = engine_b.dedup_b(states, gens)
rep = reconcile.reconcile(states, gens, resA, resB, tags, labels)
assert rep["status"] == "GREEN-AGREE" and not rep["halt"], rep["status"]
assert rep["plant_own_class"], "plant merged into an honest class"
assert rep["N_dedup_A"] == rep["N_dedup_B"] == len(states) - 2, \
    "transported copies did not merge back (or over-merge)"
print(f"[smoke dedup] {len(states)} states -> A=B={rep['N_dedup_A']} classes "
      f"GREEN-AGREE, 2 transported merged back, plant own-class PASS")

# 3. dist_cert exact ON/OFF on a 2x2 parent (A1+A1, gram diag(2,2))
G2x2 = lx.direct_sum(lx.gram_An(1), lx.gram_An(1))
d_on, wit = lx.dist_cert_exact(G2x2, [[0, 1]])
assert d_on == 0 and wit is not None, "expected ON-locus with root witness"
d_off, cert = lx.dist_cert_exact(G2x2, [[1, 1]])
assert d_off == "OFF" and cert["roots_examined"] == 0, "expected exact OFF"
print("[smoke dist_cert] flux [0,1] ON (root witness) + flux [1,1] "
      "certified OFF PASS")

# 4. cross-wing receipt via common/receipts.py (sha-stamped, validated)
pins = {"code_sha256": rc.code_sha(
    [os.path.join(HERE, "fold_hnf", "fold_hnf.py"),
     os.path.join(HERE, "dedup", "reconcile.py"), __file__])}
rec = rc.make_receipt(
    "census-smoke-2x2/v1", "COMPLETE", pins,
    window=len(W), fold_classes=F["N_classes"],
    dedup_classes=rep["N_dedup_A"], dedup_status=rep["status"],
    dist_cert={"on": str(d_on), "off": str(d_off)})
out = os.path.join(HERE, "SMOKE_RECEIPT.json")
rc.write_receipt(out, rec)
rc.validate_receipt(__import__("json").load(open(out)))
print(f"[smoke receipt] {out} written + re-validated PASS")
print("SELFTEST census_smoke: ALL PASS (fold -> dedup -> reconcile -> "
      "dist_cert -> receipt)")
