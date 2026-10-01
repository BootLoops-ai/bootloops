# fold_hnf — complete SL2(Z)-orbit fold for flux doublets (HNF + det-sign key)

SCOPE: the SL2(Z) canonicalizer used to fold a raw flux window into orbit
classes before any counting. It REPLACES `engine_a.sl2_reduce` in that role
(the dedup engines themselves are unchanged; engine_a/b monodromy dedup
still uses its own reduction internally, where an incomplete SL2 fold can
only under-merge).

## Why not engine_a.sl2_reduce
engine_a.sl2_reduce's post-descent candidate closure omits the tie-shear
identification on the HEXAGONAL CORNER of the Gauss domain
(|f|^2 = |h|^2 = 2|f.h|), so each corner orbit splits into 2-3
descent-stable reduced forms: on the G4 anchor set at B=1, 416 sl2_reduce
classes collapse to 160 correct classes (fibers 64x2 + 96x3; excess
exactly 256; all edges witness-verified). That split is in the
over-count-safe direction; this module removes it.

## Invariant (complete)
For st = (f, h) with A = [f h] of rank 2 (guaranteed at budget >= 1):
H = column Hermite normal form of A under the right GL2(Z) action,
A*U = H with U UNIQUE (A injective). key(st) = (H, det U) is a COMPLETE
SL2(Z)-orbit invariant — same key <=> same orbit; the det-sign separates
the two SL2 classes inside each GL2 class. Canonical representative:
rep_of(st) = A*V with V = U (det +1) or U*J, J = diag(1,-1) (det -1) —
an orbit MEMBER, identical for every member, idempotent. witness_word(st)
emits an S/s/T/t token word st -> rep, checkable on the referee
`dedup/witness.py` (budget asserted per step). Exact integers only; no imports.

## Public API (fold_hnf.py)
- key(st) -> ((c1, c2, i, j), det_sign) — complete orbit key
- rep_of(st) -> canonical orbit member
- witness_word(st) / inv_word(word) — referee-checkable token words
- reduce(st) -> (rep, toks) — engine_a.sl2_reduce call shape
- fold(states) -> {N_in, N_classes, keys, class_keys, reps, members}
  (deterministic; class_keys/reps input-ORDER invariant)

## Battery (tests/test_fold_hnf.py; fixed SEED 20260716)
Needs the G4 anchor-set reference data (the 256 recorded corner edges and
the reference per-B counts; not included in the package — supply the
directory with `TERRIER_G4_DIR`).
- T1 hexagonal-tie corners: all 256 recorded corner edges merge;
  416 forms -> 160 classes (fibers 64x2 + 96x3); recorded labels + SL2
  words re-verified on witness.py.
- T2 ACCEPTANCE (anchor-set): G4 raw box (Linf<=2, B=1..16) fold gives
  B=1 bin = 6,612 and full-window total = 49,544 EXACTLY; raw per-B
  histogram reproduces the reference counts; N_hnf(B) <= N_sl2(B) (the
  sl2_reduce class count) per B. Mismatches are RECORDED (measured +
  expected printed), never forced.
- T3 no over-merge: 30 fixed-seed merged pairs carry explicit SL2
  witness words verified on witness.py (+60 rep-membership witnesses,
  5 tampered-word must-fail controls).
- T4 determinism + fixed-shuffle input-order invariance; idempotence
  + key stability; reduce() words verified; refinement spot-check vs
  engine_a.sl2_reduce (public API, read-only).
Run: python3 tests/test_fold_hnf.py (single core; asserts < 25 CPU-min,
finishes far inside).

## Validation of the invariant
Independently re-derived and checked: 400/400 random SL2 words leave the
key invariant; the J-twist det-sign law; 256/256 corner merges witnessed
by explicit SL2 words; two independent folds (a closure-canonical fold and
a separate HNF fold) agree on the B=1 anchor count 6,612.
