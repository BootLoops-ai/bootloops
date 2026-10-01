# dedup — two independent flux-orbit deduplication engines (A + B), witness referee, reconciler

Two INDEPENDENT flux-orbit deduplication engines for pairs (f,h) under
monodromy x SL2(Z), a third engine-independent witness checker
(`witness.py`) that verifies every merge, and a reconciler
(`reconcile.py`) that adjudicates the two partitions.

## Group action (both engines implement this spec verbatim)
- monodromy generator M in Sp(2k,Z) (validated exactly): (f,h) -> (Mf,Mh)
- SL2(Z) doublet: S:(f,h)->(-h,f), T:(f,h)->(f+h,h), + exact inverses
- budget |f^T Sigma h| is an exact invariant (asserted per step by the
  referee); exact integer arithmetic everywhere.
- FACTORIZATION both engines exploit: monodromy (left mult on f,h) and
  SL2(Z) (column ops on [f h]) COMMUTE, so the SL2 factor is reduced
  EXACTLY (Lagrange-Gauss, token-emitting) and the bounded search runs
  over monodromy words only (4-ary, depth-capped) — no composite-ball
  blowup (the naive composite ball is intractable).

## Two-engine law
| | Engine A (engine_a.py) | Engine B (engine_b.py) |
|---|---|---|
| traversal | BFS queue | depth-capped DFS stack |
| visited container | hash dict | sorted list + bisect |
| SL2 reduction tie law | nearest, ties half-UP | floor+fix, ties DOWN |
| prune | L2^2 cap + depth | L1 cap + depth |
| canonical key | (L2^2, lex f‖h) | (L1, Linf, lex h‖f) |
No shared code. Both emit generator-word WITNESSES (start -> canonical
rep); witness.py is a third, referee implementation that checks witnesses
exactly. Bounded search can only UNDER-merge, never over-merge: every
merge carries a verifiable witness.

## Disagreement path (reconcile.py)
1. Partitions equal -> GREEN-AGREE (two-stack control GREEN).
2. Unequal, all extra-merge witnesses verify -> RED-RECONCILED-MISS: the
   engine missing a witnessed merge under-merged; reconciled partition
   (closure of witnessed merges) is reported for the audit, but a
   two-stack count mismatch is a RED control: halt=True, reported, never
   silently consumed.
3. Any witness fails -> RED-WITNESS-FAIL (integrity halt).
4. Plant sharing a class with an honest candidate -> RED-PLANT-MERGED
   (must-fail law). Plants are inserted PRE-dedup and must survive
   to the distance gate as their own orbit class.
RAW and DEDUPED counts are first-class in every report; the orbit-size
distribution is reported split on/off locus (diagnostic column).

## Tests — synthetic battery (no census data)
`python3 test_dedup.py` (single core, minutes). Fixed seeds:
TEST_SEEDS.json. Ground truth is proof-level: families seeded at pairwise
distinct budgets (exact invariant => cross-family merges impossible);
membership witnessed by construction; generation paths capped inside the
engines' search balls. T1-T8 cover generator validation, invariance,
both engines vs ground truth, witness audit, GREEN reconciliation with
plant survival + histogram columns, the RED-RECONCILED-MISS drill
(crippled engine), and the RED-WITNESS-FAIL tamper drill.
Receipts: test_receipts.json (regenerated on every run). The synthetic
battery exercises every engine path (both engines, witnesses,
reconciliation, both RED drills).
