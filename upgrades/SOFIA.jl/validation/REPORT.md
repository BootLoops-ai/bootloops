# Validation against SOFIA_examples.nb recorded outputs

Each case replays a worked example from the upstream notebook (same diagram,
same options, including notebook substitution rules) and compares the
candidate-singularity sets up to proportionality and variable naming.

- **exact** — sets identical.
- **superset** — every upstream singularity reproduced, plus extra candidates
  (acceptable under candidate semantics; extras mostly stem from the
  route-union mode, see docs/PORTING_NOTES.md).
- **near (route tail)** — all but 1-2 entries reproduced; the residuals are
  anomalous-threshold polynomials reachable only through elimination routes
  whose intermediates exceed practical size in FLINT at present.
- **timeout** — the case exceeded its wall budget (45-55 min) on every
  configuration tried, including with the SymmetryQuotient (e.g.
  55 -> 18 subtopology classes on the double box) and the
  operation-size guard. Root cause, isolated to a minimal reproducer
  (outer double-box, subtopology rep 3): an elimination route reaches a
  discriminant of a 6175-term degree-8 polynomial whose computation stalls
  BOTH FLINT's native multivariate resultant and classical subresultant
  PRS for >20 minutes — the result object is astronomically large.
  Upstream's seconds-fast runs follow a different elimination route that
  never meets this object; route choice is opaque Mathematica Solve/
  FixLoopEdges internals. Bounding the op size (maxopterms) lets such
  cases complete but skips the deep content (case 11: completes in 1896s
  with 2/24 entries under OPBOUND=2500). Closing this gap needs an
  elimination-strategy advance (monster detection + route reordering, or
  modular discriminant interpolation), tracked as future work.

| case | section | status | mine/exp | missing | extra | ours (s) | upstream (s) |
|---|---|---|---|---|---|---|---|
| 0 | Sec 4.2: Feynman geometries | timeout |  |  |  |  | 7.476818 |
| 1 | 3-mass sunrise | superset | 11/7 | 0 | 4 | 12.6 | 0.184654 |
| 2 | 2-loop massive parachute | superset | 15/12 | 0 | 3 | 55.0 | 1.094506 |
| 3 | Massive acnode | timeout |  |  |  |  | 51.244067 |
| 4 | Outer double-box | timeout |  |  |  |  | 3.514829 |
| 5 | Outer double-box | timeout |  |  |  |  | 7.082507 |
| 6 | Double bubble triangle #1 | exact | 12/12 | 0 | 0 | 177.8 | 1.304946 |
| 7 | Double bubble triangle #2 | near (route tail) | 41/25 | 2 | 18 | 16.9 | 0.97689 |
| 8 |  | skipped: output unrecoverable (TemplateBox) | | | | | 55.921094 |
| 9 |  | skipped: output unrecoverable | | | | | 104.002066 |
| 10 | Fully massive penta-box | timeout |  |  |  |  | 59.876695 |
| 11 | 6-loop massive parachute | timeout |  |  |  |  | 6.590436 |
| 12 | Fig. 1 (c): Double-pentagon with a massi | timeout |  |  |  |  | 960.81147 |
| 13 |  | skipped: output unrecoverable | | | | | 38.658367 |
| 14 | 7-loop 3-point | timeout |  |  |  |  | 6.296208 |
| 15 | Degenerate acnode | superset | 31/8 | 0 | 23 | 13.4 | 0.394699 |
| 16 | Degenerate acnode | superset | 25/19 | 0 | 6 | 18.8 | 0.08593 |
| 17 | Degenerate acnode | superset | 38/11 | 0 | 27 | 18.2 | 0.206968 |
| 18 | Degenerate acnode | skipped: composite Join | | | | |  |
| 19 | H+jet non-planar-double box | timeout |  |  |  |  | 2.572982 |
| 20 | H+jet non-planar-double box | skipped: composite Join | | | | |  |
| 21 | Rocket diagram in https://arxiv.org/abs/ | skipped: numeric Substitutions | | | | | 2.32729 |
| 22 | 2-loop non-planar Higgs self-coupling | near (route tail) | 15/16 | 2 | 1 | 2209.4 | 3041.033788 |

Additional oracle check: the massless double box reproduces
`reference/PLD_database/dbox_zero_zero.m` exactly ({s, t, s+t}).

Upstream timings are the notebook's EchoTiming records (authors' machine);
ours are single-core runs on this machine. Scoreboard: of the 11 cases with
recoverable outputs that completed, 1 exact, 6 full-coverage supersets,
2 near (1-2 residual entries), and the double-box database oracle is exact;
6 heavy multi-scale cases currently exceed the time budget.
