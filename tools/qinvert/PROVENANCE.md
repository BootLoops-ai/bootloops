# qinvert PROVENANCE — every convention → its doc cite or receipt

Rule of the kit: no convention without a source.  "SAS doc" cites name the
Base SAS 9.4 Procedures Guide section; "origin" cites name a file of the
originating analysis (not included in this package; read-only, where each
convention was validated); "receipt" cites name an executable artifact
(planted truth or negative control in `tests/test_all.py`, or a run record of
the originating analysis, not included).

## targets.py

| convention | source |
|---|---|
| SAS QNTLDEF=1..5 formulas (np=j+g; def-1 edge convention x_0:=x_1, x_{n+1}:=x_n; def-2 half-case even/odd rule; def-4 h=(n+1)p; def-5 EDF-with-averaging) | Base SAS 9.4 Procedures Guide: Statistical Procedures, "The UNIVARIATE Procedure", *Calculating Percentiles* (QNTLDEF=/PCTLDEF=; PROC MEANS shares the definitions, default 5).  Receipts: `test_sas_hand_values` |
| R types 1–9 (discontinuous 1–3 step rules incl. type-3 ties-to-even; continuous 4–9 h formulas with clamping) | Hyndman & Fan (1996), *Sample Quantiles in Statistical Packages*, Am. Stat. 50(4) 361–365; R `?quantile`.  External receipts: R's documented `quantile(1:10, .25, type=i)` outputs, `test_r_reference_values_1_to_10` |
| Equivalences SAS1=R4, SAS2=R3, SAS3=R1, SAS4=R6, SAS5=R2 | H&F sec. 2 + SAS doc notes; asserted executable: `test_documented_equivalences`; anti-aliasing negative control `test_definitions_not_all_aliased` |
| QNTLDEF=5 end-to-end | origin `src/tukey.py` (CMS Star Ratings Technical Notes 2024 final p.149–150), validated vs published K-5/K-6 fence tables; read-only cross-check `test_origin_qntldef5_crosscheck` |
| `value()` DERIVED from `support()` for every quantile (evaluator and solver sensitivity structure cannot disagree) | kit design rule; pinned externally by the evaluator receipts above |
| Tukey fence algebra lo=(1+m)Q1−mQ3, hi=(1+m)Q3−mQ1; inversion Q1=((1+m)lo+m·hi)/(1+2m), Q3=(m·lo+(1+m)hi)/(1+2m); capped-side ⇒ inequality-only (modes both/lo_only/hi_only/none) | origin `src/inverse_membership.py` + `inverse_membership2.required_qs`, planted-truth batteries there; receipts `TestFenceAlgebra` |
| Fence cap convention (lo maxed with cap, hi minned — CMS percent/rate caps) | origin `src/tukey.py`, CMS Tech Notes 2024 final p.149–150 |
| `as_fraction` refuses floats; decimal strings parsed exactly | kit exactness law (no float in any accept/reject decision; zero float-semantic sites in the replayed pipelines — module docstring) |
| Out-of-range OrderStat (index outside [1, n]) evaluates to (None, False) — not exact, never a crash | our convention (a ValueError escaping `solve()` would kill a `stack()` run); receipts `TestOrderStatOutOfRange` |

## solver.py

| convention | source |
|---|---|
| Index-shift lemma (delta of k moves any fixed value's rank by ≤ k), window pad k+4 | origin `src/membership_stacker.py` windowing; proof in README "The algebra" |
| Profile enumeration: blocks → realizations → left-to-right assembly with filler windows | origin `src/inverse_membership2.py` profile solver + `membership_stacker.solve_adds_k` / `_quartile_realizations` (planted-truth batteries + CMS K-table reconciliation receipts in that tree) |
| Parity/averaging constraints generalized per definition via `support(n')` (the origin's n mod 4 argument at QNTLDEF=5, p=1/4, 3/4) | origin `src/inverse_membership.py`; receipts `TestParity`, `TestSolverParity` |
| Equal-pinned-value SPAN MERGE: blocks pinned to the same v merge into one x'_(i)=v span over all covered indices; a wavg block in the run realizes only its both-equal branch | origin equal-quartile escape `src/membership_stacker.py` (`... and q1 != q3`) + order-statistic monotonicity.  Silent pruning of these scenarios (a strict chaining inequality) is the failure mode the receipts pin.  Planted-truth receipts `TestSameValueBlocks` (adjacent order stats k=2; equal sas3 quantiles on tie-heavy 2dp data k=2; Q1==Q3 k=1; non-adjacent equal blocks k=1; wavg+point merge k=3) + negative controls (non-monotone pins, unreachable span) |
| Realization copy-count e ranges 0..k (a span may need up to k fresh copies of a scarce v) | counting argument in `_real_span` docstring |
| Straddle-touching skip is LOUD (`touching-realizations-skipped`): a straddle end equal to a neighboring DIFFERENT pinned value is a structural limit, not handled by the chained count accounting | our convention; receipt `TestLoudness.test_touching_straddle_is_loud` (planted k=1 truth missed loudly, k=2 returned) |
| Out-of-range OrderStat at some n' = infeasible SIZE: skipped exactly (no cap note), k loop proceeds | our convention; receipts `TestOrderStatOutOfRange` (k found at the size that makes the index valid; rem_cap path; negative index) |
| `k_cap-exhausted` recorded on EVERY witness-free search (emptiness certified only up to the echoed k_cap/rem_cap) | our convention (k_cap/rem_cap exhaustion must never leave caps_hit empty while the README claims all caps loud); receipts `TestLoudness.test_k_cap_exhausted_note`, negative controls in `TestNegativeControls` |
| `register-none-grid-data-values-only` recorded when a register-None one-sided fence mode grids the free quantile | our convention; receipt `TestLoudness.test_register_none_fence_grid_is_loud` |
| 800-multiple grid cap per candidate window | origin `membership_stacker._q_pair_candidates` (same cap in the certified runs) |
| Verification gate: every witness materialized and re-evaluated via `targets.evaluate_targets`; failures dropped (soundness independent of enumeration shortcuts) | origin pattern (all origin witnesses exact-verified); enforced on every family in every solver test |
| Removal candidates complete up to region-equivalence for order-stat targets; dense (Mean) truncated loudly | argument in `_removal_candidates` docstring; receipts `TestPlantedRemoval`, duplicate-removal attack battery |
| spawn-only multiprocessing, module-level worker for spawn pickling | origin `src/solve_all_misses.py` (fork+BLAS deadlocks under load; spawn is mandatory); receipt `test_solve_many_spawn_matches_serial` |

## stacker.py

| convention | source |
|---|---|
| Cross-dataset stacking: shared entities, class-gated participation, coupling (sub multiset rides super) | origin `src/membership_stacker.py` (C/D complaints coupling, verified in-data there); receipts `TestStackerE2E` |
| Hard no-regression gate (StackerError rather than a regressing delta) | origin `removal_is_globally_safe` gate; receipts `test_removal_gate_blocks_shared_entity`, final battery in `stack()` |
| Greedy assembly certifies feasibility + no-regression, NOT entity-count minimality; `unsolved` = "not found within caps" | README "Honest limits"; stacker docstring |
| Regression anchor: star-year 2024 D02/MA-PD, published fences (0, 1.41), n=534, unique +1-addition family, witness 21/50 | originating analysis record `membership_families.json` (not included); double-gated through the origin tree's own `tukey.tukey_fences` in `TestOriginRegression` |

## Open limitations

- Straddle-touching scenarios (different pinned values sharing a
  boundary value) are skipped loudly, not solved — full handling needs
  joint boundary-count accounting across chained realizations.
