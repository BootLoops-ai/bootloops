# lattice/emit — stratum emission engines + generic sweep chassis

## Modules + batteries
| module | what it does | how to run / battery | output |
|---|---|---|---|
| harness.jl | glue-list emission core (tau_extract / entry_validate / canonical_key / emit_rank1 / emit_rank2_definite / completeness_check) with the I1-I4 / AppB / SOUND44 exact controls | `julia harness.jl all` emits the trace-<=48 list and runs the controls (needs an Oscar/Hecke env; outputs are regenerated, not shipped) | G48_LIST.jsonl + controls_report.json next to the script |
| s11.jl | rank-2 sig(1,1) certified enumeration engine; V1/V2/V3 gates + MUT mutation control | `julia s11.jl all` (needs an Oscar/Hecke env; V2/V3/MUT run on shipped fixtures; the V1 byte-exact comparison reads a reference list of the (U,U) cell that is not included in the package — set TERRIER_SWEEP_BANK to its root) | R2S11_*_ENGINE.jsonl + ENGINE_REPORT.json next to the script |
| predicates.jl | exact side_analysis / perp_root_check / congruence_diag predicate chain (module SweepPredicates, check_flux P1-P7) | exercised by both engines above | — |
| selftest_emit.py | the SHIPPED battery: exact-arithmetic verification of the glue predicate chain on bank/sound_minvecs.json (44-row known-answer table), a hand-derived closed-form case, and mutation must-fails | `python3 selftest_emit.py` (pure Python, seconds) | stdout gates |
| sweepchassis.py | generic sweep chassis: evaluate(req) + adjudicate(...) + controls, hash laws (CTRL-ID-V1, CTRL-POS-V1), atomic receipts, resume, void-on-miss rc=8, HALT law rc=9, replay, single-writer lock, injectable clock | `python3 selftest_sweepchassis.py` (SC1-SC6 on a toy stream; SC7 compares the hash laws against a reference shard receipt that is not included in the package — set TERRIER_SWEEP_SHARD_RECEIPT) | SELFTEST_SWEEPCHASSIS.json next to the script |
| bank/ | shipped fixtures: sound_minvecs.json (44-row known-answer glue table), k3xk3_flux_N_Q25_arxiv2010.10519_appB.json (published flux, arXiv:2010.10519 App B), s11_expected_counts.json (expected entry counts for the s11.jl V2/V3 cells) | sha256-pinned in ../../regression_manifest.json | — |

## Generalizations (each file header lists its own interface)
- harness.jl: pluggable predicate interface — define `PREDICATE_HOOKS = (side_analysis=..,
  perp_root_check=.., congruence_diag=..)` in Main BEFORE including to run the emission
  core under a different exact predicate chain. Defaults = SweepPredicates. The
  I1-I4/AppB/SOUND44 controls are the K3xK3 reference instance and double as the
  engine battery. entry_validate is BYTE-VERBATIM (s11.jl verbatim_check source-diffs
  it against this file at every launch).
- sweepchassis.py: sweep-driver machinery for ANY (entry stream, predicate chain).
  Injectable clock => replayed shards are byte-identical files.

## Delegation map (never duplicate)
../lattice.jl for BQF/theta/genus; ../nikulin.jl for disc-form embedding certs;
../../common/controls.py for the CTRL-ID-V1 / CTRL-POS-V1 hash primitives.

## Run law
julia: nice -n 5 julia +1.10 --project=<your Oscar env> -t 1 <file>
python: ulimit -v 32505856, nice >= 5.
