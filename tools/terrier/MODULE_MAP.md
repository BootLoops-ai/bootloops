# MODULE_MAP — module index: every shipped module, what it does, its gate battery

Battery = the named gate set that `selftest.py` runs in this tree. A battery
whose reference inputs are not included in the package prints a named SKIP
row with the `TERRIER_*` environment variable that supplies them (see
README.md, "Environment"); with the variable set it runs for real.

## periods/ wing — certified period suite
Certified CY3 period/transport/root-finding chassis: gated cards (toric/GV ->
effective Gamma-series/GKZ), exact flux contraction + frame dictionary
(xi-toggle as code), certified MUM+conifold transport (ball arithmetic
end-to-end — never reduce a certified ball through float64), Krawczyk
existence+uniqueness with FD-vs-AD Jacobian gate, two-dps verdict chassis
with mutation controls and code-sha receipts.

| module | what it does | gate battery |
|---|---|---|
| transport.py | MUM transport chassis + Krawczyk vac layer (`transport.vac_layer()`) | battery_dkmm.sh 13/13 on the two-modulus KKLT example (needs TERRIER_KKLT_BANK); replay balls byte-match manifest pins |
| fluxcurves.py | exact flux contraction / frame layer (xi_flipped toggle as code) | G2 exact + G4 frame/f3-mutation/towers_ext (needs TERRIER_KKLT_BANK); replay outputs byte-match pins |
| cards.py | gated card contract (fail-closed PENDING-CARD semantics) | selftest_cards.py GREEN (shipped fixture cards validate incl. slice-operator-v1; G2 exact) |
| geffseries.py | effective Gamma-series facade (delegates to pipeline/geff_series.py) | g1: 901 w0 coeffs exact, 48 gv exact, racetrack N_m match (needs TERRIER_KKLT_BANK) |
| opderive.py | CRT/Wang + annihilation-proof gate (wrapper; engine = tools/annihilator/) | g3: theta-form == operator_LS.json exactly; 838 windows annihilated, 316 held-out (needs TERRIER_KKLT_BANK) |
| envelope_certified.py | PROVEN tail majorant (lemma PROVEN_MAJORANT.md; certified subdisk 1/t_min law) | selftest_envelope.py E1-E6 (needs TERRIER_KKLT_BANK) |
| pfaffian.py + PFAFFIAN.md | matrix transport: (system samples \| GKZ ideal + curve) -> certified vector at target; route_choice front door with ScalarEliminationRefusal | selftest_pfaffian.py with the reference receipts (TERRIER_PFAFFIAN_BANK): ads-5-81 W0 ball dps 60+90 BYTE-EXACT; fresh-prime gate 0 diffs (147968 values); mutation control must be caught |
| conifold.py | coni_pfv CP0-CP10 chart + towers + fail-closed transport | selftest_conifold.py 6/6: DKMM NOT-CONI routing control; DKMM tower s^80; branch-flip/branch-swap mutations caught; transport fail-closed |
| gvderive.py | derive + certify GV working sets from toric data (GV-band labeling law) | none registered in this release (its worked example is not included in the package) |
| ffp.py | finite-field fingerprint engines (all-z Frobenius scan; int64 front door, BigPrimeRefusal) | selftest_ffp.py 12/12: tiny-q brute exact; HV p53 + BCM p61 pass-sets REGENERATED == shipped reference; mutation caught |
| sunit.py | refined S-integrality candidate generator (machine-readable truncation statements) | selftest_sunit.py 9/9: planted-point gate with HAND-DERIVED S_min exact (shipped synthetic control card); must-fail controls |
| summation/ | u1s summation-route modules (jet towers, entropy certs, Cauchy floors; entropy radius law in its README) | selftest_summation_reduced_u1s.py (needs the box-operator reference data, TERRIER_BOXOPS_BANK; named SKIP otherwise) |
| pipeline/cards/ | card bank incl. slice-operator-v1 class (hv4-diag-L5) | battery_l5.py PASS; selftest_cards.py GREEN |

## common/ chassis
| module | what it does | gate battery |
|---|---|---|
| verdict.py + periods/verdicts.py | two-dps verdict chassis, code-sha stamping, published-value gate protocol | selftest_verdict GREEN; full 13/13 scoreboard = selftest_verdicts.py (needs TERRIER_KKLT_BANK) |
| receipts.py | atomic tmp+fsync+rename, sha-stream, code-sha stamps, resume law | selftest_receipts.py 12/12 |
| controls.py | hash-placed blind controls, mutation + planted-error harnesses | selftest_controls.py 14/14 |
| frames.py | convention-registry template (evidence classes, toggles, tamper-evident freeze) | selftest_frames.py 9/9 |
| certs.py | two-dps / ball-honesty / published-value / float64-trim laws | selftest_certs.py 22/22 |
| dedup.py | known-point pullback-orbit dedup + newform backstop (KNOWN_ORBIT.jsonl packaged; explicit-receipt law) | selftest_dedup.py 5/5 |

## lattice/ wing — flux-lattice enumeration
Frozen conventions any consumer inherits (restated in lattice.jl header): root
lattices POSITIVE definite; disc form q into Q/2Z (Nikulin); form (a,b,c) <->
Gram [2a b; b 2c] (factor-2 trap); SL2 (proper) reduction; imprimitive classes
real; vector counts both signs; indefinite box convention.

| module | what it does | gate battery |
|---|---|---|
| lattice.jl + validate.jl | definite-BQF exact SL2 reduction + class enumeration (covers the Hecke d<0 NotImplemented gap; proper-vs-improper distinction kept) | 18/18 battery (validate.jl writes its results report and is re-runnable); B2 ideal-class-group cross-check |
| nikulin.jl | disc-form Nikulin certificates: Milgram self-test + glue-rigidity via image_in_Oq (no genus enumeration itself; sig(0,3) route avoids the sig(3,0) primitive_embeddings hang) | selftest_nikulin.jl: pins 3/3 pass (needs TERRIER_JULIA_PROJECT; writes its RUN_pins.json receipt at run time) |
| emit/ (harness.jl, s11.jl, predicates.jl) | emission engines with pluggable PREDICATE_HOOKS (emit/README.md) | selftest_emit.py 3/3 (44-row known-answer glue table exact in pure Python; closed-form case; mutations caught); engine regeneration needs an Oscar/Hecke env |
| emit/sweepchassis.py | generalized sweep chassis for any (entry stream, predicate chain); hid/hpos delegated to common/controls.py | selftest_sweepchassis.py 11/11 with a reference shard receipt (TERRIER_SWEEP_SHARD_RECEIPT): replay BYTE-IDENTICAL; planted-error rc=8; halt rc=9 |
| genus_enum.jl | Aut-free Kneser genus enumeration + PARI POST mass certificates (fills the Hecke/GAP OOM gap) | parse-checked; A1^6 control byte-match; Oscar control replay needs TERRIER_JULIA_PROJECT and TERRIER_GENUS_CONTROL_RECEIPT |
| lp_kill_rank_agnostic.py | Venkov LP kill, rank-agnostic (n = len(GK) from input) | A1^6 control byte-match + d4 replay JSON == reference (TERRIER_R18_DIR); contract = LP_KILL_RANK_GENERALIZATION.md |
| lattice_theta.py | Jacobi theta_{2,3,4}, Theta series for Z^n, D_n, D_n^+, E8, Leech, GKP flatness eps, vectorized CVP (standalone NumPy member — does NOT inherit this wing's frozen K3 conventions; theta argument is the nome q, not tau) | none registered (standalone module; docstring worked values) |

## census/ wing — certified flux census
| module | what it does | gate battery |
|---|---|---|
| fold_hnf/ | HNF fold stage | test_fold_hnf.py T1-T4 (needs TERRIER_G4_DIR) |
| dedup/ | two-stack A/B/witness/reconcile dedup | test_dedup.py 8/8 (synthetic, proof-level truth) |
| dist_cert/ | distance certificates | test_dist_cert.py 37/37 (needs TERRIER_PLANTS_DIR) |
| control_sampler/ | control/CP sampling | tests/run_all.py 19/19 (needs TERRIER_SEEDS_JSON) |
| isd_cutoff/ | isd cutoff formula law (formula committed BEFORE count; referee probes are instruments, not gates) | — |
| tile_engine/ | continuous-exclusion tile route | needs the summation card and external reference data (TERRIER_TILE_BANK) |
