# census/ — certified flux-census wing
Machinery for counting flux vacua by SL(2,Z) x monodromy orbit class with
certificates at every step: a complete orbit fold, two independent dedup
engines with a witness referee, certified distances to special loci, an
exact-rational control sampler with Clopper-Pearson bounds, an ISD
finiteness cutoff, and a continuous-exclusion tile engine.

Geometry labels used across this wing: G4 = the AESZ 34 one-parameter
family, G5 = the mirror of the Reye congruence, G6 = the Hulek-Verrill
fourfold-derived family HV4 (see `dist_cert/zchart.py`, `isd_cutoff/` and
`tile_engine/` respectively for the charts each uses).

## The pipeline (one line)
fold -> dedup two-stack -> adjudicate -> dist_cert -> control/CP -> verdicts;
tile route for CONTINUOUS exclusion over moduli boxes.

## Stages
1. `fold_hnf/`  — HNF+detU complete SL(2,Z)-orbit key: fold the raw flux window
   before any counting. On the G4 anchor set: 256/256 witnessed corner merges;
   416->160 classes at B=1; acceptance 6,612/49,544 EXACT; 5 tamper
   must-fails. Battery: tests/test_fold_hnf.py T1-T4 (needs the anchor-set
   reference data, `TERRIER_G4_DIR`).
2. `dedup/`     — TWO independent engines A+B + third-stack witness.py +
   reconcile.py. Rule: disagreement = RED — record it, never quote the raw
   count (an engine that under-merges by a fraction of a percent is caught
   by the witnessed reconciliation, and the reconciled partition is
   audit-only).
   Battery: test_dedup.py 8/8 (synthetic T1-T8, proof-level truth; runs as
   shipped).
3. adjudicate   — verdict composition per bin (verdicts are emitted through
   common/verdict.py).
4. `dist_cert/` — exact Fraction balls + certified sqrt/ln + z/WP charts +
   exact root/CM predicates; Fincke-Pohst exhaustive OFF certificates.
   Battery: tests/test_dist_cert.py 37/37 (needs the plant-rule seeds,
   `TERRIER_PLANTS_DIR`).
5. `control_sampler/` — exact-rational AD rejection sampling (dyadic k/2^64
   coins, zero floats) + exact CP brackets; Hoeffding p_U labels (only-widens).
   Battery: tests/run_all.py 19/19 (needs the seeds file, `TERRIER_SEEDS_JSON`).
6. `isd_cutoff/` — tau-eliminated ISD finiteness cutoff on n2 = f.f + h.h;
   rule: the formula is fixed BEFORE any count is computed under it.
   Derivation (ISD_CUTOFF_DERIVATION.md) + numeric checks (isd_ref.py,
   sign_probe.py).
7. `tile_engine/` — continuous-exclusion route (OL-2/4/5 in tile_engine.py,
   fail-closed: missing certificate => OPEN). Needs the u1s summation card
   (periods/summation/) and its reference bank (`TERRIER_TILE_BANK`; not
   included in the package). Reference run over the 64 extreme tiles:
   0 EXCLUDED, 63 OPEN, 1 must-fail control tile (its closure contains a
   cited vacuum, which the criterion correctly cannot exclude).
`selftest_census_smoke.py` runs the whole chain end-to-end on a synthetic 2x2
window (fold -> two-stack dedup -> dist_cert -> common receipt) with no
external data. The controls discipline throughout: two-stack RED on any
disagreement, blind plants inserted before dedup that must survive as their
own class and certify OFF-locus, and counts quoted only with their witness
closure.
