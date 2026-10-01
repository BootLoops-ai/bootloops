# control_sampler — index-density rejection sampler with exact Clopper-Pearson bounds

Ashok-Douglas index-density REJECTION sampler on the census fundamental
domain, plus the exact Clopper-Pearson (CP) machinery for the control
on-locus fraction F_ctl. The control null hypothesis (called H0-a in the
code) is that flux vacua are distributed on moduli space by the Ashok-Douglas
index density; the sampler draws the control population under that null and
the CP machinery bounds its on-locus fraction. This module touches NO census
data: the test battery is synthetic-only, at seeds fixed in advance.

## Definitions implemented
- Density [AD hep-th/0307049 eq (4.10), kappa = -1 pinned]: integrand
  det(-R - omega 1), evaluated EXACTLY (Fraction) from exact-rational
  R(z), omega(z) callables (`ad_density.py`). pi^-n and the coordinate
  measure cancel in rejection ratios and in F_ctl; never floated.
  Negative index density = data error -> halt (sign law), never |.|.
- SAME-fundamental-domain law: the sampler CONSUMES the domain object
  exported by the dedup module (Engines A+B contract in `domain.py`); it
  never builds its own domain for a census run.
- NO-EXCLUSION-ZONE guarantee: structural — `rejection.py` cannot even
  name locus data (AST identifier scan in tests) — and statistical
  (closed-form tube occupancy at fixed seeds).
- Seeds: per-bin stream =
  PCG64(master XOR int(sha256("geometry|bin").hexdigest()[0:16], 16));
  the reading is documented in `seeds.py`. The master seed is read from a
  seeds JSON file (key `control_sampler_master`; not included in the
  package — the seed tests need one supplied with `TERRIER_SEEDS_JSON`).
  Proposals and coins are exact
  dyadic rationals k/2^64 — the whole accept/reject path is
  exact-rational (no float anywhere in the sampler).
- CP-on-F_ctl (the pinned upper bound on F_ctl): Clopper-Pearson by exact
  binomial arithmetic — Fraction CDF + dyadic bisection returning a
  certified bracket [lo, hi]; the quoted p_U is the conservative end.
  Also: exceedance p^N (exact) and the INSUFFICIENT-N power gate
  min_n_for_significance (alpha0 = 1e-3 pinned). For production
  N_ctl = 100,000 pass a coarser tol (e.g. 2^-20): still exact and
  conservative, cheaper bigints.

## Files
- `seeds.py` — fixed-seed stream derivation (PCG64), dyadic uniforms
- `domain.py` — fundamental-domain interface + synthetic test domains
- `ad_density.py` — exact det(-R - omega 1), KAPPA = -1 pinned
- `rejection.py` — exact-rational rejection sampler (halt discipline)
- `cp_bound.py` — CP brackets, exceedance, power gate
- `tests/` — 19 tests: seed reproducibility + golden draws; CP vs
  published tables (Clopper-Pearson 1934 closed forms k=0; standard
  4-dp two-sided values; scipy beta-quantile cross-check); closed-form
  acceptance rates 9/16 (box) and 7/32 (triangle) exact + sampled;
  envelope/sign halts; no-exclusion-zone structural + tube occupancy
  at closed-form fractions 35/2304 and 5/384.

## Run (single core, ~25 s; needs TERRIER_SEEDS_JSON)
    python3 tests/run_all.py
