# frobenius_boundary — branch boundary at var=inf for cut-block DEs

## When to use
Any single-variable transport whose boundary anchor sits at infinity of a
mass-like variable and whose physical solution is NOT analytic there — the
target pattern: a cut block (maximal-cut DE) where the naive vacuum boundary
vanishes and the solution lives on non-analytic Frobenius branches. The tool
pins which branches can carry the physical solution (often cutting the free
boundary constants by more than half) before any quadrature/oracle work.

## Pipeline (cli.py, config-driven; see config_example.json)
1. `build_poly_DE` — exact `D(x)·dM/dx = G(x)·M` from a per-point Kira table
   (scalar LCD at rational d + DFT fit of G; fit relerr ~1e-77 on the reference
   family — the measured example used throughout; its table is not shipped).
2. `spectrum` — exact Laurent of `x·A` at infinity (poly division, no
   Vandermonde), `D1 = R_Delta − diag(L·d/2 + ai)`; exponents = eig(D1)
   ONLY if the N-frame is Fuchsian.  When a GENUINE nilpotent D_0 forward
   block is detected (fit-noise blocks ~1e-120 dropped; rel-1e-45 floor),
   the cli also shear-reduces per eps (`_reduce_fuchs`, ~12.6 s/eps on
   the reference family) and takes eig of the reduced residue Dt1 — the TRUE exponent
   multiset (independently verified by disjoint methods).
3. `classify` — pair spectra across ≥2 eps (verified at ALL given eps — the
   4-eps pattern), rationalize lambda = a + b·eps, Jordan via SVD nullity,
   trace rule tr D1 = rational-linear in eps (build-integrity gate).  The
   cli classifies the REDUCED spectrum when a genuine D_0 exists
   (`classification_frame` in the output names the frame; the eig(D1)
   classification is kept as `groups_eig_D1`); Fuchsian systems unchanged.
4. `exclude_strata` — survivors; killed branches hardened by a no-rational-
   within-1e-25 (q ≤ 1e6) scan at 2 eps; complex pairs by |Im| > tol.
5. `branch_series` — D_p-block Frobenius recursion seeded by eig-vectors;
   at integer resonances the backward path solves the compatible/min-norm
   system (SVD; incompatibility asserted numerically and reported as
   refine_rfinal / LOG_BRANCH_REQUIRED; homogeneous freedom stays with the
   resonant eigenvalue's seeds — span of Phi unchanged; a naive direct solve
   there is a ZeroDivisionError).  A GENUINE nilpotent D_0 forward block
   (fit-noise blocks dropped; e.g. D_0^2=0 rank 5) is removed by an exact
   shear reduction and the branch is built by the same backward recursion in
   the Fuchsian reduced frame, then mapped back (diag: fwd_reduced,
   orbit_multiset); a forward-coupled dense lstsq survives only as a
   fallback (fdepth>1), flagged reduce_failed.  `phi_matrix(Mstart)`
   gives the N×J branch-solution matrix for the anchor solve
   `R_rows·Phi·c = O_rows`.
   Accel envs (default off, numerics unchanged):
   FROBENIUS_POLYDE_CACHE=<dir> disk-caches the fitted poly-DE;
   FROBENIUS_SAMPLE_NCPU=<n> fork-parallel D*A sampling.

## Convention
`M_k(var) ~ var^(L·d/2 + ai_k + lambda_j)`, `lambda_j = eig(D1)`,
`ai = (#ISP powers) − (sum of dots)`, eps = (4−d)/2, L = n_loops.
N-frame: `N = diag(var^−(L·d/2+ai))·M ~ var^lambda`.

## PROVEN vs ASSUMED in the strata exclusion (do not overclaim)
PROVEN (exact): the DE and D1 (Kira IBP + exact poly division; trace rule);
the exponent spectrum and its rational/non-rational classification at the
given eps set; for the reference family also the vacuum-sector emptiness (kira_vacuum IBP)
that kills every integer eps-independent branch.
ASSUMED (flagged in output): (i) the rationalized (a,b) of survivors are
exact (multi-eps fit at ≤1e-26, not a symbolic derivation); (ii) the kill of
eps-incommensurate branches rests on the single-Laurent-eps argument — the
physical cut's large-var expansion contains only `var^(p+q·eps)` strata with
RATIONAL p,q (expansion-by-regions completeness + quasi-unipotency of the
twisted-period local system); (iii) `RATIONAL_EPS_INDEP__FLAG_review`
members are ALLOWED-but-unkilled: no recorded argument excludes them (on
the reference family's true spectrum this class fires ×5 — half-integer eps-indep).
`kill_integer_eps_indep` must only be set when the vacuum theorem holds for
YOUR family (verify with a kira_vacuum IBP run).
Overdetermination of the downstream anchor solve (extra moment rows beyond
the survivor count) is the empirical test of (ii) — keep those rows.
SPECTRUM LAW: with a genuine D_0 block the
TRUE exponent multiset = eig of the shear-reduced residue, NOT eig(D1); the
classify/exclusion stages consume the reduced spectrum automatically
(reference family: 30 → 22 survivors + 5 flagged + 3 killed = 27 allowed; a naive
"30→12" reading of eig(D1) is WRONG as a counting claim — the 18
complex/irrational eig(D1) values are exponents of NO solution).  Verified
independently by disjoint methods (monodromy + recursion rank-profiles).
The reduced-path seed→orbit-column fetch is injective
(seed_index; an argmax fallback can collide on a degenerate family and
duplicate a Phi column — rank 11 instead of 12); the selftest enforces
rank(Phi) = column count.  Downstream anchor solves must size against the
27 allowed directions or keep enough moment rows to MEASURE the extra
coefficients' vanishing (empirical, not a theorem).

## Run / regression
`cd tools/frobenius-boundary && python3 -m frobenius_boundary.selftest` — ACCEPTANCE: reproduces
the reference-family numbers from its kira_targets_all.m reduction table (not
included in this repo): degD=38, degG=40,
genuine D_0 (rank 5, 3 shear steps), reduced-frame counting 22 survivors
{1−eps(×3), 1−2eps(×6), 2−eps(×2), 7/4+eps/2(×3), 9/4+eps/2(×5),
11/4+eps/2(×2), 7/2+2eps} + 5 flagged {3/2(×3), 5/2(×2)} + 3 killed {2(×3)},
tr D1 = 69−10eps, tr Dt1 = 54−10eps, the 12 reference branch series (eig(D1)
request labels) + Phi(196) at FULL column rank 12, injective seed→column map.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).
