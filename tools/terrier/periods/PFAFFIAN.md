# periods/pfaffian.py — Pfaffian matrix transport (certified vector at target)

Worked example in `pipeline/pfaffian/`: `DESIGN.md` (the law) and `ads581/`
(route A: ads-5-81-3213, rank-17, W0 verdict). The route-B worked example
(a rank-16 conifold-curve card with a certified vector at the target) is not
included in this release. The full-size reference receipts the examples
produce are not included in the package; the battery reads them from
`TERRIER_PFAFFIAN_BANK`.

## Front door — route choice BEFORE any elimination
`route_choice(target, probe=None, budget=20000)`
* `target="point-value"` -> always `matrix-transport`: keep v' = A(s) v
  (rank = module rank; entries measured height-tame: ads-5-81 deg<=574/97-102
  bit heights at rank 17) and transport
  the SYSTEM. The restricted scalar operator is the ELIMINATION and blows up
  (measured: order 13, z-degree > 4052 from that same tame connection).
* `target="scalar-operator"` WITHOUT a mod-p (order, degree) probe ->
  `ScalarEliminationRefusal` (price with tools/pf_rank.jl / Krylov FIRST).
  With a probe: refused when order*(degree+1) > budget, and `.advice` spells
  out the matrix route. Budget default 2e4 slots = above every scalar we
  have actually landed, below the measured ads-5-81 blow-up (5.3e4).
* Stored exact series are HELD-OUT gates, never derivation inputs.

## The certified chain (generalized: system samples | GKZ ideal + curve)
1. SAMPLE mod p per node — closure solve is ALGEBRAIC in s0, no series
   enumeration (the shipped sampler: `ads581/a1_perprime.py`; route B uses
   a two-pass conifold-frame sampler with a boundary-leak assert and a basis
   from closure certificates, DEGMAX ladder, rho <= 30 — its scripts are not
   included in this release).
2. `lift_connection(banks, primes)` — per-coeff CRT + Wang; acceptance =
   0 Wang fails AND k-1 stability gate (all-but-last-prime lift agrees on a
   5% subsample). Adaptive prime schedule; stop at k > 32 (heights are
   then scalar-class -> re-examine the basis, don't buy primes).
3. `fresh_prime_gate(A_exact, P, sig, A)` — 0 diffs at a prime NEVER used
   in the lift (PA-1/PB-2). The mutation control in the battery proves this
   gate catches a single perturbed connection coefficient.
4. `series_annihilation_gate(A_exact, rows, M)` — exact annihilation over Q
   of the held-out stored series, all rows, all orders <= M (PA-2/PB-3).
   TRANSPOSE LAW: theta_s v_j = sum_i A_ij v_i (module-side connection; the
   function vector contracts through A^T). Probe the orientation on a short
   prefix first — untransposed fails loudly at orders 2-5.
5. LAND inside a proven disk (both routes land by direct summation —
   marching is used only for ordinary-point legs):
   * Route A `direct_sum_periods` + `contract_w0`: towers from the stored
     F3 dictionary fit (12-block, rank-saturating, held-out-gated;
     `frame_towers` rebuilds them from the frame-coefficient file), closed-form tail
     balls (`bound_level_ads581`, `tail_ball`), flux-frame symplectic
     contraction, two-dps digit gate.
   * Route B `entropy_landing`: entropy-majorant radius (concave f + Euler
     identity => zero-constant gradient majorant, exact vertex max over the
     level-1 support polytope -> e^{-Phi1}), geometric tail, and the in-ball
     ODE residual binding the landed vector to the gated connection.
   Ball law: arb END-TO-END; never reduce a ball through float64.

## Known structural facts (details in `pipeline/pfaffian/DESIGN.md`)
* Frame poles at MUM: the exact A(s) can have a high-order apparent pole at
  s=0 in the closure basis (order 19 for ads-5-81) — Frobenius seeding there
  is frame-blocked (Moser reduction = named tool gap); hence direct summation.
* An apparent singularity ON the transport segment makes a naive fixed-eps
  marcher Zeno-stall; the convergence-disk certificate sidesteps it.
* Fail-closed scope: no verdict without EVERY gate (PA-1..6 / PB-1..6);
  ads-5-81 W0 verdict inherits the C2 half-D3 scope string.

## Battery (`selftest_pfaffian.py`, listed in regression_manifest.json)
Needs the reference receipts (not included in the package; set
`TERRIER_PFAFFIAN_BANK`). P0 route front door (refusals + acceptances); P2
fresh-prime gate replay (ads-5-81 p=2147483489, 512 nodes, 0 diffs); P3
MUTATION control — one perturbed coefficient must be caught; P6 ads-5-81
certified W0 ball dps 60+90 == stored verdict byte-exact (midpoint digits
AND radius). The route-B legs (lift replay, annihilation replay, landing
balls) are not part of this release.
