# Counterweight — guide

Counterweight — `tools/counterweight` — https://bootloops.ai/tools/counterweight.html

The `canonical_form/` member (Canonify) has its own page:
https://bootloops.ai/tools/canonical-form.html

KIND: julia (engine) + python front-end member (canonical_form/)

REQUIREMENTS: Julia 1.11 or newer. The Julia deps (Nemo 0.56 + JSON) are pinned by the
committed `Manifest.toml`, generated on Julia 1.11; one-time setup from this directory:
`julia --project=. -e "using Pkg; Pkg.instantiate()"`. On a newer Julia minor it still
instantiates (Pkg then recommends a `Pkg.resolve()` to re-pin for that version). The
Python member (canonical_form/) needs only the repo's common core.

LAYOUT: this directory IS the Julia engine (`Project.toml` + `src/{Counterweight,RatFunc,
Fuchsia,Normalize,Multivar,KiraDE,MpolyFeed}.jl`, Manifest pinned in-tree).
`canonical_form/` is the Python LS/UT front-end member (Canonify): `canonical_form.py` +
`counterweight_bridge.jl`, defaulting its engine root to this directory (override with env
`COUNTERWEIGHT_ROOT`).

PURPOSE: No-Mathematica ε-factorization engine: given a Kira-built connection A(ε,x),
find a rotation T with T·A·T⁻¹+(dT)·T⁻¹ = ε·Ã, Ã ε-free dlog. Basis construction
(preprocessing for the bootstrap), NOT DE solving across kinematics. Doubles as an
impossibility certifier: a proven stop-class = no rational ε-gauge exists.

USE-WHEN:
- eps-factorization stalls / connection not Fuchsian → `moser_reduce` first
  (higher-order poles → first order; on failure returns an irregular-pole diagnosis =
  masters missing an LS prefactor).
- Balances stall early with a quadratic-letter alphabet, or exhaust with all
  ε⁰-residues nilpotent → R∞ residue (computed as R∞=−lim x·A, factorization-agnostic)
  + `_nilpotent_eps_grading` (Lee §5: U=diag(ε^{−Lᵢ}), Lᵢ = longest-path depth in the
  A0′ support DAG; fires automatically after `_constant_epsdecouple` fails). Reference
  scale: 13×13 block, 4-balance stall → certified dlog ε-form, 18+3 letters, ~22 s.
- Σa≠0 over rational poles (missing balance partners at an irreducible ε-free
  quadratic/cubic letter) → Case A″ Galois-orbit balance (`Normalize.jl`):
  T=(I−P)+f·P, f=q̂/Π(x−b), exact ℚ(ε) orbit residues via gcdx mod q̂; P = mixed
  rank-1 (joint rational orbit eigvec u, simultaneous dump left-eigvec v, 2nd-order ∞
  refinement) or rank-2 Tr-pair 2Π₀−pΠ₁; a ± ε⁰ pair sitting AT the roots themselves
  takes the same-orbit shear T_q=M(x)/q̂ (no dump, det T=1 — see FOOTGUNS).
- Connection arrives as uncancelled P/Q towers (IBP / on-curve kernels) →
  `MpolyFeed.feed_matrix_from_chunks` ALWAYS (raw string+Horner feed is dead: >40 min
  CPU on an n=4 block vs ~3 s mpoly-fed with the no-gcd Qex construction; `cache=`
  caches reduced pairs, selftest-on-load).
- Tower connection conjugated by stored per-block rotations T carrying 1/ε
  (dotted-master norms) → conjugate_tower driver pattern: apply_gauge(diag T_blocks) →
  block ε-grading U=diag(ε^{k_b}) (valuation longest-path) → `try_epsfactor`.
  `_eval_eps0` clears common ε-content of num/den first; genuine ε-poles still error;
  regular-vanishing returns 0.

NOT-FOR: solving DEs; half-integer-orbit exponent blocks (need a √q̂ prefactor in the
basis — an alphabet fact, not a tool gap; route to vopclose); apparent-orbit and n=2-K-eigvec
classes (provably beyond constant-P); irregular points after moser_reduce (fix the
basis/LS prefactor, not the engine).

BALANCE LADDER (`Normalize.jl` — how `try_epsfactor` removes the ε⁰ obstruction):
- Case A (primary, cross-pole): Lee rank-1 balance. Eigen-data at EVERY singular point
  via `factor(charpoly(R))` over ℚ(ε) — complete for Feynman residues (eigvals affine
  in ε). P = u·vᵀ/(vᵀu) with vᵀ a LEFT eigvec of the full ℚ(ε) residue at x₁ and u a
  RIGHT eigvec at x₂. The one-point spectral projector is NOT this P — it introduces
  an O(ε) double pole at x₂. Each accepted balance drops Σ|a| by 2 (Lee §4).
- Case A′ (same-pole ±, routed): when the remaining ± ε⁰-eigenvalue pairs all sit at
  one pole p — structurally unavoidable whenever Σa⁺ at p exceeds the available
  cross-pole negative units — a 2-step composite B(P₂,q,p)∘B(P₁,p,q) through any
  partner q carrying an a=0 eigval. The per-step projectors are generically ε-singular
  (a=0 eigvecs at q coalesce as ε→0, so vᵀu = O(ε)) but the composite gauge T₂T₁ is
  ε-regular; acceptance therefore gates the COMPOSITE only (Fuchsian + ε-regular +
  Σ|a| drop), never the individual steps.
- factor-ε (ℚ(ε)-eigenbasis): with all a=0 the ε⁰ residues can still be
  nonzero-nilpotent, where the rational coboundary fails (∫A₀ produces logs).
  `_factor_eps_eigenbasis` applies the constant-in-x gauge G = S⁻¹, S = column
  eigvecs of one residue over ℚ(ε), trying each pole and accepting on `is_epsform`.
  det S = O(ε^k) is fine — that is the canonical UT weight-grading; only G·A·G⁻¹
  must be ε-regular. Runs before `_constant_epsdecouple` in the driver.
Reference scale: a 3×3 block with an unavoidable same-pole surplus (Σa⁺ = +5 at one
pole vs 4 cross-pole negative units): 8 balance rounds (7 cross-pole + 1 routed) +
the factor-ε gauge, ~9 s, three-letter alphabet {y, y−1, y+1} → `certify_epsform`
Fuchsian + ε-form + dlog all true, post-transform residue eigvals all ∝ ε, and the
total gauge satisfies T·A·T⁻¹ + (dT)·T⁻¹ = ε·Ã exactly. The ladder has no
case-specific branches; cost per round is O(#poles²·n²) eigvec/nullspace solves —
≲10⁴ nullspace calls per round even at n≈50 with ~7 letters. `verbose=true` prints
one line per round (obstruction norm + the accepted move).

INVOKE (public, from this directory): `julia --project=. ` then `using Counterweight`;
`moser_reduce(ctx, A; max_iter=200, verbose=false)` (Fuchsia.jl);
`try_epsfactor(ctx, A; max_balance_rounds=400, verbose=false, ...)` (Normalize.jl);
`feed_matrix_from_chunks(ctx, n, chunks; cache=path)` + `feed_cache_save/load`
(MpolyFeed.jl); KiraDE parses Kira derive_dgl output → A(ε,s,t); Multivar handles
(s,t,ε) via one ratio variable + flatness check. From python: `canonical_form/`
`fuchsify()`/`factor_epsilon()` via counterweight_bridge.jl. Examples:
`examples/henn_dbox.jl` (self-contained); `examples/coupled_tower_counterweight.jl`
(driver pattern for a Kira derive_dgl tower — point `TOWER_DE_JSON` at your own
export) and `examples/sec121_baikov.py` (Baikov connection construction).
Tests: `test/{test_moser.jl,validate_henn.jl,test_orbit_shear.jl,test_sec53_qorbit.jl}` (the last with its
fixture set `test/fixtures/sec53_qorbit/`, see its README).

INPUTS: connection matrix over ℚ(ε)(x) (Nemo/FLINT); Kira derive_dgl files; raw P/Q
chunk lists as FLINT ℚ[eps,x] mpolys (MpolyFeed); per-block gauge rotations possibly
carrying 1/ε.

OUTPUTS: gauge T + ε-form Ã (dlog letters); Moser failure → irregular-pole diagnosis;
Case A″ → orbit-balance report with q-roots.

GATES: Fuchsia.jl certification (Fuchsian + ε-form checks); Case A″ acceptance = P²=P
+ Fuchsian + ε-regular + strict drop of the algebra-aware norm; MpolyFeed selftest =
exact BigInt-rational chunk evaluation at random (ε,x); timed pilot first: price every launch of
unknown cost with a short probe run BEFORE launching.

FOOTGUNS:
- Literal P_r+P_r̄ is NOT idempotent on real data — use the mixed rank-1 form.
- Case A″ projector moves (rank-1 joint-eigvec, rank-2 Tr-pair) transfer a-units between
  the orbit and OUTSIDE dump points only; a ± ε⁰ pair sitting AT the conjugate roots of
  an irreducible quadratic letter is invisible to them, and sequential per-root conjugate
  balances compose to an IRRATIONAL gauge. Cure (shipped: `Normalize.jl` same-orbit ±
  shear, tried automatically after the projector variants): single rational orbit move
  T_q=M(x)/q (deg M≤2; det T=1 so det M=q̂ⁿ; C(y)=γ(y)·u(y)v(y)ᵀ from cross-orthogonal
  K-eigvecs, adj(M) killing the balanced line to 2nd order at both roots; γ from the
  double-pole cancellation γ=(λ_hi−λ_lo−1)/(vᵀBu) with B the constant orbit Laurent
  term). Regression: `test/test_orbit_shear.jl`.
  Companion Liouville test: nonzero odd part of the q-residue ⇒ IRREMOVABLE — the
  alphabet needs the algebraic-root letter; pure rational dlog impossible (the shear's
  gates refuse the odd class; the engine stops with the residual obstruction).
- ε-factorized ≠ canonical: a D-rescale + row-0 VoP can leave transcendental T₀'s in
  an off-diag block (mixed-weight ε⁻²). Stationary-point test:
  d(m_UT[0,ε⁻²])|_basepoint=0 in all directions ⇒ NOT weight-1. Im=const on a
  Euclidean path is kernel-reality, not a UT certificate.
- Decouple-loop guard is in place (bounded rounds); `//` parser round-trips.

BATTERY (from this directory):
`julia --project=. test/validate_henn.jl` → VALIDATION GATE (a) PASS (Henn
massless double box ε-form from a scrambled basis);
`julia --project=. test/test_moser.jl` → 4/4 PASS;
`julia --project=. test/test_orbit_shear.jl` → 5/5 PASS (same-orbit ± shear at an
irreducible quadratic letter: K-arithmetic pin, move fires + invariants, public-driver
wiring, odd-class refusal control, wrong-γ mutation control);
`julia --project=. test/test_sec53_qorbit.jl` → 12/12 legs as expected (a 2×2 block with
an irreducible quadratic letter: the engine reaches a rational ε-form in the rationalizing
chart and STOPS by name in two charts with half-integer exponents; the vendored exact point
verifier passes the fresh and the two pinned ε-forms, fails the three planted controls and
the row-major reading; ~20 s standalone);
`python3 canonical_form/test_canonical_form.py` → 5/5 PASS
(T3's kite cross-check and T4's oracle wire skip cleanly when their optional
data is absent).
Registered battery (tools/BATTERIES.json): `julia --project=. test/runtests.jl` —
the shear, Henn and quadratic-letter legs in ONE Julia process (startup + Nemo load +
JIT paid once; measured ~30 s wall on a many-core Linux host, well inside run_selftests.py's
per-attempt budget).

RELATED: upstream Kira (derive_dgl → KiraDE); canonical_form/ (LS/UT rotation
front-end, bridges here for fuchsify/factor_epsilon); MpolyFeed is the Julia sibling
of `wayfinder.flintexport` (tools/wayfinder/flintexport.py); blocked classes route to vopclose (function-level closure).
External tools in the same niche — CANONICA and Libra — are Mathematica-only with no
Julia/Python port; the public fuchsia (Gituliar–Magerya, Python/Maxima) is the nearest
open-source counterpart. None is a dependency here.

CREDIT: the canonical ε-form is Henn's [Henn]; the reduction algorithm (fuchsification,
balances, ε-factorization; §5 nilpotent grading) is Roman N. Lee's [Lee], with Moser's
and Barkatou's rank-reduction theory beneath it (classical); `src/Fuchsia.jl` is our
independent Julia/Nemo implementation of Lee's algorithm. It shares its name — and its
algorithm, though no code — with the public Fuchsia of Oleksandr Gituliar and Vitaly
Magerya [Fuchsia], which we acknowledge as the first open implementation; related public
implementations are epsilon (Mario Prausa [epsilon]), CANONICA (Christoph Meyer
[CANONICA]), Libra (R. N. Lee [Libra]) and INITIAL (Christoph Dlapa, Johannes Henn, Kai
Yan [INITIAL]). The LS/UT front end follows Arkani-Hamed–Bourjaily–Cachazo–Trnka [ABCT]
and Henn–Mistlberger–Smirnov–Wasser [DlogBasis]; connections from Kira [Kira3].
Bracketed keys resolve in REFERENCES.md at the repository root.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).
