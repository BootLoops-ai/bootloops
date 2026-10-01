# Adversarial rigor audit — record

Adversarial audit of every certified-arithmetic component: 7 independent
adversarial lemma audits + a conventions audit (equation by equation against the published papers)
+ an adversarial verification pass for every rigor-relevant claim, plus two independent
oracles (GiNaC 1.8.7, Zenodo-2502.00118 ancillaries). All confirmed findings are fixed in the
shipped code, and the full test suite passes.

Two further review passes ran to completion: the Zenodo report (findings
extracted independently — see vendor/zenodo-2502.00118/REPORT.md) and the iterated.jl path-lemma referee
(result appended below). Independent of that referee, the iterated-integral code path is validated
end-to-end by the GiNaC oracle (depth-4 words with log kernels, ~152-digit agreement at
three points) and by the AMFlow ε⁷ anchor (deep log-stack words, 128-digit agreement).

## Confirmed findings (all fixed)

| # | severity | component | finding | fix |
|---|---|---|---|---|
| 1 | rigor-breaking | modular.jl `eta_quotient_qseries` | tail bound double-counted x0^w (installed the shifted majorant mass on the unshifted series, then `qshift` multiplied by x0^w again) → final tailF = x0^{2w}F instead of x0^w F: under-estimate by 4^w, demonstrably non-enclosing in adversarial setups | majorant mass computed for the unshifted Euler product; `qshift` applies the shift once |
| 2 | rigor-breaking | qseries.jl `qdq` | stored derivative radius x0′ = up(θ·x0) can exceed θ·x0 (Mag rounds up), while the tail bound used θ^{N+1}: invariant strictly violated at the ulp level (compounds as ((1+ε)θ)^N) | tail bound recomputed from the up-rounded ratio θ_ub = up(x0′/x0), with certified θ_ub < 0.91 and N ≥ 11 monotonicity check |
| 3 | rigor-breaking | qseries.jl `qsubst_pow` | silently assumed x0 ≤ 1 (uses x0^{dn} ≤ x0ⁿ); non-containing balls demonstrated for contexts with x0 > 1 | SeriesContext inner constructor now rejects x0 ≥ 1 (also load-bearing for `evaluate`) |
| 4 | wrong-result | curve.jl `incomplete_edge_quadrature` | returns the negative of the documented +i0-edge integral when the pair is passed as (larger, smaller) | certified orientation contract: requires real(a) < real(b), else error |
| 5 | wrong-result | puncture.jl `abel_image` | silently delegated to the t+iδ deformation when the marked point lies on the arcsin cut — returned ball does not contain the documented (+i0 limit) target | now errors with instructions; the deformation is explicit opt-in (`abel_image_deformed`), whose docstring and dictionary entry name record δ |
| 6 | wrong-result | dictionaries.jl `interior_dictionary` | `ψ0p, _ = _acm_psi0_with_derivative(s)` bound the period, not the derivative (destructuring order) — "dpsi0" entry stored ψ0 | take the second slot; verified numerically against the certified dK/dm chain |
| 7 | wrong-result | gamma16.jl `tau_from_t` | principal-branch K-ratio silently leaves the cusp-connected Γ₁(6) orbit for t ≲ −0.47 (most of the Euclidean axis) with only an Im τ > 0 guard | mandatory certified Hauptmodul self-check (was already the reason `sunrise2` uses the interval-Newton inversion instead) |
| 8 | wrong-result | docs/conventions.md | boundary-value equation transcribed with the γ-stripped LHS mislabelled as S₁₁₁ itself (the implementation was already correct — this same convention caused the ε¹ bug found and fixed against the AMFlow anchor during development) | doc corrected; γ-stripping now stated explicitly |
| 9 | wrong-result (stale) | puncture.jl `third_kind_G` (quadrature version) | the literal quadrature of the paper's definition of G provably never converges (P4(x−t) has a root at the endpoint; analyticity hypothesis false on all real Euclidean inputs) | superseded before audit completion by the Zenodo-PDE definition (∂ₜG one-form, anchored at G(s,0)=0, desingularized endpoint); validated against the printed series and threshold constants |

## Minor / accepted-risk findings

- `polygeom_tail` precondition k ≥ 0 now asserted (bound invalid for k < 0).
- `qshift` rejects negative shifts (previously out-of-bounds writes); `f^1` fixed via `Base.copy`;
  `qcoeff` now errors for n > N (stored zero is not an enclosure there).
- `transport` gained a step limit (no termination guarantee for wide input balls).
- `third_kind_G` entry hypotheses (s>0, −s<t<0, t(s+t)<4s) now checked at entry; inside the
  integrand they were already enforced boxwise via the analytic flag.
- `curve_periods` uses sqrt(lead) (phase kept) instead of sqrt(|lead|).
- `hyp2f1_at_r3` default series order raised so the start series does not cap output precision.
- conventions.md: generalized-Eisenstein q_M footnote (the paper's printed q_M = e^{2πiτ/M} is
  inconsistent with the printed Γ₁(6) identities; integer-q convention is operative, as in
  GiNaC); §10 note distinguishing the Bessel-moment normalization (3/4)L(χ₋₃,2) from
  S₁₁₁(2,0) = 3L(χ₋₃,2).
- `abel_image_deformed`: the O(δ·∂Z) distance to the +i0 limit is documented, not enclosed —
  a deliberate, explicit design choice recorded in the entry names of the dictionaries.
- Audit-noted heuristics that are NOT load-bearing for rigor (correctness enforced elsewhere):
  contour-radius selection in `period_quadrature` (boxwise analytic flags carry the rigor);
  Float64 conversions of certified bounds used only to pick parameters.

## Components verified sound (no rigor findings)

- `hypgeom.jl` ₂F₁ ε-series start: every ℓ₁ inequality re-derived, truncated-inverse geometric
  bound confirmed, path distances to {0,1} certified, conjugate trick sound.
- `frobenius.jl` geometric-envelope transport: induction, window width, factorial-ratio bounds,
  ℓ₁ submultiplicativity (with exact order-by-order ε-truncation), localize() re-expansion,
  evaluation tails — could not be refuted analytically or empirically.
- `modular.jl` Eisenstein series: divisor-sum loops, constant terms, coefficient bounds, tails.
- `qseries.jl` core ops (+, −, scalar, mul half-mass lemma, qshift m≥0, evaluate): sound.
- `curve.jl`: root realness-tightening (conjugation-isolation), `acm_periods`/derivative chains/
  `picard_fuchs_residual_acm` (full symbolic re-derivation), branch constructions `_pair_sqrt`,
  `_pair_sqrt_outside` (cut characterizations proved).
- `sunrise.jl`/`epsseries.jl`: assembly verified equation by equation against the published AW/ABW results; `_align`
  truncation conservative in all lead/length combinations; independent reproduction of the
  ε-expansion structure.
- conventions.md §§1–9: the transcribed equations agree with the published papers (one transcription issue, #8).

## Oracles

- GiNaC 1.8.7 (`ginac-oracle/`): all 10 kernel/iterated-integral comparisons agree to ~152
  digits; Sebbar-polynomial identity for f₃ confirmed through q³⁸ by independent machinery.
- Zenodo 10.5281/zenodo.14733100 (`vendor/zenodo-2502.00118/`): exact G differential relations
  adopted as the G definition; threshold constants a₀ = ϖ₀(4m²,m²), a₁ = −∂ₛϖ₀(4m²,m²)
  reproduced to all 68 printed digits; ω₀ mass-ODE and s-relation verified.

## Appended: iterated.jl path-lemma referee

**Verdict: SOUND — enclosure guarantee stands.** Re-derived by hand: the exact monomial/log
integration formula, the i-sum level distribution, the regularized constant term, the path
lemma (including the X ≥ 1 safeguard and the complex-q triangle-inequality direction), the
remainder-record exponent j (constant terms integrate exactly and never enter records),
the fold order vs the AW convention, and all Mag/Arb rounding directions. Empirical
containment: log-stack words I(1)…I(1⁵), geometric kernels with exact tails at
|q|/x0 = 0.88, pure-tail kernels (records as sole truth carrier), shuffle identities
within balls, N-convergence, and the GiNaC oracle word I(1,f₄,1,f₃;q_t) reproduced to
3.06e-155.

Minor findings (no action required): U_j coefficient tails are double-counted (per-level
tail AND a remainder record — widens balls only, ~8% at adversarial parameters); the
per-level j! factor is loose (j!/((j−i)!(N+1)^{i+1}) would suffice); and
`Mag(factorial(big(j)))` throws (fail-safe, loud) for j ≥ 21, capping iterated-integral
word depth at 21 — far above the depth-11 production need; documented here as the limit.
