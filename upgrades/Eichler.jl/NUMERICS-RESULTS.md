# Eichler.jl — Numerics results

**VERDICT: ALL FIVE GAUNTLET ITEMS PASS WITH CERTIFIED ENCLOSURES; 500+ certified digits
delivered for every object class; library ready as the pipeline's certified-arithmetic layer.**

Stack: Julia 1.12.6, Arblib 1.7 (FLINT/Arb), Nemo 0.54.
Ball discipline: every public function returns Acb/Arb enclosures; all tail bounds
(q-series majorants, iterated-integral path remainders, ODE Taylor envelopes) are proved
inequalities evaluated in Mag/Arb arithmetic — no point arithmetic anywhere. The exact
conventions are recorded in `docs/conventions.md` (equations restated from
arXiv:1704.08895, 1504.03255, 2502.00118).

## Validation gauntlet

### (a) Lattice-reduction anchor: sunrise t=0 ↔ (3/4)·L(χ₋₃,2) — PASS

The modular-representation boundary constant B⁽⁰⁾ (computed through the certified
ε-expansion of ₂F₁(−2ε,−ε;1−ε;r₃) by validated ODE transport along 1/8 → r₃, no
literature value consumed) satisfies, with certified overlap,

    B⁽⁰⁾ = (3√3/2)·L(χ₋₃,2),   i.e.   S₁₁₁(2,0) = 3·L(χ₋₃,2) = 4·(3/4)L(χ₋₃,2)·…

against Arb's independent Hurwitz-zeta evaluation L(χ₋₃,2) = 3⁻²(ζ(2,1/3)−ζ(2,2/3)).
Certified overlap at 102 / 170 / 236 / **744** decimal digits (targets 100/200/300/500;
2496-bit production run).
This is the lattice-reduction anchor (I = (3/4)L(χ₋₃,2) = S₁₁₁(2,0)/4) recovered with
rigorous error bars.

### (b) GiNaC cross-check of every modular form — PASS (≈152 digits)

Independent C++ oracle (`ginac-oracle/`, GiNaC 1.8.7 Eisenstein/modular-form kernels of
Walden–Weinzierl 2010.05271, CLN at 140 digits). Ten comparisons, all agreeing to the
oracle's full precision (~1.5×10⁻¹⁵²…10⁻¹⁵⁵ absolute):

| object | point | |Eichler − GiNaC| |
|---|---|---|
| e₁ = E₁(τ;χ̄₀,χ̄₁) | τ = 1/10 + i | 3.9e-155 |
| e₂ = E₁(2τ;χ̄₀,χ̄₁) | τ = 1/10 + i | 9.8e-155 |
| f₃(τ) | τ = 1/10 + i | 1.6e-153 |
| f₄(τ) | τ = 1/10 + i | 6.3e-154 |
| I(f₃;q), I(1,f₃;q), I(1,f₄,1,f₃;q) | q = nome(t=−1) ≈ −0.0797 (exact rational) | ≤1.9e-152 |
| same three words | q = 1/20 | ≤8.6e-153 |

GiNaC's independent Eisenstein machinery also confirms the Sebbar-polynomial identity
f₃ = 36√3(e₁³−e₁²e₂−4e₁e₂²+4e₂³) = −3√3[E₃(τ;χ̄₁,χ̄₀)−8E₃(2τ;χ̄₁,χ̄₀)] through q³⁸ exactly.
Additionally every kernel q-series is checked at runtime against its independent
eta-quotient pointwise evaluation (certified overlap; `test/runtests.jl`).

### (c) Picard–Fuchs annihilation by the printed L₀ of 2402.07311 — PASS

L₀ = ∂²ₛ + (4/s + 2/(16+s))∂ₛ + 6(6+s)/(s²(16+s)) applied to the printed period
ψ₀(s) = 32E(−s/16)/(πs^{3/2}(s+16)) through certified dK/dm, dE/dm chains:
residual balls **contain 0** at s ∈ {1/2, 3, 25/7, 50, 7/3} with radii ≤ 1.4e-151
(≤ 6e-151 across the precision sweep). Independently, the certified Frobenius transport
of (ψ₀, ψ₀′) with L₀ from s=2 to s=5 (geometric-envelope-validated Taylor steps) lands on
the independent closed-form value with certified overlap at 150 digits. The operator is
confirmed good (the appendix matrix of that paper has known typos and was never used).

### (d) AGM periods vs direct certified quadrature — PASS (5/5 points)

Parent curve y² = P(P+s)(P²+sP−4s) at s ∈ {7/3, 15/7, 101/13, 3/11, 89/2}: the certified
contour quadrature ∮dx/y (circle around the middle root pair; branch fixed by the
midpoint form u_ab with cut on [a,b] and the outer pair's cut pushed through infinity;
analyticity enforced boxwise via acb_calc) equals the Byrd–Friedman K-formula and the
explicit AGM evaluation K(m) = π/(2·AGM(1,√(1−m))):

    quad/K-formula ratio = 1 within ~1e-151 at all five points (certified overlap);
    K vs AGM overlap at all five points; at 1856 bits: 555 certified digits.

### (e) Audited amflow-cpp sunrise values through the modular representation — PASS

Equal-mass sunrise at p² = −1 (AMFlow convention = −S₁₁₁(4−2ε), m²=μ²=1), computed as:
Γ₁(6) iterated integrals to depth 11 (words {1,f₄}ᵏ,1,f₃,{f₂}ʲ of AW's all-orders sunrise result) ×
boundary constants × the exact equal-mass dimension shift (ABW), with the nome from
certified interval-Newton inversion of the Hauptmodul (branch: cusp-connected, L =
log|q|+iπ on the Euclidean axis):

    max over ε⁻²..ε⁷ of |Eichler − amflow-M240| = 2.0e-127

i.e. agreement to every honest digit of the audited anchor, with our enclosures carrying
~135 certified digits at 768 bits and **715–750 certified digits per ε-order** at
2496 bits (2.9 bits per certified digit). ε⁻², ε⁻¹ are reproduced exactly (3/2 and the closed-form tadpole value).

**Side-finding for the AMFlow-side anchors:** the `anchor_sunrise_dps130.json` (cauchy_M=120) values are
honest only to ~72 digits (its ε⁻² entry deviates from exactly 3/2 at digit 73); the
M240 file is good to its full ~128 digits. The staggered-goal protocol caught this
correctly; the same caveat is recorded in the AMFlow.cpp fork's PATCHES.md
(the sibling repository amflow-cpp).

## Beyond the gauntlet: Stage 2–3 objects pinned to the Zenodo ancillaries

Stage 2–3 validation runs against the Zenodo ancillaries of
2502.00118 (10.5281/zenodo.14733100, vendored in `vendor/zenodo-2502.00118/`):

- **Threshold period ring (s = 4m², the interior dictionary point):** our ϖ₀(4m²,m²)
  and −∂ₛϖ₀(4m²,m²) reproduce the ancillary constants a₀, a₁ to all 68 printed digits
  (|diff| ≈ 3.5e-69, 2.7e-69). The ancillary ω₀ mass-ODE and s-relation
  ∂ₛω₀ = −ω₀/s − (m²/s)∂_{m²}ω₀ are verified numerically.
- **Third-kind period G(s,t,m²):** the third-kind period definition's unspecified lower limit is resolved by the
  ancillary differential relations (`DiffRelationsEllipticFunctions.m`): we define G by
  certified quadrature of the exact one-form ∂ₜG dt (with r₉² ≡ P₄(m²−t), proven
  identically), anchored at G(s,0)=0, desingularized at t=0 (the integrand is made
  manifestly analytic at the endpoint). Validation: agreement with the printed (x₁,x₂)
  double series at exactly the x₂³ truncation order, with the correct 2³ scaling under
  x₂ → x₂/2 (rel. diff 1.7e-6 → 2.2e-7). Production: G(2,−1/2) to 555 certified digits
  in 0.32 s.
- **Abel images:** the ACKM closed form Z_{x_p} = F(arcsin u|k²)/2K(k²) (with the
  documented Feynman t+iδ deformation, δ exact dyadic) and the general-quartic
  edge-branch incomplete integral agree up to the expected 2-torsion origin shift:
  Z_quartic + Z_closed − (1+τ)/2 ∋ 0 at the deformation order (|d| = 3.4e-212 ≤ 2⁻⁶⁹⁰
  allowance). The general route works on the deformed Zγ curve
  y² = P(P+s−M)(P²+(s−M)P−4s) unchanged (certified roots → quadrature).

## Constant dictionaries (PSLQ-ready)

- `gauntlet/cusp-dictionary-500.txt`: 46 generators of the cusp ring (π powers, log 2,
  log 3, ζ(k), L(χ₋₃,k) k≤8, Re/Im Li_k at 3rd/6th roots of unity), 540-digit midpoints
  + radii; built in 0.34 s at 1856 bits.
- `gauntlet/interior-dictionary-s2-500.txt`: interior point (s,t) = (2,−1/2): ψ₀, ψ₁,
  ψ₀′, G, Z_{x_p}, π. Threshold values ϖ₀, ∂ϖ₀ at s=4m² to 500 digits in
  `results-production2.txt`.

## Cost vs digits (single core, Julia 1.12, this machine)

| object class | 100d | 200d | 300d | 500d | scaling |
|---|---|---|---|---|---|
| Γ₁(6) layer build (all kernels, N≈3·digits) | 0.05 s | 0.11 s | 0.23 s | 0.9–2.6 s | ~d² |
| boundary constants (ε¹²-deep ₂F₁ transport) | 10.1 s | 13.2 s | 18.2 s | 33–65 s | ~d^1.3 |
| sunrise S(4−2ε) through ε⁷ (depth-11 words) | 11.3 s | 20.4 s | 29.0 s | 62–96 s | ~d^1.4 |
| period, certified contour quadrature | 0.66 s | 1.08 s | 2.25 s | 4.1 s | ~d^1.3 |
| period, K/AGM closed form | <1 ms | <1 ms | <1 ms | ~1 ms | ~M(d) |
| PF Frobenius transport (s: 2→5) | 0.45 s | 0.24 s | 0.7 s | 1.1 s | ~d^1.2 |
| third-kind G (certified 1-form quadrature) | — | — | — | 0.32 s | ~d^1.3 |
| cusp dictionary (46 entries) | — | — | — | 0.34 s | ~M(d) |

Precision headroom calibration (working bits per certified output digit): kernels/periods
~3.4; sunrise pipeline ~4.2 (boundary-transport bound); rule of thumb adopted in
`gauntlet/production2.jl`: 2496 bits → ≥715 certified digits at every ε-order. Honest radii throughout;
derivative-series tail loss — the dominant term — is held negligible by the
θ=0.9 radius lemma in `qdq`.

## Rigor architecture (what makes the enclosures honest)

1. **QSeries**: ball coefficients + proved tail invariant Σ_{n>N}|aₙ|x₀ⁿ ≤ tailF
   (eta-quotient majorants via certified Euler factors; divisor-sum bounds for
   Eisenstein series; product tails by the half-mass H·F+F·H lemma).
2. **Iterated integrals**: exact monomial/log integration plus remainder records
   (A,j): |R(q)| ≤ A(|q|/x₀)^{N+1}(1+|L|)^j with the path lemma
   ∫₀¹σᴺ(X+ln 1/σ)ʲdσ ≤ Σᵢ C(j,i)X^{j−i}i!/(N+1)^{i+1}.
3. **ODE transport**: Taylor coefficients by exact recurrence in truncated-ε ℓ₁ algebra;
   tails by the verified geometric-envelope induction (bracket(m) ≤ 1 check in Mag).
4. **Quadrature**: acb_calc with boxwise analyticity certification of every branch cut
   (segment cuts via the midpoint form, ray cuts via principal factors).
5. **Branch bookkeeping**: cusp-connected nome by certified interval Newton with the
   Hauptmodul self-check; Euclidean log branch L = log|q| + iπ set exactly; realness
   tightenings by certified conjugation-isolation arguments; Feynman deformations by
   exact dyadic iδ, documented in the value's definition.

Independent adversarial audits of items 1–5 (independent adversarial review + GiNaC/Zenodo oracles)
are recorded in `audits/AUDIT.md`: 3 rigor-breaking and 6 wrong-result/contract findings
were confirmed and fixed (tail double-count, Mag-rounding gap, x0 validation, edge
orientation, deformation contract, dictionary destructuring, τ-branch guard, doc
transcription, stale G quadrature); the ₂F₁ series, the envelope transport, the
Eisenstein layer, the core QSeries algebra and the sunrise assembly survived dedicated
refutation attempts. All gauntlet numbers in this file were re-generated AFTER the fixes.

## Reproduction

```
julia --project=. -e 'using Pkg; Pkg.test()'        # full validation suite (~3 min)
julia --project=. gauntlet/production.jl             # cost curves + gauntlets c,d at 500d
julia --project=. gauntlet/production2.jl            # 500-digit deliverables + Zenodo checks
julia --project=. gauntlet/ginac_mirror.jl           # gauntlet b (needs ginac-oracle/results.txt)
ginac-oracle/build.sh                                # rebuild the independent oracle
```
