# Zenodo 10.5281/zenodo.14733100 — ancillaries of arXiv:2502.00118 — usage report

Reference data audit (`record.json` has the manifest). The files are
fetched on demand (`./fetch.sh`); the two large helicity-coefficient files
(`HelCoeffs2lbareExpTt.m`, 198 MB; `HelCoeffs2lbareExpLM.m`, 3.9 MB) are only needed
for the constants quoted verbatim below.

## What Eichler.jl uses (and verified)

1. **`DiffRelationsEllipticFunctions.m`** — exact differential relations. Key entries:
   - mass-ODE: ω₀''(s,m₂) = −4ω₀/(m₂(16m₂+s)) − (32m₂+s)/(16m₂²+m₂s)·ω₀'
   - homogeneity: ∂ₛω₀ = −ω₀/s − (m₂/s)∂_{m₂}ω₀   [verified numerically, src/puncture.jl]
   - **G one-form** (adopted as the *definition* of `third_kind_G`):
     ∂_{m²}G = s(s+2t)·ω₀·r₉/(t(s+t)−4m²s)²,
     ∂ₜG = −m²s²(16m²+s)·ω₀·r₉/(t(s+t)(t(s+t)−4m²s)²) + m²s(16m²+s)·r₉·∂_{m²}ω₀/(2t(s+t)(t(s+t)−4m²s)),
     ∂ₛG also given. This resolves the unspecified lower limit of the paper's definition of G:
     G is anchored by G(s, t→0) = 0 (consistent with the printed −√(x₁x₂³)[1+…] series,
     validated to the series' truncation order with correct x₂³ scaling).
2. **`Roots.m`** — r[9][s,t,m2] = Sqrt[t(s+t)(t(s+t)−4m²s)]; we proved r₉² ≡ P₄(m²−t)
   identically (P₄ the maximal-cut quartic), positive real on s>0, −s<t<0, t(s+t)<4m²s.
3. **`HelCoeffs2lbareExpTt.m`** — replacement-rule constants at (s,t) = (4m²,0):
   - a0 = 0.47250316546487902700221989863146743636324935057104918758306141468051
     = **ϖ₀(4m², m²)** (our normalization, exactly) — reproduced to all 68 printed digits.
   - a1 = 0.06513382235683734253462525034065980272226207906515132998094255579545
     = **−∂ₛϖ₀(4m², m²)** — reproduced to all 68 printed digits.
   (README claims 100 digits; the file carries 68.)
   - PLC1/PLC2/ELC1/ELC2: further boundary constants (not used by the test suite;
     candidates for a PSLQ closure at the threshold point).
   - Analytic continuation prescription recorded: y₁ → y₁ − i0⁺ in the physical region.
4. **`CanonicalBasis.m` / `CanonicalDiffEqs.m`** — 165 canonical masters and the three
   PDE matrices {DEs, DEt, DEm2} (ε-factorized). Downstream target: certified solution
   of this system with Eichler.jl periods/G as input transcendentals.

## Record 10.5281/zenodo.17141555 (2509.15315 crossings)

Not fetched (not needed for the current gauntlet).
