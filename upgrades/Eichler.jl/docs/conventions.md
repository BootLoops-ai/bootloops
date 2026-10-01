# Eichler.jl — mathematical conventions (equations restated from the cited papers)

All formulas below are restated from the cited papers and are the normative
conventions for this library. Sources: arXiv:1704.08895 (Adams–Weinzierl, AW),
arXiv:1504.03255 (Adams–Bogner–Weinzierl, ABW), arXiv:2502.00118 (diphoton/NPA, "BCNTW"),
arXiv:2402.07311 (ACKM).

## 1. Sunrise integral (AW, definition of the sunrise integral)

S_{ν1ν2ν3}(D, t) = (μ²)^{ν−D} ∫ ∏ d^Dk_i/(iπ^{D/2}) δ(p−k1−k2−k3) / ∏(−k_j²+m²)^{νj},
t = p². Euclidean region t < 0. AMFlow-port anchors use propagators (k²−m²+i0) ⇒
AMFlow sunrise = (−1)³ S_111 = −S_111 (three propagators), μ²=1, m²=1, t = p² = −1.

## 2. Elliptic curve of the second graph polynomial (AW sec. 5)

F = 0, roots e1,e2,e3 (AW, definition of the roots) with D̃ = (t−m²)³(t−9m²),
k² = (e3−e2)/(e1−e2). Periods (AW def_periods):
  ψ1 = 4μ²/D̃^{1/4} K(k),  ψ2 = 4iμ²/D̃^{1/4} K(k'),
  φ1 = 4μ²/D̃^{1/4} (K(k)−E(k)),  φ2 = 4iμ²/D̃^{1/4} E(k')
K(x), E(x) = complete elliptic integrals, **modulus convention** K(x)=∫₀¹dt/√((1−t²)(1−x²t²))
(Arb's acb_elliptic_k(m) uses m = x², i.e. K_AW(k) = elliptic_k(k²)).
Wronskian: W = ψ1 ψ2' − ψ2 ψ1' = −12πiμ⁴/(t(t−m²)(t−9m²)).
τ = ψ2/ψ1, q2 = e^{iπτ}. L2^{(0)} ψ_i = 0 with
L2^{(0)} = d²/dt² + (1/t + 1/(t−m²) + 1/(t−9m²)) d/dt
         + (1/m²)(−1/(3t) + 1/(4(t−m²)) + 1/(12(t−9m²))).

## 3. Maximal-cut curve and Γ₁(6) (AW, section on Γ₁(6)) — PROGRAM STANDARD

Maximal-cut quartic E_C: w² = (z−t/μ²)(z−(t−4m²)/μ²)(z² + 2m²z/μ² + (m⁴−4m²t)/μ⁴),
roots z1 = (t−4m²)/μ², z2 = (−m²−2m√t)/μ², z3 = (−m²+2m√t)/μ², z4 = t/μ²
(branch: t → t+iδ, √t cut along negative axis).
k_C² = (z3−z2)(z4−z1)/((z3−z1)(z4−z2)), k_C'² = (z2−z1)(z4−z3)/((z3−z1)(z4−z2)).
ψ_{1,C} = 4μ² K(k_C)/((m+√t)^{3/2}(3m−√t)^{1/2}),
ψ_{2,C} = 4iμ² K(k_C')/((m+√t)^{3/2}(3m−√t)^{1/2}).
τ_C = ψ_{2,C}/ψ_{1,C}, q_C = e^{2πiτ_C}. Near t=0: ψ_{1,C}=ψ1, 2ψ_{2,C}=ψ2+ψ1,
τ_C = (τ+1)/2, q_C = q_{2,C}² = −q2. W_C = −6πiμ⁴/(t(t−m²)(t−9m²)) = W/2.

Hauptmodul (AW eq. hauptmodul):
  t = 9m² η(6τ_C)⁸ η(τ_C)⁴ / (η(2τ_C)⁸ η(3τ_C)⁴)     [q-expansion 9m²(q_C − 4q_C² + ...)]

Sebbar generators of M_*(Γ₁(6)) (polynomial ring):
  e1 = E1(τ_C; χ̄0, χ̄1),  e2 = E1(2τ_C; χ̄0, χ̄1),
χ̄0 = trivial, χ̄1 = Kronecker (−3/n). Generalised Eisenstein series (AW app. B):
  E_k(τ; φ, ψ) = a0 + Σ_{m≥1} (Σ_{d|m} ψ(d) φ(m/d) d^{k−1}) q_M^m,  q_M = e^{2πiτ/M},
  a0 = −B_{k,ψ}/(2k) if cond(φ)=1 else 0.  B_{1,χ̄1} = −1/3 ⇒ a0(E1(·;χ̄0,χ̄1)) = 1/6.
  NOTE: the paper's printed q_M = e^{2πiτ/M} is inconsistent with the printed
  Γ₁(6) identities for e1 = E1(τ;χ̄0,χ̄1) (M = 3 would give q^{1/3} powers); the
  operative convention — used by the implementation, by GiNaC 1.8, and required by
  every η-quotient identity — is integer q-powers, q = e^{2πiτ}.
  B_{2,K}(τ) = E2(τ) − K E2(Kτ), E2 = −1/24 + Σ σ1(m) q^m.

Kernel forms (AW, section on Γ₁(6); all in τ_C, μ²=m²=1 internally):
  ψ1/π     = 2√3 (μ²/m²)(e1+e2)  = (2μ²/√3m²) η(3τ_C)η(2τ_C)⁶/(η(τ_C)³η(6τ_C)²)
  f_{1,C}  = 3√2 e1              = (t+3m²)/(2√6 μ²) · ψ1/π
  f_{2,C}  = −6(e1² + 6e1e2 − 4e2²)
           = (1/2iπ)(ψ1C²/W_C)(3t²−10m²t−9m⁴)/(2t(t−m²)(t−9m²)) = −10B_{2,2}+4B_{2,3}−2B_{2,6}
  f_{3,C}  = 36√3 (e1³ − e1²e2 − 4e1e2² + 4e2³)
           = (μ²ψ1C³/(4πW_C²))·6/(t(t−m²)(t−9m²)) = −3√3 η(τ_C)⁵η(3τ_C)η(6τ_C)⁴/η(2τ_C)⁴
  f_{4,C}  = f_{1,C}⁴
  g_{2,0,C} = −12(e1²−4e2²) = (1/2iπ)(ψ1C²/W_C)/t = η(τ)⁴η(3τ)⁴/(η(2τ)²η(6τ)²)
  g_{2,1,C} = −18(e1²+e1e2−2e2²) = (1/2iπ)(ψ1C²/W_C)/(t−m²) = −9η(3τ)³η(6τ)³/(η(τ)η(2τ))
  g_{2,9,C} = 6(e1²−3e1e2+2e2²) = (1/2iπ)(ψ1C²/W_C)/(t−9m²) = −η(τ)⁷η(6τ)⁷/(η(2τ)⁵η(3τ)⁵)
  g_{3,0,C} = −2 f_{3,C}
  g_{3,1,C} = −108√3 (e1³ − 3e1e2² + 2e2³) = −54√3 η(6τ_C)⁹/η(2τ_C)³
Cusp values (q_C=0): ψ1/π → 2μ²/(√3m²), f1 → √2/2, f2 → −1/2, f3 → 0, g_{2,0} → 1,
others → 0. (f2 cusp value: −6(1/36+6/36−4/36)=−1/2 ✓.)

## 4. Iterated integrals (AW sect. 2.4)

I(f1,...,fn; q) = ∫₀^q dq1/q1 ... with tangential base point at q=0: constant terms a0
integrate to a0·ln q (operator R removes ln(q0) terms). Single non-cusp form:
I(f;q) = a0 ln q + Σ_{j≥1} (a_j/j) q^j. Shuffle algebra holds. For q = q_C complex,
ln q_C := 2πi τ_C (NOT principal branch).

## 5. All-orders sunrise (AW, main all-orders result, Γ₁(6) version)

S_111(2−2ε, t) = (ψ_{1,C}/π) · exp(−ε I(f_{2,C};q_C) − 2εL − 2γ_E ε + 2Σ_{n≥2}(−1)ⁿζ_n εⁿ/n) ×
 { [Σ_j (ε^{2j} I({1,f_{4,C}}^j; q_C) − ½ ε^{2j+1} I({1,f_{4,C}}^j,1; q_C))] · Σ_k ε^k B^{(k)}(2,0)
   + Σ_j ε^j Σ_{k≤⌊j/2⌋} I({1,f_{4,C}}^k, 1, f_{3,C}, {f_{2,C}}^{j−2k}; q_C) },
L = ln(m²/μ²). Boundary values B^{(k)}: S_111(2−2ε,0) = (ψ1(0)/π) Γ(1+ε)² e^{−2εL} Σ ε^k B^{(k)},
where the printed AW boundary equation gives the γ-STRIPPED series (AW define
S = e^{−2γε} Σ εʲ S^(j), so the actual integral is e^{−2γε} times the RHS below):
  Σ_j ε^j S^(j)(2,0) = e^{2γε} Γ(1+2ε) (m²√3/μ²)^{−1−2ε} [3/(2ε²)·Γ(1+ε)²/Γ(1+2ε)·h − π/ε],
h = (1/i)[(−r3)^{−ε} ₂F₁(−2ε,−ε;1−ε;r3) − (−r3^{−1})^{−ε} ₂F₁(−2ε,−ε;1−ε;r3^{−1})],
r3 = e^{2πi/3}, arg(−r3) = −π/3, arg(−r3^{−1}) = +π/3.

## 6. Dimension shift, equal mass (ABW, dimension-shift relation)

S_111(4−2ε,t) = 1/(6(1−2ε)(1−3ε)(2−3ε)) ×
 { (t+3m²)(t−m²)(t−9m²)/μ⁴ · d/dt S_111(2−2ε,t)
   + [ (t−m²)(t−9m²)/μ⁴ + ε (t²+22m²t−87m⁴)/μ⁴ ] S_111(2−2ε,t)
   + 3(1−ε)² [ −6μ²/m² + ε μ²(t+21m²)/m⁴ ] [T1(4−2ε)]² },
T1(D) = tadpole = Γ(1−D/2)(m²/μ²)^{D/2−1} (check: ε⁻² coefficient −3m²/2μ² ✓).
d/dt via dq_C/dt: dt = ψ_{1,C}²/(2πi W_C) dq_C/q_C  [AW: dt = ψ1²/(iπW) dq2/q2, W_C=W/2].

## 7. 2402.07311 (ACKM) parent curve, Picard–Fuchs, periods (m_t²=1)

Curve y² = P(P+s)(P²+sP−4s) [mt²=1], roots r1=(−s−√s√(s+16))/2, r2=−s, r3=0,
r4=(√s√(s+16)−s)/2; marked point x_p = −s−t.
k² = 2√s√(s+16)/(s+√s√(s+16)+8),  u_{xp} = √(1+(t−8)√s/(t√(s+16)))/√2,
Z_{xp} = F(arcsin(u_{xp})|k²)/(2K(k²))  [Mathematica convention F(φ|m), K(m), m=k²].
Printed PF operator (gauntlet (c); operator verified good, appendix matrix NOT):
  L0 = ∂²_s + (4/s + 2/(16+s)) ∂_s + 6(6+s)/(s²(16+s))
Annihilates ψ0 = 32 E(−s/16)/(π s^{3/2}(s+16)) and
ψ1 = 32(E(s/16+1) − K(s/16+1))/(s^{3/2}(s+16))  [K(m), E(m) with m-argument convention].

## 8. 2502.00118 (BCNTW) curve, Frobenius periods, third-kind period G

Curve Y² = P4(X) = (m²−X)(m²+s−X)(m⁴−3m²s−X(2m²+s)+X²); X = m²−P maps to ACKM curve.
Singular points s ∈ {0, −16m², ∞}. Frobenius periods at s=0:
  ϖ0^[0] = (1/√(sm²))(1 − s/64m² + 9s²/16384m⁴ − 25s³/1048576m⁶ + ...)
  ϖ1^[0] = ϖ0 log(s/m²) + (1/√(sm²))(−s/32m² + 21s²/16384m⁴ − 185s³/3145728m⁶ + ...)
at s=∞: ϖ0^[∞] = (1/s)(1 − 4m²/s + 36m⁴/s² − 400m⁶/s³ + ...),
  ϖ1^[∞] = ϖ0^[∞] log(m²/s) + (1/s)(−8m²/s + 84m⁴/s² − 2960m⁶/(3s³) + ...)
at s=−16m² (v = s+16m²): ϖ0^[−16] = (1/m²)(1 + 3v/64m² + 41v²/16384m⁴ + 147v³/1048576m⁶+...),
  ϖ1^[−16] = ϖ0^[−16] log(v/m²) + (1/m²)(v/32m² + 37v²/16384m⁴ + 455v³/3145728m⁶ + ...).
Third-kind period (definition of G):
  G(s,t,m²) = ∫^{m²} dx s(s+2t)√(P4(x−t)) / (t(s+t)−4sx)² · ϖ0(s,x)
(ϖ0(s,x) = holomorphic period with mass argument x). Series check, x1 = −tm²/s², x2 = s/4m²:
  G = −√(x1 x2³)[1 − x2/32 − 10x1x2 + 3x2²/1024 − 43x1x2²/16 + 14x1²x2² + O(x³)].
Validation target: Zenodo ancillaries 10.5281/zenodo.14733100 (2502.00118) and
10.5281/zenodo.17141555 (2509.15315), NOT the printed 2402.07311 appendix.

## 9. Zγ deformed curve (m_t²=1, M = m_Z²)

y² = P(P+s−M)(P²+(s−M)P−4s); marked point P = M−s−t; five singular fibers:
s=0 (I4), s=M (I2), (s−M)²+16s=0 (I1,I1), s=∞ (I4).

## 10. Boundary-constant rings

Cusps: Dirichlet L-values of conductor | 6 (+ π, logs): L(χ−3, k). The "t=0 sunrise
constant = (3/4)L(χ−3,2)" anchor refers to the Bessel moment I = ∫₀^∞ x K₀(x)³ dx = S₁₁₁(2,0)/4;
the AW boundary value itself is S₁₁₁(2,0) = 3L(χ−3,2), i.e. B^(0) = (3√3/2)L(χ−3,2). L(χ,s) = q^{−s} Σ_{a=1}^{q} χ(a) ζ(s, a/q) via certified Hurwitz zeta.
Interior points: period/quasi-period/G rings (AGM/Frobenius-computable).
