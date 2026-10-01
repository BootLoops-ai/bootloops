# PROVEN_MAJORANT — certified tail envelope for Frobenius log-tower series

Status: proved; independently re-verified (Section 7). Implementation:
`envelope_certified.py` (this directory); battery: `selftest_envelope.py`.
Replaces the EMPIRICAL last-window ×4 tail heuristic of
`upgrades/Eichler.jl/src/cy_transport.jl` (`mum_frobenius_basis`), which is
unsound: two exact counterexamples, see §5.

## 1. Setting

L Fuchsian at z = 0, order r, θ-form  z^{-smin} L = Σ_{s=0}^{S} z^s R_s(θ),
R_s ∈ ℚ[θ], deg R_s ≤ r, R_0 = indicial polynomial with deg R_0 = r (regular
singular point), leading coefficient lc, and factorization over ℂ
R_0(θ) = lc·Π_i (θ − λ_i)^{μ_i}, Σ_i μ_i = r.

Jet space V = ℂ^J with the ℓ1 norm ‖v‖ = Σ_j |v_j|. Let D : V → V be ANY nilpotent
operator with D^J = 0 and induced norm ‖D‖ ≤ 1. The two instances used downstream:
  (i)  multiplication by ρ on ℂ[ρ]/(ρ^J)   (cy_transport.jl ρ-jets),
  (ii) the divided-power shift (Dγ)_j = γ_{j+1} (the DKMM certified-W0 log-layer towers).
Both have induced ℓ1 norm exactly 1.

DEFINITION (jet tower past horizon N0, exponent e ∈ ℝ). A sequence b_m ∈ V
(m ≥ 0) such that for every m > N0:
    R_0(e+m+D)·b_m = − Σ_{s=1}^{min(S,m)} R_s(e+m−s+D)·b_{m−s}.          (REC)
The MUM Frobenius jets satisfy (REC) with e = α, J = r, N0 = 0, D = (i);
Frobenius towers of a resonant operator satisfy it past the largest resonance;
the DKMM order-6 towers satisfy it with D = (ii), e = 0, N0 = 3.

## 2. Lemma (certified geometric tail envelope)

Let (b_m) be a jet tower past horizon N0 at exponent e. Put
  β  = max( 0,  max_i Re λ_i − e ),
  rbar_{s,k} = |[θ^k] R_s|,   c_s = max(0, e − s + 1),
  G(m)   = Π_i Σ_{j=0}^{J−1} C(μ_i+j−1, j) (m−β)^{−j},
  Φ(m,t) = G(m) / ( |lc| (m−β)^r ) · Σ_{s=1}^{S} t^{−s} Σ_{k=0}^{r} rbar_{s,k} (m+c_s)^k.

Choose N and t > 0 with
  (H1) N ≥ N0,  N > β,  N ≥ S + max(0, ⌈−e⌉);
  (H2) Φ(N+1, t) ≤ 1.
Set K = max_{N−S < m′ ≤ N} ‖b_{m′}‖ t^{−m′}. Then

  (C1)  ‖b_m‖ ≤ K t^m   for ALL m > N;
  (C2)  for 0 < x < 1/t and each jet level j:
        | Σ_{m>N} (b_m)_j x^m |  ≤  Σ_{m>N} ‖b_m‖ x^m  ≤  K (tx)^{N+1} / (1 − tx);
  (C3)  (derivative weights) for any d ≥ 0, if q := t·x·e^{d/(N+1)} < 1 then
        Σ_{m>N} ‖b_m‖ m^d x^m  ≤  K (N+1)^d (tx)^{N+1} / (1 − q).

All quantities in (H2), K, (C2), (C3) are exactly computable in ℚ when the
coefficients, e, t, x are rational and the λ_i are rational (else replace
Re λ_i by any certified upper bound; the lemma holds a fortiori).

## 3. Proof

Throughout ‖·‖ is the ℓ1 vector norm and its induced operator norm;
submultiplicativity and the triangle inequality are used silently. m > N.

(P1) Polynomial evaluation bound. For u ∈ ℝ and P ∈ ℝ[θ]:
‖P(u+D)‖ ≤ Σ_k |P_k| ‖(u+D)^k‖ ≤ Σ_k |P_k| (|u| + ‖D‖)^k ≤ Σ_k |P_k| (|u|+1)^k.

(P2) Numerator bound. By (H1), e+m−s ≥ e+N+1−S > ... ≥ 0 for 1 ≤ s ≤ S (if e ≥ 0
this is m−s > N−S ≥ 0; if e < 0 it is m ≥ N+1 ≥ S+⌈−e⌉+1 > s−e). Hence
|e+m−s| + 1 = m + (e−s+1) ≤ m + c_s, and by (P1)
  ‖R_s(e+m−s+D)‖ ≤ Σ_k rbar_{s,k} (m + c_s)^k.

(P3) Inverse bound. For each root, a_i := e+m−λ_i has Re a_i = m − (Re λ_i − e)
≥ m − β > 0, so |a_i| ≥ m − β > 0 and a_i I + D is invertible with
  (a_i I + D)^{−1} = Σ_{j=0}^{J−1} (−1)^j D^j a_i^{−j−1}
(check: multiply by a_i I + D; the sum telescopes and D^J = 0). Taking the μ_i-th
power and collecting D^j (negative-binomial coefficients, truncated exactly by
D^J = 0):
  (a_i I + D)^{−μ_i} = Σ_{j=0}^{J−1} C(μ_i+j−1, j) (−D)^j a_i^{−μ_i−j},
  ‖(a_i I + D)^{−μ_i}‖ ≤ Σ_{j=0}^{J−1} C(μ_i+j−1, j) (m−β)^{−μ_i−j}.
Multiplying over i (Σ μ_i = r):
  ‖R_0(e+m+D)^{−1}‖ ≤ G(m) / ( |lc| (m−β)^r ).

(P4) Monotonicity of Φ(·, t) on (β, ∞). Each factor (m−β)^{−j} in G is
nonincreasing. Each summand carries (m+c_s)^k / (m−β)^r with 0 ≤ k ≤ r,
c_s ≥ 0 ≥ −β; since m + c_s ≥ m ≥ m − β > 0,
  d/dm log[ (m+c_s)^k (m−β)^{−r} ] = k/(m+c_s) − r/(m−β) ≤ (k−r)/(m−β) ≤ 0.
So Φ(m, t) ≤ Φ(N+1, t) ≤ 1 for all m ≥ N+1 by (H2).

(P5) Induction for (C1). Claim ‖b_m‖ ≤ K t^m for m > N−S (strong induction).
Base: m ∈ (N−S, N] is the definition of K. Step m > N: (REC) applies since
m > N ≥ N0 and m > S (so min(S,m) = S); the used indices m−s ∈ [m−S, m−1] all
lie in (N−S, m), covered by the hypothesis. R_0(e+m+D) is invertible by (P3), so
  ‖b_m‖ ≤ ‖R_0(e+m+D)^{−1}‖ Σ_{s=1}^S ‖R_s(e+m−s+D)‖ ‖b_{m−s}‖
        ≤ [G(m)/(|lc|(m−β)^r)] Σ_s (Σ_k rbar_{s,k}(m+c_s)^k) · K t^{m−s}
        = K t^m · Φ(m, t) ≤ K t^m.                                       ∎(C1)

(P6) (C2): |(b_m)_j| ≤ ‖b_m‖ (one ℓ1 component), and Σ_{m>N} K (tx)^m is the
geometric series with ratio tx < 1.                                     ∎(C2)

(P7) (C3): the terms u_m = ‖b_m‖ m^d x^m ≤ K m^d (tx)^m have
  K(m+1)^d (tx)^{m+1} / (K m^d (tx)^m) = tx (1+1/m)^d ≤ tx e^{d/m} ≤ q < 1
for m ≥ N+1, so the sum is ≤ first-term/(1−q) = K (N+1)^d (tx)^{N+1}/(1−q). ∎(C3)

Remark (zero seed window). If all S seeds in (N−S, N] vanish then K = 0 and (C1)
asserts b_m = 0 for m > N — correct: (REC) is an S-step determination with
invertible leading operator past the horizon, so S consecutive zero jets force
all later jets to vanish identically.

## 4. Honest scope — what is NOT claimed (obstruction for coefficientwise bounds)

The ℓ1 route is blind to sign cancellation among θ-form coefficients. The
smallest admissible t (root of Φ(N+1, t) = 1) converges, as N → ∞, to the root
t*∞ of Σ_s rbar_{s,r} t^{−s} = |lc|, and 1/t*∞ can be strictly smaller than the
true convergence radius R. Measured on the reference operators (exact arithmetic):
  K3 banana   (R = 1/16, 1/x0 = 32): t_min = 24.02 (N=80), 23.29 (N=200); t*∞ ≈ 22.81
  CY3 banana  (R = 1/25, 1/x0 = 50): t_min = 45.66 (N=80), 43.05 (N=200)
An impossibility argument (route C, §7) shows this gap is INTRINSIC to any bound built
from coefficientwise absolute values: for K3, Σ_s rbar_{s,r} w^s at w = R = 1/16
equals 20/16 + 64/256 = 1.5 > 1, so no coefficientwise-monotone argument can
certify radius R. The literal target "proved envelope on the whole disk avoiding
singularities" is therefore UNPROVABLE in this class; the lemma above is the
explicit-constant repair: a certified envelope on |z| < 1/t_min(N), with
1/t_min(N) ↑ 1/t*∞ and x0 = R/2 admissible on every reference case (§6 table).
Consumers wanting x0 closer to R must either raise N (helps until t*∞) or switch
to a cancellation-aware method (certified transport re-anchoring, already
rigorous in frobenius.jl).

## 5. Counterexample search (exact arithmetic; permanent regressions in
##    selftest_envelope.py)

Against the lemma (m tested to 240-400, all Fractions):
  Legendre-type (t=1.1, 1.9), sunrise, K3, CY3, cancellation op, forced-log
  resonant tower (beta=2), negative shift e=-2 with synthetic seeds: no
  violation found.
  sup_m ||b_m||/(K t^m) in [0.53, 0.96] — scales at the claimed rate; the
  envelope is tight to within a factor ~1.04-2 on coefficient norms.
Against the EMPIRICAL x4 scheme of cy_transport.jl — FALSIFIED twice:
  D1 lacunary: h = 1/(1-z^2), N=81 (odd), x=1/2: c_N = 0 so the heuristic
     returns Mag(0); true tail 2.8e-25 > 0.  UNSOUND.
  D2 oscillatory: nearest singularities a complex pair ((3+-4i)/5, |.|=1,
     cos-modulated coefficients): last-window ratio > 0.95/x at EVERY window
     tested (180/180), so the heuristic returns Mag(0); true tails ~1e-3..1e-9
     > 0.  UNSOUND.  (The real-singularity reference cases do not trigger
     these branches, which is why the x4 legs *happened* to work — the DKMM
     certified-W0 towers inherit no error from this, but the flag was
     correct: it was never a bound.)

## 6. Comparison table — cost of rigor (generated by selftest_envelope.py E6,
##    written to envelope_comparison.md)

measured = exact |c_{j,m}| x^m summed to the stored depth + certified remainder.
  K3 banana (r=3, MUM, R=1/16, N=128, x=R/2): proven/measured 9.1-25.2 per log
    level; empirical/measured = 4.00 (the x4, sound here by luck).
  Resonant forced-log op (r=2, beta=2, N=64, x=R/2): proven/measured 5.9 (j=0),
    1.55 (j=1) — the proven bound BEATS the empirical 4.02 on the log level.
  DKMM order-6 L_s, reference tower slice (P,Q)=(0,0) (J=4, beta=3, N=240,
    x=1/112 ~ R/2, R ~ 1/55.6): proven tail 5.5e-249 absolute — sound but
    ~1e77 x looser than this slice's measured 8e-327..2e-325: the envelope
    prices the OPERATOR's worst-case growth (t_min = 80.9) while this
    particular solution slice grows like ~5^m. Certificates remain far below
    any downstream gate tolerance; cost quantified, not hidden.
  CY3 banana (e=1): admissible at x=R/2 (t_min=43.05 at N=200 < 50); in-table
    runs omitted (minutes-budget), numbers in section 4.

## 7. Verification

The lemma was checked independently of the implementation, as follows.
  1 Counterexample search (section 5): exact-arithmetic attempts to violate
    (C1) on Legendre-type, sunrise, K3, CY3, cancellation, forced-log
    resonant and negative-shift towers found no violation; the same search
    falsified the empirical x4 scheme twice.
  2 Step-by-step check of the proof (sections 2-3), which surfaced two
    hypotheses now stated explicitly:
    (i) (P2) needs N >= S + ceil(-e) (included in (H1));
    (ii) deg R_0 = r must be checked by implementations (stated in section 1).
  3 Planted-error test of the battery: an end-to-end numeric acceptance
    test FAILED to catch a deliberately dropped G(m) factor (the cases
    have slack) — the bar was tightened to STEP-LEVEL exact induced-norm
    gates (selftest E1); the planted error is then caught 12/12 on theta^r,
    where the (P3) bound is exactly attained (equality) — tight.
  4 Alternative routes: route B (a scalar majorant sequence B_m) dominates
    and agrees; route C (impossibility): no coefficientwise-absolute-value
    bound can certify the true radius R when sign cancellation is present —
    exact obstruction sum_s rbar_{s,r} R^s = 3/2 > 1 for K3; the literal
    full-disk envelope is UNPROVABLE in this class; this lemma is the
    explicit-constant repair (certified subdisk 1/t_min, finite-N sharpened).
  5 Independent re-derivation and numeric check: (P3) re-derived by direct
    matrix inversion (exact identity), (P4) by consecutive-value comparison
    m=4..200; the E2 test instance recomputed exactly and MATCHING the E2
    gate digit-for-digit: Phi(49,11/10)=0.965605317512261,
    K=2.551203884e-4, sup ratio=0.890732579021090; no circular imports;
    112 induction test points checked without error.
  6 Interface trace: section 8.

## 8. How a consumer applies the lemma

  1 Frobenius-basis truncation (a MUM-point series constructor such as
    mum_frobenius_basis in cy_transport.jl): take t from the Phi bisection
    (exact rationals fit an exact-Q coefficient path), K from the last r+S
    window, and use the tail bound Mag = K(tx)^{N+1}/(1-tx) per solution h_j.
    This is valid provided x0 < 1/t_min; at x0frac = 1/2 all four reference
    banana/K3/CY3 operators are admissible (sections 4, 6). When (H2) fails,
    shrink x0frac (x0 <= 1/(2 t_min)) or raise nterms; when no admissible t
    exists the bound does not apply and a constructor must raise an error
    rather than return a zero radius.
  2 Derivative tails: corollary (C3) with d <= r-1 replaces an N^d ratio
    factor. It needs q = t x e^{d/(N+1)} < 1, which holds at x0/2-type
    interior points for the reference cases; otherwise raise N.
  3 The DKMM certified-W0 towers (reference data, not included in the package;
    TERRIER_KKLT_BANK): the recurrence is verified in the divided-power picture
    directly on the stored JSON (E5a), and the envelope certifies the m > 240
    extension tail.
  4 Transport legs that used a stepped (empirical) majorant can be upgraded to
    the proven envelope by passing their tail gates through certify(); where a
    leg already passed with margin no numerical change results, since the
    proven bound is larger than the empirical one but within that margin.

Inputs: the restricted operator theta-form (operator_LS.json) and the DKMM
certified-W0 extended towers of the reference data, read-only; banana
theta-forms transcribed from Eichler.jl banana_pf; no duplication of
tools/annihilator / tools/pf_rank.jl / tools/wayfinder functionality (this
module bounds tails of towers others produce; it never derives operators or
towers itself).
