# B7 — box-floor Cauchy method: certified period floors over a polydisc tile

Method doc for the Cauchy-floor layer of the direct-summation route
(`b7_checks.py` verifies the inequalities below numerically — the shipped
battery requires all 4 PASS lines; `b7_numbers.py` prints the closed-form
majorant/tail tables). The derivation is self-contained: everything is an
explicit expression in the summation-layer data (center balls, truncation
order, entropy scales).

LAW: derivation route, fail-closed. A hypothesis that cannot be verified in
ball arithmetic yields NO bound (UNDECIDED => FAIL); a failed tile never
admits and never contributes a floor.

## 1. Setting + data interface

Tile T = prod_I clD(c_I, h_I) in C^n (closed polydisc), c_I != 0; optional
fattening radii r_I > h_I; sig_I = |c_I| + r_I. Fundamental period
Pi^0 = sum_k multinom(|k|;k)^2 phi^k — ALL coefficients >= 0. Jet components
a_alpha(phi) = theta^alpha Pi^0(phi), depth d = |alpha| bounded (<= 4 here).
Data consumed per tile center c: certified balls (m_alpha, r_alpha) for
a_alpha(c) (truncation at total order M + tail + rounding, outward-rounded);
the order M; ball values of the entropy scales S_h = sum_I sqrt(|c_I| + h_I)
and S_r = sum_I sqrt(sig_I).
HYPOTHESES (each checked per tile in ball arithmetic):
  H1: h_I < |c_I| for all I (tile clears the coordinate hyperplanes).
  H2: S_h < 1 (tail route); additionally S_r < 1 for the Cauchy route L2.
  H3: (closed-form tail corollary only) x*e^(d/(M+1)) < 1 at x = S_h^2.

## 2. L0/L1 — weighted entropy majorant + tails

L0 (majorant with jet weights). For |alpha| = d <= 5 and any phi with
S := sum_I sqrt|phi_I| < 1:
  |theta^alpha Pi^0(phi)| <= sum_{m>=1} m^d S^{2m} (+1 iff d=0) = Li_{-d}(x),
x = S^2.
PROOF. Termwise, theta^alpha multiplies the k-term by prod_I k_I^{alpha_I}
<= m^d (each k_I <= m = |k|). Positivity + sum-of-squares <= square-of-sum on
the multinomials: sum_{|k|=m} multinom(m;k)^2 prod|phi_I|^{k_I} <=
(sum_{|k|=m} multinom(m;k) prod sqrt|phi_I|^{k_I})^2 = S^{2m}. Sum over m.
QED (numeric spot-check 20/20 PASS, b7_checks.py check 1).
Closed forms (exact; x = S^2): Li_0 = 1/(1-x); Li_{-1} = x/(1-x)^2; Li_{-2} =
x(1+x)/(1-x)^3; Li_{-3} = x(1+4x+x^2)/(1-x)^4; Li_{-4} =
x(1+11x+11x^2+x^3)/(1-x)^5; Li_{-5} = x(1+26x+66x^2+26x^3+x^4)/(1-x)^6
(Eulerian polynomials; spot-check vs direct sums PASS, b7_checks.py check 2).
CORNER-MONOTONE: the majorant depends only on (|phi_1|,..,|phi_n|),
increasing in each; its sup over ANY polydisc is its value at the
positive-real corner — one evaluation, no n-dim optimization. This holds for
the MAJORANT only; the jets themselves are complex-valued and are never
evaluated at a "corner".
L1 (tails). tail_d(M, x) := sum_{m>M} m^d x^m. Two lawful evaluations:
 (i) EXACT: Li_{-d}(x) minus the degree-M partial sum, both in ball
     arithmetic with outward rounding;
 (ii) COROLLARY (closed form, needs H3): tail_d(M,x) <= (M+1)^d x^{M+1} /
     (1 - x e^{d/(M+1)}). Proof: m^d = (M+1)^d e^{d ln(m/(M+1))} <=
     (M+1)^d e^{d(m-M-1)/(M+1)} (ln y <= y-1), then geometric.
     Spot-check PASS (b7_checks.py check 3).

## 3. L2-L5 — box oscillation + certified period floors

L2 (Cauchy route). Under H2 (S_r < 1) each a_alpha is holomorphic on the
fattened polydisc P_r = prod_I D(c_I, r_I). For z in T the coordinate disc
D(z_I, r_I - h_I) stays inside P_r, so the one-variable Cauchy estimate gives
sup_T |d/dphi_I a_alpha| <= Li_{-d}(S_r^2)/(r_I - h_I). Integrating along the
straight segment c -> z (T convex):
  osc_T(a_alpha) := sup_{z in T} |a_alpha(z) - a_alpha(c)|
                 <= [sum_I h_I/(r_I - h_I)] * Li_{-d}(S_r^2)  =: W^C_alpha(T).
L3 (theta route, alternative). d/dphi_I a_alpha = a_{alpha+e_I}/phi_I, a
depth-(d+1) jet, so under H1+H2: osc_T(a_alpha) <=
[sum_I h_I/(|c_I| - h_I)] * Li_{-(d+1)}(S_h^2) =: W^T_alpha(T). No fattening
needed. Lawful bound: min(W^C, W^T).
L4 (split-tail workhorse). Split a_alpha = P_M + R_M (truncation + tail).
osc_T(P_M) by direct interval/ball evaluation of the FINITE sum over the tile
box (tight: scales with h * true derivative, not with the majorant);
box-uniformly |R_M| <= tail_d(M, S_h^2) (corner-monotone), hence
osc_T(R_M) <= 2 tail_d(M, S_h^2).
  omega_alpha(T) := min{ W^C_alpha, W^T_alpha,
                         osc-interval(P_M) + 2 tail_d(M, S_h^2) }.
L5 (tile-uniform enclosures + PERIOD FLOORS). For every phi in T:
  a_alpha(phi) in X_alpha(T) := ball(m_alpha, r_alpha + omega_alpha(T)), and
  min_{phi in T} |a_alpha(phi)| >= |m_alpha| - r_alpha - omega_alpha(T),
a certified LOWER bound on |period components| over the box iff the RHS is
ball-positive; RHS <= 0 or UNDECIDED => NO floor (fail-closed).
PROOF of L5: triangle inequality on a_alpha(phi) = a_alpha(c) +
(a_alpha(phi) - a_alpha(c)) with |a_alpha(c) - m_alpha| <= r_alpha and
L2/L3/L4. QED.

## 4. L6-L7 — transport factor + Gram floor

L6 (exact transport polynomial over the tile). When the monodromy frame
T_1..T_n is pairwise commuting and index-2 unipotent ((T_K - 1)^2 = 0),
N_K := log T_K = T_K - 1 EXACTLY (the log series terminates), N_K^2 = 0, all
N_K commute. Hence
  E(t) = exp(sum_I t_I N_I) = prod_I (1 + t_I N_I)
— exact multilinear, NO truncation. Tile t-box: t_I = log(phi_I)/(2 pi i);
Delta t_I = log(phi_I/c_I)/(2 pi i), principal branch lawful under H1
(|phi_I/c_I - 1| <= h_I/|c_I| < 1, no branch crossing), and
  |Delta t_I| <= -log(1 - h_I/|c_I|) / (2 pi)
(|log(1+w)| <= -log(1-|w|); spot-check 2000/2000 PASS, b7_checks.py check 4).
E over T = ball evaluation of the exact polynomial on the Delta-t box.
GATE (fail-closed): before ANY use of a truncated-degree exp form, verify the
index-2 identity (T_K - 1)^2 = 0 for all K on the frame actually in use; a
gate failure => HALT.
L7 (ball-matrix eigenvalue floor; Weyl). If H(T) = Hhat +- R is a Hermitian
ball matrix containing H(phi) for EVERY phi in T, then
  min_{phi in T} lam_min H(phi) >= lam_min(Hhat) - ||R||_F .
PROOF: Weyl, |lam_min(A+E) - lam_min(A)| <= ||E||_2 <= ||E||_F, applied at
E = H(phi) - Hhat, ||E||_F <= ||R||_F entrywise. QED.
BOUND FORM: per tile T, enclose the frame vector and its flag jets by
Pi(T) = E(Delta-t box) * X(T) (ball matrix-vector, L5+L6); run the Gram
construction end-to-end in ball arithmetic on these tile-uniform enclosures
to get H(T); then
  floor >= min_{T in tiling} [ lam_min(Hhat(T)) - ||R_H(T)||_F ]
with every ingredient computable from (m_alpha, r_alpha, M, S_h, S_r, h, r,
frame N_K) plus exact linear algebra. Fail-closed: a non-positive floor or
any UNDECIDED comparison => tile FAIL => h-halving fallback; a failed tile
never admits.

## 5. Practical sizing notes

- The pure-majorant routes W^C/W^T are the LEMMA (and become tight as
  h -> 0); at practical tile widths the WORKHORSE is L4 (interval truncation
  + tail), whose width scales with h times the TRUE derivative.
- Depth-d tails grow steeply with the entropy scale: choose the truncation
  order M per tile, driven by the ball value of S_h — only near-corner tiles
  need large M. This is the d-weighted tail law; time a pilot before any
  batch.
- The majorant layer is symmetric under permutations of the coordinates
  (S depends only on the multiset {|phi_I|}); no such collapse is claimed
  for the Gram layer — evaluate per tile.

--- B7_CAUCHY.md ends ---
