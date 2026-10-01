# abcount — manual

abcount is the computational core of the ABACUS package: certified point
counts, Frobenius characteristic polynomials, and L-factors for abelian
fourfolds presented as a period lattice with polarization and, optionally, a
Z[sqrt(-5)]-action. It combines certified theta constants (Arb/acb_theta,
with a proved tail bound that dominates the library's claimed radii), exact
integer Weil-box isolation of Frobenius coefficients, and exact receipts for
the O_K-action. This file is the mathematical authority for the package: it
states and proves the inequalities the code relies on, records the exact
semantics of the isolation and receipt layers, and explains how to run the
validation battery. Layout and invocation are in the package README and
GUIDE (`tools/abacus/`).

Throughout, "desk" means an independent in-house computation (a direct
partial sum or exact integer arithmetic) used to check a library result,
and a "receipt" is a machine-readable record of a named check in the run
output.

## 1. The theta tail bound (genus-g box truncation, Gaussian comparison)

**Setting.** tau in H_g (Im tau = Y positive definite), lambda =
lambda_min(Y) > 0 the smallest eigenvalue of Y. Half-integer characteristic
[alpha; beta], alpha, beta in {0, 1/2}^g (FLINT encoding: 2g-bit integer
a_1..a_g b_1..b_g, alpha = a/2, beta = b/2). Theta constant

  theta[alpha;beta](0, tau) = sum_{n in Z^g} exp( pi*i*(n+alpha)^T tau (n+alpha)
                                                  + 2*pi*i*(n+alpha)^T beta ).

Desk partial sum S_N = the sum restricted to the box ||n||_inf <= N.

**Theorem (tail bound).** For every tau in H_g, every half-integer
characteristic, and every integer N >= 0, with lambda = lambda_min(Im tau):

  | theta[alpha;beta](0,tau) - S_N |  <=  TailBound(g, N, lambda)
      :=  g * S_tail(N, lambda) * M(lambda)^(g-1),

  S_tail(N, lambda) = 2 * exp(-pi*lambda*(N + 1/2)^2) / (1 - exp(-pi*lambda*(2N+1))),
  M(lambda)         = 1 + 2 * exp(-pi*lambda/4) / (1 - exp(-pi*lambda)).

**Proof.**
(i) *Term modulus (Rayleigh).* |exp(pi*i*(n+alpha)^T tau (n+alpha))| =
exp(-pi*(n+alpha)^T Y (n+alpha)) and (n+alpha)^T Y (n+alpha) >=
lambda*||n+alpha||_2^2, since lambda is the smallest eigenvalue of the
symmetric positive-definite Y. The beta term is a pure phase. Hence every
term of the tail is bounded by exp(-pi*lambda*||n+alpha||^2).

(ii) *Union bound.* If ||n||_inf > N then |n_k| >= N+1 for at least one
coordinate k, so the tail sum is bounded by the sum over k = 1..g of the sum
over {n : |n_k| >= N+1}. Each such sum factorizes over coordinates
(separability of exp(-pi*lambda*sum_j (n_j+alpha_j)^2)) into
S_tail(N, lambda, alpha_k) * prod_{j != k} S_full(lambda, alpha_j).

(iii) *1-D tail lemma.* For alpha in {0, 1/2},
S_tail(N, lambda, alpha) := sum_{|m| >= N+1} exp(-pi*lambda*(m+alpha)^2)
<= S_tail(N, lambda). Proof: split into the branches m = N+1+j and
m = -(N+1+j), j >= 0. Positive branch: (m+alpha)^2 >= (N+1+j)^2 >=
(N+1/2)^2 + j*(2N+1), because (N+1+j)^2 - (N+1/2)^2 - j*(2N+1) =
N + 3/4 + j + j^2 >= 0. Negative branch: (m+alpha)^2 = (N+1+j-alpha)^2 >=
(N+1/2+j)^2 >= (N+1/2)^2 + j*(2N+1). So each branch is dominated by the
geometric series exp(-pi*lambda*(N+1/2)^2) * sum_{j>=0} r^j with
r = exp(-pi*lambda*(2N+1)) < 1; two branches give the factor 2.

(iv) *1-D full-sum lemma.* For alpha in {0, 1/2},
S_full(lambda, alpha) := sum_{m in Z} exp(-pi*lambda*(m+alpha)^2) <=
M(lambda). Proof: at most one m has |m+alpha| < 1/2 (its term is <= 1); the
remaining m split into two monotone branches with values v_j = v_0 + j,
v_0 >= 1/2, and v_j^2 >= v_0^2 + 2*j*v_0 >= 1/4 + j, so each branch is
<= exp(-pi*lambda/4) * sum_{j>=0} exp(-pi*lambda*j) =
exp(-pi*lambda/4)/(1 - exp(-pi*lambda)).

(v) Combining (i)-(iv) with the triangle inequality over the tail proves the
theorem. QED.

The bound is implemented as `abcount_s1.tail_bound_c7`; the desk partial
sums are `abcount_s1.desk_theta_g4_all` (block-diagonal genus 4) and
`abcount_s2.desk_theta_g2_all` (genus 2). All bound evaluations are
performed in Arb upper-bound (ball) arithmetic — the bound quoted in a
receipt is a certified upper bound of the mathematical expression above.

**Block-diagonal specialization.** For tau = diag(tau_1, ..., tau_4) in H_4
(four genus-1 blocks), Y is diagonal, so lambda_min(Y) = min_k Im(tau_k)
exactly. After SL2(Z) reduction of each block, Im(tau_k) >= sqrt(3)/2.

**The dominance gate (why acb_theta radii are trusted).** One pilot genus-4
acb_theta call runs on an exact dyadic pilot tau, pinned in
`battery/SEEDS.json` (key `S1_C7_pilot`), at matched precision
prec = ceil(-log2(TailBound(4, N, lambda))) + 64 guard bits with the
pre-committed desk truncation N = 12, so the wrapper and the desk target the
same truncation depth. Gate, per theta constant:

  acb_theta claimed radius  <=  TailBound(4, 12, lambda)

— the proved desk bound must DOMINATE the wrapper's claimed radius at
matched truncation; any violation stops the run. The desk enclosure
[S_12 (ball) +/- TailBound] must additionally OVERLAP the wrapper enclosure
(consistency; disjoint enclosures also stop the run). Only after this
receipt passes are the wrapper's radii trusted for the rest of the battery.

**Genus 2, non-diagonal tau.** The inequality above is genus-general and
needs no new proof. The desk side is the direct 2-D box partial sum (no
factorization shortcut off block-diagonal), with lambda_min(Im tau)
lower-bounded by the closed-form 2x2 eigenvalue in ball arithmetic
(`abcount_s2.lambda_min_lower_2x2`). The calibrated gates:
matched-truncation dominance on the 10 EVEN characteristics, enclosure
overlap on ALL 16, and exact-zero containment on the 6 ODD characteristics.
The odd thetas are exactly 0 at z = 0, and acb_theta_all reaches them by
sqrt-of-square, so its honest claimed radius there is ~2^(-prec/2) — a
rounding artifact, not a truncation tail; that is why the odd gate is
containment of the exact zero, not dominance. Any violation of the
calibrated gates stops the run.

## 2. The Weil-box inequality set at degree 8 (abelian fourfolds over Q)

For an abelian fourfold B/Q with good reduction at p, the charpoly of Frob_p
is P(T) = prod_{i=1..8} (T - alpha_i) = T^8 + a_1 T^7 + ... + a_8 (integer
coefficients; L-factor convention as fixed in the input contract,
`M3T_REGISTRATION.md` sec 6), with |alpha_i| = sqrt(p) (Weil, proved for
abelian varieties — the Riemann-hypothesis part of the Weil conjectures).

**Bounds (used by the isolation layer):**
- |a_k| <= C(8,k) * p^(k/2) for k = 1..4: a_k is (+/-) the k-th elementary
  symmetric function of the alpha_i; it is a sum of C(8,k) monomials each of
  modulus p^(k/2). QED.
- Functional equation a_(8-k) = p^(4-k) * a_k for k = 0..3 (so a_5..a_8 are
  determined and a_8 = p^4): the Weil pairing gives alpha -> p/alpha as a
  permutation of the roots, i.e. T^8 * P(p/T) = p^4 * P(T); comparing
  coefficients yields the identity. QED.

**Isolation semantics** (`abcount_s1.weil_box_isolate_deg8`): FULL
ENUMERATION of integer tuples (a_1, ..., a_4) inside the box intersected
with the enclosure, never per-coefficient rounding. Exactly one tuple ->
emit; zero -> `FAIL-WEIL-BOX-EMPTY`; more than one -> `UNDECIDED-PRECISION`
with a needed-dps estimate. Integer-exact assembly inputs give zero-radius
enclosures (isolation trivially unique; the enumeration machinery still
runs).

## 3. The Weil-box inequality set at degree 4 (abelian surfaces over Q)

For an abelian surface A/Q with good reduction at p, the charpoly of Frob_p
is P(T) = prod_{i=1..4} (T - alpha_i) = T^4 + a_1 T^3 + a_2 T^2 + a_3 T +
a_4, integer coefficients, |alpha_i| = sqrt(p) (Weil).

**Bounds (used by the degree-4 isolation layer):**
- |a_k| <= C(4,k) * p^(k/2) for k = 1..4: same elementary-symmetric argument
  as sec 2 — a_k is a sum of C(4,k) monomials each of modulus p^(k/2). QED.
- Functional equation a_(4-k) = p^(2-k) * a_k for k = 0..1 (so a_3 = p a_1
  and a_4 = p^2): alpha -> p/alpha permutes the roots, T^4 P(p/T) = p^2 P(T);
  compare coefficients. QED.

**Isolation semantics** (`abcount_s2.weil_box_isolate_deg4`): identical
discipline to sec 2 on the free pair (a_1, a_2): FULL ENUMERATION of integer
pairs inside the box intersected with the enclosures, a_3 and a_4 then
FORCED by the functional equation and verified against their own
enclosures; never per-coefficient rounding. Degree-8 objects assembled as
squares/products of genus-2 blocks use the sec-2 box unchanged
(`abcount_s1.weil_box_isolate_deg8`).

## 4. What the routes compute (the battery members)

The validation battery (`battery/SEEDS.json`; each run banks its results to
`battery/BATTERY_RESULT.json`) exercises the tool on positive members B1-B5
and negative controls N1-N5. Three route families appear throughout:
Route T (certified theta receipts — emits NO Frobenius data), Route D
(exact integer assembly + Weil-box isolation), and PARI-based paths used as
an independent cross-check ("truth side"); every cross-method comparison is
between exact integers and must match exactly.

### B1 — product of four elliptic curves (block-diagonal tau in H_4)

- Route D: per-factor a_p from in-house naive point counts (no PARI) ->
  product assembly of the four degree-2 L-factors into the degree-8 L_p by
  exact integer polynomial arithmetic -> Weil-box isolation ->
  #B1(F_p) = P(1). Cross-check path: PARI ellap per factor. Compared
  integers: a_p(E_i), i = 1..4, and the assembled degree-8 L_p.
- Route T: certified genus-4 theta receipts on the block-diagonal tau —
  odd-characteristic vanishing (120 balls contain 0), product factorization
  against the per-factor genus-1 routes (Arb + desk), per-factor Jacobi
  identity, positivity, radius discipline; plus one PSLQ recognition
  positive control (theta -> j on one block, held out from everything else).

### B2, B3 — simple genus-2 Jacobians

- **Member instantiation (analytic Jacobian, NON-real-ordered branch
  points).** The Eichler wrapper's real-ordered period path does not apply
  to x^5 - x + 1 (the B2 factor curve) or x^5 + 1 (B3) — one real root plus
  two conjugate pairs each; its docstring exposes the primitives for complex
  chains. `period_g2_bridge.jl` assembles them (the Eichler files themselves
  stay untouched): certified chord integrals over the chain committed in
  SEEDS.json (trig desingularisation, per-factor half-plane-certified
  principal square roots, Arblib.integrate with the analytic checker), then
  homology recovery per the SEEDS-committed rule (LLL on the antisymmetric
  Riemann-relation kernel; Pfaffian +-1 gate; exact integer symplectic
  reduction verified by S W S^T = J; tau = A^(-1)B with certified symmetry
  and Im tau > 0). The construction is tied to the member curve by a tie
  receipt: Moebius branch-matching of thomae_sextic(tau) against the curve's
  branch set, chi5 certified nonzero, plus anchored covariant cross-products
  where available; cross-member controls must FAIL.
- **Counting routes.** B2: Path 1 = PARI hyperellcharpoly on C (independent
  cross-check); Path 2 = in-house desk counts #C(F_p), #C(F_p^2) (no PARI)
  -> power-sum quartic -> Route-D SQUARING assembly (exact integer, deg 8)
  -> sec-2 Weil box; compared integers = the genus-2 quartic AND the
  assembled deg-8 L_p; #B2(F_p) = P(1)^2. B3: Path 1 = PARI
  hyperellcharpoly; Path 2 = the `hecke_desk.py` Jacobi-sum closed form
  (exact Z[zeta_5] desk arithmetic, no PARI code path) under the committed
  sign pin alpha_a = -J(phi, chi^a), with the +J discrimination receipt (the
  odd coefficients of the two candidate signs must differ); deg-4 Weil box
  (sec 3).
- **Route T on both members:** genus-2 theta receipts (16 characteristics,
  6 odd vanishings, radius discipline, genus-2 dominance gates of sec 1)
  and, for B2, the genus-4 block diag(tau, tau) receipts (256 values, 120
  odd vanishings, product factorization theta4 = theta2 x theta2). Route T
  emits NO Frobenius data.

### The O_K-action layer (K = Q(sqrt(-5)); module `abcount_s3`)

All predicates are named, machine-readable, and refuse rather than guess.

Input: integer 8x8 matrix rho on the product-lattice basis (per factor
{1, rho}), alternating integer polarization E, period matrix Pi (balls).
Exact integer predicates, in firing order:
- `OK-action-rho-square-eq-minus5`: rho^2 = -5*Id.
- `OK-action-lattice-preservation`: integrality on the declared basis.
- `OK-action-rosati-scaling`: rho^T E rho = 5 E.
- `OK-action-rosati-antiadjoint`: rho^T E + E rho = 0 (Rosati(rho) = -rho).

Ball-certified analytic predicates:
- `rho-C-linearity(A*Pi=Pi*M)`: the analytic representation A solved from
  the identity columns, verified on the remaining columns.
- `rho-analytic-square-eq-minus5`: A^2 = -5.
- `rho-holomorphic-eigenvalue-scalar` plus the SIGN PIN: A = s*i*sqrt(5)*Id,
  and the computed s must equal the sign declared with the input; a mismatch
  is `FAIL-SIGN-PIN`. The pin check VERIFIES or REFUSES only — it never
  converts; the single sign-conversion site is
  `abcount_s0.convert_sign_convention` (sec 6).

Polarization typing: `polarization-type-frobenius-normal-form` — Smith
normal form of E, symplectic pair structure, elementary divisors checked
against the declared type; wrong divisors fail HERE, before any tau is
emitted.

### B4, B4s — CM elliptic curve over Q(sqrt5), deg-16 assembly

E20: j = 632000 + 282880*sqrt5 (root of the Hilbert class polynomial
H_-20), pinned standard model y^2 = x^3 + 3j(1728-j)x + 2j(1728-j)^2 over
F = Q(sqrt5); twist sign u = -1 committed with the model in SEEDS.json
(a_P = -6 at p = 29, #E20(F_P) = 36 = N(pi+1)).
- Path 1 (independent cross-check): PARI ellap on the reduced model per
  prime P of F; at inert-in-K primes where the standard model degenerates
  (j = 0 or 1728 mod P — both primes above 11), ellfromj under the
  twist-invariance license (supersingular => a_P = 0 for every twist).
  At p = 29 additionally the E[3] torsion congruence: the explicit
  Frobenius matrix F on E[3] (division-polynomial basis over F_29^2) must
  satisfy charpoly(F) = T^2 - a_P T + p mod 3, and R := 2^{-1}(u*F - 3*Id)
  must satisfy R^2 = -5, tr R = 0, det R = 5 mod 3, R nonscalar; the
  claimed per-prime b_P signs must satisfy b_pin = -b_conj, |b| = 2.
- Path 2 (`hecke_desk` O_K extension, PARI-free): inert => closed form
  a_P = 0 (supersingular + Weil bound). Split => Cornacchia
  a^2 + 5b^2 = p, Frob_P = u*(a + b*rho) at the committed K-prime
  ((29, sqrt(-5) - 13) contains pi = 3 + 2sqrt(-5)), u*(a - b*rho) at the
  conjugate; output pair (a_P, b_P) = (2ua, +/-ub) in the committed basis.
- Emission rule (the rho/-rho DETECTOR, member B4s): the emitted b_P is
  s * b_P^CM — negating the input rho flips s and therefore flips every
  emitted b_P, while a_P and all outputs of the rho-free routes are
  invariant. Check s before consuming b_P.
- Assembly: L_p(B4/F, T) = prod_{P|p} L_P(E20, T^{f(P)})^4 with explicit
  degree bookkeeping; deg-16 functional equation c_{16-k} = p^{8-k} c_k;
  per-prime Weil boxes by full enumeration; at 29 the pair box
  (a/2)^2 + 5 b^2 = 29 enumerates to exactly {(+-6, +-2)}.

### B5 — Weil restriction of a genus-2 Jacobian over K

C_K: y^2 = x^5 + sqrt(-5)x + 1 over K, good reduction at {11, 29, 31}
(discriminant-norm receipts in SEEDS.json). Precondition: PARI
hyperellcharpoly over F_p AND F_p^2 must be available; the harness
re-probes this precondition in-run and refuses rather than skipping.
- Path 1: hyperellcharpoly over the residue field (F_p at split 29 for
  both conjugate reductions, sqrt(-5) -> 13 and the conjugate 16; F_p^2 at
  inert 11, 31 via ffgen towers).
- Path 2 (PARI-free): pure-python tower counts (squares-set character,
  power sums as in B2) over F_q and F_{q^2}, q = p or p^2; the F_11^4
  enumeration is timed before the F_31^4 sweep so the larger sweep is
  priced by measurement before it starts.
- Assembly: L_p(Res_{K/Q} B', T) = prod_{P|p} L_P(B', T^{f(P)}) — inert:
  quartic in T^2; split: product of the two conjugate quartics; deg-8
  functional equation + Weil box + count receipts L_p(1).

### N1-N5 — negative controls (each must FAIL by its named predicate)

- N1: tau_(1,2) += 1e-6 on B1 periods: `riemann-relation-tau-symmetry`
  fires (diagonal entries are excluded: a diagonal perturbation would be
  consistency-preserving).
- N2: a polarization with elementary divisors [1,1,2,2] against a declared
  principal type: `polarization-type-frobenius-normal-form`; no tau emitted.
- N3: rho[1][2] += 1: `OK-action-rho-square-eq-minus5`, with the residual
  entries quoted.
- N4: enclosure radii inflated to 1.5 (a 0.6 inflation would admit only one
  integer and defeat the control): `UNDECIDED-PRECISION` with a needed-dps
  estimate, NO integer emitted.
- N5: conjugated complex structure submitted under the committed sign:
  `FAIL-SIGN-PIN` (a silent same-answer would fail the whole battery; not
  observed).

## 5. Refusal semantics

The tool never guesses. Every check is a named predicate and every
non-success is a named verdict in the run output: `FAIL-<named-check>`,
`FAIL-WEIL-BOX-EMPTY`, `FAIL-SIGN-PIN`, `UNDECIDED-PRECISION` (always with
a needed-dps estimate), `UNDECIDED-PRECONDITION`. Every receipt carries the
actual input's own invariants, and every receipt carries
`pinned_lattice_claim: NONE` — no identification of the input against any
externally fixed lattice is made or implied.

## 6. Conventions

- `sign_convention` is mandatory on every artifact (including explicit
  "N/A" values such as "N/A-B2"/"N/A-B3" — carried, never dropped).
- The SINGLE sign-conversion site is `abcount_s0.convert_sign_convention`.
  Never add a second conversion site; every other sign check verifies or
  refuses only.

## 7. Running the checks

Dependencies: `python-flint` (Arb/acb_theta); PARI/GP (`gp`) on the
cross-check paths; Julia with the Eichler.jl project for the theta/period
bridges (set `ABACUS_EICHLER_PROJECT` to its root — the package ships in
this repository under `upgrades/Eichler.jl/`; the bridges refuse loudly
when it is unset); the in-repo Lockpick toolkit (`tools/lockpick/`, added
to the path automatically by the harnesses) for the j-recognition control.

From the repository root:

- `python3 tools/abacus/build/abcount/run_s0.py` — the genus-1 core checks;
  expected verdict `S0-PASS` (7/7 checks).
- `python3 tools/abacus/build/abcount/run_s3.py` — the harness of record:
  the FULL battery (B1-B4, B4s, B5, N1-N5). Expected verdict
  `BATTERY-PASS(11/11, B5 RUNS+PASS)`, 22/22 checks. Each run banks its
  result receipt to `battery/BATTERY_RESULT.json` (measured 45.4 s wall /
  0.104 GB RSS against the budget of 3600 s / 4 GB declared in
  `battery/SEEDS.json` (`S3_walls`) before the run); expect a similar
  footprint.
- `run_s1.py` / `run_s2.py` pin earlier revisions of SEEDS.json and REFUSE
  against the shipped one; their logic is exercised inside `run_s3.py`
  (native B1 re-exercise; run_s2 replay for B2/B3). This refusal is
  deliberate, not a bug.

Binding law (also in the README): abcount outputs are EVIDENCE only behind
the shipped battery. Any change to the tool re-runs the FULL battery
(`run_s3.py`) before its outputs are used.
