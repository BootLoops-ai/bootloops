# M3T VALIDATION BATTERY — abcount (spec companion)

Synthetic-truth controls bind: the tool's outputs are evidence only
after planted-truth positives + negative controls pass in its own validation
battery. The battery lives with the tool; results bank beside the code
(`battery/BATTERY_RESULT.json`).

## Scope

Every battery member is an INPUT INSTANCE in the tool's declared input contract
(M3T_REGISTRATION.md sec 6, format A or B) with a PLANTED TRUTH computed by an
independent in-house method, plus a required verdict. Seeds: every numerically
generated object (period matrices, perturbations, prime choices) derives from a
committed seed file `battery/SEEDS.json` (integers, written before any run;
sha-pinned). No member's truth may be computed by abcount itself.

Pilot primes per member: base set
p in {11, 31}, with the following PRE-COMMITTED corrections. (i) Cremona 11a1
(B1 factor) has BAD reduction at p = 11 (11 | disc = -11^5), so B1's
pilot-prime set is {13, 31}, pre-committed HERE and echoed in
`battery/SEEDS.json` — never discovered at run time. Desk truths at p = 13
(naive counts, independent of `ellap`, banked by pure-integer desk
arithmetic): a_13(11a1) = 4, a_13(14a1) = -4,
a_13(15a1) = -2, a_13(17a1) = -2; #B1(F_13) = 10*18*16*16 = 46080; 13 divides
none of the four discriminants. (ii) Both base pilot primes are INERT in
K = Q(sqrt(-5)) ((-20|11) = (-20|31) = -1), so the K-action members carry the
additional SPLIT pilot prime p = 29 ((5|29) = +1 AND (-5|29) = +1: 29 splits
in the base field Q(sqrt(5)) — residue fields F_29, Frob_P of degree 1 — and
in K, making the O_K-action per-prime-visible; member B4s below). Any residual
bad-reduction collision found at instantiation still follows the fallback:
the seed file replaces the prime with the next good one satisfying the SAME
splitting conditions, receipt filed.

## (a) Planted-truth positives (each seed-committed)

| ID | Member | Planted truth (in-house, independent) | Modules exercised |
|---|---|---|---|
| B1 | E1 x E2 x E3 x E4, four pairwise non-isogenous non-CM elliptic curves over Q (Cremona 11a1, 14a1, 15a1, 17a1; period lattices to declared dps; product principal polarization); pilot primes {13, 31} (p = 11 excluded: 11a1 bad; desk truths at 13 banked in Scope) | a_p(Ei) via PARI `ellap`; L_p(B1) = prod L_p(Ei); base_field = Q | tau constructor (block-diagonal), genus-4 theta consistency, Weil-box integer isolation, route-D trivial split |
| B2 | J(C)^2, C: y^2 = x^5 - x + 1 (genus 2, generic Jacobian; analytic Jacobian to declared dps; product polarization) | charpoly of Frob_p on J(C) via PARI `hyperellcharpoly`; square it | non-diagonal genus-2 tau block, theta vs curve counts, simple-factor handling |
| B3 | J(C5), C5: y^2 = x^5 + 1 (genus 2, CM by Q(zeta_5)) — the CM abelian variety with known L-factor | TWO independent truths: (i) `hyperellcharpoly`; (ii) Jacobi-sum / Hecke-character closed form (desk arithmetic, exact) with the SIGN PINNED: Frobenius eigenvalues alpha_a = -J(phi, chi^a) (verified at p = 11, 31: the -J polynomial matches the count-derived charpoly exactly; the +J polynomial differs in the odd coefficients — an instantiation against +J would falsely FAIL a correct tool) | CM route, L-factor assembly from characters, cross-method member par excellence |
| B4 | E20^4, E20 with CM by O_K = Z[sqrt(-5)] (disc -20, class number 2; j = 632000 + 282880*sqrt(5), root of H_{-20} = x^2 - 1264000x - 681472000; curve over base_field F = Q(sqrt(5)), the ring class field's real subfield), O_K acting DIAGONALLY; explicit rho(sqrt(-5)) as 8x8 integer matrix on the product lattice | PLANTED TRUTH stated in the base_field convention: base_field = F = Q(sqrt(5)); for each prime P of F above p (both degree 1 when (5\|p) = +1), a_P(E20) via `ellap` over F_P cross-checked by the in-house `hecke_desk` Hecke-character evaluation (registration sec 1); planted truth = the induced L_p(B4/F, T) = prod over P\|p of L_P(E20, T^f(P))^4 with the degree bookkeeping explicit | THE O_K-action verification module (exact rho receipts: rho^2 = -5, rho.Lambda = Lambda, Rosati compatibility), sig bookkeeping, sign-convention pin |
| B4s (the embedding-sensitive member) | Same object as B4 (E20^4, same rho, same polarization) at the SPLIT pilot prime p = 29: (5\|29) = +1 (29 splits in F = Q(sqrt(5)), residue fields F_29, sqrt(5) -> {11, 18}, j mod P in {12, 23}) AND (-5\|29) = +1 (29 splits in K); 29 = 3^2 + 5*2^2, so the O_K-primes above 29 are PRINCIPAL: pi = 3 + 2*sqrt(-5), pibar = 3 - 2*sqrt(-5) | REQUIRED OUTPUT (embedding-sensitive, per prime P of F above 29): the pair (a_P, b_P) with Frob_P = (a_P/2) + b_P * rho as an element of O_K in the PINNED basis {1, rho}; conjugation-aware desk truth: Frob_P = u*(3 + 2*rho) at the O_K-prime selected by the declared K-embedding and u*(3 - 2*rho) at its conjugate (twist sign u in {+1, -1} committed with the model in SEEDS.json; #E20(F_P) = 24 or 36 resp., \|a_P\| = 6 <= 2*sqrt(29)); rho -> -rho flips the SIGN of b_P (4u*rho != 0) while every B1-B4 legacy output is invariant — this member is the rho-embedding detector; b_P verified on E[m]-torsion (m odd, coprime to 10p) against 2^(-1)*(u*Frob_P - 3) mod m; desk derivation banked at seed time | K-EMBEDDING orientation: a tool that silently uses -rho for rho passes B1-B4 and N1-N5 but FAILS B4s; sign_convention K-embedding pin exercised |
| B5 (stretch, stage-3 optional) | Res_{K/Q}(J(C_K)) for a seed-committed genus-2 curve C_K over K = Q(sqrt(-5)): an abelian FOURFOLD with genuine K-action in the target class | count C_K over residue fields F_p / F_p^2 via `hyperellcharpoly` over the extension (in-house support confirmed at build time; if absent, B5 is dropped with a receipt and Q4 stands sharper) | the only member in the actual target class (simple-over-Q possible, K-action non-diagonal); answers Q4 if it lands |

Positives pass criterion (per member, per prime): abcount's output L-factor /
point-count integers EQUAL the planted truth EXACTLY (integer equality, no
tolerance), with the certified enclosure having isolated a UNIQUE integer
tuple inside the Weil box. Any non-integer, any multi-candidate box, any
mismatch = member FAIL.

## (b) Negative controls (perturbed inputs that MUST FAIL loudly)

| ID | Perturbation (seed-committed) | Required verdict |
|---|---|---|
| N1 | B1 period matrix, ONE entry perturbed by an amount strictly exceeding the declared per-entry radius, the entry PINNED PREDICATE-VISIBLE in `battery/SEEDS.json`: either an OFF-DIAGONAL tau entry of the normalized (I\|tau) form, or any single entry of a general-position (A\|A*tau) matrix — both verified detectable (residual ~1e-6 vs radius 1e-20); a DIAGONAL tau entry of the normalized form is EXCLUDED (it yields the consistent period matrix of a different abelian variety — all named predicates pass, verified) | FAIL at Riemann-relation / polarization-integrality consistency check, named check in the receipt; AND (second channel, for any consistency-preserving perturbation that slips the pin) cross-method mismatch vs the planted truth = FAIL — silent Weil-legal wrong counts are a battery FAIL of the whole registration |
| N2 | B2 with polarization E replaced by an alternating integer matrix of wrong elementary divisors | FAIL at Frobenius-normal-form / polarization-type check; no tau emitted |
| N3 | B4 with rho(sqrt(-5)) tampered in one entry (+1) | FAIL at exact O_K-action receipt (rho^2 = -5 or lattice-preservation or Rosati check); named failing predicate |
| N4 | Precision-starvation control: B1 input radii inflated until the Weil box contains >1 integer candidate | verdict UNDECIDED-PRECISION with a needed-dps estimate; ANY emitted integer = battery FAIL (the tool must refuse, never guess) |
| N5 | Sign-convention flip: B4 with the complex structure conjugated (the P vs P(-1) convention pin) | either a COHERENT declared-convention output flagged as the conjugate convention, or FAIL at the sign-pin check — WHICH one is PRE-DECLARED as a field in the sha-pinned `battery/SEEDS.json` WRITTEN BEFORE ANY RUN (the pre-declaration rides the seeds sha, no post-hoc choice); silent same-answer = battery FAIL |

Negative pass criterion: the required verdict class EXACTLY, with the failing
predicate NAMED in machine-readable output. A negative control "passing" as a
positive is a battery FAIL of the whole registration.

## (c) Cross-method agreement (per member, WHICH integer each route emits
## and its FULL code path; route T scoped to receipts)

Route T CANNOT emit an L-factor integer except through recognition -> model ->
the same PARI counting code as the truths. Its battery role is
therefore SCOPED to what it actually delivers: certified theta-consistency
receipts (Riemann relations, positivity, radius discipline) on every member,
plus ONE recognition-of-model positive control for `pslq_gate` (B1 block: theta
-> j-invariant recognition, held-out gate, never an L-factor emitter). The
integer-vs-integer cross-method pairs are:

| Member | Integer compared | Path 1 (full code path) | Path 2 (full code path) |
|---|---|---|---|
| B1 | a_p(Ei), i = 1..4, and assembled deg-8 L_p | PARI `ellap` per factor (planted truth) | abcount route-D product assembly from per-factor inputs + Weil-box isolation (abcount-internal, no PARI call in the assembly layer) |
| B2 | genus-2 charpoly of Frob_p, then its square | PARI `hyperellcharpoly` on C (planted truth) | abcount route-D squaring/product bookkeeping + Weil-box isolation (abcount-internal) |
| B3 | genus-2 charpoly (CM) | PARI `hyperellcharpoly` | `hecke_desk` Jacobi-sum closed form (in-house desk engine, registration sec 1 Builds item 3a — NO PARI code path; sign pin alpha = -J(phi, chi^a)) |
| B4 | per-prime a_P over F = Q(sqrt(5)) and induced L_p(B4/F, T) | PARI `ellap` over the residue fields | `hecke_desk` Hecke-character evaluation for K (in-house, PARI-independent — the named route-D implementation) |
| B4s | the pair (a_P, b_P) per prime P above 29 (b_P = rho-coefficient of Frob_P) | E[m]-torsion congruence check of rho vs 2^(-1)(u*Frob - 3) computed with `ellap`-side reductions | `hecke_desk` conjugation-aware desk truth u*(3 +/- 2*rho) per O_K-prime (banked at seed time) |
| B5 | deg-8 charpoly of Res_{K/Q} | PARI `hyperellcharpoly` over residue fields F_p / F_p^2 | abcount route-D Weil-restriction identity L_p(Res B') = prod over P\|p of L_P(B') assembled from `hecke_desk`-independent per-prime inputs |

Agreement = integer equality of the compared integers. One path may be the
planted-truth generator for a member ONLY if the other path is abcount-internal
or `hecke_desk` (the point stands: two computation paths sharing no code beyond
exact integer arithmetic; theta receipts ride every member but are never the
compared integer).

## (d) Acceptance criterion (registration-binding)

- abcount outputs are EVIDENCE only after: ALL positives pass at ALL declared
  primes + ALL negative controls return their required verdicts + ALL
  cross-method pairs agree. B4s is CORE, not optional (without it
  the K-action arithmetic is unexercised anywhere). Partial battery = tool
  stays NON-EVIDENTIAL (outputs consumable as observation-grade diagnostics
  only, labeled so).
- Battery results bank WITH the tool: `battery/BATTERY_RESULT.json`
  (per-member verdicts, enclosure radii, seeds sha, wall/mem actuals).
- EVERY subsequent tool change re-runs the FULL battery before the change's
  outputs are evidence (every-tool-change law rider; eras-style
  VERIFIER-RERUN discipline — the verifier is re-run, result.ok checked,
  never quoted from memory).
- B5 dropped => acceptance still available, but the registration carries a
  TARGET-CLASS-UNEXERCISED flag that rides every abcount output until a
  target-class member passes (Q4 decides whether that flag blocks
  downstream gate-use).
