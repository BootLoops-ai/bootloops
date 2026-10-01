# abcount — spec of record (theta / point-count for abelian fourfolds)

Synthetic-truth controls bind: the point-count/theta tool's outputs
are evidence only after planted-truth positives + negative controls pass in
its own validation battery; the battery is designed into this spec.
Companion (same directory): M3T_VALIDATION_BATTERY.md.

## 0. What the tool is for

Isogeny-explicit abelian data for B: an explicit period lattice with
Z[sqrt(-5)]-action + polarization for an abelian fourfold B (sig (2,2)),
plus an executable in-house theta/point-count algorithm — nothing else in
the toolkit counts abelian fourfolds.
Branch-agnostic requirement: the spec accepts both input shapes
(sec 6) — the analytic shape (period balls) and the hermitian-
lattice shape (exact O_K-module data).

## 1. TOOL SPEC

**Name: `abcount`** (shipped as the `tools/abacus` package; the `abcount_*`
file stems are the sha-pinned module/harness surface).

**Inputs** (full contract in sec 6; summary):
- a period lattice for an abelian fourfold with O_K = Z[sqrt(-5)]-action
  (the action matrix rho is INPUT, verified not discovered — evidential mode)
  + polarization data (integer alternating matrix E with declared elementary
  divisors) + mandatory sign-convention pin covering BOTH signs: the
  complex-structure convention (P vs P(-1); Q6)
  AND the K-embedding orientation (which square root of -5 the input matrix
  rho represents); the full pin lives in sec 6;
- a prime p (good-reduction declaration by the caller; the tool checks what
  it can and flags the rest);
- precision params: input dps + per-entry radii; working dps; a wall
  (wallclock/mem) per call — priced-pilot style, no unpriced calls.

**Outputs** (all with error discipline — ball arithmetic throughout,
result.ok-style verdict fields, eras-style verifier re-run law):
- certified theta constants (genus-4 theta-null values with radii) for the
  constructed tau in H_4;
- Frobenius data at p: the degree-8 integer charpoly of Frob_p over the
  declared `base_field` (equivalently #B and the L-factor in the sec-6
  convention) — emitted ONLY when the certified enclosure isolates a UNIQUE
  integer tuple inside the Weil box. The box is pinned
  exactly: (i) ENCLOSURE SOURCE — the enclosure is fed by route-D assembly
  evaluated in ball arithmetic (Hecke-character archimedean evaluation via
  `hecke_desk` and any non-exact intermediate of the product/induction
  bookkeeping); route C and integer-exact route-D paths emit EXACT integers =
  zero-radius enclosures (isolation trivially unique); route T feeds NO
  enclosure (it emits no Frobenius data — see route re-scope). (ii) INEQUALITY
  SET — |a_k| <= C(8,k) p^(k/2) for k = 1..4, functional equation
  a_(8-k) = p^(4-k) a_k for k = 0..3 (so a_5..a_8 are determined, a_8 = p^4),
  all root moduli = sqrt(p); the bounds are stated and proved in
  `manuals/abcount.md` at first build. (iii) ISOLATION SEMANTICS — FULL
  ENUMERATION of integer tuples (a_1..a_4) in the box intersected with the
  enclosure, NEVER per-coefficient rounding; exactly one tuple -> emit; zero
  tuples -> FAIL-WEIL-BOX-EMPTY; more than one -> UNDECIDED-PRECISION with a
  needed-dps estimate;
- consistency receipts: Riemann relations, polarization positivity (Im tau >
  0), O_K-action predicate (rho^2 = -5, lattice preservation, Rosati
  compatibility — Q2), sign-pin echo (BOTH pins, sec 6).
  The Q2(ii) certificate: for FORMAT A inputs the
  O_K-action receipt is EXACT — rho, J, E all derive from algebraic data and
  predicates (i)-(iii) are verified in exact integer/number-field arithmetic
  (the certificate object IS that exact verification transcript). For FORMAT B
  inputs exact J-commutation is impossible in principle (J comes from period
  balls); the receipt is the DECLARED CONDITIONAL class
  `O_K-ACTION-WITHIN-RADIUS` (commutator norm <= certified radius bound
  quoted), which RIDES every downstream output; evidential claims that need
  the exact action require format-A-grade rho provenance.
- Verdict classes (closed list, TOTAL): OK /
  FAIL-<named-check> / UNDECIDED-PRECISION / UNDECIDED-PRECONDITION /
  STOP-WALL / FAIL-INTERNAL. Any crash, exception, or fault inside a called
  library (e.g. a segfault in acb_theta) maps to FAIL-INTERNAL with the raw
  fault captured in the receipt — every possible outcome has a named class.
  Never a silent number.

**Algorithm routes (mathematical options stated honestly, with
preconditions):**
- **Route T — certified theta summation.** From (Lambda, E, rho): Frobenius
  normal form of E -> symplectic basis -> big-period matrix -> tau in H_4;
  theta-nulls via the EXISTING `acb_theta` wrapper (the Eichler package's
  `src/siegel.jl`, shipped in this repository under
  `upgrades/Eichler.jl`) with certified tails (tail-bound
  authority: Q1). Delivers theta constants
  unconditionally (given polarization positivity); delivers ARITHMETIC only
  through model recognition (theta-null ratios -> algebraic numbers via
  pslq_gate with mandatory positive controls; recognitions are HELD-OUT
  GATES, never evidence by themselves — the lockpick positive-control
  discipline).
  Route T emits NO Frobenius data and feeds no Weil-box enclosure; its
  battery role is theta-consistency receipts + one recognition positive
  control (battery sec c). Precondition: positive
  polarization; if NO positive polarization is compatible with the banked
  data, the verdict is UNDECIDED-PRECONDITION (the closed list covers this
  branch explicitly).
- **Route D — CM / isogeny decomposition.** If End^0(B) strictly contains K
  and splits B up to isogeny into factors with in-house-countable truth
  (elliptic factors via `ellap`; genus-2 Jacobian factors via
  `hyperellcharpoly`; CM factors via Hecke characters evaluated by the
  in-house `hecke_desk` engine — the NAMED independent implementation for
  B4/B4s/B5 route D, Builds item 3a), assemble L_p from
  the factors + the kernel receipt of the isogeny. Precondition: the extra
  endomorphisms are SUPPLIED with exact integer receipts (same predicate
  discipline as rho); abcount verifies and assembles, it does not conjure
  splittings in evidential mode. (A numerical endomorphism-hunt mode exists
  as observation-grade diagnostics only, clearly labeled.)
- **Route C — arithmetic-model counting.** When an algebraic model (curve
  whose Jacobian is B up to declared isogeny, or a projective model) is part
  of the input provenance, count directly: `ellap` / `hyperellcharpoly`
  (PARI) over F_p and extensions.
  Precondition: the model exists on the record. For analytic-period inputs
  route C is honest about being unavailable until the upstream computation
  yields a model or an isogeny (Q5).

**Calls (existing software consumed — none of this is rebuilt):**
- the Eichler package's `src/siegel.jl` — acb_theta genus-g theta wrapper
  (certified theta VALUES; the package ships in this repository under
  `upgrades/Eichler.jl` and the harnesses read its project root from
  `ABACUS_EICHLER_PROJECT`);
- PARI/GP `hyperellcharpoly` + `ellap` — ALL curve-side planted truths and
  route-C/D factor counts; abcount never re-implements curve counting;
- `tools/lockpick/pslq_gate.py` — algebraic recognition with mandatory
  positive controls (route T recognition layer);
- eras verifier pattern (`tools/eras/`) — battery verifier discipline
  (VERIFIER-RERUN rule, result.ok mandatory); certified-arithmetic substrate
  rule: DO NOT rebuild certified arithmetic.
- Named non-call: `posq` (certified Gauss-Jacobi evidence quadrature —
  different domain, no role here).

**Builds NEW (the actual tool body):**
1. hermitian/period front door: O_K-hermitian (Lambda, h) [format A] or
   period-matrix balls [format B] -> symplectic normal form -> tau in H_4,
   with the consistency-receipt suite (Riemann, positivity, rho predicate,
   sign pin);
2. Weil-box Frobenius-integer isolation layer (enclosure -> unique integer
   tuple or UNDECIDED-PRECISION);
3. route-D assembly (kernel receipts, L-factor product/induction
   bookkeeping, Weil-restriction identity);
3a. `hecke_desk`: in-house Hecke-character / Jacobi-sum
   evaluation engine (exact desk arithmetic over O_K and Z[zeta_5] residue
   fields, PARI-independent — shares no code with `ellap` /
   `hyperellcharpoly` beyond integer arithmetic) — the named route-D
   implementation for B3/B4/B4s/B5 cross-method truths; sign convention
   alpha = -J(phi, chi^a) pinned;
4. battery harness (`battery/` — seeds, members, verifier, BATTERY_RESULT
   banking).

## 2. VALIDATION BATTERY (designed in)

Full design: M3T_VALIDATION_BATTERY.md (same directory).
Summary, binding here: (a) planted-truth positives B1 (E1xE2xE3xE4, `ellap`
truths; pilot primes {13, 31} — p = 11 excluded, 11a1 bad there; a_13 desk
truths banked), B2 (J(C)^2 genus-2, `hyperellcharpoly`
truth), B3 (J(y^2 = x^5+1), CM, dual truths hyperellcharpoly + Jacobi-sum
closed form, sign pinned alpha = -J(phi, chi^a)), B4 (E20^4
with diagonal Z[sqrt(-5)]-action over base_field Q(sqrt(5)), `ellap` +
`hecke_desk` truths in the sec-6 base-field convention — exercises the rho
predicate), B4s (CORE; same object at the SPLIT prime p = 29, required
output the embedding-sensitive per-prime pair (a_P, b_P) — the rho-vs-(-rho)
detector), B5 stretch (Res_{K/Q} of a genus-2 Jacobian over
K — the target class; precondition-flagged); each seed-committed
(`battery/SEEDS.json`, sha-pinned, written before any run). (b) negative
controls N1-N5 (perturbed period entry — predicate-visible entry pinned in
SEEDS.json; wrong polarization type; tampered rho; precision
starvation -> must refuse; sign flip — either/or outcome pre-declared in
sha-pinned SEEDS.json) — each must FAIL with the named
predicate; a negative passing as positive fails the whole battery. (c)
cross-method: two independent INTEGER paths on every member with full code
paths tabled in the battery file; route T scoped to theta receipts + one
pslq_gate recognition positive control. (d) acceptance: abcount
outputs are EVIDENCE only after the FULL battery passes (B4s core); results
bank with the tool; every
subsequent tool change re-runs the full battery. B1-B4 sufficiency
vs B5 necessity: Q4.

## 3. EVERY-TOOL-CHANGE LAW (binding)

Any abcount change re-runs the full battery in the same cycle, and
`manuals/abcount.md` (the current-state manual) is updated with it. The
tool's one-line contract: certified theta constants / Frobenius charpoly /
L-factor for an abelian FOURFOLD given a period lattice with
Z[sqrt(-5)]-action + polarization — hermitian/period front door -> tau in
H_4, acb_theta certified values, Weil-box integer isolation, exact
rho/Rosati receipts, routes T/D/C with declared preconditions; EVIDENCE only
behind its banked battery (B1-B4 + B4s[+B5] + N1-N5).

## 4. STAGED BUILD PLAN (smallest correct core first, priced pilots)

| Stage | Content | Wall / bound | STOP semantics |
|---|---|---|---|
| S0-PRE | precondition pins: the Q3 (isogeny-class vs pinned lattice) and Q6 (the pin VALUE, both signs) resolutions stand as binding input-contract pins (S0PRE_RECEIPT.md) | none (a pin, not a run) | S0 does not start without them — else the input contract is two contracts |
| S0 | dimension-1 sanity: ONE elliptic curve — period lattice -> tau in H_1 -> theta nulls -> j recognition -> `ellap` comparison at ONE prime (p=11-class) | desk-class: minutes wallclock, <=2 GB, dps <= 60 | any theta-vs-ellap mismatch or radius blow-up = STOP, receipt, no stage 1 |
| S1 | products: battery member B1 (block-diagonal genus-4 tau); measure ONE acb_theta genus-4 call cost FIRST (priced pilot: one member, one prime) before the member sweep; independent tail-bound cross-check on that one call (Q1). The desk inequality (Gaussian-comparison bound on the ellipsoid tail, driven by the smallest eigenvalue of Im tau) is STATED AND PROVED in `manuals/abcount.md`, and the S1 receipt QUOTES the inequality; pass criterion PINNED: acb_theta's claimed radius <= the desk tail bound at matched truncation — ANY violation = STOP (never a shrug), no battery-wide trust in the wrapper's radii before this receipt | <=1 h, <=8 GB per member-prime | wall breach = STOP-WALL receipt (no silent retry); tail-bound dominance violation = STOP; cost actuals bank for S2 pricing |
| S2 | simple genus-2 blocks: B2, B3 (non-diagonal tau; cross-method vs hyperellcharpoly; Jacobi-sum desk truth for B3) | <=2 h, <=16 GB per member-prime, priced from S1 actuals | any cross-method disagreement = STOP; the disagreement is itself a banked receipt |
| S3 | the general (2,2) machinery: B4 (rho predicate wired), negative controls N1-N5, B5 if its precondition holds; FULL battery run + BATTERY_RESULT banking + the sec-3 compliance files | priced from S2 actuals before launch | full-battery pass = tool becomes evidential; any member fail = tool stays non-evidential, fix + FULL re-run |

No stage starts before the previous stage's receipt is banked. Every stage's
wall is declared BEFORE launch. UNDECIDED-* is a
verdict, never a retry loop.

## 5. DEPENDENCY NOTE (upstream input contract)

The tool consumes an upstream period computation's output (an
isogeny-explicit abelian fourfold as analytic period balls), OR
(branch-agnostic clause) a hermitian lattice. The contract below is stated
in full so producer and consumer can be checked for compatibility.

## 6. INPUT CONTRACT (what abcount accepts)

`INPUT.json`, one of two formats, both carrying the common block:
- COMMON: `sign_convention` — pins BOTH signs, (s1) the
  complex-structure convention P vs P(-1) (mandatory, Q6) AND (s2) the
  K-embedding orientation: which embedding iota: K -> C has iota(sqrt(-5))
  the value that `okaction` rho represents (equivalently: rho, not -rho, is
  declared to act as +sqrt(-5) under iota) — the two signs are independent
  for contracted inputs (no proof that one determines the other is on
  record, so BOTH are pinned; battery member B4s detects an s2 violation);
  `base_field` — the number field F over which the
  abelian variety and its L-factor are taken; "prime p" means the rational
  prime, and the emitted Frobenius data is the induced
  L_p(B/F, T) = prod over P|p of L_P(B, T^f(P)) with per-prime factors and
  residue degrees f(P) itemized in the receipt (for F = Q this collapses to
  the single degree-8 charpoly of Frob_p; B4/B4s declare F = Q(sqrt(5)));
  `provenance_sha` (upstream receipt); `precision` {dps, per-entry radius};
  `prime p`; `polarization` E as an 8x8 alternating integer matrix on the
  declared basis + its elementary divisors; `okaction` rho as an 8x8 integer
  matrix (the Z[sqrt(-5)] generator on the lattice basis, oriented per s2);
  optional `isogeny_receipt` (integer matrix + kernel data, Q3) and optional
  `model_provenance` (route-C enabler).
- FORMAT A (exact hermitian-module shape): exact O_K-hermitian module data — pseudo-basis
  (ideals a_1..a_4 + vectors) + hermitian Gram h in M_4(K) entries +
  signature declaration (2,2) + a declared complex embedding; abcount
  derives the period matrix from the CM-type/embedding data it is given
  (never guesses an embedding).
- FORMAT B (analytic-period shape): analytic big-period matrix — 4x8 complex
  balls (mid + radius per entry) at the declared dps, basis matching E and
  rho above. NOTE (C5): format-B O_K-action receipts are the conditional
  class O_K-ACTION-WITHIN-RADIUS (sec 1) — exact-action evidential claims
  need format A.
Precision policy (so the upstream side knows the ask): radii must admit (i)
Riemann + positivity receipts and (ii) Weil-box unique-integer isolation at
the pilot prime; otherwise abcount returns UNDECIDED-PRECISION carrying a
needed-dps estimate — that estimate IS the upstream re-ask,
machine-readable, no human-loop guessing.

## 7. Named design questions

Q1 tail-bound
certification authority; Q2 the O_K-action verification predicate; Q3 what
"isogeny-explicit" precisely requires (isogeny-class inputs vs one pinned
lattice); Q4 battery sufficiency (the honest
target-class gap, B5); Q5 the route-C precondition; Q6 the
global sign-convention pin. Q3/Q6 are resolved as binding pins in
S0PRE_RECEIPT.md.

## Files (all writes under the tool directory only)

M3T_REGISTRATION.md (this file), M3T_VALIDATION_BATTERY.md,
S0PRE_RECEIPT.md. Nothing here self-executes.
