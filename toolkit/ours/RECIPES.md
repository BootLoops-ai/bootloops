# Recipes

Alongside the packaged tools, the toolkit carries *recipes*: working patterns proven
in real calculations and written down with the acceptance checks and pitfalls that
made them reliable. Each entry below says what the recipe does, when to reach for it,
what must be verified before trusting a result, and where the implementation lives.
Several recipes are backed by packages shipped in this repository under
[`../../tools/`](../../tools/README.md); where no package ships, the entry says so —
the recipe text is then the reference, detailed enough to reimplement.

Contents:

- [Spectral decomposition of tabulated kernels (the FFTLog route)](#spectral-decomposition-of-tabulated-kernels-the-fftlog-route)
- [Sizing spectral/FFTLog farms: the amplification pre-check](#sizing-spectralfftlog-farms-the-amplification-pre-check)
- [Two-fold integrals: conic direction goes outside](#two-fold-integrals-conic-direction-goes-outside)
- [The five-stage close chain](#the-five-stage-close-chain)
- [Weight-graded PSLQ closure discipline](#weight-graded-pslq-closure-discipline)
- [Point-count shortcuts: one charpoly call replaces brute enumeration](#point-count-shortcuts-one-charpoly-call-replaces-brute-enumeration)
- [Boundary constants from the operator's arithmetic mod p](#boundary-constants-from-the-operators-arithmetic-mod-p)
- [Spectral densities by DE transport](#spectral-densities-by-de-transport)

---

## Spectral decomposition of tabulated kernels (the FFTLog route)

**What it does:** turns an integral over a tabulated function into a linear
combination of analytic loop-integral masters, so the whole amplitude-class machinery
(IBP, closed forms, certified evaluation) applies to data-driven integrands.

**Use when** your integral convolves a tabulated or measured function (a linear power
spectrum, a structure function, a spectral density) with rational or analytic
kernels — e.g. cosmology loop integrals over P_lin, or any loop-like integral whose
input is a table rather than a formula.

**Route:** expand the table in complex power laws via FFTLog on a log-spaced grid
(~20 lines of numpy: FFT of the log-sampled table → coefficients c_m with complex
exponents ν_m, bias chosen inside the table's convergence strip). Each power-law term
turns the integral into an analytic amplitude-class master; reassemble the answer as
the linear combination.

**Checks (mandatory):** BEFORE sizing any master-table farm from this decomposition,
run the [amplification pre-check](#sizing-spectralfftlog-farms-the-amplification-pre-check)
below — four fail-closed checks. NEVER size a farm from window self-convergence
alone: in the calculation this recipe comes from, two independent failure layers were
shared systematics invisible to internal convergence ladders — the internal spread
was 2.6e-6 while the result was wrong by 9.5e-2. On any check failure, evaluate the
TOTAL kernel directly in momentum space against the real table: subset remainders
depend on continuation conventions, and only the total has convergent tails.

**Measured performance** (two-loop matter power spectrum): the P33-II contribution by
direct total in 3.3 core-hours, P42 in 7.19 core-hours; the FFTLog master tables were
retained as mid-band cross-checks and agree with the direct route at the 5.4e-5
level.

**Implementation:** the FFTLog expansion itself is the ~20-line numpy sketch above;
no package ships for the assembly pipeline. The recipe and its pre-check are the
reference.

---

## Sizing spectral/FFTLog farms: the amplification pre-check

**What it does:** a cheap probe run before committing to any spectral/FFTLog
master-table farm. Measure the cancellation amplification
A(k) = Σ|terms| / |result| at 2–3 PRODUCTION k-points using cheap float64 masters,
and apply four fail-closed checks. If the probe fails, do not build the table:
evaluate the total kernel directly in momentum space against the real input table
(the "direct-total" pattern) instead of assembling per-monomial master tables.

**Use when:**

- you are about to size a spectral/FFTLog master-table farm from window
  self-convergence or any other proxy — run this probe FIRST;
- an FFTLog/spectral master-table assembly is losing digits, a squared-kernel block
  goes negative, or self-convergence looks clean while a cross-method check fails.

**The four checks:**

1. **Bar** — A(k_max) × the measured window-truncation error must be at or below the
   precision the downstream consumer needs.
2. **Growth** — A(k) must not grow with k.
3. **Positivity** — squared-kernel blocks must be positive at every probe point.
4. **Anchor** — at least one probe point must agree with a
   representation-independent evaluation.

Any failure means: go direct-total, and demote the table route to a cross-check.
Accounting rule: bounds decide pass/fail, residuals are what you report.

**Scope limits:** a table route whose A(k) grows with k should never be the shipped
route — it survives only as a mid-band cross-check. This entry describes the probe
and the decision rule; the assembly-side positivity checker used alongside it is not
included in this repository.

**Pitfalls:**

- Subset remainders depend on the analytic-continuation convention — in the
  originating calculation, off-strip continuation choices moved one subset remainder
  by factors of 26–330. Only the total kernel has convergent tails; always evaluate
  the total.

**Implementation:** not included in this repository. The probe is a few dozen lines
over your own master table; the checks above are the specification.

---

## Two-fold integrals: conic direction goes outside

**What it does:** a Fubini ordering rule for two-fold integrals with one elliptic and
one conic direction: put the conic direction on the OUTSIDE, and use a ratio-form log
primitive for the inner integral. The worked evaluator classifies the integrand F
into two classes:

- **F_jump** — ratio-logs log((x−1/2)/(z−1/2)) and bare log(x−1/2): these jump
  across x = 1/2 but are pole-free, so a cheap seam-split handles them;
- **F_rest** — rational terms, pole × regular-log, and polylogs: continuous along
  the +iδ path, so endpoint-only evaluations suffice.

**Use when** a two-fold integral is elliptic in the inner variable and conic in the
other.

**Pitfalls:**

- δ MUST stay ≪ s: taking δ ≫ s kills the seam jump silently.
- Working precision must be ≥ eexp+90 digits — there is a hidden cancellation depth
  of eexp+40 (measured); add adaptive extra digits near m = 1/2 (distance < 0.05).
- `sp.expand(F)` before classification — composite `Mul`s must be split so pole
  pieces separate from jump logs.
- The ratio-form log primitive is MANDATORY, not an optimization.

**Implementation:** not included in this repository; the worked evaluator exists only
inside the calculation it was built for, and the classification rules above are the
reference. (Not to be confused with the third-kind conic connection rows in
[`tools/wayfinder`](../../tools/wayfinder/GUIDE.md) — that is a different object
that happens to share the word "conic".)

---

## The five-stage close chain

**What it does:** takes a fresh SIMPLE target from graph spec to a certified,
ready-to-close bundle through five staged steps, each backed by a shipped tool.

**Use when** a new integral family fits the classic
alphabet → oracle → held-out-certify pattern and you want the stage order and
its gates pinned before any fitting happens.

**Route:** (1) ALPHABET — Landau letters from the graph spec
([`tools/landau-alphabet`](../../tools/landau-alphabet/)); (2) ORACLE —
high-precision values at rational kinematic points, farmed to a staggered digit
goal (`solve_integrals` mode of `amflow_cli`, built in the `amflow-cpp` checkout —
see [`upgrades/ENGINES.md`](../../upgrades/ENGINES.md)); (3) INTEGRITY — hash-dedup +
contamination scan of the oracle dumps before anything is fitted
([`tools/gatekeeper`](../../tools/gatekeeper/) `scan_integrity`); (4) BASIS —
basis-value columns at the sample points
([`tools/gpl-eval`](../../tools/gpl-eval/)); (5) CV — held-out certification
against oracle points never used in the fit (gatekeeper `heldout_certify`),
summarized in a READY verdict.

**Checks (mandatory):** stage order is the contract — integrity gates BEFORE
fitting, CV gates AFTER, and no stage consumes an ungated predecessor.

**Implementation:** no single package ships for the chain; the five stages are
the shipped tools above, invoked directly, and this entry is the reference for
the order and the gates.

---

## Weight-graded PSLQ closure discipline

**What it does:** a weight-graded PSLQ closure and single-valued (sv) adjudication
ladder for closing MZV-ring Laurent coefficients at scale — the engine behind a
49-function modular-graph-function atlas (results from the four-edge MGF paper,
bootloops.ai/diagrams/mgf.html). It builds stuffle-checked constant baskets,
closes every slot with a basis-invariant verdict, escalates the nulls in height and
precision, adjudicates sv-alignment, and measures quotient ranks against
Broadhurst–Kreimer dimensions.

**The ladder's components** (module names as used in the original implementation):

- **Per-slot closure with an honest verdict taxonomy** —
  `close_slot(coeff_str, tw, build_dps, basket_dps)` returns one of: ZERO /
  CLOSED_RATIONAL / CLOSED_PRODUCTS / CLOSED_DEPTH3 / CLOSED_WITH_EVEN (an sv
  violation; filed identically — a negative result is the bigger result) /
  NULL_BOUNDED (with pool, height, and precision legs stated) / NONZERO_EMPTY_POOL
  (an out-of-ring candidate). The ladder runs L1 (products) → L2 (+depth-3) → L3
  (+even probes).
- **Weight-graded baskets with receipts** — `build_basket(tw, dps)`: odd-zeta
  monomials + all-odd depth-3 MZVs (Richardson single-sum evaluation) + ζ₂^k even
  probes. Every depth-3 value is stuffle-checked at build time (two identities per
  multiset, agreeing to ≥ dps−15 digits); baskets are cached per (tw, dps).
- **Null escalation** — a PSLQ null at weight ≥ 5 is a HEIGHT failure first, not a
  span failure. Escalate maxcoeff 1e8 → 1e9, with a leading-rational rung to 1e18
  (this closed one coefficient as c₁₀ = 676694/974482387003125 with 204.6 held-out
  digits), then escalate precision and legs: dps = 560, legs at 430/530 digits,
  heights ≤ 1e12.
- **sv adjudication** — PSLQ against an sv-aligned basis (verdicts SV_CONFIRMED /
  SV_VIOLATION / UNRESOLVED), or — cheaper and exact — ALGEBRAIC sv-adjudication:
  invert the sv-formula matrix, Euler-reduce each zB(a,b) via tiny two-precision
  PSLQs (n ≤ 7, ~400 digits, cached per weight), and take the Fraction-exact
  complement — a verdict with no big scan at all.
- **Quotient ranks** (dimension denominators for values mod products) — a
  sector-restricted pool, two independent lindep legs, and a full-precision residual
  arbiter (see checks below).

**Checks (two independent engines throughout):**

- Per-slot closure inherits the full PSLQ discipline: two-precision
  canonicalize-before-compare with identical vectors, a capacity refusal (the scan
  refuses when the pool cannot resolve the requested height), a synthetic positive
  control planted IN the pool at ≥10× below maxcoeff, r[0] == 0 treated as
  DEGENERATE (never as a null), and reverification at build precision beyond the top
  leg.
- Basket values are checked against an independent FORM-based MZV engine plus
  stuffle identities; the sv adjudication was validated 21/21 against the Dorigoni
  et al. ancillary files.
- Rank law for the quotient: a candidate is DEPENDENT iff either lindep leg's vector
  passes a ~700-digit residual arbiter (> 652 digits of cancellation). Result on the
  atlas: 63 relations, residuals 695.5–700.2 digits, zero spurious and zero internal;
  independence certified to heights ≤ 1e15 in the ζ₂-free sector. The occupancy
  denominator matches the Broadhurst–Kreimer curve 1, 2, 2, 4, 5, 6 at
  w = 11..21 — the naive Lyndon depth-3 curve is WRONG from w = 15, and sv-alignment
  ranks are a spanning artifact, never a denominator.

**Honest caveats:** the ℓ=5 slots carry no independent oracle (a gap inherited
from the originating calculation). The sv formulas at w = 15–21 are this project's own derivation
(Brown's §5 fixed-point recipe, certified by a 1e−685 residual), beyond the
tabulated weights we could find in the literature (see
bootloops.ai/diagrams/mgf.html).

**Pitfalls:**

- Depth-3 + even-probe pools are ℚ-DEPENDENT by construction (stuffle relations) —
  take the independent-subset quotient FIRST, else PSLQ returns r[0] = 0 internal
  relations instead of the closure (observed in practice).
- The independence filter is a basket-level fact: persist it to disk keyed by
  (tw, dps, tier). Process-local caching made every farm worker repay a ~1-hour
  dimension-40 greedy cascade.
- Never put the constant 1 in a pool of weight > 0; single-member pools close via a
  direct-ratio check, not the scan.
- sv relation heights grow roughly ×100 per weight step (1.4e3 → 7.8e12 from w = 11
  to 21) — size maxcoeff ladders accordingly. Two coexisting relations make the two
  lindep legs return DIFFERENT genuine vectors, so leg agreement alone is not a
  criterion; the residual arbiter decides.
- The ℓ=5 even-weight ladder pools products + depth-2 + depth-4 + even probes, but
  NOT the ζ₂-tower × (ζ_odd · depth-3) product class (e.g. ζ₂^n · ζ₃ · ζ(5,3,3)).
  This is a pool-design gap, not a capacity wall: slots whose exact closure needs
  exactly that class return NULL_BOUNDED. Symptom: NULL_BOUNDED at slots the exact
  symbolic route closes as depth-3 mixed. Fix: extend the pools with that class, or
  route persistent nulls to the exact symbolic close.

**Implementation:** the ladder modules themselves are not included in this
repository. The PSLQ layer they are built on ships as
[`tools/lockpick`](../../tools/lockpick/GUIDE.md) — `pslq_gate` is the consolidated
closure harness (two-precision stability, capacity refusal, planted controls) and
`mplll` the lattice value-fitter; use those rather than raw `mp.pslq` when rebuilding
the ladder.

**Credit:** Broadhurst–Kreimer dimensions; single-valued MZVs and the f-alphabet —
F. Brown; sv formulas cross-validated against the ancillary files of Dorigoni et al.
(arXiv:2403.14816).

---

## Point-count shortcuts: one charpoly call replaces brute enumeration

**What it does:** replaces O(p⁴) brute-force point-count fingerprints with p-adic
computations:

- **genus-2 side:** the Weil pair (a_P, b_P) of a genus-2 quartic from a single
  PARI `hyperellcharpoly` call, giving N₄ = q² + 1 − (a² − 2b);
- **K3 side:** the Watson (2,2,2) K3 count fibered into p² genus-1 quartic traces
  via `ellap` on the classical-invariant Jacobian Y² = X³ − 27·I·X − 27·J.

Complexity drops from O(p⁴) to O(p² polylog) for the two-parameter side and O(p) for
the one-parameter side — measured ~400× at p = 431, and growing as p².

**Use when** an arithmetic point-count fingerprint battery is walled at O(p⁴) on K3
or genus-2 sides (curve #C(F_{p⁴}) via the norm trick, or a brute affine K3 count),
or when extending inertness certificates to fresh large primes.

**Checks:** exact-integer match at 17 reference primes (11 inert primes 7..157 and 6
split primes 23..131, both sides, every variant); cross-check of engine A
(`hyperellcharpoly`) against engine B (`ellap`); a degenerate-fiber census; the Weil
functional equation asserted on every charpoly; and two fresh inert end-to-end
certificates at p = 431 and 461 (including a dent control) — 15.4 s and 15.3 s
against a ~1.7 h projection for the brute route.

**Scope limits:** keep an independent brute-force numpy leg as the non-circular
check — do not replace it with the fast engines it is meant to check. This is not a
Costa–Harvey average-polynomial-time port; it is a direct PARI/GP pattern.

**Pitfalls:**

- PARI `hyperellcharpoly` over NON-PRIME fields (t_FFELT) stack-overflows from
  p = 37 even at 4 GB parisizemax — use `ellap` on the Jacobian as the fiber engine,
  and keep `hyperellcharpoly` cross-checks to k = 1 or p ≤ 29.
- `gp -q FILE` drops into the REPL after the script (hangs pipelines) — append
  `quit` or feed the script via stdin.
- Multi-line `for(...)` in a `.gp` file needs `{}` wrapping (a newline is a
  statement end outside braces).

**Implementation:** not included in this repository — the recipe rides directly on
PARI/GP (`hyperellcharpoly`, `ellap`) and the checks above are the specification.
For certified counting on abelian FOURFOLDS (a different object), see
[`tools/abacus`](../../tools/abacus/GUIDE.md).

**Credit:** PARI/GP (`hyperellcharpoly`, `ellap`).

---

## Boundary constants from the operator's arithmetic mod p

**What it does:** extracts the boundary constant of a class-C (Eichler/elliptic)
sector from the ARITHMETIC of its Picard–Fuchs operator alone — no high-precision
integral evaluation, no auxiliary-mass flow, no master solve. Reduce the operator
mod p per singular fiber, compute the p-adic Frobenius trace (= the Hecke eigenvalue
a_p mod p of the associated newform), rebuild {a_n}, sum the L-series, and PSLQ the
L-value to closed form.

**Use when:**

- a class-C sector needs its DE-transport starting value and you want it without any
  high-precision integral evaluation;
- the operator is already in hand from GeoTriage — the Frobenius-trace step
  then runs essentially for free.

**Scope limits:** on the K3 half-integer rung the method is honestly PARTIAL: the K3
period is Gevrey-0, the target there is an Eichler cusp EXPANSION rather than one
constant, and three singular fibers make light Stokes extraction path-dependent —
DE transport remains mandatory, and the method reports the limitation rather than a
number. It is also not a generic point counter (GeoTriage owns that role).

**Checks:** demand a 100% trace-vs-known-eigenvalue match across the declared
prime × fiber grid BEFORE trusting any reconstruction; run PSLQ with two-precision
stability and a positive control; hold out digits against an independent evaluation
where one exists (the K3 banana closure was checked to 70 digits against a constant
obtained this way).

**Pitfalls:**

- Do not conflate this with the var=∞ Frobenius-BRANCH boundary machinery
  (exponent-spectrum/branch analysis) — a different object that happens to share the
  word "Frobenius"; that machinery ships separately as
  [`tools/frobenius-boundary`](../../tools/frobenius-boundary/GUIDE.md).
- The a_p background subtraction on toric point counts is polynomial in p — remove
  it before reading off the centred a_p. The shipped classifier does this; any
  hand-rolled count must too.

**Note on feasibility:** `pari` `lfun` at LARGE conductor is CHEAP — L(E,2) at
conductor N = 23941299360 (2.39e10) runs in 28.7 s at 135 digits and 33.5 s at 175
digits, with 134.7 digits of cross-precision agreement. Do not assume the L-series
leg is infeasible at large level.

**Implementation:**
[`tools/geotriage`](../../tools/geotriage/GUIDE.md) is the reference
implementation and code anchor of this recipe — its classifier performs the
GeoTriage, the CM/a_p test, and the background subtraction. The one-off probe
scripts from the original derivation are not included. For the PSLQ leg, use
[`tools/lockpick`](../../tools/lockpick/GUIDE.md).

**Credit:** the Frobenius/point-counting mechanism follows A. Lauder (Found. Comput.
Math. 4 (2004) 221) and K. Kedlaya (arXiv:math/0105031); Chowla–Selberg; PARI/GP.
The operator → validated-boundary-constant assembly is this project's.

---

## Spectral densities by DE transport

**What it does:** evaluates a spectral density ρ(w′) = Im Σ(w′) to ~100 digits at
many points from ONE boundary value and one exact differential equation, instead of
per-point high-precision solves. The worked instance is the kite self-energy: an
exact-rational 9×9 (working 8×8) DE in w, ε-expanded to orders A0/A1/A2, walked as a
24-dimensional Laurent-block Taylor march from a single AMFlow boundary at w = 5.

**Use when:**

- ρ(w′) via parametric quadrature is slow or imprecise across many quadrature nodes;
- you need high-precision ρ at many w′ and would otherwise run a per-point solve at
  each node.

**Scope limits:** the worked instance is hardwired to the kite family — 9 masters
(master 4 dropped as a decoupled sink), singularities at w ∈ {0, 1, 9}, boundary at
w = 5. It is a pattern to re-instantiate per family, not a general spectral-density
engine.

**Checks:** demand ≥80 digits of agreement against independent per-point oracle
solves. The worked instance passed at 99.95d (w = 12), 99.90d (w = 50), and 99.89d
(w = 100); minimum checkpoint 99.89d.

**Pitfalls:**

- Arc direction around w = 9: Im(w) > 0 (the Feynman s+i0 prescription) is correct;
  an arc BELOW lands on the complex-conjugate sunrise branch (only 1.1 digits of
  agreement at w = 12). Im(w) > 0 is the validated direction — trust it over any
  docstring that names a negative-imaginary semicircle.
- Laurent-fit reference values can be self-stable to 134–217 digits and still WRONG
  at ~30 digits (an η-truncation artifact in the ε⁰ order was measured doing exactly
  this). Fit stability is not correctness — never use such fits as oracles.

**Implementation:** the kite-specific DE system and boundary are not included in this
repository. The general transport machinery the recipe runs on ships as
[`tools/wayfinder`](../../tools/wayfinder/GUIDE.md) (exact-rational DE loading,
fixed-ε Taylor march, Frobenius landings, ε-Laurent extraction, two-precision
gates), and the dispersion-quadrature consumer it pairs with ships as the
`nestor.dispersion` member of [`tools/nestor`](../../tools/nestor/GUIDE.md), whose
guide records the paired pattern (1669 nodes at 99.9 digits, cross-validated at 3
points).

---

[← back to the toolkit README](README.md)
