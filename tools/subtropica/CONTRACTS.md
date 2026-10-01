# subtropica — PINNED INTERFACE CONTRACTS

Citations `SubTropica.wl:NNNN` refer to the VENDORED copy
`reference/SubTropica.wl` (upstream SHA adac2f72; see
reference/PROVENANCE.md — note the binary is git ead8c6e, adjacent commit;
behavior mismatches are resolved by LIVE probe of the binary and recorded
here). `main.cpp:NNNN` / `handlers.cpp:NNNN` refer to
`reference/hyperflint_cli_main.cpp` / `reference/hyperflint_bridge_handlers.cpp`.
Line numbers are verified against the vendored copies
(STSubtractionFormula 10723, STfastEpsSeries 17380, STGetFaces 10255,
STProduceUs 10587, STExpandIntegral 10906/11029, stHyperFlintBuildRequest
12357); the vendored numbers are authoritative.

====================================================================
## (a) HyperFLINT eval-json wire contract (ops we actually use)
====================================================================
Binary: your `hyperflint.sh` build (an LD_LIBRARY_PATH wrapper; READ-ONLY —
set `SUBTROPICA_HF_BIN` to its path). Worked request/response examples for
the ops below appear in `test/test_bridge.jl`; the convergent smoke
integrals referred to as T1/T2/T3 in the comments (T1: ∫ 1/((1+x)(2+x)) =
log 2; T2/T3: the ζ-valued `fixtures/synth_*` integrands) were verified to
60 digits against independent mpmath quadrature.

TRANSPORT RULES (main.cpp dispatch 3162-3244; wl-side
stHyperFlintBuildRequest at SubTropica.wl:12357):
- ONE flat JSON object per invocation on stdin. The C++ JSON parser is
  REGEX-BASED AND FLAT: never nest (beyond the documented array fields),
  never batch.
- Parse stdout's LAST line only; 3 diagnostic counter lines
  (`hf_rat_split_verify` etc.) go to stderr on every run and are normal.
- Exit code is 0 even on handled failure. NEVER gate on rc for op
  results. Errors are IN-BAND: `{"error":..}`, `{"failed":true}`,
  `{"divergent":true,..}`, `{"narrow_ctx_insufficient":true}` (retry ONCE
  with narrow-ctx env off, no recursion). Gate on presence of the result
  field AND absence of failed/divergent/error keys (print-loud-rc0-gates
  footgun class). Unknown op → `{"op":..,"error":"unknown op"}`.
- `"vars"` omitted ⇒ identifiers autoscanned; always pass it explicitly.
- Always pass `"mzv_data_path"` explicitly
  (the `data/mzv_reductions.json` of your SubTropica checkout).
- Env contract: PATH, HF_MSOLVE_PATH, HF_EULER_FILTER, OMP_NUM_THREADS,
  HF_NARROW_CTX, HF_LR_TIME_BUDGET_S (default 180), HF_MAX_THREADS_PER_CALL
  (RSS lever: 1 ≈ 970 MB/call). HF_LAZY_SUM=1 is REFUSED in combination
  with check_divergences. Every invocation under `ulimit -v 32505856`.
- Version gate (DESIGN.md RT-vi): stamp (`hyperflint.sh --version` →
  `HF_VERSION: 1.2.8`; response envelopes also carry `hf_version`) AND
  behavioral canary, live-verified:
  `{"op":"eval","a":"-x^2","vars":["x"],"values":["3"]}` → `"result":"-9"`.
  A source build stamps 1.2.0 by default — stamp alone is NOT a gate.
  NOTE: `{"op":"version"}` is NOT an op (live probe: `"error":"unknown
  op"`, rc=0) — do not use it.

OP: `hyperflint` (top integration driver; handlers.cpp
hyperflint_sym 1646-2235; request built per SubTropica.wl:12357):
```json
{"op":"hyperflint",
 "expr":"<Mma-grammar string, see (d)>",      // priority expr > f > wordlist
 "f":"<bare rational — fast path>",           // alternative to expr
 "vars_int":["x4","x3"],                      // ORDERED = THE integration order, each 0→∞
 "vars_int_from":[..],"vars_int_to":[..],     // optional finite intervals (rationals, Infinity/-Infinity/oo)
 "vars":["x4","x3","mm","s12","s23"],         // int vars + kinematic symbols
 "mzv_data_path":"<abs>",
 "check_divergences":true|false,              // POLICY below
 "algebraic_letters":true|false,              // FindRoots tier: mint Wm_i/Wp_i
 "carry_discharge":true|false,
 "canonical_emission":true}                   // deterministic bytes
```
Response:
```json
{"op":"hyperflint",
 "result":[{"coef":"<SymCoef Mma-string: rational fn of kinematics + Pi, I,
             Log[n], delta[var], mzv_* atoms, Wm_i/Wp_i/sqrt_disc_i>",
            "key":[["letter",...],...]}, ...],
 "vars":[AUGMENTED ctx: kinematics + ~700 MZV atom pool + Wm/Wp pool],
 "algebraic_letters":[{"idx":..,"poly":..,"var":..,"lc":..,"sum":..,
                       "product":..,"disc":..}, ...],   // when armed
 "timing_compute_s":s}
```
Semantics (SubTropica.wl:12930-13005): each term = coef · ∏_j
ZeroInfPeriod[word_j] (boundary periods commute; shuffle product =
ordinary product); empty key ⇒ pure constant/kinematic term. Response
`"vars"` is the AUGMENTED ctx — NEVER round-trip it as user vars.
`algebraic_letters:false` failure is LOUD `{"failed":true}` — escalate to
true (FindRoots tier), never silently zero roots.
check_divergences POLICY (verified SubTropica.wl:12374-12384; parity pinned
at wl:12589, the upstream HF-DIVCHECK-PARITY fix present in the pinned
commit): hard-ON for ALL standalone calls (with the flag off a divergent
input silently returns an empty result — confirmed on the binary; control
fixture `fixtures/divergent_xinv_1var.json`); OFF only for face-sum
counter-terms (individually divergent by construction, see (b)).
`{"divergent":true}` on a standalone call = FATAL, surface it.

OP: `evaluate_periods` — reduce residual words in the last integration
variable to MZV symbols:
```json
{"op":"evaluate_periods","regulator":[{"coef":"..","key":[[..],..]},..],
 "vars":[..],"mzv_data_path":"<abs>"}
```

OP: `zero_inf_period` — single word → MZV-symbol string.
LIVE v1.2.8 SEMANTICS (the
original `["0","-1"] → "1/2*mzv_2"` example was STALE vs the live binary;
per reference/PROVENANCE.md the live probe wins, and the build-time live
pins below are gated ≥60 dps against the independent ginac large-z
extrapolation in `zero_inf_period_ginac`, hlogexpr.jl):
  ["0","-1"]      → "mzv_2"
  ["0","-1","-1"] → "mzv_3"
  ["-1"]          → "0"
  ["-2"]          → "-Log2"
  ["0","-2"]      → "1/2*Log2^2+mzv_2"
  ["0","-1","-2"] → "-Log2*mzv_2+7/8*mzv_3"
Semantics: ZeroInfPeriod[{l1,...,ln}] = the shuffle-regularized constant
term (log z → 0) of G(l1,...,ln; z) as z → +∞, same letter order. Used
for residual-key evaluation; every such evaluation is cross-checked vs
ginac ≥60 dps (DESIGN.md gate A(v)).

OPS: `simplify_with_vieta`, `combine_wm_wp_ratios`, `back_substitute` —
the Wm/Wp → Vieta/Sqrt[disc] rewrite (DESIGN.md RT-v: use THESE, do not
reimplement wl:12930-13005 locally; only token→HlogExpr parsing is local).
CHAIN-ORDER PIN (live, hlogexpr.jl `normalize_wmwp`): run
`combine_wm_wp_ratios` BEFORE `back_substitute` — combine only rewrites
LITERAL `Wm_i/Wp_i` ratios to `WmOverWp_i`, while back_substitute smears
ratios into sqrt-disc-denominator shapes (an atom-in-denominator RATIO the
HlogExpr model cannot represent) that neither op can then recover. Default
chain: simplify_with_vieta → combine_wm_wp_ratios → back_substitute.
Wm/Wp letter state (AlgebraicLetterTable) is PER-PROCESS: cleared at each
`hyperflint` entry; `algebraic_letters_*` ops mutate it only within one
process — do not assume persistence across CLI invocations.

OP: `apply_mzv_reductions` — VARS-POOL RULE (live, assemble.jl
`_mzv_symbol_pool`): the op builds its PolyCtx from the request's `vars`
ONLY (no build_mzv_var_list — hyperflint_cli_main.cpp:2117ff), so a
reduction whose OUTPUT symbol is missing from `vars` dies with an in-band
Poly-parse error (live probe: f=mzv_2_1 with vars=[mzv_2_1] → error
"Poly: parse error: (-2*mzv_3)"; adding mzv_3 → result "-2*mzv_3").
CURE: always pass vars = the full table symbol pool (every lhs token +
every rhs symbol + the basis) ∪ the monomial's own tokens.

OP: `find_lr_orders` (LR-order search; schema_version 2 'carry-aware'):
```json
{"op":"find_lr_orders",
 "groups":[[<poly strings>],[..]],   // MULTI-GROUP per-ADDEND semantics —
                                     // MANDATORY for counter-term sums;
                                     // "polys":[..] single-group form exists
                                     // but a ct-sum via polys ⇒ NOLR by
                                     // construction (assert groups for sums)
 "xvars":[..],"coeff_vars":[..],
 "verify_order":["v1","v2",...]}     // OPTIONAL FIELD (handlers.cpp:741,
                                     // VERIFY-ORDER mode): check
                                     // THIS order is LR (O(n) st_fubini_lr
                                     // calls, no search). PIN CORRECTION:
                                     // there is NO standalone "verify_order"
                                     // op in the CLI dispatch — the source
                                     // map's op list was wrong; it is this
                                     // request field. Used for best_order
                                     // replay on resume.
```
Response: `{"best_order":[..],"score":..,"nolr":bool,"strategy":..,
"timing_compute_s":..}`. Time-budgeted branch-and-bound
(HF_LR_TIME_BUDGET_S): order may vary run-to-run; answer is
order-independent by theorem; persist + replay best_order (DESIGN.md R3).

OP: `find_lr_orders_scan` (Cheng-Wu gauge scan; no wrapper in this release):
`{"op":"find_lr_orders_scan","exps":[[[a,b],...]],...}` — `exps` twist
pairs REQUIRED; `keep_rule`: "Strict"|"FindRoots" (FindRoots tier is
NECESSARY-only/speculative); `euler_filter`, `max_orders`. Response:
`{"projective":..,"truncated":..,"orders":[{"order":..,"gauge":..,
"score":..,"carried_sqrts":..,"kin_sqrts":..,"terminal_quads":..},..]}`
score-ascending.

OP: `factor_table` — single-chain replay along an order,
factor-prediction table; loud guard errors, never truncation.

OP: `parse_expr` — round-trip diagnostic; the serializer's unit test
round-trips every emitted string through it (live: `-x^2` → canonical
`-x^2`).

Support ops used by the harness: `factor`, `eval`, `partial_fractions`,
`series_expansion`, `linear_factors`, `discriminant`, `resultant`.
Full op list (69 ops) with per-op schemas: inline in
`reference/hyperflint_cli_main.cpp` (grep `// Request:`; dispatch at
3162-3244).

HF CANNOT DO (front-end = OUR responsibilities): NO ε anywhere in
integration (per-ε-order Log-integrands prepared by us); NO
counterterm/subtraction construction; divergence analysis is
detection-only; residual ZeroInfPeriod keys not auto-evaluated; NO
physical-region continuation (Euclidean-positive assumed, emits
I·Pi·delta[var] residues); NO graph→Symanzik; NO final assembly.

====================================================================
## (b) Locally-finite integrand handoff (C8 → C9/C11)
====================================================================
Pinned from STSubtractionFormula (SubTropica.wl:10723-10850),
STfastEpsSeries (wl:17380), minOrder (wl:11188), directory-layout writes
(wl:17510-17532, 17788, 17821), flat-dx comment (wl:~10742):

Per face σ ∈ Σ_div, per ε-order o ∈ [minOrder..maxOrder], per counter-term
n: ONE integrand =
  rational function of the face variables (vars ∉ J gauged to 1 — J as in
  WARNING-2 below) and kinematic coefficients
  × integer powers of Log[P_j].
Measure = FLAT dx_i on [0,∞)^k — NOT dlog; the /∏x_J division is already
applied. minOrder = |vars(last face)| − |vars(first face)| (wl:11188);
maxOrder = requested order + prefactor pole offset. Each ct individually
may be log-divergent at boundaries; divergences cancel ONLY in the face
sum ⇒ the integrator must shuffle-regularize ⇒ `check_divergences` OFF
per-ct inside the farm, hard-ON everywhere else.
The in-code representation is `LogIntegrand` (src/types.jl). Jacobian
factors (1−u_j)^(1−TropI(ρ_j)) carry ε-DEPENDENT exponents (TropI with
regulators) — they enter the poly list as (poly, EpsExp) pairs like any
P_j (DESIGN.md RT-extra-jac; unit fixture mandatory).

Persisted layout (mirrors .wl for auditability), under
`<run-dir>/integrands/<id>/ord_<o>_face_<i>/`:
`polys.json` (LR letter groups per ct — poly bases with ε-dependent or
negative exponents), `vars.json`, `best_order.json`, `ct_<n>.json`,
`partial_results/result_ct_<n>.json`.
Assembly = per-face per-order ct sums × NP prefactors × Normalization
(e^(L·γ_E·ε) type), ε-Series to output order.

====================================================================
## (c) CORRECTED PINS — WARNINGS (wrong reading spelled out)
====================================================================
**WARNING-1 — counter-term NET sign is (−1)^|subface|.**
Source: STSubtractionFormula, SubTropica.wl:10723ff. The code builds
`Times[-a[subface], ...]` where `a[s] = (−1)^(|s|+1)`; the NET assembled
sign is therefore −(−1)^(|s|+1) = **(−1)^|s|**: empty subface → +1 (the
identity term), single facet → −1 (subtracted, as it must be).
WRONG READING: "sign (−1)^(|subface|+1)" — reading `a[s]` alone
as the term sign. A transliteration of the wrong reading ADDS single-facet
counter-terms instead of subtracting them. Detection: Gate B(iii)
pole-cancellation probe — which MUST run before any production use.

**WARNING-2 — I vs J: the completion/division/gauge subset is the
COMPLEMENT.**
Source: STSubtractionFormula, SubTropica.wl:10723ff (eqs 3.34-3.35).
Correct: **I** = the FIRST index subset with det(rays[fc, I]) ≠ 0 (the
det≠0 ray-support subset); **J = Complement(1:n, I)**. The
Kronecker-delta completion rows in Vol = |det(rays|e_J)|/∏(−TropI) are
indexed by **J**; the output division /∏x_J uses **J**; the gauge-fix
sets vars ∉ J to 1 — all three use J, the complement.
WRONG READING: "J-subset = first det≠0 basis completion" —
i.e. J = I. This yields FINITE, numerically-plausible, wrong-volume,
wrong-gauge results per counter-term — silent until the face-sum
comparator. Any implementation of C7 must carry a unit test that
distinguishes I from J on an example where I ≠ {1..|I|}.

(Related conventions pinned while we are here: u-functions stored as
u = 1 − 1/(1+x^w), Factor-ed — the CODE convention (SubTropica.wl:10636,
10669), NOT the paper's v-function of eq 3.28; w normalized by
w·ρ_f = −1; re-derive once against the worked eq 3.37 example. "NotFound"
⇒ geometric property violated ⇒ B1 LOUD refusal / B2 NP
continuation.)

====================================================================
## (d) HlogExpr serialization grammar
====================================================================
The serializer (`src/serialize.jl`, unit-1) emits Mma-grammar STRINGS —
the exact input grammar of HF `parse_expr`
(include/hyperflint/convert/parse.hpp; vendored main.cpp inline docs):

- Atoms: integers, rationals `p/q` (CONTENT-CANONICAL: gcd(p,q)=1, q>0,
  sign on numerator — Rat::parse does NOT normalize non-coprime content,
  handlers.cpp T3 caveat; feeding non-canonical rationals is a silent
  wrongness source), symbols `[A-Za-z][A-Za-z0-9_]*`. NO backticks/
  contexts ever (Mma side strips them; we never emit them).
- Operators: `+ - * / ^` with Mathematica precedence; `^` binds tighter
  than unary minus (`-x^2` ≡ `-(x^2)`; guaranteed by the v1.2.8 canary,
  see (a)); integer exponents may be negative; explicit `*` always
  (never juxtaposition).
- Functions: `Log[expr]` (rewritten internally to Hlog[arg,{0}]);
  `Hlog[z,{l1,l2,...}]` with letters as rational-function strings;
  `PolyLog` accepted-but-flagged upstream — we never emit it.
- Special symbols (parse side, from HF responses): `Pi`, `I`, `Log[n]`
  (incl. `Log2` token = Log[2]), `delta[var]` (contour residues),
  `mzv_a_b_c` (ζ(a,b,c); `m` prefix on an index = negative/alternating,
  e.g. `mzv_1_m3`), `Wm_i`/`Wp_i`/`WmOverWp_i`/`sqrt_disc_i` (algebraic
  letters, resolved via the response's algebraic_letters table
  {idx,poly,var,lc,sum,product,disc}). Our own symbolic atoms
  (`EulerGamma`, zeta tokens from loggamma_series) exist only INSIDE
  HlogExpr — they are never sent to HF.
- Round-trip law: every string the serializer emits MUST round-trip
  through `{"op":"parse_expr"}` unchanged-modulo-canonicalization; this
  is a standing unit test (unit-1, test_bridge.jl).
- Parsing direction (HF coef strings → HlogExpr): tokens map to
  `src/types.jl` atoms — mzv_* → HZeta (m→negative index), Log[n] →
  HLogA, delta[v] → HDelta, Wm_i/Wp_i/sqrt_disc_i → HAlgLetter (+ table
  entry), Hlog[z,{w}] → HHlog, residual keys → HPeriod; everything
  rational stays in the term's canonical coef string.
  CANONICAL Log2 ATOM :
  log 2 has exactly ONE atom spelling, `HConst(:Log2)` — both the `Log2`
  token and `Log[2]` parse to it, and `HLogA("2")` folds to it at every
  canonicalization chokepoint (`_norm_atom`, hlogexpr.jl), so the two
  spellings can never coexist in one term list (they previously did not
  merge in term arithmetic — unit-3 hazard). The serializer emits
  `HConst(:Log2)` toward HF as `Log[2]` (grammar-safe input form).
- Canonical text form of HlogExpr (fixtures, persisted artifacts): sum of
  terms, each `coef * atom1^k1 * atom2^k2 * ...` with atoms in the fixed
  order (HConst by name, HZeta by index vector, HLogA by arg string,
  HHlog, HPeriod, HAlgLetter by (kind,idx), HDelta by var), terms sorted
  lexicographically by their atom part; deterministic bytes end-to-end
  (pairs with HF's `canonical_emission:true`).

--------------------------------------------------------------------
### §d addendum — Wm/Wp algebraic letters: engine-grammar facts + transport law
--------------------------------------------------------------------
Three facts about the engine's coefficient grammar, each established by a
raw request/response exchange with the v1.2.8 binary, govern how
coefficients carrying algebraic-letter tokens are transported into the
RT-v ops:

1. **The letter's own variable must be in `vars`.** The RT-v ops'
   allocation replay hard-fails when the letter's OWN variable is missing
   from the request `vars`: the handler prints
   `simplify_with_vieta: allocation var 'x2' not in PolyCtx` to STDERR and
   returns rc=1 with EMPTY stdout (the only known rc≠0/empty-stdout path
   in the op surface; the bridge types it `HFTransportError`). RULE (in
   `normalize_wmwp`): the engine-side vars ALWAYS include each
   AlgLetter's own `var`, even when the coefficient is constant in the
   letters.

2. **The divergence pre-check false-positives on convergent Wm/Wp-bearing
   integrands.** The per-log-atom rational pre-check (DESIGN_B1 §4b)
   reports `Log[var]^2 / var^0 at var=x1 (infinity boundary)` at every
   integration order while the counterterm tails cancel across DIFFERENT
   log atoms (Log[P1] vs Log[2x1+x2+1]) and independent dps-55 tanh-sinh
   quadrature converges. Such an integrand needs the stamped checks-OFF
   retry of the tiered driver; it is not an engine op-chain limitation.

3. **Rat::parse ≠ the coefficient emission grammar.** (i) The engine
   IN-BAND-REFUSES its own parenthesized ratio emission —
   `{"op":"simplify_with_vieta","expr":"(2/(Wm_1 - Wp_1))",...}` →
   `error "Poly: parse error"`; (ii) WORSE, it SILENTLY MIS-PARSES
   trailing factors after a division — `"N/(D)*mzv_2"` parses as
   `N/(D*mzv_2)` (mzv_2 lands in the DENOMINATOR; Mathematica precedence
   says `(N/D)*mzv_2`). Whole-coefficient-string transport into the RT-v
   ops is therefore FORBIDDEN except for the `_ratparse_safe` class:
   '/'-free polynomial bodies, or exactly ONE depth-0 '/' whose
   denominator factor reaches end-of-string, no depth-0 +/- (leading
   unary minus excepted). (iii) The chain NEVER substitutes `disc`:
   outputs carry `sqrt_disc_i^2` (and ^3, ^4, ...) verbatim.

**TRANSPORT LAW.** `hf_integrate`'s coefficient path: `parse_coef`
first; a SubTropicaParseError on a coefficient CARRYING algebraic-letter
tokens routes through `normalize_wmwp`; the refusal stays typed only if
the NORMALIZED form still fails. `normalize_wmwp` (RT-v chain order:
simplify_with_vieta → combine_wm_wp_ratios → back_substitute):
  * leg 1 — whole-string chain, `_ratparse_safe` shapes only;
  * leg 2 — structural routing: the coefficient is parsed locally with
    the chain armed as the DENOMINATOR ROUTER (`_PCtx.wmwp_router`): each
    irreducible multi-term Wm/Wp denominator D is emitted '/'-free
    (kinematic + rational content cleared exactly) and routed as the
    single-division body `"1/(D)"`; per-term Wm/Wp sub-monomials
    (Vieta products, lone roots) get the same chain pass, cached;
  * local exact closures (NOT a Vieta reimplementation — table-driven
    definitions, mirroring upstream's Mathematica auto-simplification of
    the substituted Sqrt[disc] literal): `sqrt_disc_i^2 → disc` fold at
    every parse exit (floor-mod: `1/sqrt_disc → sqrt_disc/disc`), and
    conjugate rationalization `1/(a+b·sqrt_disc) = (a−b·sqrt_disc)/
    (a²−b²·disc)`.
Letters-tier boundary KEYS can carry Wm_i/Wp_i letters (eq 4.11 k0);
`zero_inf_period_value` refuses them typed under strict (an HPeriod atom
under strict=false) — the rational closed-form branches never see atoms.
End-to-end acceptance (massive-bubble √5 fixture + eq 4.11 k0/R'1 vs
independent dps-55 tanh-sinh quadrature, ≥30 floored digits + mutation
control): test/test_wmwp_integration.jl. PRECISION NOTE: k0/R'1 match
the dps-55 quadrature values at 30-31 digits because that is the
quadrature's own error floor; a dps-70 quadrature agrees with the exact
symbolic values to 39/36 digits (k0/R'1), so the symbolic route is the
better value on both.

--------------------------------------------------------------------
### §d addendum 2 — word-level Wm/Wp periods
### (deeper gap layer) + canonical-atom law resolution
--------------------------------------------------------------------
1. **Word-level algebraic-letter periods (numeric-leg law).** The letters
   tier can mint Wm_i/Wp_i tokens INSIDE ZeroInfPeriod WORDS (live: the
   eq 4.11 driver route, word ["-1","0","Wm_1"] over the Gaussian
   quadratic 2x1²+2x1+1) — beyond the §d-addendum-1 COEF-side transport
   law. These periods are algebraic-NUMBER-valued and generally COMPLEX
   per word; there is no exact atom form for a complex-letter G in the
   HlogExpr model. LAW (wired):
   * `zero_inf_period_value` keeps the TYPED refusal under strict=true;
     strict=false keeps the explicit HPeriod atom (never silent) — the
     symbolic result then carries the atoms honestly (NOT closed-form).
   * `eval_symbolic` (verify.jl) owns the NUMERIC leg: HPeriod words with
     rational letters → `zero_inf_period_ginac`; with lone Wm_i/Wp_i root
     tokens → `zero_inf_period_alg` (hlogexpr.jl — same large-z Lagrange
     extrapolation at complex letters via ginac_gpl's re/im argument
     pairs; conjugate-symmetry cache; truncation floor ~1e-110);
     HAlgLetter atoms → `alg_letter_roots` (MONIC law
     Wm/Wp = (sum ∓ √disc)/2, root-checked against x² − sum·x + product;
     typed refusal on lc ≠ 1 or stub tables). Both are complex per TERM;
     the real-line law is enforced per COEFFICIENT SUM: residual
     imaginary part above 10^-(dps−8)·max(|re|,1) is a typed refusal —
     never a silent real().
   * The reference letter table is `LaurentSeries.letters` (response
     table, disc populated). `_letters_of` stubs never feed the leg.
2. **Canonical-atom law.** `fold_log_atoms`, its helpers and
   `SubTropicaSerializeError` live in hlogexpr.jl; `canonical_text` folds
   at ENTRY, so BOTH emission chokepoints (`canonical_text`, `mma_string`)
   emit identically-zero Log[rational] combinations as the literal "0"
   (checked on the eq 4.3 ε⁻³ coefficient, which evaluates to "0", with a
   held mutation control). Standalone contexts including hlogexpr.jl
   without serialize.jl carry the same law.
3. **Evaluation of rational-letter HPeriod atoms (test-pinned):**
   eval_symbolic on a rational-letter HPeriod atom EVALUATES it through
   the independent ginac route rather than refusing, so the strict=false
   explicit-atom path is verifiable end-to-end. Still-refused classes
   (typed, pinned in test_verify.jl): algebraic-letter words/atoms with NO
   table entry, non-monic letters, compound word letters, HDelta,
   HConst(:I).
