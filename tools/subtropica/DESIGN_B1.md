# subtropica — DESIGN_B1 (log-divergent subtraction + continuation routes)

Companion to `DESIGN.md` for the divergent routes: B1 (tropical
subtraction, C4/C5/C7) and B2 (Nilsson–Passare continuation, C6). Source
comments cite this file by section number (`DESIGN_B1 §2`, `§3`, `§4`,
`§4b`, `§4c`); the numbering below is stable for that reason.

Precedence where documents overlap: `CONTRACTS.md` (wire contracts and the
two WARNING pins, §c) first, then `DESIGN.md` (components, hazard
register), then this file (route map and gates). Source comments also
cite `SPEC_B1 §n`: those are section numbers of the detailed,
transliteration-grade specification the B1 modules were written from
(§2 polytope layer, §3 tropical layer and its ordering law, §4 the Möbius
assembly steps, §5 the code-vs-paper ledger — where the upstream `.wl` and
the paper disagree the code wins, §6 the typed-refusal contract, §7 the
gate skeleton, §8 the micro-fixtures F1/F2/F3). That specification is not
shipped as a separate document; its binding content is restated in this
file, in `CONTRACTS.md`, and executed by the tests, and the section
numbers are kept in the comments as stable cross-references.

Behavioral reference: `reference/SubTropica.wl` (read-only; line numbers
cited in comments). Shared types: `src/types.jl` and `src/b1_types.jl` —
every module codes against them as they are; a struct change is a
package-level decision, never a local edit inside one module.

====================================================================
## 1. THE TWO WARNING PINS — restated
====================================================================
**WARNING-1 — counterterm NET sign is (−1)^|subface|.**
`STSubtractionFormula` builds `Times[-a[subface],…]` with
`a[s] = (−1)^(|s|+1)` (SubTropica.wl:10823/10833), so the NET assembled
sign is −(−1)^(|s|+1) = **(−1)^|s|**: empty subface → +1 (identity term),
single facet → −1 (subtracted), facet pair → +1. The wrong reading
("(−1)^(|s|+1)") ADDS single-facet counterterms. Use
`warning1_sign(subface)` from `b1_types.jl`; the pole-cancellation gate
(§4, G3) must be able to FAIL when a sign is flipped (negative control).

**WARNING-2 — I vs J: the completion/division/gauge subset is the
COMPLEMENT.** **I** = the FIRST (lexicographic) index subset with
`det(rays[face][:,I]) ≠ 0` (wl:10777-10782; eq 3.34);
**J = Complement(1:n, I)**. The Kronecker-delta completion rows in
Vol = |det(rays|e_J)|/∏(−TropI), the output division /∏x_J, AND the
gauge-fix (vars ∉ J → 1, wl:10834) ALL use **J**. The wrong reading J = I
yields finite, numerically plausible, wrong-volume wrong-gauge results —
silent until the face-sum comparator. Micro-fixture F2 (I = [2], J = [1])
is the mandatory discriminator (gate G4).

====================================================================
## 2. COMPONENT LIST + FILE MAP
====================================================================
| Module | Files | Component / scope |
|---|---|---|
| polytope | `src/polytope.jl`, `test/test_b1_polytope.jl` | C4 — Newton-polytope layer via Oscar's first-class API (`newton_polytope`, `minkowski_sum`, `facets`, `vertices`, `affine_hull`; RT-iv). Emits `TropicalData`. 1-variable analytic branch; cache keying. Oscar in-process, NOT thread-safe — serial. |
| tropical | `src/tropical.jl`, `test/test_b1_tropical.jl` | C3+C5 — TropI max-plus evaluation, ray classification (`RayClass`/predicates from b1_types.jl), Σ_div via vertex-set argmax intersection (`SigmaFace`, ascending ordering law), w-search (`DivergentFacet`; typed refusals of kinds 1/3/4/5/6). Also owns the initial-form restriction `restrict_poly` (eq 3.30, wl:10506-10534) as an exported utility for subtract. |
| subtract | `src/subtract.jl`, `test/test_b1_subtract.jl` | C7 — Möbius counterterm assembly (steps 0-9). Emits `Vector{CounterTerm}` grouped by face in sigma_div order. WARNING-1/WARNING-2 live here; both micro-fixtures' negative controls are mandatory in its test file. |
| continue | `src/continue.jl`, `test/test_b2_continue.jl` | C6 — Nilsson–Passare continuation, the B2 route (§4c). |
| graph | `src/graph.jl`, `test/test_graph.jl` | Graph front end: propagator list → U,F Symanzik (matrix-tree, eU = ν−(L+1)(D−2ε)/2, eF = −(ν−L(D−2ε)/2), prefactor Γ(−eF)/∏Γ(ν_i), (−F)^eF sign tracked) → `EulerIntegrand`. Cross-checked vs the pySecDec `LoopIntegralFromGraph` glue. |
| fibrate | `scripts/fibrate/` (whole dir), `test/test_fibrate.jl` | Sound Goncharov fibration engine (Python/sympy: `fibrate.py`, `zip_num.py`, `gpl_num.py`): fibrates ZIP_reg words with ze-dependent letters into G({letters}; ze) × level-2 constants, with mandatory per-word numeric validation — replaces the engine's unsound `fibration_basis` op (see its README). best_order persistence+replay is `src/lr_dispatch.jl`. |
| lr_refine | `src/lr_refine.jl`, `test/test_lr_refine.jl` | Compatibility-graph LR refinement: a sharper provable upper bound on the letter set than the engine's simple Fubini reduction, certifying spurious blocking letters. The multi-group `find_lr_orders` dispatch (groups ALWAYS for sums), best_order persistence and verify_order-FIELD replay (R3) live in `src/lr_dispatch.jl`. No `find_lr_orders_scan` (Cheng-Wu gauge scan) wrapper ships — it is an engine op (CONTRACTS.md). |
| driver | `src/b1_driver.jl`, `test/test_tworeg_driver.jl` | Tier routing for census-certified divergent inputs (§4b item 3(iv), §4c). |
| shared | `src/types.jl`, `src/b1_types.jl`; `DESIGN.md`, `DESIGN_B1.md`, `CONTRACTS.md`; `reference/*`; `Project.toml`/`Manifest.toml`; `src/SubTropica.jl`, `test/runtests.jl` | — |

Test files label the modules by unit: unit-1 = hf_bridge/hlogexpr/serialize,
unit-2 = eps_expand, unit-3 = lr_dispatch/assemble, unit-4 = verify,
unit-P = polytope, unit-T = tropical, unit-S = subtract, unit-G = graph,
unit-F = fibrate, unit-L = lr_refine; "OWNER: unit-X" in a test header
names the unit whose code that file tests, and a unit's tests never depend
on another unit's code except through the shared types.

`include`/`export` wiring lives ONLY in `SubTropica.jl` (§4).

====================================================================
## 3. STANDALONE-TEST PATTERN
====================================================================
Every B1/B2/C test file is self-contained both ways: it can run on its
own under a plain Julia environment that provides Nemo and JSON (and
Oscar, for the polytope tests), and it can run inside the package suite.
To make that possible the file loads its dependencies directly and
includes the source files it needs by path, in dependency order:
```julia
using Test, Nemo, JSON            # polytope tests additionally: using Oscar
include(joinpath(@__DIR__, "..", "src", "types.jl"))
include(joinpath(@__DIR__, "..", "src", "b1_types.jl"))
include(joinpath(@__DIR__, "..", "src", "<module>.jl"))
```
Standalone invocation: `ulimit -v 32505856; julia --project=@. -e
'include("test/test_<module>.jl")'` (the address-space limit is the
convention used for every engine or helper process in this package).

`test/runtests.jl` auto-includes every `test/test_*.jl` and sandboxes
each file in its own module, so a standalone-first file joins the suite
with no edits and cannot collide with `using SubTropica` in `Main`. A
test that needs Oscar, the engine, or an external fixture directory and
does not find it SKIPS LOUDLY (`@warn` + `@test_skip`), never passes
silently.

====================================================================
## 4. INTEGRATION STRUCTURE
====================================================================

1. **Oscar → Project.toml** of the subtropica environment (Manifest
   resolved once; in-process Polymake verified under the 31 GiB
   address-space limit).
2. **includes → SubTropica.jl**: `b1_types.jl` after `types.jl`, then
   `polytope.jl`, `tropical.jl`, `subtract.jl`, `graph.jl`, `lr_refine.jl`;
   exports for the public surface (`tropical_data`, `classify_ray`,
   `sigma_div`, `produce_ws`, `subtraction`, `warning1_sign`, refusal
   types). On the B1 route C6 stays the loud typed refusal (RT-ix); the
   continuation itself ships as the B2 route, `src/continue.jl` (item 5).
3. **Driver** — the `STExpandIntegral`-equivalent (wl:11029-11145):
   quadruple → C4 → C5 (refusal ladder in the order ill-defined →
   degenerate → power → geometric property → w-norm → size guard) → C7 →
   C8 (full: per-face per-order counterterm expansion) → C10 → C12
   assembly. check_divergences: hard-ON standalone, OFF per counterterm
   inside face sums (CONTRACTS §a policy).
4. **B1 acceptance gates** (G1-G9; all must pass before any production
   use):
   - the **3 micro-fixtures** F1/F2/F3: exact FaceTerm lists, closed
     Beta/Γ totals, Laurent ledger;
   - **eq 4.3 pole tower** (eikonal, 1/ε⁴): all pole orders present with
     correct multiplicity;
   - **pole-cancellation crosscheck** (Gate B(iii) = G3): per fixture, the
     sum of counterterm pole parts cancels the raw divergence exactly
     (Γ-identity + per-Laurent-order), a local-finiteness probe along the
     divergent rays, AND the flipped-sign negative control that MUST fail
     by a strictly-greater-than-slack margin (a check that only prints and
     returns success is not a gate);
   - **a further divergent target** (Gate B(ii)) of the class where
     sector decomposition stalls, cross-checked ≥30 digits through ε⁰
     against an independent evaluator, the held-out-oracle rule B(v)
     respected;
   - **typed-refusal tests** (G6): eq 3.37 raw integrand ⇒
     `GeometricPropertyViolated`; a > 0 ray ⇒ `PowerDivergentRefusal`;
     ungauged homogeneous integrand ⇒ `B1Refusal(:TropIllDefined)`;
     segment-in-2d ⇒ `B1Refusal(:DegeneratePolytope)`; plus the
     `WNotUnitNormalized` D2 pin case.
5. **B2 = C6 continuation, a SEPARATE route**:
   Nilsson–Passare continuation + STfindFirstNPContinuation subset search,
   eq 3.37 as the u-convention unit fixture, eq 4.3/4.11 where they need
   continuation, relaxation of the `w·ρ = −1` pin (D2) by proof or by
   continuation. Nothing in B1 silently anticipates B2: power-divergent
   and geometric-property inputs REFUSE, typed, always.

====================================================================
## 4b. INTEGRATION NOTES
====================================================================
1. **Oscar → Project.toml**: Oscar 1.7 (Nemo 0.54 pin kept), Manifest
   resolved once. `import Oscar` with a selective `using Oscar: ...` sits
   in SubTropica.jl — a blanket `using Oscar` makes names such as
   `coefficients` ambiguous against Nemo's; only polytope.jl consumes the
   Oscar names.
2. **Includes/exports**: b1_types.jl right after types.jl; then (after the
   tier-A includes) polytope → tropical → subtract → graph → lr_refine →
   b1_driver. Full public surface exported (SubTropica.jl B1 export block).
3. **Integration items**:
   (i) types.jl `nu` measure LAW — the quadruple boundary is FLAT-measure
       ν; the B1 tropical layer reads DLOG ν = ν+1; the bridge is
       b1_driver `_b1_dlog`. Front ends (graph.jl) emit FLAT ν (as the
       fixtures do).
   (ii) subtract.jl `_c7_trop_poly/_c7_trop_on_ray/_c7_restrict_poly`
       DELEGATE to tropical.jl's exports; parity + hand-pin tests in
       test/test_b1_parity.jl; test_b1_subtract.jl's standalone header
       includes tropical.jl.
   (iii) **EpsExp Base-method extension check: NO COLLISION.** b1_types.jl
       is the SOLE definer of Base.zero/one/iszero/+/−/* for EpsExp; no
       other src file defines any EpsExp Base method. EpsExp*EpsExp stays
       deliberately undefined (an ε² bug detector).
   (iv) **C9/C11 slot: src/b1_driver.jl** (divergent driver).
       Census routing lives in subtropica_integrate: the tier-A path is
       unchanged for convergent/uncertifiable inputs; census-certified
       log-divergent inputs run C4→C5→C7→per-counterterm C8→per-FACE LR +
       hf_integrate on the face's per-order counterterm SUM (eq 3.33 ⇒
       order-by-order locally finite ⇒ **check_divergences stays HARD-ON
       everywhere** — the engine's detector doubles as a live negative
       control) → assemble_faces with cross-face pole alignment
       (align_laurent general fold).
       Correction found at the eq 4.3 acceptance gate: hard-ON is the
       armed FIRST LINE; on HFDivergent the driver retries that one order
       with the check OFF, stamped in facemeta["divergence_check_fallbacks"].
       The engine pre-check is per-log-atom rational and false-positives
       on transcendental infinity-tail cancellation across different log
       atoms (eq 4.3 face [2,4,10] k=1: Log[x6] vs Log[x6+1]; the
       counterterm sum decays as x6^-2·log). This is the CONTRACTS §a
       face-sum check-OFF situation; fallback values are verified by the
       gate's published-value + quadrature legs.
       CAVEAT: the fallback also converts a GENUINE non-cancelling face
       sum (e.g. a sign error) from a thrown HFDivergent into a stamped
       check-OFF value — at the engine-call level the detector is first
       line, not final. The loud negative control therefore lives in
       (a) the evidence stamp (any fallback order must be externally
       verified) and (b) the gate's value legs (published closed form,
       mpmath quadrature of the subtracted integrands, local-finiteness
       λ-scan probe) — the eq 4.3 gate demonstrates a flipped-sign
       counterterm shifting a face value by O(1) and being caught there.
       A sharper design would make the fallback conditional on an
       in-driver numeric decay probe of the summed integrand along the
       flagged ray.
       FARM DECISION: the face loop is SERIAL in-process, engine CLI
       subprocess transport (one call per face×order). The scale-out seam
       is `_b1_face_result` (a self-contained per-face body) plus
       scripts/fibrate-style job JSONs plus the patched C-ABI worker when
       a production target needs it. The `workers` kwarg is accepted and
       recorded; its semantics are unchanged in this release.
4. **runtests.jl auto-include** sandboxes each test file in its own
   module (standalone-first files include() src by path and cannot share
   `Main` with `using SubTropica`). test_b1_polytope.jl loads Oscar into
   `@__MODULE__`.

====================================================================
## 4c. B2 INTEGRATION NOTES
====================================================================
1. **Wiring**: src/continue.jl is included in SubTropica.jl after
   subtract.jl and before b1_driver.jl (it consumes the tropical.jl +
   subtract.jl surfaces; b1_driver's B2 route consumes it). Full B2 surface
   exported (continue_ray/continue_rays/find_first_np_continuation/
   expand_integral_b2/Eps2/Euler2Integrand/promote_second_regulator +
   B2Refusal/B2Prefactor carriers).
2. **Driver routing**: `subtropica_integrate` has `allow_continuation`
   (default FALSE — the typed refusals REMAIN the default surface) and
   `continuation_order` kwargs. When true, PowerDivergentRefusal /
   GeometricPropertyViolated from the census route into
   `subtropica_integrate_b2` (b1_driver.jl): geometry rebuilt from
   `_b1_census_polys`, continuation in the DLOG chart (`_b1_dlog` bridge),
   `expand_integral_b2`, then each continued integrand bridges back to the
   FLAT boundary (`_b1_undlog`) and re-enters `subtropica_integrate` with
   allow_continuation=false (ONE pass — the wl NP block runs once and its
   step ledger covers GP-removal AND power-shaving; a continued output
   that still refuses propagates typed). The exact `B2Prefactor` per term
   is convolved via `_b2pref_rational_series` (exact rational Laurent; the
   ε-proportional denominator factors are the continuation poles), terms
   are letter-unioned + pole-aligned (align_laurent) and summed. Evidence
   carries the full continuation ledger under "b2_continuation" (trigger
   refusal, np_search, per-term steps/prefactor/pole_shift/sub-route/
   sub-evidence). WNotUnitNormalized does NOT route (only the two
   continuation-shaped kinds); inside find_first_np_continuation it still
   counts as subset failure.
3. **`order` kwarg semantics**: `expand_integral_b2(; order=16)` is the
   TOTAL continuation STEP cap (guard 6 in its docstring) — NOT a Laurent
   order. The per-ray step multiplicities are always DERIVED from the
   census (wl:11110-11122), never guessed; the cap only bounds the
   exponential blow-up (paper Sec 3.2.2, bottleneck 1) and refuses typed
   (:ContinuationDepthExceeded). To avoid collision with the driver's
   Laurent `order` kwarg, subtropica_integrate surfaces it as
   `continuation_order`.
4. **find_first_np_continuation signature note**: the second-argument
   NamedTuple overload (`dd` from `divergence_data`) is only REACHABLE for
   GP-clean inputs — divergence_data throws GP before returning — so the
   driver route always uses the Vector{Int} form (`div_facets(trvals)`
   inside expand_integral_b2). The overload is kept for API symmetry (its
   docstring says so); the div_facets input must ascend (the Σ_div
   ordering law, checked).
5. **Message law**: b1_types `_NP_SUFFIX` names the route ("B2: rerun
   subtropica_integrate with allow_continuation=true"); message-law tests
   in test_b1_{subtract,parity,tropical}.jl assert the wording.
6. **B2 acceptance gate**: the mini power fixture runs end-to-end through
   subtropica_integrate(allow_continuation=true) — continued → subtracted
   → integrated → assembled Laurent — against the sympy-verified closed
   form Γ(ε−1)/Γ(1+ε) = −1/ε −1 −ε −ε² (test_b2_continue.jl,
   test_tworeg_driver.jl).
