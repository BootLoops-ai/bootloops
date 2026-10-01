# Porting notes: SOFIA (Correia, Giroux and Mizera; written in the Wolfram Language) → SOFIA.jl

The SOFIA authors' sources (scr.m, scrj.m) are pinned under `reference/` (commit in
`reference/UPSTREAM_COMMIT.txt`). Line numbers below refer to
`reference/scr.m`. This file records every deliberate behavioral divergence
and the architectural translations, so that validation discrepancies can be
traced to decisions rather than accidents.

## Architectural translations (same semantics, different mechanics)

1. **Momenta as rational vectors** (`src/baikov.jl`). scr.m manipulates
   symbolic momenta through `CenterDot` pattern rules and compresses repeated
   subexpressions into `Λ[a][i]` symbols (scr.m:714-736) — devices for
   Mathematica's pattern matcher. We represent momenta as explicit vectors
   over `[l_1..l_L, p_1..p_{n-1}]` and expand Gram matrices bilinearly once.
2. **Loop-momentum shift** (SOFIA's `FindLoopMomShift`, scr.m:649): replaced
   by a direct bounded search over candidate shifts (negated rest-vectors),
   optimizing the same objective (number of distinct loop dot products).
3. **`FixLoopEdges`** (scr.m:518-596): the GF(2) cycle-matrix heuristic is
   replaced by spanning-forest co-trees, which always give a valid loop-edge
   set; `lbl` retries alternates and validates against the variable bound.
4. **`rankDropRule`** (scr.m:626): same algorithm; random prime evaluation
   points replaced by deterministic-seeded rational points, checked twice.
5. **Resultant compression** (SOFIA's `generateAnsatzDISPATCH`, scr.m:420):
   same trick — generic resultants cached by degree pair over
   `ZZ[a_i, b_j]`, evaluated at actual coefficient polynomials. High-degree
   pairs fall back to subresultant PRS.
6. **Kinematics** (`GenerateKinematics`, scr.m:172): same cyclic Mandelstam
   basis; the linear solve is exact Gaussian elimination over QQ. Squared
   masses are ring generators from the start (scr.m homogenizes late via
   `homogeneizeKin`); node and edge labels share one generator namespace, so
   a node mass equal by label to an internal mass is the same scale — as in
   scr.m after homogenization.

## Deliberate divergences (different, and why)

7. **Exact proportionality** (`isproportional`): SOFIA's
   `ProportionalPolynomialsQ` tests at random integer points; ours is exact.
   Strictly fewer false identifications.
8. **Discriminants** are computed as `res(f, ∂f)` without the leading
   coefficient division; the lc is added to the system separately anyway
   (Fubini's `term1`) and everything is factored downstream — the
   singularity set is unchanged.
9. **Effortless** (`src/effortless.jl`): perfect-square tests are exact
   (`Nemo.is_square`) instead of 4 random integer points; the degree<5
   zero-locus pre-filter (`AllowedEvenLetters`) is skipped — it only prunes
   candidates the exact test rejects; independence of odd letters is decided
   by exact linear algebra on dlog forms at deterministic rational points
   instead of `NSolve` at 50 random real configurations. The FiniteFlow
   variants are not translated (the plain fallbacks in Effortless.m by Antonela
   Matijašić and Julian Miczajka are what is translated).
10. **Subtopology dedup** is by canonical edge order; SOFIA's
    `DeleteDuplicates` is order-sensitive (weaker). The symmetry quotient
    would merge such duplicates anyway.
11. **The GUI** (`FeynmanDraw`) and plotting are not translated.

## Route dependence and the `routes` option

The candidate-singularity set on the maximal cut is **route-dependent** in
both implementations: it depends on the loop-edge choice and order (which
momentum is `l_1` vs `l_2`), the pivot order of the dot-product → Baikov
change of variables, and the equation order (in scr.m: whatever
`FixLoopEdges` and Mathematica's `Solve` pick). The SOFIA authors' example notebook itself
demonstrates this: the "Degenerate acnode" example shows automatic
`LoopEdges` giving an incomplete ("BAD") set that pinning `LoopEdges->{1,5}`
fixes.

SOFIA.jl exposes this honestly: `routes=N` unions the candidates over up to
N deterministic route variants (loop-edge choices × both orders × pivot
orders × equation orders), producing a robust superset of any single route —
including, empirically, the one scr.m takes. With `routes=1` you get one fixed route,
comparable to (but not bit-identical with) SOFIA's default.

## Symmetries

`symmetry_map`/`symmetry_quotient` (`src/symmetry.jl`, on by default with
`symmetries=true`) implement `SymmetryQuotient` (scr.m:1107-1226) as a
performance-only feature: subtopology classes are computed once and results
transported through the kinematic maps (e.g. double box: 55 subtopologies ->
18 classes, identical output). Analysis of `leadingTerm` (scr.m:1260) shows
the `Symmetries->True` pipeline applies the quotient maps and keeps the full
δ⁰-mapped polynomials — its output equals the direct per-subtopology
computation SOFIA.jl performs, so the quotient can only affect speed, never
results.

## The PLD.jl backend bridge (`src/pld.jl`)

SOFIA's `Solver->momentumPLD` option lives on the Julia side of scrj.m
(lines 1234-1317): `prepareVariables` classifies the Gram-determinant
variables into active (integration) variables and kinematics,
`defineCoefRing`/`defineLaurentRing` build R = QQ[kinematics] and a
Laurent ring S over R in the active variables plus one α per Gram
determinant, `sendSystemToJulia` assembles Δ = Σ αᵢ Gᵢ with the LAST α
set to 1, and `ComputeDiscriminants` calls PLD.jl's `getSpecializedPAD`
with the six `PLD*` options. The port maps these one-to-one onto
`pld_delta` (Δ assembly, backend-independent, Nemo side), `pld_system`
(the diagram front end via `prepare_landau_system`, Gram list sorted by
term count as scrj.m does), and `pld_discriminants` /
`pld_singularities` (the backend call); keyword defaults are SOFIA's
option defaults (`PLDMethod->sym`, `PLDHomogeneous->true`,
`PLDHighPrecision->false`, `PLDCodimStart->-1`, `PLDFaceStart->1`,
`PLDRunASingleFace->false`, scr.m:1231). Faithfully to scrj.m's
`allActive`, ALL α's become Laurent-ring variables even though the last
one no longer appears in Δ.

Deliberate divergences: scrj.m launches a separate Julia session,
`Pkg.add`s the Oscar stack, activates a hard-wired local PLD checkout at
package load (scrj.m:1236-1250), and returns only the backend's printed
output; here the caller passes the loaded `PLD` module explicitly (so
PLD.jl never becomes a dependency and the bridge costs nothing when
unused) and the specialized discriminants come back as Nemo polynomials
in QQ[kinematics], deduplicated and (with `factor_result=true`, SOFIA's
`FactorResult->True`) factored — directly comparable with the FastFubini
lane's output. Polynomial transport between the Nemo and Oscar worlds is
term-by-term by variable name (`transplant`), replacing scrj.m's string
round trip through `ExternalEvaluate`; it is duck-typed against the
AbstractAlgebra generic-ring interface, which is how `test/pld_tests.jl`
can drive the full marshalling round trip against a mock backend when
PLD.jl is not installed (the live leg runs when it is). The
`EulerDiscriminantQ` cross-check (scrj.m:1328) is not wrapped.

## Validation

`validation/` replays the worked examples recorded in the SOFIA authors'
`SOFIA_examples.nb` (diagram + options + recorded output, extracted from the
notebook's box forms) and compares singularity sets up to proportionality
and variable naming: see `validation/REPORT.md`. The 96-file
`reference/PLD_database` doubles as an end-to-end oracle (`dbox_zero_zero`
is reproduced exactly).

## Performance frontier (heavy multi-scale topologies)

Six notebook validation cases (multi-scale acnode/outer-dbox/penta-box/
double-pentagon class) exceed practical time on this engine while the SOFIA
authors' Mathematica code completes in seconds-to-minutes. Minimal reproducer: the outer double-box's
subtopology rep 3 reaches, on our default route, a discriminant of a
6175-term degree-8 (total degree 15) polynomial; both FLINT's native
multivariate resultant and classical subresultant PRS stall on it (>20 min,
uninterruptible mid-call). The SOFIA authors' route — fixed by Mathematica `Solve`
pivoting and `FixLoopEdges` internals we deliberately did not replicate
bit-for-bit — never encounters this object.

Mitigations in place: route deadlines (`time_budget`, RouteTimeout),
UnclogTime-style subtopology skipping, and the `maxopterms` op-size guard
(SolverBound-like skip semantics). These make every case terminate but trade
away deep content on high-loop cases.

Future directions: (a) monster detection with route reordering (try other
routes when predicted op cost explodes — the union machinery already
exists); (b) per-subtopology worker processes with hard kills, so a stalled
route cannot block the union; (c) modular/evaluation-interpolation
discriminants specialized to few kinematic variables.
