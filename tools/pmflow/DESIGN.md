# gravityFlow — design document

**Scope:** self-consistent auxiliary-mass-flow for cut-eikonal (post-Minkowskian) integral
families, where the standard AMFlow recursion is a fixed point.

---

## 1. The problem

Cut-eikonal PM families carry linearized matter-line propagators `1/(2k·u_i)` that are
forbidden η-positions (matter lines are on-shell cuts).  η placed on the remaining
mass-prop slots sends η→∞ and reproduces the same propagator structure exactly: the
boundary recursion `b0→b1→b2→…` descends to byte-identical families at every level.
How the fixed point shows itself (the diagnostic evidence chain, reproducible with
`AMFLOW_DEBUG_SCHEME=1`):

- the scheme trace fires Tradition on the same mass props at every boundary node;
- `integralfamilies.yaml` diffs across levels (e.g. b3↔b4) are byte-identical
  propagator lists;
- the Cutkosky EndingScheme declines: `−2k·u` cut components are not a phase volume;
- all three scheme variants (Tradition-only, MMA order, Cutkosky-first) fire Tradition.

**Consequence:** vanilla AMFlow cannot produce boundary data for these families.  The
standard PM literature handles this via y-flow from static (y→1) boundaries.  gravityFlow
implements the self-consistent alternative: treat the fixed point AS the boundary condition
by solving the linear-response eigenproblem, then injecting via `explicit_boundary`.

The reference family used throughout this document is a 4PM cut-eikonal family with
3 matter-line cut propagators, 15 propagator slots, and 43 masters; the anchor sector
is a factorized cut tadpole with the closed form `J = −(i/8)π^{−3/2}Γ(ε−½)³` (called
J63 below, after its sector index; a second sector J111 satisfies the parameter-free
identity J111 = J63).

---

## 2. Pipeline

```
detect            probe (parallel)        solve                inject
  |                    |                    |                    |
FIXEDPOINT   →   ETAC_OVERRIDE     →   linear-response   →   FORCE_ENDING
  PROBE            (N flows α,β,…)      assembly + anchor      + explicit_boundary
                                        (mpmath, high dps)     → parent AMFlow run
```

The CLI (`pmflow.py`) drives all four stages; the engine-side hooks are four small
patches in amflow-cpp's `src/pipeline/amfsystem.cpp` (the amflow-cpp fork lives in
the sibling repository `amflow-cpp`): a fixed-point probe
mode (`AMFLOW_FIXEDPOINT_PROBE`), a forced
ending depth (`AMFLOW_FORCE_ENDING_DEPTH`), an η-position override
(`AMFLOW_ETAC_OVERRIDE`), and the self-similarity guard (§4b) that aborts a plain
run on a fixed-point family with the structured `GRAVITYFLOW_FIXEDPOINT` code
instead of recursing to the depth cap.

### Stage 1 — detect

Run the b0 family with `AMFLOW_FIXEDPOINT_PROBE=1`.  amflow-cpp skips recursive
sub-system setup and runs unit-vector transports (masters × ε-grid points; 43×12 in
~4 min on the reference family), writing the response matrix `fp_A_<tag>.json`.

The single-flow solve then:
1. drops phantom preheat-only columns (zero-norm check);
2. computes `ker(A(ε)−I)` via SVD at each ε;
3. tests the parameter-free ratio `J111/J63 == 1` and anchors via the analytic J63
   closed form.

On the reference family: `dim ker = 28`; 15 coupled rows, 28 free (pure-(1+η)-power
cut-vacuum eikonal masters).  The fixed point alone is rank-15 / 43 — one flow does
not determine the boundary.  Multi-flow is mandatory.

### Stage 2 — probe (parallel)

Run N independent flow variants by setting `AMFLOW_ETAC_OVERRIDE` to different
non-cut-prop η-subsets.  Each variant α yields a matrix A_α(ε) and a constraint set
`(A_α−I)M = 0`.  The joint rank of the stacked constraint system must reach the master
count (one analytic anchor saturates the remaining free direction).

`pmflow respond` generates one unit-injection JSON per key and runs them N-way
parallel, one subprocess per injection unit, with `AMFLOW_DUMP_EPS_GRID=1` dumping
per-eps master values to the probe logs.  `pmflow discover` auto-discovers the
required injection key set by iteratively probing with zero-injections and collecting
`no Vacuum entry` MISS keys.

The probe harness is multi-process (one child AMFlow invocation per probe; ~10 s of
startup overhead each).  An in-process batched mode is a possible extension (§4a),
not implemented.

### Stage 3 — solve

The multi-flow assembly (`pmflow solve`, delegating to `gf_solve.py`):
1. loads all probe-log EPS_GRID outputs (certified-ball parser);
2. builds the linear-response matrix C (n × n_keys) column-by-column;
3. builds `(I−CP)` for the etaC constraints and stacks the `(A_α−I)` fixed-point
   constraints into an overdetermined system (86×43 on the reference family);
4. equilibrates rows AND columns (the raw stack spans ~300 orders of magnitude);
5. pivoted row-selection forward-elimination → square system → `mp.lu_solve`;
6. validates: joint residual, the anchored masters against the analytic J63 value,
   and the J111 = J63 identity;
7. writes the per-ε explicit-boundary table keyed `<family>|i1|…|i15`.

The single-flow eigensolve path (no etaC ring) and the multi-flow + etaC assembly
share the J63 anchor formula.

### Stage 4 — inject

Set `AMFLOW_FORCE_ENDING_DEPTH=1` so that any boundary family at depth ≥ 1 takes the
Trivial ending scheme and reads values from `amf_options.explicit_boundary`
(populated from the solve-stage table).  The parent run proceeds through all 4 reduce
stages (ibp-cache HITs) + the η-DE solve with the now-anchored boundary.

---

## 3. Where the pieces live

| component | where |
|---|---|
| fp-probe hook (`AMFLOW_FIXEDPOINT_PROBE`) | amflow-cpp `src/pipeline/amfsystem.cpp` |
| force-ending hook (`AMFLOW_FORCE_ENDING_DEPTH`) | amflow-cpp `src/pipeline/amfsystem.cpp` |
| etac-override hook (`AMFLOW_ETAC_OVERRIDE`) | amflow-cpp `src/pipeline/amfsystem.cpp` |
| self-similarity guard (`GRAVITYFLOW_FIXEDPOINT`, §4b) | amflow-cpp `src/pipeline/amfsystem.cpp` |
| detect / discover / respond / inject / map drivers | `pmflow.py` (this dir) |
| closure solver (assembly, anchor, validation) | `gf_solve.py` (this dir) |
| ε⁰ ratio salvage + rational identification | `gf_eps0_ratios.py` (this dir) |

---

## 4. Possible extensions and the self-similarity guard

### 4a. In-process batched injections (would replace the subprocess probe farm) — not implemented

The per-probe subprocess model pays AMFlow startup and ibp-cache cold-MISS overhead on
every first-time key.  Replacement: extend the fp-probe code path inside
`amfsystem.cpp` to run all N unit-vector injections in a single process via a column
loop (the current code runs one column at a time), amortizing the 4 reduce stages
across all masters × ε columns.  Estimated speedup: 5–10× for a warm cache, far more
for a cold one (no per-probe startup).  Target: a `mode: fixedpoint_batch` JSON field
routing through the existing `AMFLOW_DUMP_EPS_GRID` machinery without spawning child
processes.

### 4b. Automatic self-similarity guard — implemented

The engine detects the fixed-point condition automatically.  `build_boundary`
records the parent family's canonical propagator signature (sorted
propagators-after-conservation) in the sub-system options; the child's `setup()`
— the recursion continuation point, which endings never reach — compares its own
signature against it.  On a match at any depth it raises a structured error
carrying the depth and the `GRAVITYFLOW_FIXEDPOINT` diagnostic code, rather than
silently recursing to `max_recursion_depth`.  This makes gravityFlow the mandatory
entry point instead of a manual workaround: `pmflow detect` parses the code as the
engine-certified FIXED_POINT verdict.

Scoping: the guard stands down when `AMFLOW_FORCE_ENDING_DEPTH` or
`AMFLOW_FIXEDPOINT_PROBE` is set — those are exactly the pipeline stages that must
traverse self-similar levels (forced endings terminate the descent; the probe
bypasses it) — and `AMFLOW_ALLOW_FIXEDPOINT_RECURSION=1` restores the silent
recursion for debugging.  Factorized (SingleMass) sub-families drop the inherited
signature: the guard tracks only the η→∞ boundary chain.

### 4c. Symbolic cut-tadpole anchor calculus — not implemented

The J63 closed form `−(i/8)π^{−3/2}Γ(ε−½)³` is the only built-in anchor value: the
solver's `--anchors "SECTOR:FORM,..."` flag places it on any corner master of the
fpA basis (sector index decoded per family), but `J63` is the only FORM it accepts.
A sympy module computing the analytic value of any factorized cut-tadpole sector
`(cut δ(2k·u))^L` via dimensional regularization would validate the anchor
identities symbolically (not just numerically), generalize to 3PM/5PM families with
different cut multiplicities, feed `e^{γε}` convention flags automatically, and
populate the FORM namespace beyond J63.

---

## 5. Two-flow contour-convention consistency (settled subtlety)

Two independent η-flow probes on the same b0 basis can return O(1)-disagreeing
boundary values before anchoring.  In principle this could be physical: if two flows
deform different branches of the cut measure `d^dk δ(2k·u)` (different η→∞ contour
conventions), the multi-flow solution would be convention-dependent and the injected
boundary wrong in sign/phase for one flow direction.

The decisive check: after anchoring, census EVERY redundant row.  The solve selects a
square subsystem from the overdetermined stack; the non-pivot rows are the consistency
check.  If their residuals scale as `ε⁰` (not `ε^{−k}`) and are at noise level after
anchoring, the flows deform the same contour and are truly redundant.

Measured resolution on the reference family (full closure at one kinematic point:
etaC linear response, rank 41 after two Γ³ anchors, remaining 2-dim null space solved
by least squares against the second flow's rows): **43/43 A_α rows consistent,
residual quartiles 0 / 1e-169 / 1e-131 / 3e-53**.  An apparent O(1) "disagreement"
arises when comparing against an INCOMPLETE solution (top-sector free directions set
to zero) — it is an artifact, not a contour effect.  The solve step therefore always:
anchor → eliminate → solve the null space against the second flow → census ALL rows
(the residual census IS the validation).
