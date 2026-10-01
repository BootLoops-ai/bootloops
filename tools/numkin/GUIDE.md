# numkin — guide

Tool page: https://bootloops.ai/tools/numkin.html

KIND: package (bash stage/run driver + python reconstruction + python shift optimizer;
needs a working kira + FireFly + Fermat stack to be useful — see the Kira and
AMFlow.cpp forks from the sibling repositories kira and
amflow-cpp, cloned beside this one, or obtain upstream)

PURPOSE: freeze-one-scale rescue for stuck multi-var symbolic FireFly reductions:
sed-substitute one kinematic invariant to numeric values in already-generated SYSTEM
files, run N cheap 1-var (d-symbolic) FireFly solves, reconstruct the frozen variable
exactly by Padé/Thiele with held-out verification. Plus `shift_opt.py`: per-sector
unimodular loop-shift probing that shrinks η-DE reductions. Members: the eta-rerun verb
(`etarerun/`) and basisland (`basisland.py`, the basis-solve backend).

USE-WHEN:
- Kira multi-var symbolic FireFly (2+ invariants) is "hours out" / stalls on a subset
  of high-degree entries while single-var d-only FireFly on the same system is fast.
- Need missing symbolic-DE entries when full multi-var FireFly is hours-days out
  (sed-sub route: ~10–40 min/point + exact bivariate Padé).
- An η-DE reduction is hours-long on a subsector while a sibling was fast under a loop
  shift → `shift_opt` FIRST (direct 2-var Fermat on the winning shift may bypass the
  numeric-d farm entirely).
- amflow η-DE reduction stuck/decelerating in Fermat back-sub and the SYSTEM alphabet
  is exactly {Eq,d,eta} → `etarerun/etanumd_convert.sh` (numeric-d FireFly bypass).
- Target = element of a finite-dim function space with an exactly-evaluable spanning
  set → `basisland.py` (see its section below).

NOT-FOR: many probes at ONE kinematic point; regenerating systems (the whole point is
reusing the sunk gen cost); `shift_opt` is zero-ISP targets only (zero-ISP values are
shift-invariant — ISP-carrying targets are not covered).

INVOKE:
- `numkin_sweep.sh stage <FAM> <SYS_DIR> <VAR> <VALUE> <OUT_DIR> <PT_ID>
  <PREFERRED_MASTERS> <TARGETS_FILE> <REDUCE_YAML_FRAGMENT>` — VAR e.g. m2 (NOT d; d
  stays symbolic), VALUE integer or "p/q"; SYS_DIR = completed `run_initiate:true` dir
  with `tmp/<FAM>/SYSTEM_<FAM>_*.gz` + config/ + sectormappings/.
- `numkin_sweep.sh run <OUT_DIR> <PT_ID> <NTHREADS> [KIRA_BIN] [FERMATPATH]`.
- `python3 tools/numkin/numkin_harvest.py --family FAM --targets-file F --pts-dir DIR
  --var m2 --held-out 2 --out OUT.json [--out-m OUT.m] [--min-fit-pts N]`.
- `python3 tools/numkin/shift_opt.py {build|probe|harvest} ... [--catalog]` — catalog =
  standard shifts {kᵢ→l+P−kᵢ, l→l−P, l→l+kᵢ−P, composites}; auto-completes ISPs to full
  bilinear rank; probes under Fermat env (DCAP=1, RUN_FIREFLY=0, PREHEAT_DOT=0).
- etarerun verbs: `etanumd_convert.sh SRC DST P/Q` (d-sub + relaunch);
  `etanumd_2varff.sh` (no substitution, 2-var FF); `check_degd.py` (deg_d fit from ≥2
  slices; ETANUMD_DIR); `inject_ibpcache.sh` (AMFLOW_IBP_CACHE_DIR required);
  inject/ (residual-only slice splicing). The FARM_* launcher scripts are
  TEMPLATES: their input dirs are env-required (ETANUMD_* vars); no defaults are shipped.
- Battery (no kira stack needed): `python3 tools/numkin/selftest.py` — basisland
  rank/solve on a synthetic-truth system + the rank-deficit refusal.

INPUTS: a COMPLETED kira gen pass (SYSTEM_*.gz reused verbatim, `run_initiate:false`);
sweep grid of rational points — exclude degenerate slices (a degenerate slice =
deterministic wrong-but-valid reduction, byte-identical on rerun).

OUTPUTS: per-point staged dirs `pt_<ID>/`; harvest emits a documented JSON schema AND a
synthetic kira2math-format `.m`. shift_opt: shifted family definitions + per-sector
probe verdicts (each sector may want a DIFFERENT shift).

ENV: `NUMKIN_STAGE_PAR` — parallel sed staging via xargs -P (15–20 min → ~2 min/point);
`KIRA_BIN`, `FERMATPATH` pass-through on `run`. The etarerun launch verbs preload
jemalloc from `$JEMALLOC` (fallback: the Debian/Ubuntu libjemalloc.so.2 path) and
find Fermat via `$FERMATPATH` (fallback: `fer64` on PATH, then /usr/share/Ferl7/fer64).

GATES: harvest requires ≥2 held-out points per entry, verified exactly; Padé rank-sweep
certifies point-count sufficiency; `DGridInsufficient` CONTRACT — raises (never
silently pads) on nullity-0 fits, denominator vanishing at a d-node, or degD
non-uniform across the d-grid: the fix is more points / a denser d-grid, not
tolerance. shift_opt gated vs independent oracle (168.6d). Measured shift wins that
justify the probe: 5.4×/N_eqns on one sector, >200× on a Fermat wall.

FOOTGUNS:
- Staging behavior is load-bearing: substitution is WORD-BOUNDARY
  (`sed s|\bVAR\b|(VALUE)|g`, value parenthesized) AND stage() copies config/ then
  strips the swept var from BOTH `kinematics.yaml` (kinematic_invariants) and
  `integralfamilies.yaml` (propagator masses) — kira builds FireFly's variable registry
  from ALL config symbols; leaving the var anywhere keeps FireFly in slower multi-var
  mode on a provably 1-var system.
- NEVER copy kira.db between points (stale-export footgun); each point is a fresh
  staged dir.
- Mass-slice variant (multi-mass symbolic): FRESH-GEN per slice — staged-SYSTEM reuse
  segfaults kira write-out; `preferred_masters` + `integral_ordering: 8` fixes FireFly
  "Validation failed: Entry 0" + export crash; exclude mass-space Landau-degenerate
  pairs (λ=0 → zero-poly divmod).
- Known-Q window is n−8, not the blind ceiling (n−6)/2 — budget points BEFORE farming.
- shift_opt: some sectors still hit the asymmetric-η Kira cycle bug under a shift
  (`iteration>500, items left<30` → kill by PID).
- eta-rerun: eps~1e-4 slices are ~1.3× slower than eps=1/prime — pick small-numerator
  d. FF still prints "Reconstructing in: d, eta" after the d-sub but discovers
  deg_d=0 in scan — harmless. inject without rebase_invert.py loses 73% of shared rows
  when source masters ∉ target masters.

## basisland member (the basis-solve backend)

PURPOSE: basis-constrained exact function landing — when an unknown function lives in a
KNOWN finite-dimensional space, land it exactly from #exact-samples ~ dimension,
INDEPENDENT of #variables (measured: a 14-variable dim-543 function pair landed from
596 exact rows in 1.2 min wall, 596/596 exact-verified). Consumes MIXED constraint
data: exact point values AND exact linear functionals as rows of ONE linear system.

INVOKE: `python3 tools/numkin/basisland.py rank IN.json OUT.json [nprimes=3]` (cheap
dimension discovery + pivot subset, multi-prime agreement asserted, BEFORE expensive
exact data collection); `python3 tools/numkin/basisland.py solve IN.json OUT.json
[nprimes=96]` (multi-RHS; exact x-vectors + exact_verified counts + measured heights).

GATES (all fail-closed): 25-bit primes, overflow-safe int64 elimination; rank deficit →
exit 2 with missing columns; ratrec overflow → exit 3 (more primes; 96 ≈ heights to
1e180); surplus-row inconsistency at ANY prime → assert (the overdetermination IS the
certificate; never trim surplus rows); mandatory exact Fraction verification of EVERY
row → exit 4 on any miss. Held-out battery is the CALLER's leg. NOT a replacement for
Thiele/Padé in the no-basis case: no basis, no basisland.

FOOTGUNS: the pivot basis matters for HEIGHT, not correctness (same function: ~1e300
heights in one basis, 2.29e8 in a well-chosen pivot basis — suspect the basis before the
function); functional rows fitted numerically need their own positive control at
production parameters; exact rationals only — never float rows.

RELATED: tools/kira-stack (gen pass provider; η route selector); tools/ratfit
(reconstruction backend + grid design). Flat `tools/numkin_sweep.sh`,
`tools/numkin_harvest.py`, `tools/shift_opt.py` are compat shims/wrappers.

CREDIT: numkin drives Kira (Maierhöfer, Usovitsch, Uwer; Klappert, Lange; Wu [Kira1,
Kira2, Kira3]) with FireFly (Klappert, Klein, Lange [FF1, FF2]) and Fermat (R. H.
Lewis), and AMFlow's η-differential-equation systems (Liu, Ma [AMFlow]). Freezing one
scale and reconstructing it by Thiele/Padé interpolation with held-out verification is
the classical rational-interpolation strategy that FireFly [FF1] and the finite-field
programme of von Manteuffel–Schabinger [vMS] and Peraro [Per16] made standard; the
sed-substitution driver, shift optimizer and basisland are ours. Bracketed keys resolve
in REFERENCES.md at the repository root.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).
