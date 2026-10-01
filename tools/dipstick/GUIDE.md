# dipstick — guide

Tool page: https://bootloops.ai/tools/dipstick.html

KIND: package (julia + python mix; thin `dipstick` dispatcher execs the members:
pf_rank.jl + pf_rank_validate.jl + regress_order.jl, critcount.py, cycletype.py,
nonres.py, resaudit2.py + resaudit_fast.py (the `audit` verb), flatchi.py +
flatchi_aff.py (the `flatchi` verb), regions/regions.jl (the `regions` verb, the
expansion-by-regions completeness certifier — member section below, with its
inputs under regions/tests/ and its control battery regions/selftest.jl);
`run_msolve_capped.sh`, the mandatory msolve cap wrapper (member section below);
`battery.py` + `fixtures/` vendored battery fixtures; GKZ36.json +
w5_seed_base.json vendored data for the nonres verb; fixtures/audit_toy.ms toy
fixture for the audit/flatchi legs).

PURPOSE: The pre-compute triage tool, unified — check the integral BEFORE you commit
the farm. Four verbs, one regime split: `order` = Picard–Fuchs order + GKZ holonomic-
rank probe mod p from the (Baikov/LS) integrand polynomial (Griffiths–Dwork at random
s₀∈GF(p), 2-prime cross-checked, n≤4); `count` = holonomic rank of a many-factor
Euler integral ∏pᵢ^{sᵢ} WITHOUT symbolic expansion (Lagrange/ML-degree critical
system → msolve solution count; n≥5 — closes order's NOT-FOR wall); `cycletype` =
Frobenius cycle-type statistics of the critical points (Galois/monodromy evidence);
`nonres` = GKZ facet + β-pairing nonresonance/irreducibility certificate;
`audit` = msolve-INDEPENDENT exact resultant/Bezout audit of a 2-roamer-variable
.ms critical system (deg Res_y(P1,P2) minus arrangement-point valuations = crit
count; `--fast` = evaluate-interpolate F_p variant for n≥20 lines); `flatchi` =
exact flat-lattice Möbius χ for hyperplane arrangements in any ambient dimension
(characteristic polynomial, χ_top of the projective complement, Zaslavsky
bounded-region count; a .ms file argument takes the affine route); `regions` =
expansion-by-regions completeness certificate — is this region list ALL the
regions? — checked BEFORE any region-by-region series is built or resummed
(member section below).

WHY `audit` IS MANDATORY for fiber counts: msolve (0.6.5) does not range-check
characteristic-p coefficients in .ms input; early versions of our emitters wrote
coefficients outside [0,p) and obtained wrong counts with rc=0 (see FOOTGUNS
below). Every emitter now reduces mod p, and `audit` remains the
msolve-independent cross-check; 2-variable `count` results auto-chain into
`audit`.

USE-WHEN:
- Route triage BEFORE compute: order 2 → elliptic/Eichler, 3 → K3, 4 → CY₃/LCLM;
  feeds the A/B/C value-fittability decision.
- Sizing a DE-transport or ε-factorization run (order + rank bound the system before
  IBP reduction runs).
- Series-in-jet cases with deg_s(G)≥2 — the exact jet recursion
  P_{k+1}=∂ₛ(P_k)·G−(k+1)·P_k·∂ₛG handles them (`pf_probe` auto-routes on an exact ℚ
  nonlinearity check; linear-in-s G takes the fast path).
- Many-factor Euler/GKZ integrands where symbolic expansion dies: `count` builds the
  sparse Lagrange system (degree ≤ 1+deg p) and msolve counts — X(3,6)→26 in seconds
  where the product form (deg 31) killed sympy at X(3,7).
- Galois/monodromy evidence on a 0-dimensional critical system: `cycletype` at ~6
  primes (transitivity via common subset-sums, parity, Jordan criterion).
- Irreducibility of the GKZ module at a concrete β: `nonres` (Gelfand–KZ criterion;
  facets of cone(A) + exact pairings).
- Before resumming ANY expansion-by-regions / moment / threshold / large-mass
  series, and when auditing a hand-built or literature region list: `regions`.
  A region-incomplete series still converges — to a wrong value, with no
  diagnostic from the resummation itself.

NOT-FOR:
- `order`: n=5+ variables (comfort ceiling n=4) — THAT regime is `count`'s job: same
  function (holonomic rank), complementary regimes. The order-vs-count split is
  deliberate, not a workaround.
- `count` is a RANK probe, NOT a solver; "random integer exponents = generic-weights
  rank; msolve run + degree parse are external steps".
- `cycletype` "assumes squarefree reduction at sampled primes (generic; a
  non-squarefree prime shows as a repeated factor — discard that prime)".
- `nonres`: "certificate is per tested beta; genericity by sampling; qhull float
  facets are exactness-verified before use, unverifiable facets dropped
  (conservative)".
- LCLM factorization (`reducible = pf_order < holonomic_rank` is the ONLY
  reducibility statement) — stays out.
- Series-only inputs (no integrand polynomial) — use the series-route tool
  (`series_pf.py`, tools/annihilator) instead.
- `regions` certifies completeness only; it does not construct the series.
  Hidden (Glauber-type) regions are NOT guaranteed to be enumerated from U+F
  alone — for those the χ mismatch is the tell, not the facet diff.

INVOKE:
- CLI: `dipstick order '<Gstr>' <kin_csv> <x_csv> <svar> [Lmax]` · `dipstick count
  <npts>` (emits crit3⟨n⟩.ms in CWD; then run msolve through the capped wrapper
  member — section below — and parse
  the parametrization degree) · `dipstick cycletype [file.ms]` · `dipstick nonres` ·
  `dipstick audit [--fast] <file.ms>` · `dipstick flatchi [<file.ms> | --quick]` ·
  `dipstick regions <in.json> [more.json ...]` (several inputs share one Oscar
  load; equivalently `julia --project=<your Oscar project>
  regions/regions.jl <in.json>`).
- julia-side (order): `include("pf_rank.jl"); pf_probe(Gstr, kin, xvars, "s";
  Lmax=6)`. Direct entry points: `pf_probe_gen` (Oscar-poly interface; jet carried in
  GF(p)[x,s], svar evaluated AFTER differentiation) and `pf_probe_esc` (distrusts
  order==Lmax+1 saturation / prime mismatch, escalates Lmax automatically).
- Regression: `julia --project=<your Oscar project> regress_order.jl pf_rank.jl`.
- Battery: `python3 battery.py [--legs count,cycletype,nonres,regions]` — compares
  against the vendored `fixtures/` receipts; the `regions` leg runs
  `regions/selftest.jl` (also runnable alone: `python3 regions/selftest.py`).

OUTPUTS: `order` → pf_order per variable, holonomic_rank, reducible flag; reference
walls: sunrise→order 2 (rank 7), 3L-banana→order 3 (rank 15). `count` → msolve
parametrization degree = rank (X(3,6)→26 receipt). `cycletype` → per-prime
factor-degree multisets + transitivity/parity/Jordan summary (S26 receipt). `nonres`
→ facet count + resonant pairings + verdict (X(3,6): 43 exact-verified facets; 0/43
resonant at 3 generic rational points ⇒ irreducible; 10/43 AT THE LIGHT POINT —
light points sit on module resonance loci). `audit` → per-file dict: line census,
deg Res, arrangement valuations, crit count. `flatchi` → characteristic
polynomial + χ_top (projective) or bounded/region counts (affine). `regions` →
the enumerated regions (scaling vector v_R, n!Vol, χ_reg and G_R each), the two
sums vs the whole (`[TILES OK]` / `[CHI ADDITIVITY OK]` or the mismatch), and for
a user list `VERDICT: PASS (region list complete)` or `VERDICT: FAIL (region
list INCOMPLETE)` with the MISSING facets named by v_R plus the volume and χ
deficits.

ENV: Singular backend requires p<2^29 (order). `DIPSTICK_JULIA` — the julia
invocation (default: `julia` on PATH; e.g. `julia +1.10` to pin a juliaup
channel). `DIPSTICK_JULIA_PROJECT` — a Julia
project providing Oscar/Singular (the `order` verb needs it). `DIPSTICK_LEVIATHAN` —
path to Leviathan.jl (defaults: sibling tools/leviathan, then upgrades/leviathan in
this repo; used by `order` and `regions`). `DIPSTICK_SKIP_CHI=1` — `regions`
only: skip all χ_reg Gröbner work (volume certificate unchanged, χ printed NA;
skipped χ counts as unverified). For one release the `regions` member also
accepts the old names `EBR_LEVIATHAN` and `EBR_SKIP_CHI` as fallback aliases.
`DIPSTICK_MSOLVE` — the msolve wrapper (default: the in-package
`run_msolve_capped.sh` member). python: sympy/numpy/scipy +
python-flint (cycletype). `regions` needs julia ≥ 1.10 with the Oscar and JSON
packages (Oscar.jl is GPL-3.0-or-later and bundles GAP, Singular and polymake; it
is user-installed, not shipped here, and the MIT member only imports it at run
time).

GATES: 2-prime cross-check
built in and MANDATORY; Lmax saturated (order==Lmax+1) ⇒ rerun larger or use
pf_probe_esc — never report a saturated order. Battery (from a scratch cwd):
order (regress_order.jl vs pf_rank.jl, resolving Leviathan from
upgrades/leviathan) 7/7 PASS; count PASS; cycletype PASS; nonres PASS (stdout
byte-identical to the reference log); audit PASS (both routes crit=3 on the
vendored 4-generic-lines toy, field-for-field agreement); flatchi PASS (affine
χ=q²−4q+6 / bounded 3 / regions 11 on the same toy + the projective 42-cage
calibration); regions PASS (three shipped inputs in one julia process, ~50 s:
bubble large-Q → TILES OK / CHI ADDITIVITY OK / VERDICT PASS; large-mass toy
complete list → TILES OK / VERDICT PASS; large-mass toy region-I-only list →
VERDICT FAIL naming the missing v=[1,1,1] with vol deficit 3, χ deficit 2);
reference receipts ship in `fixtures/`;
mutation-control one leg per member. Before reading any `cycletype`
statistics: at each prime, the brute-force count of F_p-rational solutions
must equal the number of linear factors of the eliminant mod p, at ≥2 primes
(the per-prime honesty gate — the only defense against a garbled eliminant
upstream). The `order_shim` leg runs the same battery through the compat shim at
tools/pf_rank.jl. A battery leg whose external engine is absent (julia + Oscar
for order and regions, msolve for count/cycletype, scipy for nonres) SKIPs with
a named reason — by design from a cold clone, not breakage. For `regions` the
named skip prints the one-time setup command
(`julia --project=upgrades/leviathan -e 'import Pkg; Pkg.instantiate()'`,
which installs Oscar from the General registry; large).

FOOTGUNS:
- A hardcoded P_k=±k!·G_sᵏ jet is WRONG for deg_s(G)≥2 — the exact recursion is the
  only valid path there (auto-routed; do not force the fast path).
- Order==Lmax+1 or prime mismatch = distrust (escalate), not a result.
- Do not read `count`'s solution count as solving anything — it is the rank, at
  generic weights.
- Do not confuse the `cycletype` verb with the galois tool (the motivic coaction
  engine) — they are unrelated.
- F_q point-count interpolation FAILS for k=4-scale hyperplane arrangements even
  at audited-generic primes (mod-q flat degenerations) — the exact flat-lattice
  route (`flatchi`) is the only honest tool there. Small-prime rows (q≤7) are
  systematically degenerate in any point-count interpolation — a fitter must
  support exclusion + held-out validation.
- `regions`: a χ mismatch WITH a complete region list means some G_R is
  Newton-degenerate (factored / Landau-singular) — a separate diagnostic, NOT
  incompleteness; check degeneracy before hunting phantom regions (the shipped
  large-mass toy is exactly this case: TILES OK, VERDICT PASS, χ 4 vs 5). Treat
  any χ mismatch, and any skipped χ, as "unverified — value-gate before use".

## The regions member (`regions/regions.jl`, expansion-by-regions completeness)

PURPOSE (member): given the Symanzik data of a family (G = U + F as a string,
the kinematic symbols, the integration variables, and an integer weight saying
which kinematic symbol is the small parameter), enumerate the lower facets of
the Newton polytope of G in one dimension higher (last coordinate = the
small-parameter weight). Those facets ARE the regions of the expansion by
regions, each with its scaling vector v_R. Two independent sums certify a
region list: (vol) Σ_R n!·Vol(G_R) = n!·Vol(G), a rigorous tiling identity, and
(χ) Σ_R χ_reg(G_R) = χ_reg(G), the master-count additivity computed with the
Leviathan engine. A user-supplied list is diffed against the enumeration and
the MISSING region is named by its v_R. Scaling vectors may have negative
components (linear-programming convention, no δ shift).

INVOKE (member): `dipstick regions in.json [more.json ...]`, or directly
`julia --project=<your Oscar project> tools/dipstick/regions/regions.jl in.json`.
Input JSON: `{"G": str, "kin": [...], "x": [...], "kin_weight": {kin: int, ...},
"regions": [[v...], ...]?}` (omit `regions` to enumerate only). About 45 s per
call, dominated by loading Oscar; extra input files in the same call add
seconds, not another load. Shipped inputs: `regions/tests/` (bubble large-Q
three-region textbook case, large-mass toy complete / incomplete pair, sunrise
threshold enumeration). Control battery: `python3 regions/selftest.py` or the
package battery's `regions` leg.

LIMITS (member): certifies completeness, does not build the series; hidden
Glauber-type regions may not be enumerated from G alone (χ mismatch is the
only tell — do not certify such families on the volume sum alone); a χ
mismatch with a complete list = Newton-degenerate G_R (see FOOTGUNS);
`DIPSTICK_SKIP_CHI=1` drops the χ sum for very large families where the
global Gröbner basis dominates, leaving the volume certificate only. Formerly
the standalone `ebr-complete` package; env names `EBR_LEVIATHAN` /
`EBR_SKIP_CHI` are still accepted for one release.

## The capped-msolve wrapper member (`run_msolve_capped.sh`)

Runs msolve (Gröbner / F4) under a MANDATORY address-space cap. msolve's F4
matrix allocation can balloon to hundreds of GB with no internal limiter; the
cap makes it die alone (malloc failure) instead of taking the machine. There
is deliberately NO default-uncapped path. msolve itself is third-party —
obtain upstream (https://msolve.lip6.fr/). The wrapper is the ONE cap door for
every msolve exec in this repository: dipstick's `count`/`cycletype` verbs and
the landau-alphabet Gröbner escalation backend both route through it.

INVOKE (member): `MSOLVE_CAP_GB=<int> tools/dipstick/run_msolve_capped.sh
<msolve args...>` — e.g. `MSOLVE_CAP_GB=31 run_msolve_capped.sh -f sys.ms -o
sys.out -t 8`. msolve args pass verbatim; stderr carries the provenance lines
`[msolve-cap] ulimit -v ... / exec: ...`.

ENV (member): `MSOLVE_CAP_GB` — MANDATORY, plain positive integer GB. Suffixed
literals ("31G") are REJECTED. `ulimit -v` takes KILOBYTES; the script
converts KB = GB·1024·1024 and execs in the SAME shell so the limit binds
msolve itself. `MSOLVE_BIN` — optional binary path (default: msolve on PATH).

rc contract (member): 2 = refusal (no/invalid cap, missing binary); otherwise
msolve's rc verbatim; rc=139 (SIGSEGV) = HIT THE CAP — raise MSOLVE_CAP_GB,
the machine is fine.

GATES (member): calibrated on katsura-13: cap=1G dies in seconds with the
machine untouched; cap=16G completes.

FOOTGUNS (member):
- CHAR-P COEFFICIENT GARBLE: msolve (0.6.5) silently garbles char-p input
  coefficients larger than 2^31 — wrong parse, rc 0, no warning; the solved
  system is then NOT the emitted system. EVERY .ms emitter MUST reduce every
  coefficient mod p to [0, p) at emission and refuse out-of-range bytes — the
  wrapper does NOT re-check them. Integer-cleared exact-arithmetic systems are
  exactly the exposed class. After wiring the fix in any emitter, validate with
  one known-answer wrong-then-right probe pair plus one emitting-scale control
  before trusting any new verdict. (The msolve-independent `audit` verb is the
  cross-check for fiber counts.)
- msolve's fglm stage has segfaulted at primes p > 2^30 — prefer primes below
  that for fglm runs (msolve's own support ceiling is 31-bit primes, < 2^31).
- rc=139 is the cap, not a bug in your system — raise and rerun; don't debug
  the input first.
- In-process Oscar/Julia F4 (`groebner_basis_f4` = msolve F4): same law —
  apply `ulimit -v` to the julia process itself. msolve F4 / Singular nf are
  NOT thread-safe — never Threads; serial per process, parallelize with
  `julia -p N` + pmap only. Singular ops cap primes p<2²⁹.
- Uncapped python/julia wrappers around msolve re-create the
  runaway-allocation hazard — cap the OUTERMOST process that mallocs.

RELATED: upgrades/leviathan in this repo (χ_regulated engine shared by `order`
and `regions`; the main in-process msolve consumer — same cap law applies);
routes into GeoTriage / A-B-C value-fittability; `regions` sits upstream of any
region-based boundary/moment resummation, and its lower-facet enumeration
parallels the tropical/Newton machinery in tools/tropical-sampler. The wrapper mirrors the jemalloc memory-guard
pattern of the AMFlow.cpp fork (the sibling repository amflow-cpp).
CREDIT: `order` — Griffiths–Dwork reduction (Griffiths [Griffiths]; Lairez [Lairez]) and
the Lee–Pomeransky count of master integrals as critical points [LP]; `count` —
holonomic rank of Euler integrals as a signed Euler characteristic / ML degree after Huh
[Huh], Catanese–Hoşten–Khetan–Sturmfels [CHKS], Bitoun–Bogner–Klausen–Panzer [BBKP] and
Agostini–Fevola–Sattelberger–Telen [AFST], with χ-drop counting following Chestnov,
Crisanti and Giroux [Leviathans] (see upgrades/leviathan); `flatchi`/`audit` —
Orlik–Terao [OT] and Zaslavsky [Zas] for arrangements; `nonres` — GKZ systems
(Gelfand–Kapranov–Zelevinsky; Saito–Sturmfels–Takayama, classical); the X(3,6) controls
are CEGM integrals of Cachazo, Early, Guevara and Mizera [CEGM]; `regions` — expansion
by regions (Beneke–Smirnov school) and Newton-polytope region geometry, with χ-drop
counting per the Leviathan engine. Gröbner engines: msolve by Berthomieu, Eder and Safey
El Din [msolve], Singular. Bracketed keys resolve in REFERENCES.md at the repository
root.