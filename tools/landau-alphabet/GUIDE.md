# Landau Alphabet — guide

Tool page: https://bootloops.ai/tools/landau-alphabet.html

KIND: package

PURPOSE: implements the graph side of the Landau bootstrap method — bounded,
face-by-face Landau / symbol-alphabet derivation from a graph spec, with a
load-bearing completeness-honesty contract.
`linear_landau.py` extends the Landau equations to mixed quadratic + linear
(worldline/HQET/eikonal) propagators. The `LandauAlphabet.jl/` member is the package's
Julia bootstrap/regression engine — ansatz construction → exact ℚ constraints →
high-precision sampling → LLL/PSLQ fit with acceptance gates (see the section below).

USE-WHEN:
- Need a graph's symbol alphabet / Landau letters, bounded and auditable →
  `landau_alphabet.py` (per-face resultant chain, hard `--face-timeout`, incremental
  flush).
- A face times out under sympy → hand exactly that face to `pld_bridge.jl` (same input,
  same output schema; needs PLD.jl, obtain upstream).
- Want an INDEPENDENT numerical cross-check of a sympy-backend alphabet without
  PLD.jl/Oscar → `julia pld_bridge.jl SPEC.json OUT.json --fallback --check
  SYMPY_OUT.json` (HomotopyContinuation-only fallback; auto-engaged when PLD/Oscar
  fail to load, e.g. any Julia >= 1.11): per-face torus critical systems of the
  Lee-Pomeransky G=U+F on a seeded generic line; exit 2 + CROSS-VALIDATION CATCH
  when a witness matches no candidate letter; writes complete=false always
  (witnesses confirm letters, never certify completeness). Witnesses are interior
  pinches — they include second-type surfaces the DE alphabet prunes (box1l: s+t)
  and exclude endpoint-class letters (box1l: s, t), so pool candidates accordingly.
- Graph has linear/eikonal propagators → `linear_landau.py`: `symanzik_mixed`,
  `landau_locus_mixed`, `pm_alphabet_x`, `verify_letter` (saturated-Groebner spurious
  filter). Lemma used: linear-prop α's are absent from U and F stays ≤ quadratic in
  each.
- Any alphabet/analysis JSON claiming completeness → `audit_faces.py` before trusting it.

NOT-FOR: PLD.jl backend is heavy and unbounded — not the default. `linear_landau.py` is
the route for mixed propagators because PLD.jl's graph-builder can't build them.

INVOKE: `python landau_alphabet.py graph_specs/<g>.json --out <g>_alphabet.json
--face-timeout 45` · `PLD_ENV=<your PLD julia env> julia pld_bridge.jl
graph_specs/<g>.json <g>_alphabet_pld.json --method sym` · `python audit_faces.py
OUT.json` (or `--spec graph_specs/<g>.json`; `--glob 'results/**/COMPLETE_*.json'
--strict` sweeps a tree). API: `from landau_alphabet import load_spec, run;
log = run(load_spec(spec), out_path=..., face_timeout=45)`.
Heavy faces escalate automatically to the in-package msolve/Singular Groebner
member `pld_groebner_backend.py` when msolve or Singular is installed
(`PLD_GROEBNER=off` to disable; see the member section below).

INPUTS: graph-spec JSON: `edges` = [[v1,v2,mass²],...] (mass² number or symbol string;
bare [[v1,v2],...] = all massless), `nodes` (vertex → external legs), `kinematics`
{invariants, external_masses, channels (one of {S, complement} suffices)}.

OUTPUTS: `{complete, alphabet (canonical, spurious-pruned), first_entry, last_entry,
spurious, faces[...face_status enum...], honesty {expected_faces, attempted_faces,
missing_faces, resolved/bounded_partial/bounded_timeout, unresolved_faces,
global_deadline_hit}}`. Top-level `complete` is true ONLY if every expected face was
attempted AND none is bounded_partial/bounded_timeout AND no global deadline hit — read
it FIRST.

ENV: `--face-timeout` per face (hard bound); `PLD_ENV`/`PLD_SRC` for the PLD.jl bridge;
`PLD_GROEBNER` (auto/off/force) for the Groebner escalation.

GATES: `audit_faces.py` hard-fails exit 2 on any silent drop: complete-claim with
timed-out/pending faces, attempted<expected, tally mismatch, or a legacy COMPLETE_*
whose strata count is far below the true face lattice. Exit codes 0 clean / 2 silent
drop / 3 unverifiable under `--strict`. Status strings `resolved-groebner:*` and
`no-free-params` count as resolved.

FOOTGUNS:
- Disconnected-cut gap (KNOWN LIMITATION): the face lattice enumerates CONNECTED subgraphs only
  — on graphs with >4 vertices a 2-particle cut can be a disconnected subgraph and is
  silently skipped, dropping canonical letters of the `s−4m²` class. K₄-type graphs
  unaffected. Workaround: seed with the validated sub-family alphabet and read additions
  from the DE-connection denominators.
- Symbolic external masses: specs are sympified with a SHARED locals dict — an invariant
  name reused inside a mass string must resolve to the same Symbol. Current code does
  this; do not bypass `load_spec` with hand sympify.
- A short-timeout run honestly reports `complete: false` — do not "fix" that by dropping
  the hard faces; the hard faces are where the missing letters live.
- `verify_letter` (saturated Groebner) is the spurious filter for the mixed engine — a
  pm_alphabet_x letter that fails it does not go in the alphabet. It is LOAD-BEARING:
  the mixed engine's `_discriminant_letters` deliberately over-generates, so unverified
  harvest letters are CANDIDATES, not alphabet members.
- The bounded resultant chain harvests kinematic factors from every intermediate
  generator set and strips pure-parameter monomial content each round, so deep-face
  letters cannot be lost to resultant self-annihilation. This guard lives in BOTH
  engines: mixed (`linear_landau._discriminant_letters`; regression
  `test_letterloss_alpha_content_face`, policed by `verify_letter`) and main
  (`landau_alphabet.discriminant_letters`; selftest leg "letter-loss (main engine)",
  policed by `canonical_alphabet`/`prune_spurious` + the Groebner escalation).
- Known test flakiness: `test_alphabet_3PM_IY` is sensitive to Python hash randomization
  (symbol-ordering-dependent spurious deg-2 letter on some orderings). Run the test file with `PYTHONHASHSEED=1` for a reproducible 8/8.

TESTS: `python3 selftest.py` — the default battery, fast code checks with known
answers in a few seconds: full engine runs on the 1-loop box (exact alphabet
{s,t} + honesty contract) and the 2-loop sunrise (normal + pseudo thresholds),
incremental flush, the auditor's exit-2 silent-drop gates (clean log passes;
tally-lie and a legacy COMPLETE_* against the 54-face K4 census both exit 2 —
census from the spec alone, no K4 engine run), the forced Groebner escalation
on the sunrise leading face vs its known thresholds (named SKIP without
msolve/Singular), the mixed-engine 2PM-box exact alphabet with the
spurious-y refutation and the letter-loss guard, and the MAIN engine's own
letter-loss guard (box1l leading face must keep both s and t with escalation
masked off, plus the self-annihilation reproducer face — pure-sympy chain).
`python3 selftest.py --full` — the original engine-driving suites plus the
bridge fallback: `test_landau_alphabet.py` (K4 face engine + auditor; minutes
— the K4 graph lives in this tier), `PYTHONHASHSEED=1 test_linear_landau.py`
(mixed engine, 8/8; 2PM box, 3PM IY/H, Apéry γ=3 anomalous-threshold
reproduction), the pld_bridge.jl HC-only fallback on box1l (positive control:
candidates = the sympy letters plus the second-type s+t, every witness matches,
exit 0; negative control: the sympy letters alone, the s+t witness must be
CAUGHT, exit 2; named SKIP without
julia or HomotopyContinuation.jl), and the Julia engine suite below (named
SKIP if julia is not on PATH). The individual suites still run standalone
exactly as before.

## The Groebner escalation member (`pld_groebner_backend.py`)

Fast Gröbner backend for the heavy face discriminants the bounded sympy
resultant chain cannot do: per-face Landau letters via Rabinowitsch-saturated
Jacobian-ideal elimination (saturate against the Feynman-parameter torus,
eliminate the x_i, factor the elimination-ideal generators over the kinematic
ring). Benchmark (a 3-loop K4 degree-267 face): sympy >900 s timeout / PLD.jl
homotopy error at 837 s → msolve 32 s; 5-edge faces 1-5 s; finds letters the
sympy chain misses.

INVOKE (member): automatic — `landau_alphabet.py` escalates to it on
FaceTimeout or a PARTIAL/EMPTY chain result (`PLD_GROEBNER=off` disables,
`=force` uses it for every face; status strings become
`resolved-groebner:msolve` / `:singular`). Direct: `from pld_groebner_backend
import face_letters; r = face_letters(F, fvars=[...], kinvars=[s,t,m2],
backend="auto", threads=8, timeout=900)` → r["letters"], r["status"],
r["backend"], r["secs"], r["elim_gens"]. CLI: `python pld_groebner_backend.py
FACE_UF.json --field s,t,m2 --backend msolve --out OUT.json`. Auto-fallback
order: msolve (`-e ELIM -g 2` modular grevlex, PREFERRED) → Singular
eliminate() over Q (correct, slow — independent cross-check) → sympy groebner
lex (last resort). External deps: msolve and/or Singular; the sympy last-resort
leg needs only sympy.

CAP (member): every msolve exec routes through the capped wrapper
`../dipstick/run_msolve_capped.sh` — one cap door, never a raw exec.
`MSOLVE_CAP_GB` defaults to 30 (an explicit environment value wins); a wrapper
refusal surfaces as backend-unavailable and the Singular leg takes over.

GATES (member): the Singular leg is kept precisely to cross-check msolve on new
territory — run it once per new face class.

FOOTGUNS (member):
- msolve expects integer coefficients; the writer clears denominators before
  emission (a rational input such as "36*a1^2/77-1" is otherwise read as
  {36*a1^2-1}, denominator dropped, rc=0).
- msolve `-S` f4sat does not emit a printable elimination GB — the `-e ELIM -g
  2` form is load-bearing, don't "simplify" it.

## The Julia engine (`LandauAlphabet.jl/`)

The package's bootstrap/regression engine (Julia; MIT): ansatz construction as
weight-graded iterated-integral words over a pluggable kernel alphabet → exact
constraint assembly over ℚ (Nemo) → high-precision sampling (Arb balls; the
amflow-cpp CLI oracle under the staggered-goal consistency protocol — a digit is
accepted only if runs at goal G and G+30 agree on it) → LLL/PSLQ fit with
acceptance gates. Two kernel instances: MPL (dlog letters) and elliptic (Γ₁(6)
modular kernels), including exact Eisenstein cusp values over Q(√−3)
(`src/cuspvals.jl`) and integrable-symbol assembly for multivariate dlog
alphabets (`src/symbols.jl`).

USE-WHEN (engine): fitting an integral to a weight-graded ansatz against exact
constraints plus high-precision samples; fast cached GPL/iterated-integral
evaluation at scale (`ItIntEvaluator`/`GPLContext` — the series path
`tools/gpl-eval` re-exports); exact level-6 modular data (Eisenstein bases, cusp
values, Landau/cusp constraints).

DEPS (engine): Julia ≥ 1.11 with the registered packages Nemo, Arblib, JSON3
(the committed `Manifest.toml` is generated on Julia 1.11; on a newer Julia
minor, run `Pkg.resolve()` to re-pin before `Pkg.instantiate()`).
The `amflow.jl` sampling path shells out to `amflow_cli` — build it from
the AMFlow.cpp fork (the sibling repository amflow-cpp) and put it on PATH,
or point `AMFLOW_CLI` at the binary;
everything else is self-contained.

TESTS (engine): first run `julia --project=LandauAlphabet.jl -e 'using Pkg;
Pkg.instantiate()'`, then `julia --project=LandauAlphabet.jl
LandauAlphabet.jl/test/runtests.jl` (GPL evaluator vs classical polylogs,
basis-tower counts, exact modular identities, integrable-symbol counts, lattice
fit + integer relation on planted truths, Cauchy boundary constants) and
`julia --project=LandauAlphabet.jl LandauAlphabet.jl/test/test_cuspvals.jl`
(exact cusp-value acceptance, weights 2-4 at all four cusps). Both are also
run by `python3 selftest.py --full`.

RELATED: PLD.jl
(cross-validation backend, obtain upstream); upgrades/SOFIA.jl (independent blind
alphabet route, in this repo); upgrades/leviathan (χ-drop letters — independent third
route, in this repo); `tools/gpl-eval` (GPLEval.jl, built on the Julia engine's
series path).

CREDIT: the Landau equations are L. D. Landau's [Landau]; this package implements the
graph side of the Landau bootstrap developed by H. S. Hannesdottir, A. J. McLeod, M. D.
Schwartz and C. Vergu [LB1, LB2, LB3]. Principal Landau determinants and the PLD.jl
cross-validation backend are by Claudia Fevola, Sebastian Mizera and Simon Telen [PLD,
FMT]; the numerical fallback uses HomotopyContinuation.jl by Paul Breiding and Sascha
Timme [HCjl]; graph polynomials in Lee–Pomeransky form [LP]; Gröbner escalation by
msolve (J. Berthomieu, C. Eder, M. Safey El Din [msolve]) and Singular; the
blind-alphabet sibling is SOFIA by M. Correia, M. Giroux and S. Mizera [SOFIA] (see
upgrades/SOFIA.jl). LandauAlphabet.jl follows the modular-form and iterated-integral
conventions of Luise Adams and Stefan Weinzierl [AW] for the Γ₁(6) sunrise/kite kernels,
evaluates Goncharov polylogarithms [Gon] in the Vollinga–Weinzierl conventions [VW], and
fits with PSLQ [PSLQ] and LLL [LLL] on Nemo/Arblib [Nemo, Arb]. Bracketed keys resolve
in REFERENCES.md at the repository root.