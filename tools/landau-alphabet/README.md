# Landau Alphabet

Standalone, bounded, face-by-face Landau / symbol-alphabet engine.  A graph-agnostic engine with an optional PLD.jl cross-validation backend; the shipped graph specs (1-loop box through the 3-loop K4) double as its regression set.  The `LandauAlphabet.jl/` member is the package's Julia bootstrap/regression engine (see below).

## Two backends, one input

| backend | file | deps | exhaustive? | bounded? |
|---|---|---|---|---|
| **sympy** (default) | `landau_alphabet.py` | `sympy` only | per-face resultant chain | **yes** — hard `--face-timeout` per face, incremental flush |
| **PLD.jl** (optional) | `pld_bridge.jl` | Oscar + HomotopyContinuation | full principal A-determinant | no (heavy) |
| **HC fallback** (in `pld_bridge.jl`) | `pld_bridge.jl` (auto when PLD/Oscar fail to load, or `--fallback`) | HomotopyContinuation only | numerical witnesses, per-face | no |

The sympy backend works standalone; PLD.jl is **not** a dependency.  Both read
the same `graph_specs/*.json` and write the same output schema, so a face the
sympy engine times out on can be handed to PLD.jl for closure.

When PLD.jl/Oscar cannot load (the default state on Julia ≥ 1.11), the bridge
runs its **HomotopyContinuation-only numerical fallback**: per face, the
Lee-Pomeransky polynomial G = U + F is restricted to a seeded generic rational
line in invariant space and the torus critical system {G = 0, ∂G/∂x_e = 0} is
solved by homotopy continuation; the surviving t-values are numerical
discriminant **witnesses**.  With `--check SYMPY_OUT.json` each witness is
matched against the candidate letters (the file's `alphabet`, `spurious` and
per-face letter lists are all pooled) — an unmatched witness means **no
candidate covers that discriminant surface** and the bridge exits 2.  Two
scope notes.  (1) Witnesses are INTERIOR torus pinches, so they include
*second-type components* a DE-level alphabet legitimately prunes: on the
massless box the leading face witnesses the u-channel surface `s + t`,
which the sympy backend's `{s, t}` does not carry — add such letters to the
candidate file when the catch is expected.  (2) Endpoint-class letters
(overall kinematic coefficients of a face polynomial, like `s` and `t` on
the box's codim-1 faces) have no interior pinch and produce NO witness —
the sympy chain's `kin_factors` harvest covers those; the fallback checks
witnesses against candidates, never candidates against witnesses.  The
fallback always writes `complete: false`: it confirms letters, it never
certifies completeness.

## Input schema

```jsonc
{
  "name": "k4box3l",
  "edges": [[1,2,"m2"], [2,3,"m2"], [3,4,"m2"], [4,1,"m2"], [1,3,0], [2,4,0]],
  // nodes: vertex -> list of external-leg indices attached there
  "nodes": {"1": [1], "2": [2], "3": [3], "4": [4]},
  "kinematics": {
    "invariants": ["s", "t", "m2"],
    "external_masses": {"1": 0, "2": 0, "3": 0, "4": 0},
    // channels: (sum of legs in S)^2 as an expression in `invariants`.
    // Only one of {S, complement(S)} is needed; the other is auto-added.
    "channels": {"1,2": "s", "2,3": "t", "1,3": "-s-t"}
  }
}
```

`edges[i] = [v1, v2, mass²]` with `mass²` either a number (`0`) or a symbol
string (`"m2"`).  This is the same shape used by the amflow-cpp family JSON
(edge list + replacement table); a bare `[[v1,v2],...]` list is also accepted
(all-massless).

## Output schema

```jsonc
{
  "complete":    false,               // TOP-LEVEL honesty flag (see below) — read this FIRST
  "alphabet":    ["s", "t", "s + t", "s - 4*m2", ...],   // canonical, spurious-pruned
  "first_entry": ["s", "t", ...],     // codim-1 F-face / bare-invariant letters
  "last_entry":  [...],               // U-only (second-type / at-infinity) letters
  "spurious":    [...],               // pruned: U-only, F/U-cancelling, param-monomial
  "faces": [
    { "face_edges": [...], "n_edges": k, "loops_L": L,
      "letters_F": [...], "letters_U": [...],
      "status_F": "resolved|bounded-partial|bounded-timeout",
      "status_U": "...", "status": "...",
      "face_status": "done|timed_out|over_budget|scaleless",   // normalized enum
      "secs": ... },
    ...
  ],
  "honesty": {
    "complete": false,                  // mirror of top-level flag
    "expected_faces": E,                // TRUE face-lattice size (recomputed independently)
    "attempted_faces": A,               // A < E  =>  enumeration was TRUNCATED
    "missing_faces": E - A,             // never-attempted faces (silent-drop class)
    "n_faces": N, "resolved": r, "bounded_partial": p, "bounded_timeout": t,
    "unresolved_faces": [[...], ...],   // re-run these with PLD.jl or longer timeout
    "global_deadline_hit": false,
    "note": "..."
  }
}
```

### Honesty fields (load-bearing)

The top-level `"complete"` flag is **`true` only if EVERY expected face was
attempted (`missing_faces == 0`) AND none is `bounded_partial`/`bounded_timeout`
AND the global deadline was not hit.** It is the first thing a reader (or a
downstream script) sees, so a truncated run can no longer *read* as "covered
everything." `honesty.unresolved_faces` lists the exact sub-sectors that need
either a longer `--face-timeout` or a `pld_bridge.jl` run; `missing_faces`
catches the orthogonal failure mode where enumeration stopped early. Each face
also carries a normalized `face_status` enum so the auditor can key off it
without parsing prose. This is the contract that makes the bounded engine safe
inside the bootstrap: it never silently drops a Landau surface.

> **Why this exists.** Without the contract, a time-pressed run can record only a fraction of the kinematic faces and silently drop the hard faces that timed out — and the hard faces are exactly where the missing letters live.

## Guard: `audit_faces.py`

Re-runs the face audit on any output/analysis JSON and **hard-fails (exit 2) on
any silent drop**:

```bash
# audit an engine output (face census taken from the embedded graph spec)
python audit_faces.py k4box3l_alphabet.json

# audit a legacy hand-written 'COMPLETE_*' file against a graph spec
python audit_faces.py path/to/COMPLETE_LANDAU_VARIETY.json --spec graph_specs/k4box3l.json

# sweep a tree, fail on any unverifiable completeness claim
python audit_faces.py --glob 'results/**/COMPLETE_*.json' --strict
```

It detects: (a) a file claiming `complete` while a face record is
`timed_out`/`over_budget`/`pending`; (b) `attempted < expected` (enumeration
truncated); (c) a tally mismatch (`honesty.bounded_timeout` ≠ the actual count
of timed-out face records); and (d) a legacy hand-written `COMPLETE_*` whose
recorded strata count is far below the true face lattice. Exit codes: `0` clean,
`2` silent drop, `3` unverifiable-under-`--strict`.

## CLI

```bash
# sympy backend
python landau_alphabet.py graph_specs/k4box3l.json --out k4box3l_alphabet.json --face-timeout 45

# PLD.jl backend (same input, same output schema)
PLD_ENV=<your PLD julia env> julia pld_bridge.jl \
      graph_specs/k4box3l.json k4box3l_alphabet_pld.json --method sym

# HC-only numerical fallback (no PLD/Oscar needed): witness the sympy
# backend's letters on a generic line; exit 2 if any witness matches none
julia pld_bridge.jl graph_specs/box1l.json box1l_hc.json \
      --fallback --check box1l_alphabet.json --seed 0

# audit any output for silently-dropped faces (exit 2 on a silent drop)
python audit_faces.py k4box3l_alphabet.json
```

> `k4box3l.json` has **54** faces; the 4 hard faces (leading K4 + three 5-edge
> sub-boxes) need a long `--face-timeout` or the PLD backend, so a short-timeout
> run honestly reports `complete: false` rather than truncating.

## API

```python
from landau_alphabet import load_spec, run
gs  = load_spec("graph_specs/topbox.json")          # or load_spec(dict)
log = run(gs, out_path="topbox_alphabet.json", face_timeout=45)
log["alphabet"], log["first_entry"], log["honesty"]
```

Lower-level pieces (all graph-agnostic, all importable):

- `symanzik_subgraph(gs, face_edges)` → `{U, F0, F, L, ...}`
- `discriminant_letters(gs, P, fvars, timeout)` → `(letters, status)` bounded resultant chain
- `canonical_alphabet(raw)` → deduped, unit-normalized list
- `classify_entries(gs, alphabet, faces)` → `(first_entry, last_entry)`
- `prune_spurious(gs, alphabet, faces, (U,F))` → `(kept, spurious)`

## Shipped graph specs (regression set)

| file | graph | known alphabet |
|---|---|---|
| `box1l.json` | 1-loop massless box | `{s, t}` |
| `sunrise2l.json` | 2-loop equal-mass sunrise | `{s, m2, s-m2, s-9m2}` |
| `n8_dbox.json` | non-planar massless double box (Tausk) | `{t, s+t}` |
| `ladder3.json` | 3-loop planar triple-box ladder | `{s, t}` |
| `topbox.json` | 2-loop non-planar gg→Zγ 'topbox' (massive 4-ring + massless 3-chain, arXiv:2402.07311 routing) | topology selftest-asserted (non-planar, closed massive ring); no engine-derived alphabet is bundled |
| `k4box3l.json` | 3-loop K4 light-by-light crossed massive box | `{s, t, s+t, s-4m2, t-4m2, st-4m2(s+t)}` |

## Tests

```bash
python selftest.py         # default battery: fast code checks with known answers (~seconds)
python selftest.py --full  # the engine-driving suites below + the Julia engine suite
```

The individual suites, unchanged:

```bash
python test_landau_alphabet.py                 # face engine + auditor (K4 graph; minutes)
PYTHONHASHSEED=1 python test_linear_landau.py  # mixed quadratic+linear engine, 8/8
```

## The Groebner escalation member: `pld_groebner_backend.py`

The heavy faces the bounded resultant chain cannot do (high-degree leading
faces, faces where the chain silently under-resolves to an empty or partial
letter set) go to a real Gröbner elimination engine. The member forms the
Jacobian ideal ⟨P, ∂P/∂x₁, …⟩, **saturates it against the Feynman-parameter
torus** (Rabinowitsch trick: adjoin `t` with `t·∏xᵢ − 1`, removing the spurious
`xᵢ=0` branches), eliminates the `xᵢ`, and factors the surviving
elimination-ideal generators over the kinematic ring. The irreducible factors
are the letters.

Backends, in auto-fallback order: **msolve** (`-e ELIM -g 2` modular grevlex
GB — `-S` f4sat does not emit a printable elimination GB) → **Singular**
(`eliminate()` on the saturated ideal over ℚ; correct but slow — the
independent cross-check engine) → **sympy** `groebner(order='lex')` (last
resort). Every msolve exec routes through the capped wrapper
`../dipstick/run_msolve_capped.sh` (`MSOLVE_CAP_GB`, default 30 GB) — one cap
door, never a raw exec.

`landau_alphabet.py`'s `discriminant_letters()` escalates to the member
automatically on FaceTimeout or a PARTIAL/EMPTY chain result; the status
string becomes `resolved-groebner:msolve` (or `:singular`) so the honesty
record shows which faces used the heavy engine. Controls: `PLD_GROEBNER=off`
(disable; pure-sympy fallback with honesty-flagged gaps), `=force` (every
face — cross-check), `=auto` (default).

```python
from pld_groebner_backend import face_letters
r = face_letters(F, fvars=[x1, x2, x3], kinvars=[s, t, m2], backend="auto",
                 threads=8, timeout=900)
# r["letters"], r["status"], r["backend"], r["secs"], r["elim_gens"]
```

```bash
python pld_groebner_backend.py FACE_UF.json --field s,t,m2 --backend msolve
```

Honesty caveat: the generators are the raw principal-determinant components on
the saturated torus — genuine elimination-ideal components may include
second-type / non-physical-region branches. The letters here are the COMPLETE
Landau variety; first-/last-entry and physical-region pruning remain the
caller's job (`classify_entries`, `prune_spurious` in `landau_alphabet.py`).

## The Julia engine: `LandauAlphabet.jl/`

The package's bootstrap/regression engine (Julia; MIT license in the member
directory): weight-graded iterated-integral ansätze over a pluggable kernel
alphabet, exact constraint assembly over ℚ (Nemo), high-precision sampling (Arb
balls; amflow-cpp CLI oracle under the staggered-goal consistency protocol), and
LLL/PSLQ fits with acceptance gates.  Two kernel instances — MPL (dlog letters)
and elliptic (Γ₁(6) modular kernels, with exact Eisenstein cusp values over
Q(√−3)).  Its fast suffix-cached series evaluator (`ItIntEvaluator` /
`GPLContext`) is the path `tools/gpl-eval` (GPLEval.jl) builds on.

```bash
# one-time: resolve the registered deps (Nemo, Arblib, JSON3)
julia --project=LandauAlphabet.jl -e 'using Pkg; Pkg.instantiate()'

# test suites
julia --project=LandauAlphabet.jl LandauAlphabet.jl/test/runtests.jl
julia --project=LandauAlphabet.jl LandauAlphabet.jl/test/test_cuspvals.jl
```

The `amflow.jl` sampling path shells out to `amflow_cli` (build from
the AMFlow.cpp fork, the sibling repository amflow-cpp; put it on PATH or
set `AMFLOW_CLI`); everything else is
self-contained.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).
