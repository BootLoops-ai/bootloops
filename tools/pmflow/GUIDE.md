# pmflow — GUIDE

Tool page: https://bootloops.ai/tools/pmflow.html

KIND: script. (`gravityflow.py` in this dir is
a compatibility alias for the same CLI.)

Files: `pmflow.py` CLI, `gf_solve.py` solver core (imported, not shelled out — do not
modify), `gf_eps0_ratios.py`, `DESIGN.md`.

## PURPOSE

Self-consistent auxiliary-mass-flow CLI for cut-eikonal (PM) integral families where
the standard AMFlow boundary recursion is a FIXED POINT (cut props are forbidden
η-positions; η→∞ rescaling reproduces the family; Cutkosky declines −2k·u cut
"masses" ≠ phase volume). Pipeline detect → discover → respond → solve → inject → map
closes the boundary system that the stock boundary recursion does not close for
this class. Plus ε⁰ ratio salvage +
exact-rational identification for such runs.

## USE-WHEN

- AMFlow boundary recursion descends `_b0_b1_b2…` with IDENTICAL master counts per
  level on a cut-eikonal family → η-flow fixed point, STRUCTURAL — no ending_scheme
  fixes it. Diagnose with `pmflow detect`; then run the full pipeline. The engine's
  self-similarity guard aborts such a plain run with the structured
  `GRAVITYFLOW_FIXEDPOINT` code (a boundary family reproducing its parent's
  propagator set), which `detect` reads as the engine-certified FIXED_POINT
  verdict; the name-heuristic verdict remains the fallback for logs without it.
- Fixed-point/forced-ending run needs ε⁰ physics out of wild-factored dumps →
  same-sector ratio salvage: solved ε-families carry sector-block-common
  wild factors (dlnf/dε 1e5–1e9, node-erratic); NO ε-fit works on raw values;
  same-sector master ratios cancel the factor EXACTLY per node → smooth →
  Lagrange-extrapolate to ε⁰.
- Value-fit PSLQ nulls at high precision on such a family → `gf_eps0_ratios.py`: test
  relation rank / rationality of the master ratios BEFORE concluding transcendence.
- Cut-quadratic (graviton-line-cut) family aborting `GRAVITYFLOW_CUTREGION`: there is
  no eta→inf region family to probe — `discover --etac ''` runs one normal-path probe
  (the guard only fires there; forced ending depth≥1 bypasses build_boundary),
  records the family as a kind=cutregion MISS and switches itself to the depth-0
  leaf-Trivial surface (AMFLOW_FORCE_ENDING_DEPTH=0); the keys it then finds
  (kind=cut) are the undeformed cut masters — values from closed-form providers,
  injected via `inject`.

## NOT-FOR

Ordinary families where AMFlow terminates (use amflow-cpp directly). Per-point
probes are multi-process (one subprocess per probe); there is no in-process batched
injection mode. Ratio-salvage ceiling = dump precision − ~13d extrapolation
Lebesgue — not unlimited.

## INVOKE

`python tools/pmflow/pmflow.py <subcommand>`: `detect`, `discover`, `respond`,
`solve`, `inject`, `map` — full flag documentation in the module docstring and
`--help`. `gf_eps0_ratios.py` turns AMFLOW_DUMP_EPS_GRID logs into verified rational
tables; its `--selftest` replays a reference log (not shipped — set
`$GF_SELFTEST_LOG` to a local copy).

Environment (no built-in paths): set `AMFLOW_CLI`
(amflow-cpp CLI; build it from the fork in the sibling repository `amflow-cpp`, `../amflow-cpp`), optionally
`AMFLOW_JEMALLOC_WRAP` (wrapper script; unset = run the CLI bare), `FERMATPATH`
(fer64), `AMFLOW_IBP_CACHE` (cache dir) — or pass `--amflow-cli/--jemalloc-wrap/
--cache` per call. BASE_ENV is DEFAULTS-ONLY (caller's environment wins over every
entry; empty defaults are not exported).

## GATES

validated on a 4PM cut-eikonal reference family (15 propagators, 43 masters)
at 2 kinematic points, 43/43 A_alpha-row closure, worst
residual ~5e-52 (re-validated at NEPS=18 grids); the two flows are mutually
consistent (contour-convention question resolved, DESIGN.md §5) — joint residual is
the check to reproduce on a new family. gf_eps0_ratios honesty: honest digits =
full-set vs inner-subset agreement; PSLQ cap 1e19; held-out = verified −
(2·height+2); ID flagged WEAK if verified <90d.

## FOOTGUNS

- Alias diagnostic: amflow `result` Laurent coefficients requiring
  multi-order cancellation to reproduce sampled |f| (e.g. |c0|~1e31 vs |f|~6e4) are
  NOT integral values — certified balls are fit-internal; symptom: coefficient
  magnitudes grow ~1/ε₀ per order. Use the EPS_GRID per-node samples instead.
- Never fit ε-dependence of raw fixed-point master values — the sector wild factor
  makes every fit garbage; ratios first, always.
- `map` without `--foreign-as-j63` leaves unmatched keys unmapped by design —
  silence there is a wiring error, read the printed list.
- Subset-stability of the Lagrange extrapolation = the honest-digit count; a
  single-grid extrapolation number is not a digit claim.
- Omitting `--etac` on discover/respond uses the DEPRECATED legacy 15-slot default: on
  any non-15-prop family the engine aborts with `AMFLOW_ETAC_OVERRIDE: length` —
  discover turns that into a loud rc=1 error, never a silent
  "discovery complete" with 0 keys. Pass `--etac CSV`, or `--etac ''` for the
  engine-chosen etac (cut-quadratic families).
- The discover MISS regex is dual-surface (key BEFORE and AFTER "has no
  Vacuum entry") — a single-shape regex is blind on one binary generation's
  logs; key lists from a single-shape parse are suspect.
- `_env` is passthrough-with-defaults: the caller's environment wins. A
  `{**os.environ, **BASE_ENV}` merge would silently discard shell-exported
  vars at the binary.
- `map` matches used-prop CONTENT — a match that never compares content
  cannot produce a nonempty circular_map on foreign ending families.
- There is no default IBP cache path — the cache comes from
  `AMFLOW_IBP_CACHE`/`--cache` only.
- `solve`/`gf_solve.py` default anchors are the 15-slot reference family's
  (`63:J63,111:J63`): on any other family the solve refuses with a named
  `--anchors` error instead of running — pass `--anchors "SECTOR:FORM,..."`
  with the family's own corner sectors (bit i of SECTOR = index slot i;
  built-in FORM: `J63`).

## What runs without an engine build (battery: `python3 selftest.py`)

- `python3 selftest.py` — CLI surface (`--help`, subcommand parse, the solve
  subparser's `--anchors`/`--n-eps`); ratio salvage on a synthetic wild-factor
  eps-grid log, recovering planted rational ε⁰ ratios exactly; the clean
  fail-closed FileNotFoundError on the gf_eps0_ratios selftest without
  `$GF_SELFTEST_LOG`; an anchored closure solve on a synthetic 2-master
  fixed-point fixture (`--anchors "1:J63"` reproduces the planted J63 / 2·J63
  solution with an all-consistent row census) plus the named `--anchors`
  refusals; and the `GRAVITYFLOW_FIXEDPOINT` surface (regex↔emit-literal
  coherence with the engine source, plus a `g++ -fsyntax-only` gate on the
  patched `amfsystem.cpp` — the engine tree is read from `$AMFLOW_CPP_SRC`, else
  from the sibling checkout `../amflow-cpp` next to this repository; both skip
  by name when the engine tree or toolchain is absent).
- Engine legs need an amflow-cpp build (from the sibling repository `amflow-cpp`); the validation record
  above (the 4PM reference family, 43/43 closure) is a recorded validation, not
  rerun here.


CREDIT: pmflow is a driver around the auxiliary-mass-flow method of Xiao Liu, Yan-Qing
Ma and Chen-Yu Wang [AMF1, AMF2] as implemented in AMFlow [AMFlow] and AMFlow.cpp (the
AMFlow.cpp contributors [AMFcpp]; our fork is the sibling repository `amflow-cpp`); the η→∞
boundary recursion, ending schemes and vacuum boundary conventions are theirs. The
fixed-point detection and the self-consistent respond/inject loop for cut-eikonal
(post-Minkowskian) families — the integral class of e.g. Driesse, Jakobsen, Mogull,
Nega, Plefka, Sauer and Usovitsch [DJMNPSU] — are this package's addition. Bracketed
keys resolve in REFERENCES.md at the repository root.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).
