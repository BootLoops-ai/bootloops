# SubTropica — guide

Unofficial Julia rewrite of SubTropica's front end (M. Giroux, S. Mizera, G. Salvatori,
arXiv:2604.20954, MIT); not affiliated with or endorsed by the SubTropica/HyperFLINT
authors — cite the original.

SubTropica — `tools/subtropica` — https://bootloops.ai/tools/subtropica.html

KIND: package (Julia driver `src/` incl. `hf_bridge.jl`,
`lr_refine.jl`, `graph.jl`, `types.jl`; `scripts/fibrate/fibrate.py`; specs
`DESIGN.md`/`DESIGN_B1.md`/`CONTRACTS.md`). The HyperFLINT engine (C++17/FLINT
hyperlogarithm CLI port) is built separately and is NOT vendored — point
`SUBTROPICA_HF_BIN` at your built engine wrapper and `SUBTROPICA_MZV_DATA` at its MZV
reduction table. The optional in-process C-ABI lib `deps/libhf_cabi_static_f2.so`
(a 56 MB binary) is not included, and neither is the engine patch it must be
built with (we observed heap corruption when reusing the v1.2.8 library
in-process in our build; the patch is a per-call PolyCtx hygiene fix, see
`deps/BUILD.md`; the CLI transport is unaffected and is the supported route): `transport=:cabi` is therefore available only to users who
have such a patched build (`deps/build_cabi.sh` refuses loudly until
`HF_CABI_WORK` and `HF_FLINT_PKGCONFIG` are set). The CLI transport is the
default and the supported route.

## PURPOSE

No-Wolfram SubTropica pipeline over the HyperFLINT C++17/FLINT hyperlogarithm engine:
exact ε-layer Laurent coefficients / boundary constants in the MZV ring (Log2 + 9
irreducible MZVs) for Euler/Feynman parametric integrals — convergent (tier A),
log-divergent tropical subtraction (B1), and power-divergent / multi-regulator
Nilsson–Passare continuation (B2). Replaces quadrature+PSLQ for this class; Laurent
poles are DERIVED, not hand-subtracted.

## USE-WHEN

- Polylog sector with a known parametric rep — want symbolic ε-layer constants instead
  of another quadrature+PSLQ → `subtropica_integrate` (or the bare hyperflint CLI for
  one-off convergent calls).
- DIVERGENT parametric integral (poles in ε; log-divergent facets; power-divergent
  rays; rays needing a second regulator) → subtropica census routes it: log-div → B1 facet
  counterterms; power-div/GP-violating → typed refusal naming
  `allow_continuation=true` → B2.
- Blocking letter in an LR search — real obstruction or spurious-Fubini artifact? Or an
  order-verify farm about to grind a big order space → `lr_refine` FIRST
  (subset-lattice DP, seconds, zero engine calls, adjudicates ALL integration orders at
  once; 10²–10⁴× cheaper than engine order-verifies).
- Boundary constants currently PSLQ-fitted and the integrand has a convergent
  parametric rep → derive them symbolically; PSLQ demotes to held-out gate.
- Graph in hand, not an Euler quadruple → `graph.jl` (`symanzik_UF`,
  `graph_to_quadruple`, `write_quadruple_json`; exact U,F parity vs pySecDec).
- Per-face/per-word farms with thousands of small engine calls per worker → C-ABI
  transport `transport=:cabi` (24–70× small ops, 176–4497× medium integrator ops vs CLI
  spawn floor).

## NOT-FOR

Numerators = typed refusal (graph front-end charter scope). Wm/Wp algebraic-letter
results are NOT MZV-closed (symbolic coefficients carry explicit conjugate period-pair
atoms; `strict_periods=true` default still refuses typed — flip only when the
numeric-leg verification path is the intended gate). LR-wall targets refuse/time out
honestly — that is correct behavior, run `lr_refine` first. C-ABI does not move genuine
engine-compute walls; `parse_expr`/`eval`/`period` ops are CLI-only.

## INVOKE

driver: `subtropica_integrate(integrand; order, allow_continuation, transport)` in the
`tools/subtropica` Julia project. Suite: `ulimit -v 32505856; julia --project=tools/subtropica
test/runtests.jl`. Battery entry: `python3 selftest.py` (in this dir) — reference
sha pins + the fibrate control battery run engine-free; the Julia suite runs only
when `SUBTROPICA_HF_BIN` points at a built engine, and skips loudly otherwise. Bare engine two-step: `op=hyperflint` (integrand `expr`/`f`, ordered
`vars_int`, spectator kinematic vars allowed) → regulator wordlist; `op=evaluate_periods`
→ named MZV symbols. LR search: `find_lr_orders`/`find_lr_orders_scan`. 69 ops total;
per-op schemas inline — grep `// Request:` in the engine source
`SubTropica/HyperFLINT/bridge/cli/main.cpp` (part of the engine tree,
not this one). lr_refine surface: `refine_letters`, `is_blocking_real`,
`blocking_certificate`. Sound fibration: `scripts/fibrate/fibrate.py`. C-ABI rebuild:
`deps/build_cabi.sh` + `deps/verify_cabi.jl` (all checks must pass).

INPUTS: Euler/Feynman parametric integrands (Euler quadruple JSON or graph via
front-end; front-end emits FLAT-measure ν — the driver's `_b1_dlog` bridges, nu MEASURE
LAW in `types.jl`); always pass `mzv_data_path` explicitly.
OUTPUTS: exact symbolic Laurent coefficients in the MZV ring (or typed refusal naming
the required input fix); B2 emits pole towers + finite parts with Wm/Wp atoms where algebraic.
ENV: `ulimit -v 32505856` per process (engine AND suite). `HF_PERIOD_TUPLES=0`
(defaults ON → opaque process-local `Period[id]` atoms with no name table).
`SUBTROPICA_RUNS` overrides the default `run_dir` root for all four driver entry points
(default: `./subtropica_runs` under cwd). `scripts/env.sh` is a template — edit the
paths for your install. Oscar env rule: SELECTIVE `using Oscar:` only; never unpin
Hwloc/Parsers/libcxxwrap jlls (33-min precompile cache at stake).

## GATES

Measured: tier-A synthetic-truth 59–60d + Smirnov 37d vs second-exact-reduction
oracle; B1 pole tower exact + refusal battery 5/5 + micro-fixtures 76–78d; B2 pole
tower exact symbolic + finite part 46d vs published closed form THROUGH THE DRIVER
ALONE; every gate mutation-controlled. Rule for every closure: `lr_refine` first on
every integrand + ≥30d numeric gate vs an independent evaluator. lr_refine NO verdicts
require the positive control beside them (smirnov_tst2 5/5 in 7.6 s); md=1 =
certificate class, md=2 = PROXY for the algebraic tier (strong evidence, NOT a
certificate).

## FOOTGUNS

- On hard multi-variable integrands (measured on a 7-var 3-loop benchmark)
  the LR-order search does not terminate at any measured budget; the tool
  returns a typed refusal, never a silent number (the sector-decomposition
  wall moves to the LR-search stage).
- Silent 0 on divergent input — set `"check_divergences": true` on any integrand not
  already known convergent (~2× that step).
- Gate on JSON fields, NEVER rc: rc=0 on failed/divergent/malformed JSON. Require
  presence of `"result"` AND absence of `"failed"`/`"divergent"`/`"error"`.
- NOLR input ⇒ present-"result" SILENT-WRONG: on non-linearly-reducible integrands
  (irreducible-quadratic/algebraic-root letters) the engine returns rc=0 + clean JSON +
  a result that silently DROPS the algebraic-root branches (rel err O(1) measured).
  Field-gating does NOT catch this — only the lr_refine-first +
  independent-≥30d-gate rule does. On a certified-LR locus the same kernel is a valid
  crosscheck oracle (gated 69.6–69.9d).
- `fibration_basis` is UNSOUND on on-contour and reg-divergent ZIP words (measured: 631/635
  silently wrong while ALL step-1 regulators were exact) — MANDATORY
  per-(word,monomial)-group numeric validation vs step-1 semantics at ≥2 kinematic
  points on ANY fibration_basis output; route failures through the sound `fibrate.py`
  differential-Goncharov engine.
- `evaluate_periods` cannot re-ingest Pi/I/delta[·] coefficients (on-contour
  bookkeeping) — do not round-trip those terms.
- `Rat::parse` misses top-level '/' in outer-paren coef strings (the very format the
  CLI emits) — normalize to `num/(den)` before re-feeding wordlists.
- Threads on a shared engine handle CRASH (measured; the bridge
  serializes as a guard) — farm unit = one PROCESS per worker, kill by PID; timeouts
  are per-process. C-ABI rules: CLI stays default + validation oracle, keep 5–10%
  dual-run.
- B2 refusal rule: never flip `allow_continuation` globally — refusals are the default
  surface; the flag is per-call.
- stderr `hf_rat_split_verify:` counters are normal, not errors.
- `scripts/fibrate/` is the home of the fibration trio
  (`fibrate.py`/`zip_num.py`/`gpl_num.py`) — import from here. `zip_num.py` series
  length is ADAPTIVE to min mapped-letter magnitude (a fixed length is a
  dps-independent ceiling).

## What ships / what runs

The package ships the full Julia driver `src/`, `test/` + `fixtures/` + `reference/`
in-tree, `scripts/` (env template, fixture builders, fibrate trio), specs
(`DESIGN*.md`/`CONTRACTS.md`), and the `deps/` build recipe. The 56 MB prebuilt C-ABI
binary is not included (build per above).

Nothing runs end to end without a HyperFLINT engine build — the test suite needs
`SUBTROPICA_HF_BIN`/`SUBTROPICA_MZV_DATA` plus the Julia project env. Test files'
src includes are repo-relative. Legs reading run receipts (raycheck fixtures, probe
banks, pilot regression memo) refuse or skip loudly until you point
`SUBTROPICA_RAYCHECK_DIR` / `SUBTROPICA_PROBES_DIR` / `SUBTROPICA_RUNS` at your own
copies. The `reference/` files are sha-pinned (`reference/SHA256SUMS`, checked
by `selftest.py`) and the `fixtures/` files are read by the tests as they are —
do not edit either by hand (fixtures regenerate via `scripts/make_fixtures.py`).

ENGINE: obtain HyperFLINT/SubTropica upstream (MIT) and build it; the optional
in-process library is described in `deps/BUILD.md`.

LICENSING: this package is MIT-licensed with one exception: `src/lr_refine.jl` and
`test/test_lr_refine.jl` are a Julia transliteration of the compatibility-graph
reduction procedures of Erik Panzer's HyperInt (Copyright (C) 2014 Erik Panzer,
GPL-3.0-or-later; arXiv:1403.3385), following F. Brown (arXiv:0910.0114), and are
therefore distributed under the GNU General Public License v3.0 or later
(`LICENSE-GPL-3.0`). Because the module includes that file, SubTropica as loaded is a
combined work under GPL-3.0-or-later terms; every other file stands alone under MIT,
and removing `src/lr_refine.jl` (and its include/export lines in `src/SubTropica.jl`)
gives an MIT-only build. See `NOTICE` beside this file.

DEPENDENCY LICENSES: apart from that file the package is MIT-licensed; `src/` is translated from
`SubTropica.wl` (MIT, (c) 2025-2026 S. Mizera, M. Giroux, G. Salvatori; the
upstream notice is in `NOTICE` beside this file) and `scripts/make_fixtures.py`
imports pySecDec (GPL-3.0) at run time to regenerate one fixture. It drives, at run time and
without vendoring them, Oscar (a GPL-licensed Julia package, used in-process for
the polytope layer) and the ginac-based `ginac_gpl` evaluator of tools/gpl-eval
(GiNaC/CLN are GPL); the vendored upstream sources in `reference/` are MIT (see
`reference/LICENSE`, `reference/LICENSE-HyperFLINT` and `reference/PROVENANCE.md`).
`deps/build_cabi.sh` links LGPL FLINT/GMP/MPFR statically into a library that is
not distributed here.

## CREDIT

CREDIT: SubTropica and its HyperFLINT engine are by Mathieu Giroux, Sebastian Mizera and
Giulio Salvatori [SubTropica] (MIT; github.com/SubTropica/SubTropica) — this package is
an unofficial Julia rewrite of their front end and we are grateful to them for releasing
the code; please cite their paper. HyperFLINT reimplements Erik Panzer's HyperInt
[HyperInt], whose algorithms (fibration bases, compatibility-graph linear-reducibility
analysis) this port follows, building on Francis Brown's theory of hyperlogarithm
integration and linear reducibility [Br09a, Br09b]; `src/lr_refine.jl` is transliterated
directly from HyperInt's compatibility-graph reduction procedures and carries Panzer's
copyright and the GPL-3.0-or-later license (see LICENSING above). The MZV reduction table descends
from the MZV Data Mine of Blümlein, Broadhurst and Vermaseren [BBV]. The divergent
routes follow the tropical/Newton-polytope analysis of Arkani-Hamed, Hillman and Mizera
[AHM] and the Nilsson–Passare / Berkesch–Forsgård–Passare continuation of Euler–Mellin
integrals [NP, BFP]; graph polynomials in Lee–Pomeransky form [LP]. Polytopes by OSCAR;
exact arithmetic by Nemo/FLINT [Nemo, Arb]; GPL numerics by GiNaC [GiNaC, VW]. Bracketed
keys resolve in REFERENCES.md at the repository root.

## License

MIT License, except `src/lr_refine.jl` and `test/test_lr_refine.jl`, which are GPL-3.0-or-later (derived from HyperInt, Copyright (C) 2014 Erik Panzer; see `NOTICE` and `LICENSE-GPL-3.0` in this directory). Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).
