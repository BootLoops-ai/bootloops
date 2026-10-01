# HyperFLINT

A C++17 reimplementation of **HyperIntica** (Mathematica) and
**HyperInt** (Maple), the Feynman-integral-friendly hyperlogarithm
integration packages of Panzer (arXiv:1403.3385), built on FLINT's
`fmpq_mpoly` polynomial infrastructure.

HyperFLINT integrates hyperlogarithms and searches linear-reducibility
orders with a 3–10× end-to-end speedup over the reference engines on
realistic Feynman diagrams, and cross-validates against HyperIntica
and, where available, HyperInt to 30+ digits. For performance numbers,
the correctness scope, and the algorithm, see the SubTropica paper
(arXiv:2604.20954).

## Using it from SubTropica (most users)

Install the add-on paclet once:

```mathematica
PacletInstall["https://subtropi.ca/HyperFLINT.paclet"]
```

Since SubTropica v1.2.2, an installed HyperFLINT is picked up
automatically and becomes the **default** symbolic integrator and
LR-order search engine — a plain `STIntegrate[...]` call already uses
it. Explicit selection still works:

```mathematica
Needs["SubTropica`"];

(* one-mass box through HyperFLINT (the default when installed) *)
STIntegrate[{{{{1,2},0},{{2,3},0},{{3,4},0},{{1,4},0}},
             {{1,M},{2,0},{3,0},{4,0}}}, Order -> 0]

(* force a specific engine for one call *)
STIntegrate[..., "Integrator" -> "HyperIntica"]    (* built-in engine *)
STIntegrateHF[...]                                 (* alias: HF for both jobs *)

(* check what is active *)
SubTropica`$HyperFLINTAvailable          (* True if binary + data resolve *)
SubTropica`$STHyperFlintUseLibraryLink   (* True = in-process dylib engaged *)
```

The add-on bundles its native libraries; nothing else to install.
Supported platforms: macOS arm64/x86_64 and Linux x86-64.

## Standalone build

Requires FLINT ≥ 3.4, `libomp` for OpenMP parallelism, and a C++17
compiler.

```bash
brew install flint libomp              # macOS
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release \
      -DOpenMP_ROOT=/opt/homebrew/opt/libomp
cmake --build build-release
```

For ASAN- or TSAN-instrumented builds:
```bash
cmake -S . -B build-asan  -DCMAKE_BUILD_TYPE=Debug -DHF_ASAN=ON
cmake -S . -B build-tsan  -DCMAKE_BUILD_TYPE=Debug -DHF_TSAN=ON
```

## Command-line examples

Two front-ends: a **pretty CLI** for human inspection and a **JSON
CLI** (`eval-json`) that reads a request on stdin and writes a JSON
response on stdout (the interface SubTropica and the cross-validation
harness use).

Polynomial factoring, pretty CLI:

```bash
$ build-release/hyperflint factor 'x^2 - 1'
1 * (x + 1) * (x - 1)
```

Partial fractions of a rational integrand, JSON style:

```bash
$ echo '{"op": "partial_fractions",
        "f": "1/(x*(1 + x)*(x - y))",
        "vars": ["x", "y"], "var": "x"}' \
  | build-release/hyperflint eval-json
{"op":"partial_fractions","var":"x","polynomial_part":"0",
 "poles":[{"pole":"0","multiplicity":1,"coefs":["-1/y"]},
          {"pole":"-1","multiplicity":1,"coefs":["1/(y + 1)"]},
          {"pole":"y","multiplicity":1,"coefs":["1/(y^2 + y)"]}],
 "vars":["x","y"]}
```

A full integration request (`"op": "hyperflint"`) takes the integrand
expression, the integration variables in order, and the kinematic
variables, and returns the eps-expanded hyperlogarithm result; the
per-op request schemas (some forty ops, from `factor` through
`integration_step`) are documented inline in `bridge/cli/main.cpp`:

```bash
$ OMP_NUM_THREADS=13 build-release/hyperflint eval-json < request.json
```

## Parallelism

The top-level `hyperflint` driver parallelizes the outer entry loop in
`integration_step` via OpenMP; control it with `OMP_NUM_THREADS`. At
one-less-than-core-count the speedup over serial is typically 5–17× on
Log²-class integrands. `-DHF_OPENMP=OFF` at configure time gives a
strictly-serial build.

## Divergence detection (opt-in)

HyperFLINT defaults to assuming the input integrand is convergent
(matching the fast HyperIntica path). Set `"check_divergences": true`
in an `integration_step` request to enable boundary-divergence
scanning: every `{logpower, power}` pole bin is tested symbolically,
and any non-zero residue raises `HyperFLINTDivergentIntegral` naming
the specific `Log[var]^k / var^n` combination and the boundary
(`zero` or `infinity`). The check roughly doubles the cost of that
integration step; convergent inputs pay nothing when the flag is off.
From SubTropica, the `"CheckDivergences"` option of `STIntegrate`
drives this flag across all integrator backends.

## Layout

- `include/hyperflint/` — public headers
  - `core/` — `Poly`, `Rat`, `SymCoef`
  - `algebra/` — shuffle, partial fractions, conversions, derivatives,
    algebraic letters (Wm/Wp)
  - `series/` — Hlog/Mpl series, Laurent expansions
  - `integrator/` — TransformWord, primitive, integration_step,
    hyperflint driver
  - `reduce/` — MZV reductions, periods, BreakUpContour
  - `symbols/` — Word, Letter
- `src/` — implementations (mirrors `include/`)
- `bridge/cli/` — command-line interface
- `bridge/librarylink/` — the in-process Mathematica dylib
- `data/` — MZV reduction table (`mzv_reductions.json`)

## Cross-validation

Every primitive is validated as a three-way comparison: the HF
implementation runs each fixture alongside HyperIntica (Mathematica)
and, where a direct equivalent exists, HyperInt (Maple), and the
canonicalized outputs are diffed. HF ↔ HyperIntica agreement is the
hard guarantee; the suite covers polynomial arithmetic through full
multi-variable integrations, and the heavy fixtures additionally ride
SubTropica's benchmark regression at the LibraryLink layer.

## License

MIT — see `LICENSE`.
