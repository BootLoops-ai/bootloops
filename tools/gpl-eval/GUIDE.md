# gpl-eval (GPLEval.jl) — guide

Tool page: https://bootloops.ai/tools/gpl-eval.html

KIND: julia package (`GPLEval.jl/`: src/{GPLEval,gpl_core,hpl,symbol_eval,
ginac_bridge,empl_bridge}.jl + src/ginac_gpl.cpp, compiled on first use; CLI
bin/gpl_eval.jl)

PURPOSE: Evaluator for Goncharov polylogarithms (GPLs) and harmonic polylogarithms
(HPLs) with an exact iterated-integration core that matches the
fit's basis conventions — evaluates symbol-word bases at value-fit sample points at
arbitrary precision (Arb balls).

USE-WHEN:
- HPL/GPL evaluation at fit sample points (`gpl_eval` / `hpl_series`).
- Building basis-value columns for lockpick fits or gatekeeper heldout CV.

NOT-FOR: elliptic-MPL / Kronecker eMPL — separate published software (pointers only;
the empl_bridge connects toward those evaluators).

DEPENDENCIES: Project.toml + Manifest.toml ship (Manifest generated under
julia 1.11; the in-repo packages are recorded by relative path), so a fresh
clone instantiates cold: `julia --project=GPLEval.jl -e 'import Pkg;
Pkg.instantiate()'` — network needed once for the registered deps Arblib +
JSON3. GPLEval hard-depends on two Julia packages that BOTH ship in
this repo:
- `Eichler` — `upgrades/Eichler.jl`.
- `LandauAlphabet` — `tools/landau-alphabet/LandauAlphabet.jl`, the Julia
  engine of the Landau Alphabet package (`ItIntEvaluator`/`GPLContext`, the
  fast suffix-cached series path).
On a different julia version, or with the deps relocated, re-wire instead:
`julia GPLEval.jl/deps/setup.jl` develops both in-repo packages
into the project (env-overridable: `LANDAUALPHABET_JL_PATH`, `EICHLER_JL_PATH`)
and instantiates the registered deps Arblib + JSON3. The ginac_bridge
additionally needs g++ and the GiNaC/CLN C++ libraries installed (external):
it compiles src/ginac_gpl.cpp into src/ginac_gpl on first use. No binary is
shipped; GiNaC and CLN are GPL-licensed, so the locally built helper is a GPL
combined work that GPLEval only runs as a subprocess. To build it by hand
(e.g. for other packages that call the same CLI):
`g++ -O2 -std=c++17 -o GPLEval.jl/src/ginac_gpl GPLEval.jl/src/ginac_gpl.cpp $(pkg-config --cflags --libs ginac)`.

INVOKE: `julia --project=GPLEval.jl GPLEval.jl/bin/gpl_eval.jl --alphabet alpha.json
--basis basis.json --point '{"s":-3,"t":-7,"z":0.5}' --prec 256 --out values.json`

INPUTS:
- alpha.json: `{"letters":{"l0":"0","l1":"1","ls":"s+t",...},"zvar":"z"}` — each
  letter is a Julia expression in the kinematic variables (keys of --point) and is the
  GPL INDEX aᵢ (the root of the dlog letter), NOT the dlog argument; "zvar" names the
  integration variable in --point.
- basis.json: list of symbol words `[["l0","l1"],["l0","l1","lm1"],...]`.
- --prec in BITS.

OUTPUTS: values.json `{"prec_bits":..., "values":[{"re":"...","im":"...",
"rad":"..."},...]}` — Arb balls; use the rad field, it is the certificate.

GATES: Arb ball radii bound the result; cross-check against the ginac bridge path
when a value is load-bearing.

FOOTGUNS:
- Letters are roots aᵢ; feeding dlog arguments silently evaluates the wrong
  functions.
- Precision is bits, not dps.

BATTERY: `python3 selftest.py` (canonical entry — wires the environment if
this checkout is not yet instantiated, then runs the suite; a missing `julia`,
or an environment it cannot wire, is a loud NAMED skip stating exactly what to
install — never a raw load error; class: partial).
Equivalently, `julia --project=GPLEval.jl GPLEval.jl/test/runtests.jl`
after the wiring (see DEPENDENCIES) — the CLI and the full test
suite run from this checkout. The ginac cross-check legs need GiNaC (external)
and skip with a warning when its build is unavailable; the empl leg needs the
separate eMPL evaluators (external published software).

RELATED: consumes alphabets from landau-alphabet/ansatzer pipelines; feeds
lockpick value-fits and gatekeeper CV. ginac_bridge/ginac_gpl.cpp = GiNaC
cross-check leg; empl_bridge connects toward the elliptic evaluators
(upgrades/Eichler.jl).

CREDIT: Goncharov multiple polylogarithms [Gon] and Remiddi–Vermaseren HPLs [RV];
reference numerics through GiNaC by Christian Bauer, Alexander Frink and Richard Kreckel
[GiNaC] with the polylogarithm algorithms of Jens Vollinga and Stefan Weinzierl [VW]
(GPL, compiled locally, run as a subprocess); on-path test words from Ajjath, Chaubey
and Shao [Ajjath]; ball arithmetic by Arb [Arb]. Elliptic MPLs route to tools/ellipticus /
Eichler.jl (Broedel–Duhr–Dulat–Penante–Tancredi [BDDT, BDDPT]; Walden–Weinzierl [WW]).
Bracketed keys resolve in REFERENCES.md at the repository root.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).
