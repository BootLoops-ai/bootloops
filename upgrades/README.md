# upgrades/ — engines of our own

The BootLoops harness runs on public engines — but not stock ones. Where an engine
needed a reimplementation (a hard Wolfram dependency, an original with no license file) or
did not exist at all, that work is part of what BootLoops is, and it ships here: three
engines in full source, each MIT-licensed and each with a `PATCHES.md` stating its
provenance — what it reimplements, from which papers or whose code, what (if
anything) is vendored beside it, and how it was validated.

| Directory | What it is | Origin / provenance | License |
|-----------|------------|-----------------------|---------|
| [`Eichler.jl/`](Eichler.jl/) | Original package: certified Eichler-integral / sunrise-period layer on Arb ball arithmetic (Γ₁(6) modular forms and iterated Eichler integrals, the all-orders sunrise representation, elliptic-curve periods and quasi-periods, a validated Frobenius/Picard–Fuchs transport engine reaching Calabi–Yau and genus-2 Siegel operators), with its PSLQ dictionaries and audit record; the toolkit members `mirror6` and `theta9` under `../tools/eichler/` build on it | original (conventions from the papers cited inside); optional GiNaC cross-check oracle as source only in `ginac-oracle/`, GiNaC not bundled | MIT |
| [`SOFIA.jl/`](SOFIA.jl/) | Julia translation of the SOFIA package by Miguel Correia, Mathieu Giroux and Sebastian Mizera (written in the Wolfram Language): the solver pipeline (FastFubini elimination, loop-by-loop Baikov, subtopologies, the Effortless odd-letter construction of Antonela Matijašić and Julian Miczajka; PLD.jl as an optional backend), so Landau singularity candidates come out with no Wolfram kernel anywhere; GUI/plotting not translated; case-by-case agreement in SOFIA.jl/validation/REPORT.md | [github.com/StrangeQuark007/SOFIA](https://github.com/StrangeQuark007/SOFIA) (Correia, Giroux, Mizera); their original files scr.m and scrj.m with their license, Effortless and the PLD_database oracle included verbatim under `reference/` with their MIT notices | MIT |
| [`leviathan/`](leviathan/) | A Julia implementation over OSCAR of the Euler-characteristic-drop method for Landau singularities, from the papers of Chestnov, Crisanti and Giroux (arXiv:2606.29612) and Chestnov and Crisanti (arXiv:2511.14875). No code was copied or translated; the authors' public Mathematica code (landau_codes.m by Chestnov, Crisanti and Giroux; DiscKosky's euler_characteristic.m by Crisanti, Lippstreu, McLeod and Polackova) was read as the reference for five details the papers leave to the code, listed in leviathan/PATCHES.md | papers cited inside (arXiv:2606.29612, arXiv:2511.14875); the authors' code on GitHub carries no license file | MIT |

The three patched forks of third-party engines that the toolkit drives — Kira 3.1
(IBP reduction), Blade (block-triangular IBP reduction, Wolfram-free) and AMFlow.cpp
(auxiliary-mass-flow numerics, with our validation wrapper suite) — are not part of
this repository. Each lives in its own repository (`kira`,
`blade`, `amflow-cpp`), published beside this one by the
same organization, with the original authors' history and license carried intact and a `PATCHES.md` stating exactly what changed versus
their code and why. [`ENGINES.md`](ENGINES.md) is the map: repository, original project, base
commit, license, what the patches add, the branch and tag to check out, and the
environment variables through which the packages under [`../tools/`](../tools/README.md)
find the built binaries.

Engines used **unmodified** are not vendored — install them from their authors'
projects; the shopping list with links is
[`../toolkit/external/TOOLS.md`](../toolkit/external/TOOLS.md). The license notes from
that list that bear on building the forks (Fermat, FiniteFlow) are repeated in
[`ENGINES.md`](ENGINES.md).

Each `PATCHES.md` ends with a License section, and the convention is the same
throughout: the original code (`Eichler.jl`, `SOFIA.jl` outside `reference/`,
`leviathan`) is MIT, Copyright (c) 2026 Anthropic, PBC, created by Matthew D. Schwartz
with the code written by Claude (Anthropic) under his supervision; the SOFIA authors'
files vendored under `SOFIA.jl/reference/` keep their own MIT license and copyright
notices unchanged; `leviathan` contains no code by the Leviathans, SPQR or DiscKosky authors at all. No GPL-licensed source
ships in this directory: the GPL-3.0-or-later Kira fork is in its own repository, and
the GPL-licensed software some engines call at run time (OSCAR under Leviathan; GiNaC
behind the optional Eichler.jl oracle, compiled locally from the MIT source in
`Eichler.jl/ginac-oracle/`) is installed from its own project, not bundled. `THIRD_PARTY.md`
at the repository root has the complete table.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).
