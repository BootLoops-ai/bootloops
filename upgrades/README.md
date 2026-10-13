# upgrades/: engines shipped in full source

Three engines ship here in full source. One was written for BootLoops (Eichler.jl). Two are
ports of other authors' programs: SOFIA.jl, a Julia translation of SOFIA by Miguel Correia,
Mathieu Giroux and Sebastian Mizera, with Effortless by Antonela Matijašić and Julian Miczajka;
and Leviathan, a Julia port of the method of Vsevolod Chestnov, Giulio Crisanti and Mathieu
Giroux, in part translated from their Mathematica code and from DiscKosky by Giulio Crisanti,
Luke Lippstreu, Andrew J. McLeod and Maria Polackova. Each directory has a `PATCHES.md` stating
its provenance: what it ports or follows, from which papers and whose code, what (if anything) is
vendored beside it, and how it was validated. The license of each is in the table; `Eichler.jl`
and `leviathan` also carry a `NOTICE` for the portions that the MIT grant does not cover.

| Directory | What it is | Origin / provenance | License |
|-----------|------------|-----------------------|---------|
| [`Eichler.jl/`](Eichler.jl/) | Written for BootLoops (one block of src/siegel.jl follows FLINT, LGPL-3.0-or-later; see its NOTICE): certified Eichler-integral / sunrise-period layer on Arb ball arithmetic (Γ₁(6) modular forms and iterated Eichler integrals, the all-orders sunrise representation, elliptic-curve periods and quasi-periods, a validated Frobenius/Picard–Fuchs transport engine reaching Calabi–Yau and genus-2 Siegel operators), with its PSLQ dictionaries and audit record; the toolkit members `mirror6` and `theta9` under `../tools/eichler/` build on it | written for BootLoops (conventions from the papers cited inside); six transvectant calls in `src/siegel.jl` follow FLINT `src/acb_theta/g2_covariants.c` (Jean Kieffer); optional GiNaC cross-check oracle as source only in `ginac-oracle/`, GiNaC not bundled | MIT; the FLINT-derived block LGPL-3.0-or-later (Eichler.jl/NOTICE) |
| [`SOFIA.jl/`](SOFIA.jl/) | Julia translation of the SOFIA package by Miguel Correia, Mathieu Giroux and Sebastian Mizera (written in the Wolfram Language): the solver pipeline (FastFubini elimination, loop-by-loop Baikov, subtopologies, the Effortless odd-letter construction of Antonela Matijašić and Julian Miczajka; PLD.jl as an optional backend), so Landau singularity candidates come out with no Wolfram kernel anywhere; GUI/plotting not translated; case-by-case agreement in SOFIA.jl/validation/REPORT.md | [github.com/StrangeQuark007/SOFIA](https://github.com/StrangeQuark007/SOFIA) (Correia, Giroux, Mizera); their original files scr.m and scrj.m with their license, Effortless and the PLD_database oracle included verbatim under `reference/` with their MIT notices | MIT |
| [`leviathan/`](leviathan/) | A Julia port over OSCAR of the Euler-characteristic-drop method for Landau singularities of Chestnov, Crisanti and Giroux (arXiv:2606.29612) and of the elimination object of SPQR (Chestnov and Crisanti, arXiv:2511.14875), in part translated from the authors' Mathematica code (landau_codes.m by Chestnov, Crisanti and Giroux; DiscKosky's euler_characteristic.m by Crisanti, Lippstreu, McLeod and Polackova); the translated portions and the portions that follow their code are listed with file and line in leviathan/PATCHES.md and leviathan/NOTICE | papers cited inside (arXiv:2606.29612, arXiv:2511.14875) and the authors' code on GitHub (Landau-s-Leviathans @ 47028b7, DiscKosky @ 5dca910, SPQR @ dec9b27), which carries no license file | MIT for the parts written for BootLoops; the MIT grant does not cover the translated portions (leviathan/NOTICE) |

The three patched forks of third-party engines that the toolkit drives — Kira 3.1
(IBP reduction), Blade (block-triangular IBP reduction, Wolfram-free) and AMFlow.cpp
(auxiliary-mass-flow numerics, with our validation wrapper suite) — are not part of
this repository. Each lives in its own repository (`kira`,
`blade`, `amflow-cpp`), published beside this one by the
same organization, with the original authors' license carried intact and a `PATCHES.md` describing what changed versus
their code and why. [`ENGINES.md`](ENGINES.md) is the map: repository, original project, base
commit, license, what the patches add, the branch to check out, and the
environment variables through which the packages under [`../tools/`](../tools/README.md)
find the built binaries.

Engines used **unmodified** are not vendored — install them from their authors'
projects; the shopping list with links is
[`../toolkit/external/TOOLS.md`](../toolkit/external/TOOLS.md). The license notes from
that list that bear on building the forks (Fermat, FiniteFlow) are repeated in
[`ENGINES.md`](ENGINES.md).

Each `PATCHES.md` ends with a License section, and the convention is the same
throughout: the code written for BootLoops (`Eichler.jl` except the FLINT-derived block of
`src/siegel.jl`, which is LGPL-3.0-or-later; the Julia translation in `SOFIA.jl` outside
`reference/`, which carries the SOFIA and Effortless authors' MIT notices with ours; and
`leviathan` apart from its translated portions) is MIT, Copyright (c) 2026 Anthropic, PBC,
created by Matthew D. Schwartz with the code written by Claude (Anthropic) under his
supervision; the SOFIA authors' files vendored under `SOFIA.jl/reference/` keep their own MIT
license and copyright notices unchanged; `leviathan` contains no file by the Leviathans, SPQR or
DiscKosky authors, but parts of it are translations of their code, their repositories carry no
license file, and the Anthropic copyright and the MIT grant do not cover the translated portions
(see leviathan/NOTICE and leviathan/PATCHES.md). No GPL-licensed source
ships in this directory (one block of Eichler.jl is LGPL-3.0-or-later; the GPL text beside it is
there because the LGPL supplements it): the GPL-3.0-or-later Kira fork is in its own repository, and
the GPL-licensed software some engines call at run time (OSCAR under Leviathan; GiNaC
behind the optional Eichler.jl oracle, compiled locally from the MIT source in
`Eichler.jl/ginac-oracle/`) is installed from its own project, not bundled. `THIRD_PARTY.md`
at the repository root has the complete table.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).
