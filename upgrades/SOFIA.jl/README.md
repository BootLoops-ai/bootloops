# SOFIA.jl

Unofficial Julia translation of the SOFIA package by Miguel Correia, Mathieu Giroux
and Sebastian Mizera (written in the Wolfram Language); not affiliated with or endorsed
by the SOFIA authors — please cite their paper (see "Scientific attribution").

A Julia translation of **SOFIA — Singularities of Feynman Integrals Automatized** —
free of any Wolfram-kernel dependency, built on [Nemo](https://github.com/Nemocas/Nemo.jl)/FLINT
and [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl).

> **Status: functional.** The full pipeline runs end-to-end — diagram in,
> candidate Landau singularities out: FastFubini elimination, loop-by-loop
> Baikov construction with rank-drop handling, subtopologies, and the
> Effortless odd-letter construction. Validated against the recorded outputs
> of the SOFIA authors' example notebook and the `PLD_database` oracle
> (see `validation/REPORT.md` and `docs/PORTING_NOTES.md`).
> The companion [PLD.jl](https://mathrepo.mis.mpg.de/PLD/) solver is
> integrated natively as an *optional* backend (`pld_singularities` —
> SOFIA's `Solver->momentumPLD` option); the caller passes the loaded
> `PLD` module, so PLD.jl and its Oscar stack never become dependencies.

## Scientific attribution

This package is a Julia translation of the Mathematica package
**SOFIA** by Miguel Correia, Mathieu Giroux, and Sebastian Mizera (their files
scr.m and scrj.m, included with their license under `reference/`):

- M. Correia, M. Giroux, S. Mizera, *SOFIA: Singularities of Feynman
  Integrals Automatized*, Comput. Phys. Commun. **320** (2026) 109970,
  [arXiv:2503.16601](https://arxiv.org/abs/2503.16601).
- Their repository: <https://github.com/StrangeQuark007/SOFIA> (MIT;
  the grant is carried in their README — see License below).

The odd-letter construction is a translation of **Effortless** by Antonela Matijašić and
Julian Miczajka (<https://github.com/antonela-matijasic/Effortless>, MIT).

**All algorithmic credit belongs to the original authors.** This repository
contributes only the Julia translation and its test/validation harness.
If you use this software in research, please cite the SOFIA paper above
(and Effortless, if you use the letter construction). The loop-by-loop Baikov
construction follows Frellesvig and Papadopoulos (JHEP 04 (2017) 083,
arXiv:1701.07356) on Baikov's representation (Nucl. Instrum. Meth. A 389 (1997)
347), as in the SOFIA authors' code.

## Why a translation?

SOFIA's default solver (FastFubini) currently exists only inside the Wolfram
kernel, which makes cluster/CI use license-bound. Its companion/back-end
[PLD.jl](https://mathrepo.mis.mpg.de/PLD/) already lives in Julia. This translation
brings the whole pipeline into one open ecosystem.

## Usage

```julia
using SOFIA
import Nemo

# end to end: candidate Landau singularities of the massless double box
mz = BigInt(0)
dbox = Diagram([((1,2),mz),((2,3),mz),((3,6),mz),((6,5),mz),
                ((5,4),mz),((4,1),mz),((2,5),mz)],
               [(1,mz),(3,mz),(4,mz),(6,mz)])
sing, ring = sofia_singularities(dbox)      # -> {s12, s23, s12+s23} = {s, t, s+t}

# with masses: labels are symbols/indexed atoms; same label = same scale
bub = Diagram([((1,2),:m1), ((1,2),:m2)], [(1,mz),(2,mz)])
sing, ring = sofia_singularities(bub)       # contains Källén λ(s, mm1, mm2)

# route-union mode: superset over elimination routes (see PORTING_NOTES)
sing, ring = sofia_singularities(bub; routes=12)

# odd letters (translation of Effortless)
odd = effortless_odd_letters(sing)          # OddLetter(P, Q, ...) entries

# optional PLD.jl backend (SOFIA's Solver->momentumPLD option): if PLD.jl is
# installed, pass the module in — an independent cross-check
# import PLD
# sing, ring = pld_singularities(PLD, bub)

# the building blocks are exported too
fubini, fastfubini, lbl, prepare_landau_system, subtopologies, ...

# the PLD_database (Fevola, Mizera, Telen; included with SOFIA) as a validation oracle
polys, ring, vars = loadsingularities("reference/PLD_database/dbox_zero_zero.m")
```

## Layout

- `src/wlparser.jl` — parser for the Wolfram-InputForm subset used by SOFIA's
  data files; converts to Nemo polynomials (the validation harness).
- `src/polyutils.jl` — proportionality classes, factor lists, resultants /
  discriminants w.r.t. a variable (subresultant PRS).
- `src/fastfubini.jl` — the FastFubini / Fubini elimination engine
  (scr.m:450–513, arXiv:2503.16601 eq. 2.30).
- `src/diagrams.jl` — typed diagram representation, contractions,
  subtopologies (replaces the implicit global state of scr.m).
- `src/pld.jl` — the optional PLD.jl backend bridge (SOFIA's
  `Solver->momentumPLD` option, scrj.m); the `PLD` module is passed in by
  the caller, never depended on.
- `reference/` — the SOFIA authors' pinned sources (scr.m, scrj.m, README) and the
  `PLD_database` (96 files, 1794 reference singularity polynomials) used as the test
  oracle; their commit recorded in `reference/UPSTREAM_COMMIT.txt`. The 14 MB
  worked-examples notebook `SOFIA_examples.nb` is not vendored (size);
  fetch it from the SOFIA repository at the pinned commit — the
  validation cases extracted from it ship in `validation/`. The Effortless
  `Examples/` notebook its vendored README points to is likewise not
  vendored; fetch it from the Effortless repository.

## Deliberate behavioral differences

Documented divergences from the SOFIA authors' code, all validated against the oracle:

1. **Proportionality testing is exact** (coefficient-ratio comparison) where
   SOFIA's `ProportionalPolynomialsQ` is probabilistic at random integer
   points. Strictly fewer false identifications.
2. **Discriminants** are computed as `res(f, ∂f)` without dividing by the
   leading coefficient; scr.m adds the leading coefficient to the system
   separately anyway, and everything is factored and deduplicated downstream,
   so the singularity set is unchanged.
3. **Subtopology dedup** is by canonical edge order (SOFIA's
   `DeleteDuplicates` is order-sensitive); a strictly stronger quotient.
4. **The PLD.jl bridge returns polynomials.** SOFIA's
   `Solver->momentumPLD` option hard-wires a local PLD checkout, runs
   `Pkg.add` at load time, and only echoes the backend's printed output
   (scrj.m:1236-1250); here the caller passes the loaded `PLD` module and
   gets the specialized discriminants back as Nemo polynomials in the
   kinematic ring. The `EulerDiscriminantQ` cross-check is not wrapped.
5. The interactive GUI (`FeynmanDraw`) is not translated.

## Credits and licenses

SOFIA is MIT by its authors (Miguel Correia, Mathieu Giroux and Sebastian
Mizera); their repository carries the grant in its README rather than a
standalone LICENSE file, and that README is included verbatim at
`reference/README.md`. Effortless (Antonela Matijašić and Julian Miczajka) is
MIT with a standalone LICENSE, included intact at
`reference/EffortlessMarch2025/LICENSE`. Everything under `reference/` keeps
its authors' terms; the package `LICENSE` in this directory covers the Julia
translation.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).
