# SOFIA.jl — what this is and how it differs from the SOFIA authors' code

This directory vendors **SOFIA.jl**, a Julia translation of the SOFIA package by
Miguel Correia, Mathieu Giroux and Sebastian Mizera (**SOFIA — Singularities of
Feynman Integrals Automatized**, written in the Wolfram Language; MIT); their
original files scr.m and scrj.m and their license are included under `reference/`.
It is not a patch set: none of their code is modified or reused verbatim in the
package source. Every function was translated from their Mathematica source (and
checked against their paper) onto Nemo/FLINT and Graphs.jl, so that the full
pipeline — diagram in, candidate Landau singularities out — runs with no Wolfram
kernel anywhere. That was the point of the translation: SOFIA's default solver
(FastFubini) exists only inside the Wolfram kernel, which makes cluster and CI use
license-bound, while its companion PLD.jl already lives in Julia.

## Attribution

- **SOFIA** by Miguel Correia, Mathieu Giroux, and Sebastian Mizera:
  Comput. Phys. Commun. **320** (2026) 109970,
  [arXiv:2503.16601](https://arxiv.org/abs/2503.16601).
  Repository: <https://github.com/StrangeQuark007/SOFIA>, pinned at commit
  `d4498b95ffcd9f43ba6bce3899c133d3487b5a65` (2025-12-20); see
  `reference/UPSTREAM_COMMIT.txt`.
- **Effortless** (the odd-letter construction) by Antonela Matijašić and
  Julian Miczajka: <https://github.com/antonela-matijasic/Effortless> (MIT).

All algorithmic credit belongs to the original authors. If you use this
software in research, cite the SOFIA paper (and Effortless if you use the
letter construction).

## What is different from the SOFIA authors' code, and why

Architectural translations (same semantics, different mechanics) and deliberate
behavioral divergences are recorded in full, with line references to scr.m and scrj.m, in
`docs/PORTING_NOTES.md`. The substantive ones:

1. **Exact proportionality test** in place of SOFIA's randomized test at
   integer points.
2. **Discriminants as res(f, ∂f)** without the leading-coefficient division;
   the leading coefficient enters the system separately anyway and everything
   is factored and deduplicated downstream, so the singularity set is
   unchanged.
3. **Order-independent subtopology deduplication.**
4. **Effortless made exact**: perfect-square tests via `Nemo.is_square`
   instead of sampling; odd-letter independence by exact linear algebra on
   dlog forms instead of numerical solving at random configurations.
5. **Route-union mode** (`routes=N`): the candidate set on the maximal cut is
   route-dependent in both implementations (the SOFIA authors' example notebook itself
   demonstrates this). SOFIA.jl exposes that honestly and can union candidates
   over deterministic route variants, giving a robust superset of any single
   route.
6. **SymmetryQuotient translated as a performance feature** (subtopology classes
   computed once, results transported through kinematic maps) — provably
   unable to change results, only speed.
7. **Not translated**: the interactive drawing GUI and plotting.

## Validation

382 test assertions (`test/`), plus a replay of the worked examples recorded
in the SOFIA authors' example notebook and an oracle check against the
`PLD_database` of Fevola, Mizera and Telen (96 files, 1794 reference polynomials,
included with SOFIA). Scoreboard, honestly
stated (`validation/REPORT.md`): of the 11 recoverable notebook cases that
completed — 1 exact match, 6 full-coverage supersets, 2 within 1–2 entries;
the massless double box reproduces the database oracle exactly. Six heavy
multi-scale cases exceed the time budget on this engine: the root cause is
isolated to a minimal reproducer (an elimination route reaching a discriminant
of a 6175-term polynomial that stalls both FLINT's multivariate resultant and
subresultant PRS), documented with mitigations and future directions in
`docs/PORTING_NOTES.md`. Every run terminates under the built-in guards
(`time_budget`, `maxopterms`).

## Licensing

- This package: MIT (`LICENSE`), Copyright (c) 2026 Anthropic, PBC for the
  Julia translation; created by Matthew D. Schwartz, code written by Claude
  (Anthropic) under his supervision.
- SOFIA (Correia, Giroux and Mizera): MIT. Their repository carries its MIT grant —
  copyright Correia, Giroux, Mizera, with the full permission text — in its
  README rather than in a standalone LICENSE file (which is why automated
  license detection reports none). That README is vendored verbatim, grant
  intact, at `reference/README.md`; the grant is byte-identical at the
  pinned commit.
- Effortless: MIT, standalone LICENSE vendored intact at
  `reference/EffortlessMarch2025/LICENSE`.

## What is included here

Vendored: the package source (`src/`), tests, porting notes, the validation
harness with its extracted case data, and `reference/` — the SOFIA authors' `scr.m` /
`scrj.m` and README (license carrier), the `PLD_database` oracle,
and the pinned Effortless sources, all under their MIT grants. **Not
vendored**: the 14 MB `SOFIA_examples.nb` worked-examples notebook (size
only); fetch it from the SOFIA repository at the pinned commit if you want
to re-extract the validation corpus with `validation/extract_notebook.py` —
the already-extracted cases ship in `validation/notebook_cases.json`.

## Candidates to offer the SOFIA authors

This is a translation of their package, not a fork; there is no patch to offer the SOFIA authors.
