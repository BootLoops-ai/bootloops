# Tools written for BootLoops

These tools were built mid-calculation, when no public tool covered the need as used here.
Each has a page at [bootloops.ai/tools](https://www.bootloops.ai/tools/index.html)
describing what it does, what existed before, and its verification story — that index
organizes them by theme. The package code ships in this repository under
[`../../tools/`](../../tools/README.md), one directory per package with a `GUIDE.md`
each. This file is the roster, grouped as the site index groups it.

The registered packages, grouped as the index groups them:

- **Reduction and IBP** — Seedling (the pre-reduction front door: certified predictive
  staging plus the family preflight), Dogtag (integral-family identity checks on a
  family file before compute is spent on it, with the loomcheck screen as member;
  formerly topology-audit), Winnow
  (upgraded-Laporta elimination with a per-row certificate), Maxcut, Trust (the
  reduction-verification triad), Formglue, and ffcapital (exact salvage of
  FireFly reconstruction state from dead or capped runs).
- **Singularities and geometry** — Landau Alphabet, Leviathan, Dipstick, GeoTriage,
  Coalescer, SubTropica, Surd (exact slice integration with radical last-variable letters
  over SubTropica's fibration conventions), plus the in-house Gröbner backend and the SOFIA
  port behind the adopted singular-locus engines.
- **Differential equations and transport** — Counterweight (the ε-factorizer, with
  Canonify as member), Wayfinder (with the ε→0 extrapolator and the FLINT
  rational-function exporter as members), Famhar, Ratfit (with the degree-alias
  triage as member), Vopclose, Numkin (with basisland, the basis-constrained exact
  lander), PMflow, Cosmoflow, Membound.
- **Periods and arithmetic geometry** — Eichler (with the genus-2 kit as member), Ellipticus
  (the certified eMPL evaluator, with the ellred reduction kit and the gmtel telescoper as
  members), GPLEval, Frobenius, Abacus, Terrier.
- **Numerical oracles and certified arithmetic** — Nestor (with the singularity-subtracted
  dispersion quadrature as member), Longhand (independent numerical ground truth the slow way:
  the arbitrary-precision parametric evaluator, with the pySecDec disteval driver
  `longhand.disteval` as member), amflow-kit (the house kit around the AMFlow.cpp fork:
  the output gates, with the IBP-cache key predictor `keypred` and the launch memory fence
  `memfence` as members; it rides the AMFlow page), Baller (the ball-arithmetic instruments under one
  front door, the mpmath precision linter in its hygiene wing), ERAS, Clinch, Gatekeeper
  (with the kink-splitting quadrature as member), Emitall (the per-result sha-pin
  writer/verifier pairs it standardizes stay with their result pages).
- **Number recognition and closure** — Lockpick (with the curated constant ring as its
  ring-basis member), Annihilator, Ansatzer, Galois, Rankscreen.
- **Exact statistics and evidence** — Mixalot (with Boxwalk, the exact box-moment
  contiguity walker, as member), Popcorn, POSQ, and JaCK & Jill.
- **Operations** — Turnstile (admission control for long jobs on a shared machine:
  priority token, RAM and CPU-width ledgers; formerly Bigram-queue) ships under
  [`../../ops/`](../../ops/README.md), not `tools/`, and has its own page outside the
  thematic index.

Where the code is: most rows above ship as packages under
[`../../tools/`](../../tools/README.md) — the index there gives each package's
verification class. A few roster names differ from their directory names:
Landau Alphabet is `tools/landau-alphabet`, Frobenius is
`tools/frobenius-boundary`, GPLEval is `tools/gpl-eval`, Canonify is
`tools/counterweight/canonical_form`, the genus-2 kit is `tools/eichler/genus2`,
the dispersion quadrature is `tools/nestor/dispersion`, the expansion-by-regions
certifier is `tools/dipstick/regions`, loomcheck ships inside `tools/dogtag`,
basisland ships inside `tools/numkin`, keypred and memfence import as
`amflow_kit.keypred` / `amflow_kit.memfence` from `tools/amflow-kit`,
and the in-house Gröbner backend is a `tools/landau-alphabet` member; Eichler,
Leviathan, and the SOFIA port live under `upgrades/` (`Eichler.jl`,
`leviathan`, `SOFIA.jl`). JaCK & Jill (Python package `phyloexact`) lives in its
own repository, `jackandjill`, published beside this one.
Winnow's certificate layer is the
receipt component vendored inside its package (byte-identical to Trust's receipt
member, `tools/trust/receipt`, test-enforced).

What runs today:

- the packages under [`../../tools/`](../../tools/README.md) — the index there says,
  per package, what runs green from a cold clone and what needs an external engine
  or your data first;
- the per-result reproducers and package bundles on the site
  ([bootloops.ai](https://www.bootloops.ai)), standalone beside their result pages,
  with their (PyPI) dependencies stated in each script's header — the standalone
  witness-verifier bundle and the FORM worked example among them;
- three engines of our own shipped in full under `../../upgrades/`: Eichler.jl (the certified
  Eichler-integral / sunrise-period package), Leviathan (a Julia
  implementation of the Euler-characteristic-drop Landau method of Chestnov, Crisanti and Giroux), and SOFIA.jl (a Julia
  translation of the SOFIA package by Correia, Giroux and Mizera, with no Wolfram dependency) — see [`../../upgrades/README.md`](../../upgrades/README.md);
- the patched forks of Kira, Blade and AMFlow.cpp that the reduction and oracle packages
  drive, each in its own repository with a `PATCHES.md` — listed, with base commits,
  licenses and the environment variables the packages read, in
  [`../../upgrades/ENGINES.md`](../../upgrades/ENGINES.md).

The working patterns the toolkit is used under — recipes kept with their gates and
footguns — are in [RECIPES.md](RECIPES.md). The general protocol layer (acceptance
gates, constant recognition, planted truth, independence bookkeeping, timing
discipline) is packaged as portable skill files in their own repository,
`skills`, published beside this one.

1.0 is installed by cloning (INSTALL.md); there is no pip distribution.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).
