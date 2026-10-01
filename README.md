# BootLoops

BootLoops 1.0 is a harness for large language models doing precision
quantitative science. At its core it is an extensive scientific software package (written,
ported, and upgraded to be driven by an LLM) together with the working
protocols that make the model's results checkable. The harness is independent
of the model driving it: clone it, point whatever agent you use at it, and the
agent gains instruments it can run. Each is documented in the terms an agent
needs: what the instrument does, when to reach for it, what its output means,
and what test its answer must pass before anyone believes it.

The kinds of things it does:

- **Frontier integrals of mathematical physics** — multi-dimensional
  parametric integrals whose answers live in the polylogarithmic, elliptic,
  K3, and Calabi–Yau classes, in closed form or to hundreds of certified
  digits.
- **Recurrences with proofs** — the equation a sum or integral provably
  satisfies, certificate included: strong enough to prove a closed form
  impossible, not merely fail to find one.
- **Bayesian evidence integrals in closed form** — the marginal likelihoods
  that decide model comparison, integrated analytically or with certified
  quadrature instead of sampled.
- **Monte Carlo replaced by convergent methods** — quadrature and recurrences
  that carry certified error bounds instead of statistical ones where the
  problem admits them; where it does not, the sampling tools state
  statistical error bars as such.
- **Ball arithmetic for effective error propagation** — every value carries
  a proven error radius through every step, so rounding and truncation error
  can never go unreported.
- **Exhaustive enumeration with completeness certificates** — check every
  case and prove that none was missed.
- **Open implementations of standard statistical procedures** — routines
  specified by public documentation and the underlying methods papers,
  written from scratch in open code and validated against published outputs.
- **Comprehensive field-specific codebases** — JaCKandJill (Bayesian
  phylogenetics, in its own repository), Mixalot (statistics), Terrier (the
  string landscape), Popcorn (population genetics).

One core runs through all of these: exact reduction, singularity analysis,
certified transport, evaluation to hundreds of digits — machinery that
applies wherever a computation ends in a number someone needs to trust. The
toolkit was built computing scattering amplitudes; that story, and the
results, are told at [bootloops.ai](https://www.bootloops.ai). The package
itself is general purpose.

The harness is meant to grow. Every problem it meets leaves tools behind for
the next one, and additions are welcome from anyone: submit a pull request,
suggest a tool addition or upgrade, or point us at independent code you want
ported or linked — through this repository's issues or the
[contact page](https://www.bootloops.ai/contact.html) at bootloops.ai.

The skills are the other half: the protocols, as plain-markdown instruction
files any agent framework reads, kept in their own repository,
`skills` (see Related repositories below). They
encode what "done" means:
a result reproduces an independent route at points no fit ever saw, with a
positive control proving the check can fail; integer-relation discipline with
the constant ring declared before the search and refusal over invention when
no relation is found; planted-truth controls that recover a known answer
before any real data is touched; provenance bookkeeping so no oracle that fed a fit ever
certifies the result; timing discipline (measure a small run before a big
one); and the heavier protocols for proving with agents, simulating referees,
auditing the literature behind a novelty claim, verifying bibliographies, and
editing prose for precision and honesty (no unsupported claims, no filler,
every quoted number traceable).

## Using it

Clone this repository, start your agent inside it, and put your problem to
the model with the toolkit in front of it:

> Read `tools/README.md` and the tool guides it points to, then propose a plan
> for this problem for me to approve.

Name the object if you can: an integral from a paper, a dataset, a published
number you want checked. The model reads the index, says which instruments
apply (or that none do), and you approve the plan before anything runs. The
tools are cheap to drive: a laptop and whatever model access you already have
are enough.

## Validating the install

One command verifies a fresh clone:

```
python3 run_selftests.py --par 8
```

This walks every package under `tools/` (49 of them), runs each package's own
battery on the toy data it ships with, and prints a status line per package.
The run takes a few minutes on a laptop: the packages come back green, and 3
stop with a named error because they need data the repository does not
include. If a battery needs an external engine you have not installed, it
skips and says what to install. For a first real computation, the one-loop
box through the Landau Alphabet engine:

```
python3 tools/landau-alphabet/test_landau_alphabet.py
```

reproduces the solver's reference results in about a minute and a half.

## Repository map

- [`tools/`](tools/README.md) — the toolkit packages, one directory per package,
  each with a `GUIDE.md` or `README.md` (purpose, acceptance gates, verification class). The
  per-package index is [`tools/README.md`](tools/README.md).
- [`toolkit/`](toolkit/) — the rosters and recipes:
  [`ours/README.md`](toolkit/ours/README.md) (the tool roster by theme),
  [`ours/RECIPES.md`](toolkit/ours/RECIPES.md) (working recipes with their
  checks), [`external/TOOLS.md`](toolkit/external/TOOLS.md) (the external
  engines, with licenses and where to obtain each).
- [`upgrades/`](upgrades/README.md) — the house engines, shipped in full
  source: SOFIA.jl (a Julia translation of the SOFIA package by Correia, Giroux
  and Mizera, with no Wolfram dependency and their original files under
  `reference/`), Eichler.jl, and the Leviathan
  Landau engine. Each directory carries `PATCHES.md` and its license. The
  patched forks of Kira, Blade, and AMFlow.cpp live in their own repositories
  (below).
- [`ops/`](ops/README.md) — operations packages: infrastructure for running
  the toolkit on a shared machine rather than scientific instruments.
  Currently Turnstile (`ops/turnstile`), admission control for long jobs —
  a priority token plus RAM and CPU-width ledgers — with its own self-test.

## Related repositories

BootLoops is six repositories published side by side at
[github.com/BootLoops-ai](https://github.com/BootLoops-ai). This
one, `bootloops`, holds the toolkit and the house engines; the other five
are listed below. Where the code or the guides assume a location for a sibling
repository, it is a checkout beside this one (`../<name>` relative to this
repository's root).

- `jackandjill` —
  JaCKandJill, exact and certified Bayesian evidence for phylogenetic models
  (Python package `phyloexact`, with its own self-certification suite). MIT.
- `amflow-cpp` — the
  patched fork of AMFlow.cpp (precision fixes, Kira/FireFly robustness,
  checkpoint/resume) with the validation wrapper suite in `bootloops-wrappers/`;
  `PATCHES.md` records every change against the AMFlow.cpp authors' code. MIT.
- `kira` — the patched fork of
  Kira 3.1 (parallel back-substitution, SQLite autocommit guards, the
  128-bit-weight build variant); `PATCHES.md` and `PATCHES.diff` record the
  delta, and `bootloops-tools/` holds two GPL Python utilities derived from
  Kira's and FireFly's file layouts that `kira-stack` uses optionally.
  GPL-3.0-or-later.
- `blade` — the patched fork
  of Blade rebuilt to run without a Wolfram kernel, the Blade authors' MIT license intact
  (plain-C FiniteFlow wrappers, a CLI, a sample-point dumper); `PATCHES.md` there.
  The Python port of Blade's search logic that drives those binaries is in
  this repository, under `tools/blade/`.
- `skills` — the working
  protocols as plain-markdown Agent Skills readable by Claude Code, Codex,
  Cursor, Copilot, and most other agent frameworks, with their packaging as
  installable Claude Code and Codex plugins.

## Requirements and tested platforms

The toolkit is plain Python 3.12 (a few packages carry Julia components) with
a small dependency core: `mpmath`, `sympy`, `numpy`, `python-flint`, and
`pytest` for the batteries. Julia packages ship manifests generated on Julia
1.11; Julia 1.12 works after a `Pkg.resolve()`. The release was tested on:

- **Linux x86_64**: Ubuntu 24.04 class (the development platform).
- **Linux containers, x86_64 and arm64**: Debian 12 with Python 3.12 and
  Julia 1.11: the full validation run passes in a fresh minimal
  container on both architectures. One known arm64 limit: the Blade fork's
  from-source installer (`blade`) stops inside MPFR's own test
  suite on arm64; the full Blade build is verified on x86_64.

macOS is not part of the release testing. The pure-Python packages have no
platform-specific code and the Blade fork ships a macOS installer, but run the
validation command above before relying on anything there.

## External engines

Several packages drive external programs (Kira, FireFly, FORM, AMFlow,
FLINT, msolve, Singular, OSCAR among them). The house engines ship in this
repository in full source under `upgrades/`; the patched forks of Kira, Blade,
and AMFlow.cpp build from their own repositories (see Related repositories
above and [INSTALL.md](INSTALL.md)); the rest are **not** downloaded
automatically: you (or your agent) install them when a tool needs one.
Nothing breaks in the meantime: a battery whose engine is
absent skips with a message naming exactly what to install, and
[`toolkit/external/TOOLS.md`](toolkit/external/TOOLS.md) lists every engine
with its license and where to obtain it. A practical route is to
let the model handle it: *"the selftest says FORM is missing; install what
it asks for."*

## Installing

See [INSTALL.md](INSTALL.md): clone the repository for the toolkit and the
house engines; the skills install from the sibling repository
`skills` by clone or through the plugin routes described there.
The toolkit packages are plain
Python (plus Julia for a few); path setup and the external-engine
prerequisites are in the same file.

The papers, the result pages, and the per-problem code live at
[bootloops.ai](https://www.bootloops.ai); this repository is the harness
itself.

## Status, provenance and responsible use

The code was written by Claude working under the author's direction: the
author set every problem, approved every plan, and checked the reported
results against independent routes; each package ships the acceptance battery
described in its GUIDE so users can re-verify.

These are research instruments. Validate outputs before relying on them;
nothing here is intended or fit for clinical, actuarial, payment, regulatory
or public-safety decisions.

## Acknowledgments

BootLoops is a harness around other people's mathematics and software, and it is a pleasure to say whose.

**Exact arithmetic and computer algebra.** Every certified digit here passes through FLINT and Arb (William Hart, Fredrik Johansson, Albin Ahlbäck and the FLINT developers; `acb_theta` by Jean Kieffer) via python-flint (Fredrik Johansson, Oscar Benjamin), Nemo/Hecke and Arblib.jl, alongside mpmath, SymPy, SageMath, `ore_algebra`, Singular, msolve, OSCAR, PARI/GP, GiNaC, fplll and Julia.

**Loop-integral engines.** Kira (Philipp Maierhöfer, Johann Usovitsch, Peter Uwer, Jonas Klappert, Fabian Lange, Zihao Wu) with FireFly (Jonas Klappert, Sven Yannick Klein, Fabian Lange) and Robert H. Lewis's Fermat; the auxiliary-mass-flow method and AMFlow (Xiao Liu, Yan-Qing Ma and collaborators) and AMFlow.cpp (its contributors, maintainer @chang18); Blade (Xin Guan, Xiao Liu, Yan-Qing Ma, Wen-Hao Wu) on Tiziano Peraro's FiniteFlow; FORM (Jos Vermaseren and the FORM developers), HyperFORM (Adam Kardos, Sven-Olaf Moch, Oliver Schnetz) and Erik Panzer's HyperInt; pySecDec (Sophia Borowka, Gudrun Heinrich, Stephen Jones et al.); FIRE and FIESTA (Alexander Smirnov et al.); SOFIA, PLD and SubTropica (Miguel Correia, Claudia Fevola, Mathieu Giroux, Sebastian Mizera, Giulio Salvatori, Simon Telen) with Effortless (Antonela Matijašić, Julian Miczajka); Landau's Leviathans and SPQR (Vsevolod Chestnov, Giulio Crisanti, Mathieu Giroux); PentagonFunctions (Dmitry Chicherin, Vasily Sotnikov, Simone Zoia). The Landau-bootstrap and sequential-discontinuity tools grew out of work with Holmfridur Hannesdottir, Andrew McLeod, Cristian Vergu and Jacob Bourjaily; the lattice-reduction regression out of work with Oscar Barrera, Aurélien Dersy, Rabia Husain and Xiaoyuan Zhang.

**Applied packages and data.** The population-genetics, ecology, seismology, phylogenetics and public-records tools stand on dadi, fitdadi, polyDFE, fastDFE, moments, msprime/tskit, `etas`, matPTF, SageMath's genus-2 modules and the papers named in each `GUIDE.md`; the datasets of Dmitry Chicherin, Ian Moult, Emery Sokatchev, Kai Yan and Yunyue Zhu and of the DravLex team (Vishnupriya Kolipakam and co-authors) are used under CC BY 4.0 with thanks.

Full author lists and references are in [`toolkit/external/TOOLS.md`](toolkit/external/TOOLS.md), [`REFERENCES.md`](REFERENCES.md) and each package's CREDIT paragraph; please cite those authors, not only BootLoops. Omissions are ours to fix: open an issue and they will be.

## Maintenance, reporting and security

**Maintenance.** This repository is maintained by Matthew D. Schwartz, not by
Anthropic. It is not an officially supported Anthropic product, and Anthropic
does not provide support, updates or fixes for it.

**Reporting issues.** Please report bugs and security problems through this
repository's GitHub issues.

**Security considerations.** Treat input files from others as code. These are
research tools meant to be run locally on inputs you trust. Many of them
evaluate the contents of their input files (JSON, YAML, `.m`, `.ms`, `.jl`,
pickle and similar), so a file received from someone else can run arbitrary
commands on your machine. Only run files you wrote yourself or got from a
source you trust, or run them in a sandbox or container. The integrity checks
and certificates in this repository guard against accidents. They are not a
security boundary.

## License and attribution

BootLoops 1.0 is released under the MIT License ([LICENSE](LICENSE)),
Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written
by Claude (Anthropic) under his supervision. This is
not an officially supported Anthropic product; it is maintained by Matthew D. Schwartz
(https://www.bootloops.ai). The prose and figures written for this repository
are released under CC BY 4.0 ([LICENSE-CONTENT](LICENSE-CONTENT)). The reference copies of third-party
sources under `upgrades/SOFIA.jl/reference/` and `tools/subtropica/reference/`, and the
Python port of Blade's search logic in `tools/blade/` and the Julia port of
SubTropica's front end in `tools/subtropica/`, keep their original
licenses and copyright notices (MIT under their authors' copyright); three files are
GPL: `tools/subtropica/src/lr_refine.jl` (with its test), derived from Erik Panzer's
HyperInt, GPL-3.0-or-later (see `tools/subtropica/NOTICE`);
`tools/eichler/genus2/mestre_port.py`, a sympy translation of SageMath's `mestre.py`
and `invariants.py` by Florian Bouyer, Marco Streng and Nick Alexander,
GPL-2.0-or-later (see `tools/eichler/NOTICE`); and `tools/formglue/form_hyper.py`,
whose HyperFORM driver template follows the example drivers shipped with HyperFORM
by Adam Kardos, Sven-Olaf Moch and Oliver Schnetz, GPL-3.0-only (see
`tools/formglue/NOTICE`). The
patched engine forks are distributed in their own repositories under the same
organization, each under the license of its original authors (Kira GPL-3.0-or-later; Blade
and AMFlow.cpp MIT), and are not part of this repository.
[THIRD_PARTY.md](THIRD_PARTY.md) has the full list and [NOTICE](NOTICE) the
attribution.

To cite BootLoops: M. D. Schwartz, *BootLoops 1.0* (2026),
[bootloops.ai](https://www.bootloops.ai).
