# Changelog

All notable changes to BootLoops are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/) (1.0.x = fixes and tool additions that change nothing existing, 1.x.0 = changes to existing tools).

## Unreleased

### Changed
- License of six files and one block. `tools/winnow/ibplapper/adapters/user_system.py`, `tools/kira-stack/parallel_kira_gen.py`, `tools/kira-stack/parallel_kira_gen_uds.py`, `tools/ffcapital/harvest_ffsave_eta.py` and `tools/ffcapital/harvest_ffsave_2var.py` move from MIT to GPL-3.0-or-later; `tools/eichler/genus2/conic_fast.py` moves to GPL-2.0-or-later; one block of `upgrades/Eichler.jl/src/siegel.jl` (the transvectant calls of `igusa_clebsch`) is LGPL-3.0-or-later. Because the package `ibplapper` imports `user_system.py`, `ibplapper` as imported is a combined work under GPL-3.0-or-later (see `tools/winnow/NOTICE`, which also says how to make an MIT-only build).

### Fixed
Credits, notices and license headers for code and data derived from other authors' work. The edits are
comments, headers, documentation and notice files, with the hash tables of the files whose headers changed
regenerated. No computed result changes; the only changes to program output are printed messages and help
texts: two messages in `counterweight`, one in `pmflow detect`, three error messages of
`wayfinder/boundary.py`, the first line of the usage text of `leviathan.jl`, and the help texts of `pmflow`,
of `seedling/audit.py` and of `emit_receipt.py` in `frobenius-boundary`. The license changes are listed under
Changed above.

- `leviathan`: the package is a Julia port of the Euler-characteristic-drop method of Chestnov, Crisanti and Giroux (arXiv:2606.29612), in part translated from their Mathematica code and from DiscKosky by Crisanti, Lippstreu, McLeod and Polackova; the statement that no code was copied or translated is withdrawn. README, `PATCHES.md` and source headers corrected, equation numbers corrected, a package `NOTICE` added.
- `dipstick`, `maxcut`: both name the Leviathan engine they call and its authors; `dipstick` credits Pak and Smirnov for the region enumeration.
- `SOFIA.jl`: the header of `src/effortless.jl` cites the Effortless thesis and states the differences from Effortless.m.
- `Eichler.jl`: one block of `src/siegel.jl` follows FLINT's `g2_covariants.c` (Jean Kieffer) and is marked LGPL-3.0-or-later, with `NOTICE` and the license text added.
- `eichler`: `genus2/conic_fast.py` carries six statements of SageMath's `mestre.py` and now has the GPL-2.0-or-later header that `mestre_port.py` has.
- `baller`: `vendor/certlane` is a port, into ball arithmetic, of the alert layer of INGV's matPTF (CC BY 4.0); `vendor/geo/mc.py` follows pyCSEP's Poisson L-test simulation (BSD-3-Clause); `pvand` is the Björck-Pereyra solver as coded by John Burkardt (MIT). File headers and `tools/baller/NOTICE` added.
- `wayfinder`: `boundary.py` states that its vacuum table, SingleMass prefactor and ending-value rule are a Python port of AMFlow.cpp routines (after AMFlow by Liu and Ma); the package ships their MIT notices and the mpmath notice in `tools/wayfinder/NOTICE`.
- `pmflow`: cites Z.-F. Liu and Y.-Q. Ma, arXiv:2201.11636, and says its fixed-point observation was made on our AMFlow.cpp build.
- `winnow`: the reader `ibplapper/adapters/user_system.py` follows Kira's parser and is now distributed under GPL-3.0-or-later with the Kira Developers' notice (`NOTICE`, `LICENSE-GPL-3.0`).
- `kira-stack`: `parallel_kira_gen.py` and `parallel_kira_gen_uds.py` follow Kira's source and are now distributed under GPL-3.0-or-later with the Kira Developers' notice (`NOTICE`, `LICENSE-GPL-3.0`).
- `ffcapital`: `harvest_ffsave_eta.py` and `harvest_ffsave_2var.py` carry the first 40 entries of FireFly's prime table and follow its state-file writer; they are now distributed under GPL-3.0-or-later with the FireFly authors' notice (`NOTICE`, `LICENSE-GPL-3.0`).
- `seedling`, `ratfit`: credits and wording corrected.
- `terrier`: the genus traversal in `lattice/genus_enum.jl` is credited as adapted from Hecke.jl (BSD-2-Clause), with the notice in `tools/terrier/NOTICE`; the license statement for the arXiv:2107.09064 data corrected.
- `popcorn`: `twowindow_sim.py` credits twosfs (Rice, Fenton, Novembre, Desai; MIT) for its seed-hash function and branch-spectrum call, and Fenton et al. for the two-site test; `tools/popcorn/NOTICE` added.
- `clinch`: the certpass adapter says which lines are copied from, and which follow, the etas package of Mizrahi, Nandan and Wiemer, and carries their copyright line.
- `posq`: credits OPQ (Gautschi; Laurie) for the Jacobi recurrence coefficients.
- `nestor`: `ladder.py` names mpmath as the source of its depth loop, and the package ships the mpmath notice in `tools/nestor/NOTICE`.
- `surd`: the cell data are stated to derive from the file P_qbp_qp_g_q.m on the authors' web page (Del Duca, Duhr, Haindl, Lazopoulos, Michel), which states no license; the CC BY 4.0 label is withdrawn.
- `landau-alphabet`: credits Fevola, Mizera and Telen and Mizera and Telen for the method of the numerical fallback; a statement about PLD.jl's output in one graph specification corrected.
- `blade`: lists all nine modules derived from Blade and describes `probe.py` as a counterpart of Blade's `ssolve`.
- `subtropica`: names the SubTropica authors in `Project.toml` and completes its `NOTICE`.
- `counterweight`: credits Lee and the public implementations of his algorithm, and states that its module `Fuchsia` is not the Fuchsia program of Gituliar and Magerya; two printed messages of `canonical_form.py` that named the wrong program and repository are corrected.
- `cosmoflow`: states that it is unrelated to the CosmoFlow package of Werth, Pinol and Renaux-Petel and lists its references.
- `annihilator`, `ansatzer`, `ellipticus`, `eras`, `frobenius-boundary`, `galois`, `gpl-eval`, `mixalot`, `numkin`, `tropical-sampler`, `trust`: credit paragraphs name the earlier programs and papers for each method, and claims of novelty are narrowed to what each package adds.
- Repository files (`THIRD_PARTY.md`, `NOTICE`, `REFERENCES.md`, `README.md`, the rosters under `toolkit/`, `tools/README.md`, `upgrades/README.md`, `upgrades/ENGINES.md`): rows and notices added for the items above; the ports under `upgrades/` are no longer described as engines of our own; the fork repositories are described as snapshots on the branch `main`, without upstream's commit history.
- Sibling repositories: `jackandjill` ships the PSF notice for its big-integer printer (after CPython's int_to_decimal by Tim Peters) and credits Sebastiano Vigna's generators; `amflow-cpp` names T. Peraro's LiteIBP as the origin of the derivative routines and lists the additions that `PATCHES.md` had left out; `blade` corrects two statements in `PATCHES.md` and adds a note on git authorship.

## 1.0.1 (2026-10-09)

### Added
- `cutbasis`: geometric master-basis selector on the maximal cut (arXiv:2511.15381
  Algorithms 1-3) writing Kira `preferred_masters` lists in three forms, with a certified
  reducer comparison (`evaluate`) and a self-test; worked example: the two-loop Møller
  double box with a massive loop.

### Fixed
- `run_selftests.py`: a crash in a data-gated package was reported as
  "REFUSED (by design)" and left the run green. A data-gated battery now
  signals a refusal only by exit code 4 (`DATA_GATED_REFUSAL_RC`) with a
  printed reason, kept in `selftest_results.json`; any other nonzero exit,
  an exit 4 without a reason, or one whose last line is an uncaught
  exception, is a FAIL. `ffcapital`, `frobenius-boundary` (was exit 3) and
  `galois` follow the contract; `galois` resolves its reference data before
  any check runs. Documented in `tools/README.md` (Verification classes)
  and `INSTALL.md`. Reported by @Aspect5 (BootLoops-ai/bootloops#3).
- `run_selftests.py`: a timeout in the Julia `Pkg.instantiate()` pre-step
  aborted the whole run and `selftest_results.json` was never written. The
  timeout now fails that one package (with the cause and the remedy in the
  record) and the run continues, exiting 1. Fix by @hzshen88
  (BootLoops-ai/bootloops#2), reported in BootLoops-ai/feedback#1.
- `nestor`: the farm's orphan guard compared each worker's parent pid with
  the pool parent's pid, which is wrong under the `forkserver` start method
  (Python 3.14's default on Linux): every worker exited at its first chunk
  and the farm self-test (P5) never finished. The guard now uses the parent
  pid seen from inside the worker and also checks that the pool parent is
  alive, so it behaves the same under `fork`, `forkserver` and `spawn`.
  Reported by @drdavidglew-afk (BootLoops-ai/feedback#2).
