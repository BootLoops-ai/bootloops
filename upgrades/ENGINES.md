# Patched engine forks

Three of the external engines the toolkit drives run as patched forks rather than stock
builds: Kira (IBP reduction), Blade (block-triangular IBP reduction) and AMFlow.cpp
(auxiliary-mass-flow numerics). Each fork is its own repository, started from a snapshot
of upstream at the base commit in the table (upstream's commit history stays in upstream's
repository), with upstream's license unchanged, upstream's notices kept and added to where
the table says so, and a `PATCHES.md` at the repository root that describes the main changes
versus upstream and why, which changes are candidates to offer upstream, and how the fork was
verified (the Kira fork also carries the complete unified diff as `PATCHES.diff`). None of
their source is in this repository. The engines that ship in this repository are a different
matter: Eichler.jl, SOFIA.jl (a Julia translation of SOFIA and of Effortless) and Leviathan (a
Julia port, in part translated from the code of Chestnov, Crisanti and Giroux and of DiscKosky)
are under [`upgrades/`](README.md).

The engine repositories are published beside this one by the same organization; the
commands and defaults below assume they are cloned side by side with this repository
(`../<name>` relative to this repository's root).

| Repository | Fork of | Base | License | What our patches add |
|------------|---------|------|---------|----------------------|
| `kira` | **Kira** by the Kira developers — created by P. Maierhöfer, J. Usovitsch and P. Uwer (Comput. Phys. Commun. 230 (2018) 99, arXiv:1705.05610); developed today by F. Lange, J. Usovitsch and Z. Wu (Kira 3, Comput. Phys. Commun. 322 (2026) 109999, arXiv:2505.20197), with J. Klappert on Kira 2.0 (arXiv:2008.06494) — [gitlab.com/kira-pyred/kira](https://gitlab.com/kira-pyred/kira) | Kira 3.1: tag `3.1` = commit `5760f42` | GPL-3.0-or-later (`COPYING` intact; our delta under the same terms, with a modification notice in each touched file) | parallel `treatcoeff` — each worker thread gets its own Fermat pool slot and the Fermat call runs outside the cache lock, removing the single-Fermat back-substitution bottleneck — plus idempotent SQLite `BEGIN`/`COMMIT` guards; +65/−19 lines over 5 files, reduction output byte-identical to stock. The jemalloc (`-Djemalloc=true`) and 128-bit-weight (`-Dweight_width=128`, binary `kira128`) build variants the harness relies on are upstream meson options, documented there. FireFly is the stock subproject. The same repository carries `bootloops-tools/`: two GPL Python utilities derived from Kira's weight layout and FireFly's save format (`pyred_weight.py`, `ffsave_degree_census.py`) that `kira-stack` and `ffcapital` here use optionally through `tools/kira-stack/kira_gpl_tools.py` (side-by-side checkout, or env `BOOTLOOPS_KIRA_TOOLS`). |
| `blade` | **Blade** by X. Guan, X. Liu, Y.-Q. Ma, W.-H. Wu (Comput. Phys. Commun. 310 (2025) 109538, arXiv:2405.14621), running on FiniteFlow by T. Peraro (JHEP 07 (2019) 031, arXiv:1905.08019) — [gitee.com/multiloop-pku/blade](https://gitee.com/multiloop-pku/blade) | commit `c0a06cb` ("Propagate MaxDegree to recmod") | MIT (`LICENSE.md` and upstream `THIRD_PARTY.md` intact) | Mathematica removed from the compiled toolchain: a no-Wolfram CMake build against a local FiniteFlow; plain-C FiniteFlow wrappers (`fflowC`, transcribed from FiniteFlow's MathLink layer with a checked failure contract) and the `fflowcli` and `dumppoints` executables in place of MathLink; the `fflowcli_stream` sibling with streaming/resumable `evalmany` and subgraph build modes (output bytes identical to stock); a hardened `auto_install`. The optional notebook layer (`Blade.wl`, `BL*/`) is upstream's, unmodified. |
| `amflow-cpp` | **AMFlow.cpp** by the AMFlow.cpp contributors (maintainer GitHub @chang18; doi:10.5281/zenodo.20087172), [github.com/chang18/amflow-cpp](https://github.com/chang18/amflow-cpp), the C++17 port of **AMFlow** by X. Liu and Y.-Q. Ma (AMFlow 2.0 with R.-J. Huang, arXiv:2607.08477; Comput. Phys. Commun. 283 (2023) 108565, arXiv:2201.11669; auxiliary-mass-flow method with C.-Y. Wang, arXiv:1711.09572; [gitlab.com/multiloop-pku/amflow](https://gitlab.com/multiloop-pku/amflow); the Mathematica original is not bundled) | upstream `main` at commit `006a275` (version 1.1.0) | MIT (`LICENSE` unchanged; `NOTICE` is upstream's text with sections appended by the fork; SPDX headers kept). The derivative routines (`src/ibp/libp_deriv.cpp`) descend from LiteIBP by T. Peraro (github.com/peraro/finiteflow-mathtools, which carries no license file) through AMFlow's Kira interface; see the fork's `NOTICE` | precision correctness headed by the exponent-separation solver fix (chop/rationalization tolerances that scale with the ε grid, so `η^(n+aε)` and `η^n` behaviors never merge silently); robustness of the embedded Kira/FireFly reduce (SIG6/SIG11 → FireFly auto-staging with a dot-cap, a CPU-progress stall watchdog, real child-error surfacing, resumable `ff_save` reduces, jemalloc auto-preload); boundary and vacuum extensions (massless-vacuum auto-complete to every loop order, the product-of-tadpoles vacuum for linear propagators, eikonal regions, `Bubble` mode, an env-gated Mercedes maximal-cut closed form); scale and operations (DE checkpoint/resume, ε-grid and node-tree parallelism, a load-adaptive IBP cache, strict JSON config). 579/579 tests pass. `bootloops-wrappers/` in the same repository is our validation harness: the 228-case Mathematica-reference parity run, closed-form and Bessel-moment anchors at 30/100/300 digits, order-by-order comparison and cost-scaling fits. The additions were written with the Mathematica source of AMFlow 2.0 open; `PATCHES.md` says which follow it. |

**What to check out.** In all three repositories the BootLoops changes are on the branch
`main`. Clone and build as each repository's own instructions
say (`README.rst`/meson for Kira, `README.md` + `install.in.txt`/CMake for Blade,
`BUILD.md`/CMake for AMFlow.cpp).

**Two license notes that bear on the builds.** **Fermat** (Robert H. Lewis, Fordham University), the
polynomial-arithmetic backend Kira's back-substitution requires, is freeware distributed
only from its author's site (the source of recent versions is available under the GNU
GPL on request to the author); it is in none of these repositories — obtain it from
<https://home.bway.net/lewis/>. **FiniteFlow**
(T. Peraro, MIT) is the substrate the Blade port links against; it is not bundled with
`blade` either — build it from <https://github.com/peraro/finiteflow> (in-source; see
`tools/blade/GUIDE.md`) and point Blade's CMake at it. The full shopping list of
unmodified engines is [`../toolkit/external/TOOLS.md`](../toolkit/external/TOOLS.md).

## How the packages here find the built engines

Nothing under `../tools/` hard-codes a fork location; the packages read the following
environment variables (or the equivalent command-line flags named in each `GUIDE.md`)
and otherwise fall back to the executable name on `PATH`.

| Engine | Variable | Meaning and default | Read by |
|--------|----------|---------------------|---------|
| AMFlow.cpp | `AMFLOW_CLI` | path to the `amflow_cli` binary built in the `amflow-cpp` checkout (`build/src/cli/amflow_cli`, or on `PATH` after `cmake --install`); default `amflow_cli` on `PATH` | `wayfinder` (sampler wing), `pmflow`, `numkin` (`shift_opt`), `landau-alphabet` (oracle stage) |
| AMFlow.cpp | `AMFLOW_JEMALLOC_WRAP` | optional launch wrapper that preloads jemalloc around `amflow_cli`; unset = run the binary bare | `pmflow`; `amflow-kit`'s smoke gate reads the analogous `AMFLOW_SMOKE_WRAP`, and `AMFLOW_SMOKE_VENDOR_BIN` for the reference binary the first time its triangle reference is generated |
| AMFlow.cpp | `AMFLOW_IBP_CACHE` | shared IBP-cache directory the fork's load-adaptive cache reads and writes (`AMFLOW_IBP_CACHE_DIR` in `numkin`'s injector) | `pmflow`, `wayfinder`, `numkin`, `amflow-kit` (its keypred member predicts the cache keys offline and does not read the variable) |
| AMFlow.cpp | `AMFLOW_CPP_SRC` | path of an `amflow-cpp` checkout (the engine SOURCE tree, read by `pmflow`'s selftest legs that inspect the patched sources; default `<repository root>/../amflow-cpp`) | `pmflow` (selftest legs 6-7 skip by name without it) |
| Kira | `kira` on `PATH`; `KIRA_BIN` | the `kira` (or `kira128`) binary from the `kira` build, found on `PATH` unless a launcher is told otherwise: `kira-stack`'s farm launcher takes `--kira <path>` (default `kira`), `seedling run` takes `--kira-cmd` (default `kira --parallel=4 jobs.yaml`), `numkin_sweep.sh run` takes the binary as its `KIRA_BIN` positional argument, `trust`'s strata legs call `kira` directly; AMFlow.cpp's job JSON names it as `kira_executable` | `kira-stack`, `seedling`, `numkin`, `trust`; `winnow` and `ffcapital` read Kira's and FireFly's files rather than launch them |
| Kira-derived GPL tools | `BOOTLOOPS_KIRA_TOOLS` | directory holding `pyred_weight.py` and `ffsave_degree_census.py` (default `<repository root>/../kira/bootloops-tools`); the loader `tools/kira-stack/kira_gpl_tools.py` raises a named ImportError when they are absent and the battery legs that need them skip by name | `kira-stack` (top-emit, census legs), `ffcapital` (`--census` input) |
| Fermat | `FERMATPATH` | path to the `fer64` binary Kira's back-substitution calls (Kira reads this variable itself; `numkin` and `pmflow` fall back to `fer64` on `PATH`); AMFlow.cpp's job JSON names it as `fermat_executable` and exports it to its Kira child as `FERMATPATH` | every package that launches Kira or `amflow_cli`: `kira-stack`, `seedling`, `numkin`, `trust`, `pmflow`, `amflow-kit`, `wayfinder` |
| Blade | `BLADE_BIN_DIR` | the `bin/` directory of the built `blade` checkout (its CMake install step puts `fflowcli`, `fflowcli_stream`, `dumppoints`, `recmod`, `dynamicrr` there); default `../blade/bin` relative to this repository's root, i.e. a sibling checkout | `blade` (the Python driver layer, `tools/blade`) |
| Blade | `BLADE_PHASE0_DUMPPOINTS` | optional separate location of the `dumppoints` binary; the `BLADE_BIN_DIR` copy is preferred whenever it exists | `blade` |
| Blade | `BLADE_PORT_ROOT` | not an engine path: the read-only bank of sha256-pinned regression artifacts `python3 -m blade.selftest --full` replays, which is not distributed with this repository; unset = the public battery legs run and the replay legs are named and skipped | `blade` |

The other run-time knobs of the AMFlow.cpp fork (`AMFLOW_KIRA_STALL_SECS`,
`AMFLOW_KIRA_FFSAVE`, `AMFLOW_NODE_PARALLEL`, `AMFLOW_DIFFEQ_NUMD`,
`AMFLOW_BOUNDARY_DUMP`, `AMFLOW_KIRA_NO_JEMALLOC`, ...) are documented in that
repository's `PATCHES.md`, `BUILD.md` and `docs/JSON_SCHEMA.md`; the packages that set
them say so in their `GUIDE.md`.

**Licenses of the forks.** The three fork repositories are distributed under their
upstream licenses (GPL-3.0-or-later for `kira`, MIT for `blade`
and `amflow-cpp`; for the derivative routines of `amflow-cpp` that descend from LiteIBP see
the table above),
with our modifications contributed under those same licenses (Copyright (c) 2026
Anthropic, PBC; created by Matthew D. Schwartz, code written by Claude (Anthropic)
under his supervision) and upstream's copyright notices unchanged; each `PATCHES.md`
ends with the full statement. Nothing in `kira` is covered by this repository's MIT
license, and this repository relicenses none of the forks: the packages here run the
built engines as separate programs.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).
