# Installing BootLoops

BootLoops has two parts: the **toolkit** (the computing packages and the
house engines, in this repository) and the **skills** (the working protocols,
plain-markdown instruction files in the
[Agent Skills](https://agentskills.io) format, readable by Claude Code, Codex,
Copilot, Cursor, and most other agent frameworks, in their own repository,
`skills`, published beside this one by the same organization). The
patched engine forks the toolkit drives (Kira, Blade, AMFlow.cpp) build from
their own sibling repositories (`kira`, `blade`,
`amflow-cpp`). Nothing activates without
your say-so on any route.

## Clone the repository (the toolkit)

```
git clone https://github.com/BootLoops-ai/bootloops
```

No GitHub account is needed: public repositories clone anonymously, and the
same tree downloads as a plain tarball from the repository page.

You get the full toolkit under `tools/` and the house engines under
`upgrades/` (SOFIA.jl, Eichler.jl, Leviathan). (The papers and the per-problem
result files live on [bootloops.ai](https://bootloops.ai), not in the
repository.)

## The toolkit packages (`tools/`)

The packages are plain Python (a few carry Julia components); there is no
installer — a clone is the installation. Per-package setup:

- **Python path.** Each package imports from its own directory; siblings find
  each other as directories under `tools/`. For scripting, put the repo's
  `tools/` on `PYTHONPATH` (or `sys.path`) — the flat shim files at the `tools/`
  root provide flat import names for the packages. Every package's
  `GUIDE.md` names its entry points.
- **Common Python dependencies** (PyPI): `mpmath`, `sympy`, `numpy`,
  `python-flint` — `pip install mpmath sympy numpy python-flint`. A few
  packages need one extra beyond this core (`scipy`, `pyyaml`, `gmpy2`,
  `networkx`, `dynesty` among them); each names its own in the Requirements
  line of its GUIDE, and all are pip-installable. `popcorn`'s simulation
  members (`popcorn.planted` regeneration, `popcorn.twowindow`) need
  `msprime` + `tskit`, and its transient cross-check optionally uses
  `moments` (PyPI `moments-popgen`); all three are optional and the battery
  legs that need them skip by name without them. When one is missing, the
  selftest runner (below) names the missing module in its per-package log.
- **External engines.** Several packages drive external programs — Kira,
  FireFly, FORM, AMFlow, Blade, FLINT, msolve, Singular, OSCAR/Julia among them.
  The catalog, with licenses and where to obtain each, is
  [`toolkit/external/TOOLS.md`](toolkit/external/TOOLS.md). The house engines
  ship in full source under [`upgrades/`](upgrades/README.md) with
  per-directory build notes; the patched forks of Kira, Blade, and AMFlow.cpp
  build from their own repositories (next section); packages locate the
  binaries via `PATH` or the environment variables named in their GUIDEs.
- **Verification.** Most packages ship a selftest or battery; the verification
  class of each (what runs green from a cold clone, what needs an engine or your
  data first) is in the [`tools/README.md`](tools/README.md) index. Legs that
  need reference data the repository does not ship skip or fail closed with a
  named error — that is by design, not breakage.

## Building the engines from their repositories

Three engines the toolkit drives are patched forks, each distributed in its own
repository beside this one (same organization) with the upstream license intact
and a `PATCHES.md` stating what changed and why. Clone the ones your problem
needs side by side with this repository, so that each sits at `../<name>`
relative to this repository's root (the one default path the code assumes, see
below); a battery whose engine is absent skips with a message naming it.

- **Kira 3.1 fork** (GPL-3.0-or-later): the sibling repository
  `kira` (clone it beside this one: `../kira`) —
  see its `PATCHES.md`
  (and `PATCHES.diff`) and the build notes there; Kira's build fetches FireFly,
  and Fermat comes from its own site (proprietary freeware). Its `bootloops-tools/`
  directory holds two GPL Python utilities (Kira weight codec, `ff_save` degree
  census) that `kira-stack` and `ffcapital` pick up from that side-by-side checkout,
  or from the directory in `BOOTLOOPS_KIRA_TOOLS`; their battery legs skip by name
  without them.
- **Blade fork, Mathematica-free toolchain** (MIT): the sibling repository
  `blade` (`../blade`) — see its `PATCHES.md`
  and the `auto_install` / `auto_install_macos` scripts; it links FiniteFlow,
  obtained upstream.
- **AMFlow.cpp fork** (MIT): the sibling repository
  `amflow-cpp` (`../amflow-cpp`) — see its
  `PATCHES.md` and `BUILD.md`; the validation wrapper suite is in
  `bootloops-wrappers/` there.

How the packages find the built binaries (only what the code reads):

- `kira` and `fer64` on `PATH`. `FERMATPATH` points Kira at the Fermat binary
  and is read (defaulting to `fer64`) by `tools/pmflow`, `tools/kira-stack`,
  and the `tools/numkin` launch scripts; `tools/kira-stack/parallel_kira_gen.py`
  takes `--kira PATH` and `tools/numkin/numkin_sweep.sh` a positional Kira path
  when the binary is elsewhere.
- `AMFLOW_CLI` — the `amflow_cli` binary (read by `tools/pmflow`,
  `tools/wayfinder/sampler.py`, `tools/numkin/shift_opt.py`, and
  `tools/landau-alphabet/LandauAlphabet.jl`; unset means `amflow_cli` on
  `PATH`). `tools/pmflow` also reads `AMFLOW_JEMALLOC_WRAP` (an optional
  wrapper script; unset runs the CLI bare) and `AMFLOW_IBP_CACHE` (its IBP
  cache directory).
- `AMFLOW_CPP_SRC` — an `amflow-cpp` checkout (source tree), read only by
  `pmflow`'s selftest legs that inspect the patched sources (default
  `../amflow-cpp` beside this repository; the legs skip by name without it).
- `BLADE_BIN_DIR` — the directory holding the built Blade binaries (`fflowcli`,
  `dumppoints`, ...), read by `tools/blade/pipeline.py`;
  `BLADE_PHASE0_DUMPPOINTS` names an alternate `dumppoints` binary.
  (`BLADE_PORT_ROOT`, read by `tools/blade/selftest.py`, is a reference-data
  root for the battery's data-gated legs, not a binary location.)

## Verifying an installation

Every package's battery has one canonical entry, pinned in
[`tools/BATTERIES.json`](tools/BATTERIES.json) — command and working directory.
`run_selftests.py` runs exactly those; nothing is guessed. To verify a clone:

```
cd bootloops
python3 -m venv .venv && source .venv/bin/activate     # macOS/Linux
pip install mpmath sympy numpy python-flint pytest     # the core
pip install scipy pyyaml gmpy2 networkx dynesty cypari2 msprime tskit   # per-package extras
python3 run_selftests.py --par 8
```

The venv activation matters: batteries invoke `python3` from your PATH, and a
bare system python (macOS ships 3.9 with nothing installed) fails every
import. For the Julia-component packages install Julia ≥1.11 (macOS:
`brew install julia` or juliaup; Linux: juliaup or your distribution's
package; on managed machines whose endpoint security blocks juliaup's
downloaded binaries, prefer Homebrew or the distribution package); the
runner instantiates each package's committed Julia environment
automatically. `abacus` additionally wants the
PARI/GP binary `gp` on PATH.

How to read the table:
- **PASS** — the battery ran green.
- **REFUSED (by design)** — data-gated packages refuse loudly without your
  data or a reference-data root; that is their contract, not breakage.
- **FAIL** — a real problem; the exit code is then 1 (a REFUSED row never sets
  it). Per-package commands, exit codes, and output tails land in
  `selftest_results.json`.

The runner walks `tools/` only. The operations package under `ops/`
(Turnstile, admission control for long jobs on a shared Linux machine) has its
own self-test, `bash ops/turnstile/selftest.sh` from the repository root; it
reports a named SKIP on macOS.

Name packages to check just those (`python3 run_selftests.py mixalot eras`);
`--timeout` raises the per-battery limit (300 s by default) on slow machines. A
single package's battery can also be run by hand: look up its line in
`tools/BATTERIES.json` and run that command from the stated directory (`cwd`
"package" means the package's own directory under `tools/`, "root" the
repository root; the runner puts `tools/`, the package directory and the root on
`PYTHONPATH`, so do the same by hand).

## The skills

The working protocols live in their own repository, `skills`
(published beside this one by the same organization), which also carries
their packaging as Claude Code and Codex plugins and the install
routes for other agents (Cursor, Copilot, the skills.sh installer, the web
apps). Two routes cover most setups: clone that repository, open Claude Code
in the clone, and run `/bootloops-setup` to pick which skills to activate and
at what scope; or, without cloning, `/plugin marketplace add
BootLoops-ai/skills` in Claude Code and install the plugins it lists. The
skills arrive inert on every route; in any other agent, point it at a skill's
`SKILL.md` directly ("read the acceptance-gate skill and hold my result to
it") or copy chosen skill folders into the agent's skills directory. That
repository's README has the details.
