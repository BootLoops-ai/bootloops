# blade — the Wolfram-free Blade block-triangular IBP pipeline

Tool page: https://bootloops.ai/tools/blade.html · pairs with the blade fork,
which lives in the sibling repository `blade` (clone it beside
this one: `../blade`; FiniteFlow comes from upstream).

KIND: package + fork (Wolfram-free C pipeline redg1→fitrel→dumppoints→ssolve→
recmod→dynamicrr). Full module docs: `README.md` in this directory.

PURPOSE: Blade block-triangular (BT) IBP port. One-time BT search per family, then:
per-point reduction probes at ~0.1–0.3 ms/pt (BTOracle), full symbolic reduction
tables (semibl), and standalone multi-prime multivariate rational reconstruction
over cached evaluations (ratrec). BT form collapses per-probe reduction cost
100–400× (the Blade paper); measured here ~10⁴–10⁵× vs per-probe Kira.

USE-WHEN:

- MANY cheap reduction probes at one kinematic point (η/ε grids, reconstruction
  sampling) where per-probe Kira/Fermat is the wall. Measured (db family): 0.286
  ms/pt python exact backend, 0.095 ms/pt ssolve-batched, vs 19.7 s real numkin Kira
  and 86 ms/pt fflow-native; one-time BT search ≈27 s ⇒ break-even at the 2nd point.
- Full symbolic reduction table from a closed BT form: `blade.semibl` (gated:
  dbox table == independent Kira 1536/1536 mod-p + 128/128 exact-ℚ).
- Standalone ratrec over a CACHED evaluation table: `blade.ratrec` — complement
  to FireFly for numeric harvests over cached tables.
- New family: START with the turnkey CLI `blade_family.py prepare` (one command
  spec→stored family directory).
- No Kira at all: give `prepare` a raw sparse IBP(-like) system in the fflow
  JSON grammar (`fflow_system` + `fflow_columns` spec keys) — the database is
  built by per-point fflow sparse elimination (`fflowcli learn`/`evalmany`),
  and the held-out relation gate replays the system's own solve at a fresh
  prime.

NOT-FOR: it does not launch Kira (`kira_config_dir` only DISCOVERS an existing
table); numeric probes are cross-checks at production scale — an exact-row claim
still needs semibl + a verified basis rotation.

INVOKE:

- Turnkey: `python3 tools/blade/blade_family.py prepare <spec.json> --workdir W
  [--max-round-seconds 1800] [--bank-dir DIR] [--bisect] [--nthreads 6]`, then
  `probe <family_dir> --points pts.json [--prime P] [--batch]` and
  `status <family_dir>`. Exit codes: 0 all blocks closed; 2 partial closure stored
  honestly; 1/3 typed failure (ERROR.json names the artifact). FRESH workdir
  required (stale-state refusal — deliberate).
- API: `from blade.probe import BTOracle; o = BTOracle(family_dir);
  o.reduce_point(kin, eps, prime)` / `o.reduce_many(...)` /
  `o.reduce_many_ssolve(...)`.
- Symbolic table: `blade.semibl`. Ratrec: `blade.ratrec`
  (`dump_points` / `scan_degrees` / `reconstruct` / `ratfun_evaluator`).

Binaries: build the blade fork from `../blade` (and
FiniteFlow — build it IN-SOURCE; out-of-source cmake misses the generated
`config.hh`), then set `BLADE_BIN_DIR` to the bin dir (default: `../blade/bin`
relative to this repository's root, i.e. the sibling checkout).
`BLADE_PHASE0_DUMPPOINTS` optionally names a separate dumppoints.

## What runs from a checkout, and what needs the reference bank

- The library (formats/pipeline/probe/ratrec/semibl/search/blade_family) is the
  shipped product; it imports clean with plain python3 + sympy (lazy) + python-flint
  and drives the binaries above. Battery (import leg, from a scratch cwd):
  package import OK, relative BLADE_BIN_DIR default resolves, fflow prime table
  regenerates.
- The regression battery `python3 -m blade.selftest --full` replays every
  phase's decisive gate from a bank of 242 sha256-pinned artifacts
  (`regression_manifest.json`) that is not distributed with this repository
  (point `BLADE_PORT_ROOT` at it) — missing/changed artifacts
  are NAMED loudly and the replay aborts, no silent skips (verified live: an
  incomplete bank yields named MISSING failures and 0/244, by design). Public
  clones get the SELF-CONTAINED raw-system section of `python3 -m
  blade.selftest` (typed lane refusals; then, with built binaries, turnkey
  `prepare` end-to-end on a planted-truth toy sparse system, probe vs the
  planted coefficients at a fresh prime, mutated-truth control — binaries
  absent → the end-to-end legs SKIP by name), plus `blade_family.py
  prepare`'s own end-to-end gates (held-out-prime oracle compare) on real
  families.

## Port vs fork — what replays upstream, what extends it

Blade here is a port plus a fork, and the boundary runs by function, not
by directory. The REPLAY half reproduces the upstream engine's behavior
and is judged by agreement with it: the compiled search and solve stages
ship as upstream wrote them (the Mathematica package layer is upstream's,
unmodified — the port removes Wolfram from the compiled toolchain, not
from the optional notebook interface); the C wrapper layer is adapted
function by function from FiniteFlow's MathLink wrappers (MIT, T. Peraro;
copyright notice retained) with the MathLink traffic replaced by plain
arguments; and the Python stage modules — formats, extension, ansatz,
scheme, search (including its reduced-analytic mode), adaptive, semibl —
port the corresponding Blade `.wl` sources (MIT) with per-function line
citations, to the point of replicating one upstream off-by-one comparison
verbatim. In this half, a
disagreement with upstream is a port bug, and the fork is pinned to the
upstream commit it replays. The FORK half extends the engine: the
no-Wolfram build itself; the streaming sibling driver (`fflowcli_stream`
— streaming evalmany writes plus subgraph build/load, compiled from its
own translation units so the stock driver stays byte-identical); the
per-point numeric reduction oracle (`probe`/BTOracle), which has no
upstream analogue; standalone rational reconstruction over cached
evaluation tables (`ratrec`); the measured artifact-gate discipline in
`pipeline` (upstream's binaries can exit 0 on failure paths, so byte
sizes, table dimensions, and held-out compares gate instead of return
codes); and a small set of declared divergences inside the ported
modules — deterministic block placement, two-point symmetry-map
confirmation, typed refusals where upstream would proceed — each
disclosed at its definition. In this half, new behavior is a feature
with its own gate, not a deviation to reconcile.

GATES (the load-bearing ones):

- Run the selftest after ANY blade/fork edit (the `--full` legs need the reference bank).
- Verification-grade probes: pass a FRESH prime — `probe` defaults to fit prime
  BIG_UINT[0].
- Gate ONLY on fit-table dims + eval byte sizes + dynamicrr state==2 + held-out
  compare — never on rc.
- Maximal-cut extension MUST declare `zero_sectors`: without it hundreds of
  thousands of cut-invalid integrals leak in.
- `banked-db` lane (`database_dir` only) SKIPS the held-out oracle unless
  `kira_table` is also given — supply both for a gated bank.
- Raw-system lane (`fflow_system`): the held-out gate replays the input
  system's OWN sparse solve at a fresh prime — it gates the search/fit/export
  chain, not the input system; validate the system itself upstream.  Masters
  must occupy the HIGHEST column ids (the eliminator pivots left to right),
  and all lane primes are fflow BIG_UINT primes (the binary refuses others).
- Respect the search STOP-line: projected-next-round refusal is deterministic.

FOOTGUNS (measured, all of them):

- `fitrel` writes flag=1 over an EMPTY fit nullspace (rc=0). Gate on fit-table dims.
- `ssolve` exits rc=0 while writing a wrong-SIZE eval. Gate on eval byte sizes.
- `ssolve` ACCEPTS freehand points — the structured-points constraint is recmod's
  only; feeding random points poisons downstream silently.
- `recmod` fails on OVER-estimated total degrees — learn degrees or use exact ones.
- `scan_degrees` ABSORBS corrupt scan points — keep held-out points independent.
- `evalmany` needs `--expect-nout`. depvars must be DESCENDING; nparsin≤9 via the
  `ptsord` digit encoding (use the evaluator lane above 9).
- `fflowcli` REFUSES non-BIG_UINT row primes (by design).
- Fresh-workdir refusal on `prepare` is deliberate — do not "resume" into a stale
  workdir.
- finiteflow MUST be built in-source; the blade fork source is the
  gitee-pinned `no-wolfram-port` branch in `blade` —
  never the deprecated gitlab repo.

The fork source lives in the sibling repository `blade`
(including the `fflowcli_stream` streaming sibling — see its `PATCHES.md`); the
reference bank of search results is not distributed.

CREDIT: Blade — X. Guan, X. Liu, Y.-Q. Ma, W.-H. Wu (arXiv:2405.14621);
block-triangular method lineage arXiv:1801.10523 + arXiv:1912.09294. Our port
removes the Wolfram dependency; algorithmic credit theirs. The fork replays
upstream commit c0a06cb (+2707/−8 across eight files; the full diff record is
`PATCHES.md` in `blade`); the C wrapper layer is adapted from
FiniteFlow's MathLink wrappers at FiniteFlow commit 278344c (T. Peraro, arXiv:1905.08019, MIT —
see `THIRD_PARTY.md` in `blade`).
