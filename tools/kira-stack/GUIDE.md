# kira-stack — guide

kira-stack — `tools/kira-stack` — https://bootloops.ai/tools/kira.html

KIND: external-wrapper + recipe-collection

WHAT SHIPS HERE: the kira extension scripts — `parallel_kira_gen.py` (Route-B sharded
IBP generation + engine-free top-sector emit), `parallel_kira_gen_uds.py` (targeted
low-r single-sector user-defined system emits), `build_restage.py` (per-sector
seed-pin restager) and `kira_gpl_tools.py` (loader for the two GPL tools below).
Flat `tools/parallel_kira_gen.py` / `tools/parallel_kira_gen_uds.py` are aliases of the same scripts.
NOT here, GPL: `pyred_weight.py` (pyred integral-weight codec, integral_ordering
1..8, used by top-emit and the census) and `ffsave_degree_census.py` (ff_save
degree-census decoder, zero solver compute) were written by reading Kira's and
FireFly's GPL sources, so they are GPL-3.0-or-later and ship in the sibling
repository kira under `bootloops-tools/`, with the Kira/FireFly
attribution in their headers. `kira_gpl_tools.py` finds them in
`../kira/bootloops-tools` (side-by-side checkouts) or in the directory
named by env `BOOTLOOPS_KIRA_TOOLS`; without them top-emit stops with a named
ImportError and the battery legs that need them skip by name.
The ENGINE is not here: use the Kira fork, which lives in its own repository,
the sibling repository kira (parallel-treatcoeff + SQLite guards) or
upstream Kira + FireFly + Fermat, obtained upstream.

PURPOSE: IBP reduction → masters + raw DE (kira driven via YAML / kira2math,
Fermat/fer64 back-sub, FireFly finite-field reconstruction), plus the operational
recipes that keep it alive at program scale: parallel generation, η-deformed SIG6 fix,
seed-pin staging, checkpoint-resume, basis-integrity caveats.

USE-WHEN:
- Any family needs a reduction table / DE rows and no precomputed table exists.
- `run_initiate:true` pegged at ~1.3 cores despite `--parallel=N` (N idle fer64) →
  parallel_kira_gen Route-B.
- Staging OOM/bad_alloc during Select/Generate on a deep-target job → seed-pin recipe
  (build_restage.py mechanizes it).
- Monolithic η-deformed reduce dies SIG6 → dot-cap + lower seed + FireFly (the fix
  lives in the amflow-cpp fork, the sibling repository `amflow-cpp`).
- kira generation censored/crawling on a monolithic reduce row on a ~1000+-sector tree →
  build_restage.py (pinning is cheap — seconds — and never hurts).
- kira+FireFly stuck on huge-degree entries and you need WHICH masters/targets carry the
  degree explosion → ffsave_degree_census.py (GPL, kira/bootloops-tools).
- Machine/process died mid back-sub with hours at stake → kira.db checkpoint-resume.

NOT-FOR: symbolic multi-var FireFly walls (→ tools/numkin); one-sector rows without
full-family IBP; independent verification of a table you must trust (→ tools/trust, receipt member).
uds route: not usable as a general reducer (13× slower, +85 spurious masters) — chunked
single-sector/gen-shard use only.

INVOKE:
- `python3 tools/parallel_kira_gen.py {setup|launch|status|top-emit|merge|probe}`
  with `--no-pm --targets-mode own-sectors --manifest`; the all-lines top sector
  is never shard-covered — `merge --top-from ORACLE_DIR` copies it from a
  completed run, or `top-emit OUT` generates SYSTEM_<fam>_<TOP>.gz engine-free
  (route-A operator templates + the pyred_weight codec; raw emit = no mod-p
  selection, no symmetry rows — gate the merged solve as always) and
  `merge --top-from OUT/top_emit` consumes it.
- `python3 tools/parallel_kira_gen_uds.py` (`--io/--no-pm`) for targeted low-r
  single-sector emits.
- `python3 tools/kira-stack/build_restage.py --base DIR --fam FAM TOP [TOP...]` or
  `--dir DIR --fam FAM`; `--nprops N` MUST match the family.
- `python3 ../kira/bootloops-tools/ffsave_degree_census.py <rundir>
  [--vars d,s2,s12] [--workers N] [--out census.json]` (or `$BOOTLOOPS_KIRA_TOOLS/...`)
  — decodes a stuck/live kira+FireFly `ff_save/states` into a per-(target,master)
  max-deg num/den census with ZERO solver compute. GPL tool, sibling repository.
- Resume: back up `results/`, rerun kira MANUALLY in the same dir; proof =
  `Equations reduced: 0 / (TOTAL−committed)`.

ENV: `BOOTLOOPS_KIRA_TOOLS` = directory holding the GPL `pyred_weight.py` /
`ffsave_degree_census.py` when `../kira/bootloops-tools` is not where they
are; `-p`/`--parallel` ≤8; jemalloc preload recommended (see BUILD.md in the sibling repository `amflow-cpp`);
memory-cap kira, no swap; `FERMATPATH` for fer64.

BATTERY (smoke, engine-free): `python3 selftest.py` (in this dir) — exercises the
shipped scripts on synthetic fixtures in a scratch directory: build_restage census /
pinned rows / receipt vs hand computation + the positive-ISP fail-closed control;
parallel_kira_gen helpers (cover_tops, parse_targets, synthetic own-sector targets);
ffsave_degree_census layout gate + planted-state decode (file + inverse naming, ZERO
states) swept over integral_ordering 1..8 + pyred_weight codec pins (hand-computed
weight, dotsp separation, sector ranks, fail-loud decode) + the out-of-range-ordering
fail-closed control; top-emit on a toy operator-template family (emitted SYSTEM
blocks vs hand computation, merge --top-from consumption) + the ordering-mismatch
fail-closed control. The census/codec/top-emit legs (S3, M2, S4, M3) need the GPL
tools from kira/bootloops-tools and SKIP BY NAME when they are absent
(the battery then also checks that the loader refuses with a named ImportError);
exit 0 either way. The engine itself is external (see WHAT SHIPS HERE).

GATES: never consume a fresh gen/merge without a body-SET compare vs a trusted serial
reference (dict-by-head is WRONG — SYSTEM has multiple eqs per head). Solve-gate after
seed-pin: Regenerate loop / leftover targets = seed clipping ⇒ rerun r+1 or d+1 (both
measured benign, eqn count moves ≤0.03%) and re-gate. Independent oracle on any table
that feeds downstream.
build_restage completeness arbiter: kira's own "unreduced integrals: 0" + (shard
farms) a dedup-UNION coverage arbiter.
ffsave_degree_census: LAYOUT round-trip (inverse(encode(name)) == name, assert) + the
exact-set file gate; the full ordering range 1..8 is implemented (5 and 8
production-validated; the rest follow the same layout and are held by the same
round-trip gate).
top-emit: raw emit carries no mod-p selection and no symmetry rows — always solve-gate
the merged system (kira "unreduced integrals: 0" + masters_check + body-SET compare
where an oracle exists); a rank deficit shows up as leftover targets, never silently.

FOOTGUNS (the operational core):
- η-deformed SIG6 = uncaught `std::bad_alloc` (Fermat back-sub swell ~115 GB, self-caps
  at RLIMIT_AS → SIG6 not SIG9). The fix is the COMBINATION, default-ON for η families
  in the amflow-cpp fork: (1) per-subsector dot-cap emitting `d: ibp_dot` (the
  load-bearing half), (2) auto-floored lower seed, (3) FireFly. FireFly alone is NOT
  the fix.
  Non-η families: strict no-op.
- Seed-pin staging: **s (ISPs) is the explosive dimension** (N_σ =
  C(min(d,r−t)+t,t)·C(s+m,m)) — never widen s above target support; the Select bulk-load
  floor is target-INDEPENDENT so chunking targets can NEVER fit under it (measured
  NO-GO). Pinned monolith datum: 313 G OOM → 7.42 G / 334 s. A pinned jobs.yaml CHANGES
  the ibp-cache key — rekey (see the keypred member of tools/amflow-kit) or use the staged dir directly.
- parallel_kira_gen as-shipped pm-copy defeats sharding — always `--no-pm` (measured
  2.64×/-p1, 1.49×/-p8 real). But no-pm gen SHIFTS the solve basis vs pm-generated artifacts (masters swap; rank preserved) — do not mix shard output with pm-generated
  coefficient artifacts without a basis reconciliation. Synthetic own-sector targets
  over-generate ~10–27× eqs vs real-target selection.
- uds: omit `integral_ordering` in UDS_JOBS = weight corruption; the >6 h wall is kira's
  INTERNAL forward elimination (superlinear in rows) — chunk to ~200-eq encodes; uds
  does NOT substitute the 2-term symmetry rows; no per-sector lock (disjoint worklists
  only); `cfg=0` unions pass `-s` tests; `export -f` snapshots defeat live patching.
- build_restage: --nprops wrong ⇒ sector bitmask and ISP assert are garbage; it
  fail-closes on positive "ISP" indices but a too-large nprops silently misbins — pin it
  from the family definition. Restaged jobs.yaml KEYS DIFFERENTLY by construction —
  inject cache products at the ORIGINAL MISS key. Naive relaunch of a censored monolith
  restarts the sweep from sector 1 — restage FIRST, then fire once.
- ffsave_degree_census: live kira rotates state files — snapshot `ff_save/states` FIRST.
  Ordering-5 layout: d1 = dots (not dots+sps as in ordering 8). Do NOT assume tags ⊆ run
  name files — DE-seed runs tag dotted/corner integrals that live in no file. The
  trivial-sector auto-search is exponential in distinct preferred sectors — capped at
  r<=2; pass --trivial-sectors explicitly for larger trivial sets. Custom-rank weights
  (small ints) are preferred-file ORDER ranks, not encodings — never re-encode them.
  ZERO-state files are counted and skipped.
- Asymmetric-η Kira cycle bug: `iteration: >500, items left: <30` frozen for 1000s of
  iterations = will NEVER finish — kill by PID.
- Family embeddings can carry a symmetry NOT realizable as momentum shifts — kira's
  preferred basis is then rank-deficient and reductions onto it are NON-UNIQUE; QUOTIENT
  the basis before ANY coefficient-level cross-run consumption.
- masters_preheat OVERCOUNTS masters — take `masters.final`; kira2math refuses rules for
  ever-preferred integrals (validate via basis rotation).
- Mass-slice farms: staged-SYSTEM reuse segfaults kira write-out — FRESH-GEN per slice;
  `preferred_masters` + `integral_ordering: 8` fixes FireFly "Validation failed: Entry 0"
  + export crash.
- Two-primes-agree is NOT a correctness gate (internally consistent, CRT-stable, wrong)
  — use an independent exact oracle.
- Never blanket-pkill kira/fer64 — kill only by verified PID/cwd.
- FF-SAVE RESUME GATE: a kira/FireFly `ff_save` WITHOUT the `shift` file is NOT
  resumable ("Shift file not found!") — the shift persists only after FireFly's first
  full checkpoint. Gate every resume on `ff_save/shift` presence BEFORE firing;
  engagement verdicts by STATES-READBACK (checkpoint counter movement), never Done-count
  (Done stays 0 through whole primes on big systems). FireFly re-checkpoints only at
  wave/prime boundaries — type tranche caps accordingly.
- GENERATION HAS NO CHECKPOINT: a kill mid sector-job generation loses everything since
  launch — tranche-resume covers FF phases only. Cap generation to complete the phase or
  don't fire it (shard generation across a farm instead).

Alternate 128-bit-weight build (kira128): a build flag of the same fork —
reach for it ONLY on the abort `Integral weight representation exceeds 64 bits`; never
as default.

CREDIT: Kira was created by Philipp Maierhöfer, Johann Usovitsch and Peter Uwer [Kira1];
Kira 2 with Jonas Klappert and Fabian Lange [Kira2]; Kira 3 by Fabian Lange, Johann
Usovitsch and Zihao Wu [Kira3] (GPL-3.0-or-later; we thank the Kira developers — our
fork's two patches are documented in `kira/PATCHES.md` and we intend to propose them
upstream). FireFly is by Jonas Klappert and Fabian Lange [FF1] with Sven Yannick Klein
[FF2]; Fermat by Robert H. Lewis (home.bway.net/lewis). The reduction strategy is
Laporta's [Lap] on Chetyrkin–Tkachov integration by parts [CT] with finite-field
reconstruction after von Manteuffel–Schabinger [vMS] and Peraro [Per16].
`pyred_weight.py` and `ffsave_degree_census.py` derive from Kira's pyRed weight layout
(pyred/integrals.cpp, kira/ReadYamlFiles.cpp) and FireFly's save-state format; they
are therefore GPL-3.0-or-later, credited to the Kira Developers and the FireFly
authors in their headers, and live in kira/bootloops-tools, not in this
MIT tree. The recipes and the parallel-generation scripts are ours on top.
Bracketed keys resolve in REFERENCES.md at the repository root.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).
