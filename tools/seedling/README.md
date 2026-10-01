# seedling — adaptive certified seeding for Kira IBP reduction

Wraps any Kira job (family definition + target list) in an adaptive seed-pin loop:

    measure target support -> pin seeds (s chain floor) -> stage (cgroup, telemetry)
    -> rank check -> two-slice value gate + DE-census gate -> escalate on failure

Motivation (measured, 3-loop 15-propagator production family): staging memory
scales with seed-box bulk, including a target-independent Select-phase load floor
(~156 GB); the unpinned staging OOMed >313 GB, while the support-pinned staging of
the identical 161-sector system completed in 334 s at 7.42 GB (>=43x censored vs a
dedicated 320 GiB-capped reference run), value-validated 442/442 against two
independent references.

Safety model (all measured exhibits, not hypotheticals):
- Pinned systems are SUBSETS of the full Laporta system: sound by construction.
- Incompleteness is NOT self-signaling: an s:0 pin solved cleanly onto a false
  master set (wrong values, no solver error). Hence the s chain floor in `pinner`
  and the MANDATORY two-slice differential gate.
- Overlap gates miss DE-feeding rows: the DE-census gate certifies the uncovered
  subset with a from-scratch mini-kira (measured 451 s / 4.7 GB class).
- Silent basis defects exist in both directions (short by 2 from a missed
  symmetry; long by 1 from selection-dependent redundancy): `rankcheck` verdicts
  are SHORT/LONG/OK/SUSPECT — never a silent OK without independent input.

Family preflight, both run BEFORE the first Kira job on a family:
- `seedling audit CFG` — full discrete autopermutation group (external-leg
  crossings ∘ charge conjugation ∘ affine unimodular loop maps). Kira's own
  finder never composes a loop map with a leg crossing; families with such a
  symmetry silently carry redundant masters. Emits `magic_relations` drop-ins.
- `seedling identity --ours ... --theirs ...` — proves (or refutes, loudly,
  naming the inconsistent monomial) that a transcribed propagator list is the
  same family as a paper's, via the affine scalar-product dictionary on the
  second Symanzik polynomial.

## CLI (examples; `pin`/`run` are dry by default)

    # emit a pinned jobs.yaml from a kira config dir + target list
    seedling pin --config CONFIGDIR --targets target --out RUNDIR \
                 [--margin r:0,s:0,d:0] [--preferred preferred] \
                 [--execute --cgroup NAME,12G,"400000 100000"]

    # full pipeline: support -> census -> margins -> schedule -> stage -> gates
    # (--sector-maps closes the target support under family symmetry maps)
    seedling run --config CONFIGDIR --targets target --out RUNDIR --margins const \
                 [--sector-maps sectormappings/FAM/sectorRelations]

    # predict staged bulk + RSS for a cell (two-phase memory model)
    seedling preflight --census sectormappings/<fam>/nonTrivialSector \
                       --top 255 --nsp 15 --nops 18 --cell 13,1,6

    # two-slice differential value gate on JSON row dicts (zero-normalized)
    seedling certify --fresh s1.json --ref banked1.json \
                     [--fresh2 s2.json --ref2 banked2.json]

    # family preflight
    seedling audit CONFIGDIR
    seedling identity --ours fam.json --theirs Propagators.m --swap "l->k,k->p"

`--margins m2` (the learned margin model) wraps a frozen ridge model:
set `SEEDLING_ML_DIR` and `SEEDLING_LABELS` to your model directory
and labels JSON; without them m2 mode fails loudly and `--margins const` is the
self-contained baseline.

Nothing executes Kira unless `--execute` and a `--cgroup NAME,MEM_MAX,CPU_MAX`
spec are both given.

The cgroup leaf launch obeys the cgroup.procs write discipline (documented
clause by clause in `runner.build_cgroup_script`'s docstring): the
spec is validated before anything is built (NAME must match
`^[A-Za-z0-9_.-]+$` and may not be a dot-only component; MEM_MAX and CPU_MAX
must be non-empty cgroup v2 literals — otherwise rc 2 with a one-line reason);
the leaf script runs `set -eu`, aborts on an empty `$CG` and on any path that is
not a proper child of `/sys/fs/cgroup` (exit 97), writes only its own `$$` to
`cgroup.procs`, verifies it is in the leaf before `exec` (else moves itself
back and exits 98), and prints a `LEAF <cgroup> pid <pid> cmd=<cmd>` launch
line; the launcher then reads the leaf while the worker runs and requires every
pid's ppid chain to end at the child it forked — a foreign pid means SIGKILL
the child, move every leaf pid back to the parent, a `phase=abort` ledger
record, and exit status 99. `sudo -n` only: a missing sudo grant fails, it
never prompts.

## Tests

    PYTHONPATH=tools python3 -m pytest tools/seedling/tests -q

52 tests are self-contained (11 of them the controls of the cgroup leaf law
in `tests/test_cgroup_law.py`: ten touch no real cgroup or sudo; the eleventh
makes a real leaf under the test process's own cgroup — no sudo — holding 3
`sleep` children it spawned, and pins that `cgroup.procs` stats at size 0
while the read lists exactly those 3 pids, so a `test -s cgroup.procs`
liveness check is inert while the read-based one is not; it skips by name
with the errno where the mkdir or the write is refused); 2 frozen-model legs
skip unless `SEEDLING_ML_DIR`/`SEEDLING_LABELS` are provided. Plain pytest
leaves a skip at exit status 0 and names it only under `-rs`; running
`tests/test_cgroup_law.py` directly exits 2 on a skip and prints ALL PASS
only at 0 skipped.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).
