# seedling — guide

Tool page: https://bootloops.ai/tools/seedling.html

KIND: package (python package `seedling`: cli, support, pinner, margins,
bulkmodel, runner, certify, receipts, escalate, rankcheck, audit; console
script `seedling` via pyproject, or `python3 -m seedling.cli` with this
directory's parent on PYTHONPATH)

REQUIREMENTS: `pyyaml` beyond the common core (`sympy` is already in it) — the
cli and audit members import yaml. The pyproject declares it, so a pip install
of the package carries it; on the plain-PYTHONPATH route, `pip install pyyaml`.

PURPOSE: the pre-reduction front door for kira — certified predictive staging
plus family preflight. Staging pipeline: targets -> support -> census ->
margin prediction -> per-sector pinned schedule -> staged kira in a fresh
cgroup -> two-slice + DE-census gates -> receipts -> ledger. ONE command,
fail-closed. Family preflight (`audit` / `identity`): discrete-symmetry and
transcription checks BEFORE the first kira run on a family. Division of
labor with its siblings: seedling audits the inputs to a reduction; trust
audits reduction output.

USE-WHEN:
- Any deep-target kira staging that would otherwise be a hand-pinned seed-pin
  run — seedling is that recipe productized, with predicted margins and
  certification. Measured motivation: an unpinned staging of a 3-loop
  15-propagator system OOMed >313 GB; the support-pinned staging of the
  identical 161-sector system completed in 334 s at 7.42 GB, value-validated
  442/442 against two independent references.
- BEFORE the first kira run on any NEW family: `seedling audit` — kira's own
  symmetry finder never composes a loop map with an external-leg crossing, so
  it silently carries redundant masters (measured on every audited production
  family that had such a symmetry).
- On EVERY paper-transcribed family BEFORE probing: `seedling identity` — a
  wrong momentum convention (outgoing vs all-incoming) with the paper's
  propagators is a WRONG family (a 76-vs-98-masters class error). A
  dictionary hit is mandatory before probing.
- Need a pin-memory/bulk prediction before staging (`preflight`).
- Need a differential correctness gate between a fresh staged run and a
  reference (`certify`), or a DE-census coverage verdict + mini-kira
  escalation spec (`census`).
- Need the exact F_p elimination engine as a LIBRARY -> tools/winnow (ships
  in this repo).

NOT-FOR: does not choose targets or basis; no symbolic reconstruction; the
bulk predictor's constants are family-portable HYPOTHESES calibrated on one
production family (treat predictions as estimates with ~2x safety factor).

INVOKE:
- `seedling pin --config <kira dir> --targets F --out DIR [--margin r:0,s:0,d:0]
  [--preferred F] [--unsafe-s0] [--execute --cgroup NAME,MEM_MAX,CPU_MAX]
  [--kira-cmd "kira --parallel=4 jobs.yaml"]` (dry by default)
- `seedling run --config <kira dir> --targets F --out DIR [--margins m2|const]
  [--labels F] [--sector-maps F] [--preferred F] [--unsafe-s0]
  [--execute --cgroup ...]
  [--kira-cmd ...] [--fresh/--ref/--fresh2/--ref2 slice-row JSONs]
  [--eta-props i,j,...] [--receipt-hook MOD.py] [--slice-hook MOD.py]`
- `seedling preflight --census <kira nonTrivialSector file> --top ID --nsp N
  --nops N --cell r,s,d`
- `seedling certify --fresh F --ref F [--fresh2 F --ref2 F]` (two-slice
  differential gate, zero-normalized rows)
- `seedling census --masters JSON --eta-props i,j,... [--gated JSON]
  [--config DIR] [--emit-minikira DIR]`
- `seedling ledger --rundir DIR`
- `seedling audit CFG [--out DIR] [--max-shift 1] [--n-denom N] [--no-cache]
  [--report-crossings]` — CFG = kira config dir OR amflow JSON; full
  autopermutation group (external-leg crossings ∘ C ∘ affine unimodular loop
  maps that permute the denominator set), anchor-solved search, every
  generator re-verified by explicit sympy substitution
- `seedling identity --ours FAMILY_DEF.json --theirs Propagators.m
  [--theirs-format mathematica|json] [--theirs-key KEY] [--swap "l->k,k->p"]
  [--masses "p1^2=M2,..."] [--x-perm j1,..,jN] [--out F]` — same-family check
  via the affine scalar-product-dictionary / second-Symanzik-F route,
  arbitrary loop/prop count, full symbolic external Gram + mass monomials;
  rc 0 PASS / rc 1 FAIL-loud with the inconsistent monomial NAMED. No-arg
  `seedling identity` = the built-in wpair regression fixture.

INPUTS: kira config dir + target list. `--sector-maps` (a kira
sectorSymmetries/sectorRelations file) closes the target support under the
family symmetry maps before pinning (support.symmetry_closure): symmetry
images land in their mapped-to sectors with target-grade cells; images the
map file leaves undecidable are ledgered, never silently dropped. Slice-row
JSONs and the receipt
emitter are INJECTED (JSON paths / hook modules) — no caller-specific paths
inside the package. Margins: `--margins m2` wraps a frozen ridge
model — set SEEDLING_ML_DIR (directory with train_eval.py + features.py) and
SEEDLING_LABELS (labels JSON); these artifacts are not part of this package
and m2 mode fails loudly without them. `--margins const` (the constant
baseline, margin 6) is fully self-contained.

OUTPUTS: pinned `jobs.yaml` + run-dir; preflight bulk/RSS prediction;
SEL/masters extraction; two-slice + DE-census verdicts; receipts
(lambda-witnesses, via tools/winnow when used as the engine); JSONL telemetry
ledger. Audit: `<fam>_AUT.json` (|Aut|, generators, D-perms, ISP action),
`<fam>_aut_relations.kira` seed file, `<fam>_magic_relations.yaml` drop-in;
cached per family-hash beside the module (override with FAMILY_AUT_CACHE);
WARNs if |Aut|>1 and the kira config lacks `magic_relations: true`. Identity:
rc 0 PASS / rc 1 FAIL-loud, verdict JSON via `--out`. Every margin emission
carries provenance (model_hash, labels path+sha256, floor flags, mode).

ENV: nothing executes kira unless `--execute` AND a `--cgroup
NAME,MEM_MAX,CPU_MAX` spec are both given. Hard floors: meff >= 1; s >= chain
floor — `--unsafe-s0` is the only override, logged. Receipt backend routing
is budget/VA-aware (never dense-OOMs a receipt pass; one loud dense->cpu
fallback). Cgroup leaf launch (the cgroup.procs write discipline, documented
clause by clause in `runner.build_cgroup_script`'s docstring): the `--cgroup` spec is validated
BEFORE anything is built — NAME `^[A-Za-z0-9_.-]+$`, never a dot-only
component, MEM_MAX/CPU_MAX non-empty cgroup v2 literals — else rc 2 with a
one-line reason; the leaf script runs `set -eu`, `: "${CG:?}"` and a
`/sys/fs/cgroup/?*` + no-dot-component path guard before its first write (exit
97), writes only its own `$$`, self-checks `grep -qx "$$" cgroup.procs` and
`/proc/self/cgroup` before `exec` (else moves back, exit 98), prints the
`LEAF <cgroup> pid <pid> cmd=<cmd>` line; the launcher (Popen, never
timeout/setsid/nohup) reads the leaf while the worker runs and requires every
pid's ppid chain to end at the child it forked — a foreign pid -> SIGKILL the
child, every leaf pid moved back to the parent, ledger `phase=abort`, exit
status 99. `sudo -n` only (never prompts).

GATES: FAIL-CLOSED contract — any stage exception => ledger `phase=failed`
(stage+traceback) + rc != 0; two-slice verdict != PASS => rc != 0; gate steps
without inputs are ledgered SKIPPED with a reason, never silently omitted;
the ledger rejects gate emissions naming no evidence file. Reference record:
golden replays byte-exact; live staged run 20 s / 0.6 GB with masters ==
reference, SEL == reference, two-slice PASS, receipts CERTIFIED. Audit member:
<= ~11 s on all validated families; a 10-propagator production family's Z2
(external swap ∘ C ∘ loop map) found exactly; equal-mass sunrise S3;
asymmetric box trivial negative; one-mass mutation shrinks the group (never
silently passes). Identity member: an 18-propagator family vs its published
ancillary list — rank 10/10, det(L)=1, dictionary = identity, polynomial
verification == 0; one flipped propagator sign -> rc=1 at the named monomial.

FOOTGUNS:
- Gate hooks must ride proven oracle machinery — hand-keyed per-target gate
  rows are UNSOUND (survivor-master labeling is slice-set-sensitive).
- Never project pass2 receipt walls from pass1 (measured 2-3x apart;
  recorded separately by design).
- Any model/feature change to margins requires re-validation; margins.py never
  retrains silently.
- Degenerate inputs raise loudly — treat a quiet run with empty gates as
  suspect; check the ledger.
- Skipping `seedling audit` before kira is how redundant masters enter
  silently — kira never tries crossing ∘ C. Wire the emitted
  `magic_relations.yaml` / `.kira` seeds into the config; reading |Aut| and
  moving on buys nothing. A trivial Aut result does NOT certify basis
  integrity — embedding-unrealizable symmetries (e.g. S4 on an equal-mass
  subsector at a special kinematic point) are invisible to this scan.
- An identity FAIL means the transcription (or the convention) is wrong — do
  NOT proceed to probing and "fix it later"; the wrong family reduces fine
  and poisons everything downstream. The no-arg fixture ENUMERATES slot
  permutations; at 18 props that is not an option — use `--x-perm`.

BATTERY (ships, self-contained): `tests/` — 52 tests pass with no external
data (core pipeline units, symmetry-closure support, the self-contained
const margin arm end-to-end, fail-closed contract, receipt routing + gate
bookkeeping + degenerate-input guards, audit/identity on the shipped sunrise
and asym_box fixtures including a mutation control, and the cgroup leaf law's
planted controls in `test_cgroup_law.py`: empty/`../x` NAME -> rc 2 with
nothing built, guard tokens in the dry-run text, the shell path guard on
forced `$CG` values, a fixture leaf holding a foreign pid -> the abort path
with mocked SIGKILL/move and the ledger `abort` record, an own-descendants
positive control — those ten touch no real cgroup or sudo — and ONE live
control, `test_positive_control_live_leaf_size_zero_with_N_pids`: a real
cgroup v2 leaf made under the calling process's own cgroup, holding 3 `sleep`
children the test spawned, reports `os.path.getsize(cgroup.procs) == 0` while
the read lists exactly those 3 pids — cgroup v2 interface files stat at size
0 whatever they hold, so a `test -s cgroup.procs` liveness check is shown
inert beside the read-based one, and a kernel that ever reports a non-zero
size fails the test by name; no sudo, the leaf is torn down in `finally`
(children terminated by pid, rmdir); where the kernel refuses the mkdir or the
cgroup.procs write the test SKIPS by name with the errno); 2 further legs
exercise the frozen margin model and SKIP unless SEEDLING_ML_DIR +
SEEDLING_LABELS are provided. Run:
`PYTHONPATH=tools python3 -m pytest tools/seedling/tests -q`. Skip contract:
that line and the registered BATTERIES.json line are plain pytest, so a skip
leaves the exit status 0 and is named only under `-rs` (add it: `-q -rs`);
run_selftests.py therefore reports the package PASS with the live leg skipped.
`python3 tools/seedling/tests/test_cgroup_law.py` runs that module directly:
ALL PASS with the count only at 0 skipped, otherwise the skipped test is
named and the exit status is 2.

RELATED: tools/winnow (the exact F_p elimination engine seedling's receipts
ride, with the standalone witness verifier); tools/trust (certifies
reduction OUTPUT — the complementary gate); kira itself (the Kira fork lives in the
sibling repository kira, cloned beside this one as ../kira).

CREDIT: Kira — P. Maierhöfer, J. Usovitsch, P. Uwer [Kira1]; J. Klappert, F. Lange, P.
Maierhöfer, J. Usovitsch [Kira2]; F. Lange, J. Usovitsch, Z. Wu [Kira3] — is the engine
being staged, with FireFly [FF1, FF2] and Fermat (R. H. Lewis) beneath it.
Support-pinned seeding is in the same spirit as Kira 3's improved seeding and equation
selection [Kira3]; this package only predicts, pins, gates and certifies what Kira runs.
Exact elimination by winnow (Laporta [Lap]). Bracketed keys resolve in REFERENCES.md at
the repository root.