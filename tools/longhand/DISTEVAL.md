# longhand.disteval — one pySecDec disteval point as a checkpointed evaluation chain

Route B of Longhand (manual; the package guide is `GUIDE.md`). `longhand.disteval` evaluates a compiled
pySecDec `loop_package` at one Euclidean kinematic point as a chain of three disteval stages with per-sector
checkpoints, a three-way lattice compare and a conservative acceptance check, and prints the memory figure a
caller declares for the run.  It works on any compiled family (`--family NAME --pins FILE`); what it adds to
a bare `python3 -m pySecDec.disteval` call is resumability (a chunk of sectors is the checkpoint unit), a
second integrator on the same kernel library, a digit count that is a minimum over pairs and bars rather than
a lattice setting, and two controls that prove the compare can see a wrong digit.  All numerical work is
pySecDec's (GPL-3.0, installed by the user: `pip install pySecDec`; never bundled here); the staging, the
compare and the check are ours.  Command line: `tools/longhand/disteval/pysecdec_point.py`; helper command
lines `stage_chunks.py`, `assemble_compare.py`, `gate.py`, `memory_basis.py` in the same directory.

**Worked example.** The records shipped under `../examples/ndpent_top/` (sha256-pinned in `PINS.json`
there) are one complete run of this chain on the generic-mass five-point top-sector scalar
(`ndpent[1,1,1,1,1,1,1,1,0,0,0]`, the double pentagon with all eight propagators; 332 sectors, 968
kernels) at `s12=-7 s23=-9 s34=-13 s45=-10 s15=-5 mm3=-12 mm4=-8`:

    eps^0 = 8.246341248953e-01 +- 7.3e-08   (stage B, the value of record; 7 digits at the weakest order)

with the per-order table `eps^-4 = 9, eps^-3 = 7, eps^-2 = 7, eps^-1 = 7, eps^0 = 7` digits, the
conservative digit class 7, the leading pole to 9 digits — a double-precision quasi-Monte-Carlo figure,
not a 30-digit certificate.  That run's stage results, compare, acceptance check and cgroup readbacks are
the records the battery replays: the check over them must reproduce the record (`--stage check --example`:
21 checks by name).  Every measured figure quoted below ("on the example") is read from those records.

## What it does

| stage | setting | form |
|---|---|---|
| `A` | `--points 1e4 --shifts 32 --epsrel 1e-4`, standard rank-1 lattices | unchunked; the coarse consistency check; its per-order magnitudes set the per-chunk `epsabs` of B and C |
| `B` | `--points 1e5 --shifts 32 --epsrel 1e-7`, standard rank-1 lattices | chunked by `--chunk-sectors 12` (28 chunks of the example's 332 sectors); the value of record |
| `C` | the SECOND INTEGRATOR: disteval's median-lattice rule `--lattice-candidates 11` at the B targets | chunked; a different lattice construction on the same kernel library |
| `compare` | pairwise per order: sigma = \|diff\| / sqrt(err_a^2 + err_b^2), agreed digits = floor(-log10(\|diff\| / \|value\|)) | `assemble_compare.py compare` |
| `check` (older spelling `gate`, kept) | the acceptance check: three stages landed; B-C sigma > 3 on any order = FAIL line; digits(order) = min over the B-C pair and the B, C bars; the digit class = min over orders; A's sigma vs B and C listed (a > 3 sigma line is MARGINAL, not a FAIL); PLANTED-FAIL and FOREIGN controls | `gate.py` (the module, its `run_gate` and the output's `gate` key keep the historical spelling so the pinned record `GATE_record.json` reproduces byte for byte) |

A **chunk** is a copy of the disteval data dir whose integral json lists only the kernels of a sector
group (`stage_chunks.py`; the `.so` files hard-linked or symlinked, `coefficients/` symlinked): disteval
writes nothing per sector, so the chunk is the checkpoint unit.  Each chunk leaves
`chunk_NNN/{result.json, disteval.log, time.txt, rc, DONE}`; `--resume` skips DONE chunks and DONE
stages; `--max-chunks N` evaluates N chunks and stops (a pilot of the per-chunk wall; the stage stays
INCOMPLETE, rc 1 by name, and resumes); the assembly adds the chunk results per order with the bars in
quadrature (`assemble_compare.py assemble`).  The per-chunk `epsabs` of B and C is
`epsrel x min_order |value_A| / sqrt(n_chunks)` (the target of the whole sum split evenly in quadrature),
floored at `--epsabs` (default 1e-12): on the example `1.706e-10`.  There is NO clock anywhere (no
`timeout`, no disteval `--timeout`): a stage runs until its chunks land or it is stopped by pid, and
resumes.

`--dry` is the stand-in form: every leg runs except the disteval call, which becomes an
import/library-presence check under GNU time plus a copy of a canned real result
(`--dry-result`, default the example's smoke); it rehearses the chunking, the checkpoints, the assembly
and the resume in seconds.

### Controls (by name, inside the acceptance check)

- **PLANTED-FAIL**: a copy of RESULT_B with the eps^0 value shifted by 1e-3 relative is compared to A and
  C and MUST raise a FAIL line (sigma > 3); if the compare does not see it, it is blind and the check FAILS.
  On the example: sigma vs A 9.99707344790171, vs C 8133.10753740778.  `--planted-fail-control` prints the line.
- **FOREIGN**: coarse independent evaluations of the same library (lattice 1e3, 4 shifts, 2 workers,
  epsrel 1e-2; one on the build host, one on a second host, the same `.so` sha) vs stage A: every order
  within 3x the summed bars (max sigmas 1.717724232013687 and 2.9485284820009676 on the example).  With no
  smoke given the check FAILS on `foreign_ok` — a control is never assumed.
- a planted digit in an example RESULT re-pinned makes the check FAIL by name (rc 1); un-re-pinned it is a
  pin mismatch (rc 3); a missing pinned record is rc 4.

## The compiled package (not delivered)

The compiled package of the example — `loop_package` (pySecDec 1.6.6, geometric decomposition, no contour
deformation, requested order eps^0) + `make disteval`: `disteval/{builtin.so, coefficients/, ndpent_top.json,
ndpent_top_integral.json, ndpent_top_integral.so}` (13 MB) with 77 MB of generated sources — is a build
of the host it was compiled on and is NOT part of this package.  Two of the five pinned files, the descriptor
JSONs `ndpent_top.json` and `ndpent_top_integral.json` (pure metadata: parameter names, prefactor expansion,
sector/order kernel lists), are shipped under `../examples/ndpent_top/` as `package_ndpent_top*.json` so the
split/dry/pins tests run without a compiled package; they are output generated by pySecDec 1.6.6
(GPL-3.0; [pySecDec, disteval] in REFERENCES.md) for our integral and contain no pySecDec code.  The `.so`
libraries and `coefficients/` are not shipped.  `disteval/PACKAGE_PINS.json` pins the example's five disteval
files by sha256 plus the generate receipt; `--disteval-dir` must match them or the run is REFUSED by name
(rc 3); `--unpinned-package` runs a rebuilt or different package with one loud `UNPINNED PACKAGE` line
and every result labelled `package_pinned false` (the example smokes are then not a FOREIGN control:
give `--smoke NAME=smoke.json` of that library).  Measured compile class of the example package:
generate 39.6 s (python) + `make -j2 disteval` 7:57.84 wall (GNU time, 139% CPU) on the build host;
198.3 s on the second host.

### Your own family: `--family NAME` and `--pins FILE`

The package names are derived from one family name (default `ndpent_top`, the example's family; `--name`
is the alias): the sum file `NAME.json`, the
integral `NAME_integral.json` / `NAME_integral.so`, the `sums` key `NAME` of every disteval result and
smoke, the integral name `NAME_integral` on disteval's per-integral statistics line.  `--family NAME` sets
them on `pysecdec_point.py`; `stage_chunks.py`, `assemble_compare.py assemble` and `gate.py` take the same
flag (their library units `load(D, name)`, `write_chunk(..., name)`, `read_result(p, name)`,
`smoke_orders(path, name)` and `run_gate(..., family=)` default to the example's family, so every existing
call is unchanged).  A family NAME is valid iff it matches `[A-Za-z0-9_]+` (letters, digits, underscore;
no separators, no dots, no empty token): anything else is refused by name (rc 2) before any file is read,
on every arm and on the three helper CLIs.  A `--family` whose sum file is absent from `--disteval-dir` is
refused by name (rc 2, the families present listed); a sum file whose `integrals` entry is not
`[NAME_integral]`, an integral file whose `name` is not `NAME_integral`, or either not JSON, is refused by
name (rc 2) on every arm that reads the package (the stage lines, `--check-package`, `--emit-pins`,
`stage_chunks.py`; `stage_chunks.family_problems(D, NAME)` is the check).  A smoke that carries no
`sums.NAME` is a FOREIGN control that FAILS by name (an `error` row), never a skipped one; the assemble CLI
refuses by name (rc 2, nothing written) a landed chunk result without `sums.NAME` and a chunks dir without
its `CHUNKS_<stage>.json`.

`--pins FILE` names the PACKAGE_PINS of that family (the same JSON shape as `PACKAGE_PINS.json`, which
stays the default and the example family's): a mismatch is rc 3 by name (a pins file of another family
counts as a mismatch, naming both families), `--unpinned-package` runs labelled as before, a `--pins FILE`
that does not exist is rc 4 by name, and a `--pins FILE` that is not a JSON object with a `files` object
of relative path -> 64-hex sha256 strings (and a `family` string when present) is refused by name (rc 2;
`--unpinned-package` does not bypass the shape check).  `--emit-pins OUT --disteval-dir D [--family NAME]
[--generate-receipt FILE]` writes that shape from the five disteval files of `D` (sha256 by hashlib at
emit; the family, the sector and kernel counts, the generate receipt's sha256 when given, a PRODUCER
block): on the example package it reproduces the five shas of `PACKAGE_PINS.json`; a `--generate-receipt`
FILE that does not exist is rc 4 by name and nothing is written.  The `generate_receipt_sha256` field is
recorded, not checked: `--check-package` reads only `family` and `files`, so pins emitted with the example
copy `../examples/ndpent_top/GENERATE_RECEIPT_record.json` (sha256 `85125efcf4787c4e...`) carry that sha where
`PACKAGE_PINS.json` carries `907489c509404e6f...` (the receipt the pins of record were cut from is not the
example copy byte for byte) and still print PACKAGE PINNED on the package of record.

The refusals of the family and pins paths, by name:

| input | rc | the line |
|---|---|---|
| `--family` not of the form `[A-Za-z0-9_]+` (empty, `a/b`, `../x`, `x/`, `a.b`), any arm or helper CLI | 2 | `REFUSED: --family '...' is not a family name (the form is [A-Za-z0-9_]+: ...)` |
| `--family NAME` with no `NAME.json` in the package dir | 2 | `REFUSED: --family NAME names no files in D (NAME.json missing); families present: ...` |
| a sum file whose `integrals` is not `[NAME_integral]`; an integral file whose `name` is not `NAME_integral`; either not JSON | 2 | `REFUSED: --family NAME is not the family of the package in D: <file> <entry> is ..., not ...` |
| `--pins FILE` not a JSON object with a non-empty `files` object of relative path -> 64-hex strings, or `family` not a string | 2 | `REFUSED: --pins FILE is not a PACKAGE_PINS file: <what is wrong> (the shape: ...)` |
| a pinned package file differing or MISSING; a pins file of another family | 3 | `REFUSED: package pin mismatch in D on <files>[, pins family X != --family Y] (...)` |
| `--pins FILE` that does not exist | 4 | `REFUSED: pins file missing: FILE (...)` |
| `--emit-pins ... --generate-receipt FILE` that does not exist | 4 | `REFUSED: --generate-receipt FILE missing (...)` |
| `assemble` over a chunk result without `sums.NAME`, or a chunks dir without `CHUNKS_<stage>.json` | 2 | `REFUSED: assemble --family NAME over O: <file> carries no sums.NAME (keys: ...)` / `CHUNKS_<stage>.json missing in O` |

## The memory figure (cgroup) and its basis

The tool runs inside whatever cgroup the caller gives it and never creates one.  `--basis W`
prints the figure to declare:

    memory.max(W) = page_up(1.3 x (driver_maxrss + W x per_worker_maxrss))      cpu.max = W x 100000      swap 0

from the measured figures in `../examples/ndpent_top/BASIS_record.json`: per-worker maxrss 6500352 B (the max
VmHWM over the `pysecdec_cpuworker` pids of a two-worker smoke sampled every 0.2 s) and the python
driver's VmHWM 93089792 B.  The nominal per-worker form `page_up(1.3 x W x per_worker)` alone leaves no
room for the resident driver (headroom 15601664 B at W 8 for a 93089792 B driver) and would OOM at
spawn, so the corrected form above, which includes the driver, is the one used; both figures are listed
per W in the basis file.  W 8 carries the PROBE label: memory.max is 2x the basis (377241600 B on the
record run; measured peak 108417024 B, oom_kill 0); the raise-once knob is 2x memory.max.  Creating the
cgroup leaf itself (memory.max by value, swap 0, cpu.max = W x 100000, a cpuset) is the caller's job, and
so is reading it back: `FENCE_AT_SPAWN_record.json` is the readback the record run's wrapper made at spawn
(memory.max equal to the declared figure, swap 0, cpu.max = the quota, the body inside the leaf), and
`--leaf-stamp FILE` hands the caller's end-of-run readback (`key=value` lines: memory.peak, memory.max,
oom_kill, wall_s, maxrss_kb_over_chunks, stages_rc, cpu.stat) to the check, which fails on a non-zero
`oom_kill` and records the rest.

## Measured walls

| tier | width | wall | where |
|---|---|---|---|
| smoke (lattice 1e3, 4 shifts, epsrel 1e-2) | W 2 | 4.03 s / 1.82 s | the build host / the second host (the FOREIGN controls) |
| stage A | W 8 | 24.68 s | the record (unchunked) |
| stage B | W 8 | 2688.68 s of chunk wall (28 chunks; max 319.47 s) | the record |
| stage C | W 8 | 3886.26 s of chunk wall (28 chunks; max 395.81 s) | the record |
| the whole chain A -> B -> C -> compare | W 8 | 6603 s inside a 377241600 B memory.max | the record |
| smoke setting, this tool | W 2 | 0:03.53 (min:s) | the build host (within the smoke's bars; eps^-4 byte-equal) |
| stage A, this tool | W 4 | 1:00.96 (min:s) | the build host (every order within the record's bars; eps^-4 byte-equal) |
| stage B chunk 0 only (`--max-chunks 1`), this tool | W 4 | 6:06.26 (min:s) | the build host |
| `--dry` B / C (28 chunks each) | W 2 | 0:30.91 / 0:34.30 (min:s) | the stand-in form |

A stage is priced from its own progress lines (`progress.jsonl`: one record per chunk with the wall,
the maxrss and disteval's per-integral statistics), never from a lattice setting: a lattice setting
alone is never a digit count.

## Usage

    cd tools/longhand/disteval      # or give the path; the records are found relative to the script
    python3 pysecdec_point.py --disteval-dir PKG/disteval --stage A --workers 8 --work WORK
    python3 pysecdec_point.py --disteval-dir PKG/disteval --stage B --workers 8 --work WORK        # epsabs from WORK/stage_A/RESULT_A.json
    python3 pysecdec_point.py --disteval-dir PKG/disteval --stage C --workers 8 --work WORK
    python3 pysecdec_point.py --stage compare --results A=WORK/stage_A/RESULT_A.json B=... C=... --out WORK/COMPARE.json
    python3 pysecdec_point.py --stage check --work WORK --compare WORK/COMPARE.json --leaf-stamp LEAF.stamp --out WORK/GATE.json
    python3 pysecdec_point.py --stage check --example               # the worked-example leg: must reproduce GATE_record.json (older spelling: --stage gate --fixtures)
    python3 pysecdec_point.py --planted-fail-control | --basis 8 | --check-package --disteval-dir PKG/disteval
    python3 pysecdec_point.py --emit-pins PINS_NAME.json --disteval-dir PKG2/disteval --family NAME    # another family's pins
    python3 pysecdec_point.py --disteval-dir PKG2/disteval --family NAME --pins PINS_NAME.json --stage A|B|C --workers 8 --work WORK2
    python3 pysecdec_point.py --stage check --work WORK2 --family NAME --smoke host1=smoke1.json --smoke host2=smoke2.json --leaf-stamp LEAF.stamp --compare WORK2/COMPARE.json
    python3 -m pytest -q -rs -p no:cacheprovider ../tests           # this route's battery: 31 passed, 1 skipped in about a minute with pySecDec installed;
                                                                    # 29 passed, 3 skipped by name in ~8 s without it (../selftest.py runs it as route B)

Exit codes: 0 PASS; 1 FAIL by name (a stage incomplete, a check FAIL, the record not reproduced); 2 usage
or a refusal by name (a bad `--point`, a missing `--work`, a `--family` of the wrong form or naming no
files, a sum or integral file of another family, a `--pins` file of the wrong shape); 3 a package or
record pin mismatch; 4 a pinned record, a `--pins` file or a `--generate-receipt` file missing.  `--point 'name=value ...'` must name
every real parameter of the package's sum file exactly once (default: the point of record).  The test
subprocesses carry the interpreter's user site (`PYTHONUSERBASE`), so a scratch `HOME` never hides the
user-site pySecDec from the stand-in check.  From Python: `from longhand.disteval import gate, stage_chunks,
assemble_compare, memory_basis` (with `tools/` on `sys.path`; `check` is accepted as an alias of `gate`).

## Not established (for the worked example)

- a value at any Minkowski point (contour deformation not built into the package);
- the other 8 top-sector masters and orders eps^1..eps^4 (separate packages);
- digits beyond the double-precision ceiling of disteval (a 30-digit closure needs another route: an
  auxiliary-mass-flow evaluator, or route A of this package for finite Euclidean integrands of low
  effective dimension);
- the value has no independent reference on the record: its bar is its own three-way agreement
  (no independent published value to compare against);
- a second point / the remaining masters (each a separate `loop_package` + make at the measured compile class).

## Files

`disteval/pysecdec_point.py` (the command line), `disteval/stage_chunks.py` (the chunk splitter),
`disteval/assemble_compare.py` (the assembly and the compare), `disteval/gate.py` (the acceptance check and its
controls), `disteval/memory_basis.py` (the basis), `disteval/PACKAGE_PINS.json` (the example package's pins),
`tests/test_disteval.py` (the battery), `examples/ndpent_top/` + `PINS.json` (the worked example's records).  The
numeric units (the kernel filtering, the add-per-order / variances-in-quadrature assembly, the sigma and agreed-digits
rule, the min-over-pairs-and-bars digit class, the planted and foreign controls, the basis formula) are the
units of the evaluation of record, unchanged but for the family name the file-naming units take.

CREDIT: this route orchestrates pySecDec (S. Borowka, G. Heinrich, S. Jahn, S. P.
Jones, M. Kerner, J. Schlenk, T. Zirke [pySecDec]) and its disteval evaluator (G.
Heinrich, S. P. Jones, M. Kerner, V. Magerya, A. Olsson, J. Schlenk [disteval]); the QMC
lattice rules are those of Borowka et al. [pySecDecQMC] and the median-lattice
construction of Goda and L'Ecuyer [GL]. All numerical work is pySecDec's; the staging,
compare and check are ours. pySecDec is GPL-3.0 and installed separately. Bracketed keys
resolve in REFERENCES.md at the repository root.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).
