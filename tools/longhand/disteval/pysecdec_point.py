#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision.
"""pysecdec_point.py -- the command line of Longhand's disteval route (longhand.disteval): one pySecDec disteval evaluation of a compiled
loop_package at one Euclidean point as an evaluation chain with per-sector checkpoints, a three-way lattice compare and a conservative
acceptance check.  The generic-mass five-point double-pentagon top sector (family ndpent_top) is the WORKED EXAMPLE whose records ship under
../examples/ndpent_top/; any compiled family runs the same way (--family NAME --pins FILE).

  STAGES (each a --stage run over the same --disteval-dir and --point; each resumable from its chunk checkpoints):
    A   lattice setting 1: --points 1e4 --shifts 32 --epsrel 1e-4, UNCHUNKED (one result; its per-order magnitudes set the per-chunk epsabs of B and C)
    B   lattice setting 2: --points 1e5 --shifts 32 --epsrel 1e-7, chunked by --chunk-sectors (12 sectors -> 28 chunks of the fixture's 332 sectors)
    C   the SECOND INTEGRATOR: disteval's median-lattice rule (--lattice-candidates 11) at the B targets, chunked
    compare   the pairwise per-order compare of RESULT jsons (assemble_compare.compare)
    check     the acceptance check by object (alias: gate) over a work dir (--work) or the worked example's records (--example, alias --fixtures):
              three stages landed, B-C sigma > 3 = FAIL line, digit class = min over the B-C pair and the B, C bars per order, the PLANTED-FAIL
              control, the FOREIGN control (gate.run_gate; the output file keeps its historical name GATE.json, matching GATE_record.json)
  A chunk = a copy of the disteval data dir whose integral json lists only the kernels of a sector group (stage_chunks.py); its outputs
  chunk_NNN/{result.json, disteval.log, time.txt, rc, DONE} are the checkpoint; --resume skips DONE chunks; the assembly adds per order with the
  bars in quadrature (assemble_compare.assemble); NO clock anywhere (no timeout, no --timeout).  --dry = the stand-in form: every leg runs except
  the disteval call, which becomes an import/library-presence check under GNU time plus a copy of a canned real result.
  THE PACKAGE is not part of this tool: --disteval-dir names a compiled package (loop_package + make disteval); it is REFUSED by name (rc 3) when its
  five files do not match the sha256 pins of the package of record in PACKAGE_PINS.json (beside this file: the worked example's family), unless
  --unpinned-package is given (then one loud UNPINNED line and every result is labelled package_pinned false).
  THE FAMILY: --family NAME (default ndpent_top, the worked example's family; --name is its alias) derives every package name: the sum file NAME.json, the
  integral NAME_integral.json / NAME_integral.so, the sums key NAME, the integral name NAME_integral (the progress statistics line).  A NAME is
  valid iff it matches [A-Za-z0-9_]+ (no separators, no dots, no empty token): anything else is REFUSED by name (rc 2) before any file is read.
  A NAME whose sum file is absent from --disteval-dir is REFUSED by name (rc 2, the families present listed); a sum file whose 'integrals' entry
  is not [NAME_integral], or an integral file whose 'name' is not NAME_integral, or either not JSON, is REFUSED by name (rc 2) on every arm that
  reads the package (the stage lines, --check-package, --emit-pins; the helper CLIs alike).  --pins FILE names the PACKAGE_PINS of that family
  (the same JSON shape; default PACKAGE_PINS.json beside this tool = the worked example's family): a pins file whose family is not --family is a pin
  mismatch (rc 3, or the labelled UNPINNED run); a --pins FILE that does not exist is rc 4 by name; a --pins FILE that is not a JSON object with
  a "files" object of relative path -> 64-hex sha256 strings (and a "family" string when present) is REFUSED by name (rc 2).  --emit-pins OUT
  writes that shape from the five disteval files of --disteval-dir (sha256 by hashlib at emit; the family recorded; --generate-receipt FILE's
  sha256 recorded when given, else null; a --generate-receipt FILE that does not exist is rc 4 by name, nothing written).
  --basis W prints the memory figure the CALLER declares on its fence (memory_basis.py): page_up(1.3 x (driver + W x per_worker)) from the
  measured figures of record; this tool runs inside whatever fence it is given and never creates one.
  exit codes: 0 PASS / 1 FAIL by name (a stage incomplete, an acceptance-check FAIL, a record not reproduced) / 2 usage or refusal by name (a bad --point, a
  --family of the wrong form or naming no package, a sum or integral file of another family, a --pins file of the wrong shape) / 3 a package or
  fixture pin mismatch / 4 a pinned fixture, a --pins file or a --generate-receipt file missing.
usage: pysecdec_point.py --disteval-dir D --stage A [--workers W] [--work DIR] [--point 's12=-7 ...'] [--resume] [--dry] [--max-chunks N]
       pysecdec_point.py --disteval-dir D --stage B|C --workers W --work DIR [--chunk-sectors 12] [--epsabs-from DIR/stage_A/RESULT_A.json]
       pysecdec_point.py --stage compare --results A=RESULT_A.json B=RESULT_B.json C=RESULT_C.json --out COMPARE.json
       pysecdec_point.py --stage check (--example | --work DIR [--smoke NAME=smoke.json ...] [--leaf-stamp FILE]) [--out GATE.json]   (aliases: --stage gate, --fixtures)
       pysecdec_point.py --planted-fail-control | --basis W | --check-package --disteval-dir D [--pins FILE] | --emit-pins OUT --disteval-dir D
       every line above takes [--family NAME] [--pins FILE] (defaults: ndpent_top, PACKAGE_PINS.json beside the tool)"""
import os, sys, json, argparse, hashlib, subprocess, math, shutil, types, re

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import stage_chunks, assemble_compare, gate as gate_mod, memory_basis  # noqa: E402

# the worked example's records (sha256-pinned in PINS.json there): the package root is one level up from this member directory
FIX = os.path.join(os.path.dirname(HERE), 'examples', 'ndpent_top')
STAGE_ALIASES = {'check': 'gate'}   # --stage check is the acceptance check; 'gate' is its older spelling, kept working
PACKAGE_PINS = os.path.join(HERE, 'PACKAGE_PINS.json')
FIXTURE_PINS = os.path.join(FIX, 'PINS.json')
DEFAULT_FAMILY, BUILTIN_SO = 'ndpent_top', 'builtin.so'


class Family:
    """the family-derived names of a compiled package: the sum file NAME.json, the integral NAME_integral.json / NAME_integral.so, the
    sums key NAME, the integral name NAME_integral, the coefficient file, the progress literal disteval prints per integral."""
    def __init__(self, name):
        self.name = name
        self.integral = name + '_integral'
        self.sum_json = name + '.json'
        self.integral_json = name + '_integral.json'
        self.integral_so = name + '_integral.so'
        self.coefficient0 = 'coefficients/' + name + '_integral_coefficient0.txt'
        self.progress_literal = '- ' + name + '_integral:'

    def pinned_files(self):
        return [BUILTIN_SO, self.sum_json, self.integral_json, self.integral_so, self.coefficient0]


FAMILY_OF_RECORD = Family(DEFAULT_FAMILY)
SUM_JSON, INTEGRAL_JSON, INTEGRAL_SO = FAMILY_OF_RECORD.sum_json, FAMILY_OF_RECORD.integral_json, FAMILY_OF_RECORD.integral_so
POINT_OF_RECORD = 's12=-7 s23=-9 s34=-13 s45=-10 s15=-5 mm3=-12 mm4=-8'
STAGE_DEFAULTS = {'A': {'points': '1e4', 'epsrel': '1e-4', 'lattice_candidates': 0, 'chunk_sectors': 0},
                  'B': {'points': '1e5', 'epsrel': '1e-7', 'lattice_candidates': 0, 'chunk_sectors': 12},
                  'C': {'points': '1e5', 'epsrel': '1e-7', 'lattice_candidates': 11, 'chunk_sectors': 12}}
FOREIGN_SMOKES = [('the smoke on the build host (lattice 1e3, 4 shifts, 2 workers, epsrel 1e-2)', os.path.join(FIX, 'smoke_build_host.json')),
                  ('the smoke on the second host (the same setting, the same .so sha)', os.path.join(FIX, 'smoke_second_host.json'))]
RC_PASS, RC_FAIL, RC_USAGE, RC_PIN, RC_MISSING = 0, 1, 2, 3, 4


def utc():
    return subprocess.run(['date', '-u', '+%Y-%m-%dT%H:%M:%SZ'], capture_output=True, text=True).stdout.strip()


def sha256(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def producer_block():
    return {'tool': 'longhand.disteval', 'files_sha256': {f: sha256(os.path.join(HERE, f))[:16] for f in ('pysecdec_point.py', 'stage_chunks.py', 'assemble_compare.py', 'gate.py', 'memory_basis.py')}, 'stamp_source': 'date -u'}


class Refusal(Exception):
    def __init__(self, rc, msg):
        super().__init__(msg)
        self.rc = rc


# ---------------------------------------------------------------- the family name and the pins file, by name
def validate_family(name):
    """--family NAME is valid iff it matches [A-Za-z0-9_]+ (no separators, no dots, no empty token): anything else is REFUSED by name (rc 2)
    before any file is read.  Returns the Family."""
    if not isinstance(name, str) or not re.fullmatch(r'[A-Za-z0-9_]+', name):
        raise Refusal(RC_USAGE, 'REFUSED: --family %r is not a family name (the form is [A-Za-z0-9_]+: letters, digits, underscore; no separators, no dots, no empty token)' % (name,))
    return Family(name)


def validate_pins(pins_path):
    """--pins FILE must exist (else rc 4 by name), be a regular readable file and a JSON object with a "files" object of relative path -> 64-hex
    sha256 strings (and a "family" string when present): a directory, an unreadable file or any other shape is REFUSED by name (rc 2).  Returns the parsed pins."""
    if not os.path.exists(pins_path):
        raise Refusal(RC_MISSING, 'REFUSED: pins file missing: %s (--pins FILE names the PACKAGE_PINS of the family; --emit-pins OUT writes one from a disteval dir)' % pins_path)

    def bad(why):
        return Refusal(RC_USAGE, 'REFUSED: --pins %s is not a PACKAGE_PINS file: %s (the shape: a JSON object with "files" {relative path: 64-hex sha256} and "family"; --emit-pins OUT writes one)' % (pins_path, why))
    if not os.path.isfile(pins_path):
        raise bad('not a regular file (%s)' % ('a directory' if os.path.isdir(pins_path) else 'not a file'))
    try:
        pins = json.load(open(pins_path))
    except (ValueError, UnicodeDecodeError) as e:
        raise bad('not JSON (%s)' % e.__class__.__name__)
    except OSError as e:
        raise bad('not readable (%s)' % e.__class__.__name__)
    if not isinstance(pins, dict):
        raise bad('the top level is a JSON %s, not an object' % type(pins).__name__)
    if 'files' not in pins:
        raise bad('no "files" key')
    if not isinstance(pins['files'], dict) or not pins['files']:
        raise bad('"files" is a JSON %s, not a non-empty object' % type(pins['files']).__name__)
    for rel, want in pins['files'].items():
        if not rel or os.path.isabs(rel) or '..' in rel.split('/'):
            raise bad('"files" key %r is not a relative path' % (rel,))
        if not isinstance(want, str) or not re.fullmatch(r'[0-9a-f]{64}', want):
            raise bad('"files" entry %s is not a 64-hex sha256 string' % rel)
    if 'family' in pins and not isinstance(pins['family'], str):
        raise bad('"family" is a JSON %s, not a string' % type(pins['family']).__name__)
    return pins


# ---------------------------------------------------------------- the package pin
def check_package(D, unpinned=False, quiet=False, pins_path=None, fam=None):
    """the five disteval files vs the pins file (--pins FILE; default PACKAGE_PINS.json beside the tool): a mismatch or a missing file is REFUSED
    by name (rc 3) unless --unpinned-package; a pins file of another family than --family counts as a mismatch by name; a missing pins file rc 4;
    a pins file of the wrong shape rc 2 (validate_pins)."""
    pins_path = pins_path or PACKAGE_PINS
    fam = fam or FAMILY_OF_RECORD
    pins = validate_pins(pins_path)
    rows = []
    ok = True
    family_ok = (pins.get('family') == fam.name)
    for rel, want in pins['files'].items():
        p = os.path.join(D, rel)
        got = sha256(p) if os.path.exists(p) else None
        rows.append((rel, want[:16], (got or 'MISSING')[:16], got == want))
        ok = ok and got == want
    if not quiet:
        for rel, w, g, m in rows:
            print('  package %-52s pin %s  dir %s  %s' % (rel, w, g, 'OK' if m else 'MISMATCH'))
        if not family_ok:
            print('  pins family %s != --family %s (%s)' % (pins.get('family'), fam.name, os.path.basename(pins_path)))
    if ok and family_ok:
        if not quiet:
            print('PACKAGE PINNED: %s matches the fixture of record (%s; %d files; generate receipt %s)' % (D, pins['family'], len(rows), (pins.get('generate_receipt_sha256') or 'none')[:16]))
        return True
    on = ', '.join(r for r, _, _, m in rows if not m) + ('' if family_ok else '%spins family %s != --family %s' % (', ' if not ok else '', pins.get('family'), fam.name))
    if unpinned:
        print('UNPINNED PACKAGE: %s does NOT match the fixture of record on %s -- every result of this run is labelled package_pinned false; the fixture smokes are NOT a valid FOREIGN control for it (give --smoke of this library)' % (D, on))
        return False
    raise Refusal(RC_PIN, 'REFUSED: package pin mismatch in %s on %s (the fixture of record is pinned in %s; a rebuilt or different package runs only with --unpinned-package, labelled, or with the pins of its own family: --pins FILE)' % (D, on, os.path.basename(pins_path)))


def check_family_files(D, fam):
    """--family NAME must name the package in --disteval-dir: the sum file NAME.json (JSON; 'type' sum; 'integrals' == [NAME_integral]) and the
    integral file NAME_integral.json (JSON; 'type' integral; 'name' == NAME_integral) -- stage_chunks.family_problems is the check.  The sum file
    absent = REFUSED by name (rc 2) listing the families present; any other defect = REFUSED by name (rc 2) naming the file and the entry."""
    probs = stage_chunks.family_problems(D, fam.name)
    if not probs:
        return True
    if probs[0] == fam.sum_json + ' missing':
        present = sorted(f[:-len('_integral.json')] for f in os.listdir(D) if f.endswith('_integral.json')) if os.path.isdir(D) else []
        raise Refusal(RC_USAGE, 'REFUSED: --family %s names no files in %s (%s missing); families present: %s' % (fam.name, D, fam.sum_json, ', '.join(present) or 'none'))
    raise Refusal(RC_USAGE, 'REFUSED: --family %s is not the family of the package in %s: %s' % (fam.name, D, '; '.join(probs)))


def emit_pins(D, fam, out, generate_receipt=None, point=None):
    """--emit-pins OUT: the PACKAGE_PINS shape (family, integral, point_of_record, n_sectors, n_kernels, files {rel: sha256}, generate_receipt_sha256,
    generator, note, stamp_utc) from the five disteval files of --disteval-dir, hashed at emit; a missing file is REFUSED by name (rc 2); a
    --generate-receipt FILE that does not exist is REFUSED by name (rc 4); nothing is written before the checks pass."""
    check_family_files(D, fam)
    if generate_receipt and not os.path.isfile(generate_receipt):
        raise Refusal(RC_MISSING, 'REFUSED: --generate-receipt %s missing (its sha256 would be recorded in the pins: give the receipt file or omit the option)' % generate_receipt)
    files = {}
    for rel in fam.pinned_files():
        p = os.path.join(D, rel)
        if not os.path.exists(p):
            raise Refusal(RC_USAGE, 'REFUSED: --emit-pins: %s missing in %s' % (rel, D))
        files[rel] = sha256(p)
    I = json.load(open(os.path.join(D, fam.integral_json)))
    secs = stage_chunks.sectors_in_order(I)
    pins = {'family': fam.name, 'integral': fam.integral, 'point_of_record': point, 'n_sectors': len(secs), 'n_kernels': len(I['kernels']), 'files': files,
            'generate_receipt_sha256': sha256(generate_receipt) if generate_receipt else None,
            'generator': I.get('generated_by'), 'emitted_from': os.path.abspath(D),
            'note': 'emitted by pysecdec_point.py --emit-pins from the five disteval files of the dir named (sha256 by hashlib at emit); the package is not delivered with the tool',
            'stamp_utc': utc(), 'PRODUCER': producer_block()}
    json.dump(pins, open(out, 'w'), indent=1)
    print('PINS EMITTED: %s (family %s; %d files; %d sectors, %d kernels; generate receipt %s) -> %s (%s)' % (D, fam.name, len(files), len(secs), len(I['kernels']), (pins['generate_receipt_sha256'] or 'none')[:16], out, sha256(out)[:16]))
    return RC_PASS


def check_fixtures(names):
    pins = json.load(open(FIXTURE_PINS))
    for n in names:
        p = os.path.join(FIX, n)
        if n not in pins:
            raise Refusal(RC_MISSING, 'REFUSED: fixture %s is not pinned in %s' % (n, os.path.join(FIX, 'PINS.json')))
        if not os.path.exists(p):
            raise Refusal(RC_MISSING, 'REFUSED: pinned fixture missing: %s' % p)
        if sha256(p) != pins[n]:
            raise Refusal(RC_PIN, 'REFUSED: fixture pin mismatch: %s (sha256 %s, pinned %s)' % (n, sha256(p)[:16], pins[n][:16]))


# ---------------------------------------------------------------- the point
def parse_point(s, D, fam=None):
    """'name=value ...' against the sum file's realp list: every real parameter exactly once, each value a number; anything else rc 2 by name."""
    S = json.load(open(os.path.join(D, (fam or FAMILY_OF_RECORD).sum_json)))
    realp = list(S.get('realp', []))
    toks = s.split()
    seen = {}
    for t in toks:
        if '=' not in t:
            raise Refusal(RC_USAGE, 'REFUSED: --point token %r is not name=value' % t)
        k, v = t.split('=', 1)
        if k not in realp:
            raise Refusal(RC_USAGE, 'REFUSED: --point names %r which is not a real parameter of the package (%s)' % (k, ', '.join(realp)))
        try:
            float(v)
        except ValueError:
            raise Refusal(RC_USAGE, 'REFUSED: --point value %r for %s is not a number' % (v, k))
        if k in seen:
            raise Refusal(RC_USAGE, 'REFUSED: --point gives %s twice' % k)
        seen[k] = v
    missing = [k for k in realp if k not in seen]
    if missing:
        raise Refusal(RC_USAGE, 'REFUSED: --point lacks %s' % ', '.join(missing))
    return ['%s=%s' % (k, seen[k]) for k in realp]


# ---------------------------------------------------------------- the per-chunk epsabs (the chain's rule, in the chain's format)
def epsabs_per_chunk(result_json, epsrel, n_chunks, floor):
    """epsabs_chunk = epsrel x min_order |value| / sqrt(n_chunks) (the target of the whole sum, split evenly in quadrature); the floor when larger."""
    r = json.load(open(result_json)); eps = float(epsrel); n = int(n_chunks)
    vals = [abs(complex(o['value_re'], o['value_im'])) for o in r['orders'] if abs(complex(o['value_re'], o['value_im'])) > 0]
    used = ('%.3e' % (eps * min(vals) / math.sqrt(n)) if vals else str(epsrel))
    return '%.3e' % max(float(used), float(floor))


# ---------------------------------------------------------------- one stage
def hms_to_s(hms):
    parts = [float(x) for x in hms.split(':')]
    return sum(x * 60 ** i for i, x in enumerate(reversed(parts)))


def time_fields(timefile):
    wall = rss = None
    if os.path.exists(timefile):
        for ln in open(timefile).read().splitlines():
            if 'Elapsed (wall clock)' in ln:
                wall = ln.rsplit(' ', 1)[1]
            if 'Maximum resident' in ln:
                rss = int(ln.rsplit(' ', 1)[1])
    return wall, rss


def run_stage(a, pinned, fam=None):
    fam = fam or FAMILY_OF_RECORD
    S = a.stage
    d = STAGE_DEFAULTS[S]
    W = int(a.workers)
    points = a.points or d['points']; epsrel = a.epsrel or d['epsrel']; presamples = a.presamples; shifts = int(a.shifts)
    latcand = d['lattice_candidates'] if a.lattice_candidates is None else int(a.lattice_candidates)
    chunk_sectors = d['chunk_sectors'] if a.chunk_sectors is None else int(a.chunk_sectors)
    D = os.path.abspath(a.disteval_dir)
    point_args = parse_point(a.point, D, fam)
    WORK = os.path.abspath(a.work); SD = os.path.join(WORK, 'stage_%s' % S); os.makedirs(SD, exist_ok=True)
    resume = bool(a.resume); dry = bool(a.dry)
    print('STAGE %s START %s W=%d points=%s presamples=%s shifts=%d epsrel=%s epsabs=%s lattice-candidates=%d chunk_sectors=%d resume=%s dry=%s package_pinned=%s family=%s dir=%s' % (S, utc(), W, points, presamples, shifts, epsrel, a.epsabs, latcand, chunk_sectors, int(resume), int(dry), pinned, fam.name, SD))
    if os.path.exists(os.path.join(SD, 'STAGE_DONE')) and resume:
        print('STAGE %s already DONE (%s); skipped' % (S, open(os.path.join(SD, 'STAGE_DONE')).read().strip()))
        return RC_PASS
    open(os.path.join(SD, 'stage.start'), 'w').write(utc() + '\n')
    # chunks (kept on --resume; rebuilt otherwise)
    chunks_json = os.path.join(SD, 'chunks', 'CHUNKS_%s.json' % S)
    if not (os.path.exists(chunks_json) and os.path.getsize(chunks_json) > 0) or not resume:
        Ssum, I = stage_chunks.load(D, fam.name)
        secs = stage_chunks.sectors_in_order(I)
        groups = [secs] if chunk_sectors <= 0 else [secs[i:i + chunk_sectors] for i in range(0, len(secs), chunk_sectors)]
        os.makedirs(os.path.join(SD, 'chunks'), exist_ok=True)
        man = {'receipt': 'CHUNKS %s: per-sector checkpoint chunks of the disteval data dir (kernels filtered per sector group; prefactor/coefficient/sum unchanged; results add per order)' % S,
               'stamp_utc': utc(), 'stage': S, 'disteval_dir': D, 'name': fam.name,
               'source_sha256': {fam.sum_json: sha256(os.path.join(D, fam.sum_json)), fam.integral_json: sha256(os.path.join(D, fam.integral_json)), fam.integral_so: sha256(os.path.join(D, fam.integral_so)), BUILTIN_SO: sha256(os.path.join(D, BUILTIN_SO))},
               'n_sectors': len(secs), 'n_kernels': len(I['kernels']), 'chunk_sectors': chunk_sectors, 'n_chunks': len(groups), 'chunks': []}
        tot = 0
        for ci, g in enumerate(groups):
            cdir = os.path.join(SD, 'chunks', 'chunk_%03d' % ci)
            kern = stage_chunks.write_chunk(D, Ssum, I, cdir, g, fam.name)
            tot += len(kern)
            man['chunks'].append({'id': ci, 'dir': cdir, 'sectors': g, 'n_sectors': len(g), 'n_kernels': len(kern), 'integral_json_sha256': sha256(os.path.join(cdir, fam.integral_json))})
        man['n_kernels_in_chunks'] = tot
        man['kernels_cover_package'] = (tot == len(I['kernels']))
        json.dump(man, open(chunks_json, 'w'), indent=1)
        print('CHUNKS %s: %d chunks, %d sectors, %d kernels (package %d), cover=%s -> %s' % (S, len(groups), sum(len(g) for g in groups), tot, len(I['kernels']), man['kernels_cover_package'], os.path.join(SD, 'chunks')))
    man = json.load(open(chunks_json))
    NCH = man['n_chunks']
    # per-chunk epsabs from a previous stage's per-order magnitudes: epsabs_chunk = epsrel x min_order |value| / sqrt(n_chunks)
    epsabs_used = a.epsabs
    src = a.epsabs_from
    if src is None and S in ('B', 'C'):
        cand = os.path.join(WORK, 'stage_A', 'RESULT_A.json')
        src = cand if os.path.exists(cand) else None
    if src and os.path.exists(src) and os.path.getsize(src) > 0:
        epsabs_used = epsabs_per_chunk(src, epsrel, NCH, a.epsabs)
        print('per-chunk epsabs := %s (= epsrel %s x the smallest non-zero order magnitude of %s / sqrt(%d); --epsabs %s is the floor when larger)' % (epsabs_used, epsrel, src, NCH, a.epsabs))
    json.dump({'_note': 'disteval worker list: W local CPU workers; no nice wrapper (the workers inherit the fence this tool runs in)', 'cluster': [{'command': [sys.executable, '-m', 'pySecDecContrib', 'pysecdec_cpuworker'], 'count': W}]}, open(os.path.join(SD, 'cluster.json'), 'w'), indent=1)
    settings = {'stage': S, 'W': str(W), 'points': str(points), 'presamples': str(presamples), 'shifts': str(shifts), 'epsrel': str(epsrel), 'epsabs_per_chunk': str(epsabs_used), 'lattice_candidates': str(latcand), 'chunk_sectors': str(chunk_sectors), 'n_chunks': str(NCH), 'standin': '1' if dry else '0', 'stamp_utc': utc(),
                'integrator': ('disteval median-lattice rule (--lattice-candidates %d): the SECOND integrator = a different lattice construction on the same kernel library' % latcand) if latcand != 0 else 'disteval standard (tabulated Korobov generating vectors) rank-1 lattices',
                'point': ' '.join(point_args), 'package': {'disteval_dir': D, 'family': fam.name, 'pinned_to_the_fixture_of_record': pinned, 'source_sha256': man['source_sha256']}, 'PRODUCER': producer_block()}
    json.dump(settings, open(os.path.join(SD, 'SETTINGS_%s.json' % S), 'w'), indent=1)
    ndone = nfail = nrun = 0
    T0 = __import__('time').time()
    progress = os.path.join(WORK, 'progress.jsonl')
    LC = ['--lattice-candidates=%d' % latcand] if latcand != 0 else []
    for ci in range(NCH):
        CH = os.path.join(SD, 'chunks', 'chunk_%03d' % ci)
        if os.path.exists(os.path.join(CH, 'DONE')) and resume:
            ndone += 1; print('chunk %d DONE already (%s); skipped' % (ci, open(os.path.join(CH, 'DONE')).read().strip())); continue
        if a.max_chunks is not None and nrun >= int(a.max_chunks):
            print('STOPPED by --max-chunks %d after %d chunk(s) this run; the stage stays INCOMPLETE -- RESUME: the identical line with --resume runs the missing chunks' % (int(a.max_chunks), nrun))
            break
        if os.path.exists(os.path.join(CH, 'DONE')):
            os.replace(os.path.join(CH, 'DONE'), os.path.join(CH, 'DONE.superseded'))
        open(os.path.join(CH, 'start'), 'w').write(utc() + '\n')
        timefile = os.path.join(CH, 'time.txt')
        if dry:
            chk = "import pySecDec.disteval, pySecDecContrib, os, sys, json; d=sys.argv[1]; so=sys.argv[2]; ij=sys.argv[3]; assert os.path.exists(os.path.join(d,so)) and os.path.exists(os.path.join(d,'builtin.so')) and os.path.isdir(os.path.join(d,'coefficients')); j=json.load(open(os.path.join(d,ij))); assert j['kernels']; wb=os.path.join(os.path.dirname(pySecDecContrib.__file__),'bin','pysecdec_cpuworker'); assert os.access(wb, os.X_OK), wb; print('standin ok', pySecDec.__version__, len(j['kernels']), 'kernels', wb)"
            with open(os.path.join(CH, 'disteval.log'), 'w') as lg:
                RC = subprocess.run(['/usr/bin/time', '-v', '-o', timefile, sys.executable, '-c', chk, CH, fam.integral_so, fam.integral_json], stdout=lg, stderr=subprocess.STDOUT).returncode
            if RC == 0 and a.dry_result and os.path.exists(a.dry_result) and os.path.getsize(a.dry_result) > 0:
                shutil.copyfile(a.dry_result, os.path.join(CH, 'result.json'))
            elif RC == 0:
                RC = 91
        else:
            cmd = ['/usr/bin/time', '-v', '-o', timefile, sys.executable, '-m', 'pySecDec.disteval', os.path.join(CH, fam.sum_json), '--cluster', os.path.join(SD, 'cluster.json'), '--epsrel=%s' % epsrel, '--epsabs=%s' % epsabs_used, '--points=%s' % points, '--presamples=%s' % presamples, '--shifts=%d' % shifts] + LC + ['--format=json'] + point_args
            with open(os.path.join(CH, 'result.json'), 'w') as ro, open(os.path.join(CH, 'disteval.log'), 'w') as lg:
                RC = subprocess.run(cmd, stdout=ro, stderr=lg).returncode
        nrun += 1
        open(os.path.join(CH, 'rc'), 'w').write('rc=%d\n' % RC); open(os.path.join(CH, 'end'), 'w').write(utc() + '\n')
        OK = 0
        if RC == 0:
            try:
                r = json.load(open(os.path.join(CH, 'result.json'))); assert 'sums' in r and fam.name in r['sums']; OK = 1
            except Exception:
                OK = 0
        if OK:
            open(os.path.join(CH, 'DONE'), 'w').write(utc() + '\n'); ndone += 1
        else:
            nfail += 1
        wall, rss = time_fields(timefile)
        stats = ''
        for ln in open(os.path.join(CH, 'disteval.log'), errors='replace').read().splitlines():
            if fam.progress_literal in ln:
                stats = ln.split(']', 1)[-1].strip()
        with open(progress, 'a') as pf:
            pf.write(json.dumps({'utc': utc(), 'stage': S, 'chunk': ci, 'n_chunks': NCH, 'n_done': ndone, 'n_fail': nfail, 'rc': RC, 'wall': wall, 'maxrss_kb': rss, 'elapsed_stage_s': int(__import__('time').time() - T0), 'stats': stats}) + '\n')
        print('chunk %d/%d rc=%d ok=%d wall=%s maxrss_kb=%s done=%d fail=%d %s %s' % (ci, NCH, RC, OK, wall, rss, ndone, nfail, utc(), stats), flush=True)
    ns = types.SimpleNamespace(chunks_dir=os.path.join(SD, 'chunks'), stage=S, out=os.path.join(SD, 'RESULT_%s.json' % S), settings=os.path.join(SD, 'SETTINGS_%s.json' % S), name=fam.name)
    ARC = assemble_compare.assemble(ns)
    open(os.path.join(SD, 'assemble.rc'), 'w').write('rc=%d\n' % ARC)
    open(os.path.join(SD, 'stage.end'), 'w').write(utc() + '\n')
    if ARC == 0 and nfail == 0 and ndone == NCH:
        open(os.path.join(SD, 'STAGE_DONE'), 'w').write('%s chunks=%d done=%d\n' % (utc(), NCH, ndone)); open(os.path.join(SD, 'stage.rc'), 'w').write('rc=0\n')
        print('STAGE %s DONE %s: %d/%d chunks, assemble rc 0 -> %s' % (S, utc(), ndone, NCH, ns.out))
        return RC_PASS
    open(os.path.join(SD, 'stage.rc'), 'w').write('rc=1\n')
    print('STAGE %s INCOMPLETE %s: done=%d fail=%d of %d, assemble rc %d -- the identical line with --resume reruns the missing chunks (FAIL by name: incomplete)' % (S, utc(), ndone, nfail, NCH, ARC))
    return RC_FAIL


# ---------------------------------------------------------------- compare / gate / basis
def do_compare(a):
    if not a.results:
        raise Refusal(RC_USAGE, 'REFUSED: --stage compare needs --results A=RESULT_A.json B=... [C=...]')
    for spec in a.results:
        if '=' not in spec or not os.path.exists(spec.split('=', 1)[1]):
            raise Refusal(RC_USAGE, 'REFUSED: --results entry %r is not NAME=existing-file' % spec)
    ns = types.SimpleNamespace(results=a.results, out=a.out or 'COMPARE.json')
    return assemble_compare.compare(ns)


def record_check(out, record_path):
    """the worked-example leg: the produced acceptance check vs the record's digit class, per-order digits, the controls' verdicts, by name."""
    rec = json.load(open(record_path))
    checks = []
    checks.append(('digit class', out['digit_class']['min_digits_over_orders(conservative)'], rec['digit_class']['min_digits_over_orders(conservative)']))
    key = 'digits(min over the B-C pair and the B, C bars)'
    checks.append(('per-order digits', [r[key] for r in out['digit_class']['per_order']], [r[key] for r in rec['digit_class']['per_order']]))
    for k in ('A-B_sigma', 'A-C_sigma', 'B-C_sigma', 'A-B_agreed_digits', 'A-C_agreed_digits', 'B-C_agreed_digits', 'A_bar_digits', 'B_bar_digits', 'C_bar_digits'):
        checks.append((k, [r.get(k) for r in out['digit_class']['per_order']], [r.get(k) for r in rec['digit_class']['per_order']]))
    checks.append(('B-C fail lines', out['fail_lines(B-C sigma > 3 on an order)'], rec['fail_lines(B-C sigma > 3 on an order)']))
    checks.append(('A consistency verdict', out['A_consistency']['verdict'], rec['A_consistency']['verdict']))
    checks.append(('A max sigma', out['A_consistency']['max_sigma_A_vs_B_or_C'], rec['A_consistency']['max_sigma_A_vs_B_or_C']))
    checks.append(('planted raised', out['controls']['planted_fail']['raised_fail_line'], rec['controls']['planted_fail']['raised_fail_line']))
    checks.append(('planted sigma vs A', out['controls']['planted_fail']['sigma_vs_A'], rec['controls']['planted_fail']['sigma_vs_A']))
    checks.append(('planted sigma vs C', out['controls']['planted_fail']['sigma_vs_C'], rec['controls']['planted_fail']['sigma_vs_C']))
    checks.append(('foreign ok', out['controls']['foreign_ok'], rec['controls']['foreign_ok']))
    checks.append(('foreign max sigmas', [f['max_sigma'] for f in out['controls']['foreign']], [f['max_sigma'] for f in rec['controls']['foreign']]))
    checks.append(('foreign rows', [f['rows'] for f in out['controls']['foreign']], [f['rows'] for f in rec['controls']['foreign']]))
    checks.append(('gate', out['gate'], rec['gate']))
    ok = True
    for name, got, want in checks:
        same = json.dumps(got, sort_keys=True) == json.dumps(want, sort_keys=True)
        ok = ok and same
        print('  RECORD %-22s %s  (got %s)' % (name, 'REPRODUCED' if same else 'DIFFERS: record %s' % json.dumps(want)[:120], json.dumps(got)[:120]))
    print('RECORD %s: %d checks vs %s (%s)' % ('REPRODUCED' if ok else 'NOT REPRODUCED (FAIL by name)', len(checks), os.path.basename(record_path), sha256(record_path)[:16]))
    return ok


def do_gate(a, pinned, fam=None):
    fam = fam or FAMILY_OF_RECORD
    if a.fixtures:
        names = ['RESULT_A.json', 'RESULT_B.json', 'RESULT_C.json', 'smoke_build_host.json', 'smoke_second_host.json', 'LEAF_STAMP_record.txt', 'COMPARE_record.json', 'GATE_record.json']
        check_fixtures(names)
        results = {s: os.path.join(FIX, 'RESULT_%s.json' % s) for s in 'ABC'}
        stamp = gate_mod.read_stamp(os.path.join(FIX, 'LEAF_STAMP_record.txt'))
        rcs = {}
        for tok in (stamp.get('stages_rc') or '').split():
            if '=' in tok:
                k, v = tok.split('=', 1)
                if k in 'ABC':
                    rcs[k] = 'rc=%s' % v
        smokes = FOREIGN_SMOKES
        out_path = a.out or 'GATE_fixtures.json'
        out = gate_mod.run_gate(results, rcs, smokes, stamp, out_path, chunks_done=None, label='(the fixture of record)', package_pinned=True, compare_path=os.path.join(FIX, 'COMPARE_record.json'), family=fam.name)
        ok = record_check(out, os.path.join(FIX, 'GATE_record.json'))
        return RC_PASS if (out['gate'] == 'PASS' and ok) else RC_FAIL
    if not a.work:
        raise Refusal(RC_USAGE, 'REFUSED: --stage check (alias gate) needs --example/--fixtures or --work DIR')
    WORK = os.path.abspath(a.work)
    results, rcs, done = {}, {}, {}
    for s in 'ABC':
        p = os.path.join(WORK, 'stage_%s' % s, 'RESULT_%s.json' % s)
        if os.path.exists(p):
            results[s] = p
        rcp = os.path.join(WORK, 'stage_%s' % s, 'stage.rc')
        rcs[s] = open(rcp).read().strip() if os.path.exists(rcp) else None
        cd = os.path.join(WORK, 'stage_%s' % s, 'chunks')
        done[s] = len([d for d in os.listdir(cd) if d.startswith('chunk_') and os.path.exists(os.path.join(cd, d, 'DONE'))]) if os.path.isdir(cd) else 0
    smokes = [tuple(x.split('=', 1)) for x in (a.smoke or [])]
    if not smokes and pinned:
        smokes = FOREIGN_SMOKES
        print('FOREIGN control: the two fixture smokes (the package is pinned to the fixture of record, the same .so sha)')
    elif not smokes:
        print('FOREIGN control: none given and the package is unpinned -- the fixture smokes are not a control for a different library; the check will FAIL on foreign_ok unless --smoke NAME=smoke.json of this library is given')
    stamp = gate_mod.read_stamp(a.leaf_stamp)
    out_path = a.out or os.path.join(WORK, 'GATE.json')
    out = gate_mod.run_gate(results, rcs, smokes, stamp, out_path, chunks_done=done, label='(work %s)' % os.path.basename(WORK), package_pinned=pinned, compare_path=a.compare, family=fam.name)
    return RC_PASS if out['gate'] == 'PASS' else RC_FAIL


def do_planted_control(a):
    if a.work:
        results = {s: os.path.join(os.path.abspath(a.work), 'stage_%s' % s, 'RESULT_%s.json' % s) for s in 'ABC'}
    else:
        check_fixtures(['RESULT_A.json', 'RESULT_B.json', 'RESULT_C.json'])
        results = {s: os.path.join(FIX, 'RESULT_%s.json' % s) for s in 'ABC'}
    O = {s: gate_mod.orders_of(json.load(open(p))) for s, p in results.items() if os.path.exists(p)}
    pl = gate_mod.planted_fail_control(O)
    if pl is None:
        raise Refusal(RC_USAGE, 'REFUSED: the planted control needs RESULT_A and RESULT_B')
    print('CONTROL %s: sigma vs A = %s, sigma vs C = %s -> %s' % (pl['name'], pl['sigma_vs_A'], pl['sigma_vs_C'], 'FAIL line RAISED (the compare sees a 1e-3 shift; the control passes)' if pl['raised_fail_line'] else 'NO FAIL line (the compare is blind; the check FAILS by name)'))
    print('  rule: ' + pl['rule'])
    return RC_PASS if pl['raised_fail_line'] else RC_FAIL


def do_basis(a):
    bf = a.basis_file or os.path.join(FIX, 'BASIS_record.json')
    if not os.path.exists(bf):
        raise Refusal(RC_MISSING, 'REFUSED: basis file missing: %s' % bf)
    worker_B, driver_B, _ = memory_basis.read_basis_file(bf)
    e = memory_basis.basis_entry(int(a.basis), worker_B, driver_B, int(a.probe_W))
    print('BASIS W=%d %s: memory.max %d B (%.3f GiB) = page_up(1.3 x (driver %d B + %d x per_worker %d B)) [the nominal per-worker formula %d B leaves headroom %d B: driver fits %s]; raise-once %d B; basis %s (%s)' % (e['W'], e['label'], e['memory_max_bytes'], e['memory_max_GiB'], driver_B, e['W'], worker_B, e['formula_memory_max_B(per_worker x W x 1.3, page_up)'], e['formula_headroom_over_W_x_worker_B'], e['driver_fits_in_formula_headroom'], e['raise_once_bytes'], os.path.basename(bf), sha256(bf)[:16]))
    print('DECLARE on the fence: memory.max = %d B; cpu.max = %d 100000 (W x 100000); swap 0; the tool runs inside it and never creates it' % (e['memory_max_bytes'], e['W'] * 100000))
    return RC_PASS


def main():
    if '-h' in sys.argv[1:] or '--help' in sys.argv[1:]:
        print(__doc__)
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--disteval-dir', default=None, help='the compiled package disteval dir (loop_package + make disteval); pinned to the fixture of record unless --unpinned-package')
    ap.add_argument('--family', '--name', dest='family', default=DEFAULT_FAMILY, help='the package family NAME: the sum file NAME.json, NAME_integral.json / .so, the sums key NAME (default ndpent_top, the worked example\'s family; --name is the alias)')
    ap.add_argument('--pins', default=None, help='the PACKAGE_PINS file of the family (the same JSON shape; default PACKAGE_PINS.json beside the tool); missing = rc 4 by name')
    ap.add_argument('--emit-pins', default=None, metavar='OUT', help='write the PACKAGE_PINS shape of --family from the five disteval files of --disteval-dir (sha256 at emit) and exit')
    ap.add_argument('--generate-receipt', default=None, help='--emit-pins: the generate receipt whose sha256 is recorded (else null)')
    ap.add_argument('--unpinned-package', action='store_true', help='run a package whose files do not match the pins file (one loud line; results labelled package_pinned false)')
    ap.add_argument('--check-package', action='store_true', help='verify the package pins and exit 0 / 3')
    ap.add_argument('--stage', choices=['A', 'B', 'C', 'compare', 'check', 'gate'], default=None, help='A | B | C | compare | check (the acceptance check; gate is the older spelling, kept)')
    ap.add_argument('--point', default=POINT_OF_RECORD, help="'name=value ...' over the package's real parameters (default: the point of the fixture of record)")
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--points', default=None); ap.add_argument('--presamples', default='1e4'); ap.add_argument('--shifts', type=int, default=32)
    ap.add_argument('--epsrel', default=None); ap.add_argument('--epsabs', default='1e-12', help='the epsabs floor (per chunk when --epsabs-from applies)')
    ap.add_argument('--epsabs-from', default=None, help="a RESULT json whose per-order magnitudes set the per-chunk epsabs (default for B/C: --work/stage_A/RESULT_A.json when present)")
    ap.add_argument('--lattice-candidates', type=int, default=None); ap.add_argument('--chunk-sectors', type=int, default=None)
    ap.add_argument('--work', default=None, help='the work dir (stage_<S>/chunks/chunk_NNN/{result.json,time.txt,rc,DONE}, RESULT_<S>.json, STAGE_DONE; progress.jsonl)')
    ap.add_argument('--resume', action='store_true', help='skip DONE chunks and DONE stages')
    ap.add_argument('--max-chunks', type=int, default=None, help='evaluate at most N chunks this run (a pilot of the per-chunk wall); the stage stays INCOMPLETE and resumes')
    ap.add_argument('--dry', action='store_true', help='the stand-in form: the disteval call becomes an import/library-presence check under GNU time + a copy of --dry-result')
    ap.add_argument('--dry-result', default=os.path.join(FIX, 'smoke_build_host.json'))
    ap.add_argument('--results', nargs='*', default=None, help='compare: NAME=RESULT.json ...')
    ap.add_argument('--example', '--fixtures', dest='fixtures', action='store_true', help='check: over the worked example\'s records in ../examples/ndpent_top (must reproduce GATE_record.json); --fixtures is the alias')
    ap.add_argument('--smoke', action='append', default=None, help='check: NAME=smoke.json FOREIGN control(s)')
    ap.add_argument('--leaf-stamp', default=None, help='check: the leaf-at-end stamp (key=value lines) of the cgroup the stages ran in')
    ap.add_argument('--compare', default=None, help='check: the COMPARE json to cite')
    ap.add_argument('--out', default=None)
    ap.add_argument('--planted-fail-control', action='store_true', help='print the PLANTED-FAIL control line over the worked example\'s records (or --work) and exit by its verdict')
    ap.add_argument('--basis', type=int, default=None, help='W: print the memory figure to declare on the cgroup (memory.max), by the recorded formula')
    ap.add_argument('--basis-file', default=None); ap.add_argument('--probe-W', type=int, default=8)
    a = ap.parse_args()
    a.stage = STAGE_ALIASES.get(a.stage, a.stage)
    pins_path = os.path.abspath(a.pins) if a.pins else PACKAGE_PINS
    try:
        fam = validate_family(a.family)
        if a.basis is not None:
            return do_basis(a)
        if a.planted_fail_control:
            return do_planted_control(a)
        if a.emit_pins:
            if not a.disteval_dir or not os.path.isdir(a.disteval_dir):
                raise Refusal(RC_USAGE, 'REFUSED: --emit-pins needs --disteval-dir D (a directory)')
            return emit_pins(a.disteval_dir, fam, a.emit_pins, a.generate_receipt, a.point)
        if a.check_package:
            if not a.disteval_dir:
                raise Refusal(RC_USAGE, 'REFUSED: --check-package needs --disteval-dir')
            check_family_files(a.disteval_dir, fam)
            check_package(a.disteval_dir, a.unpinned_package, pins_path=pins_path, fam=fam)
            return RC_PASS
        if a.stage == 'compare':
            return do_compare(a)
        if a.stage == 'gate':
            pinned = True
            if a.disteval_dir:
                pinned = check_package(a.disteval_dir, a.unpinned_package, quiet=True, pins_path=pins_path, fam=fam)
            elif a.work and not a.fixtures:
                sj = [os.path.join(os.path.abspath(a.work), 'stage_%s' % s, 'SETTINGS_%s.json' % s) for s in 'ABC']
                sj = [p for p in sj if os.path.exists(p)]
                pinned = all(json.load(open(p)).get('package', {}).get('pinned_to_the_fixture_of_record', False) for p in sj) if sj else False
            return do_gate(a, pinned, fam)
        if a.stage in ('A', 'B', 'C'):
            if not a.disteval_dir:
                raise Refusal(RC_USAGE, 'REFUSED: --stage %s needs --disteval-dir' % a.stage)
            if not os.path.isdir(a.disteval_dir):
                raise Refusal(RC_USAGE, 'REFUSED: --disteval-dir %s is not a directory' % a.disteval_dir)
            if not a.work:
                raise Refusal(RC_USAGE, 'REFUSED: --stage %s needs --work DIR (the checkpoint root)' % a.stage)
            if a.dry and not (os.path.exists(a.dry_result) and os.path.getsize(a.dry_result) > 0):
                raise Refusal(RC_MISSING, 'REFUSED: --dry needs a canned result: %s missing' % a.dry_result)
            check_family_files(a.disteval_dir, fam)
            pinned = check_package(a.disteval_dir, a.unpinned_package, pins_path=pins_path, fam=fam)
            return run_stage(a, pinned, fam)
        ap.print_help()
        return RC_USAGE
    except Refusal as e:
        print(str(e))
        return e.rc


if __name__ == '__main__':
    sys.dont_write_bytecode = True
    raise SystemExit(main())
