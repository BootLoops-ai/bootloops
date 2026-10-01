"""blade.selftest — run the whole library against the reference artifact bank.

    ulimit -v 32505856; python3 -m blade.selftest [--phase0 DIR] [--scratch DIR]
    ulimit -v 32505856; python3 -m blade.selftest --full [--root DIR]

Quick sections (default, phase0 bank only):
  A. round-trip: byte-identical read->write on EVERY contract file in the
     phase0 tree (never in place -- rewrites go to scratch).
  B. pipeline: copy the phase0 INPUTS into scratch, pipeline.run_all() on the
     copy, and require rrres / rec_0_mono / rec_0_part byte-identical to the
     reference phase0 run (which independently passed the 18/18 kira_table
     gate).
  C. ratrec: plant two known 2-var rational functions, reconstruct through
     recmod+dynamicrr from a modular evaluator callback, compare exactly in
     sympy; then a mutation control (corrupt ONE evaluation -> BladeGateError,
     never a silently-returned clean or wrong answer).

Section I (raw-system fflow lane) is SELF-CONTAINED -- no banked artifacts:
typed lane-resolution refusals, then (where the redg1/fitrel/fflowcli
binaries are built; otherwise a named SKIP) blade_family `prepare` end-to-end
on a planted-truth toy sparse system, banked-probe values vs the planted
coefficients at a fresh prime, and a mutated-truth control.  It runs even
when BLADE_PORT_ROOT is unset, so a public clone gets a real battery.

--full REGRESSION BATTERY: additionally replays every phase's decisive gate
from the banked artifact trees (no new searches; everything deterministic):
  M. manifest: sha256 of every banked artifact in regression_manifest.json --
     FAILS LOUDLY listing missing/changed artifacts, then aborts the replays
     (never skips silently).
  D. gate-mutation trio: truncated eval / empty fit
     table (flag left 1) / dynamicrr state tamper, each via a wrapped binary
     on a fresh input copy -> the NAMED BladeGateError must fire.
  E. three-route database md5 agreement (synthetic == Kira-feed ==
     fflow-native RREF) vs the banked md5_all.txt.
  F. phase3: db self-generated search REPLAY (scheme+ansatz+escalation on the
     KIRA-lane database copy, seed 20260710), sch_* byte-identity vs the
     phase0 bank, exported relations content-compare vs the banked
     db_relations.json, annihilation spot-check at 2 fresh points.
  G. phase3.5: BTOracle probe at 12 points (6 db dense + 6 dbox-w1
     compositional) vs the banked Kira reduction tables; typed
     ReconstructionOverflow refusal on the banked w1 tree at 3 primes
     (fail-closed, no artifact emitted).
  H. phase4: semibl identity run on a fresh phase0 copy -> rec_*/rrres
     byte-identical to the banked reconstruction; banked symbolic dbox table
     spot-eval (4 closed targets x 32 masters x 2 points) vs the banked
     Kira table; w0 rows still explicit PENDING.

Prints a PASS/FAIL table; exits nonzero on any FAIL.
Battery contract: the reference artifact bank (BLADE_PORT_ROOT) is READ-ONLY
here -- every rerun happens on a scratch copy.
"""

from __future__ import annotations

import argparse
import filecmp
import hashlib
import json
import os
import random
import re
import shutil
import stat
import sys
import tempfile
import time
import traceback
from fractions import Fraction

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from blade import formats, pipeline, ratrec           # noqa: E402
    from blade.formats import BladeFormatError            # noqa: E402
    from blade.pipeline import BladeGateError, PipelineConfig  # noqa: E402
else:
    from . import formats, pipeline, ratrec
    from .formats import BladeFormatError
    from .pipeline import BladeGateError, PipelineConfig

# Root of the reference artifact bank the --full battery replays against
# (242 sha-pinned files; the bank is not distributed with this repository —
# missing artifacts are NAMED loudly and the replay aborts).
# env-only root: no machine-local default; the --full battery leg SKIPs with
# a named notice when the bank is absent.
BLADE_PORT_ROOT = os.environ.get("BLADE_PORT_ROOT", "")
PHASE0_DEFAULT = os.path.join(BLADE_PORT_ROOT, "phase0")
def _manifest_path(root):
    # the manifest lives beside the artifacts it pins (in the reference bank,
    # which is not distributed with this repository)
    return os.path.join(root, "regression_manifest.json")
FRESH_PRIME = (1 << 61) - 1          # never a fit/verify/oracle prime anywhere
BATTERY_SEED = 20260709

RESULTS = []  # (section, name, ok, detail)


def record(section, name, ok, detail=""):
    RESULTS.append((section, name, ok, detail))
    mark = "PASS" if ok else "FAIL <<<<<<<<<<<<<<<<"
    print(f"  [{mark}] {section} :: {name}" + (f" -- {detail}" if detail else ""))


# --------------------------------------------------------------------------
# A. round trips
# --------------------------------------------------------------------------

_SKIP_DIRS = {"logs", "verifier", "shim", "__pycache__"}
_SKIP_FILES = {"average_time"}


def contract_files(phase0):
    for root, dirs, files in os.walk(phase0):
        dirs[:] = [d for d in dirs if d not in _SKIP_DIRS]
        for f in sorted(files):
            if f in _SKIP_FILES or f.endswith((".md", ".py", ".log", ".txt.bak")):
                continue
            p = os.path.join(root, f)
            rel = os.path.relpath(p, phase0)
            if rel == "reconstructed_table.txt":   # report, not a contract file
                continue
            yield p, rel


def section_roundtrip(phase0, scratch):
    rt_scratch = os.path.join(scratch, "roundtrip")
    n_ok = n_fail = 0
    fails = []
    unclassified = []
    for path, rel in contract_files(phase0):
        spec = formats.classify(path)
        if spec is None:
            unclassified.append(rel)
            continue
        try:
            ok, detail = formats.round_trip(path, rt_scratch)
        except Exception as e:                     # loud, counted as FAIL
            ok, detail = False, f"exception: {e}"
        if ok:
            n_ok += 1
        else:
            n_fail += 1
            fails.append((rel, detail))
    record("A.roundtrip", f"{n_ok} files byte-identical", n_fail == 0,
           "; ".join(f"{r}: {d}" for r, d in fails) if fails else "")
    record("A.roundtrip", "no unclassified contract files", not unclassified,
           ", ".join(unclassified))


# --------------------------------------------------------------------------
# B. pipeline rerun on a copy
# --------------------------------------------------------------------------

INPUT_WORK_FILES = ("sch_g1", "sch_nint", "sch_intid", "multi_sch_intid")


def copy_inputs(phase0, tree):
    """Copy ONLY the pipeline inputs (same list as the phase0 verifier's
    setup_mutation.py): databases, kinematics, scheme, config/*/in_*,
    degrees.fflow.  system/evallist files are REGENERATED by run_all."""
    src_search = os.path.join(phase0, "search")
    dst_search = os.path.join(tree, "search")
    shutil.copytree(os.path.join(src_search, "database"),
                    os.path.join(dst_search, "database"))
    shutil.copytree(os.path.join(src_search, "kinematics"),
                    os.path.join(dst_search, "kinematics"))
    src_job = os.path.join(src_search, "job0")
    dst_job = os.path.join(dst_search, "job0")
    for wid in sorted(d for d in os.listdir(src_job) if d.isdigit()):
        sw, dw = os.path.join(src_job, wid), os.path.join(dst_job, wid)
        os.makedirs(dw, exist_ok=True)
        for f in INPUT_WORK_FILES:
            shutil.copy2(os.path.join(sw, f), os.path.join(dw, f))
        src_cfg = os.path.join(sw, "config")
        for k in sorted(os.listdir(src_cfg)):
            os.makedirs(os.path.join(dw, "config", k), exist_ok=True)
            for f in ("in_nvar", "in_var"):
                shutil.copy2(os.path.join(src_cfg, k, f),
                             os.path.join(dw, "config", k, f))
    os.makedirs(os.path.join(tree, "recmod"), exist_ok=True)
    shutil.copy2(os.path.join(phase0, "recmod", "degrees.fflow"),
                 os.path.join(tree, "recmod", "degrees.fflow"))


def section_pipeline(phase0, scratch, nthreads):
    tree = os.path.join(scratch, "pipeline_copy")
    if os.path.exists(tree):
        shutil.rmtree(tree)
    try:
        copy_inputs(phase0, tree)
        cfg = PipelineConfig(
            searchdir=os.path.join(tree, "search"),
            job="job0",
            recmod_dir=os.path.join(tree, "recmod"),
            datanames=["data_p0", "data_p1", "data_p2"],
            nworks=1, nthreads=nthreads)
        t0 = time.time()
        results = pipeline.run_all(cfg)
        wall = time.time() - t0
        record("B.pipeline", "run_all completed (all artifact gates)", True,
               f"{len(results)} stages, {wall:.1f}s wall")
    except (BladeGateError, BladeFormatError, OSError) as e:
        record("B.pipeline", "run_all completed (all artifact gates)", False, str(e))
        return
    for rel in ("recmod/rrres", "recmod/rec_0_mono", "recmod/rec_0_part"):
        a = open(os.path.join(phase0, rel), "rb").read()
        b = open(os.path.join(tree, rel.replace("recmod/", "recmod/", 1)), "rb").read()
        record("B.pipeline", f"{rel} byte-identical to phase0", a == b,
               f"{len(b)}B vs {len(a)}B")
    # cross-check: assembled functions from the rerun match phase0's rrres
    try:
        part = formats.RecPart.read(os.path.join(tree, "recmod", "rec_0_part"))
        mono = formats.RecMono.read(os.path.join(tree, "recmod", "rec_0_mono"), 2)
        rr = formats.RRRes.read(os.path.join(tree, "recmod", "rrres"))
        funcs = formats.assemble_rational_functions(part, mono, rr.values)
        record("B.pipeline", "assembly of rerun artifacts",
               len(funcs) == len(part.pairs), f"{len(funcs)} rational functions")
    except Exception as e:
        record("B.pipeline", "assembly of rerun artifacts", False, str(e))


# --------------------------------------------------------------------------
# C. ratrec plant + mutation control
# --------------------------------------------------------------------------

def section_ratrec(scratch, nthreads):
    import sympy as sp
    x, y = sp.symbols("x1 x2")
    f1_num = [(3, (2, 1)), (-7, (0, 1)), (Fraction(5, 2), (0, 0))]
    f1_den = [(1, (1, 0)), (-2, (0, 1)), (1, (0, 0))]
    f2_num = [(Fraction(1, 2), (1, 0)), (Fraction(-1, 3), (0, 0))]
    f2_den = [(1, (1, 1)), (4, (0, 0))]
    expected = [
        (3 * x**2 * y - 7 * y + sp.Rational(5, 2)) / (x - 2 * y + 1),
        (x / 2 - sp.Rational(1, 3)) / (x * y + 4),
    ]
    ev1 = ratrec.ratfun_evaluator(f1_num, f1_den)
    ev2 = ratrec.ratfun_evaluator(f2_num, f2_den)

    def evaluator(coords, prime):
        return [ev1(coords, prime), ev2(coords, prime)]

    # plant test: DEFAULT path (degree scan learned from the evaluator)
    wd = os.path.join(scratch, "ratrec_plant")
    try:
        funcs = ratrec.reconstruct(evaluator, nvars=2, maxdeg=4, nfuns=2,
                                   primes=(0, 1), workdir=wd,
                                   nthreads=nthreads, symbols=(x, y))
        for i, (got, want) in enumerate(zip(funcs, expected)):
            ok = sp.simplify(sp.cancel(got - want)) == 0
            record("C.ratrec", f"plant function {i+1} exact (scanned degrees)",
                   bool(ok), f"got {sp.sstr(got)}" if not ok else "")
    except (BladeGateError, BladeFormatError) as e:
        record("C.ratrec", "plant reconstruction", False, str(e))
        return

    exact_deg = formats.DegreesFile(2, 2, [
        formats.DegreeInfo(3, 1, [(2, 0, 1, 0), (1, 0, 1, 0)]),
        formats.DegreeInfo(1, 2, [(1, 0, 1, 0), (0, 0, 1, 0)]),
    ])

    def run_mutation(name, corrupt_call, degrees, allow_corrected):
        """Corrupt exactly ONE evaluator call.  The control excludes the
        failure mode 'WRONG answer without an error'.  DETECTED
        (BladeGateError) always passes; CORRECTED (returns and the result
        equals the clean plant, all downstream gates green) passes only where
        measured to be sound (degree scan: the corrupt point is absorbed as a
        spurious (t-t_bad) num/den common factor, gcd-stripped, and the
        4-point held-out validation certifies the result -- the analogue of
        the banked mutation variant C)."""
        state = {"calls": 0}

        def corrupt_evaluator(coords, prime):
            vals = evaluator(coords, prime)
            state["calls"] += 1
            if state["calls"] == corrupt_call and \
                    prime == formats.big_uint_primes()[0]:
                vals = [(vals[0] + 1) % prime, vals[1]]
            return vals

        wd2 = os.path.join(scratch, f"ratrec_mut_{name}")
        try:
            got = ratrec.reconstruct(corrupt_evaluator, nvars=2, maxdeg=4,
                                     nfuns=2, primes=(0, 1), workdir=wd2,
                                     nthreads=nthreads, degrees=degrees,
                                     symbols=(x, y))
            clean = all(sp.simplify(sp.cancel(g - w)) == 0
                        for g, w in zip(got, expected))
            if clean and allow_corrected:
                record("C.ratrec", f"mutation({name}) never wrong silently",
                       True, "CORRECTED: corruption absorbed+validated, "
                       "result equals clean plant")
            else:
                record("C.ratrec", f"mutation({name}) never wrong silently",
                       False,
                       "returned the CLEAN answer where detection was "
                       "required -- gate hole" if clean else
                       "returned a WRONG answer without raising -- gate hole")
        except BladeGateError as e:
            record("C.ratrec", f"mutation({name}) never wrong silently", True,
                   f"DETECTED: {str(e)[:100]}")

    # (a) corruption lands in the RECONSTRUCTION evals (exact degrees given,
    #     so every evaluator call feeds the eval files).  recmod consumes the
    #     poisoned point with NO redundancy -> must be DETECTED
    #     (recmod/dynamicrr gates); a clean return here is a gate hole.
    run_mutation("recon-eval", 3, exact_deg, allow_corrected=False)
    # (b) corruption lands in the DEGREE SCAN (default path): DETECTED or
    #     CORRECTED both sound (see docstring).
    run_mutation("degree-scan", 3, None, allow_corrected=True)


# --------------------------------------------------------------------------
# I. raw-system fflow lane (self-contained: toy planted-truth system; needs
#    no banked artifacts.  End-to-end legs run only where the redg1/fitrel/
#    fflowcli binaries are built -- otherwise they SKIP by name.)
# --------------------------------------------------------------------------

def _rawsys_truth(i, j, t, e, p):
    """Planted reduction coefficients of the toy raw system, mod p:
    T0 = ((t+2*eps)*M0 + 3*M1)/(t+1),  T1 = (eps*M0 + (t^2-5)*M1)/2."""
    if (i, j) == (0, 0):
        num, den = (t + 2 * e), (t + 1)
    elif (i, j) == (0, 1):
        num, den = 3, (t + 1)
    elif (i, j) == (1, 0):
        num, den = e, 2
    else:
        num, den = (t * t - 5), 2
    return num % p * pow(den % p, -1, p) % p


def _rawsys_fixture(dirpath):
    """Write a 2-parameter toy sparse system in the fflow JSON grammar:
    6 unknowns (2 auxiliaries, 2 targets, 2 masters -- most complex first,
    masters at the highest column ids), 4 homogeneous equations, needed =
    the 2 targets.  The targets reduce onto the masters with the known
    rational coefficients of _rawsys_truth; the auxiliaries are eliminated.
    Returns the family-spec path."""
    def rf(terms):        # numerator terms [(coeff_str, (e_t, e_eps)), ...]
        return [[len(terms), [[c, list(ex)] for c, ex in terms]],
                [1, [["1", [0, 0]]]]]

    eqs = [
        # (t+1)*T0 - (t+2*eps)*M0 - 3*M1 = 0
        [3, [2, 4, 5], [rf([("1", (1, 0)), ("1", (0, 0))]),
                        rf([("-1", (1, 0)), ("-2", (0, 1))]),
                        rf([("-3", (0, 0))])]],
        # 2*T1 - eps*M0 - (t^2-5)*M1 = 0
        [3, [3, 4, 5], [rf([("2", (0, 0))]),
                        rf([("-1", (0, 1))]),
                        rf([("-1", (2, 0)), ("5", (0, 0))])]],
        # A0 - T0 - t*M1 = 0        (auxiliary, eliminated)
        [3, [0, 2, 5], [rf([("1", (0, 0))]),
                        rf([("-1", (0, 0))]),
                        rf([("-1", (1, 0))])]],
        # A1 - eps*T1 - M0 = 0      (auxiliary, eliminated)
        [3, [1, 3, 4], [rf([("1", (0, 0))]),
                        rf([("-1", (0, 1))]),
                        rf([("-1", (0, 0))])]],
    ]
    data_path = os.path.join(dirpath, "sys_data.json")
    with open(data_path, "w") as f:
        json.dump([len(eqs), eqs], f)
    sys_path = os.path.join(dirpath, "sys.json")
    with open(sys_path, "w") as f:
        # [neqs, nvars, npars, n_needed, [needed cols], nfiles, [files]]
        json.dump([len(eqs), 6, 2, 2, [2, 3], 1, [data_path]], f)
    cols_path = os.path.join(dirpath, "columns.json")
    with open(cols_path, "w") as f:
        json.dump({"columns": [[3, 1, 1], [1, 3, 1], [2, 1, 1],
                               [1, 2, 1], [1, 1, 1], [1, 1, 0]]}, f)
    spec_path = os.path.join(dirpath, "spec.json")
    with open(spec_path, "w") as f:
        json.dump({"family": "toy",
                   "search_params": ["t", "eps"],
                   "targets": [[2, 1, 1], [1, 2, 1]],
                   "masters": [[1, 1, 1], [1, 1, 0]],
                   "fflow_system": sys_path,
                   "fflow_columns": cols_path,
                   "seed": 20260901}, f)
    return spec_path


def section_rawsys(scratch, nthreads):
    from . import blade_family as bf
    from .probe import BTOracle

    # typed lane-resolution refusals (pure python; no binaries involved)
    try:
        bf.resolve_lane({"fflow_system": __file__})
        record("I.rawsys", "fflow_system without fflow_columns refused",
               False, "no exception raised")
    except bf.SpecError as e:
        record("I.rawsys", "fflow_system without fflow_columns refused",
               "fflow_columns" in str(e), str(e)[:90])
    try:
        bf.resolve_lane({})
        record("I.rawsys", "no lane keys -> typed refusal names all lanes",
               False, "no exception raised")
    except bf.LaneUnavailable as e:
        record("I.rawsys", "no lane keys -> typed refusal names all lanes",
               "fflow_system" in str(e), str(e)[:90])

    need = ("redg1", "fitrel", "fflowcli")
    missing = [b for b in need
               if not os.path.exists(os.path.join(pipeline.DEFAULT_BIN_DIR,
                                                  b))]
    if missing:
        print(f"  [SKIP] I.rawsys end-to-end: binaries {missing} not built "
              f"under BLADE_BIN_DIR={pipeline.DEFAULT_BIN_DIR} -- build the "
              f"blade fork to arm this leg", flush=True)
        return

    root = os.path.join(scratch, "rawsys")
    shutil.rmtree(root, ignore_errors=True)
    os.makedirs(root)
    spec_path = _rawsys_fixture(root)
    wd = os.path.join(root, "work")
    t0 = time.time()
    rc = bf.main(["prepare", spec_path, "--workdir", wd,
                  "--nthreads", str(nthreads),
                  "--max-round-seconds", "300"])
    record("I.rawsys", "turnkey prepare on the raw fflow system exits 0 "
           "(incl. held-out gate)", rc == 0,
           f"rc={rc}, {time.time() - t0:.1f}s wall")
    if rc != 0:
        return
    try:
        ora = BTOracle(os.path.join(wd, "family", "toy"))
    except BladeGateError as e:
        record("I.rawsys", "banked family loads through BTOracle", False,
               str(e)[:150])
        return
    p = FRESH_PRIME                     # exact backend: never a fit prime
    pts = _points("rawsys-probe", 2, 2, p)
    res = ora.reduce_many([(pt[:1], pt[1]) for pt in pts], p)
    n_bad = n_chk = 0
    for i, pt in enumerate(pts):
        trow_by_t = res[i]
        for t in range(2):
            trow = trow_by_t.get(ora.form.label(t), {})
            for m in range(2):
                got = trow.get(ora.form.label(2 + m), 0)
                want = _rawsys_truth(t, m, pt[0], pt[1], p)
                n_chk += 1
                n_bad += (got != want)
    record("I.rawsys", "probe == planted truth (2 pts x 4 entries, "
           "fresh prime)", n_bad == 0, f"{n_chk - n_bad}/{n_chk}")
    # synthetic-truth control: a mutated truth value MUST mismatch
    pt = pts[0]
    got = res[0].get(ora.form.label(0), {}).get(ora.form.label(2), 0)
    mutated = (_rawsys_truth(0, 0, pt[0], pt[1], p) + 1) % p
    record("I.rawsys", "control: mutated truth detected as mismatch",
           got != mutated, "compare has no teeth" if got == mutated else "")


# ==========================================================================
# --full regression battery (sections M, D, E, F, G, H)
# ==========================================================================

def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def _points(tag, nvars, npts, prime):
    rng = random.Random(f"{BATTERY_SEED}:{tag}")
    return [tuple(rng.randrange(2, prime - 1) for _ in range(nvars))
            for _ in range(npts)]


# ---- M. manifest -----------------------------------------------------------

def section_manifest(root):
    """sha256-verify every banked artifact.  Returns True iff the bank is
    intact; on failure every missing/changed artifact gets its own loud FAIL
    row (never a silent skip)."""
    mpath = _manifest_path(root)
    if not os.path.exists(mpath):
        record("M.manifest", "regression_manifest.json present", False,
               f"{mpath} missing -- regenerate with the bank-side "
               f"build_manifest.py")
        return False
    man = json.load(open(mpath))
    missing, changed = [], []
    for rel, want in sorted(man["artifacts"].items()):
        full = os.path.join(root, rel)
        if not os.path.exists(full):
            missing.append(rel)
        elif _sha256(full) != want:
            changed.append(rel)
    ok = not missing and not changed
    record("M.manifest", f"{len(man['artifacts'])} banked artifacts sha256-ok",
           ok, "" if ok else f"{len(missing)} missing, {len(changed)} changed")
    for rel in missing:
        record("M.manifest", f"MISSING banked artifact: {rel}", False, rel)
    for rel in changed:
        record("M.manifest", f"CHANGED banked artifact: {rel}", False, rel)
    return ok


# ---- D. gate-mutation trio (replayed) ----------------------------------------

_TRIO_BINS = ("redg1", "fitrel", "ssolve", "recmod", "dynamicrr", "dumppoints")


def _trio_wrappers(real_bin):
    return {
        "trunc-eval": ("ssolve", "ssolve", "bytes != expected", """#!/bin/bash
"%s/ssolve" "$@"
rc=$?
for out; do :; done
truncate -s 4096 "$out"
exit $rc
""" % real_bin),
        "empty-fit": ("fitrel", "fitrel", "rows < tmp_nsol", """#!/bin/bash
"%s/fitrel" "$@"
rc=$?
searchdir="$1"; db="$2"; job="$3"
: > "$searchdir/$job/0/fit/$db/3"
exit $rc
""" % real_bin),
        "state-tamper": ("dynamicrr", "dynamicrr", "!= 2", """#!/bin/bash
"%s/dynamicrr" "$@"
rc=$?
printf '1' > state
exit $rc
""" % real_bin),
    }


def section_gatetrio(phase0, scratch, nthreads):
    real_bin = pipeline.DEFAULT_BIN_DIR
    for variant, (wname, err_stage, err_frag, script) in \
            _trio_wrappers(real_bin).items():
        vroot = os.path.join(scratch, "gatetrio", variant)
        shutil.rmtree(vroot, ignore_errors=True)
        os.makedirs(vroot)
        tree = os.path.join(vroot, "tree")
        copy_inputs(phase0, tree)
        bd = os.path.join(vroot, "bin")
        os.makedirs(bd)
        for b in _TRIO_BINS:
            if b == wname:
                p = os.path.join(bd, b)
                open(p, "w").write(script)
                os.chmod(p, os.stat(p).st_mode | stat.S_IXUSR | stat.S_IXGRP)
            else:
                src = os.path.join(real_bin, b)
                if os.path.exists(src):
                    os.symlink(src, os.path.join(bd, b))
        cfg = PipelineConfig(searchdir=os.path.join(tree, "search"), job="job0",
                             recmod_dir=os.path.join(tree, "recmod"),
                             datanames=["data_p0", "data_p1", "data_p2"],
                             nworks=1, nthreads=nthreads, bin_dir=bd)
        try:
            pipeline.run_all(cfg)
            record("D.gatetrio", f"{variant}: named gate fires", False,
                   "SILENT SUCCESS -- gate hole")
        except BladeGateError as e:
            named = err_stage in str(e) and err_frag in str(e)
            record("D.gatetrio", f"{variant}: named gate fires", named,
                   f"BladeGateError: {str(e)[:120]}" if named else
                   f"WRONG gate: {str(e)[:160]}")
        except Exception as e:                                 # noqa: BLE001
            record("D.gatetrio", f"{variant}: named gate fires", False,
                   f"non-gate {type(e).__name__}: {e}")


# ---- E. three-route database md5 agreement ----------------------------------

def section_phase2(root):
    lanes = {
        "phase0": os.path.join(root, "phase0", "search", "database"),
        "kira": os.path.join(root, "phase2", "kirafeed", "database"),
        "fflow": os.path.join(root, "phase2", "fflownative", "out", "database"),
    }
    banked = {}
    for ln in open(os.path.join(root, "phase2", "verify", "md5_all.txt")):
        parts = ln.split()
        if len(parts) == 3:
            banked[(parts[0], parts[2])] = parts[1]
    n_ok = n_bad = 0
    bad = []
    sets = {lane: {} for lane in lanes}
    for (lane, rel), want in sorted(banked.items()):
        full = os.path.join(lanes[lane], rel[2:] if rel.startswith("./") else rel)
        if not os.path.exists(full):
            n_bad += 1
            bad.append(f"{lane}:{rel} MISSING")
            continue
        got = hashlib.md5(open(full, "rb").read()).hexdigest()
        sets[lane][rel] = got
        if got == want:
            n_ok += 1
        else:
            n_bad += 1
            bad.append(f"{lane}:{rel} md5 {got[:8]} != banked {want[:8]}")
    record("E.md5", f"{n_ok} route md5s == banked md5_all.txt",
           n_bad == 0, "; ".join(bad[:4]))
    three_way = (sets["phase0"] == sets["kira"] == sets["fflow"]
                 and len(sets["phase0"]) == 21)
    record("E.md5", "three routes byte-identical (21 files)",
           three_way, f"{len(sets['phase0'])} files/lane")


# ---- shared Kira-table oracles (banked .m exports, our Kira 3.1) ------------

def _kc_db(root, bs, man):
    return bs.KiraCoeffs(
        os.path.join(root, "phase2", "kirafeed", "kira", "results", "db",
                     "kira_targets_ext.m"), "db",
        man["targets"], man["masters"],
        env_fn=lambda coords, p: {
            "t": bs._GF(coords[0], p),
            "d": bs._GF(4, p) - bs._GF(2, p) * bs._GF(coords[1], p)})


def _kc_dbox(root, bs, man2):
    return bs.KiraCoeffs(
        os.path.join(root, "phase2", "verify", "dbox", "kira", "results",
                     "dbox", "kira_targets_dbox.m"), "dbox",
        man2["targets"], man2["masters"],
        env_fn=lambda coords, p: {
            "msq": bs._GF(coords[0], p), "t": bs._GF(coords[1], p),
            "d": bs._GF(4, p) - bs._GF(2, p) * bs._GF(coords[2], p)})


# ---- F. phase3 db self-generated search replay ------------------------------

def section_phase3(root, scratch, nthreads):
    from . import search as bs
    from .ansatz import FamilyAnsatz, integral_dimension
    from .scheme import RedTable, single_block_scheme, write_scheme

    p0 = os.path.join(root, "phase0")
    man = json.load(open(os.path.join(p0, "manifest.json")))
    all_ints = [tuple(t) for t in man["targets"]] + \
               [tuple(m) for m in man["masters"]]
    nm = len(man["masters"])
    labels = ["db[" + ",".join(map(str, t)) + "]" for t in all_ints]

    tree = os.path.join(scratch, "phase3_db")
    shutil.rmtree(tree, ignore_errors=True)
    search_dir = os.path.join(tree, "search")
    os.makedirs(os.path.join(search_dir, "database"))
    for pid in range(3):
        shutil.copytree(
            os.path.join(root, "phase2", "kirafeed", "database", f"data_p{pid}"),
            os.path.join(search_dir, "database", f"data_p{pid}"))

    # scheme self-generation from the numeric database (catches scheme.py rot)
    try:
        rt, rt2 = RedTable.pair(os.path.join(search_dir, "database", "data_p0"))
        res = single_block_scheme(all_ints, nm, rt, rt2, order="given")
        write_scheme(res, os.path.join(search_dir, "job0"))
    except Exception as e:                                     # noqa: BLE001
        record("F.phase3", "scheme self-generation", False,
               f"{type(e).__name__}: {e}")
        return
    record("F.phase3", "single-block scheme g1 == [0,1,2,3]",
           len(res.blocks) == 1 and res.blocks[0].g1 == [0, 1, 2, 3],
           f"g1={res.blocks[0].g1}")
    n_id = 0
    diffs = []
    sch_files = ["job0/0/sch_g1", "job0/0/sch_nint", "job0/0/sch_intid",
                 "job0/0/multi_sch_intid"]
    for rel in sch_files:
        same = filecmp.cmp(os.path.join(p0, "search", rel),
                           os.path.join(search_dir, rel), shallow=False)
        n_id += same
        if not same:
            diffs.append(rel)
    record("F.phase3", "sch_* byte-identical to phase0 bank (4 files)",
           n_id == 4, ", ".join(diffs))

    # escalation replay: identical recipe + grower seed as the banked run
    kc = _kc_db(root, bs, man)
    fam = FamilyAnsatz(search_params=("t", "eps"), var_groups=[["t", "eps"]],
                       var_weights=[[1, 1]], intmode=["dimension"], cut=[],
                       nmax=20)
    cfg = bs.SearchFamily(
        searchdir=search_dir, job="job0", fam=fam,
        work_weights=[[integral_dimension(x) for x in all_ints]],
        datanames=["data_p0", "data_p1", "data_p2"], nthreads=nthreads,
        ansatz_size=7, level_base=0, level_step=1, data_size=30,
        datasize_policy="reactive", log=lambda *_: None)
    cfg.grower = bs.make_kira_grower(cfg, kc, seed=20260710)
    try:
        out = bs.escalate_work(cfg, 0)
        record("F.phase3", "db escalation closes (flag==1)", out.closed,
               out.stop_reason if not out.closed
               else f"{len(out.rounds)} rounds")
        if not out.closed:
            return
        for pid in range(3):
            bs.run_fitrel_works(cfg, cfg.datanames[pid], 0, 1)
        doc = bs.export_relations(
            cfg, [0], labels, os.path.join(tree, "db_relations.json"),
            extra_meta={"family": "db",
                        "provenance": "self-generated scheme+ansatz"
                                      " on phase2 KIRA-lane database"})
    except BladeGateError as e:
        record("F.phase3", "search replay + export", False, str(e)[:200])
        return
    banked = json.load(open(os.path.join(root, "phase3", "search",
                                         "db_relations.json")))
    same = (doc["relations"] == banked["relations"]
            and doc["n_relations"] == banked["n_relations"] == 4
            and doc["primes_fit"] == banked["primes_fit"]
            and doc["primes_verify"] == banked["primes_verify"])
    record("F.phase3", "exported relations CONTENT == banked db_relations",
           same, f"{doc['n_relations']} relations, "
                 f"{sum(r['n_terms'] for r in doc['relations'])} terms")

    # annihilation spot-check: substitute banked Kira reduction rows at 2
    # fresh points mod a fresh prime -- residual per master must vanish
    p = FRESH_PRIME
    n_bad = n_chk = 0
    for pt in _points("phase3-annih", 2, 2, p):
        rows = {}
        for rel in doc["relations"]:
            resid = [0] * nm
            for term in rel["terms"]:
                gid = term["global_id"]
                c = Fraction(term["coeff"])
                v = c.numerator % p * pow(c.denominator % p, -1, p) % p
                for pn, e in term["monomial"].items():
                    x = pt[0] if pn == "t" else pt[1]
                    v = v * pow(x, e, p) % p
                if gid not in rows:
                    rows[gid] = kc.red_row(gid, list(pt), p, nm)
                for j in range(nm):
                    if rows[gid][j]:
                        resid[j] = (resid[j] + v * rows[gid][j]) % p
            n_chk += 1
            n_bad += any(resid)
    record("F.phase3", "annihilation vs Kira rows (2 fresh pts, fresh prime)",
           n_bad == 0, f"{n_chk - n_bad}/{n_chk} relation-residuals zero")


# ---- G. phase3.5 probe + typed export refusal -------------------------------

def section_phase35(root, scratch, nthreads):
    from . import search as bs
    from .ansatz import FamilyAnsatz, integral_dimension
    from .formats import Scheme
    from .probe import BTOracle

    p = FRESH_PRIME
    p0man = json.load(open(os.path.join(root, "phase0", "manifest.json")))
    kc_db = _kc_db(root, bs, p0man)
    fam_root = os.path.join(root, "phase35", "firstblood", "families")

    # db: dense 4 targets x 8 masters at 6 fresh points
    try:
        ora = BTOracle(os.path.join(fam_root, "db"))
        pts = _points("probe-db", 2, 6, p)
        res = ora.reduce_many([(pt[:1], pt[1]) for pt in pts], p)
        n_bad = n_chk = 0
        for i, pt in enumerate(pts):
            for t in range(4):
                trow = res[i].get(ora.form.label(t), {})
                for m in range(8):
                    want = kc_db.value(t, m, list(pt), p)
                    got = trow.get(ora.form.label(4 + m), 0)
                    n_chk += 1
                    n_bad += (got != want)
        record("G.phase35", "db probe == Kira table (6 pts x 32 entries)",
               n_bad == 0, f"{n_chk - n_bad}/{n_chk}")
    except BladeGateError as e:
        record("G.phase35", "db probe == Kira table", False, str(e)[:200])

    # dbox w1: compositional check at 6 fresh points
    man2 = json.load(open(os.path.join(root, "phase2", "verify", "dbox",
                                       "manifest_dbox.json")))
    kc_dbox = _kc_dbox(root, bs, man2)
    try:
        ora1 = BTOracle(os.path.join(fam_root, "dboxw1"))
        basis = ora1.form.independent_gids          # 32 masters then extras
        closed = list(ora1.form.closed_gids)
        pts = _points("probe-dboxw1", 3, 6, p)
        res = ora1.reduce_many([(pt[:2], pt[2]) for pt in pts], p,
                               want_gids=closed)
        n_bad = n_chk = 0
        for i, pt in enumerate(pts):
            kvec = {b: kc_dbox.red_row(b, list(pt), p, 32) for b in basis}
            for g in closed:
                want = kc_dbox.red_row(g, list(pt), p, 32)
                grow = res[i].get(ora1.form.label(g), {})
                got = [0] * 32
                for b in basis:
                    c = grow.get(ora1.form.label(b), 0)
                    if c:
                        kb = kvec[b]
                        for j in range(32):
                            if kb[j]:
                                got[j] = (got[j] + c * kb[j]) % p
                n_chk += 32
                n_bad += sum(int(a != b) for a, b in zip(got, want))
        record("G.phase35", "dbox w1 probe o Kira == Kira (6 pts, "
               "compositional)", n_bad == 0, f"{n_chk - n_bad}/{n_chk}")
    except BladeGateError as e:
        record("G.phase35", "dbox w1 probe compositional", False, str(e)[:200])

    # typed refusal: banked w1 tree at 3 primes CANNOT represent the 490-bit
    # heights -> ReconstructionOverflow, fail-closed, nothing emitted
    dbox = os.path.join(root, "phase35", "dbox_run")
    all_ints = [tuple(t) for t in man2["targets"]] + \
               [tuple(m) for m in man2["masters"]]
    labels = ["dbox[" + ",".join(map(str, t)) + "]" for t in all_ints]
    dims = [integral_dimension(x) for x in all_ints]
    sch1 = Scheme.read(os.path.join(dbox, "search", "job0", "1"))
    fam = FamilyAnsatz(search_params=tuple(man2["params"]),
                       var_groups=[list(man2["params"])],
                       var_weights=[[1, 1, 1]], intmode=["dimension"],
                       cut=[], nmax=28)
    cfg = bs.SearchFamily(
        searchdir=os.path.join(dbox, "search"), job="job0", fam=fam,
        work_weights=[[0], [dims[g] for g in sch1.intid], [0]],
        datanames=["data_p0", "data_p1", "data_p2"], log=lambda *_: None)
    out_refuse = os.path.join(scratch, "w1_relations_SHOULD_NOT_EXIST.json")
    if os.path.exists(out_refuse):
        os.unlink(out_refuse)
    try:
        bs.export_relations(cfg, [1], labels, out_refuse)
        record("G.phase35", "w1 export at 3 primes REFUSES (typed)", False,
               "export SUCCEEDED where overflow refusal was required")
    except bs.ReconstructionOverflow as e:
        typed = isinstance(e, BladeGateError)
        record("G.phase35", "w1 export at 3 primes REFUSES (typed)",
               typed and not os.path.exists(out_refuse),
               f"ReconstructionOverflow, no artifact emitted"
               if not os.path.exists(out_refuse) else
               "REFUSED but artifact was written -- fail-closed hole")
    except Exception as e:                                     # noqa: BLE001
        record("G.phase35", "w1 export at 3 primes REFUSES (typed)", False,
               f"WRONG exception {type(e).__name__}: {e}")


# ---- H. phase4 semibl db regression + banked dbox table spot-eval -----------

def section_phase4(root, scratch, nthreads):
    from . import search as bs
    from .formats import DegreesFile
    from .semibl import Orderings, SemiBLJob

    p0 = os.path.join(root, "phase0")
    man = json.load(open(os.path.join(p0, "manifest.json")))
    nt = man["nint"] - man["nmaster"]
    elements = [(i, nt + j) for i, j in man["elements"]]
    labels = [f"T{i}" for i in range(nt)] + \
             [f"M{j}" for j in range(man["nmaster"])]

    work = os.path.join(scratch, "phase4_semibl")
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work)
    shutil.copytree(os.path.join(p0, "search"), os.path.join(work, "search"))
    for pid in range(3):
        for f in (f"points_p{pid}.fflow", f"evallist_p{pid}.txt"):
            fp = os.path.join(work, "search", "job0", f)
            if os.path.exists(fp):
                os.unlink(fp)
    shutil.rmtree(os.path.join(work, "search", "job0", "eval"),
                  ignore_errors=True)
    deg0 = DegreesFile.read(os.path.join(p0, "recmod", "degrees.fflow"))
    try:
        job = SemiBLJob(
            searchdir=os.path.join(work, "search"), job="job0",
            key_dir=os.path.join(work, "key"),
            datanames=["data_p0", "data_p1", "data_p2"],
            orderings=Orderings.identity(["t", "eps"]),
            system_targets=[0, 1, 2, 3], elements=elements, degrees=deg0,
            labels=labels, nthreads=nthreads, log=lambda *_: None)
        job.run_reconstruction(min_primes=2, max_primes=3)
        record("H.phase4", "semibl grow-on-demand consumed 2 primes",
               job.meta["primes_consumed"] == 2,
               f"consumed={job.meta['primes_consumed']}")
        job.stage_prime(2)
        job.run_recmod(2)
    except BladeGateError as e:
        record("H.phase4", "semibl db identity run", False, str(e)[:200])
        return
    n_id = 0
    diffs = []
    names = [f"rec_{pid}_{s}" for pid in range(3)
             for s in ("coeff", "mono", "part")] + ["rrres"]
    for name in names:
        same = filecmp.cmp(os.path.join(work, "key", "recmod", name),
                           os.path.join(p0, "recmod", name), shallow=False)
        n_id += same
        if not same:
            diffs.append(name)
    record("H.phase4", "semibl rec_*/rrres byte-identical to banked (10)",
           n_id == len(names), ", ".join(diffs))

    # banked symbolic dbox table: spot-eval 4 closed targets vs Kira table
    man2 = json.load(open(os.path.join(root, "phase2", "verify", "dbox",
                                       "manifest_dbox.json")))
    kc_dbox = _kc_dbox(root, bs, man2)
    tab = json.load(open(os.path.join(root, "phase4", "semibl", "dbox",
                                      "dbox_table.json")))
    record("H.phase4", "dbox table: w0 rows explicit PENDING, w1/w2 closed",
           sorted(tab["rows"]) == ["4", "5", "6", "7"]
           and sorted(tab["pending"]) == ["0", "1", "2", "3"],
           f"rows={sorted(tab['rows'])} pending={sorted(tab['pending'])}")
    p = FRESH_PRIME
    G = bs._GF
    compiled = {}
    for t, row in tab["rows"].items():
        for m, expr in row.items():
            src = re.sub(r"(?<!\*\*)\b(\d+)\b", r"G(\1)",
                         expr.replace("^", "**"))
            compiled[(int(t), int(m))] = compile(src, f"<tab{t},{m}>", "eval")
    n_bad = n_chk = 0
    for pt in _points("dbox-table", 3, 2, p):     # (msq, t, eps)
        env = {"G": lambda v: G(v, p), "msq": G(pt[0], p), "t": G(pt[1], p),
               "eps": G(pt[2], p), "__builtins__": {}}
        for tgid in (4, 5, 6, 7):
            for m in range(32):
                want = kc_dbox.value(tgid, m, list(pt), p)
                code = compiled.get((tgid, 8 + m))
                got = eval(code, dict(env)).v if code is not None else 0
                n_chk += 1
                n_bad += (got != want)
    record("H.phase4", "dbox table spot-eval == Kira (4 tgts x 32 x 2 pts)",
           n_bad == 0, f"{n_chk - n_bad}/{n_chk}")


# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase0", default=None,
                    help="phase0 bank (default: <root>/phase0)")
    ap.add_argument("--root", default=BLADE_PORT_ROOT,
                    help="reference artifact bank root (mirrors allowed; the "
                         "manifest is verified against it)")
    ap.add_argument("--full", action="store_true",
                    help="regression battery: replay every phase's decisive "
                         "gate from the banked artifacts (manifest-pinned)")
    ap.add_argument("--scratch", default=None,
                    help="scratch dir (default: mkdtemp under $TMPDIR)")
    ap.add_argument("--nthreads", type=int, default=6)
    ap.add_argument("--keep", action="store_true", help="keep scratch tree")
    args = ap.parse_args()
    if not args.root:
        print("blade selftest: BLADE_PORT_ROOT unset and no --root given "
              "— the reference-bank legs are SKIPPED (the bank is not "
              "distributed with this repository); running the self-contained "
              "raw-system section only", flush=True)
    phase0 = args.phase0 or (os.path.join(args.root, "phase0")
                             if args.root else "")

    scratch = args.scratch or tempfile.mkdtemp(prefix="blade_selftest_")
    os.makedirs(scratch, exist_ok=True)
    print(f"blade selftest{' --full' if args.full else ''}: "
          f"root={args.root or '(none)'}\n  phase0={phase0 or '(none)'}\n"
          f"  scratch={scratch}")

    t0 = time.time()
    try:
        section_rawsys(scratch, args.nthreads)
    except Exception as e:                                     # noqa: BLE001
        record("I.rawsys", "section completed", False,
               f"UNCAUGHT {type(e).__name__}: {e} -- "
               + traceback.format_exc(limit=3).replace("\n", " | "))
    bank_ok = bool(args.root)
    if args.root:
        if args.full:
            bank_ok = section_manifest(args.root)
        if bank_ok:
            section_roundtrip(phase0, scratch)
            section_pipeline(phase0, scratch, args.nthreads)
            section_ratrec(scratch, args.nthreads)
            if args.full:
                for name, fn in (("D.gatetrio", section_gatetrio),
                                 ("E.md5", section_phase2),
                                 ("F.phase3", section_phase3),
                                 ("G.phase35", section_phase35),
                                 ("H.phase4", section_phase4)):
                    try:
                        if fn is section_gatetrio:
                            fn(phase0, scratch, args.nthreads)
                        elif fn is section_phase2:
                            fn(args.root)
                        else:
                            fn(args.root, scratch, args.nthreads)
                    except Exception as e:                     # noqa: BLE001
                        record(name, "section completed", False,
                               f"UNCAUGHT {type(e).__name__}: {e} -- "
                               + traceback.format_exc(limit=3).replace("\n",
                                                                      " | "))
        else:
            record("M.manifest", "replays SKIPPED: bank not intact "
                   "(fix/rebuild the manifest first)", False,
                   "battery refuses to replay against a corrupt bank")
    wall = time.time() - t0

    nfail = sum(1 for _s, _n, ok, _d in RESULTS if not ok)
    print("\n" + "=" * 72)
    print(f"{'section':<14}{'test':<44}{'result'}")
    print("-" * 72)
    for s, n, ok, d in RESULTS:
        print(f"{s:<14}{n[:43]:<44}{'PASS' if ok else 'FAIL  <<<<<<'}")
    print("-" * 72)
    print(f"{len(RESULTS) - nfail}/{len(RESULTS)} passed in {wall:.1f}s wall"
          + ("" if nfail == 0 else f"  --  {nfail} FAILURE(S)"))
    if not args.keep and args.scratch is None and nfail == 0:
        shutil.rmtree(scratch, ignore_errors=True)
    sys.exit(1 if nfail else 0)


if __name__ == "__main__":
    main()
