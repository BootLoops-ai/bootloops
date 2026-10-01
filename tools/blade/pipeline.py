"""blade.pipeline — gated stage runner for the Wolfram-free Blade C pipeline.

Stages: redg1 -> fitrel (per prime) -> dumppoints -> ssolve -> recmod (per
prime) -> dynamicrr.  Every stage = run + ARTIFACT GATE.

MEASURED GATE DISCIPLINE:
  * rc is RECORDED but NEVER gates: redg1 always exits 0; fitrel writes
    flag=1 over an EMPTY fit nullspace (mutation variant E); ssolve exits 0
    writing a wrong-size eval; recmod prints "failed" and exits 0; fflowC has
    ~15 rc-0 failure paths.
  * Gates that DO fire: fit-table DIMENSIONS, eval BYTE SIZES, dynamicrr
    state==2 + rrres coefficient count, and the caller's final held-out
    compare.  flag==0 is trusted as a FAILURE signal only (self-check fired).

All failures raise BladeGateError naming the offending artifact.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

from . import formats
from .formats import (BladeFormatError, Database, DegreesFile, EvalFile,
                      EvalList, FitTable, Kinematics, PointsFile, RecCoeff,
                      RecMono, RecPart, RRRes, Scheme, SystemFile, TmpTemplate,
                      big_uint_primes, flags_size, read_ints)

# Defaults: the blade fork lives in the sibling repository
# blade (build it; its CMake install step puts the binaries in
# <fork>/bin). BLADE_BIN_DIR names the bin dir; when unset, the default is a
# checkout of blade beside this repository:
# <repo root>/../blade/bin.
import os as _os
DEFAULT_BIN_DIR = _os.environ.get("BLADE_BIN_DIR") or _os.path.normpath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)),
                  "..", "..", "..", "blade", "bin"))
# Alternate location of the dumppoints binary (env BLADE_PHASE0_DUMPPOINTS);
# the bin-dir copy is preferred whenever it exists.
PHASE0_DUMPPOINTS = _os.environ.get(
    "BLADE_PHASE0_DUMPPOINTS",
    _os.path.join(DEFAULT_BIN_DIR, "dumppoints"))


class BladeGateError(Exception):
    """An artifact gate failed. Message names the failing artifact."""


@dataclass
class StageResult:
    stage: str
    argv: List[str]
    rc: int                 # recorded, NEVER gating
    wall: float
    gates: List[str] = field(default_factory=list)   # human-readable gate lines
    log: str = ""


@dataclass
class PipelineConfig:
    """One reduction job on one search tree (phase0 layout).

    searchdir/            database/<datanames[pid]>/red_*
                          kinematics/<job>/kin_*
                          <job>/<workid>/{sch_*, multi_sch_intid, config/<k>/in_*}
    recmod_dir/           degrees.fflow, rec_<pid>_*, rrres, state
    """
    searchdir: str
    job: str
    recmod_dir: str
    datanames: List[str]                  # index = fflow prime id (0,1,2,...)
    nworks: int = 1
    nthreads: int = 6
    bin_dir: str = DEFAULT_BIN_DIR
    dumppoints_bin: Optional[str] = None  # default: bin_dir/dumppoints, else phase0 shim
    degrees_file: Optional[str] = None    # default: recmod_dir/degrees.fflow
    maxdeg: int = 1000
    rec_prefix: str = "rec"
    rrres_name: str = "rrres"
    stage_timeout: float = 3600.0
    dynamicrr_timeout: float = 600.0

    # -- derived paths ------------------------------------------------------
    @property
    def jobdir(self) -> str:
        return os.path.join(self.searchdir, self.job)

    def workdir(self, workid: int) -> str:
        return os.path.join(self.jobdir, str(workid))

    def dbdir(self, pid: int) -> str:
        return os.path.join(self.searchdir, "database", self.datanames[pid])

    @property
    def degrees_path(self) -> str:
        return self.degrees_file or os.path.join(self.recmod_dir, "degrees.fflow")

    def system_path(self, pid: int) -> str:
        return os.path.join(self.jobdir, f"system_p{pid}.txt")

    def evallist_path(self, pid: int) -> str:
        return os.path.join(self.jobdir, f"evallist_p{pid}.txt")

    def points_path(self, pid: int) -> str:
        return os.path.join(self.jobdir, f"points_p{pid}.fflow")

    def eval_path(self, pid: int) -> str:
        return os.path.join(self.jobdir, "eval", f"eval_p{pid}_1.fflow")

    def binpath(self, name: str) -> str:
        if name == "dumppoints":
            if self.dumppoints_bin:
                return self.dumppoints_bin
            cand = os.path.join(self.bin_dir, "dumppoints")
            return cand if os.path.exists(cand) else PHASE0_DUMPPOINTS
        return os.path.join(self.bin_dir, name)

    def dims(self):
        """(nparsin, nparsout) from the degrees file (authoritative for the
        reconstruction stages)."""
        deg = DegreesFile.read(self.degrees_path)
        return deg.nparsin, deg.nparsout

    @property
    def ptsord(self) -> str:
        nparsin, _ = self.dims()
        if nparsin > 9:
            raise BladeGateError(
                "ptsord digit encoding breaks for >=10 parameters (iofflow.cpp:106)")
        return "".join(str(i) for i in range(1, nparsin + 1))


def _run(argv: Sequence[str], cwd: Optional[str], timeout: float,
         stage: str) -> StageResult:
    t0 = time.time()
    try:
        cp = subprocess.run(list(argv), cwd=cwd, timeout=timeout,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        rc, log = cp.returncode, cp.stdout.decode(errors="replace")
    except subprocess.TimeoutExpired as e:
        raise BladeGateError(f"{stage}: timeout after {timeout}s: {' '.join(argv)}") from e
    except OSError as e:
        raise BladeGateError(f"{stage}: could not exec {argv[0]}: {e}") from e
    return StageResult(stage, list(argv), rc, time.time() - t0, log=log)


# --------------------------------------------------------------------------
# input generation helpers
# --------------------------------------------------------------------------

def write_system_and_evallists(cfg: PipelineConfig) -> None:
    """Generate system_p<pid>.txt + evallist_p<pid>.txt from the tree itself
    (scheme + red_common), with ABSOLUTE paths so copies of a tree work.
    Single-block scope: targets = the G1 ids of every work, in work order."""
    targets: List[int] = []
    for w in range(cfg.nworks):
        sch = Scheme.read(cfg.workdir(w))
        targets += sch.intid[:sch.g1]
    for pid in range(len(cfg.datanames)):
        common = read_ints(os.path.join(cfg.dbdir(pid), "red_common"))
        prime, npara, nint, nmaster = common[0], common[1], common[2], common[3]
        if prime != big_uint_primes()[pid]:
            raise BladeGateError(
                f"{cfg.dbdir(pid)}/red_common: prime {prime} != BIG_UINT_PRIMES[{pid}] "
                f"= {big_uint_primes()[pid]} (dynamicrr indexes primes by id)")
        sysf = SystemFile(prime, npara, nint, nmaster, targets,
                          [os.path.abspath(cfg.jobdir)], cfg.datanames[pid])
        sysf.write(cfg.system_path(pid))
        EvalList([os.path.abspath(cfg.eval_path(pid))]).write(cfg.evallist_path(pid))


# --------------------------------------------------------------------------
# stage: redg1
# --------------------------------------------------------------------------

def run_redg1(cfg: PipelineConfig) -> StageResult:
    """redg1 <searchdir> <dataname0> <job> 0 <nworks> <nthreads> (redg1.c:16).
    Gates (rc always 0, searchalg.c:529-618): per work, flag content == 1;
    every processed config has consistent tmp_*; sum(tmp_nsol) == g1
    (ssolve asserts this as 'wrong relations', block.c:174)."""
    res = _run([cfg.binpath("redg1"), cfg.searchdir, cfg.datanames[0], cfg.job,
                "0", str(cfg.nworks), str(cfg.nthreads)],
               cwd=cfg.searchdir, timeout=cfg.stage_timeout, stage="redg1")
    for w in range(cfg.nworks):
        wd = cfg.workdir(w)
        flagp = os.path.join(wd, "flag")
        if not os.path.exists(flagp):
            raise BladeGateError(f"redg1: {flagp} missing")
        flag = open(flagp).read().strip()
        if flag != "1":
            raise BladeGateError(
                f"redg1: {flagp} == {flag!r} (g1 not reduced; grow the database "
                f"or add a higher-degree config). log tail: {res.log[-400:]}")
        sch = Scheme.read(wd)
        tot = 0
        k = 0
        while os.path.exists(os.path.join(wd, "config", str(k), "tmp_nvar")):
            tmp = TmpTemplate.read(os.path.join(wd, "config", str(k)))  # dims gate inside
            tot += tmp.nsol
            k += 1
        if k == 0:
            raise BladeGateError(f"redg1: {wd}: flag=1 but no config has tmp_nvar")
        if tot != sch.g1:
            raise BladeGateError(
                f"redg1: {wd}: sum(tmp_nsol)={tot} != g1={sch.g1} "
                f"('wrong relations' would kill ssolve)")
        res.gates.append(f"work {w}: flag=1, {k} configs, sum(tmp_nsol)={tot}==g1")
    return res


# --------------------------------------------------------------------------
# stage: fitrel
# --------------------------------------------------------------------------

def run_fitrel(cfg: PipelineConfig, pid: int) -> StageResult:
    """fitrel <searchdir> <dataname_pid> <job> 0 <nworks> <nthreads>.
    Gates: for every config k with tmp_nvar: fit table <work>/fit/<db>/<k>
    exists with word count == nrows*tmp_nvar and nrows >= tmp_nsol (ssolve
    reads the FIRST tmp_nsol rows, block.c:122-124); tmp_nsol==0 => empty
    file REQUIRED to exist.  MEASURED FOOTGUN: flag can be 1 over an empty
    table with tmp_nsol>0 -- the dimension gate is the real gate; flag==0 is
    still trusted as a failure signal (self_check fired, searchalg.c:667)."""
    db = cfg.datanames[pid]
    res = _run([cfg.binpath("fitrel"), cfg.searchdir, db, cfg.job,
                "0", str(cfg.nworks), str(cfg.nthreads)],
               cwd=cfg.searchdir, timeout=cfg.stage_timeout, stage=f"fitrel[{db}]")
    for w in range(cfg.nworks):
        wd = cfg.workdir(w)
        fitdir = os.path.join(wd, "fit", db)
        flagp = os.path.join(fitdir, "flag")
        flag = open(flagp).read().strip() if os.path.exists(flagp) else "<missing>"
        if flag == "0":
            raise BladeGateError(f"fitrel: {flagp} == 0 (self_check failed on {db})")
        k = 0
        while os.path.exists(os.path.join(wd, "config", str(k), "tmp_nvar")):
            tmp = TmpTemplate.read(os.path.join(wd, "config", str(k)))
            tab_path = os.path.join(fitdir, str(k))
            if not os.path.exists(tab_path):
                raise BladeGateError(f"fitrel: fit table {tab_path} missing")
            tab = FitTable.read(tab_path, tmp.nvar)
            if tab.nrow < tmp.nsol:
                raise BladeGateError(
                    f"fitrel: {tab_path}: {tab.nrow} rows < tmp_nsol={tmp.nsol} "
                    f"(x tmp_nvar={tmp.nvar}) -- EMPTY/short fit nullspace; "
                    f"flag was {flag!r} (flags lie, dimensions do not)")
            res.gates.append(f"work {w} cfg {k}: fit {tab.nrow}x{tmp.nvar} "
                             f">= {tmp.nsol}x{tmp.nvar}")
            k += 1
    return res


# --------------------------------------------------------------------------
# stage: dumppoints (shim)
# --------------------------------------------------------------------------

def run_dumppoints(cfg: PipelineConfig, pid: int) -> StageResult:
    """dumppoints <nparsin> <nparsout> <pid> <maxdeg> <degrees> <points>.
    Sample points MUST come from fflow's own code path (structured
    Newton/Thiele grids; opt 7-tuple identical to recmod's).  Gates: header
    n == nparsin + flags_size, nsamples >= 1, every row prime ==
    BIG_UINT_PRIMES[pid]."""
    nparsin, nparsout = cfg.dims()
    pts_path = cfg.points_path(pid)
    res = _run([cfg.binpath("dumppoints"), str(nparsin), str(nparsout), str(pid),
                str(cfg.maxdeg), cfg.degrees_path, pts_path],
               cwd=cfg.jobdir, timeout=cfg.stage_timeout, stage=f"dumppoints[p{pid}]")
    if not os.path.exists(pts_path):
        raise BladeGateError(f"dumppoints: {pts_path} missing")
    pts = PointsFile.read(pts_path)   # word-count gate inside
    fs = flags_size(nparsout)
    if pts.n != nparsin + fs:
        raise BladeGateError(f"dumppoints: {pts_path}: header n={pts.n} != "
                             f"nparsin+flags_size={nparsin}+{fs}")
    if len(pts.rows) < 1:
        raise BladeGateError(f"dumppoints: {pts_path}: zero samples")
    p = big_uint_primes()[pid]
    for i in range(len(pts.rows)):
        if pts.rows[i][nparsin] != p:
            raise BladeGateError(f"dumppoints: {pts_path}: row {i} prime "
                                 f"{pts.rows[i][nparsin]} != BIG_UINT_PRIMES[{pid}]={p}")
    res.gates.append(f"{len(pts.rows)} samples, n={pts.n}, prime ok")
    return res


# --------------------------------------------------------------------------
# stage: ssolve
# --------------------------------------------------------------------------

def run_ssolve(cfg: PipelineConfig, pid: int) -> StageResult:
    """ssolve <nparsin> <pid> 0 <nsamples> <nthreads> <system> <ptsord>
    <points> <eval> (ssolve.cpp:37), CWD = job dir (average_time lands there).
    Gates (rc 0 even on undersized output -- mutation variant E): eval file
    size == 8*(4 + nrows*(nparsin+1+flags_size+nparsout)) with nparsout from
    the degrees file; header fields consistent."""
    nparsin, nparsout = cfg.dims()
    pts = PointsFile.read(cfg.points_path(pid))
    nsamples = len(pts.rows)
    eval_path = cfg.eval_path(pid)
    os.makedirs(os.path.dirname(eval_path), exist_ok=True)
    res = _run([cfg.binpath("ssolve"), str(nparsin), str(pid), "0", str(nsamples),
                str(cfg.nthreads), cfg.system_path(pid), cfg.ptsord,
                cfg.points_path(pid), eval_path],
               cwd=cfg.jobdir, timeout=cfg.stage_timeout, stage=f"ssolve[p{pid}]")
    if not os.path.exists(eval_path):
        raise BladeGateError(f"ssolve: {eval_path} missing. log tail: {res.log[-400:]}")
    want = EvalFile.expected_bytes(nparsin, nparsout, nsamples)
    got = os.path.getsize(eval_path)
    if got != want:
        raise BladeGateError(
            f"ssolve: {eval_path}: {got} bytes != expected {want} "
            f"(nparsin={nparsin}, nparsout={nparsout}, nrows={nsamples}) -- "
            f"learned nparsout differs from degrees file; rc was {res.rc} (rc lies)")
    ev = EvalFile.read(eval_path)
    if ev.nparsin != nparsin or ev.nout_plus_flags != nparsout + flags_size(nparsout) \
            or len(ev.rows) != nsamples:
        raise BladeGateError(f"ssolve: {eval_path}: header inconsistent "
                             f"({ev.nparsin},{ev.nout_plus_flags},{len(ev.rows)})")
    res.gates.append(f"eval {got}B == expected, {nsamples} rows x "
                     f"{ev.row_width} words")
    return res


# --------------------------------------------------------------------------
# stage: recmod  (generic runner reusable by ratrec)
# --------------------------------------------------------------------------

def run_recmod_generic(recmod_bin: str, nparsin: int, nparsout: int, pid: int,
                       nthreads: int, degrees_path: str, evallist_path: str,
                       workdir: str, prefix: str, maxdeg: int,
                       timeout: float = 3600.0) -> StageResult:
    """recmod <nparsin> <nparsout> <pid> <nthreads> <degrees> <evallist>
    <prefix>_<pid> <maxdeg>, CWD = workdir (recmod.cc:12).
    Gates (rc 0 even on 'failed'): rec_<pid>_{coeff,mono,part} exist;
    part count == 2*nparsout; len(coeff) == sum(part); len(mono words) ==
    nparsin*sum(part); for pid>0, _mono and _part byte-identical to pid 0
    (fflow normalization makes the pattern prime-independent)."""
    res = _run([recmod_bin, str(nparsin), str(nparsout), str(pid), str(nthreads),
                degrees_path, evallist_path, f"{prefix}_{pid}", str(maxdeg)],
               cwd=workdir, timeout=timeout, stage=f"recmod[p{pid}]")
    paths = {s: os.path.join(workdir, f"{prefix}_{pid}_{s}")
             for s in ("coeff", "mono", "part")}
    for s, p in paths.items():
        if not os.path.exists(p):
            raise BladeGateError(f"recmod: {p} missing (recmod prints 'failed' "
                                 f"with rc 0). log tail: {res.log[-400:]}")
    part = RecPart.read(paths["part"])
    if len(part.pairs) != nparsout:
        raise BladeGateError(f"recmod: {paths['part']}: {len(part.pairs)} functions "
                             f"!= nparsout={nparsout}")
    coeff = RecCoeff.read(paths["coeff"])
    if len(coeff.tokens) != part.total_terms:
        raise BladeGateError(f"recmod: {paths['coeff']}: {len(coeff.tokens)} coeffs "
                             f"!= sum(part)={part.total_terms}")
    mono = RecMono.read(paths["mono"], nparsin)
    if len(mono.exps) != part.total_terms:
        raise BladeGateError(f"recmod: {paths['mono']}: {len(mono.exps)} monomials "
                             f"!= sum(part)={part.total_terms}")
    if pid > 0:
        for s in ("mono", "part"):
            p0 = os.path.join(workdir, f"{prefix}_0_{s}")
            if open(paths[s], "rb").read() != open(p0, "rb").read():
                raise BladeGateError(
                    f"recmod: {paths[s]} differs from {p0} -- monomial pattern "
                    f"not prime-independent (corrupt evaluations?)")
    res.gates.append(f"{len(part.pairs)} funcs, {part.total_terms} terms, "
                     f"mono/part cross-prime ok" if pid > 0 else
                     f"{len(part.pairs)} funcs, {part.total_terms} terms")
    return res


def run_recmod(cfg: PipelineConfig, pid: int) -> StageResult:
    nparsin, nparsout = cfg.dims()
    return run_recmod_generic(cfg.binpath("recmod"), nparsin, nparsout, pid,
                              cfg.nthreads, cfg.degrees_path,
                              cfg.evallist_path(pid), cfg.recmod_dir,
                              cfg.rec_prefix, cfg.maxdeg, cfg.stage_timeout)


# --------------------------------------------------------------------------
# stage: dynamicrr  (generic runner reusable by ratrec)
# --------------------------------------------------------------------------

def run_dynamicrr_generic(dynamicrr_bin: str, ncoeffs: int, nthreads: int,
                          workdir: str, prefix: str, out_name: str,
                          nprimes_available: int,
                          timeout: float = 600.0) -> StageResult:
    """dynamicrr 0 <ncoeffs> <nthreads> <prefix> <out> as a state-file
    coprocess in workdir (dynamicrr.cc:17).  Driver: whenever state==0,
    confirm <prefix>_<ip>_coeff exists for the next prime ip and set state=1;
    if dynamicrr asks for a prime beyond what exists, that IS the cross-prime
    check failing -> gate.  Gates: state file == 2, out file parses as
    exactly ncoeffs rationals.  Needs >= 2 primes (prime 0 fit + fresh-prime
    verify)."""
    if nprimes_available < 2:
        raise BladeGateError("dynamicrr: needs >= 2 primes (fit + verify)")
    argv = [dynamicrr_bin, "0", str(ncoeffs), str(nthreads), prefix, out_name]
    t0 = time.time()
    state_path = os.path.join(workdir, "state")
    if os.path.exists(state_path):
        os.unlink(state_path)   # stale state from a previous run would desync
    proc = subprocess.Popen(argv, cwd=workdir, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT)
    ip = 0
    log = b""
    try:
        while True:
            if proc.poll() is not None:
                break
            if time.time() - t0 > timeout:
                proc.kill()
                raise BladeGateError(f"dynamicrr: timeout after {timeout}s "
                                     f"(state driver stuck at prime {ip})")
            try:
                state = open(state_path).read().strip()
            except FileNotFoundError:
                state = None
            if state == "0":
                if ip >= nprimes_available:
                    proc.kill()
                    raise BladeGateError(
                        f"dynamicrr: requested prime {ip} but only "
                        f"{nprimes_available} primes staged -- cross-prime "
                        f"rational-reconstruction check FAILED (corrupt or "
                        f"insufficient evaluations)")
                coeff = os.path.join(workdir, f"{prefix}_{ip}_coeff")
                if not os.path.exists(coeff):
                    proc.kill()
                    raise BladeGateError(f"dynamicrr: {coeff} missing at handshake")
                with open(state_path, "w") as f:
                    f.write("1")
                ip += 1
            time.sleep(0.05)
        log, _ = proc.communicate(timeout=10)
    finally:
        if proc.poll() is None:
            proc.kill()
    res = StageResult("dynamicrr", argv, proc.returncode, time.time() - t0,
                      log=log.decode(errors="replace"))
    state = open(state_path).read().strip() if os.path.exists(state_path) else "<missing>"
    if state != "2":
        raise BladeGateError(f"dynamicrr: {state_path} == {state!r} != 2 "
                             f"(no certified rational reconstruction). "
                             f"log tail: {res.log[-400:]}")
    out_path = os.path.join(workdir, out_name)
    if not os.path.exists(out_path):
        raise BladeGateError(f"dynamicrr: {out_path} missing despite state==2")
    rr = RRRes.read(out_path)
    if len(rr.values) != ncoeffs:
        raise BladeGateError(f"dynamicrr: {out_path}: {len(rr.values)} rationals "
                             f"!= ncoeffs={ncoeffs}")
    res.gates.append(f"state==2, {ncoeffs} rationals, {ip} primes consumed")
    return res


def run_dynamicrr(cfg: PipelineConfig) -> StageResult:
    coeff0 = RecCoeff.read(os.path.join(cfg.recmod_dir, f"{cfg.rec_prefix}_0_coeff"))
    return run_dynamicrr_generic(cfg.binpath("dynamicrr"), len(coeff0.tokens),
                                 cfg.nthreads, cfg.recmod_dir, cfg.rec_prefix,
                                 cfg.rrres_name, len(cfg.datanames),
                                 cfg.dynamicrr_timeout)


# --------------------------------------------------------------------------
# run_all
# --------------------------------------------------------------------------

def run_all(cfg: PipelineConfig, regenerate_system: bool = True
            ) -> Dict[str, StageResult]:
    """Full invocation order (FORMATS.md 'Invocation order'): redg1 on prime 0,
    then per prime: fitrel, dumppoints, ssolve, recmod; then dynamicrr.
    Inputs (database/kinematics/scheme/configs/degrees) must already exist.
    Returns {stage_name: StageResult}; raises BladeGateError on any gate."""
    results: Dict[str, StageResult] = {}
    if regenerate_system:
        write_system_and_evallists(cfg)
    results["redg1"] = run_redg1(cfg)
    for pid in range(len(cfg.datanames)):
        results[f"fitrel_p{pid}"] = run_fitrel(cfg, pid)
        results[f"dumppoints_p{pid}"] = run_dumppoints(cfg, pid)
        results[f"ssolve_p{pid}"] = run_ssolve(cfg, pid)
        results[f"recmod_p{pid}"] = run_recmod(cfg, pid)
    results["dynamicrr"] = run_dynamicrr(cfg)
    return results
