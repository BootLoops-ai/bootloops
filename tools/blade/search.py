"""blade.search -- SEARCH-LOOP driver: Wolfram-free port of
BLSearch/SearchRelation.wl (BLSearchRelationsBasic escalation, BLRedG1,
BLFitRelations, fitnumber bisection, SearchRelation.wl:139-264) orchestrating
scheme.py (block scheme) + ansatz.py (polynomial ansatz) + pipeline.py
(gated C stages) into `blade_search(...)` -> exact-rational BT relations.

Upstream semantics ported (line cites):
  * round loop (SearchRelation.wl:181-197): redg1; if jobFinishedQ stop; else
    classify each work:
      - TERMINATED (all in_ configs have tmp_, flag==0, wl:67-73): the ansatz
        is exhausted -> raise the composite degree by one step (wl:188 -- the
        +1 is on the LAST group digit; for single-group families the port
        exposes `level_step` because the banked recipes escalate the composite
        level directly);
      - MISSING DATABASE (fewer tmp_ than in_, flag==0, wl:78-84: redg1's
        search_relation self-check failed mid-config) -> double the database
        (wl:190,195: BLDatabase[dataname, BLDataNumber[dataname]]).
    Round-limit Break (default 10, wl:153) -> $Failed (here: SearchFailed).
  * ansatzSize default = 7 on the last group (wl:140); dataSize default =
    10/30/100/300/1000 by #vars, rounded up to the thread multiple (wl:139).
  * fitnumber bisection (wl:219-246): shrink red_common's nps by bisection
    until the relative step < FitDataPercent (1/5) or the thread-block
    quotient stops moving; restore searchnumber at the end (wl:252).
  * BLGenerateAnsatzForJob breakpoint (PolynomialAnsatz.wl:195): a work whose
    flag is already 1 gets no new configs.

Documented deviations (each disclosed at its definition):
  * per-work escalation: works are driven independently (redg1 accepts a
    workid range); upstream escalates the whole job with one counter.  For
    multi-block schemes this lets cheap blocks close while an expensive block
    is still escalating, and lets a wall STOP one block only.
  * `datasize_policy="structural"`: instead of reactive doubling only, the
    driver may pre-grow the database to a structural lower bound computed
    from the config support (sum_m min(nps-1, S_m) >= nvar - g1 is NECESSARY
    for the constraint matrix to have full column rank outside the nullspace;
    S_m = #ansatz vars whose integral projects on master m).  Reactive
    doubling stays as the fallback.  The faithful default remains "reactive".
  * databases are GROWN by an injected `grower` callback (our Kira-expression
    GF(p) evaluator) -- upstream calls Wolfram BLDatabase.  No Wolfram here.

Gate discipline: artifact gates only (file content / dimensions / byte
compares); rc recorded, never gating; flag==0 trusted only as failure signal.
Every redg1/fitrel subprocess is launched with a hard timeout and killed by
its exact PID on expiry (resource rule), and the event is reported in the
round ledger, never swallowed.

Rational export: fit tables at >=3 primes -> row-wise zero-pattern identity
gate -> CRT over a GROWING prime prefix + Wang from the joint modulus ->
MANDATORY verification of EVERY entry at >=1 held-out prime never in the CRT
set (upstream growth loop: BLSolve/ReconstructRelation.wl:61-92, see
reconstruct_fit_row) -> relations as exact Fractions keyed by (global
integral id, monomial exponent vector).  When all available primes are
exhausted the export raises the TYPED ReconstructionOverflow (fail-closed,
nothing emitted); export_relations can grow the prime set via a prime_adder
callback (add_fit_prime: new BIG_UINT prime database + fitrel).  The heavy
fflow machinery (ratrec.py) is not needed because BT coefficients are
rational NUMBERS, not functions.

Reduced-analytic parameter search ("nana" mode, bottom of this module):
upstream BLAdaptiveSearch searches the BT system analytic in only
nana < n_params parameters, with the remaining parameters numerically pinned
per instance (AdaptiveSearch.wl:39-63; Database.wl:114-141;
SolveSemiBL.wl:296-304,404-410).  See the section banner for the full
semantics, the FAIL-CLOSED consumer contract (NanaContractError,
nana_pins.json, doc["reduced_analytic"]) and the per-pin instance logic
(plan_pin_instances / fit_at_pin).
"""
from __future__ import annotations

import json
import math
import os
import re
import shutil
import subprocess
import time
from dataclasses import dataclass, field, replace as _dc_replace
from fractions import Fraction
from typing import Callable, Dict, List, Optional, Sequence, Tuple

from . import formats
from .formats import (BladeFormatError, Database, FitTable, Kinematics,
                      Scheme, TmpTemplate, read_ints, write_fab_rows)
from .ansatz import FamilyAnsatz
from .pipeline import BladeGateError, StageResult, DEFAULT_BIN_DIR


# --------------------------------------------------------------------------
# red_common data-number handles (BLSetDataNumber/BLDataNumber, wl:48-56)
# --------------------------------------------------------------------------

def data_number(dbdir: str) -> int:
    return read_ints(os.path.join(dbdir, "red_common"))[-1]


def set_data_number(dbdir: str, n: int) -> None:
    """Patch ONLY the last red_common field (wl:48-53).  The C readers load
    exactly nps rows from red_ps/red_data, so extra trailing rows are inert --
    this is how the bisection shrinks the effective database."""
    common = read_ints(os.path.join(dbdir, "red_common"))
    have = _file_rows(dbdir)
    if n > have:
        raise BladeGateError(f"set_data_number({n}) > rows on disk ({have}) "
                             f"in {dbdir} -- grow the database first")
    common[-1] = n
    write_fab_rows(os.path.join(dbdir, "red_common"), [common])


def _file_rows(dbdir: str) -> int:
    """Actual point rows present in red_ps (the physical pool)."""
    common = read_ints(os.path.join(dbdir, "red_common"))
    npara = common[1]
    words = len(read_ints(os.path.join(dbdir, "red_ps")))
    if words % npara:
        raise BladeFormatError(f"{dbdir}/red_ps: {words} words % npara={npara}")
    return words // npara


# --------------------------------------------------------------------------
# job state (workTerminatedQ / workMissingDatabaseQ / jobFinishedQ, wl:67-95)
# --------------------------------------------------------------------------

def _count_seq(pattern: Callable[[int], str]) -> int:
    k = 0
    while os.path.exists(pattern(k)):
        k += 1
    return k


def n_in_configs(workdir: str) -> int:
    return _count_seq(lambda k: os.path.join(workdir, "config", str(k), "in_nvar"))


def n_tmp_configs(workdir: str) -> int:
    return _count_seq(lambda k: os.path.join(workdir, "config", str(k), "tmp_nvar"))


def work_flag(workdir: str) -> Optional[int]:
    p = os.path.join(workdir, "flag")
    if not os.path.exists(p):
        return None
    return read_ints(p)[0]


def work_state(workdir: str) -> str:
    """'finished' | 'terminated' | 'missing' | 'virgin' (wl:67-95)."""
    flag = work_flag(workdir)
    if flag == 1:
        return "finished"
    if flag is None:
        return "virgin"
    ins, tmps = n_in_configs(workdir), n_tmp_configs(workdir)
    return "terminated" if ins == tmps else "missing"


# --------------------------------------------------------------------------
# escalation config
# --------------------------------------------------------------------------

@dataclass
class SearchFamily:
    """Family recipe + tree layout for one search job."""
    searchdir: str                       # contains database/, kinematics/, <job>/
    job: str
    fam: FamilyAnsatz                    # the polynomial-ansatz recipe (nmax = ceiling)
    work_weights: List[List[int]]        # RAW integral dims per work (block order)
    datanames: List[str]                 # index == fflow prime id
    bin_dir: str = DEFAULT_BIN_DIR
    nthreads: int = 6
    # escalation knobs (upstream options, SearchRelation.wl:153)
    ansatz_size: int = 7                 # "AnsatzSize" (initial composite level)
    level_base: int = 0                  # first config's level
    level_step: int = 1                  # composite-level increment per round
    data_size: int = 30                  # "DataSize" (initial #points)
    break_rounds: int = 10               # "Break"
    fit_data_percent: Fraction = Fraction(1, 5)   # "FitDataPercent"
    datasize_policy: str = "reactive"    # "reactive" (faithful) | "structural"
    structural_margin: float = 1.10      # structural bound safety factor
    round_timeout: float = 1200.0        # 20 min hard cap per redg1 round
    mem_cap_bytes: int = 24 << 30        # refuse a round projected above this
    # STOP-line (blade_family turnkey): when set, REFUSE to launch
    # a round whose projected wall exceeds this line.  Projection is measured,
    # never guessed: last round's wall x growth, growth = ratio of the last
    # two measured walls (>= 1) once two rounds exist, else
    # stop_growth_default (conservative prior from the dbox L6-L20 ladder,
    # measured per-level ratios 1.7-8.6x).  The first round has no
    # measurement and is bounded by round_timeout (kill by exact PID).
    stop_line_seconds: Optional[float] = None
    stop_growth_default: float = 4.0
    grower: Optional[Callable[[int], None]] = None  # grow ALL prime pools to n pts
    log: Callable[[str], None] = print

    @property
    def jobdir(self) -> str:
        return os.path.join(self.searchdir, self.job)

    def workdir(self, w: int) -> str:
        return os.path.join(self.jobdir, str(w))

    def dbdir(self, pid: int) -> str:
        return os.path.join(self.searchdir, "database", self.datanames[pid])

    def levels(self, a: int) -> List[int]:
        """Config levels for escalation counter a (config index == position)."""
        return list(range(self.level_base, a + 1, self.level_step))


@dataclass
class RoundRecord:
    work: int
    rnd: int
    level: int
    nvar_top: int
    npoints: int
    proj_mat_bytes: int
    wall: float
    rc: int
    state_after: str
    note: str = ""


# --------------------------------------------------------------------------
# ansatz emission for one work at counter a (BLGenerateAnsatzForJob semantics)
# --------------------------------------------------------------------------

def write_kinematics(cfg: SearchFamily, a: int) -> None:
    """writeKinematics at the current ceiling (PolynomialAnsatz.wl:226-230).
    kin_prefix keeps row ids prefix-stable (ansatz.py O1).  NEVER shrinks:
    works are driven independently here (port deviation), and a lower-level
    work must not truncate rows an earlier work's tmp/fit already reference."""
    kndir = os.path.join(cfg.searchdir, "kinematics", cfg.job)
    kn = cfg.fam.kin_prefix(a) if a < cfg.fam.nmax else cfg.fam.kinematics()
    common = os.path.join(kndir, "kin_common")
    if os.path.exists(common) and read_ints(common)[0] >= kn.nmono:
        return
    kn.write(kndir)


def write_ansatz_for_work(cfg: SearchFamily, w: int, a: int) -> int:
    """Append missing configs of work w up to counter a; returns top in_nvar.
    Breakpoint: flag==1 -> nothing written (PolynomialAnsatz.wl:195)."""
    wd = cfg.workdir(w)
    levels = cfg.levels(a)
    if work_flag(wd) == 1:
        return _top_nvar(wd)
    have = n_in_configs(wd)
    cfgs = cfg.fam.block_configs(cfg.work_weights[w], levels)
    if have > len(cfgs):
        raise BladeGateError(f"{wd}: {have} configs on disk > {len(cfgs)} planned "
                             f"(level ledger out of sync)")
    for k in range(have, len(cfgs)):
        cfgs[k].write(os.path.join(wd, "config", str(k)))
    return cfgs[-1].nvar if cfgs else 0


def _top_nvar(workdir: str) -> int:
    n = n_in_configs(workdir)
    if n == 0:
        return 0
    return read_ints(os.path.join(workdir, "config", str(n - 1), "in_nvar"))[0]


# --------------------------------------------------------------------------
# structural sizing (support per master) -- for projections & "structural"
# datasize policy
# --------------------------------------------------------------------------

def master_supports(db: Database, sch: Scheme, workdir: str) -> List[int]:
    """S_m = #ansatz variables (over ALL in_ configs' union, i.e. the top
    cumulative config) whose integral projects on master m (red_posi support,
    generate_constraint searchalg.c:68-86)."""
    k = n_in_configs(workdir)
    if k == 0:
        return [0] * db.nmaster
    top = formats.AnsatzConfig.read(os.path.join(workdir, "config", str(k - 1)))
    per_local: Dict[int, int] = {}
    for lid, _ in top.var:
        per_local[lid] = per_local.get(lid, 0) + 1
    # posi(m, localid) != -1  <=>  global id sch.intid[localid] in master m's group
    supports = []
    for m in range(db.nmaster):
        grp = set(db.posi[db.part[m]:db.part[m + 1]])
        supports.append(sum(cnt for lid, cnt in per_local.items()
                            if sch.intid[lid] in grp))
    return supports


def structural_min_points(db: Database, sch: Scheme, workdir: str,
                          margin: float) -> int:
    """Smallest N with sum_m min(N-1, S_m) >= nvar - g1 (necessary for the
    relation search to pin all non-null directions), x margin."""
    nvar = _top_nvar(workdir)
    need = max(0, nvar - sch.g1)
    S = sorted(master_supports(db, sch, workdir), reverse=True)
    lo, hi = 2, max(2, need + 2)
    while lo < hi:
        mid = (lo + hi) // 2
        if sum(min(mid - 1, s) for s in S) >= need:
            hi = mid
        else:
            lo = mid + 1
    return int(math.ceil(lo * margin))


def projected_mat_bytes(db_npoints: int, supports: Sequence[int]) -> int:
    """Worst-case redg1 solve_master matrix: (nps-1) x max_m S_m llu words
    (searchalg.c:151 table_init_llu(mat, db->nps-1, size))."""
    return 8 * max(1, db_npoints - 1) * (max(supports) if supports else 1)


# --------------------------------------------------------------------------
# subprocess with exact-PID kill
# --------------------------------------------------------------------------

def _run_binary(argv: Sequence[str], cwd: str, timeout: float, stage: str,
                log: Callable[[str], None]) -> StageResult:
    t0 = time.time()
    proc = subprocess.Popen(list(argv), cwd=cwd, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT)
    try:
        out, _ = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        pid = proc.pid
        proc.kill()
        out, _ = proc.communicate()
        log(f"[{stage}] TIMEOUT after {timeout:.0f}s -- killed PID {pid} exactly")
        raise BladeGateError(f"{stage}: timeout after {timeout:.0f}s "
                             f"(killed PID {pid}); log tail: "
                             f"{out.decode(errors='replace')[-300:]}")
    return StageResult(stage, list(argv), proc.returncode, time.time() - t0,
                       log=out.decode(errors="replace"))


def run_redg1_works(cfg: SearchFamily, w_ini: int, w_fin: int,
                    timeout: Optional[float] = None) -> StageResult:
    """redg1 <searchdir> <dataname0> <job> <w_ini> <w_fin> <nthreads>
    (redg1.c:17).  NO flag gate here -- the escalation loop classifies."""
    return _run_binary([os.path.join(cfg.bin_dir, "redg1"), cfg.searchdir,
                        cfg.datanames[0], cfg.job, str(w_ini), str(w_fin),
                        str(cfg.nthreads)],
                       cwd=cfg.searchdir,
                       timeout=timeout or cfg.round_timeout,
                       stage=f"redg1[w{w_ini}:{w_fin}]", log=cfg.log)


def run_fitrel_works(cfg: SearchFamily, dataname: str, w_ini: int, w_fin: int,
                     timeout: Optional[float] = None) -> StageResult:
    res = _run_binary([os.path.join(cfg.bin_dir, "fitrel"), cfg.searchdir,
                       dataname, cfg.job, str(w_ini), str(w_fin),
                       str(cfg.nthreads)],
                      cwd=cfg.searchdir,
                      timeout=timeout or cfg.round_timeout,
                      stage=f"fitrel[{dataname},w{w_ini}:{w_fin}]", log=cfg.log)
    for w in range(w_ini, w_fin):
        gate_fit_work(cfg, dataname, w, res)
    return res


def gate_fit_work(cfg: SearchFamily, dataname: str, w: int,
                  res: Optional[StageResult] = None) -> bool:
    """Artifact gate for one work's fit dir (dimensions, never flags) --
    pipeline.run_fitrel semantics, parametrized.  Returns True/raises."""
    wd = cfg.workdir(w)
    fitdir = os.path.join(wd, "fit", dataname)
    flagp = os.path.join(fitdir, "flag")
    flag = read_ints(flagp)[0] if os.path.exists(flagp) else None
    if flag == 0:
        raise BladeGateError(f"fitrel: {flagp} == 0 (self_check failed)")
    k = 0
    while os.path.exists(os.path.join(wd, "config", str(k), "tmp_nvar")):
        tmp = TmpTemplate.read(os.path.join(wd, "config", str(k)))
        tab_path = os.path.join(fitdir, str(k))
        if not os.path.exists(tab_path):
            raise BladeGateError(f"fitrel: fit table {tab_path} missing")
        tab = FitTable.read(tab_path, tmp.nvar)
        if tab.nrow < tmp.nsol:
            raise BladeGateError(f"fitrel: {tab_path}: {tab.nrow} rows < "
                                 f"tmp_nsol={tmp.nsol} (flags lie, dimensions "
                                 f"do not; flag was {flag!r})")
        if res is not None:
            res.gates.append(f"w{w} cfg{k}: fit {tab.nrow}x{tmp.nvar} >= "
                             f"{tmp.nsol}x{tmp.nvar}")
        k += 1
    return True


# --------------------------------------------------------------------------
# the escalation loop (BLSearchRelationsBasic core, wl:154-202), per work
# --------------------------------------------------------------------------

@dataclass
class WorkOutcome:
    work: int
    closed: bool
    rounds: List[RoundRecord] = field(default_factory=list)
    stop_reason: str = ""


def escalate_work(cfg: SearchFamily, w: int) -> WorkOutcome:
    """Drive ONE work to flag==1 or a measured STOP.  Round semantics ==
    upstream Which[] (wl:186-196); Break counts redg1 attempts."""
    out = WorkOutcome(work=w, closed=False)
    if work_state(cfg.workdir(w)) == "finished":
        out.closed = True
        out.stop_reason = "already finished (flag==1 on disk)"
        return out
    sch = Scheme.read(cfg.workdir(w))
    a = cfg.ansatz_size
    have = n_in_configs(cfg.workdir(w))
    if have:                      # resume: counter follows configs on disk
        a = max(a, cfg.level_base + (have - 1) * cfg.level_step)
    write_kinematics(cfg, a)
    nvar = write_ansatz_for_work(cfg, w, a)
    pool0 = _file_rows(cfg.dbdir(0))
    want = min(cfg.data_size, pool0) if cfg.grower is None else cfg.data_size

    for rnd in range(cfg.break_rounds):
        # -- data sizing for this round
        if cfg.datasize_policy == "structural":
            dbh = _db_header(cfg.dbdir(0))
            smin = structural_min_points(dbh, sch, cfg.workdir(w),
                                         cfg.structural_margin)
            want = max(want, smin)
        want = _ensure_points(cfg, want)
        set_data_number(cfg.dbdir(0), want)

        # -- memory projection gate (never launch an unmeasured 30 GB round)
        dbh = _db_header(cfg.dbdir(0))
        supports = master_supports(dbh, sch, cfg.workdir(w))
        proj = projected_mat_bytes(want, supports)
        if proj > cfg.mem_cap_bytes:
            out.stop_reason = (f"memory projection {proj/2**30:.1f} GiB > cap "
                               f"{cfg.mem_cap_bytes/2**30:.0f} GiB at level "
                               f"{cfg.levels(a)[-1]} (nvar={nvar}, N={want}) "
                               f"-- STOP, re-plan trigger")
            cfg.log(f"[work {w}] {out.stop_reason}")
            return out

        # -- STOP-line: projected-next-round refusal (measured, never guessed)
        if cfg.stop_line_seconds is not None and out.rounds:
            walls = [r.wall for r in out.rounds]
            if len(walls) >= 2 and walls[-2] > 0:
                growth = max(1.0, walls[-1] / walls[-2])
                basis = f"measured ratio of last two rounds {growth:.2f}"
            else:
                growth = cfg.stop_growth_default
                basis = f"prior growth factor {growth:.1f} (one round measured)"
            projected = walls[-1] * growth
            if projected > cfg.stop_line_seconds:
                out.stop_reason = (
                    f"STOP line: projected next round {projected:.0f}s "
                    f"(last wall {walls[-1]:.1f}s x {basis}) > "
                    f"max-round-seconds {cfg.stop_line_seconds:.0f} -- "
                    f"refusing to launch (level {cfg.levels(a)[-1]}, "
                    f"nvar={nvar}, N={want})")
                cfg.log(f"[work {w}] {out.stop_reason}")
                return out

        try:
            res = run_redg1_works(cfg, w, w + 1)
        except BladeGateError as e:
            out.stop_reason = f"round {rnd}: {e}"
            return out
        st = work_state(cfg.workdir(w))
        out.rounds.append(RoundRecord(w, rnd, cfg.levels(a)[-1], nvar, want,
                                      proj, res.wall, res.rc, st))
        cfg.log(f"[work {w}] round {rnd}: level={cfg.levels(a)[-1]} "
                f"nvar={nvar} N={want} wall={res.wall:.1f}s -> {st}")
        if st == "finished":
            out.closed = True
            return out

        # -- upstream Which[] escalation (wl:186-196); one work is either
        #    terminated or missing, never both.
        if st == "terminated":
            a += cfg.level_step
            if a > cfg.fam.nmax:
                out.stop_reason = (f"level {a} beyond recipe nmax={cfg.fam.nmax}"
                                   f" -- STOP (kin_table ceiling)")
                return out
            write_kinematics(cfg, a)
            nvar = write_ansatz_for_work(cfg, w, a)
        elif st == "missing":
            want = want * 2
        else:
            out.stop_reason = f"round {rnd}: unknown pattern state={st} (wl:196)"
            return out

    out.stop_reason = (f"round-limit Break={cfg.break_rounds} reached "
                       f"(wl:201) -- STOP, re-plan trigger")
    cfg.log(f"[work {w}] {out.stop_reason}")
    return out


def _db_header(dbdir: str) -> Database:
    """Database with header + part/posi only (no ps/data load -- red_data can
    be 100s of MB; sizing needs structure only)."""
    common = read_ints(os.path.join(dbdir, "red_common"))
    prime, npara, nint, nmaster, nentry, nps = common
    db = Database(prime, npara, nint, nmaster, nentry, nps)
    db.part = read_ints(os.path.join(dbdir, "red_part"))
    db.posi = read_ints(os.path.join(dbdir, "red_posi"))
    return db


def _ensure_points(cfg: SearchFamily, want: int) -> int:
    """Grow the physical pools to >= want points (all primes, via grower);
    without a grower, cap at the pool (upstream would call BLDatabase)."""
    pool = _file_rows(cfg.dbdir(0))
    if want <= pool:
        return want
    if cfg.grower is None:
        cfg.log(f"[data] want {want} > pool {pool} and no grower -- capped")
        return pool
    t0 = time.time()
    cfg.grower(want)
    got = _file_rows(cfg.dbdir(0))
    if got < want:
        raise BladeGateError(f"grower returned {got} < requested {want}")
    cfg.log(f"[data] pool {pool} -> {got} points in {time.time()-t0:.1f}s")
    return want


# --------------------------------------------------------------------------
# fitnumber bisection (wl:219-246)
# --------------------------------------------------------------------------

def job_fit_q(cfg: SearchFamily, dataname: str, works: Sequence[int]) -> bool:
    try:
        for w in works:
            gate_fit_work(cfg, dataname, w)
        return True
    except BladeGateError:
        return False


def fitnumber_bisection(cfg: SearchFamily, works: Sequence[int],
                        pid: int = 0) -> int:
    """Find the minimal #points for which fitrel still fits (wl:219-246).
    Restores the full data number AND a valid full-data fit dir at the end."""
    dataname = cfg.datanames[pid]
    dbdir = cfg.dbdir(pid)
    searchnumber = data_number(dbdir)
    # full fit (deletes any stale dir first, wl:208)
    for w in works:
        shutil.rmtree(os.path.join(cfg.workdir(w), "fit", dataname),
                      ignore_errors=True)
    for w in works:
        run_fitrel_works(cfg, dataname, w, w + 1)
    # backup (wl:216)
    for w in works:
        d = os.path.join(cfg.workdir(w), "fit", dataname)
        shutil.rmtree(d + "_bak", ignore_errors=True)
        shutil.copytree(d, d + "_bak")

    lo, hi = 0, searchnumber
    nth = cfg.nthreads
    while True:
        mid = (lo + hi + 1) // 2                      # Ceiling[Mean[range]]
        if not (Fraction(hi - mid, mid) > cfg.fit_data_percent
                and (hi - 1) // nth > (mid - 1) // nth):
            break
        set_data_number(dbdir, mid)
        cfg.log(f"[bisect {dataname}] trying N={mid} in [{lo},{hi}]")
        for w in works:
            shutil.rmtree(os.path.join(cfg.workdir(w), "fit", dataname),
                          ignore_errors=True)
        try:
            for w in works:
                run_fitrel_works(cfg, dataname, w, w + 1)
            ok = job_fit_q(cfg, dataname, works)
        except BladeGateError:
            ok = False
        if ok:
            hi = mid
        else:
            lo = mid
    fitnumber = hi
    # restore a valid fit dir + the search data number (wl:233-252)
    if not job_fit_q(cfg, dataname, works):
        for w in works:
            d = os.path.join(cfg.workdir(w), "fit", dataname)
            shutil.rmtree(d, ignore_errors=True)
            shutil.copytree(d + "_bak", d)
    set_data_number(dbdir, searchnumber)
    with open(os.path.join(cfg.jobdir, "fitnumber"), "w") as f:
        f.write(f"{fitnumber}\n")
    cfg.log(f"[bisect {dataname}] fitnumber={fitnumber} (search N={searchnumber})")
    return fitnumber


# --------------------------------------------------------------------------
# CRT + Wang rational reconstruction (rational export)
# --------------------------------------------------------------------------

def crt_pair(a0: int, m0: int, a1: int, m1: int) -> Tuple[int, int]:
    inv = pow(m0, -1, m1)
    x = (a0 + ((a1 - a0) * inv % m1) * m0) % (m0 * m1)
    return x, m0 * m1


def wang_ratrec(a: int, m: int) -> Optional[Fraction]:
    """Standard Wang half-Euclid: a mod m -> n/d with |n|,d <= sqrt(m/2)."""
    bound = math.isqrt(m // 2)
    r0, r1 = m, a % m
    s0, s1 = 0, 1
    while r1 > bound:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound or math.gcd(r1, abs(s1)) != 1:
        return None
    return Fraction(r1, s1) if s1 > 0 else Fraction(-r1, -s1)


class ReconstructionOverflow(BladeGateError):
    """Typed refusal: the joint modulus of every admissible CRT prime prefix
    is too small (or a held-out prime rejected every candidate) -- the caller
    must FIT AT AN ADDITIONAL PRIME and retry.  Nothing is ever emitted on
    this path (fail-closed)."""

    def __init__(self, msg: str, col: Optional[int] = None,
                 n_primes: int = 0, modulus_bits: int = 0):
        super().__init__(msg)
        self.col = col
        self.n_primes = n_primes
        self.modulus_bits = modulus_bits


def reconstruct_fit_row(vals_per_prime: Sequence[Sequence[int]],
                        primes: Sequence[int], min_crt: int = 2,
                        n_verify_min: int = 1) -> Tuple[List[Fraction], int]:
    """One relation row given at >= min_crt + n_verify_min primes.

    Upstream semantics (BLSolve/ReconstructRelation.wl, BLReconstructRelations):
    the CRT prime set GROWS one fitted prime per failed round
    (`AppendTo[primeset, pid]` inside `While[rest =!= {}]`, wl:77-91), the
    rationals come from FFRatRec over the JOINT modulus `Times@@primes`
    (wl:66-67 via chineseRemainderTheorm, wl:22-24), and acceptance is
    MANDATORY agreement with the benchmark fit at a prime held OUT of the CRT
    set (wl:61-62 benchmark/testprime = prime 0; wl:69
    `If[FFRatMod[fit, testprime] =!= benchmark, Return[0], ...]`).

    Port mapping (documented deviations; semantics preserved):
      * upstream holds out prime id 0 and grows the set 1,2,...; the port
        grows the PREFIX primes[:k] and holds out primes[k:] (always
        >= n_verify_min primes, so the newest prime is NEVER in the CRT set).
        This keeps the banked phase-3 db export byte-identical
        (k=2: CRT(p0,p1), verify p2).
      * upstream's first round uses a single prime; the port starts at
        min_crt=2 (stricter: any k=1-representable coefficient is also found
        identically at k=2).
      * verification is per-entry at EVERY held-out prime -- a superset of
        the upstream full-row compare at the single testprime.

    Returns (fractions, k) with k = #CRT primes actually used.  Raises the
    TYPED ReconstructionOverflow when no admissible k verifies -- the caller
    must fit at one more prime (export_relations grows via prime_adder);
    NOTHING is emitted on failure."""
    if len(primes) != len(vals_per_prime):
        raise BladeGateError("fit rows: #primes != #value rows")
    if len(primes) < min_crt + n_verify_min:
        raise BladeGateError(f"reconstruct_fit_row needs >= {min_crt} CRT + "
                             f"{n_verify_min} held-out primes; got {len(primes)}")
    ncol = len(vals_per_prime[0])
    pat0 = [v == 0 for v in vals_per_prime[0]]
    for vp in vals_per_prime[1:]:
        if len(vp) != ncol:
            raise BladeGateError("fit rows: column count differs across primes")
        if [v == 0 for v in vp] != pat0:
            raise BladeGateError("fit rows: zero pattern differs across primes "
                                 "(pivot structure unstable -- cannot CRT)")
    last_err = "no CRT attempt made"
    for k in range(min_crt, len(primes) - n_verify_min + 1):
        out: List[Fraction] = []
        ok = True
        for j in range(ncol):
            x, m = vals_per_prime[0][j] % primes[0], primes[0]
            for i in range(1, k):
                x, m = crt_pair(x, m, vals_per_prime[i][j], primes[i])
            fr = wang_ratrec(x, m)
            if fr is None:
                last_err = (f"col {j}: Wang reconstruction failed at k={k} "
                            f"(joint modulus {m.bit_length()} bits)")
                ok = False
                break
            for extra in range(k, len(primes)):
                p = primes[extra]
                if fr.denominator % p == 0:   # unlucky held-out prime: the
                    ok = False                # candidate's mod-p image is
                    last_err = (f"col {j}: candidate {fr} (k={k}) has "
                                f"denominator divisible by held-out prime "
                                f"p[{extra}]={p}")   # undefined -- grow k
                    break
                if fr.numerator % p != fr.denominator * vals_per_prime[extra][j] % p:
                    last_err = (f"col {j}: candidate {fr} (k={k}) FAILS "
                                f"verification at held-out prime p[{extra}]={p}")
                    ok = False
                    break
            if not ok:
                break
            out.append(fr)
        if ok:
            return out, k
    raise ReconstructionOverflow(
        f"rational reconstruction UNVERIFIED with all {len(primes)} available "
        f"primes (max CRT set {len(primes) - n_verify_min}, >= {n_verify_min} "
        f"held out): {last_err} -- fit at an additional prime and retry; "
        f"refusing to emit", n_primes=len(primes),
        modulus_bits=sum(p.bit_length() for p in primes[:len(primes) - n_verify_min]))


def _collect_relations(cfg: SearchFamily, works: Sequence[int],
                       labels: Sequence[str]
                       ) -> Tuple[list, List[int], int, int]:
    """One reconstruction pass over the CURRENT cfg.datanames.  Returns
    (relations, primes, kmax, height_bits); raises ReconstructionOverflow if
    any row cannot be verified with the available primes."""
    primes = [read_ints(os.path.join(cfg.dbdir(pid), "red_common"))[0]
              for pid in range(len(cfg.datanames))]
    kn = Kinematics.read(os.path.join(cfg.searchdir, "kinematics", cfg.job))
    params = list(cfg.fam.search_params)
    relations, kmax, height_bits = [], 2, 0
    for w in works:
        wd = cfg.workdir(w)
        sch = Scheme.read(wd)
        k = 0
        while os.path.exists(os.path.join(wd, "config", str(k), "tmp_nvar")):
            tmp = TmpTemplate.read(os.path.join(wd, "config", str(k)))
            tabs = [FitTable.read(os.path.join(wd, "fit", dn, str(k)), tmp.nvar)
                    for dn in cfg.datanames]
            for tab, dn in zip(tabs, cfg.datanames):
                if tab.nrow < tmp.nsol:
                    raise BladeGateError(f"{wd}/fit/{dn}/{k}: {tab.nrow} rows "
                                         f"< tmp_nsol={tmp.nsol}")
            for r in range(tmp.nsol):   # ssolve consumes rows 0..tmp_nsol-1
                fr_row, k_used = reconstruct_fit_row(
                    [t.rows[r] for t in tabs], primes)
                kmax = max(kmax, k_used)
                terms = []
                for (lid, kinid), c in zip(tmp.var, fr_row):
                    if c == 0:
                        continue
                    height_bits = max(height_bits, abs(c.numerator).bit_length(),
                                      c.denominator.bit_length())
                    gid = sch.intid[lid]
                    terms.append({
                        "integral": labels[gid],
                        "global_id": gid,
                        "coeff": str(c),
                        "monomial": {p: e for p, e in
                                     zip(params, kn.table[kinid]) if e},
                    })
                relations.append({"work": w, "config": k, "row": r,
                                  "n_terms": len(terms), "terms": terms})
            k += 1
    return relations, primes, kmax, height_bits


def export_relations(cfg: SearchFamily, works: Sequence[int],
                     labels: Sequence[str], out_path: str,
                     extra_meta: Optional[dict] = None,
                     prime_adder: Optional[Callable[[], str]] = None,
                     max_primes: int = 8) -> dict:
    """Read fit tables of all datanames, reconstruct exact rationals, and bank
    a human-readable JSON:  relation = sum_terms coeff * prod(params^exps) *
    G[label] == 0.  Returns the dict (also written to out_path).

    Growth semantics (upstream BLReconstructRelations While-loop,
    ReconstructRelation.wl:79-91): on ReconstructionOverflow, `prime_adder()`
    must fit the SAME works at ONE more prime (database + fitrel; see
    add_fit_prime) and append the new dataname to cfg.datanames; the pass is
    then retried over the enlarged set.  Without a prime_adder, or once
    len(cfg.datanames) == max_primes, the TYPED error propagates -- REFUSE,
    never emit.  Doc metadata reports the largest CRT prefix any row needed
    (primes_fit) and the held-out primes (primes_verify); the extra
    "crt_growth" key appears only when the set actually grew beyond the
    classic 2+1 layout, keeping earlier banked exports byte-identical.

    Reduced-analytic (nana) contract, FAIL-CLOSED at the producer: if ANY
    consumed database pool carries nana_pins.json (a pinned slice pool made
    by init_pinned_database), the export must declare and match the marker
    (extra_meta["reduced_analytic"], see NanaSpec.marker) or the TYPED
    NanaContractError is raised and nothing is emitted."""
    _gate_nana_export_marker(cfg, extra_meta)
    while True:
        try:
            relations, primes, kmax, height_bits = \
                _collect_relations(cfg, works, labels)
            break
        except ReconstructionOverflow as e:
            if prime_adder is None:
                raise
            if len(cfg.datanames) >= max_primes:
                raise ReconstructionOverflow(
                    f"UNVERIFIED at max prime count {max_primes}: {e}",
                    n_primes=e.n_primes, modulus_bits=e.modulus_bits) from e
            cfg.log(f"[export] {e}")
            n0 = len(cfg.datanames)
            new = prime_adder()
            if len(cfg.datanames) != n0 + 1:
                raise BladeGateError("prime_adder must append exactly one new "
                                     "dataname to cfg.datanames")
            cfg.log(f"[export] CRT prime set grown: +{new} "
                    f"({n0} -> {len(cfg.datanames)} primes); retrying")
    doc = {
        "family_job": cfg.job,
        "search_params": list(cfg.fam.search_params),
        "primes_fit": primes[:kmax],
        "primes_verify": primes[kmax:],
        "reconstruction": ("CRT(" + ",".join(f"p{i}" for i in range(kmax))
                           + ") + Wang; every entry verified mod "
                           + ",".join(f"p{i}" for i in range(kmax, len(primes)))),
        "n_relations": len(relations),
        "relations": relations,
    }
    if kmax > 2 or len(primes) != 3:
        doc["crt_growth"] = {
            "n_primes_fitted": len(primes),
            "crt_primes_used_max": kmax,
            "max_coeff_height_bits": height_bits,
            "joint_modulus_bits": sum(p.bit_length() for p in primes[:kmax]),
        }
    if extra_meta:
        doc.update(extra_meta)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(doc, f, indent=1)
    return doc


def add_fit_prime(cfg: SearchFamily, kc: Optional["KiraCoeffs"],
                  works: Sequence[int],
                  seed: int, reserved: Sequence[int] = (),
                  timeout: Optional[float] = None,
                  grower_factory: Optional[
                      Callable[[Sequence[int]], Callable[[int], None]]]
                  = None) -> str:
    """Fit `works` at ONE new fflow BIG_UINT prime (the export growth step;
    upstream: BLDatabase+BLFitRelations at the next pid,
    ReconstructRelation.wl:83-87).

    * chooses the smallest-index BIG_UINT prime that is neither already a fit
      prime nor in `reserved` (held-out ORACLE primes must never become fit
      primes -- pass them explicitly);
    * creates database/data_p<idx> with header/part/posi copied from prime 0
      and an empty pool, grows it to prime-0's physical pool size with the
      Kira grower -- or, when `grower_factory` is given (e.g. the raw-system
      fflow lane), with grower_factory(pids) and `kc` may be None (fresh
      deterministic points either way: growers seed per (seed, pid, idx,
      attempt) and the new pid differs from every existing pool);
    * runs fitrel for the works (artifact-gated), appends the dataname to
      cfg.datanames and returns it."""
    from .formats import big_uint_primes
    P = big_uint_primes()
    used = {read_ints(os.path.join(cfg.dbdir(pid), "red_common"))[0]
            for pid in range(len(cfg.datanames))}
    blocked = used | set(reserved)
    idx = next(i for i, p in enumerate(P) if p not in blocked)
    prime = P[idx]
    name = f"data_p{idx}"
    dbdir = os.path.join(cfg.searchdir, "database", name)
    if os.path.exists(dbdir):
        raise BladeGateError(f"add_fit_prime: {dbdir} already exists but its "
                             f"prime is not in cfg.datanames -- refusing")
    src = cfg.dbdir(0)
    common = read_ints(os.path.join(src, "red_common"))
    n_pool = _file_rows(src)
    os.makedirs(dbdir)
    shutil.copy(os.path.join(src, "red_part"), os.path.join(dbdir, "red_part"))
    shutil.copy(os.path.join(src, "red_posi"), os.path.join(dbdir, "red_posi"))
    write_fab_rows(os.path.join(dbdir, "red_common"),
                   [[prime] + list(common[1:-1]) + [0]])
    write_fab_rows(os.path.join(dbdir, "red_ps"), [])
    write_fab_rows(os.path.join(dbdir, "red_data"), [])
    cfg.datanames.append(name)
    new_pid = len(cfg.datanames) - 1
    t0 = time.time()
    if grower_factory is not None:
        grower_factory((new_pid,))(n_pool)
    else:
        make_kira_grower(cfg, kc, seed=seed, pids=(new_pid,))(n_pool)
    cfg.log(f"[add_fit_prime] {name} (prime {prime}): pool 0 -> {n_pool} "
            f"points in {time.time()-t0:.1f}s")
    for w in works:
        shutil.rmtree(os.path.join(cfg.workdir(w), "fit", name),
                      ignore_errors=True)
        t0 = time.time()
        run_fitrel_works(cfg, name, w, w + 1, timeout=timeout)
        cfg.log(f"[add_fit_prime] fitrel[{name}, w{w}] {time.time()-t0:.1f}s")
    return name


# --------------------------------------------------------------------------
# Kira-expression GF(p) machinery: database GROWER + independent held-out
# relation check (generalizes phase2 build_db_dbox2.py; no sympy cancel)
# --------------------------------------------------------------------------

class _GF:
    __slots__ = ("v", "p")
    def __init__(self, v, p): self.v = v % p; self.p = p
    def _c(self, o): return o if isinstance(o, _GF) else _GF(o, self.p)
    def __add__(self, o): o = self._c(o); return _GF(self.v + o.v, self.p)
    __radd__ = __add__
    def __sub__(self, o): o = self._c(o); return _GF(self.v - o.v, self.p)
    def __rsub__(self, o): o = self._c(o); return _GF(o.v - self.v, self.p)
    def __mul__(self, o): o = self._c(o); return _GF(self.v * o.v, self.p)
    __rmul__ = __mul__
    def __truediv__(self, o):
        o = self._c(o); return _GF(self.v * pow(o.v, -1, self.p), self.p)
    def __rtruediv__(self, o):
        o = self._c(o); return _GF(o.v * pow(self.v, -1, self.p), self.p)
    def __pow__(self, e):
        return _GF(pow(self.v, e, self.p) if e >= 0
                   else pow(pow(self.v, -1, self.p), -e, self.p), self.p)
    def __neg__(self): return _GF(-self.v, self.p)
    def __pos__(self): return self


class KiraCoeffs:
    """Parse a kira2math .m table: '<fam>[i1,...] -> sum_j <fam>[m...]*(coef)'.
    Coefficients are compiled once; evaluated in GF(p).  `env_fn(coords, p)`
    must return the eval namespace (e.g. d = 4-2*eps for (msq,t,eps))."""

    def __init__(self, path: str, family: str, targets: Sequence[Tuple[int, ...]],
                 masters: Sequence[Tuple[int, ...]],
                 env_fn: Callable[[Sequence[int], int], dict]):
        self.env_fn = env_fn
        self.targets = [tuple(t) for t in targets]
        self.masters = [tuple(m) for m in masters]
        tidx = {t: i for i, t in enumerate(self.targets)}
        midx = {m: j for j, m in enumerate(self.masters)}
        jb = re.compile(re.escape(family) + r"\[([0-9,\- ]+)\]")
        txt = open(path).read().strip()
        if not (txt.startswith("{") and txt.endswith("}")):
            raise BladeFormatError(f"{path}: not a kira2math list")
        chunks = [c.strip() for c in
                  re.split(r",\s*(?=" + re.escape(family) + r"\[)", txt[1:-1])
                  if c.strip()]
        self.cexpr: Dict[Tuple[int, int], object] = {}
        for c in chunks:
            lhs_txt, rhs_txt = c.split("->", 1)
            lhs = tuple(int(x) for x in jb.match(lhs_txt.strip()).group(1).split(","))
            if lhs not in tidx:
                continue
            i = tidx[lhs]
            rhs = re.sub(r"\s+", "", rhs_txt)
            for term in _split_terms(rhs):
                mm = jb.findall(term)
                if len(mm) != 1:
                    raise BladeFormatError(f"term with {len(mm)} master refs")
                mtup = tuple(int(x) for x in mm[0].split(","))
                if mtup not in midx:
                    raise BladeFormatError(f"unknown master {mtup}")
                coef = jb.sub("1", term).replace("^", "**")
                coef = re.sub(r"(?<!\*\*)\b(\d+)\b", r"G(\1)", coef)
                key = (i, midx[mtup])
                if key in self.cexpr:
                    raise BladeFormatError(f"duplicate coefficient {key}")
                self.cexpr[key] = compile(coef, f"<c{key}>", "eval")
        if not self.cexpr:
            raise BladeFormatError(f"{path}: no target rules matched")

    def value(self, i: int, j: int, coords: Sequence[int], p: int) -> int:
        code = self.cexpr.get((i, j))
        if code is None:
            return 0
        env = dict(self.env_fn(coords, p))
        env["G"] = lambda v, _p=p: _GF(v, _p)
        env["__builtins__"] = {}
        return eval(code, env).v

    def red_row(self, gid: int, coords: Sequence[int], p: int,
                nmaster: int) -> List[int]:
        """Reduction of global integral gid onto the masters, mod p.  Global
        order contract: targets first, masters last (database convention)."""
        nt = len(self.targets)
        if gid >= nt:
            return [1 if j == gid - nt else 0 for j in range(nmaster)]
        return [self.value(gid, j, coords, p) for j in range(nmaster)]


def _split_terms(rhs: str) -> List[str]:
    out, depth, cur = [], 0, []
    for ch in rhs:
        if ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth -= 1
        if ch in "+-" and depth == 0 and cur and cur[-1] not in "*/^(+-eE":
            out.append("".join(cur)); cur = [ch]
        else:
            cur.append(ch)
    if cur:
        out.append("".join(cur))
    return [t for t in out if t.strip()]


def make_kira_grower(cfg: SearchFamily, kc: KiraCoeffs, seed: int,
                     pids: Optional[Sequence[int]] = None
                     ) -> Callable[[int], None]:
    """Grower callback: extend red_ps/red_data of the selected prime pools to
    n points (default: all).  During escalation only prime 0 is consumed
    (BLRedG1 runs on dataname 0, wl:184); fit-time pools can be grown lazily
    to the bisected fitnumber.  Deterministic per (seed, pid, point index); a
    point where any stored coefficient value-vanishes is resampled (stored
    entries must be nonzero, FORMATS.md contract)."""
    import random as _random

    def grow(n_total: int) -> None:
        for pid in (pids if pids is not None else range(len(cfg.datanames))):
            dbdir = cfg.dbdir(pid)
            # read the PHYSICAL pool (red_common nps may be logically shrunk
            # by set_data_number -- Database.read would refuse the mismatch)
            common = read_ints(os.path.join(dbdir, "red_common"))
            prime, npara, nint, nmaster, nentry, _nps_log = common
            part = read_ints(os.path.join(dbdir, "red_part"))
            posi = read_ints(os.path.join(dbdir, "red_posi"))
            ps = formats._rows(read_ints(os.path.join(dbdir, "red_ps")),
                               npara, "red_ps")
            data = formats._rows(read_ints(os.path.join(dbdir, "red_data")),
                                 nentry, "red_data")
            if len(data) != len(ps):
                raise BladeGateError(f"{dbdir}: red_ps {len(ps)} rows != "
                                     f"red_data {len(data)} rows")
            phys = len(ps)
            if phys >= n_total:
                continue
            nt = len(kc.targets)
            keys = sorted(kc.cexpr)
            for idx in range(phys, n_total):
                for attempt in range(50):
                    rng = _random.Random(f"{seed}:{pid}:{idx}:{attempt}")
                    coord = [rng.randrange(1, prime) for _ in range(npara)]
                    vals = {}
                    ok = True
                    for (i, j) in keys:
                        v = kc.value(i, j, coord, prime)
                        if v == 0:
                            ok = False
                            break
                        vals[(i, j)] = v
                    if ok:
                        break
                else:
                    raise BladeGateError("grower: 50 resamples all value-vanish")
                row = []
                for m in range(nmaster):
                    for gi in posi[part[m]:part[m + 1]]:
                        row.append(vals[(gi, m)] if gi < nt else 1)
                ps.append(coord)
                data.append(row)
            write_fab_rows(os.path.join(dbdir, "red_ps"), ps)
            write_fab_rows(os.path.join(dbdir, "red_data"), data)
            common[-1] = len(ps)
            write_fab_rows(os.path.join(dbdir, "red_common"), [common])

    return grow


# --------------------------------------------------------------------------
# raw-system fflow machinery: database GROWER over a caller-supplied JSON
# sparse system, reduced per point by the fflow sparse eliminator (fflowcli
# learn/evalmany) -- the no-Kira database lane.  FflowSystem duck-types the
# held-out oracle interface (red_row), so heldout_relation_check runs
# unchanged against the system's own solve.
# --------------------------------------------------------------------------

class FflowSystem:
    """A caller-supplied sparse linear system in the fflow JSON grammar
    (header [neqs, nvars, npars, n_needed, [needed cols...]|0, nfiles,
    [files...]]; per-file equation data), reduced per point by
    `fflowcli learn` / `fflowcli evalmany`.

    `columns` maps system column id -> integral index tuple.  The sparse
    eliminator pivots columns LEFT TO RIGHT, so the masters must occupy the
    HIGHEST column ids (order the columns most complex first).

    Construction learns the system ONCE and gates the layout:
    * the learn payload is self-consistent and dense
      (nparsout == n_dep x n_indep; a sparse layout is a typed refusal);
    * the needed dependent rows cover every target;
    * the independent set EQUALS the masters -- with masters=None the master
      list is DERIVED from the learn payload (sorted ascending, the same
      deterministic convention the kira-table lane uses).

    Values follow the solver convention dep = sum_m c_m * indep_m; `value`
    and `red_row` expose them KiraCoeffs-style (per-point cache), so the
    held-out relation gate replays the system's OWN sparse solve at a fresh
    prime -- independent of the fit/export machinery, NOT of the input
    system.  Primes must be fflow BIG_UINT primes (the binary refuses
    others, by design)."""

    def __init__(self, path: str, columns: Sequence[Sequence[int]],
                 targets: Sequence[Tuple[int, ...]],
                 masters: Optional[Sequence[Tuple[int, ...]]] = None,
                 homog: bool = True, bin_dir: str = DEFAULT_BIN_DIR,
                 nthreads: int = 1, scratch: Optional[str] = None,
                 timeout: float = 600.0,
                 log: Callable[[str], None] = print):
        if not os.path.exists(path):
            raise BladeGateError(f"fflow system not found: {path}")
        self.path = os.path.abspath(path)
        self.homog = bool(homog)
        self.bin_dir = bin_dir
        self.nthreads = int(nthreads)
        self.timeout = float(timeout)
        self.log = log
        self._scratch = scratch
        self._neval = 0
        self._cache: Dict[Tuple[int, Tuple[int, ...]],
                          Dict[Tuple[int, int], int]] = {}
        self.fflowcli = os.path.join(bin_dir, "fflowcli")
        if not os.path.exists(self.fflowcli):
            raise BladeGateError(
                f"fflowcli not found at {self.fflowcli} -- the raw-system "
                f"lane reduces per point with the fflow sparse eliminator; "
                f"build the blade fork and point BLADE_BIN_DIR at its bin "
                f"dir")
        header = json.load(open(self.path))
        if not (isinstance(header, list) and len(header) == 7):
            raise BladeFormatError(
                f"{path}: not a 7-entry fflow sparse-system header "
                f"[neqs, nvars, npars, n_needed, needed, nfiles, files]")
        self.neqs, nvars, self.npara = (int(header[0]), int(header[1]),
                                        int(header[2]))
        self.columns = [tuple(int(x) for x in c) for c in columns]
        if len(self.columns) != nvars:
            raise BladeGateError(
                f"columns map has {len(self.columns)} entries but the "
                f"system header declares nvars={nvars}")
        sdir = os.path.dirname(self.path)
        for f in header[6]:
            full = f if os.path.isabs(f) else os.path.join(sdir, f)
            if not os.path.exists(full):
                raise BladeGateError(f"system data file not found: {full} "
                                     f"(referenced by {path})")
        self.targets = [tuple(t) for t in targets]
        self._learn(masters)

    # -- subprocess (rc recorded; the FFLOWCLI-ERROR marker + artifact
    #    byte-size/echo checks gate, per the standing rc-lies discipline)
    def _run(self, argv: List[str], stage: str) -> StageResult:
        res = _run_binary(argv, cwd=os.path.dirname(self.path),
                          timeout=self.timeout, stage=stage, log=self.log)
        if "FFLOWCLI-ERROR" in res.log:
            raise BladeGateError(f"{stage}: {res.log.strip()[:300]}")
        return res

    def _learn(self, masters: Optional[Sequence[Tuple[int, ...]]]) -> None:
        argv = [self.fflowcli, "learn", str(self.npara), self.path]
        if self.homog:
            argv.append("--homog")
        res = self._run(argv, "fflowcli-learn")
        lines = res.log.splitlines()

        def ids(prefix: str) -> List[int]:
            for ln in lines:
                if ln.startswith(prefix + " "):
                    v = [int(x) for x in ln.split()[1:]]
                    if not v or v[0] != len(v) - 1:
                        raise BladeGateError(
                            f"learn payload: malformed {prefix} line {ln!r}")
                    return v[1:]
            raise BladeGateError(f"learn payload: no {prefix!r} line "
                                 f"(log tail: {res.log[-200:]!r})")

        dep_cols, ind_cols = ids("depvars"), ids("indepvars")
        for ln in lines:
            if ln.startswith("sparse_output ") and ln.split()[1] != "0":
                raise BladeGateError(
                    "learn payload: sparse output layout unsupported by "
                    "this lane (dense dep x indep blocks expected)")
        npo = None
        for ln in lines:
            if ln.startswith("nparsout "):
                npo = int(ln.split()[1])
        if npo is None:
            raise BladeGateError("learn payload: no nparsout line")
        if npo != len(dep_cols) * len(ind_cols):
            raise BladeGateError(
                f"learn payload: nparsout {npo} != n_dep x n_indep = "
                f"{len(dep_cols)}x{len(ind_cols)} (dense layout expected; "
                f"for a homogeneous system keep homog=True)")
        if ind_cols != sorted(ind_cols):
            raise BladeGateError("learn payload: indepvars not ascending")
        bad = [c for c in dep_cols + ind_cols if c >= len(self.columns)]
        if bad:
            raise BladeGateError(f"learn payload: column id {bad[0]} outside "
                                 f"the columns map (len {len(self.columns)})")
        dep_tuples = [self.columns[c] for c in dep_cols]
        ind_tuples = [self.columns[c] for c in ind_cols]
        missing = [t for t in self.targets if t not in set(dep_tuples)]
        if missing:
            raise BladeGateError(
                f"learned needed rows do not cover {len(missing)} target(s) "
                f"(e.g. {missing[0]}) -- declare them in the system's "
                f"needed list and check the columns map")
        if masters is None:
            self.masters = sorted(set(ind_tuples))
        else:
            self.masters = [tuple(m) for m in masters]
            if set(self.masters) != set(ind_tuples):
                diff = sorted(set(self.masters) ^ set(ind_tuples))
                raise BladeGateError(
                    f"learned independent set != declared masters "
                    f"(symmetric difference {diff[:4]}) -- the eliminator's "
                    f"surviving independents define the masters; fix the "
                    f"master list or the column order (masters must sit at "
                    f"the highest column ids)")
        self.nparsout = npo
        tidx = {t: i for i, t in enumerate(self.targets)}
        midx = {m: j for j, m in enumerate(self.masters)}
        # output slot j -> (target index | None, master index): dep blocks in
        # payload depvars order, within a block indepvars ascending col id
        self._slots: List[Tuple[Optional[int], int]] = []
        for dt in dep_tuples:
            ti = tidx.get(dt)           # None: a needed row that is no target
            for it in ind_tuples:
                self._slots.append((ti, midx[it]))

    def _scratch_dir(self) -> str:
        if self._scratch is None:
            import tempfile
            self._scratch = tempfile.mkdtemp(prefix="fflow_system_")
        os.makedirs(self._scratch, exist_ok=True)
        return self._scratch

    def values_at(self, points: Sequence[Sequence[int]],
                  prime: int) -> List[Dict[Tuple[int, int], int]]:
        """Batched per-point reduction: one {(target_i, master_j): value}
        dict per point (dense over the learned target slots, zeros
        included).  Convention: value = coefficient in
        target = sum_j c_j * master_j."""
        if prime not in set(formats.big_uint_primes()):
            raise BladeGateError(f"prime {prime} is not an fflow BIG_UINT "
                                 f"prime -- fflowcli refuses others")
        npo, fs = self.nparsout, formats.flags_size(self.nparsout)
        mask = (1 << npo) - 1
        flags = [(mask >> (64 * k)) & ((1 << 64) - 1) for k in range(fs)]
        rows = []
        for c in points:
            if len(c) != self.npara:
                raise BladeGateError(f"point has {len(c)} coords; the "
                                     f"system has npars={self.npara}")
            rows.append([int(x) % prime for x in c] + [prime] + flags)
        sd = self._scratch_dir()
        k = self._neval
        self._neval += 1
        ppath = os.path.join(sd, f"points_{k}.fflow")
        epath = os.path.join(sd, f"eval_{k}.fflow")
        formats.PointsFile(n=self.npara + fs, rows=rows).write(ppath)
        argv = [self.fflowcli, "evalmany", str(self.npara), self.path,
                ppath, epath, "--nthreads", str(self.nthreads),
                "--expect-nout", str(npo)]
        if self.homog:
            argv.append("--homog")
        self._run(argv, f"fflowcli-evalmany[{len(points)}pts]")
        want = formats.EvalFile.expected_bytes(self.npara, npo, len(points))
        got = os.path.getsize(epath) if os.path.exists(epath) else -1
        if got != want:
            raise BladeGateError(f"evalmany: {epath} bytes {got} != "
                                 f"expected {want} (artifact gate, never rc)")
        ev = formats.EvalFile.read(epath)
        out = []
        for i, r in enumerate(ev.rows):
            if r[:self.npara] != rows[i][:self.npara] or \
                    r[self.npara] != prime:
                raise BladeGateError(f"evalmany: row {i} coords/prime echo "
                                     f"mismatch in {epath}")
            outs = r[self.npara + 1 + fs:]
            d = {}
            for j, (ti, mj) in enumerate(self._slots):
                if ti is not None:
                    d[(ti, mj)] = outs[j]
            out.append(d)
        os.unlink(ppath)
        os.unlink(epath)
        return out

    def value(self, i: int, j: int, coords: Sequence[int], p: int) -> int:
        key = (p, tuple(int(c) % p for c in coords))
        if key not in self._cache:
            if len(self._cache) > 64:
                self._cache.pop(next(iter(self._cache)))
            self._cache[key] = self.values_at([list(key[1])], p)[0]
        return self._cache[key].get((i, j), 0)

    def red_row(self, gid: int, coords: Sequence[int], p: int,
                nmaster: int) -> List[int]:
        """KiraCoeffs.red_row contract: reduction of global integral gid
        onto the masters mod p (targets first, masters last)."""
        nt = len(self.targets)
        if gid >= nt:
            return [1 if j == gid - nt else 0 for j in range(nmaster)]
        return [self.value(gid, j, coords, p) for j in range(nmaster)]

    def measure_pattern(self, prime: int, npoints: int,
                        seed: int) -> set:
        """Measured (target, master) nonzero pattern at deterministic probe
        points -- the SAME rng stream the pid-0 grower banks first, so the
        probed pattern and the first banked points agree.  Stored = nonzero
        at >= 1 probe point; a stored pair that value-vanishes at a later
        point is resampled by the grower (stored entries must be nonzero,
        red_* contract)."""
        import random as _random
        pts = []
        for idx in range(npoints):
            rng = _random.Random(f"{seed}:0:{idx}:0")
            pts.append([rng.randrange(1, prime) for _ in range(self.npara)])
        stored = set()
        for d in self.values_at(pts, prime):
            stored |= {k for k, v in d.items() if v != 0}
        if not stored:
            raise BladeGateError("pattern probe: every (target, master) "
                                 "coefficient vanished at the probe points")
        return stored


def make_fflow_grower(cfg: SearchFamily, fsys: FflowSystem, seed: int,
                      pids: Optional[Sequence[int]] = None
                      ) -> Callable[[int], None]:
    """Grower callback over a raw fflow system (the make_kira_grower
    contract): extend red_ps/red_data of the selected prime pools to n
    points via batched evalmany.  Deterministic per (seed, pid, point
    index); a point where any STORED coefficient value-vanishes is
    resampled; a nonzero value at a pair OUTSIDE the stored pattern is a
    typed refusal (red_part/red_posi are structural -- the sparsity pattern
    must be point-uniform)."""
    import random as _random

    def grow(n_total: int) -> None:
        for pid in (pids if pids is not None else range(len(cfg.datanames))):
            dbdir = cfg.dbdir(pid)
            common = read_ints(os.path.join(dbdir, "red_common"))
            prime, npara, nint, nmaster, nentry, _nps_log = common
            if npara != fsys.npara:
                raise BladeGateError(f"{dbdir}: red_common npara={npara} != "
                                     f"system npars={fsys.npara}")
            part = read_ints(os.path.join(dbdir, "red_part"))
            posi = read_ints(os.path.join(dbdir, "red_posi"))
            ps = formats._rows(read_ints(os.path.join(dbdir, "red_ps")),
                               npara, "red_ps")
            data = formats._rows(read_ints(os.path.join(dbdir, "red_data")),
                                 nentry, "red_data")
            if len(data) != len(ps):
                raise BladeGateError(f"{dbdir}: red_ps {len(ps)} rows != "
                                     f"red_data {len(data)} rows")
            phys = len(ps)
            if phys >= n_total:
                continue
            nt = len(fsys.targets)
            stored = set()
            for m in range(nmaster):
                for gi in posi[part[m]:part[m + 1]]:
                    if gi < nt:
                        stored.add((gi, m))
            pending = {idx: 0 for idx in range(phys, n_total)}
            accepted: Dict[int, Tuple[List[int],
                                      Dict[Tuple[int, int], int]]] = {}
            while pending:
                batch = sorted(pending)
                coords = []
                for idx in batch:
                    rng = _random.Random(f"{seed}:{pid}:{idx}:{pending[idx]}")
                    coords.append([rng.randrange(1, prime)
                                   for _ in range(npara)])
                for idx, coord, d in zip(batch, coords,
                                         fsys.values_at(coords, prime)):
                    extra = sorted(k for k, v in d.items()
                                   if v != 0 and k not in stored)
                    if extra:
                        raise BladeGateError(
                            f"grower: nonzero value at (target,master) "
                            f"{extra[:3]} outside the stored pattern -- "
                            f"the system's sparsity is not point-uniform "
                            f"({dbdir})")
                    if any(d.get(k, 0) == 0 for k in stored):
                        pending[idx] += 1
                        if pending[idx] >= 50:
                            raise BladeGateError(
                                "grower: 50 resamples all value-vanish")
                        continue
                    accepted[idx] = (coord, d)
                    del pending[idx]
            for idx in sorted(accepted):
                coord, d = accepted[idx]
                row = []
                for m in range(nmaster):
                    for gi in posi[part[m]:part[m + 1]]:
                        row.append(d[(gi, m)] if gi < nt else 1)
                ps.append(coord)
                data.append(row)
            write_fab_rows(os.path.join(dbdir, "red_ps"), ps)
            write_fab_rows(os.path.join(dbdir, "red_data"), data)
            common[-1] = len(ps)
            write_fab_rows(os.path.join(dbdir, "red_common"), [common])

    return grow


def heldout_relation_check(doc: dict, kc, labels: Sequence[str],
                           nmaster: int, held_prime: int, npoints: int,
                           seed: int, log: Callable[[str], None] = print,
                           pinned_residues: Optional[Dict[int, int]] = None
                           ) -> Tuple[int, int]:
    """INDEPENDENT oracle gate: every exported rational relation must vanish
    identically when each integral is replaced by its reduction row from the
    oracle `kc` (KiraCoeffs, or any red_row duck-type such as FflowSystem),
    checked mod a prime NEVER used in search/fit/verify, at fresh points.
    Plus a synthetic-truth control (coeff+1 must be caught).  Returns
    (checks, fails).

    pinned_residues (reduced-analytic mode): {S-position: residue mod
    held_prime}.  Sample points are random ONLY in the analytic parameters;
    the pinned coordinates are fixed to the given residues -- the relations
    of a nana export hold ONLY on that slice (annihilation gate:
    relations at the pins must annihilate the banked reduction evaluated at
    the same pins).  Default None = classic full-symbolic behavior,
    byte-identical to before (the rng stream is unchanged)."""
    import random as _random
    rng = _random.Random(seed)
    label_gid = {lab: i for i, lab in enumerate(labels)}
    params = doc["search_params"]
    pts = [[rng.randrange(1, held_prime) for _ in params]
           for _ in range(npoints)]
    if pinned_residues:
        for i, r in pinned_residues.items():
            if not (0 <= i < len(params)):
                raise BladeGateError(f"pinned_residues position {i} outside "
                                     f"0..{len(params) - 1}")
            for c in pts:
                c[i] = r % held_prime
    checks = fails = 0

    def residuals(rel, mutate_first=False):
        nonlocal checks, fails
        bad = 0
        for coords in pts:
            red_cache = {}
            acc = [0] * nmaster
            for it, term in enumerate(rel["terms"]):
                gid = label_gid[term["integral"]]
                if gid not in red_cache:
                    red_cache[gid] = kc.red_row(gid, coords, held_prime, nmaster)
                mono = 1
                for pnam, e in term["monomial"].items():
                    mono = mono * pow(coords[params.index(pnam)], e,
                                      held_prime) % held_prime
                c = Fraction(term["coeff"])
                cv = c.numerator * pow(c.denominator, -1, held_prime) % held_prime
                if mutate_first and it == 0:
                    cv = (cv + 1) % held_prime
                w = cv * mono % held_prime
                for j in range(nmaster):
                    if red_cache[gid][j]:
                        acc[j] = (acc[j] + w * red_cache[gid][j]) % held_prime
            if any(acc):
                bad += 1
        return bad

    for rel in doc["relations"]:
        bad = residuals(rel)
        checks += npoints
        fails += bad
        if bad:
            log(f"[heldout] relation w{rel['work']}/c{rel['config']}/r{rel['row']}"
                f": {bad}/{npoints} residuals NONZERO")
    # synthetic-truth control on the first relation
    ctrl_bad = residuals(doc["relations"][0], mutate_first=True)
    if ctrl_bad != npoints:
        raise BladeGateError(f"heldout control: mutated coeff detected only "
                             f"{ctrl_bad}/{npoints} -- oracle has no teeth")
    log(f"[heldout] {checks} residual checks, {fails} failures; "
        f"control caught {ctrl_bad}/{npoints}")
    return checks, fails


# --------------------------------------------------------------------------
# top-level convenience
# --------------------------------------------------------------------------

def blade_search(cfg: SearchFamily, works: Sequence[int], labels: Sequence[str],
                 out_json: str, bisect: bool = True,
                 extra_meta: Optional[dict] = None
                 ) -> Tuple[List[WorkOutcome], Optional[dict]]:
    """Escalate every work; for the CLOSED subset: fitrel at every prime,
    optional fitnumber bisection (prime 0), rational export.  Returns
    (outcomes, relations_doc_or_None)."""
    outcomes = [escalate_work(cfg, w) for w in works]
    closed = [o.work for o in outcomes if o.closed]
    if not closed:
        return outcomes, None
    for pid in range(len(cfg.datanames)):
        for w in closed:
            shutil.rmtree(os.path.join(cfg.workdir(w), "fit",
                                       cfg.datanames[pid]), ignore_errors=True)
            run_fitrel_works(cfg, cfg.datanames[pid], w, w + 1)
    doc = export_relations(cfg, closed, labels, out_json, extra_meta)
    if bisect:
        fn = fitnumber_bisection(cfg, closed, pid=0)
        doc["fitnumber"] = fn
        with open(out_json, "w") as f:
            json.dump(doc, f, indent=1)
    return outcomes, doc


# ==========================================================================
# REDUCED-ANALYTIC PARAMETER SEARCH ("nana" mode) -- upstream BLAdaptiveSearch
# semantics: search relations analytic in nana < n_params parameters, with
# the remaining ("conservative") parameters numerically PINNED per instance.
#
# Upstream anatomy (line cites):
#   * nana = #analytic params; the BT system at level nana is searched with
#     VariableGroup = BLSearchParameter[[-nana;;]] -- analytic params are the
#     SUFFIX of the search-parameter order (SearchOptions.wl:229, enforced by
#     CheckSearchOptions, SearchOptions.wl:192-193); the pinned params are
#     the PREFIX: $Conservative = BLSearchParameter[[1;;ncons]]
#     (SolveSemiBL.wl:409; during search: SearchRelation.wl:174).
#   * pinning mechanism: the numeric IBP database freezes the conservative
#     coordinates to ONE shared value across the whole pool
#     (Database.wl:131-141: pslist[[All,cons]] = oldps[[1,cons]] / fixps),
#     and records their positions in red_cons (Database.wl:114-117).
#   * instance-count logic: consumer sample points are grouped by their
#     conservative coordinates -- GatherBy[pts, #[[1;;ncons]]&]
#     (gatherSamplePoints, SolveSemiBL.wl:296-304; the SAME grouping drives
#     the timeEstimate instance count, AdaptiveSearch.wl:69-73).  N search/
#     fit instances = number of DISTINCT pinned tuples; per instance
#     upstream re-fits the SEARCHED template at that pin:
#     $Conservative=...[[1;;ncons]]; BLDatabase[name, fitnumber, pid, pin];
#     BLFitRelations[name, job]  (blEvaluateAndDump, SolveSemiBL.wl:404-410),
#     amortized: buckets with <= fitnumber points go through plain numeric
#     IBP instead (SolveSemiBL.wl:399-401; AdaptiveSearch.wl:73).
#   * nana selection: BLAdaptiveSearch escalates nana = 1,2,... and picks the
#     wall-time optimum (AdaptiveSearch.wl:39-63); this port takes nana from
#     the caller (measured timed-pilot ladders replace the BLTimeIBP model) and
#     chooses WHICH params stay analytic via choose_analytic_params.
#
# FAIL-CLOSED consumer contract: relations fitted on a pinned pool are valid
# ONLY at the pins.  Every pinned pool carries nana_pins.json; the export
# must carry the matching "reduced_analytic" marker (typed NanaContractError
# otherwise); probe.BTForm and semibl.SemiBLJob REFUSE such artifacts for
# full-symbolic use.  Per-point fitting in the dropped parameters is the
# consumer's job: plan_pin_instances + fit_at_pin implement the upstream
# per-pin refit step.
# ==========================================================================

NANA_PINS_FILE = "nana_pins.json"


class NanaContractError(BladeGateError):
    """Typed refusal: a pinned (reduced-analytic) artifact was about to be
    used, or produced, as if it were full-symbolic.  Nothing is emitted."""


@dataclass(frozen=True)
class NanaSpec:
    """Which parameters stay analytic and where the rest are pinned.

    search_params : full BLSearchParameter order (S) of the family.
    analytic      : subset kept analytic, stored in S order.
    pins          : dropped param -> EXACT rational pin value (per search
                    instance; the GF(p) image is prime-dependent).
    origin        : provenance of the analytic choice ("override" |
                    "scores" | "positional-default")."""
    search_params: Tuple[str, ...]
    analytic: Tuple[str, ...]
    pins: Dict[str, Fraction] = field(default_factory=dict)
    origin: str = "override"

    def __post_init__(self):
        sp = self.search_params
        if len(set(sp)) != len(sp):
            raise BladeGateError(f"NanaSpec: duplicate search params {sp}")
        if not set(self.analytic) <= set(sp):
            raise BladeGateError(f"NanaSpec: analytic {self.analytic} not a "
                                 f"subset of {sp}")
        if len(self.analytic) < 1:
            raise BladeGateError("NanaSpec: nana must be >= 1 (at least one "
                                 "analytic parameter)")
        in_order = tuple(p for p in sp if p in set(self.analytic))
        object.__setattr__(self, "analytic", in_order)
        dropped = set(sp) - set(in_order)
        if set(self.pins) != dropped:
            raise BladeGateError(f"NanaSpec: pins keys {sorted(self.pins)} != "
                                 f"dropped params {sorted(dropped)}")
        for p, v in self.pins.items():
            if not isinstance(v, Fraction) or v == 0:
                raise BladeGateError(f"NanaSpec: pin {p}={v!r} must be a "
                                     f"NONZERO Fraction (0 is a degenerate "
                                     f"slice; see degenerate-slice footgun)")

    # ---- derived ----------------------------------------------------------
    @property
    def nana(self) -> int:
        return len(self.analytic)

    @property
    def ncons(self) -> int:
        return len(self.search_params) - self.nana

    @property
    def pinned(self) -> Tuple[str, ...]:
        return tuple(p for p in self.search_params if p not in set(self.analytic))

    @property
    def analytic_idx(self) -> List[int]:
        return [self.search_params.index(p) for p in self.analytic]

    @property
    def pinned_idx(self) -> List[int]:
        return [self.search_params.index(p) for p in self.pinned]

    def pin_residues(self, prime: int) -> Dict[int, int]:
        """{S-position: GF(prime) image of the pin} -- Database.wl fixps."""
        out = {}
        for p in self.pinned:
            fr = self.pins[p]
            if fr.denominator % prime == 0:
                raise BladeGateError(f"pin {p}={fr}: denominator divisible by "
                                     f"prime {prime}")
            out[self.search_params.index(p)] = (
                fr.numerator * pow(fr.denominator, -1, prime) % prime)
        return out

    def marker(self) -> dict:
        """The export-level FAIL-CLOSED marker (doc['reduced_analytic'])."""
        return {
            "nana": self.nana,
            "analytic_params": list(self.analytic),
            "pinned_params": {p: str(self.pins[p]) for p in self.pinned},
            "origin": self.origin,
            "valid_only_at_pins": True,
            "full_symbolic_use": (
                "FORBIDDEN -- these relations hold ONLY on the pinned slice; "
                "per-point fitting in the pinned parameters (one instance "
                "per distinct pin tuple, SolveSemiBL.wl:296-304,404-410; "
                "plan_pin_instances/fit_at_pin) is the consumer's job"),
        }


def choose_analytic_params(search_params: Sequence[str], nana: int,
                           scores: Optional[Dict[str, float]] = None,
                           override: Optional[Sequence[str]] = None
                           ) -> Tuple[Tuple[str, ...], str]:
    """Pick WHICH nana parameters stay analytic.  Returns (analytic, origin).

    Precedence:
      1. override -- explicit caller choice (validated).
      2. scores   -- per-parameter degree/complexity heuristic: PIN the
         highest-score (hardest) directions, keep the nana lowest analytic.
         This is the port's use of upstream's own complexity measure -- the
         maximal power of single-variable-analytic relations that ranks
         variable groups simple->complex (SearchRelation.wl:288-296,313-326;
         adaptive.py probes measure exactly these powerlimits).  Mass-
         dimension scoring (dimensionless eps -> 0, dimensionful params ->
         |dim|) is the cheap stand-in when no probes were run.  Ties break
         toward the upstream positional default (later S position stays
         analytic).
      3. neither  -- upstream positional default: analytic = the SUFFIX
         BLSearchParameter[[-nana;;]] (SearchOptions.wl:229), pinned = the
         prefix (SolveSemiBL.wl:409)."""
    sp = list(search_params)
    if not (1 <= nana <= len(sp)):
        raise BladeGateError(f"nana={nana} outside 1..{len(sp)}")
    if override is not None:
        ov = list(override)
        if len(ov) != nana or not set(ov) <= set(sp) or len(set(ov)) != nana:
            raise BladeGateError(f"override {ov} is not a {nana}-subset of {sp}")
        return tuple(p for p in sp if p in set(ov)), "override"
    if scores is not None:
        missing = set(sp) - set(scores)
        if missing:
            raise BladeGateError(f"scores missing for {sorted(missing)}")
        ranked = sorted(range(len(sp)), key=lambda i: (scores[sp[i]], -i))
        keep = set(ranked[:nana])
        return tuple(p for i, p in enumerate(sp) if i in keep), "scores"
    return tuple(sp[-nana:]), "positional-default"


def default_pins(params: Sequence[str], seed: int) -> Dict[str, Fraction]:
    """Deterministic generic rational pins (small height, never 0/±1, distinct
    per parameter).  Prefer EXPLICIT pins in production runs; these exist so
    scripted gates are reproducible.

    MEASURED guidance (dbox-w0 gate 3): the slice coefficients
    contain the pin raised to high powers, so EXPORT prime count scales with
    (effective pin degree) x (pin bit height).  msq=41/23 (~10 bits/power)
    drove w0's export past 24 primes; a low-height pin (2, 3, 3/2) costs
    ~1 bit/power and reconstructs in a handful of primes.  Pick the smallest
    pin that is not a degenerate slice for the family."""
    import random as _random
    small_primes = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
    out: Dict[str, Fraction] = {}
    for p in params:
        rng = _random.Random(f"{seed}:nana-pin:{p}")
        while True:
            num, den = rng.sample(small_primes, 2)
            fr = Fraction(num, den)
            if fr not in out.values() and fr != 1:
                out[p] = fr
                break
    return out


# --------------------------------------------------------------------------
# pinned database pools (Database.wl:114-141 semantics)
# --------------------------------------------------------------------------

def read_nana_pins(dbdir: str) -> Optional[dict]:
    path = os.path.join(dbdir, NANA_PINS_FILE)
    if not os.path.exists(path):
        return None
    return json.load(open(path))


def init_pinned_database(searchdir: str, dataname: str, src_dataname: str,
                         spec: NanaSpec, prime: Optional[int] = None) -> str:
    """Create database/<dataname> with header/part/posi copied from
    src_dataname, an EMPTY pool, the red_cons record (0-based S-positions of
    the pinned params -- port encoding of Database.wl:114-117's Put[posi])
    and the nana_pins.json marker.  Refuses to overwrite."""
    dbdir = os.path.join(searchdir, "database", dataname)
    src = os.path.join(searchdir, "database", src_dataname)
    if os.path.exists(dbdir):
        raise BladeGateError(f"init_pinned_database: {dbdir} exists -- refusing")
    common = read_ints(os.path.join(src, "red_common"))
    if prime is not None:
        common = [prime] + list(common[1:])
    if common[1] != len(spec.search_params):
        raise BladeGateError(f"{src}/red_common: npara {common[1]} != "
                             f"len(search_params) {len(spec.search_params)}")
    residues = spec.pin_residues(common[0])
    os.makedirs(dbdir)
    shutil.copy(os.path.join(src, "red_part"), os.path.join(dbdir, "red_part"))
    shutil.copy(os.path.join(src, "red_posi"), os.path.join(dbdir, "red_posi"))
    write_fab_rows(os.path.join(dbdir, "red_common"),
                   [list(common[:-1]) + [0]])
    write_fab_rows(os.path.join(dbdir, "red_ps"), [])
    write_fab_rows(os.path.join(dbdir, "red_data"), [])
    write_fab_rows(os.path.join(dbdir, "red_cons"), [spec.pinned_idx])
    with open(os.path.join(dbdir, NANA_PINS_FILE), "w") as f:
        json.dump({
            "search_params": list(spec.search_params),
            "analytic_params": list(spec.analytic),
            "pins": {p: str(spec.pins[p]) for p in spec.pinned},
            "pinned_positions": spec.pinned_idx,
            "prime": common[0],
            "residues": {p: residues[spec.search_params.index(p)]
                         for p in spec.pinned},
        }, f, indent=1)
    return dbdir


def _gate_pinned_pool(dbdir: str, spec: NanaSpec) -> Dict[int, int]:
    """Verify dbdir's nana_pins.json matches spec at the pool's prime; return
    the pin residues by S-position.  Typed refusal on any mismatch."""
    info = read_nana_pins(dbdir)
    if info is None:
        raise NanaContractError(f"{dbdir}: no {NANA_PINS_FILE} -- refusing to "
                                f"treat a generic pool as pinned")
    prime = read_ints(os.path.join(dbdir, "red_common"))[0]
    residues = spec.pin_residues(prime)
    want = {
        "search_params": list(spec.search_params),
        "analytic_params": list(spec.analytic),
        "pins": {p: str(spec.pins[p]) for p in spec.pinned},
        "pinned_positions": spec.pinned_idx,
        "prime": prime,
        "residues": {p: residues[spec.search_params.index(p)]
                     for p in spec.pinned},
    }
    got = {k: info.get(k) for k in want}
    if got != want:
        raise NanaContractError(f"{dbdir}: {NANA_PINS_FILE} mismatch\n"
                                f"  banked: {got}\n  spec:   {want}")
    return residues


def make_pinned_grower(cfg: SearchFamily, kc: KiraCoeffs, seed: int,
                       spec: NanaSpec,
                       pids: Optional[Sequence[int]] = None
                       ) -> Callable[[int], None]:
    """Grower with upstream fixps semantics (Database.wl:131-141): the pinned
    coordinates are FIXED to the GF(p) images of spec.pins in every point of
    every selected pool; only the analytic coordinates are sampled.  On a
    value-vanishing coefficient the ANALYTIC coordinates are resampled while
    the pins stay put (upstream keeps the shared conservative row).  Gated
    against each pool's nana_pins.json (never grows a generic pool)."""
    import random as _random
    free_idx = spec.analytic_idx

    def grow(n_total: int) -> None:
        for pid in (pids if pids is not None else range(len(cfg.datanames))):
            dbdir = cfg.dbdir(pid)
            common = read_ints(os.path.join(dbdir, "red_common"))
            prime, npara, nint, nmaster, nentry, _nps_log = common
            residues = _gate_pinned_pool(dbdir, spec)
            part = read_ints(os.path.join(dbdir, "red_part"))
            posi = read_ints(os.path.join(dbdir, "red_posi"))
            ps = formats._rows(read_ints(os.path.join(dbdir, "red_ps")),
                               npara, "red_ps")
            data = formats._rows(read_ints(os.path.join(dbdir, "red_data")),
                                 nentry, "red_data")
            if len(data) != len(ps):
                raise BladeGateError(f"{dbdir}: red_ps {len(ps)} rows != "
                                     f"red_data {len(data)} rows")
            base = [0] * npara
            for i, r in residues.items():
                base[i] = r
            for row in ps[:1]:                    # existing pool must be pinned
                for i, r in residues.items():
                    if row[i] != r:
                        raise NanaContractError(
                            f"{dbdir}: existing point has coord[{i}]={row[i]} "
                            f"!= pin residue {r} -- pool is not this slice")
            phys = len(ps)
            if phys >= n_total:
                continue
            nt = len(kc.targets)
            keys = sorted(kc.cexpr)
            for idx in range(phys, n_total):
                for attempt in range(200):
                    rng = _random.Random(f"{seed}:{pid}:{idx}:{attempt}")
                    coord = list(base)
                    for i in free_idx:
                        coord[i] = rng.randrange(1, prime)
                    vals = {}
                    ok = True
                    for (i, j) in keys:
                        v = kc.value(i, j, coord, prime)
                        if v == 0:
                            ok = False
                            break
                        vals[(i, j)] = v
                    if ok:
                        break
                else:
                    raise BladeGateError(
                        f"pinned grower: 200 resamples of the analytic "
                        f"coordinates all value-vanish at pins "
                        f"{ {p: str(spec.pins[p]) for p in spec.pinned} } -- "
                        f"degenerate slice (footgun); choose different pins")
                row = []
                for m in range(nmaster):
                    for gi in posi[part[m]:part[m + 1]]:
                        row.append(vals[(gi, m)] if gi < nt else 1)
                ps.append(coord)
                data.append(row)
            write_fab_rows(os.path.join(dbdir, "red_ps"), ps)
            write_fab_rows(os.path.join(dbdir, "red_data"), data)
            common[-1] = len(ps)
            write_fab_rows(os.path.join(dbdir, "red_common"), [common])

    return grow


# --------------------------------------------------------------------------
# restricted ansatz + config assembly
# --------------------------------------------------------------------------

def nana_family(fam: FamilyAnsatz, spec: NanaSpec) -> FamilyAnsatz:
    """Ansatz recipe restricted to the analytic parameters.  kin_table keeps
    the FULL search-parameter columns (ansatz O4: parameters outside the
    variable groups get all-zero columns -- exactly upstream's conservative
    encoding, PolynomialAnsatz.wl:62-72), so databases/points stay S-ordered.
    Feature off (ncons == 0): returns fam UNCHANGED -- the nana path is
    byte-identical to the classic search by construction."""
    if spec.ncons == 0:
        return fam
    if tuple(fam.search_params) != tuple(spec.search_params):
        raise BladeGateError(f"nana_family: fam.search_params "
                             f"{fam.search_params} != spec {spec.search_params}")
    keep = set(spec.analytic)
    groups, weights, intmode = [], [], []
    for g, w, m in zip(fam.var_groups, fam.var_weights, fam.intmode):
        gg = [v for v in g if v in keep]
        ww = [wi for v, wi in zip(g, w) if v in keep]
        if gg:
            groups.append(gg)
            weights.append(ww)
            intmode.append(m)
    if not groups:
        raise BladeGateError("nana_family: no analytic parameter survives")
    if fam.cut and len(groups) != len(fam.var_groups):
        raise BladeGateError(
            "nana_family: dropping a whole variable group of a multi-group "
            "(cut != []) recipe changes the composite-level radix -- "
            "restructure the recipe explicitly instead")
    return FamilyAnsatz(search_params=tuple(fam.search_params),
                        var_groups=groups, var_weights=weights,
                        intmode=intmode, cut=list(fam.cut), nmax=fam.nmax)


def make_nana_config(cfg: SearchFamily, spec: NanaSpec, kc: KiraCoeffs,
                     seed: int, works: Sequence[int],
                     job: Optional[str] = None,
                     datanames: Optional[List[str]] = None,
                     src_datanames: Optional[Sequence[str]] = None
                     ) -> SearchFamily:
    """Derive a reduced-analytic SearchFamily from a full one:
      * feature OFF (spec.ncons == 0): returns cfg UNCHANGED (gate-1 path);
      * else: fresh job dir (schemes byte-copied from cfg's works -- the
        scheme is parameter-independent), restricted ansatz, NEW pinned pools
        (one per src dataname, same primes) and the pinned grower.
    The returned config writes NOTHING into cfg's job/kinematics/databases."""
    if spec.ncons == 0:
        return cfg
    job = job or f"{cfg.job}_nana{spec.nana}"
    jobdir = os.path.join(cfg.searchdir, job)
    if os.path.exists(jobdir):
        raise BladeGateError(f"make_nana_config: {jobdir} exists -- refusing")
    for w in works:
        Scheme.read(cfg.workdir(w)).write(os.path.join(jobdir, str(w)))
    srcs = list(src_datanames if src_datanames is not None else cfg.datanames)
    names = list(datanames if datanames is not None
                 else [f"{job}_p{i}" for i in range(len(srcs))])
    if len(names) != len(srcs):
        raise BladeGateError("make_nana_config: datanames/src_datanames "
                             "length mismatch")
    for name, src in zip(names, srcs):
        init_pinned_database(cfg.searchdir, name, src, spec)
    ncfg = _dc_replace(cfg, job=job, fam=nana_family(cfg.fam, spec),
                       datanames=names, grower=None)
    ncfg.grower = make_pinned_grower(ncfg, kc, seed=seed, spec=spec)
    return ncfg


# --------------------------------------------------------------------------
# producer-side FAIL-CLOSED export marker gate
# --------------------------------------------------------------------------

def _gate_nana_export_marker(cfg: SearchFamily,
                             extra_meta: Optional[dict]) -> None:
    """If any consumed pool is pinned, the export MUST declare the matching
    reduced_analytic marker; if none is pinned, the marker must be absent."""
    pinned_pools = {}
    for pid in range(len(cfg.datanames)):
        info = read_nana_pins(cfg.dbdir(pid))
        if info is not None:
            pinned_pools[cfg.datanames[pid]] = info
    marker = (extra_meta or {}).get("reduced_analytic")
    if pinned_pools and marker is None:
        raise NanaContractError(
            f"export_relations: pools {sorted(pinned_pools)} are PINNED "
            f"(reduced-analytic slices) but extra_meta carries no "
            f"'reduced_analytic' marker -- refusing to emit relations that "
            f"would masquerade as full-symbolic")
    if marker is not None:
        if not pinned_pools:
            raise NanaContractError(
                "export_relations: 'reduced_analytic' marker supplied but no "
                "consumed pool is pinned -- marker would be a lie")
        for name, info in pinned_pools.items():
            if (info.get("pins") != marker.get("pinned_params")
                    or info.get("analytic_params") != marker.get("analytic_params")):
                raise NanaContractError(
                    f"export_relations: pool {name} pins {info.get('pins')} / "
                    f"analytic {info.get('analytic_params')} do not match the "
                    f"marker {marker.get('pinned_params')} / "
                    f"{marker.get('analytic_params')}")


# --------------------------------------------------------------------------
# per-pin consumer instances (SolveSemiBL.wl:296-304, 404-410)
# --------------------------------------------------------------------------

def plan_pin_instances(points: Sequence[Sequence[int]],
                       pinned_positions: Sequence[int]
                       ) -> List[Tuple[Tuple[int, ...], List[int]]]:
    """Upstream instance-count logic: group consumer sample points by their
    pinned coordinates (GatherBy[ptsp, #[[1;;ncons]]&], SolveSemiBL.wl:302;
    the identical grouping drives the timeEstimate instance count,
    AdaptiveSearch.wl:70).  Returns [(pin_tuple, [point row ids])] in
    first-occurrence order; N instances = len(result).  One pinned database +
    template fit per group (fit_at_pin); upstream amortization: a group with
    <= fitnumber points is cheaper by plain numeric IBP (SolveSemiBL.wl:
    399-401, AdaptiveSearch.wl:73) -- that routing stays with the caller."""
    order: List[Tuple[int, ...]] = []
    groups: Dict[Tuple[int, ...], List[int]] = {}
    for r, pt in enumerate(points):
        key = tuple(pt[i] for i in pinned_positions)
        if key not in groups:
            groups[key] = []
            order.append(key)
        groups[key].append(r)
    return [(k, groups[k]) for k in order]


def fit_at_pin(cfg: SearchFamily, spec: NanaSpec, kc: KiraCoeffs,
               works: Sequence[int], seed: int,
               dataname: Optional[str] = None, src_pid: int = 0,
               npoints: Optional[int] = None,
               timeout: Optional[float] = None) -> str:
    """ONE consumer instance: re-fit the already-SEARCHED template of `cfg`
    (a nana config whose works closed) at spec's pin values -- the port of
    blEvaluateAndDump's per-bucket refit ($Conservative=...; BLDatabase[name,
    fitnumber, pid, pin]; BLFitRelations[name, job], SolveSemiBL.wl:404-410).
    Creates database/<dataname> pinned at spec.pins (prime of src_pid's
    pool), grows it to `npoints` (default: the job's banked fitnumber, else
    the current data number of src_pid), runs fitrel for `works` (artifact-
    gated) and returns the dataname.  The search template is NOT re-searched:
    that is the whole point of the amortization."""
    src = cfg.datanames[src_pid]
    dataname = dataname or f"{cfg.job}_pin_" + "_".join(
        f"{p}{spec.pins[p].numerator}d{spec.pins[p].denominator}"
        for p in spec.pinned)
    if npoints is None:
        fnp = os.path.join(cfg.jobdir, "fitnumber")
        npoints = (int(open(fnp).read().strip()) if os.path.exists(fnp)
                   else data_number(cfg.dbdir(src_pid)))
    init_pinned_database(cfg.searchdir, dataname, src, spec)
    tmp_names = cfg.datanames
    try:
        cfg.datanames = list(tmp_names) + [dataname]
        pid_new = len(cfg.datanames) - 1
        make_pinned_grower(cfg, kc, seed=seed, spec=spec,
                           pids=(pid_new,))(npoints)
        for w in works:
            shutil.rmtree(os.path.join(cfg.workdir(w), "fit", dataname),
                          ignore_errors=True)
            run_fitrel_works(cfg, dataname, w, w + 1, timeout=timeout)
    finally:
        cfg.datanames = tmp_names
    return dataname


def add_fit_prime_pinned(cfg: SearchFamily, kc: KiraCoeffs, spec: NanaSpec,
                         works: Sequence[int], seed: int,
                         reserved: Sequence[int] = (),
                         timeout: Optional[float] = None) -> str:
    """add_fit_prime for nana configs: the new pool is PINNED (same pins,
    new prime), never generic -- growing a generic pool under a nana job
    would poison the CRT set with off-slice fits."""
    from .formats import big_uint_primes
    P = big_uint_primes()
    used = {read_ints(os.path.join(cfg.dbdir(pid), "red_common"))[0]
            for pid in range(len(cfg.datanames))}
    blocked = used | set(reserved)
    idx = next(i for i, p in enumerate(P) if p not in blocked)
    prime = P[idx]
    name = f"{cfg.job}_p{idx}"
    if os.path.exists(os.path.join(cfg.searchdir, "database", name)):
        raise BladeGateError(f"add_fit_prime_pinned: database {name} exists "
                             f"but is not in cfg.datanames -- refusing")
    n_pool = _file_rows(cfg.dbdir(0))
    init_pinned_database(cfg.searchdir, name, cfg.datanames[0], spec,
                         prime=prime)
    cfg.datanames.append(name)
    new_pid = len(cfg.datanames) - 1
    t0 = time.time()
    make_pinned_grower(cfg, kc, seed=seed, spec=spec, pids=(new_pid,))(n_pool)
    cfg.log(f"[add_fit_prime_pinned] {name} (prime {prime}): 0 -> {n_pool} "
            f"pinned points in {time.time()-t0:.1f}s")
    for w in works:
        shutil.rmtree(os.path.join(cfg.workdir(w), "fit", name),
                      ignore_errors=True)
        run_fitrel_works(cfg, name, w, w + 1, timeout=timeout)
    return name


def nana_search(cfg: SearchFamily, spec: NanaSpec, kc: KiraCoeffs,
                works: Sequence[int], labels: Sequence[str], out_json: str,
                seed: int, extra_meta: Optional[dict] = None,
                prime_adder: Optional[Callable[[], str]] = None,
                max_primes: int = 8
                ) -> Tuple[SearchFamily, List[WorkOutcome], Optional[dict]]:
    """Reduced-analytic blade_search.  Feature OFF (spec.ncons == 0) is the
    UNMODIFIED classic path (blade_search without bisection) on cfg itself.
    Feature ON: derive the nana config (make_nana_config), escalate, fit at
    every pinned prime, export WITH the reduced_analytic marker (fail-closed
    both ways).  Returns (nana_cfg, outcomes, doc_or_None)."""
    if spec.ncons == 0:
        outcomes, doc = blade_search(cfg, works, labels, out_json,
                                     bisect=False, extra_meta=extra_meta)
        return cfg, outcomes, doc
    ncfg = make_nana_config(cfg, spec, kc, seed, works)
    outcomes = [escalate_work(ncfg, w) for w in works]
    closed = [o.work for o in outcomes if o.closed]
    if not closed:
        return ncfg, outcomes, None
    for pid in range(len(ncfg.datanames)):
        for w in closed:
            shutil.rmtree(os.path.join(ncfg.workdir(w), "fit",
                                       ncfg.datanames[pid]), ignore_errors=True)
            run_fitrel_works(ncfg, ncfg.datanames[pid], w, w + 1)
    meta = dict(extra_meta or {})
    meta["reduced_analytic"] = spec.marker()
    if prime_adder is None:
        prime_adder = lambda: add_fit_prime_pinned(ncfg, kc, spec, closed, seed)
    doc = export_relations(ncfg, closed, labels, out_json, meta,
                           prime_adder=prime_adder, max_primes=max_primes)
    return ncfg, outcomes, doc
