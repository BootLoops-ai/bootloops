"""blade.semibl -- port of BLSolve/SolveSemiBL.wl (+ the
CollectRelation/ReconstructRelation choreography): turn a CLOSED
block-triangular form into a FULL SYMBOLIC reduction table
(target -> sum_j c_j(params) * master_j, c_j exact rational functions),
in a kira_table-comparable form.

Upstream (SemiBL is the production path; FullBL is abandoned upstream):

    BLSolveSemiBL[job, dataprefix, paras, nana, nthreads, key]:
      per prime pid = 0,1,...:
        BLDumpSamplePoints   -> points_<pid>.fflow   (R ordering)
        gatherSamplePoints   -> points_<pid>.bt      (slice-bucketed)
        blEvaluateAndDump    -> eval/<name>_<i>.fflow (BLSSolve per bucket,
                                chunked by MaxChunkMemory), eval/<name>.txt
        BLReconstructMod     -> recmod/<name>_{coeff,mono,part}
      dynamicrr coprocess (1 s state-file polling upstream)
      toRationalMPI + nzlearn scatter -> results/table

Port mapping (this module):
      SemiBLJob.stage_prime(pid)      = dumppoints + bucket + chunked ssolve
      SemiBLJob.run_recmod(pid)       = explicit per-prime recmod (gated)
      SemiBLJob.run_reconstruction()  = dynamicrr with GROW-ON-DEMAND priming
                                        (50 ms handshake loop, no 1 s polling;
                                        pipeline.run_dynamicrr_generic)
      SemiBLJob.assemble_table()      = toRationalMPI + element scatter,
                                        PENDING markers for open blocks
    BLCollectRelations / BLReconstructRelations (CollectRelation.wl,
    ReconstructRelation.wl) are already ported as search.export_relations
    (fit tables -> CRT+Wang -> exact BT relations); semibl CONSUMES that
    layer's closed works and does not duplicate it.

THE THREE PARAMETER ORDERINGS (upstream SolveSemiBL.wl:44-73; plan risk #1)

    I = BLIBPParameter        ordering of the numeric IBP database
                              (kin_table columns == red_ps coordinates)
    S = BLSearchParameter     ordering of the polynomial ansatz / fit tables
    R = blReconstructParameter ordering of degrees.fflow, points.fflow and
                              the reconstructed monomials (may be
                              complexity-sorted, blSortDegrees)

    In every tree this port produces, S == I by construction (the kinematics
    feed is written by blade.ansatz in search order); Orderings gates
    this instead of assuming it.  R may be any permutation of S.

    ptsord (upstream: an *unsigned* whose DIGITS are positions,
    iofflow.cpp:105-111 -- ordlist[i] = position of R[i] in S, 1-based;
    xin[ordlist[i]-1] = point[i]): REPLACED internally by an explicit
    permutation list (Orderings.ptsord_positions()).  The digit string is
    only formed at the C boundary and is REFUSED (typed error) for
    nparsin > 9 -- the C parser reads single digits, so >= 10 parameters
    cannot be expressed AT ALL for the shipped ssolve.  For such systems use
    evaluator= (the pure-python fallback lane, no ptsord involved).
    run_permutation_tests() is the executable spec of these semantics.

SLICE BUCKETING (upstream gatherSamplePoints, SolveSemiBL.wl:296-305):
    points sharing the values of the ncons = nparsin - nana conservative
    parameters form a slice; upstream re-fits the database per slice when a
    slice holds more points than fitnumber.  This port buckets and solves
    per slice (start/size windows into the gathered .bt file) but REFUSES
    ncons > 0 when the tree's fit tables were made at generic kinematics
    (per-slice database refit is search-lane machinery, not re-implemented
    here).  Both real families (db, dbox) run with nana == nparsin.

Gate discipline: blade.pipeline (rc never gates; artifact bytes,
dimensions, cross-prime identities, and dynamicrr state==2 do).  Failures
raise BladeGateError; nothing partial is ever written silently -- open
blocks appear in the final table ONLY as explicit PENDING markers.
"""

from __future__ import annotations

import json
import os
import time
from dataclasses import dataclass, field
from fractions import Fraction
from typing import Callable, Dict, List, Optional, Sequence, Tuple

from . import formats
from .formats import (DegreeInfo, DegreesFile, EvalFile, EvalList, PointsFile,
                      RecCoeff, RecMono, RecPart, RRRes, Scheme, SystemFile,
                      assemble_rational_functions, big_uint_primes, flags_size,
                      read_ints)
from .pipeline import (DEFAULT_BIN_DIR, BladeGateError, StageResult, _run,
                       run_dynamicrr_generic, run_recmod_generic)
from .ratrec import _evals_from_callback, scan_degrees

__all__ = ["Orderings", "SemiBLJob", "SemiBLTable", "sort_degrees_by_complexity",
           "gather_sample_points", "run_permutation_tests", "BladeGateError",
           "refuse_pinned_pool"]


def refuse_pinned_pool(dbdir: str) -> None:
    """FAIL-CLOSED nana contract (search.py reduced-analytic mode): a database
    dir carrying nana_pins.json was sampled on a PINNED slice -- fit tables /
    evals built on it hold only at the pins.  Full-symbolic reconstruction
    must refuse; per-point fitting in the pinned parameters is the consumer's
    job, one instance per distinct pin tuple (upstream SolveSemiBL.wl:296-304,
    404-410; port: search.plan_pin_instances / search.fit_at_pin).  A missing
    directory is NOT an error here (SemiBLJob bridges missing pools via the
    evaluator lane); a pinned one is."""
    marker = os.path.join(dbdir, "nana_pins.json")
    if os.path.exists(marker):
        pins = json.load(open(marker)).get("pins")
        raise BladeGateError(
            f"{dbdir}: reduced-analytic PINNED pool (pins={pins}) -- "
            f"full-symbolic SemiBL reconstruction is invalid on this data; "
            f"per-point fitting in the pinned parameters is the consumer's "
            f"job (search.plan_pin_instances / search.fit_at_pin)")


# ==========================================================================
# parameter orderings I / S / R
# ==========================================================================

@dataclass
class Orderings:
    """The three parameter orderings.  All are lists of parameter NAMES over
    the same set; validation is structural (same multiset, no duplicates)."""
    ibp: List[str]        # I
    search: List[str]     # S
    recon: List[str]      # R

    def __post_init__(self):
        for name, lst in (("I", self.ibp), ("S", self.search), ("R", self.recon)):
            if len(set(lst)) != len(lst):
                raise BladeGateError(f"orderings: duplicate parameter in {name}: {lst}")
        if set(self.ibp) != set(self.search) or set(self.search) != set(self.recon):
            raise BladeGateError(
                f"orderings: I/S/R are not permutations of one another: "
                f"I={self.ibp} S={self.search} R={self.recon}")

    @property
    def nparsin(self) -> int:
        return len(self.search)

    # ---- ptsord: R -> S (ssolve consumes S-ordered coords, points are R) ---
    def ptsord_positions(self) -> List[int]:
        """ordlist (1-based): ordlist[i] = position of R[i] in S.  C semantics
        (iofflow.cpp:130): xin[ordlist[i]-1] = point[i], i.e. the i-th
        R-ordered coordinate lands at its S position.  Upstream construction:
        FromDigits[FirstPosition[BLSearchParameter,#]&/@blReconstructParameter]
        (SolveSemiBL.wl:370)."""
        return [self.search.index(p) + 1 for p in self.recon]

    def ptsord_digits(self) -> str:
        """C-boundary digit encoding.  REFUSED (typed) beyond 9 parameters:
        the C parser peels single decimal digits (iofflow.cpp:108-111), so a
        position >= 10 is unrepresentable for the shipped binaries."""
        pos = self.ptsord_positions()
        if self.nparsin > 9 or any(p > 9 for p in pos):
            raise BladeGateError(
                f"ptsord digit encoding impossible for {self.nparsin} parameters "
                f"(C reads single digits, iofflow.cpp:105-111); use the "
                f"evaluator= fallback lane for >= 10-parameter systems")
        return "".join(str(p) for p in pos)

    # ---- coordinate permutations ------------------------------------------
    def r_to_s(self, coords: Sequence) -> List:
        """Reorder an R-ordered coordinate list into S ordering (reference
        implementation of what ssolve does internally with ordlist)."""
        out = [None] * self.nparsin
        for i, pos in enumerate(self.ptsord_positions()):
            out[pos - 1] = coords[i]
        return out

    def s_to_r(self, coords: Sequence) -> List:
        return [coords[pos - 1] for pos in self.ptsord_positions()]

    def s_to_i(self, coords: Sequence) -> List:
        """S-ordered coords -> I ordering (upstream ffEvaluateAndDump:
        posi = FirstPosition[blReconstructParameter,#]&/@BLIBPParameter,
        applied to R-ordered points; here factored through S)."""
        return [coords[self.search.index(p)] for p in self.ibp]

    def r_to_i(self, coords: Sequence) -> List:
        return self.s_to_i(self.r_to_s(coords))

    @classmethod
    def identity(cls, params: Sequence[str]) -> "Orderings":
        return cls(list(params), list(params), list(params))


def run_permutation_tests(seed: int = 20260709, rounds: int = 200) -> int:
    """Executable spec for the ordering semantics (plan risk #1).  Checks,
    for random permutations at every nparsin 1..9 (and refusal at 10+):
      1. r_to_s composed with s_to_r is the identity (both ways);
      2. r_to_s reproduces the literal C loop (xin[ordlist[i]-1]=pt[i]) with
         ordlist decoded from the DIGIT string exactly as iofflow.cpp does;
      3. semantic anchor: evaluating a random monomial with S-bound values
         equals evaluating it after shuffling into R and mapping back;
      4. s_to_i is the inverse permutation of i_to_s implied by names.
    Returns the number of checks performed; raises BladeGateError on any
    mismatch."""
    import random
    rng = random.Random(seed)
    checks = 0
    for _ in range(rounds):
        n = rng.randrange(1, 10)
        params = [f"z{k}" for k in range(n)]
        S = params[:]
        R = params[:]
        I = params[:]
        rng.shuffle(R)
        rng.shuffle(I)
        o = Orderings(I, S, R)
        vals = {p: rng.randrange(1, 1 << 20) for p in params}
        r_coords = [vals[p] for p in R]
        s_coords = [vals[p] for p in S]
        i_coords = [vals[p] for p in I]
        # 1. round trips
        if o.r_to_s(r_coords) != s_coords or o.s_to_r(s_coords) != r_coords:
            raise BladeGateError(f"permutation test: R<->S round trip failed "
                                 f"(S={S}, R={R})")
        # 2. literal C decode of the digit string
        digits = o.ptsord_digits()
        ptsord = int(digits)
        ordlist = [0] * n
        cp = ptsord
        for i in range(n, 0, -1):          # iofflow.cpp:108-111
            ordlist[i - 1] = cp % 10
            cp //= 10
        xin = [None] * n
        for i in range(n):                 # iofflow.cpp:130
            xin[ordlist[i] - 1] = r_coords[i]
        if xin != s_coords:
            raise BladeGateError(f"permutation test: C digit decode differs "
                                 f"from r_to_s (S={S}, R={R}, ptsord={digits})")
        # 3. semantic anchor: name-bound evaluation is ordering-free
        exps = [rng.randrange(0, 3) for _ in range(n)]
        p = (1 << 61) - 1
        def ev(coords, order):
            acc = 1
            for x, pn in zip(coords, order):
                acc = acc * pow(x, exps[params.index(pn)], p) % p
            return acc
        if ev(s_coords, S) != ev(o.r_to_s(r_coords), S) != ev(r_coords, R):
            raise BladeGateError("permutation test: semantic anchor failed")
        # 4. I lane
        if o.s_to_i(s_coords) != i_coords or o.r_to_i(r_coords) != i_coords:
            raise BladeGateError(f"permutation test: S/R->I failed (I={I})")
        checks += 4
    # refusal beyond 9 parameters
    big = [f"z{k}" for k in range(10)]
    try:
        Orderings.identity(big).ptsord_digits()
    except BladeGateError:
        checks += 1
    else:
        raise BladeGateError("permutation test: 10-parameter digit encoding "
                             "was NOT refused")
    return checks


# ==========================================================================
# blSortDegrees port (complexity-sorted R ordering)
# ==========================================================================

def sort_degrees_by_complexity(params: Sequence[str], deg: DegreesFile
                               ) -> Tuple[List[str], DegreesFile]:
    """Port of blSortDegrees (SolveSemiBL.wl:132-153): sort parameters by
    descending per-variable complexity max_out(num_max - num_min + den_max -
    den_min) and permute every output's per-var degree quads accordingly.
    Returns (new R parameter list, reordered DegreesFile).  Ties keep the
    original relative order (Mathematica Ordering is stable)."""
    if len(params) != deg.nparsin:
        raise BladeGateError("sort_degrees: params length != degrees nparsin")
    compl = []
    for i in range(deg.nparsin):
        m = 0
        for d in deg.info:
            nmax, nmin, dmax, dmin = d.var[i]
            m = max(m, nmax - nmin + dmax - dmin)
        compl.append(m)
    order = sorted(range(deg.nparsin), key=lambda i: (-compl[i], i))
    info = [DegreeInfo(d.numdeg, d.dendeg, [d.var[i] for i in order])
            for d in deg.info]
    return [params[i] for i in order], DegreesFile(deg.nparsin, deg.nparsout, info)


# ==========================================================================
# slice bucketing (gatherSamplePoints port)
# ==========================================================================

def gather_sample_points(pts: PointsFile, nparsin: int, prime: int,
                         const_positions: Sequence[int]
                         ) -> Tuple[PointsFile, List[Tuple[int, int]]]:
    """Port of gatherSamplePoints (SolveSemiBL.wl:296-305): keep rows of this
    prime (gate: EVERY row must carry it -- one dump is one prime), group by
    the coordinates at const_positions (0-based positions in the R-ordered
    row; upstream uses the leading ncons positions), first-seen order
    (GatherBy semantics).  Returns (gathered PointsFile, [(start, size)] per
    bucket, windows into the gathered file for ssolve start/size)."""
    for i, row in enumerate(pts.rows):
        if row[nparsin] != prime:
            raise BladeGateError(f"gather: row {i} prime {row[nparsin]} != {prime}")
    buckets: Dict[Tuple[int, ...], List[int]] = {}
    order: List[Tuple[int, ...]] = []
    for i, row in enumerate(pts.rows):
        key = tuple(row[p] for p in const_positions)
        if key not in buckets:
            buckets[key] = []
            order.append(key)
        buckets[key].append(i)
    gathered = PointsFile(pts.n, [])
    windows: List[Tuple[int, int]] = []
    at = 0
    for key in order:
        idx = buckets[key]
        gathered.rows += [pts.rows[i] for i in idx]
        windows.append((at, len(idx)))
        at += len(idx)
    return gathered, windows


def chunk_windows(windows: Sequence[Tuple[int, int]], row_words: int,
                  max_chunk_bytes: int) -> List[Tuple[int, int]]:
    """Split bucket windows into ssolve (start,size) chunks so one chunk's
    eval stays under max_chunk_bytes (upstream MaxChunkMemory -> 10^9 default,
    blEvaluateAndDump).  Chunks never span buckets (upstream solves per
    slice)."""
    cap = max(1, max_chunk_bytes // (8 * row_words))
    out: List[Tuple[int, int]] = []
    for start, size in windows:
        at = start
        while size > 0:
            take = min(size, cap)
            out.append((at, take))
            at += take
            size -= take
    return out


# ==========================================================================
# the job
# ==========================================================================

@dataclass
class SemiBLJob:
    """One SemiBL symbolic-reconstruction job over a closed BT search tree.

    searchdir/            database/<datanames[pid]>/red_*   (fit pools)
                          kinematics/<job>/kin_*            (S-ordered)
                          <job>/<w>/{sch_*, config/, fit/<dataname>/}
    key_dir/              degrees.fflow (R-ordered), points/, eval/, recmod/

    targets: global ids the table must cover; the CLOSED subset (works with
    fit tables) is solved, the rest become explicit PENDING rows.
    elements: ssolve output layout [(target_gid, master_gid)] -- target-major
    in `system_targets` order, masters ascending (the learned support,
    upstream nzlearn); gated against the exact evaluator when provided.
    """
    searchdir: str
    job: str
    key_dir: str
    datanames: List[str]                 # index == fflow prime id
    orderings: Orderings
    system_targets: List[int]            # tarid, order fixes output layout
    elements: List[Tuple[int, int]]      # (target_gid, master_gid)
    degrees: DegreesFile                 # R-ordered, EXACT totals
    labels: Optional[List[str]] = None   # gid -> printable name
    pending_targets: List[int] = field(default_factory=list)
    pending_reason: str = ""
    nana: Optional[int] = None           # default: all parameters analytic
    nthreads: int = 8
    bin_dir: str = DEFAULT_BIN_DIR
    dumppoints_bin: Optional[str] = None
    maxdeg_cli: int = 1000
    max_chunk_bytes: int = 10 ** 9       # upstream MaxChunkMemory default
    stage_timeout: float = 9000.0
    dynamicrr_timeout: float = 3600.0
    evaluator: Optional[Callable] = None # fallback lane: (S-ordered coords, prime)->vals
    force_evaluator: bool = False        # True: never invoke ssolve (>9-param lane)
    log: Callable[[str], None] = print

    # ---- derived ------------------------------------------------------------
    def __post_init__(self):
        self.searchdir = os.path.abspath(self.searchdir)
        self.key_dir = os.path.abspath(self.key_dir)
        n = self.orderings.nparsin
        if self.nana is None:
            self.nana = n
        if not (0 <= self.nana <= n):
            raise BladeGateError(f"nana={self.nana} outside 0..{n}")
        if self.nana != n and self.evaluator is None:
            raise BladeGateError(
                "ncons > 0 with the ssolve lane needs per-slice database "
                "refits (upstream BLDatabase per slice) -- not ported; run "
                "with nana == nparsin or supply evaluator=")
        # FAIL-CLOSED nana contract: fit tables made on a PINNED pool
        # (search.init_pinned_database -> nana_pins.json) are valid only at
        # the pins; a full-symbolic SemiBL table must never consume them.
        for dn in self.datanames:
            refuse_pinned_pool(os.path.join(self.searchdir, "database", dn))
        if self.degrees.nparsin != n:
            raise BladeGateError("degrees nparsin != orderings nparsin")
        if self.degrees.nparsout != len(self.elements):
            raise BladeGateError(
                f"degrees nparsout {self.degrees.nparsout} != "
                f"#elements {len(self.elements)}")
        tset = {t for t, _ in self.elements}
        if tset != set(self.system_targets):
            raise BladeGateError(
                f"elements cover targets {sorted(tset)} != system_targets "
                f"{sorted(set(self.system_targets))}")
        # element order gate: target-major in system_targets order, masters
        # ascending (the layout cache_eval emits, block.c:1067-1082)
        want = [(t, m) for t in self.system_targets
                for m in sorted(mm for tt, mm in self.elements if tt == t)]
        if list(self.elements) != want:
            raise BladeGateError(
                "elements not target-major/master-ascending in system_targets "
                "order -- would misalign with ssolve output (block.c cache_eval)")
        os.makedirs(self.key_dir, exist_ok=True)
        for sub in ("points", "eval", "recmod", "system"):
            os.makedirs(os.path.join(self.key_dir, sub), exist_ok=True)
        self.degrees_path = os.path.join(self.key_dir, "degrees.fflow")
        self.degrees.write(self.degrees_path)
        self.walls: Dict[str, float] = {}
        self.meta: Dict[str, object] = {"npoints": {}, "chunks": {}}

    def dbdir(self, pid: int) -> str:
        return os.path.join(self.searchdir, "database", self.datanames[pid])

    def jobdir(self) -> str:
        return os.path.join(self.searchdir, self.job)

    def label(self, gid: int) -> str:
        return self.labels[gid] if self.labels else f"I{gid}"

    @property
    def nparsout(self) -> int:
        return len(self.elements)

    def pid_uses_ssolve(self, pid: int) -> bool:
        """The C ssolve lane needs a database pool + fit tables at this prime
        id.  Pools may legitimately SKIP ids (e.g. the dbox bank reserves
        BIG_UINT[3]/[4] as never-fitted oracle primes, so data_p3/data_p4 do
        not exist) while recmod/dynamicrr still index primes consecutively --
        such gaps are BRIDGED by the evaluator callback (identical contract
        to the >9-parameter lane).  The bridge is recorded in
        meta['fallback_pids'] so the caller can disclose which primes were
        consumed outside the fit-table lane."""
        if self.force_evaluator:
            return False
        return pid < len(self.datanames) and os.path.isdir(self.dbdir(pid))

    # ---- inputs for the C stages ---------------------------------------------
    def _write_system(self, pid: int) -> str:
        common = read_ints(os.path.join(self.dbdir(pid), "red_common"))
        prime, npara, nint, nmaster = common[0], common[1], common[2], common[3]
        if prime != big_uint_primes()[pid]:
            raise BladeGateError(
                f"{self.dbdir(pid)}/red_common: prime {prime} != "
                f"BIG_UINT_PRIMES[{pid}] (recmod/dynamicrr index primes by id)")
        if npara != self.orderings.nparsin:
            raise BladeGateError(f"red_common npara {npara} != nparsin")
        sysf = SystemFile(prime, npara, nint, nmaster, list(self.system_targets),
                          [os.path.abspath(self.jobdir())], self.datanames[pid])
        path = os.path.join(self.key_dir, "system", f"system_p{pid}.txt")
        sysf.write(path)
        return path

    # ---- per-prime staging: points -> buckets -> chunked evals ---------------
    def stage_prime(self, pid: int) -> List[str]:
        """BLDumpSamplePoints + gatherSamplePoints + blEvaluateAndDump for one
        prime.  Returns the chunk eval paths (also recorded in
        eval/eval_p<pid>.txt).  Idempotent: existing complete chunks are kept
        (upstream break-point semantics)."""
        t0 = time.time()
        n = self.orderings.nparsin
        prime = big_uint_primes()[pid]
        fs = flags_size(self.nparsout)
        pts_path = os.path.join(self.key_dir, "points", f"points_p{pid}.fflow")
        if not os.path.exists(pts_path):
            from .ratrec import dump_points
            dump_points(n, self.nparsout, pid, self.maxdeg_cli,
                        self.degrees_path, pts_path, self.bin_dir,
                        self.dumppoints_bin)
        pts = PointsFile.read(pts_path)
        if pts.n != n + fs:
            raise BladeGateError(f"{pts_path}: header n={pts.n} != {n}+{fs}")
        # slice bucketing: conservative params = S[:ncons] at their R positions
        ncons = n - self.nana
        const_pos = [self.orderings.recon.index(p)
                     for p in self.orderings.search[:ncons]]
        gathered, windows = gather_sample_points(pts, n, prime, const_pos)
        bt_path = os.path.join(self.key_dir, "points", f"points_p{pid}.bt")
        gathered.write(bt_path)
        row_words = n + 1 + fs + self.nparsout
        chunks = chunk_windows(windows, row_words, self.max_chunk_bytes)
        self.meta["npoints"][pid] = len(gathered.rows)
        self.meta["chunks"][pid] = len(chunks)

        evalname = f"eval_p{pid}"
        paths: List[str] = []
        for ci, (start, size) in enumerate(chunks, 1):
            epath = os.path.join(self.key_dir, "eval", f"{evalname}_{ci}.fflow")
            paths.append(epath)
            want = EvalFile.expected_bytes(n, self.nparsout, size)
            if os.path.exists(epath) and os.path.getsize(epath) == want:
                continue                       # upstream break-point semantics
            if self.pid_uses_ssolve(pid):
                self._eval_chunk_ssolve(pid, bt_path, start, size, epath)
            else:
                if self.evaluator is None:
                    raise BladeGateError(
                        f"prime id {pid}: no database pool "
                        f"({self.dbdir(pid) if pid < len(self.datanames) else '<beyond datanames>'})"
                        f" and no evaluator= fallback supplied")
                fb = self.meta.setdefault("fallback_pids", [])
                if pid not in fb:
                    fb.append(pid)
                self._eval_chunk_callback(gathered, start, size, prime, epath)
            got = os.path.getsize(epath)
            if got != want:
                raise BladeGateError(
                    f"{epath}: {got} bytes != expected {want} (chunk {ci}, "
                    f"{size} rows x {self.nparsout} outputs) -- ssolve's "
                    f"LEARNED output support differs from the declared "
                    f"elements (system_learn samples a random point per run; "
                    f"re-check the element list / rerun)")
            ev = EvalFile.read(epath)
            if ev.nparsin != n or ev.nout_plus_flags != self.nparsout + fs \
                    or len(ev.rows) != size:
                raise BladeGateError(f"{epath}: header inconsistent")
        EvalList([os.path.abspath(p) for p in paths]).write(
            os.path.join(self.key_dir, "eval", f"{evalname}.txt"))
        self.walls[f"stage_p{pid}"] = time.time() - t0
        return paths

    def _eval_chunk_ssolve(self, pid: int, bt_path: str, start: int, size: int,
                           epath: str) -> None:
        """BLSSolve port: ssolve <nparsin> <pid> <start> <size> <nthreads>
        <system> <ptsord> <points.bt> <eval> (ssolve.cpp:37; start/size window
        into the gathered points file, SolveSemiBL.wl:416)."""
        spath = self._write_system(pid)
        ptsord = self.orderings.ptsord_digits()   # typed refusal > 9 params
        res = _run([os.path.join(self.bin_dir, "ssolve"),
                    str(self.orderings.nparsin), str(pid), str(start),
                    str(size), str(self.nthreads), spath, ptsord,
                    bt_path, epath],
                   cwd=self.jobdir(), timeout=self.stage_timeout,
                   stage=f"ssolve[p{pid}@{start}+{size}]")
        if not os.path.exists(epath):
            raise BladeGateError(f"ssolve: {epath} missing (rc={res.rc}, rc "
                                 f"lies anyway). log tail: {res.log[-400:]}")

    def _eval_chunk_callback(self, gathered: PointsFile, start: int, size: int,
                             prime: int, epath: str) -> None:
        """Fallback lane (>9 params, or no fit tables at this prime): the
        caller's evaluator( S-ordered coords, prime ) -> values in ELEMENT
        order.  Points are R-ordered on disk; the eval rows keep the R-ordered
        coords (recmod matches them against degrees.fflow) while the
        evaluator sees S ordering -- the same contract ssolve implements with
        ordlist."""
        window = PointsFile(gathered.n, gathered.rows[start:start + size])
        o = self.orderings

        def ev_r(coords_r, p):
            return self.evaluator(tuple(o.r_to_s(list(coords_r))), p)

        ev = _evals_from_callback(window, o.nparsin, self.nparsout, prime, ev_r)
        ev.write(epath)

    # ---- explicit per-prime recmod -------------------------------------------
    def run_recmod(self, pid: int) -> StageResult:
        t0 = time.time()
        res = run_recmod_generic(
            os.path.join(self.bin_dir, "recmod"), self.orderings.nparsin,
            self.nparsout, pid, self.nthreads, self.degrees_path,
            os.path.join(self.key_dir, "eval", f"eval_p{pid}.txt"),
            os.path.join(self.key_dir, "recmod"), "rec", self.maxdeg_cli,
            self.stage_timeout)
        self.walls[f"recmod_p{pid}"] = time.time() - t0
        return res

    def _staged_primes(self) -> int:
        k = 0
        while os.path.exists(os.path.join(self.key_dir, "recmod",
                                          f"rec_{k}_coeff")):
            k += 1
        return k

    # ---- dynamicrr with grow-on-demand priming --------------------------------
    def run_reconstruction(self, min_primes: int = 2,
                           max_primes: Optional[int] = None) -> RRRes:
        """Prime the recmod outputs for min_primes, run dynamicrr (state-file
        HANDSHAKE at 50 ms via pipeline.run_dynamicrr_generic -- the upstream
        1 s Pause loop, SolveSemiBL.wl:562-584, is replaced); whenever
        dynamicrr wants a prime beyond those staged, stage exactly one more
        (dumppoints+ssolve+recmod) and rerun.  Typed failure when max_primes
        is exhausted; nothing is emitted on failure."""
        t0 = time.time()
        if max_primes is None:
            max_primes = len(self.datanames)
        if max_primes < 2:
            raise BladeGateError("reconstruction needs >= 2 primes")
        k = self._staged_primes()
        while k < min_primes:
            self.stage_prime(k)
            self.run_recmod(k)
            k += 1
        while True:
            try:
                coeff0 = RecCoeff.read(os.path.join(self.key_dir, "recmod",
                                                    "rec_0_coeff"))
                res = run_dynamicrr_generic(
                    os.path.join(self.bin_dir, "dynamicrr"),
                    len(coeff0.tokens), self.nthreads,
                    os.path.join(self.key_dir, "recmod"), "rec", "rrres",
                    nprimes_available=k, timeout=self.dynamicrr_timeout)
                break
            except BladeGateError as e:
                if "requested prime" not in str(e):
                    raise
                if k >= max_primes:
                    raise BladeGateError(
                        f"reconstruction needs > {max_primes} primes "
                        f"(coefficient heights exceed the staged budget); "
                        f"NOT emitting a table") from e
                self.log(f"[semibl] dynamicrr wants prime {k}; staging it")
                self.stage_prime(k)
                self.run_recmod(k)
                k += 1
        self.meta["primes_consumed"] = k
        self.meta["dynamicrr_gates"] = res.gates
        self.walls["reconstruction_total"] = time.time() - t0
        rr = RRRes.read(os.path.join(self.key_dir, "recmod", "rrres"))
        return rr

    # ---- final table -----------------------------------------------------------
    def assemble_table(self, symbols=None) -> "SemiBLTable":
        """toRationalMPI + element scatter (SolveSemiBL.wl:601-611): rrres +
        rec_0_{mono,part} -> rational functions in the R-ordered parameters,
        scattered onto (target, master) elements; open targets appear as
        explicit PENDING rows and NEVER as numbers."""
        import sympy as sp
        t0 = time.time()
        part = RecPart.read(os.path.join(self.key_dir, "recmod", "rec_0_part"))
        mono = RecMono.read(os.path.join(self.key_dir, "recmod", "rec_0_mono"),
                            self.orderings.nparsin)
        rr = RRRes.read(os.path.join(self.key_dir, "recmod", "rrres"))
        if symbols is None:
            symbols = sp.symbols(self.orderings.recon)
            if self.orderings.nparsin == 1:
                symbols = (symbols,)
        funcs = assemble_rational_functions(part, mono, rr.values, symbols)
        if len(funcs) != self.nparsout:
            raise BladeGateError("assembly: function count != elements")
        rows: Dict[int, Dict[int, object]] = {t: {} for t in self.system_targets}
        for (tg, mg), f in zip(self.elements, funcs):
            rows[tg][mg] = f
        self.walls["assemble"] = time.time() - t0
        return SemiBLTable(
            params_recon=list(self.orderings.recon),
            rows=rows,
            pending={t: self.pending_reason or "block open"
                     for t in self.pending_targets},
            labels=self.labels,
            meta={"walls": dict(self.walls), **self.meta})


# ==========================================================================
# the table
# ==========================================================================

@dataclass
class SemiBLTable:
    """kira_table-comparable symbolic reduction table.  rows[target_gid] =
    {master_gid: sympy expr in params_recon}; pending[target_gid] = reason
    (fail-closed: pending targets have NO numeric/symbolic entries)."""
    params_recon: List[str]
    rows: Dict[int, Dict[int, object]]
    pending: Dict[int, str]
    labels: Optional[List[str]] = None
    meta: Dict[str, object] = field(default_factory=dict)

    def label(self, gid: int) -> str:
        return self.labels[gid] if self.labels else f"I{gid}"

    def write_text(self, path: str) -> None:
        """Same register as phase0/reconstructed_table.txt, plus explicit
        PENDING rows."""
        import sympy as sp
        with open(path, "w") as f:
            for t in sorted(self.rows):
                for m in sorted(self.rows[t]):
                    f.write(f"target {self.label(t)} -> master "
                            f"{self.label(m)} :\n  {sp.sstr(self.rows[t][m])}\n")
            for t in sorted(self.pending):
                f.write(f"target {self.label(t)} : PENDING -- "
                        f"{self.pending[t]} (no coefficients emitted)\n")

    def write_json(self, path: str) -> None:
        import sympy as sp
        doc = {
            "params": self.params_recon,
            "rows": {str(t): {str(m): sp.sstr(e)
                              for m, e in sorted(self.rows[t].items())}
                     for t in sorted(self.rows)},
            "pending": {str(t): self.pending[t] for t in sorted(self.pending)},
            "labels": self.labels,
            "meta": self.meta,
        }
        with open(path, "w") as f:
            json.dump(doc, f, indent=1)


# ==========================================================================
# element / support learning (nzlearn analogue)
# ==========================================================================

def elements_from_database(searchdir: str, dataname: str,
                           system_targets: Sequence[int]
                           ) -> List[Tuple[int, int]]:
    """The (target, master) support from the numeric database's stored
    pattern (red_part/red_posi -- Kira-provenance nonzero pattern), laid out
    target-major in system_targets order, masters ascending: the layout
    ssolve's learned cache emits (block.c system_sol_cache_info/cache_eval).
    Master gid = nint - nmaster + m."""
    dbdir = os.path.join(searchdir, "database", dataname)
    common = read_ints(os.path.join(dbdir, "red_common"))
    nint, nmaster = common[2], common[3]
    part = read_ints(os.path.join(dbdir, "red_part"))
    posi = read_ints(os.path.join(dbdir, "red_posi"))
    nt = nint - nmaster
    support: Dict[int, List[int]] = {t: [] for t in system_targets}
    for m in range(nmaster):
        for gi in posi[part[m]:part[m + 1]]:
            if gi < nt and gi in support:
                support[gi].append(nt + m)
    out = []
    for t in system_targets:
        if not support[t]:
            raise BladeGateError(f"database pattern: target {t} has empty "
                                 f"master support")
        out += [(t, m) for m in sorted(support[t])]
    return out


def gate_elements_against_evaluator(elements: Sequence[Tuple[int, int]],
                                    reduce_point: Callable, params_s: int,
                                    prime: int, npoints: int = 2,
                                    seed: int = 20260709) -> None:
    """Support gate: at `npoints` random S-ordered points, the evaluator's
    nonzero support must EQUAL the declared elements' support union (both
    directions -- a missing element misaligns every later output; an extra
    one triggers the byte gate only after an expensive run).  reduce_point
    (S-coords, prime) -> {(target_gid, master_gid): value}."""
    import random
    rng = random.Random(seed)
    seen = set()
    per_point = []
    for _ in range(npoints):
        coords = tuple(rng.randrange(1, prime) for _ in range(params_s))
        vals = reduce_point(coords, prime)
        nz = {k for k, v in vals.items() if v}
        per_point.append(nz)
        seen |= nz
    declared = set(elements)
    if seen - declared:
        raise BladeGateError(f"support gate: evaluator nonzero at UNDECLARED "
                             f"elements {sorted(seen - declared)[:8]}")
    if declared - seen:
        raise BladeGateError(
            f"support gate: declared elements never nonzero at {npoints} "
            f"random points: {sorted(declared - seen)[:8]} (accidental-zero "
            f"probability ~ deg/p per point -- if genuinely nonzero, "
            f"resample; otherwise drop them)")
