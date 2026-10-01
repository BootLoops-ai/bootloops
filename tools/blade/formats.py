"""blade.formats — readers/writers for every Blade C-pipeline contract file.

Format library (BootLoops Blade port).  Layouts are taken from the
authoritative upstream C readers and cross-checked against sample trees;
where documentation and C disagreed, C won (citations inline).

Authoritative C sources:
  blade/src/search/file.c         file_read_int/file_write_int ("%d", NO newline)
  blade/src/search/list.c         list_write_int ("%d\\n" each, :46-53),
                                  table_write_int ("%d " each + "\\n" per row, :190-199),
                                  table_write_llu ("%llu " + "\\n", :411-420);
                                  all readers are fscanf("%d"/"%llu") --
                                  whitespace-insensitive, NO error on short files.
  blade/src/search/database.c:14  database_init: red_common/part/posi/ps/data
  blade/src/search/kinematics.c:14 kin_common/kin_table
  blade/src/search/template.c     sch_*/in_*/out_* (:15,:105,:229-278)
  blade/src/search/searchalg.c    tmp_* (write_template_g1 :474-510),
                                  fit tables (fit_relations :635-686)
  blade/src/ssolve/block.c        template_info_init (:7, tmp_*),
                                  block_init (:80, fit tables, multi_sch_intid),
                                  system_init_from_txt_parallel (:778, system txt;
                                  kinematics dir derived :849 as
                                  blockdir[0]/../kinematics/<basename(blockdir[0])>)
  blade/src/ssolve/iofflow.cpp    system_dump_evaluations (:205, eval file)
  blade/src/fflow_interface/fflowC.cc  put_sparse_poly_mono/coeff/part and
                                  put_mprat (:67-165): every token "%s "/"%u "
                                  followed by ONE space, no newlines.
  finiteflow/src/alg_reconstruction.cc load_samples (:331, points file)
  finiteflow/src/alg_mp_reconstruction.cc algorithm_dump_degree_info (:99-131)
  finiteflow/include/fflow/primes.hh   BIG_UINT_PRIMES = the 201 largest
                                  63-bit primes (Select[2^63 - Range[9425], PrimeQ])

Writer styles (byte-exact per producer):
  FAB   : " ".join(row) + "\\n" per row      (phase0 fabricator inputs)
  CSCAL : "%d" single int, NO trailing newline (file_write_int, dynamicrr state)
  CLIST : "%d\\n" per entry                   (list_write_int: tmp_indep)
  CTAB  : "%d "/"%llu " per entry + "\\n" per row (table_write_int/llu:
          tmp_var, out_var, out_sol, out_rel, fit tables)
  TOK   : token + " " each, nothing else     (rec_*_coeff/_mono/_part, rrres)
  binary: flat little-endian uint64          (points/eval/degrees .fflow)

All plain python3 + fractions + json; sympy imported lazily only for assembly.
"""

from __future__ import annotations

import json
import os
import struct
from dataclasses import dataclass, field
from fractions import Fraction
from typing import Callable, List, Optional, Sequence, Tuple


class BladeFormatError(Exception):
    """A contract file does not parse / has inconsistent dimensions."""


# --------------------------------------------------------------------------
# fflow primes: the 201 largest 63-bit primes (finiteflow primes.hh:18).
# Generated with a deterministic Miller-Rabin (bases below are proven
# deterministic for n < 3.3e24) instead of hard-coding 201 literals.
# --------------------------------------------------------------------------

_MR_BASES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)


def _is_prime_u64(n: int) -> bool:
    if n < 2:
        return False
    for p in _MR_BASES:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in _MR_BASES:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


_BIG_UINT_PRIMES: Optional[List[int]] = None


def big_uint_primes() -> List[int]:
    """The fflow BIG_UINT_PRIMES list (descending). len == 201, validated."""
    global _BIG_UINT_PRIMES
    if _BIG_UINT_PRIMES is None:
        primes = [2**63 - k for k in range(1, 9426) if _is_prime_u64(2**63 - k)]
        if len(primes) != 201 or primes[0] != 9223372036854775783 \
                or primes[1] != 9223372036854775643 or primes[2] != 9223372036854775549:
            raise BladeFormatError("BIG_UINT_PRIMES generation failed self-check")
        _BIG_UINT_PRIMES = primes
    return _BIG_UINT_PRIMES


def prime_id(prime: int) -> int:
    """Index of `prime` in BIG_UINT_PRIMES; raises if absent (contract:
    database/recmod/dynamicrr primes MUST be fflow primes)."""
    try:
        return big_uint_primes().index(prime)
    except ValueError:
        raise BladeFormatError(f"prime {prime} is not an fflow BIG_UINT_PRIMES entry")


# --------------------------------------------------------------------------
# low-level text/binary helpers
# --------------------------------------------------------------------------

def read_tokens(path: str) -> List[str]:
    with open(path, "r") as f:
        return f.read().split()


def read_ints(path: str) -> List[int]:
    return [int(t) for t in read_tokens(path)]


def _rows(flat: Sequence, ncol: int, what: str) -> List[list]:
    if ncol <= 0:
        if flat:
            raise BladeFormatError(f"{what}: ncol={ncol} but {len(flat)} tokens")
        return []
    if len(flat) % ncol:
        raise BladeFormatError(f"{what}: {len(flat)} tokens not divisible by ncol={ncol}")
    return [list(flat[i:i + ncol]) for i in range(0, len(flat), ncol)]


def write_fab_rows(path: str, rows: Sequence[Sequence]) -> None:
    """Fabricator style: ' '.join per row + newline (phase0 wtext)."""
    with open(path, "w") as f:
        for r in rows:
            f.write(" ".join(str(x) for x in r))
            f.write("\n")


def write_c_scalar(path: str, v: int) -> None:
    """C file_write_int (file.c:37): '%d', no newline."""
    with open(path, "w") as f:
        f.write(str(int(v)))


def write_c_list(path: str, vals: Sequence[int]) -> None:
    """C list_write_int (list.c:46): '%d\\n' per entry."""
    with open(path, "w") as f:
        for v in vals:
            f.write(f"{int(v)}\n")


def write_c_table(path: str, rows: Sequence[Sequence[int]]) -> None:
    """C table_write_int/llu (list.c:190/:411): '%d ' per entry + '\\n' per row."""
    with open(path, "w") as f:
        for r in rows:
            for v in r:
                f.write(f"{int(v)} ")
            f.write("\n")


def write_tokens(path: str, tokens: Sequence) -> None:
    """fflow put_* style (fflowC.cc:67-165): every token followed by one space."""
    with open(path, "w") as f:
        for t in tokens:
            f.write(f"{t} ")


def read_u64_words(path: str) -> List[int]:
    data = open(path, "rb").read()
    if len(data) % 8:
        raise BladeFormatError(f"{path}: size {len(data)} not a multiple of 8")
    return list(struct.unpack(f"<{len(data)//8}Q", data))


def write_u64_words(path: str, words: Sequence[int]) -> None:
    with open(path, "wb") as f:
        f.write(struct.pack(f"<{len(words)}Q", *words))


# --------------------------------------------------------------------------
# Database: database/<dataname>/red_*   (reader: search/database.c:14)
# --------------------------------------------------------------------------

@dataclass
class Database:
    """red_common: 'prime npara nint nmaster nentry nps' (database.c:20).
    red_part: nmaster+1 prefix sums; red_posi: nentry global IDs grouped by
    master; red_ps: nps x npara sample coords; red_data: nps x nentry values.
    red_pid/red_cons are WL-only (no C binary reads them) and carried as raw
    text.  Masters are the LAST nmaster ids of 0..nint-1.  Last point row is
    the reserved self-check point (searchalg.c:135,:343)."""
    prime: int
    npara: int
    nint: int
    nmaster: int
    nentry: int
    nps: int
    part: List[int] = field(default_factory=list)
    posi: List[int] = field(default_factory=list)
    ps: List[List[int]] = field(default_factory=list)
    data: List[List[int]] = field(default_factory=list)
    pid_raw: Optional[str] = None    # red_pid, verbatim
    cons_raw: Optional[str] = None   # red_cons, verbatim

    def validate(self) -> None:
        if len(self.part) != self.nmaster + 1:
            raise BladeFormatError(f"red_part: {len(self.part)} != nmaster+1={self.nmaster+1}")
        if self.part[0] != 0 or self.part[-1] != self.nentry:
            raise BladeFormatError("red_part: not a prefix-sum list ending at nentry")
        if any(b < a for a, b in zip(self.part, self.part[1:])):
            raise BladeFormatError("red_part: not non-decreasing")
        if len(self.posi) != self.nentry:
            raise BladeFormatError(f"red_posi: {len(self.posi)} != nentry={self.nentry}")
        if len(self.ps) != self.nps or (self.ps and len(self.ps[0]) != self.npara):
            raise BladeFormatError("red_ps: wrong shape")
        if len(self.data) != self.nps or (self.data and len(self.data[0]) != self.nentry):
            raise BladeFormatError("red_data: wrong shape")
        prime_id(self.prime)  # must be an fflow prime

    @classmethod
    def read(cls, dirpath: str) -> "Database":
        common = read_ints(os.path.join(dirpath, "red_common"))
        if len(common) != 6:
            raise BladeFormatError(f"{dirpath}/red_common: expected 6 fields, got {len(common)}")
        prime, npara, nint, nmaster, nentry, nps = common
        db = cls(prime, npara, nint, nmaster, nentry, nps)
        db.part = read_ints(os.path.join(dirpath, "red_part"))
        db.posi = read_ints(os.path.join(dirpath, "red_posi"))
        db.ps = _rows(read_ints(os.path.join(dirpath, "red_ps")), npara, "red_ps")
        db.data = _rows(read_ints(os.path.join(dirpath, "red_data")), nentry, "red_data")
        for name, attr in (("red_pid", "pid_raw"), ("red_cons", "cons_raw")):
            p = os.path.join(dirpath, name)
            if os.path.exists(p):
                setattr(db, attr, open(p).read())
        db.validate()
        return db

    def write(self, dirpath: str) -> None:
        self.validate()
        os.makedirs(dirpath, exist_ok=True)
        write_fab_rows(os.path.join(dirpath, "red_common"),
                       [[self.prime, self.npara, self.nint, self.nmaster, self.nentry, self.nps]])
        write_fab_rows(os.path.join(dirpath, "red_part"), [self.part])
        write_fab_rows(os.path.join(dirpath, "red_posi"), [self.posi])
        write_fab_rows(os.path.join(dirpath, "red_ps"), self.ps)
        write_fab_rows(os.path.join(dirpath, "red_data"), self.data)
        if self.pid_raw is not None:
            open(os.path.join(dirpath, "red_pid"), "w").write(self.pid_raw)
        if self.cons_raw is not None:
            open(os.path.join(dirpath, "red_cons"), "w").write(self.cons_raw)


# --------------------------------------------------------------------------
# Kinematics: kinematics/<job>/  (reader: search/kinematics.c:14)
# --------------------------------------------------------------------------

@dataclass
class Kinematics:
    """kin_common: 'nmono npara'; kin_table: nmono x npara exponent matrix.
    Column order = BLSearchParameter = red_ps column order."""
    nmono: int
    npara: int
    table: List[List[int]] = field(default_factory=list)

    def validate(self) -> None:
        if len(self.table) != self.nmono:
            raise BladeFormatError(f"kin_table: {len(self.table)} rows != nmono={self.nmono}")
        if self.table and len(self.table[0]) != self.npara:
            raise BladeFormatError("kin_table: wrong column count")

    @classmethod
    def read(cls, dirpath: str) -> "Kinematics":
        nmono, npara = read_ints(os.path.join(dirpath, "kin_common"))[:2]
        kn = cls(nmono, npara,
                 _rows(read_ints(os.path.join(dirpath, "kin_table")), npara, "kin_table"))
        kn.validate()
        return kn

    def write(self, dirpath: str) -> None:
        self.validate()
        os.makedirs(dirpath, exist_ok=True)
        write_fab_rows(os.path.join(dirpath, "kin_common"), [[self.nmono, self.npara]])
        write_fab_rows(os.path.join(dirpath, "kin_table"), self.table)


# --------------------------------------------------------------------------
# Block scheme: <job>/<workid>/sch_*  (readers: search/template.c:15,
# ssolve/block.c:80 -- ssolve reads multi_sch_intid, not sch_intid)
# --------------------------------------------------------------------------

@dataclass
class Scheme:
    """sch_g1 / sch_nint / sch_intid (+ byte-copy multi_sch_intid for ssolve).
    intid: nint 0-based GLOBAL integral ids, G1 first then G2."""
    g1: int
    nint: int
    intid: List[int] = field(default_factory=list)

    def validate(self) -> None:
        if len(self.intid) != self.nint:
            raise BladeFormatError(f"sch_intid: {len(self.intid)} != nint={self.nint}")
        if not (0 < self.g1 <= self.nint):
            raise BladeFormatError(f"sch_g1={self.g1} outside (0, nint={self.nint}]")

    @classmethod
    def read(cls, dirpath: str) -> "Scheme":
        g1 = read_ints(os.path.join(dirpath, "sch_g1"))[0]
        nint = read_ints(os.path.join(dirpath, "sch_nint"))[0]
        sch = cls(g1, nint, read_ints(os.path.join(dirpath, "sch_intid")))
        multi = os.path.join(dirpath, "multi_sch_intid")
        if os.path.exists(multi) and read_ints(multi) != sch.intid:
            raise BladeFormatError(f"{multi} disagrees with sch_intid")
        sch.validate()
        return sch

    def write(self, dirpath: str) -> None:
        self.validate()
        os.makedirs(dirpath, exist_ok=True)
        write_fab_rows(os.path.join(dirpath, "sch_g1"), [[self.g1]])
        write_fab_rows(os.path.join(dirpath, "sch_nint"), [[self.nint]])
        write_fab_rows(os.path.join(dirpath, "sch_intid"), [self.intid])
        write_fab_rows(os.path.join(dirpath, "multi_sch_intid"), [self.intid])


# --------------------------------------------------------------------------
# Ansatz config (input): config/<k>/in_nvar,in_var  (template.c:105)
# --------------------------------------------------------------------------

@dataclass
class AnsatzConfig:
    """in_var: nvar x 2 rows {localIntegralID, kinTableRowID}, both 0-based."""
    var: List[Tuple[int, int]] = field(default_factory=list)

    @property
    def nvar(self) -> int:
        return len(self.var)

    @classmethod
    def read(cls, dirpath: str) -> "AnsatzConfig":
        nvar = read_ints(os.path.join(dirpath, "in_nvar"))[0]
        rows = _rows(read_ints(os.path.join(dirpath, "in_var")), 2, "in_var")
        if len(rows) != nvar:
            raise BladeFormatError(f"{dirpath}: in_var rows {len(rows)} != in_nvar {nvar}")
        return cls([tuple(r) for r in rows])

    def write(self, dirpath: str) -> None:
        os.makedirs(dirpath, exist_ok=True)
        write_fab_rows(os.path.join(dirpath, "in_nvar"), [[self.nvar]])
        write_fab_rows(os.path.join(dirpath, "in_var"), [list(r) for r in self.var])


# --------------------------------------------------------------------------
# redg1 outputs: config/<k>/out_*, tmp_*  (searchalg.c:474-510, template.c:229)
# --------------------------------------------------------------------------

@dataclass
class OutSolution:
    """out_nsol/out_nvar/out_sol/out_rel/out_var (template_solution_write,
    template.c:229-249).  out_sol: nsol x nvar nullspace; out_rel: nsol x nint
    relation values at the RESERVED point; out_var: nvar x 2 var list after
    absorb.  Written by C: scalars via file_write_int, tables via
    table_write_llu/int."""
    nsol: int
    nvar: int
    sol: List[List[int]] = field(default_factory=list)
    rel: List[List[int]] = field(default_factory=list)
    var: List[Tuple[int, int]] = field(default_factory=list)

    @classmethod
    def read(cls, dirpath: str, nint: int) -> "OutSolution":
        nsol = read_ints(os.path.join(dirpath, "out_nsol"))[0]
        nvar = read_ints(os.path.join(dirpath, "out_nvar"))[0]
        sol = _rows(read_ints(os.path.join(dirpath, "out_sol")), nvar, "out_sol")
        rel = _rows(read_ints(os.path.join(dirpath, "out_rel")), nint, "out_rel")
        var = [tuple(r) for r in _rows(read_ints(os.path.join(dirpath, "out_var")), 2, "out_var")]
        if len(sol) != nsol or len(rel) != nsol or len(var) != nvar:
            raise BladeFormatError(f"{dirpath}: out_* dimensions inconsistent "
                                   f"(nsol={nsol}, nvar={nvar}, sol={len(sol)}, "
                                   f"rel={len(rel)}, var={len(var)})")
        return cls(nsol, nvar, sol, rel, var)

    def write(self, dirpath: str) -> None:
        os.makedirs(dirpath, exist_ok=True)
        write_c_scalar(os.path.join(dirpath, "out_nsol"), self.nsol)
        write_c_scalar(os.path.join(dirpath, "out_nvar"), self.nvar)
        write_c_table(os.path.join(dirpath, "out_sol"), self.sol)
        write_c_table(os.path.join(dirpath, "out_rel"), self.rel)
        write_c_table(os.path.join(dirpath, "out_var"), [list(r) for r in self.var])


@dataclass
class TmpTemplate:
    """tmp_nvar/tmp_nsol/tmp_var/tmp_indep (write_template_g1, searchalg.c:474).
    tmp_indep has OUT_NSOL entries (0/1 flag per solution of the config's full
    nullspace, list_write_int); tmp_nsol = count of 1s; tmp_var = nvar' x 2
    pruned var list (table_write_int)."""
    nvar: int
    nsol: int
    var: List[Tuple[int, int]] = field(default_factory=list)
    indep: List[int] = field(default_factory=list)

    @classmethod
    def read(cls, dirpath: str) -> "TmpTemplate":
        nvar = read_ints(os.path.join(dirpath, "tmp_nvar"))[0]
        nsol = read_ints(os.path.join(dirpath, "tmp_nsol"))[0]
        var = [tuple(r) for r in _rows(read_ints(os.path.join(dirpath, "tmp_var")), 2, "tmp_var")]
        indep = read_ints(os.path.join(dirpath, "tmp_indep"))
        if len(var) != nvar:
            raise BladeFormatError(f"{dirpath}: tmp_var rows {len(var)} != tmp_nvar {nvar}")
        if sum(indep) != nsol:
            raise BladeFormatError(f"{dirpath}: sum(tmp_indep)={sum(indep)} != tmp_nsol={nsol}")
        return cls(nvar, nsol, var, indep)

    def write(self, dirpath: str) -> None:
        os.makedirs(dirpath, exist_ok=True)
        write_c_scalar(os.path.join(dirpath, "tmp_nvar"), self.nvar)
        write_c_scalar(os.path.join(dirpath, "tmp_nsol"), self.nsol)
        write_c_table(os.path.join(dirpath, "tmp_var"), [list(r) for r in self.var])
        write_c_list(os.path.join(dirpath, "tmp_indep"), self.indep)


# --------------------------------------------------------------------------
# fit tables: <work>/fit/<dataname>/<k>  (writer: searchalg.c:666 via
# table_write_llu; reader: ssolve/block.c:122 reads the FIRST tmp_nsol rows)
# --------------------------------------------------------------------------

@dataclass
class FitTable:
    """Nullspace of the pruned template on one database: rows x tmp_nvar llu.
    A config with tmp_nsol==0 legitimately has an EMPTY (0-byte) table.
    MEASURED FOOTGUN: fitrel can write flag=1 next to an EMPTY table whose
    tmp_nsol>0 -- always gate on dimensions, never the flag."""
    ncol: int
    rows: List[List[int]] = field(default_factory=list)

    @property
    def nrow(self) -> int:
        return len(self.rows)

    @classmethod
    def read(cls, path: str, ncol: int) -> "FitTable":
        flat = read_ints(path)
        if ncol == 0:
            if flat:
                raise BladeFormatError(f"{path}: ncol=0 but nonempty")
            return cls(0, [])
        return cls(ncol, _rows(flat, ncol, path))

    def write(self, path: str) -> None:
        write_c_table(path, self.rows)


# --------------------------------------------------------------------------
# ssolve system txt (writer: SolveSemiBL.wl collectSystemInfo:325;
# reader: ssolve/block.c:778 system_init_from_txt_parallel)
# --------------------------------------------------------------------------

@dataclass
class SystemFile:
    """'prime npara nint nmaster ntarget t_1..t_n nfiles' + ONE whitespace-free
    token 'blockdir_1,...,blockdir_n,dataname'.  Kinematics are found at
    blockdir_1/../kinematics/<basename(blockdir_1)> (block.c:849), so the job
    dir basename must equal the kinematics subdir name.  NOTE (block.c:788):
    field 2 is read into sys->kin_row and later overwritten by kin_common's
    nmono -- the system-file npara is effectively ignored."""
    prime: int
    npara: int
    nint: int
    nmaster: int
    targets: List[int] = field(default_factory=list)
    blockdirs: List[str] = field(default_factory=list)
    dataname: str = ""

    @classmethod
    def read(cls, path: str) -> "SystemFile":
        toks = read_tokens(path)
        prime, npara, nint, nmaster, ntar = (int(t) for t in toks[:5])
        targets = [int(t) for t in toks[5:5 + ntar]]
        nfiles = int(toks[5 + ntar])
        token = toks[6 + ntar]
        if len(toks) != 7 + ntar:
            raise BladeFormatError(f"{path}: trailing tokens beyond the dir token")
        parts = token.split(",")
        if len(parts) != nfiles + 1:
            raise BladeFormatError(f"{path}: token has {len(parts)-1} dirs, header says {nfiles}")
        return cls(prime, npara, nint, nmaster, targets, parts[:-1], parts[-1])

    def write(self, path: str) -> None:
        head = [self.prime, self.npara, self.nint, self.nmaster,
                len(self.targets)] + list(self.targets) + [len(self.blockdirs)]
        with open(path, "w") as f:
            f.write(" ".join(str(x) for x in head))
            f.write("\n" + ",".join(list(self.blockdirs) + [self.dataname]) + "\n")

    def kinematics_dir(self) -> str:
        b = self.blockdirs[0]
        return os.path.join(os.path.dirname(b), "kinematics", os.path.basename(b))


@dataclass
class EvalList:
    """recmod eval_list_file: ONE whitespace-free comma-separated token of
    evaluation-file paths (fflowC.cc:527-545, 1MB cap).  Trailing whitespace
    preserved verbatim for byte-identical round trips."""
    paths: List[str] = field(default_factory=list)
    tail: str = ""

    @classmethod
    def read(cls, path: str) -> "EvalList":
        text = open(path).read()
        body = text.rstrip()
        return cls(body.split(",") if body else [], text[len(body):])

    def write(self, path: str) -> None:
        with open(path, "w") as f:
            f.write(",".join(self.paths))
            f.write(self.tail)


# --------------------------------------------------------------------------
# points file (binary; fflow load_samples, alg_reconstruction.cc:331)
# --------------------------------------------------------------------------

def flags_size(nparsout: int) -> int:
    return (nparsout + 63) // 64


@dataclass
class PointsFile:
    """[n, nsamples] then nsamples rows of n+1 words, n = nparsin + flags_size:
    row = x_1..x_nparsin, prime, flags words (LSB-first bit j = 'output j
    needed').  Points are fflow-internal STRUCTURED sample points -- never
    fabricate them freehand; generate via the dumppoints shim."""
    n: int
    rows: List[List[int]] = field(default_factory=list)

    @classmethod
    def read(cls, path: str) -> "PointsFile":
        w = read_u64_words(path)
        if len(w) < 2:
            raise BladeFormatError(f"{path}: too short for a points header")
        n, nsamples = w[0], w[1]
        body = w[2:]
        if len(body) != nsamples * (n + 1):
            raise BladeFormatError(f"{path}: {len(body)} words != nsamples*(n+1) "
                                   f"= {nsamples}*({n}+1)")
        return cls(n, _rows(body, n + 1, path))

    def write(self, path: str) -> None:
        words = [self.n, len(self.rows)]
        for r in self.rows:
            if len(r) != self.n + 1:
                raise BladeFormatError("points row width != n+1")
            words += list(r)
        write_u64_words(path, words)

    def split_row(self, i: int, nparsin: int) -> Tuple[List[int], int, List[int]]:
        """-> (coords, prime, flags_words) for sample i."""
        r = self.rows[i]
        return r[:nparsin], r[nparsin], r[nparsin + 1:]


# --------------------------------------------------------------------------
# eval file (binary; system_dump_evaluations, iofflow.cpp:205; also
# finiteflow's evaluation-cache dump format)
# --------------------------------------------------------------------------

@dataclass
class EvalFile:
    """Header [nparsin, nparsout+flags_size, nrows, 0]; per row:
    x_1..x_nparsin, prime, flags words, then the NEEDED outputs compacted in
    order, zero-padded so the row width is nparsin+1+flags_size+nparsout."""
    nparsin: int
    nout_plus_flags: int   # header word 2 = nparsout + flags_size
    rows: List[List[int]] = field(default_factory=list)

    @property
    def row_width(self) -> int:
        return self.nparsin + 1 + self.nout_plus_flags

    @classmethod
    def expected_bytes(cls, nparsin: int, nparsout: int, nrows: int) -> int:
        fs = flags_size(nparsout)
        return 8 * (4 + nrows * (nparsin + 1 + fs + nparsout))

    @classmethod
    def read(cls, path: str) -> "EvalFile":
        w = read_u64_words(path)
        if len(w) < 4:
            raise BladeFormatError(f"{path}: too short for an eval header")
        nparsin, nout_pf, nrows, zero = w[:4]
        if zero != 0:
            raise BladeFormatError(f"{path}: header word 4 = {zero} != 0")
        width = nparsin + 1 + nout_pf
        body = w[4:]
        if len(body) != nrows * width:
            raise BladeFormatError(f"{path}: {len(body)} words != nrows*width = "
                                   f"{nrows}*{width}")
        return cls(nparsin, nout_pf, _rows(body, width, path))

    def write(self, path: str) -> None:
        words = [self.nparsin, self.nout_plus_flags, len(self.rows), 0]
        for r in self.rows:
            if len(r) != self.row_width:
                raise BladeFormatError("eval row width mismatch")
            words += list(r)
        write_u64_words(path, words)


# --------------------------------------------------------------------------
# degrees file (binary; algorithm_dump_degree_info,
# finiteflow/src/alg_mp_reconstruction.cc:99-131)
# --------------------------------------------------------------------------

@dataclass
class DegreeInfo:
    numdeg: int                                  # numerator total degree
    dendeg: int                                  # denominator total degree
    # per var (points-file coordinate order): (num_max, num_min, den_max, den_min)
    var: List[Tuple[int, int, int, int]] = field(default_factory=list)


@dataclass
class DegreesFile:
    """[nparsin, nparsout] then per output: numdeg_tot, dendeg_tot, then per
    var 4 words (num_maxdeg, num_mindeg, den_maxdeg, den_mindeg) -- dump order
    per alg_mp_reconstruction.cc:119-127."""
    nparsin: int
    nparsout: int
    info: List[DegreeInfo] = field(default_factory=list)

    @classmethod
    def read(cls, path: str) -> "DegreesFile":
        w = read_u64_words(path)
        if len(w) < 2:
            raise BladeFormatError(f"{path}: too short for a degrees header")
        nin, nout = w[0], w[1]
        need = 2 + nout * (2 + 4 * nin)
        if len(w) != need:
            raise BladeFormatError(f"{path}: {len(w)} words != {need} "
                                   f"(nparsin={nin}, nparsout={nout})")
        info, k = [], 2
        for _ in range(nout):
            numdeg, dendeg = w[k], w[k + 1]
            k += 2
            var = []
            for _ in range(nin):
                var.append(tuple(w[k:k + 4]))
                k += 4
            info.append(DegreeInfo(numdeg, dendeg, var))
        return cls(nin, nout, info)

    def write(self, path: str) -> None:
        if len(self.info) != self.nparsout:
            raise BladeFormatError("degrees: len(info) != nparsout")
        words = [self.nparsin, self.nparsout]
        for d in self.info:
            words += [d.numdeg, d.dendeg]
            if len(d.var) != self.nparsin:
                raise BladeFormatError("degrees: per-var record count != nparsin")
            for v in d.var:
                words += list(v)
        write_u64_words(path, words)

    @classmethod
    def uniform(cls, nparsin: int, nparsout: int, maxdeg: int) -> "DegreesFile":
        """Uniform bounds: total and per-var max degrees = maxdeg, mins = 0.
        WARNING (measured): NOT a valid recmod input -- recmod
        fails when the declared TOTAL degrees over-estimate the true ones
        (per-var over-estimates are tolerated).  Use exact degrees or
        ratrec.scan_degrees().  Kept for tests/fabrication only."""
        info = [DegreeInfo(maxdeg, maxdeg, [(maxdeg, 0, maxdeg, 0)] * nparsin)
                for _ in range(nparsout)]
        return cls(nparsin, nparsout, info)


# --------------------------------------------------------------------------
# recmod outputs rec_*_{part,mono,coeff} and dynamicrr rrres
# (writers: fflowC.cc put_sparse_ratfun_* :119-153, put_mprat :155-165)
# --------------------------------------------------------------------------

@dataclass
class RecPart:
    """<prefix>_part: per function 'nnum nden' (>=1 each; an empty poly is
    written as ONE placeholder term: exponents all 0, coefficient '0')."""
    pairs: List[Tuple[int, int]] = field(default_factory=list)

    @property
    def total_terms(self) -> int:
        return sum(a + b for a, b in self.pairs)

    @classmethod
    def read(cls, path: str) -> "RecPart":
        flat = read_ints(path)
        if len(flat) % 2:
            raise BladeFormatError(f"{path}: odd token count {len(flat)}")
        pairs = [(flat[i], flat[i + 1]) for i in range(0, len(flat), 2)]
        if any(a < 1 or b < 1 for a, b in pairs):
            raise BladeFormatError(f"{path}: part entry < 1")
        return cls(pairs)

    def write(self, path: str) -> None:
        write_tokens(path, [x for p in self.pairs for x in p])


@dataclass
class RecMono:
    """<prefix>_mono: per function, numerator then denominator monomial
    exponent tuples (nparsin uints each).  MUST be identical across primes."""
    nparsin: int
    exps: List[Tuple[int, ...]] = field(default_factory=list)

    @classmethod
    def read(cls, path: str, nparsin: int) -> "RecMono":
        flat = read_ints(path)
        rows = _rows(flat, nparsin, path)
        return cls(nparsin, [tuple(r) for r in rows])

    def write(self, path: str) -> None:
        write_tokens(path, [e for t in self.exps for e in t])


@dataclass
class RecCoeff:
    """<prefix>_coeff: one decimal token per monomial (residues mod the prime;
    '0' for an empty poly's placeholder).  Kept as strings for byte fidelity;
    values() parses them (mpq_get_str may emit p/q in other fflow contexts)."""
    tokens: List[str] = field(default_factory=list)

    @classmethod
    def read(cls, path: str) -> "RecCoeff":
        return cls(read_tokens(path))

    def write(self, path: str) -> None:
        write_tokens(path, self.tokens)

    def values(self) -> List[Fraction]:
        return [Fraction(t) for t in self.tokens]


@dataclass
class RRRes:
    """dynamicrr out_file: space-separated mpq decimals 'p/q' or 'p'
    (put_mprat writes '0 ' for zero), one per coefficient, in rec_*_coeff
    order."""
    values: List[Fraction] = field(default_factory=list)

    @classmethod
    def read(cls, path: str) -> "RRRes":
        return cls([Fraction(t) for t in read_tokens(path)])

    def write(self, path: str) -> None:
        write_tokens(path, [_frac_str(v) for v in self.values])


def _frac_str(v: Fraction) -> str:
    return str(v.numerator) if v.denominator == 1 else f"{v.numerator}/{v.denominator}"


# --------------------------------------------------------------------------
# assembly: rrres/coeff + mono + part -> rational functions
# --------------------------------------------------------------------------

def assemble_rational_functions(part: RecPart, mono: RecMono,
                                values: Sequence, symbols=None):
    """Rebuild the reconstructed rational functions as sympy expressions.

    values: one coefficient per monomial in part/mono order -- Fractions from
    RRRes (exact) or residues from RecCoeff (single-prime image).
    Returns a list of sympy expressions (len == len(part.pairs))."""
    import sympy as sp
    n = mono.nparsin
    if symbols is None:
        symbols = sp.symbols(f"x1:{n + 1}") if n > 1 else (sp.Symbol("x1"),)
    if len(symbols) != n:
        raise BladeFormatError(f"need {n} symbols, got {len(symbols)}")
    if len(values) != part.total_terms or len(mono.exps) != part.total_terms:
        raise BladeFormatError(
            f"assembly: {len(values)} values / {len(mono.exps)} monomials "
            f"!= sum(part) = {part.total_terms}")
    out, k = [], 0

    def poly(nterms):
        nonlocal k
        acc = sp.Integer(0)
        for _ in range(nterms):
            c = sp.Rational(Fraction(values[k]))
            term = c
            for s, e in zip(symbols, mono.exps[k]):
                if e:
                    term *= s ** int(e)
            acc += term
            k += 1
        return acc

    for nnum, nden in part.pairs:
        num = poly(nnum)
        den = poly(nden)
        if den == 0:
            raise BladeFormatError("assembly: zero denominator polynomial")
        out.append(num / den)
    return out


# --------------------------------------------------------------------------
# manifest.json (phase0 bookkeeping; json, indent=1 as the fabricator wrote it)
# --------------------------------------------------------------------------

@dataclass
class Manifest:
    data: dict = field(default_factory=dict)

    @classmethod
    def read(cls, path: str) -> "Manifest":
        return cls(json.load(open(path)))

    def write(self, path: str) -> None:
        with open(path, "w") as f:
            json.dump(self.data, f, indent=1)


# --------------------------------------------------------------------------
# generic byte-identical round-trip machinery
# --------------------------------------------------------------------------

@dataclass
class _Spec:
    kind: str
    parse: Callable[[str], object]     # path -> obj
    emit: Callable[[object, str], None]  # (obj, path) -> None


def _fab_onerow(path):
    return read_ints(path)


def _emit_fab_onerow(obj, path):
    write_fab_rows(path, [obj])


def _fab_table(ncol_of: Callable[[str], int]):
    def parse(path):
        return _rows(read_ints(path), ncol_of(path), path)
    def emit(obj, path):
        write_fab_rows(path, obj)
    return parse, emit


def _c_table(ncol_of: Callable[[str], int]):
    def parse(path):
        return _rows(read_ints(path), ncol_of(path), path)
    def emit(obj, path):
        write_c_table(path, obj)
    return parse, emit


def _sibling_int(path: str, *rel: str) -> int:
    return read_ints(os.path.join(os.path.dirname(path), *rel))[0]


def _dbdir_field(path: str, idx: int) -> int:
    return read_ints(os.path.join(os.path.dirname(path), "red_common"))[idx]


def _config_workdir(path: str) -> str:
    # <work>/config/<k>/xxx -> <work>
    return os.path.dirname(os.path.dirname(os.path.dirname(path)))


def _fit_tmp_nvar(path: str) -> int:
    # <work>/fit/<db>/<k> -> <work>/config/<k>/tmp_nvar
    k = os.path.basename(path)
    work = os.path.dirname(os.path.dirname(os.path.dirname(path)))
    return read_ints(os.path.join(work, "config", k, "tmp_nvar"))[0]


def _raw_parse(path):
    return open(path, "rb").read()


def _raw_emit(obj, path):
    open(path, "wb").write(obj)


def classify(path: str) -> Optional[_Spec]:
    """Map a phase0/production artifact path to its format spec, or None if
    the file is not a pipeline contract file (logs, scripts, reports)."""
    base = os.path.basename(path)
    parent = os.path.basename(os.path.dirname(path))
    gparent = os.path.basename(os.path.dirname(os.path.dirname(path)))

    one_row = {"red_common", "red_part", "red_posi", "red_pid", "kin_common",
               "sch_g1", "sch_nint", "sch_intid", "multi_sch_intid", "in_nvar"}
    if base in one_row:
        return _Spec(base, _fab_onerow, _emit_fab_onerow)
    if base == "red_cons":
        return _Spec(base, _raw_parse, _raw_emit)
    if base == "red_ps":
        p, e = _fab_table(lambda q: _dbdir_field(q, 1))
        return _Spec(base, p, e)
    if base == "red_data":
        p, e = _fab_table(lambda q: _dbdir_field(q, 4))
        return _Spec(base, p, e)
    if base == "kin_table":
        p, e = _fab_table(
            lambda q: read_ints(os.path.join(os.path.dirname(q), "kin_common"))[1])
        return _Spec(base, p, e)
    if base == "in_var":
        p, e = _fab_table(lambda q: 2)
        return _Spec(base, p, e)
    if base in {"tmp_nvar", "tmp_nsol", "out_nsol", "out_nvar", "flag", "state"}:
        return _Spec(base, lambda q: read_ints(q)[0],
                     lambda o, q: write_c_scalar(q, o))
    if base == "tmp_indep":
        return _Spec(base, read_ints, lambda o, q: write_c_list(q, o))
    if base in {"tmp_var", "out_var"}:
        p, e = _c_table(lambda q: 2)
        return _Spec(base, p, e)
    if base == "out_sol":
        p, e = _c_table(lambda q: _sibling_int(q, "out_nvar"))
        return _Spec(base, p, e)
    if base == "out_rel":
        p, e = _c_table(lambda q: read_ints(
            os.path.join(_config_workdir(q), "sch_nint"))[0])
        return _Spec(base, p, e)
    if base.isdigit() and gparent == "fit":
        p, e = _c_table(_fit_tmp_nvar)
        return _Spec("fit_table", p, e)
    if base.startswith("system_") and base.endswith(".txt"):
        return _Spec("system", SystemFile.read, lambda o, q: o.write(q))
    if base.startswith("evallist") and base.endswith(".txt"):
        return _Spec("evallist", EvalList.read, lambda o, q: o.write(q))
    if base.endswith(".fflow"):
        if "points" in base:
            return _Spec("points", PointsFile.read, lambda o, q: o.write(q))
        if "eval" in base or parent == "eval":
            return _Spec("eval", EvalFile.read, lambda o, q: o.write(q))
        if "degree" in base:
            return _Spec("degrees", DegreesFile.read, lambda o, q: o.write(q))
        return None
    if base.endswith("_part"):
        return _Spec("rec_part", RecPart.read, lambda o, q: o.write(q))
    if base.endswith("_mono"):
        def parse_mono(q):
            deg = os.path.join(os.path.dirname(q), "degrees.fflow")
            nin = read_u64_words(deg)[0] if os.path.exists(deg) else 1
            return RecMono.read(q, nin)
        return _Spec("rec_mono", parse_mono, lambda o, q: o.write(q))
    if base.endswith("_coeff"):
        return _Spec("rec_coeff", RecCoeff.read, lambda o, q: o.write(q))
    if base == "rrres" or base.startswith("rrres"):
        return _Spec("rrres", RRRes.read, lambda o, q: o.write(q))
    if base == "manifest.json":
        return _Spec("manifest", Manifest.read, lambda o, q: o.write(q))
    return None


def round_trip(path: str, scratch_dir: str) -> Tuple[bool, str]:
    """Parse `path` with its format reader, re-write it with the format
    writer into scratch_dir (NEVER in place), byte-compare.
    -> (ok, detail).  Raises BladeFormatError if the file is unclassifiable."""
    spec = classify(path)
    if spec is None:
        raise BladeFormatError(f"{path}: not a known contract file")
    obj = spec.parse(path)
    os.makedirs(scratch_dir, exist_ok=True)
    out = os.path.join(scratch_dir, "rt_" + os.path.basename(path))
    spec.emit(obj, out)
    a = open(path, "rb").read()
    b = open(out, "rb").read()
    os.unlink(out)
    if a == b:
        return True, f"{spec.kind}: {len(a)}B byte-identical"
    return False, (f"{spec.kind}: MISMATCH orig {len(a)}B vs rewrite {len(b)}B")
