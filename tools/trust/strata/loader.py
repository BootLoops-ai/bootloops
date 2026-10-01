#!/usr/bin/env python3
"""loader.py — SYSTEM_*.gz -> F_p loader at numeric (d, eta).

Parses kira's tmp/<family>/SYSTEM_*.gz selected-equation dumps (format observed
on vac3/k2sh Generates, kira 3.x, integral_ordering 5):

    Eq
    <eq id (weight-like)>
    <n_terms>
    <n_terms> <coeff-expr> <weight> <sector> 0 0     (x n_terms)

Coefficient expressions are polynomials/rationals in (d, eta) with ^ powers.
They are evaluated mod p at a NUMERIC slice (d0, eta0) via a restricted eval
into F_p elements (compiled once per distinct string; kira emits heavy reuse).

Weight semantics per the banked empirical decode (weights.py, validated on
production k2sh artifacts and re-checked against vac3 kira.db
WEIGHTBITS=(5,4,17,13), INTEGRALORDERING=5):
  - w >= 2^30 : ordinary integral; (dots, s) = decode(w); sector = explicit
                sector field of the term line.
  - w <  2^30 : master-class ordinal; mapped to an integral via the run's
                results/<family>/masters file (matched BY SECTOR; asserted
                unique per sector).
Terms whose sector is in sectormappings/<family>/trivialsector are ZERO
integrals and are dropped (with a counter).

Output of load_system(): dict with
  rows        list[dict weight->F_p int]  (nonempty rows only)
  order       dict weight->weight  (kira weight IS the Laporta rank)
  stratum_of  dict weight->t (popcount of sector)
  sector_of   dict weight->sector
  class_of    dict weight->(dots, s) or None for master-class
  masters     dict weight->{"sector": int, "indices": [..] or None}
  counters    parse statistics (eqs, terms, dropped_trivial, zero_coeff, ...)
"""
import glob
import gzip
import os
import re

from weights import MASTER_WEIGHT_CEILING, decode, popcount


class Fp:
    """Minimal F_p element supporting the operators kira coeff strings use."""
    __slots__ = ("v", "p")

    def __init__(self, v, p):
        self.v = v % p
        self.p = p

    def _c(self, o):
        return o.v if isinstance(o, Fp) else o % self.p

    def __add__(self, o):
        return Fp(self.v + self._c(o), self.p)
    __radd__ = __add__

    def __sub__(self, o):
        return Fp(self.v - self._c(o), self.p)

    def __rsub__(self, o):
        return Fp(self._c(o) - self.v, self.p)

    def __mul__(self, o):
        return Fp(self.v * self._c(o), self.p)
    __rmul__ = __mul__

    def __truediv__(self, o):
        ov = self._c(o)
        if ov % self.p == 0:
            raise ZeroDivisionError("denominator hit mod p")
        return Fp(self.v * pow(ov, self.p - 2, self.p), self.p)

    def __rtruediv__(self, o):
        if self.v == 0:
            raise ZeroDivisionError("denominator hit mod p")
        return Fp(self._c(o) * pow(self.v, self.p - 2, self.p), self.p)

    def __pow__(self, n):
        return Fp(pow(self.v, int(n), self.p), self.p)

    def __neg__(self):
        return Fp(-self.v, self.p)

    def __int__(self):          # exponents arrive as Fp after literal-wrapping
        return self.v


_ALLOWED = re.compile(r"^[0-9deta+\-*/^() ]+$")


class CoeffEvaluator:
    """Compiled-eval cache for kira coefficient strings at one (p, d0, eta0)."""

    def __init__(self, p, d0, eta0):
        self.p = p
        self.env = {"d": Fp(d0, p), "eta": Fp(eta0, p),
                    "F": lambda v, _p=p: Fp(v, _p), "__builtins__": {}}
        self.cache = {}

    def __call__(self, expr):
        code = self.cache.get(expr)
        if code is None:
            if not _ALLOWED.match(expr):
                raise ValueError(f"unexpected coeff chars: {expr!r}")
            # Wrap every integer literal as an F_p element BEFORE compiling:
            # pure-numeric coefficients like (-1)/(2) otherwise evaluate with
            # Python float division (k2disp finding — exactness bug;
            # vac3 coeffs always carried d/eta, masking it).
            wrapped = re.sub(r"(\d+)", r"F(\1)", expr.replace("^", "**"))
            code = compile(wrapped, "<coeff>", "eval")
            self.cache[expr] = code
        val = eval(code, self.env)  # noqa: S307 (restricted charset + no builtins)
        assert isinstance(val, Fp), f"non-Fp coefficient value for {expr!r}"
        return val.v


def parse_masters_file(path):
    """masters/preferred file -> [(indices tuple, sector), ...] (file order)."""
    out = []
    for line in open(path):
        m = re.search(r"\[([-0-9, ]+)\]", line)
        if not m:
            continue
        idx = tuple(int(x) for x in m.group(1).split(","))
        sec = sum(1 << i for i, a in enumerate(idx) if a > 0)
        out.append((idx, sec))
    return out


def load_trivial_sectors(art_dir, family):
    path = os.path.join(art_dir, "sectormappings", family, "trivialsector")
    triv = set()
    if os.path.exists(path):
        for tok in re.findall(r"\d+", open(path).read()):
            triv.add(int(tok))
    return triv


def load_system(art_dir, family, p, d0, eta0, drop_trivial=True):
    ev = CoeffEvaluator(p, d0, eta0)
    trivial = load_trivial_sectors(art_dir, family) if drop_trivial else set()
    files = sorted(glob.glob(os.path.join(art_dir, "tmp", family, "SYSTEM_*.gz")))
    if not files:
        raise FileNotFoundError(f"no SYSTEM_*.gz under {art_dir}/tmp/{family}")

    rows = []
    sector_of, class_of = {}, {}
    counters = {"eqs": 0, "terms": 0, "dropped_trivial": 0, "zero_coeff": 0,
                "cancelled_to_empty": 0, "files": len(files)}
    for f in files:
        with gzip.open(f, "rt") as fh:
            lines = fh.read().split("\n")
        i = 0
        n = len(lines)
        while i < n:
            if lines[i] != "Eq":
                i += 1
                continue
            n_terms = int(lines[i + 2])
            row = {}
            for j in range(n_terms):
                parts = lines[i + 3 + j].split()
                assert len(parts) == 6, f"bad term line: {lines[i + 3 + j]!r}"
                coeff, w, sec = parts[1], int(parts[2]), int(parts[3])
                counters["terms"] += 1
                if sec in trivial:
                    counters["dropped_trivial"] += 1
                    continue
                cv = ev(coeff)
                if cv == 0:
                    counters["zero_coeff"] += 1
                    continue
                row[w] = (row.get(w, 0) + cv) % p
                if row[w] == 0:
                    del row[w]
                sector_of[w] = sec
                class_of[w] = decode(w)
            i += 3 + n_terms
            counters["eqs"] += 1
            row = {w: v for w, v in row.items() if v}
            if row:
                rows.append(row)
            else:
                counters["cancelled_to_empty"] += 1

    # masters: master-class weights PLUS file-masters resolved to ordinary
    # weights. Empirics (measured on vac3-family sweep Generates):
    # kira grants small ordinals to the PREFERRED masters, in preferred-file
    # order (ordinal rank k <-> k-th preferred entry; sector asserted);
    # every other masters-file entry keeps an ordinary weight, resolved here
    # by (sector, (dots,s)) class. Any wrong assignment is caught by the
    # held-out dictionary match / full-table cross-check downstream.
    masters = {}
    file_masters = []
    mfile = os.path.join(art_dir, "results", family, "masters")
    if os.path.exists(mfile):
        file_masters = parse_masters_file(mfile)
    preferred = []
    pfile = os.path.join(art_dir, "preferred")
    if os.path.exists(pfile):
        preferred = parse_masters_file(pfile)
    mclass_ws = sorted(w for w in sector_of if w < MASTER_WEIGHT_CEILING)
    if preferred and len(mclass_ws) == len(preferred):
        for w, (idx, sec) in zip(mclass_ws, preferred):
            assert sector_of[w] == sec, \
                f"preferred-order master map broken: w={w} sec {sector_of[w]} != {sec}"
            masters[w] = {"sector": sec, "indices": idx}
    else:               # fallback: unique-sector master-class weights
        by_sec = {}
        for idx, sec in file_masters:
            by_sec.setdefault(sec, []).append(idx)
        for w in mclass_ws:
            sec = sector_of[w]
            idxs = by_sec.get(sec, [])
            masters[w] = {"sector": sec,
                          "indices": idxs[0] if len(idxs) == 1 else None}
    # Non-preferred masters keep ORDINARY weights and are never resolved
    # here: a true master can never be pivoted on (a row leading on a master
    # would be a relation among masters), so they emerge as elimination
    # SURVIVORS and are labeled post-solve (corpus_bench) by (sector, class).

    # Initiate-census pseudo-master demotion (the
    # "+N false masters" class seen from the STAGING side): kira grants
    # small ordinals to EVERY preferred-file entry, but the run's own
    # results/<fam>/masters is the FINAL basis — preferred entries absent
    # from it were reduced by kira during initiate (relations exist in the
    # dumped system). Keeping them as forbidden masters would leave their
    # columns un-eliminable and poison closed-row tails (measured:
    # preferred 87 vs final 54). Demote them to ordinary columns;
    # the run's final masters file is the authority. No-op whenever
    # preferred is a subset of the final basis (vac3,
    # k2disp, lp1disp, sh: preferred==final or preferred subset).
    demoted = []
    if file_masters:
        final_idx = {tuple(idx) for idx, _ in file_masters}
        demoted = [w for w, m in masters.items()
                   if m["indices"] is not None
                   and tuple(m["indices"]) not in final_idx]
        for w in demoted:
            del masters[w]
        counters["pseudo_masters_demoted"] = len(demoted)

    order = {w: w for w in sector_of if w not in masters}
    # Demoted pseudo-masters keep their natural TINY kira ordinals as the
    # elimination order (probe-measured): rows that reduce to
    # {pseudo, masters} pass DOWN to the pseudo's low stratum and pivot on
    # it there (measured: the stratified solve completes this way). The one
    # invariant they break is WEIGHT-monotone closure order — a pseudo
    # pivot row may reference lower-stratum ordinary pivots that close
    # later in a plain ascending sweep. The cure lives in the closure
    # sweeps (coverage.assemble_closed / certify.close_all), which are
    # dependency-driven, not order-driven. (An order-remap above 2^62 was
    # tried and REJECTED: it diverts mixed rows into low strata and breaks
    # the pass-down direction — loud "pass-down row not below current
    # stratum" assert.)
    stratum_of = {w: popcount(s) for w, s in sector_of.items()}
    return {"rows": rows, "order": order, "stratum_of": stratum_of,
            "sector_of": sector_of, "class_of": class_of, "masters": masters,
            "file_masters": file_masters, "counters": counters}


def load_rows_filtered(art_dir, family, p, d0, eta0, include_sigs):
    """Parse SYSTEM_*.gz keeping ONLY equations whose content signature
    (frozenset of (coeff-str, weight)) is in include_sigs; evaluate at the
    slice like load_system. Returns list of dict rows (trivial dropped)."""
    ev = CoeffEvaluator(p, d0, eta0)
    trivial = load_trivial_sectors(art_dir, family)
    rows = []
    for f in sorted(glob.glob(os.path.join(art_dir, "tmp", family,
                                           "SYSTEM_*.gz"))):
        with gzip.open(f, "rt") as fh:
            lines = fh.read().split("\n")
        i = 0
        while i < len(lines):
            if lines[i] != "Eq":
                i += 1
                continue
            n_terms = int(lines[i + 2])
            sig, terms = [], []
            for j in range(n_terms):
                parts = lines[i + 3 + j].split()
                sig.append((parts[1], int(parts[2])))
                terms.append((parts[1], int(parts[2]), int(parts[3])))
            i += 3 + n_terms
            if frozenset(sig) not in include_sigs:
                continue
            row = {}
            for coeff, w, sec in terms:
                if sec in trivial:
                    continue
                cv = ev(coeff)
                if cv:
                    row[w] = (row.get(w, 0) + cv) % p
                    if row[w] == 0:
                        del row[w]
            if row:
                rows.append(row)
    return rows


def eq_signature(art_dir, family):
    """Set of content signatures (frozenset of (coeff-str, weight)) of all
    selected equations — used to identify FRESH identities in a bigger box."""
    sigs = set()
    for f in sorted(glob.glob(os.path.join(art_dir, "tmp", family, "SYSTEM_*.gz"))):
        with gzip.open(f, "rt") as fh:
            lines = fh.read().split("\n")
        i = 0
        while i < len(lines):
            if lines[i] != "Eq":
                i += 1
                continue
            n_terms = int(lines[i + 2])
            sig = []
            for j in range(n_terms):
                parts = lines[i + 3 + j].split()
                sig.append((parts[1], int(parts[2])))
            sigs.add(frozenset(sig))
            i += 3 + n_terms
    return sigs
