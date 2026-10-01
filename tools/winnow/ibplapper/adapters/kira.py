"""adapters/kira.py — Kira SYSTEM_*.gz sidecar adapter (Route A):
read a Kira Generate+Select artifact tree, evaluate at
one numeric (p, d0, eta0) slice, hand the core an opaque-column System.

ALL kira specifics live here, never in the core: artifact-tree layout
(tmp/<fam>/SYSTEM_*.gz, results/<fam>/masters, preferred,
sectormappings/<fam>/trivialsector), the weight decode, the ordinal-master
mapping, and the popcount(sector) default stratification.

HONESTY FLAGS:
  - SYSTEM_*.gz is an UNDOCUMENTED internal kira dump. The weight decode is
    EMPIRICAL and version-pinned: kira 3.1 / kira.db VERSION '2.1' /
    INTEGRALORDERING 5 / WEIGHTBITS (5,4,17,13) (see adapters/weights.py
    for the validated layout). A kira
    release can break it silently -> the version guard below REFUSES a
    kira.db whose WEIGHTBITS/ordering differ, and records "ABSENT-UNGUARDED"
    when no kira.db is present.
  - Kira never exports its weight<->integral dictionary. The ordinal-master
    mapping (small ordinals <-> preferred-file order) is asserted empirics;
    meta labels are "best-effort, receipt-checked", not guaranteed. The
    downstream cross-check (oracle table / receipt verify) is MANDATORY
    before any label is trusted.
  - Provenance gate: pass expect_counters={"eqs":..., "terms":...} to
    refuse a silently-different Generate (eq/term counts are the cheap
    invariant of a staged box).

Integer literals in coefficient strings are wrapped as F_p elements
BEFORE compile — pure-numeric coefficients like (-1)/(2) otherwise
evaluate with Python float division (a silent exactness bug).
"""
import glob
import gzip
import os
import re
import sqlite3

from .weights import MASTER_WEIGHT_CEILING, decode, popcount

# ------------------------------------------------------------ version guard
PINNED = {"WEIGHTBITS": (5, 4, 17, 13), "INTEGRALORDERING": 5,
          "DB_VERSION": "2.1"}


class KiraVersionError(RuntimeError):
    """kira.db declares a weight packing this adapter's EMPIRICAL decode was
    never validated on — refusing to guess (design §3 guard)."""


class ProvenanceError(RuntimeError):
    """Loaded eq/term counts differ from the declared expectation."""


def version_guard(art_dir):
    """Check results/kira.db against the pinned empirical-decode constants.

    Returns a guard record (stored in counters); RAISES KiraVersionError on
    a mismatching db. A missing kira.db is recorded as ABSENT-UNGUARDED —
    loud in the record, not fatal (some stagepacks ship without the db)."""
    db = os.path.join(art_dir, "results", "kira.db")
    if not os.path.exists(db):
        return {"status": "ABSENT-UNGUARDED", "db": db,
                "note": "no kira.db; empirical decode UNVERIFIED for this "
                        "artifact tree"}
    con = sqlite3.connect(db)
    try:
        wb = tuple(con.execute("SELECT A,B,C,D FROM WEIGHTBITS").fetchone())
        io = con.execute("SELECT ID FROM INTEGRALORDERING").fetchone()[0]
        ver = con.execute("SELECT number FROM VERSION").fetchone()[0]
    except sqlite3.Error as e:
        raise KiraVersionError(f"{db}: cannot read guard tables ({e})")
    finally:
        con.close()
    if wb != PINNED["WEIGHTBITS"] or int(io) != PINNED["INTEGRALORDERING"]:
        raise KiraVersionError(
            f"{db}: WEIGHTBITS={wb} INTEGRALORDERING={io} != pinned "
            f"{PINNED['WEIGHTBITS']}/{PINNED['INTEGRALORDERING']} — the "
            "empirical weight decode is NOT validated for this kira; "
            "refusing (re-run the probe battery before unpinning)")
    rec = {"status": "GUARDED", "WEIGHTBITS": list(wb),
           "INTEGRALORDERING": int(io), "db_version": str(ver)}
    if str(ver) != PINNED["DB_VERSION"]:
        rec["note"] = (f"kira.db VERSION {ver!r} != pinned "
                       f"{PINNED['DB_VERSION']!r} (WEIGHTBITS/ordering "
                       "match; recorded, not fatal)")
    return rec


# ---------------------------------------------------- F_p coeff evaluation
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
# m2 path (for ordering-5 banks symbolic in (d, m2)): adds
# ONLY the char 'm' — an expression without 'm' can never enter this path,
# so m2-free inputs keep the pre-change code path verbatim (charset guard).
_ALLOWED_M2 = re.compile(r"^[0-9detam+\-*/^() ]+$")
# literal-wrap for the m2 path: a digit run inside a symbol name ("m2")
# must NOT be wrapped — wrap only integer literals not preceded by a name
# character (the old regex would split "m2" into m F(2)).
_INT_LITERAL_M2 = re.compile(r"(?<![a-zA-Z_])(\d+)")


class CoeffEvaluator:
    """Compiled-eval cache for kira coefficient strings at one (p, d0, eta0)
    [+ optional m20 for (d, m2)-symbolic banks; m2-free inputs are untouched
    by the m2 path — see the charset-guard note at _ALLOWED_M2]."""

    def __init__(self, p, d0, eta0, m20=None):
        self.p = p
        self.env = {"d": Fp(d0, p), "eta": Fp(eta0, p),
                    "F": lambda v, _p=p: Fp(v, _p), "__builtins__": {}}
        if m20 is not None:
            self.env["m2"] = Fp(m20, p)
        self.cache = {}

    def __call__(self, expr):
        code = self.cache.get(expr)
        if code is None:
            if "m" in expr:
                # m2 path (charset-guarded: unreachable for m2-free input)
                if "m2" not in self.env:
                    raise ValueError(
                        f"coeff uses m2 but no m2 value was supplied "
                        f"(pass m20=...): {expr!r}")
                if not _ALLOWED_M2.match(expr):
                    raise ValueError(f"unexpected coeff chars: {expr!r}")
                wrapped = _INT_LITERAL_M2.sub(r"F(\1)",
                                              expr.replace("^", "**"))
            else:
                if not _ALLOWED.match(expr):
                    raise ValueError(f"unexpected coeff chars: {expr!r}")
                # Wrap every integer literal as an F_p element BEFORE
                # compiling: pure-numeric coefficients like (-1)/(2)
                # otherwise evaluate with Python float division — a
                # silent exactness bug (coefficients carrying d/eta
                # mask it).
                wrapped = re.sub(r"(\d+)", r"F(\1)", expr.replace("^", "**"))
            code = compile(wrapped, "<coeff>", "eval")
            self.cache[expr] = code
        val = eval(code, self.env)  # noqa: S307 (restricted charset + no builtins)
        assert isinstance(val, Fp), f"non-Fp coefficient value for {expr!r}"
        return val.v


# ------------------------------------------------------------- file parsing
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


def load_system(art_dir, family, p, d0, eta0, drop_trivial=True,
                expect_counters=None, m20=None):
    """Parse tmp/<family>/SYSTEM_*.gz at one numeric slice -> sysd dict
    (rows/order/stratum_of/sector_of/class_of/masters/counters).

    Reads the artifact tree and evaluates one numeric slice, plus:
      - the version guard record in counters["version_guard"] (raises
        KiraVersionError on a mismatching kira.db);
      - the optional provenance gate expect_counters={"eqs":..,"terms":..}
        (raises ProvenanceError on mismatch);
      - optional m20 for (d, m2)-symbolic banks; omitted =
        pre-change behavior, and an m2-bearing coefficient then raises."""
    guard = version_guard(art_dir)
    ev = CoeffEvaluator(p, d0, eta0, m20)
    trivial = load_trivial_sectors(art_dir, family) if drop_trivial else set()
    files = sorted(glob.glob(os.path.join(art_dir, "tmp", family,
                                          "SYSTEM_*.gz")))
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
    # weights. Empirically (measured on two reference Generate trees):
    # kira grants small ordinals to the PREFERRED masters, in preferred-file
    # order (ordinal rank k <-> k-th preferred entry; sector asserted);
    # every other masters-file entry keeps an ordinary weight. Any wrong
    # assignment is caught by the downstream cross-check (MANDATORY).
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
    # here: a true master can never be pivoted on, so they emerge as
    # elimination SURVIVORS and are labeled post-solve by (sector, class).

    counters["version_guard"] = guard
    if expect_counters:
        for k, v in expect_counters.items():
            if counters.get(k) != v:
                raise ProvenanceError(
                    f"{family}: loaded {k}={counters.get(k)} != declared "
                    f"{v} — artifact tree is not the staged box")

    order = {w: w for w in sector_of if w not in masters}
    stratum_of = {w: popcount(s) for w, s in sector_of.items()}
    return {"rows": rows, "order": order, "stratum_of": stratum_of,
            "sector_of": sector_of, "class_of": class_of, "masters": masters,
            "file_masters": file_masters, "counters": counters}


def to_system(sysd, p):
    """sysd (load_system output) -> ibplapper.System, with the adapter
    default stratification stratum_of = popcount(sector) and meta carrying
    the best-effort labels (opaque to the core, receipt-checked downstream).

    Row-index integrity gate (checker lineage): refuses a load where empty
    rows were dropped, so witness lambda indices always match the caller's
    independent re-parse of the same artifact tree."""
    from ..api import System
    assert sysd["counters"]["cancelled_to_empty"] == 0, \
        "row indexing mismatch risk: loader dropped cancelled-to-empty rows"
    return System(sysd["rows"], p, sysd["order"],
                  forbid=set(sysd["masters"]),
                  stratum_of=sysd["stratum_of"],
                  meta={"sector_of": sysd["sector_of"],
                        "class_of": sysd["class_of"],
                        "masters": sysd["masters"],
                        "counters": sysd["counters"]})


def load(art_dir, family, p, d0, eta0, expect_counters=None, m20=None):
    """One-call convenience: artifact tree -> (System, sysd)."""
    sysd = load_system(art_dir, family, p, d0, eta0,
                       expect_counters=expect_counters, m20=m20)
    return to_system(sysd, p), sysd
