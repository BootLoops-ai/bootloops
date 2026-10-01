"""adapters/user_system.py — kira user_defined_system reader (Route B).

Reads the DOCUMENTED, stable kira input format for user-defined systems of
equations — the same files a `reduce_user_defined_system` job consumes,
i.e. what a USER hands kira, not an internal dump. Because the format is
documented and stable there is nothing empirical to pin: Route B carries
NO version guard and NO weight decode (contrast Route A, adapters/kira.py,
which parses an undocumented internal artifact behind a version guard and
can refuse with KiraVersionError / record ABSENT-UNGUARDED).

Format (per the kira documentation and shipped examples):
  - `src` is one equation file or a directory; a directory contributes its
    non-hidden regular files with extension .kira or .kira.gz, sorted by
    name. Files may be gzip-compressed or plain text either way.
  - An equation is a run of consecutive non-blank lines; a blank (or
    whitespace-only) line ends it. Equations never span files.
  - Each line is ONE term,  <integral>*<coefficient>,  split at the FIRST
    `*`. The integral is either
      * an integer WEIGHT  (weight notation: the integer is the column id
        AND its Laporta rank — higher weight eliminated first), or
      * name[i1,i2,...,in] (integral notation: indices = propagator
        powers; every integral in one load must carry the same index
        count).
    Mixing the two notations in one load is refused loudly.
  - The coefficient is a rational expression in d and the kinematic
    invariants. It is evaluated here at ONE numeric (point, prime) slice —
    rows hand the core plain F_p values, per the library contract. The
    caller declares every symbol with its slice value (`values=`); an
    undeclared symbol in a coefficient is a loud error, never a guess.

Ordering (integral notation): each distinct integral gets an opaque column
id equal to its rank in the ascending sort by

    (t, r, s, sector, family, indices)

with t = number of positive indices, r = sum of positive indices, s = sum
of |negative indices|, sector = bitmask of positive indices. Higher rank =
eliminated first — the standard Laporta direction (more lines, then more
dots, then more numerator first; corner integrals at the bottom).
stratum_of defaults to t, matching the Route A default. In weight notation
the file's integers are used verbatim (column id = rank = weight) and the
default is a single stratum. Either way the order/stratum maps returned
here are only the adapter DEFAULT: order, strata, and masters (`forbid`)
remain the caller's physics surface, exactly as for any System.

Quickstart:

    from ibplapper.adapters import user_system
    system, sysd = user_system.load(
        "eqs/", p, values={"d": 1234577, "s": 87654321, "m2": 424243},
        forbid=[("T", (1, 0, 0))],       # or raw column ids
        expect_counters={"eqs": 7135, "terms": 42286})

`sysd["col_of"]` maps (family, indices) -> column id, `sysd["labels"]`
inverts it; witness lambda indices match this loader's row order, so a
downstream verifier re-parses the same files independently and calls
`ibplapper.receipt.verify` — mandatory, as for every route.
"""
import gzip
import os
import re

from .kira import Fp, ProvenanceError

_EXTS = (".kira", ".kira.gz")
_INTEGRAL = re.compile(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*"
                       r"\[\s*([-0-9,\s]+)\]\s*$")
_WEIGHT = re.compile(r"^\s*(\d+)\s*$")
_IDENT = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")
# wrap integer literals as F_p elements BEFORE compile (pure-numeric
# coefficients otherwise evaluate with float division — the silent
# exactness bug Route A documents); a digit run inside a symbol name
# ("m2") must NOT be wrapped.
_INT_LITERAL = re.compile(r"(?<![A-Za-z0-9_])(\d+)")
_ALLOWED = re.compile(r"^[0-9A-Za-z_+\-*/^() ]+$")


class UserSystemFormatError(ValueError):
    """A line/file does not follow the documented user-system format."""


class SliceEvaluator:
    """Compiled-eval cache for coefficient strings at one numeric slice.

    values: dict symbol -> integer slice value (e.g. {"d": d0, "s": s0}).
    Every identifier appearing in a coefficient must be declared here —
    an unknown symbol raises, never guesses. Same restricted-charset +
    no-builtins + F_p-literal-wrap discipline as the Route A evaluator.
    """

    def __init__(self, p, values):
        self.p = int(p)
        self.env = {"__builtins__": {}, "F": lambda v, _p=p: Fp(v, _p)}
        for name, v in dict(values).items():
            if not _IDENT.fullmatch(name) or name == "F":
                raise ValueError(f"bad symbol name {name!r}")
            self.env[name] = Fp(int(v), p)
        self.cache = {}

    def __call__(self, expr, where=""):
        code = self.cache.get(expr)
        if code is None:
            if not _ALLOWED.match(expr):
                raise UserSystemFormatError(
                    f"{where}: unexpected coefficient chars: {expr!r}")
            unknown = sorted({s for s in _IDENT.findall(expr)
                              if s not in self.env})
            if unknown:
                raise UserSystemFormatError(
                    f"{where}: coefficient uses undeclared symbols "
                    f"{unknown} (declare them in values=): {expr!r}")
            wrapped = _INT_LITERAL.sub(r"F(\1)", expr.replace("^", "**"))
            code = compile(wrapped, "<coeff>", "eval")
            self.cache[expr] = code
        val = eval(code, self.env)  # noqa: S307 (restricted charset, no builtins)
        assert isinstance(val, Fp), f"non-Fp coefficient value for {expr!r}"
        return val.v


def _source_files(src):
    """One file, or a directory's non-hidden *.kira / *.kira.gz files
    sorted by name (the documented directory convention)."""
    if os.path.isfile(src):
        return [src]
    if not os.path.isdir(src):
        raise FileNotFoundError(f"no such file or directory: {src}")
    names = sorted(
        n for n in os.listdir(src)
        if not n.startswith(".") and n.endswith(_EXTS)
        and os.path.isfile(os.path.join(src, n)))
    if not names:
        raise FileNotFoundError(
            f"{src}: no {'/'.join(_EXTS)} files in directory")
    return [os.path.join(src, n) for n in names]


def _open_text(path):
    """gzip or plain text, decided by content (both are legal for either
    extension — kira reads compressed and uncompressed alike)."""
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_terms(src):
    """Parse the equation files into raw term lists (no evaluation).

    Returns (eqs, notation, np, n_files) where eqs is a list of equations,
    each a list of (key, coeff_string) with key = int weight (weight
    notation) or (family, indices tuple) (integral notation); notation is
    "weight" or "integral"; np is the index count (integral notation) or
    None; n_files counts the source files read.
    """
    eqs = []
    notation = None
    np = None
    n_files = 0
    for path in _source_files(src):
        n_files += 1
        cur = []
        with _open_text(path) as fh:
            for ln, raw in enumerate(fh, 1):
                line = raw.rstrip("\r\n")
                if not line.strip():
                    if cur:
                        eqs.append(cur)
                        cur = []
                    continue
                where = f"{path}:{ln}"
                if "*" not in line:
                    raise UserSystemFormatError(
                        f"{where}: ill-formed integral*coefficient: "
                        f"{line!r}")
                lhs, coeff = line.split("*", 1)
                m = _WEIGHT.match(lhs)
                if m:
                    kind, key = "weight", int(m.group(1))
                else:
                    m = _INTEGRAL.match(lhs)
                    if not m:
                        raise UserSystemFormatError(
                            f"{where}: integral is neither an integer "
                            f"weight nor name[indices]: {lhs!r}")
                    fam = m.group(1)
                    idx = tuple(int(x) for x in m.group(2).split(","))
                    kind, key = "integral", (fam, idx)
                    if np is None:
                        np = len(idx)
                    elif len(idx) != np:
                        raise UserSystemFormatError(
                            f"{where}: {len(idx)} indices, expected {np}")
                if notation is None:
                    notation = kind
                elif notation != kind:
                    raise UserSystemFormatError(
                        f"{where}: mixed weight/integral notation in one "
                        "load is refused (ordering semantics differ)")
                cur.append((key, coeff.strip()))
        if cur:                     # equations never span files
            eqs.append(cur)
    if not eqs:
        raise UserSystemFormatError(f"{src}: no equations found")
    return eqs, notation, np, n_files


def _integral_rank_key(fam_idx):
    """Ascending Laporta key: corner integrals first, deeper/dressed last."""
    fam, idx = fam_idx
    t = sum(1 for a in idx if a > 0)
    r = sum(a for a in idx if a > 0)
    s = sum(-a for a in idx if a < 0)
    sector = sum(1 << i for i, a in enumerate(idx) if a > 0)
    return (t, r, s, sector, fam, idx)


def load_system(src, p, values, expect_counters=None):
    """Parse + evaluate at one numeric slice -> sysd dict.

    Returns {"rows", "order", "stratum_of", "notation", "labels",
    "col_of", "counters"}; order covers EVERY column (restrict via
    forbid= at to_system). The optional provenance gate expect_counters=
    {"eqs":..., "terms":...} raises ProvenanceError on mismatch, same
    contract as Route A.
    """
    eqs, notation, np_, n_files = parse_terms(src)
    ev = SliceEvaluator(p, values)

    labels, col_of = None, None
    if notation == "integral":
        distinct = sorted({key for eq in eqs for key, _ in eq},
                          key=_integral_rank_key)
        col_of = {key: rank for rank, key in enumerate(distinct)}
        labels = {rank: key for key, rank in col_of.items()}

    rows = []
    counters = {"eqs": 0, "terms": 0, "zero_coeff": 0,
                "cancelled_to_empty": 0, "files": n_files,
                "notation": notation}
    for eq in eqs:
        row = {}
        for key, coeff in eq:
            counters["terms"] += 1
            cv = ev(coeff, where=f"eq {counters['eqs']}")
            if cv == 0:
                counters["zero_coeff"] += 1
                continue
            c = key if notation == "weight" else col_of[key]
            row[c] = (row.get(c, 0) + cv) % p
            if row[c] == 0:
                del row[c]
        counters["eqs"] += 1
        if row:
            rows.append(row)
        else:
            counters["cancelled_to_empty"] += 1

    if notation == "weight":
        cols = sorted({c for r in rows for c in r})
        order = {c: c for c in cols}
        stratum_of = {c: 0 for c in cols}       # no structure to stratify on
    else:
        order = {c: c for c in labels}
        stratum_of = {c: _integral_rank_key(labels[c])[0] for c in labels}

    if expect_counters:
        for k, v in expect_counters.items():
            if counters.get(k) != v:
                raise ProvenanceError(
                    f"user system {src}: loaded {k}={counters.get(k)} != "
                    f"declared {v} — input is not the declared system")

    return {"rows": rows, "order": order, "stratum_of": stratum_of,
            "notation": notation, "labels": labels, "col_of": col_of,
            "counters": counters}


def resolve_cols(sysd, cols):
    """Map a mixed list of column ids / (family, indices) labels to ids."""
    out = set()
    for c in cols:
        if isinstance(c, int):
            out.add(c)
            continue
        fam, idx = c
        key = (fam, tuple(idx))
        if not sysd["col_of"] or key not in sysd["col_of"]:
            raise KeyError(
                f"unknown integral label {key} (labels resolve only in "
                "integral notation, against integrals present in the load)")
        out.add(sysd["col_of"][key])
    return out


def to_system(sysd, p, forbid=()):
    """sysd -> ibplapper.System. forbid: master columns, as raw ids or
    (family, indices) labels; forbidden columns leave `order` (a master
    never carries an elimination rank).

    Row-index integrity gate (same contract as Route A): refuses a load
    where cancelled-to-empty rows were dropped, so witness lambda indices
    always match an independent re-parse of the same files."""
    from ..api import System
    assert sysd["counters"]["cancelled_to_empty"] == 0, \
        "row indexing mismatch risk: loader dropped cancelled-to-empty rows"
    fb = resolve_cols(sysd, forbid)
    order = {c: r for c, r in sysd["order"].items() if c not in fb}
    return System(sysd["rows"], p, order, forbid=fb,
                  stratum_of=sysd["stratum_of"],
                  meta={"labels": sysd["labels"],
                        "notation": sysd["notation"],
                        "counters": sysd["counters"]})


def load(src, p, values, forbid=(), expect_counters=None):
    """One-call convenience: equation files -> (System, sysd)."""
    sysd = load_system(src, p, values, expect_counters=expect_counters)
    return to_system(sysd, p, forbid=forbid), sysd
