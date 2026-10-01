#!/usr/bin/env python3
"""weightdict.py — weight<->integral dictionary via a kira full-box
kira2math run, matched against the Route-B closed rows.

kira never exports its weight->integral map; we recover it by running kira's
OWN solve (triangular + back substitution + kira2math) on an explicit
mandatory list of EVERY box integral, then matching index-keyed reduction
rows to our weight-keyed closed rows by reduction-vector equality at the
ANCHOR slice, held out on the remaining slices. The kira table doubles as a
full-table cross-check oracle (kira's algebra vs Route-B, per slice).

Matching honesty: identical reduction vectors (symmetric partners) form an
ambiguity CLASS; assignment inside a class is arbitrary-but-fixed and
counted. A fault cannot hide there: all class members share the same row at
every slice by construction of the match.
"""
import os
import re
import subprocess
import sys
import time

from loader import _ALLOWED, Fp

# Kept in sync with corpus_bench.parse_oracle_table (same kira2math grammar).
_ROW_HEAD = re.compile(r"^\s*(\w+)\[([-\d, ]+)\]\s*->\s*(.*)$")
_TERM = re.compile(r"([+-])?\s*(\w+)\[([-\d, ]+)\](?:\*\((.*)\))?\s*,?\s*$")

DICT_JOBS = """jobs:
 - reduce_sectors:
    reduce:
     - {{topologies: [{fam}], sectors: [{secs}], r: {r}, s: {s}}}
    select_integrals:
      select_mandatory_list:
        - [{fam}, target]
    preferred_masters: preferred
    integral_ordering: 5
    run_initiate: true
    run_triangular: true
    run_back_substitution: true
    run_firefly: false
 - kira2math:
    target:
     - [{fam}, target]
"""


def enumerate_box(sectors, n_idx, rmax, smax, dmax=None):
    """All integrals of the given (canonical, nontrivial) sectors inside the
    (rmax, smax) selection box; per sector the dots window is rmax - t.

    dmax=None (default): dots run to the full r-budget
    rmax - t per sector. Measured falsification (x3): the recursive solve selection pivots deep-dot integrals at low
    sectors (dots 6-7 at t=2 — the r-budget allows it, the d-budget does not
    bind selection), so capping dots at solve.d+2 left those
    (sector,(dots,s)) classes ABSENT from the dict kira table (k2disp: 232
    clean rows unmatchable, heldout_mismatch=0). Passing an integer dmax
    reproduces the legacy capped window min(dmax, rmax - t) (regression /
    reproduction use only)."""
    import itertools
    out = []
    for sec in sectors:
        pos = [i for i in range(n_idx) if sec >> i & 1]
        neg = [i for i in range(n_idx) if not sec >> i & 1]
        t = len(pos)
        dot_hi = rmax - t if dmax is None else min(dmax, rmax - t)
        for dots in range(0, dot_hi + 1):
            for dist in itertools.combinations_with_replacement(pos, dots):
                for stot in range(0, smax + 1):
                    for sdist in itertools.combinations_with_replacement(
                            neg, stot):
                        a = [0] * n_idx
                        for i in pos:
                            a[i] = 1
                        for i in dist:
                            a[i] += 1
                        for i in sdist:
                            a[i] -= 1
                        out.append(tuple(a))
    return sorted(set(out))


def run_dict_kira(workdir, fam, staging, box, sectors_seen, n_idx,
                  timeout=1800):
    """Stage and run the dictionary solve; returns (wall_s, n_targets)."""
    import shutil
    os.makedirs(workdir, exist_ok=True)
    shutil.copytree(os.path.join(staging, "config"),
                    os.path.join(workdir, "config"), dirs_exist_ok=True)
    shutil.copy(os.path.join(staging, "preferred"),
                os.path.join(workdir, "preferred"))
    tuples = enumerate_box(sectors_seen, n_idx, box["r"], box["s"],
                           box.get("d"))
    with open(os.path.join(workdir, "target"), "w") as fh:
        for a in tuples:
            fh.write(f"{fam}[{','.join(map(str, a))}]\n")
    with open(os.path.join(workdir, "jobs.yaml"), "w") as fh:
        fh.write(DICT_JOBS.format(fam=fam,
                                  secs=",".join(map(str, box["sectors"])),
                                  r=box["r"], s=box["s"]))
    t0 = time.time()
    with open(os.path.join(workdir, "kira_dict.log"), "w") as lg:
        rc = subprocess.call(["nice", "-n", "10", "kira", "--parallel=4",
                              "jobs.yaml"], cwd=workdir, stdout=lg,
                             stderr=subprocess.STDOUT, timeout=timeout)
    assert rc == 0, f"dictionary kira rc={rc} in {workdir}"
    return round(time.time() - t0, 2), len(tuples)


# ---------------------------------------------------------------------------
# Eta-dict mode of record: the dict table is consumed STREAMING
# (single pass, one compile per distinct expr, no materialization) and the
# dict kira run is NUMERIC-SLICE by default (eta pinned per distinct slice
# eta; measured on a production dict run: total 1493 s vs stream 2107 s vs batch ~2800 s /
# 17.3 GB; all correctness gates green, zero row diffs over 77,347 tuples x
# 4 slices x 2 primes). Stream mode retained for slice-hungry runs
# (>~3 eta points: numeric pays ~one kira run per eta, stream pays one
# symbolic run flat).

def config_has_eta(staging):
    """True iff eta is USED in the config (propagators / scalarproduct
    rules), not merely declared as an invariant. Declared-but-unused eta
    (vac3-class eta-stage configs) means the table coefficients cannot
    depend on eta, so numeric-slice pinning is a no-op and the caller falls
    back to the single symbolic run consumed streaming."""
    decl = re.compile(r"(?m)^\s*- \[eta, *\d+\]\s*$")
    for f in ("kinematics.yaml", "integralfamilies.yaml"):
        txt = decl.sub("", open(os.path.join(staging, "config", f)).read())
        if re.search(r"\beta\b", txt):
            return True
    return False


def stage_numeric_config(staging, workdir, eta_val):
    """Write an eta-pinned scratch COPY of the config into workdir/config
    (never in place). Returns True if eta was pinned, False if the config
    has no eta (caller falls back to the symbolic run). Fails loudly on an
    unrecognized eta declaration (degradation 3a catches it upstream)."""
    kin = open(os.path.join(staging, "config", "kinematics.yaml")).read()
    fam_y = open(os.path.join(staging, "config",
                              "integralfamilies.yaml")).read()
    if not (re.search(r"\beta\b", kin) or re.search(r"\beta\b", fam_y)):
        return False
    old = "  kinematic_invariants:\n    - [eta, 2]"
    if old not in kin:
        raise RuntimeError("eta present but invariant block not in the "
                           "validated form; numeric-slice staging refused")
    kin = kin.replace(old, "  kinematic_invariants: []")
    rep = f"({int(eta_val)})"
    kin, n1 = re.subn(r"\beta\b", rep, kin)
    fam_y, n2 = re.subn(r"\beta\b", rep, fam_y)
    assert n1 + n2 >= 1, "no eta occurrence substituted"
    os.makedirs(os.path.join(workdir, "config"), exist_ok=True)
    with open(os.path.join(workdir, "config", "kinematics.yaml"), "w") as fh:
        fh.write(kin)
    with open(os.path.join(workdir, "config",
                           "integralfamilies.yaml"), "w") as fh:
        fh.write(fam_y)
    return True


def run_dict_kira_numeric(workdir, fam, staging, box, sectors_seen, n_idx,
                          eta_val, timeout=1800, target_tuples=None):
    """run_dict_kira with eta pinned NUMERIC in a scratch config copy
    (numkin pattern; Fermat backend, so numeric kinematics is exact — the
    --set_value FireFly footgun does not apply because we stage the config).
    target_tuples: optional explicit tuple list (--dict-scope
    core/targets; None = enumerate_box over the dict box). Writes
    DICT_ETA.json (eta + n_targets stamp) so numeric-slice tables are
    reusable across legs (--reuse-dict numeric-slice layout).
    Returns (wall_s, n_targets)."""
    import json as _json
    import shutil
    os.makedirs(workdir, exist_ok=True)
    if not stage_numeric_config(staging, workdir, eta_val):
        raise RuntimeError("config has no eta; use run_dict_kira")
    shutil.copy(os.path.join(staging, "preferred"),
                os.path.join(workdir, "preferred"))
    tuples = (sorted(set(map(tuple, target_tuples)))
              if target_tuples is not None else
              enumerate_box(sectors_seen, n_idx, box["r"], box["s"],
                            box.get("d")))
    with open(os.path.join(workdir, "target"), "w") as fh:
        for a in tuples:
            fh.write(f"{fam}[{','.join(map(str, a))}]\n")
    with open(os.path.join(workdir, "DICT_ETA.json"), "w") as fh:
        _json.dump({"eta": eta_val, "n_targets": len(tuples),
                    "scope": "explicit" if target_tuples is not None
                    else "box"}, fh)
    with open(os.path.join(workdir, "jobs.yaml"), "w") as fh:
        fh.write(DICT_JOBS.format(fam=fam,
                                  secs=",".join(map(str, box["sectors"])),
                                  r=box["r"], s=box["s"]))
    t0 = time.time()
    with open(os.path.join(workdir, "kira_dict.log"), "w") as lg:
        rc = subprocess.call(["nice", "-n", "10", "kira", "--parallel=4",
                              "jobs.yaml"], cwd=workdir, stdout=lg,
                             stderr=subprocess.STDOUT, timeout=timeout)
    assert rc == 0, f"numeric dictionary kira rc={rc} in {workdir}"
    return round(time.time() - t0, 2), len(tuples)


def run_dict_kira_numeric_persector(workdir, fam, staging, n_idx, eta_val,
                                    target_tuples, box_r, box_s=0,
                                    timeout_per_leg=1500, r_budget="trim",
                                    on_timeout="raise", sector_subset=None):
    """Per-sector dict legs (measured finding: the eta-pinned mandatory-list dict
    kira is CELL-REDUCE-BOUND, not target-count-bound: 56.9 s @10 tuples /
    >=2160 s @1000 / >=6900 s @631 on the sh sec1018 cell).

    One eta-pinned dict kira run PER TARGET SECTOR instead of one cell run:
    each leg's reduce entry is scoped to that sector (kira closes the
    subsector chain itself) with the r-budget trimmed to the leg's own
    targets (r_budget="trim": r_leg = max target r; "box": legacy cell r).
    Correctness is NOT assumed from rc=0 (sec1015 lesson: incompleteness is
    not self-signaling): callers must gate leg tables against banked
    references / the independent oracle rows.

    Legs run smallest-sector-first (ladder order) so a caller-imposed
    cumulative cap censors the expensive tail, not the cheap head.
    on_timeout="continue": a timed-out/failed leg is recorded
    {"failed": ...} and later legs still run (driver mode); "raise" keeps
    corpus_bench degradation-3a semantics. Timed-out kira children are
    killed strictly cwd-scoped to the leg dir.

    Returns list of leg dicts {"sector","t","r","s","n_targets","wall_s",
    "workdir","table"} (failed legs carry "failed" instead of "table")."""
    import json as _json
    import shutil
    import signal
    os.makedirs(workdir, exist_ok=True)
    by_sec = {}
    for tup in sorted(set(map(tuple, target_tuples))):
        sec = sum(1 << i for i, a in enumerate(tup) if a > 0)
        by_sec.setdefault(sec, []).append(tup)
    if sector_subset is not None:
        by_sec = {s: t for s, t in by_sec.items() if s in set(sector_subset)}
    def _t(sec):
        return bin(sec).count("1")
    legs = []
    for sec in sorted(by_sec, key=lambda s: (_t(s), len(by_sec[s]), s)):
        tups = by_sec[sec]
        r_leg = (max(sum(a for a in tup if a > 0) for tup in tups)
                 if r_budget == "trim" else box_r)
        s_leg = max([box_s] + [-sum(a for a in tup if a < 0)
                               for tup in tups]) if r_budget == "trim" \
            else box_s
        legdir = os.path.join(workdir, f"leg_sec{sec}")
        os.makedirs(legdir, exist_ok=True)
        if not stage_numeric_config(staging, legdir, eta_val):
            raise RuntimeError("config has no eta; persector needs numeric")
        shutil.copy(os.path.join(staging, "preferred"),
                    os.path.join(legdir, "preferred"))
        with open(os.path.join(legdir, "target"), "w") as fh:
            for a in tups:
                fh.write(f"{fam}[{','.join(map(str, a))}]\n")
        with open(os.path.join(legdir, "DICT_ETA.json"), "w") as fh:
            _json.dump({"eta": eta_val, "n_targets": len(tups),
                        "scope": "persector", "sector": sec,
                        "r": r_leg, "s": s_leg}, fh)
        with open(os.path.join(legdir, "jobs.yaml"), "w") as fh:
            fh.write(DICT_JOBS.format(fam=fam, secs=str(sec),
                                      r=r_leg, s=s_leg))
        leg = {"sector": sec, "t": _t(sec), "r": r_leg, "s": s_leg,
               "n_targets": len(tups), "workdir": legdir}
        t0 = time.time()
        try:
            with open(os.path.join(legdir, "kira_dict.log"), "w") as lg:
                rc = subprocess.call(
                    ["nice", "-n", "10", "kira", "--parallel=4", "jobs.yaml"],
                    cwd=legdir, stdout=lg, stderr=subprocess.STDOUT,
                    timeout=timeout_per_leg)
            leg["wall_s"] = round(time.time() - t0, 2)
            if rc != 0:
                raise RuntimeError(f"persector dict kira rc={rc} in {legdir}")
            leg["table"] = kira_table_path(legdir, fam)
        except Exception as exc:  # noqa: BLE001 — timeout/rc/no-table
            leg["wall_s"] = round(time.time() - t0, 2)
            leg["failed"] = f"{type(exc).__name__}: {exc}"[:200]
            # cwd-scoped orphan sweep (fer64 children of a killed kira),
            # strictly under THIS leg dir — never touches foreign jobs.
            # /proc is Linux-only: without it there is no cwd-scoped sweep,
            # so skip LOUDLY rather than crash inside this handler.
            legreal = os.path.realpath(legdir)
            if os.path.isdir("/proc"):
                for pd in os.listdir("/proc"):
                    if not pd.isdigit():
                        continue
                    try:
                        cwd = os.readlink(f"/proc/{pd}/cwd")
                    except OSError:
                        continue
                    if cwd == legreal or cwd.startswith(legreal + os.sep):
                        try:
                            os.kill(int(pd), signal.SIGKILL)
                        except OSError:
                            pass
            else:
                print(f"[weightdict] WARN: no /proc on this platform — "
                      f"skipping the cwd-scoped orphan sweep under {legdir}; "
                      f"stray engine children may survive", file=sys.stderr)
            if on_timeout == "raise":
                legs.append(leg)
                raise RuntimeError(
                    f"persector leg sec{sec} failed: {leg['failed']} "
                    f"(legs so far: {len(legs)})") from exc
        legs.append(leg)
    return legs


class MultiEvalSome:
    """Per-slice-TOLERANT variant of MultiEval for the streamed oracle
    comparator: one compile per distinct expr, but a
    denominator hit at one slice must NOT censor the other slices (the batch
    reference — corpus_bench.eval_vec — reports 'denhit' per slice, first
    failing term wins at that slice only). Bounded cache like MultiEval."""
    CACHE_MAXLEN = 64

    def __init__(self, slices):
        self.envs = [{"d": Fp(d0, p), "eta": Fp(e0, p),
                      "F": (lambda v, _p=p: Fp(v, _p)), "__builtins__": {}}
                     for (p, d0, e0) in slices]
        self.cache = {}

    def eval_some(self, expr, active):
        """-> {si: value or None-on-denhit} for si in active."""
        code = self.cache.get(expr)
        if code is None:
            if not _ALLOWED.match(expr):
                raise ValueError(f"unexpected coeff chars: {expr!r}")
            code = compile(re.sub(r"(\d+)", r"F(\1)",
                                  expr.replace("^", "**")), "<coeff>", "eval")
            if len(expr) <= self.CACHE_MAXLEN:
                self.cache[expr] = code
        out = {}
        for si in active:
            try:
                out[si] = eval(code, self.envs[si]).v  # noqa: S307 restricted
            except ZeroDivisionError:
                out[si] = None
        return out


def stream_rows(path):
    """Yield (tgt_tuple, [(sign, midx, expr), ...]) one row at a time with
    the EXACT line discipline of corpus_bench.parse_oracle_table (gated
    byte-equivalent on the banked lp1disp table: same 77,347 rows, identical
    dict stats)."""
    cur = None
    with open(path) as fh:
        for raw in fh:
            line = raw.strip()
            if line in ("{", "}", ""):
                continue
            if line == ",":
                if cur is not None:
                    yield cur
                cur = None
                continue
            m = _ROW_HEAD.match(line)
            if m:
                if cur is not None:
                    yield cur
                cur = (tuple(int(x) for x in m.group(2).split(",")), [])
                rest = m.group(3).strip().rstrip(",")
                if rest in ("", "+", "0"):
                    continue
                line = rest
            if cur is None:
                raise ValueError(f"orphan line in {path}: {raw!r}")
            body = line.rstrip(",").strip()
            if body in ("0", "+ 0", "+0", "- 0", "-0") or not body:
                continue
            t = _TERM.match(body)
            if not t:
                raise ValueError(f"unparsed term in {path}: {raw!r}")
            cur[1].append((-1 if t.group(1) == "-" else 1,
                           tuple(int(x) for x in t.group(3).split(",")),
                           t.group(4) or "1"))
    if cur is not None:
        yield cur


class MultiEval:
    """One compile per distinct expr across N slice envs. Only short exprs
    are cached (kira2math back-sub coefficients are mostly unique big
    rationals — caching them would rebuild the batch-mode RSS wall)."""
    CACHE_MAXLEN = 64

    def __init__(self, slices):
        self.envs = [{"d": Fp(d0, p), "eta": Fp(e0, p),
                      "F": (lambda v, _p=p: Fp(v, _p)), "__builtins__": {}}
                     for (p, d0, e0) in slices]
        self.cache = {}

    def eval_all(self, expr):
        code = self.cache.get(expr)
        if code is None:
            if not _ALLOWED.match(expr):
                raise ValueError(f"unexpected coeff chars: {expr!r}")
            code = compile(re.sub(r"(\d+)", r"F(\1)",
                                  expr.replace("^", "**")), "<coeff>", "eval")
            if len(expr) <= self.CACHE_MAXLEN:
                self.cache[expr] = code
        out = []
        for env in self.envs:
            v = eval(code, env)   # noqa: S307 restricted charset, no builtins
            out.append(v.v)
        return out


def stream_eval_table(path, slices, midx2w):
    """Single-pass streaming parse + evaluation of a kira2math table at ALL
    given (p, d0, eta0) slices. Returns (kira_rows per slice, n_rows,
    n_skipped). Skip semantics identical to the batch path: a tuple with a
    denominator hit or an unmapped master at ANY slice is dropped at ALL
    slices."""
    me = MultiEval(slices)
    ps = [s[0] for s in slices]
    krs = [{} for _ in slices]
    n_rows = n_skip = 0
    for tgt, terms in stream_rows(path):
        n_rows += 1
        vecs = [{} for _ in slices]
        ok = True
        for sign, midx, expr in terms:
            mw = midx2w.get(midx)
            if mw is None:
                ok = False
                break
            try:
                vals = me.eval_all(expr)
            except ZeroDivisionError:
                ok = False
                break
            for k, v in enumerate(vals):
                vecs[k][mw] = (vecs[k].get(mw, 0) - sign * v) % ps[k]
        if not ok:
            n_skip += 1
            continue
        for k in range(len(slices)):
            krs[k][tgt] = {c: v for c, v in vecs[k].items() if v}
    return krs, n_rows, n_skip


def kira_table_path(workdir, fam):
    for cand in (os.path.join(workdir, "results", fam, "kira_target.m"),):
        if os.path.exists(cand):
            return cand
    hits = []
    for root, _d, files in os.walk(os.path.join(workdir, "results")):
        hits += [os.path.join(root, f) for f in files if f.endswith(".m")]
    assert hits, f"no kira2math output under {workdir}/results"
    return hits[0]


def match_dictionary(kira_rows_by_slice, closed_by_slice, midx2w, p_of_slice):
    """kira_rows_by_slice: per slice {tuple: {master_weight: NEGATED val}}
    (same convention as closed rows). closed_by_slice: per slice
    {weight: row}. Returns (dict tuple->weight, stats)."""
    n_slices = len(closed_by_slice)
    stats = {"n_tuples": 0, "matched": 0, "ambiguous_class": 0,
             "unmatched": 0, "heldout_mismatch": 0, "master_tuples": 0,
             "zero_tuples": 0}
    # index closed rows at anchor slice by canonical signature
    def sig(row):
        return tuple(sorted(row.items()))
    anchor_index = {}
    for w, row in closed_by_slice[0].items():
        anchor_index.setdefault(sig(row), []).append(w)
    used = set()
    wdict = {}
    p0 = p_of_slice[0]
    for a in sorted(kira_rows_by_slice[0]):
        stats["n_tuples"] += 1
        vec0 = kira_rows_by_slice[0][a]
        sec = sum(1 << i for i, x in enumerate(a) if x > 0)
        if not vec0:                    # reduces to zero at anchor
            stats["zero_tuples"] += 1
            continue
        if a in midx2w and vec0 == {midx2w[a]: p0 - 1}:
            wdict[a] = midx2w[a]        # the tuple IS a master
            stats["master_tuples"] += 1
            continue
        cands = [w for w in anchor_index.get(sig(vec0), []) if w not in used]
        if not cands:
            stats["unmatched"] += 1
            continue
        if len(anchor_index.get(sig(vec0), [])) > 1:
            stats["ambiguous_class"] += 1
        w = min(cands)
        ok = all(closed_by_slice[s].get(w, {}) == kira_rows_by_slice[s][a]
                 for s in range(1, n_slices))
        if ok:
            used.add(w)
            wdict[a] = w
            stats["matched"] += 1
        else:
            stats["heldout_mismatch"] += 1
    return wdict, stats
