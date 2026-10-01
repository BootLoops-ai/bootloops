#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Part of BootLoops 1.0 (repository-root
# LICENSE and NOTICE). Reads the plain-text state files FireFly writes; contains
# no FireFly code (THIRD_PARTY.md B2).
# ffcapital member — 2-var salvage emitter.
r"""harvest_ffsave_2var.py — ffcapital salvage emitter: extract every FULLY-
RECONSTRUCTED 2-var (d, eta) rational function from a dead kira+FireFly
ff_save (+ snapshot saves), assign each to its (target, master) slot, gate
fail-closed, and bank the results. Zero solver compute — pure decode + exact
arithmetic.

extends: harvest_ffsave_eta (the 1-var member beside this file) — the
'2-var emitter' generalization of the 1-var member.
Monomial exponent layout per
ffsave_degree_census.py pervar_from_monomials (the GPL census tool in the
sibling repository kira, bootloops-tools/; g_ni/g_di lines
"e1 e2 num den"; default FireFly variable order => position 1 = d, position
2 = eta). Ground truth: RatReconst::save_state (RatReconst.cpp:2537-2716);
validation.gz written by Reconstructor.hpp:2157ff (line 1 = the black-box
input values mod the write-time prime; line i+2 = the black-box result of
fn i, tag order).

NAME MAPPING (why not a weight decode): ffsave_degree_census's weight layout
does not close on this configuration (sectormapped standalone job; its file
gate fails and inverse-decode names disagree with kira's own kira_target.m
names — measured, 0/9153 pair overlap). The mapping here is
instead DISCOVERED: FireFly tag order gives consecutive target-weight groups;
kira_target.m rows give named slots; fns are matched to slots by modular
fingerprints at >=2 slice d-values, hardened by group->row constraint
propagation and a global master-weight consistency check, and finally
CERTIFIED per fn by exact cross-multiplication against every slice (gate 2).
Any ambiguity that survives is a typed refusal, never a guess (identical
exact duplicate classes are placed arbitrarily ONLY after proving all class
members exactly equal in every slice — placement then provably immaterial).

GATES (all mandatory, fail-closed; any failure = typed quarantine record,
never a silent drop):
  G1  per-fn black-box validation against EACH save's own validation.gz
      (prime auto-detected per save, unanimity over a sample, assert unique).
  G1b cross-save agreement: a fn done in several saves must be the SAME
      rational function in all of them (digest fast path; cross-multiplication
      slow path — never coefficient-vector diff).
  A   slot assignment: fingerprint match + constraints as above; no-match,
      assignment conflict, or master-weight inconsistency = typed refusal.
  G2  cross-multiplication against the banked exact rational-d slices
      (kira2math kira_target.m): N(d_s,eta)*Q(eta) == P(eta)*D(d_s,eta) as
      exact polynomial identities over Q, on every slice; fails==0 and
      passes >= --min-slice-pass required. Slice-local num/den cancellation
      is immaterial by construction (function-level equality).
  KILL --kill-rate (default 0.1%): if contamination-evidence failures
      (G1 + G1b + no-match/conflict/inconsistency + G2 mismatch) / |union
      done| exceeds it, the save is declared CONTAMINATED: no bank is
      emitted, receipt verdict = CONTAMINATED_STOP, exit rc 3.

usage:
  harvest_ffsave_2var.py salvage \
    --primary <ff_save dir> [--extra label=<ff_save dir> ...] \
    --slice label=<kira_target.m>=<d_num>/<d_den> ... (>=2) \
    [--census <ffsave_degree_census json>]   (degree stats for RESIDUAL_TODO) \
    [--residual-target label=<target file> ...] (staged slices to shrink) \
    --bank <outdir> [--family myfam] [--workers 24] \
    [--min-slice-pass 2] [--kill-rate 0.001]

outputs (in --bank): funcs_2var.jsonl.gz (exact N/D + slot + gate provenance),
extracted_<family>.m (kira2math-compatible, complete+partial rows, banked
terms only), index.json, complete_targets.txt, partial_targets.txt,
EXTRACTION_RECEIPT.json, RESIDUAL_TODO.json, QUARANTINE.jsonl,
residual_targets/<label>.target.
"""
import argparse, gzip, hashlib, json, os, sys, time
from collections import Counter
from fractions import Fraction
from multiprocessing import Pool

sys.setrecursionlimit(200000)

# FireFly 63-bit prime table head (ReconstHelper.cpp:24ff)
PRIMES = [
    9223372036854775783, 9223372036854775643, 9223372036854775549,
    9223372036854775507, 9223372036854775433, 9223372036854775421,
    9223372036854775417, 9223372036854775399, 9223372036854775351,
    9223372036854775337, 9223372036854775291, 9223372036854775279,
    9223372036854775259, 9223372036854775181, 9223372036854775159,
    9223372036854775139, 9223372036854775097, 9223372036854775073,
    9223372036854775057, 9223372036854774959, 9223372036854774937,
    9223372036854774917, 9223372036854774893, 9223372036854774797,
    9223372036854774739, 9223372036854774713, 9223372036854774679,
    9223372036854774629, 9223372036854774587, 9223372036854774571,
    9223372036854774559, 9223372036854774511, 9223372036854774509,
    9223372036854774499, 9223372036854774451, 9223372036854774413,
    9223372036854774341, 9223372036854774319, 9223372036854774307,
    9223372036854774277]

KEYWORDS = {"combined_prime", "tag_name", "is_done", "max_deg_num", "max_deg_den",
            "individual_degrees_num", "individual_degrees_den", "need_prime_shift",
            "normalizer_deg", "normalize_to_den", "normalizer_den_num",
            "shifted_max_num_eqn", "shift", "sub_num", "sub_den", "zero_degs_num",
            "zero_degs_den", "g_ni", "g_di", "combined_ni", "combined_di",
            "combined_primes_ni", "combined_primes_di", "interpolations"}

FP_ETAS = [1234577, 7654321, 33550337, 87178291, 479001599]  # fingerprint eta points
FP_PRIME = PRIMES[0]


def n_primes_exact(cp):
    k = 0
    while cp > 1:
        if k >= len(PRIMES) or cp % PRIMES[k]:
            return -1
        cp //= PRIMES[k]
        k += 1
    return k


def parse_state(path):
    sections, cur = {}, None
    with gzip.open(path, "rt") as f:
        first = f.readline().rstrip("\n")
        if first == "ZERO":
            return {"zero": True}
        cur = first
        sections[cur] = []
        for line in f:
            line = line.rstrip("\n")
            if line in KEYWORDS:
                cur = line
                sections[cur] = []
            elif line.strip():
                sections[cur].append(line)
    return sections


def poly2_from_g(lines, fn, side):
    """g_ni/g_di lines 'e1 e2 num den' -> {(e1,e2): Fraction} (2-var)."""
    P = {}
    for ln in lines:
        p = ln.split()
        assert len(p) == 4, (f"fn {fn} {side}: monomial line has {len(p)} fields "
                             f"(want 4 = 2 exponents + num + den): {ln!r}")
        key = (int(p[0]), int(p[1]))
        assert key not in P, f"fn {fn} {side}: duplicate exponent {key}"
        P[key] = Fraction(int(p[2]), int(p[3]))
    return P


def digest_poly2(N, D):
    h = hashlib.sha256()
    for tag, P in (("N", N), ("D", D)):
        h.update(tag.encode())
        for e, c in sorted(P.items()):
            h.update(f"{e[0]},{e[1]}:{c.numerator}/{c.denominator};".encode())
    return h.hexdigest()


# ---------------- per-save extraction worker ----------------
W = {}  # worker context (fork-inherited)


def _extract_one(fn):
    """states/<fn>_*.gz -> record with 2-var exact N/D + G1 verdict."""
    fname = W["fnmap"].get(fn)
    if fname is None:
        return {"fn": fn, "error": "no unique state file (missing or duplicate)"}
    sec = parse_state(os.path.join(W["ffdir"], "states", fname))
    if sec.get("zero"):
        rec = {"fn": fn, "zero": True, "done": 1, "np": None,
               "N": {}, "D": {(0, 0): Fraction(1)}}
    else:
        done = int(sec["is_done"][0])
        if not done:
            return {"fn": fn, "done": 0}
        tag = sec["tag_name"][0]
        if W["tags"] is not None and tag != W["tags"][fn]:
            return {"fn": fn, "error":
                    f"tag_name {tag} != tags[{fn}] {W['tags'][fn]} — index/tag order broken"}
        rec = {"fn": fn, "done": 1, "tag": tag,
               "np": n_primes_exact(int(sec["combined_prime"][0])),
               "N": poly2_from_g(sec.get("g_ni", []), fn, "g_ni"),
               "D": poly2_from_g(sec.get("g_di", []), fn, "g_di")}
    p, (x1, x2), vals = W["prime"], W["pt"], W["vals"]
    ok, err = _validate_mod(rec["N"], rec["D"], x1, x2, p, vals[fn])
    rec["g1"] = ok
    if not ok:
        rec["g1_error"] = err
    rec["digest"] = digest_poly2(rec["N"], rec["D"])
    rec["degd_num"] = max((e[0] for e in rec["N"]), default=0)
    rec["degd_den"] = max((e[0] for e in rec["D"]), default=0)
    rec["dege_num"] = max((e[1] for e in rec["N"]), default=0)
    rec["dege_den"] = max((e[1] for e in rec["D"]), default=0)
    return rec


def _probe_done(fn):
    fname = W["fnmap"].get(fn)
    if fname is None:
        return (fn, False)
    try:
        with gzip.open(os.path.join(W["ffdir"], "states", fname), "rt") as f:
            lines = [f.readline().rstrip("\n") for _ in range(8)]
    except OSError:
        return (fn, False)
    if lines and lines[0] == "ZERO":
        return (fn, True)
    for i, l in enumerate(lines):
        if l == "is_done" and i + 1 < len(lines):
            return (fn, lines[i + 1].strip() == "1")
    return (fn, False)


def eval2_mod(P, x1, x2, p):
    r = 0
    for (e1, e2), c in P.items():
        dm = c.denominator % p
        if dm == 0:
            return None
        t = (c.numerator % p) * pow(dm, p - 2, p) % p
        r = (r + t * pow(x1, e1, p) * pow(x2, e2, p)) % p
    return r


def _validate_mod(N, D, x1, x2, p, want):
    dv = eval2_mod(D, x1, x2, p)
    if dv is None:
        return False, "coefficient denominator divisible by validation prime"
    if dv == 0:
        return False, "denominator vanishes at validation point"
    nv = eval2_mod(N, x1, x2, p)
    if nv is None:
        return False, "coefficient denominator divisible by validation prime"
    if nv * pow(dv, p - 2, p) % p != want % p:
        return False, "blackbox mismatch"
    return True, None


def load_validation(ffdir):
    vals = []
    with gzip.open(os.path.join(ffdir, "validation.gz"), "rt") as f:
        pt = [int(x) for x in f.readline().split()]
        for ln in f:
            ln = ln.strip()
            if ln:
                vals.append(int(ln))
    assert len(pt) == 2, f"{ffdir}: validation point has {len(pt)} coords (want 2 = [d, eta])"
    return pt, vals


def detect_prime(ffdir, fnmap, pt, vals, done_fns):
    """Unanimous unique prime index under which sampled done fns reproduce
    validation.gz. Sample up to 8 nonzero done fns (lowest fn index)."""
    sample = []
    for fn in done_fns:
        fname = fnmap.get(fn)
        if fname is None:
            continue
        sec = parse_state(os.path.join(ffdir, "states", fname))
        if sec.get("zero") or not int(sec["is_done"][0]):
            continue
        N = poly2_from_g(sec.get("g_ni", []), fn, "g_ni")
        D = poly2_from_g(sec.get("g_di", []), fn, "g_di")
        if not N:
            continue
        sample.append((fn, N, D))
        if len(sample) >= 8:
            break
    assert sample, f"{ffdir}: no nonzero done fn available for prime detection"
    cands = []
    for pi, p in enumerate(PRIMES):
        x1, x2 = pt[0] % p, pt[1] % p
        if all(_validate_mod(N, D, x1, x2, p, vals[fn])[0] for fn, N, D in sample):
            cands.append(pi)
    assert len(cands) == 1, (f"{ffdir}: prime autodetect ambiguous/failed: {cands} "
                             f"(sample fns {[s[0] for s in sample]})")
    return cands[0]


# ---------------- kira_target.m slice parser ----------------
def parse_name(s):
    """'<family>[0,2,-1,...]' -> tuple of ints."""
    i, j = s.index("["), s.index("]")
    return tuple(int(x) for x in s[i + 1:j].split(","))


class ExprParser:
    """Recursive-descent parser for kira2math coefficient expressions:
    integers, 'eta', + - * / ^ ( ). Produces a rational function in eta as
    (num, den) dense Fraction coefficient lists."""
    __slots__ = ("s", "i", "n")

    def __init__(self, s):
        self.s = s
        self.i = 0
        self.n = len(s)

    def parse(self):
        r = self.expr()
        assert self.i == self.n, f"trailing junk at {self.i}: {self.s[self.i:self.i+30]!r}"
        return r

    def peek(self):
        # "\0" sentinel at end: never a member of any operator set ("" would
        # be — the empty string is a substring of every string)
        return self.s[self.i] if self.i < self.n else "\0"

    def expr(self):
        if self.peek() == "-":
            self.i += 1
            r = rneg(self.term())
        else:
            if self.peek() == "+":
                self.i += 1
            r = self.term()
        while self.peek() in "+-":
            op = self.s[self.i]
            self.i += 1
            t = self.term()
            r = radd(r, rneg(t) if op == "-" else t)
        return r

    def term(self):
        r = self.factor()
        while self.peek() in "*/":
            op = self.s[self.i]
            self.i += 1
            f = self.factor()
            r = rmul(r, f) if op == "*" else rdiv(r, f)
        return r

    def factor(self):
        b = self.base()
        if self.peek() == "^":
            self.i += 1
            j = self.i
            while j < self.n and self.s[j].isdigit():
                j += 1
            k = int(self.s[self.i:j])
            self.i = j
            b = rpow(b, k)
        return b

    def base(self):
        c = self.peek()
        if c == "(":
            self.i += 1
            r = self.expr()
            assert self.peek() == ")", f"missing ) at {self.i}"
            self.i += 1
            return r
        if c == "e":
            assert self.s[self.i:self.i + 3] == "eta", f"bad token at {self.i}"
            self.i += 3
            return ([Fraction(0), Fraction(1)], [Fraction(1)])
        assert c.isdigit(), f"bad char {c!r} at {self.i} in {self.s[max(0,self.i-20):self.i+20]!r}"
        j = self.i
        while j < self.n and self.s[j].isdigit():
            j += 1
        v = int(self.s[self.i:j])
        self.i = j
        return ([Fraction(v)], [Fraction(1)])


def ptrim(a):
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def pmul(a, b):
    if (len(a) == 1 and a[0] == 0) or (len(b) == 1 and b[0] == 0):
        return [Fraction(0)]
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    out[i + j] += x * y
    return ptrim(out)


def padd(a, b):
    out = list(a) + [Fraction(0)] * (len(b) - len(a)) if len(b) > len(a) else list(a)
    for j, y in enumerate(b):
        out[j] += y
    return ptrim(out)


def radd(a, b):
    return (padd(pmul(a[0], b[1]), pmul(b[0], a[1])), pmul(a[1], b[1]))


def rneg(a):
    return ([-c for c in a[0]], a[1])


def rmul(a, b):
    return (pmul(a[0], b[0]), pmul(a[1], b[1]))


def rdiv(a, b):
    assert not (len(b[0]) == 1 and b[0][0] == 0), "division by zero expression"
    return (pmul(a[0], b[1]), pmul(a[1], b[0]))


def rpow(a, k):
    r = ([Fraction(1)], [Fraction(1)])
    b = a
    while k:
        if k & 1:
            r = rmul(r, b)
        b = rmul(b, b)
        k >>= 1
    return r


def _parse_term_line(args):
    idx, line = args
    fam = W["family"]
    body = line.strip()
    assert body.startswith("+ "), f"term line {idx}: {body[:40]!r}"
    body = body[2:]
    assert body.startswith(fam + "["), f"term line {idx}: family mismatch {body[:40]!r}"
    close = body.index("]")
    master = parse_name(body[:close + 1])
    rest = body[close + 1:].replace(" ", "")
    assert rest.startswith("*(") and rest.endswith(")"), f"term line {idx}: shape {rest[:30]!r}"
    num, den = ExprParser(rest[2:-1]).parse()
    return idx, master, num, den


def parse_slice_m(path, family, workers):
    """kira_target.m -> ({target: {master: (num, den)}}, [target row order])."""
    rows = {}
    order = []
    cur = None
    term_lines = []
    with open(path) as f:
        for raw in f:
            line = raw.rstrip("\n")
            if not line or line in ("{", "}", ","):
                continue
            if line.startswith(" + "):
                assert cur is not None, "term before any row head"
                term_lines.append((cur, line))
            else:
                assert line.startswith(family + "[") and line.rstrip().endswith("->"), \
                    f"unexpected row head {line[:60]!r}"
                cur = parse_name(line)
                assert cur not in rows, f"duplicate target row {cur}"
                rows[cur] = {}
                order.append(cur)
    global W
    W = {"family": family}
    with Pool(workers) as pool:
        parsed = pool.map(_parse_term_line,
                          [(i, l) for i, (_, l) in enumerate(term_lines)], chunksize=8)
    for (i, master, num, den), (tgt, _) in zip(sorted(parsed), term_lines):
        assert master not in rows[tgt], f"duplicate term {tgt} <- {master}"
        rows[tgt][master] = (num, den)
    return rows, order


# ---------------- fingerprints ----------------
def eval1_mod(coeffs, x, p):
    r = 0
    xp = 1
    for c in coeffs:
        dm = c.denominator % p
        if dm == 0:
            return None
        r = (r + (c.numerator % p) * pow(dm, p - 2, p) % p * xp) % p
        xp = xp * x % p
    return r


def fp_ratfun(num, den, p, etas):
    out = []
    for x in etas:
        nv = eval1_mod(num, x, p)
        dv = eval1_mod(den, x, p)
        if nv is None or dv is None or dv == 0:
            out.append(None)
        else:
            out.append(nv * pow(dv, p - 2, p) % p)
    return tuple(out)


def subst_d(P2, dfrac):
    """{(e1,e2): c} at d=dfrac -> dense eta coefficient list."""
    if not P2:
        return [Fraction(0)]
    dege = max(e[1] for e in P2)
    maxd = max(e[0] for e in P2)
    pw = [Fraction(1)]
    for _ in range(maxd):
        pw.append(pw[-1] * dfrac)
    out = [Fraction(0)] * (dege + 1)
    for (e1, e2), c in P2.items():
        out[e2] += c * pw[e1]
    return ptrim(out)


def _slot_fp_one(args):
    tgt, m = args
    fps = []
    for label, dfrac, rows in W["fpslices"]:
        num, den = rows[tgt][m]
        fps.append(fp_ratfun(num, den, FP_PRIME, FP_ETAS))
    return (tgt, m), tuple(fps)


def _fn_fp_one(args):
    fn, N2, D2 = args
    fps = []
    for label, dfrac, rows in W["fpslices"]:
        fps.append(fp_ratfun(subst_d(N2, dfrac), subst_d(D2, dfrac),
                             FP_PRIME, FP_ETAS))
    return fn, tuple(fps)


# ---------------- gate 2 ----------------
def _gate2_one(args):
    fn, N2, D2, tgt, m = args
    passes, fails, skips = [], [], []
    for label, dfrac, rows in W["slices"]:
        pq = rows.get(tgt, {}).get(m)
        if pq is None:
            skips.append((label, "term absent from slice"))
            continue
        P, Q = pq
        N1 = subst_d(N2, dfrac)
        D1 = subst_d(D2, dfrac)
        if len(D1) == 1 and D1[0] == 0:
            skips.append((label, "denominator vanishes identically at slice d"))
            continue
        if pmul(N1, Q) == pmul(P, D1):
            passes.append(label)
        else:
            fails.append(label)
    return fn, passes, fails, skips


def cross_eq_1var(a, b):
    """(num,den) == (num,den) as rational functions (cross-multiplication)."""
    return pmul(a[0], b[1]) == pmul(b[0], a[1])


def cross_eq_2var(Na, Da, Nb, Db):
    """Na/Da == Nb/Db as rational functions (sparse 2-var cross-multiplication)."""
    def smul(A, B):
        out = {}
        for ea, ca in A.items():
            for eb, cb in B.items():
                k = (ea[0] + eb[0], ea[1] + eb[1])
                out[k] = out.get(k, Fraction(0)) + ca * cb
        return {k: v for k, v in out.items() if v}
    return smul(Na, Db) == smul(Nb, Da)


# ---------------- slot assignment (typed, fail-closed) ----------------
def assign_slots(tags, chosen, slice_data, order, quarantine, notes):
    """Discover fn -> (target, master). Returns dict fn -> slot.
    Typed refusals appended to quarantine (contamination-evidence) or notes
    (identifiability limits, not contamination)."""
    # tag groups (consecutive target-weight runs)
    groups = []
    for i, t in enumerate(tags):
        w1 = t.split("_", 1)[0]
        if not groups or groups[-1][0] != w1:
            groups.append([w1, []])
        groups[-1][1].append(i)
    assert len({g[0] for g in groups}) == len(groups), \
        "target weights repeat non-consecutively in tags — grouping invalid"
    rows0 = slice_data[0][2]
    assert len(groups) == len(order), \
        f"{len(groups)} tag groups != {len(order)} .m rows"
    gs = Counter(len(g[1]) for g in groups)
    rs = Counter(len(rows0[t]) for t in order)
    assert gs == rs, f"group-size multiset != row-size multiset: {gs} vs {rs}"

    # fingerprints (first two slices; 5 eta points each)
    global W
    W = {"fpslices": slice_data[:2]}
    slot_list = [(tgt, m) for tgt in order for m in rows0[tgt]]
    with Pool(W_POOL) as pool:
        slot_fps = dict(pool.map(_slot_fp_one, slot_list, chunksize=64))
        fn_fps = dict(pool.map(
            _fn_fp_one, [(fn, chosen[fn]["N"], chosen[fn]["D"]) for fn in sorted(chosen)],
            chunksize=32))
    fp2slots = {}
    for s, fp in slot_fps.items():
        fp2slots.setdefault(fp, []).append(s)
    fp2fns = {}
    for fn, fp in fn_fps.items():
        fp2fns.setdefault(fp, []).append(fn)

    fn_group = {}
    for gi, (w1, fns) in enumerate(groups):
        for fn in fns:
            fn_group[fn] = gi
    row_of_slot = {s: s[0] for s in slot_list}

    assign = {}
    assign_src = {}   # fn -> "unique" | "class" | "pass2"; only non-"class"
    # assignments carry master-weight information (a "class" placement pairs
    # exactly-identical functions arbitrarily — correct as functions, but its
    # fn<->master-name pairing is arbitrary and must never teach w2 names)
    # candidate rows per group, by size
    size_rows = {}
    for t in order:
        size_rows.setdefault(len(rows0[t]), []).append(t)
    cand_rows = {gi: set(size_rows[len(fns)]) for gi, (w1, fns) in enumerate(groups)}
    # prune by done-fp multiset containment
    row_fp_count = {t: Counter(slot_fps[(t, m)] for m in rows0[t]) for t in order}
    for gi, (w1, fns) in enumerate(groups):
        need = Counter(fn_fps[fn] for fn in fns if fn in fn_fps)
        cand_rows[gi] = {t for t in cand_rows[gi]
                        if all(row_fp_count[t][fp] >= k for fp, k in need.items())}
        if not cand_rows[gi]:
            for fn in fns:
                if fn in fn_fps:
                    quarantine.append({"fn": fn, "gate": "A", "error":
                                       "no candidate row for tag group (fp multiset)"})
    # iterate: unique-candidate fixing + direct fn matches fixing rows
    fixed = {}
    changed = True
    while changed:
        changed = False
        taken = set(fixed.values())
        for gi in list(cand_rows):
            if gi in fixed:
                continue
            c = cand_rows[gi] - taken
            if len(c) == 1:
                r = next(iter(c))
                fixed[gi] = r
                taken.add(r)
                changed = True
        # globally-unique fn<->slot matches pin their group's row
        for fn, fp in fn_fps.items():
            gi = fn_group[fn]
            if gi in fixed:
                continue
            slots = fp2slots.get(fp, [])
            if len(fp2fns[fp]) == 1 and len(slots) == 1:
                newc = cand_rows[gi] & {row_of_slot[slots[0]]}
                if newc != cand_rows[gi]:
                    cand_rows[gi] = newc
                    changed = True
    unresolved_groups = [gi for gi in cand_rows if gi not in fixed]

    # within fixed (group,row): assign fns to slots
    for gi, t in sorted(fixed.items()):
        fns = [fn for fn in groups[gi][1] if fn in fn_fps]
        free = {m: slot_fps[(t, m)] for m in rows0[t]}
        # unique-fp direct
        byfp = {}
        for m, fp in free.items():
            byfp.setdefault(fp, []).append(m)
        fnbyfp = {}
        for fn in fns:
            fnbyfp.setdefault(fn_fps[fn], []).append(fn)
        for fp, fnlist in fnbyfp.items():
            slots = byfp.get(fp, [])
            if not slots:
                for fn in fnlist:
                    quarantine.append({"fn": fn, "gate": "A", "error":
                                       "no slot with matching fp in assigned row"})
                continue
            if len(fnlist) == 1 and len(slots) == 1:
                assign[fnlist[0]] = (t, slots[0])
                assign_src[fnlist[0]] = "unique"
                continue
            # degenerate class: prove all slots exactly identical in EVERY
            # slice and all fns exactly identical; then placement immaterial
            ok = True
            for label, dfrac, rows in slice_data:
                ref = rows[t][slots[0]]
                if not all(cross_eq_1var(rows[t][m], ref) for m in slots[1:]):
                    ok = False
                    break
            d0 = chosen[fnlist[0]]["digest"]
            if ok and not all(chosen[fn]["digest"] == d0 for fn in fnlist[1:]):
                ok = False
            if not ok:
                for fn in fnlist:
                    quarantine.append({"fn": fn, "gate": "A", "error":
                                       "degenerate fp class not exactly identical"})
                continue
            if len(fnlist) > len(slots):
                for fn in fnlist:
                    quarantine.append({"fn": fn, "gate": "A", "error":
                                       f"{len(fnlist)} identical fns > {len(slots)} slots"})
                continue
            # place k fns into first k of j identical slots (immaterial)
            for fn, m in zip(fnlist, slots):
                assign[fn] = (t, m)
                assign_src[fn] = "class"
            if len(fnlist) < len(slots):
                notes.append({"row": list(t), "class_slots": [list(m) for m in slots],
                              "placed": len(fnlist),
                              "note": "identical-coefficient class partially done; "
                                      "uncovered members stay residual"})
    # ---- pass 2: learned master-weight names pin slots in unresolved groups.
    # w2 knowledge comes ONLY from "unique"/"pass2" assignments (never "class").
    # A group is fixed when exactly ONE candidate row admits a placement and
    # that placement is forced (every done fn has exactly one compatible slot).
    changed = True
    while changed:
        changed = False
        w2learn = {}
        for fn, (t, m) in assign.items():
            if assign_src.get(fn) != "class":
                w2learn[tags[fn].split("_", 1)[1]] = m
        taken_rows = set(fixed.values())
        assigned_slots = set(assign.values())
        for gi in [g for g in cand_rows if g not in fixed]:
            donefns = [fn for fn in groups[gi][1] if fn in fn_fps]
            if not donefns:
                continue
            feas = []
            for t in cand_rows[gi] - taken_rows:
                mapping = {}
                used = set()
                ok = True
                for fn in donefns:
                    w2 = tags[fn].split("_", 1)[1]
                    cand = [m for m in rows0[t]
                            if (t, m) not in assigned_slots and m not in used
                            and slot_fps[(t, m)] == fn_fps[fn]
                            and (w2 not in w2learn or w2learn[w2] == m)]
                    if len(cand) != 1:
                        ok = False
                        break
                    mapping[fn] = cand[0]
                    used.add(cand[0])
                if ok:
                    feas.append((t, mapping))
            if len(feas) == 1:
                t, mapping = feas[0]
                fixed[gi] = t
                taken_rows.add(t)
                for fn, m in mapping.items():
                    assign[fn] = (t, m)
                    assign_src[fn] = "pass2"
                    assigned_slots.add((t, m))
                changed = True
    unresolved_groups = [gi for gi in cand_rows if gi not in fixed]
    for gi in unresolved_groups:
        for fn in groups[gi][1]:
            if fn in fn_fps:
                notes.append({"fn": fn, "note": "row unresolved for tag group — "
                              "fn banked without slot (excluded from .m/coverage)"})

    # global master-weight consistency check (strong crosscheck, zero cost):
    # one master weight must map to one master name over all NON-class
    # assignments; conflicts void EVERY fn carrying that weight (fail-closed —
    # cannot tell which side is wrong). "class" placements are exempt: their
    # fn<->name pairing is arbitrary by construction (functions identical).
    w2map = {}
    for fn, (t, m) in assign.items():
        if assign_src.get(fn) != "class":
            w2map.setdefault(tags[fn].split("_", 1)[1], {}).setdefault(m, []).append(fn)
    for w2, byname in w2map.items():
        if len(byname) > 1:
            for m, fns in byname.items():
                for fn in fns:
                    del assign[fn]
                    quarantine.append({"fn": fn, "gate": "A", "error":
                                       f"master weight {w2} maps to "
                                       f"{len(byname)} distinct master names"})
    # target-weight consistency is structural (group->row is a function)
    return assign, assign_src, fixed, unresolved_groups


# ---------------- emission helpers ----------------
def poly2_str(P):
    if not P:
        return "0"
    parts = []
    for (e1, e2), c in sorted(P.items(), key=lambda kv: (-(kv[0][0] + kv[0][1]), kv[0])):
        if c == 0:
            continue
        mono = "*".join(([f"d^{e1}" if e1 > 1 else "d"] if e1 else []) +
                        ([f"eta^{e2}" if e2 > 1 else "eta"] if e2 else []))
        cs = (f"{c.numerator}" if c.denominator == 1 else
              f"{c.numerator}/{c.denominator}")
        if c < 0 or c.denominator != 1:
            cs = f"({cs})"
        parts.append(f"{cs}*{mono}" if mono else cs)
    return "+".join(parts) if parts else "0"


def name_str(family, tup):
    return f"{family}[{','.join(str(x) for x in tup)}]"


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _hist(it):
    h = {}
    for x in it:
        h[str(x)] = h.get(str(x), 0) + 1
    return h


W_POOL = 24


# ---------------- main ----------------
def main():
    global W, W_POOL
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["salvage"])
    ap.add_argument("--primary", required=True, help="main ff_save dir")
    ap.add_argument("--extra", action="append", default=[],
                    help="label=<ff_save dir> snapshot saves (repeatable)")
    ap.add_argument("--census", default=None,
                    help="ffsave_degree_census.py json (degree stats for RESIDUAL_TODO; "
                         "the census tool ships in kira/bootloops-tools)")
    ap.add_argument("--slice", action="append", default=[], dest="slices",
                    help="label=<kira_target.m>=<d_num>/<d_den> (repeatable, >=2)")
    ap.add_argument("--residual-target", action="append", default=[],
                    help="label=<target file> staged slice target lists to residualize")
    ap.add_argument("--bank", required=True)
    ap.add_argument("--family", default="fam",
                    help="integral-family prefix as it appears in kira_target.m names")
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--min-slice-pass", type=int, default=2)
    ap.add_argument("--kill-rate", type=float, default=0.001)
    a = ap.parse_args()
    assert len(a.slices) >= 2, "need >= 2 slices for gate 2"
    W_POOL = a.workers
    t0 = time.time()
    os.makedirs(a.bank, exist_ok=True)
    os.makedirs(os.path.join(a.bank, "residual_targets"), exist_ok=True)
    receipt = {"tool": "harvest_ffsave_2var.py salvage",
               "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
               "extends": "harvest_ffsave_eta (1-var); slot mapping discovered "
                          "(fingerprint + tag-group constraints), not weight-decoded",
               "primary": a.primary, "extra": a.extra, "census": a.census,
               "slices": a.slices, "gates": {},
               "params": {"min_slice_pass": a.min_slice_pass,
                          "kill_rate": a.kill_rate, "workers": a.workers,
                          "fp_etas": FP_ETAS, "fp_prime": FP_PRIME}}

    tags = [ln.strip() for ln in open(os.path.join(a.primary, "tags")) if ln.strip()]
    nfn = len(tags)
    saves = [("primary", a.primary)] + [tuple(x.split("=", 1)) for x in a.extra]

    # ---- per-save: done scan, prime detect, extract + G1
    per_save = {}
    quarantine = []
    notes = []
    for label, ffdir in saves:
        sfiles = os.listdir(os.path.join(ffdir, "states"))
        assert len(sfiles) == nfn, f"{label}: {len(sfiles)} state files != {nfn} tags"
        fnmap, dups = {}, set()
        for fname in sfiles:
            fn = int(fname.split("_")[0])
            if fn in fnmap:
                dups.add(fn)
            fnmap[fn] = fname
        for fn in dups:
            del fnmap[fn]
            quarantine.append({"fn": fn, "save": label, "gate": "G1",
                               "error": "duplicate state files (stale twin — rsync --delete law)"})
        stags = [ln.strip() for ln in open(os.path.join(ffdir, "tags")) if ln.strip()]
        assert stags == tags, f"{label}: tags file differs from primary — fn order broken"
        pt, vals = load_validation(ffdir)
        assert len(vals) == nfn, f"{label}: validation.gz has {len(vals)} values != {nfn}"
        W = {"ffdir": ffdir, "fnmap": fnmap}
        with Pool(a.workers) as pool:
            heads = pool.map(_probe_done, sorted(fnmap), chunksize=128)
        done_fns = sorted(fn for fn, d in heads if d)
        if not done_fns:
            per_save[label] = {"dir": ffdir, "done": 0, "validated": 0,
                               "prime_idx": None, "g1_failures": []}
            print(f"[salvage] {label}: 0 done fns — skipped", file=sys.stderr)
            continue
        pi = detect_prime(ffdir, fnmap, pt, vals, done_fns)
        p = PRIMES[pi]
        W = {"ffdir": ffdir, "fnmap": fnmap, "tags": tags, "prime": p,
             "pt": (pt[0] % p, pt[1] % p), "vals": vals}
        with Pool(a.workers) as pool:
            recs = []
            for k, r in enumerate(pool.imap_unordered(_extract_one, done_fns, chunksize=32)):
                recs.append(r)
                if (k + 1) % 2000 == 0:
                    print(f"[salvage] {label}: extracted {k + 1}/{len(done_fns)}",
                          file=sys.stderr)
        g1fail = []
        for r in recs:
            if r.get("error"):
                g1fail.append({"fn": r["fn"], "save": label, "gate": "G1",
                               "error": r["error"]})
            elif r.get("done") and not r.get("g1"):
                g1fail.append({"fn": r["fn"], "save": label, "gate": "G1",
                               "error": r.get("g1_error", "unknown")})
        quarantine.extend(g1fail)
        per_save[label] = {"dir": ffdir, "done": len(done_fns),
                           "validated": sum(1 for r in recs if r.get("g1")),
                           "prime_idx": pi, "prime": p,
                           "g1_failures": [q["fn"] for q in g1fail],
                           "recs": {r["fn"]: r for r in recs if r.get("g1")}}
        print(f"[salvage] {label}: done={len(done_fns)} validated="
              f"{per_save[label]['validated']} prime_idx={pi} g1_fail={len(g1fail)}",
              file=sys.stderr)

    # ---- union + G1b cross-save agreement
    union = sorted(set().union(*[set(v.get("recs", {})) for v in per_save.values()]))
    disagreements = []
    chosen = {}
    for fn in union:
        cands = [(lbl, v["recs"][fn]) for lbl, v in per_save.items()
                 if fn in v.get("recs", {})]
        digs = {r["digest"] for _, r in cands}
        agree = True
        if len(digs) > 1:
            base = cands[0][1]
            for lbl, r in cands[1:]:
                if r["digest"] == base["digest"]:
                    continue
                if not cross_eq_2var(base["N"], base["D"], r["N"], r["D"]):
                    agree = False
                    disagreements.append({"fn": fn, "gate": "G1b",
                                          "saves": [c[0] for c in cands],
                                          "error": "inter-save disagreement"})
                    break
        if agree:
            cands.sort(key=lambda c: (-(c[1]["np"] or 0), c[0] != "primary"))
            rec = dict(cands[0][1])
            rec["saves_agreeing"] = [c[0] for c in cands]
            chosen[fn] = rec
    quarantine.extend(disagreements)
    print(f"[salvage] union done={len(union)} agree={len(chosen)} "
          f"disagreements={len(disagreements)}", file=sys.stderr)

    # ---- slices: parse; verify identical row/term structure
    slice_specs = []
    for spec in a.slices:
        label, path, dv = spec.split("=")
        dn, dd = dv.split("/")
        slice_specs.append((label, path, Fraction(int(dn), int(dd))))
    slice_data = []
    order = None
    for label, path, dfrac in slice_specs:
        rows, o = parse_slice_m(path, a.family, a.workers)
        if order is None:
            order = o
        else:
            assert o == order, f"slice {label}: row order differs"
            assert all(set(rows[t]) == set(slice_data[0][2][t]) for t in order), \
                f"slice {label}: term support differs"
        slice_data.append((label, dfrac, rows))
        print(f"[salvage] slice {label}: {len(rows)} targets, "
              f"{sum(len(m) for m in rows.values())} terms, d={dfrac}", file=sys.stderr)
    nslots = sum(len(slice_data[0][2][t]) for t in order)
    assert nslots == nfn, f"{nslots} slice slots != {nfn} tags"

    # ---- assignment
    assign, assign_src, fixed_rows, unresolved = assign_slots(
        tags, chosen, slice_data, order, quarantine, notes)
    print(f"[salvage] assignment: {len(assign)} fns -> slots "
          f"(src {_hist(assign_src.get(fn) for fn in assign)}); "
          f"{len(fixed_rows)}/{len(order)} rows fixed; "
          f"{len(unresolved)} groups unresolved", file=sys.stderr)

    # ---- gate 2 (full cross-multiplication on every slice)
    W = {"slices": slice_data}
    jobs = [(fn, chosen[fn]["N"], chosen[fn]["D"], assign[fn][0], assign[fn][1])
            for fn in sorted(assign)]
    with Pool(a.workers) as pool:
        g2 = []
        for k, r in enumerate(pool.imap_unordered(_gate2_one, jobs, chunksize=16)):
            g2.append(r)
            if (k + 1) % 1000 == 0:
                print(f"[salvage] gate2: {k + 1}/{len(jobs)}", file=sys.stderr)
    banked = {}
    g2_mismatch, g2_insufficient = [], []
    for fn, passes, fails, skips in g2:
        rec = chosen[fn]
        rec["g2_pass"], rec["g2_fail"], rec["g2_skip"] = passes, fails, skips
        if fails:
            g2_mismatch.append({"fn": fn, "gate": "G2", "fails": fails,
                                "passes": passes, "skips": skips,
                                "slot": [list(assign[fn][0]), list(assign[fn][1])],
                                "error": "slice cross-multiplication mismatch"})
        elif len(passes) < a.min_slice_pass:
            g2_insufficient.append({"fn": fn, "gate": "G2", "passes": passes,
                                    "skips": skips,
                                    "error": f"insufficient slice evidence (<{a.min_slice_pass})"})
        else:
            banked[fn] = rec
    quarantine.extend(g2_mismatch)
    quarantine.extend(g2_insufficient)
    print(f"[salvage] gate2: banked={len(banked)} mismatch={len(g2_mismatch)} "
          f"insufficient={len(g2_insufficient)}", file=sys.stderr)

    # ---- kill criterion (contamination-evidence failures only)
    hard_failures = (sum(len(v["g1_failures"]) for v in per_save.values())
                     + len(disagreements) + len(g2_mismatch)
                     + sum(1 for q in quarantine if q.get("gate") == "A"))
    rate = hard_failures / max(1, len(union))
    receipt["gates"] = {
        "per_save": {k: {kk: vv for kk, vv in v.items() if kk != "recs"}
                     for k, v in per_save.items()},
        "union_done": len(union),
        "g1b_disagreements": len(disagreements),
        "assignment": {"assigned": len(assign), "rows_fixed": len(fixed_rows),
                       "groups_unresolved": len(unresolved),
                       "src_histogram": _hist(assign_src.get(fn) for fn in assign),
                       "notes": len(notes)},
        "g2": {"banked": len(banked), "mismatch": len(g2_mismatch),
               "insufficient": len(g2_insufficient),
               "pass_histogram": _hist(len(r["g2_pass"]) for r in banked.values())},
        "hard_failures": hard_failures, "failure_rate": rate,
        "kill_rate": a.kill_rate}
    qpath = os.path.join(a.bank, "QUARANTINE.jsonl")
    with open(qpath, "w") as f:
        for q in quarantine:
            f.write(json.dumps(q) + "\n")
        for n in notes:
            f.write(json.dumps({"note_only": True, **n}) + "\n")
    if rate > a.kill_rate:
        receipt["verdict"] = "CONTAMINATED_STOP"
        json.dump(receipt, open(os.path.join(a.bank, "EXTRACTION_RECEIPT.json"), "w"),
                  indent=1)
        print(f"[salvage] KILL: failure rate {rate:.4%} > {a.kill_rate:.2%} — "
              f"save declared CONTAMINATED; no bank emitted", file=sys.stderr)
        sys.exit(3)

    # ---- emit bank
    fpath = os.path.join(a.bank, "funcs_2var.jsonl.gz")
    with gzip.open(fpath, "wt") as f:
        for fn in sorted(banked):
            r = banked[fn]
            tgt, m = assign[fn]
            f.write(json.dumps({
                "fn": fn, "tag": tags[fn], "target": list(tgt), "master": list(m),
                "np": r["np"], "zero": bool(r.get("zero")),
                "degd_num": r["degd_num"], "degd_den": r["degd_den"],
                "dege_num": r["dege_num"], "dege_den": r["dege_den"],
                "N": {f"{e[0]},{e[1]}": f"{c.numerator}/{c.denominator}"
                      for e, c in sorted(r["N"].items())},
                "D": {f"{e[0]},{e[1]}": f"{c.numerator}/{c.denominator}"
                      for e, c in sorted(r["D"].items())},
                "gates": {"g1_saves": r["saves_agreeing"],
                          "g2_slices_pass": r["g2_pass"],
                          "g2_slices_skip": r["g2_skip"]}}) + "\n")

    by_target = {}
    for fn in sorted(banked):
        tgt, m = assign[fn]
        by_target.setdefault(tgt, []).append((m, banked[fn]))
    rows0 = slice_data[0][2]
    complete = sorted(t for t in by_target
                      if {m for m, _ in by_target[t]} == set(rows0[t]))
    partial = sorted(t for t in by_target if t not in set(complete))
    mpath = os.path.join(a.bank, f"extracted_{a.family}.m")
    with open(mpath, "w") as f:
        f.write("{\n")
        first = True
        for tgt in order:
            if tgt not in by_target:
                continue
            if not first:
                f.write(",\n")
            first = False
            f.write(f"{name_str(a.family, tgt)} -> \n")
            for m, r in by_target[tgt]:
                f.write(f" + {name_str(a.family, m)}*(({poly2_str(r['N'])})/"
                        f"({poly2_str(r['D'])}))\n")
        f.write("}\n")
    with open(os.path.join(a.bank, "complete_targets.txt"), "w") as f:
        for t in complete:
            f.write(name_str(a.family, t) + "\n")
    with open(os.path.join(a.bank, "partial_targets.txt"), "w") as f:
        for t in partial:
            f.write(name_str(a.family, t) + "\n")

    json.dump({"family": a.family, "n_banked": len(banked),
               "complete_targets": len(complete), "partial_targets": len(partial),
               "fns": {str(fn): {"tag": tags[fn],
                                 "target": list(assign[fn][0]),
                                 "master": list(assign[fn][1]),
                                 "np": banked[fn]["np"],
                                 "dege_num": banked[fn]["dege_num"],
                                 "dege_den": banked[fn]["dege_den"],
                                 "degd_num": banked[fn]["degd_num"],
                                 "degd_den": banked[fn]["degd_den"]}
                       for fn in sorted(banked)}},
              open(os.path.join(a.bank, "index.json"), "w"), indent=1)

    # RESIDUAL_TODO: residual slots (names) + not-banked fns (census classes)
    covered = {assign[fn] for fn in banked}
    residual_slots = [[list(t), list(m)] for t in order for m in rows0[t]
                      if (t, m) not in covered]
    tag2meta = {}
    if a.census:
        cen = json.load(open(a.census))
        tag2meta = {e["tag"]: e for e in cen["entries"]}
    qfns = {q["fn"] for q in quarantine if "fn" in q}
    residual_fns = []
    for fn in range(nfn):
        if fn in banked:
            continue
        e = tag2meta.get(tags[fn], {})
        nprec = None
        for lbl in per_save:
            r = per_save[lbl].get("recs", {}).get(fn)
            if r is not None:
                nprec = r["np"]
        cls = ("quarantined" if fn in qfns else
               "done_unplaced" if fn in chosen else "not_done")
        residual_fns.append({"fn": fn, "tag": tags[fn], "class": cls,
                             "deg_num": e.get("deg_num"), "deg_den": e.get("deg_den"),
                             "ind_num": e.get("ind_num"), "ind_den": e.get("ind_den"),
                             "np": nprec})
    classes = {}
    for r in residual_fns:
        b = max(r["deg_num"] or 0, r["deg_den"] or 0)
        key = f"{r['class']}|deg<={20 * (b // 20 + 1)}"
        classes[key] = classes.get(key, 0) + 1
    json.dump({"n_residual_fns": len(residual_fns),
               "n_residual_slots": len(residual_slots),
               "classes": classes, "fns": residual_fns,
               "slots": residual_slots},
              open(os.path.join(a.bank, "RESIDUAL_TODO.json"), "w"), indent=1)

    # residualized target lists (drop targets whose EVERY slot is banked)
    complete_set = set(complete)
    resid_written = {}
    for spec in a.residual_target:
        label, path = spec.split("=", 1)
        kept, dropped = [], 0
        block = []

        def flush(block, kept, dropped):
            if not block:
                return dropped
            nm = parse_name(block[0])
            if nm in complete_set:
                return dropped + 1
            kept.extend(block + ["\n"])
            return dropped

        for ln in open(path):
            if ln.strip():
                block.append(ln)
            else:
                dropped = flush(block, kept, dropped)
                block = []
        dropped = flush(block, kept, dropped)
        rp = os.path.join(a.bank, "residual_targets", f"{label}.target")
        with open(rp, "w") as f:
            f.writelines(kept)
        resid_written[label] = {"path": rp, "dropped_targets": dropped,
                                "kept_targets": sum(1 for l in kept if l.strip())}
    receipt["residual_targets"] = resid_written

    receipt["bank"] = {"n_banked": len(banked), "n_union_done": len(union),
                       "n_states": nfn, "n_residual_fns": len(residual_fns),
                       "n_residual_slots": len(residual_slots),
                       "complete_targets": len(complete),
                       "partial_targets": len(partial),
                       "files": {os.path.basename(p): sha256_file(p)
                                 for p in [fpath, mpath, qpath,
                                           os.path.join(a.bank, "index.json"),
                                           os.path.join(a.bank, "RESIDUAL_TODO.json")]}}
    receipt["verdict"] = "BANKED"
    receipt["wall_s"] = round(time.time() - t0, 1)
    json.dump(receipt, open(os.path.join(a.bank, "EXTRACTION_RECEIPT.json"), "w"),
              indent=1)
    print(f"[salvage] BANKED {len(banked)}/{nfn} (union done {len(union)}) in "
          f"{receipt['wall_s']}s; receipt + quarantine written", file=sys.stderr)


if __name__ == "__main__":
    main()
