"""trust.oracle_k2disp — lp_syz check-3 adapter for the lbl3m2L_k2disp
fixture family (2-loop).

This is TRUST package code (NOT vendored bytes): it loads a PRIVATE pinned
instance of the vendored lp_syz engine (the certifier variant with
reduce_target + kira-table compare) and configures the FAMILY surface only:

  - PROPS: lbl3m2L_k2disp — 2 loops (l,k2), 9 propagators; masses: props 1-6
    m2=1, prop 7 m2=eta+5 (the dispersion line), props 8-9 massless ISPs
    (suite/T2_standard gen config, verbatim);
  - GRAM: p1p2=-1/2, p1p3=2/3, p2p3=-1/6, p_i^2=0 (numeric, = the
    kinematic point this family was derived at);
  - build_G: L-loop generalization (LxL Symanzik A-matrix; the vendored
    engines hardcode 3x3);
  - conv: 2-loop Lee-Pomeransky Gamma prefactor — prod_{j=1..w}((L+1)d/2 - j)
    with L=2, i.e. 3d/2, replacing the 3-loop 2d (the vendored conv).

Everything downstream of build_G/conv — Singular syz, tower descent, the
MANDATORY Aut(G) quotient, per-node fmpq RREF, violation rows, target
reduction — runs on the PINNED engine bytes unchanged.
"""
import ast
import itertools
import re
import time
from fractions import Fraction

import sympy as sp

from . import _core

FAMILY = "lbl3m2L_k2disp"

# prop id -> (loop coeffs over (l,k2), ext coeffs over (p1,p2,p3), mass kind)
# mass kind: 0 = massless, 1 = m2=1, "eta5" = m2 = eta+5
PROPS = {
    1: ((1, 0), (1, 0, 0), 1),    # (l+p1)^2 - 1
    2: ((1, -1), (1, 0, 0), 1),   # (l-k2+p1)^2 - 1
    3: ((1, -1), (1, 1, 0), 1),   # (l-k2+p1+p2)^2 - 1
    4: ((1, -1), (1, 1, 1), 1),   # (l-k2+p1+p2+p3)^2 - 1
    5: ((1, 0), (1, 1, 1), 1),    # (l+p1+p2+p3)^2 - 1
    6: ((1, 0), (0, 0, 0), 1),    # l^2 - 1
    7: ((0, 1), (0, 0, 0), "eta5"),  # k2^2 - eta - 5
    8: ((0, 1), (1, 0, 0), 0),    # (k2+p1)^2
    9: ((0, 1), (0, 1, 0), 0),    # (k2+p2)^2
}
NPROPS = 9
NLOOPS = 2

GRAM = [[sp.Integer(0), sp.Rational(-1, 2), sp.Rational(2, 3)],
        [sp.Rational(-1, 2), sp.Integer(0), sp.Rational(-1, 6)],
        [sp.Rational(2, 3), sp.Rational(-1, 6), sp.Integer(0)]]


def _mass(kind, eta0):
    if kind == "eta5":
        return eta0 + 5
    return sp.Integer(kind)


def make_build_G(eta0):
    """L-loop Lee-Pomeransky G = U + F builder for this family."""
    def build_G(node, _eta0_ignored=None):
        n = len(node)
        ys = sp.symbols(f"y1:{n + 1}")
        L = NLOOPS
        A = sp.zeros(L, L)
        Bp = [[sp.Integer(0)] * 3 for _ in range(L)]
        C = sp.Integer(0)
        for k, pid in enumerate(node):
            cL, cP, mm = PROPS[pid]
            y = ys[k]
            for a in range(L):
                for b in range(L):
                    A[a, b] += y * cL[a] * cL[b]
                for e in range(3):
                    Bp[a][e] += y * cL[a] * cP[e]
            pp = sum(cP[e] * GRAM[e][f] * cP[f]
                     for e in range(3) for f in range(3))
            C += y * (pp - _mass(mm, eta0))
        U = sp.expand(A.det())
        if U == 0:
            return None, True
        adj = A.adjugate()
        F = sp.Integer(0)
        for a in range(L):
            for b in range(L):
                BB = sum(Bp[a][e] * GRAM[e][f] * Bp[b][f]
                         for e in range(3) for f in range(3))
                F += adj[a, b] * BB
        F = sp.expand(F - U * C)
        G = sp.expand(U + F)
        P = sp.Poly(G, *ys, domain="QQ")
        Gd = {tuple(e): Fraction(int(c.p), int(c.q)) for e, c in P.terms()}
        return Gd, False
    return build_G


def load_engine(eta0_frac, tag="k2disp"):
    """Private pinned lp_syz instance configured for this family.
    eta0_frac: Fraction. Returns the configured module."""
    m = _core.load_vendored_instance("lp_syz", tag)
    eta0 = sp.Rational(eta0_frac.numerator, eta0_frac.denominator)
    m.build_G = make_build_G(eta0)

    base_tower = m.Tower

    class Tower2(base_tower):
        """2-loop conv: prod_{j=1..w}((L+1)d/2 - j), L=2 -> 3d/2."""

        def conv(self, col, sigma):
            R, mm = self.cols[col] if isinstance(col, int) else col
            w = sum(mm) + len(R)
            c = Fraction(1)
            for mi in mm:
                for k in range(2, mi + 1):
                    c *= k
            half3d = Fraction(3, 2) * self.d0
            for j in range(1, w + 1):
                c /= (half3d - j)
            if sigma == -1 and w % 2:
                c = -c
            return c

    m.Tower = Tower2
    return m


# ------------------------------------------------------------- kira table IO
# The table is SEMI-TRUSTED input (check 3's whole premise is "independent
# oracle when kira itself is in question") and sp.sympify is eval-class — a
# crafted coefficient would execute arbitrary code inside the certifier
# BEFORE any value-level check. Coefficients are rational
# functions in (d, eta): a locked ast grammar (int literals, d, eta,
# + - * / unary +- and integer-literal powers, NOTHING else) is built
# directly into sympy objects. Anything outside the grammar refuses TYPED
# (ValueError) WITHOUT evaluation; sympify is never called on table text.
_D, _ETA = sp.symbols("d eta")
_SYMS = {"d": _D, "eta": _ETA}

# IN-grammar resource bombs: a Pow with a huge integer-literal exponent is
# grammatical but evaluates eagerly at tree build (a 60-byte coefficient
# like '2^10^18' would chew CPU/RAM unbounded), and a literal-zero
# denominator would build sympy zoo silently, dying later UNTYPED in
# eval_coeff. Grammar-time budgets, all refusing typed ValueError:
#   - integer literals longer than 10^4 digits refuse before ast.parse;
#   - |exponent| > 10^4 refuses (_lit_int);
#   - numeric-base powers whose result estimate exceeds 2*10^6 bits refuse
#     (closes nested-pow amplification: each stage's actual bit length feeds
#     the next stage's estimate);
#   - identically-zero denominators (and 0**negative) refuse at build —
#     zoo/nan are never constructed.
# Single-Pow budgets are NOT enough — Mult/Div chains can rebuild refused
# values from in-cap pieces ('2^9999*...*2^9999' x2000 = the same integer
# (2^9999)^2000 refuses as a Pow), same-base symbolic powers auto-MERGE
# their exponents past the literal cap ('(d+1)^9999*...' x2000 -> exponent
# 19,998,000, detonating at eval), sympy EAGERLY distributes an integer
# power over a Mul's numeric head ('(2^9999*d)^9999' computes a 1e8-bit
# integer at construction), and long Add chains build quadratically
# (measured DoS band: 4-44 KB tables -> 6-60+ s, multi-GiB RSS). So:
# AGGREGATE budgets across the WHOLE
# coefficient — ast node count, cumulative numeric bits, cumulative op-units
# (each op charges 1 + the flattened arg count), ANY intermediate's numeric
# content (Number results AND Mul numeric heads, not just Pow), a pre-
# construction estimate for powers of numeric-headed bases, and a post-build
# walk refusing ACCUMULATED (merged) exponents past _MAX_POW_EXP. Evaluation
# rides its own Fraction tree-walk under the SAME budgets (see eval_coeff) —
# typed refusal, never a hang. Caps sized >=10x the pinned real table's
# measured maxima (4405 nodes / 5652 op-units / 1.9e5 cumulative eval bits).
_MAX_POW_EXP = 10 ** 4
_MAX_LIT_DIGITS = 10 ** 4
_MAX_NUM_BITS = 2 * 10 ** 6      # GUARD NUM-BITS (battery mutation anchor)
_MAX_TOTAL_BITS = 2 * 10 ** 6    # GUARD TOTAL-BITS (battery mutation anchor)
_MAX_TOTAL_OPS = 10 ** 5         # GUARD TOTAL-OPS (battery mutation anchor)
_MAX_AST_NODES = 5 * 10 ** 4     # GUARD AST-NODES (battery mutation anchor)
_HUGE_LIT = re.compile(r"[0-9]{%d}" % (_MAX_LIT_DIGITS + 1))
# The raw TABLE byte alphabet — ASCII-printable + newline only
# (the kira .m format is ASCII); any other byte refuses BEFORE parsing.
_TABLE_BAD_BYTE = re.compile(rb"[^\x20-\x7e\r\n]")


def _bits_of(q):
    """Bit size (numerator + denominator) of a sympy Rational/Integer."""
    return int(q.p).bit_length() + int(q.q).bit_length()


def _charge(r, budget):
    """Charge ONE built intermediate against the whole-coefficient
    aggregate budget (budget = [total_bits, total_ops]). Numeric
    content is charged for Number results AND for the Rational head sympy
    keeps as a Mul's first arg (the '2^9999*d'-chain class); op cost adds the
    flattened arg count (the quadratic Add-flood class). Typed refusal."""
    b = 0
    if r.is_Number:
        if r.is_Rational:
            b = _bits_of(r)
    elif r.is_Mul and r.args and r.args[0].is_Rational:
        b = _bits_of(r.args[0])
    if b > _MAX_NUM_BITS:
        raise ValueError(
            f"kira-table coefficient intermediate carries {b} numeric bits — "
            f"exceeds the {_MAX_NUM_BITS}-bit cap (ANY op, not just Pow) — "
            f"refusing (aggregate resource guard)")
    budget[0] += b
    budget[1] += 1 + (len(r.args) if (r.is_Add or r.is_Mul) else 0)
    if budget[0] > _MAX_TOTAL_BITS:
        raise ValueError(
            f"kira-table coefficient CUMULATIVE numeric size {budget[0]} bits "
            f"exceeds the {_MAX_TOTAL_BITS}-bit whole-coefficient budget — "
            f"refusing (rebuild chains inside the per-Pow caps "
            f"must refuse, not evaluate)")
    if budget[1] > _MAX_TOTAL_OPS:
        raise ValueError(
            f"kira-table coefficient build cost {budget[1]} op-units exceeds "
            f"the {_MAX_TOTAL_OPS} whole-coefficient budget — refusing "
            f"(quadratic Add-flood class)")


def _lit_int(nd):
    """Integer-literal (optionally unary-signed) exponent, or typed refusal."""
    if isinstance(nd, ast.UnaryOp) and isinstance(nd.op, (ast.UAdd, ast.USub)):
        v = _lit_int(nd.operand)
        return -v if isinstance(nd.op, ast.USub) else v
    if isinstance(nd, ast.Constant) and type(nd.value) is int:
        if abs(nd.value) > _MAX_POW_EXP:
            raise ValueError(
                f"kira-table coefficient exponent {nd.value} exceeds the "
                f"grammar cap {_MAX_POW_EXP} — refusing (resource-bomb "
                f"guard; in-grammar bombs must refuse, not evaluate)")
        return nd.value
    raise ValueError(
        "kira-table coefficient has a non-integer-literal exponent — "
        "refusing (rational-function grammar)")


def _safe_rational_expr(s):
    """ONE coefficient body -> sympy expr under the locked rational-function
    grammar. Iterative post-order build (no eval/sympify anywhere; long sum
    chains must not hit the recursion limit). Literal/exponent/numeric-power
    budgets refuse typed BEFORE any expensive evaluation, and AGGREGATE
    budgets run across the whole coefficient — ast node
    count, per-intermediate + cumulative numeric bits (any op, incl. Mul
    numeric heads), cumulative op-units, and a post-build walk over merged
    Pow exponents — all typed refusals (see the budget block above)."""
    if _HUGE_LIT.search(s):
        raise ValueError(
            f"kira-table coefficient contains an integer literal longer than "
            f"{_MAX_LIT_DIGITS} digits — refusing (resource-bomb "
            f"guard): {s[:60]!r}")
    try:
        tree = ast.parse(s, mode="eval")
    except (SyntaxError, ValueError, MemoryError, RecursionError) as e:
        raise ValueError(
            f"kira-table coefficient is not parseable arithmetic "
            f"({e.__class__.__name__}) — refusing: {s[:60]!r}")
    n_nodes = sum(1 for _ in ast.walk(tree))
    if n_nodes > _MAX_AST_NODES:
        raise ValueError(
            f"kira-table coefficient has {n_nodes} ast nodes — exceeds the "
            f"{_MAX_AST_NODES}-node whole-coefficient budget — refusing "
            f"(aggregate resource guard; the pinned real table "
            f"maxes at ~4.4e3 nodes)")
    budget = [0, 0]   # cumulative [bits, op-units], whole build
    val = {}
    stack = [(tree.body, False)]
    while stack:
        nd, done = stack.pop()
        if done:
            if isinstance(nd, ast.UnaryOp):
                v = val.pop(id(nd.operand))
                val[id(nd)] = -v if isinstance(nd.op, ast.USub) else v
            elif isinstance(nd.op, ast.Pow):
                base = val.pop(id(nd.left))
                exp = _lit_int(nd.right)
                head = None
                if base.is_Number:
                    if not base.is_Rational:
                        raise ValueError(
                            "kira-table coefficient has a non-rational "
                            "numeric power base — refusing typed")
                    if base.is_zero and exp < 0:
                        raise ValueError(
                            "kira-table coefficient raises zero to a "
                            "negative power — refusing typed at build "
                            "(sympy zoo is never constructed)")
                    head = base
                elif base.is_Mul and base.args and base.args[0].is_Rational:
                    # sympy distributes an integer power over a Mul's
                    # numeric head EAGERLY — (2^9999*d)^9999 computes
                    # a 1e8-bit integer at construction; estimate FIRST.
                    head = base.args[0]
                if head is not None:
                    bits = _bits_of(head)
                    if bits * abs(exp) > _MAX_NUM_BITS:
                        raise ValueError(
                            f"kira-table coefficient numeric power result "
                            f"estimate {bits}*{abs(exp)} bits exceeds the "
                            f"{_MAX_NUM_BITS}-bit budget — refusing "
                            f"(resource-bomb guard; Mul-head "
                            f"bases are estimated too)")
                val[id(nd)] = base ** exp
                _charge(val[id(nd)], budget)
            else:
                a = val.pop(id(nd.left))
                b = val.pop(id(nd.right))
                if isinstance(nd.op, ast.Add):
                    val[id(nd)] = a + b
                elif isinstance(nd.op, ast.Sub):
                    val[id(nd)] = a - b
                elif isinstance(nd.op, ast.Mult):
                    val[id(nd)] = a * b
                else:
                    if b.is_zero:   # GUARD DIV0 (battery mutation anchor)
                        raise ValueError(
                            "kira-table coefficient divides by an "
                            "identically zero denominator — refusing typed "
                            "at build (sympy zoo is never "
                            "constructed)")
                    val[id(nd)] = a / b
                _charge(val[id(nd)], budget)
            continue
        if isinstance(nd, ast.Constant) and type(nd.value) is int:
            val[id(nd)] = sp.Integer(nd.value)
        elif (isinstance(nd, ast.Name) and isinstance(nd.ctx, ast.Load)
              and nd.id in _SYMS):
            val[id(nd)] = _SYMS[nd.id]
        elif isinstance(nd, ast.UnaryOp) and isinstance(nd.op,
                                                        (ast.UAdd, ast.USub)):
            stack.append((nd, True))
            stack.append((nd.operand, False))
        elif isinstance(nd, ast.BinOp) and isinstance(
                nd.op, (ast.Add, ast.Sub, ast.Mult, ast.Div)):
            stack.append((nd, True))
            stack.append((nd.right, False))
            stack.append((nd.left, False))
        elif isinstance(nd, ast.BinOp) and isinstance(nd.op, ast.Pow):
            _lit_int(nd.right)          # grammar-gate the exponent up front
            stack.append((nd, True))
            stack.append((nd.left, False))
        else:
            raise ValueError(
                f"kira-table coefficient contains a non-rational-function "
                f"construct ({type(nd).__name__}) — REFUSING WITHOUT "
                f"EVALUATION (semi-trusted input; coefficients must be "
                f"rational functions in (d, eta))")
    ex = val[id(tree.body)]
    # Same-base symbolic powers auto-MERGE their exponents past
    # the literal cap ((d+1)^9999 * ... * (d+1)^9999 -> exponent 19,998,000,
    # detonating later at eval) — merged exponents route through the same
    # cap the literals face.
    for _pw in ex.atoms(sp.Pow):   # GUARD EXP-WALK (battery mutation anchor)
        _pe = _pw.exp
        if _pe.is_Integer and abs(int(_pe)) > _MAX_POW_EXP:
            raise ValueError(
                f"kira-table coefficient carries an ACCUMULATED power "
                f"exponent {int(_pe)} past the {_MAX_POW_EXP} cap (same-base "
                f"merges route through the same estimate as literals) — "
                f"refusing (aggregate resource guard)")
    return ex


def parse_kira_table(path):
    """Parse ALL entries of a kira_target.m-style table for this family.
    Returns dict target_nu -> dict master_nu -> sympy coeff expr.
    Coefficient bodies go through the LOCKED grammar (_safe_rational_expr) —
    never sympify/eval; out-of-grammar bodies refuse typed.
    PARSE COVERAGE is asserted — every raw '->' record must
    become a parsed entry, and every non-blank line of every entry body must
    be consumed by exactly one term of the exact ' + fam[...]*(...)' form
    (duplicate masters/targets refuse too). Any residue refuses typed,
    NAMING the offending line numbers: the 'refuses WITHOUT evaluation'
    contract binds the whole byte stream, never just the regex-matching
    subset — a table is never certified against a silently parsed SUBSET of
    its bytes.
    BYTE ALPHABET asserted on the RAW bytes BEFORE parsing —
    the table path accepts ASCII-printable + newline (LF/CR) only; any other
    byte refuses typed naming offset and line. Regex-based coverage counters
    alone are blind to a lookalike record (U+2192 arrow, fullwidth digit,
    lookalike family name) — it matches NEITHER counter and would be
    silently dropped; the alphabet gate kills every lookalike class at once
    and makes the whole-byte-stream contract hold byte-for-byte."""
    raw = open(path, "rb").read()
    _mba = _TABLE_BAD_BYTE.search(raw)   # GUARD ALPHABET (battery mutation anchor)
    if _mba is not None:
        _off = _mba.start()
        _ln = raw.count(b"\n", 0, _off) + 1
        raise ValueError(
            f"kira-table byte-alphabet violation: byte 0x{raw[_off]:02x} at "
            f"offset {_off} (line {_ln}) of {path} — the table path accepts "
            f"ASCII-printable + newline bytes only; refusing typed BEFORE "
            f"parsing (lookalike bytes match no grammar and "
            f"are invisible to regex-based coverage counters)")
    txt = (raw.decode("ascii", "surrogateescape")
           .replace("\r\n", "\n").replace("\r", "\n"))
    out = {}
    fam = re.escape(FAMILY)
    ent_pat = re.compile(fam + r"\[([0-9,\s\-]+)\] -> \n(.*?)(?:\n,\n|\n\}\Z|\n\}\n)",
                         re.S)
    term_line_pat = re.compile(
        r"\s*\+\s*" + fam + r"\[([0-9,\s\-]+)\]\*(\(.*\))\s*")
    rec_pat = re.compile(fam + r"\[[0-9,\s\-]*\]\s*->")
    ents = list(ent_pat.finditer(txt))
    raw_recs = list(rec_pat.finditer(txt))
    if len(ents) != len(raw_recs):   # GUARD COVERAGE (battery anchor A)
        matched = {m.start() for m in ents}
        lines = [txt.count("\n", 0, m.start()) + 1
                 for m in raw_recs if m.start() not in matched]
        raise ValueError(
            f"kira-table parse coverage FAILED: {len(raw_recs)} raw '->' "
            f"records but {len(ents)} parsed entries — unconsumed record(s) "
            f"at line(s) {lines[:8]} of {path} — refusing typed "
            f"(out-of-grammar content refuses, never silently dropped)")
    for mm in ents:
        nu = tuple(int(x) for x in mm.group(1).split(","))
        if nu in out:
            raise ValueError(
                f"kira-table duplicate entry target {nu} at line "
                f"{txt.count(chr(10), 0, mm.start()) + 1} of {path} — "
                f"refusing typed (silent overwrite forbidden)")
        body = mm.group(2)
        body_line0 = txt.count("\n", 0, mm.start(2)) + 1
        terms = {}
        bad_lines = []
        for i, ln in enumerate(body.split("\n")):
            if not ln.strip():
                continue
            m2 = term_line_pat.fullmatch(ln)
            nu2 = (tuple(int(x) for x in m2.group(1).split(","))
                   if m2 else None)
            if m2 is None or nu2 in terms:
                bad_lines.append(body_line0 + i)
                continue
            terms[nu2] = _safe_rational_expr(m2.group(2).replace("^", "**"))
        if bad_lines:                # GUARD COVERAGE (battery anchor B)
            raise ValueError(
                f"kira-table entry {nu}: {len(bad_lines)} body line(s) not "
                f"consumed by the term grammar (or duplicate master) at "
                f"line(s) {bad_lines[:8]} of {path} — refusing typed "
                f"(parse coverage: certifying a parsed SUBSET "
                f"of the table bytes is forbidden)")
        out[nu] = terms
    return out


def _charge_frac(r, budget):
    """Eval-side twin of _charge: one Fraction intermediate
    against the whole-coefficient budget. Typed refusal."""
    b = r.numerator.bit_length() + r.denominator.bit_length()
    if b > _MAX_NUM_BITS:
        raise ValueError(
            f"kira-table coefficient evaluation intermediate carries {b} "
            f"bits — exceeds the {_MAX_NUM_BITS}-bit cap — refusing "
            f"(eval soft budget)")
    budget[0] += b
    budget[1] += 1
    if budget[0] > _MAX_TOTAL_BITS:
        raise ValueError(
            f"kira-table coefficient CUMULATIVE evaluation size {budget[0]} "
            f"bits exceeds the {_MAX_TOTAL_BITS}-bit whole-coefficient "
            f"budget — refusing (eval soft budget)")
    if budget[1] > _MAX_TOTAL_OPS:
        raise ValueError(
            f"kira-table coefficient evaluation cost {budget[1]} ops exceeds "
            f"the {_MAX_TOTAL_OPS} whole-coefficient budget — refusing "
            f"(eval soft budget)")


def _eval_fraction(ex, dv, ev, budget):
    """Evaluator: stdlib-Fraction tree-walk over the LOCKED
    grammar's node set (Integer/Rational/Symbol d,eta/Add/Mul/Pow) with the
    aggregate budgets checked at every step — an in-grammar bomb refuses
    typed BEFORE the expensive operation runs; sympy .subs (which evaluates
    first and can hang unbounded) is off the table path. Any node outside
    the grammar set refuses typed (grammar-built exprs never contain one).
    A vanishing denominator raises ZeroDivisionError for the caller to type."""
    if ex is _D:
        return dv
    if ex is _ETA:
        return ev
    if ex.is_Integer:
        return Fraction(int(ex))
    if ex.is_Rational:
        return Fraction(int(ex.p), int(ex.q))
    if ex.is_Add or ex.is_Mul:
        r = _eval_fraction(ex.args[0], dv, ev, budget)
        for a in ex.args[1:]:
            v = _eval_fraction(a, dv, ev, budget)
            r = (r + v) if ex.is_Add else (r * v)
            _charge_frac(r, budget)
        return r
    if ex.is_Pow:
        pe = ex.exp
        if not pe.is_Integer:
            raise ValueError(
                "kira-table coefficient has a non-integer power exponent at "
                "evaluation — refusing (locked grammar)")
        e = int(pe)
        if abs(e) > _MAX_POW_EXP:
            raise ValueError(
                f"kira-table coefficient carries an ACCUMULATED power "
                f"exponent {e} past the {_MAX_POW_EXP} cap at evaluation — "
                f"refusing (eval soft budget; never a hang)")
        b = _eval_fraction(ex.base, dv, ev, budget)
        bits = b.numerator.bit_length() + b.denominator.bit_length()
        if bits * max(1, abs(e)) > _MAX_NUM_BITS:
            raise ValueError(
                f"kira-table coefficient power result estimate "
                f"{bits}*{abs(e)} bits exceeds the {_MAX_NUM_BITS}-bit "
                f"budget at evaluation — refusing (eval estimate guard)")
        if b == 0 and e < 0:
            raise ZeroDivisionError("pole at evaluation point")
        r = b ** e
        _charge_frac(r, budget)
        return r
    raise ValueError(
        f"kira-table coefficient contains a non-grammar node "
        f"({type(ex).__name__}) at evaluation — refusing (locked grammar)")


def eval_coeff(ex, d0, eta0):
    """Exact-QQ evaluation of ONE table coefficient at rational (d0, eta0).
    Accepts the sympy exprs parse_kira_table built, or a raw string — a
    string goes through the SAME locked grammar (never sympify).
    Evaluation is a stdlib-Fraction tree-walk under an INTERNAL
    soft budget (the same aggregate caps as the build: per-intermediate and
    cumulative bits, op count, merged-exponent cap) — every over-budget path
    refuses typed ValueError BEFORE the expensive operation, never a hang;
    sympy .subs is retired from the table path. A denominator
    vanishing AT the point refuses typed."""
    if isinstance(ex, (str, bytes)):
        ex = _safe_rational_expr(
            (ex.decode() if isinstance(ex, bytes) else ex).replace("^", "**"))
    dv = Fraction(d0.numerator, d0.denominator)
    ev = Fraction(eta0.numerator, eta0.denominator)
    try:
        return _eval_fraction(ex, dv, ev, [0, 0])
    except ZeroDivisionError:
        # a denominator vanishing AT THE POINT must refuse typed here — it
        # must never die later as an untyped TypeError ('invalid input: zoo').
        raise ValueError(
            f"kira-table coefficient is non-finite at (d0,eta0)=({d0},{eta0})"
            f" — a denominator vanishes at this evaluation point; refusing "
            f"typed: pick another rational point")
    except RecursionError:
        raise ValueError(
            "kira-table coefficient expr too deeply nested to evaluate — "
            "refusing typed (eval soft budget)")


def support(nu):
    return tuple(i + 1 for i, v in enumerate(nu) if v > 0)


# ------------------------------------------------------------------ driver
def reduce_and_compare(table_path, d0, eta0, targets=None, D0=None,
                       verbose=True, engine_tag="k2disp"):
    """Build ONE tower over the union support of the chosen targets (+ their
    kira masters), reduce every target on the pinned engine, compare exact-QQ
    vs the kira table at rational (d0, eta0). Returns result dict."""
    t_all = time.time()
    ents = parse_kira_table(table_path)
    # fail CLOSED on an empty/junk table and on requested-but-absent
    # targets — a silent filter would leave a success-shaped dict (sigma=1,
    # n_violations=0, n_match=0) indistinguishable from a clean certification.
    if not ents:
        raise ValueError(
            f"kira table {table_path!r} parsed to 0 {FAMILY} entries — not a "
            f"table for this family; refusing a vacuous certification")
    if targets is None:
        targets = sorted(ents)
    else:
        targets = [tuple(t) for t in targets]
        missing = [t for t in targets if t not in ents]
        if missing:
            raise KeyError(
                f"{len(missing)} requested target(s) absent from the kira "
                f"table (e.g. {missing[0]}) — refusing a vacuous certification")
    if not targets:
        raise ValueError(
            "empty target list — refusing a vacuous certification")
    union = set()
    maxdot = 0
    masters = set()
    for t in targets:
        union.update(support(t))
        maxdot = max(maxdot, sum(v - 1 for v in t if v > 0))
        for k in ents[t]:
            union.update(support(k))
            masters.add(k)
            maxdot = max(maxdot, sum(v - 1 for v in k if v > 0))
    top = tuple(sorted(union))
    if D0 is None:
        D0 = maxdot + 1
    eng = load_engine(eta0, tag=engine_tag)
    res = {"family": FAMILY, "top": list(top), "n_targets": len(targets),
           "n_forced_masters": len(masters), "D0": D0,
           "d0": str(d0), "eta0": str(eta0)}
    tw = eng.Tower(top, sp.Rational(eta0.numerator, eta0.denominator),
                   d0, D0, sorted(masters), verbose=verbose)
    tw.build_relations()
    tw.solve()
    res.update({k: tw.stats[k] for k in
                ("nodes_nonzero", "n_repclasses", "n_columns", "n_relations",
                 "n_dropped", "t_singular", "t_build", "t_solve",
                 "n_violations", "n_free_nonmaster")})

    details_by_sigma = {}
    chosen = None
    for sigma in (1, -1):
        allok = True
        details = []
        for nu in targets:
            mine = tw.reduce_target(nu, sigma)
            mine2 = {k: v for k, v in mine.items() if v != 0}
            kmap = {}
            for knu, kex in ents[nu].items():
                K = support(knu)
                mK = tuple(knu[i - 1] - 1 for i in K)
                kc = tw.canon(K, mK) if K in tw.rep else None
                if kc is None:
                    GK, _ = eng.build_G(K, None)
                    for R in tw.reps:
                        if len(R) != len(K):
                            continue
                        sig = tw._iso_find(GK, tw.G[R], len(K))
                        if sig is not None:
                            m2 = [0] * len(K)
                            for i, v in enumerate(mK):
                                m2[sig[i]] = v
                            kc = (R, tw._aut_min(R, tuple(m2)))
                            break
                if kc is None:
                    details.append({"target": list(nu),
                                    "match": False,
                                    "reason": f"kira master {knu} unmapped"})
                    allok = False
                    break
                kmap[kc] = kmap.get(kc, Fraction(0)) + eval_coeff(kex, d0, eta0)
            else:
                diffs = []
                for k in set(mine2) | set(kmap):
                    a = mine2.get(k, Fraction(0))
                    b = kmap.get(k, Fraction(0))
                    if a != b:
                        diffs.append((str(k), str(a), str(b)))
                ok = not diffs
                details.append({"target": list(nu), "match": ok,
                                "n_masters_kira": len(kmap),
                                "n_masters_mine": len(mine2),
                                "diffs": diffs[:6]})
                if not ok:
                    allok = False
        details_by_sigma[sigma] = details
        if allok:
            chosen = sigma
            break
    res["sigma"] = chosen
    res["details"] = details_by_sigma.get(chosen if chosen is not None else -1,
                                          details_by_sigma.get(-1))
    if chosen is None:
        res["details_sigma_1"] = details_by_sigma.get(1)
    res["n_match"] = sum(1 for dd in (res["details"] or []) if dd.get("match"))
    res["wall_s"] = round(time.time() - t_all, 1)
    res["_tower"] = tw          # for downstream legs (witness emission)
    res["_engine"] = eng
    res["_targets"] = targets
    res["_entries"] = ents
    return res


def tower_relation_rows_modp(tw, p):
    """The tower's exact relation rows evaluated mod p (row order = node
    order then insertion order — the order lam indexes). Rows whose
    denominators vanish mod p are REFUSED loudly (caller picks another p)."""
    rows = []
    for R in tw.node_order:
        for r in tw.rel_by_node[R]:
            row = {}
            for c, v in r.items():
                den = v.denominator % p
                if den == 0:
                    raise ZeroDivisionError(
                        f"denominator vanishes mod {p} — pick another prime")
                row[c] = (v.numerator * pow(den, p - 2, p)) % p
            row = {c: v for c, v in row.items() if v}
            if row:
                rows.append(row)
    return rows
