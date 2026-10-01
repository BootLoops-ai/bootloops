"""baller.hygiene_checks — the ball footgun catalog as ENFORCED checks
(each check cites its source discipline).

Four source LINTS (scan .py text/AST; findings -> nonzero return like
dps_lint) and two RUNTIME helpers:

  lint_fixed_wdps       workdps(<literal>) / mpf('1e-K') guards must scale
                        with dps (dps_lint FOOTGUNS in GUIDE.md; cancellation
                        guards need ~2*dps)
  lint_ambient_conv     mpf()/mpf(str()) conversions outside any workdps
                        block blind two-precision checks (FIXES.md:31 trap;
                        the popgen/mplll ambient-dps law)
  lint_unkeyed_cache    lru_cache on precision-reading functions: an un-keyed
                        cache lets the first low-dps call cap all later
                        inputs
  lint_ctx_restore      flint ctx.prec/ctx.cap are process-global — every
                        mutation must be save/restore-wrapped (acbfast law)

  ftrim(vec)            ball midpoints once per ACCEPTED step — radii never
                        accumulate across steps. LAW: call ONLY after the
                        step was validated/accepted; trimming an unvalidated
                        ball discards rigor.
  RadiusWatch(bar)      halts (RadiusBlowup) when a ball's relative radius
                        log10(rad/|mid|) crosses -bar digits — the alert
                        lever that catches radius blowup DURING a run.
  bridge_mpf(x)         re-round-safe mpf bridge: mp.mpf(x) on an existing
                        mpf re-rounds to AMBIENT dps; the bridge goes through
                        ._mpf_ raw (eval_int_J2L law: 'mantissas are dyadic;
                        no re-rounding').
"""
import ast
import re


# ---------- lints ----------

def _findings_from_ast(path, text, visitor_cls):
    try:
        tree = ast.parse(text)
    except SyntaxError as e:
        return [(path, 0, f"SYNTAX (unlintable): {e}")]
    v = visitor_cls()
    v.visit(tree)
    return [(path, ln, msg) for ln, msg in v.findings]


def lint_fixed_wdps(path, text=None):
    """Flag fixed-precision guards: workdps with an int literal, and mpf
    small-literal cutoffs in quoted (any case/mantissa) or bare-float form.
    NOT FLAGGED (documented reach limits): guards routed through
    named constants; non-mpf comparison literals; call spellings the regex
    cannot see (alias-renamed workdps calls, space before the paren) —
    dataflow and token games are beyond a lint; reviewers own that class.
    (Examples are deliberately NOT written literally here — the lint lints
    its own docstring.)"""
    text = open(path).read() if text is None else text
    out = []
    for i, line in enumerate(text.split("\n"), 1):
        if re.search(r"workdps\(\s*\d+\s*\)", line):
            out.append((path, i, "workdps(<literal>): fixed working precision"
                        " — must scale with the caller's dps"))
        if re.search(r"""mpf\(\s*['"][0-9.]*[eE]-\d+['"]\s*\)""", line):
            out.append((path, i, "mpf('<small literal>') cutoff — a fixed"
                        " guard ceilings held-out digits; scale with dps"))
        if re.search(r"mpf\(\s*[0-9.]*[eE]-\d+\s*\)", line):
            out.append((path, i, "mpf(<float literal>) cutoff — fixed guard"
                        " AND a float-literal rounding trap; scale with dps"))
    return out


class _AmbientConv(ast.NodeVisitor):
    """Recognizes workdps under IMPORT ALIASES (`from mpmath import
    workdps as w`) and any attribute spelling (`mp.workdps`/`ctx.workdps`);
    tracks mpf ALIASES (`M = mp.mpf`) so aliased conversions are flagged too;
    the aliased-with recognition kills the false positive on properly guarded
    aliased blocks."""

    def __init__(self):
        self.findings = []
        self._wdps_depth = 0
        self._wdps_names = {"workdps"}
        self._mpf_names = {"mpf"}

    def visit_ImportFrom(self, node):
        if node.module == "mpmath":
            for a in node.names:
                if a.name == "workdps":
                    self._wdps_names.add(a.asname or a.name)
                if a.name == "mpf":
                    self._mpf_names.add(a.asname or a.name)
        self.generic_visit(node)

    def visit_Assign(self, node):
        # track `M = mp.mpf` / `w = mp.workdps` style aliases
        if isinstance(node.value, (ast.Attribute, ast.Name)):
            tail = getattr(node.value, "attr", getattr(node.value, "id", ""))
            for t in node.targets:
                if isinstance(t, ast.Name):
                    if tail == "mpf":
                        self._mpf_names.add(t.id)
                    if tail == "workdps":
                        self._wdps_names.add(t.id)
        self.generic_visit(node)

    def _is_wdps_call(self, expr):
        if not isinstance(expr, ast.Call):
            return False
        name = getattr(expr.func, "attr", getattr(expr.func, "id", ""))
        return name in self._wdps_names or name == "workdps"

    def visit_With(self, node):
        is_wdps = any(self._is_wdps_call(it.context_expr) for it in node.items)
        self._wdps_depth += is_wdps
        self.generic_visit(node)
        self._wdps_depth -= is_wdps

    def visit_Call(self, node):
        name = getattr(node.func, "attr", getattr(node.func, "id", ""))
        if name in self._mpf_names and self._wdps_depth == 0 and node.args:
            suppressed = (0 < node.lineno <= len(getattr(self, "_lines", []))
                          and "baller:ambient-ok" in self._lines[node.lineno - 1])
            if not isinstance(node.args[0], ast.Constant) and not suppressed:
                self.findings.append(
                    (node.lineno, "mpf() conversion OUTSIDE any workdps block"
                     " — rounds the input at ambient dps; both two-precision"
                     " legs then share the rounded input"))
        self.generic_visit(node)


def lint_ambient_conv(path, text=None):
    """Flag non-literal mpf() conversions outside every workdps block. Scope:
    files with two-precision structure — either workdps usage OR the
    `mp.dps = ...` assignment spelling (R1: the most common form of the
    FIXES.md:31 trap has no workdps at all). `# baller:ambient-ok` on the
    line suppresses (documented-deliberate ambient conversions)."""
    text = open(path).read() if text is None else text
    if "workdps" not in text and not re.search(r"mp\.dps\s*=", text):
        return []
    try:
        tree = ast.parse(text)
    except SyntaxError as e:
        return [(path, 0, f"SYNTAX (unlintable): {e}")]
    v = _AmbientConv()
    v._lines = text.split("\n")
    v.visit(tree)
    return [(path, ln, msg) for ln, msg in v.findings]


class _UnkeyedCache(ast.NodeVisitor):
    AMBIENT = ("mp.dps", "ctx.prec", "getprec", "prec()", "mp.prec",
               "workdps", "mp.workprec")

    def __init__(self):
        self.findings = []

    def visit_FunctionDef(self, node):
        deco = [ast.unparse(d) for d in node.decorator_list]
        # functools.cache (the modern lru_cache(maxsize=None)
        # spelling) must be caught alongside lru_cache
        if any("lru_cache" in d or d == "cache" or d.endswith(".cache")
               or d.startswith("cache(") for d in deco):
            body_src = ast.unparse(node)
            args = {a.arg for a in node.args.args}
            keyed = any(k in args for k in ("dps", "prec", "wp", "digits"))
            reads_ambient = any(t in body_src for t in self.AMBIENT)
            if reads_ambient and not keyed:
                self.findings.append(
                    (node.lineno, f"memo-cached function {node.name}() reads"
                     " ambient precision but is not keyed by dps/prec — the"
                     " first low-dps call silently caps every later result"))
        self.generic_visit(node)

    visit_AsyncFunctionDef = visit_FunctionDef


def lint_unkeyed_cache(path, text=None):
    text = open(path).read() if text is None else text
    if "lru_cache" not in text and "cache" not in text:
        return []
    return _findings_from_ast(path, text, _UnkeyedCache)


class _CtxRestore(ast.NodeVisitor):
    """Proper scope tracking (a ctx.prec/ctx.cap
    write is RESTORED iff it sits in a finalbody, or in the body of a try
    whose finalbody also writes ctx state); AugAssign (ctx.prec += 64,
    the precision-doubling idiom) and tuple-form restores
    (`ctx.prec, ctx.cap = old` — ctx_guard's own idiom) are now seen; the
    base name must be 'ctx' or a tracked `from flint import ctx as X` alias
    (the endswith('c') heuristic flagged any `vec.cap`); a line-comment
    `# baller:ctx-ok` suppresses (for guard classes that restore across
    method boundaries, e.g. __enter__/__exit__ pairs)."""

    def __init__(self, src_lines=None):
        self.findings = []
        self._protected = 0
        self._in_finally = 0
        self._ctx_names = {"ctx"}
        self._lines = src_lines or []

    def visit_ImportFrom(self, node):
        if node.module == "flint":
            for a in node.names:
                if a.name == "ctx":
                    self._ctx_names.add(a.asname or a.name)
        self.generic_visit(node)

    def _is_ctx_target(self, t):
        if isinstance(t, ast.Tuple):
            return any(self._is_ctx_target(e) for e in t.elts)
        return (isinstance(t, ast.Attribute) and t.attr in ("prec", "cap")
                and (getattr(t.value, "id", getattr(t.value, "attr", ""))
                     in self._ctx_names))

    def _writes_ctx(self, node):
        if isinstance(node, ast.Assign):
            return any(self._is_ctx_target(t) for t in node.targets)
        if isinstance(node, ast.AugAssign):
            return self._is_ctx_target(node.target)
        return False

    def _finalbody_restores(self, node):
        for stmt in node.finalbody:
            for n in ast.walk(stmt):
                if self._writes_ctx(n):
                    return True
        return False

    def visit_Try(self, node):
        restores = self._finalbody_restores(node)
        self._protected += restores
        for n in node.body + node.orelse:
            self.visit(n)
        self._protected -= restores
        for h in node.handlers:
            self.visit(h)
        self._in_finally += 1
        for n in node.finalbody:
            self.visit(n)
        self._in_finally -= 1

    def _flag(self, node):
        if self._protected or self._in_finally:
            return
        if 0 < node.lineno <= len(self._lines) and \
                "baller:ctx-ok" in self._lines[node.lineno - 1]:
            return
        t = node.targets[0] if isinstance(node, ast.Assign) else node.target
        self.findings.append(
            (node.lineno, f"{ast.unparse(t)} written with NO restoring "
             "try/finally in scope — process-global flint state; use "
             "baller.hygiene.ctx_guard or save/restore (acbfast law)"))

    def visit_Assign(self, node):
        if self._writes_ctx(node):
            self._flag(node)
        self.generic_visit(node)

    def visit_AugAssign(self, node):
        if self._writes_ctx(node):
            self._flag(node)
        self.generic_visit(node)


def lint_ctx_restore(path, text=None):
    """Flag ctx.prec/ctx.cap writes (incl. augmented) outside any restoring
    try/finally scope. `# baller:ctx-ok` on the line suppresses (guard
    classes restoring across method boundaries). NOT SEEN (documented reach
    limit): setattr(ctx, 'prec', v) — dataflow beyond a lint."""
    text = open(path).read() if text is None else text
    if ".prec" not in text and ".cap" not in text:
        return []
    try:
        tree = ast.parse(text)
    except SyntaxError as e:
        return [(path, 0, f"SYNTAX (unlintable): {e}")]
    v = _CtxRestore(src_lines=text.split("\n"))
    v.visit(tree)
    return [(path, ln, msg) for ln, msg in v.findings]


LINTS = {"fixed_wdps": lint_fixed_wdps, "ambient_conv": lint_ambient_conv,
         "unkeyed_cache": lint_unkeyed_cache, "ctx_restore": lint_ctx_restore}


def lint_files(paths, which=None):
    """Run the selected lints (default all) over paths; returns findings
    list; empty = clean (the dps_lint exit-code convention is the caller's)."""
    out = []
    for p in paths:
        text = open(p).read()
        for name, fn in LINTS.items():
            if which and name not in which:
                continue
            out.extend((name,) + f for f in fn(p, text))
    return out


# ---------- runtime helpers ----------

class RadiusBlowup(RuntimeError):
    pass


def ftrim(vec):
    """Ball midpoints once per ACCEPTED step (radii never accumulate).
    LAW: only after the step was validated/accepted — trimming an
    unvalidated ball discards rigor. Accepts arb/acb scalars or lists;
    non-ball scalars (int/float/Fraction — already zero-radius semantically)
    pass through unchanged (R1: the old fallback crashed inside its own
    except handler)."""
    def one(x):
        mid_attr = getattr(x, "mid", None)
        if mid_attr is None:
            return x                      # exact scalar: nothing to trim
        m = mid_attr() if callable(mid_attr) else mid_attr   # iv.mpf has
        try:                              # .mid as a PROPERTY, not a method
            return type(x)(m)
        except (AttributeError, TypeError):
            try:
                return type(x)(m.real, m.imag)
            except (AttributeError, TypeError):
                return m                  # reconstruction impossible: the mid
    return [one(x) for x in vec] if isinstance(vec, (list, tuple)) else one(vec)


class RadiusWatch:
    """Halts when relative radius crosses the bar: log10(rad/|mid|) > -bar
    digits. Use as the in-run alert lever against radius blowup."""

    def __init__(self, bar_digits, label="radius-watch"):
        self.bar, self.label = float(bar_digits), label
        self.worst = None

    def check(self, ball, where=""):
        import math
        rad = float(ball.rad()) if hasattr(ball, "rad") else float(abs(ball.imag.rad()) + abs(ball.real.rad()))
        mid = ball.mid() if hasattr(ball, "mid") else ball
        m = abs(complex(mid)) if hasattr(mid, "imag") else abs(float(mid))
        if not math.isfinite(m):
            # arb('inf') has rad EXACTLY 0 — the rad==0 shortcut must
            # exclude nonfinite midpoints (it would certify an overflowed
            # value as 'exact'); non-finite mid halts BEFORE any shortcut
            digits = -math.inf
            self.worst = digits if self.worst is None else min(self.worst, digits)
            raise RadiusBlowup(
                f"{self.label}{' @' + where if where else ''}: NON-FINITE "
                f"midpoint ({mid}) — HALT")
        if rad == 0:
            return 9999.0
        if not math.isfinite(rad) or not math.isfinite(m):
            digits = -math.inf            # inf/nan anywhere: always halt
        elif m == 0:
            # zero midpoint = TOTAL cancellation = infinite relative
            # radius — the exact case this watch exists for (treating mid=0
            # as mid=1 would report a healthy number)
            digits = -math.inf
        else:
            digits = -(math.log10(rad) - math.log10(m))
        self.worst = digits if self.worst is None else min(self.worst, digits)
        if not (digits >= self.bar):
            raise RadiusBlowup(
                f"{self.label}{' @' + where if where else ''}: relative "
                f"radius at {digits:.1f} digits < bar {self.bar} — HALT "
                f"(ball={ball})")
        return digits


class ctx_guard:
    """Save/restore python-flint's process-global ctx.prec/ctx.cap around a
    block (the vendored kklt
    engines set ctx.prec at function tops with no restore — callers wrap
    engine calls in this guard; the acbfast law as a context manager)."""

    def __init__(self, prec=None, cap=None):
        self._prec, self._cap = prec, cap

    def __enter__(self):
        from flint import ctx
        self._old = (ctx.prec, ctx.cap)
        if self._prec is not None:
            ctx.prec = self._prec         # baller:ctx-ok (restored in __exit__)
        if self._cap is not None:
            ctx.cap = self._cap           # baller:ctx-ok (restored in __exit__)
        return ctx

    def __exit__(self, *exc):
        from flint import ctx
        ctx.prec, ctx.cap = self._old     # baller:ctx-ok (the restore itself)
        return False


def bridge_mpf(x, mp):
    """Re-round-safe bridge for an EXISTING mpf. BOTH mp.mpf(x) and
    mp.mpf(x._mpf_) re-round to ambient dps (mpmath normalizes tuples to the
    current precision) — the raw-wrap constructor make_mpf is the only
    no-reround path ('mantissas are dyadic; no re-rounding', eval_int_J2L
    law; measured: a 50-dps sqrt2 survives make_mpf at dps 15 bit-exact and
    loses 115 bits through mpf())."""
    if hasattr(x, "_mpf_"):
        return mp.make_mpf(x._mpf_)
    return mp.mpf(x)   # baller:ambient-ok (non-mpf input: rounding is the
    #                     caller's stated intent at this point)
