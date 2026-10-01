#!/usr/bin/env python3
"""emitall.lint — template lint for emission scripts.

    python3 lint.py <script-or-dir> [more paths...] [--json report.json] [--quiet]

Enforces the emission-integrity rule: digit literals in emission prose are
findings, sampled fields must be sorted/seeded — it turns this error class
into a machine-caught one for every project.

Failure shape: an emitted headline carries the hardcoded prose
"(max 5 cents)" while the stored receipt's maximum is $0.09 — every OTHER
number in that string was interpolated from the receipt, so the one typed
digit is invisible to eyeball review ("typed-not-emitted"; such defects can
UNDERSTATE the receipt). Second class: emission order depending on
PYTHONHASHSEED (set/dict-hash iteration feeding emitted lists) — a double-run
byte-stability gate catches it downstream; this lint catches the pattern at
the source.

What it checks
--------------
(a) digit-literal: any digit run in the LITERAL text of a string that feeds
    an emission sink — f-string literal parts (interpolation placeholders are
    exempt, including format specs), ``.format()`` templates (``{...}`` fields
    stripped, ``{{``/``}}`` escapes kept as literal text), %-templates
    (placeholders incl. width/precision stripped), and plain string literals.
    An emission sink is, heuristically:
      - assignment / augmented assignment to a name, attribute, or string
        subscript key matching /emit|headline|summar/i;
      - .append/.extend/.insert/.add on such a name;
      - any argument to a call whose function name matches the pattern
        (a "registered emit function", e.g. emit(...), emit_headline(...));
      - a keyword argument or dict key matching the pattern;
      - the return value of a function whose own name matches the pattern.
    Escape hatch: ``# emitall-lint: allow-literal <reason>`` on any line of
    the flagged string expression. The reason is MANDATORY — a bare pragma is
    itself a finding (pragma-missing-reason) and suppresses nothing.

(b) unordered-iteration: set iteration (set()/frozenset(), set literals and
    comprehensions, &|^- combinations, .union()-family calls) or dict-view
    iteration (.keys()/.values()/.items()) feeding an emission sink — in the
    sink expression itself, inside its f-string interpolations, or a
    statement-level ``for`` over such an iterable whose body writes a sink —
    unless wrapped in sorted(). Order-insensitive reducers
    (sum/min/max/len/any/all/set/frozenset/sorted) are exempt consumers.

(c) unseeded-sampling: random.sample/choice/choices/shuffle (module-level
    ``random`` or ``*.random``, e.g. np.random) with no preceding .seed()
    call in the file, or the same methods on an instance constructed with no
    seed (random.Random(), np.random.default_rng(), RandomState()). Seeded
    instances and post-seed calls pass. Applies file-wide: emission scripts
    are the intended scan target, and an unseeded sample anywhere in one is
    the leak class.

Exit codes: 0 clean, 1 any finding, 2 usage error. Findings are listed,
never smoothed; pragma-allowed literals are counted and reported separately.

Honest limits: sink detection is by NAME pattern — an emission list called
``lines`` is invisible, and a quote-transcription argument to an emit-named
function is flagged even though it is a comparison target, not prose (fix:
pragma with reason, or move the quotes into an emitall claims spec — spec
JSON is data, not linted Python). Templates held in variables before
.format()/% are not resolved. Iterating a bare dict NAME cannot be
classified. An unresolvable ``x.sample()`` receiver (e.g. a pandas frame) is
not flagged. Seeding is detected per file, not per RNG stream.
"""
import argparse
import ast
import json
import os
import re
import subprocess
import sys

EMIT_RE = re.compile(r"emit|headline|summar", re.I)
SAMPLERS = {"sample", "choice", "choices", "shuffle"}
RNG_CTORS = {"Random", "default_rng", "RandomState", "SystemRandom", "Generator"}
ORDER_SAFE_CALLS = {"sorted", "sum", "min", "max", "len", "any", "all",
                    "set", "frozenset"}
ITER_UNWRAP = {"list", "tuple", "reversed", "enumerate", "iter"}
SETOPS = (ast.BitAnd, ast.BitOr, ast.BitXor, ast.Sub)
SET_METHODS = {"union", "intersection", "difference", "symmetric_difference"}
DICT_VIEWS = {"keys", "values", "items"}

DIGIT_RE = re.compile(r"\d+")
PRAGMA_RE = re.compile(r"#\s*emitall-lint:\s*allow-literal\b[:\s]*(.*)")
_FMT_FIELD = re.compile(r"\{[^{}]*\}")
_PCT_FIELD = re.compile(
    r"%(?:\([^)]*\))?[-+ #0]*(?:\*|\d+)?(?:\.(?:\*|\d+))?[hlL]?"
    r"[diouxXeEfFgGcrsa%]")


def strip_format_placeholders(text):
    """Remove {field} placeholders; keep {{ }} escapes as literal braces."""
    t = text.replace("{{", "\x00").replace("}}", "\x01")
    t = _FMT_FIELD.sub("", t)
    return t.replace("\x00", "{").replace("\x01", "}")


def strip_percent_placeholders(text):
    """Remove %-placeholders (incl. width/precision digits); %% -> literal."""
    return _PCT_FIELD.sub(lambda m: "" if m.group() != "%%" else "%", text)


def _name_of(node):
    """Best-effort name of an assignment target / call receiver."""
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        return node.attr
    if isinstance(node, ast.Subscript):
        s = node.slice
        if isinstance(s, ast.Constant) and isinstance(s.value, str):
            return s.value
        return _name_of(node.value)
    return None


def _flat_targets(t):
    if isinstance(t, (ast.Tuple, ast.List)):
        for e in t.elts:
            yield from _flat_targets(e)
    else:
        yield t


def classify_iter(node):
    """'set' / 'dict-view' when the iterable's order is hash-dependent or
    insertion-opaque; None when acceptable (sorted, lists, unknown names)."""
    while (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
           and node.func.id in ITER_UNWRAP and node.args):
        node = node.args[0]
    if isinstance(node, (ast.Set, ast.SetComp)):
        return "set"
    if isinstance(node, ast.Call):
        if isinstance(node.func, ast.Name):
            if node.func.id in {"set", "frozenset"}:
                return "set"
            return None                      # sorted(...) and friends
        if isinstance(node.func, ast.Attribute):
            if node.func.attr in DICT_VIEWS:
                return "dict-view"
            if node.func.attr in SET_METHODS:
                return "set"
        return None
    if isinstance(node, ast.BinOp) and isinstance(node.op, SETOPS):
        if classify_iter(node.left) == "set" or classify_iter(node.right) == "set":
            return "set"
    return None


class FileLinter(ast.NodeVisitor):
    def __init__(self, path, source):
        self.path = path
        self.findings = []
        self.allowed = []
        self._seen = set()
        self._func_stack = []
        self.pragmas = {}                     # lineno -> reason text
        for i, line in enumerate(source.splitlines(), 1):
            m = PRAGMA_RE.search(line)
            if m:
                self.pragmas[i] = m.group(1).strip()
        self.tree = ast.parse(source, filename=path)
        self.seed_lines = []
        self.rng_bindings = {}                # name -> bool(seeded)
        for n in ast.walk(self.tree):
            if (isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute)
                    and n.func.attr == "seed"):
                self.seed_lines.append(n.lineno)
            if isinstance(n, ast.Assign) and len(n.targets) == 1 \
                    and isinstance(n.targets[0], ast.Name) \
                    and isinstance(n.value, ast.Call):
                fn = n.value.func
                ctor = fn.attr if isinstance(fn, ast.Attribute) else (
                    fn.id if isinstance(fn, ast.Name) else None)
                if ctor in RNG_CTORS:
                    self.rng_bindings[n.targets[0].id] = bool(
                        n.value.args or n.value.keywords)

    # ------------------------------------------------------------ findings

    def finding(self, kind, lineno, detail):
        key = (kind, lineno, detail)
        if key in self._seen:
            return
        self._seen.add(key)
        self.findings.append({"file": self.path, "line": lineno,
                              "kind": kind, "detail": detail})

    def check_text(self, text, tplnode, what):
        runs = list(DIGIT_RE.finditer(text))
        if not runs:
            return
        span = range(tplnode.lineno,
                     getattr(tplnode, "end_lineno", tplnode.lineno) + 1)
        pragma_lines = [ln for ln in span if ln in self.pragmas]
        if pragma_lines:
            reason = self.pragmas[pragma_lines[0]]
            if reason:
                self.allowed.append({
                    "file": self.path, "line": tplnode.lineno,
                    "runs": [m.group() for m in runs], "reason": reason})
                return
            self.finding("pragma-missing-reason", pragma_lines[0],
                         "allow-literal pragma without a reason — reason is "
                         "mandatory; pragma suppresses nothing")
        for m in runs:
            ctx = text[max(0, m.start() - 24):m.end() + 24].replace("\n", " ")
            self.finding("digit-literal", tplnode.lineno,
                         f"digit run {m.group()!r} in {what}: ...{ctx}...")

    def flag_if_unordered(self, iter_node):
        cls = classify_iter(iter_node)
        if cls:
            self.finding("unordered-iteration", iter_node.lineno,
                         f"{cls} iteration feeds an emission — order is "
                         "hash/insertion-dependent; wrap in sorted() (emission "
                         "order: sampled/enumerated fields sorted or seeded)")

    # -------------------------------------------------- sink value scanning

    def scan_value(self, node):
        if node is None:
            return
        if isinstance(node, ast.JoinedStr):
            for part in node.values:
                if isinstance(part, ast.Constant) and isinstance(part.value, str):
                    self.check_text(part.value, node, "f-string literal text")
                elif isinstance(part, ast.FormattedValue):
                    self.scan_order(part.value)      # format_spec exempt
            return
        if isinstance(node, ast.Constant):
            if isinstance(node.value, str):
                self.check_text(node.value, node, "string literal")
            return
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                and node.func.attr == "format"
                and isinstance(node.func.value, ast.Constant)
                and isinstance(node.func.value.value, str)):
            self.check_text(strip_format_placeholders(node.func.value.value),
                            node, ".format() template")
            for a in node.args:
                self.scan_order(a)
            for kw in node.keywords:
                self.scan_order(kw.value)
            return
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Mod) \
                and isinstance(node.left, ast.Constant) \
                and isinstance(node.left.value, str):
            self.check_text(strip_percent_placeholders(node.left.value),
                            node, "%-template")
            self.scan_order(node.right)
            return
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
            self.scan_value(node.left)
            self.scan_value(node.right)
            return
        if isinstance(node, (ast.List, ast.Tuple, ast.Set)):
            for e in node.elts:
                self.scan_value(e)
            return
        if isinstance(node, ast.Dict):
            for v in node.values:
                self.scan_value(v)
            return
        if isinstance(node, ast.IfExp):
            self.scan_value(node.body)
            self.scan_value(node.orelse)
            self.scan_order(node.test)
            return
        if isinstance(node, (ast.ListComp, ast.GeneratorExp, ast.SetComp)):
            self.scan_value(node.elt)
            for gen in node.generators:
                self.flag_if_unordered(gen.iter)
                for cond in gen.ifs:
                    self.scan_order(cond)
            return
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                and node.func.attr == "join"):
            if isinstance(node.func.value, ast.Constant) \
                    and isinstance(node.func.value.value, str):
                self.check_text(node.func.value.value, node, "join separator")
            for a in node.args:
                self.flag_if_unordered(a)
                self.scan_value(a) if isinstance(
                    a, (ast.GeneratorExp, ast.ListComp)) else self.scan_order(a)
            return
        self.scan_order(node)

    def scan_order(self, node):
        """Order-only scan (rule b) of an arbitrary expression subtree."""
        if node is None or isinstance(node, ast.Constant):
            return
        if isinstance(node, ast.Call):
            if isinstance(node.func, ast.Name) \
                    and node.func.id in ORDER_SAFE_CALLS:
                return                        # order-insensitive consumer
            if isinstance(node.func, ast.Attribute) and node.func.attr == "join":
                for a in node.args:
                    self.flag_if_unordered(a)
        if isinstance(node, (ast.ListComp, ast.GeneratorExp,
                             ast.SetComp, ast.DictComp)):
            for gen in node.generators:
                self.flag_if_unordered(gen.iter)
        for child in ast.iter_child_nodes(node):
            self.scan_order(child)

    # ------------------------------------------------------------- visitors

    def _is_sink_name(self, name):
        return bool(name and EMIT_RE.search(name))

    def visit_Assign(self, node):
        if any(self._is_sink_name(_name_of(t))
               for tgt in node.targets for t in _flat_targets(tgt)):
            self.scan_value(node.value)
        self.generic_visit(node)

    def visit_AnnAssign(self, node):
        if self._is_sink_name(_name_of(node.target)):
            self.scan_value(node.value)
        self.generic_visit(node)

    def visit_AugAssign(self, node):
        if self._is_sink_name(_name_of(node.target)):
            self.scan_value(node.value)
        self.generic_visit(node)

    def visit_Dict(self, node):
        for k, v in zip(node.keys, node.values):
            if isinstance(k, ast.Constant) and isinstance(k.value, str) \
                    and self._is_sink_name(k.value):
                self.scan_value(v)
        self.generic_visit(node)

    def visit_Call(self, node):
        fn = node.func
        if isinstance(fn, ast.Attribute) \
                and fn.attr in {"append", "extend", "insert", "add",
                                "appendleft"} \
                and self._is_sink_name(_name_of(fn.value)):
            for a in node.args:
                self.scan_value(a)
        elif self._is_sink_name(_name_of(fn)):
            for a in node.args:
                self.scan_value(a)
            for kw in node.keywords:
                self.scan_value(kw.value)
        else:
            for kw in node.keywords:
                if kw.arg and self._is_sink_name(kw.arg):
                    self.scan_value(kw.value)
        if isinstance(fn, ast.Attribute) and fn.attr in SAMPLERS:
            base = fn.value
            is_random_module = (
                (isinstance(base, ast.Name) and base.id == "random")
                or (isinstance(base, ast.Attribute) and base.attr == "random"))
            if is_random_module:
                if not any(sl < node.lineno for sl in self.seed_lines):
                    self.finding(
                        "unseeded-sampling", node.lineno,
                        f"random.{fn.attr}() with no preceding .seed() call "
                        "in this file — emission order/content varies per "
                        "run (emission rule: sampled fields sorted or seeded)")
            elif isinstance(base, ast.Name) and base.id in self.rng_bindings:
                if not self.rng_bindings[base.id]:
                    self.finding(
                        "unseeded-sampling", node.lineno,
                        f"{base.id}.{fn.attr}() on an RNG constructed with "
                        "no seed — pass an explicit seed")
        self.generic_visit(node)

    def visit_For(self, node):
        cls = classify_iter(node.iter)
        if cls and self._body_has_sink(node.body):
            self.finding("unordered-iteration", node.iter.lineno,
                         f"for-loop over a {cls} writes an emission sink — "
                         "order is hash/insertion-dependent; wrap the "
                         "iterable in sorted()")
        self.generic_visit(node)

    def _body_has_sink(self, stmts):
        for stmt in stmts:
            for n in ast.walk(stmt):
                if isinstance(n, (ast.Assign, ast.AnnAssign, ast.AugAssign)):
                    tgts = n.targets if isinstance(n, ast.Assign) else [n.target]
                    for t in tgts:
                        for tt in _flat_targets(t):
                            if self._is_sink_name(_name_of(tt)):
                                return True
                elif isinstance(n, ast.Call):
                    f = n.func
                    if isinstance(f, ast.Attribute) \
                            and f.attr in {"append", "extend", "insert",
                                           "add", "appendleft"} \
                            and self._is_sink_name(_name_of(f.value)):
                        return True
                    if self._is_sink_name(_name_of(f)):
                        return True
                    if any(kw.arg and self._is_sink_name(kw.arg)
                           for kw in n.keywords):
                        return True
                elif isinstance(n, ast.Dict):
                    if any(isinstance(k, ast.Constant)
                           and isinstance(k.value, str)
                           and self._is_sink_name(k.value) for k in n.keys):
                        return True
        return False

    def visit_FunctionDef(self, node):
        self._func_stack.append(node.name)
        self.generic_visit(node)
        self._func_stack.pop()

    visit_AsyncFunctionDef = visit_FunctionDef

    def visit_Return(self, node):
        if self._func_stack and self._is_sink_name(self._func_stack[-1]):
            self.scan_value(node.value)
        self.generic_visit(node)


# -------------------------------------------------------------------- API

def lint_source(source, path="<string>"):
    """Returns (findings, allowed) for one source text."""
    try:
        linter = FileLinter(path, source)
    except SyntaxError as e:
        return ([{"file": path, "line": e.lineno or 0, "kind": "parse-error",
                  "detail": f"cannot parse: {e.msg}"}], [])
    linter.visit(linter.tree)
    linter.findings.sort(key=lambda f: (f["line"], f["kind"], f["detail"]))
    return linter.findings, linter.allowed


def lint_file(path):
    with open(path, encoding="utf-8", errors="replace") as fh:
        return lint_source(fh.read(), path)


def iter_py_files(target):
    if os.path.isfile(target):
        yield target
        return
    for root, dirs, files in os.walk(target):
        dirs[:] = sorted(d for d in dirs
                         if d != "__pycache__" and not d.startswith("."))
        for name in sorted(files):
            if name.endswith(".py"):
                yield os.path.join(root, name)


def lint_paths(targets):
    """Returns (findings, allowed, files_scanned)."""
    findings, allowed, n_files = [], [], 0
    for target in targets:
        if not os.path.exists(target):
            findings.append({"file": target, "line": 0, "kind": "parse-error",
                             "detail": "path does not exist"})
            continue
        for path in iter_py_files(target):
            n_files += 1
            f, a = lint_file(path)
            findings.extend(f)
            allowed.extend(a)
    return findings, allowed, n_files


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="emitall template lint — digit literals in emission "
                    "prose and unsorted/unseeded sampled fields")
    ap.add_argument("targets", nargs="+", help="python script(s) or dir(s)")
    ap.add_argument("--json", default=None, help="write JSON report here")
    ap.add_argument("--quiet", action="store_true",
                    help="counts only, no per-finding lines")
    args = ap.parse_args(argv)
    findings, allowed, n_files = lint_paths(args.targets)
    counts = {}
    for f in findings:
        counts[f["kind"]] = counts.get(f["kind"], 0) + 1
    if not args.quiet:
        cur = None
        for f in findings:
            if f["file"] != cur:
                cur = f["file"]
                print(f"\n{cur}")
            print(f"  L{f['line']:<5d} {f['kind']:22s} {f['detail']}")
    summary = ", ".join(f"{v} {k}" for k, v in sorted(counts.items())) or "clean"
    print(f"\n{len(findings)} findings ({summary}) in {n_files} files scanned;"
          f" {len(allowed)} literal(s) allowed by pragma")
    for a in allowed:
        print(f"  allowed {a['file']}:L{a['line']} runs={a['runs']} "
              f"reason: {a['reason']}")
    if args.json:
        report = {
            "run": "emitall-lint",
            "date_utc": subprocess.run(
                ["date", "-u", "+%Y-%m-%dT%H:%M:%SZ"],
                capture_output=True, text=True).stdout.strip(),
            "targets": [os.path.abspath(t) for t in args.targets],
            "files_scanned": n_files,
            "counts": counts,
            "findings": findings,
            "allowed_by_pragma": allowed,
        }
        with open(args.json, "w") as fh:
            json.dump(report, fh, indent=1)
        print("report ->", args.json)
    return 1 if findings else 0


if __name__ == "__main__":
    sys.exit(main())
