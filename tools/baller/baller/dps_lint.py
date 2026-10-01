#!/usr/bin/env python3
"""
baller.dps_lint — static linter for the mpmath import-time-dps footgun
(member module of BALLER's hygiene wing; stdlib-only, runs standalone).

The bug (it has capped a production evaluator at 17.6 digits, twice):
    import mpmath as mp
    T = mp.mpf(-1)/3          # ← evaluated at import, dps=15 → 15-digit constant
    def f(x=mp.mpf('0.02')):  # ← default arg evaluated at def-time = import time
        mp.mp.dps = 100       # too late for T and x's default
        ...

This linter walks the module AST and flags any mpmath constructor / numeric
function / lazy-constant attribute that is *evaluated at import time* and
appears textually BEFORE the first module-level assignment to mp.dps /
mp.mp.dps / mpmath.mp.dps. Function/method *bodies* are NOT walked (they run
at call time), but function *default arguments*, *decorators*, and *class-body
assignments* ARE (they run at import time).

Cross-module rules (the import-resets-dps
footgun: module A sets mp.mp.dps at import time, module B imports it ->
B's dps silently reset mid-run, and the false bias can cost hours to find):
  [imported-sets-dps]   module-level `mp.mp.dps =` in a module IMPORTED by
                        another linted file (re-set dps after EVERY import).
  [import-frozen-const] module-level mpf/mpc/... in an IMPORTED module: it is
                        evaluated ONCE at import time at whatever dps holds
                        then, even though the module sets dps itself; a main
                        that raises dps later still inherits the frozen value.
  [pre-main-const]      module-level literal whose only protecting dps
                        assignment sits under `if __name__ == '__main__'` -
                        as an import the guard never runs (default dps=15).

Usage (any of the three; the module is stdlib-only):
    python3 tools/baller/baller/dps_lint.py FILE_OR_DIR [FILE_OR_DIR ...]
    python3 -m baller.hygiene lint FILE_OR_DIR [...]   (tools/baller on sys.path)
    from baller import hygiene; hygiene.dps_lint.main([...])
    ... --selftest   (built-in battery; also leg L21 of battery/battery.py)
Exit 1 if any findings.  Cross-module rules run over the whole file set given.
Exit 2 on a nonexistent path or a non-.py file argument.
--selftest plants one finding per rule in toy fixtures under a temp dir
(honors TMPDIR) and asserts the expected counts and exit codes (exit 0
all-pass, 1 otherwise).
"""
import ast, os, sys

# mpmath callables whose return value's precision is fixed at call time
MP_FUNCS = {
    'mpf', 'mpc', 'mpmathify', 'matrix', 'mpi',
    'log', 'ln', 'log10', 'exp', 'sqrt', 'cbrt', 'power', 'root',
    'sin', 'cos', 'tan', 'asin', 'acos', 'atan', 'atan2',
    'sinh', 'cosh', 'tanh', 'gamma', 'loggamma', 'zeta', 'polylog',
    'ellipk', 'ellipe', 'ellipf', 'hyp2f1', 'besselj', 'besselk',
    'fac', 'factorial', 'binomial', 'bernoulli',
}
# mpmath lazy-constant attributes (mpmath properties that compute at access time)
MP_CONSTS = {
    'pi', 'e', 'euler', 'catalan', 'apery', 'phi', 'glaisher',
    'khinchin', 'mertens', 'twinprime', 'ln2', 'ln10', 'eps', 'inf',
}


def _is_main_guard(test):
    """True for `__name__ == '__main__'` (either operand order)."""
    if not (isinstance(test, ast.Compare) and len(test.comparators) == 1):
        return False
    nm = lambda n: isinstance(n, ast.Name) and n.id == '__name__'
    mn = lambda n: isinstance(n, ast.Constant) and n.value == '__main__'
    l, r = test.left, test.comparators[0]
    return (nm(l) and mn(r)) or (nm(r) and mn(l))


def _attr_chain(node):
    """Return dotted-name list for an Attribute/Name chain, or None."""
    parts = []
    while isinstance(node, ast.Attribute):
        parts.append(node.attr)
        node = node.value
    if isinstance(node, ast.Name):
        parts.append(node.id)
        return list(reversed(parts))
    return None


class DpsLinter(ast.NodeVisitor):
    def __init__(self, src_lines, fname, collect_all=False):
        self.src = src_lines
        self.fname = fname
        # names that refer to the mpmath module or its mp context
        self.mp_aliases = set()
        # bare names imported from mpmath (mpf, mpc, ...)
        self.mp_bare = set()
        self.dps_set_line = None  # first module-level dps assignment line
        self.findings = []
        # cross-module extension: keep walking past dps sets,
        # tag each finding with its protection state, track imports
        self.collect_all = collect_all
        self._protect = None          # None | 'main' | 'module'
        self.dps_module_lines = []    # true module-level dps assignments
        self.dps_main_lines = []      # dps assignments under __main__ guard
        self.imported = set()         # top-level names of imported modules

    # ── import tracking ──────────────────────────────────────────────
    def visit_Import(self, node):
        for a in node.names:
            self.imported.add(a.name.split('.')[0])
            if a.name == 'mpmath':
                self.mp_aliases.add(a.asname or 'mpmath')

    def visit_ImportFrom(self, node):
        if node.module:
            self.imported.add(node.module.split('.')[0])
        if node.module and node.module.split('.')[0] == 'mpmath':
            for a in node.names:
                nm = a.asname or a.name
                if a.name == 'mp':
                    self.mp_aliases.add(nm)
                else:
                    self.mp_bare.add(nm)

    # ── dps assignment detection ─────────────────────────────────────
    def _is_dps_target(self, tgt):
        ch = _attr_chain(tgt)
        if not ch or ch[-1] != 'dps':
            return False
        # mp.dps  /  mp.mp.dps  /  mpmath.mp.dps  /  <alias>.dps
        return ch[0] in self.mp_aliases or ch[0] == 'mpmath'

    # ── mpmath-expression detection ──────────────────────────────────
    def _flag_if_mpmath(self, node, ctx):
        if isinstance(node, ast.Call):
            f = node.func
            ch = _attr_chain(f)
            if ch:
                if len(ch) == 1 and ch[0] in self.mp_bare and ch[0] in MP_FUNCS:
                    self._emit(node, ctx, '.'.join(ch))
                    return
                if ch[0] in self.mp_aliases and ch[-1] in MP_FUNCS:
                    self._emit(node, ctx, '.'.join(ch))
                    return
        elif isinstance(node, ast.Attribute):
            ch = _attr_chain(node)
            if ch and ch[0] in self.mp_aliases and ch[-1] in MP_CONSTS:
                # skip if this Attribute is the .func of a Call (handled above)
                self._emit(node, ctx, '.'.join(ch))

    def _emit(self, node, ctx, what):
        ln = node.lineno
        snip = self.src[ln - 1].rstrip() if ln - 1 < len(self.src) else ''
        self.findings.append((self.fname, ln, ctx, what, snip, self._protect))

    def _note_dps(self, lineno, ctx):
        """Record a dps assignment; return True if the (legacy) walk stops."""
        in_main = '__main__-guard' in ctx
        (self.dps_main_lines if in_main else self.dps_module_lines).append(lineno)
        if self.dps_set_line is None:
            self.dps_set_line = lineno
        if self._protect is None or (self._protect == 'main' and not in_main):
            self._protect = 'main' if in_main else 'module'
        return not self.collect_all

    def _scan_expr(self, node, ctx):
        """Walk an expression tree, pruning lambda bodies."""
        self._scan_expr_prune_lambda(node, ctx)

    # ── top-level driver ─────────────────────────────────────────────
    def lint(self, tree):
        self._walk_block(tree.body, 'module-level')

    def _walk_block(self, stmts, ctx):
        """Sequential walk of a statement list at import/script time.
        Descends into If/For/While/With/Try bodies (those execute in order),
        but NOT into FunctionDef/Lambda bodies (those execute at call time).
        Stops flagging once a dps assignment is seen in this linear path.
        Returns True if a dps assignment was encountered (so caller can stop too).
        """
        for stmt in stmts:
            # imports
            if isinstance(stmt, (ast.Import, ast.ImportFrom)):
                self.visit(stmt)
                continue
            # dps assignment?  (Assign / AugAssign; also handle `mp.mp.dps = N; X` via Expr-of-Tuple? no — semicolons split into separate stmts)
            if isinstance(stmt, (ast.Assign, ast.AugAssign)):
                tgts = stmt.targets if isinstance(stmt, ast.Assign) else [stmt.target]
                if any(self._is_dps_target(t) for t in tgts):
                    if self._note_dps(stmt.lineno, ctx):
                        return True
                    continue
                # otherwise scan RHS
                if getattr(stmt, 'value', None) is not None:
                    self._scan_expr(stmt.value, ctx)
                continue
            if isinstance(stmt, ast.AnnAssign):
                if self._is_dps_target(stmt.target):
                    if self._note_dps(stmt.lineno, ctx):
                        return True
                    continue
                if stmt.value is not None:
                    self._scan_expr(stmt.value, ctx)
                continue
            # function def: only defaults+decorators are import-time
            if isinstance(stmt, (ast.FunctionDef, ast.AsyncFunctionDef)):
                for d in (stmt.args.defaults + stmt.args.kw_defaults):
                    if d is not None:
                        self._scan_expr(d, f'default-arg in {stmt.name}()')
                for dec in stmt.decorator_list:
                    self._scan_expr(dec, f'decorator on {stmt.name}()')
                continue
            if isinstance(stmt, ast.ClassDef):
                for dec in stmt.decorator_list:
                    self._scan_expr(dec, f'decorator on class {stmt.name}')
                if self._walk_block(stmt.body, f'class-body {stmt.name}'):
                    return True
                continue
            # compound statements: descend into bodies in order
            if isinstance(stmt, ast.If):
                self._scan_expr(stmt.test, ctx)
                # We treat both branches conservatively: if EITHER branch sets dps,
                # we consider dps set after the If (optimistic for `if __name__` guards;
                # if you want strict, change to `body_set and else_set`).
                # __main__-guard bodies are tagged so a dps set there counts as
                # 'main' protection only (rule [pre-main-const]).
                bctx = (ctx + ' __main__-guard') if _is_main_guard(stmt.test) else ctx
                b = self._walk_block(stmt.body, bctx)
                e = self._walk_block(stmt.orelse, ctx)
                if b or e:
                    return True
                continue
            if isinstance(stmt, (ast.For, ast.While)):
                for f in ast.iter_fields(stmt):
                    if f[0] in ('iter', 'test') and f[1] is not None:
                        self._scan_expr(f[1], ctx)
                if self._walk_block(stmt.body, ctx):
                    return True
                self._walk_block(stmt.orelse, ctx)
                continue
            if isinstance(stmt, ast.With):
                for it in stmt.items:
                    self._scan_expr(it.context_expr, ctx)
                if self._walk_block(stmt.body, ctx):
                    return True
                continue
            if isinstance(stmt, ast.Try):
                if self._walk_block(stmt.body, ctx):
                    return True
                for h in stmt.handlers:
                    self._walk_block(h.body, ctx)
                self._walk_block(stmt.orelse, ctx)
                self._walk_block(stmt.finalbody, ctx)
                continue
            # simple statement (Expr, Return, Raise, ...): scan its expressions,
            # skipping any nested lambda bodies.
            self._scan_stmt_shallow(stmt, ctx)
        return False

    def _scan_stmt_shallow(self, stmt, ctx):
        for field, val in ast.iter_fields(stmt):
            if val is None:
                continue
            for node in (val if isinstance(val, list) else [val]):
                if isinstance(node, ast.expr):
                    self._scan_expr_prune_lambda(node, ctx)

    def _scan_expr_prune_lambda(self, node, ctx):
        class V(ast.NodeVisitor):
            def __init__(s, o, c): s.o=o; s.c=c
            def visit_Lambda(s, n):
                for d in (n.args.defaults + n.args.kw_defaults):
                    if d is not None: s.o._scan_expr(d, 'default-arg in lambda')
                # skip body
            def generic_visit(s, n):
                s.o._flag_if_mpmath(n, s.c)
                ast.NodeVisitor.generic_visit(s, n)
        V(self, ctx).visit(node)


# ── rule [mpdiff-real-branched] ──────────────────────────────────────────
# mp.diff(f, x) uses complex-step differentiation.  If f (or something f
# calls) branches on real order (assert, min/max, sorted, `if a < b:`),
# the complex step silently takes the wrong branch or raises TypeError.
# Heuristic AST check; false positives OK.
_BRANCH_BUILTINS = {'min', 'max', 'sorted', 'abs'}


def _has_real_branch(fn_node):
    """True if fn_node's body contains assert / min / max / sorted / an
    `if <cmp>:` on a Compare with </<=/>/>=."""
    for n in ast.walk(fn_node):
        if isinstance(n, ast.Assert):
            return True
        if isinstance(n, ast.Call):
            fid = getattr(n.func, 'id', None)
            if fid in _BRANCH_BUILTINS:
                return True
        if isinstance(n, ast.If) and isinstance(n.test, ast.Compare):
            if any(isinstance(op, (ast.Lt, ast.LtE, ast.Gt, ast.GtE))
                   for op in n.test.ops):
                return True
    return False


def _lint_mpdiff(tree, src_lines, fname, mp_aliases, mp_bare):
    """Scan whole tree for mp.diff(f,...) where f is a Name/Lambda that
    (heuristically) branches on real order.  False positives OK."""
    # index all module-local function defs by name (including nested)
    defs = {n.name: n for n in ast.walk(tree)
            if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))}
    branching = {nm for nm, fn in defs.items() if _has_real_branch(fn)}
    out = []
    for n in ast.walk(tree):
        if not isinstance(n, ast.Call):
            continue
        ch = _attr_chain(n.func)
        is_diff = (ch and ch[-1] == 'diff'
                   and (ch[0] in mp_aliases
                        or (len(ch) == 1 and 'diff' in mp_bare)))
        if not is_diff or not n.args:
            continue
        arg0 = n.args[0]
        hit = None
        if isinstance(arg0, ast.Lambda):
            # lambda body: branch in body OR calls a branching local def
            if _has_real_branch(arg0):
                hit = 'lambda body'
            else:
                for c in ast.walk(arg0.body):
                    if isinstance(c, ast.Call):
                        cid = getattr(c.func, 'id', None)
                        if cid in branching:
                            hit = f'calls {cid}()'
                            break
        elif isinstance(arg0, ast.Name) and arg0.id in branching:
            hit = f'{arg0.id}()'
        if hit:
            snip = src_lines[n.lineno - 1].rstrip() \
                if n.lineno - 1 < len(src_lines) else ''
            out.append((fname, n.lineno, 'mpdiff-real-branched',
                        f'mp.diff on real-branched callable ({hit}) — '
                        'complex-step will misbranch', snip))
    return out


def lint_file(path):
    try:
        src = open(path, encoding='utf-8', errors='replace').read()
        tree = ast.parse(src, filename=path)
    except SyntaxError as e:
        return [(path, e.lineno or 0, 'SYNTAX', '', str(e))]
    L = DpsLinter(src.splitlines(), path)
    L.lint(tree)
    if not L.mp_aliases and not L.mp_bare:
        return []  # no mpmath import → nothing to flag
    out = [f[:5] for f in L.findings]
    out += _lint_mpdiff(tree, L.src, path, L.mp_aliases, L.mp_bare)
    return out


def lint_cross(files):
    """Cross-module rules over a file SET: [imported-sets-dps],
    [import-frozen-const] (imported modules), [pre-main-const] (any file).
    Returns findings as (file, line, rule, what, snippet)."""
    facts = {}
    for f in files:
        try:
            src = open(f, encoding='utf-8', errors='replace').read()
            tree = ast.parse(src, filename=f)
        except SyntaxError:
            continue
        L = DpsLinter(src.splitlines(), f, collect_all=True)
        L.lint(tree)
        if L.mp_aliases or L.mp_bare or L.imported:
            facts[f] = L
    byname = {os.path.splitext(os.path.basename(f))[0]: f for f in facts}
    importers = {}
    for f, L in facts.items():
        for mod in L.imported:
            tgt = byname.get(mod)
            if tgt and tgt != f:
                importers.setdefault(tgt, set()).add(os.path.basename(f))
    out = []
    for f, L in facts.items():
        imps = sorted(importers.get(f, ()))
        snip = lambda ln: L.src[ln - 1].strip() if ln - 1 < len(L.src) else ''
        if imps:
            for ln in L.dps_module_lines:
                out.append((f, ln, 'imported-sets-dps',
                            'imported by: ' + ', '.join(imps), snip(ln)))
            for (_fn, ln, ctx, what, sn, prot) in L.findings:
                if prot == 'module':  # prot None already a per-file finding
                    out.append((f, ln, 'import-frozen-const',
                                f'{what} frozen at import dps (imported by: '
                                + ', '.join(imps) + ')', sn))
        for (_fn, ln, ctx, what, sn, prot) in L.findings:
            if prot == 'main':
                out.append((f, ln, 'pre-main-const',
                            f'{what} — dps only set under __main__ guard', sn))
    return out


def gather(paths):
    files = []
    for p in paths:
        if os.path.isfile(p) and p.endswith('.py'):
            files.append(p)
        elif os.path.isdir(p):
            for root, _, fs in os.walk(p):
                for f in fs:
                    if f.endswith('.py'):
                        files.append(os.path.join(root, f))
    return sorted(set(files))


# ── built-in battery ─────────────────────────────────────────────────────
def selftest():
    """Plant exactly one finding per rule in toy fixture files under a temp
    dir, then assert: per-rule counts land in the expected files, the CLI
    exits 1 on the planted set, 0 on a clean set, and 2 on bad path
    arguments.  A fixed copy of the planted constant (dps set first) must
    lint clean, so the pass state is reachable.  Returns 0 all-pass, 1 else."""
    import io, shutil, tempfile
    from contextlib import redirect_stdout, redirect_stderr
    planted_src = {
        # [module-level]: constructor evaluated before the first dps set
        'const_predps.py':
            "import mpmath as mp\n"
            "T = mp.mpf(-1) / 3\n"
            "mp.mp.dps = 50\n",
        # [imported-sets-dps] + [import-frozen-const]: module sets dps and
        # freezes a constant, and another linted file imports it
        'helper_mod.py':
            "import mpmath as mp\n"
            "mp.mp.dps = 40\n"
            "C = mp.mpf('0.1')\n",
        'uses_helper.py':
            "import helper_mod\n"
            "import mpmath as mp\n"
            "mp.mp.dps = 60\n"
            "X = mp.mpf('2')\n",
        # [pre-main-const]: only protecting dps set sits under __main__
        'main_guarded.py':
            "import mpmath as mp\n"
            "if __name__ == '__main__':\n"
            "    mp.mp.dps = 50\n"
            "T = mp.pi\n",
        # [mpdiff-real-branched]: mp.diff on a callable that branches on
        # real order
        'diff_branch.py':
            "import mpmath as mp\n"
            "mp.mp.dps = 30\n"
            "def f(x):\n"
            "    if x < 0:\n"
            "        return -x\n"
            "    return x\n"
            "d = mp.diff(f, 1)\n",
    }
    clean_src = (
        "import mpmath as mp\n"
        "mp.mp.dps = 50\n"
        "Z = mp.mpf('2') / 3\n"
        "def g(x):\n"
        "    return mp.sqrt(x)\n")
    expect = {
        'module-level': 'const_predps.py',
        'imported-sets-dps': 'helper_mod.py',
        'import-frozen-const': 'helper_mod.py',
        'pre-main-const': 'main_guarded.py',
        'mpdiff-real-branched': 'diff_branch.py',
    }
    ok, fails = [], []
    tmp = tempfile.mkdtemp(prefix='dps_lint_selftest_')
    try:
        planted = os.path.join(tmp, 'planted')
        clean = os.path.join(tmp, 'clean')
        os.makedirs(planted); os.makedirs(clean)
        for name, src in planted_src.items():
            with open(os.path.join(planted, name), 'w') as fh:
                fh.write(src)
        with open(os.path.join(clean, 'clean_mod.py'), 'w') as fh:
            fh.write(clean_src)

        # per-rule counts over the planted set (per-file + cross rules)
        files = gather([planted])
        counts, where = {}, {}
        for f in files:
            for (fn, ln, rule, what, snip) in lint_file(f):
                counts[rule] = counts.get(rule, 0) + 1
                where.setdefault(rule, set()).add(os.path.basename(fn))
        for (fn, ln, rule, what, snip) in lint_cross(files):
            counts[rule] = counts.get(rule, 0) + 1
            where.setdefault(rule, set()).add(os.path.basename(fn))
        for rule, fname in expect.items():
            n = counts.get(rule, 0)
            if n != 1:
                fails.append(f"rule [{rule}]: {n} finding(s), expected 1")
            elif where[rule] != {fname}:
                fails.append(f"rule [{rule}]: found in "
                             f"{sorted(where[rule])}, expected {fname}")
            else:
                ok.append(f"rule [{rule}]: 1 finding in {fname}")
        total = sum(counts.values())
        if total != len(expect):
            fails.append(f"total findings {total} != {len(expect)} "
                         f"(per-rule: {counts})")
        else:
            ok.append(f"total findings == {len(expect)}")

        # exit codes through the real CLI entry point (output swallowed)
        notpy = os.path.join(tmp, 'notes.txt')
        with open(notpy, 'w') as fh:
            fh.write('not python\n')
        buf = io.StringIO()
        with redirect_stdout(buf), redirect_stderr(buf):
            rcs = {
                'planted dir -> 1': (main([planted]), 1),
                'clean dir -> 0': (main([clean]), 0),
                'nonexistent path -> 2':
                    (main([os.path.join(tmp, 'no_such_path')]), 2),
                'non-.py file -> 2': (main([notpy]), 2),
            }
        for name, (got, want) in rcs.items():
            if got != want:
                fails.append(f"exit code {name}: got {got}")
            else:
                ok.append(f"exit code {name}")

        # fixed control: dps set BEFORE the constant -> finding disappears
        fixed = os.path.join(tmp, 'fixed_ok.py')
        with open(fixed, 'w') as fh:
            fh.write("import mpmath as mp\n"
                     "mp.mp.dps = 50\n"
                     "T = mp.mpf(-1) / 3\n")
        got = lint_file(fixed)
        if got:
            fails.append(f"fixed fixture: {len(got)} finding(s), expected 0")
        else:
            ok.append("fixed fixture: 0 findings")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    for line in ok:
        print(f"  PASS {line}", flush=True)
    for line in fails:
        print(f"  FAIL {line}", flush=True)
    print(f"dps_lint selftest: {len(ok)} pass, {len(fails)} fail", flush=True)
    return 1 if fails else 0


def main(argv):
    if '--selftest' in argv:
        return selftest()
    if not argv:
        print(__doc__); return 2
    bad = 0
    for p in argv:
        if p.startswith('-'):
            print(f"dps_lint: unknown option: {p}", file=sys.stderr, flush=True)
            bad += 1
        elif os.path.isfile(p) and not p.endswith('.py'):
            print(f"dps_lint: not a .py file: {p}", file=sys.stderr, flush=True)
            bad += 1
        elif not os.path.isfile(p) and not os.path.isdir(p):
            print(f"dps_lint: no such file or directory: {p}",
                  file=sys.stderr, flush=True)
            bad += 1
    if bad:
        return 2
    files = gather(argv)
    total = 0
    for f in files:
        for (fn, ln, ctx, what, snip) in lint_file(f):
            total += 1
            print(f"{fn}:{ln}: [{ctx}] {what}  ::  {snip.strip()}")
    for (fn, ln, rule, what, snip) in lint_cross(files):
        total += 1
        print(f"{fn}:{ln}: [{rule}] {what}  ::  {snip.strip()}")
    if total:
        print(f"\ndps_lint: {total} finding(s) in {len(files)} file(s)", file=sys.stderr)
        return 1
    print(f"dps_lint: clean ({len(files)} file(s))", file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
