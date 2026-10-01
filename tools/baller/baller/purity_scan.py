"""Mechanical purity scanner for a numerics purity bar.

Scans python (AST walk + regex) and C (regex) sources for finite-precision
literals, classifying each hit:

  PURE-INT-RAT   : no float involved (ints, mpf(p)/mpf(q), sp.Rational) —
                   not reported as hits.
  PURE-DPSREL    : dps-relative bound  10**(-dps+K), 10**(-mp.dps+K), etc.
  PURE-CONV      : float literal that can only affect CONVERGENCE/gating,
                   not the defined value: comparison thresholds, atol/rtol/
                   tol/bar/eps keywords, maxterms, error prints, plot code.
  RAT-STRING     : exact rational-valued decimal string ('0.5', '0.25') fed
                   to mpf/Float/Rational — exact, PURE.
  DEC-STRING     : long decimal string literal (>=8 sig digits) — banked
                   value string.  VIOLATION if on a value path (adjudicate).
  FLOAT-EXPR     : float literal in an arithmetic/value expression —
                   VIOLATION candidate (adjudicate).
  MPF-OF-FLOAT   : mp.mpf(<float literal>) / mpf(0.28...) — the classic
                   disease signature.  VIOLATION candidate.

Output: JSON hits table per file + a summary; every FLOAT-EXPR / DEC-STRING
/ MPF-OF-FLOAT hit carries its source line for manual adjudication.
Exit code 0 always (report, don't gate — adjudication is human).
"""
import ast
import json
import os
import re
import sys
from fractions import Fraction

# keywords whose keyword-argument float values are convergence-only
CONV_KEYWORDS = {
    "atol", "rtol", "tol", "eps", "abs_err", "maxterms", "maxdegree",
    "bar", "thresh", "threshold", "cutoff", "window", "pad", "margin",
    "alpha", "figsize", "dpi", "lw", "ms", "s", "fontsize", "timeout",
}
# names that mark a line as convergence/gate context
CONV_NAME_RE = re.compile(
    r"(atol|rtol|\btol\b|bar|thresh|gate|assert|abs_err|err_est|"
    r"maxterms|maxdegree|timeout|fuse|budget|plot|plt\.|figsize|"
    r"legend|xlim|ylim|title|label)", re.I)

DEC_STRING_RE = re.compile(r"^-?\d*\.\d{8,}([eE][+-]?\d+)?$")
SHORT_DEC_RE = re.compile(r"^-?\d*\.?\d+([eE][+-]?\d+)?$")


def _is_exact_rational_string(s):
    """decimal string whose value is an exact dyadic-free rational with a
    SHORT literal (0.5, 0.25, 1e-6 as 10^-6 etc.) — exact iff it round-trips
    through Fraction without precision loss (always true for decimal
    strings) AND is short enough to be a coefficient, not a banked value."""
    try:
        Fraction(s)
    except (ValueError, ZeroDivisionError):
        return False
    digs = re.sub(r"[^0-9]", "", s.split("e")[0].split("E")[0]).lstrip("0")
    return len(digs) < 8


class Hit:
    __slots__ = ("path", "line", "col", "kind", "literal", "src")

    def __init__(self, path, line, col, kind, literal, src):
        self.path, self.line, self.col = path, line, col
        self.kind, self.literal, self.src = kind, literal, src

    def to_dict(self):
        return {"file": self.path, "line": self.line, "col": self.col,
                "kind": self.kind, "literal": self.literal,
                "src": self.src.strip()[:180]}


def scan_python(path):
    hits = []
    with open(path, "rb") as f:
        src_bytes = f.read()
    try:
        src = src_bytes.decode("utf-8")
    except UnicodeDecodeError:
        src = src_bytes.decode("latin-1")
    lines = src.splitlines()

    def L(n):
        return lines[n - 1] if 0 < n <= len(lines) else ""

    try:
        tree = ast.parse(src)
    except SyntaxError as e:
        hits.append(Hit(path, e.lineno or 0, 0, "PARSE-ERROR", str(e), ""))
        return hits

    # map: (lineno,col) of float constants inside dps-relative power exprs
    dpsrel = set()
    conv_kw = set()
    mpf_float = set()
    call_ctx = {}

    class V(ast.NodeVisitor):
        def visit_BinOp(self, node):
            # 10 ** (-dps + K) / 10.0 ** ( ... dps ... )
            if isinstance(node.op, ast.Pow):
                seg = ast.get_source_segment(src, node) or ""
                if re.search(r"\bdps\b|\bprec\b|\bqdps\b|\bwp\b", seg):
                    for sub in ast.walk(node):
                        if (isinstance(sub, ast.Constant)
                                and isinstance(sub.value, float)):
                            dpsrel.add((sub.lineno, sub.col_offset))
            self.generic_visit(node)

        def visit_Call(self, node):
            fn = ast.get_source_segment(src, node.func) or ""
            for kw in node.keywords:
                if kw.arg in CONV_KEYWORDS:
                    for sub in ast.walk(kw.value):
                        if (isinstance(sub, ast.Constant)
                                and isinstance(sub.value, float)):
                            conv_kw.add((sub.lineno, sub.col_offset))
            if re.search(r"(^|\.)mpf$|(^|\.)mpmathify$|(^|\.)Float$", fn):
                for a in node.args[:1]:
                    if (isinstance(a, ast.Constant)
                            and isinstance(a.value, float)):
                        mpf_float.add((a.lineno, a.col_offset))
                    if (isinstance(a, ast.Constant)
                            and isinstance(a.value, str)):
                        call_ctx[(a.lineno, a.col_offset)] = (fn, a.value)
            self.generic_visit(node)

    V().visit(tree)

    for node in ast.walk(tree):
        if not isinstance(node, ast.Constant):
            continue
        ln, co = node.lineno, node.col_offset
        v = node.value
        if isinstance(v, float):
            key = (ln, co)
            src_line = L(ln)
            if key in mpf_float:
                hits.append(Hit(path, ln, co, "MPF-OF-FLOAT", repr(v),
                                src_line))
            elif key in dpsrel:
                hits.append(Hit(path, ln, co, "PURE-DPSREL", repr(v),
                                src_line))
            elif key in conv_kw or CONV_NAME_RE.search(src_line):
                hits.append(Hit(path, ln, co, "PURE-CONV", repr(v),
                                src_line))
            else:
                hits.append(Hit(path, ln, co, "FLOAT-EXPR", repr(v),
                                src_line))
        elif isinstance(v, str):
            sv = v.strip()
            if DEC_STRING_RE.match(sv):
                hits.append(Hit(path, ln, co, "DEC-STRING", sv, L(ln)))
            elif ((ln, co) in call_ctx and SHORT_DEC_RE.match(sv)
                  and "." in sv and not _is_exact_rational_string(sv)):
                hits.append(Hit(path, ln, co, "DEC-STRING", sv, L(ln)))
    return hits


C_FLOAT_RE = re.compile(
    r"(?<![\w.])(\d+\.\d*([eE][+-]?\d+)?|\.\d+([eE][+-]?\d+)?"
    r"|\d+[eE][+-]?\d+)(?![\w.])")


def scan_c(path):
    hits = []
    with open(path, "r", errors="replace") as f:
        for ln, line in enumerate(f, 1):
            stripped = line.split("//")[0]
            if stripped.lstrip().startswith(("*", "/*")):
                continue
            for m in C_FLOAT_RE.finditer(stripped):
                lit = m.group(0)
                if _is_exact_rational_string(lit):
                    kind = ("PURE-CONV" if CONV_NAME_RE.search(line)
                            else "FLOAT-EXPR")
                    # short exact-rational C literals still float-typed:
                    # exact only if double covers them — flag short ones
                    # as FLOAT-EXPR-SHORT for adjudication
                    kind = kind if kind == "PURE-CONV" else "FLOAT-EXPR-SHORT"
                else:
                    kind = "DEC-STRING"
                hits.append(Hit(path, ln, m.start(), kind, lit, line))
    return hits


def main(list_file, out_json):
    files = [l.strip() for l in open(list_file)
             if l.strip() and not l.startswith("#")]
    all_hits, table = [], {}
    for p in files:
        if not os.path.isfile(p):
            table[p] = {"status": "MISSING", "hits": []}
            continue
        hits = (scan_c(p) if p.endswith((".c", ".h"))
                else scan_python(p))
        all_hits += hits
        kinds = {}
        for h in hits:
            kinds[h.kind] = kinds.get(h.kind, 0) + 1
        table[p] = {"status": "scanned", "kind_counts": kinds,
                    "hits": [h.to_dict() for h in hits]}
    summary = {}
    for h in all_hits:
        summary[h.kind] = summary.get(h.kind, 0) + 1
    out = {"n_files": len(files), "summary_kind_counts": summary,
           "files": table}
    with open(out_json, "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(summary, indent=1))
    print(f"wrote {out_json}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
