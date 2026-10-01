#!/usr/bin/env python3
"""battery.py — emitall's document-side verb: the paper claims battery.

Canonical core of the FROZEN/NEVER/NEG/ABSTRACT assertion battery that a
shipped paper or report carries: instead of hand-rebuilding the same
comment-stripping / normalization / wrap-proof-join harness per project,
a project writes ONE declarative rows spec and runs this module.

Row taxonomy:
  FROZEN   — pinned headline number/phrase verbatim (whitespace-normalized)
             in its OWNING section file (or the rendered pdf text), optionally
             at an exact count.
  NEVER    — forbidden phrasing appears NOWHERE in rendered prose (comments
             stripped); optional context whitelist for scoped uses.
  NEG      — pins an ABSENCE: the string must land exactly 0.
  ABSTRACT — sentinel phrases present in (or dropped from) the abstract block.

Exit codes (emitall convention): 0 all rows pass; 1 any FAIL; 2 refusal
(stale pdf vs tex — fail closed) or malformed spec.

Invoke:
  python3 tools/emitall/battery.py SPEC.json   [--quiet]
  python3 tools/emitall/battery.py --selftest        # gated fixture battery

Spec shape (JSON):
{
  "label": "example-paper",
  "base": "/abs/path/of/paper_build",       # sources resolved against this
  "tex":  ["sec1.tex", "sec2.tex"],          # prose = comment-stripped tex
  "pdf":  "paper.pdf",                       # rendered text via pdftotext -> d()
  "txt":  "render.txt",                      # OR pre-rendered text (no pdftotext)
  "pdf_fresh_vs": "paper.tex",               # REFUSE (rc 2) if pdf older (fail closed)
  "ligature_fold": false,                    # fold fi/fl/ff ligatures + space-free compare
  "abstract": {"tex": "main.tex"}            # or {"span": ["Abstract", "Contents"]} on rendered text
  "rows": [
    {"class":"FROZEN","name":"catalog-318","file":"sec2.tex","needle":"318 entries"},
    {"class":"FROZEN","name":"pdf-pin","where":"pdf","needle":"318 entries","count":2},
    {"class":"FROZEN","name":"regex-pin","where":"pdf","pattern":"p\\\\s*<\\\\s*10"},
    {"class":"NEVER","name":"no-refuted","pattern":"\\\\brefuted\\\\b",
     "unless_context":["scoped whitelisted phrase"]},
    {"class":"NEG","name":"stale-274","needle":"274"},
    {"class":"ABSTRACT","name":"abs-open","needle":"We grade every catalog entry"},
    {"class":"ABSTRACT","name":"abs-dropped","needle":"retired sentinel clause","absent":true}
  ]
}
"""
import json
import os
import re
import subprocess
import sys

# ---------------------------------------------------------------- primitives


def strip_comments(text):
    """Drop unescaped %-comments per line. Prose rows must never be
    satisfied — or tripped — by commented-out tex."""
    out = []
    for ln in text.split("\n"):
        m = re.search(r"(?<!\\)%", ln)
        if m:
            ln = ln[: m.start()]
        out.append(ln)
    return "\n".join(out)


def norm(s, ligature_fold=False):
    """Whitespace-normalize + curly-quote fold. With ligature_fold, also
    drop fi/fl/ff and ALL spaces (pdftotext drops ligatures, sometimes
    leaving a space)."""
    s = re.sub(r"\s+", " ", s)
    s = s.replace("\u2019", "'").replace("\u201c", '"').replace("\u201d", '"')
    if ligature_fold:
        s = s.replace("fi", "").replace("fl", "").replace("ff", "").replace(" ", "")
    return s


def d(text):
    """Wrap-proof join of rendered (pdftotext) text: form feeds and
    pure-digit page-number lines dropped, hyphen line-wraps rejoined, lines
    joined, whitespace squeezed. A pinned row that straddles a page turn
    still counts as contiguous in reading order."""
    text = text.replace("-\n", "")  # rejoin hyphen wraps
    lines = [ln for ln in text.split("\n") if not re.fullmatch(r"\d{1,4}", ln.strip())]
    text = " ".join(lines).replace("\f", " ")
    return re.sub(r"\s+", " ", text)


def pdf_text(pdf_path, fresh_vs=None):
    """pdftotext -> d(). Fail closed on a stale render: if the pdf is older
    than the tex it claims to render, REFUSE (BatteryRefuse) rather than
    pass against a pre-splice dump ("green != gated": never PASS against a
    stale render).
    """
    if fresh_vs is not None:
        if os.path.getmtime(pdf_path) < os.path.getmtime(fresh_vs):
            raise BatteryRefuse(
                f"{os.path.basename(pdf_path)} older than "
                f"{os.path.basename(fresh_vs)} - rebuild first"
            )
    t = subprocess.run(["pdftotext", pdf_path, "-"], capture_output=True,
                       text=True, check=True).stdout
    return d(t)


class BatteryRefuse(Exception):
    pass


class SpecError(Exception):
    pass


# ---------------------------------------------------------------- the battery


class Battery:
    """Assembles the per-section prose map + rendered text, then runs rows."""

    def __init__(self, spec, base=None):
        self.spec = spec
        self.base = base or spec.get("base") or "."
        self.lig = bool(spec.get("ligature_fold", False))
        self.rows = []  # (class, name, ok, detail)
        # per-section ownership: FROZEN rows anchor to the
        # OWNING file, never to the concatenation.
        self.nprose = {}
        for f in spec.get("tex", []):
            raw = open(os.path.join(self.base, f)).read()
            self.nprose[f] = norm(strip_comments(raw), self.lig)
        self.all_prose = " \n ".join(self.nprose.values())
        # rendered text: pdf (via pdftotext + d()) or pre-rendered txt (via d())
        self.flat = None
        if spec.get("pdf"):
            fv = spec.get("pdf_fresh_vs")
            self.flat = norm(pdf_text(
                os.path.join(self.base, spec["pdf"]),
                os.path.join(self.base, fv) if fv else None), self.lig)
        elif spec.get("txt"):
            tp = os.path.join(self.base, spec["txt"])
            fv = spec.get("pdf_fresh_vs")
            if fv is not None and os.path.getmtime(tp) < os.path.getmtime(
                    os.path.join(self.base, fv)):
                raise BatteryRefuse(
                    f"{spec['txt']} older than {fv} - rebuild first")
            self.flat = norm(d(open(tp).read()), self.lig)
        self.abstract = self._abstract(spec.get("abstract"))

    def _abstract(self, ab):
        if not ab:
            return None
        if "tex" in ab:  # the \begin{abstract} block, comments stripped
            raw = open(os.path.join(self.base, ab["tex"])).read()
            m = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", raw, re.S)
            return norm(strip_comments(m.group(1)), self.lig) if m else ""
        if "span" in ab:  # rendered text between two markers
            if self.flat is None:
                raise SpecError("abstract span needs pdf/txt rendered text")
            a, b = ab["span"]
            m = re.search(re.escape(norm(a, self.lig)) + r"(.*?)" +
                          re.escape(norm(b, self.lig)), self.flat, re.S)
            return m.group(1) if m else ""
        raise SpecError("abstract needs 'tex' or 'span'")

    # -- scopes -------------------------------------------------------------
    def _scope(self, row):
        where = row.get("where")
        if row.get("file"):
            f = row["file"]
            if f not in self.nprose:
                raise SpecError(f"row {row.get('name')}: unknown file {f}")
            return self.nprose[f]
        if where == "pdf":
            if self.flat is None:
                raise SpecError(f"row {row.get('name')}: no rendered text loaded")
            return self.flat
        # default "all": every tex prose + rendered text if present
        return self.all_prose + (" \n " + self.flat if self.flat else "")

    def _hits(self, row, scope):
        if "pattern" in row:
            return [m for m in re.finditer(row["pattern"], scope, re.I if row.get("i") else 0)]
        needle = norm(row["needle"], self.lig)
        return [m for m in re.finditer(re.escape(needle), scope)]

    # -- row runners --------------------------------------------------------
    def run(self):
        for row in self.spec.get("rows", []):
            cls = row.get("class")
            name = row.get("name", "?")
            try:
                if cls == "FROZEN":
                    scope = self._scope(row)
                    hits = self._hits(row, scope)
                    want = row.get("count")
                    ok = (len(hits) == want) if want is not None else bool(hits)
                    det = (f"want {want} got {len(hits)}" if want is not None
                           else ("present" if ok else "MISSING")) + \
                          ": " + (row.get("needle") or row.get("pattern"))[:60]
                elif cls == "NEVER":
                    scope = self._scope(row)
                    wl = [norm(w, self.lig) for w in row.get("unless_context", [])]
                    bad = []
                    for m in self._hits(row, scope):
                        ctx = scope[max(0, m.start() - 120):m.end() + 120]
                        if not any(w in ctx for w in wl):
                            bad.append(ctx[:80])
                    ok = not bad
                    det = bad[0] if bad else ("clean" + (f" ({len(wl)} whitelisted ctx)" if wl else ""))
                elif cls == "NEG":
                    scope = self._scope(row)
                    n = len(self._hits(row, scope))
                    ok = (n == 0)
                    det = f"want 0 got {n}: " + (row.get("needle") or row.get("pattern"))[:60]
                elif cls == "ABSTRACT":
                    if self.abstract is None:
                        raise SpecError("ABSTRACT row without abstract source")
                    if self.abstract == "":
                        ok, det = False, "abstract block not located"
                    else:
                        present = bool(self._hits(row, self.abstract))
                        ok = (not present) if row.get("absent") else present
                        det = ("absent as pinned" if row.get("absent") else
                               "sentinel " + ("present" if ok else "MISSING")) + \
                              ": " + (row.get("needle") or row.get("pattern"))[:50]
                else:
                    raise SpecError(f"row {name}: unknown class {cls!r}")
            except SpecError:
                raise
            except BatteryRefuse:
                raise
            except Exception as e:  # a row error is a FAIL, never a pass
                ok, det = False, f"ROW-ERROR {e}"
            self.rows.append((cls, name, ok, det))
        return self.rows

    def report(self, quiet=False):
        fails = sum(1 for _, _, ok, _ in self.rows if not ok)
        if not quiet:
            for cls, name, ok, det in self.rows:
                print(f"[{'PASS' if ok else 'FAIL'}] {cls:8s} {name}  ({det})")
            print(f"\n{len(self.rows) - fails}/{len(self.rows)} rows pass")
        return 0 if fails == 0 else 1


def run_spec(spec_path, quiet=False):
    try:
        spec = json.load(open(spec_path))
    except Exception as e:
        raise SpecError(f"cannot parse {spec_path}: {e}")
    b = Battery(spec, base=spec.get("base") or os.path.dirname(os.path.abspath(spec_path)))
    b.run()
    return b.report(quiet=quiet)


# ---------------------------------------------------------------- selftest


def selftest():
    """Gated fixture battery: every class must catch its planted failure and
    pass its clean twin; comment-stripping, wrap-proof join, and the stale-
    render refusal are each exercised. Prints PASS/FAIL per gate; rc 0 iff
    all gates hold."""
    import tempfile
    gates = []

    def gate(name, ok):
        gates.append((name, ok))
        print(f"[{'PASS' if ok else 'FAIL'}] selftest: {name}")

    with tempfile.TemporaryDirectory() as td:
        s1 = os.path.join(td, "s1.tex")
        with open(s1, "w") as f:
            f.write("\\begin{abstract}\nWe grade 318 catalog entries under 12 checks.\n"
                    "\\end{abstract}\n"
                    "The catalog holds 318   entries.\n"
                    "% commented: the result is refuted and 274 appears here\n"
                    "Body prose, unverified statements about that estimate.\n")
        # rendered text with a hyphen line-wrap + a page turn INSIDE a pinned
        # phrase + a page-number line
        txt = os.path.join(td, "render.txt")
        with open(txt, "w") as f:
            f.write("Abstract We grade 318 catalog entries under 12 checks. Contents\n"
                    "the combined sc-\nore comes to\n17\n\f41.25 points\n")
        base = {"tex": ["s1.tex"], "txt": "render.txt",
                "abstract": {"tex": "s1.tex"}}

        def run_rows(rows, extra=None):
            spec = dict(base)
            spec.update(extra or {})
            spec["rows"] = rows
            b = Battery(spec, base=td)
            b.run()
            return b

        # 1 clean spec -> rc 0 (whitespace-normalized FROZEN; scoped NEVER;
        #   NEG clean because 274 lives only in a comment; ABSTRACT sentinel)
        b = run_rows([
            {"class": "FROZEN", "name": "catalog", "file": "s1.tex",
             "needle": "The catalog holds 318 entries."},
            {"class": "FROZEN", "name": "wrap-proof", "where": "pdf",
             "needle": "the combined score comes to 41.25 points"},
            {"class": "NEVER", "name": "scoped-unverified", "pattern": r"unverified\w*",
             "unless_context": ["statements about that estimate"]},
            {"class": "NEG", "name": "stale-274", "needle": "274"},
            {"class": "ABSTRACT", "name": "abs-open", "needle": "We grade 318 catalog entries"},
        ])
        gate("clean spec rc 0", b.report(quiet=True) == 0)
        # 2 planted FROZEN break -> rc 1
        b = run_rows([{"class": "FROZEN", "name": "broken", "file": "s1.tex",
                       "needle": "The catalog holds 999 entries."}])
        gate("FROZEN catches planted break", b.report(quiet=True) == 1)
        # 3 NEVER catches an unwhitelisted hit -> rc 1
        b = run_rows([{"class": "NEVER", "name": "unverified-nowl",
                       "pattern": r"unverified\w*"}])
        gate("NEVER catches unscoped hit", b.report(quiet=True) == 1)
        # 4 NEVER/NEG never satisfied — and never tripped — by comments:
        #   "refuted" sits ONLY in a % comment, so NEVER stays clean...
        b = run_rows([{"class": "NEVER", "name": "no-refuted", "pattern": r"\brefuted\b"}])
        ok_a = b.report(quiet=True) == 0
        #   ...and a FROZEN row must NOT be satisfied by commented text.
        b = run_rows([{"class": "FROZEN", "name": "comment-not-prose", "file": "s1.tex",
                       "needle": "the result is refuted"}])
        gate("comment stripping (both directions)", ok_a and b.report(quiet=True) == 1)
        # 5 NEG catches a real occurrence -> rc 1
        b = run_rows([{"class": "NEG", "name": "neg-hit", "needle": "318"}])
        gate("NEG catches presence", b.report(quiet=True) == 1)
        # 6 ABSTRACT missing sentinel -> rc 1; absent-row (drop pin) passes
        b = run_rows([
            {"class": "ABSTRACT", "name": "abs-missing", "needle": "nonexistent sentinel"},
            {"class": "ABSTRACT", "name": "abs-drop", "needle": "old dropped clause",
             "absent": True}])
        gate("ABSTRACT sentinel + drop pin",
             b.report(quiet=True) == 1 and b.rows[0][2] is False and b.rows[1][2] is True)
        # 7 abstract via rendered span
        b = run_rows([{"class": "ABSTRACT", "name": "abs-span",
                       "needle": "We grade 318 catalog entries"}],
                     extra={"abstract": {"span": ["Abstract", "Contents"]}})
        gate("ABSTRACT from rendered span", b.report(quiet=True) == 0)
        # 8 stale-render refusal (fail closed, rc 2 semantics)
        os.utime(txt, (1, 1))  # make the render ancient
        try:
            run_rows([], extra={"pdf_fresh_vs": "s1.tex"})
            gate("stale render refused", False)
        except BatteryRefuse:
            gate("stale render refused", True)

    fails = [n for n, ok in gates if not ok]
    print(f"\nselftest: {len(gates) - len(fails)}/{len(gates)} gates green")
    return 0 if not fails else 1


# ---------------------------------------------------------------- CLI


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    quiet = "--quiet" in argv
    argv = [a for a in argv if a != "--quiet"]
    if argv and argv[0] == "--selftest":
        return selftest()
    if len(argv) != 1:
        print(__doc__.split("\n\n")[0])
        print("usage: battery.py SPEC.json [--quiet] | battery.py --selftest")
        return 2
    try:
        return run_spec(argv[0], quiet=quiet)
    except SpecError as e:
        print(f"SPEC ERROR: {e}", file=sys.stderr)
        return 2
    except BatteryRefuse as e:
        print(f"BATTERY REFUSE: {e}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
