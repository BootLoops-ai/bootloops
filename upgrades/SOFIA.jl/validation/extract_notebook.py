#!/usr/bin/env python3
"""Extract (diagram, recorded singularities) pairs from SOFIA_examples.nb.

Strategy: every `candidateSingularities = SOFIA[diag...]` input cell is
followed by an Output cell whose BoxData encodes the singularity list
(possibly via TableForm). We convert box forms to linear Wolfram InputForm
text and emit JSON: [{"pos", "diag", "output", "tableform"}...].
The nearest preceding `diag... = {...}` assignment supplies the diagram.
"""
import re
import json
import sys

SRC = "reference/SOFIA_examples.nb"
text = open(SRC, encoding="utf-8", errors="ignore").read()

# ---------- generic box-expression parser ----------
# Grammar (textual WL in .nb files): Head[arg, ...], "string", {a, b}, atom

class P:
    def __init__(self, s, i=0):
        self.s = s
        self.i = i

    def ws(self):
        while self.i < len(self.s) and self.s[self.i] in " \t\r\n":
            self.i += 1

    def parse(self):
        e = self.primary()
        self.ws()
        # infix rules / settings in option position: ignore semantics
        while self.i + 1 < len(self.s) and self.s[self.i:self.i+2] in ("->", ":>"):
            self.i += 2
            rhs = self.primary()
            e = ("rule", e, rhs)
            self.ws()
        return e

    def primary(self):
        self.ws()
        c = self.s[self.i]
        if c == '"':
            return self.string()
        if c == "{":
            return self.list_()
        m = re.match(r"[-+]?[0-9][0-9.]*(`+[0-9.]*)?(\*\^[-+0-9]+)?", self.s[self.i:])
        if m:
            self.i += len(m.group(0))
            return ("sym", m.group(0))
        return self.symbol_or_call()

    def string(self):
        assert self.s[self.i] == '"'
        self.i += 1
        out = []
        while True:
            c = self.s[self.i]
            if c == "\\":
                out.append(self.s[self.i:self.i+2])
                self.i += 2
            elif c == '"':
                self.i += 1
                break
            else:
                out.append(c)
                self.i += 1
        return ("str", "".join(out))

    def list_(self):
        assert self.s[self.i] == "{"
        self.i += 1
        items = []
        self.ws()
        if self.s[self.i] == "}":
            self.i += 1
            return ("list", items)
        while True:
            items.append(self.parse())
            self.ws()
            if self.s[self.i] == ",":
                self.i += 1
                continue
            if self.s[self.i] == "}":
                self.i += 1
                return ("list", items)
            raise ValueError(f"bad list at {self.i}: {self.s[self.i:self.i+40]!r}")

    def symbol_or_call(self):
        m = re.match(r"[\w$`]+", self.s[self.i:])
        if not m:
            # operators like -> inside option rules; treat single char as atom
            ch = self.s[self.i]
            self.i += 1
            return ("sym", ch)
        name = m.group(0)
        self.i += len(name)
        self.ws()
        if self.i < len(self.s) and self.s[self.i] == "[":
            self.i += 1
            args = []
            self.ws()
            if self.s[self.i] == "]":
                self.i += 1
                return ("call", name, args)
            while True:
                args.append(self.parse())
                self.ws()
                if self.s[self.i] == ",":
                    self.i += 1
                    continue
                if self.s[self.i] == "]":
                    self.i += 1
                    return ("call", name, args)
                raise ValueError(f"bad call {name} at {self.i}: {self.s[self.i:self.i+40]!r}")
        return ("sym", name)


IGNORED_STR = {
    "\\[IndentingNewLine]", "\\[InvisibleSpace]", "\\[NoBreak]", "\\[NonBreakingSpace]",
}

def boxes_to_text(node):
    """Convert a parsed box expression to linear WL input text."""
    kind = node[0]
    if kind == "str":
        s = node[1]
        if s in IGNORED_STR:
            return ""
        s = s.replace("\\[IndentingNewLine]", " ").replace("\\n", " ")
        s = s.replace("\\[InvisibleSpace]", "")
        return s
    if kind == "sym":
        return ""          # bare symbols in box context are formatting atoms
    if kind == "rule":
        return ""          # option rules are formatting, not content
    if kind == "list":
        return "".join(boxes_to_text(x) for x in node[1])
    # calls
    name, args = node[1], node[2]
    if name == "RowBox":
        return "".join(boxes_to_text(x) for x in args[0][1])
    if name == "SuperscriptBox":
        return f"({boxes_to_text(args[0])})^({boxes_to_text(args[1])})"
    if name == "SubscriptBox":
        return f"{boxes_to_text(args[0])}{boxes_to_text(args[1])}"
    if name == "FractionBox":
        return f"(({boxes_to_text(args[0])})/({boxes_to_text(args[1])}))"
    if name == "SqrtBox":
        return f"Sqrt[{boxes_to_text(args[0])}]"
    if name in ("TagBox", "StyleBox", "InterpretationBox", "FormBox",
                "ItemBox", "PaneBox", "AdjustmentBox", "TooltipBox"):
        return boxes_to_text(args[0])
    if name == "GridBox":
        rows = args[0]
        cells = []
        for row in rows[1]:
            t = ",".join(boxes_to_text(c) for c in row[1])
            cells.append(t)
        return "{" + ",".join(cells) + "}"
    if name == "Cell" or name == "BoxData":
        return boxes_to_text(args[0])
    return "".join(boxes_to_text(a) for a in args)


def cell_region(start):
    """Return (text, end) of the bracket-balanced Cell[...] starting at start."""
    i = text.index("[", start)
    depth = 0
    j = i
    instr = False
    while j < len(text):
        c = text[j]
        if instr:
            if c == "\\":
                j += 2
                continue
            if c == '"':
                instr = False
        else:
            if c == '"':
                instr = True
            elif c == "[":
                depth += 1
            elif c == "]":
                depth -= 1
                if depth == 0:
                    return text[start:j+1], j + 1
        j += 1
    raise ValueError("unbalanced cell")


# ---------- section structure for robust pairing ----------
sections = []
for m in re.finditer(r'Cell\[(?:TextData\[)?"((?:[^"\\]|\\.)*)"(?:\])?,\s*"(Section|Subsection|Subsubsection|Title|Chapter)"', text):
    if m.group(2) in ("Section", "Subsection"):
        sections.append((m.start(), m.group(1)[:60]))
sections.sort()

def section_of(pos):
    cur = (0, "preamble")
    for s in sections:
        if s[0] < pos:
            cur = s
        else:
            break
    return cur

# ---------- diagram corpus (built here, with substitution-rule capture) ----------
diags = []
for m in re.finditer(r'RowBox\[\{"(diag\w*)", "="', text):
    start = m.end()
    chunk = text[start:start+9000]
    toks = re.findall(r'"((?:[^"\\]|\\.)*)"', chunk)
    expr, depth, started = '', 0, False
    rules = ''
    it = iter(range(len(toks)))
    k = 0
    while k < len(toks):
        tok = toks[k]; k += 1
        if tok.startswith('\\['):
            continue
        expr += tok
        depth += tok.count('{') - tok.count('}')
        if '{' in tok:
            started = True
        if started and depth == 0:
            break
    if not (started and expr.startswith('{')):
        continue
    # capture a following rule chain: /. or //. up to ; or a non-rule token
    while k < len(toks):
        tok = toks[k]
        if tok.startswith('\\['):
            k += 1; continue
        if tok in ('/.', '//.', ':>', '->', '|', '_', '__') or re.fullmatch(r'[\w|_]+', tok or '') :
            rules += tok
            k += 1
        else:
            break
    rules = rules.split(';')[0]
    if not re.match(r'^(/\.|//\.)', rules):
        rules = ''
    diags.append((m.start(), m.group(1), ' '.join(expr.split()), rules))
diags.sort()
with open("validation/notebook_diagrams.txt", "w") as f:
    for pos, name, e, r in diags:
        f.write(f"{name}\t{pos}\t{e}\t{r}\n")

cases = []
for m in re.finditer(r'RowBox\[\{"(candidateSingularities\w*|sings?)", "="', text):
    callpos = m.start()
    # is this a plain singularity call (not FindLetters)?
    head = text[callpos:callpos+1200]
    if "FindLetters" in head and "True" in head.split("FindLetters")[1][:40]:
        continue
    # nearest preceding diagram assignment WITHIN the same section
    secstart, sectitle = section_of(callpos)
    prior = [d for d in diags if secstart <= d[0] < callpos]
    if not prior:
        prior = [d for d in diags if d[0] < callpos]   # fallback: any prior
        sectitle += " [cross-section pairing]"
    if not prior:
        continue
    dpos, dname, dexpr, drules = prior[-1]
    # find the next Output cell after the call
    out_m = re.compile(r'Cell\[BoxData\[').search(text, callpos)
    found = None
    search_from = callpos
    for _ in range(8):
        out_m = re.compile(r'Cell\[BoxData\[').search(text, search_from)
        if out_m is None:
            break
        cell, end = cell_region(out_m.start())
        search_from = end
        if '"Output"' in cell:
            found = cell
            break
    if found is None:
        continue
    try:
        node = P(found).parse()
        lin = boxes_to_text(node)
    except Exception as e:
        cases.append({"pos": callpos, "diag": dexpr, "error": str(e)})
        continue
    # linearize the input cell containing the call (for options/substitutions)
    cellstart = text.rfind("Cell[BoxData[", 0, callpos)
    try:
        incell, _ = cell_region(cellstart)
        intext = boxes_to_text(P(incell).parse())
    except Exception:
        intext = ""
    cases.append({"pos": callpos, "diag_name": dname, "diag": dexpr,
                  "rules": drules, "section": sectitle, "input": intext,
                  "output": lin})

json.dump(cases, open("validation/notebook_cases.json", "w"), indent=1)
with open("validation/notebook_cases.tsv", "w") as f:
    for c in cases:
        clean = lambda s: ' '.join(str(s).split()).replace('\t', ' ')
        f.write("\t".join([str(c.get("pos", "")), clean(c.get("diag", "")),
                           clean(c.get("input", c.get("error", ""))),
                           clean(c.get("output", "")),
                           clean(c.get("section", "")),
                           clean(c.get("rules", ""))]) + "\n")
print(f"extracted {len(cases)} cases")
for c in cases[:6]:
    print("DIAG:", c["diag"][:90])
    print("OUT :", c.get("output", c.get("error", ""))[:160])
    print()
