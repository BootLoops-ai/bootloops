#!/usr/bin/env python3
"""form_route — momentum-routing glue for FORM 5 diagram-generator output.

FORM 5's integrated generator (diagrams_) emits diagrams as terms

    coeff * topo_(t) * node_(id,coupling,field(mom),...) [* edge_(id,field(mom),v1,v2) ...]

where EVERY internal edge carries an independent momentum label and NO
momentum conservation is imposed (port battery C finding). This tool is the missing
spanning-tree glue between generator output and integrand programs.

INPUT (file path or string; unwrapped via tools/formglue/form_io.read_form_output):
  - full graph output (node_ form, with or without WithEdges_),
  - the census form (TopologiesOnly_ + WithEdges_, bare momenta),
  - vacuum graphs (no external legs) and self-loops (tadpole edges) included.

CONVENTIONS (pinned on real FORM 5.0.1 output,
tests/fixtures/gen/ + the FORM 5 manual "Diagram generation"):
  - node_ field arguments are functions of the momenta flowing INTO the node;
    momentum conservation at a node <=> the printed arguments sum to zero.
  - external legs are node_ functions with a single field and coupling 1.
  - every momentum label appears exactly twice, once with each sign (a
    self-loop carries both signs on the same node).
  - edge_(id, field(mom), v1, v2): mom flows v1 -> v2 (printed -mom at v1,
    +mom at v2 in the node_ list).

ROUTING: per diagram, build the multigraph on internal vertices, pick a
deterministic BFS spanning tree, keep each non-tree (chord) edge's own
generator label as an independent loop momentum (L = E_int - V_int + 1 of
them; self-loops are always chords), and solve every tree-edge momentum from
vertex conservation by leaf-stripping — exact integer/Fraction linear
algebra, no floats. External momenta keep the generator's leg labels (or are
renamed/ordered via externals=[...], leg-node-id order).

MANDATORY SELF-CHECKS (run on EVERY invocation, rc-coded — RouteError and
nonzero exit on any failure):
  1. momentum conservation at EVERY internal vertex after assignment
     (residual reduced modulo overall external conservation sum p_i = 0);
  2. Euler count E_int = V_int + L - 1 (connected graph);
  3. every external appears exactly once (one leg node, one internal vertex);
  4. loop-momentum count == first Betti number, recomputed INDEPENDENTLY via
     union-find (E_int - V_int + n_components, n_components must be 1);
  5. every chord carries exactly its own loop label, unit coefficient.
The --mutate-edge control hook corrupts one routed momentum AFTER solving
and BEFORE the checks: it exists so mutation controls can prove check (1)
is live (a gate that cannot fail is void). Never use it for physics runs.

OUTPUT:
  (a) routing table (edge label -> momentum expression),
  (b) a FORM substitution block: one `#procedure route<i>()` of `id` (and
      `#define ROUTED<i>LOOPS "..."`) statements per diagram, ready to
      #include into an integrand program (RHS contain only externals and
      chord labels, so statement order is immaterial),
  (c) optional JSON with the full per-diagram routing structure.

Acceptance tests: tools/formglue/tests/test_form_route.py (demo gg->qqbar box
equivalence vs hand-wiring, 30/30 on-shell set, 74/74 QCD 3-loop vacuum
set (L=3, self-loops), census form, mutation controls).
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from fractions import Fraction

# Package-relative import — never routed through the
# flat tools/ shims. A dirname(__file__) sys.path.insert is deliberately
# absent: inside the package it would alias formglue modules as top-level
# duplicates. Script mode (python3 tools/formglue/form_route.py) anchors
# the package PARENT (tools/) instead.
try:
    from .form_io import read_form_output  # noqa: E402
except ImportError:  # direct script execution: no parent package
    _PKG_PARENT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if _PKG_PARENT not in sys.path:
        sys.path.insert(0, _PKG_PARENT)
    from formglue.form_io import read_form_output  # noqa: E402

__all__ = ["RouteError", "route_text", "route_file", "RoutedDiagram",
           "mom_str", "form_block", "routing_table"]

_SYM_RE = re.compile(r"^([+-]?)([A-Za-z][A-Za-z0-9]*)$")
_FIELD_RE = re.compile(r"^([A-Za-z][A-Za-z0-9]*)\(([+-]?[A-Za-z][A-Za-z0-9]*)\)$")
_COEF_RE = re.compile(r"^\d+(/\d+)?$")
_NAME_EQ_RE = re.compile(r"^[A-Za-z][A-Za-z0-9]*=")


class RouteError(RuntimeError):
    """Any parse/graph/self-check failure. Always fatal, always loud."""


# --------------------------------------------------------------- tokenization

def _split_top(s: str, seps: str):
    """Split s at depth-0 (w.r.t. parentheses) occurrences of chars in seps.
    Separators are dropped; for +/- splitting use _split_terms instead."""
    parts, depth, cur = [], 0, []
    for ch in s:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth < 0:
                raise RouteError(f"unbalanced ')' in {s!r}")
        if depth == 0 and ch in seps:
            parts.append("".join(cur))
            cur = []
        else:
            cur.append(ch)
    if depth != 0:
        raise RouteError(f"unbalanced '(' in {s!r}")
    parts.append("".join(cur))
    return parts


def _split_terms(s: str):
    """Split an expression string into signed top-level terms."""
    terms, depth, cur = [], 0, []
    for ch in s:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth < 0:
                raise RouteError(f"unbalanced ')' in expression")
        if depth == 0 and ch in "+-" and cur and cur[-1] not in "+-^*/(":
            terms.append("".join(cur))
            cur = [ch]
        else:
            cur.append(ch)
    if depth != 0:
        raise RouteError("unbalanced '(' in expression")
    if cur:
        terms.append("".join(cur))
    return [t for t in terms if t.strip()]


# --------------------------------------------------------------------- momenta

def _mvec(items=()):
    v = {}
    for sym, c in items:
        c = Fraction(c)
        if c:
            v[sym] = v.get(sym, Fraction(0)) + c
            if not v[sym]:
                del v[sym]
    return v


def _madd(a, b, cb=Fraction(1)):
    out = dict(a)
    for s, c in b.items():
        out[s] = out.get(s, Fraction(0)) + cb * c
        if not out[s]:
            del out[s]
    return out


def mom_str(v, order=()):
    """Deterministic FORM-ready string for a momentum vector."""
    if not v:
        return "0"
    keys = [s for s in order if s in v] + sorted(s for s in v if s not in order)
    out = []
    for s in keys:
        c = v[s]
        sign = "-" if c < 0 else "+"
        mag = abs(c)
        out.append(f"{sign}{s}" if mag == 1 else f"{sign}{mag}*{s}")
    txt = "".join(out)
    return txt[1:] if txt.startswith("+") else txt


# ---------------------------------------------------------------- data model

class _Edge:
    __slots__ = ("key", "eid", "field", "sym", "v1", "v2", "kind", "mom")

    def __init__(self, key, eid, field, sym, v1, v2, kind):
        self.key, self.eid, self.field, self.sym = key, eid, field, sym
        self.v1, self.v2, self.kind = v1, v2, kind  # kind: ext | tree | chord
        self.mom = None  # momentum vector flowing v1 -> v2


class RoutedDiagram:
    """One routed diagram: graph metadata + exact routing + check record."""

    def __init__(self, idx, topo, coeff):
        self.idx, self.topo, self.coeff = idx, topo, coeff
        self.nodes = {}          # node id -> [(field, sign, sym), ...]
        self.leg_nodes = []      # ascending node ids
        self.externals = []      # external syms in leg-node order
        self.edges = {}          # sym -> _Edge
        self.loop_syms = []      # chord labels (= loop momenta) in edge order
        self.V_int = 0
        self.E_int = 0
        self.L = 0
        self.checks = []         # (name, "ok") — failures raise instead

    def routing(self):
        """dict label -> momentum vector (flowing v1->v2 for that edge)."""
        return {e.sym: dict(e.mom) for e in self.edges.values()}

    def sym_order(self):
        return list(self.externals) + list(self.loop_syms)

    def to_json(self):
        order = self.sym_order()
        return {
            "diagram": self.idx, "topo": self.topo, "coeff": self.coeff,
            "externals": list(self.externals),
            "loop_momenta": list(self.loop_syms),
            "V_int": self.V_int, "E_int": self.E_int, "L": self.L,
            "edges": [
                {"label": e.sym, "id": e.eid, "field": e.field,
                 "v1": e.v1, "v2": e.v2, "kind": e.kind,
                 "momentum": mom_str(e.mom, order)}
                for e in self._edge_list()],
            "checks": [n for n, _ in self.checks],
        }

    def _edge_list(self):
        return sorted(self.edges.values(),
                      key=lambda e: (e.eid if e.eid is not None else 0, e.sym))


# ------------------------------------------------------------------- parsing

def _parse_fieldspec(arg, where):
    m = _FIELD_RE.match(arg)
    if m:
        field, inner = m.group(1), m.group(2)
        sm = _SYM_RE.match(inner)
        return field, (-1 if sm.group(1) == "-" else +1), sm.group(2)
    m = _SYM_RE.match(arg)
    if m:
        return None, (-1 if m.group(1) == "-" else +1), m.group(2)
    raise RouteError(f"{where}: cannot parse field/momentum spec {arg!r} "
                     "(only single signed labels are generator output; "
                     "routed expressions are not valid router input)")


def _parse_term(term, idx):
    term = term.strip()
    sign = "+"
    if term[:1] in "+-":
        sign, term = term[0], term[1:]
    coeff_parts, topo = [], None
    dia = None
    factors = _split_top(term, "*")
    nodes_raw, edges_raw = [], []
    for f in factors:
        f = f.strip()
        if not f:
            raise RouteError(f"diagram {idx}: empty factor in term")
        if _COEF_RE.match(f):
            coeff_parts.append(f)
        elif f.startswith("topo_(") and f.endswith(")"):
            topo = int(f[len("topo_("):-1])
        elif f.startswith("node_(") and f.endswith(")"):
            nodes_raw.append(_split_top(f[len("node_("):-1], ","))
        elif f.startswith("edge_(") and f.endswith(")"):
            edges_raw.append(_split_top(f[len("edge_("):-1], ","))
        elif (f.startswith("block_(") or f.startswith("onepi_(")) and f.endswith(")"):
            pass  # tagging functions: tolerated, not needed for routing
        else:
            raise RouteError(f"diagram {idx}: unrecognized factor {f!r}")
    if not nodes_raw:
        raise RouteError(f"diagram {idx}: no node_ functions "
                         "(WithoutNodes_ output is not supported: node "
                         "occurrences are the routing ground truth)")
    coeff = sign + ("*".join(coeff_parts) if coeff_parts else "1")
    dia = RoutedDiagram(idx, topo, coeff)

    for args in nodes_raw:
        if len(args) < 3:
            raise RouteError(f"diagram {idx}: node_ with <3 args: {args}")
        nid = int(args[0])
        if nid in dia.nodes:
            raise RouteError(f"diagram {idx}: duplicate node id {nid}")
        fields = [_parse_fieldspec(a, f"diagram {idx} node {nid}")
                  for a in args[2:]]
        dia.nodes[nid] = fields
    dia.leg_nodes = sorted(n for n, ff in dia.nodes.items() if len(ff) == 1)

    # momentum-label occurrence census from node_ functions (ground truth)
    occ = {}
    for nid, fields in dia.nodes.items():
        for _, s, sym in fields:
            occ.setdefault(sym, []).append((nid, s))
    for sym, oc in occ.items():
        if len(oc) != 2 or oc[0][1] * oc[1][1] != -1:
            raise RouteError(
                f"diagram {idx}: label {sym} must appear exactly twice with "
                f"opposite signs, got {oc} (corrupt/edited input?)")

    def _mk_edge(eid, field, sym, v1, v2):
        if sym in dia.edges:
            raise RouteError(f"diagram {idx}: duplicate edge label {sym}")
        kind = "ext" if (v1 in dia.leg_nodes or v2 in dia.leg_nodes) else "int"
        dia.edges[sym] = _Edge(sym, eid, field, sym, v1, v2, kind)

    if edges_raw:
        for args in edges_raw:
            if len(args) != 4:
                raise RouteError(f"diagram {idx}: edge_ with !=4 args: {args}")
            eid = int(args[0])
            field, s, sym = _parse_fieldspec(args[1], f"diagram {idx} edge {args[0]}")
            if s != +1:
                raise RouteError(f"diagram {idx}: edge_ {eid} momentum "
                                 f"printed with a sign: {args[1]!r}")
            v1, v2 = int(args[2]), int(args[3])
            # cross-check vs node occurrences: mom flows v1 -> v2
            want = sorted([(v1, -1), (v2, +1)])
            if sorted(occ[sym]) != want:
                raise RouteError(
                    f"diagram {idx}: edge_({eid},{sym},{v1},{v2}) is "
                    f"inconsistent with node_ occurrences {occ[sym]} "
                    "(corrupt/edited input?)")
            _mk_edge(eid, field, sym, v1, v2)
        if set(dia.edges) != set(occ):
            raise RouteError(f"diagram {idx}: edge_ labels "
                             f"{sorted(dia.edges)} != node labels {sorted(occ)}")
    else:
        for sym in sorted(occ):
            oc = occ[sym]
            v_minus = [n for n, s in oc if s == -1][0]
            v_plus = [n for n, s in oc if s == +1][0]
            _mk_edge(None, None, sym, v_minus, v_plus)
    return dia


# ------------------------------------------------------- routing (per diagram)

def _route_diagram(dia, externals=None, loops=None, mutate_edge=None):
    idx = dia.idx
    # ---- externals: leg-node order; optional user ordering / renaming
    det_ext = []
    for leg in dia.leg_nodes:
        (field, s, sym) = dia.nodes[leg][0]
        det_ext.append(sym)
    if len(set(det_ext)) != len(det_ext):
        raise RouteError(f"diagram {idx}: repeated external label in {det_ext}")
    rename = {}
    if externals is not None:
        externals = list(externals)
        if len(externals) != len(det_ext):
            raise RouteError(
                f"diagram {idx}: {len(externals)} external names supplied "
                f"for {len(det_ext)} legs {det_ext}")
        if set(externals) == set(det_ext):
            dia.externals = externals          # user ORDER, generator names
        else:
            all_syms = set(dia.edges)
            for new in externals:
                if new in all_syms and new not in det_ext:
                    raise RouteError(
                        f"diagram {idx}: external rename {new!r} collides "
                        f"with an internal label")
            rename = dict(zip(det_ext, externals))
            dia.externals = externals
    else:
        dia.externals = det_ext

    ext_edges = [e for e in dia.edges.values() if e.kind == "ext"]
    int_edges = [e for e in dia.edges.values() if e.kind == "int"]
    int_vertices = sorted(n for n in dia.nodes if n not in dia.leg_nodes)
    dia.V_int, dia.E_int = len(int_vertices), len(int_edges)

    # self-check 3: every external exactly once; legs are 1-valent
    if len(ext_edges) != len(dia.leg_nodes):
        raise RouteError(f"diagram {idx}: {len(ext_edges)} external edges for "
                         f"{len(dia.leg_nodes)} legs")
    leg_hits = {}
    for e in ext_edges:
        for v in (e.v1, e.v2):
            if v in dia.leg_nodes:
                leg_hits[v] = leg_hits.get(v, 0) + 1
    if any(c != 1 for c in leg_hits.values()) or set(leg_hits) != set(dia.leg_nodes):
        raise RouteError(f"diagram {idx}: external legs not 1-valent: {leg_hits}")
    dia.checks.append(("externals-once", "ok"))

    # ---- spanning tree (deterministic BFS) over internal vertices
    def ekey(e):
        return (e.eid if e.eid is not None else 0, e.sym)

    adj = {v: [] for v in int_vertices}
    for e in sorted(int_edges, key=ekey):
        if e.v1 != e.v2:
            adj[e.v1].append((e.v2, e))
            adj[e.v2].append((e.v1, e))
    tree_parent = {}   # vertex -> (parent_vertex, edge)
    depth = {}
    order = []
    if int_vertices:
        root = int_vertices[0]
        depth[root] = 0
        queue = [root]
        seen = {root}
        while queue:
            v = queue.pop(0)
            order.append(v)
            for w, e in adj[v]:
                if w not in seen:
                    seen.add(w)
                    depth[w] = depth[v] + 1
                    tree_parent[w] = (v, e)
                    queue.append(w)
        if len(seen) != len(int_vertices):
            raise RouteError(f"diagram {idx}: internal graph disconnected "
                             f"({len(seen)}/{len(int_vertices)} vertices reached)")
    tree_edges = {id(e) for _, e in tree_parent.values()}
    chords = [e for e in sorted(int_edges, key=ekey) if id(e) not in tree_edges]
    dia.L = len(chords)

    # ---- loop momentum names: chord's own label by default, else rename
    if loops is not None:
        loops = list(loops)
        if len(loops) < dia.L:
            raise RouteError(f"diagram {idx}: {len(loops)} loop names for "
                             f"L={dia.L}")
        all_syms = set(dia.edges) | set(dia.externals)
        for e, new in zip(chords, loops):
            if new != e.sym and new in all_syms:
                raise RouteError(f"diagram {idx}: loop rename {new!r} "
                                 f"collides with an existing label")
        dia.loop_syms = list(loops[:dia.L])
    else:
        dia.loop_syms = [e.sym for e in chords]

    # ---- momentum assignment
    for e in ext_edges:
        sym = rename.get(e.sym, e.sym)
        e.mom = _mvec([(sym, 1)])
        e.kind = "ext"
    for e, lname in zip(chords, dia.loop_syms):
        e.mom = _mvec([(lname, 1)])
        e.kind = "chord"
    for e in int_edges:
        if e.kind == "int":
            e.kind = "tree"

    # incident edge-ends per internal vertex: (edge, sign of inflow at v)
    ends = {v: [] for v in int_vertices}
    for e in dia.edges.values():
        if e.v1 in ends:
            ends[e.v1].append((e, -1))
        if e.v2 in ends:
            ends[e.v2].append((e, +1))

    # leaf-strip: deepest vertices first; parent edge solved from conservation
    for v in sorted(order[1:], key=lambda w: -depth[w]):
        pv, pe = tree_parent[v]
        if pe.mom is not None:
            raise RouteError(f"diagram {idx}: tree edge {pe.sym} solved twice")
        acc = {}
        psign = None
        seen_parent = False
        for e, s in ends[v]:
            if e is pe and not seen_parent:
                psign, seen_parent = s, True
                continue
            if e.mom is None:
                raise RouteError(f"diagram {idx}: vertex {v}: unsolved "
                                 f"non-parent edge {e.sym} (solver order bug)")
            acc = _madd(acc, e.mom, Fraction(s))
        pe.mom = _madd({}, acc, Fraction(-psign))

    missing = [e.sym for e in int_edges if e.mom is None]
    if missing:
        raise RouteError(f"diagram {idx}: unrouted edges {missing}")

    # ---- mutation-control hook (BEFORE the self-checks; controls only)
    if mutate_edge is not None:
        if mutate_edge not in dia.edges:
            raise RouteError(f"diagram {idx}: --mutate-edge {mutate_edge!r} "
                             f"is not an edge label")
        bump = dia.externals[0] if dia.externals else dia.loop_syms[0]
        e = dia.edges[mutate_edge]
        e.mom = _madd(e.mom, _mvec([(bump, 1)]))

    # ---- self-check 1: conservation at EVERY internal vertex,
    #      modulo overall external conservation sum(externals) = 0
    def _reduce(vec):
        if not dia.externals:
            return dict(vec)
        last = dia.externals[-1]
        out = dict(vec)
        c = out.pop(last, Fraction(0))
        if c:
            for s in dia.externals[:-1]:
                out[s] = out.get(s, Fraction(0)) - c
                if not out[s]:
                    del out[s]
        return out

    for v in int_vertices:
        acc = {}
        for e, s in ends[v]:
            acc = _madd(acc, e.mom, Fraction(s))
        red = _reduce(acc)
        if red:
            raise RouteError(
                f"diagram {idx}: MOMENTUM CONSERVATION BROKEN at vertex {v}: "
                f"residual { mom_str(red, dia.sym_order()) }")
    dia.checks.append(("vertex-conservation", "ok"))

    # ---- self-check 2: Euler count
    if dia.E_int != dia.V_int + dia.L - 1:
        raise RouteError(f"diagram {idx}: Euler count E_int={dia.E_int} != "
                         f"V_int+L-1={dia.V_int + dia.L - 1}")
    dia.checks.append(("euler-count", "ok"))

    # ---- self-check 4: Betti number, independent union-find recount
    parent = {v: v for v in int_vertices}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for e in int_edges:
        ra, rb = find(e.v1), find(e.v2)
        if ra != rb:
            parent[ra] = rb
    ncomp = len({find(v) for v in int_vertices}) if int_vertices else 0
    if int_vertices and ncomp != 1:
        raise RouteError(f"diagram {idx}: {ncomp} internal components")
    betti = dia.E_int - dia.V_int + ncomp
    if betti != dia.L:
        raise RouteError(f"diagram {idx}: chord count L={dia.L} != first "
                         f"Betti number {betti}")
    dia.checks.append(("betti-count", "ok"))

    # ---- self-check 5: chords carry exactly their own loop label
    for e, lname in zip(chords, dia.loop_syms):
        if e.mom != _mvec([(lname, 1)]) and mutate_edge != e.sym:
            raise RouteError(f"diagram {idx}: chord {e.sym} momentum "
                             f"{mom_str(e.mom)} != {lname}")
    dia.checks.append(("chord-unit", "ok"))
    return dia


# ----------------------------------------------------------------- public API

def route_text(text, externals=None, loops=None, mutate_edge=None):
    """Route every diagram term in a diagrams_ output string/file content.
    Returns [RoutedDiagram]; raises RouteError on ANY parse or check failure."""
    clean = read_form_output(text)
    lines = [ln.strip() for ln in clean.splitlines() if ln.strip()]
    if not lines:
        raise RouteError("empty input")
    terms = []
    for ln in lines:
        ln = _NAME_EQ_RE.sub("", ln, count=1).rstrip(";")
        terms.extend(_split_terms(ln))
    out = []
    for i, t in enumerate(terms, 1):
        dia = _parse_term(t, i)
        out.append(_route_diagram(dia, externals=externals, loops=loops,
                                  mutate_edge=mutate_edge))
    return out


def route_file(path, **kw):
    with open(path) as f:
        return route_text(f.read(), **kw)


def routing_table(diagrams):
    rows = []
    for d in diagrams:
        order = d.sym_order()
        rows.append(f"# diagram {d.idx} (topo {d.topo}, coeff {d.coeff}): "
                    f"V_int={d.V_int} E_int={d.E_int} L={d.L} "
                    f"externals={','.join(d.externals) or '-'} "
                    f"loops={','.join(d.loop_syms) or '-'}")
        for e in d._edge_list():
            fld = f" {e.field}" if e.field else ""
            rows.append(f"#   {e.sym:>6}{fld:>5} [{e.v1}->{e.v2} {e.kind:>5}]"
                        f"  ->  {mom_str(e.mom, order)}")
    return "\n".join(rows) + "\n"


def form_block(diagrams):
    """FORM `id`-substitution block: one #procedure route<i>() per diagram,
    ready to #include; RHS contain only externals + chord labels."""
    out = ["* momentum routing generated by tools/formglue/form_route.py",
           "* call route<i>() inside the module processing diagram <i>", ""]
    for d in diagrams:
        order = d.sym_order()
        out.append(f"#procedure route{d.idx}()")
        out.append(f"* diagram {d.idx} (topo {d.topo}): loop momenta "
                   f"{','.join(d.loop_syms) or '(none: tree)'}"
                   f"; externals {','.join(d.externals) or '(vacuum)'}")
        for e in d._edge_list():
            if e.kind == "ext":
                tgt = mom_str(e.mom, order)
                if tgt != e.sym:
                    out.append(f"id {e.sym} = {tgt};")
            elif e.kind == "chord":
                tgt = mom_str(e.mom, order)
                if tgt != e.sym:
                    out.append(f"id {e.sym} = {tgt};")
            else:
                out.append(f"id {e.sym} = {mom_str(e.mom, order)};")
        out.append("#endprocedure")
        out.append("")
    return "\n".join(out)


# ------------------------------------------------------------------------ CLI

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Momentum-routing glue for FORM 5 diagrams_ output "
                    "(spanning-tree loop/tree assignment + rc-coded self-checks)")
    ap.add_argument("input", help="diagrams_ #write dump (file path)")
    ap.add_argument("--externals", default=None,
                    help="comma list: external names in leg-node order "
                         "(reorders if same set as generator's, renames if not)")
    ap.add_argument("--loops", default=None,
                    help="comma list: loop-momentum names for the chords "
                         "(default: keep each chord's own generator label)")
    ap.add_argument("--json", default=None, help="write JSON routing here")
    ap.add_argument("--form", default=None, help="write FORM id-block here")
    ap.add_argument("--table", default=None,
                    help="write routing table here ('-' = stdout)")
    ap.add_argument("--mutate-edge", default=None, metavar="LABEL",
                    help="CONTROL ONLY: corrupt LABEL's routed momentum "
                         "before the self-checks (must make the run FAIL)")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args(argv)

    ext = args.externals.split(",") if args.externals else None
    lps = args.loops.split(",") if args.loops else None
    try:
        dias = route_file(args.input, externals=ext, loops=lps,
                          mutate_edge=args.mutate_edge)
    except RouteError as e:
        print(f"FAIL form_route: {e}", file=sys.stderr)
        return 2
    if args.json:
        with open(args.json, "w") as f:
            json.dump([d.to_json() for d in dias], f, indent=1)
    if args.form:
        with open(args.form, "w") as f:
            f.write(form_block(dias))
    if args.table:
        tbl = routing_table(dias)
        if args.table == "-":
            sys.stdout.write(tbl)
        else:
            with open(args.table, "w") as f:
                f.write(tbl)
    if not args.quiet:
        ls = sorted({d.L for d in dias})
        print(f"ROUTED {len(dias)} diagrams, all self-checks PASS "
              f"(L values: {ls})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
