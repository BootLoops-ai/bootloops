#!/usr/bin/env python3
"""Render the gg -> q qbar Feynman diagrams from the generator's own data.

Input:  gg2qq_topologies.json (same directory as this script) — the
        structured edge lists parsed from FORM 5 diagrams_ output
        (node_ terms) by the momentum-routing glue: per diagram, every
        edge with its endpoints, kind (external / tree / loop chord),
        field (gluon / quark / ghost), and fermion-flow arrow.
Output: form-gg2qq-tree.png   (the 3 tree diagrams, one row)
        form-gg2qq-1loop.png  (all 30 one-loop diagrams, 6 x 5 grid)

Dependencies: matplotlib only.  Nothing here is drawn by hand: the
layout is computed from the edge lists by fixed rules (externals pinned
to the frame corners, loop vertices on a circle oriented by where their
external legs pull, hanging vertices stepped outward toward their
externals), and every line style is keyed off the field label the
generator emitted for that edge.

The JSON is produced from the generator output (the shipped
enum_1loop.frm run, plus the same enumeration at loop count zero for
the tree set) by the spanning-tree momentum router; it ships alongside
so this script runs without FORM installed.
"""
import json
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Polygon

# deterministic output: fixed fonts/params, no metadata that varies
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "mathtext.fontset": "dejavusans",
    "savefig.facecolor": "white",
    "figure.facecolor": "white",
})
PNG_META = {"Software": None}   # strip version-dependent PNG metadata

HERE = os.path.dirname(os.path.abspath(__file__))
LW = 1.5
COL = "black"
HILITE = "#b02020"

# pinned external positions: incoming gluons left, quark pair right
PIN = {"p1": (-1.55, 0.95), "p2": (-1.55, -0.95),
       "q1": (1.55, 0.95), "q2": (1.55, -0.95)}


# ---------------------------------------------------------------- geometry

def _norm(v):
    h = math.hypot(*v)
    return (v[0] / h, v[1] / h) if h > 1e-12 else (1.0, 0.0)


def _bezier(p0, p1, sag, t):
    """Quadratic Bezier from p0 to p1 bowed by sag (perpendicular), plus
    its unit tangent, at parameter t."""
    mx, my = (p0[0] + p1[0]) / 2, (p0[1] + p1[1]) / 2
    d = (p1[0] - p0[0], p1[1] - p0[1])
    n = _norm((-d[1], d[0]))
    c = (mx + sag * n[0], my + sag * n[1])
    u = 1 - t
    x = u * u * p0[0] + 2 * u * t * c[0] + t * t * p1[0]
    y = u * u * p0[1] + 2 * u * t * c[1] + t * t * p1[1]
    dx = 2 * u * (c[0] - p0[0]) + 2 * t * (p1[0] - c[0])
    dy = 2 * u * (c[1] - p0[1]) + 2 * t * (p1[1] - c[1])
    return (x, y), _norm((dx, dy))


def _pathlen(p0, p1, sag):
    L, prev = 0.0, None
    for i in range(33):
        pt, _ = _bezier(p0, p1, sag, i / 32)
        if prev:
            L += math.hypot(pt[0] - prev[0], pt[1] - prev[1])
        prev = pt
    return L


# ------------------------------------------------------------- line styles

def draw_gluon(ax, p0, p1, sag=0.0):
    """Curly line: rolling-circle curls along the (possibly bowed) path,
    with an end envelope so the curve starts and ends on the vertices."""
    L = _pathlen(p0, p1, sag)
    ncurl = max(3, int(round(L / 0.17)))
    amp = min(0.052, 0.16 * L / ncurl)
    xs, ys = [], []
    nsamp = 24 * ncurl
    for i in range(nsamp + 1):
        t = i / nsamp
        (bx, by), (tx, ty) = _bezier(p0, p1, sag, t)
        nx, ny = -ty, tx
        th = 2 * math.pi * ncurl * t
        env = min(1.0, t / 0.06, (1 - t) / 0.06)
        ax_off = amp * env * (math.sin(th) * nx + math.cos(th) * tx)
        ay_off = amp * env * (math.sin(th) * ny + math.cos(th) * ty)
        xs.append(bx + ax_off)
        ys.append(by + ay_off)
    ax.plot(xs, ys, color=COL, lw=LW * 0.85, solid_capstyle="round")


def draw_line(ax, p0, p1, sag=0.0, dashed=False):
    xs, ys = [], []
    for i in range(41):
        (x, y), _ = _bezier(p0, p1, sag, i / 40)
        xs.append(x)
        ys.append(y)
    ax.plot(xs, ys, color=COL, lw=LW,
            ls=(0, (4, 2.6)) if dashed else "-", solid_capstyle="round")


def draw_arrowhead(ax, p0, p1, sag, flip):
    """Filled triangle at the path midpoint, pointing along the fermion
    flow (flip=True when flow runs p1 -> p0)."""
    (mx, my), (tx, ty) = _bezier(p0, p1, sag, 0.5)
    if flip:
        tx, ty = -tx, -ty
    nx, ny = -ty, tx
    s = 0.085
    tip = (mx + s * 0.75 * tx, my + s * 0.75 * ty)
    b1 = (mx - s * 0.55 * tx + s * 0.5 * nx, my - s * 0.55 * ty + s * 0.5 * ny)
    b2 = (mx - s * 0.55 * tx - s * 0.5 * nx, my - s * 0.55 * ty - s * 0.5 * ny)
    ax.add_patch(Polygon([tip, b1, b2], closed=True, fc=COL, ec="none",
                         zorder=5))


def draw_edge(ax, pos, e, sag):
    p0, p1 = pos[e["v1"]], pos[e["v2"]]
    if e["field"] == "gluon":
        draw_gluon(ax, p0, p1, sag)
    else:
        draw_line(ax, p0, p1, sag, dashed=(e["field"] == "ghost"))
        if e.get("arrow"):
            draw_arrowhead(ax, p0, p1, sag, flip=(e["arrow"][0] == e["v2"]))


# ----------------------------------------------------------------- layout

def internal_edges(d):
    return [e for e in d["edges"] if e["kind"] != "ext"]


def tree_layout(d):
    """4-point tree, two internal vertices: canonical channel layout
    read off the external partition."""
    ints = sorted({v for e in internal_edges(d)
                   for v in (e["v1"], e["v2"])})
    att = {}   # internal vertex -> set of external labels
    for e in d["edges"]:
        if e["kind"] == "ext":
            v = e["v2"] if e["v2"] in ints else e["v1"]
            att.setdefault(v, set()).add(e["label"])
    a = next(v for v in ints if "p1" in att[v])
    b = next(v for v in ints if v != a)
    if "p2" in att[a]:            # p1 with p2: s-channel, horizontal
        pos = {a: (-0.55, 0.0), b: (0.55, 0.0)}
    else:                         # p1 with q1 (t) or q2 (u): vertical
        pos = {a: (0.0, 0.6), b: (0.0, -0.6)}
    return pos


def channel_name(d):
    att = {}
    for e in d["edges"]:
        if e["kind"] == "ext":
            att.setdefault(max(e["v1"], e["v2"]), set()).add(e["label"])
    for legs in att.values():
        if "p1" in legs:
            if "p2" in legs:
                return "s-channel"
            return "t-channel" if "q1" in legs else "u-channel"
    raise ValueError("p1 attachment not found")


def _seg_cross(a, b, c, d):
    """Crossing point of the open interiors of ab and cd, else None."""
    eps = 1e-9
    for p, q in ((a, c), (a, d), (b, c), (b, d)):
        if abs(p[0] - q[0]) < eps and abs(p[1] - q[1]) < eps:
            return None         # shared endpoint: not a crossing
    d1 = (b[0] - a[0], b[1] - a[1])
    d2 = (d[0] - c[0], d[1] - c[1])
    den = d1[0] * d2[1] - d1[1] * d2[0]
    if abs(den) < eps:
        return None
    t = ((c[0] - a[0]) * d2[1] - (c[1] - a[1]) * d2[0]) / den
    u = ((c[0] - a[0]) * d1[1] - (c[1] - a[1]) * d1[0]) / den
    if eps < t < 1 - eps and eps < u < 1 - eps:
        return (a[0] + t * d1[0], a[1] + t * d1[1])
    return None


def _pt_seg_dist(p, a, b):
    dx, dy = b[0] - a[0], b[1] - a[1]
    L2 = dx * dx + dy * dy
    if L2 < 1e-18:
        return math.hypot(p[0] - a[0], p[1] - a[1])
    t = max(0.0, min(1.0, ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / L2))
    return math.hypot(p[0] - a[0] - t * dx, p[1] - a[1] - t * dy)


def _polyline(p0, p1, sag, nseg=10):
    if sag == 0.0:
        return [p0, p1]
    return [_bezier(p0, p1, sag, i / nseg)[0] for i in range(nseg + 1)]


def _score(d, pos):
    """(crossings, clearance violations, total leg length): lower is
    better.  Elements are the internal edges (with their bows) and the
    straight external legs, all as polylines.  Clearance violations
    count: a vertex within 0.19 of a non-incident element, two vertices
    within 0.28, a crossing within 0.32 of a vertex, two non-adjacent
    elements within 0.15 anywhere, and two lines meeting at a vertex at
    a shallow (<12 degree) angle."""
    sag = sagittas(d)
    elems = []                            # (set of incident vertices, pts)
    for e in internal_edges(d):
        elems.append(({e["v1"], e["v2"]},
                      _polyline(pos[e["v1"]], pos[e["v2"]], sag[id(e)])))
    leglen = 0.0
    for e in d["edges"]:
        if e["kind"] != "ext":
            continue
        v = e["v2"] if e["v2"] in pos else e["v1"]
        elems.append(({v}, [PIN[e["label"]], pos[v]]))
        leglen += math.dist(PIN[e["label"]], pos[v])
    ncross = 0
    nclear = 0
    for i in range(len(elems)):
        for j in range(i + 1, len(elems)):
            pi, pj = elems[i][1], elems[j][1]
            adjacent = bool(elems[i][0] & elems[j][0])
            mind = float("inf")
            for s in range(len(pi) - 1):
                for t in range(len(pj) - 1):
                    x = _seg_cross(pi[s], pi[s + 1], pj[t], pj[t + 1])
                    if x is not None:
                        ncross += 1
                        if min(math.dist(x, p) for p in pos.values()) < 0.32:
                            nclear += 1
                    if not adjacent:
                        mind = min(mind,
                                   _pt_seg_dist(pj[t], pi[s], pi[s + 1]),
                                   _pt_seg_dist(pj[t + 1], pi[s], pi[s + 1]),
                                   _pt_seg_dist(pi[s], pj[t], pj[t + 1]),
                                   _pt_seg_dist(pi[s + 1], pj[t], pj[t + 1]))
            if not adjacent and mind < 0.15:
                nclear += 1
    verts = list(pos.items())
    for v, p in verts:
        for inc, pts in elems:
            if v in inc:
                continue
            if min(_pt_seg_dist(p, pts[s], pts[s + 1])
                   for s in range(len(pts) - 1)) < 0.19:
                nclear += 1
    for i in range(len(verts)):
        for j in range(i + 1, len(verts)):
            if math.dist(verts[i][1], verts[j][1]) < 0.28:
                nclear += 1
    # shallow angles between the lines meeting at each vertex
    for v, p in verts:
        dirs = []
        for inc, pts in elems:
            if v not in inc:
                continue
            ends = [pts[0], pts[-1]]
            other = max(ends, key=lambda q: math.dist(p, q))
            near = pts[1] if math.dist(pts[0], p) < math.dist(pts[-1], p) \
                else pts[-2]
            q = near if len(pts) > 2 else other
            dirs.append(math.atan2(q[1] - p[1], q[0] - p[0]))
        for i in range(len(dirs)):
            for j in range(i + 1, len(dirs)):
                dphi = abs((dirs[i] - dirs[j] + math.pi) % (2 * math.pi)
                           - math.pi)
                if dphi < math.radians(12):
                    nclear += 1
    return (ncross, nclear, round(leglen, 6))


def loop_layout(d):
    """One-loop layout as a deterministic search over canonical
    configurations: cycle vertices on a circle (24 offsets x 2
    traversal directions), each hanging vertex stepped 0.5 from its
    cycle attachment along one of a fixed candidate list of directions
    (toward each of its externals' pins, toward their mean, or radially
    outward).  The configuration minimizing (line crossings, clearance
    violations, total external-leg length) wins; ties break on
    enumeration order, so the result is reproducible."""
    edges = internal_edges(d)
    ints = sorted({v for e in edges for v in (e["v1"], e["v2"])})
    chord = next(e for e in edges if e["kind"] == "chord")
    adj = {v: [] for v in ints}
    for e in edges:
        if e["kind"] == "tree":
            adj[e["v1"]].append(e["v2"])
            adj[e["v2"]].append(e["v1"])
    # cycle = chord + unique tree path between its endpoints
    a, b = chord["v1"], chord["v2"]
    parent = {a: None}
    stack = [a]
    while stack:
        v = stack.pop()
        for w in adj[v]:
            if w not in parent:
                parent[w] = v
                stack.append(w)
    path = [b]
    while path[-1] != a:
        path.append(parent[path[-1]])
    cycle = list(reversed(path))

    # external labels carried by each hanging vertex; hang parents
    att_v = {}
    for e in d["edges"]:
        if e["kind"] == "ext":
            att_v[e["label"]] = e["v2"] if e["v2"] in ints else e["v1"]
    cyc = set(cycle)
    hang_ext = {}                         # hanging vertex -> [labels]
    for lab, v in att_v.items():
        if v not in cyc:
            hang_ext.setdefault(v, []).append(lab)
    hangs = []                            # (vertex, cycle parent)
    for v in sorted(set(ints) - cyc):
        par = [w for w in adj[v] if w in cyc]
        if len(par) != 1 or v not in hang_ext:
            raise ValueError(f"diagram {d['diagram']}: unsupported "
                             f"hanging vertex {v}")
        hangs.append((v, par[0]))

    n = len(cycle)
    R = {2: 0.40, 3: 0.48, 4: 0.55}.get(n, 0.55)
    best = None
    for rev in (0, 1):
        seq = list(reversed(cycle)) if rev else cycle
        for k in range(24):
            off = k * math.pi / 12
            base = {v: (R * math.cos(off + 2 * math.pi * j / n),
                        R * math.sin(off + 2 * math.pi * j / n))
                    for j, v in enumerate(seq)}
            # hang placement candidates for this slot assignment:
            # direction (toward each pin, their mean, radially outward,
            # radial rotated +-45 degrees) x step length
            cands = []
            for v, par in hangs:
                pp = base[par]
                dirs = [_norm((PIN[l][0] - pp[0], PIN[l][1] - pp[1]))
                        for l in sorted(hang_ext[v])]
                mx = sum(PIN[l][0] for l in hang_ext[v]) / len(hang_ext[v])
                my = sum(PIN[l][1] for l in hang_ext[v]) / len(hang_ext[v])
                if math.hypot(mx - pp[0], my - pp[1]) > 0.15:
                    dirs.append(_norm((mx - pp[0], my - pp[1])))
                rx, ry = _norm(pp)        # radially outward + rotations
                c45 = math.cos(math.pi / 4)
                dirs += [(rx, ry),
                         (c45 * (rx - ry), c45 * (rx + ry)),
                         (c45 * (rx + ry), c45 * (ry - rx))]
                cands.append([(u, step) for u in dirs
                              for step in (0.5, 0.75)])
            idx = [0] * len(hangs)
            while True:
                pos = dict(base)
                for (v, par), opts, i in zip(hangs, cands, idx):
                    u, step = opts[i]
                    pos[v] = (base[par][0] + step * u[0],
                              base[par][1] + step * u[1])
                key = _score(d, pos) + (rev, k, tuple(idx))
                if best is None or key < best[0]:
                    best = (key, pos)
                # advance the mixed-radix counter
                for slot in range(len(idx) - 1, -1, -1):
                    idx[slot] += 1
                    if idx[slot] < len(cands[slot]):
                        break
                    idx[slot] = 0
                else:
                    break
                if not hangs:
                    break
    return best[1]


def sagittas(d):
    """Bow parallel internal edges apart; everything else stays straight."""
    groups = {}
    for e in internal_edges(d):
        groups.setdefault(frozenset((e["v1"], e["v2"])), []).append(e)
    sag = {}
    for g in groups.values():
        if len(g) == 1:
            sag[id(g[0])] = 0.0
        else:
            for e, s in zip(sorted(g, key=lambda e: e["label"]),
                            (0.24, -0.24)):
                sag[id(e)] = s
    return sag


# ----------------------------------------------------------------- tiles

def draw_diagram(ax, d, pos):
    sag = sagittas(d)
    for e in d["edges"]:
        if e["kind"] == "ext":
            p_ext = PIN[e["label"]]
            v_int = e["v2"] if e["v2"] in pos else e["v1"]
            epos = dict(pos)
            epos[e["v1"] if e["v1"] not in pos else e["v2"]] = p_ext
            # externals drawn pinned-point -> vertex, straight
            draw_edge(ax, {e["v1"]: epos[e["v1"]], e["v2"]: epos[e["v2"]]},
                      e, 0.0)
        else:
            draw_edge(ax, pos, e, sag[id(e)])
    for v, (x, y) in pos.items():
        ax.plot([x], [y], marker="o", ms=3.4, mfc=COL, mec=COL, zorder=6)
    ax.set_xlim(-1.75, 1.75)
    ax.set_ylim(-1.35, 1.35)
    ax.set_aspect("equal")
    ax.axis("off")


def render_tree(doc, out):
    fig, axes = plt.subplots(1, 3, figsize=(9.2, 2.9))
    order = {"t-channel": 0, "u-channel": 1, "s-channel": 2}
    tiles = sorted(doc["tree"], key=lambda d: order[channel_name(d)])
    for ax, d in zip(axes, tiles):
        draw_diagram(ax, d, tree_layout(d))
        ax.set_title(channel_name(d), fontsize=11, pad=4)
        for lab, (x, y) in PIN.items():
            ax.text(x * 1.07, y * 1.12, f"${lab[0]}_{lab[1]}$",
                    ha="center", va="center", fontsize=9)
        ax.set_xlim(-1.95, 1.95)
        ax.set_ylim(-1.45, 1.45)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.86, bottom=0.03,
                        wspace=0.05)
    fig.savefig(out, dpi=170, metadata=PNG_META)
    plt.close(fig)


def render_oneloop(doc, out, highlight=21):
    rows, cols = 5, 6
    fig, axes = plt.subplots(rows, cols, figsize=(13.2, 10.6))
    dias = sorted(doc["oneloop"], key=lambda d: d["diagram"])
    assert len(dias) == rows * cols
    for ax, d in zip(axes.flat, dias):
        draw_diagram(ax, d, loop_layout(d))
        ax.text(-1.68, 1.22, str(d["diagram"]), fontsize=10,
                ha="left", va="top", color="#555555")
        if d["diagram"] == highlight:
            ax.add_patch(FancyBboxPatch(
                (-1.7, -1.3), 3.4, 2.6,
                boxstyle="round,pad=0.02,rounding_size=0.08",
                fill=False, ec=HILITE, lw=1.6, zorder=7))
    fig.subplots_adjust(left=0.005, right=0.995, top=0.995, bottom=0.005,
                        wspace=0.04, hspace=0.04)
    fig.savefig(out, dpi=150, metadata=PNG_META)
    plt.close(fig)


def main():
    with open(os.path.join(HERE, "gg2qq_topologies.json")) as f:
        doc = json.load(f)
    assert len(doc["tree"]) == 3 and len(doc["oneloop"]) == 30
    render_tree(doc, "form-gg2qq-tree.png")
    render_oneloop(doc, "form-gg2qq-1loop.png")
    print("wrote form-gg2qq-tree.png (3 tiles), "
          "form-gg2qq-1loop.png (30 tiles)")


if __name__ == "__main__":
    main()
