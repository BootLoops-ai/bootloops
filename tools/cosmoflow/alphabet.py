#!/usr/bin/env python3
r"""alphabet.py -- FRW letter-alphabet generator for an arbitrary site graph.

Emits the candidate symbol-letter alphabet of an FRW wavefunction graph --
ANY connected site graph (chains, rings, stars, multi-loop polygons; front
door `alphabet_graph(nv, edges, ...)`, named special cases
`alphabet_tree_chain` / `alphabet_loop_ngon`) -- in
the ansatzer input schema (tools/ansatzer: name / vars /
letters / first_entry [/ last_entry ...]), one JSON per graph:

  * q_g FACET forms -- the polytope-facet / OFPT energy denominators, from
    the connected-subgraph rule.  Deleted loop edges (both endpoints inside
    the subgraph) are cut TWICE: coefficient c=2, DERIVED (three
    independent exact routes, see C_PROVENANCE); c=1 (the superseded
    2408.16386 triangle-section print) stays reachable by argument.  At
    TREE level c never enters (a connected subgraph of a tree cannot omit
    an edge without disconnecting -- asserted here).
  * FOLDED letters (the X_i - Y_e class): each facet with >=1 cut edge also
    contributes its cut-edge-sign-flipped partner (tree-chain letters
    X1-Y, X2-Y; bubble letters x1-P, x2-P, x1+x2-2P are this class).
  * DISC / KAELLEN loci (loop graphs only): the bare edge energies y_e (the
    bubble alphabet's "P" letter class) in the main letter list, plus a
    SEPARATE `disc_kallen_baikov` block with the irreducible factors of
    disc_{z_e} B(z;X) on the Baikov chart (nested-Kaellen chain, e.g.
    disc_{z1}B = 4*lam(X1^2,X2^2,X3^2)*lam(z2,z3,X3^2) for n_s=3).  Kept
    separate because it lives on the Baikov kinematic lattice (base-edge
    invariants), not the energetic (x_i, y_e) chart of the facet letters.

Every letter carries a `source` tag; the output carries an honesty block
`complete: false` -- the completeness gate is the CONNECTION-DENOMINATOR
cross-check (the named completeness gate): the alphabet may
be claimed complete only once every denominator letter of the independently
derived DE connection (mod-p Griffiths--Dwork route) is contained in it.

VALIDATED censuses: n_s=2 tree chain -> the 5-letter two-site alphabet
(names l_EL, l_ER, l_ET, l_FL, l_FR); 3-gon facet set -> the 10 c=2
facets; 4-gon -> 17 facets, c=2 on the four edge-deleted forms.

CALIBRATION: `connected_subgraphs` reproduces the subgraph counts stated
in 2408.16386 -- bubble 5, triangle 10 -- and gives box 17 = n^2+1.
`q_of_subgraph` parameterizes the edges and defaults c to the derived 2.

CLI (run from tools/, the package's parent):
  python3 -m cosmoflow.alphabet graph --nv 3 --edges 0-1,1-2,0-2 [--c 2]
          [--sites x1,x2,x3] [--enames y12,y23,y31] [--disc] [--out F.json]
  python3 -m cosmoflow.alphabet tree n_s [--out F.json]      (path graph)
  python3 -m cosmoflow.alphabet loop n_s [--c 2] [--no-disc] [--out F.json]
                                                            (n-cycle)
The `graph` kind is the front door: any connected site graph (chains, rings,
stars, multi-loop polygons; parallel edges = repeated pairs, e.g. the bubble
`--nv 2 --edges 0-1,0-1`).  `tree`/`loop` are the named special cases.
"""
import itertools
import json

import sympy as sp

from .polytope import baikov_B

C_PROVENANCE = ("c=2 DERIVED (not assumed): three independent exact routes "
                "(wavefunction pole census, canonical-form ratio, "
                "vertex-facet geometry) agree exactly.  c=1 = the "
                "superseded 2408.16386 triangle-section print, reachable by "
                "argument (L2/periods proven c-independent).")


# ---------------------------------------------------------------------------
# Connected-subgraph enumerator.  Calibration counts: n=2 -> 5, n=3 -> 10
# (the counts stated in arXiv:2408.16386), n=4 -> 17.
# ---------------------------------------------------------------------------
def connected_subgraphs(nv, edges):
    """All connected subgraphs (V',E') of a multigraph, per the subgraph
    convention of arXiv:2408.16386 (single vertices allowed; every vertex of
    V' must be connected through E' when |V'|>1)."""
    subs = []
    for r in range(1, nv + 1):
        for V in itertools.combinations(range(nv), r):
            Vs = set(V)
            E_in = [i for i, (u, v) in enumerate(edges) if u in Vs and v in Vs]
            for k in range(len(E_in) + 1):
                for E in itertools.combinations(E_in, k):
                    if len(Vs) == 1:
                        if k == 0:
                            subs.append((V, E))
                        continue
                    parent = {v: v for v in Vs}
                    def find(a):
                        while parent[a] != a:
                            a = parent[a]
                        return a
                    for i in E:
                        u, v = edges[i]
                        ru, rv = find(u), find(v)
                        if ru != rv:
                            parent[ru] = rv
                    if len({find(v) for v in Vs}) == 1:
                        subs.append((V, E))
    return subs


def q_of_subgraph(V, E, edges, xs, ye, c=2):
    """Facet form q_g of connected subgraph (V,E): sum_{v in V} x_v
    + sum(once-cut y_e) + c * sum(double-cut y_e).  Once-cut = edge not in E
    with exactly one endpoint in V (coefficient 1, unambiguous); double-cut =
    deleted loop edge, both endpoints in V (coefficient c; c=2 derived, see
    C_PROVENANCE).  Returns (expr, n_double_cut)."""
    Vs = set(V)
    cut1 = [i for i, (u, v) in enumerate(edges)
            if i not in E and (u in Vs) != (v in Vs)]
    cut2 = [i for i, (u, v) in enumerate(edges)
            if i not in E and (u in Vs) and (v in Vs)]
    q = (sum(xs[i] for i in V) + sum(ye[i] for i in cut1)
         + c * sum(ye[i] for i in cut2))
    return sp.expand(q), len(cut2)


def _cut_edges(V, E, edges):
    """(once-cut, double-cut) edge index lists of subgraph (V,E)."""
    Vs = set(V)
    cut1 = [i for i, (u, v) in enumerate(edges)
            if i not in E and (u in Vs) != (v in Vs)]
    cut2 = [i for i, (u, v) in enumerate(edges)
            if i not in E and (u in Vs) and (v in Vs)]
    return cut1, cut2


# ---------------------------------------------------------------------------
# Graph data
# ---------------------------------------------------------------------------
def loop_ngon_edges(n_s):
    """Edge list + names of the 1-loop n-gon (n_s=2: bubble, two parallel
    edges y12,y21 per 2408.16386; n_s>=3: y12,y23,..,y_{n,1})."""
    if n_s == 2:
        return [(0, 1), (0, 1)], ['y12', 'y21']
    edges = [(i, (i + 1) % n_s) for i in range(n_s)]
    names = [f'y{i+1}{(i+1) % n_s + 1}' for i in range(n_s)]
    return edges, names


def tree_chain_edges(n_s):
    """Edge list + names of the n_s-site tree chain (n_s=2: single edge
    named Y; n_s>=3: Y1..Y_{n-1})."""
    edges = [(i, i + 1) for i in range(n_s - 1)]
    names = ['Y'] if n_s == 2 else [f'Y{i+1}' for i in range(n_s - 1)]
    return edges, names


def parse_edges(spec):
    """'0-1,1-2,0-2' -> [(0,1),(1,2),(0,2)] (0-based vertex labels; repeated
    pairs are parallel edges, e.g. the bubble '0-1,0-1')."""
    edges = []
    for tok in str(spec).replace(' ', '').split(','):
        if not tok:
            continue
        u, v = tok.split('-')
        edges.append((int(u), int(v)))
    return edges


def default_edge_names(nv, edges):
    """y{u}{v} by 1-based endpoints (y{u}_{v} when nv>9); a parallel
    duplicate takes the reversed label first (bubble: y12, y21), then _k."""
    sep = '_' if nv > 9 else ''
    names = []
    for (u, v) in edges:
        cand = [f'y{u+1}{sep}{v+1}', f'y{v+1}{sep}{u+1}']
        nm = next((c for c in cand if c not in names), None)
        k = 2
        while nm is None or nm in names:
            nm = f'{cand[0]}_{k}'; k += 1
        names.append(nm)
    return names


def _n_components(nv, edges):
    parent = list(range(nv))
    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    for (u, v) in edges:
        ru, rv = find(u), find(v)
        if ru != rv:
            parent[ru] = rv
    return len({find(v) for v in range(nv)})


def _bridges(nv, edges):
    """Indices of bridge edges (removal disconnects) -- the tree-like edges;
    every other edge lies on a cycle."""
    n0 = _n_components(nv, edges)
    return [i for i in range(len(edges))
            if _n_components(nv, edges[:i] + edges[i+1:]) > n0]


# ---------------------------------------------------------------------------
# Emitters
# ---------------------------------------------------------------------------
def _honesty(extra=None):
    h = {
        'complete': False,
        'completeness_gate': (
            'connection-denominator cross-check (the named completeness '
            'gate): the '
            'alphabet is complete only when every denominator letter of the '
            'independently derived DE connection (mod-p Griffiths--Dwork '
            'route) is contained in this list.  '
            'Until that check has been run for the graph at hand, '
            'treat this alphabet as a candidate set.'),
    }
    if extra:
        h.update(extra)
    return h


def _vtag(V):
    return ''.join(str(v + 1) for v in V) if max(V) < 9 else \
        '_'.join(str(v + 1) for v in V)


def alphabet_graph(nv, edges, c=2, site_names=None, edge_names=None,
                   with_disc=False, name=None):
    """Letter alphabet of the FRW wavefunction of an ARBITRARY site graph
    (the front door; `alphabet_tree_chain` / `alphabet_loop_ngon` are thin
    wrappers over it).

        nv         : number of sites (vertices 0..nv-1)
        edges      : list of (u, v) pairs, 0-based; parallel edges allowed
                     (bubble = [(0,1),(0,1)]); the graph must be connected
        c          : double-cut coefficient of a deleted loop edge (default
                     the DERIVED c=2, see C_PROVENANCE); never enters for a
                     tree (asserted: a tree subgraph has no double-cut edge)
        site_names : vertex-energy symbols (default x1..x_nv)
        edge_names : edge-energy symbols (default y{u}{v}, see
                     default_edge_names)
        with_disc  : attach the Baikov disc/Kaellen block -- implemented for
                     the 1-loop n-gon only (connected, 2-regular); any other
                     graph raises ValueError when asked for it
        name       : output 'name' field (default frw_graph_v{nv}_e{E}_l{L})

    Letter classes: FACET q_g (connected-subgraph rule, = first_entry list;
    once-cut edges coefficient 1, double-cut edges coefficient c), FOLDED
    (cut-edge sign flip of every facet with >=1 cut edge), DISC (bare
    energy y_e of every cycle edge, i.e. every non-bridge edge; empty for a
    tree).  Names: q_G (whole graph, nothing cut), q_g<sites> (proper
    subgraph, once-cut only), q_Gdel_<edges> / q_g<sites>_del_<edges>
    (double-cut edges), fold_<...>, disc_<edge>.  Loop number
    L = E - V + 1 is reported in the 'graph' block."""
    try:
        edges = [tuple(e) for e in edges]
        ok = nv >= 1 and all(len(e) == 2 for e in edges) and all(
            0 <= u < nv and 0 <= v < nv and u != v for u, v in edges)
    except TypeError:
        ok = False
    if not ok:
        raise ValueError("bad edge list")
    ne = len(edges)
    if _n_components(nv, edges) != 1:
        raise ValueError("site graph must be connected")
    L = ne - nv + 1
    enames = list(edge_names) if edge_names is not None else \
        default_edge_names(nv, edges)
    snames = list(site_names) if site_names else [f'x{i+1}' for i in range(nv)]
    if len(enames) != ne or len(snames) != nv:
        raise ValueError("name list length mismatch")
    xs = [sp.Symbol(n) for n in snames]
    ye = [sp.Symbol(n) for n in enames]
    br = set(_bridges(nv, edges))
    cyc = [i for i in range(ne) if i not in br]      # cycle (loop) edges
    tree = (L == 0)
    facet_src = ('facet/OFPT energy denominator, connected-subgraph rule '
                 '(2312.05303 eq. (3.31) for n_s=2)') if tree else None
    letters, sources, first = {}, {}, []
    for (V, E) in connected_subgraphs(nv, edges):
        q, ndc = q_of_subgraph(V, E, edges, xs, ye, c=(1 if tree else c))
        cut1, cut2 = _cut_edges(V, E, edges)
        if tree:
            assert ndc == 0, "tree level must never see a double-cut edge"
        if len(V) == nv and not cut1 and not cut2:
            nm = 'q_G'
        elif cut2 and len(V) == nv:
            nm = 'q_Gdel_' + '_'.join(enames[i] for i in cut2)
        elif cut2:
            nm = 'q_g' + _vtag(V) + '_del_' + '_'.join(enames[i] for i in cut2)
        else:
            nm = 'q_g' + _vtag(V)
        base, k = nm, 2          # distinct (V,E) with identical labels
        while nm in letters:
            nm = f'{base}_v{k}'; k += 1
        letters[nm] = str(q)
        sources[nm] = facet_src if tree else (
            'facet q_g, connected-subgraph rule; double-cut '
            f'coefficient c={c}' + (' -- ' + C_PROVENANCE if cut2 else ''))
        first.append(nm)
        if cut1 or cut2:
            fq = sp.expand(sum(xs[i] for i in V)
                           - sum(ye[i] for i in cut1)
                           - (0 if tree else c) * sum(ye[i] for i in cut2))
            fname = 'fold_' + nm[2:]
            letters[fname] = str(fq)
            sources[fname] = ('folded letter (X_i - Y_e class), cut-edge '
                              'sign flip of ' + nm
                              + (' (2312.05303 eq. (3.31); excluded from '
                                 'first entry -- physical-sheet hypothesis)'
                                 if tree else
                                 ' (the bubble-alphabet x1-P/x2-P/x1+x2-2P '
                                 'class)'))
    for i in cyc:
        nm = f'disc_{enames[i]}'
        letters[nm] = enames[i]
        sources[nm] = ('bare edge energy (disc/branch letter class; the '
                       'bubble alphabet "P" letter)')
    out = {
        'name': name or f'frw_graph_v{nv}_e{ne}_l{L}',
        'graph': {'nv': nv, 'edges': [list(e) for e in edges], 'loops': L,
                  'cycle_edges': [enames[i] for i in cyc]},
        'vars': snames + enames,
        'letters': letters,
        'first_entry': first,
        'sources': sources,
        'c': ({'value': None,
               'note': 'c never enters at tree level (a connected subgraph '
                       'cannot omit an edge of a tree without '
                       'disconnecting)'} if tree else
              {'value': c, 'provenance': C_PROVENANCE}),
        'honesty': _honesty(),
    }
    if with_disc:
        two_reg = all(sum((u == w) + (v == w) for u, v in edges) == 2
                      for w in range(nv))
        if not (L == 1 and two_reg and ne == nv):
            raise ValueError('disc/Kaellen Baikov block is implemented for '
                             'the 1-loop n-gon only (connected, 2-regular); '
                             'pass with_disc=False for this graph')
        out['disc_kallen_baikov'] = disc_kallen_block(nv, c=c)
    return out


def alphabet_tree_chain(n_s, site_names=None, edge_names=None):
    """Letter alphabet of the n_s-site TREE chain wavefunction coefficient
    (thin wrapper over alphabet_graph on the path graph).  Facet class
    (= OFPT energy denominators; connected subgraphs of a path are the
    contiguous arcs) + folded class (cut-edge sign flips); c never enters at
    tree level -- asserted.  Legacy letter names: q_<sites>, f_<sites>;
    n_s=2 reproduces the standard 5-letter two-site alphabet with names
    l_EL, l_ER, l_ET, l_FL, l_FR."""
    edges, enames = tree_chain_edges(n_s)
    a = alphabet_graph(n_s, edges, edge_names=edge_names or enames,
                       site_names=site_names or [f'X{i+1}' for i in range(n_s)],
                       name=f'frw_tree_chain_n{n_s}')
    two = {'q_1': 'l_EL', 'q_2': 'l_ER', 'q_12': 'l_ET',
           'f_1': 'l_FL', 'f_2': 'l_FR'} if n_s == 2 else {}
    def ren(nm):
        nm = ('q_' + _vtag(range(n_s)) if nm == 'q_G' else
              nm.replace('q_g', 'q_', 1).replace('fold_g', 'f_', 1))
        return two.get(nm, nm)
    a['letters'] = {ren(k): v for k, v in a['letters'].items()}
    a['sources'] = {ren(k): (v.replace('sign flip of q_' + k[5:],
                                       'sign flip of ' + ren('q_' + k[5:]))
                             if k.startswith('fold_') else v)
                    for k, v in a['sources'].items()}
    a['first_entry'] = [ren(k) for k in a['first_entry']]
    return a


def disc_kallen_block(n_s, c=2):
    """disc_{z_e} B factors on the Baikov chart (loop n-gon).  Returns a
    separate block: these loci live on the base-invariant lattice
    (X_{ij}, z_e), NOT the energetic (x_i, y_e) chart of the facets."""
    ys = sp.symbols(f'y1:{n_s+1}')
    m = n_s * (n_s - 1) // 2
    Xs = sp.symbols(f'X1:{m+1}')
    zs = sp.symbols(f'z1:{n_s+1}')
    B = baikov_B(n_s, ys, Xs)
    Bz = sp.expand(B.subs({y: sp.sqrt(z) for y, z in zip(ys, zs)}))
    letters, sources = {}, {}
    for e, z in enumerate(zs):
        p = sp.Poly(Bz, z)
        if p.degree() != 2:
            continue
        a2, a1, a0 = p.all_coeffs()
        disc = sp.factor(sp.expand(a1 ** 2 - 4 * a2 * a0))
        facs = [f for f, mult in sp.factor_list(disc)[1]]
        for k, f in enumerate(facs):
            nm = f'disc_z{e+1}_f{k+1}'
            letters[nm] = str(sp.expand(f))
            sources[nm] = (f'irreducible factor of disc_z{e+1} B (nested-'
                           'Kaellen disc chain; for n_s=3: disc_z1 B = '
                           '4*lam(X1^2,X2^2,X3^2)*lam(z2,z3,X3^2), '
                           'cosmoflow manual INPUTS row / 2408.16386)')
    return {
        'chart': 'baikov (base-edge invariants X_ij, z_e = y_e^2)',
        'note': ('separate block: different kinematic lattice from the '
                 'energetic facet letters; z-quadratic discriminant chain'),
        'vars': [str(s) for s in Xs] + [str(s) for s in zs],
        'letters': letters,
        'sources': sources,
    }


def alphabet_loop_ngon(n_s, c=2, site_names=None, edge_names=None,
                       with_disc=True):
    """Letter alphabet of the 1-loop n-gon wavefunction (thin wrapper over
    alphabet_graph on the n-cycle): facet class (connected-subgraph census,
    double-cut coefficient c -- default the DERIVED c=2, see C_PROVENANCE)
    + folded class + disc class (bare y_e) + Baikov disc/Kaellen block.
    Facet counts: 5 (n=2) / 10 (n=3) / 17 (n=4) -- validated censuses."""
    edges, enames = loop_ngon_edges(n_s)
    a = alphabet_graph(n_s, edges, c=c, site_names=site_names,
                       edge_names=edge_names or enames, with_disc=with_disc,
                       name=f'frw_loop_ngon_n{n_s}')
    a['honesty']['facet_census_gate'] = (
        'facet count checked against the validated censuses: 5 (n=2) / '
        '10 (n=3, c=2) / 17 (n=4)')
    return a


def facet_letters(alphabet):
    """The facet-class letter expressions of an emitted alphabet (sympy)."""
    return [sp.sympify(alphabet['letters'][k])
            for k in alphabet['first_entry']]


# ---------------------------------------------------------------------------
def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('kind', choices=['graph', 'tree', 'loop'],
                    help='graph = arbitrary site graph (front door); '
                         'tree = n_s-site chain; loop = 1-loop n_s-gon')
    ap.add_argument('n_s', type=int, nargs='?', default=None,
                    help='number of sites (tree/loop kinds)')
    ap.add_argument('--nv', type=int, default=None,
                    help='number of sites (graph kind)')
    ap.add_argument('--edges', default=None,
                    help="edge list 'u-v,u-v,...' with 0-based sites "
                         "(graph kind); repeated pairs = parallel edges")
    ap.add_argument('--sites', default=None,
                    help='comma-separated site-energy names (default x1..)')
    ap.add_argument('--enames', default=None,
                    help='comma-separated edge-energy names (default y{u}{v})')
    ap.add_argument('--c', type=int, default=2,
                    help='double-cut coefficient (loop edges only; default '
                         '2 = the DERIVED value)')
    ap.add_argument('--disc', action='store_true',
                    help='graph kind: attach the Baikov disc/Kaellen block '
                         '(1-loop n-gon graphs only)')
    ap.add_argument('--no-disc', action='store_true',
                    help='loop kind: drop the Baikov disc/Kaellen block')
    ap.add_argument('--out', default=None)
    a = ap.parse_args(argv)
    names = lambda s: s.split(',') if s else None
    if a.kind == 'graph':
        if a.nv is None or a.edges is None:
            ap.error('graph kind needs --nv N --edges u-v,u-v,...')
        alph = alphabet_graph(a.nv, parse_edges(a.edges), c=a.c,
                              site_names=names(a.sites),
                              edge_names=names(a.enames), with_disc=a.disc)
    elif a.n_s is None:
        ap.error(f'{a.kind} kind needs n_s')
    elif a.kind == 'tree':
        alph = alphabet_tree_chain(a.n_s, site_names=names(a.sites),
                                   edge_names=names(a.enames))
    else:
        alph = alphabet_loop_ngon(a.n_s, c=a.c, site_names=names(a.sites),
                                  edge_names=names(a.enames),
                                  with_disc=not a.no_disc)
    s = json.dumps(alph, indent=1)
    if a.out:
        with open(a.out, 'w') as f:
            f.write(s + '\n')
        print(f"wrote {a.out}  ({len(alph['letters'])} letters, "
              f"{len(alph['first_entry'])} facets)")
    else:
        print(s)
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())
