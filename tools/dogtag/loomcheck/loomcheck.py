#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision.
"""loomcheck — member of the Dogtag package (tools/dogtag/loomcheck/; the
package was formerly topology-audit): mechanical Yangian/loom/fishnet
applicability screen for position-space conformal integrals given as
x_ij^2 propagator/numerator power dictionaries (or parsed from a
Mathematica basis file in the Jiang arXiv:2607.11645 anc format).

Checks, per graph, against the published Yangian applicability conditions:
  (1) position-space conformal weight at each internal vertex:
      sum_e a_e = D (numerator x_ij^2 counted as power -1)
  (2) planarity of the FULL position-space graph including numerator edges
      (required by every Yangian theorem: CKLMZ 1708.00007, loom 2304.04654,
      planar-P-hat PRL 2505.05550)
  (3) the level-one momentum (P-hat) face condition of Loebbert-Ruenaufer-
      Stawinski 2505.05550: every position-space face k-gon must satisfy
      sum_e a_e = (k-2) D/2 (with unit denominators in D=4 only
      quadrilateral faces qualify = the fishnet/loom structure)

SCREEN-ONLY: a PASS row means "not excluded by these mechanical conditions",
it does NOT construct Yangian PDEs or SoV determinants.  A FAIL row is a
proof of non-applicability of the cited theorems as published.
Reference verdict (the nine default targets): 0/9 in any proven Yangian/SoV class
(4/9 non-planar incl. numerator edges; all planar ones fail the face
condition).  WATCH: the 2505.05550 lineage may weaken these conditions —
re-run the screen against any successor theorem before citing the NO.

INVOKE (cwd = tools/dogtag, the package directory):
  python3 -m loomcheck --selftest            # structural battery, rc 0/1
  python3 loomcheck/loomcheck.py --selftest   # the same, as a script
  python3 -m loomcheck --basis 4Loop_int_basis.m [--targets 144,163,...]
                       [--D 4] [--json RECEIPT.json]
  from loomcheck import screen_graph          # API (package dir on sys.path)
The whole-package self-test (python3 topology_audit.py --selftest) runs
this battery as its [loomcheck] leg.
"""
import argparse, json, re, sys

import networkx as nx

# default basis: the arXiv:2607.11645 ancillary basis file (download the anc
# archive from arXiv); point --basis (or LOOMCHECK_BASIS) at 4Loop_int_basis.m.
import os
DEFAULT_BASIS = os.environ.get('LOOMCHECK_BASIS', '4Loop_int_basis.m')
NINE = [144, 163, 167, 168, 169, 212, 214, 233, 260]


def load_entries(basis_path):
    with open(basis_path) as fh:
        s = fh.read().replace('\n', '').replace(' ', '')
    s = s[s.index('*)') + 2:]
    depth = 0; entries = []; cur = ''
    for ch in s:
        if ch == '{':
            depth += 1
            if depth == 2: cur = ''
            elif depth == 1: continue
        elif ch == '}':
            depth -= 1
            if depth == 1: entries.append(cur); continue
            elif depth == 0: break
        if depth >= 2: cur += ch
    return entries


def parse(entries, t):
    """Entry t (1-indexed) -> {(i,j): net power} with numerators negative."""
    e = entries[t - 1]; d = 0; pos = None
    for i, ch in enumerate(e):
        if ch == '[': d += 1
        elif ch == ']': d -= 1
        elif ch == ',' and d == 0: pos = i
    el1 = e[:pos]
    m = re.match(r'^(.*?)/\((.*)\)$', el1)
    numint, den = m.group(1), m.group(2)
    powers = {}
    for a, b, p in re.findall(r'x\[(\d),(\d)\](?:\^(\d))?', den):
        powers[(int(a), int(b))] = powers.get((int(a), int(b)), 0) + (int(p) if p else 1)
    for a, b, p in re.findall(r'x\[(\d),(\d)\](?:\^(\d))?', numint):
        powers[(int(a), int(b))] = powers.get((int(a), int(b)), 0) - (int(p) if p else 1)
    return powers


def screen_graph(powers, D=4, internal=None):
    """Screen one graph.  powers: {(i,j): net power}; internal: iterable of
    internal vertex labels (default: Jiang convention, vertices >= 5)."""
    if internal is None:
        internal = sorted({v for e in powers for v in e if v >= 5})
    w = {v: 0 for v in internal}
    for (a, b), p in powers.items():
        for v in (a, b):
            if v in w: w[v] += p
    conf = all(x == D for x in w.values())
    G = nx.MultiGraph(); G.add_edges_from(powers.keys())
    planar, emb = nx.check_planarity(G)
    numerator_edges = sorted(e for e, p in powers.items() if p < 0)
    res = {
        'conformal_internal_weight': conf,
        'planar_full_graph': bool(planar),
        'numerator_edges': [list(e) for e in numerator_edges],
        'faces': None, 'phat_face_violations': None,
    }
    if planar:
        faces = set()
        for u in emb.nodes():
            for v in emb[u]:
                f = emb.traverse_face(u, v)
                faces.add(tuple(sorted(zip(f, f[1:] + [f[0]]))))
        bad = sum(1 for f in faces
                  if sum(powers.get((min(a, b), max(a, b)), 0) for a, b in f)
                  != (len(f) - 2) * D // 2)
        res['faces'] = len(faces)
        res['phat_face_violations'] = bad
    res['yangian_class'] = bool(
        conf and planar and res['phat_face_violations'] == 0)
    return res


def selftest_cases(D=4):
    """The built-in battery as data: structural known-answer checks, no data
    files needed.  One case per screened condition, asserted from the
    published definitions — K5 is non-planar (Kuratowski) even when one edge
    is a numerator edge; the 4-point star's internal vertex has unit-power
    weight 4 = D and 5 != D with one propagator squared; a unit 4-cycle's two
    quadrilateral faces each sum to 4 = (4-2)*4/2 (0 violations) while a unit
    5-cycle's two pentagon faces sum to 5 != 6 (2 violations).  Returns a list
    of (name, got, expected) triples; a case passes iff got == expected."""
    cases = []

    def check(name, got, want):
        cases.append((name, got, want))

    # condition (2): planarity of the FULL graph including numerator edges
    k5 = {(i, j): 1 for i in range(1, 6) for j in range(i + 1, 6)}
    k5[(1, 2)] = -1                       # numerator edge still counts
    r = screen_graph(k5, D=D)
    check("K5+numerator: non-planar", r['planar_full_graph'], False)
    check("K5+numerator: numerator-edge list", r['numerator_edges'], [[1, 2]])
    check("K5+numerator: internal weight 4=D", r['conformal_internal_weight'], True)
    check("K5+numerator: excluded", r['yangian_class'], False)

    # condition (1): internal-vertex conformal weight on the 4-point star
    star = {(1, 5): 1, (2, 5): 1, (3, 5): 1, (4, 5): 1}
    check("star: unit powers weigh 4=D",
          screen_graph(star, D=D)['conformal_internal_weight'], True)
    check("star: squared propagator weighs 5!=D",
          screen_graph({**star, (1, 5): 2}, D=D)['conformal_internal_weight'], False)

    # condition (3): P-hat face sums, quadrilateral vs pentagon
    sq = screen_graph({(1, 2): 1, (2, 3): 1, (3, 4): 1, (1, 4): 1}, D=D)
    check("4-cycle: two faces", sq['faces'], 2)
    check("4-cycle: quadrilateral faces satisfy sum=(k-2)D/2",
          sq['phat_face_violations'], 0)
    pent = screen_graph({(1, 2): 1, (2, 3): 1, (3, 4): 1, (4, 5): 1, (1, 5): 1}, D=D)
    check("5-cycle: both pentagon faces violate", pent['phat_face_violations'], 2)
    check("5-cycle: excluded", pent['yangian_class'], False)
    return cases


def selftest():
    """Run the built-in battery (selftest_cases), print one PASS / FAIL line
    per check and a `SELFTEST PASS|FAIL: n/N checks passed` summary.
    Returns 0 iff every check passed, else 1."""
    cases = selftest_cases(D=4)
    fails = []
    for name, got, want in cases:
        if got == want:
            print(f"PASS  {name}")
        else:
            print(f"FAIL  {name}: got {got!r}, expected {want!r}")
            fails.append(name)
    print(f"SELFTEST {'FAIL' if fails else 'PASS'}: "
          f"{len(cases) - len(fails)}/{len(cases)} checks passed")
    return 1 if fails else 0


def main(argv=None):
    ap = argparse.ArgumentParser(prog='loomcheck',
                                 description=__doc__.split('\n')[0])
    ap.add_argument('--basis', default=DEFAULT_BASIS)
    ap.add_argument('--targets', default=','.join(map(str, NINE)),
                    help='comma-separated 1-indexed entries (default: the nine default targets)')
    ap.add_argument('--D', type=int, default=4)
    ap.add_argument('--json', metavar='OUT', help='write per-target JSON receipt')
    ap.add_argument('--selftest', action='store_true',
                    help='run the built-in battery (structural known-answer checks, '
                         'no data files needed) and exit 0/1')
    a = ap.parse_args(argv)
    if a.selftest:
        return selftest()
    entries = load_entries(a.basis)
    out = {}
    for t in (int(x) for x in a.targets.split(',')):
        r = screen_graph(parse(entries, t), D=a.D)
        out[f'I{t}'] = r
        msg = (f"I{t}: conformal(int wt={a.D})={r['conformal_internal_weight']}"
               f"  planar(full)={r['planar_full_graph']}")
        if r['planar_full_graph']:
            msg += (f"  faces={r['faces']}"
                    f"  P-hat-face-violations={r['phat_face_violations']}")
        msg += f"  yangian_class={r['yangian_class']}"
        print(msg)
    if a.json:
        with open(a.json, 'w') as fh:
            json.dump(out, fh, indent=1)
    n_in = sum(1 for r in out.values() if r['yangian_class'])
    print(f"SUMMARY: {n_in}/{len(out)} in a proven Yangian/SoV class")
    return 0


if __name__ == '__main__':
    sys.exit(main())
