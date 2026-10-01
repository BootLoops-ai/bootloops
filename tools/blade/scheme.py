"""blade.scheme -- block-scheme construction (Wolfram-free port of Blade BLSearch/Job.wl).

Builds the `<job>/<workid>/sch_{g1,nint,intid}` block
partition consumed by redg1/fitrel (src/search/template.c:15 template_integrals_init)
and ssolve (src/ssolve/block.c:80, reads multi_sch_intid) from a NUMERIC reduction
table, replacing the Wolfram package BLSearch/Job.wl.

Semantics ported (authoritative .wl line cites):

* G1/G2 split -- Job.wl:71-96 `G1G2[G, mode]`.  The nullspace of the projection
  matrix is taken with the REVERSED-COLUMN pivot convention (Job.wl:99-101
  `reverseColumn` / `evaluateSol`): columns of Transpose[RedTable[[ids]]] are
  reversed before NullSpace and the result is reversed back, so RREF pivots land
  on the LAST-listed (simplest / master) integrals and the free columns -- the
  integrals marked "solved" -- are the EARLIEST-listed (most complex) ones.
  Equivalent characterization used here (proved in `solved_positions` docstring):
  g[i] is solved  <=>  proj(g[i]) is in the span of {proj(g[j]) : j > i}.
  Getting the convention backwards marks the MASTERS as solved (the unit rows
  are always pivot candidates); the db byte-gate discriminates this.

* Block-per-sector scheme -- Job.wl:248-317 `BLGenerateJob`:
  intsToScheme (Job.wl:215), sectorDown/sectorDownCollect
  (IntegralExtension.wl:34-43, frobeniusSolve2 = AuxiliaryFunctions.wl:92-95),
  selectSubScheme (Job.wl:207-212), the sector-first / non-master-first /
  MIs-last reordering (Job.wl:267-268, note line 268 appends the FULL MIs list
  in its STORED order to every block), the skip rule (Job.wl:269), PartitionDot
  (Job.wl:281-286), PartitionRankThreshold (Job.wl:289-295), and the
  user-defined block path userDefinedG1G2 (Job.wl:181-196).

* Magic-relation loops -- upstream collectMap/backTrack/pinchMap
  (Job.wl:130-169) repeatedly finds ONE directed cycle by depth-first
  backtracking (O(N^3)) and contracts ("pinches") it until the graph is
  acyclic.  Here this is replaced by Tarjan's strongly-connected-components
  algorithm.  EQUIVALENCE: contracting any directed cycle merges only vertices
  of one SCC (a cycle is a closed walk, and all vertices on a closed walk are
  strongly connected); contraction preserves reachability between the remaining
  vertices, so iterating until no cycle remains contracts vertices u,v into the
  same super-node IFF u and v lie on some closed walk, i.e. IFF they are in the
  same SCC.  The fixed point of collectMap is therefore exactly the SCC
  condensation, independent of the order in which cycles are found (confluence).
  Self-loops cannot occur because the map stores Complement[projsec,{sec}]
  (Job.wl:269).  `collect_map_reference` is a faithful port of the upstream
  contraction loop and is property-tested against `tarjan_scc` in the phase-3
  gate driver.  ONE DELIBERATE DEVIATION, documented: upstream places a merged
  block at position First[Flatten[loop]] -- the DFS cycle entry point, which
  depends on backTrack's search order (Job.wl:117,142-151); this port places it
  at the MINIMUM member index (deterministic).  Block CONTENTS are unaffected
  (newG1G2 sorts the merged G1 and the merged integral set, Job.wl:115-116);
  only the workid assignment can permute when a run has >=2 multi-sector loops.

* SymMap dedup -- Job.wl:53-62 builds mapped->unique integral rules by comparing
  RedTable rows at ONE numeric point; a false collision (two distinct integrals
  agreeing at the test point, probability ~ nint^2/prime per pair but certain
  for e.g. int pairs equal on a slice) silently merges distinct integrals.
  IMPROVEMENT (mandatory here): `symmap_dedup` requires >= 2 points; a pair that
  matches at point 0 but differs at any confirmation point is NOT merged and is
  reported loudly in `SymMapResult.collisions` (and to stderr).

* Scheme files -- writeScheme (Job.wl:223-229) via blade.formats.Scheme
  (byte-exact writers, gated in the phase-0 roundtrip suite).  sch_intid holds
  0-based GLOBAL ids, G1 first then G2 (IntegralID/@Join[G1,G2]-1).

Integral representation: a bare tuple of propagator powers (the second argument
of BL[fam,{...}]); one family per scheme, ids are positions in `all_ints`
(0-based; upstream IntegralID is 1-based, Job.wl:46).  Masters are the LAST
`nmaster` entries of `all_ints` in their stored order (FORMATS.md, red_common).

Integral ordering (Integrals.wl:57-62, $IntegralOrdering=1 default,
SearchOptions.wl Options[BLSetReducerOptions]): the WL weight is
{-props, -dots, -rank, !ExtMappedQ[int], !MappedQ[int], Min[powers], int}.
In the BLSearch context ExtMappedQ/MappedQ have NO definition (Integrals.wl:65
"useless"), so entries 4-5 stay unevaluated and WL canonical order compares
them by their BL argument -- i.e. the tie-break after (-props,-dots,-rank) is
the canonical (lexicographic, numeric-ascending) order of the power lists, and
Min[powers]/int are never reached for distinct integrals.  Effective key
implemented: (-props, -dots, -rank, powers).  Orderings 2-4 are also provided.

Second-point discipline (this port, not upstream): every `solved_positions`
verdict that feeds a block split is recomputed at a confirmation point when a
second RedTable is supplied; a mismatch raises SchemeInstabilityError (a
reduction table whose rank structure differs between sample points cannot
define a scheme).
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

from flint import nmod_mat

try:  # package-relative when imported as blade.scheme
    from . import formats
except ImportError:  # pragma: no cover - direct script use
    import formats  # type: ignore

Powers = Tuple[int, ...]


class SchemeError(Exception):
    """Semantic failure while building a block scheme."""


class SchemeInstabilityError(SchemeError):
    """A rank/solved verdict differed between two sample points."""


# ---------------------------------------------------------------------------
# Integral helpers (BLSearch/Integrals.wl:31-45)
# ---------------------------------------------------------------------------

def props(p: Powers) -> int:
    """BLIntPropagators, Integrals.wl:31."""
    return sum(1 for a in p if a > 0)


def dots(p: Powers) -> int:
    """BLIntDots, Integrals.wl:32."""
    return sum(a - 1 for a in p if a > 1)


def rank(p: Powers) -> int:
    """BLIntRank, Integrals.wl:33."""
    return -sum(a for a in p if a < 0)


def sector(p: Powers) -> Powers:
    """BLSector, Integrals.wl:40."""
    return tuple(1 if a > 0 else 0 for a in p)


def denominator(p: Powers) -> Powers:
    """BLDenominator, Integrals.wl:41."""
    return tuple(a if a > 0 else 0 for a in p)


def weight_key(p: Powers, ordering: int = 1):
    """IntegralWeight, Integrals.wl:57-62; see module docstring for why the
    effective tie-break after the numeric entries is the powers tuple itself."""
    pr, dt, rk = props(p), dots(p), rank(p)
    if ordering == 1:
        return (-pr, -dt, -rk, p)
    if ordering == 2:
        return (-pr, -rk, -dt, p)
    if ordering == 3:
        return (-pr, -dt - rk, -dt, -rk, p)
    if ordering == 4:
        return (-pr, -dt - rk, -rk, -dt, p)
    raise SchemeError(f"unknown integral ordering {ordering}")


def sort_integrals(ints: Iterable[Powers], ordering: int = 1) -> List[Powers]:
    """BLSortIntegrals, Integrals.wl:66."""
    return sorted(ints, key=lambda p: weight_key(p, ordering))


def subsector_or_equal(s1: Powers, s2: Powers) -> bool:
    """BLSectorOrSubsectorQ, Integrals.wl:36."""
    return all(a - b <= 0 for a, b in zip(s1, s2))


def frobenius_solve(nvars: int, degree: int) -> List[Powers]:
    """FrobeniusSolve[{1,..,1}, degree] with WL's ascending-lex output order
    (AuxiliaryFunctions.wl:92-95 frobeniusSolve2 with unit weights)."""
    if nvars == 0:
        return [()] if degree == 0 else []
    if nvars == 1:
        return [(degree,)]
    out: List[Powers] = []
    for first in range(degree + 1):
        out.extend((first,) + rest for rest in frobenius_solve(nvars - 1, degree - first))
    return out


def sector_down(sect: Powers, n: int) -> List[Powers]:
    """sectorDown, IntegralExtension.wl:34-39: degree-n ISP extensions of a
    (generalized) sector; ISP slots are the ZERO positions of `sect`."""
    posi = [i for i, a in enumerate(sect) if a == 0]
    out: List[Powers] = []
    for sol in frobenius_solve(len(posi), n):
        tmp = list(sect)
        for slot, e in zip(posi, sol):
            tmp[slot] = -e
        out.append(tuple(tmp))
    return out


def sector_down_collect(sect: Powers, n: int) -> List[Powers]:
    """sectorDownCollect, IntegralExtension.wl:41-43 (single-integrand family:
    the BLExtraIntDerivDen branch is not ported).  Degrees n..0, descending."""
    out: List[Powers] = []
    for d in range(n, -1, -1):
        out.extend(sector_down(sect, d))
    return out


def _dedup(seq: Iterable) -> List:
    """DeleteDuplicates: order-preserving."""
    return list(dict.fromkeys(seq))


# ---------------------------------------------------------------------------
# Numeric reduction table
# ---------------------------------------------------------------------------

class RedTable:
    """Dense nint x nmaster projection table at one sample point of a red_*
    database (upstream: results/data {TestPrime, RedTable}, Job.wl:50).

    Structural support comes from red_posi/red_part (FORMATS.md guarantees a
    stored entry never VALUE-vanishes at a sample point, so the value pattern
    at any point equals the structural pattern)."""

    def __init__(self, db: "formats.Database", point: int = 0):
        if not 0 <= point < db.nps:
            raise SchemeError(f"point {point} outside 0..{db.nps - 1}")
        self.prime: int = db.prime
        self.nint: int = db.nint
        self.nmaster: int = db.nmaster
        self.point = point
        rows = [[0] * db.nmaster for _ in range(db.nint)]
        support: List[Set[int]] = [set() for _ in range(db.nint)]
        for m in range(db.nmaster):
            for k in range(db.part[m], db.part[m + 1]):
                gid = db.posi[k]
                rows[gid][m] = db.data[point][k] % db.prime
                support[gid].add(m)
        self.rows = rows
        self.support = support

    @classmethod
    def pair(cls, dirpath: str, points: Tuple[int, int] = (0, 1)) -> Tuple["RedTable", "RedTable"]:
        """Primary + confirmation tables from one database directory."""
        db = formats.Database.read(dirpath)
        if db.nps < 2:
            raise SchemeError("second-point discipline requires nps >= 2")
        return cls(db, points[0]), cls(db, points[1])


def solved_positions(ids: Sequence[int], rt: RedTable, reversed_pivot: bool = True) -> List[int]:
    """Positions i in `ids` whose integral is 'solved' by the reduction data.

    Faithful to evaluateSol + FirstPosition (Job.wl:74-75, 99-101):
    NullSpace[reverseColumn@Transpose@RedTable[[ids]]] puts RREF pivots on the
    reversed-leading = ORIGINAL-TRAILING columns; the free (non-pivot) columns
    are the solved integrals, and after reversing back each nullspace vector's
    first nonzero entry (normalized 1) sits at its free column, which is what
    FirstPosition[#,1] picks.  Free columns of an RREF are exactly the columns
    lying in the span of the pivot columns to their LEFT; reading the reversal
    back, that is:  ids[i] solved  <=>  proj(ids[i]) in span{proj(ids[j]): j>i}.
    Implemented as flint nmod RREF of the column-reversed matrix; the free
    columns are returned mapped back to ascending original positions (which is
    also the row order Mathematica's NullSpace induces on `solved`: free
    columns are emitted last-free-first in reversed coordinates = ascending in
    original coordinates).

    reversed_pivot=False deliberately drops the reversal (control only): pivots
    then prefer the EARLIEST-listed integrals, which mis-marks masters as
    solved; the phase-3 gate uses this as a convention-discriminating control.
    """
    k = len(ids)
    if k == 0:
        return []
    order = list(range(k - 1, -1, -1)) if reversed_pivot else list(range(k))
    entries: List[int] = []
    for m in range(rt.nmaster):
        row = [rt.rows[i][m] for i in ids]
        entries.extend(row[j] for j in order)
    mat = nmod_mat(rt.nmaster, k, entries, rt.prime)
    ref, rnk = mat.rref()
    pivots: Set[int] = set()
    for r in range(rnk):
        for c in range(k):
            if ref[r, c] != 0:
                pivots.add(c)
                break
    free = [c for c in range(k) if c not in pivots]
    if reversed_pivot:
        return sorted(k - 1 - f for f in free)
    return sorted(free)


def _confirm_solved(ids: Sequence[int], rt: RedTable, rt2: Optional[RedTable],
                    context: str) -> List[int]:
    sol = solved_positions(ids, rt)
    if rt2 is not None:
        sol2 = solved_positions(ids, rt2)
        if sol2 != sol:
            raise SchemeInstabilityError(
                f"{context}: solved set differs between point {rt.point} and "
                f"point {rt2.point}: {sol} vs {sol2}")
    return sol


# ---------------------------------------------------------------------------
# G1G2 (Job.wl:71-96)
# ---------------------------------------------------------------------------

def g1g2(g_ids: Sequence[int], mode: str, all_ints: Sequence[Powers], rt: RedTable,
         rt_confirm: Optional[RedTable] = None) -> Tuple[List[int], List[int]]:
    """Port of G1G2[G, mode] (Job.wl:71-96).  `g_ids` should be ordered
    (complex first, masters last) exactly as upstream expects G sorted."""
    if not g_ids:
        return [], []
    sol = _confirm_solved(g_ids, rt, rt_confirm, "g1g2")
    solved = [g_ids[i] for i in sol]
    top = all_ints[g_ids[0]]
    if mode == "sector":
        g1 = [i for i in solved if sector(all_ints[i]) == sector(top)]
    elif mode == "rank":
        g1 = [i for i in solved
              if sector(all_ints[i]) == sector(top) and rank(all_ints[i]) == rank(top)]
    elif mode == "dot":
        g1 = [i for i in solved
              if sector(all_ints[i]) == sector(top) and dots(all_ints[i]) == dots(top)]
    elif mode == "all":
        g1 = list(solved)
    else:
        raise SchemeError(f"undefined mode for G1G2: {mode!r}")  # Job.wl:91
    g1set = set(g1)
    g2 = [i for i in g_ids if i not in g1set]
    return g1, g2


# ---------------------------------------------------------------------------
# Magic-relation loops: Tarjan SCC (replaces collectMap/backTrack, Job.wl:130-169)
# ---------------------------------------------------------------------------

def tarjan_scc(nnodes: int, edges: Dict[int, Sequence[int]]) -> List[List[int]]:
    """SCCs with >= 2 members, each sorted ascending, ordered by min member.

    Deterministic iterative Tarjan; equivalence with upstream's
    contract-cycles-until-acyclic loop is argued in the module docstring and
    property-tested against `collect_map_reference` in the gate driver.
    Self-loops are ignored (upstream cannot produce them, Job.wl:269)."""
    index = [0] * nnodes
    low = [0] * nnodes
    on = [False] * nnodes
    idx = [1]
    stack: List[int] = []
    comps: List[List[int]] = []

    def strongconnect(v0: int) -> None:
        work = [(v0, 0)]
        while work:
            v, pi = work.pop()
            if pi == 0:
                index[v] = low[v] = idx[0]
                idx[0] += 1
                stack.append(v)
                on[v] = True
            recurse = False
            nbrs = edges.get(v, ())
            for j in range(pi, len(nbrs)):
                w = nbrs[j]
                if w == v:
                    continue
                if index[w] == 0:
                    work.append((v, j + 1))
                    work.append((w, 0))
                    recurse = True
                    break
                if on[w]:
                    low[v] = min(low[v], index[w])
            if not recurse:
                if low[v] == index[v]:
                    comp = []
                    while True:
                        w = stack.pop()
                        on[w] = False
                        comp.append(w)
                        if w == v:
                            break
                    comps.append(sorted(comp))
                if work:
                    parent = work[-1][0]
                    low[parent] = min(low[parent], low[v])

    for v in range(nnodes):
        if index[v] == 0:
            strongconnect(v)
    multi = [c for c in comps if len(c) > 1]
    multi.sort(key=lambda c: c[0])
    return multi


def collect_map_reference(nnodes: int, edges: Dict[int, Sequence[int]]) -> List[List[int]]:
    """Faithful (python) port of the UPSTREAM contraction loop
    collectMap/backTrack/pinchMap (Job.wl:130-169), used only as a test ORACLE
    for tarjan_scc.  Nodes are kept as frozensets of original indices; a DFS
    from a virtual root (Job.wl:131) finds one cycle, which is pinched into a
    combined node (Job.wl:155-161); repeat until acyclic.  Returns the combined
    nodes with >= 2 members, sorted like tarjan_scc."""
    node_of = {i: frozenset([i]) for i in range(nnodes)}
    adj: Dict[frozenset, List[frozenset]] = {
        frozenset([i]): [frozenset([j]) for j in edges.get(i, ()) if j != i]
        for i in range(nnodes)}

    def find_cycle() -> Optional[List[frozenset]]:
        # backTrack (Job.wl:142-151): DFS over nodes in key order from a root
        # connected to every node; returns the path with the repeated node.
        keys = list(adj.keys())

        def dfs(path: List[frozenset], point: Optional[frozenset]) -> Optional[List[frozenset]]:
            targets = keys if point is None else adj.get(point, [])
            for ii in targets:
                values = adj.get(ii, [])
                if not values:
                    continue
                hit = [p for p in path if p in values]
                if hit:
                    return path + [ii, hit[0]]
                res = dfs(path + [ii], ii)
                if res is not None:
                    return res
            return None

        return dfs([], None)

    while True:
        res = find_cycle()
        if res is None:
            break
        # pinchMap (Job.wl:155-161)
        last = res[-1]
        start = res.index(last)
        com_members = _dedup(res[start:])
        com = frozenset().union(*com_members)
        new_adj: Dict[frozenset, List[frozenset]] = {}
        out: List[frozenset] = []
        for k, vals in adj.items():
            if k in com_members:
                out.extend(vals)
            else:
                new_adj[k] = _dedup(com if v in com_members else v for v in vals)
        new_adj[com] = [v if v not in com_members else com for v in out]
        new_adj[com] = _dedup(v for v in new_adj[com] if v != com)
        adj = new_adj

    multi = [sorted(k) for k in adj if len(k) > 1]
    multi.sort(key=lambda c: c[0])
    return multi


# ---------------------------------------------------------------------------
# SymMap dedup with mandatory second-point confirmation (Job.wl:53-62 improved)
# ---------------------------------------------------------------------------

@dataclass
class SymMapResult:
    rules: Dict[int, int]                 # mapped id -> unique id (first match, like x/.SymMap)
    collisions: List[Tuple[int, int]]     # matched at point 0 but NOT at a confirm point
    symflag: bool                         # ContainsAll[Keys[SymMap], mapped] (Job.wl:62)


def symmap_dedup(mapped_ids: Sequence[int], unique_ids: Sequence[int],
                 rts: Sequence[RedTable]) -> SymMapResult:
    """Port of the SymMap construction (Job.wl:60-62) with a MANDATORY
    second-point confirmation.  Upstream compares RedTable rows at ONE point;
    a false collision silently merges distinct integrals into one unique
    sector.  Here >= 2 RedTables (distinct points of the same database) are
    required; any pair equal at rts[0] but different at a later point is
    flagged LOUDLY (stderr + SymMapResult.collisions) and NOT merged."""
    if len(rts) < 2:
        raise SchemeError("symmap_dedup: >= 2 sample points are mandatory")
    rules: Dict[int, int] = {}
    collisions: List[Tuple[int, int]] = []
    for i in mapped_ids:
        for j in unique_ids:
            if rts[0].rows[i] != rts[0].rows[j]:
                continue
            if all(rt.rows[i] == rt.rows[j] for rt in rts[1:]):
                rules.setdefault(i, j)   # first matching rule wins, like /.SymMap
            else:
                collisions.append((i, j))
                print(f"scheme.symmap_dedup: LOUD COLLISION ids {i}->{j}: rows equal "
                      f"at point {rts[0].point} but differ at a confirmation point; "
                      f"NOT merged", file=sys.stderr, flush=True)
    symflag = set(mapped_ids) <= set(rules)
    return SymMapResult(rules, collisions, symflag)


# ---------------------------------------------------------------------------
# Scheme construction (BLGenerateJob, Job.wl:248-317)
# ---------------------------------------------------------------------------

@dataclass
class Block:
    g1: List[int]
    g2: List[int]

    @property
    def ids(self) -> List[int]:
        return self.g1 + self.g2


@dataclass
class SchemeResult:
    blocks: List[Block]
    allsec: List[Powers]
    magic_map: Dict[int, List[int]]
    loops: List[List[int]]
    nmaster: int
    all_ints: List[Powers]
    notes: List[str] = field(default_factory=list)

    def sizes(self) -> List[Tuple[int, int]]:
        """(len G1, nint) per block."""
        return [(len(b.g1), len(b.g1) + len(b.g2)) for b in self.blocks]

    def summary(self) -> str:
        alli = len(_dedup(i for b in self.blocks for i in b.ids))
        red = sum(len(b.g1) for b in self.blocks)
        # Job.wl:303-305 message
        return (f"block number: {len(self.blocks)}, all integrals: {alli}, "
                f"reduced integrals: {red}, masters: {alli - red}; "
                f"size of blocks -> {[len(b.g1) for b in self.blocks]}")


def _scheme_entries(all_ints: Sequence[Powers], id_of: Dict[Powers, int],
                    ordering: int) -> List[Tuple[Powers, int]]:
    """intsToScheme (Job.wl:215): one (denominator, rank-of-heaviest) entry per
    denominator group of the weight-sorted family integrals."""
    fam = sort_integrals(all_ints, ordering)
    entries: List[Tuple[Powers, int]] = []
    seen: Set[Powers] = set()
    for p in fam:
        d = denominator(p)
        if d not in seen:
            seen.add(d)
            entries.append((d, rank(p)))
    return entries


def _select_sub_scheme(scheme: Sequence[Tuple[Powers, int]], sec: Powers,
                       sector_rules: Dict[Powers, Powers],
                       unique_subsector: bool, symflag: bool) -> List[Tuple[Powers, int]]:
    """selectSubScheme (Job.wl:207-212).  `sec` is a 0/1 sector, so
    BLLowerDenom[sec] = all 0/1 subsector masks; sub2 = proper subsectors,
    optionally SymMap-replaced (sector-level keys only)."""
    npr = props(sec)
    ones = [i for i, a in enumerate(sec) if a > 0]
    sub2: List[Powers] = []
    for mask in range(1 << len(ones)):
        t = [0] * len(sec)
        nset = 0
        for b, pos in enumerate(ones):
            if mask >> b & 1:
                t[pos] = 1
                nset += 1
        if nset < npr:
            sub2.append(tuple(t))
    if unique_subsector and symflag:
        sub2 = _dedup(sector_rules.get(s, s) for s in sub2)
    keep = set(sub2) | {sec}   # Complement[sub,sub2] == {sec} at equal props
    return [e for e in scheme if sector(e[0]) in keep]


def _family_sorted_known(gen: Iterable[Powers], id_of: Dict[Powers, int],
                         ordering: int) -> List[Powers]:
    """FamilyIntegrals (Job.wl:320): keep known integrals, weight-sort."""
    return sort_integrals((p for p in _dedup(gen) if p in id_of), ordering)


def build_scheme(all_ints: Sequence[Powers], nmaster: int,
                 rt: RedTable, rt_confirm: Optional[RedTable] = None, *,
                 usints: Sequence[Tuple[int, Sequence[Powers]]] = (),
                 exints: Optional[Dict[int, Sequence[Powers]]] = None,
                 symmap: Optional[SymMapResult] = None,
                 partition_dot: bool = True,
                 partition_rank_threshold: int = 100,
                 unique_subsector: bool = True,
                 integral_ordering: int = 1) -> SchemeResult:
    """Full port of BLGenerateJob (Job.wl:248-317).

    all_ints  : global integral list; index == 0-based global id == database id;
                the LAST `nmaster` entries are the masters in database order.
    rt        : RedTable at the primary point; rt_confirm at a second point
                (recommended; every solved verdict is re-checked there).
    usints    : user-defined blocks, list of (global id of the us-integral,
                [BL integrals occurring in its definition]) -- userDefinedG1G2.
    exints    : preferred-master definitions, master GLOBAL id -> definition
                integrals (MIs[[proj]]/.EXInts, Job.wl:263); default none.
    symmap    : result of symmap_dedup (SymFlag semantics, Job.wl:62); default
                empty rules with symflag=True (ContainsAll[{}, {}] is True, so
                upstream's $UniqueSubsectorQ branch is a no-op then).

    Defaults mirror Options[BLSetSchemeOptions] (SearchOptions.wl:
    PartitionDot->True, PartitionRankThreshold->100, UniqueSubsectorQ->True)
    and IntegralOrdering->1.
    """
    all_ints = [tuple(p) for p in all_ints]
    nint = len(all_ints)
    if len(set(all_ints)) != nint:
        raise SchemeError("duplicate integrals in all_ints")
    if rt.nint != nint or rt.nmaster != nmaster:
        raise SchemeError("RedTable shape disagrees with all_ints/nmaster")
    id_of = {p: i for i, p in enumerate(all_ints)}
    mis = list(range(nint - nmaster, nint))
    mis_set = set(mis)
    exints = exints or {}
    if symmap is None:
        symmap = SymMapResult(rules={}, collisions=[], symflag=True)
    sector_rules = {sector(all_ints[k]): sector(all_ints[v])
                    for k, v in symmap.rules.items()
                    if dots(all_ints[k]) == 0 and rank(all_ints[k]) == 0}

    notes: List[str] = []
    scheme = _scheme_entries(all_ints, id_of, integral_ordering)
    allsec = _dedup(sector(d) for d, _ in scheme)

    # --- per-sector blocks + magic map (Job.wl:258-270) ---
    blocks: List[Optional[Block]] = []
    magic_map: Dict[int, List[int]] = {}
    secindex = {s: i for i, s in enumerate(allsec)}
    for sec in allsec:
        # projections of this sector's scheme integrals (Job.wl:261-262);
        # structural support is used (FORMATS.md: stored entries never
        # value-vanish at sample points, so this equals the value pattern).
        gen = []
        for d, r in scheme:
            if sector(d) == sec:
                gen.extend(sector_down_collect(d, r))
        proj: Set[int] = set()
        for p in _dedup(gen):
            gid = id_of.get(p)
            if gid is not None:
                proj |= rt.support[gid]
        # projsec (Job.wl:263): sectors of the (EXInts-replaced) masters at the
        # SAME propagator count as sec.
        projsec: List[Powers] = []
        for m in sorted(proj):
            mid = mis[m]
            defs = exints.get(mid, [all_ints[mid]])
            projsec.extend(sector(q) for q in defs)
        projsec = [s for s in _dedup(projsec) if props(s) == props(sec)]
        # sch (Job.wl:264): Union sorts canonically = ascending lex on tuples.
        secs_for_sub = sorted(set(projsec) | {sec})
        sch: List[Tuple[Powers, int]] = []
        for s in secs_for_sub:
            sch.extend(_select_sub_scheme(scheme, s, sector_rules,
                                          unique_subsector, symmap.symflag))
        sch = _dedup(sch)
        # G (Job.wl:266-268): known+sorted, sector-sec first, then non-masters
        # first with the FULL MIs list appended in stored order.
        gen2: List[Powers] = []
        for d, r in sch:
            gen2.extend(sector_down_collect(d, r))
        fam = _family_sorted_known(gen2, id_of, integral_ordering)
        fam = ([p for p in fam if sector(p) == sec]
               + [p for p in fam if sector(p) != sec])
        g_ids = [id_of[p] for p in fam if id_of[p] not in mis_set] + mis
        # skip rule (Job.wl:269)
        if not g_ids or g_ids[0] in mis_set or sector(all_ints[g_ids[0]]) != sec:
            blocks.append(None)
        else:
            g1, g2 = g1g2(g_ids, "sector", all_ints, rt, rt_confirm)
            if not g1:
                # Upstream would keep {{},G2} and later mangle it via the
                # /.{}->Nothing replacement (Job.wl:278) -- a latent bug path;
                # a sector with an empty G1 has nothing reducible, drop it.
                notes.append(f"sector {sec}: empty G1, block dropped "
                             f"(upstream {{}}->Nothing edge, Job.wl:278)")
                blocks.append(None)
            else:
                blocks.append(Block(g1, g2))
        magic_map[secindex[sec]] = sorted(secindex[s] for s in projsec if s != sec)

    # --- magic loops (Job.wl:272-278) via Tarjan SCC ---
    loops = tarjan_scc(len(allsec), magic_map)
    drop: Set[int] = set()
    for comp in loops:                       # newG1G2 (Job.wl:110-123)
        member_blocks = [blocks[i] for i in comp if blocks[i] is not None]
        ints_union: List[int] = []
        g1_union: List[int] = []
        for b in member_blocks:
            g1_union.extend(b.g1)
            ints_union.extend(b.ids)
        ints_sorted = sort_integrals([all_ints[i] for i in _dedup(ints_union)],
                                     integral_ordering)
        g1_sorted = sort_integrals([all_ints[i] for i in _dedup(g1_union)],
                                   integral_ordering)
        g1_ids = [id_of[p] for p in g1_sorted]
        g1_set = set(g1_ids)
        g2_ids = [id_of[p] for p in ints_sorted if id_of[p] not in g1_set]
        home = comp[0]   # min member (deterministic; see docstring DEVIATION)
        blocks[home] = Block(g1_ids, g2_ids) if g1_ids else None
        drop |= set(comp[1:])
    merged = [b for i, b in enumerate(blocks) if i not in drop and b is not None]

    # --- PartitionDot (Job.wl:281-286) ---
    if partition_dot:
        out: List[Block] = []
        for b in merged:
            maxd = max(dots(all_ints[i]) for i in b.g1)
            for jj in range(0, maxd + 1):
                cur = [i for i in b.g1 if dots(all_ints[i]) == jj]
                if not cur:
                    continue
                ge = {i for i in b.g1 if dots(all_ints[i]) >= jj}
                out.append(Block(cur, [i for i in b.ids if i not in ge]))
        merged = out

    # --- PartitionRankThreshold (Job.wl:289-295) ---
    out2: List[Block] = []
    for b in merged:
        if len(b.g1) <= partition_rank_threshold:
            out2.append(b)
            continue
        ranks = _dedup(rank(all_ints[i]) for i in b.g1)
        if len(ranks) == 1:
            out2.append(b)
            continue
        # upstream: GroupBy[...,BLIntRank]//ReverseSort, First = s_max group.
        # ReverseSort on an Association sorts by VALUES (canonically: length
        # first), which is ambiguous when group sizes differ; the documented
        # INTENT ("rank equal to s_max", Job.wl:11,288) is implemented: the
        # top group is the maximal rank.  Noted in SchemeResult.notes.
        top_rank = max(ranks)
        top = [i for i in b.g1 if rank(all_ints[i]) == top_rank]
        rest = [i for i in b.g1 if rank(all_ints[i]) != top_rank]
        out2.append(Block(top, rest + b.g2))
        out2.append(Block(rest, list(b.g2)))
        notes.append(f"PartitionRankThreshold split applied (intent reading of "
                     f"Job.wl:291 ReverseSort; see scheme.py)")
    merged = out2

    # --- user-defined blocks (Job.wl:181-196, 298) ---
    for usblock in _user_defined_blocks(usints, scheme, all_ints, id_of, mis,
                                        mis_set, sector_rules, unique_subsector,
                                        symmap.symflag, integral_ordering):
        merged.append(usblock)

    # --- final consistency check (Job.wl:300-312) ---
    alli = _dedup(i for b in merged for i in b.ids)
    reduced = sum(len(b.g1) for b in merged)
    if alli and len(alli) - reduced != nmaster:
        raise SchemeError(
            f"block-triangular masters ({len(alli) - reduced}) are not the same "
            f"as database masters ({nmaster})")   # Job.wl:311

    return SchemeResult(blocks=merged, allsec=allsec, magic_map=magic_map,
                        loops=loops, nmaster=nmaster, all_ints=list(all_ints),
                        notes=notes)


def _user_defined_blocks(usints, scheme, all_ints, id_of, mis, mis_set,
                         sector_rules, unique_subsector, symflag,
                         ordering) -> List[Block]:
    """userDefinedG1G2 (Job.wl:181-196).  usints: [(us global id, def ints)]."""
    if not usints:
        return []
    us = [(uid, _dedup(sector(p) for p in defs)) for uid, defs in usints]
    # GatherBy key (Job.wl:183): the max-propagator sectors of the definition.
    groups: Dict[Tuple[Powers, ...], List[Tuple[int, List[Powers]]]] = {}
    order: List[Tuple[Powers, ...]] = []
    for uid, secs in us:
        mx = max(props(s) for s in secs)
        key = tuple(s for s in secs if props(s) == mx)
        if key not in groups:
            groups[key] = []
            order.append(key)
        groups[key].append((uid, secs))
    blocks: List[Block] = []
    for key in order:
        members = groups[key]
        projsec = _dedup(s for _, secs in members for s in secs)
        sch: List[Tuple[Powers, int]] = []
        for s in projsec:            # Job.wl:190: plain Join, NO Union sort
            sch.extend(_select_sub_scheme(scheme, s, sector_rules,
                                          unique_subsector, symflag))
        sch = _dedup(sch)
        gen: List[Powers] = []
        for d, r in sch:
            gen.extend(sector_down_collect(d, r))
        fam = _family_sorted_known(gen, id_of, ordering)
        g2 = [id_of[p] for p in fam if id_of[p] not in mis_set] + mis
        blocks.append(Block([uid for uid, _ in members], g2))
    return blocks


# ---------------------------------------------------------------------------
# Single user-target block (degenerate scheme used by the phase-0/2 banks)
# ---------------------------------------------------------------------------

def single_block_scheme(all_ints: Sequence[Powers], nmaster: int, rt: RedTable,
                        rt_confirm: Optional[RedTable] = None, *,
                        order: str = "given",
                        integral_ordering: int = 1) -> SchemeResult:
    """One block reducing ALL non-master integrals at once: G1G2 with mode
    "all" (Job.wl:87-88) on G = non-masters + MIs.  This is the degenerate
    scheme the phase-0 db bank and the phase-2 dbox bank were fabricated with
    (a user-target amplitude reduced in one work; cf. userDefinedG1G2 which
    would put the single us-symbol in G1 instead -- see NOTES.md).

    order="given"  : non-masters in all_ints order (the fabricators' order);
    order="weight" : non-masters weight-sorted (what FamilyIntegrals would do,
                     Job.wl:320).  For the db family the two coincide; for the
                     dbox bank only "given" matches (targets were emitted in
                     run.wl order, which is the REVERSE of the weight order)."""
    all_ints = [tuple(p) for p in all_ints]
    nint = len(all_ints)
    mis = list(range(nint - nmaster, nint))
    non = list(range(nint - nmaster))
    if order == "weight":
        non = [i for _, i in sorted(((weight_key(all_ints[i], integral_ordering), i)
                                     for i in non))]
    elif order != "given":
        raise SchemeError(f"unknown order {order!r}")
    g_ids = non + mis
    g1, g2 = g1g2(g_ids, "all", all_ints, rt, rt_confirm)
    return SchemeResult(blocks=[Block(g1, g2)], allsec=[], magic_map={},
                        loops=[], nmaster=nmaster, all_ints=list(all_ints))


# ---------------------------------------------------------------------------
# Emission (writeScheme, Job.wl:223-229)
# ---------------------------------------------------------------------------

def write_scheme(result: SchemeResult, jobdir: str) -> List[str]:
    """Write <jobdir>/<workid>/sch_{g1,nint,intid}+multi_sch_intid for every
    block, workid = 0..nblocks-1 (Job.wl:314 writes block i to workid i-1).
    Byte-exact via blade.formats.Scheme.  Returns the work dirs."""
    import os
    dirs: List[str] = []
    for wid, b in enumerate(result.blocks):
        d = os.path.join(jobdir, str(wid))
        sch = formats.Scheme(g1=len(b.g1), nint=len(b.g1) + len(b.g2),
                             intid=b.ids)
        sch.write(d)
        dirs.append(d)
    return dirs
