#!/usr/bin/env python3
"""Dedup Engine B: flux-orbit canonicalizer (DFS ball, L1 key).

INDEPENDENT second dedup stack (two-engine law). Same SPEC as Engine A —
orbit canonicalization of flux pairs (f,h)
under monodromy x SL2(Z) with generator-word witnesses — but shares NO
code and NO design choices with engine_a.py:
  traversal   : depth-first stack walk (A: breadth-first queue)
  container   : sorted flat-key list + bisect  (A: hash dict)
  ball metric : L1 norm                        (A: L2^2)
  canonical key_B = (L1, Linf, lex on h||f)    (A: (L2^2, lex on f||h))
  pre-phase   : greedy steepest-descent        (A: none)
Own matrix algebra. Exact integer arithmetic only. Witnesses checkable by
witness.py. Bounded search can only UNDER-merge, never over-merge.
"""
import bisect

def _sigma_b(n):
    k = n // 2
    rows = []
    for i in range(n):
        row = [0] * n
        if i < k:
            row[k + i] = 1
        else:
            row[i - k] = -1
        rows.append(row)
    return rows

def _mvb(M, v):
    out = []
    for row in M:
        acc = 0
        for a, b in zip(row, v):
            acc += a * b
        out.append(acc)
    return tuple(out)

def _mmb(A, B):
    Bt = list(zip(*B))
    return [[sum(a * b for a, b in zip(row, col)) for col in Bt] for row in A]

def _symplb_ok(M, S):
    Mt = [list(c) for c in zip(*M)]
    return _mmb(Mt, _mmb(S, M)) == S

def _invb(M, S):
    Mt = [list(c) for c in zip(*M)]
    P = _mmb(S, _mmb(Mt, S))
    return [[-x for x in row] for row in P]

def l1(st):
    f, h = st
    return sum(abs(x) for x in f) + sum(abs(x) for x in h)

def linf(st):
    f, h = st
    return max(max(abs(x) for x in f), max(abs(x) for x in h))

def key_b(st):
    f, h = st
    return (l1(st), linf(st), h + f)

def _enc(st):
    f, h = st
    return h + f

def moves_b(gens, n):
    S = _sigma_b(n)
    tab = {}
    for i, M in enumerate(gens):
        if not _symplb_ok(M, S):
            raise ValueError("engine B: generator %d not symplectic" % i)
        tab["M%d" % i] = M
        tab["m%d" % i] = _invb(M, S)
    return tab

def _dotb(u, v):
    acc = 0
    for a, b in zip(u, v):
        acc += a * b
    return acc

def sl2_reduce_b(st):
    """Engine-B exact SL2(Z) doublet reduction (Lagrange-Gauss).

    Same commuting-action factorization as Engine A's spec, but an
    INDEPENDENT implementation: floor-division quotient with a
    round-up-only-if-strictly-past-half fix (ties round DOWN — a different
    tie law than engine A's half-up), and L1-based key_b normalization.
    Emits the realizing token word (witness-checkable)."""
    f, h = st
    toks = []
    while True:
        hh = _dotb(h, h)
        if hh == 0:
            break
        q, r = divmod(_dotb(f, h), hh)
        if 2 * r > hh:
            q += 1
        if q:
            f = tuple(a - q * b for a, b in zip(f, h))
            toks += ["t"] * q if q > 0 else ["T"] * (-q)
        if _dotb(f, f) < hh:
            f, h = tuple(-y for y in h), f
            toks.append("S")
        else:
            break
    cands = [((f, h), ()),
             ((tuple(-x for x in f), tuple(-x for x in h)), ("S", "S"))]
    if _dotb(f, f) == _dotb(h, h):
        cands.append(((tuple(-x for x in h), f), ("S",)))
        cands.append(((h, tuple(-x for x in f)), ("s",)))
    kb = min(cands, key=lambda c: key_b(c[0]))
    return kb[0], tuple(toks) + kb[1]

class _Seen(object):
    """Sorted parallel lists + bisect (engine-B container law).

    Stores (word, depth) per encoded state; a re-visit at strictly
    shallower depth updates in place and re-opens the state."""
    def __init__(self):
        self.keys = []
        self.vals = []
    def add(self, enc, val):
        i = bisect.bisect_left(self.keys, enc)
        if i < len(self.keys) and self.keys[i] == enc:
            if val[1] < self.vals[i][1]:
                self.vals[i] = val
                return True
            return False
        self.keys.insert(i, enc)
        self.vals.insert(i, val)
        return True
    def val(self, enc):
        return self.vals[bisect.bisect_left(self.keys, enc)]

def _mono_dfs(start, mono, ball, floor, depth_max, max_nodes):
    """Depth-capped DFS over monodromy moves + exact SL2 reduction."""
    cap = max(ball * l1(start), floor)
    seen = _Seen()
    seen.add(_enc(start), ((), 0))
    stack = [start]
    nodes = 0
    best, bkey = start, key_b(start)
    while stack and nodes < max_nodes:
        cur = stack.pop()
        wcur, dcur = seen.val(_enc(cur))
        if dcur >= depth_max:
            continue
        for tok, M in mono:
            raw = (_mvb(M, cur[0]), _mvb(M, cur[1]))
            nx, rtoks = sl2_reduce_b(raw)
            if l1(nx) > cap:
                continue
            if not seen.add(_enc(nx), (wcur + (tok,) + rtoks, dcur + 1)):
                continue
            nodes += 1
            stack.append(nx)
            k = key_b(nx)
            if k < bkey:
                best, bkey = nx, k
    return best, seen.val(_enc(best))[0], nodes

def canonicalize_b(st, mono, ball=8, floor=48, depth_max=9, max_nodes=30000):
    cur, word = sl2_reduce_b(st)
    tot = 0
    while True:
        b, w, nd = _mono_dfs(cur, mono, ball, floor, depth_max,
                             max_nodes - tot)
        tot += nd
        if b == cur:
            return cur, word, tot
        cur, word = b, word + w

def dedup_b(states, gens, ball=8, floor=48, depth_max=9, max_nodes=30000):
    """RAW candidates in -> witnessed canonical partition out.

    Independent Engine-B result dict, same schema as engine A."""
    n = len(states[0][0])
    tab = moves_b(gens, n)
    mono = [(t, tab[t]) for t in sorted(tab, reverse=True)]
    reps, words, nodes, trunc = [], [], 0, False
    for st in states:
        r, w, nd = canonicalize_b(st, mono, ball, floor, depth_max,
                                  max_nodes)
        if nd >= max_nodes:
            trunc = True
        reps.append(r)
        words.append(w)
        nodes += nd
    cls = {}
    for i, r in enumerate(reps):
        cls.setdefault(r, []).append(i)
    part = sorted(cls.values())
    return {"engine": "B", "partition": part, "reps": reps, "words": words,
            "nodes": nodes, "truncated": trunc,
            "N_raw": len(states), "N_dedup": len(part)}
