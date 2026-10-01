#!/usr/bin/env python3
"""Dedup Engine A: flux-orbit canonicalizer (BFS ball, L2 key).

Flux-orbit canonicalizer for IIB pairs (f,h) under monodromy x SL2(Z):
  monodromy gen M in Sp(2k,Z):  (f,h) -> (Mf, Mh)
  SL2(Z) doublet:  S:(f,h)->(-h,f)   T:(f,h)->(f+h,h)   (+ inverses s,t)
Both preserve the budget |f^T Sigma h| exactly.

Algorithm: BFS-ball fixpoint canonicalization. Visited = hash dict;
canonical key_A = (L2^2 norm, lex on f||h); ball cap = ball * L2^2(center),
iterated to a fixpoint whenever a lower-key state is found.
Every candidate carries a generator-word WITNESS (start -> canonical rep),
checkable by witness.py; merges are therefore sound by construction —
bounded search can only UNDER-merge (split), never over-merge.

TWO-ENGINE LAW: this file shares NO code with engine_b.py — different
traversal (BFS vs DFS), container (dict vs sorted list + bisect), ball
metric (L2 vs L1), key (f-major vs h-major), no greedy pre-phase. Exact
integer arithmetic only.
"""
from collections import deque

def sigma(n):
    k = n // 2
    S = [[0] * n for _ in range(n)]
    for i in range(k):
        S[i][k + i] = 1
        S[k + i][i] = -1
    return S

def _mv(M, v):
    return tuple(sum(r[j] * v[j] for j in range(len(v))) for r in M)

def _mm(A, B):
    n = len(A)
    return [[sum(A[i][k] * B[k][j] for k in range(n)) for j in range(n)]
            for i in range(n)]

def _sympl_ok(M, S):
    MT = [[M[j][i] for j in range(len(M))] for i in range(len(M))]
    return _mm(MT, _mm(S, M)) == S

def _inv_sympl(M, S):
    MT = [[M[j][i] for j in range(len(M))] for i in range(len(M))]
    A = _mm(S, _mm(MT, S))
    return [[-x for x in r] for r in A]

def n2(st):
    f, h = st
    return sum(x * x for x in f) + sum(x * x for x in h)

def key_a(st):
    f, h = st
    return (n2(st), f + h)

def build_moves(gens, n):
    S = sigma(n)
    mbt = {}
    for i, M in enumerate(gens):
        if not _sympl_ok(M, S):
            raise ValueError("engine A: generator %d not symplectic" % i)
        mbt["M%d" % i] = M
        mbt["m%d" % i] = _inv_sympl(M, S)
    return mbt

def _dot(u, v):
    return sum(a * b for a, b in zip(u, v))

def sl2_reduce(st):
    """Exact Lagrange-Gauss SL2(Z) reduction of the doublet (f,h).

    The SL2(Z) action (column ops on [f h]) COMMUTES with monodromy (left
    multiplication), so it is reduced EXACTLY here — no search — and the
    orbit search runs over monodromy words only. Emits the token word
    realizing the reduction (witness-checkable by witness.py).
    Rounding law (engine A): nearest integer, ties round HALF-UP."""
    f, h = st
    toks = []
    while True:
        hh = _dot(h, h)
        if hh == 0:
            break
        k = (2 * _dot(f, h) + hh) // (2 * hh)
        if k:
            f = tuple(a - k * b for a, b in zip(f, h))
            toks += ["t"] * k if k > 0 else ["T"] * (-k)
        if _dot(f, f) < hh:
            f, h = tuple(-x for x in h), f
            toks.append("S")
        else:
            break
    cands = [((f, h), ()),
             ((tuple(-x for x in f), tuple(-x for x in h)), ("S", "S"))]
    if _dot(f, f) == _dot(h, h):
        cands.append(((tuple(-x for x in h), f), ("S",)))
        cands.append(((h, tuple(-x for x in f)), ("s",)))
    best = min(cands, key=lambda c: key_a(c[0]))
    return best[0], tuple(toks) + best[1]

def _mono_ball(start, mono, ball, floor, depth_max, max_nodes):
    """BFS over monodromy moves, each composed with exact SL2 reduction."""
    cap = max(ball * n2(start), floor)
    words = {start: ((), 0)}
    q = deque([start])
    nodes = 0
    best, bkey = start, key_a(start)
    while q and nodes < max_nodes:
        cur = q.popleft()
        wcur, dcur = words[cur]
        if dcur >= depth_max:
            continue
        for tok, M in mono:
            raw = (_mv(M, cur[0]), _mv(M, cur[1]))
            nx, rtoks = sl2_reduce(raw)
            if nx in words or n2(nx) > cap:
                continue
            words[nx] = (wcur + (tok,) + rtoks, dcur + 1)
            nodes += 1
            q.append(nx)
            k = key_a(nx)
            if k < bkey:
                best, bkey = nx, k
    return best, words[best][0], nodes

def canonicalize(st, mono, ball=25, floor=64, depth_max=9, max_nodes=30000):
    cur, word = sl2_reduce(st)
    tot = 0
    while True:
        b, w, nd = _mono_ball(cur, mono, ball, floor, depth_max,
                              max_nodes - tot)
        tot += nd
        if b == cur:
            return cur, word, tot
        cur, word = b, word + w

def dedup(states, gens, ball=25, floor=64, depth_max=9, max_nodes=30000):
    """RAW candidates in -> witnessed canonical partition out."""
    n = len(states[0][0])
    mbt = build_moves(gens, n)
    mono = [(t, mbt[t]) for t in sorted(mbt)]
    reps, words, nodes, trunc = [], [], 0, False
    for st in states:
        r, w, nd = canonicalize(st, mono, ball, floor, depth_max, max_nodes)
        if nd >= max_nodes:
            trunc = True
        reps.append(r)
        words.append(w)
        nodes += nd
    cls = {}
    for i, r in enumerate(reps):
        cls.setdefault(r, []).append(i)
    part = sorted(cls.values())
    return {"engine": "A", "partition": part, "reps": reps, "words": words,
            "nodes": nodes, "truncated": trunc,
            "N_raw": len(states), "N_dedup": len(part)}
