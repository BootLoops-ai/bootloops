#!/usr/bin/env python3
"""Referee stack: exact generator-word WITNESS checker.

Third, engine-independent implementation of the group action
(monodromy x SL2(Z) on flux pairs (f,h)) — used ONLY to verify the merge
witnesses emitted by engine_a/engine_b and to build synthetic test inputs.
Exact integer arithmetic; the budget |f^T Sigma h| is asserted invariant at
every step (violation = ArithmeticError, never a silent pass).
Token alphabet: S s T t (SL2 doublet + inverses), M<i> m<i> (monodromy gen
i and its exact symplectic inverse M^{-1} = -Sigma M^T Sigma).
"""

def sig(n):
    k = n // 2
    S = [[0] * n for _ in range(n)]
    for i in range(k):
        S[i][k + i], S[k + i][i] = 1, -1
    return S

def mv(M, v):
    return tuple(sum(a * b for a, b in zip(row, v)) for row in M)

def mm(A, B):
    n = len(A)
    return [[sum(A[i][k] * B[k][j] for k in range(n)) for j in range(n)]
            for i in range(n)]

def sympl_ok(M, S):
    MT = [list(r) for r in zip(*M)]
    return mm(MT, mm(S, M)) == S

def inv_sympl(M, S):
    MT = [list(r) for r in zip(*M)]
    A = mm(S, mm(MT, S))
    return [[-x for x in r] for r in A]

def make_table(gens, n):
    S = sig(n)
    tab = {}
    for i, M in enumerate(gens):
        if not sympl_ok(M, S):
            raise ValueError("witness: generator %d not symplectic" % i)
        tab["M%d" % i] = M
        tab["m%d" % i] = inv_sympl(M, S)
    return tab

def step(st, tok, tab):
    f, h = st
    if tok == "S":
        return (tuple(-x for x in h), f)
    if tok == "s":
        return (h, tuple(-x for x in f))
    if tok == "T":
        return (tuple(a + b for a, b in zip(f, h)), h)
    if tok == "t":
        return (tuple(a - b for a, b in zip(f, h)), h)
    M = tab[tok]
    return (mv(M, f), mv(M, h))

def budget(st, S=None):
    f, h = st
    if S is None:
        S = sig(len(f))
    return abs(sum(a * b for a, b in zip(f, mv(S, h))))

def apply(st, word, gens, tab=None):
    """Apply a token word exactly; enforce budget invariance per step."""
    n = len(st[0])
    if tab is None:
        tab = make_table(gens, n)
    S = sig(n)
    b0 = budget(st, S)
    for tok in word:
        st = step(st, tok, tab)
        if budget(st, S) != b0:
            raise ArithmeticError("budget invariant broken at token " + tok)
    return st

def verify(src, word, dst, gens, tab=None):
    """True iff word maps src -> dst exactly (the merge witness check)."""
    try:
        return apply(src, word, gens, tab) == dst
    except Exception:
        return False
