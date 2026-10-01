#!/usr/bin/env python3
r"""
sparse_cascade.py — reusable GF(p) symbol-cascade engine.

Replaces the dense `kron(WORD[w-1], I_A)` pattern (which blows up to A^w
columns) with a *sparse-support* word matrix: we only ever store the columns
(word indices) that are actually nonzero in the current basis.  All linear
algebra is over GF(p) using float64-BLAS-backed matmul with periodic mod-p
reduction; the chunk length is clamped so every partial sum stays ≤ 2^53
(exactly representable in float64), which admits any p with
(p-1)² + (p-1) ≤ 2^53 (p up to 94 906 266); larger p raises ValueError.
Boundary proof: `python3 -m maxcut.sparse_cascade` runs the exactness probe.

Core objects
------------
SparseWord
    A basis of N vectors expanded on a sparse set of word-columns:
    `mat` is (N, ncol) int64, `cols` is the length-ncol array of global word
    indices (big-endian base-A encoding: idx = d0·A^{w-1}+…+d_{w-1}).

incremental_nullspace(C, n, p)
    Nullspace of an (nr × n) constraint matrix over GF(p), processed in
    row-batches so nr ≫ n costs O(nr·rank·n / batch) BLAS flops, never
    materialising more than O((batch+rank)·n) at once.

Typical use:
    W1 = SparseWord(eye, first_entry_letters, p)
    for w in 2..wmax:
        Bw = nullspace_mod(integ_rows, p)          # (Nw, Nm1, A) recursive basis
        Ww = W_{w-1}.extend(Bw, A)                 # no kron — chunked over A
    project_sparse(Ww, constraint_row_iter)        # crossing / LE / LOCKFE / drop
"""
import numpy as np

# ---------------------------------------------------------------- GF(p) dense
def mod(M, p):
    return np.rint(M).astype(np.int64) % p

def matmul_mod(X, Y, p, ch=2048):
    """(X @ Y) mod p via float64 BLAS, chunked along the inner dim so that
    every partial sum is an exactly-representable float64 integer.

    After each chunk the accumulator is reduced to [0, p), so exactness needs
    (p-1) + ch·(p-1)² ≤ 2^53; `ch` is clamped to that bound (ch=2048 is kept
    up to p ≤ 2 097 152).  Admits any p ≤ 94 906 266 — beyond that even a
    single product overflows 2^53 and a ValueError is raised."""
    X = np.ascontiguousarray(X); Y = np.ascontiguousarray(Y)
    m, k = X.shape; _, n = Y.shape
    pm1 = int(p) - 1
    if pm1 > 0:
        cap = ((1 << 53) - pm1) // (pm1 * pm1)
        if cap < 1:
            raise ValueError(f"matmul_mod: p={p} exceeds the float64-exact "
                             "envelope ((p-1)^2 + (p-1) > 2^53)")
        ch = min(int(ch), cap)
    Xf = X.astype(np.float64); Yf = Y.astype(np.float64)
    out = np.zeros((m, n), dtype=np.float64)
    for i in range(0, k, ch):
        out = (out + Xf[:, i:i+ch] @ Yf[i:i+ch, :]) % p
    return mod(out, p)

def matmul_max_p():
    """Largest modulus matmul_mod admits: max p with (p-1)² + (p-1) ≤ 2^53."""
    from math import isqrt
    pm1 = isqrt(1 << 53)
    while pm1 * pm1 + pm1 > (1 << 53):
        pm1 -= 1
    return pm1 + 1

def matmul_boundary_probe():
    """Prove matmul_mod exactness at the envelope boundary.  Returns (ok, checks).

    Worst-case legs: all entries p-1 (the largest possible addends), where the
    exact answer is k·(p-1)² ≡ k (mod p) — run at p=2,999,999 (a modulus the
    pre-clamp code got wrong at ch=2048) and at the largest admissible p.
    Cross-check leg: random entries in [p/2, p) against the exact big-integer
    product.  Guard leg: p above the envelope must raise ValueError."""
    checks = []
    for p in (2_999_999, matmul_max_p()):
        k = 4096
        X = np.full((4, k), p - 1, dtype=np.int64)
        Y = np.full((k, 4), p - 1, dtype=np.int64)
        got = matmul_mod(X, Y, p)
        checks.append((f"allmax p={p} k={k}", bool(np.all(got == k % p))))
    rng = np.random.default_rng(20260901)
    p = 2_999_999
    X = rng.integers(p // 2, p, size=(8, 2048)).astype(np.int64)
    Y = rng.integers(p // 2, p, size=(2048, 8)).astype(np.int64)
    exact = ((X.astype(object) @ Y.astype(object)) % p).astype(np.int64)
    checks.append((f"random-vs-exact p={p}",
                   bool(np.all(matmul_mod(X, Y, p) == exact))))
    try:
        matmul_mod(np.ones((1, 1), dtype=np.int64),
                   np.ones((1, 1), dtype=np.int64), matmul_max_p() + 1)
        checks.append(("overlarge-p-refused", False))
    except ValueError:
        checks.append(("overlarge-p-refused", True))
    return all(ok for _, ok in checks), checks

def rref_mod(Min, p):
    """Reduced row-echelon form over GF(p).  Returns (R[:rank], pivot_cols)."""
    M = (np.asarray(Min, dtype=np.int64) % p).astype(np.float64)
    m, n = M.shape; piv = []; r = 0
    for c in range(n):
        if r >= m: break
        col = mod(M[r:, c], p); nz = np.nonzero(col)[0]
        if len(nz) == 0: continue
        pr = r + int(nz[0])
        if pr != r: M[[r, pr]] = M[[pr, r]]
        ip = pow(int(mod(M[r, c], p)), p - 2, p)
        M[r] = mod(M[r] * ip, p).astype(np.float64)
        col2 = mod(M[:, c], p).astype(np.float64); col2[r] = 0
        idx = np.nonzero(col2)[0]
        if len(idx):
            M[idx] = mod(M[idx] - np.outer(col2[idx], M[r]), p).astype(np.float64)
        piv.append(c); r += 1
    return mod(M[:r], p), piv

def nullspace_mod(M, p):
    m, n = M.shape
    if m == 0: return np.eye(n, dtype=np.int64)
    R, piv = rref_mod(M, p)
    return _nullspace_from_reduced(R, piv, n, p)

def _nullspace_from_reduced(R, piv, n, p):
    """R: (r,n) with R[i,piv[i]]=1 and R[i,piv[j]]=0 for j!=i (pivots need not
    be sorted).  Returns NS (n, n-r)."""
    pivset = set(piv)
    free = [c for c in range(n) if c not in pivset]
    NS = np.zeros((n, len(free)), dtype=np.int64)
    if free:
        Rf = R[:, free]  # (r, nfree)
        for j, f in enumerate(free):
            NS[f, j] = 1
        for i, pc in enumerate(piv):
            NS[pc, :] = (-Rf[i, :]) % p
    return NS

# ------------------------------------------------------ incremental nullspace
def incremental_nullspace(C_iter, n, p, batch=512):
    """Nullspace over GF(p) of the row-stack of C_iter (each row length n),
    processed in batches so memory is O((batch+rank)·n).  C_iter may be a
    2-D ndarray, a list of 1-D arrays, or any iterable of 1-D arrays."""
    if isinstance(C_iter, np.ndarray) and C_iter.ndim == 2:
        src = C_iter
        def chunks():
            for i in range(0, src.shape[0], batch):
                yield src[i:i+batch]
        gen = chunks()
    else:
        rows = list(C_iter)
        if not rows:
            return np.eye(n, dtype=np.int64)
        src = np.asarray(rows, dtype=np.int64)
        def chunks():
            for i in range(0, src.shape[0], batch):
                yield src[i:i+batch]
        gen = chunks()
    R = np.zeros((0, n), dtype=np.int64)
    piv = []
    piv_arr = np.zeros(0, dtype=np.int64)
    for Mb in gen:
        Mb = np.asarray(Mb, dtype=np.int64) % p
        if R.shape[0]:
            # forward-reduce batch against existing pivots (BLAS)
            coef = Mb[:, piv_arr]                        # (b, r)
            Mb = (Mb - matmul_mod(coef, R, p)) % p
        # drop zero rows early
        nz = np.any(Mb, axis=1)
        if not np.any(nz): continue
        Mb = Mb[nz]
        Rb, pivb = rref_mod(Mb, p)
        if Rb.shape[0] == 0: continue
        pivb_arr = np.asarray(pivb, dtype=np.int64)
        if R.shape[0]:
            # back-reduce existing rows against the new pivots
            coef2 = R[:, pivb_arr]                       # (r, rb)
            R = (R - matmul_mod(coef2, Rb, p)) % p
        R = np.vstack([R, Rb])
        piv.extend(pivb)
        piv_arr = np.asarray(piv, dtype=np.int64)
        if R.shape[0] >= n:       # full rank ⇒ nullspace = {0}
            return np.zeros((n, 0), dtype=np.int64)
    return _nullspace_from_reduced(R, piv, n, p)

# ------------------------------------------------------------- SparseWord
class SparseWord:
    """N basis vectors expanded on a sparse set of word-columns.

    mat  : (N, ncol) int64, entries in [0,p)
    cols : (ncol,)   int64, global word indices (big-endian base-A)
    pos  : dict word_idx -> local column position
    """
    __slots__ = ('mat', 'cols', 'pos', 'p')

    def __init__(self, mat, cols, p, prune=True):
        mat = np.asarray(mat, dtype=np.int64) % p
        cols = np.asarray(cols, dtype=np.int64)
        if prune and mat.size:
            nz = np.any(mat, axis=0)
            mat = np.ascontiguousarray(mat[:, nz])
            cols = cols[nz]
        self.mat = mat; self.cols = cols; self.p = p
        self.pos = {int(c): i for i, c in enumerate(cols)}

    @property
    def N(self):    return self.mat.shape[0]
    @property
    def ncol(self): return self.mat.shape[1]

    @classmethod
    def from_letters(cls, letters, p):
        """Weight-1 basis: one basis vector per allowed first-entry letter."""
        letters = sorted(letters)
        N = len(letters)
        return cls(np.eye(N, dtype=np.int64), np.asarray(letters, dtype=np.int64), p)

    def extend(self, Bw, A):
        """Tensor with the alphabet and rotate into the new recursive basis.

        Bw : (Nw, N_self, A)  — nullspace of the weight-w integrability(+extra)
             constraints in the (N_{w-1} × A) extension basis.
        Returns SparseWord with N = Nw and cols ⊆ {c·A+l : c∈self.cols, l<A}.

        This is exactly  Bw.reshape(Nw,-1) @ kron(self.mat, I_A)  but never
        materialises the kron: it loops over the A last-letter slots, each a
        (Nw × N_{w-1}) @ (N_{w-1} × ncol) BLAS matmul.
        """
        Nw, Nm1, A_ = Bw.shape
        assert Nm1 == self.N and A_ == A
        p = self.p
        # new_mat[:, j*A + l] = Bw[:,:,l] @ self.mat[:, j]
        chunks = [matmul_mod(Bw[:, :, l], self.mat, p) for l in range(A)]  # each (Nw, ncol)
        new_mat = np.stack(chunks, axis=2).reshape(Nw, self.ncol * A)
        new_cols = (self.cols[:, None] * A + np.arange(A, dtype=np.int64)).ravel()
        return SparseWord(new_mat, new_cols, p)

    def left_mul(self, T):
        """T : (M, N).  Returns SparseWord with mat = T @ self.mat (mod p)."""
        return SparseWord(matmul_mod(T, self.mat, self.p), self.cols.copy(), self.p)

    def constraint_matrix(self, rows):
        """rows : iterable of dict{word_idx: coef}.  Returns dense (nr, N)
        with entry r = Σ coef · self.mat[:, word].  Words outside the support
        contribute 0 (column absent ⇒ already zero)."""
        rows = list(rows)
        nr = len(rows)
        C = np.zeros((nr, self.N), dtype=np.int64)
        p = self.p
        for r, row in enumerate(rows):
            acc = None
            for c, v in row.items():
                j = self.pos.get(int(c))
                if j is None: continue
                term = (int(v) % p) * self.mat[:, j]
                acc = term if acc is None else acc + term
            if acc is not None:
                C[r, :] = acc % p
        return C

    def support_words(self):
        return [int(c) for c in self.cols]

# ---------------------------------------------------------------- projection
def project_sparse(W, rows, batch=512):
    """Impose linear constraints `rows` (iterable of dict{word:coef}) on the
    basis carried by SparseWord W.  Returns (new_dim, new_W)."""
    rows = list(rows)
    if not rows:
        return W.N, W
    C = W.constraint_matrix(rows)             # (nr, N)
    # drop trivially-zero constraint rows before reduction
    nz = np.any(C, axis=1)
    C = C[nz]
    if C.shape[0] == 0:
        return W.N, W
    NS = incremental_nullspace(C, W.N, W.p, batch=batch)   # (N, nfree)
    nfree = NS.shape[1]
    if nfree == 0:
        empty = SparseWord(np.zeros((0, 0), dtype=np.int64),
                           np.zeros(0, dtype=np.int64), W.p, prune=False)
        return 0, empty
    return nfree, W.left_mul(NS.T)

# ----------------------------------------------- support-restricted row-gens
def crossing_rows_on_support(cols, swap_word, p):
    """Yield s↔t crossing constraints {w:1, σ(w):-1} only for words touching
    the support (others are 0=0).  σ must be an involution on word indices."""
    sup = set(int(c) for c in cols)
    seen = set()
    for c in sup:
        j = swap_word(c)
        if c == j: continue
        key = (c, j) if c < j else (j, c)
        if key in seen: continue
        seen.add(key)
        yield {c: 1, j: p - 1}

def kill_rows_on_support(cols, predicate):
    """Yield {c:1} for each support word c with predicate(c)==True."""
    for c in cols:
        c = int(c)
        if predicate(c):
            yield {c: 1}

# ---------------------------------------------------------------- utilities
def word_digits(idx, w, A):
    k = idx; d = []
    for _ in range(w): d.append(k % A); k //= A
    return d[::-1]

def word_index(digits, A):
    o = 0
    for x in digits: o = o * A + x
    return o

# ---------------------------------------------------------------- probe CLI
if __name__ == '__main__':
    import sys
    _ok, _checks = matmul_boundary_probe()
    for _name, _o in _checks:
        print(f"[{'PASS' if _o else 'FAIL'}] {_name}")
    print(f"matmul_mod boundary probe: {'PASS' if _ok else 'FAIL'} "
          f"(envelope max p = {matmul_max_p()})")
    sys.exit(0 if _ok else 1)
