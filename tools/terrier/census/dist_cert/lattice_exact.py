"""dist_cert.lattice_exact -- exact lattice predicates for G2/G3 (eps0 = 0).

G2 metric = exact lattice predicate (CM + quotient-surviving roots); G3
metric = exact lattice predicate, the unfrozen-sector root locus, defined
per row by the NAMED frozen sublattice (rk 17/14/8). Special set S = exact
root/CM membership, closed-form. Membership is EXACT integer arithmetic:
dist_cert = 0 iff the predicate fires; a certified OFF verdict = exhaustive
Fincke-Pohst root enumeration returning empty (a finite exact certificate).
CM stratum: explicit finite certified list ONLY (else dropped) --
membership = exact equality against that list.
NO census data is touched here; sign convention: even NEGATIVE-definite
sectors, roots r have r.G.r = -2.
"""
from fractions import Fraction

F = Fraction


def mat_vec(G, v):
    return [sum(G[i][j] * v[j] for j in range(len(v))) for i in range(len(G))]


def inner(G, u, v):
    return sum(u[i] * x for i, x in enumerate(mat_vec(G, v)))


def gram_of_basis(G, B):
    """Gram matrix of the sublattice spanned by rows of B (parent Gram G)."""
    return [[inner(G, bi, bj) for bj in B] for bi in B]


def gram_An(n):
    """A_n root lattice, NEGATIVE definite (diag -2, offdiag +1)."""
    return [[-2 if i == j else (1 if abs(i - j) == 1 else 0)
             for j in range(n)] for i in range(n)]


def gram_Dn(n):
    """D_n root lattice (n >= 3), NEGATIVE definite."""
    G = [[0] * n for _ in range(n)]
    for i in range(n):
        G[i][i] = -2
    for i in range(n - 2):
        G[i][i + 1] = G[i + 1][i] = 1
    G[n - 3][n - 1] = G[n - 1][n - 3] = 1
    return G


def direct_sum(A, B):
    na, nb = len(A), len(B)
    G = [[0] * (na + nb) for _ in range(na + nb)]
    for i in range(na):
        G[i][:na] = list(A[i])
    for i in range(nb):
        for j in range(nb):
            G[na + i][na + j] = B[i][j]
    return G


# G3 frozen sublattices (named realizations):
#   G3-1 [1 1 1]: rk 17, realized A17; G3-2 [2 0 2]: rk 14, realized A7+A7;
#   G3-3 [6 6 6]: rk 8, realized D4+D4.
def frozen_gram_G3(row):
    if row == 1:
        return gram_An(17)
    if row == 2:
        return direct_sum(gram_An(7), gram_An(7))
    if row == 3:
        return direct_sum(gram_Dn(4), gram_Dn(4))
    raise ValueError("G3 row must be 1, 2 or 3")


def integer_kernel(A):
    """Saturated Z-basis of {x in Z^n : A x = 0} for integer matrix A (m x n).

    Column-HNF style: unimodular column ops over an identity tracker; columns
    zeroed in A give a saturated kernel basis (exact integer arithmetic).
    """
    m, n = len(A), len(A[0])
    W = [list(r) for r in A]
    U = [[1 if i == j else 0 for j in range(n)] for i in range(n)]

    def colop(j, k, q):  # col_k -= q * col_j
        for i in range(m):
            W[i][k] -= q * W[i][j]
        for i in range(n):
            U[i][k] -= q * U[i][j]

    def colswap(j, k):
        for i in range(m):
            W[i][j], W[i][k] = W[i][k], W[i][j]
        for i in range(n):
            U[i][j], U[i][k] = U[i][k], U[i][j]

    r = 0
    for i in range(m):
        piv = None
        while True:
            nz = [j for j in range(r, n) if W[i][j] != 0]
            if not nz:
                piv = None
                break
            j0 = min(nz, key=lambda j: abs(W[i][j]))
            colswap(r, j0)
            done = True
            for k in range(r + 1, n):
                if W[i][k] != 0:
                    colop(r, k, W[i][k] // W[i][r])
                    if W[i][k] != 0:
                        done = False
            if done:
                piv = r
                break
        if piv is not None:
            r += 1
    return [[U[i][j] for i in range(n)] for j in range(r, n)]


def _ldl(P):
    """Exact LDL^T of a symmetric positive-definite rational matrix P."""
    n = len(P)
    L = [[F(1) if i == j else F(0) for j in range(n)] for i in range(n)]
    D = [F(0)] * n
    A = [[F(x) for x in row] for row in P]
    for j in range(n):
        D[j] = A[j][j] - sum(L[j][k] * L[j][k] * D[k] for k in range(j))
        if D[j] <= 0:
            raise ValueError("matrix not positive definite")
        for i in range(j + 1, n):
            L[i][j] = (A[i][j] - sum(L[i][k] * L[j][k] * D[k]
                                     for k in range(j))) / D[j]
    return L, D


def short_vectors(P, target=2, node_cap=2_000_000):
    """ALL x in Z^n \\ {0} with x^T P x == target (P positive definite, exact).

    Complete Fincke-Pohst enumeration; exact Fraction arithmetic throughout.
    Returns the full list (both signs). node_cap guards runaway inputs.
    """
    n = len(P)
    L, D = _ldl(P)
    from math import floor
    out, nodes = [], [0]
    x = [0] * n

    def rec(i, rem):
        nodes[0] += 1
        if nodes[0] > node_cap:
            raise RuntimeError("short_vectors: node cap exceeded")
        if i < 0:
            if rem == 0 and any(v != 0 for v in x):
                out.append(list(x))
            return
        c = sum(L[j][i] * x[j] for j in range(i + 1, n))
        # need integers m with D[i]*(m + c)^2 <= rem; scan out from floor(-c)
        m = floor(-c)
        while D[i] * (m + c) * (m + c) <= rem:
            x[i] = m
            rec(i - 1, rem - D[i] * (m + c) * (m + c))
            m -= 1
        m = floor(-c) + 1
        while D[i] * (m + c) * (m + c) <= rem:
            x[i] = m
            rec(i - 1, rem - D[i] * (m + c) * (m + c))
            m += 1
        x[i] = 0

    rec(n - 1, F(target))
    return out


def roots(neg_gram, node_cap=2_000_000):
    """All roots (norm -2 vectors) of an even NEGATIVE-definite Gram matrix.
    Complete exact enumeration -- empty return is a certified OFF certificate."""
    if not neg_gram:
        return []
    P = [[-x for x in row] for row in neg_gram]
    return short_vectors(P, 2, node_cap)


def rank_Q(vecs):
    """Exact rank over Q of a list of integer/rational vectors."""
    M = [[F(x) for x in v] for v in vecs if any(x != 0 for x in v)]
    r = 0
    for col in range(len(M[0]) if M else 0):
        piv = next((i for i in range(r, len(M)) if M[i][col] != 0), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        for i in range(len(M)):
            if i != r and M[i][col] != 0:
                f = M[i][col] / M[r][col]
                M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        r += 1
    return r


def on_locus_root(parent_gram, flux_vecs, sector_basis=None, frozen_basis=None,
                  node_cap=2_000_000):
    """G2/G3 root predicate. ON-locus iff the
    flux-orthogonal part of the censused sector contains a root NOT in the
    Q-span of the frozen sublattice (G3: named W_root row; G2: quotient-
    surviving roots via sector_basis = invariant sublattice basis).

    Returns (on: bool, witness_parent_coords or None, n_roots_examined).
    Exact integer/rational arithmetic; complete enumeration certificate.
    """
    G = parent_gram
    n = len(G)
    B = sector_basis if sector_basis is not None else \
        [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    if flux_vecs:
        A = [[inner(G, b, f) for b in B] for f in flux_vecs]
        K = integer_kernel(A)
    else:
        K = [[1 if i == j else 0 for j in range(len(B))] for i in range(len(B))]
    if not K:
        return False, None, 0
    Bp = [[sum(k[a] * B[a][j] for a in range(len(B))) for j in range(n)]
          for k in K]
    Gp = gram_of_basis(G, Bp)
    rts = roots(Gp, node_cap)
    base = list(frozen_basis) if frozen_basis else []
    r0 = rank_Q(base) if base else 0
    for y in rts:
        r_par = [sum(y[a] * Bp[a][j] for a in range(len(Bp))) for j in range(n)]
        if not base or rank_Q(base + [r_par]) > r0:
            return True, r_par, len(rts)
    return False, None, len(rts)


def cm_member(candidate_key, frozen_cm_list):
    """CM stratum membership: EXACT equality against the explicit finite
    certified list (empty list = CM dropped)."""
    return candidate_key in frozen_cm_list


def dist_cert_exact(parent_gram, flux_vecs, sector_basis=None,
                    frozen_basis=None, cm_key=None, frozen_cm_list=()):
    """G2/G3 dist_cert: 0 (exact) iff on-locus, else certified OFF (eps0=0).
    Returns (0, witness) or ('OFF', {'roots_examined': k})."""
    on, wit, k = on_locus_root(parent_gram, flux_vecs, sector_basis,
                               frozen_basis)
    if on:
        return 0, wit
    if cm_key is not None and cm_member(cm_key, frozen_cm_list):
        return 0, ("CM", cm_key)
    return "OFF", {"roots_examined": k}
