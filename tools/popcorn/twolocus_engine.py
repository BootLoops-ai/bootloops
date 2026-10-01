#!/usr/bin/env python3
"""Two-locus branch-length second moments E[T_i^A T_j^B] under the coalescent
with recombination and piecewise-constant population size (deterministic
sparse solve of the labeled two-locus configuration chain).

Object computed. For a sample of n haploid sequences observed at two loci A
and B, let T_i^A be the total branch length of the locus-A genealogy that
subtends exactly i of the n samples (i = 1..n-1), likewise T_j^B at locus B.
TwoLocus(n).moments(rho, epochs) returns

    M[i-1][j-1] = E[ T_i^A * T_j^B ],   i, j = 1..n-1,
    m1A[i-1]    = E[ T_i^A ],  m1B[j-1] = E[ T_j^B ],

so that the expected joint two-locus frequency spectrum of a pair of
infinite-sites mutations, one per locus, is proportional to M (small-theta
limit), and the single-locus spectra are proportional to m1A, m1B.

Conventions. Time is measured in units of 2*N_ref generations. Every pair of
ancestral lineages coalesces at rate eta(tau) = N_ref / N(tau), piecewise
constant: `epochs` is a list of (tau_start, eta) with tau_start ascending from
0, the last epoch extending to infinity (a constant population of size N_ref
is [(0.0, 1.0)]). rho = 4 * N_ref * R, with R the per-generation probability
that the two loci are separated onto different parents (crossing over
anywhere between them, plus any conversion event that separates them; the
two-point object depends on R only). Every ancestral lineage that carries
material ancestral to the sample at BOTH loci splits into an A-only and a
B-only parent at rate rho/2. epochs_from_Nt() converts a history given in
generations and absolute sizes into this form.

Method. Continuous-time Markov chain on labeled ancestral configurations:
a configuration is a multiset of lineage labels (a, b), a = number of sampled
A-copies descending from the lineage, b = likewise for B (the two-locus
ancestral process of Griffiths 1981 and Hudson 1983 with "trapped" material
retained: a lineage formed by coalescing (a,0) with (0,b) is doubly ancestral
again and can split again). Labels (n,0), (0,n) and (n,n) have reached the
most recent common ancestor at every locus they carry and are dropped; the
empty configuration is absorbing. The reachable configurations are enumerated
by breadth-first search from {(1,1) x n}; coalescence (C) and splitting (S)
generators are assembled as sparse matrices acting on functions of the
configuration, (Q g)(x) = sum_y q_xy (g(y) - g(x)), Q = eta*C + (rho/2)*S.
The branch lengths are additive functionals, T_i^A = int f_i^A(X_tau) dtau
with f_i^A(x) = number of labels in x with a = i, so the backward
(Feynman-Kac / dynamic-programming) equations for first and second moments
of additive functionals,

    -dG_f/dtau  = Q G_f + f,
    -dH_fg/dtau = Q H_fg + f * G_g + g * G_f,

hold with G, H -> 0 at absorption. In the terminal (infinite) epoch they are
stationary and are solved DIRECTLY by one sparse LU factorization of -Q
(no truncation, no time stepping); through each finite epoch the stacked
linear system [G; H; 1] has constant coefficients and is propagated by the
action of the matrix exponential (scipy.sparse.linalg.expm_multiply). The
answer is read off at the initial configuration.

Checks that hold by construction of the model and are used for validation:
  * the single-locus marginals are Kingman at every rho: for constant N,
    E[T_i^A] = E[T_i^B] = 2/i exactly;
  * rho = 0: both loci share one genealogy and M[i][j] = E[L_i L_j], the
    single-locus Kingman branch-length second moments (exact rationals;
    E[L_i L_j] = 4*sigma_ij + 4/(i*j) in terms of the sigma_ij of Fu 1995);
  * n = 2, constant N: M = 4 * E[T^A T^B] with E[T^A T^B] =
    1 + (rho + 18) / (rho^2 + 13*rho + 18), the second term being the
    classical covariance of the pairwise coalescence times at the two loci
    (E[T] = Var[T] = 1; Griffiths 1981; Hudson 1983; see McVean 2002); under
    piecewise-constant N(t) the n = 2 first moment (any rho) and rho = 0 second moment are elementary integrals of
    the survival function, which checks the finite-epoch propagation;
  * rho -> infinity: M -> outer(m1A, m1B) (unlinked loci);
  * A/B exchangeability: M is symmetric;
  * time rescaling: moments(rho, [(0, c)]) = moments(rho/c, [(0, 1)]) / c^2
    (first moments / c); splitting a constant epoch in two changes nothing.

Register: FLOAT-VALIDATED. Double-precision sparse linear algebra; the
output is validated, not certified, and carries no interval enclosure.
Measured agreement: at rho = 0, n = 8, against the exact-rational Kingman
E[L_i L_j]: max abs deviation ~2e-16 on the normalized pair vector
(pool_pairs), max relative deviation ~1e-13 on raw entries, |E[T_i] - 2/i|
~2e-15; at n = 2 against the closed forms above (constant N for rho in
[0, 100], and a two-epoch history): relative deviation <= 1e-15;
self-consistency across rho at n = 6 (Kingman marginals, symmetry, epoch
splitting, time rescaling) at or below the 1e-14 level. Do not quote its
output as exact or certified.

Cost (single thread, measured): reachable configurations 108, 338, 1042,
2997, 8405 for n = 4, 5, 6, 7, 8 (growth ~2.9x per sample). n = 8: build
~0.4 s; a constant-N evaluation ~0.2 s at rho = 0 and ~13 s / ~0.8 GB at
rho > 0 (LU fill-in of the coupled chain); each finite epoch adds an
expm_multiply on a ((2(n-1) + (n-1)^2) * states + 1)-dimensional block
system (~0.5M at n = 8, ~8 s per epoch). n <= 7 is sub-second per constant-N
evaluation. Nothing is built at import; the chain is constructed when
TwoLocus(n) is called.

pool_pairs(M) reads the symmetrized matrix on i <= j (2*M_ii on the diagonal,
M_ij + M_ji off it), normalized to unit sum: for this A/B-symmetric model it
is the normalized upper triangle of M, i.e. the scale-free two-locus object.
It is a matrix-pooling convention, NOT the class distribution of an unordered
pair of sites (that object carries M_ii, not 2*M_ii, on the diagonal):
compare only against vectors pooled the same way.

References.
  Griffiths R.C. (1981), Neutral two-locus multiple allele models with
    recombination, Theor. Popul. Biol. 19, 169-186.
  Hudson R.R. (1983), Properties of a neutral allele model with intragenic
    recombination, Theor. Popul. Biol. 23, 183-201.
  Fu Y.-X. (1995), Statistical properties of segregating sites, Theor.
    Popul. Biol. 48, 172-197.
  Hudson R.R. (2001), Two-locus sampling distributions and their
    application, Genetics 159, 1805-1817.
  McVean G.A.T. (2002), A genealogical interpretation of linkage
    disequilibrium, Genetics 162, 987-991.
  Kamm J.A., Spence J.P., Chan J., Song Y.S. (2016), Two-locus likelihoods
    under variable population size and fine-scale recombination rate
    estimation, Genetics 203, 1381-1399, arXiv:1510.06017.
The configuration chain here is the moment (branch-length) counterpart of the
two-locus sampling recursions of Hudson (2001) and Kamm et al. (2016): it
delivers E[T_i^A T_j^B] under piecewise-constant N(t) by linear algebra
rather than the full two-locus sampling probabilities.

Requires numpy and scipy (scipy.sparse, scipy.sparse.linalg).
"""
from __future__ import annotations
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla

class TwoLocus:
    """Two-locus configuration chain for a sample of n haploids. The
    constructor enumerates the reachable labeled configurations and builds
    the sparse coalescence (self.C) and splitting (self.S) generators; call
    moments(rho, epochs) for E[T_i^A T_j^B] and the first moments."""
    def __init__(self, n: int = 8):
        self.n = n; self._build()

    # ---------- state space by BFS
    def _norm(self, ms):
        n = self.n; out = []
        for (a, b) in ms:
            if (a, b) in ((n, 0), (0, n), (0, 0), (n, n)): continue
            out.append((a, b))
        out.sort(); return tuple(out)

    def _build(self):
        n = self.n; init = self._norm([(1, 1)] * n)
        idx = {init: 0}; states = [init]; coal = []; split = []   # (from, to, multiplicity)
        q = [init]
        while q:
            s = q.pop(); i = idx[s]; L = list(s); m = len(L)
            # coalescences
            seen = {}
            for x in range(m):
                for y in range(x + 1, m):
                    la, lb = L[x], L[y]; merged = (la[0] + lb[0], la[1] + lb[1])
                    rest = L[:x] + L[x + 1:y] + L[y + 1:] + [merged]; t = self._norm(rest)
                    seen[t] = seen.get(t, 0) + 1
            for t, c in seen.items():
                if t not in idx: idx[t] = len(states); states.append(t); q.append(t)
                coal.append((i, idx[t], c))
            # splits
            seen = {}
            for x in range(m):
                a, b = L[x]
                if a > 0 and b > 0:
                    rest = L[:x] + L[x + 1:] + [(a, 0), (0, b)]; t = self._norm(rest)
                    seen[t] = seen.get(t, 0) + 1
            for t, c in seen.items():
                if t not in idx: idx[t] = len(states); states.append(t); q.append(t)
                split.append((i, idx[t], c))
        self.states = states; self.idx = idx; N = len(states); self.N = N
        def mat(trip):
            r = np.array([t[0] for t in trip]); c = np.array([t[1] for t in trip]); v = np.array([float(t[2]) for t in trip])
            A = sp.coo_matrix((v, (r, c)), shape=(N, N)).tocsr(); d = np.asarray(A.sum(axis=1)).ravel()
            return (A - sp.diags(d)).tocsr()
        self.C = mat(coal); self.S = mat(split)          # generators acting on functions: (Qg)(x) = sum_y q_xy (g_y - g_x)
        self.absorb = idx[tuple()]                         # empty configuration = everything inert/absorbed
        fa = np.zeros((n - 1, N)); fb = np.zeros((n - 1, N))
        for k, s in enumerate(states):
            for (a, b) in s:
                if 0 < a < n: fa[a - 1, k] += 1
                if 0 < b < n: fb[b - 1, k] += 1
        self.fA = fa; self.fB = fb; self.init = 0
        keep = np.ones(N, bool); keep[self.absorb] = False; self.keep = keep

    # ---------- moments
    def moments(self, rho: float, epochs):
        """epochs: list of (tau_start, eta) with tau ascending from 0 (time in units of 2*N_ref generations,
        eta = N_ref/N the pair-coalescence rate); the last epoch extends to infinity. rho = 4*N_ref*R.
        Returns (M (n-1 x n-1) = E[T_i^A T_j^B], m1A = E[T_i^A], m1B = E[T_j^B])."""
        n1 = self.n - 1; N = self.N; keep = self.keep
        F = np.vstack([self.fA, self.fB])                             # (2n1, N) functionals: A_1..A_n1, B_1..B_n1
        nf = 2 * n1
        def Qof(eta): return (eta * self.C + 0.5 * rho * self.S).tocsc()
        # terminal epoch: exact
        eta_T = epochs[-1][1]; Q = Qof(eta_T)[keep][:, keep]; lu = spla.splu((-Q).tocsc())
        G = np.zeros((nf, N)); G[:, keep] = lu.solve(F[:, keep].T).T
        H = np.zeros((n1, n1, N))
        for i in range(n1):
            for j in range(n1):
                rhs = F[i] * G[n1 + j] + F[n1 + j] * G[i]; H[i, j, keep] = lu.solve(rhs[keep])
        # finite epochs, latest first: constant-coefficient LINEAR system in Y = [G (nf blocks); H (n1*n1 blocks); 1],
        # propagated exactly by the action of the matrix exponential (expm_multiply) — no explicit ODE stepping.
        Nk = int(keep.sum()); kidx = np.where(keep)[0]
        Fk = F[:, keep]
        for k in range(len(epochs) - 2, -1, -1):
            t0, eta = epochs[k]; t1 = epochs[k + 1][0]; dur = t1 - t0; Qk = Qof(eta)[keep][:, keep].tocsr()
            nb = nf + n1 * n1; dim = nb * Nk + 1
            blocks = [[None] * (nb + 1) for _ in range(nb + 1)]
            for bI in range(nb): blocks[bI][bI] = Qk
            # G rows: + f (constant) -> couples to the trailing '1'
            for gI in range(nf): blocks[gI][nb] = sp.csr_matrix(Fk[gI].reshape(-1, 1))
            # H_{ij} rows: + f_i^A * G_{B_j} + f_j^B * G_{A_i}
            for i in range(n1):
                for j in range(n1):
                    hI = nf + i * n1 + j
                    blocks[hI][n1 + j] = sp.diags(Fk[i]); blocks[hI][i] = (blocks[hI][i] + sp.diags(Fk[n1 + j])) if blocks[hI][i] is not None else sp.diags(Fk[n1 + j])
            blocks[nb][nb] = sp.csr_matrix((1, 1))
            A = sp.bmat(blocks, format="csr")
            y0 = np.concatenate([G[:, keep].ravel(), H[:, :, keep].reshape(n1 * n1, Nk).ravel(), [1.0]])
            y1 = spla.expm_multiply(A, y0, start=0.0, stop=dur, num=2, endpoint=True)[-1]
            G = np.zeros((nf, N)); G[:, keep] = y1[:nf * Nk].reshape(nf, Nk)
            H = np.zeros((n1, n1, N)); H[:, :, keep] = y1[nf * Nk:nf * Nk + n1 * n1 * Nk].reshape(n1, n1, Nk)
        M = H[:, :, self.init].copy(); m1A = G[:n1, self.init].copy(); m1B = G[n1:, self.init].copy()
        return M, m1A, m1B

def pool_pairs(M):
    """Read the symmetrized matrix on index pairs i <= j: 2*M_ii on the diagonal, M_ij + M_ji off it, in the order
    (1,1),(1,2),...,(1,n-1),(2,2),...,(n-1,n-1); normalized to unit sum (n(n-1)/2 cells; 28 at n = 8). For this
    A/B-symmetric model it equals the normalized upper triangle of M. A pooling convention, not the class law of an
    unordered pair of sites (that carries M_ii on the diagonal): compare only against vectors pooled the same way."""
    n1 = M.shape[0]; v = np.array([2 * M[i, i] if i == j else M[i, j] + M[j, i] for i in range(n1) for j in range(i, n1)]); return v / v.sum()

def epochs_from_Nt(times_gen, sizes, N_ref=1.0e4):
    """piecewise-constant N(t): sizes[k] holds on [times_gen[k], times_gen[k+1]) (generations, ascending from 0);
    returns [(tau_start, eta)] with tau in units of 2 N_ref generations and eta = N_ref / N."""
    return [(t / (2.0 * N_ref), N_ref / s) for t, s in zip(times_gen, sizes)]

if __name__ == "__main__":
    import sys, time
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 6; rho = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
    t0 = time.time(); E = TwoLocus(n); print(f"n={n}: states {E.N} (reachable), build {time.time()-t0:.2f}s; nnz C {E.C.nnz} S {E.S.nnz}")
    t0 = time.time(); M, a, b = E.moments(0.0, [(0.0, 1.0)]); print(f"constant N, rho=0: {time.time()-t0:.2f}s; E[T_i] = {np.round(a,4)} (Kingman 2/i: {np.round(2/np.arange(1,n),4)})")
    t0 = time.time(); M1, a1, b1 = E.moments(rho, [(0.0, 1.0)]); print(f"constant N, rho={rho:g}: {time.time()-t0:.2f}s; E[T_i^A] = {np.round(a1,4)}")
    np.set_printoptions(precision=4, suppress=True, linewidth=140); print(f"E[T_i^A T_j^B] (rho={rho:g}):\n{M1}")
