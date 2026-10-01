#!/usr/bin/env python3
"""twosfs_engine — exact second moments E[L_i L_j] of the branch-length spectrum
under Lambda-coalescents, and the two-time Kingman block-count kernel.

Register: EXACT. Every moment, probability and polynomial coefficient is a
fractions.Fraction (the two-time kernel is a sympy polynomial with exact
rational coefficients); no float touches a returned value or a sign decision.
Standard library only, except twotime_kernel_poly, which imports sympy when
called. Rate families come from lambda_exact (shipped beside this file).

Object
------
For an exchangeable coalescent on a sample of n, let L_i (i = 1..n-1) be the
total branch length subtending exactly i leaves, and L_tot = sum_i L_i. Two
completely linked sites (no recombination between them) sit on ONE genealogy,
so under the infinite-sites model with mutation rate theta/2 per unit branch
length the expected numbers of segregating sites satisfy
    E[xi_i]          = (theta/2)   E[L_i]
    E[xi_i xi_j]     = (theta/2)^2 E[L_i L_j]            (i != j)
    E[xi_i (xi_i-1)] = (theta/2)^2 E[L_i^2],
i.e. the second factorial moments of the site-frequency spectrum of a
non-recombining locus are the matrix E[L_i L_j]: the expected two-site
frequency spectrum (2-SFS) of linked sites, up to scale. The scale-free object
    q_ij = E[L_i L_j] / E[L_tot^2]        (sum over ordered (i,j) equals 1)
is invariant under any global rescaling of time and is the distribution of
the ordered class pair of the two segregating sites conditional on exactly
two segregating sites, in the limit theta -> 0. First moments E[L_i] are the
expected (unfolded) SFS up to theta/2; they are returned alongside.

What is computed, and how (all exact rationals)
-----------------------------------------------
(1) lambda_moments(n, lam) -> (u, v), u[i] = E[L_i], v[(i,j)] = E[L_i L_j]
    (i <= j), for the Lambda-coalescent with merger rates lam(b, k) (any
    callable returning Fractions; lambda_exact provides kingman_rate(),
    beta_rate(alpha), dirac_rate(psi), msprime_dirac_rate(psi, c)).
    Lumped chain on block-size multisets (integer partitions of n): from a
    state with b blocks each specific set of k blocks merges at rate lam(b,k);
    a size-multiset choice {j_c copies of size c} has prod_c C(mult_c, j_c)
    realizations. Reward r_i(s) = number of blocks of size i accrues L_i at
    rate r_i during the holding time. First-step second moments: in state s
    with total rate q, T ~ Exp(q) independent of the next state S', and
    Y_i := (L_i accrued from s on) = r_i T + Y_i'(S'), hence
        u_i(s)  = r_i/q + sum_s' P(s,s') u_i(s')
        v_ij(s) = 2 r_i r_j / q^2
                + (1/q) sum_s' P(s,s') [ r_i u_j(s') + r_j u_i(s') ]
                + sum_s' P(s,s') v_ij(s')
    using E[T] = 1/q, E[T^2] = 2/q^2. The chain only coarsens, so this is a
    single pass over states ordered by block count (no linear solve). The
    total rate out of each state is asserted equal to sum_k C(b,k) lam(b,k).
    Hand check, Kingman n = 3 (L_1 = 3 T_3 + T_2, L_2 = T_2, T_k ~ Exp(C(k,2))):
    u = (2, 1), v = [[6, 3], [3, 2]].
(2) setpartition_moments(n, lam): the same recursion on raw set partitions of
    {0..n-1}, an independent state space; must agree entrywise with (1).
(3) kingman_moments_route2(n): Kingman only. Holding times T_k ~ Exp(C(k,2))
    are independent of the jump chain (uniform pair merges), so
        E[L_i L_j] = sum_{k,l} E[T_k T_l] E[N_i^(k) N_j^(l)],
    N_i^(k) = number of size-i blocks while k blocks remain; same-level and
    cross-level joint moments of the jump chain by exact enumeration
    (jump_levels, cond_level_dist). An independent decomposition that must
    agree entrywise with (1); the Kingman values are the branch-length
    covariances of Fu (1995).
(4) Kingman block-counting process and single-time moments. The number of
    ancestral lineages K(g) at coalescent intensity g is the pure-death chain
    b -> b-1 at rate C(b,2) (Tavare 1984). With the level t = e^{-g},
        P_{b,k}(t) := P(K(g) = k | K(0) = b) = sum_j c_{bkj} t^{C(j,2)},
    obtained here by exact spectral decomposition of the bidiagonal generator
    (eigenvalues -C(j,2), j = 1..n, all distinct): death_semigroup_coeffs(n),
    death_prob(n, b, k, t). Since the partition at intensity g is the jump
    chain stopped at K(g) blocks, P(partition = a) = P_{n,|a|}(t) x
    jump_levels(n)[|a|][a], and component_moments(n, t) returns
        m1[i] = E[N_i(g)],  m2[(i,j)] = E[N_i(g) N_j(g)],   t = e^{-g},
    exactly for rational t in [0, 1] (t = 0 is the directional limit t -> 0+,
    where the last two lineages dominate; t = 1 is the star, all singletons).
(5) twotime_kernel_poly(n) (needs sympy): the symmetrized two-time kernel
        Csym_ij(s, u) = ( E[N_i(g_s) N_j(g_su)] + E[N_j(g_s) N_i(g_su)] ) / 2,
    s = e^{-g_early} in [0,1], u = e^{-(g_late - g_early)} in [0,1], as sympy
    polynomials in (s, u) with rational coefficients (i <= j), from (4) and
    the Markov property: sum_a P(Pi(g_s) = a) N_i(a) sum_l P_{|a|,l}(u)
    E[N_j at l blocks | start a]. The (s, u) square maps onto the triangle
    y <= x of the (x, y) = (s, s u) unit square; the symmetrization covers
    the other triangle.

Why (4)-(5) matter: the variable-population-size Kingman class as a kernel
--------------------------------------------------------------------------
A single-population Kingman coalescent with any deterministic size history
is a time change of the standard coalescent (Griffiths & Tavare 1994): with
cumulative intensity R(tau) and level s(tau) = e^{-R(tau)}, let nu be the
push-forward of Lebesgue measure d tau under s (a measure on (0, 1]). Since
L_i = int N_i(R(tau)) d tau, Fubini and the Markov property give
    E[L_i]     = int      m1_i(x)       nu(dx)
    E[L_i L_j] = int int  Csym_ij(x, y) nu(dx) nu(dy),
so the expected 1-SFS of the whole class is the cone spanned by the curve
t -> m1(t), and the expected 2-SFS of the class, together with every mixture
across loci of such histories ("pooling": sums of nu (x) nu terms), lies in
the cone of the kernel. Constant size N corresponds to nu(dx) = dx/x, and the
kernel then reproduces lambda_moments(n, kingman_rate()) exactly (every
monomial s^a u^b of Csym has a, b >= 1 and contributes 2 c_ab/(a b); this
identity is one of the shipped self-checks). A single atom nu = D delta_t is
the limit history that coalesces instantly to intensity g = -log t, holds
with no coalescence for a duration D, then coalesces instantly to the root;
its moments are D m1(t) and D^2 m2(t). These representations are what
twosfs_certificates uses to certify that a given exact 2-SFS (for example a
Beta- or Dirac-coalescent's) lies outside the 2-SFS set of every pooled
variable-size Kingman history, or inside the set of atom poolings.

Helpers: int_partitions, multiset_choices, lumped_transitions, jump_levels,
cond_level_dist, normalize_matrix(v, n) -> (q, E[L_tot^2]) with the full
(i,j)+(j,i) sum in the denominator, flatten(q, n) -> list over i <= j.

Cost (pure-Python Fractions, single core, indicative): lambda_moments is
milliseconds at n <= 8 and grows with the number of integer partitions of n
and with the operand size of the rate family (Beta(3/2): n = 10, 42 states,
0.06 s; n = 15, 176 states, 1.3 s; n = 20, 627 states, ~16 s; Kingman
n = 20 ~2.5 s); setpartition_moments enumerates all Bell(n) set partitions
(n = 6: 203, 0.1 s; a small-n cross-check only); kingman_moments_route2 is
a cross-check too (n = 15: ~13 s); component_moments n = 15: 0.6 s;
twotime_kernel_poly ~0.03-0.2 s for n = 3..6 and ~1 s at n = 8. Memory
tens of MB. Scope: neutral, single population, no
recombination between the two sites, constant size on the Lambda side
(lambda_moments), arbitrary deterministic size history on the Kingman side
through the kernel; expected values (first and second moments) only.

CLI:  python3 twosfs_engine.py [--n N] [--family kingman|beta|dirac]
                               [--param P] [--out FILE]
  runs the built-in exact identity checks at sample size n (exit 1 on any
  failure), prints E[L_i], E[L_i L_j] and q for the requested coalescent
  (parameter as a rational string, e.g. 3/2 or 1/4), and writes a JSON with
  exact num/den strings only if --out is given.

References
  Kingman, J. F. C. (1982). The coalescent. Stochastic Process. Appl. 13,
    235-248.
  Tavare, S. (1984). Line-of-descent and genealogical processes, and their
    applications in population genetics models. Theor. Popul. Biol. 26,
    119-164.
  Griffiths, R. C. & Tavare, S. (1994). Sampling theory for neutral alleles
    in a varying environment. Phil. Trans. R. Soc. Lond. B 344, 403-410.
  Fu, Y.-X. (1995). Statistical properties of segregating sites. Theor.
    Popul. Biol. 48, 172-197.
  Pitman, J. (1999). Coalescents with multiple collisions. Ann. Probab. 27,
    1870-1902.
  Sagitov, S. (1999). The general coalescent with asynchronous mergers of
    ancestral lines. J. Appl. Probab. 36, 1116-1125.
  Birkner, M., Blath, J. & Eldon, B. (2013). Statistical properties of the
    site-frequency spectrum associated with Lambda-coalescents. Genetics 195,
    1037-1053.  (Second moments of the SFS under Lambda-coalescents.)
  Ferretti, L., Klassmann, A., Raineri, E., Ramos-Onsins, S. E., Wiehe, T. &
    Achaz, G. (2018). The neutral frequency spectrum of linked sites. Theor.
    Popul. Biol. 123, 70-79.  (Joint spectrum of two linked sites, Kingman.)
  Hobolth, A., Siri-Jegousse, A. & Bladt, M. (2019). Phase-type distributions
    in population genetics. Theor. Popul. Biol. 127, 16-32.  (SFS moments
    under Lambda-coalescents by phase-type methods.)
  Rice, D. P., Novembre, J. & Desai, M. M. (2018). Distinguishing
    multiple-merger from Kingman coalescence using two-site frequency
    spectra. bioRxiv 461517; journal version: Fenton, E. F., Rice, D. P.,
    Novembre, J. & Desai, M. M. (2025). Detecting deviations from Kingman
    coalescence using 2-site frequency spectra. Genetics 229(4), iyaf023.
"""
import argparse
import json
import sys
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from math import comb

from lambda_exact import (beta_rate, dirac_rate, kingman_rate,          # noqa
                          msprime_dirac_rate, expected_lengths)

# ---------------------------------------------------------------- states

@lru_cache(maxsize=None)
def int_partitions(n):
    """All integer partitions of n as sorted-descending tuples."""
    out = []
    def rec(rem, mx, acc):
        if rem == 0:
            out.append(tuple(acc))
            return
        for p in range(min(rem, mx), 0, -1):
            rec(rem - p, p, acc + [p])
    rec(n, n, [])
    return out

def multiset_choices(state, k):
    """Distinct k-multisets of block sizes from `state`, with counts.

    Yields (chosen_sizes_tuple, ways) where ways = prod_c C(mult_c, j_c)."""
    from collections import Counter
    cnt = Counter(state)
    sizes = sorted(cnt)
    out = []
    def rec(idx, left, acc, ways):
        if left == 0:
            out.append((tuple(acc), ways))
            return
        if idx == len(sizes):
            return
        c = sizes[idx]
        for j in range(0, min(left, cnt[c]) + 1):
            rec(idx + 1, left - j, acc + [c] * j, ways * comb(cnt[c], j))
    rec(0, k, [], 1)
    return out

def lumped_transitions(state, lam):
    """[(rate, new_state)] for a k-merger Lambda chain on block-size multisets."""
    b = len(state)
    out = []
    for k in range(2, b + 1):
        r = lam(b, k)
        if r == 0:
            continue
        for chosen, ways in multiset_choices(state, k):
            new = list(state)
            for c in chosen:
                new.remove(c)
            new.append(sum(chosen))
            out.append((ways * r, tuple(sorted(new, reverse=True))))
    return out

# ------------------------------------------------- (1) lumped moment DP

def lambda_moments(n, lam):
    """Exact (u, v): u[i] = E[L_i], v[(i,j)] = E[L_i L_j], i<=j, i,j=1..n-1.

    First-step DP over block-size multisets (see module docstring)."""
    states = sorted(int_partitions(n), key=len)          # ascending #blocks
    U = {}
    V = {}
    for s in states:
        b = len(s)
        if b == 1:
            U[s] = {i: F(0) for i in range(1, n)}
            V[s] = {(i, j): F(0) for i in range(1, n) for j in range(i, n)}
            continue
        trans = lumped_transitions(s, lam)
        q = sum(r for r, _ in trans)
        # check: the total rate must equal sum_k C(b,k) lam(b,k)
        qref = sum(comb(b, k) * lam(b, k) for k in range(2, b + 1))
        assert q == qref, (s, q, qref)
        assert q > 0, s
        r_ = {i: F(s.count(i)) for i in range(1, n)}
        u = {i: r_[i] / q for i in range(1, n)}
        v = {(i, j): 2 * r_[i] * r_[j] / q ** 2
             for i in range(1, n) for j in range(i, n)}
        for rate, s2 in trans:
            P = rate / q
            for i in range(1, n):
                u[i] += P * U[s2][i]
            for i in range(1, n):
                for j in range(i, n):
                    v[(i, j)] += P * ((r_[i] * U[s2][j] + r_[j] * U[s2][i]) / q
                                      + V[s2][(i, j)])
        U[s], V[s] = u, v
    init = tuple([1] * n)
    return U[init], V[init]

# --------------------------------------- (2) set-partition cross-check DP

def setpartition_moments(n, lam):
    """Same DP on raw set partitions of [n] (independent state space)."""
    def parts_of(elems):
        if not elems:
            yield []
            return
        first, rest = elems[0], elems[1:]
        for p in parts_of(rest):
            for i in range(len(p)):
                yield p[:i] + [p[i] | {first}] + p[i + 1:]
            yield p + [[first] if isinstance(first, set) else {first}]
    allp = [tuple(sorted((frozenset(b) for b in p), key=lambda x: sorted(x)))
            for p in parts_of(list(range(n)))]
    allp = sorted(set(allp), key=len)
    U, V = {}, {}
    for s in allp:
        b = len(s)
        if b == 1:
            U[s] = {i: F(0) for i in range(1, n)}
            V[s] = {(i, j): F(0) for i in range(1, n) for j in range(i, n)}
            continue
        trans = []
        for k in range(2, b + 1):
            r = lam(b, k)
            if r == 0:
                continue
            for idxs in combinations(range(b), k):
                merged = frozenset().union(*(s[x] for x in idxs))
                new = [blk for x, blk in enumerate(s) if x not in idxs]
                new.append(merged)
                trans.append((r, tuple(sorted(new, key=lambda x: sorted(x)))))
        q = sum(r for r, _ in trans)
        r_ = {i: F(sum(1 for blk in s if len(blk) == i)) for i in range(1, n)}
        u = {i: r_[i] / q for i in range(1, n)}
        v = {(i, j): 2 * r_[i] * r_[j] / q ** 2
             for i in range(1, n) for j in range(i, n)}
        for rate, s2 in trans:
            P = rate / q
            for i in range(1, n):
                u[i] += P * U[s2][i]
            for i in range(1, n):
                for j in range(i, n):
                    v[(i, j)] += P * ((r_[i] * U[s2][j] + r_[j] * U[s2][i]) / q
                                      + V[s2][(i, j)])
        U[s], V[s] = u, v
    init = tuple(sorted((frozenset([x]) for x in range(n)),
                        key=lambda x: sorted(x)))
    return U[init], V[init]

# ----------------------------- Kingman jump chain (uniform pair merges)

@lru_cache(maxsize=None)
def jump_levels(n):
    """levels[k] = {state: prob} of the Kingman jump chain at k blocks."""
    levels = {n: {tuple([1] * n): F(1)}}
    for b in range(n, 1, -1):
        nxt = {}
        for s, p in levels[b].items():
            tot = comb(b, 2)
            for chosen, ways in multiset_choices(s, 2):
                new = list(s)
                for c in chosen:
                    new.remove(c)
                new.append(sum(chosen))
                s2 = tuple(sorted(new, reverse=True))
                nxt[s2] = nxt.get(s2, F(0)) + p * F(ways, tot)
        levels[b - 1] = nxt
    return levels

def cond_level_dist(start, l):
    """{state: prob} of the jump chain at l blocks given current state `start`."""
    cur = {start: F(1)}
    for b in range(len(start), l, -1):
        nxt = {}
        for s, p in cur.items():
            tot = comb(b, 2)
            for chosen, ways in multiset_choices(s, 2):
                new = list(s)
                for c in chosen:
                    new.remove(c)
                new.append(sum(chosen))
                s2 = tuple(sorted(new, reverse=True))
                nxt[s2] = nxt.get(s2, F(0)) + p * F(ways, tot)
        cur = nxt
    return cur

def _N(s, i):
    return s.count(i)

# ------------------------------- (3) Kingman route-2 second moments

def kingman_moments_route2(n):
    """E[L_i], E[L_i L_j] via T_k x jump-chain decomposition (Kingman only)."""
    levels = jump_levels(n)
    ET = {k: F(1, comb(k, 2)) for k in range(2, n + 1)}
    ET2 = {k: F(2, comb(k, 2) ** 2) for k in range(2, n + 1)}
    u = {i: sum(p * _N(s, i) * ET[k]
                for k in range(2, n + 1) for s, p in levels[k].items())
         for i in range(1, n)}
    v = {}
    # cross-level joint moments: for k > l,
    #   E[N_i^{(k)} N_j^{(l)}] = sum_{s at k} P(s) N_i(s) E[N_j | jump to l from s]
    for i in range(1, n):
        for j in range(i, n):
            acc = F(0)
            for k in range(2, n + 1):
                # same level
                m = sum(p * _N(s, i) * _N(s, j) for s, p in levels[k].items())
                acc += ET2[k] * m
                for l in range(2, k):
                    cij = F(0)
                    cji = F(0)
                    for s, p in levels[k].items():
                        if _N(s, i) == 0 and _N(s, j) == 0:
                            continue
                        dl = cond_level_dist(s, l)
                        ei = sum(pp * _N(s2, i) for s2, pp in dl.items())
                        ej = sum(pp * _N(s2, j) for s2, pp in dl.items())
                        cij += p * _N(s, i) * ej
                        cji += p * _N(s, j) * ei
                    acc += ET[k] * ET[l] * (cij + cji)
            v[(i, j)] = acc
    return u, v

# --------------------- (4) death-chain semigroup + single-time moments

@lru_cache(maxsize=None)
def death_semigroup_coeffs(n):
    """coeffs[b][k][j]: P_{b,k}(t) = sum_j coeffs[b][k][j] * t^C(j,2).

    Exact spectral decomposition of the pure-death chain b -> b-1 at rate
    C(b,2) on {n,...,1}: P(g) = exp(gQ); with t = e^{-g},
    exp(gQ) = sum_j t^{C(j,2)} R_j, R_j = prod_{m != j} (Q - lam_m I)/(lam_j - lam_m),
    lam_j = -C(j,2) all distinct for j = 1..n."""
    idx = {b: n - b for b in range(n, 0, -1)}     # matrix index
    dim = n
    Q = [[F(0)] * dim for _ in range(dim)]
    for b in range(n, 1, -1):
        Q[idx[b]][idx[b]] = -F(comb(b, 2))
        Q[idx[b]][idx[b - 1]] = F(comb(b, 2))
    lams = {j: -F(comb(j, 2)) for j in range(1, n + 1)}
    def matmul(A, B):
        return [[sum(A[r][m] * B[m][c] for m in range(dim))
                 for c in range(dim)] for r in range(dim)]
    R = {}
    for j in range(1, n + 1):
        M = [[F(1) if r == c else F(0) for c in range(dim)] for r in range(dim)]
        for m in range(1, n + 1):
            if m == j:
                continue
            A = [[(Q[r][c] - (lams[m] if r == c else 0)) / (lams[j] - lams[m])
                  for c in range(dim)] for r in range(dim)]
            M = matmul(M, A)
        R[j] = M
    coeffs = {}
    for b in range(1, n + 1):
        coeffs[b] = {}
        for k in range(1, b + 1):
            coeffs[b][k] = {j: R[j][idx[b]][idx[k]] for j in range(1, n + 1)
                            if R[j][idx[b]][idx[k]] != 0}
    return coeffs

def death_prob(n, b, k, t):
    """P(block count k at intensity g | started at b), t = e^{-g} rational."""
    t = F(t)
    co = death_semigroup_coeffs(n)[b].get(k, {})
    return sum(c * t ** comb(j, 2) for j, c in co.items())

def component_moments(n, t):
    """Single-time Kingman block-count moments at level t: (m1, m2) with
    m1[i] = E[N_i(g)], m2[(i,j)] = E[N_i N_j](g), t = e^{-g} rational,
    N_i(g) = number of blocks of size i of the Kingman n-coalescent at
    intensity g. These are the unnormalized first and second moments
    (D = 1) of the single-atom limit history nu = D delta_t (see the module
    docstring); m1 normalized to sum 1 is the extreme point at level t of
    the normalized expected 1-SFS set of the variable-size Kingman class.

    t = 0 is the DIRECTIONAL limit t -> 0+: every P_{n,k}(t) has leading
    order t^{C(k,2)} (slowest surviving eigenvalue is -C(k,2)), so the
    normalized moments are dominated by the level-2 (last-two-lineages)
    conditioned moments; we return those (the exact limit of the normalized
    moments). t = 1 (g = 0) is the star: all n singletons."""
    t = F(t)
    if t == 0:
        levels = jump_levels(n)
        m1 = {i: F(0) for i in range(1, n)}
        m2 = {(i, j): F(0) for i in range(1, n) for j in range(i, n)}
        for s, p in levels[2].items():
            for i in range(1, n):
                Ni = _N(s, i)
                if Ni:
                    m1[i] += p * Ni
                    for j in range(i, n):
                        m2[(i, j)] += p * Ni * _N(s, j)
        return m1, m2
    levels = jump_levels(n)
    m1 = {i: F(0) for i in range(1, n)}
    m2 = {(i, j): F(0) for i in range(1, n) for j in range(i, n)}
    for k in range(2, n + 1):
        pk = death_prob(n, n, k, t)
        if pk == 0:
            continue
        for s, p in levels[k].items():
            w = pk * p
            for i in range(1, n):
                Ni = _N(s, i)
                if Ni:
                    m1[i] += w * Ni
                    for j in range(i, n):
                        m2[(i, j)] += w * Ni * _N(s, j)
    return m1, m2

# ------------------------------------------- (5) two-time kernel (sympy)

def twotime_kernel_poly(n):
    """Symbolic Csym_ij(s, u): E-sym[N_i(g_s) N_j(g_{su})] as sympy
    expressions polynomial in (s, u), s = e^{-g_early}, su = e^{-g_late},
    (s,u) in [0,1]^2, symmetrized over which index is read at the earlier
    time. Returns ({(i,j): expr} for i <= j, (s, u)) with s, u the sympy
    symbols. Setting u = 1 gives the single-time second moments of
    component_moments as polynomials in s. Requires sympy."""
    import sympy as sp
    s_, u_ = sp.symbols('s u', nonnegative=True)
    levels = jump_levels(n)
    co = death_semigroup_coeffs(n)
    # P_{n,k}(s): polynomial in s
    def Pnk(k):
        return sum(sp.Rational(c) * s_ ** comb(j, 2)
                   for j, c in co[n].get(k, {}).items())
    # P_{k,l}(u)
    def Pkl(k, l):
        return sum(sp.Rational(c) * u_ ** comb(j, 2)
                   for j, c in co[k].get(l, {}).items())
    out = {}
    for i in range(1, n):
        for j in range(i, n):
            acc = sp.Integer(0)
            for k in range(2, n + 1):
                pk = Pnk(k)
                for a, pa in levels[k].items():
                    Ni, Nj = _N(a, i), _N(a, j)
                    if Ni == 0 and Nj == 0:
                        continue
                    for l in range(2, k + 1):
                        pl = Pkl(k, l)
                        dl = cond_level_dist(a, l)
                        ei = sum(pp * _N(a2, i) for a2, pp in dl.items())
                        ej = sum(pp * _N(a2, j) for a2, pp in dl.items())
                        # symmetrized: (N_i(early) N_j(late) + N_j(early) N_i(late))/2
                        acc += pk * sp.Rational(pa) * pl * \
                            (sp.Rational(Ni) * sp.Rational(ej)
                             + sp.Rational(Nj) * sp.Rational(ei)) / 2
            out[(i, j)] = sp.expand(acc)
    return out, (s_, u_)

# ------------------------------------------------------------- helpers

def normalize_matrix(v, n):
    """q_ij = v_ij / sum_full where sum_full counts (i,j) and (j,i)."""
    tot = sum((1 if i == j else 2) * v[(i, j)]
              for i in range(1, n) for j in range(i, n))
    return {k: val / tot for k, val in v.items()}, tot

def flatten(q, n):
    return [q[(i, j)] for i in range(1, n) for j in range(i, n)]

# ------------------------------------------------------ identity checks

CHECK_FAMILIES = [('kingman', None), ('beta', '3/2'), ('beta', '1'),
                  ('beta', '1/2'), ('dirac', '1/4'), ('dirac', '3/4'),
                  ('msprime_dirac', '1/2,10')]

def rate_function(family, param=None):
    """lam(b, k) for family in {kingman, beta, dirac, msprime_dirac}; param a
    rational string ('3/2'), a Fraction, or 'psi,c' for msprime_dirac."""
    if family == 'kingman':
        return kingman_rate()
    if family == 'beta':
        return beta_rate(F(param))
    if family == 'dirac':
        return dirac_rate(F(param))
    if family == 'msprime_dirac':
        psi, c = (param.split(',') if isinstance(param, str) else param)
        return msprime_dirac_rate(F(psi), F(c))
    raise ValueError(f"unknown family {family!r}")

def identity_checks(n, setpartition_max_n=6):
    """Built-in exact identities at sample size n (dict name -> bool; a check
    that does not apply at this n maps to None). All comparisons are exact
    Fraction equalities.
      check1  first moments of lambda_moments == lambda_exact.expected_lengths
              for the seven CHECK_FAMILIES
      check2  lambda_moments == setpartition_moments (u and v) for Kingman,
              Beta(3/2), Dirac(1/4)  (n <= setpartition_max_n only: Bell(n)
              states)
      check3  Kingman: lambda_moments == kingman_moments_route2 (u and v)
      check4  Kingman n = 3 hand values u = (2, 1), v = [[6, 3], [3, 2]]
      check5  family limits: Beta(2) == Kingman, Dirac(0) == Kingman (u and v);
              Dirac(1) == star tree (L_1 = n T, T ~ Exp(1): u_1 = n,
              v_11 = 2 n^2, all other entries 0)
      check6  death semigroup rows are probability vectors: sum_k P_{n,k}(t)
              == 1 and every P_{n,k}(t) >= 0 at t in {1/7, 3/5, 255/256}
      check7  single-time moments: t = 1 is the star (m1_1 = n, m2_11 = n^2);
              leaf conservation on {K(g) >= 2} (the root block of size n is
              not an i < n class): sum_i i m1_i(t) == n (1 - P_{n,1}(t)) and
              sum_{i,j} i j m2_ij(t) == n^2 (1 - P_{n,1}(t)) at t in {1/512,
              3/8, 373/512, 1}, and == n, n^2 at the directional limit t = 0
              (conditioned on two blocks)"""
    out = {}
    ok = True
    for fam, par in CHECK_FAMILIES:
        lam = rate_function(fam, par)
        u, _v = lambda_moments(n, lam)
        h = expected_lengths(n, lam)[n]
        ok &= all(u[i] == h[i] for i in range(1, n))
    out['check1_first_moments_vs_recursion'] = bool(ok)
    if n <= setpartition_max_n:
        ok = True
        for fam, par in (('kingman', None), ('beta', '3/2'), ('dirac', '1/4')):
            lam = rate_function(fam, par)
            ok &= lambda_moments(n, lam) == setpartition_moments(n, lam)
        out['check2_lumped_vs_setpartition'] = bool(ok)
    else:
        out['check2_lumped_vs_setpartition'] = None
    out['check3_kingman_two_routes'] = bool(
        lambda_moments(n, kingman_rate()) == kingman_moments_route2(n))
    u3, v3 = lambda_moments(3, kingman_rate())
    out['check4_kingman_n3_hand_values'] = bool(
        u3 == {1: F(2), 2: F(1)} and v3 == {(1, 1): F(6), (1, 2): F(3),
                                            (2, 2): F(2)})
    king = lambda_moments(n, kingman_rate())
    ok = (lambda_moments(n, beta_rate(F(2))) == king
          and lambda_moments(n, dirac_rate(F(0))) == king)
    us, vs = lambda_moments(n, dirac_rate(F(1)))
    ok &= (us[1] == n and vs[(1, 1)] == 2 * n * n
           and all(us[i] == 0 for i in range(2, n))
           and all(vs[(i, j)] == 0 for i in range(1, n) for j in range(i, n)
                   if (i, j) != (1, 1)))
    out['check5_family_limits'] = bool(ok)
    ok = True
    for t in (F(1, 7), F(3, 5), F(255, 256)):
        row = [death_prob(n, n, k, t) for k in range(1, n + 1)]
        ok &= (sum(row) == 1 and all(p >= 0 for p in row))
    out['check6_death_semigroup_stochastic'] = bool(ok)
    m1, m2 = component_moments(n, F(1))
    ok = (m1[1] == n and m2[(1, 1)] == n * n
          and all(m1[i] == 0 for i in range(2, n)))
    for t in (F(0), F(1, 512), F(3, 8), F(373, 512), F(1)):
        m1, m2 = component_moments(n, t)
        live = 1 if t == 0 else 1 - death_prob(n, n, 1, t)
        ok &= sum(i * m1[i] for i in range(1, n)) == n * live
        ok &= sum((1 if i == j else 2) * i * j * m2[(i, j)]
                  for i in range(1, n) for j in range(i, n)) == n * n * live
    out['check7_component_identities'] = bool(ok)
    return out

# ---------------------------------------------------------------- main

def main(argv=None):
    ap = argparse.ArgumentParser(
        description="exact E[L_i], E[L_i L_j] and normalized 2-SFS q under a "
                    "Lambda-coalescent, plus the built-in identity checks")
    ap.add_argument('--n', type=int, default=4, help="sample size (default 4)")
    ap.add_argument('--family', default='kingman',
                    choices=['kingman', 'beta', 'dirac', 'msprime_dirac'])
    ap.add_argument('--param', default=None,
                    help="rational parameter: alpha for beta, psi for dirac, "
                         "'psi,c' for msprime_dirac (e.g. 3/2, 1/4, 1/2,10)")
    ap.add_argument('--out', default=None,
                    help="write the results JSON (exact num/den strings) here")
    a = ap.parse_args(argv)
    n = a.n
    if a.family != 'kingman' and a.param is None:
        ap.error(f"--family {a.family} needs --param")
    checks = identity_checks(n)
    failed = [k for k, v in checks.items() if v is False]
    for k, v in checks.items():
        print(f"{k:38s} {v}")
    if failed:
        print('IDENTITY CHECK FAILURE:', failed)
        return 1
    lam = rate_function(a.family, a.param)
    u, v = lambda_moments(n, lam)
    q, tot = normalize_matrix(v, n)
    tag = a.family + (f"({a.param})" if a.param is not None else '')
    print(f"n = {n}, {tag}:  E[L_tot^2] = {tot}")
    for i in range(1, n):
        print(f"  E[L_{i}] = {u[i]}  ({float(u[i]):.6f})")
    for i in range(1, n):
        for j in range(i, n):
            print(f"  E[L_{i} L_{j}] = {v[(i, j)]}   q_{i}{j} = {q[(i, j)]}"
                  f"  ({float(q[(i, j)]):.6f})")
    if a.out:
        doc = {'n': n, 'family': a.family, 'param': a.param, 'checks': checks,
               'E_L': {str(i): str(u[i]) for i in range(1, n)},
               'E_LL': {f'{i},{j}': str(v[(i, j)])
                        for i in range(1, n) for j in range(i, n)},
               'q': {f'{i},{j}': str(q[(i, j)])
                     for i in range(1, n) for j in range(i, n)},
               'E_Ltot2': str(tot)}
        with open(a.out, 'w') as fh:
            json.dump(doc, fh, indent=1)
        print('wrote', a.out)
    return 0

if __name__ == '__main__':
    sys.exit(main())
