"""Three-level neutral sampling formula (NEW) + nested process-DP oracle.

Model (3-level HDP / nested Hoppe urns): global metacommunity CRP(theta);
regions r = 1..R each with pool CRP(rho_r) whose base is the global urn;
plots p (each in region reg(p)) with local CRP(I_p) whose base is region
reg(p)'s pool.

Sequential construction, adding an individual to plot p (region r):
  local copy:          j_p/(I_p+j_p) * n_{ps}/j_p
  immigrate to plot:   I_p/(I_p+j_p), then within region r's pool:
     existing regional lineage of species s:  abar_{rs}/(rho_r + Abar_r)
        [abar_{rs} = # plot-level immigration events in region r, species s]
     new regional lineage: rho_r/(rho_r + Abar_r), then global urn:
        existing global lineage of s: g_s/(theta + G)
        new species:                  theta/(theta + G)
  [g_s = # regional-lineage creation events of species s, G = sum g_s]

CONJECTURED FORMULA (this file verifies it exactly at small sizes):

P[D] = [ prod_p J_p! / prod_{p,s} n_{ps}! / prod_vec Phi_vec! ]
     x sum over ancestry arrays {a_{ps}}, {b_{rs}}, consistent supports:
       prod_p [ I_p^{A_p} / (I_p)_{J_p} ]  * prod_{p,s} sbar(n_{ps}, a_{ps})
     x prod_r [ rho_r^{B_r} / (rho_r)_{Abar_r} ] * prod_{r,s} sbar(abar_{rs}, b_{rs})
     x theta^S / (theta)_G * prod_s (g_s - 1)!
where A_p = sum_s a_{ps}; abar_{rs} = sum_{p in r} a_{ps}; Abar_r = sum_s abar_{rs};
B_r = sum_s b_{rs}; g_s = sum_r b_{rs}; G = sum_s g_s;
1 <= a_{ps} <= n_{ps} (0 iff n_{ps}=0), 1 <= b_{rs} <= abar_{rs} (0 iff 0).

All exact Fractions; brute-force enumeration (small J only).
"""
from fractions import Fraction
from functools import lru_cache
from math import factorial
from collections import defaultdict, Counter
from itertools import product as iproduct


@lru_cache(maxsize=None)
def sbar_row(n):
    if n == 0:
        return (1,)
    prev = sbar_row(n - 1)
    row = [0] * (n + 1)
    for a in range(n):
        row[a] += (n - 1) * prev[a]
        row[a + 1] += prev[a]
    return tuple(row)


def sbar(n, a):
    if a < 0 or a > n:
        return 0
    return sbar_row(n)[a]


def pochhammer(x, k):
    x = Fraction(x)
    out = Fraction(1)
    for i in range(k):
        out *= x + i
    return out


def ranges(n):
    """valid ancestry counts for count n: [0] if n==0 else 1..n."""
    return [0] if n == 0 else list(range(1, n + 1))


def formula_P3(Dmat, reg_of, theta, rhos, Ivec):
    """Brute-force 3-level formula. Dmat[s][p] = n_{ps} (species rows).

    reg_of[p] = region index of plot p.
    """
    theta = Fraction(theta)
    rhos = [Fraction(x) for x in rhos]
    Ivec = [Fraction(x) for x in Ivec]
    S = len(Dmat)
    P_ = len(Dmat[0])
    R = len(rhos)
    Js = [sum(Dmat[s][p] for s in range(S)) for p in range(P_)]

    pref = Fraction(1)
    for p in range(P_):
        pref *= factorial(Js[p])
        for s in range(S):
            pref /= factorial(Dmat[s][p])
    for cnt in Counter(tuple(row) for row in Dmat).values():
        pref /= factorial(cnt)

    total = Fraction(0)
    # enumerate a_{ps}
    a_choices = [[ranges(Dmat[s][p]) for p in range(P_)] for s in range(S)]
    flat_a = [a_choices[s][p] for s in range(S) for p in range(P_)]
    for a_flat in iproduct(*flat_a):
        a = [[a_flat[s * P_ + p] for p in range(P_)] for s in range(S)]
        # per-plot factors
        w = Fraction(1)
        for p in range(P_):
            A_p = sum(a[s][p] for s in range(S))
            w *= Ivec[p] ** A_p / pochhammer(Ivec[p], Js[p])
            for s in range(S):
                w *= sbar(Dmat[s][p], a[s][p])
        # abar_{rs}
        abar = [[0] * S for _ in range(R)]
        for s in range(S):
            for p in range(P_):
                abar[reg_of[p]][s] += a[s][p]
        # enumerate b_{rs}
        b_lists = [ranges(abar[r][s]) for r in range(R) for s in range(S)]
        for b_flat in iproduct(*b_lists):
            b = [[b_flat[r * S + s] for s in range(S)] for r in range(R)]
            w2 = Fraction(1)
            ok = True
            for r in range(R):
                Abar_r = sum(abar[r])
                B_r = sum(b[r])
                w2 *= rhos[r] ** B_r / pochhammer(rhos[r], Abar_r) if Abar_r else 1
                for s in range(S):
                    w2 *= sbar(abar[r][s], b[r][s])
            g = [sum(b[r][s] for r in range(R)) for s in range(S)]
            if any(gs == 0 for gs in g):
                ok = False
            if not ok:
                continue
            G = sum(g)
            w3 = theta ** S / pochhammer(theta, G)
            for s in range(S):
                w3 *= factorial(g[s] - 1)
            total += w * w2 * w3
    return pref * total


def process_distribution3(Jvec, reg_of, theta, rhos, Ivec):
    """Exact nested-urn DP. State per species: (n_1..n_P, abar_1..abar_R, g)."""
    theta = Fraction(theta)
    rhos = [Fraction(x) for x in rhos]
    Ivec = [Fraction(x) for x in Ivec]
    P_ = len(Jvec)
    R = len(rhos)
    states = {(): Fraction(1)}
    for p in range(P_):
        r = reg_of[p]
        I = Ivec[p]
        rho = rhos[r]
        for j in range(Jvec[p]):
            new = defaultdict(Fraction)
            for st, prob in states.items():
                sp = [list(x) for x in st]
                G = sum(x[-1] for x in sp)
                Abar_r = sum(x[P_ + r] for x in sp)
                # local copy
                if j > 0:
                    for c, x in enumerate(sp):
                        if x[p] == 0:
                            continue
                        q = prob * Fraction(j, 1) / (I + j) * Fraction(x[p], j)
                        ns = [tuple(y) for y in sp]
                        y = list(x); y[p] += 1
                        ns[c] = tuple(y)
                        new[tuple(sorted(ns))] += q
                # immigrant -> existing regional lineage of species c
                for c, x in enumerate(sp):
                    if x[P_ + r] == 0:
                        continue
                    q = prob * I / (I + j) * Fraction(x[P_ + r]) / (rho + Abar_r)
                    ns = [tuple(y) for y in sp]
                    y = list(x); y[p] += 1; y[P_ + r] += 1
                    ns[c] = tuple(y)
                    new[tuple(sorted(ns))] += q
                # immigrant -> new regional lineage -> existing global species c
                for c, x in enumerate(sp):
                    q = (prob * I / (I + j) * rho / (rho + Abar_r)
                         * Fraction(x[-1]) / (theta + G))
                    ns = [tuple(y) for y in sp]
                    y = list(x); y[p] += 1; y[P_ + r] += 1; y[-1] += 1
                    ns[c] = tuple(y)
                    new[tuple(sorted(ns))] += q
                # immigrant -> new regional lineage -> NEW species
                q = (prob * I / (I + j) * rho / (rho + Abar_r)
                     * theta / (theta + G))
                fresh = [0] * (P_ + R + 1)
                fresh[p] = 1; fresh[P_ + r] = 1; fresh[-1] = 1
                new[tuple(sorted(st + (tuple(fresh),)))] += q
            states = dict(new)
    out = defaultdict(Fraction)
    for st, prob in states.items():
        key = tuple(sorted(tuple(x[:P_]) for x in st))
        out[key] += prob
    return dict(out)
