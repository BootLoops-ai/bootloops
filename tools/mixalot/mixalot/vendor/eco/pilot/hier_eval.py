"""Lattice-DP evaluator for the 3-level formula (scalable route).

DP over species on the coupling lattice (Abar_1..Abar_R, G); per-species
transfer built from within-region across-plot convolutions (u_rs) and the
level-2 Stirling transition; (g_s-1)! twist applied at species close; final
weights prod_r 1/(rho_r)_{Abar_r} * theta^S/(theta)_G at the end.

This exact-Fraction version (no pruning) is gate G1 vs the brute-force
formula_P3. The ball version with certified pruning (NOTES_HIERARCHICAL.md)
reuses the same structure.
"""
from fractions import Fraction
from math import factorial
from collections import defaultdict, Counter
from itertools import product as iproduct

from hierarchical3 import sbar_row, pochhammer


def u_poly(nvec_plots, Ivec_plots):
    """Within-region across-plot convolution for one species:
    u(abar) = [x^abar] prod_p sum_a sbar(n_p,a) (I_p x)^a.
    Returns list of Fractions indexed by abar (0..sum n_p)."""
    u = [Fraction(1)]
    for n, I in zip(nvec_plots, Ivec_plots):
        if n == 0:
            continue
        row = sbar_row(n)
        v = [row[a] * I ** a for a in range(n + 1)]
        # note a=0 coeff is 0 for n>=1
        new = [Fraction(0)] * (len(u) + len(v) - 1)
        for i, x in enumerate(u):
            if x:
                for j, y in enumerate(v):
                    if y:
                        new[i + j] += x * y
        u = new
    return u


def hier_P3(Dmat, reg_of, theta, rhos, Ivec):
    """Exact 3-level P[D] via species-DP on the (Abar_r, G) lattice."""
    theta = Fraction(theta)
    rhos = [Fraction(x) for x in rhos]
    Ivec = [Fraction(x) for x in Ivec]
    S = len(Dmat)
    P_ = len(Dmat[0])
    R = len(rhos)
    Js = [sum(Dmat[s][p] for s in range(S)) for p in range(P_)]
    plots_of = [[p for p in range(P_) if reg_of[p] == r] for r in range(R)]

    # prefactor (I folded into u): prod_p J_p!/prod_s n_ps! / (I_p)_{J_p} / prod Phi!
    pref = Fraction(1)
    for p in range(P_):
        pref *= factorial(Js[p])
        pref /= pochhammer(Ivec[p], Js[p])
        for s in range(S):
            pref /= factorial(Dmat[s][p])
    for cnt in Counter(tuple(row) for row in Dmat).values():
        pref /= factorial(cnt)

    # per-species per-region increment lists: [(abar, b, weight)]
    def species_moves(s):
        per_region = []
        for r in range(R):
            nvec = [Dmat[s][p] for p in plots_of[r]]
            if sum(nvec) == 0:
                per_region.append([(0, 0, Fraction(1))])
                continue
            u = u_poly(nvec, [Ivec[p] for p in plots_of[r]])
            moves = []
            for abar in range(1, len(u)):
                if u[abar] == 0:
                    continue
                row = sbar_row(abar)
                for b in range(1, abar + 1):
                    moves.append((abar, b, u[abar] * row[b] * rhos[r] ** b))
            per_region.append(moves)
        return per_region

    # DP
    state = {((0,) * R, 0): Fraction(1)}
    for s in range(S):
        pr = species_moves(s)
        new = defaultdict(Fraction)
        for (Abar, G), val in state.items():
            for combo in iproduct(*pr):
                g = sum(b for _, b, _ in combo)
                if g == 0:
                    continue
                w = val * factorial(g - 1)
                for _, _, wr in combo:
                    w *= wr
                A2 = tuple(Abar[r] + combo[r][0] for r in range(R))
                new[(A2, G + g)] += w
        state = dict(new)

    total = Fraction(0)
    for (Abar, G), val in state.items():
        w = theta ** S / pochhammer(theta, G)
        for r in range(R):
            w /= pochhammer(rhos[r], Abar[r])
        total += val * w
    return pref * total
