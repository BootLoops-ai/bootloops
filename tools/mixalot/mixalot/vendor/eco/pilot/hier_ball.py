"""Ball + certified-pruning evaluator for the 3-level formula.

Same DP as hier_eval.hier_P3 but in arb balls with rigorous cell pruning:
a cell's future contribution is bounded by
   val * finalweight(current coords) * prod_{remaining species} W_s
(final weights 1/(rho)_Abar, 1/(theta)_G are decreasing in Abar,G; W_s = total
move weight of species s incl. its (g-1)! twist — all positive). Pruned mass
is accumulated into `slack` and added to the final ball radius.

ln-space output for scale safety: returns (lnP ball, ncells stats).
"""
import math
from collections import Counter
from math import factorial

from flint import arb, ctx

from hierarchical3 import sbar_row
from itertools import product as iproduct


def hier_lnP_ball(Dmat, reg_of, theta_s, rhos_s, Ivec_s, prec=96, prune_rel=1e-28):
    ctx.prec = prec
    theta = arb(theta_s)
    rhos = [arb(x) for x in rhos_s]
    Ivec = [arb(x) for x in Ivec_s]
    S = len(Dmat)
    P_ = len(Dmat[0])
    R = len(rhos)
    Js = [sum(Dmat[s][p] for s in range(S)) for p in range(P_)]
    plots_of = [[p for p in range(P_) if reg_of[p] == r] for r in range(R)]

    # ln prefactor
    lnpref = arb(0)
    for p in range(P_):
        lnpref += arb(Js[p] + 1).lgamma()
        lnpref -= (Ivec[p] + Js[p]).lgamma() - Ivec[p].lgamma()
        for s in range(S):
            lnpref -= arb(Dmat[s][p] + 1).lgamma()
    for cnt in Counter(tuple(row) for row in Dmat).values():
        lnpref -= arb(cnt + 1).lgamma()

    # per-species moves and their total weights W_s (for pruning bounds)
    def u_poly_ball(nvec, Ivs):
        u = [arb(1)]
        for n, I in zip(nvec, Ivs):
            if n == 0:
                continue
            row = sbar_row(n)
            Ia = arb(1)
            v = []
            for a in range(n + 1):
                v.append(arb(row[a]) * Ia)
                Ia *= I
            new = [arb(0)] * (len(u) + len(v) - 1)
            for i, x in enumerate(u):
                for j, y in enumerate(v):
                    new[i + j] += x * y
            u = new
        return u

    species = []
    for s in range(S):
        per_region = []
        for r in range(R):
            nvec = [Dmat[s][p] for p in plots_of[r]]
            if sum(nvec) == 0:
                per_region.append([(0, 0, arb(1))])
                continue
            u = u_poly_ball(nvec, [Ivec[p] for p in plots_of[r]])
            moves = []
            for abar in range(1, len(u)):
                row = sbar_row(abar)
                for b in range(1, abar + 1):
                    moves.append((abar, b, u[abar] * arb(row[b]) * rhos[r] ** b))
            per_region.append(moves)
        species.append(per_region)

    # W_s: total weight incl. (g-1)! (upper bound for pruning)
    W = []
    for pr in species:
        tot = arb(0)
        for combo in iproduct(*pr):
            g = sum(b for _, b, _ in combo)
            if g == 0:
                continue
            w = arb(factorial(g - 1))
            for _, _, wr in combo:
                w *= wr
            tot += w
        W.append(tot)
    # suffix products of W upper bounds (floats of upper endpoints, ln-space)
    lnWtail = [0.0] * (S + 1)
    for s in range(S - 1, -1, -1):
        hi = float(W[s].mid()) + float(W[s].rad())
        lnWtail[s] = lnWtail[s + 1] + math.log(max(hi, 1e-300))

    def ln_finalweight(Abar, G):
        # ln of prod_r 1/(rho_r)_{Abar_r} * theta^S/(theta)_G  (decreasing in coords)
        w = S * theta.log() - ((theta + G).lgamma() - theta.lgamma())
        for r in range(R):
            w -= (rhos[r] + Abar[r]).lgamma() - rhos[r].lgamma()
        return float(w.mid())

    # DP with pruning
    state = {((0,) * R, 0): arb(1)}
    slack = 0.0  # accumulated ln-space-safe absolute bound on discarded mass
    for s in range(S):
        pr = species[s]
        new = {}
        for (Abar, G), val in state.items():
            for combo in iproduct(*pr):
                g = sum(b for _, b, _ in combo)
                if g == 0:
                    continue
                w = val * arb(factorial(g - 1))
                for _, _, wr in combo:
                    w *= wr
                A2 = tuple(Abar[r] + combo[r][0] for r in range(R))
                key = (A2, G + g)
                if key in new:
                    new[key] += w
                else:
                    new[key] = w
        # prune
        if new:
            # ln of the largest cell*finalweight as scale reference
            best = max(math.log(max(float(v.mid()), 1e-300)) + ln_finalweight(k[0], k[1])
                       for k, v in new.items())
            thresh = best + math.log(prune_rel)
            kept = {}
            for k, v in new.items():
                bound = (math.log(max(float(v.mid()) + float(v.rad()), 1e-300))
                         + ln_finalweight(k[0], k[1]) + lnWtail[s + 1])
                if bound < thresh:
                    slack += math.exp(bound)
                else:
                    kept[k] = v
            state = kept
        else:
            state = new
    # final assembly (work relative to lnpref + best scale)
    tot = arb(0)
    for (Abar, G), val in state.items():
        w = theta ** S / _poch(theta, G)
        for r in range(R):
            w /= _poch(rhos[r], Abar[r])
        tot += val * w
    tot += arb(0, slack)  # certified slack for pruned mass
    return lnpref + tot.log(), {'cells_final': len(state)}


def _poch(x, k):
    out = arb(1)
    for i in range(k):
        out *= x + i
    return out
