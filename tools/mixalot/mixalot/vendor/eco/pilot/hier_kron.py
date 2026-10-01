"""Kronecker-encoded certified evaluator for the 3-level formula (R=2).

Lattice (Abar_1, Abar_2, G) embedded in one variable:
   index = Abar_1 + B1*Abar_2 + B1*B2*G,  B1 = J1+1, B2 = J2+1, BG = J1+J2+1.
Per-species transfer = sparse arb_poly; the whole DP = balanced product tree
(FLINT C-speed, certified balls, no pruning). Final: coefficient sweep with
precomputed Pochhammer weight ladders.

Feasible for region totals J_r up to a few hundred (memory ~ B1*B2*BG coeffs).
No wraparound: Abar_r <= J_r, G <= J1+J2 by construction.
"""
import math
from collections import Counter
from math import factorial

from flint import arb, arb_poly, ctx

from hierarchical3 import sbar_row


def _u_poly_ball(nvec, Ivs, prec):
    ctx.prec = prec
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


def hier_lnP_kron(Dmat, reg_of, theta_s, rhos_s, Ivec_s, prec=96):
    """Certified ln P for R=2 regions. Returns (arb, info)."""
    ctx.prec = prec
    theta = arb(theta_s)
    rhos = [arb(x) for x in rhos_s]
    Ivec = [arb(x) for x in Ivec_s]
    S = len(Dmat)
    P_ = len(Dmat[0])
    R = 2
    assert len(rhos) == 2
    Js = [sum(Dmat[s][p] for s in range(S)) for p in range(P_)]
    plots_of = [[p for p in range(P_) if reg_of[p] == r] for r in range(R)]
    Jr = [sum(Js[p] for p in plots_of[r]) for r in range(R)]
    B1, B2 = Jr[0] + 1, Jr[1] + 1
    BG = Jr[0] + Jr[1] + 1
    length = B1 * B2 * BG

    # ln prefactor
    lnpref = arb(0)
    for p in range(P_):
        lnpref += arb(Js[p] + 1).lgamma()
        lnpref -= (Ivec[p] + Js[p]).lgamma() - Ivec[p].lgamma()
        for s in range(S):
            lnpref -= arb(Dmat[s][p] + 1).lgamma()
    for cnt in Counter(tuple(row) for row in Dmat).values():
        lnpref -= arb(cnt + 1).lgamma()

    # scale guard: work with weights scaled by exp(-mu) per species to keep
    # coefficients in arb-friendly range (they can be enormous); we scale each
    # species poly by 1/max_coeff and track ln(scale) exactly.
    ln_scale = arb(0)
    polys = []
    for s in range(S):
        # per-region (abar, b) move lists
        moves_r = []
        for r in range(R):
            nvec = [Dmat[s][p] for p in plots_of[r]]
            if sum(nvec) == 0:
                moves_r.append([(0, 0, arb(1))])
                continue
            u = _u_poly_ball(nvec, [Ivec[p] for p in plots_of[r]], prec)
            mv = []
            for abar in range(1, len(u)):
                row = sbar_row(abar)
                for b in range(1, abar + 1):
                    mv.append((abar, b, u[abar] * arb(row[b]) * rhos[r] ** b))
            moves_r.append(mv)
        # combine regions with (g-1)! twist -> sparse dict index -> weight
        entries = {}
        for a1, b1_, w1 in moves_r[0]:
            for a2, b2_, w2 in moves_r[1]:
                g = b1_ + b2_
                if g == 0:
                    continue
                w = w1 * w2 * arb(factorial(g - 1))
                idx = a1 + B1 * a2 + B1 * B2 * g
                entries[idx] = entries.get(idx, arb(0)) + w
        # scale
        mx = max(float(v.mid()) for v in entries.values())
        sc = arb(mx)
        ln_scale += sc.log()
        coeffs = [arb(0)] * (max(entries) + 1)
        for idx, v in entries.items():
            coeffs[idx] = v / sc
        polys.append(arb_poly(coeffs))

    # balanced product tree
    while len(polys) > 1:
        polys = [polys[i] * polys[i + 1] for i in range(0, len(polys) - 1, 2)] \
                + ([polys[-1]] if len(polys) % 2 else [])
    M = polys[0]
    assert M.length() <= length + 1, (M.length(), length)

    # weight ladders: 1/(rho1)_{A1}, 1/(rho2)_{A2}, theta^S/(theta)_G
    inv_p1 = [arb(1)]
    for k in range(Jr[0]):
        inv_p1.append(inv_p1[-1] / (rhos[0] + k))
    inv_p2 = [arb(1)]
    for k in range(Jr[1]):
        inv_p2.append(inv_p2[-1] / (rhos[1] + k))
    inv_pg = [arb(1)]
    for k in range(Jr[0] + Jr[1]):
        inv_pg.append(inv_pg[-1] / (theta + k))

    # separable-weight contraction via convolutions:
    # sum_idx c[A1 + B1 A2 + B1B2 G] w1(A1) w2(A2) wg(G)
    #   stage 1: multiply by reversed-w1 poly, read stride-B1 at offset B1-1
    W1rev = arb_poly(list(reversed(inv_p1[:B1])))
    M1 = M * W1rev
    L = M.length()
    n2g = B2 * BG
    vals2 = [arb(0)] * n2g
    for j in range(n2g):
        pos = j * B1 + (B1 - 1)
        if pos < M1.length():
            vals2[j] = M1[pos]
    # stage 2: contract A2 (vals2 indexed j = A2 + B2*G)
    W2rev = arb_poly(list(reversed(inv_p2[:B2])))
    M2 = arb_poly(vals2) * W2rev
    tot = arb(0)
    for G in range(BG):
        pos = G * B2 + (B2 - 1)
        if pos < M2.length():
            tot += M2[pos] * inv_pg[G]
    tot *= theta ** S
    return lnpref + ln_scale + tot.log(), {'poly_len': L, 'S': S, 'Jr': Jr}
