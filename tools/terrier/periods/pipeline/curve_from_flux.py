#!/usr/bin/env python3
"""
curve_from_flux.py — PFFV flat direction -> monomial restriction curve,
general h_eff (exact fractions).

  N_ab = kappa_abc M^c        (PFFV lemma matrix; det N != 0 required)
  p    = N^{-1} K             (flat direction; K.p = 0 and p > 0 asserted)
  d    = lcm of denominators of p;  nu = d p  (integer, primitive)
  curve of record:  z~_a = s^{nu_a}   (tilde = sigma-twisted GV-positive frame)
  flat coordinate:  T = p tau,  q_s = e^{2 pi i tau / d} = prod q_a^{c_a},
                    c.nu = 1 (Bezout)
  flux vectors (convention C1: Pi-ordering (F_0, F_a, X^0, X^a)):
      F = (M.b, (a.M)^T, 0, M^T),   H = (0, K^T, 0, 0)
  tadpole:  N_flux = -M.K/2 <= Q_D3
  Half-integral tadpole policy (convention C2):
      default = integer-tadpole gate.  Card flag allow_half_tadpole
      (the paper's true -M.K/2 value as a pin + a provenance string)
      enables CONDITIONAL-C2 mode: tad in (1/2)Z accepted iff
      0 <= tad <= Q_D3, tad == pin, and N_D3 = Q_D3 - tad - 1/2 in Z>=0;
      the gate record then carries c2_live: half-D3 and a verdict-scope
      string, inherited by every downstream verdict.
GATE (run_pipe G2): DKMM card -> p = (2/5, 3/10), nu = (4,3) (curve z = (s^4,
s^3)), F = (7,3,-24,0,-16,50), H = (0,3,-4,0,0,0), tadpole 124 <= 138.
"""
from fractions import Fraction as Fr
from math import gcd


def solve_exact(A, b):
    """exact Gaussian solve A x = b (A nonsingular, Fractions)."""
    n = len(A)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[p] = M[p], M[c]
        pv = M[c][c]
        M[c] = [x / pv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [M[r][t] - f * M[c][t] for t in range(n + 1)]
    return [M[i][n] for i in range(n)]


def pffv_curve(card):
    h = card.h
    M, K = card.M, card.K
    N = [[sum(card.kap(i, j, c) * M[c] for c in range(h)) for j in range(h)]
         for i in range(h)]
    p = solve_exact([row[:] for row in N], K[:])
    assert sum(K[a] * p[a] for a in range(h)) == 0, "K.N^-1.K != 0"
    assert all(x > 0 for x in p), "flat direction p not in the Kahler cone"
    d = 1
    for x in p:
        d = d * x.denominator // gcd(d, x.denominator)
    nu = [int(x * d) for x in p]
    g = 0
    for x in nu:
        g = gcd(g, x)
    assert g == 1, "nu not primitive (d not minimal)"
    tad = Fr(-1, 2) * sum(M[a] * K[a] for a in range(h))
    c2 = getattr(card, "allow_half_tadpole", None)
    c2_live = None
    if c2 is None:
        # stock gate (integer tadpole bookkeeping).  Integrality and budget
        # are asserted separately so the message names the failing branch.
        assert tad.denominator == 1, \
            f"tadpole {tad} half-integral (convention C2: needs the card's " \
            f"allow_half_tadpole pin of the paper's true -M.K/2)"
        assert tad <= card.Q_D3, f"tadpole {tad} > Q_D3 {card.Q_D3}"
    else:
        # CONDITIONAL-C2 mode (silent acceptance prohibited — c2_live is
        # written into the gate record, and every downstream verdict row
        # inherits the verdict_scope wording).
        assert tad.denominator in (1, 2), f"tadpole {tad} not in (1/2)Z"
        assert Fr(0) <= tad <= card.Q_D3, \
            f"tadpole {tad} outside [0, Q_D3={card.Q_D3}]"
        assert tad == c2["pin"], \
            f"tadpole {tad} != card pin {c2['pin']} ({c2['provenance']})"
        if tad.denominator == 2:
            nd3 = card.Q_D3 - tad - Fr(1, 2)
            assert nd3.denominator == 1 and nd3 >= 0, \
                f"N_D3 = Q_D3 - tad - 1/2 = {nd3} not in Z>=0"
            c2_live = {"c2_live": "half-D3", "N_D3": int(nd3),
                       "provenance": c2["provenance"],
                       "verdict_scope":
                       "CONDITIONAL on C2: source papers' flux-quantization "
                       "convention (odd integer quanta; additionally one "
                       "stuck half D3-brane for this vacuum, Q_flux = "
                       f"{tad} in Z+1/2)"}
    aM = [sum(card.a_mat[i][j] * M[j] for j in range(h)) for i in range(h)]
    bM = sum(card.b_vec[i] * M[i] for i in range(h))
    assert all(x.denominator == 1 for x in aM) and bM.denominator == 1, \
        "flux vector integrality fails (a.M or b.M non-integer)"
    F = [bM] + aM + [Fr(0)] + list(M)
    H = [Fr(0)] + list(K) + [Fr(0)] * (h + 1)
    assert all(x.denominator == 1 for x in F + H)
    out = dict(N=N, p=p, d=d, nu=tuple(nu),
               tadpole=int(tad) if tad.denominator == 1 else str(tad),
               F=[int(x) for x in F], H=[int(x) for x in H])
    if c2_live is not None:
        out.update(c2_live)
    return out


if __name__ == "__main__":
    import sys, os, json
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from family import load_card
    card = load_card(sys.argv[1])
    out = pffv_curve(card)
    print("N  =", out["N"])
    print("p  =", out["p"], " d =", out["d"], " nu =", out["nu"])
    print("F  =", out["F"], "\nH  =", out["H"], "\ntadpole =", out["tadpole"],
          "<= Q_D3 =", card.Q_D3)
