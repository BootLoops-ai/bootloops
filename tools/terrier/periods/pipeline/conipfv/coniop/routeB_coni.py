#!/usr/bin/env python3
"""routeB_coni.py — independent D-module route to L_s (DKMM routeB pattern,
automated): theta-algebra left GB (Singular PLURAL) of the box ideal built
from the exact Gamma-shift law (shiftbox.py — independent of the card's op
strings AND of the series), staircase/rank over Q(z), normal-form connection,
monomial-curve pullback Theta = sum nu_a theta_a, first Q(s)-dependency.
h_eff = 5 card only.
Timebox: caller enforces wall cap; this script prints stage walls.
"""
import json, os, sys, subprocess, time, itertools
from fractions import Fraction as Fr
_PIPE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, _PIPE)
from family import load_card
HERE = os.path.dirname(os.path.abspath(__file__))


def shift_box_poly(card, u):
    """theta-form box for class shift u as Singular string:
    P_L(t) - z^{u+} * (P_R at t shifted) ... built EXACTLY like the verified
    coefficient identity P_L(q) c(q) = P_R(q) c(q-u):
    on F = sum c(q) z^q:  P_L(theta) F = z^u P_R(theta - u)?? — index care:
    sum_q P_L(q) c(q) z^q  and  z^u sum P_R(q) c(q-u) z^{q-u}·z^u ... both
    equal termwise, so OP = P_L(theta) - z^u P_R(theta) with P_R INDEXED AT
    the SAME q as the identity: term q of z^u{P_R(theta)F} = P_R(q-u)c(q-u)
    -> need P_R(theta) replaced by P_R evaluated at theta ... = q on z^{q-u}
    is q-u, so use P_R(theta + u) inside: OP = P_L(theta) - z^u P_R(theta+u).
    Verified numerically below before use (fail-closed)."""
    import sympy as sp
    h = card.h
    th = sp.symbols(f"t1:{h+1}")
    PL, PR = sp.Integer(1), sp.Integer(1)
    for (v, m, isn) in card.factors:
        d = sum(v[i] * u[i] for i in range(h))
        if d == 0:
            continue
        vq = sum(v[i] * th[i] for i in range(h))
        if isn:
            if d > 0:
                for j in range(d):
                    PR *= (vq - j) ** m if m > 1 else (vq - j)
            else:
                for j in range(-d):
                    PL *= (vq - d - j) ** m if m > 1 else (vq - d - j)
        else:
            if d > 0:
                for j in range(d):
                    PL *= (vq - j) ** m if m > 1 else (vq - j)
            else:
                for j in range(-d):
                    PR *= (vq - d - j) ** m if m > 1 else (vq - d - j)
    # OP acts on F: coefficient of z^q: PL(q)c(q) - PR(q)c(q-u) = 0
    shift = {th[i]: th[i] + u[i] for i in range(h)}
    PRs = sp.expand(PR.subs(shift, simultaneous=True))
    return sp.expand(PL), PRs


def op_check(card, u, PL, PRs, cfrac_ext, box):
    """numeric fail-closed check of the OPERATOR form on random q."""
    import sympy as sp
    h = card.h
    th = sp.symbols(f"t1:{h+1}")
    for q in box:
        c = cfrac_ext(card, q)
        qu = tuple(q[i] - u[i] for i in range(h))
        cu = cfrac_ext(card, qu)
        sub = dict(zip(th, q))
        subu = dict(zip(th, qu))
        lhs = Fr(str(PL.subs(sub))) * c - Fr(str(PRs.subs(subu))) * cu
        assert lhs == 0, f"operator form check fails u={u} q={q}"
