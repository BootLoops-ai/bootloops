#!/usr/bin/env python3
"""shiftbox.py — exact lattice-shift (box) identities from the Gamma law.

For any class shift u in Z^h the Gamma-series coefficients satisfy the
polynomial identity  P_L(q) c(q) = P_R(q) c(q-u)  with
  P_L(q) = prod_{num,d<0} ff(v.q-d,-d)^m * prod_{den,d>0} ff(v.q,d)^m
  P_R(q) = prod_{num,d>0} ff(v.q,d)^m * prod_{den,d<0} ff(v.q-d,-d)^m
(d = v.u, ff = falling factorial).  For u = unit class shifts these are the
GKZ boxes in theta form: OP_u = P_L(theta) - z^u P_R(theta) acting on
sum c(q) z^q  (P_R evaluated at q, the identity index).  Used to (i) verify
the card's raw boxes, (ii) pin card errata (ds-lorien op for the l0=+2
vector), (iii) provide the Route-B ideal generators independently of the
card strings.
"""
import sys, os
from fractions import Fraction as Fr
_PIPE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, _PIPE)


def ff(x, k):
    """falling factorial x(x-1)...(x-k+1) as coefficient list (ascending)?
    Here numeric eval only."""
    out = 1
    for j in range(k):
        out *= (x - j)
    return out


def shift_polys_eval(card, u, q):
    """(P_L(q), P_R(q)) exact ints for shift u at lattice point q."""
    h = card.h
    PL = PR = 1
    for (v, m, isn) in card.factors:
        d = sum(v[i] * u[i] for i in range(h))
        if d == 0:
            continue
        vq = sum(v[i] * q[i] for i in range(h))
        if isn:
            if d > 0:
                PR *= ff(vq, d) ** m
            else:
                PL *= ff(vq - d, -d) ** m
        else:
            if d > 0:
                PL *= ff(vq, d) ** m
            else:
                PR *= ff(vq - d, -d) ** m
    return PL, PR


def check_shift(card, u, box, cfrac_ext):
    nid = nontriv = 0
    h = card.h
    for q in box:
        cq = cfrac_ext(card, q)
        qu = tuple(q[i] - u[i] for i in range(h))
        cqu = cfrac_ext(card, qu)
        PL, PR = shift_polys_eval(card, u, q)
        assert PL * cq == PR * cqu, f"shift law fails u={u} q={q}"
        nid += 1
        nontriv += (cq != 0 or cqu != 0)
    return nid, nontriv


if __name__ == "__main__":
    import random
    from family import load_card
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from boxcheck import cfrac_ext
    from gseries import frame, tvec, cone_points, _one_cert
    from conipfv.coni_curve import coni_pfv_curve
    card = load_card(sys.argv[1])
    Msm = int(sys.argv[2]) if len(sys.argv) > 2 else 60
    h = card.h
    cur = coni_pfv_curve(card)
    nu = [int(x) for x in cur["nu"]]
    facs, Sinv = frame(card)
    V = [list(v) for (v, m, isn) in facs]
    nup = [int(x) for x in tvec(nu, Sinv, h)]
    certs, caps = _one_cert(V, nup, h)
    pts = cone_points(V, nup, Msm, certs, caps, h)
    sup = [tuple(sum(Sinv[i][j]*n[j] for j in range(h)) for i in range(h))
           for n, m in pts]
    random.seed(3)
    box = set()
    for q in sup:
        for _ in range(6):
            box.add(tuple(q[i] + random.randint(-1, 2) for i in range(h)))
    box = sorted(box)
    for b in range(h):
        u = tuple(1 if i == b else 0 for i in range(h))
        nid, ntv = check_shift(card, u, box, cfrac_ext)
        print(f"{card.name}: shift e_{b+1}: {nid} identities exact "
              f"({ntv} nontrivial) — PASS")
