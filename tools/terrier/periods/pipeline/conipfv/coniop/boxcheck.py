#!/usr/bin/env python3
"""boxcheck.py — r1-analog: the card's RAW GKZ box operators annihilate the
Gamma coefficients c(q) EXACTLY on a lattice box (coefficient identities
sum_mu P_mu(q - mu) c(q - mu) = 0 per q).  Validates the coniop coefficient
law + frame against the card's own (independently derived) GKZ ideal.
Operators parsed from card raw: gkz_ideal_raw_theta (lorien schema) or
gkz_box_ideal_UNGATED (the other conifold-PFV card schema) — sympy in
th1..thh, z1..zh.
"""
import sys, os, json, random
from fractions import Fraction as Fr
_PIPE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, _PIPE)
from family import load_card
from geff_series import cfrac


def parse_ops(card):
    import sympy as sp
    h = card.h
    th = sp.symbols(f"th1:{h+1}")
    zz = sp.symbols(f"z1:{h+1}")
    raw = card.raw.get("gkz_ideal_raw_theta") or \
        card.raw.get("gkz_box_ideal_UNGATED")
    assert raw, "card carries no GKZ box ideal"
    ops = []
    for s in raw:
        e = sp.expand(sp.sympify(s, dict([(f"th{i+1}", th[i]) for i in range(h)]
                                         + [(f"z{i+1}", zz[i]) for i in range(h)])))
        terms = {}                      # mu (z-exponent) -> poly-in-theta dict
        for mono, coef in e.as_poly(*zz, domain=f"QQ[{','.join(map(str,th))}]"
                                    ).terms():
            terms[tuple(int(x) for x in mono)] = sp.Poly(coef, *th)
        ops.append(terms)
    return ops, th


def cfrac_ext(card, q):
    """Gamma coefficient on the FULL class lattice, cfrac negative-arg law:
    den zeros strictly dominating num poles -> exact 0; undominated num pole
    -> fail-closed; else num!/prod den!  (matches gseries.coeff on support)."""
    from math import factorial
    num, den, npole, nzero = 1, 1, 0, 0
    for (v, m, isn) in card.factors:
        a = sum(v[i] * q[i] for i in range(card.h))
        if a < 0:
            if isn:
                npole += m
            else:
                nzero += m
        elif isn:
            num *= factorial(a) ** m
        else:
            den *= factorial(a) ** m
    if nzero > npole:
        return Fr(0)
    assert npole == 0, f"undominated numerator Gamma pole at q={q}"
    return Fr(num, den)


def apply_ops(card, ops, box, tag):
    """check sum_mu P_mu(q-mu) c(q-mu) = 0 for all q in box, per op (exact)."""
    h = card.h
    nid = nontriv = 0
    for q in box:
        for k, terms in enumerate(ops):
            tot = Fr(0); seen = 0
            for mu, P in terms.items():
                qq = tuple(q[i] - mu[i] for i in range(h))
                c = cfrac_ext(card, qq)
                if c == 0:
                    continue
                seen += 1
                val = P.eval(dict(zip(P.gens, qq)))
                tot += c * Fr(val)
            assert tot == 0, f"{tag}: box identity FAILS at q={q}, op {k}"
            nid += 1
            nontriv += (seen > 0)
    return nid, nontriv


if __name__ == "__main__":
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from gseries import frame, tvec, cone_points, _one_cert
    from conipfv.coni_curve import coni_pfv_curve
    card = load_card(sys.argv[1])
    Msm = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    cur = coni_pfv_curve(card)
    nu = [int(x) for x in cur["nu"]]
    h = card.h
    facs, Sinv = frame(card)
    V = [list(v) for (v, m, isn) in facs]
    nup = [int(x) for x in tvec(nu, Sinv, h)]
    certs, caps = _one_cert(V, nup, h)
    pts = cone_points(V, nup, Msm, certs, caps, h)
    # card-frame support points q = Sinv.n'
    sup = [tuple(sum(Sinv[i][j]*n[j] for j in range(h)) for i in range(h))
           for n, m in pts]
    random.seed(11)
    box = set()
    for q in sup:
        for _ in range(6):
            box.add(tuple(q[i] + random.randint(-1, 2) for i in range(h)))
    ops, th = parse_ops(card)
    nid, nontriv = apply_ops(card, ops, sorted(box), card.name)
    print(f"{card.name}: BOXCHECK PASS — {nid} exact identities "
          f"({nontriv} nontrivial) at {len(box)} lattice points x {len(ops)} ops "
          f"(support seed M<={Msm}, {len(sup)} pts)")
