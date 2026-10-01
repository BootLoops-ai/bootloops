#!/usr/bin/env python3
"""
coni_curve.py — coni-PFV restriction chart + g2' gate (CP0-CP6), exact.

Lemma (2406.13751 SS3.2, eqs eq:detN_conition..eq:integrality_condition_II):
  N = kappa.M;  p = N^-1 K + lam N^-1 q_cf with lam fixed by q_cf.p = 0:
      lam = -(q_cf.N^-1 K)/(q_cf.N^-1 q_cf)
  gates: K.p = 0 exact (Diophantine); N.p - K = lam q_cf exact re-check;
  facet interiority vs the pinned GV window: q.nu >= 0, = 0 iff q || q_cf.
All arithmetic Fraction-exact; see conipfv/DESIGN.md SS1/SS5.
"""
import sys, os
from fractions import Fraction as Fr
from math import gcd

_PIPE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _PIPE)
from curve_from_flux import solve_exact  # exact Gaussian solve


def _is_parallel(q, qcf):
    """q = k*qcf for positive integer k?"""
    ks = set()
    for a, b in zip(q, qcf):
        if b == 0:
            if a != 0:
                return False
        else:
            ks.add(Fr(a, b))
    return len(ks) == 1 and ks.pop() >= 1


def _primitive(q):
    g = 0
    for x in q:
        g = gcd(g, abs(int(x)))
    return g == 1


def _hnf_kernel_basis(qcf):
    """integer basis of {v in Z^h : qcf.v = 0} (h-1 rows), via column ops."""
    h = len(qcf)
    U = [[1 if i == j else 0 for j in range(h)] for i in range(h)]
    r = [int(x) for x in qcf]
    # reduce r to (g,0,...,0) tracking unimodular column ops on U
    while True:
        nz = [i for i in range(h) if r[i] != 0]
        if len(nz) <= 1:
            break
        nz.sort(key=lambda i: abs(r[i]))
        i0 = nz[0]
        for i in nz[1:]:
            f = r[i] // r[i0]
            r[i] -= f * r[i0]
            for a in range(h):
                U[a][i] -= f * U[a][i0]
    assert any(r), "qcf = 0"
    piv = next(i for i in range(h) if r[i] != 0)
    return [[U[a][i] for a in range(h)] for i in range(h) if i != piv]


def _solve_candidate(N, K, q, h):
    """p(lam) for candidate conifold class q; None if degenerate."""
    try:
        NiK = solve_exact([r[:] for r in N], [Fr(x) for x in K])
        Niq = solve_exact([r[:] for r in N], [Fr(x) for x in q])
    except StopIteration:
        return None                                   # det N = 0
    den = sum(Fr(q[a]) * Niq[a] for a in range(h))    # q.N^-1 q
    if den == 0:
        return None
    lam = -sum(Fr(q[a]) * NiK[a] for a in range(h)) / den
    p = [NiK[a] + lam * Niq[a] for a in range(h)]
    return p, lam


def coni_pfv_curve(card):
    """g2': CP0-CP6 (DESIGN SS5). Returns dict; AssertionError = gate FAIL."""
    h, M, K = card.h, card.M, card.K
    N = [[sum(card.kap(i, j, c) * M[c] for c in range(h)) for j in range(h)]
         for i in range(h)]
    # CP0 routing: PFFV cards are NOT ours
    try:
        NiK = solve_exact([r[:] for r in N], K[:])
    except StopIteration:
        raise AssertionError("CP0: det N = 0 (eq:detN violated)")
    KNiK = sum(K[a] * NiK[a] for a in range(h))
    if KNiK == 0:
        return dict(routing="NOT-CONI",
                    note="K.N^-1.K = 0: PFFV card, defer to g2/pffv_curve")
    # CP1: unique primitive pinned class surviving lemma + CP4
    survivors = []
    for q in sorted(card.gv_pinned):
        if not _primitive(q):
            continue
        sol = _solve_candidate(N, K, q, h)
        if sol is None:
            continue
        p, lam = sol
        if sum(K[a] * p[a] for a in range(h)) != 0:
            continue                                  # Diophantine gate
        degs = {qq: sum(Fr(qq[a]) * p[a] for a in range(h))
                for qq in card.gv_pinned}
        if any(d < 0 for d in degs.values()):
            continue
        if any(d == 0 and not _is_parallel(qq, q) for qq, d in degs.items()):
            continue
        if degs[q] != 0 or all(d == 0 for d in degs.values()):
            continue
        survivors.append((q, p, lam))
    assert len(survivors) == 1, \
        f"CP1: {len(survivors)} conifold candidates survive (need exactly 1)"
    q_cf, p, lam = survivors[0]
    # CP2: exact lemma residual re-check
    Np = [sum(N[a][b] * p[b] for b in range(h)) for a in range(h)]
    for a in range(h):
        assert Np[a] - K[a] == lam * q_cf[a], "CP2: N.p - K != lam*q_cf"
    # CP3: conifold multiplicities + D4 validity flags
    n_cf = card.gv_pinned[q_cf]
    assert n_cf >= 1, "CP3: GV(q_cf) < 1"
    M_cf = sum(Fr(q_cf[a]) * M[a] for a in range(h))
    assert M_cf != 0, "CP3: M_cf = q_cf.M = 0 (no flux on shrinking cycle)"
    K_cf = sum(Fr(q_cf[a]) * K[a] for a in range(h))
    # CP4: nu = r p primitive integer, degrees integral
    r = 1
    for x in p:
        r = r * x.denominator // gcd(r, x.denominator)
    nu = [int(x * r) for x in p]
    g = 0
    for x in nu:
        g = gcd(g, abs(x))
    assert g == 1, "CP4: nu not primitive"
    degs = {q: sum(q[a] * nu[a] for a in range(h)) for q in card.gv_pinned}
    assert all(isinstance(d, int) or d.denominator == 1
               for d in degs.values()), "CP4: non-integer degree"
    # CP5: flux integrality on the bulk projection (existence of P)
    aM = [sum(card.a_mat[i][j] * M[j] for j in range(h)) for i in range(h)]
    assert all(x.denominator == 1 for x in aM), "CP5: a.M not integral"
    # evenness of bulk (A.M) (eq:integrality I) is REPORTED, not fail-closed:
    # the a-representative is frame-dependent mod 2 (freeze B2 residue,
    # PIN-REQ); the card frame need not be the 2406 coni frame.  Parity is
    # lattice-linear, so any integer kernel basis gives one verdict per frame.
    W = _hnf_kernel_basis(q_cf)          # integer basis of q_cf-perp
    bulk_even = all(
        (lambda v: v.denominator == 1 and int(v) % 2 == 0)(
            sum(Fr(w[a]) * aM[a] for a in range(h))) for w in W)
    c2M = sum(Fr(card.c2D[a]) * M[a] for a in range(h))
    assert c2M % 24 == 0, "CP5: c2'.M not in 24Z (eq:integrality II; " \
        "card c2D must be the conifold-SHIFTED c2', freeze B3)"
    # CP6: tadpole + symplectic flux vectors (C1 ordering)
    tad = Fr(-1, 2) * sum(M[a] * K[a] for a in range(h))
    assert tad.denominator == 1 and tad <= card.Q_D3, \
        f"CP6: tadpole {tad} > Q_D3 {card.Q_D3}"
    bM = sum(card.b_vec[i] * M[i] for i in range(h))
    assert bM.denominator == 1, "CP6: b.M not integral"
    F = [int(bM)] + [int(x) for x in aM] + [0] + [int(x) for x in M]
    H = [0] + [int(x) for x in K] + [0] * (h + 1)
    out = dict(routing="CONI-PFV", N=N, p=p, lam=lam, q_cf=list(q_cf),
               n_cf=int(n_cf), M_cf=M_cf, K_cf=K_cf, r=r, nu=nu,
               degs={",".join(map(str, q)): int(d) for q, d in degs.items()},
               KNiK=KNiK, tadpole=int(tad), F=F, H=H,
               bulk_aM_even=bulk_even,
               d4_flags=dict(same_sign=(M_cf * K_cf > 0),
                             K_over_M=str(Fr(K_cf, M_cf) if M_cf else None)))
    _crosswalk_gate(card, out)
    return out


def _crosswalk_gate(card, out):
    """When the card banks coni_basis_crosswalk, match EXACTLY through T."""
    cw = card.raw.get("coni_basis_crosswalk")
    out["crosswalk"] = None
    if not cw:
        return
    h = card.h
    chk = {}
    T = cw.get("T_ours_to_coni")
    fd = cw.get("flat_direction_exact") or {}
    cc = cw.get("conifold_curve") or {}
    if T and cc.get("ours"):
        qc = [sum(T[i][j] * cc["ours"][j] for j in range(h)) for i in range(h)]
        chk["q_cf_coni_map"] = (qc == list(cc.get("coni", [])))
        chk["q_cf_match"] = (list(cc["ours"]) == out["q_cf"])
    if cc.get("n_cf") is not None:
        chk["n_cf_match"] = (int(cc["n_cf"]) == out["n_cf"])
    if T and fd.get("nu_coni") is not None:
        # z-type (contravariant) vectors: nu_ours = T^T nu_coni
        nuc = [int(x) for x in fd["nu_coni"]]
        back = [sum(T[i][j] * nuc[i] for i in range(h)) for j in range(h)]
        chk["nu_coni_match"] = (back == out["nu"])
    if fd.get("d") is not None:
        chk["r_match"] = (int(fd["d"]) == out["r"])
    if cw.get("M_coni") and T:
        # curve classes map by T; flux M is divisor-indexed: pair-invariance
        # gate instead: q_cf.M equal in both frames
        Mc = [Fr(x) for x in cw["M_coni"]]
        qcfc = [Fr(x) for x in cc.get("coni", [0] * h)]
        chk["M_cf_pairing"] = (sum(a * b for a, b in zip(qcfc, Mc))
                               == out["M_cf"])
    mism = fd.get("conifold_flux_mismatch", "")
    if mism:
        chk["lam_note"] = f"card note: {mism}; derived lam = {out['lam']}"
    bad = [k for k, v in chk.items() if v is False]
    assert not bad, f"crosswalk gate FAIL: {bad}"
    out["crosswalk"] = chk


if __name__ == "__main__":
    import json
    from family import load_card
    card = load_card(sys.argv[1])
    out = coni_pfv_curve(card)
    if out["routing"] == "NOT-CONI":
        print(f"{card.name}: NOT-CONI (PFFV card, defer to g2)")
    else:
        print(f"{card.name}: CONI-PFV  K.N^-1.K = {out['KNiK']}")
        print(f"  q_cf = {out['q_cf']}  n_cf = {out['n_cf']}  "
              f"M_cf = {out['M_cf']}  K_cf = {out['K_cf']}  lam = {out['lam']}")
        print(f"  p = {[str(x) for x in out['p']]}  r = {out['r']}")
        print(f"  nu = {out['nu']}  tadpole = {out['tadpole']} <= {card.Q_D3}")
        print(f"  crosswalk = {out['crosswalk']}")
