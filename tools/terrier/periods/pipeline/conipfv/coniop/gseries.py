#!/usr/bin/env python3
"""gseries.py — class-graded restricted fundamental period for CONI-PFV cards.

Frame-free support law: the Gamma-series support is {q in Z^h : ALL factor
rows (num and den) >= 0} — the card-basis positive orthant is NOT part of the
definition (ds-lorien's support cone exits the orthant in every coordinate;
the orthant restriction silently truncates the series).  We enumerate in a
UNIMODULAR ROW FRAME: h factor rows S with det +-1; n' = S q >= 0 on the
support, all other rows transform integrally, coefficients depend on row
VALUES only.  Level grading m = q.nu is proven finite-per-level by exact
LP certificates nu' = alpha + beta.V' (alpha,beta >= 0), which also cap the
DFS.  Coefficient law = geff_series.cfrac's (num!/prod den!, all rows >= 0
on the enumerated set by construction).  Conifold ray: q_cf support-checked
OUT (D5 resummation never enters the z_cf = 0 graded series).
"""
import sys, os, json, time
from math import factorial, lcm
from fractions import Fraction as Fr
_PIPE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, _PIPE)


def frame(card):
    """(facs', Sinv): factor rows in a unimodular row frame (support in n'>=0).
    facs' = list of (row', mult, is_num); q = Sinv.n' (integer matrix)."""
    import itertools, numpy as np
    from sympy import Matrix
    h = card.h
    rows = [list(v) for (v, m, isn) in card.factors]
    A = np.array(rows)
    S = None
    for sub in itertools.combinations(range(len(rows)), h):
        if abs(round(np.linalg.det(A[list(sub)].astype(float)))) == 1:
            S = Matrix([rows[i] for i in sub])
            if abs(S.det()) == 1:
                break
            S = None
    assert S is not None, "no unimodular row frame (support lattice not Z^h)"
    Sinv = S.inv()
    assert all(x == int(x) for x in Sinv), "Sinv not integral"
    Sinv = [[int(Sinv[i, j]) for j in range(h)] for i in range(h)]
    facs = []
    for (v, m, isn) in card.factors:
        vp = tuple(sum(v[i] * Sinv[i][j] for i in range(h)) for j in range(h))
        facs.append((vp, m, isn))
    return facs, Sinv


def tvec(vec, Sinv, h):
    """transform a card-frame covector to the row frame: vec.q = vec'.n'."""
    return [sum(Fr(vec[i]) * Sinv[i][j] for i in range(h)) for j in range(h)]


def coeff(facs, vals):
    """c at row values vals (aligned with facs); all vals >= 0 guaranteed."""
    num, den = 1, 1
    for (v, m, isn), a in zip(facs, vals):
        if isn:
            num *= factorial(a) ** m
        else:
            den *= factorial(a) ** m
    return Fr(num, den)


def _one_cert(V, nu, h):
    """exact alpha >= 0: nu' = alpha + beta.V (beta >= 0), maximize alpha.obj
    for each obj in {1} u {e_a}; returns list of exact alpha certs."""
    from scipy.optimize import linprog
    nv = len(V)
    outs = []
    for obj in [[1]*h] + [[1 if i == a else 0 for i in range(h)] for a in range(h)]:
        res = linprog([float(sum(obj[a]*V[j][a] for a in range(h)))
                       for j in range(nv)],
                      A_ub=[[float(V[j][a]) for j in range(nv)] for a in range(h)],
                      b_ub=[float(x) for x in nu], bounds=[(0, 1e5)]*nv,
                      method="highs")
        if not res.success:
            continue
        beta = [Fr(max(0.0, x)).limit_denominator(10**7) for x in res.x]
        for _ in range(80):
            alpha = [Fr(nu[a]) - sum(beta[j]*V[j][a] for j in range(nv))
                     for a in range(h)]
            if all(x >= 0 for x in alpha):
                outs.append(alpha); break
            beta = [b * Fr(999, 1000) for b in beta]
    assert outs, "no grading certificate: nu not in dual support cone"
    caps_a = [max(c[a] for c in outs) for a in range(h)]
    assert all(x > 0 for x in caps_a), f"non-strict grading: caps {caps_a}"
    return outs, caps_a


def cone_points(V, nu, M, certs, caps_a, h, visit=None):
    """all n' >= 0 with V.n' >= 0 and n'.nu <= M; exact, certificate-capped.
    visit(n_tuple, m): streaming mode (no point list kept; returns count)."""
    order = sorted(range(h), key=lambda a: -caps_a[a])
    NU = [nu[a] for a in order]
    VR = [[v[a] for a in order] for v in V]
    capo = [int(Fr(M) / caps_a[order[d]]) for d in range(h)]
    L = lcm(*[c[a].denominator for c in certs for a in range(h)])
    W = [[int(certs[c][order[d]] * L) for d in range(h)]
         for c in range(len(certs))]
    ML, nc, nV = M * L, len(certs), len(V)
    maxrow = [[sum(max(0, VR[j][i]) * capo[i] for i in range(d, h))
               for d in range(h + 1)] for j in range(nV)]
    mindeg = [sum(min(0, NU[i]) * capo[i] for i in range(d, h))
              for d in range(h + 1)]
    out = []
    cnt = [0]
    def rec(d, n, pc, pr, pm):
        if d == h:
            if all(x >= 0 for x in pr) and 0 <= pm <= M:
                nn = [0] * h
                for i in range(h):
                    nn[order[i]] = n[i]
                cnt[0] += 1
                if visit is None:
                    out.append((tuple(nn), pm))
                else:
                    visit(nn, pm)
            return
        for na in range(capo[d] + 1):
            npc = [pc[c] + W[c][d] * na for c in range(nc)]
            if any(x > ML for x in npc):
                break                      # cert partials nondecreasing
            bad = next((j for j in range(nV)
                        if pr[j] + VR[j][d]*na + maxrow[j][d+1] < 0), None)
            if bad is not None:
                if VR[bad][d] < 0:
                    break
                continue
            npm = pm + NU[d] * na
            if npm + mindeg[d + 1] > M:
                continue
            n[d] = na
            rec(d + 1, n, npc,
                [pr[j] + VR[j][d]*na for j in range(nV)], npm)
        n[d] = 0
    rec(0, [0]*h, [0]*nc, [0]*nV, 0)
    return out if visit is None else cnt[0]


def graded_series(card, nu, M, sigma=None):
    """a_m, J1_m, J2_m (m = 0..M, exact): graded period + z_cf-jet gradings
    J_k = sum (q.M)^k c(q)  (q.xi = (q.M)/M_cf, DESIGN SS1)."""
    h = card.h
    facs, Sinv = frame(card)
    V = [list(v) for (v, m, isn) in facs]         # ALL rows >= 0 = support
    nup = [int(x) for x in tvec(nu, Sinv, h)]
    Mfp = [Fr(x) for x in tvec(card.M, Sinv, h)]
    sig = [int(x) for x in tvec(sigma, Sinv, h)] if sigma and \
        any(s != 1 for s in sigma) else None
    assert sig is None, "sign-twisted coni frame not implemented (A7 PIN-REQ)"
    certs, caps_a = _one_cert(V, nup, h)
    a = [Fr(0)]*(M+1); J1 = [Fr(0)]*(M+1); J2 = [Fr(0)]*(M+1)
    skip_jets = bool(os.environ.get("TERRIER_SKIP_JETS"))
    def visit(n, m):                      # STREAMING (no point list: 31G law)
        vals = [sum(v[i]*n[i] for i in range(h)) for v in V]
        c = coeff(facs, vals)
        a[m] += c
        if not skip_jets:
            nM = sum(Mfp[i]*n[i] for i in range(h))
            J1[m] += c*nM; J2[m] += c*nM*nM
    t0 = time.time()
    npts = cone_points(V, nup, M, certs, caps_a, h, visit=visit)
    t_all = time.time() - t0
    for x in a:
        assert x.denominator == 1, "non-integer graded coefficient"
    return a, J1, J2, dict(npoints=npts, t_enum=t_all, t_coef=0.0,
                           caps_alpha=[str(x) for x in caps_a],
                           Sinv=Sinv)


if __name__ == "__main__":
    from family import load_card
    from conipfv.coni_curve import coni_pfv_curve
    path, M = sys.argv[1], int(sys.argv[2])
    card = load_card(path)
    cur = coni_pfv_curve(card)
    assert cur["routing"] == "CONI-PFV"
    nu = [int(x) for x in cur["nu"]]
    qcf = [int(x) for x in cur["q_cf"]]
    # conifold ray must be OUTSIDE the support cone (else D5 resummation)
    rows = [v for (v, m, isn) in card.factors]
    assert any(sum(v[i]*qcf[i] for i in range(card.h)) < 0 for v in rows), \
        "q_cf ray inside support cone: graded series needs D5 resummation"
    t0 = time.time()
    a, J1, J2, info = graded_series(card, nu, M, card.sigma)
    dt = time.time() - t0
    nz = sum(1 for x in a if x)
    print(f"{card.name}: M={M} pts={info['npoints']} nonzero_a={nz} "
          f"wall={dt:.1f}s (enum {info['t_enum']:.1f} coef {info['t_coef']:.1f})")
    print("  head:", [str(x) for x in a[:12]])
    out = dict(card=card.name, M=M, nu=nu, q_cf=qcf, r=int(cur["r"]),
               lam=str(cur["lam"]), wall_s=dt,
               info={k: info[k] for k in ("npoints", "caps_alpha", "Sinv")},
               a=[str(x) for x in a], J1=[str(x) for x in J1],
               J2=[str(x) for x in J2])
    fn = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      f"{card.name}.series.json")
    json.dump(out, open(fn, "w"))
    print("wrote", fn)
