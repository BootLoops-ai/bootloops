#!/usr/bin/env python3
"""
pipe_vac.py — vacuum layer (jet ideal + certified vacuum), parametrized
over (PF ideal, module basis, curve, flux, conventions).

EXACT part (sympy over Q; generic normal-form closure, no hand-ordered
elimination):
  * parse the effective PF ideal (normal-ordered sympy strings from the card);
  * generic module closure: singles cascade + coupled Gaussian solve over
    Q(z_1..z_h) for the boundary monomials' normal forms; leftover relations
    must reduce to 0 (consistency gate, generalizes th2*L2 -> 0);
  * connection matrices A_a (rank x rank, rational in z) — gate vs reference;
  * curve matrix M(s) (Theta = sum nu_a th_a, z_a = sigma_a s^{nu_a}),
    det M(s*) != 0 exact — gate vs reference;
  * generalized r3: first Q(s)-dependency of Theta^k[1] at k = r must
    reproduce the restrict_op operator EXACTLY (independent D-module route).
CERTIFIED part (arb/acb):
  * V0 = transverse jets at s* from the certified companion (M(s*)^T solve);
  * off-curve leg dV/dt = B(t)V with Groenwall+Cauchy tail (h_a may be balls);
  * F-terms + Jacobian by forward AD over acb (2h+2 real coordinates);
  * FD-vs-AD Jacobian gate;  Krawczyk interval-Newton (existence+uniqueness).
Resource caps: single thread; ulimit -v 32505856, nice 5.
"""
from fractions import Fraction as Fr
from math import comb, factorial
import json, os, sys, time

import pipe_lib as PL
import pipe_transport as PT

# ----------------------------------------------------- exact layer (sympy)
def parse_ideal(card):
    import sympy as sp
    h = card.h
    ths = sp.symbols(" ".join(f"th{i+1}" for i in range(h)))
    zs = sp.symbols(" ".join(f"z{i+1}" for i in range(h)))
    if h == 1:
        ths, zs = (ths,), (zs,)
    loc = {f"th{i+1}": ths[i] for i in range(h)}
    loc.update({f"z{i+1}": zs[i] for i in range(h)})
    ops = []
    for s_expr in card.ideal:
        e = sp.expand(sp.sympify(s_expr, locals=loc))
        p = sp.Poly(e, *ths)
        op = {}
        for mono, c in zip(p.monoms(), p.coeffs()):
            op[tuple(mono)] = sp.cancel(op.get(tuple(mono), 0) + c)
        ops.append({k: v for k, v in op.items() if v != 0})
    return ops, list(ths), list(zs)


def opmul_theta(op, a, zs):
    """left-multiply sum c_beta(z) th^beta by th_a."""
    import sympy as sp
    out = {}
    for mono, c in op.items():
        key = tuple(mono[i] + (1 if i == a else 0) for i in range(len(mono)))
        out[key] = sp.cancel(out.get(key, 0) + c)
        dc = sp.cancel(zs[a] * sp.diff(c, zs[a]))
        if dc != 0:
            out[mono] = sp.cancel(out.get(mono, 0) + dc)
    return {k: v for k, v in out.items() if v != 0}


def reduce_rel(rel, NF, B):
    import sympy as sp
    vec = {b: sp.Integer(0) for b in B}
    unk = {}
    for m, c in rel.items():
        if m in NF:
            for b, w in NF[m].items():
                vec[b] = sp.cancel(vec[b] + c * w)
        else:
            unk[m] = sp.cancel(unk.get(m, 0) + c)
    return vec, {k: v for k, v in unk.items() if v != 0}

def solve_module(gens, B, h, zs):
    """generic normal-form closure: NF for all boundary monomials b+e_a.
    Returns (NF, leftovers) — leftovers = generated relations not used in the
    solve (consistency material)."""
    import sympy as sp
    NF = {b: {b: sp.Integer(1)} for b in B}
    needed = set()
    for b in B:
        for a in range(h):
            m = tuple(b[i] + (1 if i == a else 0) for i in range(h))
            if m not in NF:
                needed.add(m)
    dmax = max(sum(m) for m in needed)
    rounds = [list(gens)]
    while len(rounds) <= dmax:
        nxt = []
        for rel in rounds[-1]:
            for a in range(h):
                nxt.append(opmul_theta(rel, a, zs))
        rounds.append(nxt)
    allrels = [r for rnd in rounds for r in rnd]
    used = [False] * len(allrels)
    progress = True
    while progress:                               # singles cascade
        progress = False
        for i, rel in enumerate(allrels):
            if used[i]:
                continue
            vec, unk = reduce_rel(rel, NF, B)
            if len(unk) == 1:
                (tgt, alpha), = unk.items()
                if sum(tgt) <= dmax and tgt not in NF:
                    NF[tgt] = {b: sp.cancel(-v / alpha)
                               for b, v in vec.items() if v != 0}
                    used[i] = True
                    progress = True
    if not needed <= set(NF):                     # coupled Gaussian solve
        sel = []
        U = set()
        for i, rel in enumerate(allrels):
            if used[i]:
                continue
            vec, unk = reduce_rel(rel, NF, B)
            if unk and all(sum(m) <= dmax for m in unk):
                sel.append((i, vec, unk))
                U |= set(unk)
        U = sorted(U)
        assert len(sel) >= len(U), \
            f"cannot close module: {len(sel)} rels for {len(U)} unknowns"
        Ar = [[unk.get(m, sp.Integer(0)) for m in U] for i, vec, unk in sel]
        Rr = [[vec[b] for b in B] for i, vec, unk in sel]
        rl = 0
        piv = []
        for c in range(len(U)):
            pr = next((r for r in range(rl, len(Ar))
                       if sp.cancel(Ar[r][c]) != 0), None)
            assert pr is not None, "coupled system rank-deficient"
            Ar[rl], Ar[pr] = Ar[pr], Ar[rl]
            Rr[rl], Rr[pr] = Rr[pr], Rr[rl]
            pv = Ar[rl][c]
            Ar[rl] = [sp.cancel(x / pv) for x in Ar[rl]]
            Rr[rl] = [sp.cancel(x / pv) for x in Rr[rl]]
            for r in range(len(Ar)):
                if r != rl and Ar[r][c] != 0:
                    f = Ar[r][c]
                    Ar[r] = [sp.cancel(Ar[r][t] - f * Ar[rl][t])
                             for t in range(len(U))]
                    Rr[r] = [sp.cancel(Rr[r][t] - f * Rr[rl][t])
                             for t in range(len(B))]
            piv.append(c)
            rl += 1
        for i, c in enumerate(piv):               # x_c = -R_i
            NF[U[c]] = {B[j]: sp.cancel(-Rr[i][j]) for j in range(len(B))
                        if sp.cancel(Rr[i][j]) != 0}
        for r in range(rl, len(Ar)):              # leftover rows: consistency
            assert all(sp.cancel(x) == 0 for x in Ar[r])
            assert all(sp.cancel(x) == 0 for x in Rr[r]), \
                "over-determined module system INCONSISTENT"
        for i, vec, unk in sel:
            used[i] = True
    assert needed <= set(NF), f"unclosed: {needed - set(NF)}"
    leftovers = [allrels[i] for i in range(len(allrels)) if not used[i]]
    return NF, leftovers

def nf_theta_times(NF, mono_nf, a, zs):
    """theta_a . (sum c_b b) reduced onto the basis (needs NF of b+e_a)."""
    import sympy as sp
    out = {}
    for b, c in mono_nf.items():
        key = tuple(b[i] + (1 if i == a else 0) for i in range(len(b)))
        tgt = NF[key] if key in NF else {key: sp.Integer(1)}
        for bb, w in tgt.items():
            out[bb] = sp.cancel(out.get(bb, 0) + c * w)
        dc = sp.cancel(zs[a] * sp.diff(c, zs[a]))
        if dc != 0:
            for bb, w in (NF[b] if b in NF else {b: sp.Integer(1)}).items():
                out[bb] = sp.cancel(out.get(bb, 0) + dc * w)
    return {k: v for k, v in out.items() if v != 0}


def nf_extend(NF, m, h, zs):
    """recursively derive NF(m) = theta_a NF(m - e_a) for higher monomials."""
    if m in NF:
        return
    for a in range(h):
        if m[a] > 0:
            prev = tuple(m[i] - (1 if i == a else 0) for i in range(h))
            nf_extend(NF, prev, h, zs)
            sub = {b: NF[prev][b] for b in NF[prev]}
            # ensure all b+e_a of the support are known
            for b in list(sub):
                key = tuple(b[i] + (1 if i == a else 0) for i in range(h))
                if key not in NF and key != m:
                    nf_extend(NF, key, h, zs)
            NF[m] = nf_theta_times(NF, sub, a, zs)
            return
    raise AssertionError(f"cannot extend NF to {m}")


def consistency_gate(NF, leftovers, B, h, zs, nmax=4):
    """the first nmax leftover relations must reduce to 0 in the module."""
    import sympy as sp
    nchk = 0
    for rel in leftovers:
        if nchk >= nmax:
            break
        for m in rel:
            if m not in NF:
                nf_extend(NF, m, h, zs)
        vec, unk = reduce_rel(rel, NF, B)
        assert not unk
        assert all(sp.cancel(v) == 0 for v in vec.values()), \
            "leftover relation does NOT reduce to 0 — module inconsistent"
        nchk += 1
    return nchk


def connection_matrices(NF, B, h, zs):
    """(A_a)_{lj}: th_a b_j = sum_l (A_a)_{lj} b_l."""
    import sympy as sp
    A = []
    for a in range(h):
        Aa = [[sp.Integer(0)] * len(B) for _ in range(len(B))]
        for j, bj in enumerate(B):
            key = tuple(bj[i] + (1 if i == a else 0) for i in range(h))
            nf = NF[key] if key in NF else {key: sp.Integer(1)}
            for bb, c in nf.items():
                assert bb in B, (a, bj, bb)
                Aa[B.index(bb)][j] = sp.cancel(c)
        A.append(Aa)
    return A

def curve_matrix(NF, B, h, zs, nu, sigma, r_ord, s_star):
    """M[j,k] = coeff of b_j in Theta^k[1], Theta = sum nu_a th_a, on the
    curve z_a = sigma_a s^{nu_a}.  Also returns the k=r_ord column for the
    generalized r3 (Route-B) gate and det M(s*) != 0 (exact)."""
    import sympy as sp
    s = sp.symbols("s")
    sub = {zs[a]: sigma[a] * s ** nu[a] for a in range(h)}

    def theta_action(vecdict):
        out = {b: sp.Integer(0) for b in B}
        for b, c in vecdict.items():
            dc = sp.cancel(s * sp.diff(c, s))
            if dc != 0:
                out[b] = sp.cancel(out[b] + dc)
            for a in range(h):
                key = tuple(b[i] + (1 if i == a else 0) for i in range(h))
                tgt = NF[key] if key in NF else {key: sp.Integer(1)}
                for bb, cc in tgt.items():
                    assert bb in B
                    out[bb] = sp.cancel(out[bb] + nu[a] * c * cc.subs(sub))
        return {k: sp.cancel(v) for k, v in out.items() if v != 0}

    V = [{B[0]: sp.Integer(1)}]
    assert B[0] == (0,) * h
    for k in range(r_ord):
        V.append(theta_action(V[-1]))
    M = sp.Matrix([[sp.together(V[k].get(b, 0)) for k in range(r_ord)]
                   for b in B])
    detM = sp.cancel(sp.together(M.det()))
    assert detM != 0, "M singular over Q(s)"
    sstar = sp.Rational(s_star.numerator, s_star.denominator)
    dnum, dden = sp.fraction(detM)
    assert dnum.subs(s, sstar) != 0, "det M vanishes at s*"
    return M, detM, V, s


def routeB_gate(V, B, r_ord, s):
    """first Q(s)-dependency of Theta^k[1] must land at k = r_ord; returns
    the integer-normalized operator coefficients b_k(s) (ascending lists)."""
    import sympy as sp
    Mfull = sp.Matrix([[sp.together(V[k].get(b, 0)) for k in range(r_ord + 1)]
                       for b in B])
    for k in range(1, r_ord):
        assert Mfull[:, :k + 1].rank() == k + 1, f"early dependency at k={k}"
    ns = Mfull.nullspace()
    assert len(ns) == 1, "dependency space not 1-dimensional at k=r"
    coefv = [sp.cancel(x) for x in ns[0]]
    lcm = sp.Integer(1)
    for x in coefv:
        lcm = sp.lcm(lcm, sp.denom(sp.together(x)))
    bpoly = [sp.expand(sp.cancel(x * lcm)) for x in coefv]
    g = sp.Integer(0)
    mmin = None
    for x in bpoly:
        if x == 0:
            continue
        P = sp.Poly(x, s)
        for c in P.all_coeffs():
            g = sp.gcd(g, c)
        low = min(mm[0] for mm in P.monoms())
        mmin = low if mmin is None else min(mmin, low)
    bpoly = [sp.expand(sp.cancel(x / (g * s ** mmin))) for x in bpoly]
    out = []
    for x in bpoly:
        if x == 0:
            out.append([Fr(0)])
            continue
        P = sp.Poly(x, s)
        cs = [Fr(0)] * (P.degree() + 1)
        for mm, c in zip(P.monoms(), P.coeffs()):
            cs[mm[0]] = Fr(int(sp.Integer(c)))
        out.append(cs)
    return out

def ser_rath(expr, zs):
    """rational function in z -> {num:{'d1,..,dh': str}, den:{...}}."""
    import sympy as sp
    n, d = sp.fraction(sp.cancel(sp.together(expr)))
    out = {}
    for pol, tag in ((sp.expand(n), "num"), (sp.expand(d), "den")):
        P = sp.Poly(pol, *zs)
        out[tag] = {",".join(map(str, m)): str(c)
                    for m, c in zip(P.monoms(), P.coeffs())}
    return out


def ser_rat1(expr, s):
    import sympy as sp
    n, d = sp.fraction(sp.cancel(sp.together(expr)))
    out = {}
    for pol, tag in ((sp.expand(n), "num"), (sp.expand(d), "den")):
        P = sp.Poly(pol, s)
        out[tag] = {str(m[0]): str(c) for m, c in zip(P.monoms(), P.coeffs())}
    return out


def jet_layer(card, curve, op, outdir):
    """full exact vacuum-layer data build + banking; returns file paths."""
    import sympy as sp
    t0 = time.time()
    h = card.h
    gens, ths, zs = parse_ideal(card)
    B = card.module_basis
    assert B is not None and B[0] == (0,) * h
    r_ord = max(int(k.split(",")[1]) for k in op["theta_form"])
    assert len(B) == r_ord, "module basis size != operator order (order drop?)"
    NF, leftovers = solve_module(gens, B, h, zs)
    print(f"[jets] normal forms closed ({time.time()-t0:.0f}s)", flush=True)
    nchk = consistency_gate(NF, leftovers, B, h, zs)
    print(f"[jets] consistency: {nchk} leftover relations reduce to 0  [PASS]",
          flush=True)
    A = connection_matrices(NF, B, h, zs)
    nu = tuple(curve["nu"])
    sigma = card.sigma
    M, detM, V, s = curve_matrix(NF, B, h, zs, nu, sigma, r_ord, card.s_star)
    print(f"[jets] connection + curve matrix done ({time.time()-t0:.0f}s); "
          f"det M(s*) != 0 exact", flush=True)
    bB = routeB_gate(V, B, r_ord, s)
    # generalized r3: cross-products of route-B polys vs the operator
    PA = {k: [Fr(c) for c in PL.PK[k]] for k in range(r_ord + 1)}
    def polymul(a, b):
        out = [Fr(0)] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b):
                    if y:
                        out[i + j] += x * y
        return out
    def poltrim(a):
        while a and a[-1] == 0:
            a = a[:-1]
        return a
    for k in range(r_ord + 1):
        for j in range(k):
            lhs = poltrim(polymul(PA[k], bB[j]))
            rhs = poltrim(polymul(PA[j], bB[k]))
            assert lhs == rhs, f"r3 mismatch at (k,j)=({k},{j})"
    print(f"[jets] r3 GATE PASS — D-module route == annihilator route exactly "
          f"(all {r_ord*(r_ord+1)//2} cross-products)", flush=True)
    conn = {"B": [list(b) for b in B],
            "A": [[[ser_rath(A[a][l][j], zs) for j in range(len(B))]
                   for l in range(len(B))] for a in range(h)],
            "note": "th_a (b_j Pi) = sum_l (A_a)_{lj}(z) (b_l Pi)"}
    with open(os.path.join(outdir, "connection_z.json"), "w") as f:
        json.dump(conn, f)
    curveT = {"B": [list(b) for b in B], "s_star": str(card.s_star),
              "M": [[ser_rat1(M[i, j], s) for j in range(len(B))]
                    for i in range(len(B))],
              "detM": ser_rat1(detM, s)}
    with open(os.path.join(outdir, "curve_T.json"), "w") as f:
        json.dump(curveT, f)
    print("[jets] wrote connection matrices + curve matrix to outdir", flush=True)
    return conn, curveT

# ===================================================== certified layer
from flint import arb, acb, fmpq, ctx

CARDG = None; RK = None; NAD = None; ACONN = None; M_AT_S = None
ZS = None; S_STAR = None; TAU_PIN = None; BASIS = None; IDXE = None
STATE = {}


def fr2arb(x):
    return arb(fmpq(x.numerator, x.denominator))


def fr2acb(x):
    return acb(fr2arb(x))


def parse_ph(d):
    return {tuple(map(int, k.split(","))): Fr(v) for k, v in d.items()}


def parse_p1(d):
    return {int(k): Fr(v) for k, v in d.items()}


def configure_cert(card, curve, conn, curveT):
    """call after PL.configure() + PT.configure()."""
    global CARDG, RK, NAD, ACONN, M_AT_S, ZS, S_STAR, TAU_PIN, BASIS, IDXE
    CARDG = card
    h = card.h
    BASIS = [tuple(b) for b in conn["B"]]
    RK = len(BASIS)
    assert RK == PT.R, "rank != operator order"
    NAD = 2 * h + 2
    ACONN = [[[ (parse_ph(conn["A"][a][l][j]["num"]),
                 parse_ph(conn["A"][a][l][j]["den"]))
               for j in range(RK)] for l in range(RK)] for a in range(h)]
    S_STAR = card.s_star
    TAU_PIN = card.tau_pin
    assert curveT["s_star"] == str(S_STAR)
    def rat1_at(nd, x):
        num = sum(c * x ** d for d, c in parse_p1(nd["num"]).items()) \
            if nd["num"] else Fr(0)
        den = sum(c * x ** d for d, c in parse_p1(nd["den"]).items())
        return num / den
    global M_AT_S
    M_AT_S = [[rat1_at(curveT["M"][j][k], S_STAR) for k in range(RK)]
              for j in range(RK)]
    detm = rat1_at(curveT["detM"], S_STAR)
    assert detm != 0
    nu = tuple(curve["nu"])
    ZS = [Fr(card.sigma[a]) * S_STAR ** nu[a] for a in range(h)]
    IDXE = [BASIS.index(tuple(1 if i == a else 0 for i in range(h)))
            for a in range(h)]


def V0_at_star(comp, prec):
    """V0[j][i] = (b_j Pi_i)(z*) (v^3-scaled) from the certified companion."""
    ctx.prec = prec
    D = []
    for k in range(RK):
        row = []
        for i in range(PL.NPER):
            v = acb(0)
            for j in range(k + 1):
                if PT.S2[k][j]:
                    v += fr2acb(Fr(PT.S2[k][j] * factorial(j))
                                * S_STAR ** j) * comp[i][j]
            row.append(v)
        D.append(row)
    A = [[fr2acb(M_AT_S[j][k]) for j in range(RK)] for k in range(RK)]  # M^T
    X = [row[:] for row in D]
    for c in range(RK):
        p = max(range(c, RK), key=lambda r: float(abs(A[r][c]).mid()))
        A[c], A[p] = A[p], A[c]; X[c], X[p] = X[p], X[c]
        for r in range(RK):
            if r != c:
                f = A[r][c] / A[c][c]
                A[r] = [A[r][t] - f * A[c][t] for t in range(RK)]
                X[r] = [X[r][t] - f * X[c][t] for t in range(PL.NPER)]
    V0 = [[X[c][i] / A[c][c] for i in range(PL.NPER)] for c in range(RK)]
    i0 = BASIS.index((0,) * CARDG.h)
    for i in range(PL.NPER):
        assert (V0[i0][i] - comp[i][0]).contains(acb(0)), "V0 row-0 mismatch"
    return V0

# ---------------- certified off-curve leg (tail proofs: cert_vac docstring)
def ser_inv(p, N):
    out = [acb(0)] * (N + 1); out[0] = 1 / p[0]
    for m in range(1, N + 1):
        s = acb(0)
        for k in range(1, min(m, len(p) - 1) + 1):
            s += p[k] * out[m - k]
        out[m] = -s * out[0]
    return out


def ser_mul(a, b, N):
    out = [acb(0)] * (N + 1)
    for i, x in enumerate(a[:N + 1]):
        for j in range(0, min(N - i, len(b) - 1) + 1):
            out[i + j] += x * b[j]
    return out


def polyh_tser(P, zc, hv):
    """P(z_c + t h) as acb t-poly list; zc, hv = lists of acb."""
    h = len(zc)
    dm = [max((k[a] for k in P), default=0) for a in range(h)]
    pw = []
    for a in range(h):
        tab = {0: [acb(1)]}
        for d in range(1, dm[a] + 1):
            tab[d] = [fr2acb(Fr(comb(d, k))) * zc[a] ** (d - k) * hv[a] ** k
                      for k in range(d + 1)]
        pw.append(tab)
    out = [acb(0)] * (sum(dm) + 1)
    for mono, c in P.items():
        pr = [acb(1)]
        for a in range(h):
            if mono[a]:
                pr = ser_mul(pr, pw[a][mono[a]], sum(mono))
        ca = fr2acb(c)
        for i, x in enumerate(pr):
            out[i] += ca * x
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def polyh_at(P, zz):
    v = acb(0)
    for mono, c in P.items():
        t = fr2acb(c)
        for a in range(len(zz)):
            if mono[a]:
                t *= zz[a] ** mono[a]
        v += t
    return v


def abs_lb(x):
    a = abs(x)
    lb = a - arb(0, a.rad())
    return lb if lb > 0 else arb(0)


def leg_transport(V0, hvec, prec, digits):
    """V(1) enclosure for dV/dt = B(t)V, z(t) = z* + t h (h may be balls)."""
    ctx.prec = prec
    h = CARDG.h
    zc = [fr2acb(ZS[a]) for a in range(h)]
    if all(x == 0 for x in hvec):
        return [row[:] for row in V0]
    R0 = 2000
    import math
    while True:
        tb = acb(arb(0, R0), arb(0, R0))
        zb = [zc[a] + tb * hvec[a] for a in range(h)]
        ok = all(abs_lb(z) > 0 for z in zb)
        alpha = arb(0)
        if ok:
            rows = [arb(0)] * RK
            for a in range(h):
                for l in range(RK):
                    for j in range(RK):
                        num, den = ACONN[a][l][j]
                        if not num:
                            continue
                        dl = abs_lb(polyh_at(den, zb)) * abs_lb(zb[a])
                        if not dl > 0:
                            ok = False; break
                        rows[j] += abs(hvec[a]) * \
                            abs(polyh_at(num, zb)).upper() / dl
                    if not ok:
                        break
                if not ok:
                    break
            if ok:
                alpha = arb(0)
                for r in rows:
                    alpha = alpha.max(r.upper())
        if ok and float((alpha * R0).upper()) <= 6.0:
            break
        R0 //= 4
        assert R0 >= 4, "no valid analyticity radius for leg"
    N = int(digits * math.log(10) / math.log(R0)) + 15
    # B(t) Taylor coefficients: RK x RK lists;  C[j][l] <- (A_a)_{lj} h_a/z_a
    Bser = [[[acb(0)] * (N + 1) for _ in range(RK)] for _ in range(RK)]
    for a in range(h):
        ha = hvec[a]
        if ha == 0:
            continue
        zt = [zc[a], ha]
        for l in range(RK):
            for j in range(RK):
                num, den = ACONN[a][l][j]
                if not num:
                    continue
                nt = polyh_tser(num, zc, hvec)
                dt = polyh_tser(den, zc, hvec)
                assert len(dt) + 1 <= N, "leg N too small for denominator"
                dt = ser_mul(dt, zt, len(dt) + 1)
                ent = ser_mul(nt, ser_inv(dt, N), N)
                dst = Bser[j][l]
                for m in range(N + 1):
                    dst[m] += ha * ent[m]
    Vm = [row[:] for row in V0]
    S = [row[:] for row in V0]
    hist = [Vm]
    for m in range(N):
        nxt = [[acb(0)] * PL.NPER for _ in range(RK)]
        for k in range(m + 1):
            Vk = hist[m - k]
            for l in range(RK):
                for j in range(RK):
                    b = Bser[l][j][k]
                    if b != 0:
                        for i in range(PL.NPER):
                            nxt[l][i] += b * Vk[j][i]
        inv = arb(1) / (m + 1)
        Vm = [[x * inv for x in row] for row in nxt]
        hist.append(Vm)
        for l in range(RK):
            for i in range(PL.NPER):
                S[l][i] += Vm[l][i]
    ee = (alpha * R0).exp()
    geo = (arb(1) / R0) ** (N + 1) / (1 - arb(1) / R0)
    for i in range(PL.NPER):
        Mc = arb(0)
        for l in range(RK):
            Mc = Mc.max(abs(V0[l][i]).upper())
        tb = (Mc * ee * geo).upper()
        pm = arb(0, tb)
        for l in range(RK):
            S[l][i] += acb(pm, pm)
    return S


# ---------------- forward AD over acb (NAD real coordinates)
class AD:
    __slots__ = ("v", "p")
    def __init__(self, v, p=None):
        self.v = v if isinstance(v, acb) else acb(v)
        self.p = p if p is not None else [acb(0)] * NAD
    def __add__(s, o):
        o = o if isinstance(o, AD) else AD(o)
        return AD(s.v + o.v, [s.p[i] + o.p[i] for i in range(NAD)])
    __radd__ = __add__
    def __sub__(s, o):
        o = o if isinstance(o, AD) else AD(o)
        return AD(s.v - o.v, [s.p[i] - o.p[i] for i in range(NAD)])
    def __rsub__(s, o):
        return (AD(o) if not isinstance(o, AD) else o).__sub__(s)
    def __mul__(s, o):
        o = o if isinstance(o, AD) else AD(o)
        return AD(s.v * o.v, [s.p[i] * o.v + s.v * o.p[i] for i in range(NAD)])
    __rmul__ = __mul__
    def __truediv__(s, o):
        o = o if isinstance(o, AD) else AD(o)
        iv = 1 / o.v
        return AD(s.v * iv, [(s.p[i] * o.v - s.v * o.p[i]) * iv * iv
                             for i in range(NAD)])
    def __rtruediv__(s, o):
        return (AD(o) if not isinstance(o, AD) else o).__truediv__(s)
    def __neg__(s):
        return AD(-s.v, [-x for x in s.p])
    def conj(s):
        return AD(s.v.conjugate(), [x.conjugate() for x in s.p])

def sympl(a, b):
    n = CARDG.h + 1
    return sum(a[i] * b[n + i] for i in range(n)) \
        - sum(a[n + i] * b[i] for i in range(n))


def eval_state(x, prec, digits):
    """x = NAD arb balls (Re w_1, Im w_1, ..., Re w_h, Im w_h, Re tau, Im tau).
    Returns dict: E (h+1 AD nodes: E_tau, E_a) + extras at the point."""
    ctx.prec = prec
    h = CARDG.h
    I = acb(0, 1)
    wv = [acb(x[2*a], x[2*a + 1]) for a in range(h)]
    tauv = acb(x[NAD - 2], x[NAD - 1])
    w = []
    for a in range(h):
        p = [acb(0)] * NAD
        p[2*a], p[2*a + 1] = acb(1), I
        w.append(AD(wv[a], p))
    ptau = [acb(0)] * NAD
    ptau[NAD - 2], ptau[NAD - 1] = acb(1), I
    tau = AD(tauv, ptau)
    hvec = [fr2acb(ZS[a]) * wv[a] for a in range(h)]
    V = leg_transport(STATE["V0"], hvec, prec, digits)
    v3 = acb(0, 2 * arb.pi()) ** 3
    zz = [fr2acb(ZS[a]) * (1 + wv[a]) for a in range(h)]
    zn = [1 + wv[a] for a in range(h)]
    i0 = BASIS.index((0,) * h)
    Pi, Th = [], [[] for _ in range(h)]
    for i in range(PL.NPER):
        pv = V[i0][i] / v3
        t = [V[IDXE[a]][i] / v3 for a in range(h)]
        tt = [[None] * h for _ in range(h)]
        for a in range(h):
            for b in range(a, h):
                val = acb(0)
                for l in range(RK):
                    num, den = ACONN[a][l][IDXE[b]]
                    if num:
                        val += polyh_at(num, zz) / polyh_at(den, zz) * V[l][i]
                tt[a][b] = tt[b][a] = val / v3
        pp = [acb(0)] * NAD
        for a in range(h):
            pp[2*a] = t[a] / zn[a]
            pp[2*a + 1] = I * t[a] / zn[a]
        Pi.append(AD(pv, pp))
        for a in range(h):
            pa = [acb(0)] * NAD
            for b in range(h):
                pa[2*b] = tt[b][a] / zn[b]
                pa[2*b + 1] = I * tt[b][a] / zn[b]
            Th[a].append(AD(t[a], pa))
    one = AD(acb(1))
    dPi = [[Th[a][i] / (one + w[a]) for i in range(PL.NPER)] for a in range(h)]
    sq2pi = AD((2 / arb.pi()).sqrt())
    FmtH = [AD(acb(PL.FFLUX[i])) - tau * acb(PL.HFLUX[i])
            for i in range(PL.NPER)]
    Hc = [AD(acb(hh)) for hh in PL.HFLUX]
    W = sq2pi * sympl(FmtH, Pi)
    Bal = sympl(Hc, Pi)
    PiC = [p.conj() for p in Pi]
    Ncs = AD(acb(0, -1)) * sympl(PiC, Pi)
    Ktau = AD(acb(-1)) / (tau - tau.conj())
    Etau = sq2pi * (-Bal) + Ktau * W
    E = [Etau]
    for a in range(h):
        Na = AD(acb(0, -1)) * sympl(PiC, dPi[a])
        Ka = -(Na / Ncs)
        E.append(sq2pi * sympl(FmtH, dPi[a]) + Ka * W)
    emK = Ncs.v.real * 2 * x[NAD - 1]
    X0 = Pi[h + 1].v
    return {"E": E, "W": W.v, "emKcs": Ncs.v, "emK": emK,
            "W0": abs(W.v) / emK.sqrt(),
            "U": [Pi[h + 2 + a].v / X0 for a in range(h)],
            "B": Bal.v, "Pi": [p.v for p in Pi], "tau": tauv}

# ---------------- Krawczyk interval-Newton (+ FD-Jacobian gate)
def Fvec(E):
    out = []
    for e in range(len(E)):
        out += [E[e].v.real, E[e].v.imag]
    return out


def jac_real(E):
    J = [[None] * NAD for _ in range(NAD)]
    for e in range(len(E)):
        for j in range(NAD):
            J[2 * e][j] = E[e].p[j].real
            J[2 * e + 1][j] = E[e].p[j].imag
    return J


def inv_mid(J):
    n = NAD
    A = [[arb(J[i][j].mid()) for j in range(n)]
         + [arb(1 if k == i else 0) for k in range(n)] for i in range(n)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: float(abs(A[r][c]).mid()))
        A[c], A[p] = A[p], A[c]
        pv = A[c][c]
        A[c] = [x / pv for x in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]
                A[r] = [A[r][t] - f * A[c][t] for t in range(2 * n)]
    return [[arb(A[i][n + j].mid()) for j in range(n)] for i in range(n)]


def newton_center(x0, prec, digits, iters=30):
    x = [arb(v) for v in x0]
    res = None
    for it in range(iters):
        st = eval_state(x, prec, digits)
        F = Fvec(st["E"])
        res = max(float(abs(f).upper()) for f in F)
        J = jac_real(st["E"])
        Y = inv_mid(J)
        dx = [sum(Y[i][j] * F[j] for j in range(NAD)) for i in range(NAD)]
        x = [arb((x[i] - dx[i]).mid()) for i in range(NAD)]
        if res < 1e-300:
            break
    return x, res


def krawczyk(y, r, prec, digits):
    stF = eval_state(y, prec, digits)
    F = Fvec(stF["E"])
    xb = [y[i] + arb(0, 1) * r[i] for i in range(NAD)]
    stJ = eval_state(xb, prec, digits)
    J = jac_real(stJ["E"])
    Y = inv_mid(J)
    K = []
    for i in range(NAD):
        v = y[i] - sum(Y[i][j] * F[j] for j in range(NAD))
        for j in range(NAD):
            c = (arb(1) if i == j else arb(0)) - \
                sum(Y[i][t] * J[t][j] for t in range(NAD))
            v += c * (arb(0, 1) * r[j])
        K.append(v)
    contained = all((abs(arb(K[i].mid()) - y[i]) + arb(0, K[i].rad())).upper()
                    < r[i] for i in range(NAD))
    return K, contained, stF


def fd_jacobian_gate(y, prec, digits, hstep=None, bar=1e-6):
    """AD Jacobian vs central finite differences at a point (development
    gate: catches e.g. a transposed connection matrix)."""
    hstep = hstep or arb(10) ** (-(digits // 3))
    st = eval_state(y, prec, digits)
    J = jac_real(st["E"])
    worst = 0.0
    for j in range(NAD):
        yp = list(y); yp[j] = y[j] + hstep
        ym = list(y); ym[j] = y[j] - hstep
        Fp = Fvec(eval_state(yp, prec, digits)["E"])
        Fm = Fvec(eval_state(ym, prec, digits)["E"])
        for i in range(NAD):
            fd = (Fp[i] - Fm[i]) / (2 * hstep)
            ad = J[i][j]
            sc = max(float(abs(ad).mid()), float(abs(fd).mid()), 1e-30)
            rel = float(abs(fd - ad).mid()) / sc
            if sc > 1e-12:                     # entries with signal
                worst = max(worst, rel)
    assert worst < bar, f"FD-vs-AD Jacobian mismatch: {worst:.3e}"
    return worst
