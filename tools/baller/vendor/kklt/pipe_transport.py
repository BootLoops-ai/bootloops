# baller transport engine — certified numeric (arb/acb ball) layer.
#!/usr/bin/env python3
"""
pipe_transport.py — CERTIFIED numeric layer (arb/acb balls), parametrized
over the operator (order r, degree D) and the family card; generalized
from a fixed-instance pilot.  Every inequality is the pilot's, with 6 -> r
and 63 -> D parametric; for the DKMM card the numbers must reproduce the
pilot's banked values.

MUM leg (s = 0 -> s0, rational, single certified summation):
  tail by recurrence-majorant induction |a_m| <= C rho^{-m}; envelopes
  (exact Fractions, monotone decreasing in m, evaluated at m0 = N+1):
    phi(m)  <= [sum_d (|P_{r,d}| + sum_{k<r}|P_{k,d}|/m) rho^d]
               / [|lead| prod_i (1 - e_i/m)^{mult_i}]
    psi_i(m)<= [sum_d rho^d sum_k C(k,i)|P_{k,d}| m^{k-i-r}] / [same]
  using ind(m) = lead prod (m - e_i)^{mult_i} = lead m^r prod (1-e_i/m)^mult
  (exact integer-root factorization asserted in pipe_lib.ind_m).

Ordinary-point legs: D-form Q_j = s^j sum_k S2(k,j) P_k;  local relation
  sum_delta A_delta(m) c_{m-delta} = 0, A_delta(m) = sum_j Q_{j,delta-r+j}
  (m-delta)_j, A_0 = Q_{r,0}(m)_r;  majorant phi_loc(m) <= [sum rho^delta
  sum_j |Q_j,.| m^{j-r}] / [|Q_{r,0}| (1-(r-1)/m)^r];  window-C ball union.
Tails of derivative sums: descriptor algebra as pilot (u_e = ((N+2)/(N+1))^e).
Apparent singularities are never crossed; they enter only through rho.
"""
from fractions import Fraction as Fr
from math import comb, factorial
import json, os
from flint import arb, acb, fmpq, ctx

import pipe_lib as PL

# ------------------------------------------------ configured globals
R = None; DEG = None; S2 = None; QPOLY = None; ABSR = None; ABSLOW = None


def configure():
    """call after PL.configure()."""
    global R, DEG, S2, QPOLY, ABSR, ABSLOW
    R = PL.R_ORD
    DEG = PL.DEG
    S2 = [[0] * (R + 1) for _ in range(R + 1)]
    S2[0][0] = 1
    for k in range(1, R + 1):
        for j in range(1, k + 1):
            S2[k][j] = j * S2[k - 1][j] + S2[k - 1][j - 1]
    QPOLY = {}
    for j in range(0, R + 1):
        base = [0] * (DEG + 1)
        for k in range(j, R + 1):
            if S2[k][j]:
                for d in range(DEG + 1):
                    base[d] += S2[k][j] * PL.PK[k][d]
        QPOLY[j] = [0] * j + base
    ABSR = [abs(PL.PK[R][d]) for d in range(DEG + 1)]
    ABSLOW = [sum(abs(PL.PK[k][d]) for k in range(R)) for d in range(DEG + 1)]


def fr2arb(x):
    return arb(fmpq(x.numerator, x.denominator))


def fr2acb(z):
    if isinstance(z, tuple):
        return acb(fr2arb(z[0]), fr2arb(z[1]))
    return acb(fr2arb(z))


def cxq(z):
    return z if isinstance(z, tuple) else (Fr(z), Fr(0))


def cxq_sub(a, b):
    a, b = cxq(a), cxq(b)
    return (a[0] - b[0], a[1] - b[1])


def cxq_mul(a, b):
    a, b = cxq(a), cxq(b)
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def cxq_abs_ub(a):
    a = cxq(a)
    return abs(a[0]) + abs(a[1])


def cxq_abs_lb(a):
    a = cxq(a)
    return max(abs(a[0]), abs(a[1]))

# ---------------------------------------------------- MUM majorant (exact)
def ind_den(m0):
    """|lead| * prod (1 - e_i/m0)^mult  ==  |ind(m0)| / m0^r, exact."""
    den = Fr(abs(PL.IND_LEAD))
    for e, mu in PL.IND_FAC:
        if e:
            den *= (1 - Fr(e, m0)) ** mu
    return den


def phi_env(m0, rho):
    num = sum((Fr(ABSR[d]) + Fr(ABSLOW[d], m0)) * rho ** d
              for d in range(1, DEG + 1))
    return num / ind_den(m0)


def psi_env(m0, rho, i):
    num = Fr(0)
    for d in range(DEG + 1):
        s = Fr(0)
        for k in range(i, R + 1):
            c = comb(k, i) * abs(PL.PK[k][d])
            if c:
                s += Fr(c, m0 ** (R + i - k))     # m^{k-i-r}
        num += s * rho ** d
    return num / ind_den(m0)


def mum_layer_consts(F, layers, N, rho):
    phibar = phi_env(N + 1, rho)
    assert phibar < 1, f"MUM majorant fails: phibar={float(phibar):.3f}"
    Cs = {}
    for lp in layers:
        base = max(abs(F[lp][j]) * rho ** j for j in range(N - DEG, N + 1))
        src = Fr(0)
        for l2 in layers:
            if l2 > lp:
                i = l2 - lp
                src += PL.falling(l2, i) * Cs[l2] * psi_env(N + 1, rho, i)
        Cs[lp] = max(base, src / (1 - phibar))
    return Cs


def tail_desc(desc, N, x_ub, rho):
    y = fr2arb(x_ub) / fr2arb(rho)
    tot = arb(0)
    for e, K in desc:
        u = (fr2arb(Fr(N + 2, N + 1))) ** e
        r = u * y
        assert r < 1, "tail ratio >= 1"
        tot += fr2arb(K) * fr2arb(Fr(N + 1)) ** e * y ** (N + 1) / (1 - r)
    return tot.upper()


def pm_ball(ub):
    return arb(0, 1) * ub


# --------------------------------------------- MUM evaluation of a period
def mum_eval(T, N, s0, rho, nth=None):
    """theta^i y (i < nth) at rational s0 for an extended tower T."""
    nth = nth or R
    grades = sorted(set((P, Q) for (lp, P, Q) in T))
    ls0 = fr2arb(s0).log()
    v = acb(0, 2 * arb.pi())
    z3 = arb.zeta(arb(3))
    out = []
    diag = {"tail_worst": arb(0)}
    state = {}
    for g in grades:
        layers = sorted([lp for (lp, P, Q) in T if (P, Q) == g], reverse=True)
        F = {lp: T[(lp,) + g] for lp in layers}
        Cs = mum_layer_consts(F, layers, N, rho)
        state[g] = {lp: (F[lp][:], [(0, Cs[lp])]) for lp in layers}
    for i in range(nth):
        val = acb(0)
        for g in grades:
            P, Q = g
            gval = arb(0)
            for lp, (coef, desc) in state[g].items():
                s = arb(0)
                for m in range(N, -1, -1):
                    s = s * fr2arb(s0) + fr2arb(coef[m])
                tb = tail_desc(desc, N, s0, rho)
                diag["tail_worst"] = diag["tail_worst"].max(tb)
                gval += (s + pm_ball(tb)) * ls0 ** lp
            val += acb(gval) * v ** P * z3 ** Q
        out.append(val)
        for g in grades:
            st = state[g]
            keys = set(st) | {lp - 1 for lp in st if lp >= 1}
            new = {}
            for lp in keys:
                if lp in st:
                    coef, desc = st[lp]
                    nc = [Fr(m) * coef[m] for m in range(N + 1)]
                    nd = [(e + 1, K) for (e, K) in desc]
                else:
                    nc = [Fr(0)] * (N + 1)
                    nd = []
                if lp + 1 in st:
                    c2, d2 = st[lp + 1]
                    for m in range(N + 1):
                        if c2[m]:
                            nc[m] += (lp + 1) * c2[m]
                    nd += [(e, (lp + 1) * K) for (e, K) in d2]
                new[lp] = (nc, nd)
            state[g] = new
    return out, diag

# ------------------------------------------------- ordinary-point local leg
M0_LOC = 128


def local_shift(x1):
    """Q_{j,d} = t-coefficients of Q_j(x1+t); x1 Fr or (Fr,Fr).  Exact."""
    x1 = cxq(x1)
    Qs = {}
    for j in range(R + 1):
        q = QPOLY[j]
        deg = max(d for d in range(len(q)) if q[d]) if any(q) else 0
        pw = [(Fr(1), Fr(0))]
        for _ in range(deg):
            pw.append(cxq_mul(pw[-1], x1))
        out = [(Fr(0), Fr(0))] * (deg + 1)
        for e in range(deg + 1):
            if q[e]:
                for d in range(e + 1):
                    c = comb(e, d) * q[e]
                    re, im = pw[e - d]
                    out[d] = (out[d][0] + c * re, out[d][1] + c * im)
        Qs[j] = out
    return Qs


def phi_loc(Qs, rho, m0=M0_LOC):
    lb = cxq_abs_lb(Qs[R][0])
    assert lb > 0, "center is a root of Q_r (singular or s=0)"
    num = Fr(0)
    dmax = max(len(Qs[j]) for j in range(R + 1)) - 1
    for delta in range(1, dmax + R + 1):
        s = Fr(0)
        for j in range(R + 1):
            d = delta - R + j
            if 0 <= d < len(Qs[j]):
                ub = cxq_abs_ub(Qs[j][d])
                if ub:
                    s += ub * Fr(1, m0 ** (R - j))
        if s:
            num += s * rho ** delta
    return num / (lb * (1 - Fr(R - 1, m0)) ** R)


def pick_rho(Qs, hint=Fr(1, 50)):
    hi = hint
    while phi_loc(Qs, hi) > Fr(3, 4):
        hi = hi / 2
    lo, up = hi, hi * 2
    for _ in range(30):
        mid = (lo + up) / 2
        if phi_loc(Qs, mid) <= Fr(3, 4):
            lo = mid
        else:
            up = mid
    return lo


def local_leg(x1, seeds_list, h, N, prec, deriv_out=None):
    """Taylor-march all periods from center x1 to x1+h; certified."""
    deriv_out = deriv_out or R
    ctx.prec = prec
    DMAX = DEG + R
    Qs = local_shift(x1)
    rho = pick_rho(Qs)
    hub = cxq_abs_ub(h)
    assert hub < rho * Fr(9, 10), f"step too long: |h|={float(hub)}"
    assert N >= max(M0_LOC, DMAX + R)
    Qa = {j: [fr2acb(z) for z in Qs[j]] for j in range(R + 1)}
    ha = fr2acb(h)
    y = fr2arb(hub) / fr2arb(rho)
    diag = {"rho": float(rho), "phi": float(phi_loc(Qs, rho)), "N": N,
            "x1": (float(cxq(x1)[0]), float(cxq(x1)[1]))}
    cs = [list(seeds) + [acb(0)] * (N - (R - 1)) for seeds in seeds_list]
    for m in range(R, N + 1):
        Arow = []
        for delta in range(1, min(m, DMAX) + 1):
            s = acb(0)
            for j in range(max(0, R - delta), R + 1):
                d = delta - R + j
                if 0 <= d < len(Qs[j]) and Qs[j][d] != (Fr(0), Fr(0)):
                    s += Qa[j][d] * PL.falling(m - delta, j)
            Arow.append(s)
        lead = Qa[R][0] * PL.falling(m, R)
        for c in cs:
            acc = acb(0)
            for delta in range(1, min(m, DMAX) + 1):
                acc += Arow[delta - 1] * c[m - delta]
            c[m] = -acc / lead
    new_seeds, radmax = [], arb(0)
    rr = fr2arb(rho)
    for c in cs:
        Cu = arb(0)
        for jw in range(N - (DMAX + R), N + 1):
            Cu = Cu.max(abs(c[jw]) * rr ** jw)
        Cu = Cu.upper()
        outs = []
        for k in range(deriv_out):
            s = acb(0)
            for m in range(N, k - 1, -1):
                s = s * ha + c[m] * comb(m, k)
            tail_k = Cu * tail_geom_poly(N, k, y, hub)
            outs.append(s + acb(pm_ball(tail_k), pm_ball(tail_k)))
            radmax = radmax.max(arb(outs[-1].real.rad()).max(
                arb(outs[-1].imag.rad())))
        new_seeds.append(outs)
    diag["rad_max"] = float(radmax.upper()) if radmax != 0 else 0.0
    diag["Cu_last"] = float(Cu)
    return new_seeds, diag

def tail_geom_poly(N, k, y, hub):
    """arb ub of sum_{m>N} C(m,k) |h|^{m-k} rho^{-m} (without the C factor)."""
    u = fr2arb(Fr(N + 2, N + 1)) ** k
    r = u * y
    assert r < 1, "landing tail ratio >= 1"
    fk = 1
    for t in range(2, k + 1):
        fk *= t
    return (fr2arb(Fr(N + 1)) ** k * y ** (N + 1) / (1 - r)
            / (fr2arb(Fr(fk)) * fr2arb(hub) ** k)).upper()


# ------------------------------------------------- seeds and route runner
def mum_seeds(thvals, s0):
    """theta^j y(s0) -> Taylor seeds c_k = D^k y(s0)/k!."""
    fp = [[Fr(1)]]
    for k in range(R - 1):
        prev = fp[-1]
        nxt = [Fr(0)] * (len(prev) + 1)
        for j, cj in enumerate(prev):
            nxt[j + 1] += cj
            nxt[j] -= k * cj
        fp.append(nxt)
    s0a = fr2acb(s0)
    seeds = []
    for k in range(R):
        val = acb(0)
        for j, cj in enumerate(fp[k]):
            if cj:
                val += fr2acb(cj) * thvals[j]
        seeds.append(val / (s0a ** k * factorial(k)))
    return seeds


def leg_N(hub, rho, prec):
    y = hub / rho
    import math
    digits = prec * 0.30103 + 12
    N = int(digits * math.log(10) / math.log(float(1 / y))) + 80
    return max(N, 140)


def schedule(s0, hints, s_star):
    """auto-refine waypoints: insert midpoints until every leg fits inside
    0.55 * local majorant radius (the pilot's schedule routine, verbatim)."""
    path = [s0] + hints + [s_star]
    for _ in range(25):
        newpath, changed = [path[0]], False
        for i in range(len(path) - 1):
            x, t = path[i], path[i + 1]
            rho = pick_rho(local_shift(x))
            h = cxq_abs_ub(cxq_sub(t, x))
            if h > Fr(55, 100) * rho:
                mid = cxq(t)
                xx = cxq(x)
                m = ((xx[0] + mid[0]) / 2, (xx[1] + mid[1]) / 2)
                m = (Fr(round(m[0] * 40000), 40000),
                     Fr(round(m[1] * 40000), 40000))
                newpath.append(m); changed = True
            newpath.append(t)
        path = newpath
        if not changed:
            break
    assert not changed, "schedule did not converge"
    return path[1:-1]


def run_route(towers, N_EXT, s0, rho0, waypoints, targets_final, prec):
    ctx.prec = prec
    diags = []
    seeds = []
    for T in towers:
        th, dg = mum_eval(T, N_EXT, s0, rho0, nth=R)
        seeds.append(mum_seeds(th, s0))
        diags.append({"mum_tail": float(dg["tail_worst"])})
    centers = [s0] + waypoints
    for i, ctr in enumerate(centers):
        if i < len(centers) - 1:
            tgt = centers[i + 1]
            h = cxq_sub(tgt, ctr)
            N = leg_N(cxq_abs_ub(h), pick_rho(local_shift(ctr)), prec)
            seeds, dg = local_leg(ctr, seeds, h, N, prec)
            diags.append(dg)
    ctr = centers[-1]
    out = {}
    for tgt in targets_final:
        h = cxq_sub(tgt, ctr)
        N = leg_N(cxq_abs_ub(h), pick_rho(local_shift(ctr)), prec)
        vals, dg = local_leg(ctr, seeds, h, N, prec)
        out[tgt] = vals
        diags.append(dg)
    return out, diags

# ------------------------------------------------------- flux contraction
def sympl(a, b):
    """v^T Sigma w, Sigma = [[0,I],[-I,0]], (h+1)+(h+1) blocks."""
    n = PL.H + 1
    return sum(a[i] * b[n + i] for i in range(n)) \
        - sum(a[n + i] * b[i] for i in range(n))


def contract(vals, tau_im):
    """vals: list over 2h+2 periods of [v_0..v_{r-1}] (v^3-scaled homogeneous
    Pi in ordering (F_0, F_a, X^0, X^a)).  tau pinned externally (FC2:
    curve-restricted stationarity cannot recover tau; tau_hat diagnostic)."""
    h = PL.H
    v3 = acb(0, 2 * arb.pi()) ** 3
    Pi = [vals[i][0] / v3 for i in range(PL.NPER)]
    dPi = [vals[i][1] / v3 for i in range(PL.NPER)]
    Ff = [acb(x) for x in PL.FFLUX]
    Hf = [acb(x) for x in PL.HFLUX]
    A = sympl(Ff, Pi)
    B = sympl(Hf, Pi)
    tau = acb(0, fr2arb(tau_im))
    taubar = tau.conjugate()
    sq2pi = (2 / arb.pi()).sqrt()
    W = sq2pi * (A - tau * B)
    DtauW = sq2pi * (-B) + W * (-1 / (tau - taubar))
    PiC = [z.conjugate() for z in Pi]
    z = sympl(PiC, Pi) * acb(0, -1)
    emKcs = z.real
    assert z.imag.contains(arb(0)), "Pi^dag Sigma Pi not i*real"
    emK = emKcs * 2 * fr2arb(tau_im)
    W0 = abs(W) / emK.sqrt()
    dW = sq2pi * (sympl(Ff, dPi) - tau * sympl(Hf, dPi))
    X0 = Pi[h + 1]
    return {"tau_hat": (A / B).conjugate(),
            "U": [Pi[h + 2 + a] / X0 for a in range(h)],
            "emKcs_over_X0sq": emKcs / abs(X0) ** 2,
            "W": W, "absW": abs(W), "emK": emK, "W0": W0, "dW_ds": dW,
            "DtauW_rel": abs(DtauW) / abs(W), "A": A, "B": B, "Pi": Pi}
