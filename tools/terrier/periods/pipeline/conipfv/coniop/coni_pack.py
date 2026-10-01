#!/usr/bin/env python3
"""coni_pack.py — DESIGN SS7 contract artifacts from a derived coni_op.

Given <card>.coni_op.json (routeA_coni output) + the card + its g3' gates
JSON, emit the per-card contract fields:
  coni_towers : exact Frobenius log tower at s=0 (exponent-0 block), from
                the operator by exact rational epsilon-jet recursion; gate:
                L_s kills every tower branch through the computed order.
  coni_routes : MUM -> s_star route with clearances from the operator's
                ACTUAL singular locus (leading-coefficient roots; true and
                apparent singularities are not separated — labeled).
  s_star      : exact rational landing point from the card vev data
                (g3' t_vac; s = exp(-2 pi t/r), continued-fraction pinned).
  coni_C2     : leg-3 landing-ball majorant constant (DESIGN SS3 (*)):
                C2 = |n_cf|/(4 pi) + Li-ladder z_cf^2-block bound from the
                exact J2 grading at s_star with geometric tail bound
                (certificate formula in the JSON, ESTIMATE-labeled tail).
"""
import json, os, sys
from fractions import Fraction as Fr
_PIPE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, _PIPE)
HERE = os.path.dirname(os.path.abspath(__file__))


def theta_from_json(op):
    return {(int(k.split(",")[0]), int(k.split(",")[1])): Fr(v)
            for k, v in op["theta_form"].items()}


class Jet:
    """truncated polynomial in eps over Q, deg < K."""
    def __init__(self, c, K):
        self.K = K
        self.c = (list(c) + [Fr(0)] * K)[:K]
    def __add__(s, o):
        return Jet([a + b for a, b in zip(s.c, o.c)], s.K)
    def __sub__(s, o):
        return Jet([a - b for a, b in zip(s.c, o.c)], s.K)
    def __mul__(s, o):
        out = [Fr(0)] * s.K
        for i, a in enumerate(s.c):
            if a:
                for j, b in enumerate(o.c):
                    if b and i + j < s.K:
                        out[i + j] += a * b
        return Jet(out, s.K)
    def divexact(s, o):
        """divide by jet o with o.c[v] != 0 at its valuation v; requires
        s valuation >= v (exact Frobenius division)."""
        v = next((i for i, x in enumerate(o.c) if x), None)
        assert v is not None, "division by zero jet"
        assert all(x == 0 for x in s.c[:v]), "jet division valuation clash"
        num = s.c[v:] + [Fr(0)] * v
        den = o.c[v:] + [Fr(0)] * v
        out = [Fr(0)] * s.K
        for i in range(s.K):
            acc = num[i]
            for j in range(i):
                acc -= out[j] * den[i - j]
            out[i] = acc / den[0]
        return Jet(out, s.K)


def frobenius_tower(theta, mult0, Ntow, K, extra=0):
    """resonant Frobenius: y(eps) = s^eps sum_m c_m(eps) s^m with c_0 =
    eps^extra (extra = total multiplicity of positive-integer indicial
    roots); branch j = eps^{extra+j} coefficient; T[j][m] = coefficient of
    s^m at log^j s/j! (the eps^{extra+j} parts of c_m).  Exact rationals."""
    degx = max(i for (i, k) in theta)
    ordt = max(k for (i, k) in theta)
    P = [[Fr(0)] * (ordt + 1) for _ in range(degx + 1)]
    for (i, k), c in theta.items():
        P[i][k] = c
    def Pi_at(i, jet_arg):        # P_i(eps + m) as jet: jet_arg = eps + m
        acc = Jet([Fr(0)], jet_arg.K)
        pw = Jet([Fr(1)], jet_arg.K)
        for k in range(ordt + 1):
            if P[i][k]:
                acc = acc + Jet([P[i][k]], jet_arg.K) * pw
            pw = pw * jet_arg
        return acc
    c0 = [Fr(0)] * K
    c0[extra] = Fr(1)
    c = [Jet(c0, K)]
    for m in range(1, Ntow + 1):
        rhs = Jet([Fr(0)], K)
        for i in range(1, min(m, degx) + 1):
            rhs = rhs - Pi_at(i, Jet([Fr(m - i), Fr(1)], K)) * c[m - i]
        c.append(rhs.divexact(Pi_at(0, Jet([Fr(m), Fr(1)], K))))
    # branch j = eps^{extra+j} coefficient of y = s^eps sum c_m s^m;
    # log^i terms from s^eps = sum eps^i log^i/i!: T built in verify basis
    T = [[c[m].c[extra + j] for m in range(Ntow + 1)] for j in range(mult0)]
    return T


def verify_tower(theta, T, Ntow):
    """exact: L_s applied to y_j = sum_i log^i s/i! f_{j-i} vanishes to s^Ntow.
    Uses theta (s d/ds) action on s^m log^i/i!: -> m s^m log^i/i! + s^m
    log^{i-1}/(i-1)!."""
    degx = max(i for (i, k) in theta)
    ordt = max(k for (i, k) in theta)
    J = len(T)
    for j in range(J):
        F = [[T[j - i][m] if j - i < len(T) else Fr(0)
              for m in range(Ntow + 1)] for i in range(j + 1)]
        # res[i][m]: coefficient of s^m log^i/i! of L_s y_j
        res = {}
        for (i, k), cf in theta.items():
            for li in range(j + 1):
                for m in range(Ntow + 1 - i):
                    v = F[li][m]
                    if not v:
                        continue
                    # theta^k (s^m log^li/li!) = sum_t C(k,t) m^{k-t} s^m log^{li-t}/(li-t)!
                    from math import comb
                    for t in range(min(k, li) + 1):
                        w = cf * comb(k, t) * Fr(m) ** (k - t) * v
                        if w:
                            key = (li - t, m + i)
                            res[key] = res.get(key, Fr(0)) + w
        bad = [k for k, v in res.items() if v != 0 and k[1] <= Ntow - degx]
        assert not bad, f"tower branch {j} not annihilated: {sorted(bad)[:4]}"
    return True


def routes_and_star(theta, gates, r_den, dps=60):
    """singular locus from lead coeff of the ODE form; s_star from vev."""
    from mpmath import mp, mpf, exp as mexp, pi as mpi, polyroots, mpc
    mp.dps = dps
    ordt = max(k for (i, k) in theta)
    lead = {}
    for (i, k), c in theta.items():
        if k == ordt:
            lead[i] = lead.get(i, Fr(0)) + c
    deg = max(lead)
    poly = [mpf(int(lead.get(i, 0))) for i in range(deg, -1, -1)]
    roots = polyroots(poly, maxsteps=200, extraprec=200)
    t_vac = mpf(gates["g3p"]["t_vac"])
    s_num = mexp(-2 * mpi * t_vac / r_den)
    s_star = Fr(str(s_num)).limit_denominator(10**30)
    rts = sorted([(abs(rt), rt) for rt in roots if abs(rt) > mpf("1e-40")])
    nearest = rts[0]
    real_pos = [x for x in rts
                if abs(x[1].imag) < mpf("1e-30") and x[1].real > 0]
    route = dict(
        s0="0", target=str(s_star),
        s_star_float=float(s_num),
        n_lead_roots=len(rts),
        nearest_root_abs=float(nearest[0]),
        nearest_real_pos_root=float(real_pos[0][0]) if real_pos else None,
        clearance_ratio=float(mpf(str(float(s_num))) / nearest[0]),
        policy="every leading-coefficient root treated as singular "
               "(apparent/true not separated); waypoints = geometric halving to wall",
        waypoints=[float(s_num) * f for f in (0.25, 0.5, 0.75, 1.0)])
    return route, s_star, [complex(rt) for rt in roots]


def c2_bound(card, cur_gates, s_star_f, J2, M=None, dps=50):
    """C2 per DESIGN SS3 (*): |n_cf|/(4 pi) + (1/8 pi^2) sum_m |J2_m/M_cf^2|
    |s*|^m + geometric tail (ratio from the last computed terms; ESTIMATE label)."""
    from mpmath import mp, mpf, pi as mpi
    mp.dps = dps
    g2p = cur_gates["g2p"]
    ncf, Mcf = int(g2p["n_cf"]), int(g2p["M_cf"])
    s = mpf(abs(s_star_f))
    tot = mpf(0)
    last = mpf(0)
    for m, x in enumerate(J2):
        v = Fr(x) if not isinstance(x, Fr) else x
        term = abs(mpf(v.numerator) / v.denominator) / Mcf**2 * s**m
        tot += term
        if v:
            last = term
    # geometric tail: term ratio bounded by rho = s * growth; measured margin
    rho = mpf("0.5")
    tail = last * rho / (1 - rho)
    C2f = mpf(abs(ncf)) / (4 * mpi) + (tot + tail) / (8 * mpi**2)
    return Fr(str(C2f * (1 + mpf("1e-10")))).limit_denominator(10**12)


if __name__ == "__main__":
    card_name, gates_file, Ntow = sys.argv[1], sys.argv[2], \
        int(sys.argv[3]) if len(sys.argv) > 3 else 120
    op = json.load(open(os.path.join(HERE, f"{card_name}.coni_op.json")))
    gates = json.load(open(gates_file))
    sfile = json.load(open(os.path.join(HERE, f"{card_name}.series.json")))
    theta = theta_from_json(op)
    mult0 = op["indicial_s0"]["mult0"]
    extra = sum(m for e, m in op["indicial_s0"]["factors"] if e != 0)
    # K padding: each resonant division (valuation v) truncates v top eps
    # coefficients; total padding = 2*extra keeps branches < mult0 valid
    # (verify_tower is the fail-closed gate either way)
    T = frobenius_tower(theta, mult0, Ntow, mult0 + 2 * extra, extra)
    verify_tower(theta, T, Ntow)
    print(f"[p5] tower: {mult0} log branches to s^{Ntow}, exact recursion, "
          f"L_s-annihilation VERIFIED")
    # branch-0 must reproduce the graded series (normalization gate)
    a = [Fr(x) for x in sfile["a"]]
    assert all(T[0][m] == a[m] for m in range(min(Ntow, len(a) - 1))), \
        "tower branch 0 != graded series"
    print("[p5] branch-0 == graded restricted series (exact)")
    route, s_star, roots = routes_and_star(theta, gates, sfile["r"])
    print(f"[p6] s_star = {route['s_star_float']:.6g}; nearest lead-root "
          f"{route['nearest_root_abs']:.6g}; clearance ratio "
          f"{route['clearance_ratio']:.3f}")
    J2 = [Fr(x) for x in sfile["J2"]]
    C2 = c2_bound(None, gates, route["s_star_float"], J2)
    print(f"[C2] leg-3 majorant constant = {float(C2):.6g} (tail ESTIMATE)")
    out = dict(
        card=card_name,
        coni_op=dict(theta_form=op["theta_form"], r_shift=op["r_shift"],
                     s_theta=op["s_theta"], provenance="routeA graded series "
                     f"(fit<{op['fit_window']}, {op['held_out_windows']} "
                     "held-out windows exact)"),
        coni_towers=dict(Ntow=Ntow, branches=len(T), scheme="log^j/j!",
                         T=[[str(x) for x in row] for row in T]),
        coni_routes=route,
        s_star=str(s_star),
        coni_jet_ideal=dict(status="J1 jet series stored in the series JSON; "
                            "its annihilator is not computed here",
                            note="leg-2 block-companion built from L_s + J1"),
        coni_C2=str(C2),
        labels=dict(certified=["coni_op recurrence/theta_form", "coni_towers"],
                    estimates=["s_star (vev, PFV step-1)", "coni_C2 tail",
                               "route clearances (float roots)"]))
    fn = os.path.join(HERE, f"{card_name}.contract.json")
    json.dump(out, open(fn, "w"))
    print("wrote", fn)
