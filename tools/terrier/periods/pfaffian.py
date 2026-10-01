#!/usr/bin/env python3
"""pfaffian.py (periods wing) — Pfaffian first-order matrix transport:
(system samples | GKZ ideal + curve) -> CERTIFIED vector at target.

FRONT DOOR — route-choice checklist (pitfall: scalar-elimination blowup; run
it BEFORE any elimination / long series run; pfaffian/DESIGN.md §0 is law):
  1. Is the target a POINT-VALUE / vector-at-a-point?  -> matrix transport.
     Keep the D-module as v' = A(s) v (rank = module rank, entries stay
     height-tame) and transport THAT.  NEVER form the cyclic-vector scalar:
     elimination is the measured blow-up (ads-5-81: rank-17 deg<=574
     connection -> scalar ORDER 13, z-degree > 4052).
  2. Is the scalar operator itself the deliverable?  Only then eliminate —
     and price it FIRST with a mod-p order/degree probe (tools/pf_rank.jl /
     Krylov).  route_choice() REFUSES (ScalarEliminationRefusal, with the
     matrix-route advice) when order*(degree+1) exceeds the budget.
  3. Stored exact series are HELD-OUT gates, never derivation inputs.

Chain (worked per-example scripts in pipeline/pfaffian/ads581/; the
full-size reference receipts are not included in the package —
TERRIER_PFAFFIAN_BANK):
  lift_connection   exact height-tame matrix lift: per-coeff CRT + Wang with
                    the k-1 stability gate (a1_lift pattern)
  fresh_prime_gate  exact A reproduces every sampled entry at a prime NEVER
                    used in the lift (PA-1/PB-2; 0 diffs or FAIL)
  series_annihilation_gate  exact annihilation certificate over Q against a
                    held-out stored series (PA-2/PB-3; transpose law)
  direct_sum_periods + contract_w0  route-A proven-disk direct summation
                    with closed-form tail balls + F3 flux-frame contraction
                    (12-block dictionary fit w/ rank saturation, stored as
                    the F3 frame-coefficient file; verdict.py stage 3)
  entropy_landing   route-B proven-disk landing: entropy-majorant radius
                    (concave f, Euler identity => zero-constant gradient
                    majorant, exact vertex max) + geometric tail + in-ball
                    ODE residual binding the vector to the gated connection
Ball law: arb END-TO-END, no float64 reduction (known pitfall).
Battery: selftest_pfaffian.py — ads-5-81 W0 ball reproduced
byte-comparably from the reference receipts + mutation control."""
import json
from fractions import Fraction as Fr
from math import isqrt, gcd

# largest scalar-elimination we have actually landed is O(10^3-10^4)
# coefficient slots; the ads-5-81 refusal case is 13*(4052+1) ~ 5.3e4.
DEFAULT_SCALAR_BUDGET = 20_000

MATRIX_ADVICE = (
    "REFUSED: scalar elimination prices at {slots} coefficient slots "
    "(order {order} x (degree {degree} + 1) > budget {budget}). The scalar "
    "ODE is the ELIMINATION of the first-order matrix system and generically "
    "explodes (measured: ads-5-81 rank-17/deg-574 connection -> order-13 "
    "scalar with z-degree > 4052). If your target is a certified value at a "
    "point you do not need it: (1) sample the connection A(s) mod p per node "
    "(closure solve is ALGEBRAIC in s0 — no series enumeration); (2) "
    "lift_connection() with the fresh-prime + k-1 stability gates; (3) gate "
    "with series_annihilation_gate() against stored series HELD OUT from the "
    "lift; (4) land with direct_sum_periods()/entropy_landing() inside a "
    "proven disk. See periods/PFAFFIAN.md + pipeline/pfaffian/DESIGN.md.")


class ScalarEliminationRefusal(RuntimeError):
    """Raised by route_choice when a scalar-operator request prices above
    budget (or arrives without a probe). Carries .advice — the matrix
    route."""

    def __init__(self, advice):
        super().__init__(advice)
        self.advice = advice


def route_choice(target, probe=None, budget=DEFAULT_SCALAR_BUDGET):
    """Front door. target: 'point-value' (certified value/vector at a point)
    or 'scalar-operator' (the eliminated ODE itself is the deliverable).
    probe: (order, degree) from a mod-p probe (pf_rank/Krylov) — MANDATORY
    for scalar requests; a scalar request without a probe is refused
    outright."""
    if target in ("point-value", "vector-at-point"):
        return {"route": "matrix-transport",
                "law": "transport v' = A(s) v; never eliminate (DESIGN §0)"}
    if target != "scalar-operator":
        raise ValueError(f"unknown target {target!r}")
    if probe is None:
        raise ScalarEliminationRefusal(
            "REFUSED: scalar-operator request without a mod-p (order, "
            "degree) probe. Price it first (tools/pf_rank.jl / Krylov mod p) "
            "— elimination blow-up is only visible in the probe. If you "
            "actually need a point-value, use route_choice('point-value').")
    order, degree = probe
    slots = order * (degree + 1)
    if slots > budget:
        raise ScalarEliminationRefusal(MATRIX_ADVICE.format(
            slots=slots, order=order, degree=degree, budget=budget))
    return {"route": "scalar-elimination-ok", "slots": slots,
            "budget": budget, "probe": {"order": order, "degree": degree}}


# ------------------- exact height-tame matrix lift (a1_lift law, verbatim)
def crt_pair(a1, m1, a2, m2):
    d = pow(m1, -1, m2)
    return (a1 + ((a2 - a1) * d % m2) * m1) % (m1 * m2), m1 * m2


def wang(r, M):
    """Wang rational reconstruction of r mod M (balanced, B = isqrt(M/2))."""
    B = isqrt(M // 2)
    r0, r1 = M, r % M
    t0, t1 = 0, 1
    while r1 > B:
        q = r0 // r1
        r0, r1, t0, t1 = r1, r0 - q * r1, t1, t0 - q * t1
    if abs(t1) > B or t1 == 0 or gcd(r1, abs(t1)) != 1 or gcd(t1, M) != 1:
        return None
    return Fr(r1, t1) if t1 > 0 else Fr(-r1, -t1)


def lift_arr(vals_by_prime, plist):
    """per-coeff CRT over plist then Wang. -> (fracs-with-None, nfail)."""
    n = len(vals_by_prime[0])
    out, nfail = [], 0
    for k in range(n):
        a, m = int(vals_by_prime[0][k]) % plist[0], plist[0]
        for t in range(1, len(plist)):
            a, m = crt_pair(a, m, int(vals_by_prime[t][k]) % plist[t],
                            plist[t])
        f = wang(a, m)
        if f is None:
            nfail += 1
            out.append(None)
        else:
            out.append(f)
    return out, nfail


def lift_connection(banks, primes, seed=20260711, subsample=20):
    """Exact lift of a sampled connection: banks = per-prime dicts with keys
    DEN (len dden+1 int array), Anum (RK x RK x maxdeg+1), degs (RK x RK,
    -1 = zero entry), dden. Adaptive-schedule acceptance = 0 Wang fails AND
    the k-1 stability gate (all-but-last-prime lift agrees on a 1/subsample
    coeff subsample). Returns the A_exact dict (a1_lift schema)."""
    import random
    degs = banks[0]["degs"]
    dden = int(banks[0]["dden"])
    for b in banks[1:]:
        assert (b["degs"] == degs).all() if hasattr(degs, "all") \
            else b["degs"] == degs, "deg table mismatch"
        assert int(b["dden"]) == dden
    RK = len(degs)
    DENf, dfail = lift_arr([b["DEN"] for b in banks], primes)
    entries, nfail_tot, ncoef = {}, dfail, len(DENf)
    hmax = max((abs(f.numerator).bit_length() + f.denominator.bit_length()
                for f in DENf if f is not None), default=0)
    for i in range(RK):
        for j in range(RK):
            d = int(degs[i][j]) if not hasattr(degs, "shape") \
                else int(degs[i, j])
            if d < 0:
                continue
            fr, nf = lift_arr([b["Anum"][i, j, :d + 1] for b in banks],
                              primes)
            nfail_tot += nf
            ncoef += d + 1
            if nf == 0:
                hmax = max(hmax, max(abs(f.numerator).bit_length()
                                     + f.denominator.bit_length()
                                     for f in fr))
            entries[f"{i},{j}"] = fr
    if nfail_tot:
        return {"incomplete": True, "wang_fails": nfail_tot,
                "ncoef": ncoef, "advice": "harvest more primes (A-1 law)"}
    random.seed(seed)
    idx = [(i, j, k) for i in range(RK) for j in range(RK)
           if int(degs[i, j]) >= 0 for k in range(int(degs[i, j]) + 1)]
    sub = random.sample(idx, max(1, ncoef // subsample))
    mism = 0
    if len(primes) - 1 >= 2:
        for (i, j, k) in sub:
            a, m = None, None
            for t, p in enumerate(primes[:-1]):
                v = int(banks[t]["Anum"][i, j, k]) % p
                a, m = (v, p) if a is None else crt_pair(a, m, v, p)
            if wang(a, m) != entries[f"{i},{j}"][k]:
                mism += 1
    else:
        mism = -1                                  # k-1 gate needs k >= 3
    return {"primes": list(primes), "dden": dden, "rank": RK,
            "den": [str(f) for f in DENf],
            "deg_table": [list(map(int, r)) for r in degs],
            "entries": {k: [str(f) for f in v] for k, v in entries.items()},
            "max_height_bits": hmax, "stability_pass": (mism == 0)}


# ---------------- fresh-prime gate (PA-1 / PB-2 core, generalized verbatim)
def fresh_prime_gate(AX, P, sig, A, require_fresh=True):
    """Exact A(s) must reproduce EVERY sampled entry value at nodes sig of a
    prime P never used in the lift. sig: int array of nodes; A: RK x RK x N
    mod-P samples. Returns receipt dict; pass iff diffs == 0."""
    import numpy as np
    if require_fresh:
        assert P not in AX["primes"], "fresh prime was used in the lift!"

    def red(fs):
        out = np.zeros(len(fs), dtype=np.int64)
        for k, s in enumerate(fs):
            f = Fr(s)
            out[k] = f.numerator % P * pow(f.denominator % P, P - 2, P) % P
        return out

    def pev(pl, xs):
        acc = np.zeros(len(xs), dtype=np.int64)
        for cf in pl[::-1]:
            acc = (acc * xs + int(cf)) % P
        return acc

    sig = sig.astype(np.int64)
    A = A.astype(np.int64)
    N, RK = len(sig), A.shape[0]
    dv = pev(red(AX["den"]), sig)
    assert np.all(dv % P != 0), "den vanishes at a fresh node"
    dinv = np.array([pow(int(v), P - 2, P) for v in dv], dtype=np.int64)
    ndiff = 0
    for i in range(RK):
        for j in range(RK):
            key = f"{i},{j}"
            if key in AX["entries"]:
                v = pev(red(AX["entries"][key]), sig) * dinv % P
            else:
                v = np.zeros(N, dtype=np.int64)
            ndiff += int(np.count_nonzero((v - A[i, j]) % P))
    return {"gate": "fresh-prime", "fresh_prime": P, "nodes": N,
            "entry_values_checked": RK * RK * N, "diffs": ndiff,
            "pass": ndiff == 0}


def series_annihilation_gate(AX, rows, M, transpose=True):
    """Exact annihilation certificate over Q against a HELD-OUT stored
    series: den(s)*theta_s v - A^T v == 0 through s^M, all rows (PB-3/PA-2
    law; transpose = module-side convention, measured on a short prefix
    first — the untransposed residual fails loudly at orders 2-5)."""
    from flint import fmpq_poly, fmpq
    RK = AX["rank"]

    def q2f(s):
        f = Fr(s)
        return fmpq(f.numerator, f.denominator)

    den = fmpq_poly([q2f(x) for x in AX["den"]])
    vs, tv = [], []
    for j in range(RK):
        c = [q2f(x) for x in rows[j]]
        vs.append(fmpq_poly(c))
        tv.append(fmpq_poly([fmpq(m) * c[m] for m in range(M + 1)]))
    bad = []
    for i in range(RK):
        r = den * tv[i]
        for j in range(RK):
            e = AX["entries"].get(f"{j},{i}" if transpose else f"{i},{j}")
            if e is None:
                continue
            r = r - fmpq_poly([q2f(x) for x in e]) * vs[j]
        cs = r.coeffs()
        if [m for m in range(min(M, len(cs) - 1) + 1) if cs[m] != 0]:
            bad.append(i)
    return {"gate": "series-annihilation", "M": M, "rank": RK,
            "rows_clean": RK - len(bad), "bad_rows": bad,
            "pass": not bad}


# --------- route A: proven-disk direct summation + flux-frame contraction
# (verdict.py stages 3-5 verbatim, parameterized; arb end-to-end)
def fr2arb(x):
    from flint import arb, fmpq
    return arb(fmpq(x.numerator, x.denominator))


def col_tower(phi_al, p, q, M):
    """F3 dictionary column: Phi_alpha tower shifted into the (p, q) slot,
    key (lpow, p+2*e2, q+e3), e2-sector weight (-1/24)^e2 (verdict stage 3)."""
    out = {}
    for (lp, e2, e3), v in phi_al.items():
        key = (lp, p + 2 * e2, q + e3)
        dst = out.setdefault(key, [Fr(0)] * (M + 1))
        f = Fr(-1, 24) ** e2
        for m in range(M + 1):
            if v[m]:
                dst[m] += f * v[m]
    return out


def frame_towers(coeffs, phi, M):
    """Rebuild the fitted period towers from stored frame coefficients
    (F3 frame schema: per-period dict (alpha, p, q) -> Fr), preserving
    the stored coefficient order (byte-parity law)."""
    towers = []
    for cd in coeffs:
        T = {}
        for (al, p, q), cv in cd.items():
            for key, v in col_tower(phi[al], p, q, M).items():
                dst = T.setdefault(key, [Fr(0)] * (M + 1))
                for m in range(M + 1):
                    if v[m]:
                        dst[m] += cv * v[m]
        towers.append({k: v for k, v in T.items()
                       if any(x != 0 for x in v)})
    return towers


def tail_ball(cmax_abs, s_abs, M, bound_level):
    """sum_{m>M} bound_level(m)*prefmax*s^m, geometric domination with the
    conservative ratio r = 2*s_abs*1.5 (valid M >= 20; verdict law)."""
    from flint import arb
    r = 2 * fr2arb(s_abs) * arb(1.5)
    assert r < 1
    t = bound_level(M + 1) * fr2arb(s_abs) ** (M + 1) / (1 - r)
    return (t * cmax_abs).upper()


def direct_sum_periods(towers, coeffs, svac, mtow, dps, bound_level,
                       pref_bound=8000, prec_guard=20):
    """Certified period vector by DIRECT SUMMATION inside a proven disk:
    every tower summed to s^mtow with the closed-form tail ball; branch
    log s = ln|s| + i*pi for svac < 0 (eps = +1). Returns (Pi, v=2*pi*i)."""
    from flint import arb, acb, ctx
    ctx.prec = int(dps * 3.33) + prec_guard
    sA = fr2arb(abs(svac))
    L = acb(sA.log(), arb.pi()) if svac < 0 else acb(sA.log())
    v = acb(0, 2 * arb.pi())
    z3 = arb.zeta(arb(3))
    sv = fr2arb(svac)
    Pi = []
    for i in range(len(towers)):
        val = acb(0)
        cmax = arb(0)
        for u, cv in coeffs[i].items():
            cmax += fr2arb(abs(cv)) * pref_bound
        tb = tail_ball(cmax, abs(svac), mtow, bound_level)
        for (lp, P, Q), cf in towers[i].items():
            s = arb(0)
            for m in range(mtow, -1, -1):
                s = s * sv + fr2arb(cf[m])
            pm = arb(0, tb)
            val += acb(s + pm, pm) * L ** lp * v ** P * z3 ** Q
        Pi.append(val)
    return Pi, v


def bound_level_ads581(m):
    """Rigorous |C_beta[m]| bound for the ads-5-81 jet towers (verdict.py
    law, proven elementary: multinomial 2^m, point count, psi-jet Faa di
    Bruno majorants B1..B3, x8 sector count)."""
    from flint import arb
    npts = (fr2arb(Fr(m, 9)) + 1) * (fr2arb(Fr(m, 10)) + 1) ** 3
    B1 = 26 * (1 + (arb(1 + 4 * m)).log())
    B2 = B1 ** 2 + 264
    B3 = B1 ** 3 + 3 * B1 * 264 + 2880
    return npts * arb(2) ** m * (1 + B1 + B2 + B3) * 8


def sympl(a, b, n):
    """standard symplectic pairing <a,b> on 2n-vectors (PT.contract law)."""
    return sum(a[i] * b[n + i] for i in range(n)) \
        - sum(a[n + i] * b[i] for i in range(n))


def contract_w0(Pi, v, F, H, tau_re, tau_im, n=6):
    """Flux-frame contraction (verdict stage 5 verbatim): Pi/v^3, symplectic
    A = <F,Pi>, B = <H,Pi>, W = sqrt(2/pi) |A - tau B| (UNDRESSED, paper
    09064 convention). tau_im: exact Fraction pin. Returns dict."""
    from flint import arb, acb
    v3 = v ** 3
    Piv = [p / v3 for p in Pi]
    Ff = [acb(x) for x in F]
    Hf = [acb(x) for x in H]
    A = sympl(Ff, Piv, n)
    B = sympl(Hf, Piv, n)
    tau = acb(arb(tau_re), fr2arb(tau_im))
    sq2pi = (2 / arb.pi()).sqrt()
    W = sq2pi * (A - tau * B)
    return {"W0": abs(W), "A": A, "B": B, "tau_hat": (A / B)}


# --------- route B: entropy-majorant certified landing
def entropy_landing(V, VERTS, YT, s_star, M, rows, dps, point_count,
                    growth_pow, A_exact=None, prec_guard=40,
                    tail_slack=1.001):
    """Certified ball vector at s_star by direct summation of the stored
    exact graded series with the ENTROPY-MAJORANT tail: c(q) <= e^{f(a)},
    f concave + 1-homogeneous => zero-constant gradient majorant at interior
    YT (Euler identity, live receipt c0 ~ 0), maximized EXACTLY over the
    level-1 support-polytope vertices => radius e^{-Phi1}; point count +
    weight polynomial point_count(M+1) (arb callable); geometric tail with
    ratio qhat*((M+2)/(M+1))^growth_pow. If A_exact given, the landing ODE
    residual den*theta_v - A^T v must contain 0 in-ball, every row."""
    from flint import arb, ctx
    ctx.prec = int(dps * 3.33) + prec_guard
    f2a = fr2arb
    h = len(YT)
    at = [sum(Fr(V[i][j]) * YT[j] for j in range(h)) for i in range(len(V))]
    assert all(x > 0 for x in at), "y~ not interior"
    Aq = sum(at)
    g = [(f2a(Aq) / f2a(x)).log() for x in at]
    c0 = (f2a(Aq) * f2a(Aq).log()
          - sum(f2a(x) * f2a(x).log() for x in at)
          - sum(g[i] * f2a(at[i]) for i in range(len(V))))
    Phi1 = None
    for vtx in VERTS:
        val = sum(g[i] * sum(f2a(V[i][j]) * f2a(vtx[j]) for j in range(h))
                  for i in range(len(V)))
        Phi1 = val if Phi1 is None else arb.max(Phi1, val)
    sst = f2a(s_star)
    qhat = (Phi1.exp() * sst).upper()
    assert float(qhat) < 1
    mp1 = M + 1
    Pm = point_count(mp1)
    ratio = arb(qhat) * (f2a(M + 2) / f2a(M + 1)) ** growth_pow
    assert ratio < 1, "tail ratio >= 1"
    T = (c0.exp() * Pm * arb(qhat) ** mp1 / (1 - ratio)).upper()
    tb = arb(0, float(T) * tail_slack)
    out_v, out_tv = [], []
    RK = len(rows)
    for b in range(RK):
        cs = [f2a(Fr(x)) for x in rows[b]]
        acc, acct = arb(0), arb(0)
        for m in range(M, -1, -1):
            acc = acc * sst + cs[m]
            acct = acct * sst + m * cs[m]
        out_v.append(acc + tb)
        out_tv.append(acct + tb)
    res = {"v": out_v, "theta_v": out_tv, "tail": float(T), "c0": c0,
           "Phi1": Phi1, "radius": float((-Phi1).exp())}
    if A_exact is not None:
        denv = arb(0)
        for c in [f2a(Fr(x)) for x in A_exact["den"]][::-1]:
            denv = denv * sst + c
        worst = 0.0
        for i in range(RK):
            r = denv * out_tv[i]
            for j in range(RK):
                e = A_exact["entries"].get(f"{j},{i}")
                if e is None:
                    continue
                ev = arb(0)
                for c in [f2a(Fr(x)) for x in e][::-1]:
                    ev = ev * sst + c
                r -= ev * out_v[j]
            assert r.contains(arb(0)), f"ODE residual row {i} excludes 0"
            worst = max(worst, float(r.rad()))
        res["ode_residual_worst_rad"] = worst
    return res
