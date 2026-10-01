#!/usr/bin/env python3
"""W4 blind leg: finite-g mixture partition function Z_g and its DPM limit,
implemented INDEPENDENTLY from the statements only (blind protocol: no W4
files read). Generating-function route, exact Fractions throughout.

Notation: k states, g components, U=(u_1..u_k) state counts, N=sum U,
NORMALIZED Dir(1) state priors, Z(U) with NO multinomial coefficient
(fixed observation sequence).

(A) Z_g(U;a) = [prod u_v!]/(a)_N * sum_{weak comps m of N, len g}
        CT_U(m) * prod_a (a/g)_{m_a} * (k-1)!/(m_a+k-1)!
(B) Derived here from EPPF x Dir(1)-multinomial (blocks B_i, |B_i|=m_i,
    EPPF weight a^t prod (m_i-1)!/(a)_N; block state counts c_i get
    (k-1)! prod_v c_vi! / (m_i+k-1)!; collapsing observation->block
    assignments prod_v u_v!/prod_i c_vi! per table, /t! for block order,
    the c! cancel):
    Z_DPM(U;a) = [prod u_v!]/(a)_N * sum_{t=1}^{N} (a^t/t!) *
        sum_{comps m of N, m_i>=1} CT_U(m) * prod_i (m_i-1)!(k-1)!/(m_i+k-1)!
(C) Claim: Z_g -> Z_DPM and g*(Z_g - Z_DPM) -> constant.

GF mechanism: CT_U(m) = [prod_a z_a^{m_a}] prod_v h_{u_v}(z_1..z_g).
Contracting each column against its weight w(m) gives one k-variate
column polynomial F(x) = sum_{0<=c<=U} w(|c|) x^c; then
  Z_g/pref  = [x^U] F_g(x)^g            (binary powering, truncated to box U)
  Z_DPM/pref= [x^U] sum_{t=1}^N (a F1)^t/t!   (truncated exponential)
"""
from fractions import Fraction as Fr
from itertools import product
from math import factorial as fact
import json, os, sys, time


def rising(x, n):
    r = Fr(1)
    for i in range(n):
        r *= (x + i)
    return r


# ---------- truncated multivariate polynomial arithmetic (box 0<=c<=U) ----
def pmul(A, B, U):
    out = {}
    for ca, va in A.items():
        for cb, vb in B.items():
            c = tuple(x + y for x, y in zip(ca, cb))
            ok = True
            for ci, ui in zip(c, U):
                if ci > ui:
                    ok = False
                    break
            if ok:
                out[c] = out.get(c, Fr(0)) + va * vb
    return out


def ppow(A, n, U):
    k = len(U)
    R = {tuple([0] * k): Fr(1)}
    base = A
    while n:
        if n & 1:
            R = pmul(R, base, U)
        n >>= 1
        if n:
            base = pmul(base, base, U)
    return R


def boxpoints(U):
    return product(*[range(u + 1) for u in U])


def prefactor(U, alpha):
    p = Fr(1)
    for u in U:
        p *= fact(u)
    return p / rising(alpha, sum(U))


# ---------- (A) finite-g via GF ----------
def Zg_gf(U, alpha, g):
    k = len(U)
    F = {}
    for c in boxpoints(U):
        m = sum(c)
        F[c] = rising(Fr(alpha, g) if isinstance(alpha, int) else alpha / g, m) \
            * Fr(fact(k - 1), fact(m + k - 1))
    P = ppow(F, g, U)
    return prefactor(U, alpha) * P.get(tuple(U), Fr(0))


# ---------- (B) DPM via truncated exp GF ----------
def Zdpm_gf(U, alpha, drop_tfact=False):
    """drop_tfact=True is the G4 planted error: omit the 1/t!."""
    k = len(U)
    N = sum(U)
    F1 = {}
    for c in boxpoints(U):
        m = sum(c)
        if m >= 1:
            F1[c] = Fr(fact(m - 1) * fact(k - 1), fact(m + k - 1))
    total = {}
    term = {tuple([0] * k): Fr(1)}
    for t in range(1, N + 1):
        term = pmul(term, F1, U)
        w = alpha ** t if drop_tfact else (alpha ** t) / fact(t)
        for c, v in term.items():
            total[c] = total.get(c, Fr(0)) + w * v
    return prefactor(U, alpha) * total.get(tuple(U), Fr(0))


# ---------- G1 reference: direct allocation brute force ----------
def Zg_brute(U, alpha, g):
    """Straight from the model: allocations z in {1..g}^N with
    Dirichlet-multinomial(alpha/g) sequence probability, per-component
    Dir(1)-multinomial marginal, fixed observation sequence."""
    k = len(U)
    N = sum(U)
    x = []
    for v, u in enumerate(U):
        x += [v] * u
    ag = alpha / g if not isinstance(alpha, int) else Fr(alpha, g)
    total = Fr(0)
    for z in product(range(g), repeat=N):
        m = [0] * g
        c = [[0] * k for _ in range(g)]
        for obs, a in zip(x, z):
            m[a] += 1
            c[a][obs] += 1
        p = Fr(1)
        for a in range(g):
            p *= rising(ag, m[a])
            num = fact(k - 1)
            for v in range(k):
                num *= fact(c[a][v])
            p *= Fr(num, fact(m[a] + k - 1))
        total += p
    return total / rising(alpha, N)


# ---------- G2 reference: CRP sequential-predictive brute force ----------
def set_partitions(items):
    if not items:
        yield []
        return
    first, rest = items[0], items[1:]
    for part in set_partitions(rest):
        for i in range(len(part)):
            yield part[:i] + [[first] + part[i]] + part[i + 1:]
        yield [[first]] + part


def Zdpm_brute(U, alpha):
    """Enumerate set partitions of {1..N}; probability of each computed by
    the CRP sequential predictive x Dir(1) posterior predictive, item by
    item -- a mechanism disjoint from the GF route."""
    k = len(U)
    N = sum(U)
    x = []
    for v, u in enumerate(U):
        x += [v] * u
    total = Fr(0)
    for part in set_partitions(list(range(N))):
        assign = {}
        for bi, B in enumerate(part):
            for item in B:
                assign[item] = bi
        p = Fr(1)
        seen = {}  # block -> [m, counts]
        for n_prev in range(N):
            b = assign[n_prev]
            v = x[n_prev]
            if b in seen:
                mi, ci = seen[b]
                p *= Fr(mi) / (alpha + n_prev)
                p *= Fr(ci[v] + 1, mi + k)
                seen[b][0] += 1
                ci[v] += 1
            else:
                p *= alpha / (alpha + n_prev)
                p *= Fr(1, k)
                seen[b] = [1, [0] * k]
                seen[b][1][v] += 1
        total += p
    return total


# ---------- gates ----------
def run_gates():
    t0 = time.time()
    gate = {"gate": "GATE_W4_DPM", "run_id": "w4-dpm-blind",
            "date": "2026-07-18", "blind_protocol": True,
            "implementation": "w4_blind_gf.py (generating functions, exact Fractions)",
            "formula_disagreements": []}

    a13, a32, a227 = Fr(1, 3), Fr(3, 2), Fr(22, 7)

    # --- G1: Z_g GF vs direct allocation brute force ---
    g1_cases = [
        ((2, 1), Fr(1), 2), ((1, 1, 1), a13, 3), ((3, 2), a32, 2),
        ((2, 2), a227, 3), ((2, 1, 1), a13, 4), ((3, 1), Fr(2), 5),
        ((2, 2, 1), a32, 3), ((1, 2, 3), a227, 2),
    ]
    g1 = []
    for U, al, g in g1_cases:
        zg = Zg_gf(U, al, g)
        zb = Zg_brute(U, al, g)
        g1.append({"U": list(U), "alpha": str(al), "g": g,
                   "Z_gf": str(zg), "Z_brute": str(zb), "match": zg == zb})
    gate["G1"] = {"desc": "Z_g GF vs fresh direct-allocation brute force",
                  "cases": g1, "pass": all(c["match"] for c in g1)}

    # --- G2: Z_DPM GF vs CRP set-partition brute force ---
    g2_cases = [
        ((2, 1), Fr(1)), ((1, 1, 1), a13), ((3, 2), a32), ((2, 2), a227),
        ((2, 2, 2), Fr(1, 2)), ((4, 3), a227), ((3, 3, 2), a13), ((5, 1), a32),
    ]
    g2 = []
    for U, al in g2_cases:
        zd = Zdpm_gf(U, al)
        zb = Zdpm_brute(U, al)
        g2.append({"U": list(U), "alpha": str(al),
                   "Z_gf": str(zd), "Z_brute": str(zb), "match": zd == zb})
    gate["G2"] = {"desc": "Z_DPM GF vs fresh CRP sequential-predictive brute force",
                  "cases": g2, "pass": all(c["match"] for c in g2)}

    # --- G3: convergence g*(Z_g - Z_DPM) at g=64,128,256 ---
    g3_cases = [((2, 1), Fr(1)), ((2, 2), a13), ((3, 2), a32), ((2, 1, 1), a227)]
    g3 = []
    for U, al in g3_cases:
        zd = Zdpm_gf(U, al)
        D = {}
        for g in (64, 128, 256):
            D[g] = g * (Zg_gf(U, al, g) - zd)
        d1 = D[64] - D[128]
        d2 = D[128] - D[256]
        # D_g = C + c1/g + O(1/g^2) => successive differences halve;
        # d1 == d2 == 0 means D_g is exactly constant (strongest stability)
        ratio = float(d1 / d2) if d2 != 0 else None
        stable = (d1 == 0 and d2 == 0) or \
            ((ratio is not None) and (1.7 < ratio < 2.3))
        # Richardson extrapolation for the constant
        C_extrap = D[256] + (D[256] - D[128])
        g3.append({"U": list(U), "alpha": str(al),
                   "gD_64": float(D[64]), "gD_128": float(D[128]),
                   "gD_256": float(D[256]),
                   "diff_ratio_(64-128)/(128-256)": ratio,
                   "C_extrapolated": float(C_extrap),
                   "doubling_stable": stable})
    gate["G3"] = {"desc": "g*(Z_g - Z_DPM) doubling-stability, g=64,128,256",
                  "cases": g3, "pass": all(c["doubling_stable"] for c in g3)}

    # --- G4: planted error (drop 1/t!) must be caught by G2 machinery ---
    g4 = []
    for U, al in [((3, 2), a32), ((2, 2), a227)]:
        zbad = Zdpm_gf(U, al, drop_tfact=True)
        zb = Zdpm_brute(U, al)
        zgood = Zdpm_gf(U, al)
        g4.append({"U": list(U), "alpha": str(al),
                   "Z_planted(no 1/t!)": str(zbad), "Z_brute": str(zb),
                   "planted_detected": zbad != zb,
                   "good_still_matches": zgood == zb})
    gate["G4"] = {"desc": "planted error in own impl (dropped 1/t!) caught by G2 brute",
                  "cases": g4,
                  "pass": all(c["planted_detected"] and c["good_still_matches"] for c in g4)}

    gate["all_pass"] = all(gate[G]["pass"] for G in ("G1", "G2", "G3", "G4"))
    gate["runtime_sec"] = round(time.time() - t0, 2)
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "gates", "GATE_W4_DPM.json")
    with open(out, "w") as f:
        json.dump(gate, f, indent=1)
    print(json.dumps({k: (v if not isinstance(v, dict) else v.get("pass"))
                      for k, v in gate.items() if k.startswith("G") or k == "all_pass"}))
    print("recorded:", out)
    return gate


if __name__ == "__main__":
    run_gates()
