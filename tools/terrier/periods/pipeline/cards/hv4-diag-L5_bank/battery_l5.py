#!/usr/bin/env python3
# Part of the hv4-diag-L5 card battery: periods/selftest_slice_operator_card.py replays it and compares the regenerated receipt with the sha-pinned reference receipt, so edits that change its output make that check fail.
# battery_l5.py (DIAG-SLICE) — HV4 diagonal-slice order-5 operator card
# battery: (F) unique-annihilator FIT from the two-engine series + exact
# match vs the printed JKK/GvdH operator (arXiv:2404.12422 ~eq (5.60));
# (M) minimality probes; (S) symbol vs the reference diagonal conifold
# points {1/36, 1/16, 1/4};
# (E) local exponents at {0, 1/36, 1/16, 1/4, inf}; (MUM) recurrence
# regeneration; (P) pfaffian theta-companion payload + jet gate;
# (K3) L3/Domb subsector; (RC) pfaffian.route_choice conformance.
import json, sys, time
from fractions import Fraction as Fr
T0 = time.time()
def pmul(a, b):
    r = [Fr(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b): r[i+j] += x*y
    return r
def padd(a, b):
    r = [Fr(0)]*max(len(a), len(b))
    for i, x in enumerate(a): r[i] += x
    for i, x in enumerate(b): r[i] += x
    return r
def pscal(a, s): return [x*s for x in a]
def peval(a, x): return sum(c*x**i for i, c in enumerate(a))
def ptrim(a):
    while len(a) > 1 and a[-1] == 0: a = a[:-1]
    return a
def pdiv(a, b):          # exact division, returns (q, r)
    a = [Fr(x) for x in a]; b = ptrim([Fr(x) for x in b]); q = [Fr(0)]*max(1, len(a)-len(b)+1)
    while len(ptrim(a)) >= len(b) and any(a):
        a = ptrim(a); d = len(a)-len(b); c = a[-1]/b[-1]; q[d] = c
        for i, y in enumerate(b): a[i+d] -= c*y
        a = ptrim(a)
        if len(a) < len(b): break
    return ptrim(q), ptrim(a)
TH = [Fr(0), Fr(1)]      # theta
def C1(c): return [Fr(c)]
# --- printed operators, expanded exactly from the tex factored forms -----
# L5 (arXiv:2404.12422 ~eq (5.60); cited to JKK 2312.07611)
inner = padd(pscal(pmul(pmul(TH, padd(TH, C1(1))),
                        padd(pmul(TH, padd(TH, C1(1))), C1(1))), 14), C1(3))
Q1 = pscal(pmul(padd(pscal(TH, 2), C1(1)), inner), -2)
th1 = padd(TH, C1(1)); th2 = padd(TH, C1(2))
Q2 = pscal(pmul(pmul(pmul(th1, th1), th1),
                padd(pscal(pmul(TH, th2), 196), C1(255))), 4)
Q3 = pscal(pmul(pmul(pmul(th1, th1), pmul(th2, th2)),
                padd(pscal(TH, 2), C1(3))), -1152)
Q0 = [Fr(0)]*5 + [Fr(1)]
QL5 = [Q0, Q1, Q2, Q3]                       # L5 = sum_j phi^j Q_j(theta)
# L3 (same tex block; cited Verrill 1996 + JKK): theta^3 + 64 phi^2(th+1)^3
#   - 2 phi (2th+1)(5 th(th+1) + 2)
Q3_0 = [Fr(0)]*3 + [Fr(1)]
Q3_1 = pscal(pmul(padd(pscal(TH, 2), C1(1)),
                  padd(pscal(pmul(TH, th1), 5), C1(2))), -2)
Q3_2 = pscal(pmul(pmul(th1, th1), th1), 64)
QL3 = [Q3_0, Q3_1, Q3_2]
# --- series (two-engine gated) ------------------------------------------
SA = json.load(open("SERIES_A.json")); SB = json.load(open("SERIES_B.json"))
assert all(v for k, v in SB.items() if k.startswith("agree")), "engine gate"
c6 = [int(x) for x in SA["c6"]]; domb = [int(x) for x in SA["c4_domb"]]
N = SA["N"]
def rec_ok(Q, c):        # sum_j Q_j(m-j) c_{m-j} == 0 for all m
    for m in range(len(c)):
        s = sum(peval(Q[j], m-j)*c[m-j] for j in range(len(Q)) if m-j >= 0)
        if s != 0: return False, m
    return True, None
ok5, bad5 = rec_ok(QL5, c6); ok3, bad3 = rec_ok(QL3, domb)
assert ok5 and ok3, (bad5, bad3)
# MUM regeneration (Frobenius at 0, analytic solution, y(0)=1)
cre = [Fr(1)]
for m in range(1, N+1):
    s = sum(peval(QL5[j], m-j)*cre[m-j] for j in range(1, 4) if m-j >= 0)
    cre.append(Fr(-s, m**5))
mum_regen_ok = all(cre[m] == c6[m] for m in range(N+1))
# --- FIT: nullspace of the (order, phi-deg) ansatz over Q ----------------
def nullspace(ordr, deg, c, M):
    # unknowns q[j][i], j=0..deg, i=0..ordr ; rows m=0..M
    nu = (deg+1)*(ordr+1); rows = []
    for m in range(M+1):
        row = []
        for j in range(deg+1):
            for i in range(ordr+1):
                row.append(Fr((m-j)**i * c[m-j]) if m-j >= 0 else Fr(0))
        rows.append(row)
    r = 0
    for col in range(nu):
        p = next((i for i in range(r, len(rows)) if rows[i][col] != 0), None)
        if p is None: continue
        rows[r], rows[p] = rows[p], rows[r]
        rows[r] = [x/rows[r][col] for x in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][col] != 0:
                rows[i] = [a - rows[i][col]*b for a, b in zip(rows[i], rows[r])]
        r += 1
    # free columns -> kernel basis
    piv = []; rr = 0
    for col in range(nu):
        if rr < r and rows[rr][col] == 1 and all(rows[k][col] == 0 for k in range(len(rows)) if k != rr):
            piv.append(col); rr += 1
    free = [c_ for c_ in range(nu) if c_ not in piv]
    basis = []
    for f in free:
        v = [Fr(0)]*nu; v[f] = Fr(1)
        for k, col in enumerate(piv): v[col] = -rows[k][f]
        basis.append(v)
    return basis
ker53 = nullspace(5, 3, c6, 60)
fit_unique = (len(ker53) == 1)
v = ker53[0]; v = [x / v[5] for x in v]            # normalize q[0][5] = 1
Qfit = [[v[j*6+i] for i in range(6)] for j in range(4)]
fit_matches_printed = all(ptrim(Qfit[j]) == ptrim(QL5[j]) for j in range(4))
min_ord4 = (len(nullspace(4, 10, c6, 90)) == 0)    # no order<=4, deg<=10
min_deg2 = (len(nullspace(5, 2, c6, 60)) == 0)     # no order-5, deg<=2
ker32 = nullspace(3, 2, domb, 30)                  # K3 sector fit
fitK3_unique = (len(ker32) == 1)
w = ker32[0]; w = [x / w[3] for x in w]
fitK3_matches = all(ptrim([w[j*4+i] for i in range(4)]) == ptrim(QL3[j]) for j in range(3))
# --- theta-normal form R_k(phi), symbol, D-form, exponents ---------------
def Rform(Q, ordr):
    return [[Q[j][k] if k < len(Q[j]) else Fr(0) for j in range(len(Q))]
            for k in range(ordr+1)]                 # R[k] = coeffs in phi
R5 = Rform(QL5, 5); R3 = Rform(QL3, 3)
sym5 = ptrim(R5[5]); sym3 = ptrim(R3[3])
f4, f16, f36 = [Fr(1), Fr(-4)], [Fr(1), Fr(-16)], [Fr(1), Fr(-36)]
symbol_ok = (sym5 == ptrim(pmul(pmul(f4, f16), f36)) and
             sym5 == [Fr(1), Fr(-56), Fr(784), Fr(-2304)] and
             sym3 == ptrim(pmul(f4, f16)) and
             sym5 == ptrim(pmul(sym3, f36)))        # L5 symbol = L3 sym*(1-36phi)
# jet gate (pfaffian R-form annihilation): sum_k R_k(phi) * (theta^k y) = 0
def jet_gate(R, c, ordr):
    vk = [[Fr(n)**k * c[n] for n in range(len(c))] for k in range(ordr+1)]
    for m in range(len(c)):
        s = Fr(0)
        for k in range(ordr+1):
            for j, rj in enumerate(R[k]):
                if rj and m-j >= 0: s += rj*vk[k][m-j]
        if s != 0: return False
    return True
jet5_ok = jet_gate(R5, c6, 5); jet3_ok = jet_gate(R3, domb, 3)
# D-form a_j(phi) = phi^j sum_k S2(k,j) R_k(phi);  a_5 = phi^5 * symbol
S2 = [[Fr(1)]]
for k in range(1, 6):
    S2.append([Fr(0)]*(k+1))
    for j in range(k+1):
        S2[k][j] = (S2[k-1][j-1] if j >= 1 and j-1 < len(S2[k-1]) else Fr(0)) \
                   + (Fr(j)*S2[k-1][j] if j < len(S2[k-1]) else Fr(0))
def Dform(R, ordr):
    a = []
    for j in range(ordr+1):
        s = [Fr(0)]
        for k in range(j, ordr+1): s = padd(s, pscal(R[k], S2[k][j]))
        a.append(ptrim(pmul([Fr(0)]*j + [Fr(1)], s)))
    return a
a5 = Dform(R5, 5); a3 = Dform(R3, 3)
lead_ok = (a5[5] == ptrim(pmul([Fr(0)]*5+[Fr(1)], sym5)) and
           a3[3] == ptrim(pmul([Fr(0)]*3+[Fr(1)], sym3)))  # no apparent sing.
from math import comb
def indicial_at(a, c, ordr):     # exponents at finite regular singular c
    A = []                       # A[j][t] = coeff of x^t in a_j(c+x)
    for j in range(ordr+1):
        aj = a[j]
        A.append([sum(aj[m]*comb(m, t)*Fr(c)**(m-t) for m in range(t, len(aj)))
                  for t in range(len(aj))])
    dmin = min(t-j for j in range(ordr+1) for t in range(len(A[j])) if A[j][t] != 0)
    I = [Fr(0)]
    for j in range(ordr+1):
        t = dmin + j
        if 0 <= t < len(A[j]) and A[j][t] != 0:
            ff = [Fr(1)]
            for i in range(j): ff = pmul(ff, [Fr(-i), Fr(1)])
            I = padd(I, pscal(ff, A[j][t]))
    return dmin, ptrim(I)
def rat_roots(I):                # exact rational roots w/ multiplicity
    roots = []; P = I[:]
    cands = [Fr(p, q) for q in (1, 2, 3, 4, 6) for p in range(-24, 49)]
    prog = True
    while len(P) > 1 and prog:
        prog = False
        for r in cands:
            if peval(P, r) == 0:
                P, rem = pdiv(P, [-r, Fr(1)]); assert rem in ([], [Fr(0)])
                roots.append(r); prog = True; break
    return sorted(roots), ptrim(P)      # residual must be constant
def exps_at(a, c, ordr):
    dmin, I = indicial_at(a, c, ordr)
    roots, resid = rat_roots(I)
    return {"dmin": dmin, "exponents": [str(r) for r in roots],
            "residual_deg": len(resid)-1, "complete": len(roots) == ordr and len(resid) == 1}
E5 = {"0": {"exponents": [str(Fr(0))]*5, "note": "R_k(0) -> nu^5 MUM"},
      "1/36": exps_at(a5, Fr(1, 36), 5), "1/16": exps_at(a5, Fr(1, 16), 5),
      "1/4": exps_at(a5, Fr(1, 4), 5)}
assert ptrim([R5[k][0] for k in range(6)]) == [Fr(0)]*5+[Fr(1)]  # nu^5 at 0
Iinf5 = ptrim([R5[k][3]*(-1)**k for k in range(6)])   # I_inf(nu) = Q3(-nu)
rts_inf5, res5 = rat_roots(Iinf5)
E5["inf"] = {"exponents": [str(r) for r in rts_inf5], "residual_deg": len(res5)-1}
E3 = {"0": {"exponents": [str(Fr(0))]*3, "note": "nu^3 MUM"},
      "1/16": exps_at(a3, Fr(1, 16), 3), "1/4": exps_at(a3, Fr(1, 4), 3)}
Iinf3 = ptrim([R3[k][2]*(-1)**k for k in range(4)])
rts_inf3, res3 = rat_roots(Iinf3)
E3["inf"] = {"exponents": [str(r) for r in rts_inf3], "residual_deg": len(res3)-1}
# --- pfaffian payload (matrix-transport route; scalar-refusal conformance)
import os as _os
sys.path.insert(0, _os.path.abspath(_os.path.join(
    _os.path.dirname(_os.path.abspath(__file__)), *[_os.pardir] * 3)))
from pfaffian import route_choice
rc_point = route_choice("point-value")
rc_scalar = route_choice("scalar-operator", probe=(5, 3))   # 20 slots: OK
den = [int(x) for x in sym5]
Anum = [[[0] for _ in range(5)] for _ in range(5)]
for k in range(4): Anum[k][k+1] = den[:]
for j in range(5): Anum[4][j] = [int(-x) for x in ptrim(R5[j])]
deg_table = [[len(ptrim([Fr(x) for x in e]))-1 for e in row] for row in Anum]
payload = {"basis": ["1", "th", "th^2", "th^3", "th^4"],
 "frame": ("theta-frame: theta v = (A_num/den) v on v=(y,th y,..,th^4 y); "
           "d/dphi frame: dv/dphi = A_num/(phi*den) v; den = symbol "
           "(1-4phi)(1-16phi)(1-36phi); integer numerators"),
 "A_num": Anum, "den": den, "deg_table": deg_table,
 "route_choice_point": rc_point, "route_choice_scalar": rc_scalar}
rec = {"card": "hv4-diag-L5", "scope": "DIAG-SLICE",
 "series_gate": {"two_engine_N": N, "direct6_anchor_N": 40, "pass": True},
 "fit": {"ansatz": "(order 5, phi-deg 3, theta-deg<=5), rows m=0..60",
   "nullspace_dim_1": fit_unique, "matches_printed_L5": fit_matches_printed},
 "minimality": {"no_order_le4_deg_le10": min_ord4,
   "no_order5_deg_le2": min_deg2},
 "recurrence_full_depth": {"L5_n400": ok5, "L3_domb_n400": ok3},
 "mum_regeneration_n400": mum_regen_ok,
 "jet_gate_Rform": {"L5": jet5_ok, "L3": jet3_ok},
 "symbol": {"coeff_line": [str(x) for x in sym5],
   "factors": ["1-4*phi", "1-16*phi", "1-36*phi"],
   "matches_reference_discriminant_points": symbol_ok,
   "leading_Dform_no_apparent_sing": lead_ok},
 "exponents_L5": E5, "exponents_L3": E3,
 "k3_fit": {"nullspace_dim_1": fitK3_unique, "matches_printed_L3": fitK3_matches},
 "pfaffian_payload": payload}
elapsed_s = round(time.time()-T0, 3)   # wall-clock: printed, not written (the receipt is deterministic)
json.dump(rec, open("BATTERY_RECEIPT.json", "w"), indent=1)
gates = [fit_unique, fit_matches_printed, min_ord4, min_deg2, ok5, ok3,
         mum_regen_ok, jet5_ok, jet3_ok, symbol_ok, lead_ok, fitK3_unique,
         fitK3_matches, E5["1/36"]["complete"], E5["1/16"]["complete"],
         E5["1/4"]["complete"]]
print("BATTERY", "PASS" if all(gates) else "FAIL", elapsed_s, "s")
print("E5:", {k: (v["exponents"] if "exponents" in v else v) for k, v in E5.items()})
print("E3:", {k: v["exponents"] for k, v in E3.items()})
assert all(gates)
