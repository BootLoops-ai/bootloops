#!/usr/bin/env python3
"""selftest_envelope.py — validation battery for envelope_certified.py.

Gates (all exact unless noted):
  E1 step-level norm gates: P1/P3 bounds vs exact induced l1 operator norms,
     both pictures; REGRESSION: the drop-G plant MUST be violated on theta^r
     (step-level bar: an end-to-end numeric acceptance test alone did not
     catch a planted dropped-G error, see PROVEN_MAJORANT.md section 7).
  E2 envelope (C1) exact on: Legendre-type, sunrise banana, cancellation op.
  E3 resonant forced-log tower: REC sanity + envelope past horizon.
  E4 unsoundness exhibits for the empirical x4 scheme (lacunary + oscillatory).
  E5 DKMM order-6 reference tower (read-only, from TERRIER_KKLT_BANK): REC in
     the divided-power picture + certified envelope at N=240.
  E6 comparison table: proven vs empirical-x4 vs measured on 3 reference cases.
Needs the reference data of the two-modulus KKLT example (restricted
operator + certified-W0 towers; not included in the package): set
TERRIER_KKLT_BANK.

Resource caps: single-core, minutes-class; ulimit -v 32505856, nice 5."""
import json, os, sys
from fractions import Fraction as Fr
from math import comb, factorial

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from envelope_certified import (theta_form, mum_jets, rational_roots, phi,
                                certify, deriv_tail, rec_check,
                                empirical_tail_x4, fmt, l1, _jmul, _jpoly, _jinv)

KKLT = os.environ.get("TERRIER_KKLT_BANK")      # reference data root (read-only)
if not KKLT:
    sys.exit("REFUSE: reference data not included in the package (restricted "
             "operator + certified-W0 towers of the two-modulus KKLT "
             "example) -- set TERRIER_KKLT_BANK to run")
FAIL = []

def gate(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name} {detail}")
    if not ok:
        FAIL.append(name)

# ---------------------------------------------------------------- E1
def opmat(P, u, J, picture):
    D = [[Fr(0)] * J for _ in range(J)]
    for j in range(J - 1):
        if picture == "mult":
            D[j + 1][j] = Fr(1)
        else:
            D[j][j + 1] = Fr(1)
    M = [[Fr(0)] * J for _ in range(J)]
    for c in reversed(P):
        T = [[sum(M[i][k] * ((Fr(u) if k == jj else 0) + D[k][jj])
                  for k in range(J)) for jj in range(J)] for i in range(J)]
        M = [[T[i][jj] + (Fr(c) if i == jj else 0) for jj in range(J)]
             for i in range(J)]
    return M

def l1ind(M):
    return max(sum(abs(M[i][j]) for i in range(len(M))) for j in range(len(M)))

def matinv(M):
    n = len(M)
    A = [row[:] + [Fr(1) if i == j else Fr(0) for j in range(n)]
         for i, row in enumerate(M)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[p] = A[p], A[c]
        A[c] = [x / A[c][c] for x in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                A[r] = [x - A[r][c] * y for x, y in zip(A[r], A[c])]
    return [row[n:] for row in A]

def run_E1():
    ok, plant_caught = True, 0
    cases = [([Fr(0), Fr(0), Fr(1)], 0, 2),                  # theta^2 MUM
             ([Fr(0), Fr(0), Fr(0), Fr(1)], 0, 3),           # theta^3 MUM
             ([Fr(0), Fr(-2), Fr(1)], 2, 2)]                 # theta(theta-2), beta=2
    for R0, beta, J in cases:
        lc, roots = rational_roots(R0)
        r = sum(roots.values())
        for pic in ("mult", "div"):
            for m in (beta + 2, beta + 5, beta + 20):
                exact = l1ind(matinv(opmat(R0, m, J, pic)))
                mb = Fr(m) - beta
                G = Fr(1)
                for lam, mu in roots.items():
                    G *= sum(Fr(comb(mu + j - 1, j)) / mb ** j for j in range(J))
                honest = G / (abs(lc) * mb ** r)
                if exact > honest:
                    ok = False
                if R0[-1] == 1 and roots == {Fr(0): r}:      # theta^r: plant check
                    if exact > 1 / (abs(lc) * mb ** r):
                        plant_caught += 1
    gate("E1a inverse-norm bound (P3) both pictures", ok)
    gate("E1b drop-G plant violated on theta^r (bar regression)", plant_caught == 12,
         f"caught {plant_caught}/12")
    # numerator bound (P1/P2)
    ok = True
    Rs = [Fr(-1, 4), Fr(-1), Fr(-1)]                          # -(theta+1/2)^2
    for pic in ("mult", "div"):
        for m in (5, 9, 33):
            exact = l1ind(opmat(Rs, m, 2, pic))
            bnd = sum(abs(c) * (Fr(m) + 1) ** k for k, c in enumerate(Rs))
            if exact > bnd:
                ok = False
    gate("E1c numerator bound (P1/P2) both pictures", ok)

# ---------------------------------------------------------------- E2
def sup_ratio(R, e, J, b, N, t, mmax):
    S = len(R) - 1
    K = max(l1(b[m]) / Fr(t) ** m for m in range(N - S + 1, N + 1))
    if K == 0:
        return Fr(0), K
    return max(l1(b[m]) / (K * Fr(t) ** m) for m in range(N + 1, mmax + 1)), K

def run_E2():
    R_leg = [[Fr(0), Fr(0), Fr(1)], [Fr(-1, 4), Fr(-1), Fr(-1)]]
    R_sun, _ = theta_form([[-3, 1], [9, -20, 3], [0, 9, -10, 1]])
    R_can = [[Fr(0), Fr(1)], [Fr(1), Fr(3)], [Fr(1), Fr(-2)]]
    for tag, R, J, t in (("legendre", R_leg, 2, Fr(11, 10)),
                         ("sunrise", R_sun, 2, Fr(13, 10)),
                         ("cancel", R_can, 1, Fr(36, 10))):
        b = mum_jets(R, 0, J, 160)
        lc, roots = rational_roots(R[0])
        ph = phi(R, 0, J, lc, roots, t, 49)
        w, _ = sup_ratio(R, 0, J, b, 48, t, 160)
        gate(f"E2 {tag}: Phi<=1 and (C1) exact to m=160",
             ph <= 1 and w <= 1, f"Phi={float(ph):.4f} sup={float(w):.4f}")

# ---------------------------------------------------------------- E3
def resonant_jets(N2, J=2, Jb=4):
    """Forced-log tower of L = theta(theta-2) - z(theta+1/2)^2 (beta=2)."""
    R = [[Fr(0), Fr(-2), Fr(1)], [Fr(-1, 4), Fr(-1), Fr(-1)]]
    b = [[Fr(0)] * Jb for _ in range(N2 + 1)]
    b[0][1] = Fr(1)                                   # b_0 = rho
    for m in range(1, N2 + 1):
        acc = _jmul(_jpoly(R[1], m - 1, Jb), b[m - 1], Jb)
        R0j = _jpoly(R[0], m, Jb)
        if m == 2:                                    # divide once by rho
            assert R0j[0] == 0 and acc[0] == 0
            b[m] = [-x for x in _jmul(_jinv(R0j[1:] + [Fr(0)], Jb),
                                      acc[1:] + [Fr(0)], Jb)]
        else:
            b[m] = [-x for x in _jmul(_jinv(R0j, Jb), acc, Jb)]
    return R, [row[:J] for row in b]

def run_E3():
    R, bt = resonant_jets(200)
    bad = rec_check(R, 0, 2, bt, 4, 60, picture="mult")
    gate("E3a resonant tower REC (mult picture, m=4..60)", bad is None)
    cert = certify(R, 0, 2, bt, 64, Fr(1, 2))
    w, _ = sup_ratio(R, 0, 2, bt, 64, cert["t"], 200)
    gate("E3b resonant envelope certify + (C1) to m=200",
         cert["ok"] and w <= 1, f"t={float(cert['t']):.4f} sup={float(w):.4f}")
    return R, bt, cert

# ---------------------------------------------------------------- E4
def run_E4():
    c = [1 if m % 2 == 0 else 0 for m in range(400)]  # 1/(1-z^2)
    emp = empirical_tail_x4(c, 81, Fr(1, 2), 8)
    true = sum(c[m] * 0.5 ** m for m in range(82, 400))
    gate("E4a empirical x4 UNSOUND on lacunary (claims 0)", float(emp) < true,
         f"claimed {float(emp):.1e} true {true:.1e}")
    c = [Fr(1), Fr(6, 5)]
    for m in range(2, 700):
        c.append(Fr(6, 5) * c[-1] - c[-2])            # poles (3+-4i)/5, |.|=1
    emp = empirical_tail_x4(c, 60, Fr(7, 10), 8)
    true = sum(abs(float(ci)) * 0.7 ** m for m, ci in enumerate(c) if m > 60)
    gate("E4b empirical x4 UNSOUND on oscillatory (theta>=0.95 -> 0)", float(emp) < true,
         f"claimed {float(emp):.1e} true {true:.1e}")

# ---------------------------------------------------------------- E5
def load_dkmm():
    """Order-6 L_s theta form + tower-0 (P,Q)=(0,0) slice, divided-power jets,
    read from the reference data root (TERRIER_KKLT_BANK)."""
    OA = json.load(open(os.path.join(KKLT, "restrict", "operator_LS.json")))
    P = {k: [0] * 64 for k in range(7)}
    for key, cval in OA["theta_form"].items():
        d, k = map(int, key.split(","))
        P[k][d] = int(cval)
    R = [[Fr(P[k][d]) for k in range(7)] for d in range(64)]
    T = json.load(open(os.path.join(KKLT, "cert_w0", "towers_ext.json")))
    tw = T["periods"][0]
    J = 4
    layers = {int(k.split(",")[0]): [Fr(v) for v in tw[k]]
              for k in tw if k.endswith(",0,0")}
    n = len(layers[0])
    jets = [[factorial(j) * layers.get(j, [Fr(0)] * n)[m] for j in range(J)]
            for m in range(n)]
    return R, jets

def run_E5():
    R, jets = load_dkmm()
    bad = rec_check(R, 0, 4, jets, 300, 316, picture="div")
    gate("E5a DKMM tower REC (divided picture, m=300..316)", bad is None,
         f"first fail: {bad}")
    cert = certify(R, 0, 4, jets, 240, Fr(1, 112))   # ~R/2; R~1/55.6, l1-cert radius 1/80.9
    ok = cert.get("ok", False)
    w = None
    if ok:
        w, _ = sup_ratio(R, 0, 4, jets, 240, cert["t"], 420)
        ok = ok and w <= 1
    gate("E5b DKMM certify N=240 x=1/112 + (C1) to m=420", ok,
         f"t_min={float(cert.get('t_min', 0)):.3f} sup={float(w) if w is not None else -1:.4f}")
    return R, jets, cert

# ---------------------------------------------------------------- E6
def measured_tail(jets, N, x, Mtop, level):
    return sum(abs(jets[m][level]) * Fr(x) ** m for m in range(N + 1, Mtop + 1))

def table_case(tag, R, e, J, jets, N, x, win, Mtop, f):
    cert = certify(R, e, J, jets, N, Fr(x))
    assert cert["ok"], (tag, cert)
    remc = certify(R, e, J, jets, Mtop - 1, Fr(x))     # certified remainder pad
    rem = remc["tail"] if remc.get("ok") else Fr(0)
    for j in range(J):
        meas = measured_tail(jets, N, x, Mtop, j)
        emp = empirical_tail_x4([jets[m][j] for m in range(N + 1)], N, x, win)
        pv, mv, ev = cert["tail"], meas + rem, emp
        er = fmt(ev / mv) if mv else "n/a"
        pr = fmt(pv / mv) if mv else "n/a"
        f.write(f"| {tag} | {j} | {fmt(mv)} | {fmt(ev)} | {fmt(pv)} | "
                f"{er} | {pr} | {'YES' if ev < mv else 'no'} |\n")
        f.flush()

def run_E6(res_pack, dkmm_pack):
    R_k3, _ = theta_form([[0, -4, 64], [0, 1, -68, 448], [0, 0, 3, -90, 384],
                          [0, 0, 0, 1, -20, 64]])
    b_k3 = mum_jets(R_k3, 0, 3, 512)
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "envelope_comparison.md")
    with open(out, "w") as f:
        f.write("| case | j | measured | empirical x4 | proven | emp/meas | "
                "prov/meas | emp unsound? |\n|---|---|---|---|---|---|---|---|\n")
        table_case("K3-banana N=128 x=1/32", R_k3, 0, 3, b_k3, 128, Fr(1, 32),
                   8, 512, f)
        R_res, bt, _ = res_pack
        table_case("resonant N=64 x=1/2", R_res, 0, 2, bt, 64, Fr(1, 2),
                   8, 200, f)
        R_d, jets_d, _ = dkmm_pack
        table_case("DKMM-L_s N=240 x=1/112", R_d, 0, 4, jets_d, 240, Fr(1, 112),
                   69, 420, f)
    print("E6 comparison table ->", out)
    # C3 corollary sanity on K3, d=2
    cert = certify(R_k3, 0, 3, b_k3, 128, Fr(1, 32))
    bnd = deriv_tail(cert, 2)
    meas = sum(l1(b_k3[m]) * m ** 2 * Fr(1, 32) ** m for m in range(129, 513))
    gate("E6b (C3) derivative tail bound >= exact partial sum", bnd is not None
         and meas <= bnd, f"meas={float(meas):.3e} bound={float(bnd):.3e}")

if __name__ == "__main__":
    run_E1(); run_E2()
    res_pack = run_E3()
    run_E4()
    dkmm_pack = run_E5()
    run_E6(res_pack, dkmm_pack)
    print("SELFTEST:", "FAIL " + ",".join(FAIL) if FAIL else "ALL GATES PASS")
    sys.exit(1 if FAIL else 0)
