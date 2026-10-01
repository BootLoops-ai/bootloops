#!/usr/bin/env python3
"""Swap route, k=2 ONLY: log-tilted float64 evaluation of ln Z via the
substitution [s^u1 t^N] F(s t, t)^g (finite g, Dir(1) mixing = W1 register)
and ln Z_DPM via series exp of alpha*F1 (EPPF register, alpha=1 conventions).

Registers (pre-verified):
  finite g, Dir(1): column weight w(m) = m!(k-1)!/(m+k-1)! = 1/(m+1) at k=2;
      exact reference Zg_gf(U, alpha=g, g) from w4_blind_gf.py (Dir(alpha/g)
      with alpha=g equals Dir(1), verified identity).
  DPM: per-block weight (m-1)!(k-1)!/(m+k-1)! = 1/(m(m+1)) at k=2;
      exact reference Zdpm_gf(U, alpha).
  Z has NO multinomial coefficient; prefactor = u1! u2! / rising(alpha, N).

Mechanics: 2D coefficient array A[i,j] = coeff of s^i t^j (i<=u1, j<=N,
band mask i<=j, j-i<=u2); exponential tilting t -> e^theta t chosen by
bisection so the corner (u1, N) sits at the tilted mode; binary powering
with per-squaring max-normalization (log scale accumulated); FFT
convolution (numpy.fft.rfft2) for products; DPM via scaled Taylor exp
(Horner, T terms, l1-norm <= 0.5 after 2^d scaling) then d squarings.

CLI (all results recorded incrementally):
  selftest                          tiny-N float-vs-exact smoke test
  gate_g N g [label]                finite-g gate vs Zg_gf(alpha=g)
  gate_dpm N alpha [label]          DPM gate vs Zdpm_gf
  gate_planted N g ref_label        planted error (skip one renorm log)
  row_g N g [label]                 timed row, appended to SWAP_ROWS.md
  row_dpm N alpha [label]           timed row, appended to SWAP_ROWS.md
  probe                             squaring walls + N=50000 projection
"""
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
GATE_PATH = os.path.join(HERE, "gates", "GATE_SWAP_ROUTE.json")
ROWS_PATH = os.path.join(HERE, "SWAP_ROWS.md")
LN2 = math.log(2.0)


# ---------------- utilities ----------------
def next_fast_len(n):
    best = None
    p2 = 1
    while p2 < 8 * n:
        p3 = p2
        while p3 < 8 * n:
            p5 = p3
            while p5 < n:
                p5 *= 5
            if p5 >= n and (best is None or p5 < best):
                best = p5
            p3 *= 3
        p2 *= 2
    return best


def lnint(n):
    b = n.bit_length()
    if b <= 900:
        return math.log(n)
    sh = b - 900
    return math.log(n >> sh) + sh * LN2


def lnfr(fr):
    return lnint(fr.numerator) - lnint(fr.denominator)


def bank_gate(entry):
    if os.path.exists(GATE_PATH):
        with open(GATE_PATH) as f:
            doc = json.load(f)
    else:
        doc = {"gate": "GATE_SWAP_ROUTE", "run_id": "swap-route-w4",
               "engine": "swap_route.py (k=2, log-tilted float64, FFT powering)",
               "exact_twins": "w4_blind_gf.py Zg_gf(U, alpha=g, g) [Dir(1) identity], Zdpm_gf(U, alpha)",
               "criterion": "|lnZ_float - lnZ_exact|/|lnZ_exact| < 1e-10",
               "entries": []}
    entry["utc"] = time.strftime("%Y-%m-%d %H:%M:%S")
    doc["entries"].append(entry)
    with open(GATE_PATH, "w") as f:
        json.dump(doc, f, indent=1)
    print("recorded gate entry:", entry.get("label"))


def bank_row(lines):
    new = not os.path.exists(ROWS_PATH)
    with open(ROWS_PATH, "a") as f:
        if new:
            f.write("# Swap-route rows (k=2, balanced U, log-tilted float64 + FFT powering)\n\n"
                    "Engine: swap_route.py. Every row k=2, balanced U (u_v = N/k). "
                    "Finite-g rows use Dir(1,..,1) mixing (W1-gated register: weight per "
                    "column-sum m is m!(k-1)!/(m+k-1)!); DPM rows use the EPPF register "
                    "(per-block (m-1)!(k-1)!/(m+k-1)!). Z has no multinomial coefficient. "
                    "Gate receipts: gates/GATE_SWAP_ROUTE.json.\n")
        f.write("\n" + "\n".join(lines) + "\n")
    print("recorded row lines ->", ROWS_PATH)


# ---------------- float engine ----------------
def build_base(N, kind):
    u1 = N // 2
    u2 = N - u1
    i = np.arange(u1 + 1)[:, None]
    j = np.arange(N + 1)[None, :]
    mask = (i <= j) & ((j - i) <= u2)
    w = np.zeros(N + 1)
    if kind == "g":
        w = 1.0 / (np.arange(N + 1) + 1.0)
    else:
        jj = np.arange(1, N + 1, dtype=float)
        w[1:] = 1.0 / (jj * (jj + 1.0))
    A = mask * w[None, :]
    return A, mask.astype(float), u1, u2


def solve_tilt_g(A, N, g):
    """theta with single-factor tilted mean t-degree = N/g."""
    W = A.sum(axis=0)
    j = np.arange(N + 1, dtype=float)
    ok = W > 0
    lnW = np.full(N + 1, -np.inf)
    lnW[ok] = np.log(W[ok])
    target = N / float(g)
    lo, hi = -60.0, 5.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        a = lnW + mid * j
        M = a.max()
        e = np.exp(a - M)
        mean = float((j * e).sum() / e.sum())
        if mean < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def solve_tilt_dpm(A, N, alpha):
    """theta with alpha * sum_j j W_j e^{theta j} = N (total mean t-degree)."""
    W = A.sum(axis=0)
    j = np.arange(N + 1, dtype=float)
    ok = W > 0
    lnW = np.full(N + 1, -np.inf)
    lnW[ok] = np.log(W[ok])
    lnT = math.log(N) - math.log(alpha)
    lo, hi = -60.0, 5.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        a = lnW + mid * j
        M = a.max()
        s = float((j * np.exp(a - M)).sum())
        val = M + math.log(s) if s > 0 else -np.inf
        if val < lnT:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


class Conv2:
    def __init__(self, u1, N, mask):
        self.u1, self.N, self.mask = u1, N, mask
        self.P = next_fast_len(2 * u1 + 1)
        self.Q = next_fast_len(2 * N + 1)

    def fwd(self, a):
        return np.fft.rfft2(a, s=(self.P, self.Q))

    def back(self, fa):
        c = np.fft.irfft2(fa, s=(self.P, self.Q))[: self.u1 + 1, : self.N + 1]
        c *= self.mask
        np.maximum(c, 0.0, out=c)
        return c


def poly_pow_corner(B, g, cv, planted=False):
    """ln coeff of s^u1 t^N in B^g (B >= 0). Per-squaring max-normalization;
    planted=True drops the log-scale accumulation of the LAST squaring (the
    first squarings of a heavily damped base have max exactly 1, so a plant
    there is vacuous; the last one is load-bearing)."""
    m0 = float(B.max())
    b = B / m0
    lb = math.log(m0)
    R, lR = None, 0.0
    nsq = g.bit_length() - 1
    isq = 0
    while g:
        if g & 1:
            if R is None:
                R, lR = b.copy(), lb
            else:
                # deferred-spectrum R-update: never hold two spectra past the
                # multiply, free everything before the inverse (peak ~3 grids)
                fR = cv.fwd(R)
                del R
                Fb_t = cv.fwd(b)
                fR *= Fb_t
                del Fb_t
                R = cv.back(fR)
                del fR
                m = float(R.max())
                R /= m
                lR += lb + math.log(m)
        g >>= 1
        if g:
            Fb = cv.fwd(b)
            del b
            np.multiply(Fb, Fb, out=Fb)
            b = cv.back(Fb)
            del Fb
            m = float(b.max())
            b /= m
            isq += 1
            if planted and isq == nsq:
                lb = 2 * lb          # PLANTED ERROR: ln(m) dropped once
            else:
                lb = 2 * lb + math.log(m)
    return math.log(float(R[cv.u1, cv.N])) + lR


def lnZ_g_float(N, g, planted=False):
    A, mask, u1, u2 = build_base(N, "g")
    theta = solve_tilt_g(A, N, g)
    j = np.arange(N + 1, dtype=float)
    B = A * np.exp(theta * j)[None, :]
    cv = Conv2(u1, N, mask)
    lncorner = poly_pow_corner(B, g, cv, planted=planted) - theta * N
    lnpref = (math.lgamma(u1 + 1) + math.lgamma(u2 + 1)
              - (math.lgamma(g + N) - math.lgamma(g)))
    return lnpref + lncorner, theta


def lnZ_dpm_float(N, alpha, T=22):
    A, mask, u1, u2 = build_base(N, "dpm")
    theta = solve_tilt_dpm(A, N, alpha)
    j = np.arange(N + 1, dtype=float)
    X = alpha * A * np.exp(theta * j)[None, :]
    S = float(X.sum())
    d = max(0, math.ceil(math.log2(S / 0.5))) if S > 0.5 else 0
    X /= 2.0 ** d
    cv = Conv2(u1, N, mask)
    Xf = cv.fwd(X)
    E = X / T
    E[0, 0] += 1.0
    for p in range(T - 1, 0, -1):     # Horner: E <- I + (X/p) conv E
        fE = cv.fwd(E)
        del E
        fE *= Xf
        E = cv.back(fE)
        del fE
        E /= p
        E[0, 0] += 1.0
    ls = 0.0
    for _ in range(d):
        f = cv.fwd(E)
        del E
        np.multiply(f, f, out=f)
        E = cv.back(f)
        del f
        m = float(E.max())
        E /= m
        ls = 2 * ls + math.log(m)
    lncorner = math.log(float(E[u1, N])) + ls - theta * N
    lnpref = (math.lgamma(u1 + 1) + math.lgamma(u2 + 1)
              - (math.lgamma(alpha + N) - math.lgamma(alpha)))
    return lnpref + lncorner, theta, d


# ---------------- commands ----------------
def cmd_selftest():
    from w4_blind_gf import Zg_gf, Zdpm_gf
    out = []
    for N, g in ((8, 3), (12, 7), (10, 64)):
        u = (N // 2, N - N // 2)
        lf, _ = lnZ_g_float(N, g)
        le = lnfr(Zg_gf(u, g, g))
        out.append(("g", N, g, lf, le, abs(lf - le) / abs(le)))
    from fractions import Fraction as Fr
    for N, al in ((8, 1), (10, 2)):
        u = (N // 2, N - N // 2)
        lf, _, _ = lnZ_dpm_float(N, al)
        le = lnfr(Zdpm_gf(u, Fr(al)))
        out.append(("dpm", N, al, lf, le, abs(lf - le) / abs(le)))
    for r in out:
        print(f"{r[0]} N={r[1]} p={r[2]}: float={r[3]:.12f} exact={r[4]:.12f} rel={r[5]:.2e}")
    print("SELFTEST", "PASS" if all(r[5] < 1e-10 for r in out) else "FAIL")


def cmd_gate_g(N, g, label):
    from w4_blind_gf import Zg_gf
    u = (N // 2, N - N // 2)
    t0 = time.perf_counter()
    lf, theta = lnZ_g_float(N, g)
    tf = time.perf_counter() - t0
    t0 = time.perf_counter()
    le = lnfr(Zg_gf(u, g, g))
    te = time.perf_counter() - t0
    rel = abs(lf - le) / abs(le)
    entry = {"label": label, "kind": "finite_g_dir1", "N": N, "U": list(u), "g": g,
             "lnZ_float": lf, "lnZ_exact": le, "abs_diff": abs(lf - le),
             "rel_err": rel, "rel_err_str": f"{rel:.1e}", "theta_t": theta,
             "wall_float_s": round(tf, 3), "wall_exact_s": round(te, 1),
             "exact_ref": "Zg_gf(U, alpha=g, g) [Dir(1) identity]",
             "pass": bool(rel < 1e-10)}
    bank_gate(entry)
    print(json.dumps(entry, indent=1))


def cmd_gate_dpm(N, alpha, label):
    from fractions import Fraction as Fr
    from w4_blind_gf import Zdpm_gf
    u = (N // 2, N - N // 2)
    t0 = time.perf_counter()
    lf, theta, d = lnZ_dpm_float(N, alpha)
    tf = time.perf_counter() - t0
    t0 = time.perf_counter()
    le = lnfr(Zdpm_gf(u, Fr(alpha)))
    te = time.perf_counter() - t0
    rel = abs(lf - le) / abs(le)
    entry = {"label": label, "kind": "dpm_eppf", "N": N, "U": list(u), "alpha": alpha,
             "lnZ_float": lf, "lnZ_exact": le, "abs_diff": abs(lf - le),
             "rel_err": rel, "rel_err_str": f"{rel:.1e}", "theta_t": theta,
             "exp_scaling_d": d, "wall_float_s": round(tf, 3),
             "wall_exact_s": round(te, 1), "exact_ref": "Zdpm_gf(U, alpha)",
             "pass": bool(rel < 1e-10)}
    bank_gate(entry)
    print(json.dumps(entry, indent=1))


def cmd_gate_planted(N, g, ref_label):
    with open(GATE_PATH) as f:
        doc = json.load(f)
    ref = [e for e in doc["entries"] if e["label"] == ref_label][-1]
    le = ref["lnZ_exact"]
    lp, _ = lnZ_g_float(N, g, planted=True)
    rel = abs(lp - le) / abs(le)
    entry = {"label": f"planted_{ref_label}_final", "kind": "planted_error",
             "N": N, "g": g,
             "description": "skip one renormalization (log-scale of last squaring dropped)",
             "lnZ_planted": lp, "lnZ_exact": le, "rel_err_planted": rel,
             "clean_rel_err": ref["rel_err"],
             "caught": bool(rel > 1e-10),
             "pass": bool(rel > 1e-10)}
    bank_gate(entry)
    print(json.dumps(entry, indent=1))


def cmd_row_g(N, g, glabel):
    t0 = time.perf_counter()
    lf, theta = lnZ_g_float(N, g)
    wall = time.perf_counter() - t0
    ws = f"{wall:.1f} s" if wall < 100 else f"{wall:.0f} s"
    lines = [f"## row g={glabel} k=2 N={N}",
             f"- balanced U = ({N // 2}, {N - N // 2}); Dir(1,..,1) mixing "
             "(weight per column-sum m: m!(k-1)!/(m+k-1)!); no multinomial coefficient",
             f"- ln Z = {lf:.6f}",
             f"- wall {ws} (tilt theta_t = {theta:.6f}, binary powering, "
             "per-squaring max-normalization, rfft2 convolution)",
             "- verification: floating point with exponential tilting; same engine "
             "gated against exact-rational evaluation in gates/GATE_SWAP_ROUTE.json"]
    bank_row(lines)
    print(f"ROW g={glabel} N={N}: lnZ={lf:.6f} wall={ws}")


def cmd_row_dpm(N, alpha):
    t0 = time.perf_counter()
    lf, theta, d = lnZ_dpm_float(N, alpha)
    wall = time.perf_counter() - t0
    ws = f"{wall:.1f} s" if wall < 100 else f"{wall:.0f} s"
    lines = [f"## row DPM alpha={alpha} k=2 N={N}",
             f"- balanced U = ({N // 2}, {N - N // 2}); EPPF register "
             "(per-block (m-1)!(k-1)!/(m+k-1)!); no multinomial coefficient",
             f"- ln Z = {lf:.6f}",
             f"- wall {ws} (tilt theta_t = {theta:.6f}, series exp: scaled Taylor "
             f"T=22 then {d} squarings, rfft2 convolution)",
             "- verification: floating point with exponential tilting; same engine "
             "gated against exact-rational evaluation in gates/GATE_SWAP_ROUTE.json"]
    bank_row(lines)
    print(f"ROW dpm alpha={alpha} N={N}: lnZ={lf:.6f} wall={ws}")


def cmd_probe():
    walls = {}
    for N in (5000, 10000):
        A, mask, u1, u2 = build_base(N, "g")
        cv = Conv2(u1, N, mask)
        t0 = time.perf_counter()
        f = cv.fwd(A)
        _ = cv.back(f * f)
        walls[N] = time.perf_counter() - t0
        print(f"probe N={N}: padded ({cv.P},{cv.Q}) squaring wall {walls[N]:.2f} s")
    u1_50k = 25000
    P50, Q50 = next_fast_len(2 * u1_50k + 1), next_fast_len(2 * 50000 + 1)
    pts10 = next_fast_len(10001) * next_fast_len(20001)
    pts50 = P50 * Q50
    scale = (pts50 / pts10) * (math.log(pts50) / math.log(pts10))
    sq50 = walls[10000] * scale
    mem_real = pts50 * 8 / 2 ** 30
    mem_cplx = P50 * (Q50 // 2 + 1) * 16 / 2 ** 30
    # transform counts: g=1e6 -> 19 squarings + 7 result-multiplies (popcount);
    # DPM alpha=1 -> T=22 Horner multiplies + ~5 squarings
    proj_g = sq50 * (19 + 7 * 1.5)
    proj_dpm = sq50 * (22 + 5)
    lines = ["## projection receipt: N=50000 cells (measured scaling basis)",
             f"- measured squaring wall: {walls[5000]:.2f} s at N=5000, "
             f"{walls[10000]:.2f} s at N=10000 (single core)",
             f"- N=50000 padded FFT grid ({P50},{Q50}): {mem_real:.1f} GiB real + "
             f"{mem_cplx:.1f} GiB complex per transform -> exceeds the 31 GiB "
             "address-space cap",
             f"- projected squaring wall at N=50000: {sq50:.0f} s; "
             f"g=10^6 needs ~26 products -> ~{proj_g / 60:.0f} min; "
             f"DPM alpha=1 needs ~27 products -> ~{proj_dpm / 60:.0f} min; "
             "both far beyond the 10-min single-cell limit",
             "- VERDICT: rows (g=10^6, N=50000) and (DPM alpha=1, N=50000) SKIPPED "
             "(projected > 10 min AND over the memory cap); receipts recorded here "
             "and in gates/GATE_SWAP_ROUTE.json"]
    bank_row(lines)
    bank_gate({"label": "skip_projection_N50000", "kind": "projection",
               "sq_wall_5000_s": round(walls[5000], 2),
               "sq_wall_10000_s": round(walls[10000], 2),
               "padded_grid_N50000": [P50, Q50],
               "mem_real_GiB": round(mem_real, 1), "mem_cplx_GiB": round(mem_cplx, 1),
               "proj_squaring_s": round(sq50), "proj_g1e6_min": round(proj_g / 60),
               "proj_dpm_min": round(proj_dpm / 60),
               "verdict": "SKIP both N=50000 rows (>10 min projected; >31 GiB)"})


def main():
    cmd = sys.argv[1]
    if cmd == "selftest":
        cmd_selftest()
    elif cmd == "gate_g":
        cmd_gate_g(int(sys.argv[2]), int(sys.argv[3]), sys.argv[4])
    elif cmd == "gate_dpm":
        cmd_gate_dpm(int(sys.argv[2]), int(sys.argv[3]), sys.argv[4])
    elif cmd == "gate_planted":
        cmd_gate_planted(int(sys.argv[2]), int(sys.argv[3]), sys.argv[4])
    elif cmd == "row_g":
        cmd_row_g(int(sys.argv[2]), int(sys.argv[3]), sys.argv[4])
    elif cmd == "row_dpm":
        cmd_row_dpm(int(sys.argv[2]), int(sys.argv[3]))
    elif cmd == "probe":
        cmd_probe()
    else:
        raise SystemExit(f"unknown cmd {cmd}")


if __name__ == "__main__":
    main()
