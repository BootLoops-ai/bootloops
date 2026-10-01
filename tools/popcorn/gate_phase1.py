#!/usr/bin/env python3
"""Held-out two-route gate: recurrence vs direct, >=40 matched digits at
held-out points across the S range. Writes GATE_PHASE1.json into the cwd —
run it from a scratch directory.

Route 1 (recurrence): exact-Q tridiagonal engine (n<=1024) / arb stable-
  direction recursion (n>1024) — sfs_engine.py + arb_route.py.
Route 2 (direct): per-i mpmath hyp1f1 — independent of the tridiagonal
  system and of arb; re-anchored here at 3 points against an exact-rational
  certified positive-term series.

Held-out points: fixed fresh seed 77190708. Coverage: the moderate-|S| leak
band, near-neutral, hostile |S|=1000-2000, n up to 10007. Vector-level
(ALL i) for n<=1024; 12 spot i for n=10007. Gradient gate via the
parameter-shift identity at 3 i per point. Bar: >=40 matched digits.
"""
import json, os, random, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction
from mpmath import mp, mpf, exp, hyp1f1, fabs, log10
from sfs_engine import solve_M_exact, solve_dM_exact, assemble_M
from arb_route import solve_vector_arb

TARGET = 40
CMP_DPS = 60   # compare at 60d so a 40d bar has headroom


def digits(a, b):
    if a == b:
        return 9999.0
    if b == 0:
        return 9999.0 if a == 0 else 0.0
    return float(-log10(fabs((a - b) / b)))


def series_anchor(n, i, S):
    """Exact-rational certified positive-series 1F1 anchor."""
    zfr = Fraction(-S) if S < 0 else Fraction(S)
    a, b = (n - i, n) if S < 0 else (i, n)
    t = Fraction(1); s = Fraction(1); k = 0
    tol = Fraction(1, 10 ** (CMP_DPS + 30))
    while True:
        t *= Fraction(zfr.numerator * (a + k), zfr.denominator * (b + k) * (k + 1))
        s += t
        k += 1
        r = Fraction(zfr.numerator * (a + k), zfr.denominator * (b + k) * (k + 1))
        if r < Fraction(1, 2) and t < s * tol:
            break
    with mp.workdps(CMP_DPS + 25):
        R = mpf(s.numerator) / mpf(s.denominator)
        return R if S < 0 else exp(-mpf(S)) * R  # = 1F1(n-i;n;-S)


def gate_point(n, S):
    rec = {"n": n, "S": S}
    t0 = time.time()
    # ONE exact rational S threaded through BOTH routes (decimal-string exact;
    # a binary-float/rational mismatch is amplified by the cancellation factor
    # — measured -7000d on the first gate run)
    Sfr = Fraction(str(S))
    if n <= 1024:
        alpha, beta = solve_M_exact(n, Sfr)
        M1 = assemble_M(alpha, beta, Sfr, CMP_DPS + 20)
        idxs = list(range(1, n))
        ders = solve_dM_exact(n, Sfr, alpha, beta, 1)
        a1, b1 = ders[0]
        G1 = assemble_M(a1, b1, Sfr, CMP_DPS + 20)
        route1 = "exact-tridiag"
    else:
        sol = solve_vector_arb(n, S, CMP_DPS + 10, want_grad=True, S_frac=Sfr)
        M1 = [None] * (n + 1)
        G1 = [None] * (n + 1)
        idxs = sorted(set([1, 2, 3, n // 7, n // 3, n // 2, (2 * n) // 3,
                           (9 * n) // 10, n - 3, n - 2, n - 1,
                           max(1, n // 1000)]))
        with mp.workdps(CMP_DPS + 20):
            for i in idxs:
                M1[i] = mp.mpmathify(sol["M"][i].mid().str(CMP_DPS + 15, radius=False))
                G1[i] = mp.mpmathify(sol["Mp"][i].mid().str(CMP_DPS + 15, radius=False))
        route1 = f"arb-directional(prec={sol['prec']})"
    # route 2: per-i direct hyp1f1 at the SAME exact rational S
    worstM = 9999.0
    worstG = 9999.0
    with mp.workdps(CMP_DPS + 20):
        Sm = mpf(Sfr.numerator) / Sfr.denominator
        for i in idxs:
            truth = hyp1f1(n - i, n, -Sm)
            if truth != 0:
                worstM = min(worstM, digits(M1[i], truth))
        for i in (idxs[0], idxs[len(idxs) // 2], idxs[-1]):
            truthG = -(mpf(n - i) / n) * hyp1f1(n - i + 1, n + 1, -Sm)
            if truthG != 0:
                worstG = min(worstG, digits(G1[i], truthG))
    rec.update({"route1": route1, "n_indices": len(idxs),
                "worst_M_digits": round(worstM, 1),
                "worst_grad_digits": round(worstG, 1),
                "wall_s": round(time.time() - t0, 2)})
    return rec


def main():
    rng = random.Random(77190708)
    pts = []
    for n in (23, 137, 1024):
        pts += [(n, round(rng.choice([-1, 1]) * 10 ** rng.uniform(-2, 3.3), 4))
                for _ in range(4)]
        pts += [(n, round(rng.choice([-1, 1]) * 10 ** rng.uniform(0.5, 1.5), 4))]  # leak band
        pts += [(n, rng.choice([-1, 1]) * 1e-4)]                                   # near-neutral
    pts += [(10007, -1000.0), (10007, 913.0), (10007, -7.3), (10007, 2000.0)]
    recs = [gate_point(n, S) for (n, S) in pts]
    # re-anchor route 2 against the exact-rational series at 3 points
    anchors = []
    for (n, i, S) in [(137, 68, -1000.0), (1024, 512, 555.0), (23, 1, -0.0001)]:
        with mp.workdps(CMP_DPS + 20):
            d = digits(hyp1f1(n - i, n, -mpf(S)), series_anchor(n, i, S))
        anchors.append({"n": n, "i": i, "S": S, "hyp1f1_vs_series_digits": round(d, 1)})
    worst = min(r["worst_M_digits"] for r in recs)
    worstg = min(r["worst_grad_digits"] for r in recs)
    worsta = min(a["hyp1f1_vs_series_digits"] for a in anchors)
    verdict = {"worst_M_digits": worst, "worst_grad_digits": worstg,
               "worst_anchor_digits": worsta, "bar": TARGET,
               "n_points": len(recs),
               "PASS": bool(worst >= TARGET and worstg >= TARGET and worsta >= 55)}
    out = {"verdict": verdict, "points": recs, "anchors": anchors}
    with open("GATE_PHASE1.json", "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(verdict, indent=1))
    for r in recs:
        print(f"n={r['n']:6d} S={r['S']:+10.4f} [{r['route1']:>28s}] "
              f"M {r['worst_M_digits']:6.1f}d grad {r['worst_grad_digits']:6.1f}d "
              f"({r['n_indices']} idx, {r['wall_s']}s)")
    for a in anchors:
        print(f"anchor n={a['n']} i={a['i']} S={a['S']}: hyp1f1 vs exact series {a['hyp1f1_vs_series_digits']}d")
    return 0 if verdict["PASS"] else 1


if __name__ == "__main__":
    sys.exit(main())
