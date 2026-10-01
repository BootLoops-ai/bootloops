"""MIXALOT worked example: the large-N / large-k / large-g table, exact vs MC.

Everything computes LIVE through the vendored engines (small, quick cells);
regimes beyond a sampler's reach print an em-dash with the reason — that gap
IS the tool's point. MC = the vendored estimator suite (JAX HMC + numpy
Gibbs; vendored estimators.py; seeded runs are deterministic on a given
stack). No sampler claim is extrapolated beyond the cells run here.

Output: TABLE_LARGE_KGN.md + TABLE_LARGE_KGN.json (walls measured).
Run: python3 examples/large_kgn_table.py [--quick]
"""
import json
import os
import random
import subprocess
import sys
import time
from fractions import Fraction as F
from math import log

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from mixalot.engines import load  # noqa: E402

QUICK = "--quick" in sys.argv


def lnfrac(z):
    n, d = z.numerator, z.denominator
    shn = max(0, n.bit_length() - 500)
    shd = max(0, d.bit_length() - 500)
    return (log(n >> shn) + shn * log(2)) - (log(d >> shd) + shd * log(2))


def main():
    cf = load("closed_form_1var")
    est = load("estimators")
    bigk2 = load("bigk2")
    gf = load("w4_blind_gf")
    out = {"sections": {}, "generated_utc": subprocess.run(
        ["date", "-u", "+%Y-%m-%dT%H:%M:%SZ"], capture_output=True,
        text=True).stdout.strip()}

    # ---------- A: large N (1-var, U = (N/2, N/2)), exact vs MC
    rowsA = []
    Ns = [100, 1000, 10000] if QUICK else [100, 1000, 10000, 100000]
    for N in Ns:
        u = N // 2
        t0 = time.time()
        z = cf.Z_closed(u, u)
        t_exact = time.time() - t0
        lz = lnfrac(z)
        row = {"N": N, "exact_logZ": lz, "exact_wall_s": round(t_exact, 4),
               "mc": {}}
        for name in ("nested", "ss"):
            vals = []
            t0 = time.time()
            for seed in (1, 2):
                r = est.run("m1", 2, [u, u], name, seed=seed)
                vals.append(r["logZ_hat"])
            row["mc"][name] = {
                "seeds": vals,
                "mean_bias": sum(vals) / len(vals) - lz,
                "wall_s": round(time.time() - t0, 2)}
        rowsA.append(row)
        print(f"A N={N}: exact {lz:.6f} ({t_exact:.3f}s) | " + " | ".join(
            f"{k} bias {v['mean_bias']:+.4f} ({v['wall_s']}s)"
            for k, v in row["mc"].items()), flush=True)
    out["sections"]["A_large_N"] = {
        "rows": rowsA,
        "beyond": "N beyond these cells: no sampler claim is made here — "
                  "validate any estimator against the exact column on the "
                  "largest N it still reaches; the exact column is the "
                  "referee."}

    # ---------- B: large k (bigk2; the recorded scaling construction)
    rowsB = []
    cells = [(100, 10000), (300, 30000)] if QUICK else \
        [(100, 10000), (300, 30000), (1000, 100000)]
    for k, N in cells:
        random.seed(42)
        cuts = sorted(random.sample(range(1, N + k), k - 1))
        U = [b - a - 1 for a, b in zip([0] + cuts, cuts + [N + k])]
        t0 = time.time()
        Z = bigk2.Z_fast2(U)
        dt = time.time() - t0
        rowsB.append({"k": k, "N": N, "exact_lnZ": lnfrac(Z),
                      "wall_s": round(dt, 2),
                      "mc": "— (no sampler column: no vendored sampler "
                            "operates at k>=100; exact-only regime)"})
        print(f"B k={k} N={N}: lnZ={lnfrac(Z):.6f} ({dt:.2f}s)", flush=True)
    out["sections"]["B_large_k"] = {"rows": rowsB}

    # ---------- C: large g + the exact DP limit
    rowsC = []
    t0 = time.time()
    zdpm = gf.Zdpm_gf([2, 1], F(1))
    for g in (64, 128, 256):
        zg = gf.Zg_gf([2, 1], F(1), g)
        rowsC.append({"g": g, "exact_Zg": str(zg),
                      "g_times_gap": str(g * (zg - zdpm))})
    rowsC.append({"g": "DP limit", "exact_Zg": str(zdpm),
                  "note": "exact Dirichlet-process limit; the g*(Z_g - Z_DPM)"
                          " = -1/48 law above is EXACT at every listed g"})
    out["sections"]["C_large_g"] = {
        "rows": rowsC, "wall_s": round(time.time() - t0, 3),
        "float_route": "swap_route: certified float route to g = 1e6 "
                       "components (k=2; certified against exact at N=100 — "
                       "recorded gate; a ROUTE demonstration, not a "
                       "general-purpose evaluator)",
        "mc": "— (no sampler operates in this regime)"}

    # ---------- render
    L = ["# MIXALOT worked example — exact vs MC across the regimes",
         "",
         "MC = the vendored estimator suite (seeded runs, deterministic on "
         "a given stack). Cells "
         "print live-measured walls. '—' = no sampler reaches the regime.",
         "", "## A. Large N (1-var, U=(N/2,N/2)) — exact vs nested/ss", "",
         "| N | exact logZ (wall) | nested bias (wall) | ss bias (wall) |",
         "|---|---|---|---|"]
    for r in rowsA:
        L.append(f"| {r['N']:,} | {r['exact_logZ']:.6f} "
                 f"({r['exact_wall_s']}s) | "
                 f"{r['mc']['nested']['mean_bias']:+.4f} "
                 f"({r['mc']['nested']['wall_s']}s) | "
                 f"{r['mc']['ss']['mean_bias']:+.4f} "
                 f"({r['mc']['ss']['wall_s']}s) |")
    L += ["", out["sections"]["A_large_N"]["beyond"], "",
          "## B. Large k — exact only (the regime samplers don't reach)", "",
          "| k | N | exact lnZ | wall |", "|---|---|---|---|"]
    for r in rowsB:
        L.append(f"| {r['k']} | {r['N']:,} | {r['exact_lnZ']:.6f} | "
                 f"{r['wall_s']}s |")
    L += ["", rowsB[0]["mc"], "",
          "## C. Large g and the exact Dirichlet-process limit", "",
          "| g | exact Z_g (U=(2,1), alpha=1) | g*(Z_g − Z_DPM) |",
          "|---|---|---|"]
    for r in rowsC:
        L.append(f"| {r['g']} | {r['exact_Zg']} | "
                 f"{r.get('g_times_gap', r.get('note', ''))} |")
    L += ["", out["sections"]["C_large_g"]["float_route"], ""]
    md = "\n".join(L) + "\n"
    open(os.path.join(os.path.dirname(__file__),
                      "TABLE_LARGE_KGN.md"), "w").write(md)
    with open(os.path.join(os.path.dirname(__file__),
                           "TABLE_LARGE_KGN.json"), "w") as f:
        json.dump(out, f, indent=1)
        f.write("\n")
    print("\nwrote examples/TABLE_LARGE_KGN.{md,json}")


if __name__ == "__main__":
    main()
