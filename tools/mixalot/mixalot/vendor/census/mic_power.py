#!/usr/bin/env python3
"""Micah synthetic-truth power probe.
Planted-contrast construction (fixed profiles, Micah per-chapter
closed-class token sizes, cut after ch4, fixed seeds and sampling loop),
scored with the validated seg_v1 float engine (uniform prior over g=1..6,
alpha=1/2). Detection criterion: P(g=1) < 0.5.

Arms:
  A. Genesis-strength closed-class contrast at Micah size (pooled Gen
     chs 1-11 vs 12-50 profiles), 30 trials, seeds 1000+i  — gate: 0/30.
  B. Half-strength variant (process B = 50/50 profile mix), 30 trials,
     seeds 2000+i — gate: 0/30.
  C. Language-swap contrast at Micah size (pooled Dan Hebrew chs
     1,8-12 vs Aramaic chs 3-7 closed profiles; the Daniel-seam
     mechanism), 10 trials, seeds 3000+i — gate: 10/10.

Output: raw/mic_power.json (numbers only; consumed by assemble_sweep).
"""
import json, os, random, sys
from fractions import Fraction
from math import log, exp, comb

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# public port: pilot engine dir env-pointed (no machine-local default)
PILOT = os.environ.get("MIXALOT_PILOT_DIR", "")
if not PILOT:
    raise RuntimeError("MIXALOT_PILOT_DIR unset — path to the pilot engine "
                       "dir (seg_v1.py) is required")
sys.path.insert(0, PILOT)
sys.path.insert(0, f"{BASE}/code")
import seg_v1 as sv
import production_sweep as pr

G = 6
ALPHA = Fraction(1, 2)


def p_g1(X):
    """P(g=1 | X), uniform prior over g=1..6, seg_v1 float engine."""
    LW = sv._logW_table(X, ALPHA)
    n = len(X)
    lz = [sv.dp_float(X, g, ALPHA, LW) - log(comb(n - 1, g - 1))
          for g in range(1, G + 1)]
    m = max(lz)
    w = [exp(v - m) for v in lz]
    return w[0] / sum(w)


def sample_units(sizes, cut, wA, wB, k, seed):
    """Gate's sampling loop verbatim: sequential acc-scan per token."""
    rng = random.Random(seed)
    units = []
    for i, s in enumerate(sizes):
        w = wA if i < cut else wB
        row = [0] * k
        for _ in range(s):
            r = rng.random()
            acc = 0.0
            for c, wc in enumerate(w):
                acc += wc
                if r <= acc:
                    row[c] += 1
                    break
            else:
                row[-1] += 1
        units.append(row)
    return units


def arm(name, sizes, wA, wB, k, trials, seed0):
    det = 0
    for i in range(trials):
        X = sample_units(sizes, 4, wA, wB, k, seed0 + i)
        if p_g1(X) < 0.5:
            det += 1
    print(f"POWER {name}: detected (P(g=1)<0.5) in {det}/{trials} trials",
          flush=True)
    return {"detected": det, "trials": trials, "seed0": seed0}


def pooled(X, rows):
    return [sum(X[i][c] for i in rows) for c in range(len(X[0]))]


def main():
    Xm, _, _ = pr.chapter_matrix("Mic", "closed")
    sizes = [sum(u) for u in Xm]
    print("Micah chapter closed token sizes:", sizes, "total", sum(sizes))

    # arms A, B: Genesis-strength profiles (gate construction)
    Xg, gv, _ = pr.chapter_matrix("Gen", "closed")
    kg = len(gv)
    pA = pooled(Xg, range(0, 11))
    pB = pooled(Xg, range(11, 50))
    wA = [x / sum(pA) for x in pA]
    wB = [x / sum(pB) for x in pB]
    gen = arm("Genesis-strength contrast at Micah size",
              sizes, wA, wB, kg, 30, 1000)
    wB2 = [(x + y) / 2 for x, y in zip(wA, wB)]
    half = arm("half-strength contrast at Micah size",
               sizes, wA, wB2, kg, 30, 2000)

    # arm C: language-swap profiles (Daniel seam; ch2 mixed, excluded)
    Xd, dv, chd = pr.chapter_matrix("Dan", "closed")
    kd = len(dv)
    heb = pooled(Xd, [chd.index(f"ch{c}") for c in (1, 8, 9, 10, 11, 12)])
    ara = pooled(Xd, [chd.index(f"ch{c}") for c in (3, 4, 5, 6, 7)])
    wH = [x / sum(heb) for x in heb]
    wR = [x / sum(ara) for x in ara]
    swap = arm("language-swap contrast at Micah size",
               sizes, wH, wR, kd, 10, 3000)

    out = {"spec": "planted-contrast power probe at Micah's exact shape "
                   "(same per-chapter closed token counts, cut after ch4); "
                   "detection = P(g=1) < 0.5, uniform prior g=1..6, "
                   "alpha=1/2, seg_v1 float engine (gated)",
           "adopted_from": "the gate's power probe "
                           "(Genesis-strength 0/30 reproduced)",
           "mic_sizes": sizes, "total_closed_tokens": sum(sizes),
           "n_chapters": len(sizes),
           "genesis_strength": gen, "half_strength": half,
           "language_swap": swap}
    with open(f"{BASE}/raw/mic_power.json", "w") as f:
        json.dump(out, f, indent=1, sort_keys=True)
    print("wrote raw/mic_power.json")


if __name__ == "__main__":
    main()
