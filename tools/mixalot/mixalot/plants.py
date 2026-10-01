"""Synthetic-truth plants + blind recovery — the house discipline, in the box.

plant(): generate allocation counts from a KNOWN g-component mixture.
recover_blind(): run the exact gstar posterior on planted data and score
recovery. Every MIXALOT deployment should run its battery before trusting a
verdict on real data (Rule: plants before claims)."""
import random
from fractions import Fraction

from .core import gstar


def plant(k, g_true, N, seed, separation="strong"):
    """Counts U from a known g_true-component mixture of k-state categoricals.

    separation: 'strong' = components concentrated on disjoint state blocks;
    'weak' = components drawn iid uniform on the simplex (harder recovery,
    honest INDET expected at small N). Deterministic given seed."""
    rng = random.Random(seed)
    if k < g_true:
        raise ValueError("need k >= g_true for block separation")
    # component weights ~ uniform simplex (exchangeable)
    w = [rng.random() for _ in range(g_true)]
    tot = sum(w)
    w = [x / tot for x in w]
    # component-state profiles
    profiles = []
    if separation == "strong":
        blocks = [list(range(i, k, g_true)) for i in range(g_true)]
        for b in blocks:
            p = [0.02 / k] * k
            for s in b:
                p[s] = 1.0
            t = sum(p)
            profiles.append([x / t for x in p])
    else:
        for _ in range(g_true):
            p = [rng.random() for _ in range(k)]
            t = sum(p)
            profiles.append([x / t for x in p])
    U = [0] * k
    for _ in range(N):
        c = rng.choices(range(g_true), weights=w)[0]
        s = rng.choices(range(k), weights=profiles[c])[0]
        U[s] += 1
    return U


def recover_blind(k=4, N=60, gmax=3, seeds=range(5), separation="strong",
                  p_threshold=Fraction(7, 10)):
    """Plant-and-recover harness. Returns per-truth results; the caller
    decides pass/fail. MEASURED TRUTH at the defaults (k=4, N=60, strong
    separation): MAP hits 0/5, 0/5, 4/5 for g_true=1,2,3 and 0/15 confident —
    iid category-count g-inference is prior-geometry dominated at this size
    (see MANUAL, weak-identifiability fence). This harness exists so YOU
    measure the power at YOUR shape before trusting a verdict, not because
    recovery is expected at small N."""
    results = {}
    for g_true in range(1, gmax + 1):
        rows = []
        for seed in seeds:
            U = plant(k, g_true, N, seed=1000 * g_true + seed,
                      separation=separation)
            r = gstar(U, gmax=gmax)
            rows.append({"seed": seed, "U": U, "map_g": r["map_g"],
                         "hit": r["map_g"] == g_true,
                         "p_true": r["posterior"][g_true],
                         "confident": r["posterior"][g_true] >= p_threshold})
        results[g_true] = rows
    return results
