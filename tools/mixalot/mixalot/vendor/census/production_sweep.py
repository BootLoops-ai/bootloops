#!/usr/bin/env python3
"""Production grid for the first census sweep wave. Engine: seg_v1
float mode, gated by GATE_SEGMENT.json PASS (the Isaiah pilot); grid settings
IDENTICAL to the pilot production grid:

Grid: book in {Dan, Gen, Mic, Esth} x granularity in {chapter, block} x
channel in {closed, full} x alpha in {1/k, 1/2, 1}; g = 1..6; priors on g:
uniform and geometric(1/2), both declared. Pericope-block rule (mechanical,
identical to pilot): close a block once its cumulative FULL-segment count
reaches >= 400; remainder after the last close joins as the final block.

Esther rides as the wave's negative control under these identical settings
(its count files are byte-identical to the pilot's).

Output: raw/prod_grid.json (numbers only; no prose).
"""
import json, sys, time, os
from fractions import Fraction
from math import log, exp, comb

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# public port: pilot engine dir env-pointed (no machine-local default)
PILOT = os.environ.get("MIXALOT_PILOT_DIR", "")
if not PILOT:
    raise RuntimeError("MIXALOT_PILOT_DIR unset — path to the pilot engine "
                       "dir (seg_v1.py) is required")
sys.path.insert(0, PILOT)
import seg_v1 as sv

G, BLOCK_TARGET = 6, 400
BOOKS = ["Dan", "Gen", "Mic", "Esth"]


def chapter_matrix(book, channel):
    X, vocab, chs = sv.load_matrix(f"{BASE}/data/counts_{book}_A.json", channel)
    labels = [f"ch{c}" for c in chs]
    return X, vocab, labels


def block_matrix(book, channel):
    obj = json.load(open(f"{BASE}/data/verses_{book}.json"))
    verses = obj["verses"]
    keys = set()
    for v in verses:
        keys.update(v[channel])
    vocab = sorted(keys)
    idx = {v: i for i, v in enumerate(vocab)}
    blocks, cur, cur_full, last_verse = [], [0] * len(vocab), 0, None
    labels = []
    for v in verses:
        for kk, n in v[channel].items():
            cur[idx[kk]] += n
        cur_full += sum(v["full"].values())
        last_verse = v["verse"]
        if cur_full >= BLOCK_TARGET:
            blocks.append(cur)
            labels.append(last_verse)          # block ENDS after this verse
            cur, cur_full = [0] * len(vocab), 0
    if cur_full:
        blocks.append(cur)
        labels.append(last_verse)
    return blocks, vocab, labels


def posteriors(logZ, prior):
    lw = [logZ[g] + log(prior[g]) for g in range(1, G + 1)]
    m = max(lw)
    w = [exp(v - m) for v in lw]
    s = sum(w)
    return {g: w[g - 1] / s for g in range(1, G + 1)}


def run_config(X, a):
    out = {"logZ": {}, "boundary": {}}
    LW = sv._logW_table(X, a)
    n = len(X)
    for g in range(1, G + 1):
        out["logZ"][g] = sv.dp_float(X, g, a, LW) - log(comb(n - 1, g - 1))
    for g in range(2, G + 1):
        bp = sv.boundary_posterior_float(X, g, a)
        out["boundary"][g] = {t: round(p, 8) for t, p in bp.items()}
    return out


def main():
    uni = {g: 1.0 / G for g in range(1, G + 1)}
    geo = {g: 2.0 ** -g for g in range(1, G + 1)}
    gs = sum(geo.values())
    geo = {g: v / gs for g, v in geo.items()}
    grid = {"block_target_full_segments": BLOCK_TARGET, "G": G,
            "priors": {"uniform": uni, "geometric_half": geo}, "configs": {}}
    for book in BOOKS:
        for gran in ("chapter", "block"):
            for channel in ("closed", "full"):
                if gran == "chapter":
                    X, vocab, labels = chapter_matrix(book, channel)
                else:
                    X, vocab, labels = block_matrix(book, channel)
                k = len(vocab)
                N = sum(map(sum, X))
                for aname, a in (("1/k", Fraction(1, k)),
                                 ("1/2", Fraction(1, 2)), ("1", Fraction(1))):
                    t0 = time.time()
                    r = run_config(X, a)
                    key = f"{book}|{gran}|{channel}|{aname}"
                    pg_u = posteriors(r["logZ"], uni)
                    pg_g = posteriors(r["logZ"], geo)
                    gstar = max(pg_u, key=pg_u.get)
                    entry = {
                        "n": len(X), "k": k, "N": N,
                        "unit_labels": labels,
                        "logZ": {g: round(v, 6) for g, v in r["logZ"].items()},
                        "P_g_uniform": {g: round(v, 6) for g, v in pg_u.items()},
                        "P_g_geometric_half":
                            {g: round(v, 6) for g, v in pg_g.items()},
                        "g_star_uniform": gstar,
                        "boundary": r["boundary"],
                        "secs": round(time.time() - t0, 2)}
                    grid["configs"][key] = entry
                    print(key, "n=%d k=%d" % (len(X), k),
                          "P(g)=", {g: round(v, 3) for g, v in pg_u.items()},
                          "%.1fs" % entry["secs"], flush=True)
    with open(f"{BASE}/raw/prod_grid.json", "w") as f:
        json.dump(grid, f, sort_keys=True)
    print("wrote raw/prod_grid.json")


if __name__ == "__main__":
    main()
