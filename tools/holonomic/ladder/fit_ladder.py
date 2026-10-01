#!/usr/bin/env python3
"""fit_ladder.py — least-squares power-law fit of certified-transport cost vs
holonomic rank from the ladder_R*.json receipts written by rank_ladder.py.
Fits log t = a + b log R per (path, eps) stage over the ranks where the stage
was run; extrapolates to the ranks you name.

    python3 fit_ladder.py [--dir RECEIPT_DIR] [--min-rank 24]
                          [--extrapolate 96,120,720]
"""
import argparse
import glob
import json
import math
import os

ap = argparse.ArgumentParser()
ap.add_argument("--dir", default=os.path.dirname(os.path.abspath(__file__)),
                help="directory holding ladder_R*.json receipts")
ap.add_argument("--min-rank", type=int, default=24,
                help="drop warm-up ranks below this from the fit")
ap.add_argument("--extrapolate", default="96,120,720",
                help="comma list of ranks to price from the fit")
args = ap.parse_args()

data = {}
for fn in sorted(glob.glob(f"{args.dir}/ladder_R*.json")):
    d = json.load(open(fn))
    R = int(d["rank"])
    for st, s in d.get("stages", {}).items():
        data.setdefault(st, []).append(
            (R, float(s["wall_s"]), float(s["max_entry_rad"]),
             float(s["maxrss_mb_after"])))
if not data:
    raise SystemExit(f"no ladder_R*.json receipts in {args.dir} — "
                     f"run rank_ladder.py first (see its header)")


def fit(pairs):
    xs = [math.log(r) for r, t, *_ in pairs]
    ys = [math.log(t) for r, t, *_ in pairs]
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) \
        / sum((x - mx) ** 2 for x in xs)
    a = my - b * mx
    resid = max(abs(y - (a + b * x)) for x, y in zip(xs, ys))
    return a, b, resid


def human(sec):
    if sec < 120:
        return f"{sec:.1f} s"
    if sec < 7200:
        return f"{sec/60:.1f} min"
    if sec < 48 * 3600:
        return f"{sec/3600:.1f} h"
    if sec < 2 * 365.25 * 86400:
        return f"{sec/86400:.1f} d"
    return f"{sec/(365.25*86400):.2f} yr"


print(f"power-law fits over ranks >= {args.min_rank} (t = C * R^b):")
fits = {}
for st in sorted(data):
    pairs = [p for p in sorted(data[st]) if p[0] >= args.min_rank]
    if len(pairs) < 2:
        print(f"  {st}: <2 points, skipped")
        continue
    a, b, resid = fit(pairs)
    fits[st] = (a, b)
    pts = ", ".join(f"R{r}:{t:g}s" for r, t, *_ in pairs)
    print(f"  {st}: exponent b = {b:.2f}  (C = {math.exp(a):.3g}; "
          f"max log-resid {resid:.2f})  [{pts}]")

targets = [int(t) for t in args.extrapolate.split(",") if t.strip()]
if targets and fits:
    print("\nextrapolation (single process, same machine, naive path):")
    for st, (a, b) in sorted(fits.items()):
        row = "  ".join(f"R={R}: {human(math.exp(a + b*math.log(R)))}"
                        for R in targets)
        print(f"  {st}:  {row}")
