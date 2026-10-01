#!/usr/bin/env python3
"""Truth-checked march endpoint case as a standalone gate over an arbitrary
march_lib file ({FILE} target for the mutation harness). rc=0 iff the
complex-step ball contains the 256-bit point-accurate e^h AND the real-h
control passes. Used by battery leg L5 (pristine PASS + mutants caught)."""
import importlib.util
import sys

from flint import acb, acb_poly, ctx


def main(path):
    spec = importlib.util.spec_from_file_location("march_target", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    DENP = acb_poly([0, 1]); NUMP = {(0, 0): acb_poly([0, 1])}
    old = ctx.prec
    try:
        ctx.prec = 256
        truths = {0.05: acb(0.05).exp(), (0.05, 0.03): acb(0.05, 0.03).exp()}
    finally:
        ctx.prec = old
    # K=30: truncation ~1e-40 (rounding-dominated); K=6: truncation ~5e-13 —
    # at K=6 the rigorous tail is LOAD-BEARING, so a broken/zeroed tail term
    # breaks containment there (the mutation-visibility case).
    for K in (30, 6):
        for h, truth in ((0.05, truths[0.05]), (acb(0.05, 0.03), truths[(0.05, 0.03)])):
            out = mod.step(NUMP, DENP, 1, acb(1), [acb(1)], h, 0.5, K)
            if not out[0].contains(truth):
                print(f"FAIL K={K} h={h}: ball {out[0]} does not contain exact {truth}")
                return 1
            # tail honesty: finite and not absurdly wide
            rad = float(out[0].real.rad())
            if not (rad < (1e-10 if K == 30 else 1e-4)):
                print(f"FAIL K={K} h={h}: radius {rad} unreasonably wide")
                return 1
    # contract refusals: |h|>=r, r<0, and wrong v0 length must all raise —
    # accepting them silently would emit too-small balls via negative-tail
    # absolutization
    cases = [
        dict(h=acb(0.6), r=0.5, v0=[acb(1)], why="|h| >= r"),
        dict(h=0.05, r=-0.5, v0=[acb(1)], why="negative r"),
        dict(h=0.05, r=0.5, v0=[acb(1), acb(2)], why="v0 length != RK"),
    ]
    for c in cases:
        try:
            mod.step(NUMP, DENP, 1, acb(1), c["v0"], c["h"], c["r"], 6)
            print(f"FAIL contract: {c['why']} was ACCEPTED (must refuse)")
            return 1
        except ValueError:
            pass
    print("PASS march endpoint (K=30 + tail-load-bearing K=6, 256-bit truth,"
          " 3 contract refusals)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
