#!/usr/bin/env python3
"""assemble_r4b — R4B tune-era design assembly, BY IDENTITY.

Rebuilds (y, logref, Z) exactly as r4b/code/fit_multipliers.main() built
them for the banked fits (cube_B0 path, TUNE = 2010-2012), using the
frozen module's OWN functions (load_cube, ref_forecasts, playoff
load_truth) — imported read-only, never transliterated. Saves one npz per
floor + sha16 receipts of every input consumed.

Reads the originating study's un-shipped tree (round4/, stage2/) and its
helper modules (fit_multipliers, playoff) through the
CERTPASS_REFERENCE_DIR environment variable; refuses loudly when it is unset.
"""
import hashlib
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
HR = (os.environ.get("CERTPASS_REFERENCE_DIR")
      or os.environ.get("CERTPASS_LANE_ROOT", ""))
if not HR or not os.path.isdir(HR):
    raise SystemExit(
        "certpass: CERTPASS_REFERENCE_DIR is not set (or not a directory); this "
        "script reads the originating study's un-shipped artifacts (round4/, "
        "stage2/) and its helper modules fit_multipliers/playoff. Point "
        "CERTPASS_REFERENCE_DIR at that tree.")
RB = f"{HR}/round4/r4b"
S2 = f"{HR}/stage2"
sys.path.insert(0, f"{RB}/code")
sys.path.insert(0, f"{S2}/build/playoff")

import fit_multipliers as FM                                    # noqa
from playoff import load_truth, FLOORS                          # noqa


def sha16(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()[:16]


def main():
    cube_p = f"{RB}/features/cube_B0.npz"
    Z, span0, meta = FM.load_cube(cube_p)
    groups = meta["groups"]
    cols_of = {}
    cols = []
    for g in FM.GROUPS_ORDER:
        cols = cols + groups[g]
        cols_of["G" + str(len(cols_of) + 1)] = list(cols)
    out_meta = dict(
        cube_sha16=sha16(cube_p),
        meta_sha16=sha16(f"{RB}/features/FEATURE_META.json"),
        train_sha_file=FM.TRAIN_SHA,
        span0=int(span0), tune=list(FM.TUNE),
        clip_eta=float(FM.CLIP_ETA))
    for floor, mag in (("hi", 2.0), ("lo", 3.0)):
        truth = load_truth(mag)
        ytr, ltr, Ztr = [], [], []
        n_fit = 0
        for Y in FM.TUNE:
            days, lam = FM.ref_forecasts(floor, Y)
            if Y in (2010, 2011):
                n_fit += len(days) * lam.shape[1]
            ytr.append(truth[days].ravel())
            ltr.append(np.log(np.maximum(lam, 1e-12)).ravel())
            Ztr.append(Z[days - span0].reshape(-1, Z.shape[2]))
        ytr = np.concatenate(ytr)
        ltr = np.concatenate(ltr)
        Ztr = np.concatenate(Ztr).astype(np.float64)
        for cd in (("G1",) if floor == "hi" else ("G2",)):
            fit = json.load(open(f"{RB}/fits/{floor}_{cd}.json"))
            cc = cols_of[cd]
            assert cc == fit["cols"], (cc, fit["cols"])
            path = os.path.join(HERE, "data", f"r4b_{floor}_{cd}.npz")
            np.savez_compressed(
                path, y=ytr, logref=ltr, Zg=Ztr[:, cc],
                cols=np.asarray(cc, np.int64),
                x_banked=np.array([fit["b0"]] + fit["beta"], np.float64),
                l2=np.float64(fit["l2"]),
                tune_nll_banked=np.float64(fit["tune_nll"]),
                n_fit=np.int64(n_fit))
            print(floor, cd, "n", len(ytr), "p", len(cc) + 1,
                  "l2", fit["l2"], "->", path)
            out_meta[f"{floor}_{cd}_fit_sha16"] = sha16(
                f"{RB}/fits/{floor}_{cd}.json")
    out_meta["assembler"] = "fit_multipliers by identity (load_cube/"\
        "ref_forecasts) + playoff.load_truth by identity"
    json.dump(out_meta, open(os.path.join(HERE, "data",
                                          "r4b_assembly_meta.json"), "w"),
              indent=1)
    print("meta ->", os.path.join(HERE, "data", "r4b_assembly_meta.json"))


if __name__ == "__main__":
    main()
