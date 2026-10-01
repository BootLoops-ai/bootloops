#!/usr/bin/env python3
"""assemble_r4a — R4A sealed-fit design assembly + frozen-code identity.

The banked sealed optimum (round4/r4a/build/params_sealed.json,
F1F2_F6_{hi,lo}: the SEALED_VERDICT entrant) is variant-B fit_cond on
ref-emission windows 2010-2014. This script:

1. rebuilds the fit design exactly as the frozen step4_sealed_eval.py
   built it (cube.npz + F1/F6 log1p + ref-emission day stack + full-cube
   standardize_stats), and verifies the banked per-column mu/sd match the
   recomputation EXACTLY (byte-equal floats) — the design-pipeline
   identity receipt;
2. extracts fit_cond/standardize_stats FROM THE FROZEN SCRIPT'S OWN TEXT
   (ast extraction, sha16 receipted — step4 is a script, not importable)
   and receipts the frozen closure's (value, grad) at the banked beta;
3. saves one npz per floor with everything the oracle needs.

Reads the originating study's un-shipped tree (round4/, stage2/) and its
helper module (playoff) through the CERTPASS_REFERENCE_DIR environment
variable; refuses loudly when it is unset.
"""
import ast
import hashlib
import json
import os
import sys

import numpy as np
from scipy.optimize import minimize            # noqa (frozen code needs it)
from scipy.special import logsumexp, softmax   # noqa

HERE = os.path.dirname(os.path.abspath(__file__))
HR = (os.environ.get("CERTPASS_REFERENCE_DIR")
      or os.environ.get("CERTPASS_LANE_ROOT", ""))
if not HR or not os.path.isdir(HR):
    raise SystemExit(
        "certpass: CERTPASS_REFERENCE_DIR is not set (or not a directory); this "
        "script reads the originating study's un-shipped artifacts (round4/, "
        "stage2/) and its helper module playoff. Point CERTPASS_REFERENCE_DIR at "
        "that tree.")
R4A = f"{HR}/round4/r4a"
S2 = f"{HR}/stage2"
sys.path.insert(0, f"{S2}/build/playoff")
import playoff                                  # noqa

STEP4 = f"{R4A}/build/step4_sealed_eval.py"
FIT_YEARS = [2010, 2011, 2012, 2013, 2014]


def sha16(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()[:16]


def frozen_functions():
    """Compile fit_cond + standardize_stats from the frozen step4 text."""
    src = open(STEP4).read()
    tree = ast.parse(src)
    keep = [n for n in tree.body if isinstance(n, ast.FunctionDef)
            and n.name in ("fit_cond", "standardize_stats")]
    mod = ast.Module(body=keep, type_ignores=[])
    ns = dict(np=np, minimize=minimize, logsumexp=logsumexp,
              softmax=softmax)
    exec(compile(mod, STEP4, "exec"), ns)
    return ns["fit_cond"], ns["standardize_stats"], sha16(STEP4)


def main():
    fit_cond, standardize_stats, step4_sha = frozen_functions()
    cube = np.load(f"{R4A}/build/cube.npz")
    days_cube = cube["days"]
    X = cube["X"]
    names = [str(s) for s in cube["names"]]
    day2row = {int(d): i for i, d in enumerate(days_cube)}
    Xt = X.astype(np.float64).copy()
    for j, nm in enumerate(names):
        if nm.startswith(("F1", "F6")):
            Xt[:, :, j] = np.log1p(Xt[:, :, j])

    def load_ref(fl, year):
        z = np.load(f"{S2}/components/ref_etas/forecasts/{fl}_{year}.npz")
        return z["days"].astype(np.int64), z["lam"]

    def stack(years, fl):
        dd, ll = [], []
        for y in years:
            d, lam = load_ref(fl, y)
            dd.append(d)
            ll.append(lam)
        d = np.concatenate(dd)
        lam = np.concatenate(ll)
        rows = np.array([day2row[int(x)] for x in d])
        return d, lam, Xt[rows]

    sealed = json.load(open(f"{R4A}/build/params_sealed.json"))
    truths = {fl: playoff.load_truth(mag)
              for fl, mag in playoff.FLOORS.items()}
    meta = dict(step4_sha16=step4_sha,
                cube_sha16=sha16(f"{R4A}/build/cube.npz"),
                params_sealed_sha16=sha16(
                    f"{R4A}/build/params_sealed.json"),
                fit_years=FIT_YEARS, group="F1F2_F6")
    for fl in ("hi", "lo"):
        rec = sealed[f"F1F2_F6_{fl}"]
        assert rec["fit_windows"] == FIT_YEARS
        cols = [j for j, nm in enumerate(names)
                if nm.startswith(("F1", "F2", "F6"))]
        assert [names[j] for j in cols] == rec["cols"]
        d_f, lam_f, Xf_f = stack(FIT_YEARS, fl)
        mu, sd = standardize_stats(Xf_f)
        b_mu = np.asarray(rec["mu"], float)
        b_sd = np.asarray(rec["sd"], float)
        mu_exact = bool(np.array_equal(mu[cols], b_mu))
        sd_exact = bool(np.array_equal(sd[cols], b_sd))
        Xs = (Xf_f[:, :, cols] - b_mu) / b_sd    # banked constants of record
        y = truths[fl][d_f]
        beta = np.asarray(rec["beta"], float)
        ridge = float(rec["ridge"])
        assert rec["weight"] == "equal"
        # frozen-closure receipt at the banked beta
        o = np.log(np.maximum(lam_f, 1e-9))
        yy = y.astype(np.float64)
        Yd = yy.sum(1)
        wd = 1.0 / np.maximum(Yd, 1.0)
        gconst = -((wd[:, None] * yy)[:, :, None] * Xs).sum((0, 1))
        wY = wd * Yd

        def nll_frozen(b):
            z = o + Xs @ b
            lse = logsumexp(z, axis=1)
            f = (-((wd[:, None] * yy) * z).sum() + (wY * lse).sum()
                 + ridge * (b ** 2).sum())
            w = softmax(z, axis=1)
            g = (gconst + ((wY[:, None] * w)[:, :, None] * Xs).sum((0, 1))
                 + 2 * ridge * b)
            return f, g

        f0, g0 = nll_frozen(beta)
        # independent same-data refit through the FROZEN fit_cond
        b_ref, ok_ref = fit_cond(Xs, lam_f, y, ridge, "equal")
        dev_refit = float(np.max(np.abs(b_ref - beta)))
        print(fl, "n_days", len(d_f), "p", len(cols),
              "mu_exact", mu_exact, "sd_exact", sd_exact,
              "f(banked)", f0, "gmax", np.abs(g0).max(),
              "refit_dev", dev_refit, "refit_ok", ok_ref)
        path = os.path.join(HERE, "data", f"r4a_{fl}_F1F2F6.npz")
        np.savez_compressed(
            path, Xs=Xs, o=o, y=y, wd=wd, wY=wY,
            beta_banked=beta, ridge=np.float64(ridge),
            f_at_banked=np.float64(f0), gmax_at_banked=np.float64(
                np.abs(g0).max()),
            refit_dev=np.float64(dev_refit))
        meta[fl] = dict(n_days=int(len(d_f)), k=int(Xs.shape[1]),
                        p=int(Xs.shape[2]), mu_exact=mu_exact,
                        sd_exact=sd_exact, f_at_banked=float(f0),
                        gmax_at_banked=float(np.abs(g0).max()),
                        refit_dev=dev_refit, refit_converged=bool(ok_ref))
    json.dump(meta, open(os.path.join(HERE, "data",
                                      "r4a_assembly_meta.json"), "w"),
              indent=1)
    print("meta ->", os.path.join(HERE, "data", "r4a_assembly_meta.json"))


if __name__ == "__main__":
    main()
