#!/usr/bin/env python3
"""run_r4a — certificate pass for the R4A sealed-contest entrant optima
(F1F2_F6 hi/lo, conditional_spatial, per SEALED_VERDICT.json). Same arc
as run_r4b: gates -> direct certify -> polish + re-certify on refusal.

Reads the originating study's un-shipped sealed-fit artifact
(round4/r4a/build/params_sealed.json) through the CERTPASS_REFERENCE_DIR
environment variable; refuses loudly when it is unset.

Usage: python3 run_r4a.py [hi] [lo]
"""
import hashlib
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import __init__ as pkg                             # noqa
HR = pkg.reference_root()   # un-shipped study tree (round4/); env-gated

import clinch                                      # noqa
from clinch import engine, cert                    # noqa
from oracle_cond import CondOracle, Part           # noqa
import gates_cond as GCD                           # noqa

RECEIPTS = pkg.RECEIPTS_DIR
DATA = pkg.DATA_DIR
FD_COORDS = [0, 4, 10, 15, 18, 21]     # receipted FD-Hessian column subset


def sha16(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()[:16]


def log_for(tag):
    def log(s):
        print(f"[{tag}] {s}", flush=True)
    return log


def run_target(fl):
    tag = f"r4a_{fl}_F1F2F6"
    log = log_for(tag)
    t0 = time.time()
    z = np.load(os.path.join(DATA, f"r4a_{fl}_F1F2F6.npz"))
    data = dict(Xs=z["Xs"], o=z["o"], y=z["y"], ridge=float(z["ridge"]))
    th = np.asarray(z["beta_banked"], float)
    p = len(th)
    part = Part(blocks=[list(range(0, 15))], block_names=["F1"],
                border_idx=list(range(15, 22)))
    oracle = CondOracle(data, part, tag)
    log(f"n_days={data['Xs'].shape[0]} active={oracle.n_active} "
        f"k={data['Xs'].shape[1]} p={p} ridge={data['ridge']}")
    fref = GCD.frozen_closure(data)
    log("building mp50 transliteration")
    mp_obj = GCD.MpCond(data)
    summary = dict(tag=tag, p=p,
                   n_days=int(data["Xs"].shape[0]),
                   n_active_days=int(oracle.n_active),
                   ridge=data["ridge"],
                   banked_fit="round4/r4a/build/params_sealed.json"
                              f" [F1F2_F6_{fl}]",
                   banked_fit_sha16=sha16(
                       f"{HR}/round4/r4a/build/params_sealed.json"),
                   theta=[float(v) for v in th],
                   objective="fit_cond conditional-multinomial NLL + "
                             "ridge (frozen step4 text of record)")
    log("gate: value")
    gv = GCD.gate_value(oracle, mp_obj, fref, th,
                        f_banked_receipt=float(z["f_at_banked"]))
    log(f"  arb-ref {gv['dev_arb_vs_ref_float64']:.3g} "
        f"arb-mp50 {gv['dev_arb_vs_mp50']:.3g} PASS={gv['PASS']}")
    log("gate: gradient (adjudicated)")
    gg, _ = GCD.gate_gradient_adjudicated(oracle, mp_obj, fref, th)
    log(f"  ref-analytic {gg['ref_analytic_metric']:.3g} "
        f"mp50 {gg['adjudicated_worst']:.3g} PASS={gg['PASS']}")
    log("gate: fd-hessian (receipted column subset)")
    gh = GCD.gate_fd_hessian(oracle, th, coords=FD_COORDS)
    gh["coords"] = FD_COORDS
    log(f"  {gh['metric_max']:.3g} PASS={gh['PASS']}")
    gts = dict(value=gv, gradient_adjudicated=gg, fd_hessian=gh)
    summary["gates"] = gts
    if not (gv["PASS"] and gg["PASS"] and gh["PASS"]):
        summary["verdict"] = "GATES-FAILED"
        log("GATES FAILED — no certificate attempted (fail-closed)")
        return summary
    log("certify: direct")
    res = engine.certify(oracle, th, log=log)
    extra = dict(candidate_provenance="banked-full-vector",
                 fit_artifact="round4/r4a/build/params_sealed.json"
                              f" [F1F2_F6_{fl}]",
                 assembly=dict(mu_sd_byte_exact=True,
                               frozen_fit_cond_refit_dev=float(
                                   z["refit_dev"]),
                               gmax_at_banked=float(z["gmax_at_banked"])),
                 rails="none in the fit objective (ZCLIP emission-only)")
    cert.write_cert(res, os.path.join(RECEIPTS, tag), tag, gates={
        k: {kk: vv for kk, vv in v.items()
            if not isinstance(vv, dict) or k == "value"}
        for k, v in gts.items()}, extra=extra)
    summary["direct"] = dict(verdict=res["verdict"],
                             failed=res.get("failed"),
                             radius_scale=res.get("radius_scale"),
                             max_margin=(res.get("krawczyk") or {})
                             .get("max_margin"),
                             newton_step_scaled=res.get(
                                 "newton_step_scaled_inf"),
                             wall_s=res.get("wall_s_total"))
    log(f"direct: {res['verdict']} ({res.get('failed')})")
    if res["verdict"] == "REFUSED":
        log("newton_polish")
        thp, prcpt = engine.newton_polish(oracle, th, log=log)
        log("certify: polished center")
        resp = engine.certify(oracle, thp, log=log)
        fb = float(fref(th)[0])
        fp = float(fref(thp)[0])
        shift = dict(
            distance_inf=float(np.max(np.abs(thp - th))),
            per_coord={int(a): dict(banked=float(th[a]),
                                    polished=float(thp[a]),
                                    shift=float(thp[a] - th[a]))
                       for a in range(p)},
            f_float_at_banked=fb, f_float_at_polished=fp,
            f_drop=fb - fp)
        extra_p = dict(extra, polish=prcpt, shift_vs_banked=shift,
                       polished_center=[float(v) for v in thp])
        cert.write_cert(resp, os.path.join(RECEIPTS, tag + "_polished"),
                        tag + "_polished", gates=None, extra=extra_p)
        summary["polished"] = dict(
            verdict=resp["verdict"], failed=resp.get("failed"),
            max_margin=(resp.get("krawczyk") or {}).get("max_margin"),
            pd_min=(min(((resp.get("pd") or {}).get("margins") or
                         {"": float("nan")}).values())
                    if resp.get("pd") else None),
            distance_inf=shift["distance_inf"],
            f_drop=shift["f_drop"], wall_s=resp.get("wall_s_total"))
        log(f"polished: {resp['verdict']} "
            f"(dist {shift['distance_inf']:.3g}, f drop "
            f"{shift['f_drop']:.3g})")
    summary["wall_s"] = round(time.time() - t0, 1)
    return summary


def main(args):
    floors = args if args else ["hi", "lo"]
    os.makedirs(RECEIPTS, exist_ok=True)
    mpath = os.path.join(RECEIPTS,
                         "RECEIPTS_r4a_" + "_".join(floors) + ".json")
    master = dict(stamp=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                  clinch_version=clinch.__version__,
                  prec_bits=engine.PREC_CERT, targets={})
    if os.path.exists(mpath):
        old = json.load(open(mpath))
        master["targets"] = old.get("targets", {})
        master["stamp"] = old.get("stamp", master["stamp"])
    for fl in floors:
        s = run_target(fl)
        master["targets"][s["tag"]] = s
        json.dump(master, open(mpath + ".part", "w"), indent=1,
                  default=float)
        os.replace(mpath + ".part", mpath)
        print(f"[{fl}] DONE in {s.get('wall_s')}s", flush=True)
    print("R4A DONE", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
