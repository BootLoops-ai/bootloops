#!/usr/bin/env python3
"""run_r4b — certificate pass for the R4B sealed-contest primaries
(hi_G1, lo_G2). Per target: identity gates (fail-closed) -> direct
engine.certify at the banked vector AS GIVEN -> on refusal newton_polish
-> certify at the polished center (separate, labelled product) + shift
ledger. cert.write_cert per product. Banked numbers never modified.

Reads the originating study's un-shipped fit artifacts and helper module
(fit_multipliers) through the CERTPASS_REFERENCE_DIR environment variable;
refuses loudly when it is unset (see the package __init__).

Usage: python3 run_r4b.py [hi_G1] [lo_G2]
"""
import hashlib
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))          # certpass package importable
sys.path.insert(0, HERE)
import __init__ as pkg                             # noqa (path setup)
HR = pkg.reference_root()   # un-shipped study tree (round4/, stage2/); env-gated
sys.path.insert(0, f"{HR}/round4/r4b/code")
sys.path.insert(0, f"{HR}/stage2/build/playoff")

import fit_multipliers as FM                       # noqa
import clinch                                      # noqa
from clinch import engine, cert                    # noqa
from oracle_coded import CodedOracle, Part         # noqa
import gates_coded as GC                           # noqa

RECEIPTS = pkg.RECEIPTS_DIR
DATA = pkg.DATA_DIR


def sha16(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()[:16]


def log_for(tag):
    def log(s):
        print(f"[{tag}] {s}", flush=True)
    return log


def run_target(tag):
    log = log_for(tag)
    t0 = time.time()
    z = np.load(os.path.join(DATA, f"r4b_{tag}.npz"))
    data = dict(y=z["y"], logref=z["logref"], Zg=z["Zg"],
                l2=float(z["l2"]))
    th = np.asarray(z["x_banked"], float)
    p = len(th)
    part = Part(blocks=[list(range(1, p))], block_names=["betas"],
                border_idx=[0])
    log(f"n={len(data['y'])} p={p} l2={data['l2']}")
    oracle = CodedOracle(data, part, tag)
    nll = FM.coded_nll_factory(data["y"], data["logref"], data["Zg"],
                               data["l2"])
    mp_obj = GC.MpCoded(data)
    summary = dict(tag=f"r4b_{tag}", p=p, n=int(len(data["y"])),
                   l2=data["l2"],
                   banked_fit=f"round4/r4b/fits/{tag}.json",
                   banked_fit_sha16=sha16(f"{HR}/round4/r4b/fits/{tag}.json"),
                   theta=[float(v) for v in th],
                   objective="capped coded log-lik mean + l2*||beta||^2 "
                             "(fit_multipliers.coded_nll_factory of record)")
    log("gate: value")
    gv = GC.gate_value(oracle, mp_obj, nll, float(z["tune_nll_banked"]), th)
    log(f"  arb-ref {gv['dev_arb_vs_ref_float64']:.3g} "
        f"arb-mp50 {gv['dev_arb_vs_mp50']:.3g} "
        f"ref-banked {gv['dev_ref_vs_banked']:.3g} PASS={gv['PASS']}")
    log("gate: gradient (adjudicated)")
    gg, _ = GC.gate_gradient_adjudicated(oracle, mp_obj, nll, th)
    log(f"  ref-analytic {gg['ref_analytic_metric']:.3g} "
        f"mp50 {gg['adjudicated_worst']:.3g} PASS={gg['PASS']}")
    log("gate: fd-hessian")
    gh = GC.gate_fd_hessian(oracle, th)
    log(f"  {gh['metric_max']:.3g} PASS={gh['PASS']}")
    gts = dict(value=gv, gradient_adjudicated=gg, fd_hessian=gh)
    summary["gates"] = gts
    if not (gv["PASS"] and gg["PASS"] and gh["PASS"]):
        summary["verdict"] = "GATES-FAILED"
        log("GATES FAILED — no certificate attempted (fail-closed)")
        return summary
    log("certify: direct")
    res = engine.certify(oracle, th, log=log)
    out_dir = os.path.join(RECEIPTS, f"r4b_{tag}")
    extra = dict(candidate_provenance="banked-full-vector",
                 fit_artifact=f"round4/r4b/fits/{tag}.json",
                 banked_tune_nll=float(z["tune_nll_banked"]),
                 rail_census="clip |eta|<4 and p3-floor 1e-9 provably "
                             "inactive over every evaluated box "
                             "(fail-closed oracle rails; zero crossings)")
    cert.write_cert(res, out_dir, f"r4b_{tag}", gates={
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
        fb = float(nll(th)[0])
        fp = float(nll(thp)[0])
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
        cert.write_cert(resp, os.path.join(RECEIPTS,
                                           f"r4b_{tag}_polished"),
                        f"r4b_{tag}_polished", gates=None, extra=extra_p)
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
    tags = args if args else ["hi_G1", "lo_G2"]
    os.makedirs(RECEIPTS, exist_ok=True)
    mpath = os.path.join(RECEIPTS, "RECEIPTS_certpass.json")
    master = dict(stamp=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                  clinch_version=clinch.__version__,
                  prec_bits=engine.PREC_CERT, targets={})
    if os.path.exists(mpath):
        master["targets"] = json.load(open(mpath)).get("targets", {})
    for tag in tags:
        s = run_target(tag)
        master["targets"][f"r4b_{tag}"] = s
        json.dump(master, open(mpath + ".part", "w"), indent=1,
                  default=float)
        os.replace(mpath + ".part", mpath)
        print(f"[{tag}] DONE in {s.get('wall_s')}s", flush=True)
    print("R4B DONE", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
