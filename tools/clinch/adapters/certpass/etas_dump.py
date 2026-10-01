#!/usr/bin/env python3
"""etas_dump — regenerate the ETAS fit state BY THE MACHINERY OF RECORD.

RUNS UNDER build/etas_env/bin/python3 ONLY (the pinned install,
commit 097f08b6...; same env the banked fits ran under). For one banked
fit (e2_params.json + catalog csv):

1. rebuild meta exactly as build/etas_e2_fit.mode_fit does (from the
   banked config block), instantiate ETASParameterCalculation, prepare();
2. verify n_target_events / beta / area match the banked e2_params.json;
3. run ONE expectation_step at the banked theta (the package's own code)
   and dump: pair table (source pos, target pos, dt, r2, m_source),
   source integrals G + l_hat, per-target P_background + tot_rates,
   n_hat — the identity-gate state of record;
4. verify the regenerated state matches the BANKED sources_*.csv /
   trig_and_bg_probs_*.csv of the original fit (max abs dev receipts).

Usage: etas_env/bin/python3 etas_dump.py <e2_params.json> <catalog.csv>
       <banked_sources.csv> <banked_trigbg.csv> <out.npz>
"""
import json
import sys

import numpy as np
import pandas as pd

from etas.inversion import (ETASParameterCalculation,
                            parameter_dict2array,
                            expected_aftershocks)


def main():
    e2_p, cat_p, src_p, trg_p, out_p = sys.argv[1:6]
    e2 = json.load(open(e2_p))
    cfg = e2["config"]
    meta = {
        "name": cfg["name"], "id": cfg["name"], "fn_catalog": cat_p,
        "auxiliary_start": cfg["auxiliary_start"],
        "timewindow_start": cfg["timewindow_start"],
        "timewindow_end": cfg["timewindow_end"],
        "mc": cfg["mc"], "delta_m": cfg["delta_m"],
        "coppersmith_multiplier": cfg["coppersmith_multiplier"],
        "shape_coords": cfg["shape_coords"],
        "theta_0": cfg["theta_0"],
    }
    calc = ETASParameterCalculation(meta)
    calc.prepare()
    assert len(calc.target_events) == e2["n_target_events"], \
        (len(calc.target_events), e2["n_target_events"])
    assert abs(float(calc.beta) - e2["beta"]) < 1e-12
    assert abs(float(calc.area) - e2["area_km2"]) < 1e-6
    assert not getattr(calc, "constraints", []), "constraints not empty"
    theta = parameter_dict2array(e2["final_parameters"])
    mc_min = calc.m_ref - calc.delta_m / 2
    pij, targ, srce, n_hat, i_hat = calc.expectation_step(theta, mc_min)
    # G at the banked theta (the package's own integral)
    srce = srce.copy()
    srce["G"] = expected_aftershocks(
        [srce["source_magnitude"],
         srce["pos_source_to_start_time_distance"],
         srce["source_to_end_time_distance"]],
        [theta[2:], mc_min])
    # ---- reproduction vs the banked fit artifacts -------------------
    bs = pd.read_csv(src_p).set_index("source_id")
    bt = pd.read_csv(trg_p).set_index("target_id")
    a = srce.join(bs, rsuffix="_banked")
    dev_G = float((a["G"] - a["G_banked"]).abs().max())
    dev_lhat = float((a["l_hat"] - a["l_hat_banked"]).abs().max())
    b = targ.join(bt, rsuffix="_banked")
    dev_pbg = float((b["P_background"] - b["P_background_banked"])
                    .abs().max())
    comp_max = float(np.abs(np.asarray(
        pij.reset_index()["source_completeness_above_ref"],
        float)).max())
    zeta_dev = float((targ["zeta_plus_1"] - 1.0).abs().max())
    # ---- dump -------------------------------------------------------
    src_ids = srce.index.to_numpy()
    trg_ids = targ.index.to_numpy()
    spos = {int(v): i for i, v in enumerate(src_ids)}
    tpos = {int(v): i for i, v in enumerate(trg_ids)}
    pr = pij.reset_index()
    pair_s = np.array([spos[int(v)] for v in pr["source_id"]], np.int32)
    pair_t = np.array([tpos[int(v)] for v in pr["target_id"]], np.int32)
    np.savez_compressed(
        out_p,
        pair_s=pair_s, pair_t=pair_t,
        pair_dt=pr["time_distance"].to_numpy(np.float64),
        pair_r2=pr["spatial_distance_squared"].to_numpy(np.float64),
        pair_m=pr["source_magnitude"].to_numpy(np.float64),
        src_m=srce["source_magnitude"].to_numpy(np.float64),
        src_ts=srce["pos_source_to_start_time_distance"]
        .to_numpy(np.float64),
        src_te=srce["source_to_end_time_distance"].to_numpy(np.float64),
        src_G=srce["G"].to_numpy(np.float64),
        src_lhat=srce["l_hat"].to_numpy(np.float64),
        trg_pbg=targ["P_background"].to_numpy(np.float64),
        trg_totrate=(targ["mu"] / targ["P_background"])
        .to_numpy(np.float64),
        theta=np.array([float(v) if v is not None else np.nan
                        for v in theta], np.float64),
        beta=np.float64(calc.beta), area=np.float64(calc.area),
        tw_length=np.float64(calc.timewindow_length),
        m_ref=np.float64(calc.m_ref), delta_m=np.float64(calc.delta_m),
        mc_min=np.float64(mc_min), n_hat=np.float64(n_hat),
        n_hat_banked=np.float64(e2["n_hat"]))
    print(json.dumps(dict(
        n_pairs=int(len(pair_s)), n_sources=int(len(src_ids)),
        n_targets=int(len(trg_ids)),
        dev_G=dev_G, dev_lhat=dev_lhat, dev_pbg=dev_pbg,
        dev_nhat=abs(float(n_hat) - e2["n_hat"]),
        completeness_above_ref_max=comp_max, zeta_dev=zeta_dev,
        out=out_p)))


if __name__ == "__main__":
    main()
