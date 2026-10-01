"""Certificate run, tier D — the flagship: the reference Panama-65
different-I fit (66 parameters; PANAMA65_DIFFI_FIT.json; the
study's named hole "optimum at search precision, not certified").

Direct certificate at the published point AS GIVEN first; the polished-
center certificate is the separate labelled product if direct refuses
(the fit13/fit17 flagship pattern). Battery (incl. E7 real-scale identity)
must be OVERALL PASS on disk before any certificate is attempted.
"""
import json
import math
import os
import sys
import time

sys.dont_write_bytecode = True

import numpy as np

from . import DATA_DIR, PILOT_DIR, RECEIPTS_DIR, ident
from .oracle_diffi import EtienneDiffIOracle

from clinch import engine, cert

STAMP = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
BUDGET_S = 7200.0
POLISH_PREC = 160
POLISH_TOL_SIGMA = 1e-8
POLISH_MAXIT = 16


def load():
    import csv
    p65 = os.path.join(
        DATA_DIR,
        "panama65/derived/panama65_recent_abundance_matrix_dbh100.csv")
    ident.log_read(p65, "panama-65 abundance matrix (fit data of record)")
    rows = []
    with open(p65) as f:
        rd = csv.reader(f)
        next(rd)
        for line in rd:
            rows.append(tuple(int(x) for x in line[1:]))
    r_path = os.path.join(PILOT_DIR, "PANAMA65_DIFFI_FIT.json")
    ident.log_read(r_path, "diffI fit receipt (candidate of record)")
    rec = json.load(open(r_path))
    return rows, rec, r_path, p65


def x_at(u, prec):
    from flint import arb, ctx
    ctx.prec = prec
    return ([[arb(v) for v in u[1:]]], [arb(u[0])])


def frozen_h_polish(o, u0, H_mid, log):
    """Heuristic frozen-Hessian Newton polish: fresh gradient each step,
    ONE fixed midpoint NLL Hessian (rigor lives only in the subsequent
    certificate; this is the labelled polish product's producer)."""
    u = np.array(u0, dtype=float)
    d = np.sqrt(np.maximum(np.diag(H_mid), 1e-6))
    hist = []
    for it in range(POLISH_MAXIT):
        t0 = time.time()
        r = o.evaluate(x_at(list(u), POLISH_PREC), order=1)
        g = np.array([float(v.mid()) for v in r["grad"]])
        g = -g                                  # NLL gradient
        step = np.linalg.solve(H_mid, g)
        scaled = float(np.max(np.abs(step) * d))
        hist.append(dict(it=it, step_inf=float(np.max(np.abs(step))),
                         step_scaled_inf=scaled,
                         wall_s=round(time.time() - t0, 1)))
        log(f"polish it{it}: scaled step {scaled:.3g} "
            f"({hist[-1]['wall_s']}s)")
        u = u - step
        if scaled <= POLISH_TOL_SIGMA:
            break
    dist = float(np.max(np.abs(u - np.array(u0))))
    return list(u), dict(iterations=hist, distance_inf=dist,
                         frozen_hessian=True, prec_bits=POLISH_PREC,
                         tol_sigma=POLISH_TOL_SIGMA,
                         note="frozen-midpoint-Hessian Newton (heuristic; "
                              "the certificate at the polished center "
                              "carries the rigor)")


def main():
    out_root = os.path.join(RECEIPTS_DIR, "diffi_panama65")
    os.makedirs(out_root, exist_ok=True)
    logf = open(os.path.join(out_root, "LOG.txt"), "a")

    def log(s):
        logf.write(time.strftime("[%H:%M:%SZ] ", time.gmtime()) + s + "\n")
        logf.flush()

    log(f"== diffi_panama65 @ {STAMP}")
    # battery must be green (incl. E7) before certificates
    bat_path = os.path.join(RECEIPTS_DIR, "ADAPTER_BATTERY.json")
    bat = json.load(open(bat_path))
    assert bat["overall"] == "PASS", "adapter battery not green"
    assert "E7_diffi_real_scale_identity" in bat["legs"], "E7 missing"

    rows, rec, r_path, p65 = load()
    o = EtienneDiffIOracle(rows, tag="diffi_panama65")
    u_pub = [math.log(rec["theta"])] + [math.log(v) for v in rec["I"]]

    # value identity gate at the published point (receipt register)
    from flint import arb, ctx
    ctx.prec = 224
    lnP = o.evaluate(x_at(u_pub, 224), order=0)["lnP"]
    ref = arb(rec["lnL_certified"])
    vg = dict(register="receipt",
              my_lnP=str(lnP), receipt=rec["lnL_certified"],
              mid_abs_diff=abs(float((lnP - ref).mid())),
              bar=1e-9,
              lnL_search_prec_receipt=rec["lnL_search_prec"],
              search_vs_certified_gap=abs(rec["lnL_search_prec"]
                                          - float(ref.mid())))
    vg["ok"] = bool(vg["mid_abs_diff"] <= vg["bar"])
    gates = dict(value_receipt_register=vg,
                 battery=dict(path=bat_path, overall=bat["overall"],
                              e7=bat["legs"]
                              ["E7_diffi_real_scale_identity"]))
    log(f"value gate: diff {vg['mid_abs_diff']:.3g} ok={vg['ok']}")

    prov = dict(receipt=r_path, receipt_sha=ident.sha256(r_path),
                data=p65, data_sha=ident.sha256(p65),
                theta=rec["theta"], n_plots=len(rec["I"]),
                published_register="coordinate-ascent optimum at search "
                                   "precision (the study's own caveat: "
                                   "'optimum at search precision, not "
                                   "certified')")
    extra = dict(adapter="clinch.adapters.etienne",
                 coordinates="(lntheta | lnI_1..lnI_65)",
                 structure="ONE dense 65-dim lnI block + 1-dim lntheta "
                           "border (I_p x I_q couplings are nonzero; 65 "
                           "independent dim-1 blocks would be unsound)",
                 objective="NLL = -lnP",
                 candidate_provenance=prov, run_stamp_utc=STAMP)

    t0 = time.time()
    res = engine.certify(o, u_pub, budget_s=BUDGET_S, log=log)
    log(f"DIRECT: {res['verdict']} ({time.time()-t0:.0f}s)")
    cert.write_cert(res, out_root, "diffi_panama65", gates=gates,
                    extra=dict(extra, product="direct"))

    if res["verdict"].startswith("REFUSED"):
        # optional probe-cache reuse (CLINCH_DIFFI_H_MID points at a saved
        # H midpoint .npy; unset => recompute below)
        H_path = os.environ.get("CLINCH_DIFFI_H_MID", "")
        if H_path and os.path.exists(H_path):
            H_mid = np.load(H_path)
            log("polish: frozen H from probe cache")
        else:
            r2 = o.evaluate(x_at(u_pub, 128), order=2)
            H_mid = -np.array([[float(r2["hess"][i][j].mid())
                                for j in range(66)] for i in range(66)])
            log("polish: frozen H computed at 128 bits")
        thp, prcpt = frozen_h_polish(o, u_pub, H_mid, log)
        t0 = time.time()
        resp = engine.certify(o, thp, budget_s=BUDGET_S, log=log)
        log(f"POLISHED: {resp['verdict']} ({time.time()-t0:.0f}s)")
        pdir = os.path.join(out_root, "polished")
        # distances in original parameters
        th_new = math.exp(thp[0])
        dI = [abs(math.exp(thp[1 + p]) - rec["I"][p])
              for p in range(len(rec["I"]))]
        cert.write_cert(resp, pdir, "diffi_panama65_polished",
                        gates=gates,
                        extra=dict(extra, product="polished_center",
                                   polish=prcpt,
                                   theta_polished=th_new,
                                   theta_published=rec["theta"],
                                   max_abs_I_move=max(dI)))
    json.dump(dict(stamp=STAMP, direct=res["verdict"],
                   polished=(None if not res["verdict"].startswith("REF")
                             else resp["verdict"])),
              open(os.path.join(RECEIPTS_DIR, "MANIFEST_TIER_D.json"),
                   "w"), indent=1)
    ident.flush_read_log("tier_d")
    log("TIER D DONE")
    print("TIER D DONE", flush=True)


if __name__ == "__main__":
    main()
