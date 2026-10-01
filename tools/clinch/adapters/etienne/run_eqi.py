"""Certificate rerun: the two equal-I multisample fits (tier A crashed at
the MPEqI structural-zero assertion before these ran; ss receipts landed).
Also reconstructs the tier-A read log (same deterministic path list)."""
import csv
import json
import math
import os
import sys
import time

sys.dont_write_bytecode = True

from . import ETIENNE_DIR, DATA_DIR, PILOT_DIR, RECEIPTS_DIR, ident
from . import gates as GT
from .oracle_eqi import EtienneEqIOracle
from .run_ss import run_one, u_from_theta_I, FITS_SS, SIXFOREST

from flint import arb, ctx


def main():
    # reconstruct the tier-A read audit (deterministic list, same files)
    for _, rel, _ in FITS_SS:
        ident.log_read(os.path.join(ETIENNE_DIR, rel),
                       "abundance vector (tier A, reconstructed log)")
    ident.log_read(os.path.join(PILOT_DIR, "PHASE2_CERTIFIED_FITS.json"),
                   "study fit receipts (tier A)")
    ident.log_read(os.path.join(DATA_DIR, "packages/bci_dbh1cm_c1.txt"),
                   "BCI dbh>=1cm abundance vector (tier A)")
    ident.log_read(os.path.join(PILOT_DIR, "BCI_235K_CERTIFIED.json"),
                   "235k fit receipt (tier A)")
    ident.log_read(os.path.join(PILOT_DIR, "SIXFOREST_PM_RF.json"),
                   "six-forest receipt (tier A)")
    for f in SIXFOREST:
        ident.log_read(os.path.join(DATA_DIR,
                                    f"volkov2005/derived/{f}_abundvec.txt"),
                       "forest abundance vector (tier A)")
    manifest = []

    # Panama 3-plot equal-I
    ms_path = os.path.join(DATA_DIR, "packages/panama3_full.txt")
    ident.log_read(ms_path, "panama 3-plot abundance matrix")
    rows = [tuple(int(x) for x in line.split()) for line in open(ms_path)]
    r_path = os.path.join(PILOT_DIR, "PHASE2_MS_EQI_PANAMA.json")
    ident.log_read(r_path, "eqI panama3 receipt (candidate)")
    r = json.load(open(r_path))
    o = EtienneEqIOracle(rows, tag="eqi_panama3")
    u = u_from_theta_I(r["theta_mid"], r["I_mid"])
    mpm = GT.MPEqI(rows)
    old = ctx.prec
    try:
        ctx.prec = 256
        x = ([[arb(u[1])]], [arb(u[0])])
        nll = o.nll(x)
    finally:
        ctx.prec = old
    gp = dict(value_receipt_register=GT.value_gate(nll, r["lnP"]),
              a5_adjudication_register=GT.adjudication_gate(
                  o, mpm, [(u[0] + 0.08, u[1] - 0.1)]),
              fd_hessian_register=GT.fd_hessian_gate(o, u))
    prov = dict(receipt=r_path, receipt_sha=ident.sha256(r_path),
                theta=r["theta_mid"], I=r["I_mid"],
                published_register="float64 midpoints of the study's "
                                   "certified equal-I fit (R4)")
    pr = run_one("eqi_panama3", o, u, prov, gp)
    manifest.append(dict(tag="eqi_panama3",
                         products={k: v["verdict"] for k, v in pr.items()
                                   if k != "polish"}))
    print("eqi_panama3", manifest[-1]["products"], flush=True)

    # Panama-65 equal-I
    p65_path = os.path.join(
        DATA_DIR,
        "panama65/derived/panama65_recent_abundance_matrix_dbh100.csv")
    ident.log_read(p65_path, "panama-65 abundance matrix (dbh>=100mm)")
    rows65 = []
    with open(p65_path) as f:
        rd = csv.reader(f)
        next(rd)
        for line in rd:
            rows65.append(tuple(int(x) for x in line[1:]))
    r_path = os.path.join(PILOT_DIR, "PANAMA65_EQI_CERTIFIED.json")
    ident.log_read(r_path, "eqI panama65 receipt (candidate)")
    r = json.load(open(r_path))
    o = EtienneEqIOracle(rows65, tag="eqi_panama65")
    assert o.S == r["S"] and o.N == r["N"] and sum(o.Js) == r["Jtot"]
    u = u_from_theta_I(r["theta_mid"], r["I_mid"])
    old = ctx.prec
    try:
        ctx.prec = 256
        x = ([[arb(u[1])]], [arb(u[0])])
        nll = o.nll(x)
    finally:
        ctx.prec = old
    gp = dict(value_receipt_register=GT.value_gate(nll, r["lnP"]),
              fd_hessian_register=GT.fd_hessian_gate(o, u))
    prov = dict(receipt=r_path, receipt_sha=ident.sha256(r_path),
                theta=r["theta_mid"], I=r["I_mid"],
                published_register="float64 midpoints of the study's "
                                   "certified equal-I network fit (R6)")
    pr = run_one("eqi_panama65", o, u, prov, gp)
    manifest.append(dict(tag="eqi_panama65",
                         products={k: v["verdict"] for k, v in pr.items()
                                   if k != "polish"}))
    print("eqi_panama65", manifest[-1]["products"], flush=True)

    json.dump(manifest, open(os.path.join(RECEIPTS_DIR,
                                          "MANIFEST_TIER_A_EQI.json"),
                             "w"), indent=1)
    ident.flush_read_log("tier_a")
    print("EQI DONE", flush=True)


if __name__ == "__main__":
    main()
