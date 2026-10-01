"""Certificate run, tier A: single-sample interior fits of record
(+ the two equal-I multisample fits, same 2-dim machinery).

Per fit: identity gates -> DIRECT certificate at the candidate AS GIVEN
(the study's receipted floats / printed table row) -> polished-center
certificate as a separate labelled product if direct refuses. Receipts to
v3/receipts/clinch_etienne/<tag>/.
"""
import json
import math
import os
import sys
import time

sys.dont_write_bytecode = True

from . import ETIENNE_DIR, DATA_DIR, PILOT_DIR, RECEIPTS_DIR, ident
from . import gates as GT
from .oracle_ss import EtienneSSOracle
from .oracle_eqi import EtienneEqIOracle

from clinch import engine, cert

STAMP = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def load_vec(path):
    ident.log_read(path, "abundance vector (fit data of record)")
    return sorted(int(x) for x in open(path).read().split())


def load_receipt(path):
    ident.log_read(path, "study fit receipt (candidate of record)")
    return json.load(open(path))


def u_from_theta_m(theta, m, J):
    I = m * (J - 1) / (1.0 - m)
    return [math.log(theta), math.log(I)]


def u_from_theta_I(theta, I):
    return [math.log(theta), math.log(I)]


def run_one(tag, oracle, u_cand, cand_prov, gate_pack, budget=1800.0):
    out_dir = os.path.join(RECEIPTS_DIR, tag)
    os.makedirs(out_dir, exist_ok=True)
    logf = open(os.path.join(out_dir, "LOG.txt"), "a")

    def log(s):
        logf.write(s + "\n")
        logf.flush()

    log(f"== {tag} @ {STAMP} (date -u at run start)")
    t0 = time.time()
    res = engine.certify(oracle, u_cand, budget_s=budget, log=log)
    res_direct = res
    products = dict(direct=res)
    polish_rec = None
    if res["verdict"].startswith("REFUSED"):
        thp, polish_rec = engine.newton_polish(oracle, u_cand, log=log)
        resp = engine.certify(oracle, list(thp), budget_s=budget, log=log)
        products["polished"] = resp
        products["polish"] = polish_rec
    extra = dict(
        adapter="clinch.adapters.etienne", coordinates="(lntheta, lnI)",
        objective="NLL = -lnP (certified strict local min == strict "
                  "local max of the likelihood)",
        candidate_u=list(map(float, u_cand)),
        candidate_provenance=cand_prov,
        run_stamp_utc=STAMP)
    c = cert.write_cert(res_direct, out_dir, tag, gates=gate_pack,
                        extra=dict(extra, product="direct"))
    if "polished" in products:
        pdir = os.path.join(out_dir, "polished")
        cert.write_cert(products["polished"], pdir, tag + "_polished",
                        gates=gate_pack,
                        extra=dict(extra, product="polished_center",
                                   polish=polish_rec))
    log(f"total {time.time() - t0:.1f}s; direct={res_direct['verdict']}"
        + (f"; polished={products['polished']['verdict']}"
           if "polished" in products else ""))
    logf.close()
    return products


def gate_pack_ss(oracle, u_cand, receipt_lnp=None, do_a5=False, mp_model=None,
                 fd=True):
    from flint import arb, ctx
    pack = {}
    old = ctx.prec
    try:
        ctx.prec = 256
        x = ([[arb(u_cand[1])]], [arb(u_cand[0])])
        nll = oracle.nll(x)
    finally:
        ctx.prec = old
    if receipt_lnp is not None:
        pack["value_receipt_register"] = GT.value_gate(nll, receipt_lnp)
    else:
        pack["value_receipt_register"] = dict(
            register="receipt", ok=None,
            note="no receipted lnP ball for this fit; A5 register carries "
                 "the identity claim", oracle_nll=str(nll))
    if do_a5 and mp_model is not None:
        pts = [(u_cand[0] + 0.1, u_cand[1] - 0.15),
               (u_cand[0] - 0.07, u_cand[1] + 0.12)]
        pack["a5_adjudication_register"] = GT.adjudication_gate(
            oracle, mp_model, pts)
    if fd:
        pack["fd_hessian_register"] = GT.fd_hessian_gate(oracle, u_cand)
    return pack


FITS_SS = [
    # tag, data path, receipt source
    ("ss_bci1982_N100", "data/ctfs_bci_panama/derived/bci1982_N100_abundvec.txt", 0),
    ("ss_bci1985_N100", "data/ctfs_bci_panama/derived/bci1985_N100_abundvec.txt", 1),
    ("ss_bci1990_N100", "data/ctfs_bci_panama/derived/bci1990_N100_abundvec.txt", 2),
    ("ss_bci1995_N100", "data/ctfs_bci_panama/derived/bci1995_N100_abundvec.txt", 3),
    ("ss_bci2000_N100", "data/ctfs_bci_panama/derived/bci2000_N100_abundvec.txt", 4),
    ("ss_cm1984", "data/winemiller1990/derived/cm1984_abundvec_sorted.txt", 5),
    ("ss_bci_volkov_21457_225", "data/packages/bci_volkov_21457_225.txt", 6),
    ("ss_bci1982_untb_vintage", "data/ctfs_bci_panama/derived/bci1982_N100_20071118_abundvec.txt", 7),
]

SIXFOREST = ["bci", "korup", "pasoh", "lambir", "yasuni", "sinharaja"]
A5_TAGS = {"ss_bci1982_N100", "ss_cm1984", "ss_worked_example"}


def main():
    manifest = []
    p2_path = os.path.join(PILOT_DIR, "PHASE2_CERTIFIED_FITS.json")
    p2 = load_receipt(p2_path)

    for tag, rel, idx in FITS_SS:
        path = os.path.join(ETIENNE_DIR, rel)
        D = load_vec(path)
        r = p2[idx]
        assert r["J"] == sum(D) and r["S"] == len(D), (tag, "data mismatch")
        o = EtienneSSOracle(D, tag=tag)
        u = u_from_theta_m(r["theta_mid"], r["m_mid"], r["J"])
        mpm = GT.MPSingleSample(D) if tag in A5_TAGS else None
        gp = gate_pack_ss(o, u, receipt_lnp=r["lnP"],
                          do_a5=tag in A5_TAGS, mp_model=mpm)
        prov = dict(receipt=p2_path, receipt_sha=ident.sha256(p2_path),
                    fields=["theta_mid", "m_mid"], tag_in_receipt=r["tag"],
                    theta=r["theta_mid"], m=r["m_mid"],
                    published_register="float64 midpoints of the "
                                       "study's certified fit")
        pr = run_one(tag, o, u, prov, gp)
        manifest.append(dict(tag=tag, products={k: v["verdict"]
                                                for k, v in pr.items()
                                                if k != "polish"}))
        print(tag, manifest[-1]["products"], flush=True)

    # worked example: published record is the certified 22-digit theta
    D = [1, 1, 2, 3, 5, 8]
    o = EtienneSSOracle(D, tag="ss_worked_example")
    p = ident.pilot()
    thf, If, _ = p["phase2_certify"].float_opt(D, o.KP, th0=7.0, I0=6.0)
    thp, _ = engine.newton_polish(o, [math.log(thf), math.log(If)])
    u = list(map(float, thp))
    mpm = GT.MPSingleSample(D)
    gp = gate_pack_ss(o, u, receipt_lnp=None, do_a5=True, mp_model=mpm)
    prov = dict(
        receipt=os.path.join(ETIENNE_DIR, "RESULT.md"),
        published="theta_hat = 7.047969433001057383852(+/-5e-22) "
                  "(R5, corrected Etienne worked example, D=(1,1,2,3,5,8))",
        candidate="self-polished center; the certified ball must contain "
                  "the published interval (uniqueness transfers the "
                  "statement to the published value)")
    pr = run_one("ss_worked_example", o, u, prov, gp)
    # containment check of the published interval
    from flint import arb, ctx
    old = ctx.prec
    try:
        ctx.prec = 256
        res = pr.get("polished") or pr["direct"]
        if res["verdict"].startswith("CERTIFIED"):
            th_ball = arb(u[0]).exp()
            pub = arb("[7.047969433001057383852 +/- 5e-22]")
            d = abs(float((th_ball - pub).mid()))
            note = dict(theta_at_center=str(th_ball),
                        published_interval="7.047969433001057383852+/-5e-22",
                        center_distance=d,
                        radius_lntheta=res.get("radius"),
                        contained=bool(d < res.get("radius", 0.0) * 7.0))
            json.dump(note, open(os.path.join(
                RECEIPTS_DIR, "ss_worked_example",
                "PUBLISHED_CONTAINMENT.json"), "w"), indent=1)
    finally:
        ctx.prec = old
    manifest.append(dict(tag="ss_worked_example",
                         products={k: v["verdict"] for k, v in pr.items()
                                   if k != "polish"}))
    print("ss_worked_example", manifest[-1]["products"], flush=True)

    # BCI dbh>=1cm, J=235360
    path = os.path.join(DATA_DIR, "packages/bci_dbh1cm_c1.txt")
    D = load_vec(path)
    r_path = os.path.join(PILOT_DIR, "BCI_235K_CERTIFIED.json")
    r = load_receipt(r_path)
    assert r["J"] == sum(D) and r["S"] == len(D)
    o = EtienneSSOracle(D, tag="ss_bci235k")
    u = u_from_theta_m(r["theta_mid"], r["m_mid"], r["J"])
    gp = gate_pack_ss(o, u, receipt_lnp=r["lnP"], do_a5=False)
    prov = dict(receipt=r_path, receipt_sha=ident.sha256(r_path),
                theta=r["theta_mid"], m=r["m_mid"],
                published_register="float64 midpoints of the study's "
                                   "certified fit (R6, J=235,360)")
    pr = run_one("ss_bci235k", o, u, prov, gp)
    manifest.append(dict(tag="ss_bci235k",
                         products={k: v["verdict"] for k, v in pr.items()
                                   if k != "polish"}))
    print("ss_bci235k", manifest[-1]["products"], flush=True)

    # six-forest pm rows (published table register, 4-6 printed digits)
    sf_path = os.path.join(PILOT_DIR, "SIXFOREST_PM_RF.json")
    sf = load_receipt(sf_path)
    for forest in SIXFOREST:
        tag = f"ss_sixforest_{forest}_pm"
        path = os.path.join(DATA_DIR, f"volkov2005/derived/{forest}_abundvec.txt")
        D = load_vec(path)
        row = sf[forest]
        assert row["J"] == sum(D) and row["S"] == len(D)
        o = EtienneSSOracle(D, tag=tag)
        u = u_from_theta_m(row["pm_theta"], row["pm_m"], row["J"])
        gp = gate_pack_ss(o, u, receipt_lnp=None, do_a5=False)
        gp["value_receipt_register"] = dict(
            register="receipt_float",
            note="receipt stores pm_LL float + truncated ball string; "
                 "float register compared",
            oracle_nll=gp["value_receipt_register"]["oracle_nll"],
            receipt_lnp_float=row["pm_LL"])
        prov = dict(receipt=sf_path, receipt_sha=ident.sha256(sf_path),
                    theta=row["pm_theta"], m=row["pm_m"],
                    published_register="printed table row (4-6 digits) — "
                                       "the EH2011-format published point")
        pr = run_one(tag, o, u, prov, gp)
        manifest.append(dict(tag=tag, products={k: v["verdict"]
                                                for k, v in pr.items()
                                                if k != "polish"}))
        print(tag, manifest[-1]["products"], flush=True)

    # equal-I multisample: Panama 3-plot + Panama-65
    ms_path = os.path.join(DATA_DIR, "packages/panama3_full.txt")
    ident.log_read(ms_path, "panama 3-plot abundance matrix")
    rows = [tuple(int(x) for x in line.split())
            for line in open(ms_path)]
    r_path = os.path.join(PILOT_DIR, "PHASE2_MS_EQI_PANAMA.json")
    r = load_receipt(r_path)
    o = EtienneEqIOracle(rows, tag="eqi_panama3")
    u = u_from_theta_I(r["theta_mid"], r["I_mid"])
    mpm = GT.MPEqI(rows)
    from flint import arb, ctx
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

    import csv
    p65_path = os.path.join(
        DATA_DIR, "panama65/derived/panama65_recent_abundance_matrix_dbh100.csv")
    ident.log_read(p65_path, "panama-65 abundance matrix (dbh>=100mm)")
    rows65 = []
    with open(p65_path) as f:
        rd = csv.reader(f)
        next(rd)
        for line in rd:
            rows65.append(tuple(int(x) for x in line[1:]))
    r_path = os.path.join(PILOT_DIR, "PANAMA65_EQI_CERTIFIED.json")
    r = load_receipt(r_path)
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
                                          "MANIFEST_TIER_A.json"), "w"),
              indent=1)
    ident.flush_read_log("tier_a")
    print("TIER A DONE", flush=True)


if __name__ == "__main__":
    main()
