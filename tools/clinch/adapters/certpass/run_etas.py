#!/usr/bin/env python3
"""run_etas — certificate pass for the ETAS reference-fit optima
(T1: q3_hi/q3_lo live-chain refits; T4: q0_hi/q0_lo retro-seal refits).

Certified objective: the incomplete-data (observed) ETAS NLL of record
(oracle_etas module doc). The banked candidates are EM fixed points
converged to param-diff < 1e-3 — the DIRECT certificate at the banked
vector is expected to REFUSE (fat-ball law); the polished-center product
+ the per-coordinate shift ledger vs the banked vector is the
deliverable. Banked numbers are never modified.

Reads the originating study's un-shipped banked fit artifacts (stack/)
through the CERTPASS_REFERENCE_DIR environment variable; refuses loudly when
it is unset.

Usage: python3 run_etas.py [q3_hi] [q3_lo] [q0_hi] [q0_lo]
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
HR = pkg.reference_root()   # un-shipped study tree (stack/); env-gated

import clinch                                      # noqa
from clinch import engine, cert                    # noqa
from oracle_etas import EtasOracle, Part           # noqa
from twin_etas import TwinEtas                     # noqa
import gates_etas as GE                            # noqa

RECEIPTS = pkg.RECEIPTS_DIR
DATA = pkg.DATA_DIR
BUDGET_S = 10800.0
NW_FULL = int(os.environ.get("CERTPASS_NW", "10"))

FIT_ART = {
    "q3_hi": "stack/runner/ref_fit_hi/e2_params.json",
    "q3_lo": "stack/runner/ref_fit_lo/e2_params.json",
    "q0_hi": "stack/q0/ref_fit_hi/e2_params.json",
    "q0_lo": "stack/q0/ref_fit_lo/e2_params.json",
}

# the FIT'S OWN parameter bounds (etas.inversion.RANGES, iota dropped) —
# the M-step L-BFGS-B optimized inside these; the bounded-problem
# certificate uses the same box
RANGES_LO = np.array([-10.0, -20.0, 0.01, -8.0, -0.99, 0.01, -4.0,
                      -1.0, 0.01])
RANGES_HI = np.array([0.0, 10.0, 20.0, 0.0, 1.0, 12.26, 3.0, 5.0, 5.0])
NAMES9 = ["log10_mu", "log10_k0", "a", "log10_c", "omega", "log10_tau",
          "log10_d", "gamma", "rho"]


def sha16(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()[:16]


def log_for(tag):
    def log(s):
        print(f"[{tag}] {s}", flush=True)
    return log


def theta9(dump):
    """banked 10-vector (iota NaN) -> 9-vector of record."""
    th = np.asarray(dump["theta"], float)
    return np.r_[th[0], th[2:]]


def run_target(tag):
    log = log_for(tag)
    t0 = time.time()
    z = dict(np.load(os.path.join(DATA, f"etas_{tag}.npz")))
    th = theta9(z)
    dump_rec = json.load(open(os.path.join(
        HERE, f"logs_dump_{tag.replace('_', '')}.out")))
    part = Part(blocks=[list(range(1, 9))], block_names=["triggering"],
                border_idx=[0])
    log(f"n_pairs={len(z['pair_t'])} n_t={len(z['trg_pbg'])} "
        f"n_s={len(z['src_m'])}")
    oracle = EtasOracle(z, part, tag, nworkers=NW_FULL)
    twin = TwinEtas(z)
    sdata, srec = GE.slice_data(z)
    oracle_s = EtasOracle(sdata, part, tag + "_slice", nworkers=2)
    twin_s = TwinEtas(sdata)
    log(f"slice: {srec}")
    log("building mp50 slice transliteration")
    mp_obj = GE.MpEtas(sdata)
    e2 = FIT_ART[tag]
    summary = dict(tag=f"etas_{tag}", p=9,
                   n_pairs=int(len(z["pair_t"])),
                   n_targets=int(len(z["trg_pbg"])),
                   n_sources=int(len(z["src_m"])),
                   banked_fit=e2, banked_fit_sha16=sha16(f"{HR}/{e2}"),
                   theta=[float(v) for v in th],
                   dump_reproduction=dump_rec, slice_receipt=srec,
                   objective="incomplete-data (observed) ETAS NLL of "
                             "the pinned Mizrahi EM (oracle_etas doc)")
    # mu self-consistency receipt (EM: mu = n_hat/(A*T))
    mu = 10.0 ** th[0]
    summary["mu_nhat_selfconsistency"] = dict(
        mu_AT=float(mu * twin.A * twin.T),
        n_hat_banked=float(z["n_hat_banked"]),
        rel_dev=float(abs(mu * twin.A * twin.T - float(z["n_hat"]))
                      / float(z["n_hat"])))
    log("gate: EM identities (twin, full data)")
    ge = GE.gate_em_identities(twin, th, z)
    log(f"  pbg {ge['dev_pbg']:.3g} lhat {ge['dev_lhat']:.3g} "
        f"G_rel {ge['dev_G_rel']:.3g} nhat {ge['dev_nhat']:.3g} "
        f"PASS={ge['PASS']}")
    log("gate: value")
    gv = GE.gate_value(oracle, twin, oracle_s, mp_obj, twin_s, th)
    log(f"  arb-twin(rel) {gv['dev_arb_vs_twin_rel']:.3g} "
        f"slice arb-mp50 {gv['dev_slice_arb_vs_mp50_rel']:.3g} "
        f"PASS={gv['PASS']}")
    log("gate: gradient (mp50 slice adjudication, all 9 coords)")
    gg, _ = GE.gate_gradient_adjudicated(oracle_s, mp_obj, twin, th,
                                         oracle_full=oracle)
    log(f"  full twin: analytic {gg['full_arb_vs_twin_analytic_cols']:.3g}"
        f" fd-cols {gg['full_arb_vs_twin_fd_cols']:.3g} "
        f"mp50 {gg['adjudicated_worst']:.3g} PASS={gg['PASS']}")
    log("gate: fd-hessian (slice, all 9 cols)")
    gh = GE.gate_fd_hessian(oracle_s, th)
    log(f"  {gh['metric_max']:.3g} PASS={gh['PASS']}")
    gts = dict(em_identities=ge, value=gv, gradient_adjudicated=gg,
               fd_hessian=gh)
    summary["gates"] = gts
    if not (ge["PASS"] and gv["PASS"] and gg["PASS"] and gh["PASS"]):
        summary["verdict"] = "GATES-FAILED"
        log("GATES FAILED — no certificate attempted (fail-closed)")
        oracle.close(); oracle_s.close()
        return summary
    log("certify: direct (banked EM point AS GIVEN)")
    res = engine.certify(oracle, th, budget_s=BUDGET_S, log=log)
    out_dir = os.path.join(RECEIPTS, f"etas_{tag}")
    extra = dict(candidate_provenance="banked-full-vector (EM fixed "
                                      "point, param-diff<1e-3 stop)",
                 fit_artifact=e2,
                 em_register="direct refusal expected: the EM stop "
                             "tolerance leaves a nonzero observed-NLL "
                             "gradient; the polished-center product is "
                             "the certificate deliverable")
    cert.write_cert(res, out_dir, f"etas_{tag}", gates={
        k: {kk: vv for kk, vv in v.items() if not isinstance(vv, dict)}
        for k, v in gts.items()}, extra=extra)
    summary["direct"] = dict(verdict=res["verdict"],
                             failed=res.get("failed"),
                             newton_step_scaled=res.get(
                                 "newton_step_scaled_inf"),
                             max_margin=(res.get("krawczyk") or {})
                             .get("max_margin"),
                             wall_s=res.get("wall_s_total"))
    log(f"direct: {res['verdict']} ({res.get('failed')})")
    if res["verdict"] == "REFUSED":
        if os.environ.get("CERTPASS_PIN", "") == "ltau_hi":
            # corner-first arc (hi targets): the taper-ridge mechanism is
            # receipted on the lo twins' full free-walk diagnostics; here
            # the corner is proven by the certificate's own KKT leg
            # (strict sign of the full gradient at the bound) — the long
            # free walk is skipped for compute honesty, not rigor.
            th0 = th.copy()
            pins = {5: float(RANGES_HI[5])}
            pre = [dict(note="free walk skipped (CERTPASS_PIN=ltau_hi)")]
        else:
            log("hybrid damped Newton pre-polish (analytic pairs + "
                "arb-mid sources), FREE")
            th0, pre = twin.hybrid_polish(th, oracle, log=log)
            pins = {}
            for k in range(9):
                if th0[k] >= RANGES_HI[k]:
                    pins[k] = float(RANGES_HI[k])
                elif th0[k] <= RANGES_LO[k]:
                    pins[k] = float(RANGES_LO[k])
        summary["prepolish_free"] = pre[-3:]
        if pins:
            # free walk left the FIT'S OWN bounds — the bounded-problem
            # optimum is a corner; free-walk diagnostic receipted
            _, gfree, Hfree = (None, None, None)
            sv, sg, sH = oracle.sources_f64(th0)
            pv, pg, pH = twin.pairs_val_grad_hess(th0)
            ev = np.linalg.eigvalsh(pH + sH)
            free_diag = dict(
                endpoint={NAMES9[a]: float(th0[a]) for a in range(9)},
                nll_twin_at_endpoint=float(pv + sv),
                hess_eigs=[float(v) for v in ev],
                pinned={NAMES9[k]: v for k, v in pins.items()},
                note=("corner-first arc (CERTPASS_PIN): mechanism "
                      "receipted on the lo twins' free-walk "
                      "diagnostics; the corner is proven by this "
                      "certificate's KKT leg"
                      if os.environ.get("CERTPASS_PIN") else
                      "free polish crossed the fit's own L-BFGS-B "
                      "bounds (etas.inversion.RANGES); the taper "
                      "direction is a noncompact likelihood ridge — "
                      "the bounded-problem optimum is a CORNER"))
            log(f"free walk crossed bounds: {free_diag['pinned']} — "
                f"corner product")
            log("hybrid pinned polish")
            thc, prec_ = twin.hybrid_polish(th0, oracle, log=log,
                                            pinned=pins)
            from clinch.mask import MaskedOracle
            masked = MaskedOracle(oracle, pins)
            log("arb newton_polish (masked/pinned)")
            thp, prcpt = engine.newton_polish(masked, thc, log=log)
            for k, v in pins.items():
                thp[k] = v
            log("certify: corner (bounds = the fit's RANGES)")
            resp = engine.certify(oracle, thp,
                                  bounds=(RANGES_LO, RANGES_HI),
                                  budget_s=BUDGET_S, log=log)
            product = f"etas_{tag}_corner"
        else:
            free_diag = None
            log("arb newton_polish")
            thp, prcpt = engine.newton_polish(oracle, th0, log=log)
            log("certify: polished center")
            resp = engine.certify(oracle, thp, budget_s=BUDGET_S,
                                  log=log)
            product = f"etas_{tag}_polished"
        shift = dict(
            distance_inf=float(np.max(np.abs(thp - th))),
            per_coord={NAMES9[a]: dict(banked=float(th[a]),
                                       polished=float(thp[a]),
                                       shift=float(thp[a] - th[a]))
                       for a in range(9)},
            nll_twin_at_banked=float(twin.value(th)),
            nll_twin_at_polished=float(twin.value(thp)))
        shift["nll_drop"] = (shift["nll_twin_at_banked"]
                             - shift["nll_twin_at_polished"])
        extra_p = dict(extra, polish=prcpt, shift_vs_banked=shift,
                       free_walk_diagnostic=free_diag,
                       polished_center=[float(v) for v in thp])
        cert.write_cert(resp, os.path.join(RECEIPTS, product),
                        product, gates=None, extra=extra_p)
        summary["polished"] = dict(
            product=product,
            verdict=resp["verdict"], failed=resp.get("failed"),
            max_margin=(resp.get("krawczyk") or {}).get("max_margin"),
            pd_min=(min(((resp.get("pd") or {}).get("margins") or
                         {"": float("nan")}).values())
                    if resp.get("pd") else None),
            kkt=(resp.get("kkt") or {}).get("margins"),
            distance_inf=shift["distance_inf"],
            nll_drop=shift["nll_drop"], wall_s=resp.get("wall_s_total"))
        log(f"{product}: {resp['verdict']} "
            f"(dist {shift['distance_inf']:.3g}, nll drop "
            f"{shift['nll_drop']:.4g})")
    oracle.close(); oracle_s.close()
    summary["wall_s"] = round(time.time() - t0, 1)
    return summary


def main(args):
    tags = args if args else ["q3_lo", "q0_lo", "q3_hi", "q0_hi"]
    os.makedirs(RECEIPTS, exist_ok=True)
    mpath = os.path.join(RECEIPTS, f"RECEIPTS_etas_{'_'.join(tags)}.json")
    master = dict(stamp=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                  clinch_version=clinch.__version__,
                  prec_bits=engine.PREC_CERT, targets={})
    for tag in tags:
        s = run_target(tag)
        master["targets"][s["tag"]] = s
        json.dump(master, open(mpath + ".part", "w"), indent=1,
                  default=float)
        os.replace(mpath + ".part", mpath)
        print(f"[{tag}] DONE in {s.get('wall_s')}s", flush=True)
    print("ETAS PASS DONE", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
