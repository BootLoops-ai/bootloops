#!/usr/bin/env python3
"""
run_conipfv.py — orchestrator: g2' (CP0-CP6) + g3' (CP7-CP10) + transport
routing for coni-PFV cards.  Slots into the card battery next to run_pipe:
  python3 run_conipfv.py cards/<card>.json     one card, transcript + gates
  python3 run_conipfv.py --all cards/          suite sweep + scoreboard
Verdicts (exit 0 unless FAIL):
  NOT-CONI      PFFV card; g2/g3 of run_pipe own it (CP0 routing)
  CONI-PENDING  g2' + g3' PASS; transport awaits coni_op (fail-closed)
  CONI-READY    g2' + g3' PASS and transport inputs present
  CONI-FAIL     a CP gate failed (assertion text recorded)
Writes conipfv/out/<card>.gates.json per card and a scoreboard JSON in out/
for the cards given on the command line (a single-card call rewrites the
scoreboard with that card only); run_pipe outputs and reference data are
never touched.
"""
import sys, os, json, traceback

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.dirname(_HERE))

from family import load_card
from coni_curve import coni_pfv_curve
from coni_frobenius import g3_series
from coni_transport import transport_ready
from mpmath import mp, nstr


def _ser(x):
    from fractions import Fraction
    if isinstance(x, Fraction):
        return str(x)
    if isinstance(x, (list, tuple)):
        return [_ser(v) for v in x]
    if isinstance(x, dict):
        return {str(k): _ser(v) for k, v in x.items()}
    if type(x).__module__.startswith("mpmath"):
        return nstr(x, 20)
    return x


def run_card(path, outdir):
    card = load_card(path)
    rec = dict(card=card.name, path=os.path.basename(path))
    try:
        cur = coni_pfv_curve(card)
        if cur["routing"] == "NOT-CONI":
            rec.update(verdict="NOT-CONI", note=cur["note"])
            return rec
        rec["g2p"] = dict(PASS=True, q_cf=cur["q_cf"], n_cf=cur["n_cf"],
                          M_cf=str(cur["M_cf"]), K_cf=str(cur["K_cf"]),
                          lam=str(cur["lam"]), r=cur["r"], nu=cur["nu"],
                          KNiK=str(cur["KNiK"]), tadpole=cur["tadpole"],
                          bulk_aM_even=cur["bulk_aM_even"],
                          crosswalk=_ser(cur["crosswalk"]))
        g3 = g3_series(card, cur, dps=50)
        res = g3["racetrack"]
        rec["g3p"] = dict(
            PASS=True,
            S_ladder={str(N): str(v) for N, v in sorted(g3["ladder"]["S"].items())},
            t_vac=nstr(res["t"], 15), roots=res["joint"]["roots"],
            q_eff_times_r=res["joint"]["q_eff"],
            g_s=nstr(res["z_cf"]["g_s"], 10),
            zcf_Fterm=nstr(res["zmag_F"], 10),
            zcf_Qthroat=nstr(res["z_cf"]["zvev"], 10),
            vev_gap_efolds=nstr(res["vev_form_gap_efolds"], 6),
            zcf_floor=nstr(res["z_cf"]["floor"], 6),
            W0_est=nstr(res["W0"], 12), W0_floor=nstr(res["W0_floor"], 6),
            c_log=str(g3["frobenius"]["c_log"]),
            c_tau_eq_lam=str(g3["frobenius"]["c_tau"]),
            pin_checks=_ser(g3["pin_checks"]),
            label="ESTIMATE (pinned window); certified = transport legs 1-3")
        ready = transport_ready(card)
        rec.update(verdict="CONI-READY" if ready else "CONI-PENDING",
                   transport_ready=ready)
    except AssertionError as e:
        cw = card.raw.get("coni_basis_crosswalk") or {}
        if cw.get("dictionary_blocked_gates") or card.raw.get("blocker"):
            rec.update(verdict="CONI-BLOCKED", gate_fail=str(e),
                       note="card-documented structural blocker (frame "
                            "dictionary not closed); fail-closed as designed")
        else:
            rec.update(verdict="CONI-FAIL", gate_fail=str(e))
    except Exception as e:
        rec.update(verdict="CONI-FAIL", error=traceback.format_exc(limit=3))
    return rec


def main():
    outdir = os.path.join(_HERE, "out")
    os.makedirs(outdir, exist_ok=True)
    args = sys.argv[1:]
    paths = []
    if args and args[0] == "--all":
        d = args[1]
        paths = sorted(os.path.join(d, f) for f in os.listdir(d)
                       if f.endswith(".json"))
    else:
        paths = args
    board = []
    nfail = 0
    for p in paths:
        rec = run_card(p, outdir)
        with open(os.path.join(outdir, rec["card"] + ".gates.json"), "w") as f:
            json.dump(rec, f, indent=1)
        line = f"{rec['card']:<34} {rec['verdict']}"
        if rec["verdict"] == "CONI-PENDING":
            line += (f"  lam={rec['g2p']['lam']}  n_cf={rec['g2p']['n_cf']}"
                     f"  M_cf={rec['g2p']['M_cf']}  t={rec['g3p']['t_vac']}"
                     f"  W0_est={rec['g3p']['W0_est']}")
        elif rec["verdict"] == "CONI-BLOCKED":
            line += "  " + rec.get("gate_fail", "")[:80]
        elif rec["verdict"] == "CONI-FAIL":
            line += "  " + rec.get("gate_fail", rec.get("error", ""))[:110]
            nfail += 1
        print(line)
        board.append(rec)
    with open(os.path.join(outdir, "scoreboard.json"), "w") as f:
        json.dump(board, f, indent=1)
    sys.exit(1 if nfail else 0)


if __name__ == "__main__":
    main()
