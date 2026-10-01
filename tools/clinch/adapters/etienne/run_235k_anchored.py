"""235k rerun with the ANCHORED accumulator oracle (adapter-owned;
lnP_derivs2_anchored). The reference-route certificate refused on
enclosure width: at the 1e-9-sigma box the H_lnI,lnI enclosure has
relative width ~0.28 (measured: mid 3.803, rad 1.081 at lnI box width
2.2e-10) because the reference loop anchors its running product at A=S
(~2.3e5 ball multiplies before the mass of the A-distribution). The
anchored variant is the same mathematics re-anchored at a0 = E[A];
identity gates: agreement with the reference route at tight centers to
ball-radius level (receipted), plus the standard value/FD gates."""
import json
import math
import os
import sys
import time

sys.dont_write_bytecode = True

from . import DATA_DIR, PILOT_DIR, RECEIPTS_DIR, ident
from . import gates as GT
from .oracle_ss import EtienneSSOracle, lnP_derivs2_anchored
from .run_ss import run_one, u_from_theta_m

from flint import arb, ctx


def main():
    path = os.path.join(DATA_DIR, "packages/bci_dbh1cm_c1.txt")
    D = sorted(int(x) for x in open(path).read().split())
    ident.log_read(path, "BCI dbh>=1cm abundance vector (anchored rerun)")
    r_path = os.path.join(PILOT_DIR, "BCI_235K_CERTIFIED.json")
    ident.log_read(r_path, "235k fit receipt (anchored rerun)")
    r = json.load(open(r_path))
    o = EtienneSSOracle(D, tag="ss_bci235k_anchored", anchored=True)
    u = u_from_theta_m(r["theta_mid"], r["m_mid"], r["J"])
    # identity gate: anchored vs reference route at the tight center
    p2 = ident.pilot()["phase2_certify"]
    old = ctx.prec
    try:
        ctx.prec = 192
        th, I = arb(u[0]).exp(), arb(u[1]).exp()
        l1, g1, H1 = p2.lnP_derivs2(D, o.KP, th, I, 192)
        l2, g2, H2 = lnP_derivs2_anchored(D, o.KP, th, I, 192)
        route_gate = dict(
            register="anchored_vs_reference_identity",
            lnP_mid_diff=float((l1 - l2).mid()),
            g_mid_diffs=[float((g1[k] - g2[k]).mid()) for k in (0, 1)],
            H_mid_diffs=[[float((H1[i][j] - H2[i][j]).mid())
                          for j in (0, 1)] for i in (0, 1)],
            bar=1e-30)
        route_gate["ok"] = bool(
            abs(route_gate["lnP_mid_diff"]) <= 1e-30
            and all(abs(v) <= 1e-30 for v in route_gate["g_mid_diffs"])
            and all(abs(v) <= 1e-20 for row in route_gate["H_mid_diffs"]
                    for v in row))
        ctx.prec = 256
        x = ([[arb(u[1])]], [arb(u[0])])
        nll = o.nll(x)
    finally:
        ctx.prec = old
    gp = dict(value_receipt_register=GT.value_gate(nll, r["lnP"]),
              anchored_route_identity=route_gate,
              fd_hessian_register=GT.fd_hessian_gate(o, u))
    prov = dict(receipt=r_path, receipt_sha=ident.sha256(r_path),
                theta=r["theta_mid"], m=r["m_mid"],
                oracle_variant="lnP_derivs2_anchored (adapter-owned; "
                               "reference-route certificate refused on "
                               "enclosure width, receipt ss_bci235k/)",
                published_register="float64 midpoints of the study's "
                                   "certified fit (R6, J=235,360)")
    pr = run_one("ss_bci235k_anchored", o, u, prov, gp)
    print("ss_bci235k_anchored",
          {k: v["verdict"] for k, v in pr.items() if k != "polish"},
          flush=True)
    ident.flush_read_log("t235k_anchored")


if __name__ == "__main__":
    main()
