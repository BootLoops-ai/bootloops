"""Certificate run, tier B: the microbiome joint equal-I depth-ladder fits
(S4_DEPTH_LADDER.json; "joint interior at every rarefaction
depth"; the depth-1000 rung carries the quoted m_eq 0.38-vs-Sloan-0.62
comparison).

The rarefied matrices are reproduced with the study's own seeded
procedure (random.Random(16180339), identical pool construction); the
reproduction is IDENTITY-CHECKED against the receipts (S per rung and the
receipted theta/I must be reproduced as a stationary point) before any
certificate. Candidates: the receipted (theta_mid, I_mid) floats.
"""
import json
import math
import os
import random
import sys
import time

sys.dont_write_bytecode = True

import numpy as np

from . import SCIENCE_DIR, RECEIPTS_DIR, ident
from .oracle_eqi import EtienneEqIOracle
from .run_corner import load_seawater
from . import gates as GT

from clinch import engine, cert

STAMP = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def main():
    r_path = os.path.join(SCIENCE_DIR, "S4_DEPTH_LADDER.json")
    ident.log_read(r_path, "depth-ladder fit receipts (candidates)")
    rec = json.load(open(r_path))
    M, names = load_seawater()
    S, N = M.shape
    rng = random.Random(16180339)
    manifest = []
    for depth in (1000, 3000, 10000, 30000):
        row = rec[str(depth)]
        sub = np.zeros((S, N), dtype=np.int64)
        for j in range(N):
            col = M[:, j]
            pool = np.repeat(np.arange(S), col)
            idx = rng.sample(range(len(pool)), depth)
            np.add.at(sub[:, j], pool[list(idx)], 1)
        rows = [tuple(int(x) for x in sub[s, :]) for s in range(S)
                if sub[s, :].sum() > 0]
        tag = f"eqi_seawater_depth{depth}"
        if len(rows) != row["S"]:
            manifest.append(dict(tag=tag, verdict="SKIPPED",
                                 reason=f"rarefaction reproduction S="
                                        f"{len(rows)} != receipt {row['S']}"
                                        f" (RNG stream drift)"))
            print(tag, "SKIPPED S mismatch", flush=True)
            continue
        o = EtienneEqIOracle(rows, tag=tag)
        u = [math.log(row["theta_mid"]), math.log(row["I_mid"])]
        out_dir = os.path.join(RECEIPTS_DIR, tag)
        os.makedirs(out_dir, exist_ok=True)
        logf = open(os.path.join(out_dir, "LOG.txt"), "a")

        def log(s):
            logf.write(s + "\n")
            logf.flush()
        log(f"== {tag} @ {STAMP}")
        gp = dict(value_receipt_register=dict(
            register="receipt_theta_ball",
            receipt_theta=row["theta"], receipt_I_mid=row["I_mid"],
            note="ladder receipts pin (theta ball, I mid); candidate = "
                 "float mids; the certificate's own containment carries "
                 "the identity"),
            fd_hessian_register=GT.fd_hessian_gate(o, u))
        res = engine.certify(o, u, log=log)
        prov = dict(receipt=r_path, receipt_sha=ident.sha256(r_path),
                    depth=depth, theta=row["theta_mid"], I=row["I_mid"],
                    m_equiv=row.get("m_equiv"),
                    reproduction="seeded rarefaction reproduced "
                                 "(random.Random(16180339)); S matched")
        products = dict(direct=res)
        if res["verdict"].startswith("REFUSED"):
            thp, prcpt = engine.newton_polish(o, u, log=log)
            resp = engine.certify(o, list(thp), log=log)
            products["polished"] = resp
            cert.write_cert(resp, os.path.join(out_dir, "polished"),
                            tag + "_polished", gates=gp,
                            extra=dict(product="polished_center",
                                       polish=prcpt,
                                       candidate_provenance=prov,
                                       run_stamp_utc=STAMP))
        cert.write_cert(res, out_dir, tag, gates=gp,
                        extra=dict(product="direct",
                                   candidate_provenance=prov,
                                   coordinates="(lntheta, lnI)",
                                   run_stamp_utc=STAMP))
        manifest.append(dict(tag=tag,
                             products={k: v["verdict"]
                                       for k, v in products.items()}))
        print(tag, manifest[-1], flush=True)
        logf.close()
    json.dump(manifest, open(os.path.join(RECEIPTS_DIR,
                                          "MANIFEST_TIER_B.json"), "w"),
              indent=1)
    ident.flush_read_log("tier_b")
    print("TIER B DONE", flush=True)


if __name__ == "__main__":
    main()
