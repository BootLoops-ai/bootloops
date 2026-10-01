"""clinch.cert — CERT.json + the one-paragraph human statement."""
import json
import os
import time

__all__ = ["write_cert", "human_statement"]


def human_statement(res, model_tag):
    d = res.get("dims", {})
    v = res["verdict"]
    if v == "OUT-OF-SCOPE":
        return (f"OUT-OF-SCOPE. The {model_tag} candidate lives in a model "
                f"space this oracle has not ported and identity-gated: "
                f"{res['reason']} Certifying an unported space would be "
                f"untestable code; port and gate first.")
    if v.startswith("CERTIFIED"):
        kw = res["krawczyk"]; pd = res["pd"]
        pdm = pd.get("margins", {})
        rtxt = (f"per-coordinate radii in [{res.get('radius_min'):.3g}, "
                f"{res.get('radius'):.3g}] (curvature-scaled box)"
                if res.get("radius_min") not in (None, res.get("radius"))
                else f"radius {res['radius']:.3g}")
        corner = ""
        if v == "CERTIFIED-CORNER":
            ks = sorted((res.get("active_set") or {}).items())
            corner = (
                f" The candidate sits EXACTLY on {len(ks)} declared "
                f"parameter bound(s) ({', '.join(f'coord {k} at {s} bound'
                for k, (s, _) in ks)}); the statement is the BOUNDED "
                f"problem's: the reduced (free-coordinate) system "
                f"certifies as above and every active bound's KKT "
                f"multiplier sign is strictly certified over the box "
                f"(strict complementarity).")
        return (
            f"{v}. Inside the box of {rtxt} over the "
            f"{d.get('n_total')} certified coordinates around "
            f"the {model_tag} candidate optimum, a stationary point of the "
            f"penalized objective EXISTS, is UNIQUE, and every Hessian in "
            f"the box is positive definite — a certified strict local "
            f"minimum. Proof: hierarchical block-arrow Krawczyk "
            f"contraction (worst margin {kw['max_margin']:.3g} < 1 over "
            f"{d.get('n_blocks')} latent blocks + the {d.get('border')}-"
            f"dim border, interval Schur complement on the border) and "
            f"interval block-arrow Cholesky "
            f"(worst PD margin {min(pdm.values()):.3g} > 0), in arb "
            f"ball arithmetic at {res['prec_dps']} digits.{corner} This "
            f"is a LOCAL statement only; no global-optimality claim is "
            f"made.")
    return (
        f"REFUSED. The {model_tag} candidate could not be certified: "
        f"{res['reason']} (failing subspace: {res.get('failed')}). "
        f"A refusal names the failed contraction and its measured defect; "
        f"it is the honest outcome — a fat ball is a refusal, never a "
        f"shrug.")


def write_cert(res, out_dir, model_tag, gates=None, extra=None):
    """CERT.json + STATEMENT.txt in out_dir; returns the CERT dict."""
    os.makedirs(out_dir, exist_ok=True)
    cert = dict(
        tool="clinch",
        model=model_tag,
        stamp=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        verdict=res["verdict"],
        ball_center_sha16=res.get("theta_sha16"),
        radius=res.get("radius"),
        prec_bits=res.get("prec_bits"), prec_dps=res.get("prec_dps"),
        dims=res.get("dims"),
        contraction_margins=(res.get("krawczyk") or {}).get("margins_named"),
        newton_step_inf=(res.get("krawczyk") or {}).get("newton_step_inf"),
        pd_margins=(res.get("pd") or {}).get("margins_named")
        or (res.get("pd") or {}).get("margins"),
        failed=res.get("failed"),
        reason=res.get("reason"),
        oracle_identity_gates=gates,
        rungs=[dict(radius=r["radius"],
                    verdict=r["krawczyk"]["verdict"],
                    max_margin=r["krawczyk"].get("max_margin"),
                    failed=r["krawczyk"].get("failed_named",
                                             r["krawczyk"].get("failed")),
                    wall_s=r["krawczyk"].get("wall_s"))
               for r in res.get("rungs", [])],
        wall_s_total=res.get("wall_s_total"),
        statement=human_statement(res, model_tag))
    if extra:
        cert.update(extra)
    p = os.path.join(out_dir, "CERT.json")
    json.dump(cert, open(p + ".part", "w"), indent=1, default=float)
    os.replace(p + ".part", p)
    open(os.path.join(out_dir, "STATEMENT.txt"), "w").write(
        cert["statement"] + "\n")
    return cert
