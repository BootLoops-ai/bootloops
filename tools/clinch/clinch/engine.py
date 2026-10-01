"""clinch.engine — the certificate engine.

Consumes baller.certify.block_krawczyk BY IDENTITY (link, never copy):
the structured Krawczyk existence+uniqueness leg and the block-arrow PD
leg both live in BALLER; CLINCH owns the oracle, the radius policy, and
the products. An oracle that declares kkt_multiplier_idx (the
extended-KKT centered spaces, oracle_v36) routes the PD leg to
clinch.kkt_pd.pd_reduced_congruence — the reduced-Hessian
signed-congruence inertia certificate — since the extended Jacobian is
indefinite by construction and the plain PD target does not apply.

RADIUS POLICY (all constants CODE-PINNED, none user-tunable):
  PURE CURVATURE SCALING:  r_k = r_unit * s_k,
       s_k = 1/sqrt(max(mid H_kk, DIAG_FLOOR)) — a UNIFORM box in the
       Newton metric. Measured on fit13, the two wrong geometries both
       fail: a raw-uniform box suffers Bauer-Skeel width amplification
       ~cond(H) (~6 decades of curvature; margin 7.7e6), and mixing
       per-coordinate Newton-step floors into the radii creates radius
       RATIOS r_j/r_k spanning decades that amplify the cross-block terms
       by 1/r_k (margin 8.9e3 with the own-row width term at 1.0 exactly).
       Only the equilibrated geometry keeps row sums ~ r_unit. Scales are
       HEURISTIC input (box geometry only); rigor never depends on them.
  r_unit0 = STEP_COVER * ||scaled Newton step||_inf ("the optimizer's
       final step scale" in the Newton metric, measured from the oracle's
       own midpoint structured Schur solve): the box must COVER the
       residual step to the true optimum, so r_unit can never usefully go
       below the scaled residual — the ladder only climbs.
  ladder: r_unit0 * LADDER[k] (increasing), first CERTIFIED rung wins
       (the PD leg then runs on the SAME H(box) enclosure). All rungs
       refused -> REFUSED with the best rung's named margins.
  FAT-BALL LAW (spec: "a fat ball is a refusal, never a shrug"): a
       certificate is only issued when every radius is at most
       FAT_SIGMA_MAX curvature-sigmas (r_k <= FAT_SIGMA_MAX * s_k). A
       candidate whose scaled Newton residual already forces a fatter box
       (STEP_COVER * ||c/s||_inf + r_unit > FAT_SIGMA_MAX) REFUSES with
       that finding — a big certified ball around a poor candidate would
       be a true but WEAK statement, and CLINCH does not emit weak
       statements (the corner-without-declared-bounds control measured
       exactly this failure shape).
  WALL_BUDGET_S: if the next rung would start after the budget, stop and
       refuse honestly ("budget exhausted") — checkpoint-and-report, never
       silent grinding. No third outcome: the result is CERTIFIED or a
       named REFUSAL either way.
"""
import os
import sys
import time

import numpy as np

# baller lives beside this package in the tools tree (env-overridable)
BALLER_DIR = os.environ.get("CLINCH_BALLER_DIR") or os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "baller")

PREC_CERT = 192          # bits (~57 decimal digits) — pinned
STEP_COVER = 2.0         # pinned Newton-step coverage factor
R_UNIT_MIN = 1e-9        # pinned floor on r_unit (sigma units)
R_FLOOR = 1e-13          # pinned
DIAG_FLOOR = 1e-6        # pinned curvature floor for the radius scaling
FAT_SIGMA_MAX = 1.0      # pinned fat-ball law: max radius in sigma units
LADDER = (1.0, 2.0, 4.0, 8.0)        # pinned multipliers (climb only)
WALL_BUDGET_S = 1800.0   # pinned checkpoint bar (the spec's 30 min)
# accepted PD-leg verdicts: plain block-arrow PD, or the reduced-Hessian
# inertia certificate (extended-KKT centered spaces, clinch.kkt_pd)
_PD_OK = ("PD_CERTIFIED", "PD_REDUCED_CERTIFIED")

__all__ = ["certify", "certify_guarded", "PREC_CERT", "STEP_COVER",
           "R_UNIT0", "FAT_SIGMA_MAX", "LADDER", "WALL_BUDGET_S"]


def _baller():
    if BALLER_DIR not in sys.path:
        sys.path.insert(0, BALLER_DIR)
    import baller
    baller.verify(quiet=True)                 # fail-closed, run first (law)
    from baller.certify import block_krawczyk as bk
    return bk


def _newton_step(bk, Ft, Ht):
    """(step_blocks, step_border): float approximate Newton step at the
    candidate (midpoint structured Schur solve; heuristic seed only)."""
    (Fz, Fg), (D, B, G) = Ft, Ht
    Y, P = [], []
    S = bk._mid_np(G)
    for i in range(len(D)):
        Dm = bk._mid_np(D[i]); Bm = bk._mid_np(B[i])
        Yi = np.linalg.inv(Dm)
        Pi = Yi @ Bm
        S -= Bm.T @ Pi
        Y.append(Yi); P.append(Pi)
    Ys = np.linalg.inv(S)
    fz = [np.array([float(v.mid()) for v in blk]) for blk in Fz]
    fg = np.array([float(v.mid()) for v in Fg])
    acc = fg.copy()
    for j in range(len(D)):
        acc -= P[j].T @ fz[j]
    cg = Ys @ acc
    cz = [Y[i] @ fz[i] - P[i] @ cg for i in range(len(D))]
    return cz, cg


def _radius_scales(bk, Ht):
    """Curvature scales s_i = 1/sqrt(max(mid H_ii, DIAG_FLOOR)) per
    coordinate, block-structured (heuristic box geometry only)."""
    (D, B, G) = Ht
    sz = []
    for i in range(len(D)):
        d = np.maximum(np.diag(bk._mid_np(D[i])), DIAG_FLOOR)
        sz.append(1.0 / np.sqrt(d))
    sg = 1.0 / np.sqrt(np.maximum(np.diag(bk._mid_np(G)), DIAG_FLOOR))
    return sz, sg


def certify(oracle, theta, bounds=None, prec=PREC_CERT,
            budget_s=WALL_BUDGET_S, log=None):
    """Run the full CLINCH certificate on candidate theta (floats).

    VERDICT SET (A2):
      CERTIFIED-INTERIOR  no active bound; unique stationary point + PD
                          Hessian in the box (certified strict local min).
                          For an extended-KKT oracle (kkt_multiplier_idx
                          declared): unique KKT point + reduced-Hessian
                          PD via the signed-congruence inertia leg
                          (clinch.kkt_pd) — a certified strict local
                          minimum of the equality-constrained problem.
      CERTIFIED-CORNER    coordinates EXACTLY at a declared bound: the
                          reduced (free-coordinate) system certifies, the
                          reduced Hessian is PD, and every active bound's
                          KKT multiplier sign is STRICTLY certified over
                          the box (strict complementarity) — a certified
                          strict local minimum of the BOUNDED problem.
                          (an AT-CORNER audit flag alone is out of scope;
                          CLINCH adjudicates the KKT.)
      REFUSED(named)      the failed contraction/PD/KKT leg named with its
                          measured defect.
      OUT-OF-SCOPE        unported model space (use certify_guarded).
    bounds: None (unconstrained fit — the home v3.1 case) or (lo, hi)
    float arrays; a coordinate is ACTIVE iff theta_k == bound_k EXACTLY
    (published optima are clipped exactly; near-misses stay interior)."""
    log = log or (lambda s: None)
    bk = _baller()
    from .oracle_v31 import OracleRefusal
    active = {}
    if bounds is not None:
        lo, hi = bounds
        th = np.asarray(theta, dtype=float)
        for k in range(len(th)):
            if lo is not None and np.isfinite(lo[k]) and th[k] == lo[k]:
                active[k] = ("lower", float(lo[k]))
            elif hi is not None and np.isfinite(hi[k]) and th[k] == hi[k]:
                active[k] = ("upper", float(hi[k]))
    if active:
        from .mask import MaskedOracle
        masked = MaskedOracle(oracle, {k: v for k, (_, v) in
                                       active.items()})
        res = _certify_system(masked, theta, prec, budget_s, log, bk,
                              OracleRefusal, kkt_active=active,
                              base_oracle=oracle)
        return res
    return _certify_system(oracle, theta, prec, budget_s, log, bk,
                           OracleRefusal, kkt_active=None,
                           base_oracle=oracle)


def certify_guarded(make_oracle, theta, **kw):
    """certify() with the OUT-OF-SCOPE verdict: make_oracle() constructing
    an oracle for an unported space (check_space refusal) returns the
    typed OUT-OF-SCOPE verdict instead of raising."""
    from .oracle_v31 import OracleRefusal
    try:
        oracle = make_oracle()
    except OracleRefusal as e:
        if e.where == "check_space":
            return dict(verdict="OUT-OF-SCOPE", failed="check_space",
                        reason=str(e),
                        note="unported model space — certifying it would "
                             "be untestable code; port + gate first")
        raise
    return certify(oracle, theta, **kw)


def _certify_system(oracle, theta, prec, budget_s, log, bk, OracleRefusal,
                    kkt_active=None, base_oracle=None):
    part = oracle.part
    interior_verdict = ("CERTIFIED-CORNER" if kkt_active
                        else "CERTIFIED-INTERIOR")
    names = ["border"] + list(part.block_names)

    def named(margins):
        out = {}
        for k, v in margins.items():
            if k == "border":
                out["border"] = v
            else:
                out[part.block_names[int(k.split("_")[1])]] = v
        return out

    t_all = time.time()
    border_idx = getattr(part, "border_idx", None)
    if border_idx is None:
        border_idx = list(range(oracle.dims[1]))
    center = ([[float(theta[t]) for t in idx] for idx in part.blocks],
              [float(theta[t]) for t in border_idx])
    rungs = []
    try:
        t0 = time.time()
        Ft = None
        with __import__("baller").hygiene.ctx_guard(prec=prec):
            box0 = bk.box_around(center, 0.0)
            Ft = oracle.F(box0)
            Ht = oracle.H(box0)
        wall_center = time.time() - t0
        log(f"center eval {wall_center:.1f}s")
        cz, cg = _newton_step(bk, Ft, Ht)
        sz, sg = _radius_scales(bk, Ht)
        s_inf = max([float(np.max(np.abs(c))) for c in cz if len(c)]
                    + [float(np.max(np.abs(cg)))])
        scaled = [float(np.max(np.abs(c) / s)) for c, s in zip(cz, sz)
                  if len(c)] + [float(np.max(np.abs(cg) / sg))]
        s_inf_scaled = max(scaled)
        log(f"newton step inf ~ {s_inf:.3g} (scaled {s_inf_scaled:.3g}); "
            f"curvature scales "
            f"[{min(np.min(s) for s in sz + [sg]):.3g}, "
            f"{max(np.max(s) for s in sz + [sg]):.3g}]")
        ru0 = max(STEP_COVER * s_inf_scaled, R_UNIT_MIN)
        if ru0 > FAT_SIGMA_MAX:
            res = dict(
                verdict="REFUSED", failed="fat-ball",
                newton_step_inf_est=s_inf,
                newton_step_scaled_inf=s_inf_scaled, rungs=rungs,
                reason=(f"fat-ball law: the candidate's scaled Newton "
                        f"residual ({s_inf_scaled:.3g} curvature-sigma) "
                        f"already forces radii beyond "
                        f"{FAT_SIGMA_MAX:g} sigma — a certificate here "
                        f"would be a true but WEAK statement about a "
                        f"poorly-converged candidate; polish further or "
                        f"declare the active bounds"))
            return _finish(res, oracle, theta, prec, t_all)
        log(f"r_unit0 = {ru0:.3g} sigma (step cover)")
        best = None
        for mult in LADDER:
            if time.time() - t_all > budget_s:
                reason = (f"wall budget {budget_s:.0f}s exhausted after "
                          f"{len(rungs)} rungs — honest refusal, rerun "
                          f"with a larger budget to continue the ladder")
                res = dict(verdict="REFUSED", failed="budget",
                           reason=reason, rungs=rungs,
                           newton_step_inf_est=s_inf)
                if best is not None:
                    res["best_rung"] = best
                return _finish(res, oracle, theta, prec, t_all)
            ru = ru0 * mult
            if ru > FAT_SIGMA_MAX:
                log(f"rung r_unit={ru:.3g}: skipped (fat-ball law)")
                continue
            rs = ru
            r = ([[ru * float(v) for v in s] for s in sz],
                 [ru * float(v) for v in sg])
            t1 = time.time()
            try:
                with __import__("baller").hygiene.ctx_guard(prec=prec):
                    box = bk.box_around(center, r)
                    Hb = oracle.H(box)
            except OracleRefusal as e:
                # a clip guard tripping at THIS radius is this rung's
                # refusal, not the ladder's death — record and continue
                rungs.append(dict(radius=rs, krawczyk=dict(
                    verdict="REFUSED", failed=f"oracle:{e.where}",
                    reason=str(e), margins={}, margins_named={},
                    wall_s=round(time.time() - t1, 1))))
                log(f"rung r_unit={ru:.3g}: REFUSED oracle:{e.where}")
                continue
            fh = (Ft, Hb)
            kw = bk.block_krawczyk(oracle, center, r, prec, fh=fh)
            kw["radius_scale"] = rs
            kw["wall_s"] = round(time.time() - t1, 1)
            kw["margins_named"] = named(kw.get("margins", {}))
            if kw.get("failed") and kw["failed"].startswith("block_"):
                kw["failed_named"] = part.block_names[
                    int(kw["failed"].split("_")[1])]
            rungs.append(dict(radius=rs, krawczyk=kw))
            log(f"rung r_unit={rs:.3g}: {kw['verdict']} "
                f"max_margin={kw.get('max_margin')} "
                f"failed={kw.get('failed_named', kw.get('failed'))}")
            if kw["verdict"] == "CERTIFIED":
                t2 = time.time()
                mult_idx = getattr(oracle, "kkt_multiplier_idx", None)
                if mult_idx:
                    from .kkt_pd import pd_reduced_congruence
                    pd = pd_reduced_congruence(*Hb, mult_idx=mult_idx,
                                               prec=prec)
                else:
                    pd = bk.pd_block_arrow(*Hb, prec)
                pd["wall_s"] = round(time.time() - t2, 1)
                if "margins" in pd and pd.get("verdict") in _PD_OK:
                    pd["margins_named"] = named(
                        {k.replace("border_schur", "border"): v
                         for k, v in pd["margins"].items()})
                rungs[-1]["pd"] = pd
                if pd["verdict"] in _PD_OK:
                    kkt = None
                    if kkt_active:
                        kkt = _kkt_leg(bk, oracle, base_oracle, box,
                                       kkt_active, prec)
                        rungs[-1]["kkt"] = kkt
                        if not kkt["ok"]:
                            res = dict(verdict="REFUSED",
                                       failed=kkt["failed"],
                                       reason=kkt["reason"],
                                       krawczyk=kw, pd=pd, kkt=kkt,
                                       rungs=rungs,
                                       newton_step_inf_est=s_inf)
                            return _finish(res, oracle, theta, prec, t_all)
                    res = dict(verdict=interior_verdict,
                               radius=kw["radius"],
                               radius_min=kw["radius_min"],
                               radius_scale=rs,
                               krawczyk=kw, pd=pd, kkt=kkt, rungs=rungs,
                               failed=None,
                               active_set=(None if not kkt_active else
                                           {int(k): v for k, v in
                                            kkt_active.items()}),
                               newton_step_inf_est=s_inf,
                               reason=("reduced system certified + PD + "
                                       "strict KKT multiplier signs: "
                                       "certified strict local minimum "
                                       "of the bounded problem (corner)"
                                       if kkt_active else
                                       "unique KKT point in the ball + "
                                       "reduced Hessian PD throughout "
                                       "(signed-congruence inertia "
                                       "certificate): certified strict "
                                       "local minimum of the equality-"
                                       "constrained problem"
                                       if mult_idx else
                                       "unique stationary point in the "
                                       "ball + PD Hessian throughout: "
                                       "certified strict local minimum"))
                    return _finish(res, oracle, theta, prec, t_all)
                res = dict(verdict="REFUSED", failed=pd["failed"],
                           reason=pd["reason"], krawczyk=kw, pd=pd,
                           rungs=rungs, newton_step_inf_est=s_inf)
                return _finish(res, oracle, theta, prec, t_all)
            if best is None or kw.get("max_margin", float("inf")) < \
                    best["krawczyk"].get("max_margin", float("inf")):
                best = rungs[-1]
        if best is None:
            # ALL rungs failed before any margin could be measured (every
            # rung raised an oracle refusal, or every rung was skipped by
            # the fat-ball law) — a TYPED refusal, never a crash (the
            # etienne-deployment F6 finding: best=None dereference).
            classes = [r["krawczyk"].get("failed") for r in rungs]
            res = dict(
                verdict="REFUSED", failed="all-rungs-failed",
                reason=("every radius rung failed before a contraction "
                        "margin could be measured: "
                        + (f"per-rung failure classes {classes}"
                           if classes else
                           "all rungs skipped by the fat-ball law")),
                rungs=rungs, newton_step_inf_est=s_inf)
            return _finish(res, oracle, theta, prec, t_all)
        bkw = best["krawczyk"]
        res = dict(verdict="REFUSED",
                   failed=bkw.get("failed_named", bkw.get("failed")),
                   reason=(f"no ladder rung contracted; best margin "
                           f"{bkw.get('max_margin', float('nan')):.4g} at "
                           f"r_scale={best['radius']:.3g} on "
                           f"{bkw.get('failed_named', bkw.get('failed'))}"
                           + (f" ({bkw.get('reason')})"
                              if bkw.get("max_margin") is None else "")),
                   rungs=rungs, best_rung=best, newton_step_inf_est=s_inf)
        return _finish(res, oracle, theta, prec, t_all)
    except OracleRefusal as e:
        res = dict(verdict="REFUSED", failed=f"oracle:{e.where}",
                   reason=str(e), rungs=rungs)
        return _finish(res, oracle, theta, prec, t_all)


POLISH_MAXIT = 12        # pinned
POLISH_TOL_SIGMA = 1e-7  # pinned: stop when scaled step < this


def newton_polish(oracle, theta, prec=PREC_CERT, log=None):
    """Float structured-Newton polish of the candidate (midpoint Schur
    solves on the arb oracle's midpoints). Returns (theta_polished,
    receipt). NOT part of the direct certificate — the polished-center
    certificate is a SEPARATE, labeled product (certify at the polished
    cluster center, receipt the distance to the published point)."""
    log = log or (lambda s: None)
    bk = _baller()
    part = oracle.part
    border_idx = getattr(part, "border_idx", None)
    if border_idx is None:
        border_idx = list(range(oracle.dims[1]))
    th = np.array(theta, dtype=float)
    hist = []
    for it in range(POLISH_MAXIT):
        center = ([[float(th[t]) for t in idx] for idx in part.blocks],
                  [float(th[t]) for t in border_idx])
        with __import__("baller").hygiene.ctx_guard(prec=prec):
            box0 = bk.box_around(center, 0.0)
            Ft = oracle.F(box0)
            Ht = oracle.H(box0)
        cz, cg = _newton_step(bk, Ft, Ht)
        sz, sg = _radius_scales(bk, Ht)
        scaled = max([float(np.max(np.abs(c) / s))
                      for c, s in zip(cz, sz) if len(c)]
                     + [float(np.max(np.abs(cg) / sg))])
        raw = max([float(np.max(np.abs(c))) for c in cz if len(c)]
                  + [float(np.max(np.abs(cg)))])
        hist.append(dict(it=it, step_inf=raw, step_scaled_inf=scaled))
        log(f"polish it{it}: step {raw:.3g} (scaled {scaled:.3g})")
        if scaled <= POLISH_TOL_SIGMA:
            break
        for bi, idx in enumerate(part.blocks):
            th[idx] -= cz[bi]
        th[np.asarray(border_idx)] -= cg
    d_inf = float(np.max(np.abs(th - np.asarray(theta, dtype=float))))
    return th, dict(iterations=hist, distance_inf=d_inf,
                    tol_sigma=POLISH_TOL_SIGMA,
                    note="float64 structured Newton on the oracle "
                         "midpoints; heuristic only — the certificate "
                         "at the polished center carries the rigor")


def _kkt_leg(bk, masked, base_oracle, box, active, prec):
    """Certified KKT multiplier signs at a corner: the FULL-system gradient
    over the certified reduced box (active coords pinned exact) must be
    STRICTLY > 0 at a lower-active coordinate (< 0 at upper-active) for
    every point of the box — strict complementarity, arb-certified."""
    from flint import arb
    bp = base_oracle.part
    where = {}
    for bi, idx in enumerate(bp.blocks):
        for pos, t in enumerate(idx):
            where[int(t)] = ("z", bi, pos)
    for pos, t in enumerate(bp.border_idx):
        where[int(t)] = ("g", None, pos)
    with __import__("baller").hygiene.ctx_guard(prec=prec):
        Fz, Fg = masked.F_full(box)
        margins = {}
        for k, (side, bound) in active.items():
            reg, bi, pos = where[int(k)]
            val = Fz[bi][pos] if reg == "z" else Fg[pos]
            ok = (val > 0) if side == "lower" else (val < 0)
            margins[int(k)] = dict(side=side, bound=bound,
                                   grad_lower=float(val.lower()),
                                   grad_upper=float(val.upper()),
                                   ok=bool(ok))
            if not ok:
                return dict(ok=False, failed=f"kkt:coord_{k}",
                            margins=margins,
                            reason=(f"KKT multiplier sign at coord {k} "
                                    f"({side} bound {bound:g}) not "
                                    f"strictly certified: gradient "
                                    f"enclosure [{float(val.lower()):.3g},"
                                    f" {float(val.upper()):.3g}]"))
    return dict(ok=True, failed=None, margins=margins,
                reason="all active-bound multiplier signs strictly "
                       "certified over the box")


def _finish(res, oracle, theta, prec, t_all):
    import hashlib
    res["prec_bits"] = int(prec)
    res["prec_dps"] = round(prec * 0.30103, 1)
    res["dims"] = dict(blocks=oracle.dims[0], border=oracle.dims[1],
                       n_blocks=len(oracle.dims[0]),
                       n_total=int(sum(oracle.dims[0]) + oracle.dims[1]))
    res["theta_sha16"] = hashlib.sha256(
        np.ascontiguousarray(np.asarray(theta, float)).tobytes()
    ).hexdigest()[:16]
    res["wall_s_total"] = round(time.time() - t_all, 1)
    return res
