#!/usr/bin/env python3
r"""frobenius_d3.py — G-block RE-LAND DESIGN #3 side copy of frobenius.land()
(frobenius.py itself is NOT touched — kept frozen by design).

WHY (measured evidence):
  design #1 (5-pt ray):            cond 6.05e21-5.29e22 (150x log spread)
  design #2 (col-equil, 2ring x 5ph, 250x25, match_dps 90): cond 9.08e12
    (5 orders short of the 1e8 bar)
  branch census (gbranches_r2048_j4): post-shear clusters = lambda=-3 with
    log tower l<=2  AND  lambda=-3+2eps_n, 2eps_n = 2(a+ib)/(c*rden),
    |2eps_n| = 2/2048.  The cond thief is the angle between u^{2eps_n} and
    the log tower span{1, L, L^2}: best quadratic-in-L approximation of
    e^{2eps_n L} over the match-point L-cloud leaves (2eps_n * dL)^3/3!
    => cond ~ 6/(2eps_n * dL)^3 * tower prefactor,  dL = log-u spread.
    Design #2's dL ~= 1  =>  ~1e13.  NO diagonal rescaling can fix this
    (van der Sluis: max-abs col equil is within sqrt(n_col) of ANY column
    scaling), and land() already solves by QR with cond_est = single-power
    SVD of the equilibrated matrix (no normal-equations squaring exists).

DESIGN #3 (this module): grow dL with
  (a) more rings across a wider radius span (Re L), and
  (b) MONODROMY SHEETS (Im L): match points on analytic continuations of
      the data around u=0 (theta -> theta + 2*pi*s).  Each loop adds 2*pi
      of Im-L spread at O(ring-transport) cost; the log tower is polynomial
      in s while the eps-branch is exponential e^{4*pi*i*eps_n*s} — the
      lp1family monodromy-separation pattern ("tail u^{j+m eps} branch
      separation by MONODROMY").
  plus (c) ROW equilibration = the GLS weighting for per-datum RELATIVE
      accuracy (march data is relative-accurate; with rings spanning
      decades the u^-3 row profile otherwise dominates the LS).  Exact,
      disclosed; kappa unscaled after solve.  BOTH conds are reported:
      cond_est_cols (design-#2's exact observable, column-equil only) and
      cond_est (the doubly-equilibrated solve matrix).

The Frobenius tower is IMPORTED UNCHANGED from production frobenius.py
(frobenius_basis — byte-identical numerics); this module only re-implements
the matching leg.  Sheet-aware basis eval: Y(u, theta_cum) built from the
basis dict with L = ln|u| + i*theta_cum (theta_cum = CUMULATIVE angle, NOT
reduced mod 2pi; integer u-powers in shear factors D_i(u) are single-valued
so only exp(M L) sees the sheet).

Transport route (full mode): ray -> outer ring, then per sheet s_min..s_max
ring sweeps with angular waypoints every <=pi/4 (chord min distance
>= cos(pi/8)*ring) and radial hops at fixed angle.  cond_only mode does NO
transports (pilot observable = exact production matrix cond).
"""

import sys
sys.dont_write_bytecode = True
from mpmath import mp, mpf, mpc

try:
    from .frobenius import frobenius_basis
    from .transport import _to_mpc, transport_fixed_eps
except ImportError:
    from frobenius import frobenius_basis
    from transport import _to_mpc, transport_fixed_eps

__all__ = ["land3", "eval_sheet"]


def eval_sheet(basis, u, theta_cum, wp):
    """Y(u) on the sheet where arg-history = theta_cum (L = ln|u| + i*theta).
    Built from the basis dict (M, P, shear_factors at basis['dps'])."""
    n = basis["n"]
    kmax = basis["kmax"]
    with mp.workdps(wp):
        uv = _to_mpc(u)
        if uv == 0:
            raise ValueError("eval_sheet at u=0")
        L = mp.log(abs(uv)) + mpc(0, 1) * mpf(theta_cum)
        E = mp.expm(mp.matrix(basis["M"]) * L)
        P = basis["P"]
        acc = [[mpc(v) for v in row] for row in P[kmax]]
        for q in range(kmax - 1, -1, -1):
            Pq = P[q]
            for i in range(n):
                ai = acc[i]
                Pqi = Pq[i]
                for j in range(n):
                    ai[j] = ai[j] * uv + Pqi[j]
        Y = [[mp.fsum(acc[i][k] * E[k, j] for k in range(n))
              for j in range(n)] for i in range(n)]
        for T, mtop in reversed(basis["shear_factors"]):
            for i in range(mtop):
                for j in range(n):
                    Y[i][j] = Y[i][j] * uv
            Y = [[mp.fsum(T[i][k] * Y[k][j] for k in range(n))
                  for j in range(n)] for i in range(n)]
        return Y


def _design_points(base_arg, r_f, rings, phases, sheets):
    """[(R_frac, theta_cum)] in TRANSPORT order: sheets ascending, rings
    outer->inner->outer alternating per sheet, phases swept monotonically.
    theta_cum is the cumulative angle (base_arg + phase + 2*pi*sheet)."""
    pts = []
    o_sorted = sorted(phases)
    for si, s in enumerate(sorted(sheets)):
        ring_seq = list(rings) if si % 2 == 0 else list(reversed(rings))
        for ri, R in enumerate(ring_seq):
            seq = o_sorted if ri % 2 == 0 else list(reversed(o_sorted))
            for o in seq:
                th = base_arg + mpf(str(o)) + 2 * mp.pi * s
                pts.append((mpf(str(R)), th))
    return pts


def land3(desys, eps0, x_from, y_from, x_sing, dps, kmax,
          match_rings=(0.5, 0.125), match_phase_offsets=None,
          match_sheets=(0,), match_dps=None, cluster_tol=None,
          equilibrate_rows=True, cond_only=False, basis=None,
          max_leg_angle=None):
    """Design-#3 Frobenius landing.  Returns the frobenius.land()-style dict
    plus: cond_est (row+col equilibrated solve matrix), cond_est_cols
    (column-equil only — design-#2's observable), match_design='rings3'.
    Solver is mp.qr_solve ONLY (no normal-equations fallback at all)."""
    wp = dps + 30
    if basis is None:
        basis = frobenius_basis(desys, eps0, x_sing, dps, kmax,
                                cluster_tol=cluster_tol)
    with mp.workdps(wp):
        eps_v = _to_mpc(eps0)
        xs = _to_mpc(x_sing)
        xf = _to_mpc(x_from)
        L_from = abs(xf - xs)
        if L_from == 0:
            raise ValueError("x_from coincides with x_sing")
        r_f = mpf(basis["r_circle"])
        offs = (match_phase_offsets if match_phase_offsets is not None
                else [-mp.pi/3, -mp.pi/6, mpf(0), mp.pi/6, mp.pi/3])
        for R in match_rings:
            if not (0 < mpf(str(R)) < 1):
                raise ValueError(f"match_rings fraction {R} not in (0,1)")
        base_arg = mp.arg((xf - xs) / L_from)
        pts = _design_points(base_arg, r_f, match_rings, offs, match_sheets)
        n = desys.n
        # ---- transports (skipped in cond_only pilot mode) ----
        ys = []
        if not cond_only:
            leg_cap = mp.pi / 4 if max_leg_angle is None else mpf(str(max_leg_angle))
            t_dps = match_dps if match_dps is not None else wp
            cur_x = xf
            cur_y = [_to_mpc(v) for v in y_from]
            cur_R, cur_th = None, None
            for (R, th) in pts:
                if cur_R is None:
                    # ray descent to the first point's angle/radius: first go
                    # radially (same arg as incoming ray), then swing.
                    tgt0 = xs + R * r_f * mp.exp(mpc(0, 1) * base_arg)
                    cur_y = transport_fixed_eps(desys, eps_v, cur_x, tgt0,
                                                cur_y, t_dps)
                    cur_x, cur_R, cur_th = tgt0, R, base_arg
                # angular waypoints at the CURRENT radius, then radial hop
                dth = th - cur_th
                nleg = max(1, int(mp.ceil(abs(dth) / leg_cap)))
                for k in range(1, nleg + 1):
                    thk = cur_th + dth * mpf(k) / nleg
                    tgt = xs + cur_R * r_f * mp.exp(mpc(0, 1) * thk)
                    cur_y = transport_fixed_eps(desys, eps_v, cur_x, tgt,
                                                cur_y, t_dps)
                    cur_x = tgt
                cur_th = th
                if R != cur_R:
                    tgt = xs + R * r_f * mp.exp(mpc(0, 1) * th)
                    cur_y = transport_fixed_eps(desys, eps_v, cur_x, tgt,
                                                cur_y, t_dps)
                    cur_x, cur_R = tgt, R
                with mp.workdps(wp):
                    cur_y = [mpc(v) for v in cur_y]
                ys.append(list(cur_y))
        # ---- stack the match system ----
        rows = []
        rhs = []
        for pi_, (R, th) in enumerate(pts):
            u = R * r_f * mp.exp(mpc(0, 1) * th)   # value of u (sheet in th)
            Yp = eval_sheet(basis, u, th, wp)
            for i in range(n):
                rows.append([mpc(v) for v in Yp[i]])
                rhs.append(mpc(ys[pi_][i]) if not cond_only else mpc(0))
        m = len(rows)
        # row equilibration (GLS relative-noise weighting; rhs scaled too)
        row_scale = None
        if equilibrate_rows:
            row_scale = []
            for i in range(m):
                mx = max(abs(v) for v in rows[i])
                mx = mx if mx > 0 else mpf(1)
                row_scale.append(mx)
                rows[i] = [v / mx for v in rows[i]]
                rhs[i] = rhs[i] / mx
        # column equilibration (design #2, unchanged)
        col_scale = []
        for j in range(n):
            mx = max(abs(rows[i][j]) for i in range(m))
            col_scale.append(mx if mx > 0 else mpf(1))
        for i in range(m):
            rows[i] = [rows[i][j] / col_scale[j] for j in range(n)]
        Amat = mp.matrix(rows)
        bvec = mp.matrix(rhs)

        def _cond(Mx):
            try:
                sv = mp.svd_c(Mx, compute_uv=False)
                smax, smin = max(sv), min(sv)
                return mp.nstr(smax / smin, 8) if smin > 0 else "inf"
            except (ValueError, ZeroDivisionError, ArithmeticError):
                return None
        cond_est = _cond(Amat)
        cond_cols = None
        if equilibrate_rows:
            # design-#2 observable: column-equil ONLY (rebuild w/o row scale)
            rows2 = []
            for i in range(m):
                rows2.append([rows[i][j] * col_scale[j] * row_scale[i]
                              for j in range(n)])
            cs2 = []
            for j in range(n):
                mx = max(abs(rows2[i][j]) for i in range(m))
                cs2.append(mx if mx > 0 else mpf(1))
            for i in range(m):
                rows2[i] = [rows2[i][j] / cs2[j] for j in range(n)]
            cond_cols = _cond(mp.matrix(rows2))
        else:
            cond_cols = cond_est
        u_tail = max(mpf(str(R)) for R in match_rings) * r_f
        tail = mpf(basis["p_last_norm"]) * u_tail ** basis["kmax"]
        if cond_only:
            with mp.workdps(dps):
                return {"cond_est": cond_est, "cond_est_cols": cond_cols,
                        "match_design": "rings3",
                        "equilibrated": True,
                        "equilibrated_rows": bool(equilibrate_rows),
                        "n_points": len(pts),
                        "rings": [str(R) for R in match_rings],
                        "sheets": [int(s) for s in match_sheets],
                        "phases": [mp.nstr(mpf(str(o)), 8) for o in offs],
                        "series_tail": +tail, "basis": basis,
                        "dps": dps, "kmax": kmax, "cond_only": True}
        kappa_m, resnorm = mp.qr_solve(Amat, bvec)   # QR ONLY — loud on fail
        kappa = [mpc(kappa_m[j]) / col_scale[j] for j in range(n)]
        bnorm = mp.norm(bvec)
        res_rel = mpf(resnorm) / max(mpf(1), bnorm)
    with mp.workdps(dps):
        return {"kappa": [+v for v in kappa],
                "residual": +mpf(resnorm), "residual_rel": +res_rel,
                "solver": "qr",
                "cond_est": cond_est, "cond_est_cols": cond_cols,
                "match_design": "rings3", "equilibrated": True,
                "equilibrated_rows": bool(equilibrate_rows),
                "n_points": len(pts),
                "sheets": [int(s) for s in match_sheets],
                "series_tail": +tail, "basis": basis,
                "dps": dps, "kmax": kmax}
