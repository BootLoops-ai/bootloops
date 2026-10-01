"""gates_etas — oracle identity gates for the ETAS reference-fit targets.

Register design (stated in the addendum): the mp50 independent-
transliteration adjudication runs on a deterministic REDUCED SLICE (the
same objective formula over receipted reduced index sets, identical for
all three implementations — the mp gate certifies the FORMULA and its
derivatives, incl. the omega-derivative of the incomplete-gamma factor
via mpmath's own generalized gammainc, an implementation independent of
the oracle's Kummer-jet route). The FULL-data identities ride:
  (a) arb-vs-float64-twin value/gradient receipts (full data);
  (b) the EM identities against the pinned-env state of record:
      P_background / l_hat / G / n_hat reproduction (twin, full data);
  (c) FD-of-arb-gradient vs assembled arb Hessian on the slice oracle
      (all 9 columns, exact 2^-16 step).
"""
import numpy as np

from baller.hygiene import ctx_guard

PREC = 192
MP_DPS = 50
BAR_VALUE_TWIN = 1e-7          # relative (float64 sum register)
BAR_VALUE_MP = 1e-9            # slice, absolute/relative hybrid
BAR_GRAD_MP = 1e-10
BAR_GRAD_TWIN_ANALYTIC = 1e-7  # max(1,|g|) register
BAR_GRAD_TWIN_FD = 1e-4        # (lc,om,ltau) source-FD columns
BAR_HESS_FD = 1e-6
BAR_EM_PBG = 1e-12
BAR_EM_LHAT = 1e-9
BAR_EM_G_REL = 1e-9
FD_HMP = 1e-9
FD_HESS_STEP = 2.0 ** -16


def _center(part, th):
    return ([[float(th[t]) for t in blk] for blk in part.blocks],
            [float(th[t]) for t in part.border_idx])


def tight_box(part, th, prec=PREC):
    import baller.certify as bc
    bk = bc.block_krawczyk
    with ctx_guard(prec=prec):
        return bk.box_around(_center(part, th), 0.0)


def arb_val(oracle, th, prec=PREC):
    with ctx_guard(prec=prec):
        return oracle.nll(tight_box(oracle.part, th, prec))


def arb_grad(oracle, th, prec=PREC):
    with ctx_guard(prec=prec):
        g, _ = oracle._grad_full(tight_box(oracle.part, th, prec))
    return g


def slice_data(data, n_t_red=120, n_s_red=300):
    """Deterministic reduced index sets: first n_t_red targets (with all
    their pairs), first n_s_red//2 sources + n_s_red//2 from the middle
    (covers both the ts>0 aux branch and the ts=0 branch)."""
    pt = np.asarray(data["pair_t"], np.int64)
    keep = pt < n_t_red
    n_s = len(data["src_m"])
    half = n_s_red // 2
    sidx = np.r_[0:half, n_s // 2:n_s // 2 + half]
    out = dict(
        pair_t=pt[keep], pair_s=np.asarray(data["pair_s"])[keep],
        pair_dt=np.asarray(data["pair_dt"])[keep],
        pair_r2=np.asarray(data["pair_r2"])[keep],
        pair_m=np.asarray(data["pair_m"])[keep],
        src_m=np.asarray(data["src_m"])[sidx],
        src_ts=np.asarray(data["src_ts"])[sidx],
        src_te=np.asarray(data["src_te"])[sidx],
        trg_pbg=np.asarray(data["trg_pbg"])[:n_t_red],
        src_lhat=np.zeros(len(sidx)), src_G=np.zeros(len(sidx)),
        n_hat=np.float64(0),
        mc_min=data["mc_min"], area=data["area"],
        tw_length=data["tw_length"])
    rec = dict(n_t_red=int(n_t_red), n_pairs_red=int(keep.sum()),
               src_idx=[int(v) for v in sidx[:3]] + ["..."],
               n_s_red=int(len(sidx)),
               n_aux=(int((out["src_ts"] > 0).sum())))
    return out, rec


class MpEtas:
    """Independent mpmath 50-dps transliteration of the SLICE objective.
    Time factor via mpmath.gammainc(s, a=x1, b=x2) — mpmath's own
    generalized incomplete gamma, independent of the Kummer-jet route."""

    def __init__(self, sdata, dps=MP_DPS):
        from mpmath import mp, mpf
        self.dps = dps
        old = mp.dps
        mp.dps = dps
        try:
            self.pt = np.asarray(sdata["pair_t"], np.int64)
            self.pdt = [mpf(float(v)) for v in sdata["pair_dt"]]
            self.pr2 = [mpf(float(v)) for v in sdata["pair_r2"]]
            self.pdm = [mpf(float(v) - float(sdata["mc_min"]))
                        for v in sdata["pair_m"]]
            self.sdm = [mpf(float(v) - float(sdata["mc_min"]))
                        for v in sdata["src_m"]]
            self.sts = [mpf(float(v)) for v in sdata["src_ts"]]
            self.ste = [mpf(float(v)) for v in sdata["src_te"]]
            self.n_t = int(len(sdata["trg_pbg"]))
            self.A = mpf(float(sdata["area"]))
            self.T = mpf(float(sdata["tw_length"]))
        finally:
            mp.dps = old

    def f(self, th):
        from mpmath import mp, mpf, exp, log, pi, power, gammainc
        old = mp.dps
        mp.dps = self.dps
        try:
            W = log(mpf(10))
            lmu, lk0, a, lc, om, ltau, ld, gam, rho = [
                mpf(float(v)) for v in th]
            mu = exp(lmu * W)
            c = exp(lc * W)
            tau = exp(ltau * W)
            k0 = exp(lk0 * W)
            d = exp(ld * W)
            lamt = [mu for _ in range(self.n_t)]
            for i in range(len(self.pt)):
                dm = self.pdm[i]
                kap = k0 * exp(a * dm)
                P = self.pdt[i] + c
                Tv = exp(-self.pdt[i] / tau) * power(P, -(1 + om))
                Db = d * exp(gam * dm)
                Sv = power(self.pr2[i] + Db, -(1 + rho))
                lamt[self.pt[i]] += kap * Tv * Sv
            tot = mpf(0)
            for t in range(self.n_t):
                tot -= log(lamt[t])
            tot += mu * self.A * self.T
            s = -om
            for j in range(len(self.sdm)):
                dm = self.sdm[j]
                kap = k0 * exp(a * dm)
                Db = d * exp(gam * dm)
                area = pi * power(Db, -rho) / rho
                x1 = (self.sts[j] + c) / tau
                x2 = (self.ste[j] + c) / tau
                tf = exp(c / tau) * power(tau, -om) \
                    * gammainc(s, a=x1, b=x2)
                tot += kap * area * tf
            return tot
        finally:
            mp.dps = old


def gate_value(oracle_full, twin, oracle_slice, mp_obj, twin_slice, th):
    th = np.asarray(th, float)
    va = arb_val(oracle_full, th)
    amid = float(va.mid())
    vtwin = float(twin.value(th))
    vs = arb_val(oracle_slice, th)
    vsm = float(vs.mid())
    vmp = float(mp_obj.f(th))
    vts = float(twin_slice.value(th))
    scale = max(1.0, abs(vtwin))
    out = dict(
        f_arb_mid=amid, f_arb_rad=float(va.rad()), f_twin=vtwin,
        dev_arb_vs_twin_rel=abs(amid - vtwin) / scale,
        f_slice_arb=vsm, f_slice_arb_rad=float(vs.rad()),
        f_slice_mp50=vmp, f_slice_twin=vts,
        dev_slice_arb_vs_mp50_rel=abs(vsm - vmp) / max(1.0, abs(vmp)),
        dev_slice_twin_vs_mp50_rel=abs(vts - vmp) / max(1.0, abs(vmp)),
        bars=dict(twin=BAR_VALUE_TWIN, mp=BAR_VALUE_MP))
    out["PASS"] = (out["dev_arb_vs_twin_rel"] <= BAR_VALUE_TWIN
                   and out["dev_slice_arb_vs_mp50_rel"] <= BAR_VALUE_MP)
    return out


def gate_gradient_adjudicated(oracle_slice, mp_obj, twin, th,
                              oracle_full=None):
    from mpmath import mp, mpf
    th = np.asarray(th, float)
    # full-data arb vs twin (receipt registers per column class)
    rec = {}
    if oracle_full is not None:
        ga = arb_grad(oracle_full, th)
        ga_mid = np.array([float(v.mid()) for v in ga])
        gt = twin.grad(th)
        m = np.abs(ga_mid - gt) / np.maximum(1.0, np.abs(gt))
        rec["full_arb_vs_twin_analytic_cols"] = float(
            np.max(m[[0, 1, 2, 6, 7, 8]]))
        rec["full_arb_vs_twin_fd_cols"] = float(np.max(m[[3, 4, 5]]))
        rec["full_grad_arb"] = [float(v) for v in ga_mid]
        rec["full_grad_arb_max_rad"] = float(max(v.rad() for v in ga))
        pass_full = (rec["full_arb_vs_twin_analytic_cols"]
                     <= BAR_GRAD_TWIN_ANALYTIC
                     and rec["full_arb_vs_twin_fd_cols"]
                     <= BAR_GRAD_TWIN_FD)
    else:
        pass_full = True
    # slice mp50 central-FD adjudication, ALL 9 coords
    gs = arb_grad(oracle_slice, th)
    gs_mid = np.array([float(v.mid()) for v in gs])
    adj = {}
    worst = 0.0
    for a in range(9):
        h = FD_HMP * max(1.0, abs(th[a]))
        tp = th.copy(); tp[a] += h
        tm = th.copy(); tm[a] -= h
        old = mp.dps
        mp.dps = MP_DPS
        try:
            step = mpf(float(tp[a])) - mpf(float(tm[a]))
            gm = (mp_obj.f(tp) - mp_obj.f(tm)) / step
        finally:
            mp.dps = old
        gm = float(gm)
        m = abs(gs_mid[a] - gm) / max(1.0, abs(gm))
        adj[a] = dict(g_arb_slice=float(gs_mid[a]), g_mp50=gm, metric=m)
        worst = max(worst, m)
    out = dict(rec, adjudicated_coords=adj, adjudicated_worst=worst,
               bar_mp50=BAR_GRAD_MP,
               note=("mp50 slice adjudication covers all 9 coords incl. "
                     "the omega/incomplete-gamma route; full-data twin "
                     "receipts verbatim (FD columns at the FD register)"))
    out["PASS"] = bool(pass_full and worst <= BAR_GRAD_MP)
    return out, gs_mid


def gate_fd_hessian(oracle_slice, th, prec=PREC):
    """FD of the slice-oracle arb gradient vs its assembled H (all 9)."""
    th = np.asarray(th, float)
    p = 9
    part = oracle_slice.part
    with ctx_guard(prec=prec):
        D, B, G = oracle_slice.H(tight_box(part, th, prec))
    Hm = np.zeros((p, p))
    bidx = [int(t) for t in part.border_idx]
    for bi, blk in enumerate(part.blocks):
        blk = [int(t) for t in blk]
        for i, s in enumerate(blk):
            for j, t in enumerate(blk):
                Hm[s, t] = float(D[bi][i, j].mid())
            for j, t in enumerate(bidx):
                Hm[s, t] = float(B[bi][i, j].mid())
                Hm[t, s] = Hm[s, t]
    for i, s in enumerate(bidx):
        for j, t in enumerate(bidx):
            Hm[s, t] = float(G[i, j].mid())
    h = FD_HESS_STEP
    worst = 0.0
    worst_at = None
    for a in range(p):
        tp = th.copy(); tp[a] += h
        tm = th.copy(); tm[a] -= h
        gp = np.array([float(v.mid())
                       for v in arb_grad(oracle_slice, tp, prec)])
        gm = np.array([float(v.mid())
                       for v in arb_grad(oracle_slice, tm, prec)])
        col = (gp - gm) / (2 * h)
        m = np.abs(col - Hm[:, a]) / np.maximum(1.0, np.abs(Hm[:, a]))
        j = int(np.argmax(m))
        if m[j] > worst:
            worst = float(m[j]); worst_at = (j, int(a))
    return dict(step=h, metric_max=worst, worst_entry=worst_at,
                bar=BAR_HESS_FD, PASS=worst <= BAR_HESS_FD,
                note="slice oracle, all 9 columns; same code path as the "
                     "full-data H")


def gate_em_identities(twin, th, dump):
    ids = twin.em_identities(th, dump)
    g_reg = dump.get("g_register")
    if g_reg is not None and b"penultimate" in bytes(g_reg):
        # np-assembled table: the banked CSV G column is the
        # PENULTIMATE-EM-iteration register (receipted); the final-theta
        # G identity instead rides the l_hat/P_bg/n_hat bars (final
        # register) + the package-dump diff receipt
        bar_G = 1e-4
        ids["G_register"] = "banked-penultimate (receipted bar 1e-4)"
    else:
        bar_G = BAR_EM_G_REL
    ids["bars"] = dict(pbg=BAR_EM_PBG, lhat=BAR_EM_LHAT, G_rel=bar_G)
    ids["PASS"] = (ids["dev_pbg"] <= BAR_EM_PBG
                   and ids["dev_lhat"] <= BAR_EM_LHAT
                   and ids["dev_G_rel"] <= bar_G)
    return ids
