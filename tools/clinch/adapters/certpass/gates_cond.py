"""gates_cond — oracle identity gates for the R4A conditional-spatial
targets (run BEFORE any certificate; CLINCH law).

  gate_value      arb value (tight center) vs (a) the frozen-text float64
                  closure (assemble_r4a's nll_frozen receipt, re-derived
                  here from the SAME banked constants) vs (b) mpmath
                  50-dps independent transliteration.
  gate_gradient_adjudicated
                  arb analytic gradient vs the frozen closure's analytic
                  float64 gradient (verbatim) and vs mp50 central FD of
                  the independent transliteration (top-K + max-|g|).
  gate_fd_hessian arb assembled H vs central FD of the arb gradient at
                  the exact 2^-20 step.
"""
import numpy as np
from scipy.special import logsumexp, softmax

from baller.hygiene import ctx_guard

PREC = 192
MP_DPS = 50
BAR_VALUE_F64 = 1e-6
BAR_VALUE_MP = 1e-9
BAR_GRAD_MP = 1e-10
BAR_HESS_FD = 1e-6
FD_HMP = 1e-9
FD_HESS_STEP = 2.0 ** -20


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


def frozen_closure(data):
    """The float64 reference closure (identical arithmetic to the frozen
    fit_cond inner nll, on the banked standardized design)."""
    Xs = np.asarray(data["Xs"], float)
    o = np.log(np.maximum(np.asarray(data["lam_f"], float), 1e-9)) \
        if "lam_f" in data else np.asarray(data["o"], float)
    y = np.asarray(data["y"], float)
    ridge = float(data["ridge"])
    yy = y.astype(np.float64)
    Yd = yy.sum(1)
    wd = 1.0 / np.maximum(Yd, 1.0)
    gconst = -((wd[:, None] * yy)[:, :, None] * Xs).sum((0, 1))
    wY = wd * Yd

    def nll(b):
        z = o + Xs @ b
        lse = logsumexp(z, axis=1)
        f = (-((wd[:, None] * yy) * z).sum() + (wY * lse).sum()
             + ridge * (b ** 2).sum())
        w = softmax(z, axis=1)
        g = (gconst + ((wY[:, None] * w)[:, :, None] * Xs).sum((0, 1))
             + 2 * ridge * b)
        return f, g

    return nll


class MpCond:
    """Independent mpmath 50-dps transliteration (active days only —
    inactive days contribute exactly zero)."""

    def __init__(self, data, dps=MP_DPS):
        from mpmath import mp, mpf
        self.dps = dps
        Xs = np.asarray(data["Xs"], float)
        o = np.asarray(data["o"], float)
        y = np.asarray(data["y"], float)
        Yd = y.sum(1)
        act = np.where(Yd > 0)[0]
        old = mp.dps
        mp.dps = dps
        try:
            self.X = [[[mpf(float(v)) for v in Xs[d, c]]
                       for c in range(Xs.shape[1])] for d in act]
            self.o = [[mpf(float(v)) for v in o[d]] for d in act]
            self.y = [[mpf(float(v)) for v in y[d]] for d in act]
            self.wd = [mpf(1.0) / mpf(float(max(Yd[d], 1.0)))
                       for d in act]
            self.Yd = [mpf(float(Yd[d])) for d in act]
            self.ridge = mpf(float(data["ridge"]))
            self.k = Xs.shape[1]
            self.p = Xs.shape[2]
        finally:
            mp.dps = old

    def f(self, b):
        from mpmath import mp, mpf, exp, log
        old = mp.dps
        mp.dps = self.dps
        try:
            bb = [mpf(float(v)) for v in b]
            tot = mpf(0)
            for i in range(len(self.X)):
                zs = []
                for c in range(self.k):
                    row = self.X[i][c]
                    z = self.o[i][c]
                    for a in range(self.p):
                        z += row[a] * bb[a]
                    zs.append(z)
                M = max(zs)
                S = mpf(0)
                for z in zs:
                    S += exp(z - M)
                lse = M + log(S)
                acc = mpf(0)
                for c in range(self.k):
                    if self.y[i][c] != 0:
                        acc -= self.y[i][c] * zs[c]
                tot += self.wd[i] * (acc + self.Yd[i] * lse)
            rid = mpf(0)
            for v in bb:
                rid += v * v
            return tot + self.ridge * rid
        finally:
            mp.dps = old


def gate_value(oracle, mp_obj, fref, th, f_banked_receipt=None):
    th = np.asarray(th, float)
    f64, _ = fref(th)
    va = arb_val(oracle, th)
    amid = float(va.mid())
    vmp = float(mp_obj.f(th))
    out = dict(
        f_ref_float64=float(f64),
        f_arb_mid=amid, f_arb_rad=float(va.rad()),
        dev_arb_vs_ref_float64=abs(amid - float(f64)),
        f_mp50=vmp, dev_arb_vs_mp50=abs(amid - vmp),
        bars=dict(arb_vs_ref=BAR_VALUE_F64, arb_vs_mp50=BAR_VALUE_MP))
    if f_banked_receipt is not None:
        out["f_assembly_receipt"] = float(f_banked_receipt)
        out["dev_ref_vs_assembly_receipt"] = abs(
            float(f64) - float(f_banked_receipt))
    out["PASS"] = (out["dev_arb_vs_ref_float64"] <= BAR_VALUE_F64
                   and out["dev_arb_vs_mp50"] <= BAR_VALUE_MP)
    return out


def gate_gradient_adjudicated(oracle, mp_obj, fref, th, K=4):
    from mpmath import mp, mpf
    th = np.asarray(th, float)
    p = len(th)
    ga = arb_grad(oracle, th)
    ga_mid = np.array([float(v.mid()) for v in ga])
    _, g64 = fref(th)
    g64 = np.asarray(g64, float)
    dis = np.abs(ga_mid - g64) / np.maximum(1.0, np.abs(g64))
    idx = list(np.argsort(-dis)[:K]) + [int(np.argmax(np.abs(ga_mid)))]
    idx = sorted(set(int(i) for i in idx))
    adj = {}
    worst = 0.0
    for a in idx:
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
        m = abs(ga_mid[a] - gm) / max(1.0, abs(gm))
        adj[int(a)] = dict(g_arb=float(ga_mid[a]), g_mp50=gm, metric=m,
                           g_ref64=float(g64[a]))
        worst = max(worst, m)
    out = dict(
        ref_analytic_metric=float(np.max(dis)),
        adjudicated_coords=adj, adjudicated_worst=worst,
        bar_mp50=BAR_GRAD_MP,
        grad_arb_max_rad=float(max(v.rad() for v in ga)))
    out["PASS"] = worst <= BAR_GRAD_MP
    return out, ga_mid


def gate_fd_hessian(oracle, th, prec=PREC, coords=None):
    th = np.asarray(th, float)
    p = len(th)
    part = oracle.part
    with ctx_guard(prec=prec):
        D, B, G = oracle.H(tight_box(part, th, prec))
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
    for a in (range(p) if coords is None else coords):
        tp = th.copy(); tp[a] += h
        tm = th.copy(); tm[a] -= h
        gp = np.array([float(v.mid()) for v in arb_grad(oracle, tp, prec)])
        gm = np.array([float(v.mid()) for v in arb_grad(oracle, tm, prec)])
        col = (gp - gm) / (2 * h)
        m = np.abs(col - Hm[:, a]) / np.maximum(1.0, np.abs(Hm[:, a]))
        j = int(np.argmax(m))
        if m[j] > worst:
            worst = float(m[j]); worst_at = (j, int(a))
    return dict(step=h, metric_max=worst, worst_entry=worst_at,
                bar=BAR_HESS_FD, PASS=worst <= BAR_HESS_FD,
                note="FD of the arb gradient: truncation-only error")
