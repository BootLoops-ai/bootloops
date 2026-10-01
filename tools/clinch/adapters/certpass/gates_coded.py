"""gates_coded — oracle identity gates for the R4B coded-Poisson targets
(run BEFORE any certificate; CLINCH law; sunspot-adapter pattern).

  gate_value      arb oracle f (tight center) vs (a) the FROZEN float64
                  reference closure fit_multipliers.coded_nll_factory
                  (imported by identity from round4/r4b/code) vs (b) the
                  banked tune_nll vs (c) an independent mpmath 50-dps
                  transliteration using the reference's own DIRECT p3
                  form (1 - e^-lam (1+lam+lam^2/2)).
  gate_gradient_adjudicated
                  arb analytic gradient vs the reference closure's OWN
                  analytic float64 gradient (verbatim receipt) and vs
                  mpmath 50-dps central FD of the independent
                  transliteration at the top-K disagreement coords + the
                  intercept + the max-|g| coord. Bar 1e-10 in the
                  max(1,|g|) register.
  gate_fd_hessian arb assembled Hessian vs central FD of the arb gradient
                  at the exact power-of-two step 2^-20.
"""
import numpy as np
from flint import arb

from baller.hygiene import ctx_guard

PREC = 192
MP_DPS = 50
BAR_VALUE_F64 = 1e-9
BAR_VALUE_MP = 1e-12
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


# ---------------------------------------------------------------------
# independent mpmath transliteration (direct reference forms)
# ---------------------------------------------------------------------
class MpCoded:
    def __init__(self, data, dps=MP_DPS):
        from mpmath import mp, mpf
        self.dps = dps
        old = mp.dps
        mp.dps = dps
        try:
            self.X = np.minimum(np.asarray(data["y"]), 3).astype(int)
            self.o = [mpf(float(v)) for v in data["logref"]]
            Z = np.asarray(data["Zg"], float)
            self.Z = [[mpf(float(v)) for v in row] for row in Z]
            self.n, self.q = Z.shape
            self.l2 = mpf(float(data["l2"]))
        finally:
            mp.dps = old

    def f(self, x):
        from mpmath import mp, mpf, exp, log
        old = mp.dps
        mp.dps = self.dps
        try:
            b0 = mpf(float(x[0]))
            beta = [mpf(float(v)) for v in x[1:]]
            L2C = log(mpf(2))
            tot = mpf(0)
            for r in range(self.n):
                row = self.Z[r]
                eta = b0
                for a in range(self.q):
                    eta += row[a] * beta[a]
                if eta > 4:
                    eta = mpf(4)
                elif eta < -4:
                    eta = mpf(-4)
                loglam = self.o[r] + eta
                lam = exp(loglam)
                X = self.X[r]
                if X == 0:
                    tot += lam
                elif X == 1:
                    tot += lam - loglam
                elif X == 2:
                    tot += lam + L2C - 2 * loglam
                else:
                    p3 = 1 - exp(-lam) * (1 + lam + lam * lam / 2)
                    if p3 < mpf("1e-9"):
                        p3 = mpf("1e-9")
                    tot += -log(p3)
            rid = mpf(0)
            for v in beta:
                rid += v * v
            return tot / self.n + self.l2 * rid
        finally:
            mp.dps = old


def gate_value(oracle, mp_obj, fref, banked_nll, th):
    th = np.asarray(th, float)
    f64, _ = fref(th)
    va = arb_val(oracle, th)
    amid = float(va.mid())
    vmp = float(mp_obj.f(th))
    out = dict(
        f_ref_float64=float(f64), f_banked=float(banked_nll),
        dev_ref_vs_banked=abs(float(f64) - float(banked_nll)),
        f_arb_mid=amid, f_arb_rad=float(va.rad()),
        dev_arb_vs_ref_float64=abs(amid - float(f64)),
        f_mp50=vmp, dev_arb_vs_mp50=abs(amid - vmp),
        bars=dict(arb_vs_ref=BAR_VALUE_F64, arb_vs_mp50=BAR_VALUE_MP))
    out["PASS"] = (out["dev_arb_vs_ref_float64"] <= BAR_VALUE_F64
                   and out["dev_arb_vs_mp50"] <= BAR_VALUE_MP)
    return out


def gate_gradient_adjudicated(oracle, mp_obj, fref, th, K=5):
    from mpmath import mp, mpf
    th = np.asarray(th, float)
    p = len(th)
    ga = arb_grad(oracle, th)
    ga_mid = np.array([float(v.mid()) for v in ga])
    _, g64 = fref(th)          # the reference's OWN analytic gradient
    g64 = np.asarray(g64, float)
    dis = np.abs(ga_mid - g64) / np.maximum(1.0, np.abs(g64))
    idx = list(np.argsort(-dis)[:K]) + [0, int(np.argmax(np.abs(ga_mid)))]
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
        grad_arb_max_rad=float(max(v.rad() for v in ga)),
        note=("ref analytic float64 metric receipted verbatim; the claim "
              "rides the mpmath 50-dps central-FD adjudication register "
              "(truncation ~h^2 f''' ~ 1e-18 scale)"))
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
