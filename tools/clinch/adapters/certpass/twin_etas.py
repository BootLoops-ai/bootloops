"""twin_etas — vectorized float64 twin of the ETAS observed-data NLL.

Same formulas as oracle_etas (module doc there), numpy/scipy register.
Used for: (a) full-data value/gradient identity receipts vs the arb
oracle, (b) the EM-identity gates against the pinned-env state of record
(P_background / l_hat / G / n_hat), (c) the float64 Newton pre-polish.

Source-side time-factor TF uses upper_gamma_ext transliterated VERBATIM
from the pinned etas.inversion (sha-receipted by the dump); its
(lc, om, ltau) derivative columns are central-FD of the source sum
(receipted at the FD register — the rigorous claim rides the arb oracle
+ the mp50 slice adjudication, never this twin).
"""
import numpy as np
from scipy.special import gammaincc, gamma as gamma_func, exp1

LN10 = np.log(10.0)
IDX = dict(lmu=0, lk0=1, a=2, lc=3, om=4, ltau=5, ld=6, gam=7, rho=8)


def upper_gamma_ext(a, x):
    """VERBATIM transliteration of the pinned etas.inversion form."""
    if a > 0:
        return gammaincc(a, x) * gamma_func(a)
    elif a == 0:
        return exp1(x)
    else:
        return (upper_gamma_ext(a + 1, x) - np.power(x, a) * np.exp(-x)) / a


class TwinEtas:
    def __init__(self, data):
        self.pt = np.asarray(data["pair_t"], np.int64)
        self.pdt = np.asarray(data["pair_dt"], float)
        self.pr2 = np.asarray(data["pair_r2"], float)
        self.pdm = np.asarray(data["pair_m"], float) - float(data["mc_min"])
        self.ps = np.asarray(data["pair_s"], np.int64)
        self.src_dm = np.asarray(data["src_m"], float) - float(
            data["mc_min"])
        self.src_ts = np.asarray(data["src_ts"], float)
        self.src_te = np.asarray(data["src_te"], float)
        self.n_t = int(len(data["trg_pbg"]))
        self.n_s = int(len(data["src_m"]))
        self.A = float(data["area"])
        self.T = float(data["tw_length"])

    # ------------- pieces --------------------------------------------
    def _pair_g(self, th):
        lmu, lk0, a, lc, om, ltau, ld, gam, rho = th
        c = 10.0 ** lc
        tau = 10.0 ** ltau
        k0 = 10.0 ** lk0
        d = 10.0 ** ld
        P = self.pdt + c
        L = np.log(P)
        Tv = np.exp(-self.pdt / tau - (1 + om) * L)
        Db = d * np.exp(gam * self.pdm)
        R = self.pr2 + Db
        M = np.log(R)
        q = Db / R
        Sv = np.exp(-(1 + rho) * M)
        kap = k0 * np.exp(a * self.pdm)
        g = kap * Tv * Sv
        return g, dict(P=P, L=L, M=M, q=q, c=c, tau=tau)

    def lam(self, th):
        g, aux = self._pair_g(th)
        mu = 10.0 ** th[0]
        S0 = np.bincount(self.pt, weights=g, minlength=self.n_t)
        return mu + S0, g, aux

    def G_sum(self, th, sl=None):
        """Vector of G_j (source integrals)."""
        lmu, lk0, a, lc, om, ltau, ld, gam, rho = th
        c = 10.0 ** lc
        tau = 10.0 ** ltau
        k0 = 10.0 ** lk0
        d = 10.0 ** ld
        dm = self.src_dm if sl is None else self.src_dm[sl]
        ts = self.src_ts if sl is None else self.src_ts[sl]
        te = self.src_te if sl is None else self.src_te[sl]
        kap = k0 * np.exp(a * dm)
        Db = d * np.exp(gam * dm)
        area = np.pi * np.power(Db, -rho) / rho
        x1 = (ts + c) / tau
        x2 = (te + c) / tau
        tf = (np.exp(c / tau) * np.power(tau, -om)
              * (upper_gamma_ext(-om, x1) - upper_gamma_ext(-om, x2)))
        return kap * area * tf

    def value(self, th):
        lamv, _, _ = self.lam(th)
        mu = 10.0 ** th[0]
        return (-np.log(lamv).sum() + mu * self.A * self.T
                + self.G_sum(th).sum())

    # ------------- EM identity receipts ------------------------------
    def em_identities(self, th, dump):
        lamv, g, _ = self.lam(th)
        mu = 10.0 ** th[0]
        pbg = mu / lamv
        lhat = np.bincount(self.ps, weights=g / lamv[self.pt],
                           minlength=self.n_s)
        G = self.G_sum(th)
        return dict(
            dev_pbg=float(np.abs(pbg - dump["trg_pbg"]).max()),
            dev_lhat=float(np.abs(lhat - dump["src_lhat"]).max()),
            dev_G_rel=float((np.abs(G - dump["src_G"])
                             / np.maximum(np.abs(dump["src_G"]),
                                          1e-300)).max()),
            dev_nhat=float(abs(pbg.sum() - float(dump["n_hat"]))),
            n_hat=float(pbg.sum()))

    # ------------- gradient (hybrid) ---------------------------------
    def grad(self, th):
        """9-vector: pairs+mu analytic; source (lc,om,ltau) via central
        FD of the G-sum (h = 1e-7); everything else analytic."""
        th = np.asarray(th, float)
        lmu, lk0, a, lc, om, ltau, ld, gam, rho = th
        lamv, g, aux = self.lam(th)
        mu = 10.0 ** lmu
        il = 1.0 / lamv
        ilp = il[self.pt]
        W = LN10
        gil = g * ilp
        c, tau = aux["c"], aux["tau"]
        iP = 1.0 / aux["P"]
        L, M, q = aux["L"], aux["M"], aux["q"]
        gr = np.zeros(9)
        gr[0] = -(mu * W * il).sum() + mu * W * self.A * self.T
        gr[1] = -W * gil.sum()
        gr[2] = -(gil * self.pdm).sum()
        gr[3] = ((1 + om) * c * W) * (gil * iP).sum()
        gr[4] = (gil * L).sum()
        gr[5] = -(W / tau) * (gil * self.pdt).sum()
        gr[6] = (1 + rho) * W * (gil * q).sum()
        gr[7] = (1 + rho) * (gil * self.pdm * q).sum()
        gr[8] = (gil * M).sum()
        # source side, analytic block
        G = self.G_sum(th)
        dm = self.src_dm
        lnD = ld * W + gam * dm
        gr[1] += W * G.sum()
        gr[2] += (G * dm).sum()
        gr[6] += -rho * W * G.sum()
        gr[7] += -rho * (G * dm).sum()
        gr[8] += (G * (-lnD - 1.0 / rho)).sum()
        # source side, (lc, om, ltau) via central FD of sum G
        for k in (3, 4, 5):
            h = 1e-7 * max(1.0, abs(th[k]))
            tp = th.copy(); tp[k] += h
            tm = th.copy(); tm[k] -= h
            gr[k] += (self.G_sum(tp).sum() - self.G_sum(tm).sum()) \
                / (tp[k] - tm[k])
        return gr

    def pairs_val_grad_hess(self, th):
        """Analytic float64 (val, grad, hess) of the PAIRS + mu part of
        the NLL (everything except sum_j G_j), vectorized einsum over the
        same hand formulas as the arb oracle (oracle_etas module doc)."""
        th = np.asarray(th, float)
        lmu, lk0, a, lc, om, ltau, ld, gam, rho = th
        W = LN10
        mu = 10.0 ** lmu
        c = 10.0 ** lc
        tau = 10.0 ** ltau
        g, aux = self._pair_g(th)
        S0 = np.bincount(self.pt, weights=g, minlength=self.n_t)
        lam = mu + S0
        il = 1.0 / lam
        val = float(-np.log(lam).sum() + mu * self.A * self.T)
        iP = 1.0 / aux["P"]
        L, M, q = aux["L"], aux["M"], aux["q"]
        dt = self.pdt
        dm = self.pdm
        # d-vector components (indices 1..8), per pair
        dvec = [None] * 9
        dvec[1] = np.full(len(g), W)
        dvec[2] = dm
        dvec[3] = -(1 + om) * c * W * iP
        dvec[4] = -L
        dvec[5] = (W / tau) * dt
        dvec[6] = -(1 + rho) * W * q
        dvec[7] = -(1 + rho) * dm * q
        dvec[8] = -M
        q1q = q * (1 - q)
        hmap = {(3, 3): -(1 + om) * W * W * c * dt * iP * iP,
                (3, 4): -c * W * iP,
                (5, 5): -(W * W / tau) * dt,
                (6, 6): -(1 + rho) * W * W * q1q,
                (6, 7): -(1 + rho) * W * dm * q1q,
                (7, 7): -(1 + rho) * dm * dm * q1q,
                (6, 8): -W * q,
                (7, 8): -dm * q}
        nt = self.n_t
        S1 = np.zeros((9, nt))
        for k in range(1, 9):
            S1[k] = np.bincount(self.pt, weights=g * dvec[k],
                                minlength=nt)
        gl = np.zeros((9, nt))
        gl[0] = mu * W
        gl[1:] = S1[1:]
        grad = np.zeros(9)
        grad[0] = -(mu * W * il).sum() + mu * W * self.A * self.T
        for k in range(1, 9):
            grad[k] = -(S1[k] * il).sum()
        hess = np.zeros((9, 9))
        il2 = il * il
        for k in range(9):
            for l in range(k, 9):
                if k >= 1 and l >= 1:
                    e = dvec[k] * dvec[l]
                    hkl = hmap.get((k, l))
                    if hkl is not None:
                        e = e + hkl
                    S2kl = np.bincount(self.pt, weights=g * e,
                                       minlength=nt)
                else:
                    S2kl = np.zeros(nt)
                Hl = S2kl
                if k == 0 and l == 0:
                    Hl = Hl + mu * W * W
                hess[k, l] = (gl[k] * gl[l] * il2 - Hl * il).sum()
        hess[0, 0] += mu * W * W * self.A * self.T
        hess = np.triu(hess) + np.triu(hess, 1).T
        return val, grad, hess

    def hybrid_polish(self, th, oracle, maxit=40, tol_g=1e-6, log=None,
                      pinned=None):
        """Damped (Levenberg) Newton on the EXACT-register hybrid
        derivatives: analytic pairs (this twin) + arb-mid sources
        (oracle.sources_f64). Acceptance: gradient-norm decrease on the
        FREE coordinates. pinned: {coord: exact value} (corner
        candidates). Heuristic only — the certificate carries the rigor."""
        log = log or (lambda s: None)
        th = np.asarray(th, float).copy()
        pinned = pinned or {}
        for k, v in pinned.items():
            th[k] = v
        free = np.array([k for k in range(len(th)) if k not in pinned])

        def eval_at(t):
            sv, sg, sH = oracle.sources_f64(t)
            pv, pg, pH = self.pairs_val_grad_hess(t)
            return pg + sg, pH + sH

        g, H = eval_at(th)
        gn = float(np.abs(g[free]).max())
        mu_lm = 0.0
        hist = [dict(it=-1, gmax=gn)]
        step = np.zeros(len(free))
        for it in range(maxit):
            if gn < tol_g:
                break
            ok = False
            for _ in range(12):
                Hf = H[np.ix_(free, free)]
                Hd = Hf + mu_lm * np.diag(np.maximum(np.diag(Hf), 1e-12))
                try:
                    step = np.linalg.solve(Hd, g[free])
                except np.linalg.LinAlgError:
                    mu_lm = max(mu_lm * 10, 1e-6)
                    continue
                th_try = th.copy()
                th_try[free] -= step
                g2, H2 = eval_at(th_try)
                gn2 = float(np.abs(g2[free]).max())
                if gn2 < gn:
                    th, g, H, gn = th_try, g2, H2, gn2
                    mu_lm = mu_lm / 3 if mu_lm > 1e-14 else 0.0
                    ok = True
                    break
                mu_lm = max(mu_lm * 10, 1e-6)
            hist.append(dict(it=it, gmax=gn,
                             step_inf=float(np.abs(step).max()),
                             mu_lm=mu_lm))
            log(f"hybrid it{it}: gmax {gn:.3g} "
                f"step {np.abs(step).max():.3g} mu {mu_lm:.1g}")
            if not ok:
                log("hybrid: no gradient decrease — stopping")
                break
        return th, hist

    def hess_fd(self, th, h=2.0 ** -16):
        p = len(th)
        H = np.zeros((p, p))
        for a_ in range(p):
            tp = np.asarray(th, float).copy(); tp[a_] += h
            tm = np.asarray(th, float).copy(); tm[a_] -= h
            H[:, a_] = (self.grad(tp) - self.grad(tm)) / (2 * h)
        return 0.5 * (H + H.T)

    def newton_prepolish(self, th, maxit=60, tol=1e-9, log=None):
        """Backtracking line-search Newton on the twin NLL (full steps
        along flat valleys; value-decrease acceptance). Heuristic only —
        rigor rides the arb polish + certificate downstream."""
        log = log or (lambda s: None)
        th = np.asarray(th, float).copy()
        hist = []
        f0 = self.value(th)
        H = None
        for it in range(maxit):
            g = self.grad(th)
            if H is None or it % 2 == 0:
                H = self.hess_fd(th)
            try:
                step = np.linalg.solve(H, g)
            except np.linalg.LinAlgError:
                step = g / np.maximum(np.diag(H), 1e-9)
            si = float(np.abs(step).max())
            hist.append(dict(it=it, gmax=float(np.abs(g).max()),
                             step_inf=si, f=float(f0)))
            log(f"prepolish it{it}: gmax {np.abs(g).max():.3g} "
                f"step {si:.3g} f {f0:.6f}")
            if si < tol:
                break
            accepted = False
            for s in (1.0, 0.5, 0.25, 0.1, 0.03, 0.01):
                th_try = th - s * step
                f_try = self.value(th_try)
                if f_try < f0:
                    th, f0 = th_try, f_try
                    accepted = True
                    break
            if not accepted:
                log("prepolish: no descent step accepted — stopping")
                break
            if accepted and s != 1.0:
                H = None          # force H refresh after damped move
        return th, hist
