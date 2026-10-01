"""clinch.oracle_v36 — the CENTERED-space oracle (v3.6 unpinned / v3.7
pinned) via the AUXILIARY-BORDER KKT EXTENSION — the L3 design of record
(MANUAL.md; battery leg text).

THE PROBLEM. v3.6/v3.7 spaces carry HARD MEAN-ZERO CENTERING on the surv
and adv species effects: eta uses u_raw[sp] - ubar with ubar =
mean_species(u_raw) — every latent of a channel couples to every other
through the mean, so the latent Hessian is family-block + rank-2 and the
plain block-arrow form is lost.

THE EXTENSION. Introduce ubar as an auxiliary BORDER variable per centered
channel with its constraint multiplier:
    variables  x = (theta, m_s, m_a, mu_s, mu_a[, lam])
    L = Phi(theta, m) + mu_s (m_s - ubar_s(theta))
                      + mu_a (m_a - ubar_a(theta))
                      [+ lam (T(border, m_a) - T_ref)]        (v3.7 pin)
where Phi is the penalized NLL with eta using u_raw[sp] - m_channel (m an
INDEPENDENT variable — the latent-latent coupling through the mean is
gone: the theta-theta block of grad L's Jacobian is family-block-arrow
again), and T is the latents-zero advance total of fit18_common (the pin;
in extended variables it touches ONLY border coordinates + m_a). Zeros of
F = grad L correspond exactly to stationary points of the true composed
objective (m = ubar, mu = -dPhi/dm resolve them; with the pin, to
KKT-stationary points of the pinned problem, lam_ext = -lambda in
fit18_common's sign convention). The extended system is a block-arrow with
border = 53 model globals + 4 (or 5) auxiliaries; baller's structured
Krawczyk applies UNCHANGED.

SEMANTICS NOTE (recorded honestly): the extended Jacobian is a KKT matrix
— indefinite by construction (mu/lam rows) — so the KRAWCZYK leg
(existence + uniqueness of the KKT point in the box) certifies the
stationarity half; the PD statement for centered spaces is the REDUCED
Hessian's and rides the rank-4 congruence extension,
clinch.kkt_pd.pd_reduced_congruence (the signed-congruence inertia
certificate; battery leg L8 proves it on the toy centered fixture). This
oracle declares its multiplier rows via kkt_multiplier_idx, which routes
the engine's PD leg there. The L3 battery statistic is the
contraction-margin comparison, which is Krawczyk-side.

Survival climate surface (v3.6): clim_s = g1*MAT + g2*fs + g4*seas +
g3b*MAT*fs with fs = f_sat(MAP_mm; m0), m0 = exp(log_m0b_surv), per-site
fs/dfs/d2fs enclosed from the log_m0b ball (dfs/dlm0 = -fs(1-fs),
d2fs/dlm02 = (2 fs - 1) dfs). ADVANCE keeps the fit13 linear plane.
Fecundity untouched (uncentered). Priors per the reference v3.6 register
(g_surv sd [0.5, 4.9, 0.5]; g3b sd 4.9; log_m0b ~ N(7.0901, 1.0986^2)).
All other math identical to oracle_v31 (log-space forms, endpoint hulls,
clip guards; the reference is likelihood.py's v36 branch).
"""
import numpy as np
from flint import arb, arb_mat

from .jets import Jet, dig_ball, trig_ball, lgam_ball
from .oracle_v31 import (OracleRefusal, Partition, _Acc, V31_GLOBAL_SLOTS)

__all__ = ["check_space_v36", "ModelV36Oracle", "V36_GLOBAL_SLOTS",
           "N_GLOBAL36"]

V36_GLOBAL_SLOTS = V31_GLOBAL_SLOTS + [("g3b_surv", 1), ("log_m0b_surv", 1)]
N_GLOBAL36 = sum(w for _, w in V36_GLOBAL_SLOTS)          # 53
M0B_LOG_CENTER = 7.090076835776092         # v3lib.M0B_PRIOR_LOG_CENTER
M0B_LOG_SD = 1.0986122886681098            # v3lib.M0B_PRIOR_LOG_SD
SD_G_SAT = 4.9                             # likelihood.PRIOR_SD['g_sat']
CLIMATE_MAP_CENTER = 2000.0                # v3lib.CLIMATE_CENTER[1]
CLIMATE_MAP_SCALE = 1000.0                 # v3lib.CLIMATE_SCALE[1]


def check_space_v36(slice_bounds, n_total):
    """Scope guard: refuse any dial space that is not EXACTLY v3.6."""
    names = [nm for nm, _ in V36_GLOBAL_SLOTS]
    got = {nm: b for nm, b in slice_bounds.items()
           if not nm.startswith("z_")}
    if sorted(got) != sorted(names):
        raise OracleRefusal(
            "check_space", f"not the v3.6 space: extra="
            f"{sorted(set(got) - set(names))} missing="
            f"{sorted(set(names) - set(got))}")
    off = 0
    for nm, w in V36_GLOBAL_SLOTS:
        if tuple(slice_bounds[nm]) != (off, off + w):
            raise OracleRefusal("check_space",
                                f"slot {nm} at {slice_bounds[nm]}, "
                                f"expected ({off},{off + w})")
        off += w
    return True


class ModelV36Oracle:
    """Extended-KKT block-arrow oracle for the v3.6 space (pinned=v3.7).

    arrs: the SAME 10-unit cell arrays as fit13 (registry equality is the
    ADAPTER's assertion); consts: v3.6 slice_bounds + the v3.1 prior
    constants. Extended layout: [53 globals][m_s, m_a, mu_s, mu_a(, lam)]
    [latents shifted by naux]."""

    def __init__(self, arrs, consts, pinned=False, T_ref=None):
        check_space_v36(consts["slice_bounds"], consts["n_total"])
        self.arrs, self.consts = arrs, consts
        self.pinned = bool(pinned)
        if self.pinned:
            assert T_ref is not None, "pinned oracle needs T_ref"
        self.T_ref = None if T_ref is None else float(T_ref)
        self.naux = 5 if self.pinned else 4
        sb = consts["slice_bounds"]
        self.sl = {nm: sb[nm][0] for nm in sb}     # REAL indices (globals)
        self.zoff = self.naux                      # latent ext shift
        self.ix_ms = N_GLOBAL36
        self.ix_ma = N_GLOBAL36 + 1
        self.ix_mus = N_GLOBAL36 + 2
        self.ix_mua = N_GLOBAL36 + 3
        self.ix_lam = N_GLOBAL36 + 4 if self.pinned else None
        self.n_ext = consts["n_total"] + self.naux
        self.part = Partition(
            sb, self.n_ext, consts["S"], consts["NG"], consts["NF"],
            arrs["gen_of"], arrs["fam_of"], consts["families"],
            n_border=N_GLOBAL36 + self.naux, z_shift=self.naux)
        S = consts["S"]
        self.cnt_gen = np.bincount(np.asarray(arrs["gen_of"]),
                                   minlength=consts["NG"]).astype(float)
        self.cnt_fam = np.bincount(np.asarray(arrs["fam_of"]),
                                   minlength=consts["NF"]).astype(float)
        self._const_cache = {}

    @property
    def dims(self):
        return self.part.dims

    @property
    def kkt_multiplier_idx(self):
        """Border offsets of the constraint-multiplier rows (mu_s, mu_a
        and, under the pin, lam) — the declaration clinch.kkt_pd's
        reduced-Hessian inertia certificate consumes; the engine routes
        its PD leg there when this is present."""
        out = [self.ix_mus, self.ix_mua]
        if self.pinned:
            out.append(self.ix_lam)
        return tuple(int(a) for a in out)

    def F(self, x):
        th = self.part.flatten(x)
        _, grad, _ = self.evaluate(th, order2=False)
        return self.part.split_grad(grad)

    def H(self, x):
        th = self.part.flatten(x)
        _, _, acc = self.evaluate(th, order2=True)
        return acc.to_mats()

    # ------------------------------------------------------------------
    def ubar_np(self, theta_real):
        """(ubar_surv, ubar_adv) at a float REAL-layout theta (center
        initialization; heuristic floats — rigor lives in the residual)."""
        sb = self.consts["slice_bounds"]
        v = {nm: np.asarray(theta_real[s[0]:s[1]]) for nm, s in sb.items()}
        sig = np.exp(v["log_sig"])
        S = self.consts["S"]
        out = []
        for p, tag in ((0, "surv"), (1, "adv")):
            u = (sig[3 * p] * v[f"z_sp_{tag}"].sum()
                 + sig[3 * p + 1] * (self.cnt_gen @ v[f"z_gen_{tag}"])
                 + sig[3 * p + 2] * (self.cnt_fam @ v[f"z_fam_{tag}"])) / S
            out.append(float(u))
        return out[0], out[1]

    def ext_center(self, theta_real, mu_s=0.0, mu_a=0.0, lam=0.0):
        """Extended float center [globals, aux, latents]."""
        th = np.asarray(theta_real, dtype=float)
        ms, ma = self.ubar_np(th)
        aux = [ms, ma, float(mu_s), float(mu_a)]
        if self.pinned:
            aux.append(float(lam))
        return np.concatenate([th[:N_GLOBAL36], np.array(aux),
                               th[N_GLOBAL36:]])

    def init_multipliers(self, theta_real, lam=0.0, prec=192):
        """Run one gradient pass at mu=0 (lam fixed) and return the ext
        center with mu = -mid(F_m rows) — the values resolving the m-rows."""
        from baller.hygiene import ctx_guard
        th0 = self.ext_center(theta_real, 0.0, 0.0, lam)
        with ctx_guard(prec=prec):
            thb = [arb(float(t)) for t in th0]
            _, grad, _ = self.evaluate(thb, order2=False)
            mu_s = -float(grad[self.ix_ms].mid())
            mu_a = -float(grad[self.ix_ma].mid())
        return self.ext_center(theta_real, mu_s, mu_a, lam)

    # ------------------------------------------------------------------
    def _nll_const(self, prec_key):
        if prec_key in self._const_cache:
            return self._const_cache[prec_key]
        a = self.arrs
        s = arb(0)
        for n, y in zip(a["sv_n"], a["sv_die"]):
            s -= (arb(float(n) + 1).lgamma() - arb(float(y) + 1).lgamma()
                  - arb(float(n) - float(y) + 1).lgamma())
        for R in a["rc_count"]:
            s += arb(float(R) + 1).lgamma()
        self._const_cache[prec_key] = s
        return s

    def evaluate(self, th, order2=True):
        from flint import ctx
        a, sl = self.arrs, self.sl
        part = self.part
        zoff = self.zoff
        n_ext = self.n_ext
        grad = [arb(0)] * n_ext
        cen = {}
        acc = _Acc(part) if order2 else None
        nll = self._nll_const(ctx.prec)

        def gv(name, i=0):
            return th[sl[name] + i]

        sig = [gv("log_sig", j).exp() for j in range(9)]
        c19 = arb("1e-9"); c110 = arb("1e-10"); c112 = arb("1e-12")
        ln19 = c19.log(); ln110 = c110.log()

        def logsigmoid(x):
            if float(x.v.mid()) >= 0.0:
                return -((-x).exp().log1p())
            return x - x.exp().log1p()

        # per-site saturating water surface from the log_m0b ball
        m0 = gv("log_m0b_surv").exp()
        C = a["C"]
        n_site = C.shape[0]
        fs_site, dfs_site, d2fs_site = [], [], []
        for p in range(n_site):
            mm = arb(float(C[p, 1]) * CLIMATE_MAP_SCALE + CLIMATE_MAP_CENTER)
            fs = mm / (mm + m0)
            dfs = -fs * (1 - fs)
            fs_site.append(fs)
            dfs_site.append(dfs)
            d2fs_site.append((2 * fs - 1) * dfs)
        ig1 = sl["g_surv"] + 1
        ig3b = sl["g3b_surv"]
        ilm0 = sl["log_m0b_surv"]

        # ---------------- survival (jets in eta, log_nu) ----------------
        ilognu = sl["log_nu"]
        lognu_j = Jet.var(gv("log_nu"), 1, 2, order2)
        nu_j = lognu_j.exp()
        lg_nu_j = nu_j._chain(lgam_ball(nu_j.v), dig_ball(nu_j.v),
                              None if not order2 else trig_ball(nu_j.v))
        dig_cache = {}

        def lg_of_nu_plus(nfl):
            got = dig_cache.get(nfl)
            if got is None:
                w = nu_j.v + nfl
                got = (lgam_ball(w), dig_ball(w),
                       None if not order2 else trig_ball(w))
                dig_cache[nfl] = got
            lv, dg, tg = got
            return nu_j._chain(lv, dg, tg)

        zs_b = sl["z_sp_surv"] + zoff
        zg_b = sl["z_gen_surv"] + zoff
        zf_b = sl["z_fam_surv"] + zoff
        for c in range(len(a["sv_n"])):
            sp = int(a["sv_sp"][c]); k = int(a["sv_k"][c])
            st = int(a["sv_site"][c])
            n = float(a["sv_n"][c]); y = float(a["sv_die"][c])
            dt5 = float(a["sv_dt5"][c])
            x0 = float(a["X"][sp, 0]); x1 = float(a["X"][sp, 1])
            kc = (k - 3.0) / 2.0
            zb = float(a["sv_zbar"][c]); xd = float(a["sv_xdry"][c])
            Dz = float(a["D_site"][st]) * zb
            gen = int(a["gen_of"][sp]); fam = int(a["fam_of"][sp])
            izs = zs_b + sp; izg = zg_b + gen; izf = zf_b + fam
            zsv, zgv, zfv = th[izs], th[izg], th[izf]
            C0 = float(C[st, 0]); C2 = float(C[st, 2])
            fs = fs_site[st]; dfs = dfs_site[st]; d2fs = d2fs_site[st]
            g1 = gv("g_surv", 1); g3b = gv("g3b_surv")
            eta = (gv("a_surv", k) + x0 * gv("b_surv", 0)
                   + x1 * gv("b_surv", 1) + (x0 * kc) * gv("b_sxk")
                   + C0 * gv("g_surv", 0) + fs * g1
                   + C2 * gv("g_surv", 2) + (C0 * fs) * g3b
                   + zb * gv("w_surv", 0) + Dz * gv("w_surv", 1)
                   + xd * gv("w_surv", 2)
                   + sig[0] * zsv + sig[1] * zgv + sig[2] * zfv
                   - th[self.ix_ms])
            ej = Jet.var(eta, 0, 2, order2)
            lp5 = logsigmoid(ej); lm5 = logsigmoid(-ej)
            # clip CENSUS, not refusal — unclipped-objective target (see
            # oracle_v31 docstring/census note)
            if not (lp5.v > ln19 and lm5.v > ln19):
                cen["sv_p5"] = cen.get("sv_p5", 0) + 1
            lpsur = lp5 * dt5
            psur = lpsur.exp()
            pdie = -(lpsur.expm1())
            if not (pdie.v > c110 and lpsur.v > ln110):
                cen["sv_pdie"] = cen.get("sv_pdie", 0) + 1
            al = pdie * nu_j
            be = psur * nu_j
            ll = ((al + y).lgamma() + (be + (n - y)).lgamma()
                  - lg_of_nu_plus(n) + lg_nu_j - al.lgamma() - be.lgamma())
            nll -= ll.v
            ge = -ll.g[0]; gn = -ll.g[1]
            coords = ((sl["a_surv"] + k, 1.0), (sl["b_surv"], x0),
                      (sl["b_surv"] + 1, x1), (sl["b_sxk"], x0 * kc),
                      (sl["g_surv"], C0), (sl["g_surv"] + 2, C2),
                      (sl["w_surv"], zb), (sl["w_surv"] + 1, Dz),
                      (sl["w_surv"] + 2, xd), (self.ix_ms, -1.0))
            bcoords = ((izs, sig[0]), (izg, sig[1]), (izf, sig[2]),
                       (sl["log_sig"], sig[0] * zsv),
                       (sl["log_sig"] + 1, sig[1] * zgv),
                       (sl["log_sig"] + 2, sig[2] * zfv),
                       (ig1, fs), (ig3b, C0 * fs),
                       (ilm0, (g1 + C0 * g3b) * dfs))
            for idx, cf in coords:
                grad[idx] += ge * cf
            for idx, cf in bcoords:
                grad[idx] += ge * cf
            grad[ilognu] += gn
            if order2:
                hee = -ll.h[0]; hen = -ll.h[1]; hnn = -ll.h[2]
                allc = coords + bcoords
                for i in range(len(allc)):
                    ia, ca = allc[i]
                    hca = hee * ca
                    for j in range(i, len(allc)):
                        ib, cb = allc[j]
                        acc.add(ia, ib, hca * cb)
                    acc.add(ia, ilognu, hen * ca)
                acc.add(ilognu, ilognu, hnn)
                acc.add(izs, sl["log_sig"], ge * sig[0])
                acc.add(izg, sl["log_sig"] + 1, ge * sig[1])
                acc.add(izf, sl["log_sig"] + 2, ge * sig[2])
                acc.add(sl["log_sig"], sl["log_sig"], ge * (sig[0] * zsv))
                acc.add(sl["log_sig"] + 1, sl["log_sig"] + 1,
                        ge * (sig[1] * zgv))
                acc.add(sl["log_sig"] + 2, sl["log_sig"] + 2,
                        ge * (sig[2] * zfv))
                # d2eta from the saturating surface
                acc.add(ig1, ilm0, ge * dfs)
                acc.add(ig3b, ilm0, ge * (C0 * dfs))
                acc.add(ilm0, ilm0, ge * ((g1 + C0 * g3b) * d2fs))

        # ---------------- advance (jet in eta) ---------------------------
        zs_a = sl["z_sp_adv"] + zoff
        zg_a = sl["z_gen_adv"] + zoff
        zf_a = sl["z_fam_adv"] + zoff
        for c in range(len(a["av_n"])):
            sp = int(a["av_sp"][c]); k = int(a["av_k"][c])
            st = int(a["av_site"][c])
            n = float(a["av_n"][c]); y = float(a["av_up"][c])
            dt5 = float(a["av_dt5"][c])
            x0 = float(a["X"][sp, 0]); x1 = float(a["X"][sp, 1])
            zb = float(a["av_zbar"][c]); xd = float(a["av_xdry"][c])
            Dz = float(a["D_site"][st]) * zb
            gen = int(a["gen_of"][sp]); fam = int(a["fam_of"][sp])
            izs = zs_a + sp; izg = zg_a + gen; izf = zf_a + fam
            zsv, zgv, zfv = th[izs], th[izg], th[izf]
            C0 = float(C[st, 0]); C1 = float(C[st, 1]); C2 = float(C[st, 2])
            eta = (gv("a_adv", k) + x0 * gv("b_gro", 0)
                   + x1 * gv("b_gro", 1)
                   + C0 * gv("g_gro", 0) + C1 * gv("g_gro", 1)
                   + C2 * gv("g_gro", 2)
                   + zb * gv("w_gro", 0) + Dz * gv("w_gro", 1)
                   + xd * gv("w_gro", 2)
                   + sig[3] * zsv + sig[4] * zgv + sig[5] * zfv
                   - th[self.ix_ma])
            ej = Jet.var(eta, 0, 1, order2)
            lq5 = logsigmoid(ej); lm5a = logsigmoid(-ej)
            if not (lq5.v > ln19 and lm5a.v > ln19):
                cen["av_q5"] = cen.get("av_q5", 0) + 1
            lstay = lm5a * dt5
            q = -(lstay.expm1())
            if not (lstay.v > ln110 and q.v > c110):
                cen["av_q_stay5"] = cen.get("av_q_stay5", 0) + 1
            ll = q.log() * y + lstay * (n - y)
            nll -= ll.v
            ge = -ll.g[0]
            coords = ((sl["a_adv"] + k, 1.0), (sl["b_gro"], x0),
                      (sl["b_gro"] + 1, x1),
                      (sl["g_gro"], C0), (sl["g_gro"] + 1, C1),
                      (sl["g_gro"] + 2, C2), (sl["w_gro"], zb),
                      (sl["w_gro"] + 1, Dz), (sl["w_gro"] + 2, xd),
                      (self.ix_ma, -1.0))
            bcoords = ((izs, sig[3]), (izg, sig[4]), (izf, sig[5]),
                       (sl["log_sig"] + 3, sig[3] * zsv),
                       (sl["log_sig"] + 4, sig[4] * zgv),
                       (sl["log_sig"] + 5, sig[5] * zfv))
            for idx, cf in coords:
                grad[idx] += ge * cf
            for idx, cf in bcoords:
                grad[idx] += ge * cf
            if order2:
                hee = -ll.h[0]
                allc = coords + bcoords
                for i in range(len(allc)):
                    ia, ca = allc[i]
                    hca = hee * ca
                    for j in range(i, len(allc)):
                        ib, cb = allc[j]
                        acc.add(ia, ib, hca * cb)
                acc.add(izs, sl["log_sig"] + 3, ge * sig[3])
                acc.add(izg, sl["log_sig"] + 4, ge * sig[4])
                acc.add(izf, sl["log_sig"] + 5, ge * sig[5])
                acc.add(sl["log_sig"] + 3, sl["log_sig"] + 3,
                        ge * (sig[3] * zsv))
                acc.add(sl["log_sig"] + 4, sl["log_sig"] + 4,
                        ge * (sig[4] * zgv))
                acc.add(sl["log_sig"] + 5, sl["log_sig"] + 5,
                        ge * (sig[5] * zfv))

        # ------------- fecundity (uncentered; v3.1 form verbatim) --------
        offs = self.consts["fec_offsets"]
        ilogphi = sl["log_phi"]; ia1 = sl["a_fec"] + 1
        logphi_j = Jet.var(gv("log_phi"), 2, 3, order2)
        phi_j = logphi_j.exp()
        lgphi_j = phi_j._chain(lgam_ball(phi_j.v), dig_ball(phi_j.v),
                               None if not order2 else trig_ball(phi_j.v))
        zs_f = sl["z_sp_fec"] + zoff
        zg_f = sl["z_gen_fec"] + zoff
        zf_f = sl["z_fam_fec"] + zoff
        for c in range(len(a["rc_count"])):
            sp = int(a["rc_sp"][c]); st = int(a["rc_site"][c])
            R = float(a["rc_count"][c]); dt5 = float(a["rc_dt5"][c])
            zb = float(a["rc_zbar"][c])
            x0 = float(a["X"][sp, 0]); x1 = float(a["X"][sp, 1])
            grp = int(a["group_of"][sp])
            gen = int(a["gen_of"][sp]); fam = int(a["fam_of"][sp])
            izs = zs_f + sp; izg = zg_f + gen; izf = zf_f + fam
            zsv, zgv, zfv = th[izs], th[izg], th[izf]
            C0 = float(C[st, 0]); C1 = float(C[st, 1]); C2 = float(C[st, 2])
            eta = (gv("a_fec", 0) + x0 * gv("b_fec", 0)
                   + x1 * gv("b_fec", 1)
                   + C0 * gv("g_fec", 0) + C1 * gv("g_fec", 1)
                   + C2 * gv("g_fec", 2)
                   + zb * gv("lam", grp)
                   + sig[6] * zsv + sig[7] * zgv + sig[8] * zfv)
            ej = Jet.var(eta, 0, 3, order2)
            a1j = Jet.var(gv("a_fec", 1), 1, 3, order2)
            O_SHIFT = 1.5
            Bj = Jet.const(arb(0), 3, order2)
            for j, ad in enumerate(a["rc_adults"][c]):
                adf = float(ad)
                if adf != 0.0:
                    Bj = Bj + (a1j * (float(offs[j]) - O_SHIFT)).exp() * adf
            lmu = (ej + a1j * O_SHIFT + Bj.log()) + arb(dt5).log()
            if not (lmu.v > c112.log()):
                cen["rc_mu"] = cen.get("rc_mu", 0) + 1
            dlm = lmu - logphi_j
            if float(dlm.v.mid()) <= 0.0:
                lpm = logphi_j + (dlm.exp()).log1p()
            else:
                lpm = lmu + ((-dlm).exp()).log1p()
            ll = ((phi_j + R).lgamma() - lgphi_j
                  + phi_j * (logphi_j - lpm) + (lmu - lpm) * R)
            nll -= ll.v
            ge = -ll.g[0]; ga1 = -ll.g[1]; gph = -ll.g[2]
            coords = ((sl["a_fec"], 1.0), (sl["b_fec"], x0),
                      (sl["b_fec"] + 1, x1),
                      (sl["g_fec"], C0), (sl["g_fec"] + 1, C1),
                      (sl["g_fec"] + 2, C2), (sl["lam"] + grp, zb))
            bcoords = ((izs, sig[6]), (izg, sig[7]), (izf, sig[8]),
                       (sl["log_sig"] + 6, sig[6] * zsv),
                       (sl["log_sig"] + 7, sig[7] * zgv),
                       (sl["log_sig"] + 8, sig[8] * zfv))
            from .jets import hidx
            for idx, cf in coords:
                grad[idx] += ge * cf
            for idx, cf in bcoords:
                grad[idx] += ge * cf
            grad[ia1] += ga1
            grad[ilogphi] += gph
            if order2:
                hee = -ll.h[hidx(3, 0, 0)]
                hea = -ll.h[hidx(3, 0, 1)]
                hep = -ll.h[hidx(3, 0, 2)]
                haa = -ll.h[hidx(3, 1, 1)]
                hap = -ll.h[hidx(3, 1, 2)]
                hpp = -ll.h[hidx(3, 2, 2)]
                allc = coords + bcoords
                for i in range(len(allc)):
                    ia, ca = allc[i]
                    hca = hee * ca
                    for j in range(i, len(allc)):
                        ib, cb = allc[j]
                        acc.add(ia, ib, hca * cb)
                    acc.add(ia, ia1, hea * ca)
                    acc.add(ia, ilogphi, hep * ca)
                acc.add(ia1, ia1, haa)
                acc.add(ia1, ilogphi, hap)
                acc.add(ilogphi, ilogphi, hpp)
                acc.add(izs, sl["log_sig"] + 6, ge * sig[6])
                acc.add(izg, sl["log_sig"] + 7, ge * sig[7])
                acc.add(izf, sl["log_sig"] + 8, ge * sig[8])
                acc.add(sl["log_sig"] + 6, sl["log_sig"] + 6,
                        ge * (sig[6] * zsv))
                acc.add(sl["log_sig"] + 7, sl["log_sig"] + 7,
                        ge * (sig[7] * zgv))
                acc.add(sl["log_sig"] + 8, sl["log_sig"] + 8,
                        ge * (sig[8] * zfv))

        # ---------------- priors (v3.6 register) -------------------------
        co = self.consts
        psd = co["prior_sd"]
        nll_box = [nll]

        def gauss(name, width, center, sd):
            base = sl[name]
            for i in range(width):
                cen = center[i] if isinstance(center, (list, tuple,
                                                       np.ndarray)) else center
                sdi = sd[i] if isinstance(sd, (list, tuple,
                                               np.ndarray)) else sd
                r = (th[base + i] - float(cen)) / float(sdi)
                nll_box[0] += 0.5 * (r * r)
                grad[base + i] += r / float(sdi)
                if order2:
                    acc.add(base + i, base + i,
                            arb(1) / (float(sdi) * float(sdi)))

        gauss("a_surv", 7, co["pc_surv_logit"], psd["a_surv"])
        gauss("a_adv", 6, 0.0, psd["a_adv"])
        gauss("a_fec", 2, [co["pc_fec_log"], 0.0],
              [psd["a_fec0"], psd["a_fec1"]])
        for nm in ("b_surv", "b_gro", "b_fec"):
            gauss(nm, 2, 0.0, psd["b"])
        gauss("b_sxk", 1, 0.0, psd["b_sxk"])
        gauss("g_surv", 3, 0.0, [psd["g"], SD_G_SAT, psd["g"]])
        gauss("g_gro", 3, 0.0, psd["g"])
        gauss("g_fec", 3, 0.0, psd["g"])
        gauss("g3b_surv", 1, 0.0, SD_G_SAT)
        gauss("log_m0b_surv", 1, M0B_LOG_CENTER, M0B_LOG_SD)
        gauss("log_sig", 9, co["pc_log_sig"], psd["log_sig"])
        gauss("lam", 3, 0.0, psd["lam"])
        gauss("w_surv", 3, 0.0, psd["w"])
        gauss("w_gro", 3, 0.0, psd["w"])
        gauss("log_phi", 1, co["pc_log_phi"], psd["log_phi"])
        gauss("log_nu", 1, co["pc_log_nu"], psd["log_nu"])
        nll = nll_box[0]
        for p in ("surv", "adv", "fec"):
            for lev in ("sp", "gen", "fam"):
                base, hi = co["slice_bounds"][f"z_{lev}_{p}"]
                for t in range(base, hi):
                    z = th[t + zoff]
                    nll += 0.5 * (z * z)
                    grad[t + zoff] += z
                    if order2:
                        acc.add(t + zoff, t + zoff, arb(1))

        # -------- centering constraints + multiplier couplings -----------
        S = float(co["S"])
        for chan, (p0, ix_m, ix_mu) in (("surv", (0, self.ix_ms,
                                                  self.ix_mus)),
                                        ("adv", (1, self.ix_ma,
                                                 self.ix_mua))):
            p, ix_m_, ix_mu_ = p0, ix_m, ix_mu
            mu = th[ix_mu_]
            zsl = co["slice_bounds"][f"z_sp_{chan}"]
            gsl = co["slice_bounds"][f"z_gen_{chan}"]
            fsl = co["slice_bounds"][f"z_fam_{chan}"]
            s0, s1, s2 = sig[3 * p], sig[3 * p + 1], sig[3 * p + 2]
            ls0, ls1, ls2 = (sl["log_sig"] + 3 * p, sl["log_sig"] + 3 * p
                             + 1, sl["log_sig"] + 3 * p + 2)
            # zbar sums (balls)
            sum_sp = arb(0)
            for t in range(zsl[0], zsl[1]):
                sum_sp += th[t + zoff]
            sum_gen = arb(0)
            for j, t in enumerate(range(gsl[0], gsl[1])):
                sum_gen += th[t + zoff] * float(self.cnt_gen[j])
            sum_fam = arb(0)
            for j, t in enumerate(range(fsl[0], fsl[1])):
                sum_fam += th[t + zoff] * float(self.cnt_fam[j])
            zbar0 = sum_sp / S
            zbar1 = sum_gen / S
            zbar2 = sum_fam / S
            ubar = s0 * zbar0 + s1 * zbar1 + s2 * zbar2
            # F_m += mu ; F_mu = m - ubar ; grad_z/ls -= mu * d ubar
            grad[ix_m_] += mu
            grad[ix_mu_] += th[ix_m_] - ubar
            for j, t in enumerate(range(zsl[0], zsl[1])):
                grad[t + zoff] += -mu * (s0 / S)
                if order2:
                    acc.add(t + zoff, ix_mu_, -(s0 / S))
                    acc.add(t + zoff, ls0, -mu * (s0 / S))
            for j, t in enumerate(range(gsl[0], gsl[1])):
                cj = float(self.cnt_gen[j])
                if cj:
                    grad[t + zoff] += -mu * (s1 * (cj / S))
                    if order2:
                        acc.add(t + zoff, ix_mu_, -(s1 * (cj / S)))
                        acc.add(t + zoff, ls1, -mu * (s1 * (cj / S)))
            for j, t in enumerate(range(fsl[0], fsl[1])):
                cj = float(self.cnt_fam[j])
                if cj:
                    grad[t + zoff] += -mu * (s2 * (cj / S))
                    if order2:
                        acc.add(t + zoff, ix_mu_, -(s2 * (cj / S)))
                        acc.add(t + zoff, ls2, -mu * (s2 * (cj / S)))
            grad[ls0] += -mu * (s0 * zbar0)
            grad[ls1] += -mu * (s1 * zbar1)
            grad[ls2] += -mu * (s2 * zbar2)
            if order2:
                acc.add(ix_m_, ix_mu_, arb(1))
                acc.add(ls0, ix_mu_, -(s0 * zbar0))
                acc.add(ls1, ix_mu_, -(s1 * zbar1))
                acc.add(ls2, ix_mu_, -(s2 * zbar2))
                acc.add(ls0, ls0, -mu * (s0 * zbar0))
                acc.add(ls1, ls1, -mu * (s1 * zbar1))
                acc.add(ls2, ls2, -mu * (s2 * zbar2))

        # ---------------- the v3.7 pin (KKT border row) ------------------
        if self.pinned:
            lam_v = th[self.ix_lam]
            Tsum = arb(0)
            for c in range(len(a["av_n"])):
                k = int(a["av_k"][c]); st = int(a["av_site"][c])
                sp = int(a["av_sp"][c])
                n = float(a["av_n"][c]); dt5 = float(a["av_dt5"][c])
                x0 = float(a["X"][sp, 0]); x1 = float(a["X"][sp, 1])
                zb = float(a["av_zbar"][c]); xd = float(a["av_xdry"][c])
                Dz = float(a["D_site"][st]) * zb
                C0 = float(C[st, 0]); C1 = float(C[st, 1])
                C2 = float(C[st, 2])
                eta0 = (gv("a_adv", k) + x0 * gv("b_gro", 0)
                        + x1 * gv("b_gro", 1)
                        + C0 * gv("g_gro", 0) + C1 * gv("g_gro", 1)
                        + C2 * gv("g_gro", 2)
                        + zb * gv("w_gro", 0) + Dz * gv("w_gro", 1)
                        + xd * gv("w_gro", 2) - th[self.ix_ma])
                ej = Jet.var(eta0, 0, 1, order2)
                lq5 = logsigmoid(ej); lm5p = logsigmoid(-ej)
                if not (lq5.v > ln19 and lm5p.v > ln19):
                    cen["pin_q5"] = cen.get("pin_q5", 0) + 1
                padv = -( (lm5p * dt5).expm1() )     # 1 - stay5
                Tc = padv * n
                Tsum += Tc.v
                w = Tc.g[0]                          # dT_c/deta0
                coords = ((sl["a_adv"] + k, 1.0), (sl["b_gro"], x0),
                          (sl["b_gro"] + 1, x1),
                          (sl["g_gro"], C0), (sl["g_gro"] + 1, C1),
                          (sl["g_gro"] + 2, C2), (sl["w_gro"], zb),
                          (sl["w_gro"] + 1, Dz), (sl["w_gro"] + 2, xd),
                          (self.ix_ma, -1.0))
                for idx, cf in coords:
                    grad[idx] += (lam_v * w) * cf
                    if order2:
                        acc.add(idx, self.ix_lam, w * cf)
                if order2:
                    h0 = Tc.h[0]
                    lh = lam_v * h0
                    for i in range(len(coords)):
                        ia, ca = coords[i]
                        for j in range(i, len(coords)):
                            ib, cb = coords[j]
                            acc.add(ia, ib, lh * (ca * cb))
            grad[self.ix_lam] += Tsum - self.T_ref
        self.last_clip_census = dict(cen)
        return nll, grad, acc
