"""clinch.adapter_v31 — the CLINCH model adapter for the v3.1
dial layout (fit13/fit17/fit18 model of record).

Data flows through the reference pipeline's OWN code by identity import (refbridge ->
the fit runners' load_data — every receipted data-identity gate runs
untouched; REGISTRY keys asserted against the receipts). The adapter then
extracts the flat cell arrays the oracle consumes and caches them as an
npz (sha-pinned provenance) so certificate iteration does not re-run the
half-hour assembly; the FLAGSHIP gate path records whether it ran fresh
or from cache.

ORACLE IDENTITY GATES (issued BEFORE any certificate):
  gate_gradient  arb oracle (tight balls, high prec) vs the receipted
                 numpy analytic gradient at the receipted theta:
                 max_i |mid(g_arb) - g_ref| / max(1, |g_ref_i|) <= 1e-12,
                 plus the nll absolute deviation. The reference carries its
                 own float64 rounding — the measured number is receipted
                 verbatim either way.
  gate_fd_hessian  the assembled arb Hessian column vs central FD of the
                 arb gradient (exact power-of-two step, matched precision)
                 on sampled columns — the Hessian has NO numpy reference;
                 this is the internal honesty gate, mutation-tested in the
                 battery.
"""
import hashlib
import json
import os
import sys
import time

import numpy as np

# Reference-side inputs: the originating study's reference code and fits
# (not distributed with this package). Set CLINCH_CODE_DIR (the reference
# model code) and CLINCH_REFERENCE_DIR (the reference fits and their stored
# results) to your copies; loaders refuse loudly when unset.


def _env_first(*names):
    """First non-empty value among the named environment variables."""
    for nm in names:
        v = os.environ.get(nm, "")
        if v:
            return v
    return ""


CODE_DIR = _env_first("CLINCH_CODE_DIR", "CLINCH_LANE_CODE")
REFERENCE_DIR = _env_first("CLINCH_REFERENCE_DIR", "CLINCH_RECEIPTS_ROOT",
                           "CLINCH_LANE_RECEIPTS")
DEFAULT_CACHE = os.path.join(REFERENCE_DIR, "clinch_cache") if REFERENCE_DIR else ""


def _require_reference(what):
    if not (CODE_DIR and REFERENCE_DIR):
        raise RuntimeError(
            f"SKIP {what}: CLINCH_CODE_DIR / CLINCH_REFERENCE_DIR unset — "
            "these legs read the originating study's reference code and "
            "fits (not distributed with this package); point the env vars "
            "at your copies")

__all__ = ["load", "extract_arrays", "constants", "gate_gradient",
           "gate_fd_hessian", "CODE_DIR", "REFERENCE_DIR"]


def _reference_code():
    for p in (CODE_DIR, os.path.join(CODE_DIR, "jaxlike")):
        if p not in sys.path:
            sys.path.insert(0, p)
    import refbridge
    import likelihood as L
    import v3lib as V
    return refbridge, L, V


def extract_arrays(data, V):
    """The flat numpy arrays the oracle consumes (v3.1 cell grain)."""
    return dict(
        X=np.asarray(data.X, dtype=float),
        C=np.asarray(data.C, dtype=float),
        D_site=np.asarray(V.site_dryness(data.C), dtype=float),
        group_of=np.asarray(data.group_of, dtype=int),
        gen_of=np.asarray(data.species["genus_of"], dtype=int),
        fam_of=np.asarray(data.species["family_of"], dtype=int),
        sv_sp=data.sv_sp, sv_k=data.sv_k, sv_site=data.sv_site,
        sv_n=data.sv_n, sv_die=data.sv_die, sv_dt5=data.sv_dt5,
        sv_zbar=data.sv_zbar, sv_xdry=data.sv_xdry,
        av_sp=data.av_sp, av_k=data.av_k, av_site=data.av_site,
        av_n=data.av_n, av_up=data.av_up, av_dt5=data.av_dt5,
        av_zbar=data.av_zbar, av_xdry=data.av_xdry,
        rc_sp=data.rc_sp, rc_site=data.rc_site, rc_count=data.rc_count,
        rc_dt5=data.rc_dt5, rc_zbar=data.rc_zbar,
        rc_adults=np.asarray(data.rc_adults, dtype=float))


def constants(data, dials, L, V):
    pc = L.compadre_centers_cached()
    return dict(
        slice_bounds={nm: (sl.start, sl.stop)
                      for nm, sl in dials.slices.items()},
        n_total=int(dials.n_total),
        S=int(data.S), NG=int(data.NG), NF=int(data.NF),
        families=list(data.species["families"]),
        prior_sd={k: float(v) for k, v in L.PRIOR_SD.items()},
        pc_surv_logit=float(pc["surv_logit"]),
        pc_fec_log=float(pc["fec_log"]),
        pc_log_sig=float(L.PRIOR_CENTER["log_sig"]),
        pc_log_phi=float(L.PRIOR_CENTER["log_phi"]),
        pc_log_nu=float(L.PRIOR_CENTER["log_nu"]),
        fec_offsets=[float(x) for x in V.FECUNDITY_STAGE_OFFSETS])


def _sha(x):
    return hashlib.sha256(x).hexdigest()[:16]


def load(tag, cache_dir=None, use_cache=True, verbose=False):
    """(arrs, consts, theta, provenance) for 'fit13' or 'fit17'.

    Fresh path: the reference runners' load_data + registry assertion via
    refbridge.load_assembly (minutes). Cache path: the sha-pinned npz this
    function wrote on a previous fresh run — provenance records which."""
    _require_reference(f"load({tag!r})")
    if cache_dir is None:
        cache_dir = DEFAULT_CACHE
    assert tag in ("fit13", "fit17"), tag
    os.makedirs(cache_dir, exist_ok=True)
    npz_path = os.path.join(cache_dir, f"{tag}_arrays.npz")
    js_path = os.path.join(cache_dir, f"{tag}_consts.json")
    if use_cache and os.path.exists(npz_path) and os.path.exists(js_path):
        blob = open(npz_path, "rb").read()
        meta = json.load(open(js_path))
        assert meta["npz_sha16"] == _sha(blob), \
            f"{tag} cache drifted — delete {npz_path} and re-run fresh"
        z = np.load(npz_path)
        arrs = {k: z[k] for k in z.files if k != "theta"}
        theta = z["theta"]
        prov = dict(meta["provenance"], path="cache", npz=npz_path,
                    npz_sha16=meta["npz_sha16"])
        return arrs, meta["consts"], theta, prov
    t0 = time.time()
    refbridge, L, V = _reference_code()
    data, dials, theta, nll_rec = refbridge.load_assembly(tag,
                                                          verbose=verbose)
    arrs = extract_arrays(data, V)
    consts = constants(data, dials, L, V)
    prov = dict(path="fresh_load", tag=tag,
                loader="refbridge.load_assembly (reference runners' load_data; "
                       "registry keys asserted vs the receipt)",
                nll_record=float(nll_rec),
                theta_sha16=_sha(np.ascontiguousarray(theta).tobytes()),
                wall_s=round(time.time() - t0, 1),
                stamp=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
    np.savez_compressed(npz_path + ".part.npz", theta=theta, **arrs)
    os.replace(npz_path + ".part.npz", npz_path)
    blob = open(npz_path, "rb").read()
    json.dump(dict(consts=consts, provenance=prov,
                   npz_sha16=_sha(blob)),
              open(js_path + ".part", "w"), indent=1)
    os.replace(js_path + ".part", js_path)
    # keep the live reference objects for gate use on the fresh path
    prov["_live"] = (data, dials)
    return arrs, consts, theta, prov


# ---------------------------------------------------------------------------
# identity gates
# ---------------------------------------------------------------------------
def gate_gradient(oracle, theta, tag, prec=256, live=None):
    """arb oracle vs the receipted numpy analytic gradient at theta."""
    from flint import arb
    from baller.hygiene import ctx_guard
    if live is None:
        refbridge, L, V = _reference_code()
        data, dials, th_r, _ = refbridge.load_assembly(tag)
        assert np.allclose(th_r, theta), "theta drifted vs receipt"
    else:
        data, dials = live
        refbridge, L, V = _reference_code()
    nll_ref, g_ref = L.nll_grad(np.asarray(theta, dtype=float), data, dials)
    t0 = time.time()
    with ctx_guard(prec=prec):
        th = [arb(float(t)) for t in theta]
        nll_b, grad, _ = oracle.evaluate(th, order2=False)
        nll_mid = float(nll_b.mid())
        diffs = np.array([abs(float(grad[i].mid()) - g_ref[i])
                          / max(1.0, abs(g_ref[i]))
                          for i in range(len(g_ref))])
    return dict(gate="gradient_identity", tag=tag, prec_bits=prec,
                bar=1e-12, max_rel=float(diffs.max()),
                argmax=int(diffs.argmax()),
                nll_ref=float(nll_ref), nll_arb_mid=nll_mid,
                nll_abs_dev=abs(nll_mid - float(nll_ref)),
                ok=bool(diffs.max() <= 1e-12),
                wall_s=round(time.time() - t0, 1),
                metric="max_i |mid(g_arb)-g_ref| / max(1,|g_ref_i|)")


def _mp_grad_coord(data, dials, V, theta, j, dps=50):
    """Independent high-precision reference for gradient coordinate j:
    mpmath transliteration of likelihood.py's HAND-DERIVED gradient chain
    (dll_dpdie/dq/dlogmu etc.) — a derivation-independent cross-check of
    the jet-based arb oracle (the two share only the model definition)."""
    import mpmath as mp
    mp.mp.dps = dps
    v = dials.view(np.asarray(theta, dtype=float))
    sl = dials.slices
    gen_of = data.species["genus_of"]
    fam_of = data.species["family_of"]
    D_site = V.site_dryness(data.C)
    F = lambda x: mp.mpf(float(x))
    sig = [mp.e ** F(v["log_sig"][i]) for i in range(9)]

    def name_of(j):
        for nm, s in sl.items():
            if s.start <= j < s.stop:
                return nm, j - s.start
        raise KeyError(j)

    nm, off = name_of(j)
    tot = mp.mpf(0)

    def eta_cell(proc, c):
        if proc == 0:
            sp, k, st = int(data.sv_sp[c]), int(data.sv_k[c]), int(data.sv_site[c])
            u = (sig[2] * F(v["z_fam_surv"][fam_of[sp]])
                 + sig[1] * F(v["z_gen_surv"][gen_of[sp]])
                 + sig[0] * F(v["z_sp_surv"][sp]))
            return (F(v["a_surv"][k]) + F(data.X[sp, 0]) * F(v["b_surv"][0])
                    + F(data.X[sp, 1]) * F(v["b_surv"][1])
                    + F(v["b_sxk"][0]) * F(data.X[sp, 0]) * F((k - 3.0) / 2.0)
                    + sum(F(data.C[st, i]) * F(v["g_surv"][i]) for i in range(3))
                    + F(v["w_surv"][0]) * F(data.sv_zbar[c])
                    + F(v["w_surv"][1]) * F(D_site[st]) * F(data.sv_zbar[c])
                    + F(v["w_surv"][2]) * F(data.sv_xdry[c]) + u)
        if proc == 1:
            sp, k, st = int(data.av_sp[c]), int(data.av_k[c]), int(data.av_site[c])
            u = (sig[5] * F(v["z_fam_adv"][fam_of[sp]])
                 + sig[4] * F(v["z_gen_adv"][gen_of[sp]])
                 + sig[3] * F(v["z_sp_adv"][sp]))
            return (F(v["a_adv"][k]) + F(data.X[sp, 0]) * F(v["b_gro"][0])
                    + F(data.X[sp, 1]) * F(v["b_gro"][1])
                    + sum(F(data.C[st, i]) * F(v["g_gro"][i]) for i in range(3))
                    + F(v["w_gro"][0]) * F(data.av_zbar[c])
                    + F(v["w_gro"][1]) * F(D_site[st]) * F(data.av_zbar[c])
                    + F(v["w_gro"][2]) * F(data.av_xdry[c]) + u)
        sp, st = int(data.rc_sp[c]), int(data.rc_site[c])
        u = (sig[8] * F(v["z_fam_fec"][fam_of[sp]])
             + sig[7] * F(v["z_gen_fec"][gen_of[sp]])
             + sig[6] * F(v["z_sp_fec"][sp]))
        return (F(v["a_fec"][0]) + F(data.X[sp, 0]) * F(v["b_fec"][0])
                + F(data.X[sp, 1]) * F(v["b_fec"][1])
                + sum(F(data.C[st, i]) * F(v["g_fec"][i]) for i in range(3))
                + F(v["lam"][data.group_of[sp]]) * F(data.rc_zbar[c]) + u)

    def coef(proc, c, sp):
        """deta/dtheta_j for this cell (0 -> skip); reference chain."""
        if proc == 0:
            if nm == "a_surv":
                return 1.0 if int(data.sv_k[c]) == off else 0.0
            if nm == "b_surv":
                return float(data.X[sp, off])
            if nm == "b_sxk":
                return float(data.X[sp, 0]) * (int(data.sv_k[c]) - 3.0) / 2.0
            if nm == "g_surv":
                return float(data.C[int(data.sv_site[c]), off])
            if nm == "w_surv":
                st = int(data.sv_site[c])
                return [float(data.sv_zbar[c]),
                        float(D_site[st]) * float(data.sv_zbar[c]),
                        float(data.sv_xdry[c])][off]
            if nm == "z_sp_surv":
                return sig[0] if sp == off else 0.0
            if nm == "z_gen_surv":
                return sig[1] if gen_of[sp] == off else 0.0
            if nm == "z_fam_surv":
                return sig[2] if fam_of[sp] == off else 0.0
            if nm == "log_sig" and off == 0:
                return sig[0] * F(v["z_sp_surv"][sp])
            if nm == "log_sig" and off == 1:
                return sig[1] * F(v["z_gen_surv"][gen_of[sp]])
            if nm == "log_sig" and off == 2:
                return sig[2] * F(v["z_fam_surv"][fam_of[sp]])
            return 0.0
        if proc == 1:
            if nm == "a_adv":
                return 1.0 if int(data.av_k[c]) == off else 0.0
            if nm == "b_gro":
                return float(data.X[sp, off])
            if nm == "g_gro":
                return float(data.C[int(data.av_site[c]), off])
            if nm == "w_gro":
                st = int(data.av_site[c])
                return [float(data.av_zbar[c]),
                        float(D_site[st]) * float(data.av_zbar[c]),
                        float(data.av_xdry[c])][off]
            if nm == "z_sp_adv":
                return sig[3] if sp == off else 0.0
            if nm == "z_gen_adv":
                return sig[4] if gen_of[sp] == off else 0.0
            if nm == "z_fam_adv":
                return sig[5] if fam_of[sp] == off else 0.0
            if nm == "log_sig" and off == 3:
                return sig[3] * F(v["z_sp_adv"][sp])
            if nm == "log_sig" and off == 4:
                return sig[4] * F(v["z_gen_adv"][gen_of[sp]])
            if nm == "log_sig" and off == 5:
                return sig[5] * F(v["z_fam_adv"][fam_of[sp]])
            return 0.0
        if nm == "a_fec" and off == 0:
            return 1.0
        if nm == "b_fec":
            return float(data.X[sp, off])
        if nm == "g_fec":
            return float(data.C[int(data.rc_site[c]), off])
        if nm == "lam":
            return (float(data.rc_zbar[c])
                    if int(data.group_of[sp]) == off else 0.0)
        if nm == "z_sp_fec":
            return sig[6] if sp == off else 0.0
        if nm == "z_gen_fec":
            return sig[7] if gen_of[sp] == off else 0.0
        if nm == "z_fam_fec":
            return sig[8] if fam_of[sp] == off else 0.0
        if nm == "log_sig" and off == 6:
            return sig[6] * F(v["z_sp_fec"][sp])
        if nm == "log_sig" and off == 7:
            return sig[7] * F(v["z_gen_fec"][gen_of[sp]])
        if nm == "log_sig" and off == 8:
            return sig[8] * F(v["z_fam_fec"][fam_of[sp]])
        return 0.0

    nu = mp.e ** F(v["log_nu"][0])
    phi = mp.e ** F(v["log_phi"][0])
    # survival cells
    if nm in ("a_surv", "b_surv", "b_sxk", "g_surv", "w_surv", "log_nu",
              "z_sp_surv", "z_gen_surv", "z_fam_surv") \
            or (nm == "log_sig" and off <= 2):
        for c in range(len(data.sv_n)):
            sp = int(data.sv_sp[c])
            cf = None if nm == "log_nu" else coef(0, c, sp)
            if cf == 0.0:
                continue
            n = F(data.sv_n[c]); y = F(data.sv_die[c])
            dt5 = F(data.sv_dt5[c])
            eta = eta_cell(0, c)
            p5 = 1 / (1 + mp.e ** (-eta))
            psur = p5 ** dt5
            pdie = 1 - psur
            al, be = pdie * nu, (1 - pdie) * nu
            dal = mp.digamma(y + al) - mp.digamma(al)
            dbe = mp.digamma(n - y + be) - mp.digamma(be)
            if nm == "log_nu":
                dll_dnu = (pdie * dal + (1 - pdie) * dbe
                           - mp.digamma(n + nu) + mp.digamma(nu))
                tot += -dll_dnu * nu
            else:
                dll_dpdie = nu * (dal - dbe)
                dpdie_deta = -psur * dt5 * (1 - p5)
                tot += -(dll_dpdie * dpdie_deta) * cf
    # advance cells
    if nm in ("a_adv", "b_gro", "g_gro", "w_gro",
              "z_sp_adv", "z_gen_adv", "z_fam_adv") \
            or (nm == "log_sig" and 3 <= off <= 5):
        for c in range(len(data.av_n)):
            sp = int(data.av_sp[c])
            cf = coef(1, c, sp)
            if cf == 0.0:
                continue
            n = F(data.av_n[c]); y = F(data.av_up[c])
            dt5 = F(data.av_dt5[c])
            eta = eta_cell(1, c)
            q5 = 1 / (1 + mp.e ** (-eta))
            stay5 = mp.e ** (dt5 * mp.log(1 - q5))
            q = 1 - stay5
            dll_dq = y / q - (n - y) / (1 - q)
            dq_deta = stay5 * dt5 * q5
            tot += -(dll_dq * dq_deta) * cf
    # fecundity cells
    if nm in ("a_fec", "b_fec", "g_fec", "lam", "log_phi",
              "z_sp_fec", "z_gen_fec", "z_fam_fec") \
            or (nm == "log_sig" and off >= 6):
        offs = [mp.mpf(i) for i in range(data.rc_adults.shape[1])]
        a1 = F(v["a_fec"][1])
        for c in range(len(data.rc_count)):
            sp = int(data.rc_sp[c])
            R = F(data.rc_count[c]); dt5 = F(data.rc_dt5[c])
            A_ = sum(F(data.rc_adults[c][i]) * mp.e ** (a1 * offs[i])
                     for i in range(len(offs)))
            eta = eta_cell(2, c)
            mu = dt5 * mp.e ** eta * A_
            dll_dlogmu = R - mu * (R + phi) / (phi + mu)
            gc = -dll_dlogmu
            if nm == "log_phi":
                dll_dphi = (mp.digamma(R + phi) - mp.digamma(phi)
                            + mp.log(phi / (phi + mu)) + 1
                            - (R + phi) / (phi + mu))
                tot += -dll_dphi * phi
            elif nm == "a_fec" and off == 1:
                dA = sum(F(data.rc_adults[c][i]) * offs[i]
                         * mp.e ** (a1 * offs[i]) for i in range(len(offs)))
                tot += gc * dA / A_
            else:
                cf = coef(2, c, sp)
                if cf != 0.0:
                    tot += gc * cf
    # priors (reference formulas)
    import likelihood as L
    pc = L.compadre_centers_cached()
    th = np.asarray(theta, dtype=float)
    if nm.startswith("z_"):
        tot += F(th[j])
    elif nm == "a_surv":
        tot += (F(th[j]) - F(pc["surv_logit"])) / F(L.PRIOR_SD["a_surv"]) ** 2
    elif nm == "a_adv":
        tot += F(th[j]) / F(L.PRIOR_SD["a_adv"]) ** 2
    elif nm == "a_fec":
        c0 = pc["fec_log"] if off == 0 else 0.0
        sd = L.PRIOR_SD["a_fec0"] if off == 0 else L.PRIOR_SD["a_fec1"]
        tot += (F(th[j]) - F(c0)) / F(sd) ** 2
    elif nm in ("b_surv", "b_gro", "b_fec"):
        tot += F(th[j]) / F(L.PRIOR_SD["b"]) ** 2
    elif nm == "b_sxk":
        tot += F(th[j]) / F(L.PRIOR_SD["b_sxk"]) ** 2
    elif nm in ("g_surv", "g_gro", "g_fec"):
        tot += F(th[j]) / F(L.PRIOR_SD["g"]) ** 2
    elif nm == "log_sig":
        tot += (F(th[j]) - F(L.PRIOR_CENTER["log_sig"])) \
            / F(L.PRIOR_SD["log_sig"]) ** 2
    elif nm == "lam":
        tot += F(th[j]) / F(L.PRIOR_SD["lam"]) ** 2
    elif nm in ("w_surv", "w_gro"):
        tot += F(th[j]) / F(L.PRIOR_SD["w"]) ** 2
    elif nm == "log_phi":
        tot += (F(th[j]) - F(L.PRIOR_CENTER["log_phi"])) \
            / F(L.PRIOR_SD["log_phi"]) ** 2
    elif nm == "log_nu":
        tot += (F(th[j]) - F(L.PRIOR_CENTER["log_nu"])) \
            / F(L.PRIOR_SD["log_nu"]) ** 2
    return tot, f"{nm}[{off}]"


def gate_gradient_adjudicated(oracle, theta, tag, live=None, prec=256,
                              top_k=5, dps=50):
    """The full identity gate, both registers receipted:
      (i) the registered fit-theta metric vs the float64 reference,
          number receipted VERBATIM (its bar is float64-noise-limited);
      (ii) ADJUDICATION at matched precision: the top-K disagreement
          coordinates re-derived in mpmath at `dps` digits via the
          reference's own hand-derived chain — |arb - mp|/max(1,|mp|)
          must meet the 1e-12 bar (this register carries the identity
          claim); |ref_float64 - mp| is receipted alongside (the
          reference's own rounding, adjudicated)."""
    import mpmath as mp
    if live is None:
        refbridge, L, V = _reference_code()
        data, dials, th_r, _ = refbridge.load_assembly(tag)
        assert np.allclose(th_r, theta), "theta drifted vs receipt"
        live = (data, dials)
    else:
        refbridge, L, V = _reference_code()
        data, dials = live
    base = gate_gradient(oracle, theta, tag, prec=prec, live=live)
    from flint import arb
    from baller.hygiene import ctx_guard
    with ctx_guard(prec=prec):
        th = [arb(float(t)) for t in theta]
        _, grad, _ = oracle.evaluate(th, order2=False)
        _, g_ref = L.nll_grad(np.asarray(theta, dtype=float), data, dials)
        diffs = np.array([abs(float(grad[i].mid()) - g_ref[i])
                          / max(1.0, abs(g_ref[i]))
                          for i in range(len(g_ref))])
        top = np.argsort(diffs)[::-1][:top_k]
        adj = []
        for j in top:
            mpv, label = _mp_grad_coord(data, dials, V, theta, int(j),
                                        dps=dps)
            am = mp.mpf(grad[int(j)].mid().str(40, radius=False))
            arb_vs_mp = float(abs(am - mpv) / max(1, abs(mpv)))
            ref_vs_mp = float(abs(mp.mpf(float(g_ref[int(j)])) - mpv))
            adj.append(dict(coord=int(j), name=label,
                            arb_vs_mp_rel=arb_vs_mp,
                            ref_float64_vs_mp_abs=ref_vs_mp,
                            registered_metric_here=float(diffs[int(j)])))
    ok_adj = all(a["arb_vs_mp_rel"] <= 1e-12 for a in adj)
    return dict(registered=base, adjudication=adj,
                adjudication_bar=1e-12, adjudication_ok=bool(ok_adj),
                verdict=("IDENTITY-ESTABLISHED (adjudicated at matched "
                         "precision)" if ok_adj else "IDENTITY-FAILED"),
                note=("registered float64-register number receipted "
                      "verbatim; the identity claim rides the matched-"
                      "precision adjudication register"))


def gate_fd_hessian(oracle, theta, cols, prec=256, h_pow=-20,
                    hess_mats=None):
    """Assembled arb Hessian columns vs central FD of the arb gradient
    (exact step 2**h_pow). cols: theta indices to check. Returns the
    receipt with per-column max rel diffs (bar 1e-6: FD truncation)."""
    from flint import arb
    from baller.hygiene import ctx_guard
    h = 2.0 ** h_pow
    part = oracle.part
    t0 = time.time()
    with ctx_guard(prec=prec):
        if hess_mats is None:
            th = [arb(float(t)) for t in theta]
            _, _, acc = oracle.evaluate(th, order2=True)
            D, B, G = acc.to_mats()
        else:
            D, B, G = hess_mats
        out = {}
        for j in cols:
            thp = [arb(float(t)) for t in theta]
            thm = [arb(float(t)) for t in theta]
            thp[j] = thp[j] + h
            thm[j] = thm[j] - h
            _, gp, _ = oracle.evaluate(thp, order2=False)
            _, gm, _ = oracle.evaluate(thm, order2=False)
            fd = [(gp[i] - gm[i]) / (2 * h) for i in range(len(theta))]
            # assembled column j
            bj, oj = part.block_of[j], part.off_of[j]
            worst = 0.0
            for i in range(len(theta)):
                bi, oi = part.block_of[i], part.off_of[i]
                if bi < 0 and bj < 0:
                    hij = G[oi, oj]
                elif bi >= 0 and bj >= 0:
                    hij = D[bi][oi, oj] if bi == bj else arb(0)
                elif bi >= 0:
                    hij = B[bi][oi, oj]
                else:
                    hij = B[bj][oj, oi]
                d = abs(float(hij.mid()) - float(fd[i].mid())) \
                    / max(1.0, abs(float(hij.mid())))
                worst = max(worst, d)
            out[int(j)] = worst
    return dict(gate="fd_hessian", cols={int(k): float(v)
                                         for k, v in out.items()},
                h=h, prec_bits=prec, bar=1e-6,
                ok=bool(max(out.values()) <= 1e-6),
                wall_s=round(time.time() - t0, 1))
