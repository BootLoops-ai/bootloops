"""L3: the flat-valley control, run on the reference fit artifacts (not
shipped; CLINCH_REFERENCE_DIR + CLINCH_CODE_DIR, the adapter's env gates).

fit16 stall optimum = CKPT_STAGE1.npz (gmax 0.026999, the stage-1
checkpoint); THETA_STAGE1.npz is a later tail-polished point, receipted
here as an addendum."""
import sys, time, json
import numpy as np
import os as _os
_TOOLS = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, _os.path.join(_TOOLS, "clinch"))
sys.path.insert(0, _os.path.join(_TOOLS, "baller"))
import clinch; clinch.verify(quiet=True)
from clinch import adapter_v31 as A, engine
from clinch.oracle_v36 import ModelV36Oracle, N_GLOBAL36
import baller
from baller.certify import block_krawczyk as bk
from baller.hygiene import ctx_guard
from flint import arb

RCPT = A.REFERENCE_DIR
if not RCPT:
    raise SystemExit('CLINCH_REFERENCE_DIR unset — path to the reference results tree is required')
PREC = 192
arrs, consts13, th13, prov = A.load('fit13')
refbridge, L, V = A._reference_code()
data, dials13, _, _ = refbridge.load_assembly('fit13')
d36 = L.DialSpace(data.S, data.NG, data.NF, n_storm=0,
                  n_weather=V.V31_WEATHER_DIALS, n_v36=V.V36_SURV_DIALS)
consts36 = dict(consts13)
consts36['slice_bounds'] = {nm: (sl.start, sl.stop)
                            for nm, sl in d36.slices.items()}
consts36['n_total'] = int(d36.n_total)
th16_stall = np.load(f'{RCPT}/fit16_null/CKPT_STAGE1.npz')['theta']
th16_rec = np.load(f'{RCPT}/fit16_null/THETA_STAGE1.npz')['theta']
th18 = np.load(f'{RCPT}/fit18a_null/THETA_STAGE1.npz')['theta']
T_ref = json.load(open(f'{RCPT}/fit18a_null/FIT18A_SUMMARY.json'))['pin']['T_ref_fit13']
import fit18_common as F18
_, g_full = L.nll_grad(th18.astype(float), data, d36)
_, gT = F18.latents_zero_adv_total(th18.astype(float), data, d36,
                                   want_grad=True)
slz = d36.slices['z_sp_adv']
lam_ext = -(float(g_full[slz].sum()) / float(gT[slz].sum()))

def name_of(j_real):
    for nm, s in d36.slices.items():
        if s.start <= j_real < s.stop:
            return f'{nm}[{j_real - s.start}]'
    return f'?{j_real}'

def rownorm_scales(Ht):
    D, B, G = Ht
    Dm = [np.abs(bk._mid_np(M)) for M in D]
    Bm = [np.abs(bk._mid_np(M)) for M in B]
    Gm = np.abs(bk._mid_np(G))
    border = Gm.max(axis=1)
    for Bi in Bm:
        if Bi.size:
            border = np.maximum(border, Bi.max(axis=0))
    sz = [1.0 / np.sqrt(np.maximum(
        np.maximum(Di.max(axis=1), Bi.max(axis=1) if Bi.size else 0),
        1e-6)) for Di, Bi in zip(Dm, Bm)]
    return sz, 1.0 / np.sqrt(np.maximum(border, 1e-6))

def analyze(tag, oracle, ext):
    part = oracle.part
    center = ([[float(ext[t]) for t in i] for i in part.blocks],
              [float(ext[t]) for t in part.border_idx])
    with ctx_guard(prec=PREC):
        Ft = oracle.F(bk.box_around(center, 0.0))
        Ht = oracle.H(bk.box_around(center, 0.0))
    cz, cg = engine._newton_step(bk, Ft, Ht)
    sz, sg = rownorm_scales(Ht)
    cflat = np.zeros(oracle.n_ext); sflat = np.ones(oracle.n_ext)
    for bi, idx in enumerate(part.blocks):
        cflat[idx] = cz[bi]; sflat[idx] = sz[bi]
    cflat[part.border_idx] = cg; sflat[part.border_idx] = sg
    scaled = np.abs(cflat) / sflat
    top = np.argsort(scaled)[::-1][:5]
    naux = oracle.naux
    names = []
    for j in top:
        if j < N_GLOBAL36:
            names.append((name_of(int(j)), float(scaled[j])))
        elif j < N_GLOBAL36 + naux:
            names.append((f'aux[{int(j) - N_GLOBAL36}]', float(scaled[j])))
        else:
            names.append((name_of(int(j) - naux), float(scaled[j])))
    gz = d36.slices['z_sp_adv']
    gauge_ids = np.arange(gz.start + naux, gz.stop + naux)
    gvec = cflat[gauge_ids]
    proj = abs(gvec.sum()) / np.sqrt(len(gvec)) / max(
        np.linalg.norm(cflat), 1e-300)
    frac_zspadv = float((np.abs(cflat[gauge_ids]) ** 2).sum()
                        / max((np.abs(cflat) ** 2).sum(), 1e-300))
    return dict(tag=tag, center=center, Ft=Ft, sz=sz, sg=sg,
                scaled_inf=float(scaled.max()),
                top_scaled=names, gauge_proj=float(proj),
                step_mass_frac_z_sp_adv=frac_zspadv)

def margins_at(oracle, info, ru):
    part = oracle.part
    r = ([[ru * float(v) for v in s] for s in info['sz']],
         [ru * float(v) for v in info['sg']])
    with ctx_guard(prec=PREC):
        Hb = oracle.H(bk.box_around(info['center'], r))
    kw = bk.block_krawczyk(oracle, info['center'], r, PREC,
                           fh=(info['Ft'], Hb))
    out = dict(r_unit=ru, verdict=kw['verdict'],
               max_margin=kw.get('max_margin'),
               border_margin=kw.get('margins', {}).get('border'),
               failed=kw.get('failed'), reason=kw.get('reason', '')[:160],
               clip_census=dict(getattr(oracle, 'last_clip_census', {})))
    if (kw.get('failed') or '').startswith('block_'):
        out['failed_named'] = part.block_names[int(kw['failed'].split('_')[1])]
    return out

R = {}
# ---- fit18a: deep polish then certify ---------------------------------
o18 = ModelV36Oracle(arrs, consts36, pinned=True, T_ref=T_ref)
ext = o18.init_multipliers(th18, lam=lam_ext, prec=256)
th = np.array(ext, dtype=float)
part = o18.part
for it in range(5):
    center = ([[float(th[t]) for t in i] for i in part.blocks],
              [float(th[t]) for t in part.border_idx])
    with ctx_guard(prec=PREC):
        Ft = o18.F(bk.box_around(center, 0.0))
        Ht = o18.H(bk.box_around(center, 0.0))
    cz, cg = engine._newton_step(bk, Ft, Ht)
    sz, sg = rownorm_scales(Ht)
    sc = max([float(np.max(np.abs(c)/s)) for c, s in zip(cz, sz) if len(c)]
             + [float(np.max(np.abs(cg)/sg))])
    print(f'fit18a polish it{it}: scaled {sc:.3g}', flush=True)
    if sc < 1e-12:
        break
    for bi, idx in enumerate(part.blocks):
        th[idx] -= cz[bi]
    th[part.border_idx] -= cg
R['fit18a_polish_scaled_final'] = sc
i18 = analyze('fit18a_pinned_polished', o18, th)
R['fit18a_analysis'] = {k: v for k, v in i18.items()
                        if k not in ('center', 'Ft', 'sz', 'sg')}
ru18 = max(2.0 * i18['scaled_inf'], 3e-10)
r18 = None
for ru in (ru18, 2 * ru18, 4 * ru18):
    r18 = margins_at(o18, i18, ru)
    print('fit18a rung', ru, r18['verdict'], r18['max_margin'], flush=True)
    if r18['verdict'] == 'CERTIFIED':
        break
R['fit18a_cert'] = r18
R['fit18a_dist_from_receipted'] = float(np.max(np.abs(
    th - np.array(o18.init_multipliers(th18, lam=lam_ext, prec=256)))))

# ---- fit16 STALL optimum (the registered artifact) --------------------
o16 = ModelV36Oracle(arrs, consts36, pinned=False)
ext16 = o16.init_multipliers(th16_stall, prec=256)
i16 = analyze('fit16_stall_summary_optimum', o16, ext16)
R['fit16_stall_analysis'] = {k: v for k, v in i16.items()
                             if k not in ('center', 'Ft', 'sz', 'sg')}
ru16 = 2.0 * i16['scaled_inf']
R['fit16_stall_at_coverage'] = margins_at(o16, i16, ru16)
print('fit16 stall @ coverage', ru16, R['fit16_stall_at_coverage'], flush=True)
R['fit16_stall_at_common'] = margins_at(o16, i16, r18['r_unit'])
print('fit16 stall @ common ru', R['fit16_stall_at_common'], flush=True)

# ---- fit16 record theta (addendum) ------------------------------------
ext16r = o16.init_multipliers(th16_rec, prec=256)
i16r = analyze('fit16_record_theta', o16, ext16r)
R['fit16_record_analysis'] = {k: v for k, v in i16r.items()
                              if k not in ('center', 'Ft', 'sz', 'sg')}
R['fit16_record_at_coverage'] = margins_at(
    o16, i16r, max(2.0 * i16r['scaled_inf'], 3e-10))
print('fit16 record @ coverage', R['fit16_record_at_coverage'], flush=True)

# ---- the statistic ----------------------------------------------------
m16c = R['fit16_stall_at_common']['max_margin']
R['comparison'] = dict(
    fit16_stall_at_coverage=dict(ru=ru16,
                                 **{k: v for k, v in
                                    R['fit16_stall_at_coverage'].items()
                                    if k != 'clip_census'}),
    fit18a_margin=r18['max_margin'], fit18a_ru=r18['r_unit'],
    defect_ratio_common_radius=(None if m16c is None else
                                m16c / max(r18['max_margin'], 1e-300)),
    common_ru_pair=(m16c, r18['max_margin']),
    flat_valley_geometry=dict(
        gauge_gradient_at_record=7.974078416022947e-11,
        gauge_gradient_at_stall=7.573e-05,
        pin_displacement_a_adv=0.7692, delta_nll=0.65094,
        note='two receipted optima 0.77 apart along the advance-level '
             'direction at 0.65 nll cost — the near-flat direction, '
             'measured'))
json.dump(R, open('L3_RESULTS.json', 'w'), indent=1, default=float)
print(json.dumps(R['comparison'], indent=1, default=float))
print('fit16 stall top scaled-step coords:', R['fit16_stall_analysis']['top_scaled'])
print('fit16 stall z_sp_adv step-mass fraction:',
      R['fit16_stall_analysis']['step_mass_frac_z_sp_adv'])
