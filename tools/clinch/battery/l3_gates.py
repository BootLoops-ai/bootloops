"""L3 preflight: v36 oracle identity gates at fit16 + fit18a thetas."""
import sys, time, json
import numpy as np
import os as _os
_TOOLS = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, _os.path.join(_TOOLS, "clinch"))
sys.path.insert(0, _os.path.join(_TOOLS, "baller"))
import clinch; clinch.verify(quiet=True)
from clinch import adapter_v31 as A
from clinch.oracle_v36 import ModelV36Oracle, N_GLOBAL36
from baller.hygiene import ctx_guard
from flint import arb

RCPT = A.REFERENCE_DIR
if not RCPT:
    raise SystemExit('CLINCH_REFERENCE_DIR unset — path to the reference results tree is required')
arrs, consts13, th13, prov = A.load('fit13')

# reference objects (identity gate reference + registry assertions)
refbridge, L, V = A._reference_code()
data, dials13, _, _ = refbridge.load_assembly('fit13')
d36 = L.DialSpace(data.S, data.NG, data.NF, n_storm=0,
                  n_weather=V.V31_WEATHER_DIALS, n_v36=V.V36_SURV_DIALS)
assert d36.n_global == 53 and d36.n_total == 5594
r13 = json.load(open(f'{RCPT}/fit13_null/REGISTRY_FIT13.json'))
for tag, reg in (('fit16', 'REGISTRY_FIT16'), ('fit18a', 'REGISTRY_FIT18A')):
    rr = json.load(open(f'{RCPT}/{tag}_null/{reg}.json'))
    assert rr['keys'] == r13['keys'], f'{tag} registry != fit13 — STOP'
print('registries identical to fit13 assembly: OK', flush=True)

consts36 = dict(consts13)
consts36['slice_bounds'] = {nm: (sl.start, sl.stop)
                            for nm, sl in d36.slices.items()}
consts36['n_total'] = int(d36.n_total)

th16 = np.load(f'{RCPT}/fit16_null/THETA_STAGE1.npz')['theta']
th18 = np.load(f'{RCPT}/fit18a_null/THETA_STAGE1.npz')['theta']

out = {}
for tag, th in (('fit16', th16), ('fit18a', th18)):
    o = ModelV36Oracle(arrs, consts36, pinned=False)
    ext = o.init_multipliers(th, prec=256)
    nll_ref, g_ref = L.nll_grad(th.astype(float), data, d36)
    t0 = time.time()
    with ctx_guard(prec=256):
        thb = [arb(float(t)) for t in ext]
        nll_b, grad, _ = o.evaluate(thb, order2=False)
    # theta rows of F == true composed gradient (m=ubar, mu=-Phi_m)
    diffs = []
    for j in range(d36.n_total):
        je = j if j < N_GLOBAL36 else j + o.naux
        diffs.append(abs(float(grad[je].mid()) - g_ref[j])
                     / max(1.0, abs(g_ref[j])))
    diffs = np.array(diffs)
    aux_res = [abs(float(grad[N_GLOBAL36 + i].mid())) for i in range(4)]
    rec = dict(tag=tag, max_rel=float(diffs.max()),
               argmax=int(diffs.argmax()),
               nll_ref=float(nll_ref), nll_arb_mid=float(nll_b.mid()),
               nll_abs_dev=abs(float(nll_b.mid()) - float(nll_ref)),
               aux_row_residuals=aux_res,
               wall_s=round(time.time() - t0, 1))
    print(json.dumps(rec), flush=True)
    out[tag] = rec

# FD gate on the EXTENDED system at fit16 (cols: a global, log_m0b, a
# latent, m_s, mu_a) — self-consistency of F vs H incl. aux rows
o = ModelV36Oracle(arrs, consts36, pinned=False)
ext = o.init_multipliers(th16, prec=256)
h = 2.0 ** -20
cols = [0, 52, o.ix_ms, o.ix_mua, N_GLOBAL36 + o.naux + 100]
with ctx_guard(prec=256):
    thb = [arb(float(t)) for t in ext]
    _, _, acc = o.evaluate(thb, order2=True)
    D, B, G = acc.to_mats()
    part = o.part
    worst_all = {}
    for j in cols:
        thp = [arb(float(t)) for t in ext]; thm = [arb(float(t)) for t in ext]
        thp[j] = thp[j] + h; thm[j] = thm[j] - h
        _, gp, _ = o.evaluate(thp, order2=False)
        _, gm, _ = o.evaluate(thm, order2=False)
        bj, oj = part.block_of[j], part.off_of[j]
        worst = 0.0
        for i in range(o.n_ext):
            bi, oi = part.block_of[i], part.off_of[i]
            if bi < 0 and bj < 0: hij = G[oi, oj]
            elif bi >= 0 and bj >= 0:
                hij = D[bi][oi, oj] if bi == bj else arb(0)
            elif bi >= 0: hij = B[bi][oi, oj]
            else: hij = B[bj][oj, oi]
            fd = (gp[i] - gm[i]) / (2 * h)
            dv = abs(float(hij.mid()) - float(fd.mid())) \
                / max(1.0, abs(float(hij.mid())))
            worst = max(worst, dv)
        worst_all[int(j)] = worst
        print('FD col', j, 'worst rel', worst, flush=True)
out['fd'] = worst_all
json.dump(out, open('L3_GATES.json', 'w'), indent=1, default=float)
print('L3 gates done')
