# BCM production pass-sets for AESZ4 (h33) and AESZ11 (h34) — the two carded
# hypergeometric ops with known rank-2 points (positive controls, BKSZ 2203.09426):
#   AESZ4: z* = -1/5832 (wt-4 partner 54.4.a.c); AESZ11: z* = -1/432 (180.4.a.e).
# Law (validated, bcm_validate.py): a=-H_p(Mz), s2=+H_{p^2}(Mz), b=(a^2-s2)/(2p),
# FIRE := b integral & weil_circle & split_12 nonempty. Unit-root congruence asserted
# at every z as in-production integrity check. Primes = 1 mod 12 (both ops served).
import json, os, sys, time
import numpy as np
from bcm import GaussField, h_coeffs, h_all_fp
from bcm_validate import OPS, series_trunc_mod_p, horner_all

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from weil import weil_circle
from quartic import split_12

PRIMES = [61, 73, 97, 109, 157, 181, 193, 229, 241, 277, 313, 337]
JOBS = {'AESZ4': ('h33', (-1, 5832)), 'AESZ11': ('h34', (-1, 432))}

def run(op, hname, ctrl, p):
    t0 = time.time()
    alpha, M = OPS[hname]
    H1r = h_all_fp(GaussField(p, 1), h_coeffs(GaussField(p, 1), alpha))
    H2r = h_all_fp(GaussField(p, 2), h_coeffs(GaussField(p, 2), alpha))
    resid = max(float(np.max(np.abs(H1r[1:] - np.round(H1r.real[1:])))),
                float(np.max(np.abs(H2r[1:] - np.round(H2r.real[1:])))))
    assert resid < 0.05, (op, p, 'float certification fail', resid)
    H1 = np.round(H1r.real).astype(np.int64); H2 = np.round(H2r.real).astype(np.int64)
    S = horner_all(series_trunc_mod_p(alpha, M, p), p)
    Minv = pow(M % p, p - 2, p)
    fires, alphas = [], {}
    for z in range(1, p):
        if z == Minv:
            continue
        t = M * z % p
        a, s2 = -int(H1[t]), int(H2[t])
        assert (-a - S[z]) % p == 0, (op, p, z, 'unit-root congruence FAIL')
        if (a * a - s2) % (2 * p):
            continue
        b = (a * a - s2) // (2 * p)
        if not weil_circle(a, b, p):
            continue
        hits = split_12(a, b, p)
        if hits:
            fires.append(z)
            alphas[str(z)] = sorted({al for al, be in hits})
    czr = ctrl[0] * pow(ctrl[1], p - 2, p) % p
    ctrl_hit = dict(z_residue=int(czr), fired=czr in fires)
    if ctrl_hit['fired']:
        t = M * czr % p
        a = -int(H1[t]); b = (a * a - int(H2[t])) // (2 * p)
        ctrl_hit['splits'] = split_12(a, b, p)
    n = p - 2
    return dict(op=op, hname=hname, p=p, status='OK',
                oracle='BCM finite-hypergeometric H_q (complex-exact Gauss sums)',
                label='CONTROLLED', control_point=f'{ctrl[0]}/{ctrl[1]}',
                control=ctrl_hit, n_smooth=n, sing=[int(Minv)],
                fire_count=len(fires), rate=round(len(fires) / n, 4),
                fires=fires, alphas=alphas, includes_fp2_pass=True,
                float_resid=resid, compute_seconds=round(time.time() - t0, 1))

def main():
    out = os.path.join(HERE, 'passsets')
    for op, (hname, ctrl) in JOBS.items():
        for p in PRIMES:
            fn = os.path.join(out, f'{op}.p{p}.json')
            if os.path.exists(fn):
                continue
            res = run(op, hname, ctrl, p)
            with open(fn, 'w') as fh:
                json.dump(res, fh)
            print(op, p, 'fires', res['fire_count'], 'rate', res['rate'],
                  'ctrl', res['control']['fired'], f"{res['compute_seconds']}s",
                  f"resid {res['float_resid']:.1e}", flush=True)

if __name__ == '__main__':
    main()
