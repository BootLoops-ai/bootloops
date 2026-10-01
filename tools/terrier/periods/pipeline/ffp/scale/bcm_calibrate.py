# Calibrate H_q convention on the Legendre family (alpha=(1/2,1/2), beta=(1,1)):
# brute a_p of y^2 = x(x-1)(x-lambda) vs H variants (sign, conj, lambda<->1/lambda,
# chi(-1) twist). Requires a UNIQUE variant matching ALL lambda at ALL test primes.
import json, numpy as np
from bcm import GaussField, h_coeffs, h_all_fp

def brute_ap_legendre(p):
    chi = np.zeros(p, dtype=np.int64)
    for x in range(1, p):
        chi[x] = 1 if pow(x, (p - 1) // 2, p) == 1 else -1
    ap = {}
    for lam in range(2, p):
        if lam == 1:
            continue
        s = 0
        for x in range(p):
            s += chi[x * (x - 1) % p * ((x - lam) % p) % p]
        ap[lam] = -s
    return ap

def variants(H, p):
    inv = [0] + [pow(t, p - 2, p) for t in range(1, p)]
    chim1 = 1 if p % 4 == 1 else -1
    out = {}
    for s in (1, -1):
        for tw in (1, chim1):
            for u in ('t', '1/t'):
                key = (s, tw, u)
                out[key] = {t: s * tw * H[t if u == 't' else inv[t]]
                            for t in range(2, p)}
    return out

def main():
    # integrality diagnostic of the H table (max distance of the float H
    # values from integers); the match search itself is the two-pass block
    # below, which recomputes H from scratch
    resid_max = 0.0
    for p in (13, 19, 23, 29):
        gf = GaussField(p, 1)
        C = h_coeffs(gf, ['1/2', '1/2'])
        for conj in (False, True):
            Hc = h_all_fp(gf, C, conj=conj)
            resid_max = max(resid_max, float(np.max(np.abs(Hc[2:] - np.round(Hc.real[2:])))))
    # two-pass: collect matches per prime then intersect
    match_sets = []
    for p in (13, 19, 23, 29):
        gf = GaussField(p, 1)
        C = h_coeffs(gf, ['1/2', '1/2'])
        ap = brute_ap_legendre(p)
        ms = set()
        for conj in (False, True):
            H = np.round(h_all_fp(gf, C, conj=conj).real).astype(np.int64)
            for key, vals in variants(H, p).items():
                if all(vals[t] == ap[t] for t in ap if t != 1):
                    ms.add(('conj' if conj else 'noconj',) + key)
        match_sets.append((p, ms))
        print(p, 'matching variants:', sorted(ms))
    inter = set.intersection(*[m for _, m in match_sets])
    print('intersection:', sorted(inter), 'max int residual:', resid_max)
    # calibration record, written to the current directory (an output of this
    # script, not a packaged input; the frozen convention lives in bcm_passsets.py)
    with open('BCM_CALIB.json', 'w') as fh:
        json.dump(dict(primes=[13, 19, 23, 29], intersection=[list(x) for x in inter],
                       max_residual=resid_max,
                       per_prime={str(p): [list(x) for x in sorted(m)]
                                  for p, m in match_sets}), fh, indent=1)
    return inter

if __name__ == '__main__':
    main()
