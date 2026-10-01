#!/usr/bin/env python3
"""build_sym5.py — CAL-B instrument builder: the exact scalar Sym^5-Legendre
operator (order 6), constructed by exact series elimination from the HYP2F1
letter (f0 = 2F1(1/2,1/2;1;lam), c_n = binom(2n,n)^2/16^n), and VERIFIED:
(i) annihilates f0^5 to depth, (ii) indicial roots at 0 are {0,0,0,0,0,0}
(full MUM J(6) candidate; jordan measured downstream by the transport's own
gates).  Emits work/SYM5_OP.json (producer block, optional lint)."""
import json, os, sys, math, subprocess, datetime, hashlib, time
from fractions import Fraction

MEMBER_DIR = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import deform_kit as dk

NTER = 420
t0 = time.time()
# f0^5 exact integer-scaled series: c_n * 16^n integer; y = f0^5: coeff * 16^n
c = [math.comb(2 * n, n) ** 2 for n in range(NTER)]          # c_n * 16^n
def mul16(a, b):
    out = [0] * NTER
    for i, ai in enumerate(a):
        if ai:
            for k in range(NTER - i):
                if b[k]:
                    out[i + k] += ai * b[k]
    return out
f2_ = mul16(c, c)
f4_ = mul16(f2_, f2_)
y16 = mul16(f4_, c)          # y_n * 16^n exact
y = [Fraction(v, 16 ** n) for n, v in enumerate(y16)]

# fit sum_{i<=6, d<=dcap} b_{i,d} * n^i * lam^d annihilating y:
# coefficient of lam^m: sum_{i,d} b_{i,d} (m-d)^i y_{m-d} = 0
best = None
for dcap in range(2, 12):
    nun = 7 * (dcap + 1)
    rows = []
    for m in range(0, min(NTER - 2, nun + 60)):
        row = []
        for i in range(7):
            for d in range(dcap + 1):
                row.append(Fraction((m - d) ** i) * y[m - d] if m - d >= 0 else Fraction(0))
        rows.append(row)
    A = [r[:] for r in rows]
    piv = {}
    rr = 0
    for cix in range(nun):
        pr = next((r_ for r_ in range(rr, len(A)) if A[r_][cix] != 0), None)
        if pr is None:
            continue
        A[rr], A[pr] = A[pr], A[rr]
        inv = 1 / A[rr][cix]
        A[rr] = [x * inv for x in A[rr]]
        for r_ in range(len(A)):
            if r_ != rr and A[r_][cix] != 0:
                f_ = A[r_][cix]
                A[r_] = [x - f_ * z for x, z in zip(A[r_], A[rr])]
        piv[cix] = rr
        rr += 1
    free = [cix for cix in range(nun) if cix not in piv]
    if free:
        fc = free[0]
        vec = [Fraction(0)] * nun
        vec[fc] = Fraction(1)
        for cix, r_ in piv.items():
            vec[cix] = -A[r_][fc]
        den = 1
        for x in vec:
            den = den * x.denominator // math.gcd(den, x.denominator)
        bint = [int(x * den) for x in vec]
        g = 0
        for x in bint:
            g = math.gcd(g, x)
        if g:
            bint = [x // g for x in bint]
        best = (dcap, bint)
        break
assert best, 'no operator found up to dcap 11'
dcap, bint = best
btheta = [[bint[i * (dcap + 1) + d] for d in range(dcap + 1)] for i in range(7)]

# verify annihilation to full depth on exact series
def check_annihilate():
    for m in range(NTER - 1):
        s = Fraction(0)
        for i in range(7):
            for d in range(dcap + 1):
                if m - d >= 0 and btheta[i][d]:
                    s += btheta[i][d] * Fraction((m - d) ** i) * y[m - d]
        if s != 0:
            return m
    return None
bad = check_annihilate()
assert bad is None, ('annihilation fails at', bad)

# theta -> d/dx rows: a_k(x) = x^k * sum_i S(i,k) b_i(x)   (Stirling 2nd kind)
S2 = [[0] * 8 for _ in range(8)]
S2[0][0] = 1
for i in range(1, 8):
    for k in range(1, i + 1):
        S2[i][k] = S2[i - 1][k - 1] + k * S2[i - 1][k]
maxdeg = dcap + 6
arows = []
for k in range(7):
    poly = [0] * (maxdeg + 1)
    for i in range(k, 7):
        if S2[i][k]:
            for d in range(dcap + 1):
                if btheta[i][d]:
                    poly[d + k] += S2[i][k] * btheta[i][d]
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    arows.append(poly)

# indicial check via OpData
class _A:
    pass
op = dk.OpData(arows, 0, 'Sym5Legendre')
W0 = op.Wk_poly[0]
res = op.resonances
rec = {'producer': {'component': 'frobenius-boundary/deform_transport',
                    'script': 'code/build_sym5.py',
                    'script_sha256': dk.sha256_file(os.path.abspath(__file__)),
                    'stamp_utc': dk.utc()},
       'spec': 'CAL-B operator (exact series elimination; deterministic)',
       'construction': {'series': 'f0^5, f0 = 2F1(1/2,1/2;1;lam), c_n = C(2n,n)^2/16^n exact',
                        'n_terms_used': NTER, 'theta_deg_cap_found': dcap,
                        'order': 6},
       'verification': {'annihilates_f0pow5_to': NTER - 2,
                        'resonances_at_0': res,
                        'W0_poly': [str(x) for x in W0],
                        'expected': 'single resonance m=0 with multiplicity 6 (full MUM)'},
       'coeffs': arows,
       'theta_form_coeffs': btheta,
       'wall_s': round(time.time() - t0, 2)}
assert res == {0: 6}, ('unexpected indicial structure', res)
out = os.path.join(MEMBER_DIR, 'work', 'SYM5_OP.json')
os.makedirs(os.path.dirname(out), exist_ok=True)
with open(out, 'w') as f:
    json.dump(rec, f, indent=1)
print('WROTE', out)
dk.lint(out)
print('dcap:', dcap, 'resonances:', res, 'order-6 rows deg:', [len(r_) - 1 for r_ in arows])
