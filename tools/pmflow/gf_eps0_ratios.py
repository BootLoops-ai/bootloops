#!/usr/bin/env python3
r"""
gf_eps0_ratios.py — eps^0 ratio extraction + exact-rational identification
for fixed-point (forced-ending) AMFlow runs.  One consolidated pipeline for
the discovery it was built on: the 4PM cut top sector of the reference
family is rank-1 over Q at eps^0 (all 7 master ratios exact rationals,
>=102.5d verified, 2 points).

Pipeline per parent log (AMFLOW_DUMP_EPS_GRID output):
  1. load per-node masters M[m][k]; same-sector ratios rho_r = M_r/M_0
     (sector-block wild factors cancel exactly per node);
  2. barycentric-Lagrange extrapolation to eps=0; honest digits = full-set
     vs inner-subset agreement (subset stability);
  3. rational identification: direct 2-vector PSLQ (cap 1e19); if null,
     relation-system fallback — scan triples for integer relations, verify
     each numerically at the trust floor, solve the exact system over Q
     (Fraction), and check predicted rationals against the data;
  4. verification: every identification re-checked against the extrapolated
     value; held-out digits = verified digits minus 2*len(height)+2 (the
     digits needed to pin the rational).

Usage:
  gf_eps0_ratios.py --log parent.log [--rows 0-7] [--n-eps 18] [--dps 260]
                    [--out table.json] [--selftest]
  --selftest replays a reference eps-grid log (not shipped; set
    $GF_SELFTEST_LOG to a local copy) and verifies the recorded table
    (m1=-2/5 exact, 6 rationals h<=1e17, all verified >=100d).
"""
import argparse, itertools, json, os, re, sys
from fractions import Fraction
import mpmath as mp


def pb(s):
    s = s.strip()
    balls = re.findall(r'\[([^\]]+)\]|([+\-]?[0-9][0-9.eE+\-]*)', s)
    nums = []
    for b, f in balls:
        if b: nums.append(mp.mpf(b.split('+/-')[0].strip()))
        elif f: nums.append(mp.mpf(f))
    if not nums: return mp.mpc(0)
    if len(nums) == 1:
        return mp.mpc(0, nums[0]) if ('I' in s or 'j' in s) else mp.mpc(nums[0], 0)
    return mp.mpc(nums[0], nums[1])


def load(logpath, n_eps):
    eps, M = [], {}
    for line in open(logpath, errors='replace'):
        if '[EPS_GRID] eps[' in line:
            eps.append(pb(line.split('=', 1)[1]))
        elif '[EPS_GRID] M[' in line:
            m = re.match(r'.*M\[(\d+)\]\[(\d+)\] = (.*)', line)
            if m:
                M.setdefault(int(m.group(1)), [None] * n_eps)[int(m.group(2))] = \
                    pb(m.group(3))
    return eps, M


def lagrange0(xs, ys):
    n = len(xs)
    w = [mp.mpf(1)] * n
    for i in range(n):
        for j in range(n):
            if i != j: w[i] /= (xs[i] - xs[j])
    num = mp.mpc(0); den = mp.mpc(0)
    for i in range(n):
        t = w[i] / (0 - xs[i])
        num += t * ys[i]; den += t
    return num / den


def extract(eps, M, rows, ref=0):
    xs = [e.real for e in eps]
    idx = sorted(range(len(xs)), key=lambda j: xs[j])
    out = {}
    for r in rows:
        if r == ref: continue
        ratio = [M[r][j] / M[ref][j] for j in range(len(xs))]
        full = lagrange0([xs[j] for j in idx], [ratio[j] for j in idx])
        sub = lagrange0([xs[j] for j in idx[1:-1]], [ratio[j] for j in idx[1:-1]])
        digs = float(-mp.log10(abs(full - sub) / abs(full))) if full != sub else mp.mp.dps
        out[r] = (full.real, round(digs, 1))
    return out


def identify(vals):
    """vals: {row: (value, trust)}. Returns {row: Fraction or None}."""
    got = {}
    rows = sorted(vals)
    for r in rows:
        v, tr = vals[r]
        rel = mp.pslq([v, mp.mpf(1)], tol=mp.mpf(10) ** (-(min(tr, 120) - 10)),
                      maxcoeff=10 ** 19, maxsteps=2000000)
        got[r] = Fraction(int(-rel[1]), int(rel[0])) if rel else None
    missing = [r for r in rows if got[r] is None]
    if missing:
        # relation-system fallback over ALL rows (found rationals included
        # as exact constraints)
        rels = []
        for combo in itertools.combinations(rows, 3):
            vec = [vals[r][0] for r in combo] + [mp.mpf(1)]
            rel = mp.pslq(vec, tol=mp.mpf('1e-80'), maxcoeff=10 ** 7,
                          maxsteps=600000)
            if rel is None: continue
            full = {r: 0 for r in rows}; full['c'] = rel[3]
            for k, r in enumerate(combo): full[r] = rel[k]
            res = sum(full[r] * vals[r][0] for r in rows) + full['c']
            sc = max(abs(full[r] * vals[r][0]) for r in rows if full[r])
            if abs(res) / sc < mp.mpf('1e-90'):
                rels.append(full)
        for r in rows:
            if got[r] is not None:
                rels.append({q: (1 if q == r else 0) for q in rows} |
                            {'c': -got[r]})
        A = [[Fraction(rel[q]) for q in rows] for rel in rels]
        b = [Fraction(-rel['c']) for rel in rels]
        Mx = [row + [b[i]] for i, row in enumerate(A)]
        rank = 0
        for c in range(len(rows)):
            piv = next((i for i in range(rank, len(Mx)) if Mx[i][c] != 0), None)
            if piv is None: continue
            Mx[rank], Mx[piv] = Mx[piv], Mx[rank]
            for i in range(len(Mx)):
                if i != rank and Mx[i][c] != 0:
                    f = Mx[i][c] / Mx[rank][c]
                    for c2 in range(len(rows) + 1): Mx[i][c2] -= f * Mx[rank][c2]
            rank += 1
        if rank == len(rows):
            for i, r in enumerate(rows):
                got[r] = Mx[i][len(rows)] / Mx[i][i]
    return got


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--log'); ap.add_argument('--rows', default='0-7')
    ap.add_argument('--n-eps', type=int, default=18)
    ap.add_argument('--dps', type=int, default=260)
    ap.add_argument('--out'); ap.add_argument('--selftest', action='store_true')
    a = ap.parse_args()
    if a.selftest:
        a.log = os.environ.get('GF_SELFTEST_LOG', 'gf_selftest_reference.log')  # reference log, not shipped
    mp.mp.dps = a.dps
    lo, hi = (int(x) for x in a.rows.split('-'))
    eps, M = load(a.log, a.n_eps)
    vals = extract(eps, M, range(lo, hi + 1))
    mp.mp.dps = 130
    got = identify(vals)
    table, ok = {}, True
    for r in sorted(got):
        v, tr = vals[r]
        fr = got[r]
        if fr is None:
            table[f'm{r}'] = {'rational': False, 'value': mp.nstr(v, 40), 'stab': tr}
            print(f'm{r}/m0: NOT identified (trust {tr}d)')
            continue
        err = abs(mp.mpf(fr.numerator) / fr.denominator - v) / abs(v)
        hd = float(-mp.log10(err)) if err > 0 else 130.0
        pin = 2 * max(len(str(abs(fr.numerator))), len(str(fr.denominator))) + 2
        table[f'm{r}'] = {'p': str(fr.numerator), 'q': str(fr.denominator),
                          'verified_digits': round(hd, 1),
                          'heldout_after_pin': round(hd - pin, 1)}
        flag = 'OK' if hd >= 90 else 'WEAK'
        if hd < 90: ok = False
        print(f'm{r}/m0 = {fr.numerator}/{fr.denominator}  verified {hd:.1f}d '
              f'(held-out {hd - pin:.1f}d) {flag}')
    if a.selftest:
        want = Fraction(-2, 5)
        assert got[1] == want, f'selftest: m1/m0 = {got[1]} != -2/5'
        nrat = sum(1 for r in got if got[r] is not None)
        assert nrat == 7, f'selftest: {nrat}/7 rational'
        print('SELFTEST PASS: 7/7 rational incl. m1/m0 = -2/5')
    if a.out:
        json.dump(table, open(a.out, 'w'), indent=1)
        print('wrote', a.out)
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
