#!/usr/bin/env python3
"""ctrl_per_order.py — per-epsilon-order control comparator for m2-DE tower
boundaries vs an AMFlow oracle.

Why this exists: a scalar "min over all orders" metric can report ~16d and be
mistaken for a boundary-precision problem.  The real structure is a jet
truncation-alias ladder: roughly
|log10 eps0| digits are lost per epsilon-order above kmin, so the min-digit
floor is set by the HIGHEST requested order, not by boundary precision.
This tool prints the per-(master,order) breakdown, fits the digit-loss slope,
and emits an ALIAS-LIMITED verdict when the slope matches |log10 eps0|.

Formats:
  towers : json['towers'][str(j)][str(order)] = [re_str, im_str]
  oracle : json[key]['masters'][str(j)][str(order)] = decimal string
"""
import argparse, json, os, sys
from mpmath import mp, mpf, mpc, fabs, log10

CAP = 99.0
SKIP_ABS = mpf('1e-30')


def rel_digits(val, ora):
    d = fabs(val - ora)
    if d == 0:
        return CAP
    return min(CAP, float(-log10(d / fabs(ora))))


def compare(towers_path, oracle_path, oracle_key, m_lo, m_hi, eps0_log10):
    tw = json.load(open(towers_path))
    orc = json.load(open(oracle_path))[oracle_key]['masters']
    if eps0_log10 is None:
        eps0_log10 = tw.get('eps0_log10')
    cells = {}  # (master, order) -> digits
    for j in range(m_lo, m_hi):
        tj, oj = tw['towers'].get(str(j)), orc.get(str(j))
        if tj is None or oj is None:
            continue
        for k_s, ora_s in oj.items():
            if k_s not in tj:
                continue
            ora = mpf(ora_s)
            if fabs(ora) < SKIP_ABS:
                continue
            re_s, im_s = tj[k_s]
            cells[(j, int(k_s))] = rel_digits(mpc(mpf(re_s), mpf(im_s)), ora)

    orders = sorted({k for _, k in cells})
    per_order = {}
    for k in orders:
        col = {j: d for (j, kk), d in cells.items() if kk == k}
        jmin = min(col, key=col.get)
        per_order[k] = dict(min=col[jmin], argmin=jmin, count=len(col))
    per_master = {}
    for j in sorted({j for j, _ in cells}):
        row = {k: d for (jj, k), d in cells.items() if jj == j}
        kw = min(row, key=row.get)
        per_master[j] = dict(worst_order=kw, digits=row[kw])

    # least-squares slope of min-digits vs order
    xs = orders; ys = [per_order[k]['min'] for k in orders]
    n = len(xs); xb = sum(xs) / n; yb = sum(ys) / n
    sxx = sum((x - xb) ** 2 for x in xs)
    slope = sum((x - xb) * (y - yb) for x, y in zip(xs, ys)) / sxx if sxx else 0.0

    verdict, K, D = "NOT alias-shaped — investigate", None, None
    if eps0_log10 is not None and abs(eps0_log10) > 0:
        e = abs(float(eps0_log10))
        if abs(abs(slope) - e) <= 0.35 * e:
            ok = [k for k in orders if per_order[k]['min'] >= 30.0]
            if ok:
                K = max(ok); D = per_order[K]['min']
                verdict = "ALIAS-LIMITED above order %+d" % K
    return dict(towers=towers_path, oracle=oracle_path, oracle_key=oracle_key,
                eps0_log10=eps0_log10, per_order=per_order, per_master=per_master,
                slope=slope, verdict=verdict, usable_max_order=K, usable_digits=D,
                cells={"%d,%d" % jk: d for jk, d in cells.items()})


def report(r):
    print("== per-order (min over masters) ==")
    print("%6s %10s %7s %6s" % ("order", "min_d", "argmin", "count"))
    for k in sorted(r['per_order']):
        p = r['per_order'][k]
        print("%+6d %10.2f %7d %6d" % (k, p['min'], p['argmin'], p['count']))
    print("== per-master worst order ==")
    for j in sorted(r['per_master']):
        p = r['per_master'][j]
        print("  master %2d : worst at order %+d (%.2f d)" % (j, p['worst_order'], p['digits']))
    print("slope of min-digits vs order: %.3f d/order (|eps0_log10|=%s)"
          % (r['slope'], r['eps0_log10']))
    print("VERDICT:", r['verdict'])
    if r['usable_max_order'] is not None:
        print("usable orders <= %+d at %.2f d" % (r['usable_max_order'], r['usable_digits']))


def selftest():
    # selftest fixtures (reference tower pairs) are not shipped in this repo.
    # Point CTRL_TOWERS/CTRL_ORACLE at your own towers/oracle pair.
    TW = os.environ.get('CTRL_TOWERS', 'towers_P0b_K8.json')
    OR = os.environ.get('CTRL_ORACLE', 'master_towers_P0hi.json')
    if not (os.path.exists(TW) and os.path.exists(OR)):
        print('SKIP selftest: reference fixture towers not present '
              '(set CTRL_TOWERS/CTRL_ORACLE to your own pair)')
        return
    REF = {-3: 67.08, -2: 57.95, -1: 49.06, 0: 40.19, 1: 31.52, 2: 23.33, 3: 16.36}
    r = compare(TW, OR, 'P0hi', 0, 13, None)
    report(r)
    for k, v in REF.items():
        got = r['per_order'][k]['min']
        assert abs(got - v) <= 0.3, "order %+d: %.2f vs ref %.2f" % (k, got, v)
    assert 7.5 <= abs(r['slope']) <= 10.5, "slope %.3f out of [7.5,10.5]" % r['slope']
    assert r['verdict'].startswith("ALIAS-LIMITED"), r['verdict']
    assert r['usable_max_order'] <= 1, "usable max order %s > +1" % r['usable_max_order']
    assert r['per_order'][r['usable_max_order']]['min'] >= 30.0
    print("SELFTEST PASS")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--towers'); ap.add_argument('--oracle')
    ap.add_argument('--oracle-key', default='P0hi')
    ap.add_argument('--masters', default='0:13', help='half-open range lo:hi')
    ap.add_argument('--eps0', type=float, default=None, help='log10(eps0); default from towers file')
    ap.add_argument('--json', dest='json_out', default=None)
    ap.add_argument('--selftest', action='store_true')
    a = ap.parse_args()
    mp.dps = max(mp.dps, 80)
    if a.selftest:
        selftest(); return
    if not (a.towers and a.oracle):
        ap.error('--towers and --oracle are required (or use --selftest)')
    lo, hi = (int(x) for x in a.masters.split(':'))
    r = compare(a.towers, a.oracle, a.oracle_key, lo, hi, a.eps0)
    report(r)
    if a.json_out:
        json.dump(r, open(a.json_out, 'w'), indent=1)
        print("wrote", a.json_out)


if __name__ == '__main__':
    main()
