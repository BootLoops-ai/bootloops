#!/usr/bin/env python3
"""gmtel route B (DATA route): exact moment series of the Watson/HLY K3
pencil at rational rates -> minimal recurrence via annihilator.pf_from_series
(mod-p pin at TWO primes, exact Q-nullspace, held-out terms), then the
KNOWN-comparator gate vs the reference L5 (L5_theta.txt, cited by sha).

All verdict fields SCRIPT-EMITTED.  Usage:
  python3 famgen_pf_routeB.py --abc 1 4 16 --mmax 140 --out CONTROL_PF_ROUTEB.json
"""
import argparse, hashlib, json, os, subprocess, sys, time
from fractions import Fraction as F
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
# Dir holding the annihilator engine module (annihilator.py); the packaged
# default is the sibling tools/annihilator/ directory in this tree.
TOOLS = os.environ.get('FAMGEN_TOOLS',
                       os.path.join(HERE, '..', '..', 'annihilator'))
sys.path.insert(0, TOOLS)

from famgen_common import (moments, load_L5, L5_moment_rec, rec_to_op_x,
                           op_primitive, op_equal, L5_PATH)
from annihilator import pf_from_series, find_recurrence_fast  # house engine

P1 = (1 << 31) - 1          # 2147483647 (annihilator default P31)
P2 = 2147483629             # second pin prime


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--abc', nargs=3, type=int, required=True,
                    help='(A,B,C) = squared rates, integers')
    ap.add_argument('--mmax', type=int, default=140)
    ap.add_argument('--nverify', type=int, default=15)
    ap.add_argument('--out', default='CONTROL_PF_ROUTEB.json')
    args = ap.parse_args()
    Av, Bv, Cv = args.abc
    t0 = time.time()

    rep = {'leg': 'gmtel', 'route': 'B (data: exact moments -> annihilator fit)',
           'ABC': [Av, Bv, Cv], 'mmax': args.mmax,
           'stamp_utc': subprocess.run(['date', '-u', '+%Y-%m-%dT%H:%M:%SZ'],
                                       capture_output=True, text=True).stdout.strip(),
           'engine': 'tools/annihilator pf_from_series'}

    # ---- exact moments
    a = moments(Av, Bv, Cv, args.mmax)
    rep['moments_first5'] = [str(x) for x in a[:5]]
    rep['t_moments_s'] = round(time.time() - t0, 2)

    # ---- mod-p pin at two primes (independent scans), then exact fit
    hit1 = find_recurrence_fast([int(x % P1) if x.denominator == 1 else
                                 int((x.numerator * pow(x.denominator, -1, P1)) % P1)
                                 for x in a], P1, rmax=10, dmax=8, overdet=5,
                                prefer='ode_order')
    hit2 = find_recurrence_fast([int(x % P2) if x.denominator == 1 else
                                 int((x.numerator * pow(x.denominator, -1, P2)) % P2)
                                 for x in a], P2, rmax=10, dmax=8, overdet=5,
                                prefer='ode_order')
    rep['modp_pin'] = {'p1': P1, 'rs_p1': hit1[:2] if hit1 else None,
                       'p2': P2, 'rs_p2': hit2[:2] if hit2 else None}
    if hit1 is None or hit2 is None or hit1[:2] != hit2[:2]:
        rep['verdict'] = 'FAIL: two-prime (r,s) pin disagreement or no hit'
        json.dump(rep, open(args.out, 'w'), indent=1)
        print(json.dumps(rep, indent=1)); return 1
    r, s = hit1[:2]

    res = pf_from_series(a, rmax=r + 1, smax=s + 1, nverify=args.nverify, mode='auto')
    if res is None or res[2] is None:
        rep['verdict'] = 'FAIL: exact Q-reconstruction did not verify'
        json.dump(rep, open(args.out, 'w'), indent=1)
        print(json.dumps(rep, indent=1)); return 1
    r, s, coeffs, nheld = res
    rep['fit'] = {'r': r, 's': s, 'n_heldout_verified': nheld}
    rep['t_fit_s'] = round(time.time() - t0, 2)

    # recurrence as sympy polys c_j(m)
    m = sp.symbols('m')
    rec = {}
    for (j, k), val in coeffs.items():
        rec[j] = sp.expand(rec.get(j, 0) + sp.Rational(val.numerator, val.denominator) * m ** k)
    rec = {j: sp.Poly(c, m) for j, c in rec.items() if sp.expand(c) != 0}

    # full-series annihilation re-check (independent of the fit's own internals)
    ok_all = True
    for n in range(0, args.mmax - r):
        v = F(0)
        for j, cp in rec.items():
            v += F(str(cp.eval(n))) * a[n + j]
        if v != 0:
            ok_all = False
            break
    rep['recurrence_annihilates_all_terms'] = bool(ok_all)

    # ---- KNOWN comparator: the reference L5 at (A,B,C)
    sha = hashlib.sha256(open(L5_PATH, 'rb').read()).hexdigest()[:16]
    rep['known_comparator'] = {'file': L5_PATH, 'sha256_16': sha}
    P = load_L5()
    l5rec = L5_moment_rec(P, Av, Bv, Cv)
    # does the L5 recurrence annihilate the moments? (consistency vs the reference operator)
    ok_l5 = True
    for n in range(0, args.mmax - 6):
        v = sp.Integer(0)
        for j, cp in l5rec.items():
            v += cp.eval(n) * sp.Rational(a[n + j].numerator, a[n + j].denominator)
        if sp.simplify(v) != 0:
            ok_l5 = False
            break
    rep['L5_annihilates_moments'] = bool(ok_l5)

    # operator-level comparison in x-space (canonical primitive D-forms)
    x = sp.symbols('x')
    op_fit = op_primitive(rec_to_op_x(rec, x), x)
    op_l5 = op_primitive(rec_to_op_x(l5rec, x), x)
    same = op_equal(op_fit, op_l5, x)
    rep['fit_equals_L5_exactly'] = bool(same)
    if not same:
        # minimal-order sub-case (symmetric locus): report right-factor status honestly
        rep['note'] = ('fitted minimal operator differs from L5 specialization; '
                       'expected only on degeneration loci (e.g. A=B=C) where the '
                       'minimal order drops — both-annihilate facts above still hold')
    rep['op_fit_order'] = max(op_fit)
    rep['op_L5_order'] = max(op_l5)
    rep['op_fit_xdeg'] = max(sp.Poly(c, x).degree() for c in op_fit.values())
    rep['op_fit'] = {str(j): str(c) for j, c in sorted(op_fit.items())}

    verdict = 'PASS' if (ok_all and ok_l5 and (same or max(op_fit) < max(op_l5))) else 'FAIL'
    if same:
        verdict_detail = 'fitted minimal recurrence == L5 specialization EXACTLY (primitive forms)'
    elif max(op_fit) < max(op_l5) and ok_l5:
        verdict_detail = ('fitted minimal order < L5 order (degeneration locus); both '
                          'annihilate all computed terms')
    else:
        verdict_detail = 'mismatch'
    rep['verdict'] = verdict
    rep['verdict_detail'] = verdict_detail
    rep['t_total_s'] = round(time.time() - t0, 2)
    json.dump(rep, open(args.out, 'w'), indent=1)
    print(json.dumps({k: rep[k] for k in ('ABC', 'modp_pin', 'fit',
                                          'recurrence_annihilates_all_terms',
                                          'L5_annihilates_moments',
                                          'fit_equals_L5_exactly', 'verdict',
                                          'verdict_detail', 't_total_s')}, indent=1))
    return 0 if verdict == 'PASS' else 1


if __name__ == '__main__':
    sys.exit(main())
