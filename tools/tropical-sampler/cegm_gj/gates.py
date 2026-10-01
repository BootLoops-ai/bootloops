#!/usr/bin/env python3
"""Two-run gate arithmetic: floored-integer digit agreement between run values.

Usage: gates.py A B [--bar D]   where A,B are run JSON files (field 'value')
or literal decimal strings. Prints rel diff and floored agreed digits;
rc=0 iff floored digits >= bar (default 30).
"""
import sys, json, argparse
import mpmath as mp

def load(x):
    if x.endswith('.json'):
        return json.load(open(x))['value']
    return x

ap = argparse.ArgumentParser()
ap.add_argument('a'); ap.add_argument('b')
ap.add_argument('--bar', type=float, default=30.0)
ap.add_argument('--relmax', type=float, default=None,
                help='optional strict rel-diff gate, e.g. 1e-30')
args = ap.parse_args()
mp.mp.dps = 120
va, vb = mp.mpf(load(args.a)), mp.mpf(load(args.b))
rel = abs(va - vb)/abs(va)
digits = float(-mp.log10(rel)) if rel > 0 else 120.0
fl = int(mp.floor(digits))
res = dict(rel=mp.nstr(rel, 4), digits=round(digits, 2), floored=fl,
           bar=args.bar, gate='PASS' if fl >= args.bar else 'FAIL')
if args.relmax is not None:
    res['relmax'] = args.relmax
    res['rel_gate'] = 'PASS' if rel < mp.mpf(args.relmax) else 'FAIL'
print(json.dumps(res))
ok = fl >= args.bar and (args.relmax is None or rel < mp.mpf(args.relmax))
sys.exit(0 if ok else 1)
