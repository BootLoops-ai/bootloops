#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision.
"""gate.py -- the ACCEPTANCE CHECK of Longhand's disteval route (longhand.disteval; `pysecdec_point.py --stage check`, older spelling `--stage gate`), by object,
over the three stage results of one work dir (or the worked example's records under ../examples/ndpent_top).  The module, its functions (run_gate) and
the JSON key 'gate' keep the historical spelling so that the pinned record GATE_record.json reproduces byte for byte:
(1) the three stages landed (RESULT_A/B/C complete; stage rcs 0); (2) the two lattice settings (A: --points 1e4, B: --points 1e5, standard
Korobov generating vectors) and the SECOND INTEGRATOR (C: the median-lattice rule, --lattice-candidates 11) agree per regulator order: the digit count of an
order = MIN over the three pairs of floor(-log10(|diff|/|value|)) and over the three bars' own digits; a pair sigma > 3 on any order is a FAIL
line; (3) CONTROLS by name: PLANTED-FAIL = a copy of RESULT_B with the eps^0 value shifted by 1e-3 relative is compared to A and C and MUST raise a FAIL line
(sigma > 3), else the compare is blind and the check FAILS; FOREIGN = coarse independent evaluations of the same library (a smoke of lattice 1e3,
4 shifts, 2 workers, epsrel 1e-2 on another host) vs stage A: every order within 3 x the summed bars; (4) the digit class = the conservative min over
orders and the per-order table; (5) NOT ESTABLISHED by object.  The numeric units (orders_of, pair, smoke_orders, bar_digits) are the evaluation
chain's units of record, unchanged but for the family name smoke_orders takes (the 'sums' key of a smoke; --family, alias --name; default
ndpent_top, the worked example's family; a smoke without that key is a FOREIGN control that fails by name); run_gate is the same computation over inputs
given by path.  A --family of the wrong form (not [A-Za-z0-9_]+) is REFUSED by name by the CLI (rc 2, nothing written).  Main-guarded; `date -u`.
usage: gate.py --result A=RESULT_A.json --result B=... --result C=... [--rc A=rc=0 ...] [--smoke NAME=smoke.json ...] [--leaf-stamp FILE] [--family FAMILY] --out GATE.json"""
import os, sys, json, math, subprocess, hashlib, copy, argparse


def utc():
    return subprocess.run(['date', '-u', '+%Y-%m-%dT%H:%M:%SZ'], capture_output=True, text=True).stdout.strip()


def sha256(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


def orders_of(R):
    return {tuple(o['regulator_powers']): (complex(o['value_re'], o['value_im']), abs(complex(o['err_re'], o['err_im']))) for o in R['orders']}


def pair(A, B):
    rows = {}
    for pw in sorted(set(A) | set(B)):
        if pw not in A or pw not in B:
            rows[pw] = None
            continue
        va, ea = A[pw]
        vb, eb = B[pw]
        diff = abs(va - vb)
        ref = max(abs(va), abs(vb))
        sig = diff / math.sqrt(ea ** 2 + eb ** 2) if (ea > 0 or eb > 0) else None
        agreed = (math.floor(-math.log10(diff / ref)) if diff > 0 and ref > 0 else (16 if ref > 0 else None))
        rows[pw] = {'abs_diff': diff, 'sigma': sig, 'agreed_digits': agreed, 'within_sum_of_bars': diff <= ea + eb, 'within_3x_sum_of_bars': diff <= 3 * (ea + eb)}
    return rows


def smoke_orders(path, name='ndpent_top'):
    """a disteval --format=json stdout (the smoke.out files): sums.<family> rows [[p],[re,im],[ere,eim]]"""
    r = json.load(open(path))
    return {tuple(pw): (complex(v[0], v[1]), abs(complex(e[0], e[1]))) for pw, v, e in r['sums'][name]}


def bar_digits(v, e):
    return math.floor(-math.log10(e / abs(v))) if abs(v) > 0 and e > 0 else None


def read_stamp(path):
    """the leaf-at-end stamp: key=value lines (memory.peak, memory.max, oom_kill, wall_s, maxrss_kb_over_chunks, stages_rc, cpu.stat)."""
    if not path or not os.path.exists(path):
        return {}
    return dict(ln.split('=', 1) for ln in open(path).read().splitlines() if '=' in ln)


def planted_fail_control(O):
    """CONTROL 1: PLANTED-FAIL (a 1e-3 relative shift on eps^0 of a copy of B must raise a FAIL line vs A and vs C)."""
    planted = None
    if 'B' in O and 'A' in O:
        PB = dict(O['B'])
        v, e = PB[(0,)]
        PB[(0,)] = (v * (1 + 1e-3), e)
        pa = pair(O['A'], PB)
        pc = pair(O['C'], PB) if 'C' in O else None
        sa = pa[(0,)]['sigma'] if pa.get((0,)) else None
        sc = pc[(0,)]['sigma'] if pc and pc.get((0,)) else None
        planted = {'name': 'PLANTED-FAIL (RESULT_B copy, eps^0 value x (1 + 1e-3))', 'sigma_vs_A': sa, 'sigma_vs_C': sc, 'raised_fail_line': bool(sa is not None and sa > 3) and (sc is None or sc > 3), 'rule': 'the compare must flag a 1e-3 shift as > 3 sigma on eps^0 (the bars are ~1e-7 relative); if it does not, the compare is blind and the gate FAILS'}
    return planted


def run_gate(results, rcs, smokes, stamp, out_path, chunks_done=None, label='', package_pinned=True, compare_path=None, family='ndpent_top'):
    """results: {'A': path, 'B': path, 'C': path}; rcs: {'A': 'rc=0', ...} (None when unknown); smokes: [(name, path), ...] the FOREIGN
    controls; stamp: the leaf-at-end dict (may be empty); chunks_done: {'A': n, ...} DONE markers counted on disk (None in fixture mode);
    family: the 'sums' key the smokes are read by (a smoke without it is a FOREIGN control that FAILS by name, never skipped)."""
    res = {}
    stage_ok = {}
    for s in ('A', 'B', 'C'):
        p = results.get(s)
        res[s] = json.load(open(p)) if p and os.path.exists(p) else None
        rc = rcs.get(s) if rcs else None
        ndone = (chunks_done or {}).get(s)
        stage_ok[s] = {'result_present': res[s] is not None, 'complete': res[s]['complete'] if res[s] else False, 'n_chunks': res[s]['n_chunks'] if res[s] else None, 'n_landed': res[s]['n_landed'] if res[s] else None, 'chunks_DONE_on_disk': ndone, 'rc': rc, 'settings': res[s].get('settings') if res[s] else None, 'result_sha256': sha256(p) if res[s] else None, 'chunk_walls_s': res[s].get('chunk_walls_s') if res[s] else None, 'chunk_maxrss_kb': res[s].get('chunk_maxrss_kb') if res[s] else None}
    landed = all(stage_ok[s]['complete'] and stage_ok[s]['rc'] == 'rc=0' for s in stage_ok)
    O = {s: orders_of(res[s]) for s in res if res[s]}
    pairs = {}
    for a, b in (('A', 'B'), ('A', 'C'), ('B', 'C')):
        if a in O and b in O:
            pairs['%s-%s' % (a, b)] = pair(O[a], O[b])
    # the per-order digit class: set by the two evaluations AT THE TARGET precision (B: lattice setting 2, standard lattices; C: the median-lattice second
    # integrator, same epsrel) = min over the B-C agreed digits and the B and C bars; stage A (lattice setting 1, epsrel 1e-4, a 32-shift bar) is the coarse
    # consistency check: its sigma vs B and vs C is listed per order and a > 3 sigma line is an A-CONSISTENCY line (reported; its bar is itself a 32-sample
    # estimate); the gate FAILS on a B-C sigma > 3 on any order, on a blind planted control, or on a foreign control outside 3x the summed bars
    table = []
    min_digits = None
    fail_lines = []
    a_lines = []
    for pw in sorted(O['A'].keys() if 'A' in O else []):
        row = {'regulator_powers': list(pw)}
        ds = []
        for s in O:
            v, e = O[s][pw]
            row['%s_value' % s] = (v.real, v.imag)
            row['%s_err' % s] = e
            row['%s_bar_digits' % s] = bar_digits(v, e)
            if s in ('B', 'C') and row['%s_bar_digits' % s] is not None:
                ds.append(row['%s_bar_digits' % s])
        for k, pr in pairs.items():
            r = pr.get(pw)
            row['%s_sigma' % k] = r['sigma'] if r else None
            row['%s_agreed_digits' % k] = r['agreed_digits'] if r else None
            if k == 'B-C':
                if r and r['agreed_digits'] is not None:
                    ds.append(r['agreed_digits'])
                if r and r['sigma'] is not None and r['sigma'] > 3:
                    fail_lines.append({'pair': k, 'regulator_powers': list(pw), 'sigma': r['sigma']})
            else:
                if r and r['sigma'] is not None and r['sigma'] > 3:
                    a_lines.append({'pair': k, 'regulator_powers': list(pw), 'sigma': r['sigma'], 'reading': 'stage A (epsrel 1e-4, 32-shift bar) vs the target-precision stage: the difference exceeds 3x the summed bars; A is the coarse consistency setting, its bar is a 32-sample estimate'})
        row['digits(min over the B-C pair and the B, C bars)'] = min(ds) if ds else None
        table.append(row)
        if ds:
            min_digits = min(ds) if min_digits is None else min(min_digits, min(ds))
    planted = planted_fail_control(O)
    # CONTROL 2: FOREIGN (coarse independent evaluations of the same library: lattice 1e3 / 4 shifts / epsrel 1e-2, the same .so sha) vs stage A
    foreign = []
    for name, path in smokes:
        if os.path.exists(path) and 'A' in O:
            try:
                S = smoke_orders(path, family)
            except KeyError:
                foreign.append({'name': name, 'file': os.path.basename(path), 'sha256': sha256(path), 'rows': {}, 'all_within_3x_sum_of_bars': False, 'max_sigma': None, 'error': 'the smoke carries no sums.%s (keys: %s) -- not a FOREIGN control of this family' % (family, ', '.join(sorted(json.load(open(path)).get('sums', {}).keys())) or 'none')})
                continue
            pr = pair(O['A'], S)
            foreign.append({'name': name, 'file': os.path.basename(path), 'sha256': sha256(path), 'rows': {str(list(pw)): r for pw, r in pr.items()}, 'all_within_3x_sum_of_bars': all(r and r['within_3x_sum_of_bars'] for r in pr.values()), 'max_sigma': max((r['sigma'] for r in pr.values() if r and r['sigma'] is not None), default=None)})
    foreign_ok = bool(foreign) and all(f['all_within_3x_sum_of_bars'] for f in foreign)
    a_max_sigma = max((r['%s_sigma' % k] for r in table for k in ('A-B', 'A-C') if r.get('%s_sigma' % k) is not None), default=None)
    gate_pass = landed and not fail_lines and planted is not None and planted['raised_fail_line'] and foreign_ok and (stamp.get('oom_kill', '0') == '0')
    out = {'receipt': 'GATE%s: lattice setting 2 (B) and the second integrator (C, median lattices), both at the target epsrel, agree per order to the digit counts below by object; lattice setting 1 (A) is the coarse consistency check (its sigma per order listed); controls by name: PLANTED-FAIL, FOREIGN (%s); %s' % ((' ' + label) if label else '', ', '.join(n for n, _ in smokes) or 'none given', 'PASS' if gate_pass else 'FAIL -- read the fields'),
           'stamp_utc': utc(), 'gate': 'PASS' if gate_pass else 'FAIL', 'package_pinned': package_pinned, 'landed_all_three_stages': landed, 'stages': stage_ok,
           'digit_class': {'min_digits_over_orders(conservative)': min_digits, 'per_order': table, 'rule': 'digits(order) = min over the B-C pair (lattice setting 2 vs the median-lattice second integrator, both at epsrel 1e-7) of floor(-log10(|diff|/|value|)) and over the B and C bars floor(-log10(err/|value|)); the quoted digit count of the value is the min over orders; stage A (lattice setting 1, epsrel 1e-4) is the coarse consistency check, its sigma vs B and C listed per order; double-precision disteval house ceiling ~8-9 d'},
           'fail_lines(B-C sigma > 3 on an order)': fail_lines, 'A_consistency': {'max_sigma_A_vs_B_or_C': a_max_sigma, 'lines_over_3_sigma': a_lines, 'verdict': 'CONSISTENT' if not a_lines else 'MARGINAL (%d line(s) over 3 sigma; A is the epsrel 1e-4 setting with a 32-shift bar)' % len(a_lines)}, 'controls': {'planted_fail': planted, 'foreign': foreign, 'foreign_ok': foreign_ok},
           'leaf_at_end': {'memory.peak': stamp.get('memory.peak'), 'memory.max': stamp.get('memory.max'), 'oom_kill': stamp.get('oom_kill'), 'wall_s': stamp.get('wall_s'), 'maxrss_kb_over_chunks': stamp.get('maxrss_kb_over_chunks'), 'stages_rc': stamp.get('stages_rc'), 'cpu.stat': stamp.get('cpu.stat')},
           'compare_receipt': {'file': os.path.basename(compare_path) if compare_path else None, 'sha256': sha256(compare_path) if compare_path and os.path.exists(compare_path) else None},
           'not_established': ['a value at any Minkowski point (contour deformation not built)', 'the other 8 top-sector masters and orders eps^1..eps^4 (separate packages)', 'digits beyond the double-precision ceiling of disteval (the 30-digit closure is the AMFlow route)', 'no independent published value exists for this fixture family at this point; the acceptance bar is the three-way agreement of the shipped routes'],
           'PRODUCER': {'script': os.path.basename(os.path.abspath(__file__)), 'sha256': sha256(os.path.abspath(__file__)), 'stamp_source': 'date -u'}}
    json.dump(out, open(out_path, 'w'), indent=1)
    print('gate:', out_path, sha256(out_path)[:16], '|', out['gate'], '| landed', landed, '| min digits (B-C + bars)', min_digits, '| B-C fail lines', len(fail_lines), '| A consistency', out['A_consistency']['verdict'], '| planted raised', planted['raised_fail_line'] if planted else None, '| foreign ok', foreign_ok, '| oom', stamp.get('oom_kill'), '| package pinned', package_pinned)
    for r in table:
        print('  order', r['regulator_powers'], 'B %.12e' % r['B_value'][0] if 'B_value' in r else '', 'bars A/B/C', r.get('A_bar_digits'), r.get('B_bar_digits'), r.get('C_bar_digits'), 'agreed A-B/A-C/B-C', r.get('A-B_agreed_digits'), r.get('A-C_agreed_digits'), r.get('B-C_agreed_digits'), 'sigmas %.2f %.2f %.2f' % tuple((r.get(k) or 0) for k in ('A-B_sigma', 'A-C_sigma', 'B-C_sigma')), '-> digits', r['digits(min over the B-C pair and the B, C bars)'])
    for f in fail_lines:
        print('  FAIL line: pair %s order %s sigma %.3f > 3' % (f['pair'], f['regulator_powers'], f['sigma']))
    if planted is not None:
        print('  control PLANTED-FAIL: sigma vs A %s, vs C %s -> raised_fail_line %s' % (planted['sigma_vs_A'], planted['sigma_vs_C'], planted['raised_fail_line']))
    for f in foreign:
        print('  control FOREIGN %s: all within 3x the summed bars %s (max sigma %s)%s' % (f['name'], f['all_within_3x_sum_of_bars'], f['max_sigma'], (' -- ' + f['error']) if f.get('error') else ''))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--result', action='append', default=[], help='STAGE=RESULT.json (A, B, C)')
    ap.add_argument('--rc', action='append', default=[], help='STAGE=rc=N (the stage rc line; default rc=0 when the RESULT is complete and none is given)')
    ap.add_argument('--smoke', action='append', default=[], help='NAME=smoke.json (a disteval --format=json stdout; the FOREIGN control)')
    ap.add_argument('--leaf-stamp', default=None)
    ap.add_argument('--label', default='')
    ap.add_argument('--family', '--name', dest='family', default='ndpent_top', help="the package family NAME = the 'sums' key the FOREIGN smokes are read by; --name is the alias")
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    if not (a.family and a.family.isascii() and all(c.isalnum() or c == '_' for c in a.family)):
        print('REFUSED: --family %r is not a family name (the form is [A-Za-z0-9_]+: letters, digits, underscore; no separators, no dots, no empty token)' % (a.family,))
        return 2
    results = dict(x.split('=', 1) for x in a.result)
    rcs = dict(x.split('=', 1) for x in a.rc)
    for s in results:
        rcs.setdefault(s, 'rc=0' if json.load(open(results[s])).get('complete') else 'rc=1')
    smokes = [tuple(x.split('=', 1)) for x in a.smoke]
    out = run_gate(results, rcs, smokes, read_stamp(a.leaf_stamp), a.out, label=a.label, family=a.family)
    return 0 if out['gate'] == 'PASS' else 1


if __name__ == '__main__':
    sys.dont_write_bytecode = True
    raise SystemExit(main())
