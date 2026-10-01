#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision.
"""assemble_compare.py -- (1) ASSEMBLE the per-chunk disteval results of one stage into RESULT_<stage>.json (values add per regulator
order; variances add, i.e. the bars add in quadrature; partial when a chunk has no DONE marker: rc 3); (2) COMPARE the stages' results per order
(pairwise |diff|, sigma = |diff| / sqrt(err_a^2 + err_b^2), agreed digits = floor(-log10(|diff| / |value|))) into a COMPARE json.
NO digit is quoted here: the file states the rule (a digit count is the MIN over the pairs of agreed digits AND over each bar's own digits, per
order; gate.py reads it).  Main-guarded; stamps from `date -u`.  The units below are the evaluation chain's units of record, unchanged but for
the family name the assembly takes (--family, alias --name; default ndpent_top, the worked example's family = the 'sums' key of a disteval result).
The assemble CLI REFUSES by name (rc 2, nothing written) a --family of the wrong form (not [A-Za-z0-9_]+), a chunks dir without its
CHUNKS_<stage>.json, and a landed chunk result that carries no sums.<family> (chunk_results_problems lists what is wrong).
usage: assemble_compare.py assemble --chunks-dir O --stage NAME --out RESULT.json [--family FAMILY]
       assemble_compare.py compare --out COMPARE.json A=RESULT_A.json B=RESULT_B.json [C=RESULT_C.json ...]"""
import os, sys, json, argparse, hashlib, subprocess, math, glob


def utc():
    return subprocess.run(['date', '-u', '+%Y-%m-%dT%H:%M:%SZ'], capture_output=True, text=True).stdout.strip()


def sha256(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


def read_result(p, name='ndpent_top'):
    """disteval --format=json stdout: {'regulators': [...], 'sums': {'<family>': [[[p], [re, im], [ere, eim]], ...]}, 'integrals': {...}}"""
    txt = open(p).read()
    r = json.loads(txt)
    terms = r['sums'][name]
    out = {}
    for pw, (vre, vim), (ere, eim) in terms:
        out[tuple(pw)] = (complex(vre, vim), complex(ere, eim))
    return out


def assemble(a):
    man = json.load(open(os.path.join(a.chunks_dir, 'CHUNKS_%s.json' % a.stage)))
    val = {}
    var = {}
    landed = []
    missing = []
    walls = []
    rss = []
    for c in man['chunks']:
        d = c['dir']
        done = os.path.exists(os.path.join(d, 'DONE'))
        rp = os.path.join(d, 'result.json')
        if not (done and os.path.exists(rp)):
            missing.append(c['id'])
            continue
        r = read_result(rp, a.name)
        for pw, (v, e) in r.items():
            val[pw] = val.get(pw, 0) + v
            var[pw] = var.get(pw, 0) + complex(e.real ** 2, e.imag ** 2)
        landed.append(c['id'])
        tv = os.path.join(d, 'time.txt')
        if os.path.exists(tv):
            t = open(tv).read()
            for ln in t.splitlines():
                if 'Elapsed (wall clock)' in ln:
                    hms = ln.split(':', 1)[1].strip()
                    hms = ln.rsplit(' ', 1)[1]
                    parts = [float(x) for x in hms.split(':')]
                    walls.append(sum(x * 60 ** i for i, x in enumerate(reversed(parts))))
                if 'Maximum resident set size' in ln:
                    rss.append(int(ln.rsplit(' ', 1)[1]))
    orders = []
    for pw in sorted(val.keys()):
        v = val[pw]
        e = complex(math.sqrt(var[pw].real), math.sqrt(var[pw].imag))
        absv = abs(v)
        abse = abs(e)
        orders.append({'regulator_powers': list(pw), 'value_re': v.real, 'value_im': v.imag, 'err_re': e.real, 'err_im': e.imag,
                       'rel_err': (abse / absv) if absv > 0 else None, 'bar_digits': (math.floor(-math.log10(abse / absv)) if absv > 0 and abse > 0 else None)})
    out = {'receipt': 'RESULT stage %s: the sum of %d landed chunks of %d (values add per order; bars add in quadrature); NOT a quoted value -- the digits come from the stage-to-stage COMPARE' % (a.stage, len(landed), man['n_chunks']),
           'stamp_utc': utc(), 'stage': a.stage, 'chunks_manifest': os.path.join(a.chunks_dir, 'CHUNKS_%s.json' % a.stage), 'chunks_manifest_sha256': sha256(os.path.join(a.chunks_dir, 'CHUNKS_%s.json' % a.stage)),
           'name': a.name, 'n_chunks': man['n_chunks'], 'n_landed': len(landed), 'landed': landed, 'missing': missing, 'complete': not missing,
           'kernels_cover_package': man.get('kernels_cover_package'), 'n_kernels': man.get('n_kernels_in_chunks'),
           'orders': orders, 'chunk_walls_s': {'n': len(walls), 'sum': sum(walls) if walls else None, 'max': max(walls) if walls else None},
           'chunk_maxrss_kb': {'n': len(rss), 'max': max(rss) if rss else None}, 'settings': (json.load(open(a.settings)) if a.settings and os.path.exists(a.settings) else None)}
    json.dump(out, open(a.out, 'w'), indent=1)
    print('RESULT %s: landed %d/%d complete=%s orders=%s -> %s' % (a.stage, len(landed), man['n_chunks'], out['complete'], [(o['regulator_powers'], o['bar_digits']) for o in orders], a.out))
    return 0 if not missing else 3


def compare(a):
    R = {}
    for spec in a.results:
        name, path = spec.split('=', 1)
        j = json.load(open(path))
        R[name] = {'path': path, 'sha256': sha256(path), 'complete': j['complete'], 'orders': {tuple(o['regulator_powers']): o for o in j['orders']}, 'settings': j.get('settings')}
    names = list(R.keys())
    pairs = []
    orders = sorted(set().union(*[set(R[n]['orders'].keys()) for n in names]))
    min_agreed = None
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            A, B = R[names[i]], R[names[j]]
            rows = []
            for pw in orders:
                oa = A['orders'].get(pw)
                ob = B['orders'].get(pw)
                if not oa or not ob:
                    rows.append({'regulator_powers': list(pw), 'status': 'missing in one stage'})
                    continue
                va = complex(oa['value_re'], oa['value_im'])
                vb = complex(ob['value_re'], ob['value_im'])
                ea = abs(complex(oa['err_re'], oa['err_im']))
                eb = abs(complex(ob['err_re'], ob['err_im']))
                diff = abs(va - vb)
                sig = diff / math.sqrt(ea ** 2 + eb ** 2) if (ea > 0 or eb > 0) else None
                ref = max(abs(va), abs(vb))
                agreed = (math.floor(-math.log10(diff / ref)) if diff > 0 and ref > 0 else (None if ref == 0 else 16))
                bard = min([d for d in (oa.get('bar_digits'), ob.get('bar_digits')) if d is not None] or [None])
                rows.append({'regulator_powers': list(pw), 'a': (va.real, va.imag), 'b': (vb.real, vb.imag), 'err_a': ea, 'err_b': eb, 'abs_diff': diff, 'sigma': sig, 'agreed_digits': agreed, 'bar_digits_min': bard,
                             'within_sum_of_bars': (diff <= ea + eb)})
                if agreed is not None:
                    dd = min([x for x in (agreed, bard) if x is not None])
                    min_agreed = dd if min_agreed is None else min(min_agreed, dd)
            pairs.append({'pair': [names[i], names[j]], 'rows': rows})
    out = {'receipt': 'COMPARE: stage-to-stage agreement per regulator order (two lattice settings + the median-lattice second rule); the quoted digit count of a value is MIN over the pairs of agreed_digits and over the bars (bar_digits) for that order -- read, never assumed; a sigma > 3 on any pair is a FAIL line for that order',
           'stamp_utc': utc(), 'stages': {n: {k: v for k, v in R[n].items() if k != 'orders'} for n in names}, 'all_complete': all(R[n]['complete'] for n in names),
           'pairs': pairs, 'min_digits_over_all_orders_and_pairs(the conservative digit count)': min_agreed,
           'fail_lines': [{'pair': p['pair'], 'regulator_powers': r['regulator_powers'], 'sigma': r['sigma']} for p in pairs for r in p['rows'] if r.get('sigma') is not None and r['sigma'] > 3]}
    json.dump(out, open(a.out, 'w'), indent=1)
    print('COMPARE: stages %s all_complete=%s min_digits=%s fail_lines=%d -> %s' % (names, out['all_complete'], min_agreed, len(out['fail_lines']), a.out))
    return 0


def chunk_results_problems(chunks_dir, stage, name):
    """the by-name checks of the assemble CLI, in order: the family NAME form ([A-Za-z0-9_]+), the manifest CHUNKS_<stage>.json (present; JSON;
    an object with a 'chunks' list), and every landed chunk (DONE + result.json): result.json is JSON, an object whose 'sums' object carries
    the key NAME.  Returns the list of problems found (empty = the assembly can run); reads only, raises nothing."""
    if not (isinstance(name, str) and name and name.isascii() and all(c.isalnum() or c == '_' for c in name)):
        return ['--family %r is not a family name (the form is [A-Za-z0-9_]+: letters, digits, underscore; no separators, no dots, no empty token)' % (name,)]
    mp = os.path.join(chunks_dir, 'CHUNKS_%s.json' % stage)
    if not os.path.isfile(mp):
        return ['%s missing in %s' % (os.path.basename(mp), chunks_dir)]
    try:
        man = json.load(open(mp))
    except (ValueError, UnicodeDecodeError) as e:
        return ['%s is not JSON (%s)' % (os.path.basename(mp), e.__class__.__name__)]
    if not isinstance(man, dict) or not isinstance(man.get('chunks'), list):
        return ['%s carries no chunks list' % os.path.basename(mp)]
    out = []
    for c in man['chunks']:
        d = c.get('dir') if isinstance(c, dict) else None
        rp = os.path.join(d, 'result.json') if d else None
        if not (d and os.path.exists(os.path.join(d, 'DONE')) and os.path.exists(rp)):
            continue
        try:
            r = json.load(open(rp))
        except (ValueError, UnicodeDecodeError) as e:
            out.append('%s is not JSON (%s)' % (rp, e.__class__.__name__))
            continue
        sums = r.get('sums') if isinstance(r, dict) else None
        if not isinstance(sums, dict) or name not in sums:
            out.append('%s carries no sums.%s (keys: %s)' % (rp, name, ', '.join(sorted(sums.keys())) if isinstance(sums, dict) and sums else 'none'))
    return out


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    p1 = sub.add_parser('assemble')
    p1.add_argument('--chunks-dir', required=True)
    p1.add_argument('--stage', required=True)
    p1.add_argument('--out', required=True)
    p1.add_argument('--settings', default=None)
    p1.add_argument('--family', '--name', dest='name', default='ndpent_top', help="the package family NAME = the 'sums' key of the chunk results; --name is the alias")
    p2 = sub.add_parser('compare')
    p2.add_argument('--out', required=True)
    p2.add_argument('results', nargs='+')
    a = ap.parse_args()
    if a.cmd == 'assemble':
        probs = chunk_results_problems(a.chunks_dir, a.stage, a.name)
        if probs:
            print('REFUSED: assemble --family %s over %s: %s' % (a.name, a.chunks_dir, '; '.join(probs)))
            return 2
        return assemble(a)
    return compare(a)


if __name__ == '__main__':
    sys.dont_write_bytecode = True
    raise SystemExit(main())
