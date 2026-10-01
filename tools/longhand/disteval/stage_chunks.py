#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision.
"""stage_chunks.py -- per-sector CHECKPOINT chunks for one disteval stage of Longhand's disteval route (longhand.disteval).
disteval evaluates one integral JSON at a time and writes nothing per sector (read: pySecDec/disteval.py 1.6.6, prepare_eval / do_eval /
the result dict) -- so the checkpoint unit is a CHUNK: a copy of the disteval data dir whose <family>_integral.json lists only the kernels of a
group of sectors (the kernels are named sector_<n>_order_<k>; the 'orders' lists are filtered the same way; the prefactor, the coefficient and
the sum file are unchanged, so the chunk results ADD per regulator order and their variances add).  The worker dlopens ./<family>.so from the
data dir it is started in (strings of pysecdec_cpuworker: './%s.so'), so every chunk dir carries hard links to <family>_integral.so and
builtin.so and a symlink to coefficients/.  The family NAME (--family, alias --name; default ndpent_top, the worked example's family) derives the file
names: NAME.json the sum file, NAME_integral.json / NAME_integral.so the integral; load(D, name) and write_chunk(..., name) take it.  A NAME of
the wrong form (not [A-Za-z0-9_]+), a sum file absent or not JSON or whose 'integrals' is not [NAME_integral], an integral file absent or not
JSON or whose 'name' is not NAME_integral: the CLI REFUSES by name, rc 2 (family_problems lists what is wrong; the families present are listed).
usage: stage_chunks.py --disteval-dir D --out O (--chunk-sectors N | --sectors 5,17,...) [--stage NAME] [--family FAMILY]
       stage_chunks.py --disteval-dir D --choose-typical [--family FAMILY]     (prints the typical sector: median index among the sectors with the modal kernel count)
Main-guarded; stamps from `date -u`; nothing written outside --out.  The units below are the evaluation chain's units of record, unchanged
but for the family name they take (family_problems is the by-name check the CLIs run before load)."""
import os, sys, json, argparse, hashlib, subprocess, collections, shutil


def utc():
    return subprocess.run(['date', '-u', '+%Y-%m-%dT%H:%M:%SZ'], capture_output=True, text=True).stdout.strip()


def sha256(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def sector_of(kernel):
    # kernel names: sector_<n>_order_<k>
    assert kernel.startswith('sector_'), kernel
    return int(kernel.split('_')[1])


def family_problems(D, name):
    """the by-name checks of a family in a disteval dir, in order: the NAME form ([A-Za-z0-9_]+: letters, digits, underscore; no separators, no
    dots, no empty token), the sum file NAME.json (present; JSON; an object of 'type' sum with 'integrals' == [NAME_integral]), the integral file
    NAME_integral.json (present; JSON; an object of 'type' integral with 'name' == NAME_integral).  Returns the list of problems found, the first
    one '<NAME>.json missing' when the sum file is absent (empty = the family names this package); reads only, raises nothing."""
    if not (isinstance(name, str) and name and name.isascii() and all(c.isalnum() or c == '_' for c in name)):
        return ['--family %r is not a family name (the form is [A-Za-z0-9_]+: letters, digits, underscore; no separators, no dots, no empty token)' % (name,)]
    out = []
    for rel, kind, key, want in ((name + '.json', 'sum', 'integrals', [name + '_integral']), (name + '_integral.json', 'integral', 'name', name + '_integral')):
        p = os.path.join(D, rel)
        if not os.path.isfile(p):
            out.append('%s missing' % rel)
            continue
        try:
            J = json.load(open(p))
        except (ValueError, UnicodeDecodeError) as e:
            out.append('%s is not JSON (%s)' % (rel, e.__class__.__name__))
            continue
        if not isinstance(J, dict) or J.get('type') != kind:
            out.append('%s is not a disteval %s file (type %r)' % (rel, kind, J.get('type') if isinstance(J, dict) else type(J).__name__))
            continue
        if J.get(key) != want:
            out.append('%s %r is %r, not %r' % (rel, key, J.get(key), want))
    return out


def load(D, name='ndpent_top'):
    S = json.load(open(os.path.join(D, name + '.json')))
    assert S['type'] == 'sum' and S['integrals'] == [name + '_integral'], 'not the %s sum file' % name
    I = json.load(open(os.path.join(D, name + '_integral.json')))
    assert I['type'] == 'integral' and I['name'] == name + '_integral'
    return S, I


def sectors_in_order(I):
    seen = []
    for k in I['kernels']:
        s = sector_of(k)
        if s not in seen:
            seen.append(s)
    return sorted(seen)


def choose_typical(I):
    cnt = collections.Counter(sector_of(k) for k in I['kernels'])
    mode = collections.Counter(cnt.values()).most_common(1)[0][0]
    cands = sorted(s for s, n in cnt.items() if n == mode)
    return cands[len(cands) // 2], mode, len(cands), dict(collections.Counter(cnt.values()))


def write_chunk(D, S, I, cdir, sectors, name='ndpent_top'):
    os.makedirs(cdir, exist_ok=True)
    keep = set(sectors)
    kernels = [k for k in I['kernels'] if sector_of(k) in keep]
    J = dict(I)
    J['kernels'] = kernels
    J['orders'] = [dict(o, kernels=[k for k in o['kernels'] if sector_of(k) in keep]) for o in I['orders']]
    json.dump(S, open(os.path.join(cdir, name + '.json'), 'w'), indent=1)
    json.dump(J, open(os.path.join(cdir, name + '_integral.json'), 'w'), indent=1)
    for so in (name + '_integral.so', 'builtin.so'):
        dst = os.path.join(cdir, so)
        if os.path.exists(dst) or os.path.islink(dst):
            os.remove(dst)
        try:
            os.link(os.path.join(D, so), dst)
        except OSError:
            os.symlink(os.path.abspath(os.path.join(D, so)), dst)
    co = os.path.join(cdir, 'coefficients')
    if os.path.islink(co) or os.path.exists(co):
        if os.path.islink(co):
            os.remove(co)
        else:
            shutil.rmtree(co)
    os.symlink(os.path.abspath(os.path.join(D, 'coefficients')), co)
    return kernels


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--disteval-dir', required=True)
    ap.add_argument('--out')
    ap.add_argument('--chunk-sectors', type=int, default=None, help='sectors per chunk; 0 = one chunk with every sector (unchunked)')
    ap.add_argument('--sectors', default=None, help='comma list: exactly these sectors as ONE chunk (the smoke)')
    ap.add_argument('--stage', default='stage')
    ap.add_argument('--choose-typical', action='store_true')
    ap.add_argument('--family', '--name', dest='name', default='ndpent_top', help='the package family NAME (sum file NAME.json, integral NAME_integral.json / .so); --name is the alias')
    a = ap.parse_args()
    D = a.disteval_dir
    probs = family_problems(D, a.name)
    if probs:
        present = sorted(f[:-len('_integral.json')] for f in os.listdir(D) if f.endswith('_integral.json')) if os.path.isdir(D) else []
        print('REFUSED: --family %s is not the family of the package in %s: %s; families present: %s' % (a.name, D, '; '.join(probs), ', '.join(present) or 'none'))
        return 2
    S, I = load(D, a.name)
    secs = sectors_in_order(I)
    if a.choose_typical:
        s, mode, ncand, hist = choose_typical(I)
        print(json.dumps({'sector': s, 'modal_kernel_count': mode, 'n_sectors_with_modal_count': ncand, 'kernels_per_sector_histogram': hist,
                          'kernels': [k for k in I['kernels'] if sector_of(k) == s], 'n_sectors': len(secs), 'n_kernels': len(I['kernels'])}))
        return 0
    assert a.out, '--out required'
    if a.sectors:
        groups = [[int(x) for x in a.sectors.split(',') if x.strip()]]
        for g in groups:
            for s in g:
                assert s in secs, 'sector %d not in the package' % s
    else:
        assert a.chunk_sectors is not None, '--chunk-sectors or --sectors'
        n = a.chunk_sectors
        groups = [secs] if n <= 0 else [secs[i:i + n] for i in range(0, len(secs), n)]
    os.makedirs(a.out, exist_ok=True)
    man = {'receipt': 'CHUNKS %s: per-sector checkpoint chunks of the disteval data dir (kernels filtered per sector group; prefactor/coefficient/sum unchanged; results add per order)' % a.stage,
           'stamp_utc': utc(), 'stage': a.stage, 'disteval_dir': os.path.abspath(D),
           'name': a.name, 'source_sha256': {a.name + '.json': sha256(os.path.join(D, a.name + '.json')), a.name + '_integral.json': sha256(os.path.join(D, a.name + '_integral.json')),
                             a.name + '_integral.so': sha256(os.path.join(D, a.name + '_integral.so')), 'builtin.so': sha256(os.path.join(D, 'builtin.so'))},
           'n_sectors': len(secs), 'n_kernels': len(I['kernels']), 'chunk_sectors': (a.chunk_sectors if not a.sectors else 'explicit'), 'n_chunks': len(groups), 'chunks': []}
    tot = 0
    for ci, g in enumerate(groups):
        cdir = os.path.join(a.out, 'chunk_%03d' % ci)
        kern = write_chunk(D, S, I, cdir, g, a.name)
        tot += len(kern)
        man['chunks'].append({'id': ci, 'dir': cdir, 'sectors': g, 'n_sectors': len(g), 'n_kernels': len(kern), 'integral_json_sha256': sha256(os.path.join(cdir, a.name + '_integral.json'))})
    man['n_kernels_in_chunks'] = tot
    man['kernels_cover_package'] = (tot == len(I['kernels'])) if not a.sectors else False
    json.dump(man, open(os.path.join(a.out, 'CHUNKS_%s.json' % a.stage), 'w'), indent=1)
    print('CHUNKS %s: %d chunks, %d sectors, %d kernels (package %d), cover=%s -> %s' % (a.stage, len(groups), sum(len(g) for g in groups), tot, len(I['kernels']), man['kernels_cover_package'], a.out))
    return 0


if __name__ == '__main__':
    sys.dont_write_bytecode = True
    raise SystemExit(main())
