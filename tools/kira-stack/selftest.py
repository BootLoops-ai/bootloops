#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision.
"""kira-stack smoke battery — run: python3 selftest.py

The engine (Kira/FireFly/Fermat) is external; this battery exercises the
shipped extension scripts on synthetic fixtures in a scratch directory, with
hand-computed expectations and fail-closed controls:

  S1  build_restage --dir: per-sector seed-pin census on a synthetic target
      list — sector bitmasks, (r,s,d) boxes, pinned-row emission, global box,
      receipt fields, and the jobs.yaml backup all checked against hand
      computation.
  M1  build_restage fail-closed: a positive ISP index in the target list must
      abort with the named assert (never a silently mis-binned census).
  S2  parallel_kira_gen helpers: cover_tops bit-drop subsectors, parse_targets
      sector extraction, and synthetic own-sector target generation (corner
      present, sector preserved, r/s bounds respected).
  S3  ffsave_degree_census on a synthetic run dir (tiny kira.db + name files):
      the layout gate (encode->decode round-trip on every name) must pass; a
      planted ff_save state tagged with the tool's own weight encoding must
      decode back to the same (target, master) names and degrees, a tag naming
      an integral outside the name files must resolve via inverse decode, and
      a ZERO state must be counted and skipped. The gate is then swept over
      the full engine ordering range (integral_ordering 1..8), and the
      pyred_weight codec is pinned against hand computation: a hand-derived
      ordering-8 weight, the four dotsp layouts separating on a dots+sps
      name, sector_ordering 1 vs 2 ranks, out-of-tree sector ranking under a
      declared top sector, and a fail-loud decode of an inconsistent weight.
  M2  ffsave_degree_census fail-closed: an out-of-range integral_ordering in
      kira.db must abort with the named assert.
  S4  parallel_kira_gen top-emit on a toy family (operator templates + a
      planted shard kira.db, no engine): the emitted SYSTEM_<fam>_<top>.gz
      must parse with read_system_eqs, contain exactly the hand-computed
      equations (head = highest weight, terms descending, engine 6-token
      term lines), and merge --top-from must consume it (has_top, eq count).
  M3  top-emit fail-closed: a kira.db integral ordering that contradicts the
      setup io must abort with the named error.

S3, M2, S4 and M3 need the two Kira-derived GPL tools (pyred_weight.py,
ffsave_degree_census.py), which ship in the sibling repository
kira under bootloops-tools/ and are found through
kira_gpl_tools.py (../kira beside this checkout, or the
directory in env BOOTLOOPS_KIRA_TOOLS). When they are absent those four legs
SKIP BY NAME and the battery still exits 0 on the legs that ran; nothing is
silently substituted.

Exit 0 all PASS (or PASS + named SKIPs); assertion failure (nonzero exit)
otherwise.
"""
import gzip
import json
import os
import re
import runpy
import sqlite3
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import kira_gpl_tools  # noqa: E402  (loader for the GPL tools in kira)

GPL_LEGS = ('S3', 'M2', 'S4', 'M3')


def s1_m1_build_restage():
    with tempfile.TemporaryDirectory() as td:
        # fam[i1..i4, i5, i6]: 4 propagators + 2 ISPs
        with open(os.path.join(td, 'de_targets'), 'w') as f:
            f.write('fam[1,1,0,0,0,0]\n'      # sector 3,  r=2 s=0 d=0
                    'fam[2,1,0,0,-1,0]\n'     # sector 3,  r=3 s=1 d=1
                    'fam[1,0,1,1,0,-2]\n')    # sector 13, r=3 s=2 d=0
        with open(os.path.join(td, 'jobs.yaml'), 'w') as f:
            f.write('old monolith\n')
        r = subprocess.run([sys.executable, os.path.join(HERE, 'build_restage.py'),
                            '--dir', td, '--fam', 'fam', '--nprops', '4'],
                           capture_output=True, text=True)
        assert r.returncode == 0, r.stdout + r.stderr
        jy = open(os.path.join(td, 'jobs.yaml')).read()
        assert '- {topologies: [fam], sectors: [13], r: 3, s: 2, d: 0}' in jy
        assert '- {topologies: [fam], sectors: [3], r: 3, s: 1, d: 1}' in jy
        assert jy.index('sectors: [13]') < jy.index('sectors: [3]'), \
            'pinned rows not in descending sector order'
        assert open(os.path.join(td, 'jobs.yaml.pre_restage')).read() == 'old monolith\n'
        rcpt = json.load(open(os.path.join(td, 'RESTAGE_STAGING_RECEIPT.json')))
        (st,) = rcpt['staged'].values()
        assert st['n_targets'] == 3 and st['n_pinned_sectors'] == 2
        assert st['global_box_rsd'] == [3, 2, 1]
        print('S1 PASS  build_restage census + pinned rows + receipt vs hand computation')

        # M1: positive ISP index must abort with the named assert
        with open(os.path.join(td, 'de_targets'), 'w') as f:
            f.write('fam[1,1,0,0,2,0]\n')
        r = subprocess.run([sys.executable, os.path.join(HERE, 'build_restage.py'),
                            '--dir', td, '--fam', 'fam', '--nprops', '4'],
                           capture_output=True, text=True)
        assert r.returncode != 0 and 'positive ISP index' in r.stdout + r.stderr, \
            'positive-ISP control did not fire'
        print('M1 PASS  build_restage fail-closed on a positive ISP index')


def s2_parallel_kira_gen():
    import parallel_kira_gen as pkg
    assert pkg.bits(0b1011) == 3
    assert pkg.cover_tops(0b1011, 4) == [3, 9, 10], 'bit-drop subsectors wrong'
    with tempfile.TemporaryDirectory() as td:
        tf = os.path.join(td, 'de_targets')
        with open(tf, 'w') as f:
            f.write('fam[1,1,0,0,-1,0]\n')
        (line, idx, sec), = pkg.parse_targets(tf)
        assert idx == (1, 1, 0, 0, -1, 0) and sec == 3
    lines = pkg.synthetic_targets('fam', [3], r=3, s=1, nidx=6, nprops_top=4)
    assert 'fam[1,1,0,0,0,0]' in lines, 'corner rep missing'
    for ln in lines:
        idx = tuple(int(x) for x in re.search(r'\[(.*)\]', ln).group(1).split(','))
        assert sum(1 << i for i in range(4) if idx[i] >= 1) == 3, 'sector drifted'
        assert sum(a for a in idx if a > 0) <= 3, 'r bound exceeded'
        assert -min(idx) <= 1, 's bound exceeded'
    print('S2 PASS  parallel_kira_gen helpers (cover_tops, parse_targets, '
          'synthetic own-sector targets)')


def _make_rundir(td, iord=5):
    d = os.path.join(td, 'run')
    os.makedirs(os.path.join(d, 'results', 'fam'))
    os.makedirs(os.path.join(d, 'ff_save', 'states'))
    con = sqlite3.connect(os.path.join(d, 'results', 'kira.db'))
    con.execute('CREATE TABLE WEIGHTBITS (b_d1 INT, b_d2 INT, b_pp INT, b_sp INT)')
    con.execute('INSERT INTO WEIGHTBITS VALUES (8, 8, 12, 12)')
    con.execute('CREATE TABLE INTEGRALORDERING (iord INT)')
    con.execute(f'INSERT INTO INTEGRALORDERING VALUES ({iord})')
    con.commit()
    con.close()
    with open(os.path.join(d, 'results', 'fam', 'masters'), 'w') as f:
        f.write('fam[1,0,1,0]\nfam[0,1,1,-1]\n')
    with open(os.path.join(d, 'de_targets'), 'w') as f:
        f.write('fam[1,1,0,0]\nfam[1,1,1,-1]\nfam[3,1,0,-1]\n')
        # fam[3,1,0,-1] carries dots=2 AND sps=1 so the round-trip gate
        # separates all four dotsp layouts
    open(os.path.join(d, 'preferred'), 'w').close()
    return d


def s3_m2_ffsave_census():
    script = kira_gpl_tools.script_path('ffsave_degree_census')
    with tempfile.TemporaryDirectory() as td:
        d = _make_rundir(td)
        # first pass (no states): runs the layout gate on every name and
        # hands us the tool's own encoder for tag construction
        argv0 = sys.argv
        sys.argv = ['ffsave_degree_census.py', d, '--family', 'fam',
                    '--out', os.path.join(td, 'c0.json')]
        try:
            g = runpy.run_path(script)
        finally:
            sys.argv = argv0
        w = g['default_weight']
        # plant states: one tag inside the name files, one target outside them
        # (inverse-decode path), one ZERO state
        recs = [(f"{w((1, 1, 0, 0))}_{w((1, 0, 1, 0))}", 37, 12),
                (f"{w((2, 1, 1, 0))}_{w((0, 1, 1, -1))}", 5, 8)]
        sdir = os.path.join(d, 'ff_save', 'states')
        for tag, dn, dd in recs:
            with gzip.open(os.path.join(sdir, f'{tag}_1.gz'), 'wt') as f:
                f.write(f'tag_name\n{tag}\nmax_deg_num\n{dn}\n'
                        f'max_deg_den\n{dd}\nis_done\n1\n')
        with gzip.open(os.path.join(sdir, 'z_1.gz'), 'wt') as f:
            f.write('ZERO\n')
        out = os.path.join(td, 'census.json')
        r = subprocess.run([sys.executable, script, d, '--family', 'fam',
                            '--out', out], capture_output=True, text=True)
        assert r.returncode == 0 and 'GATE PASS' in r.stderr, r.stdout + r.stderr
        c = json.load(open(out))
        assert c['zero_states'] == 1
        e = {t['tag']: t for t in c['entries']}
        assert len(e) == 2
        e1 = e[recs[0][0]]
        assert (e1['target'], e1['master']) == ([1, 1, 0, 0], [1, 0, 1, 0])
        assert (e1['deg_num'], e1['deg_den'], e1['done']) == (37, 12, 1)
        assert e1['target_src'] == 'file' and e1['master_src'] == 'file'
        e2 = e[recs[1][0]]
        assert e2['target'] == [2, 1, 1, 0] and e2['target_src'] == 'inverse', \
            'inverse weight decode failed to name the out-of-file target'
        print('S3 PASS  ffsave_degree_census layout gate + planted-state decode '
              '(file + inverse naming, ZERO counted)')

    # layout gate across the full engine ordering range
    for iord in range(1, 9):
        with tempfile.TemporaryDirectory() as td:
            d = _make_rundir(td, iord=iord)
            r = subprocess.run([sys.executable, script, d, '--family', 'fam',
                                '--out', os.path.join(td, 'c.json')],
                               capture_output=True, text=True)
            assert r.returncode == 0 and 'GATE PASS' in r.stderr, \
                (iord, r.stdout + r.stderr)

    # pyred_weight codec pins (hand computation; bits b_d1=8 b_d2=8 b_pp=12
    # b_sp=12 as planted in kira.db above)
    pw = kira_gpl_tools.load('pyred_weight')
    c8 = pw.WeightCodec(4, (8, 8, 12, 12), 8)
    # fam[2,1,0,-1]: sector 3 (rank 5 under (#lines, sector)), dots=1, sps=1,
    # d1=dots+sps=2, d2=sps=1, ppw=idx of (1,0) in comps(1,2)=1,
    # spw=(2-1)-idx of (0,1)=1
    assert c8.encode((2, 1, 0, -1)) == \
        ((((5 << 8 | 2) << 8 | 1) << 12 | 1) << 12 | 1), \
        'hand-computed ordering-8 weight mismatch'
    ws = set()
    for iord in (1, 2, 3, 4):
        c = pw.WeightCodec(4, (8, 8, 12, 12), iord)
        nu = (3, 1, 0, -1)          # dots=2, sps=1: separates all four layouts
        w = c.encode(nu)
        assert c.decode(w) == nu, (iord, 'round trip failed')
        ws.add(w)
    assert len(ws) == 4, 'dotsp branches do not separate'
    # sector_ordering 1 ranks by sector number; 2 by (#lines, sector number)
    assert pw.WeightCodec(4, (8, 8, 12, 12), 1).encode((1, 1, 0, 0)) >> 40 == 3
    assert c8.encode((1, 1, 0, 0)) >> 40 == 5
    # out-of-tree sectors rank after in-tree sectors of the same line count
    ct = pw.WeightCodec(4, (8, 8, 12, 12), 8, topsectors=[7])
    assert (ct.rank[8], ct.rank[3]) == (4, 5), 'out-of-tree ranking wrong'
    try:
        pw.WeightCodec(4, (8, 8, 12, 12), 3).decode((5 << 8 | 1) << 32 | 2 << 24)
        assert False, 'negative-sps decode did not fail loud'
    except ValueError:
        pass
    print('S3 PASS  ordering sweep 1..8 + pyred_weight codec pins '
          '(hand weight, dotsp separation, sector ranks, fail-loud decode)')

    with tempfile.TemporaryDirectory() as td:
        d = _make_rundir(td, iord=9)
        r = subprocess.run([sys.executable, script, d, '--family', 'fam'],
                           capture_output=True, text=True)
        assert r.returncode != 0 and 'engine range is 1..8' in r.stdout + r.stderr, \
            'out-of-range-ordering control did not fire'
        print('M2 PASS  ffsave_degree_census fail-closed on an out-of-range ordering')


def s4_m3_top_emit():
    import parallel_kira_gen as pkg
    pw = kira_gpl_tools.load('pyred_weight')
    with tempfile.TemporaryDirectory() as td:
        # toy 4-index family: 3-line top sector 7 with one ISP (index 3)
        ref = os.path.join(td, 'ref')
        smdir = os.path.join(ref, 'sectormappings', 'fam')
        os.makedirs(os.path.join(ref, 'config'))
        os.makedirs(smdir)
        with open(os.path.join(smdir, 'IBP'), 'w') as f:
            f.write('fam[0,0,0,0]*(d-2*a0)\n'    # op1: (d-2n0) I[n] - I[n+e0]
                    'fam[1,0,0,0]*(-1)\n'
                    '\n'
                    'fam[0,0,0,0]*(a3)\n'        # op2: n3 I[n] + m2 n3 I[n-e3]
                    'fam[0,0,0,-1]*(m2*a3)\n')
        with open(os.path.join(smdir, 'trivialsector'), 'w') as f:
            f.write('0\n')
        with open(os.path.join(ref, 'de_targets'), 'w') as f:
            f.write('fam[1,1,1,0]\n')
        out = os.path.join(td, 'out')
        pkg.setup(ref, out, 'fam', top=7, r=3, s=1, tops=[3, 5, 6],
                  targets_file='de_targets', io=8, no_pm=True,
                  targets_mode='own-sectors', nidx=4)
        # plant a completed shard's kira.db (the weight-bits source)
        dbp = os.path.join(out, 'shard_3', 'results', 'kira.db')
        con = sqlite3.connect(dbp)
        con.execute('CREATE TABLE WEIGHTBITS (b_d1 INT, b_d2 INT, b_pp INT, b_sp INT)')
        con.execute('INSERT INTO WEIGHTBITS VALUES (8, 8, 12, 12)')
        con.execute('CREATE TABLE INTEGRALORDERING (iord INT)')
        con.execute('INSERT INTO INTEGRALORDERING VALUES (8)')
        con.commit()
        con.close()
        rcpt = pkg.emit_top(out)
        path = os.path.join(out, 'top_emit', 'tmp', 'fam', 'SYSTEM_fam_7.gz')
        assert rcpt['n_eqs'] == 3 and os.path.abspath(path) == rcpt['file']
        eqs = pkg.read_system_eqs(path)
        # hand computation: seeds (1,1,1,0) and (1,1,1,-1); op2 vanishes on
        # the n3=0 seed -> 3 equations, each block in the engine layout
        # (head = highest weight, terms descending, 6-token term lines)
        w = pw.WeightCodec(4, (8, 8, 12, 12), 8, topsectors=[7]).encode
        exp = {
            # op1 @ (1,1,1,0):  -I[2,1,1,0] + (d-2) I[1,1,1,0]
            w((2, 1, 1, 0)): (f'Eq\n{w((2,1,1,0))}\n2\n'
                              f'2 -1 {w((2,1,1,0))} 7 0 0\n'
                              f'2 1*d-2 {w((1,1,1,0))} 7 0 0'),
            # op1 @ (1,1,1,-1): -I[2,1,1,-1] + (d-2) I[1,1,1,-1]
            w((2, 1, 1, -1)): (f'Eq\n{w((2,1,1,-1))}\n2\n'
                               f'2 -1 {w((2,1,1,-1))} 7 0 0\n'
                               f'2 1*d-2 {w((1,1,1,-1))} 7 0 0'),
            # op2 @ (1,1,1,-1): -m2 I[1,1,1,-2] - I[1,1,1,-1]
            w((1, 1, 1, -2)): (f'Eq\n{w((1,1,1,-2))}\n2\n'
                               f'2 -1*m2 {w((1,1,1,-2))} 7 0 0\n'
                               f'2 -1 {w((1,1,1,-1))} 7 0 0'),
        }
        got = {h: b for h, b in eqs}
        assert got == exp, f'emitted blocks differ from hand computation:' \
                           f'\n{got}\nvs\n{exp}'
        # merge consumes the emit dir through the documented --top-from route
        target = os.path.join(td, 'target')
        m = pkg.merge(out, target, top_from=os.path.join(out, 'top_emit'),
                      dedupe=True)
        assert m['has_top'] and m['n_eqs'] == 3
        assert open(os.path.join(target, 'tmp', 'fam', 'SYSTEMconfig')
                    ).read().strip() == '3'
        print('S4 PASS  top-emit engine-format SYSTEM file vs hand computation '
              '+ merge --top-from consumption')

        # M3: kira.db ordering contradicting the setup io must abort
        con = sqlite3.connect(dbp)
        con.execute('UPDATE INTEGRALORDERING SET iord = 5')
        con.commit()
        con.close()
        try:
            pkg.emit_top(out)
            assert False, 'io-mismatch control did not fire'
        except SystemExit as e:
            assert 'integral_ordering' in str(e), e
        print('M3 PASS  top-emit fail-closed on a kira.db/setup ordering mismatch')


if __name__ == '__main__':
    s1_m1_build_restage()
    s2_parallel_kira_gen()
    if kira_gpl_tools.available():
        s3_m2_ffsave_census()
        s4_m3_top_emit()
        print('kira-stack selftest: ALL PASS (7/7)')
    else:
        # the top-emit route must refuse by name, never fall back silently
        try:
            kira_gpl_tools.load('pyred_weight')
            assert False, 'loader did not raise on an absent tool'
        except ImportError as e:
            assert 'kira' in str(e) and kira_gpl_tools.ENV in str(e), e
        for leg in GPL_LEGS:
            print(f'{leg} SKIP  Kira-derived GPL tools (pyred_weight, '
                  f'ffsave_degree_census) not found in {kira_gpl_tools.tools_dir()} '
                  f'— check out ../kira or set {kira_gpl_tools.ENV}')
        print('kira-stack selftest: ALL PASS (3/3 run; S3 M2 S4 M3 SKIPPED by '
              'name, named ImportError verified)')
