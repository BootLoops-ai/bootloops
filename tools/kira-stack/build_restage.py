#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision.
"""build_restage.py — per-sector seed-pin RESTAGE stager for kira jobs.yaml
(kira-stack member).

SYMPTOM this fixes: kira equation GENERATION censored/slow on a MONOLITHIC
reduce row {sectors: [TOP], r,s,d = full-shard box} — the full-box seed sweep
applies the top box to every sector of a ~2,500-3,000-sector tree.
MECHANISM: emit ONE reduce row PER TARGET SECTOR, each with the (r,s,d) box
MEASURED from the targets in that sector (s is the explosive dimension, d
next; measured pin collapse 313G->7.4G, 40x). Sectors
outside the target support get no row; tree sectors below a pinned row
inherit only that row's small box. select_mandatory_list stays the target
list — sweep COUNT is unchanged (closure invariance), the per-sector seed box
collapse owns the win.

VALIDATION (measured, same shard / same box):
- Monolithic staging: on a ~2,800-sector tree, all 4 heavy shards censored at
  a 10,800 s (3 h) generation cap with 0 SYSTEM_*.gz written and ~22 h sunk;
  one shard censored at sector 2634/2826.
- Restaged (SAME shard, SAME box): the full 2826/2826 tree sweep completed
  inside a 7,200 s window; measured early sweep rate ~2.9 sect/s vs the
  monolith's 0.24 avg — a ~12x datum (~46 min sweep vs the censored 3 h cap).
- --dir mode: a monolithic target_reduce jobs.yaml, censored in generation at
  an r15 s2 d2 full-tree sweep, restaged to 591 pinned rows with the global
  box collapsed r15s2d2->r14s2d2.

USAGE (fail-closed: no silent defaults):
  build_restage.py --base DIR --fam FAMILY TOP [TOP ...]      # shard layout:
      stages DIR/shard_<TOP>/jobs.yaml from DIR/shard_<TOP>/<targets>
  build_restage.py --dir DIR --fam FAMILY                     # single dir:
      stages DIR/jobs.yaml from DIR/<targets> (monolithic-run shape)
Options:
  --targets NAME   targets filename inside each staged dir (default: de_targets)
  --nprops N       number of propagator indices (sector bitmask + ISP split;
                   default 14; MUST match your family or the ISP assert fires)
  --head-file F / --tail-file F   override the jobs.yaml head/tail blocks
                   (defaults emit: select_mandatory_list = [FAM, <targets>],
                   integral_ordering 5, run_initiate true,
                   run_triangular/run_firefly false). A run that needs
                   preferred_masters + run_firefly: true is a --tail-file job,
                   not an edit here.
  --backup-suffix S  (default .pre_restage; existing jobs.yaml renamed ONCE)
Writes RESTAGE_STAGING_RECEIPT.json (per-dir pin census + jobs.yaml sha256)
into --base (or --dir). Readback-only otherwise; touches nothing outside the
staged dirs. Targets format: kira integral lines `Fam[..., i1, ..., iN]`;
positive ISP indices are asserted absent (they would break the r/s/d census).
"""
import argparse
import hashlib
import json
import os
import re
import sys
import time

HEAD_DEFAULT = """jobs:
  - reduce_sectors:
      reduce:
"""
TAIL_DEFAULT = """      select_integrals:
        select_mandatory_list:
          - [{fam}, {targets}]
      integral_ordering: 5
      run_initiate: true
      run_triangular: false
      run_firefly: false
"""


def census(path, nprops):
    """sector -> [n, maxr, maxs, maxd] over the target list at `path`."""
    secs = {}
    for line in open(path):
        m = re.search(r'\[(.*)\]', line.strip())
        if not m:
            continue
        idx = [int(x) for x in m.group(1).split(',')]
        assert not any(a > 0 for a in idx[nprops:]), \
            f'positive ISP index in {path} (check --nprops): {line}'
        sec = sum(1 << i for i, a in enumerate(idx[:nprops]) if a > 0)
        r = sum(a for a in idx if a > 0)
        s = -sum(a for a in idx if a < 0)
        d = sum(a - 1 for a in idx if a > 1)
        e = secs.setdefault(sec, [0, 0, 0, 0])
        e[0] += 1
        e[1] = max(e[1], r)
        e[2] = max(e[2], s)
        e[3] = max(e[3], d)
    return secs


def stage(sd, fam, targets, head, tail, nprops, backup_suffix):
    tpath = os.path.join(sd, targets)
    secs = census(tpath, nprops)
    rows = []
    for sec in sorted(secs, reverse=True):
        n, r, s, d = secs[sec]
        rows.append(f'        - {{topologies: [{fam}], sectors: [{sec}], '
                    f'r: {r}, s: {s}, d: {d}}}\n')
    jy = head + ''.join(rows) + tail
    jp = os.path.join(sd, 'jobs.yaml')
    bak = jp + backup_suffix
    if os.path.exists(jp) and not os.path.exists(bak):
        os.rename(jp, bak)
    with open(jp, 'w') as f:
        f.write(jy)
    gbox = (max(v[1] for v in secs.values()),
            max(v[2] for v in secs.values()),
            max(v[3] for v in secs.values()))
    return {
        'n_targets': sum(v[0] for v in secs.values()),
        'n_pinned_sectors': len(secs),
        'global_box_rsd': list(gbox),
        'jobs_sha256': hashlib.sha256(jy.encode()).hexdigest(),
        'jobs_yaml': jp,
        'backup': bak if os.path.exists(bak) else None,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--base', help='shard layout root (stages shard_<TOP>/)')
    ap.add_argument('--dir', dest='onedir', help='single staged dir (monolith shape)')
    ap.add_argument('--fam', required=True)
    ap.add_argument('--targets', default='de_targets')
    ap.add_argument('--nprops', type=int, default=14)
    ap.add_argument('--head-file')
    ap.add_argument('--tail-file')
    ap.add_argument('--backup-suffix', default='.pre_restage')
    ap.add_argument('tops', nargs='*', type=int)
    a = ap.parse_args()
    if bool(a.base) == bool(a.onedir):
        ap.error('exactly one of --base (with TOP args) or --dir')
    if a.base and not a.tops:
        ap.error('--base needs TOP sector numbers')
    head = open(a.head_file).read() if a.head_file else HEAD_DEFAULT
    tail = (open(a.tail_file).read() if a.tail_file
            else TAIL_DEFAULT.format(fam=a.fam, targets=a.targets))
    root = a.base or a.onedir
    rcpt_path = os.path.join(root, 'RESTAGE_STAGING_RECEIPT.json')
    rcpt = json.load(open(rcpt_path)) if os.path.exists(rcpt_path) else {}
    rcpt.setdefault('tool', 'tools/kira-stack/build_restage.py (seed-pin restage)')
    rcpt.setdefault('recipe', 'one reduce row per target sector, (r,s,d) measured '
                    'from that sector\'s targets')
    rcpt.setdefault('staged', {})
    units = ([(str(t), os.path.join(a.base, f'shard_{t}')) for t in a.tops]
             if a.base else [(os.path.basename(root.rstrip('/')), root)])
    for name, sd in units:
        rcpt['staged'][name] = stage(sd, a.fam, a.targets, head, tail,
                                     a.nprops, a.backup_suffix)
        rcpt['staged'][name]['staged_utc'] = time.strftime(
            '%Y-%m-%dT%H:%M:%SZ', time.gmtime())
        print(f"[restage] {name}: {rcpt['staged'][name]['n_pinned_sectors']} "
              f"pinned sectors, {rcpt['staged'][name]['n_targets']} targets, "
              f"box {rcpt['staged'][name]['global_box_rsd']}")
    rcpt['updated_utc'] = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
    json.dump(rcpt, open(rcpt_path, 'w'), indent=1)
    print(f'[restage] receipt -> {rcpt_path}')


if __name__ == '__main__':
    main()
