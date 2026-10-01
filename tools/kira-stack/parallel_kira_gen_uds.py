#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision.
"""parallel_kira_gen --mode=uds : generate IBP equations per-sector directly
from Kira's own sectormappings/{IBP,LI} operator templates, write in
`TOPO[idx,...]*(coeff)` blank-line-separated format, and feed to
`reduce_user_defined_system + input_system:` — Kira does the pyred weight
encoding + Phase-C SYSTEM_*.gz emission itself. Bypasses Phase-A entirely
AND covers the top sector.

Operator template semantics (Kira's sectormappings/<fam>/IBP, blank-line
separated ops; each line `TOPO[s1,...,s15]*(poly(a1..a15,d,m2))`):
applied to a seed n = (n1..n15), substitute a_i -> n_i (integer), add shift
s to the seed indices, and the equation is Σ_terms coeff * I[n+s] = 0.
This is the same convention pyred uses (relations.cpp generate_system).

Coeffs after substitution are polys in {d, m2} with RATIONAL constants —
exactly what Kira's read_equations / FireFly parser accept.

CLI:
  gen  REF OUT [--fam F --top T --r R --s S --sectors S1,.. --nproc N]
  uds  OUT [--fam F]     # write jobs_uds.yaml + de_targets, ready for kira
"""
import argparse, gzip, json, multiprocessing as mp, os, re, sys, time
from fractions import Fraction as F

FAM_RE = re.compile(r'^(\w+)\[([\d, +-]+)\]\*\((.+)\)$')
AI_RE = re.compile(r'\ba(\d+)\b')


# ----------------------------------------------------------------- template

def load_ops(ref, fam):
    """→ list of ops; each op = list of (shift:tuple[int], coeff_terms).
    coeff_terms = list of (const:Fraction, dpow:int, m2pow:int, aidx:int|None)
    from a factored monomial (Kira's IBP/LI coeffs are sums of such)."""
    ops = []
    for name in ('IBP', 'LI'):
        path = f'{ref}/sectormappings/{fam}/{name}'
        if not os.path.exists(path): continue
        raw = open(path).read().split('\n\n')
        for block in raw:
            block = block.strip()
            if not block: continue
            terms = []
            for line in block.split('\n'):
                m = FAM_RE.match(line.strip())
                if not m: continue
                shift = tuple(int(x) for x in m.group(2).split(','))
                coeff = _parse_coeff(m.group(3))
                terms.append((shift, coeff))
            if terms: ops.append(terms)
    return ops


def _parse_coeff(s):
    """Kira IBP coeffs are sums of monomials in {a1..a15, d, m2} with small
    integer/rational constants, e.g. '-2*a8-a12+d-a4', '2*m2*a1', '3/2*a11'.
    → list of (const, dpow, m2pow, aidx|None). Exactly one a-factor per
    monomial (or none) — verified against a production Kira IBP template."""
    out = []
    # split into signed monomials
    s = s.replace(' ', '')
    parts = re.findall(r'[+-]?[^+-]+', s)
    for part in parts:
        sign = F(1)
        if part[0] == '-': sign = F(-1); part = part[1:]
        elif part[0] == '+': part = part[1:]
        facs = part.split('*')
        const = sign; dpow = 0; m2pow = 0; aidx = None
        for f in facs:
            if f == 'd': dpow += 1
            elif f == 'm2': m2pow += 1
            elif f.startswith('a') and f[1:].isdigit():
                assert aidx is None, f'>1 a-factor in {part!r}'
                aidx = int(f[1:])  # Kira IBP uses a0..a{np-1} (0-indexed)
            else:
                const *= F(f)
        out.append((const, dpow, m2pow, aidx))
    return out


def _fmt(c, dpow, m2pow):
    """rational const * d^dpow * m2^m2pow → Kira/FireFly-parseable string."""
    num, den = c.numerator, c.denominator
    s = f'{num}' if den == 1 else f'{num}/{den}'
    for _ in range(dpow): s += '*d'
    for _ in range(m2pow): s += '*m2'
    return s


# ----------------------------------------------------------------- seeds

def seeds_for_sector(S, nprops, r, s):
    """all seeds in sector S with sum(pos) <= r, sum(-neg) <= s, in-sector
    props >= 1, others <= 0. Matches Kira's SeedSpecSector::expand corner
    (dots = r - lines, sps = s)."""
    lines = [i for i in range(nprops) if S >> i & 1]
    isp = [i for i in range(nprops) if not (S >> i & 1)]
    dots = r - len(lines)
    if dots < 0: return
    base = [1 if S >> i & 1 else 0 for i in range(nprops)]
    # enumerate dot compositions
    dot_pats = [[]]
    def _dots(rem, i, cur):
        if i == len(lines): dot_pats.append(list(cur)); return
        for e in range(rem + 1):
            _dots(rem - e, i + 1, cur + [e])
    dot_pats.clear(); _dots(dots, 0, [])
    # sps compositions over ISP slots
    sps_pats = [[]]
    def _sps(rem, i, cur):
        if i == len(isp): sps_pats.append(list(cur)); return
        for e in range(rem + 1):
            _sps(rem - e, i + 1, cur + [e])
    sps_pats.clear(); _sps(s, 0, [])
    for dp in dot_pats:
        for sp in sps_pats:
            n = list(base)
            for k, i in enumerate(lines): n[i] += dp[k]
            for k, i in enumerate(isp): n[i] -= sp[k]
            yield tuple(n)


# ----------------------------------------------------------------- gen

_G = {}  # per-worker: ops, fam, nprops, triv


def _init(ref, fam, nprops):
    _G['ops'] = load_ops(ref, fam)
    _G['fam'] = fam; _G['np'] = nprops
    tp = f'{ref}/sectormappings/{fam}/trivialsector'
    _G['triv'] = set(int(x) for x in open(tp).read().replace('\n', ',').split(',')
                     if x.strip()) if os.path.exists(tp) else set()


def _sector_of(t, nprops):
    return sum(1 << i for i in range(nprops) if t[i] >= 1)


def gen_sector(args):
    """→ (S, n_eqs, path). Writes eq_S.txt.gz in TOPO[idx]*(coeff)\\n\\n."""
    S, r, s, out = args
    ops, fam, nprops, triv = _G['ops'], _G['fam'], _G['np'], _G['triv']
    path = f'{out}/eq_{S}.txt.gz'
    n_eqs = 0
    with gzip.open(path, 'wt') as f:
        for seed in seeds_for_sector(S, nprops, r, s):
            for op in ops:
                # accumulate coeff by target tuple
                acc = {}
                for shift, cterms in op:
                    # substitute a_i -> seed[i]; drop monomial if seed[aidx]==0
                    for const, dpow, m2pow, aidx in cterms:
                        c = const if aidx is None else const * seed[aidx]
                        if c == 0: continue
                        tgt = tuple(seed[i] + shift[i] for i in range(nprops))
                        key = (tgt, dpow, m2pow)
                        acc[key] = acc.get(key, F(0)) + c
                # group by tgt
                by_tgt = {}
                for (tgt, dpow, m2pow), c in acc.items():
                    if c == 0: continue
                    # drop trivial-sector targets (they're exact zeros)
                    if _sector_of(tgt, nprops) in triv: continue
                    by_tgt.setdefault(tgt, []).append((c, dpow, m2pow))
                if len(by_tgt) < 2: continue
                for tgt, mons in by_tgt.items():
                    coeff = ''.join(('+' if m[0] > 0 and i > 0 else '') + _fmt(*m)
                                    for i, m in enumerate(mons))
                    f.write(f'{fam}[{",".join(str(x) for x in tgt)}]*({coeff})\n')
                f.write('\n')
                n_eqs += 1
    return S, n_eqs, path


# ----------------------------------------------------------------- driver

def gen(ref, out, fam, top, r, s, sectors, nproc):
    os.makedirs(out, exist_ok=True)
    # nprops = the IBP first line's tuple length
    with open(f'{ref}/sectormappings/{fam}/IBP') as fh:
        first = fh.readline()
    nprops = len(FAM_RE.match(first).group(2).split(','))
    if sectors is None:
        # all non-trivial subsectors of top (Kira's own list)
        nts = set()
        for line in open(f'{ref}/sectormappings/{fam}/nonTrivialSector'):
            p = line.split()
            if len(p) >= 1:
                sn = int(p[0])
                if sn <= top and (sn & top) == sn:
                    nts.add(sn)
        sectors = sorted(nts)
    _init(ref, fam, nprops)  # main proc for meta
    t0 = time.time()
    tasks = [(S, r, s, out) for S in sectors]
    results = []
    if nproc == 1:
        for t in tasks: results.append(gen_sector(t))
    else:
        with mp.Pool(nproc, initializer=_init, initargs=(ref, fam, nprops)) as p:
            for r_ in p.imap_unordered(gen_sector, tasks):
                results.append(r_)
                if len(results) % 20 == 0:
                    print(f'  [{len(results)}/{len(tasks)}] '
                          f'{time.time()-t0:.1f}s', flush=True)
    dt = time.time() - t0
    tot = sum(n for _, n, _ in results)
    meta = {'ref': os.path.abspath(ref), 'fam': fam, 'top': top, 'r': r, 's': s,
            'nprops': nprops, 'n_ops': len(_G['ops']),
            'n_sectors': len(sectors), 'n_eqs': tot, 'wall_s': round(dt, 2),
            'files': [p for _, _, p in sorted(results)], 'nproc': nproc}
    json.dump(meta, open(f'{out}/UDS_META.json', 'w'), indent=1)
    print(f'[gen] {len(sectors)} sectors, {tot} eqs, {dt:.1f}s ({nproc} procs)')
    return meta


UDS_JOBS = """jobs:
  - reduce_user_defined_system:
      input_system:
        files: {files}
        otf: true
      select_integrals:
        select_mandatory_list:
          - [{fam}, de_targets]
{pm_line}      integral_ordering: {io}
      run_initiate: true
      run_triangular: false
      run_firefly: false
"""


def uds(out, fam, io=8, no_pm=False):
    # io + no-pm params: merged SYSTEM weight IDs are only deterministic if
    # integral_ordering (and the pm choice) match the shard gens — the
    # template must carry io explicitly.
    meta = json.load(open(f'{out}/UDS_META.json'))
    files = json.dumps([os.path.relpath(f, out) for f in meta['files']])
    pm = '' if no_pm else '      preferred_masters: preferred_masters\n'
    with open(f'{out}/jobs_uds.yaml', 'w') as f:
        f.write(UDS_JOBS.format(files=files, fam=fam, io=io, pm_line=pm))
    print(f'[uds] wrote {out}/jobs_uds.yaml ({len(meta["files"])} files, '
          f'io={io}, no_pm={no_pm})')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('gen')
    p.add_argument('ref'); p.add_argument('out')
    p.add_argument('--fam', required=True, help='family name'); p.add_argument('--top', type=int, default=1023)
    p.add_argument('--r', type=int, required=True); p.add_argument('--s', type=int, required=True)
    p.add_argument('--sectors'); p.add_argument('--nproc', type=int, default=1)
    p = sub.add_parser('uds'); p.add_argument('out'); p.add_argument('--fam', required=True, help='family name')
    p.add_argument('--io', type=int, default=8); p.add_argument('--no-pm', action='store_true')
    a = ap.parse_args()
    if a.cmd == 'gen':
        secs = [int(x) for x in a.sectors.split(',')] if a.sectors else None
        gen(a.ref, a.out, a.fam, a.top, a.r, a.s, secs, a.nproc)
    else:
        uds(a.out, a.fam, a.io, a.no_pm)
