#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision.
"""parallel_kira_gen — parallelize Kira's initiate (IBP-generation) phase by
sharding the sector tree across N independent kira instances, then merging
their SYSTEM_*.gz files for a single `run_initiate:false, run_firefly:true`
solve.

Root cause of the serial bottleneck (measured on a 15-index 4-point family,
r=13 s=1 generation, 2.6 h):
  Phase A (pyred generate_solve): worker threads generate per-sector mod-p
    equations in parallel, but the on-the-fly forward SOLVE runs in the master
    thread (relations.h:984-991 GeneratorParallelMasterSO::operator()) — serial.
  Phase C (Regenerate → SYSTEM_*.gz): treatcoeff = kira.cpp treatcoeff2 uses
    fermat[0] and is wrapped in a global mutex s_treatcoeff_mtx
    (relations.h:1242) — serial by design.
  Net: `--parallel=N` spawns N idle fer64's; kira runs at ~1.3 cores.
  (Verified on a live production run: 133% CPU, 28× fer64 @ 0.0%.)

Route B (this tool): shard by (nprops-1)-line "cover" tops.
  - pyred weight IDs are DETERMINISTIC given identical config/, sectormappings/,
    integral_ordering, preferred_masters (verified on two shards of the same
    reference run: 89 shared SYSTEM sectors, common eqs byte-identical).
  - reductionFF.cpp:load_ff_system() scans tmp/<fam>/SYSTEM_<fam>_<S>[_<k>].gz
    for every S in 0..2^jule; SYSTEMconfig is only bb->reserve() capacity.
    → merged shard files with `_k` suffixes are consumed transparently.
  - each shard = full run_initiate on its subtree only → both Phase-A and
    Phase-C work is subtree-local and runs in parallel across shards.
  Measured prior: a 9-line top at r=12 s=0 initiates in 97 s on one shard.

Known gap: the true top sector (all lines) is not a subsector of any cover top.
  Its SYSTEM file must come from ONE of:
    (a) --top-from DIR   copy SYSTEM_<fam>_<TOP>.gz from a completed run
                         (e.g. a completed full-tree oracle run);
    (b) --include-top    run one extra full-tree shard (long pole ≈ serial time
                         — defeats the point except as a correctness fallback);
    (c) top-emit         route-A single-sector emit via the pyred_weight codec:
                         apply the sectormappings IBP/LI operator templates to
                         every top-sector seed in the (r, s) box (the machinery
                         of parallel_kira_gen_uds), weight-encode heads/terms
                         with pyred_weight.py, and write SYSTEM_<fam>_<TOP>.gz
                         directly — no engine run at emit time. The codec is a
                         GPL-3.0-or-later tool (derived from Kira's sources)
                         shipped in the sibling repository kira,
                         bootloops-tools/; it is found through
                         kira_gpl_tools.py (../kira beside this
                         checkout, or env BOOTLOOPS_KIRA_TOOLS) and this route
                         stops with a named ImportError without it. Consume via
                         `merge OUT TARGET --top-from OUT/top_emit`. Raw emit =
                         no mod-p selection and no symmetry-row substitution
                         (over-generation, same direction as own-sectors mode);
                         gate the merged solve as always.

CLI:
  setup    REF OUT [--fam F --top T --r R --s S --tops S1,..]
  launch   OUT [--parallel P --nice N --host HOST]
  status   OUT
  top-emit OUT [--db RUNDIR --out-dir DIR]      # option (c), engine-free
  merge    OUT TARGET [--top-from DIR] [--dedupe]
  probe    REF [--r R --s S --tops S1,..]       # setup+launch+wait+verify, small
"""
import argparse, gzip, hashlib, json, os, re, shutil, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))


# ----------------------------------------------------------------- helpers

def bits(s): return bin(s).count('1')

def read_nontrivial(ref, fam):
    p = f'{ref}/sectormappings/{fam}/nonTrivialSector'
    out = []
    for line in open(p):
        parts = line.split()
        if len(parts) >= 2:
            out.append((int(parts[0]), int(parts[1])))
    return out

def parse_targets(path):
    """FAM[a1,...,a15] → (indices tuple, sector)."""
    out = []
    for line in open(path):
        m = re.search(r'\[([\d, -]+)\]', line)
        if not m: continue
        idx = tuple(int(x) for x in m.group(1).split(','))
        s = sum(1 << i for i in range(len(idx)) if idx[i] >= 1)
        out.append((line.rstrip('\n'), idx, s))
    return out

def cover_tops(top, nprops):
    """all (nprops-1)-line subsectors of top (bit-drop each set bit)."""
    return sorted(top ^ (1 << i) for i in range(nprops) if top >> i & 1)

def stamp(): return time.strftime('%Y-%m-%d %H:%M:%S %Z')

def _load_sibling(name):
    """Import a module from this script's own directory by path — immune to
    the flat compat shims (importing a shim would exec it and clobber
    sys.argv[0])."""
    import importlib.util
    p = os.path.join(HERE, name + '.py')
    spec = importlib.util.spec_from_file_location(f'_kira_stack_{name}', p)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m

def read_names(path, fam):
    """FAM[a1,...,aN] lines -> list of index tuples (file order)."""
    out = []
    for ln in open(path):
        m = re.search(re.escape(fam) + r'\[([\d, -]+)\]', ln)
        if m:
            out.append(tuple(int(x) for x in m.group(1).split(',')))
    return out


# ----------------------------------------------------------------- setup

JOBS_TMPL = """jobs:
  - reduce_sectors:
      reduce:
        - {{topologies: [{fam}], sectors: [{secs}], r: {r}, s: {s}{d_part}}}
      select_integrals:
        select_mandatory_list:
          - [{fam}, de_targets]
{pm_line}      integral_ordering: {io}
      run_initiate: true
      run_triangular: false
      run_firefly: false
"""
PM_LINE = '      preferred_masters: preferred_masters\n'


def synthetic_targets(fam, secs, r, s, nidx, nprops_top):
    """own-sector de_targets: for each
    ASSIGNED sector, corner + max-dot rep at the r bound + single-ISP reps at
    the s bound. Cures the containment-filter zero-output bug when all real
    targets live in the top sector. Over-generation is correctness-safe
    (gate compare direction is ref ⊆ merged); under-generation is detectable
    (masters_check + row oracle)."""
    lines = []
    isps = list(range(nprops_top, nidx))
    for S in secs:
        pbits = [i for i in range(nprops_top) if S >> i & 1]
        corner = [0] * nidx
        for i in pbits:
            corner[i] = 1
        variants = [corner]
        extra = r - len(pbits)
        if extra > 0:
            maxdot = list(corner)
            maxdot[pbits[0]] += extra          # deterministic: dots on lowest line
            variants.append(maxdot)
        for base in list(variants):
            for k in isps[:s and len(isps)]:   # s=0 → no ISP reps
                v = list(base)
                v[k] = -s
                variants.append(v)
        for v in variants:
            lines.append(f'{fam}[{",".join(str(x) for x in v)}]')
    return lines


def setup(ref, out, fam, top, r, s, tops, targets_file, io,
          no_pm=False, targets_mode='filter', manifest=None, nidx=15, d=None):
    # d: optional dot bound for the reduce box — the template previously had
    # no d key, which silently dropped dotted seeds for d>=1 references.
    # None = omit (byte-identical to prior behavior). NB captured as d_part
    # IMMEDIATELY: the shard loop reuses the name `d` for the shard dir
    # (shadowing hazard — d would render as the shard path).
    d_part = f', d: {d}' if d is not None else ''
    os.makedirs(out, exist_ok=True)
    all_targets = parse_targets(f'{ref}/{targets_file}')
    # shard spec: manifest = [{'shard_id': .., 'tops': [..]}, ..] (multi-top);
    # else legacy one-shard-per-top from tops list.
    if manifest is not None:
        shard_list = [(m['shard_id'], list(m['tops'])) for m in manifest]
    else:
        if tops is None:
            nprops = bits(top)
            tops = cover_tops(top, nprops)
        shard_list = [(str(S), [S]) for S in tops]
    meta = {'ref': os.path.abspath(ref), 'fam': fam, 'top': top, 'r': r, 's': s,
            'io': io, 'targets_file': targets_file,
            'tops': [S for _, ts in shard_list for S in ts],
            'no_pm': no_pm, 'targets_mode': targets_mode,
            'n_targets': len(all_targets), 'shards': {}, 'created': stamp()}
    for sid, secs in shard_list:
        d = f'{out}/shard_{sid}'
        os.makedirs(d, exist_ok=True)
        # config, sectormappings: MUST be identical across shards → symlink
        for sub in ('config', 'sectormappings'):
            dst = f'{d}/{sub}'
            if os.path.islink(dst) or os.path.exists(dst): continue
            os.symlink(os.path.abspath(f'{ref}/{sub}'), dst)
        if not no_pm:   # pm at gen is the wrong default — see the no-pm footgun
            shutil.copy(f'{ref}/preferred_masters', f'{d}/preferred_masters')
        if targets_mode == 'own-sectors':
            sub_t = synthetic_targets(fam, secs, r, s, nidx, bits(top))
        else:
            # de_targets: subset in this shard's subtree (empty is fine — kira
            # then selects nothing and writes 0 eqs; but we filter to nonempty)
            sub_t = [line for (line, idx, sec) in all_targets
                     if any((sec & S) == sec for S in secs)]
        with open(f'{d}/de_targets', 'w') as f:
            f.write('\n'.join(sub_t) + ('\n' if sub_t else ''))
        with open(f'{d}/jobs.yaml', 'w') as f:
            f.write(JOBS_TMPL.format(fam=fam, secs=','.join(str(S) for S in secs),
                                     r=r, s=s, io=io, d_part=d_part,
                                     pm_line='' if no_pm else PM_LINE))
        os.makedirs(f'{d}/tmp', exist_ok=True)
        os.makedirs(f'{d}/results', exist_ok=True)
        meta['shards'][sid] = {'dir': os.path.abspath(d), 'tops': secs,
                               'n_targets': len(sub_t)}
    json.dump(meta, open(f'{out}/PGEN_META.json', 'w'), indent=1)
    print(f'[setup] {len(shard_list)} shards in {out} (targets/shard: '
          f'{[meta["shards"][sid]["n_targets"] for sid, _ in shard_list]})')
    return meta


# ----------------------------------------------------------------- launch

def pid_alive(pid):
    """True if pid is a live process. Signal-0 probe — portable (/proc is
    Linux-only); PermissionError means alive but foreign-owned."""
    try:
        pid = int(pid)
    except (TypeError, ValueError):
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def launch(out, parallel, nice_lvl, host, env_extra, kira_bin='kira',
           as_cap_gb=None, max_concurrent=None):
    """as_cap_gb: prlimit --as per shard (Phase-A r=13 s=1 peaks ~100GB VSZ;
    set the HARD limit generously — raising later without root is impossible).
    max_concurrent: launch in batches of N, wait for each batch's SYSTEMconfig."""
    meta = json.load(open(f'{out}/PGEN_META.json'))
    fam = meta['fam']
    env = os.environ.copy()
    env.setdefault('FERMATPATH', 'fer64')  # resolve fer64 on PATH; env overrides
    for kv in env_extra:
        k, v = kv.split('=', 1); env[k] = v
    if not host and shutil.which('setsid') is None:
        sys.exit('FATAL: launch needs setsid(1) (util-linux; absent on stock '
                 'macOS) — its failure would be SILENT here (the shell error '
                 'lands in the redirected log while echo $! still prints a '
                 'pid). Launch on a Linux host (--host) or install setsid.')
    pids = {}
    shard_items = list(meta['shards'].items())
    batches = ([shard_items[i:i+max_concurrent]
                for i in range(0, len(shard_items), max_concurrent)]
               if max_concurrent else [shard_items])
    for bi, batch in enumerate(batches):
        for S, sh in batch:
            d = sh['dir']
            log = f'{d}/kira.log'
            cmd = f'cd {d} && exec {kira_bin} jobs.yaml --parallel={parallel}'
            if host:
                full = (f"ssh {host} 'setsid nohup bash -lc "
                        f"\"{cmd}\" </dev/null > {log} 2>&1 & echo $!'")
            else:
                full = (f'setsid nohup nice -n {nice_lvl} bash -c '
                        f'"{cmd}" </dev/null > {log} 2>&1 & echo $!')
            p = subprocess.run(full, shell=True, capture_output=True, text=True, env=env)
            pid = p.stdout.strip()
            pids[S] = pid
            if as_cap_gb and pid.isdigit():
                if sys.platform != 'linux':
                    sys.exit(f'FATAL: the address-space cap (--as) needs '
                             f'Linux (/proc cwd match + prlimit), absent on '
                             f'{sys.platform}; relaunch without --as')
                time.sleep(0.5)
                for cand in subprocess.getoutput('pgrep -f "jobs.yaml"').split():
                    if subprocess.getoutput(f'readlink /proc/{cand}/cwd') == d:
                        subprocess.run(['prlimit', f'--pid={cand}',
                                        f'--as={as_cap_gb*1024**3}'])
                        break
            print(f'[launch] shard {S}: pid {pid} → {log}'
                  + (f' (as={as_cap_gb}GB)' if as_cap_gb else ''))
            time.sleep(0.15)
        if max_concurrent and bi < len(batches) - 1:
            print(f'[launch] batch {bi+1}/{len(batches)} launched; waiting...')
            while any(not os.path.exists(f'{sh["dir"]}/tmp/{fam}/SYSTEMconfig')
                      and pid_alive(pids[S])
                      for S, sh in batch):
                time.sleep(30)
    meta['pids'] = pids; meta['launched'] = stamp()
    meta['host'] = host or subprocess.getoutput('hostname')
    meta['as_cap_gb'] = as_cap_gb; meta['kira_bin'] = kira_bin
    json.dump(meta, open(f'{out}/PGEN_META.json', 'w'), indent=1)
    return pids


def status(out):
    meta = json.load(open(f'{out}/PGEN_META.json'))
    fam = meta['fam']
    rows = []
    for S, sh in meta['shards'].items():
        d = sh['dir']
        cfg = f'{d}/tmp/{fam}/SYSTEMconfig'
        neq = int(open(cfg).read().strip()) if os.path.exists(cfg) else None
        nfiles = len([f for f in os.listdir(f'{d}/tmp/{fam}')
                      if f.startswith('SYSTEM_')]) if os.path.isdir(f'{d}/tmp/{fam}') else 0
        pid = meta.get('pids', {}).get(str(S)) or meta.get('pids', {}).get(S)
        alive = bool(pid) and pid_alive(pid)
        # progress from log
        prog = ''
        log = f'{d}/kira.log'
        if os.path.exists(log):
            for line in reversed(open(log).readlines()[-200:]):
                m = re.search(r'sector \d+ \((\d+) of (\d+)\)', line)
                if m: prog = f'{m.group(1)}/{m.group(2)}'; break
        rows.append((S, 'DONE' if neq is not None else ('RUN' if alive else 'DEAD'),
                     neq, nfiles, prog, pid))
    for r in rows:
        print(f'  shard {r[0]:>4}  {r[1]:<4}  eqs={r[2]!s:<8}  '
              f'files={r[3]:<4}  {r[4]:<10}  pid={r[5]}')
    done = all(r[1] == 'DONE' for r in rows)
    return rows, done


# ----------------------------------------------------------------- merge

def read_system_eqs(path):
    """→ list of (head_weight:int, body:bytes) per Eq block."""
    out = []
    with gzip.open(path, 'rt') as f:
        lines = f.read().split('\n')
    i = 0
    while i < len(lines):
        if lines[i] != 'Eq':
            i += 1; continue
        head = int(lines[i+1]); n = int(lines[i+2])
        body = '\n'.join(lines[i:i+3+n])
        out.append((head, body))
        i += 3 + n
    return out


def merge(out, target, top_from, dedupe):
    meta = json.load(open(f'{out}/PGEN_META.json'))
    fam, top = meta['fam'], meta['top']
    tdir = f'{target}/tmp/{fam}'
    os.makedirs(tdir, exist_ok=True)
    # gather: sector → list of source files
    by_sec = {}
    for S, sh in meta['shards'].items():
        sd = f'{sh["dir"]}/tmp/{fam}'
        if not os.path.isdir(sd): continue
        for f in os.listdir(sd):
            m = re.match(rf'SYSTEM_{fam}_(\d+)(?:_\d+)?\.gz$', f)
            if m:
                by_sec.setdefault(int(m.group(1)), []).append(f'{sd}/{f}')
    if top_from:
        src = f'{top_from}/tmp/{fam}/SYSTEM_{fam}_{top}.gz'
        if os.path.exists(src):
            by_sec.setdefault(top, []).append(src)
            print(f'[merge] injected top-sector SYSTEM from {src}')
        else:
            print(f'[merge] WARN --top-from {src} not found')
    total = 0
    n_dup = 0
    masters_union = set()
    for sec in sorted(by_sec):
        srcs = by_sec[sec]
        if len(srcs) == 1 and not dedupe:
            shutil.copy(srcs[0], f'{tdir}/SYSTEM_{fam}_{sec}.gz')
            with gzip.open(srcs[0], 'rt') as g:
                total += g.read().count('Eq\n')
            continue
        # union + dedup by full body (identical eqs from overlapping subtrees)
        seen = {}
        for src in srcs:
            for head, body in read_system_eqs(src):
                h = hashlib.md5(body.encode()).digest()
                if h in seen: n_dup += 1; continue
                seen[h] = body
        with gzip.open(f'{tdir}/SYSTEM_{fam}_{sec}.gz', 'wt') as g:
            for body in seen.values():
                g.write(body + '\n')
        total += len(seen)
    with open(f'{tdir}/SYSTEMconfig', 'w') as f:
        f.write(f'{total}\n')
    # masters: union of shard masters
    for S, sh in meta['shards'].items():
        mf = f'{sh["dir"]}/tmp/{fam}/masters'
        if os.path.exists(mf):
            masters_union.update(l.rstrip('\n') for l in open(mf) if l.strip())
    with open(f'{tdir}/masters', 'w') as f:
        f.write('\n'.join(sorted(masters_union)) + '\n')
    # results/masters (some kira paths read this)
    os.makedirs(f'{target}/results/{fam}', exist_ok=True)
    shutil.copy(f'{tdir}/masters', f'{target}/results/{fam}/masters')
    print(f'[merge] {len(by_sec)} sectors, {total} eqs '
          f'({n_dup} duplicates dropped), {len(masters_union)} masters → {tdir}')
    if top not in by_sec:
        print(f'[merge] WARN top sector {top} has NO SYSTEM file — '
              f'supply via --top-from (an oracle run, or a top-emit dir) '
              f'or --include-top')
    meta['merge'] = {'target': os.path.abspath(target), 'n_sectors': len(by_sec),
                     'n_eqs': total, 'n_dup': n_dup, 'has_top': top in by_sec,
                     'merged': stamp()}
    json.dump(meta, open(f'{out}/PGEN_META.json', 'w'), indent=1)
    return meta['merge']


# ----------------------------------------------------------------- top-emit

def emit_top(out, db=None, out_dir=None):
    """Option (c): emit the top sector's SYSTEM_<fam>_<TOP>.gz directly.

    Route-A raw generation (the IBP/LI operator templates in
    ref/sectormappings applied to every top-sector seed in the (r, s) box —
    the machinery of parallel_kira_gen_uds) plus the pyred_weight codec for
    head/term weights (the GPL tool from kira/bootloops-tools,
    loaded through kira_gpl_tools), written in the engine's SYSTEM block
    format (kira.cpp writer):

        Eq
        <head weight>
        <n terms>
        <n> <coeff> <weight> <sector> <topology> <flag2>   x n, descending
                                                           weight, head first

    Weight bits + integral ordering come from a completed run's
    results/kira.db (--db RUNDIR, else every shard that has one — all must
    agree and match the setup io); weight IDs are only defined relative to
    that run configuration. Preferred-master custom weights (pm mode) are
    applied in file order, trivial sectors skipped, exactly as the engine
    assigns them. Single-family runs assumed (topology id 0). Coefficients
    use the operator-template grammar of the uds route (monomials in a_i,
    d, m2 with rational constants).

    Equations whose head carries a positive ISP power land in a sector
    outside the family tree; they are emitted too (the codec ranks
    out-of-tree sectors after in-tree ones of the same line count, as
    sector_weight_table does) and counted in the receipt. The whole file is
    consumed by `merge --top-from <out>/top_emit`.
    """
    meta = json.load(open(f'{out}/PGEN_META.json'))
    fam, top, r, s, io = meta['fam'], meta['top'], meta['r'], meta['s'], meta['io']
    ref = meta['ref']
    uds = _load_sibling('parallel_kira_gen_uds')
    pw = _load_sibling('kira_gpl_tools').load('pyred_weight')
    # weight bits: from a completed run's kira.db, fail-closed
    dbs = []
    if db:
        dbs.append(os.path.join(db, 'results', 'kira.db'))
        if not os.path.exists(dbs[0]):
            sys.exit(f'FATAL: --db {db} has no results/kira.db')
    else:
        for S, sh in meta['shards'].items():
            p = os.path.join(sh['dir'], 'results', 'kira.db')
            if os.path.exists(p):
                dbs.append(p)
    if not dbs:
        sys.exit('FATAL: top-emit needs WEIGHTBITS from a completed run '
                 '(results/kira.db in a shard, or --db RUNDIR) — weight IDs '
                 'are only defined relative to the run configuration')
    got = {p: pw.from_kiradb(p) for p in dbs}
    vals = set(got.values())
    if len(vals) != 1:
        sys.exit(f'FATAL: kira.db files disagree on WEIGHTBITS/ordering: {got}')
    ((bits, iord),) = vals
    if iord != io:
        sys.exit(f'FATAL: kira.db integral_ordering {iord} != setup io {io} — '
                 f'weights would not match the merged solve')
    # route-A machinery: np from the operator template, ops, trivial sectors
    with open(f'{ref}/sectormappings/{fam}/IBP') as fh:
        first = fh.readline()
    np_ = len(uds.FAM_RE.match(first).group(2).split(','))
    codec = pw.WeightCodec(np_, bits, iord, topsectors=[top])
    ops = uds.load_ops(ref, fam)
    tp = f'{ref}/sectormappings/{fam}/trivialsector'
    triv = set(int(x) for x in open(tp).read().replace('\n', ',').split(',')
               if x.strip()) if os.path.exists(tp) else set()
    # preferred-master custom weights (file order, trivial-sector entries
    # skipped) — only when the shards generated WITH preferred_masters
    custom = {}
    if not meta.get('no_pm'):
        k = 0
        for nu in read_names(f'{ref}/preferred_masters', fam):
            if pw.sec_of(nu) in triv:
                continue
            k += 1
            custom[nu] = k
    wof = lambda nu: custom.get(nu) or codec.encode(nu)
    from fractions import Fraction as F
    tdir = out_dir or f'{out}/top_emit'
    fdir = f'{tdir}/tmp/{fam}'
    os.makedirs(fdir, exist_ok=True)
    path = f'{fdir}/SYSTEM_{fam}_{top}.gz'
    n_eq = n_seed = heads_out = n_terms = 0
    t0 = time.time()
    with gzip.open(path, 'wt') as fz:
        for seed in uds.seeds_for_sector(top, np_, r, s):
            n_seed += 1
            if n_seed % 2000 == 0:
                print(f'[top-emit] {n_seed} seeds, {n_eq} eqs, '
                      f'{time.time()-t0:.1f}s', flush=True)
            for op in ops:
                # accumulate coefficient by (target, d-power, m2-power) —
                # the gen_sector algorithm of parallel_kira_gen_uds
                acc = {}
                for shift, cterms in op:
                    for const, dpow, m2pow, aidx in cterms:
                        c = const if aidx is None else const * seed[aidx]
                        if c == 0:
                            continue
                        tgt = tuple(seed[i] + shift[i] for i in range(np_))
                        key = (tgt, dpow, m2pow)
                        acc[key] = acc.get(key, F(0)) + c
                by_tgt = {}
                for (tgt, dpow, m2pow), c in acc.items():
                    if c == 0:
                        continue
                    if pw.sec_of(tgt) in triv:
                        continue   # trivial-sector targets are exact zeros
                    by_tgt.setdefault(tgt, []).append((c, dpow, m2pow))
                # single-term relations stay: c(d,m2)*I = 0 pins I to zero,
                # and the engine writer emits any block with >= 1 nonzero term
                if not by_tgt:
                    continue
                terms = sorted(((wof(t), t, mons) for t, mons in by_tgt.items()),
                               key=lambda x: -x[0])
                n = len(terms)
                if pw.sec_of(terms[0][1]) & ~top:
                    heads_out += 1
                fz.write(f'Eq\n{terms[0][0]}\n{n}\n')
                for w, tgt, mons in terms:
                    coeff = ''.join(('+' if m[0] > 0 and i > 0 else '')
                                    + uds._fmt(*m) for i, m in enumerate(mons))
                    fz.write(f'{n} {coeff} {w} {pw.sec_of(tgt)} 0 0\n')
                n_eq += 1
                n_terms += n
    dt = time.time() - t0
    receipt = {'out': os.path.abspath(out), 'ref': ref, 'fam': fam, 'top': top,
               'r': r, 's': s, 'io': io, 'np': np_, 'weight_bits': list(bits),
               'db_files': dbs, 'n_ops': len(ops), 'n_seeds': n_seed,
               'n_eqs': n_eq, 'n_terms': n_terms,
               'heads_outside_tree': heads_out,
               'n_custom_weights': len(custom), 'wall_s': round(dt, 2),
               'file': os.path.abspath(path), 'created': stamp()}
    json.dump(receipt, open(f'{tdir}/TOP_EMIT_META.json', 'w'), indent=1)
    print(f'[top-emit] sector {top}: {n_seed} seeds x {len(ops)} ops -> '
          f'{n_eq} eqs ({n_terms} terms, {heads_out} out-of-tree heads) '
          f'in {dt:.1f}s -> {path}', flush=True)
    print(f'[top-emit] consume: merge {out} TARGET --top-from {tdir}',
          flush=True)
    return receipt


# ----------------------------------------------------------------- probe

def probe(ref, fam, top, r, s, tops, io):
    """setup a small probe, launch here, wait, verify weight-compatibility."""
    out = f'{ref}/../pgen_probe_{int(time.time())}'
    out = os.path.abspath(out)
    if tops is None:
        tops = cover_tops(top, bits(top))[:4]
    meta = setup(ref, out, fam, top, r, s, tops,
                 targets_file='de_targets_m2', io=io)
    t0 = time.time()
    launch(out, parallel=4, nice_lvl=10, host=None, env_extra=[])
    # wait
    while True:
        rows, done = status(out)
        if done: break
        if any(r[1] == 'DEAD' and r[2] is None for r in rows):
            print('[probe] shard died, abort'); break
        if time.time() - t0 > 3600:
            print('[probe] 1h timeout'); break
        time.sleep(15)
    dt = time.time() - t0
    # weight-compat verify: for each pair of shards sharing a sector, check
    # that identical head-weights map to identical bodies.
    verdict = {'shards': tops, 'r': r, 's': s, 'wall_s': round(dt, 1),
               'per_shard': {r[0]: {'eqs': r[2], 'files': r[3]} for r in rows}}
    compat_ok = True; checked = 0
    for i in range(len(tops)):
        for j in range(i+1, len(tops)):
            di = meta['shards'][tops[i]]['dir']; dj = meta['shards'][tops[j]]['dir']
            fdi = f'{di}/tmp/{fam}'; fdj = f'{dj}/tmp/{fam}'
            if not (os.path.isdir(fdi) and os.path.isdir(fdj)): continue
            common = set(os.listdir(fdi)) & set(os.listdir(fdj))
            common = [c for c in common if c.startswith('SYSTEM_') and c.endswith('.gz')]
            for c in common[:3]:
                ea = {h: b for h, b in read_system_eqs(f'{fdi}/{c}')}
                eb = {h: b for h, b in read_system_eqs(f'{fdj}/{c}')}
                for h in set(ea) & set(eb):
                    checked += 1
                    if ea[h] != eb[h]: compat_ok = False
    verdict['weight_compat'] = {'ok': compat_ok, 'eqs_checked': checked}
    verdict['probe_dir'] = out
    return verdict


# ----------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)

    p = sub.add_parser('setup')
    p.add_argument('ref'); p.add_argument('out')
    p.add_argument('--fam', required=True, help='family name'); p.add_argument('--top', type=int, default=1023)
    p.add_argument('--r', type=int, required=True); p.add_argument('--s', type=int, required=True)
    p.add_argument('--d', type=int, help='dot bound for the reduce box (omit = no d key, prior behavior)')
    p.add_argument('--tops'); p.add_argument('--targets', default='de_targets_m2')
    p.add_argument('--io', type=int, default=8)
    p.add_argument('--no-pm', action='store_true',
                   help='omit preferred_masters at gen (gated Route-B mode; '
                        'measured: pm at gen re-inflates seeds ~90x)')
    p.add_argument('--targets-mode', choices=['filter', 'own-sectors'], default='filter',
                   help='own-sectors: synthetic per-assigned-sector de_targets '
                        '(corner+max-dot+ISP reps at r,s bound)')
    p.add_argument('--manifest', help='JSON [{shard_id,tops:[..]},..] multi-top shards')
    p.add_argument('--nidx', type=int, default=15)

    p = sub.add_parser('launch')
    p.add_argument('out'); p.add_argument('--parallel', type=int, default=8)
    p.add_argument('--nice', type=int, default=5); p.add_argument('--host')
    p.add_argument('--env', action='append', default=[])
    p.add_argument('--kira', default='kira')
    p.add_argument('--as-cap-gb', type=int)
    p.add_argument('--max-concurrent', type=int)

    p = sub.add_parser('status'); p.add_argument('out')

    p = sub.add_parser('top-emit')
    p.add_argument('out')
    p.add_argument('--db', help='run dir whose results/kira.db supplies '
                                'WEIGHTBITS + integral ordering (default: '
                                'every shard that has one; all must agree)')
    p.add_argument('--out-dir', help='emit dir (default OUT/top_emit)')

    p = sub.add_parser('merge')
    p.add_argument('out'); p.add_argument('target')
    p.add_argument('--top-from'); p.add_argument('--dedupe', action='store_true', default=True)

    p = sub.add_parser('probe')
    p.add_argument('ref'); p.add_argument('--fam', required=True, help='family name')
    p.add_argument('--top', type=int, default=1023)
    p.add_argument('--r', type=int, default=11); p.add_argument('--s', type=int, default=0)
    p.add_argument('--tops'); p.add_argument('--io', type=int, default=8)

    a = ap.parse_args()
    tops = [int(x) for x in a.tops.split(',')] if getattr(a, 'tops', None) else None

    if a.cmd == 'setup':
        manifest = None
        if a.manifest:
            mj = json.load(open(a.manifest))
            manifest = mj['shards'] if isinstance(mj, dict) else mj
        setup(a.ref, a.out, a.fam, a.top, a.r, a.s, tops, a.targets, a.io,
              no_pm=a.no_pm, targets_mode=a.targets_mode, manifest=manifest,
              nidx=a.nidx, d=a.d)
    elif a.cmd == 'launch':
        launch(a.out, a.parallel, a.nice, a.host, a.env, a.kira,
               a.as_cap_gb, a.max_concurrent)
    elif a.cmd == 'status':
        status(a.out)
    elif a.cmd == 'top-emit':
        emit_top(a.out, a.db, a.out_dir)
    elif a.cmd == 'merge':
        merge(a.out, a.target, a.top_from, a.dedupe)
    elif a.cmd == 'probe':
        v = probe(a.ref, a.fam, a.top, a.r, a.s, tops, a.io)
        print(json.dumps(v, indent=1))


if __name__ == '__main__':
    main()
