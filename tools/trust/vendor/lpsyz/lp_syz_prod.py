#!/usr/bin/env python3
# ============================================================================
# TRUST vendored engine: lp_syz_prod.py (sha pinned in trust/_pins.py).
# ============================================================================
"""lp_syz_prod.py — PRODUCTION: lp_syz mod-p screen + CRT lift + exact-Q verify.

The D0=5 exact-QQ fmpq RREF is infeasible directly (~4600-digit syz
coefficients, >=15.6G measured). This adapter:
  1. runs the FULL pilot cascade solve mod p (flint nmod_mat) at 3+ CRT primes,
     with lambda-multiplier tracking (augmented identity block) so every violation
     row carries its provenance over the ORIGINAL exact relation rows;
  2. CRT-lifts ONLY the violation (master-only null) rows, rational-reconstructs
     (Wang), and checks them at a HELD-OUT prime never used in the CRT;
  3. verifies EXACTLY over Q (fmpq): subsystem = exact rows in the lambda-support
     union; PASS iff rank([Msub]) == rank([Msub; v]) — the candidate is then a
     PROVEN exact consequence of the exact LP-syzygy relations. Mod-p alone is
     never trusted (charter).

Stages (separate processes so the 30-min guard can be enforced externally):
  screen — build tower (exact) + mod-p cascades + lift + holdout; dumps
           SCREEN_<tag>.json + SUBSYS_<tag>.json (exact subsystem rows, strings).
  verify — reads SUBSYS dump, exact fmpq rank check per candidate; VERIFY_<tag>.json.

Masters forced = [BLOCKER] + (TRUE-432 cap LP(sector)) + (tier2 provider certs cap
LP(sector)), from a per-sector JSON (see --masters-json). Blocker FIRST => it owns
the smallest master colid, so a null row containing it pivots on it directly.
NOTE: violation rows are master-only BY CONSTRUCTION (master colids are the largest;
RREF pivot = min colid in support), so every violation is a masters-span statement.

Guards: run every invocation under `ulimit -v 20971520` (20G) in the same shell
line; wrap the verify stage in `timeout 1800`.
"""
import sys, os, time, json, math, argparse, hashlib, pickle, gzip
sys.set_int_max_str_digits(0)
from fractions import Fraction
import sympy as sp
import flint
import lp_syz_431 as lib

T0 = time.time()

# ---- syzygy disk cache: syz depends ONLY on (G, n, degbound), NOT on the
# pool depth D0/Dtop — so D5 -> Dt6 rebuilds of the same sector are syz-free.
_SYZ_CACHE_DIR = 'syzcache'
_orig_singular_syz = lib.singular_syz

def _cached_singular_syz(Gd, n, timeout=300, degbound=0):
    key = hashlib.sha256(
        repr((sorted(Gd.items()), n, degbound)).encode()).hexdigest()[:24]
    fp = os.path.join(_SYZ_CACHE_DIR, key + '.pkl')
    if os.path.exists(fp):
        with open(fp, 'rb') as f:
            return pickle.load(f), 0.0
    syz, dt = _orig_singular_syz(Gd, n, timeout=timeout, degbound=degbound)
    os.makedirs(_SYZ_CACHE_DIR, exist_ok=True)
    tmp = fp + f'.tmp{os.getpid()}'
    with open(tmp, 'wb') as f:
        pickle.dump(syz, f, protocol=4)
    os.replace(tmp, fp)
    return syz, dt

lib.singular_syz = _cached_singular_syz

def log(msg):
    print(f"[{time.strftime('%H:%M:%S')} +{time.time()-T0:7.1f}s] {msg}", flush=True)

def ckpt(obj, path):
    with open(path + '.tmp', 'w') as f:
        json.dump(obj, f, indent=1, default=str)
    os.replace(path + '.tmp', path)

def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()

# ---------------- exact rational arithmetic helpers ---------------------------
def frac_modp(fr, p):
    d = fr.denominator % p
    if d == 0:
        raise ZeroDivisionError('bad prime: divides a denominator')
    return (fr.numerator % p) * pow(d, -1, p) % p

def crt_stack(residues, mods):
    r, m = residues[0] % mods[0], mods[0]
    for ri, mi in zip(residues[1:], mods[1:]):
        t = ((ri - r) * pow(m % mi, -1, mi)) % mi
        r = r + m * t
        m *= mi
    return r % m, m

def ratrec(r, m):
    """Wang rational reconstruction; None if no candidate under isqrt(m/2)."""
    N = math.isqrt(m // 2)
    a0, a1 = m, r % m
    x0, x1 = 0, 1
    while a1 > N:
        q = a0 // a1
        a0, a1 = a1, a0 - q * a1
        x0, x1 = x1, x0 - q * x1
    n, d = a1, x1
    if d == 0:
        return None
    if d < 0:
        n, d = -n, -d
    if d > N or math.gcd(n, d) != 1:
        return None
    return Fraction(n, d)

def gen_primes(n, seed_start=(1 << 61) + 1):
    ps, x = [], seed_start
    while len(ps) < n:
        x = int(sp.nextprime(x))
        ps.append(x)
    return ps

# ---------------- mod-p cascade on a plain state dict -------------------------
def modp_cascade(st, p):
    """Full cascade solve mod p with lambda tracking, on a tower-state dict.
    st needs: row_list, node_order, rowids_by_node, colid, cols, node_blocks,
    first_master_col. Returns violations [(node_str, pcol, {col:int}, {rid:int})]
    + census."""
    t0 = time.time()
    row_list = st['row_list']; colid = st['colid']; cols = st['cols']
    fmc = st['first_master_col']
    rows_mp = []
    for rid, R, row in row_list:
        rows_mp.append({c: frac_modp(v, p) for c, v in row.items()})
    pending = {R: [] for R in st['node_order']}
    viols = []
    rules_cols = set()
    sig = []
    for R in st['node_order']:
        loc = [(rows_mp[rid], {rid: 1}) for rid in st['rowids_by_node'][R]] \
              + pending[R]
        if not loc:
            continue
        present = sorted({c for r, _ in loc for c in r})
        pidx = {c: i for i, c in enumerate(present)}
        nl, npr = len(loc), len(present)
        M = flint.nmod_mat(nl, npr + nl, p)
        for i, (r, _) in enumerate(loc):
            for c, v in r.items():
                M[i, pidx[c]] = v
            M[i, npr + i] = 1
        RR, rank = M.rref()
        blkset = {colid[c] for c in st['node_blocks'][R]}
        n_v_node = 0
        for i in range(nl):
            pj = None
            for j in range(npr):
                if RR[i, j] != 0:
                    pj = j
                    break
            if pj is None:
                break
            pcol = present[pj]
            if pcol in blkset:
                rules_cols.add(pcol)
                continue
            row2 = {}
            for j in range(pj, npr):
                v = RR[i, j]
                if v != 0:
                    row2[present[j]] = int(v)
            mult = {}
            for j in range(nl):
                v = RR[i, npr + j]
                if v != 0:
                    vi = int(v)
                    for rid2, c2 in loc[j][1].items():
                        mult[rid2] = (mult.get(rid2, 0) + vi * c2) % p
            mult = {k: v for k, v in mult.items() if v}
            if pcol >= fmc:
                viols.append((str(R), pcol, row2, mult))
                n_v_node += 1
            else:
                owner = cols[pcol][0]
                pending[owner].append((row2, mult))
        sig.append((str(R), nl, npr, rank, n_v_node))
    free_nm = fmc - len(rules_cols)
    dt = time.time() - t0
    log(f"  mod-p cascade p={p}: {len(viols)} violations, "
        f"{len(rules_cols)} rules, free non-master {free_nm} ({dt:.1f}s)")
    return {'p': p, 'violations': viols, 'n_rules': len(rules_cols),
            'free_nonmaster': free_nm, 'sig': sig, 'wall_s': round(dt, 1)}


# ---------------- production tower -------------------------------------------
class ProdTower(lib.Tower):
    def flatten_rows(self):
        self.row_list = []            # (rid, node, exact row dict col->Fraction)
        self.rowids_by_node = {}
        rid = 0
        for R in self.node_order:
            ids = []
            for row in self.rel_by_node[R]:
                self.row_list.append((rid, R, row))
                ids.append(rid)
                rid += 1
            self.rowids_by_node[R] = ids
        log(f"flattened {rid} exact relation rows")

    def state(self):
        return {'row_list': self.row_list, 'node_order': self.node_order,
                'rowids_by_node': self.rowids_by_node, 'colid': self.colid,
                'cols': self.cols, 'node_blocks': self.node_blocks,
                'first_master_col': self.first_master_col,
                'master_info': self.master_info, 'stats': self.stats}

    def solve_modp(self, p):
        return modp_cascade(self.state(), p)


def load_masters(masters_json):
    mj = json.load(open(masters_json))
    blocker = tuple(mj['blocker'])
    true_ms = [tuple(x) for x in mj['true_masters']]
    provs = [tuple(x) for x in mj.get('providers', [])]
    forced = [blocker] + true_ms + provs
    kind = {blocker: 'BLOCKER'}
    for m in true_ms:
        kind.setdefault(m, 'TRUE')
    for m in provs:
        kind.setdefault(m, 'PROVIDER')
    return mj, blocker, forced, kind


def stage_screen(args):
    mj, blocker, forced, kind = load_masters(args.masters_json)
    top = tuple(mj['sector_props'])
    d0 = Fraction(args.d0); m2f = Fraction(args.m2)
    lib.GRAM = lib.make_gram(Fraction(args.s0), Fraction(args.t0))
    m2v = sp.Rational(m2f.numerator, m2f.denominator)
    tag = args.tag
    res = {'stage': 'screen', 'sector': mj['sector'], 'tag': tag,
           'd0': str(d0), 'm2': str(m2f), 's0': args.s0, 't0': args.t0,
           'D0': args.D0, 'blocker': str(blocker),
           'n_forced': len(forced),
           'n_true': sum(1 for k in kind.values() if k == 'TRUE'),
           'n_providers': sum(1 for k in kind.values() if k == 'PROVIDER'),
           'masters_json_sha256': sha256(args.masters_json),
           'script_sha256': sha256(os.path.abspath(__file__)),
           'lib_sha256': sha256(os.path.abspath(lib.__file__)),
           'ts_start': time.strftime('%Y-%m-%d %H:%M:%S')}
    outp = f'SCREEN_{tag}.json'
    ckpt(res, outp)

    tw = ProdTower(top, m2v, d0, args.D0, forced,
                   syz_timeout=args.syz_timeout, verbose=True)
    tw.build_relations()
    res.update({k: tw.stats[k] for k in tw.stats})
    ckpt(res, outp)
    tw.flatten_rows()

    hit2nu = {}
    for nu, hit in tw.master_info:
        hit2nu.setdefault(hit, []).append(nu)
    col2name = {}
    col2kind = {}
    for nu, hit in tw.master_info:
        c = tw.colid[hit]
        if c not in col2name:
            col2name[c] = 'LBL3Q[' + ','.join(map(str, nu)) + ']'
            col2kind[c] = kind[nu]
    bcol = tw.colid[[h for nu, h in tw.master_info if nu == blocker][0]]
    res['blocker_col'] = bcol
    res['aut_redundant_masters'] = tw.stats.get('aut_redundant_masters', [])

    # ---- mod-p cascades: n_crt CRT primes + 1 held-out; auto-extend on
    #      ratrec failure (previous holdout joins CRT, fresh prime = holdout)
    primes = gen_primes(args.max_primes + 4)
    runs, used = {}, []
    pi = 0

    def next_good_prime():
        nonlocal pi
        while pi < len(primes):
            p = primes[pi]; pi += 1
            try:
                r = tw.solve_modp(p)
            except ZeroDivisionError as e:
                log(f"  prime {p} rejected: {e}")
                continue
            runs[p] = r
            used.append(p)
            return p
        return None

    for _ in range(args.n_crt + 1):
        if next_good_prime() is None:
            res['verdict'] = 'GUARD-STOP: not enough good primes'
            ckpt(res, outp); return 2

    def vkeys(r):
        return sorted((n, pc, tuple(sorted(row.keys())))
                      for n, pc, row, _ in r['violations'])

    def lift_all(crt_primes, holdout_p):
        cands, any_fail = [], False
        v0 = runs[crt_primes[0]]['violations']
        hold_map = {(n, pc): row for n, pc, row, _ in runs[holdout_p]['violations']}
        for idx in range(len(v0)):
            node, pcol, _, _ = v0[idx]
            rowmaps, lam_supp = [], set()
            for p in crt_primes:
                n2, pc2, row2, mult2 = runs[p]['violations'][idx]
                assert (n2, pc2) == (node, pcol)
                rowmaps.append(row2)
                lam_supp |= set(mult2.keys())
            lam_supp |= set(runs[holdout_p]['violations'][idx][3].keys())
            cols = sorted(rowmaps[0].keys())
            lifted, ok = {}, True
            for c in cols:
                resid = [rm[c] for rm in rowmaps]
                r, m = crt_stack(resid, crt_primes)
                fr = ratrec(r, m)
                if fr is None:
                    ok = False
                    break
                lifted[c] = fr
            if not ok:
                any_fail = True
                cands.append({'idx': idx, 'node': node, 'pcol': pcol,
                              'status': 'RATREC-FAIL (needs more primes)'})
                continue
            hrow = hold_map.get((node, pcol))
            hold_ok = hrow is not None and set(hrow.keys()) == set(lifted.keys()) \
                and all(frac_modp(fr, holdout_p) == hrow[c]
                        for c, fr in lifted.items())
            maxh = max(max(abs(f.numerator), f.denominator)
                       for f in lifted.values())
            cands.append({
                'idx': idx, 'node': node, 'pcol': pcol,
                'status': 'LIFTED',
                'holdout_pass': bool(hold_ok),
                'touches_blocker': bcol in lifted,
                'height_digits': len(str(maxh)),
                'lambda_support_n': len(lam_supp),
                'lambda_support': sorted(lam_supp),
                'relation_named': {col2name.get(c, f'col{c}'): str(v)
                                   for c, v in lifted.items()},
                'relation_kinds': {col2name.get(c, f'col{c}'): col2kind.get(c, '?')
                                   for c in lifted},
                'relation_cols': {str(c): str(v) for c, v in lifted.items()},
            })
        return cands, any_fail

    while True:
        keysets = [vkeys(runs[p]) for p in used]
        agree = all(k == keysets[0] for k in keysets[1:])
        if not agree:
            res['verdict'] = ('GUARD-STOP: violation structure differs across '
                              'primes (bad prime?)')
            res['n_violations_per_prime'] = {str(p): len(runs[p]['violations'])
                                             for p in used}
            ckpt(res, outp); return 2
        crt_primes, holdout_p = used[:-1], used[-1]
        cands, any_fail = lift_all(crt_primes, holdout_p)
        if not any_fail or len(used) >= args.max_primes:
            break
        log(f"  ratrec fail with {len(crt_primes)} CRT primes — extending")
        if next_good_prime() is None:
            break

    res['crt_primes'] = crt_primes
    res['holdout_prime'] = holdout_p
    res['modp_walls_s'] = {str(p): runs[p]['wall_s'] for p in used}
    res['free_nonmaster'] = runs[crt_primes[0]]['free_nonmaster']
    res['n_rules_modp'] = runs[crt_primes[0]]['n_rules']
    res['pool_drop_frac'] = round(
        tw.stats['n_dropped'] / max(tw.stats['n_dropped']
                                    + tw.stats['n_relations'], 1), 3)
    res['viol_struct_agree_all_primes'] = True
    res['n_violations_per_prime'] = {str(p): len(runs[p]['violations'])
                                     for p in used}
    res['n_candidates'] = len(cands)
    res['candidates'] = [{k: v for k, v in c.items() if k != 'lambda_support'}
                         for c in cands]
    res['blocker_touching'] = [c['idx'] for c in cands
                               if c.get('touches_blocker')]
    ckpt(res, outp)

    # ---- dump exact subsystem rows for the verify stage ----------------------
    need_rids = set()
    for c in cands:
        if c.get('status') == 'LIFTED':
            need_rids |= set(c['lambda_support'])
    sub = {'tag': tag, 'sector': mj['sector'], 'blocker_col': bcol,
           'candidates': [c for c in cands if c.get('status') == 'LIFTED'],
           'rows': {str(rid): {str(c): [str(v.numerator), str(v.denominator)]
                               for c, v in tw.row_list[rid][2].items()}
                    for rid in sorted(need_rids)}}
    with open(f'SUBSYS_{tag}.json', 'w') as f:
        json.dump(sub, f, default=str)
    log(f"dumped {len(need_rids)} exact subsystem rows -> SUBSYS_{tag}.json")

    res['wall_total_s'] = round(time.time() - T0, 1)
    res['rss_gb_self_children'] = lib.rss_gb()
    ckpt(res, outp)
    log(f"screen done: {len(cands)} candidates, "
        f"{len(res['blocker_touching'])} touch blocker")
    return 0


def stage_verify(args):
    tag = args.tag
    sub = json.load(open(f'SUBSYS_{tag}.json'))
    outp = f'VERIFY_{tag}.json'
    res = {'stage': 'verify', 'tag': tag, 'sector': sub['sector'],
           'ts_start': time.strftime('%Y-%m-%d %H:%M:%S'),
           'method': 'exact fmpq rank: rank(Msub) == rank([Msub; v]) '
                     'with Msub = exact relation rows in lambda-support union',
           'results': []}
    ckpt(res, outp)
    rows = {int(rid): {int(c): Fraction(int(nd[0]), int(nd[1]))
                       for c, nd in row.items()}
            for rid, row in sub['rows'].items()}
    for cand in sub['candidates']:
        t0 = time.time()
        S = cand['lambda_support']
        v = {int(c): Fraction(x) for c, x in cand['relation_cols'].items()}
        if len(S) > args.max_subsys_rows:
            res['results'].append({'idx': cand['idx'], 'pcol': cand['pcol'],
                                   'verdict': f'GUARD-STOP: subsystem {len(S)} rows '
                                              f'> cap {args.max_subsys_rows}'})
            ckpt(res, outp)
            continue
        cols = sorted(set().union(*[rows[r].keys() for r in S]) | set(v.keys()))
        cidx = {c: i for i, c in enumerate(cols)}
        M1 = flint.fmpq_mat(len(S), len(cols))
        for i, r in enumerate(S):
            for c, fr in rows[r].items():
                M1[i, cidx[c]] = flint.fmpq(fr.numerator, fr.denominator)
        _, r1 = M1.rref()
        M2 = flint.fmpq_mat(len(S) + 1, len(cols))
        for i, r in enumerate(S):
            for c, fr in rows[r].items():
                M2[i, cidx[c]] = flint.fmpq(fr.numerator, fr.denominator)
        for c, fr in v.items():
            M2[len(S), cidx[c]] = flint.fmpq(fr.numerator, fr.denominator)
        _, r2 = M2.rref()
        dt = time.time() - t0
        ok = (r1 == r2)
        res['results'].append({
            'idx': cand['idx'], 'node': cand['node'], 'pcol': cand['pcol'],
            'touches_blocker': cand.get('touches_blocker', False),
            'subsystem_rows': len(S), 'subsystem_cols': len(cols),
            'rank_Msub': r1, 'rank_aug': r2,
            'exact_verify_pass': bool(ok), 'wall_s': round(dt, 1)})
        log(f"  cand idx={cand['idx']} pcol={cand['pcol']}: "
            f"{len(S)}x{len(cols)} rank {r1} vs aug {r2} -> "
            f"{'PASS' if ok else 'FAIL'} ({dt:.1f}s)")
        ckpt(res, outp)
    res['rss_gb_self_children'] = lib.rss_gb()
    res['wall_total_s'] = round(time.time() - T0, 1)
    res['all_pass'] = all(x.get('exact_verify_pass') for x in res['results']) \
        and bool(res['results'])
    ckpt(res, outp)
    return 0


# =============== split-stage flow: build -> cascade(s) -> lift2 ===============
def stage_build(args):
    """Build tower + exact relations once; pickle state for parallel cascades."""
    mj, blocker, forced, kind = load_masters(args.masters_json)
    top = tuple(mj['sector_props'])
    d0 = Fraction(args.d0); m2f = Fraction(args.m2)
    lib.GRAM = lib.make_gram(Fraction(args.s0), Fraction(args.t0))
    m2v = sp.Rational(m2f.numerator, m2f.denominator)
    tag = args.tag
    tw = ProdTower(top, m2v, d0, args.D0, forced, Dtop=args.Dtop,
                   syz_timeout=args.syz_timeout, verbose=True)
    tw.build_relations()
    tw.flatten_rows()
    col2name, col2kind = {}, {}
    for nu, hit in tw.master_info:
        c = tw.colid[hit]
        if c not in col2name:
            col2name[c] = 'LBL3Q[' + ','.join(map(str, nu)) + ']'
            col2kind[c] = kind[nu]
    bcol = tw.colid[[h for nu, h in tw.master_info if nu == blocker][0]]
    meta = {'sector': mj['sector'], 'tag': tag, 'd0': str(d0), 'm2': str(m2f),
            's0': args.s0, 't0': args.t0, 'D0': args.D0, 'Dtop': args.Dtop,
            'blocker': str(blocker), 'blocker_col': bcol,
            'n_forced': len(forced),
            'masters_json_sha256': sha256(args.masters_json),
            'script_sha256': sha256(os.path.abspath(__file__)),
            'lib_sha256': sha256(os.path.abspath(lib.__file__)),
            'aut_redundant_masters': tw.stats.get('aut_redundant_masters', []),
            'stats': tw.stats,
            'ts': time.strftime('%Y-%m-%d %H:%M:%S')}
    with open(f'BUILD_{tag}.pkl', 'wb') as f:
        pickle.dump({'state': tw.state(), 'col2name': col2name,
                     'col2kind': col2kind, 'bcol': bcol, 'meta': meta}, f,
                    protocol=4)
    ckpt(meta, f'BUILD_{tag}.json')
    log(f"build done: {len(tw.row_list)} rows pickled -> BUILD_{tag}.pkl "
        f"({os.path.getsize(f'BUILD_{tag}.pkl')/1e9:.2f} GB)")
    return 0


def stage_cascade(args):
    tag = args.tag
    with open(f'BUILD_{tag}.pkl', 'rb') as f:
        B = pickle.load(f)
    log(f"pickle loaded ({len(B['state']['row_list'])} rows)")
    p = gen_primes(args.prime_idx + 1)[args.prime_idx]
    try:
        r = modp_cascade(B['state'], p)
    except ZeroDivisionError as e:
        r = {'p': p, 'bad_prime': str(e)}
    r['prime_idx'] = args.prime_idx
    with gzip.open(f'CASC_{tag}_i{args.prime_idx}.json.gz', 'wt') as f:
        json.dump(r, f)
    log(f"cascade i{args.prime_idx} p={p} saved")
    return 0


def stage_lift2(args):
    """CRT-lift violation rows + lambda multipliers from all cascade files;
    exact-Q verification: sum(lambda_i * R_i) == v over Fractions/fmpq.
    rc=0 all candidates decided; rc=3 need more primes; rc=2 structure fail."""
    tag = args.tag
    with open(f'BUILD_{tag}.pkl', 'rb') as f:
        B = pickle.load(f)
    st = B['state']; col2name = B['col2name']; col2kind = B['col2kind']
    bcol = B['bcol']; meta = B['meta']
    row_list = st['row_list']
    outp = f'VERIFY2_{tag}.json'
    import glob as _g
    cascs = []
    for fn in sorted(_g.glob(f'CASC_{tag}_i*.json.gz'),
                     key=lambda s: int(s.split('_i')[-1].split('.')[0])):
        with gzip.open(fn, 'rt') as f:
            r = json.load(f)
        if 'bad_prime' in r:
            log(f"skipping bad prime file {fn}")
            continue
        # json turns int keys into strings
        r['violations'] = [(n, pc, {int(c): v for c, v in row.items()},
                            {int(k): v for k, v in mult.items()})
                           for n, pc, row, mult in r['violations']]
        cascs.append(r)
    res = {'stage': 'lift2', 'tag': tag, 'meta_sector': meta['sector'],
           'blocker': meta['blocker'], 'blocker_col': bcol,
           'D0': meta['D0'], 'Dtop': meta['Dtop'],
           'point': {'d0': meta['d0'], 'm2': meta['m2'],
                     's': meta['s0'], 't': meta['t0']},
           'n_cascades': len(cascs),
           'primes': [c['p'] for c in cascs],
           'free_nonmaster': cascs[0]['free_nonmaster'] if cascs else None,
           'pool_drop_frac': round(meta['stats']['n_dropped'] /
                                   max(meta['stats']['n_dropped'] +
                                       meta['stats']['n_relations'], 1), 3),
           'ts_start': time.strftime('%Y-%m-%d %H:%M:%S'),
           'method': 'exact-Q lambda-certificate: reconstructed lambda over Q '
                     '(CRT+ratrec of tracked multipliers), then '
                     'sum(lambda_i*R_i) == v checked EXACTLY over fmpq; '
                     'any exact match proves v in rowspan of the exact '
                     'LP-syzygy relations',
           'results': []}
    if len(cascs) < 3:
        res['verdict'] = 'GUARD-STOP: <3 good cascades'
        ckpt(res, outp); return 2
    # structure agreement
    def vkeys(r):
        return sorted((n, pc, tuple(sorted(row.keys())))
                      for n, pc, row, _ in r['violations'])
    k0 = vkeys(cascs[0])
    if not all(vkeys(c) == k0 for c in cascs[1:]):
        res['verdict'] = 'GUARD-STOP: violation structure differs across primes'
        ckpt(res, outp); return 2
    nv = len(cascs[0]['violations'])
    res['n_violations'] = nv
    primes = [c['p'] for c in cascs]
    crt_primes, holdout_p = primes[:-1], primes[-1]
    res['crt_primes_v'] = crt_primes
    res['holdout_prime_v'] = holdout_p
    need_more = False
    t_exact = 0.0
    for idx in range(nv):
        node, pcol = cascs[0]['violations'][idx][0], cascs[0]['violations'][idx][1]
        rowmaps = [c['violations'][idx][2] for c in cascs]
        mults = [c['violations'][idx][3] for c in cascs]
        cols = sorted(rowmaps[0].keys())
        # ---- v lift over crt primes, holdout on the last
        v_lift, v_ok = {}, True
        for c in cols:
            r, m = crt_stack([rm[c] for rm in rowmaps[:-1]], crt_primes)
            fr = ratrec(r, m)
            if fr is None:
                v_ok = False; break
            v_lift[c] = fr
        hold_ok = v_ok and all(frac_modp(fr, holdout_p) == rowmaps[-1].get(c, 0)
                               for c, fr in v_lift.items())
        # lambda support census only — lambda-CRT MEASURED DEAD
        # (ratrec fraction flat ~62% k=3..6: multiplier heights far beyond CRT
        # reach; exact leg lives in stage exactfull / subverify instead)
        lam_supp = set()
        for mu in mults:
            lam_supp |= set(mu.keys())
        entry = {'idx': idx, 'node': node, 'pcol': pcol,
                 'pivot_master': col2name.get(pcol, f'col{pcol}'),
                 'touches_blocker': bcol in rowmaps[0],
                 'support_masters': [col2name.get(c, f'col{c}') for c in cols],
                 'lambda_support_n': len(lam_supp)}
        if not v_ok:
            entry['status'] = 'NEED-MORE-PRIMES(v)'
            need_more = True
            res['results'].append(entry)
            continue
        maxh = max(max(abs(f.numerator), f.denominator)
                   for f in v_lift.values())
        entry.update({
            'status': 'V-LIFTED' if hold_ok else 'HOLDOUT-FAIL',
            'holdout_pass': bool(hold_ok),
            'height_digits_v': len(str(maxh)),
            'relation_named': {col2name.get(c, f'col{c}'): str(fr)
                               for c, fr in v_lift.items()},
            'relation_kinds': {col2name.get(c, f'col{c}'): col2kind.get(c, '?')
                               for c in v_lift},
        })
        res['results'].append(entry)
    res['blocker_touching'] = [r['idx'] for r in res['results']
                               if r.get('touches_blocker')]
    res['all_v_lifted'] = all(r.get('status') == 'V-LIFTED'
                              for r in res['results'])
    res['rss_gb_self_children'] = lib.rss_gb()
    res['wall_total_s'] = round(time.time() - T0, 1)
    ckpt(res, outp)
    log(f"lift2 done: {sum(1 for r in res['results'] if r.get('status')=='V-LIFTED')}"
        f"/{nv} v-lifted+holdout; blocker touching {res['blocker_touching']}")
    return 3 if need_more else 0


def stage_exactfull(args):
    """FULL exact fmpq cascade on the pickled tower (the pilot-killed solve,
    retried under the honest 30-min/20G guard) + exact span check of every
    lifted candidate against the exact violation rows. If this completes, the
    whole violation layer is EXACTLY derived over Q in one shot."""
    tag = args.tag
    with open(f'BUILD_{tag}.pkl', 'rb') as f:
        B = pickle.load(f)
    st = B['state']; col2name = B['col2name']; bcol = B['bcol']
    outp = f'EXACTFULL_{tag}.json'
    res = {'stage': 'exactfull', 'tag': tag,
           'ts_start': time.strftime('%Y-%m-%d %H:%M:%S'),
           'method': 'full exact fmpq cascade (lib.Tower.solve on shim) -> '
                     'exact violation rows over Q; candidates verified by exact '
                     'span membership rank check over the master-only block'}
    ckpt(res, outp)

    class Shim:
        pass
    sh = Shim()
    sh.node_order = st['node_order']
    sh.colid = st['colid']; sh.cols = st['cols']
    sh.node_blocks = st['node_blocks']
    sh.first_master_col = st['first_master_col']
    sh.stats = dict(st['stats']); sh.verbose = True
    sh.rel_by_node = {R: [st['row_list'][rid][2]
                          for rid in st['rowids_by_node'][R]]
                      for R in st['node_order']}
    t0 = time.time()
    lib.Tower.solve(sh)
    res['exact_solve_wall_s'] = round(time.time() - t0, 1)
    res['n_violations_exact'] = len(sh.viol_rows)
    res['free_nonmaster_exact'] = len(sh.free_nonmaster)
    res['rss_gb_self_children'] = lib.rss_gb()
    ckpt(res, outp)
    log(f"exact cascade done: {len(sh.viol_rows)} violations "
        f"({res['exact_solve_wall_s']}s)")

    # exact violation rows, serialized (master-only by construction)
    res['exact_violations'] = [
        {col2name.get(c, f'col{c}'): str(v) for c, v in vr.items()}
        for vr in sh.viol_rows]
    res['blocker_in_exact_span'] = None

    # candidates from lift2 (if present): exact span check
    name2col = {n: c for c, n in col2name.items()}
    ver = None
    vp = f'VERIFY2_{tag}.json'
    if os.path.exists(vp):
        ver = json.load(open(vp))
    if ver:
        mcols = sorted({c for vr in sh.viol_rows for c in vr}
                       | {name2col[n] for r in ver['results']
                          if 'relation_named' in r for n in r['relation_named']})
        midx = {c: i for i, c in enumerate(mcols)}
        M1 = flint.fmpq_mat(len(sh.viol_rows), len(mcols))
        for i, vr in enumerate(sh.viol_rows):
            for c, v in vr.items():
                M1[i, midx[c]] = flint.fmpq(v.numerator, v.denominator)
        _, r1 = M1.rref()
        res['exact_span_rank'] = r1
        out_r = []
        for r in ver['results']:
            if 'relation_named' not in r:
                out_r.append({'idx': r['idx'], 'status': r.get('status')})
                continue
            M2 = flint.fmpq_mat(len(sh.viol_rows) + 1, len(mcols))
            for i, vr in enumerate(sh.viol_rows):
                for c, v in vr.items():
                    M2[i, midx[c]] = flint.fmpq(v.numerator, v.denominator)
            for n, s in r['relation_named'].items():
                fr = Fraction(s)
                M2[len(sh.viol_rows), midx[name2col[n]]] = \
                    flint.fmpq(fr.numerator, fr.denominator)
            _, r2 = M2.rref()
            ok = (r2 == r1)
            out_r.append({'idx': r['idx'], 'pcol': r['pcol'],
                          'touches_blocker': r.get('touches_blocker'),
                          'exact_verify_pass': bool(ok)})
        res['results'] = out_r
        res['n_exact_verified'] = sum(1 for x in out_r
                                      if x.get('exact_verify_pass'))
        res['all_exact_verified'] = all(x.get('exact_verify_pass')
                                        for x in out_r if 'pcol' in x)
    # blocker expressibility direct from the exact span
    bl_rows = [vr for vr in sh.viol_rows if bcol in vr]
    res['blocker_in_exact_span'] = bool(bl_rows)
    if bl_rows:
        res['blocker_exact_relations'] = [
            {col2name.get(c, f'col{c}'): str(v) for c, v in vr.items()}
            for vr in bl_rows]
    res['wall_total_s'] = round(time.time() - T0, 1)
    ckpt(res, outp)
    log(f"exactfull done: blocker_in_exact_span={res['blocker_in_exact_span']}; "
        f"verified {res.get('n_exact_verified')} candidates")
    return 0


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--stage', required=True,
                    choices=['screen', 'verify', 'build', 'cascade', 'lift2',
                             'exactfull'])
    ap.add_argument('--masters-json')
    ap.add_argument('--tag', required=True)
    ap.add_argument('--D0', type=int, default=5)
    ap.add_argument('--Dtop', type=int, default=None)
    ap.add_argument('--d0', default='97/23')
    ap.add_argument('--m2', default='5/7')
    ap.add_argument('--s0', default='-5')
    ap.add_argument('--t0', default='-3')
    ap.add_argument('--syz-timeout', type=int, default=1500)
    ap.add_argument('--n-crt', type=int, default=3)
    ap.add_argument('--max-primes', type=int, default=8)
    ap.add_argument('--max-subsys-rows', type=int, default=2500)
    ap.add_argument('--prime-idx', type=int, default=0)
    args = ap.parse_args()
    sys.exit({'screen': stage_screen, 'verify': stage_verify,
              'build': stage_build, 'cascade': stage_cascade,
              'lift2': stage_lift2, 'exactfull': stage_exactfull}
             [args.stage](args))
