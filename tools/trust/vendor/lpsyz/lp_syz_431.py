#!/usr/bin/env python3
# ============================================================================
# TRUST vendored engine: lp_syz_431.py (sha pinned in trust/_pins.py).
# ============================================================================
"""lp_syz_431.py — timed pilot: lp_syz port to LBL3Q sector 431.

Port of the probe2_syzygy lp_syz.py
(validated lbl3mopp certifier, sec39-class detector) to the LBL3Q family.

Family: LBL3Q (3-loop QED LbL non-planar), loops (k1,k2,k3), externals (p1,p2,p3),
kinematics p_i^2=0, p1.p2=s/2, p2.p3=t/2, p1.p3=(-s-t)/2; production slice
(s,t)=(-5,-3) (matches the production rowfarm cone gens), m2 numeric rational here.

Sector 431 = props {1,2,3,4,6,8,9} (7 vars). Blocker (off-basis dotted own-target
master, sec39/938 class): LBL3Q[1,1,1,3,0,1,0,1,1,0,0,0,0,0,0].
Question: does the LP syzygy tower yield a NEW relation (violation row) with
nonzero support on the blocker => provider row expressing it over other masters.

Stages:
  control  — sector-45 tower ({1,3,4,6}, K3-banana-like, 5 kira/TRUE masters):
             expect 0 violations at matched D (positive control of the port).
  probe    — smallest-first timing basis (timed pilot): P0 = 6-var child full syz;
             P1 = 7-var top syz with degBound cap; P2 = 7-var full syz.
             + tower node/iso census (no Singular) => cost projection.
  tower    — full sector-431 tower, forced masters = LP-representable cone
             masters (support ⊆ 431, all nu>=0), verdict on the blocker.

All arithmetic exact QQ at numeric rational (d0, m2, s, t). Checkpoint-first:
every stage writes its JSON as soon as it lands.
"""
import sys, os, time, json, re, itertools, subprocess, tempfile, resource
sys.set_int_max_str_digits(0)   # 7-var syz coefficients exceed 4300 digits
from fractions import Fraction
import sympy as sp
import flint

T0 = time.time()

# ---------------- family: LBL3Q ---------------------------------------------
# prop id -> (loop coeffs over (k1,k2,k3), ext coeffs over (p1,p2,p3), massive?)
PROPS = {
 1: ((1, 0, 0), (0, 0, 0), 1),   # k1^2 - m2
 2: ((1,-1, 0), (0, 0, 0), 1),   # (k1-k2)^2 - m2
 3: ((1,-1, 0), (0, 1, 0), 1),   # (k1-k2+p2)^2 - m2
 4: ((1,-1,-1), (0, 1, 0), 1),   # (k1-k2-k3+p2)^2 - m2
 5: ((1,-1,-1), (0, 1, 1), 1),   # (k1-k2-k3+p2+p3)^2 - m2
 6: ((1, 0,-1), (0, 1, 1), 1),   # (k1-k3+p2+p3)^2 - m2
 7: ((1, 0,-1), (-1, 0, 0), 1),  # (k1-k3-p1)^2 - m2
 8: ((1, 0, 0), (-1, 0, 0), 1),  # (k1-p1)^2 - m2
 9: ((0, 1, 0), (0, 0, 0), 0),   # k2^2
10: ((0, 0, 1), (0, 0, 0), 0),   # k3^2
}

def make_gram(s0, t0):
    s0 = sp.Rational(s0); t0 = sp.Rational(t0)
    return [[sp.Integer(0), s0/2, (-s0-t0)/2],
            [s0/2, sp.Integer(0), t0/2],
            [(-s0-t0)/2, t0/2, sp.Integer(0)]]

GRAM = None   # set in main from --s0/--t0

def build_G(node, m2v):
    """node = sorted tuple of prop ids. Returns (Gdict, iszero).
    Gdict: {exp_tuple: Fraction} over vars in node order. G = U + F."""
    n = len(node); ys = sp.symbols(f'y1:{n+1}')
    A = sp.zeros(3,3); Bp = [[sp.Integer(0)]*3 for _ in range(3)]; C = sp.Integer(0)
    for k,pid in enumerate(node):
        cL, cP, mm = PROPS[pid]; y = ys[k]
        for a in range(3):
            for b in range(3): A[a,b] += y*cL[a]*cL[b]
            for e in range(3): Bp[a][e] += y*cL[a]*cP[e]
        pp = sum(cP[e]*GRAM[e][f]*cP[f] for e in range(3) for f in range(3))
        C += y*(pp - (m2v if mm else 0))
    U = sp.expand(A.det())
    if U == 0: return None, True
    adj = A.adjugate()
    F = sp.Integer(0)
    for a in range(3):
        for b in range(3):
            BB = sum(Bp[a][e]*GRAM[e][f]*Bp[b][f] for e in range(3) for f in range(3))
            F += adj[a,b]*BB
    F = sp.expand(F - U*C)
    G = sp.expand(U + F)
    if G == 0: return None, True
    P = sp.Poly(G, *ys, domain='QQ')
    Gd = {tuple(e): Fraction(int(c.p), int(c.q)) for e,c in P.terms()}
    return Gd, False

# ---------------- dict-poly ops (verbatim from lp_syz.py) --------------------
def dmul_mono(p, e0, c0=Fraction(1)):
    return {tuple(a+b for a,b in zip(e,e0)): c*c0 for e,c in p.items()}
def dadd(p, q, cq=Fraction(1)):
    r = dict(p)
    for e,c in q.items():
        v = r.get(e, Fraction(0)) + cq*c
        if v: r[e] = v
        elif e in r: del r[e]
    return r
def dderiv(p, i):
    r = {}
    for e,c in p.items():
        if e[i]:
            e2 = list(e); e2[i] -= 1
            r[tuple(e2)] = c*e[i]
    return r
def drestrict(p, i):
    r = {}
    for e,c in p.items():
        if e[i]==0:
            r[tuple(v for k,v in enumerate(e) if k!=i)] = c
    return r

# ---------------- Singular syz (degBound knob added) --------------------------
def poly_str(p, n):
    if not p: return '0'
    ts = []
    for e,c in sorted(p.items()):
        mono = '*'.join(f'y{i+1}^{ei}' for i,ei in enumerate(e) if ei)
        cs = f'({c.numerator}/{c.denominator})' if c.denominator!=1 else f'({c.numerator})'
        ts.append(cs + ('*'+mono if mono else ''))
    return '+'.join(ts)

def singular_syz(Gd, n, timeout=300, degbound=0):
    gens = [poly_str(dderiv(Gd,i), n) for i in range(n)]
    gens.append(poly_str({e: -c for e,c in Gd.items()}, n))
    vec = ',\n  '.join(f'[{g}]' for g in gens)
    vlist = ','.join(f'y{i+1}' for i in range(n))
    db = f'degBound = {degbound};\n' if degbound else ''
    script = f"""ring R = 0,({vlist}),dp;
{db}module m1 = {vec};
module S = syz(m1);
"---BEGIN_SYZ---";
int i; int k;
for(i=1;i<=ncols(S);i++){{
  "---GEN---";
  for(k=1;k<=nrows(S);k++){{ string(S[k,i]); }}
}}
"---END_SYZ---";
exit;
"""
    with tempfile.NamedTemporaryFile('w', suffix='.sing', delete=False) as f:
        f.write(script); fn = f.name
    t = time.time()
    try:
        out = subprocess.run(['Singular','-q',fn], capture_output=True, text=True, timeout=timeout)
    finally:
        os.unlink(fn)
    dt = time.time()-t
    ms = re.search(r'---BEGIN_SYZ---\n(.*?)\n---END_SYZ---', out.stdout, re.S)
    if not ms: raise RuntimeError('singular fail: '+out.stdout[-500:]+out.stderr[-500:])
    syz = []
    for blk in ms.group(1).split('---GEN---'):
        lines = [l.strip() for l in blk.splitlines() if l.strip()]
        if len(lines) != n+1: continue
        gen = []
        for ln in lines:
            gen.append(parse_sing_poly(ln, n))
        syz.append((gen[:n], gen[n]))
    return syz, dt

def parse_sing_poly(s, n):
    if s == '0': return {}
    p = {}
    s = s.replace('-', '+-')
    for term in s.split('+'):
        term = term.strip()
        if not term: continue
        c = Fraction(1); e = [0]*n
        neg = term.startswith('-')
        if neg: term = term[1:]
        for fac in term.split('*'):
            if not fac: continue
            if fac[0] == 'y':
                if '^' in fac:
                    v,po = fac[1:].split('^'); e[int(v)-1] += int(po)
                else: e[int(fac[1:])-1] += 1
            else:
                c *= Fraction(fac)
        if neg: c = -c
        key = tuple(e); p[key] = p.get(key, Fraction(0)) + c
        if p[key] == 0: del p[key]
    return p

# ---------------- tower (verbatim core, degbound plumbed) ---------------------
def monos(n, dmax, _c={}):
    key = (n, dmax)
    if key not in _c:
        _c[key] = [e for e in itertools.product(range(dmax+1), repeat=n) if sum(e) <= dmax]
    return _c[key]

class Tower:
    def __init__(self, top, m2v, d0, D0, forced_masters, verbose=True,
                 Dtop=None, Dlow=None, degbound=0, syz_timeout=1200, census_only=False):
        self.top = tuple(sorted(top)); self.m2v = m2v; self.d0 = d0
        self.D0 = D0; self.Dtop = Dtop if Dtop is not None else D0
        self.Dlow = Dlow if Dlow is not None else D0
        self.degbound = degbound; self.syz_timeout = syz_timeout
        self.verbose = verbose
        self.stats = {}
        self._build_nodes()
        self._iso_classes()
        if census_only: return
        self._map_masters(forced_masters)
        self._columns()

    def node_D(self, R):
        n = len(R)
        if R == self.top: return self.Dtop
        if n <= 4: return self.Dlow
        return self.D0

    def _build_nodes(self):
        t = time.time(); self.G = {}; self.zero = set()
        for r in range(len(self.top), 0, -1):
            for S in itertools.combinations(self.top, r):
                Gd, isz = build_G(S, self.m2v)
                if isz: self.zero.add(S)
                else: self.G[S] = Gd
        self.stats['nodes_nonzero'] = len(self.G); self.stats['nodes_zero'] = len(self.zero)
        self.stats['t_symanzik'] = round(time.time()-t,2)
        if self.verbose:
            print(f"[tower] nodes: {len(self.G)} nonzero, {len(self.zero)} zero ({self.stats['t_symanzik']}s)", flush=True)

    @staticmethod
    def _iso_find(G1, G2, n):
        if len(G1) != len(G2): return None
        for sig in itertools.permutations(range(n)):
            ok = True
            for e,c in G1.items():
                e2 = [0]*n
                for i,v in enumerate(e): e2[sig[i]] = v
                if G2.get(tuple(e2)) != c: ok = False; break
            if ok: return sig
        return None

    def _iso_classes(self):
        t = time.time()
        self.rep = {}; self.perm = {}
        reps = []
        for S in sorted(self.G, key=lambda s:(-len(s), s)):
            n = len(S); found = False
            for R in reps:
                if len(R) != n: continue
                sig = self._iso_find(self.G[S], self.G[R], n)
                if sig is not None:
                    self.rep[S] = R; self.perm[S] = sig; found = True; break
            if not found:
                reps.append(S); self.rep[S] = S; self.perm[S] = tuple(range(n))
        self.reps = reps
        self.auts = {}
        for R in reps:
            n = len(R); auts = []
            for sig in itertools.permutations(range(n)):
                if sig == tuple(range(n)): continue
                ok = True
                for e,c in self.G[R].items():
                    e2 = [0]*n
                    for i,v in enumerate(e): e2[sig[i]] = v
                    if self.G[R].get(tuple(e2)) != c: ok = False; break
                if ok: auts.append(sig)
            self.auts[R] = auts
        self.stats['n_repclasses'] = len(reps)
        self.stats['repclass_sizes'] = {}
        for R in reps:
            self.stats['repclass_sizes'].setdefault(len(R), 0)
            self.stats['repclass_sizes'][len(R)] += 1
        self.stats['t_iso'] = round(time.time()-t,2)
        if self.verbose:
            print(f"[tower] {len(reps)} iso-classes of {len(self.G)} nodes; per-level {self.stats['repclass_sizes']} ({self.stats['t_iso']}s)", flush=True)

    def _aut_min(self, R, m):
        best = tuple(m)
        for sig in self.auts[R]:
            m2 = [0]*len(m)
            for i,v in enumerate(m): m2[sig[i]] = v
            t2 = tuple(m2)
            if t2 < best: best = t2
        return best

    def canon(self, S, m):
        R = self.rep[S]; sig = self.perm[S]
        m2 = [0]*len(S)
        for i,v in enumerate(m): m2[sig[i]] = v
        return R, self._aut_min(R, m2)

    def _map_masters(self, masters_nu):
        self.master_cols = []; self.master_info = []
        for nu in masters_nu:
            K = tuple(i+1 for i,v in enumerate(nu) if v > 0)
            mK = tuple(nu[i-1]-1 for i in K)
            GK, isz = build_G(K, self.m2v)
            assert not isz, f"master {nu} zero?"
            hit = None
            if K in self.rep:
                hit = self.canon(K, mK)
            else:
                for R in self.reps:
                    if len(R) != len(K): continue
                    sig = self._iso_find(GK, self.G[R], len(K))
                    if sig is not None:
                        m2 = [0]*len(K)
                        for i,v in enumerate(mK): m2[sig[i]] = v
                        hit = (R, self._aut_min(R, tuple(m2))); break
            if hit is None: raise RuntimeError(f"master {nu}: no iso class in tower")
            self.master_cols.append(hit); self.master_info.append((tuple(nu), hit))
        ncanon = len(set(self.master_cols))
        if ncanon != len(self.master_cols):
            # aut-redundant masters (family_aut class) — report, do not die
            from collections import Counter
            dup = [k for k,v in Counter(self.master_cols).items() if v>1]
            self.stats['aut_redundant_masters'] = [str(k) for k in dup]
            print(f"[tower] WARNING: {len(self.master_cols)-ncanon} aut-redundant forced masters: {dup}", flush=True)
        if self.verbose:
            print(f"[tower] mapped {len(self.master_cols)} forced masters ({ncanon} canonical)", flush=True)

    def _columns(self):
        mset = set(self.master_cols)
        self.colid = {}; self.cols = []
        self.node_order = sorted(self.reps, key=lambda s:(-len(s), s))
        self.node_blocks = {}
        for R in self.node_order:
            n = len(R); DR = self.node_D(R)
            blk = []
            for m in sorted(monos(n, DR),
                key=lambda e:(-sum(e), tuple(-v for v in e))):
                if self._aut_min(R, m) != m: continue
                if (R,m) in mset: continue
                self.colid[(R,m)] = len(self.cols); self.cols.append((R,m)); blk.append((R,m))
            self.node_blocks[R] = blk
        self.first_master_col = len(self.cols)
        for c in self.master_cols:
            if c not in self.colid:
                self.colid[c] = len(self.cols); self.cols.append(c)
        self.stats['n_columns'] = len(self.cols)
        if self.verbose: print(f"[tower] {len(self.cols)} columns (masters from {self.first_master_col})", flush=True)

    def build_relations(self):
        t = time.time(); self.rel_by_node = {R: [] for R in self.node_order}
        self.syz = {}; tsing = 0.0; nrel = 0; ndrop = 0
        for R in self.node_order:
            n = len(R)
            db = self.degbound if n >= 7 else 0
            syz, dt = singular_syz(self.G[R], n, timeout=self.syz_timeout, degbound=db)
            self.syz[R] = syz; tsing += dt
            if self.verbose:
                print(f"[syz] {R}: {len(syz)} gens {dt:.1f}s (degbound={db})", flush=True)
            for (alist, b) in syz:
                da = max((max((sum(e) for e in a), default=0) for a in alist), default=0)
                db2 = max((sum(e) for e in b), default=0)
                drel = max(da-1, db2, 0)
                if da == 0 and db2 == 0: continue
                DR = self.node_D(R)
                for m in monos(n, max(DR - drel, 0)):
                    row = {}; ok = True
                    Rm = {}
                    for i in range(n):
                        if alist[i]: Rm = dadd(Rm, dderiv(dmul_mono(alist[i], m), i))
                    if b: Rm = dadd(Rm, dmul_mono(b, m), cq=-self.d0/2)
                    for e,c in Rm.items():
                        key = (R, self._aut_min(R, e))
                        if key not in self.colid: ok = False; break
                        row[self.colid[key]] = row.get(self.colid[key], Fraction(0)) + c
                    if not ok: ndrop += 1; continue
                    for i in range(n):
                        if m[i] != 0 or not alist[i]: continue
                        bp = drestrict(dmul_mono(alist[i], m), i)
                        if not bp: continue
                        Sc = tuple(v for k,v in enumerate(R) if k != i)
                        if Sc in self.zero: continue
                        for e,c in bp.items():
                            key = self.canon(Sc, e)
                            if key not in self.colid: ok = False; break
                            row[self.colid[key]] = row.get(self.colid[key], Fraction(0)) + c
                        if not ok: break
                    if not ok: ndrop += 1; continue
                    row = {k:v for k,v in row.items() if v != 0}
                    if row:
                        self.rel_by_node[R].append(row); nrel += 1
        self.stats.update(t_singular=round(tsing,2), n_relations=nrel, n_dropped=ndrop,
                          t_build=round(time.time()-t,2))
        if self.verbose:
            print(f"[rel] {nrel} relations ({ndrop} dropped out-of-pool) "
                  f"singular {tsing:.1f}s build {self.stats['t_build']}s", flush=True)

    def solve(self):
        t = time.time(); self.rules = {}; self.violations = []; self.viol_rows = []
        pending = {R: [] for R in self.node_order}
        peak_local = (0,0)
        for R in self.node_order:
            rows = self.rel_by_node[R] + pending[R]
            if not rows: continue
            present = sorted({c for r in rows for c in r})
            pidx = {c:i for i,c in enumerate(present)}
            M = flint.fmpq_mat(len(rows), len(present))
            for ri,r in enumerate(rows):
                for c,v in r.items(): M[ri, pidx[c]] = flint.fmpq(v.numerator, v.denominator)
            if len(rows)*len(present) > peak_local[0]*peak_local[1]: peak_local = (len(rows), len(present))
            RR, rank = M.rref()
            blkset = set(self.colid[c] for c in self.node_blocks[R] if self.colid.get(c) is not None)
            for ri in range(rank):
                pj = None
                for j in range(len(present)):
                    if RR[ri, j] != 0: pj = j; break
                if pj is None: continue
                pcol = present[pj]
                tail = []
                pv = Fraction(int(RR[ri,pj].p), int(RR[ri,pj].q))
                for j in range(pj+1, len(present)):
                    v = RR[ri, j]
                    if v != 0:
                        tail.append((present[j], -Fraction(int(v.p), int(v.q))/pv))
                if pcol in blkset:
                    self.rules[pcol] = tail
                elif pcol >= self.first_master_col:
                    self.violations.append((R, pcol, len(tail)))
                    vr = {pcol: Fraction(1)}
                    for c,v in tail: vr[c] = -v
                    self.viol_rows.append(vr)
                else:
                    owner = self.cols[pcol][0]
                    row2 = {pcol: Fraction(1)}
                    for c,v in tail: row2[c] = -v
                    pending[owner].append(row2)
            if self.verbose:
                nfree_blk = sum(1 for c in self.node_blocks[R] if self.colid[c] not in self.rules)
                print(f"[solve] {R}: {len(rows)}x{len(present)} rank {rank}; unpivoted in-node {nfree_blk}", flush=True)
        self.stats['t_solve'] = round(time.time()-t,2)
        self.stats['peak_local_mat'] = peak_local
        self.stats['n_violations'] = len(self.violations)
        free_nonmaster = [self.cols[i] for i in range(self.first_master_col) if i not in self.rules]
        self.stats['n_free_nonmaster'] = len(free_nonmaster)
        self.free_nonmaster = free_nonmaster
        if self.verbose:
            print(f"[solve] done {self.stats['t_solve']}s; free non-master {len(free_nonmaster)}; violations {len(self.violations)}", flush=True)

    def conv(self, col, sigma):
        R, m = self.cols[col] if isinstance(col, int) else col
        w = sum(m) + len(R)
        c = Fraction(1)
        for mi in m:
            for k in range(2, mi+1): c *= k
        for j in range(1, w+1): c /= (2*self.d0 - j)
        if sigma == -1 and w % 2: c = -c
        return c

# ---------------- masters parsing --------------------------------------------
def parse_masters_file(path, support):
    """LP-representable masters: all indices >=0, support subset of `support`."""
    keep = []; skipped = []
    supp = set(support)
    for ln in open(path):
        m = re.match(r'\s*LBL3Q\[([0-9,\-\s]+)\]', ln)
        if not m: continue
        nu = tuple(int(x) for x in m.group(1).split(','))
        pos = {i+1 for i,v in enumerate(nu) if v > 0}
        if any(v < 0 for v in nu) or not pos <= supp:
            skipped.append(nu); continue
        keep.append(nu)
    return keep, skipped

def rss_gb():
    self_r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    ch = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    if sys.platform == 'darwin':  # ru_maxrss: bytes on macOS, KiB on Linux
        self_r /= 1024; ch /= 1024
    return round(self_r/1048576,2), round(ch/1048576,2)

def ckpt(obj, path):
    with open(path+'.tmp','w') as f: json.dump(obj, f, indent=1, default=str)
    os.replace(path+'.tmp', path)

BLOCKER = (1,1,1,3,0,1,0,1,1,0,0,0,0,0,0)
SEC431 = (1,2,3,4,6,8,9)

# ---------------- main ---------------------------------------------------------
if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--stage', required=True, choices=['control','control452','probe','tower'])
    ap.add_argument('--D0', type=int, default=4)
    ap.add_argument('--Dlow', type=int, default=None)
    ap.add_argument('--Dtop', type=int, default=None)
    ap.add_argument('--d0', default='97/23')
    ap.add_argument('--m2', default='5/7')
    ap.add_argument('--s0', default='-5')
    ap.add_argument('--t0', default='-3')
    ap.add_argument('--degbound', type=int, default=0)
    ap.add_argument('--syz-timeout', type=int, default=900)
    ap.add_argument('--masters-file', default='cone447_masters.txt')
    ap.add_argument('--out', default=None)
    args = ap.parse_args()
    d0 = Fraction(args.d0); m2f = Fraction(args.m2)
    GRAM = make_gram(Fraction(args.s0), Fraction(args.t0))
    globals()['GRAM'] = GRAM
    m2v = sp.Rational(m2f.numerator, m2f.denominator)
    res = {'stage': args.stage, 'd0': str(d0), 'm2': str(m2f),
           's0': args.s0, 't0': args.t0, 'D0': args.D0, 'degbound': args.degbound,
           'ts_start': time.strftime('%Y-%m-%d %H:%M:%S')}
    outp = args.out or f'{args.stage}.json'

    if args.stage == 'probe':
        # P0: 6-var child (drop prop 9 -> massive hexagon-like node), full syz
        node6 = (1,2,3,4,6,8)
        t = time.time(); Gd6, isz = build_G(node6, m2v)
        res['P0_node'] = str(node6); res['P0_G_terms'] = len(Gd6)
        res['P0_t_symanzik'] = round(time.time()-t,2)
        syz6, dt6 = singular_syz(Gd6, 6, timeout=args.syz_timeout)
        degs6 = [max(max((sum(e) for a in al for e in a), default=0),
                     max((sum(e) for e in b), default=0)) for al,b in syz6]
        res.update(P0_wall_s=round(dt6,2), P0_n_syz=len(syz6),
                   P0_syz_degs=sorted(set(degs6)))
        res['rss_gb_self_children'] = rss_gb()
        ckpt(res, outp); print(f"[probe] P0 done {dt6:.1f}s n_syz={len(syz6)}", flush=True)

        # census of the full 431 tower (no Singular) for projection
        tw = Tower(SEC431, m2v, d0, args.D0, [], verbose=True, census_only=True)
        res['census'] = {k: tw.stats[k] for k in
                         ('nodes_nonzero','nodes_zero','n_repclasses','repclass_sizes','t_symanzik','t_iso')}
        ckpt(res, outp)

        # P1: 7-var top node, degBound-capped syz
        t = time.time(); Gd7, isz = build_G(SEC431, m2v)
        res['P1_G_terms'] = len(Gd7); res['P1_t_symanzik'] = round(time.time()-t,2)
        db = args.degbound or 6
        try:
            syz7, dt7 = singular_syz(Gd7, 7, timeout=args.syz_timeout, degbound=db)
            degs7 = [max(max((sum(e) for a in al for e in a), default=0),
                         max((sum(e) for e in b), default=0)) for al,b in syz7]
            res.update(P1_wall_s=round(dt7,2), P1_n_syz=len(syz7), P1_degbound=db,
                       P1_syz_degs=sorted(set(degs7)))
        except subprocess.TimeoutExpired:
            res.update(P1_wall_s=None, P1_degbound=db, P1_verdict='TIMEOUT@%ds' % args.syz_timeout)
        res['rss_gb_self_children'] = rss_gb()
        ckpt(res, outp); print(f"[probe] P1 done: {res.get('P1_wall_s')}s", flush=True)

        # P2: 7-var full syz (no degBound)
        try:
            syz7f, dt7f = singular_syz(Gd7, 7, timeout=args.syz_timeout)
            degs7f = [max(max((sum(e) for a in al for e in a), default=0),
                          max((sum(e) for e in b), default=0)) for al,b in syz7f]
            res.update(P2_wall_s=round(dt7f,2), P2_n_syz=len(syz7f),
                       P2_syz_degs=sorted(set(degs7f)))
        except subprocess.TimeoutExpired:
            res.update(P2_wall_s=None, P2_verdict='TIMEOUT@%ds' % args.syz_timeout)
        res['rss_gb_self_children'] = rss_gb()
        res['wall'] = round(time.time()-T0,1)
        ckpt(res, outp)
        print(json.dumps(res, indent=1, default=str), flush=True)
        sys.exit(0)

    if args.stage == 'control':
        # sector 45 = props {1,3,4,6}; kira/TRUE masters: corner + 4 dotted (5 total)
        # NOTE: this is a 4-line equal-mass banana at P^2=t -> sec39/938 CLASS:
        # a violation here may be a GENUINE kira rank-deficiency, not a port bug.
        top = (1,3,4,6)
        masters = [(1,0,1,1,0,1,0,0,0,0,0,0,0,0,0),
                   (1,0,1,2,0,1,0,0,0,0,0,0,0,0,0),
                   (1,0,1,3,0,1,0,0,0,0,0,0,0,0,0),
                   (1,0,1,2,0,2,0,0,0,0,0,0,0,0,0),
                   (1,0,1,4,0,1,0,0,0,0,0,0,0,0,0)]
    elif args.stage == 'control452':
        # sector 452 = props {3,7,8,9}; sector-restricted reduction gated 8/8 exact vs kira;
        # kira masters: corner + dot-3 (2 total). Clean port control: expect 0 violations.
        top = (3,7,8,9)
        masters = [(0,0,1,0,0,0,1,1,1,0,0,0,0,0,0),
                   (0,0,2,0,0,0,1,1,1,0,0,0,0,0,0)]
    else:
        top = SEC431
        masters, skipped = parse_masters_file(args.masters_file, SEC431)
        res['n_masters_forced'] = len(masters)
        res['n_masters_skipped_nonLP'] = len(skipped)
        res['skipped'] = [str(x) for x in skipped]
        assert BLOCKER in masters, 'blocker not in forced masters!'
        print(f"[main] {len(masters)} LP-representable forced masters; {len(skipped)} skipped (numerator/off-support)", flush=True)

    tw = Tower(top, m2v, d0, args.D0, masters, Dtop=args.Dtop, Dlow=args.Dlow,
               degbound=args.degbound, syz_timeout=args.syz_timeout)
    tw.build_relations()
    tw.solve()
    res.update(tw.stats)
    res['free_nonmaster_sample'] = [f"{c}" for c in tw.free_nonmaster[:20]]
    res['rss_gb_self_children'] = rss_gb()
    ckpt(res, outp)

    # ---- verdict on the blocker (tower stage) --------------------------------
    if args.stage == 'tower':
        hit2nu = {}
        for nu, hit in tw.master_info:
            hit2nu.setdefault(hit, []).append(nu)
        bcol = tw.colid[[h for nu,h in tw.master_info if nu==BLOCKER][0]]
        res['blocker_col'] = bcol
        touching = []
        for vr in tw.viol_rows:
            if bcol in vr:
                named = {}
                clean = True
                for c,v in vr.items():
                    cc = tw.cols[c]
                    if c >= tw.first_master_col or cc in hit2nu:
                        nm = 'M' + str(hit2nu.get(cc, [cc])[0])
                    else:
                        nm = 'FREE' + str(cc); clean = False
                    named[nm] = str(v)
                touching.append({'clean_master_only': clean, 'row': named})
        res['n_viol_rows'] = len(tw.viol_rows)
        res['viol_rows_touching_blocker'] = touching
        res['blocker_expressible'] = bool(touching) and any(t['clean_master_only'] for t in touching)
        # all violation rows named (they are new null-vectors regardless of blocker)
        allv = []
        for vr in tw.viol_rows[:40]:
            named = {}
            for c,v in vr.items():
                cc = tw.cols[c]
                nm = ('M'+str(hit2nu[cc][0])) if cc in hit2nu else ('X'+str(cc))
                named[nm] = str(v)
            allv.append(named)
        res['viol_rows_named'] = allv
    elif args.stage in ('control','control452'):
        res['control_pass'] = (len(tw.viol_rows) == 0)
        # name every violation row for audit
        hit2nu = {}
        for nu, hit in tw.master_info:
            hit2nu.setdefault(hit, []).append(nu)
        allv = []
        for vr in tw.viol_rows[:20]:
            named = {}
            for c,v in vr.items():
                cc = tw.cols[c]
                nm = ('M'+str(hit2nu[cc][0])) if cc in hit2nu else ('X'+str(cc))
                named[nm] = str(v)
            allv.append(named)
        res['viol_rows_named'] = allv

    res['rss_gb_self_children'] = rss_gb()
    res['wall'] = round(time.time()-T0,1)
    ckpt(res, outp)
    print(json.dumps({k:v for k,v in res.items() if k not in ('viol_rows_named','skipped')}, indent=1, default=str), flush=True)
