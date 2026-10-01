#!/usr/bin/env python3
# ============================================================================
# TRUST vendored engine: lp_syz.py (sha pinned in trust/_pins.py).
# ============================================================================
"""lp_syz.py — PROBE 2: E4C DE20-port syzygy-IBP adapted to lbl3mopp_k2sh.

Lee-Pomeransky per sector-tower node S: J_S[m] = int_{x>0} x^m G_S^{-d/2},
G_S = U_S + F_S. Singular syz (a,b): sum a_i dG/dx_i = b G  ==>
  J_S[ sum_i d_i(a_i x^m) - (d/2) b x^m ] = sum_{i: m_i=0} J_{S\\i}[(a_i x^m)|x_i=0]
Conversion: I[nu] with nu=m+1 on S: J = prod(m_i!) * (1/prod_{j=1..w}(2d-j)) * sigma^w * I,
w = |m|+|S| (up to a common Gamma(2d)Gamma(d/2) factor that cancels in ratios).
All arithmetic exact QQ at numeric rational (d0, eta0). Kinematics numeric
(s=-1, t=-1/3, msq=1), eta symbolic->numeric.
"""
import sys, os, time, json, re, itertools, subprocess, tempfile, resource
sys.set_int_max_str_digits(0)   # py3.12's 4300-digit cap crashes Fraction() on big syz coefficients
from fractions import Fraction
import sympy as sp
import flint

T0 = time.time()
TMAX = float(os.environ.get('LPSYZ_TMAX', '5400'))

# ---------------- family ----------------------------------------------------
# prop id -> (loop coeffs over (l,k1,k2), ext coeffs over (p1,p2,p3), massive?)
PROPS = {
 1: ((1,0,0),(1,0,0),1), 2: ((1,-1,0),(1,0,0),1), 3: ((0,1,-1),(0,0,0),1),
 4: ((0,0,1),(0,0,0),1), 5: ((0,0,1),(0,1,0),1), 6: ((0,0,1),(0,1,1),1),
 7: ((1,0,0),(1,1,1),1), 8: ((1,0,0),(0,0,0),1),
 9: ((0,1,0),(0,0,0),0),10: ((1,0,-1),(1,0,0),0),11: ((1,0,0),(0,1,0),0),
12: ((0,1,0),(1,0,0),0),13: ((0,1,0),(0,1,0),0),14: ((0,1,0),(0,0,1),0),
15: ((0,0,1),(1,0,0),0)}
GRAM = [[sp.Integer(0), sp.Rational(-1,2), sp.Rational(2,3)],
        [sp.Rational(-1,2), sp.Integer(0), sp.Rational(-1,6)],
        [sp.Rational(2,3), sp.Rational(-1,6), sp.Integer(0)]]

def build_G(node, eta0):
    """node = sorted tuple of prop ids. Returns (Gdict, iszero).
    Gdict: {exp_tuple: Fraction} over vars in node order."""
    n = len(node); ys = sp.symbols(f'y1:{n+1}')
    m2 = sp.Integer(1) + eta0
    A = sp.zeros(3,3); Bp = [[sp.Integer(0)]*3 for _ in range(3)]; C = sp.Integer(0)
    for k,pid in enumerate(node):
        cL, cP, mm = PROPS[pid]; y = ys[k]
        for a in range(3):
            for b in range(3): A[a,b] += y*cL[a]*cL[b]
            for e in range(3): Bp[a][e] += y*cL[a]*cP[e]
        pp = sum(cP[e]*GRAM[e][f]*cP[f] for e in range(3) for f in range(3))
        C += y*(pp - (m2 if mm else 0))
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
    P = sp.Poly(G, *ys, domain='QQ')
    Gd = {tuple(e): Fraction(int(c.p), int(c.q)) for e,c in P.terms()}
    return Gd, False

# ---------------- dict-poly ops ---------------------------------------------
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
    """set x_i=0 and drop coordinate i"""
    r = {}
    for e,c in p.items():
        if e[i]==0:
            r[tuple(v for k,v in enumerate(e) if k!=i)] = c
    return r
def dmaxdeg(p): return max((sum(e) for e in p.items() and p), default=0) if p else 0

# ---------------- Singular syz ----------------------------------------------
def poly_str(p, n):
    if not p: return '0'
    ts = []
    for e,c in sorted(p.items()):
        mono = '*'.join(f'y{i+1}^{ei}' for i,ei in enumerate(e) if ei)
        cs = f'({c.numerator}/{c.denominator})' if c.denominator!=1 else f'({c.numerator})'
        ts.append(cs + ('*'+mono if mono else ''))
    return '+'.join(ts)

def singular_syz(Gd, n, timeout=300):
    gens = [poly_str(dderiv(Gd,i), n) for i in range(n)]
    gens.append(poly_str({e: -c for e,c in Gd.items()}, n))
    vec = ',\n  '.join(f'[{g}]' for g in gens)
    vlist = ','.join(f'y{i+1}' for i in range(n))
    script = f"""ring R = 0,({vlist}),dp;
module m1 = {vec};
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
        syz.append((gen[:n], gen[n]))   # (a_list, b)
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

# ---------------- tower ------------------------------------------------------
def monos(n, dmax, _c={}):
    key = (n, dmax)
    if key not in _c:
        _c[key] = [e for e in itertools.product(range(dmax+1), repeat=n) if sum(e) <= dmax]
    return _c[key]

class Tower:
    def __init__(self, top, eta0, d0, D0, forced_masters, verbose=True, Dtop=None, Dlow=None):
        """top: tuple of prop ids. forced_masters: list of 15-index nu vectors."""
        self.top = tuple(sorted(top)); self.eta0 = eta0; self.d0 = d0
        self.D0 = D0; self.Dtop = Dtop if Dtop is not None else D0
        self.Dlow = Dlow if Dlow is not None else D0
        self.verbose = verbose
        self.stats = {}
        self._build_nodes()
        self._iso_classes()
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
                Gd, isz = build_G(S, self.eta0)
                if isz: self.zero.add(S)
                else: self.G[S] = Gd
        self.stats['nodes_nonzero'] = len(self.G); self.stats['nodes_zero'] = len(self.zero)
        self.stats['t_symanzik'] = round(time.time()-t,2)
        if self.verbose:
            print(f"[tower] nodes: {len(self.G)} nonzero, {len(self.zero)} zero ({self.stats['t_symanzik']}s)", flush=True)

    @staticmethod
    def _iso_find(G1, G2, n):
        """perm sigma with G1 relabeled by sigma == G2; sigma[i]=j means var i of node1 -> var j of node2"""
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
        self.rep = {}; self.perm = {}   # node -> repnode, node->perm (node coords -> rep coords)
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
        # automorphism groups of the representatives (nontrivial perms fixing G)
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
        self.stats['aut_sizes'] = {str(R): len(self.auts[R])+1 for R in reps if self.auts[R]}
        self.stats['t_iso'] = round(time.time()-t,2)
        if self.verbose:
            print(f"[tower] {len(reps)} iso-classes of {len(self.G)} nodes; "
                  f"nontrivial Aut on {sum(1 for R in reps if self.auts[R])} reps ({self.stats['t_iso']}s)", flush=True)

    def _aut_min(self, R, m):
        best = tuple(m)
        for sig in self.auts[R]:
            m2 = [0]*len(m)
            for i,v in enumerate(m): m2[sig[i]] = v
            t2 = tuple(m2)
            if t2 < best: best = t2
        return best

    def canon(self, S, m):
        """(node, exps in node order) -> Aut-canonical (repnode, exps in rep order)"""
        R = self.rep[S]; sig = self.perm[S]
        m2 = [0]*len(S)
        for i,v in enumerate(m): m2[sig[i]] = v
        return R, self._aut_min(R, m2)

    def _map_masters(self, masters_nu):
        """kira master nu-vectors (len 15) -> canonical (rep, m) columns; the master's
        prop-subset may lie outside the tower — match by G-iso against reps."""
        self.master_cols = []; self.master_info = []
        for nu in masters_nu:
            K = tuple(i+1 for i,v in enumerate(nu) if v > 0)
            mK = tuple(nu[i-1]-1 for i in K)
            GK, isz = build_G(K, self.eta0)
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
        if self.verbose:
            print(f"[tower] mapped {len(self.master_cols)} forced masters", flush=True)

    def _columns(self):
        """global column order: per rep (topo: |S| desc, lex), nonmaster monos
        (|m| desc then lex), then ALL master cols at the end."""
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

    # ---------- relations ----------
    def build_relations(self):
        t = time.time(); self.rel_by_node = {R: [] for R in self.node_order}
        self.syz = {}; tsing = 0.0; nrel = 0; ndrop = 0
        for R in self.node_order:
            n = len(R)
            syz, dt = singular_syz(self.G[R], n); tsing += dt
            self.syz[R] = syz
            children = {}
            for i in range(n):
                Sc = tuple(v for k,v in enumerate(R) if k != i)
                if Sc in self.zero: children[i] = None          # zero sector
                elif Sc in self.rep: children[i] = Sc
                else: children[i] = Sc   # subsets of a rep are nodes of the full tower iff rep sub top...
            for (alist, b) in syz:
                da = max((max((sum(e) for e in a), default=0) for a in alist), default=0)
                db = max((sum(e) for e in b), default=0)
                drel = max(da-1, db, 0)
                if da == 0 and db == 0: continue
                DR = self.node_D(R)
                for m in monos(n, max(DR - drel, 0)):
                    row = {}; ok = True
                    # S-part
                    Rm = {}
                    for i in range(n):
                        if alist[i]: Rm = dadd(Rm, dderiv(dmul_mono(alist[i], m), i))
                    if b: Rm = dadd(Rm, dmul_mono(b, m), cq=-self.d0/2)
                    for e,c in Rm.items():
                        key = (R, self._aut_min(R, e))
                        if key not in self.colid: ok = False; break
                        row[self.colid[key]] = row.get(self.colid[key], Fraction(0)) + c
                    if not ok: ndrop += 1; continue
                    # children (boundary)
                    for i in range(n):
                        if m[i] != 0 or not alist[i]: continue
                        bp = drestrict(dmul_mono(alist[i], m), i)
                        if not bp: continue
                        Sc = tuple(v for k,v in enumerate(R) if k != i)
                        if Sc in self.zero: continue
                        for e,c in bp.items():
                            # J_S[R] = -sum_child J_child[bp]  (int_0^inf d_i F = -F|_0)
                            key = self.canon(Sc, e)
                            if key not in self.colid: ok = False; break
                            row[self.colid[key]] = row.get(self.colid[key], Fraction(0)) + c
                        if not ok: break
                    if not ok: ndrop += 1; continue
                    row = {k:v for k,v in row.items() if v != 0}
                    if row:
                        self.rel_by_node[R].append(row); nrel += 1
        self.stats.update(t_singular=round(tsing,2), n_relations=nrel, n_dropped=ndrop,
                          t_build=round(time.time()-t,2),
                          syz_sizes={str(R): len(self.syz[R]) for R in self.node_order})
        if self.verbose:
            print(f"[rel] {nrel} relations ({ndrop} dropped out-of-pool) "
                  f"singular {tsing:.1f}s build {self.stats['t_build']}s", flush=True)

    # ---------- per-node elimination ----------
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
            rr = 0
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

    def expand(self, colid, _memo=None):
        if _memo is None: _memo = self._memo = getattr(self, '_memo', {})
        if colid in _memo: return _memo[colid]
        if colid not in self.rules:
            _memo[colid] = {colid: Fraction(1)}; return _memo[colid]
        out = {}
        for c, v in self.rules[colid]:
            for cc, vv in self.expand(c, _memo).items():
                w = out.get(cc, Fraction(0)) + v*vv
                if w: out[cc] = w
                elif cc in out: del out[cc]
        _memo[colid] = out; return out

    # ---------- I <-> J conversion ----------
    def conv(self, col, sigma):
        R, m = self.cols[col] if isinstance(col, int) else col
        w = sum(m) + len(R)
        c = Fraction(1)
        for mi in m:
            for k in range(2, mi+1): c *= k
        for j in range(1, w+1): c /= (2*self.d0 - j)
        if sigma == -1 and w % 2: c = -c
        return c

    def reduce_target(self, nu, sigma):
        """nu: 15-vector. Returns dict {(rep,m): Fraction} of I-space coefficients."""
        K = tuple(i+1 for i,v in enumerate(nu) if v > 0)
        mK = tuple(nu[i-1]-1 for i in K)
        assert K in self.rep, f"target sector {K} not in tower"
        tc = self.canon(K, mK)
        cid = self.colid[tc]
        sys.setrecursionlimit(100000)
        Jred = self.expand(cid)
        ct = self.conv(tc, sigma)
        out = {}
        for c, v in Jred.items():
            out[self.cols[c]] = v * self.conv(c, sigma) / ct
        return out

# ---------------- kira table parse -------------------------------------------
def parse_kira_entry(path, nu):
    txt = open(path).read()
    pat = r'lbl3mopp_k2sh\[' + r',\s*'.join(str(v) for v in nu) + r'\] -> \n(.*?)(?:\n,\n|\n\})'
    m = re.search(pat, txt, re.S)
    if not m: return None
    body = m.group(1)
    terms = re.findall(r'lbl3mopp_k2sh\[([0-9,\s\-]+)\]\*(\(.*?\))\s*(?=$|\n)', body)
    d, eta = sp.symbols('d eta')
    out = {}
    for idx, coeff in terms:
        nu2 = tuple(int(x) for x in idx.split(','))
        ex = sp.sympify(coeff.replace('^','**'))
        out[nu2] = ex
    return out

def eval_coeff(ex, d0, eta0):
    d, eta = sp.symbols('d eta')
    v = ex.subs({d: sp.Rational(d0.numerator, d0.denominator),
                 eta: sp.Rational(eta0.numerator if hasattr(eta0,'numerator') else eta0, getattr(eta0,'denominator',1))})
    v = sp.nsimplify(v, rational=True)
    v = sp.Rational(v)
    return Fraction(int(v.p), int(v.q))

# ---------------- main -------------------------------------------------------
def rss_gb():
    kib = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform == 'darwin': kib /= 1024  # ru_maxrss: bytes on macOS, KiB on Linux
    return round(kib/1048576, 2)

if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--stage', required=True, choices=['control','sec63','syz255'])
    ap.add_argument('--D0', type=int, default=6)
    ap.add_argument('--Dlow', type=int, default=None)
    ap.add_argument('--Dtop', type=int, default=None)
    ap.add_argument('--d0', default='97/23')
    ap.add_argument('--eta0', default='5/7')
    ap.add_argument('--out', default=None)
    args = ap.parse_args()
    d0 = Fraction(args.d0); eta0f = Fraction(args.eta0)
    eta0 = sp.Rational(eta0f.numerator, eta0f.denominator)
    # KTAB is env-only (TRUST_LPSYZ_KTAB, REQUIRED) — no machine-local
    # default; refuse loudly when unset.
    KTAB = os.environ.get('TRUST_LPSYZ_KTAB')
    if KTAB is None:
        raise SystemExit('TRUST_LPSYZ_KTAB unset — path to the kira_target.m '
                         'truth table is required (no default)')
    res = {'stage': args.stage, 'd0': str(d0), 'eta0': str(eta0f), 'D0': args.D0}

    if args.stage == 'syz255':
        t = time.time()
        Gd, isz = build_G(tuple(range(1,9)), eta0)
        res['t_symanzik_255'] = round(time.time()-t,2)
        res['G_terms'] = len(Gd)
        syz, dt = singular_syz(Gd, 8, timeout=3000)
        res['t_singular_255'] = round(dt,2)
        res['n_syz_255'] = len(syz)
        degs = [max(max((sum(e) for a in al for e in a), default=0), max((sum(e) for e in b), default=0)) for al,b in syz]
        res['syz_degs'] = sorted(set(degs))
        res['rss_gb'] = rss_gb()
        print(json.dumps(res, indent=1))
        if args.out: json.dump(res, open(args.out,'w'), indent=1)
        sys.exit(0)

    if args.stage == 'control':
        top = (3,4,8)
        targets = [(0,0,1,1,0,0,0,2,0,0,0,0,0,0,0),
                   (0,0,1,2,0,0,0,1,0,0,0,0,0,0,0),
                   (0,0,2,1,0,0,0,1,0,0,0,0,0,0,0)]
        masters = [(0,0,1,1,0,0,0,1,0,0,0,0,0,0,0)]
    else:
        top = (1,2,3,4,5,6)
        targets = [(2,1,1,1,1,1,0,0,0,0,0,0,0,0,0),
                   (1,1,2,1,1,1,0,0,0,0,0,0,0,0,0),
                   (1,1,1,1,1,2,0,0,0,0,0,0,0,0,0),
                   (1,2,1,2,1,1,0,0,0,0,0,0,0,0,0)]
        txt = open(KTAB).read()
        masters = set()
        for nu in targets:
            ent = parse_kira_entry(KTAB, nu)
            assert ent, f'no kira entry for {nu}'
            for k in ent: masters.add(k)
        masters = sorted(masters)
        print(f"[main] {len(masters)} kira masters over the 3 targets", flush=True)

    tw = Tower(top, eta0, d0, args.D0, masters, Dtop=args.Dtop, Dlow=args.Dlow)
    tw.build_relations()
    tw.solve()
    res.update(tw.stats)
    res['free_nonmaster_sample'] = [f"{c}" for c in tw.free_nonmaster[:20]]
    res['viol_rows_named'] = [{f"{tw.cols[c]}": str(v) for c,v in r.items()} for r in tw.viol_rows]

    ok_all = {}
    for sigma in (1,-1):
        allok = True; details = []
        for nu in targets:
            ent = parse_kira_entry(KTAB, nu)
            mine = tw.reduce_target(nu, sigma)
            # map kira masters to canonical cols
            kmap = {}
            okt = True; diffs = []
            for knu, kex in ent.items():
                kc = None
                K = tuple(i+1 for i,v in enumerate(knu) if v>0)
                mK = tuple(knu[i-1]-1 for i in K)
                GK,_ = build_G(K, eta0)
                if K in tw.rep: kc = tw.canon(K, mK)
                else:
                    for R in tw.reps:
                        if len(R)!=len(K): continue
                        sig = Tower._iso_find(GK, tw.G[R], len(K))
                        if sig is not None:
                            m2=[0]*len(K)
                            for i,v in enumerate(mK): m2[sig[i]]=v
                            kc=(R,tw._aut_min(R,tuple(m2))); break
                assert kc, f'kira master {knu} unmatched'
                kmap[kc] = kmap.get(kc, Fraction(0)) + eval_coeff(kex, d0, eta0f)
            mine2 = {k:v for k,v in mine.items() if v != 0}
            forced = set(tw.master_cols)
            leak = [k for k in mine2 if k not in forced]
            keys = set(mine2) | set(kmap)
            for k in keys:
                a = mine2.get(k, Fraction(0)); b = kmap.get(k, Fraction(0))
                if a != b:
                    okt = False; diffs.append((str(k), str(a), str(b)))
            # if plain mismatch: test Delta in rowspan(violation rows) over masters
            mod_ok = okt
            if not okt and tw.viol_rows and not leak:
                mcols = sorted({c for r in tw.viol_rows for c in r} |
                               {tw.colid[k] for k in keys if k in tw.colid})
                mi = {c:i for i,c in enumerate(mcols)}
                V = flint.fmpq_mat(len(tw.viol_rows)+1, len(mcols))
                for ri,r in enumerate(tw.viol_rows):
                    for c,v in r.items(): V[ri, mi[c]] = flint.fmpq(v.numerator, v.denominator)
                r1 = flint.fmpq_mat([[V[i,j] for j in range(len(mcols))] for i in range(len(tw.viol_rows))]).rank()
                for k in keys:
                    dv = mine2.get(k, Fraction(0)) - kmap.get(k, Fraction(0))
                    if dv: V[len(tw.viol_rows), mi[tw.colid[k]]] = flint.fmpq(dv.numerator, dv.denominator)
                mod_ok = (V.rank() == r1)
            details.append({'target': list(nu), 'match': okt, 'match_mod_degeneracy': mod_ok,
                            'n_masters_kira': len(kmap), 'n_masters_mine': len(mine2),
                            'n_leak': len(leak), 'leak_sample': [str(k) for k in leak[:5]],
                            'diffs': diffs[:6]})
            if not mod_ok: allok = False
        ok_all[sigma] = (allok, details)
        if allok:
            res['sigma'] = sigma; break
    for sigma,(allok,details) in ok_all.items():
        res[f'match_sigma_{sigma}'] = allok
        res[f'details_sigma_{sigma}'] = details
    res['rss_gb'] = rss_gb(); res['wall'] = round(time.time()-T0,1)
    print(json.dumps({k:v for k,v in res.items() if not k.startswith('details')}, indent=1, default=str))
    if args.out: json.dump(res, open(args.out,'w'), indent=1, default=str)
