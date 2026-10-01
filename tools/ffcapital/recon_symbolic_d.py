#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Part of BootLoops 1.0 (repository-root
# LICENSE and NOTICE). Reads the plain-text state files FireFly writes; contains
# no FireFly code (THIRD_PARTY.md B2).
# ffcapital member — cross-slice symbolic-d salvage.
# Siblings: harvest_ffsave_eta.py (1-var), harvest_ffsave_2var.py (2-var).
# All roots via CLI/--config; fail-closed on unset CONFIG.
r"""recon_symbolic_d.py — symbolic-(d,eta) reconstruction from a farm's
numeric-d eta-slices, with the full gate battery.

INPUT (all exact rationals, no floats anywhere on the value path):
  - per-slice harvest jsonl (harvest_ffsave_eta.py output): per fn,
    N/D = {eta_exp: "p/q"} — the 1-var eta rational function at that slice's
    numeric d, defined up to an overall per-slice scalar (FireFly's per-slice
    normalization sets SOME coefficient to 1/1; which one varies by fn/slice).
  - wave1_perfn_coverage.jsonl.gz: per fn: tag, nd, dd (deg_d of num/den from
    the upstream degree census), need_fit=nd+dd+5,
    need_total=nd+dd+9, cov, slice list.
  - ffsave_9prime_bank/ff_save: a (2-var, symbolic-d) FireFly state of the
    same reduction at 9 primes = the NEVER-IN-FIT oracle.

NORMALIZATION CONVENTION (the cross-slice pin, established before any fitting;
verified by `probe3` and the synthetic control):
  Per-slice scalars are eliminated by fitting ONLY scalar-free ratios.
  Per fn choose anchor positions aN (num) and aD (den) = the lowest eta-exponent
  nonzero on ALL covered slices (fallback: max-count position; slices where an
  anchor vanishes are dropped from that fn's usable set).
    num position k: G_k(d) = A_k(d)/B_aD(d)   — rational of type (nd, dd)
    den position j: H_j(d) = B_j(d)/A_aN(d)   — rational of type (dd, nd)
  Every fitted object has num+den d-degree <= nd+dd, so the coverage law
  need_fit = nd+dd+5 (nd+dd+1 interpolation + 4 guard) applies uniformly.
  Internal identity G_aN * H_aD == 1 is checked exactly per fn.
  Reassembly: F = sum_k G_k eta^k / sum_j (H_j/H_aD) eta^j, cleared to
  P(d,eta), Q(d,eta) in Z[d,eta], jointly content-primitive, sign fixed by
  leading d-coefficient of Q's aD eta-coefficient > 0.

FIT / HELD-OUT SPLIT: usable slices sorted by d ascending; fit = need_fit
slices spread across the d-range (deterministic even-stride pick); held-out =
ALL remaining (>= 4 by the coverage bar). Held-out slices never enter any fit.

CANCELLATION SLICES: at special d a
fn's P(d0,eta), Q(d0,eta) can acquire a common eta-factor; FireFly's slice then
carries the REDUCED form — the same rational function, but its coefficient
vector is NOT on the generic coefficient curves (observed in production: two
fns at d=118/29, common factor eta+22/5). Detection: slice (degN,degD) BOTH below the
fn's cross-slice generic degrees (min joint drop >= 1). Such slices are
auto-EXCLUDED from all coefficient fits and counted toward validation only,
via the cancellation-invariant cross-multiplication check below.

GATE BATTERY (ALL must pass to bank; any failure voids the fn -> RETRY list):
  (fit)      overdetermined exact nullspace fit (the 4 guard points make the
             (nd+dd+5)x(nd+dd+2) system overdetermined by >=3; a nontrivial
             nullspace exists only if the data is truly rational of the type)
  (internal) G_aN * H_aD == 1 exactly; assembled deg_d(P)<=nd, deg_d(Q)<=dd
  (a) heldout: for EVERY validation slice (held-out + anchor-dropped +
      cancellation slices), rational-function equality by cross-multiplication:
      P(d_s,eta)*D_slice(eta) == lambda * Q(d_s,eta)*N_slice(eta) exactly,
      single lambda (cancellation-invariant; a plain coefficient-vector
      comparison fails on reduced slices)
  (b) loo: for every fit slice, every position refit without it must return
      the IDENTICAL (N,D)
  (c) oracle: vs the 9-prime bank, never in any fit. Bank-done fn: exact
      2-var cross-multiplication identity P*bankD == Q*bankN. Bank-active fn:
      per prime p_i of the bank's own combined_prime chain (factored against
      the FireFly 63-bit table), the bank's committed residues/exact entries
      must equal lambda_i * (our integer coefficients) mod p_i for a single
      lambda_i per prime across num AND den.
  (d) synthetic-truth control (subcommand `control`): a KNOWN (d,eta) rational
      function of representative degrees, sliced at the run's OWN d-grid (the
      coverage file's slice names) with FireFly-style per-slice normalization,
      pushed through THIS code path, must be recovered exactly and pass every
      gate — plus negative controls (corrupted guard point, corrupted held-out
      point, corrupted oracle residue must each FAIL loudly). The planted
      vanishing/cancellation d-points are picked FROM the grid
      (control_plants); grids too small to satisfy the coverage law for the
      synthetic degrees are refused loudly, never fitted thin.

Output serialization (one json line per banked fn):
  {"fn", "tag", "nd", "dd", "N": {eta_exp: [c0,...]}, "D": {eta_exp: [c0,...]},
   "fit_slices", "heldout_slices", "anchors": [aN, aD], "dropped_slices",
   "gates": {...}}
  N/D: integer d-coefficient lists (ascending powers, decimal strings),
  jointly content-primitive, sign convention as above. F = N/D at (d, eta).
"""
import argparse
import glob
import gzip
import json
import math
import os
import re
import sys
import time
from fractions import Fraction

from flint import fmpq, fmpq_poly, fmpz, fmpz_mat, fmpz_poly

# --- wiring: the 6 module constants are CLI/config-driven (flag > --config
# JSON key > derived-from-BASE); nothing is silently defaulted here — an
# unset/missing path fails LOUDLY at first use (fail-closed).
BASE = None          # run root                     (--base / config "base")
WAVE1 = None         # per-slice harvest dir        (--wave1; derives <BASE>/wave1_harvest)
P2_FILE = None       # extra single-slice harvest   (--p2-file; derives <BASE>/p2_harvest_d86_23.jsonl)
BANK_STATES = None   # never-in-fit oracle bank     (--bank-states; derives <BASE>/ffsave_9prime_bank/ff_save/states)
FF_HELPER = None     # FireFly ReconstHelper.cpp    (--ff-helper; the prime-table source)
COVERAGE = None      # per-fn coverage jsonl.gz     (--coverage; derives <WAVE1>/wave1_perfn_coverage.jsonl.gz)


def _require(val, name, flag):
    if val is None:
        raise SystemExit(f'CONFIG: {name} unset — pass {flag} (or --config JSON'
                         f' / --base derivation); no silent defaults')
    if not os.path.exists(val):
        raise SystemExit(f'CONFIG: {name} path does not exist: {val}')
    return val


def resolve_config(args):
    """Flag > --config JSON key > derived-from-BASE. Sets the module globals
    (mp workers inherit them via fork; Pool is created after this runs)."""
    global BASE, WAVE1, P2_FILE, BANK_STATES, FF_HELPER, COVERAGE
    cfg = {}
    if args.config:
        with open(args.config) as f:
            cfg = json.load(f)

    def pick(cli, key):
        return cli if cli is not None else cfg.get(key)

    BASE = pick(args.base, 'base')
    WAVE1 = pick(args.wave1, 'wave1') or (
        BASE and os.path.join(BASE, 'wave1_harvest'))
    P2_FILE = pick(args.p2_file, 'p2_file') or (
        BASE and os.path.join(BASE, 'p2_harvest_d86_23.jsonl'))
    BANK_STATES = pick(args.bank_states, 'bank_states') or (
        BASE and os.path.join(BASE, 'ffsave_9prime_bank', 'ff_save', 'states'))
    FF_HELPER = pick(args.ff_helper, 'ff_helper')
    COVERAGE = pick(args.coverage, 'coverage') or (
        WAVE1 and os.path.join(WAVE1, 'wave1_perfn_coverage.jsonl.gz'))

t0 = time.time()


def log(msg):
    print(f'[{time.time()-t0:7.1f}s] {msg}', flush=True)


# ---------------------------------------------------------------- primes
def firefly_primes(n=40):
    txt = open(_require(FF_HELPER, 'FF_HELPER', '--ff-helper')).read()
    m = re.search(r'the 300 largest 63-bit primes.*?\{(.*?)\}', txt, re.S)
    ps = [int(x) for x in re.findall(r'(\d+)uLL', m.group(1))]
    assert ps[0] == 9223372036854775783
    return ps[:n]


def factor_combined_prime(M, ptab):
    """combined_prime must be a product of a prefix of the FireFly table."""
    primes = []
    for p in ptab:
        if M == 1:
            break
        if M % p:
            raise ValueError('combined_prime not a prefix product')
        M //= p
        primes.append(p)
    if M != 1:
        raise ValueError('combined_prime not exhausted by table prefix')
    return primes


# ---------------------------------------------------------------- slices
def slice_d(name):
    if name == 'p2_d86_23':
        return Fraction(86, 23)
    m = re.match(r'node_eps_(-?\d+)_(\d+)(_HELDOUT)?$', name)
    if not m:
        raise ValueError(f'unknown slice name {name}')
    eps = Fraction(int(m.group(1)), int(m.group(2)))
    return 4 - 2 * eps


# NEVER-FIT slices (validation-only):
#  - *_HELDOUT nodes: NEVER enter fits, ever.
#  - Slices once quarantined and later exonerated (e.g. a common-factor
#    cancellation point misread as corruption) can be listed in
#    QUARANTINE_VALIDATION_ONLY: they then validate but never enter fits.
#    Cancellation slices are in any case auto-excluded from coefficient fits
#    per-fn by the shape detector, which is the correct guard.
QUARANTINE_VALIDATION_ONLY = set()


def never_fit(name):
    return name.endswith('_HELDOUT') or name in QUARANTINE_VALIDATION_ONLY


def slice_path(name):
    if name == 'p2_d86_23':
        return _require(P2_FILE, 'P2_FILE', '--p2-file')
    _require(WAVE1, 'WAVE1', '--wave1')
    return os.path.join(WAVE1, f'{name}_harvest.jsonl')


def build_offset_index(path):
    """fn -> byte offset of its jsonl line."""
    idx = {}
    off = 0
    with open(path, 'rb') as f:
        for line in f:
            # cheap fn extraction: line starts {"fn": N,
            b = line.find(b':')
            e = line.find(b',', b)
            fn = int(line[b + 1:e])
            idx[fn] = off
            off += len(line)
    return idx


class SliceStore:
    """Lazy random access to per-slice harvest rows by fn index."""

    def __init__(self, names):
        self.names = names
        self.d = {n: slice_d(n) for n in names}
        self.paths = {n: slice_path(n) for n in names}
        self.idx = {}
        self.fh = {}

    def build_indices(self):
        for n in self.names:
            self.idx[n] = build_offset_index(self.paths[n])
            log(f'indexed {n}: {len(self.idx[n])} fns')

    def row(self, name, fn):
        off = self.idx[name].get(fn)
        if off is None:
            return None
        fh = self.fh.get(name)
        if fh is None:
            fh = self.fh[name] = open(self.paths[name], 'rb')
        fh.seek(off)
        r = json.loads(fh.readline())
        assert r['fn'] == fn
        return r


def parse_frac(s):
    """'p/q' -> (int p, int q>0)."""
    if '/' in s:
        a, b = s.split('/')
        a, b = int(a), int(b)
    else:
        a, b = int(s), 1
    if b < 0:
        a, b = -a, -b
    return a, b


# ---------------------------------------------------------------- exact fit
class FitFail(Exception):
    pass


def ratfit_pq(xs_uv, ys_nd, p, q):
    """Exact rational fit y(x)=N(x)/D(x), deg N<=p, deg D<=q, from ALL points
    (overdetermined when len>p+q+1). Returns (N, D) ascending int-coeff lists,
    primitive, gcd-reduced, D-lead>0. Raises FitFail if no consistent rational
    function of this type exists (nullity 0) or verification fails."""
    K = len(xs_uv)
    P = max(p, q)
    Mat = fmpz_mat(K, p + q + 2)
    for k in range(K):
        A, B = ys_nd[k]
        u, v = xs_uv[k]
        vp = [1] * (P + 1)
        for b in range(1, P + 1):
            vp[b] = vp[b - 1] * v
        pu = 1
        pw = [0] * (P + 1)
        for b in range(P + 1):
            pw[b] = pu * vp[P - b]
            pu *= u
        for b in range(p + 1):
            Mat[k, b] = B * pw[b]
        for b in range(q + 1):
            Mat[k, p + 1 + b] = -A * pw[b]
    NS, nullity = Mat.nullspace()
    if nullity == 0:
        raise FitFail(f'nullity 0 (K={K}, p={p}, q={q})')
    col = None
    for c in range(nullity):
        if any(NS[p + 1 + b, c] != 0 for b in range(q + 1)):
            col = c
            break
    if col is None:
        raise FitFail('nullspace has no nonzero denominator column')
    Np = fmpz_poly([int(NS[b, col]) for b in range(p + 1)])
    Dp = fmpz_poly([int(NS[p + 1 + b, col]) for b in range(q + 1)])
    g = Np.gcd(Dp)
    if g.degree() > 0:
        Np //= g
        Dp //= g
    N = [int(Np[i]) for i in range(Np.degree() + 1)] if Np != 0 else [0]
    D = [int(Dp[i]) for i in range(Dp.degree() + 1)]
    g = 0
    for c in N + D:
        g = math.gcd(g, c)
    if g > 1:
        N = [c // g for c in N]
        D = [c // g for c in D]
    if D[-1] < 0:
        N = [-c for c in N]
        D = [-c for c in D]
    # verify against every input point (protects nullity>1 column choice)
    if not verify_points(N, D, xs_uv, ys_nd):
        raise FitFail('fit does not reproduce its own input points')
    return N, D


def poly_eval_uv(P, u, v):
    """ascending int coeffs at x=u/v -> (num, den) unreduced ints."""
    num, den = 0, 1
    for i, c in enumerate(P):
        if c:
            tn, td = c * (u ** i), v ** i
            num = num * td + tn * den
            den = den * td
    return num, den


def verify_points(N, D, xs_uv, ys_nd):
    for (u, v), (A, B) in zip(xs_uv, ys_nd):
        nn, nd_ = poly_eval_uv(N, u, v)
        dn, dd_ = poly_eval_uv(D, u, v)
        if nn * dd_ * B != dn * nd_ * A:
            return False
    return True


def poly_lcm(a, b):
    if a == 0:
        return b
    if b == 0:
        return a
    g = a.gcd(b)
    return a * (b // g)


# ---------------------------------------------------------------- per-fn core
def choose_anchor(rows, side):
    """positions nonzero on the most slices; prefer lowest eta exponent."""
    count = {}
    for r in rows.values():
        for k, v in r[side].items():
            if parse_frac(v)[0] != 0:
                count[int(k)] = count.get(int(k), 0) + 1
    if not count:
        return None, 0
    best = max(count.values())
    anchor = min(k for k, c in count.items() if c == best)
    return anchor, best


def spread_pick(n_avail, n_fit):
    """deterministic even-spread index pick, then fill low-to-high."""
    if n_fit >= n_avail:
        return list(range(n_avail))
    sel = sorted({(i * (n_avail - 1)) // (n_fit - 1) for i in range(n_fit)}
                 if n_fit > 1 else {0})
    i = 0
    while len(sel) < n_fit:
        if i not in sel:
            sel.append(i)
        i += 1
    return sorted(sel)


def get_val(row, side, k):
    v = row[side].get(str(k))
    return parse_frac(v) if v is not None else (0, 1)


def reconstruct_fn(fn, cov_row, store, oracle_reader, ptab, loo=True):
    """Full pipeline for one fn. Returns dict with status + gates + result."""
    nd, dd = cov_row['nd'], cov_row['dd']
    need_fit = nd + dd + 5
    res = {'fn': fn, 'tag': cov_row['tag'], 'nd': nd, 'dd': dd,
           'status': None, 'gates': {}}
    rows = {}
    for s in cov_row['slices']:
        r = store.row(s, fn)
        if r is None or not r.get('done') or r.get('validated') is False:
            continue
        if r['tag'] != cov_row['tag']:
            res['status'] = 'TAG_MISMATCH'
            return res
        rows[s] = r
    # cancellation-slice detection: a common eta-factor of
    # P,Q at that slice's d drops BOTH eta-degrees below the generic shape;
    # such slices carry the reduced form and must not enter coefficient fits
    def degpair(r):
        dN = max((int(k) for k, v in r['N'].items()
                  if parse_frac(v)[0] != 0), default=-1)
        dD = max((int(k) for k, v in r['D'].items()
                  if parse_frac(v)[0] != 0), default=-1)
        return dN, dD
    gen_N = max((degpair(r)[0] for r in rows.values()), default=-1)
    gen_D = max((degpair(r)[1] for r in rows.values()), default=-1)
    cancel = [s for s, r in rows.items()
              if min(gen_N - degpair(r)[0], gen_D - degpair(r)[1]) >= 1]
    rows_clean = {s: r for s, r in rows.items() if s not in cancel}

    aN, cN = choose_anchor(rows_clean, 'N')
    aD, cD = choose_anchor(rows_clean, 'D')
    if aN is None or aD is None:
        res['status'] = 'NO_ANCHOR'
        return res
    usable, dropped, neverfit = [], [], []
    for s, r in rows_clean.items():
        if get_val(r, 'N', aN)[0] != 0 and get_val(r, 'D', aD)[0] != 0:
            (neverfit if never_fit(s) else usable).append(s)
        else:
            dropped.append(s)
    if len(usable) < need_fit or \
            (len(usable) - need_fit + len(dropped) + len(cancel)
             + len(neverfit)) < 4:
        res['status'] = 'INSUFFICIENT_AFTER_ANCHOR'
        res['usable'] = len(usable)
        return res
    usable.sort(key=lambda s: store.d[s])
    sel = spread_pick(len(usable), need_fit)
    fit_slices = [usable[i] for i in sel]
    held = [usable[i] for i in range(len(usable)) if i not in set(sel)]
    validation = held + dropped + cancel + neverfit
    res['anchors'] = [aN, aD]
    res['fit_slices'] = fit_slices
    res['heldout_slices'] = validation
    res['dropped_slices'] = dropped
    res['cancel_slices'] = cancel
    res['neverfit_slices'] = neverfit

    # support union over usable slices
    SN = sorted({int(k) for s in usable for k in rows[s]['N']})
    SD = sorted({int(k) for s in usable for k in rows[s]['D']})

    def ratio(r, side, k, aside, a):
        A1, B1 = get_val(r, side, k)
        A2, B2 = get_val(r, aside, a)
        # (A1/B1)/(A2/B2) = A1*B2 / (B1*A2)
        return (A1 * B2, B1 * A2)

    xs_fit = [(store.d[s].numerator, store.d[s].denominator)
              for s in fit_slices]

    def fit_all(slc_names, xs):
        G, H = {}, {}
        for k in SN:
            ys = [ratio(rows[s], 'N', k, 'D', aD) for s in slc_names]
            G[k] = ratfit_pq(xs, ys, nd, dd)
        for j in SD:
            ys = [ratio(rows[s], 'D', j, 'N', aN) for s in slc_names]
            H[j] = ratfit_pq(xs, ys, dd, nd)
        return G, H

    try:
        G, H = fit_all(fit_slices, xs_fit)
    except FitFail as e:
        res['status'] = 'FIT_FAIL'
        res['gates']['fit'] = f'FAIL: {e}'
        return res
    res['gates']['fit'] = 'pass'

    # internal identity G_aN * H_aD == 1
    NG, DG = (fmpz_poly(G[aN][0]), fmpz_poly(G[aN][1]))
    NH, DH = (fmpz_poly(H[aD][0]), fmpz_poly(H[aD][1]))
    if NG * NH != DG * DH:
        res['status'] = 'INTERNAL_FAIL'
        res['gates']['internal'] = 'FAIL: G_aN*H_aD != 1'
        return res
    res['gates']['internal'] = 'pass'

    # assemble P(d,eta), Q(d,eta): W_j = H_j / H_aD
    W = {}
    for j in SD:
        num = fmpz_poly(H[j][0]) * DH
        den = fmpz_poly(H[j][1]) * NH
        g = num.gcd(den)
        if g.degree() > 0 or g != 1:
            num //= g
            den //= g
        W[j] = (num, den)
    L = fmpz_poly([1])
    for k in SN:
        L = poly_lcm(L, fmpz_poly(G[k][1]))
    for j in SD:
        L = poly_lcm(L, W[j][1])
    Pk, Qj = {}, {}
    for k in SN:
        Nk, Dk = fmpz_poly(G[k][0]), fmpz_poly(G[k][1])
        Pk[k] = Nk * (L // Dk)
    for j in SD:
        Qj[j] = W[j][0] * (L // W[j][1])
    # joint integer content + sign normalization
    g = fmpz(0)
    for p in list(Pk.values()) + list(Qj.values()):
        for i in range(p.degree() + 1):
            g = g.gcd(p[i])
    if g > 1:
        Pk = {k: p // fmpz_poly([g]) for k, p in Pk.items()}
        Qj = {j: p // fmpz_poly([g]) for j, p in Qj.items()}
    qa = Qj[aD]
    if qa[qa.degree()] < 0:
        Pk = {k: -p for k, p in Pk.items()}
        Qj = {j: -p for j, p in Qj.items()}
    degP = max((p.degree() for p in Pk.values() if p != 0), default=0)
    degQ = max((p.degree() for p in Qj.values() if p != 0), default=0)
    if degP > nd or degQ > dd:
        res['status'] = 'DEGCK_FAIL'
        res['gates']['degck'] = f'FAIL: deg_d P={degP}(<= {nd}?) Q={degQ}(<= {dd}?)'
        return res
    res['gates']['degck'] = 'pass'

    # ---- gate (a): every validation slice (held-out + dropped + cancel),
    # rational-function equality by cross-multiplication (single lambda;
    # cancellation-invariant)
    for s in validation:
        if not crossmult_ok(Pk, Qj, rows[s], store.d[s]):
            res['status'] = 'HELDOUT_FAIL'
            res['gates']['heldout'] = f'FAIL@{s}: cross-mult function equality'
            return res
    res['gates']['heldout'] = f'pass({len(validation)})'

    # ---- gate (b): leave-one-out on the fit set
    if loo:
        for i in range(len(fit_slices)):
            sub = [fit_slices[x] for x in range(len(fit_slices)) if x != i]
            xs_sub = [xs_fit[x] for x in range(len(fit_slices)) if x != i]
            try:
                G2, H2 = fit_all(sub, xs_sub)
            except FitFail as e:
                res['status'] = 'LOO_FAIL'
                res['gates']['loo'] = f'FAIL drop {fit_slices[i]}: {e}'
                return res
            if G2 != G or H2 != H:
                res['status'] = 'LOO_FAIL'
                res['gates']['loo'] = f'FAIL drop {fit_slices[i]}: fit moved'
                return res
        res['gates']['loo'] = f'pass({len(fit_slices)})'
    else:
        res['gates']['loo'] = 'skipped'

    # ---- gate (c): 9-prime bank oracle (never in any fit)
    oracle = oracle_reader(fn)
    if oracle is None:
        res['status'] = 'ORACLE_MISSING'
        res['gates']['oracle'] = 'FAIL: no bank state'
        return res
    ok, msg = oracle_check(Pk, Qj, oracle, ptab)
    res['gates']['oracle'] = msg
    if not ok:
        res['status'] = 'ORACLE_FAIL'
        return res

    res['status'] = 'BANKED'
    res['N'] = {str(k): [str(int(p[i])) for i in range(p.degree() + 1)]
                for k, p in Pk.items() if p != 0}
    res['D'] = {str(j): [str(int(p[i])) for i in range(p.degree() + 1)]
                for j, p in Qj.items() if p != 0}
    return res


def conv1(a, b):
    """1-var coefficient-dict product (Fraction values)."""
    out = {}
    for k1, c1 in a.items():
        for k2, c2 in b.items():
            k = k1 + k2
            out[k] = out.get(k, 0) + c1 * c2
    return {k: v for k, v in out.items() if v}


def crossmult_ok(Pk, Qj, row, dfr):
    """Rational-function equality of assembled F(d_s, .) vs the slice's N/D:
    P(d_s,eta)*D_slice(eta) == lambda * Q(d_s,eta)*N_slice(eta), one lambda.
    Invariant under slice-local num/den common-factor cancellation."""
    u, v = dfr.numerator, dfr.denominator

    def at_d(table):
        out = {}
        for k, p in table.items():
            q = Fraction(*poly_eval_uv(
                [int(p[i]) for i in range(p.degree() + 1)], u, v))
            if q:
                out[k] = q
        return out

    A, B = at_d(Pk), at_d(Qj)
    Ns = {int(k): Fraction(*parse_frac(x)) for k, x in row['N'].items()}
    Ds = {int(k): Fraction(*parse_frac(x)) for k, x in row['D'].items()}
    Ns = {k: x for k, x in Ns.items() if x}
    Ds = {k: x for k, x in Ds.items() if x}
    lhs = conv1(A, Ds)
    rhs = conv1(B, Ns)
    lam = None
    for k in sorted(set(lhs) | set(rhs)):
        a, b = lhs.get(k, 0), rhs.get(k, 0)
        if b == 0:
            if a != 0:
                return False
        else:
            l = a / b
            if lam is None:
                lam = l
            elif l != lam:
                return False
    return lam is not None and lam != 0


# ---------------------------------------------------------------- oracle
STATE_SECTIONS = ('combined_prime', 'tag_name', 'is_done', 'max_deg_num',
                  'max_deg_den', 'individual_degrees_num',
                  'individual_degrees_den', 'need_prime_shift',
                  'normalizer_deg', 'normalize_to_den', 'normalizer_den_num',
                  'shifted_max_num_eqn', 'shift', 'sub_num', 'sub_den',
                  'zero_degs_num', 'zero_degs_den', 'g_ni', 'g_di',
                  'combined_ni', 'combined_di', 'combined_primes_ni',
                  'combined_primes_di', 'interpolations')


def read_bank_state(fn, states_dir=None):
    if states_dir is None:
        states_dir = _require(BANK_STATES, 'BANK_STATES', '--bank-states')
    paths = glob.glob(os.path.join(states_dir, f'{fn}_*.gz'))
    if not paths:
        return None
    if len(paths) > 1:  # rename-forward should leave exactly one
        paths.sort(key=lambda p: int(p.rsplit('_', 1)[1][:-3]))
    txt = gzip.open(paths[-1], 'rt').read()
    sec = {}
    cur = None
    for line in txt.splitlines():
        if line in STATE_SECTIONS:
            cur = line
            sec[cur] = []
        elif cur is not None:
            sec[cur].append(line)
    out = {'is_done': int(sec['is_done'][0]),
           'tag_name': sec['tag_name'][0].strip(),
           'combined_prime': int(sec['combined_prime'][0]),
           'shift': [int(x) for x in sec['shift'][0].split()],
           'g_n': {}, 'g_d': {}, 'r_n': {}, 'r_d': {}}
    for key, dst in (('g_ni', 'g_n'), ('g_di', 'g_d')):
        for line in sec.get(key, []):
            p = line.split()
            if len(p) == 4:
                out[dst][(int(p[0]), int(p[1]))] = (int(p[2]), int(p[3]))
    for key, dst in (('combined_ni', 'r_n'), ('combined_di', 'r_d')):
        for line in sec.get(key, []):
            p = line.split()
            if len(p) == 3:
                out[dst][(int(p[0]), int(p[1]))] = int(p[2])
    return out


def _monomials(table):
    """{eta_exp: fmpz_poly in d} -> {(e_d, e_eta): int}."""
    out = {}
    for e_eta, poly in table.items():
        for e_d in range(poly.degree() + 1):
            c = int(poly[e_d])
            if c:
                out[(e_d, e_eta)] = c
    return out


def oracle_check(Pk, Qj, oracle, ptab):
    """gate (c). Pk/Qj: {eta_exp: fmpz_poly}. oracle: read_bank_state dict."""
    if oracle['shift'] and any(oracle['shift']):
        return False, 'FAIL: bank state has nonzero shift (unhandled)'
    ours_n = _monomials(Pk)
    ours_d = _monomials(Qj)
    if oracle['is_done']:
        # exact cross-multiplication: P * bankD == Q * bankN  (2-var, exact)
        bn, bd = oracle['g_n'], oracle['g_d']
        lhs = conv2(ours_n, frac_clear(bd))
        rhs = conv2(ours_d, frac_clear(bn))
        # frac_clear scales bn and bd by DIFFERENT integers; compare as
        # rational multiples: cross-multiply with the two clear factors
        lam1 = frac_clear_den(bd)
        lam2 = frac_clear_den(bn)
        lhs = {k: v * lam2 for k, v in lhs.items()}
        rhs = {k: v * lam1 for k, v in rhs.items()}
        if lhs != rhs:
            return False, 'FAIL: exact 2-var identity vs bank-done state'
        return True, 'pass(exact-vs-bank-done)'
    # active: residue check per prime, single lambda_i across num+den
    M = oracle['combined_prime']
    primes = factor_combined_prime(M, ptab)
    committed = {}
    for (mono, val) in oracle['g_n'].items():
        committed[('N', mono)] = ('g', val)
    for (mono, val) in oracle['g_d'].items():
        committed[('D', mono)] = ('g', val)
    for (mono, r) in oracle['r_n'].items():
        committed[('N', mono)] = ('r', r)
    for (mono, r) in oracle['r_d'].items():
        committed[('D', mono)] = ('r', r)
    ours = {('N', m): c for m, c in ours_n.items()}
    ours.update({('D', m): c for m, c in ours_d.items()})
    # support: every nonzero coefficient of ours must be committed in bank
    missing = [k for k in ours if k not in committed]
    if missing:
        return False, f'FAIL: {len(missing)} fitted monomials absent from bank commit (e.g. {missing[0]})'
    nch = 0
    for p in primes:
        # bank value mod p per committed key
        lam = None
        ref = None
        for key, c in ours.items():
            if c % p:
                ref = key
                break
        if ref is None:
            return False, 'FAIL: all fitted coefficients vanish mod p (degenerate)'
        bt, bv = committed[ref]
        bref = (bv[0] * pow(bv[1], -1, p)) % p if bt == 'g' else bv % p
        lam = (bref * pow(ours[ref] % p, -1, p)) % p
        if lam == 0:
            return False, f'FAIL: lambda=0 mod {p}'
        for key, (bt, bv) in committed.items():
            if bt == 'g':
                if bv[1] % p == 0:
                    continue  # denominator not invertible mod this prime
                bval = (bv[0] * pow(bv[1], -1, p)) % p
            else:
                bval = bv % p
            oval = (lam * (ours.get(key, 0) % p)) % p
            if oval != bval:
                return False, f'FAIL: residue mismatch at {key} mod {p}'
            nch += 1
    return True, f'pass(residues:{len(primes)}p x {len(committed)}c)'


def conv2(a, b):
    """2-var monomial-dict product, integer coefficients."""
    out = {}
    for (d1, e1), c1 in a.items():
        for (d2, e2), c2 in b.items():
            k = (d1 + d2, e1 + e2)
            out[k] = out.get(k, 0) + c1 * c2
    return {k: v for k, v in out.items() if v}


def frac_clear(tab):
    """{mono: (num,den)} -> {mono: int} scaled by lcm of dens."""
    L = 1
    for (_, den) in tab.values():
        L = L * den // math.gcd(L, den)
    return {m: num * (L // den) for m, (num, den) in tab.items()}


def frac_clear_den(tab):
    L = 1
    for (_, den) in tab.values():
        L = L * den // math.gcd(L, den)
    return L


# ---------------------------------------------------------------- coverage
def load_coverage(path=None):
    if path is None:
        path = _require(COVERAGE, 'COVERAGE', '--coverage')
    rows = {}
    for line in gzip.open(path, 'rt'):
        r = json.loads(line)
        rows[r['fn']] = r
    return rows


# ---------------------------------------------------------------- control
def control_plants(ds):
    """Grid-parametrized plant points for the synthetic control.

    The planted cancellation must sit ON a grid slice (the whole point is
    that reconstruct_fn sees a genuinely degree-dropped slice), so both
    plant d-points are picked FROM the run's own grid instead of being
    hard-wired. The reference grid (which contains d=118/29) keeps its
    banked plant values so the byte-identical-control law still holds; any
    other grid gets deterministic picks: the cancellation at an unrepeated
    interior d-value, the vanishing factor at the smallest other distinct d.
    Refuses loudly (SystemExit) when no valid plant exists — never a bare
    crash, never a silently moved plant.
    """
    legacy_cancel, legacy_vanish = Fraction(118, 29), Fraction(55, 13)
    if legacy_cancel in ds:
        return legacy_cancel, legacy_vanish
    uniq = sorted(set(ds))
    once = [d for d in uniq if ds.count(d) == 1]
    if len(uniq) < 2 or not once:
        raise SystemExit(
            'CONTROL: cannot plant the cancellation on this grid — need >=2 '
            'distinct slice d-values with >=1 of them unrepeated (got '
            f'{len(uniq)} distinct over {len(ds)} slices); extend the grid')
    cancel_d = once[len(once) // 2]
    vanish_d = next(d for d in uniq if d != cancel_d)
    return cancel_d, vanish_d


def run_control(store_names, outdir):
    """Gate (d): synthetic-truth control + negative controls, through the
    IDENTICAL pipeline code path (reconstruct_fn with a fake store/oracle)."""
    import random
    rnd = random.Random(20260801)
    ptab = firefly_primes()
    ds = sorted(slice_d(n) for n in store_names)
    # F* engineered with (i) a coefficient factor VANISHING at a grid slice d
    # and (ii) a COMMON eta-factor (eta+2) of P*,Q* exactly at another grid
    # slice d — the cancellation-point scenario. With cancel_d = u/v,
    # vanish_d = u'/v' (both picked from the grid by control_plants):
    #   P* = (v*d-u)*p1 + (eta+2)*p2,  Q* = (v*d-u)*q1 + (eta+2)*q2
    cancel_d, vanish_d = control_plants(ds)

    def rp(de_eta, de_d):
        return {e: fmpz_poly([rnd.randint(-99, 99) for _ in range(de_d)]
                             + [rnd.randint(1, 99)])
                for e in range(de_eta + 1)}

    p1, p2 = rp(4, 3), rp(4, 3)
    q1, q2 = rp(6, 2), rp(6, 2)
    # planted vanishing coefficient factor (v'*d - u')
    p2[2] = p2[2] * fmpz_poly([-vanish_d.numerator, vanish_d.denominator])
    ufac = fmpz_poly([-cancel_d.numerator, cancel_d.denominator])

    def add2(x, y):
        out = dict(x)
        for e, p in y.items():
            out[e] = out.get(e, fmpz_poly([0])) + p
        return {e: p for e, p in out.items() if p != 0}

    def mul_c(t):  # multiply by (eta + 2)
        return add2({e: 2 * p for e, p in t.items()},
                    {e + 1: p for e, p in t.items()})

    A = add2({e: ufac * p for e, p in p1.items()}, mul_c(p2))
    B = add2({e: ufac * p for e, p in q1.items()}, mul_c(q2))
    nd = max(p.degree() for p in A.values())
    dd = max(p.degree() for p in B.values())
    SN, SD = sorted(A), sorted(B)
    # coverage law for the synthetic fn itself: need_fit = nd+dd+5 fit slices
    # + >=3 more held-out + the cancellation slice (fit-excluded) = nd+dd+9
    if len(ds) < nd + dd + 9:
        raise SystemExit(
            f'CONTROL: grid has {len(ds)} slices; the synthetic control '
            f'(nd={nd}, dd={dd}) needs >= {nd + dd + 9} '
            f'(need_fit={nd + dd + 5} fit slices + 3 held-out + the planted '
            'cancellation slice) — extend the grid')
    cancel_name = f's{ds.index(cancel_d)}'

    class FakeStore:
        def __init__(self, corrupt=None):
            self.names = [f's{i}' for i in range(len(ds))]
            self.d = {f's{i}': ds[i] for i in range(len(ds))}
            self.corrupt = corrupt or {}

        def row(self, name, fn):
            dv = self.d[name]
            u, v = dv.numerator, dv.denominator
            Nv = {k: Fraction(*poly_eval_uv(
                [int(p[i]) for i in range(p.degree() + 1)], u, v))
                for k, p in A.items()}
            Dv = {j: Fraction(*poly_eval_uv(
                [int(p[i]) for i in range(p.degree() + 1)], u, v))
                for j, p in B.items()}
            # FireFly returns the REDUCED form: divide out any common factor
            # (this is what makes cancellation-point slices degree-dropped)
            FQ = lambda q: fmpq(q.numerator, q.denominator)  # noqa: E731
            Np = fmpq_poly([FQ(Nv.get(k, Fraction(0)))
                            for k in range(max(SN) + 1)])
            Dp = fmpq_poly([FQ(Dv.get(j, Fraction(0)))
                            for j in range(max(SD) + 1)])
            g = Np.gcd(Dp)
            if g.degree() > 0:
                Np, Dp = Np // g, Dp // g
            Nv = {i: Fraction(str(Np[i])) for i in range(Np.degree() + 1)
                  if Np[i] != 0}
            Dv = {i: Fraction(str(Dp[i])) for i in range(Dp.degree() + 1)
                  if Dp[i] != 0}
            # FireFly-style per-slice normalization: divide by den's lowest
            # nonzero coeff if any, else num's lowest nonzero
            c = None
            for j in sorted(Dv):
                if Dv[j] != 0:
                    c = Dv[j]
                    break
            if c is None:
                for k in sorted(Nv):
                    if Nv[k] != 0:
                        c = Nv[k]
                        break
            Nn = {str(k): f'{(x/c).numerator}/{(x/c).denominator}'
                  for k, x in Nv.items() if x != 0}
            Dn = {str(j): f'{(x/c).numerator}/{(x/c).denominator}'
                  for j, x in Dv.items() if x != 0}
            key = (name,)
            if key in self.corrupt:
                side, pos, delta = self.corrupt[key]
                tgt = Nn if side == 'N' else Dn
                pk = str(pos) if str(pos) in tgt else sorted(tgt, key=int)[0]
                a, b = parse_frac(tgt[pk])
                tgt[pk] = f'{a + delta * b}/{b}'
            return {'fn': fn, 'tag': 'SYNTH', 'done': 1, 'np': 1,
                    'N': Nn, 'D': Dn, 'validated': 'synthetic'}

    def make_oracle(mode, corrupt_residue=False):
        # bank normalization: scale by an arbitrary rational lambda
        lam = Fraction(rnd.randint(1, 10 ** 6), rnd.randint(1, 10 ** 6))
        gn = {}
        gd = {}
        for k, p in A.items():
            for e in range(p.degree() + 1):
                if int(p[e]):
                    q = lam * int(p[e])
                    gn[(e, k)] = (q.numerator, q.denominator)
        for j, p in B.items():
            for e in range(p.degree() + 1):
                if int(p[e]):
                    q = lam * int(p[e])
                    gd[(e, j)] = (q.numerator, q.denominator)
        if mode == 'done':
            return {'is_done': 1, 'combined_prime': ptab[0], 'shift': [0, 0],
                    'g_n': gn, 'g_d': gd, 'r_n': {}, 'r_d': {}}
        # active: commit everything as residues mod prod(primes[:8])
        M = 1
        for p in ptab[:8]:
            M *= p
        rn = {m: (a * pow(b, -1, M)) % M for m, (a, b) in gn.items()}
        rd = {m: (a * pow(b, -1, M)) % M for m, (a, b) in gd.items()}
        if corrupt_residue:
            m0 = sorted(rn)[0]
            rn[m0] = (rn[m0] + 12345) % M
        return {'is_done': 0, 'combined_prime': M, 'shift': [0, 0],
                'g_n': {}, 'g_d': {}, 'r_n': rn, 'r_d': rd}

    cov_row = {'fn': 0, 'tag': 'SYNTH', 'nd': nd, 'dd': dd,
               'slices': [f's{i}' for i in range(len(ds))]}
    report = []

    def proportional2(x, y):
        lam = None
        for k in set(x) | set(y):
            a, b = x.get(k, 0), y.get(k, 0)
            if b == 0:
                if a != 0:
                    return False
            else:
                l = Fraction(a, b)
                if lam is None:
                    lam = l
                elif l != lam:
                    return False
        return lam not in (None, 0)

    # positive control, bank-done oracle path; requires (i) BANKED, (ii) exact
    # recovery of the known truth, (iii) the planted cancellation slice
    # (d=cancel_d) auto-detected and excluded from fits
    st = FakeStore()
    r = reconstruct_fn(0, cov_row, st, lambda fn: make_oracle('done'), ptab)
    r0 = r
    ok = r['status'] == 'BANKED'
    if ok:
        Pk = {int(k): fmpz_poly([int(c) for c in v]) for k, v in r['N'].items()}
        Qj = {int(j): fmpz_poly([int(c) for c in v]) for j, v in r['D'].items()}
        lhs = conv2(_monomials(Pk), _monomials(B))
        rhs = conv2(_monomials(Qj), _monomials(A))
        ok = (proportional2(lhs, rhs)
              and r.get('cancel_slices') == [cancel_name]
              and cancel_name not in r['fit_slices'])
    report.append(('positive+bank-done+exact-recovery+cancel-detect', ok,
                   {k: v for k, v in r['gates'].items()}, r['status']))

    # positive control, bank-ACTIVE (residue) oracle path
    r = reconstruct_fn(0, cov_row, st, lambda fn: make_oracle('active'), ptab)
    report.append(('positive+bank-active-residues', r['status'] == 'BANKED',
                   dict(r['gates']), r['status']))

    # negative 1: corrupt one FIT (guard) slice coefficient -> must fail
    fitslice = r0['fit_slices'][2] if r0.get('fit_slices') else 's7'
    stc = FakeStore(corrupt={(fitslice,): ('N', 1, 1)})
    r = reconstruct_fn(0, cov_row, stc, lambda fn: make_oracle('done'), ptab)
    report.append(('negative:corrupt-fit-guard-slice',
                   r['status'] in ('FIT_FAIL', 'LOO_FAIL', 'HELDOUT_FAIL',
                                   'INTERNAL_FAIL', 'DEGCK_FAIL'),
                   dict(r['gates']), r['status']))

    # negative 2: corrupt one clean HELD-OUT slice -> heldout gate fails
    ho = next(s for s in r0['heldout_slices']
              if s not in r0.get('cancel_slices', [])
              and s not in r0.get('dropped_slices', []))
    stc = FakeStore(corrupt={(ho,): ('D', 2, 3)})
    r = reconstruct_fn(0, cov_row, stc, lambda fn: make_oracle('done'), ptab)
    report.append(('negative:corrupt-heldout-slice',
                   r['status'] == 'HELDOUT_FAIL', dict(r['gates']),
                   r['status']))

    # negative 2b: corrupt the CANCELLATION slice -> cross-mult check on the
    # excluded slice must still catch it (heldout gate)
    stc = FakeStore(corrupt={(cancel_name,): ('D', 1, 5)})
    r = reconstruct_fn(0, cov_row, stc, lambda fn: make_oracle('done'), ptab)
    report.append(('negative:corrupt-cancel-slice',
                   r['status'] == 'HELDOUT_FAIL', dict(r['gates']),
                   r['status']))

    # negative 3: corrupt one bank residue -> oracle gate fails
    r = reconstruct_fn(0, cov_row, st,
                       lambda fn: make_oracle('active', corrupt_residue=True),
                       ptab)
    report.append(('negative:corrupt-oracle-residue',
                   r['status'] == 'ORACLE_FAIL', dict(r['gates']),
                   r['status']))

    allok = all(x[1] for x in report)
    out = {'control': 'PASS' if allok else 'FAIL',
           'degrees': {'nd': nd, 'dd': dd, 'SN': SN, 'SD': SD,
                       'planted_vanishing':
                           f'p2[2] ~ ({vanish_d.denominator}d'
                           f'{-vanish_d.numerator:+d}) vanishes '
                           f'at d={vanish_d}',
                       'planted_cancellation':
                           f'P*,Q* share (eta+2) at d={cancel_d} '
                           f'({cancel_name})'},
           'cases': [{'name': n, 'ok': o, 'gates': g, 'status': s}
                     for n, o, g, s in report]}
    with open(os.path.join(outdir, 'control_result.json'), 'w') as f:
        json.dump(out, f, indent=1)
    for n, o, g, s in report:
        log(f'CONTROL {n}: {"OK" if o else "**FAIL**"} (status={s}, gates={g})')
    log(f'CONTROL VERDICT: {out["control"]}')
    return allok


# ---------------------------------------------------------------- workers
def worker_main(args_tuple):
    (wid, fns, cov_rows, slice_names, outdir, loo) = args_tuple
    os.nice(10)
    store = SliceStore(slice_names)
    # load pre-built indices
    with open(os.path.join(outdir, 'slice_index.json')) as f:
        raw = json.load(f)
    store.idx = {s: {int(k): v for k, v in d.items()} for s, d in raw.items()}
    ptab = firefly_primes()
    outp = os.path.join(outdir, f'shard_{wid:03d}.jsonl')
    done = 0
    with open(outp, 'w') as f:
        for fn in fns:
            try:
                r = reconstruct_fn(fn, cov_rows[fn], store, read_bank_state,
                                   ptab, loo=loo)
            except Exception as e:  # loud, never silent
                r = {'fn': fn, 'tag': cov_rows[fn]['tag'],
                     'status': 'EXCEPTION', 'error': repr(e)}
            f.write(json.dumps(r) + '\n')
            done += 1
            if done % 200 == 0:
                f.flush()
    return wid, len(fns)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=['control', 'probe3', 'index', 'run'])
    ap.add_argument('--fns', help='comma list or file of fn indices')
    ap.add_argument('--limit', type=int, default=0)
    ap.add_argument('--workers', type=int, default=24)
    ap.add_argument('--outdir', default=None,
                    help='output dir (was <BASE>/reconstruction; now explicit)')
    ap.add_argument('--batchdir', default=None)
    ap.add_argument('--no-loo', action='store_true')
    ap.add_argument('--skip-banked', default=None,
                    help='jsonl(.gz) of prior results; BANKED fns are skipped')
    # wiring: flag > --config > --base
    ap.add_argument('--config', default=None,
                    help='JSON with keys base/wave1/p2_file/bank_states/'
                         'ff_helper/coverage (CLI flags override)')
    ap.add_argument('--base', default=None, help='run root dir')
    ap.add_argument('--wave1', default=None, help='per-slice harvest dir')
    ap.add_argument('--p2-file', default=None, help='extra single-slice harvest jsonl')
    ap.add_argument('--bank-states', default=None,
                    help='never-in-fit oracle bank ff_save/states dir')
    ap.add_argument('--ff-helper', default=None,
                    help='FireFly ReconstHelper.cpp (prime-table source)')
    ap.add_argument('--coverage', default=None, help='per-fn coverage jsonl.gz')
    args = ap.parse_args()

    resolve_config(args)
    outdir = args.batchdir or args.outdir
    if not outdir:
        ap.error('need --outdir or --batchdir (was <BASE>/reconstruction)')
    os.makedirs(outdir, exist_ok=True)

    cov = load_coverage()
    if args.cmd in ('probe3', 'run'):   # fail-closed BEFORE any worker launch
        _require(BANK_STATES, 'BANK_STATES', '--bank-states')
        _require(FF_HELPER, 'FF_HELPER', '--ff-helper')
    slice_names = sorted({s for r in cov.values() for s in r.get('slices', [])},
                         key=lambda n: slice_d(n))
    log(f'coverage: {len(cov)} fns, {len(slice_names)} slices')

    if args.cmd == 'control':
        ok = run_control(slice_names, outdir)
        sys.exit(0 if ok else 3)

    if args.cmd == 'index':
        store = SliceStore(slice_names)
        store.build_indices()
        with open(os.path.join(outdir, 'slice_index.json'), 'w') as f:
            json.dump({s: {str(k): v for k, v in d.items()}
                       for s, d in store.idx.items()}, f)
        log('index written')
        sys.exit(0)

    # eligible = covered fns
    eligible = [fn for fn, r in sorted(cov.items())
                if not r.get('excluded_p2open')
                and r.get('cov', 0) >= r.get('need_total', 10 ** 9)]
    if args.skip_banked:
        op = gzip.open if args.skip_banked.endswith('.gz') else open
        prior = {json.loads(l)['fn'] for l in op(args.skip_banked, 'rt')
                 if '"BANKED"' in l}
        eligible = [f for f in eligible if f not in prior]
        log(f'skip-banked: {len(prior)} prior; {len(eligible)} remain')
    if args.fns:
        if os.path.exists(args.fns):
            want = {int(x) for x in open(args.fns).read().split()}
        else:
            want = {int(x) for x in args.fns.split(',')}
        eligible = [f for f in eligible if f in want]
    if args.limit:
        eligible = eligible[:args.limit]
    log(f'eligible fns this run: {len(eligible)}')

    if args.cmd == 'probe3':
        # normalization-convention probe: report which position FireFly set
        # to 1/1 per slice, per fn, and fit end-to-end
        store = SliceStore(slice_names)
        with open(os.path.join(outdir, 'slice_index.json')) as f:
            raw = json.load(f)
        store.idx = {s: {int(k): v for k, v in d.items()}
                     for s, d in raw.items()}
        ptab = firefly_primes()
        for fn in eligible[:3] if not args.fns else eligible:
            r = cov[fn]
            ones = {}
            for s in r['slices']:
                row = store.row(s, fn)
                if row is None:
                    continue
                one = [(side, k) for side in 'ND'
                       for k, v in row[side].items() if v == '1/1']
                ones[s] = one
            log(f'fn {fn} (nd={r["nd"]} dd={r["dd"]}): FireFly-normalized-to-1 '
                f'positions per slice: {ones}')
            rr = reconstruct_fn(fn, r, store, read_bank_state, ptab)
            log(f'fn {fn}: status={rr["status"]} gates={rr.get("gates")}')
        sys.exit(0)

    # run
    import multiprocessing as mp
    nw = args.workers
    shards = [eligible[i::nw] for i in range(nw)]
    work = [(i, shards[i], {f: cov[f] for f in shards[i]}, slice_names,
             outdir, not args.no_loo) for i in range(nw) if shards[i]]
    log(f'launching {len(work)} workers (nice 10)')
    with mp.Pool(len(work)) as pool:
        for wid, n in pool.imap_unordered(worker_main, work):
            log(f'worker {wid} done ({n} fns)')
    log('all workers done')


if __name__ == '__main__':
    main()
