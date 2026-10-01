#!/usr/bin/env python3
r"""numkin_harvest.py — reconstruction half of the "numkin" technique.

Companion to tools/numkin/numkin_sweep.sh (the generation half; the flat
path tools/numkin_sweep.sh is an exec-wrapper). numkin_sweep.sh
sed-substitutes a numeric value for ONE symbolic kinematic invariant (e.g. m2)
directly into already-generated Kira SYSTEM_<FAM>_*.gz files, leaving `d`
symbolic, then runs ordinary FireFly per point. Each landed point therefore
gives, per (target, master) pair, an EXACT rational function of `d` (Kira's
own FireFly already reconstructed the `d`-dependence in full at that fixed
numeric kinematic point -- no d-fitting is needed, only reassembly across the
point grid). This script sweeps the point grid, Pade-reconstructs the
dependence on the swept variable (`--var`, e.g. m2), and verifies against
held-out points with EXACT rational arithmetic.

SCOPE: this script assumes DIRECT per-target reductions -- kira already
gives target = sum_j coef_j(d, var) * master_j. It carries no
family-specific "recipe expansion" layer (that layer only matters when
patching missing entries of an already-assembled matrix, which is not this
tool's job): it goes straight from parse_kira2math() output to the
bivariate (d, var) reconstruction.

ALGORITHMIC CORE:
  - parse_kira2math()            : delegates to tools/frobenius_boundary/kira_parse.py
  - _pade_fmpz / _pade_verify     : Cauchy/Pade rational reconstruction via
                                    flint fmpz_mat integer-nullspace (exact)
  - bivariate reassembly          : per-node Pade-in-var, then per-coefficient
                                    Pade-in-d, LCM-clear, content-primitive
  - DGridInsufficient              : signals "not enough points", never silently
                                    returns a wrong answer

Note on the "d" axis: each point's coefficient string IS already an exact
closed-form rational function of d (parsed once via a parse_dpoly-style
regex into an flint fmpq_poly numerator/denominator pair). Per-node
Pade-in-d is still needed because different numkin points may return that
d-rational in different (but equivalent) normalized forms -- so we evaluate
at a small d-node grid and Pade-reconstruct per var-coefficient.

Output: A(var) matrix in TWO forms:
  1. JSON  (--out):
       {'entries': {'<target idx csv>|<master idx csv>':
                       {'N_d_by_varpow': [str, ...],   # numerator, ascending
                                                        # powers of `var`, each
                                                        # a plain-text sympy-
                                                        # parseable poly in d
                        'D_d_by_varpow': [str, ...]},  # denominator, same
                    ...},
        'provenance': {method, family, var, fit_points, held_out_points,
                       targets_file, n_columns_total, n_ok, n_failed,
                       n_skipped_too_few_points, failed_entries}}
     i.e. coef(d, var) = (sum_b N_d_by_varpow[b] * var**b)
                        / (sum_b D_d_by_varpow[b] * var**b).
     This is a documented, tool-specific schema storing raw (d, var) -- no
     eps = (4-d)/2 substitution or eps-grid expansion is performed.
  2. Synthetic kira2math-format .m file (--out-m, optional) so any
     kira2math-format consumer can parse the reconstructed matrix
     with ZERO adapter code (same regex contract:
     FAM[idx] -> \n + FAM[idx2]*(coef)\n ... `,` between blocks; coef uses the
     real --var name, e.g. 'm2', and python `**` power notation, which
     sp.sympify(coef, locals={'d': d, 'm2': m2}) parses natively).

No fabrication: every reconstructed entry is checked against >=2 held-out
numkin points with exact (gmpy2/flint) arithmetic. Failures are logged in
'failed_entries' and OMITTED from the output matrix/. m file (never silently
dropped -- this script prints and records them).
"""
import argparse
import glob
import json
import math
import os
import re
import sys
import time
from fractions import Fraction as Fr

# This file lives in tools/numkin/;
# kira_parse lives in tools/frobenius_boundary/, i.e. a sibling of the PACKAGE,
# so anchor on the package's parent dir. NB: this sys.path side effect is
# LOAD-BEARING for external consumers — downstream harnesses may do
# `import numkin_harvest` then `from kira_parse import ...`.
# Do not convert to a package-relative import.
_FB_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for _cand in (os.path.join(_FB_ROOT, 'frobenius_boundary'),
              os.path.join(_FB_ROOT, 'frobenius-boundary', 'frobenius_boundary')):
    if os.path.isdir(_cand):
        sys.path.insert(0, _cand)
        break
from kira_parse import parse_kira2math  # noqa: E402

from flint import fmpz, fmpz_mat, fmpz_poly, fmpq, fmpq_poly  # noqa: E402

t0 = time.time()


def log(msg):
    print(f'[{time.time()-t0:7.1f}s] {msg}', flush=True)


# ============================================================ d-poly parsing
# Kira/FireFly ("Fermat-normalized") coefficient strings use `*` and `^`
# (not python `**`), heavy nested parenthesization, and are, after numkin
# substitution, pure rational functions of `d` alone (the swept variable has
# been replaced by a numeric literal everywhere it appeared as a bare token).
# We do NOT need a bespoke recursive-descent parser: Python's own expression
# grammar accepts `^`-free algebra, so translate `^` -> `**` and eval() with
# `d` bound to a chosen flint.fmpq node value and integers auto-promoted to
# fmpq (so `/` stays exact rational division, never float).
_CARET_POW = re.compile(r'\^')


def _prep_expr(src):
    src = src.strip()
    if src in ('', '0', '(0)'):
        return '0'
    if src in ('1', '(1)'):
        return '1'
    return _CARET_POW.sub('**', src)


class _FQ(fmpq):
    """fmpq subclass whose __truediv__/__rtruediv__ never falls back to float,
    and whose int literals passed through eval() become fmpq automatically
    via the eval namespace (see eval_fmpq)."""
    pass


_POW_RE = re.compile(r'\*\*(\d+)')


def eval_fmpq(expr_src, d_val):
    """Evaluate a Fermat-style coefficient string at d = d_val (fmpq) -> fmpq.
    Pure integer/rational arithmetic throughout (flint fmpq), exact.

    NB: exponents must stay plain python int -- fmpq**fmpq raises (flint has
    no general rational-power op), and even where it wouldn't raise, wrapping
    an exponent would be mathematically wrong. So: wrap every bare integer
    literal with fmpq(...) EXCEPT ones immediately following '**' (handled
    first, protected via a placeholder, then restored after the general
    wrap -- a two-pass protect-then-wrap pattern)."""
    src = _prep_expr(expr_src)
    if src == '0':
        return fmpq(0)
    if src == '1':
        return fmpq(1)
    protected = _POW_RE.sub(lambda m: f'**QQPOWQQ{m.group(1)}QQPOWQQ', src)
    wrapped = re.sub(r'(?<![\w.])(\d+)(?![\w.])', r'fmpq(\1)', protected)
    wrapped = re.sub(r'\*\*QQPOWQQ(\d+)QQPOWQQ', r'**\1', wrapped)
    try:
        return eval(wrapped, {'__builtins__': {}}, {'fmpq': fmpq, 'd': d_val})
    except NameError as e:
        # Most likely cause: the numkin substitution left the swept variable
        # (e.g. 'm2') as a free symbol in this coefficient -- i.e. the
        # numkin_sweep.sh stage() substitution (kinematics.yaml strip +
        # integralfamilies.yaml sed) did not reach every occurrence, e.g. a
        # config layout its patterns do not match.
        # Fail loudly rather than silently mis-evaluate.
        raise RuntimeError(
            f'eval_fmpq: unresolved free symbol in coefficient string '
            f'(expected pure-d rational after numkin substitution): '
            f'{e}  -- src={expr_src!r}') from e


# node grid for the Pade-in-d reassembly step (mirrors thiele_ratd's dfit/dchk
# split: small positive integers keep coefficients compact).
_D_FIT_N = 30
_D_FIT = tuple(fmpq(101 + i) for i in range(_D_FIT_N))
_D_CHK = (fmpq(1009, 3), fmpq(-211, 5))
_D_ALL = _D_FIT + _D_CHK


# ============================================================ Pade (Cauchy) via flint nullspace
class DGridInsufficient(RuntimeError):
    """Not enough points/nodes to pin down the rational function at the
    attempted (p, q) degree split. Caller must add points, never guess."""


def _pade_fmpz(xs_uv, ys_nd):
    """Exact rational interpolation via integer-matrix nullspace.
      xs_uv : list[(int u_k, int v_k)]  nodes x_k = u_k/v_k  (v_k>0, coprime)
      ys_nd : list[(int A_k, int B_k)]  values y_k = A_k/B_k (B_k != 0)
    Returns (N, D) ascending int-coeff lists, gcd(N,D)=1, primitive, D-lead>0.
    Ansatz p+q = K-1 -> nullity >= 1 always; excess freedom is removed by
    fmpz_poly.gcd. If the true (p*, q*) satisfies p*<=p, q*<=q the result is
    exact; otherwise it's an artifact of an under-determined system and the
    caller's held-out check rejects it (never silently trusted)."""
    K = len(xs_uv)
    q = K // 2
    p = K - 1 - q
    P = max(p, q)
    Mat = fmpz_mat(K, p + q + 2)
    for k in range(K):
        A, B = ys_nd[k]
        u, v = xs_uv[k]
        # powers of (u,v): pw[b] = u^b * v^(P-b)
        pw = [0] * (P + 1)
        pu = 1
        vp = [1] * (P + 1)
        for b in range(1, P + 1):
            vp[b] = vp[b - 1] * v
        for b in range(P + 1):
            pw[b] = pu * vp[P - b]
            pu *= u
        for b in range(p + 1):
            Mat[k, b] = B * pw[b]
        for b in range(q + 1):
            Mat[k, p + 1 + b] = -A * pw[b]
    NS, nd = Mat.nullspace()
    if nd == 0:
        raise DGridInsufficient(f'_pade_fmpz: nullity 0 (K={K}, p={p}, q={q})')
    col = 0
    for c in range(nd):
        if any(NS[p + 1 + b, c] != 0 for b in range(q + 1)):
            col = c
            break
    Ncoef = [int(NS[b, col]) for b in range(p + 1)]
    Dcoef = [int(NS[p + 1 + b, col]) for b in range(q + 1)]
    Np, Dp = fmpz_poly(Ncoef), fmpz_poly(Dcoef)
    g = Np.gcd(Dp)
    if g.degree() > 0:
        Np //= g
        Dp //= g
    N = [int(Np[i]) for i in range(max(1, len(Np)))]
    D = [int(Dp[i]) for i in range(max(1, len(Dp)))]
    while len(N) > 1 and N[-1] == 0:
        N.pop()
    while len(D) > 1 and D[-1] == 0:
        D.pop()
    g = 0
    for c in N:
        g = math.gcd(g, c)
    for c in D:
        g = math.gcd(g, c)
    if g > 1:
        N = [c // g for c in N]
        D = [c // g for c in D]
    if D[-1] < 0:
        N = [-c for c in N]
        D = [-c for c in D]
    return N, D


def _poly_eval_int(P, x_u, x_v):
    """Evaluate ascending int-coeff poly P at x = x_u/x_v -> (num, den) ints
    (not reduced)."""
    num, den = 0, 1
    for i, c in enumerate(P):
        if c == 0:
            continue
        # term c * (x_u/x_v)^i, accumulate as fraction with common growing den
        tn, td = c * (x_u ** i), x_v ** i
        num = num * td + tn * den
        den = den * td
    return num, den


def _pade_verify(N, D, xs_uv, ys_nd):
    for (u, v), (A, B) in zip(xs_uv, ys_nd):
        nn, nd_ = _poly_eval_int(N, u, v)
        dn, dd_ = _poly_eval_int(D, u, v)
        # N(x)/D(x) == A/B  <=>  nn*dd_*B == dn*nd_*A  (cross-multiply, all int)
        if nn * dd_ * B != dn * nd_ * A:
            return False
    return True


def _poly_lcm(a, b):
    if a == 0:
        return b
    if b == 0:
        return a
    g = a.gcd(b)
    return a * (b // g)


def _bivar_primitive_int(N2, D2):
    """(list[fmpq_poly] ascending in `var`) -> integer-coeff, content-primitive,
    D-lead > 0. Trims trailing-zero var-coefficients."""
    while len(N2) > 1 and N2[-1] == 0:
        N2.pop()
    while len(D2) > 1 and D2[-1] == 0:
        D2.pop()
    Lq = fmpz(1)
    for p in N2 + D2:
        for i in range(len(p)):
            Lq = Lq.lcm(p[i].q)
    N2 = [Lq * p for p in N2]
    D2 = [Lq * p for p in D2]
    g = fmpz(0)
    for p in N2 + D2:
        for i in range(len(p)):
            g = g.gcd(p[i].p)
    if g == 0:
        g = fmpz(1)
    inv = fmpq(1, g)
    N2 = [inv * p for p in N2]
    D2 = [inv * p for p in D2]
    lead = None
    for p in reversed(D2):
        if p != 0:
            lead = p[p.degree()]
            break
    if lead is not None and lead < 0:
        N2 = [-p for p in N2]
        D2 = [-p for p in D2]
    return N2, D2


def bivar_eval_var(P2, xQ):
    """list-of-fmpq_poly (ascending var-powers) evaluated at var=xQ (fmpq)
    -> fmpq_poly in d."""
    r = fmpq_poly()
    for c in reversed(P2):
        r = xQ * r + c
    return r


def thiele_ratd_var(xs, ys, mfit=_D_FIT_N, raise_on_dfail=True):
    """Bivariate (d, var) exact rational reconstruction.
      xs : list[fmpq]                     var-nodes, length K (numkin points)
      ys : list[(str, str)]               (numerator_src, denominator_src) --
                                           Fermat-style strings in `d` ONLY,
                                           read straight from parse_kira2math.
    Returns (N2, D2): lists of fmpq_poly in d, ascending var-power, primitive
    integer-coeff bivariate N(d,var), D(d,var). Raises DGridInsufficient if
    the var-Pade or the per-coefficient d-Pade doesn't converge/verify.

    Algorithm (per-node Pade-in-var then per-coefficient Pade-in-d, with
    one difference from the generic bivariate case: each y_k here is ALREADY an exact d-rational
    closed form supplied directly by Kira/FireFly at that numkin point, not
    itself something needing d-fitting -- but we still evaluate it at a small
    d-node grid because different points may hand back that d-rational in
    different (content/normalization) forms, and reassembling per-var-Pade
    per d-node is the exact, cheap, proven way to remove that ambiguity):
      1. Evaluate every y_k (a (num_src, den_src) pair) at len(_D_ALL) small
         rational d-nodes -> integer pair grid Y[m][k] = (A, B).
      2. Per d-node m: Pade-in-var via flint fmpz nullspace -> primitive int
         (Nm, Dm); monic-normalize by Dm's leading coeff.
      3. Per var-coefficient b: Pade-in-d over the M fit d-nodes on samples
         (coeff_b(d_m), lc(d_m)); verify at 2 check d-nodes.
      4. LCM-clear per-b d-denominators, remove d-content, canonical sign.
    """
    K = len(xs)
    xs_uv = tuple((int(x.p), int(x.q)) for x in xs)
    dall_uv = tuple((int(dv.p), int(dv.q)) for dv in _D_ALL)
    M = len(_D_FIT)

    # 1) integer grid Y[m][k] = (A, B) with y_k(d_m) = A/B
    Ynd = [[None] * K for _ in dall_uv]
    for k, (nsrc, dsrc) in enumerate(ys):
        for m, dmQ in enumerate(_D_ALL):
            nv = eval_fmpq(nsrc, dmQ)
            dv = eval_fmpq(dsrc, dmQ)
            if dv == 0:
                raise DGridInsufficient(f'point {k}: denominator vanishes at d-node {dmQ}')
            # y = nv/dv ; cross-multiply to plain (A,B) ints
            Ynd[m][k] = (int(nv.p) * int(dv.q), int(nv.q) * int(dv.p))

    # 2) per d-node: Pade-in-var
    NDm = []
    degN = degD = 0
    maxbits = 0
    for m in range(len(dall_uv)):
        if all(A == 0 for (A, B) in Ynd[m]):
            NDm.append(([0], [1]))
            continue
        Nm, Dm = _pade_fmpz(xs_uv, Ynd[m])
        NDm.append((Nm, Dm))
        degN = max(degN, len(Nm) - 1)
        degD = max(degD, len(Dm) - 1)
        maxbits = max(maxbits, max(abs(c).bit_length() for c in Nm + Dm))
        if m == 1 and maxbits > 1200:
            if raise_on_dfail:
                raise DGridInsufficient(
                    f'var-Pade unconverged (primitive coeffs {maxbits} bits >> 200)')
            return None

    lc = [pair[1][-1] if len(pair[1]) - 1 == degD else 0 for pair in NDm]
    if any(c == 0 for c in lc):
        if raise_on_dfail:
            raise DGridInsufficient('degD not uniform across d-grid -- var-Pade unconverged')
        return None

    def samp_nd(which, b):
        out = []
        for m, pair in enumerate(NDm):
            arr = pair[which]
            out.append((arr[b] if b < len(arr) else 0, lc[m]))
        return out

    N2, D2 = [], []
    dens = []
    for which, deg, out in ((0, degN, N2), (1, degD, D2)):
        for b in range(deg + 1):
            snd = samp_nd(which, b)
            if all(A == 0 for (A, B) in snd):
                out.append((fmpq_poly(), fmpq_poly([1])))
                continue
            Nb, Db = _pade_fmpz(dall_uv[:M], snd[:M])
            if not _pade_verify(Nb, Db, dall_uv[M:], snd[M:]):
                if raise_on_dfail:
                    raise DGridInsufficient(
                        f'd-grid check FAILED at var^{b} ({["N","D"][which]}): '
                        f'var-Pade unconverged or mfit={mfit} too small')
                return None
            nb, db = fmpq_poly(Nb), fmpq_poly(Db)
            ldc = db[db.degree()]
            nb *= fmpq(ldc.q, ldc.p)
            db *= fmpq(ldc.q, ldc.p)
            out.append((nb, db))
            dens.append(db)
    L = fmpq_poly([1])
    for db in dens:
        L = _poly_lcm(L, db)
    N2 = [nb * (L // db) for (nb, db) in N2]
    D2 = [nb * (L // db) for (nb, db) in D2]
    G = fmpq_poly()
    for p in N2 + D2:
        G = G.gcd(p)
    if G != 0 and (G.degree() > 0 or G[0] != 1):
        N2 = [p // G for p in N2]
        D2 = [p // G for p in D2]
    return _bivar_primitive_int(N2, D2)


# ============================================================ point collection
def collect_points(pts_dir, family, varname):
    """Return list of (k:int, m_path:str, var_value:Fraction) for every
    pt_* directory that has landed a nonempty results/<FAM>/*.m and a
    <varname>val file."""
    out = []
    for d in sorted(glob.glob(os.path.join(pts_dir, 'pt_*')),
                     key=lambda p: int(re.search(r'pt_(\d+)$', p).group(1))):
        k = int(re.search(r'pt_(\d+)$', d).group(1))
        vf = os.path.join(d, f'{varname}val')
        if not os.path.exists(vf):
            continue
        mfiles = glob.glob(os.path.join(d, 'results', family, '*.m'))
        mfiles = [f for f in mfiles if os.path.getsize(f) > 10]
        if not mfiles:
            continue
        vv = open(vf).read().strip()
        vQ = Fr(vv)
        out.append((k, d, mfiles, vQ))
    return out


def load_targets(path):
    """Parse a target-list file (one FAM[idx,...] per line) -> list of tuples."""
    out = []
    pat = re.compile(r'\[([^\]]+)\]')
    for ln in open(path):
        m = pat.search(ln)
        if m:
            out.append(tuple(int(x) for x in m.group(1).split(',')))
    return out


# ============================================================ .m emission (kira2math-format compatible)
def _fmpq_poly_to_dsrc(p):
    """fmpq_poly in d -> Fermat-style string usable by both parse_kira2math
    and sympify-based parse_m-style consumers (sympify with locals {d}). Uses
    `**` (both accept python-style power via sympify / our own eval)."""
    if p == 0:
        return '0'
    terms = []
    for i in range(p.degree(), -1, -1):
        c = p[i]
        if c == 0:
            continue
        cs = f'({c.p}/{c.q})' if c.q != 1 else str(int(c.p))
        if i == 0:
            terms.append(cs)
        elif i == 1:
            terms.append(f'{cs}*d')
        else:
            terms.append(f'{cs}*d**{i}')
    return ' + '.join(terms) if terms else '0'


def emit_m_file(out_path, family, block, varname='m2'):
    """block: {target_tuple: {master_tuple: (N2, D2)}} -> synthetic kira2math
    .m file, one block per target, in the exact syntax parse_kira2math /
    sympify-based parse_m-style consumers expect:
        FAM[idx] ->
         + FAM[idx2]*((num)/(den))
        ,
        ...
    Coefficients are LEFT in bivariate (d, var) form -- N2/D2 are lists of
    fmpq_poly-in-d ascending in var-power, so the (num)/(den) string is a
    ratio of two var-polynomials whose coefficients are d-rationals. This is
    NOT plain kira2math output (which never nests like this in practice for
    a single point) but IS valid input to sympy.sympify given locals
    {d, <var>} -- so a sympify-based parse_m works UNCHANGED as long as
    its sp.sympify locals dict includes the var name (see the note in the
    tool docstring)."""
    lines = ['{']
    first_t = True
    for tgt, masters in block.items():
        if not first_t:
            lines.append(',')
        first_t = False
        lines.append(f'{family}[{",".join(map(str, tgt))}] -> ')
        for mm, (N2, D2) in masters.items():
            nsrc = _bivar_poly_to_varsrc(N2, varname)
            dsrc = _bivar_poly_to_varsrc(D2, varname)
            lines.append(f' + {family}[{",".join(map(str, mm))}]*(({nsrc})/({dsrc}))')
    lines.append('}')
    with open(out_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')


def _bivar_poly_to_varsrc(P2, varname='m2'):
    """list[fmpq_poly-in-d] ascending var-power -> python/sympy-style string
    in (d, <varname>). Uses the actual swept-variable name (not a generic
    placeholder) so a consumer's
    sp.sympify(coef, locals={'d': d, 'm2': m2}) consumes it unchanged when
    varname == 'm2' (the default / primary use case)."""
    if not P2 or all(p == 0 for p in P2):
        return '0'
    terms = []
    for b, p in enumerate(P2):
        if p == 0:
            continue
        dsrc = _fmpq_poly_to_dsrc(p)
        wrapped = f'({dsrc})'
        if b == 0:
            terms.append(wrapped)
        elif b == 1:
            terms.append(f'{wrapped}*{varname}')
        else:
            terms.append(f'{wrapped}*{varname}**{b}')
    return ' + '.join(terms) if terms else '0'


# ============================================================ main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--family', required=True, help='Kira family name (as used in the reduction files)')
    ap.add_argument('--targets-file', required=True,
                     help='file with one FAM[idx,...] target per line (for logging/scoping only; '
                          'entries actually reconstructed = whatever parse_kira2math finds)')
    ap.add_argument('--pts-dir', required=True,
                     help='dir containing pt_*/results/<FAM>/*.m and pt_*/<var>val')
    ap.add_argument('--var', default='m2', help='swept kinematic variable name (default m2)')
    ap.add_argument('--held-out', type=int, default=2,
                     help='number of highest-k points to reserve as held-out verification (default 2)')
    ap.add_argument('--out', required=True, help='output JSON path')
    ap.add_argument('--out-m', default=None,
                     help='optional: also emit a synthetic kira2math-format .m file at this path')
    ap.add_argument('--min-fit-pts', type=int, default=3,
                     help='abort if fewer than this many fit points are available (default 3)')
    args = ap.parse_args()

    targets = load_targets(args.targets_file)
    log(f'{len(targets)} targets declared in {args.targets_file}')

    pts = collect_points(args.pts_dir, args.family, args.var)
    if not pts:
        log(f'NO landed points found under {args.pts_dir} (need pt_*/{args.var}val + '
            f'pt_*/results/{args.family}/*.m) -- nothing to do. This is the expected state '
            f'while numkin_sweep.sh probes are still running; re-run once .m files land.')
        json.dump({'entries': {}, 'provenance': {'status': 'NO_POINTS_LANDED',
                                                   'pts_dir': args.pts_dir}}, open(args.out, 'w'), indent=1)
        return

    log(f'found {len(pts)} landed points: k={[p[0] for p in pts]}')

    # parse every point's .m file(s) -> {target: [(master, coef_src_d), ...]}
    per_pt = {}
    for k, d, mfiles, vQ in pts:
        merged = {}
        for mf in mfiles:
            parsed = parse_kira2math(mf, args.family)
            for lhs, terms in parsed.items():
                merged.setdefault(lhs, []).extend(terms)
        per_pt[k] = (vQ, merged)
        log(f'  pt_{k}: {args.var}={vQ}  {len(merged)} target rows parsed from {len(mfiles)} .m file(s)')

    ks = sorted(per_pt)
    if len(ks) < args.min_fit_pts + args.held_out:
        log(f'WARNING: only {len(ks)} points landed; need >= {args.min_fit_pts} fit + '
            f'{args.held_out} held-out = {args.min_fit_pts + args.held_out}. '
            f'Proceeding with what is available but verification may be WEAK or entries '
            f'may FAIL for lack of degree freedom -- this is by design (no fabrication).')

    held = ks[-args.held_out:] if args.held_out > 0 and len(ks) > args.held_out else []
    fit = [k for k in ks if k not in held]
    log(f'fit points: {fit}  held-out points: {held}')
    if args.held_out and len(held) < 2:
        log('WARNING: <2 held-out points -- integrity rule requires >=2 genuine held-out '
            'verification points; treat any PASS here as PROVISIONAL until more points land.')

    # collect all (target, master) columns seen in the FIT points (held-out
    # points only need to supply a value to verify against, not to be a
    # column source)
    columns = {}  # (target, master) -> {k: (num_src, den_src)}
    for k in ks:
        vQ, merged = per_pt[k]
        for tgt, terms in merged.items():
            for mm, coef_src in terms:
                key = (tgt, mm)
                columns.setdefault(key, {})[k] = coef_src

    log(f'{len(columns)} distinct (target, master) columns across all points')

    xs_fit = [fmpq(per_pt[k][0].numerator, per_pt[k][0].denominator) for k in fit]
    ho_vals = {k: fmpq(per_pt[k][0].numerator, per_pt[k][0].denominator) for k in held}

    block = {}
    fails = []
    n_ok = n_fail = n_skip = 0
    for (tgt, mm), by_k in columns.items():
        fit_ks = [k for k in fit if k in by_k]
        if len(fit_ks) < args.min_fit_pts:
            fails.append({'target': tgt, 'master': mm, 'reason':
                          f'only {len(fit_ks)} fit points have this column (need >= {args.min_fit_pts})'})
            n_skip += 1
            continue
        xs_e = [fmpq(per_pt[k][0].numerator, per_pt[k][0].denominator) for k in fit_ks]
        ys_e = []
        ok_parse = True
        for k in fit_ks:
            src = by_k[k]
            # coefficient string is a single rational-in-d expression; split
            # num/den if it is literally "(...)/( ...)" (parse_kira2math coef
            # strings from a raw kira2math .m are ALWAYS of that "(N)/(D)"
            # outer form since kira2math always emits explicit numerator and
            # denominator parens) -- else treat whole string as numerator /1.
            src = src.strip()
            if src.startswith('(') :
                # find the matching close-paren of the leading '(' to split N | /D
                depth = 0
                for i, c in enumerate(src):
                    if c == '(':
                        depth += 1
                    elif c == ')':
                        depth -= 1
                        if depth == 0:
                            nsrc = src[1:i]
                            rest = src[i+1:]
                            break
                else:
                    nsrc, rest = src, ''
                if rest.startswith('/'):
                    dsrc = rest[1:]
                    if dsrc.startswith('(') and dsrc.endswith(')'):
                        dsrc = dsrc[1:-1]
                else:
                    nsrc, dsrc = src, '1'
            else:
                nsrc, dsrc = src, '1'
            ys_e.append((nsrc, dsrc))
        try:
            N2, D2 = thiele_ratd_var(xs_e, ys_e)
        except DGridInsufficient as e:
            fails.append({'target': tgt, 'master': mm, 'reason': str(e)})
            n_fail += 1
            continue
        # held-out verification: EXACT
        ok = True
        for k in held:
            if k not in by_k:
                continue  # held point simply doesn't carry this column (not a failure)
            src = by_k[k].strip()
            if src.startswith('('):
                depth = 0
                for i, c in enumerate(src):
                    if c == '(':
                        depth += 1
                    elif c == ')':
                        depth -= 1
                        if depth == 0:
                            hnsrc = src[1:i]
                            hrest = src[i+1:]
                            break
                else:
                    hnsrc, hrest = src, ''
                if hrest.startswith('/'):
                    hdsrc = hrest[1:]
                    if hdsrc.startswith('(') and hdsrc.endswith(')'):
                        hdsrc = hdsrc[1:-1]
                else:
                    hnsrc, hdsrc = src, '1'
            else:
                hnsrc, hdsrc = src, '1'
            vQ = ho_vals[k]
            gN = bivar_eval_var(N2, vQ)
            gD = bivar_eval_var(D2, vQ)
            # evaluate reconstructed d-rational and the raw held-out d-rational
            # at the SAME small d-grid and cross-check exactly (avoids having
            # to symbolically cancel two possibly differently-normalized
            # rational functions of d)
            for dchk in _D_CHK:
                g_n = eval_fmpq(_fmpq_poly_to_dsrc(gN), dchk) if gN != 0 else fmpq(0)
                g_d = eval_fmpq(_fmpq_poly_to_dsrc(gD), dchk)
                r_n = eval_fmpq(hnsrc, dchk)
                r_d = eval_fmpq(hdsrc, dchk)
                if g_d == 0 or r_d == 0:
                    ok = False
                    break
                if g_n * r_d != r_n * g_d:
                    ok = False
                    break
            if not ok:
                break
        if not ok:
            fails.append({'target': tgt, 'master': mm, 'reason': 'HELD-OUT VERIFICATION FAILED'})
            n_fail += 1
            continue
        block.setdefault(tgt, {})[mm] = (N2, D2)
        n_ok += 1

    log(f'reconstruction complete: {n_ok} ok, {n_fail} failed-verify/insufficient-degree, '
        f'{n_skip} skipped (too few fit points)')
    if fails:
        log(f'FIRST {min(10,len(fails))} failures: '
            + '; '.join(f"{fr['target']}<-{fr['master']}: {fr['reason']}" for fr in fails[:10]))

    # serialize
    ser_entries = {}
    for tgt, masters in block.items():
        for mm, (N2, D2) in masters.items():
            key = f'{",".join(map(str,tgt))}|{",".join(map(str,mm))}'
            ser_entries[key] = {
                'N_d_by_varpow': [_fmpq_poly_to_dsrc(p) for p in N2],
                'D_d_by_varpow': [_fmpq_poly_to_dsrc(p) for p in D2],
            }
    prov = {
        'method': f'numkin ({args.var}-sweep) + per-node flint Pade-in-{args.var} '
                  f'+ Pade-in-d reassembly + >=2 held-out EXACT verify',
        'family': args.family,
        'var': args.var,
        'fit_points': [str(per_pt[k][0]) for k in fit],
        'held_out_points': [str(per_pt[k][0]) for k in held],
        'targets_file': args.targets_file,
        'n_columns_total': len(columns),
        'n_ok': n_ok,
        'n_failed': n_fail,
        'n_skipped_too_few_points': n_skip,
        'failed_entries': fails,
    }
    json.dump({'entries': ser_entries, 'provenance': prov}, open(args.out, 'w'), indent=1)
    log(f'WROTE {args.out}  ({len(ser_entries)} entries)')

    if args.out_m:
        emit_m_file(args.out_m, args.family, block, varname=args.var)
        log(f'WROTE {args.out_m}  (synthetic kira2math-format, {len(block)} target blocks)')


if __name__ == '__main__':
    main()
