"""boxwalk.syz — tangency-module generators: Singular `syz` driver + parser
+ exact flint certification; sampled-scan fallback when Singular times out.

Module: kernel of R^{n+k} -> R^k,
  v = (a_1..a_n, b_1..b_k) |-> sum_e a_e*X_e*dP_j/dx_e - b_j*P_j  (row j).
Every module element gives a field a with r_j = 0 (tangent) for ALL j in the
set — boundary-safe relations with no down-shift into settled classes.

Certification is unconditional: exact flint divisibility P_j | G_j over Z.
The sampled fallback (pilot p1+p2 pattern) is degree-truncated and
LIGHTLY TESTED — prefer Singular.
"""
import os, re, json, time, subprocess, tempfile, itertools
import numpy as np

from . import core
from .relations import RelationSystem, RelationFamily

INT = np.int64
TERM = re.compile(
    r'([+-]?)(\d*)((?:\*?x\d+(?:\^\d+)?)*)\*?gen\((\d+)\)')


def _poly_str(pd, n):
    terms = []
    for mo, c in pd.items():
        mon = '*'.join(f'x{f+1}' + (f'^{mo[f]}' if mo[f] > 1 else '')
                       for f in range(n) if mo[f])
        terms.append(f'({c})' + (f'*{mon}' if mon else ''))
    return '+'.join(terms) or '0'


def singular_syz(spec, labs, timeout=600, workdir=None, keep=False):
    """Run Singular syz for tangency set labs. Returns list of a-dict tuples
    (one per module generator) or None on timeout/absence."""
    n = spec.n
    k = len(labs)
    wd = workdir or tempfile.mkdtemp(prefix='boxwalk_syz_')
    tag = '_'.join(labs)[:60]
    outtxt = os.path.join(wd, f'SYZ_{tag}.txt')
    lines = [f'ring r = 0, ({",".join(f"x{i+1}" for i in range(n))}), dp;']
    for e in range(n):
        lines.append(f'poly X{e+1} = x{e+1}-x{e+1}^2;')
    for j, lab in enumerate(labs):
        lines.append(f'poly P{j+1} = {_poly_str(spec.polys[lab], n)};')
    cols = []
    for e in range(n):
        ent = ','.join(f'X{e+1}*diff(P{j+1},x{e+1})' for j in range(k))
        cols.append(f'[{ent}]')
    for j in range(k):
        ent = ','.join(f'-P{j+1}' if i == j else '0' for i in range(k))
        cols.append(f'[{ent}]')
    lines += [f'module M = {",".join(cols)};',
              'module S = syz(M);',
              f'write(":w {outtxt}", string(S));', 'quit;']
    src = os.path.join(wd, f'syz_{tag}.sing')
    open(src, 'w').write('\n'.join(lines))
    try:
        subprocess.run(['Singular', '-q', src], capture_output=True,
                       text=True, timeout=timeout, check=False)
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return None
    if not os.path.exists(outtxt):
        return None
    gens = parse_syz(open(outtxt).read(), n)
    if not keep:
        for f in (src, outtxt):
            try:
                os.remove(f)
            except OSError:
                pass
    return gens


def parse_syz(txt, n):
    """Parse string(S) -> list of a-dict lists [{mono: coeff}]*n."""
    elems = txt.replace('\n', '').split(',')
    out = []
    for el in elems:
        comps = {}
        for sign, coef, monos, gi in TERM.findall(el):
            c = int(coef) if coef else 1
            if sign == '-':
                c = -c
            mono = [0] * n
            for mm in re.findall(r'x(\d+)(?:\^(\d+))?', monos):
                mono[int(mm[0]) - 1] += int(mm[1]) if mm[1] else 1
            gi = int(gi)
            d = comps.setdefault(gi, {})
            key = tuple(mono)
            d[key] = d.get(key, 0) + c
        a_dicts = [dict(comps.get(e + 1, {})) for e in range(n)]
        a_dicts = [{m: c for m, c in d.items() if c} for d in a_dicts]
        if any(a_dicts):
            out.append(a_dicts)
    return out


def certify(system, a_dicts, labs):
    """Exact flint certification: P_j | sum_e a_e X_e d_e P_j for all j in
    labs. Returns True/False."""
    n = system.n
    a = [system.ctx.from_dict({m: c for m, c in d.items() if c})
         for d in a_dicts]
    for lab in labs:
        G = system.zero()
        for e in range(n):
            G += a[e] * system.X[e] * system.P[lab].derivative(e)
        q, r = divmod(G, system.P[lab])
        if r != system.zero():
            return False
    return True


def certified_families(spec, system, settled, current, timeout=600,
                       dmax_fallback=4, log=print):
    """Generators tangent to `settled`, packaged as RelationFamily objects
    with classes = settled + [current]. Singular first; sampled-scan
    fallback on timeout. Returns (families, meta)."""
    meta = {'settled': list(settled), 'current': current, 'source': None,
            'n_raw': 0, 'n_certified': 0}
    gens = singular_syz(spec, settled, timeout=timeout)
    if gens is not None:
        meta['source'] = 'singular_syz'
    else:
        log(f"[syz] Singular timeout/missing for {settled}; "
            f"sampled-scan fallback (dmax={dmax_fallback})")
        gens = sampled_scan(spec, system, settled, dmax_fallback)
        meta['source'] = 'sampled_scan'
    meta['n_raw'] = len(gens)
    fams = []
    for a_dicts in gens:
        if certify(system, a_dicts, settled):
            fams.append(RelationFamily(system, a_dicts,
                                       list(settled) + [current]))
    meta['n_certified'] = len(fams)
    # mutation control: corrupt first certified generator -> must FAIL
    if fams:
        bad = [dict(d) for d in fams[0].a_dicts]
        k0 = next(iter(bad[0])) if bad[0] else (0,) * spec.n
        bad[0][k0] = bad[0].get(k0, 0) + 1
        meta['mutation_control_fails'] = not certify(system, bad, settled)
    return fams, meta


# ------------------------------------------------- sampled-scan fallback --
def _monomials(dmax, n):
    return [g for g in itertools.product(range(dmax + 1), repeat=n)
            if sum(g) <= dmax]


def _eval_poly(pd, X, p):
    """pd at points X (S,n) mod p -> (S,)."""
    S, n = X.shape
    acc = np.zeros(S, dtype=INT)
    for mo, c in pd.items():
        t = np.full(S, int(c) % p, dtype=INT)
        for f in range(n):
            for _ in range(mo[f]):
                t = (t * X[:, f]) % p
        acc = (acc + t) % p
    return acc


def _dpoly(pd, e):
    out = {}
    for mo, c in pd.items():
        if mo[e]:
            m2 = tuple(v - 1 if f == e else v for f, v in enumerate(mo))
            out[m2] = out.get(m2, 0) + c * mo[e]
    return {m: c for m, c in out.items() if c}


def _sample_on_divisor(pd, n, nn, p, rng):
    """nn points on {P=0} mod p: fix all coords random, solve the rotating
    coordinate as a univariate root (flint nmod_poly)."""
    import flint
    pts = []
    tries = 0
    while len(pts) < nn and tries < 200 * nn:
        tries += 1
        s = tries % n
        x = rng.integers(1, p, size=n).astype(INT)
        dmax_s = max(mo[s] for mo in pd)
        coeffs = [0] * (dmax_s + 1)
        for mo, c in pd.items():
            t = int(c) % p
            for f in range(n):
                if f != s:
                    t = t * pow(int(x[f]), mo[f], p) % p
            coeffs[mo[s]] = (coeffs[mo[s]] + t) % p
        if all(c == 0 for c in coeffs[1:]):
            continue
        f = flint.nmod_poly(coeffs, p)
        for root, _ in f.roots():
            y = x.copy()
            y[s] = int(root) % p
            pts.append(y)
            break
    return np.array(pts[:nn], dtype=INT)


def _build_rows(spec, labs, mons, nrows, p, rng):
    n = spec.n
    nmon = len(mons)
    per_j = max(40, nrows // max(1, len(labs)) + 1)
    mon_arr = np.array(mons, dtype=INT)
    rows = []
    for lab in labs:
        pd = spec.polys[lab]
        X = _sample_on_divisor(pd, n, per_j, p, rng)
        if len(X) == 0:
            continue
        S = X.shape[0]
        assert int(_eval_poly(pd, X, p).max()) == 0
        dmax = int(mon_arr.max())
        powt = np.ones((S, n, dmax + 1), dtype=INT)
        for k in range(1, dmax + 1):
            powt[:, :, k] = (powt[:, :, k - 1] * X) % p
        xg = np.ones((S, nmon), dtype=INT)
        for f in range(n):
            xg = (xg * powt[:, f, mon_arr[:, f]]) % p
        blk = np.empty((S, n * nmon), dtype=INT)
        for e in range(n):
            Xe = (X[:, e] * (1 - X[:, e])) % p
            sc = (Xe * _eval_poly(_dpoly(pd, e), X, p)) % p
            blk[:, e * nmon:(e + 1) * nmon] = (xg * sc[:, None]) % p
        rows.append(blk)
    return np.vstack(rows) if rows else np.zeros((0, n * nmon), dtype=INT)


def _nullspace_basis(A, p):
    """Canonical-RREF nullspace basis mod p (pilot p2 pattern)."""
    A = A % p
    m, nc = A.shape
    piv, row = [], 0
    for c in range(nc):
        if row == m:
            break
        col = A[row:, c]
        nz = np.nonzero(col)[0]
        if nz.size == 0:
            continue
        i = row + int(nz[0])
        if i != row:
            A[[row, i]] = A[[i, row]]
        inv = pow(int(A[row, c]), p - 2, p)
        A[row] = (A[row] * inv) % p
        f = A[:, c].copy()
        f[row] = 0
        nzr = np.nonzero(f)[0]
        if nzr.size:
            A[nzr] = (A[nzr] - f[nzr, None] * A[row]) % p
        piv.append(c)
        row += 1
    free = [c for c in range(nc) if c not in piv]
    basis = []
    for fcol in free:
        v = np.zeros(nc, dtype=INT)
        v[fcol] = 1
        for r_, pc in enumerate(piv):
            v[pc] = (-int(A[r_, fcol])) % p
        basis.append(v)
    return tuple(free), basis


def sampled_scan(spec, system, labs, dmax, primes=None, seed=20260707):
    """Degree-truncated sampled lift (pilot p1+p2): nullspace mod 2 primes,
    pair by free column, CRT + Wang ratrec, clear denominators.
    Returns list of a-dict lists (uncertified — caller certifies)."""
    from fractions import Fraction
    primes = primes or core.primes31(2)
    n = spec.n
    mons = _monomials(dmax, n)
    nmon = len(mons)
    ncols = n * nmon
    rng = np.random.default_rng(seed)
    bases = {}
    for p in primes[:2]:
        A = _build_rows(spec, labs, mons, int(1.6 * ncols) + 80, p, rng)
        if A.shape[0] < ncols:
            return []
        bases[p] = _nullspace_basis(A, p)
    pA, pB = primes[:2]
    if bases[pA][0] != bases[pB][0]:
        return []
    out = []
    M = pA * pB
    for vA, vB in zip(bases[pA][1], bases[pB][1]):
        fracs, ok = [], True
        for a1, a2 in zip(vA.tolist(), vB.tolist()):
            r, _ = core.crt_list([a1, a2], [pA, pB])
            fr = core.ratrec(r, M)
            if fr is None:
                ok = False
                break
            fracs.append(fr)
        if not ok:
            continue
        lcm = 1
        for fr in fracs:
            lcm = lcm * fr.denominator // np.gcd(lcm, fr.denominator)
        ints = [int(fr * lcm) for fr in fracs]
        g = int(np.gcd.reduce([abs(v) for v in ints if v] or [1]))
        ints = [v // g for v in ints]
        a_dicts = []
        for e in range(n):
            dd = {}
            for i, mo in enumerate(mons):
                c = ints[e * nmon + i]
                if c:
                    dd[mo] = c
            a_dicts.append(dd)
        if any(a_dicts):
            out.append(a_dicts)
    return out
