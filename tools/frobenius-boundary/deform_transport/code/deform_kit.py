#!/usr/bin/env python3
"""deform_kit.py — deformation-method matrix Frobenius transport: exact local
factors of operator motives (matrix wing of the deform_transport member).

Instrument:
 (i)  log-aware Frobenius basis at x=0 over Q: solutions as J-vector series
      (component k = coefficient of log^k / k!), delta-normalized at the
      resonance free positions, then echelonized by L=0 valuation; exact
      nilpotent N_struct with (d/dlog) Y = Y N_struct; Y(x,L) = Y0(x) e^{LN}.
 (ii) F(x) = Y0(x) E Y0(x^p)^{-1} with E in {E : N E = p E N} (exact basis).
      GATE GD: theta F = A_th F - p F A_th(x^p) holds iff the whole chain
      (series, jets, N_struct, commutant, assembly) is coherent.
 (iii) E pinned by REGULARITY (pole cancellation, the Laurent coefficients
      below grade 0 vanish) + p-integrality; linear digit solve; the residual
      unit scalar gamma from det/FE emergence; the last +-1 by a declared
      containment bit.
 (iv) evaluation at Teichmuller points (empirical tail certificate + the
      two-truncation stability law), division-free charpoly, Weil-box + FE
      integer isolation.

Scale algebra: stored = p^scale * true, residues mod p^modk; true values
certified mod p^{modk - scale}.  Divisions assert exactness (undersized
scale fails loudly, never silently).
Vendored engines consumed by sibling IMPORT (shas in receipts).
Failure-class vocabulary used in verdict strings: K1 = fixable bug (one
debug cycle), K2 = compute-priced wall, K3 = structural design wall.
Requires python-flint (series algebra) and gmpy2 (fast modular inverses);
both are guarded — a missing dependency refuses by name at first use.
"""
import json, os, sys, time, hashlib, datetime, subprocess, math
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lfact_kit as lk
r2, f2 = lk.r2, lk.f2
try:
    import flint
except ImportError:                                  # guarded: refuse by name at first use
    flint = None

MEMBER_DIR = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

def utc():
    return datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')

def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for ch in iter(lambda: f.read(1 << 20), b''):
            h.update(ch)
    return h.hexdigest()

KIT_PATH = os.path.abspath(__file__)

def snapshot_kit():
    """Content-addressed snapshot of this kit's source under
    work/kit_snapshots/, with a ledger work/KIT_SNAPSHOTS.json (both under the
    member's work dir, created on demand).  Returns the kit sha256."""
    sha = sha256_file(KIT_PATH)
    bdir = os.path.join(MEMBER_DIR, 'work', 'kit_snapshots')
    os.makedirs(bdir, exist_ok=True)
    dst = os.path.join(bdir, 'deform_kit_%s.py' % sha[:8])
    if not os.path.exists(dst):
        with open(KIT_PATH, 'rb') as fsrc, open(dst, 'wb') as fdst:
            fdst.write(fsrc.read())
        led = os.path.join(MEMBER_DIR, 'work', 'KIT_SNAPSHOTS.json')
        try:
            L = json.load(open(led))
        except Exception:
            L = {'producer_script': 'code/deform_kit.py (snapshot_kit)', 'producer_sha256': sha,
                 'stamp_utc': utc(), 'component': 'frobenius-boundary/deform_transport',
                 'snapshots': []}
        L['snapshots'].append({'sha256': sha, 'file': os.path.basename(dst), 'stamp_utc': utc()})
        L['producer_sha256'] = sha
        L['stamp_utc'] = utc()
        with open(led, 'w') as f:
            json.dump(L, f, indent=1)
    return sha

def producer(script_relpath):
    eng = {'deform_kit.py': snapshot_kit(),
           'lfact_kit.py': sha256_file(lk.__file__),
           'r2_kit.py': lk.R2_SHA, 'f2_lib.py': lk.F2_SHA}
    if os.environ.get('FROB_L6_OPERATOR') and os.path.exists(os.environ['FROB_L6_OPERATOR']):
        eng['l6_operator'] = lk.l6_sha()
    return {'component': 'frobenius-boundary/deform_transport',
            'script': script_relpath,
            'script_sha256': sha256_file(os.path.join(MEMBER_DIR, script_relpath)),
            'stamp_utc': utc(),
            'engine_shas': eng}

def lint(path):
    """Optional producer-lint hook: set FROB_PRODUCER_LINT to a lint script to
    check emitted receipts; skipped (with a printed note) when unset."""
    linter = os.environ.get('FROB_PRODUCER_LINT')
    if not linter or not os.path.exists(linter):
        line = 'producer_lint: SKIPPED (FROB_PRODUCER_LINT unset)'
        print(line)
        return line
    r = subprocess.run([sys.executable, linter, '--check', path],
                       capture_output=True, text=True)
    line = (r.stdout.strip() or r.stderr.strip())
    print(line)
    return line

def emit(rec, relpath):
    out = os.path.join(MEMBER_DIR, relpath)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, 'w') as f:
        json.dump(rec, f, indent=1, default=str)
    print('WROTE', out)
    return out, lint(out)

def valp(n, p):
    if n == 0:
        return None
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v

def centered(r, M):
    r %= M
    return r - M if r > M // 2 else r

# ===================================================================== OpData

class OpData:
    """Exact operator scaffold in the W_k(n) frame (r2_kit letter)."""
    def __init__(self, A, e, name):
        self.A, self.e, self.name = A, e, name
        A7, D, r_, s0, dmax, Wk = r2.op_scaffold(A, e)
        self.order, self.dmax, self.s0, self.D = r_, dmax, s0, D
        def ff_poly(i):
            c = [1]
            for t in range(i):
                a = e - t
                new = [0] * (len(c) + 1)
                for j, cj in enumerate(c):
                    new[j + 1] += cj
                    new[j] += a * cj
                c = new
            return c
        ffs = [ff_poly(i) for i in range(r_ + 1)]
        self.Wk_poly = []
        for k in range(dmax + 1):
            poly = [0] * (r_ + 2)
            for i in range(r_ + 1):
                j = i + s0 + k
                if 0 <= j <= D and A7[i][j]:
                    for t, ft in enumerate(ffs[i]):
                        poly[t] += A7[i][j] * ft
            while len(poly) > 1 and poly[-1] == 0:
                poly.pop()
            self.Wk_poly.append(poly)
        self.Wk = Wk
        for k in (0, 1, min(dmax, 7)):
            for n in (0, 3, 11, 29):
                v = sum(c * n ** t for t, c in enumerate(self.Wk_poly[k]))
                assert v == Wk(k, n), (k, n, 'Wk poly mismatch')
        W0 = self.Wk_poly[0]
        self.resonances = {}
        for m in range(0, 200):
            if sum(c * m ** t for t, c in enumerate(W0)) == 0:
                mult, d1 = 1, [t * c for t, c in enumerate(W0)][1:]
                while sum(c * m ** t for t, c in enumerate(d1)) == 0:
                    mult += 1
                    d1 = [t * c for t, c in enumerate(d1)][1:]
                self.resonances[m] = mult

    def W_jets(self, k, m, J):
        """[U^j] W_k(m + U) for j = 0..J-1, exact ints."""
        poly = self.Wk_poly[k]
        out = []
        for j in range(J):
            s = 0
            for t in range(j, len(poly)):
                b = 1
                for u in range(j):
                    b = b * (t - u) // (u + 1)
                s += poly[t] * b * (m ** (t - j))
            out.append(s)
        return out

    def theta_form(self):
        """c_i(x) coefficient lists with sum_i c_i(x) theta^i = the operator
        L = sum_k x^k W_k(theta) (the r2_kit source-index convention):
        c_i(x) = sum_k x^k [n^i] W_k(n)."""
        ci = [[0] * (self.dmax + 1) for _ in range(self.order + 1)]
        for k in range(self.dmax + 1):
            poly = self.Wk_poly[k]
            for i in range(min(len(poly), self.order + 1)):
                ci[i][k] = poly[i]
        return ci

# ====================================================== seed basis over Q

def toeplitz_solve_Q(a, b, freevals=None):
    """Solve sum_j a_j x_{i+j} = b_i over Q. v = mult of leading zeros of a.
    x_0..x_{v-1} free (freevals), x_v..x_{J-1} determined; rows J-1..J-v are
    consistency (must have b=0). Returns (x, fails)."""
    J = len(b)
    v = 0
    while v < len(a) and a[v] == 0:
        v += 1
    v = min(v, J)
    x = [Fraction(0)] * J
    fails = []
    for i in range(v):
        x[i] = Fraction(freevals[i]) if freevals and i < len(freevals) else Fraction(0)
    for i in range(J - 1, J - 1 - v, -1):
        s = b[i]
        for j in range(v, len(a)):
            if i + j < J and a[j] != 0:
                s -= a[j] * x[i + j]
        if s != 0:
            fails.append(i)
    for i in range(J - 1 - v, -1, -1):
        s = b[i]
        for j in range(v + 1, len(a)):
            if i + j < J and a[j] != 0:
                s -= a[j] * x[i + j]
        x[i + v] = s / a[v]
    return x, fails

def seed_basis(op, J, seedN):
    """Log-aware Frobenius basis over Q to m <= seedN (delta-normalized)."""
    basis, freepos = [], []
    for m in range(seedN + 1):
        a = op.W_jets(0, m, J)
        v = 0
        while v < len(a) and a[v] == 0:
            v += 1
        rhss = []
        for col in basis:
            rhs = [Fraction(0)] * J
            for k in range(1, min(m, op.dmax) + 1):
                cm = col.get(m - k)
                if cm is None:
                    continue
                wj = op.W_jets(k, m - k, J)
                for i in range(J):
                    s = Fraction(0)
                    for j in range(J - i):
                        if wj[j] and cm[i + j]:
                            s += wj[j] * cm[i + j]
                    rhs[i] += s
            rhss.append([-t for t in rhs])
        for col, rhs in zip(basis, rhss):
            x, fails = toeplitz_solve_Q(a, rhs)
            assert not fails, (m, 'log-depth consistency failure (J too small?)', fails)
            if any(x):
                col[m] = x
        if v > 0 and m in op.resonances:
            for comp in range(min(v, J)):
                fv = [Fraction(0)] * J
                fv[comp] = Fraction(1)
                x, fails = toeplitz_solve_Q(a, [Fraction(0)] * J, freevals=fv)
                assert not fails
                basis.append({m: x})
                freepos.append((m, comp))
    assert len(basis) == op.order, (len(basis), op.order, 'basis count != order')
    maxdepth = max((i for col in basis for vv in col.values()
                    for i in range(J) if vv[i] != 0), default=0)
    for i, (m, k) in enumerate(freepos):     # delta normalization check
        for jb, col in enumerate(basis):
            got = col.get(m, [Fraction(0)] * J)[k]
            assert got == (1 if jb == i else 0), 'delta normalization broken'
    N = [[Fraction(0)] * op.order for _ in range(op.order)]
    for j, col in enumerate(basis):
        for i, (m, k) in enumerate(freepos):
            vv = col.get(m)
            N[i][j] = vv[k + 1] if (vv is not None and k + 1 < J) else Fraction(0)
    return {'basis': basis, 'freepos': freepos, 'maxdepth': maxdepth,
            'N_struct': N, 'J': J, 'seedN': seedN}

def matmulQ(X, Y):
    n = len(X)
    return [[sum(X[i][t] * Y[t][j] for t in range(len(Y))) for j in range(len(Y[0]))]
            for i in range(n)]

def matinvQ(M):
    n = len(M)
    A_ = [ [Fraction(x) for x in row] + [Fraction(1 if i == j else 0) for j in range(n)]
           for i, row in enumerate(M)]
    for c in range(n):
        piv = next(r_ for r_ in range(c, n) if A_[r_][c] != 0)
        A_[c], A_[piv] = A_[piv], A_[c]
        inv = 1 / A_[c][c]
        A_[c] = [x * inv for x in A_[c]]
        for r_ in range(n):
            if r_ != c and A_[r_][c] != 0:
                f_ = A_[r_][c]
                A_[r_] = [x - f_ * y for x, y in zip(A_[r_], A_[c])]
    return [row[n:] for row in A_]

def jordan_of_nilpotent(N):
    n = len(N)
    def rankQ(M):
        M = [row[:] for row in M]
        rr = 0
        for c in range(n):
            piv = next((r_ for r_ in range(rr, n) if M[r_][c] != 0), None)
            if piv is None:
                continue
            M[rr], M[piv] = M[piv], M[rr]
            inv = 1 / M[rr][c]
            M[rr] = [x * inv for x in M[rr]]
            for r_ in range(n):
                if r_ != rr and M[r_][c] != 0:
                    f_ = M[r_][c]
                    M[r_] = [x - f_ * y for x, y in zip(M[r_], M[rr])]
            rr += 1
        return rr
    ranks = []
    P = [row[:] for row in N]
    for _ in range(n):
        ranks.append(rankQ(P))
        P = matmulQ(P, N)
    nilp = ranks[-1] == 0
    rseq = [n] + ranks
    ge = [rseq[k - 1] - rseq[k] for k in range(1, n + 1)]
    blocks = []
    for k in range(1, n + 1):
        cnt = ge[k - 1] - (ge[k] if k < n else 0)
        blocks += [k] * max(0, cnt)
    return sorted(blocks, reverse=True), nilp

def echelonize_basis(sb, op):
    """Grade columns by the JET-column valuation dgrade = min m with a nonzero
    coefficient VECTOR (the theta-jets of a log column reach x^dgrade through
    the log parts).  Columns sharing a grade with PARALLEL lead vectors are
    subtracted to raise the grade; independent leads at a shared grade are
    fine (Lambda handles them).  Updates basis + N_struct; adds deff."""
    J = sb['J']
    n = op.order
    def dgrade(col):
        for m in sorted(col):
            if any(col[m]):
                return m
        return None
    basis = sb['basis']
    Tmat = [[Fraction(1 if i == j else 0) for j in range(n)] for i in range(n)]
    guard = 0
    while True:
        guard += 1
        assert guard < 120, 'echelonization loop guard'
        vals = [dgrade(c) for c in basis]
        assert all(v is not None for v in vals), 'zero column: basis broken'
        assert max(vals) <= sb['seedN'] - 2, ('echelon needs longer seed range', vals)
        clash = None
        for j in range(n):
            for i in range(j):
                if vals[i] != vals[j]:
                    continue
                a_ = basis[i][vals[i]]
                b_ = basis[j][vals[j]]
                # parallel test
                ratio = None
                par = True
                for x, y in zip(a_, b_):
                    if x == 0 and y == 0:
                        continue
                    if x == 0 or y == 0:
                        par = False
                        break
                    r_ = y / x
                    if ratio is None:
                        ratio = r_
                    elif r_ != ratio:
                        par = False
                        break
                if par and ratio is not None:
                    clash = (i, j, ratio)
                    break
            if clash:
                break
        if clash is None:
            break
        i, j, f_ = clash
        newcol = {}
        for m in set(basis[j]) | set(basis[i]):
            a_ = basis[j].get(m, [Fraction(0)] * J)
            b_ = basis[i].get(m, [Fraction(0)] * J)
            vec = [x - f_ * y for x, y in zip(a_, b_)]
            if any(vec):
                newcol[m] = vec
        basis[j] = newcol
        for t in range(n):
            Tmat[t][j] -= f_ * Tmat[t][i]
    sb['deff'] = [dgrade(c) for c in basis]
    sb['N_struct'] = matmulQ(matmulQ(matinvQ(Tmat), sb['N_struct']), Tmat)
    sb['echelon_T'] = [[str(x) for x in row] for row in Tmat]
    jordan, nilp = jordan_of_nilpotent(sb['N_struct'])
    sb['jordan'], sb['nilpotent'] = jordan, nilp
    return sb

def basis_maxdenval(sb, p):
    d = 0
    for col in sb['basis']:
        for vv in col.values():
            for x in vv:
                dv = valp(x.denominator, p)
                if dv:
                    d = max(d, dv)
    return d

# ====================================================== ramped-scale march

def march_ramped(op, sb, p, Ntot, keep_true, sigma, Jm=None, vpad=None,
                 progress=None):
    """RAMPED-SCALE continuation: coefficient m stored as p^{v(m)} * c_m with
    v(m) = ledger-cumsum(m) + vpad (v covers the measured true denominator
    growth = the ledger law).  NO p-divisions occur (pivot valuations are
    absorbed by the ramp; only unit inverses), so residues stay SMALL:
    stored mod p^{Sr}, Sr = keep_true + vpad + 16.  Archive: sols[j][m] =
    Jm-vector of ints at scale v(m), plus the v-array.  Downstream converts
    to the uniform sigma scale (exact, multiplications only)."""
    t0 = time.time()
    J = sb['J']
    Jm = Jm or (sb['maxdepth'] + 1)
    basis, seedN = sb['basis'], sb['seedN']
    m_start = seedN + 1
    Ltot, sing = r2.loss_ledger(op.A, op.e, p, Ntot, m_start)
    vpad = vpad if vpad is not None else (J + 6)
    Sr = keep_true + vpad + 16
    M = p ** Sr
    Mz = r2.mpz(M)
    singv = dict(sing)
    v = [0] * (Ntot + 1)
    run = 0
    for m in range(Ntot + 1):
        run += singv.get(m, 0)
        v[m] = run + vpad
    rec = {'p': p, 'Ntot': Ntot, 'Sr': Sr, 'vpad': vpad, 'ramp_final': v[-1],
           'Jm': Jm, 'L_ledger_scalar': Ltot, 'n_singular': len(sing),
           'mode': 'ramped'}
    nsol = len(basis)
    sols = [[None] * (Ntot + 1) for _ in range(nsol)]
    window = {}
    for m in range(seedN + 1):
        row = []
        for jcol in range(nsol):
            vec = basis[jcol].get(m)
            out = []
            for i in range(Jm):
                fr = vec[i] if vec else Fraction(0)
                num, den = fr.numerator, fr.denominator
                dv = valp(den, p) or 0
                assert dv <= v[m], ('seed denominator exceeds ramp', m, jcol, i, dv)
                out.append(int(r2.mpz(num) * (r2.mpz(p) ** (v[m] - dv)) % M
                               * r2.INV(den // p ** dv, M) % M))
            row.append(out)
        window[m] = row
        for jcol in range(nsol):
            sols[jcol][m] = row[jcol]
    pw_cache = {}
    for m in range(m_start, Ntot + 1):
        kmax = min(m, op.dmax)
        a = op.W_jets(0, m, Jm)
        pivv = valp(a[0], p) or 0
        assert v[m] - v[m - 1] == pivv if m > 0 else True
        rhs = [[r2.mpz(0)] * Jm for _ in range(nsol)]
        for k in range(1, kmax + 1):
            row = window.get(m - k)
            if row is None:
                continue
            wj = op.W_jets(k, m - k, Jm)
            if not any(wj):
                continue
            dlt = v[m] - v[m - k]
            pd = pw_cache.get(dlt)
            if pd is None:
                pd = r2.mpz(p) ** dlt % Mz if dlt < Sr else 0
                pw_cache[dlt] = pd
            if pd == 0:
                continue
            wjm = [r2.mpz(w % M) * pd % Mz for w in wj]
            nz = [j for j in range(Jm) if wj[j]]
            for s_ in range(nsol):
                cm = row[s_]
                rv = rhs[s_]
                for i in range(Jm):
                    acc = rv[i]
                    for j in nz:
                        ij = i + j
                        if ij < Jm and cm[ij]:
                            acc += wjm[j] * cm[ij]
                    rv[i] = acc % Mz
        # Toeplitz solve with pivot a0 = p^pivv * unit: the p^pivv cancels the
        # ramp increment (v[m] = v[m-1] + pivv), so we divide by p^pivv AND
        # rescale to v[m]: net = divide by unit only.
        unit = a[0] // p ** pivv
        uinv = r2.INV(r2.mpz(unit % M), Mz)
        # for the off-diagonal jets a_j (j>=1): they multiply SAME-step values
        # (already at v[m]); the equation was scaled by p^{v[m]}: consistent.
        am = [r2.mpz(x % M) for x in a]
        newrow = []
        for s_ in range(nsol):
            x = [0] * Jm
            for i in range(Jm - 1, -1, -1):
                sacc = (-rhs[s_][i]) % Mz
                for j in range(1, Jm - i):
                    if am[j] and x[i + j]:
                        sacc = (sacc - am[j] * x[i + j]) % Mz
                # solve p^pivv*unit * x_i = sacc at scale v[m]: the stored
                # unknown x̃_i = p^{v[m]} x_i: p^pivv*unit*x_i = sacc/p^{v[m]}
                # => x̃_i = sacc * uinv / p^{pivv} * ... : rescaling algebra:
                # sacc is p^{v[m]}*(true rhs); true x_i = true_rhs/(p^pivv u):
                # x̃_i = p^{v[m]} x_i = sacc/(p^pivv u) — need divisibility:
                if pivv:
                    assert sacc % (p ** pivv) == 0, ('ramp shortfall', m, s_, i)
                    sacc //= p ** pivv
                    sacc = sacc * (r2.mpz(p) ** pivv) % Mz if False else sacc
                    # NOTE: dividing by p^pivv then... true x̃ = p^{v[m]}x =
                    # p^{v[m]} * truerhs/(p^pivv u) and sacc = p^{v[m]}truerhs
                    # => x̃ = (sacc/p^pivv)*u^{-1} * p^{... 0}: wait:
                    # x̃ = sacc * u^{-1} / p^{pivv}: yes as done.
                x[i] = int(sacc * uinv % Mz)
            newrow.append(x)
        window[m] = newrow
        for s_ in range(nsol):
            sols[s_][m] = newrow[s_]
        window.pop(m - op.dmax - 1, None)
        if progress and m % progress == 0:
            print('  ramped m=%d/%d t=%.1fs' % (m, Ntot, time.time() - t0), flush=True)
    rec['t_wall_s'] = round(time.time() - t0, 2)
    return {'sols': sols, 'rec': rec, 'v': v, 'Sr': Sr, 'sigma': sigma,
            'Jm': Jm, 'p': p, 'Ntot': Ntot, 'mode': 'ramped',
            'keep_true': keep_true}

# ====================================================== mod-p^S march

def march(op, sb, p, Ntot, keep, sigma, Jm=None, progress=None):
    """Continue the echelonized basis mod p^{S} (S = keep + Jm*L + 8) with
    fixed-point scale p^sigma; archive coefficients reduced mod p^keep.
    keep should be >= (true abs target) + sigma."""
    t0 = time.time()
    J = sb['J']
    Jm = Jm or (sb['maxdepth'] + 1)
    basis, seedN = sb['basis'], sb['seedN']
    m_start = seedN + 1
    Ltot, sing = r2.loss_ledger(op.A, op.e, p, Ntot, m_start)
    # cascade allocation: worst case is Jm*Ltot but the realized cascade is
    # ~Ltot (derivative-valuation compensation); allocate 2*Ltot + 96 and let
    # the exact-divisibility asserts catch any shortfall loudly.
    S = keep + 2 * Ltot + 96
    M = int(r2.mpz(p) ** S)
    Mkeep = p ** keep
    rec = {'p': p, 'Ntot': Ntot, 'S_alloc': S, 'keep': keep, 'sigma': sigma,
           'Jm': Jm, 'L_ledger_scalar': Ltot, 'n_singular': len(sing)}
    nsol = len(basis)
    sols = [[None] * (Ntot + 1) for _ in range(nsol)]
    window = {}
    for m in range(seedN + 1):
        row = []
        for jcol in range(nsol):
            vec = basis[jcol].get(m)
            out = []
            for i in range(Jm):
                fr = vec[i] if vec else Fraction(0)
                num, den = fr.numerator, fr.denominator
                dv = 0
                while den % p == 0:
                    den //= p
                    dv += 1
                assert dv <= sigma, ('seed denominator exceeds sigma', m, jcol, i, dv)
                out.append(int(r2.mpz(num) * (r2.mpz(p) ** (sigma - dv)) % M * r2.INV(den, M) % M))
            row.append(out)
        window[m] = row
        for jcol in range(nsol):
            sols[jcol][m] = [v % Mkeep for v in row[jcol]]
    loss = 0
    shrink_at = 0
    suffix_after = {}
    run = 0
    for i in range(len(sing) - 1, -1, -1):
        suffix_after[sing[i][0]] = run
        run += sing[i][1] * Jm
    Mz = r2.mpz(M)
    for m in range(m_start, Ntot + 1):
        kmax = min(m, op.dmax)
        a = op.W_jets(0, m, Jm)
        assert a[0] != 0, (m, 'char-0 resonance above seed range')
        v = valp(a[0], p) or 0
        rhs = [[r2.mpz(0)] * Jm for _ in range(nsol)]
        for k in range(1, kmax + 1):
            row = window.get(m - k)
            if row is None:
                continue
            wj = op.W_jets(k, m - k, Jm)
            if not any(wj):
                continue
            wjm = [r2.mpz(w % M) for w in wj]
            nz = [j for j in range(Jm) if wj[j]]
            for s_ in range(nsol):
                cm = row[s_]
                rv = rhs[s_]
                for i in range(Jm):
                    acc = rv[i]
                    for j in nz:
                        ij = i + j
                        if ij < Jm and cm[ij]:
                            acc += wjm[j] * cm[ij]
                    rv[i] = acc % Mz
        pv = p ** v
        unit = a[0] // pv
        uinv = r2.INV(r2.mpz(unit % M), Mz)
        am = [r2.mpz(x % M) for x in a]
        newrow = []
        for s_ in range(nsol):
            x = [0] * Jm
            for i in range(Jm - 1, -1, -1):
                sacc = (-rhs[s_][i]) % Mz
                for j in range(1, Jm - i):
                    if am[j] and x[i + j]:
                        sacc = (sacc - am[j] * x[i + j]) % Mz
                if v:
                    assert sacc % pv == 0, ('sigma exceeded (denominator beyond allocation)', m, s_, i)
                    sacc //= pv
                x[i] = sacc * uinv % Mz
            newrow.append(x)
        loss += v * Jm
        window[m] = newrow
        for s_ in range(nsol):
            sols[s_][m] = [int(t_ % Mkeep) for t_ in newrow[s_]]
        window.pop(m - op.dmax - 1, None)
        if v and loss - shrink_at >= 64:
            Srem = keep + suffix_after.get(m, 0) + 8
            Mnew = int(r2.mpz(p) ** Srem)
            if Mnew < M:
                M, Mz = Mnew, r2.mpz(Mnew)
                for mm in list(window):
                    window[mm] = [[int(t_ % M) for t_ in vec] for vec in window[mm]]
            shrink_at = loss
        if progress and m % progress == 0:
            print('  march m=%d/%d loss=%d t=%.1fs' % (m, Ntot, loss, time.time() - t0), flush=True)
    rec['loss_realized_scaled'] = loss
    rec['ledger_expected_scaled_max'] = Ltot * Jm
    rec['t_wall_s'] = round(time.time() - t0, 2)
    return {'sols': sols, 'rec': rec, 'keep': keep, 'sigma': sigma, 'Jm': Jm, 'p': p,
            'Ntot': Ntot}

# ====================================================== Y0 + series algebra

def y0_columns(op, mr):
    """Y0[i][j] coefficient lists (ints scaled p^sigma mod p^keep):
    [x^n] Y0[i][j] = sum_t C(i,t) n^{i-t} sols[j][n][t]."""
    sols, keep, Jm, p = mr['sols'], mr['keep'], mr['Jm'], mr['p']
    Mk = p ** keep
    n_ = op.order
    Ntot = mr['Ntot']
    binom = [[math.comb(i, t) for t in range(n_ + 1)] for i in range(n_ + 1)]
    Y = [[None] * n_ for _ in range(n_)]
    for j in range(n_):
        col = sols[j]
        for i in range(n_):
            out = [0] * (Ntot + 1)
            for n in range(Ntot + 1):
                vec = col[n]
                if vec is None:
                    continue
                acc = 0
                for t in range(min(i, Jm - 1), -1, -1):
                    if vec[t]:
                        acc += binom[i][t] * pow(n, i - t) * vec[t]
                if acc:
                    out[n] = acc % Mk
            Y[i][j] = out
    return Y

class Ser:
    """Series toolbox over Z/p^modk via flint fmpz_mod_poly."""
    def __init__(self, p, modk):
        if flint is None:
            raise RuntimeError('deform_kit.Ser: python-flint is required for the '
                               'series algebra (pip install python-flint)')
        self.p, self.modk = p, modk
        self.mod = p ** modk
        self.ctx = flint.fmpz_mod_poly_ctx(flint.fmpz(p) ** modk)
    def mk(self, coeffs):
        return self.ctx(coeffs)
    def mul(self, a, b, T):
        c = a * b
        return self.ctx([c[i] for i in range(min(c.degree() + 1, T))]) if c.degree() >= 0 else c
    def coeffs(self, a, T):
        return [int(a[i]) for i in range(T)]
    def subs_xp(self, coeffs, T):
        out = [0] * T
        for n, c in enumerate(coeffs):
            if n * self.p >= T:
                break
            out[n * self.p] = c
        return out

def mat_mul_ser(S, A, B, T):
    n_ = len(A)
    m_ = len(B[0])
    K = len(B)
    out = [[None] * m_ for _ in range(n_)]
    for i in range(n_):
        for j in range(m_):
            acc = S.mk(0)
            for k in range(K):
                acc = acc + A[i][k] * B[k][j]
            d = acc.degree()
            out[i][j] = S.mk([acc[t] for t in range(min(d + 1, T))]) if d >= 0 else acc
    return out

def mat_newton_inv_ser(S, M, lead_inv_int, T):
    """Inverse of series matrix M (flint polys) whose x^0 matrix has inverse
    lead_inv_int (int matrix mod p^modk). Newton doubling to length T."""
    n_ = len(M)
    X = [[S.mk(lead_inv_int[i][j]) for j in range(n_)] for i in range(n_)]
    t = 1
    while t < T:
        t = min(2 * t, T)
        MX = mat_mul_ser(S, M, X, t)
        R = [[(S.mk(2 if i == j else 0) - MX[i][j]) for j in range(n_)] for i in range(n_)]
        X = mat_mul_ser(S, X, R, t)
    return X

# ====================================================== commutant + pinning

def commutant_basis(N, p, order):
    """Exact basis over Q of {E : N E = p E N} as integer matrices."""
    n_ = order
    rows = []
    for i in range(n_):
        for j in range(n_):
            row = [Fraction(0)] * (n_ * n_)
            for t in range(n_):
                row[t * n_ + j] += N[i][t]
                row[i * n_ + t] -= p * N[t][j]
            rows.append(row)
    A_ = [r_[:] for r_ in rows]
    nv = n_ * n_
    piv_of_col = {}
    rr = 0
    for c in range(nv):
        piv = next((r_ for r_ in range(rr, len(A_)) if A_[r_][c] != 0), None)
        if piv is None:
            continue
        A_[rr], A_[piv] = A_[piv], A_[rr]
        inv = 1 / A_[rr][c]
        A_[rr] = [x * inv for x in A_[rr]]
        for r_ in range(len(A_)):
            if r_ != rr and A_[r_][c] != 0:
                f_ = A_[r_][c]
                A_[r_] = [x - f_ * y for x, y in zip(A_[r_], A_[rr])]
        piv_of_col[c] = rr
        rr += 1
    free = [c for c in range(nv) if c not in piv_of_col]
    out = []
    for fc in free:
        vec = [Fraction(0)] * nv
        vec[fc] = Fraction(1)
        for c, r_ in piv_of_col.items():
            vec[c] = -A_[r_][fc]
        den = 1
        for x in vec:
            den = den * x.denominator // math.gcd(den, x.denominator)
        out.append([[int(vec[i * n_ + j] * den) for j in range(n_)] for i in range(n_)])
    return out

def padic_kernel(rows, p, K):
    """Solution structure of rows . e == 0 (mod p^K), e in (Z/p^K)^nv, by
    valuation-aware elimination.  Returns dict:
      pivots: list of (col, val) — column col is determined mod p^{K-val}
      free_cols: columns with full (unit) freedom
      gens: one generator per free col (kernel vectors mod p^K), built by
            back-substitution; a free col whose back-substitution requires a
            non-exact division is reported in 'obstructed' with the reachable
            p-power scaling instead.
    """
    mod = p ** K
    nv = len(rows[0]) if rows else 0
    A = []
    for r_ in rows:
        rr = [x % mod for x in r_]
        if any(rr):
            A.append(rr)
    piv = {}                      # col -> normalized row (pivot entry = p^v * 1)
    pivval = {}
    for r_ in A:
        r_ = r_[:]
        for _ in range(nv + 2):
            best = None
            for c in range(nv):
                if r_[c]:
                    v = valp(r_[c], p)
                    if best is None or v < best[1]:
                        best = (c, v)
            if best is None:
                break
            c, v = best
            if v >= K:
                break
            if c not in piv:
                unit = r_[c] // p ** v
                uinv = pow(unit % mod, -1, mod)
                r_ = [x * uinv % mod for x in r_]
                piv[c] = r_
                pivval[c] = v
                break
            pv = pivval[c]
            if v >= pv:
                q = (r_[c] // p ** pv) % mod
                pr = piv[c]
                r_ = [(x - q * y) % mod for x, y in zip(r_, pr)]
            else:
                unit = r_[c] // p ** v
                uinv = pow(unit % mod, -1, mod)
                rnew = [x * uinv % mod for x in r_]
                old = piv[c]
                piv[c] = rnew
                pivval[c] = v
                q = (old[c] // p ** v) % mod
                r_ = [(x - q * y) % mod for x, y in zip(old, rnew)]
    # full inter-reduction of pivot rows (bounded passes)
    for _ in range(len(piv) + 2):
        changed = False
        for c in list(piv):
            r_ = piv[c]
            for c2 in list(piv):
                if c2 == c:
                    continue
                v2 = pivval[c2]
                if r_[c2] and (valp(r_[c2], p) or 0) >= v2:
                    q = (r_[c2] // p ** v2) % mod
                    if q:
                        r_ = [(x - q * y) % mod for x, y in zip(r_, piv[c2])]
                        changed = True
            piv[c] = r_
        if not changed:
            break
    free = [c for c in range(nv) if c not in piv]
    gens, obstructed = [], []
    for fc in free:
        vec = [0] * nv
        vec[fc] = 1
        ok = True
        # back-substitute in decreasing pivot dependence (iterate to fixpoint)
        for _ in range(len(piv) + 2):
            done = True
            for c in sorted(piv):
                r_ = piv[c]
                v = pivval[c]
                s = sum(r_[c2] * vec[c2] for c2 in range(nv) if c2 != c) % mod
                # need p^v * e_c == -s  => e_c = -s / p^v (exact division required)
                if s % (p ** v) != 0:
                    ok = False
                    break
                want = (-(s // p ** v)) % (p ** (K - v))
                if vec[c] % (p ** (K - v)) != want:
                    vec[c] = want
                    done = False
            if not ok or done:
                break
        if ok and all(sum(r * e for r, e in zip(row, vec)) % mod == 0 for row in A[:40]):
            gens.append(vec)
        else:
            obstructed.append(fc)
    return {'pivots': sorted((c, pivval[c]) for c in piv), 'free_cols': free,
            'gens': gens, 'obstructed': obstructed, 'n_rows': len(A), 'nv': nv}

def nullspace_modpk(rows, p, K):
    """Kernel generators of rows.e == 0 mod p^K (unit-pivot Gaussian);
    returns (gens, n_unit_pivots, n_free)."""
    mod = p ** K
    A_ = [[x % mod for x in r_] for r_ in rows if any(x % mod for x in r_)]
    if not A_:
        nv = len(rows[0])
        return [[1 if i == j else 0 for j in range(nv)] for i in range(nv)], 0, len(rows[0])
    nv = len(A_[0])
    pivcols = []
    rr = 0
    for c in range(nv):
        piv = next((r_ for r_ in range(rr, len(A_)) if A_[r_][c] % p != 0), None)
        if piv is None:
            continue
        A_[rr], A_[piv] = A_[piv], A_[rr]
        inv = pow(A_[rr][c], -1, mod)
        A_[rr] = [x * inv % mod for x in A_[rr]]
        for r_ in range(len(A_)):
            if r_ != rr and A_[r_][c]:
                f_ = A_[r_][c]
                A_[r_] = [(x - f_ * y) % mod for x, y in zip(A_[r_], A_[rr])]
        pivcols.append(c)
        rr += 1
    free = [c for c in range(nv) if c not in pivcols]
    gens = []
    for fc in free:
        vec = [0] * nv
        vec[fc] = 1
        for idx, c in enumerate(pivcols):
            vec[c] = (-A_[idx][fc]) % mod
        gens.append(vec)
    return gens, len(pivcols), len(free)

# ====================================================== charpoly + isolation

def berkowitz_charpoly(Mmat, mod):
    """Division-free charpoly det(xI - M) mod `mod`; ascending coeff list,
    length n+1, leading coefficient 1."""
    n_ = len(Mmat)
    C = [1, (-Mmat[0][0]) % mod]
    for k in range(1, n_):
        R = [Mmat[k][j] % mod for j in range(k)]
        Ccol = [Mmat[i][k] % mod for i in range(k)]
        A_ = [row[:k] for row in Mmat[:k]]
        s = [1, (-Mmat[k][k]) % mod]
        v = Ccol[:]
        for t in range(1, k + 1):
            dot = sum(R[i] * v[i] for i in range(k)) % mod
            s.append((-dot) % mod)
            if t < k:
                v = [sum(A_[i][j] * v[j] for j in range(k)) % mod for i in range(k)]
        newC = [0] * (len(C) + len(s) - 1)
        for i, si in enumerate(s):
            if si:
                for j, cj in enumerate(C):
                    newC[i + j] = (newC[i + j] + si * cj) % mod
        C = newC[:k + 2]
    return [c % mod for c in reversed(C)]

def rp_coeffs(cp_asc, p, K):
    """det(1 - T F) coefficients from charpoly det(xI - F) (ascending, monic):
    det(xI-F) = sum_j a_j x^j = prod(x - l_i); det(1-TF) = prod(1 - T l_i)
    = T^n sum_j a_j T^{-j} = sum_k a_{n-k} T^k  =>  c_k = a_{n-k} mod p^K."""
    mod = p ** K
    n_ = len(cp_asc) - 1
    return [cp_asc[n_ - k] % mod for k in range(n_ + 1)]

def isqrt_pow(p, wt, k):
    """floor(p^{wt*k/2}) exactly."""
    if (wt * k) % 2 == 0:
        return p ** (wt * k // 2)
    return math.isqrt(p ** (wt * k))

def isolate_factor(c, p, K, wt=4, n=6):
    """Integer isolation from c_k = [T^k] det(1-TF) mod
    p^K; weight wt (|roots| = p^{wt/2}, FE pairs to p^wt); n even.
    Isolates c_1..c_{n/2}; FE forces the top half.  Emergence residuals of the
    top half at depth K are the GF gate."""
    mod = p ** K
    half = n // 2
    out = {'K': K, 'n': n, 'wt': wt, 'c_mod': [int(x % mod) for x in c]}
    ints = {}
    ok = True
    for k in range(1, half + 1):
        bound = math.comb(n, k) * isqrt_pow(p, wt, k)
        Kk = 1
        while p ** Kk <= 2 * bound:
            Kk += 1
        if Kk > K:
            ints[k] = None
            ok = False
            continue
        modk = p ** min(K, Kk + 1)
        val = centered(c[k] % modk, modk)
        inbox = abs(val) <= bound
        ints[k] = val if inbox else None
        ok = ok and inbox
    out['c_low'] = {str(k): ints[k] for k in ints}
    out['weil_box_ok'] = ok
    if ok:
        R = [1] + [ints[k] for k in range(1, half + 1)]
        for j in range(half + 1, n + 1):
            R.append(p ** (wt * (2 * j - n) // 2) * (ints[n - j] if n - j >= 1 else 1))
        fe = {}
        for j in range(half + 1, n + 1):
            diff = (c[j] - R[j]) % mod
            fe['c%d' % j] = None if diff == 0 else (valp(diff, p) or 0)
        out['FE_residual_vals_at_K'] = fe
        out['FE_emerges_full_at_K'] = all(v is None for v in fe.values())
        out['Rp_integer'] = R
        out['weil_circle_ok'] = bool(lk.weil_circle_ok(R, p, wt=wt))
        out['ordinary_c1_unit'] = (ints[1] % p != 0)
    return out

def elliptic_containment(R, p, ap):
    """C2: exact divisibility by the declared elliptic factor shapes over Z."""
    ap2 = ap * ap - 2 * p
    def divides(quad):
        q = [Fraction(x) for x in R]
        quot = [Fraction(0)] * (len(R) - len(quad) + 1)
        d = quad
        for i in range(len(quot) - 1, -1, -1):
            c_ = q[i + len(quad) - 1] / d[-1]
            quot[i] = c_
            for j in range(len(quad)):
                q[i + j] -= c_ * d[j]
        return all(x == 0 for x in q[:len(quad) - 1])
    return {'ap': ap, 'ap2': ap2,
            'S1_div_1_papt2_p4': divides([1, -p * ap2, p ** 4]),
            'S2_div_1_pap_p3': divides([1, -p * ap, p ** 3]),
            'S2b_div_1_p2ap_p5': divides([1, -p * p * ap, p ** 5])}

def teich(x0, p, K):
    """Teichmuller lift of fiber x0 (int or (num,den)) mod p^K."""
    mod = p ** K
    if isinstance(x0, tuple):
        n, d = x0
        return pow(n, p ** K, mod) * pow(pow(d, p ** K, mod), -1, mod) % mod
    return pow(x0, p ** K, mod)

# ====================================================== operator loaders

def op_L6():
    return OpData(r2.load_L6(lk.l6_path()), -2, 'L6(e=-2)')

def op_legendre():
    return OpData(r2.HYP2F1, 0, 'HYP2F1')

def op_sym5(path=None):
    path = path or os.path.join(MEMBER_DIR, 'work', 'SYM5_OP.json')
    d = json.load(open(path))
    return OpData([[int(c) for c in row] for row in d['coeffs']], 0, 'Sym5Legendre')

# ====================================================== full transport driver

DEFAULT_TWISTS = {
    # (coeffs, name, mult): twist power = mult * s_twist; mult follows the
    # factor's multiplicity in the theta-leading coefficient (measured law:
    # Sym^5 c6 = -64(1-x)^5 needed 5x the Legendre rate).
    'HYP2F1': [([1, -1], '1-x', 1)],
    'Sym5Legendre': [([1, -1], '1-x', 5)],
    'L6(e=-2)': [([1, -1], '1-x', 3), ([1, 1], '1+x', 5),
                 ([1, 1, 1], 'Phi3', 1), ([1, 1, 1, 1, 1], 'Phi5', 1)],
}

def transport(op, p, K_true, fibers=None, NF=None, sigma=None, sigma2=48, extraK=4,
              progress=None, twist=None, s_twist=None, win_per_digit=3.0):
    """End-to-end matrix transport at prime p.
    K_true: target certified depth for F(xhat) charpoly digits.
    twist: list of (coeff-list, name) singular factors Delta_i(x); the pinning
    and evaluation ride Ptw = prod Delta_i(x^p)^{s_twist} (the measured
    denominator law).  Sizes: T ~ p*maxdeff +
    p*s*degDelta + p*win_per_digit*(dtarget) + slack, all recorded.
    Honest failures raise AssertionError (caller receipts them)."""
    t0 = time.time()
    n_ = op.order
    J = n_
    sb = seed_basis(op, J, max(op.resonances) + 10 if op.resonances else 10)
    sb = echelonize_basis(sb, op)
    maxdeff = max(sb['deff'])
    twist = twist if twist is not None else DEFAULT_TWISTS.get(op.name, [([1, -1], '1-x')])
    dtarget = K_true + 6
    s_twist = s_twist if s_twist is not None else dtarget + 2
    twist = [(t[0], t[1], t[2] if len(t) > 2 else 1) for t in twist]
    degD = sum((len(t[0]) - 1) * t[2] for t in twist)
    if NF is None:
        NF = int(p * (s_twist * degD + win_per_digit * dtarget + 8))
    Ntot = NF + p * maxdeff + 16
    Jm = sb['maxdepth'] + 1
    # sigma auto-law: the TRUE denominator growth of the log-basis follows the
    # ledger (measured on Sym^5: exactly the W0-root-multiplicity sum per
    # p-window); allocate sigma = L(Ntot) + J + 12.
    Ltot_pre, _sing = r2.loss_ledger(op.A, op.e, p, Ntot, sb['seedN'] + 1)
    if sigma is None:
        sigma = Ltot_pre + J + 12
    dbar = basis_maxdenval(sb, p)
    assert dbar <= sigma - 8, ('sigma too small for measured seed denominators', dbar)
    keep_true = K_true + extraK + 14        # margin for det-val extraction + FE
    keep = keep_true + sigma + sigma2
    # Newton budget: realized divide-back loss measured ~50 digits across the
    # controls (worst case would be rounds*(sigma+sigma2)); allocate 192 and
    # let the post-hoc certificate assert (inv_residual_certified_depth >=
    # keep) catch any shortfall loudly.
    keepN = keep + 192
    res = {'p': p, 'op': op.name, 'K_true': K_true, 'keep_true': keep_true,
           'sigma': sigma, 'sigma2': sigma2, 'NF': NF, 'Ntot': Ntot,
           'keepN_budget': keepN,
           'deff': sb['deff'], 'jordan': sb['jordan'], 'nilpotent': sb['nilpotent'],
           'maxdepth_log': sb['maxdepth'], 'dbar_seed': dbar}
    # march at the budgeted modulus
    mr = march(op, sb, p, Ntot, keepN, sigma, Jm=Jm, progress=progress)
    res['march'] = mr['rec']
    Y = y0_columns(op, mr)
    S = Ser(p, keepN)
    T = Ntot + 1
    # commutant
    Nfrac = sb['N_struct']
    den = 1
    for row in Nfrac:
        for x in row:
            den = den * x.denominator // math.gcd(den, x.denominator)
    Nint = [[int(x * den) for x in row] for row in Nfrac]   # N scaled by den (nilpotent scaling ok)
    Ebasis = commutant_basis([[Fraction(x, den) for x in row] for row in Nint], p, n_)
    res['E_space_dim'] = len(Ebasis)
    # M~ = Y0(x^p) * diag(x^{-p deff}); leading matrix Lambda exact over Q
    Lam = [[None] * n_ for _ in range(n_)]
    for j in range(n_):
        dj = sb['deff'][j]
        for i in range(n_):
            vec = sb['basis'][j].get(dj, [Fraction(0)] * J)
            acc = Fraction(0)
            for t in range(min(i, J - 1) + 1):
                acc += math.comb(i, t) * Fraction(dj) ** (i - t) * vec[t]
            Lam[i][j] = acc
    detLam = None
    LamInv = matinvQ(Lam)
    mod = S.mod
    LamInv_int = [[int(r2.mpz(x.numerator % mod) * r2.INV(r2.mpz(x.denominator % mod), mod) % mod)
                   for x in row] for row in LamInv]
    for row in LamInv:
        for x in row:
            dv = valp(x.denominator, p)
            assert not dv or dv <= sigma2 // 2, ('Lambda inverse denominator too deep', dv)
    # build M~ columns: (Y0[i][j] shifted down by deff_j) composed x->x^p
    Mt = [[None] * n_ for _ in range(n_)]
    for j in range(n_):
        dj = sb['deff'][j]
        for i in range(n_):
            base = Y[i][j][dj:]
            Mt[i][j] = S.mk(S.subs_xp(base, T))
    # Newton inversion (stored scale sigma; lead stored = p^sigma*Lam):
    # feed lead_inv = p^{-sigma}*LamInv ... work with TRUE-scale objects instead:
    # Divide M~ stored by p^sigma is not exact entrywise; instead invert scaled:
    # X_true = M~_true^{-1}; stored X = p^{sigma2} X_true. Newton on stored with
    # correction: X <- (X(2*p^{sigma+sigma2} I - M~ X)) / p^{sigma+sigma2}, exact-divide.
    ps12 = p ** (sigma + sigma2)
    X = [[S.mk(int(r2.mpz(LamInv_int[i][j]) * (p ** sigma2) % mod)) for j in range(n_)]
         for i in range(n_)]
    t = 1
    while t < T:
        t = min(2 * t, T)
        MX = mat_mul_ser(S, Mt, X, t)
        R_ = [[(S.mk(2 * ps12 if i == j else 0) - MX[i][j]) for j in range(n_)] for i in range(n_)]
        XR = mat_mul_ser(S, X, R_, t)
        # exact divide by p^{sigma+sigma2}
        Xn = []
        for i in range(n_):
            rown = []
            for j in range(n_):
                cs = [int(XR[i][j][u]) for u in range(min(XR[i][j].degree() + 1, t))]
                out = []
                for cval in cs:
                    assert cval % ps12 == 0, 'sigma2 undersized in Newton divide-back'
                    out.append(cval // ps12)
                rown.append(S.mk(out))
            Xn.append(rown)
        X = Xn
    # CERTIFY the inverse: X*Mt == p^{s12} I; the realized min-valuation of the
    # residual is the certificate consumed downstream (budget != certificate).
    chk = mat_mul_ser(S, X, Mt, T)
    cert_inv = None
    for i in range(n_):
        for j in range(n_):
            for u in range(T):
                cv = int(chk[i][j][u]) - (ps12 if (i == j and u == 0) else 0)
                if cv % S.mod:
                    v = valp(cv % S.mod, p) or 0
                    cert_inv = v if cert_inv is None else min(cert_inv, v)
    if cert_inv is None:
        cert_inv = keepN
    res['inv_residual_certified_depth'] = cert_inv
    assert cert_inv >= keep, ('Newton certificate below target', cert_inv, keep)
    # F_b rows: G_b = Y0*E_b ; P_b col k shifted by p*(maxdeff - deff_k); F_b = P_b X / x^{p maxdeff}
    shift = [p * (maxdeff - d) for d in sb['deff']]
    offset = p * maxdeff       # stored index n corresponds to true exponent n - offset
    def assemble_rows(E, rows):
        G = [[None] * n_ for _ in range(len(rows))]
        for ri, i in enumerate(rows):
            for k in range(n_):
                coeffs = [0] * T
                for jj in range(n_):
                    e_ = E[jj][k]
                    if e_:
                        yv = Y[i][jj]
                        for u in range(T):
                            if yv[u]:
                                coeffs[u] = (coeffs[u] + e_ * yv[u]) % S.mod
                # apply shift for column k
                sh = shift[k]
                coeffs = [0] * sh + coeffs[:T - sh]
                G[ri][k] = S.mk(coeffs)
        return mat_mul_ser(S, G, X, T)
    # ---- E pinning: the true Frobenius ray is selected by
    # DECAY + regularity + integrality.  A wrong commutant direction stays
    # bounded-but-non-decaying (its F does not converge at |x|=1); the true F
    # is overconvergent: val(c_n) grows ~ (n - n0)/(C p).  Depth ladder:
    # rows demand val(true c_idx) >= d over a tail window, d = 1,2,...; the
    # largest d with a surviving ray gives the pinned E and the MEASURED
    # discrimination depth (gate GE data).  Regularity rows (true exponent
    # < 0 must vanish) and integrality rows (d = 0) always included.
    # ---- twist product Ptw = prod Delta_i(x^p)^{mult_i * s_twist}
    Ptw = S.mk(1)
    for coeffs, name, mult in twist:
        fac = [0] * ((len(coeffs) - 1) * p + 1)
        for ii, cc in enumerate(coeffs):
            fac[ii * p] = cc % S.mod
        fp = S.mk(fac)
        acc = S.mk(1)
        e_ = s_twist * mult
        base = fp
        while e_:
            if e_ & 1:
                acc = acc * base
            e_ >>= 1
            if e_:
                base = base * base
        Ptw = Ptw * acc
    Ptw = S.mk([Ptw[u] for u in range(min(Ptw.degree() + 1, T))])
    # ---- E pinning: depth ladder on the TWISTED tails (the singular-locus
    # denominator law; s_twist and the window are the measured N ~ C*p*m
    # constants).  Regularity rows (true exponent < 0)
    # always included.
    pinrows_idx = list(range(min(n_, 3)))
    Frows = [assemble_rows(E, pinrows_idx) for E in Ebasis]
    KE = keep - sigma - sigma2
    nb = len(Ebasis)
    FrowsT = [[[Frows[b][ri][j] * Ptw for j in range(n_)]
               for ri in range(len(pinrows_idx))] for b in range(nb)]
    # rows fed RAW (stored scale s12 = sigma+sigma2); the valuation-aware
    # kernel handles per-row scales.  Ladder modulus for depth d: p^{s12+d}
    # (condition: val(true combo) >= d).  Deep negative valuations of the
    # basis elements (the wrong directions diverge at the ledger rate) are
    # measured and reported — they are the discrimination mechanism.
    s12v = sigma + sigma2
    minv = None
    for b in range(nb):
        for ri in range(len(pinrows_idx)):
            for j in range(n_):
                poly = FrowsT[b][ri][j]
                for idx in range(T - 220, T, 4):
                    cv = int(poly[idx])
                    if cv:
                        v = valp(cv, p) or 0
                        minv = v if minv is None else min(minv, v)
    dmax_try = KE - 4
    modcap = p ** (s12v + dmax_try)
    base_rows = []
    stride_neg = max(1, offset // 48) if offset else 1
    for ri in range(len(pinrows_idx)):
        for j in range(n_):
            for idx in range(0, offset, stride_neg):
                row = [int(Frows[b][ri][j][idx]) % modcap for b in range(nb)]
                if any(row):
                    base_rows.append(row)
    wlen = max(24, min(220, (T - offset) // 6))
    tail_rows = []
    for ri in range(len(pinrows_idx)):
        for j in range(n_):
            for idx in range(T - wlen, T - 4, 2):
                row = [int(FrowsT[b][ri][j][idx]) % modcap for b in range(nb)]
                if any(row):
                    tail_rows.append(row)
    allrows = base_rows + tail_rows
    ker = None
    d_star = None
    lo, hi = -(s12v - 32), dmax_try      # deep negative floor: the ray may sit
    while lo <= hi:                      # well below the jet-lattice baseline
        mid = (lo + hi) // 2
        k_ = padic_kernel(allrows, p, s12v + mid)
        if k_['gens']:
            ker, d_star = k_, mid
            lo = mid + 1
        else:
            hi = mid - 1
    assert ker is not None, 'pinning: no ray survives even at the negative floor — design wall'
    gens = ker['gens']
    res['twist'] = {'factors': [(nm, mu) for _, nm, mu in twist], 's_twist': s_twist,
                    'degD_weighted': degD, 'NF': NF}
    res['pin'] = {'n_base_rows': len(base_rows), 'n_tail_rows': len(tail_rows),
                  'tail_window': [T - wlen, T - 4], 'depth_ladder_dstar': d_star,
                  'stored_scale_s12': s12v, 'min_stored_val_measured': minv,
                  'pivots': ker['pivots'], 'free_cols': ker['free_cols'],
                  'obstructed': ker['obstructed'], 'n_gens': len(gens),
                  'KE_available': KE, 'ray_unique_up_to_scalar': len(gens) == 1}
    return {'res': res, 'sb': sb, 'mr': mr, 'S': S, 'Y': Y, 'X': X, 'Ebasis': Ebasis,
            'gens': gens, 'offset': offset, 'shift': shift, 'KE': KE, 'T': T,
            'ps12': ps12, 'assemble_rows': assemble_rows, 'p': p, 'n': n_,
            'keep': keepN, 'keep_target': keep, 'sigma': sigma, 'sigma2': sigma2,
            'op': op, 'Ptw': Ptw, 'twist': twist, 's_twist': s_twist,
            'd_star': d_star, 't_wall': time.time() - t0}

def pinned_E(tr):
    """The pinned integer-combination E matrix (ray representative) from the
    nullspace generator; asserts a single free direction."""
    gens, Eb = tr['gens'], tr['Ebasis']
    assert len(gens) >= 1
    e = gens[0]
    n_ = tr['n']
    E = [[sum(e[b] * Eb[b][i][j] for b in range(len(Eb))) for j in range(n_)]
         for i in range(n_)]
    return E, e

def normalize_ray(tr, E, F):
    """Rescale the ray representative so the assembled F is p-integral at the
    stored scale: measure baseline = min(true val) over F entries; if
    negative, scale E and F by p^{-baseline}.  gamma absorbs the change
    (gamma_branches carries valuation).  Returns (E', F', baseline)."""
    S, T, n_, p = tr['S'], tr['T'], tr['n'], tr['p']
    s12 = tr['sigma'] + tr['sigma2']
    vmin = None
    for i in range(n_):
        for j in range(n_):
            poly = F[i][j]
            for idx in range(0, T, max(1, T // 400)):
                cv = int(poly[idx])
                if cv:
                    v = (valp(cv, p) or 0) - s12
                    vmin = v if vmin is None else min(vmin, v)
    baseline = min(0, vmin if vmin is not None else 0)
    if baseline < 0:
        sc = p ** (-baseline)
        E = [[x * sc for x in row] for row in E]
        F = [[F[i][j] * S.mk(sc) for j in range(n_)] for i in range(n_)]
    return E, F, baseline

def assemble_full_F(tr, E):
    """Full 6x6 (n x n) F series for the pinned E: stored scale s12, Laurent
    offset tr['offset'].  Returns list-of-lists of flint polys."""
    S, T, n_ = tr['S'], tr['T'], tr['n']
    Y, X, shift = tr['Y'], tr['X'], tr['shift']
    G = [[None] * n_ for _ in range(n_)]
    for i in range(n_):
        for k in range(n_):
            coeffs = [0] * T
            for jj in range(n_):
                e_ = E[jj][k]
                if e_:
                    yv = Y[i][jj]
                    for u in range(T):
                        if yv[u]:
                            coeffs[u] = (coeffs[u] + e_ * yv[u]) % S.mod
            sh = shift[k]
            coeffs = [0] * sh + coeffs[:T - sh]
            G[i][k] = S.mk(coeffs)
    return mat_mul_ser(S, G, X, T)

def gd_residual(tr, F, check_len=None):
    """GATE GD: c6(x) c6(x^p) theta F - c6(x^p) C(x) F + p F C(x^p) c6(x) == 0
    (C = c6 * A_th companion).  Verifies the whole chain.  Returns max
    valuation-deficit found (None = clean at available precision)."""
    op, S, T, n_, p = tr['op'], tr['S'], tr['T'], tr['n'], tr['p']
    offset = tr['offset']
    ci = op.theta_form()
    c6 = ci[op.order]
    c6f = S.mk([x % S.mod for x in c6])
    c6p = S.mk(S.subs_xp([x % S.mod for x in c6], T))
    check_len = check_len or min(T, 4 * p)
    # theta F on stored (true exponent = idx - offset): coeff *= (idx - offset)
    thF = [[S.mk([(idx - offset) * int(F[i][j][idx]) % S.mod for idx in range(T)])
            for j in range(n_)] for i in range(n_)]
    # C(x) F: rows 0..n-2 of C: c6 * shift; row n-1: -c_i
    CF = [[None] * n_ for _ in range(n_)]
    for i in range(n_ - 1):
        for j in range(n_):
            CF[i][j] = c6f * F[i + 1][j]
    cif = [S.mk([x % S.mod for x in ci[k]]) for k in range(op.order)]
    for j in range(n_):
        acc = S.mk(0)
        for k in range(n_):
            acc = acc - cif[k] * F[k][j]
        CF[n_ - 1][j] = acc
    # F C(x^p): columns: (F C)_{ij} = F_{i,j-1}*c6(x^p) [j>=1] + F_{i,n-1} * (-c_j(x^p))
    cipf = [S.mk(S.subs_xp([x % S.mod for x in ci[k]], T)) for k in range(op.order)]
    FC = [[None] * n_ for _ in range(n_)]
    for i in range(n_):
        for j in range(n_):
            acc = S.mk(0)
            if j >= 1:
                acc = acc + F[i][j - 1] * c6p
            acc = acc - F[i][n_ - 1] * cipf[j]
            FC[i][j] = acc
    worst = None
    modv = S.mod
    for i in range(n_):
        for j in range(n_):
            R_ = c6f * (c6p * thF[i][j]) - c6p * CF[i][j] + p * FC[i][j] * c6f
            for idx in range(min(check_len, T)):
                cv = int(R_[idx]) % modv
                if cv:
                    v = valp(cv, tr['p']) or 0
                    deficit = tr['keep'] - v
                    worst = max(worst or 0, deficit)
    return worst

def eval_F_at(tr, F, x0, tail_chunk=24):
    """TWISTED evaluation of stored F at the Teichmuller point of fiber x0:
    sum the decaying series Ptw*F, divide by prod Delta_i(xhat)^{s} (a unit at
    good fibers; Teichmuller makes Delta(xhat^p) = Delta(xhat)).  Returns
    (Fmat ints true-value mod p^{cert}, cert, tailval): cert = min(KE, d_star,
    measured tail valuation of the twisted series)."""
    S, T, n_, p = tr['S'], tr['T'], tr['n'], tr['p']
    offset, ps12 = tr['offset'], tr['ps12']
    keep = tr['keep']
    Ptw, twist, s_twist = tr['Ptw'], tr['twist'], tr['s_twist']
    modv = p ** keep
    xh = teich(x0, p, keep)
    xinv = pow(xh, -1, modv)
    den = 1
    for coeffs, name, mult in twist:
        val = sum(c * pow(xh, ii, modv) for ii, c in enumerate(coeffs)) % modv
        assert val % p != 0, ('twist factor %s vanishes at fiber — bad fiber' % name, x0)
        den = den * pow(val, s_twist * mult, modv) % modv
    dinv = pow(den, -1, modv)
    tailv = None
    Fm = [[0] * n_ for _ in range(n_)]
    xh_powers = [1] * T
    for idx in range(1, T):
        xh_powers[idx] = xh_powers[idx - 1] * xh % modv
    for i in range(n_):
        for j in range(n_):
            G = F[i][j] * Ptw
            acc = 0
            for idx in range(T):
                cv = int(G[idx])
                if cv:
                    acc = (acc + cv * xh_powers[idx]) % modv
            for idx in range(T - tail_chunk, T):
                cv = int(G[idx])
                if cv:
                    v = (valp(cv, p) or 0) - (tr['sigma'] + tr['sigma2'])
                    tailv = v if tailv is None else min(tailv, v)
            acc = acc * dinv % modv
            acc = acc * pow(xinv, offset, modv) % modv
            assert acc % ps12 == 0, ('F(xhat) not integral at stored scale', i, j,
                                     valp(acc, p))
            Fm[i][j] = acc // ps12
    if tailv is None:
        tailv = tr['KE']
    cert = max(1, min(tr['KE'], tr.get('d_star', tr['KE']), tailv))
    Fm = [[x % (p ** cert) for x in row] for row in Fm]
    return Fm, cert, tailv

def gamma_branches(E0, p, wt, n, Kg):
    """Scalar normalization from the INTEGER det of the pinned ray matrix
    (det F(xhat) = det E exactly at Teichmuller points since
    the Wronskian-monomial x^{sum deff (1-p)} evaluates to 1):
    gamma^n * det(E0) = sgn * p^{n*wt/2}.  Returns [(gamma mod p^Kg, tag)...]
    over sgn in {+1,-1} and the n-th roots of unity in Z_p."""
    # integer determinant of E0
    import itertools
    n_ = len(E0)
    det = 0
    for perm in itertools.permutations(range(n_)):
        sgn = 1
        seen = list(perm)
        # parity via inversion count
        inv = sum(1 for a in range(n_) for b in range(a + 1, n_) if perm[a] > perm[b])
        term = (-1) ** inv
        for i2 in range(n_):
            term *= E0[i2][perm[i2]]
        det += term
    vd = valp(det, p)
    out = {'detE0': det if abs(det) < 10 ** 40 else 'big(val %s)' % vd,
           'val_detE0': vd, 'target_val': n * wt // 2, 'branches': []}
    if vd is None:
        return out
    # gamma^n det = sgn p^{n wt/2}: val(gamma) = (n wt/2 - vd)/n, must be integral
    dv = n * wt // 2 - vd
    if dv % n != 0:
        out['gamma_val_nonintegral'] = dv
        return out
    vgam = dv // n
    out['gamma_valuation'] = vgam
    modg = p ** Kg
    u = (det // p ** vd) % modg
    for sgn in (1, -1):
        rhs = sgn * pow(u, -1, modg) % modg
        r0 = rhs % p
        roots = [r for r in range(1, p) if pow(r, n, p) == r0]
        for r in roots:
            x = r
            pk = p
            while pk < modg:
                pk = min(pk * pk, modg)
                fx = (pow(x, n, pk) - rhs) % pk
                fpx = n * pow(x, n - 1, pk) % pk
                x = (x - fx * pow(fpx, -1, pk)) % pk
            out['branches'].append({'gamma_unit': int(x), 'gamma_val': vgam,
                                    'sign': sgn, 'root_mod_p': r})
    return out

def apply_gamma(c, br, p, K):
    """Normalized R_p coefficients: c_k * gamma^k with gamma = p^{v} * unit.
    Returns list (ints mod p^{K - max(0, -v*n)}), or None if a negative-val
    gamma makes a coefficient non-integral (branch dies)."""
    mod = p ** K
    g, v = br['gamma_unit'], br['gamma_val']
    out = []
    for k, ck in enumerate(c):
        val = ck * pow(g, k, mod) % mod
        if v >= 0:
            val = val * pow(p, v * k, mod) % mod
        else:
            pv = p ** (-v * k)
            if val % pv:
                return None
            val //= pv
        out.append(val)
    return out

def gamma_normalize(cp_asc, p, K, wt, n):
    """Solve the unit scalar gamma from det + FE emergence:
    true charpoly of gamma*F.  Returns list of (gamma mod p^Kg, branch tag);
    empty list = normalization impossible (gate fail).  n even."""
    mod = p ** K
    c = rp_coeffs(cp_asc, p, K)     # c_k of det(1-TF)
    half = n // 2
    out = []
    # det F = (-1)^n * cp_asc[0] = c_n; want gamma^n c_n = ±p^{n wt/2}
    cn = c[n] % mod
    vdet = valp(cn, p)
    if vdet is None or vdet != n * wt // 2 or K <= n * wt // 2:
        return out
    u6 = (cn // p ** vdet) % (p ** (K - vdet))
    # FE pair: c_{half+1} = p^{wt} c_{half-1}: gamma^{half+1} c_{half+1} = p^wt gamma^{half-1} c_{half-1}
    # => gamma^2 = p^wt c_{half-1} / c_{half+1}
    ch1 = c[half + 1] % mod
    chm = c[half - 1] % mod if half >= 1 else 1
    v1 = valp(ch1, p)
    if v1 is None or v1 < wt or (chm % p == 0 and half - 1 >= 1):
        return out
    Kg = K - max(v1, wt) - 1
    if Kg < 3:
        return out
    modg = p ** Kg
    rhs = p ** wt * chm
    g2 = (rhs // p ** v1) * pow(ch1 // p ** v1, -1, modg) % modg if v1 <= wt and rhs % p ** v1 == 0 else None
    if g2 is None or g2 % p == 0:
        return out
    # sqrt mod p^Kg (p odd)
    a = g2 % p
    if pow(a, (p - 1) // 2, p) != 1:
        return out
    r = pow(a, (p + 1) // 4, p) if p % 4 == 3 else None
    if r is None:
        # Tonelli-Shanks mod p
        q, s = p - 1, 0
        while q % 2 == 0:
            q //= 2
            s += 1
        z = 2
        while pow(z, (p - 1) // 2, p) != p - 1:
            z += 1
        mfac, cq, tq, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
        while tq != 1:
            i2, t2 = 0, tq
            while t2 != 1:
                t2 = t2 * t2 % p
                i2 += 1
            b = pow(cq, 1 << (mfac - i2 - 1), p)
            mfac, cq, tq, r = i2, b * b % p, tq * b * b % p, r * b % p
    # Hensel lift sqrt to modg
    x = r
    pk = p
    while pk < modg:
        pk = min(pk * pk, modg)
        x = (x - (x * x - g2) * pow(2 * x, -1, pk)) % pk
    for gam in (x % modg, (-x) % modg):
        # consistency with det: gamma^n * u6-part
        gn = pow(gam, n, modg)
        lhsu = gn * (u6 % modg) % modg
        for sgn in (1, -1):
            if (lhsu - sgn) % modg == 0 or (lhsu - sgn) % (p ** max(1, Kg - 2)) == 0:
                out.append((int(gam), 'sqrt%+d_det%+d' % (1 if gam == x % modg else -1, sgn)))
    return out
