#!/usr/bin/env python3
"""lfact_kit.py — instrument kit for exact local factors of operator motives
(single-series wing of the deform_transport member).

Consumes the vendored r2_kit engine by sibling import (sha pinned in every
receipt).  Instrument members, each carrying its own control:
  - multi_level_lift: unit_root_lift extended to snapshots at
    m = s + p^k - 1 for k = 1..kmax (accumulators mod p^{kmax+1}, Teichmuller
    at matching depth, guard S = kmax + L + 8).  Engine letter otherwise the
    single-series one (FD registers, seed law, loss ledger, shrink, events).
  - sym5_control: exact Sym^5-Legendre rank-6 control for the truncation-
    snapshot algebra (closed-form 2F1 coefficients; NO division steps), with
    exact Frobenius eigenvalues from the character-sum a_p.
  - hankel_reads: u_k = T_{k+1}/T_k and v1 = (T_3 T_1 - T_2^2)/(p (T_2 - T_1^2))
    normalized by u — the slope-1 unit read (EMPIRICAL register; licensed
    only by the sym5 control's measured precision).
  - pcurv_jordan: independent jet-based p-curvature Jordan profile at one
    fiber (an independent code path from any series route).
  - enumerate_sextics: Weil-box + p^4-functional-equation integer isolation
    under declared factor shapes S0/S1/S2.

The order-6 reference operator (L6-class) is NOT distributed with the repo:
set FROB_L6_OPERATOR to a JSON operator file ({"coeffs": [[...], ...]}) to
use the L6-dependent members (l6_path() refuses loudly otherwise).
"""
import importlib.util, json, hashlib, os, sys, time
from fractions import Fraction

_HERE = os.path.dirname(os.path.abspath(__file__))

def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for ch in iter(lambda: f.read(1 << 20), b''):
            h.update(ch)
    return h.hexdigest()

def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod

# vendored engines, imported from beside this file (shas in every receipt)
r2 = load_module('r2_kit', os.path.join(_HERE, 'r2_kit.py'))
f2 = load_module('f2_lib', os.path.join(_HERE, 'f2_lib.py'))
R2_SHA = sha256_file(os.path.join(_HERE, 'r2_kit.py'))
F2_SHA = sha256_file(os.path.join(_HERE, 'f2_lib.py'))

def l6_path():
    """Path to the order-6 reference operator JSON (env-gated; loud refusal)."""
    path = os.environ.get('FROB_L6_OPERATOR')
    if not path or not os.path.exists(path):
        raise RuntimeError(
            'lfact_kit: the order-6 reference operator is not distributed with '
            'the repo; set FROB_L6_OPERATOR to a JSON operator file '
            '({"coeffs": [[...], ...]}, 7 rows) to use the L6-dependent members')
    return path

def l6_sha():
    return sha256_file(l6_path())

mpz, INV = r2.mpz, r2.INV  # exact backend (gmpy2 when available)

def utc():
    import datetime
    return datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')

def producer(script):
    eng = {'r2_kit.py': R2_SHA, 'f2_lib.py': F2_SHA,
           'lfact_kit.py': sha256_file(os.path.abspath(__file__))}
    if os.environ.get('FROB_L6_OPERATOR') and os.path.exists(os.environ['FROB_L6_OPERATOR']):
        eng['l6_operator'] = l6_sha()
    return {'component': 'frobenius-boundary/deform_transport (lfact_kit)',
            'script': script,
            'script_sha256': sha256_file(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                      os.path.basename(script))),
            'stamp_utc': utc(),
            'engine_shas': eng}

# ------------------------------------------------------------------ engine

def multi_level_lift(A, e, p, seeds, seed_upto, xgrid, lead_exps, kmax,
                     progress=None):
    """unit_root_lift letter, extended: snapshots at m = s + p^k - 1,
    k = 1..kmax; accumulators and Teichmuller mod p^{kmax+1}; guard
    S = kmax + L(N) + 8; N = p^kmax - 1 + max(lead_exps).
    Returns receipt with per-fiber per-column T-list [T_1..T_kmax] mod
    p^{kmax+1} (G-normalized: xhat^{-s} applied), plus raw_modp, events,
    ledger fields — same law as the single-series engine."""
    t0 = time.time()
    A7, D, r, s0, dmax, Wk = r2.op_scaffold(A, e)
    ncol = len(seeds)
    smax = max(lead_exps)
    N = p ** kmax - 1 + smax
    m_start = seed_upto + 1
    Ltot, sing = r2.loss_ledger(A, e, p, N, m_start)
    S = kmax + Ltot + 8
    pk1 = p ** (kmax + 1)
    rec = {"p": p, "N": N, "kmax": kmax, "S_alloc": S, "L_ledger": Ltot,
           "n_singular_steps": len(sing), "e": e, "backend": r2.BACKEND,
           "seed_upto": seed_upto, "lead_exps": lead_exps}
    seedvals = []
    for col in seeds:
        vals = []
        for m in range(seed_upto + 1):
            fr = col.get(m, Fraction(0))
            if fr.denominator % p == 0:
                rec["seed_denominator_hit"] = [m, str(fr)]
                return rec
            vals.append(fr)
        seedvals.append(vals)
    suffix_after = {}
    run = 0
    for i in range(len(sing) - 1, -1, -1):
        suffix_after[sing[i][0]] = run
        run += sing[i][1]
    M = mpz(p ** S)
    regs = None
    def fd_init(m0, mod):
        regs = []
        for k in range(dmax + 1):
            exact = [Wk(k, m0 - k + t) for t in range(r + 1)]
            for lev in range(1, r + 1):
                for t in range(r, lev - 1, -1):
                    exact[t] = exact[t] - exact[t - 1]
            top = r
            while top > 0 and exact[top] == 0:
                top -= 1
            regs.append([mpz(z % mod) for z in exact[:top + 1]])
        return regs
    regs = fd_init(m_start - 1, M)
    window = []
    for m in range(seed_upto + 1):
        rowv = []
        for c in range(ncol):
            fr = seedvals[c][m]
            rowv.append(mpz(fr.numerator) * INV(fr.denominator, M) % M)
        window.append(rowv)
    def fib_parts(x0):
        return (x0, 1) if isinstance(x0, int) else x0
    def fib_key(x0):
        n, d = fib_parts(x0)
        return str(n) if d == 1 else "%d/%d" % (n, d)
    # Teichmuller at depth kmax+1: Teich(x) = x^{p^kmax} mod p^{kmax+1}
    xhat, xres = {}, {}
    for x0 in xgrid:
        n, d = fib_parts(x0)
        tn = pow(n, p ** kmax, pk1)
        td = pow(d, p ** kmax, pk1)
        xhat[x0] = tn * pow(td, -1, pk1) % pk1
        xres[x0] = n * pow(d, -1, p) % p
        assert pow(xhat[x0], p, pk1) == xhat[x0]  # fixed point law
    accs = {x0: [0] * ncol for x0 in xgrid}
    snaps = {x0: [[None] * (kmax + 1) for _ in range(ncol)] for x0 in xgrid}
    xpow = {x0: 1 for x0 in xgrid}
    raw_modp = {x0: [0] * ncol for x0 in xgrid}
    x0pow = {x0: 1 for x0 in xgrid}
    snap_at = {}
    for c in range(ncol):
        for k in range(1, kmax + 1):
            snap_at.setdefault(lead_exps[c] + p ** k - 1, []).append((c, k))
    for m in range(seed_upto + 1):
        for x0 in xgrid:
            for c in range(ncol):
                accs[x0][c] = (accs[x0][c] + int(window[m][c] % pk1) * xpow[x0]) % pk1
                if m < p:
                    raw_modp[x0][c] = (raw_modp[x0][c] + int(window[m][c] % p) * x0pow[x0]) % p
        for x0 in xgrid:
            xpow[x0] = xpow[x0] * xhat[x0] % pk1
            if m < p:
                x0pow[x0] = x0pow[x0] * xres[x0] % p
    loss = 0
    events = []
    alive = [True] * ncol
    last_shrink_loss = 0
    check_every = 4096
    for m in range(m_start, N + 1):
        kmaxd = min(m, dmax)
        for reg in regs:
            for lev in range(len(reg) - 1):
                reg[lev] = (reg[lev] + reg[lev + 1]) % M
        if m % check_every == 0:
            w0d = Wk(0, m) % M
            assert regs[0][0] == w0d, (m, 'FD drift')
        pivot = regs[0][0]
        rhs = [0] * ncol
        for k in range(1, kmaxd + 1):
            wkv = regs[k][0]
            if wkv:
                row = window[m - k]
                for c in range(ncol):
                    if alive[c] and row[c] is not None:
                        rhs[c] += wkv * row[c]
        piv = pivot % M
        v = 0
        pw = piv
        while pw % p == 0 and v < 8:
            pw //= p
            v += 1
        if pw % p == 0:
            events.append([m, "pivot_val_gt8"])
            rec["events"] = events
            rec["ABORT"] = "pivot valuation > 8"
            return rec
        if v == 0:
            inv = INV(piv, M)
            newrow = [None if not alive[c] else (-rhs[c] % M) * inv % M
                      for c in range(ncol)]
        else:
            pv = p ** v
            unit = piv // pv
            inv = INV(unit, M)
            newrow = [None] * ncol
            for c in range(ncol):
                if not alive[c]:
                    continue
                rc = (-rhs[c]) % M
                if rc % pv != 0:
                    events.append({"m": m, "col_lead_exp": lead_exps[c],
                                   "event": "integrality_violation",
                                   "pivot_val": v, "rhs_val": r2.valp(rc, p)})
                    alive[c] = False
                    continue
                newrow[c] = (rc // pv) * inv % M
            loss += v
            if loss - last_shrink_loss >= 48:
                Srem = kmax + suffix_after.get(m, 0) + 8
                Mnew = p ** Srem
                if Mnew < M:
                    M = Mnew
                    for row in window:
                        if row is not None:
                            for c in range(ncol):
                                if row[c] is not None:
                                    row[c] %= M
                    for k in range(dmax + 1):
                        regs[k] = [z % M for z in regs[k]]
                    newrow = [None if z is None else z % M for z in newrow]
                last_shrink_loss = loss
        window.append(newrow)
        for x0 in xgrid:
            xp = xpow[x0]
            av = accs[x0]
            for c in range(ncol):
                if newrow[c]:
                    av[c] = int(av[c] + int(newrow[c] % pk1) * xp) % pk1
            xpow[x0] = xp * xhat[x0] % pk1
            if m < p:
                for c in range(ncol):
                    if newrow[c] is not None:
                        raw_modp[x0][c] = int(raw_modp[x0][c]
                                              + int(newrow[c] % p) * x0pow[x0]) % p
                x0pow[x0] = x0pow[x0] * xres[x0] % p
        if m in snap_at:
            for (c, k) in snap_at[m]:
                if alive[c]:
                    for x0 in xgrid:
                        snaps[x0][c][k] = accs[x0][c]
        if len(window) > dmax + 1 and m > m_start + dmax + 2:
            window[m - dmax - 1] = None
        if progress and m % progress == 0:
            print("  m=%d/%d loss=%d/%d t=%.1fs" % (m, N, loss, Ltot,
                                                    time.time() - t0), flush=True)
    rec["alive_at_end"] = {("s%d" % lead_exps[c]): alive[c] for c in range(ncol)}
    if events:
        rec["events"] = events
    rec["loss_realized"] = loss
    rec["final_precision_digits"] = S - loss
    rec["precision_ok"] = (S - loss) >= kmax + 1
    assert loss == Ltot, (loss, Ltot, "ledger mismatch")
    fibers = {}
    for x0 in xgrid:
        fx = {"xhat_mod_pk1": int(xhat[x0])}
        xr = xres[x0]
        badf = [q for q in (xr, xr - 1, xr + 1, xr**2 + xr + 1,
                            xr**4 + xr**3 + xr**2 + xr + 1) if q % p == 0]
        if badf:
            fx["skip"] = "singular fiber mod p"
            fibers[fib_key(x0)] = fx
            continue
        fx["raw_modp"] = [int(z) for z in raw_modp[x0]]
        cols = {}
        for c in range(ncol):
            s = lead_exps[c]
            xs_inv = pow(pow(xhat[x0], s, pk1), -1, pk1)
            Ts = snaps[x0][c]
            if any(Ts[k] is None for k in range(1, kmax + 1)):
                cols["s%d" % s] = {"non_integral": True,
                                   "T_partial": [None if t is None else int(t)
                                                 for t in Ts[1:]]}
                continue
            T = [1] + [int(Ts[k]) * xs_inv % pk1 for k in range(1, kmax + 1)]
            d = {"T": [int(t) for t in T]}   # T[0]=1, T[k] = [G]_{<p^k}(xhat)
            d["T_vals"] = [(r2.valp(t, p) or 0) if t else None for t in T]
            if T[1] % p == 0:
                d["ordinary"] = False
            else:
                d["ordinary"] = True
                # u_k = T_{k+1}/T_k mod p^{k+1}, k = 1..kmax-1 (unit T_k only)
                for k in range(1, kmax):
                    mod = p ** (k + 1)
                    if T[k] % p == 0:
                        d["u_level%d" % k] = None
                        d["u_level%d_note" % k] = "T_%d nonunit" % k
                    else:
                        d["u_level%d" % k] = int(T[k + 1] * pow(T[k], -1, mod) % mod)
            cols["s%d" % s] = d
        fx["cols"] = cols
        fibers[fib_key(x0)] = fx
    rec["fibers"] = fibers
    rec["t_wall_s"] = round(time.time() - t0, 2)
    return rec

# --------------------------------------------------------------- sym5 control

def sym5_series_modpk(p, kmax, Nlen, mod):
    """f0 = 2F1(1/2,1/2;1;lam) coefficients c_n = binom(2n,n)^2/16^n mod `mod`
    (EXACT integers via math.comb; no divisions by p anywhere: 16^{-1} exists),
    then d = f0^5 truncated to Nlen terms mod `mod` (int convolution)."""
    import numpy as np
    from scipy.signal import fftconvolve
    assert mod < 1 << 24, 'limb-split convolution guard (mod = p^{kmax+1} <= 61^4)'
    inv16 = pow(16, -1, mod)
    c = np.zeros(Nlen, dtype=np.int64)
    b = 1  # binom(2n, n), exact integer recurrence
    ipow = 1
    for n in range(Nlen):
        if n:
            b = b * (2 * (2 * n - 1)) // n     # binom(2n,n) = binom(2n-2,n-1)*2(2n-1)/n
            ipow = ipow * inv16 % mod
        c[n] = (b % mod) * ipow % mod * (b % mod) % mod
    B13 = (1 << 13) % mod
    B26 = (1 << 26) % mod
    def mul(a, b2):
        # exact mod-mul via base-2^13 limb split + float FFT (error << 1, values <= 2^44)
        alo, ahi = (a & 8191).astype(np.float64), (a >> 13).astype(np.float64)
        blo, bhi = (b2 & 8191).astype(np.float64), (b2 >> 13).astype(np.float64)
        ll = np.rint(fftconvolve(alo, blo)[:Nlen]).astype(np.int64)
        cr = np.rint(fftconvolve(alo, bhi)[:Nlen] + fftconvolve(ahi, blo)[:Nlen]).astype(np.int64)
        hh = np.rint(fftconvolve(ahi, bhi)[:Nlen]).astype(np.int64)
        out = (ll % mod + (cr % mod) * B13 % mod + (hh % mod) * B26 % mod) % mod
        return out.astype(np.int64)
    d2 = mul(c, c)      # f^2
    d4 = mul(d2, d2)    # f^4
    d = mul(d4, c)      # f^5
    return [int(v) for v in d]

def sym5_control(p, fibers, kmax=3):
    """Exact rank-6 control at prime p: truncation snapshots of f0^5 at
    Teichmuller points vs the EXACT Sym^5 Frobenius eigenvalues alpha^i beta^{5-i}
    (alpha = unit root of T^2 - a_p T + p; here directly the Legendre-lambda
    charsum via f2.legendre_ap_fast).
    Measures: u-read precision per level; Hankel v1-read precision.
    Returns receipt rows; EXPECTED (empirical model): u_k == alpha^5 mod p^{k+1};
    v1 = Hankel/(p*u) == alpha^3 mod p^{m1} with m1 measured here."""
    pk1 = p ** (kmax + 1)
    Nlen = p ** kmax
    d = sym5_series_modpk(p, kmax, Nlen, pk1)
    chi = f2.chi_table(p)
    rows = {}
    for lam0 in fibers:
        lam0 = lam0 % p
        if lam0 in (0, 1):
            continue
        ap = f2.legendre_ap_fast(lam0, p, chi)
        if ap % p == 0:
            rows[str(lam0)] = {'skip': 'supersingular'}
            continue
        # exact unit root alpha of E_lam: root of T^2 - ap T + p with T=ap mod p
        mod = pk1
        a = ap % mod
        for _ in range(kmax + 3):   # Hensel to full modulus
            a = (a - (a * a - ap * a + p) * pow(2 * a - ap, -1, mod)) % mod
        alpha = a
        lamhat = pow(lam0, p ** kmax, pk1)
        T = [1]
        for k in range(1, kmax + 1):
            acc = 0
            po = 1
            for m in range(p ** k):
                acc = (acc + d[m] * po) % pk1
                po = po * lamhat % pk1
            T.append(acc)
        row = {'ap': ap, 'alpha_mod': int(alpha), 'T1_unit': T[1] % p != 0}
        if T[1] % p == 0:
            row['skip'] = 'nonordinary read'
            rows[str(lam0)] = row
            continue
        a5 = pow(alpha, 5, pk1)
        a3 = pow(alpha, 3, pk1)
        for k in range(1, kmax):
            mod = p ** (k + 1)
            u = T[k + 1] * pow(T[k], -1, mod) % mod
            # measured agreement depth vs alpha^5
            diff = (u - a5) % mod
            depth = k + 1 if diff == 0 else (r2.valp(diff, p) or 0)
            row['u_level%d' % k] = int(u)
            row['u_level%d_match_depth_vs_alpha5' % k] = depth
        # Hankel v1
        num = (T[3] * T[1] - T[2] * T[2]) % pk1
        den = (T[2] - T[1] * T[1]) % pk1
        row['hankel_num_val'] = r2.valp(num, p) if num else None
        row['hankel_den_val'] = r2.valp(den, p) if den else None
        if den % p != 0 and num and (r2.valp(num, p) or 0) >= 1:
            H = (num // p) * pow(den, -1, p ** kmax) % (p ** kmax)
            u2 = T[3] * pow(T[2], -1, p ** kmax) % (p ** kmax)
            v1 = H * pow(u2, -1, p ** kmax) % (p ** kmax)
            diff = (v1 - a3) % (p ** kmax)
            depth = kmax if diff == 0 else (r2.valp(diff, p) or 0)
            row['v1_read'] = int(v1)
            row['v1_match_depth_vs_alpha3'] = depth
        rows[str(lam0)] = row
    return rows

# --------------------------------------------------------------- p-curvature

def pcurv_jordan(p, x0, e=-2):
    """Jordan profile of the p-curvature of the L6-class operator at fiber x0
    mod p, by the (x - x0)-jet route (an independent code path): companion A(x) of the
    monic L6 (y^(6) = -sum a_i y^(i)), jets of order p; M_1 = A,
    M_{k+1} = M_k' + M_k A; Jordan type from ranks of powers of M_p(x0)."""
    import numpy as np
    A = r2.load_L6(l6_path())
    # polynomial coefficient rows a_i(x), i=0..6; monic-ize at the fiber jet
    T = p + 2  # jet length
    def poly_jet(row):
        # Taylor shift by repeated synthetic division: jet coeffs of the poly at x0
        work = [c % p for c in row]
        out = []
        for _ in range(T):
            if not work:
                break
            # synthetic division of work by (x - x0), high-to-low
            q = [0] * (len(work) - 1)
            acc = 0
            for i in range(len(work) - 1, 0, -1):
                acc = (acc * x0 + work[i]) % p
                q[i - 1] = acc
            rem = (acc * x0 + work[0]) % p
            out.append(rem)
            work = q
        return (out + [0] * T)[:T]
    jets = [poly_jet(row) for row in A]
    import numpy as _np
    def jmul(a, b):
        # int64-safe: values < p <= 200, length <= p+2: max sum ~ p^2*(p+2) << 2^63
        out = _np.convolve(_np.asarray(a, dtype=_np.int64),
                           _np.asarray(b, dtype=_np.int64))[:T] % p
        return out.tolist() + [0] * (T - min(T, len(a) + len(b) - 1))
    def jder(a):
        return [(k + 1) * a[k + 1] % p for k in range(T - 1)] + [0]
    lead = jets[6]
    if lead[0] % p == 0:
        return {'skip': 'leading coeff vanishes at fiber (pole of presentation)'}
    # invert lead as jet
    inv0 = pow(lead[0], -1, p)
    linv = [inv0] + [0] * (T - 1)
    for _ in range(T.bit_length() + 1):
        e2 = jmul(lead, linv)
        e2[0] = (e2[0] - 2) % p
        linv = [(-x) % p for x in jmul(linv, e2)]
    ai = [jmul(jets[i], linv) for i in range(6)]
    # companion: A[j][j+1] = 1; A[5][i] = -ai[i]
    Amat = [[[0] * T for _ in range(6)] for _ in range(6)]
    one = [1] + [0] * (T - 1)
    for j in range(5):
        Amat[j][j + 1] = one[:]
    for i in range(6):
        Amat[5][i] = [(-c) % p for c in ai[i]]
    M = [row[:] for row in Amat]
    for step in range(p - 1):
        new = [[jder(M[i][j]) for j in range(6)] for i in range(6)]
        for i in range(6):
            for j in range(6):
                acc = new[i][j]
                for k in range(6):
                    if any(M[i][k]) and any(Amat[k][j]):
                        pr = jmul(M[i][k], Amat[k][j])
                        acc = [(a + b) % p for a, b in zip(acc, pr)]
                new[i][j] = acc
        M = new
    import numpy as np
    Mp = np.array([[M[i][j][0] % p for j in range(6)] for i in range(6)], dtype=np.int64)
    def rank_modp(mat):
        m = mat.copy() % p
        rk = 0
        rows, cols = m.shape
        rr = 0
        for c in range(cols):
            piv = None
            for r_ in range(rr, rows):
                if m[r_, c] % p:
                    piv = r_
                    break
            if piv is None:
                continue
            m[[rr, piv]] = m[[piv, rr]]
            inv = pow(int(m[rr, c]), -1, p)
            m[rr] = m[rr] * inv % p
            for r_ in range(rows):
                if r_ != rr and m[r_, c] % p:
                    m[r_] = (m[r_] - m[r_, c] * m[rr]) % p
            rr += 1
            rk += 1
        return rk
    ranks = []
    P_ = Mp.copy()
    for _ in range(6):
        ranks.append(rank_modp(P_))
        P_ = P_ @ Mp % p
    nilpotent = ranks[-1] == 0
    # jordan type from the rank sequence: #blocks of size >= k = r_{k-1} - r_k
    prof = tuple(ranks[:5])
    rseq = [6] + ranks          # rseq[k] = rank(N^k)
    ge = [rseq[k - 1] - rseq[k] for k in range(1, 8 - 1)]  # k = 1..6
    blocks = []
    for k in range(1, 7):
        exact = ge[k - 1] - (ge[k] if k < 6 else 0)
        blocks += [k] * max(0, exact)
    blocks = sorted(blocks, reverse=True)
    return {'p': p, 'x0': x0, 'rank_sequence': ranks, 'nilpotent': bool(nilpotent),
            'profile': str((bool(nilpotent), prof)), 'jordan_type': blocks}

# --------------------------------------------------------------- enumeration

def enumerate_sextics(p, u3, u3_mod, v1=None, v1_depth=0, ap_ell=None):
    """Weil-box + p^4-FE enumeration under declared factor shapes.
    u3: unit root mod p^{u3_mod}; v1: slope-1 unit read mod p^{v1_depth} or None;
    ap_ell: exact elliptic a_p at the fiber (control-2 datum).
    Returns dict with S1/S2 candidate lists (each candidate = dict of
    quadratic traces + the sextic coefficients) and adjudication fields."""
    out = {'p': p, 'u3_mod': u3_mod, 'v1_depth': v1_depth}
    pm = p ** u3_mod
    # outer trace A = alpha0 + p^4/alpha0 mod p^{u3_mod}
    A_mod = (u3 + p ** 4 * pow(u3, -1, pm)) % pm
    Amax = 2 * p * p
    cands = [a for a in range(-Amax, Amax + 1) if (a - A_mod) % pm == 0]
    out['outer_A_candidates'] = cands
    ap2 = None if ap_ell is None else ap_ell * ap_ell - 2 * p
    out['ap_ell'] = ap_ell
    out['ap2_ell'] = ap2
    S1, S2 = [], []
    if ap_ell is not None:
        B1 = p * ap2
        for A in cands:
            for w in (-2, -1, 0, 1, 2):
                t = w * p * p
                # v1 congruence gates (only at measured depth)
                ok = True
                if v1 is not None and v1_depth >= 1:
                    ok = (v1 - ap2) % (p ** min(v1_depth, u3_mod)) == 0
                if ok:
                    S1.append({'A': A, 'B': B1, 't': t, 'w': w})
        a2 = p * ap_ell
        for A in cands:
            for w in (-2, -1, 0, 1, 2):
                t = w * p * p
                ok = True
                if v1 is not None and v1_depth >= 1:
                    # S2 slope-1 root = p * gamma, gamma = unit root of E
                    gam = ap_ell % p  # gamma == a_p mod p
                    ok = (v1 - gam) % p == 0 if v1_depth >= 1 else True
                if ok:
                    S2.append({'a': a2, 't': t, 'w': w})
    out['S1_candidates'] = S1
    out['S2_candidates'] = S2
    return out

def sextic_from_S1(A, B, t, p):
    """(1 - A T + p^4 T^2)(1 - B T + p^4 T^2)(1 - t T + p^4 T^2) -> coeffs c0..c6."""
    P4 = p ** 4
    def q(s):
        return [1, -s, P4]
    def mul(a, b):
        out = [0] * (len(a) + len(b) - 1)
        for i, ai in enumerate(a):
            for j, bj in enumerate(b):
                out[i + j] += ai * bj
        return out
    return mul(mul(q(A), q(B)), q(t))

def weil_circle_ok(coeffs, p, wt=4):
    """All complex roots of reciprocal poly on |z| = p^{wt/2} (numeric screen)."""
    import numpy as np
    c = np.array(coeffs, dtype=float)
    rts = np.roots(c[::-1])
    tgt = p ** (wt / 2.0)
    return bool(np.all(np.abs(np.abs(rts) - tgt) / tgt < 1e-6))
