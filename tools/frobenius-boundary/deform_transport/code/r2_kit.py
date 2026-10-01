#!/usr/bin/env python3
"""r2_kit.py — guard-digit-tracked mod-p^S power-series engine (vendored
single-series engine for the deform_transport member).

A mod-p^S engine for the log-free solutions of an order-6 operator (and
control operators), pushed to N = p^2 + O(1) terms so the level-1 Dwork
ratio at Teichmuller points yields candidate unit roots mod p^2 per
(prime, fiber):

    u_s(p, x0)  =  [G_s]_{<p^2}(xhat) / [G_s]_{<p}(xhat)   mod p^2,
    xhat = Teich(x0) = x0^p mod p^2  (so xhat^p = xhat: rational gauges
    A(x) contribute A(xhat)/A(xhat^p) = 1 exactly; constants cancel in the
    ratio),  G_s = x^{-s} F_s,  with F_s the canonical log-free solutions.

CANDIDATE a_p: centered lift of u mod p^2 (for motivic weight 2 or 3 the
p^w/u antidual term vanishes mod p^2); Weil windows |a| <= 2p (wt-3)
and |a| <= 2 p^{3/2} (wt-4, needs p >= 17) are SCREENS, not claims.

PRECISION LAW: the forward recurrence divides by the indicial value W0(m);
at residues where that value has p-valuation v, each division by p^v costs
v digits of absolute precision FOREVER AFTER (the ambiguity propagates
through a homogeneous solution that later singular steps divide by p
again).  Total loss L = sum of pivot valuations over m in [m_start, N].
The engine allocates S = 2 + L + 8 p-digits, verifies rhs divisibility at
every singular step (an integrality violation is a recorded event, not a
crash), tracks realized loss against the precomputed ledger, and shrinks
the working modulus as remaining loss decreases.  Final values are EXACT
residues mod p^{S-L} >= p^2 — no float, no estimate.

SEEDING: char-0 singular positions are handled by seeding the recurrence
with the EXACT rational char-0 canonical-basis coefficients (Fraction
solver), so the p-adic engine starts where every pivot is p-adically
finite.

CONTROLS: C1 ground truth (engine coefficients vs exact char-0 Fraction
solve, residue-exact); C2 Legendre end-to-end (unit root mod p^2 vs the
EXACT character-sum a_p of the quadratic-twist curve); the raw mod-p reads
double as a cross-check surface.  The Dwork-ratio reading for
maximally-unipotent-adjacent points is EMPIRICAL; the Weil-window
statistics over a prime grid measure whether the hypothesis has content —
empirical certification, not a theorem.
"""
import json, os, sys, time
from fractions import Fraction
try:
    from gmpy2 import mpz, invert as _invert
    def INV(a, m):
        return _invert(a, m)
    BACKEND = "gmpy2"
except Exception:                                    # pragma: no cover
    def mpz(x):
        return x
    def INV(a, m):
        return pow(int(a), -1, int(m))
    BACKEND = "int"

INDICIAL_CONST = 169840080411043248019112229026157035520
HYP2F1 = [[-1], [4, -8], [0, 4, -4]]   # order 2, annihilates 2F1(1/2,1/2;1;lam)
XGRID = [2, 3, 4, 5, 7, 8, 9, 10, 11, 13, (1, 2)]

# ---------------------------------------------------------------- basics

def load_L6(path):
    d = json.load(open(path))
    A = [[int(c) for c in row] for row in d['coeffs']]
    assert len(A) == 7 and all(len(r) == 281 for r in A)
    return A

def pad(A):
    D = max(len(r) for r in A) - 1
    return [list(r) + [0] * (D + 1 - len(r)) for r in A], D

def is_split(p):
    return pow(-5 % p, (p - 1) // 2, p) == 1

def good_prime(p):
    return p >= 23 and is_split(p) and INDICIAL_CONST % p != 0

def centered(r, M):
    r %= M
    return r - M if r > M // 2 else r

def valp(n, p):
    if n == 0:
        return None
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v

# ------------------------------------------------- operator scaffolding

def op_scaffold(A, e):
    """s0, D, r, and W_k(n) as exact-integer callables (k = diagonal index)."""
    A7, D = pad(A)
    r = len(A) - 1
    val_i = [next(j for j, c in enumerate(A7[i]) if c) for i in range(r + 1)]
    s0 = min(val_i[i] - i for i in range(r + 1))
    dmax = D - s0

    def Wk(k, n):
        s = 0
        for i in range(r + 1):
            j = i + s0 + k
            if 0 <= j <= D and A7[i][j]:
                ff = 1
                for t in range(i):
                    ff *= (n + e - t)
                s += A7[i][j] * ff
        return s
    return A7, D, r, s0, dmax, Wk

# ------------------------------------------------- char-0 exact solver

def logfree_char0(A, e, N, order):
    """Log-free solutions over Q to N terms; returns (cols, freepos), cols are
    dicts m->Fraction (param columns)."""
    A7, D, r, s0, dmax, Wk = op_scaffold(A, e)
    cols, freepos = [], []
    for m in range(N + 1):
        rhs = [sum(Wk(k, m - k) * col.get(m - k, Fraction(0))
                   for k in range(1, min(m, dmax) + 1)) for col in cols]
        lead = Wk(0, m)
        if lead != 0:
            for col, rv in zip(cols, rhs):
                col[m] = Fraction(-rv, 1) / lead
        else:
            nz = [j for j, rv in enumerate(rhs) if rv != 0]
            if nz:
                j = nz[-1]
                cj = rhs[j]
                base = cols.pop(j)
                freepos.pop(j)
                rhs2 = rhs[:j] + rhs[j + 1:]
                for k2, col in enumerate(cols):
                    fac = Fraction(-rhs2[k2], 1) / cj
                    for mm, vv in base.items():
                        col[mm] = col.get(mm, Fraction(0)) + fac * vv
            cols.append({m: Fraction(1)})
            freepos.append(m)
    return cols, freepos

def canonical_seed(A, e, order, freepos_want, upto):
    """Exact rational canonical-basis coefficients c_0..c_upto:
    basis[k] = solution with 1 at freepos_want[k], 0 at the others."""
    cols, freepos = logfree_char0(A, e, upto, order)
    assert freepos == freepos_want, (freepos, freepos_want)
    P = len(cols)
    M = [[cols[b].get(fp, Fraction(0)) for b in range(P)] for fp in freepos_want]
    if P == 1:
        inv = [[1 / M[0][0]]]
    elif P == 2:
        det = M[0][0] * M[1][1] - M[0][1] * M[1][0]
        assert det != 0
        inv = [[M[1][1] / det, -M[0][1] / det], [-M[1][0] / det, M[0][0] / det]]
    else:
        raise NotImplementedError
    basis = []
    for k in range(P):
        col = {}
        for m in range(upto + 1):
            v = sum(cols[b].get(m, Fraction(0)) * inv[b][k] for b in range(P))
            if v != 0:
                col[m] = v
        basis.append(col)
    return basis

# ------------------------------------------------- loss ledger

def loss_ledger(A, e, p, N, m_start):
    """Sum of p-valuations of W0(m) over m in [m_start, N] (char-0 zeros of
    W0 must lie below m_start). Returns (L_total, suffix) with suffix[i] =
    loss at steps > singular_step_i (for modulus shrinking)."""
    A7, D, r, s0, dmax, Wk = op_scaffold(A, e)
    resid = [m for m in range(p) if Wk(0, m) % p == 0]
    sing = []
    for r0 in resid:
        m = r0 if r0 >= m_start else r0 + p * ((m_start - r0 + p - 1) // p)
        while m <= N:
            w = Wk(0, m)
            assert w != 0, (m, 'char-0 zero above m_start')
            sing.append((m, valp(w, p)))
            m += p
    sing.sort()
    L = sum(v for _, v in sing)
    return L, sing

# ------------------------------------------------- the mod-p^S engine

def unit_root_lift(A, e, p, seeds, seed_upto, xgrid, lead_exps,
                   extra_terms=None, progress=None):
    """Guard-tracked mod-p^S solve to N = p^2 - 1 + max(lead_exps) (so every
    column's [G]_{<p^2} is complete), Dwork ratios at Teichmuller points.

    seeds: list (per column) of dict m->Fraction, exact char-0 coefficients
           through m = seed_upto. Engine recurrence starts at seed_upto+1.
    lead_exps: leading exponent s per column (F_s = x^s + ...).
    Returns receipt dict."""
    t0 = time.time()
    A7, D, r, s0, dmax, Wk = op_scaffold(A, e)
    ncol = len(seeds)
    smax = max(lead_exps)
    N = p * p - 1 + (smax if extra_terms is None else extra_terms)
    m_start = seed_upto + 1
    Ltot, sing = loss_ledger(A, e, p, N, m_start)
    S = 2 + Ltot + 8
    p2 = p * p
    rec = {"p": p, "N": N, "S_alloc": S, "L_ledger": Ltot,
           "n_singular_steps": len(sing), "e": e, "backend": BACKEND,
           "seed_upto": seed_upto, "lead_exps": lead_exps}
    # seed p-integrality (fixed rationals; denominators must be prime to p)
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
    # working modulus + suffix-loss table for shrinking
    suffix_after = {}   # singular index -> loss remaining AFTER that step
    run = 0
    for i in range(len(sing) - 1, -1, -1):
        suffix_after[sing[i][0]] = run
        run += sing[i][1]
    M = mpz(p ** S)
    # FD registers: reg[k][0..r] forward differences of W_k at n = m-k
    # init at m0 = m_start: reg[k] built from W_k(m0-k .. m0-k+r)  (deg <= r)
    def fd_init(m0, mod):
        # registers at argument n = m0 - k; each advance moves n -> n+1.
        # trailing zero differences trimmed (W_k has degree <= r in n).
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
    # window list indexed by m; old entries nulled once out of the recurrence window
    window = []
    for m in range(seed_upto + 1):
        rowv = []
        for c in range(ncol):
            fr = seedvals[c][m]
            rowv.append(mpz(fr.numerator) * INV(fr.denominator, M) % M)
        window.append(rowv)
    # accumulators (mod p^2); fibers are ints or (num, den) rationals.
    # xhat = Teichmuller lift mod p^2 (multiplicative: Teich(n/d)=Teich(n)/Teich(d))
    def fib_parts(x0):
        return (x0, 1) if isinstance(x0, int) else x0
    def fib_key(x0):
        n, d = fib_parts(x0)
        return str(n) if d == 1 else "%d/%d" % (n, d)
    xhat = {}
    xres = {}   # fiber residue mod p
    for x0 in xgrid:
        n, d = fib_parts(x0)
        xhat[x0] = pow(n, p, p2) * pow(pow(d, p, p2), -1, p2) % p2
        xres[x0] = n * pow(d, -1, p) % p
    accs = {x0: [0] * ncol for x0 in xgrid}       # sum c_m xhat^m mod p^2
    snapB = {x0: [None] * ncol for x0 in xgrid}   # at m = s+p-1
    snapT = {x0: [None] * ncol for x0 in xgrid}   # at m = s+p^2-1
    xpow = {x0: 1 for x0 in xgrid}
    raw_modp = {x0: [0] * ncol for x0 in xgrid}   # sum_{m<p} c_m x0^m mod p (raw read)
    x0pow = {x0: 1 for x0 in xgrid}
    for m in range(seed_upto + 1):
        for x0 in xgrid:
            for c in range(ncol):
                accs[x0][c] = (accs[x0][c] + window[m][c] * xpow[x0]) % p2
                if m < p:
                    raw_modp[x0][c] = (raw_modp[x0][c] + window[m][c] * x0pow[x0]) % p
        # snapshots cannot occur in the seed range (Bpoint >= p-1 > seed_upto)
        for x0 in xgrid:
            xpow[x0] = xpow[x0] * xhat[x0] % p2
            if m < p:
                x0pow[x0] = x0pow[x0] * xres[x0] % p
    # snapshot checkpoints
    Bpoint = {c: lead_exps[c] + p - 1 for c in range(ncol)}
    Tpoint = {c: lead_exps[c] + p2 - 1 for c in range(ncol)}
    loss = 0
    events = []
    alive = [True] * ncol
    last_shrink_loss = 0
    check_every = 4096
    for m in range(m_start, N + 1):
        kmax = min(m, dmax)
        # advance FD registers to argument n = m-k (they sit at m-1-k after last step)
        for reg in regs:
            for lev in range(len(reg) - 1):
                reg[lev] = (reg[lev] + reg[lev + 1]) % M
        # guardrail: direct recompute of W_0(m) occasionally
        if m % check_every == 0:
            w0d = Wk(0, m) % M
            assert regs[0][0] == w0d, (m, 'FD drift')
        pivot = regs[0][0]
        rhs = [0] * ncol
        for k in range(1, kmax + 1):
            wkv = regs[k][0]
            if wkv:
                row = window[m - k]
                for c in range(ncol):
                    if alive[c] and row[c] is not None:
                        rhs[c] += wkv * row[c]
        piv = pivot % M
        v = 0
        pw = piv
        while pw % p == 0 and v < 6:
            pw //= p
            v += 1
        if pw % p == 0:
            events.append([m, "pivot_val_gt6"])
            rec["events"] = events
            rec["ABORT"] = "pivot valuation > 6"
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
                    # column's char-0 solution is NOT p-integral: the naive
                    # truncation congruence is ill-defined for it — recorded
                    # event, column retired, survivors continue.
                    events.append({"m": m, "col_lead_exp": lead_exps[c],
                                   "event": "integrality_violation",
                                   "pivot_val": v, "rhs_val": valp(rc, p)})
                    alive[c] = False
                    continue
                newrow[c] = (rc // pv) * inv % M
            loss += v
            # shrink working modulus when worthwhile
            if loss - last_shrink_loss >= 48:
                Srem = 2 + suffix_after.get(m, 0) + 8
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
        # accumulate values
        for x0 in xgrid:
            xp = xpow[x0]
            av = accs[x0]
            for c in range(ncol):
                if newrow[c]:
                    av[c] = int(av[c] + (newrow[c] % p2) * xp) % p2
            xpow[x0] = xp * xhat[x0] % p2
            if m < p:
                for c in range(ncol):
                    if newrow[c] is not None:
                        raw_modp[x0][c] = int(raw_modp[x0][c]
                                              + newrow[c] * x0pow[x0]) % p
                x0pow[x0] = x0pow[x0] * xres[x0] % p
        for c in range(ncol):
            if alive[c] and m == Bpoint[c]:
                for x0 in xgrid:
                    snapB[x0][c] = accs[x0][c]
            if alive[c] and m == Tpoint[c]:
                for x0 in xgrid:
                    snapT[x0][c] = accs[x0][c]
        # drop old window entries (keep last dmax)
        if len(window) > dmax + 1 and m > m_start + dmax + 2:
            window[m - dmax - 1] = None
        if progress and m % progress == 0:
            print("  m=%d/%d loss=%d/%d t=%.1fs" % (m, N, loss, Ltot,
                                                    time.time() - t0), flush=True)
    prec_final = S - loss
    rec["alive_at_end"] = {("s%d" % lead_exps[c]): alive[c]
                           for c in range(ncol)}
    if events:
        rec["events"] = events
    rec["loss_realized"] = loss
    rec["final_precision_digits"] = prec_final
    rec["precision_ok"] = prec_final >= 2
    assert loss == Ltot, (loss, Ltot, "ledger mismatch")
    # Dwork ratios
    fibers = {}
    for x0 in xgrid:
        fx = {"xhat_mod_p2": xhat[x0]}
        xr = xres[x0]
        # singular-fiber screen mod p (on the residue)
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
            xs_inv = pow(pow(xhat[x0], s, p2), -1, p2)
            if snapB[x0][c] is None or snapT[x0][c] is None:
                cols["s%d" % s] = {"non_integral": True}
                continue
            B = int(snapB[x0][c]) * xs_inv % p2   # [G]_{<p}(xhat) mod p^2
            T = int(snapT[x0][c]) * xs_inv % p2   # [G]_{<p^2}(xhat) mod p^2
            d = {"B": int(B), "T": int(T)}
            if B % p == 0:
                d["ordinary"] = False
            else:
                d["ordinary"] = True
                u = T * pow(B, -1, p2) % p2
                d["u_mod_p2"] = int(u)
                d["u_mod_p"] = int(u % p)
                a = int(centered(u, p2))
                d["a_cand"] = a
                d["weil_wt3"] = abs(a) <= 2 * p
                d["weil_wt4"] = abs(a) * abs(a) <= 4 * p**3
            cols["s%d" % s] = d
        fx["cols"] = cols
        fibers[fib_key(x0)] = fx
    rec["fibers"] = fibers
    rec["t_wall_s"] = round(time.time() - t0, 2)
    return rec

# ------------------------------------------------- Legendre exact a_p

def curve_ap(p, x0):
    """a_p of C_x0: v^2 = u(x0-u)(1-x0*u), exact character sum (pure python)."""
    ap = 0
    for u in range(1, p):
        f = u * ((x0 - u) % p) % p * ((1 - x0 * u) % p) % p
        if f == 0:
            continue
        ap -= 1 if pow(f, (p - 1) // 2, p) == 1 else -1
    return ap

def legendre_control(p, xs=(2, 3, 4, 5)):
    """Engine end-to-end on 2F1: candidate a_p = chi2(-1)*centered(u + p/u)
    vs exact character-sum a_p. Returns receipt with per-fiber PASS/FAIL."""
    seeds = [{0: Fraction(1)}]
    rec = unit_root_lift(HYP2F1, 0, p, seeds, 0, [x * x % p for x in xs],
                         [0])
    if "fibers" not in rec:
        return rec
    chim1 = 1 if p % 4 == 1 else -1
    p2 = p * p
    out = {"p": p, "engine": {k: rec[k] for k in
                              ("S_alloc", "L_ledger", "loss_realized",
                               "final_precision_digits", "t_wall_s")},
           "fibers": {}}
    allpass = True
    for x in xs:
        lam = x * x % p
        fx = rec["fibers"].get(str(lam))
        row = {"lam": lam}
        if fx is None or "skip" in fx:
            row["skip"] = "singular"
        elif lam in (0, 1):
            row["skip"] = "degenerate lam"
        else:
            d = fx["cols"]["s0"]
            apx = curve_ap(p, x)
            row["ap_exact"] = apx
            if not d["ordinary"]:
                row["ordinary"] = False
                row["ap_div_p"] = (apx % p == 0)
                allpass = allpass and row["ap_div_p"]
            else:
                u = d["u_mod_p2"]
                acand = chim1 * centered((u + p * pow(u, -1, p2)) % p2, p2)
                # note: u + p*u^{-1} mod p^2 (weight-1 motive: antidual term p/u)
                row["a_cand"] = acand
                row["match"] = (acand == apx)
                allpass = allpass and row["match"]
        out["fibers"][str(x)] = row
    out["ALL_PASS"] = allpass
    return out

# ------------------------------------------------- L6 wrappers

def l6_seeds(A):
    """Canonical char-0 seeds for the L6-class operator at e=-2 (freepos {2,7}), to m=8."""
    return canonical_seed(A, -2, 6, [2, 7], 8)

def l6_lift(A, p, seeds, progress=None):
    if not good_prime(p):
        return {"p": p, "skip": "not a good split prime"}
    return unit_root_lift(A, -2, p, seeds, 8, XGRID, [2, 7],
                          progress=progress)

# ------------------------------------------------- ground truth gate

def ground_truth_gate(A, p=23, upto=70):
    """Engine window coefficients vs exact char-0 Fractions, residue-exact."""
    seeds = l6_seeds(A)
    basis = canonical_seed(A, -2, 6, [2, 7], upto)
    # run a tiny engine to N=upto by reusing unit_root_lift internals is
    # overkill; instead run the real lift at p (N=p^2+6 >= upto for p>=9)
    # and have it export its early window — simpler: rerun logic inline.
    A7, D, r, s0, dmax, Wk = op_scaffold(A, -2)
    Ltot, sing = loss_ledger(A, -2, p, upto, 9)
    S = 2 + Ltot + 8
    M = mpz(p ** S)
    window = []
    for m in range(9):
        window.append([basis[c].get(m, Fraction(0)) for c in range(2)])
    win = [[fr.numerator * pow(fr.denominator, -1, M) % M for fr in row]
           for row in window]
    loss = 0
    alive = [True, True]
    report = {"p": p, "upto": upto, "S": S, "L": Ltot, "mismatches": [],
              "col_events": []}
    for m in range(9, upto + 1):
        rhs = [0, 0]
        for k in range(1, min(m, dmax) + 1):
            w = Wk(k, m - k) % M
            if w:
                for c in range(2):
                    if alive[c] and win[m - k][c] is not None:
                        rhs[c] += w * win[m - k][c]
        v = valp(Wk(0, m), p) or 0
        pv = p ** v
        unit = (Wk(0, m) // pv) % M
        inv = pow(unit, -1, M)
        row = [None, None]
        for c in range(2):
            if not alive[c]:
                continue
            rc = (-rhs[c]) % M
            if rc % pv != 0:
                # cross-check the drop against exact valuations
                fr = basis[c].get(m, Fraction(0))
                report["col_events"].append(
                    [m, c, "engine_drop", "exact_den_val_p",
                     valp(fr.denominator, p)])
                alive[c] = False
                continue
            row[c] = (rc // pv) * inv % M
        loss += v
        win.append(row)
        prec = S - loss
        Mp = p ** prec
        for c in range(2):
            if not alive[c] or row[c] is None:
                continue
            fr = basis[c].get(m, Fraction(0))
            if fr.denominator % p == 0:
                report["mismatches"].append([m, c, "exact_nonintegral_but_engine_alive"])
                continue
            exact = fr.numerator * pow(fr.denominator, -1, Mp) % Mp
            if exact != row[c] % Mp:
                report["mismatches"].append([m, c])
    report["loss"] = loss
    report["alive"] = alive
    # PASS: no value mismatches, and every engine drop coincides with a
    # genuine exact-denominator p (col_events all carry exact_den_val >= 1)
    drops_ok = all(ev[4] and ev[4] >= 1 for ev in report["col_events"])
    report["drops_match_exact"] = drops_ok
    report["PASS"] = (not report["mismatches"]) and drops_ok
    return report

if __name__ == "__main__":
    print("r2_kit: library module (imported by lfact_kit / deform_kit)")
