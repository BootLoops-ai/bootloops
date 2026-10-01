#!/usr/bin/env python3
"""mirror6 — exact mirror-map / instanton-type fingerprint engine for MUM-type
points of D-finite operators.

WHAT IT COMPUTES, given an operator with a MUM-type point (a point whose local
monodromy has a unipotent Jordan block of size k >= 2 on an integer exponent
class):
  1. The exact log-Frobenius CHAIN y_0 .. y_{k-1} of the maximal block
     (parameterized tower solve; every resonance handled by exact elimination;
     leftover free parameters frozen to 0 = the CANONICAL-ZERO GAUGE,
     DEFINITION D3 below; the freeze list is part of the receipt).
  2. The canonical q-coordinate (mirror map), DEFINITION D1:
         log q = y_1 / y_0   (normalized so q = s(1 + O(s)))
     and its inverse series s(q) exactly over Q.
  3. The structure-series ladder (couplings), DEFINITION D2: with
     R_j = y_j/y_0 rewritten in t = log q,
         alpha_j = theta(...theta(theta(R_j)/alpha_1)...)/alpha_{j-1})
     (j-fold theta = q d/dq, dividing by the previously found alpha's);
     alpha_1 = 1 by construction; alpha_2 = the Yukawa-type coupling
     normalized alpha_2 = 1 + O(q).  For an order-4 CY operator
     n0 * alpha_2 is the standard normalized Yukawa coupling K(q)
     (van Straten 1704.00164: K(q) = n0 + sum n_d d^3 q^d/(1-q^d)).
  4. Instanton-type numbers, DEFINITION D4: Lambert-kernel inversion of
     alpha_2 - 1 (and higher alphas) at kernel weights w in {2,3,4}:
         alpha - 1 = sum_d m_d d^w q^d/(1-q^d),  m_d in Q exactly.
     The CY3 convention is w=3 (the control's published convention).
  5. Script-emitted integrality verdicts (per series, per kernel):
     "integral" / "integral-after-rescaling-by-N=<lcm of denominators>"
     / "non-integral: <failure pattern>", plus the q-rescale
     (kappa) integrality mode for the mirror map.
  6. Arithmetic fingerprint: p-adic valuation tables of the integralized
     numbers, growth/ratio-test radius (conifold-distance read, numeric
     measured-grade), and a rational-q-point value screen (data only, no
     constant-recognition verdicts — sealed scans belong to a separate
     PSLQ gate layer).

DEFINITIONS (the WARNING-LAW register; every generalization is pinned):
  D1 q-coordinate: q := s * exp(h/y_0) where y_1 = y_0 log s + h is the
     chain log-partner; the s^rho coefficient of h is gauged to 0 so that
     q = s(1+O(s)).  At a TRUE MUM (single block, one exponent) this is the
     standard mirror map.  At a partial-MUM point ((5,1) or (3,1,1,1)) the
     chain members are ambiguous by log-free solutions whose admixture
     parameters the solve either FIXES (obstruction-cured) or freezes to 0
     (canonical-zero gauge); the receipt lists both sets.  The FIXED ones are
     forced by the operator; only the frozen ones are conventions.
  D2 couplings: the iterated theta-quotient ladder above — for order-4 CY
     operators it reproduces the classical normal form
     L ~ theta^2 (1/K) theta^2 (self-duality check alpha_3 == alpha_2 is a
     control row); for a general order-6 weight-4 shape it is taken as the
     DEFINITION of the coupling tower.
  D3 canonical-zero gauge: all resonance parameters still free after the
     chain solve are set to 0.  A sensitivity block re-runs the fingerprint
     with each frozen parameter set to 1 to measure which outputs are
     gauge-robust.
  D4 kernel weights: w=3 is the CY3 (genus-0 GV / Lambert d^3) convention,
     applied verbatim to the weight-4 wall as a DEFINITION; w=2 and w=4
     columns are emitted alongside so the weight-4 reading is not presumed.
  D5 normalization n0: the ladder normalizes alpha_2(0)=1 ("n0=1 frame").
     Multiplying by an external n0 (published H^3 for controls) is an INPUT,
     never derived here; downstream receipts carry the n0=1 numbers plus
     the lcm-rescale N.

No unit-root, CY, or geometric-mirror claim is made for any input operator
anywhere in this module; the outputs are exact invariants of the operator's
local structure under the definitions above.  No Pi/Hodge-functional claim
(the chain data is local solution structure, not a period functional).

Exact arithmetic: gmpy2.mpq end-to-end; numerics only in the labeled
measured-grade screens (mpmath, dps stated).
"""
import hashlib
import json
import math
import os
import sys
import time

from gmpy2 import mpq, mpz

MPQ0 = mpq(0)
MPQ1 = mpq(1)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for blk in iter(lambda: f.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


# ---------------------------------------------------------------- polynomials
# poly = list of mpq, low -> high.

def peval(poly, x):
    acc = MPQ0
    for c in reversed(poly):
        acc = acc * x + c
    return acc


def pderiv(poly):
    return [mpq(i) * poly[i] for i in range(1, len(poly))]


def padd(a, b):
    n = max(len(a), len(b))
    out = [MPQ0] * n
    for i, c in enumerate(a):
        out[i] += c
    for i, c in enumerate(b):
        out[i] += c
    return out


def pmul(a, b):
    if not a or not b:
        return []
    out = [MPQ0] * (len(a) + len(b) - 1)
    for i, ca in enumerate(a):
        if ca:
            for j, cb in enumerate(b):
                if cb:
                    out[i + j] += ca * cb
    return out


def pscal(k, a):
    return [k * c for c in a]


def taylor_shift(poly, c):
    """p(x) -> p(x + c) exactly (Horner-style synthetic shift)."""
    out = list(poly)
    n = len(out)
    for i in range(n - 1):
        for j in range(n - 2, i - 1, -1):
            out[j] = out[j] + c * out[j + 1]
    return out


# ------------------------------------------------- operator -> theta form at s
# Operator given as L = sum_i p_i(x) D^i (D = d/dx), p_i coeff lists (mpq).
# At point x0 (s = x - x0): compute the LEFT-NORMAL theta form
#     s^order * L = sum_m s^m R_m(theta_s),   theta_s = s d/ds.
# Using s^i D^i = theta(theta-1)...(theta-i+1) =: FF_i(theta) and
# s^m FF(theta) already left-normal.

def falling_factorial_poly(i):
    """FF_i(u) = u(u-1)...(u-i+1) as coeff list."""
    poly = [MPQ1]
    for r in range(i):
        poly = padd(pmul(poly, [mpq(-r), MPQ1]), [])
    return poly


def _normalize_theta_form(R):
    """Drop zero blocks; left-factor the common s-power (s^v is a unit on
    Laurent-type local solutions, so the solution space is unchanged)."""
    R = {m: poly for m, poly in R.items() if any(poly)}
    if not R:
        raise ValueError("zero operator")
    v = min(R.keys())
    return {m - v: poly for m, poly in R.items()}


def theta_form_at_point(pcoeffs, x0):
    """pcoeffs: list over i of coeff lists p_i(x).  Returns dict m -> R_m(u)
    (theta-polynomial coeff lists) for the left-normal local form
    sum_m s^m R_m(theta), s = x - x0, normalized so R_0 != 0 (the indicial
    block)."""
    order = len(pcoeffs) - 1
    R = {}
    for i, p in enumerate(pcoeffs):
        if not any(p):
            continue
        psh = taylor_shift([mpq(c) for c in p], mpq(x0))  # p_i(s + x0)
        ff = falling_factorial_poly(i)
        # term: p_i(s) * s^{order-i} * FF_i(theta)
        for m_rel, c in enumerate(psh):
            if c == 0:
                continue
            m = m_rel + (order - i)
            R[m] = padd(R.get(m, []), pscal(c, ff))
    return _normalize_theta_form(R)


def theta_form_from_thetaop(zblocks):
    """Operator given directly as sum_m z^m Q_m(theta): zblocks = dict
    m -> theta-poly coeff list (ints/mpq).  Returns normalized dict."""
    return {m: [mpq(c) for c in poly] for m, poly in zblocks.items()}


def theta_form_at_infinity(R):
    """Given left-normal sum_m s^m R_m(theta_s), substitute s = 1/u:
    theta_s = -theta_u; s^m = u^{-m}; multiply left by u^{Mmax}:
        sum_m u^{Mmax-m} R_m(-theta_u).
    Returns the new dict (left-normal in u automatically)."""
    Mmax = max(R.keys())
    out = {}
    for m, poly in R.items():
        neg = [c if i % 2 == 0 else -c for i, c in enumerate(poly)]
        out.setdefault(Mmax - m, [])
        out[Mmax - m] = padd(out[Mmax - m], neg)
    return _normalize_theta_form(out)


# ------------------------------------------------------- parameterized scalars
# A PVal is a linear form c0 + sum_i c_i * P_i over active parameters:
# dict {0: const, param_id: coeff}.  Param ids start at 1.

def pv(c=MPQ0):
    return {0: mpq(c)} if c != 0 else {}


def pv_add(a, b):
    out = dict(a)
    for k, v in b.items():
        nv = out.get(k, MPQ0) + v
        if nv:
            out[k] = nv
        else:
            out.pop(k, None)
    return out


def pv_scal(k, a):
    if k == 0:
        return {}
    return {i: k * v for i, v in a.items()}


def pv_is_zero(a):
    return not a


def pv_const_part(a):
    return a.get(0, MPQ0)


def pv_subst(a, subs):
    """Apply substitution dict param_id -> PVal (repeatedly resolved)."""
    out = {}
    for k, v in a.items():
        if k != 0 and k in subs:
            for k2, v2 in subs[k].items():
                nv = out.get(k2, MPQ0) + v * v2
                if nv:
                    out[k2] = nv
                else:
                    out.pop(k2, None)
        else:
            nv = out.get(k, MPQ0) + v
            if nv:
                out[k] = nv
            else:
                out.pop(k, None)
    return out


# ------------------------------------------------------------ the chain solve

class ChainSolveError(Exception):
    pass


def solve_chain(R, rho, k, N, label="", overrides=None):
    """Exact parameterized log-tower Frobenius solve.

    Ansatz  y_top = sum_{j=0}^{k-1} u_j(s) log(s)^j / j!,
            u_j = s^rho sum_{n=0}^{N} a[j][n] s^n,
    rho = the MINIMAL exponent of the integer class (partners may sit at any
    class exponent; leading zeros are free).  The system is HOMOGENEOUS: every
    resonant order introduces free parameters, every resonant row constrains
    them (exact elimination).  The solution space of the ansatz = all
    solutions of log-degree <= k-1.

    Chain extraction: y_m = sum_{r=0}^{m} u_{(k-1)-m+r} log^r/r!  (m=0..k-1),
    y_0 = u_{k-1} = the chain bottom.  Normalization post-pass: the lowest
    order of u_{k-1} whose param-form is nonzero pins the chain (that form is
    set to 1); all remaining free parameters are frozen to 0 (canonical-zero
    gauge, D3).  If u_{k-1} vanishes identically, there is NO length-k chain
    at this point: ChainSolveError with mechanism (a valid structural
    landing, per the warning law).

    Equations: for each power N', tower row j (top j first):
       sum_{m} sum_{r>=0} R_m^{(r)}(rho+N'-m)/r! * a[j+r][N'-m] = 0.
    R_0(rho+N') multiplies a[j][N']; at resonances (R_0(rho+N') = 0) the row
    is a CONSTRAINT (eliminate a parameter) and a[j][N'] becomes fresh free.
    """
    derivs = {}
    for m, poly in R.items():
        ds = [poly]
        for _ in range(k):
            ds.append(pderiv(ds[-1]))
        derivs[m] = ds
    R0 = R.get(0)
    if R0 is None:
        raise ChainSolveError("no R_0 block (bad theta form)")

    rho = mpq(rho)
    a = [[{} for _ in range(N + 1)] for _ in range(k)]
    subs = {}          # param -> PVal substitution (fully resolved forms)
    next_param = [1]
    introduced = []    # (param_id, j, n)
    eliminated = []    # (param_id, j, n, at_row)
    resonances = []

    # factorial denominators
    invfact = [MPQ1]
    for r in range(1, k + 1):
        invfact.append(invfact[-1] / r)

    def np_table(Np):
        """Precompute [(m, n, r, cval)] for this order (shared across rows)."""
        tab = []
        for m, ds in derivs.items():
            n = Np - m
            if n < 0 or n > N:
                continue
            base = rho + n
            for r in range(0, k):
                poly = ds[r]
                if not poly:
                    continue
                cval = peval(poly, base) * invfact[r]
                if cval:
                    tab.append((m, n, r, cval))
        return tab

    def row_value(j, tab):
        """The (j, N') equation's KNOWN part: everything except the
        R_0(rho+N') * a[j][N'] pivot term."""
        acc = {}
        for (m, n, r, cval) in tab:
            if r == 0 and m == 0:
                continue  # pivot
            jr = j + r
            if jr >= k:
                continue
            term = a[jr][n]
            if term:
                acc = pv_add(acc, pv_scal(cval, term))
        return acc

    for Np in range(0, N + 1):
        tab_np = np_table(Np)
        pivot = peval(R0, rho + Np)
        res = (pivot == 0)
        if res:
            resonances.append(Np)
        for j in range(k - 1, -1, -1):
            known = row_value(j, tab_np)
            known = pv_subst(known, subs)
            if not res:
                a[j][Np] = pv_scal(MPQ1 / pivot, pv_scal(mpq(-1), known))
            else:
                # constraint: known must vanish (homogeneous system: forms
                # carry no constant part, so a nonzero form has parameters)
                if not pv_is_zero(known):
                    pids = [p for p in known if p != 0]
                    if not pids:
                        raise ChainSolveError(
                            f"[{label}] inhomogeneous obstruction at relative "
                            f"order {Np}, row j={j}: {known.get(0)} != 0 "
                            f"(should be impossible in homogeneous solve)")
                    p_el = max(pids)
                    coeff = known[p_el]
                    rest = {kk: vv for kk, vv in known.items() if kk != p_el}
                    sol = pv_scal(mpq(-1) / coeff, rest)
                    sol = pv_subst(sol, subs)
                    subs[p_el] = sol
                    for pk in list(subs):
                        subs[pk] = pv_subst(subs[pk], subs)
                    eliminated.append((p_el, j, Np))
                # fresh free parameter for a[j][Np]
                pid = next_param[0]
                next_param[0] += 1
                a[j][Np] = {pid: MPQ1}
                introduced.append((pid, j, Np))

    # resolve all coefficient forms against the final substitutions
    for j in range(k):
        for n in range(N + 1):
            a[j][n] = pv_subst(a[j][n], subs)

    # chain normalization: lowest order of u_{k-1} with a nonzero form
    pivot_n, pivot_form = None, None
    for n in range(N + 1):
        if a[k - 1][n]:
            pivot_n, pivot_form = n, a[k - 1][n]
            break
    if pivot_form is None:
        raise ChainSolveError(
            f"[{label}] NO length-{k} chain: the log^{k-1} tower row "
            f"vanishes identically after elimination (structural landing)")
    p_star = min(p for p in pivot_form if p != 0)
    c_star = pivot_form[p_star]
    rest = {kk: vv for kk, vv in pivot_form.items() if kk != p_star}
    # p_star = (1 - rest)/c_star
    norm = pv_add(pv(MPQ1), pv_scal(mpq(-1), rest))
    subs2 = {p_star: pv_scal(MPQ1 / c_star, norm)}
    for j in range(k):
        for n in range(N + 1):
            a[j][n] = pv_subst(a[j][n], subs2)
    b = pivot_n  # chain-bottom relative valuation

    # ---- canonical q-GAUGE pass (pole-kill): use the remaining free
    # parameters (log-free admixtures into the chain representatives) to
    # kill every u_j (j <= k-2) coefficient at relative orders < b, plus the
    # u_{k-2}[b] constant (the q-normalization).  Equations processed
    # (n ascending, j descending); an equation reducible to a nonzero
    # CONSTANT is an UN-GAUGEABLE component — recorded; a genuine pole left
    # below the chain bottom is a structural q-gauge obstruction.
    subs3 = {}
    gauge_killed = []
    gauge_unkillable = []
    eqs = []
    for n in range(0, b + 1):
        for j in range(k - 2, -1, -1):
            if n == b and j != k - 2:
                continue
            eqs.append((j, n))
    for (j, n) in eqs:
        form = pv_subst(a[j][n], subs3)
        if pv_is_zero(form):
            continue
        pids = [p for p in form if p != 0]
        if not pids:
            gauge_unkillable.append((int(j), int(n), str(form.get(0))))
            continue
        p_el = max(pids)
        coeff = form[p_el]
        restf = {kk: vv for kk, vv in form.items() if kk != p_el}
        sol = pv_subst(pv_scal(mpq(-1) / coeff, restf), subs3)
        subs3[p_el] = sol
        for pk in list(subs3):
            subs3[pk] = pv_subst(subs3[pk], subs3)
        gauge_killed.append((int(p_el), int(j), int(n)))
    pole_rows = [(j, n, r) for (j, n, r) in gauge_unkillable if n < b]
    if pole_rows:
        raise ChainSolveError(
            f"[{label}] NO CONSISTENT q-GAUGE (structural): log-row "
            f"components below the chain bottom (rel val {b}) cannot be "
            f"gauged away by the log-free space; un-killable (row j, rel "
            f"order n, residual): {pole_rows}.  Mechanism: the chain's "
            f"log-partners carry genuine s^(rho+n) poles relative to y0 — "
            f"the Bessel-class obstruction; a mirror-type q-coordinate "
            f"does not exist at this point under D1.")

    frozen = []
    live = set()
    ov = overrides or {}
    for j in range(k):
        for n in range(N + 1):
            val = pv_subst(a[j][n], subs3)
            a[j][n] = val
            for p in val:
                if p != 0:
                    live.add(p)
    for (pid, j, n) in introduced:
        if pid in live:
            frozen.append((pid, j, n))
    for j in range(k):
        for n in range(N + 1):
            val = a[j][n]
            out = pv_const_part(val)
            for p, coeff in val.items():
                if p != 0 and p in ov:
                    out += coeff * mpq(ov[p])
            a[j][n] = out

    # rebase by b: drop the (now all-zero below b) leading orders so the
    # chain bottom is monic at index 0; effective rho = rho + b
    if b > 0:
        for j in range(k):
            a[j] = a[j][b:] + [MPQ0] * b

    return {
        "a": a,
        "resonances": resonances,
        "chain_bottom_relative_valuation": int(pivot_n),
        "chain_normalization_param": int(p_star),
        "gauge_killed": gauge_killed,
        "gauge_unkillable_constants": gauge_unkillable,
        "params_introduced": [(int(p), int(j), int(n)) for p, j, n in introduced],
        "params_eliminated": [(int(p), int(j), int(n)) for p, j, n in eliminated],
        "params_frozen_to_zero": [(int(p), int(j), int(n)) for p, j, n in frozen],
    }


# ------------------------------------------------------------- series algebra
# series = list of mpq, index = power of q (or s), truncated at length N+1.

def smul(a, b, N):
    out = [MPQ0] * (N + 1)
    for i, ca in enumerate(a[:N + 1]):
        if ca:
            top = min(len(b), N + 1 - i)
            for j in range(top):
                if b[j]:
                    out[i + j] += ca * b[j]
    return out


def sinv(a, N):
    """1/a for a[0] != 0."""
    if a[0] == 0:
        raise ZeroDivisionError("series not a unit")
    out = [MPQ0] * (N + 1)
    out[0] = MPQ1 / a[0]
    for n in range(1, N + 1):
        acc = MPQ0
        for i in range(1, min(n, len(a) - 1) + 1):
            if i < len(a) and a[i]:
                acc += a[i] * out[n - i]
        out[n] = -acc / a[0]
    return out


def sdiv(a, b, N):
    return smul(a, sinv(b, N), N)


def sexp(a, N):
    """exp of series with a[0] = 0 via ODE: E' = a' E."""
    if a[0] != 0:
        raise ValueError("sexp needs valuation >= 1")
    out = [MPQ0] * (N + 1)
    out[0] = MPQ1
    for n in range(1, N + 1):
        acc = MPQ0
        for i in range(1, n + 1):
            if i < len(a) and a[i]:
                acc += mpq(i) * a[i] * out[n - i]
        out[n] = acc / n
    return out


def slog(a, N):
    """log of series with a[0] = 1: L' = a'/a."""
    if a[0] != 1:
        raise ValueError("slog needs a[0]=1")
    ai = sinv(a, N)
    da = [mpq(i) * a[i] for i in range(1, min(len(a), N + 1))]
    integ = smul(da, ai, N)  # (a'/a) shifted: index i -> power i (of q^{i})
    out = [MPQ0] * (N + 1)
    for i, c in enumerate(integ[:N]):
        out[i + 1] = c / (i + 1)
    return out


def scompose(f, g, N):
    """f(g(q)) with g[0] = 0; Horner."""
    if g and g[0] != 0:
        raise ValueError("compose needs g(0)=0")
    out = [MPQ0] * (N + 1)
    for c in reversed(f[:N + 1]):
        out = smul(out, g, N)
        out[0] += c
    return out


def srevert(f, N):
    """Compositional inverse of f = q + O(q^2) by Newton iteration."""
    if f[0] != 0 or f[1] != 1:
        raise ValueError("revert needs f = q + O(q^2)")
    g = [MPQ0, MPQ1] + [MPQ0] * (N - 1)
    df = [mpq(i) * f[i] for i in range(1, len(f))]
    prec = 2
    while prec <= N:
        prec = min(2 * prec, N + 1)
        # g <- g - (f(g) - q)/f'(g), computed to order prec-1
        fg = scompose(f, g, prec - 1)
        err = [fg[i] if i < len(fg) else MPQ0 for i in range(prec)]
        err[1] -= 1
        dfg = scompose(df, g, prec - 1)
        corr = sdiv(err, dfg, prec - 1)
        g = [(g[i] if i < len(g) else MPQ0) -
             (corr[i] if i < len(corr) else MPQ0) for i in range(prec)]
    return [g[i] if i < len(g) else MPQ0 for i in range(N + 1)]


def stheta(a, N):
    return [mpq(i) * a[i] for i in range(min(len(a), N + 1))] + \
        [MPQ0] * max(0, N + 1 - len(a))


# t-polynomials: dict m -> q-series (coefficient of t^m/m! is NOT used here;
# plain t^m coefficients).

def tp_theta(tp, N):
    """theta_q on sum_m t^m f_m(q) with t = log q: theta(t^m f) =
    t^m theta f + m t^{m-1} f."""
    out = {}
    for m, f in tp.items():
        tf = stheta(f, N)
        if any(tf):
            out[m] = padd_series(out.get(m), tf, N)
        if m >= 1:
            mf = pscal_series(mpq(m), f, N)
            if any(mf):
                out[m - 1] = padd_series(out.get(m - 1), mf, N)
    return {m: f for m, f in out.items() if any(f)}


def padd_series(a, b, N):
    if a is None:
        a = [MPQ0] * (N + 1)
    return [a[i] + (b[i] if i < len(b) else MPQ0) for i in range(N + 1)]


def pscal_series(k, a, N):
    return [k * (a[i] if i < len(a) else MPQ0) for i in range(N + 1)]


def tp_div_scalar_series(tp, den, N):
    di = sinv(den, N)
    return {m: smul(f, di, N) for m, f in tp.items()}


# -------------------------------------------------------------- kernel invert

def lambert_invert(series, w, N):
    """Given F = 1 + sum c_n q^n, solve F - 1 = sum_d m_d d^w q^d/(1-q^d):
    c_n = sum_{d|n} m_d d^w  =>  m_n exactly (mpq)."""
    m = {}
    for n in range(1, N + 1):
        acc = series[n] if n < len(series) else MPQ0
        for d in range(1, n):
            if n % d == 0 and d in m:
                acc -= m[d] * mpq(d) ** w
        m[n] = acc / (mpq(n) ** w)
    return m


# ---------------------------------------------------------- integrality reads

def integrality_verdict(numbers):
    """numbers: dict d -> mpq.  Script-emitted verdict string + data."""
    dens = {d: int(v.denominator) for d, v in numbers.items()}
    if all(x == 1 for x in dens.values()):
        return {"verdict": "integral", "rescale_N": 1,
                "failure_pattern": None}
    lcm = 1
    for x in dens.values():
        lcm = lcm * x // math.gcd(lcm, x)
    # failure pattern: prime-valuation growth of denominators vs d
    primes = set()
    for x in dens.values():
        y = x
        for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43):
            while y % p == 0:
                primes.add(p)
                y //= p
        if y > 1:
            primes.add(("large", y))
    pat = {}
    for p in [p for p in primes if isinstance(p, int)]:
        vals = []
        for d in sorted(dens):
            v, y = 0, dens[d]
            while y % p == 0:
                v += 1
                y //= p
            vals.append(v)
        pat[str(p)] = vals
    growth = "bounded" if all(
        isinstance(p, int) for p in primes) and lcm < 10 ** 12 else "unbounded-or-large"
    if lcm < 10 ** 12:
        return {"verdict": f"integral-after-rescaling-by-N={lcm}",
                "rescale_N": int(lcm),
                "failure_pattern": {"denominator_primes_valuations_by_d": pat,
                                    "growth": growth}}
    return {"verdict": "non-integral",
            "rescale_N": None,
            "failure_pattern": {"denominator_primes_valuations_by_d": pat,
                                "lcm_digits": len(str(lcm)),
                                "growth": growth}}


def q_rescale_kappa(series_coeffs):
    """Mirror-map N-integrality mode: minimal kappa in Z>0 with
    b_d * kappa^(d-1) in Z for all d (per-prime ceil of v_p(den)/(d-1)).
    series_coeffs: dict d -> mpq with d >= 1, b_1 = 1."""
    from collections import defaultdict
    need = defaultdict(int)
    for d, b in series_coeffs.items():
        if d < 2:
            continue
        den = int(b.denominator)
        y = den
        p = 2
        while p * p <= y:
            if y % p == 0:
                v = 0
                while y % p == 0:
                    v += 1
                    y //= p
                need[p] = max(need[p], -(-v // (d - 1)))
            p += 1
        if y > 1:
            need[y] = max(need[y], -(-1 // (d - 1)) if d > 1 else 1)
    kappa = 1
    for p, e in need.items():
        kappa *= p ** e
    return int(kappa)


def valuation_table(numbers, primes=(2, 3, 5, 7, 11, 13)):
    out = {}
    for p in primes:
        row = []
        for d in sorted(numbers):
            v = numbers[d]
            num = int(v.numerator)
            if num == 0:
                row.append(None)
                continue
            val, y = 0, abs(num)
            while y % p == 0:
                val += 1
                y //= p
            val_den, y = 0, int(v.denominator)
            while y % p == 0:
                val_den += 1
                y //= p
            row.append(val - val_den)
        out[str(p)] = row
    return out


# -------------------------------------------------------------- the main run

def fingerprint(R, rho, k, N, label, n0=1, overrides=None):
    """Full fingerprint at one MUM-type point.  R = theta-form dict,
    rho = chain bottom exponent (minimal class exponent), k = maximal block
    size, N = depth.  overrides: {param_id: value} un-freezes canonical-gauge
    parameters (sensitivity runs).  Returns exact data (mpq stringified)."""
    t0 = time.time()
    N_req = N
    sol = solve_chain(R, rho, k, N, label=label, overrides=overrides)
    a = sol["a"]
    b = sol["chain_bottom_relative_valuation"]
    if b > 0:
        N = N - b  # rebased series lose b top orders; work at N_eff
        a = [row[:N + 1] for row in a]
        rho = mpq(rho) + b

    # u_j series (relative to s^rho_eff); chain: y_m tower rows u_{k-1-m+r}
    u = [a[j] for j in range(k)]  # u[j][n]
    y0 = u[k - 1]
    if y0[0] == 0:
        raise ChainSolveError(f"[{label}] chain bottom has zero lead")

    # ratios v_{j,r} = u_{k-1-j+r}/y0 for chain member y_j
    y0inv = sinv(y0, N)

    # canonical chain re-gauge: if the solve FORCED a nonzero y0-admixture
    # constant c in y1 (v10(0) = c != 0), apply exp(-c N): u[j] <-
    # sum_i (-c)^i/i! u[j+i] — preserves the chain, zeroes the constant.
    regauge_c = smul(u[k - 2], y0inv, N)[0]
    if regauge_c != 0:
        factr = [MPQ1]
        for i in range(1, k + 1):
            factr.append(factr[-1] * i)
        unew = []
        for j in range(k):
            acc = [MPQ0] * (N + 1)
            for i in range(0, k - j):
                coef = (-regauge_c) ** i / factr[i]
                acc = [acc[n] + coef * u[j + i][n] for n in range(N + 1)]
            unew.append(acc)
        u = unew
        y0 = u[k - 1]
        y0inv = sinv(y0, N)

    # mirror map: y1/y0 = log s + v10, v10 = u_{k-2}/y0 (val >= 1 by gauge)
    v10 = smul(u[k - 2], y0inv, N)
    if v10[0] != 0:
        raise ChainSolveError(
            f"[{label}] v10(0) != 0 after re-gauge: {v10[0]}")
    qofs = smul([MPQ0, MPQ1] + [MPQ0] * (N - 1), sexp(v10, N), N)  # q(s)
    sofq = srevert(qofs, N)                                        # s(q)

    # A~(q) = v10(s(q)); log s = t - A~
    Atil = scompose(v10, sofq, N)

    # chain ratios in q with t-decomposition:
    # R_j = sum_r (t - A~)^r / r! * vt_{j,r}(q),  vt = v(s(q))
    # collect true t-poly; verify translation property; extract g_j.
    g = {0: [MPQ1] + [MPQ0] * N, 1: [MPQ0] * (N + 1)}
    check_translation = {}
    fact = [1] * (k + 1)
    for i in range(1, k + 1):
        fact[i] = fact[i - 1] * i
    tpolys = {}
    negA_pows = [[MPQ1] + [MPQ0] * N]
    for _ in range(k):
        negA_pows.append(smul(negA_pows[-1], [-c for c in Atil], N))
    for j in range(1, k):
        tp = {}
        for r in range(0, j + 1):
            vjr = smul(u[k - 1 - j + r], y0inv, N)
            vq = scompose(vjr, sofq, N)
            for i in range(0, r + 1):
                coef = MPQ1 / (fact[i] * fact[r - i])
                term = pscal_series(coef, smul(negA_pows[r - i], vq, N), N)
                tp[i] = padd_series(tp.get(i), term, N)
        # tp[i] = coefficient of t^i (true coefficient, includes 1/i!)
        tpolys[j] = {m: f for m, f in tp.items() if any(f)}
        g[j] = pscal_series(MPQ1, tpolys[j].get(0, [MPQ0] * (N + 1)), N)
        # translation check: coefficient of t^m must equal g_{j-m}/m!
        for m, f in list(tpolys[j].items()):
            if m >= 1 and (j - m) in g:
                want = pscal_series(MPQ1 / fact[m], g[j - m], N)
                diff = [f[i] - want[i] for i in range(N + 1)]
                check_translation[f"j{j}_t{m}"] = not any(diff)

    # alpha ladder
    alphas = {1: [MPQ1] + [MPQ0] * N}
    ladder_ok = True
    ladder_note = []
    for j in range(2, k):
        X = dict(tpolys[j])  # t-poly of R_j
        ok = True
        for i in range(1, j + 1):
            X = tp_theta(X, N)
            if i < j:
                X = tp_div_scalar_series(X, alphas[i], N)
        # X should now be t-free
        tfree = {m: f for m, f in X.items() if m != 0 and any(f)}
        if tfree:
            ok = False
            ladder_ok = False
            ladder_note.append(
                f"alpha_{j}: residual t-powers {sorted(tfree)} nonzero "
                f"(first nonzero orders: "
                f"{ {m: next((i for i, c in enumerate(f) if c), None) for m, f in tfree.items()} })")
        alphas[j] = X.get(0, [MPQ0] * (N + 1))
        if alphas[j][0] != 1:
            ladder_note.append(
                f"alpha_{j}(0) = {alphas[j][0]} (expected 1)")

    # instanton-type numbers from alpha_2 (and higher), kernels w = 2,3,4
    inst = {}
    for j in sorted(alphas):
        if j == 1:
            continue
        row = {}
        for w in (2, 3, 4):
            nums = lambert_invert(alphas[j], w, N)
            nums = {d: mpq(n0) * v for d, v in nums.items()}
            row[f"kernel_d^{w}"] = {
                "numbers_first": {str(d): str(nums[d])
                                  for d in sorted(nums) if d <= min(N, 24)},
                "integrality": integrality_verdict(nums),
            }
        inst[f"alpha_{j}"] = row

    # mirror map data
    mm = {str(d): str(sofq[d]) for d in range(1, min(N, 24) + 1)}
    mm_nums = {d: sofq[d] for d in range(1, N + 1)}
    mm_int = integrality_verdict(mm_nums)
    mm_kappa = q_rescale_kappa(mm_nums)

    dt = time.time() - t0
    return {
        "label": label,
        "rho_effective": str(rho),
        "chain_length_k": k,
        "depth_N_requested": N_req,
        "depth_N": N,
        "chain_bottom_relative_valuation":
            sol["chain_bottom_relative_valuation"],
        "gauge_killed": sol["gauge_killed"],
        "gauge_unkillable_constants": sol["gauge_unkillable_constants"],
        "resonances_relative_orders": sol["resonances"],
        "params_introduced": sol["params_introduced"],
        "params_eliminated_by_obstructions": sol["params_eliminated"],
        "params_frozen_to_zero_canonical_gauge": sol["params_frozen_to_zero"],
        "y0_head": {str(n): str(y0[n]) for n in range(min(N, 12) + 1)},
        "mirror_map_s_of_q_head": mm,
        "mirror_map_integrality": mm_int,
        "mirror_map_q_rescale_kappa": mm_kappa,
        "q_of_s_head": {str(n): str(qofs[n]) for n in range(1, min(N, 12) + 1)},
        "translation_property_checks": check_translation,
        "alpha_ladder_log_free": ladder_ok,
        "alpha_ladder_notes": ladder_note,
        "alphas_head": {f"alpha_{j}": {str(n): str(alphas[j][n])
                                       for n in range(min(N, 16) + 1)}
                        for j in alphas},
        "instanton_type": inst,
        "n0_frame": n0,
        "wall_seconds": round(dt, 3),
        "_series": {  # exact full series for downstream use (strings)
            "s_of_q": [str(c) for c in sofq],
            "q_of_s": [str(c) for c in qofs],
            "alphas": {str(j): [str(c) for c in alphas[j]] for j in alphas},
            "g": {str(j): [str(c) for c in g[j]] for j in g if j >= 2},
            "y0": [str(c) for c in y0],
        },
    }


# ------------------------------------------------------------- receipt helper

def producer_block(script_path):
    return {
        "tool": "eichler.mirror6",
        "script": script_path,
        "script_sha256": sha256_file(script_path),
        "stamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }


def emit_receipt(path, payload, script_path):
    payload = dict(payload)
    payload["producer"] = producer_block(script_path)
    with open(path, "w") as f:
        json.dump(payload, f, indent=1, default=str)
    return path
