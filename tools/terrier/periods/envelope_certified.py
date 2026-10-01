#!/usr/bin/env python3
"""
envelope_certified.py — PROVEN tail envelope for Frobenius log-tower series
of a Fuchsian operator (lemma + proof: PROVEN_MAJORANT.md, this directory).

Replaces the EMPIRICAL last-window x4 tail of
upgrades/Eichler.jl/src/cy_transport.jl (mum_frobenius_basis), which is
unsound (two exact counterexamples in PROVEN_MAJORANT.md section 5).
Everything here is exact Fraction arithmetic; floats appear only in
__repr__/reporting.

API:
  theta_form(coeffs)                p_i(z) coeff lists -> (R, smin)
  mum_jets(R, e, J, N)              exact rho-jet recursion (MUM/nonresonant m>=1)
  rational_roots(R0)                exact factorization (lc, {root: mult}) or None
  phi(R, e, J, lc, roots, t, m)     the certified Phi(m, t) of the lemma (exact)
  certify(R, e, J, jets, N, x, ...) -> dict: proven tail bound per jet level
  deriv_tail(cert, d)               corollary (C3) weight-m^d tail

Resource caps: single-core, exact; run under ulimit -v 32505856, nice 5.
Inputs: theta-forms from a stored operator JSON (e.g. pipeline/dkmm/
operator_LS.json), stored log-layer towers (the DKMM certified-W0 towers of
the reference data, TERRIER_KKLT_BANK, not included in the package), banana
operators transcribed from Eichler.jl banana_pf (cy_transport.jl).
"""
from fractions import Fraction as Fr
from math import comb

# ---------------------------------------------------------------- theta form
def theta_form(coeffs):
    """coeffs[i][j] = [z^j] p_i(z) (exact ints/Fractions). Returns (R, smin) with
    R[s][k] = [theta^k] R_s, z^{-smin} L = sum_s z^s R_s(theta)."""
    r = len(coeffs) - 1
    pairs = [(i, j) for i in range(r + 1) for j in range(len(coeffs[i]))
             if coeffs[i][j]]
    if not pairs:
        raise ValueError("zero operator")
    smin = min(j - i for i, j in pairs)
    smax = max(j - i for i, j in pairs)
    ff = [[Fr(1)]]                       # falling factorials theta^(i)
    for i in range(1, r + 1):
        prev, cur = ff[-1], [Fr(0)] * (i + 1)
        for k, c in enumerate(prev):
            cur[k + 1] += c
            cur[k] -= (i - 1) * c
        ff.append(cur)
    R = [[Fr(0)] * (r + 1) for _ in range(smax - smin + 1)]
    for i, j in pairs:
        for k, fk in enumerate(ff[i]):
            R[j - i - smin][k] += Fr(coeffs[i][j]) * fk
    return R, smin

# ---------------------------------------------------------------- jet helpers
def _jmul(a, b, J):
    c = [Fr(0)] * J
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b[:J - i]):
                c[i + j] += ai * bj
    return c

def _jinv(a, J):
    if a[0] == 0:
        raise ZeroDivisionError("jet not invertible (resonance)")
    c = [Fr(0)] * J
    c[0] = 1 / a[0]
    for k in range(1, J):
        c[k] = -c[0] * sum(a[j] * c[k - j] for j in range(1, k + 1))
    return c

def _jpoly(P, u, J):
    """P(u + rho) as a length-J jet (exact)."""
    res, cur = [Fr(0)] * J, [Fr(0)] * J
    cur[0] = Fr(1)
    for k in range(len(P)):
        if P[k]:
            for j in range(J):
                res[j] += Fr(P[k]) * cur[j]
        nxt = [Fr(0)] * J
        for j in range(J):
            if cur[j]:
                nxt[j] += Fr(u) * cur[j]
                if j + 1 < J:
                    nxt[j + 1] += cur[j]
        cur = nxt
    return res

def l1(v):
    return sum(abs(x) for x in v)

def mum_jets(R, e, J, N):
    """b_m, m=0..N, from REC with b_0 = 1 (requires R_0(e+m) != 0 for m>=1)."""
    S = len(R) - 1
    b = [[Fr(0)] * J for _ in range(N + 1)]
    b[0][0] = Fr(1)
    for m in range(1, N + 1):
        acc = [Fr(0)] * J
        for s in range(1, min(S, m) + 1):
            if any(R[s]):
                t = _jmul(_jpoly(R[s], Fr(e) + m - s, J), b[m - s], J)
                acc = [x + y for x, y in zip(acc, t)]
        d = _jinv(_jpoly(R[0], Fr(e) + m, J), J)
        b[m] = [-x for x in _jmul(d, acc, J)]
    return b

# ---------------------------------------------------------------- roots of R0
def rational_roots(R0, num_range=2000, dens=(1, 2, 3, 4, 5, 6, 8, 10, 12)):
    """Exact rational factorization of R0. Returns (lc, {root: mult}) or None."""
    P = [Fr(c) for c in R0]
    while P and P[-1] == 0:
        P.pop()
    lc, roots = P[-1], {}
    def val(Q, x):
        v = Fr(0)
        for c in reversed(Q):
            v = v * x + c
        return v
    while len(P) > 1:
        found = None
        for num in range(-num_range, num_range + 1):
            for den in dens:
                x = Fr(num, den)
                if val(P, x) == 0:
                    found = x
                    break
            if found is not None:
                break
        if found is None:
            return None
        roots[found] = roots.get(found, 0) + 1
        q, acc = [], Fr(0)
        for c in reversed(P):
            acc = c + acc * found
            q.append(acc)
        assert q[-1] == 0
        P = list(reversed(q[:-1]))
    return lc, roots

# ---------------------------------------------------------------- the lemma
def phi(R, e, J, lc, roots, t, m):
    """Phi(m, t) of PROVEN_MAJORANT.md section 2 — exact Fraction."""
    S = len(R) - 1
    r = sum(roots.values())
    beta = max([Fr(0)] + [Fr(lam) - Fr(e) for lam in roots])
    mb = Fr(m) - beta
    if mb <= 0:
        raise ValueError("m <= beta: below the certified horizon")
    G = Fr(1)
    for lam, mu in roots.items():
        G *= sum(Fr(comb(mu + j - 1, j)) / mb ** j for j in range(J))
    tot = Fr(0)
    for s in range(1, S + 1):
        cs = max(Fr(0), Fr(e) - s + 1)
        inner = sum(abs(Fr(R[s][k])) * (Fr(m) + cs) ** k for k in range(len(R[s])))
        if inner:
            tot += inner / Fr(t) ** s
    return G * tot / (abs(lc) * mb ** r)

def _t_min(R, e, J, lc, roots, N, iters=60, thi=None):
    """Smallest t with Phi(N+1, t) <= 1, by exact bisection. None if none < thi."""
    lo, hi = Fr(0), None
    probe = Fr(1)
    for _ in range(80):                      # find any admissible upper end
        if phi(R, e, J, lc, roots, probe, N + 1) <= 1:
            hi = probe
            break
        probe *= 2
    if hi is None:
        return None
    for _ in range(iters):
        mid = (lo + hi) / 2
        if mid > 0 and phi(R, e, J, lc, roots, mid, N + 1) <= 1:
            hi = mid
        else:
            lo = mid
    return hi

def certify(R, e, J, jets, N, x, t=None, n_grid=8):
    """PROVEN tail envelope at 0 < x < 1/t for the jet tower `jets` (exact seeds
    b_m, m <= N, satisfying REC past their horizon; caller vouches for REC —
    use `rec_check` below to verify a stretch).

    Returns dict with exact Fractions:
      t, K, Phi(N+1,t), tail  = sum_{m>N} ||b_m|| x^m bound (C2),
      per_level[j]            = same bound (valid for each jet component),
      ok = False (with reason) if no admissible t < 1/x exists at this N."""
    S = len(R) - 1
    rr = rational_roots(R[0])
    if rr is None:
        return {"ok": False, "reason": "indicial roots not rational — supply beta bound"}
    lc, roots = rr
    r = sum(roots.values())
    if len(R[0]) - 1 > r or all(c == 0 for c in R[0]):
        return {"ok": False, "reason": "deg R_0 != r (not regular singular form)"}
    beta = max([Fr(0)] + [Fr(lam) - Fr(e) for lam in roots])
    from math import ceil
    if not (N > beta and N >= S + max(0, ceil(-Fr(e)))):
        return {"ok": False, "reason": f"(H1) fails: need N > {beta}, N >= S + ceil(-e)"}
    x = Fr(x)
    cands = []
    tmin = _t_min(R, e, J, lc, roots, N) if t is None else Fr(t)
    if tmin is None or tmin * x >= 1:
        return {"ok": False, "reason": f"no admissible t < 1/x at N={N} "
                f"(t_min={float(tmin) if tmin else 'inf'}); raise N or shrink x"}
    if t is not None:
        cands = [Fr(t)]
    else:
        eps = (1 / x - tmin) / (n_grid + 1)
        cands = [tmin * Fr(2 ** 20 + 1, 2 ** 20)] + \
                [tmin + k * eps for k in range(1, n_grid)]
        cands = [c for c in cands if c * x < 1]
    best = None
    for tc in cands:
        if phi(R, e, J, lc, roots, tc, N + 1) > 1:
            continue
        K = max(l1(jets[m]) / tc ** m for m in range(N - S + 1, N + 1))
        tail = K * (tc * x) ** (N + 1) / (1 - tc * x)
        if best is None or tail < best["tail"]:
            best = {"ok": True, "t": tc, "K": K, "N": N, "x": x, "J": J,
                    "beta": beta, "t_min": tmin, "tail": tail,
                    "phi": phi(R, e, J, lc, roots, tc, N + 1)}
    if best is None:
        return {"ok": False, "reason": "no candidate t passed (H2)"}
    return best

def deriv_tail(cert, d):
    """(C3): proven bound on sum_{m>N} ||b_m|| m^d x^m. None if q >= 1."""
    from fractions import Fraction as F
    t, x, N, K = cert["t"], cert["x"], cert["N"], cert["K"]
    # rational upper bound on e^{d/(N+1)}: sum_{i<=8} (d/(N+1))^i/i! * (1 + tiny pad)
    u = F(d, N + 1)
    epad = sum(u ** i / __import__("math").factorial(i) for i in range(9))
    epad = epad * (1 - u / 9) ** -1 if u < 9 else None   # geometric remainder pad
    if epad is None:
        return None
    q = t * x * epad
    if q >= 1:
        return None
    return K * F(N + 1) ** d * (t * x) ** (N + 1) / (1 - q)

def rec_check(R, e, J, jets, m_from, m_to, picture="mult"):
    """Exact check that jets satisfy (REC) on [m_from, m_to].
    picture='mult': jet-product form (cy_transport). picture='div': divided-power
    operator form (the cert_w0 log layers, gamma_j = j! * f_j)."""
    S = len(R) - 1
    for m in range(m_from, m_to + 1):
        if picture == "mult":
            lhs = _jmul(_jpoly(R[0], Fr(e) + m, J), jets[m], J)
            rhs = [Fr(0)] * J
            for s in range(1, min(S, m) + 1):
                t = _jmul(_jpoly(R[s], Fr(e) + m - s, J), jets[m - s], J)
                rhs = [a - b for a, b in zip(rhs, t)]
        else:
            def apply(P, u, v):
                out = [Fr(0)] * J
                for k in range(len(P)):
                    if P[k]:
                        for j in range(J):
                            acc = sum(Fr(comb(k, i)) * Fr(u) ** (k - i) * v[j + i]
                                      for i in range(0, min(k, J - 1 - j) + 1))
                            out[j] += Fr(P[k]) * acc
                return out
            lhs = apply(R[0], Fr(e) + m, jets[m])
            rhs = [Fr(0)] * J
            for s in range(1, min(S, m) + 1):
                t = apply(R[s], Fr(e) + m - s, jets[m - s])
                rhs = [a - b for a, b in zip(rhs, t)]
        if lhs != rhs:
            return m
    return None

def empirical_tail_x4(coeffs, N, x, win):
    """EXACT-rational emulation of cy_transport.jl lines 278-288 (the unsound
    empirical scheme), for one jet level's coefficient list."""
    ratio = Fr(0)
    for m in range(max(1, N - win), N):
        if coeffs[m] != 0:
            ratio = max(ratio, abs(Fr(coeffs[m + 1]) / Fr(coeffs[m])))
    th = ratio * Fr(x)
    last = abs(Fr(coeffs[N])) * Fr(x) ** N
    return 4 * last * th / (1 - th) if (th < Fr(95, 100) and last > 0) else Fr(0)

def fmt(v):
    """Format a Fraction (or float) as a.bbe+-xx without float under/overflow."""
    from math import log10, floor
    v = Fr(v)
    if v == 0:
        return "0"
    lg = log10(abs(v.numerator)) - log10(v.denominator)
    ex = floor(lg)
    mant = 10 ** (lg - ex)
    return f"{'-' if v < 0 else ''}{mant:.2f}e{ex:+03d}"
