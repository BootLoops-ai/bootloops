#!/usr/bin/env python3
"""asdlib.py — Atkin–Swinnerton-Dyer congruence miner core (member of the
frobenius-boundary package).

Law mined:
    a(np) - gamma_p a(n) + p^(k-1) a(n/p) == 0 (mod p^(s(n))),
    a(n/p) := 0 if p !| n;  floor law s(n) = 1 + ord_p(n).
Everything here is deterministic; all verdict strings are computed, never set.

GAUGE LESSON (load-bearing; measured before this member shipped): family-MUM
mp^r-1 coefficient grids mine the TRIVIAL unit root — CONSISTENT-TRIVIAL is
not "modularity found"; proven nontrivial three-term gamma_p lives on
fixed-variety formal groups (Honda class).  classify_gamma() separates the
classes mechanically.
"""
import json, hashlib
from fractions import Fraction
from math import isqrt

# ---------------------------------------------------------------- valuations

def vp_int(x, p):
    """p-adic valuation of a nonzero int; None for 0 (=+infinity)."""
    if x == 0:
        return None
    v = 0
    while x % p == 0:
        x //= p
        v += 1
    return v

def vp_of(x, p):
    if isinstance(x, Fraction):
        if x == 0:
            return None
        return vp_int(x.numerator, p) - vp_int(x.denominator, p)
    return vp_int(x, p)

def red_mod(x, m, p):
    """Reduce int/Fraction mod m = p^t. Raises ValueError if p divides the
    denominator (P-DENOMINATOR event; caller converts to the wall verdict)."""
    if isinstance(x, Fraction):
        num, den = x.numerator, x.denominator
        if den % p == 0:
            raise ValueError("p-denominator")
        return (num % m) * pow(den % m, -1, m) % m
    return x % m

# ---------------------------------------------------------------- the miner

def _rows(p, N, k, aval, tier_shift):
    """Yield per-row constraints (n, e, residue) with gamma == residue mod p^e,
    plus degenerate/inconsistent row records. aval(n) -> exact int/Fraction.
    tier_shift: 0 = guaranteed floor law, 1 = exploratory (+1)."""
    out = {"cons": [], "degen": [], "incons": [], "pden": []}
    pk1 = p ** (k - 1)
    for n in range(1, N // p + 1):
        vn = 0
        nn = n
        while nn % p == 0:
            nn //= p
            vn += 1
        t = 1 + vn + tier_shift
        m = p ** t
        try:
            an = aval(n); anp = aval(n * p)
            A = red_mod(anp, m, p)
            if n % p == 0:
                A = (A + pk1 * red_mod(aval(n // p), m, p)) % m
            B = red_mod(an, m, p)
        except ValueError:
            out["pden"].append(n)
            continue
        if B == 0:
            # gamma*0 == A (mod m): constraint only if A != 0
            if A % m != 0:
                out["incons"].append({"n": n, "why": "a(n)==0 mod p^t but rhs!=0",
                                      "t": t})
            else:
                out["degen"].append(n)
            continue
        v = vp_int(B, p)  # 0 <= v < t
        if v:
            if A % p ** v != 0:
                out["incons"].append({"n": n, "why": "val(rhs)<val(a(n))",
                                      "t": t})
                continue
            A //= p ** v
            B //= p ** v
            m = p ** (t - v)
            A %= m
            B %= m
        g = A * pow(B, -1, m) % m
        out["cons"].append((n, t - v, g))
    return out

def _merge(cons, p):
    """CRT-merge single-prime constraints [(n, e, residue)].
    Returns (D, gamma, witnesses): if witnesses nonempty the system is
    inconsistent; gamma is the deepest row's value (mod p^D)."""
    if not cons:
        return 0, None, []
    best = max(cons, key=lambda r: r[1])
    D, g = best[1], best[2]
    wit = []
    for (n, e, r) in cons:
        if (g - r) % p ** e != 0:
            wit.append({"n": n, "e": e, "row_gamma": r,
                        "deep_gamma_mod_pe": g % p ** e, "deep_row_n": best[0]})
    return D, g, wit

def weil_candidates(g, p, D, k):
    """Integers c == g mod p^D with |c| <= floor(2 p^((k-1)/2))."""
    if D == 0 or g is None:
        return []
    m = p ** D
    b = isqrt(4 * p ** (k - 1))
    c0 = g % m
    cands = []
    c = c0 - ((c0 + b) // m) * m  # smallest rep >= -b
    while c <= b:
        if c >= -b:
            cands.append(c)
        c += m
    return cands

def drift_fit(cons, p):
    """Gauge-drift detector: on plain rows with e>=2 fit
    gamma_n == g0*(1 + n*eps*p) mod p^2. Returns dict."""
    rows = [(n, r % p * 0 + r) for (n, e, r) in cons if e >= 2]
    rows = [(n, r) for (n, r) in rows]
    if len(rows) < 3:
        return {"fit": "insufficient rows"}
    m2 = p * p
    (n1, g1), (n2, g2) = rows[0], rows[1]
    if (g1 - g2) % p != 0:
        return {"fit": "no mod-p base agreement"}
    # gamma_n = g0 + g0*eps*p*n: mod p: g0 == g1 mod p. Write g_n = g1row...
    # Solve: g_n - g_m == g0*eps*p*(n-m) mod p^2 -> eps0 := g0*eps mod p
    dn = (n1 - n2) % p
    if dn == 0:
        return {"fit": "degenerate n spacing"}
    eps0 = ((g1 - g2) // p) % p * pow(dn, -1, p) % p
    g0 = (g1 - eps0 * p * n1) % m2
    ok, bad = 0, 0
    for (n, g) in rows:
        if (g0 + eps0 * p * n - g) % m2 == 0:
            ok += 1
        else:
            bad += 1
    return {"fit": "solved", "eps0_mod_p": eps0, "g0_mod_p2": g0,
            "rows_ok": ok, "rows_bad": bad,
            "verdict_drift": ("DRIFT-CONSISTENT" if bad == 0 and eps0 != 0 else
                              "NO-DRIFT" if bad == 0 else "DRIFT-FAIL")}

def emp_depth_profile(aval, p, N, k, c, cap=None):
    """Empirical depth d(n) = ord_p(a(np) - c a(n) + p^(k-1) a(n/p)) per row,
    integer candidate c. cap = precision floor (modular series) or None."""
    pk1 = p ** (k - 1)
    plain, deep = [], []
    for n in range(1, N // p + 1):
        r = aval(n * p) - c * aval(n)
        if n % p == 0:
            r = r + pk1 * aval(n // p)
        d = vp_of(r, p)
        d = (cap if d is None else min(d, cap)) if cap is not None else d
        (deep if n % p == 0 else plain).append((n, d))
    def summ(rows):
        if not rows:
            return {"rows": 0}
        ds = [d for (_, d) in rows]
        fin = [d for d in ds if d is not None]
        return {"rows": len(rows),
                "min": (min(fin) if fin else "exact-all"),
                "n_exact_zero_rows": sum(1 for d in ds if d is None),
                "max_capped": cap}
    return {"plain": summ(plain), "deep": summ(deep)}

def mine(aval, p, N, k, want_profile=True):
    """Full mine at one prime. aval(n) exact int/Fraction, 1<=n<=N."""
    res = {"p": p, "N": N, "k": k, "rows_max": N // p}
    g_final = None
    for tier, shift in (("guaranteed", 0), ("exploratory", 1)):
        rr = _rows(p, N, k, aval, shift)
        D, g, wit = _merge(rr["cons"], p)
        v = ("P-DENOMINATOR-WALL" if (rr["pden"] and not rr["cons"]) else
             "UNDERDETERMINED" if not rr["cons"] else
             "INCONSISTENT" if (wit or rr["incons"]) else
             "CONSISTENT-WITH-PDEN-ROWS" if rr["pden"] else
             "CONSISTENT")
        entry = {"verdict": v, "depth_D": D, "rows_used": len(rr["cons"]),
                 "rows_degenerate": len(rr["degen"]),
                 "gamma_mod_pD": (g if v.startswith("CONSISTENT") else None),
                 "witnesses": (wit + rr["incons"])[:8],
                 "n_witnesses": len(wit) + len(rr["incons"]),
                 "pden_rows": rr["pden"][:4]}
        if v.startswith("CONSISTENT"):
            entry["weil_box_k"] = weil_candidates(g, p, D, k)
            entry["weil_unique"] = (len(entry["weil_box_k"]) == 1)
            m = p ** D
            eis = (1 + p ** (k - 1)) % m
            entry["ordinary"] = (g % p != 0)
            entry["pins"] = {"one": g % m == 1 % m, "zero": g % m == 0,
                             "eis_1_plus_pk1": g % m == eis,
                             "neg_eis": g % m == (-(1 + p ** (k - 1))) % m}
            entry["gamma_class"] = classify_gamma(entry, D, k)
        if tier == "exploratory" and v == "INCONSISTENT":
            entry["drift_detector"] = drift_fit(rr["cons"], p)
        res[tier] = entry
        if tier == "guaranteed" and v.startswith("CONSISTENT"):
            g_final = (g, D)
    if want_profile and g_final:
        c = profile_candidate(res["guaranteed"], p, k)
        if c is not None:
            res["integer_candidate_k"] = c
            res["empirical_depth"] = emp_depth_profile(aval, p, N, k, c)
    return res

def classify_gamma(entry, D, k):
    """Script-computed gamma classification (fixed decision grid)."""
    pins = entry["pins"]
    if pins["zero"] and not entry["ordinary"]:
        return "ZERO/NON-ORDINARY"
    if pins["one"] and pins["eis_1_plus_pk1"]:
        return f"TRIVIAL-CLASS(1==1+p^{k-1} at depth {D})"
    if pins["eis_1_plus_pk1"]:
        return f"TRIVIAL-EISENSTEIN(1+p^{k-1})"
    if pins["one"]:
        return "TRIVIAL(1)"
    if entry.get("weil_unique"):
        return f"NONTRIVIAL-CUSPIDAL({entry['weil_box_k'][0]})"
    if not entry.get("weil_box_k"):
        return "OUTSIDE-CUSPIDAL-BOX"
    return f"AMBIGUOUS({len(entry['weil_box_k'])} candidates)"

def profile_candidate(entry, p, k):
    """Integer gamma used for the empirical depth profile."""
    if entry.get("weil_unique"):
        return entry["weil_box_k"][0]
    pins = entry.get("pins", {})
    if pins.get("eis_1_plus_pk1"):
        return 1 + p ** (k - 1)
    if pins.get("one"):
        return 1
    if pins.get("zero"):
        return 0
    return None

# ------------------------------------------------------- control generators

def apery_zeta3(nmax):
    """A(0..nmax) exact; recurrence (n+1)^3 A(n+1)=(2n+1)(17n^2+17n+5)A(n)-n^3 A(n-1),
    self-checked against the double sum for n<=40."""
    A = [1, 5]
    for n in range(1, nmax):
        num = (2 * n + 1) * (17 * n * n + 17 * n + 5) * A[n] - n ** 3 * A[n - 1]
        d = (n + 1) ** 3
        assert num % d == 0, f"Apery recurrence non-integral at n={n}"
        A.append(num // d)
    # double-sum spot check
    from math import comb
    for n in range(0, min(41, nmax + 1)):
        s = sum(comb(n, kk) ** 2 * comb(n + kk, kk) ** 2 for kk in range(n + 1))
        assert s == A[n], f"Apery self-check fail at n={n}"
    return A[:nmax + 1]

def eta_product_8_4(nmax):
    """q-expansion of eta(2t)^4 eta(4t)^4 to q^nmax (list a[0..nmax], a[0]=0)."""
    def eta_pow4(step, N):
        # (prod (1-q^(step*n)))^4 truncated at N: repeated sparse mult
        e = [0] * (N + 1); e[0] = 1
        # Euler pentagonal for prod(1-x^n), x = q^step
        f = [0] * (N + 1); f[0] = 1
        kk = 1
        while True:
            g1 = kk * (3 * kk - 1) // 2 * step
            g2 = kk * (3 * kk + 1) // 2 * step
            if g1 > N and g2 > N:
                break
            s = -1 if kk % 2 else 1
            if g1 <= N:
                f[g1] += s
            if g2 <= N:
                f[g2] += s
            kk += 1
        # f = prod(1-q^{step n}); raise to 4th by two squarings
        def mul(u, v):
            w = [0] * (N + 1)
            for i, ui in enumerate(u):
                if ui:
                    for j in range(0, N + 1 - i):
                        if v[j]:
                            w[i + j] += ui * v[j]
            return w
        f2 = mul(f, f)
        return mul(f2, f2)
    N = nmax - 1  # after multiplying by q
    A = eta_pow4(2, N)
    B = eta_pow4(4, N)
    C = [0] * (N + 1)
    for i, ai in enumerate(A):
        if ai:
            for j in range(0, N + 1 - i):
                if B[j]:
                    C[i + j] += ai * B[j]
    out = [0] * (nmax + 1)
    for i, ci in enumerate(C):
        out[i + 1] = ci
    return out

def apery_zeta2(nmax):
    """b(0..nmax): zeta(2)-Apery numbers (probe/report block only)."""
    from math import comb
    return [sum(comb(n, kk) ** 2 * comb(n + kk, kk) for kk in range(n + 1))
            for n in range(nmax + 1)]

def ell_omega_series(ainvs, N, m):
    """c(1..N) mod m (m odd): coefficients of the invariant differential
    omega = dx/(2y+a1x+a3) of y^2+a1xy+a3y = x^3+a2x^2+a4x+a6 in the
    Z-integral formal parameter t = -x/y (Silverman IV.1). Honda's theorem
    gives the ASD law with gamma_p = a_p(E), k=2, at good p."""
    a1, a2, a3, a4, a6 = [a % m for a in ainvs]
    NN = N + 6
    s = [0] * (NN + 1)      # w(t) = sum s_k t^k, k>=3
    w2 = [0] * (NN + 1)     # w^2 coefficients
    w3 = [0] * (NN + 1)     # w^3 coefficients
    for k in range(3, NN + 1):
        # w2[k], w3[k] involve only s_j, j <= k-3 (all set)
        v2 = 0
        for i in range(3, k - 2):
            if s[i]:
                j = k - i
                if 3 <= j <= NN and s[j]:
                    v2 += s[i] * s[j]
        w2[k] = v2 % m
        v3 = 0
        for i in range(3, k - 5):
            if s[i] and w2[k - i]:
                v3 += s[i] * w2[k - i]
        # NOTE w2 here is complete for indices <= k-3 (needed j=k-i>=6 ok)
        w3[k] = v3 % m
        val = (1 if k == 3 else 0)
        val += a1 * s[k - 1] + a2 * s[k - 2] + a3 * w2[k]
        val += a4 * (w2[k - 1] if k >= 1 else 0) + a6 * w3[k]
        s[k] = val % m
    # u = w/t^3, unit; uinv
    u = [s[k + 3] for k in range(0, N + 3)]
    u[0] = 1
    ui = [0] * len(u)
    ui[0] = 1
    for k in range(1, len(u)):
        acc = 0
        for j in range(1, k + 1):
            if u[j] and ui[k - j]:
                acc += u[j] * ui[k - j]
        ui[k] = (-acc) % m
    # num = -2*ui + t*ui'; g = t^3*(2y+a1x+a3) = -2*ui + a1*t*ui + a3 t^3
    num = [0] * len(ui)
    g = [0] * len(ui)
    for k in range(len(ui)):
        num[k] = (-2 * ui[k] + (k * ui[k])) % m  # t*ui' has coeff k*ui[k] at t^k
        g[k] = (-(2 * ui[k]) + (a1 * ui[k - 1] if k >= 1 else 0)) % m
    if len(g) > 3:
        g[3] = (g[3] + a3) % m
    # omega = num * g^{-1}
    gi = [0] * len(g)
    g0inv = pow(g[0], -1, m)
    gi[0] = g0inv
    for k in range(1, len(g)):
        acc = 0
        for j in range(1, k + 1):
            if g[j] and gi[k - j]:
                acc += g[j] * gi[k - j]
        gi[k] = (-acc * g0inv) % m
    om = [0] * (N + 1)
    for k in range(N + 1):
        acc = 0
        for j in range(k + 1):
            if num[j] and gi[k - j]:
                acc += num[j] * gi[k - j]
        om[k] = acc % m
    # c(n) = [t^(n-1)] omega  -> return list c[1..N]
    return [None] + om[:N]

# ------------------------------------------------------- operator machinery

def load_operator(path):
    d = json.load(open(path))
    coeffs = [[int(v) for v in ci] for ci in d["coeffs"]]
    return d, coeffs

def shift_operator(coeffs, x0):
    """P_i(x) -> P_i(u + x0) exact (binomial transform); x0 int."""
    if x0 == 0:
        return [list(ci) for ci in coeffs]
    out = []
    for ci in coeffs:
        n = len(ci)
        sh = [0] * n
        # Horner in (u + x0): P(u+x0) = (((c_{n-1})(u+x0) + c_{n-2})...)
        acc = [0] * n
        deg = 0
        for c in reversed(ci):
            # acc = acc*(u+x0) + c
            new = [0] * (deg + 2)
            for i in range(deg + 1):
                if acc[i]:
                    new[i + 1] += acc[i]
                    new[i] += acc[i] * x0
            new[0] += c
            deg = deg + 1
            acc = new[:deg + 1]
        acc = acc[:n] + [0] * max(0, n - len(acc))
        out.append(acc[:n])
    return out

def ffall(s, order):
    """[ff(s,0..order)] falling factorials, s int."""
    out = [1]
    for i in range(order):
        out.append(out[-1] * (s - i))
    return out

def op_local_data(coeffs):
    """d0, dmax, indicial polynomial I(s) as coefficient list (exact ints),
    and its integer roots."""
    order = len(coeffs) - 1
    d0 = min(min(j for j, v in enumerate(ci) if v) - i
             for i, ci in enumerate(coeffs) if any(ci))
    dmax = max(max(j for j, v in enumerate(ci) if v) - i
               for i, ci in enumerate(coeffs) if any(ci))
    # I(s) = sum_i ff(s,i)*c_i[i+d0]  — build as polynomial in s
    Ipoly = [0] * (order + 2)
    for i, ci in enumerate(coeffs):
        j = i + d0
        if 0 <= j < len(ci) and ci[j]:
            # ff poly: s(s-1)...(s-i+1)
            fp = [1]
            for t in range(i):
                fp = [a - t * b for a, b in
                      zip([0] + fp, fp + [0])]  # fp*(s-t): shift minus t*fp
            for e, cf in enumerate(fp):
                Ipoly[e] += cf * ci[j]
    while len(Ipoly) > 1 and Ipoly[-1] == 0:
        Ipoly.pop()
    # integer roots via rational root theorem (monic-insensitive scan)
    def ival(s):
        v = 0
        for cf in reversed(Ipoly):
            v = v * s + cf
        return v
    roots = [s for s in range(-30, 31) if ival(s) == 0]
    return {"d0": d0, "dmax": dmax, "Ipoly": Ipoly, "int_roots_pm30": roots,
            "ival": ival, "order": order}

def frobenius_exact(coeffs, rho, N):
    """Exact exponent-rho Frobenius solution b_0..b_N (Fractions), b_0=1.
    Resonances: RHS must vanish exactly, coefficient set to 0 (gauge law);
    else raises RuntimeError('LOG-OBSTRUCTED at n=..')."""
    loc = op_local_data(coeffs)
    d0, dmax, ival, order = loc["d0"], loc["dmax"], loc["ival"], loc["order"]
    BW = dmax - d0
    b = [Fraction(1)]
    ffc = {0: ffall(rho, order)}
    resonances = []
    for n in range(1, N + 1):
        rhs = Fraction(0)
        for m in range(max(0, n - BW), n):
            if b[m] == 0:
                continue
            dd = d0 + n - m
            w = 0
            fm = ffc[m]
            for i, ci in enumerate(coeffs):
                j = i + dd
                if 0 <= j < len(ci) and ci[j]:
                    w += fm[i] * ci[j]
            if w:
                rhs += b[m] * w
        Iv = ival(rho + n)
        if Iv == 0:
            if rhs != 0:
                raise RuntimeError(f"LOG-OBSTRUCTED at n={n}")
            resonances.append(n)
            b.append(Fraction(0))
        else:
            b.append(-rhs / Iv)
        ffc[n] = ffall(rho + n, order)
    return b, resonances, loc

def frobenius_modular(coeffs, rho, N, p, S_target, resonances_cert):
    """Exponent-rho solution mod p^F via working modulus p^(S_target+V_total).
    resonances_cert: list of resonant n (from the exact route) where b_n=0.
    Returns (b_res list mod p^F, F, W, V_total). Detects p-denominators."""
    loc = op_local_data(coeffs)
    d0, dmax, ival, order = loc["d0"], loc["dmax"], loc["ival"], loc["order"]
    BW = dmax - d0
    vs = []
    for n in range(1, N + 1):
        Iv = ival(rho + n)
        vs.append(0 if Iv == 0 else vp_int(Iv, p))
    V_total = sum(vs)
    W = S_target + V_total
    M = p ** W
    cmod = [[c % M for c in ci] for ci in coeffs]
    b = [1]
    ffs = [[f % M for f in ffall(rho, order)]]
    floor = W
    res_set = set(resonances_cert)
    for n in range(1, N + 1):
        rhs = 0
        for m in range(max(0, n - BW), n):
            bm = b[m]
            if bm == 0:
                continue
            dd = d0 + n - m
            w = 0
            fm = ffs[m]
            for i in range(order + 1):
                j = i + dd
                if 0 <= j < len(cmod[i]) and cmod[i][j]:
                    w += fm[i] * cmod[i][j]
            if w:
                rhs = (rhs + bm * w) % M
        Iv = ival(rho + n)
        if Iv == 0:
            if n not in res_set:
                raise RuntimeError(f"uncertified resonance at n={n}")
            if rhs % p ** min(floor, 40):
                raise RuntimeError(f"resonance rhs !=0 mod p^40 at n={n}")
            b.append(0)
        else:
            v = vs[n - 1]
            if v:
                if rhs % p ** v:
                    raise RuntimeError(f"P-DENOMINATOR at step n={n}")
                rhs //= p ** v
                floor -= v
            u = Iv // p ** v if v else Iv
            b.append(-rhs * pow(u % M, -1, M) % M)
        ffs.append([f % M for f in ffall(rho + n, order)])
    F = floor
    pf = p ** F
    return [x % pf for x in b], F, W, V_total

def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for ch in iter(lambda: f.read(1 << 20), b""):
            h.update(ch)
    return h.hexdigest()

def sha256_self():
    return sha256_file(__file__)

def lint_check(path):
    """Optional producer-lint hook: set FROB_PRODUCER_LINT to a lint script;
    skipped (returns a note string) when unset."""
    import os, subprocess, sys
    linter = os.environ.get("FROB_PRODUCER_LINT")
    if not linter or not os.path.exists(linter):
        return "producer_lint: SKIPPED (FROB_PRODUCER_LINT unset)"
    r = subprocess.run([sys.executable, linter, "--check", path],
                       capture_output=True, text=True)
    return (r.stdout.strip() or r.stderr.strip())
