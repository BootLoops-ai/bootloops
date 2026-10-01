"""abcount — S2 simple-genus-2-blocks extension (stage S2 of the sha-pinned
spec M3T_REGISTRATION.md sec 4): battery members B2, B3 — NON-DIAGONAL genus-2 tau,
cross-method vs hyperellcharpoly, Jacobi-sum desk truth for B3.

EXTENDS abcount_s0 (S0-banked, UNMODIFIED) and abcount_s1 (S1-banked,
UNMODIFIED); the Q6 single conversion site REMAINS
abcount_s0.convert_sign_convention — this module introduces NO second site.
`hecke_desk` (Builds item 3a) lives in hecke_desk.py.

New in S2 (all seed-committed in battery/SEEDS.json S2_* blocks):
  * member instantiation of the analytic Jacobian for NON-real-ordered branch
    points (complex-chain chord integrals through the NEW bridge
    period_g2_bridge.jl, which assembles the Eichler wrapper's exposed-but-
    unassembled primitives; the Eichler files themselves stay untouched);
  * homology recovery: LLL on the antisymmetric Riemann-relation kernel,
    Pfaffian gate, exact integer symplectic reduction, tau = A^{-1}B with
    certified symmetry + positivity (SEEDS homology_recovery_rule);
  * desk genus-2 theta box partial sums for NON-diagonal tau + the C7 tail
    bound (reused verbatim from abcount_s1.tail_bound_c7 — the inequality is
    genus-general) — dominance gate per SEEDS S2_theta_receipts;
  * degree-4 Weil box (manuals/abcount.md sec 2b) + quartic-from-counts
    assembly + deg-8 squaring assembly for the B2 product member;
  * Moebius branch-matching tie receipt (SEEDS S2_tie_receipt).
"""

import math
import os
import subprocess
import time
from fractions import Fraction

from flint import acb, arb, ctx, fmpz_mat

import abcount_s0 as ab
import abcount_s1 as s1

PERIOD_BRIDGE_JL = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "period_g2_bridge.jl")

J_STD_4 = [[0, 0, 1, 0], [0, 0, 0, 1], [-1, 0, 0, 0], [0, -1, 0, 0]]


def normalize_dyadic(man, e):
    """Strip trailing zero bits so python-side dyadics match the bridge's
    arf_dump normalization (exact value unchanged; enables raw echo compare)."""
    man = int(man)
    e = int(e)
    if man == 0:
        return 0, 0
    while man % 2 == 0:
        man //= 2
        e += 1
    return man, e


# ----------------------------------------------------------------------------
# bridge I/O (exact dyadic both ways — the S1 bridge discipline, new modes)
# ----------------------------------------------------------------------------

def _parse_ball(w, off):
    """8 ints at w[off:off+8] -> acb ball (mid man/exp pairs + rad man/exp)."""
    re_m = s1._arb_exact_from_man_exp(w[off], w[off + 1])
    im_m = s1._arb_exact_from_man_exp(w[off + 2], w[off + 3])
    re_r = s1._arb_exact_from_man_exp(w[off + 4], w[off + 5])
    im_r = s1._arb_exact_from_man_exp(w[off + 6], w[off + 7])
    return acb(re_m + arb("0 +/- 1") * re_r, im_m + arb("0 +/- 1") * im_r), re_r.max(im_r)


def bridge_call(mode_lines, workdir, tag, timeout=3300):
    inp = os.path.join(workdir, "s2_in_%s.txt" % tag)
    outp = os.path.join(workdir, "s2_out_%s.txt" % tag)
    with open(inp, "w") as f:
        f.write("\n".join(mode_lines) + "\n")
    t0 = time.time()
    s1.require_eichler()
    r = subprocess.run(
        [s1.JULIA, "--project=%s" % s1.EICHLER_PROJECT, PERIOD_BRIDGE_JL, inp, outp],
        capture_output=True, text=True, timeout=timeout,
        env=dict(os.environ, ABACUS_SIEGEL_JL=s1.SIEGEL_JL))
    wall = time.time() - t0
    if r.returncode != 0:
        raise RuntimeError("FAIL-INTERNAL: period bridge rc=%d stderr=%s"
                           % (r.returncode, r.stderr[-2000:]))
    out = {"lines": [], "subprocess_wall_seconds": wall}
    with open(outp) as f:
        for line in f:
            w = line.split()
            if w:
                out["lines"].append(w)
    return out


def periods_from_bridge(coeffs, chain, prec, workdir, tag, budget_ns=0):
    """MODE periods: certified chord integrals for the four chain cycles.
    Timed pilot: the bridge times ONE pilot integral and STOP-WALLs before the
    sweep if 15x the pilot exceeds budget_ns.  Returns (Pi, roots, receipt)."""
    lines = ["mode periods", "prec %d" % prec,
             "coeffs " + " ".join(str(c) for c in coeffs),
             "chain " + " ".join(str(c) for c in chain)]
    if budget_ns:
        lines.append("budget_ns %d" % budget_ns)
    res = bridge_call(lines, workdir, tag)
    for w in res["lines"]:
        if w[0] == "stop_wall":
            raise RuntimeError("STOP-WALL: period pilot pricing breach: %s"
                               % " ".join(w[1:]))
    Pi = [[None] * 4 for _ in range(2)]
    roots = [None] * 5
    pilot_seconds = sweep_seconds = load_seconds = maxrss = None
    max_rad = arb(0)
    for w in res["lines"]:
        if w[0] == "pi":
            k, j = int(w[1]) - 1, int(w[2]) - 1
            ball, rad = _parse_ball(w, 3)
            Pi[k][j] = ball
            max_rad = max_rad.max(rad)
        elif w[0] == "root":
            roots[int(w[1]) - 1], _ = _parse_ball(w, 2)
        elif w[0] == "pilot_integral_ns":
            pilot_seconds = int(w[1]) / 1e9
        elif w[0] == "sweep_ns":
            sweep_seconds = int(w[1]) / 1e9
        elif w[0] == "load_ns":
            load_seconds = int(w[1]) / 1e9
        elif w[0] == "maxrss_bytes":
            maxrss = int(w[1]) / (1024 ** 3)
    if any(v is None for row in Pi for v in row) or any(r is None for r in roots):
        raise RuntimeError("FAIL-INTERNAL: periods bridge output incomplete")
    receipt = {"prec_bits": prec, "max_period_radius": max_rad.str(6),
               "pilot_integral_seconds": pilot_seconds,
               "sweep_seconds": sweep_seconds,
               "julia_load_seconds": load_seconds, "maxrss_gb": maxrss,
               "subprocess_wall_seconds": res["subprocess_wall_seconds"],
               "scheme": "SEEDS S2_period_construction.segment_integral_scheme (certified; analytic checker + half-plane certificates)"}
    return Pi, roots, receipt


def theta_from_bridge(g, tau_entries, prec, workdir, tag):
    """MODE theta: full symmetric tau (upper-triangle dyadic entries) ->
    all 2^(2g) theta values. tau_entries: dict (i,j)->(man_re, exp_re,
    man_im, exp_im) 0-indexed upper triangle."""
    lines = ["mode theta", "prec %d" % prec, "g %d" % g]
    for (i, j), (rm, re_, im, ime) in sorted(tau_entries.items()):
        lines.append("entry %d %d %d %d %d %d" % (i + 1, j + 1, rm, re_, im, ime))
    res = bridge_call(lines, workdir, tag)
    n = 1 << (2 * g)
    values, claimed = [None] * n, [None] * n
    echo = {}
    call_seconds = maxrss = None
    for w in res["lines"]:
        if w[0] == "theta":
            idx = int(w[1])
            values[idx], claimed[idx] = _parse_ball(w, 2)
        elif w[0] == "tau_echo":
            i, j = int(w[1]) - 1, int(w[2]) - 1
            echo[(i, j)] = (int(w[3]), int(w[4]), int(w[5]), int(w[6]))
        elif w[0] == "call_ns":
            call_seconds = int(w[1]) / 1e9
        elif w[0] == "maxrss_bytes":
            maxrss = int(w[1]) / (1024 ** 3)
    if any(v is None for v in values):
        raise RuntimeError("FAIL-INTERNAL: theta bridge output incomplete")
    echo_ok = all(echo.get((i, j)) == (rm, re_, im, ime)
                  for (i, j), (rm, re_, im, ime) in tau_entries.items())
    return {"values": values, "claimed_radii": claimed, "echo_exact": echo_ok,
            "call_seconds": call_seconds, "maxrss_gb": maxrss,
            "subprocess_wall_seconds": res["subprocess_wall_seconds"], "prec": prec}


def thomae_roots_from_bridge(tau, prec, workdir, tag):
    """MODE thomae: thomae_sextic(tau) -> chi5 nonzero cert + 6 certified
    roots + igusa_clebsch covariants of the Thomae sextic.

    Takes the FULL-precision member tau BALLS (mid + radius) — NOT the 64-bit
    dyadic truncation: the tie receipt compares branch cross-ratios at the
    certified-period accuracy, and a 2^-64 truncation would shift them far
    above the ball radii (the dyadic rule is a THETA-input rule only)."""
    lines = ["mode thomae", "prec %d" % prec]
    for (i, j) in ((0, 0), (0, 1), (1, 1)):
        z = tau[i][j]
        rm, re_ = s1._man_exp_of_mid(z.real)
        im, ime = s1._man_exp_of_mid(z.imag)
        rad = z.real.rad().abs_upper().max(z.imag.rad().abs_upper())
        rman, rexp = (rad.man_exp() if not rad.is_zero() else (0, 0))
        lines.append("entry %d %d %d %d %d %d %d %d"
                     % (i + 1, j + 1, rm, re_, im, ime, int(rman), int(rexp)))
    res = bridge_call(lines, workdir, tag)
    roots = [None] * 6
    cov = [None] * 4
    chi5_nonzero = None
    for w in res["lines"]:
        if w[0] == "troot":
            roots[int(w[1]) - 1], _ = _parse_ball(w, 2)
        elif w[0] == "cov":
            cov[int(w[1]) - 1], _ = _parse_ball(w, 2)
        elif w[0] == "chi5_nonzero":
            chi5_nonzero = (w[1] == "true")
    if any(r is None for r in roots) or chi5_nonzero is None:
        raise RuntimeError("FAIL-INTERNAL: thomae bridge output incomplete")
    return roots, cov, chi5_nonzero


def curve_covariants_from_bridge(coeffs, prec, workdir, tag):
    """MODE covariants: igusa_clebsch of the input curve (deg-5 model
    sextic-normalized inside the unmodified Eichler
    igusa_clebsch(HyperellipticCurve))."""
    lines = ["mode covariants", "prec %d" % prec,
             "coeffs " + " ".join(str(c) for c in coeffs)]
    res = bridge_call(lines, workdir, tag)
    cov = [None] * 4
    for w in res["lines"]:
        if w[0] == "cov":
            cov[int(w[1]) - 1], _ = _parse_ball(w, 2)
    if any(c is None for c in cov):
        raise RuntimeError("FAIL-INTERNAL: covariants bridge output incomplete")
    return cov


# ----------------------------------------------------------------------------
# homology recovery (SEEDS S2_period_construction.homology_recovery_rule)
# ----------------------------------------------------------------------------

PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def riemann_relation_vectors(Pi, scale_exp=80):
    """LLL relation vectors for (Pi N Pi^T)_{12} = 0 over the 6 antisym unknowns.

    Local precision boost (S1 'desk sums boost locally' pattern): the period
    balls carry ~1e-100 mids, and the products must be formed at matching
    working precision or the 10^scale_exp-scaled residual columns drown in
    arithmetic rounding (member arithmetic runs at dps 60)."""
    old_dps = ctx.dps
    ctx.dps = 130
    try:
        cvec = [Pi[0][i] * Pi[1][j] - Pi[0][j] * Pi[1][i] for (i, j) in PAIRS]
    finally:
        ctx.dps = old_dps
    K = 10 ** scale_exp
    rows = []
    for t in range(6):
        row = [1 if s == t else 0 for s in range(6)]
        for part in (cvec[t].real, cvec[t].imag):
            man, e = s1._man_exp_of_mid(part)
            val = Fraction(man) * (Fraction(2) ** e) * K
            row.append(int(val))
        rows.append(row)
    L = fmpz_mat(rows).lll()
    sols = []
    for r in range(6):
        v = [int(L[r, c]) for c in range(8)]
        if abs(v[6]) < 10 ** 20 and abs(v[7]) < 10 ** 20 and any(v[:6]):
            sols.append(v[:6])
    return sols


def antisym(nv):
    n12, n13, n14, n23, n24, n34 = nv
    return [[0, n12, n13, n14], [-n12, 0, n23, n24],
            [-n13, -n23, 0, n34], [-n14, -n24, -n34, 0]]


def pfaffian(A):
    return A[0][1] * A[2][3] - A[0][2] * A[1][3] + A[0][3] * A[1][2]


def inv4_unimodular(A):
    n = 4
    M = [[Fraction(A[i][j]) for j in range(n)] for i in range(n)]
    I = [[Fraction(int(i == j)) for j in range(n)] for i in range(n)]
    for col in range(n):
        piv = next(r for r in range(col, n) if M[r][col] != 0)
        M[col], M[piv] = M[piv], M[col]
        I[col], I[piv] = I[piv], I[col]
        d = M[col][col]
        M[col] = [x / d for x in M[col]]
        I[col] = [x / d for x in I[col]]
        for r in range(n):
            if r != col and M[r][col] != 0:
                f = M[r][col]
                M[r] = [a - f * b for a, b in zip(M[r], M[col])]
                I[r] = [a - f * b for a, b in zip(I[r], I[col])]
    assert all(x.denominator == 1 for row in I for x in row), "not unimodular"
    return [[int(x) for x in row] for row in I]


def _ext_gcd(a, b):
    old_r, r = a, b
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r != 0:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t
    return old_s, old_t, old_r


def symplectic_reduce(W):
    """Exact integer symplectic Gram-Schmidt: unimodular antisymmetric 4x4 W
    -> S with S W S^T = J_STD_4 (rows of S = new basis in old coordinates).
    Verified by the caller via the identity before use."""
    n = 4

    def form(u, v):
        return sum(u[i] * W[i][j] * v[j] for i in range(n) for j in range(n))

    avail = [[int(i == j) for j in range(n)] for i in range(n)]
    pairs = []
    for _ in range(2):
        e = f = None
        for i in range(len(avail)):
            for j in range(len(avail)):
                if i != j and abs(form(avail[i], avail[j])) == 1:
                    e, f = avail[i], avail[j]
                    if form(e, f) == -1:
                        e, f = f, e
                    break
            if e is not None:
                break
        if e is None:
            for i in range(len(avail)):
                vals = [form(avail[i], avail[j]) for j in range(len(avail))]
                g = 0
                for v in vals:
                    g = math.gcd(g, abs(v))
                if g == 1:
                    e = avail[i]
                    cur_g, cur_vec = 0, [0] * n
                    for j, v in enumerate(vals):
                        if v == 0:
                            continue
                        if cur_g == 0:
                            cur_g, cur_vec = v, list(avail[j])
                        else:
                            a, b, g2 = _ext_gcd(cur_g, v)
                            cur_vec = [a * x + b * y for x, y in zip(cur_vec, avail[j])]
                            cur_g = g2
                        if abs(cur_g) == 1:
                            break
                    assert abs(cur_g) == 1
                    f = cur_vec if cur_g == 1 else [-x for x in cur_vec]
                    break
            assert e is not None, "no unimodular symplectic pair (W not unimodular?)"
        pairs.append((list(e), list(f)))
        new_avail = []
        for v in avail:
            if v is e or v is f:
                continue
            a = form(v, f)
            b = form(e, v)
            vv = [x - a * y - b * z for x, y, z in zip(v, e, f)]
            if any(vv):
                new_avail.append(vv)
        avail = new_avail
    (e1, f1), (e2, f2) = pairs
    return [e1, e2, f1, f2]


def tau_from_periods_g2(Pi, nv, sign):
    """Candidate acceptance step (iv) of the SEEDS homology_recovery_rule.
    Returns (tau 2x2 acb or None, receipt)."""
    N = antisym([sign * x for x in nv])
    pf = pfaffian(N)
    receipt = {"candidate_nv": nv, "orientation_sign": sign, "pfaffian": pf}
    if abs(pf) != 1:
        receipt["reject"] = "Pfaffian not +-1"
        return None, receipt
    resid = acb(0)
    for a in range(4):
        for b in range(4):
            if N[a][b]:
                resid += N[a][b] * Pi[0][a] * Pi[1][b]
    if not bool(resid.contains(acb(0))):
        receipt["reject"] = "certified Riemann-relation ball check failed"
        return None, receipt
    receipt["riemann_relation_ball_contains_0"] = True
    W = inv4_unimodular(N)
    S = symplectic_reduce(W)
    chk = [[sum(S[i][a] * W[a][b] * S[j][b] for a in range(4) for b in range(4))
            for j in range(4)] for i in range(4)]
    if chk != J_STD_4:
        raise RuntimeError("FAIL-INTERNAL: symplectic reduction identity failed")
    receipt["symplectic_identity_SWS^T=J"] = True
    Pip = [[sum(Pi[k][j] * S[i][j] for j in range(4)) for i in range(4)]
           for k in range(2)]
    A = [[Pip[0][0], Pip[0][1]], [Pip[1][0], Pip[1][1]]]
    B = [[Pip[0][2], Pip[0][3]], [Pip[1][2], Pip[1][3]]]
    detA = A[0][0] * A[1][1] - A[0][1] * A[1][0]
    if bool(detA.contains(acb(0))):
        receipt["reject"] = "det A contains 0 (cannot certify inversion)"
        return None, receipt
    Ainv = [[A[1][1] / detA, -A[0][1] / detA], [-A[1][0] / detA, A[0][0] / detA]]
    tau = [[sum(Ainv[i][k] * B[k][j] for k in range(2)) for j in range(2)]
           for i in range(2)]
    if not bool((tau[0][1] - tau[1][0]).contains(acb(0))):
        receipt["reject"] = "tau symmetry ball check failed"
        return None, receipt
    receipt["tau_symmetry_residual_rad"] = (tau[0][1] - tau[1][0]).abs_upper().str(6)
    off = tau[0][1].union(tau[1][0])     # ball intersection is ideal; union is safe
    tau[0][1] = off
    tau[1][0] = off
    im11 = tau[0][0].imag
    det_im = tau[0][0].imag * tau[1][1].imag - tau[0][1].imag * tau[1][0].imag
    if not (bool(im11 > 0) and bool(det_im > 0)):
        receipt["reject"] = "Im tau not certified positive-definite"
        return None, receipt
    receipt["im_tau_positive_definite_certified"] = True
    receipt["symplectic_S_rows"] = S
    return tau, receipt


def recover_tau(Pi):
    """Full SEEDS homology_recovery_rule: LLL -> deterministic candidate order
    -> first accepted candidate. Returns (tau, receipt).

    Runs under a local precision boost (S1 pattern) so the certified tau balls
    inherit the period accuracy rather than dps-60 arithmetic rounding."""
    old_dps = ctx.dps
    ctx.dps = 130
    try:
        return _recover_tau_inner(Pi)
    finally:
        ctx.dps = old_dps


def _recover_tau_inner(Pi):
    sols = riemann_relation_vectors(Pi)
    trace = {"relation_vectors": sols, "candidates_tried": []}
    if not len(sols):
        raise RuntimeError("FAIL-PERIOD-HOMOLOGY: no Riemann relation vectors found")
    combos = []
    if len(sols) == 1:
        combos = [sols[0]]
    else:
        v1, v2 = sols[0], sols[1]
        for x in range(-3, 4):
            for y in range(-3, 4):
                if (x, y) != (0, 0):
                    combos.append([x * a + y * b for a, b in zip(v1, v2)])
    for nv in combos:
        for sign in (1, -1):
            tau, rec = tau_from_periods_g2(Pi, nv, sign)
            trace["candidates_tried"].append(rec)
            if tau is not None:
                trace["accepted"] = {"nv": nv, "sign": sign}
                return tau, trace
    raise RuntimeError("FAIL-PERIOD-HOMOLOGY: no candidate produced a certified principal tau")


# ----------------------------------------------------------------------------
# desk genus-2 theta box sums (non-diagonal tau) — C7 desk side at genus 2
# ----------------------------------------------------------------------------

def desk_theta_g2_all(tau, N):
    """Direct 2-D box partial sums S_N for all 16 characteristics of the
    (possibly non-diagonal) 2x2 tau: sum over ||n||_inf <= N of
    exp(i pi ((n+a/2)^T tau (n+a/2) + (n+a/2).b)).  Ball arithmetic;
    independent of acb_theta's code path (substrate law)."""
    ipi = acb(0, arb.pi())
    out = []
    for idx in range(16):
        a1, a2 = (idx >> 3) & 1, (idx >> 2) & 1
        b1, b2 = (idx >> 1) & 1, idx & 1
        s = acb(0)
        for n1 in range(-N, N + 1):
            h1 = acb(2 * n1 + a1) / 2
            for n2 in range(-N, N + 1):
                h2 = acb(2 * n2 + a2) / 2
                q = (h1 * h1 * tau[0][0] + 2 * h1 * h2 * tau[0][1]
                     + h2 * h2 * tau[1][1] + h1 * b1 + h2 * b2)
                s += (ipi * q).exp()
        out.append(s)
    return out


def lambda_min_lower_2x2(im_tau):
    """Certified lower bound for the smallest eigenvalue of the symmetric 2x2
    real ball matrix Im tau: (tr - sqrt(tr^2 - 4 det))/2 evaluated in balls."""
    t = im_tau[0][0] + im_tau[1][1]
    d = im_tau[0][0] * im_tau[1][1] - im_tau[0][1] * im_tau[1][0]
    disc = t * t - 4 * d
    if not bool(disc >= 0):
        disc = disc.union(arb(0))
    lam = (t - disc.sqrt()) / 2
    return lam.lower()


# ----------------------------------------------------------------------------
# integer routes: quartic from counts, deg-8 squaring, deg-4 Weil box
# ----------------------------------------------------------------------------

def quartic_from_counts(n1, n2, p):
    """SEEDS derivation: charpoly [1, -e1, e2, -p*e1, p^2] (descending) from
    #C(F_p) = n1, #C(F_p^2) = n2."""
    t1 = p + 1 - n1
    t2 = p * p + 1 - n2
    if (t1 * t1 - t2) % 2 != 0:
        raise ValueError("FAIL-COUNT-PARITY: t1^2 - t2 odd")
    e1 = t1
    e2 = (t1 * t1 - t2) // 2
    return [1, -e1, e2, -p * e1, p * p]


def count_curve_Fp(fcoeffs, p):
    """#C(F_p) for y^2 = f(x), deg f = 5 (one point at infinity); pure python."""
    n = 1
    for x in range(p):
        v = 0
        for c in reversed(fcoeffs):
            v = (v * x + c) % p
        if v == 0:
            n += 1
        elif pow(v, (p - 1) // 2, p) == 1:
            n += 2
    return n


def count_curve_Fp2(fcoeffs, p):
    """#C(F_{p^2}); F_{p^2} = F_p[s]/(s^2 - ns), ns the least nonresidue."""
    ns = next(n for n in range(2, p) if pow(n, (p - 1) // 2, p) == p - 1)
    q = p * p

    def mul(a, b):
        return ((a[0] * b[0] + ns * a[1] * b[1]) % p,
                (a[0] * b[1] + a[1] * b[0]) % p)

    def powe(a, e):
        r = (1, 0)
        while e:
            if e & 1:
                r = mul(r, a)
            a = mul(a, a)
            e >>= 1
        return r

    n = 1
    for a0 in range(p):
        for a1 in range(p):
            x = (a0, a1)
            v = (0, 0)
            for c in reversed(fcoeffs):
                v = mul(v, x)
                v = ((v[0] + c) % p, v[1])
            if v == (0, 0):
                n += 1
            elif powe(v, (q - 1) // 2) == (1, 0):
                n += 2
    return n


def hyperellcharpoly_pari(fcoeffs, p):
    """Path 1 planted truth: PARI hyperellcharpoly (padic_counters row).
    Returns monic descending coefficient list [1, c3, c2, c1, c0]."""
    import cypari2
    pari = cypari2.Pari()
    x = pari("x")
    f = sum(int(c) * x ** i for i, c in enumerate(fcoeffs))
    cp = pari.hyperellcharpoly(pari.Mod(f, p))
    return [int(cp.polcoef(k)) for k in range(4, -1, -1)]


def square_quartic_deg8(q4):
    """Exact integer squaring of the monic quartic (descending) -> the deg-8
    charpoly of the product member J(C)^2, ASCENDING c_0..c_8 with c_0 = 1
    monic-reversed to match abcount_s1's assembly convention (c[k] = coefficient
    of T^(8-k) of the monic deg-8 charpoly)."""
    a = list(reversed(q4))       # ascending quartic
    c = [0] * 9
    for i in range(5):
        for j in range(5):
            c[i + j] += a[i] * a[j]
    c = list(reversed(c))        # descending deg-8: c[0] = 1 leading
    assert c[0] == 1
    return c                     # c[k] = coeff of T^(8-k); matches s1 layer


def functional_equation_receipt_deg4(cp, p):
    """manuals sec 2b: quartic a_(4-k) = p^(2-k) a_k, k = 0, 1 (cp descending:
    [1, a1, a2, a3, a4])."""
    a = cp
    rows = {"k=0": a[4] == p ** 2 * a[0], "k=1": a[3] == p * a[1]}
    return all(rows.values()), rows


def weil_box_isolate_deg4(coeff_enclosures, p):
    """Degree-4 Weil box (manuals sec 2b): FULL ENUMERATION of integer pairs
    (a_1, a_2) with |a_k| <= C(4,k) p^(k/2) intersected with the enclosures;
    a_3 = p a_1 and a_4 = p^2 are then forced by the functional equation and
    verified against their enclosures.  Same isolation semantics as deg-8 (C3):
    never per-coefficient rounding."""
    per_k = []
    for k in (1, 2):
        A = math.isqrt(math.comb(4, k) ** 2 * p ** k)
        ball = coeff_enclosures[k - 1]
        cands = [a for a in range(-A, A + 1) if ball.contains(acb(a))]
        per_k.append({"k": k, "box_halfwidth": A, "n_candidates": len(cands),
                      "cands": cands})
    tuples = [(a1, a2) for a1 in per_k[0]["cands"] for a2 in per_k[1]["cands"]
              if coeff_enclosures[2].contains(acb(p * a1))
              and coeff_enclosures[3].contains(acb(p * p))]
    receipt = {"per_k": [{kk: vv for kk, vv in d.items() if kk != "cands"}
                         for d in per_k],
               "n_tuples_enumerated": len(tuples),
               "semantics": "FULL enumeration box-intersect-enclosure; functional equation forces a3, a4; never rounding"}
    if len(tuples) == 1:
        return "OK", tuples, receipt
    if len(tuples) == 0:
        return "FAIL-WEIL-BOX-EMPTY", tuples, receipt
    return "UNDECIDED-PRECISION", tuples, receipt


# ----------------------------------------------------------------------------
# Moebius branch-matching tie receipt (SEEDS S2_tie_receipt)
# ----------------------------------------------------------------------------

def moebius_branch_match(curve_roots, thomae_roots):
    """curve side: 5 finite certified roots + infinity; canonical normalization
    M(z) = (z - r1)/(r2 - r1) maps r1 -> 0, r2 -> 1, infinity -> infinity;
    targets = images of r3, r4, r5.  Thomae side: all 120 ordered triples
    (i, j, k) of the 6 roots; N(z) = ((z - s_i)(s_j - s_k)) / ((z - s_k)(s_j - s_i));
    remaining 3 images must ball-match the 3 targets as a multiset."""
    r = curve_roots
    targets = [(r[m] - r[0]) / (r[1] - r[0]) for m in (2, 3, 4)]
    s = thomae_roots
    matches = []
    for i in range(6):
        for j in range(6):
            if j == i:
                continue
            for k in range(6):
                if k == i or k == j:
                    continue
                rest = [m for m in range(6) if m not in (i, j, k)]
                denom_ji = s[j] - s[i]
                if bool(denom_ji.contains(acb(0))):
                    continue
                imgs = []
                ok_map = True
                for m in rest:
                    den = (s[m] - s[k]) * denom_ji
                    if bool(den.contains(acb(0))):
                        ok_map = False
                        break
                    imgs.append((s[m] - s[i]) * (s[j] - s[k]) / den)
                if not ok_map:
                    continue
                used = [False] * 3
                good = True
                for t in targets:
                    hit = next((q for q in range(3) if not used[q]
                                and bool((imgs[q] - t).contains(acb(0)))), None)
                    if hit is None:
                        good = False
                        break
                    used[hit] = True
                if good:
                    matches.append((i, j, k))
    return matches


def covariant_cross_products(cov_curve, cov_thomae):
    """Secondary tie: weight-balanced cross-products (SEEDS secondary_covariant_check).
    Returns list of rows {check, overlap, anchored}."""
    w20, w40, w60, w100 = cov_curve
    t20, t40, t60, t100 = cov_thomae
    rows = []
    for name, lhs, rhs, anchors in (
            ("C20^5*C100' = C20'^5*C100", w20 ** 5 * t100, t20 ** 5 * w100, (w20, t20, w100, t100)),
            ("C40*C20'^2 = C40'*C20^2", w40 * t20 ** 2, t40 * w20 ** 2, (w40, t40, w20, t20)),
            ("C60*C20'^3 = C60'*C20^3", w60 * t20 ** 3, t60 * w20 ** 3, (w60, t60, w20, t20))):
        anchored = all(not bool(a.contains(acb(0))) for a in anchors)
        rows.append({"check": name,
                     "overlap": bool((lhs - rhs).contains(acb(0))),
                     "anchored": anchored})
    return rows
