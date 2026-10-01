#!/usr/bin/env python3
r"""
frobenius.py — local Frobenius solution basis at a regular-singular point,
and Frobenius LANDING of a transported vector onto that basis.

DESIGN
======
* Per-node Frobenius landing at a regular-singular point,
  matched to the regular backward transport. The ancestor pattern built only
  the analytic (exponent-0, resonance-free) tower seeded in ker(Mres); this
  module builds the FULL fundamental matrix Y(u) = P(u) * u^M (matrix power
  via expm(M log u)), so equal-eigenvalue/Jordan resonances and their log
  layers are included.
* Why Frobenius landing exists at all: large production fixed-eps
  transports need a Frobenius landing at points where some rows of the
  system are singular; plain Taylor transport asserts out there.
* REAL-eps0-only policy: the upstream complex-eps Frobenius matching
  path carries a FLAGGED Im-linear defect — matching residuals
  grow linearly in Im(eps0), root cause not yet closed. Until that flag is
  lifted, frobenius_basis()/land() raise ValueError on complex eps0. Use
  real-eps0 nodes + epsfan.cauchy_laurent on a real segment instead... note:
  the eps-fan circle itself feeds transport_fixed_eps (which IS complex-eps
  safe); only the Frobenius MATCHING is restricted.
* Import-dps footgun: all
  string parsing inside mp.workdps; global mp.dps never touched.

MATH (matrix Frobenius, regular-singular point)
===============================================
Shift u = x - x_sing. Regular-singular means A(u) = M/u + B(u) with B
analytic at 0. Seek the fundamental matrix

    Y(u) = P(u) * exp(M log u),   P(u) = I + sum_{q>=1} P_q u^q .

Y' = AY with E' = (M/u)E gives u P' = [M, P] + u B P, i.e. the Sylvester
recursion (solved by n^2 x n^2 LU):

    q P_q - M P_q + P_q M = C_q,   C_q = sum_{k=1..q} B_{k-1} P_{q-k} .

The operator q - ad_M has spectrum {q - (lambda_i - lambda_j)}, so it is
invertible iff no two eigenvalues of M differ by the nonzero integer q. Equal eigenvalues (Jordan blocks) are fine — the log
layers live inside exp(M log u).

SHEARING (nonzero-integer eigenvalue gaps):
while some pair of residue eigenvalues differs by a nonzero integer, the
TOP cluster lambda* of the gapped congruence class is split off via its
resolvent-contour spectral projector

    Pi = (1/2 pi i) oint_{|z-lambda*|=sep/2} (z I - M)^{-1} dz

(trapezoid rule on the circle — geometric convergence; projector QUALITY IS
CHECKED: ||Pi^2-Pi||, ||[M,Pi]||, tr Pi = mult, else RuntimeError), a basis
T = [range(Pi) | range(I-Pi)] block-diagonalizes M (off-blocks checked ~0,
then zeroed), and the Moser/Turrittin gauge Y = T diag(u I_G, I) Z lowers
the lambda* eigengroup by exactly 1:

    M -> [[M_11 - I, Btil_12(0)], [0, M_22]],
    B -> [[Btil_11, (Btil_12(u)-Btil_12(0))/u], [u Btil_21, Btil_22]] .

Each shear consumes one Taylor order of the 12-block, so the Cauchy ring is
over-sampled upfront by the shear-count bound (sum over congruence classes
of their integer span — each shear lowers the class top by 1 until merge).
The chain repeats until gap-free, then the Sylvester recursion applies to
the reduced (M, B) and the returned basis is
Y(u) = [prod_i T_i D_i(u)] P(u) exp(M_red log u) (gap-free systems take the
IDENTICAL old path: no factors, byte-same numerics). Verified against
hand-derived 2x2 gap-1 closed forms (diagonal, resonant-log, conjugated
resonant-log) in tests/test_frobenius.py. Any shear quality-gate failure
raises RuntimeError naming the alternatives: wayfinder.frob_scalar for
SCALAR operators with resonant towers (validated at 260d two-precision),
the driver-side resonance-tolerant landing pattern (integer-gap
eigenvalue classes clustered, LS-matched), and
tools/frobenius_boundary (via wayfinder.boundary_branches) for branch
classification at var=infinity.
One large production case verified det(qI - Mres) != 0 for
q >= 1, so it never enters the shearing path.

M and the B_k are extracted from pointwise A via a Cauchy circle around
x_sing (shared machinery from transport.py), with a double-pole tripwire:
a non-negligible u^-2 Laurent coefficient => not regular-singular =>
ValueError. Matching in land() is least-squares over several matching
points (overdetermined), residual reported — no silent trust.
"""

from mpmath import mp, mpf, mpc

try:  # optional acceleration hook — NEVER required
    import flint as _flint  # noqa: F401
    HAVE_FLINT = True
except Exception:
    HAVE_FLINT = False

try:
    from .transport import (_to_mpc, _circle_tables, _sample_circle, _fft,
                            _next_pow2, _declared_sings, transport_fixed_eps)
except ImportError:  # script-style use without package __init__
    from transport import (_to_mpc, _circle_tables, _sample_circle, _fft,
                           _next_pow2, _declared_sings, transport_fixed_eps)

__all__ = ["frobenius_basis", "land", "HAVE_FLINT"]


def _require_real_eps(eps_v, dps, who):
    """REAL eps0 nodes only: complex-eps Frobenius matching has a flagged
    Im-linear defect upstream."""
    tol = mpf(10) ** (-max(10, dps // 2))
    if abs(eps_v.imag) > tol * max(mpf(1), abs(eps_v.real)):
        raise ValueError(
            f"{who}: complex eps0 = {mp.nstr(eps_v, 8)} refused by policy — "
            f"the complex-eps Frobenius matching path carries a FLAGGED "
            f"Im-linear defect upstream. Use real eps0 nodes; complex eps is fine "
            f"in transport_fixed_eps, only the Frobenius MATCHING is "
            f"restricted until the flag is lifted.")


def _laurent_ring(desys, center, eps_v, r, kmax, wp):
    """Laurent coefficients of A about center from one Cauchy circle:
    returns (Cm2, M, B) with Cm2 = u^-2 coefficient matrix (regular-singular
    tripwire), M = residue (u^-1), B = [B_0..B_{kmax}] Taylor part.
    All entries mpc at working dps wp (call inside mp.workdps(wp))."""
    n = desys.n
    N = _next_pow2(max(2 * (kmax + 4), int(3.5 * wp) + 16))
    _, wroots_inv = _circle_tables(N, wp)
    samples = _sample_circle(desys, center, eps_v, r, N, wp)
    rinv = 1 / mpf(r)
    Cm2 = [[mpc(0)] * n for _ in range(n)]
    M = [[mpc(0)] * n for _ in range(n)]
    B = [[[mpc(0)] * n for _ in range(n)] for _ in range(kmax + 1)]
    for i in range(n):
        for j in range(n):
            col = [mpc(samples[k][i][j]) for k in range(N)]
            if all(v == 0 for v in col):
                continue
            hat = _fft(col, wroots_inv)
            Cm2[i][j] = hat[(-2) % N] / N * r ** 2
            M[i][j] = hat[(-1) % N] / N * r
            rp = mpf(1)
            for m in range(kmax + 1):
                B[m][i][j] = hat[m] / N * rp
                rp *= rinv
    scale = max(max(abs(mpc(samples[k][i][j])) for i in range(n)
                    for j in range(n)) for k in range(N))
    return Cm2, M, B, scale


def _mat_mul(A, B):
    n = len(A)
    m = len(B[0])
    K = len(B)
    out = [[mpc(0)] * m for _ in range(n)]
    for i in range(n):
        Ai = A[i]
        for k in range(K):
            a = Ai[k]
            if a == 0:
                continue
            Bk = B[k]
            oi = out[i]
            for j in range(m):
                if Bk[j] != 0:
                    oi[j] += a * Bk[j]
    return out


def _solve_sylvester(q, M, C, lu_cache):
    """Solve q*P - M P + P M = C for P (n x n).

    n <= 6: original dense n^2 x n^2 LU path (byte-identical legacy behavior;
    verified against a component-wise ODE plug-in of the resonant 2x2 test).
    n > 6:  Bartels-Stewart via ONE cached complex Schur M = Q U Q^H
    (the dense path is O(n^6) bignum per Taylor order — measured
    ~40 min/order at n=27, wp 250; Schur sweep is O(n^3)/order, seconds):
        q P' - U P' + P' U = C',  P' = Q^H P Q,  C' = Q^H C Q
        (q - U_ii + U_jj) P'_ij = C'_ij + sum_{c>i} U_ic P'_cj
                                        - sum_{c<j} P'_ic U_cj
    swept j ascending / i descending. Post-shear M is gap-free, so the
    denominators are bounded away from 0 for all q >= 1; a loud tripwire
    guards the claim (never silent).
    lu_cache: dict; key q -> dense matrix (legacy), key '_schur' -> (Q, U)."""
    n = len(M)
    if n <= 6:
        if q not in lu_cache:
            big = mp.zeros(n * n)
            for a in range(n):        # row block index (row of P)
                for b in range(n):    # col block index (col of P)
                    r0 = b * n + a    # vec index of P[a][b]
                    big[r0, r0] += q
                    for c in range(n):
                        # -(M P)[a][b] = -sum_c M[a][c] P[c][b]
                        big[r0, b * n + c] -= M[a][c]
                        # +(P M)[a][b] = +sum_c P[a][c] M[c][b]
                        big[r0, c * n + a] += M[c][b]
            lu_cache[q] = big
        rhs = mp.matrix(n * n, 1)
        for a in range(n):
            for b in range(n):
                rhs[b * n + a] = C[a][b]
        sol = mp.lu_solve(lu_cache[q], rhs)
        P = [[mpc(0)] * n for _ in range(n)]
        for a in range(n):
            for b in range(n):
                P[a][b] = mpc(sol[b * n + a])
        return P
    # ---- Bartels-Stewart (n > 6) ----
    if '_schur' not in lu_cache:
        Q, U = mp.schur(mp.matrix(M))
        lu_cache['_schur'] = (Q, U)
    Q, U = lu_cache['_schur']
    Cm = Q.H * mp.matrix(C) * Q
    Pp = mp.zeros(n)
    dmin = None
    for j in range(n):
        for i in range(n - 1, -1, -1):
            acc = Cm[i, j]
            for c in range(i + 1, n):
                if U[i, c] != 0 and Pp[c, j] != 0:
                    acc += U[i, c] * Pp[c, j]
            for c in range(j):
                if Pp[i, c] != 0 and U[c, j] != 0:
                    acc -= Pp[i, c] * U[c, j]
            den = q - U[i, i] + U[j, j]
            ad = abs(den)
            if dmin is None or ad < dmin:
                dmin = ad
            if ad < mpf(10) ** (-(mp.dps // 3)):
                raise RuntimeError(
                    f"_solve_sylvester: near-singular Bartels-Stewart pivot "
                    f"|q - lam_i + lam_j| = {mp.nstr(ad, 3)} at q={q} — "
                    f"residue not gap-free (shearing incomplete?)")
            Pp[i, j] = acc / den
    Pm = Q * Pp * Q.H
    return [[mpc(Pm[i, j]) for j in range(n)] for i in range(n)]


def _log_layers(Mmat, eigs, wp, ctol=None):
    """Max log power in exp(M log u): for each eigenvalue cluster lambda,
    the longest Jordan chain length minus 1 = smallest k with
    rank((M - lambda I)^k) == n - alg_mult. Numerical rank via svd_c."""
    n = Mmat.rows
    tol_c = ctol if ctol is not None else mpf(10) ** (-(wp // 2))
    clusters = []
    for lam in eigs:
        for cl in clusters:
            if abs(lam - cl[0]) < tol_c * max(mpf(1), abs(cl[0])):
                cl[1] += 1
                break
        else:
            clusters.append([lam, 1])
    log_max = 0
    for lam, mult in clusters:
        if mult == 1:
            continue
        Bm = Mmat - lam * mp.eye(n)
        Pk = mp.eye(n)
        for k in range(1, mult + 1):
            Pk = Pk * Bm
            sv = mp.svd_c(Pk, compute_uv=False)
            smax = max(sv) if len(sv) else mpf(0)
            rk = sum(1 for s in sv if s > tol_c * max(smax, mpf(1)))
            if rk <= n - mult:
                log_max = max(log_max, k - 1)
                break
        else:
            log_max = max(log_max, mult - 1)
    return log_max, clusters


_SHEAR_ALTERNATIVES = (
    "Alternatives if shearing cannot proceed: (a) SCALAR operators / "
    "P-recurrences with resonant integer towers -> wayfinder.frob_scalar "
    "(per-root eps-jet Frobenius, validated at 260d two-precision); (b) a "
    "driver-side resonance-tolerant landing (integer-gap eigenvalue classes "
    "clustered, LS-matched); (c) branch classification at var=infinity -> "
    "tools/frobenius_boundary via wayfinder.boundary_branches.")


def _eig_clusters(eigs, wp, ctol=None):
    """Cluster numerically-equal eigenvalues: [[value, mult], ...]
    (same tolerance discipline as _log_layers). ctol overrides the default
    10^-(wp/2) RELATIVE tolerance — needed when Jordan pseudospectrum rings
    (radius ~10^-(wp/m) for depth-m chains, a measured footgun) sit AT the
    default tolerance and would shatter a defective cluster into singletons
    whose 'spectral projectors' then rightly fail the quality gates
    (measured on a production defective-cluster landing)."""
    tol_c = ctol if ctol is not None else mpf(10) ** (-(wp // 2))
    clusters = []
    for lam in eigs:
        for cl in clusters:
            if abs(lam - cl[0]) < tol_c * max(mpf(1), abs(cl[0])):
                cl[1] += 1
                break
        else:
            clusters.append([mpc(lam), 1])
    return clusters


def _integer_gaps(clusters, wp):
    """[(i_hi, i_lo, k)]: cluster pairs whose values differ by the nonzero
    POSITIVE integer k (within tol) — the Sylvester-singular resonances."""
    tol_int = mpf(10) ** (-(wp // 3))
    gaps = []
    for i, (li, _mi) in enumerate(clusters):
        for j, (lj, _mj) in enumerate(clusters):
            if i == j:
                continue
            d = li - lj
            k = int(mp.nint(d.real))
            if k >= 1 and abs(d - k) < tol_int:
                gaps.append((i, j, k))
    return gaps


def _shear_bound(clusters, wp):
    """Upper bound on the number of shears: sum over integer-gap-connected
    cluster components of the max integer gap (= the class span; each shear
    lowers the class top by exactly 1 until it merges downward)."""
    gaps = _integer_gaps(clusters, wp)
    if not gaps:
        return 0
    parent = list(range(len(clusters)))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for i, j, _k in gaps:
        parent[find(i)] = find(j)
    span = {}
    for i, _j, k in gaps:
        r = find(i)
        span[r] = max(span.get(r, 0), k)
    return sum(span.values())


def _spectral_projector_hermite(Mlst, clusters, ci, wp):
    """Pi_c = f_c(M) by Hermite-interpolation holomorphic calculus:
    f_c == 1 + O((z-lam_c)^{m_c}) at the target cluster, == O((z-lam_j)^{m_j})
    at every other cluster — a degree<n polynomial evaluated in M by Horner.
    Ported from a validated production spectral-projector implementation
    (validated to working precision incl. Jordan + shearing +
    tiny-exponent clusters).
    Cost ~n matmuls — replaces the ~wp*3.3-node resolvent contour (measured
    ~19 min/cluster at n=27, wp 220 -> seconds). Defective-cluster means
    carry only O(ring^2) pseudospectrum error (symmetric ring), far below
    wp; the shear quality gates remain the arbiter and the contour is the
    fallback."""
    n = len(Mlst)
    lam, mc = clusters[ci]

    def _shift(coef, x0):
        work = list(coef)
        out = []
        for _ in range(len(coef)):
            acc = mpc(0)
            for c in reversed(work):
                acc = acc * x0 + c
            out.append(acc)
            if len(work) <= 1:
                work = [mpc(0)]
            else:
                q = [mpc(0)] * (len(work) - 1)
                acc2 = mpc(0)
                for i in range(len(work) - 1, 0, -1):
                    acc2 = work[i] + acc2 * x0
                    q[i - 1] = acc2
                work = q
        return out

    with mp.workdps(wp):
        # g(z) = prod_{j != ci} (z - lam_j)^{m_j}
        gcoef = [mpc(1)]
        for j, (l2, m2) in enumerate(clusters):
            if j == ci:
                continue
            for _ in range(m2):
                new = [mpc(0)] * (len(gcoef) + 1)
                for k, c in enumerate(gcoef):
                    new[k + 1] += c
                    new[k] -= l2 * c
                gcoef = new
        gt = _shift(gcoef, lam)         # g about lam
        h = [mpc(0)] * mc               # Taylor inverse of g to order mc
        h[0] = 1 / gt[0]
        for q in range(1, mc):
            acc = mpc(0)
            for k in range(1, q + 1):
                gk = gt[k] if k < len(gt) else mpc(0)
                acc += gk * h[q - k]
            h[q] = -acc / gt[0]
        hz = _shift(h, -lam)            # back to z-coordinates
        # f = g * hz; Pi = f(M) by Horner
        fz = [mpc(0)] * (len(gcoef) + len(hz) - 1)
        for i, va in enumerate(gcoef):
            if va == 0:
                continue
            for j, vb in enumerate(hz):
                if vb:
                    fz[i + j] += va * vb
        res = [[mpc(0)] * n for _ in range(n)]
        for c in reversed(fz):
            res = _mat_mul(res, Mlst)
            for i in range(n):
                res[i][i] += c
        return res


def _spectral_projector(Mlst, center, r, wp):
    """Spectral projector Pi = (1/2 pi i) oint (zI-M)^{-1} dz on the circle
    |z-center| = r via the trapezoid rule (geometric convergence: nearest
    excluded eigenvalue sits at 2r by construction, ~log10(2) digits/node).
    Working precision is bumped by the resolvent leverage log10(1/r)."""
    n = len(Mlst)
    wp2 = wp + 15 + max(0, int(mp.ceil(-mp.log10(r))))
    N = int(mpf(wp) * mp.log(10) / mp.log(2)) + 32
    with mp.workdps(wp2):
        P = [[mpc(0)] * n for _ in range(n)]
        for kk in range(N):
            w = mpf(r) * mp.expjpi(mpf(2 * kk) / N)
            z = center + w
            Az = mp.matrix(n)
            for i in range(n):
                for j in range(n):
                    Az[i, j] = -mpc(Mlst[i][j])
                Az[i, i] += z
            R = mp.inverse(Az)
            for i in range(n):
                for j in range(n):
                    P[i][j] += w * R[i, j] / N
        return [[mpc(v) for v in row] for row in P]


def _shear_reduce(Mlst, B, kmax_valid, wp, ctol=None):
    """Moser/Turrittin shearing chain: repeat top-cluster split + gauge
    Y = T diag(u I_G, I) Z until the residue has no nonzero-integer
    eigenvalue gaps. Returns (M_fin, B_fin, factors) with factors =
    [(T, m_top), ...] in application order (Y = T_1 D_1 ... T_s D_s Z).
    All quality gates raise RuntimeError (with alternatives) on failure."""
    n = len(Mlst)
    tol_gate = mpf(10) ** (-(wp - 18))
    clusters0 = _eig_clusters(mp.eig(mp.matrix(Mlst), left=False, right=False),
                              wp, ctol)
    bound = _shear_bound(clusters0, wp)
    factors = []
    for _it in range(bound + n + 2):
        eigs = mp.eig(mp.matrix(Mlst), left=False, right=False)
        clusters = _eig_clusters(eigs, wp, ctol)
        gaps = _integer_gaps(clusters, wp)
        if not gaps:
            return Mlst, B, factors
        if kmax_valid < 1:
            raise RuntimeError(
                "shearing: ran out of valid Taylor orders (each shear "
                "consumes one 12-block order) — the upfront shear bound was "
                "violated; this is a bug, report it. " + _SHEAR_ALTERNATIVES)
        # top cluster of a gapped class: maximal 'hi' member (it cannot have
        # anything integer-above it, else that would be a gap pair too)
        hi_idx = max((i for i, _j, _k in gaps),
                     key=lambda i: (clusters[i][0].real, clusters[i][0].imag))
        lam, mtop = clusters[hi_idx]
        sep = min(abs(cl[0] - lam) for q, cl in enumerate(clusters)
                  if q != hi_idx)
        # quality gates: projector idempotent, commutes with M, rank = mult.
        # Hermite holomorphic-calculus projector FIRST (seconds); resolvent
        # contour (~19 min/cluster at n=27, wp 220) only as fallback. The
        # gates below arbitrate either way.
        Mm = mp.matrix(Mlst)
        nrm = lambda X: max(abs(X[i, j]) for i in range(n) for j in range(n))

        def _gate(Pm):
            q_idem = nrm(Pm * Pm - Pm)
            q_comm = nrm(Mm * Pm - Pm * Mm) / max(mpf(1), nrm(Mm))
            q_tr = abs(mp.fsum(Pm[i, i] for i in range(n)) - mtop)
            bad = (max(q_idem, q_comm) > tol_gate * 1000
                   or q_tr > mpf(10) ** (-(wp // 2)))
            return bad, q_idem, q_comm, q_tr

        Pi = _spectral_projector_hermite(Mlst, clusters, hi_idx, wp)
        Pm = mp.matrix(Pi)
        bad, q_idem, q_comm, q_tr = _gate(Pm)
        if bad:
            Pi = _spectral_projector(Mlst, lam, sep / 2, wp)
            Pm = mp.matrix(Pi)
            bad, q_idem, q_comm, q_tr = _gate(Pm)
        if bad:
            raise RuntimeError(
                f"shearing: spectral projector for cluster {mp.nstr(lam, 8)} "
                f"(mult {mtop}) failed quality gates (hermite AND contour): "
                f"||Pi^2-Pi||="
                f"{mp.nstr(q_idem, 3)}, ||[M,Pi]||/||M||={mp.nstr(q_comm, 3)}, "
                f"|tr Pi - m|={mp.nstr(q_tr, 3)} (wp={wp}). "
                + _SHEAR_ALTERNATIVES)
        # invariant bases: range(Pi) (m columns), range(I-Pi) (n-m columns)
        U1, S1, _V1 = mp.svd_c(Pm)
        rank1 = sum(1 for s in S1 if s > mpf(1) / 2)
        Qm = mp.eye(n) - Pm
        U2, S2, _V2 = mp.svd_c(Qm)
        rank2 = sum(1 for s in S2 if s > mpf(1) / 2)
        if rank1 != mtop or rank2 != n - mtop:
            raise RuntimeError(
                f"shearing: projector ranks ({rank1}, {rank2}) != "
                f"({mtop}, {n - mtop}) for cluster {mp.nstr(lam, 8)}. "
                + _SHEAR_ALTERNATIVES)
        T = mp.matrix(n)
        for i in range(n):
            for j in range(mtop):
                T[i, j] = U1[i, j]
            for j in range(n - mtop):
                T[i, mtop + j] = U2[i, j]
        sv = mp.svd_c(T, compute_uv=False)
        if min(sv) < mpf(10) ** (-(wp // 3)) * max(sv):
            raise RuntimeError(
                f"shearing: basis T ill-conditioned (smin/smax = "
                f"{mp.nstr(min(sv) / max(sv), 3)}) — skewed spectral split. "
                + _SHEAR_ALTERNATIVES)
        # Mt = T^-1 M T must be block-diagonal (both off-blocks ~ 0)
        Ti = mp.inverse(T)
        Mt = Ti * Mm * T
        off = max([abs(Mt[i, j]) for i in range(mtop) for j in range(mtop, n)]
                  + [abs(Mt[i, j]) for i in range(mtop, n) for j in range(mtop)]
                  + [mpf(0)])
        if off > tol_gate * 1000 * max(mpf(1), nrm(Mm)):
            raise RuntimeError(
                f"shearing: T^-1 M T off-block norm {mp.nstr(off, 3)} not "
                f"negligible — invariant-subspace split failed. "
                + _SHEAR_ALTERNATIVES)
        Bt = [Ti * mp.matrix(Bk) * T for Bk in B]
        # sheared residue and Taylor blocks (see module docstring MATH)
        Mn = [[mpc(0)] * n for _ in range(n)]
        for i in range(mtop):
            for j in range(mtop):
                Mn[i][j] = mpc(Mt[i, j])
            Mn[i][i] -= 1
        for i in range(mtop, n):
            for j in range(mtop, n):
                Mn[i][j] = mpc(Mt[i, j])
        for i in range(mtop):
            for j in range(mtop, n):
                Mn[i][j] = mpc(Bt[0][i, j])       # Btil_12(0) enters M
        Bn = []
        K = len(B) - 1
        for k in range(K + 1):
            Bk = [[mpc(0)] * n for _ in range(n)]
            for i in range(mtop):
                for j in range(mtop):
                    Bk[i][j] = mpc(Bt[k][i, j])
                for j in range(mtop, n):          # (Btil_12(u)-Btil_12(0))/u
                    Bk[i][j] = mpc(Bt[k + 1][i, j]) if k + 1 <= K else mpc(0)
            for i in range(mtop, n):
                for j in range(mtop):             # u * Btil_21(u)
                    Bk[i][j] = mpc(Bt[k - 1][i, j]) if k >= 1 else mpc(0)
                for j in range(mtop, n):
                    Bk[i][j] = mpc(Bt[k][i, j])
            Bn.append(Bk)
        Mlst, B = Mn, Bn
        factors.append(([[mpc(T[i, j]) for j in range(n)] for i in range(n)],
                        mtop))
        kmax_valid -= 1
    raise RuntimeError(
        "shearing: iteration bound exceeded without reaching a gap-free "
        "residue — non-terminating shear chain (should be impossible; "
        "report). " + _SHEAR_ALTERNATIVES)


def frobenius_basis(desys, eps0, x_sing, dps, kmax, direction=+1,
                    cluster_tol=None):
    """Frobenius fundamental basis Y(u) = P(u) exp(M log u) at a
    regular-singular point x_sing (u = x - x_sing), truncated at u^kmax.

    REAL eps0 only (policy — see the module docstring's Im-linear flag).

    Returns dict:
        n, var, x_sing, eps0, dps, kmax, direction
        indicial     eigenvalues of the ORIGINAL residue matrix M
                     (list[mpc], @dps)
        log_max      max power of log u appearing (Jordan resonance layers
                     of the REDUCED residue — the actual log structure)
        clusters     [(eigenvalue, algebraic multiplicity), ...] of the
                     reduced residue
        M            residue matrix of the (shear-)REDUCED system (== the
                     original when sheared is False) (list of list mpc, @dps)
        P            [P_0..P_kmax] Taylor matrices of the reduced P(u) (@dps)
        sheared      bool — True iff nonzero-integer indicial gaps were
                     removed by a Moser/Turrittin shearing chain
        n_shears     number of shears applied (0 when gap-free)
        shear_factors [(T, m_top), ...] @dps; the basis is
                     Y(u) = T_1 D_1(u) ... T_s D_s(u) P(u) exp(M log u),
                     D_i(u) = diag(u I_{m_top_i}, I) — empty when gap-free,
                     in which case Y(u) = P(u) exp(M log u) exactly as before
        indicial_reduced  eigenvalues of the reduced residue (u-powers in
                     exp(M log u); == indicial when not sheared)
        r_circle     Cauchy sampling radius used (mpf)
        p_last_norm  max-entry norm of P_kmax (truncation honesty handle:
                     series tail at |u| is roughly p_last_norm*|u|^kmax,
                     times the bounded shear-polynomial factors if sheared)
        eval(u, dps_out=None)    -> n x n fundamental matrix at u (complex
                     u allowed; principal branch of log u — approach side
                     for direction=-1 means negative real u, log picks +i*pi)
        eval_t(t, dps_out=None)  -> eval(direction * t) for real t > 0

    Raises:
        ValueError    complex eps0 (policy), or double pole at x_sing
                      (not regular-singular).
        RuntimeError  a shearing quality gate failed (projector idempotency/
                      commutation/rank, T conditioning, off-block residual)
                      — message names the frob_scalar / driver-side /
                      frobenius_boundary alternatives.
    """
    wp = dps + 30
    with mp.workdps(wp):
        eps_v = _to_mpc(eps0)
        _require_real_eps(eps_v, dps, "frobenius_basis")
        xs = _to_mpc(x_sing)
        sings = _declared_sings(desys, wp)
        others = [abs(s - xs) for s in sings
                  if abs(s - xs) > mpf(10) ** (-(wp // 2))]
        r_f = mpf(1) / 8
        if others:
            r_f = min(r_f, min(others) / 2)
        n = desys.n
        import time as _time
        timings = {}
        # SLACK Taylor orders sampled upfront: the circle
        # SAMPLING dominates (N is set by wp, not kmax) and shearing used to
        # force a full RE-sample; with slack >= shear_bound the one ring
        # serves both. Gap-free systems: byte-identical result, tiny extra
        # FFT post-processing only.
        _SLACK = 8
        _t0 = _time.time()
        Cm2, Mlst, B, scale = _laurent_ring(desys, xs, eps_v, r_f,
                                            kmax + _SLACK, wp)
        timings['ring_s'] = round(_time.time() - _t0, 1)
        # regular-singular tripwire: u^-2 coefficient must vanish
        m2 = max(abs(Cm2[i][j]) for i in range(n) for j in range(n))
        if m2 > scale * r_f ** 2 * mpf(10) ** (-(wp - 20)):
            raise ValueError(
                f"x_sing={mp.nstr(xs, 8)} is NOT regular-singular: u^-2 "
                f"Laurent coefficient of A has norm {mp.nstr(m2, 3)} "
                f"(circle scale {mp.nstr(scale, 3)}) — Frobenius basis "
                f"undefined; check the point / the DE loader")
        Mmat = mp.matrix(Mlst)
        eigs_orig = mp.eig(Mmat, left=False, right=False)
        # nonzero-integer indicial gaps => Moser/Turrittin shearing chain.
        # Each shear consumes one Taylor order of the 12-block; re-sample
        # ONLY if the bound exceeds the upfront slack.
        ctol = None if cluster_tol is None else mpf(str(cluster_tol))
        shear_bound = _shear_bound(_eig_clusters(eigs_orig, wp, ctol), wp)
        factors = []
        if shear_bound > 0:
            _t0 = _time.time()
            if shear_bound > _SLACK:
                _cm2, Mlst, B, _scale = _laurent_ring(desys, xs, eps_v, r_f,
                                                      kmax + shear_bound, wp)
            Mlst, B, factors = _shear_reduce(Mlst, B, kmax + shear_bound, wp,
                                             ctol)
            Mmat = mp.matrix(Mlst)
            timings['shear_s'] = round(_time.time() - _t0, 1)
        eigs = mp.eig(Mmat, left=False, right=False) if factors else eigs_orig
        log_max, clusters = _log_layers(Mmat, eigs, wp, ctol)
        # Sylvester recursion for P_q (residue is gap-free here)
        _t0 = _time.time()
        Plist = [[[mpc(1) if i == j else mpc(0) for j in range(n)]
                  for i in range(n)]]
        lu_cache = {}
        for q in range(1, kmax + 1):
            C = [[mpc(0)] * n for _ in range(n)]
            for k in range(1, q + 1):
                if k - 1 <= kmax:
                    T = _mat_mul(B[k - 1], Plist[q - k])
                    for i in range(n):
                        for j in range(n):
                            C[i][j] += T[i][j]
            Plist.append(_solve_sylvester(q, Mlst, C, lu_cache))
        timings['sylvester_s'] = round(_time.time() - _t0, 1)
        p_last_norm = max(abs(Plist[kmax][i][j])
                          for i in range(n) for j in range(n))
        # freeze copies for the closures (working precision, rounded on exit)
        M_wp = [[mpc(v) for v in row] for row in Mlst]
        P_wp = [[[mpc(v) for v in row] for row in Pq] for Pq in Plist]
        facs_wp = [([[mpc(v) for v in row] for row in T], m)
                   for T, m in factors]
        var = getattr(desys, "var", "x")

    def _eval(u, dps_out=None):
        d_out = dps_out or dps
        with mp.workdps(wp):
            uv = _to_mpc(u)
            if uv == 0:
                raise ValueError("eval at u=0: fundamental matrix is "
                                 "singular/log-divergent at the Frobenius "
                                 "point itself")
            E = mp.expm(mp.matrix(M_wp) * mp.log(uv))
            # P(u) via Horner on matrices
            acc = [[mpc(v) for v in row] for row in P_wp[kmax]]
            for q in range(kmax - 1, -1, -1):
                for i in range(n):
                    for j in range(n):
                        acc[i][j] = acc[i][j] * uv + P_wp[q][i][j]
            Y = [[sum(acc[i][k] * E[k, j] for k in range(n))
                  for j in range(n)] for i in range(n)]
            # shear factors (empty when gap-free): Y = T_1 D_1 ... T_s D_s Y
            for T, mtop in reversed(facs_wp):
                for i in range(mtop):          # D_i(u) = diag(u I_m, I)
                    for j in range(n):
                        Y[i][j] = Y[i][j] * uv
                Y = [[sum(T[i][k] * Y[k][j] for k in range(n))
                      for j in range(n)] for i in range(n)]
        with mp.workdps(d_out):
            return [[+v for v in row] for row in Y]

    def _eval_t(t, dps_out=None):
        with mp.workdps(wp):
            uv = _to_mpc(t) * direction
        return _eval(uv, dps_out)

    with mp.workdps(dps):
        out = {
            "n": n, "var": var,
            "x_sing": +xs, "eps0": +eps_v, "dps": dps, "kmax": kmax,
            "direction": direction,
            "indicial": [+mpc(e) for e in eigs_orig],
            "log_max": log_max,
            "clusters": [(+mpc(l), m) for l, m in clusters],
            "M": [[+v for v in row] for row in Mlst],
            "P": [[[+v for v in row] for row in Pq] for Pq in Plist],
            "sheared": bool(factors),
            "n_shears": len(factors),
            "shear_factors": [([[+v for v in row] for row in T], m)
                              for T, m in factors],
            "indicial_reduced": [+mpc(e) for e in eigs],
            "r_circle": +r_f,
            "p_last_norm": +p_last_norm,
            "timings": timings,
            "eval": _eval,
            "eval_t": _eval_t,
        }
    return out


def land(desys, eps0, x_from, y_from, x_sing, dps, kmax,
         allow_normal_eq=False, cluster_tol=None, match_dps=None,
         match_fracs=None, match_rings=None, match_phase_offsets=None,
         equilibrate_cols=True, cond_only=False):
    """Match a transported vector onto the Frobenius basis at x_sing.

    match_dps: precision for the internal matching-leg
    transports only (default None = dps, legacy). Set to the DATA's own
    accuracy when y_from is far less accurate than the tower precision —
    marching low-accuracy data at wp pays the quadratic bignum cost for
    zero information. The basis and the LS match stay at wp.

    match_rings / match_phase_offsets (the G-block conditioning
    fix): two-ring x multi-phase OVERDETERMINED match design.
    match_rings = list of |u| as FRACTIONS of the local Cauchy radius
    (e.g. [0.5, 0.125] — dynamic range 4, not the 5-point ray's 150x log
    spread); match_phase_offsets = phase offsets (radians) relative to the
    incoming direction (default [-pi/3, -pi/6, 0, +pi/6, +pi/3]; all
    points must stay off the log branch cut). Phase separates distinct
    lambda columns, the two moduli separate log towers. Overrides
    match_fracs when given. Transport route: down the ray to the outer
    ring, sweep phases, radial hop, sweep the inner ring (all chords stay
    >= cos(pi/4)*ring radius from the singularity).

    equilibrate_cols (default ON, disclosed tool change):
    scale each column of the matching matrix by its own max-abs before
    qr_solve and unscale kappa after — removes the dominant diagonal
    u^{Re lambda_i} scaling part of the measured cond 6.2e21/5.3e22
    (G-block ledger). Exact reparametrization: same LS problem, same
    residual; cond_est is reported for the EQUILIBRATED matrix.

    Pattern: least-squares matching of the regular backward transport onto
    the local Frobenius tower, upgraded to an
    OVERDETERMINED match: the vector y_from at x_from is transported
    (transport_fixed_eps, chained legs) to three matching points
    u_p = f_p * (x_from - x_sing) on the ray toward x_sing (f pinned so the
    points sit inside the Cauchy circle), and

        Y(u_p) kappa = y(x_p),   p = 1..3     (3n equations, n unknowns)

    is solved by least squares (mp.qr_solve).  The residual is REPORTED,
    never silently trusted.  If qr_solve fails, land() FAILS LOUDLY by
    default (F4 fix): the normal-equations fallback (A^H A)
    SQUARES the condition number of the matching matrix — whose columns
    scale as u^lambda_i and are genuinely ill-conditioned for spread
    indicial exponents — and a small residual does NOT bound the kappa
    error in near-degenerate directions.  Pass allow_normal_eq=True to
    accept the fallback anyway; it is then recorded as solver='normal_eq'.

    Returns dict:
        kappa          list[mpc] @dps — coefficients on the basis columns
        residual       absolute LS residual 2-norm
        residual_rel   residual / max(1, ||rhs||)
        solver         'qr' | 'normal_eq' (which solver produced kappa)
        cond_est       condition-number estimate smax/smin of the stacked
                       matching matrix (svd; str) — None if the svd failed
        match_points   the x-points used
        series_tail    p_last_norm * |u_1|^kmax truncation estimate
        basis          the frobenius_basis dict (with eval closures)
        dps, kmax

    Raises ValueError on complex eps0 (real-eps0 policy, via
    frobenius_basis), RuntimeError if qr_solve fails and
    allow_normal_eq=False, and everything transport_fixed_eps raises.
    """
    wp = dps + 30
    basis = frobenius_basis(desys, eps0, x_sing, dps, kmax,
                            cluster_tol=cluster_tol)  # eps policy here
    with mp.workdps(wp):
        eps_v = _to_mpc(eps0)
        xs = _to_mpc(x_sing)
        xf = _to_mpc(x_from)
        y = [_to_mpc(v) for v in y_from]
        L = abs(xf - xs)
        if L == 0:
            raise ValueError("x_from coincides with x_sing")
        r_f = mpf(basis["r_circle"])
        # matching radii: inside the sampling circle AND on the segment.
        # match_fracs: caller-supplied multipliers of
        # (x_from - x_sing). RATIONALE: the legacy 3 points span only a
        # factor 2 in |u| (log u nearly constant), so log^l towers and
        # u^{small-exponent} branch columns are close to collinear — the
        # matching cannot separate them (measured: cond 5.9e13, log-content
        # leak ~1e-10 on a production defective-cluster landing). A log-wide spread
        # (factor >=100 across the points, still inside the Cauchy disc and
        # the series-tail budget) restores the separation.
        match_design = "ray"
        if match_rings is not None:
            # two-ring x multi-phase design (G-block conditioning fix).
            match_design = "rings"
            offs = (match_phase_offsets if match_phase_offsets is not None
                    else [-mp.pi/3, -mp.pi/6, mpf(0), mp.pi/6, mp.pi/3])
            offs = [mpf(str(o)) if not isinstance(o, (type(mpf(0)),))
                    else o for o in offs]
            rims = sorted((mpf(str(R)) for R in match_rings), reverse=True)
            for R in rims:
                if R * r_f >= r_f:
                    raise ValueError(
                        f"match_rings: fraction {mp.nstr(R, 6)} not < 1 "
                        f"(points must sit inside the Cauchy radius)")
            diru = (xf - xs) / L
            base_arg = mp.arg(diru)
            # branch-cut guard: every point must keep its argument strictly
            # inside (-pi, pi) — the tower's log(u) is principal-branch.
            for R in rims:
                for o in offs:
                    a = base_arg + o
                    if not (-mp.pi + mpf('0.05') < a < mp.pi - mpf('0.05')):
                        raise ValueError(
                            f"match_phase_offsets: point argument "
                            f"{mp.nstr(a, 6)} too close to the log branch "
                            f"cut (incoming arg {mp.nstr(base_arg, 6)})")
            # transport route: enter outer ring at offset 0, sweep +side,
            # cross to -side (chords subtend <= pi/2 => distance from the
            # singularity >= cos(pi/4)*ring radius), radial hop at the last
            # phase, sweep the inner ring back.
            o_sorted = sorted(offs)
            i0 = min(range(len(o_sorted)), key=lambda i: abs(o_sorted[i]))
            route_offs = ([o_sorted[i0]] + o_sorted[i0+1:]
                          + list(reversed(o_sorted[:i0])))
            pts = []
            for ri, R in enumerate(rims):
                seq = route_offs if ri == 0 else list(reversed(route_offs))
                for o in seq:
                    pts.append(xs + R * r_f * diru * mp.exp(mpc(0, 1) * o))
            fracs = [abs(p - xs) / L for p in pts]   # for the tail estimate
        elif match_fracs is not None:
            fracs = [mpf(str(f)) for f in match_fracs]
            for f in fracs:
                if f * L >= r_f:
                    raise ValueError(
                        f"match_fracs: point |u|={mp.nstr(f*L, 6)} outside "
                        f"the Cauchy radius {mp.nstr(r_f, 6)}")
            pts = [xs + f * (xf - xs) for f in fracs]
        else:
            u1 = min(L / 16, r_f / 2)
            fracs = [u1 / L, u1 / L * mpf(2) / 3, u1 / L / 2]
            pts = [xs + f * (xf - xs) for f in fracs]
        n = desys.n
        # chain transports x_from -> p1 -> p2 -> ... (route sing-free).
        # cond_only (pilot mode): skip the transports and the
        # solve — build the EXACT production matching matrix (basis evals
        # at the design points + equilibration) and report its cond only.
        ys = []
        if not cond_only:
            cur_x, cur_y = xf, [mpc(v) for v in y]
            t_dps = match_dps if match_dps is not None else wp
            for p in pts:
                cur_y = transport_fixed_eps(desys, eps_v, cur_x, p,
                                            cur_y, t_dps)
                with mp.workdps(wp):
                    cur_y = [mpc(v) for v in cur_y]
                ys.append(list(cur_y))
                cur_x = p
        # stack Y(u_p) kappa = y_p
        rows = []
        rhs = []
        for pi_, p in enumerate(pts):
            Yp = basis["eval"](p - xs, wp)
            for i in range(n):
                rows.append([mpc(v) for v in Yp[i]])
                rhs.append(mpc(ys[pi_][i]) if not cond_only else mpc(0))
        Amat = mp.matrix(rows)
        bvec = mp.matrix(rhs)
        # column equilibration: exact reparametrization —
        # solve (A D^-1)(D kappa) = b, D = diag(max-abs of each column).
        col_scale = None
        if equilibrate_cols:
            col_scale = []
            for j in range(Amat.cols):
                m = max(abs(Amat[i, j]) for i in range(Amat.rows))
                col_scale.append(m if m > 0 else mpf(1))
            for i in range(Amat.rows):
                for j in range(Amat.cols):
                    Amat[i, j] = Amat[i, j] / col_scale[j]
        # condition-number estimate of the matching matrix (columns scale as
        # u^lambda_i — genuinely ill-conditioned for spread indicial
        # exponents); recorded so a clean-looking residual cannot hide a
        # near-degenerate direction.
        cond_est = None
        try:
            sv = mp.svd_c(Amat, compute_uv=False)
            smax = max(sv) if len(sv) else mpf(0)
            smin = min(sv) if len(sv) else mpf(0)
            cond_est = mp.nstr(smax / smin, 8) if smin > 0 else "inf"
        except (ValueError, ZeroDivisionError, ArithmeticError):
            cond_est = None
        if cond_only:
            u_tail_c = max(abs(p - xs) for p in pts)
            tail_c = mpf(basis["p_last_norm"]) * u_tail_c ** basis["kmax"]
            with mp.workdps(dps):
                return {
                    "cond_est": cond_est,
                    "match_design": match_design,
                    "equilibrated": bool(col_scale is not None),
                    "match_points": [+p for p in pts],
                    "series_tail": +tail_c,
                    "basis": basis,
                    "dps": dps, "kmax": kmax,
                    "cond_only": True,
                }
        solver = "qr"
        try:
            kappa_m, resnorm = mp.qr_solve(Amat, bvec)
        except (ValueError, ZeroDivisionError, ArithmeticError) as exc:
            # F4: narrow except — programming errors propagate;
            # only a genuine solver failure reaches here, and the
            # normal-equations fallback must be OPTED INTO (it squares the
            # condition number; wrong-in-last-digits kappa can emerge with a
            # clean-looking residual at these digit targets).
            if not allow_normal_eq:
                raise RuntimeError(
                    f"land: mp.qr_solve failed ({type(exc).__name__}: {exc}) "
                    f"on the {Amat.rows}x{Amat.cols} matching matrix "
                    f"(cond_est={cond_est}). Refusing the normal-equations "
                    f"fallback by default — A^H A squares the condition "
                    f"number and a small residual does not bound the kappa "
                    f"error in near-degenerate directions. Pass "
                    f"allow_normal_eq=True to accept it (recorded as "
                    f"solver='normal_eq'), or increase dps / move x_from.")
            solver = "normal_eq"
            AH = Amat.H
            kappa_m = mp.lu_solve(AH * Amat, AH * bvec)
            resnorm = mp.norm(Amat * kappa_m - bvec)
        if col_scale is not None:
            kappa = [mpc(kappa_m[i]) / col_scale[i] for i in range(n)]
        else:
            kappa = [mpc(kappa_m[i]) for i in range(n)]
        bnorm = mp.norm(bvec)
        res_rel = mpf(resnorm) / max(mpf(1), bnorm)
        u_tail = max(f * L for f in fracs)     # worst (largest) matching |u|
        tail = mpf(basis["p_last_norm"]) * u_tail ** basis["kmax"]
        pts_out = list(pts)
    with mp.workdps(dps):
        return {
            "kappa": [+v for v in kappa],
            "residual": +mpf(resnorm),
            "residual_rel": +res_rel,
            "solver": solver,
            "cond_est": cond_est,
            "match_design": match_design,
            "equilibrated": bool(col_scale is not None),
            "match_points": [+p for p in pts_out],
            "series_tail": +tail,
            "basis": basis,
            "dps": dps,
            "kmax": kmax,
        }
