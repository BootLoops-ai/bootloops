"""baller.block_arrow — generic block-arrow structured Krawczyk certification
+ block-arrow interval PD verification.

THE SHAPE. A hierarchical-MAP stationarity system has the arrow form: the
symmetric Jacobian (Hessian of the objective) is

    H = [ D_1              B_1 ]
        [      D_2         B_2 ]
        [           ...    ... ]
        [ B_1^T B_2^T ...  G   ]

with B independent diagonal blocks D_i (the latent blocks) and a dense
border G (the globals). Nothing here is model-specific: the ORACLE owns the
model; this module owns the structured certification.

ORACLE CONTRACT (fail-closed; violations raise BlockArrowContractError):
  oracle.dims  -> ([n_1, ..., n_B], n_g)
  oracle.F(x)  -> (Fz, Fg): Fz = list of B lists of flint.arb, Fg = list of
                  arb — the system value (gradient) enclosed over x.
  oracle.H(x)  -> (D, B, G): D = list of B arb_mat (n_i x n_i),
                  B = list of B arb_mat (n_i x n_g), G = arb_mat (n_g x n_g)
                  — the symmetric block-arrow Jacobian enclosed over x.
  x is (zs, g): zs = list of B lists of arb, g = list of arb. The oracle
  MUST return rigorous enclosures valid for every point of the input box
  (ball arithmetic end to end; python-flint arb only).

SEMANTICS (standard Krawczyk, blockwise assembled; Rump-2010 semantics):
with x~ the candidate center, X = x~ + [-r, r]^n, C the structured
approximate inverse of mid H (blockwise approximate inverse + interval
Schur complement on the border; the float preconditioner is heuristic, the
RIGOR never depends on it), the Krawczyk image is

    K(X) = x~ - C F(x~) + (I - C H(X)) (X - x~).

K(X) strictly interior to X proves a zero of F EXISTS in X and is UNIQUE
there. The per-coordinate test: mag(c_k) + (|E| r)_k < r with c = C F(x~),
E = I - C H(X); every magnitude is a certified upper bound and the final
comparison is arb-certified — an undecidable comparison REFUSES (a fat
ball is a refusal, never a shrug). The structured assembly never forms the
dense n x n E: with Y_i ~ inv(mid D_i), P_i = Y_i mid(B_i),
S^ = mid(G) - sum_i mid(B_i)^T Y_i mid(B_i), Ys ~ inv(S^),

    E_zz[i][j] = delta_ij (I - Y_i D_i) - P_i Q_j
    E_zg[i]    = P_i (Ys W) - Y_i B_i
    E_gz[j]    = Q_j,          Q_j = Ys (P_j^T D_j - B_j^T)
    E_gg       = I - Ys W,     W   = G - sum_j P_j^T B_j

and the cross-block row sums are bounded via
sum_j |P_i Q_j| <= |P_i| q,  q_m = sum_j sum_cols |Q_j|[m, :] —
cost LINEAR in the number of blocks, never quadratic.

REFUSALS name the failing subspace ('block_<i>' or 'border') and the
measured contraction defect (the margin that had to be provably < 1).
No knobs: no tolerance, no damping, no retry lives here — radius policy
belongs to the CALLER (CLINCH pins its own).

PD LEG. pd_block_arrow certifies that EVERY symmetric matrix in the
block-arrow enclosure is positive definite, by the same structure: each
D_i via certified congruence (float Cholesky R_i of the midpoint, X_i the
float triangular inverse — triangular with nonzero diagonal, hence EXACTLY
nonsingular, so the congruence preserves definiteness — then
M_i = X_i^T D_i X_i enclosed in arb; eps_i = max row sum of |M_i - I|
provably < 1 certifies M_i PD hence D_i PD, and gives the Neumann bound
|M_i^{-1} - I| <= eps_i/(1-eps_i) entrywise). The border Schur complement
S = G - sum_i B_i^T D_i^{-1} B_i is enclosed via
D_i^{-1} = X_i M_i^{-1} X_i^T (the Neumann term entering as a certified
entrywise ball), and S gets the same congruence certificate. All D_i PD
+ S PD <=> H PD (symmetric Schur). interval_cholesky is the scalar
reference implementation (small n), cross-checked in the battery.

Symmetry: both legs symmetrize defensively (A -> (A + A^T)/2 — an exact
no-op for a correctly assembled Hessian) so the PD statement is airtight
for the symmetric part the semantics require.
"""
import math

import numpy as np

__all__ = [
    "BlockArrowContractError", "block_krawczyk", "pd_block_arrow",
    "certify_local_min", "interval_cholesky", "box_around", "eval_oracle",
]


class BlockArrowContractError(TypeError):
    """Oracle/shape contract violation (typed; distinct from a REFUSAL)."""


def _flint():
    from flint import arb, arb_mat, ctx
    return arb, arb_mat, ctx


def _ctx_guard(prec):
    from baller.hygiene import ctx_guard
    return ctx_guard(prec=prec)


# ---------------------------------------------------------------------------
# input plumbing
# ---------------------------------------------------------------------------
def _check_dims(oracle):
    try:
        bl, ng = oracle.dims
        bl = [int(n) for n in bl]
        ng = int(ng)
    except Exception as e:
        raise BlockArrowContractError(f"oracle.dims unreadable: {e}")
    if not bl or any(n <= 0 for n in bl) or ng <= 0:
        raise BlockArrowContractError(f"bad dims: blocks={bl} border={ng}")
    return bl, ng


def _radius_struct(center, radius):
    """Normalize radius to per-coordinate structure ((rz, rg) floats).
    Accepts a scalar (uniform box) or an (rz, rg) structure matching
    center — the ANISOTROPIC box a stiff hierarchical Hessian needs
    (curvature-scaled radii; the caller owns the scaling policy)."""
    zs, g = center
    if isinstance(radius, (tuple, list)):
        rz, rg = radius
        rz = [[float(v) for v in blk] for blk in rz]
        rg = [float(v) for v in rg]
        if len(rz) != len(zs) or any(len(a) != len(b)
                                     for a, b in zip(rz, zs)) \
                or len(rg) != len(g):
            raise BlockArrowContractError("radius structure mismatch")
    else:
        r = float(radius)
        rz = [[r] * len(blk) for blk in zs]
        rg = [r] * len(g)
    for v in [x for blk in rz for x in blk] + list(rg):
        if not (v >= 0.0 and math.isfinite(v)):
            raise BlockArrowContractError(f"radius must be finite >= 0: {v}")
    return rz, rg


def box_around(center, radius):
    """(zs, g) box: coordinate k gets center_k + [-r_k, r_k] as an arb
    ball. center = (zs, g) with arb-convertible entries; radius = scalar
    (uniform) or (rz, rg) per-coordinate structure (0 -> tight)."""
    arb, _, _ = _flint()
    rz, rg = _radius_struct(center, radius)
    zs, g = center
    zbox = [[arb(v) + (arb(0, r) if r else arb(0))
             for v, r in zip(blk, rblk)] for blk, rblk in zip(zs, rz)]
    gbox = [arb(v) + (arb(0, r) if r else arb(0)) for v, r in zip(g, rg)]
    return zbox, gbox


def _mid_np(M):
    """Float midpoint matrix of an arb_mat (heuristic input only)."""
    n, m = M.nrows(), M.ncols()
    Mm = M.mid()
    return np.array([[float(Mm[i, j]) for j in range(m)] for i in range(n)])


def _np_to_mat(A):
    arb, arb_mat, _ = _flint()
    return arb_mat([[arb(float(v)) for v in row] for row in A])


def _colvec(v):
    arb, arb_mat, _ = _flint()
    return arb_mat([[x if isinstance(x, arb) else arb(x)] for x in v])


def _sym(M):
    arb, _, _ = _flint()
    return (M + M.transpose()) * arb("0.5")


def _ub_max(a, b):
    """arb upper bound for max(a, b): (a+b+|a-b|)/2 with |.| an upper
    bound — rigorous, never undecidable."""
    arb, _, _ = _flint()
    return (a + b + (a - b).abs_upper()) * arb("0.5")


def _row_abs_sums(M):
    """Certified row sums of entry magnitudes as arb balls (each entry's
    abs_upper is exact; the sum stays a rigorous upper enclosure)."""
    arb, _, _ = _flint()
    n, m = M.nrows(), M.ncols()
    out = []
    for i in range(n):
        s = arb(0)
        for j in range(m):
            s += M[i, j].abs_upper()
        out.append(s)
    return out


def _row_weighted(M, w):
    """Certified per-row sums sum_j |M[i,j]| * w_j (w: list of arb >= 0)."""
    arb, _, _ = _flint()
    n, m = M.nrows(), M.ncols()
    out = []
    for i in range(n):
        s = arb(0)
        for j in range(m):
            s += M[i, j].abs_upper() * w[j]
        out.append(s)
    return out


def _col_abs_sums(M):
    arb, _, _ = _flint()
    n, m = M.nrows(), M.ncols()
    out = [arb(0) for _ in range(m)]
    for i in range(n):
        for j in range(m):
            out[j] += M[i, j].abs_upper()
    return out


def _all_finite(M):
    n, m = M.nrows(), M.ncols()
    return all(M[i, j].is_finite() for i in range(n) for j in range(m))


def _fup(x):
    """Float upper bound of an arb (reporting only; may round)."""
    try:
        return float(x.abs_upper())
    except Exception:
        return float("inf")


def _refuse(failed, defect, why, extra=None):
    out = {"verdict": "REFUSED", "failed": failed,
           "defect": (None if defect is None else float(defect)),
           "reason": why}
    if extra:
        out.update(extra)
    return out


def eval_oracle(oracle, center, radius, prec):
    """One (F at the tight center, H over the box) evaluation pair, shape-
    checked. Both certificate legs consume this; call it once."""
    bl, ng = _check_dims(oracle)
    with _ctx_guard(prec):
        tight = box_around(center, 0.0)
        box = box_around(center, radius)
        Fz, Fg = oracle.F(tight)
        if len(Fz) != len(bl) or any(len(f) != n for f, n in zip(Fz, bl)) \
                or len(Fg) != ng:
            raise BlockArrowContractError("oracle.F shape mismatch")
        D, B, G = oracle.H(box)
        if len(D) != len(bl) or len(B) != len(bl):
            raise BlockArrowContractError("oracle.H block-count mismatch")
        for i, n in enumerate(bl):
            if D[i].nrows() != n or D[i].ncols() != n \
                    or B[i].nrows() != n or B[i].ncols() != ng:
                raise BlockArrowContractError(f"oracle.H block {i} shape")
        if G.nrows() != ng or G.ncols() != ng:
            raise BlockArrowContractError("oracle.H border shape")
    return (Fz, Fg), (D, B, G)


# ---------------------------------------------------------------------------
# the structured Krawczyk certificate
# ---------------------------------------------------------------------------
def block_krawczyk(oracle, center, radius, prec, fh=None):
    """Structured Krawczyk existence+uniqueness certificate on the box
    center +/- radius (uniform infinity-norm ball). Returns a verdict dict;
    REFUSED verdicts name the failing block/border and the contraction
    defect. Raises BlockArrowContractError only for contract violations.
    fh: optional precomputed eval_oracle(...) result (plumbing, not a knob)."""
    arb, arb_mat, _ = _flint()
    bl, ng = _check_dims(oracle)
    rz_f, rg_f = _radius_struct(center, radius)
    if not all(v > 0.0 for blk in rz_f for v in blk) \
            or not all(v > 0.0 for v in rg_f):
        raise BlockArrowContractError("radius must be > 0 per coordinate")
    if fh is None:
        fh = eval_oracle(oracle, center, radius, prec)
    (Fz, Fg), (D, B, G) = fh
    with _ctx_guard(prec):
        rzb = [[arb(v) for v in blk] for blk in rz_f]
        rgb = [arb(v) for v in rg_f]
        if not all(f.is_finite() for blk in Fz for f in blk) \
                or not all(f.is_finite() for f in Fg):
            return _refuse("nonfinite_oracle", None,
                           "F enclosure not finite",
                           {"tool": "baller.certify.block_krawczyk"})
        if not all(_all_finite(M) for M in D) or \
                not all(_all_finite(M) for M in B) or not _all_finite(G):
            return _refuse("nonfinite_oracle", None,
                           "H enclosure not finite",
                           {"tool": "baller.certify.block_krawczyk"})

        # --- structured approximate inverse of mid H (floats; heuristic) --
        Y_np, P_np = [], []
        S_np = _mid_np(G)
        for i in range(len(bl)):
            Dm = _mid_np(D[i])
            Bm = _mid_np(B[i])
            try:
                Yi = np.linalg.inv(Dm)
            except np.linalg.LinAlgError:
                return _refuse(f"block_{i:02d}", None,
                               "midpoint block not invertible "
                               "(approximate-inverse stage)",
                               {"tool": "baller.certify.block_krawczyk"})
            if not np.all(np.isfinite(Yi)):
                return _refuse(f"block_{i:02d}", None,
                               "midpoint inverse not finite",
                               {"tool": "baller.certify.block_krawczyk"})
            Pi = Yi @ Bm
            S_np -= Bm.T @ Pi
            Y_np.append(Yi)
            P_np.append(Pi)
        try:
            Ys_np = np.linalg.inv(S_np)
        except np.linalg.LinAlgError:
            return _refuse("border", None,
                           "midpoint border Schur not invertible",
                           {"tool": "baller.certify.block_krawczyk"})
        if not np.all(np.isfinite(Ys_np)):
            return _refuse("border", None,
                           "midpoint Schur inverse not finite",
                           {"tool": "baller.certify.block_krawczyk"})

        Y = [_np_to_mat(A) for A in Y_np]
        P = [_np_to_mat(A) for A in P_np]
        Ys = _np_to_mat(Ys_np)

        # --- c = C F(x~), blockwise (the SAME C the E-assembly uses) ------
        u = [Y[i] * _colvec(Fz[i]) for i in range(len(bl))]
        acc = _colvec(Fg)
        for j in range(len(bl)):
            acc = acc - P[j].transpose() * _colvec(Fz[j])
        c_g = Ys * acc
        c_z = [u[i] - P[i] * c_g for i in range(len(bl))]

        # --- border helpers ----------------------------------------------
        W = G
        for j in range(len(bl)):
            W = W - P[j].transpose() * B[j]
        YsW = Ys * W
        # q_m = sum_j sum_l |Q_j[m,l]| * rz_j[l]  (radius-weighted)
        q = [arb(0) for _ in range(ng)]
        for j in range(len(bl)):
            Qj = Ys * (P[j].transpose() * D[j] - B[j].transpose())
            rs = _row_weighted(Qj, rzb[j])
            for m in range(ng):
                q[m] += rs[m]

        margins = {}
        newton_inf = arb(0)

        # border rows: mag(c_g_k) + q_k + sum_j |I - YsW|[k,j] rg_j < rg_k
        m_border = arb(0)
        border_ok = True
        for k in range(ng):
            s = arb(0)
            for j in range(ng):
                e = (arb(1) if k == j else arb(0)) - YsW[k, j]
                s += e.abs_upper() * rgb[j]
            tot = c_g[k, 0].abs_upper() + q[k] + s
            if not (tot < rgb[k]):
                border_ok = False
            m_border = _ub_max(m_border, tot / rgb[k])
            newton_inf = _ub_max(newton_inf, c_g[k, 0].abs_upper())
        margins["border"] = (border_ok, m_border)

        # latent rows, block i
        for i in range(len(bl)):
            n = bl[i]
            A = Y[i] * D[i]                          # interval, n x n
            Ezg = P[i] * YsW - Y[i] * B[i]           # n x ng interval
            ezg_rs = _row_weighted(Ezg, rgb)
            ok_i = True
            m_i = arb(0)
            for k in range(n):
                s = arb(0)
                for j in range(n):
                    e = (arb(1) if k == j else arb(0)) - A[k, j]
                    s += e.abs_upper() * rzb[i][j]
                cross = arb(0)
                for m in range(ng):
                    cross += P[i][k, m].abs_upper() * q[m]
                tot = (c_z[i][k, 0].abs_upper() + s + cross + ezg_rs[k])
                if not (tot < rzb[i][k]):
                    ok_i = False
                m_i = _ub_max(m_i, tot / rzb[i][k])
                newton_inf = _ub_max(newton_inf, c_z[i][k, 0].abs_upper())
            margins[f"block_{i:02d}"] = (ok_i, m_i)

        certified = all(ok for ok, _ in margins.values())
        pool = ([(k, m) for k, (ok, m) in margins.items() if not ok]
                if not certified else list(
                    (k, m) for k, (ok, m) in margins.items()))
        worst_name, worst_val = max(
            ((k, _fup(m)) for k, m in pool), key=lambda t: t[1])
        r_all = [v for blk in rz_f for v in blk] + list(rg_f)
        out = {
            "tool": "baller.certify.block_krawczyk",
            "prec_bits": int(prec),
            "radius": max(r_all),
            "radius_min": min(r_all),
            "radius_anisotropic": bool(min(r_all) != max(r_all)),
            "dims": {"blocks": bl, "border": ng},
            "margins": {k: _fup(m) for k, (ok, m) in margins.items()},
            "margin_ok": {k: bool(ok) for k, (ok, _) in margins.items()},
            "max_margin": max(_fup(m) for _, m in margins.values()),
            "newton_step_inf": _fup(newton_inf),
        }
        if certified:
            out.update({"verdict": "CERTIFIED", "failed": None,
                        "defect": None,
                        "reason": "K(X) strictly interior to X: a zero of F "
                                  "EXISTS in the ball and is UNIQUE there"})
        else:
            out.update({"verdict": "REFUSED", "failed": worst_name,
                        "defect": worst_val,
                        "reason": (f"contraction not proven on {worst_name}: "
                                   f"margin {worst_val:.6g} not provably < 1 "
                                   f"(fat ball or true defect — both "
                                   f"refuse)")})
        return out


# ---------------------------------------------------------------------------
# PD verification
# ---------------------------------------------------------------------------
def _congruence_eps(M_int, tag):
    """Certified congruence PD step with EXACT power-of-two equilibration
    (ill-conditioned blocks: the float Cholesky residual scales with the
    condition number; the stiffness is mostly diagonal scale, and dividing
    row/column i by 2^round(log2 sqrt(A_ii)) is EXACT in both float and
    arb, an exact congruence — PD is preserved both ways). Then:
    R = mid Cholesky of the equilibrated A' -> X' its float triangular
    inverse (exactly nonsingular) -> M = X'^T A' X' in arb -> eps =
    certified upper bound of max row sum |M - I|. Returns
    (eps: arb, X_eff: np array with A^{-1} = X_eff M^{-1} X_eff^T) or a
    refusal dict."""
    arb, arb_mat, _ = _flint()
    n = M_int.nrows()
    A = _sym(M_int)
    Am0 = _mid_np(A)
    Am0 = 0.5 * (Am0 + Am0.T)
    d = np.abs(np.diag(Am0))
    d = np.where(d > 0, d, 1.0)
    inv_s = 2.0 ** (-np.round(0.5 * np.log2(d)))     # exact powers of two
    A2 = arb_mat(n, n)
    for i in range(n):
        for j in range(n):
            A2[i, j] = A[i, j] * (inv_s[i] * inv_s[j])
    Am = _mid_np(A2)
    Am = 0.5 * (Am + Am.T)
    try:
        R = np.linalg.cholesky(Am).T          # upper triangular
    except np.linalg.LinAlgError:
        w, v = np.linalg.eigh(Am)
        return _refuse(f"midpoint_cholesky:{tag}", float(w[0]),
                       f"{tag}: midpoint not numerically PD "
                       f"(min eig {w[0]:.6g} after equilibration; "
                       f"direction is a float diagnostic, not a "
                       f"certificate)",
                       {"diagnostic_negative_direction":
                        [float(x) for x in (inv_s * v[:, 0])]})
    Rinv = np.linalg.solve(R, np.eye(n))      # triangular, diag != 0
    if not np.all(np.isfinite(Rinv)):
        return _refuse(f"midpoint_cholesky:{tag}", None,
                       f"{tag}: triangular inverse not finite")
    X = _np_to_mat(Rinv)
    M = X.transpose() * A2 * X
    eps = arb(0)
    for i in range(n):
        s = arb(0)
        for j in range(n):
            e = M[i, j] - (arb(1) if i == j else arb(0))
            s += e.abs_upper()
        eps = _ub_max(eps, s)
    # A^{-1} = (S^{-1} X') M^{-1} (S^{-1} X')^T; S^{-1} exact -> exact rows
    X_eff = inv_s[:, None] * Rinv
    return (eps, X_eff)


def _equilibrated_chol(Am, tag):
    """Float Cholesky of the exact-power-of-two equilibrated midpoint.
    Returns (X = S^{-1} R^{-1}: np upper-triangular with nonzero diagonal
    — EXACTLY nonsingular) or a refusal dict with the float diagnostic
    direction. A' = S^{-1} A S^{-1} ~ R^T R  =>  A ~ (RS)^T (RS)."""
    d = np.abs(np.diag(Am))
    d = np.where(d > 0, d, 1.0)
    inv_s = 2.0 ** (-np.round(0.5 * np.log2(d)))
    A2 = (inv_s[:, None] * Am) * inv_s[None, :]
    A2 = 0.5 * (A2 + A2.T)
    try:
        R = np.linalg.cholesky(A2).T
    except np.linalg.LinAlgError:
        w, v = np.linalg.eigh(A2)
        return _refuse(f"midpoint_cholesky:{tag}", float(w[0]),
                       f"{tag}: midpoint not numerically PD "
                       f"(min eig {w[0]:.6g} after equilibration; "
                       f"direction is a float diagnostic, not a "
                       f"certificate)",
                       {"diagnostic_negative_direction":
                        [float(x) for x in (inv_s * v[:, 0])]})
    Rinv = np.linalg.solve(R, np.eye(len(d)))
    if not np.all(np.isfinite(Rinv)):
        return _refuse(f"midpoint_cholesky:{tag}", None,
                       f"{tag}: triangular inverse not finite")
    return inv_s[:, None] * Rinv


def pd_block_arrow(D, B, G, prec):
    """Certify that EVERY symmetric matrix in the block-arrow interval
    enclosure (D_i, B_i, G) is positive definite, by ONE global certified
    congruence with the block-arrow approximate Cholesky factor:
    X = [[Xz (block-diag), C], [0, Xg]] (float, block-triangular with
    triangular nonsingular diagonal blocks — EXACTLY nonsingular however
    the floats round), M = X^T H X enclosed blockwise in arb; every row
    sum of |M - I| provably < 1 certifies M PD hence H PD. No rigorous
    inverse anywhere: ill-conditioning only degrades how close M is to I
    (measured on a production-scale hierarchical system: the
    per-block-Schur routes fail — the Neumann
    bound rank-collapses and the certified interval solve fattens by
    ~1e2 on cond ~ 1e13 collinear-latent blocks; the global congruence
    certifies the same enclosure with eps ~ 1e-3). Refusals name the
    failing block/border row group."""
    arb, arb_mat, _ = _flint()
    nb = len(D)
    with _ctx_guard(prec):
        ng = G.nrows()
        one = arb(1)
        margins = {}
        # ---- float factors (heuristic; rigor lives in the congruence) --
        Dm = [0.5 * (_mid_np(M) + _mid_np(M).T) for M in D]
        Bm = [_mid_np(M) for M in B]
        Gm = 0.5 * (_mid_np(G) + _mid_np(G).T)
        Xz_np, U_np = [], []
        Sm = Gm.copy()
        for i in range(nb):
            got = _equilibrated_chol(Dm[i], f"block_{i:02d}")
            if isinstance(got, dict):
                got["tool"] = "baller.certify.pd_block_arrow"
                return got
            Xz_np.append(got)
            Ui = Bm[i].T @ got                       # ng x n_i
            U_np.append(Ui)
            Sm -= Ui @ Ui.T
        got = _equilibrated_chol(Sm, "border_schur")
        if isinstance(got, dict):
            got["tool"] = "baller.certify.pd_block_arrow"
            return got
        Xg_np = got
        Xg = _np_to_mat(Xg_np)
        # ---- M = X^T H X blockwise in arb ------------------------------
        Gsym = _sym(G)
        Mgg = Xg.transpose() * Gsym * Xg
        border_extra = [arb(0) for _ in range(ng)]   # sum_i |M_zg,i| cols
        worst_block = None
        for i in range(nb):
            n = D[i].nrows()
            Xz = _np_to_mat(Xz_np[i])
            Ci_np = -(Xz_np[i] @ (U_np[i].T @ Xg_np))
            Ci = _np_to_mat(Ci_np)
            Dsym = _sym(D[i])
            Mzz = Xz.transpose() * Dsym * Xz
            Ti = Dsym * Ci + B[i] * Xg
            Mzg = Xz.transpose() * Ti
            # M_gg accumulation: C^T D C + C^T B Xg + (C^T B Xg)^T
            CB = Ci.transpose() * (B[i] * Xg)
            Mgg = Mgg + Ci.transpose() * (Dsym * Ci) + CB + CB.transpose()
            ok_i = True
            eps_i = arb(0)
            for k in range(n):
                s = arb(0)
                for j in range(n):
                    e = Mzz[k, j] - (arb(1) if k == j else arb(0))
                    s += e.abs_upper()
                for l in range(ng):
                    a = Mzg[k, l].abs_upper()
                    s += a
                if not (s < one):
                    ok_i = False
                eps_i = _ub_max(eps_i, s)
            for l in range(ng):
                col = arb(0)
                for k in range(n):
                    col += Mzg[k, l].abs_upper()
                border_extra[l] += col
            m_i = float((one - eps_i).lower())
            margins[f"block_{i:02d}"] = m_i
            if not ok_i:
                return _refuse(
                    f"block_{i:02d}", _fup(eps_i),
                    f"block_{i:02d}: congruence row sum {_fup(eps_i):.6g} "
                    f"not provably < 1 (fat ball or true indefiniteness)",
                    {"tool": "baller.certify.pd_block_arrow",
                     "margins": margins})
        eps_g = arb(0)
        ok_g = True
        for k in range(ng):
            s = border_extra[k]
            for j in range(ng):
                e = Mgg[k, j] - (arb(1) if k == j else arb(0))
                s += e.abs_upper()
            if not (s < one):
                ok_g = False
            eps_g = _ub_max(eps_g, s)
        margins["border_schur"] = float((one - eps_g).lower())
        if not ok_g:
            return _refuse(
                "border_schur", _fup(eps_g),
                f"border congruence row sum {_fup(eps_g):.6g} not "
                f"provably < 1 (fat ball or true indefiniteness)",
                {"tool": "baller.certify.pd_block_arrow",
                 "margins": margins})
        return {"tool": "baller.certify.pd_block_arrow",
                "verdict": "PD_CERTIFIED", "failed": None, "defect": None,
                "prec_bits": int(prec), "margins": margins,
                "reason": "every symmetric matrix in the enclosure is "
                          "positive definite (global block-arrow "
                          "congruence, all rows contract)"}


def certify_local_min(oracle, center, radius, prec):
    """The combined certificate CLINCH consumes: ONE oracle evaluation pair
    feeds both legs. CERTIFIED means: inside the ball a stationary point
    EXISTS, is UNIQUE, and every Hessian in the enclosure is PD (a
    certified strict local minimum). Any leg's refusal is the verdict."""
    fh = eval_oracle(oracle, center, radius, prec)
    kw = block_krawczyk(oracle, center, radius, prec, fh=fh)
    if kw["verdict"] != "CERTIFIED":
        kw["leg"] = "krawczyk"
        return {"verdict": "REFUSED", "krawczyk": kw, "pd": None,
                "failed": kw["failed"], "reason": kw["reason"]}
    D, B, G = fh[1]
    pd = pd_block_arrow(D, B, G, prec)
    if pd["verdict"] != "PD_CERTIFIED":
        pd["leg"] = "pd"
        return {"verdict": "REFUSED", "krawczyk": kw, "pd": pd,
                "failed": pd["failed"], "reason": pd["reason"]}
    return {"verdict": "CERTIFIED", "krawczyk": kw, "pd": pd,
            "failed": None,
            "reason": "unique stationary point in the ball + PD Hessian "
                      "throughout the ball: certified strict local minimum"}


# ---------------------------------------------------------------------------
# scalar reference: interval Cholesky (small n; battery cross-check)
# ---------------------------------------------------------------------------
def interval_cholesky(A):
    """Scalar interval Cholesky on an arb_mat enclosure. Succeeds iff every
    pivot is provably > 0, which certifies PD for every symmetric matrix in
    the enclosure. Returns (ok, min_pivot_lower: float|None, fail_index)."""
    arb, arb_mat, _ = _flint()
    A = _sym(A)
    n = A.nrows()
    L = [[arb(0)] * n for _ in range(n)]
    minpiv = None
    for j in range(n):
        s = A[j, j]
        for k in range(j):
            s -= L[j][k] * L[j][k]
        if not (s > 0):
            return (False, minpiv, j)
        piv = s.sqrt()
        lo = float(s.lower())
        minpiv = lo if minpiv is None else min(minpiv, lo)
        L[j][j] = piv
        for i in range(j + 1, n):
            t = A[i, j]
            for k in range(j):
                t -= L[i][k] * L[j][k]
            L[i][j] = t / piv
    return (True, minpiv, None)
