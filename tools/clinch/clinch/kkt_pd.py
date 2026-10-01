"""clinch.kkt_pd — the reduced-Hessian PD certificate for CENTERED
(extended-KKT) spaces: the rank-4 congruence extension named in
oracle_v36's SEMANTICS NOTE, completing the Krawczyk leg's certificate.

THE PROBLEM. The extended system of a centered space is a KKT system:
variables x = (theta, m, mu) with primal theta+m and one multiplier row
per constraint (mu_s, mu_a, and lam under the pin). Its Jacobian is the
symmetric bordered matrix

    K = [ H_L  A^T ]      H_L = Lagrangian Hessian over the PRIMAL
        [ A    0   ]            variables, A = constraint Jacobian,

indefinite BY CONSTRUCTION — with q constraints K carries at least q
negative eigenvalues, so the plain PD leg (pd_block_arrow, target I)
does not apply. The statement that certifies a strict local minimum of
the CONSTRAINED problem is positive definiteness of the REDUCED Hessian
Z^T H_L Z, Z a basis of the constraint tangent space null(A).

THE CERTIFICATE (signed congruence -> inertia -> reduced PD):

1. SIGNED CONGRUENCE. Build ONE float block-triangular factor
   X = [[Xz (block-diag latent Cholesky factors), C], [0, Xg]] exactly as
   the PD leg does — except the border factor Xg comes from an
   equilibrated float EIGENDECOMPOSITION of the midpoint border Schur
   complement (which is indefinite), scaled to |eigenvalue|^{-1/2}, so
   X^T mid(K) X ~ T with the SIGNED target T = diag(+1 on every latent
   row, sign pattern of the border eigenvalues on the border rows).
   The construction refuses unless the float sign pattern shows EXACTLY
   n_mult minus signs (n_mult = the number of multiplier rows the oracle
   declares) — the 'kkt-signature' refusal, the honest verdict when the
   reduced Hessian is indefinite.
2. RIGOR. M = X^T K X is enclosed blockwise in arb over the ENTIRE
   H-enclosure; eps = certified upper bound of the max row sum of
   |M - T| must be provably < 1. Then for every real symmetric K0 in the
   enclosure, ||X^T K0 X - T||_2 <= ||.||_inf < 1 (symmetric), so by
   Weyl every eigenvalue of X^T K0 X lies within distance < 1 of the
   matching +/-1 eigenvalue of T: X^T K0 X is NONSINGULAR with inertia
   In(T) = (p, q, 0). M nonsingular forces X nonsingular (X v = 0 would
   give M v = 0), so Sylvester's law of inertia gives In(K0) = (p, q, 0)
   — the float factor is heuristic, the inertia claim is certified.
3. SEMANTICS. K0 nonsingular forces A to have full row rank (a
   dependence A^T y = 0, y != 0 would make K0 [0; y] = 0). The classical
   KKT inertia identity (Gould 1985; Haynsworth inertia additivity) then
   gives  In(K0) = In(Z^T H_L Z) + (q, q, 0)  with q = n_mult, so
   q certified negative eigenvalues and NONE anywhere else means
   In(Z^T H_L Z) = (p - q, 0, 0): the reduced Hessian is positive
   definite for EVERY matrix in the enclosure. Combined with the
   Krawczyk leg (a KKT point exists in the box and is unique), that is
   second-order sufficiency: a certified strict local minimum of the
   equality-constrained problem.

CONTRACT (fail-closed): the oracle declares its multiplier rows as
border offsets (kkt_multiplier_idx); the corresponding diagonal block of
G must be EXACTLY zero balls (the defining KKT structure — anything else
raises the baller contract error, typed, distinct from a refusal).

'RANK-4': the unpinned centered space borders the model Hessian by the 4
auxiliary rows (m_s, m_a, mu_s, mu_a) — the reduced-Hessian statement
rides a rank-4 bordering (rank-5 under the pin). Nothing here is
model-specific: n_mult is read from the declaration, the toy battery
fixture runs the same code at rank-2.

REFUSAL CLASSES: midpoint_cholesky:block_<i> (a latent block's midpoint
is not numerically PD — float diagnostic direction attached);
midpoint_singular:border_schur; kkt-signature (midpoint negative count
!= n_mult — the planted-indefinite control lands here); block_<i> /
border_schur (congruence row sum not provably < 1: fat ball or true
defect — both refuse).

Kernel policy: the generic congruence/row-sum machinery is consumed from
baller.certify.block_krawczyk BY IDENTITY (helpers, contract error,
refusal shape); this module owns only the KKT semantics — the same
split the corner KKT leg uses.
"""
import numpy as np

__all__ = ["pd_reduced_congruence"]


def _bk():
    from .engine import _baller
    return _baller()


def _signed_eig_factor(bk, Sm, n_mult, tag):
    """Equilibrated float eigenfactor of the (indefinite) midpoint border
    Schur complement: X with X^T Sm X ~ diag(signs), signs ascending
    (negatives first). Returns (X_eff, signs) or a refusal dict."""
    d = np.abs(np.diag(Sm))
    d = np.where(d > 0, d, 1.0)
    inv_s = 2.0 ** (-np.round(0.5 * np.log2(d)))     # exact powers of two
    S2 = (inv_s[:, None] * Sm) * inv_s[None, :]
    S2 = 0.5 * (S2 + S2.T)
    w, V = np.linalg.eigh(S2)
    if np.any(w == 0.0) or not np.all(np.isfinite(w)):
        return bk._refuse(
            f"midpoint_singular:{tag}", None,
            f"{tag}: midpoint eigenvalue exactly zero/nonfinite — no "
            f"signed factor (float diagnostic, not a certificate)")
    neg = int(np.sum(w < 0.0))
    if neg != n_mult:
        return bk._refuse(
            "kkt-signature", float(w[min(neg, len(w) - 1)]),
            f"{tag}: midpoint inertia has {neg} negative "
            f"direction(s) where the KKT structure requires exactly "
            f"{n_mult} (one per multiplier row) — the reduced Hessian "
            f"is not numerically PD (float diagnostic; the certificate "
            f"is withheld, never forced)",
            {"negative_count": neg, "expected_negative_count": n_mult,
             "diagnostic_eigenvalues": [float(x) for x in w]})
    X = V * (np.abs(w) ** -0.5)[None, :]
    return (inv_s[:, None] * X, np.sign(w))


def pd_reduced_congruence(D, B, G, mult_idx, prec):
    """Certify that the REDUCED Hessian Z^T H_L Z is positive definite
    for EVERY symmetric matrix in the extended-KKT block-arrow enclosure
    (D_i, B_i, G) — the signed-congruence inertia certificate documented
    in the module header. mult_idx: border offsets of the multiplier
    rows (the oracle's kkt_multiplier_idx). Returns a verdict dict:
    PD_REDUCED_CERTIFIED with per-block margins and the certified
    inertia, or a named refusal."""
    bk = _bk()
    from flint import arb
    nb = len(D)
    ng = G.nrows()
    mult = sorted(set(int(a) for a in mult_idx))
    if not mult or len(mult) >= ng or any(a < 0 or a >= ng for a in mult):
        raise bk.BlockArrowContractError(
            f"kkt multiplier rows {mult} not a proper subset of the "
            f"{ng}-row border")
    n_mult = len(mult)
    with bk._ctx_guard(prec):
        one = arb(1)
        # ---- KKT structure contract: multiplier block EXACTLY zero ----
        for a in mult:
            for b in mult:
                if not G[a, b].is_zero():
                    raise bk.BlockArrowContractError(
                        f"multiplier block G[{a},{b}] is not the exact "
                        f"zero ball — not a KKT system (the multiplier "
                        f"rows must carry a zero diagonal block)")
        # ---- float factors (heuristic; rigor lives in the congruence) --
        Dm = [0.5 * (bk._mid_np(M) + bk._mid_np(M).T) for M in D]
        Bm = [bk._mid_np(M) for M in B]
        Gm = 0.5 * (bk._mid_np(G) + bk._mid_np(G).T)
        Xz_np, U_np = [], []
        Sm = Gm.copy()
        for i in range(nb):
            got = bk._equilibrated_chol(Dm[i], f"block_{i:02d}")
            if isinstance(got, dict):
                got["tool"] = "clinch.kkt_pd.pd_reduced_congruence"
                return got
            Xz_np.append(got)
            Ui = Bm[i].T @ got
            U_np.append(Ui)
            Sm -= Ui @ Ui.T
        got = _signed_eig_factor(bk, Sm, n_mult, "border_schur")
        if isinstance(got, dict):
            got["tool"] = "clinch.kkt_pd.pd_reduced_congruence"
            return got
        Xg_np, signs = got
        Xg = bk._np_to_mat(Xg_np)
        tgt = [arb(-1) if s < 0 else arb(1) for s in signs]
        # ---- M = X^T K X blockwise in arb vs the SIGNED target --------
        Gsym = bk._sym(G)
        Mgg = Xg.transpose() * Gsym * Xg
        border_extra = [arb(0) for _ in range(ng)]
        margins = {}
        for i in range(nb):
            n = D[i].nrows()
            Xz = bk._np_to_mat(Xz_np[i])
            Ci = bk._np_to_mat(-(Xz_np[i] @ (U_np[i].T @ Xg_np)))
            Dsym = bk._sym(D[i])
            Mzz = Xz.transpose() * Dsym * Xz
            Ti = Dsym * Ci + B[i] * Xg
            Mzg = Xz.transpose() * Ti
            CB = Ci.transpose() * (B[i] * Xg)
            Mgg = Mgg + Ci.transpose() * (Dsym * Ci) + CB + CB.transpose()
            ok_i = True
            eps_i = arb(0)
            for k in range(n):
                s = arb(0)
                for j in range(n):
                    e = Mzz[k, j] - (one if k == j else arb(0))
                    s += e.abs_upper()
                for l in range(ng):
                    s += Mzg[k, l].abs_upper()
                if not (s < one):
                    ok_i = False
                eps_i = bk._ub_max(eps_i, s)
            for l in range(ng):
                col = arb(0)
                for k in range(n):
                    col += Mzg[k, l].abs_upper()
                border_extra[l] += col
            margins[f"block_{i:02d}"] = float((one - eps_i).lower())
            if not ok_i:
                return bk._refuse(
                    f"block_{i:02d}", bk._fup(eps_i),
                    f"block_{i:02d}: signed-congruence row sum "
                    f"{bk._fup(eps_i):.6g} not provably < 1 (fat ball or "
                    f"true defect)",
                    {"tool": "clinch.kkt_pd.pd_reduced_congruence",
                     "margins": margins})
        eps_g = arb(0)
        ok_g = True
        for k in range(ng):
            s = border_extra[k]
            for j in range(ng):
                e = Mgg[k, j] - (tgt[k] if k == j else arb(0))
                s += e.abs_upper()
            if not (s < one):
                ok_g = False
            eps_g = bk._ub_max(eps_g, s)
        margins["border_schur"] = float((one - eps_g).lower())
        if not ok_g:
            return bk._refuse(
                "border_schur", bk._fup(eps_g),
                f"border signed-congruence row sum {bk._fup(eps_g):.6g} "
                f"not provably < 1 (fat ball or true defect)",
                {"tool": "clinch.kkt_pd.pd_reduced_congruence",
                 "margins": margins})
        n_lat = sum(M.nrows() for M in D)
        p = n_lat + int(np.sum(signs > 0))
        return {"tool": "clinch.kkt_pd.pd_reduced_congruence",
                "verdict": "PD_REDUCED_CERTIFIED", "failed": None,
                "defect": None, "prec_bits": int(prec),
                "margins": margins,
                "inertia": [int(p), int(n_mult), 0],
                "n_mult": int(n_mult),
                "reduced_dim": int(p - n_mult),
                "reason": (f"certified inertia ({p}, {n_mult}, 0) of the "
                           f"extended KKT Jacobian over the whole "
                           f"enclosure (signed congruence, all rows "
                           f"contract): exactly one negative direction "
                           f"per multiplier row, so by the KKT inertia "
                           f"identity the reduced Hessian "
                           f"({p - n_mult} coordinates on the constraint "
                           f"tangent space) is positive definite "
                           f"throughout")}
