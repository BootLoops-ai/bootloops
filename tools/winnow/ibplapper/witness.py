"""witness.py — retrofit lambda-witness emission (v1 witnesses).

v1 witnesses = RETROFIT mode — per requested pivot row, a Wiedemann
lambda-solve against the caller's rows (measured economics: ~2.2-2.3
s/row/prime at the 7,135-eq reference family, min-poly reuse 5.5x,
KB-class certificates, density mean 3.3%), then an INDEPENDENT verify
through the vendored RECEIPT core (~1 ms/row). This module is glue around
the vendored emitter (ibplapper.receipt.emitter — imported, never forked).

Square systems (n_eqs == n_pivots) take the direct Abar^T Wiedemann path
(emitter.solve_target). Rectangular systems (redundant rows present) take
the Kaltofen-Saunders normal-equation path only — same emitter primitives,
every candidate exactly verified before acceptance.

Convention: a closed row r means  t + sum_m r[m]*m = 0, so the certified
coefficients are c_m = (-r[m]) mod p and the witness identity is
sum_i lam_i R_i == e_t - sum_m c_m e_m (WITNESS_FORMAT v1.0).
"""
import os
import time

from .receipt import core
from .receipt.emitter import FpSystem, emit_for_target, wiedemann_solve


def _solve_rectangular(S, t_pidx, rng, master_cols, claimed, p):
    """KS normal-equation path for n_eqs != n_pivots (emitter primitives;
    the direct path assumes a square partition)."""
    import numpy as np
    b_vec = np.zeros(S.n, np.int64)
    b_vec[t_pidx] = 1
    D = rng.integers(1, p, S.n_eqs, dtype=np.int64)
    mvG = lambda x: S.mv_AT((D * S.mv_A(x)) % p)            # noqa: E731
    t0 = time.time()
    z, st = wiedemann_solve(mvG, S.n, b_vec, p, rng)
    rec = {"p": p, "solve_wall_s": round(time.time() - t0, 3),
           "matvecs": st.get("matvecs"), "mode": st.get("mode"),
           "path": "KS_rectangular"}
    if z is None:
        rec["status"] = "SOLVE_FAILED"
        return rec, None
    lam = (D * S.mv_A(z)) % p
    # exact residual check on the pivot block
    if ((S.mv_AT(lam) - b_vec) % p).any():
        rec["status"] = "SOLVE_FAILED"
        return rec, None
    r_out = S.mv_out_T(lam)
    c, shell_hit = {}, 0
    for j, w in enumerate(S.others):
        if r_out[j]:
            if w in master_cols:
                c[w] = int((-int(r_out[j])) % p)
            else:
                shell_hit += 1
    if shell_hit:
        rec["status"] = "SHELL_RESIDUAL"
        rec["shell_hit_cols"] = shell_hit
        return rec, None
    rec["c"] = c
    rec["lambda_nnz"] = int((lam != 0).sum())
    if claimed is not None and claimed != c:
        rec["status"] = "CLAIM_MISMATCH"
    else:
        rec["status"] = "CERTIFIED"
    return rec, lam


def retrofit(system, closed_subs, targets, masters=None, write_dir=None,
             seed=20260706, family=None, point=None):
    """Emit + independently verify witnesses for `targets` (pivot columns of
    `closed_subs`). Returns dict target -> record.

    masters: allowed coefficient columns; DEFAULT = every non-pivot column.
    If the system carries shell columns (dump-closure gap), pass an explicit
    master set so shell-touching rows are refused certificates
    (SHELL_RESIDUAL) instead of certified against shell columns.
    """
    import numpy as np
    p = system.p
    rows = system.rows
    pivots = sorted(closed_subs)
    pivot_set = set(pivots)
    for t in targets:
        if t not in pivot_set:
            raise ValueError(f"witness target {t} is not a solved pivot")
    if masters is None:
        masters = set(system.all_cols) - pivot_set
    else:
        masters = set(masters)
    rows_terms = [list(r.items()) for r in rows]
    S = FpSystem(rows_terms, pivots, system.all_cols, p)
    square = (S.n_eqs == S.n)
    rng = np.random.default_rng(seed)
    fpr = core.system_fingerprint(rows, p)
    if write_dir:
        os.makedirs(write_dir, exist_ok=True)

    out, minpoly = {}, None
    for t in targets:
        claimed = {m: (-v) % p for m, v in closed_subs[t].items()}
        if square:
            rec, lam, mp_new = emit_for_target(S, t, masters, rng,
                                               minpoly=minpoly,
                                               claimed=claimed)
            if mp_new is not None:
                minpoly = mp_new
        else:
            rec, lam = _solve_rectangular(S, S.pidx[t], rng, masters,
                                          claimed, p)
        rec["target_col"] = int(t)
        if lam is not None and rec["status"] == "CERTIFIED":
            nz = np.flatnonzero(lam)
            lam_idx = [int(i) for i in nz]
            lam_val = [int(lam[i]) for i in nz]
            wit = {"p": p, "target_col": int(t), "c": rec["c"],
                   "lam_idx": lam_idx, "lam_val": lam_val,
                   "n_rows": S.n_eqs, "meta": {}, "path": "<in-memory>"}
            ok, det = core.verify_row(rows, wit)
            rec["verify_pass"] = bool(ok)
            if not ok:                          # loud, never silent
                rec["status"] = "VERIFY_FAILED"
                rec["verify_detail"] = det
            rec["lam"] = dict(zip(lam_idx, lam_val))
            if write_dir and ok:
                wp = os.path.join(write_dir, f"w_{p}_{t}.json")
                core.write_witness(
                    wp, p=p, target_col=t, c=rec["c"], lam_idx=lam_idx,
                    lam_val=lam_val, n_rows=S.n_eqs, family=family,
                    point=point,
                    system={"n_rows": S.n_eqs, "fingerprint": fpr})
                rec["witness_path"] = wp
        out[t] = rec
    return out
