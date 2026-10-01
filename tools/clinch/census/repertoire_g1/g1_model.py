"""Gate-1 — model family, anchored cell grid, float pre-pass.

Family of record (pinned parameter box): De Boer-Perelson clone-niche
competition, dx_i/dt = x_i (sum_j b_ij phi_j(y_j) - d_i) + s_i,
y_j = sum_k c_kj x_k, phi(y) = 1/(1+(y/K)^h) (saturating consumption class).

Census object = the s=0 (peripheral-maintenance) competition skeleton, the
support-pattern object of the De Boer-Perelson
competitive-exclusion structure (den Braber B5 lower edge; thymic share is a
box coordinate whose s>0 interior is OUT OF SCOPE for Gate-1 and fenced in the
receipt).

Structure b_ij = p_i * ub_ij (ub rows sum to 1; p_i = max per-capita division
rate, B1/B2), c_ij = chi * uc_ij with uc = (1-a)*ub + a*(1/m): 'a' is the
benefit-crowding pattern asymmetry (a=0 aligned self-pMHC-like crowding; a>0
adds the global IL-7-pool crowding component, Napolitano B9 reading).

Exact-reduction lemmas:
 (L1) at s=0 the equilibrium set depends on (ub, uc, rho_i=d_i/p_i, h, K*)
      only; per-clone division scale p_i enters stability alone.
 (L2) K is a pure unit of x at s=0 (x -> (K'/K) x maps the K-cell onto the
      K'-cell, count-preserving) — K-control cell pairs must match exactly.
 (L3) aligned cells (a=0) admit the concave Lyapunov W(x) =
      sum_j Phi_j(y_j) - sum_i (d_i chi/p_i) x_i, so their stable state is
      unique up to degenerate flat faces (the MacArthur/gradient structure of
      the pitch's lemma family; the census verifies it instance by instance).

All parameters are float64 = exact dyadics (arb reads them exactly).
Seed law: SEED_BASE=20260725; per-cell rng = PCG64(sha256(key)) with key
excluding K and the division band (so K-controls and scale pairs share draws).
"""
import hashlib
import itertools
import json
import math

import numpy as np

SEED_BASE = 20260725
N_CLONES = 20
CHI = 0.5  # crowding intensity (dyadic); c rows sum to chi

DIV_BANDS = {  # log-uniform bands, /day (B1/B2 envelope upper halves; L1: only
    # rho=d/p and patterns move the equilibrium set, receipt Sec on coverage)
    "hi": ((5e-4, 1.5e-3), (3.3e-4, 1.0e-3)),  # (CD4 rows 0-9, CD8 rows 10-19)
}
LOSS_MODES = {  # B3 register
    "slice": None,                 # d_i = rho_i p_i, rho log-uniform [0.8,1.2]
    "lo": (1e-4, 3e-4),            # B3 low edge
    "hi": (4.4e-3, 1.32e-2),       # B3 top third (labeled-loss upper bracket)
}


def _rng(key: str):
    h = hashlib.sha256(key.encode()).digest()
    return np.random.default_rng(int.from_bytes(h[:8], "big"))


def build_cell(m, h, w, sp, a, K, loss, div, seed_tag, tag=""):
    """Return cell dict with exact float64 parameter arrays."""
    key = f"g1|{SEED_BASE}|m={m}|h={h}|w={w}|sp={sp}|a={a}|loss={loss}|s={seed_tag}"
    rng = _rng(key)
    n = N_CLONES
    # usage pattern ub: dominant private-ish niche + shared mix
    mix = rng.uniform(0.5, 1.5, (n, m))
    mix /= mix.sum(1, keepdims=True)
    ub = np.zeros((n, m))
    for i in range(n):
        ub[i, i % m] += (1.0 - w)
        ub[i] += w * mix[i]
    if sp:  # row-sparsity: keep top min(3,m) niches per clone, renormalize
        keep = min(3, m)
        for i in range(n):
            idx = np.argsort(ub[i])[::-1][:keep]
            row = np.zeros(m)
            row[idx] = ub[i, idx]
            ub[i] = row / row.sum()
    uc = (1.0 - a) * ub + a / m
    (c4lo, c4hi), (c8lo, c8hi) = DIV_BANDS[div]
    p = np.empty(n)
    p[:10] = np.exp(rng.uniform(math.log(c4lo), math.log(c4hi), 10))
    p[10:] = np.exp(rng.uniform(math.log(c8lo), math.log(c8hi), 10))
    if loss == "slice":
        rho = np.exp(rng.uniform(math.log(0.8), math.log(1.2), n))
        d = rho * p
    else:
        lo, hi = LOSS_MODES[loss]
        d = np.exp(rng.uniform(math.log(lo), math.log(hi), n))
    b = p[:, None] * ub
    c = CHI * uc
    cid = f"m{m:02d}_h{h}_w{w}_sp{sp}_a{a}_K{K}_L{loss}_D{div}_{seed_tag}{tag}"
    return {"cell_id": cid, "m": m, "h": h, "w": w, "sp": sp, "a": a, "K": K,
            "loss": loss, "div": div, "seed_tag": seed_tag, "n": n,
            "b": b, "c": c, "d": d, "p": p}


def build_grid():
    """The pinned Gate-1 grid: 100 cells at m=8 (>=64 bar) + 40 m-sweep cells.
    K-controls: 4 cells duplicating two m=8 slice cells at K in {0.1, 10}
    (L2 exact-copy law: counts must match their K=1 partners)."""
    cells = []
    for h in (1, 2):
        for w in (0.1, 1.0):
            for sp in (0, 1):
                for a in (0.0, 0.5):
                    for st in ("A", "B"):
                        for loss in ("slice", "lo", "hi"):
                            cells.append(build_cell(8, h, w, sp, a, 1.0,
                                                    loss, "hi", st))
    for K in (0.1, 10.0):  # K-control partners of two slice cells (seed A)
        cells.append(build_cell(8, 1, 0.1, 0, 0.5, K, "slice", "hi", "A", "_KC"))
        cells.append(build_cell(8, 2, 1.0, 0, 0.0, K, "slice", "hi", "A", "_KC"))
    for m in (2, 4, 6, 10, 12):
        for h in (1, 2):
            for w in (0.1, 1.0):
                for a in (0.0, 0.5):
                    cells.append(build_cell(m, h, w, 0, a, 1.0,
                                            "slice", "hi", "A"))
    return cells


def build_plants():
    """Planted instrument controls (bench-only, pilot pattern):
    PLANT_BISTABLE: n=2, m=2, crossed benefit/crowding patterns — exactly two
    stable singleton states (mutual inhibition; analytic truth).
    PLANT_MAC: aligned m=8 seed-A slice cell (L3: must return one stable
    coexistence state)."""
    n, m = 2, 2
    p = np.array([1e-3, 1e-3])
    b = p[:, None] * np.array([[1.0, 0.0], [0.0, 1.0]])
    # each clone crowds the OTHER's niche 3x harder than its own =>
    # two stable singletons + interior saddle (analytic truth, receipted):
    # singleton x*=32, invader margin G_l = p/13 - d < 0 certified.
    c = np.array([[0.125, 0.375], [0.375, 0.125]])
    d = np.array([2e-4, 2e-4])
    plant1 = {"cell_id": "PLANT_BISTABLE", "m": m, "h": 1, "w": 0.0, "sp": 0,
              "a": -1, "K": 1.0, "loss": "plant", "div": "plant",
              "seed_tag": "P", "n": n, "b": b, "c": c, "d": d, "p": p,
              "truth_stable": 2}
    plant2 = build_cell(8, 1, 1.0, 0, 0.0, 1.0, "slice", "hi", "A", "_PLM")
    plant2["cell_id"] = "PLANT_MAC"
    plant2["truth_stable"] = 1
    return [plant1, plant2]


# ---------------------------------------------------------------- float layer
def f_phi(y, K, h):
    t = y / K
    return 1.0 / (1.0 + t) if h == 1 else 1.0 / (1.0 + t * t)


def f_dphi(y, K, h):
    t = y / K
    if h == 1:
        f = 1.0 / (1.0 + t)
        return -(f * f) / K
    f = 1.0 / (1.0 + t * t)
    return -(2.0 * t / K) * (f * f)


def f_G(cell, S, x):
    y = x @ cell["c"][S]
    ph = f_phi(y, cell["K"], cell["h"])
    return cell["b"][S] @ ph - cell["d"][S]


def f_DG(cell, S, x):
    y = x @ cell["c"][S]
    dp = f_dphi(y, cell["K"], cell["h"])
    return cell["b"][S] @ (dp[:, None] * cell["c"][S].T)


def f_xmax(cell, i, hint=None):
    """Float upper bound on x_i at any equilibrium (certified redone in arb)."""
    b, c, d = cell["b"][i], cell["c"][i], cell["d"][i]

    def h_i(t):
        return b @ f_phi(c * t, cell["K"], cell["h"]) - d
    hi = hint or 1.0
    for _ in range(400):
        if h_i(hi) < 0:
            break
        hi *= 2.0
    else:
        return None
    lo = 0.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if h_i(mid) < 0:
            hi = mid
        else:
            lo = mid
    return hi


def float_roots(cell, S, xmax_S, key):
    """Deterministic multi-start Newton; seeds for the certified layer only."""
    k = len(S)
    xm = np.asarray(xmax_S)
    rng = _rng(key)
    starts = [0.5 * xm, 0.15 * xm, 0.85 * xm]
    starts += [xm * rng.uniform(0.05, 0.95, k) for _ in range(2)]
    roots = []
    scale = np.maximum(cell["d"][S], 1e-12)
    for x0 in starts:
        x = x0.copy()
        ok = False
        for _ in range(80):
            g = f_G(cell, S, x)
            if np.max(np.abs(g) / scale) < 1e-12:
                ok = True
                break
            try:
                step = np.linalg.solve(f_DG(cell, S, x), g)
            except np.linalg.LinAlgError:
                break
            t = 1.0
            while np.any(x - t * step <= 0) and t > 1e-12:
                t *= 0.5
            if t <= 1e-12:
                break
            x = x - t * step
        if ok and np.all(x > 0) and np.all(x < 1.05 * xm):
            if not any(np.max(np.abs(x - r) / xm) < 1e-5 for r in roots):
                roots.append(x)
    return roots


def cells_manifest(cells):
    rows = []
    for c in cells:
        hsh = hashlib.sha256()
        for a in ("b", "c", "d", "p"):
            hsh.update(np.ascontiguousarray(c[a]).tobytes())
        rows.append({"cell_id": c["cell_id"], "m": c["m"], "h": c["h"],
                     "w": c["w"], "sp": c["sp"], "a": c["a"], "K": c["K"],
                     "loss": c["loss"], "div": c["div"], "seed": c["seed_tag"],
                     "param_sha256": hsh.hexdigest()[:16]})
    return rows


if __name__ == "__main__":
    cs = build_grid()
    man = cells_manifest(cs)
    print(json.dumps({"n_cells": len(cs),
                      "n_m8": sum(1 for c in cs if c["m"] == 8),
                      "first": man[0], "last": man[-1]}, indent=1))
