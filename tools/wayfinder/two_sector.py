#!/usr/bin/env python3
r"""
two_sector.py — two-branch (s^n, s^{n-eps}) Frobenius series at a point with
APPARENT higher poles (entries up to s^-pmax), fixed eps.

PIECES (package-coherent — none depend on an eps-Laurent-window
representation):

  * scc_blocks            (pure graph code;
                           reverse-topological SCC blocks, dependencies first)
  * block_delays          (Bellman fixed
                           point phi_i = max(0, max_j(phi_j + p_ij - 1));
                           raises on a positive pole cycle, i.e. within-block
                           upward coupling must be one-directional)
  * apparent-pole s^-k entry tables  (the
                           per-entry local Laurent coefficients come from one
                           Cauchy circle on the pointwise DESystem contract,
                           the same shared machinery as frobenius.py, with an
                           on-circle significance test for the pole order and
                           a pole-order > pmax tripwire)
  * two-branch slot schedule  (slot
                           tau = n + phi_i; per-slot small dense solve; the
                           level n <= 0 rows are CONSTRAINTS, evaluated and
                           REPORTED, never fudged)
  * single-pin autopin    (solve_sector_autopin, k=1 — resonant
                           coordinate DETERMINED by the constraint rows via a
                           particular + homogeneous pair, t = -r_p/r_h on the
                           best-coupled row; asserts k == 1, same as the reference driver)

NOT ported (driver-specific; a copy would be dishonest here):
  * the eps-Laurent WINDOW arithmetic (engine._lmul/_ladd/_linv flat windows)
    and lsolve_small's eps-offset resonance pivoting: this package is
    FIXED-eps by design (the windowed Laurent transport is measured-UNSOUND on
    a pointwise basis — see the no-Laurent-window design rule in
    README.md). At fixed eps the coefficients are mpc scalars; resonance shows up
    as a numerically singular slot matrix (SVD-detected) and must be pinned
    (explicitly or via autopin) — a weaker but honest discipline, reported.
  * the original driver's parallel sources service and its own JSON loading
    (de_load owns loading; systems here are any wayfinder DESystem).

MATH
====
Ansatz (region-analysis-derived in the reference driver; generic here):
    M(s) = sum_n a_n s^n + (-s)^{-eps} sum_n b_n s^n,  branch sigma in {0,eps}
with the DE dM/ds = A(s) M and A_ij(s) = sum_{k >= -p_ij} A_ij^{(k)} s^k
(p_ij <= pmax APPARENT poles allowed). Matching s^{n-1-sigma}:
    (n - sigma) x_n = sum_j sum_k A_ij^{(k)} x_{j, n-1-k}
n >= 1 rows are solved (slot schedule closes them despite the k < -1 terms);
rows with n <= 0 are constraints on the seeds, reported as residuals.
References x_{j, l} with l > NL are truncated to 0 — the top ~pmax levels of
downstream components carry truncation-order error, same as the reference driver's
engine's SRC windows; keep NL comfortably above the order you consume.

Import-dps discipline: all parsing inside mp.workdps; global mp.dps never
touched. No sympy.
"""

from mpmath import mp, mpf, mpc

try:
    from .transport import (_to_mpc, _circle_tables, _sample_circle, _fft,
                            _next_pow2, _declared_sings)
except ImportError:  # flat mode (tests, scripts)
    from transport import (_to_mpc, _circle_tables, _sample_circle, _fft,
                           _next_pow2, _declared_sings)

__all__ = ["scc_blocks", "block_delays", "two_sector_series",
           "eval_two_sector"]


class _ResonantSlot(AssertionError):
    """Numerically singular slot system — the fixed-eps face of the reference driver's
    'RESONANT slot needs pin'. .coords = [(component, level), ...]."""

    def __init__(self, msg, coords):
        super().__init__(msg)
        self.coords = coords


# ---------------------------------------------------------------------------
# graph structure (VERBATIM ports from the reference driver)
# ---------------------------------------------------------------------------

def scc_blocks(entry_keys, n):
    """SCC blocks of the dependency graph (i depends on j for each stored
    (i, j)), in reverse-topological order (dependencies first).
    VERBATIM port of the reference driver's scc_blocks (input generalized from the
    dict to any iterable of (i, j) keys)."""
    g = {v: set() for v in range(n)}
    for (i, j) in entry_keys:
        if i != j:
            g[i].add(j)   # i depends on j
    idx_c = [0]
    stack, low, idx, onst, out = [], {}, {}, {}, []

    def strong(v):
        idx[v] = low[v] = idx_c[0]
        idx_c[0] += 1
        stack.append(v)
        onst[v] = True
        for w in g[v]:
            if w not in idx:
                strong(w)
                low[v] = min(low[v], low[w])
            elif onst.get(w):
                low[v] = min(low[v], idx[w])
        if low[v] == idx[v]:
            comp = []
            while True:
                w = stack.pop()
                onst[w] = False
                comp.append(w)
                if w == v:
                    break
            out.append(sorted(comp))
    for v in range(n):
        if v not in idx:
            strong(v)
    return out    # already reverse-topological (dependencies first)


def block_delays(block, poles):
    """Slot delays phi_i for one block: Bellman fixed point
    phi_i = max(0, max_j (phi_j + p_ij - 1)). VERBATIM port of
    the reference driver's block_delays; raises on a positive pole cycle (the certified
    reference condition: within-block upward coupling one-directional)."""
    phi = {i: 0 for i in block}
    for _ in range(len(block) + 2):
        ch = False
        for i in block:
            for j in block:
                p = poles.get((i, j))
                if p is None:
                    continue
                need = phi[j] + p - 1
                if need > phi[i]:
                    phi[i] = need
                    ch = True
        if not ch:
            return phi
    raise AssertionError("positive pole cycle in block %s" % block)


# ---------------------------------------------------------------------------
# local Laurent tables from the pointwise DESystem contract
# ---------------------------------------------------------------------------

def _laurent_tables(desys, x_sing, eps_v, NL, pmax, wp):
    """dict (i,j) -> dict k -> mpc coefficient (k = -p_ij..NL) of A about
    x_sing, from ONE Cauchy circle (shared transport machinery; the role of
    engine.entry_slaurent). Pole order per entry by an on-circle significance
    test; tripwire if anything significant sits below s^-pmax."""
    n = desys.n
    sings = _declared_sings(desys, wp)
    others = [abs(s - x_sing) for s in sings
              if abs(s - x_sing) > mpf(10) ** (-(wp // 2))]
    r = mpf(1) / 8
    if others:
        r = min(r, min(others) / 2)
    N = _next_pow2(max(2 * (NL + pmax + 8), int(3.5 * wp) + 16))
    _, wroots_inv = _circle_tables(N, wp)
    samples = _sample_circle(desys, x_sing, eps_v, r, N, wp)
    tol = mpf(10) ** (-(wp - 20))
    tab = {}
    poles = {}
    for i in range(n):
        for j in range(n):
            col = [mpc(samples[k][i][j]) for k in range(N)]
            scale_e = max(abs(v) for v in col)
            if scale_e == 0:
                continue
            hat = _fft(col, wroots_inv)
            # tripwire: nothing significant below s^-pmax
            for q in (pmax + 1, pmax + 2):
                cq = hat[(-q) % N] / N * r ** q
                if abs(cq) > scale_e * r ** q * tol:
                    raise ValueError(
                        f"A[{i}][{j}]: significant s^-{q} Laurent "
                        f"coefficient at x_sing (|c|={mp.nstr(abs(cq), 3)}) "
                        f"— pole order exceeds pmax={pmax}; raise pmax or "
                        f"check the point")
            ent = {}
            p_ord = 0
            for q in range(pmax, 0, -1):
                cq = hat[(-q) % N] / N * r ** q
                if abs(cq) > scale_e * r ** q * tol:
                    ent[-q] = cq
                    if p_ord == 0:
                        p_ord = q
            rp = mpf(1)
            for m in range(NL + 1):
                cm = hat[m] / N * rp
                if cm != 0:
                    ent[m] = cm
                rp /= r
            if ent:
                tab[(i, j)] = ent
                # negative p (leading zero at x_sing) relaxes the delays,
                # exactly as the reference driver's split_entry qd-qn < 0 did
                if p_ord == 0:
                    p_ord = -min(ent)   # all kept keys >= 0 here
                poles[(i, j)] = p_ord
    return tab, poles, r


# ---------------------------------------------------------------------------
# one branch: slot-scheduled recursion
# ---------------------------------------------------------------------------

def _solve_branch(tab, poles, n, blocks, sigma_v, seed, NL, wp, pinned):
    """Solve one sector (sigma_v = 0 for branch A, eps for branch B) given
    level-0 seeds. pinned: dict (i, level) -> mpc for resonant coordinates.
    Returns (X, resid, sdiag) with X[j] = [x_{j,0}..x_{j,NL}], resid the
    constraint/consistency residual dict {(i, level<=0): mpc}, and sdiag
    honesty diagnostics {'lsq_worst': mpf} (worst least-squares residual of
    pinned rectangular slots — the reference driver recorded its leftover rows too).
    Faithful fixed-eps transcription of the reference driver's LocalSolver.solve_sector
    (simplified: direct RHS convolution instead of SRC/ACC accumulators —
    package clarity over driver throughput; identical equations)."""
    X = [[None] * (NL + 1) for _ in range(n)]
    pmax_eff = max([p for p in poles.values() if p > 0], default=0)
    sdiag = {"lsq_worst": mpf(0)}

    def rhs_row(i, nlev, exclude):
        """sum_j sum_k A_ij^{(k)} x_{j, nlev-1-k} over FINALIZED levels,
        skipping (j, l) in exclude (the slot's own unknowns). References
        beyond NL truncate to 0 (documented)."""
        acc = mpc(0)
        for (ii, j), ent in tab.items():
            if ii != i:
                continue
            for k, c in ent.items():
                l = nlev - 1 - k
                if l < 0 or l > NL or (j, l) in exclude:
                    continue
                v = X[j][l]
                if v is None:
                    raise AssertionError(
                        f"slot schedule violation: row ({i},{nlev}) needs "
                        f"unfinalized x[{j}][{l}] — delays did not close "
                        f"(pole table inconsistent?)")
                if v != 0:
                    acc += c * v
        return acc

    for block in blocks:
        phi = block_delays(block, {ij: p for ij, p in poles.items()
                                   if ij[0] in block and ij[1] in block})
        maxphi = max(phi.values()) if phi else 0
        for j in block:
            X[j][0] = mpc(seed[j])
        for tau in range(1, NL + maxphi + 1):
            slot = [(i, tau - phi[i]) for i in block
                    if 1 <= tau - phi[i] <= NL]
            if not slot:
                continue
            pin_here = [(i, nl) for (i, nl) in slot if (i, nl) in pinned]
            free = [(i, nl) for (i, nl) in slot if (i, nl) not in pinned]
            slotset = set(slot)
            # rows: all slot rows; unknowns: free coordinates
            rows = []
            rhs = []
            for (i, ni) in slot:
                base = rhs_row(i, ni, slotset)
                row = []
                for (jj, nj) in free:
                    coef = mpc(ni) - sigma_v if (jj, nj) == (i, ni) else mpc(0)
                    ent = tab.get((i, jj))
                    if ent is not None:
                        c = ent.get(ni - 1 - nj)
                        if c is not None:
                            coef -= c
                    row.append(coef)
                for (jj, nj) in pin_here:
                    coef = mpc(ni) - sigma_v if (jj, nj) == (i, ni) else mpc(0)
                    ent = tab.get((i, jj))
                    if ent is not None:
                        c = ent.get(ni - 1 - nj)
                        if c is not None:
                            coef -= c
                    base -= coef * pinned[(jj, nj)]
                rows.append(row)
                rhs.append(base)
            if free:
                Amat = mp.matrix(len(slot), len(free))
                for a, row in enumerate(rows):
                    for b, v in enumerate(row):
                        Amat[a, b] = v
                bvec = mp.matrix(len(slot), 1)
                for a, v in enumerate(rhs):
                    bvec[a] = v
                # resonance detection: numerically singular slot system
                sv_m = mp.svd_c(Amat, compute_uv=False)
                svs = [mpf(sv_m[k]) for k in range(sv_m.rows)]
                smax = max(svs) if svs else mpf(0)
                smin = min(svs) if svs else mpf(0)
                if smin <= max(smax, mpf(1)) * mpf(10) ** (-(wp - 25)):
                    # identify the resonant coordinate(s) via the null space
                    reso = list(free)
                    if len(free) > 1:
                        try:
                            U, S, V = mp.svd_c(Amat)
                            nullv = [abs(V[len(free) - 1, b])
                                     for b in range(len(free))]
                            mx = max(nullv)
                            reso = [free[b] for b in range(len(free))
                                    if nullv[b] > mx / 2]
                        except (ValueError, ArithmeticError):
                            pass
                    raise _ResonantSlot(
                        f"RESONANT slot (numerically singular, smin/smax="
                        f"{mp.nstr(smin / max(smax, mpf(1)), 3)}): "
                        f"coordinates {reso} need a pin — pass pins= or "
                        f"autopin=True (reference autopin pattern: "
                        f"solve_sector_autopin)", reso)
                if len(slot) == len(free):
                    xs = mp.lu_solve(Amat, bvec)
                else:
                    xs, lres = mp.qr_solve(Amat, bvec)
                    if mpf(lres) > sdiag["lsq_worst"]:
                        sdiag["lsq_worst"] = mpf(lres)
                for b, (jj, nj) in enumerate(free):
                    X[jj][nj] = mpc(xs[b])
            for (jj, nj) in pin_here:
                X[jj][nj] = mpc(pinned[(jj, nj)])
    # constraint rows n = -pmax+1 .. 0 (evaluated and REPORTED, never fudged)
    resid = {}
    for i in range(n):
        for nlev in range(-pmax_eff + 1, 1):
            lhs = mpc(0)
            if nlev == 0:
                lhs = (mpc(0) - sigma_v) * X[i][0]
            resid[(i, nlev)] = lhs - rhs_row(i, nlev, set())
    return X, resid, sdiag


def _autopin_branch(tab, poles, n, blocks, sigma_v, seed, NL, wp, coord):
    """Fixed-eps scalar transcription of the reference solve_sector_autopin
    (k=1): the resonant coordinate is DETERMINED by the constraint rows via
    particular (pin=0) + homogeneous (zero seed, pin=1) runs and
    t = -r_part/r_hom on the best-coupled constraint row."""
    zero_seed = [mpc(0)] * n
    Xp, rp, dp = _solve_branch(tab, poles, n, blocks, sigma_v, seed, NL, wp,
                               {coord: mpc(0)})
    Xh, rh, dh = _solve_branch(tab, poles, n, blocks, sigma_v, zero_seed, NL,
                               wp, {coord: mpc(1)})
    scale = max([abs(v) for row in Xp for v in row if v is not None]
                + [mpf(1)])
    rows = [k for k, v in rh.items()
            if abs(v) > scale * mpf(10) ** (-(wp - 25))]
    if not rows:
        raise AssertionError(
            "autopin: no constraint row couples to the resonant pin "
            f"{coord} — the coordinate is not determined by the level<=0 "
            f"equations; supply it via pins= (external data needed, same "
            f"contract as the reference driver)")
    best = max(rows, key=lambda k: abs(rh[k]))
    t = -rp[best] / rh[best]
    X = [[(a + t * b if a is not None else None)
          for a, b in zip(Xp[i], Xh[i])] for i in range(n)]
    resid = {k: rp[k] + t * rh.get(k, mpc(0)) for k in rp}
    sdiag = {"lsq_worst": max(dp["lsq_worst"], dh["lsq_worst"])}
    return X, resid, sdiag, t, best


# ---------------------------------------------------------------------------
# public API
# ---------------------------------------------------------------------------

def two_sector_series(desys, eps0, x_sing, seed_a, seed_b, NL, dps,
                      pmax=4, pins_a=None, pins_b=None, autopin=False):
    """Two-branch Frobenius series of y' = A y at x_sing (apparent poles up
    to s^-pmax allowed), fixed eps:

        y(x_sing + s) ~ sum_n a_n s^n + (-s)^{-eps0} sum_n b_n s^n

    Args:
        desys: wayfinder DESystem (pointwise .A contract; declared
            singular points cap the Cauchy radius).
        eps0: FIXED eps (real or complex; the recursion is algebraic in eps).
            Resonance risk grows as eps0 approaches rationals with small
            denominator — the slot SVD test catches it (raise, never silent).
        seed_a, seed_b: level-0 seeds (length n; region-derived — this
            module does NOT derive seeds, same split as the reference driver, where they
            came from the region analysis).
        NL: series depth. Top ~pmax levels of downstream components carry
            truncation-order error (module docstring); keep NL > needed+pmax.
        pins_a/pins_b: dict {(component, level): value} for resonant
            coordinates known externally.
        autopin: if True, a SINGLE resonant coordinate per branch is
            determined from the constraint rows (particular+homogeneous,
            t = -r_p/r_h autopin pattern). >1 resonant coordinate
            raises (generalize when a real case needs it — the reference driver asserted
            the same).

    Returns dict:
        'A', 'B'         list[n] of list[NL+1] mpc coefficients (@dps)
        'residuals_A/B'  {(i, level<=0): mpc} constraint residuals (@dps) —
                         REPORTED, never fudged; gate on them
        'pins_A/B'       pins used (given or autopinned, @dps)
        'autopin_rows'   {'A'/'B': constraint row used} when autopinned
        'blocks', 'poles', 'r_circle', 'x_sing', 'eps0', 'NL', 'dps'
    """
    wp = dps + 30
    with mp.workdps(wp):
        eps_v = _to_mpc(eps0)
        xs = _to_mpc(x_sing)
        n = desys.n
        if len(seed_a) != n or len(seed_b) != n:
            raise ValueError(f"seeds must have length n={n}")
        sa = [_to_mpc(v) for v in seed_a]
        sb = [_to_mpc(v) for v in seed_b]
        tab, poles, r = _laurent_tables(desys, xs, eps_v, NL, pmax, wp)
        blocks = scc_blocks(tab.keys(), n)
        out = {"blocks": blocks,
               "poles": {f"{i},{j}": p for (i, j), p in poles.items()},
               "autopin_rows": {}}
        for name, sigma_v, seed, pins in (
                ("A", mpc(0), sa, dict(pins_a or {})),
                ("B", eps_v, sb, dict(pins_b or {}))):
            pins_v = {k: _to_mpc(v) for k, v in pins.items()}
            try:
                X, resid, sdiag = _solve_branch(tab, poles, n, blocks,
                                                sigma_v, seed, NL, wp, pins_v)
            except _ResonantSlot as exc:
                if not autopin:
                    raise
                if pins_v:
                    raise AssertionError(
                        f"branch {name}: resonant slot {exc.coords} on top "
                        f"of explicit pins — mixed mode unsupported; pin "
                        f"everything explicitly") from exc
                if len(exc.coords) != 1:
                    raise AssertionError(
                        f"branch {name}: {len(exc.coords)} resonant "
                        f"coordinates {exc.coords} — autopin implements "
                        f"k=1 only (the reference driver asserted the same; generalize "
                        f"when a real case needs it)") from exc
                X, resid, sdiag, t, row_used = _autopin_branch(
                    tab, poles, n, blocks, sigma_v, seed, NL, wp,
                    exc.coords[0])
                pins_v = {exc.coords[0]: t}
                out["autopin_rows"][name] = row_used
            with mp.workdps(dps):
                out[name] = [[+mpc(v) if v is not None else mpc(0)
                              for v in X[i]] for i in range(n)]
                out["residuals_" + name] = {k: +mpc(v)
                                            for k, v in resid.items()}
                out["pins_" + name] = {k: +mpc(v)
                                       for k, v in pins_v.items()}
                out["slot_lsq_worst_" + name] = +sdiag["lsq_worst"]
        with mp.workdps(dps):
            out.update({"r_circle": +r, "x_sing": +xs, "eps0": +eps_v,
                        "NL": NL, "dps": dps})
    return out


def eval_two_sector(result, s, dps=None):
    """Evaluate A(s) + (-s)^{-eps0} B(s) from a two_sector_series result
    (principal branch of log(-s), matching frobenius.py / the reference driver
    minus_s_pow_minus_eps). Truncation is the caller's business — check the
    last retained coefficients like any Frobenius series."""
    dps = dps or result["dps"]
    wp = dps + 15
    with mp.workdps(wp):
        sv = _to_mpc(s)
        eps_v = mpc(result["eps0"])
        pref = mp.exp(-eps_v * mp.log(-sv))
        vals = []
        for Ai, Bi in zip(result["A"], result["B"]):
            va = mpc(0)
            for c in reversed(Ai):
                va = va * sv + c
            vb = mpc(0)
            for c in reversed(Bi):
                vb = vb * sv + c
            vals.append(va + pref * vb)
    with mp.workdps(dps):
        return [+v for v in vals]
