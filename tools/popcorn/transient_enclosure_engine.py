#!/usr/bin/env python3
"""transient_enclosure_engine — certified enclosures of the transient selected SFS.

WHAT IT COMPUTES. Rigorous two-sided interval enclosures of the expected
unfolded sample site-frequency spectrum E_n(i; t), i = 1..n-1 (theta = 1
units), under genic (additive, h = 1/2) selection of scaled strength S <= 0
and a PIECEWISE-CONSTANT population-size history that follows an ancestral
Wright equilibrium at the same S: the transient (non-equilibrium) selected
spectrum. The moment system of the Poisson-random-field forward diffusion
is truncated at M >= n with a proved two-sided tail closure, integrated
through the epochs in ball arithmetic (Arb, via python-flint: rigorous
matrix exponential and linear solve), and projected to sample size n with
exact integer coefficients. Every reported entry is a ball [mid +/- rad]
guaranteed to contain the value of the model below; the width certificate
is self-reporting, so a truncation M that is too small shows up as a wide
ball, never as a silently wrong number.

REGISTER: CERTIFIED-ENCLOSURE, MODEL-CONDITIONAL. The enclosure is rigorous
for the stated diffusion model, its moment hierarchy and the stated
truncation brackets; it is not a bound on model error (it does not certify
the diffusion against a discrete Wright-Fisher chain, and it knows nothing
about linkage or dominance). Three other registers live in this file and
must never be mixed with it in one sentence:
  FLOAT      float_cell(): scipy dense expm of the same truncated LOWER
             system in double precision. An implementation check only:
             at the shipped cells it deviates from the enclosure midpoint
             by up to 3e-8 (relative), i.e. some 60 orders of magnitude
             outside the certified width and inside its own rounding
             error, so the certified value adjudicates the float, not the
             reverse.
  FLOAT-HP   xcheck(): an independent mpmath vector-Taylor substep
             integrator of the LOWER system (different arithmetic, different
             algorithm, different code path) compared with the ball route;
             an implementation cross-check. Its worst_rel_dev is |mpmath
             LOWER value - hull midpoint| / |midpoint|; at the shipped
             instance (M = 40) that equals the half LOWER/UPPER bracket at
             j = 20 (1.5e-32, a truncation quantity, hence substep-
             independent), while against the LOWER-bracket midpoint the
             mpmath route agrees to 2e-50 (its working precision).
  FIELD-FLOAT  moments_cell(): the moment-closure solver `moments`
             (Jouganous et al. 2017), an OPTIONAL dependency imported only
             inside that function, run at n directly and at n_big projected
             down; run_cell() reports its deviation from the enclosure when
             it is installed and says so by name when it is not.

MODEL AND CONVENTIONS (the certificate's conditionality). Density f(x, t)
of derived-allele frequencies under the infinite-sites forward diffusion
with genic selection and piecewise-constant size,

    df/dt = (1/(2 rho(t))) d^2/dx^2 [x(1-x) f] - gamma d/dx [x(1-x) f]

on (0, 1), with S = 4 N_ref s (S < 0 deleterious; the same S as
popcorn.sfs; dadi/moments gamma = S/2), gamma = S/2, time t in units of
2 N_ref generations, rho(t) = N(t)/N_ref piecewise constant, and new
mutations entering at x = 0+ at rate theta/2 per unit time with
theta = 4 N_ref u = 1 (spectra scale linearly in theta). Regularity assumed
of the model: f >= 0 with finite moments, x f(x, t) -> rho(t) theta as
x -> 0 (the mutational baseline) and (1-x) f -> 0 as x -> 1. The history
`epochs` = [(rho_1, T_1), (rho_2, T_2), ...] runs past -> present AFTER an
ancestral equilibrium epoch at rho = 1 (ANC_RHO); each rho_k, T_k is an
exact rational (flint.fmpq). The module constant EPOCHS is the reference
history of the shipped cells: rho = 1/5 for T = 1/10, then rho = 1 for
T = 1/20.

THE MOMENT SYSTEM. With w_j(t) = int_0^1 x^j (1-x) f(x, t) dx, j >= 1, the
test functions x^j (1-x) close under the generator (both boundary terms
vanish at x = 1; at x = 0 only j = 1 picks up the influx, at the
N-independent rate theta/2):

    dw_j/dt = [ j(j-1) w_{j-1} - j(j+1) w_j ] / (2 rho(t))
              + gamma [ j w_j - (j+1) w_{j+1} ] + (theta/2) delta_{j,1}.

Tridiagonal in structure: drift couples down (sub-diagonal), selection
couples up (super-diagonal). On each epoch the system is linear
time-invariant: w(t_k + T) = e^{A_k T} (w(t_k) - p_k) + p_k with
A_k p_k = -b, b = (theta/2) e_1. This is the class of moment equations for
the non-equilibrium frequency spectrum of Evans, Shvets & Slatkin (2007).
Equilibrium (initial condition, ancestral rho = 1), by the 1F1 integral
representation:

    w_j^eq = theta [1 - 1F1(1; j+1; -S)] / (j (1 - e^{-S})),
    S = 0:  w_j^eq = theta / (j (j+1)).

Projection to the sample (exact integer coefficients), i = 1..n-1:

    E_n(i) = C(n, i) sum_{k=0}^{n-i-1} C(n-i-1, k) (-1)^k w_{i+k},

which needs only w_1..w_{n-1}; the alternating sum amplifies ball width by
up to ~2^(n-2), a precision-padding rule (below), not an obstacle. The
stationary sample entry it must reproduce at t = 0 is the closed form
E_n(i; S) = n/(i(n-i)) (1 - 1F1(n-i; n; -S))/(1 - e^{-S}) (S = 0: 1/i), the
same object popcorn.sfs certifies (equilibrium_entry()).

THE TRUNCATION-TAIL CLOSURE (LEMMA; the enclosure's teeth). Truncate at
M >= n; only row M references the discarded tail, through
-gamma (M+1) w_{M+1}. Fix gamma <= 0. Let A_lo be the M x M truncated
matrix with w_{M+1} := 0 and A_hi = A_lo + |gamma| (M+1) e_M e_M^T (i.e.
w_{M+1} := w_M). Let v, W solve the corresponding affine systems with
v(0) <= w(0) <= W(0) entrywise, w the true moments. Then for all t >= 0
and j = 1..M:  v_j(t) <= w_j(t) <= W_j(t).
  Proof. (i) 0 <= w_{M+1}(t) <= w_M(t) exactly, since
  0 <= x^{M+1} (1-x) <= x^M (1-x) on [0, 1] and f >= 0. (ii) A_lo and
  A_hi are Metzler (off-diagonal entries j(j-1)/(2 rho) >= 0 and
  -gamma (j+1) >= 0 for gamma <= 0), so their propagators are entrywise
  nonnegative (e^{At} = lim_k (I + tA/k)^k with I + tA/k >= 0 for k
  large). (iii) z = (w_1..w_M) - v satisfies z' = A_lo z + g(t) e_M with
  g = |gamma| (M+1) w_{M+1} >= 0 and z(0) >= 0, so by Duhamel
  z(t) = Phi_lo(t, 0) z(0) + int_0^t Phi_lo(t, s) g(s) e_M ds >= 0
  (composed across epochs, every factor entrywise nonnegative). Likewise
  W - (w_1..w_M) is driven by |gamma| (M+1) (w_M - w_{M+1}) >= 0 through
  A_hi from nonnegative initial data, hence stays >= 0. QED.
  Corollaries. (a) S = 0: the super-diagonal vanishes, the hierarchy is
  lower triangular and truncation is EXACT (v = W = w up to arithmetic;
  the UPPER pass is skipped). (b) The lemma composes across epochs and
  accepts interval initial data, so ball-valued equilibrium initial
  conditions propagate rigorously. (c) Beneficial S > 0 is NOT covered
  (A is no longer Metzler); certified_cell() refuses S > 0 (ValueError).
  Tightness (informal, not load-bearing; the enclosure is the bracket
  regardless): the closure gap forces row M only and reaches w_j only
  through selection entries against drift damping, one level per factor
  ~ 2 rho |gamma| / j, so it is factorially small for M >> 2 rho |gamma|.
  Measured on the reference history at S = -5, n = 20, prec = 256: the
  relative width (2 rad/|mid|) of the sample entries is 5e-15 at M = 40
  (truncation-dominated: the width IS the bracket), 9e-37 at M = 60,
  4e-61 at M = 80 and ~4e-69 from M = 100 on (precision floor); at
  M = 200 the LOWER/UPPER separation is ~7e-8 at j = M and below the
  ball radius (~1e-77) at j <= 20.

THE CERTIFIED PIPELINE (certified_cell).
  1. Exact rational matrices: A_lo, A_hi and A*T entries as fmpq (S, rho, T
     rational); one rounding into arb at `prec` bits.
  2. Initial condition: w_j^eq balls via arb hypgeom_1f1, j = 1..M (the
     exact infinite-system equilibrium; no truncation enters).
  3. Per epoch, both brackets: E = (A T).exp() (rigorous arb_mat matrix
     exponential); p = A.solve(-b) (rigorous); w <- p + E (w - p).
  4. Hull per moment: the union ball of the LOWER and UPPER results
     (contains [inf v, sup W], hence the true w_j).
  5. Projection with exact integer coefficients -> enclosure balls E_n(i).
  Two certified checks run inside every cell and are reported under
  rec["gates"]: 'stationarity_interior_contains_zero' (the ball residual
  A_anc w^eq + b contains 0 in every interior row: rows 1..M-1 of the
  truncated system coincide with the infinite system, and w^eq is its
  exact equilibrium; catches a wrong matrix, initial condition or forcing)
  and 'eq_projection_overlaps_closed_form' (the projected initial balls
  overlap the closed-form stationary sample SFS for every i; ties this
  engine to the stationary object of popcorn.sfs).
  Precision rule (measured): prec >= 64 + 3.33 (0.30 n + D) bits for D
  target digits (projection amplification plus scaling-and-squaring loss);
  prec = 256 at n = 20 lands 60+ digit enclosures. Sizing rule: validity
  needs only M >= n; tightness wants M >> 2 rho_max |gamma|; run an M
  ladder and read the widths. Cost is dominated by the matrix exponentials
  (~M^3 each, 2 per epoch per bracket): single-core, prec 256, two epochs,
  both brackets: ~0.1 s at M = 40, ~1 s at M = 80, ~25 s at M = 200,
  ~6 min at M = 500.

REFERENCE DATA (reference/enclosure/, all at n = 20, prec = 256 on EPOCHS):
  cell_M{200,500}_S{0,-5,-50}_n20.json  six certified cells (checks, moment
      diagnostics, per-entry mid to 30 digits + radius + relative width,
      and the FLOAT route's values under 'float_scipy'). For each S the
      M = 200 and M = 500 enclosures agree to all 30 recorded digits (two
      truncations, one value); relative widths 1e-76 .. 1e-64.
  xcheck_M40_S-5.json  the FLOAT-HP cross-check instance (M = 40, S = -5,
      first reference epoch, 512 substeps, 50 digits): worst relative
      deviation 1.55e-32 on w_1..w_20 (= the half LOWER/UPPER bracket at
      j = 20 for M = 40; the mpmath-vs-LOWER agreement itself is 2e-50).
  load_reference_cell(M, S), load_reference_xcheck() and entry_ball(entry)
  read them back (entry_ball rebuilds an arb ball from 'mid30' and 'rad').

SIDE EFFECTS AND TYPES. certified_cell() and xcheck() set the process-wide
flint precision ctx.prec = prec and leave it there (the returned balls are
meant to be used at that precision); reset it yourself when mixing with
other flint code. S is a non-positive integer (an fmpq rational S is
accepted by certified_cell but not by float_cell/moments_cell or the CLI);
epoch entries are fmpq. Importing this module imports python-flint;
numpy/scipy (float_cell), mpmath (float_cell, xcheck) and moments
(moments_cell) are imported inside the functions that need them.

CLI (prints; writes JSON only where told):
  transient_enclosure_engine.py cell   --M 200 --S -5 [--n 20 --prec 256] [--json OUT.json]
  transient_enclosure_engine.py xcheck [--outdir DIR]      # FLOAT-HP cross-check
  transient_enclosure_engine.py all    --outdir DIR        # the six reference cells (~15 min)
`cell` and `all` also run the FLOAT route and, if installed, the moments
cross-check, and report containment/deviation per entry.

REFERENCES.
  Kimura M. (1964) Diffusion models in population genetics. J. Appl.
    Probab. 1:177-232. doi:10.2307/3211856
  Sawyer S.A., Hartl D.L. (1992) Population genetics of polymorphism and
    divergence. Genetics 132:1161-1176. doi:10.1093/genetics/132.4.1161
  Evans S.N., Shvets Y., Slatkin M. (2007) Non-equilibrium theory of the
    allele frequency spectrum. Theor. Popul. Biol. 71:109-119.
    doi:10.1016/j.tpb.2006.06.005  (moment equations for the transient SFS)
  Zivkovic D., Steinrucken M., Song Y.S., Stephan W. (2015) Transition
    densities and sample frequency spectra of diffusion processes with
    selection and variable population size. Genetics 200:601-617.
    doi:10.1534/genetics.115.175265  (spectral route to the same object)
  Jouganous J., Long W., Ragsdale A.P., Gravel S. (2017) Inferring the joint
    demographic history of multiple populations: beyond the diffusion
    approximation. Genetics 206:1549-1567. doi:10.1534/genetics.117.200493
    (the jackknife moment closure of `moments`, the optional cross-check)
  Johansson F. (2017) Arb: efficient arbitrary-precision midpoint-radius
    interval arithmetic. IEEE Trans. Comput. 66:1281-1292.
    doi:10.1109/TC.2017.2690633  (the ball arithmetic behind python-flint)
"""
import argparse, json, os, sys, time
from datetime import datetime, timezone
from fractions import Fraction
from math import comb

from flint import arb, arb_mat, ctx, fmpq

# ---------------------------------------------------------------- conventions
THETA = fmpq(1)                       # theta = 4 N_ref u = 1
ANC_RHO = fmpq(1)                     # ancestral size = N_ref (equilibrium IC)
# reference history of the shipped cells (past -> present, after the
# ancestral equilibrium at rho=1):  rho=1/5 for T=1/10  ->  rho=1 for T=1/20
EPOCHS = [(fmpq(1, 5), fmpq(1, 10)), (fmpq(1), fmpq(1, 20))]

HERE = os.path.dirname(os.path.abspath(__file__))
REFERENCE_DIR = os.path.join(HERE, "reference", "enclosure")
REFERENCE_CELLS = tuple((M, S) for M in (200, 500) for S in (0, -5, -50))   # n = 20


def qarb(q):
    """fmpq/int/Fraction -> arb (ball rounding at current prec, tracked)."""
    if isinstance(q, Fraction):
        q = fmpq(q.numerator, q.denominator)
    q = fmpq(q)
    return arb(q.p) / arb(q.q)


# ---------------------------------------------------- exact rational matrices
def rat_rows(M, S, rho, dt=None, upper=False):
    """Exact fmpq rows of A (or A*dt if dt given) for the truncated system.
    Row r=0..M-1 is moment j=r+1.  upper=True adds the rank-one tail term
    -gamma*(M+1) at (M-1,M-1)  [w_{M+1} := w_M closure]."""
    g = fmpq(S) / 2
    rows = []
    for r in range(M):
        j = r + 1
        row = [fmpq(0)] * M
        if j >= 2:
            row[r - 1] = fmpq(j * (j - 1)) / (2 * rho)
        row[r] = fmpq(-j * (j + 1)) / (2 * rho) + g * j
        if j < M:
            row[r + 1] = -g * (j + 1)
        rows.append(row)
    if upper:
        rows[M - 1][M - 1] += -g * (M + 1)
    if dt is not None:
        rows = [[e * dt for e in row] for row in rows]
    return rows


def to_arb_mat(rat):
    return arb_mat([[qarb(e) for e in row] for row in rat])


def bvec(M):
    """forcing b = (theta/2) e_1 as arb column."""
    col = [[arb(0)] for _ in range(M)]
    col[0][0] = qarb(THETA / 2)
    return arb_mat(col)


# ------------------------------------------------------- equilibrium moments
def weq_ball(j, S):
    """w_j at Wright equilibrium (rho=1), theta=1, as an arb ball.
    w_j = 1/(j(1-e^{-S})) * (1 - 1F1(1; j+1; -S));  S=0: 1/(j(j+1)).
    (closed form via the 1F1 integral representation; certified by arb)"""
    if S == 0:
        return qarb(fmpq(1, j * (j + 1)))
    z = arb(-S)                                    # exact integer here
    one_f1 = z.hypgeom_1f1(arb(1), arb(j + 1))    # 1F1(1; j+1; -S)
    den = arb(j) * (-(arb(-S).expm1()))           # j * (1 - e^{-S})
    return (arb(1) - one_f1) / den


def equilibrium_entry(n, i, S):
    """Certified stationary (Wright-equilibrium, rho=1) sample-SFS entry as
    an arb ball, theta=1: E_n(i;S) = n/(i(n-i)) (1 - 1F1(n-i; n; -S)) /
    (1 - e^{-S}); S=0 -> 1/i. The same closed form popcorn.sfs certifies."""
    if S == 0:
        return qarb(fmpq(1, i))
    z = arb(-S)
    num = arb(1) - z.hypgeom_1f1(arb(n - i), arb(n))
    den = -(arb(-S).expm1())
    return qarb(fmpq(n, i * (n - i))) * num / den


# ------------------------------------------------------------------ certified
def certified_cell(M, S, n, prec, epochs=EPOCHS, verbose=True):
    """One certified cell: truncation M >= n, integer S <= 0, sample size n,
    `prec` bits, history `epochs` [(fmpq rho, fmpq T), ...] after the
    ancestral rho=1 equilibrium. Sets flint ctx.prec = prec (left set).
    Returns a dict: M, S, n, prec_bits, epochs, theta, anc_rho, walls_s,
    gates (the two certified in-cell checks), moment_diag (LOWER/UPPER
    separations), enclosure (per i: mid30, rad, rel_width, str), hull (arb
    ball per moment) and entry_balls (arb ball per sample entry E_n(i),
    i = 1..n-1). Register CERTIFIED-ENCLOSURE, model-conditional."""
    if not S <= 0:
        raise ValueError("the two-sided closure lemma requires deleterious S <= 0")
    if not M >= n:
        raise ValueError("need M >= n for the projection")
    ctx.prec = prec
    t0 = time.time()
    rec = {"M": M, "S": S, "n": n, "prec_bits": prec,
           "epochs": [[str(r), str(d)] for r, d in epochs],
           "theta": "1", "anc_rho": "1", "walls_s": {}, "gates": {}}

    # --- initial condition: exact-equilibrium balls (infinite system, no truncation)
    t = time.time()
    w0 = [weq_ball(j, S) for j in range(1, M + 1)]
    rec["walls_s"]["ic_1f1"] = round(time.time() - t, 3)

    # --- check 1: certified stationarity residual of the IC
    A_anc = to_arb_mat(rat_rows(M, S, ANC_RHO))
    w0col = arb_mat([[x] for x in w0])
    res = A_anc * w0col + bvec(M)
    interior_ok = all(res[r, 0].contains(arb(0)) for r in range(M - 1))
    max_rad = max(float(res[r, 0].rad()) for r in range(M - 1))
    rec["gates"]["stationarity_interior_contains_zero"] = bool(interior_ok)
    rec["gates"]["stationarity_max_residual_radius"] = max_rad
    rec["gates"]["row_M_tail_residual_bound"] = str(abs(res[M - 1, 0]).str(5))

    # --- check 2: certified equilibrium projection vs the closed form
    eq_proj_ok = True
    for i in range(1, n):
        e = project_entry(n, i, w0)
        if not e.overlaps(equilibrium_entry(n, i, S)):
            eq_proj_ok = False
    rec["gates"]["eq_projection_overlaps_closed_form"] = bool(eq_proj_ok)

    # --- propagate bracket systems through the epochs
    v = arb_mat([[x] for x in w0])     # LOWER (w_{M+1} := 0)
    W = arb_mat([[x] for x in w0])     # UPPER (w_{M+1} := w_M)
    b = bvec(M)
    for k, (rho, dt) in enumerate(epochs):
        for tag, upper in (("lo", False), ("hi", True)):
            if S == 0 and upper:
                continue
            A = to_arb_mat(rat_rows(M, S, rho, upper=upper))
            Adt = to_arb_mat(rat_rows(M, S, rho, dt=dt, upper=upper))
            t = time.time()
            E = Adt.exp()
            wall_exp = time.time() - t
            t = time.time()
            wp = A.solve(-b)           # particular (equilibrium of truncated affine)
            wall_solve = time.time() - t
            if upper:
                W = wp + E * (W - wp)
            else:
                v = wp + E * (v - wp)
            rec["walls_s"][f"epoch{k}_{tag}_exp"] = round(wall_exp, 3)
            rec["walls_s"][f"epoch{k}_{tag}_solve"] = round(wall_solve, 3)
    if S == 0:
        W = v

    # --- enclosure hulls per moment + diagnostics
    hull = [v[r, 0].union(W[r, 0]) for r in range(M)]
    rec["moment_diag"] = {
        "bracket_sep_mid_j1": (W[0, 0] - v[0, 0]).mid().str(5),
        "bracket_sep_mid_jn": (W[n - 1, 0] - v[n - 1, 0]).mid().str(5),
        "bracket_sep_mid_jM": (W[M - 1, 0] - v[M - 1, 0]).mid().str(5),
        "max_ball_rad_j_le_n": max(float(hull[r].rad()) for r in range(n)),
    }

    # --- projection to sample n
    entries = []
    for i in range(1, n):
        e = project_entry(n, i, hull)
        mid, rad = float(e.mid()), float(e.rad())
        entries.append({"i": i, "mid30": e.mid().str(30), "rad": rad,
                        "rel_width": (2 * rad / abs(mid)) if mid else None,
                        "str": e.str(20, radius=True)})
    rec["enclosure"] = entries
    rec["walls_s"]["total_certified"] = round(time.time() - t0, 3)
    rec["hull"] = hull            # arb objects (in-process use; stripped for JSON)
    rec["entry_balls"] = [project_entry(n, i, hull) for i in range(1, n)]
    if verbose:
        g = rec["gates"]
        print(f"[cell M={M} S={S} prec={prec}] gates: "
              f"stationarity={g['stationarity_interior_contains_zero']} "
              f"eqproj={g['eq_projection_overlaps_closed_form']} | certified wall "
              f"{rec['walls_s']['total_certified']}s")
    return rec


def project_entry(n, i, w):
    """E_n(i) = C(n,i) sum_{k=0}^{n-i-1} C(n-i-1,k) (-1)^k w_{i+k}  (exact ints).
    w: list of arb (index j-1) or arb_mat column."""
    get = (lambda j: w[j - 1, 0]) if isinstance(w, arb_mat) else (lambda j: w[j - 1])
    s = arb(0)
    for k in range(0, n - i):
        c = comb(n - i - 1, k) * (-1 if k % 2 else 1)
        s += arb(c) * get(i + k)
    return arb(comb(n, i)) * s


# --------------------------------------------------------------------- floats
def float_cell(M, S, n, epochs=EPOCHS):
    """FLOAT register: scipy dense expm of the same LOWER-closure system."""
    import numpy as np
    from scipy.linalg import expm, solve as fsolve
    t0 = time.time()

    def fA(rho):
        A = np.zeros((M, M))
        g = S / 2.0
        for r in range(M):
            j = r + 1
            if j >= 2:
                A[r, r - 1] = j * (j - 1) / (2 * rho)
            A[r, r] = -j * (j + 1) / (2 * rho) + g * j
            if j < M:
                A[r, r + 1] = -g * (j + 1)
        return A

    from mpmath import mp, mpf, hyp1f1, expm1
    mp.dps = 30
    if S == 0:
        w = np.array([1.0 / (j * (j + 1)) for j in range(1, M + 1)])
    else:
        w = np.array([float((1 - hyp1f1(1, j + 1, -S)) / (j * (-expm1(mpf(-S)))))
                      for j in range(1, M + 1)])
    b = np.zeros(M); b[0] = 0.5
    for rho, dt in epochs:
        A = fA(int(rho.p) / int(rho.q))
        wp = fsolve(A, -b)
        w = wp + expm(A * (int(dt.p) / int(dt.q))) @ (w - wp)
    out = [float(sum(comb(n - i - 1, k) * (-1) ** k * w[i + k - 1]
                     for k in range(0, n - i)) * comb(n, i)) for i in range(1, n)]
    return {"sfs": out, "wall_s": round(time.time() - t0, 3)}


def moments_cell(S, n, n_big, epochs=EPOCHS):
    """FIELD-FLOAT register: the moment-closure solver `moments` (optional
    dependency, imported here only), run at n directly (default and fine
    time step) and at n_big projected down to n. Integer S."""
    import moments
    g = S / 2.0
    out = {}
    t0 = time.time()
    for lbl, nn, dtf in (("direct", n, 0.02), ("direct_dtfine", n, 0.001),
                         ("projected", n_big, 0.02)):
        fs = moments.Spectrum(moments.LinearSystem_1D.steady_state_1D(nn, gamma=g))
        for rho, dt in epochs:
            fs.integrate([int(rho.p) / int(rho.q)], int(dt.p) / int(dt.q),
                         dt_fac=dtf, gamma=g, h=0.5, theta=1.0)
        if nn != n:
            fs = fs.project([n])
        out[lbl] = [float(fs[i]) for i in range(1, n)]
    out["n_big"] = n_big
    out["wall_s"] = round(time.time() - t0, 3)
    return out


# --------------------------------------------------------- independent xcheck
def xcheck(M=40, S=-5, n=20, prec=256, dps=50, nsub=512, verbose=True):
    """FLOAT-HP independent route: mpmath vector-Taylor substep integrator of
    the LOWER system vs the arb route, single epoch (the first of EPOCHS:
    rho=1/5, T=1/10). Different arithmetic (mpf vs arb balls), different
    algorithm (substepped vector Taylor vs arb_mat scaling-and-squaring exp),
    different code path. Returns dict M, S, nsub, dps, worst_rel_dev (over
    w_1..w_n: |mpmath LOWER value - hull midpoint| / |midpoint|, where the
    hull midpoint sits half the LOWER/UPPER bracket away from LOWER, so at
    small M this reports the half bracket, 1.5e-32 at M = 40, j = 20; the
    mpmath route agrees with the LOWER-bracket midpoint to ~1e-dps) and
    agree_digits. The defaults are the shipped reference instance
    (reference/enclosure/xcheck_M40_S-5.json)."""
    from mpmath import mp, mpf, hyp1f1, expm1, lu_solve, matrix
    ep = [EPOCHS[0]]
    # --- arb route
    rec = certified_cell(M, S, n, prec, epochs=ep, verbose=False)
    v_arb = rec["hull"]
    # --- mpmath route (independent implementation)
    mp.dps = dps
    rho, T = mpf(1) / 5, mpf(1) / 10
    g = mpf(S) / 2
    if S == 0:
        w = [mpf(1) / (j * (j + 1)) for j in range(1, M + 1)]
    else:
        w = [(1 - hyp1f1(1, j + 1, -S)) / (j * (-expm1(mpf(-S)))) for j in range(1, M + 1)]

    def matvec(y):
        out = [mpf(0)] * M
        for r in range(M):
            j = r + 1
            acc = (-j * (j + 1) / (2 * rho) + g * j) * y[r]
            if j >= 2:
                acc += (j * (j - 1) / (2 * rho)) * y[r - 1]
            if j < M:
                acc += (-g * (j + 1)) * y[r + 1]
            out[r] = acc
        return out

    b = [mpf(0)] * M; b[0] = mpf(1) / 2
    h = T / nsub
    for _ in range(nsub):
        # affine step: w <- w + sum_{k>=1} h^k/k! (A^k w + A^{k-1} b)
        term_w = list(w); term_b = list(b)
        acc = list(w)
        k = 1
        while True:
            new_w = matvec(term_w)
            fac = h / k
            term = [fac * (nw + tb) for nw, tb in zip(new_w, term_b)]
            acc = [a + t for a, t in zip(acc, term)]
            term_b = [fac * x for x in matvec(term_b)]
            term_w = [fac * x for x in new_w]
            k += 1
            tn = max(abs(x) for x in term)
            an = max(abs(x) for x in acc)
            if tn < an * mpf(10) ** (-(dps + 8)):
                break
        w = acc
    # --- agreement (digits) on the first n moments, compared AT dps precision
    mp.dps = dps + 10
    worst = mpf(0)
    for j in range(1, n + 1):
        a_hp = mpf(v_arb[j - 1].mid().str(dps + 10, radius=False))
        d = abs(a_hp - w[j - 1]) / abs(a_hp)
        worst = max(worst, d)
    worst = float(worst)
    import math as _m
    digits = -_m.log10(worst) if worst > 0 else dps
    if verbose:
        print(f"[xcheck M={M} S={S} nsub={nsub} dps={dps}] worst rel dev arb-vs-mpmath "
              f"on w_1..w_{n}: {worst:.3e}  (~{digits:.1f} digits)")
    return {"M": M, "S": S, "nsub": nsub, "dps": dps,
            "worst_rel_dev": worst, "agree_digits": digits}


# ------------------------------------------------------------- reference data
def reference_cell_path(M, S, n=20):
    return os.path.join(REFERENCE_DIR, f"cell_M{M}_S{S}_n{n}.json")


def load_reference_cell(M, S, n=20):
    """A shipped certified cell (dict: M, S, n, prec_bits, epochs, theta,
    anc_rho, gates, moment_diag, enclosure[{i, mid30, rad, rel_width}],
    float_scipy {sfs}). (M, S) in REFERENCE_CELLS, n = 20."""
    with open(reference_cell_path(M, S, n)) as f:
        return json.load(f)


def load_reference_xcheck():
    """The shipped FLOAT-HP cross-check instance (dict, xcheck() schema)."""
    with open(os.path.join(REFERENCE_DIR, "xcheck_M40_S-5.json")) as f:
        return json.load(f)


def entry_ball(entry):
    """arb ball from a stored enclosure entry {"mid30": "[m +/- r]", "rad": r}:
    the 30-digit printed midpoint interval widened by the stored radius, so
    it contains the ball that was computed. Uses the current ctx.prec."""
    return arb(entry["mid30"]) + arb(0, entry["rad"])


# ----------------------------------------------------------------------- main
def run_cell(M, S, n, prec, jsonpath=None):
    """CLI driver: certified cell + FLOAT route + (if installed) the moments
    cross-check; per-entry containment/deviation report; optional JSON."""
    rec = certified_cell(M, S, n, prec)
    hull = rec.pop("hull"); balls = rec.pop("entry_balls")
    fl = float_cell(M, S, n)
    try:
        mo = moments_cell(S, n, n_big=100)
    except ImportError:
        mo = None
        print("  SKIP moments cross-check: moments not installed (optional)")
    rec["float_scipy"] = fl
    rec["field_moments"] = mo
    # containment + deviation report
    comp = []
    for idx, i in enumerate(range(1, n)):
        ball = balls[idx]
        mid = float(ball.mid())
        row = {"i": i}
        routes = [("scipy", fl["sfs"][idx])]
        if mo is not None:
            routes += [("moments_n20", mo["direct"][idx]),
                       ("moments_n20_dtfine", mo["direct_dtfine"][idx]),
                       ("moments_n100proj", mo["projected"][idx])]
        for lbl, val in routes:
            row[lbl] = {"value": val,
                        "inside": bool(ball.contains(arb(val))),
                        "rel_dev_from_mid": abs(val - mid) / abs(mid) if mid else None}
        comp.append(row)
    rec["comparison"] = comp
    rec["date_utc"] = datetime.now(timezone.utc).isoformat()
    if jsonpath:
        if os.path.dirname(jsonpath):
            os.makedirs(os.path.dirname(jsonpath), exist_ok=True)
        with open(jsonpath, "w") as f:
            json.dump(rec, f, indent=1)
        print(f"  wrote {jsonpath}")
    return rec


if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="certified enclosures of the transient selected SFS")
    ap.add_argument("mode", choices=["cell", "all", "xcheck"])
    ap.add_argument("--M", type=int, default=200)
    ap.add_argument("--S", type=int, default=-5)
    ap.add_argument("--n", type=int, default=20)
    ap.add_argument("--prec", type=int, default=256)
    ap.add_argument("--json", default=None, help="cell: write the report here")
    ap.add_argument("--outdir", default=None,
                    help="all (required) / xcheck (optional): output directory")
    a = ap.parse_args()
    if a.mode == "cell":
        run_cell(a.M, a.S, a.n, a.prec, a.json)
    elif a.mode == "xcheck":
        r = xcheck()
        if a.outdir:
            out = os.path.join(a.outdir, f"xcheck_M{r['M']}_S{r['S']}.json")
            os.makedirs(a.outdir, exist_ok=True)
            json.dump(r, open(out, "w"), indent=1)
            print(f"  wrote {out}")
    else:
        if not a.outdir:
            ap.error("mode 'all' needs --outdir")
        for M, S in REFERENCE_CELLS:
            run_cell(M, S, a.n, a.prec,
                     os.path.join(a.outdir, f"cell_M{M}_S{S}_n{a.n}.json"))
