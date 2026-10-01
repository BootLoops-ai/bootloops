#!/usr/bin/env python3
"""Transient selected site-frequency spectrum for large samples (float PDE route).

WHAT IT COMPUTES. The expected unfolded site-frequency spectrum E_i,
i = 1..n-1 (theta = 1 units), of a sample of size n (n ~ 1000-2500 is the
design range) under genic (additive, h = 1/2) selection of scaled strength S
and a piecewise-constant population-size history epochs = [(nu, T), ...]
(past -> present; nu = N/N_ref, T in units of 2*N_ref generations) that
follows an ancestral Wright equilibrium at the same S. Route: direct
numerical integration of the Kimura forward diffusion (Kimura 1964) in the
Poisson-random-field setting (Sawyer & Hartl 1992) with time-dependent
population size (the non-equilibrium frequency-spectrum problem of Evans,
Shvets & Slatkin 2007), NOT a moment closure: no jackknife, no n_large; the
sample size enters only through the final binomial projection.

THE SCHEME. Working variable u(x,t) = x(1-x) f(x,t), f the density of the
derived-allele frequency. With gamma = S/2, time in 2*N_ref generations and
rho = nu the relative size of the current epoch, the forward equation

    df/dt = D d^2/dx^2 [x(1-x) f] - gamma d/dx [x(1-x) f],   D = 1/(2 rho)

becomes, in u, a CONSTANT-COEFFICIENT advection-diffusion flux law

    du/dt = -V(x) dJ/dx,   J = gamma u - D du/dx,   V = x(1-x).

Boundary conditions (the model's own, not numerical devices):
  u -> rho * theta  at x -> 0   (mutational-baseline influx; this replaces
                                 delta-function injection of new mutations),
  u -> 0            at x -> 1   ((1-x) f -> 0).

Space: Scharfetter-Gummel exponentially fitted interface fluxes
(Scharfetter & Gummel 1969) on a fixed log/lin/log composite grid (default
G = 3273 nodes, x down to 1e-8 at both ends, seams at the spacing-continuity
point). Time: theta-method, Crank-Nicolson by default with a Rannacher
backward-Euler startup (Rannacher 1984) to damp the initial-layer ringing;
banded (tridiagonal) implicit solves through LAPACK dgttrf/dgttrs with one
factorization per epoch. Because J has constant coefficients, the DISCRETE
stationary state of the scheme satisfies the exact two-point exponential
relation: the scheme's equilibrium is the exact Wright equilibrium
u = rho*theta*(1-e^{-S rho (1-x)})/(1-e^{-S rho}) sampled on the grid, up to
the x_min/x_max truncation, so the stationary error is dominated by the
projection quadrature, not by the PDE discretization.

Projection to sample size n: E_i = INT Bin(i; n, x) u(x)/(x(1-x)) dx by a
precomputed (n-1, G) kernel matrix (trapezoid weights; the binomial pmf is
<= 1, so no overflow). Utilities: hypergeometric down-sampling (formula-exact, float arithmetic) of an
unfolded expected SFS from n to m < n (project_down), and unfolded -> folded
minor-allele classes (fold_sfs).

CONVENTIONS. S = 4*N_ref*s with s the selection coefficient of the
heterozygote under additive selection; S < 0 deleterious, S > 0
advantageous (the same S as popcorn.sfs; dadi's gamma = S/2). theta = 1
ancestral (4*N_ref*mu = 1); spectra scale linearly in theta. Epoch list
[(nu, T), ...] runs past -> present after an ancestral equilibrium epoch at
nu = 1; an epoch with T <= 0 is an instantaneous size change (only the x -> 0
boundary value changes). Output E_i, i = 1..n-1, unfolded, theta = 1 units.
expected_sfs returns (E, meta); meta["negative_entries"] counts entries
below zero (they can appear in the exponentially suppressed high-frequency
tail at very strong negative selection) and callers must check it: the
engine does not raise.

REGISTER: FLOAT-VALIDATED. Double precision throughout, with a MEASURED
trust radius; nothing here is certified or exact. Measured against the
certified stationary references shipped under reference/transient/
(n = 1000, S in {0, -1, -5, -20, -100}; every entry i = 1..999 computed by
two independent arbitrary-precision routes of the package's certified
stationary evaluator, agreeing to >= 47 digits at generation; shipped to 20 significant digits): the projected analytic
equilibrium AND the equilibrium held for T = 0.5 through the integrator both
agree with the reference to <= 9.4e-5 max-relative in every frequency band
(i = 1, 2-5, 6-20, 21-100, 101-500, 501-900, 901-999) at every S in that set,
including S = -100 where the entries span 44 orders of magnitude; the error
is quadrature-dominated and concentrated in the log-spaced bands (mid
classes i = 101-900: <= 1e-9 projected, ~1e-6 held). For transient
histories no certified reference exists; the shipped self-convergence table
(two-epoch bottleneck-then-growth history [(0.33, 0.47), (3.38, 0.017)],
n = 1000, entries with E_i > 1e-5) gives time-step refinement
dt0 4e-4 -> 1e-4: <= 3.9e-4 for |S| <= 20, 4.0e-3 at S = -100, 1.1e-2 at
S = -200; grid refinement G 3273 -> 4611: <= 4.6e-5 for |S| <= 20, 1.2e-3
at S = -100, 5.2e-3 at S = -200. The exponentially suppressed tail entries
(E_i < 1e-5, zero data weight) show up to 9e-2 unfiltered dt sensitivity at
S = -200. Any other regime (other n, custom grids that break spacing
continuity, larger dt0, |S| > 200, sharper histories) is UNMEASURED. Where
the question is stationary, use the certified engine (popcorn.sfs /
certsfs.py) instead of this module.

THREADING. scipy's LAPACK honors OMP_NUM_THREADS / OPENBLAS_NUM_THREADS /
MKL_NUM_THREADS; the reference numbers above were produced single-threaded.

REFERENCES.
  Kimura M. (1964) Diffusion models in population genetics. J. Appl.
    Probab. 1:177-232. doi:10.2307/3211856
  Sawyer S.A., Hartl D.L. (1992) Population genetics of polymorphism and
    divergence. Genetics 132:1161-1176. doi:10.1093/genetics/132.4.1161
  Evans S.N., Shvets Y., Slatkin M. (2007) Non-equilibrium theory of the
    allele frequency spectrum. Theor. Popul. Biol. 71:109-119.
    doi:10.1016/j.tpb.2006.06.005
  Gutenkunst R.N., Hernandez R.D., Williamson S.H., Bustamante C.D. (2009)
    Inferring the joint demographic history of multiple populations from
    multidimensional SNP frequency data. PLoS Genet. 5:e1000695.
    doi:10.1371/journal.pgen.1000695  (finite-difference diffusion route)
  Jouganous J., Long W., Ragsdale A.P., Gravel S. (2017) Inferring the joint
    demographic history of multiple populations: beyond the diffusion
    approximation. Genetics 206:1549-1567. doi:10.1534/genetics.117.200493
    (the moment-closure route this engine deliberately avoids at large n)
  Scharfetter D.L., Gummel H.K. (1969) Large-signal analysis of a silicon
    Read diode oscillator. IEEE Trans. Electron Devices 16:64-77.
    doi:10.1109/T-ED.1969.16566
  Rannacher R. (1984) Finite element solution of diffusion problems with
    irregular data. Numer. Math. 43:309-327. doi:10.1007/BF01390130
"""
import numpy as np
from scipy.linalg import lapack
from scipy.special import gammaln

ENGINE_LABEL = "transient_sfs_sg_cn"


# ----------------------------------------------------------------- grid
def make_grid(x_min=1e-8, x_lin_lo=None, x_lin_hi=None, one_minus_max=1e-8,
              per_decade=160, dx_mid=8e-4):
    """Composite grid: log in x on [x_min, x_lin_lo], uniform on
    [x_lin_lo, x_lin_hi], log in (1-x) on [x_lin_hi, 1-one_minus_max].
    Seams default to the spacing-continuity point (log-spacing == dx_mid),
    which restores second-order trapezoid cancellation there (an abrupt
    spacing jump at a seam costs about 2.5e-4 relative on the classes
    that straddle it)."""
    if x_lin_lo is None:
        x_lin_lo = dx_mid * per_decade / np.log(10.0)
    if x_lin_hi is None:
        x_lin_hi = 1.0 - dx_mid * per_decade / np.log(10.0)
    nA = int(np.ceil(np.log10(x_lin_lo / x_min) * per_decade))
    gA = np.geomspace(x_min, x_lin_lo, nA + 1)
    nB = int(np.ceil((x_lin_hi - x_lin_lo) / dx_mid))
    gB = np.linspace(x_lin_lo, x_lin_hi, nB + 1)[1:]
    nC = int(np.ceil(np.log10((1 - x_lin_hi) / one_minus_max) * per_decade))
    gC = 1.0 - np.geomspace(1 - x_lin_hi, one_minus_max, nC + 1)[1:]
    return np.concatenate([gA, gB, gC])


def _bernoulli(z):
    """B(z) = z / (e^z - 1), stable; B(0) = 1."""
    z = np.asarray(z, float)
    out = np.empty_like(z)
    small = np.abs(z) < 1e-8
    out[small] = 1.0 - 0.5 * z[small]
    zb = np.clip(z[~small], -700.0, 700.0)
    out[~small] = zb / np.expm1(zb)
    return out


def ln_g_ratio(x, S):
    """ln g, g = (1-e^{-S(1-x)})/(1-e^{-S}); the Wright equilibrium is
    u = theta*g. S<0 deleterious; S=0 limit g = 1-x."""
    x = np.asarray(x, float)
    if S == 0.0:
        return np.log1p(-x)
    a = abs(S)

    def _ln_expm1(t):
        t = np.asarray(t, float)
        return np.where(t > 30.0, t, np.log(np.maximum(np.expm1(
            np.minimum(t, 30.0)), 1e-300)))
    return _ln_expm1(a * (1.0 - x)) - _ln_expm1(np.array(a))


class TransientSFSEngine:
    """Grid, projection kernel and stepping machinery for one sample size n.
    Build once per n (the (n-1, G) kernel costs ~0.3 s / ~80 MB at n=1000,
    default grid) and reuse across S and histories."""

    def __init__(self, n=1000, theta=1.0, grid_kw=None):
        self.n = int(n)
        self.theta = float(theta)
        self.x = make_grid(**(grid_kw or {}))
        x = self.x
        self.G = len(x)
        self.V = x * (1.0 - x)
        self.h = np.diff(x)                          # (G-1,) interface gaps
        self.hbar = np.empty(self.G)
        self.hbar[1:-1] = 0.5 * (x[2:] - x[:-2])
        self.hbar[0] = 0.5 * (x[1] - x[0])
        self.hbar[-1] = 0.5 * (x[-1] - x[-2])
        # projection kernel: E_i = K @ u  (trapezoid weights / (x(1-x)))
        i_arr = np.arange(1, self.n)
        lnC = (gammaln(self.n + 1) - gammaln(i_arr + 1)
               - gammaln(self.n - i_arr + 1))
        lnpmf = (lnC[:, None] + i_arr[:, None] * np.log(x)[None, :]
                 + (self.n - i_arr)[:, None] * np.log1p(-x)[None, :])
        w = self.hbar / self.V                       # du-quadrature weight
        self.K = np.exp(np.clip(lnpmf, -745.0, 30.0)) * w[None, :]

    # -------------------------------------------------------- equilibrium
    def eq_u(self, S, rho=1.0):
        """Wright equilibrium u = rho*theta*g(x; S*rho) at size rho."""
        return rho * self.theta * np.exp(ln_g_ratio(self.x, S * rho))

    # ------------------------------------------------------ epoch stepping
    def _epoch_operator(self, S, nu):
        """Tridiagonal L (interior nodes 1..G-2) and boundary source c for
        du/dt = L u + c with Dirichlet u_0 = nu*theta, u_{G-1} = 0."""
        Dc = 1.0 / (2.0 * nu)
        gam = S / 2.0
        P = gam * self.h / Dc                        # interface Peclet
        a = (Dc / self.h) * _bernoulli(-P)           # coeff of left node
        b = (Dc / self.h) * _bernoulli(P)            # coeff of right node
        # J_{j+1/2} = a_j u_j - b_j u_{j+1}
        m = self.G - 2                               # interior count
        fac = self.V[1:-1] / self.hbar[1:-1]
        # du_j/dt = -V_j (J_{j+1/2} - J_{j-1/2}) / hbar_j, j = 1..G-2
        #         = fac * ( a_{j-1} u_{j-1} - (b_{j-1} + a_j) u_j
        #                   + b_j u_{j+1} )
        sub = fac[1:] * a[1:-1]                      # couples u_{j-1}, j>=2
        dia = -fac * (b[:-1] + a[1:])
        sup = fac[:-1] * b[1:-1]                     # couples u_{j+1}
        c = np.zeros(m)
        c[0] = fac[0] * a[0] * (nu * self.theta)     # left Dirichlet source
        # right Dirichlet u_{G-1} = 0 contributes nothing
        return sub, dia, sup, c

    def _theta_step_factor(self, sub, dia, sup, dt, th):
        """Factor (I - th*dt*L) once; returns solver closure."""
        m = len(dia)
        dl = -th * dt * sub                          # (m-1,)
        d = 1.0 - th * dt * dia
        du = -th * dt * sup
        dlf, df, duf, du2, ipiv, info = lapack.dgttrf(dl, d, du)
        assert info == 0, f"dgttrf info={info}"

        def solve(rhs):
            xs, info2 = lapack.dgttrs(dlf, df, duf, du2, ipiv, rhs)
            assert info2 == 0
            return xs
        return solve

    def run_epochs(self, S, epochs, dt0=4e-4, min_steps=24, rannacher=4):
        """Integrate from ancestral equilibrium through epochs; returns
        interior+boundary u on the full grid."""
        u = self.eq_u(S, 1.0).copy()
        for nu, T in epochs:
            nu, T = float(nu), float(T)
            if T <= 0:
                # instantaneous size change: only BC value changes; state
                # carries over (density continuous)
                u[0] = nu * self.theta
                continue
            sub, dia, sup, c = self._epoch_operator(S, nu)
            nst = max(min_steps, int(np.ceil(T / dt0)))
            dt = T / nst
            ui = u[1:-1].copy()

            def _Lu(v):
                r = dia * v + c
                r[1:] += sub * v[:-1]
                r[:-1] += sup * v[1:]
                return r
            # Rannacher: damp CN ringing with a few backward-Euler halves
            nran = min(rannacher, nst)
            if nran:
                sBE = self._theta_step_factor(sub, dia, sup, 0.5 * dt, 1.0)
                for _ in range(2 * nran):
                    ui = sBE(ui + 0.5 * dt * c)
            sCN = self._theta_step_factor(sub, dia, sup, dt, 0.5)
            for _ in range(nst - nran):
                ui = sCN(ui + 0.5 * dt * _Lu(ui) + 0.5 * dt * c)
            u = np.empty(self.G)
            u[0] = nu * self.theta
            u[-1] = 0.0
            u[1:-1] = ui
        return u

    # ---------------------------------------------------------- projection
    def project(self, u):
        """u on grid -> E_i, i = 1..n-1 (unfolded, theta=1 units)."""
        return self.K @ u

    def expected_sfs(self, S, epochs, dt0=4e-4):
        """Transient expected SFS after the epoch list; returns (E, meta).
        Callers must check meta["negative_entries"] == 0."""
        u = self.run_epochs(S, epochs, dt0=dt0)
        E = self.project(u)
        return E, {"engine": ENGINE_LABEL, "G": self.G, "dt0": dt0,
                   "negative_entries": int((E < 0).sum())}

    def equilibrium_sfs(self, S):
        """Analytic Wright equilibrium at nu = 1 projected to the sample
        (float; for certified stationary values use popcorn.sfs)."""
        return self.project(self.eq_u(S, 1.0))


# ---------------------------------------------------- hypergeometric down-proj
def project_down(E, n, m):
    """Hypergeometric projection (formula-exact, float arithmetic) of an unfolded
    expected SFS from sample n to m < n: E_m(k) = sum_i E_n(i) * P[hyp(k; n, i, m)],
    segregating classes only (k = 1..m-1)."""
    E = np.asarray(E, float)
    i = np.arange(1, n)
    out = np.empty(m - 1)
    lnCnm = gammaln(n + 1) - gammaln(m + 1) - gammaln(n - m + 1)
    for k in range(1, m):
        ln = (gammaln(i + 1) - gammaln(k + 1) - gammaln(i - k + 1)
              + gammaln(n - i + 1) - gammaln(m - k + 1)
              - gammaln(n - i - (m - k) + 1) - lnCnm)
        ok = (i >= k) & (n - i >= m - k)
        t = np.where(ok, np.exp(np.where(ok, ln, -np.inf)), 0.0)
        out[k - 1] = float(t @ E)
    return out


def fold_sfs(E):
    """Unfolded E_i (i=1..n-1) -> folded minor classes k=1..n//2."""
    E = np.asarray(E, float)
    n = len(E) + 1
    K = n // 2
    F = np.empty(K)
    for k in range(1, K + 1):
        F[k - 1] = E[k - 1] + (E[n - k - 1] if k < n - k else 0.0)
    return F


def _selftest():
    eng = TransientSFSEngine(n=1000)
    i = np.arange(1, 1000)
    # 1. neutral equilibrium projection: 1/i
    E = eng.equilibrium_sfs(0.0)
    r1 = np.max(np.abs(E * i - 1.0))
    # 2. neutral nu=2 long hold -> 2/i (dynamics + BC test)
    E2, _ = eng.expected_sfs(0.0, [(2.0, 30.0)], dt0=2e-3)
    r2 = np.max(np.abs(E2 * i / 2.0 - 1.0))
    # 3. S=-5 equilibrium held T=0.5 at nu=1 must stay put (stationarity)
    E3a = eng.equilibrium_sfs(-5.0)
    E3b, _ = eng.expected_sfs(-5.0, [(1.0, 0.5)])
    r3 = np.max(np.abs(E3b / E3a - 1.0))
    print(f"selftest: eq-neutral maxrel {r1:.3e}; nu2-hold maxrel {r2:.3e}; "
          f"S-5 hold-still maxrel {r3:.3e}; G={eng.G}")
    assert r1 < 2e-4 and r2 < 5e-3 and r3 < 5e-5
    return r1, r2, r3


if __name__ == "__main__":
    _selftest()
