#!/usr/bin/env python3
r"""
kronecker.py -- Kronecker-Eisenstein kernel evaluator (high-precision, mpmath).

Supplies the marked-point
Kronecker-Eisenstein coefficients g^(n)(z, tau) and the AGW differential-form
COEFFICIENTS omega_k(z_j, N tau) that a fixed-base {dlog, dtau}-only eMPL
evaluator cannot represent.  These are exactly the
mass-dependent marked-point structures (z_i(t) moving on the torus) that the
(1,1,2) unequal-mass sunrise needs and that the equal-mass z_i = 1/3
degeneration switches off.

Primary reference: arXiv:1907.01251 (Bogner, Mueller-Stach, Weinzierl,
"The unequal mass sunrise integral expressed through iterated integrals on
Mbar_{1,3}", Nucl.Phys.B 954 (2020) 114991).  All equation numbers below refer
to that paper.

==========================================================================
DEFINITIONS (1907.01251)
==========================================================================
Kronecker function (eq.42), q = exp(i pi tau):

    F(x, y, tau) = pi * theta1'(0,q) * theta1(pi(x+y), q)
                                     / [ theta1(pi x, q) * theta1(pi y, q) ]

Generating function (eq.43):

    F(z, alpha, tau) = (1/alpha) sum_{n>=0} g^(n)(z, tau) alpha^n .

So the g^(n) are the Taylor coefficients of  alpha |-> alpha * F(z, alpha, tau).

q-bar series (eq.44-46), q_bar = exp(2 pi i tau), w_bar = exp(2 pi i z),
ELi_{n;m}(x;y;q) = sum_{j,k>=1} x^j y^k q^{jk}/(j^n k^m),
E_{n;m} = ELi_{n;m}(x;y) - (-1)^{n+m} ELi_{n;m}(x^-1;y^-1):

    g^(0) = 1
    g^(1) = -2 pi i [ (1+w_bar)/(2(1-w_bar)) + E_{0,0}(w_bar;1;q_bar) ]
    g^(n) = -(2 pi i)^n/(n-1)! [ -B_n/n + E_{0,1-n}(w_bar;1;q_bar) ]   (n>1)

Eisenstein representation (eq.50-54,57):
    g^(1) = E1(z,tau)                                      ( = pi theta1'/theta1 )
    g^(2) = -1/2 [ E2 - e2 - E1^2 ]                        ( E2 - e2 = Weierstrass-P )
    g^(3) =  1/6 [ 2 E3 - 3 (E2 - e2) E1 + E1^3 ]
with  E_k(z,tau) = sum_e 1/(z + n1 + n2 tau)^k  (Eisenstein summation),
      E2 = -E1' , E3 = +E1''/2 (since E_k' = -k E_{k+1}),
      e2(tau) = (pi^2/3) * Ehat2(tau),  Ehat2 = 1 - 24 sum_{n>=1} sigma1(n) q_bar^n.

Quasi-periodicity (eq.48):
    g^(n)(z+1,   tau) = g^(n)(z, tau)
    g^(n)(z+tau, tau) = sum_{j=0}^n ((-2 pi i)^j / j!) g^(n-j)(z, tau)

Poles (eq. text after 48): g^(n) has only SIMPLE poles in z.  g^(1) has a
simple pole of unit residue at every lattice point.  For n>1 the pole sits
only at lattice points NOT on the real axis.

Differential-form coefficients (eq.60, g^(-1):=0):
    omega_k(z,tau) = (2 pi)^{2-k} [ g^(k-1)(z,tau) dz
                                    + (k-1) g^(k)(z,tau) dtau/(2 pi i) ]
so the dz- and dtau-COEFFICIENTS of omega_k are returned separately.  With the
N-rescaling tau -> N tau the dtau coefficient carries an extra factor N
(eq.61: omega_2(z,2tau) has  2 g^(2)(z,2tau) dtau/(2 pi i)).

Pure-tau forms (eq.63,65):
    eta_2(tau) = b2(tau) dtau/(2 pi i),  b2 = e2(tau) - 2 e2(2 tau)  in M2(Gamma0(2))
    eta_4(tau) = (2 pi)^{-4} e4(tau) dtau/(2 pi i)   in M4(SL2(Z))

==========================================================================
API
==========================================================================
    nome(tau)                         q   = exp(i pi tau)
    qbar(tau)                         q_bar = exp(2 pi i tau)
    theta1(u, tau)                    Jacobi theta_1 (argument u, nome exp(i pi tau))
    theta1_prime0(tau)                d/du theta_1(u,tau)|_{u=0}
    kronecker_F(z, alpha, tau)        F(z, alpha, tau)  (eq.42)
    g_coeffs(z, tau, nmax=4, route="qbar")
                                      list [g^(0),...,g^(nmax)] ; route in
                                      {"qbar","theta","eisenstein"}
    g_n(n, z, tau, route="qbar")      single g^(n)
    g_diff(z_i, z_j, tau, n, ...)     g^(n)(z_i - z_j, tau)  (marked-point pair)
    omega_k(k, z, tau, N=1, route)    -> dict {"dz":..., "dtau":...} coefficients
    e2(tau), e4(tau)                  Eisenstein series e_k(tau) (eq.52)
    eta2_coeff(tau), eta4_coeff(tau)  dtau-coefficient of eta_2, eta_4 (eq.63,65)

All functions accept mpf/mpc and respect the ambient mp.mp.dps.  The qbar route
is the fast/default path; theta and eisenstein are independent cross-checks.
"""
import mpmath as mp

__all__ = [
    "nome", "qbar", "theta1", "theta1_prime0",
    "kronecker_F", "g_coeffs", "g_n", "g_diff",
    "omega_k", "e2", "e4", "eta2_coeff", "eta4_coeff",
]


def _twopii():
    return 2 * mp.pi * mp.mpc(0, 1)


def _sigma(n, k):
    """Divisor sigma function sigma_k(n) = sum_{d|n} d^k  (exact integer)."""
    s = 0
    d = 1
    while d * d <= n:
        if n % d == 0:
            s += d ** k
            e = n // d
            if e != d:
                s += e ** k
        d += 1
    return s


# -------------------------------------------------------------------------
# Theta layer
# -------------------------------------------------------------------------
def nome(tau):
    """q = exp(i pi tau)."""
    return mp.exp(mp.mpc(0, 1) * mp.pi * mp.mpc(tau))


def qbar(tau):
    """q_bar = exp(2 pi i tau) = q^2."""
    return mp.exp(_twopii() * mp.mpc(tau))


def theta1(u, tau):
    """Jacobi theta_1 with argument u and nome q = exp(i pi tau).

    Uses mpmath's converged jtheta q-series (Im tau > 0 required for |q|<1)."""
    return mp.jtheta(1, u, nome(tau))


def theta1_prime0(tau):
    """theta_1'(0, q) = d/du theta_1(u, q)|_{u=0}."""
    return mp.jtheta(1, 0, nome(tau), 1)


def kronecker_F(z, alpha, tau):
    """Kronecker function F(z, alpha, tau) of eq.42."""
    q = nome(tau)
    return (mp.pi * mp.jtheta(1, 0, q, 1) * mp.jtheta(1, mp.pi * (z + alpha), q)
            / (mp.jtheta(1, mp.pi * z, q) * mp.jtheta(1, mp.pi * alpha, q)))


# -------------------------------------------------------------------------
# g^(n) -- three independent routes
# -------------------------------------------------------------------------
def _g_theta(z, tau, nmax, r=None, N=None):
    """g^(n) from the Kronecker function F (eq.42,43) by a Cauchy contour integral
    in alpha (numerically stable, valid for any z, unlike a naive Taylor series):

        phi(a) = a * F(z,a,tau) = sum_{n>=0} g^(n)(z,tau) a^n
        g^(n)  = (1/2 pi i) oint phi(a)/a^{n+1} da ,   a = r e^{i theta}.

    The radius r is chosen smaller than the nearest alpha-singularity of F
    (alpha in the lattice Z + tau Z; nearest nonzero point has |.| >= min(1,Im tau)),
    so r = min(1, Im tau)/3 is safe."""
    z = mp.mpc(z)
    tau = mp.mpc(tau)
    if r is None:
        r = min(mp.mpf(1), abs(mp.im(tau))) / 3
    if N is None:
        # trapezoid on a circle is spectrally accurate; scale with precision
        N = max(8 * (nmax + 1), int(2.2 * mp.mp.dps) + 40)
    phi = lambda a: a * kronecker_F(z, a, tau)
    out = []
    twopi = 2 * mp.pi
    for n in range(nmax + 1):
        s = mp.mpc(0)
        for kk in range(N):
            th = twopi * kk / N
            a = r * mp.exp(mp.mpc(0, 1) * th)
            # 1/(2 pi i) * phi(a)/a^{n+1} * (i a dth) = phi(a)/(2 pi a^n) dth
            s += phi(a) / a ** n
        out.append(s / N)
    out[0] = mp.mpc(1)   # g^(0) = 1 exactly (eq.43); contour value confirms to ~dps
    return out


def _ELi(nn, mm, x, y, qb, J):
    """ELi_{nn;mm}(x;y;q_bar) = sum_{j,k>=1} x^j y^k q^{jk}/(j^nn k^mm) (eq.44)."""
    s = mp.mpc(0)
    for j in range(1, J + 1):
        xj = x ** j
        qbj = qb ** j
        # inner geometric-ish sum over k; |qb^j|<1 so k-tail decays fast
        qjk = qbj
        for k in range(1, J + 1):
            s += xj * (y ** k) * qjk / (mp.mpf(j) ** nn * mp.mpf(k) ** mm)
            qjk *= qbj
    return s


def _Ecombo(nn, mm, x, y, qb, J):
    """E_{nn;mm} = ELi(x;y) - (-1)^{nn+mm} ELi(x^-1;y^-1) (eq.45)."""
    return (_ELi(nn, mm, x, y, qb, J)
            - (-1) ** (nn + mm) * _ELi(nn, mm, 1 / x, 1 / y, qb, J))


def _g_qbar_cell(z, tau, nmax, J):
    """q-bar series of eq.46 -- VALID ONLY for z in the convergence strip
    0 < Im(z) < Im(tau).  Use _g_qbar for the all-domain wrapper."""
    tpi = _twopii()
    qb = qbar(tau)
    wb = mp.exp(tpi * mp.mpc(z))
    out = [mp.mpc(1)]                                              # g^(0)
    out.append(-tpi * ((1 + wb) / (2 * (1 - wb))
                       + _Ecombo(0, 0, wb, mp.mpf(1), qb, J)))     # g^(1)
    for n in range(2, nmax + 1):
        out.append(-(tpi ** n) / mp.factorial(n - 1)
                   * (-mp.bernoulli(n) / n
                      + _Ecombo(0, 1 - n, wb, mp.mpf(1), qb, J)))  # g^(n)
    return out


def _reduce_z(z, tau):
    """Reduce z to a representative z0 with 0 <= Re(z0) < 1 and
    0 <= Im(z0) < Im(tau), returning (z0, n_tau) where z = z0 + n_tau*tau (mod 1).
    The q-bar series converges for z0 in this half-open strip (the real axis
    Im(z0)=0 IS included: the (1+wbar)/(2(1-wbar)) term of g^(1) carries the
    z=0 lattice pole correctly there).  Points are NOT pushed off the real axis,
    so the pole structure is preserved."""
    z = mp.mpc(z)
    tau = mp.mpc(tau)
    imt = mp.im(tau)
    # subtract tau-multiples to land Im in [0, Im tau)
    n_tau = int(mp.floor(mp.im(z) / imt))
    z0 = z - n_tau * tau
    # guard against rounding pushing Im just below 0 or up to/over Im(tau)
    if mp.im(z0) < 0:
        z0 += tau
        n_tau -= 1
    elif mp.im(z0) >= imt:
        z0 -= tau
        n_tau += 1
    # remove integer part of Re (periodic in 1)
    z0 = z0 - mp.floor(mp.re(z0))
    return z0, n_tau


def _shift_tau_once(g, tau, forward):
    """One step of the eq.48 quasi-periodicity, returning g^(k)(z +/- tau) from g^(k)(z).
    forward=True : g^(k)(z+tau) = sum_{j=0}^k (-2pi i)^j/j! g^(k-j)(z)   (eq.48).
    forward=False: g^(k)(z-tau) = sum_{j=0}^k (+2pi i)^j/j! g^(k-j)(z)
                   (the inverse unipotent map; verified numerically)."""
    tpi = _twopii()
    base = -tpi if forward else tpi
    nmax = len(g) - 1
    out = []
    for k in range(nmax + 1):
        s = mp.mpc(0)
        for j in range(0, k + 1):
            s += base ** j / mp.factorial(j) * g[k - j]
        out.append(s)
    return out


def _shift_tau(g, tau, n):
    """Apply the tau-shift (eq.48) |n| times; n>0 shifts z by +tau, n<0 by -tau."""
    if n == 0:
        return list(g)
    forward = n > 0
    g = list(g)
    for _ in range(abs(n)):
        g = _shift_tau_once(g, tau, forward)
    return g


def _qbar_J_for(z0, tau, extra=30):
    """Pick a q-bar truncation J that reaches ~mp.mp.dps digits for the cell point
    z0 (with 0 <= Im(z0) < Im(tau)).  The two ELi pieces of E_{n;m} decay as
    |wbar|^j and |wbar|^{-j}|qbar|^{jk} respectively, i.e. the controlling per-index
    decay rate is min(Im(z0), Im(tau)-Im(z0)) (in units where the per-step base is
    exp(-2pi * that)).  J ~ dps / (that decay) padded."""
    imt = mp.im(tau)
    d = min(mp.im(z0), imt - mp.im(z0))
    rate = 2 * mp.pi * max(d, mp.mpf("1e-6")) / mp.log(10)   # decimal digits per index
    J = int(mp.mp.dps / max(rate, mp.mpf("0.01"))) + extra
    return max(J, 40)


# If the q-bar series would need more than this many terms to reach mp.mp.dps
# (i.e. the reduced point sits within ~ this distance of a horizontal lattice line,
# where exp(-2pi*dist) ~ 1 and the double sum stalls), fall back to the theta route,
# which is exact there (explicit theta log-derivative / Cauchy contour).
_QBAR_J_CAP = 600


def _g_qbar_strip(z0, tau, nmax, J):
    """g^(n)(z0) for z0 in the strip 0<=Im(z0)<Im(tau).  Folds via parity
    g^(n)(z0)=(-1)^n g^(n)(-z0) so the q-bar double series is always evaluated in
    the lower half-strip Im<=Im(tau)/2 (where the ELi(x;y) term dominates), and
    sizes J from the controlling decay rate.  -z0 (Im<0) reduces by +1*tau.

    Near a horizontal lattice line (Im(z0)->0 or ->Im(tau)) the q-bar series needs
    an impractically large J; there it transparently defers to the theta route
    (which represents g^(n) exactly via the theta functions, no real-axis issue)."""
    imt = mp.im(tau)
    if J is None and _qbar_J_for(z0, tau) > _QBAR_J_CAP:
        return _g_theta(z0, tau, nmax)
    if mp.im(z0) > imt / 2:
        # z0' = -z0 + tau  has Im = Im(tau)-Im(z0) < Im(tau)/2, in the strip.
        zp = -z0 + tau
        zp = zp - mp.floor(mp.re(zp))
        Jp = J if J is not None else _qbar_J_for(zp, tau)
        # g^(n)(z0) = (-1)^n g^(n)(-z0) ; and g^(n)(-z0) = g^(n)(zp - tau)
        #           = _shift_tau(g(zp), -1).
        g_zp = _g_qbar_cell(zp, tau, nmax, Jp)
        g_mz0 = _shift_tau(g_zp, tau, -1)                  # g^(n)(-z0)
        return [(-1) ** n * g_mz0[n] for n in range(len(g_mz0))]
    J0 = J if J is not None else _qbar_J_for(z0, tau)
    return _g_qbar_cell(z0, tau, nmax, J0)


def _g_qbar(z, tau, nmax, J):
    """All-domain g^(n) via the q-bar series (eq.46) + fundamental-cell reduction
    and the eq.48 quasi-periodicity.  Periodic in 1; tau-shifts are undone
    analytically; valid for ANY z in the upper-half-plane lattice.  J=None lets
    the routine auto-size the truncation to reach ~mp.mp.dps digits."""
    z0, n_tau = _reduce_z(z, tau)
    g0 = _g_qbar_strip(z0, tau, nmax, J)
    return _shift_tau(g0, tau, n_tau)


def e2(tau):
    """Eisenstein series e_2(tau) = sum_e' 1/(n1+n2 tau)^2 (Eisenstein summation,
    eq.52).  Value: e2 = (pi^2/3) * Ehat2(tau), Ehat2 = 1 - 24 sum sigma1(n) qbar^n."""
    qb = qbar(tau)
    # truncation set by ambient precision
    nterms = int(mp.mp.dps / (-mp.log10(abs(qb)) + mp.mpf("1e-30"))) + 20
    nterms = max(nterms, 40)
    s = mp.mpf(0)
    qn = mp.mpc(1)
    for n in range(1, nterms + 1):
        qn *= qb
        s += _sigma(n, 1) * qn
    Ehat2 = 1 - 24 * s
    return mp.pi ** 2 / 3 * Ehat2


def e4(tau):
    """Eisenstein series e_4(tau) = sum_e' 1/(n1+n2 tau)^4 (eq.52).
    Value: e4 = (pi^4/45)(1 + 240 sum sigma3(n) qbar^n)."""
    qb = qbar(tau)
    nterms = int(mp.mp.dps / (-mp.log10(abs(qb)) + mp.mpf("1e-30"))) + 20
    nterms = max(nterms, 40)
    s = mp.mpf(0)
    qn = mp.mpc(1)
    for n in range(1, nterms + 1):
        qn *= qb
        s += _sigma(n, 3) * qn
    return mp.pi ** 4 / 45 * (1 + 240 * s)


def _g_eisenstein(z, tau, nmax):
    """g^(1),g^(2),g^(3) via the Weierstrass/Eisenstein representation (eq.54).
    Independent of the q-bar and theta routes (uses analytic z-derivatives of
    E1 = pi theta1'(pi z)/theta1(pi z) and the weight-2 Eisenstein e2)."""
    if nmax > 3:
        raise ValueError("Eisenstein route implemented for n<=3 (eq.54).")
    z = mp.mpc(z)
    q = nome(tau)

    def E1f(zz):
        return mp.pi * mp.jtheta(1, mp.pi * zz, q, 1) / mp.jtheta(1, mp.pi * zz, q)

    E1 = E1f(z)
    out = [mp.mpc(1), E1]
    if nmax >= 2:
        E2 = -mp.diff(E1f, z)               # E_k' = -k E_{k+1}  =>  E2 = -E1'
        e2v = e2(tau)
        out.append(-mp.mpf(1) / 2 * (E2 - e2v - E1 ** 2))
    if nmax >= 3:
        E3 = mp.diff(E1f, z, 2) / 2          # E2' = -2 E3  =>  E3 = E1''/2
        out.append(mp.mpf(1) / 6 * (2 * E3 - 3 * (E2 - e2v) * E1 + E1 ** 3))
    return out


def g_coeffs(z, tau, nmax=4, route="qbar", J=None):
    """List [g^(0), g^(1), ..., g^(nmax)] at (z, tau).

    route: "qbar" (eq.46, fast, default), "theta" (eq.42/43 alpha-expansion),
           "eisenstein" (eq.54, n<=3).
    J:     q-bar series truncation (default scales with mp.mp.dps and Im tau)."""
    if route == "qbar":
        # J=None -> auto-size per reduced representative inside _g_qbar (recommended).
        return _g_qbar(z, tau, nmax, None if J is None else int(J))
    if route == "theta":
        return _g_theta(z, tau, nmax)
    if route == "eisenstein":
        return _g_eisenstein(z, tau, nmax)
    raise ValueError(f"unknown route {route!r}")


def g_n(n, z, tau, route="qbar", J=None):
    """Single coefficient g^(n)(z, tau)."""
    return g_coeffs(z, tau, nmax=n, route=route, J=J)[n]


def g_at_zero(n, tau, J=None):
    """g^(n)(0, tau) for n>=2 (the FINITE second-kind values at a coincident pair,
    e.g. the m1=m2 punctures z1=z2 of the (1,1,2) sunrise -> J6=0 structure).

    From eq.46 at wbar=1:  g^(n)(0) = -(2pi i)^n/(n-1)! [ -B_n/n + E_{0,1-n}(1;1;qbar) ].
    With x=y=1, E_{0,1-n}(1;1) = ELi_{0,1-n}(1;1)[1-(-1)^{1-n}], and
    ELi_{0,1-n}(1;1;qbar) = sum_{j,k>=1} k^{n-1} qbar^{jk}.  For odd n>=3 this gives 0
    (consistent with parity g^(odd)(0)=0); g^(1)(0) is a POLE (raises)."""
    if n == 1:
        raise ValueError("g^(1)(0) is a simple pole (no finite value at z=0).")
    if n < 0:
        raise ValueError("n must be >= 0")
    if n == 0:
        return mp.mpc(1)
    tpi = _twopii()
    qb = qbar(tau)
    if n % 2 == 1:
        return mp.mpc(0)            # B_n=0 (odd) and E-combo vanishes
    # n even: E_{0,1-n}(1;1) = 2 * sum_{j,k>=1} k^{n-1} qbar^{jk}
    if J is None:
        lq = -mp.log10(abs(qb))
        J = max(int(mp.mp.dps / max(lq, mp.mpf("0.05"))) + 30, 40)
    s = mp.mpc(0)
    for j in range(1, J + 1):
        qbj = qb ** j
        for k in range(1, J + 1):
            s += mp.mpf(k) ** (n - 1) * qbj ** k
    E = 2 * s
    return -(tpi ** n) / mp.factorial(n - 1) * (-mp.bernoulli(n) / n + E)


def g_diff(z_i, z_j, tau, n, route="qbar", J=None, coincide_tol=None):
    """g^(n)(z_i - z_j, tau): the marked-point pair kernel argument z_i - z_j.

    If z_i and z_j coincide (mod the lattice) within coincide_tol, the argument is
    a lattice point: g^(1) is then a pole (raises ValueError) while g^(n>=2) returns
    the finite g_at_zero value (the second-kind objects that survive at z1=z2)."""
    zz = mp.mpc(z_i) - mp.mpc(z_j)
    if coincide_tol is None:
        coincide_tol = mp.mpf(10) ** (-int(mp.mp.dps * 0.6))
    z0, _ = _reduce_z(zz, tau)
    # distance to nearest lattice point (0 or 1 in Re, 0 in Im within the strip)
    re0 = min(mp.re(z0), abs(1 - mp.re(z0)))
    if abs(mp.im(z0)) < coincide_tol and re0 < coincide_tol:
        return g_at_zero(n, tau, J=J)
    return g_n(n, zz, tau, route=route, J=J)


# -------------------------------------------------------------------------
# Differential-form coefficients omega_k(z, N tau)  (eq.60-61)
# -------------------------------------------------------------------------
def omega_k(k, z, tau, N=1, route="qbar", J=None):
    """Differential-form coefficient omega_k(z, N tau) (eq.60).

    Returns dict {"dz": c_dz, "dtau": c_dtau} so that
        omega_k = c_dz * dz + c_dtau * dtau,
    with
        c_dz   = (2 pi)^{2-k} g^(k-1)(z, N tau)
        c_dtau = (2 pi)^{2-k} (k-1) g^(k)(z, N tau) * N / (2 pi i)
    (the N factor from tau -> N tau on the dtau leg, eq.61; for k=1 the dtau
    coefficient vanishes since the (k-1) prefactor is 0; g^(-1):=0 so omega_0
    has c_dz=0 and c_dtau = 2 pi i * N, i.e. omega_0(.,N tau)=2 pi i N dtau)."""
    Ntau = N * mp.mpc(tau)
    pref = (2 * mp.pi) ** (2 - k)
    # g^(k-1) and g^(k); g^(-1):=0 handled below
    nmax = max(k, 1)
    gs = g_coeffs(z, Ntau, nmax=nmax, route=route, J=J)
    g_km1 = mp.mpc(0) if k - 1 < 0 else gs[k - 1]
    g_k = gs[k] if k <= nmax else mp.mpc(0)
    c_dz = pref * g_km1
    c_dtau = pref * (k - 1) * g_k * N / _twopii()
    return {"dz": c_dz, "dtau": c_dtau}


def eta2_coeff(tau):
    """dtau-coefficient of eta_2(tau) = b2(tau) dtau/(2 pi i), b2=e2(tau)-2 e2(2tau)
    (eq.63).  Returns b2(tau)/(2 pi i)."""
    b2 = e2(tau) - 2 * e2(2 * mp.mpc(tau))
    return b2 / _twopii()


def eta4_coeff(tau):
    """dtau-coefficient of eta_4(tau) = (2 pi)^{-4} e4(tau) dtau/(2 pi i) (eq.65)."""
    return (2 * mp.pi) ** (-4) * e4(tau) / _twopii()


if __name__ == "__main__":
    mp.mp.dps = 50
    tau = mp.mpc(0.3, 1.1)
    z = mp.mpc(0.22, 0.37)
    print("# Kronecker-Eisenstein g^(n)(z,tau), three routes")
    print(f"tau = {tau},  z = {z}")
    gq = g_coeffs(z, tau, 4, "qbar")
    gt = g_coeffs(z, tau, 4, "theta")
    ge = g_coeffs(z, tau, 3, "eisenstein")
    for n in range(5):
        print(f"g^{n} = {mp.nstr(gq[n], 20)}")
        print(f"     theta agree to     {mp.nstr(abs(gq[n]-gt[n]), 3)}")
        if n <= 3:
            print(f"     eisenstein agree to {mp.nstr(abs(gq[n]-ge[n]), 3)}")
