#!/usr/bin/env python3
r"""rowscripts_empl.py -- SHARED self-contained Kronecker-eMPL machinery for the
unequal-mass two-loop sunrise at eps^0, mass configurations (1,1,2), (1,2,3) and
(1,1,4) (labeled rows 9, 10, 11 in the data files row09/row10/row11_data.json).

Implements the uniform per-puncture closed form for the sunrise family:

    E^(0)(t) = J[1,1,1]^(0)(t) / psihat1(t)
             = - C_{4,2}(t) - (6/(2 pi)^3) * (1/3) * sum_{j=1}^{3} [ W(z_j,1) - 8 W(z_j,2) ]

with
  * psihat1(t) = |psi1_F(t)| / pi, the holomorphic period of the BMSW Feynman
    curve E_{3,F} (arXiv:1907.01251), psi1_F = (2/sqrt(Z_{3,F})) K(k_F^2);
  * z_1, z_2, z_3 the Abel-Jacobi images of the three Feynman punctures on the
    torus (incomplete elliptic integral F(asin u', k'_F^2) / (2 K(k'_F^2)), BMSW
    def_coordinate_torus); z_j real for Euclidean t < 0, z1+z2+z3 = 1, freezing
    to the Gamma_1(6) 3-torsion z* = 1/3 at equal mass;
  * W(z,N) = I(1, g^(3)(z, N tau_C); q_C) = sum_{n>=1} b_n(z,N) q_C^n / n^2, the
    depth-two Kronecker-Eisenstein eMPL word; the kernel q-series is
      g^(3)(z, N tau_C) = -(2 pi i)^3/2! * sum_{j,k>=1} (w^j - w^{-j}) k^2 q_C^{N j k},
    w = exp(2 pi i z)  (Kronecker-coefficient qbar-series, eq. (46) of
    arXiv:1907.01251, in the ELi notation of Adams-Weinzierl arXiv:1704.08895; qbar = exp(2 pi i * N tau_C) = q_C^N);
  * tau_C = (tau_F + 1)/2, q_C = -exp(i pi tau_F): NEGATIVE REAL for t < 0;
  * C_{4,2}(t) = sum_j (1/2i) [ Li2(w_j) - Li2(1/w_j) ], w_j = exp(2 pi i z_j),
    the elliptic-dilogarithm boundary term (equal-mass limit: (3/2) sqrt3 L(chi_{-3},2)).

The three PSLQ-exact rational coefficients [-1, -1/3, +8/3] were established
offline by integer-relation (PSLQ) closures at 149-160 digits against held-out
oracles; nothing here is fit at runtime -- the formula is evaluated with
ZERO free parameters.

DEGENERACY: for m1 = m2 (rows 9 and 11) the punctures satisfy z1 = z2 exactly.
marked_z_mults() groups coincident punctures and returns [(z, multiplicity)], so
the per-puncture sum handles the degenerate and generic cases uniformly (and
computes each distinct W only once).

DOMAIN: Euclidean t < 0 (below all thresholds; q_C negative real, all z_j real).
The q-series is auto-sized from the proven tail bound |b_n| <= C n^2 with the
dps-independent prefactor C = (2 pi)^3 zeta(2) = 408.03 < 10^3 budgeted
EXPLICITLY (see Nq_for) => |tail| <= C |q_C|^(Nq+1) / (1 - |q_C|); Nq_for
certifies that bound at runtime and GROWS Nq until it clears 10^-(dps+extra).
Points with |q_C| >= QMAX = 0.8995 (chosen so 1/(1-|q_C|) < 10 strictly) are
refused explicitly rather than silently truncated.  |q_C| -> 0 at the soft
point t -> 0^- (harmless; control passes at full precision at t = -1e-3) and
creeps toward 1 only at astronomically deep |t| (|q_C| = 0.89 at t = -1e12 for
(1,1,2); measured control at t = -1e4: 42/50d).  Gate-verified window:
-7 <= t <= -1/3 (|q_C| ~ 0.03 - 0.2); practical domain: any Euclidean
t in [-1e6, 0) at the few-digit-loss level documented above.
fcurve carries the Abel-Jacobi BRANCH REPAIR (z3 -> 1 - z1 - z2 when the
ellipf/asin principal branch reflects across the half-period); without it
(1,1,4) fails on -2 < t < 0.

PORTABILITY: pure mpmath, no file/network access, no import-time mpf constants
(everything precision-dependent is a zero-arg function evaluated at the ACTIVE
mp.mp.dps -- the import-time dps=15 footgun).

Positive control: equal_mass_control(t) evaluates the SAME machinery at
(m1,m2,m3) = (1,1,1) and compares against the independent classical Gamma_1(6)
route (Eisenstein q-series e1, e2 -> newform f3 -> I(1,f3;q_C), boundary
(3/2) sqrt3 L(chi_{-3},2) via Hurwitz zeta) -- two constructions that share no
code path beyond the curve itself.

Value-level certified-bound REPORTING:
qtail_bound() evaluates the same proven per-word tail bound Nq_for already
enforces; predict_E/predict_J propagate it to value level through the l1-norm
of the exact prefactors (diag keys err_E / err_J); run_row prints the certified
BOUND lines (gate table summary + --point).  Print-only -- the bound never
enters the value math; the runtime enforcement itself is the fail-closed
ValueError at the Nq cap.
"""
import mpmath as mp

# QMAX = 0.8995 (exact rational 1799/2000; not precision-dependent), NOT 0.9:
# the tail bound carries a 1/(1-|q_C|) factor whose supremum on |q_C| < QMAX is
# 1/(1-0.8995) = 9.9503 < 10 -- strictly ONE decimal digit, matching the single
# QLOG guard digit budgeted in Nq_for.  At the old 0.9 boundary the supremum is
# exactly 10, so the factor could shave the last guard digit.
QMAX = mp.mpf(1799) / 2000   # refuse |q_C| >= 0.8995


# --------------------------------------------------------------------------
# precision-correct lazy constants (never frozen at import-time dps)
# --------------------------------------------------------------------------
def two_pi_i():
    return 2 * mp.pi * mp.mpc(0, 1)


def a_norm():
    """A_NORM = 6/(2 pi)^3: f3 = A_NORM [g3(z*,tau_C) - 8 g3(z*,2 tau_C)] at z*=1/3."""
    return 6 / (2 * mp.pi) ** 3


def g3_pref():
    """-(2 pi i)^3 / 2!  (B_3 = 0, so no constant term)."""
    return -(two_pi_i()) ** 3 / 2


def _csqrt(x):
    return mp.sqrt(mp.mpc(x))


# --------------------------------------------------------------------------
# 1. BMSW Feynman curve E_{3,F}: roots, moduli, period, Abel-Jacobi punctures
# --------------------------------------------------------------------------
def fcurve(t, m1, m2, m3, mu=1):
    r"""All E_{3,F} data at t (arXiv:1907.01251 def_roots / points_on_E_hat /
    def_coordinate_torus).  Masses are the ACTUAL masses (not squared).
    Returns dict with e1F, e2F, e3F, kF2, kFp2, tauF, psi1F, z1F, z2F, z3F."""
    t = mp.mpc(t)
    m1, m2, m3, mu = mp.mpf(m1), mp.mpf(m2), mp.mpf(m3), mp.mpf(mu)
    m1s, m2s, m3s = m1 * m1, m2 * m2, m3 * m3
    M100 = m1s + m2s + m3s
    mu1 = -m1 + m2 + m3          # pseudo-threshold masses
    mu2 = m1 - m2 + m3
    mu3 = m1 + m2 - m3
    mu4 = m1 + m2 + m3
    Delta = mu1 * mu2 * mu3 * mu4
    mu4p = mu ** 4
    rad = 3 * (_csqrt(mu1 * mu1 - t) * _csqrt(mu2 * mu2 - t)
               * _csqrt(mu3 * mu3 - t) * _csqrt(mu4 * mu4 - t))
    e1F = (-t * t + 2 * M100 * t + Delta + rad) / (24 * mu4p)
    e2F = (-t * t + 2 * M100 * t + Delta - rad) / (24 * mu4p)
    e3F = (2 * t * t - 4 * M100 * t - 2 * Delta) / (24 * mu4p)
    Z1F = e3F - e2F
    Z2F = e1F - e3F
    Z3F = e1F - e2F
    kF2 = Z1F / Z3F
    kFp2 = -Z1F / Z2F            # puncture-map modulus k'_F^2
    tauF = mp.mpc(0, 1) * mp.ellipk(1 - kF2) / mp.ellipk(kF2)
    psi1F = 2 / _csqrt(Z3F) * mp.ellipk(kF2)
    Kp = mp.ellipk(kFp2)

    def xhat(mi2, mj2):          # xhat_{i,j,F} = e3F + m_i^2 m_j^2 / mu^4
        return e3F + mi2 * mj2 / mu4p

    def zF(xjk):                 # Abel-Jacobi image (incomplete elliptic integral)
        up = _csqrt((e1F - e3F) / (xjk - e3F))
        return mp.ellipf(mp.asin(up), kFp2) / (2 * Kp)

    z1F = zF(xhat(m2s, m3s))     # i=1, (j,k)=(2,3)
    z2F = zF(xhat(m3s, m1s))     # i=2, (j,k)=(3,1)
    z3F = zF(xhat(m1s, m2s))     # i=3, (j,k)=(1,2)
    # Abel-Jacobi BRANCH REPAIR (verified to 70 digits):
    # the lattice relation z1+z2+z3 = 1 must hold exactly (AGW: z3 = 1 - z1 - z2).
    # The ellipf/asin principal branch reflects z3 -> 1 - z3 once z3 crosses the
    # half-period (e.g. (1,1,4) for -2 < t < 0: raw sum = 0.903 at t = -1).
    # Enforcing the relation repairs the branch; validated against the held-out
    # AMFlow oracle at (1,1,4) t = -1 (stored 288-digit sample) to full working
    # precision.  z1, z2 stay far from 1/2 on the documented domain.
    ssum = z1F + z2F + z3F
    if abs(ssum - 1) > mp.mpf(10) ** (-mp.mp.dps // 2):
        z3F = 1 - z1F - z2F
    return dict(e1F=e1F, e2F=e2F, e3F=e3F, kF2=kF2, kFp2=kFp2, tauF=tauF,
                psi1F=psi1F, z1F=z1F, z2F=z2F, z3F=z3F)


def psihat1(t, masses, mu=1):
    """AW-normalized dressing period psihat1 = |psi1_F|/pi (real for t < 0)."""
    return abs(fcurve(t, *masses, mu=mu)["psi1F"]) / mp.pi


def nome_qC(tauF):
    """q_C = -exp(i pi tau_F) = exp(2 pi i tau_C), tau_C = (tau_F+1)/2.
    NEGATIVE REAL for Euclidean t < 0."""
    return -mp.e ** (mp.mpc(0, 1) * mp.pi * tauF)


def check_domain(t, qC):
    """Refuse points outside the documented domain instead of silently truncating."""
    if not (mp.im(mp.mpc(t)) == 0 and mp.re(mp.mpc(t)) < 0):
        raise ValueError(f"t = {t}: only Euclidean t < 0 is supported by this script")
    if abs(qC) >= QMAX:
        raise ValueError(f"t = {t}: |q_C| = {mp.nstr(abs(qC), 6)} >= {mp.nstr(QMAX, 4)} "
                         "(series-domain limit; t too close to 0^- or too deep)")


def marked_z_mults(cv, tol_digits=None):
    """Group the three punctures into [(z, multiplicity)] (handles the z1 = z2
    degeneracy of m1 = m2 rows exactly; generic rows give three mult-1 entries)."""
    if tol_digits is None:
        tol_digits = mp.mp.dps // 2
    tol = mp.mpf(10) ** (-tol_digits)
    zs = [cv["z1F"], cv["z2F"], cv["z3F"]]
    groups = []
    for z in zs:
        for g in groups:
            if abs(z - g[0]) < tol:
                g[1] += 1
                break
        else:
            groups.append([z, 1])
    return [(g[0], g[1]) for g in groups]


# --------------------------------------------------------------------------
# 2. Kronecker-Eisenstein depth-two words W(z,N) = I(1, g^(3)(z, N tau_C); q_C)
# --------------------------------------------------------------------------
def Nq_for(qC, extra=12):
    r"""Series length CERTIFIED at runtime from the PROVEN tail bound.

    Derivation of the dps-INDEPENDENT prefactor C (previously this constant was
    silently absorbed by the +extra guard digits; now it is budgeted explicitly):
    the coefficients returned by g3_qC_coeffs are
        b_n = -(2 pi i)^3/2! * S_n,  S_n = sum_{j,k>=1, N j k = n} (w^j - w^{-j}) k^2,
    and on the wired domain z is real, so |w| = |e^{2 pi i z}| = 1 and
    |w^j - w^{-j}| <= 2, giving
        |S_n| <= 2 sum_{k | n/N} k^2 = 2 sigma_2(n/N) <= 2 zeta(2) (n/N)^2 <= 2 zeta(2) n^2
    (sigma_2(m) = sum_{d|m} d^2 = m^2 sum_{d|m} d^-2 < zeta(2) m^2).  Hence
        |b_n| <= ((2 pi)^3 / 2) * 2 zeta(2) n^2 = C n^2,   C = (2 pi)^3 zeta(2) = 408.03,
    and with the 1/n^2 of the depth-two word the truncation tail obeys
        |tail(Nq)| <= C |q|^(Nq+1) / (1 - |q|).
    C < 10^3 costs CLOG = ceil(log10 C) = 3 digits; |q| < QMAX = 0.8995 keeps
    1/(1-|q|) <= 9.9503 < 10, i.e. QLOG = 1 digit.  Both enter the guess below.

    The closed-form solve is only the STARTING guess (fast path: in the wired
    domain it already over-satisfies the bound, so no escalation happens).  The
    gate then ACTS: it evaluates the bound and grows Nq geometrically until
        C |q|^(Nq+1) / (1 - |q|) < 10^-(dps+extra)
    holds (with a factor-1/2 slack absorbing the dps+10 evaluation rounding),
    raising with a precise message if no admissible depth exists."""
    aq = abs(qC)
    if aq <= 0:
        return 60
    if aq >= 1:
        raise ValueError(f"Nq_for: |q_C| = {mp.nstr(aq, 6)} >= 1 -- q-series divergent, "
                         "no truncation depth certifies the tail")
    dps = mp.mp.dps
    CLOG, QLOG = 3, 1   # ceil(log10 C) and ceil(log10 sup 1/(1-|q|)); see docstring
    need = (dps + extra + CLOG + QLOG) / max(-mp.log10(aq), mp.mpf("0.046"))
    Nq = max(int(need) + 8, 60)
    # ACTING gate: certify the runtime bound by direct evaluation, escalating Nq
    # until it clears the target (mp.workdps context: no global dps mutation).
    NQ_CAP = 500000
    with mp.workdps(dps + 10):
        C = (2 * mp.pi) ** 3 * mp.zeta(2)
        target = mp.mpf(10) ** (-(dps + extra)) / 2   # /2 = directed-rounding slack
        while C * aq ** (Nq + 1) / (1 - aq) >= target:
            if Nq >= NQ_CAP:
                # Name the ACHIEVED bound as a number (not only the target), so
                # the raise states the actual shortfall instead of an
                # unquantified inequality.
                achieved = C * aq ** (Nq + 1) / (1 - aq)
                raise ValueError(
                    f"Nq_for: proven tail bound C|q_C|^(Nq+1)/(1-|q_C|) = "
                    f"{mp.nstr(achieved, 4)} at Nq = {NQ_CAP} still >= target "
                    f"{mp.nstr(target, 4)} (= 10^-({dps}+{extra})/2) for "
                    f"|q_C| = {mp.nstr(aq, 6)} "
                    "-- point too deep for the q-series route at this precision")
            Nq = min(Nq + max(Nq // 4, 16), NQ_CAP)
    return Nq


def qtail_bound(qC, Nq):
    r"""The PROVEN per-word truncation bound actually achieved at depth Nq:
        |tail(Nq)| <= C |q_C|^(Nq+1) / (1 - |q_C|),   C = (2 pi)^3 zeta(2),
    i.e. the same certified expression Nq_for enforces, EVALUATED at the depth
    in use so it can be propagated to value level and printed (a BOUND, not an
    estimate; see Nq_for's docstring for the coefficient-bound derivation).
    Evaluated under mp.workdps(dps+10); returns an mpf."""
    aq = abs(qC)
    if aq <= 0:
        return mp.mpf(0)
    with mp.workdps(mp.mp.dps + 10):
        C = (2 * mp.pi) ** 3 * mp.zeta(2)
        return +(C * aq ** (Nq + 1) / (1 - aq))


def g3_qC_coeffs(zj, N, Nq):
    r"""Coefficients b_n (n = 0..Nq) of g^(3)(z_j, N tau_C) as a power series in q_C:
    g^(3)(z, N tau_C) = g3_pref() * sum_{j,k>=1} (w^j - w^{-j}) k^2 q_C^{N j k}."""
    wbar = mp.exp(two_pi_i() * zj)
    b = [mp.mpc(0)] * (Nq + 1)
    j = 1
    while N * j <= Nq:
        wj = wbar ** j - wbar ** (-j)
        k = 1
        while N * j * k <= Nq:
            b[N * j * k] += wj * (k * k)
            k += 1
        j += 1
    pref = g3_pref()
    return [pref * x for x in b]


def I1_g3(zj, N, qC, Nq):
    """Depth-two word W(z_j, N) = sum_{n>=1} b_n(z_j, N) q_C^n / n^2."""
    b = g3_qC_coeffs(zj, N, Nq)
    s = mp.mpc(0)
    qn = qC
    for n in range(1, Nq + 1):
        s += b[n] * qn / (n * n)
        qn *= qC
    return s


# --------------------------------------------------------------------------
# 3. Elliptic-dilogarithm boundary C_{4,2} and the assembled closed form
# --------------------------------------------------------------------------
def C42_from_z(z_mults):
    r"""C_{4,2} = sum_j (1/2i)[Li2(w_j) - Li2(1/w_j)], w_j = exp(2 pi i z_j)."""
    tot = mp.mpc(0)
    for z, mult in z_mults:
        w = mp.exp(two_pi_i() * z)
        tot += mult * (mp.polylog(2, w) - mp.polylog(2, 1 / w)) / (2 * mp.mpc(0, 1))
    return tot


def predict_E(t, masses, mu=1, diag=None):
    r"""The printed closed form (zero free parameters):
        E^(0)(t) = -C_{4,2}(t) - (6/(2 pi)^3) * (1/3) sum_j [ W(z_j,1) - 8 W(z_j,2) ].
    Returns a real mpf for Euclidean t < 0.  diag (optional dict) receives the
    frame data (tauF, qC, z-multiplicities, Nq, psihat1) plus the value-level
    certified truncation bound err_E: each word W(z_j,N) is truncated at the
    SAME depth Nq, so its tail obeys the proven per-word bound qtail_bound
    (enforced < 10^-(dps+12) by Nq_for); the per-word bounds propagate through
    the l1-norm of the EXACT prefactors, |a_norm()/3| * sum_j mult_j * (1 + 8).
    C_{4,2} and the curve data are mpmath builtins evaluated at working
    precision (no truncation on this code path), so err_E is a certified BOUND
    on the q-series truncation error of E (a bound, not an estimate)."""
    cv = fcurve(t, *masses, mu=mu)
    qC = nome_qC(cv["tauF"])
    check_domain(t, qC)
    zm = marked_z_mults(cv)
    Nq = Nq_for(qC)
    ell = mp.mpc(0)
    for z, mult in zm:
        ell += mult * (I1_g3(z, 1, qC, Nq) - 8 * I1_g3(z, 2, qC, Nq))
    ell = a_norm() * ell / 3
    c42 = C42_from_z(zm)
    if diag is not None:
        # value-level certified bound (print-only; never enters the value math)
        bW = qtail_bound(qC, Nq)                       # per-word proven tail bound
        l1 = mp.fsum(mult * (1 + 8) for _, mult in zm)  # exact |1| + |-8| per puncture
        err_E = abs(a_norm()) / 3 * l1 * bW
        diag.update(tauF=cv["tauF"], qC=qC, z_mults=zm, Nq=Nq,
                    psihat1=abs(cv["psi1F"]) / mp.pi, C42=c42, ell_block=ell,
                    qtail_word=bW, err_E=err_E)
    return (-c42 - ell).real


def predict_J(t, masses, mu=1, diag=None):
    """J[1,1,1]^(0)(t) = psihat1(t) * E^(0)(t): directly comparable to the raw
    AMFlow eps^0 master (indices (1,1,1,0,0), d = 2 - 2 eps normalization of the
    stored sample files).  diag gains err_J = psihat1 * err_E: the value-level
    certified q-series truncation bound of J (exact positive prefactor psihat1;
    BOUND, not estimate)."""
    d = {} if diag is None else diag
    E = predict_E(t, masses, mu=mu, diag=d)
    d["err_J"] = d["psihat1"] * d["err_E"]
    return d["psihat1"] * E


# --------------------------------------------------------------------------
# 4. Positive control: equal-mass (1,1,1) vs the independent Gamma_1(6) route
# --------------------------------------------------------------------------
def _chi_m3(n):
    r = n % 3
    return 1 if r == 1 else (-1 if r == 2 else 0)


def _f3_qcoeffs(N):
    """q-coefficients of f3 = 36 sqrt3 (e1^3 - e1^2 e2 - 4 e1 e2^2 + 4 e2^3),
    e1 = 1/6 + sum_{m>=1} (sum_{d|m} chi_{-3}(d)) q^m, e2(q) = e1(q^2)
    (conventions of arXiv:1704.08895; same route as the row-8 gate script)."""
    e1 = [mp.mpf(1) / 6] + [mp.mpf(0)] * N
    for d in range(1, N + 1):
        c = _chi_m3(d)
        if c:
            for m in range(d, N + 1, d):
                e1[m] += c
    e2 = [mp.mpf(0)] * (N + 1)
    for m in range(0, N // 2 + 1):
        e2[2 * m] = e1[m]
    e2[0] = mp.mpf(1) / 6

    def pmul(a, b):
        r = [mp.mpf(0)] * (N + 1)
        for i, ai in enumerate(a):
            if ai:
                for j in range(0, N + 1 - i):
                    if b[j]:
                        r[i + j] += ai * b[j]
        return r

    e1e1 = pmul(e1, e1)
    e1e1e1 = pmul(e1e1, e1)
    e1e1e2 = pmul(e1e1, e2)
    e2e2 = pmul(e2, e2)
    e1e2e2 = pmul(e1, e2e2)
    e2e2e2 = pmul(e2e2, e2)
    s3 = 36 * mp.sqrt(3)
    return [s3 * (e1e1e1[n] - e1e1e2[n] - 4 * e1e2e2[n] + 4 * e2e2e2[n])
            for n in range(N + 1)]


def L_chi_m3_2():
    """L(chi_{-3}, 2) = (zeta(2,1/3) - zeta(2,2/3)) / 9 (Hurwitz-zeta closed form)."""
    return (mp.zeta(2, mp.mpf(1) / 3) - mp.zeta(2, mp.mpf(2) / 3)) / 9


def equal_mass_control(t):
    """At (m1,m2,m3) = (1,1,1) compare THIS module's marked-point evaluation of
    E^(0)(t) against the independent classical Gamma_1(6) route
        E_ref = -( (3/2) sqrt3 L(chi_{-3},2) + I(1, f3; q_C) ).
    Returns (E_marked, E_ref, agree_digits)."""
    ones = (mp.mpf(1), mp.mpf(1), mp.mpf(1))
    diag = {}
    E_marked = predict_E(t, ones, diag=diag)
    qC = diag["qC"]
    Nq = Nq_for(qC, extra=18)
    f3 = _f3_qcoeffs(Nq)
    If3 = mp.mpc(0)
    qn = qC
    for n in range(1, Nq + 1):
        If3 += f3[n] * qn / (n * n)
        qn *= qC
    B0 = mp.mpf(3) / 2 * mp.sqrt(3) * L_chi_m3_2()
    E_ref = (-(B0 + If3)).real
    dd = abs(E_marked - E_ref)
    agree = int(-mp.log10(dd / max(abs(E_ref), mp.mpf(1)))) if dd > 0 else mp.mp.dps
    return E_marked, E_ref, agree


# --------------------------------------------------------------------------
# 5. Shared gate/CLI driver for the three mass configurations
# --------------------------------------------------------------------------
def masses_from_sq(masses_sq):
    """sqrt of the squared masses at the CURRENT dps (call again after any dps
    change -- an sqrt(2) frozen at low dps poisons the --check rerun)."""
    return tuple(mp.sqrt(mp.mpf(m)) for m in masses_sq)


def run_row(row, data, args):
    """Standard gate table + positive control + optional --point/--check.
    data: {"masses_sq":[...], "gate_points":[{"t_num","t_den","oracle_mid",
    "oracle_digits","source",...}], "paper_literals":[...]}.
    Called by the row scripts AFTER mp.mp.dps has been set in main()."""
    import time
    t0 = time.time()
    masses_sq = data["masses_sq"]
    masses = masses_from_sq(masses_sq)
    label = ",".join(str(m) for m in masses_sq)
    print("=" * 78)
    print(f"ROW {row}: two-loop sunrise (m^2 = {label}) eps^0 -- "
          f"uniform per-puncture Kronecker-eMPL closed form")
    print(f"dps = {mp.mp.dps}")
    print("=" * 78)

    # positive control (independent classical route, equal mass)
    tc = mp.mpf(-3)
    gate_failures = []   # fail-closed: below-bar => raise at end
    Em, Er, agree = equal_mass_control(tc)
    ctrl_target = min(mp.mp.dps - 8, 999)
    ctrl_ok = agree >= min(30, ctrl_target)
    if not ctrl_ok:
        gate_failures.append(f"positive control t=-3 agree {agree}d < bar "
                             f"{min(30, ctrl_target)}d")
    print(f"positive control @ t={float(tc)} (1,1,1): marked-point vs Gamma_1(6) "
          f"classical route agree {agree} digits "
          f"[{'PASS' if ctrl_ok else 'FAIL'}]")

    worst = None
    worst_rel = mp.mpf(0)
    results = []
    maxcert = mp.mpf(0)
    for gp in data["gate_points"]:
        t = mp.mpf(gp["t_num"]) / mp.mpf(gp["t_den"])
        diag = {}
        J_pred = predict_J(t, masses, diag=diag)
        oracle = mp.mpf(gp["oracle_mid"])
        rel = abs(J_pred - oracle) / abs(oracle)
        d = int(-mp.log10(rel)) if rel > 0 else mp.mp.dps
        cap = min(gp["oracle_digits"], mp.mp.dps)
        results.append((t, J_pred, oracle, d, cap, gp))
        worst = d if worst is None else min(worst, d)
        worst_rel = max(worst_rel, rel)
        maxcert = max(maxcert, diag["err_J"])
        print(f"  t = {mp.nstr(t, 8):>10}  J_pred = {mp.nstr(J_pred, 24)}")
        print(f"    vs AMFlow oracle ({gp['source']}): agree {d} digits "
              f"(cap: {cap} = min(oracle {gp['oracle_digits']}d, dps))  "
              f"|q_C| = {mp.nstr(abs(diag['qC']), 5)}  Nq = {diag['Nq']}")
    # Fail-closed gate: STRICTLY-GREATER bar (worst > 30, i.e.
    # >= 31 truncated digits) AND the /2 directed-rounding slack at rel level
    # (worst rel < 10^-30/2, the same slack Nq_for budgets).  A plain
    # 'worst >= 30' integer compare would let a 1e-30 mutation land at EXACTLY
    # 30 agree-digits and pass with rc=0.
    bar_rel = mp.mpf(10) ** (-30) / 2
    gate_ok = (worst > 30) and (worst_rel < bar_rel)
    if not gate_ok:
        gate_failures.append(
            f"WORST gate agree {worst}d not STRICTLY > 30d, or rel error "
            f"{mp.nstr(worst_rel, 3)} not < 10^-30/2 = {mp.nstr(bar_rel, 3)}")
    print(f"\n  WORST gate agreement over {len(results)} stored held-out points: "
          f"{worst} digits  [{'PASS' if gate_ok else 'FAIL'}: strict bar "
          f"> 30 d AND worst-rel {mp.nstr(worst_rel, 3)} < 10^-30/2, "
          f"directed-rounding slack]")
    print(f"  [certified] q-series truncation error |Delta J| <= "
          f"{mp.nstr(maxcert, 3)} at every gate point (per-word proven tail "
          f"bound C|q_C|^(Nq+1)/(1-|q_C|), C = (2 pi)^3 zeta(2), enforced "
          f"< 10^-(dps+12) by Nq_for, propagated by the l1-norm of the exact "
          f"prefactors; BOUND, not estimate)")

    for lit in data.get("paper_literals", []):
        t = mp.mpf(lit["t_num"]) / mp.mpf(lit["t_den"])
        E = predict_E(t, masses)
        ref = mp.mpf(lit["E0"])
        nprint = len(lit["E0"].replace("-", "").replace(".", ""))
        rel = abs(E - ref) / abs(ref)
        d = min(int(-mp.log10(rel)) if rel > 0 else mp.mp.dps, nprint, mp.mp.dps)
        lit_ok = d >= min(nprint, mp.mp.dps) - 2
        if not lit_ok:
            gate_failures.append(
                f"paper literal E^(0) @ t={mp.nstr(t, 6)} agree {d}d < bar "
                f"{min(nprint, mp.mp.dps) - 2}d")
        print(f"  paper literal E^(0) @ t={mp.nstr(t, 6)} [{nprint} printed digits, "
          f"{lit['where'].split('(')[0].strip()}]: agree {d}d "
          f"[{'PASS' if lit_ok else 'FAIL'}]")

    for ps in args.point:
        t = mp.mpf(mp.mpf(ps.split('/')[0]) / mp.mpf(ps.split('/')[1])) if '/' in ps \
            else mp.mpf(ps)
        diag = {}
        try:
            E = predict_E(t, masses, diag=diag)
            print(f"\n--point t = {mp.nstr(t, 12)}:")
            print(f"  E^(0) = J/psihat1 = {mp.nstr(E, mp.mp.dps)}")
            print(f"  psihat1 = {mp.nstr(diag['psihat1'], 30)}")
            print(f"  J^(0)  = {mp.nstr(diag['psihat1'] * E, mp.mp.dps)}")
            print(f"  certified bound: |Delta E| <= {mp.nstr(diag['err_E'], 3)}, "
                  f"|Delta J| <= {mp.nstr(diag['psihat1'] * diag['err_E'], 3)}  "
                  f"(q-series truncation; BOUND, not estimate)")
        except ValueError as e:
            print(f"\n--point t = {ps}: REFUSED ({e})")
            raise SystemExit(2)  # a refusal must be rc=2, never rc=0

    wall = time.time() - t0
    print(f"\nmeasured wall time: {wall:.2f} s")

    if args.check:
        print(f"\n--check: rerunning all gate points at dps {mp.mp.dps + 60} ...")
        lo = {mp.nstr(r[0], 10): r[1] for r in results}
        old_dps = mp.mp.dps
        mp.mp.dps = old_dps + 60
        masses_hi = masses_from_sq(masses_sq)   # re-sqrt at the new dps
        okall = True
        t1 = time.time()
        for gp in data["gate_points"]:
            t = mp.mpf(gp["t_num"]) / mp.mpf(gp["t_den"])
            J_hi = predict_J(t, masses_hi)
            dd = abs(J_hi - lo[mp.nstr(t, 10)]) / abs(J_hi)
            stab = int(-mp.log10(dd)) if dd > 0 else mp.mp.dps
            ok = stab >= old_dps - 12
            okall = okall and ok
            print(f"  t = {mp.nstr(t, 8):>10}: low-dps value stable to {stab} digits "
                  f"(target ~{old_dps})  [{'ok' if ok else 'DRIFT'}]")
        mp.mp.dps = old_dps
        print(f"--check wall time: {time.time() - t1:.2f} s")
        print(f"--check verdict: {'PASS' if okall else 'FAIL'}")
        if not okall:
            gate_failures.append("--check two-precision DRIFT (see table)")
    # Fail-closed: any FAIL printed above must also set a nonzero exit code.
    if gate_failures:
        raise RuntimeError(
            "gate FAIL (below bar => nonzero exit): " + "; ".join(gate_failures))
    return worst, wall


# --------------------------------------------------------------------------
# 6. Module self-test:  python3 rowscripts_empl.py  (no data files needed)
# --------------------------------------------------------------------------
def _selftest():
    import time
    t0 = time.time()
    mp.mp.dps = 50
    ok = True

    # (a) positive control: marked-point machinery vs independent Gamma_1(6) route
    Em, Er, agree = equal_mass_control(mp.mpf(-3))
    good = agree >= mp.mp.dps - 10
    ok &= good
    print(f"[{'PASS' if good else 'FAIL'}] equal-mass control t=-3: "
          f"marked-point vs classical Gamma_1(6) agree {agree}d (dps {mp.mp.dps})")

    # (b) z1=z2 degeneracy grouping (rows 9/11: m1=m2)
    cv = fcurve(-3, *masses_from_sq([1, 1, 4]))
    zm = marked_z_mults(cv)
    good = sorted(m for _, m in zm) == [1, 2]
    ok &= good
    print(f"[{'PASS' if good else 'FAIL'}] (1,1,4) degeneracy: multiplicities "
          f"{sorted(m for _, m in zm)} == [1, 2]")

    # (c) Abel-Jacobi branch repair active: lattice relation z1+z2+z3 = 1
    for msq, t in ([1, 1, 4], -1), ([1, 1, 2], -1), ([1, 2, 3], -2):
        cv = fcurve(t, *masses_from_sq(msq))
        s = cv["z1F"] + cv["z2F"] + cv["z3F"]
        good = abs(s - 1) < mp.mpf(10) ** (-(mp.mp.dps - 8))
        ok &= good
        print(f"[{'PASS' if good else 'FAIL'}] z1+z2+z3 = 1 at m^2={msq}, t={t}: "
              f"|sum-1| = {mp.nstr(abs(s - 1), 3)}")

    # (d) generic point all-distinct
    cv = fcurve(-3, *masses_from_sq([1, 2, 3]))
    good = len(marked_z_mults(cv)) == 3
    ok &= good
    print(f"[{'PASS' if good else 'FAIL'}] (1,2,3) generic: 3 distinct marked points")

    # (e) domain refusals: t > 0 rejected; |q_C| >= QMAX rejected
    try:
        predict_E(mp.mpf(1), masses_from_sq([1, 1, 2]))
        good = False
    except ValueError:
        good = True
    ok &= good
    print(f"[{'PASS' if good else 'FAIL'}] refusal at t=+1 (only Euclidean t<0 wired)")
    try:
        check_domain(mp.mpf(-1), mp.mpf("0.95"))
        good = False
    except ValueError:
        good = True
    ok &= good
    print(f"[{'PASS' if good else 'FAIL'}] refusal at |q_C| = 0.95 >= QMAX = {mp.nstr(QMAX, 2)}")

    # (f) control near the soft point and deep Euclidean (domain edges)
    for tt, floor in (("-0.001", mp.mp.dps - 10), ("-100", mp.mp.dps - 12)):
        _, _, ag = equal_mass_control(mp.mpf(tt))
        good = ag >= floor
        ok &= good
        print(f"[{'PASS' if good else 'FAIL'}] equal-mass control t={tt}: {ag}d (floor {floor})")

    print(f"selftest wall time: {time.time() - t0:.2f} s")
    print(f"SELFTEST {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(_selftest())
