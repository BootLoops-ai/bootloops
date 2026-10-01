# lockpick ring15 member — characters, Dirichlet L, Chowla-Selberg periods, CM periods.
"""ring15_core: conductor-15 / Q(sqrt(-15)) CM constant ring -- core pieces.

Derivations (every identity below is re-verified numerically at runtime):

[Lerch]  For primitive odd quadratic chi mod q (q=|D|, D<0 fundamental):
  L(s,chi) = q^{-s} sum_a chi(a) zeta_H(s, a/q);
  zeta_H'(0,x) = logGamma(x) - (1/2)log(2pi);  zeta_H(0,x) = 1/2 - x;
  L(0,chi) = -(1/q) sum_a a*chi(a) = 2h/w  (Dirichlet class number formula).
  =>  sum_a chi(a) logGamma(a/q) = L'(0,chi) + (2h/w) log q.            (*)

[CS period -- our convention, exponent w/(4h) from (*)]
  Omega_D := (2 pi |D|)^{-1/2} * [ prod_a Gamma(a/|D|)^{chi_D(a)} ]^{w/(4h)}
  D=-4 (h=1,w=4):  Omega_{-4} = Gamma(1/4)^2/(4 pi^{3/2}) = |eta(i)|^2.
  Independent control: lemniscate constant varpi = pi/agm(1,sqrt2) (Gauss AGM,
  no Gamma functions): (Omega_{-4}/(varpi/pi))^2 = 1/2 exactly.

[h=2 pinning, D=-15]  chi_{-15} = chi_{-3}*chi_5; +1 on {1,2,4,8}, -1 on
  {7,11,13,14} mod 15.  h=2, w=2 => exponent 1/4. Reduced forms [1,1,4] and
  [2,1,2]; CM points tau_1=(-1+sqrt(-15))/2, tau_2=(-1+sqrt(-15))/4.
  PINNED normalization (PSLQ on logs, height cap 10^4, dps 120 AND 210,
  identical hit [2,2,-1,-1,0] on {log r, log2, log3, log5, log pi}):
      prod_{j=1,2} sqrt(Im tau_j)*|eta(tau_j)|^2 = (sqrt(15)/2) * Omega_{-15}^2.
  Unified empirical CS normalization, verified on ALL of D=-3,-4,-15:
      prod_{j=1}^{h} sqrt(y_j)|eta(tau_j)|^2 = (sqrt|D|/2)^{h/2} * Omega_D^h
  (D=-4: factor=1, Omega_{-4}=|eta(i)|^2;  D=-3: PSLQ hit [-4,-2,1] => factor
   (3/4)^{1/4};  D=-15: the relation above).
  Gross unit factor between the two individual CM periods (PSLQ, both dps,
  hit [-3,-1,0,0,0] on {log(cm1/cm2), log phi, ...}):
      cm1/cm2 = phi^{-1/3},  phi = (1+sqrt(5))/2.
"""
import mpmath as mp

CHI15 = {1: 1, 2: 1, 4: 1, 8: 1, 7: -1, 11: -1, 13: -1, 14: -1}


def chi15(n): return CHI15.get(n % 15, 0)
def chi3(n):  return [0, 1, -1][n % 3]
def chi5(n):  return {1: 1, 4: 1, 2: -1, 3: -1}.get(n % 5, 0)
def chi4(n):  return {1: 1, 3: -1}.get(n % 4, 0)


# D -> (chi, modulus q, class number h, number of units w); h,w for D<0 only.
CHARS = {-15: (chi15, 15, 2, 2), -3: (chi3, 3, 1, 6),
         -4: (chi4, 4, 1, 4), 5: (chi5, 5, None, None)}


def dirichlet_L(s, D):
    """L(s, chi_D) via Hurwitz zeta; at s=1 (termwise pole, sum chi = 0 cancels
    it) use the digamma formula L(1,chi) = -(1/q) sum_a chi(a) psi(a/q)."""
    chi, q, _, _ = CHARS[D]
    if s == 1:
        return -mp.fsum(chi(a) * mp.digamma(mp.mpf(a) / q)
                        for a in range(1, q) if chi(a)) / q
    return mp.power(q, -s) * mp.fsum(
        chi(a) * mp.zeta(s, mp.mpf(a) / q) for a in range(1, q) if chi(a))


def dirichlet_Lprime0(D):
    """L'(0, chi_D) = sum chi(a) zeta_H'(0,a/q) - log(q) sum chi(a)(1/2 - a/q)."""
    chi, q, _, _ = CHARS[D]
    t = mp.fsum(chi(a) * mp.zeta(0, mp.mpf(a) / q, 1)
                for a in range(1, q) if chi(a))
    t -= mp.log(q) * mp.fsum(chi(a) * (mp.mpf(1) / 2 - mp.mpf(a) / q)
                             for a in range(1, q) if chi(a))
    return t


def lerch_residual(D):
    """Residual of (*): must be ~0 at working dps. Validates the CS exponent."""
    chi, q, h, w = CHARS[D]
    s = mp.fsum(chi(a) * mp.loggamma(mp.mpf(a) / q)
                for a in range(1, q) if chi(a))
    return s - dirichlet_Lprime0(D) - mp.mpf(2 * h) / w * mp.log(q)


def gamma_quotient(D):
    """prod_{a=1}^{q-1} Gamma(a/q)^{chi_D(a)}."""
    chi, q, _, _ = CHARS[D]
    return mp.exp(mp.fsum(chi(a) * mp.loggamma(mp.mpf(a) / q)
                          for a in range(1, q) if chi(a)))


def Omega(D):
    """Chowla-Selberg period, convention in module docstring."""
    chi, q, h, w = CHARS[D]
    return mp.power(2 * mp.pi * q, mp.mpf(-1) / 2) * \
        mp.power(gamma_quotient(D), mp.mpf(w) / (4 * h))


def eta(tau):
    """Dedekind eta, pentagonal-number series, exact to working dps."""
    q = mp.exp(2j * mp.pi * tau)
    s, n = mp.mpc(1), 1
    tol = mp.mpf(10) ** (-mp.mp.dps - 10)
    while True:
        t1 = (-1) ** n * q ** (n * (3 * n - 1) // 2)
        t2 = (-1) ** n * q ** (n * (3 * n + 1) // 2)
        s += t1 + t2
        if abs(t1) < tol:
            break
        n += 1
    return mp.exp(mp.pi * 1j * tau / 12) * s


def cm_period_eta(a, b, c):
    """sqrt(Im tau)*|eta(tau)|^2 at the CM point of reduced form [a,b,c]."""
    D = b * b - 4 * a * c
    tau = (-b + mp.sqrt(mp.mpc(D))) / (2 * a)
    return mp.sqrt(tau.imag) * abs(eta(tau)) ** 2


def disc4_controls():
    """Positive controls for the whole construction (disc -4, lemniscatic)."""
    out = {}
    out['catalan_diff'] = dirichlet_L(2, -4) - mp.catalan
    varpi = mp.pi / mp.agm(1, mp.sqrt(2))      # AGM: Gamma-free lemniscate const
    out['varpi_agm'] = varpi
    out['lemn_ratio2_minus_half'] = (Omega(-4) / (varpi / mp.pi)) ** 2 - mp.mpf(1) / 2
    out['omega4_minus_eta_at_i'] = Omega(-4) - cm_period_eta(1, 0, 1)
    out['lerch_resid_m4'] = lerch_residual(-4)
    out['lerch_resid_m3'] = lerch_residual(-3)
    out['lerch_resid_m15'] = lerch_residual(-15)
    # real-character control: L(1,chi_5) = 2*log(phi)/sqrt(5)
    phi = (1 + mp.sqrt(5)) / 2
    out['L5_1_minus_2logphi_sqrt5'] = dirichlet_L(1, 5) - 2 * mp.log(phi) / mp.sqrt(5)
    return out
