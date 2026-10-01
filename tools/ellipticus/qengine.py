"""ellipticus.qengine — Kronecker-Eisenstein / tau-word q-series engine
(engine q: the GiNaC-class fast path) with certified tail bounds.

COMPOSES (bundled engines, path-imported, NOT rewritten):
  * vendor/kronecker.py (Kronecker-Eisenstein engine) — g^(n)(z,tau) (three independent
    routes: qbar / theta / eisenstein), omega_k, Eisenstein e2/e4, eta-forms;
  * vendor/itint.py — AGW tau-iterated integrals with
    shuffle / tangential-base-point regularization at the cusp (the L-letter
    q-poly algebra);
  * vendor/rowscripts_empl.py — BMSW sunrise rows machinery with
    the PROVEN q-tail bound (|b_n| <= (2pi)^3 zeta(2) n^2), runtime-certified
    truncation Nq_for, QMAX refusal, value-level bound propagation
    (predict_E / predict_J / equal_mass_control).

Certification layer added here:
  * g_tail_bound(n, z0, tau, J): PROVEN geometric-polynomial tail bound for
    the eq.46 qbar series of g^(n) after fundamental-cell reduction —
    for the reduced point z0 (0 <= Im z0 <= Im tau / 2, the strip the
    packaged _g_qbar_strip always folds into) both ELi pieces have
    |term(j,k)| <= max(|w|, |w^-1 qbar^k|)-controlled decay with per-index
    base  rho = exp(-2 pi d),  d = min(Im z0, Im tau - Im z0)  (Im z0 = 0:
    the series is still summed against |qbar|^(jk); base |qbar|);
    |coeff of the (j,k) block| <= (2 pi)^n / (n-1)! * k^(n-1).  Summing the
    double tail j+k > J geometrically gives an EXPLICIT bound returned as an
    mp value (see manual for the inequality chain).  Fail-closed: if the
    bound cannot be pushed below 10^-(dps+guard) at the packaged J the API
    raises rather than truncating silently (QMAX-refusal pattern).
  * two-precision self-gate in ellipticus.evaluate wraps every verb.
"""
import sys
import mpmath as mp

# The engines ship bundled in vendor/ beside this package. Override with
# ELLIPTICUS_KRON_DIR / ELLIPTICUS_ROWS_DIR (former names EMPL_EVAL_KRON_DIR /
# EMPL_EVAL_ROWS_DIR still read) to point at live copies.
import os as _os
from . import env as _env
_VENDOR = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "vendor")
_KRON_DIR = _env("ELLIPTICUS_KRON_DIR", "EMPL_EVAL_KRON_DIR", _VENDOR)
_ROWS_DIR = _env("ELLIPTICUS_ROWS_DIR", "EMPL_EVAL_ROWS_DIR", _VENDOR)
for _d in (_KRON_DIR, _ROWS_DIR):
    if _d not in sys.path:
        sys.path.insert(0, _d)

import kronecker as K            # noqa: E402  (vendored engine)
import itint as IT               # noqa: E402  (vendored tau-words)
import rowscripts_empl as RE     # noqa: E402  (vendored sunrise rows)

from .curve import exact, sunrise_frame  # noqa: E402


# ----------------------------------------------------------- certified tails
def g_tail_bound(n, z, tau, J):
    """PROVEN tail bound for the qbar double series of g^(n) (eq. 46 of
    arXiv:1907.01251; ELi conventions of arXiv:1704.08895, as implemented in the vendored
    kronecker._g_qbar_cell) truncated at J terms per index, for the
    fundamental-cell representative z0 of z.

    Bound chain (manual has the derivation): with qb = exp(2 pi i tau),
    w = exp(2 pi i z0), a = |qb|, u = max(|w|, 1/|w|) (u <= a^(-1/2) after
    the vendored lower-half-strip fold),
      |ELi_{0,1-n} term (j,k)| <= u^j k^(n-1) a^(jk),
    tail(j > J or k > J) <= 2 * sum_{j>J} u^j a^j / (1-a)-type geometric
    majorants; we use the crude but explicit
      T(J) = 2 * (2 pi)^n / (n-1)! * [ S_j + S_k ],
      S_j = sum_{j>J} (u a)^j * C_k,   C_k = sum_{k>=1} k^(n-1) a'^k
            with a' = a (polylog closed form Li_{1-n}(a)),
      S_k = C_j' * sum_{k>J} k^(n-1) a^k   (Lerch tail, bounded by
            (J+1)^(n-1) a^(J+1) / (1-a)^n * n!),
    everything evaluated in mp at the ambient precision.  Requires
    u * a < 1 (true in the folded strip for Im tau > 0); raises otherwise."""
    n = int(n)
    if n < 1:
        return mp.mpf(0)
    z0, _ = K._reduce_z(mp.mpc(z), mp.mpc(tau))
    imt = mp.im(mp.mpc(tau))
    # fold to lower half-strip as the vendored strip evaluator does
    if mp.im(z0) > imt / 2:
        z0 = -z0 + mp.mpc(tau)
    a = abs(K.qbar(tau))
    u = mp.exp(2 * mp.pi * abs(mp.im(z0)))
    ua = u * a
    if not ua < 1:
        raise ValueError("g_tail_bound: u*a >= 1 (point too close to the "
                         "cell edge for the qbar route) — use route='theta'")
    pref = 2 * (2 * mp.pi) ** n / mp.factorial(max(n - 1, 1))
    # |term(j,k)| <= u^j k^(n-1) a^(jk)  (both ELi pieces; |w^j| <= u^j).
    # S_j (j > J, all k >= 1):
    #   sum_k k^(n-1) (a^j)^k <= a^j * D,  D = Li_{1-n}(a)/a  (x <= a since
    #   j >= 1; termwise k^(n-1) x^k <= x * k^(n-1) a^(k-1))
    #   => S_j <= D * sum_{j>J} (u a)^j = D * (ua)^(J+1)/(1-ua).
    D = mp.polylog(1 - n, a) / a
    Sj = D * ua ** (J + 1) / (1 - ua)
    # S_k (k > J, all j >= 1):
    #   sum_j u^j a^(jk) = (u a^k)/(1 - u a^k) <= u a^k / (1 - u a^(J+1)),
    #   sum_{k>J} k^(n-1) a^k <= (J+1)^(n-1) a^(J+1) / (1 - a e^((n-1)/(J+1)))
    #   [k^(n-1) = (J+1)^(n-1) e^((n-1) ln(k/(J+1))) <= (J+1)^(n-1)
    #    e^((n-1)(k-(J+1))/(J+1))], requires a e^((n-1)/(J+1)) < 1.
    ae = a * mp.exp(mp.mpf(n - 1) / (J + 1))
    if not ae < 1:
        raise ValueError("g_tail_bound: a e^((n-1)/(J+1)) >= 1 — raise J")
    Sk = (u / (1 - u * a ** (J + 1))) * mp.mpf(J + 1) ** (n - 1) \
        * a ** (J + 1) / (1 - ae)
    return pref * (Sj + Sk)


def g_n_certified(n, z, tau, dps, extra=12):
    """g^(n)(z, tau) via the vendored qbar route with a PROVEN tail bound
    pushed below 10^-(dps+extra); returns (value, bound).  Falls back to
    (and cross-checks against) the theta route when the bound cannot
    converge (cell-edge points); fail-closed on disagreement."""
    with mp.workdps(dps + 30):
        target = mp.mpf(10) ** (-(dps + extra))
        z0, _ = K._reduce_z(mp.mpc(z), mp.mpc(tau))
        Jauto = K._qbar_J_for(z0, mp.mpc(tau), extra=30)
        J = min(max(Jauto, 40), K._QBAR_J_CAP)
        try:
            bnd = g_tail_bound(n, z, tau, J)
            grow = 0
            while bnd > target and grow < 12:
                J = int(J * 1.5) + 10
                if J > 20000:
                    break
                bnd = g_tail_bound(n, z, tau, J)
                grow += 1
        except ValueError:
            bnd = None
        if bnd is not None and bnd <= target:
            val = K.g_n(n, z, tau, route="qbar", J=J)
            return val, bnd
        # theta-route fallback: no series tail; certify by route agreement
        v1 = K.g_n(n, z, tau, route="theta")
        v2 = K.g_n(n, z, tau, route="qbar")
        d = abs(v1 - v2)
        assert d < mp.mpf(10) ** (-(dps - 5)) * max(abs(v1), mp.mpf(1)), \
            f"qbar/theta route disagreement {mp.nstr(d, 4)}"
        return v1, d


# ------------------------------------------------------- sunrise front door
def sunrise_E0(t, masses_sq, dps):
    """eps^0 unequal-mass two-loop sunrise E^(0)(t) (BMSW normalization) via
    the vendored PROVEN-bound machinery (rowscripts_empl.predict_E: zero free
    parameters, runtime-certified q-truncation, QMAX refusal).
    masses_sq: exact rationals (m1^2, m2^2, m3^2).  Euclidean t < 0.
    Returns (value, certified_bound, diag)."""
    t = exact(t, 't')
    ms = [exact(m, 'mass^2') for m in masses_sq]
    with mp.workdps(dps):
        masses = RE.masses_from_sq(
            [mp.mpf(m.numerator) / mp.mpf(m.denominator) for m in ms])
        diag = {}
        val = RE.predict_E(mp.mpf(t.numerator) / mp.mpf(t.denominator),
                           masses, diag=diag)
        return val, diag.get('err_E'), diag


def sunrise_J(t, masses_sq, dps):
    """J[1,1,1]^(0)(t) = E^(0) * psihat1 (gated against the bundled
    rows-09/10/11 AMFlow reference points) via the vendored proven-bound
    machinery (rowscripts_empl.predict_J).  Returns (value, certified_bound, diag)."""
    t = exact(t, 't')
    ms = [exact(m, 'mass^2') for m in masses_sq]
    with mp.workdps(dps):
        masses = RE.masses_from_sq(
            [mp.mpf(m.numerator) / mp.mpf(m.denominator) for m in ms])
        diag = {}
        val = RE.predict_J(mp.mpf(t.numerator) / mp.mpf(t.denominator),
                           masses, diag=diag)
        return val, diag.get('err_J'), diag


def sunrise_equal_mass_control(t, dps):
    """positive control: same machinery at (1,1,1) vs the independent
    classical Gamma_1(6) route (vendored equal_mass_control)."""
    t = exact(t, 't')
    with mp.workdps(dps):
        return RE.equal_mass_control(
            mp.mpf(t.numerator) / mp.mpf(t.denominator))


# ------------------------------------------------------------ tau-words
def tau_word_frame(qC, z1, z2, z3, Nq):
    """AGW tau-word frame (vendored itint.Frame)."""
    return IT.Frame(qC, z1, z2, z3, Nq)


def tau_words(words_coeffs, frame):
    """value of sum_i coeff_i * F(word_i) via the vendored shuffle-
    regularized qbar engine (itint.itint)."""
    return IT.itint(words_coeffs, frame)
