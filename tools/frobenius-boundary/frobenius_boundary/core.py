#!/usr/bin/env python3
"""core.py — family-agnostic API for the Frobenius-branch boundary machinery.

Typical use (see cli.py for the config-driven driver):
    basis = de_core.load_masters(masters_path); ai = de_core.alpha_int(basis)
    D, G, fit_err = build_poly_DE(kira_m, basis, var='m2', eps='1/101', prec=80)
    spec = spectrum(D, G, ai, eps='1/101')          # exponents at var=inf
    fams, groups = classify([spec1['eigvals'], ...], eps_list)
    strata = exclude_strata(groups, EVs, eps_vals)  # survivors
    Nn, diag = branch_series(D, G, lam, v, 30, ai=ai, eps='1/101')
    Phi = phi_matrix(196, [(lam, Nn), ...], ai, dd)
Convention: M_k(var) ~ var^(L*d/2 + ai_k + lambda_j), lambda_j = eig(D1).
"""
import os
from fractions import Fraction
import mpmath as mp
mp.mp.dps = max(mp.mp.dps, 50)      # set FIRST
os.environ.setdefault('OPENBLAS_NUM_THREADS', '2')   # mandatory cap

from . import de_core, poly_de, infinity, classify_mod, branches_mod
from .classify_mod import exclude_strata, trace_rule, jordan_info   # re-export
from .branches_mod import phi_matrix                                # re-export

LAST_BUILD = {}      # extras of the most recent build_poly_DE (roots, degs, N)


def eps_to_d(eps):
    """'1/101' | Fraction | float -> (dd_rat sympy, dd mp, eps mp)."""
    import sympy as sp
    fe = Fraction(eps) if not isinstance(eps, Fraction) else eps
    dd_rat = sp.Rational(4) - 2 * sp.Rational(fe.numerator, fe.denominator)
    epsv = mp.mpf(fe.numerator) / fe.denominator
    return dd_rat, mp.mpf(4) - 2 * epsv, epsv


def delta_of(ai):
    N = len(ai)
    return [[ai[l] - ai[k] for l in range(N)] for k in range(N)]


def build_poly_DE(kira_targets_m_path, masters, var='m2', s_t_numeric=None,
                  eps='1/101', prec=80, mass_slots=None, family=None,
                  verbose=False):
    """Exact polynomial DE  D(var)*dM/dvar = G(var)*M  at fixed numeric
    kinematics (already baked into the Kira table; s_t_numeric is recorded
    as provenance only).  Returns (D, G, fit_err_log10); D = mp coeff list
    (ascending), G = list of mp.matrix.  Extras (roots, degD, degG, Dc_rat,
    N) in core.LAST_BUILD."""
    basis = de_core.load_masters(masters)
    if mass_slots is None:
        mass_slots = list(range(8))          # reference-family default; configure!
    dd_rat, dd, epsv = eps_to_d(eps)
    r = poly_de.build_poly_DE_full(kira_targets_m_path, basis, dd_rat, prec,
                                   mass_slots, var=var, family=family,
                                   verbose=verbose)
    LAST_BUILD.clear()
    LAST_BUILD.update(r, s_t_numeric=s_t_numeric, eps=str(eps), var=var,
                      basis=basis)
    return r['Dc'], r['Gc'], r['fit_err_log10']


def spectrum(D, G, ai, eps='1/101', n_loops=3, n_ord=10, prec=None,
             dps_eig=60):
    """Exponent spectrum at var=inf via the D1 = R_Delta - diag(L*d/2+ai)
    eigenproblem.  Returns dict(D1, eigvals, Rj, alpha, jordan-ready data).
    Exponent of master k on branch j:  L*d/2 + ai_k + eigvals[j]."""
    prec = prec or mp.mp.dps + 40
    _, dd, epsv = eps_to_d(eps)
    N = len(ai)
    Delta = delta_of(ai)
    Rj = infinity.laurent_xA(D, G, N, n_ord, prec)
    D1 = infinity.D1_matrix(Rj, ai, dd, Delta, N, n_loops)
    ev = infinity.eig_spectrum(D1, dps_eig)
    return {'D1': D1, 'eigvals': ev, 'Rj': Rj, 'eps': str(eps), 'dd': dd,
            'n_loops': n_loops, 'ai': ai, 'Delta': Delta, 'N': N}


def classify(exponents, eps_list, maxden=64, kill_integer_eps_indep=True):
    """exponents: list of eigenvalue lists (one per eps in eps_list).
    Returns (families, groups) — rational-linear families lambda=a+b*eps with
    multiplicities, verified at ALL provided eps (the 4-eps pattern)."""
    eps_vals = [eps_to_d(e)[2] for e in eps_list]
    return classify_mod.classify(exponents, eps_vals, maxden,
                                 kill_integer_eps_indep)


def branch_series(D, G, lambda_j, v_j, n_terms, ai=None, eps='1/101',
                  n_loops=3, prec=None, eigvals=None, seed_index=None):
    """Branch solution N^(j) = x^{lambda_j} sum_n N_n x^{-n} (x = var) via the
    D_p-block recursion seeded by eigenvector v_j of D1 (D1 v = lambda_j v).
    seed_index: rank of v_j in your seed_vectors basis — REQUIRED for an
    injective branch fetch on degenerate families of genuine-D_0 systems
    (see branches_mod docstring).
    Returns (Nn dict, diag with resonance report)."""
    prec = prec or mp.mp.dps + 40
    if ai is None:
        ai = de_core.alpha_int(LAST_BUILD['basis'])
    _, dd, _ = eps_to_d(eps)
    N = len(ai)
    Delta = delta_of(ai)
    Rj = infinity.laurent_xA(D, G, N, n_terms, prec)
    Dp = infinity.build_Dp(Rj, ai, dd, Delta, N, n_terms, n_loops)
    return branches_mod.branch_series(Dp, N, lambda_j, v_j, n_terms, prec,
                                      eigvals=eigvals, seed_index=seed_index)


def seed_vectors(D1, lam, tol_exp=15):
    """Nullspace of (D1 - lam I): the Frobenius seeds of branch lam."""
    return infinity.seed_vectors(D1, lam, tol_exp)
