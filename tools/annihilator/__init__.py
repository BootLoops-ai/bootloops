#!/usr/bin/env python3
r"""annihilator — PF operator / holonomic recurrence from a series alone.

Package init: `import annihilator` exposes the full core surface
(pf_from_series, find_recurrence_fast, nullspace_mod_fast, rec_to_theta,
pretty_rec, moments_sc, P31), which lives in annihilator/annihilator.py.
The flat legacy import names ship as header-marked shims (`series_pf.py`,
`nullspace_fast/`).

Members (import explicitly; see GUIDE.md):
  annihilator/ore/       — Ore-algebra mod-p operator kit: pure-sympy Ore over
                           Q(x) (ore_ops), pairwise LCLM over GF(p)[x] via
                           Sylvester nullspace (lclm_modp), Ore pseudo
                           right-division over GF(p)[x] (ore_rdiv_modp), plus
                           the transcribed arXiv:1110.1705 operator fixtures
                           (known_factors) that drive the LCLM selftest.
  annihilator/factor/    — order-4 right-factor extraction for theta-form
                           operators (rf_lib).
  annihilator/eigenring/ — eigenring/hom-space engine that splits Jordan/
                           log-entangled right factors mod p
                           (s4_eigenring_pilot, with --selftest).

LIFTED PRIMITIVE: ore_rightdiv_modp is exported at the package top level,
guarded so environments without python-flint keep the pure numpy/sympy core
surface.
"""
from .annihilator import *  # noqa: F401,F403
from .annihilator import (P31, nullspace_mod_fast, find_recurrence_fast,  # noqa: F401
                          pf_from_series, rec_to_theta, pretty_rec, moments_sc)

# Lifted Ore primitive. Requires python-flint; core surface unaffected
# when flint is absent.
try:
    from .ore.ore_rdiv_modp import ore_rightdiv_modp  # noqa: F401
except ImportError:  # python-flint not installed
    ore_rightdiv_modp = None
