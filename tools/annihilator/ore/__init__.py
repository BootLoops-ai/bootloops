"""annihilator.ore — Ore-algebra mod-p operator kit.

Import members explicitly — the GF(p) members need python-flint:

  ore_ops        — pure-sympy Ore over Q(x): mul/apply_to_series/lclm/
                   right_divide/formal_solution (no flint)
  lclm_modp      — pairwise LCLM over GF(p)[x] via Sylvester nullspace,
                   ascending-order minimality search with deg_cap
                   escalation, every stage verified by right-divisibility
                   (flint)
  ore_rdiv_modp  — Ore pseudo right-division over GF(p)[x] (flint); lifted
                   to the package top level as annihilator.ore_rightdiv_modp
  known_factors  — exact differential-operator factors transcribed from
                   arXiv:1110.1705 (selftest fixture for lclm_modp; carries
                   its own exact annihilation self-checks, no flint)

Law: validate right-division on known-good factors FIRST (S=0) before
trusting S!=0 negatives.
"""
