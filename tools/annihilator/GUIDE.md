# Annihilator — PF operator / holonomic recurrence from a series alone

Tool page: https://bootloops.ai/tools/annihilator.html

KIND: package (`annihilator.py` core; pure python/numpy — no Julia; member subdirs
`ore/`, `factor/`, `eigenring/` need python-flint, see MEMBERS). The flat import
names ship as header-marked import shims: `series_pf.py` and `nullspace_fast/`
(which also carries the 15/15 self-test).

PURPOSE: Find the annihilating (Picard–Fuchs / holonomic) operator of a series:
minimal P-finite recurrence and θ-form ODE with ℚ coefficients, from an exact
moment/period sequence alone. Data-route sibling of the integrand-route Griffiths–
Dwork rank probe (`pf_rank.jl`, shipped at the `tools/` root; that tool needs the integrand polynomial —
this needs ONLY the series a[0..N]).

USE-WHEN:

- Have a fast exact series generator (LGF moments, maxcut Taylor coefficients,
  diagonal coefficients) but no clean integrand, and need the PF operator /
  recurrence.
- Need the holonomic recurrence order of a long integer sequence mod p where a sympy
  nullspace is ETA-prohibitive → the GF(p) numpy kernel (10–100× over sympy;
  reference scale: order-24/deg-19 recurrence from 550 terms in 204 s; LGF
  orthorhombic order-5 in <5 s).
- Need true minimal ODE order → the ODE-order-first (s-outer, r-inner) search; the
  default unknowns-first sort inflates orders with apparent-singularity indicial
  factors.

NOT-FOR: Does not certify the operator geometrically (pair with an integrand-route
tool when the integrand exists). Kernel returns ONE nullvector, not a basis. Numpy
engine requires p < 2^31 (int64; constraint p²·dim < 2^63); default p = 2^31−1.

INVOKE: `import annihilator` (or legacy `import series_pf` / `from nullspace_fast
import nullspace_mod_fast` through the shims). Key surface: `nullspace_mod_fast(A,
p)` (int64 GF(p) RREF, one nullvector or None); `find_recurrence_fast` ((r,d)-ansatz
sweep, overdetermined solve, full-row verify); `pf_from_series` with `mode='auto'`
(mod-p s-first scan pins (r,s), then ONE exact ℚ-nullspace at that (r,s)
reconstructs coefficients — avoids full-ℚ Fraction blowup); `rec_to_theta`
(recurrence → θ-form ODE). Use `prefer='ode_order'` to force minimal ODE order.

INPUTS: exact integer/rational sequence a[0..N] (long enough: unknowns of the (r,d)
ansatz must be overdetermined); prime p < 2^31.

OUTPUTS: (order r, coeff-degree d) verdict; recurrence coefficients over ℚ; θ-form
PF operator.

GATES: Always 2-prime-confirm (r,d) before trusting; held-out verify (fit on prefix,
check on tail) is built into `find_recurrence_fast`'s full-row verify — keep it on.
Self-test battery: Fibonacci + Catalan P-finite plus the CLI selftest
exit-status gate (green run exits 0; a planted-corruption mutant exits 1;
unknown flags rejected), 15/15 PASS
(`nullspace_fast/test_nullspace_fast.py`, runs through the shim).

Battery (run from a scratch cwd): shim self-test 15/15 ALL PASS; `python3
annihilator.py --selftest` built-in selftest PASS (auto-mode ladder d=2..5,
each rung pinned to its known minimal (r,s), held-out 15 each; exits 0 only
if every rung passes, 1 otherwise).

MEMBERS (import explicitly; the GF(p) members need python-flint — without it the
core numpy/sympy surface still imports, with `annihilator.ore_rightdiv_modp = None`):

- `ore/` — Ore-algebra mod-p operator kit: `ore_ops` (pure-sympy Ore over ℚ(x):
  mul, apply_to_series, lclm, right_divide, formal_solution — no
  flint), `lclm_modp` (pairwise LCLM over GF(p)[x] via Sylvester nullspace;
  ascending-order minimality search with deg_cap escalation, every stage verified
  by right-divisibility; its `__main__` selftest builds the order-11 LCLM of the
  four arXiv:1110.1705 direct-summand operators from the `known_factors` fixture
  and gates on order 11 + all four right-division checks), `ore_rdiv_modp` (Ore
  pseudo right-division over GF(p)[x]; lifted to the package top level as
  `annihilator.ore_rightdiv_modp`; `__main__` runs a planted-factor S=0/S≠0
  control), `known_factors` (transcribed arXiv:1110.1705 operator fixtures with
  exact annihilation self-checks; the printed L^(5)_4 carries a transcription
  caveat — see its docstring).
- `factor/` — `rf_lib`: order-4 right-factor extraction for theta-form operators
  (mod-p MUM log-towers, degree-scan fit, Ore right-division verification,
  multi-prime CRT + rational reconstruction, exact-ℚ right-division, theta_mul).
- `eigenring/` — `s4_eigenring_pilot.py`: eigenring/hom-space engine that splits
  Jordan/log-entangled right factors mod p; COPAIR verification (candidates are
  applied to independently generated series, never accepted on prime-agreement or
  division alone). Self-contained `--selftest`.

Member law (S=0 first): validate right-division on known-good factors FIRST
(S=0) before trusting S≠0 negatives.

FOOTGUNS:

- Unknowns-first search returns inflated ODE orders (apparent-singularity indicial
  factors); s-first guarantees the FIRST hit is the true minimal
  θ-order.
- One nullvector ≠ nullspace: rank/degeneracy conclusions need more than this kernel.
- Single-prime (r,d) can alias; the 2-prime confirmation is mandatory, not advisory.

RELATED: `pf_rank.jl` (integrand-route sibling; mod-p Griffiths–Dwork) and `sparse_cascade.py` (sparse complement of the dense-ansatz kernel), both at the `tools/` root; the maxcut kit's Krylov/series route (`tools/maxcut/`).

CREDIT: recovering an annihilating operator from series coefficients is the 'guessing'
method of Salvy and Zimmermann's gfun [gfun] and Kauers' Guess, used at scale with
modular arithmetic by Boukraa, Hassani, Maillard, Zenine, McCoy, van Hoeij, Koutschan
and collaborators, from whose diagonal-Ising paper (Assis, Boukraa, Hassani, van Hoeij,
Maillard, McCoy [Ising11]) our operator fixtures are transcribed. The order-degree
trade-off behind the ODE-order-first search, and the apparent singularities it avoids,
are explained by Chen, Jaroschek, Kauers and Singer [CJKS]. Operator factorization and
the eigenring method follow Mark van Hoeij [vH97] (whose Maple DFactor is the reference
implementation) and Singer; LCLMs via Sylvester matrices are classical (Ore). Mod-p
linear algebra by numpy and python-flint [Arb]. Bracketed keys resolve in REFERENCES.md
at the repository root.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).
