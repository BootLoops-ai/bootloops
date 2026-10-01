# Leviathan — chi-drop Landau singularities (a Julia implementation of arXiv:2606.29612)

A Julia implementation (over OSCAR) of the Euler-characteristic-drop method for
Landau singularities of Feynman integrals, from the papers of Vsevolod Chestnov,
Giulio Crisanti and Mathieu Giroux, "Landau's Leviathans" (arXiv:2606.29612) — from
which this package also borrows its name, with thanks — and of Vsevolod Chestnov and
Giulio Crisanti, SPQR (arXiv:2511.14875), whose companion-matrix minimal-polynomial
construction the exact candidate route computes. All algorithmic credit is theirs;
please cite both papers. No code was copied or translated. The authors' public
Mathematica code (landau_codes.m in Landau-s-Leviathans by Chestnov, Crisanti and
Giroux, and euler_characteristic.m in DiscKosky by Giulio Crisanti, Luke Lippstreu,
Andrew J. McLeod and Maria Polackova, both on GitHub) was read as the reference for
five details the papers leave to the code: the Type-1.2 non-degeneracy test, the
"Diophantine" locus sampling (lowest-degree variable, 400 tries), the cut-edge
regulated ideal, normalization by the least common multiple of denominators, and
dehomogenization in the last kinematic variable (PATCHES.md, "What follows the
authors' code", gives our file:line and theirs). Their worked examples and published
singularity lists validate this engine.
Not affiliated with or endorsed by the authors. The χ = number-of-master-integrals
principle behind the method is due to Lee and Pomeransky (JHEP 11 (2013) 165), Bitoun,
Bogner, Klausen and Panzer (Lett. Math. Phys. 109 (2019) 497) and Agostini, Fevola,
Sattelberger and Telen (Commun. Num. Theor. Phys. 18 (2024) 327); principal Landau
determinants to Fevola, Mizera and Telen (Comput. Phys. Commun. 303 (2024) 109278;
Phys. Rev. Lett. 132 (2024) 101601); Gröbner bases by msolve (Berthomieu, Eder,
Safey El Din) inside OSCAR.

Landau singularities of a Feynman family = kinematic loci where the Euler characteristic
chi of the regulated Lee–Pomeransky critical-point ideal drops (Eq. 5 of the paper).
See PATCHES.md for the full provenance and `LICENSE` for the license (MIT).

## Stack
Julia 1.10 + Oscar 1.7.3 (Project.toml in this directory; Leviathan 1.0):

    julia +1.10 --project=. leviathan.jl <cmd> <input> [seed]

Commands: `chi` (per-sector chi + regulated chi + additivity), `candidates`
(singularity factors from elimination), `diagnose` (4 Appendix-A checks), `landau`
(full pipeline incl. Type-2.2 chi-drop verification of every candidate).

Input file (key=value): `G = ...` (or `U =`/`F =`), `kin = s, mm`, `x = x1, x2`,
optional `cut = 1,2` (cut edges: plain dG generators, excluded from saturation).
Example: examples/bubble.in.

## Method/implementation map
- chi counts (Eq. 12): sector ideal J_S = <dG_S, 1 - x0 G_S>, regulated ideal
  I = <nu_i G - (d/2) x_i dG_i, 1 - x0 (prod x) G>; random 31-bit F_p values for
  kinematics and (nu, d); GB via Oscar groebner_basis_f4 (msolve F4 in-process);
  chi = #standard monomials by own DRL staircase counter (src/staircase.jl) —
  interval-sliced recursion whose cost scales with the staircase's threshold structure,
  not box volume (counts wider than Int come back as Int128);
  `nothing` = positive-dimensional (Type-1.1 Indeterminate).
- Candidates (Eq. 8): EXACT route over K = Q(kinematics): one DRL GB of the ideal over K,
  then per variable the monic minimal polynomial of x_i in K[x]/I by normal forms +
  K-linear solve (the SPQR companion-matrix object, without FiniteFlow); candidates =
  irreducible factors of the LCM of denominators of the monic coefficients f_k.
  Routes: regulated (generic rational nu, d — catches second-type, e.g. bubble s=0)
  and nu=0 per sector. Driver unions both.
- Reconstruction route (src/reconstruct.jl `landau_reconstruct`): the same monic
  minpoly, sampled over (prime, kinematic point) and rationally reconstructed
  (src/ratrec.jl dense-ansatz nullspace + CRT) — for families where the exact
  function-field GB is out of reach. NOTE input requirement: G must be jointly
  homogeneous in the kin list and the DIMENSIONFUL scale (e.g. m2) must be the LAST kin
  entry — the route dehomogenizes kin[end]=1 and rehomogenizes factors by total degree
  (setting m2=1 in the input silently breaks rehomogenization).
- Verification (Type 2.2): chi_regulated constrained to each candidate locus
  (random F_p point on l=0, "Diophantine" mode, <=400 tries) must drop below generic.
- Diagnostics: type1.1 (sector dim), type1.2 noDegeneracyQ (kin-var/x-var role swap,
  2 random values, GCD of univariate eliminants mod p < 2^29), type2.1 (chi additivity),
  type2.2 (candidate chi-drop).

## Files
src/{Leviathan.jl,primes.jl,parse.jl,staircase.jl,family.jl,chi.jl,elimination.jl,
diagnostics.jl,landau.jl,ratrec.jl,reconstruct.jl}; leviathan.jl (CLI); chi_at.jl
(chi at specific kinematic points, JSON in/out); builder.jl (graph->G via PLD's getUF —
needs the external PLD package, set ENV["PLD_PATH"]); test_leviathan.jl
(+ test_smoke.jl, test_sunrise.jl, test_reconstruct.jl, test_staircase.jl — the last
is engine-free: it runs on bare Julia with no Oscar install).

## Footguns (measured)
- Singular-backed ops (eliminate/reduce/GB over GF(p)) need p < 2^29; f4 takes 31-bit.
- `leading_exponent_vector(g)` without explicit ordering is WRONG for f4 output; we
  compute DRL leading exponents ourselves.
- Per-variable `eliminate()` over Q(s) on regulated ideals hangs already at sunrise-top
  (>90 s); the minpoly route does the same job in seconds. Don't "fix" back.
- msolve F4 / Singular normal_form are NOT thread-safe: threaded sampling gives silently
  wrong coefficients. Sampling is serial per process; parallelize across processes.

## Validation
- Bubble: exact {s, mm, s-4mm}, chi 3=1+1+1; diagnostics 1.1/2.1/2.2 PASS, 1.2
  correctly FAILS flagging s (second-type, nu=0-degenerate).
- Sunrise: chi_top=4, chi 7=4+1+1+1, factors {s, mm, s-mm, s-9mm} (the paper's
  worked example).
- Reconstruction route: bubble and sunrise reproduce the exact route, and the sunrise
  its standard threshold/pseudo-threshold singular set (test_reconstruct.jl).
- Not implemented (documented frontier): SPQRDet (Faddeev–LeVerrier characteristic-polynomial
  route), 4d-restriction parametrisation, the SPQR/FiniteFlow Macaulay-matrix route.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).
