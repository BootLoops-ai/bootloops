# Leviathan: chi-drop Landau singularities (a Julia port of the method of arXiv:2606.29612)

A Julia port (over OSCAR) of the Euler-characteristic-drop method for Landau
singularities of Feynman integrals of Vsevolod Chestnov, Giulio Crisanti and Mathieu
Giroux, "Landau's Leviathans" (arXiv:2606.29612), from which this package also takes
its name, with thanks, and of the elimination object of SPQR by Vsevolod Chestnov and
Giulio Crisanti (arXiv:2511.14875), whose companion-matrix minimal polynomial the exact
candidate route computes. All algorithmic credit is theirs; please cite both papers.

This package was written with the authors' public Mathematica code open
(landau_codes.m in Landau-s-Leviathans by Chestnov, Crisanti and Giroux;
euler_characteristic.m and irreducible_monomials.m in DiscKosky by Giulio Crisanti, Luke
Lippstreu, Andrew J. McLeod and Maria Polackova; SPQR by Chestnov and Crisanti; all on
GitHub), and parts of it are translations of that code into Julia/Oscar: the Type-1.2
non-degeneracy test, the "Diophantine" locus sampling and the cut-propagator regulated
ideal. Other parts follow its conventions: normalization by the least common multiple of
denominators, dehomogenization in the last kinematic variable, the wrapper around the
monomial count, the retry over primes and the driver's verify switch. PATCHES.md gives
our file:line and theirs for each. The candidate generation by minimal polynomials, the
monomial counting algorithm, the rational reconstruction, the parser and the input and
output were written for this package, for the objects defined in the papers. Counted
line by line on release 1.0.1, 97 of the 885 non-blank engine lines are translations of
the authors' routines and 92 more follow their code; counted by whole functions (a
function that contains a translated block counts in full), 415 of 885.

The authors' repositories carry no license file. NOTICE in this directory lists the
translated portions and states what the LICENSE in this directory covers.

Not affiliated with or endorsed by the authors. The χ = number-of-master-integrals
principle behind the method is due to Lee and Pomeransky (JHEP 11 (2013) 165), Bitoun,
Bogner, Klausen and Panzer (Lett. Math. Phys. 109 (2019) 497) and Agostini, Fevola,
Sattelberger and Telen (Commun. Num. Theor. Phys. 18 (2024) 327); principal Landau
determinants to Fevola, Mizera and Telen (Comput. Phys. Commun. 303 (2024) 109278;
Phys. Rev. Lett. 132 (2024) 101601); Gröbner bases by msolve (Berthomieu, Eder,
Safey El Din) inside OSCAR.

Landau singularities of a Feynman family = kinematic loci where the Euler characteristic
chi of the regulated Lee–Pomeransky critical-point ideal drops (Eq. 5 of the paper).
Equation numbers here and in the source files are those of arXiv:2606.29612v1, counted in
the TeX source of version 1.
See PATCHES.md for the full provenance, and `LICENSE` and `NOTICE` for the license.

## Stack
Julia 1.10 + Oscar 1.7.3 (Project.toml in this directory; Leviathan 1.0):

    julia +1.10 --project=. leviathan.jl <cmd> <input> [seed]

Commands: `chi` (per-sector chi + regulated chi + additivity), `candidates`
(singularity factors from elimination), `diagnose` (4 Appendix-A checks), `landau`
(full pipeline incl. Type-2.2 chi-drop verification of every candidate).

Input file (key=value): `G = ...` (or `U =`/`F =`), `kin = s, mm`, `x = x1, x2`,
optional `cut = 1,2` (cut edges: plain dG generators, excluded from saturation; the
cut list and its rule are those of DiscKosky's counters).
Example: examples/bubble.in.

## Method/implementation map
- chi counts: sector ideals J and J_S = <dG_S, 1 - x0 G_S> (Eqs. 10 and 11), regulated
  ideal I = <nu_i G - (d/2) x_i dG_i, 1 - x0 (prod x) G> (Eq. 3), sector sum (Eq. 12);
  random 31-bit F_p values for kinematics and (nu, d); GB via Oscar groebner_basis_f4
  (msolve F4 in-process); chi = #standard monomials by the DRL staircase counter in
  src/staircase.jl. The counting algorithm was written for this package (an
  interval-sliced recursion whose cost scales with the staircase's threshold structure,
  not box volume; counts wider than Int come back as Int128); its guards, return
  contract and three-helper interface follow DiscKosky's irreducible_monomials.m and
  SPQR's various_functions.m. `nothing` = positive-dimensional (Type-1.1 Indeterminate).
- Candidates (Eqs. 8 and 9): EXACT route over K = Q(kinematics): one DRL GB of the ideal
  over K, then per variable the monic minimal polynomial of x_i in K[x]/I by normal
  forms + K-linear solve (the SPQR companion-matrix object, without FiniteFlow);
  candidates = irreducible factors of the LCM of denominators of the monic coefficients
  f_k (this read-off follows landau_codes.m).
  Routes: regulated (generic rational nu, d; catches second-type, e.g. bubble s=0)
  and nu=0 per sector. Driver unions both.
- Reconstruction route (src/reconstruct.jl `landau_reconstruct`): the same monic
  minpoly, sampled over (prime, kinematic point) and rationally reconstructed
  (src/ratrec.jl dense-ansatz nullspace + CRT), for families where the exact
  function-field GB is out of reach. NOTE input requirement: G must be jointly
  homogeneous in the kin list and the DIMENSIONFUL scale (e.g. m2) must be the LAST kin
  entry: the route dehomogenizes kin[end]=1 and rehomogenizes factors by total degree
  (setting m2=1 in the input silently breaks rehomogenization). The last-variable
  convention follows landau_codes.m.
- Verification (Type 2.2): chi_regulated constrained to each candidate locus
  (random F_p point on l=0, "Diophantine" mode, <=400 tries; translated from DiscKosky)
  must drop below generic.
- Diagnostics: type1.1 (sector dim), type1.2 noDegeneracyQ (kin-var/x-var role swap,
  2 random values, GCD of univariate eliminants mod p < 2^29; translated from
  landau_codes.m), type2.1 (chi additivity), type2.2 (candidate chi-drop).

## Files
src/{Leviathan.jl,primes.jl,parse.jl,staircase.jl,family.jl,chi.jl,elimination.jl,
diagnostics.jl,landau.jl,ratrec.jl,reconstruct.jl}; leviathan.jl (CLI); chi_at.jl
(chi at specific kinematic points, JSON in/out); builder.jl (graph->G via PLD's getUF;
needs the external PLD package, set ENV["PLD_PATH"]); test_leviathan.jl
(+ test_smoke.jl, test_sunrise.jl, test_reconstruct.jl, test_staircase.jl; the last
is engine-free: it runs on bare Julia with no Oscar install); NOTICE, PATCHES.md.

## Footguns (measured)
- Singular-backed ops (eliminate/reduce/GB over GF(p)) need p < 2^29; f4 takes 31-bit.
- `leading_exponent_vector(g)` without explicit ordering is WRONG for f4 output; we
  compute DRL leading exponents ourselves.
- Oscar's per-variable `eliminate()` over Q(s) on regulated ideals hangs already at
  sunrise-top (>90 s); the minpoly route does the same job in seconds. Don't "fix" back.
- msolve F4 / Singular normal_form are NOT thread-safe: threaded sampling gives silently
  wrong coefficients. Sampling is serial per process; parallelize across processes.

All run times quoted in this directory are those of this Julia program. They say nothing
about the speed of the authors' method or of their own code.

## Validation
- Bubble (the Type-1.2 example of the authors' demo notebook functions_demo.nb): exact
  {s, mm, s-4mm}, chi 3=1+1+1; diagnostics 1.1/2.1/2.2 PASS, 1.2 correctly FAILS
  flagging s (second-type, nu=0-degenerate).
- Sunrise: chi_top=4, chi 7=4+1+1+1, factors {s, mm, s-mm, s-9mm} (the paper's
  worked example).
- Reconstruction route: bubble and sunrise reproduce the exact route, and the sunrise
  its standard threshold/pseudo-threshold singular set (test_reconstruct.jl).
- The bubble and the sunrise are the smallest cases, chosen so that the output can be
  checked by hand. They are far simpler than the diagrams for which Chestnov, Crisanti
  and Giroux developed the method.
- Not implemented (documented frontier): SPQRDet (Faddeev–LeVerrier characteristic-polynomial
  route), 4d-restriction parametrisation, the SPQR/FiniteFlow Macaulay-matrix route.

## License

The parts of this package written for BootLoops are released under the MIT License (LICENSE in this directory), Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. The translated portions listed in NOTICE and PATCHES.md derive from the Mathematica code of Chestnov, Crisanti and Giroux (github.com/Giu989/Landau-s-Leviathans @ 47028b7) and of Crisanti, Lippstreu, McLeod and Polackova (github.com/Giu989/DiscKosky @ 5dca910). Those repositories carry no license file, and the Anthropic copyright and the MIT grant do not cover the translated portions. Documentation is licensed CC BY 4.0. See also LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).
