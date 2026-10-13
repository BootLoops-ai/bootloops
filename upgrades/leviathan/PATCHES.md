# PATCHES: Leviathan chi-drop engine

## What this is

A Julia port of the Euler-characteristic-drop method for Landau singularities of Feynman
integrals of Vsevolod Chestnov, Giulio Crisanti and Mathieu Giroux (arXiv:2606.29612) and
of the elimination object of SPQR by Vsevolod Chestnov and Giulio Crisanti
(arXiv:2511.14875). This is not a patched tree of the authors' code and no file of theirs
ships here. Parts of the Julia code are translations of their Mathematica routines, and
other parts follow conventions of their code; both are listed below, with our file and
line and theirs (section "What is translated from, or follows, the authors' code"). The
rest implements the equations of the papers with algorithms written for this package.
The engine was validated against the worked example of arXiv:2606.29612, the bubble
example of the authors' demo notebook and the singularity lists published in the
authors' repositories.

## Attribution

- Papers (the method; article text CC BY 4.0):
  - arXiv:2606.29612: Vsevolod Chestnov, Giulio Crisanti, Mathieu Giroux, *Landau's
    Leviathans* (the regulated ideal, Eq. 3; the chi-drop criterion, Eq. 5; the
    elimination ideal and its leading coefficient, Eqs. 8 and 9; the sector ideals and
    the sector sum, Eqs. 10 to 12; the Appendix-A diagnostics). Equation numbers are
    counted in the TeX source of version 1.
  - arXiv:2511.14875: Vsevolod Chestnov, Giulio Crisanti, *SPQR* (the companion-matrix
    minimal-polynomial object the exact candidate route computes).
- The authors' code (written in the Wolfram Language; this package was written with it
  open, and it is also the source of published validation values; the repositories
  carry no license file; we thank the authors for publishing them):
  - github.com/Giu989/Landau-s-Leviathans @ 47028b7 (Vsevolod Chestnov, Giulio Crisanti,
    Mathieu Giroux; codes_and_examples/landau_codes.m and functions_demo.nb)
  - github.com/Giu989/DiscKosky @ 5dca910 (Giulio Crisanti, Luke Lippstreu, Andrew J.
    McLeod, Maria Polackova; Kernel/euler_characteristic.m and
    Kernel/irreducible_monomials.m)
  - github.com/Giu989/SPQR @ dec9b27 (Vsevolod Chestnov, Giulio Crisanti;
    Kernel/various_functions.m; a Mathematica package on FiniteFlow)

## Why a Julia program instead of a patch

1. **License**: the authors' repositories carry no license file, so their files are not
   redistributed here. The article text of arXiv:2606.29612 is CC BY 4.0. The engine is
   written partly from the papers and partly from the authors' code, as listed below.
2. **Stack**: the authors' SPQR requires Mathematica 13.1+ with the FiniteFlow MathLink;
   this program runs entirely on the open Julia + Oscar (Singular/msolve/FLINT) stack.

## What is different from the authors' code, and why

- **Exact candidate route** (src/elimination.jl): the per-variable elimination polynomial
  is obtained as the monic minimal polynomial in K[x]/I via ONE degrevlex Groebner basis
  plus incremental normal forms and a K-linear solve. Oscar's direct per-variable
  `eliminate()` over the rational function field hangs already on the sunrise top
  sector (>90 s against seconds, both in this Julia program).
- **Reconstruction route** (src/reconstruct.jl, src/ratrec.jl): where the authors' SPQR uses
  FiniteFlow's Macaulay-matrix machinery, this engine samples the same minimal polynomial
  over (prime, kinematic point) grids and reconstructs coefficients by dense-ansatz
  nullspace + CRT + rational reconstruction, with no FiniteFlow dependency.
- **Staircase counter** (src/staircase.jl): the counting algorithm is an exact
  interval-sliced recursion over the staircase's threshold structure, written for this
  package (the authors' counter is a different algorithm), so chi is exact regardless of
  box volume; test_staircase.jl (engine-free) cross-checks it against a brute-force
  box-walk oracle and huge-box closed forms. Leading exponents are computed in this
  package (Oscar 1.7's `leading_exponent_vector` without an explicit ordering is wrong
  for f4 output; documented in README.md). The wrapper around the counter follows the
  authors' code (item 6 below).
- **Diagnostics** (src/diagnostics.jl): the four Appendix-A checks. The Type-1.2 test
  (`nodegeneracy`) is a translation of `noDegeneracyQ` / `degenerateForKinAndVar` in the
  authors' landau_codes.m, which the paper names as the algorithm (Appendix A: "included
  in the repository as noDegeneracyQ") without spelling it out: for each (kinematic,
  variable) pair the kinematic is promoted to a ring variable, the variable is demoted
  to a parameter at two random values, the x-variables are eliminated, and the GCD of the
  two eliminants must be trivial. It runs on the nu=0 ideal <dG/dx_i, 1 - x0 G> of the
  paper's Eq. 10, as the paper and the authors' code both specify. Differences from
  their routine: the two routines are fused into one, one fixed prime below 2^29 is
  used, and their options are absent.
- **Not implemented** (documented frontier, not silent gaps): SPQRDet (Faddeev–LeVerrier
  characteristic-polynomial route), the 4d-restriction parametrisation, and the
  SPQR/FiniteFlow Macaulay-matrix route itself.

## What is translated from, or follows, the authors' code

Eight points of contact. "Translation" means that our statements correspond one to one,
in the same order, to a block of the authors' code, and that the block contains at least
one choice that is not in the paper. "Follows" means that a convention, an interface or
a skeleton was taken from their code and written in different statements. Names differ
throughout and the authors' options are absent. Line numbers: ours at this revision,
theirs at the pinned commits.

1. **Type-1.2 non-degeneracy test (translation).** src/diagnostics.jl:60-99
   (`nodegeneracy`) : landau_codes.m:33-63 (`degenerateForKinAndVar`) and :66-116
   (`noDegeneracyQ`) @ 47028b7 (Chestnov, Crisanti, Giroux). Its ring helper follows the
   same routine: src/family.jl:62-74 (`spec_ring_midq`) : landau_codes.m:39-40, 73.
2. **"Diophantine" locus sampling (translation).** Solve the candidate locus for its
   lowest-degree kinematic variable at random values of the others, up to 400 attempts:
   src/chi.jl:50-77 (`point_on_locus`) : DiscKosky Kernel/euler_characteristic.m:96-131
   (`prepareKinematicSpecializationMS`, "Diophantine" branch; the 400-attempt loop is
   :118-127) @ 5dca910 (Crisanti, Lippstreu, McLeod, Polackova).
3. **Cut-propagator regulated ideal (translation).** A cut edge contributes the plain
   derivative dG/dx_i and is excluded from the regulators and the saturation product:
   src/chi.jl:27-48 (`chi_regulated`; the cut rule is at :33-46), src/elimination.jl:26-33 (`ff_ideal`) and
   src/reconstruct.jl:15-21 (`fp_ideal`) : euler_characteristic.m:157-180
   (`prepareCountRegulatedMS`) @ 5dca910. With no cut edges these blocks are Eq. 3 of
   the paper. The cut list as a third input (src/family.jl:10; leviathan.jl:13, 29-30)
   follows the third argument of `CountSectorsUnregulated` and `CountSectorsRegulated`
   (euler_characteristic.m:204, 310).
4. **Normalization by the least common multiple of denominators (follows).** Candidates
   are the irreducible factors of the LCM of the denominators of the monic eliminant
   coefficients, numeric factors dropped: src/elimination.jl:89-91, 116-124 and
   src/reconstruct.jl:107-113 : landau_codes.m:238 (`SPQRLandau`) @ 47028b7.
5. **Dehomogenization in the last kinematic variable (follows).** kin[end] is set to 1,
   factors are rehomogenized by total degree, and kin[end] is added back as a candidate:
   src/reconstruct.jl:55, 109, 140 and src/ratrec.jl:111-117 (`rehomogenize`) :
   landau_codes.m:143-147, 240-247 ("AutoDehomogenise", `extraSing`) @ 47028b7.
6. **Wrapper around the monomial count (follows).** The return contract (0 for the unit
   ideal, a sentinel for a positive-dimensional ideal, else the count), the two guards
   and the split into three helpers: src/staircase.jl:1-4, 21-29, 78-79, 89-93 :
   DiscKosky Kernel/irreducible_monomials.m:5-6, 267-271, 303-304 @ 5dca910 and SPQR
   Kernel/various_functions.m:45-53 @ dec9b27 (Chestnov, Crisanti).
7. **Retry over primes in the Type-2.2 check (follows).** The point search is repeated
   over the prime list until the locus has a point: src/diagnostics.jl:41-47 :
   euler_characteristic.m:43-52 (`countInSector`) @ 5dca910.
8. **Driver's verify switch and labels (follows).** One switch turns verification on or
   off; with the switch off every candidate is returned under the label used for
   verified ones; the labels genuine and spurious: src/landau.jl:35-45 and
   src/reconstruct.jl:142-152 : landau_codes.m:253-262 ("CheckCaveats") @ 47028b7.

Counted on release 1.0.1 (the code lines are the same at this revision; only comments
changed): of the 885 non-blank lines of src/, 97 are translations (40 from
landau_codes.m, 57 from DiscKosky), 92 follow the authors' code, 278 implement the
equations of the papers (a counterpart exists in the authors' code, but nothing that is
only in their code is reproduced), and 418 have no counterpart in the authors' code
(parser, counting algorithm, rational reconstruction, input and output).
Counted by whole functions, where a function that contains a translated block counts in
full, 415 of the 885 lines come from the authors' code.

## Test-vector provenance

The test expectations embedded in test_smoke.jl / test_sunrise.jl / test_leviathan.jl
are the sunrise worked example of arXiv:2606.29612 and the equal-mass bubble, the
Type-1.2 example of the authors' demo notebook (functions_demo.nb). The bubble
polynomial is written here with both kinematic signs opposite to the notebook's; the choice of example and the
expected outcome (the non-degeneracy test fails and flags s) are the notebook's. One
testset in test_reconstruct.jl checks the sunrise with equal internal masses (mm1) and
external invariant MM1, variable names as in the authors' data files and in SOFIA,
against its standard singular set (normal threshold MM1 = 9 mm1, pseudo-threshold
MM1 = mm1, and the vanishing loci mm1 = 0, MM1 = 0), which both routes must return in
the same test. No data files from the authors' repositories are redistributed here.

## License

Julia code, in part translated from the authors' Mathematica code (see above).

- The parts written for BootLoops: MIT License, Copyright (c) 2026 Anthropic, PBC (see
  `LICENSE` in this directory and the repository-root LICENSE and NOTICE). Created by
  Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision.
- The translated blocks of items 1 to 3 above, as listed in NOTICE: they derive from landau_codes.m in
  github.com/Giu989/Landau-s-Leviathans @ 47028b7 by Chestnov, Crisanti and Giroux and
  from Kernel/euler_characteristic.m in github.com/Giu989/DiscKosky @ 5dca910 by
  Crisanti, Lippstreu, McLeod and Polackova. Those repositories carry no license file.
  The Anthropic copyright line and the MIT grant do not cover these parts. `NOTICE` in
  this directory lists them by file and line.

Oscar (and the Singular, msolve and GAP libraries underneath it) is a GPL-licensed
dependency obtained through the Julia package manager; it is not bundled, and using it
at run time does not change the license of this source.

## Validation

Validation: bubble and sunrise exact routes, reconstruction-vs-exact
cross-checks, and sector-chi / singularity-list agreement with the singularity lists
published in the authors' repositories.
