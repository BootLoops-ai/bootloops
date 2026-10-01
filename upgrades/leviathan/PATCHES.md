# PATCHES — Leviathan chi-drop engine

## What this is

A Julia implementation of the Euler-characteristic-drop method for Landau singularities
of Feynman integrals, from the papers of Vsevolod Chestnov, Giulio Crisanti and Mathieu
Giroux (arXiv:2606.29612) and Vsevolod Chestnov and Giulio Crisanti (arXiv:2511.14875).
No code was copied or translated: this is **not a patched tree** of the authors' code,
none of their code ships here, and none was translated. The authors' public Mathematica
code (Landau-s-Leviathans/codes_and_examples/landau_codes.m and
DiscKosky/Kernel/euler_characteristic.m on GitHub) was read as the reference for five
details the papers leave to the code: the Type-1.2 non-degeneracy test, the
"Diophantine" locus sampling (lowest-degree variable, 400 tries), the cut-edge regulated
ideal, normalization by the least common multiple of denominators, and dehomogenization
in the last kinematic variable (section "What follows the authors' code" below). The
engine was validated against the papers' worked examples and the singularity lists
published in the authors' repositories.

## Attribution

- Papers (the algorithm source, CC BY 4.0):
  - arXiv:2606.29612 — Vsevolod Chestnov, Giulio Crisanti, Mathieu Giroux, *Landau's
    Leviathans* (the chi-drop criterion, Eqs. 5/8/12, the Appendix-A diagnostics).
  - arXiv:2511.14875 — Vsevolod Chestnov, Giulio Crisanti, *SPQR* (the companion-matrix
    minimal-polynomial object the exact candidate route computes).
- The authors' code (read as the reference for the five details listed below and as the
  source of published validation values; the repositories carry no license file, and no
  code from them is copied or translated; we thank the authors for publishing them):
  - github.com/Giu989/Landau-s-Leviathans @ 47028b7 (Chestnov, Crisanti, Giroux;
    codes_and_examples/landau_codes.m, written in the Wolfram Language)
  - github.com/Giu989/SPQR @ dec9b27 (Chestnov, Crisanti; Kernel/*.m, a Mathematica
    package on FiniteFlow)
  - github.com/Giu989/DiscKosky @ 5dca910 (Giulio Crisanti, Luke Lippstreu, Andrew J.
    McLeod, Maria Polackova; Kernel/euler_characteristic.m, a Mathematica package)

## Why a new implementation instead of a patch

1. **License**: the authors' repositories carry no license file, so modified
   redistribution — the normal route for the engines in this directory — is not
   available. The papers' text is CC BY 4.0, so the algorithm specification is.
   Accordingly the engine is written from the papers; the authors' code serves as the
   reference for the five details listed below and as the source of published
   validation data.
2. **Stack**: the authors' SPQR requires Mathematica 13.1+ with the FiniteFlow MathLink;
   this implementation runs entirely on the open Julia + Oscar (Singular/msolve/FLINT)
   stack.

## What is different from the authors' code, and why

- **Exact candidate route** (src/elimination.jl): the per-variable elimination polynomial
  is obtained as the monic minimal polynomial in K[x]/I via ONE degrevlex Groebner basis
  plus incremental normal forms and a K-linear solve — the direct per-variable
  `eliminate()` over the rational function field hangs already on the sunrise top
  sector (>90 s vs seconds).
- **Reconstruction route** (src/reconstruct.jl, src/ratrec.jl): where the authors' SPQR uses
  FiniteFlow's Macaulay-matrix machinery, this engine samples the same minimal polynomial
  over (prime, kinematic point) grids and reconstructs coefficients by dense-ansatz
  nullspace + CRT + rational reconstruction — no FiniteFlow dependency.
- **Own DRL staircase counter** (src/staircase.jl): standard-monomial counting from
  leading exponents computed in-house (Oscar 1.7's `leading_exponent_vector` without an
  explicit ordering is wrong for f4 output — measured footgun, documented in README.md).
  Counting is an exact interval-sliced recursion over the staircase's threshold
  structure, so chi is exact regardless of box volume; test_staircase.jl (engine-free)
  cross-checks it against a brute-force box-walk oracle and huge-box closed forms.
- **Diagnostics** (src/diagnostics.jl): the four Appendix-A checks. The Type-1.2 test
  (`nodegeneracy`) follows `noDegeneracyQ` / `degenerateForKinAndVar` in the authors'
  landau_codes.m, which the paper names as the algorithm (Appendix A: "included in the
  repository as noDegeneracyQ") without spelling it out: for each (kinematic, variable)
  pair the kinematic is promoted to a ring variable, the variable is demoted to a
  parameter at two random values, the x-variables are eliminated, and the GCD of the two
  eliminants must be trivial. It runs on the nu=0 ideal <dG/dx_i, 1 - x0 G> of the
  paper's Eq. 12; a regulated-ideal variant false-PASSes, which is why the nu=0 ideal is
  the one used. Re-expressed in Oscar (eliminate, gcd, residues mod p < 2^29); no code
  copied.
- **Not implemented** (documented frontier, not silent gaps): SPQRDet (Faddeev–LeVerrier
  characteristic-polynomial route), the 4d-restriction parametrisation, and the
  SPQR/FiniteFlow Macaulay-matrix route itself.

## What follows the authors' code

The papers leave five details to the implementation; for these the authors' Mathematica
code was read and its behavior reproduced in new Julia/Oscar code (ours : theirs, line
numbers at the pinned commits):

1. Type-1.2 non-degeneracy test — src/diagnostics.jl:56-95 (`nodegeneracy`) :
   landau_codes.m:33-63 (`degenerateForKinAndVar`) and :66-116 (`noDegeneracyQ`)
   @ 47028b7 (Chestnov, Crisanti, Giroux).
2. "Diophantine" locus sampling (solve the candidate locus for its lowest-degree
   kinematic variable at random values of the others, up to 400 attempts) —
   src/chi.jl:45-72 (`point_on_locus`) : DiscKosky Kernel/euler_characteristic.m:4, 16,
   43-52, 118-127 (`countInSector`) @ 5dca910 (Crisanti, Lippstreu, McLeod, Polackova).
3. Cut-edge regulated ideal (a cut edge contributes the plain derivative dG/dx_i and is
   excluded from the regulators and the saturation product) — src/chi.jl:17-43
   (`chi_regulated`) and src/elimination.jl:27-33 : euler_characteristic.m:157-180
   (`prepareCountRegulatedMS`) @ 5dca910.
4. Normalization by the least common multiple of denominators (candidates = irreducible
   factors of the LCM of the denominators of the monic eliminant coefficients, numeric
   factors dropped) — src/elimination.jl:3-4, 90-91, 118-121 : landau_codes.m:238
   (`SPQRLandau`) @ 47028b7.
5. Dehomogenization in the last kinematic variable (kin[end] set to 1, factors
   rehomogenized by total degree, kin[end] re-added as a candidate) —
   src/reconstruct.jl:55, 109, 140 and src/ratrec.jl:112 (`rehomogenize`) :
   landau_codes.m:143-147, 240-247 ("AutoDehomogenise", `extraSing`) @ 47028b7.

## Test-vector provenance

The test expectations embedded in test_smoke.jl / test_sunrise.jl / test_leviathan.jl
are the papers' own worked examples (equal-mass bubble and sunrise; CC-BY 4.0). One
testset in test_reconstruct.jl checks the sunrise with equal internal masses (mm1) and
external invariant MM1 against its standard singular set (normal threshold MM1 = 9 mm1,
pseudo-threshold MM1 = mm1, and the vanishing loci mm1 = 0, MM1 = 0), independently
reproduced by both routes in the same test. No data files from the authors' repositories
are redistributed here.

## License

New Julia code (nothing copied or translated). MIT License, Copyright (c) 2026 Anthropic, PBC (see `LICENSE` in
this directory and the repository-root LICENSE and NOTICE). Created by Matthew D. Schwartz;
code written by Claude (Anthropic) under his supervision. Oscar (and the
Singular, msolve and GAP libraries underneath it) is a GPL-licensed dependency obtained
through the Julia package manager; it is not bundled, and using it at run time does not
change the license of this source.

## Validation

Validation: bubble and sunrise exact routes, reconstruction-vs-exact
cross-checks, and sector-chi / singularity-list agreement with the singularity lists
published in the authors' repositories.
