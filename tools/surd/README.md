# Surd

Surd integrates three-fold Cheng-Wu (Feynman-parameter simplex) integrals exactly on a one-parameter slice t of the kinematics when the
last integration variable carries radical letters. Two Brown-linear fibration steps run exactly in Q[x,t,j]/(j^2+1); the function-level
third step carries the square- and cube-root letters of the last variable exactly, as algebraic generators over Q(i)(t), and returns
weight <= 3 hyperlogarithms with algebraic arguments and exact coefficients. Structure tools (letter table, Galois symbol test, P-adic slot
test, local behavior at letter roots, Landau-family attribution) and certified Arb reference evaluators come with it. The worked problem
is the LO QCD / N=4 collinear four-point energy correlator on the dipole slice.

The name: a surd is the old word for an irreducible root such as √2 or ∛5. The integrator carries such square- and cube-root letters
exactly through the last integration, the step at which linear reducibility fails and rational-alphabet hyperlogarithm programs
(HyperInt [HyperInt], HyperFLINT) hand over; for other treatments of square-root letters see RationalizeRoots by Besier, Wasser and
Weinzierl [RR] and Heller, von Manteuffel and Schabinger [HvMS] (keys in REFERENCES.md at the repository root).

Requirements: python3 >= 3.10 with python-flint (>= 0.6), sympy, mpmath; the sibling package `tools/subtropica` for one pure-Python
evaluator module (found automatically in the repository layout). Nothing to build.

Smoke test (one support group end to end, 3-5 min, single thread, writes only under the target directory):

    bash tests/smoke_test.sh /path/to/fresh/scratch/dir          # ends with "[smoke] RESULT: PASS"

Self-test suite: `python3 tests/run_tests.py` (a few seconds; data tests skip by name unless `--data DIR` is given).

Start with `GUIDE.md` (purpose, when to use it, invocation, checks, pitfalls), then `MANUAL.md` (installation, tests, a worked example,
setting up a new problem, method, tools, environment variables, formats, validation, limits); `FORMAT.md` specifies the output files and
`DATA.md` the package data. All outputs go under `$SURD_WORK` (default `./surd_work`), never into the package tree.

Tool page: https://bootloops.ai/tools/surd.html

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his
supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components
keep their own licenses (THIRD_PARTY.md).
