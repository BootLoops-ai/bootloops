# Surd — exact hyperlogarithm integration with radical last-variable letters (three-fold Cheng-Wu integrals on a one-parameter slice)

NAME: a surd is the old word for an irreducible root such as √2 or ∛5. The integrator carries such square- and cube-root letters
exactly through the last integration, where standard hyperlogarithm programs stop: the roots of the quadratic and cubic last-variable
letters enter the result as exact algebraic generators over Q(i)(t), never as floating-point numbers.
(Formerly `collinear-hyperlog`.)

Tool page: https://bootloops.ai/tools/surd.html

KIND: package (`slice/scripts/` drivers, stage 3, assembly and structure tools; `symbolic/scripts/` the a-priori-alphabet fibration
engine; `gate/scripts/` cell builder, numeric reference evaluators, Arb one-fold integrator; `gate/cells/`, `gate/lr/`, `data/` package
data; `tests/`). Short `README.md`; long manual `MANUAL.md`; output file format `FORMAT.md`; data `DATA.md`.

## PURPOSE

Exact closed-form integration of three-fold Cheng-Wu (Feynman-parameter / energy-fraction simplex) integrals with a rational integrand
over a KNOWN denominator alphabet, on a ONE-parameter slice t of the kinematics: two Brown-linear fibration steps exact in
Q[x,t,j]/(j^2+1), then a function-level third step with radical last-variable letters (square roots, cube roots) carried exactly in
Q(i)(t)[rho]/(L), giving weight <= 3 hyperlogarithms with algebraic arguments and exact coefficients; plus structure tools on the result
(alphabet table, square-class Galois symbol test per weight and slot, P-adic slot test, local behavior at letter roots, Landau-family
attribution, support versus an external alphabet) and the independent numeric reference evaluators that check every step. Built for and
validated on the LO QCD / N=4 collinear four-point energy correlator on the dipole slice (seven 1 -> 4 splitting channels, 342 support
groups; the N=4 result reproduces the closed form of arXiv:2401.06463 on the slice to 1e-26).

## USE-WHEN

- A parametric integral whose first eliminations are Brown-linear but whose LAST variable carries algebraic letters, so hyperlogarithm
  engines of the HyperInt class (`tools/subtropica`/HyperFLINT, `tools/formglue`'s `form_hyper`) return nothing or silently drop the
  algebraic branches, and you want the exact function on a one-parameter family of kinematics.
- You know the stage alphabets a priori (from a linear-reducibility run such as `tools/subtropica`'s `lr_refine`, or a Landau census)
  and want an integrator that never rediscovers letters: every denominator is factored and looked up, never expanded over a common
  denominator.
- Exact symbol/alphabet questions on the integrated FUNCTION (which algebraic letters survive at which weight; is log P a letter; the
  local exponent bound at a letter's root) rather than on the integrand.
- A >= 30-digit certified numeric value of the same class of integrals (Arb residues + double-exponential quadrature, two routes) as a
  reference value for any closed form: `gate/scripts/oracle.py`, `slice/scripts/slice_oracle.py`, `slice_po.py`; or the last integration
  of a HyperFLINT one-fold object in ball arithmetic: `gate/scripts/onefold_arb.py`.

## NOT-FOR

More than ONE symbolic kinematic parameter (four-variable objects reach 1e6-1e7 monomials and do not finish; measured); homotopy /
monodromy of hyperlogs (no evaluation off the principal sheet with the -i0 rule, so pinch-sheet Landau coefficients are out of reach);
Galois symbol tests for CUBIC letters in `symtest.py` (square classes only); automatic choice of the elimination order (it is fixed per
first pair by the chart table of `slice_common.py`, following the reducibility run); integrands not linear in the first two eliminated
variables; more than three integrations.

## INVOKE

One support group end to end (stages 1-2 with resume-safe checkpoints, the stage 1-2 check, stage 3 with its numeric check):

    export SURD_WORK=<scratch dir>          # every write goes here (default ./surd_work; the package tree is refused)
    cd slice/scripts
    python3 slice_run.py    --channel q_qbpqpgq --sigma 1234 --pair 34
    python3 slice_gate12.py --tag q_qbpqpgq_s1234_p34_g2
    python3 stage3.py       --tag q_qbpqpgq_s1234_p34_g2 --t 1/2,4/5

Many groups: `bash slice/scripts/s13_worker.sh LIST NAME` (one `channel sigma pair` per line; one worker per core on disjoint lists).
Classes and assembly: `slice_classes.py [--channels C]`, `assemble.py --channels C|all`. Evaluation of an assembled file:
`evaluate.py --file <label|path> --t p/q [--digits N]`. Structure tools: `letters.py --files L`, `symtest.py --files L [--data DIR]
[--tstar e/7 --tstar2 pi/9]`, `valtest.py --label L`, `tstar_local.py --label L`, `near_root.py`, `landau_family.py`, `n4_alphabet.py`.
Numeric reference evaluators: `slice_oracle.py --t LIST [--channels ...] [--nproc k]`, `file_gate.py [--receipt]`,
`gate/scripts/oracle.py --channel q_jet --shape 'u=3/7+2/5j;v=-1/3+4/9j' --digits 30`, `gate/scripts/onefold_arb.py --piece PIECE.json`.
Library: `import slice_prov, assemble; assemble.evaluate_file(path, Fraction(1,2), dps=60)` (with `slice/scripts` on `sys.path`).
Smoke test: `bash tests/smoke_test.sh [TARGET]`; self-test suite: `python3 tests/run_tests.py` (details in MANUAL.md section 3).

## INPUTS

Exact cell integrands (`gate/scripts/oracle_cells.build_cells`; included as `gate/cells/*.pkl` for the channel of the smoke test, built
from splitting-function term lists under `$SURD_TERMS` otherwise), the slice map d_ab(t) (`slice_common.py`, polynomial in t), the
elimination order per support group (chart table in `slice_common.py`) and the stage alphabets (`gate/lr/out_groups_symbolic.json`);
optionally a user-supplied `ckpt/cmsyz_cubics_slice.json` under the work directory (external cubic letters restricted to the slice,
`{"cubics": [{"sigmas": [...], "factor": "<polynomial in X, t>"}]}`) for the comparison in `letters.py`, skipped when absent. A new
problem supplies atoms, numerator, powers, external t-factor, an order that is Brown-linear in its first two variables, and its stage
alphabets (MANUAL.md section 5). Environment: `SURD_WORK` (output root), `SURD_DATA` (data root for the structure tools and tests;
default WORK), `SURD_TERMS`, `SURD_LR`, `SURD_N4SLICE`, `SURD_FIBRATE`, `SURD_E4C_ENGINE`, `SURD_VMEM_KB` (all with defaults or a clear
error when a tool needs them; table in MANUAL.md section 9).

## OUTPUTS

Under the output root: per group `ckpt/<tag>/stage{1,2,3}.pkl` and the run records `SLICE_S12_/G12_/S3_<tag>.json` (results plus a
producer block: script sha, launch line, input and module shas, host, UTC stamp; write-once, superseded copies kept); assembled
`out/*.json.gz` in the format of `FORMAT.md` (value = sum coef(t) * prod rho^e * prod const * prod Z(w); coefficients (A+iB)/D in
Q(i)(t); root-slot symmetric); structure records `SLICE_S5_*.json`; reference-value records `SLICE_S0_ORACLE_t*.json`.

## REQUIREMENTS

python3 (3.10 or newer) with python-flint >= 0.6 (developed and measured with 0.8.0), sympy (>= 1.12) and mpmath (>= 1.3); no compiled
component of its own. The hyperlog evaluators import `zip_num` from the sibling package `tools/subtropica` (`scripts/fibrate/zip_num.py`,
pure Python; found automatically in the repository layout, or point `SURD_FIBRATE` at its directory). The HyperFLINT engine itself is
NOT needed except for `hf_piece.py --run` and the piece files of `onefold_arb.py`; numpy only for the optional `--e4c` external-engine
cross-check. bash and GNU coreutils for the shell drivers; single thread throughout.

## CHECKS

Never trust a group without its stage 1-2 check (`slice_gate12.py`: symbolic one-fold object vs the independent residue-route fiber
function, >= 30 digits at >= 2 real x per t) AND its stage-3 check (`stage3.py`: closed form vs the two-chart Arb piece evaluator
`slice_po.py`, >= 30 digits at two t); an assembled file gets a file check (`file_gate.py`) against the full numeric reference values of
`slice_oracle.py` at >= 3 rational t including both sides of every letter root before any structure claim; a new integrand class needs
one end-to-end control against an independently known closed form (the N=4 pattern: 5/24 x the published function to 1e-26); symbol
verdicts only on COMPLETE files and at two point pairs with root-continuation labeling; check records carry `min_agree_digits` and a
boolean computed from the rows, never a typed status. Measured: stage 1-2 checks 33.9-51 digits and stage-3 checks 35.2-48 digits over
all 342 groups; file checks 32.8-40.6 digits over ten files.

Self-test suite: `python3 tests/run_tests.py` (hermetic, a few seconds: hyperlog evaluator vs closed forms and quadrature, the -i0 path
prescription vs an independent contour integral, exact cross-chart identity of the included cells, the N=4 alphabet data; the five data
tests are named and SKIP unless `--data`/`SURD_DATA` points at the assembled files, which are not part of the package); `--only s12`
adds stages 1-2 and the stage 1-2 check of one group (~50 s). `bash tests/smoke_test.sh [TARGET]`: one group end to end + assembly +
file evaluation against two frozen 30-digit values (stage 1-2 check 38.4 digits, stage-3 check 46.2 digits; 3-5 min, 0.2 GB).
Integrity: `sha256sum -c MANIFEST.sha256` from the package directory checks every included script and data file (the manifest is
regenerated by
`find slice symbolic gate data tests -type f \( -name '*.py' -o -name '*.sh' -o -name '*.pkl' -o -name '*.json' \) | LC_ALL=C sort | xargs sha256sum > MANIFEST.sha256`).

## FOOTGUNS

- Outputs never go into the package tree: set `SURD_WORK` (default `./surd_work`); a root inside the package is refused at the first
  write. Cells for channels other than the included one are built on first use from term lists you supply (`SURD_TERMS`, format in
  DATA.md) and cached under the output root.
- Representation-singular rational t (a letter coefficient or point denominator vanishes while F is analytic; t = 3/5 and 5/8 in the E4C
  files): `evaluate_file` raises; use `evaluate.py` (symmetric limit) or another t. Digits ~ dps minus termwise cancellation depth (10
  generic, up to ~30 near colliding roots); |Im F| is the monitor.
- Euclidean slices make fiber letters Gaussian: everything runs in Q[x,t,j]/(j^2+1) (`fibr_gi.py`, `alphabet_gi.py`); a real-rational
  shortcut silently loses terms.
- "Absent from the weight-w SYMBOL" is not "absent at weight w" (i pi and zeta terms are symbol-invisible); the local Laurent capacity m
  bounds exponents (no worse than (t - t*)^(-m)); it does not give the pinch-sheet coefficient.
- LLL relation searches in real quadratic fields with large discriminant can miss a relation: `symtest.py` verifies every relation found
  at t* again at an independent t**; treat an unverified relation as absent.
- Splitting numerators to save memory destroys residue cancellations (a 40k-term chunk -> 11M monomials where the whole group has 195k):
  split by SUPPORT group (`hf_piece.first_pair`), not by term blocks.
- In the strict engine (`symbolic/scripts/fibr.py` on `alphabet.Alphabet`) an off-alphabet denominator raises `OffAlphabet` by design:
  extend the stage-alphabet file from a reducibility run; do not bypass. The slice drivers run the auto-extending Gaussian alphabet
  (`alphabet_gi.GIAlphabet`) and list every letter it added in the run record; compare that list with your reducibility run.
- A fixed-precision real-axis quadrature of a HyperFLINT one-fold object returns garbage (termwise magnitudes exceed the sum by 100+
  digits and the representation is termwise singular on the axis): `onefold_arb.py` rotates the contour and works in balls for that reason.

## CREDIT

The fibration method (linear reducibility, integration in a fibration basis of hyperlogarithms) is F. Brown, "The massless higher-loop
two-point function", Commun. Math. Phys. 287 (2009) 925, arXiv:0804.1660; its algorithmic form and the compatibility-graph reduction
follow E. Panzer, "Algorithms for the symbolic integration of hyperlogarithms with applications to Feynman integrals" (HyperInt),
Comput. Phys. Commun. 188 (2015) 148, arXiv:1403.3385. The projective (Cheng-Wu) freedom to set one Feynman parameter to 1 is the
Cheng-Wu theorem (H. Cheng and T. T. Wu, Expanding Protons, MIT Press, 1987). The QCD integrands are the tree-level quadruple-collinear splitting functions of V. Del Duca, C. Duhr, R. Haindl,
A. Lazopoulos and M. Michel, JHEP 02 (2020) 189 (arXiv:1912.06425) and JHEP 10 (2020) 093 (arXiv:2007.05345) [DDHLM1, DDHLM2], and the
shipped cell data derived from them (gate/cells/) carry their CC BY 4.0 attribution (DATA.md); the N=4 reference function
and its alphabet are from D. Chicherin, I. Moult, E. Sokatchev, K. Yan and Y. Zhu, "The Collinear Limit of the Four-Point Energy
Correlator in N = 4 Super Yang-Mills Theory", arXiv:2401.06463. Exact arithmetic by FLINT through python-flint (Arb ball arithmetic
included); polynomial algebra by SymPy; multiprecision evaluation by mpmath. The ZIP-word conventions and the `zip_num` evaluator are
those of the sibling package `tools/subtropica`.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his
supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components
keep their own licenses (THIRD_PARTY.md).
