# Surd — manual

Surd is an a-priori-alphabet fibration integrator for three-fold Cheng-Wu (Feynman-parameter simplex) integrals whose integrand is a
rational function with a known denominator alphabet, on a one-parameter slice of the kinematics. The name: a surd is the old word for an
irreducible root such as √2 or ∛5, and the integrator carries such square- and cube-root letters of the last integration variable exactly
through the last integration, where standard hyperlogarithm programs stop. `README.md` is the one-page entry; `GUIDE.md` is the short
guide (purpose, when to use it, how to invoke, checks, hazards); `FORMAT.md` specifies the assembled closed-form files; `DATA.md`
describes the package data. This manual gives installation, the tests, a worked example, how to set up a new problem, the mathematics
and algorithms, the tools, the file and record formats, the measured costs, the validation numbers and the limits.

Contents: 1 What Surd does; 2 Installation and requirements; 3 Smoke test and self-test suite; 4 Worked example; 5 Setting up a new
problem; 6 Method in more detail; 7 Structure tools; 8 Numeric reference evaluators; 9 Environment variables; 10 Outputs, record
format, file format and package data; 11 Branches, base point, precision; 12 Resource profile; 13 Validation; 14 Limits.

## 1. What Surd does

The tool computes, exactly,

    I(t) = int_{(0,oo)^3} d^3x  N(x_a, x_b, x_c; t) / prod_k A_k(x; t)^{M_k}

where every kinematic invariant is a polynomial in one parameter t, the atoms A_k belong to a known alphabet, and the integrand is
Brown-linear (linearly reducible) in the first two variables to be eliminated. The result is a sum of weight <= 3 hyperlogarithms with
algebraic arguments and exact coefficients in Q(i)(t) extended by the roots of the last-variable letters:

* stages 1-2: two Brown-linear fibration steps, exact in Q[x, t, j]/(j^2 + 1). Gaussian letters arise on Euclidean slices, so the
  fibration arithmetic of `symbolic/scripts/fibr.py` is transported verbatim to the Gaussian ring (`slice/scripts/fibr_gi.py`,
  `alphabet_gi.py`). Output: the one-fold object, a term store of rational prefactors times ZIP words with x-dependent letters, checked
  against an independent residue-route fiber function (`slice_gate12.py`);
* stage 3 (function level): fibration of every ZIP word into G_0(points; x) over the roots of the last-variable alphabet (rational points
  in Q(i)(t), roots of quadratic and cubic X-letters) with base-point constants at x0 = 11/2, then exact integration over (0, oo) by
  algebraic partial fractions in Q(i)(t)[y]/(L) and pole/polynomial-part recursions, coefficients in the tensor ring of root generators
  (`k3field.py`, `fib3.py`, `int3.py`); checked per group against a two-chart Arb piece evaluator (`slice_po.py`) and, after assembly,
  against the full numeric reference values (`file_gate.py` vs `slice_oracle.py`);
* assembly of groups into channels and jets with exact collection of identical hyperlog monomials (`assemble.py`) and a numeric evaluator
  of the resulting files (`assemble.evaluate_file`; guarded command line `evaluate.py`);
* structure tools on the assembled files (section 7): `letters.py` (alphabet and discriminant table), `symtest.py` (square-class Galois
  content of the weight-w symbol, per slot), `valtest.py` (P-adic slot functionals: is log P a letter), `tstar_local.py` and `near_root.py`
  (local behavior at letter roots), `landau_family.py` (which denominator atoms generate a letter), `n4_alphabet.py` (support versus an
  external alphabet).

Everything is exact (flint `fmpq`/`fmpq_mpoly`, `Fraction`); no PSLQ, fit or tolerance enters the construction.

The worked problem included with the package is the LO QCD and N=4 SYM collinear four-point energy correlator (E4C) on the dipole slice:
seven 1 -> 4 splitting channels plus the N=4 control, 342 support groups, and the N=4 result reproduces the known closed form of
arXiv:2401.06463 on the slice (section 13). The stage labels used in file names are S12 (stages 1-2), G12 (the stage 1-2 check), S3
(stage 3 and its check), S0 (the full numeric reference values), S5 (structure tools).

## 2. Installation and requirements

There is nothing to build. Requirements: python3 (3.10 or newer) with python-flint >= 0.6 (developed and measured with 0.8.0), sympy
(>= 1.12) and mpmath (>= 1.3), for example `pip install python-flint sympy mpmath`; bash and GNU coreutils for the shell drivers. The
hyperlog evaluators import `zip_num` from the sibling package `tools/subtropica` (`scripts/fibrate/zip_num.py`, pure Python, mpmath
only): in the repository layout it is found automatically; otherwise set `SURD_FIBRATE` to the directory that holds `zip_num.py`. The
HyperFLINT engine of `tools/subtropica` is NOT needed except for `hf_piece.py --run` and for producing the piece files that
`onefold_arb.py` consumes (then set `SUBTROPICA_HF_BIN` and `SUBTROPICA_MZV_DATA`); numpy only for the optional `--e4c` external-engine
cross-check. Everything runs single thread.

Integrity of the package files: `sha256sum -c MANIFEST.sha256` from the package directory checks every included script and data file.

All outputs go under one output root, `$SURD_WORK` (default `./surd_work` under the current directory); a root inside the package tree is
refused at the first write. Set it before running anything:

    export SURD_WORK=/path/to/scratch/surd_work

## 3. Smoke test and self-test suite

Smoke test (one support group end to end plus assembly and file evaluation, in a fresh scratch directory, against two frozen 30-digit
values; single thread, 3-5 min, ~0.2 GB):

    bash tests/smoke_test.sh [TARGET]          # TARGET = a fresh directory (default ${TMPDIR:-/tmp}/surd_smoke.<stamp>)

It runs `slice_run.py`, `slice_gate12.py`, `stage3.py`, `slice_classes.py`, `assemble.py` and `evaluate.py` for the group
`q_qbpqpgq_s1234_p34_g2` from the included cell pickles, prints the stage 1-2 check (38.4 digits), the stage-3 check (46.2 digits), the
two frozen-value matches and the assembled-file evaluation, and ends with `[smoke] RESULT: PASS` and exit code 0. It sets `SURD_WORK`
and `SURD_DATA` to the target itself, so nothing else needs to be set.

Self-test suite:

    python3 tests/run_tests.py                       # hermetic tests, a few seconds, writes nothing
    python3 tests/run_tests.py --only s12            # + stages 1-2 and the stage 1-2 check of one group (~50 s; writes under $SURD_WORK or a temporary directory)
    python3 tests/run_tests.py --data DIR            # + the five data tests (DIR holds out/*.json.gz and SLICE_S5_LETTERS.json; ~56 min, 2.1 GiB)
    python3 tests/run_tests.py --data DIR --only b,d2   # the two quick data tests (~8 min)

The hermetic tests are (h) the hyperlog evaluator against closed forms and quadrature, (p) the -i0 prescription against an independent
contour integral, (g) the exact cross-chart identity of the included cells and the frozen census of the smoke-test group, (n) the N=4
alphabet data file. The five data tests (a, b, c, d2, e; the section 13 numbers as frozen expectations) need the assembled files, which
are not part of the package (about 0.5 GB); without `--data DIR` or `$SURD_DATA` they are SKIPPED BY NAME, not failed. The last line is
`RESULT: PASS (n ran, m skipped: names)` with exit code 0 when every test that ran passed.

## 4. Worked example: one support group end to end

The unit of work is a support group: the cells of one channel in one chart that share a first admissible vertex pair, at one label
assignment sigma (section 6). Its tag is `<channel>_s<sigma>_p<pair>_g<gauge>`. Every step is single thread and resume-safe (a rerun
continues from the last checkpoint or returns at once if its record exists).

    export SURD_WORK=<scratch dir>
    cd slice/scripts
    python3 slice_run.py    --channel q_qbpqpgq --sigma 1234 --pair 34      # stages 1-2: ckpt/<tag>/stage{1,2}.pkl, SLICE_S12_<tag>.json   (~50 s)
    python3 slice_gate12.py --tag q_qbpqpgq_s1234_p34_g2                    # stage 1-2 check: SLICE_G12_<tag>.json                        (~3 s)
    python3 stage3.py       --tag q_qbpqpgq_s1234_p34_g2 --t 1/2,4/5        # stage 3 + piece check: stage3.pkl, SLICE_S3_<tag>.json       (~90 s)
                                                                            #   --reuse re-evaluates and re-checks from the pickle; --no-po skips the piece evaluator

What to look at: `SLICE_G12_<tag>.json` and `SLICE_S3_<tag>.json` carry `min_agree_digits` and the boolean `PASS(>=30)`; the stage-3
record also lists the closed-form value divided by the external factor at each t (`gate_rows`), the census of terms, letters, root
generators and weights, wall time and peak RSS. The log of each step is under `$SURD_WORK/logs/`.

Then the sigma classes of the channel, the assembly of every finished group of the channel into one file, and the evaluation of that
file at a rational t:

    python3 slice_classes.py --channels q_qbpqpgq                           # SLICE_CLASSES.json (exact relabeling symmetries, class multiplicities)
    python3 assemble.py      --channels q_qbpqpgq                           # out/E4C_LO_QCD_dipole_slice_q_qbpqpgq.json.gz, the partial quark-jet file (+ .txt summary), SLICE_ASSEMBLY.json
    python3 evaluate.py      --file $SURD_WORK/out/E4C_LO_QCD_dipole_slice_q_qbpqpgq.json.gz --t 1/2 --digits 30

With one group done the channel file is marked incomplete (`complete: false`, `missing_groups` listed) and its value is that group's
contribution times its class multiplicity and flavor weight; the smoke test checks exactly this number. A complete channel needs all of
its groups (48 for q_qbpqpgq: 12 sigma classes x 4 first pairs).

Many groups: write a list file with one `channel sigma pair` per line and run one worker per core on disjoint lists,

    bash slice/scripts/s13_worker.sh LIST NAME [t-list]                     # skips groups whose SLICE_S3 record exists; logs under $SURD_WORK/logs

then `assemble.py --channels LIST|all` per jet (records merge across invocations), the full numeric reference values
`slice_oracle.py --t 1/2,2/3,4/5 --channels LIST --nproc k` (checkpoint per record; resume = rerun; `--assemble-only`), and the file
check `file_gate.py [--labels ...] [--s ...]` followed by `file_gate.py --receipt` for the summary record. Labels accepted by the tools
are `quark`, `gluon`, `n4`, a channel name, or a file name/path under `$SURD_DATA/out`.

From Python (with `slice/scripts` on `sys.path`):

    import slice_prov, assemble
    from fractions import Fraction
    F = assemble.evaluate_file(path, Fraction(1, 2), dps=60)      # mpmath complex; the exact value is real, |Im F| is the error monitor

## 5. Setting up a new problem

A problem is specified by five things: the atoms (polynomials in x_a, x_b, x_c and t), the numerator, the atom powers, the external
t-factor pulled out of the integrand, and an elimination order (first pair, chart) for which the first two eliminations are
Brown-linear, together with the stage alphabets of that order. In the code these enter through a few small interfaces:

* **The group integrand.** `slice_common.build_group(channel, gauge, pair, sigma)` returns the dict that `slice_run.py`, `stage3.py`
  and the structure tools consume: `N` (the numerator, an `fmpq_mpoly` in the context `(x_a, x_b, x_c, t)`), `atoms` (`{name: fmpq_mpoly}`
  in the same context), `M` (`{name: integer power}`), `xvars` (the three chart variables), `external_d_powers` and `ext_factor_poly`
  (an `fmpq_poly` in t; the true integrand is `N / prod atoms^M / ext_factor(t)`), `ctx`, and labels (`channel`, `gauge`, `pair`,
  `sigma`, `n_cells`, `dsigma_polys`). For the E4C it is built from the exact cell integrands (`oracle_cells.build_cells(channel, gauge)`;
  cell format in `DATA.md`) composed with the slice map. A new integrand inside the same atom family {x_j, Om = sum x, ell_S = partial
  sums, s_S = sum_{a<b in S} d_ab x_a x_b} needs only a new term list or cell pickle (`$SURD_TERMS`, `DATA.md`) and, if the kinematics
  differ, a new slice map; an integrand outside that family needs a replacement for `build_group` with the same return contract.
* **The slice.** `slice_common.d_of_t` / `d_polys_t` / `slice_ctx` give the six d_ab as polynomials in t (here d12 = t^2, d34 = 1,
  d13 = d24 = t^2/4 + 3t/10 + 1/4, d14 = d23 = t^2/4 - 3t/10 + 1/4: detector positions w1 = (t/2) e^{i phi}, w2 = -w1, w3 = -1/2,
  w4 = 1/2 with cos phi = 3/5, so that all d_ab are polynomial in t). The invariants must be POLYNOMIAL in t; choose the slice so.
* **The elimination order.** `slice_common.CHART` fixes the chart variable per first pair and `order_for(pair, gauge)` the order
  (x_k, x_j of the pair, then the last variable); `hf_piece.first_pair` / `group_cells` assign cells to groups. The first two variables
  must be Brown-linear for the group; take the order from a linear-reducibility run.
* **The stage alphabets.** `gate/lr/out_groups_symbolic.json` (format in `DATA.md`; produced by the `lr_refine` member of
  `tools/subtropica` on the cell integrands) lists, per (channel, chart, group, order), the letters alive after each integration. The
  strict engine (`symbolic/scripts/alphabet.Alphabet.load_lr` with `fibr.Engine`) looks every denominator factor up there and raises
  `OffAlphabet` for an unknown one; the slice drivers use the auto-extending Gaussian alphabet `alphabet_gi.GIAlphabet`, which splits
  factors over Q(i), adds what it meets and writes the full letter list into `SLICE_S12_<tag>.json` (`alphabet`,
  `nonlinear_last_variable_letters_used`), to be compared with the reducibility run. Point `SURD_LR` at your file.
* **The numeric reference.** The checks need the same integrand on an independent numeric route: `onefold_arb.fibre_residue` (stage 1-2
  check) and `oracle_core` (piece evaluator and full reference) integrate the cells numerically, so a problem given as cells gets its
  checks for free; a problem given through a custom `build_group` must also supply its cells (or its own reference values in the piece
  file format of section 8) before any result is trusted.
* **Group lists** (`channel sigma pair` per line) for `s13_worker.sh`, and the sigma classes from `slice_classes.py` (`SLICE_CLASSES.json`),
  which verifies every relabeling symmetry exactly on the cells before using it.

A new integrand class also needs one end-to-end control against an independently known closed form before its results are used (for
the E4C this is the N=4 function, section 13).

## 6. Method in more detail

**Cells.** `gate/scripts/oracle_cells.py` maps a 1 -> 4 splitting-function term list to Cheng-Wu cells in the chart x_g = 1: a cell is a
denominator multiset over the atoms {x_j, Om = sum x, ell_S = sum_{j in S} x_j, s_S = sum_{a<b in S} d_ab x_a x_b} with a numerator that is a
polynomial in the three chart variables and the six pair distances d_ab. The kinematics stays symbolic in the d_ab, so a permutation
sigma of the final-state labels and a point of the slice are both substitutions d_ab -> d_{sigma(a) sigma(b)}(t).

**Support groups.** The cells of a channel and chart are grouped by their first admissible vertex pair {k, j}: the lexicographically first
pair contained in no three-particle quadric s_S and no partial-sum plane ell_S of the cell's support (`hf_piece.first_pair`). A group is
integrated in the order x_k, x_j, then the last variable x_l, in the chart fixed per pair (`slice_common.CHART`). Splitting by support group
rather than by term blocks is essential: the residue cancellations live inside a group (a 40k-term block split reached 11M monomials where
the whole group has 195k).

**Sigma classes.** The 24 label assignments are reduced modulo the slice automorphism (12)(34) of the labeled distance matrix (identical
integrands) and the channel's exact relabeling symmetries tau (I_{sigma tau} = I_sigma), each verified exactly on the cells at random
rational points before use (`slice_classes.py`); 12 or 6 classes per channel remain, 342 (channel, class, first-pair) group objects in all.

**Stages 1-2.** For a group, the integrand N / prod A^M is lifted to Q[x_a, x_b, x_c, t, j]/(j^2 + 1) and integrated in x_k then x_j by the
fibration-basis engine: partial fractions in the active variable with every new denominator (resultants a b' - a' b, leading coefficients,
constant terms) FACTORED and each irreducible factor LOOKED UP in the stage alphabet (`symbolic/scripts/alphabet.py`; the strict table is
seeded from `gate/lr/out_groups_symbolic.json` and refuses an unknown factor with `OffAlphabet`, the slice drivers run its auto-extending
Gaussian version `alphabet_gi.GIAlphabet`, which splits factors over Q(i) and records every letter it adds); nothing is ever expanded
over a common denominator. Primitives are ZIP words (shuffle-regularized G(...; oo) from 0, the ZeroInfPeriod convention:
log-divergences at 0 and oo dropped). Checkpoints `ckpt/<tag>/stage1.pkl`, `stage2.pkl`; a rerun resumes.

**Stage 1-2 check.** The stage-2 one-fold object F(x; t), specialized to rational t and j -> i and evaluated at real x in {7/3, 11/2, 3/5,
37/4} for t in {1/2, 4/5} (hyperlogs by `hpath.py` with the -i0 rule), is compared with the residue-route fiber function of
`gate/scripts/onefold_arb.fibre_residue` (exact cells with numeric d_ab(t), inner integral by residues in acb, one-dimensional exp-sinh in
Arb; no code shared with the symbolic engine). Pass = >= 30 digits at >= 2 points per t (`slice_gate12.py`).

**Stage 3.** A. Every product of ZIP words of the one-fold object is fibrated into G_0(points; x) with constants at the base point
x0 = 11/2 (`fib3.Fibrator`); the points are 0, rational points of Q(i)(t), and the roots rho_{L,s} of the non-linear last-variable letters
L(X; t). B. Every (rational function of x) x (basis word) is integrated exactly over (0, oo) with algebraic partial fractions in
Q(i)(t)[y]/(L) and the pole / polynomial-part recursions (`int3.Integrator`); coefficients live in the tensor ring over Q(i)(t) generated by
the (letter, slot) root symbols (`k3field.TRing`), so no root is ever approximated. Output `ckpt/<tag>/stage3.pkl`: terms
coef(t, roots) x CONST x prod Z(words), with CONST = ZIP words at x0 and G_0(v; x0), and Z(u) = reg_oo G_0(u; oo) over points.

**Stage-3 check.** The closed form, evaluated at t = 1/2 and 4/5 (roots by `mpmath.polyroots`, hyperlogs by `hpath.py`) and divided by
the external factor prod d^ext(t), is compared with the piece evaluator of the group (`slice_po.py`: `oracle_core.py` integrates the
group's cells numerically by inner residues in Arb and a two-dimensional exp-sinh rule, in two charts/routes A = (x_4 = 1; outer x_1, inner
x_2, analytic x_3) and B = (x_1 = 1; x_2, x_3, x_4); A-vs-B agreement ~1e-42). Pass = >= 30 digits at both t.

**Assembly.** `assemble.py` maps every group term store to canonical keys (letter polynomials as strings in (X, t, j), points as (polynomial,
slot) or exact Q(i)(t) values, hyperlog words over points, constants) and collects identical hyperlog monomials with exact Q(i)(t)
coefficients: F_ch(t) = w_flavor(n_f = 5) S_ch sum_classes mult sum_groups I_group(t). The file format is in `FORMAT.md`.

**File check and the full numeric reference (S0).** `slice_oracle.py` computes, independently of everything above, the channel and jet
values at rational t (`gate/scripts/oracle.py` machinery: exact rational-function assembly per (channel, chart, sigma), inner residues
in Arb, two-dimensional exp-sinh, two routes; certified digits = -log10 |A - B|/|A|). `file_gate.py` evaluates the assembled files at
every t with a reference record (t > 1 through F(T) = s^6 F(s), s = 1/T) and reports the agreement.

## 7. Structure tools

All read assembled files `out/*.json.gz` (from `$SURD_WORK`, else `$SURD_DATA`) and write a structure record `SLICE_S5_*.json` under
`$SURD_WORK`. Labels: `quark`, `gluon`, `n4`, a channel name, or a file name/path.

| tool | question | command | record |
|---|---|---|---|
| `letters.py` | the letter table: algebraic X-letters per function with discriminants factored over Q[t] and matched against the twelve rational tangency letters of the dipole slice (`letters.TABLE1`); cubic letters (discriminant class, field isomorphism vs an external cubic if `ckpt/cmsyz_cubics_slice.json` is present); the rational t-letters; per quadratic letter and weight the net contribution S_L^(w)(t) of the terms carrying its roots | `letters.py --files quark,gluon,n4 [--t 1/2,4/5] [--noS]` | `SLICE_S5_LETTERS.json` |
| `symtest.py` | does the square class of Delta(t) survive in the weight-w symbol, per slot? (Galois functionals at a transcendental t* with LLL relations verified at t**; square classes only) | `symtest.py --files L [--data DIR] [--tstar e/7 --tstar2 pi/9] [--no-continuation]` (free-form points: names, `sqrt(3)/5`, `2/(3+sqrt5)`, decimals, p/q) | `SLICE_S5_SYMTEST*.json` |
| `valtest.py` | is log P(t) a letter in slot p of the weight-w symbol, for the rational tangency letters P? (P-adic slot functionals) | `valtest.py --label L [--weights 1,2,3]` | `SLICE_S5_VALTEST_<label>.json` |
| `tstar_local.py` | function-level behavior at the square-root letter classes on the principal sheet: odd part, boundedness, termwise pole order at t* -+ h | `tstar_local.py --label L [--classes table1] [--hs ...]` | `SLICE_S5_TSTAR_<label>.json` |
| `near_root.py` | analyticity of the complete jet files at every tangency root in (0, 1) from values at t* -+ 1e-3, 1e-5 (uses the `file_gate.py` cache) | `near_root.py [--labels quark,gluon]` | `SLICE_S5_NEARROOT.json` |
| `landau_family.py` | which denominator atoms (Landau family) generate each tangency letter, per group | `landau_family.py [--channels ...] [--all-sigma] [--all-letters]` | `SLICE_S5_LANDAU.json` |
| `n4_alphabet.py` | prime support of the N=4 representation versus the external slice alphabet `data/N4_SLICE.json` | `n4_alphabet.py [--labels n4]` | `SLICE_S5_N4ALPHABET.json` |

Symbol verdicts are meaningful only on COMPLETE files and at two point pairs with root-continuation labeling; "absent from the weight-w
symbol" is not "absent at weight w" (i pi and zeta terms are symbol-invisible).

## 8. Numeric reference evaluators

These share no code with the symbolic engine beyond the cell data and give certified values (two independent routes in Arb ball
arithmetic; certified digits = -log10 of their relative difference).

* `gate/scripts/oracle.py --channel q_jet|g_jet|<channel> --shape 'u=3/7+2/5j;v=-1/3+4/9j' --digits 30 [--nproc 3]`: the E4C at an
  arbitrary shape (points {0, 1, u, v} with rational real and imaginary parts): exact rational-function assembly per (channel, chart,
  sigma), inner residues in Arb, two-dimensional exp-sinh, routes A and B; per-task results are appended to
  `$SURD_WORK/oracle_ckpt/<shapekey>.jsonl` (rerun = resume; `--assemble-only`). Core: `oracle_core.py`; cells: `oracle_cells.py`
  (`--e4c` adds the optional float cross-check against `$SURD_E4C_ENGINE`).
* `slice/scripts/slice_oracle.py --t LIST [--channels LIST] [--nproc k] [--t-major] [--assemble-only]`: the same machinery on the slice
  at rational t, per channel and per jet (n_f = 5), in the normalization of the assembled files; writes `SLICE_S0_ORACLE_t<p>_<q>.json`
  per t and the table `SLICE_S0_ORACLE.json`; checkpoint per record.
* `slice/scripts/slice_po.py --channel C --sigma S --pair P --t T`: the per-group piece evaluator (two charts) used by the stage-3
  check; cache `ckpt/po/`.
* `slice/scripts/file_gate.py [--labels quark,gluon,n4] [--s LIST] [--shard k/n] [--e4c]`, then `--receipt`: assembled files versus every
  reference record; per-(file, t) cache `ckpt/file_eval/`; summary record `SLICE_FILE_GATE.json`.
* `gate/scripts/onefold_arb.py --piece PIECE.json --chart g<gauge> [--digits 32]`: the last integration of a numeric-kinematics
  HyperFLINT one-fold object in ball arithmetic on a rotated contour, closed against its two-chart piece reference; also holds
  `fibre_residue`, the residue-route fiber function of the stage 1-2 check. The piece file is JSON: `channel`, `sigma` (e.g. `"1234"`),
  `pair_group` (`"34"` or `"all"`), `d_sigma` (`{"12": "p/q", ...}`, the six labeled distances as rationals), `onefold` =
  `{"g<gauge>": {"gauge": int, "order": [k, j, l], "response": path}}` where `response` is a HyperFLINT `eval-json` response produced
  with `nint = 2` (the last variable left symbolic; `hf_piece.py --nint 2 --run` builds and runs the request with the engine of
  `tools/subtropica`), and optionally `oracle` = `{"A": {"value_mid": "..."}, ...}` with reference values.
* `gate/scripts/hf_piece.py --channel C --shape u,ub,v,vb --gauge g --order kjl --tag T [--pair kj] [--nint 2] [--run]`: HyperFLINT
  request builder per piece (library part: `first_pair`, `group_cells`); files under `$SURD_WORK/hf_piece/`.
* `gate/scripts/hlog_eval.py` (G(w; 1) and ZIP_reg words with algebraic letters; `selfgate(dps)` self-check) and `slice/scripts/hpath.py`
  (path-deformed evaluation with the -i0 rule) are the hyperlog evaluators behind every numeric comparison of the closed forms.

## 9. Environment variables

All optional except where a tool says otherwise.

| variable | meaning | default |
|---|---|---|
| `SURD_WORK` | output root: run records `SLICE_*.json` (+ `superseded/`), `ckpt/`, `logs/`, `out/`, `cells/`, `oracle_ckpt/`, `hf_piece/` | `./surd_work`; a root inside the package tree is refused |
| `SURD_DATA` | data root holding assembled files `out/*.json.gz` and the records the structure tools and tests read | WORK |
| `SURD_TERMS` | directory of splitting-function term lists `<name>.terms.pkl`, read only when a channel's cells are not cached | `<package>/data/terms` |
| `SURD_LR` | stage-alphabet file | `gate/lr/out_groups_symbolic.json` |
| `SURD_N4SLICE` | external N=4 slice alphabet | `data/N4_SLICE.json` |
| `SURD_FIBRATE` | directory holding `zip_num.py` | `<repo>/tools/subtropica/scripts/fibrate` |
| `SURD_E4C_ENGINE` | directory of an external float engine `e4c.py` for the optional float cross-checks of `oracle_cells.py --e4c` and `file_gate.py --e4c` | unset (checks skipped / refused with a message) |
| `SURD_VMEM_KB` | address-space cap applied by the shell drivers (`ulimit -v`) | 8388608 |
| `SUBTROPICA_HF_BIN`, `SUBTROPICA_MZV_DATA` | HyperFLINT engine of `tools/subtropica`, needed only by `hf_piece.py --run` and for the piece files `onefold_arb.py` consumes | unset |
| `HF_LR_TIME_BUDGET_S` | time budget in seconds that `hf_piece.py` passes to the HyperFLINT linear-reducibility step | 3600 |

The drivers set `OMP_NUM_THREADS=1` and run under `nice`; the shell drivers apply the address-space cap.

## 10. Outputs, record format, file format and package data

Per group: `ckpt/<tag>/stage{1,2,3}.pkl` (exact pickles) and the run records `SLICE_S12_<tag>.json`, `SLICE_G12_<tag>.json`,
`SLICE_S3_<tag>.json`, tag = `<channel>_s<sigma>_p<pair>_g<gauge>`. Assembled: `out/E4C_LO_*_dipole_slice_<label>.json.gz` (+ `.txt`
summaries) and `SLICE_ASSEMBLY.json`; structure tools write `SLICE_S5_*.json`; the reference evaluator writes `SLICE_S0_ORACLE_t<p>_<q>.json`;
the file check writes `SLICE_FILE_GATE.json`.

Every record is JSON with the tool's results plus a `producer` block: script path and sha256, launch line, cwd, host, pid, the sha256 of
every input file, the sha256 of the engine modules and of every `slice/scripts` module at import (and the list of modules whose hash
changed during the run), package versions and a UTC stamp. Records are write-once: writing a record whose file already exists first moves
the old one to `superseded/<name>_<sha12>.json` and appends a line to `superseded/LEDGER.txt`; writes are atomic (temporary file + rename).
Check records carry `min_agree_digits` and a boolean `PASS(>=30)` computed from the rows, never a typed status.

The format of the assembled closed-form files (value = sum coef(t) * prod rho^e * prod const * prod Z(w); coefficients (A+iB)/D in
Q(i)(t); root-slot symmetric), their evaluation and precision notes and the census of the E4C reference files are in `FORMAT.md`. The
package data (the three cell pickles of the smoke-test channel, the stage-alphabet file, the N=4 slice alphabet; formats, provenance,
and the term-list format for building cells of other channels) are in `DATA.md`.

## 11. Branches, base point, precision

* Branch rule: every hyperlog letter that lies EXACTLY on the positive real integration path is taken as `letter - i0` (the contour passes
  above it), at every stage and in the final evaluation (`hpath.py`: path concatenation 0 -> m -> 1 with m = 1/2 + ic above the near-path
  letters, segment integrals as ordinary G(.; 1)); all other hyperlogs on their principal branch. The rule is validated, not assumed:
  stage 1-2 checks (43-51 digits on groups with intrinsically on-path letters), stage-3 piece checks (>= 35 digits), the file check at
  every reference t, and reference comparisons >= 33 digits in every on-path configuration interval of (0, 1) down to t = 1/20 (the
  configurations change at the real roots of the letters' discriminants and leading coefficients).
* Base point x0 = 11/2 for the fibration constants (G_0(v; x0), ZIP words at x0); a second base point (7/2) reproduces the same function to
  42 digits; base-point polynomials are representation letters, not letters of the function.
* Precision of `assemble.evaluate_file(path, t, dps)`: the exact value is real, so |Im F| is the error monitor. Digits delivered ~ dps minus
  the termwise cancellation depth: ~10 generically, ~25 near t = 4/5, up to ~30 within ~1e-2 of the roots 0.4110, 0.8110 of 15t^2 -+ 6t - 5
  (poles of rational hyperlog points; F itself is analytic) and near the Table-1 roots 0.8178, 0.9511 (q_qbqgq at t = 163/200: 34.2 digits at
  dps 60, 37.4 at dps 80). Raise dps there, or use `evaluate.py`, which escalates dps until |Im F|/|F| meets the target.
* Representation-singular rational t: `evaluate_file` can raise `ZeroDivisionError` where a letter coefficient, a rational-point denominator
  or a coefficient denominator vanishes exactly while F is analytic (t = 3/5 for 7 of the 9 QCD files; t = 5/8 for the files carrying the
  8t -+ 5 class). `evaluate.py` returns the symmetric limit (F(t + delta) + F(t - delta))/2 there, escalating delta until the nearly colliding
  roots clear the 1e-12 path-distance guard of `hpath.py` (delta = 1e-20 at 5/8), with the O(F'' delta^2) error reported. Per-file lists are in
  `FORMAT.md`.
* Domain: files are constructed and checked on 0 < t < 1; F(1/t) = t^6 F(t) gives t > 1.
* "Table-1 letters": the twelve rational tangency letters P(t) predicted for the dipole slice by the Landau analysis of the E4C
  (`letters.TABLE1`: 8t +- 5, 5t +- 8, 45t^2 +- 6t - 35, 35t^2 +- 6t - 45, 2t^2 - 1, t^2 - 2 and two quartics). The structure tools ask
  whether they are letters of the integrated functions.

## 12. Resource profile (measured, single thread)

Per group stages 1-3 <= 36.6 min wall (the g_gggg first-pair-14 groups are the worst; median minutes), <= 0.85 GiB RSS; jet assembly
5.10 GiB (quark, 156 groups, 18.8 min) / 4.98 GiB (gluon, 114 groups, 17.5 min); `evaluate_file` at dps 60: 95-250 s and 1.0-1.7 GB per QCD
file, ~150 s / 0.45 GB for n4; full reference values for the 7 QCD channels per t ~2 CPU-h (29 min on 4 processes, 0.53 GB), n4 per t
~27 min on 3 processes; `symtest` gluon ~9 min / 1.7 GB, quark ~3.6 min / 1.15 GB; disk for the full problem 4.8 GB (stage-3 pickles
~2.3 GB, `out/` 0.5 GB). The smoke-test group (q_qbpqpgq, sigma 1234, pair 34): stages 1-2 49 s, stage 1-2 check 3 s, stage 3 with its
check 89 s, 0.16 GB.

Per channel (all groups; stage 1-2 check = min over the 8 (t, x) points, stage-3 check = min over t = 1/2, 4/5):

| channel | groups | stage-1 N terms | stage-2 monomials | nonlinear last-variable letters (deg: count) | stages 1-2 wall s (median, max) | stages 1-2 RSS MB max | stage 1-2 check digits (median, min) | stage-3 terms (median, max) | stage-3 wall s (median, max) | stage-3 RSS MB max | stage-3 check digits (median, min) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| q_qbpqpgq | 48 | 264-4176 | 8526-76300 | deg2: 4, deg3: 8 | 65, 125 | 183 | 49.5, 36.0 | 609, 1521 | 14, 368 | 236 | 48.1, 42.4 |
| q_qbqgq | 36 | 3397-7882 | 87758-183392 | deg3: 12, deg2: 30 | 174, 434 | 290 | 44.6, 36.3 | 1436, 3072 | 47, 502 | 279 | 44.0, 41.2 |
| q_gggq | 72 | 1346-11812 | 70557-176578 | deg3: 40, deg2: 20 | 166, 604 | 341 | 38.1, 35.8 | 1704, 2020 | 307, 484 | 279 | 37.3, 35.2 |
| g_qbpqpqbq | 36 | 561-7484 | 36349-108836 | deg3: 24, deg2: 12 | 121, 305 | 237 | 39.0, 33.9 | 1109, 1336 | 373, 571 | 247 | 42.4, 40.6 |
| g_qbqqbq | 18 | 2694-22998 | 91609-176757 | deg3: 12, deg2: 6 | 281, 483 | 358 | 37.8, 34.3 | 1365, 1821 | 530, 655 | 285 | 41.5, 40.2 |
| g_qbggq | 24 | 4428-30708 | 60514-184700 | deg3: 16, deg2: 24 | 278, 728 | 408 | 38.2, 36.2 | 1952, 2155 | 590, 757 | 509 | 42.0, 40.9 |
| g_gggg | 36 | 427-31551 | 22371-247668 | deg2: 84, deg3: 12 | 407, 1096 | 571 | 37.8, 35.9 | 2380, 3903 | 105, 831 | 866 | 37.4, 35.8 |
| n4 | 72 | 26110-63798 | 25349-39457 | deg3: 32, deg2: 16 | 34, 126 | 285 | 38.6, 38.0 | 1976, 2203 | 8, 48 | 360 | 39.7, 38.0 |

No stage 1-2 or stage-3 check is below 30 digits (the n4 groups are checked per sigma class against the full reference values instead
of the per-group piece evaluator).

## 13. Validation

* **Full numeric reference values** (S0, `slice_oracle.py`; certified digits from the two-route agreement in parentheses). F_ch(t) and
  the jets, n_f = 5:

  | t | q_qbpqpgq | q_qbqgq | q_gggq | g_qbpqpqbq | g_qbqqbq | g_qbggq | g_gggg | n4 | quark jet | gluon jet |
  |---|---|---|---|---|---|---|---|---|---|---|
  | 1/2 | 39.8231429717647981229 (37) | 9.85332338311198121274 (36) | 299.670398801443081256 (36) | 13.9533678453451131882 (36) | 3.45318361375936036833 (36) | 441.663260015706397713 (35) | 1446.07717306700165238 (35) | 305.864446892640099545 (37) | 349.34686515631986059239 (36) | 1905.1469845418125236546 (35) |
  | 2/3 | 31.0140931016272659854 (37) | 7.66208136245655829922 (36) | 205.796055212522436393 (36) | 12.5348510266286536938 (36) | 3.11535468308721894604 (36) | 326.957465395077078289 (35) | 986.259035164963288085 (35) | 211.834937557978234585 (37) | 244.47222967660626067863 (36) | 1328.8667062697562390154 (35) |
  | 4/5 | 24.0463148979389746626 (37) | 5.93656298612860891797 (36) | 147.709212022406279337 (36) | 10.2379963077283493558 (36) | 2.54895740141279899990 (36) | 245.400009091309171044 (35) | 706.291992478400918179 (35) | 152.979270122224247581 (37) | 177.69208990647386291832 (36) | 964.47895527885123757983 (35) |
  | 9/10 | 18.8585477179829417601 (37) | 4.65470983430885620744 (36) | 112.236820653719917376 (36) | 8.16852261432771022922 (36) | 2.03493484105280994996 (36) | 189.925270729070171228 (35) | 536.308342121248789996 (35) | - | 135.75007820601171534378 (36) | 736.43707030569948140417 (35) |
  | 3/2 | 2.72277360562982856388 (36) | 0.67266557914570607839 (36) | 18.0671433931432591621 (36) | 1.10045331372322885652 (35) | 0.27350164570313033271 (35) | 28.7040847534772743628 (35) | 86.5851553505591912722 (35) | - | 21.462582577918793804434 (36) | 116.66319506346282482440 (35) |
  | 5/2 | 0.19619820083108279816 (35) | 0.04858506739510899039 (35) | 1.56936041495046782316 (36) | 0.05947144856071119094 (34) | 0.01467029232602645624 (34) | 2.21540738409071676934 (35) | 7.63101288603668063320 (35) | - | 1.8141436831766596117170 (35) | 9.9205620110141350497361 (35) |

* **File check** (assembled files vs the reference values, minimum over all reference t; t > 1 through F(T) = s^6 F(s)): quark 35.8, gluon
  34.7, n4 40.2, q_qbpqpgq 40.6, q_qbqgq 39.1, q_gggq 37.9, g_qbpqpqbq 33.8, g_qbqqbq 32.8, g_qbggq 34.2, g_gggg 36.1 digits; every file
  complete.
* **N=4 control against a known closed form**: F_n4(t) / Re E4C_CMSYZ(t) = 5/24 (5 = recoiler sum of the squared five-point form factor,
  1/24 = S_n4) with |24 F/(5 C) - 1| = 1.8e-27, 9.9e-27, 1.0e-26 at t = 1/2, 2/3, 4/5 and <= 2e-26 at t = 3/7, C evaluated independently from
  the GPL expressions of arXiv:2401.06463 (C(1/2) = 1468.1493450846724778184508646, C(2/3) = 1016.80770027829552600981296509, C(4/5) =
  734.300496586676388389623415671, C(3/7) = 1734.992090049981175025454775084718202052); F_n4(1/2) =
  305.86444689264009954551059734352338942858597093132. The prime support of the N=4 representation equals the 34 CMSYZ slice letters
  (34 of 34, none extra, none missing).
* **Slot invariance**: the value of a file is unchanged to 55 digits under reversal or cyclic permutation of every root list (the term
  store lives in the symmetric splitting algebra); base point 7/2 instead of 11/2 reproduces the function to 42 digits.
* **Symbol tests** (`symtest.py`, weight 3, two point pairs (e/7, pi/9) and (1/sqrt7, pi/9), fixed root labeling by canonical order plus
  continuation): in the quark jet the two genuine square classes (25t^2 + 7)(t + 1)(t - 1)(t^2 + 1) and -(7t^2 + 25)(t + 1)(t - 1)(t^2 + 1)
  SURVIVE and every Table-1 class CANCELS (functionals <= 1e-58); n4 carries only the two genuine classes; identical verdicts at both point
  pairs; every LLL relation found at t* is verified at t**.
* **Analyticity at the Table-1 roots** (`near_root.py`, complete jet files at t* +- 1e-3 and +- 1e-5 around every Table-1 root in (0, 1)):
  odd-part ratio D(1e-3)/D(1e-5) = 100 to within 6e-6 relative, even-part change |E(1e-3) - E(1e-5)|/|F| = 1.5e-6 .. 3.7e-6, |Im F| <= 7e-30
  on both sides, i.e. no (t - t*)^{+-1/2} or (t - t*)^{-3/2} component on the principal sheet (a square-root cusp would give ~3e-2).
* **Representation-singular points** (`evaluate.py`): quark at t = 3/5 (symmetric limit) 282.5818179378116526279331062518124624777, n4 at
  3/5 (direct) 245.87527586150966006840229747961171011, q_qbqgq at 5/8 (symmetric limit, delta = 1e-20)
  8.184320525040340232233572530595188909, each >= 30 digits against direct evaluations nearby.

## 14. Limits

* One symbolic kinematic parameter only: four-variable intermediate objects reach 10^6-10^7 monomials and do not finish (measured).
* Three integrations; the first two eliminated variables must be Brown-linear (no change of variables or root-rationalization is
  attempted); the elimination order is fixed per first pair, not searched.
* Principal sheet only, with the -i0 rule for on-path letters: no homotopy or monodromy of the hyperlogarithms, so the coefficient of a
  Landau singularity on a pinch sheet is out of reach; the local Laurent capacity bounds exponents but does not give that coefficient.
* `symtest.py` treats square classes only; there is no Galois symbol test for the cube-root (cubic) letters.
* Kinematic invariants must be polynomial in t on the slice; Euclidean slices force the Gaussian ring Q[x,t,j]/(j^2+1), and a
  real-rational shortcut silently loses terms.
* Numeric evaluation of the assembled files loses digits to termwise cancellation (section 11) and raises at representation-singular
  rational t; `evaluate.py` handles both, at a cost in time.
* The integrand interface is the atom family of the worked problem (section 5); other families need a `build_group` replacement and
  their own cells for the numeric checks.
