# Surd — package data

Five data files come with the package. All are plain data, produced by this package's scripts from the sources named below (or, for
`N4_SLICE.json`, by an exact computer-algebra reduction described below); none contains code, and the pickles hold only built-in Python
types. The three cell pickles are derived from a file by V. Del Duca, C. Duhr, R. Haindl, A. Lazopoulos and M. Michel whose source
states no license; see "Source of the cell data and its license status" below.

| file | size | sha256[:16] | what it is |
|---|---|---|---|
| `gate/cells/q_qbpqpgq_g1.pkl` | 1.0 MB | `896d5e12ea8e8b42` | exact Cheng-Wu cells of the channel q -> q̄′q′g q (chart x_1 = 1), derived from the file `P_qbp_qp_g_q.m` of Del Duca, Duhr, Haindl, Lazopoulos and Michel (source and license status below) |
| `gate/cells/q_qbpqpgq_g2.pkl` | 1.0 MB | `a19c79c604dd8594` | same channel, chart x_2 = 1 |
| `gate/cells/q_qbpqpgq_g4.pkl` | 1.0 MB | `a74ed433e8756dcb` | same channel, chart x_4 = 1 |
| `gate/lr/out_groups_symbolic.json` | 1.2 MB | `9b492c1df59a9a98` | per-support-group stage alphabets from a linear-reducibility run (all seven QCD channels and the N=4 control) |
| `data/N4_SLICE.json` | 56 kB | `cc3bd6b61fec4a1d` | the N=4 SYM slice alphabet used as the external reference of `slice/scripts/n4_alphabet.py` and `tests/run_tests.py` |

## The cell pickles (`gate/cells/<channel>_g<gauge>.pkl`)

Produced by `gate/scripts/oracle_cells.py build_cells(channel, gauge)` from the term list of the tree-level 1 -> 4 splitting function of the
channel: every term is mapped to the Cheng-Wu chart x_gauge = 1 and collected over its denominator multiset, keeping the kinematics symbolic
through the six pair distances d_ab (a < b). Each file is a `pickle` (protocol 4) of one `dict` with keys

| key | Python type | meaning |
|---|---|---|
| `channel` | `str` | channel name (`q_qbpqpgq`) |
| `gauge` | `int` | the chart variable set to 1 |
| `xvars` | `list` of 3 `int` | the remaining chart variables, increasing |
| `dshift` | `int` | global shift of the d exponents (all stored exponents are >= 0; the evaluator divides by prod d^dshift) |
| `hz` | `int` | degree of homogeneity of the term list in the energy fractions z (1 for quark-parent files) |
| `n_terms`, `n_cells`, `n_monomials` | `int` | census: input terms, cells, (x, d)-monomials |
| `build_s` | `float` | build wall time in seconds |
| `cells_raw` | `list` of `(denominator_key, numerator)` tuples | the cells |

where `denominator_key` is a `tuple` of `(atom_name: str, power: int)` pairs over the atoms {`x1`..`x4`, `Om` = x1+x2+x3+x4,
`ell<S>` = sum_{j in S} x_j, `s<S>` = sum_{a<b in S} d_ab x_a x_b}, and `numerator` is a `dict` mapping an x-exponent `tuple` of 3 `int`
(powers of the three chart variables) to a `dict` mapping a d-exponent `tuple` of 6 `int` (powers of d12, d13, d14, d23, d24, d34, shifted
by `dshift`) to a `(numerator: int, denominator: int)` pair. The only Python types present are `dict`, `list`, `tuple`, `str`, `int`,
`float`. `oracle_cells._hydrate` turns the numerators into `fmpq_mpoly` objects at load time. The three charts of one channel describe the
same integrand; `tests/run_tests.py` (test `g`) checks that identity exactly at random rational points.

Cells for the other channels (`q_qbqgq`, `q_gggq`, `g_qbpqpqbq`, `g_qbqqbq`, `g_qbggq`, `g_gggg`, `n4`) are not included (about 40 MB).
On a cache miss `build_cells` reads the channel's term list from `$SURD_TERMS/<name>.terms.pkl` and writes the new cell pickle
under `$SURD_WORK/cells/`. A term list is a pickle of a `(terms, names)` pair: `names` is a `list` of factor names (`CF`,
`CA`, `D`, `gg`, `z<S>`, `s<S>`, `k<i>k<j>`) and `terms` a `dict` mapping a `tuple` of `(name_index: int, exponent: int)` pairs to a
rational coefficient, i.e. the splitting function as a sum of monomials in those factors. The file names expected per channel are the
`CHFILES` table in `oracle_cells.py`; the factor names are the notation of the source files.

### Source of the cell data and its license status

The QCD term lists are the quadruple-collinear splitting functions of V. Del Duca, C. Duhr, R. Haindl, A. Lazopoulos and M. Michel,
"Tree-level splitting amplitudes for a quark into four collinear partons", JHEP 02 (2020) 189 [arXiv:1912.06425], and "Tree-level
splitting amplitudes for a gluon into four collinear partons", JHEP 10 (2020) 093 [arXiv:2007.05345]. The quark-parent paper states
that the explicit results "are too long to be recorded in this paper" and "are made available in computer-readable form". The authors
provide them as Mathematica files on the web page "Quadruple Collinear Splitting Amplitudes. Supplementary material from arXiv:1912.06425 and
arXiv:2007.05345", https://people.phys.ethz.ch/~pheno/quadruple_collinear. The seven QCD file names in the `CHFILES` table are file
names on that page.

The three shipped cell pickles are derived data of one of those files: `files/P_qbp_qp_g_q.m` on that page (listed there as
"Splitting: q->qb'(1) q'(2) g(3) q(4), total, eq. (3.5), LHS"; 1,810,350 bytes; sha256
`8143dfd031a53647e0f7e9f2b223e10fe97e372beba05c1cd3dfff81786bee59` as fetched on 2026-10-13). The pickles are a mechanical
re-expression of the authors' term list (28,418 input terms): the color and dimension constants are substituted (CF = 4/3, CA = 3,
D = 4), the energy fractions and invariants are mapped to the Cheng-Wu chart x_gauge = 1 as described above, and the terms are
collected over their denominators. The atom polynomials in `gate/lr/out_groups_symbolic.json` are denominators only (the atoms x_j,
Om, ell_S, s_S of the cell integrands, for the seven QCD channels and the N=4 control).

License status, as found on 2026-10-13. The web page and the file carry no license, copyright or terms-of-use statement. The journal
versions of the two articles are published under CC BY 4.0 and the arXiv records under the arXiv non-exclusive distribution license;
the arXiv records have no ancillary files, and whether the journal records carry the Mathematica files as supplementary material has
not been established. The file is therefore not known to be covered by the articles' CC BY 4.0 license. The Anthropic copyright and the
MIT and CC BY 4.0 grants of this repository do not cover the cell pickles; the result they re-express is the work of its authors,
named above.

### The N=4 control term list (`P_N4_coll.m`, not included)

The `n4` entry of `CHFILES` names a term list `P_N4_coll.m` in the same format. It holds the N=4 1 -> 4 collinear integrand, the
collinear limit of the squared five-point form factor (normalization <P>_N4 = s1234^3 P_N4). The file is not included in the package
and neither are cells built from it. It is not a file of the web page above, and it is not one of the ancillary files of
arXiv:2401.06463. This release does not record which publication or computation the file was produced from. The N=4 reference
function against which the control is compared is that of D. Chicherin, I. Moult, E. Sokatchev, K. Yan and Y. Zhu, arXiv:2401.06463
(see the last section).

## The stage alphabets (`gate/lr/out_groups_symbolic.json`)

One JSON object: `description`, `algorithm` (compatibility-graph reduction, Fubini), `frame` (the symbolic split-signature frame
w = (u,1,0,v), w̄ = (ū,1,0,v̄) and the six d_ab as products in u, ū, v, v̄), `max_degree`, `lr_refine` (the sibling-package member that
produced it) and `runs`: a list of 276 records, one per (channel, chart `gauge`, first-pair `group`, elimination `order`/`perm`), each
with the `atoms`, the `input_polys` (atom polynomials in the cell variables and frame symbols, with a sha256 of their canonical text), the
`verdict`, and `stages`: for stage k the polynomial `letters` alive after k integrations, the `next_var`, the admissible and blocking
next variables. `symbolic/scripts/alphabet.py` seeds its letter table from `stages[k].letters` of the matching run; every denominator the
fibration engine produces is factored and looked up there, and a factor outside the table raises `OffAlphabet`. Produced by running the
linear-reducibility refinement of `tools/subtropica` (`src/lr_refine.jl`) on the cell integrands of `oracle_cells.py`.

## The N=4 slice alphabet (`data/N4_SLICE.json`)

The 63 GPL/q symbol letters of the N=4 SYM 1 -> 4 collinear limit (ancillary files of D. Chicherin, I. Moult, E. Sokatchev, K. Yan and
Y. Zhu, arXiv:2401.06463; CC BY 4.0, adapted), restricted to the dipole slice over
the 12 labeling classes and reduced to 34 irreducible Q[t]-norm factors: `letters` maps each canonical primitive polynomial in t to its
`degree`, discriminant `disc`, whether it descends from a q letter (`from_q`), the number of source letters and a sample of them;
`degenerate_on_slice` lists pre-slice letters that become 0 or constant per class; `field_discs` are the quadratic field discriminants that
occur; `gates` records the exact identity checks the reduction passed. Produced by an exact polynomial-algebra reduction (norms L·conj(L) over
Q(i)[t], root letters through the exact cubic resultant), independent of this package's integrator; it is the external reference against
which `n4_alphabet.py` compares the prime support of the N=4 closed form (34 of 34).
