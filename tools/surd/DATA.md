# Surd — package data

Five data files come with the package. All are plain data produced by this package (or, for `N4_SLICE.json`, by an exact
computer-algebra reduction described below); none contains code, and the pickles hold only built-in Python types.

| file | size | sha256[:16] | what it is |
|---|---|---|---|
| `gate/cells/q_qbpqpgq_g1.pkl` | 1.0 MB | `896d5e12ea8e8b42` | exact Cheng-Wu cells of the channel q -> q̄′q′g q (chart x_1 = 1), derived from the splitting function of Del Duca, Duhr, Haindl, Lazopoulos and Michel (CC BY 4.0, adapted; see below) |
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
rational coefficient, i.e. the splitting function as a sum of monomials in those factors (QCD: the quadruple-collinear splitting functions of
V. Del Duca, C. Duhr, R. Haindl, A. Lazopoulos and M. Michel, JHEP 02 (2020) 189 [arXiv:1912.06425] and JHEP 10 (2020) 093
[arXiv:2007.05345]; results (c) the authors, CC BY 4.0 (journal versions, JHEP / SCOAP3), adapted as described above; N=4: the 1 -> 4
collinear limit of the squared five-point form factor). The shipped cell pickles are therefore derived data of the quark-parent result and
carry its attribution (THIRD_PARTY.md section B and NOTICE at the repository root). The file names expected per channel are the `CHFILES`
table in `oracle_cells.py`.

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
