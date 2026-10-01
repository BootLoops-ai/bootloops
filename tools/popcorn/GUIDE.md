# POPCORN — certified population-genetics likelihoods

Tool page: https://bootloops.ai/tools/popcorn.html

KIND: package (`import popcorn`; CLI `certsfs.py`; battery `selftest.py`;
flat layout — the engine modules ship beside the package files and are
aliased in place by identity).
NAME NOTE: unrelated to Popcorn, the trans-ethnic genetic-correlation estimator
of Brown, Ye, Price & Zaitlen (2016, Am. J. Hum. Genet. 99:76;
github.com/brielin/Popcorn).

PURPOSE: The Poisson Random Field selection-SFS/DFE stack of the
polyDFE/fitdadi/fastDFE class made exact and certified:

- **Selection SFS, exact**: M_i = 1F1(n-i; n; -S) for the whole vector via ONE
  exact-rational tridiagonal solve — M_i = alpha_i + beta_i e^{-S} with alpha,
  beta exact Fractions (n <= 2000 practical); exact gradient stack, one extra
  solve per order, same matrix. Certified-precision assembly: the Fraction
  heights BOUND the cancellation a priori, so the assembly dps is computed,
  not guessed.
- **Selection SFS, biobank-n**: stable-direction contiguity recursion in arb
  ball arithmetic (python-flint) — ball-certified vectors + gradients to
  n = 10^5; radii track the directional contamination exactly, and precision
  doubles until the target is met.
- **Certified DFE mixing** (`popcorn.dfe.CertifiedKernel`): polyDFE model-C
  class mixtures with a certified kernel cache — likelihood evaluation is
  matrix-vector fast while every number carries a proof budget (GL degree-pair
  self-check, analytic S->0 and S->infinity tail bounds, exact analytic
  parameter gradients).
- **Dominance h != 1/2** (`popcorn.dominance`): E(S,h) is an entire
  q-deformation of the Kummer line. A certified 2-fold nested-quadrature
  oracle (degree-pair self-checks at both levels, h=1/2 collapse gate) plus a
  stable q-series evaluator, cross-checked against pinned reference values.
- **Exact certificate family** (`popcorn.certificates`): positivity of
  rational polynomials on [0,1] (two independent routes), pure-rational LP
  membership with Farkas witnesses, float-propose/exact-decide LP, exact
  integer-cone dual simplex, certify-after-float hull membership with
  certificate reuse (CertHull), and exact-rational REGION certificates
  (Bernstein box positivity, certified sup over a chart, Cauchy-Schwarz
  chi^2 lower bounds, rational reconstruction + Sturm root counting).
- **Exact Lambda-coalescent machinery** (`popcorn.lambda_coalescent`):
  merger rates lambda_{b,k} for Kingman, Beta(2-alpha, alpha) and Dirac(psi)
  coalescents (and the Kingman+Dirac two-atom mixture) as exact rationals, and
  the exact expected branch-length spectrum E[L_i], i = 1..n-1 (hence the
  expected unfolded SFS) by the standard first-transition recursion, all in
  Fractions with built-in exact identity checks (Kingman 2/i closed form,
  alpha=2 -> Kingman, alpha=1 -> Bolthausen-Sznitman, psi=0 -> Kingman,
  psi=1 -> star); ships an n = 20 exact reference grid (Kingman, 7 Beta,
  6 Dirac) and an exact linear functional separating that grid from the
  variable-population-size Kingman class.
- **Two-locus branch-length moments** (`popcorn.twolocus`): the
  second moments E[T_i^A T_j^B] (i, j = 1..n-1) of the branch lengths
  subtending i samples at locus A and j at locus B, for a sample of n under
  the coalescent with recombination (scaled rate rho = 4 N_ref R, each doubly
  ancestral lineage splitting at rho/2) and piecewise-constant N(t); i.e. the
  expected joint two-locus frequency spectrum up to scale, plus the Kingman
  first moments. Deterministic: the labeled two-locus configuration chain is
  enumerated once per n, the terminal epoch is one sparse LU solve and each
  finite epoch one expm_multiply — no simulation, no time stepping.
- **Transient selected SFS, large samples** (`popcorn.transient`): the
  expected unfolded SFS E_i (i = 1..n-1, theta = 1) under genic selection S
  after a piecewise-constant size history [(nu, T), ...] following an
  ancestral Wright equilibrium, by direct integration of the Kimura forward
  diffusion in u = x(1-x)f — Scharfetter-Gummel exponentially fitted fluxes
  on a fixed log/lin/log grid, Crank-Nicolson with Rannacher startup, banded
  LAPACK solves with one factorization per epoch, binomial projection to
  sample size n (n ~ 1000-2500 design range) — plus hypergeometric down-sampling
  n -> m (formula-exact, float) and folding. Not a moment closure: n enters only
  through the projection. FLOAT-VALIDATED register (measured trust radius,
  see FOOTGUNS), never certified.
- **Certified enclosures of the TRANSIENT selected SFS** (`popcorn.enclosure`,
  engine `transient_enclosure_engine.py`): the expected unfolded sample SFS
  E_n(i; t), i = 1..n-1 (theta = 1), under genic selection S <= 0 and a
  piecewise-constant size history [(rho, T), ...] after an ancestral Wright
  equilibrium. The moment system of the forward diffusion, truncated at M >= n,
  is integrated through the epochs in ball arithmetic (Arb via python-flint:
  exact rational matrices, rigorous matrix exponential and solve per epoch),
  and projected to sample size n with exact integer coefficients, so every
  reported entry is an interval guaranteed to contain the model's value.
  Truncation is closed two-sidedly by a proved lemma (for S <= 0 the
  truncated matrix is Metzler and 0 <= w_{M+1} <= w_M, so the systems with
  w_{M+1} := 0 and w_{M+1} := w_M bracket the untruncated moments
  entrywise and the reported ball is their hull; at S = 0 the closure is
  exact), and independently controlled by comparing truncations: the
  shipped M = 200 and M = 500 enclosures agree to all 30 recorded digits at
  every shipped cell, with relative widths 1e-76 .. 1e-64. Register:
  CERTIFIED-ENCLOSURE, model-conditional (rigorous for the stated diffusion
  model and truncation; not a bound on model error). Use it to certify a float
  transient solver (`popcorn.transient`) at spot cells, exactly as the package
  certifies its float stationary paths against the exact stationary engine.
- **Ancestral-state join** (`popcorn.ancestral`, engine `epo_join.py`):
  stream a sites table (chrom pos ref alt, or pos ref alt) against an
  Ensembl-EPO-style ancestral-allele FASTA, classify each site by EPO
  confidence (uppercase = high, lowercase = low, `N`/`-`/`.` = no call),
  orient REF/ALT to ancestral/derived over biallelic SNVs, and write the
  polarized `pos ref alt anc conf` table with counters (confidence
  histogram, anc==ref / anc==alt / third-allele mismatch by confidence,
  unpolarizable) and the summary rates. Pure Python, standard library;
  measured about 0.6 million sites per second on a 5 Mb synthetic ancestor.
- **Exact linked two-site frequency spectrum and class certificates**
  (`popcorn.twosfs`): the exact second moments E[L_i L_j] (i, j = 1..n-1)
  of the branch-length spectrum, hence the expected two-site frequency
  spectrum of completely linked sites up to (theta/2)^2 and its scale-free
  normalization q_ij = E[L_i L_j]/E[L_tot^2], for any exchangeable
  Lambda-coalescent (Kingman, Beta(2-alpha, alpha), Dirac, Kingman+Dirac
  mixture) by an exact first-step recursion on block-size multisets, with two
  independent exact cross-routes (raw set partitions; for Kingman the
  holding-time x jump-chain decomposition), all in Fractions. On the Kingman
  side, the exact block-counting semigroup and the two-time block-count
  kernel Csym_ij(x, y) (sympy, rational coefficients) represent the 1- and
  2-SFS of EVERY single-population Kingman coalescent with a deterministic
  size history as integrals against one measure; on top of that,
  `twosfs_certificates` decides, exactly, whether a given model's 2-SFS lies
  outside the 2-SFS set of every such history and every pooling of them
  (a rational witness W with sum W_ij Csym_ij >= 0 on the unit square,
  certified by exact tensor-Bernstein subdivision, and <W, M> < 0), outside
  all poolings of single-atom histories only, or inside (exact atom
  mixture). A float LP proposes witnesses; only exact arithmetic decides.
  Ships exact reference rows for n = 3..6: Beta and Dirac coalescents whose
  1-SFS is provably a variable-size-Kingman mixture yet whose 2-SFS is
  certified outside the whole class from n = 4 on, the constant-size Kingman
  controls, and a Monte Carlo cross-check of the engine against msprime.
- **Folded two-site spectra** (`popcorn.foldgate`, engine `foldgate_engine.py`):
  the minor-allele folding map Phi of a two-site frequency spectrum (what an
  unpolarized data set observes: cell (i,j) pooled with (i,n-j), (n-i,j),
  (n-i,n-j); classes C_a = {a, n-a}) on stored symmetric coordinates, its exact
  transpose (pullback) and exact ranks, plus exact deciders that place a folded
  target spectrum against the folded image of the variable-population-size
  Kingman class: an exact convex combination of folded single-atom components
  (ANNIHILATED: inside the folded class closure, nothing separates it), a folded
  functional whose pullback kernel is Bernstein-certified nonnegative on the unit
  square with an exact negative margin (SURVIVES_OUT_FULLCLASS), an exact conic
  obstruction that rules out every such functional (NO_LINEAR_WITNESS), and the
  weaker diagonal certificate (outside the hull of all single-atom poolings).
  Floats (HiGHS LPs) only propose; Fractions decide. Ships the exact Kingman-class
  polynomials (single-atom second moments and the two-time block-count kernel)
  and a sixteen-target reference table (Beta, Dirac and Kingman coalescents) at
  n = 4, 5, 6.
- **Two-site spectrum from genotypes** (`popcorn.twosite`; engines
  `twosite_projection.py`, EXACT, and `twosite_spectrum.py`, data
  utility): for a pair of biallelic sites typed in the same N diploids,
  `project_grid(cells, N, m)` is the exact probability table of the two
  alt-allele counts (i, j), i, j = 0..2m, in a uniformly random subset of
  m individuals — the expectation over ALL C(N, m) subsets, by an integer
  dynamic program over the nine genotype cells (multivariate
  hypergeometric) with one exact division, as Fractions summing to
  exactly 1; the two-site analogue of hypergeometric SFS projection, over
  individuals so both sites share the subsample, no phase needed. Around
  it: explicit-enumeration reference route (`brute_grid`, `brute_check`),
  derived-allele orientation, unfolded segregating block, minor-allele
  fold, exact accumulation and a memo by distinct table. The estimator
  `estimate()` / `twosite_spectrum.py` reads a VCF (or a `pos ref alt gts`
  dosage TSV), an optional `pos ref alt anc conf` polarization table (the
  table `popcorn.ancestral` writes) and an optional bedGraph genetic map,
  scans site pairs within `max_bp`, classifies them (candidate
  multi-nucleotide pairs, map gaps, hotspot spans, out of range, signal),
  and accumulates the exact projected pair spectra per distance bin (cM
  with a map, bp without) at each requested n — unfolded over pairs whose
  two sites orient, folded over all — with per-bin normalized spectra as
  exact rationals and delete-one-block jackknife intervals. Standard
  library only.
- **Planted truths for two-site / SFS estimators** (`popcorn.planted`, engine
  `synth_truths.py`): (i) the closed-form (float) expected branch-length
  spectrum E[L_i] under any piecewise-constant Kingman history (Tavaré's
  ancestral-process transition function integrated epoch by epoch; double
  precision, measured accurate to 1e-16 at n <= 6 and 1e-12 at n = 20
  against a 60-digit route) and an epoch-size fit that builds an SFS-matched
  pair of histories; (ii) seeded msprime generators for named truth classes
  (TK1 three-epoch bottleneck, TK2 its four-epoch SFS-matched partner, Beta
  coalescents alpha = 1.3 / 1.7, Dirac coalescents psi = 0.05 / 0.2, or
  classes declared in a pins file) that plant q_ij = E[L_i L_j]/E[L_tot^2]
  and xi_i = E[L_i]/E[L_tot] at n = 4, 6 with jackknife Monte Carlo errors,
  and write two-site-ready synthetic data per arm (recombination ON/OFF, one
  gene-conversion arm): a VCF-like unphased genotype table plus the true
  `pos ref alt anc conf` polarization table, with mu calibrated per arm to a
  fixed segregating-site density and r fixed by a map-average r/mu. One
  pins file fixes the base seed every stream derives from, so any object
  regenerates alone and, under the recorded library versions, bit for bit.
  SIMULATION / DATA-UTILITY register.
- **Two-window coalescent simulation harness** (`popcorn.twowindow`, engine
  `twowindow_sim.py`, msprime): simulate a pair of 1-bp windows a distance
  d apart on one sequence (`run_two_window`), or two loci at population-scaled
  recombination rho = 4 Ne m (`run_two_locus`), under Kingman,
  Beta-/Dirac-coalescent, selective-sweep (`sweep_model`) or gene-conversion
  (geometric tract) models, and return the ordered branch-mode second-moment
  matrix M[i-1][j-1] = sum afs_i(A) afs_j(B) (i, j = 1..n-1), its pooled pair
  vector (diagonal doubled, unit sum; 28 entries at the default 4 diploid
  samples, n = 8) with delete-one-block jackknife standard errors, and the
  measured decorrelation fraction n_afs_distinct / replicates (replicates
  whose two windows carry different genealogies); `run_probe_track` does many
  same-window probes along one long recombining sequence. SIMULATION register
  (Monte Carlo with stated errors; never exact or certified).

USE-WHEN:

- Certified selection-SFS entries/vectors/exact gradients for PRF DFE
  likelihoods; auditing a float-pipeline popgen fit; biobank-scale n where
  float pipelines are structurally unavailable.
- Certified DFE-mixing kernel; dominance h != 1/2.
- Exact hull/positivity/LP certificates, and region certificates — "the data
  point is provably >= this far from the WHOLE model region", not just from a
  grid of samples.
- Exact neutral expected SFS under any Lambda-coalescent for n up to a couple
  of hundred (constant population size); an exact reference against which to
  validate a float or Monte Carlo implementation (e.g. msprime's Beta and
  Dirac models); exact sign decisions on linear functionals of the spectrum.
- Two-locus: you need the recombining two-site spectrum (or its rho = 0 /
  rho -> infinity limits) at small n as a smooth function of rho and of a
  piecewise-constant size history — as the model side of a two-site
  composite likelihood, to calibrate LD-decay summaries, or as the
  deterministic comparand for a coalescent simulator. n <= 7 is sub-second;
  n = 8 at rho > 0 costs ~13 s and ~0.8 GB per evaluation (states grow ~2.9x
  per added sample).
- Non-equilibrium selected spectra at sample sizes where moment closures
  degrade (n ~ 10^3), e.g. a DFE likelihood under a bottleneck/growth
  history; stationary questions go to the certified engine (`popcorn.sfs`,
  `certsfs.py`) instead.
- Polarizing a VCF-derived sites list to ancestral/derived before building
  an unfolded SFS or any derived-allele statistic; measuring the
  ancestral-mismatch and low-confidence rates of a call set against the EPO
  ancestor; producing the `pos ref alt anc conf` table other tools consume.
- You need the exact expected 2-SFS (or SFS covariances) of linked sites
  under a Beta-, Dirac- or general Lambda-coalescent at small n (n <= ~20),
  e.g. as the reference for a simulator or a float implementation, or to ask
  whether a multiple-merger signal in the linked 2-SFS can be mimicked by ANY
  Kingman demography (population-size history) or mixture of demographies —
  with an exact certificate either way rather than a fit over a family of
  size histories.
- You have (or model) a two-site frequency spectrum WITHOUT ancestral states and
  need the folded object and its exact linear algebra (fold, transpose, rank,
  kernel dimension); you want to know, exactly, whether folding destroys a
  model-class separation that exists for the polarized spectrum (e.g. whether a
  folded multiple-merger two-site spectrum can still be told apart from every
  variable-size Kingman history by a linear certificate) at small n.
- You need the observed two-site (joint) frequency spectrum of a diploid
  panel at a fixed small n as a function of genetic or physical distance —
  the data side of a two-locus composite likelihood, an LD-decay summary
  that does not pick a subsample, or the comparand of `popcorn.twolocus` /
  a coalescent simulator; or the exact pair-projection operator itself
  inside your own scan. The exact work scales with the number of DISTINCT
  3x3 tables, not pairs: measured worst case (random genotypes, nearly
  every pair a new table, n = 4 and 6, one core) about 3e4 pairs/s at
  N = 12 and 3e3 pairs/s at N = 50-100; panels dominated by rare variants
  repeat tables and run faster.
- A certified (interval) value of a non-equilibrium selected spectrum at small
  n, e.g. to price the error of a float transient solver (a diffusion grid or
  a moment closure) at chosen cells of (S, history), or to carry a rigorous
  error bar on a transient SFS entry into a likelihood comparison; two cells
  at different truncation M give an independent control. n = 20 at prec = 256
  is the measured regime (60+ digit enclosures, ~25 s per M = 200 cell).
- Validating a two-site or SFS pipeline end to end: run your estimator on a
  planted arm's tables and require the planted q_ij / xi_i back within the
  quoted Monte Carlo error; checking a coalescent simulator's Beta / Dirac /
  variable-size conventions against exact spectra (the battery does exactly
  this against `popcorn.lambda_coalescent` and `popcorn.twolocus`); getting
  the closed-form (float) expected SFS of a multi-epoch Kingman history at
  small n, or a second history with the same SFS.
- Checking an exact two-site prediction by Monte Carlo (at rho = 0 or with no
  process between the windows the two-window harness estimates exactly the
  linked E[L_i L_j] of `popcorn.twosfs`; at rho > 0 the E[T_i^A T_j^B] of
  `popcorn.twolocus`); measuring how gene conversion and crossover decorrelate
  neighboring sites with distance; putting sweeps, multiple-merger
  coalescents or a demography behind a two-site statistic where no exact
  engine exists. About 0.2 ms per two-locus replicate at n = 8; at 16,000
  replicates (a few seconds) the 28 cells carry 1-5 % relative error.

NOT-FOR: certified demography or demographic inference. The certified
SFS/DFE/dominance stack is equilibrium, constant-N only; piecewise-constant
N(t) is certified only for the transient selected SFS at small n through
`popcorn.enclosure` (model-conditional enclosures, S <= 0, genic selection);
`popcorn.transient` (single-locus selected SFS at large n) and
`popcorn.twolocus` (neutral two-locus branch-length moments) remain
FLOAT-VALIDATED; still no S > 0, dominance or linkage in the enclosure; and
`popcorn.lambda_coalescent` is neutral and constant-N. No migration or multi-population spectra, no continuous size
change, no fitting driver. `popcorn.twosfs` and `popcorn.foldgate` certify,
exactly, whether a given model's linked two-site spectrum at small n lies
inside or outside the WHOLE variable-size Kingman class (polarized or
minor-allele-folded); they do not fit or infer a size history.
`popcorn.twosite` is model-free bookkeeping of observed genotypes (exact
projection arithmetic, float jackknife), not a likelihood. `popcorn.planted`
and `popcorn.twowindow` are msprime simulation harnesses (planted truths /
two-window Monte Carlo with stated errors): they plant and measure, they do
not score, fit or certify. The 1F1 closed
form for the fixed-S SFS is the classical PRF result (Sawyer & Hartl 1992;
Bustamante et al. 2001), not ours. The
package's built-in self-checks certify internal consistency; a production
claim should additionally be gated against an oracle that shares no code with
the leg being checked.

INVOKE:

```sh
python3 tools/popcorn/certsfs.py selftest --dps 40      # pins + reference values
python3 tools/popcorn/certsfs.py entry  20 3 -100 --dps 40
python3 tools/popcorn/certsfs.py check  100 50 -1000 --dps 40   # two routes + agreement
python3 tools/popcorn/certsfs.py vector 100 -1000 --dps 40 --grad
python3 tools/popcorn/lambda_exact.py                       # exact identity checks + n=20 reference grid (prints only)
python3 tools/popcorn/twolocus_engine.py 6 1.0               # two-locus demo: states, E[T_i], E[T_i^A T_j^B] at n=6, rho=1
python3 tools/popcorn/transient_sfs_engine.py                # transient engine smoke (~1 s)
python3 tools/popcorn/epo_join.py --fasta ancestor_21.fa --sites sites.tsv --out anc_chr21.tsv.gz --chrom chr21 --summary counts.json
python3 tools/popcorn/twosfs_engine.py --n 5 --family beta --param 3/2        # identity checks + exact E[L_i], E[L_i L_j], q
python3 tools/popcorn/twosfs_certificates.py --n 4 --family beta --param 1/2  # class verdict: OUT_FULLCLASS, exact margin, witness
python3 tools/popcorn/twosfs_certificates.py --n 4 --family kingman           # control: OUT_ATOMPOOLINGS only, never OUT_FULLCLASS
python3 tools/popcorn/foldgate_engine.py            # re-derives all 16 reference decisions from the shipped ingredients (a few seconds; prints only)
python3 tools/popcorn/foldgate_engine.py --n 4 --out folded_n4.json
python3 tools/popcorn/twosite_projection.py                     # DP vs enumeration on the shipped 12 x 24 panel, n = 4,6,8
python3 tools/popcorn/twosite_spectrum.py --geno calls.chr21.vcf.gz --anc anc_chr21.tsv.gz \
        --map genetic_map.bedGraph.gz --chrom chr21 --n 4,6 --out spectrum_chr21.json
python3 tools/popcorn/twosite_spectrum.py --geno panel.tsv --n 8 --bins 1,100,1000,10000,100000 --out s.json   # bp bins, folded only
python3 tools/popcorn/twosite_spectrum.py --brute-check --geno small_panel.tsv --n 4,6,8                       # <= 30 individuals
python3 tools/popcorn/transient_enclosure_engine.py cell --M 200 --S -5 --json cell.json   # one certified cell + float route (+ moments if installed), ~30 s
python3 tools/popcorn/transient_enclosure_engine.py xcheck                              # independent mpmath route vs the ball route (~8 s)
python3 tools/popcorn/synth_truths.py check                                   # analytic checks vs the reference set (no msprime)
python3 tools/popcorn/synth_truths.py pin --out sim/ --seed 7 --map genetic_map.bedGraph --chrom chr7   # or --r-per-bp 1.2e-8; --extra my_classes.json
python3 tools/popcorn/synth_truths.py truth --cls B13 --out sim/ [--pins sim/pins.json] [--reps 20000]  # default pins = the shipped reference
python3 tools/popcorn/synth_truths.py arm --arm TK1_rON_GC --out sim/ [--pins sim/pins.json]           # <ARM>.gt.tsv.gz, <ARM>.anc.tsv.gz, arm_<ARM>.json
python3 tools/popcorn/synth_truths.py selftest                                # reduced-scale end to end in a temp dir (<1 s, needs msprime)
python3 tools/popcorn/twowindow_sim.py 1.0 2000                # demo: two-locus arm at rho = 0 (vacuity check) and rho = 1, pooled vectors + jackknife SEs
python3 -c "import sys; sys.path.insert(0,'tools'); import popcorn.twowindow as W; r = W.run_two_locus(0.0, 2000, seed=1, n_blocks=10); print(r['n_afs_distinct'], r['pooled'][:3], r['pooled_se'][:3])"
python3 -c "import sys; sys.path.insert(0,'tools'); import popcorn.twowindow as W; r = W.run_two_window(1000, 2000, seed='gc', gc_rate=2e-8, gc_tract=300, n_blocks=10); print(r['params']['margin_bp'], r['afs_distinct_fraction'])"
```

```python
import popcorn                       # with tools/ on PYTHONPATH
popcorn.verify()                     # sha pins; [] iff clean — treat any hit as fatal
alpha, beta = popcorn.sfs.solve_M_exact(50, -10)         # exact Fractions
K = popcorn.dfe.CertifiedKernel(20, dps=40); K.build(-1); K.build(1)
Ev, selfcons_d, tail_rel = K.mix({'b': 0.4, 'Sd': -1000, 'pb': 0.02, 'Sb': 10})
E, sc = popcorn.dominance.E_dominant(20, 3, -100.0, 0.3, dps=50)
C = popcorn.certificates             # lazy: scipy only if cert_lp touched
status, cert = C.exact_lp.membership(points, q)          # exact Farkas witness
R = C.region                         # region certificates (see its docstring)
from fractions import Fraction
L = popcorn.lambda_coalescent        # exact, standard library
h = L.expected_lengths(50, L.beta_rate(Fraction(3, 2)))[50]   # {i: E[L_i]} exact Fractions
x = L.xi_hat(50, L.dirac_rate(Fraction(1, 10)))               # normalized expected SFS, sums to 1 exactly
W2 = popcorn.twolocus                # numpy + scipy; float-validated
M, a, b = W2.moments(6, 0.9)         # E[T_i^A T_j^B], E[T_i] at n=6, rho=0.9; W2.pool_pairs(M)
T = popcorn.transient                # numpy + scipy; float-validated
eng = T.TransientSFSEngine(n=1000)   # build once per n
E, meta = eng.expected_sfs(-5.0, [(0.33, 0.47), (3.38, 0.017)]); assert meta['negative_entries'] == 0
from popcorn.ancestral import join_sites, polarization_rates
counts = join_sites('ancestor_21.fa', 'sites.tsv', 'anc_chr21.tsv.gz', chrom='chr21'); rates = polarization_rates(counts)
T2 = popcorn.twosfs                  # engine + hull: standard library; certificates lazy (sympy)
u, v = T2.lambda_moments(5, T2.beta_rate(Fraction(3, 2)))    # exact E[L_i], E[L_i L_j] (i <= j)
q, EL2 = T2.normalize_matrix(v, 5)              # exact normalized linked 2-SFS, E[L_tot^2]
m1, m2 = T2.component_moments(5, Fraction(1, 2))              # Kingman block-count moments at level t
r = T2.class_verdict(4, T2.dirac_rate(Fraction(1, 4)))        # sympy + numpy/scipy: 'OUT_FULLCLASS', exact margin, witness
c = T2.certify_witness(4, T2.dirac_rate(Fraction(1, 4)), r['witness_W'])   # re-derive the certificate exactly, no LP
FG = popcorn.foldgate                 # folding map + loaders: standard library; hunts need numpy+scipy+sympy
Phi = FG.fold_coeffs(6)               # {(a,b): {(i,j): weight}}; FG.rank_report(6) -> rank 6, kernel 9
phi = FG.apply_fold(5, q)             # q: {(i,j): Fraction}, i <= j  ->  folded {(a,b): Fraction}
W = FG.pullback(5, vf)                # vf: {(a,b): Fraction}; Phi^T vf: <vf, Phi q> == <W, q> exactly
ing = FG.load_ingredients(4); tg = FG.tgrid_of()            # shipped Kingman-class polynomials, 83-node grid
kern, (s, u) = FG.kernel_sympy(ing); kernF = FG.fold_exprs(4, kern)
row = FG.load_reference_rows()['rows'][0]                   # n=4 Beta(1/2): exact M = E[L_i L_j]
obj = FG.make_object(4, row['M'])                           # {'n','M','tot','q','phi'}
rec = FG.decide(4, obj, FG.folded_components(4, ing, tg), tg, kernF, kern, s, u)   # -> 'NO_LINEAR_WITNESS' + exact conic certificate
TW = popcorn.twosite                 # standard library; exact core + data utility
G = TW.project_grid(TW.nine_cell(d1, d2), len(d1), 3)   # (i, j) table at n = 6, Fractions, sums to 1 exactly
assert G == TW.brute_grid(d1, d2, len(d1), 3)           # explicit enumeration route
U = TW.seg_matrix(TW.flip_grid(G, alt1_is_ancestral, alt2_is_ancestral)); F = TW.fold_matrix(G)
doc = TW.estimate('calls.chr21.vcf.gz', anc='anc_chr21.tsv.gz', gmap='genetic_map.bedGraph.gz', chrom='chr21', n=(4, 6))
doc['results']['n4']['unfolded']['bins']['3']['qhat_exact']    # normalized 3x3 spectrum of cM bin 3 as 'p/q' strings
EN = popcorn.enclosure               # lazy; python-flint on first use
rec = EN.certified_cell(200, -5, 20, 256, verbose=False)   # M, S, n, prec bits; ~25 s
rec["gates"]; rec["enclosure"][0]    # in-cell certified checks; {'i': 1, 'mid30': '[0.4994... +/- ...]', 'rad': ..., 'rel_width': ...}
balls = rec["entry_balls"]           # arb balls E_n(i), i = 1..n-1 (ctx.prec is left at 256)
hist = EN.epochs_fmpq([("1/10", "1/20"), (2, "1/40")])     # exact rationals; past -> present after the rho = 1 equilibrium
rec2 = EN.certified_cell(80, -20, 20, 256, epochs=hist, verbose=False)
EN.float_cell(80, -20, 20, epochs=hist)["sfs"]             # FLOAT register, implementation check only
ref = EN.load_reference_cell(500, -5); EN.entry_ball(ref["enclosure"][0])   # shipped cells as balls
PT = popcorn.planted                 # numpy; scipy for the fit; msprime/tskit only when generating
EL = PT.expected_branch_sfs([(0.0, 2e4), (500.0, 3e3), (1500.0, 1.6e4)], 6)   # closed-form E[L_i] in generations (ploidy-1 units, float)
sizes, rep = PT.match_epoch_sizes(PT.TK1_EPOCHS, [0.0, 150.0, 700.0, 2500.0])  # SFS-matched partner history
pins = PT.load_pins()                # shipped reference pins (or PT.make_pins("sim/", seed=7, r_per_bp=1.2e-8))
rec = PT.run_truth("B13", pins, 20000, "sim/")                # planted q_ij / xi with jackknife SE -> sim/truth_B13.json
arm = PT.run_arm("TK1_rON_GC", pins, "sim/", "sim/")          # tables + sim/arm_TK1_rON_GC.json (seeds, mu, r, counts, payload hashes)
ref = PT.load_reference_truth("D020")["mc"]["n4"]["q_ij"]     # the shipped planted value
import msprime                       # optional; only the simulation calls below need it
TW2 = popcorn.twowindow              # numpy; msprime + tskit on first simulation; SIMULATION register
r = TW2.run_two_locus(0.0, 16000, seed=1, n_blocks=40)        # rho = 0: r['n_afs_distinct'] must be 0
keys, exact = TW2.exact_linked_pooled(8)                       # exact Kingman comparand from popcorn.twosfs (Fractions); pass lam=popcorn.twosfs.beta_rate(Fraction(3,2)) for Beta
z = (r['pooled'] - [float(x) for x in exact]) / r['pooled_se']
g = TW2.run_two_window(10_000, 20000, seed='arm1', gc_rate=1e-8, gc_tract=300, n_blocks=20)   # margin law -> 3000 bp
sw = TW2.run_two_window(1000, 5000, seed=7, recomb_rate=1e-8, model=[TW2.sweep_model(500, 0.05), msprime.StandardCoalescent()])
bt = TW2.run_two_locus(0.0, 16000, seed=3, samples=8, ploidy=1, model=msprime.BetaCoalescent(alpha=1.5), n_blocks=40)   # ploidy 1 for Lambda-coalescents
```

INPUTS: S = 4*Ne*s (S>0 advantageous), unfolded derived-allele SFS, theta=1
normalization. polyDFE/fastDFE S == this S; dadi/fitdadi gamma == S/2.
Dominance conventions: fitnesses 1 : 1+2sh : 1+2s. entry/check/selftest need
mpmath only; vector needs python-flint above n=2000. Certificate polynomials:
dict {exponent: rational} on [0,1] (positivity) or {(es, eu): Fraction} on
[0,1]^2 (region). Lambda-coalescent parameters as Fractions or 'p/q' strings;
two-locus time in 2 N_ref generations, eta = N_ref/N, rho = 4 N_ref R;
transient epochs [(nu, T), ...] past -> present, genic S only; EPO FASTA one
chromosome per file, sites as `chrom pos ref alt` or `pos ref alt`. Two-site
moment matrices as dicts keyed (i, j), i <= j, on stored symmetric
coordinates with pooled-pair mass sum (2-[i=j]) q_ij = 1; folded objects on
1 <= a <= b <= floor(n/2); levels t = e^{-g} and witnesses as rationals /
'p/q' strings (msprime_dirac takes 'psi,c'); Kingman-class ingredients
shipped for n = 4, 5, 6 (other n: feed E[L_i L_j], the diagonal polynomials
and the two-time kernel from `popcorn.twosfs` in the same dict format).
Two-site genotypes as VCF (biallelic SNVs, no missing GT at a retained site;
dosage = count of `1` alleles) or `pos ref alt gts` TSV (gts = one 0/1/2
character per individual); polarization table `pos ref alt anc conf` as
written by `epo_join.py`; genetic map as bedGraph `chrom start end cM/Mb`;
n even (m = n/2 individuals <= N). Enclosure engine: S a non-positive
integer (S = 4 N_ref s, the package convention; dadi/moments gamma = S/2),
epochs [(rho, T), ...] past -> present with rho = N/N_ref and T in 2 N_ref
generations as exact rationals (`flint.fmpq`; `popcorn.enclosure.epochs_fmpq`
converts ints, Fractions and 'p/q' strings and refuses binary floats),
applied after an ancestral equilibrium at rho = 1; theta = 1; truncation
M >= n; prec in bits (rule of thumb prec >= 64 + 3.33 (0.30 n + D) for D
digits). Planted truths: ploidy-1 units — epochs [(t_start_gen,
N_chromosomes)], pair rate 1/N per generation (a diploid N_e = 1e4 is
N = 2e4); q_ij / xi are unit-free; pins.json fixes the base seed, r/mu and
the fitted TK2 sizes. Two-window harness: msprime's units — `Ne` is
msprime's population_size (diploid effective size at ploidy 2), branch
lengths in generations (pooled vectors are scale free), `samples` =
individuals of the given `ploidy` (n = samples x ploidy haploid sequences,
default 4 x 2 = 8), crossover and gene-conversion initiation rates per bp
per generation, `rho_scaled` = 4 Ne m with m the per-generation crossover
probability between the two sites (the rho of `popcorn.twolocus`); `seed`
an int or a str (hashed by blake2b to a 32-bit root seed);
`sweep_model(position, s)` takes msprime's s (fitnesses 1, 1 + s/2, 1 + s,
so S = 2 Ne s in this package's convention).

OUTPUTS: CLI prints values at the requested dps (check mode: two routes +
agreement digits; below dps-5 agreement is a named CertSFSGateError, rc=3 —
never a silent fallback). Library returns exact Fractions / certified balls /
kernel objects / exact certificates that are independently re-substituted
before being returned.

REQUIREMENTS: python3 + mpmath (core). Optional, each unlocking legs that
otherwise SKIP BY NAME: python-flint (ball route, kernel, cone/positivity
route R; `popcorn.enclosure` — legs enclosure_cell, enclosure_stationary,
enclosure_reference and the `--full` leg enclosure_full SKIP BY NAME without
it), gmpy2 (LP certificate family), numpy + scipy (float proposers;
`popcorn.twolocus`, `popcorn.transient`, the LP proposers of
`popcorn.twosfs.class_verdict` and the `popcorn.foldgate` hunts, the
enclosure's `float_cell` sub-check, `popcorn.planted`'s epoch-size fit —
legs twolocus_*, transient_*, twosfs_verdict, foldgate_deciders and
planted_vs_twolocus SKIP BY NAME without them), sympy (Sturm root counting,
positivity route S; the `popcorn.twosfs` kernel + certificate layer and the
`popcorn.foldgate` hunts — legs twosfs_kernel, twosfs_class_certs,
twosfs_verdict and foldgate_deciders SKIP BY NAME without it), msprime >= 1.0
(brings tskit; `popcorn.planted` generators / regeneration and every
`popcorn.twowindow` simulation — legs planted_sim, twowindow_linked,
twowindow_recomb and the `--full` leg planted_regen SKIP BY NAME without it;
`import popcorn.planted` / `import popcorn.twowindow` and their analytic,
geometry, pooling and jackknife helpers do not need it), moments
(moments-popgen; only the enclosure CLI's optional FIELD-FLOAT cross-check,
imported inside `moments_cell`, which prints "SKIP moments cross-check"
without it; no leg needs it). `popcorn.lambda_coalescent`,
`popcorn.ancestral`, `popcorn.twosite`, the `popcorn.twosfs` engine + hull,
the `popcorn.foldgate` folding map / exact verifiers / loaders and
`popcorn.planted`'s planted_vs_lambda leg need the standard library only;
`popcorn.twowindow`'s helpers and `popcorn.planted`'s closed form need numpy.
The two-window harness's exact comparand `exact_linked_pooled` uses
`popcorn.twosfs` (the n = 8 Kingman values also ship under
`reference/twowindow/`); its rho > 0 battery cross-check uses
`popcorn.twolocus` (scipy) and is skipped with a printed note without it.

GATES: every pinned file byte-identical (`popcorn.verify()`, and the CLI
re-verifies before any engine import); certsfs check agreement >= dps-5 or
rc=3; battery default tier = the fast legs below; `--full` adds the held-out
two-route gate (gate_phase1.py: >= 40 matched digits at held-out points to
n = 10007, gradient gate via the parameter-shift identity), the dominance
oracle's full selftest, the n=20 kernel selftest, the M = 200 certified
transient cell reproduced against the shipped one, and the CLI regeneration
of two planted truths plus the gene-conversion arm.

Battery: `python3 tools/popcorn/selftest.py` → one line per leg + OVERALL
PASS, rc=0. Measured ~35 s (43 legs) with all optional engines present; legs
missing an engine SKIP BY NAME stating what to install (skips do not fail the
battery). Includes planted negative controls (exterior LP points, negative
polynomials, the region selftest's MUST-FAIL plants, planted MUST-FAIL
polynomials and perturbed convex identities in the two-site legs, and the
constant-size Kingman control that must never certify out of its own class).
Two-site legs (a few seconds together): twosfs_identities, twosfs_hull_ref,
twosfs_kernel, twosfs_class_certs, twosfs_verdict, twosfs_montecarlo,
foldgate_map, foldgate_annihilated, foldgate_obstructions, foldgate_deciders,
twosite_projection, twosite_spectrum, twosite_cli. Enclosure legs (about
4.5 s together): enclosure_cell (live M = 40 and M = 60 cells at S = -5,
n = 20: both in-cell certified checks true, relative width 5e-15 -> 9e-37 as
M grows, the two enclosures overlap entrywise, the float route within 1e-8
of the certified midpoint, `moments` never imported), enclosure_stationary
(97 certified stationary entries and the projected equilibrium moments agree
with `popcorn.sfs`'s exact route to >= 30 digits, measured 71-75; S = 0 gives
1/i), enclosure_reference (the six shipped cells parse, checks true, M = 200
vs M = 500 identical to 30 digits and overlapping, a live M = 60 enclosure
overlaps the shipped S = -5 cells, and the independent mpmath route
reproduces the shipped worst_rel_dev of 1.55e-32 bit for bit — that number
is the half LOWER/UPPER bracket at M = 40, j = 20; the mpmath-vs-LOWER
agreement is 2e-50). Planted legs (about 2 s): planted_analytic (Tavaré
identities, constant-N 2N/i, epoch-split and rescaling laws, the shipped TK2
fit and E[L_i] reproduced, planted xi of TK1/TK2 vs the closed form
|z| <= 2.7), planted_vs_lambda (planted xi of the Beta/Dirac classes vs
`popcorn.lambda_coalescent` exact spectra |z| <= 2.3, one-atom Dirac control
rejected at |z| 41), planted_vs_twolocus (planted q_ij vs `popcorn.twolocus`
at rho = 0 |z| <= 2.4; TK2-vs-TK1 distance recomputed), planted_sim (seeded
reduced end-to-end regeneration, bit-identical under the recorded versions,
Monte Carlo agreement otherwise, mode printed). Two-window legs (about 5.5 s):
twowindow_laws (geometry, margin, seed, pooling and jackknife laws on
literals; live `popcorn.twosfs` == the shipped n = 8 reference),
twowindow_linked (rho = 0, n = 8, 16000 replicates vs exact Kingman
max|z| <= 4 with the unlinked limit > 20 SE away; n_afs_distinct == 0;
gene conversion decorrelates), twowindow_recomb (rho = 1, n = 6 vs
`popcorn.twolocus` max|z| <= 4.5, rho/2 and 2 rho rejected). `--full` is
minutes-scale (enclosure_full ~25 s, planted_regen ~10 s on top of the three
heavy gates). The battery writes nothing into the package tree (scratch goes
to a temp dir).

FOOTGUNS:

- The shipped engine bytes are sha-pinned; `popcorn.verify()` and the CLI
  fail closed on drift. Regenerate pins only when deliberately changing
  shipped bytes (`python3 -c "import popcorn._pins as p; p.regen()"`), and
  update `certsfs.PINNED_SHA256` in the same change.
- Exact assembly: alpha/beta are exact for a SPECIFIC rational S — assembling
  with any other nearby S (e.g. the binary float of a decimal) is amplified
  by the full cancellation factor (measured: thousands of digits lost).
  Thread ONE exact rational S through solve and assembly.
- Never substitute mp.quad for certquad: mp.quad's error estimate is absolute
  and it silently underconverges on exponential-ramp panels — the defect
  certquad exists to cure. A wrong lnf estimator costs panels, not
  correctness, BUT an lnf that reports the floor over live support silently
  drops mass — sanity-check lnf against f.
- Convert parameters to mpf INSIDE the workdps block: converting outside
  rounds them at ambient precision (measured: a uniform 16-digit ceiling).
- Degree-pair self-consistency between correlated legs is not proven digits;
  the h=1/2 collapse gate and the two-route CLI check are the independent
  anchors shipped here.
- `bern_nonneg` (region) is subdivision-only: a polynomial that TOUCHES zero
  on the box needs touch-root deflation — positivity route B.
- `cone_lp_b`'s one-RREF bit-sorted basis completion is load-bearing: the
  greedy per-column variant of the same simplex measured ~40x slower. Do not
  "simplify" it back.
- q-series dominance points at extreme h and large |S| cost the full computed
  precision pad (deliberate; can run minutes) — prefer the oracle there.
- Consumers fitting likelihoods built on these outputs: optimize on the
  deviance scale C*(d - log1p(d)), not raw logL — when the objective's
  dynamic range swamps the curvature, L-BFGS-B/FD-Hessian breaks on raw logL
  and the deviance form recovered ~1e-11 relative accuracy in testing.
- Register EXACT but pure-Python Fraction arithmetic: cost is O(n^3) rational
  operations whose operand size grows with n and with the family (Beta(3/2):
  ~2 s at n=100, ~17 s at n=200; Dirac(1/2): ~3 min at n=200; Kingman is
  O(n^2)). Pass parameters as Fractions or strings ('3/2'), never as binary
  floats, or the rates are exact for the wrong rational. Constant population
  size and expected values only. The shipped n=20 Kingman-class functional's
  nonnegativity over the variable-size Kingman class is carried by the
  reference file as given, not re-derived by this module; the module only
  evaluates it exactly.
- Two-locus register is FLOAT-VALIDATED, not certified — double-precision
  sparse linear algebra with no enclosure; validated at rho = 0 against the
  exact-rational Kingman E[L_i L_j] (n = 8: ~2e-16 on the normalized pair
  vector, ~1e-13 relative on raw entries), at n = 2 against the closed-form
  two-locus covariance (rho + 18)/(rho^2 + 13 rho + 18) for all rho, and by
  self-consistency across rho (Kingman marginals 2/i at any rho, A/B
  symmetry, epoch-splitting invariance, time-rescaling covariance, ~1e-15).
  Never mix its output into an exact/certified table without the label.
  Units: time in 2 N_ref generations, eta = N_ref/N, rho = 4 N_ref R;
  `epochs_from_Nt` divides generations by 2 N_ref (default N_ref = 1e4).
  `pool_pairs` doubles the diagonal (it reads M + M^T on i <= j): it is a
  pooling convention, not the class law of an unordered pair of sites —
  compare only against vectors pooled the same way.
- `popcorn.transient` is FLOAT-VALIDATED, not certified: double precision
  with a measured trust radius and nothing more. Measured (data shipped
  under `reference/transient/`): at n = 1000 on the default grid, the
  projected equilibrium and the equilibrium held T = 0.5 through the
  integrator agree with the certified two-route stationary references
  (shipped to 20 significant digits) to <= 9.4e-5 max-relative in every
  frequency band for S in {0, -1, -5, -20, -100}, quadrature-dominated and
  concentrated in the log-spaced bands; transient self-convergence on a
  two-epoch bottleneck-then-growth history, entries E_i > 1e-5: dt
  refinement 4e-4 -> 1e-4 <= 3.9e-4 for |S| <= 20, 4.0e-3 at S = -100,
  1.1e-2 at S = -200 (grid refinement <= 5.2e-3 throughout); the
  exponentially suppressed tail (E_i < 1e-5) moves by up to 9e-2 at
  S = -200. Other n, |S| > 200, custom grids that break seam spacing
  continuity (an abrupt jump costs about 2.5e-4), larger dt0 and sharper
  histories are unmeasured. `expected_sfs` never raises on negative tail
  entries — check `meta["negative_entries"]`. Selection is genic (h = 1/2)
  only; for dominance use `popcorn.dominance` (stationary). Never quote a
  number from this module beside certified ones without labeling each.
- Register is DATA-UTILITY: integer counters and string columns, no
  estimator and nothing certified; the summary rates are plain ratios of the
  counters. One chromosome per FASTA: only the first record is read (further
  records are ignored and flagged as `fasta_extra_sequences`), so always
  pass `--chrom` with a multi-chromosome sites table or rows of other
  chromosomes are joined against the wrong sequence. Positions are 1-based;
  a position past the sequence end yields anc `.` and counts in
  `out_of_fasta_range`, it is not an error. Orientation case-folds the
  ancestral base (a low-confidence `t` orients like `T`); filter on `conf`
  downstream if only high-confidence calls are wanted. Compare `.gz`
  outputs by decompressed payload, never by file hash (gzip stores an
  mtime).
- `popcorn.twosfs` is EXACT but two things in it are not certificates: the
  LP inside `class_verdict`/`hunt` only PROPOSES (a different scipy may
  propose a different, equally valid rational witness, or none: UNDECIDED is
  "not found", never "inside"), and `reference/twosfs/msprime_moments.json`
  is a float Monte Carlo CHECK of the engine (worst |z| 1.8 over 72
  statistics), not part of any proof. OUT certificates are sufficient-side:
  K_W >= 0 pointwise certifies over a cone that contains the class, so
  failure to find a witness proves nothing (copositive-type witnesses could
  exist); margins are certificate-specific lower bounds, not distances, and
  only the verdicts should be compared across n. OUT_ATOMPOOLINGS excludes
  poolings of single-atom limit histories on any support but NOT genuinely
  time-varying histories — the constant-size Kingman coalescent is itself
  OUT_ATOMPOOLINGS (and outside every finite-support hull of
  `twosfs_hull.project_hull`), which is exactly why those two registers must
  never be read as "non-Kingman". `bern_1d_nonneg`/`bern_2d_nonneg` are
  three-valued (True / False / None = undecided at max depth); a polynomial
  touching zero inside the box needs `popcorn.certificates.positivity`
  route B. `twosfs_hull.solve_exact` (square, None if singular) and
  `popcorn.certificates.region.solve_exact` (rectangular) are different
  routines with the same name — never star-import both. Scope: neutral,
  one population, no recombination between the two sites, constant size on
  the Lambda side, second moments only; cost grows with the number of
  integer partitions of n (Beta(3/2): n = 20 in ~16-20 s; the set-partition
  and second-Kingman-route cross-checks are for small n only).
- `popcorn.foldgate` is register EXACT, but read the verdicts for what they
  certify: ANNIHILATED is membership in the folded class CLOSURE (single-atom
  poolings are limits of histories); NO_LINEAR_WITNESS is an obstruction to
  every pointwise-nonnegative kernel functional, NOT a membership proof; the
  diagonal certificate excludes poolings of single-atom histories only (at
  n = 6 it excludes even the constant-size Kingman target, which is a smooth
  history, not an atom pooling). `bern_1d_nonneg` / `bern_2d_nonneg` return
  None (inconclusive), not True, for a polynomial that touches zero inside
  the domain; the hunts treat anything but True as "no certificate". LP
  proposals are rationalized with `limit_denominator(10**6)` and then judged
  exactly, so a different scipy build may propose a different (equally
  exact) certificate: compare certificates by re-verifying them, not by
  string equality. Fully linked sites only (no recombination between the
  pair); expected spectra only. Pass reference rationals as Fractions or
  'p/q' strings, never as binary floats.
- `popcorn.twosite`: the projection is EXACT but everything is projected
  over INDIVIDUALS, so n = 2m must be even and every retained site must be
  typed in all N individuals (the VCF reader drops a site with any missing
  GT and counts it; it does not impute). `project_grid` takes the 3x3
  TABLE, `brute_grid` the two dosage VECTORS; `brute_grid`/`brute_check`
  are checkers, O(C(N, m)) per pair, refused above 30 individuals. Objects
  returned by `GridCache` are shared: never mutate them or pass them as
  the first argument of `add_into`/`add_scaled`; the estimator's cache
  grows with the number of distinct tables (peak RSS 1.6 GB at 3e5
  (table, n) entries, measured) — pass `cache_max_entries` on large N,
  values are unchanged. In the estimator the ancestral base is compared
  case-insensitively, so with an EPO-derived table low-confidence
  (lowercase) calls ARE used unless `high_conf_only` is set — decide which
  you want; sites whose table REF/ALT disagree with the genotype file are
  left unoriented (`refalt_mismatch`), not an error. With a map, a site
  inside an uncovered interval has no genetic position and all its pairs
  are excluded (`na_gap_excluded`), and any pair whose span touches an
  interval above `hotspot_cm_per_mb` is set aside (`hotspot_control`) — on
  a fine-scale map with the default 10 cM/Mb this can be a large share of
  pairs; inspect `pair_class_counts`. Pairs at <= `mnv_bp` (default 2 bp)
  never enter the signal. One chromosome per run. The folded register needs
  no polarization; the unfolded one only counts pairs whose two sites both
  orient, so its `n_pairs` is smaller. phi/q_hat are exact; the jackknife
  se/ci are floats over blocks of the pair midpoint and are absent (a note
  instead) when fewer than two blocks carry pairs.
- `popcorn.enclosure` is CERTIFIED-ENCLOSURE but MODEL-CONDITIONAL: the balls
  enclose the moment system of the stated diffusion (genic selection
  h = 1/2, S <= 0, piecewise-constant N(t), theta = 1) through the proved
  two-sided truncation bracket, and nothing else — not the discrete
  Wright-Fisher chain, not linkage, not dominance, and S > 0 is refused
  (ValueError; the bracket lemma needs a Metzler matrix). Validity needs
  only M >= n, tightness needs M >> rho_max |S|: the width self-reports
  truncation (at S = -5, n = 20: 5e-15 relative at M = 40, 9e-37 at M = 60,
  precision floor ~4e-69 from M = 100), so run two M and read the widths
  rather than trusting one. The projection amplifies width by up to
  2^(n-2): raise prec with n (measured only at n = 20, prec = 256).
  `certified_cell` and `xcheck` set the process-wide flint precision
  `ctx.prec` and leave it set. Epoch entries must be exact rationals (fmpq);
  S must be an integer for the float and moments routes and the CLI. The
  engine's `float_cell` (scipy), `xcheck` (mpmath) and `moments_cell`
  (optional `moments`) are implementation checks in FLOAT / FLOAT-HP /
  FIELD-FLOAT register: at the shipped cells the float route sits up to
  3e-8 from the certified midpoint, some 60 orders outside the certified
  width, so the enclosure adjudicates the floats, never the reverse; never
  quote those routes as certified. `xcheck`'s worst_rel_dev is measured
  against the HULL midpoint, which sits half the LOWER/UPPER bracket from
  LOWER, so at small M it reports the half bracket (1.5e-32 at M = 40) and
  is substep-independent; the mpmath route's agreement with the LOWER
  system itself is ~1e-dps (2e-50 at the shipped instance). Cost is ~M^3
  per matrix exponential (two per epoch per bracket): ~25 s per M = 200
  cell, ~6 min at M = 500, single core.
- `popcorn.planted` is SIMULATION register: the planted q_ij / xi_i are Monte
  Carlo expectations known to their quoted jackknife SE (2 x 10^6 replicates
  for TK1/TK2, 20000 for the Beta/Dirac classes in the shipped set), never
  exact; the only closed form is `expected_branch_sfs` and it is double
  precision (validated against the 2N/i law, against `popcorn.twolocus` at
  rho = 0 and against the seeded simulations, |z| <= 2.7; the alternating
  sum loses digits with n: 1e-16 at n <= 6, 1e-12 at n = 20 against a
  60-digit route). Units are ploidy-1: epochs are (generations, chromosomes) with
  pair rate 1/N, so a diploid N_e = 10^4 is N = 2 x 10^4 here; the
  "diploid" individuals of an arm are exchangeable chromosome pairs from a
  haploid coalescent, not a diploid multiple-merger model. mu is
  recalibrated PER ARM to hold the site density, so absolute rates of the
  multiple-merger arms (B13: ~5.6e-7) are synthetic; compare
  population-scaled quantities. Regenerated .gz tables never match by file
  hash (gzip stores an mtime): compare `gz_payload_sha256` or the record's
  `payload_sha256`. Bit-for-bit regeneration of the reference set holds
  under the msprime/tskit/numpy versions recorded in
  `reference/planted/pins.json`; under other versions the battery falls
  back to agreement within Monte Carlo error and says so. A pins file is
  written once and never regenerated (an existing one is returned as is);
  truths are clean (no genotype error, no mis-polarization).
- `popcorn.twowindow` is SIMULATION register — every vector carries Monte
  Carlo error; quote the replicate count, the block count and the jackknife
  SE with any number, and never place one beside an exact or certified value
  without the label. The jackknife SE is itself noisy at small block counts
  (relative noise about 1/sqrt(2(B-1)): 24 % at B = 10, 11 % at B = 40), and
  cells with small expected mass (pairs of high classes) are right-skewed at
  a few thousand replicates: a max-|z| acceptance over all 28 cells needs
  >= 10^4 replicates to behave like a Gaussian tail. `pool_pairs` doubles
  the diagonal (same convention as `popcorn.twolocus.pool_pairs` and the
  i <= j order of `popcorn.twosfs.flatten`); compare only against vectors
  pooled the same way. An arm with no process between the windows MUST give
  n_afs_distinct == 0 — if it does not, the window geometry or the simulator
  call is wrong and nothing else from that run means anything. Batches in
  which msprime raises are counted in `failed` (check it is 0) and a run
  where everything failed raises. Results are reproducible from (seed,
  n_reps, batch, n_blocks) only at a fixed msprime version: changing
  `batch` or `n_blocks` changes the seed stream. In `run_two_locus` the
  crossover mass sits on the rate-map interval [1, 2) because in msprime's
  discrete genome that interval drives the breakpoint between sites 0 and 1
  (mass on [0, 1) does nothing); do not "fix" the map to [m, 0].
  `sweep_model` passes s and the frequencies straight to
  msprime.SweepGenicSelection, whose own limits apply (single population, no
  size change during the sweep; see msprime's documentation on
  start_frequency); the default margin law max(3000, 10 x tract) bp keeps
  gene-conversion tracts that cover a window from being edge-truncated. A
  Beta- or Dirac-coalescent arm that is to be compared with `popcorn.twosfs`
  must run at ploidy = 1 (samples = n): at ploidy p > 1 msprime's
  multiple-merger models merge up to 2p groups at once (a Xi-coalescent) —
  measured at n = 8, Beta(3/2): ploidy 1 matches the exact spectrum at
  max|z| 2.9 over 28 cells, ploidy 2 is 38 SE away.

CREDIT: popcorn implements, exactly and with certificates, the
Poisson-random-field selection-SFS / DFE methodology introduced by Sawyer &
Hartl (1992, Genetics 132:1161); the fixed-S closed form M_i = 1F1(n-i; n;
-S) is the classical PRF result (Sawyer & Hartl 1992; written explicitly by
Bustamante, Wakeley, Sawyer & Hartl 2001, Genetics 159:1779), and selection
with arbitrary dominance h follows Williamson, Fledel-Alon & Bustamante
(2004, Genetics 168:463). The likelihood stack it certifies is the one built
and made standard by dadi (Gutenkunst, Hernandez, Williamson & Bustamante
2009, PLoS Genet. 5:e1000695), fitdadi (Kim, Huber & Lohmueller 2017,
Genetics 206:345), polyDFE (Tataru, Mollion, Glémin & Bataillon 2017,
Genetics 207:1103; the 'model C' DFE family is theirs), fastDFE (Sendrowski
& Bataillon 2024, Mol. Biol. Evol. 41:msae070) and moments (Jouganous, Long,
Ragsdale & Gravel 2017, Genetics 206:1549); we are indebted to those
authors, whose papers and open code define the conventions (S vs gamma =
S/2) this package pins. The recurrence numerics are classical: Olver's
(1967, J. Res. NBS 71B:111) boundary-value formulation and Miller's backward
recurrence (Gautschi 1967, SIAM Rev. 9:24) for the contiguity relations;
Bernstein/de Casteljau subdivision, Farkas certificates and Sturm's theorem
for the exact certificates, arranged under a float-proposes / exact-decides
discipline.

Lambda-coalescents: Pitman (1999), Sagitov (1999); Beta(2-alpha, alpha)
family: Schweinsberg (2003); psi-coalescent: Eldon & Wakeley (2006); expected
SFS recursions under Lambda-coalescents: Birkner, Blath & Eldon (2013);
Kingman closed form E[L_i] = 2/i: Fu (1995); Bolthausen & Sznitman (1998).

Two-locus: the two-locus ancestral process with recombination is due
to Griffiths (1981) and Hudson (1983); the n = 2 covariance formula is
theirs (see McVean 2002, Genetics 162:987); the rho = 0 comparand is the
Kingman branch-length covariance of Fu (1995, Theor. Popul. Biol. 48:172);
two-locus sampling theory under constant and variable N(t) is Hudson (2001,
Genetics 159:1805) and Kamm, Spence, Chan & Song (2016, Genetics 203:1381,
arXiv:1510.06017), of which this chain is the branch-length-moment
counterpart.

Transient: integrates the Kimura (1964) forward diffusion in the
Poisson-random-field setting (Sawyer & Hartl 1992) for the non-equilibrium
frequency spectrum (Evans, Shvets & Slatkin 2007, Theor. Popul. Biol.
71:109), using Scharfetter & Gummel (1969, IEEE Trans. Electron Devices
16:64) exponential fitting and Rannacher (1984, Numer. Math. 43:309)
startup. dadi (Gutenkunst et al. 2009) is the standard finite-difference
diffusion route and moments (Jouganous et al. 2017, Genetics 206:1549) the
moment-closure route for the same problem; this module takes the PDE route
so that the sample size enters only through the final projection.

Ancestral join: The ancestral calls and their upper/lower-case confidence convention
are those of the Ensembl EPO (Enredo-Pecan-Ortheus) pipeline: Paten et al.
(2008) Genome Res. 18:1814 (doi:10.1101/gr.076554.108) and 18:1829
(doi:10.1101/gr.076521.108); Herrero et al. (2016) Database bav096
(doi:10.1093/database/bav096); the same convention annotates the ancestral
allele in the 1000 Genomes Project releases (Nature 526:68, 2015,
doi:10.1038/nature15393). This module only joins and counts.

Linked two-site spectrum: second moments / covariances of the SFS: Fu (1995,
Theor. Popul. Biol. 48:172) for Kingman; Birkner, Blath & Eldon (2013,
Genetics 195:1037) and Hobolth, Siri-Jegousse & Bladt (2019, Theor. Popul.
Biol. 127:16) under Lambda-coalescents; the joint spectrum of two linked
sites under Kingman: Ferretti, Klassmann, Raineri, Ramos-Onsins, Wiehe &
Achaz (2018, Theor. Popul. Biol. 123:70); the linked 2-SFS as a
discriminator of multiple-merger versus Kingman coalescence: Rice, Novembre
& Desai (2018, bioRxiv 461517) and its journal version Fenton, Rice,
Novembre & Desai (2025, Genetics 229:iyaf023). Block-counting process of
the coalescent: Tavare (1984, Theor. Popul. Biol. 26:119); variable
population size as a time change: Griffiths & Tavare (1994, Phil. Trans. R.
Soc. B 344:403). Lambda-coalescents and the rate families as under
`popcorn.lambda_coalescent`. Monte Carlo reference: msprime (Baumdicker et
al. 2022, Genetics 220:iyab229). Bernstein subdivision certificates are
classical (see `popcorn.certificates`).

Folded two-site spectra: folded site classes: Fu (1995, Theor. Popul. Biol.
48:172). The two-site frequency spectrum of linked sites: Ferretti,
Klassmann, Raineri, Ramos-Onsins, Wiehe & Achaz (2018, Theor. Popul. Biol.
123:70); branch-length moments under Lambda-coalescents: Birkner, Blath &
Eldon (2013, Genetics 195:1037). Bernstein-basis positivity certificates
rest on the classical convex-hull and subdivision properties (Farouki 2012,
Comput. Aided Geom. Des. 29:379). The folding map, the pullback
certificates and the conic obstruction are this package's own
constructions.

Two-site genotype spectra: the pair projection is the
multivariate-hypergeometric (subsampling without replacement) analogue,
over individuals and two sites jointly, of the frequency-spectrum
projection of Marth, Czabarka, Murvai & Sherry (2004, Genetics 166:351) as
used in dadi (Gutenkunst et al. 2009, PLoS Genet. 5:e1000695); the object
it estimates is the two-locus sampling distribution of Hudson (2001,
Genetics 159:1805). Pairs a few bp apart are set aside because of
multinucleotide mutation events (Schrider, Hourmozdi & Hahn 2011, Curr.
Biol. 21:1051). Block jackknife: Efron & Tibshirani (1993), An Introduction
to the Bootstrap, Chapman & Hall. This module counts and projects; it fits
nothing.

Transient enclosures: the moment equations are the non-equilibrium
frequency-spectrum system of Evans, Shvets & Slatkin (2007, Theor. Popul.
Biol. 71:109) for the Kimura (1964) diffusion in the Poisson-random-field
setting (Sawyer & Hartl 1992); Živković, Steinrücken, Song & Stephan (2015,
Genetics 200:601) treat the same transient selected spectrum with variable
population size by spectral methods, and moments (Jouganous et al. 2017,
Genetics 206:1549) integrates the moment system in floats with a jackknife
closure. This module instead closes the truncation two-sidedly (the Metzler
bracket lemma stated in the engine docstring) and propagates it in the ball
arithmetic of Arb (Johansson 2017, IEEE Trans. Comput. 66:1281,
doi:10.1109/TC.2017.2690633) through python-flint.

Planted truths: the expected SFS under piecewise-constant size uses
Tavaré's (1984, Theor. Popul. Biol. 26:119) ancestral-process transition
function with the branch-class weights of Fu (1995); explicit expressions of
this kind are Polanski & Kimmel (2003, Genetics 165:427); that different
histories can share a frequency spectrum is Myers, Fefferman & Patterson
(2008, Theor. Popul. Biol. 73:342); the site-density target is Watterson's
(1975) estimator. The one map-derived constant (r/mu via a mean of 1.96
cM/Mb) is the length-weighted chromosome-21 mean of the UCSC GRCh38
recombAvg track, i.e. the sex-averaged deCODE genetic map of Halldorsson et
al. (2019, Science 363:eaau1043); only that mean is used.
Simulation is msprime (Kelleher, Etheridge & McVean 2016, PLoS Comput.
Biol. 12:e1004842; Baumdicker et al. 2022, Genetics 220:iyab229) — we thank
the msprime and tskit developers (Jerome Kelleher and colleagues) for the
simulators every planted-truth leg rests on — with its
Beta (Schweinsberg 2003) and Dirac (Eldon & Wakeley 2006) coalescent
models; this module only plants, calibrates and records.

Two-window harness: the simulator is msprime (Baumdicker et al. 2022,
Genetics 220:iyab229) and the branch-mode allele frequency spectrum of a
genealogy is tskit's (Ralph, Thornton & Kelleher 2020, Genetics 215:779);
the two-locus ancestral process is Griffiths (1981, Theor. Popul. Biol.
19:169) and Hudson (1983, Theor. Popul. Biol. 23:183); gene conversion in
the coalescent is Wiuf & Hein (2000, Genetics 155:451); the sweep model
msprime implements is the structured coalescent of Braverman et al. (1995,
Genetics 140:783) as in discoal (Kern & Schrider 2016, Bioinformatics
32:3839); the rho = 0 exact comparand is the Kingman E[L_i L_j] of Fu (1995,
Theor. Popul. Biol. 48:172) as computed by `popcorn.twosfs`. This module
only arranges the windows, the seed/batch law, the pooling and the
jackknife.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).
