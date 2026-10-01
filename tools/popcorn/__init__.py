"""POPCORN — certified population-genetics likelihoods (package front door).

Certified likelihoods for the DFE: the Poisson Random Field selection-SFS/DFE
stack of the polyDFE/fitdadi/fastDFE class made exact and certified — exact
QQ(e^{-S}) SFS vectors + exact gradients, ball-certified vectors to n = 10^5,
a certified DFE-mixing kernel, dominance h != 1/2 (an entire q-deformation of
the Kummer line), plus an exact certificate-instrument family
(positivity / LP / hull / region certificates).

Modules:
  popcorn.sfs          exact/certified selection-SFS spectra (sfs_engine,
                       arb_route)
  popcorn.dfe          certified DFE-mixing kernel + certified quadrature
                       (dfe_layer, certquad)
  popcorn.dominance    dominance evaluators, h != 1/2 (dominance_oracle,
                       dominance_qseries)
  popcorn.certificates exact certificate family (positivity, exact_lp,
                       fast_lp, cone_lp_b host, cert_lp, region) — lazy
  popcorn.lambda_coalescent  exact Lambda-coalescent merger rates (Kingman,
                       Beta(2-a,a), Dirac) and expected SFS E[L_i] (lambda_exact)
  popcorn.twolocus     two-locus branch-length moments E[T_i^A T_j^B] under
                       recombination + piecewise N(t) (twolocus_engine; float)
  popcorn.transient    transient selected SFS through piecewise-constant N(t),
                       large n; float-validated PDE route (transient_sfs_engine)
  popcorn.ancestral    ancestral-state join against an EPO ancestral FASTA:
                       classify, orient, polarized table + counters (epo_join)
  popcorn.twosfs       exact linked two-site frequency spectrum E[L_i L_j]
                       under Lambda-coalescents + exact certificates against
                       the variable-size Kingman class (twosfs_engine,
                       twosfs_hull, twosfs_certificates; sympy for the
                       certificates, numpy+scipy propose only)
  popcorn.foldgate     minor-allele folding of a two-site spectrum + exact
                       deciders vs the folded variable-size Kingman class
                       (foldgate_engine; exact)
  popcorn.twosite      exact hypergeometric projection of two-site genotype
                       tables to fixed n (twosite_projection; EXACT) and the
                       two-site frequency spectrum of a diploid panel from
                       VCF/TSV genotypes, binned by cM or bp, block jackknife
                       (twosite_spectrum; data utility)
  popcorn.enclosure    certified enclosures (Arb balls) of the TRANSIENT selected
                       SFS through piecewise-constant N(t), S <= 0, small n;
                       model-conditional (transient_enclosure_engine) — lazy
  popcorn.planted      planted truths for SFS / two-site estimators: closed-form
                       (float) expected spectra under piecewise-constant N(t)
                       and seeded msprime generators (synth_truths; simulation)
  popcorn.twowindow    two-window / two-locus msprime harness for the two-site
                       spectrum: pooled pair vector + jackknife SEs, measured
                       decorrelation fraction (twowindow_sim; simulation)
CLI:
  popcorn/certsfs.py   hardened certified-SFS front door (entry/vector/check/
                       selftest; sha-pin verified, fail-closed rc=3)
Battery:
  popcorn/selftest.py  fast default legs (~40 s); heavier gates behind --full

Scope: certified likelihoods for the DFE. NOT-FOR: certified demography (the
certified stack is equilibrium constant-N; piecewise-constant N(t) certified
only for the transient selected SFS at small n via popcorn.enclosure
(model-conditional Arb enclosures, genic S <= 0), float-validated in
popcorn.transient (large n) / popcorn.twolocus; lambda_coalescent is
neutral constant-N; twosfs/foldgate certify statements about the whole
variable-size Kingman class for a model's linked 2-SFS at small n but fit no
history; twosite is model-free bookkeeping of observed genotypes (exact
projection arithmetic, not a certified likelihood); planted and twowindow
are msprime simulation harnesses (SIMULATION register: planted truths with
Monte Carlo SEs and a two-window Monte Carlo with jackknife SEs) that plant
and measure but never score, fit or certify); the 1F1 closed form for the
fixed-S SFS is the classical PRF result (Sawyer & Hartl 1992; Bustamante et al.
2001), not ours. Conventions: S = 4*Ne*s, S>0
advantageous, unfolded SFS; polyDFE/fastDFE S == this S; dadi/fitdadi
gamma == S/2.
"""
from ._pins import PKG_DIR, verify

__version__ = "1.0.0"

_SUBMODULES = ("sfs", "dfe", "dominance", "certificates",
               "lambda_coalescent", "twolocus", "transient", "ancestral",
               "twosfs", "foldgate", "twosite",
               "enclosure", "planted", "twowindow")


def __getattr__(name):
    if name in _SUBMODULES:
        import importlib
        mod = importlib.import_module(f".{name}", __name__)
        globals()[name] = mod
        return mod
    raise AttributeError(f"module 'popcorn' has no attribute {name!r}")


__all__ = list(_SUBMODULES) + ["PKG_DIR", "verify", "__version__"]
