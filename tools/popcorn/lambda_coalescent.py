"""popcorn.lambda_coalescent — exact Lambda-coalescent rates and expected SFS.

Aliases, in place and by identity (see popcorn._loader):
  lambda_exact — merger rates lambda_{b,k} for Kingman, Beta(2-alpha, alpha),
                 Dirac(psi) and the Kingman+Dirac mixture as exact rationals;
                 the exact expected branch-length spectrum E[L_i] = h(n,i),
                 i = 1..n-1 (hence the expected unfolded SFS) by the standard
                 first-transition recursion; the normalized spectrum xi_hat;
                 built-in exact identity checks; the n = 20 reference grid and
                 Kingman-class functional under reference/lambda_coalescent/.
                 Register EXACT: fractions.Fraction end to end, standard
                 library only. Constant population size.

CLI: lambda_exact.py [--n N] [--out FILE] (identity checks + reference grid).
"""
from ._loader import import_in_place

lambda_exact = import_in_place("lambda_exact")

# Convenience re-exports.
N = lambda_exact.N
beta_rate = lambda_exact.beta_rate
dirac_rate = lambda_exact.dirac_rate
kingman_rate = lambda_exact.kingman_rate
msprime_dirac_rate = lambda_exact.msprime_dirac_rate
expected_lengths = lambda_exact.expected_lengths
xi_hat = lambda_exact.xi_hat
identity_checks = lambda_exact.identity_checks
reference_grid = lambda_exact.reference_grid
load_reference_spectra = lambda_exact.load_reference_spectra
load_witness = lambda_exact.load_witness
witness_margin = lambda_exact.witness_margin

__all__ = ["lambda_exact", "N", "beta_rate", "dirac_rate", "kingman_rate",
           "msprime_dirac_rate", "expected_lengths", "xi_hat",
           "identity_checks", "reference_grid", "load_reference_spectra",
           "load_witness", "witness_margin"]
