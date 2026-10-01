"""popcorn.twosfs — exact two-site frequency spectrum of linked sites under
Lambda-coalescents, and exact certificates against the variable-size Kingman
class.

Aliases, in place and by identity (see popcorn._loader):
  twosfs_engine — E[L_i], E[L_i L_j] (i, j = 1..n-1) for any exchangeable
                  Lambda-coalescent (Kingman, Beta(2-alpha, alpha), Dirac,
                  Kingman+Dirac mixture; rates from popcorn.lambda_coalescent)
                  by an exact first-step recursion on block-size multisets,
                  with two independent exact cross-routes; the normalized
                  linked 2-SFS q_ij = E[L_i L_j]/E[L_tot^2]; the Kingman
                  block-counting semigroup, single-time block-count moments
                  m1(t), m2(t) and the symmetrized two-time kernel Csym(s, u)
                  (sympy) that represent every variable-size Kingman
                  history's 1- and 2-SFS as integrals against one measure;
                  identity_checks(n). Standard library (+ sympy for the
                  kernel). Loaded eagerly.
  twosfs_hull   — exact Euclidean projection of q onto the hull of finitely
                  many normalized atom moments (project_hull), exact square
                  solve. Standard library. Loaded eagerly.
  twosfs_certificates — exact Bernstein certificates (bern_1d_nonneg,
                  bern_2d_nonneg, poly_coeffs_2d, m2_polys), the verdict
                  engine hunt / class_verdict (OUT_FULLCLASS: outside the
                  2-SFS set of every pooled variable-size Kingman history;
                  OUT_ATOMPOOLINGS; INSIDE_2SFS; UNDECIDED), certify_witness,
                  kingman_kernel, kernel_constant_size_check, load_reference.
                  Needs sympy at import and numpy + scipy for class_verdict;
                  loaded LAZILY on first access so that `import
                  popcorn.twosfs` needs the standard library only.

Register: EXACT (Fractions / sympy rationals for every moment, certificate,
margin and verdict). The only floats: the LP inside class_verdict/hunt,
which PROPOSES a witness or a support that exact arithmetic then accepts or
rejects, and display fields named *_float. Reference data:
reference/twosfs/ (exact hull projections and 1-SFS mixture certificates,
exact class certificates with rational witnesses, the identity-check
record, and a Monte Carlo cross-check against an independent simulator).

CLI: twosfs_engine.py [--n N --family F --param P --out FILE] (identity
checks + exact moments); twosfs_certificates.py [same flags] (class verdict).
"""
from ._loader import import_in_place

lambda_exact = import_in_place("lambda_exact")      # rate families (shared)
twosfs_engine = import_in_place("twosfs_engine")
twosfs_hull = import_in_place("twosfs_hull")

# Convenience re-exports (engine).
lambda_moments = twosfs_engine.lambda_moments
setpartition_moments = twosfs_engine.setpartition_moments
kingman_moments_route2 = twosfs_engine.kingman_moments_route2
component_moments = twosfs_engine.component_moments
death_prob = twosfs_engine.death_prob
death_semigroup_coeffs = twosfs_engine.death_semigroup_coeffs
jump_levels = twosfs_engine.jump_levels
cond_level_dist = twosfs_engine.cond_level_dist
twotime_kernel_poly = twosfs_engine.twotime_kernel_poly
normalize_matrix = twosfs_engine.normalize_matrix
flatten = twosfs_engine.flatten
identity_checks = twosfs_engine.identity_checks
rate_function = twosfs_engine.rate_function
beta_rate = lambda_exact.beta_rate
dirac_rate = lambda_exact.dirac_rate
kingman_rate = lambda_exact.kingman_rate
msprime_dirac_rate = lambda_exact.msprime_dirac_rate
# Convenience re-exports (hull).
project_hull = twosfs_hull.project_hull
solve_exact = twosfs_hull.solve_exact
pairs = twosfs_hull.pairs
dot = twosfs_hull.dot
q_of_moments = twosfs_hull.q_of_moments

_LAZY = ("poly_coeffs_2d", "bern_1d_nonneg", "bern_2d_nonneg", "m2_polys",
         "hunt", "kingman_kernel", "kernel_constant_size_check",
         "class_verdict", "certify_witness", "load_reference")


def __getattr__(name):
    if name == "twosfs_certificates" or name in _LAZY:
        mod = import_in_place("twosfs_certificates")      # needs sympy
        globals()["twosfs_certificates"] = mod
        for k in _LAZY:
            globals()[k] = getattr(mod, k)
        return mod if name == "twosfs_certificates" else globals()[name]
    raise AttributeError(f"module 'popcorn.twosfs' has no attribute {name!r}")


def q_matrix(n, lam):
    """(q, E[L_tot^2]): the exact normalized linked 2-SFS of the coalescent
    with merger rates lam at sample size n, q[(i, j)] for i <= j."""
    _u, v = lambda_moments(n, lam)
    return normalize_matrix(v, n)


__all__ = ["twosfs_engine", "twosfs_hull", "twosfs_certificates",
           "lambda_exact", "lambda_moments", "setpartition_moments",
           "kingman_moments_route2", "component_moments", "death_prob",
           "death_semigroup_coeffs", "jump_levels", "cond_level_dist",
           "twotime_kernel_poly", "normalize_matrix", "flatten",
           "identity_checks", "rate_function", "beta_rate", "dirac_rate",
           "kingman_rate", "msprime_dirac_rate", "project_hull",
           "solve_exact", "pairs", "dot", "q_of_moments", "q_matrix"] + list(_LAZY)
