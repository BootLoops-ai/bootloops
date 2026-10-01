"""popcorn.twolocus — two-locus branch-length second moments E[T_i^A T_j^B]
under the coalescent with recombination and piecewise-constant N(t).

Aliases, in place and by identity (see popcorn._loader):
  twolocus_engine — TwoLocus(n): labeled two-locus configuration chain
                    (sparse generators built once per n); .moments(rho,
                    epochs) -> (M, m1A, m1B) with M[i-1][j-1] =
                    E[T_i^A T_j^B], i, j = 1..n-1, by a direct sparse LU
                    solve in the terminal epoch and expm_multiply through
                    finite epochs. Time in 2*N_ref generations, eta =
                    N_ref/N, rho = 4*N_ref*R. numpy + scipy only.
                    pool_pairs(M): normalized symmetrized i <= j vector;
                    epochs_from_Nt(times_gen, sizes, N_ref): history in
                    generations -> [(tau_start, eta)].

Register: FLOAT-VALIDATED (double precision; validated at rho = 0 against
exact-rational Kingman E[L_i L_j], at n = 2 against the closed-form
two-locus covariance, and by self-consistency across rho — not certified,
no enclosure). Reference data: reference/twolocus/. Cost grows ~2.9x in
states per added sample; n <= 7 is sub-second per constant-N evaluation,
n = 8 at rho > 0 is ~13 s and ~0.8 GB.
"""
import functools

from ._loader import import_in_place

twolocus_engine = import_in_place("twolocus_engine")

# Convenience re-exports.
TwoLocus = twolocus_engine.TwoLocus
pool_pairs = twolocus_engine.pool_pairs
epochs_from_Nt = twolocus_engine.epochs_from_Nt


@functools.lru_cache(maxsize=8)
def chain(n):
    """TwoLocus(n), built once per n and cached (the build is the only
    n-dependent setup; moments() can then be called at any rho / history)."""
    return TwoLocus(int(n))


def moments(n, rho, epochs=((0.0, 1.0),)):
    """(M, m1A, m1B) for sample size n at scaled recombination rate rho under
    the piecewise-constant history `epochs` = [(tau_start, eta), ...]
    (default: constant N). M[i-1][j-1] = E[T_i^A T_j^B]."""
    return chain(n).moments(float(rho), list(epochs))


__all__ = ["twolocus_engine", "TwoLocus", "pool_pairs", "epochs_from_Nt",
           "chain", "moments"]
