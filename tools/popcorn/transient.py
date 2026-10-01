"""popcorn.transient — transient selected SFS for large samples (float PDE route).

Aliases, in place and by identity (see popcorn._loader):
  transient_sfs_engine — expected unfolded SFS E_i (i = 1..n-1, theta = 1)
                         under genic selection S and a piecewise-constant
                         size history [(nu, T), ...] after an ancestral
                         Wright equilibrium: Scharfetter-Gummel fluxes on a
                         fixed log/lin/log grid, Crank-Nicolson with
                         Rannacher startup, banded LAPACK solves (one
                         factorization per epoch), binomial projection to
                         sample size n; exact hypergeometric down-sampling
                         n -> m and folding utilities. numpy + scipy.

REGISTER: FLOAT-VALIDATED, never certified. Measured trust radius (data
under reference/transient/): stationary spectra at n = 1000 agree with the
certified two-route references to <= 9.4e-5 max-relative in every frequency
band for S in {0, -1, -5, -20, -100}; transient self-convergence on a
two-epoch history (entries E_i > 1e-5): dt refinement <= 3.9e-4 for
|S| <= 20, 4.0e-3 at S = -100, 1.1e-2 at S = -200. Other regimes are
unmeasured. For stationary questions use the certified engine (popcorn.sfs,
certsfs.py); label every number from this module as float when it appears
beside certified ones.
"""
from ._loader import import_in_place

transient_sfs_engine = import_in_place("transient_sfs_engine")

# Convenience re-exports.
TransientSFSEngine = transient_sfs_engine.TransientSFSEngine
project_down = transient_sfs_engine.project_down
fold_sfs = transient_sfs_engine.fold_sfs
make_grid = transient_sfs_engine.make_grid
ln_g_ratio = transient_sfs_engine.ln_g_ratio

__all__ = ["transient_sfs_engine", "TransientSFSEngine", "project_down",
           "fold_sfs", "make_grid", "ln_g_ratio"]
