"""popcorn.sfs — certified selection-SFS engine core (exact + ball routes).

Aliases, in place and by identity (see popcorn._loader):
  sfs_engine — exact QQ(e^{-S}) tridiagonal routes: M_i = alpha_i + beta_i
               e^{-S} exactly (two rational tridiagonal solves); derivatives =
               one more solve each, same matrix. n <= 2000 practical.
               mpmath only.
  arb_route  — arb ball recursion (needs python-flint; loaded lazily so the
               exact route works without it); ball-certified vectors +
               gradients to n = 10^5.

CLI front door: certsfs.py (entry/vector/check/selftest). Held-out two-route
gate: gate_phase1.py (battery --full).
"""
from ._loader import import_in_place

sfs_engine = import_in_place("sfs_engine")

# Convenience re-exports.
solve_M_exact = sfs_engine.solve_M_exact
solve_dM_exact = sfs_engine.solve_dM_exact
assemble_M = sfs_engine.assemble_M

__all__ = ["sfs_engine", "arb_route", "solve_M_exact", "solve_dM_exact",
           "assemble_M"]


def __getattr__(name):
    if name == "arb_route":                      # lazy: needs python-flint
        mod = import_in_place("arb_route")
        globals()["arb_route"] = mod
        return mod
    raise AttributeError(f"module 'popcorn.sfs' has no attribute {name!r}")
