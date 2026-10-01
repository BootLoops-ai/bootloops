"""popcorn.foldgate — minor-allele folding of a two-site frequency spectrum, and
exact deciders for a folded spectrum against the folded Kingman class.

Aliases, in place and by identity (see popcorn._loader):
  foldgate_engine — the folding map Phi on stored symmetric coordinates
                    (fold_coeffs, apply_fold, pullback = Phi^T, exact_rank,
                    rank_report; standard library, exact) and the exact
                    deciders for a folded target {n, M, tot, q, phi}
                    (make_object) against the folded image of the
                    variable-population-size Kingman class: inside_hunt
                    (ANNIHILATED: exact convex combination of folded
                    single-atom components), witness_hunt
                    (SURVIVES_OUT_FULLCLASS: folded functional with a
                    Bernstein-certified nonnegative pullback kernel and an
                    exact negative margin), outer_membership
                    (NO_LINEAR_WITNESS: exact conic obstruction),
                    diagonal_witness_hunt (outside the hull of all
                    single-atom poolings), realization_hunt_n4, and the
                    chain decide(); exact building blocks solve_exact,
                    bern_1d_nonneg, bern_2d_nonneg; loaders for the shipped
                    Kingman-class polynomials and the sixteen-target
                    reference table (n = 4, 5, 6) under reference/foldgate/.

Register: EXACT (Fractions / sympy Rationals decide; scipy's HiGHS LPs on float
grids only propose, and a proposal that fails exact reconstruction or exact
certification is discarded). The folding map, the exact verifiers and the
loaders need the standard library only; the hunts import numpy, scipy and
sympy on first call (legs foldgate_deciders SKIP BY NAME without them).

CLI: foldgate_engine.py [--n {4,5,6}] [--out FILE] re-derives every reference
decision from the shipped ingredients (prints; writes only with --out).
"""
from ._loader import import_in_place

foldgate_engine = import_in_place("foldgate_engine")

# Convenience re-exports.
PAIRS = foldgate_engine.PAIRS
FPAIRS = foldgate_engine.FPAIRS
fold_coeffs = foldgate_engine.fold_coeffs
apply_fold = foldgate_engine.apply_fold
pullback = foldgate_engine.pullback
exact_rank = foldgate_engine.exact_rank
rank_report = foldgate_engine.rank_report
mults = foldgate_engine.mults
full_sum = foldgate_engine.full_sum
make_object = foldgate_engine.make_object
inside_hunt = foldgate_engine.inside_hunt
witness_hunt = foldgate_engine.witness_hunt
outer_membership = foldgate_engine.outer_membership
diagonal_witness_hunt = foldgate_engine.diagonal_witness_hunt
realization_hunt_n4 = foldgate_engine.realization_hunt_n4
decide = foldgate_engine.decide
solve_exact = foldgate_engine.solve_exact
bern_1d_nonneg = foldgate_engine.bern_1d_nonneg
bern_2d_nonneg = foldgate_engine.bern_2d_nonneg
poly_coeffs_2d = foldgate_engine.poly_coeffs_2d
load_reference_rows = foldgate_engine.load_reference_rows
load_ingredients = foldgate_engine.load_ingredients
component_q = foldgate_engine.component_q
folded_components = foldgate_engine.folded_components
m2_sympy = foldgate_engine.m2_sympy
kernel_sympy = foldgate_engine.kernel_sympy
fold_exprs = foldgate_engine.fold_exprs
tgrid_of = foldgate_engine.tgrid_of
lam_of = foldgate_engine.lam_of

__all__ = ["foldgate_engine", "PAIRS", "FPAIRS", "fold_coeffs", "apply_fold",
           "pullback", "exact_rank", "rank_report", "mults", "full_sum",
           "make_object", "inside_hunt", "witness_hunt", "outer_membership",
           "diagonal_witness_hunt", "realization_hunt_n4", "decide",
           "solve_exact", "bern_1d_nonneg", "bern_2d_nonneg", "poly_coeffs_2d",
           "load_reference_rows", "load_ingredients", "component_q",
           "folded_components", "m2_sympy", "kernel_sympy", "fold_exprs",
           "tgrid_of", "lam_of"]
