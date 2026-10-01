# Counterweight.jl — top-level module for the no-Mathematica epsilon-factorization engine.
#
# Given a coupled-DE sector's connection matrix A(eps,s,t) (built by Kira IBP), find a rotation
# T to an eps-factorized (canonical / UT) basis:
#     T A T^{-1} + (dT) T^{-1} = eps * Atilde,    Atilde eps-independent, dlog.
#
# This is BASIS CONSTRUCTION (preprocessing for the Landau bootstrap), NOT solving the DE across
# kinematics.  Engine = Lee's balance-transformation algorithm (arXiv:1411.0911) implemented over
# Q(eps)(x) using Nemo (Flint).  No Wolfram, no SageMath.
#
# Submodules:
#   RatFunc    — field tower Q(eps)(x), residue/pole analysis, gauge action helpers.
#   Fuchsia    — Fuchsian check, eps-form check, certification.
#   Normalize  — Lee balance driver: eigenvalue normalization + eps-decoupling.
#   Multivar   — multivariate (s,t,eps): eps-factorize in one ratio variable, verify flatness.
#   KiraDE     — parse Kira derive_dgl output into a connection matrix A(eps,s,t).

module Counterweight

using Nemo

include("RatFunc.jl")
include("Fuchsia.jl")
include("Normalize.jl")
include("Multivar.jl")
include("KiraDE.jl")
include("MpolyFeed.jl")
include(joinpath(@__DIR__, "..", "examples", "henn_dbox.jl"))

using .RatFunc
using .Fuchsia
using .Normalize
using .Multivar
using .KiraDE
using .MpolyFeed
using .HennDbox

# re-export the public surface
export Ctx, make_ctx, F2x, matF2x, eps_of, x_of,
       singular_points_x, residue_matrix, residue_at_infinity, is_fuchsian, eval_x_at,
       pole_profile, laurent_coeff_matrix,
       is_epsform, certify_epsform, apply_gauge, gauge_balance,
       moser_reduce, MoserResult, dlog_alphabet,
       try_epsfactor, EpsResult,
       epsfactor_ratio, check_flatness, integrability_holds,
       load_kira_de, build_A_from_rows,
       feed_matrix_from_chunks, qex_from_mpoly_pair

end # module
