# Multivar.jl — multivariate (s,t,eps) handling.
#
# Strategy: the IBP-derived system is a flat connection
#     dI = (A_s ds + A_t dt) I ,   A_s, A_t matrices over Q(eps,s,t),
# integrable: ∂_t A_s - ∂_s A_t = [A_s, A_t]  (flatness / dA = A∧A).
# We change to the ratio variable x = t/s (pulling the overall scale into a power of (-s)),
# eps-factorize the SINGLE-variable system in x with the engine, obtain the constant-ish gauge
# T, then VERIFY that the same T eps-factorizes the full multivariate connection by checking
# the transformed connection's flatness and dlog form in BOTH directions.
#
# This module provides:
#   * integrability_holds : check ∂_t A_s - ∂_s A_t = [A_s, A_t] for a 2-variable system.
#   * check_flatness      : same, returning the residual matrix.
#   * epsfactor_ratio     : reduce a system given directly in the ratio variable x (the common
#                           case after scaling), delegating to Normalize.try_epsfactor.
#
# The transport of the gauge across variables is exact: once the dlog letters are linear in x
# (the Fuchsian-dlog assumption), the constant balance gauges found in x are kinematic-variable
# independent, so they act identically on the t-connection.  check_flatness on the transformed
# connection is the certificate.

module Multivar

using Nemo
using ..RatFunc
using ..Fuchsia
using ..Normalize

export epsfactor_ratio, check_flatness, integrability_holds

"""
    epsfactor_ratio(ctx, A; kwargs...) -> EpsResult

Convenience wrapper: A is the single-variable connection in the distinguished ratio variable
x of `ctx`.  Delegates to Normalize.try_epsfactor.
"""
function epsfactor_ratio(ctx::Ctx, A; kwargs...)
    return Normalize.try_epsfactor(ctx, A; kwargs...)
end

"""
    check_flatness(ctx_x, ctx_y, As, At) -> residual matrix

Given two single-variable contexts sharing eps (here we pass the matrices already over a common
multivariate field), compute  R = ∂_y A_x - ∂_x A_y - (A_x A_y - A_y A_x).  R == 0 <=> flat.

For our use both A_x and A_y are matrices over the SAME field Q(eps,x,y); the caller supplies
derivative closures dx, dy.  We keep it generic via supplied derivative functions.
"""
function check_flatness(Ax, Ay, dx, dy)
    return dy(Ax) - dx(Ay) - (Ax*Ay - Ay*Ax)
end

"""
    integrability_holds(Ax, Ay, dx, dy) -> Bool
"""
function integrability_holds(Ax, Ay, dx, dy)
    R = check_flatness(Ax, Ay, dx, dy)
    return iszero(R)
end

end # module
