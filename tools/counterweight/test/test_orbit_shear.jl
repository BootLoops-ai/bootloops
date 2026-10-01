#!/usr/bin/env julia
# test_orbit_shear.jl — unit tests for the same-orbit ± shear (the single rational orbit
# move T_q = M(x)/q̂ at an irreducible quadratic letter, Normalize._orbit_balance!).
#
# Strategy: a 2×2 Fuchsian connection whose ε⁰ obstruction is a ± eigenvalue pair sitting
# AT the conjugate roots of q̂ = x²+1 — the configuration every projector-based orbit move
# (rank-1 joint-eigvec, rank-2 Tr-pair) is blind to, because those only transfer a-units
# between the orbit and OUTSIDE dump points.  The shear must fire, stay rational
# (T = M(x)/q̂, deg M ≤ 2, det T = 1), preserve Fuchsianity, and zero the orbit's ε⁰
# eigenvalues.  Controls: the nonzero-odd (Liouville-irremovable) class must be refused,
# and a wrong-γ mutation of the accepted gauge must fail the Fuchsian gate.
#
# Run:  julia --project=. test/test_orbit_shear.jl
#       (or both files in one process: julia --project=. test/runtests.jl — the
#        registered battery; the include below is guarded so that reuses one load)

isdefined(Main, :Counterweight) ||
    include(joinpath(@__DIR__, "..", "src", "Counterweight.jl"))
using .Counterweight
using .Counterweight.RatFunc, .Counterweight.Fuchsia, .Counterweight.Normalize
using Nemo

ctx = make_ctx()
x   = ctx.x
Qe  = ctx.Qe
Qex = ctx.Qex

println("="^70)
println("Same-orbit ± shear (T_q = M(x)/q̂) — unit tests")
println("="^70)

# Fixture: q̂ = x²+1 (irreducible), A = (2x/q̂)·R0 + (1/x)·E12 with R0 = [[0,1],[1,0]].
# Orbit residue R(y) = R0: ε⁰ eigenvalues ±1 at BOTH conjugate roots (a ± pair on the
# orbit); the pole at 0 is nilpotent (a = 0) and ∞ carries no rational eigenvalue, so no
# projector move can strictly drop the algebraic-aware norm — only the shear can.
q  = x^2 + Qex(1)
cs = [Nemo.QQ(1), Nemo.QQ(0), Nemo.QQ(1)]           # monic x² + 0·x + 1
A = zero_matrix(Qex, 2, 2)
A[1,2] = 2*x//q + Qex(1)//x
A[2,1] = 2*x//q
@assert is_fuchsian(ctx, A)                "fixture must be Fuchsian"
obs0 = Normalize._obstruction_norm_alg(ctx, A)
@assert obs0 > Nemo.QQ(4)                  "fixture must carry the ± orbit pair (Σ|a| = 4)"

# ----------------------------------------------------------------------
# (S1) K-arithmetic: the constant orbit Laurent term B(y) = B₀ + y·B₁.
#      Hand value: B(y) = −y·(R0/2 + E12), i.e. B₀ = 0, B₁ = [[0,−3/2],[−1/2,0]].
# ----------------------------------------------------------------------
println("\n[S1] constant orbit Laurent term in K components")
BK = Normalize._orbit_const_term(ctx, A, cs)
@assert BK !== nothing                     "S1: const-term extraction failed"
B0, B1 = BK
@assert iszero(B0)                         "S1: B0 must vanish for this fixture"
B1ref = matrix(Qe, 2, 2, [Qe(0) Qe(-3//2); Qe(-1//2) Qe(0)])
@assert B1 == B1ref                        "S1: B1 mismatch vs hand value"
println("    B(y) = y*", B1, "   => PASS")

# ----------------------------------------------------------------------
# (S2) The shear fires and is the promised move: T = M(x)/q̂, deg M ≤ 2, det T = 1,
#      Fuchsian + strict algebraic-norm drop + orbit ε⁰ residue nilpotent afterwards.
# ----------------------------------------------------------------------
println("\n[S2] same-orbit ± shear fires on the ± pair")
moved, Anew, Td = Normalize._orbit_balance!(ctx, A; verbose=false)
@assert moved                              "S2: no orbit move found"
T, desc = Td
println("    move: ", desc)
@assert occursin("same-orbit", desc)       "S2: a projector move fired instead of the shear"
@assert det(T) == Qex(1)                   "S2: det T must be exactly 1"
# T − I = N(x)/q̂ with N of x-degree ≤ 2: every entry's denominator divides q̂
qpoly = numerator(q)
for i in 1:2, j in 1:2
    e = T[i,j] - (i == j ? Qex(1) : Qex(0))
    iszero(e) && continue
    @assert iszero(mod(qpoly, denominator(e)))  "S2: T−I not of the M(x)/q̂ form"
    @assert degree(numerator(e)) <= 2      "S2: numerator degree > 2"
end
@assert is_fuchsian(ctx, Anew)             "S2: result not Fuchsian"
@assert Normalize._qex_mat_eps_regular(ctx, Anew)  "S2: result not ε-regular"
obs1 = Normalize._obstruction_norm_alg(ctx, Anew)
@assert obs1 < obs0                        "S2: no strict norm drop"
@assert obs1 < Nemo.QQ(1//100)             "S2: orbit ± pair not cleared (Σ|a| should be 0)"
# post-move orbit ε⁰ residue is nilpotent (all eigenvalues 0 at both roots)
Rj2 = Normalize._orbit_residue_components(ctx, Anew, cs)
@assert Rj2 !== nothing
H0 = eps0_part(ctx, Normalize._orbit_hat_matrix(ctx, Rj2, cs))
cpH = Normalize.charpoly_q(H0)
@assert cpH == gen(parent(cpH))^4          "S2: orbit ε⁰ residue not nilpotent"
@assert apply_gauge(ctx, A, T) == Anew     "S2: gauge T inconsistent"
println("    obs(alg) ", obs0, " -> ", obs1, "   det T = 1   => PASS")

# ----------------------------------------------------------------------
# (S3) Driver wiring: try_epsfactor reaches the shear through the public API.
#      (The ε-free toy stops afterwards by design — the leg checks the MOVE is wired,
#      the accepted history names it, and the recorded total gauge is exact.)
# ----------------------------------------------------------------------
println("\n[S3] public driver reaches the shear")
res = try_epsfactor(ctx, A; max_balance_rounds=6, verbose=false)
hist = get(res.report, "history", String[])
@assert !isempty(hist)                     "S3: driver accepted no move"
@assert occursin("same-orbit", hist[1])    "S3: first accepted move is not the shear"
@assert apply_gauge(ctx, A, res.T) == res.Anew  "S3: total gauge inconsistent"
@assert Normalize._obstruction_norm_alg(ctx, res.Anew) < Nemo.QQ(1//100)  "S3: orbit obstruction survived"
println("    history[1] = ", hist[1])
println("    => PASS")

# ----------------------------------------------------------------------
# (S4) Control — nonzero ODD orbit residue (Liouville-irremovable class): no move.
#      R(y) = y·[[0,1],[−1,0]] has eigvals ±1 at the roots but purely odd structure;
#      the rational shear does not exist and nothing may pretend it does.
# ----------------------------------------------------------------------
println("\n[S4] odd-class control refuses")
B = zero_matrix(Qex, 2, 2)
B[1,2] = Qex(-2)//q
B[2,1] = Qex(2)//q
@assert is_fuchsian(ctx, B)
movedb, _, _ = Normalize._orbit_balance!(ctx, B; verbose=false)
@assert !movedb                            "S4: a move fired on the irremovable odd class"
println("    moved = false   => PASS")

# ----------------------------------------------------------------------
# (S5) Mutation control — the γ condition is load-bearing: doubling the shear part
#      (wrong γ) must fail the Fuchsian gate (double pole at the roots survives).
# ----------------------------------------------------------------------
println("\n[S5] wrong-γ mutation fails the gate")
Twrong = identity_matrix(Qex, 2) + 2*(T - identity_matrix(Qex, 2))
Awrong = apply_gauge(ctx, A, Twrong)
@assert !is_fuchsian(ctx, Awrong)          "S5: wrong-γ gauge stayed Fuchsian — gate vacuous"
println("    wrong-γ result non-Fuchsian   => PASS")

println("\n" * "="^70)
println("ALL same-orbit shear unit tests PASS (5/5)")
println("="^70)
