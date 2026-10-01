#!/usr/bin/env julia
# test_moser.jl — unit tests for Barkatou–Moser pole-order reduction (Fuchsia.moser_reduce).
#
# Strategy: build a Fuchsian connection, scramble it with a rational diagonal gauge that
# introduces a higher-order pole (so the scrambled system is *guaranteed* regular-singular),
# and verify moser_reduce recovers a Fuchsian form with a consistent total gauge.
#
# Run:  julia --project=. test/test_moser.jl

include(joinpath(@__DIR__, "..", "src", "Counterweight.jl"))
using .Counterweight
using .Counterweight.RatFunc, .Counterweight.Fuchsia
using Nemo

ctx = make_ctx()
x   = ctx.x
eps = ctx.eps
Qe  = ctx.Qe
Qex = ctx.Qex

mat(n, ent) = matF2x(ctx, n, ent)
diagx(vals) = begin
    n = length(vals); D = zero_matrix(Qex, n, n)
    for i in 1:n; D[i,i] = Qex(vals[i]); end; D
end

println("="^70)
println("Moser pole-order reduction — unit tests")
println("="^70)

# ----------------------------------------------------------------------
# (T1) 2×2: Fuchsian seed scrambled by diag(x,1) → order-2 pole at 0, nilpotent A₀.
# ----------------------------------------------------------------------
println("\n[T1] 2×2, order-2 nilpotent pole at 0 (regular singular)")
R0 = matrix(Qe, 2, 2, [Qe(eps) Qe(1); Qe(1) Qe(-eps)])
R1 = matrix(Qe, 2, 2, [Qe(0) Qe(eps); Qe(eps) Qe(0)])
Afuchs = mat(2, (i,j)-> Qex(R0[i,j])//x + Qex(R1[i,j])//(x-Qex(1)))
@assert is_fuchsian(ctx, Afuchs)
S = diagx([x, Qex(1)])                       # scrambling gauge (zero at 0, pole at ∞)
A1 = apply_gauge(ctx, Afuchs, S)
prof1 = pole_profile(ctx, A1)
println("    scrambled profile: ", Fuchsia._profstr(prof1))
@assert prof1[Qe(0)] == 2          "T1 setup: expected order-2 pole at 0"
@assert !is_fuchsian(ctx, A1)

mr1 = moser_reduce(ctx, A1; max_iter=30, verbose=false)
println("    out: ", Fuchsia._profstr(mr1.profile_out), "   ok=", mr1.ok,
        "   iters=", mr1.iters)
@assert mr1.ok                          "T1: moser_reduce failed"
@assert is_fuchsian(ctx, mr1.A)         "T1: result not Fuchsian"
@assert apply_gauge(ctx, A1, mr1.T) == mr1.A   "T1: gauge T inconsistent"
println("    => PASS")

# ----------------------------------------------------------------------
# (T2) 3×3: scramble by diag(x², x, 1) → order-3 pole at 0; reduces in ≥2 steps.
# ----------------------------------------------------------------------
println("\n[T2] 3×3, order-3 nilpotent pole at 0 (regular singular)")
R = matrix(Qe, 3, 3, [Qe(eps) Qe(1) Qe(2);
                      Qe(1)   Qe(2*eps) Qe(1);
                      Qe(1)   Qe(1)   Qe(-3*eps)])
Afuchs3 = mat(3, (i,j)-> Qex(R[i,j])//x + Qex(i==j ? eps : 0)//(x-Qex(1)))
@assert is_fuchsian(ctx, Afuchs3)
S3 = diagx([x^2, x, Qex(1)])
A2 = apply_gauge(ctx, Afuchs3, S3)
prof2 = pole_profile(ctx, A2)
println("    scrambled profile: ", Fuchsia._profstr(prof2))
@assert prof2[Qe(0)] == 3          "T2 setup: expected order-3 pole at 0"

mr2 = moser_reduce(ctx, A2; max_iter=60, verbose=false)
println("    out: ", Fuchsia._profstr(mr2.profile_out), "   ok=", mr2.ok,
        "   iters=", mr2.iters)
@assert mr2.ok                          "T2: moser_reduce failed"
@assert is_fuchsian(ctx, mr2.A)         "T2: result not Fuchsian"
@assert mr2.iters >= 2                  "T2: expected ≥2 reduction steps"
@assert apply_gauge(ctx, A2, mr2.T) == mr2.A   "T2: gauge T inconsistent"
println("    => PASS")

# ----------------------------------------------------------------------
# (T3) Non-nilpotent leading matrix → returns diagnostic, does not loop.
# ----------------------------------------------------------------------
println("\n[T3] non-nilpotent leading (eigenvalue 1) → irregular diagnostic")
B = matrix(Qe, 2, 2, [Qe(1) Qe(0); Qe(0) Qe(0)])
A3 = mat(2, (i,j)-> Qex(B[i,j])//x^2 + Qex(i==j ? eps : 0)//x)
mr3 = moser_reduce(ctx, A3; max_iter=20, verbose=false)
println("    ok=", mr3.ok, "   stop=\"", mr3.stop, "\"")
println("    irregular=", mr3.irregular)
@assert !mr3.ok                         "T3: should have failed (irregular)"
@assert mr3.irregular !== nothing       "T3: missing irregular diagnostic"
@assert mr3.iters == 0                  "T3: should return immediately, no iteration"
@assert haskey(mr3.irregular, "leading_eigfactors")
println("    => PASS")

# ----------------------------------------------------------------------
# (T4) Pole at infinity (polynomial growth) → reduced via x→1/y substitution.
# ----------------------------------------------------------------------
println("\n[T4] order-3 pole at infinity")
S4 = diagx([Qex(1)//x, Qex(1)])              # zero at ∞, pole at 0
A4 = apply_gauge(ctx, Afuchs, S4)
prof4 = pole_profile(ctx, A4)
println("    scrambled profile: ", Fuchsia._profstr(prof4))
@assert prof4[:inf] >= 2           "T4 setup: expected higher-order pole at ∞"

mr4 = moser_reduce(ctx, A4; max_iter=30, verbose=false)
println("    out: ", Fuchsia._profstr(mr4.profile_out), "   ok=", mr4.ok,
        "   iters=", mr4.iters)
@assert mr4.ok                          "T4: moser_reduce at infinity failed"
@assert is_fuchsian(ctx, mr4.A)         "T4: result not Fuchsian"
@assert apply_gauge(ctx, A4, mr4.T) == mr4.A   "T4: gauge T inconsistent"
println("    => PASS")

println("\n" * "="^70)
println("ALL Moser unit tests PASS")
println("="^70)
