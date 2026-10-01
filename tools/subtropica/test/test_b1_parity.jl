# test_b1_parity.jl — INTEGRATION tests.
# (1) Parity assertion between subtract.jl's
#     _c7_* trop/restrict entry points (now DELEGATING) and unit-T's
#     exported tropical.jl semantics — with hand-pinned SPEC_B1 §8 values,
#     so a future de-delegation regression is caught, not just the aliasing.
# (2) b1_driver census routing (b1_ray_census): 1-var cases only — the
#     1-var analytic C4 branch needs NO Oscar, so this file stays cheap and
#     Oscar-free both standalone and under the suite.
#
# Standalone: ulimit -v 32505856; julia tools/subtropica/test/test_b1_parity.jl
# Suite: auto-included by runtests.jl (own sandbox module).

using Test
using Nemo
using JSON

if !@isdefined(EpsExp)
    include(joinpath(@__DIR__, "..", "src", "types.jl"))
end
if !@isdefined(RayClass)
    include(joinpath(@__DIR__, "..", "src", "b1_types.jl"))
end
if !@isdefined(trop_on_ray)
    include(joinpath(@__DIR__, "..", "src", "tropical.jl"))
end
if !@isdefined(_c7_J_of)
    include(joinpath(@__DIR__, "..", "src", "subtract.jl"))
end
if !@isdefined(tropical_data)     # 1-var analytic branch only (no Oscar call)
    include(joinpath(@__DIR__, "..", "src", "polytope.jl"))
end
if !@isdefined(_canon_atomvec)    # runtime deps of the driver's 0-var path
    include(joinpath(@__DIR__, "..", "src", "hlogexpr.jl"))
end
if !@isdefined(EpsExpandRefusal)  # _compositions + the typed refusal
    include(joinpath(@__DIR__, "..", "src", "eps_expand.jl"))
end
if !@isdefined(B2Prefactor)       # B2 continuation types — the driver's B2
    # route (subtropica_integrate_b2/_b2_apply_pref) annotates
    # signatures with them, so b1_driver.jl needs continue.jl at include time
    include(joinpath(@__DIR__, "..", "src", "continue.jl"))
end
if !@isdefined(b1_ray_census)
    include(joinpath(@__DIR__, "..", "src", "b1_driver.jl"))
end

pt_mkE(nu, polys, vars, ring) =
    EulerIntegrand(Prefactor(), nu, polys, vars, Symbol[], ring)

@testset "B1 integration: C7↔C5 parity + census routing" begin

# --------------------------------------------------------------- parity (ii)
@testset "_c7_* == tropical.jl exports (F1/F3 pins)" begin
    R2, (x1, x2) = polynomial_ring(QQ, ["x1", "x2"])
    P1 = 1 + x1 + x2                                    # F1 (SPEC_B1 §8.1)
    P3 = 1 + x2 + x1 * x2                               # F3 (SPEC_B1 §8.3)
    E1 = pt_mkE([EpsExp(0, 1), EpsExp(0, 1)],           # dlog ν = (ε, ε)
                [(P1, EpsExp(-1, -3))], [:x1, :x2], R2)
    rays = [BigInt[-1, 0], BigInt[0, -1], BigInt[1, 1], BigInt[1, -1]]
    for rho in rays
        @test _c7_trop_poly(P1, rho, 2) == trop_poly(P1, rho, 2)
        @test _c7_trop_poly(P3, rho, 2) == trop_poly(P3, rho, 2)
        @test _c7_trop_on_ray(E1, rho, 2) == trop_on_ray(E1, rho)
        @test _c7_restrict_poly(P1, rho, 2) == restrict_poly(P1, rho, 2)
        @test _c7_restrict_poly(P3, rho, 2) == restrict_poly(P3, rho, 2)
    end
    # hand pins (SPEC_B1 §8.1/§8.3 tables) — teeth against de-delegation:
    @test trop_on_ray(E1, BigInt[-1, 0]) == EpsExp(0, -1)          # TropI(ρ1) = −ε
    @test trop_on_ray(E1, BigInt[1, 1]) == EpsExp(-1, -1)          # −1−ε (convergent)
    @test _c7_restrict_poly(P1, BigInt[-1, 0], 2) == 1 + x2        # F1 restriction
    @test _c7_restrict_poly(P3, BigInt[1, -1], 2) == 1 + x1 * x2   # F3 along ρ_C
end

# ------------------------------------------------------- census routing (iv)
@testset "b1_ray_census routing — 1-var (Oscar-free analytic branch)" begin
    R1, (x,) = polynomial_ring(QQ, ["x"])

    # (a) log-divergent at x→0: FLAT ν = ε−1  ⇒ dlog ν = ε ⇒ TropI(−1) = −ε
    Ediv = pt_mkE([EpsExp(-1, 1)], [(1 + x, EpsExp(-1, -3))], [:x], R1)
    c = b1_ray_census(Ediv)
    @test c.route === :b1
    @test c.dd.div_facets == [1] && length(c.dd.sigma_div) == 2
    @test c.Ed.nu[1] == EpsExp(0, 1)                    # measure bridge ν+1
    # C7 through the census objects: 3 counterterms, WARNING-1 signs
    cts = subtraction_terms(c.Ed, c.td, c.dd.sigma_div, c.dd.w)
    @test length(cts) == 3
    @test [ct.sign for ct in cts] == [1, -1, 1]         # ∅-id, ∅-{1}, face-{1}
    @test cts[3].face == [1] && isempty(cts[3].face_vars)   # 0-var face, J = []
    @test cts[3].vol_det == 1 && cts[3].vol_trops == [EpsExp(0, 1)]  # Vol = 1/ε
    # FLAT-measure mono on the ∅ face: ν_flat = ε−1 comes back out
    @test cts[1].mono == [EpsExp(-1, 1)]

    # (b) unregulated (TropI ≡ 0 ray) ⇒ Phase-A fall-through, never silent B1
    Exinv = pt_mkE([EpsExp(-1)], [(1 + x, EpsExp(-1))], [:x], R1)
    cx = b1_ray_census(Exinv)
    @test cx.route === :phase_a
    @test occursin("TropIllDefined", cx.note)

    # (c) power divergence ⇒ TYPED refusal naming B2 (PowerDivergentRefusal)
    Epow = pt_mkE([EpsExp(-3, 1)], [(1 + x, EpsExp(-1, -3))], [:x], R1)
    @test_throws PowerDivergentRefusal b1_ray_census(Epow)
    err = try b1_ray_census(Epow); nothing catch e; e end
    @test err isa PowerDivergentRefusal
    @test occursin("Nilsson-Passare", sprint(showerror, err))  # B2 named
    @test occursin("allow_continuation", sprint(showerror, err))  # B2 route named

    # (d) no non-constant polys ⇒ no polytope ⇒ Phase-A path
    Emono = pt_mkE([EpsExp(0, 1)], Tuple{QQMPolyRingElem,EpsExp}[], [:x], R1)
    @test b1_ray_census(Emono).route === :phase_a

    # (e) convergent input ⇒ Phase-A path with census certificate
    Econv = pt_mkE([EpsExp(0, 2)], [(1 + x, EpsExp(-2, -3))], [:x], R1)
    cc = b1_ray_census(Econv)
    @test cc.route === :phase_a && occursin("0 divergent", cc.note)
end

# ------------------------------------------- 0-var face constant-series unit
@testset "0-var face exact series (_b1_const_ct_orders)" begin
    R1, (x,) = polynomial_ring(QQ, ["x"])
    # ct with one constant poly 2^(−1−3ε): series 1/2·exp(−3ε log 2)
    ct = CounterTerm(Int[1], Int[], 1, big(1), [EpsExp(0, 1)], Int[],
                     Symbol[], EpsExp[], [EpsExp(0)],
                     [(R1(QQ(2)), EpsExp(-1, -3))])
    os = _b1_const_ct_orders(ct, Rational{BigInt}(1), 2)
    @test length(os) == 3
    @test os[1].terms[1].coef == "1/2" && isempty(os[1].terms[1].atoms)
    # ε¹: −3/2·Log2 (HLogA("2") canonicalizes to HConst(:Log2))
    @test os[2].terms[1].coef == "-3/2"
    @test os[2].terms[1].atoms == [HConst(:Log2) => 1]
    # ε²: 9/4·Log2²
    @test os[3].terms[1].coef == "9/4"
    @test os[3].terms[1].atoms == [HConst(:Log2) => 2]
    # kinvar-dependent 0-var poly refuses loudly
    ctbad = CounterTerm(Int[1], Int[], 1, big(1), [EpsExp(0, 1)], Int[],
                        Symbol[], EpsExp[], [EpsExp(0)],
                        [(1 + x, EpsExp(-1, -3))])
    @test_throws EpsExpandRefusal _b1_const_ct_orders(ctbad, Rational{BigInt}(1), 1)
end

end # top testset
