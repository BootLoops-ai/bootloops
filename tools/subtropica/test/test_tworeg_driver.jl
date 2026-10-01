# test_tworeg_driver.jl — two-regulator driver tests. The
# two-regulator driver route: subtropica_integrate(::Euler2Integrand) with the
# exact [q^aux_order] convolution (b1_driver.jl §B2Q; the eq 4.11 gate's
# MG7 decision transliterated — VERDICT_EQ411.md §3).
#
# Layers:
#  (a) HF-free units: pref q-series algebra, log-insertion multinomials,
#      q0-projection scope refusal, census/refusal law (ContinuationRequired
#      message names the flag), double-embed reconstruction on the mini.
#  (b) HF end-to-end (skipped loudly without the engine): the MINI TWO-REG
#      SYNTHETIC with hand-worked closed form,
#         J(q,eps) = ∫∫ x1^{-1-q} x2^{eps+q} (1+x1+x2)^{-1-2eps} dx1 dx2
#                  = Γ(-q)·Γ(1+eps+q)·Γ(eps)/Γ(1+2eps),
#      through subtropica_integrate ALONE:
#         [q^-1] J = -Γ(1+eps)Γ(eps)/Γ(1+2eps) = -1/eps + 0 + ζ₂·eps - ...
#         [q^0]  J = -(γ+ψ(1+eps))·Γ(1+eps)Γ(eps)/Γ(1+2eps)
#                  = -ζ₂ + ζ₃·eps + (π⁴/60)·eps² + ...
#      ([q^0] exercises the FULL machinery: unregulated-by-eps ray →
#      forced continuation → q-pole prefactor series → v=1 log insertions
#      through the B1 Möbius subtraction with reconstructed q-data.)
#
# Self-contained: run via runtests.jl (auto-included, own sandbox module) or
#   ulimit -v 32505856; julia --project=tools/subtropica tools/subtropica/test/test_tworeg_driver.jl

using Test
using SubTropica
using Nemo

@testset "two-reg driver" begin

    cfg = SubTropica.HFConfig()
    have_hf = isfile(cfg.bin) && isfile(cfg.mzv_data_path)
    have_hf || @warn "HyperFLINT engine not found at $(cfg.bin) — two-reg end-to-end testsets skipped"

    R2, (x1, x2) = polynomial_ring(Nemo.QQ, ["x1", "x2"])
    P = 1 + x1 + x2
    # THE MINI (flat measure): x1^{-1-q} x2^{eps+q} P^{-1-2eps}
    mkmini() = Euler2Integrand(Prefactor(),
        [Eps2(-1//1, 0//1, -1//1), Eps2(0//1, 1//1, 1//1)],
        [(P, Eps2(-1//1, -2//1, 0//1))],
        [:x1, :x2], Symbol[], R2, [:eps, :q])

    # ---------------------------------------------------------------- (a) units
    @testset "pref q-series (exact)" begin
        # continuation-shaped pref: (-1-2eps)/(q)  [num e_P, den TropI = q]
        p = B2Prefactor{Eps2}(1//1, [Eps2(-1//1, -2//1, 0//1)],
                              [Eps2(0//1, 0//1, 1//1)])
        (qmin, cc) = SubTropica._b2q_pref_qseries(p, 0)
        @test qmin == -1 && length(cc) == 2
        @test length(cc[1]) == 1 && isempty(cc[2])     # pure 1/q pole
        t = cc[1][1]
        @test t.c == 1 && t.num == [EpsExp(-1//1, -2//1)] && isempty(t.den)
        # scaled-q negative control: den 2q halves every coefficient
        p2 = B2Prefactor{Eps2}(1//1, [Eps2(-1//1, -2//1, 0//1)],
                               [Eps2(0//1, 0//1, 2//1)])
        (qmin2, cc2) = SubTropica._b2q_pref_qseries(p2, 0)
        @test qmin2 == -1 && cc2[1][1].c == 1//2
        # mixed denominator 1/(eps - q): geometric series in EpsExp prefs
        pm = B2Prefactor{Eps2}(1//1, Eps2[], [Eps2(0//1, 1//1, -1//1)])
        (m3, c3) = SubTropica._b2q_pref_qseries(pm, 2)
        @test m3 == 0 && length(c3) == 3
        @test SubTropica.b2pref_eval(c3[1][1], big(1)//3) == 3      # 1/eps
        @test SubTropica.b2pref_eval(c3[2][1], big(1)//3) == 9      # q/eps^2
        @test SubTropica.b2pref_eval(c3[3][1], big(1)//3) == 27     # q^2/eps^3
        # exactness cross-check at a rational (eps,q): the truncated series
        # differs from b2pref_eval by EXACTLY the geometric tail
        # exact·(q/eps)^3  (1/(eps-q) = Σ q^j/eps^{j+1})
        ee = big(1)//3; qq = big(1)//17
        exact = SubTropica.b2pref_eval(pm, ee, qq)
        approx = sum(SubTropica.b2pref_eval(c3[j][1], ee) * qq^(m3 + j - 1)
                     for j in 1:3)
        @test exact - approx == exact * (qq / ee)^3
    end

    @testset "log-insertion multinomials" begin
        items = Tuple{Rational{BigInt},QQMPolyRingElem}[
            (Rational{BigInt}(-1), x1), (Rational{BigInt}(2), x2)]
        T0 = SubTropica._b2q_log_power_terms(items, 0)
        @test T0 == [(1//1, Pair{QQMPolyRingElem,Int}[])]
        T1 = SubTropica._b2q_log_power_terms(items, 1)
        @test length(T1) == 2
        T2 = SubTropica._b2q_log_power_terms(items, 2)
        @test length(T2) == 3                                     # ½L² expansion
        @test any(t -> t[1] == 1//2 && t[2] == [x1 => 2], T2)
        @test any(t -> t[1] == -2//1 && t[2] == [x1 => 1, x2 => 1], T2)
        @test any(t -> t[1] == 2//1 && t[2] == [x2 => 2], T2)
        # constants: Log[1] kills; Log[2] survives; Log[-3] is LOUD
        @test isempty(SubTropica._b2q_log_power_terms(
            [(Rational{BigInt}(5), R2(1))], 1))
        Tc = SubTropica._b2q_log_power_terms([(Rational{BigInt}(5), R2(2))], 1)
        @test length(Tc) == 1 && Tc[1][1] == 5
        @test_throws SubTropica.EpsExpandRefusal SubTropica._b2q_log_power_terms(
            [(Rational{BigInt}(1), R2(-3))], 1)
    end

    @testset "q0-projection + scope refusal" begin
        E2 = mkmini()
        Ed = SubTropica._b2q_dlog(E2)
        (E0, cvec) = SubTropica._b2q_q0_flat(Ed)
        @test E0.nu == [EpsExp(-1//1, 0//1), EpsExp(0//1, 1//1)]  # flat, q=0
        @test cvec == [-1//1, 1//1]
        @test E0.polys[1][2] == EpsExp(-1//1, -2//1)
        # q-carrying POLY exponent refuses typed (v1 scope law)
        Ebad = Euler2Integrand(Prefactor(), [Eps2(1), Eps2(1)],
            [(P, Eps2(-1//1, -2//1, 1//1))], [:x1, :x2], Symbol[], R2,
            [:eps, :q])
        err = try SubTropica._b2q_q0_flat(Ebad); nothing catch e; e end
        @test err isa B2Refusal && err.kind == :AuxRegulatorInPolyExponent
    end

    @testset "census + ContinuationRequired message law" begin
        # default surface refuses TYPED and names the flag (no HF needed —
        # the refusal fires before any engine call... except hf_gate; skip
        # without the engine)
        if have_hf
            err = try subtropica_integrate(mkmini(); order=1, aux_order=0)
                nothing catch e; e end
            @test err isa B2Refusal && err.kind == :ContinuationRequired
            msg = sprint(showerror, err)
            @test occursin("allow_continuation=true", msg)
            @test occursin("unregulated_by_eps", msg)
        end
        # two-reg census values (trop2_on_ray — the mandated census surface)
        Ed = SubTropica._b2q_dlog(mkmini())
        @test trop2_on_ray(Ed, [1, 1])  == Eps2(0//1, -1//1, 0//1)   # -eps  (log)
        @test trop2_on_ray(Ed, [-1, 0]) == Eps2(0//1, 0//1, 1//1)    #  q    (UNREG)
        @test trop2_on_ray(Ed, [0, -1]) == Eps2(-1//1, -1//1, -1//1) # conv
        @test SubTropica.unregulated_by_eps(trop2_on_ray(Ed, [-1, 0]))
    end

    @testset "double-embed reconstruction (mini, exact)" begin
        # continued-by-hand mini term: nu_dlog = (1-q, 1+eps+q),
        # polys P^{-2-2eps}  (the (-1,0)-step output; poly exps q-free)
        Ed = Euler2Integrand(Prefactor(),
            [Eps2(1//1, 0//1, -1//1), Eps2(1//1, 1//1, 1//1)],
            [(P, Eps2(-2//1, -2//1, 0//1))],
            [:x1, :x2], Symbol[], R2, [:eps, :q])
        (E0, cvec) = SubTropica._b2q_q0_flat(Ed)
        census = SubTropica.b1_ray_census(E0)
        @test census.route === :b1
        rec = SubTropica._b2q_reconstruct(E0, cvec, census)
        @test length(rec.cts) >= 2                    # identity + face cts
        # faces: [] then [(1,1)-facet]; the singleton face keeps x2 (J=[2])
        faces = unique([ct.face for ct in rec.cts])
        @test Int[] in faces && length(faces) == 2
        for (i, ct) in enumerate(rec.cts)
            if ct.face == Int[]
                @test isempty(ct.subface) ? rec.mono_c[i] == [-1//1, 1//1] : true
            else
                # vol_trop = -TropI((1,1)) = eps exactly (q-free): C = 0
                @test all(iszero, rec.voltrop_c[i])
            end
            # jacobian/restricted exponents inherit NO q here (TropI q-free)
            @test all(iszero, rec.polyexp_c[i])
        end
    end

    # ------------------------------------------------------- (b) END-TO-END
    if have_hf
        # closed-form comparators (BigFloat, dps ~ 40)
        setprecision(BigFloat, 200) do
            z2 = BigFloat(pi)^2 / 6
            z3 = BigFloat("1.2020569031595942853997381615114499907649862923405")
            pi4_60 = BigFloat(pi)^4 / 60

            @testset "mini end-to-end [q^-1] (v=0 route)" begin
                L = subtropica_integrate(mkmini(); order=2, aux_order=-1,
                                      allow_continuation=true, hf=cfg)
                # -Γ(1+eps)Γ(eps)/Γ(1+2eps) = -1/eps + 0 + ζ₂ eps - 2ζ₃ eps²
                E = SubTropica.eval_symbolic(L, NamedTuple(); dps=40)
                co(k) = E.coeffs[k - E.minorder + 1]
                @test E.minorder <= -1
                @test abs(co(-1) + 1) < 1e-30
                @test abs(co(0)) < 1e-30
                @test abs(co(1) - z2) < 1e-25
                @test abs(co(2) + 2 * z3) < 1e-25
                @test L.evidence["b2q"]["continuation"] == true
                @test L.evidence["b2q"]["aux_order"] == -1
            end

            @testset "mini end-to-end [q^0] (v=1 insertion route)" begin
                L = subtropica_integrate(mkmini(); order=2, aux_order=0,
                                      allow_continuation=true, hf=cfg)
                # -ζ₂ + ζ₃ eps + (π⁴/60) eps²
                E = SubTropica.eval_symbolic(L, NamedTuple(); dps=40)
                co(k) = E.coeffs[k - E.minorder + 1]
                for k in E.minorder:-1
                    @test abs(co(k)) < 1e-30            # explicit-zero poles only
                end
                @test abs(co(0) + z2) < 1e-25
                @test abs(co(1) - z3) < 1e-25
                @test abs(co(2) - pi4_60) < 1e-25
                # provenance ledger: exact-convolution route + census + steps
                b2q = L.evidence["b2q"]
                @test b2q["route"] == "exact-q-convolution (MG7 decision, VERDICT_EQ411 §3)"
                @test length(b2q["census"]) == 3
                @test b2q["trigger"]["unregulated_by_eps_rays"] != []
                @test length(b2q["terms"]) == 1
                @test b2q["terms"][1]["q_pole_order"] == 1
                @test b2q["terms"][1]["insertion_orders"] == [1]
            end
        end
    end
end

println("test_tworeg_driver: all testsets finished")
