# runtests.jl — SHARED SKELETON (shared structure; each unit adds
# tests by DROPPING test/test_*.jl files — this file auto-includes them,
# NO edits here needed). Run:
#   ulimit -v 32505856; julia --project=tools/subtropica \
#     tools/subtropica/test/runtests.jl

using Test
using SubTropica

@testset "SubTropica" begin

    @testset "scaffold: shared types construct" begin
        e = SubTropica.EpsExp(1//2, -3//4)
        @test e.a == 1//2 && e.b == -3//4
        @test SubTropica.iseps(e)
        @test !SubTropica.iseps(SubTropica.EpsExp(2))
        pf = SubTropica.Prefactor(1//1)
        @test pf.eps_power == 0 && isempty(pf.gammas)

        t = SubTropica.HTerm("1/2", [SubTropica.HConst(:Pi) => 2,
                                  SubTropica.HZeta([1, -3]) => 1])
        h = SubTropica.HlogExpr([t], Symbol[:s, :t])
        @test length(h.terms) == 1
        @test SubTropica.HZeta([1, -3]) == SubTropica.HZeta([1, -3])   # hash/== sanity
        @test SubTropica.HConst(:Pi) != SubTropica.HConst(:I)

        L = SubTropica.LaurentSeries(-2, [SubTropica.HlogExpr(), SubTropica.HlogExpr()])
        @test L.minorder == -2 && L.provenance == :analytic

        cfg = SubTropica.HFConfig()
        @test endswith(cfg.bin, "hyperflint.sh")
        @test endswith(cfg.mzv_data_path, "mzv_reductions.json")
    end

    @testset "scaffold: engine present (read-only)" begin
        cfg = SubTropica.HFConfig()
        # Presence-only checks; behavioral stamp+canary gate is unit-1's
        # hf_gate (test_bridge.jl). Skip cleanly on boxes without the engine.
        if isfile(cfg.bin)
            @test isfile(cfg.mzv_data_path)
        else
            @info "HyperFLINT engine not found at $(cfg.bin) — bridge tests will skip"
        end
    end

    @testset "scaffold: API stubs refuse loudly" begin
        L = SubTropica.LaurentSeries(0, [SubTropica.HlogExpr()])
        @test_throws ErrorException SubTropica.verify_laurent(L, :amflow)
        # loggamma_series is IMPLEMENTED (unit-2, eps_expand.jl); behavior
        # covered in test_expand.jl. Sanity: it returns without throwing.
        @test SubTropica.loggamma_series(SubTropica.EpsExp(1), 3) !== nothing
        @test_throws ErrorException SubTropica.hf_call("parse_expr", Dict{String,Any}())
    end

    # ---- unit test files: auto-include every test/test_*.jl ------------
    # Each file runs in its OWN sandbox module:
    # the B1/C test files are standalone-first and include() src files by
    # path (types.jl structs included), which cannot coexist with
    # `using SubTropica` in a shared Main. Test registration is unaffected
    # (the @testset stack is task-local; module boundaries don't break it).
    testdir = @__DIR__
    for f in sort(filter(f -> startswith(f, "test_") && endswith(f, ".jl"),
                         readdir(testdir)))
        modname = Symbol("SubTropicaSuite_" * replace(f, r"\.jl$" => ""))
        path = joinpath(testdir, f)
        @testset "$f" begin
            Core.eval(Main, Expr(:module, true, modname,
                                 Expr(:block, :(include($path)))))
        end
    end

end
