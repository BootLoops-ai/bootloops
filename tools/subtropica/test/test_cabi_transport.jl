# test_cabi_transport.jl — C-ABI transport tests. The OPT-IN
# in-process C-ABI transport (hf_bridge.jl C-ABI section; deps/BUILD.md):
#   (1) CLI-vs-ccall BYTE-IDENTITY on 5 identity inputs (envelope+timing
#       strip = the repo's own snapshot-ctest normalization),
#   (2) the F2 pair sequence in THIS process (the unpatched-lib killer),
#   (3) a 100-call heterogeneous in-process soak (byte-stable vs first
#       response — the F2 regression gate shape),
#   (4) hf_call end-to-end parity + typed refusals (transport law).
# Self-contained: run via test/runtests.jl (auto-included) or
#   ulimit -v 32505856; julia --project=tools/subtropica tools/subtropica/test/test_cabi_transport.jl
# Skips cleanly if the vendored lib or the CLI engine is absent.
# LAW REMINDER: the CLI stays the default transport + validation oracle;
# threads sharing the ABI handle are FORBIDDEN (measured crash inside FLINT)
# — nothing here spawns threads.

using Test
using SubTropica
using JSON

@testset "C-ABI transport (opt-in ccall farm route)" begin

    cfg = SubTropica.HFConfig()                    # transport=:cli default
    have_lib = isfile(cfg.cabi_lib)
    have_cli = isfile(cfg.bin) && isfile(cfg.mzv_data_path)
    have_lib || @warn "C-ABI lib not found at $(cfg.cabi_lib) — cabi testsets skipped"
    have_cli || @warn "HyperFLINT CLI not found at $(cfg.bin) — cabi testsets skipped"

    @testset "transport selection law" begin
        # default stays :cli — the validation oracle (never silently :cabi)
        @test cfg.transport === :cli ||
              get(ENV, "SUBTROPICA_HF_TRANSPORT", "cli") != "cli"
        # unknown transport refuses typed, pre-spawn
        @test_throws SubTropica.HFError SubTropica._hf_transport(
            Dict{String,Any}("op" => "hyperflint"), cfg;
            timeout_s=1, transport=:carrier_pigeon)
        # non-ABI op on the raw cabi leg refuses typed (G1 gap named)
        @test_throws SubTropica.HFError SubTropica._cabi_call_raw(
            "parse_expr", "{}"; cfg=cfg)
        # bad lib path refuses typed at the dlopen gate
        badcfg = SubTropica.HFConfig(cabi_lib="/nonexistent/libhf.so")
        @test_throws SubTropica.HFGateError SubTropica._cabi_call_raw(
            "partial_fractions", "{}"; cfg=badcfg)
    end

    if have_lib && have_cli
        MZV = cfg.mzv_data_path
        mk_hf(expr, vi, vs) = "{\"op\":\"hyperflint\",\"expr\":\"$expr\"," *
            "\"vars_int\":$vi,\"vars\":$vs,\"mzv_data_path\":\"$MZV\"," *
            "\"algebraic_letters\":false}"

        # THE 5 identity inputs (probe set: 3 pfrac incl. unary-minus canary,
        # 1 linear_factors, 1 full integrator op).
        IDENT = [
            ("partial_fractions",
             """{"op":"partial_fractions","f":"1/(x*(1 + x)*(x - y))","vars":["x","y"],"var":"x"}"""),
            ("partial_fractions",
             """{"op":"partial_fractions","f":"(x^2 + y)/((x + y)^2*(1 - x))","vars":["x","y"],"var":"x"}"""),
            ("partial_fractions",   # unary-minus canary through BOTH parsers
             """{"op":"partial_fractions","f":"-x^2","vars":["x"],"var":"x"}"""),
            ("linear_factors",
             """{"op":"linear_factors","poly":"2*x^2 + 3*x*y + y^2","var":"x","vars":["x","y"]}"""),
            ("hyperflint", mk_hf("1/((1+x)*(2+x))", "[\"x\"]", "[\"x\"]")),
        ]

        @testset "byte-identity CLI vs ccall (5 identity inputs)" begin
            for (op, req) in IDENT
                abi_resp = SubTropica._cabi_call_raw(op, req; cfg=cfg)
                cli_resp = SubTropica._cli_call_raw(req; cfg=cfg, timeout_s=300)
                @test SubTropica._cabi_normalize(abi_resp) ==
                      SubTropica._cabi_normalize(cli_resp)
                @test !occursin("\"error\"", abi_resp)
            end
        end

        @testset "F2 pair sequence, one process" begin
            # 1-var pfrac -> 2-var pfrac -> identical 2-var pfrac; the
            # unpatched lib asserts (Debug) / aborts-segvs (Release) here.
            p1 = """{"op":"partial_fractions","f":"1/(x*(1 + x))","vars":["x"],"var":"x"}"""
            p2 = """{"op":"partial_fractions","f":"1/(x*(1 + x)*(x - y))","vars":["x","y"],"var":"x"}"""
            r1 = SubTropica._cabi_call_raw("partial_fractions", p1; cfg=cfg)
            r2 = SubTropica._cabi_call_raw("partial_fractions", p2; cfg=cfg)
            r3 = SubTropica._cabi_call_raw("partial_fractions", p2; cfg=cfg)
            @test !occursin("\"error\"", r1)
            @test !occursin("\"error\"", r2)
            @test r2 == r3                       # byte-identical repeat
            # linear_factors face of F2 (the get_str_pretty segv face)
            lf = """{"op":"linear_factors","poly":"2*x^2 + 3*x*y + y^2","var":"x","vars":["x","y"]}"""
            l1 = SubTropica._cabi_call_raw("linear_factors", lf; cfg=cfg)
            l2 = SubTropica._cabi_call_raw("linear_factors", lf; cfg=cfg)
            @test !occursin("\"error\"", l1) && l1 == l2
        end

        @testset "100-call heterogeneous in-process soak" begin
            # Mixed op shapes interleaved (heap-churn pattern — the trigger
            # surface of the original F2 crash); every response must be
            # byte-identical (after timing strip) to its first occurrence.
            reqs = [IDENT[1], IDENT[2], IDENT[3], IDENT[4], IDENT[5],
                    ("hyperflint", mk_hf("1/((1+x)*(1+2*x))", "[\"x\"]", "[\"x\"]")),
                    ("hyperflint", mk_hf("x/((1+x)^2*(2+x))", "[\"x\"]", "[\"x\"]"))]
            ref = Dict{Int,String}()
            mismatches = 0
            for i in 1:100
                j = (i - 1) % length(reqs) + 1
                (op, req) = reqs[j]
                s = SubTropica._cabi_normalize(
                        SubTropica._cabi_call_raw(op, req; cfg=cfg))
                if haskey(ref, j)
                    s == ref[j] || (mismatches += 1)
                else
                    ref[j] = s
                    @test !occursin("\"error\"", s)
                end
            end
            @test mismatches == 0
        end

        @testset "hf_call end-to-end parity + taxonomy on :cabi" begin
            fields = Dict{String,Any}(
                "expr" => "1/((1+x)*(1+y)*(1+x+y))",
                "vars_int" => ["x", "y"], "vars" => ["x", "y"],
                "check_divergences" => true, "algebraic_letters" => false,
                "canonical_emission" => true)
            ra = SubTropica.hf_call("hyperflint", fields; cfg=cfg, transport=:cabi)
            rc = SubTropica.hf_call("hyperflint", fields; cfg=cfg, transport=:cli)
            @test JSON.json(ra["result"]) == JSON.json(rc["result"])
            @test ra["hf_version"] == "1.2.8"    # envelope double check
            # in-band divergence taxonomy IDENTICAL on the cabi transport:
            dfields = Dict{String,Any}(
                "expr" => "1/(x*(1+x))", "vars_int" => ["x"], "vars" => ["x"],
                "check_divergences" => true)
            @test_throws SubTropica.HFDivergent SubTropica.hf_call(
                "hyperflint", dfields; cfg=cfg, transport=:cabi)
        end
    end
end
