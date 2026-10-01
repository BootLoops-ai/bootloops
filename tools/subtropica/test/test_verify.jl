# test_verify.jl — OWNER: unit-4. Self-contained tests for src/verify.jl +
# the fixtures/ set (no cross-unit deps: units 1-3 surfaces are not called).
# Run standalone:
#   ulimit -v 32505856; julia --project=tools/subtropica \
#     tools/subtropica/test/test_verify.jl
# or via test/runtests.jl (auto-included).
#
# Core self-test (task spec): synthetic-truth fixture vs OWN mpmath oracle
# artifact must PASS; a mutate_coefficient scratch copy must FAIL.

using Test
using SubTropica
using Nemo
using JSON

const _VF_FIXDIR = normpath(joinpath(@__DIR__, "..", "fixtures"))
const _VF_SCRATCH = get(ENV, "SUBTROPICA_TEST_SCRATCH", "") == "" ?
    mktempdir(; prefix="subtropica_verify_test_") :
    joinpath(ENV["SUBTROPICA_TEST_SCRATCH"], "verify_test_$(getpid())")
mkpath(_VF_SCRATCH)

_msg_of(f) = try
    f()
    ""
catch e
    sprint(showerror, e)
end

@testset "verify.jl (unit-4)" begin

    # ---------------------------------------------------------------- fixtures
    @testset "fixture inventory + schema" begin
        expected = Dict(
            "moller_box_eq41.json" => false,
            "smirnov_tst2_5var.json" => false,
            "synth_log2_1var.json" => false,
            "synth_zeta2_2var.json" => false,
            "synth_zeta3_beta_1var.json" => false,
            "divergent_xinv_1var.json" => true,
        )
        for (f, div) in expected
            path = joinpath(_VF_FIXDIR, f)
            @test isfile(path)
            fx = load_fixture(path)
            @test fx["schema"] == "subtropica-euler-quad-v1"
            @test fx["divergent"] == div
            @test haskey(fx["provenance"], "generated_by")
            # fixtures carry NO result values (RT-iii, structural)
            @test !haskey(fx, "laurent") && !haskey(fx, "value") &&
                  !haskey(fx, "result")
        end
        # Moller box specifics (DESIGN.md C2 conventions)
        box = load_fixture(joinpath(_VF_FIXDIR, "moller_box_eq41.json"))
        @test box["vars"] == ["x0", "x1", "x2"] &&
              box["kinvars"] == ["sb", "tb", "msq"]
        @test box["prefactor"]["gammas"] == [Dict("arg" => ["2", "1"], "power" => 1)]
        @test box["polys"][1]["exp"] == ["0", "2"]    # U^(2 eps)
        @test box["polys"][2]["exp"] == ["-2", "-1"]  # F^(-2-eps)
        # Smirnov twist exponents (verified symbolically at generation time)
        sm = load_fixture(joinpath(_VF_FIXDIR, "smirnov_tst2_5var.json"))
        @test sm["nu"] == [["0","3"], ["2","5"], ["4","-87"], ["2","-12"], ["1","19"]]
        @test length(sm["polys"]) == 5 && isempty(sm["kinvars"])
        bad = joinpath(_VF_SCRATCH, "not_a_fixture.json")
        write(bad, "{\"schema\": \"something-else\", \"id\": \"x\"}\n")
        @test_throws ErrorException load_fixture(bad)
    end

    # ------------------------------------------------- structural enforcement
    @testset "oracle artifacts: structural RT-iii enforcement" begin
        art_path = joinpath(_VF_FIXDIR, "oracle", "synth_log2_1var.mpmath.json")
        art = load_oracle_artifact(art_path)
        @test art.oracle === :mpmath && !art.sanity_only
        @test art.fixture_id == "synth_log2_1var"
        @test art.minorder == 0 && haskey(art.values, 0)
        @test art.claimed_digits[0] >= 30
        @test haskey(art.provenance, "method")

        # a FIXTURE is refused as a comparison-target source
        fix_path = joinpath(_VF_FIXDIR, "synth_log2_1var.json")
        @test_throws ErrorException load_oracle_artifact(fix_path)
        @test occursin("FORBIDDEN", _msg_of(() -> load_oracle_artifact(fix_path)))

        # stamp-less artifact refused
        raw = JSON.parsefile(art_path)
        stripped = Dict{String,Any}(String(k) => v for (k, v) in raw)
        stripped["provenance"] = Dict{String,Any}("note" => "no stamp")
        p1 = joinpath(_VF_SCRATCH, "stampless.json")
        write(p1, JSON.json(stripped))
        @test_throws ErrorException load_oracle_artifact(p1)
        @test occursin("provenance stamp", _msg_of(() -> load_oracle_artifact(p1)))

        # unknown oracle refused (FORM5 is NOT an oracle — RT-iii)
        raw2 = Dict{String,Any}(String(k) => v for (k, v) in JSON.parsefile(art_path))
        raw2["oracle"] = "form5"
        p2 = joinpath(_VF_SCRATCH, "form5.json")
        write(p2, JSON.json(raw2))
        @test_throws ErrorException load_oracle_artifact(p2)

        # pysecdec forced sanity tier: loads, but its reports can never gate
        raw3 = Dict{String,Any}(String(k) => v for (k, v) in JSON.parsefile(art_path))
        raw3["oracle"] = "pysecdec"
        p3 = joinpath(_VF_SCRATCH, "pysecdec.json")
        write(p3, JSON.json(raw3))
        art3 = load_oracle_artifact(p3)
        @test art3.sanity_only
        ctrl = synthetic_controls()[1]
        rep3 = compare_laurent(ctrl.candidate(), art3)
        @test rep3.sanity_only && !can_gate(rep3)
    end

    # ------------------------------------------------------------ floored digits
    @testset "floored_matched_digits" begin
        setprecision(BigFloat, 384) do
            a = BigFloat(1)
            @test floored_matched_digits(a, a; cap=60) == 60
            @test floored_matched_digits(big"0.0", big"0.0"; cap=45) == 45
            @test floored_matched_digits(a, a + parse(BigFloat, "1e-10")) == 10
            @test floored_matched_digits(a, -a) == 0
            @test floored_matched_digits(a, a * (1 + parse(BigFloat, "1e-25") *
                                                 BigFloat(pi))) in 24:25
        end
    end

    # ------------------------------------------- self-test: positive controls
    @testset "synthetic-truth self-test PASSES (Arb candidate vs mpmath artifact)" begin
        for ctrl in synthetic_controls()
            cand = ctrl.candidate()
            rep = verify_fixture(ctrl.fixture, cand, ctrl.artifact; min_digits=30)
            @test rep.pass
            @test can_gate(rep)
            @test rep.oracle === :mpmath
            @test rep.details["worst_digits"] >= 30
            # honest ceiling: cannot claim more than the artifact provides
            @test all(d -> d <= 60, values(rep.digits))
        end
        # multi-order fixture really compared 4 orders (pole-guard coverage)
        rep = verify_fixture(synthetic_controls()[3].fixture,
                             synthetic_controls()[3].candidate(),
                             synthetic_controls()[3].artifact)
        @test sort(collect(keys(rep.digits))) == [0, 1, 2, 3]
    end

    # ---------------------------------------- self-test: negative control FAILS
    @testset "mutation negative control FAILS (scratch copy, never in place)" begin
        ctrl = synthetic_controls()[1]
        cand = ctrl.candidate()
        orig_coeffs = copy(cand.coeffs)

        mut = mutate_coefficient(cand)                    # in-memory scratch copy
        @test mut !== cand && mut.coeffs !== cand.coeffs
        @test cand.coeffs == orig_coeffs                  # source untouched
        @test mut.coeffs[1] != cand.coeffs[1]
        @test haskey(mut.evidence, "mutated")
        rep = verify_fixture(ctrl.fixture, mut, ctrl.artifact; min_digits=30)
        @test !rep.pass && !can_gate(rep)
        # failed for the calibrated reason: ~24 digits (1e-25*pi), not garbage
        @test 20 <= rep.details["worst_digits"] <= 27

        # multi-order flavor: mutate a HIGHER order, worst-case must still govern
        ctrl3 = synthetic_controls()[3]
        mut3 = mutate_coefficient(ctrl3.candidate(); order=2)
        rep3 = verify_fixture(ctrl3.fixture, mut3, ctrl3.artifact; min_digits=30)
        @test !rep3.pass
        @test rep3.details["worst_order"] == 2

        # file flavor: mutates a SCRATCH COPY of the artifact, source untouched
        art_path = ctrl.artifact
        bytes_before = read(art_path)
        mpath = mutate_coefficient(art_path, joinpath(_VF_SCRATCH, "mut"))
        @test isfile(mpath) && dirname(mpath) != dirname(art_path)
        @test read(art_path) == bytes_before
        mart = load_oracle_artifact(mpath)
        @test haskey(mart.provenance, "mutated")
        repf = compare_laurent(ctrl.candidate(), mart; min_digits=30)
        @test !repf.pass
        @test haskey(repf.details, "mutated_artifact")
        # in-place mutation is REFUSED
        @test_throws ErrorException mutate_coefficient(art_path, dirname(art_path))
    end

    # ------------------------------------------------------- divergent refusal
    @testset "divergent fixture refused LOUDLY (gate A iv)" begin
        div_path = joinpath(_VF_FIXDIR, "divergent_xinv_1var.json")
        @test load_fixture(div_path)["divergent"] === true
        ctrl = synthetic_controls()[1]
        @test_throws ErrorException verify_fixture(div_path, ctrl.candidate(),
                                                   ctrl.artifact)
        @test occursin("REFUSING divergent",
                       _msg_of(() -> verify_fixture(div_path, ctrl.candidate(),
                                                    ctrl.artifact)))
    end

    # -------------------------------------------------------- pole-guard strictness
    @testset "order-set mismatch fails (truncated-Laurent pole guard)" begin
        ctrl = synthetic_controls()[3]                # 4-order artifact
        art = load_oracle_artifact(ctrl.artifact)
        short = ctrl.candidate()
        truncated = LaurentSeries(0, short.coeffs[1:2])   # drops eps^2, eps^3
        rep = compare_laurent(truncated, art)
        @test !rep.pass && haskey(rep.details, "order_mismatch")
        shifted = LaurentSeries(-1, short.coeffs)         # fake pole shift
        rep2 = compare_laurent(shifted, art)
        @test !rep2.pass && haskey(rep2.details, "order_mismatch")
    end

    # ------------------------------------------------------------ NOT-RUN stubs
    @testset "oracle adapters: AMFlow/hiprec NOT-RUN stubs, pysecdec refused" begin
        fx = joinpath(_VF_FIXDIR, "synth_log2_1var.json")
        for (oracle, marker) in ((:amflow, "NOT-RUN"),
                                 (:hiprec_sectordecomp, "NOT-RUN"),
                                 (:pysecdec, "SANITY-ONLY"))
            @test_throws ErrorException run_oracle(oracle, fx)
            @test occursin(marker, _msg_of(() -> run_oracle(oracle, fx)))
        end
        @test_throws ErrorException run_oracle(:form5, fx)
        # symbol-only verify_laurent stays a loud stub (bring an artifact)
        L = LaurentSeries(0, [HlogExpr()])
        @test_throws ErrorException verify_laurent(L, :amflow)
    end

    # ----------------------------------------------------------- point picker
    @testset "pick_point: deterministic, positive, degeneracy-avoiding" begin
        R, (sb, tb) = polynomial_ring(QQ, ["sb", "tb"])
        letters = [sb - tb, sb * tb - 1, 3 * sb - 1]
        p1 = pick_point([:sb, :tb, :msq]; letters=letters, seed=42)
        p2 = pick_point([:sb, :tb, :msq]; letters=letters, seed=42)
        @test p1 == p2                                    # deterministic
        p3 = pick_point([:sb, :tb, :msq]; letters=letters, seed=43)
        @test p3 != p1                                    # seed-sensitive
        for k in (:sb, :tb, :msq)
            v = p1[k]
            @test v isa Rational{BigInt} && v > 0
            @test v != 1 && v != 1 // 3 && denominator(v) > 1
        end
        @test length(unique(collect(values(p1)))) == 3    # pairwise distinct
        for L in letters                                  # away from degeneracies
            @test !iszero(evaluate(L, [QQ(p1[:sb]), QQ(p1[:tb])]))
        end
        pneg = pick_point([:s, :t, :msq]; signs=Dict(:s => -1, :t => -1), seed=7)
        @test pneg[:s] < 0 && pneg[:t] < 0 && pneg[:msq] > 0
        @test_throws ErrorException pick_point([:z]; letters=[sb - tb])  # unassigned var
    end

    # ---------------------------------------------------------- misc hygiene
    @testset "parse_rational exactness" begin
        @test parse_rational("22/7") == 22 // 7
        @test parse_rational("-3") == -3 // 1
        @test_throws ErrorException parse_rational("0.5")     # dyadic footgun
        @test_throws ErrorException parse_rational("1e-3")
    end

    # -------------------------------------------- live mpmath oracle (guarded)
    @testset "live run_oracle(:mpmath) round trip" begin
        can_run = Sys.which("python3") !== nothing &&
                  success(`bash -c "ulimit -v 32505856; python3 -c 'import mpmath, sympy'"`)
        if can_run
            fx = joinpath(_VF_FIXDIR, "synth_log2_1var.json")
            outdir = joinpath(_VF_SCRATCH, "live_oracle")
            apath = run_oracle(:mpmath, fx; dps=40, outdir=outdir)
            @test isfile(apath) && startswith(normpath(apath), normpath(outdir))
            art = load_oracle_artifact(apath)
            @test art.fixture_id == "synth_log2_1var" && art.oracle === :mpmath
            rep = compare_laurent(synthetic_controls()[1].candidate(), art;
                                  min_digits=30)
            @test rep.pass && can_gate(rep)
            # divergent fixture refused by the live adapter too
            div = joinpath(_VF_FIXDIR, "divergent_xinv_1var.json")
            @test_throws Exception run_oracle(:mpmath, div; dps=30, outdir=outdir)
        else
            @info "python3+mpmath+sympy unavailable — live :mpmath adapter test skipped"
            @test true
        end
    end

    # ------------------------------------- eval_symbolic 
    @testset "eval_symbolic: symbolic Laurent -> BigFloat at a point" begin
        # (1) Smirnov gate-(b) exact result: 1 - 4/5*mzv_2^2 + mzv_3
        #     == 1 - 2 zeta(4) + zeta(3); reference via Nemo/Arb — a chain
        #     independent of the ginac_gpl evaluator under test.
        if isfile(SubTropica._ginac_bin())
            RR = SubTropica.ArbField(320)
            Lsm = LaurentSeries(0, [SubTropica.parse_coef(
                      "1 - 4/5*mzv_2^2 + mzv_3", Symbol[])])
            N = SubTropica.eval_symbolic(Lsm, Dict{Symbol,Rational{BigInt}}(); dps=60)
            @test N isa LaurentSeries{BigFloat} && N.minorder == 0
            refsm = setprecision(BigFloat, 340) do
                BigFloat(RR(1) - 2*SubTropica.zeta(RR(4)) + SubTropica.zeta(RR(3)))
            end
            @test SubTropica._matched_digits(N.coeffs[1], refsm) >= 55
            # (2) synth beta fixture orders 0..3 vs the Arb synthetic-truth
            #     candidates (verify.jl synthetic_controls, independent build)
            Lbeta = LaurentSeries(0, [SubTropica.parse_coef(c, Symbol[]) for c in
                        ("1", "-3", "9 - 2*mzv_2", "-27 + 6*mzv_2 + 6*mzv_3")])
            Nbeta = SubTropica.eval_symbolic(Lbeta, Dict{Symbol,Rational{BigInt}}(); dps=60)
            cand = SubTropica.synthetic_controls(prec=320)[3].candidate()
            for k in 1:4
                @test SubTropica._matched_digits(Nbeta.coeffs[k], cand.coeffs[k]) >= 55
            end
            # (3) kinematic coef + Hlog atom: (s/(s+t)) * G(0,-1;1) at
            #     s=3/7, t=2/5; G(0,-1;1) = -Li2(-1) = pi^2/12
            ek = HlogExpr([SubTropica.HTerm("s/(s + t)",
                     [SubTropica.HHlog("1", ["0", "-1"]) => 1])], Symbol[:s, :t])
            Lk = LaurentSeries(0, [ek])
            Nk = SubTropica.eval_symbolic(Lk, Dict(:s => 3//7, :t => 2//5); dps=60)
            refk = setprecision(BigFloat, 340) do
                (BigFloat(15) / 29) * BigFloat(pi)^2 / 12   # (3/7)/(3/7+2/5) = 15/29
            end
            @test SubTropica._matched_digits(Nk.coeffs[1], refk) >= 55
            # evidence stamp carries the point
            @test Nk.evidence["eval_symbolic_point"]["s"] == "3/7"
            # (3b) word-level numeric leg: a residual
            #      rational-letter HPeriod atom (resolve_periods strict=false
            #      explicit-atom path) EVALUATES via the independent
            #      zero_inf_period_ginac route — ZeroInfPeriod[{0,-1}] = zeta_2
            #      (live v1.2.8 pin, CONTRACTS.md §a). The still-refused classes
            #      are in (4).
            Lper = LaurentSeries(0, [HlogExpr([SubTropica.HTerm("1",
                       [SubTropica.HPeriod(["0", "-1"]) => 1])], Symbol[])])
            Nper = SubTropica.eval_symbolic(Lper, Dict{Symbol,Rational{BigInt}}();
                                         dps=60)
            refz2 = setprecision(BigFloat, 340) do
                BigFloat(pi)^2 / 6
            end
            @test SubTropica._matched_digits(Nper.coeffs[1], refz2) >= 55
        else
            @info "ginac_gpl CLI unavailable — eval_symbolic numeric tests skipped"
            @test true
        end
        # (4) typed refusals (need no ginac): algebraic-letter period word
        #     with NO table entry; HAlgLetter atom with no table entry;
        #     missing kinematic symbol; numeric-input misuse
        Lwm = LaurentSeries(0, [HlogExpr([SubTropica.HTerm("1",
                  [SubTropica.HPeriod(["0", "Wm_9"]) => 1])], Symbol[])])
        @test_throws SubTropica.SubTropicaEvalError SubTropica.eval_symbolic(
            Lwm, Dict{Symbol,Rational{BigInt}}())
        Lal = LaurentSeries(0, [HlogExpr([SubTropica.HTerm("1",
                  [SubTropica.HAlgLetter(:sqrt_disc, 9) => 1])], Symbol[])])
        @test_throws SubTropica.SubTropicaEvalError SubTropica.eval_symbolic(
            Lal, Dict{Symbol,Rational{BigInt}}())
        Lmiss = LaurentSeries(0, [HlogExpr([SubTropica.HTerm("s", Pair{SubTropica.HAtom,Int}[])],
                                           Symbol[:s])])
        @test_throws ErrorException SubTropica.eval_symbolic(
            Lmiss, Dict{Symbol,Rational{BigInt}}())
        @test_throws ErrorException SubTropica.eval_symbolic(
            LaurentSeries(0, [big"1.0"]), Dict{Symbol,Rational{BigInt}}())
    end

    # -------------------------------------- mzv_provenance_gate (A(iii)/R7)
    @testset "mzv_provenance_gate: refusals + live ginac cross-eval" begin
        # (a) missing table: typed refusal naming the env var (engine-free)
        nofile = joinpath(_VF_SCRATCH, "no_such_table.json")
        @test_throws ErrorException SubTropica.mzv_provenance_gate(table_path=nofile)
        @test occursin("SUBTROPICA_MZV_DATA",
                       _msg_of(() -> SubTropica.mzv_provenance_gate(table_path=nofile)))

        # (b) toy table — classical identities in the pinned ASCENDING
        # convention (mzv_bigfloat docstring): Euler zeta(1,2)=zeta(3);
        # zeta(2)=Pi^2/6; zeta(4)=2/5*zeta(2)^2; sum (-1)^n/n = -Log2.
        toy = joinpath(_VF_SCRATCH, "toy_mzv_reductions.json")
        toyrules = [Dict("lhs" => "mzv_1_2", "rhs" => "mzv_3"),
                    Dict("lhs" => "mzv_2",   "rhs" => "1/6*Pi^2"),
                    Dict("lhs" => "mzv_4",   "rhs" => "2/5*mzv_2^2"),
                    Dict("lhs" => "mzv_m1",  "rhs" => "-Log2")]
        write(toy, JSON.json(Dict("reductions" => toyrules,
                                  "basis" => ["Log2", "mzv_2", "mzv_3"])))

        # (c) pin mismatch refuses BEFORE evaluation; malformed pin refuses
        @test occursin("MISMATCH", _msg_of(() -> SubTropica.mzv_provenance_gate(
            table_path=toy, expect_sha256="0"^64)))
        @test_throws ErrorException SubTropica.mzv_provenance_gate(
            table_path=toy, expect_sha256="not-a-sha")

        # (d) dps floor: below 30d the mutation control cannot bite — refused
        @test_throws ErrorException SubTropica.mzv_provenance_gate(
            table_path=toy, dps=10)

        # (e) divergent-LHS-only table: loud refusal, not a silent no-op
        divtab = joinpath(_VF_SCRATCH, "toy_mzv_divergent.json")
        write(divtab, JSON.json(Dict(
            "reductions" => [Dict("lhs" => "mzv_2_1", "rhs" => "mzv_3")],
            "basis" => String[])))
        @test occursin("divergent",
                       _msg_of(() -> SubTropica.mzv_provenance_gate(table_path=divtab)))

        if isfile(SubTropica._ginac_bin())
            # (f) live gate on the toy table: 4/4 rules >= 40d via ginac
            rep = SubTropica.mzv_provenance_gate(table_path=toy; nrules=8, dps=40)
            @test rep isa SubTropica.GateReport
            @test rep.pass && can_gate(rep)
            @test rep.oracle === :ginac
            @test rep.details["nrules_checked"] == 4        # min(8, 4 rules)
            @test rep.details["worst_digits"] >= 40
            @test Set(rep.details["sampled_lhs"]) ==
                  Set(["mzv_1_2", "mzv_2", "mzv_4", "mzv_m1"])
            @test rep.details["table_sha256_pinned"] == "UNPINNED"
            @test occursin(r"^[0-9a-f]{64}$", rep.details["table_sha256"])
            # negative control fired at the calibrated 1e-25*pi magnitude
            @test 20 <= rep.details["negative_control_digits"] <= 27
            # correct pin (script-computed from the bytes, never hand-typed)
            pin = SubTropica._sha256_hex(read(toy))
            rep_p = SubTropica.mzv_provenance_gate(table_path=toy; nrules=2,
                                                   dps=40, expect_sha256=pin)
            @test rep_p.pass && rep_p.details["table_sha256_pinned"] == pin
            # deterministic sampling: same seed => same sample
            r1 = SubTropica.mzv_provenance_gate(table_path=toy; nrules=2,
                                                dps=40, seed=11)
            r2 = SubTropica.mzv_provenance_gate(table_path=toy; nrules=2,
                                                dps=40, seed=11)
            @test r1.details["sampled_lhs"] == r2.details["sampled_lhs"]

            # (g) doctored rule must FAIL (scratch copy — source untouched):
            # swap the Euler rhs to mzv_2, an O(1) mismatch
            bad = joinpath(_VF_SCRATCH, "toy_mzv_bad.json")
            badrules = deepcopy(toyrules)
            badrules[1]["rhs"] = "mzv_2"
            write(bad, JSON.json(Dict("reductions" => badrules,
                                      "basis" => ["Log2", "mzv_2", "mzv_3"])))
            repb = SubTropica.mzv_provenance_gate(table_path=bad; nrules=8, dps=40)
            @test !repb.pass && !can_gate(repb)
            @test repb.details["worst_lhs"] == "mzv_1_2"
            @test repb.details["rule_digits"]["mzv_1_2"] < 40
        else
            @info "ginac_gpl CLI unavailable — mzv_provenance_gate live legs skipped (SUBTROPICA_GINAC_GPL)"
            @test occursin("ginac_gpl",
                           _msg_of(() -> SubTropica.mzv_provenance_gate(table_path=toy)))
        end
    end
end
