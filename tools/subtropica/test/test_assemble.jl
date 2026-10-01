# test_assemble.jl — OWNER: unit-3. Tests for src/lr_dispatch.jl +
# src/assemble.jl. Self-contained: NO dependence on other units' code or
# tests — the HF transport is injected (canned-JSON mocks copied verbatim
# from live v1.2.8 probes), plus ONE live smoke assembly on a
# trivial 1-d integral (skipped cleanly if the engine is absent).
# Standalone run:
#   ulimit -v 32505856; julia --project=tools/subtropica \
#     tools/subtropica/test/test_assemble.jl

using Test
using SubTropica
using JSON

# --------------------------------------------------------------------------
# helpers (a3_ prefix — runtests.jl includes every unit's tests into Main)
# --------------------------------------------------------------------------

# Mock transport: dispatch function gets (op, fields) and returns the raw
# JSON string of a canned response; requests are captured for assertions.
function a3_mock(dispatch::Function; capture::Vector{Any}=Any[])
    return function (op::String, fields::Dict; cfg::SubTropica.HFConfig=SubTropica.HFConfig())
        push!(capture, (op, deepcopy(fields)))
        return JSON.parse(dispatch(op, fields))
    end
end

# Live transport (test-local; unit-1's hf_call is a stub at scaffold time
# and cross-unit test deps are forbidden). CONTRACTS.md §a transport:
# one flat JSON object on stdin, parse stdout's LAST line only, stderr
# counters are normal, never gate on rc. The `ulimit -v` prefix used
# throughout this package is applied inside the spawned shell.
function a3_live_transport(op::String, fields::Dict;
                           cfg::SubTropica.HFConfig=SubTropica.HFConfig())
    req = Dict{String,Any}("op" => op)
    merge!(req, Dict{String,Any}(String(k) => v for (k, v) in fields))
    haskey(req, "mzv_data_path") || (req["mzv_data_path"] = cfg.mzv_data_path)
    cmd = `bash -c "ulimit -v 32505856; exec \"$(cfg.bin)\" eval-json"`
    out = IOBuffer()
    p = run(pipeline(cmd; stdin=IOBuffer(JSON.json(req)), stdout=out,
                     stderr=devnull); wait=true)
    lines = filter(!isempty, split(String(take!(out)), '\n'))
    isempty(lines) && error("a3_live_transport: no stdout from HF (rc=$(p.exitcode))")
    return JSON.parse(String(last(lines)))
end

# Rational-coefficient HlogExpr builders/readers for assembly checks.
a3_hl(c::String, atoms::Vector{<:Pair}=Pair{SubTropica.HAtom,Int}[];
      vars::Vector{Symbol}=Symbol[]) =
    SubTropica.HlogExpr(
        [SubTropica.HTerm(c, Pair{SubTropica.HAtom,Int}[a => p for (a, p) in atoms])],
        vars)
a3_hl(r::Rational, args...; kw...) =
    a3_hl(SubTropica._a3_ratstr(Rational{BigInt}(r)), args...; kw...)

# Sum of plain-rational coefs of the atom-free part (fails test on opaque coef).
function a3_constval(e::SubTropica.HlogExpr)
    acc = 0 // big(1)
    for t in e.terms
        isempty(t.atoms) || continue
        r = SubTropica._a3_tryrat(t.coef)
        r === nothing && error("a3_constval: non-rational coef $(t.coef)")
        acc += r
    end
    return acc
end

# Dict atom-part-string => rational coef (for atom-carrying checks).
function a3_coefdict(e::SubTropica.HlogExpr)
    d = Dict{String,Rational{BigInt}}()
    for t in e.terms
        key = join((string(SubTropica._a3_atom_key(a), "^", p) for (a, p) in
                    SubTropica._a3_norm_atoms(t.atoms)), "*")
        r = SubTropica._a3_tryrat(t.coef)
        r === nothing && error("a3_coefdict: non-rational coef $(t.coef)")
        d[key] = get(d, key, 0 // big(1)) + r
    end
    return d
end

a3_letter(idx, poly, var; lc="1", s="-1", pr="1", disc="5") =
    SubTropica.AlgLetter(idx, poly, var, lc, s, pr, disc)

const A3_CFG = SubTropica.HFConfig()
const A3_ENGINE = isfile(A3_CFG.bin)
const A3_RUNS = joinpath(get(ENV, "SUBTROPICA_RUNS", "subtropica_runs"), "unit3_selftest")

# Canned responses — VERBATIM from live v1.2.8 probes.
const A3_R_SEARCH_X = """{"op":"find_lr_orders","schema_version":2,"hf_version":"1.2.8","best_order":["x"],"score":1,"nolr":false,"strategy":"LR_NoOpt","timing_compute_s":9.0544e-05,"nXVars":1,"nGroups":1,"nPolys":[1]}"""
const A3_R_SEARCH_XY = """{"op":"find_lr_orders","schema_version":2,"hf_version":"1.2.8","best_order":["x","y"],"score":4.21914,"nolr":false,"strategy":"LR_NoOpt","timing_compute_s":0.000183326,"nXVars":2,"nGroups":2,"nPolys":[1,1]}"""
const A3_R_NOLR = """{"op":"find_lr_orders","schema_version":2,"hf_version":"1.2.8","best_order":[],"score":null,"nolr":true,"strategy":"Fubini_Lungo","timing_compute_s":5.9094e-05,"nXVars":1,"nGroups":1,"nPolys":[1]}"""
const A3_R_VERIFY_OK = """{"op":"find_lr_orders","schema_version":2,"hf_version":"1.2.8","best_order":[],"score":null,"nolr":false,"strategy":"LR_NoOpt","timing_compute_s":1.2206e-05,"nXVars":1,"nGroups":1,"nPolys":[1],"order_is_lr":true,"verify_malformed":false,"verify_blocking_step":-1,"verify_blocking_degree":0,"verify_forbidden_dep":false,"verify_blocking_letter":""}"""
const A3_R_VERIFY_BAD = """{"op":"find_lr_orders","schema_version":2,"hf_version":"1.2.8","best_order":[],"score":null,"nolr":false,"strategy":"LR_NoOpt","timing_compute_s":2.0e-05,"nXVars":1,"nGroups":1,"nPolys":[1],"order_is_lr":false,"verify_malformed":false,"verify_blocking_step":0,"verify_blocking_degree":2,"verify_forbidden_dep":false,"verify_blocking_letter":"1+x+x^2"}"""
const A3_R_ERROR = """{"op":"find_lr_orders","error":"unknown op"}"""
# NOTE: nolr:false is INERT in verify mode (handlers.cpp:858) — the search
# envelope without order_is_lr must be rejected by verify_lr_order:
const A3_R_VERIFY_MISSING = A3_R_SEARCH_X

@testset "unit-3: lr_dispatch" begin

    @testset "find_lr_order: groups form, happy path" begin
        cap = Any[]
        t = a3_mock((op, f) -> A3_R_SEARCH_X; capture=cap)
        (order, meta) = SubTropica.find_lr_order([["1+x"]], [:x], String[];
                                              transport=t)
        @test order == [:x]
        @test meta["score"] == 1
        @test meta["strategy"] == "LR_NoOpt"
        @test meta["hf_version"] == "1.2.8"
        (op, fields) = cap[1]
        @test op == "find_lr_orders"
        @test haskey(fields, "groups")           # ALWAYS multi-group
        @test !haskey(fields, "polys")           # NEVER the single-group form
        @test fields["groups"] == [["1+x"]]
        @test fields["xvars"] == ["x"]
        @test !haskey(fields, "verify_order")
        @test !haskey(fields, "algebraic_letters")  # baseline envelope
    end

    @testset "find_lr_order: NOLR is a LOUD typed error" begin
        t = a3_mock((op, f) -> A3_R_NOLR)
        err = try
            SubTropica.find_lr_order([["1+x+x^2"]], [:x], String[]; transport=t)
            nothing
        catch e
            e
        end
        @test err isa SubTropica.NoLROrderError
        msg = sprint(showerror, err)
        @test occursin("REFUSING", msg)          # opposite of upstream silent-0
        @test occursin("group", msg)             # multi-group hint for ct-sums
    end

    @testset "find_lr_order: in-band error / malformed envelope" begin
        @test_throws SubTropica.LRDispatchError SubTropica.find_lr_order(
            [["1+x"]], [:x], String[]; transport=a3_mock((op, f) -> A3_R_ERROR))
        # caller bugs
        @test_throws SubTropica.LRDispatchError SubTropica.find_lr_order(
            Vector{String}[], [:x], String[]; transport=a3_mock((op, f) -> A3_R_SEARCH_X))
        @test_throws SubTropica.LRDispatchError SubTropica.find_lr_order(
            [String[]], [:x], String[]; transport=a3_mock((op, f) -> A3_R_SEARCH_X))
        @test_throws SubTropica.LRDispatchError SubTropica.find_lr_order(
            [["1+x"]], Symbol[], String[]; transport=a3_mock((op, f) -> A3_R_SEARCH_X))
    end

    @testset "verify_lr_order: verify_order FIELD replay" begin
        cap = Any[]
        t = a3_mock((op, f) -> A3_R_VERIFY_OK; capture=cap)
        @test SubTropica.verify_lr_order([:x], [["1+x"]], [:x], String[];
                                      transport=t)
        @test cap[1][2]["verify_order"] == ["x"]
        det = Dict{String,Any}()
        tb = a3_mock((op, f) -> A3_R_VERIFY_BAD)
        ok = SubTropica.verify_lr_order([:x], [["1+x+x^2"]], [:x], String[];
                                     transport=tb, details=det)
        @test ok == false
        @test det["verify_blocking_letter"] == "1+x+x^2"
        @test det["verify_blocking_degree"] == 2
        # response without order_is_lr (wrong binary / search envelope) throws
        @test_throws SubTropica.LRDispatchError SubTropica.verify_lr_order(
            [:x], [["1+x"]], [:x], String[];
            transport=a3_mock((op, f) -> A3_R_VERIFY_MISSING))
    end

    @testset "integrand_hash + persistence round-trip" begin
        h1 = SubTropica.integrand_hash([["1+x"]], [:x], String[])
        h2 = SubTropica.integrand_hash([["1+x"]], [:x], String[])
        @test h1 == h2                       # deterministic
        @test occursin(r"^[0-9a-f]{64}$", h1)
        @test h1 != SubTropica.integrand_hash([["1+x"], ["1+x+y"]], [:x], String[])
        @test h1 != SubTropica.integrand_hash([["1+x"]], [:y], String[])
        # xvar ORDER is part of the key
        ha = SubTropica.integrand_hash([["1+x+y"]], [:x, :y], String[])
        hb = SubTropica.integrand_hash([["1+x+y"]], [:y, :x], String[])
        @test ha != hb

        mkpath(A3_RUNS)
        dir = mktempdir(A3_RUNS)
        path = joinpath(dir, "best_order.json")
        meta = Dict{String,Any}("score" => 1, "integrand_hash" => h1)
        SubTropica.persist_best_order(path, [:x2, :x1], meta)
        (order, meta2) = SubTropica.load_best_order(path)
        @test order == [:x2, :x1]
        @test meta2["score"] == 1
        # corrupt file is a typed error
        write(path, "{not json")
        @test_throws SubTropica.LRReplayError SubTropica.load_best_order(path)
        @test_throws SubTropica.LRReplayError SubTropica.load_best_order(
            joinpath(dir, "nope.json"))
        rm(dir; recursive=true, force=true)
    end

    @testset "ensure_lr_order: search→persist, then verify-replay" begin
        mkpath(A3_RUNS)
        rundir = mktempdir(A3_RUNS)
        groups = [["1+x"]]
        # 1st call: fresh search + persist
        cap1 = Any[]
        t1 = a3_mock((op, f) -> haskey(f, "verify_order") ?
                         error("verify before any persisted state") :
                         A3_R_SEARCH_X; capture=cap1)
        (o1, m1) = SubTropica.ensure_lr_order(groups, [:x], String[];
                                           run_dir=rundir, transport=t1)
        @test o1 == [:x] && m1["replayed"] == false
        p = SubTropica.best_order_path(rundir, groups, [:x], String[])
        @test isfile(p)
        # 2nd call: MUST replay via verify_order field, never re-search
        t2 = a3_mock((op, f) -> haskey(f, "verify_order") ? A3_R_VERIFY_OK :
                         error("SEARCH dispatched on resume — replay bug"))
        (o2, m2) = SubTropica.ensure_lr_order(groups, [:x], String[];
                                           run_dir=rundir, transport=t2)
        @test o2 == [:x] && m2["replayed"] == true
        # stale persisted order: verify fails -> re-search + record change
        t3 = a3_mock((op, f) -> haskey(f, "verify_order") ? A3_R_VERIFY_BAD :
                         A3_R_SEARCH_XY)
        (o3, m3) = SubTropica.ensure_lr_order(groups, [:x], String[];
                                           run_dir=rundir, transport=t3)
        @test o3 == [:x, :y]                    # from the fresh search
        @test m3["order_changed_from"] == ["x"] # R3: gate must re-verify
        # on_stale=:error path
        SubTropica.persist_best_order(p, [:x],
            Dict{String,Any}("integrand_hash" =>
                SubTropica.integrand_hash(groups, [:x], String[])))
        t4 = a3_mock((op, f) -> A3_R_VERIFY_BAD)
        @test_throws SubTropica.LRReplayError SubTropica.ensure_lr_order(
            groups, [:x], String[]; run_dir=rundir, transport=t4,
            on_stale=:error)
        rm(rundir; recursive=true, force=true)
    end

    # ---- Phase B wrappers (find_lr_orders_scan / factor_table / gauge
    # derivation). Canned responses CONSTRUCTED from the reference handler
    # emitters (reference/hyperflint_bridge_handlers.cpp,
    # find_lr_orders_scan ~:1428-1454, factor_table ~:1146-1212) — shape-
    # pinned to the shipped source, not live probes; the live pin rides the
    # engine-armed suite.
    @testset "find_lr_orders_scan: request shape + response taxonomy" begin
        R_SCAN_OK = """{"op":"find_lr_orders_scan","schema_version":2,"hf_version":"1.2.8","projective":true,"truncated":false,"orders":[{"order":["x","y","z"],"gauge":"z","score":2.5,"carried_sqrts":0,"kin_sqrts":0,"terminal_quads":1},{"order":["y","x","z"],"gauge":"z","score":3.75,"carried_sqrts":1,"kin_sqrts":0,"terminal_quads":1}],"timing_compute_s":0.0123,"nXVars":3,"nGroups":2}"""
        R_SCAN_NONPROJ = """{"op":"find_lr_orders_scan","schema_version":2,"hf_version":"1.2.8","projective":false,"truncated":false,"orders":[],"timing_compute_s":0.001,"nXVars":3,"nGroups":2}"""
        R_SCAN_EMPTY = """{"op":"find_lr_orders_scan","schema_version":2,"hf_version":"1.2.8","projective":true,"truncated":false,"orders":[],"timing_compute_s":0.002,"nXVars":3,"nGroups":2}"""
        groups = [["1+x+y+z"], ["x*y+z"]]
        exps = [[[-2, -1]], [[1, 1]]]

        cap = Any[]
        t = a3_mock((op, f) -> R_SCAN_OK; capture=cap)
        (orders, meta) = SubTropica.find_lr_orders_scan(groups, [:x, :y, :z],
                                                        String[], exps;
                                                        transport=t)
        (op, fields) = cap[1]
        @test op == "find_lr_orders_scan"
        @test fields["groups"] == [["1+x+y+z"], ["x*y+z"]]
        @test fields["exps"] == [[[-2, -1]], [[1, 1]]]
        # baseline envelope: default knobs NOT armed
        @test !haskey(fields, "keep_rule") && !haskey(fields, "euler_filter") &&
              !haskey(fields, "max_orders")
        @test length(orders) == 2
        @test orders[1]["order"] == [:x, :y, :z] && orders[1]["gauge"] == :z
        @test orders[1]["score"] == 2.5 && orders[2]["carried_sqrts"] == 1
        @test meta["projective"] === true && meta["truncated"] === false
        @test meta["keep_rule"] == "Strict" && meta["nGroups"] == 2

        # non-default knobs armed on request
        cap2 = Any[]
        t2 = a3_mock((op, f) -> R_SCAN_OK; capture=cap2)
        SubTropica.find_lr_orders_scan(groups, [:x, :y, :z], String[], exps;
                                       keep_rule="FindRoots", euler_filter=true,
                                       max_orders=64, transport=t2)
        @test cap2[1][2]["keep_rule"] == "FindRoots"
        @test cap2[1][2]["euler_filter"] === true
        @test cap2[1][2]["max_orders"] == 64

        # client-side refusals (mirror the engine's loud checks)
        @test_throws SubTropica.LRDispatchError SubTropica.find_lr_orders_scan(
            groups, [:x, :y, :z], String[], exps; keep_rule="Speculative",
            transport=t)
        @test_throws SubTropica.LRDispatchError SubTropica.find_lr_orders_scan(
            groups, [:x, :y, :z], String[], [[[1, 1]]]; transport=t)  # count
        @test_throws SubTropica.LRDispatchError SubTropica.find_lr_orders_scan(
            groups, [:x, :y, :z], String[],
            [[[-2, -1], [0, 0]], [[1, 1]]]; transport=t)              # length
        @test_throws SubTropica.LRDispatchError SubTropica.find_lr_orders_scan(
            groups, [:x, :y, :z], String[], [[[1]], [[1, 1]]]; transport=t)

        # projective:false = INPUT error (typed, names the gate), NOT NOLR
        errp = try
            SubTropica.find_lr_orders_scan(groups, [:x, :y, :z], String[],
                exps; transport=a3_mock((op, f) -> R_SCAN_NONPROJ))
            nothing
        catch e; e end
        @test errp isa SubTropica.LRDispatchError
        @test occursin("projectivity", sprint(showerror, errp))
        @test occursin("NOT a NOLR", sprint(showerror, errp))

        # projective:true + empty orders = the scan's NOLR analog: LOUD typed
        erre = try
            SubTropica.find_lr_orders_scan(groups, [:x, :y, :z], String[],
                exps; transport=a3_mock((op, f) -> R_SCAN_EMPTY))
            nothing
        catch e; e end
        @test erre isa SubTropica.NoLROrderError
        @test occursin("keep_rule=Strict", erre.detail)

        # in-band error / malformed envelope
        @test_throws SubTropica.LRDispatchError SubTropica.find_lr_orders_scan(
            groups, [:x, :y, :z], String[], exps;
            transport=a3_mock((op, f) -> A3_R_ERROR))
        @test_throws SubTropica.LRDispatchError SubTropica.find_lr_orders_scan(
            groups, [:x, :y, :z], String[], exps;
            transport=a3_mock((op, f) -> """{"op":"find_lr_orders_scan","projective":true}"""))
    end

    @testset "factor_table: permutation guard + table passthrough" begin
        R_FT_OK = """{"op":"factor_table","schema_version":2,"hf_version":"1.2.8","order":["x","y"],"polys":["1+x","1+x+y","y"],"stages":[{"var":"x","admissible":[[0,1]],"pool":[0,1],"n_pairs":1,"n_singletons":2,"n_inadmissible":0,"t_build_s":0.001}],"pairs":[{"var":"x","f":0,"g":1,"c":"1","factors":[[2,1]],"oop":false}],"singletons":[{"var":"x","id":0,"deg":1,"coeffs":[{"power":1,"c":"1","factors":[],"oop":false}]}],"stats":{"pairs_total":1,"singletons_total":2,"oop":0,"pair_fallbacks":0,"trial_s":0.0,"fallback_s":0.0}}"""
        groups = [["1+x", "1+x+y"]]

        cap = Any[]
        t = a3_mock((op, f) -> R_FT_OK; capture=cap)
        tab = SubTropica.factor_table(groups, [:x, :y], String[], [:x, :y];
                                      max_pairs=100, transport=t)
        (op, fields) = cap[1]
        @test op == "factor_table"
        @test haskey(fields, "groups") && !haskey(fields, "polys")  # multi-group ALWAYS
        @test fields["order"] == ["x", "y"]
        @test fields["max_pairs"] == 100
        @test !haskey(fields, "max_singletons") && !haskey(fields, "max_response_mb")
        @test !haskey(fields, "algebraic_letters")                  # baseline envelope
        @test tab["order"] == [:x, :y]
        @test tab["polys"] == ["1+x", "1+x+y", "y"]
        @test tab["stats"]["pairs_total"] == 1
        @test length(tab["stages"]) == 1 && tab["stages"][1]["var"] == "x"

        # order must be a permutation of xvars (client-side, engine-mirroring)
        @test_throws SubTropica.LRDispatchError SubTropica.factor_table(
            groups, [:x, :y], String[], [:x]; transport=t)
        @test_throws SubTropica.LRDispatchError SubTropica.factor_table(
            groups, [:x, :y], String[], [:x, :z]; transport=t)
        @test_throws SubTropica.LRDispatchError SubTropica.factor_table(
            groups, [:x, :y], String[], Symbol[]; transport=t)

        # loud guard errors arrive in-band — typed, never truncation
        @test_throws SubTropica.LRDispatchError SubTropica.factor_table(
            groups, [:x, :y], String[], [:x, :y];
            transport=a3_mock((op, f) -> """{"op":"factor_table","error":"max_pairs exceeded"}"""))
        @test_throws SubTropica.LRDispatchError SubTropica.factor_table(
            groups, [:x, :y], String[], [:x, :y];
            transport=a3_mock((op, f) -> """{"op":"factor_table","polys":[]}"""))
    end

    @testset "derive_gauge_from_homogeneous_lr (stDeriveGaugeFromHomogeneousLR port)" begin
        R_SEARCH_XYZ = """{"op":"find_lr_orders","schema_version":2,"hf_version":"1.2.8","best_order":["x","y","z"],"score":3.0,"nolr":false,"strategy":"LR_NoOpt","timing_compute_s":0.0004,"nXVars":3,"nGroups":1,"nPolys":[2]}"""
        letters = ["1+x+y+z", "x*y+z^2"]

        cap = Any[]
        t = a3_mock((op, f) -> R_SEARCH_XYZ; capture=cap)
        (gauge, order, meta) = SubTropica.derive_gauge_from_homogeneous_lr(
            letters, [:x, :y, :z]; transport=t)
        # LAST var of the full order = gauge; rest = integration order (wl:17108)
        @test gauge == :z && order == [:x, :y]
        @test meta["gauge_strategy"] == "derive_from_homogeneous_lr"
        @test meta["gauge"] == "z"
        # dispatched as ONE group over ALL xvars (upstream {letters} call)
        @test cap[1][2]["groups"] == [letters]
        @test cap[1][2]["xvars"] == ["x", "y", "z"]

        # NOLR: upstream \$Failed becomes the LOUD typed refusal
        @test_throws SubTropica.NoLROrderError SubTropica.derive_gauge_from_homogeneous_lr(
            letters, [:x, :y, :z]; transport=a3_mock((op, f) -> A3_R_NOLR))
        # caller bugs: no letters / too few xvars
        @test_throws SubTropica.LRDispatchError SubTropica.derive_gauge_from_homogeneous_lr(
            String[], [:x, :y]; transport=t)
        @test_throws SubTropica.LRDispatchError SubTropica.derive_gauge_from_homogeneous_lr(
            letters, [:x]; transport=t)
        # partial best_order = malformed response (gauge split needs it all)
        @test_throws SubTropica.LRDispatchError SubTropica.derive_gauge_from_homogeneous_lr(
            letters, [:x, :y, :z]; transport=a3_mock((op, f) -> A3_R_SEARCH_X))
    end
end

@testset "unit-3: assemble" begin

    @testset "local exact HlogExpr helpers" begin
        a = a3_hl(1 // 2, [SubTropica.HConst(:Pi) => 2])
        b = a3_hl(1 // 3, [SubTropica.HConst(:Pi) => 2])
        c = SubTropica._a3_hadd(a, b)
        pikey = first(keys(a3_coefdict(a)))
        @test a3_coefdict(c) == Dict(pikey => 5 // big(6))
        # different atom parts stay separate
        d = SubTropica._a3_hadd(a, a3_hl(1 // 1, [SubTropica.HZeta([3]) => 1]))
        @test length(d.terms) == 2
        # exact product: coef and atom powers
        e = SubTropica._a3_hmul(a, a3_hl(4 // 1, [SubTropica.HConst(:Pi) => 1]))
        cd = a3_coefdict(e)
        @test length(cd) == 1 && first(values(cd)) == 2 // big(1)
        @test first(e.terms).atoms == [SubTropica.HConst(:Pi) => 3]
        # cancellation drops the term
        z = SubTropica._a3_hadd(a, SubTropica._a3_hscale(a, -1 // 1))
        @test isempty(z.terms)
    end

    @testset "union_letters + remap" begin
        L1 = a3_letter(1, "1+w+w^2", :w)
        L2 = a3_letter(2, "2+w+w^2", :w; disc="7")
        L1b = a3_letter(1, "2+w+w^2", :w; disc="7")   # same entry, other resp
        (letters, remaps) = SubTropica.union_letters([[L1, L2], [L1b]])
        @test length(letters) == 2
        @test remaps[1] == Dict(1 => 1, 2 => 2)
        @test remaps[2] == Dict(1 => 2)
        @test letters[2].poly == "2+w+w^2" && letters[2].idx == 2
        # inconsistent duplicate is an upstream bug => typed throw
        badL = a3_letter(1, "1+w+w^2", :w; disc="999")
        @test_throws SubTropica.LetterTableMismatch SubTropica.union_letters(
            [[L1], [badL]])
        # remap rewrites HAlgLetter atoms
        e = a3_hl(1 // 1, [SubTropica.HAlgLetter(:Wm, 1) => 1])
        e2 = SubTropica.remap_letters(e, Dict(1 => 7))
        @test first(e2.terms).atoms == [SubTropica.HAlgLetter(:Wm, 7) => 1]
        @test_throws SubTropica.LetterTableMismatch SubTropica.remap_letters(
            e, Dict(3 => 1))
        # remap rewrites Wm/Wp tokens INSIDE HPeriod
        # WORDS too (the letters tier mints them, eq 4.11) — single-pass,
        # so a swap remap cannot double-apply; missing idx throws typed
        ew = a3_hl(1 // 1,
                   [SubTropica.HPeriod(["-1", "0", "Wm_1"]) => 1,
                    SubTropica.HPeriod(["Wp_2", "sqrt_disc_1"]) => 1])
        ew2 = SubTropica.remap_letters(ew, Dict(1 => 2, 2 => 1))
        words = sort([a.word for (a, _) in first(ew2.terms).atoms
                      if a isa SubTropica.HPeriod])
        @test words == sort([["-1", "0", "Wm_2"], ["Wp_1", "sqrt_disc_2"]])
        @test_throws SubTropica.LetterTableMismatch SubTropica.remap_letters(
            a3_hl(1 // 1, [SubTropica.HPeriod(["Wm_5"]) => 1]), Dict(1 => 2))
        # rational-only words pass through untouched (never rewritten)
        # (compare (coef, atoms) — HTerm itself has no == overload)
        er = a3_hl(1 // 1, [SubTropica.HPeriod(["0", "-1"]) => 1])
        rp = SubTropica.remap_letters(er, Dict(1 => 7))
        @test [(t.coef, t.atoms) for t in rp.terms] ==
              [(t.coef, t.atoms) for t in er.terms]
    end

    @testset "align_laurent: general multi-face pole padding" begin
        p1 = (-2, [a3_hl(1 // 1), a3_hl(2 // 1)])          # ε^-2, ε^-1
        p2 = (0, [a3_hl(5 // 1)])                          # ε^0
        (gmin, padded) = SubTropica.align_laurent(
            Tuple{Int,Vector{SubTropica.HlogExpr}}[p1, p2])
        @test gmin == -2
        @test length(padded[1]) == 3 && length(padded[2]) == 3
        @test isempty(padded[1][3].terms)          # face1 padded at ε^0
        @test isempty(padded[2][1].terms)          # face2 padded at ε^-2
        @test a3_constval(padded[2][3]) == 5
    end

    # empty letter table shorthand
    a3_nl = SubTropica.AlgLetter[]
    a3_ev() = Dict{String,Any}("check_divergences" => true,
                               "lr_order" => ["x"])

    @testset "check_divergences policy: HARD-ON, assembling without it throws" begin
        blocks = [[(a3_hl(1 // 1), a3_nl)]]
        pref = (0, [a3_hl(1 // 1)])
        # missing flag
        @test_throws SubTropica.DivergenceCheckViolation SubTropica.assemble(
            0, blocks, pref, 0; evidence=Dict{String,Any}("lr_order" => ["x"]))
        # explicit false is just as fatal
        @test_throws SubTropica.DivergenceCheckViolation SubTropica.assemble(
            0, blocks, pref, 0;
            evidence=Dict{String,Any}("check_divergences" => false,
                                      "lr_order" => ["x"]))
        # face_sum mode is Phase B — refuse loudly
        @test_throws SubTropica.DivergenceCheckViolation SubTropica.assemble(
            0, blocks, pref, 0; evidence=a3_ev(), mode=:face_sum)
        # LR order is REQUIRED provenance
        @test_throws SubTropica.AssembleProvenanceError SubTropica.assemble(
            0, blocks, pref, 0;
            evidence=Dict{String,Any}("check_divergences" => true))
    end

    @testset "assemble: sums, convolution, provenance" begin
        # integrand series I(ε) = 1 - ε   (order 0 split into two
        # contributions 1/3 + 2/3 to exercise the per-order sum)
        blocks = [[(a3_hl(1 // 3), a3_nl), (a3_hl(2 // 3), a3_nl)],
                  [(a3_hl(-1 // 1), a3_nl)]]
        # prefactor series P(ε) = 1 + 2ε
        pref = (0, [a3_hl(1 // 1), a3_hl(2 // 1)])
        L = SubTropica.assemble(0, blocks, pref, 1; evidence=a3_ev())
        @test L isa SubTropica.LaurentSeries{SubTropica.HlogExpr}
        @test L.minorder == 0
        @test length(L.coeffs) == 2
        @test a3_constval(L.coeffs[1]) == 1          # 1*1
        @test a3_constval(L.coeffs[2]) == 1          # 1*2 + (-1)*1
        @test L.provenance == :analytic
        # provenance metadata (task spec: pins + LR order + flag + wall)
        for k in ("hf_bin_sha256", "hf_bin_path", "upstream_wl_sha",
                  "upstream_bin_git", "lr_order", "check_divergences",
                  "wall_s_assemble", "assembled_at", "mode",
                  "requested_order")
            @test haskey(L.evidence, k)
        end
        @test L.evidence["upstream_wl_sha"] ==
              "adac2f722be64337aa095b2b2e7266628b03289b"
        @test L.evidence["upstream_bin_git"] == "ead8c6e"
        @test L.evidence["check_divergences"] === true
        @test L.evidence["wall_s_assemble"] >= 0
        if A3_ENGINE
            @test occursin(r"^[0-9a-f]{64}$", L.evidence["hf_bin_sha256"])
        end
        # prefactor pole shifts the output window
        Lp = SubTropica.assemble(0, blocks, (-1, [a3_hl(1 // 1), a3_hl(1 // 1)]),
                              0; evidence=a3_ev())
        @test Lp.minorder == -1
        @test a3_constval(Lp.coeffs[1]) == 1         # ε^-1: 1*1
        @test a3_constval(Lp.coeffs[2]) == 0         # ε^0: 1*1 + (-1)*1
    end

    @testset "truncation guard (R2): missing input orders throw" begin
        blocks = [[(a3_hl(1 // 1), a3_nl)]]          # I known to ε^0 only
        pref0 = (0, [a3_hl(1 // 1)])                 # P known to ε^0 only
        @test_throws SubTropica.LaurentTruncationError SubTropica.assemble(
            0, blocks, pref0, 1; evidence=a3_ev())   # needs I_1: missing
        pref1 = (0, [a3_hl(1 // 1), a3_hl(1 // 1)])
        @test_throws SubTropica.LaurentTruncationError SubTropica.assemble(
            0, [[(a3_hl(1 // 1), a3_nl)], [(a3_hl(1 // 1), a3_nl)]],
            pref0, 1; evidence=a3_ev())              # needs P_1: missing
        # requested below the leading order
        @test_throws SubTropica.LaurentTruncationError SubTropica.assemble(
            0, blocks, pref0, -1; evidence=a3_ev())
        # explicit zeros make it legal
        Lz = SubTropica.assemble(0,
            [[(a3_hl(1 // 1), a3_nl)], Tuple{SubTropica.HlogExpr,Vector{SubTropica.AlgLetter}}[]],
            pref1, 1; evidence=a3_ev())
        @test a3_constval(Lz.coeffs[2]) == 1
    end

    @testset "multi-face fold: pole padding + letter union end-to-end" begin
        w1 = SubTropica.AlgLetter(1, "1+w+w^2", :w, "1", "-1", "1", "-3")
        w1dup = SubTropica.AlgLetter(1, "1+w+w^2", :w, "1", "-1", "1", "-3")
        f1 = SubTropica.FaceResult(-1,
            [[(a3_hl(1 // 1, [SubTropica.HAlgLetter(:Wm, 1) => 1]), [w1])],
             [(a3_hl(3 // 1), SubTropica.AlgLetter[])]])
        f2 = SubTropica.FaceResult(0,
            [[(a3_hl(4 // 1, [SubTropica.HAlgLetter(:Wm, 1) => 1]), [w1dup])]])
        # the ε^-1 pole in f1 means output ε^0 needs the prefactor through
        # ε^1 — the R2 guard demands the explicit zero:
        pref = (0, [a3_hl(1 // 1), SubTropica.HlogExpr()])
        L = SubTropica.assemble_faces([f1, f2], pref, 0; evidence=a3_ev())
        @test L.minorder == -1
        # and WITHOUT the explicit zero the guard throws (general fold)
        @test_throws SubTropica.LaurentTruncationError SubTropica.assemble_faces(
            [f1, f2], (0, [a3_hl(1 // 1)]), 0; evidence=a3_ev())
        @test length(L.letters) == 1        # unioned across faces
        d = a3_coefdict(L.coeffs[2])        # ε^0: 3 + 4*Wm_1
        @test d[""] == 3
        @test length(d) == 2 && any(v == 4 for v in values(d))
    end

    @testset "mzv_reduce (mock; canned strings from live probes)" begin
      if !isfile(A3_CFG.mzv_data_path)
        @info "mzv table not found — mzv_reduce mock tests skipped (pool needs the table)"
      else
        # live-pinned: mzv_2_1 -> -2*mzv_3 ; mzv_2_1^2 -> 4*mzv_3^2 ;
        # Log2*mzv_2_1 -> -2*Log2*mzv_3 (Log2 rides along linearly)
        canned = Dict(
            "mzv_2_1" => "-2*mzv_3",
            "mzv_2_1^2" => "4*mzv_3^2",
            "Log2*mzv_2_1" => "-2*Log2*mzv_3",
        )
        cap = Any[]
        t = a3_mock(; capture=cap) do op, f
            op == "apply_mzv_reductions" || error("unexpected op $op")
            r = canned[f["f"]]
            return """{"op":"apply_mzv_reductions","result":"$(r)","vars":[]}"""
        end
        e = SubTropica.HlogExpr([
            SubTropica.HTerm("3/4", [SubTropica.HZeta([2, 1]) => 1]),
            SubTropica.HTerm("1", [SubTropica.HZeta([2, 1]) => 2]),
            SubTropica.HTerm("1", [SubTropica.HConst(:Log2) => 1,
                                SubTropica.HZeta([2, 1]) => 1]),
            SubTropica.HTerm("5", [SubTropica.HConst(:Pi) => 1]),  # non-MZV: untouched
        ], Symbol[])
        # fresh memo for a deterministic capture count
        empty!(SubTropica._MZV_RED_MEMO)
        r = SubTropica.mzv_reduce(e; transport=t)
        d = a3_coefdict(r)
        k3 = join(("$(SubTropica._a3_atom_key(SubTropica.HZeta([3])))^1",), "*")
        @test d[k3] == -3 // big(2)                        # 3/4 * (-2)
        @test d["$(SubTropica._a3_atom_key(SubTropica.HZeta([3])))^2"] == 4
        # Log2 comes back as HConst(:Log2) (canonical Log2 atom decision
        # matches parse_coef's _norm_atom fold);
        # ordering: HConst (rank 1) before HZeta (2)
        @test d["$(SubTropica._a3_atom_key(SubTropica.HConst(:Log2)))^1*$(SubTropica._a3_atom_key(SubTropica.HZeta([3])))^1"] == -2
        @test d["$(SubTropica._a3_atom_key(SubTropica.HConst(:Pi)))^1"] == 5
        # vars pool rule: every request carries the table symbol pool
        for (op, f) in cap
            @test "mzv_3" in f["vars"]      # reduction TARGET present
            @test "mzv_3_5" in f["vars"]    # basis token present
        end
        # in-band error surfaces as typed MzvReduceError
        terr = a3_mock((op, f) -> """{"op":"apply_mzv_reductions","error":"Poly: parse error: (-2*mzv_3)"}""")
        empty!(SubTropica._MZV_RED_MEMO)
        @test_throws SubTropica.MzvReduceError SubTropica.mzv_reduce(
            a3_hl(1 // 1, [SubTropica.HZeta([2, 1]) => 1]); transport=terr)
        empty!(SubTropica._MZV_RED_MEMO)
      end
    end

    @testset "_parse_mzv_poly: pinned reduction-output grammar" begin
        P = SubTropica._parse_mzv_poly
        @test a3_coefdict(P("-2*mzv_3")) ==
              Dict("$(SubTropica._a3_atom_key(SubTropica.HZeta([3])))^1" => -2 // big(1))
        d = a3_coefdict(P("1/2*Log2*mzv_2 - 1/4*mzv_3"))
        @test d["$(SubTropica._a3_atom_key(SubTropica.HConst(:Log2)))^1*$(SubTropica._a3_atom_key(SubTropica.HZeta([2])))^1"] == 1 // big(2)
        @test d["$(SubTropica._a3_atom_key(SubTropica.HZeta([3])))^1"] == -1 // big(4)
        @test a3_coefdict(P("mzv_2^2")) ==
              Dict("$(SubTropica._a3_atom_key(SubTropica.HZeta([2])))^2" => 1 // big(1))
        @test a3_coefdict(P("-(Log2)")) ==
              Dict("$(SubTropica._a3_atom_key(SubTropica.HConst(:Log2)))^1" => -1 // big(1))
        @test a3_constval(P("7")) == 7
        @test a3_constval(P("-3/2")) == -3 // 2
        @test isempty(P("0").terms)
        # alternating index tokens round-trip through the m-prefix
        e = P("mzv_1_m3")
        @test first(e.terms).atoms == [SubTropica.HZeta([1, -3]) => 1]
        @test SubTropica._zeta_token([1, -3]) == "mzv_1_m3"
        # outside the grammar => LOUD
        @test_throws SubTropica.MzvReduceError P("mystery_token")
        @test_throws SubTropica.MzvReduceError P("(1+x)/(1+y)")
    end

    @testset "to_ginsh: pinned exports + refusals" begin
        # exports
        e = SubTropica.HlogExpr([
            SubTropica.HTerm("3/2", [SubTropica.HConst(:Pi) => 2]),
            SubTropica.HTerm("-1", [SubTropica.HZeta([3]) => 1,
                                 SubTropica.HConst(:EulerGamma) => 1]),
            SubTropica.HTerm("1", [SubTropica.HLogA("1+s") => 2]),
            SubTropica.HTerm("1", [SubTropica.HHlog("s", ["0", "-1"]) => 1]),
            SubTropica.HTerm("1", [SubTropica.HConst(:Log2) => 1]),
        ], Symbol[:s])
        g = SubTropica.to_ginsh(e)
        @test occursin("(3/2)*Pi^2", g)
        @test occursin("Euler", g) && occursin("zeta(3)", g)
        @test occursin("log(1+s)^2", g)
        @test occursin("G({0,-1},s)", g)
        @test occursin("log(2)", g)
        @test SubTropica.to_ginsh(SubTropica.HlogExpr()) == "0"
        # series form
        L = SubTropica.LaurentSeries(-1, [a3_hl(1 // 2), SubTropica.HlogExpr()])
        @test SubTropica.to_ginsh(L) == [(-1, "(1/2)"), (0, "0")]
        @test SubTropica.ginsh_series_string(L) == "((1/2))*eps^(-1)+(0)*eps^(0)"
        # refusals: no fabricated conventions/constants
        @test_throws SubTropica.GinshExportError SubTropica.to_ginsh(
            a3_hl(1 // 1, [SubTropica.HZeta([2, 1]) => 1]))    # depth-2 token
        @test_throws SubTropica.GinshExportError SubTropica.to_ginsh(
            a3_hl(1 // 1, [SubTropica.HZeta([1, -3]) => 1]))   # alternating basis
        @test_throws SubTropica.GinshExportError SubTropica.to_ginsh(
            a3_hl(1 // 1, [SubTropica.HPeriod(["0", "-1"]) => 1]))
        @test_throws SubTropica.GinshExportError SubTropica.to_ginsh(
            a3_hl(1 // 1, [SubTropica.HAlgLetter(:Wm, 1) => 1]))
        @test_throws SubTropica.GinshExportError SubTropica.to_ginsh(
            a3_hl(1 // 1, [SubTropica.HDelta(:x) => 1]))
        # coef guard: unparsed Mma/function tokens refuse
        @test_throws SubTropica.GinshExportError SubTropica.to_ginsh(
            a3_hl("Log[2]"))
        @test_throws SubTropica.GinshExportError SubTropica.to_ginsh(
            a3_hl("mzv_2"))
    end
end

# --------------------------------------------------------------------------
# LIVE smoke assembly: I(ε) = ∫₀^∞ (1+x)^(-2-ε) dx = 1/(1+ε) = 1 - ε + ε² …
# Per-order integrands (the exp(εb·logP) multinomial with b = -1, unit-2's
# job in the real pipeline) prepared by hand:
#   k=0: ∫ (1+x)^-2 = 1 ; k=1: (-1)·∫ Log[1+x](1+x)^-2 = -1 ;
#   k=2: (1/2)·∫ Log[1+x]²(1+x)^-2 = 1.
# Exercises: live find_lr_orders search + verify_order replay + persistence,
# live hyperflint (check_divergences HARD-ON), live divergent control,
# assembly with provenance, live apply_mzv_reductions, ginsh export.
# --------------------------------------------------------------------------
@testset "unit-3: LIVE smoke (trivial 1-d integral)" begin
    if !A3_ENGINE
        @info "HyperFLINT engine not found at $(A3_CFG.bin) — live smoke skipped"
    else
        mkpath(A3_RUNS)
        rundir = mktempdir(A3_RUNS)
        groups = [["1+x"]]

        # LR search + persistence + replay, live
        (order, meta) = SubTropica.ensure_lr_order(groups, [:x], String[];
            run_dir=rundir, transport=a3_live_transport)
        @test order == [:x]
        @test meta["replayed"] == false
        (order2, meta2) = SubTropica.ensure_lr_order(groups, [:x], String[];
            run_dir=rundir, transport=a3_live_transport)
        @test order2 == [:x]
        @test meta2["replayed"] == true          # verify_order-field replay

        # live divergent control: the flag catches ∫ dx/(x(1+x)) loudly
        rdiv = a3_live_transport("hyperflint",
            Dict{String,Any}("f" => "1/(x*(1+x))", "vars_int" => ["x"],
                             "vars" => ["x"], "check_divergences" => true))
        @test get(rdiv, "divergent", false) === true

        # per-ε-order live integrations (check_divergences HARD-ON)
        exprs = ["1/(1+x)^2", "Log[1+x]/(1+x)^2", "Log[1+x]^2/(1+x)^2"]
        scales = [1 // 1, -1 // 1, 1 // 2]       # b^k / k! with b = -1
        expected = [1 // 1, -1 // 1, 1 // 1]     # of 1/(1+ε)
        blocks = Vector{Tuple{SubTropica.HlogExpr,Vector{SubTropica.AlgLetter}}}[]
        hf_wall = 0.0
        for (k, ex) in enumerate(exprs)
            fields = Dict{String,Any}(
                "vars_int" => ["x"], "vars" => ["x"],
                "check_divergences" => true, "canonical_emission" => true)
            k == 1 ? (fields["f"] = ex) : (fields["expr"] = ex)
            resp = a3_live_transport("hyperflint", fields)
            @test haskey(resp, "result")
            @test !haskey(resp, "error") && get(resp, "failed", false) == false
            @test get(resp, "divergent", false) == false
            hf_wall += get(resp, "timing_compute_s", 0.0)
            acc = SubTropica.HlogExpr()
            for term in resp["result"]
                @test isempty(term["key"])       # trivial case: no residual words
                acc = SubTropica._a3_hadd(acc,
                    SubTropica._parse_mzv_poly(String(term["coef"])))
            end
            push!(blocks, [(SubTropica._a3_hscale(acc, scales[k]),
                            SubTropica.AlgLetter[])])
        end

        ev = Dict{String,Any}(
            "check_divergences" => true,
            "lr_order" => String.(order),
            "wall_s_hf" => hf_wall,
            "run_dir" => rundir)
        # prefactor ≡ 1: the R2 guard treats every series as truncated, so
        # the vanishing orders are EXPLICIT zeros
        pref = (0, [a3_hl(1 // 1), SubTropica.HlogExpr(), SubTropica.HlogExpr()])
        L = SubTropica.assemble(0, blocks, pref, 2; evidence=ev)
        @test L.minorder == 0
        @test [a3_constval(c) for c in L.coeffs] == expected
        @test occursin(r"^[0-9a-f]{64}$", L.evidence["hf_bin_sha256"])
        @test L.evidence["upstream_bin_git"] == "ead8c6e"
        @test L.evidence["lr_order"] == ["x"]
        @test L.evidence["wall_s_hf"] > 0
        @test SubTropica.to_ginsh(L) == [(0, "(1)"), (1, "(-1)"), (2, "(1)")]

        # live MZV reduction pass: ζ-token convention pinned by the table
        empty!(SubTropica._MZV_RED_MEMO)
        red = SubTropica.mzv_reduce(
            a3_hl(1 // 1, [SubTropica.HZeta([2, 1]) => 1]);
            cfg=A3_CFG, transport=a3_live_transport)
        @test a3_coefdict(red) ==
              Dict("$(SubTropica._a3_atom_key(SubTropica.HZeta([3])))^1" => -2 // big(1))

        rm(rundir; recursive=true, force=true)
    end
end
