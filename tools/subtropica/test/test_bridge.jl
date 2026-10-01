# test_bridge.jl — OWNER: unit-1. Bridge transport + version gate + error
# taxonomy + serializer round-trip law + hyperflint/period wrappers.
# Self-contained: run either via test/runtests.jl (auto-included) or
#   ulimit -v 32505856; julia --project=tools/subtropica tools/subtropica/test/test_bridge.jl
# Skips engine-dependent testsets cleanly if the READ-ONLY engine is absent.

using Test
using SubTropica
using Nemo

@testset "unit-1 bridge" begin

    cfg = SubTropica.HFConfig()
    have_hf = isfile(cfg.bin) && isfile(cfg.mzv_data_path)
    have_hf || @warn "HyperFLINT engine not found at $(cfg.bin) — engine testsets skipped"

    # ---- caller-side misuse: ErrorException BEFORE any spawn ---------------
    @testset "required-field validation" begin
        @test_throws ErrorException SubTropica.hf_call("parse_expr", Dict{String,Any}())
        @test_throws ErrorException SubTropica.hf_call("hyperflint", Dict{String,Any}("vars" => ["x"]))
        @test_throws ErrorException SubTropica.hf_call("zero_inf_period", Dict{String,Any}("vars" => ["x"]))
    end

    # ---- version gate (stamp + canary, RT-vi) ------------------------------
    @testset "version gate" begin
        # a binary with no HF_VERSION stamp must FAIL the gate loudly
        badcfg = SubTropica.HFConfig(bin="/bin/false")
        @test_throws SubTropica.HFGateError SubTropica.hf_gate(badcfg)
        misscfg = SubTropica.HFConfig(bin="/nonexistent/hyperflint.sh")
        @test_throws SubTropica.HFGateError SubTropica.hf_gate(misscfg)
        if have_hf
            @test SubTropica.hf_gate(cfg) === true
            @test SubTropica._HF_GATE_OK[cfg.bin]   # memoized
        end
    end

    # ---- Log[rational] → prime-basis canonical emission fold
    # (serialize.jl fold_log_atoms + mma_string wiring). HF-free:
    # the fold is exact HlogExpr algebra; the emitted forms are already
    # covered by the round-trip law testset (Log[2]/Log[p] parse shapes).
    @testset "fold_log_atoms canonical emission" begin
        H = SubTropica.HlogExpr
        T = SubTropica.HTerm
        LA = SubTropica.HLogA
        L2 = SubTropica.HConst(:Log2)
        one_atom(a; coef="1", k=1) = H([T(coef, [a => k])], Symbol[])

        # basis pins: Log[1/2] → −Log2; Log[4] → 2·Log2; Log[6] → Log2+Log[3]
        f = SubTropica.fold_log_atoms(one_atom(LA("1/2")))
        @test length(f.terms) == 1 && f.terms[1].coef == "-1" &&
              f.terms[1].atoms == [L2 => 1]
        f = SubTropica.fold_log_atoms(one_atom(LA("4")))
        @test f.terms[1].coef == "2" && f.terms[1].atoms == [L2 => 1]
        f = SubTropica.fold_log_atoms(one_atom(LA("6")))
        @test length(f.terms) == 2
        @test any(t -> t.atoms == [L2 => 1] && t.coef == "1", f.terms)
        @test any(t -> t.atoms == [LA("3") => 1] && t.coef == "1", f.terms)
        # composite rational: Log[9/8] → 2·Log[3] − 3·Log2
        f = SubTropica.fold_log_atoms(one_atom(LA("9/8")))
        @test any(t -> t.atoms == [L2 => 1] && t.coef == "-3", f.terms)
        @test any(t -> t.atoms == [LA("3") => 1] && t.coef == "2", f.terms)
        # power expansion: Log[6]^2 → Log2² + 2·Log2·Log[3] + Log[3]²
        f = SubTropica.fold_log_atoms(one_atom(LA("6"); k=2))
        @test length(f.terms) == 3
        @test any(t -> t.atoms == [L2 => 2] && t.coef == "1", f.terms)
        @test any(t -> t.atoms == [L2 => 1, LA("3") => 1] && t.coef == "2", f.terms)
        @test any(t -> t.atoms == [LA("3") => 2] && t.coef == "1", f.terms)
        # single-prime negative power: Log[4]^(-1) → (1/2)·Log2^(-1)
        f = SubTropica.fold_log_atoms(one_atom(LA("4"); k=-1))
        @test f.terms[1].coef == "1/2" && f.terms[1].atoms == [L2 => -1]
        # Log[1] kills the term exactly; Log[1] in a denominator is LOUD
        f = SubTropica.fold_log_atoms(H([T("7", [LA("1") => 1]), T("5", [])], Symbol[]))
        @test length(f.terms) == 1 && f.terms[1].coef == "5"
        @test_throws SubTropica.SubTropicaSerializeError SubTropica.fold_log_atoms(
            one_atom(LA("1"); k=-1))
        # pass-throughs: kinematic arg, non-positive rational, multi-prime
        # negative power — all KEPT (and the no-fold path returns e itself)
        for a in (LA("1+s"), LA("-2"), LA("0"))
            e0 = one_atom(a)
            @test SubTropica.fold_log_atoms(e0) === e0
        end
        f = SubTropica.fold_log_atoms(one_atom(LA("6"); k=-1))
        @test f.terms[1].atoms == [LA("6") => -1]

        # THE REGRESSION (mandate): eq 4.3 ε⁻³ combination
        # 3·Log2 − Log[1/2] + Log[1/4] − Log[4] canonicalizes to literal 0.
        combo = H([T("3", [L2 => 1]), T("-1", [LA("1/2") => 1]),
                   T("1", [LA("1/4") => 1]), T("-1", [LA("4") => 1])], Symbol[])
        @test isempty(SubTropica.fold_log_atoms(combo).terms)
        @test SubTropica.mma_string(combo) == "0"                      # emission side
        @test SubTropica.canonical_text(SubTropica.fold_log_atoms(combo)) == "0"
        # mutation control: breaking one coefficient must NOT cancel
        combo_bad = H([T("3", [L2 => 1]), T("-1", [LA("1/2") => 1]),
                       T("1", [LA("1/4") => 1]), T("-2", [LA("4") => 1])], Symbol[])
        @test SubTropica.mma_string(combo_bad) != "0"

        # idempotence + kinematic coefficients preserved
        e2 = H([T("(s + 1)/(s - 1)", [LA("6") => 1])], [:s])
        f2 = SubTropica.fold_log_atoms(e2)
        @test SubTropica.fold_log_atoms(f2) === f2
        @test all(occursin("s", t.coef) for t in f2.terms)
        # round-trip parity with parse: folded emission re-parses to the
        # same folded object (Log[2]→HConst(:Log2), Log[3]→HLogA("3"))
        s6 = SubTropica.mma_string(one_atom(LA("6")))
        pe = SubTropica.parse_coef(s6, Symbol[])
        @test SubTropica.canonical_text(pe) ==
              SubTropica.canonical_text(SubTropica.fold_log_atoms(one_atom(LA("6"))))
    end

    if have_hf
        # ---- unary-minus canary via the public path ------------------------
        @testset "canary behavior" begin
            r = SubTropica.hf_call("eval", Dict{String,Any}(
                "a" => "-x^2", "vars" => ["x"], "values" => ["3"]); cfg=cfg)
            @test r["result"] == "-9"      # -(x^2), CONTRACTS.md §a live pin
        end

        # ---- in-band error taxonomy (typed, never defaults; rc never gated) -
        @testset "error taxonomy" begin
            @test_throws SubTropica.HFError SubTropica.hf_call("nonsense_op",
                Dict{String,Any}("expr" => "x"); cfg=cfg)          # "unknown op", rc=2
            # non-linearly-reducible => {"failed":true}, rc=0 (control outcome, CONTRACTS.md §a)
            @test_throws SubTropica.HFFailed SubTropica.hf_call("hyperflint",
                Dict{String,Any}("expr" => "1/(1+x+x^2)", "vars_int" => ["x"],
                                 "vars" => ["x"], "check_divergences" => true); cfg=cfg)
            # divergent with checks hard-ON => {"divergent":true}, rc=0
            @test_throws SubTropica.HFDivergent SubTropica.hf_call("hyperflint",
                Dict{String,Any}("expr" => "1/(x*(1+x))", "vars_int" => ["x"],
                                 "vars" => ["x"], "check_divergences" => true); cfg=cfg)
            # parse error => in-band {"error":..}
            @test_throws SubTropica.HFError SubTropica.hf_call("parse_expr",
                Dict{String,Any}("expr" => "1+", "vars" => ["x"]); cfg=cfg)
            # flatness: nested object refused pre-flight
            @test_throws SubTropica.HFError SubTropica.hf_call("eval",
                Dict{String,Any}("a" => "x", "vars" => ["x"], "values" => ["1"],
                                 "nested" => Dict("q" => 1)); cfg=cfg)
            # timeout: killed by PID, typed
            @test_throws SubTropica.HFTimeout SubTropica.hf_call("hyperflint",
                Dict{String,Any}("expr" => "Log[1+x]^2/(x*(1+x))",
                                 "vars_int" => ["x"], "vars" => ["x"],
                                 "check_divergences" => true); cfg=cfg,
                timeout_s=0.05)
        end

        # ---- serializer ROUND-TRIP LAW (CONTRACTS.md §d) --------------------
        # parse_expr must ACCEPT every emitted string; semantic equality is
        # checked through the `eval` op at exact rational points vs Nemo
        # (canonicalization may rewrite bytes — "x^-2" -> "1/x^2" live pin).
        @testset "serializer round-trip law" begin
            R, (x, y) = polynomial_ring(Nemo.QQ, ["x", "y"])
            polys = [3 * x^2 * y - y + R(Nemo.QQ(1, 2)),
                     -x^2,                       # the canary shape
                     x - y,
                     R(Nemo.QQ(-4, 6)),          # content-canonical: -2/3
                     2 * x^3 * y^2 - 7 * x + 5]
            pts = [["2/3", "-5"], ["3", "1"], ["-7/9", "2"]]
            parseq = v -> begin
                parts = split(v, "/")
                Nemo.QQ(parse(BigInt, parts[1]),
                        length(parts) == 2 ? parse(BigInt, parts[2]) : big(1))
            end
            for p in polys
                s = SubTropica.mma_string(p)
                @test !occursin("//", s)         # Julia-ism must never leak
                rp = SubTropica.hf_call("parse_expr",
                    Dict{String,Any}("expr" => s, "vars" => ["x", "y"]); cfg=cfg)
                @test haskey(rp, "canonical")    # parse SUCCESS (no error key)
                for vals in pts
                    re = SubTropica.hf_call("eval", Dict{String,Any}(
                        "a" => s, "vars" => ["x", "y"], "values" => vals); cfg=cfg)
                    exact = evaluate(p, [parseq(v) for v in vals])
                    @test String(re["result"]) == SubTropica.mma_string(exact)
                end
            end
            # LogIntegrand emission parses
            li = SubTropica.LogIntegrand(x - y, x * y, [x + 1 => 2, y + 2 => 1],
                                      [:x, :y], Symbol[])
            sli = SubTropica.mma_string(li)
            @test occursin("Log[", sli)
            rli = SubTropica.hf_call("parse_expr",
                Dict{String,Any}("expr" => sli, "vars" => ["x", "y"]); cfg=cfg)
            @test occursin("Hlog", String(rli["canonical"]))  # Log -> Hlog[..,{0}]
            # HlogExpr emission: HF-parseable subset round-trips; internal
            # atoms REFUSE (never sent toward HF — CONTRACTS.md §d)
            eok = SubTropica.HlogExpr([SubTropica.HTerm("1/2", [SubTropica.HConst(:Pi) => 2]),
                                    SubTropica.HTerm("-1", [SubTropica.HLogA("2") => 1])],
                                   Symbol[])
            sok = SubTropica.mma_string(eok)
            rok = SubTropica.hf_call("parse_expr",
                Dict{String,Any}("expr" => sok, "vars" => ["x"]); cfg=cfg)
            @test haskey(rok, "canonical")
            for bad in (SubTropica.HZeta([2]), SubTropica.HDelta(:x4),
                        SubTropica.HPeriod(["0", "-1"]), SubTropica.HConst(:EulerGamma),
                        SubTropica.HAlgLetter(:Wm, 1))
                ebad = SubTropica.HlogExpr([SubTropica.HTerm("1", [bad => 1])], Symbol[])
                @test_throws SubTropica.SubTropicaSerializeError SubTropica.mma_string(ebad)
            end
        end

        # ---- hyperflint wrapper + period resolution (T1-T3 live pins) -------
        @testset "hf_integrate + resolve_periods" begin
            # T2: ∫ Log[1+x]/(x(1+x)) = ζ(2)
            terms, letters = SubTropica.hf_integrate("Log[1+x]/(x*(1+x))",
                                                  [:x], [:x]; cfg=cfg)
            @test isempty(letters)
            @test length(terms) == 1
            @test terms[1][2] == [["0", "-1"]]
            @test SubTropica.canonical_text(SubTropica.resolve_periods(terms; cfg=cfg)) == "mzv_2"
            # T1: ∫ 1/((1+x)(2+x)) = log 2 — canonical Log2 atom decision
            # : canonical text is "Log2", not "Log[2]"
            t1, _ = SubTropica.hf_integrate("1/((1+x)*(2+x))", [:x], [:x]; cfg=cfg)
            @test SubTropica.canonical_text(SubTropica.resolve_periods(t1; cfg=cfg)) == "Log2"
            # T3: ∫ Log[1+x]^2/(x(1+x)) = 2ζ(3)
            t3, _ = SubTropica.hf_integrate("Log[1+x]^2/(x*(1+x))", [:x], [:x]; cfg=cfg)
            @test SubTropica.canonical_text(SubTropica.resolve_periods(t3; cfg=cfg)) == "2*mzv_3"
            # kinematic single-letter word: ∫ 1/((1+x)(a+x)) = log(a)/(a-1)
            ta, _ = SubTropica.hf_integrate("1/((1+x)*(a+x))", [:x], [:x, :a]; cfg=cfg)
            @test ta[1][2] == [["-a"]]
            ra = SubTropica.resolve_periods(ta; cfg=cfg)
            @test SubTropica.canonical_text(ra) == "(1)/(a - 1)*Log[a]"
            # multi-letter kinematic word REFUSES (strict) / HPeriod (strict=false)
            tb, _ = SubTropica.hf_integrate("Log[1+x]/((1+x)*(a+x))", [:x], [:x, :a]; cfg=cfg)
            @test tb[1][2] == [["-a", "-1"]]
            @test_throws SubTropica.SubTropicaEvalError SubTropica.resolve_periods(tb; cfg=cfg)
            rb = SubTropica.resolve_periods(tb; cfg=cfg, strict=false)
            @test occursin("ZeroInfPeriod[{-a,-1}]", SubTropica.canonical_text(rb))
            # ---- Wm/Wp coef path (algebraic-letter gap closure):
            # 1/(x^2+3x+1) at the letters tier mints the quadratic letter
            # (disc 5) and returns coefs of the refused-before shape
            # (2/(Wm_1 - Wp_1)); hf_integrate now normalizes them through
            # HF's OWN RT-v chain — every parsed coef is atom-canonical
            # (sqrt_disc power ≤1, no multi-term atom denominators).
            tw, lw = SubTropica.hf_integrate("1/(x^2+3*x+1)", [:x], [:x]; cfg=cfg,
                                          algebraic_letters=true)
            @test length(lw) == 1 && lw[1].poly == "x^2 + 3*x + 1"
            @test lw[1].disc == "5"
            @test !isempty(tw)
            sawdisc = false
            for (ce, _) in tw
                for t in ce.terms, (a, k) in t.atoms
                    if a isa SubTropica.HAlgLetter
                        @test a.kind === :sqrt_disc && k == 1
                        sawdisc = true
                    end
                end
            end
            @test sawdisc          # the normalized coefs carry sqrt_disc_1
            # letter-bearing boundary words: typed refusal (strict) /
            # explicit HPeriod atom (strict=false) — never a crash in the
            # rational closed-form branches (live bug class)
            @test_throws SubTropica.SubTropicaEvalError SubTropica.zero_inf_period_value(
                ["Wm_1"]; cfg=cfg, letters=lw)
            pv = SubTropica.zero_inf_period_value(["Wm_1"]; cfg=cfg, letters=lw,
                                               strict=false)
            @test pv.terms[1].atoms == [SubTropica.HPeriod(["Wm_1"]) => 1]
            @test_throws SubTropica.SubTropicaEvalError SubTropica.zero_inf_period_value(
                ["0", "0", "Wm_1"]; cfg=cfg, letters=lw)
            rw = SubTropica.resolve_periods(tw; cfg=cfg, letters=lw, strict=false)
            @test occursin("ZeroInfPeriod", SubTropica.canonical_text(rw))
            # vars_int ⊄ vars refused
            @test_throws SubTropica.HFError SubTropica.hf_integrate("1/(1+x)", [:x], [:y]; cfg=cfg)
            # on_period verification hook fires per resolved word
            seen = Vector{Vector{String}}()
            SubTropica.resolve_periods(terms; cfg=cfg,
                                    on_period=(w, v) -> push!(seen, w))
            @test seen == [["0", "-1"]]
        end
    end
end
