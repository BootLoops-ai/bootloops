# test_hlogexpr.jl — OWNER: unit-1. Token parser, canonical text form,
# exact HlogExpr arithmetic, HF-op Wm/Wp normalization (RT-v), and THE
# PERIOD UNIT GATE: every ZeroInfPeriod residual key in the corpus is sent
# to HF's zero_inf_period op AND cross-evaluated numerically via the
# independent ginac_gpl CLI to ≥60 dps (DESIGN.md RT-extra-period /
# gate A(v)), with synthetic-truth positive controls and a non-rational
# ×(1+δπ) negative control (known pitfall: rational perturbations alias).
# Self-contained: run via runtests.jl or
#   ulimit -v 32505856; julia --project=tools/subtropica tools/subtropica/test/test_hlogexpr.jl

using Test
using SubTropica

@testset "unit-1 hlogexpr" begin

    cfg = SubTropica.HFConfig()
    have_hf = isfile(cfg.bin) && isfile(cfg.mzv_data_path)
    have_ginac = isfile(SubTropica._ginac_bin())
    have_hf || @warn "HyperFLINT engine not found — HF-dependent testsets skipped"
    have_ginac || @warn "ginac_gpl not found — numeric-gate testsets skipped"

    lt = SubTropica.AlgLetter(1, "u^2+3*u+1", :u, "1", "-3", "1", "5")

    # ---- parser: tokens, precedence, canonicalization ----------------------
    @testset "parse_coef basics" begin
        # content-canonical rationals (gcd, sign-on-numerator)
        @test SubTropica.canonical_text(SubTropica.parse_coef("4/6", Symbol[])) == "2/3"
        @test SubTropica.canonical_text(SubTropica.parse_coef("(-1)", Symbol[])) == "-1"
        @test SubTropica.canonical_text(SubTropica.parse_coef("1/2/3", Symbol[])) == "1/6"
        @test SubTropica.canonical_text(SubTropica.parse_coef("2^-2", Symbol[])) == "1/4"
        # unary minus vs ^ (the canary precedence, live-pinned)
        e1 = SubTropica.parse_coef("-x^2", [:x])
        e2 = SubTropica.parse_coef("-(x^2)", [:x])
        e3 = SubTropica.parse_coef("(-x)^2", [:x])
        @test SubTropica.canonical_text(e1) == SubTropica.canonical_text(e2)
        @test SubTropica.canonical_text(e3) == "x^2"
        # rational functions of kinematics stay in coef strings (Nemo-canonical)
        @test SubTropica.canonical_text(SubTropica.parse_coef("(-1/(a - 1))", [:a])) ==
              "(-1)/(a - 1)"
        # x^-2 == 1/x^2 (negative integer exponents, live parse_expr pin)
        @test SubTropica.canonical_text(SubTropica.parse_coef("x^-2", [:x])) ==
              SubTropica.canonical_text(SubTropica.parse_coef("1/x^2", [:x]))
        # token dictionary (SubTropica.wl:12261-12321)
        @test SubTropica.parse_coef("Log2", Symbol[]).terms[1].atoms ==
              [SubTropica.HConst(:Log2) => 1]
        # canonical Log2 atom convention (convention):
        # HLogA("2") folds to HConst(:Log2) everywhere — regression tests.
        @test SubTropica.parse_coef("Log[2]", Symbol[]).terms[1].atoms ==
              [SubTropica.HConst(:Log2) => 1]
        let e = SubTropica.parse_coef("Log2 + Log[2]", Symbol[])
            @test length(e.terms) == 1                       # the split is DEAD
            @test e.terms[1].atoms == [SubTropica.HConst(:Log2) => 1]
            @test e.terms[1].coef == "2"
            @test SubTropica.canonical_text(e) == "2*Log2"
        end
        # arithmetic chokepoint: a manually-built HLogA("2") term merges
        # with an HConst(:Log2) term through hlog_add (_canon_atomvec fold)
        let a = SubTropica.HlogExpr([SubTropica.HTerm("1", [SubTropica.HLogA("2") => 1])], Symbol[]),
            b = SubTropica.HlogExpr([SubTropica.HTerm("3", [SubTropica.HConst(:Log2) => 1])], Symbol[])
            c = SubTropica.hlog_add(a, b)
            @test length(c.terms) == 1
            @test c.terms[1].atoms == [SubTropica.HConst(:Log2) => 1]
            @test c.terms[1].coef == "4"
        end
        # numeric sanity: the canonical atom still evaluates to log(2)
        @test SubTropica._matched_digits(
                  SubTropica.hlog_const_bigfloat(SubTropica.parse_coef("Log[2]", Symbol[]); digits=60),
                  log(BigFloat(2))) >= 60
        @test SubTropica.parse_coef("mzv_1_m3", Symbol[]).terms[1].atoms ==
              [SubTropica.HZeta([1, -3]) => 1]
        @test SubTropica.parse_coef("Log[1]", Symbol[]).terms == SubTropica.HTerm[]  # Log[1]=0
        @test SubTropica.canonical_text(SubTropica.parse_coef("I*Pi*delta[x4]", Symbol[])) ==
              "I*Pi*delta[x4]"
        @test SubTropica.parse_coef("Hlog[s,{0,-1,t}]", [:s, :t]).terms[1].atoms ==
              [SubTropica.HHlog("s", ["0", "-1", "t"]) => 1]
        # algebraic letters need a table entry (wiring-bug refusal, wl:12289-12297)
        @test_throws SubTropica.SubTropicaParseError SubTropica.parse_coef("Wm_1", Symbol[])
        @test SubTropica.parse_coef("Wm_1", Symbol[]; letters=[lt]).terms[1].atoms ==
              [SubTropica.HAlgLetter(:Wm, 1) => 1]
        # WmOverWp_i = Wm_i * Wp_i^-1 (wl:12274-12285); single-term atom
        # denominators are representable
        ew = SubTropica.parse_coef("WmOverWp_1", Symbol[]; letters=[lt])
        @test ew.terms[1].atoms == [SubTropica.HAlgLetter(:Wm, 1) => 1,
                                    SubTropica.HAlgLetter(:Wp, 1) => -1]
        @test SubTropica.canonical_text(SubTropica.parse_coef("1/Wp_1", Symbol[]; letters=[lt])) ==
              SubTropica.canonical_text(SubTropica.parse_coef("Wp_1^-1", Symbol[]; letters=[lt]))
        # back_substitute sqrt-disc ratio shapes are now REPRESENTABLE
        # (Wm/Wp gap closure): sqrt_disc_1^2 folds to the table
        # disc and a + b*sqrt_disc denominators rationalize by conjugate —
        # both exact, table-driven (CONTRACTS.md §d). For u^2+3*u+1
        # (disc 5): (-d/2 - 3/2)/(d/2 - 3/2) = 7/2 + 3d/2 exactly.
        @test SubTropica.canonical_text(SubTropica.parse_coef(
            "(-1/2*sqrt_disc_1 - 3/2)/(1/2*sqrt_disc_1 - 3/2)", Symbol[]; letters=[lt])) ==
            "7/2 + 3/2*sqrt_disc_1"
        # the sqrt_disc^2 -> disc exit fold (ONE canonical spelling)
        @test SubTropica.canonical_text(SubTropica.parse_coef(
            "sqrt_disc_1^2", Symbol[]; letters=[lt])) == "5"
        @test SubTropica.canonical_text(SubTropica.parse_coef(
            "1/sqrt_disc_1", Symbol[]; letters=[lt])) == "1/5*sqrt_disc_1"
        # a disc-less table (stubs) keeps the old typed refusal
        ltnod = SubTropica.AlgLetter(1, "u^2+3*u+1", :u, "1", "-3", "1", "")
        @test_throws SubTropica.SubTropicaParseError SubTropica.parse_coef(
            "(-1/2*sqrt_disc_1 - 3/2)/(1/2*sqrt_disc_1 - 3/2)", Symbol[]; letters=[ltnod])
        # multi-term Wm/Wp denominators still refuse typed WITHOUT the
        # normalize_wmwp router (plain parse_coef is engine-free)
        @test_throws SubTropica.SubTropicaParseError SubTropica.parse_coef(
            "2/(Wm_1 - Wp_1)", Symbol[]; letters=[lt])
        # unknown symbols are contract drift — loud
        @test_throws SubTropica.SubTropicaParseError SubTropica.parse_coef("qq_bogus", Symbol[])
        @test_throws SubTropica.SubTropicaParseError SubTropica.parse_coef("Sqrt[5]", Symbol[])
        @test_throws SubTropica.SubTropicaParseError SubTropica.parse_coef("1+", Symbol[])
    end

    # ---- canonical text: determinism + idempotence --------------------------
    @testset "canonical_text determinism" begin
        a = SubTropica.parse_coef("mzv_2 + 2*Log2^2 - Pi*I", Symbol[])
        b = SubTropica.parse_coef("-I*Pi + Log2^2 + mzv_2 + Log2*Log2", Symbol[])
        @test SubTropica.canonical_text(a) == SubTropica.canonical_text(b)
        for (s, vars) in (("1/2*mzv_2", Symbol[]), ("(-1/(a - 1))", [:a]),
                          ("mzv_2 + 1/2*Log2^2", Symbol[]),
                          ("-Log2*mzv_2 + 7/8*mzv_3", Symbol[]),
                          ("I*Pi*delta[x4] + s/(s + t)", [:s, :t]))
            ct = SubTropica.canonical_text(SubTropica.parse_coef(s, vars))
            @test SubTropica.canonical_text(SubTropica.parse_coef(ct, vars)) == ct  # idempotent
        end
    end

    # ---- exact arithmetic (one owner for unit-2/unit-3 semantics) -----------
    @testset "hlog arithmetic" begin
        a = SubTropica.parse_coef("1/2*mzv_2 + Log2", Symbol[])
        b = SubTropica.parse_coef("2*mzv_2", Symbol[])
        @test SubTropica.canonical_text(SubTropica.hlog_mul(a, b)) ==
              SubTropica.canonical_text(SubTropica.parse_coef("mzv_2^2 + 2*Log2*mzv_2", Symbol[]))
        @test SubTropica.canonical_text(SubTropica.hlog_add(a, a)) ==
              SubTropica.canonical_text(SubTropica.parse_coef("mzv_2 + 2*Log2", Symbol[]))
        # cancellation to exact zero
        @test SubTropica.canonical_text(SubTropica.hlog_add(a, SubTropica.hlog_scale(a, -1))) == "0"
        # rational and string scaling (string parsed exactly over a.vars)
        @test SubTropica.canonical_text(SubTropica.hlog_scale(b, 1//4)) ==
              SubTropica.canonical_text(SubTropica.parse_coef("1/2*mzv_2", Symbol[]))
        k = SubTropica.parse_coef("s*mzv_3", [:s, :t])
        @test SubTropica.canonical_text(SubTropica.hlog_scale(k, "t/(s + t)")) ==
              SubTropica.canonical_text(SubTropica.parse_coef("s*t/(s + t)*mzv_3", [:s, :t]))
        # var-context union across operands
        u = SubTropica.hlog_mul(SubTropica.parse_coef("s", [:s]), SubTropica.parse_coef("t", [:t]))
        @test u.vars == [:s, :t]
        @test SubTropica.canonical_text(u) == "s*t"
    end

    # ---- Wm/Wp normalization through HF's OWN ops (RT-v) --------------------
    if have_hf
        @testset "normalize_wmwp (HF ops, live-pinned chain)" begin
            # Vieta product: Wm·Wp = product = 1 (simplify_with_vieta)
            e, s = SubTropica.normalize_wmwp("Wm_1*Wp_1", [lt]; vars=[:u], cfg=cfg)
            @test s == "1" && SubTropica.canonical_text(e) == "1"
            # Vieta sum: Wm+Wp = -3 (back_substitute leg)
            e, s = SubTropica.normalize_wmwp("Wm_1 + Wp_1", [lt]; vars=[:u], cfg=cfg)
            @test s == "-3" && SubTropica.canonical_text(e) == "-3"
            # difference -> -sqrt_disc_1 (back_substitute)
            e, s = SubTropica.normalize_wmwp("Wm_1 - Wp_1", [lt]; vars=[:u], cfg=cfg)
            @test s == "-sqrt_disc_1"
            @test e.terms[1].atoms == [SubTropica.HAlgLetter(:sqrt_disc, 1) => 1]
            # ratio -> WmOverWp_1: combine_wm_wp_ratios MUST run before
            # back_substitute (which smears ratios into sqrt-denominator
            # shapes neither op recovers — live pin)
            e, s = SubTropica.normalize_wmwp("Wm_1/Wp_1", [lt]; vars=[:u], cfg=cfg)
            @test s == "WmOverWp_1"
            @test e.terms[1].atoms == [SubTropica.HAlgLetter(:Wm, 1) => 1,
                                       SubTropica.HAlgLetter(:Wp, 1) => -1]
            # positional-replay validation: idx gap refused
            lt7 = SubTropica.AlgLetter(7, "u^2+3*u+1", :u, "1", "-3", "1", "5")
            @test_throws SubTropica.SubTropicaEvalError SubTropica.normalize_wmwp(
                "Wm_7", [lt7]; vars=[:u], cfg=cfg)
        end

        # ---- normalize_wmwp v2: structural routing (gap closure)
        # The engine cannot read its own parenthesized "(N/D)" coef emission
        # (in-band Poly parse error) and SILENTLY mis-parses "N/(D)*t" (t
        # lands in the denominator) — live probes, CONTRACTS.md §d. These
        # shapes now go through the parser-armed router: each multi-term
        # Wm/Wp denominator is routed as the single-division body "1/(D)".
        @testset "normalize_wmwp structural routing (eq 4.11 shapes)" begin
            # the eq 4.11 letter: x2^2+1, Gaussian roots, disc -4
            lg = SubTropica.AlgLetter(1, "x2^2 + 1", :x2, "1", "0", "1", "-4")
            # k0 refused coef, EXACTLY as returned by the letters tier —
            # including vars=Symbol[] (the reroute config whose missing
            # letter var was the HFTransportError root cause)
            e, s = SubTropica.normalize_wmwp("(2/(Wm_1 - Wp_1))", [lg];
                                          vars=Symbol[], cfg=cfg)
            @test SubTropica.canonical_text(e) == "1/2*sqrt_disc_1"
            # R'1 4-term sum coef (mzv/Log2 tokens + Vieta-collapsible ratio
            # -> ratio is exactly 0 for this letter: -p-s+1 = 0)
            r1c = "((-Wm_1*Wp_1 - Wm_1 - Wp_1 + 1)/(Wm_1*Wp_1 + Wm_1 + " *
                  "Wp_1 + 1)) + (-2)*Log2*mzv_2 + ((-Wm_1*Wp_1 - Wm_1 - " *
                  "Wp_1 + 1)/(Wm_1*Wp_1 + Wm_1 + Wp_1 + 1))*mzv_2 + (-1)*mzv_3"
            e, s = SubTropica.normalize_wmwp(r1c, [lg]; vars=Symbol[], cfg=cfg)
            @test SubTropica.canonical_text(e) == "-2*Log2*mzv_2 - mzv_3"
            # R'1 non-symmetric cubic ratio: hand value at Wm=-i, Wp=i is
            # (8-4i)/(-4) = -2+i = -2 + sqrt_disc_1/2
            e, s = SubTropica.normalize_wmwp(
                "(-3*Wm_1^3 + 4*Wm_1^2*Wp_1 - 2*Wm_1^2 + 6*Wm_1*Wp_1 - " *
                "Wm_1 + 2*Wp_1)/(Wm_1^3 - Wm_1^2*Wp_1 + 2*Wm_1^2 - " *
                "2*Wm_1*Wp_1 + Wm_1 - Wp_1)", [lg]; vars=Symbol[], cfg=cfg)
            @test SubTropica.canonical_text(e) == "-2 + 1/2*sqrt_disc_1"
        end
    end

    # ---- THE PERIOD UNIT GATE (RT-extra-period / gate A(v)) ------------------
    # Every residual ZeroInfPeriod key in the corpus: HF zero_inf_period op
    # (symbolic) vs INDEPENDENT ginac_gpl large-z Lagrange extrapolation
    # (numeric), ≥60 matched decimal digits. Semantics pin:
    # ZeroInfPeriod[w] = const term (log z -> 0) of G(w;z), z -> ∞.
    if have_hf && have_ginac
        @testset "period unit gate: HF vs ginac ≥60 dps" begin
            corpus = [["0", "-1"], ["0", "-1", "-1"], ["-1"], ["-2"],
                      ["0", "-2"], ["0", "0", "-1"], ["-1", "-2"],
                      ["0", "-1", "-2"]]
            for word in corpus
                hfe = SubTropica.zero_inf_period_value(word; cfg=cfg)
                hfv = SubTropica.hlog_const_bigfloat(hfe; digits=70)
                gv = SubTropica.zero_inf_period_ginac(word; digits=70)
                d = SubTropica._matched_digits(hfv, gv)
                @test d >= 60
                d >= 60 || @error "period gate FAIL" word SubTropica.canonical_text(hfe) d
            end
            # non-integer letters: HF refuses in-band ("Phase 6c scope")
            @test_throws SubTropica.HFError SubTropica.hf_call("zero_inf_period",
                Dict{String,Any}("word" => ["0", "-1/2"], "vars" => ["x"]); cfg=cfg)
        end

        @testset "synthetic-truth controls" begin
            # positive: extrapolator vs independent BigFloat truths (neither
            # HF nor the mzv dictionary involved)
            z2 = setprecision(BigFloat, 500) do
                BigFloat(pi)^2 / 6
            end
            lg2 = setprecision(BigFloat, 500) do
                log(BigFloat(2))
            end
            @test SubTropica._matched_digits(
                SubTropica.zero_inf_period_ginac(["0", "-1"]; digits=70), z2) >= 60
            @test SubTropica._matched_digits(
                SubTropica.zero_inf_period_ginac(["-2"]; digits=70), -lg2) >= 60
            # positive: mzv dictionary anchors
            @test SubTropica._matched_digits(SubTropica.mzv_bigfloat([2]; digits=70), z2) >= 60
            @test SubTropica._matched_digits(SubTropica.mzv_bigfloat([-1]; digits=70), -lg2) >= 60
            # ascending-convention consistency: ζ_H(1,2) = ζ(3) (Euler)
            @test SubTropica._matched_digits(SubTropica.mzv_bigfloat([1, 2]; digits=70),
                                          SubTropica.mzv_bigfloat([3]; digits=70)) >= 60
            # divergent mzv refused (outermost index +1)
            @test_throws SubTropica.SubTropicaEvalError SubTropica.mzv_bigfloat([2, 1])
            # NEGATIVE control: ×(1+δπ) perturbation (NON-rational — known
            # pitfall: rational perturbations alias into exact relations)
            # must FAIL the 60-digit bar and land near the injected 40 digits
            gv = SubTropica.zero_inf_period_ginac(["0", "-1"]; digits=70)
            pert = setprecision(BigFloat, 500) do
                gv * (1 + big(10.0)^-40 * BigFloat(pi))
            end
            dneg = SubTropica._matched_digits(pert, z2)
            @test dneg < 60
            @test 30 <= dneg <= 45
        end
    end

    # ---- comparator unit behavior -------------------------------------------
    @testset "matched-digits comparator" begin
        # precision-preserving (ambient-dps footgun): two 400-bit values
        # differing at 1e-100 must NOT count as equal at ambient 256 bits
        x = setprecision(BigFloat, 400) do
            BigFloat(pi)
        end
        y = setprecision(BigFloat, 400) do
            BigFloat(pi) * (1 + big(10.0)^-100)
        end
        @test 95 <= SubTropica._matched_digits(x, y) <= 105
        # absolute fallback against exact zero
        @test SubTropica._matched_digits(big(0.0), big(10.0)^-80) in (79, 80)
        @test SubTropica._matched_digits(big(1.0), big(1.0)) == 999
    end
end
