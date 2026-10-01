# test_wmwp_integration.jl — Wm/Wp integration tests.
# END-TO-END acceptance for the Wm/Wp algebraic-letter gap closure
# (CONTRACTS.md §d diagnosis + the RT-v chain wired into hf_integrate's
# coefficient path):
#
#   A. fresh minimal sqrt-letter fixture — the 1-loop massive-bubble class
#      at (m^2=1, s=-1): Feynman parameter x = t/(1+t) turns
#      ∫₀¹ dx/(1 - (s/m²)x(1-x)) into ∫₀^∞ dt/(t²+3t+1); the letter is
#      β = sqrt(1-4m²/s) = √5 (poly t²+3t+1, disc 5). HAND-DERIVED closed
#      form (partial fractions over the roots (-3∓√5)/2):
#         ∫₀^∞ dt/(t²+3t+1) = (1/√5)·log((3+√5)/(3-√5))
#      (sympy-checked: integrate(1/(t**2+3*t+1),(t,0,oo)) ==
#       -2*sqrt(5)*atanh(3*sqrt(5)/5)/5... == log((3+sqrt(5))/(3-sqrt(5)))/sqrt(5)).
#   B. the eq 4.11 k0 integral at the letters tier, checks hard-ON —
#      vs an independent dps-55 tanh-sinh quadrature reference value.
#   C. the eq 4.11 R'1 integral: checks-ON refusal is the DOCUMENTED
#      per-log-atom divergence-detector false positive (DESIGN_B1 §4b —
#      cross-log-atom infinity-tail cancellation; the dps-55 quadrature
#      converged), asserted typed; the stamped checks-OFF retry then goes
#      through, coefs normalize, and the value matches the dps-55
#      quadrature reference.
#
# Numeric route: coef HlogExprs (sqrt_disc atoms) evaluated exactly from the
# AlgLetter table; residual Wm/Wp boundary words (HPeriod atoms under
# strict=false) evaluated via the INDEPENDENT ginac_gpl CLI at complex
# letters with the same large-z Lagrange extrapolation as
# zero_inf_period_ginac (letters Wm/Wp = (sum ∓ sqrt(disc))/2, root-checked;
# conjugate-symmetry cache halves the calls). Bar: ≥30 floored digits, with
# a ×(1+1e-25π) mutation control (non-rational — anti-aliasing law).
#
# Self-contained: run via runtests.jl or
#   ulimit -v 32505856; julia --project=tools/subtropica tools/subtropica/test/test_wmwp_integration.jl

using Test
using SubTropica
using Nemo

@testset "wmwp integration (algebraic-letter gap closure)" begin

    cfg = SubTropica.HFConfig()
    have_hf = isfile(cfg.bin) && isfile(cfg.mzv_data_path)
    have_ginac = isfile(SubTropica._ginac_bin())
    have_hf || @warn "HyperFLINT engine not found — wmwp integration tests skipped"
    have_ginac || @warn "ginac_gpl not found — wmwp numeric gates skipped"

    DIG = 60          # comparator digits (the quadrature references are dps-55)
    GUARD = 40        # extrapolation guard (prototype: ~n·log10(15) loss)
    PREC = ceil(Int, (DIG + GUARD + 20) * 3.33) + 64

    # ---- numeric helpers (verification side; ginac_gpl is the independent
    # evaluator — same extrapolation scheme as zero_inf_period_ginac, complex
    # letters allowed by the CLI's re/im argument pairs) -----------------------

    ginacGc = (lets, zstr) -> begin
        args = String[string(DIG + GUARD), zstr, "0"]
        for (re, im) in lets
            push!(args, re, im)
        end
        r = SubTropica._spawn_capture(vcat(["bash", "-c",
            "ulimit -v 32505856; exec \"\$0\" \"\$@\"", SubTropica._ginac_bin()],
            args), nothing; timeout_s=600)
        lines = split(strip(r.stdout), '\n')
        length(lines) >= 2 || error("ginac_gpl bad output rc=$(r.rc): $(r.stderr)")
        setprecision(BigFloat, PREC) do
            Complex{BigFloat}(parse(BigFloat, replace(String(lines[1]), "E" => "e")),
                              parse(BigFloat, replace(String(lines[2]), "E" => "e")))
        end
    end

    ZCACHE = Dict{Vector{Tuple{String,String}},Complex{BigFloat}}()
    negim = lets -> [(re, startswith(im, "-") ? im[2:end] :
                          (im == "0" ? im : "-" * im)) for (re, im) in lets]
    function Zc(lets::Vector{Tuple{String,String}})
        all(l -> l == ("0", "0"), lets) && return Complex{BigFloat}(0)
        haskey(ZCACHE, lets) && return ZCACHE[lets]
        nl = negim(lets)
        # conjugate symmetry: all-real -> flipped-im letters give conj(G)
        haskey(ZCACHE, nl) && return (ZCACHE[lets] = conj(ZCACHE[nl]))
        setprecision(BigFloat, PREC) do
            n = length(lets)
            E0, S = 120, 8
            Ls = BigFloat[]
            gs = Complex{BigFloat}[]
            for j in 0:n
                push!(Ls, (E0 + S * j) * log(BigFloat(10)))
                push!(gs, ginacGc(lets, "1" * "0"^(E0 + S * j)))
            end
            c0 = Complex{BigFloat}(0)
            for j in 0:n
                w = one(BigFloat)
                for k in 0:n
                    k == j && continue
                    w *= (0 - Ls[k+1]) / (Ls[j+1] - Ls[k+1])
                end
                c0 += w * gs[j+1]
            end
            ZCACHE[lets] = c0
        end
    end

    pq = s -> begin
        p = split(s, "/")
        Rational{BigInt}(parse(BigInt, p[1]),
                         length(p) == 2 ? parse(BigInt, p[2]) : big(1))
    end

    # Wm/Wp/sqrt_disc values from the response table: monic law
    # Wm = (sum - sqrt(disc))/2, Wp = (sum + sqrt(disc))/2 (back_substitute
    # pin Wm-Wp -> -sqrt_disc); both ROOT-CHECKED against x² - sum·x + prod.
    function lettervals(l::SubTropica.AlgLetter)
        l.lc == "1" || error("non-monic letter (lc=$(l.lc)) — extend the law first")
        setprecision(BigFloat, PREC) do
            s = pq(l.sum); p = pq(l.product); D = pq(l.disc)
            bf = q -> BigFloat(numerator(q)) / BigFloat(denominator(q))
            sd = sqrt(Complex{BigFloat}(bf(D)))
            Wm = (bf(s) - sd) / 2
            Wp = (bf(s) + sd) / 2
            for W in (Wm, Wp)
                abs(W^2 - bf(s) * W + bf(p)) < big"1e-50" ||
                    error("AlgLetter root check FAILED for idx $(l.idx)")
            end
            (Wm, Wp, sd)
        end
    end

    cstr = x -> (string(real(x)), string(imag(x)))
    wordlets = (word, lv) -> Tuple{String,String}[
        haskey(lv, strip(w)) ? cstr(lv[strip(w)]) : (String(strip(w)), "0")
        for w in word]

    function atomval(a, lv, lvd)
        if a isa SubTropica.HConst
            a.name === :Pi && return Complex{BigFloat}(BigFloat(pi))
            a.name === :Log2 && return Complex{BigFloat}(log(BigFloat(2)))
            a.name === :I && return Complex{BigFloat}(0, 1)
            a.name === :EulerGamma &&
                return Complex{BigFloat}(BigFloat(Base.MathConstants.eulergamma))
            error("unhandled HConst $(a.name)")
        elseif a isa SubTropica.HZeta
            return Complex{BigFloat}(SubTropica.mzv_bigfloat(a.idx; digits=DIG + GUARD))
        elseif a isa SubTropica.HLogA
            q = SubTropica._parse_rat_literal(a.arg)
            q > 0 || error("non-positive Log arg $(a.arg)")
            return Complex{BigFloat}(log(BigFloat(numerator(q)) /
                                         BigFloat(denominator(q))))
        elseif a isa SubTropica.HAlgLetter
            return lvd[(a.kind, a.idx)]
        elseif a isa SubTropica.HPeriod
            return Zc(wordlets(a.word, lv))
        end
        error("unhandled atom $(typeof(a))")
    end

    # numeric value of a resolved constant HlogExpr; per-term values are
    # returned so the mutation control can perturb ONE term
    function evalterms(e::SubTropica.HlogExpr, lv, lvd)
        ctx = SubTropica._pctx(e.vars, SubTropica._letters_of(e))
        setprecision(BigFloat, PREC) do
            vals = Complex{BigFloat}[]
            for t in e.terms
                fr = SubTropica._parse_rat_frac(t.coef, ctx)
                nu = numerator(fr)
                de = denominator(fr)
                (is_constant(nu) && is_constant(de)) ||
                    error("non-constant coef '$(t.coef)'")
                q = constant_coefficient(nu) // constant_coefficient(de)
                v = Complex{BigFloat}(BigFloat(numerator(q)) /
                                      BigFloat(denominator(q)))
                for (a, k) in t.atoms
                    v *= atomval(a, lv, lvd)^k
                end
                push!(vals, v)
            end
            vals
        end
    end

    # full pipeline value: hf_integrate (letters tier) -> resolve_periods
    # strict=false -> numeric eval; returns (value, resolved, letters, terms)
    function pipeline_value(expr, order, vars; check_divergences=true)
        terms, letters = SubTropica.hf_integrate(expr, order, vars; cfg=cfg,
            check_divergences=check_divergences, algebraic_letters=true)
        length(letters) == 1 || error("expected exactly one letter, got $(length(letters))")
        res = SubTropica.resolve_periods(terms; cfg=cfg, letters=letters, strict=false)
        (Wm, Wp, sd) = lettervals(letters[1])
        lv = Dict("Wm_$(letters[1].idx)" => Wm, "Wp_$(letters[1].idx)" => Wp)
        lvd = Dict((:Wm, letters[1].idx) => Wm, (:Wp, letters[1].idx) => Wp,
                   (:sqrt_disc, letters[1].idx) => sd)
        vals = evalterms(res, lv, lvd)
        (sum(vals; init=Complex{BigFloat}(0)), res, letters, vals)
    end

    imag_ok = v -> abs(imag(v)) <=
        max(abs(real(v)), one(BigFloat)) * big(10.0)^(-(DIG - 8))

    if have_hf && have_ginac
        # ------------------------------------------------------------------
        # A. massive-bubble sqrt-letter fixture (fresh, minimal)
        # ------------------------------------------------------------------
        @testset "massive bubble: sqrt(1-4m2/s) letter end-to-end" begin
            (v, res, letters, vals) = pipeline_value("1/(x^2+3*x+1)", [:x], [:x])
            # the parsed artifact carries the normalized algebraic atoms
            @test occursin("sqrt_disc_1", SubTropica.canonical_text(res))
            @test occursin("ZeroInfPeriod", SubTropica.canonical_text(res))
            @test letters[1].disc == "5"
            # hand-derived closed form: (1/√5)·log((3+√5)/(3−√5))
            ref = setprecision(BigFloat, PREC) do
                s5 = sqrt(BigFloat(5))
                log((3 + s5) / (3 - s5)) / s5
            end
            @test imag_ok(v)
            d = SubTropica._matched_digits(real(v), ref)
            @test d >= 30
            @info "massive bubble vs closed form: $d floored digits (bar 30)"
            # mutation control (×(1+1e-25π), non-rational — anti-alias law):
            # perturbing ONE term must drop below the bar
            mut = setprecision(BigFloat, PREC) do
                vm = copy(vals)
                vm[1] *= (1 + big"1e-25" * BigFloat(pi))
                real(sum(vm))
            end
            @test SubTropica._matched_digits(mut, ref) < 30
        end

        # ------------------------------------------------------------------
        # B. eq 4.11 k0 (vs the dps-55 quadrature reference)
        # ------------------------------------------------------------------
        @testset "eq 4.11 k0 at the letters tier vs reference dps-55" begin
            k0expr = "-2*Log[x1]/((x1*x2+x1+x2)*(2*x1+x2+1))"
            (v, res, letters, _) = pipeline_value(k0expr, [:x1, :x2], [:x1, :x2])
            @test letters[1].poly == "x2^2 + 1" && letters[1].disc == "-4"
            @test occursin("sqrt_disc_1", SubTropica.canonical_text(res))
            ref = setprecision(BigFloat, PREC) do
                parse(BigFloat,
                      "1.2697979381877088371491554851606652594605627130925857")
            end
            @test imag_ok(v)
            d = SubTropica._matched_digits(real(v), ref)
            @test d >= 30
            @info "eq 4.11 k0 vs reference dps-55 quadrature: $d floored digits (bar 30)"
        end

        # ------------------------------------------------------------------
        # C. eq 4.11 R'1: typed divergence-detector false positive first,
        #    then the stamped checks-OFF leg vs the dps-55 reference value
        # ------------------------------------------------------------------
        @testset "eq 4.11 R'1: typed FP refusal + checks-OFF leg vs reference" begin
            r1expr = "Log[x1]*Log[x1*x2+x1+x2]*(x1+x2)/(x1*x2+x1+x2)^2" *
                     "-2*Log[x1]*Log[2*x1+x2+1]*(x1+x2)/(x1*x2+x1+x2)^2" *
                     "+2*Log[x1]*Log[2]*(x1+x2)/(x1*x2+x1+x2)^2" *
                     "+Log[x1]*Log[x2]*1/((1+x2)*(1+x1)^2)" *
                     "-2*Log[x1]*Log[2]*1/((1+x2)*(1+x1)^2)" *
                     "-Log[x1]*Log[1+x1]*1/((1+x2)*(1+x1)^2)" *
                     "+Log[x1]^2*1/((1+x1)*(1+x2)^2)" *
                     "-Log[x1]*Log[1+x2]*1/((1+x1)*(1+x2)^2)"
            # checks hard-ON refuses typed: the engine's per-log-atom rational
            # pre-check cannot see the cross-log-atom infinity-tail
            # cancellation (DESIGN_B1 §4b false-positive class; the
            # dps-55 quadrature of this integrand CONVERGED)
            @test_throws SubTropica.HFDivergent SubTropica.hf_integrate(
                r1expr, [:x1, :x2], [:x1, :x2]; cfg=cfg,
                check_divergences=true, algebraic_letters=true)
            # stamped checks-OFF retry (the CONTRACTS §a face-sum situation;
            # external verification = the reference value compared right here)
            (v, res, letters, _) = pipeline_value(r1expr, [:x1, :x2],
                [:x1, :x2]; check_divergences=false)
            @test letters[1].poly == "x2^2 + 1"
            ref = setprecision(BigFloat, PREC) do
                parse(BigFloat,
                      "-6.601432780780210667590251167035173355892626833219107")
            end
            @test imag_ok(v)
            d = SubTropica._matched_digits(real(v), ref)
            @test d >= 30
            @info "eq 4.11 R'1 vs reference dps-55 quadrature: $d floored digits (bar 30)"
        end
    end
end
