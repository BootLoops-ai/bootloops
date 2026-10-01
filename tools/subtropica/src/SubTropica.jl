# SubTropica.jl — module shell: the include/export wiring for every module.
# Shared types live in types.jl and b1_types.jl. Design: DESIGN.md;
# wire/handoff contracts: CONTRACTS.md.
#
# subtropica — no-Wolfram Julia rewrite of SubTropica's front-end driving the
# HyperFLINT C++ engine (v1.2.8 binary, READ-ONLY, CLI eval-json transport).
# Phase A scope: convergent, proven-LR Euler integrals from raw quadruples;
# check_divergences hard-ON; Euclidean region only.
# Portions translated from SubTropica.wl, (c) 2025-2026 S. Mizera, M. Giroux,
# G. Salvatori (MIT; see ../NOTICE and ../reference/LICENSE); cite arXiv:2604.20954.

module SubTropica

using JSON
using Libdl   # C-ABI in-process transport (hf_bridge.jl)
using Nemo
# DESIGN_B1 §4 step 1: the C4 polytope layer needs
# Oscar's polyhedral API. SELECTIVE import only — a blanket `using Oscar`
# makes names like `coefficients` ambiguous (Oscar/Hecke/Singular export
# different bindings than Nemo's) and broke Phase-A serialize.jl at suite
# time. Only polytope.jl consumes these; everything else stays pure Nemo.
import Oscar
using Oscar: newton_polytope, minkowski_sum, facets, vertices, affine_hull,
             normal_fan, halfspace_matrix_pair, convex_hull, normal_vector,
             negbias

# ---- shared types (changes here are package-level design decisions) ----------
include("types.jl")
include("b1_types.jl")     # B1 shared types (after types.jl — EpsExp ops)

# ---- HF boundary ---------------------------------------------------
include("hlogexpr.jl")     # token parse HF coef-string -> HlogExpr; canonical text form
include("serialize.jl")    # Julia -> Mma-grammar strings (round-trip law, CONTRACTS.md §d)
include("hf_bridge.jl")    # hf_call, version stamp+canary gate, error taxonomy, env contract

# ---- ε-expansion ----------------------------------------------------
include("eps_expand.jl")   # loggamma_series (RT-vii), exp(εΣb·logP) multinomial,
                           # sign-definiteness rule (RT-viii), Laurent pole guard (R2)
# (eps_expand.jl exports; loggamma_series is in the shared API export block below)
export gamma_prefactor_series, sign_split, eps_expand, laurent_pole_guard,
       SignSplitResult, SignIndefiniteError, EpsExpandRefusal

# ---- LR dispatch + assembly -----------------------------------------
include("lr_dispatch.jl")  # find_lr_orders (always `groups`), verify_order-FIELD replay,
                           # best_order persistence
include("assemble.jl")     # ct/order sums, letter-table union, normalization,
                           # LaurentSeries{HlogExpr} with provenance
# (lr_dispatch.jl/assemble.jl exports)
export find_lr_order, verify_lr_order, ensure_lr_order, persist_best_order,
       load_best_order, integrand_hash, best_order_path,
       find_lr_orders_scan, factor_table, derive_gauge_from_homogeneous_lr,
       NoLROrderError, LRDispatchError, LRReplayError,
       FaceResult, assemble, assemble_faces, union_letters, remap_letters,
       align_laurent, mzv_reduce, to_ginsh, ginsh_series_string,
       DivergenceCheckViolation, AssembleProvenanceError, LetterTableMismatch,
       LaurentTruncationError, GinshExportError, MzvReduceError

# ---- verification ---------------------------------------------------
include("verify.jl")       # oracle comparator {amflow, hiprec_sectordecomp, mpmath},
                           # pySecDec sanity-only, synthetic-truth controls,
                           # rational-point picker
# (verify.jl exports)
export OracleArtifact, load_fixture, load_oracle_artifact, parse_rational,
       pick_point, floored_matched_digits, compare_laurent, verify_fixture,
       run_oracle, mutate_coefficient, synthetic_controls, can_gate,
       eval_symbolic

# ---- B1/C build (DESIGN_B1 §4 step 2) ----------------
# Dependency order: b1_types.jl is included right after types.jl above;
# polytope (C4) → tropical (C3+C5) → subtract (C7, delegates trop/restrict to
# tropical.jl) → graph (C front-end) → lr_refine (gap tool) → b1_driver
# (C9/C11 divergent driver; needs eps_expand/lr_dispatch/hf_bridge/assemble,
# all included above).
include("polytope.jl")     # C4 Newton-polytope layer (Oscar first-class)
include("tropical.jl")     # TropI, ray classes, Σ_div, w-search, restrict
include("subtract.jl")     # C7 Möbius counterterm assembly (WARNING-1/2)
include("continue.jl")     # B2: C6 Nilsson-Passare continuation (consumes
                           # tropical.jl + subtract.jl surfaces; wired into the
                           # driver via subtropica_integrate allow_continuation —
                           # DESIGN_B1 §4c)
include("graph.jl")        # graph → U,F → EulerIntegrand front-end
include("lr_refine.jl")    # compatibility-graph LR letter refinement
include("b1_driver.jl")    # census routing + divergent face farm

# ---- default run-dir config ----
# The four `run_dir` keyword defaults (subtropica_integrate x2, subtropica_integrate_b1,
# subtropica_integrate_b2) route through here instead of repeating a hardcoded
# literal. Public default = ./subtropica_runs (cwd-relative);
# override with ENV["SUBTROPICA_RUNS"] — the same variable scripts/env.sh already
# exports. Resolved at CALL time (keyword-default lowering), so an env change
# takes effect without reloading the package.
const _DEFAULT_RUN_ROOT = "subtropica_runs"  # default (cwd-relative)
default_run_root() = get(ENV, "SUBTROPICA_RUNS", _DEFAULT_RUN_ROOT)
default_run_dir() = joinpath(default_run_root(), "subtropica_state")
export default_run_root, default_run_dir

# B1 shared types + classification (b1_types.jl)
export RayClass, convergent, log_divergent, power_divergent, classify_ray,
       warning1_sign, is_illdefined, is_divergent, is_power, is_log,
       TropicalData, DivergentFacet, SigmaFace, CounterTerm,
       B1RefusalException, B1Refusal, PowerDivergentRefusal,
       GeometricPropertyViolated, WNotUnitNormalized
# C4 (polytope.jl)
export tropical_data, primitive_integer, forget_coeffs,
       tropical_data_cache_key, refuse_if_degenerate, normal_fan_diagnostic
# C3+C5 (tropical.jl)
export support_x, trop_poly, trop_on_ray, trop_values, div_facets,
       lex_combinations, sigma_div, produce_ws, u_pair, one_minus_u_pair,
       restrict_poly, divergence_data
# C7 (subtract.jl)
export subtraction_terms
# B2/C6 (continue.jl) — Nilsson-Passare continuation surface
export continue_ray, continue_rays, find_first_np_continuation,
       expand_integral_b2, Eps2, Euler2Integrand, promote_second_regulator,
       B2Refusal, B2Prefactor, ContinuedIntegrand, B2Continued,
       b2pref_identity, b2pref_mul, b2pref_simplify, b2pref_eval,
       trop2_on_ray, trop2_values, is_promoted, unregulated_by_eps,
       epsexp_eval, eps2_eval
# graph front-end (graph.jl)
export symanzik_UF, graph_to_quadruple, quadruple_json, write_quadruple_json,
       load_quadruple, GraphRefusal
# LR refinement (lr_refine.jl)
export refine_letters, is_blocking_real, blocking_certificate, LRRefineRefusal
# B1 driver (b1_driver.jl) + B2 continuation route (same file)
export b1_ray_census, subtropica_integrate_b1, subtropica_integrate_b2
# two-reg driver glue (continue.jl §C6b) + the
# Log[rational]→prime-basis canonical-emission fold (now hlogexpr.jl —
# canonical_text folds at entry, canonical-atom law); the
# two-reg subtropica_integrate method itself rides the existing export.
export find_np_continuation_forced, fold_log_atoms
# word-level algebraic-letter boundary periods —
# numeric leg for HPeriod atoms over Wm/Wp roots (hlogexpr.jl), consumed
# by eval_symbolic via LaurentSeries.letters (verify.jl).
export alg_letter_roots, zero_inf_period_alg

# ---- exported API (stubs; signatures are the contract — DESIGN.md §6) ------
export EpsExp, Prefactor, EulerIntegrand,
       HAtom, HConst, HZeta, HLogA, HHlog, HPeriod, HAlgLetter, HDelta,
       HTerm, HlogExpr, AlgLetter, LaurentSeries, LogIntegrand,
       HFConfig, GateReport
export subtropica_integrate, euler_integrand, verify_laurent,
       hf_call, hf_gate, loggamma_series

"""
    euler_integrand(prefactor, polys, nu, vars; kinvars=Symbol[]) -> EulerIntegrand

Construct the raw Euler quadruple (THE Phase-A input, RT-x): integrand =
prefactor · ∏ x_i^nu_i · ∏ P_j^(a_j+b_j·ε) on [0,∞)^n, flat dx, all
variables and kinematic coefficients positive (declared-Euclidean).
`polys` = Vector of (poly, EpsExp); polys may be given as canonical
Mma-grammar strings (parsed into one Nemo ring, vars-then-kinvars order)
or as ready QQMPolyRingElem sharing one ring.
"""
function euler_integrand(prefactor::Prefactor, polys, nu, vars; kinvars=Symbol[])
    vs = Symbol.(collect(vars))
    ks = Symbol.(collect(kinvars))
    names = vcat(vs, ks)
    isempty(vs) && error("SubTropica.euler_integrand: no integration variables")
    # ring: either adopt the shared parent of supplied Nemo polys, or build one
    elems = [P for (P, _) in polys if P isa QQMPolyRingElem]
    local R::QQMPolyRing
    if !isempty(elems)
        R = parent(elems[1])
        all(parent(P) === R for P in elems) ||
            error("SubTropica.euler_integrand: supplied polys do not share one ring")
        nvars(R) == length(names) ||
            error("SubTropica.euler_integrand: ring has $(nvars(R)) vars ≠ |vars|+|kinvars|=$(length(names))")
    else
        R, _ = polynomial_ring(Nemo.QQ, String.(names))
    end
    gmap = Dict{Symbol,QQMPolyRingElem}(names[i] => gen(R, i) for i in eachindex(names))
    _eps(e) = e isa EpsExp ? e :
        EpsExp(Rational{BigInt}(parse_rational(string(e[1]))),
               Rational{BigInt}(parse_rational(string(e[2]))))
    ps = Tuple{QQMPolyRingElem,EpsExp}[]
    for (P, e) in polys
        Pp = P isa QQMPolyRingElem ? P : _ti_parse_poly(String(P), R, gmap)
        push!(ps, (Pp, _eps(e)))
    end
    nus = EpsExp[_eps(n) for n in nu]
    length(nus) == length(vs) ||
        error("SubTropica.euler_integrand: |nu|=$(length(nus)) ≠ |vars|=$(length(vs))")
    return EulerIntegrand(prefactor, nus, ps, vs, ks, R)
end

"""
    euler_integrand(fixture::AbstractDict) -> EulerIntegrand
    euler_integrand(fixture_path::AbstractString) -> EulerIntegrand

Build the quadruple from a `subtropica-euler-quad-v1` fixture (verify.jl schema,
fixtures/README.md). Carries NO result values by construction. The
`divergent` flag is NOT consulted here — divergence is detected LIVE by
check_divergences (hard-ON) in the pipeline, never assumed from metadata.
"""
function euler_integrand(fx::AbstractDict)
    pf = fx["prefactor"]
    gammas = Tuple{EpsExp,Int}[
        (EpsExp(parse_rational(string(g["arg"][1])),
                parse_rational(string(g["arg"][2]))), Int(g["power"]))
        for g in pf["gammas"]]
    prefactor = Prefactor(parse_rational(string(pf["c"])),
                          Int(get(pf, "eps_power", 0)), gammas,
                          parse_rational(string(pf["gammaE_eps"])))
    nu = [(string(n[1]), string(n[2])) for n in fx["nu"]]
    polys = [(String(p["poly"]), (string(p["exp"][1]), string(p["exp"][2])))
             for p in fx["polys"]]
    return euler_integrand(prefactor, polys, nu, Symbol.(fx["vars"]);
                           kinvars=Symbol.(fx["kinvars"]))
end
euler_integrand(path::AbstractString) = euler_integrand(load_fixture(path))

# -- fixture poly-string parser (integration glue; +,-,*,^ over ring gens) ---
function _ti_poly_eval(ex, R::QQMPolyRing, gmap::Dict{Symbol,QQMPolyRingElem})
    ex isa Integer && return R(Nemo.QQ(BigInt(ex)))
    if ex isa Symbol
        haskey(gmap, ex) || error("SubTropica.euler_integrand: unknown symbol '$ex' in fixture poly")
        return gmap[ex]
    end
    ex isa Expr || error("SubTropica.euler_integrand: unsupported token $(repr(ex)) in fixture poly")
    ex.head === :call || error("SubTropica.euler_integrand: unsupported expr head $(ex.head)")
    op = ex.args[1]
    as = ex.args[2:end]
    if op === :+
        return sum(_ti_poly_eval(a, R, gmap) for a in as)
    elseif op === :- && length(as) == 1
        return -_ti_poly_eval(as[1], R, gmap)
    elseif op === :- && length(as) == 2
        return _ti_poly_eval(as[1], R, gmap) - _ti_poly_eval(as[2], R, gmap)
    elseif op === :*
        return prod(_ti_poly_eval(a, R, gmap) for a in as)
    elseif op === :^
        (as[2] isa Integer && as[2] >= 0) ||
            error("SubTropica.euler_integrand: non-natural exponent in fixture poly")
        return _ti_poly_eval(as[1], R, gmap)^Int(as[2])
    elseif op === :/
        d = as[2]
        d isa Integer || error("SubTropica.euler_integrand: '/' only by integer constants in fixture polys")
        return _ti_poly_eval(as[1], R, gmap) * R(Nemo.QQ(1, BigInt(d)))
    end
    error("SubTropica.euler_integrand: unsupported operator '$op' in fixture poly")
end
_ti_parse_poly(s::AbstractString, R::QQMPolyRing, gmap) =
    _ti_poly_eval(Meta.parse(replace(s, r"\s+" => "")), R, gmap)

# -- factored HF expr for a LogIntegrand (integration glue) ------------------
# eps_expand builds num/den as EXPANDED Nemo products (types.jl contract);
# serializing the expanded polynomials blows up the HF CLI request (live
# the Smirnov 5-var eps^0 denominator expands to ~39 KB and the
# engine returns empty stdout rc=0). Upstream keeps products factored
# (STIntegrateSubtractionNP) — re-factor via Nemo before serialization.
# Mathematically identical to mma_string(li); bytes differ only in grouping.
function _ti_factored_poly_string(P::QQMPolyRingElem)
    is_constant(P) && return mma_string(P)
    fac = Nemo.factor(P)
    parts = String[]
    u = unit(fac)
    isone(u) || push!(parts, "(" * mma_string(u) * ")")
    for (f, m) in fac
        s = "(" * mma_string(f) * ")"
        push!(parts, m == 1 ? s : s * "^" * string(m))
    end
    return join(parts, "*")
end
function _ti_factored_expr(li::LogIntegrand)
    parts = String[]
    if isone(li.den)
        push!(parts, "(" * _ti_factored_poly_string(li.num) * ")")
    else
        push!(parts, "(" * _ti_factored_poly_string(li.num) * ")/(" *
                     _ti_factored_poly_string(li.den) * ")")
    end
    for (P, k) in li.logs
        k >= 1 || error("LogIntegrand log power must be ≥1, got $k")
        s = "Log[" * mma_string(P) * "]"
        push!(parts, k == 1 ? s : s * "^" * string(k))
    end
    return join(parts, "*")
end

"""
    subtropica_integrate(quad::EulerIntegrand; order=:auto, dps=60, gauge=:auto,
                      region=:euclidean, check_divergences=:auto,
                      findroots=:cascade, workers=1, hf=HFConfig())
        -> LaurentSeries{HlogExpr}

Symbolic Laurent-in-ε of a (Phase A: convergent, proven-LR) Euler integral
via HyperFLINT. `region=:physical` REFUSES in v1 (R4 quarantine).
`check_divergences=:auto` = hard-ON standalone (Phase A always), OFF only
inside Phase-B face-sum farms.

Continuation (B2): with the DEFAULT
`allow_continuation=false`, inputs needing Nilsson-Passare continuation
REFUSE loudly (typed PowerDivergentRefusal / GeometricPropertyViolated —
refusals stay the default surface, RT-ix). With `allow_continuation=true`
those two census refusals route into `expand_integral_b2` (dlog chart via
the `_b1_dlog` measure bridge); each continued integrand re-enters this
driver ONCE (allow_continuation=false inside — expand_integral_b2 derives
the full step ledger from the census, so a second pass is never
legitimate), its exact B2Prefactor is convolved in, and the term series
are pole-aligned and summed. `continuation_order` is the TOTAL
continuation STEP cap forwarded to `expand_integral_b2(; order=...)` —
NOT a Laurent order (DESIGN_B1 §4c).
"""
function subtropica_integrate(quad::EulerIntegrand; order=:auto, dps::Int=60,
                           gauge=:auto, region::Symbol=:euclidean,
                           check_divergences=:auto, findroots::Symbol=:cascade,
                           workers::Int=1, hf::HFConfig=HFConfig(),
                           run_dir::AbstractString=default_run_dir(),
                           domains=Dict{Symbol,Any}(), mzv::Bool=true,
                           strict_periods::Bool=true, on_period=nothing,
                           allow_continuation::Bool=false,
                           continuation_order::Integer=16)
    t0 = time()
    region === :euclidean ||
        error("SubTropica.subtropica_integrate: region=$region REFUSED in v1 — Euclidean only (DESIGN.md R4)")
    cd = check_divergences === :auto ? true : Bool(check_divergences)
    cd || throw(DivergenceCheckViolation(
        "subtropica_integrate: check_divergences=false is the Phase-B face-sum " *
        "farm ONLY (DESIGN.md C11) — standalone calls run the flag HARD-ON"))
    ord = order === :auto ? 2 : Int(order)
    ord >= 0 || error("SubTropica.subtropica_integrate: order=$ord < 0")

    # engine stamp + unary-minus canary (memoized; RT-vi)
    hf_gate(hf) || throw(HFGateError("hf_gate returned false — wrong/mis-stamped binary"))

    # ---- B1 ray census + routing (DESIGN_B1 §4) ----
    # The census is B1's OWN (SPEC_B1 §3.2; fixture metadata never consulted):
    # convergent / uncertifiable inputs continue on the Phase-A path UNCHANGED
    # below (HF check_divergences hard-ON remains that path's soundness bar);
    # census-certified log-divergent inputs take the B1 subtraction driver
    # (b1_driver.jl — check_divergences stays hard-ON there too). Power
    # divergence / GP violation / w-norm refusals propagate TYPED from here
    # — unless allow_continuation=true, in which case the two
    # continuation-shaped kinds (PowerDivergentRefusal /
    # GeometricPropertyViolated) route into the B2 continuation driver
    # (WNotUnitNormalized and all other refusals stay typed refusals).
    census = try
        b1_ray_census(quad)
    catch err
        if allow_continuation && (err isa PowerDivergentRefusal ||
                                  err isa GeometricPropertyViolated)
            return subtropica_integrate_b2(quad, err, ord;
                       hf=hf, run_dir=run_dir, domains=domains, mzv=mzv,
                       strict_periods=strict_periods, on_period=on_period,
                       findroots=findroots, workers=workers, dps=dps,
                       continuation_order=Int(continuation_order))
        end
        rethrow()
    end
    if census.route === :b1
        return subtropica_integrate_b1(quad, census.td, census.dd, census.Ed, ord;
                                    hf=hf, run_dir=run_dir, domains=domains,
                                    mzv=mzv, strict_periods=strict_periods,
                                    on_period=on_period, findroots=findroots,
                                    workers=workers, census_note=census.note)
    end

    # C8 (degenerate, Phase A): per-ε-order LogIntegrand lists + Γ-prefactor
    ex = eps_expand(quad, ord; domains=domains)

    # C10: LR order over the letter set (ensure = persisted search + replay, R3)
    letters = String[]
    for (P, _) in quad.polys
        is_constant(P) && continue
        for (f, _) in Nemo.factor(P)
            s = mma_string(f)
            s in letters || push!(letters, s)
        end
    end
    local lr_order::Vector{Symbol}
    lr_meta = Dict{String,Any}()
    alg_letters = false
    carry_flag = false
    if isempty(letters)
        lr_order = copy(quad.vars)     # monomial-only integrand: any order is LR
        lr_meta["note"] = "no non-constant letters — LR search skipped"
    else
        try
            (lr_order, lr_meta) = ensure_lr_order([letters], quad.vars,
                                                  String.(quad.kinvars);
                                                  run_dir=run_dir, cfg=hf)
        catch err
            # FindRoots-tier escalation (HFConfig: "escalation is the
            # cascade's job"): strict LR failed — retry allowing algebraic
            # letters (Wm/Wp/sqrt_disc), then + carried-sqrt discharge
            # (live pin against the v1.2.8 binary: the Moller-class quadratic F needs BOTH
            # algebraic_letters AND carry_discharge to find [x1,x0,x2]).
            (err isa NoLROrderError && findroots === :cascade) || rethrow()
            carry = false
            try
                (lr_order, lr_meta) = ensure_lr_order([letters], quad.vars,
                                                      String.(quad.kinvars);
                                                      run_dir=joinpath(run_dir, "algletters"),
                                                      cfg=hf, algebraic_letters=true)
            catch err2
                err2 isa NoLROrderError || rethrow()
                (lr_order, lr_meta) = ensure_lr_order([letters], quad.vars,
                                                      String.(quad.kinvars);
                                                      run_dir=joinpath(run_dir, "algletters_carry"),
                                                      cfg=hf, algebraic_letters=true,
                                                      carry_discharge=true)
                carry = true
                carry_flag = true
            end
            alg_letters = true
            lr_meta["findroots_tier"] = "algebraic_letters" *
                (carry ? "+carry_discharge" : "") * " (strict LR: NOLR)"
        end
    end

    # C9: HF integration per ε-order (check_divergences HARD-ON), residual
    # ZeroInfPeriod keys resolved (integer words → HF op; single kinematic
    # letter → -Log[-ℓ]; multi-letter kinematic → typed refusal)
    allvars = vcat(quad.vars, quad.kinvars)
    blocks = Vector{Tuple{HlogExpr,Vector{AlgLetter}}}[]
    hf_wall = 0.0
    for lst in ex.integrands
        entries = Tuple{HlogExpr,Vector{AlgLetter}}[]
        for li in lst
            expr = _ti_factored_expr(li)
            tI = time()
            (terms, tletters) = hf_integrate(expr, lr_order, allvars; cfg=hf,
                                             check_divergences=true,
                                             algebraic_letters=alg_letters,
                                             carry_discharge=carry_flag)
            hf_wall += time() - tI
            val = resolve_periods(terms; cfg=hf, letters=tletters,
                                  strict=strict_periods, on_period=on_period)
            push!(entries, (val, tletters))
        end
        push!(blocks, entries)
    end

    # C12: convolve with the Γ-prefactor Laurent, stamp provenance, MZV-reduce
    ev = Dict{String,Any}(
        "check_divergences" => true,
        "b1_census" => census.note,          # census certificate (route = :phase_a)
        "lr_order" => isempty(letters) ? "no_integration_required" : String.(lr_order),
        "lr_meta" => lr_meta,
        "wall_s_hf" => hf_wall,
        "wall_s_total_pre_assemble" => time() - t0,
        "run_dir" => String(run_dir),
        "requested_order" => ord,
        "prefactor_minorder" => ex.prefactor_minorder,
    )
    return assemble(ex.minorder, blocks,
                    (ex.prefactor_minorder, ex.prefactor_coeffs), ord;
                    evidence=ev, cfg=hf, mzv=mzv)
end

"""
    verify_laurent(L::LaurentSeries, oracle::Symbol; point=:auto,
                   min_digits::Int=30, sanity::Bool=false) -> GateReport

Compare a symbolic LaurentSeries against an independent numeric oracle at
an exact rational kinematic point (letter-degeneracy-avoiding picker).
Gate oracles (RT-iii): :amflow, :hiprec_sectordecomp, :mpmath.
:pysecdec is accepted ONLY with sanity=true and its GateReport can never
gate. Fixture-value comparison is FORBIDDEN in gates.
"""
function verify_laurent(L::LaurentSeries, oracle::Symbol; point=:auto,
                        min_digits::Int=30, sanity::Bool=false)
    error("SubTropica.verify_laurent: not implemented — src/verify.jl")
end

end # module SubTropica
