# b1_driver.jl — B1 DIVERGENT-PATH DRIVER (C9/C11 farm slot). Integration-
# owned (
# DESIGN_B1 §4 step 3).
# ===========================================================================
# The STExpandIntegral-equivalent for divergent inputs (wl:11029-11145):
#
#   quadruple ──census──> tropical_data (C4) ──> divergence_data (C5:
#   refusal ladder §6 — power ⇒ PowerDivergentRefusal [names B2 via the
#   _NP_SUFFIX message law], GP ⇒ GeometricPropertyViolated, w-norm ⇒
#   WNotUnitNormalized) ──> produce_ws ──> subtraction_terms (C7)
#   ──> per-counterterm eps_expand (C8) ──> per-FACE LR search +
#   hf_integrate on the face's PER-ORDER COUNTERTERM SUM (theorem eq 3.33:
#   each face's Möbius sum is locally finite order-by-order, so
#   check_divergences stays HARD-ON — HF's divergence detector is a live
#   negative control on the subtraction, not a disabled safety) ──>
#   assemble_faces with cross-face pole alignment (assemble.jl align_laurent
#   — the general fold written for exactly this).
#
# ROUTING (wired into subtropica_integrate): convergent inputs take the
# Phase-A path UNCHANGED; the B1 path only engages after ITS OWN ray
# census (b1_ray_census below — fixture metadata is never consulted).
#
# MEASURE BRIDGE (types.jl nu law, design decision): the quadruple
# boundary carries FLAT-measure nu; the B1 tropical layer reads nu in the
# DLOG chart — `_b1_dlog` shifts ν_dlog = ν_flat + 1 on entry. CounterTerm
# output is already back in FLAT measure (subtract.jl mono = amb − 1 after
# the /∏x_J shift, wl:10837), so C8 consumes ct.mono directly.
#
# FARM-SLOT DECISION (recorded): faces are independent;
# the transport is the HyperFLINT CLI (one subprocess per hf_integrate call,
# exactly the Phase-A pattern), and the driver runs the face loop SERIALLY
# in-process — Oscar/Nemo state is not thread-safe and at B1 gate scales
# (≤ ~10 faces, fixtures ≪ 1 s/face) process-farm dispatch would be pure
# overhead. The seam for scaling out is the per-face body `_b1_face_result`
# (self-contained; scripts/fibrate-style per-face job JSONs + the C-ABI
# worker are the named upgrade when a production
# target needs it). `workers` kwarg is accepted and recorded in evidence;
# values > 1 do not change semantics yet.
# ===========================================================================

# ---------------------------------------------------------------------------
# Census helpers
# ---------------------------------------------------------------------------

"""ν_dlog = ν_flat + 1 view of the quadruple (measure bridge; SPEC_B1 §0)."""
_b1_dlog(E::EulerIntegrand) =
    EulerIntegrand(E.prefactor, EpsExp[ν + 1 for ν in E.nu], E.polys,
                   E.vars, E.kinvars, E.ring)

"""Coefficient-free, kinvar-projected factor list over QQ[vars...] for the
C4 polytope build (SPEC_B1 §2.1: the polytope lives over the integration
variables ONLY — kinematic ring gens are coefficients upstream and their
axes must not enter; coefficients all → 1 so monomial collapse cannot
cancel). Pure-constant projections (kinvar-only factors) are dropped
(wl:10415). Accepts the two-regulator quadruple too
— only poly SUPPORTS are read, regulator type is irrelevant here."""
function _b1_census_polys(E::Union{EulerIntegrand,Euler2Integrand})
    n = length(E.vars)
    Rv, _ = polynomial_ring(Nemo.QQ, String.(E.vars))
    out = QQMPolyRingElem[]
    for (P, _) in E.polys
        iszero(P) && throw(ArgumentError("b1 census: zero polynomial factor"))
        sup = Set{Vector{Int}}()
        for i in 1:length(P)
            push!(sup, Int.(exponent_vector(P, i)[1:n]))
        end
        (length(sup) == 1 && all(iszero, first(sup))) && continue
        B = MPolyBuildCtx(Rv)
        for ev in sup
            push_term!(B, one(Nemo.QQ), ev)
        end
        push!(out, finish(B))
    end
    return out
end

"""
    b1_ray_census(E::EulerIntegrand) -> NamedTuple

B1's OWN ray census (SPEC_B1 §3.2) — the routing decision for
subtropica_integrate. Returns `route = :phase_a` (with `note`) when the B1
route must NOT engage:
  * no non-constant poly factors over the integration vars (no polytope —
    monomial-only integrands stay on the Phase-A path, where HF
    check_divergences hard-ON is the soundness bar), or
  * census refusal :TropIllDefined / :DegeneratePolytope (unregulated or
    lineality-degenerate inputs — the trop criterion cannot certify them;
    Phase-A HF check_divergences catches genuine divergence loudly, so the
    historical Phase-A behavior is preserved byte-for-byte), or
  * all rays convergent (locally finite — Phase-A path unchanged).
Returns `route = :b1` with (td, dd, Ed) when the census certifies ≥1
log-divergent ray and the full C5 ladder passed. PowerDivergentRefusal,
GeometricPropertyViolated, WNotUnitNormalized and
B1Refusal(:SigmaDivTooLarge/:NoBasisCompletion) PROPAGATE typed (SPEC_B1
§6 contract — never a silent fall-through on a divergence the census DID
certify).
"""
function b1_ray_census(E::EulerIntegrand)
    cps = _b1_census_polys(E)
    isempty(cps) && return (route = :phase_a, td = nothing, dd = nothing,
        Ed = nothing,
        note = "b1 census: no non-constant poly factors over the integration vars — no polytope; Phase-A path")
    Ed = _b1_dlog(E)
    local td, dd
    try
        td = tropical_data(cps, copy(E.vars))
        dd = divergence_data(Ed, td)
    catch err
        if err isa B1Refusal && err.kind in (:TropIllDefined, :DegeneratePolytope)
            return (route = :phase_a, td = nothing, dd = nothing, Ed = nothing,
                note = "b1 census refusal $(err.kind) — trop criterion cannot certify; Phase-A path (HF check_divergences hard-ON is the soundness bar)")
        end
        rethrow()
    end
    if isempty(dd.div_facets)
        return (route = :phase_a, td = td, dd = nothing, Ed = nothing,
            note = "b1 census: $(length(td.rays)) rays, 0 divergent — locally finite; Phase-A path")
    end
    return (route = :b1, td = td, dd = dd, Ed = Ed,
        note = "b1 census: $(length(td.rays)) rays, $(length(dd.div_facets)) log-divergent, $(length(dd.sigma_div)) Σ_div faces — B1 subtraction path")
end

# ---------------------------------------------------------------------------
# Per-counterterm reduction to the face chart
# ---------------------------------------------------------------------------

"""Face-chart ring + gauge substitution of one CounterTerm: reduced ring
QQ[face_vars..., kinvars...] and the ambient→reduced substitution vector
(vars ∈ J → their reduced gens, vars ∉ J → 1, kinvars → trailing gens).
Split out of `_b1_ct_quad` so the two-reg route can map
INSERTED log arguments through the identical substitution."""
function _b1_ct_chart(E::EulerIntegrand, ct::CounterTerm)
    n  = length(E.vars)
    nf = length(ct.face_vars)
    nk = length(E.kinvars)
    Rf, _ = polynomial_ring(Nemo.QQ, vcat(String.(ct.face_vars), String.(E.kinvars)))
    posJ = Dict{Int,Int}(j => k for (k, j) in enumerate(ct.J))
    subs = Vector{QQMPolyRingElem}(undef, nvars(E.ring))
    for i in 1:n
        subs[i] = haskey(posJ, i) ? gen(Rf, posJ[i]) : one(Rf)
    end
    for k in 1:nk
        subs[n+k] = gen(Rf, nf + k)
    end
    return (Rf, subs)
end

"""Reduced-ring quadruple of one CounterTerm: ring QQ[face_vars..., kinvars...],
vars = ct.face_vars, nu = ct.mono (FLAT measure — the /∏x_J shift is already
inside), polys = the gauge-fixed restricted P's + jacobian factors mapped
into the reduced ring. Prefactor(1) — sign/Vol/global prefactor are applied
at the face/assembly layer, never here (no double counting)."""
function _b1_ct_quad(E::EulerIntegrand, ct::CounterTerm)
    (Rf, subs) = _b1_ct_chart(E, ct)
    polys = Tuple{QQMPolyRingElem,EpsExp}[(evaluate(P, subs), e) for (P, e) in ct.polys]
    return EulerIntegrand(Prefactor(1//1), copy(ct.mono), polys,
                          copy(ct.face_vars), copy(E.kinvars), Rf)
end

_b1_qstr(q::Rational{BigInt}) =
    denominator(q) == 1 ? string(numerator(q)) :
                          string(numerator(q)) * "/" * string(denominator(q))

"""Zero-variable face path (|J| = 0, i.e. |face| = n): every gauge-fixed
poly is constant, the ct value is sign·∏c^(a+bε) and the ε-series is exact —
no eps_expand ring, no HF call (the wl:11047-11053 zero-var branch). Returns
`orders`: Vector of ε^0..ε^K coefficients as HlogExpr (Log[c] atoms).
Kinvar-dependent 0-var polys REFUSE loudly (not needed by any B1 gate
target; |face| = n only occurs for fixture-scale n)."""
function _b1_const_ct_orders(ct::CounterTerm, mul::Rational{BigInt}, K::Int)
    q0 = mul
    lets = Tuple{Rational{BigInt},Rational{BigInt}}[]   # (c, b) with c^(a+bε)
    for (P, e) in ct.polys
        is_constant(P) || throw(EpsExpandRefusal(
            "b1 driver: zero-variable face $(ct.face) carries a NON-CONSTANT " *
            "poly after gauge-fix ($(P)) — kinvar-dependent 0-var faces are " *
            "outside the B1 gate scope; REFUSING"))
        c = Rational{BigInt}(BigInt(numerator(leading_coefficient(P))),
                             BigInt(denominator(leading_coefficient(P))))
        c > 0 || throw(EpsExpandRefusal(
            "b1 driver: 0-var face constant $c ≤ 0 (Euclidean violation)"))
        denominator(e.a) == 1 || throw(EpsExpandRefusal(
            "b1 driver: 0-var face constant exponent $(e.a) non-integer"))
        q0 *= c^Int(numerator(e.a))
        (iszero(e.b) || isone(c)) || push!(lets, (c, e.b))
    end
    T = length(lets)
    out = HlogExpr[]
    for o in 0:K
        terms = HTerm[]
        for k in _compositions(o, T)
            (T == 0 && o > 0) && continue
            q = q0
            atoms = Pair{HAtom,Int}[]
            for t in 1:T
                k[t] == 0 && continue
                q *= lets[t][2]^k[t] // factorial(big(k[t]))
                push!(atoms, HLogA(_b1_qstr(lets[t][1])) => k[t])
            end
            iszero(q) || push!(terms, HTerm(_b1_qstr(q), _canon_atomvec(atoms)))
        end
        push!(out, isempty(terms) ? HlogExpr() : HlogExpr(terms, Symbol[]))
    end
    return out
end

# ---------------------------------------------------------------------------
# Per-face LR cascade (mirror of the Phase-A tier escalation in
# subtropica_integrate; multi-group law: per-counterterm letter GROUPS for the
# ct sum — SubTropica wl:15659-15691 initial pair sets, DESIGN_B1 §2 lr_refine.jl
# row "groups ALWAYS for sums").
# ---------------------------------------------------------------------------
function _b1_lr_cascade(groups::Vector{Vector{String}}, xvars::Vector{Symbol},
                        kin::Vector{String}; run_dir::AbstractString,
                        hf::HFConfig, findroots::Symbol)
    local lr_order::Vector{Symbol}
    lr_meta = Dict{String,Any}()
    alg = false
    carry = false
    try
        (lr_order, lr_meta) = ensure_lr_order(groups, xvars, kin;
                                              run_dir=run_dir, cfg=hf)
    catch err
        (err isa NoLROrderError && findroots === :cascade) || rethrow()
        try
            (lr_order, lr_meta) = ensure_lr_order(groups, xvars, kin;
                                                  run_dir=joinpath(run_dir, "algletters"),
                                                  cfg=hf, algebraic_letters=true)
        catch err2
            err2 isa NoLROrderError || rethrow()
            (lr_order, lr_meta) = ensure_lr_order(groups, xvars, kin;
                                                  run_dir=joinpath(run_dir, "algletters_carry"),
                                                  cfg=hf, algebraic_letters=true,
                                                  carry_discharge=true)
            carry = true
        end
        alg = true
        lr_meta["findroots_tier"] = "algebraic_letters" *
            (carry ? "+carry_discharge" : "") * " (strict LR: NOLR)"
    end
    return (lr_order, lr_meta, alg, carry)
end

# ---------------------------------------------------------------------------
# One face → FaceResult (the farm unit; self-contained)
# ---------------------------------------------------------------------------
function _b1_face_result(E::EulerIntegrand, cts::Vector{CounterTerm}, fidx::Int,
                         K::Int; hf::HFConfig, run_dir::AbstractString,
                         domains, strict_periods::Bool, on_period,
                         findroots::Symbol, facemeta::Dict{String,Any})
    ct0 = cts[1]
    d = length(ct0.face)
    # face Vol = vol_det / ∏(−TropI(ρ)) with −TropI = vt.b·ε on log rays ⇒
    # Vol = q_face · ε^(−d); q_face exact rational (b1_types CounterTerm doc)
    q_face = Rational{BigInt}(ct0.vol_det)
    for vt in ct0.vol_trops
        iszero(vt.a) || throw(AssertionError(
            "b1 driver: vol_trop with nonzero constant part on face $(ct0.face) — non-log ray leaked past the C5 ladder"))
        iszero(vt.b) && throw(AssertionError(
            "b1 driver: zero vol_trop on face $(ct0.face)"))
        q_face //= vt.b
    end
    facemeta["face"] = string(ct0.face)
    facemeta["depth"] = d
    facemeta["q_face"] = _b1_qstr(q_face)
    facemeta["n_cts"] = length(cts)

    # ---- zero-variable face: exact constant series, no HF ------------------
    if isempty(ct0.face_vars)
        acc = [HlogExpr() for _ in 0:K]
        for ct in cts
            os = _b1_const_ct_orders(ct, q_face * ct.sign, K)
            for k in 1:K+1
                acc[k] = hlog_add(acc[k], os[k])
            end
        end
        facemeta["lr_order"] = "no_integration_required"
        orders = [isempty(a.terms) ? Tuple{HlogExpr,Vector{AlgLetter}}[] :
                  Tuple{HlogExpr,Vector{AlgLetter}}[(a, AlgLetter[])] for a in acc]
        return FaceResult(-d, orders)
    end

    # ---- eps-expand every counterterm on the face chart --------------------
    quads = [_b1_ct_quad(E, ct) for ct in cts]
    exs = [eps_expand(q, K; domains=domains) for q in quads]

    # ---- per-ct letter groups (multi-group law) + LR cascade ---------------
    # LETTER RULE (eq 4.3 gate): a poly's factors are LR letters
    # iff it can appear as a DENOMINATOR or LOG argument in some ε-order —
    # e.b ≠ 0 (log letters) or e.a < 0 (denominator). Pure-numerator factors
    # (e.a ≥ 0 integer, e.b == 0, e.g. the eikonal tensor numerator Q_N with
    # exponent +1) never enter the integrand as letters; feeding them to
    # find_lr_orders NOLRs solvable faces (live: eq 4.3 ∅-face fell through
    # the whole FindRoots cascade and the algletters_carry 5-var call died,
    # while upstream reports a strict-LR order for the same counterterms).
    groups = Vector{String}[]
    for q in quads
        g = String[]
        for (P, e) in q.polys
            is_constant(P) && continue
            (iszero(e.b) && e.a >= 0) && continue   # pure numerator — no letter
            for (f, _) in Nemo.factor(P)
                s = mma_string(f)
                s in g || push!(g, s)
            end
        end
        isempty(g) || push!(groups, g)
    end
    allvars = vcat(ct0.face_vars, E.kinvars)
    local lr_order::Vector{Symbol}
    alg = false; carry = false
    if isempty(groups)
        lr_order = copy(ct0.face_vars)
        facemeta["lr_order"] = "no non-constant letters - LR search skipped"
    else
        (lr_order, lr_meta, alg, carry) = _b1_lr_cascade(groups, ct0.face_vars,
            String.(E.kinvars); run_dir=run_dir, hf=hf, findroots=findroots)
        facemeta["lr_order"] = String.(lr_order)
        facemeta["lr_meta"] = lr_meta
    end

    # ---- per-order: SUM the counterterms, ONE hf_integrate per order -------
    # (theorem eq 3.33: the face's Möbius sum is order-by-order locally
    # finite — check_divergences stays hard-ON and doubles as the live
    # negative control on the subtraction itself)
    orders = Vector{Vector{Tuple{HlogExpr,Vector{AlgLetter}}}}()
    hf_wall = 0.0
    for k in 0:K
        parts = String[]
        for (ct, ex) in zip(cts, exs)
            mul = q_face * ct.sign          # exact rational, folded into num
            for li in ex.integrands[k+1]
                sli = LogIntegrand(li.num * parent(li.num)(Nemo.QQ(numerator(mul),
                                                                   denominator(mul))),
                                   li.den, li.logs, li.vars, li.kinvars)
                push!(parts, "(" * _ti_factored_expr(sli) * ")")
            end
        end
        if isempty(parts)
            push!(orders, Tuple{HlogExpr,Vector{AlgLetter}}[])
            continue
        end
        expr = join(parts, "+")
        tI = time()
        (terms, tlet) = try
            hf_integrate(expr, lr_order, allvars; cfg=hf,
                         check_divergences=true,
                         algebraic_letters=alg,
                         carry_discharge=carry)
        catch err
            err isa HFDivergent || rethrow()
            # DETECTOR FALSE-POSITIVE CLASS (measured live, eq 4.3
            # gate, face [2,4,10] k=1): the engine's divergence pre-check is
            # per-log-atom rational and cannot certify TRANSCENDENTAL tail
            # cancellation across different log atoms (Log[x6] vs Log[x6+1]
            # infinity-boundary tails; the ct-SUM decays like x6^-2·log —
            # verified numerically on the eq 4.3 gate). Theorem
            # eq 3.33 + the census (log-only rays, GP certified) guarantee the
            # face's per-order Möbius sum is locally finite, so this is the
            # CONTRACTS §a face-sum situation the check-OFF policy was written
            # for. Fall back to check OFF for THIS order only, STAMPED in
            # facemeta/evidence (never silent); the acceptance-gate quadrature
            # and published-value legs verify the fallback values externally.
            fb = get!(facemeta, "divergence_check_fallbacks", Any[])
            push!(fb, Dict{String,Any}("order" => k,
                  "message" => sprint(showerror, err)))
            hf_integrate(expr, lr_order, allvars; cfg=hf,
                         check_divergences=false,
                         algebraic_letters=alg,
                         carry_discharge=carry)
        end
        hf_wall += time() - tI
        val = resolve_periods(terms; cfg=hf, letters=tlet,
                              strict=strict_periods, on_period=on_period)
        push!(orders, Tuple{HlogExpr,Vector{AlgLetter}}[(val, tlet)])
    end
    facemeta["wall_s_hf"] = hf_wall
    return FaceResult(-d, orders)
end

# ---------------------------------------------------------------------------
# THE DIVERGENT ENTRY POINT (called from subtropica_integrate after the census)
# ---------------------------------------------------------------------------
"""
    subtropica_integrate_b1(E, td, dd, Ed, requested_order; ...) -> LaurentSeries

B1 divergent pipeline (header comment above). `td`/`dd`/`Ed` come from
`b1_ray_census` (route = :b1). All HF calls run check_divergences hard-ON.
"""
function subtropica_integrate_b1(E::EulerIntegrand, td::TropicalData, dd, Ed,
                              requested_order::Int;
                              hf::HFConfig=HFConfig(),
                              run_dir::AbstractString=default_run_dir(),
                              domains=Dict{Symbol,Any}(), mzv::Bool=true,
                              strict_periods::Bool=true, on_period=nothing,
                              findroots::Symbol=:cascade, workers::Int=1,
                              census_note::AbstractString="")
    t0 = time()
    # C7 — Möbius counterterm assembly (dlog-chart input; FLAT-measure out)
    cts = subtraction_terms(Ed, td, dd.sigma_div, dd.w)

    # group by face, sigma_div order preserved (C7 emits in that order)
    faces_cts = Vector{CounterTerm}[]
    for ct in cts
        if !isempty(faces_cts) && faces_cts[end][1].face == ct.face
            push!(faces_cts[end], ct)
        else
            push!(faces_cts, CounterTerm[ct])
        end
    end
    @assert length(faces_cts) == length(dd.sigma_div) "b1 driver: face grouping does not match sigma_div"

    dmax = maximum(length(g[1].face) for g in faces_cts)
    gmin = -dmax
    # global prefactor series to the depth the convolution needs
    (pmin0, _) = gamma_prefactor_series(E.prefactor, 0)
    parg = max(0, requested_order - gmin - pmin0)
    (pmin_g, pcoeffs_g) = gamma_prefactor_series(E.prefactor, parg)
    # per-face integrand top order (uniform top after alignment):
    # face coeff at ε^(k−d) needs ct-integrand order k ≤ K = req − pmin + d
    faces = FaceResult[]
    facemetas = Dict{String,Any}[]
    for (fidx, g) in enumerate(faces_cts)
        d = length(g[1].face)
        K = requested_order - pmin_g + d
        K >= 0 || throw(LaurentTruncationError(
            "b1 driver: requested_order=$requested_order below the face pole depth bookkeeping (K=$K < 0)"))
        fm = Dict{String,Any}()
        push!(faces, _b1_face_result(E, g, fidx, K; hf=hf,
              run_dir=joinpath(run_dir, "b1_face_$(fidx)"), domains=domains,
              strict_periods=strict_periods, on_period=on_period,
              findroots=findroots, facemeta=fm))
        push!(facemetas, fm)
    end

    ev = Dict{String,Any}(
        "check_divergences" => true,          # hard-ON on every B1 HF call
        "lr_order" => Dict{String,Any}(fm["face"] => fm["lr_order"]
                                       for fm in facemetas),
        "b1_path" => true,
        "b1_census" => String(census_note),
        "b1_faces" => facemetas,
        "b1_n_counterterms" => length(cts),
        "b1_workers" => workers,              # serial v1 (farm-slot decision, header)
        "run_dir" => String(run_dir),
        "requested_order" => requested_order,
        "prefactor_minorder" => pmin_g,
        "wall_s_total_pre_assemble" => time() - t0,
    )
    return assemble_faces(faces, (pmin_g, pcoeffs_g), requested_order;
                          evidence=ev, cfg=hf, mzv=mzv)
end

# ===========================================================================
# B2 CONTINUATION ROUTE (DESIGN_B1 §4c).
# Engaged from subtropica_integrate ONLY when allow_continuation=true and the
# B1 census refused PowerDivergentRefusal / GeometricPropertyViolated.
# Chart law: continuation runs in the DLOG chart (`_b1_dlog` bridge, exactly
# as continue.jl's header prescribes); continued integrands bridge BACK to
# the FLAT quadruple boundary (`_b1_undlog`) and re-enter subtropica_integrate
# with allow_continuation=false — ONE continuation pass, matching upstream
# STExpandIntegral (the NP block runs once, wl:11093-11135; the step ledger
# already covers GP-removal AND power-shaving multiplicities). A continued
# output that still refuses propagates typed (honest, never looped).
# ===========================================================================

"""FLAT-measure view of a dlog-chart integrand: ν_flat = ν_dlog − 1 (inverse
of `_b1_dlog`; types.jl nu MEASURE LAW)."""
_b1_undlog(E::EulerIntegrand) =
    EulerIntegrand(E.prefactor, EpsExp[ν - 1 for ν in E.nu], E.polys,
                   E.vars, E.kinvars, E.ring)

"""ε-order shift of a B2Prefactor's Laurent expansion: +1 per ε-proportional
num factor, −1 per ε-proportional den factor (the continuation poles).
Zero factors are LOUD (a zero num factor is a vanished term upstream — a
continue.jl bug, not a value; a zero den factor is 1/0)."""
function _b2pref_poleshift(p::B2Prefactor{EpsExp})
    m = 0
    for x in p.num
        iszero(x) && throw(ArgumentError(
            "b2 prefactor: zero numerator factor — vanished term should not have been emitted"))
        iszero(x.a) && (m += 1)
    end
    for x in p.den
        iszero(x) && throw(ArgumentError("b2 prefactor: zero denominator factor (1/0)"))
        iszero(x.a) && (m -= 1)
    end
    return m
end

"""
    _b2pref_rational_series(p::B2Prefactor{EpsExp}, K) -> (m, r)

EXACT rational Laurent expansion p(ε) = Σ_{k=0}^{K} r[k+1]·ε^(m+k) + O(ε^(m+K+1))
with m = `_b2pref_poleshift(p)`. Linear factors with a ≠ 0 are convolved /
long-divided; ε-proportional factors shift m and scale by b."""
function _b2pref_rational_series(p::B2Prefactor{EpsExp}, K::Int)
    K >= 0 || throw(ArgumentError("_b2pref_rational_series: K=$K < 0"))
    m = 0
    r = zeros(Rational{BigInt}, K + 1)
    r[1] = p.c
    for x in p.num
        iszero(x) && throw(ArgumentError(
            "b2 prefactor: zero numerator factor — vanished term should not have been emitted"))
        if iszero(x.a)
            m += 1
            r = r .* x.b
        else
            prev = copy(r)
            r = r .* x.a
            for k in 2:K+1
                r[k] += x.b * prev[k-1]
            end
        end
    end
    for x in p.den
        iszero(x) && throw(ArgumentError("b2 prefactor: zero denominator factor (1/0)"))
        if iszero(x.a)
            m -= 1
            r = r .// x.b
        else
            # long division by (a + b·ε): d[k] = (c[k] − b·d[k−1]) / a
            r[1] = r[1] // x.a
            for k in 2:K+1
                r[k] = (r[k] - x.b * r[k-1]) // x.a
            end
        end
    end
    return (m, r)
end

"""Multiply a term's LaurentSeries by its exact B2Prefactor, truncated at
`requested_order`. Returns (minorder, coeffs::Vector{HlogExpr}); the inner
series must reach ε^(requested_order − m) — guaranteed by the caller's
`subord` bookkeeping, guarded here (truncated-Laurent pole guard, R2)."""
function _b2_apply_pref(L::LaurentSeries{HlogExpr}, p::B2Prefactor{EpsExp},
                        requested_order::Int)
    m = _b2pref_poleshift(p)
    outmin = L.minorder + m
    outmin > requested_order && return (outmin, HlogExpr[])  # provably above cutoff
    (m2, r) = _b2pref_rational_series(p, requested_order - outmin)
    @assert m2 == m
    Ltop = L.minorder + length(L.coeffs) - 1
    Ltop >= requested_order - m || throw(LaurentTruncationError(
        "b2 driver: prefactor convolution to ε^$(requested_order) needs the " *
        "term series up to ε^$(requested_order - m) but it stops at ε^$(Ltop)"))
    coeffs = HlogExpr[]
    for o in outmin:requested_order
        acc = HlogExpr()
        for k in 0:(o - outmin)
            j = o - m - k - L.minorder + 1          # index into L.coeffs
            (1 <= j <= length(L.coeffs)) || continue
            iszero(r[k+1]) && continue
            acc = hlog_add(acc, hlog_scale(L.coeffs[j], r[k+1]))
        end
        push!(coeffs, acc)
    end
    return (outmin, coeffs)
end

_b2_eestr(e::EpsExp) = "$(_b1_qstr(e.a)) + ($(_b1_qstr(e.b)))*eps"
# JSON-safe copy of a continuation step ledger (EpsExp/Eps2 → strings)
function _b2_step_ledger(steps)
    out = Any[]
    for st in steps
        d = Dict{String,Any}()
        for (k, v) in st
            d[k] = v isa EpsExp ? _b2_eestr(v) : v
        end
        push!(out, d)
    end
    return out
end

"""
    subtropica_integrate_b2(E, refusal, requested_order; ...) -> LaurentSeries

The B2 continuation driver (header block above): census geometry rebuilt,
`expand_integral_b2` in the dlog chart, each continued output integrated
through `subtropica_integrate` (allow_continuation=false; its own census
routes it to the B1 subtraction path or Phase A), exact B2Prefactor
convolution per term, pole-aligned cross-term sum. Provenance carries the
full continuation ledger (trigger refusal, np_search, per-term steps +
prefactors + sub-evidence)."""
function subtropica_integrate_b2(E::EulerIntegrand, refusal, requested_order::Int;
                              hf::HFConfig=HFConfig(),
                              run_dir::AbstractString=default_run_dir(),
                              domains=Dict{Symbol,Any}(), mzv::Bool=true,
                              strict_periods::Bool=true, on_period=nothing,
                              findroots::Symbol=:cascade, workers::Int=1,
                              dps::Int=60, continuation_order::Int=16)
    t0 = time()
    refusal isa B1RefusalException || throw(ArgumentError(
        "subtropica_integrate_b2: `refusal` must be the census refusal that triggered the route"))
    cps = _b1_census_polys(E)
    isempty(cps) && throw(ArgumentError(
        "subtropica_integrate_b2: no census polys — a census refusal cannot have come from here"))
    Ed = _b1_dlog(E)
    td = tropical_data(cps, copy(E.vars))
    outs = expand_integral_b2(Ed, td; order=continuation_order)

    pieces  = Tuple{Int,Vector{HlogExpr}}[]
    tables  = Vector{AlgLetter}[]
    termev  = Any[]
    skipped = Any[]
    for (i, ci) in enumerate(outs)
        m = _b2pref_poleshift(ci.pref)
        subord = max(0, requested_order - m)
        Efl = _b1_undlog(ci.integrand)
        L = subtropica_integrate(Efl; order=subord, dps=dps, hf=hf,
                              run_dir=joinpath(run_dir, "b2_term_$(i)"),
                              domains=domains, mzv=mzv,
                              strict_periods=strict_periods,
                              on_period=on_period, findroots=findroots,
                              workers=workers)   # allow_continuation=false: ONE pass
        (mo, cs) = _b2_apply_pref(L, ci.pref, requested_order)
        tev = Dict{String,Any}(
            "term" => i,
            "steps" => _b2_step_ledger(get(ci.provenance, "steps", Any[])),
            "prefactor" => Dict{String,Any}(
                "c" => _b1_qstr(ci.pref.c),
                "num" => [_b2_eestr(x) for x in ci.pref.num],
                "den" => [_b2_eestr(x) for x in ci.pref.den]),
            "pole_shift" => m,
            "sub_route" => get(L.evidence, "b1_path", false) ? "b1" : "phase_a",
            "sub_census" => get(L.evidence, "b1_census", ""),
            "sub_lr_order" => get(L.evidence, "lr_order", ""),
        )
        haskey(L.evidence, "b1_faces") &&
            (tev["sub_faces"] = L.evidence["b1_faces"])
        if isempty(cs)                       # term starts above requested_order
            tev["dropped"] = "leading order ε^$(mo) > requested ε^$(requested_order)"
            push!(skipped, tev)
            continue
        end
        push!(pieces, (mo, cs))
        push!(tables, L.letters)
        push!(termev, tev)
    end
    isempty(pieces) && throw(LaurentTruncationError(
        "subtropica_integrate_b2: every continued term starts above ε^$(requested_order)"))

    # letter union + remap, then pole-aligned fold (assemble.jl general fold)
    (letters, remaps) = union_letters(tables)
    pieces = [(mo, HlogExpr[remap_letters(c, remaps[i]) for c in cs])
              for (i, (mo, cs)) in enumerate(pieces)]
    (gmin, padded) = align_laurent(pieces)
    nC = length(padded[1])
    coeffs = [reduce(hlog_add, (p[k] for p in padded); init=HlogExpr())
              for k in 1:nC]

    ev = Dict{String,Any}(
        "check_divergences" => true,          # hard-ON on every inner call
        "lr_order" => Dict{String,Any}(string(t["term"]) => t["sub_lr_order"]
                                       for t in termev),
        "b2_continuation" => Dict{String,Any}(
            "trigger_type" => string(nameof(typeof(refusal))),
            "trigger" => sprint(showerror, refusal),
            "np_search" => get(outs[1].provenance, "np_search",
                               Dict{String,Any}()),
            "continuation_order_cap" => continuation_order,
            "n_continued_terms" => length(outs),
            "terms" => termev,
            "terms_dropped" => skipped),
        "b1_path" => true,                    # divergent route (B2 flavor)
        "run_dir" => String(run_dir),
        "requested_order" => requested_order,
        "b2_workers" => workers,
        "wall_s_total" => time() - t0,
    )
    return LaurentSeries{HlogExpr}(gmin, coeffs, letters, :analytic, ev)
end

# ===========================================================================
# TWO-REGULATOR DRIVER ROUTE (DESIGN_B1 §4c
# successor). subtropica_integrate accepts an Euler2Integrand END-TO-END:
#
#   two-reg census (trop2_on_ray, EXACT Eps2 triples) ──> triggers:
#     * unregulated-by-eps rays (a=0, b=0, c≠0)  — MUST be continued away
#       (the [q^s] extraction evaluates continued integrands at q=0, where
#       such rays are genuinely divergent);
#     * power rays (a>0) — a NP steps (a+1 when in the removal set);
#     * GP violation on the surviving divergent fan.
#   Continuation = the FROZEN B2 core (continue_ray/continue_rays over
#   Eps2, continue.jl), removal search = find_np_continuation_forced
#   (S ⊇ unregulated set).
#
#   [q^s] EXTRACTION — the eq 4.11 gate's MG7 decision transliterated
#   (VERDICT_EQ411.md §3, "exact [q¹] assembly"): the regulator q is KEPT
#   EXACT through the continuation; the target coefficient is an exact
#   convolution of three q-series, never a q → c·ε promotion run:
#     (i)  continuation prefactor B2Prefactor{Eps2} → exact q-Laurent whose
#          coefficients are sums of B2Prefactor{EpsExp} (`_b2q_pref_qseries`
#          — the 1/(−q) continuation poles give the q-pole; mixed factors
#          (a+bε)+cq expand as geometric series in EpsExp prefactors);
#     (ii) face Vol factors vol_det/∏(−TropI(ρ)) with TropI = Bε + Cq →
#          exact rational series Σ_j h_j·q^j·ε^{−d−j} (`h` convolution);
#     (iii)integrand q-dependence — POST-CONTINUATION it must sit in the
#          MONOMIAL exponents only (x^{cq}; q-carrying POLY exponents
#          refuse typed :AuxRegulatorInPolyExponent — v1 scope, exactly
#          the eq 4.10/4.11 class where "the remaining q-dependence is
#          only the x1^{−q} measure factor", FIXTURE_EQ411 §3) — expanded
#          as exact log insertions [q^k]e^{qL} = L^k/k!, L = Σ c_i log x_i
#          (+ Σ c_j Log[P_j] at the COUNTERTERM level, where jacobian
#          exponents inherit q through TropI).
#
#   The q-carrying counterterm data (mono/vol_trops/jacobian exponents as
#   Eps2) is reconstructed EXACTLY from two single-ε embeddings q = 0 and
#   q = c0·ε of the SAME geometry (`_b2q_reconstruct`): subtraction_terms'
#   exponent bookkeeping is affine in the input exponents while the
#   geometry (faces, signs, dets, restricted polys, u's) is exponent-free,
#   so c-parts = (b-parts(c0-embed) − b-parts(0-embed))/c0, with FULL
#   structural asserts between the two runs. This consumes the FROZEN
#   subtract.jl/CounterTerm surfaces untouched (types are frozen — no Eps2
#   CounterTerm exists). The sanctioned q=0-projection precedent is the
#   eq 4.11 gate itself (gate_struct.jl "single-eps projection at q=0").
#
#   The result is the ε-Laurent of  [q^aux_order]( Prefactor(ε) · J(q,ε) )
#   — J the raw two-regulator integral of the quadruple. External
#   q-prefactors that are NOT ε-linear (the eq 4.11 sin(πq)/π) cannot ride
#   in the Prefactor type; callers reduce them exactly first (for eq 4.11:
#   sin(πq)/π = q·(1+O(q²)) and J's simple q-pole give
#   [q¹]Fq = [q⁰]J·Γ-prefactor EXACTLY — VERDICT_EQ411 §3).
# ===========================================================================

"""ν_dlog = ν_flat + 1 view of the two-regulator quadruple (measure bridge,
Eps2 flavor of `_b1_dlog`)."""
_b2q_dlog(E::Euler2Integrand) =
    Euler2Integrand(E.prefactor, Eps2[ν + 1 for ν in E.nu], E.polys,
                    E.vars, E.kinvars, E.ring, E.regulators,
                    copy(E.provenance))

_b2q_eestr(e::Eps2) =
    "$(_b1_qstr(e.a)) + ($(_b1_qstr(e.b)))*eps + ($(_b1_qstr(e.c)))*q"

"""Single-ε FLAT projection of a (dlog-chart) two-reg integrand at q = 0,
plus the monomial q-coefficient vector. Poly exponents MUST be q-free
(:AuxRegulatorInPolyExponent otherwise — v1 scope, header block)."""
function _b2q_q0_flat(E::Euler2Integrand)
    for (P, e) in E.polys
        iszero(e.c) || throw(B2Refusal(:AuxRegulatorInPolyExponent,
            Dict{String,Any}("poly" => string(P), "exponent" => _b2q_eestr(e),
                "hint" => "two-reg v1 scope: post-continuation q-dependence must sit in the monomial exponents only (eq 4.10/4.11 class); rewrite the quadruple or extend the route")))
    end
    E0 = EulerIntegrand(E.prefactor,
                        EpsExp[EpsExp(ν.a - 1, ν.b) for ν in E.nu],
                        Tuple{QQMPolyRingElem,EpsExp}[(P, EpsExp(e.a, e.b))
                                                      for (P, e) in E.polys],
                        E.vars, E.kinvars, E.ring)
    cvec = Rational{BigInt}[ν.c for ν in E.nu]
    return (E0, cvec)
end

# ---------------------------------------------------------------------------
# (i) exact q-Laurent of a B2Prefactor{Eps2}
# ---------------------------------------------------------------------------
"""
    _b2q_pref_qseries(p::B2Prefactor{Eps2}, qtop::Int) -> (qmin, coeffs)

Exact q-Laurent  p(ε,q) = Σ_{j≥0} (Σ coeffs[j+1]) · q^(qmin+j) + O(q^(qtop+1)),
each coefficient a SUM of `B2Prefactor{EpsExp}` terms (kept exact — the
ε-structure is never expanded here). Pure-q factors (a=b=0) shift qmin;
mixed den factors (a+bε)+cq expand as Σ_j (−c)^j q^j/(a+bε)^{j+1}."""
function _b2q_pref_qseries(p::B2Prefactor{Eps2}, qtop::Int)
    qmin = 0
    for x in p.num
        iszero(x) && throw(ArgumentError("b2q prefactor: zero numerator factor"))
        (iszero(x.a) && iszero(x.b)) && (qmin += 1)
    end
    for x in p.den
        iszero(x) && throw(ArgumentError("b2q prefactor: zero denominator factor (1/0)"))
        (iszero(x.a) && iszero(x.b)) && (qmin -= 1)
    end
    L = qtop - qmin + 1
    L <= 0 && return (qmin, Vector{B2Prefactor{EpsExp}}[])
    coeffs = [B2Prefactor{EpsExp}[] for _ in 1:L]
    coeffs[1] = [B2Prefactor{EpsExp}(p.c, EpsExp[], EpsExp[])]
    for x in p.num
        if iszero(x.a) && iszero(x.b)                       # pure q: c·q
            coeffs = [B2Prefactor{EpsExp}[B2Prefactor{EpsExp}(t.c * x.c, t.num, t.den)
                                          for t in cc] for cc in coeffs]
        elseif iszero(x.c)                                  # ε-only factor
            f = B2Prefactor{EpsExp}(1 // 1, EpsExp[EpsExp(x.a, x.b)], EpsExp[])
            coeffs = [B2Prefactor{EpsExp}[b2pref_mul(t, f) for t in cc]
                      for cc in coeffs]
        else                                                # (a+bε) + c·q
            f0 = B2Prefactor{EpsExp}(1 // 1, EpsExp[EpsExp(x.a, x.b)], EpsExp[])
            nc = [B2Prefactor{EpsExp}[] for _ in 1:L]
            for j in 1:L
                append!(nc[j], B2Prefactor{EpsExp}[b2pref_mul(t, f0)
                                                   for t in coeffs[j]])
                j >= 2 && append!(nc[j],
                    B2Prefactor{EpsExp}[B2Prefactor{EpsExp}(t.c * x.c, t.num, t.den)
                                        for t in coeffs[j-1]])
            end
            coeffs = nc
        end
    end
    for x in p.den
        if iszero(x.a) && iszero(x.b)                       # pure q pole: 1/(c·q)
            coeffs = [B2Prefactor{EpsExp}[B2Prefactor{EpsExp}(t.c // x.c, t.num, t.den)
                                          for t in cc] for cc in coeffs]
        elseif iszero(x.c)
            f = B2Prefactor{EpsExp}(1 // 1, EpsExp[], EpsExp[EpsExp(x.a, x.b)])
            coeffs = [B2Prefactor{EpsExp}[b2pref_mul(t, f) for t in cc]
                      for cc in coeffs]
        else
            nc = [B2Prefactor{EpsExp}[] for _ in 1:L]
            for j0 in 1:L, j in 0:(L-j0)
                g = B2Prefactor{EpsExp}((-x.c)^j // 1, EpsExp[],
                                        EpsExp[EpsExp(x.a, x.b) for _ in 1:(j+1)])
                append!(nc[j0+j], B2Prefactor{EpsExp}[b2pref_mul(t, g)
                                                      for t in coeffs[j0]])
            end
            coeffs = nc
        end
    end
    return (qmin, coeffs)
end

# ---------------------------------------------------------------------------
# (iii) exact log-insertion expansion  [q^v] e^{q·L} = L^v / v!
# ---------------------------------------------------------------------------
"""
    _b2q_log_power_terms(items, v) -> Vector{(coef, logs)}

Multinomial expansion of (Σ_t c_t·Log[A_t])^v / v! over `items` =
Vector{(c_t::Rational, A_t::QQMPolyRingElem)} (zero c's pre-filtered by the
caller): each output term = exact rational coefficient ∏ c^α/α! and the log
list [(A_t, α_t)]. CONSTANT arguments: A ≡ 1 kills the term (log 1 = 0);
A ≡ const ≤ 0 refuses LOUDLY (Euclidean violation); positive constants stay
as Log[const] factors (HF grammar-legal; item-3 fold canonicalizes the
output atoms). v = 0 ⇒ the single identity term."""
function _b2q_log_power_terms(items::Vector{Tuple{Rational{BigInt},QQMPolyRingElem}},
                              v::Int)
    v >= 0 || throw(ArgumentError("_b2q_log_power_terms: v=$v < 0"))
    v == 0 && return Tuple{Rational{BigInt},Vector{Pair{QQMPolyRingElem,Int}}}[
        (1 // 1, Pair{QQMPolyRingElem,Int}[])]
    isempty(items) &&
        return Tuple{Rational{BigInt},Vector{Pair{QQMPolyRingElem,Int}}}[]
    out = Tuple{Rational{BigInt},Vector{Pair{QQMPolyRingElem,Int}}}[]
    for alpha in _compositions(v, length(items))
        q = Rational{BigInt}(1)
        logs = Pair{QQMPolyRingElem,Int}[]
        dead = false
        for (t, a) in zip(items, alpha)
            a == 0 && continue
            (c, A) = t
            if is_constant(A)
                cA = leading_coefficient(A)
                isone(cA) && (dead = true; break)          # log 1 = 0
                cA > 0 || throw(EpsExpandRefusal(
                    "b2q insertion: Log[$(cA)] with non-positive constant argument — Euclidean violation"))
            end
            q *= c^a // factorial(big(a))
            push!(logs, A => a)
        end
        dead && continue
        iszero(q) || push!(out, (q, logs))
    end
    return out
end

# ---------------------------------------------------------------------------
# double-embed reconstruction of the q-carrying counterterm data
# ---------------------------------------------------------------------------
const _B2Q_C0_CANDIDATES = (1 // big(3), 1 // big(5), 2 // big(5), 1 // big(7),
                            3 // big(7), 2 // big(9))

"""
    _b2q_reconstruct(E0, cvec, census) -> NamedTuple

Reconstruct the q-parts of the counterterm bookkeeping (header block): run
`subtraction_terms` at the q = 0 embedding (the `census` already computed)
AND at q = c0·ε (c0 from `_B2Q_C0_CANDIDATES`, chosen so no divergent-log
ray degenerates: b + tc·c0 ≠ 0 whenever a = 0), then recover every c-part
by exact linear algebra. Returns
    (cts = cts0, mono_c, polyexp_c, voltrop_c, c0)
with mono_c[i]/polyexp_c[i]/voltrop_c[i] the per-ct q-coefficient vectors
aligned to cts0[i].mono/.polys/.vol_trops. FULL structural asserts between
the two embeds (faces, subfaces, signs, dets, J, restricted polys,
exponent a-parts) — any mismatch is a LOUD internal error, never silence."""
function _b2q_reconstruct(E0::EulerIntegrand, cvec::Vector{Rational{BigInt}},
                          census)
    td = census.td
    trv0 = census.dd.trvals
    tcs = [sum(cvec[i] * Rational{BigInt}(rho[i]) for i in eachindex(cvec))
           for rho in td.rays]
    c0 = nothing
    for cand in _B2Q_C0_CANDIDATES
        ok = all(!iszero(trv0[i].a) || !iszero(trv0[i].b + tcs[i] * cand)
                 for i in eachindex(trv0))
        ok && (c0 = Rational{BigInt}(cand); break)
    end
    c0 === nothing && error("_b2q_reconstruct: no admissible c0 among " *
        "$_B2Q_C0_CANDIDATES — extend the candidate list")
    cts0 = subtraction_terms(census.Ed, td, census.dd.sigma_div, census.dd.w)
    EC = EulerIntegrand(E0.prefactor,
                        EpsExp[EpsExp(E0.nu[i].a, E0.nu[i].b + cvec[i] * c0)
                               for i in eachindex(E0.nu)],
                        E0.polys, E0.vars, E0.kinvars, E0.ring)
    EdC = _b1_dlog(EC)
    ddC = divergence_data(EdC, td)
    ddC.div_facets == census.dd.div_facets || error(
        "_b2q_reconstruct: embed q=$(c0)*eps changed the divergent fan " *
        "($(ddC.div_facets) vs $(census.dd.div_facets)) — c0 constraint bug")
    ddC.w == census.dd.w || error("_b2q_reconstruct: embed changed the w-search output")
    ctsC = subtraction_terms(EdC, td, ddC.sigma_div, ddC.w)
    length(ctsC) == length(cts0) || error(
        "_b2q_reconstruct: counterterm count mismatch between embeds")
    dq(b0::EpsExp, bC::EpsExp) = begin
        b0.a == bC.a || error("_b2q_reconstruct: exponent a-part drifted between embeds")
        (bC.b - b0.b) // c0
    end
    mono_c    = Vector{Vector{Rational{BigInt}}}(undef, length(cts0))
    polyexp_c = Vector{Vector{Rational{BigInt}}}(undef, length(cts0))
    voltrop_c = Vector{Vector{Rational{BigInt}}}(undef, length(cts0))
    for i in eachindex(cts0)
        a = cts0[i]; b = ctsC[i]
        (a.face == b.face && a.subface == b.subface && a.sign == b.sign &&
         a.vol_det == b.vol_det && a.J == b.J && a.face_vars == b.face_vars) ||
            error("_b2q_reconstruct: ct $(i) structural mismatch between embeds")
        length(a.polys) == length(b.polys) || error(
            "_b2q_reconstruct: ct $(i) poly-count mismatch between embeds")
        all(a.polys[j][1] == b.polys[j][1] for j in eachindex(a.polys)) ||
            error("_b2q_reconstruct: ct $(i) restricted-poly mismatch between embeds")
        mono_c[i]    = Rational{BigInt}[dq(a.mono[k], b.mono[k])
                                        for k in eachindex(a.mono)]
        polyexp_c[i] = Rational{BigInt}[dq(a.polys[j][2], b.polys[j][2])
                                        for j in eachindex(a.polys)]
        voltrop_c[i] = Rational{BigInt}[dq(a.vol_trops[k], b.vol_trops[k])
                                        for k in eachindex(a.vol_trops)]
    end
    return (cts = cts0, mono_c = mono_c, polyexp_c = polyexp_c,
            voltrop_c = voltrop_c, c0 = c0)
end

# ---------------------------------------------------------------------------
# inner ε-Laurent of [q^v] of one continued term (single-reg pipeline with
# exact log insertions)
# ---------------------------------------------------------------------------

# Phase-A flavor: insertions on the raw quadruple (no counterterms — the
# census certified local finiteness / could not certify and HF hard-ON is
# the soundness bar, exactly the single-reg Phase-A policy). L = Σ c_i·log x_i
# over the ORIGINAL variables (poly exponents are q-free by _b2q_q0_flat).
function _b2q_phase_a_inner(E0::EulerIntegrand, cvec::Vector{Rational{BigInt}},
                            v::Int, requested_order::Int;
                            hf::HFConfig, run_dir::AbstractString, domains,
                            mzv::Bool, strict_periods::Bool, on_period,
                            findroots::Symbol, census_note::AbstractString)
    t0 = time()
    ex = eps_expand(E0, requested_order; domains=domains)
    items = Tuple{Rational{BigInt},QQMPolyRingElem}[
        (cvec[i], gen(E0.ring, i)) for i in eachindex(cvec) if !iszero(cvec[i])]
    T = _b2q_log_power_terms(items, v)
    letters = String[]
    for (P, _) in E0.polys
        is_constant(P) && continue
        for (f, _) in Nemo.factor(P)
            s = mma_string(f)
            s in letters || push!(letters, s)
        end
    end
    local lr_order::Vector{Symbol}
    lr_meta = Dict{String,Any}()
    alg = false; carry = false
    if isempty(letters)
        lr_order = copy(E0.vars)
        lr_meta["note"] = "no non-constant letters — LR search skipped"
    else
        (lr_order, lr_meta, alg, carry) = _b1_lr_cascade([letters], E0.vars,
            String.(E0.kinvars); run_dir=run_dir, hf=hf, findroots=findroots)
    end
    allvars = vcat(E0.vars, E0.kinvars)
    blocks = Vector{Tuple{HlogExpr,Vector{AlgLetter}}}[]
    hf_wall = 0.0
    for lst in ex.integrands
        entries = Tuple{HlogExpr,Vector{AlgLetter}}[]
        for li in lst, (qc, ilogs) in T
            sli = LogIntegrand(li.num * parent(li.num)(Nemo.QQ(numerator(qc),
                                                              denominator(qc))),
                               li.den, vcat(li.logs, ilogs), li.vars, li.kinvars)
            tI = time()
            (terms, tlet) = hf_integrate(_ti_factored_expr(sli), lr_order,
                                         allvars; cfg=hf,
                                         check_divergences=true,
                                         algebraic_letters=alg,
                                         carry_discharge=carry)
            hf_wall += time() - tI
            val = resolve_periods(terms; cfg=hf, letters=tlet,
                                  strict=strict_periods, on_period=on_period)
            push!(entries, (val, tlet))
        end
        push!(blocks, entries)
    end
    ev = Dict{String,Any}(
        "check_divergences" => true,
        "b1_census" => String(census_note),
        "b2q_insertion_order" => v,
        "b2q_insertion_terms" => length(T),
        "lr_order" => isempty(letters) ? "no_integration_required" : String.(lr_order),
        "lr_meta" => lr_meta,
        "wall_s_hf" => hf_wall,
        "wall_s_total_pre_assemble" => time() - t0,
        "run_dir" => String(run_dir),
        "requested_order" => requested_order,
        "prefactor_minorder" => ex.prefactor_minorder,
    )
    return assemble(ex.minorder, blocks,
                    (ex.prefactor_minorder, ex.prefactor_coeffs),
                    requested_order; evidence=ev, cfg=hf, mzv=mzv)
end

# B1 flavor: the census certified log-divergent rays — Möbius subtraction
# with the q-family carried exactly (header block): per face,
#   [q^v] contribution = Σ_{v1+v2=v}  vol_det·h[v1+1]·ε^{−(d+v1)}
#                        · Σ_ct sign · ∫ [q^{v2}]-inserted ct integrand,
# h = the exact rational q-series of ∏ 1/(−TropI) (TropI = Bε + Cq), and
# the ct insertions L_ct = Σ mono_c·log x_fv + Σ polyexp_c·Log[P^chart]
# (jacobian exponents inherit q through TropI — reconstructed data).
function _b2q_b1_inner(E0::EulerIntegrand, census, cvec::Vector{Rational{BigInt}},
                       v::Int, requested_order::Int;
                       hf::HFConfig, run_dir::AbstractString, domains,
                       mzv::Bool, strict_periods::Bool, on_period,
                       findroots::Symbol, evsink::Dict{String,Any})
    t0 = time()
    rec = _b2q_reconstruct(E0, cvec, census)
    cts = rec.cts
    # group by face (sigma_div order — C7 emission order, as in the B1 driver)
    groups = Vector{Int}[]
    for (i, ct) in enumerate(cts)
        if !isempty(groups) && cts[groups[end][1]].face == ct.face
            push!(groups[end], i)
        else
            push!(groups, Int[i])
        end
    end
    dmax = maximum(length(cts[g[1]].face) for g in groups)
    (pmin0, _) = gamma_prefactor_series(E0.prefactor, 0)
    gmin = -(dmax + v)
    parg = max(0, requested_order - gmin - pmin0)
    (pmin_g, pcoeffs_g) = gamma_prefactor_series(E0.prefactor, parg)

    faces = FaceResult[]
    facemetas = Dict{String,Any}[]
    for (fidx, g) in enumerate(groups)
        ct1 = cts[g[1]]
        d = length(ct1.face)
        fm = Dict{String,Any}("face" => string(ct1.face), "depth" => d,
                              "n_cts" => length(g))
        # ---- face Vol q-series: ∏ 1/(B_r·ε + C_r·q), exact rationals -------
        # (vol_trops = (−TropI) per face ray; A-part must vanish — log rays)
        vt0 = ct1.vol_trops
        vtC = rec.voltrop_c[g[1]]
        for j in g
            (cts[j].vol_trops == vt0 && rec.voltrop_c[j] == vtC) || error(
                "b2q driver: vol_trops differ within face group $(ct1.face)")
        end
        h = zeros(Rational{BigInt}, v + 1); h[1] = 1
        for (r, vt) in enumerate(vt0)
            iszero(vt.a) || throw(AssertionError(
                "b2q driver: vol_trop with nonzero constant part on face $(ct1.face)"))
            iszero(vt.b) && throw(AssertionError(
                "b2q driver: zero eps-part vol_trop on face $(ct1.face)"))
            B = vt.b; C = vtC[r]
            f = Rational{BigInt}[(1 // B) * ((-C) // B)^j for j in 0:v]
            nh = zeros(Rational{BigInt}, v + 1)
            for i in 0:v, j in 0:(v-i)
                nh[i+j+1] += h[i+1] * f[j+1]
            end
            h = nh
        end
        fm["vol_qseries"] = [_b1_qstr(x) for x in h]

        # ---- zero-variable faces: constants only ---------------------------
        if isempty(ct1.face_vars)
            for v2 in 0:v
                v1 = v - v2
                if v2 == 0
                    K = requested_order - pmin_g + d + v1
                    K >= 0 || throw(LaurentTruncationError(
                        "b2q driver: requested_order below face pole bookkeeping"))
                    acc = [HlogExpr() for _ in 0:K]
                    for j in g
                        os = _b1_const_ct_orders(cts[j],
                            cts[j].vol_det * h[v1+1] * cts[j].sign, K)
                        for k in 1:K+1
                            acc[k] = hlog_add(acc[k], os[k])
                        end
                    end
                    orders = [isempty(a.terms) ? Tuple{HlogExpr,Vector{AlgLetter}}[] :
                              Tuple{HlogExpr,Vector{AlgLetter}}[(a, AlgLetter[])]
                              for a in acc]
                    push!(faces, FaceResult(-(d + v1), orders))
                else
                    all(all(iszero, rec.polyexp_c[j]) for j in g) && continue
                    throw(B2Refusal(:AuxZeroVarFaceUnsupported, Dict{String,Any}(
                        "face" => ct1.face,
                        "hint" => "zero-variable face with q-carrying constant exponents at insertion order ≥ 1 — outside the v1 two-reg scope")))
                end
            end
            fm["lr_order"] = "no_integration_required"
            push!(facemetas, fm)
            continue
        end

        # ---- charts, expansions, insertion items ---------------------------
        Kmax = requested_order - pmin_g + d + v
        Kmax >= 0 || throw(LaurentTruncationError(
            "b2q driver: requested_order below face pole bookkeeping (Kmax<0)"))
        quads = [_b1_ct_quad(E0, cts[j]) for j in g]
        exs = [eps_expand(q, Kmax; domains=domains) for q in quads]
        itemss = Vector{Tuple{Rational{BigInt},QQMPolyRingElem}}[]
        for (gi, j) in enumerate(g)
            items = Tuple{Rational{BigInt},QQMPolyRingElem}[]
            for k in eachindex(rec.mono_c[j])
                iszero(rec.mono_c[j][k]) && continue
                push!(items, (rec.mono_c[j][k], gen(quads[gi].ring, k)))
            end
            for pj in eachindex(rec.polyexp_c[j])
                iszero(rec.polyexp_c[j][pj]) && continue
                push!(items, (rec.polyexp_c[j][pj], quads[gi].polys[pj][1]))
            end
            push!(itemss, items)
        end
        # ---- letter groups (LETTER RULE + q-carrying log args) -------------
        lgroups = Vector{String}[]
        for (gi, j) in enumerate(g)
            lg = String[]
            for (pj, (P, e)) in enumerate(quads[gi].polys)
                is_constant(P) && continue
                (iszero(e.b) && e.a >= 0 && iszero(rec.polyexp_c[j][pj])) && continue
                for (f, _) in Nemo.factor(P)
                    s = mma_string(f)
                    s in lg || push!(lg, s)
                end
            end
            isempty(lg) || push!(lgroups, lg)
        end
        allvars = vcat(ct1.face_vars, E0.kinvars)
        local lr_order::Vector{Symbol}
        alg = false; carry = false
        if isempty(lgroups)
            lr_order = copy(ct1.face_vars)
            fm["lr_order"] = "no non-constant letters - LR search skipped"
        else
            (lr_order, lr_meta, alg, carry) = _b1_lr_cascade(lgroups,
                ct1.face_vars, String.(E0.kinvars);
                run_dir=joinpath(run_dir, "b2q_face_$(fidx)"), hf=hf,
                findroots=findroots)
            fm["lr_order"] = String.(lr_order)
            fm["lr_meta"] = lr_meta
        end
        # ---- per (v2, ε-order): ONE hf_integrate on the inserted ct-SUM ----
        hf_wall = 0.0
        for v2 in 0:v
            v1 = v - v2
            K2 = requested_order - pmin_g + d + v1
            Ts = [_b2q_log_power_terms(itemss[gi], v2) for gi in eachindex(g)]
            (v2 > 0 && all(isempty, Ts)) && continue
            orders = Vector{Vector{Tuple{HlogExpr,Vector{AlgLetter}}}}()
            for k in 0:K2
                parts = String[]
                for (gi, j) in enumerate(g)
                    mul0 = Rational{BigInt}(cts[j].vol_det) * h[v1+1] * cts[j].sign
                    for li in exs[gi].integrands[k+1], (qc, ilogs) in Ts[gi]
                        mul = mul0 * qc
                        iszero(mul) && continue
                        sli = LogIntegrand(
                            li.num * parent(li.num)(Nemo.QQ(numerator(mul),
                                                            denominator(mul))),
                            li.den, vcat(li.logs, ilogs), li.vars, li.kinvars)
                        push!(parts, "(" * _ti_factored_expr(sli) * ")")
                    end
                end
                if isempty(parts)
                    push!(orders, Tuple{HlogExpr,Vector{AlgLetter}}[])
                    continue
                end
                expr = join(parts, "+")
                tI = time()
                (terms, tlet) = try
                    hf_integrate(expr, lr_order, allvars; cfg=hf,
                                 check_divergences=true,
                                 algebraic_letters=alg, carry_discharge=carry)
                catch err
                    err isa HFDivergent || rethrow()
                    # same detector false-positive class as the B1 driver
                    # (per-log-atom rational pre-check vs transcendental tail
                    # cancellation); the q-family Möbius cancellation is an
                    # exponent-level identity, so its q-Taylor coefficients
                    # (our inserted sums) inherit local finiteness. STAMPED.
                    fb = get!(fm, "divergence_check_fallbacks", Any[])
                    push!(fb, Dict{String,Any}("order" => k, "v2" => v2,
                          "message" => sprint(showerror, err)))
                    hf_integrate(expr, lr_order, allvars; cfg=hf,
                                 check_divergences=false,
                                 algebraic_letters=alg, carry_discharge=carry)
                end
                hf_wall += time() - tI
                val = resolve_periods(terms; cfg=hf, letters=tlet,
                                      strict=strict_periods, on_period=on_period)
                push!(orders, Tuple{HlogExpr,Vector{AlgLetter}}[(val, tlet)])
            end
            push!(faces, FaceResult(-(d + v1), orders))
        end
        fm["wall_s_hf"] = hf_wall
        push!(facemetas, fm)
    end

    ev = Dict{String,Any}(
        "check_divergences" => true,
        "lr_order" => Dict{String,Any}(fm["face"] => fm["lr_order"]
                                       for fm in facemetas),
        "b1_path" => true,
        "b1_census" => String(census.note),
        "b2q_insertion_order" => v,
        "b2q_embed_c0" => _b1_qstr(rec.c0),
        "b1_faces" => facemetas,
        "b1_n_counterterms" => length(cts),
        "run_dir" => String(run_dir),
        "requested_order" => requested_order,
        "prefactor_minorder" => pmin_g,
        "wall_s_total_pre_assemble" => time() - t0,
    )
    evsink["b2q_embed_c0"] = _b1_qstr(rec.c0)
    return assemble_faces(faces, (pmin_g, pcoeffs_g), requested_order;
                          evidence=ev, cfg=hf, mzv=mzv)
end

"""[q^v] inner dispatch for one continued two-reg term: v = 0 is the plain
single-reg driver (its own census/routing, allow_continuation=false — one
continuation pass total); v ≥ 1 routes on the term's own census to the
insertion-carrying Phase-A or B1 flavor. Typed refusals propagate."""
function _b2q_inner(E0::EulerIntegrand, cvec::Vector{Rational{BigInt}},
                    v::Int, requested_order::Int;
                    hf::HFConfig, run_dir::AbstractString, domains, mzv::Bool,
                    strict_periods::Bool, on_period, findroots::Symbol,
                    workers::Int, dps::Int, evsink::Dict{String,Any})
    if v == 0
        return subtropica_integrate(E0; order=requested_order, dps=dps, hf=hf,
                                 run_dir=run_dir, domains=domains, mzv=mzv,
                                 strict_periods=strict_periods,
                                 on_period=on_period, findroots=findroots,
                                 workers=workers)
    end
    census = b1_ray_census(E0)
    if census.route === :b1
        return _b2q_b1_inner(E0, census, cvec, v, requested_order; hf=hf,
                             run_dir=run_dir, domains=domains, mzv=mzv,
                             strict_periods=strict_periods, on_period=on_period,
                             findroots=findroots, evsink=evsink)
    end
    return _b2q_phase_a_inner(E0, cvec, v, requested_order; hf=hf,
                              run_dir=run_dir, domains=domains, mzv=mzv,
                              strict_periods=strict_periods, on_period=on_period,
                              findroots=findroots, census_note=census.note)
end

# ---------------------------------------------------------------------------
# THE TWO-REGULATOR ENTRY POINT — subtropica_integrate(::Euler2Integrand)
# ---------------------------------------------------------------------------
"""
    subtropica_integrate(E2::Euler2Integrand; aux_order=0, order=:auto,
                      allow_continuation=false, kwargs...)
        -> LaurentSeries{HlogExpr}

END-TO-END two-regulator route (header block of the B2Q
section): the ε-Laurent of  [q^aux_order]( Prefactor(ε) · J(q,ε) )  with
J the raw two-regulator integral of the quadruple and q =
`E2.regulators[2]`. Census via `trop2_on_ray` (exact Eps2 triples);
unregulated-by-eps rays and power rays trigger Nilsson–Passare
continuation THROUGH THE FROZEN B2 CORE (kept exact in both regulators);
the [q^aux_order] coefficient is assembled by EXACT convolution (the
eq 4.11 gate's MG7 decision — never a q → c·ε promotion run).

Default surface law unchanged: with `allow_continuation=false` (default)
inputs needing continuation refuse TYPED
(`B2Refusal(:ContinuationRequired)` naming the flag). External q-dependent
prefactors that are not ε-linear (eq 4.11's sin(πq)/π) cannot ride in
`Prefactor` — reduce them exactly first (docs in the section header).
"""
function subtropica_integrate(E2::Euler2Integrand; order=:auto, dps::Int=60,
                           gauge=:auto, region::Symbol=:euclidean,
                           check_divergences=:auto, findroots::Symbol=:cascade,
                           workers::Int=1, hf::HFConfig=HFConfig(),
                           run_dir::AbstractString=default_run_dir(),
                           domains=Dict{Symbol,Any}(), mzv::Bool=true,
                           strict_periods::Bool=true, on_period=nothing,
                           allow_continuation::Bool=false,
                           continuation_order::Integer=16,
                           aux_order::Integer=0)
    t0 = time()
    region === :euclidean ||
        error("SubTropica.subtropica_integrate: region=$region REFUSED in v1 — Euclidean only (DESIGN.md R4)")
    cd = check_divergences === :auto ? true : Bool(check_divergences)
    cd || throw(DivergenceCheckViolation(
        "subtropica_integrate: check_divergences=false is the Phase-B face-sum " *
        "farm ONLY — standalone calls run the flag HARD-ON"))
    ord = order === :auto ? 2 : Int(order)
    ord >= 0 || error("SubTropica.subtropica_integrate: order=$ord < 0")
    s = Int(aux_order)
    hf_gate(hf) || throw(HFGateError("hf_gate returned false — wrong/mis-stamped binary"))

    Ed2 = _b2q_dlog(E2)
    cps = _b1_census_polys(E2)
    censustab = Any[]
    need_cont = false
    trigger = Dict{String,Any}()
    local td
    if isempty(cps)
        # monomial-only: no polytope, no rays, nothing to continue; HF
        # hard-ON in the inner runs is the soundness bar (the single-reg
        # monomial-only Phase-A policy verbatim)
        td = nothing
    else
        td = tropical_data(cps, copy(E2.vars))
        isempty(td.equations) || throw(B1Refusal(:DegeneratePolytope,
            Dict{String,Any}("equations" => td.equations,
                             "ambient_dim" => td.ambient_dim)))
        trv2 = trop2_values(Ed2, td)
        for (i, t) in enumerate(trv2)
            is_illdefined(t) && throw(B1Refusal(:TropIllDefined, Dict{String,Any}(
                "ray_index" => i, "ray" => td.rays[i],
                "hint" => "two-reg TropI == 0 identically — unfixed projective gauge")))
            push!(censustab, Dict{String,Any}("ray" => string(td.rays[i]),
                                              "tropI" => _b2q_eestr(t)))
        end
        sel   = Int[i for i in eachindex(trv2) if trv2[i].a >= 0]
        unreg = Int[i for i in sel if unregulated_by_eps(trv2[i])]
        power = Int[i for i in sel if trv2[i].a > 0]
        for i in power
            denominator(trv2[i].a) == 1 || throw(B2Refusal(:NonIntegerDivergenceOrder,
                Dict{String,Any}("ray" => i, "tropI" => trv2[i],
                                 "order" => trv2[i].a)))
        end
        gp_broken = false
        if isempty(unreg) && isempty(power) && !isempty(sel)
            try
                produce_ws(td, sigma_div(td, sel))
            catch err
                (err isa GeometricPropertyViolated || err isa WNotUnitNormalized) ||
                    rethrow()
                gp_broken = true
                trigger["gp"] = sprint(showerror, err)
            end
        end
        need_cont = !isempty(unreg) || !isempty(power) || gp_broken
        trigger["unregulated_by_eps_rays"] = unreg
        trigger["power_rays"] = power
        if need_cont && !allow_continuation
            throw(B2Refusal(:ContinuationRequired, Dict{String,Any}(
                "unregulated_by_eps" => unreg, "power" => power,
                "gp" => get(trigger, "gp", ""),
                "hint" => "two-reg census requires Nilsson-Passare continuation — rerun subtropica_integrate with allow_continuation=true")))
        end
        if need_cont
            fr = find_np_continuation_forced(td, sel, unreg)
            np_list = Int[]
            for r in fr.rays
                append!(np_list, fill(r, Int(trv2[r].a) + 1))
            end
            for r in power
                r in fr.rays && continue
                append!(np_list, fill(r, Int(trv2[r].a)))
            end
            length(np_list) > continuation_order && throw(B2Refusal(
                :ContinuationDepthExceeded, Dict{String,Any}(
                    "steps_required" => length(np_list),
                    "order" => Int(continuation_order),
                    "np_list" => copy(np_list))))
            trigger["np_search"] = Dict{String,Any}(
                "rays_removed" => collect(fr.rays),
                "np_ray_sequence" => copy(np_list),
                "subsets_scanned" => fr.scanned)
        end
    end

    outs = if need_cont
        continue_rays(Ed2, [td.rays[r] for r in trigger["np_search"]["np_ray_sequence"]])
    else
        ContinuedIntegrand{Eps2,Euler2Integrand}[
            ContinuedIntegrand{Eps2,Euler2Integrand}(b2pref_identity(Eps2), Ed2,
                Dict{String,Any}("np" => "no continuation needed"))]
    end

    # ---- exact [q^s] convolution over continued terms ----------------------
    pieces  = Tuple{Int,Vector{HlogExpr}}[]
    tables  = Vector{AlgLetter}[]
    termev  = Any[]
    for (i, ci) in enumerate(outs)
        (E0, cvec) = _b2q_q0_flat(ci.integrand)
        (qmin, prefC) = _b2q_pref_qseries(ci.pref, s)
        # per q-order plan: v -> max inner ε-order over the pref terms
        plan = Dict{Int,Int}()
        for (j, cc) in enumerate(prefC)
            isempty(cc) && continue
            v = s - (qmin + j - 1)
            v < 0 && continue
            need = maximum(max(0, ord - _b2pref_poleshift(t)) for t in cc)
            plan[v] = max(get(plan, v, 0), need)
        end
        tev = Dict{String,Any}(
            "term" => i,
            "steps" => _b2_step_ledger(get(ci.provenance, "steps", Any[])),
            "prefactor" => Dict{String,Any}(
                "c" => _b1_qstr(ci.pref.c),
                "num" => [_b2q_eestr(x) for x in ci.pref.num],
                "den" => [_b2q_eestr(x) for x in ci.pref.den]),
            "q_pole_order" => -qmin,
            "insertion_orders" => sort!(collect(keys(plan))),
        )
        inner = Dict{Int,LaurentSeries{HlogExpr}}()
        for (v, subord) in sort!(collect(plan))
            evs = Dict{String,Any}()
            inner[v] = _b2q_inner(E0, cvec, v, subord; hf=hf,
                run_dir=joinpath(run_dir, "b2q_term_$(i)_v$(v)"),
                domains=domains, mzv=mzv, strict_periods=strict_periods,
                on_period=on_period, findroots=findroots, workers=workers,
                dps=dps, evsink=evs)
            isempty(evs) || (tev["inner_v$(v)"] = evs)
        end
        nadded = 0
        for (j, cc) in enumerate(prefC)
            isempty(cc) && continue
            v = s - (qmin + j - 1)
            v < 0 && continue
            L = inner[v]
            for pt in cc
                (mo, cs) = _b2_apply_pref(L, pt, ord)
                isempty(cs) && continue
                push!(pieces, (mo, cs))
                push!(tables, L.letters)
                nadded += 1
            end
        end
        tev["pieces"] = nadded
        push!(termev, tev)
    end
    isempty(pieces) && throw(LaurentTruncationError(
        "subtropica_integrate(two-reg): no contribution at or below ε^$(ord) " *
        "for [q^$(s)] — the requested coefficient may vanish identically or " *
        "start above the requested order"))

    (letters, remaps) = union_letters(tables)
    pieces = [(mo, HlogExpr[remap_letters(c, remaps[i]) for c in cs])
              for (i, (mo, cs)) in enumerate(pieces)]
    (gmin, padded) = align_laurent(pieces)
    nC = length(padded[1])
    coeffs = [reduce(hlog_add, (p[k] for p in padded); init=HlogExpr())
              for k in 1:nC]

    ev = Dict{String,Any}(
        "check_divergences" => true,
        "lr_order" => "b2q: per-inner-run (see b2q.terms)",
        "b1_path" => need_cont,
        "b2q" => Dict{String,Any}(
            "aux_order" => s,
            "regulators" => String.(E2.regulators),
            "census" => censustab,
            "continuation" => need_cont,
            "trigger" => trigger,
            "n_continued_terms" => length(outs),
            "terms" => termev,
            "route" => "exact-q-convolution (MG7 decision, VERDICT_EQ411 §3)"),
        "run_dir" => String(run_dir),
        "requested_order" => ord,
        "wall_s_total" => time() - t0,
    )
    return LaurentSeries{HlogExpr}(gmin, coeffs, letters, :analytic, ev)
end
