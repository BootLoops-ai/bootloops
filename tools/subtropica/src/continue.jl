# continue.jl — C6: Nilsson–Passare tropical continuation (Phase B2 core).
# The B2 companion to the B1 stack (DESIGN_B1 §4 item 5).
# Reference specification for the surrounding machinery:
# SPEC_B1 (see DESIGN_B1.md); behavioral reference
# tools/subtropica/reference/SubTropica.wl (READ-ONLY, cited as wl:NNNN);
# paper citations = the SubTropica paper text (eq 3.22-3.23,
# Sec 3.2.2/3.2.3, Sec 4.3.2 eqs 4.8-4.11).
#
# Include contract: include()d AFTER src/types.jl, src/b1_types.jl,
# src/tropical.jl and src/subtract.jl. NO `using` here (types.jl header
# law); Nemo is in scope from the enclosing module / standalone header.
# This file CONSUMES tropical.jl's surfaces (trop_poly, trop_on_ray,
# trop_values, div_facets, sigma_div, produce_ws, lex_combinations) and
# subtract.jl's surface (subtraction_terms) — no private duplicates
# (DESIGN_B1 §2 ownership; same pattern as subtract.jl's delegation).
#
# MEASURE convention: like tropical.jl/subtract.jl, everything here lives
# in the DLOG chart (SPEC_B1 §0): EulerIntegrand.nu is read as the dlog
# exponent. The flat-measure quadruple boundary is bridged by b1_driver's
# `_b1_dlog` (types.jl nu MEASURE LAW) — B2 wiring into the driver happens
# at B2 integration, not here.
#
# ---------------------------------------------------------------------------
# WHAT C6 IS (paper Sec 3.2.2, eqs 3.22-3.23; wl:10543-10570):
# For a divergent ray rho, the shifted integral I(λ;ε) = ∫ dlog(x) I(x_i λ^{-ρ_i})
# is λ-independent; d/dλ at λ=1 of the STFactor-ed integrand
#     λ^{-TropI(ρ)} · x^ν · ∏_j Q_j(x,λ)^{e_j},
#     Q_j(x,λ) = Σ_m c_m x^m λ^{t_j - m·ρ},  t_j = trop_poly(P_j, ρ)
# yields  0 = -TropI(ρ)·I + Σ_j I_j,  i.e.  I = (1/TropI(ρ)) Σ_j I_j  (eq 3.23)
# with I_j = e_j · x^ν ∏_{k≠j} P_k^{e_k} · P_j^{e_j-1} · D_j,
#      D_j = dQ_j/dλ|_{λ=1} = Σ_{m: m·ρ < t_j} c_m (t_j - m·ρ) x^m.
#
# SIGN PIN (resolves the loose "-1/TropI" phrasing): the .wl builds
# newPref = -1/mons[[1]] (wl:10547) where mons[[1]] is the λ-monomial
# exponent of the factored integrand, mons[[1]] = -TropI(ρ); the NET
# per-term prefactor is therefore  e_j · (+1/TropI(ρ))  — POSITIVE control:
# eq 3.37 continued along (1,1) gives {{3, x1^ε x2^ε (1+x1+x2)^{-1-3ε}}}
# (paper eq 3.38 / paper:1329-1331) = (-3ε)/(-ε) = +3. Verbatim in tests.
#
# Term-splitting parity note: the .wl runs the substituted integrand
# through STFactor (wl:10545), so Mathematica may present a Q_j as a
# PRODUCT of factors and the d/dλ Table then emits one term per factor.
# By the product rule the resulting term LIST is a refinement of ours
# (one term per P_j entry); the SUMS are identical. B1/B2 keep the P_j
# list structure (the polys list is the identity of the integrand here);
# no re-factorization is performed.
#
# STextractCoefficient parity (wl:10554, def wl:10320-10324): upstream
# pulls the variable-free coefficient and the x-monomial content of each
# continued term out of the poly list; here `_extract_content_monomial`
# folds the rational content of D_j into the prefactor and its x-monomial
# gcd (integration-var slots only) into nu. Kinvar-only factors (e.g. a
# bare X11) STAY in the poly list — nu has no kinvar slots; the B1 census
# drops constant factors downstream (wl:10415 / b1_driver _b1_census_polys).
# ---------------------------------------------------------------------------

# ===========================================================================
# §C0. Typed refusals (B2) — same LOUD/typed/never-partial contract as
#      SPEC_B1 §6; subtyped into the B1 exception family so the farm's
#      catch surface is uniform. Kinds:
#        :ContinuationSearchTooLarge — subset-search combinatorics guard
#        :ContinuationDepthExceeded  — total NP-step cap (`order`) exceeded
#        :ContinuationIllDefined     — TropI(ρ) ≡ 0 along the requested ray
#        :NonIntegerDivergenceOrder  — power order a ∉ Z (repeat count ill-def)
#        :TooManyRegulators          — >2 regulators requested (§C5)
# ===========================================================================
struct B2Refusal <: B1RefusalException
    kind   :: Symbol
    detail :: Dict{String,Any}
end
Base.showerror(io::IO, e::B2Refusal) =
    print(io, "subtropica B2: ", e.kind, " — ", e.detail, "; REFUSING")

# ===========================================================================
# §C1. Two-regulator exponent Eps2 = a + b·ε + c·q (paper Sec 4.3.2).
# Defined HERE (b1_types.jl is frozen; EpsExp stays single-ε and this file
# is the sole definer of Eps2 Base methods). Eps2*Eps2 deliberately
# undefined, same ε²-bug-detector rationale as EpsExp (b1_types §preamble).
# ===========================================================================
struct Eps2
    a::Rational{BigInt}   # constant part
    b::Rational{BigInt}   # coefficient of ε (the primary regulator)
    c::Rational{BigInt}   # coefficient of q (the auxiliary regulator)
end
Eps2(a::Union{Integer,Rational}) = Eps2(Rational{BigInt}(a), 0//1, 0//1)
Eps2(x::EpsExp) = Eps2(x.a, x.b, 0//1)                 # single-ε embedding
Base.:(==)(x::Eps2, y::Eps2) = x.a == y.a && x.b == y.b && x.c == y.c
Base.hash(x::Eps2, h::UInt) = hash((x.a, x.b, x.c), hash(:Eps2, h))
Base.zero(::Type{Eps2}) = Eps2(0//1, 0//1, 0//1)
Base.zero(::Eps2)       = zero(Eps2)
Base.one(::Type{Eps2})  = Eps2(1//1, 0//1, 0//1)
Base.one(::Eps2)        = one(Eps2)
Base.iszero(x::Eps2)    = iszero(x.a) && iszero(x.b) && iszero(x.c)
Base.:+(x::Eps2, y::Eps2) = Eps2(x.a + y.a, x.b + y.b, x.c + y.c)
Base.:-(x::Eps2)          = Eps2(-x.a, -x.b, -x.c)
Base.:-(x::Eps2, y::Eps2) = Eps2(x.a - y.a, x.b - y.b, x.c - y.c)
Base.:+(x::Eps2, c::Union{Integer,Rational}) = Eps2(x.a + c, x.b, x.c)
Base.:+(c::Union{Integer,Rational}, x::Eps2) = x + c
Base.:-(x::Eps2, c::Union{Integer,Rational}) = Eps2(x.a - c, x.b, x.c)
Base.:-(c::Union{Integer,Rational}, x::Eps2) = Eps2(c - x.a, -x.b, -x.c)
Base.:*(x::Eps2, c::Union{Integer,Rational}) = Eps2(x.a * c, x.b * c, x.c * c)
Base.:*(c::Union{Integer,Rational}, x::Eps2) = x * c

# Exact evaluation at rational regulator values.
epsexp_eval(t::EpsExp, e::Rational{BigInt}) = t.a + t.b * e
eps2_eval(t::Eps2, e::Rational{BigInt}, q::Rational{BigInt}) =
    t.a + t.b * e + t.c * q

# ===========================================================================
# §C2. B2 prefactor carrier — c · ∏ num_k / ∏ den_k with ε-linear (EpsExp)
# or two-regulator (Eps2) factors, KEPT EXACT (regulators never expanded
# here; the .wl keeps prefactors as raw symbolics through FullSimplify,
# wl:10563 — this is the typed equivalent). The explicit continuation poles
# live in `den` (the TropI(ρ) factors), exactly like B1's vol_trops carry
# the subtraction poles (SPEC_B1 §4 step 7).
# ===========================================================================
struct B2Prefactor{T}
    c   :: Rational{BigInt}
    num :: Vector{T}          # e.g. the e_j exponent factors, one per step
    den :: Vector{T}          # the TropI(ρ) factors, one per step
end
B2Prefactor{T}() where {T} = B2Prefactor{T}(1//1, T[], T[])
b2pref_identity(::Type{T}) where {T} = B2Prefactor{T}()

# fieldwise ==/hash (same pattern/rationale as types.jl/b1_types.jl —
# default struct == is field-identity and fails on equal Vectors)
Base.:(==)(x::B2Prefactor, y::B2Prefactor) =
    x.c == y.c && x.num == y.num && x.den == y.den
Base.hash(x::B2Prefactor, h::UInt) =
    hash((x.c, x.num, x.den), hash(:B2Prefactor, h))
Base.:(==)(x::B2Refusal, y::B2Refusal) =
    x.kind == y.kind && isequal(x.detail, y.detail)
Base.hash(x::B2Refusal, h::UInt) = hash((x.kind, x.detail), hash(:B2Refusal, h))

b2pref_mul(p::B2Prefactor{T}, q::B2Prefactor{T}) where {T} =
    B2Prefactor{T}(p.c * q.c, vcat(p.num, q.num), vcat(p.den, q.den))

# proportionality r with x == r·y (exact), or nothing
function _proportional_ratio(x::EpsExp, y::EpsExp)
    (iszero(x) || iszero(y)) && return nothing     # keep zeros loud, never cancel

    r = !iszero(y.a) ? x.a // y.a : x.b // y.b
    return (x.a == r * y.a && x.b == r * y.b) ? r : nothing
end
function _proportional_ratio(x::Eps2, y::Eps2)
    (iszero(x) || iszero(y)) && return nothing     # keep zeros loud, never cancel

    r = !iszero(y.a) ? x.a // y.a : (!iszero(y.b) ? x.b // y.b : x.c // y.c)
    return (x.a == r * y.a && x.b == r * y.b && x.c == r * y.c) ? r : nothing
end

"""
    b2pref_simplify(p) -> B2Prefactor

Cancel proportional num/den factor pairs into the rational constant, and
fold regulator-free (pure-rational) factors into `c`. Exact; display/parity
convenience only (eq 3.38's prefactor collapses to the literal 3)."""
function b2pref_simplify(p::B2Prefactor{T}) where {T}
    c = p.c
    num = copy(p.num); den = copy(p.den)
    # regulator-free factors -> c
    for v in (num, den)
        keep = T[]
        for x in v
            if x == (x isa EpsExp ? EpsExp(x.a, 0//1) : Eps2(x.a, 0//1, 0//1))
                iszero(x.a) && (push!(keep, x); continue)   # keep 0 (pole/zero: loud)
                c = (v === num) ? c * x.a : c // x.a
            else
                push!(keep, x)
            end
        end
        empty!(v); append!(v, keep)
    end
    # proportional num/den pairs -> c
    i = 1
    while i <= length(num)
        hit = false
        for j in eachindex(den)
            r = _proportional_ratio(num[i], den[j])
            r === nothing && continue
            c *= r
            deleteat!(num, i); deleteat!(den, j)
            hit = true
            break
        end
        hit || (i += 1)
    end
    return B2Prefactor{T}(c, num, den)
end

"""Exact rational value of a single-ε prefactor at ε = e (throws
DivideError on a vanishing denominator factor — loud, never silent)."""
function b2pref_eval(p::B2Prefactor{EpsExp}, e::Rational{BigInt})
    v = p.c
    for x in p.num; v *= epsexp_eval(x, e); end
    for x in p.den; v //= epsexp_eval(x, e); end
    return v
end
function b2pref_eval(p::B2Prefactor{Eps2}, e::Rational{BigInt}, q::Rational{BigInt})
    v = p.c
    for x in p.num; v *= eps2_eval(x, e, q); end
    for x in p.den; v //= eps2_eval(x, e, q); end
    return v
end

# ===========================================================================
# §C3. Continued-integrand carrier — the typed {prefactor, integrand} pair
# of STContinueRays (wl:10559-10570). `provenance` records the step ledger
# (ray, poly index, exponent, TropI, extracted content/monomial) and the
# regulator-promotion record (§C5) — auditability law.
# ===========================================================================
struct ContinuedIntegrand{T,IT}
    pref       :: B2Prefactor{T}
    integrand  :: IT
    provenance :: Dict{String,Any}
end
const B2Continued = ContinuedIntegrand{EpsExp,EulerIntegrand}

# ===========================================================================
# §C4. The single Nilsson–Passare step — STcontinueRay (wl:10543-10556),
# eqs 3.22-3.23.
# ===========================================================================

# D_j = dQ_j/dλ|_{λ=1} = Σ_{m: m·ρ < t_j} c_m (t_j − m·ρ) x^m  — the
# "shaved facet" derivative (paper below eq 3.23). Exact; zero iff P_j
# equals its own initial form along ρ (all support on the exposed face) —
# such terms are dropped (wl:10553 DeleteCases[#,0]&).
function _np_deriv_poly(P::QQMPolyRingElem, rho::AbstractVector, n::Integer,
                        tj::Rational{BigInt})
    R = parent(P)
    B = MPolyBuildCtx(R)
    for i in 1:length(P)
        ev = exponent_vector(P, i)
        md = _xdot(ev, rho, n)                     # tropical.jl exact dot
        md < tj || continue                        # face terms drop
        push_term!(B, coeff(P, i) * Nemo.QQ(tj - md), ev)
    end
    return finish(B)
end

# STextractCoefficient equivalent (wl:10554/10320-10324): pull the rational
# content and the integration-var monomial gcd out of D. Returns
# (content::Rational, mu::Vector{Int} length n, Dtil) with
# D == content · x^mu · Dtil and Dtil either 1 (drop) or a genuine poly
# (possibly kinvar-only — stays in the poly list, see header note).
function _extract_content_monomial(D::QQMPolyRingElem, n::Integer)
    @assert !iszero(D) "_extract_content_monomial: zero polynomial"
    R = parent(D)
    c0 = content(D)                                # positive rational gcd
    ccont = Rational{BigInt}(BigInt(numerator(c0)), BigInt(denominator(c0)))
    D1 = divexact(D, c0)
    mu = [minimum(exponent_vector(D1, i)[k] for i in 1:length(D1)) for k in 1:n]
    if any(!iszero, mu)
        B = MPolyBuildCtx(R)
        push_term!(B, one(Nemo.QQ), vcat(Int.(mu), zeros(Int, nvars(R) - n)))
        D1 = divexact(D1, finish(B))
    end
    if length(D1) == 1 && all(iszero, exponent_vector(D1, 1))
        c1 = coeff(D1, 1)                          # ±1 by construction
        ccont *= Rational{BigInt}(BigInt(numerator(c1)), BigInt(denominator(c1)))
        return (ccont, mu, one(R))                 # pure monomial: drop
    end
    return (ccont, mu, D1)
end

# Generic step core over exponent type T (EpsExp or Eps2): one continued
# integrand per poly entry with nonzero derivative. Returns Vector of
# (j, efactor::T, ccont, nu2, polys2, dtil_kept::Bool).
function _np_step_core(nu::Vector{T}, polys::Vector{Tuple{QQMPolyRingElem,T}},
                       rho::AbstractVector, n::Integer) where {T}
    out = Tuple{Int,T,Rational{BigInt},Vector{T},
                Vector{Tuple{QQMPolyRingElem,T}},Bool}[]
    for (j, (P, e)) in enumerate(polys)
        iszero(e) && continue                      # e_j = 0 ⇒ term ≡ 0 (wl:10553)
        tj = trop_poly(P, rho, n)                  # tropical.jl (consumed)
        D = _np_deriv_poly(P, rho, n, tj)
        iszero(D) && continue                      # wl:10553
        (ccont, mu, Dtil) = _extract_content_monomial(D, n)
        nu2 = copy(nu)
        for i in 1:n
            iszero(mu[i]) || (nu2[i] += mu[i])
        end
        polys2 = copy(polys)
        polys2[j] = (P, e - 1)                     # exponent shift −1 (wl:10551)
        kept = !isone(Dtil)
        kept && push!(polys2, (Dtil, one(T)))      # derivative factor, exponent 1
        push!(out, (j, e, ccont, nu2, polys2, kept))
    end
    return out
end

"""
    continue_ray(E::EulerIntegrand, td::TropicalData, ray::Integer) -> Vector{B2Continued}
    continue_ray(E::EulerIntegrand, rho::AbstractVector{<:Integer})  -> Vector{B2Continued}

ONE Nilsson–Passare continuation step along `ray` (global index into
`td.rays`) or an explicit primitive integer ray vector — STcontinueRay
(wl:10543-10556), paper eqs 3.22-3.23. Emits one continued integrand per
poly entry with a nonzero λ-derivative, each carrying the exact prefactor
e_j / TropI(ρ) (see the §header SIGN PIN: net = e_j·(+1/TropI); positive
control eq 3.38 pins +3). Regulators are KEPT exact (EpsExp arithmetic
throughout); the input's `Prefactor` passes through UNTOUCHED (C8 global
apply — same law as subtract.jl).

Internal invariant (paper's "shaving" statement below eq 3.23, made exact):
every emitted integrand has TropI_new(ρ) = TropI(ρ) − k with k ∈ Z, k ≥ 1
(a-part strictly reduced, b-part unchanged) — asserted per term.

Refuses B2Refusal(:ContinuationIllDefined) when TropI(ρ) ≡ 0 (the eq 3.23
prefactor 1/TropI(ρ) is 1/0; upstream would emit a raw ComplexInfinity).
The caller iterates via `continue_rays`/`expand_integral_b2`.
"""
function continue_ray(E::EulerIntegrand, rho::AbstractVector{<:Integer};
                      _ray_index::Union{Int,Nothing}=nothing)
    n = length(E.vars)
    length(rho) == n || throw(ArgumentError(
        "continue_ray: ray length $(length(rho)) != n = $n"))
    T = trop_on_ray(E, rho)                        # tropical.jl (consumed)
    is_illdefined(T) && throw(B2Refusal(:ContinuationIllDefined,
        Dict{String,Any}("ray" => collect(rho), "ray_index" => _ray_index,
                         "tropI" => T,
                         "hint" => "TropI(rho) == 0: eq 3.23 prefactor 1/TropI undefined — unfixed gauge or wrong ray")))
    out = B2Continued[]
    for (j, e, ccont, nu2, polys2, kept) in _np_step_core(E.nu, E.polys, rho, n)
        E2 = EulerIntegrand(E.prefactor, nu2, polys2, E.vars, E.kinvars, E.ring)
        # shaving invariant: ΔTropI(ρ) = trop_D(ρ) − t_j ∈ Z, ≤ −1
        Tn = trop_on_ray(E2, rho)
        (Tn.b == T.b && Tn.a <= T.a - 1 && denominator(Tn.a - T.a) == 1) ||
            error("continue_ray: shaving invariant violated (TropI ",
                  "$(T.a)+($(T.b))eps -> $(Tn.a)+($(Tn.b))eps on rho=$(collect(rho)))")
        prov = Dict{String,Any}("steps" => Any[Dict{String,Any}(
            "ray" => collect(rho), "ray_index" => _ray_index, "poly" => j,
            "efactor" => e, "tropI" => T, "content" => ccont,
            "deriv_poly_kept" => kept)])
        push!(out, B2Continued(
            B2Prefactor{EpsExp}(ccont, EpsExp[e], EpsExp[T]), E2, prov))
    end
    return out
end
function continue_ray(E::EulerIntegrand, td::TropicalData, ray::Integer)
    1 <= ray <= length(td.rays) || throw(ArgumentError(
        "continue_ray: ray index $ray out of range 1:$(length(td.rays))"))
    td.ambient_dim == length(E.vars) || throw(ArgumentError(
        "continue_ray: TropicalData ambient_dim != length(vars)"))
    return continue_ray(E, td.rays[ray]; _ray_index = Int(ray))
end

"""
    continue_rays(state::Vector{B2Continued}, E-ray sequence) -> Vector{B2Continued}

STContinueRays (wl:10559-10570): apply `continue_ray` sequentially along
`rays_seq` (a Vector of ray vectors), flat-joining the per-integrand term
lists and MULTIPLYING prefactors exactly (the .wl's
FullSimplify[pri[[1]]·#[[1]]] becomes exact B2Prefactor multiplication).
Empty ray list is the identity (wl:10561). Step ledgers concatenate."""
function continue_rays(state::Vector{B2Continued},
                       rays_seq::AbstractVector{<:AbstractVector})
    for rho in rays_seq
        nxt = B2Continued[]
        for ci in state
            for t in continue_ray(ci.integrand, rho)
                prov = Dict{String,Any}(
                    "steps" => vcat(get(ci.provenance, "steps", Any[]),
                                    t.provenance["steps"]))
                for (k, v) in ci.provenance
                    k == "steps" || (prov[k] = v)
                end
                push!(nxt, B2Continued(b2pref_mul(ci.pref, t.pref),
                                       t.integrand, prov))
            end
        end
        state = nxt
    end
    return state
end
continue_rays(E::EulerIntegrand, rays_seq::AbstractVector{<:AbstractVector}) =
    continue_rays(B2Continued[B2Continued(b2pref_identity(EpsExp), E,
                                          Dict{String,Any}())], rays_seq)

# ===========================================================================
# §C5. Multi-regulator support (paper Sec 4.3.2, eq 4.10-4.11 treatment).
# Exactly TWO regulators (ε primary, one auxiliary — the paper's q);
# anything more is a typed refusal (:TooManyRegulators). The pattern:
#  1. build the two-regulator integrand (Euler2Integrand, Eps2 exponents);
#  2. the naive ε-only census is BLIND to q-regulated rays (paper: "the
#     tropical analysis only sees the ε-regulated divergence");
#  3. PROMOTE q -> q·ε (paper's `integrand /. q -> q eps`): on Eps2 data
#     (a,b,c) [= a + bε + cq] this is a REINTERPRETATION to a + (b + c·q)ε
#     with q a generic parameter — same stored triple, promoted
#     classification: regulated ⇔ (b,c) ≠ (0,0). Recorded in provenance.
#  4. continue along the q-only rays (continue_ray method below), THEN
#     extract the O(q) series coefficient (downstream job, not B2 core —
#     the paper does it with SeriesCoefficient after the continuation).
# ===========================================================================
struct Euler2Integrand
    prefactor  :: Prefactor
    nu         :: Vector{Eps2}
    polys      :: Vector{Tuple{QQMPolyRingElem,Eps2}}
    vars       :: Vector{Symbol}
    kinvars    :: Vector{Symbol}
    ring       :: QQMPolyRing
    regulators :: Vector{Symbol}          # exactly [eps_like, aux]
    provenance :: Dict{String,Any}

    function Euler2Integrand(prefactor, nu, polys, vars, kinvars, ring,
                             regulators, provenance = Dict{String,Any}())
        length(regulators) > 2 && throw(B2Refusal(:TooManyRegulators,
            Dict{String,Any}("regulators" => collect(regulators),
                             "count" => length(regulators),
                             "hint" => "B2 supports exactly 2 regulators (eps + one auxiliary, paper Sec 4.3.2); reduce or split the problem")))
        length(regulators) == 2 || throw(ArgumentError(
            "Euler2Integrand: exactly 2 regulators required (got $(length(regulators))); single-regulator inputs use EulerIntegrand"))
        length(nu) == length(vars) || throw(ArgumentError(
            "Euler2Integrand: length(nu) != length(vars)"))
        new(prefactor, nu, polys, vars, kinvars, ring,
            collect(regulators), provenance)
    end
end

"""TropI(ρ) for the two-regulator integrand — same eq 3.12 sum, Eps2-exact."""
function trop2_on_ray(E::Euler2Integrand, rho::AbstractVector)::Eps2
    n = length(E.vars)
    length(rho) == n || throw(ArgumentError(
        "trop2_on_ray: ray length $(length(rho)) != n = $n"))
    acc = zero(Eps2)
    for i in 1:n
        acc += E.nu[i] * Rational{BigInt}(rho[i])
    end
    for (P, e) in E.polys
        acc += e * trop_poly(P, rho, n)            # tropical.jl (consumed)
    end
    return acc
end
trop2_values(E::Euler2Integrand, td::TropicalData) =
    Eps2[trop2_on_ray(E, rho) for rho in td.rays]

# Classification views. UNPROMOTED (the paper's naive first call): the
# ε-only analysis — a ray with a == 0, b == 0, c ≠ 0 is divergent but
# UNREGULATED BY ε (the paper's abort case). PROMOTED (q -> q·ε):
# regulated ⇔ (b,c) ≠ (0,0) at generic q.
is_illdefined(t::Eps2)     = iszero(t.a) && iszero(t.b) && iszero(t.c)
is_divergent(t::Eps2)      = t.a >= 0
is_power(t::Eps2)          = t.a > 0
is_log(t::Eps2)            = iszero(t.a) && !(iszero(t.b) && iszero(t.c))
unregulated_by_eps(t::Eps2) = iszero(t.a) && iszero(t.b) && !iszero(t.c)

"""
    promote_second_regulator(E::Euler2Integrand) -> Euler2Integrand

The paper's q -> q·ε promotion (Sec 4.3.2: `STPreAnalysis[jac integrand /.
q -> q eps, ...]`). On Eps2 data the stored triple (a,b,c) is UNCHANGED —
the promotion reinterprets a + b·ε + c·q as a + (b + c·q)·ε with q a
generic parameter, so every q-regulated ray becomes ε-regulated (paper:
trops {ε, ε−q, −q} -> {ε, −ε(−1+q), −εq}). The promotion is RECORDED in
provenance (`"promotion" => "q -> q*eps"` + the regulator names); promoting
twice refuses (the double substitution q -> qε² has no Eps2 home)."""
function promote_second_regulator(E::Euler2Integrand)
    haskey(E.provenance, "promotion") && throw(B2Refusal(:TooManyRegulators,
        Dict{String,Any}("regulators" => E.regulators,
                         "promotion" => E.provenance["promotion"],
                         "hint" => "already promoted — a second promotion would need an eps^2 exponent slot")))
    prov = copy(E.provenance)
    prov["promotion"] = "$(E.regulators[2]) -> $(E.regulators[2])*$(E.regulators[1])"
    prov["promotion_law"] = "Eps2(a,b,c) reinterpreted: a + (b + c*q)*eps, q generic (paper Sec 4.3.2)"
    return Euler2Integrand(E.prefactor, E.nu, E.polys, E.vars, E.kinvars,
                           E.ring, E.regulators, prov)
end
is_promoted(E::Euler2Integrand) = haskey(E.provenance, "promotion")

"""continue_ray for the two-regulator integrand (the paper's
STTropicalContinuation call on eq 4.10 along {1,1}). Same NP step, Eps2
prefactors; ill-definedness = a ≡ b ≡ c ≡ 0."""
function continue_ray(E::Euler2Integrand, rho::AbstractVector{<:Integer})
    n = length(E.vars)
    T = trop2_on_ray(E, rho)
    is_illdefined(T) && throw(B2Refusal(:ContinuationIllDefined,
        Dict{String,Any}("ray" => collect(rho), "tropI" => T)))
    out = ContinuedIntegrand{Eps2,Euler2Integrand}[]
    for (j, e, ccont, nu2, polys2, kept) in _np_step_core(E.nu, E.polys, rho, n)
        prov = copy(E.provenance)
        prov["steps"] = vcat(get(E.provenance, "steps", Any[]),
            Any[Dict{String,Any}("ray" => collect(rho), "poly" => j,
                                 "efactor" => e, "tropI" => T,
                                 "content" => ccont, "deriv_poly_kept" => kept)])
        E2 = Euler2Integrand(E.prefactor, nu2, polys2, E.vars, E.kinvars,
                             E.ring, E.regulators, prov)
        push!(out, ContinuedIntegrand{Eps2,Euler2Integrand}(
            B2Prefactor{Eps2}(ccont, Eps2[e], Eps2[T]), E2, prov))
    end
    return out
end

# ===========================================================================
# §C6. Minimal-continuation search — STfindFirstNPContinuation
# (wl:10850-10895): smallest subset of divergent rays whose REMOVAL makes
# the w-search succeed with SIMPLE u's.
# ===========================================================================

# The .wl "satisfactory u" filter (wl:10873-10876), transliterated to
# w-space. Upstream rejects a hypus set when
#   Count[hypus, "NotFound" | var^(n>1), ∞] > 0
#     — some u missing (GP) or some variable appears squared+ in
#       u = x^{w⁺}/(x^{w⁻}+x^{w⁺})  ⟺  some |w_i| > 1;
#   Count[Denominator/@hypus, den with > multiDegreeNP var-occurrences] > 0
#     — the denominator x^{w⁻}+x^{w⁺} mentions more than multiDegreeNP
#       (= 4, wl:2670) variables ⟺ count(w_i ≠ 0) > 4 (given |w_i| ≤ 1).
# NOTE this is an NP-SEARCH filter only, NOT a subtraction precondition
# (SPEC_B1 §6 non-refusals note).
_w_is_simple(w::AbstractVector, multi_degree_np::Integer) =
    all(abs(x) <= 1 for x in w) && count(!iszero, w) <= multi_degree_np

"""
    find_first_np_continuation(td, div_facets; kwargs...) -> NamedTuple

STfindFirstNPContinuation (wl:10850-10895): scan subsets S of the divergent
rays in Mathematica Subsets order (size ascending from `start_from`, lex
within each size — `lex_combinations`, tropical.jl); for each S, pretend
the rays in S have been deleted by NP continuation and test whether the
remaining divergent fan admits the geometric property with SIMPLE
u-functions (`_w_is_simple`; upstream wl:10873-10876, multiDegreeNP=4 at
wl:2670). First success wins. Returns
    (rays = S, ws, sigma, scanned)
with `ws`/`sigma` the successful `produce_ws`/`sigma_div` output for the
SURVIVING divergent set (upstream's hypus, normalized — here produce_ws's
D2 pin `w·ρ = −1` stands in for normalizeQ=True; a WNotUnitNormalized on a
candidate subset counts as subset FAILURE, like GP, since upstream's
rational normalization would produce fractional exponents downstream).

Termination: S = all divergent rays always succeeds vacuously (empty fan),
so the search cannot exhaust — upstream's unreachable {NotFound, NotFound}
tail (wl:10888-10891) is NOT ported. `forcedUs` short-circuit (wl:10857-
10859) is not ported either: B1/B2 have no external-u injection channel.

Guards (typed, BEFORE the exponential scan):
  * length(div_facets) > max_div_facets (default 16)  ⇒
      B2Refusal(:ContinuationSearchTooLarge)
  * more than `subset_budget` subsets scanned (default 100_000) ⇒
      B2Refusal(:ContinuationSearchTooLarge)
sigma_div's own :SigmaDivTooLarge (>24) propagates untouched.

Second-arg forms: a Vector{Int} of divergent ray indices (ascending), or
the NamedTuple returned by `divergence_data` (its .div_facets is used —
only reachable for GP-clean inputs, kept for API symmetry)."""
function find_first_np_continuation(td::TropicalData, div_facets_::Vector{Int};
                                    start_from::Integer = 0,
                                    multi_degree_np::Integer = 4,
                                    max_div_facets::Integer = 16,
                                    subset_budget::Integer = 100_000)
    issorted(div_facets_) || throw(ArgumentError(
        "find_first_np_continuation: div_facets must ascend (SPEC_B1 §3.2 ordering law)"))
    nd = length(div_facets_)
    nd > max_div_facets && throw(B2Refusal(:ContinuationSearchTooLarge,
        Dict{String,Any}("count" => nd, "limit" => Int(max_div_facets),
                         "div_facets" => copy(div_facets_),
                         "hint" => "2^$(nd) candidate subsets — raise max_div_facets deliberately or shrink the problem")))
    scanned = 0
    for k in Int(start_from):nd                    # size ascending (wl loop)
        for S in lex_combinations(div_facets_, k)  # Subsets order within size
            scanned += 1
            scanned > subset_budget && throw(B2Refusal(:ContinuationSearchTooLarge,
                Dict{String,Any}("scanned" => scanned,
                                 "subset_budget" => Int(subset_budget),
                                 "count" => nd)))
            rest = setdiff(div_facets_, S)         # order-preserving: still ascending
            sigma = sigma_div(td, rest)            # tropical.jl (consumed)
            ws = try
                produce_ws(td, sigma)              # tropical.jl (consumed)
            catch err
                (err isa GeometricPropertyViolated ||
                 err isa WNotUnitNormalized) && continue
                rethrow()
            end
            all(_w_is_simple(w, multi_degree_np) for w in values(ws)) || continue
            return (rays = S, ws = ws, sigma = sigma, scanned = scanned)
        end
    end
    error("find_first_np_continuation: unreachable — removing all divergent rays always succeeds")
end
find_first_np_continuation(td::TropicalData, dd::NamedTuple; kwargs...) =
    find_first_np_continuation(td, collect(Int, dd.div_facets); kwargs...)

# ===========================================================================
# §C6b. TWO-REGULATOR DRIVER GLUE (core
# semantics above UNTOUCHED; these are the Eps2 driver-side companions the
# b1_driver two-reg route consumes).
# ===========================================================================

"""
    continue_rays(state::Vector{ContinuedIntegrand{Eps2,Euler2Integrand}}, rays_seq)
    continue_rays(E::Euler2Integrand, rays_seq)

Two-regulator sequential continuation — the Eps2 mirror of the single-ε
`continue_rays` above (STContinueRays, wl:10559-10570): apply the Eps2
`continue_ray` along each ray in `rays_seq`, flat-joining term lists and
multiplying `B2Prefactor{Eps2}`s exactly. Step ledgers thread through the
INTEGRAND provenance (the Eps2 `continue_ray` chains them there), so the
carrier's provenance is taken from the newest step."""
function continue_rays(state::Vector{ContinuedIntegrand{Eps2,Euler2Integrand}},
                       rays_seq::AbstractVector{<:AbstractVector})
    for rho in rays_seq
        nxt = ContinuedIntegrand{Eps2,Euler2Integrand}[]
        for ci in state
            for t in continue_ray(ci.integrand, rho)
                push!(nxt, ContinuedIntegrand{Eps2,Euler2Integrand}(
                    b2pref_mul(ci.pref, t.pref), t.integrand, t.provenance))
            end
        end
        state = nxt
    end
    return state
end
continue_rays(E::Euler2Integrand, rays_seq::AbstractVector{<:AbstractVector}) =
    continue_rays(ContinuedIntegrand{Eps2,Euler2Integrand}[
        ContinuedIntegrand{Eps2,Euler2Integrand}(b2pref_identity(Eps2), E,
                                                 copy(E.provenance))], rays_seq)

"""
    find_np_continuation_forced(td, div_facets, force; kwargs...) -> NamedTuple

Minimal-removal search CONSTRAINED to contain `force` (ascending global ray
indices): scan subsets S ⊇ force of the divergent rays, |S∖force| ascending
then lex (`lex_combinations`), same accept test as
`find_first_np_continuation` (GP + D2 pin + `_w_is_simple`). The two-reg
driver uses `force` = the rays regulated ONLY by the auxiliary regulator
(`unregulated_by_eps`) — those MUST be continued away regardless of the GP
verdict, because the [q^k]-extraction evaluates the continued integrands at
q = 0 where such rays are genuinely divergent. Same typed guards as the
core search. `force ⊆ div_facets` required."""
function find_np_continuation_forced(td::TropicalData, div_facets_::Vector{Int},
                                     force::Vector{Int};
                                     multi_degree_np::Integer = 4,
                                     max_div_facets::Integer = 16,
                                     subset_budget::Integer = 100_000)
    issorted(div_facets_) || throw(ArgumentError(
        "find_np_continuation_forced: div_facets must ascend (SPEC_B1 §3.2 ordering law)"))
    issubset(force, div_facets_) || throw(ArgumentError(
        "find_np_continuation_forced: force ⊄ div_facets"))
    nd = length(div_facets_)
    nd > max_div_facets && throw(B2Refusal(:ContinuationSearchTooLarge,
        Dict{String,Any}("count" => nd, "limit" => Int(max_div_facets),
                         "div_facets" => copy(div_facets_))))
    pool = setdiff(div_facets_, force)             # ascending (order-preserving)
    scanned = 0
    for k in 0:length(pool)
        for Sx in lex_combinations(pool, k)
            scanned += 1
            scanned > subset_budget && throw(B2Refusal(:ContinuationSearchTooLarge,
                Dict{String,Any}("scanned" => scanned,
                                 "subset_budget" => Int(subset_budget),
                                 "count" => nd)))
            S = sort!(vcat(collect(Int, force), collect(Int, Sx)))
            rest = setdiff(div_facets_, S)
            sigma = sigma_div(td, rest)
            ws = try
                produce_ws(td, sigma)
            catch err
                (err isa GeometricPropertyViolated ||
                 err isa WNotUnitNormalized) && continue
                rethrow()
            end
            all(_w_is_simple(w, multi_degree_np) for w in values(ws)) || continue
            return (rays = S, ws = ws, sigma = sigma, scanned = scanned)
        end
    end
    error("find_np_continuation_forced: unreachable — removing all divergent rays always succeeds")
end

# ===========================================================================
# §C7. Orchestration — the NP block of STExpandIntegral (wl:11093-11135).
# ===========================================================================
"""
    expand_integral_b2(E::EulerIntegrand, td::TropicalData;
                       order = 16, kwargs...) -> Vector{B2Continued}

The STExpandIntegral continuation orchestration (wl:11093-11135), stopping
at the B1 handoff: returns the list of continued integrands (exact typed
prefactors, dlog chart) READY for the B1 subtraction path — feed each
`.integrand` through the B1 census (C4 polytope rebuild + divergence_data)
and subtraction_terms; the census MUST be recomputed per output unless
`provenance["trdata_reusable"]` is true (poly support unchanged — the .wl
reuses trData only in the single-output case, wl:11140, and forces a
recompute otherwise).

Steps (regulator kept exact throughout):
 1. trop values on td.rays (tropical.jl); refusal ladder head as in C5:
    :TropIllDefined (wl:11071-11079 abort guard), :DegeneratePolytope.
 2. locally finite (no divergent ray): identity wrap, no continuation
    (wl:11090-11092 branch).
 3. per-ray divergence orders a = TropI.a ("ordersOfDivergences",
    wl:11085); non-integer positive a refuses :NonIntegerDivergenceOrder
    (the repeat count Table[ray, a(+1)] needs a ∈ Z).
 4. minimal-removal search: `find_first_np_continuation`.
 5. ray multiplicity ledger (wl:11110-11122):
      ray ∈ rays_to_continue      -> (a + 1) steps  ("made finite":
                                     deleted from the divergent fan;
                                     a = 0 log rays get exactly 1 step)
      power ray ∉ rays_to_continue-> a steps        ("made log-divergent")
 6. guard: total steps > `order` ⇒ B2Refusal(:ContinuationDepthExceeded)
    (`order` caps the TOTAL STEP COUNT length(np_list), summed over all
    rays — NOT a Laurent order and NOT per-ray; the per-ray multiplicities
    are always DERIVED from the census as above, never guessed; the paper
    warns the iteration is exponential in these counts, Sec 3.2.2. The
    driver surfaces this kwarg as `continuation_order` — DESIGN_B1 §4c.3).
 7. sequential continuation (`continue_rays`), provenance-chained.

Every output integrand keeps E.prefactor untouched and carries its own
B2Prefactor (the wl prefactor slot). The per-output ε-order bookkeeping the
.wl prints (ordersPref, wl:11126-11128) is recoverable exactly from the
B2Prefactor den factors (each log-ray step contributes one ε-linear pole
factor); it is not duplicated here."""
function expand_integral_b2(E::EulerIntegrand, td::TropicalData;
                            order::Integer = 16,
                            start_from::Integer = 0,
                            multi_degree_np::Integer = 4,
                            max_div_facets::Integer = 16,
                            subset_budget::Integer = 100_000)
    n = length(E.vars)
    td.ambient_dim == n || throw(ArgumentError(
        "expand_integral_b2: TropicalData ambient_dim != length(vars)"))
    trvals = trop_values(E, td)                    # tropical.jl (consumed)
    for (i, t) in enumerate(trvals)                # wl:11071-11079 guard
        is_illdefined(t) && throw(B1Refusal(:TropIllDefined, Dict{String,Any}(
            "ray_index" => i, "ray" => td.rays[i], "trvals" => trvals)))
    end
    isempty(td.equations) || throw(B1Refusal(:DegeneratePolytope,
        Dict{String,Any}("equations" => td.equations,
                         "ambient_dim" => td.ambient_dim)))
    sel = div_facets(trvals)                       # tropical.jl (consumed)
    if isempty(sel)                                # wl:11090-11092
        prov = Dict{String,Any}("np" => "locally finite — no continuation",
                                "trdata_reusable" => true)
        return B2Continued[B2Continued(b2pref_identity(EpsExp), E, prov)]
    end
    # divergence orders (wl:11085); integer check for the repeat counts
    orders = Dict{Int,Rational{BigInt}}(i => trvals[i].a for i in sel)
    for i in sel
        (orders[i] > 0 && denominator(orders[i]) != 1) &&
            throw(B2Refusal(:NonIntegerDivergenceOrder, Dict{String,Any}(
                "ray" => i, "tropI" => trvals[i], "order" => orders[i])))
    end
    power = [i for i in sel if trvals[i].a > 0]
    fr = find_first_np_continuation(td, sel; start_from = start_from,
        multi_degree_np = multi_degree_np, max_div_facets = max_div_facets,
        subset_budget = subset_budget)
    # ray multiplicity ledger (wl:11110-11122)
    np_list = Int[]
    for r in fr.rays                               # deleted rays: a+1 steps
        append!(np_list, fill(r, Int(orders[r]) + 1))
    end
    for r in power                                 # residual power: a steps
        r in fr.rays && continue
        append!(np_list, fill(r, Int(orders[r])))
    end
    if isempty(np_list)                            # wl:11124 NPcontinuation==={}
        prov = Dict{String,Any}("np" => "GP holds, log-only — no continuation needed",
                                "trdata_reusable" => true,
                                "np_search" => Dict{String,Any}(
                                    "rays_removed" => Int[],
                                    "subsets_scanned" => fr.scanned))
        return B2Continued[B2Continued(b2pref_identity(EpsExp), E, prov)]
    end
    length(np_list) > order && throw(B2Refusal(:ContinuationDepthExceeded,
        Dict{String,Any}("steps_required" => length(np_list),
                         "order" => Int(order),
                         "np_list" => copy(np_list),
                         "hint" => "raise `order` deliberately — each step multiplies the integrand count by up to #polys (paper Sec 3.2.2 bottleneck 1)")))
    state = continue_rays(E, [td.rays[r] for r in np_list])
    # provenance stamps: search result + trData reusability
    for ci in state
        ci.provenance["np_search"] = Dict{String,Any}(
            "rays_removed" => collect(fr.rays),
            "residual_power_rays" => [r for r in power if !(r in fr.rays)],
            "np_ray_sequence" => copy(np_list),
            "subsets_scanned" => fr.scanned)
        ci.provenance["trdata_reusable"] =
            all(!get(st, "deriv_poly_kept", true) for st in ci.provenance["steps"])
    end
    return state
end
