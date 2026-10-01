# PLD.jl backend bridge (optional).
#
# Port of upstream's Julia-side integration (scrj.m:1234-1317): the
# `Solver->momentumPLD` lane hands the Landau system to PLD.jl — the
# principal-Landau-determinant solver, SOFIA's companion back end
# (https://mathrepo.mis.mpg.de/PLD/) — as
#
#     Δ = α_0 G_0 + ... + α_{n-2} G_{n-2} + G_{n-1}
#
# (the last α set to 1, upstream sendSystemToJulia/momentumPLD), over a
# coefficient ring QQ[kinematics] and a Laurent ring in the active Baikov
# variables and the α's, then calls `getSpecializedPAD`.
#
# PLD.jl and its stack (Oscar, HomotopyContinuation) are deliberately NOT
# dependencies of this package — upstream likewise loads them into a
# separate Julia session at run time (scrj.m:1236-1250). Here the caller
# passes the loaded module:
#
#     import PLD
#     sing, ring = pld_singularities(PLD, diagram)
#
# The bridge is duck-typed against the AbstractAlgebra generic-ring
# interface that Oscar implements, so it costs nothing when unused and is
# testable against a mock backend (see test/pld_tests.jl); the live leg of
# the test suite runs whenever PLD.jl is installed.
#
# Upstream's EulerDiscriminantQ cross-check (scrj.m:1328) is not wrapped.

"""
    PLDSystem

The Landau system in the shape PLD.jl consumes (upstream `prepareVariables`
+ `sendSystemToJulia`, scrj.m:1253-1307):

- `delta`: Δ = Σ αᵢ Gᵢ with the last α set to 1, as a polynomial in
  `ring` = QQ[kin..., active..., alphas...]
- `kin`: kinematic variable names (the coefficient ring of the backend)
- `active`: active (integration) variable names
- `alphas`: the α names — all `ngrams` of them, matching upstream
  `allActive`, although the last one no longer appears in `delta`
"""
struct PLDSystem
    delta::Any
    ring::Any
    kin::Vector{String}
    active::Vector{String}
    alphas::Vector{String}
    ngrams::Int
end

# ---------- generic polynomial transport ----------

_ratcoeff(c::Rational) = Rational{BigInt}(c)
_ratcoeff(c::Integer) = Rational{BigInt}(c)
_ratcoeff(c::Nemo.ZZRingElem) = BigInt(c) // BigInt(1)
_ratcoeff(c) = BigInt(Nemo.numerator(c)) // BigInt(Nemo.denominator(c))

"""
    transplant(p, images, unit)

Rebuild polynomial `p` term by term in another ring: each variable of
`parent(p)` is looked up BY NAME in `images` (a name => generator dict) and
each coefficient is carried over as a `Rational{BigInt}`. `unit` is the
target's one, used for constant terms. Only variables actually appearing in
`p` need an image. The target generators may live in different rings (e.g.
a coefficient ring and a Laurent ring over it, as in the PLD call): the
term products promote generically. This is the native replacement for
upstream's string round trip through `ExternalEvaluate` (scrj.m:1282-1287).
"""
function transplant(p, images::AbstractDict{String,<:Any}, unit)
    names = map(string, Nemo.symbols(Nemo.parent(p)))
    imgs = Any[get(images, nm, nothing) for nm in names]
    out = nothing
    for (c, ev) in zip(Nemo.coefficients(p), Nemo.exponent_vectors(p))
        t = nothing
        for (i, e) in enumerate(ev)
            e == 0 && continue
            imgs[i] === nothing &&
                error("transplant: no image for variable '$(names[i])'")
            f = imgs[i]^e
            t = t === nothing ? f : t * f
        end
        term = (t === nothing ? unit : t) * _ratcoeff(c)
        out = out === nothing ? term : out + term
    end
    return out === nothing ? zero(unit) : out
end

# ---------- Δ construction (Nemo side, backend-independent) ----------

"""
    pld_delta(polys, active_idx) -> PLDSystem

Assemble Δ = α_0 P_1 + ... + α_{n-2} P_{n-1} + P_n from a polynomial
system, in the given order (upstream dehomogenizes by setting the LAST α
to 1; scrj.m:1306). `active_idx` are the generator indices of the active
(integration) variables in `parent(polys[1])`; every other variable
appearing in the system is kinematic. The result lives in a fresh ring
QQ[kin..., active..., α_0...α_{n-1}].
"""
function pld_delta(polys::AbstractVector, active_idx::Vector{Int})
    isempty(polys) && error("pld_delta: empty polynomial system")
    R0 = Nemo.parent(first(polys))
    names0 = map(string, Nemo.symbols(R0))
    nv = length(names0)
    all(i -> 1 <= i <= nv, active_idx) ||
        error("pld_delta: active index out of range")
    appears = falses(nv)
    for p in polys, i in 1:nv
        appears[i] = appears[i] || !free_of(p, i)
    end
    activeset = Set(active_idx)
    kin = [names0[i] for i in 1:nv if appears[i] && !(i in activeset)]
    active = names0[active_idx]
    isempty(kin) &&
        error("pld_delta: no kinematic variables in the system — " *
              "the specialized discriminants would have nowhere to live")
    n = length(polys)
    alphas = ["α$i" for i in 0:n-1]
    used = Set(vcat(kin, active))
    for a in alphas
        a in used && error("pld_delta: α name '$a' collides with a system variable")
    end
    Re, ge = Nemo.polynomial_ring(Nemo.QQ, vcat(kin, active, alphas))
    byname = Dict{String,Any}(zip(vcat(kin, active, alphas), ge))
    unit = one(Re)
    lift(p) = transplant(p, byname, unit)
    delta = lift(polys[n])
    for i in 1:n-1
        delta += byname[alphas[i]] * lift(polys[i])
    end
    return PLDSystem(delta, Re, kin, active, alphas, n)
end

"""
    pld_system(d::Diagram; maxcut=true, loopedges=nothing, dotperm=0,
               rowperm=0) -> PLDSystem

The diagram's Landau system in PLD shape: LBL Baikov Gram determinants
(`prepare_landau_system`), sorted by term count as upstream does before
dispatching to the solver (`SortBy[gramList, Length]`, scr.m:1246), then
assembled into Δ by `pld_delta`.
"""
function pld_system(d::Diagram; maxcut::Bool=true, loopedges=nothing,
                    dotperm::Int=0, rowperm::Int=0)
    gramlist, active, _ = prepare_landau_system(d; maxcut, loopedges,
                                                dotperm, rowperm)
    isempty(gramlist) &&
        error("pld_system: the Landau system is empty for this diagram")
    return pld_delta(sort(gramlist; by=nterms), active)
end

# ---------- the backend call ----------

"""Namespace holding the ring constructors: `PLD.Oscar` when PLD carries
the binding (it does — PLD.jl loads Oscar), else the module itself."""
_ring_namespace(PLD::Module) =
    isdefined(PLD, :Oscar) ? getfield(PLD, :Oscar) : PLD

function _backend_polyring(NS::Module, names::Vector{String})
    if isdefined(NS, :polynomial_ring)
        R, g = NS.polynomial_ring(NS.QQ, names)
    elseif isdefined(NS, :PolynomialRing)
        R, g = NS.PolynomialRing(NS.QQ, names)
    else
        error("PLD bridge: backend namespace $(NS) has neither " *
              "polynomial_ring nor PolynomialRing")
    end
    return R, collect(g)
end

function _backend_laurentring(NS::Module, R, names::Vector{String})
    if isdefined(NS, :LaurentPolynomialRing)
        S, g = NS.LaurentPolynomialRing(R, names)
    elseif isdefined(NS, :laurent_polynomial_ring)
        S, g = NS.laurent_polynomial_ring(R, names)
    else
        error("PLD bridge: backend namespace $(NS) has neither " *
              "LaurentPolynomialRing nor laurent_polynomial_ring")
    end
    return S, collect(g)
end

"""Collect every ring element in a (possibly nested) backend result."""
function _harvest!(out::Vector{Any}, x)
    if x isa Nemo.RingElem
        push!(out, x)
    elseif x isa Union{Tuple,AbstractArray}
        for y in x
            _harvest!(out, y)
        end
    elseif x isa AbstractDict
        for y in values(x)
            _harvest!(out, y)
        end
    end
    return out
end

"""
    pld_discriminants(PLD, sys::PLDSystem; method=:sym, homogeneous=true,
                      high_prec=false, codim_start=-1, face_start=1,
                      single_face=false) -> (discs, ring)

Hand `sys` to PLD.jl's `getSpecializedPAD` and return the specialized
discriminants as polynomials in a Nemo ring QQ[kinematics], deduplicated
up to constants. Keyword defaults mirror the upstream SOFIA option
defaults (`PLDMethod->sym, PLDHomogeneous->true, PLDHighPrecision->false,
PLDCodimStart->-1, PLDFaceStart->1, PLDRunASingleFace->false`,
scr.m:1231); the call itself is `ComputeDiscriminants` (scrj.m:1289).
Where upstream only echoes the backend's printed output, the bridge
returns the actual polynomials.
"""
function pld_discriminants(PLD::Module, sys::PLDSystem; method::Symbol=:sym,
                           homogeneous::Bool=true, high_prec::Bool=false,
                           codim_start::Integer=-1, face_start::Integer=1,
                           single_face::Bool=false)
    NS = _ring_namespace(PLD)
    Rb, rg = _backend_polyring(NS, sys.kin)
    snames = vcat(sys.active, sys.alphas)
    Sb, sg = _backend_laurentring(NS, Rb, snames)
    images = Dict{String,Any}(zip(sys.kin, rg))
    for (nm, g) in zip(snames, sg)
        images[nm] = g
    end
    Δ = transplant(sys.delta, images, one(Sb))
    raw = PLD.getSpecializedPAD(Δ, rg, sg; high_prec, method, homogeneous,
                                codim_start, face_start, single_face)
    keep = Any[]
    for p in _harvest!(Any[], raw)
        if p isa Nemo.MPolyRingElem
            is_constant(p) || push!(keep, p)
        elseif p isa Nemo.FieldElem || p isa Nemo.ZZRingElem
            # a bare numeric constant carries no singularity content
        else
            error("pld_discriminants: unexpected ring element of type " *
                  "$(typeof(p)) in the backend result")
        end
    end
    Rk, kg = Nemo.polynomial_ring(Nemo.QQ, sys.kin)
    back = Dict{String,Any}(zip(sys.kin, kg))
    unit = one(Rk)
    return dedup_proportional([transplant(p, back, unit) for p in keep]), Rk
end

"""
    pld_singularities(PLD, d::Diagram; maxcut=true, factor_result=true,
                      kwargs...) -> (sing, ring)

Candidate Landau singularities of one diagram via the PLD.jl backend —
the native counterpart of upstream `SOFIASingularitiesSEED` with
`Solver->momentumPLD` (scr.m:1245), and the independent in-ecosystem
cross-check of the FastFubini lane (`singularities_seed`). `factor_result`
mirrors upstream `FactorResult->True`: the discriminants are split into
distinct irreducible factors. Solver keywords (`method`, `homogeneous`,
`high_prec`, `codim_start`, `face_start`, `single_face`) and route
keywords (`loopedges`, `dotperm`, `rowperm`) pass through.

    import PLD
    sing, ring = pld_singularities(PLD, dbox)
"""
function pld_singularities(PLD::Module, d::Diagram; maxcut::Bool=true,
                           loopedges=nothing, dotperm::Int=0, rowperm::Int=0,
                           factor_result::Bool=true, kwargs...)
    sys = pld_system(d; maxcut, loopedges, dotperm, rowperm)
    discs, Rk = pld_discriminants(PLD, sys; kwargs...)
    out = factor_result ? factor_list_unique(discs) : discs
    return [p for p in out if !is_constant(p) && !iszero(p)], Rk
end
