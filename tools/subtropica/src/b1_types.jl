# b1_types.jl — B1 SHARED TYPES. Companion to src/types.jl. Every B1 module
# (see DESIGN_B1.md) codes against these definitions as they are; a change
# here is a package-level design decision that touches every B1 module.
#
# Include order: AFTER src/types.jl (needs EpsExp, QQMPolyRingElem in scope).
# Like types.jl, this file is include()d inside `module SubTropica` after
# `using Nemo`, `using JSON` — do not add `using` here. Standalone tests
# (default global Julia env) do:
#     using Nemo, JSON
#     include(".../src/types.jl"); include(".../src/b1_types.jl")
#
# Reference specification: SPEC_B1 (transliteration-grade; see DESIGN_B1.md).
# The TWO WARNING pins (CONTRACTS.md §c) are LAW here:
#   WARNING-1 — net counterterm sign = (−1)^|subface|
#               (SubTropica.wl:10823 `Times[-a[subface],…]`, wl:10833
#                `a[s]=(−1)^(|s|+1)` ⇒ net −(−1)^(|s|+1) = (−1)^|s|).
#   WARNING-2 — I = FIRST det≠0 ray-support index subset (wl:10777-10782);
#               J = Complement(1:n, I); Vol completion rows, gauge-fix
#               (vars ∉ J → 1, wl:10834) and /∏x_J division (wl:10837)
#               ALL use J, the complement.
# Refusal contract: SPEC_B1 §6 (six kinds, checked in order). Sentinel logic
# is the CORRECTED one — upstream's guard `Count[us,"NotFound!"]>0`
# (wl:10760) is DEAD CODE (STProduceUs emits "NotFound", no bang, wl:10638;
# SPEC_B1 D1). Here there are NO string sentinels at all: a missing w is
# `missing`, and refusals are the typed exceptions below.

# ---------------------------------------------------------------------------
# EpsExp frozen arithmetic (ε-linear ring ops shared by ALL B1 builders).
# types.jl defines only the struct + ==/hash; builders need +,−,scalar-* for
# TropI evaluation (SPEC_B1 §3.1) and u_σ/Vol assembly (§4 steps 6-7).
# Defining them ONCE here prevents per-builder method collisions at
# integration. EpsExp*EpsExp is DELIBERATELY UNDEFINED — ε² has no home in
# an ε-linear exponent; a MethodError there is a bug detector, not a gap.
# ---------------------------------------------------------------------------
Base.zero(::Type{EpsExp}) = EpsExp(0//1, 0//1)
Base.zero(::EpsExp)       = zero(EpsExp)
Base.one(::Type{EpsExp})  = EpsExp(1//1, 0//1)
Base.one(::EpsExp)        = one(EpsExp)
Base.iszero(x::EpsExp)    = iszero(x.a) && iszero(x.b)
Base.:+(x::EpsExp, y::EpsExp) = EpsExp(x.a + y.a, x.b + y.b)
Base.:-(x::EpsExp)            = EpsExp(-x.a, -x.b)
Base.:-(x::EpsExp, y::EpsExp) = EpsExp(x.a - y.a, x.b - y.b)
Base.:+(x::EpsExp, c::Union{Integer,Rational}) = EpsExp(x.a + c, x.b)
Base.:+(c::Union{Integer,Rational}, x::EpsExp) = x + c
Base.:-(x::EpsExp, c::Union{Integer,Rational}) = EpsExp(x.a - c, x.b)
Base.:-(c::Union{Integer,Rational}, x::EpsExp) = EpsExp(c - x.a, -x.b)
Base.:*(x::EpsExp, c::Union{Integer,Rational}) = EpsExp(x.a * c, x.b * c)
Base.:*(c::Union{Integer,Rational}, x::EpsExp) = x * c

# ---------------------------------------------------------------------------
# RayClass — exact classification of TropI(ρ) = a + b·ε at regulators→0
# (SPEC_B1 §3.2; paper eq 3.14, paper:941-946; .wl divFacets wl:10739 and
# power split wl:10970/11093). TropI ≡ 0 (a==b==0) is NOT a class — it is
# the :TropIllDefined refusal (.wl abort guard wl:10948-10956/11071-11079).
#   convergent      : a < 0
#   log_divergent   : a == 0, b ≠ 0
#   power_divergent : a > 0            (B1 refuses — no C6 continuation)
# ---------------------------------------------------------------------------
@enum RayClass convergent log_divergent power_divergent

# Frozen classification predicates (SPEC_B1 §3.2, verbatim semantics —
# exact Rational comparisons; nothing float, nothing Reduce-shaped):
is_illdefined(t::EpsExp) = iszero(t.a) && iszero(t.b)   # wl:10948-10956 guard
is_divergent(t::EpsExp)  = t.a >= 0                     # divFacets: v ≥ 0
is_power(t::EpsExp)      = t.a > 0                      # refusal §6.3
is_log(t::EpsExp)        = iszero(t.a) && !iszero(t.b)

"""classify_ray(t) — RayClass of an exact TropI value; throws the
:TropIllDefined refusal (SPEC_B1 §6.1) on t ≡ 0. Check ORDER in callers:
ill-defined → degenerate polytope → power → GP (SPEC_B1 §3.2)."""
function classify_ray(t::EpsExp)::RayClass
    is_illdefined(t) &&
        throw(B1Refusal(:TropIllDefined, Dict{String,Any}("tropI" => t)))
    is_power(t) && return power_divergent
    is_log(t)   && return log_divergent
    return convergent
end

"""warning1_sign(subface) — the NET counterterm sign (−1)^|subface|
(WARNING-1; wl:10823/10833). Empty subface → +1 (identity term), single
facet → −1 (subtracted), pair → +1. The wrong v1 reading (−1)^(|s|+1)
ADDS single-facet counterterms; gate G3(c)/(d) catches it."""
warning1_sign(subface::AbstractVector{<:Integer})::Int =
    iseven(length(subface)) ? 1 : -1

# ---------------------------------------------------------------------------
# TropicalData — C4 output (SPEC_B1 §1 `TropData`, §2; mirror of .wl assData,
# wl:9856). Built by src/polytope.jl from the COEFFICIENT-FREE
# factor list over QQ[vars...] (integration vars ONLY — SPEC_B1 §2.1).
# Invariants (documented, asserted in polytope.jl's tests, not enforced here):
#   • rays[i] are OUTER facet normals, PRIMITIVE integer (lcm-cleared,
#     gcd-reduced, direction preserved — SPEC_B1 §2.2 pin; the w·ρ = −1
#     refusal check §3.3 is scale-sensitive).
#     Oscar halfspace form A·x ≤ b ⇒ outer normal = row of A, NO sign flip
#     (the .wl −FACETS[:,2:] flip at wl:9853 is polymake-convention-specific).
#   • vertices are lattice points, NO homogenizing leading 1 (we store
#     stripped; .wl stores with the 1 and strips at use, wl:10257).
#   • vertex_sets[i] = sorted indices of vertices maximizing v·rays[i]
#     (STGetFaces argmax, wl:10259-10264) — the ONLY geometric input Σ_div
#     needs (normal fan requested-then-discarded upstream, SPEC_B1 D7).
#   • equations = affine hull a·x = b rows; NONEMPTY ⇒ dim Newt < n ⇒
#     :DegeneratePolytope refusal (SPEC_B1 §6.2). Diagnostic like facet_rhs.
#   • length(rays) == length(vertex_sets) == length(facet_rhs);
#     all ray/vertex vectors have length ambient_dim.
# ---------------------------------------------------------------------------
struct TropicalData
    rays        :: Vector{Vector{BigInt}}    # OUTER normals, primitive integer
    vertices    :: Vector{Vector{BigInt}}    # lattice vertices, stripped
    vertex_sets :: Vector{Vector{Int}}       # per ray: argmax vertex indices, sorted
    facet_rhs   :: Vector{Rational{BigInt}}  # b_i of rays[i]·x ≤ b_i (diagnostic)
    equations   :: Vector{Tuple{Vector{Rational{BigInt}},Rational{BigInt}}}
                                             # affine hull rows (nonempty ⇒ refusal)
    ambient_dim :: Int                       # n = length(vars)
end

# ---------------------------------------------------------------------------
# DivergentFacet — C5 per-ray record (SPEC_B1 §1 `DivData`, flattened to one
# record per DIVERGENT ray; built by src/tropical.jl).
#   ray   : GLOBAL index into TropicalData.rays (all B1 face/subface indexing
#           is by these global indices — SPEC_B1 §3.2 ordering law).
#   class : log_divergent or power_divergent (power ⇒ the builder throws
#           PowerDivergentRefusal BEFORE any CounterTerm exists; a
#           power_divergent record may still be constructed transiently for
#           the refusal's detail payload).
#   tropI : exact TropI(ρ) = a + b·ε (EpsExp; regulators KEPT — the jacobian
#           exponent 1 − TropI(ρ) is ε-dependent, RT-extra-jac).
#   w     : the sign-flipped PRIMITIVE damping vector with w·ρ == −1
#           (SPEC_B1 §3.3 pin, discrepancy D2), or `missing` when the
#           nullspace search found none — the CORRECTED "NotFound" sentinel
#           (D1: upstream's "NotFound!" guard is dead; here missing-ness is
#           typed and the builder must throw GeometricPropertyViolated).
#           u_f = 1 − 1/(1+x^w) = x^w/(1+x^w) is DERIVED from w on demand
#           (CODE convention wl:10636, NOT the paper's v = 1/(1+x^w); D3).
# ---------------------------------------------------------------------------
struct DivergentFacet
    ray   :: Int
    class :: RayClass
    tropI :: EpsExp
    w     :: Union{Vector{BigInt},Missing}
end

# ---------------------------------------------------------------------------
# SigmaFace — one face of Σ_div (SPEC_B1 §3.2; STGetFaces wl:10255-10277;
# paper eq 3.24). Built by tropical.jl.
#   facets          : ASCENDING global ray indices; Int[] = the identity
#                     (empty) face, which is ALWAYS sigma_div[1] (wl:10272).
#                     Face list order = size, then lex in div_facets order
#                     (Mathematica Subsets order, wl:10267-10268).
#   common_vertices : sorted vertex indices in ⋂_{i∈facets} vertex_sets[i] —
#                     the compatibility WITNESS (nonempty ⇔ the subset is a
#                     genuine face; wl:10272 common-vertex test). Convention
#                     for the empty face: ALL vertex indices (intersection
#                     over the empty family), i.e. collect(1:nvertices).
# ---------------------------------------------------------------------------
struct SigmaFace
    facets          :: Vector{Int}
    common_vertices :: Vector{Int}
end

# ---------------------------------------------------------------------------
# CounterTerm — ONE Möbius term of the subtraction formula (C7 output unit;
# SPEC_B1 §1 CTTerm+FaceTerm flattened self-contained; built by subtract.jl,
# src/subtract.jl; consumed by C8). subtraction(...) returns
# Vector{CounterTerm}; grouping by `face` recovers SPEC_B1's FaceTerm list,
# ordered like sigma_div, terms per face in sigma_div-excess order.
#
#   face      : element of Σ_div this term belongs to (ascending global ray
#               indices; [] = the finite remainder face).
#   subface   : the EXCESS σ'∖face of the superset face σ' ⊇ face that
#               generated this term (wl:10807-10812); [] = identity term.
#   sign      : NET sign = warning1_sign(subface) = (−1)^|subface| —
#               WARNING-1. Stored explicitly (not recomputed) so the G3(d)
#               flipped-sign NEGATIVE control can build a wrong term at the
#               evaluation layer; builders MUST fill it via warning1_sign.
#   vol_det   : |det(rows: rays[face], then e_j for j ∈ J)| ∈ BigInt, ≥ 1
#               (wl:10792-10797; paper eq 3.35 is the transpose — D4). This
#               is the EXACT RATIONAL part of Vol.
#   vol_trops : the (−TropI(ρ)) factors, ρ ∈ face, order of `face`;
#               Vol(face) = vol_det / prod(vol_trops) — kept UNEXPANDED:
#               on log rays −TropI = −b·ε, so Vol ∝ 1/ε^|face| carries the
#               ONLY explicit ε-poles of B1 (paper:1219-1221); expansion is
#               C8's job (truncated-Laurent pole guard applies there).
#               Empty face ⇒ vol_det = 1, vol_trops = [] (Vol(∅) = 1).
#   J         : WARNING-2 — J = Complement(1:n, I), I = FIRST (lex) index
#               subset with det(rays[face][:,I]) ≠ 0 (wl:10777-10782;
#               eq 3.34). Empty face ⇒ J = collect(1:n). Used by all three:
#               Vol completion rows, gauge-fix, /∏x_J.
#   face_vars : vars[J], in `vars` order — the surviving integration
#               variables (flat measure on [0,∞)^|J|).
#   mono      : exponents of face_vars (length == length(J), same order),
#               AFTER the vars∉J → 1 gauge-fix (wl:10834) and AFTER the
#               /∏x_J flat-dx shift (wl:10837) and INCLUDING the x^{u_σ}
#               extra monomial (below). ε-dependent (EpsExp).
#   u_sigma   : the x^{u_σ} twist exponents (SPEC_B1 §4 step 6; eq 3.36):
#               u_σ = Σ_{ρ∈face} TropI(ρ)·π_J(w_ρ), length n (AMBIENT
#               indexing, zeros off J), components O(ε) on log rays. Kept
#               separately for auditability even though already folded into
#               `mono`; identity face ⇒ all-zero.
#   polys     : (P, e) factors of the term — the RESTRICTED P_j's (initial
#               forms along the face+subface rays, ORIGINAL coefficients
#               kept, wl:10506-10534) AND the subface jacobian factors
#               (1+x^{w_j}, −(1 − TropI(ρ_j))) for j ∈ subface — the
#               ε-DEPENDENT exponents of RT-extra-jac. Polys live in the
#               integrand's ring, already gauge-fixed (vars∉J → 1).
# ---------------------------------------------------------------------------
struct CounterTerm
    face      :: Vector{Int}
    subface   :: Vector{Int}
    sign      :: Int
    vol_det   :: BigInt
    vol_trops :: Vector{EpsExp}
    J         :: Vector{Int}
    face_vars :: Vector{Symbol}
    mono      :: Vector{EpsExp}
    u_sigma   :: Vector{EpsExp}
    polys     :: Vector{Tuple{QQMPolyRingElem,EpsExp}}
end

# ---------------------------------------------------------------------------
# Typed refusals — SPEC_B1 §6: refuse LOUDLY, typed, BEFORE any CounterTerm
# is emitted; no partial-output mode; the farm surfaces these verbatim.
# Three named types for the D1/D2-corrected checks + one generic carrier for
# the remaining kinds (:TropIllDefined, :DegeneratePolytope,
# :SigmaDivTooLarge, :NoBasisCompletion). All carry an ALWAYS-POPULATED
# detail Dict (rays/faces/TropI values involved).
# ---------------------------------------------------------------------------
abstract type B1RefusalException <: Exception end

"""PowerDivergentRefusal — some ray has TropI leading part a > 0 (SPEC_B1
§6.3, D10). The B1 route itself has no Nilsson-Passare continuation — that is
the B2 route (src/continue.jl), engaged by subtropica_integrate with
allow_continuation=true, which consumes this refusal (DESIGN RT-ix);
upstream merely PRINTS a warning (wl:10732-10734) and proceeds to garbage —
B1 throws instead. detail: ray indices, TropI values, divergence orders."""
struct PowerDivergentRefusal <: B1RefusalException
    ray    :: Int                 # offending GLOBAL ray index (first hit)
    tropI  :: EpsExp              # its TropI (a > 0)
    detail :: Dict{String,Any}
end
PowerDivergentRefusal(ray::Integer, tropI::EpsExp) =
    PowerDivergentRefusal(Int(ray), tropI,
        Dict{String,Any}("ray" => Int(ray), "tropI" => tropI,
                         "order" => tropI.a))

"""GeometricPropertyViolated — the w-search (STProduceUs wl:10587-10648)
found NO primitive w with w·ρ_f ≠ 0 orthogonal to all compatible rays for
some divergent facet (SPEC_B1 §3.3/§6.4; paper GP, paper:1118-1121).
CORRECTED sentinel logic (D1): upstream's `Count[us,"NotFound!"]>0` guard
(wl:10760) can never fire ("NotFound" is emitted, wl:10638) — here the
missing w IS the sentinel and this exception IS the guard. detail: facet,
compatible set, sorted candidate nullspace, eq 3.37 diagnosis hint."""
struct GeometricPropertyViolated <: B1RefusalException
    facet      :: Int             # divergent facet (global ray index)
    compatible :: Vector{Int}     # its compatible divergent facets
    detail     :: Dict{String,Any}
end
GeometricPropertyViolated(facet::Integer, compatible::AbstractVector{<:Integer}) =
    GeometricPropertyViolated(Int(facet), Int.(compatible),
        Dict{String,Any}("facet" => Int(facet),
                         "compatible" => Int.(compatible),
                         "hint" => "all nullspace candidates orthogonal to rho_f (cf. paper eq 3.37: every v_rho NotFound)"))

"""WNotUnitNormalized — the selected primitive w has w·ρ_f ≠ −1 (SPEC_B1
§3.3 pin / §6.5; discrepancy D2: upstream STSubtractionFormula runs
normalizeQ=False, wl:10758, and can silently use w·ρ_f = −c, c > 1 — the
subtraction theorem is stated for the normalized case, so B1 refuses).
detail: facet, w, ρ_f, the product."""
struct WNotUnitNormalized <: B1RefusalException
    facet   :: Int                # divergent facet (global ray index)
    w       :: Vector{BigInt}     # selected primitive w (post sign-flip)
    rho     :: Vector{BigInt}     # ρ_f (primitive outer normal)
    product :: BigInt             # w·ρ_f (≠ −1)
    detail  :: Dict{String,Any}
end
WNotUnitNormalized(facet::Integer, w::AbstractVector, rho::AbstractVector,
                   product::Integer) =
    WNotUnitNormalized(Int(facet), BigInt.(w), BigInt.(rho), BigInt(product),
        Dict{String,Any}("facet" => Int(facet), "w" => BigInt.(w),
                         "rho" => BigInt.(rho), "product" => BigInt(product)))

"""B1Refusal — generic carrier for the remaining SPEC_B1 §6 kinds:
:TropIllDefined (§6.1), :DegeneratePolytope (§6.2), :SigmaDivTooLarge
(§6.6, length(div_facets) > 24), :NoBasisCompletion (§4 step 5 — cannot
happen for lin. indep. face normals; refuse loudly, never FirstCase-Missing).
detail ALWAYS populated."""
struct B1Refusal <: B1RefusalException
    kind   :: Symbol
    detail :: Dict{String,Any}
end

# Message law (SPEC_B1 §6): the two continuation-shaped kinds name the B2
# continuation route (C6: the refusal
# stays the DEFAULT surface; the message tells the caller the opt-in);
# the rest are plain.
const _NP_SUFFIX = " — Nilsson-Passare continuation required (B2: rerun " *
                   "subtropica_integrate with allow_continuation=true); REFUSING"

Base.showerror(io::IO, e::PowerDivergentRefusal) =
    print(io, "subtropica B1: PowerDivergentRay — ray ", e.ray,
              " has TropI = ", e.tropI.a, " + (", e.tropI.b, ")*eps, a > 0",
              _NP_SUFFIX)
Base.showerror(io::IO, e::GeometricPropertyViolated) =
    print(io, "subtropica B1: GeometricPropertyViolated — no valid w for facet ",
              e.facet, " (compatible: ", e.compatible, ")", _NP_SUFFIX)
Base.showerror(io::IO, e::WNotUnitNormalized) =
    print(io, "subtropica B1: WNotUnitNormalized — facet ", e.facet,
              ": w·rho = ", e.product, " ≠ -1 (w = ", e.w, ", rho = ", e.rho,
              "); REFUSING")
Base.showerror(io::IO, e::B1Refusal) =
    print(io, "subtropica B1: ", e.kind, " — ", e.detail, "; REFUSING")

# ---------------------------------------------------------------------------
# Fieldwise ==/hash (same pattern as types.jl; isequal so that two `missing`
# w's compare EQUAL — desired for frozen data records, mutation tests rely
# on it). RayClass is an @enum: ==/hash already canonical.
# ---------------------------------------------------------------------------
for T in (:TropicalData, :DivergentFacet, :SigmaFace, :CounterTerm,
          :PowerDivergentRefusal, :GeometricPropertyViolated,
          :WNotUnitNormalized, :B1Refusal)
    @eval Base.:(==)(x::$T, y::$T) =
        all(isequal(getfield(x, f), getfield(y, f)) for f in fieldnames($T))
    @eval Base.hash(x::$T, h::UInt) =
        hash(Tuple(getfield(x, f) for f in fieldnames($T)), hash($(QuoteNode(T)), h))
end
