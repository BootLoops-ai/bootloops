# polytope.jl — C4: Newton-polytope layer via Oscar FIRST-CLASS API.
#
# Reference specification: SPEC_B1 §2 (+ §6.2 refusal,
# §7 gates G1/G7). Behavioral reference: tools/subtropica/reference/SubTropica.wl
# (READ-ONLY; line numbers cited per semantic below). Governing docs:
# DESIGN_B1.md §2 (ownership), CONTRACTS.md. This file replaces the .wl
# polymake-script path (STtropicalDataOLD wl:9712-9858,
# STtropicalDataBuildScript wl:9869-9954, STTropicalDataPrecompute
# wl:9960-10054) with Oscar's first-class polyhedral API (RT-iv):
# newton_polytope / minkowski_sum / facets / vertices / affine_hull.
# One raw-Polymake escape hatch is PERMITTED by spec but none is needed —
# every query below is first-class.
#
# Include context: inside `module SubTropica` after `using Nemo`, `using JSON`,
# and (from integration on, DESIGN_B1 §4 step 1) `using Oscar` — do NOT add
# `using` here. Standalone tests load Oscar themselves and include this file
# by path. If Oscar is absent, definitions still load; any polytope call
# refuses loudly via _c4_require_oscar() (never a silent UndefVarError).
#
# SIGN PIN (SPEC_B1 §2.2): Oscar's halfspace form is A·x ≤ b, so the OUTER
# facet normal IS the row of A — NO sign flip. The .wl's rays = −FACETS[:,2:]
# (wl:9853) is polymake-path-specific (polymake stores b + a·x ≥ 0, i.e.
# INNER-pointing a). Do not add a flip here.
#
# SCALE PIN (SPEC_B1 §1 note + §2.2): rays are PRIMITIVE INTEGER vectors
# (lcm-cleared, gcd-reduced, direction preserved). Vol and the divergence
# classification are ray-scale-invariant, but the w·ρ = −1 refusal check
# (SPEC_B1 §3.3, discrepancy D2) is NOT — primitivity is load-bearing.
#
# Oscar/Polymake is in-process and NOT thread-safe (known pitfall):
# everything here runs SERIALLY; parallelism
# only ever at the C11 process farm (SPEC_B1 §2.4).
#
# Determinism: Oscar/polymake facet and vertex emission order is a backend
# detail; downstream indexing is by OUR arrays (gate G1: "up to ray ORDER").
# We canonicalize by sorting vertices and rays lexicographically (ascending)
# so TropicalData is reproducible across runs/versions. The 1-var branch's
# pinned ray order [[-1],[1]] (wl:9701) coincides with lex order.

# ---------------------------------------------------------------------------
# Oscar-presence guard: polytope calls must fail LOUDLY when Oscar is not
# loaded into the including module (Phase-A SubTropica has no Oscar; it enters
# the project env at integration, DESIGN_B1 §4 step 1).
# ---------------------------------------------------------------------------
function _c4_require_oscar()
    isdefined(@__MODULE__, :newton_polytope) && return nothing
    error("subtropica C4: Oscar is not loaded in $(@__MODULE__) — the polytope " *
          "layer needs `using Oscar` (added to the subtropica env at integration, " *
          "DESIGN_B1 §4 step 1; standalone tests use the default global env). " *
          "REFUSING")
end

# ---------------------------------------------------------------------------
# Exact conversions (everything rational-exact; nothing float anywhere).
# ---------------------------------------------------------------------------
_qq_to_rat(x)::Rational{BigInt} =
    Rational{BigInt}(BigInt(numerator(x)), BigInt(denominator(x)))
_qq_to_rat(x::Rational)::Rational{BigInt} = Rational{BigInt}(x)
_qq_to_rat(x::Integer)::Rational{BigInt} = Rational{BigInt}(x)

"""
    primitive_integer(v) -> (r::Vector{BigInt}, s::Rational{BigInt})

Primitive integer vector in the direction of `v` (SPEC_B1 §2.2): multiply by
the lcm of denominators, divide by the gcd of the resulting integers,
PRESERVE direction. Returns the primitive vector `r == s·v` and the positive
scale `s`, so a halfspace `v·x ≤ b` rescales consistently to `r·x ≤ s·b`
(the facet_rhs kept alongside the ray — "keep consistent pair", SPEC_B1
§2.2 pseudocode). Zero vector is refused (a facet normal is never 0).
"""
function primitive_integer(v::AbstractVector)
    vr = Rational{BigInt}[_qq_to_rat(x) for x in v]
    all(iszero, vr) && error("subtropica C4: primitive_integer of the zero vector")
    den = lcm([denominator(x) for x in vr])          # clear denominators
    w = BigInt[numerator(x * den) for x in vr]
    g = gcd(w)                                       # g > 0 (w ≠ 0)
    r = BigInt[div(x, g) for x in w]
    return r, Rational{BigInt}(den, g)               # s = den/g > 0
end

# ---------------------------------------------------------------------------
# Coefficient forgetting — STforgetCoeffs (wl:10325/10350): every coefficient
# → 1. Only the MONOMIAL SUPPORT survives into the polytope layer; kinematic
# data lives in coefficients and is gone by construction (SPEC_B1 §2.1).
# Coefficient-free ⇒ all coefficients +1 ⇒ no cancellation in products ⇒
# Newt(∏ factors) = Minkowski sum of the Newt(factor) EXACTLY.
# ---------------------------------------------------------------------------
function forget_coeffs(p::QQMPolyRingElem)
    R = parent(p)
    B = MPolyBuildCtx(R)
    for ev in exponent_vectors(p)
        push_term!(B, one(base_ring(R)), ev)
    end
    return finish(B)
end

# support(p) — exponent vectors as exact BigInt vectors (used by the 1-var
# branch here; tropical.jl has its own copy of this contract).
_support(p::QQMPolyRingElem) =
    Vector{BigInt}[BigInt.(ev) for ev in exponent_vectors(p)]

# ---------------------------------------------------------------------------
# Input contract checks (SPEC_B1 §2.1): factors live in QQ[vars...] — the
# integration variables ONLY, in `vars` order. Feeding the full
# QQ[vars..., kinvars...] ring is a BUG (the polytope would grow kinvar
# axes) — refuse it here, loudly, by symbol comparison.
# ---------------------------------------------------------------------------
function _check_input(polys::Vector{QQMPolyRingElem}, vars::Vector{Symbol})
    isempty(polys) &&
        throw(ArgumentError("subtropica C4: empty factor list — no ring context; " *
                            "pass the (possibly constant) polynomials"))
    R = parent(polys[1])
    for p in polys
        parent(p) === R ||
            throw(ArgumentError("subtropica C4: factors live in different rings"))
        iszero(p) &&
            throw(ArgumentError("subtropica C4: zero polynomial has no Newton polytope"))
    end
    symbols(R) == vars ||
        throw(ArgumentError("subtropica C4: ring symbols $(symbols(R)) ≠ vars $(vars) — " *
                            "SPEC_B1 §2.1: ring must be QQ[vars...] (integration " *
                            "variables ONLY, in order); kinvars in the ring is a bug"))
    return R
end

# ---------------------------------------------------------------------------
# Cache keying (SPEC_B1 §2.4; monPolsTropicalData wl:10412-10427).
# Upstream: cacheKey = {Times @@ polsNoCoeffs, vars} (wl:10419, 9975) —
# a Mathematica EXPRESSION, i.e. the UNEXPANDED product. B1 pins (G7) that
# factor-list refinement (P vs its factors) must HIT, so we canonicalize:
#   key = ( string(forget_coeffs(∏ forget_coeffs(P_j))), vars )
# — the coefficient-free EXPANDED product, coefficient-forgotten AGAIN so
# monomial multiplicities from the expansion (e.g. the 2·x1 in
# (1+x1)(1+x1+x2)) cannot leak coefficient values into the key. Nemo's
# `string` of a QQMPolyRingElem is deterministic (canonical internal
# monomial order) — the key is content-addressed by the support.
# Pins (SPEC_B1 §2.4):
#   • keyed on the PRODUCT — factor-list order/refinement cannot split entries;
#   • `vars` WITH ORDER is part of the key (permuted vars = different frame);
#   • coefficient VALUES never enter the key (wl comment 10416-10418);
#   • cache is process-local; the C11 farm gets per-worker caches (Oscar is
#     not thread-safe — no shared mutable state);
#   • STTropicalDataPrecompute's batched concurrent polymake launch
#     (wl:9999-10034) is polymake-process-specific: B1 precompute is just a
#     serial in-process loop over cache misses.
# Constant factors are dropped BEFORE keying (DeleteCases[polsNoCoeffs, 1],
# wl:10415 — a coeff^ε "polynomial" degenerates to 1 and must not enter).
# ---------------------------------------------------------------------------
if !@isdefined(TROPICAL_DATA_CACHE)
    const TROPICAL_DATA_CACHE =
        Dict{Tuple{String,Vector{Symbol}},TropicalData}()
end

clear_tropical_data_cache!() = empty!(TROPICAL_DATA_CACHE)

# coefficient-free, constant-dropped factor list + canonical coeff-free product
function _cf_factors_and_product(polys::Vector{QQMPolyRingElem})
    R = parent(polys[1])
    factors_cf = [forget_coeffs(p) for p in polys]     # wl:10325/10350
    filter!(f -> !is_constant(f), factors_cf)          # wl:10415
    prodcf = isempty(factors_cf) ? one(R) :
             forget_coeffs(prod(factors_cf))           # canonical (see key pins)
    return factors_cf, prodcf
end

"""
    tropical_data_cache_key(polys, vars) -> (String, Vector{Symbol})

The content-addressed cache key (SPEC_B1 §2.4): canonical string of the
coefficient-free product of the coefficient-free, constant-dropped factors,
paired with `vars` (order-sensitive). Exposed for the G7 keying tests.
"""
function tropical_data_cache_key(polys::Vector{QQMPolyRingElem},
                                 vars::Vector{Symbol})
    _check_input(polys, vars)
    _, prodcf = _cf_factors_and_product(polys)
    return (string(prodcf), copy(vars))
end

# ---------------------------------------------------------------------------
# vertex_sets — per-ray argmax vertex indices (STGetFaces inner loop,
# wl:10259-10264: sps = vertices . r; Position[sps, Max[sps]]), computed ONCE
# here and cached inside TropicalData (SPEC_B1 §2.2). This is the ONLY
# geometric input Σ_div needs (the normal fan is requested-then-discarded
# upstream — discrepancy D7).
# ---------------------------------------------------------------------------
function _argmax_vertex_sets(verts::Vector{Vector{BigInt}},
                             rays::Vector{Vector{BigInt}})
    sets = Vector{Vector{Int}}(undef, length(rays))
    for (k, r) in enumerate(rays)
        dots = [sum(v[i] * r[i] for i in eachindex(r)) for v in verts]
        m = maximum(dots)
        sets[k] = sort!(findall(==(m), dots))
    end
    return sets
end

# ---------------------------------------------------------------------------
# The 1-variable ANALYTIC branch — transliteration of STtropicalData
# (wl:9697-9709; duplicated in STNewton wl:9679-9691). This is what actually
# runs upstream for n == 1: polymake is NEVER called.
#
# The OTHER upstream 1-var mechanism — the polymake dummy-projection path
# (STtropicalDataBuildScript wl:9869-9954: declares a dummy second variable
# `local_var_names<Polynomial>(qw(x_0 x_1))` + `+0*x1` dependence, then
# `projection($p, [2], revert => true)` at wl:9913/9921; same logic at
# wl:9745-9754/9770 in STtropicalDataOLD) — exists only because polymake's
# newton() builds each factor's polytope in the subspace of variables that
# actually occur. Oscar's newton_polytope(f) for f ∈ QQ[x_1..x_n] always
# lives in R^n, so the dummy machinery is unnecessary AND it is dead code in
# practice (STtropicalData routes n==1 to the analytic branch above it).
# Documented per SPEC_B1 §2.3(b); NOT ported.
#
# Transliteration pins (SPEC_B1 §2.3(a)), acting on the coeff-free PRODUCT:
#   • rays = [[-1],[1]] hard-coded (wl:9701; ray 1 = lower endpoint
#     "x ≥ minExp", ray 2 = upper endpoint "x ≤ maxExp");
#   • vertices = ALL support exponents (wl:9704-9705 keeps interior lattice
#     support points, NOT just the two hull endpoints; argmax then selects
#     endpoints anyway — parity kept). Upstream order = CoefficientRules
#     order; immaterial downstream, we pin ascending for determinism;
#   • vertex_sets by the same argmax rule (shared helper above);
#   • single-monomial poly ⇒ nonempty `equations` (wl:9703) ⇒ B1
#     :DegeneratePolytope refusal downstream (SPEC_B1 §2.3a/§6.2) — a 1-var
#     integrand whose only poly is a monomial is scaleless;
#   • facet-offset sign QUIRK upstream (wl:9704 writes {minExp, +1} where
#     polymake's row convention for x ≥ min would be (−min, +1)) — harmless
#     upstream because nothing downstream ever consumes "facets"/"equations"
#     (grep in SPEC_B1 §2.3; discrepancy D8). We store canonical Oscar-style
#     values facet_rhs = (−minExp, maxExp) for rays ([-1], [1]) and canonical
#     a·x = b equations; the quirk is cited, NOT replicated.
# ---------------------------------------------------------------------------
function _tropical_data_1var(prodcf::QQMPolyRingElem, vars::Vector{Symbol})
    sup = sort!([ev[1] for ev in _support(prodcf)])      # ascending, exact
    rays = [BigInt[-1], BigInt[1]]                       # wl:9701, hard-coded
    verts = Vector{BigInt}[[e] for e in sup]
    vertex_sets = _argmax_vertex_sets(verts, rays)       # argmax endpoints
    lo, hi = sup[1], sup[end]
    facet_rhs = Rational{BigInt}[-lo, hi]                # canonical (D8 quirk
                                                         # cited above, not kept)
    equations = length(sup) == 1 ?                       # wl:9703 (single
        [(Rational{BigInt}[1], Rational{BigInt}(sup[1]))] : # monomial ⇒ point)
        Tuple{Vector{Rational{BigInt}},Rational{BigInt}}[]
    return TropicalData(rays, verts, vertex_sets, facet_rhs, equations, 1)
end

# ---------------------------------------------------------------------------
# Multi-var path: Oscar first-class (SPEC_B1 §2.2 pseudocode; replaces
# wl:9712-9858). Minkowski sum over the coeff-free factors mirrors
# `$p = minkowski_sum($p, $p_k)` (wl:9764-9784).
# ---------------------------------------------------------------------------
function _tropical_data_multivar(factors_cf::Vector{QQMPolyRingElem},
                                 R, n::Int)
    _c4_require_oscar()
    P = isempty(factors_cf) ?
        convex_hull(matrix(QQ, zeros(Int, 1, n))) :      # all-constant input:
                                                         # Newt = {0} (degenerate
                                                         # ⇒ refusal downstream)
        reduce(minkowski_sum, (newton_polytope(f) for f in factors_cf))
                                                         # wl:9764-9784

    # --- vertices: lattice points, NO homogenizing leading 1 to strip (the
    #     .wl strip at wl:10257 is polymake-format-specific; SPEC_B1 §1).
    verts = Vector{Vector{BigInt}}()
    for v in vertices(P)
        w = Vector{BigInt}(undef, n)
        for i in 1:n
            isone(denominator(v[i])) ||
                error("subtropica C4: non-lattice vertex $(v) — Newton polytopes " *
                      "of integer-exponent polynomials must be lattice polytopes")
            w[i] = BigInt(numerator(v[i]))
        end
        push!(verts, w)
    end
    sort!(verts)                                         # lex, deterministic

    # --- facets: Oscar halfspace form A·x ≤ b ⇒ OUTER normal = row of A,
    #     NO sign flip (SIGN PIN above; the .wl −FACETS[:,2:] flip at wl:9853
    #     is polymake-convention-specific).
    A, b = halfspace_matrix_pair(facets(P))
    nf = nrows(A)
    rays = Vector{Vector{BigInt}}(undef, nf)
    rhs = Vector{Rational{BigInt}}(undef, nf)
    for i in 1:nf
        r, s = primitive_integer([A[i, j] for j in 1:n])
        rays[i] = r
        rhs[i] = _qq_to_rat(b[i]) * s                    # consistent pair:
    end                                                  # r·x ≤ s·b
    perm = sortperm(rays)                                # lex, deterministic
    rays = rays[perm]
    rhs = rhs[perm]

    vertex_sets = _argmax_vertex_sets(verts, rays)       # wl:10259-10264

    # --- affine hull ("equations", wl assData data[[-1]], wl:9856):
    #     a·x = b rows; NONEMPTY ⇔ dim Newt < n ⇒ :DegeneratePolytope
    #     refusal downstream (SPEC_B1 §6.2 — upstream stores but never
    #     inspects these; B1 refuses rather than silently under-check).
    equations = Tuple{Vector{Rational{BigInt}},Rational{BigInt}}[]
    for h in affine_hull(P)
        a = normal_vector(h)
        push!(equations, (Rational{BigInt}[_qq_to_rat(a[i]) for i in 1:n],
                          _qq_to_rat(negbias(h))))
    end
    sort!(equations, by = first)                         # deterministic

    td = TropicalData(rays, verts, vertex_sets, rhs, equations, n)
    # documented TropicalData invariants (b1_types.jl) — cheap, assert here:
    @assert length(td.rays) == length(td.vertex_sets) == length(td.facet_rhs)
    @assert all(length(r) == n for r in td.rays)
    @assert all(length(v) == n for v in td.vertices)
    return td
end

# ---------------------------------------------------------------------------
# tropical_data — the C4 entry point (SPEC_B1 §2.1/§2.2).
# ---------------------------------------------------------------------------
"""
    tropical_data(polys::Vector{QQMPolyRingElem}, vars::Vector{Symbol};
                  use_cache = true) -> TropicalData

Newton-polytope data of the Euler integrand's polynomial factor list
(SPEC_B1 §2; mirror of the .wl assData, wl:9856). `polys` is the factor
list over QQ[vars...] — integration variables ONLY (SPEC_B1 §2.1);
coefficients are forgotten internally (STforgetCoeffs, wl:10325/10350) and
constant factors dropped (wl:10415), so passing either the coeff-free image
or the original coefficient-carrying factors yields the identical result
AND the identical cache key (coefficient values never matter — that is the
point of the cache, wl:10416-10418).

Output invariants (b1_types.jl TropicalData docs): OUTER primitive-integer
rays (NO sign flip vs Oscar's A·x ≤ b rows); stripped lattice vertices;
per-ray argmax vertex_sets (wl:10259-10264); facet_rhs consistent with the
primitive rays; affine-hull `equations` (nonempty ⇒ caller must refuse
:DegeneratePolytope — see `refuse_if_degenerate`). Vertices and rays are
sorted lexicographically (deterministic; downstream indexes by these arrays,
gate G1 "up to ray ORDER"). n == 1 takes the analytic branch (wl:9697-9709),
Oscar untouched.

Serial only — Oscar/Polymake is in-process and NOT thread-safe.
"""
function tropical_data(polys::Vector{QQMPolyRingElem}, vars::Vector{Symbol};
                       use_cache::Bool = true)
    R = _check_input(polys, vars)
    n = length(vars)
    factors_cf, prodcf = _cf_factors_and_product(polys)
    key = (string(prodcf), copy(vars))                   # §2.4 (wl:10419)
    if use_cache && haskey(TROPICAL_DATA_CACHE, key)
        return TROPICAL_DATA_CACHE[key]                  # hit: same object
    end
    td = n == 1 ? _tropical_data_1var(prodcf, vars) :
                  _tropical_data_multivar(factors_cf, R, n)
    use_cache && (TROPICAL_DATA_CACHE[key] = td)
    return td
end

# ---------------------------------------------------------------------------
# Degeneracy refusal helper (SPEC_B1 §6.2). C4 emits `equations` as a
# diagnostic; the B1 check ladder (ill-defined → DEGENERATE → power → GP,
# SPEC_B1 §6) refuses on it. Callers (driver / C5) use this helper so the
# refusal payload is uniform. Returns td unchanged when full-dimensional.
# ---------------------------------------------------------------------------
function refuse_if_degenerate(td::TropicalData)
    isempty(td.equations) && return td
    throw(B1Refusal(:DegeneratePolytope, Dict{String,Any}(
        "equations"   => td.equations,
        "dim"         => td.ambient_dim - length(td.equations),
        "ambient_dim" => td.ambient_dim)))
end

# ---------------------------------------------------------------------------
# normal_fan — DIAGNOSTIC ONLY, behind a flag (SPEC_B1 §2.2; discrepancy D7:
# the .wl requests MAXIMAL_CONES from polymake, wl:9828/9938, parses them,
# wl:9852/10043, and then DROPS them from the production assData,
# wl:9856/10045-10047; only the legacy fan variants wl:10061/10123 keep one).
# Σ_div is built from vertex-set intersections (SPEC_B1 §3.2), never from
# the fan — B1 correctness must not depend on this function.
# ---------------------------------------------------------------------------
function normal_fan_diagnostic(polys::Vector{QQMPolyRingElem},
                               vars::Vector{Symbol}; enable::Bool = false)
    enable || error("subtropica C4: normal_fan_diagnostic is DIAGNOSTIC-ONLY " *
                    "(SPEC_B1 §2.2, D7 — upstream discards the fan; Σ_div " *
                    "uses vertex-set intersections). Pass enable=true only " *
                    "for by-hand inspection; REFUSING")
    _c4_require_oscar()
    _check_input(polys, vars)
    length(vars) >= 2 ||
        error("subtropica C4: no fan for n == 1 (analytic branch, wl:9697-9709)")
    factors_cf, _ = _cf_factors_and_product(polys)
    isempty(factors_cf) &&
        error("subtropica C4: all factors constant — no fan")
    P = reduce(minkowski_sum, (newton_polytope(f) for f in factors_cf))
    return normal_fan(P)
end
