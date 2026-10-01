# tropical.jl — C3+C5: exact tropical evaluation, ray classification, Σ_div,
# w-vector search, u-functions, initial-form restriction.
# Ownership: DESIGN_B1 §2. Reference specification: SPEC_B1 §3 (+ §6 refusal
# ladder). Behavioral reference: tools/subtropica/reference/SubTropica.wl
# (git adac2f72; cited as wl:NNNN — line numbers re-verified in SPEC_B1).
# Paper citations = SubTropica paper line/eq numbers (reference/ vendored copy; SPEC_B1 pins).
#
# Like types.jl / b1_types.jl, this file is include()d inside `module SubTropica`
# after `using Nemo`, `using JSON` — do NOT add `using` here. Include order:
# AFTER types.jl (EulerIntegrand, EpsExp) and b1_types.jl (TropicalData,
# SigmaFace, DivergentFacet, RayClass predicates, typed refusals).
# Standalone tests (default global env): using Nemo; include types.jl,
# b1_types.jl, then this file (test/test_b1_tropical.jl pattern).
#
# The TWO WARNING pins (CONTRACTS.md §c) do not live here (C7/subtract.jl), but
# WARNING-2's scale-sensitivity is why produce_ws enforces the w·ρ = −1 pin
# on PRIMITIVE rays (SPEC_B1 §2.2 pin + §3.3/D2).

# ===========================================================================
# §T1. Exact tropical evaluation (SPEC_B1 §3.1)
# ===========================================================================

# _xdot(ev, rho, n) — exact dot of the x-part (FIRST n entries) of an exponent
# vector with a ray; kinematic-variable exponent slots (> n) never enter
# (kinematics live in coefficients; constants tropicalize to 0, wl:10363 and
# the numeric-drop at wl:10377). Exact Rational{BigInt} throughout — nothing
# float, nothing Reduce-shaped (SPEC_B1 §3.2 preamble).
_xdot(ev::AbstractVector, rho::AbstractVector, n::Integer) =
    sum(Rational{BigInt}(ev[i]) * Rational{BigInt}(rho[i]) for i in 1:n;
        init = zero(Rational{BigInt}))

"""
    support_x(P, n) -> Vector{Vector{Int}}

x-part support of `P`: the distinct exponent vectors over the FIRST n
(integration) variables. This is the coefficient-forgotten support
(STforgetCoeffs wl:10325/10350 image) — coefficients never enter, so
`support_x` of P equals the support of its coefficient-free image
(SPEC_B1 §3.1 pin on STTropicalizeIntegrand, wl:10479-10481).
"""
function support_x(P::QQMPolyRingElem, n::Integer)
    iszero(P) && throw(ArgumentError("support_x: zero polynomial"))
    unique(Vector{Int}[exponent_vector(P, i)[1:n] for i in 1:length(P)])
end

"""
    trop_poly(P, rho, n) -> Rational{BigInt}

Tropicalization of the polynomial `P` evaluated on the ray `rho`:
max_{m ∈ supp P} m[1:n]·rho — the max-plus image of STtropicalize
(wl:10362-10381: coefficients dropped, Plus→Max, monomial → m·vars)
composed with STEvalRay (wl:10354-10356); paper eq 3.12 (definition) /
eq 3.20 (worked form). Julia has no symbolic Max, so we evaluate per ray
directly and exactly (SPEC_B1 §3.1). Only the polynomial branch of the
.wl num/den split (wl:10367-10380) exists here: B1 `EulerIntegrand.polys`
entries are QQMPolyRingElem, denominator ≡ 1 by type. Constant P (support
= {0}) gives 0 — constants tropicalize to 0 (wl:10363).
"""
function trop_poly(P::QQMPolyRingElem, rho::AbstractVector, n::Integer)
    iszero(P) && throw(ArgumentError("trop_poly: zero polynomial has no Newton polytope"))
    maximum(_xdot(exponent_vector(P, i), rho, n) for i in 1:length(P))
end

"""
    trop_on_ray(E, rho) -> EpsExp

TropI(rho) = ν·rho + Σ_j e_j · trop_poly(P_j, rho), e_j = a_j + b_j·ε —
SPEC_B1 §3.1; .wl trint (wl:10730-10731, built by STTropicalizeIntegrand
wl:10479-10481); paper eq 3.12 (definition), eq 3.14 (ε-linearity — exact
here by EpsExp arithmetic; the .wl //Factor at wl:10731 is subsumed).

NOTE (measure convention): ν enters DIRECTLY — all tropical bookkeeping is
in the dlog chart ∫∏dx_i/x_i (SPEC_B1 §0; .wl converts to flat measure
only at output, wl:10740-10745/10837). The §8 fixtures pin this: F1 has
ν = (ε,ε) and TropI(ρ₁) = −ε. (types.jl's EulerIntegrand comment says
"flat measure dx"; the B1 spec + fixtures use ν as the dlog exponent —
flagged as an open REQUEST; this file follows SPEC_B1 §3.1
verbatim.)

NOTE (paper-text slip): below eq 3.21 the paper states Trop I(ρ₁) = +ε for
the eq 3.19 example; plugging ρ₁ = (−1,−2) into the paper's own eq 3.20
(and into this function) gives −ε. The a-part — hence the classification —
is identical; code wins (SPEC_B1 preamble). Pinned in the tests.
"""
function trop_on_ray(E::EulerIntegrand, rho::AbstractVector)::EpsExp
    n = length(E.vars)
    length(rho) == n || throw(ArgumentError(
        "trop_on_ray: ray length $(length(rho)) != n = $n integration variables"))
    length(E.nu) == n || throw(ArgumentError(
        "trop_on_ray: malformed EulerIntegrand — length(nu) != length(vars)"))
    acc = zero(EpsExp)
    for i in 1:n
        acc += E.nu[i] * Rational{BigInt}(rho[i])       # EpsExp * Rational (b1_types)
    end
    for (P, e) in E.polys
        acc += e * trop_poly(P, rho, n)                 # e = a + b·ε, exact
    end
    return acc
end

"""
    trop_values(E, td) -> Vector{EpsExp}

TropI on every facet ray of `td` (wl:10731 trValues). Per-ray classification
is `classify_ray` (b1_types.jl, frozen — throws :TropIllDefined on TropI ≡ 0);
the full refusal ladder lives in `divergence_data`.
"""
function trop_values(E::EulerIntegrand, td::TropicalData)
    td.ambient_dim == length(E.vars) || throw(ArgumentError(
        "trop_values: TropicalData ambient_dim $(td.ambient_dim) != length(vars) $(length(E.vars))"))
    EpsExp[trop_on_ray(E, rho) for rho in td.rays]
end

"""
    div_facets(trvals) -> Vector{Int}

Indices of divergent rays: a ≥ 0 at regulators→0 — .wl divFacets
`Position[trValues, v ≥ 0]` (wl:10739/11082); ascending global ray index
by construction (SPEC_B1 §3.2). Callers run the §6 ladder first
(ill-defined/power are refusals, not members of a "divergent" list that
C7 may consume).
"""
div_facets(trvals::AbstractVector{EpsExp}) =
    Int[i for (i, t) in enumerate(trvals) if is_divergent(t)]

# ===========================================================================
# §T2. Σ_div — faces of the divergent fan (SPEC_B1 §3.2)
# ===========================================================================

"""
    lex_combinations(v, m) -> Vector{Vector{T}}

All m-element subsets of the ordered collection `v`, in lexicographic order
of positions — exactly Mathematica `Subsets[v, {m}]` order (the order
STGetFaces inherits at wl:10267-10268; SPEC_B1 §3.2 ordering law, G2 pin:
`Subsets[Range[4],{2}]` = 12,13,14,23,24,34 — asserted in tests).
"""
function lex_combinations(v::AbstractVector{T}, m::Integer) where {T}
    k = length(v)
    m < 0 && return Vector{T}[]
    m == 0 && return [T[]]
    m > k && return Vector{T}[]
    out = Vector{T}[]
    idx = collect(1:m)
    while true
        push!(out, v[idx])
        j = m
        while j >= 1 && idx[j] == k - m + j
            j -= 1
        end
        j == 0 && break
        idx[j] += 1
        for t in (j+1):m
            idx[t] = idx[t-1] + 1
        end
    end
    return out
end

"""
    sigma_div(td, sel) -> Vector{SigmaFace}

Σ_div — transliteration of STGetFaces (wl:10255-10277); paper eq 3.24:
a subset of divergent facets is a face iff their vertex argmax sets share a
common vertex (the wl:10272 common-vertex test on `vertexInFacet`,
wl:10259-10264 — precomputed as `td.vertex_sets`). `sel` = div_facets
(ascending GLOBAL ray indices).

Ordering law (SPEC_B1 §3.2, preserved EXACTLY — C7 hard-codes face 1 = ∅,
wl:10781/10796): the empty face FIRST (always a face, wl:10272), witness =
ALL vertex indices (intersection over the empty family — b1_types SigmaFace
convention); then nonempty subsets by SIZE, then lex in `sel` order
(Mathematica Subsets order, wl:10267-10268). Elements of each face ascend in
global ray index (subsets of the ascending `sel`; the .wl posRays index map,
wl:10258/10276, is identity here — we work in indices throughout).

Parity notes: `sel == []` returns SigmaFace[] (wl:10256 `If[rays==={}, {}]`;
the locally-finite early return is C7's step 0). The upstream re-filter
`divFaces = Cases[faces, SubsetQ[divFacets,v]]` (wl:10749/11084) is a no-op
(SPEC_B1 §3.2 pin) — nothing to port.

Guard: length(sel) > 24 ⇒ B1Refusal(:SigmaDivTooLarge) (SPEC_B1 §6.6).
Ladder kind 6, but it must fire BEFORE the exponential subset scan — so
operationally it precedes the GP/w-norm checks (noted in SPEC_B1 §3.2:
"guard with a loud size check rather than silence").
"""
function sigma_div(td::TropicalData, sel::Vector{Int})
    isempty(sel) && return SigmaFace[]                          # wl:10256
    issorted(sel) || throw(ArgumentError(
        "sigma_div: sel must ascend in global ray index (SPEC_B1 §3.2 ordering law)"))
    all(i -> 1 <= i <= length(td.rays), sel) || throw(ArgumentError(
        "sigma_div: sel contains an out-of-range ray index"))
    length(sel) > 24 && throw(B1Refusal(:SigmaDivTooLarge, Dict{String,Any}(
        "count" => length(sel), "div_facets" => copy(sel), "limit" => 24)))
    nv = length(td.vertices)
    faces = SigmaFace[SigmaFace(Int[], collect(1:nv))]          # empty face ALWAYS first
    for m in 1:length(sel)                                      # size, then …
        for S in lex_combinations(sel, m)                       # … lex (wl:10267-10268)
            common = td.vertex_sets[S[1]]
            for i in 2:m
                common = intersect(common, td.vertex_sets[S[i]])
                isempty(common) && break
            end
            isempty(common) && continue                         # wl:10272 test
            push!(faces, SigmaFace(S, sort(common)))
        end
    end
    return faces
end

# ===========================================================================
# §T3. w-vector search (SPEC_B1 §3.3; STProduceUs wl:10587-10648)
# ===========================================================================

# exact integer dot (w, ρ ∈ Z^n)
_dotbig(u::AbstractVector, v::AbstractVector) =
    sum(BigInt(u[i]) * BigInt(v[i]) for i in eachindex(u); init = zero(BigInt))

# primitive integer vector: clear denominators (lcm), divide by content (gcd),
# DIRECTION PRESERVED (SPEC_B1 §2.2 primitive_integer semantics — the w·ρ=−1
# pin is scale-sensitive, §1 scale-invariance note).
function _primitive_intvec(v::AbstractVector{Rational{BigInt}})
    den = foldl(lcm, denominator.(v); init = one(BigInt))
    w = BigInt[numerator(x * den) for x in v]
    g = foldl(gcd, w; init = zero(BigInt))
    return (iszero(g) || isone(g)) ? w : BigInt[div(x, g) for x in w]
end
_primitive_intvec(v::AbstractVector{BigInt}) =
    _primitive_intvec(Rational{BigInt}[x // 1 for x in v])

"""
    _integer_nullspace_rows(rows, n) -> Vector{Vector{BigInt}}

Primitive-integer basis of {w ∈ Q^n : r·w = 0 ∀ r ∈ rows} — the exact
rational nullspace of the compatible-ray matrix (rows = rays; Mathematica
`NullSpace` at wl:10620). Deterministic presentation: Nemo rref over QQ,
then ONE standard kernel vector per free column (+1 in the free slot,
−R[i,free] in pivot slots), ascending free-column order, each primitivized
(content 1, direction kept — free slot stays +1 unless primitivization is
trivial, which preserves it too). Parity with upstream is pinned at the
level of the SELECTED w (SPEC_B1 §3.3), not the basis presentation — the
score-sort + orthogonality filter + sign flip + the w·ρ = −1 refusal (D2)
make the outcome equivalent; Mathematica's basis presentation differs and
is NOT reproduced.
"""
function _integer_nullspace_rows(rows::AbstractVector{<:AbstractVector}, n::Integer)
    k = length(rows)
    k >= 1 || throw(ArgumentError(
        "_integer_nullspace_rows: empty row set (simplest-w branch handles no-compatible facets)"))
    M = matrix(QQ, k, n, [QQ(rows[i][j]) for i in 1:k for j in 1:n])
    r, R = rref(M)
    pivots = Int[]
    for i in 1:r
        push!(pivots, findfirst(j -> !iszero(R[i, j]), 1:n))
    end
    free = setdiff(1:n, pivots)
    basis = Vector{Vector{BigInt}}()
    for j in free
        v = zeros(Rational{BigInt}, n)
        v[j] = 1
        for (i, p) in enumerate(pivots)
            v[p] = -(BigInt(numerator(R[i, j])) // BigInt(denominator(R[i, j])))
        end
        push!(basis, _primitive_intvec(v))
    end
    return basis
end

"""
STScoreWvecs (wl:10582): score {max|w_i|, #entries attaining that max,
Σ|w_i|} — lexicographically lower = simpler. Used with a STABLE sort
(Mathematica SortBy is stable; we pin MergeSort) so ties keep basis order.
"""
function _w_score(v::AbstractVector)
    a = BigInt[abs(BigInt(x)) for x in v]
    mx = maximum(a)
    return (mx, count(==(mx), a), sum(a))
end

"""
    produce_ws(td, sigma) -> Dict{Int, Vector{BigInt}}

w-vector search — transliteration of STProduceUs (wl:10587-10648) +
STScoreWvecs (wl:10582); paper geometric property (paper:1118-1121): for
each divergent facet f, a w with w·ρ_f = −1 and w·ρ' = 0 for every
divergent ρ' compatible with f.

* facets/pairs are read off `sigma` (the singleton/pair faces; wl:10593-10594);
  compatibility comes from the PAIR faces ONLY (wl:10596) — Σ_div is a
  simplicial complex (paper:1122-1123), so pairs determine it (SPEC_B1 §3.3).
* no compatible partner ⇒ simplest-w branch (wl:10603-10607):
  k = first nonzero slot of ρ_f, w = −ρ_f[k]·e_k.
* else: exact rational nullspace of the compatible rays (wl:10610-10620),
  primitive basis, score-sorted (stable), candidates orthogonal to ρ_f
  removed (wl:10622). Empty result — either trivial nullspace
  (wl:10611-10616 zero-vector placeholder, geometricQ = False) or
  filtered-to-empty (wl:10626) — ⇒ GeometricPropertyViolated. This is the
  CORRECTED sentinel (SPEC_B1 D1): upstream's guard
  `Count[us,"NotFound!"]>0` (wl:10760) is DEAD CODE — STProduceUs emits
  "NotFound", no bang (wl:10638). Here missing-ness is TYPED and thrown;
  there are no string sentinels (b1_types preamble). Upstream's
  `First[{}]` message-and-continue path (wl:10626) is NOT reproduced.
* sign flip so that w·ρ_f < 0 (wl:10633; same at 10668 in STProduceAllUs).
* D2 PIN (SPEC_B1 §3.3/§6.5): the selected PRIMITIVE w must satisfy
  w·ρ_f == −1 EXACTLY, else WNotUnitNormalized — STRICTER than upstream:
  STSubtractionFormula calls with normalizeQ = False (wl:10758) and can
  silently proceed with w·ρ_f = −c, c > 1 (incl. the simplest-w branch,
  where w·ρ_f = −ρ_f[k]²). The subtraction theorem is stated for the
  normalized case; relaxation is B2 scope.

Facets are processed in ascending global index (singleton order in `sigma`);
with several offenders the FIRST throws. Requires `sigma` to obey the
SPEC_B1 §3.2 ordering law (as produced by `sigma_div`).
"""
function produce_ws(td::TropicalData, sigma::Vector{SigmaFace})
    n = td.ambient_dim
    facs = Int[f.facets[1] for f in sigma if length(f.facets) == 1]     # wl:10593
    prs = Vector{Int}[f.facets for f in sigma if length(f.facets) == 2] # wl:10594
    ws = Dict{Int,Vector{BigInt}}()
    for f in facs                                                       # wl:10600-10642
        rho = td.rays[f]
        compat = Int[]                                                  # wl:10596
        for p in prs
            f in p || continue
            for g in p
                (g == f || g in compat) || push!(compat, g)
            end
        end
        local wf::Vector{BigInt}
        if isempty(compat)
            # simplest-w branch (wl:10603-10607): w = −ρ_f[k]·e_k
            k = findfirst(!iszero, rho)
            k === nothing && throw(ArgumentError("produce_ws: zero ray $f (rays must be primitive nonzero)"))
            wf = zeros(BigInt, n)
            wf[k] = -BigInt(rho[k])
        else
            basis = _integer_nullspace_rows([td.rays[g] for g in compat], n)
            cands = sort(basis; by = _w_score, alg = Base.Sort.MergeSort) # wl:10620 SortBy, stable
            keep = [v for v in cands if _dotbig(v, rho) != 0]             # wl:10622
            if isempty(keep)                                              # GP violated (D1-corrected)
                throw(GeometricPropertyViolated(f, copy(compat), Dict{String,Any}(
                    "facet" => f, "compatible" => copy(compat),
                    "nullspace" => cands,
                    "hint" => "all nullspace candidates orthogonal to rho_f (cf. paper eq 3.37: every v_rho NotFound)")))
            end
            wf = _primitive_intvec(keep[1])   # already primitive; kept for SPEC_B1 §3.3 pseudocode parity
        end
        _dotbig(wf, rho) > 0 && (wf = -wf)                                # wl:10633
        pr = _dotbig(wf, rho)
        pr == -1 || throw(WNotUnitNormalized(f, wf, BigInt.(rho), pr))    # D2 pin, §6.5
        ws[f] = wf
    end
    return ws
end

# ===========================================================================
# §T4. u-functions (CODE convention — wl:10636/10669; SPEC_B1 §3.3 + D3)
# ===========================================================================

"""
    u_pair(xs, w) -> (num, den)::NTuple{2, polynomial}

u_f = 1 − 1/(1+x^w) = x^w/(1+x^w) in the Factor-ed polynomial-pair form
(`1-1/(1+Times@@(xvars^nullVector))//Factor`, wl:10636; same at wl:10669 in
STProduceAllUs): with the sign split w = w⁺ − w⁻ (componentwise),
num = x^{w⁺}, den = x^{w⁻} + x^{w⁺} — exactly Mma Factor's cleared form
(gcd(num,den) = 1: the two monomials x^{w⁺}, x^{w⁻} have disjoint variable
support). CODE convention, NOT the paper's damping function
v_ρ = 1/(1+x^w) (eq 3.28) — discrepancy D3: v = 1 − u, see
`one_minus_u_pair`; every C7 formula uses (1−u). B1 stores w natively
(DivergentFacet.w) and materializes u on demand; the .wl u→w recovery
(step 4, wl:10770-10775) is a parity round-trip in the tests only.
`xs` = the integration-variable generators (first n gens of the ring).
"""
function u_pair(xs::Vector{<:MPolyRingElem}, w::AbstractVector)
    isempty(xs) && throw(ArgumentError("u_pair: empty variable list"))
    length(xs) == length(w) || throw(ArgumentError(
        "u_pair: length(xs) = $(length(xs)) != length(w) = $(length(w))"))
    R = parent(xs[1])
    wp = one(R)
    wm = one(R)
    for i in eachindex(xs)
        wi = Int(w[i])
        wi > 0 && (wp *= xs[i]^wi)
        wi < 0 && (wm *= xs[i]^(-wi))
    end
    return (wp, wm + wp)
end

"""
    one_minus_u_pair(xs, w) -> (num, den)

1 − u = v_ρ = 1/(1+x^w) (paper eq 3.28) in cleared form:
(x^{w⁻}, x^{w⁻} + x^{w⁺}) — the jacobian factor base of C7 step 8
(wl:10801-10804), exponent 1 − TropI(ρ) kept ε-dependent there
(RT-extra-jac; not this file's concern).
"""
function one_minus_u_pair(xs::Vector{<:MPolyRingElem}, w::AbstractVector)
    num, den = u_pair(xs, w)
    return (den - num, den)
end

# ===========================================================================
# §T5. Initial-form restriction (SPEC_B1 §4 step 9 operator; exported utility
#      for subtract.jl — DESIGN_B1 §2 ownership table)
# ===========================================================================

"""
    restrict_poly(P, rho, n) -> polynomial

Initial form of `P` along the ray `rho`: keep EXACTLY the monomials whose
x-part maximizes m·rho — the face of Newt P exposed by the outer normal rho —
with ORIGINAL coefficients kept (only monomial SELECTION uses the support;
wl:10510 `(cr[[1]]/.coefsRules)`). Transliterates STrestrictPoly
(wl:10506-10511 — stated there as the min of −m·ρ with depth 0; identical;
`depth` is never nonzero in C7, so it is not ported). Paper eq 3.30,
"initial form" (paper:1168-1170). Monomial prefactor x^ν and exponents e_j
pass through unchanged at the integrand level (wl:10532 — C7's concern).
"""
function restrict_poly(P::QQMPolyRingElem, rho::AbstractVector{<:Union{Integer,Rational}},
                       n::Integer)
    iszero(P) && throw(ArgumentError("restrict_poly: zero polynomial"))
    m0 = maximum(_xdot(exponent_vector(P, i), rho, n) for i in 1:length(P))
    out = zero(parent(P))
    for i in 1:length(P)
        _xdot(exponent_vector(P, i), rho, n) == m0 && (out += term(P, i))
    end
    return out
end

"""
    restrict_poly(P, rays, n) -> polynomial

Sequential multi-ray (face) restriction: iterated initial forms, applied
LEFT-TO-RIGHT. The .wl applies face[2:end] first, then face[1]
(STrestrictIntegrand recursion wl:10534; STrestrict loop wl:10519-10525);
the result is order-independent — iterated initial forms along normals of a
common face commute (SPEC_B1 §4 step 9 restriction pin) — so B1 uses plain
left-to-right and the tests assert order-independence (G8).
"""
restrict_poly(P::QQMPolyRingElem, rays::AbstractVector{<:AbstractVector}, n::Integer) =
    foldl((q, r) -> restrict_poly(q, r, n), rays; init = P)

# ===========================================================================
# §T6. C5 driver — refusal ladder + divergence data (SPEC_B1 §6, §3.2-3.3)
# ===========================================================================

"""
    divergence_data(E, td) -> NamedTuple

C5 entry point: run the §6 refusal ladder over the tropical data of `E`
(TropicalData from C4), then build Σ_div and the w-vectors. Returns

    (trvals, classes, div_facets, sigma_div, w, facets)

with `facets :: Vector{DivergentFacet}` — one record per divergent ray, all
`log_divergent` with `w` filled (power NEVER returns: it throws; a
power_divergent DivergentFacet exists only transiently inside the refusal's
detail payload, per the b1_types doc). Empty `div_facets` ⇒ locally-finite
early return with empty sigma/w/facets — the flat-measure conversion of
SPEC_B1 §4 step 0 is C7's job, not done here.

Refusal ladder (SPEC_B1 §6 — first hit throws; all values exact):
 1. B1Refusal(:TropIllDefined)   — some TropI(ρ) ≡ 0 (a = b = 0); mirror of
    the .wl abort guard `MemberQ[trValues,0] → Abort[]`
    (wl:10948-10956/11071-11079). Typical cause: unfixed projective
    (Cheng-Wu) gauge. detail: ray index, ray, all trvals.
 2. B1Refusal(:DegeneratePolytope) — td.equations nonempty ⇒ dim Newt < n;
    facet normals do not cover lineality directions, the eq 3.13 criterion
    is incomplete ⇒ refuse rather than under-check (SPEC_B1 §6.2; upstream
    stores the affine hull but never inspects it). The authoritative
    construction-side test is polytope.jl's (G1/G6-iv); re-checked here so the C5
    entry cannot be driven past an under-checked polytope. detail:
    equations, dim, ambient_dim.
 3. PowerDivergentRefusal — some ray has a > 0 (SPEC_B1 §6.3/D10; paper
    945-946). Upstream merely PRINTS "Power divergence on …"
    (wl:10732-10734) and proceeds to garbage — B1 throws (the B1-route
    contract, DESIGN RT-ix; the shipped B2 route continue.jl consumes the
    refusal). detail collects ALL power rays for B2 to consume:
    indices, TropI values, orders a ("ordersOfDivergences", wl:10962), and
    transient DivergentFacet records (w = missing).
 4./5. GeometricPropertyViolated / WNotUnitNormalized — inside produce_ws.
 6. B1Refusal(:SigmaDivTooLarge) — inside sigma_div; ladder kind 6, fired
    before the exponential subset scan (SPEC_B1 §3.2 guard note).
"""
function divergence_data(E::EulerIntegrand, td::TropicalData)
    trvals = trop_values(E, td)
    # --- ladder 1: ill-defined (wl:10948-10956/11071-11079) ----------------
    for (i, t) in enumerate(trvals)
        is_illdefined(t) && throw(B1Refusal(:TropIllDefined, Dict{String,Any}(
            "ray_index" => i, "ray" => td.rays[i], "trvals" => trvals,
            "hint" => "TropI == 0 identically on a ray - typically an unfixed projective (Cheng-Wu) gauge")))
    end
    # --- ladder 2: degenerate polytope (SPEC_B1 §6.2) -----------------------
    if !isempty(td.equations)
        eqm = matrix(QQ, length(td.equations), td.ambient_dim,
                     [QQ(q[1][j]) for q in td.equations for j in 1:td.ambient_dim])
        throw(B1Refusal(:DegeneratePolytope, Dict{String,Any}(
            "equations" => td.equations, "dim" => td.ambient_dim - rank(eqm),
            "ambient_dim" => td.ambient_dim)))
    end
    classes = RayClass[classify_ray(t) for t in trvals]   # (re-throws ill-defined; unreachable)
    # --- ladder 3: power divergence (§6.3, D10) ------------------------------
    powers = findall(==(power_divergent), classes)
    if !isempty(powers)
        i0 = powers[1]
        throw(PowerDivergentRefusal(i0, trvals[i0], Dict{String,Any}(
            "ray" => i0, "tropI" => trvals[i0], "order" => trvals[i0].a,
            "rays" => powers,
            "tropIs" => EpsExp[trvals[i] for i in powers],
            "orders" => Rational{BigInt}[trvals[i].a for i in powers],
            "facets" => DivergentFacet[DivergentFacet(i, power_divergent, trvals[i], missing)
                                       for i in powers])))
    end
    # --- Σ_div (ladder 6 inside) + w-search (ladders 4/5 inside) ------------
    sel = div_facets(trvals)
    sigma = sigma_div(td, sel)
    ws = produce_ws(td, sigma)
    facets = DivergentFacet[DivergentFacet(i, classes[i], trvals[i], ws[i]) for i in sel]
    return (trvals = trvals, classes = classes, div_facets = sel,
            sigma_div = sigma, w = ws, facets = facets)
end
