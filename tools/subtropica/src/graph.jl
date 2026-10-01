# graph.jl — (DESIGN_B1.md §2) Phase-C graph front-end:
# propagator graph -> Symanzik U,F -> raw Euler quadruple (EulerIntegrand).
#
# This file is include()d inside `module SubTropica` after `using Nemo, JSON`
# (same rule as types.jl — no `using` here). Standalone tests include it by
# path after types.jl/b1_types.jl (DESIGN_B1 §3).
#
# Semantics transliterated from (cite-as-you-go below):
#  * SubTropica.wl STGetIntegrandData (wl:2975-3011): exponent conventions
#      eU = ν_tot − (L+1)(D−2ε)/2   (wl:2997, `eU=(ν-(L+1)(dim-2 eps)/2)`)
#      eF = −(ν_tot − L(D−2ε)/2)    (wl:2996)
#    prefactor Γ(−eF)/∏_e Γ(ν_e) (wl:3010), Feynman branch (−F)^eF with
#    −F ≥ 0 Euclidean (wl:3001-3005; also wl:3583-3591 `Fphys = -F`).
#  * SubTropica.wl STSymanzik (Options wl:3104, main wl:3122): the A,B,c
#    quadratic-form data Qβ/Lβ/cβ (wl:3474-3486), U = Det[Qβ] (wl:3498),
#    F = cβ·U − Lβ·Adj(Qβ)·Lβᵀ (wl:3500-3504), zeroU guard (wl:3510-3516),
#    zeroF/vacuumPeriod guards (wl:3519-3531), prefactor Γ(dd)/∏Γ(ν_e)
#    with dd = ν_tot − L·D/2 (wl:3474, wl:3596-3598).
#    We do NOT port the momentum-space Qβ/Lβ/cβ route: the input here is a
#    GRAPH, so U and F are built by the equivalent graph-combinatorial
#    formulas (weighted matrix-tree + spanning-2-forest), reusing the
#    live-validated pattern of
#    pld_complete_c3.jl (internal reference)
#    (build_UF: reduced weighted Laplacian for U, lines 25-51; two_forest_F
#    with union-find component split, lines 82-132; F = F0 + mass·U,
#    line 68). Equality of the two routes is the standard matrix-tree /
#    all-minors theorem; we additionally cross-check U against the direct
#    spanning-tree sum on every call (crosscheck kwarg).
#  * CenterDot NUMERATOR machinery (STSymanzik step 5 wl:3546-3568, step 7
#    postprocessing wl:3606-3660; STGetIntegrandDataNumerators wl:3014-3080):
#    NOT ported in this pass — numerators REFUSE, typed (:NumeratorsNotSupported),
#    per the B1 charter ("numerators OPTIONAL for this pass, refuse typed").
#  * Paper eq 1.3-1.5 (SubTropica paper, reference/ vendored copy, lines 222-249):
#    Schwinger rep I = Γ(ω)·∏_e (−1)^{ν_e}/Γ(ν_e) · ∫ dx/GL(1)
#    ∏ x^{ν_e−1} · U^{−(D/2−ω)} (−F)^{−ω}, ω = ν_tot − LD/2 = −eF.
#    Hence Prefactor.c = (−1)^{ν_tot}; the e^{εLγ_E} of eq 1.3 is OPTIONAL
#    here (gammaE_norm kwarg, default OFF to match the frozen fixture
#    convention — fixtures/moller_box_eq41.json provenance pins the
#    pySecDec-disteval normalization with gammaE_eps = 0).
#  * Lee-Pomeransky variant (lee_pomeransky kwarg; Lee & Pomeransky
#    arXiv:1308.6676; same G = U + F assembly as pld_complete_c3.jl:148):
#    I = (−1)^{ν_tot} Γ(D/2) / (Γ((L+1)D/2 − ν_tot) ∏Γ(ν_e)) ·
#        ∫_{[0,∞)^E} ∏ dx x^{ν_e−1} · G^{−D/2},  G = U + F(Euclid-positive),
#    NO GL(1) gauge fixing (projectively complete).
#
# Conventions pinned against the FROZEN Phase-A fixtures (independent checks):
#  * fixtures/moller_box_eq41.json: eU=(0,2ε·…)=EpsExp(0,2), eF=EpsExp(−2,−1),
#    Γ(2+ε), nu(flat)=ν_e−1, vars x0..x{E-2} (Cheng-Wu on the LAST edge,
#    0-based pySecDec naming), kinvars-after-vars ring order (types.jl).
#  * a 3-loop 8-propagator graph (ν_tot=8, L=3, in test_graph.jl):
#    eU=EpsExp(0,4), eF=EpsExp(−2,−3), Γ(2+3ε) = Γ(−eF). Both cases match
#    the formulas above with zero free parameters.
#
# Kinematic input model: masses[e] = the SQUARED mass of edge e (0, exact
# rational, kinvar Symbol, or a linear combination — see _graph_kinform);
# `kinematics` maps dot products p_i·p_j (1 ≤ i ≤ j ≤ nlegs−1; the LAST leg
# is eliminated by momentum conservation p_n = −Σ_{i<n} p_i, same route as
# pld_complete_c3.jl psq(), lines 86-97) to linear forms over kinvars.
# Supply them in the DECLARED-EUCLIDEAN chart (types.jl EulerIntegrand
# region contract): graph_to_quadruple refuses :NonEuclidean if any U/F
# coefficient comes out non-positive (same check as scripts/make_fixtures.py
# mma_poly require_positive).

# ---------------------------------------------------------------------------
# Typed refusal (graph front-end scope — B1Refusal kinds are C4/C5/C7's;
# this front-end refuses BEFORE any tropical machinery runs).
# Kinds: :NumeratorsNotSupported, :NumeratorExponent, :NoLoop, :Disconnected,
#        :Scaleless, :NonEuclidean, :MissingDotProduct, :BadInput, :TooLarge,
#        :InternalCheckFailed.
# ---------------------------------------------------------------------------
struct GraphRefusal <: Exception
    kind   :: Symbol
    detail :: Dict{String,Any}
end
Base.showerror(io::IO, e::GraphRefusal) =
    print(io, "subtropica graph front-end: ", e.kind, " — ", e.detail, "; REFUSING")

_graph_refuse(kind::Symbol, kv...) =
    throw(GraphRefusal(kind, Dict{String,Any}(kv...)))

# ---------------------------------------------------------------------------
# Kinematic linear forms: Dict{Symbol,Rational{BigInt}} with the constant
# part under :__const__. Accepted inputs (recursively): Integer/Rational
# (constant), Symbol (coefficient-1 kinvar), Pair{Symbol,<:Number}
# (coefficient·kinvar), Vector/Tuple of those (sum).
# ---------------------------------------------------------------------------
const _GRAPH_CONSTKEY = :__const__
const _GraphKinForm = Dict{Symbol,Rational{BigInt}}

function _graph_kf_add!(f::_GraphKinForm, x, mult::Rational{BigInt})
    if x isa Integer || x isa Rational
        f[_GRAPH_CONSTKEY] = get(f, _GRAPH_CONSTKEY, 0 // BigInt(1)) +
                             mult * Rational{BigInt}(x)
    elseif x isa Symbol
        x == _GRAPH_CONSTKEY && _graph_refuse(:BadInput,
            "msg" => "kinvar name collides with the internal constant key")
        f[x] = get(f, x, 0 // BigInt(1)) + mult
    elseif x isa Pair && x.first isa Symbol &&
           (x.second isa Integer || x.second isa Rational)
        _graph_kf_add!(f, x.first, mult * Rational{BigInt}(x.second))
    elseif x isa AbstractVector || x isa Tuple
        for el in x
            _graph_kf_add!(f, el, mult)
        end
    else
        _graph_refuse(:BadInput, "msg" => "unrecognized kinematic entry",
                      "entry" => string(x))
    end
    return f
end
_graph_kinform(x) = _graph_kf_add!(_GraphKinForm(), x, 1 // BigInt(1))

_graph_kf_symbols(f::_GraphKinForm) =
    Symbol[k for k in keys(f) if k != _GRAPH_CONSTKEY && !iszero(f[k])]

"""Linear form -> element of ring R, using kinmap :: Dict{Symbol,gen}."""
function _graph_kf_poly(f::_GraphKinForm, R::QQMPolyRing,
                        kinmap::Dict{Symbol,QQMPolyRingElem})
    p = zero(R)
    for (k, c) in f
        iszero(c) && continue
        if k == _GRAPH_CONSTKEY
            p += R(QQ(c))
        else
            haskey(kinmap, k) || _graph_refuse(:BadInput,
                "msg" => "kinematic symbol not in kinvars", "symbol" => string(k))
            p += QQ(c) * kinmap[k]
        end
    end
    return p
end

# ---------------------------------------------------------------------------
# Small exact combinatorics (no external deps).
# ---------------------------------------------------------------------------
"""All k-subsets of 1:n, lexicographic (ascending index vectors)."""
function _graph_subsets(n::Int, k::Int)
    (k < 0 || k > n) && return Vector{Int}[]
    k == 0 && return [Int[]]
    out = Vector{Int}[]
    idx = collect(1:k)
    while true
        push!(out, copy(idx))
        j = k
        while j >= 1 && idx[j] == n - k + j
            j -= 1
        end
        j == 0 && break
        idx[j] += 1
        for l in (j+1):k
            idx[l] = idx[l-1] + 1
        end
    end
    return out
end

"""Union-find over `subset` of edges. Returns nothing if the subset contains
a cycle (incl. self-loops), else the component root id per vertex 1..nv.
(Same union-find pattern as pld_complete_c3.jl comps(), lines 103-113.)"""
function _graph_forest_components(nv::Int, pairs::Vector{Tuple{Int,Int}},
                                  subset::Vector{Int})
    parent = collect(1:nv)
    find(a) = begin
        while parent[a] != a
            parent[a] = parent[parent[a]]
            a = parent[a]
        end
        a
    end
    for e in subset
        (a, b) = pairs[e]
        a == b && return nothing
        ra, rb = find(a), find(b)
        ra == rb && return nothing
        parent[ra] = rb
    end
    return Int[find(v) for v in 1:nv]
end

"""Number of connected components of the whole graph (cycles allowed)."""
function _graph_ncomponents(nv::Int, pairs::Vector{Tuple{Int,Int}})
    parent = collect(1:nv)
    find(a) = begin
        while parent[a] != a
            parent[a] = parent[parent[a]]
            a = parent[a]
        end
        a
    end
    for (a, b) in pairs
        a == b && continue
        ra, rb = find(a), find(b)
        ra != rb && (parent[ra] = rb)
    end
    return length(unique(Int[find(v) for v in 1:nv]))
end

# ---------------------------------------------------------------------------
# U — first Symanzik.
# Primary route: direct spanning-tree sum  U = Σ_{spanning trees T} ∏_{e∉T} x_e
# (pld_complete_c3.jl:33-35 states this equivalence; pure ring ops, exact).
# Cross-check route: ∏x_e · det(reduced weighted Laplacian, conductance
# 1/x_e) — the weighted MATRIX-TREE theorem, transliterated from
# pld_complete_c3.jl build_UF lines 36-51 (fraction-field det, divexact).
# ---------------------------------------------------------------------------
function _graph_U_treesum(R::QQMPolyRing, xs::Vector{QQMPolyRingElem},
                          pairs::Vector{Tuple{Int,Int}}, nv::Int)
    ne = length(pairs)
    U = zero(R)
    for T in _graph_subsets(ne, nv - 1)
        _graph_forest_components(nv, pairs, T) === nothing && continue
        # nv-1 acyclic edges ⇒ exactly one component ⇒ spanning tree
        cut = one(R)
        inT = falses(ne)
        for e in T
            inT[e] = true
        end
        for e in 1:ne
            inT[e] || (cut *= xs[e])
        end
        U += cut
    end
    return U
end

function _graph_U_laplacian(R::QQMPolyRing, xs::Vector{QQMPolyRingElem},
                            pairs::Vector{Tuple{Int,Int}}, nv::Int)
    FR = fraction_field(R)                       # pld_complete_c3.jl:38
    L = zero_matrix(FR, nv, nv)
    for (e, (a, b)) in enumerate(pairs)
        a == b && continue                       # self-loop: no conductance
        w = FR(one(R)) // FR(xs[e])              # conductance 1/x_e (pld:40)
        L[a, a] += w; L[b, b] += w
        L[a, b] -= w; L[b, a] -= w
    end
    Lr = L[1:nv-1, 1:nv-1]                       # delete row/col nv (pld:46)
    U_frac = prod(FR.(xs)) * det(Lr)             # pld:48
    return divexact(numerator(U_frac), denominator(U_frac))   # pld:49-51
end

# ---------------------------------------------------------------------------
# F — second Symanzik (Euclidean-positive arrangement, i.e. the −F_paper /
# pySecDec-F sign: scripts/make_fixtures.py convention notes; SubTropica.wl
# Feynman branch wl:3001-3005, wl:3591 `Fphys = -F`):
#   F = Σ_{spanning 2-forests (T1,T2)} (∏_{e∉T1∪T2} x_e)·(−p(T1)²)
#       + U · Σ_e m_e² x_e
# transliterated from pld_complete_c3.jl two_forest_F (lines 82-132) +
# the mass-term assembly F = F0 + mass·U (line 68, generalized to all edges).
# p(T1)² is evaluated from the user dot-product table with the last leg
# eliminated by momentum conservation (pld psq(), lines 86-97): the value is
# component-choice independent since (Σ_{i∈S}p_i)² = (Σ_{i∉S}p_i)².
# ---------------------------------------------------------------------------
function _graph_psq(S::Set{Int}, nlegs::Int,
                    kd::Dict{Tuple{Int,Int},_GraphKinForm},
                    R::QQMPolyRing, kinmap::Dict{Symbol,QQMPolyRingElem})
    nlegs <= 1 && return zero(R)                 # 0/1 legs ⇒ zero inflow
    c = Int[(i in S ? 1 : 0) - (nlegs in S ? 1 : 0) for i in 1:(nlegs-1)]
    p = zero(R)
    for i in 1:(nlegs-1), j in i:(nlegs-1)
        w = i == j ? c[i]^2 : 2 * c[i] * c[j]
        w == 0 && continue
        haskey(kd, (i, j)) || _graph_refuse(:MissingDotProduct,
            "msg" => "kinematics must supply p_i.p_j for all 1<=i<=j<=nlegs-1 " *
                     "(last leg eliminated by momentum conservation)",
            "pair" => [i, j])
        p += w * _graph_kf_poly(kd[(i, j)], R, kinmap)
    end
    return p
end

function _graph_F0(R::QQMPolyRing, xs::Vector{QQMPolyRingElem},
                   pairs::Vector{Tuple{Int,Int}}, nv::Int,
                   legs_at::Vector{Int}, nlegs::Int,
                   kd::Dict{Tuple{Int,Int},_GraphKinForm},
                   kinmap::Dict{Symbol,QQMPolyRingElem})
    ne = length(pairs)
    k = nv - 2
    k < 0 && return zero(R)
    F0 = zero(R)
    for sub in _graph_subsets(ne, k)
        comp = _graph_forest_components(nv, pairs, sub)
        comp === nothing && continue
        roots = unique(comp)
        length(roots) == 2 || continue           # spanning 2-forest (pld:99-119)
        S1 = Set{Int}(i for i in 1:nlegs if comp[legs_at[i]] == roots[1])
        coefpoly = -_graph_psq(S1, nlegs, kd, R, kinmap)   # (pld:127)
        iszero(coefpoly) && continue
        insub = falses(ne)
        for e in sub
            insub[e] = true
        end
        cut = one(R)
        for e in 1:ne
            insub[e] || (cut *= xs[e])
        end
        F0 += coefpoly * cut                     # (pld:128-129)
    end
    return F0
end

# ---------------------------------------------------------------------------
# symanzik_UF — the raw graph->U,F builder (UNGAUGED, all edges).
# ---------------------------------------------------------------------------
"""
    symanzik_UF(edges, nodes, masses, kinematics;
                kinvars=nothing, xnames=nothing, crosscheck=true)

Build the Symanzik polynomials of a Feynman graph natively (weighted
matrix-tree + spanning-2-forest; pattern of pld_complete_c3.jl, cited at the
implementation sites). Arguments:
- `edges`   : internal edges as (v,w) vertex pairs (parallel edges and
              self-loops allowed).
- `nodes`   : external-leg attachment vertices; leg i carries momentum p_i,
              with p_nlegs = −Σ_{i<nlegs} p_i (momentum conservation).
- `masses`  : per-edge SQUARED masses (0 / Rational / Symbol / linear form).
- `kinematics` : dict-like (i,j) => linear form for p_i·p_j,
              1 ≤ i ≤ j ≤ nlegs−1 ONLY (see _graph_psq).
- `kinvars` : pinned kinvar ORDER (types.jl ring order is vars-then-kinvars);
              default = all used symbols, sorted by name.
- `xnames`  : Schwinger-parameter names; default `x0..x{E-1}` (0-based
              pySecDec/fixture naming — fixtures/README.md).
- `crosscheck` : also compute U by the reduced-Laplacian determinant and
              refuse :InternalCheckFailed on mismatch.

Returns a NamedTuple (U, F, ring, xsyms, kinsyms, xgens, kinmap, L, nv, ne).
F is the Euclidean-positive arrangement (see _graph_F0 header). Refusals:
:Disconnected (U ≡ 0 ⇔ disconnected graph — the graph-side analogue of
STSymanzik::zeroU, wl:3101/3510-3516), :NoLoop, :Scaleless (F ≡ 0,
STSymanzik::zeroF wl:3527; the vacuumPeriod F=1 shortcut wl:3519-3526 is NOT
ported — refuse loudly), :TooLarge, :BadInput, :MissingDotProduct.
"""
function symanzik_UF(edges, nodes, masses, kinematics;
                     kinvars::Union{Nothing,Vector{Symbol}}=nothing,
                     xnames::Union{Nothing,Vector{String}}=nothing,
                     crosscheck::Bool=true)
    # --- normalize + validate the graph -------------------------------------
    ne = length(edges)
    ne >= 1 || _graph_refuse(:BadInput, "msg" => "no internal edges")
    length(masses) == ne || _graph_refuse(:BadInput,
        "msg" => "masses length != number of edges",
        "n_edges" => ne, "n_masses" => length(masses))
    rawpairs = Tuple{Int,Int}[]
    for e in edges
        (length(e) == 2 && all(x -> x isa Integer, collect(e))) ||
            _graph_refuse(:BadInput, "msg" => "edge is not a vertex pair",
                          "edge" => string(e))
        push!(rawpairs, (Int(e[1]), Int(e[2])))
    end
    legs = Int[Int(v) for v in nodes]
    nlegs = length(legs)
    verts = sort(unique(vcat([v for p in rawpairs for v in p], legs)))
    vidx = Dict{Int,Int}(v => i for (i, v) in enumerate(verts))
    nv = length(verts)
    pairs = Tuple{Int,Int}[(vidx[a], vidx[b]) for (a, b) in rawpairs]
    legs_at = Int[vidx[v] for v in legs]

    ncomp = _graph_ncomponents(nv, pairs)
    ncomp == 1 || _graph_refuse(:Disconnected,
        "msg" => "graph (edges + leg attachment vertices) is disconnected; " *
                 "U would vanish identically (STSymanzik::zeroU analogue, wl:3101)",
        "n_components" => ncomp)
    L = ne - nv + 1
    L >= 1 || _graph_refuse(:NoLoop,
        "msg" => "tree-level graph (L = E - V + 1 = 0): not a loop integral",
        "E" => ne, "V" => nv)

    # loud size guard on the exponential subset scans (SPEC_B1 §3.2 style)
    work = binomial(BigInt(ne), BigInt(nv - 1)) +
           binomial(BigInt(ne), BigInt(max(nv - 2, 0)))
    work <= 500_000 || _graph_refuse(:TooLarge,
        "msg" => "spanning-tree/2-forest enumeration too large",
        "subsets" => string(work))

    # --- kinematic table + kinvar collection --------------------------------
    kd = Dict{Tuple{Int,Int},_GraphKinForm}()
    for (key, val) in _graph_pairs_iterable(kinematics)
        (length(key) == 2 && all(x -> x isa Integer, collect(key))) ||
            _graph_refuse(:BadInput, "msg" => "kinematics key is not (i,j)",
                          "key" => string(key))
        i, j = minmax(Int(key[1]), Int(key[2]))
        (1 <= i && j <= nlegs - 1) || _graph_refuse(:BadInput,
            "msg" => "dot products are declared for legs 1..nlegs-1 only " *
                     "(leg nlegs is eliminated by momentum conservation)",
            "pair" => [i, j], "nlegs" => nlegs)
        haskey(kd, (i, j)) && _graph_refuse(:BadInput,
            "msg" => "duplicate dot product", "pair" => [i, j])
        kd[(i, j)] = _graph_kinform(val)
    end
    massforms = _GraphKinForm[_graph_kinform(m) for m in masses]

    used = Set{Symbol}()
    for f in values(kd)
        union!(used, _graph_kf_symbols(f))
    end
    for f in massforms
        union!(used, _graph_kf_symbols(f))
    end
    kv = kinvars === nothing ? sort(collect(used); by=string) : copy(kinvars)
    issubset(used, Set(kv)) || _graph_refuse(:BadInput,
        "msg" => "kinvars does not cover all symbols used",
        "missing" => string.(sort(collect(setdiff(used, Set(kv))); by=string)))

    xn = xnames === nothing ? String["x$(e-1)" for e in 1:ne] : copy(xnames)
    length(xn) == ne || _graph_refuse(:BadInput, "msg" => "xnames length != E")
    names = vcat(xn, String[string(s) for s in kv])
    length(unique(names)) == length(names) || _graph_refuse(:BadInput,
        "msg" => "variable name collision (vars/kinvars)", "names" => names)

    # ring order = integration vars first, then kinvars (types.jl contract)
    R, g = polynomial_ring(QQ, names)
    xs = QQMPolyRingElem[g[i] for i in 1:ne]
    kinmap = Dict{Symbol,QQMPolyRingElem}(kv[i] => g[ne+i] for i in 1:length(kv))

    # --- U (spanning trees; matrix-tree cross-check) ------------------------
    U = _graph_U_treesum(R, xs, pairs, nv)
    iszero(U) && _graph_refuse(:Disconnected,
        "msg" => "U vanished identically (no spanning tree)")
    if crosscheck
        Ulap = _graph_U_laplacian(R, xs, pairs, nv)
        U == Ulap || _graph_refuse(:InternalCheckFailed,
            "msg" => "spanning-tree U != matrix-tree (Laplacian) U — " *
                     "internal bug, do not proceed")
    end

    # --- F = F0(2-forests) + U·Σ m_e² x_e ------------------------------------
    F = _graph_F0(R, xs, pairs, nv, legs_at, nlegs, kd, kinmap)
    for e in 1:ne
        mp = _graph_kf_poly(massforms[e], R, kinmap)
        iszero(mp) || (F += U * mp * xs[e])      # pld_complete_c3.jl:68 pattern
    end
    iszero(F) && _graph_refuse(:Scaleless,
        "msg" => "F (second Symanzik) is identically zero — scaleless " *
                 "(STSymanzik::zeroF, wl:3527; vacuumPeriod shortcut " *
                 "wl:3519-3526 deliberately NOT ported)")

    return (U=U, F=F, ring=R, xsyms=Symbol.(xn), kinsyms=kv, xgens=xs,
            kinmap=kinmap, L=L, nv=nv, ne=ne)
end

# kinematics may arrive as a Dict or as a Vector of Pairs
_graph_pairs_iterable(d::AbstractDict) = collect(d)
_graph_pairs_iterable(d::AbstractVector) = [(p.first, p.second) for p in d]
_graph_pairs_iterable(d) = _graph_refuse(:BadInput,
    "msg" => "kinematics must be Dict-like or a Vector of Pairs")

# ---------------------------------------------------------------------------
# graph_to_quadruple — graph -> EulerIntegrand (the raw Euler quadruple,
# types.jl / DESIGN.md RT-x). Two variants:
#  * Feynman/Schwinger (default): paper eq 1.4 with GL(1) fixed Cheng-Wu
#    style on ONE Schwinger parameter (kwarg `gauge`, default = the LAST
#    kept edge — the fixture convention, "Cheng-Wu x3=1 (last Feynman
#    parameter)" / the 3-loop fixture's "Cheng-Wu x7=1"):
#      polys = [(U, eU), (F, eF)],
#      eU = ν_tot − (L+1)(D−2ε)/2 = EpsExp(ν_tot − (L+1)D/2, L+1)  (wl:2997)
#      eF = −(ν_tot − L(D−2ε)/2)  = EpsExp(LD/2 − ν_tot, −L)       (wl:2996)
#      prefactor c = (−1)^{ν_tot} (eq 1.4 ∏(−1)^{ν_e}),
#      gammas = Γ(−eF) · ∏_{ν_e≠1} Γ(ν_e)^{−1}                     (wl:3010)
#  * Lee-Pomeransky (`lee_pomeransky=true`; arXiv:1308.6676): NO gauge,
#      polys = [(U+F, EpsExp(−D/2, 1))],
#      gammas = Γ(D/2−ε…) = (EpsExp(D/2,−1),1), (EpsExp((L+1)D/2−ν_tot,−(L+1)),−1),
#      ∏_{ν_e≠1} Γ(ν_e)^{−1}; same c = (−1)^{ν_tot}.
# Monomial exponents are stored FLAT-measure: nu_i = ν_e − 1
# (fixtures/README.md schema note; paper eq 1.4 ∏ x^{ν_e−1}).
# Pinched edges (ν_e = 0) are dropped from the graph but keep the 0-based
# x-naming of their original edge index (STGetIntegrandData pinched
# semantics, wl:2977-2979 + xvars naming wl:2984). Negative ν_e = numerator
# territory ⇒ typed refusal (this pass; STGetIntegrandDataNumerators
# wl:3014-3080 not ported).
# ---------------------------------------------------------------------------
"""
    graph_to_quadruple(edges, nodes, masses, kinematics;
                       nu=ones(Int,length(edges)), D=4, gauge=nothing,
                       lee_pomeransky=false, gammaE_norm=false,
                       kinvars=nothing, numerators=nothing) -> EulerIntegrand

Graph front-end (DESIGN.md C2 spec pulled forward; see file header for the
full citation map). `nu` = integer propagator powers (0 = pinched, dropped;
negative REFUSES). `D` = integer spacetime dimension, always continued as
D−2ε (STSymanzik default Dimension -> 4 − 2 eps, wl:3104-3106). `gauge` =
1-based ORIGINAL edge index whose Schwinger parameter is set to 1 (Feynman
variant only; default last kept edge). `gammaE_norm=true` adds the
e^{L·γ_E·ε} normalization of paper eq 1.3 (default OFF = frozen-fixture /
pySecDec-disteval convention). `numerators` MUST be nothing/empty
(:NumeratorsNotSupported refusal — CenterDot machinery wl:3546-3568/3606-3660
deliberately not ported this pass).

Declared-Euclidean contract (types.jl): every coefficient of U and F must be
POSITIVE, else :NonEuclidean refusal (parity with scripts/make_fixtures.py
mma_poly require_positive).
"""
function graph_to_quadruple(edges, nodes, masses, kinematics;
                            nu::AbstractVector{<:Integer}=ones(Int, length(edges)),
                            D::Integer=4,
                            gauge::Union{Nothing,Integer}=nothing,
                            lee_pomeransky::Bool=false,
                            gammaE_norm::Bool=false,
                            kinvars::Union{Nothing,Vector{Symbol}}=nothing,
                            numerators=nothing)
    if !(numerators === nothing ||
         (numerators isa AbstractVector && isempty(numerators)))
        _graph_refuse(:NumeratorsNotSupported,
            "msg" => "numerator/CenterDot handling (STSymanzik step 5 " *
                     "wl:3546-3568, step 7 wl:3606-3660; " *
                     "STGetIntegrandDataNumerators wl:3014-3080) is NOT " *
                     "ported in this pass",
            "numerators" => string(numerators))
    end
    ne = length(edges)
    length(nu) == ne || _graph_refuse(:BadInput,
        "msg" => "nu length != number of edges (upstream guard, wl:3240-3243)")
    any(v -> v < 0, nu) && _graph_refuse(:NumeratorExponent,
        "msg" => "negative propagator exponent = numerator (wl:3020-3025 " *
                 "semantics); not supported this pass",
        "nu" => collect(Int, nu))
    D >= 1 || _graph_refuse(:BadInput, "msg" => "D must be a positive integer")

    keep = Int[e for e in 1:ne if nu[e] > 0]     # pinch drop (wl:2977-2979)
    isempty(keep) && _graph_refuse(:BadInput,
        "msg" => "no propagators left after pinching (wl:3287-3290 guard)")

    kedges = [edges[e] for e in keep]
    kmasses = [masses[e] for e in keep]
    xn = String["x$(e-1)" for e in keep]         # 0-based fixture naming
    UF = symanzik_UF(kedges, nodes, kmasses, kinematics;
                     kinvars=kinvars, xnames=xn)

    # declared-Euclidean positivity (types.jl region contract)
    for (P, what) in ((UF.U, "U"), (UF.F, "F"))
        for i in 1:length(P)
            coeff(P, i) > 0 || _graph_refuse(:NonEuclidean,
                "msg" => "non-positive coefficient in $what — violates the " *
                         "declared-Euclidean positivity contract " *
                         "(types.jl EulerIntegrand; make_fixtures.py mma_poly)",
                "poly" => what, "coeff" => string(coeff(P, i)),
                "monomial" => string(exponent_vector(P, i)))
        end
    end

    L = UF.L
    nutot = sum(Int, nu)                          # ν_e = 0 contribute nothing
    D2 = Rational{BigInt}(D) // 2
    gammadenoms = Tuple{EpsExp,Int}[(EpsExp(Rational{BigInt}(nu[e])), -1)
                                    for e in keep if nu[e] != 1]
    c = Rational{BigInt}(iseven(nutot) ? 1 : -1)  # eq 1.4 ∏(−1)^{ν_e}
    gE = gammaE_norm ? Rational{BigInt}(L) : 0 // BigInt(1)   # eq 1.3 e^{εLγE}

    if lee_pomeransky
        gauge === nothing || _graph_refuse(:BadInput,
            "msg" => "lee_pomeransky variant is projectively complete — " *
                     "no GL(1) gauge fixing; pass gauge=nothing")
        G = UF.U + UF.F                           # pld_complete_c3.jl:148
        polys = Tuple{QQMPolyRingElem,EpsExp}[
            (G, EpsExp(-D2, Rational{BigInt}(1)))]
        gammas = vcat(Tuple{EpsExp,Int}[
                (EpsExp(D2, Rational{BigInt}(-1)), 1),
                (EpsExp((L + 1) * D2 - nutot, Rational{BigInt}(-(L + 1))), -1)],
            gammadenoms)
        pref = Prefactor(c, 0, gammas, gE)
        nuflat = EpsExp[EpsExp(Rational{BigInt}(nu[e] - 1)) for e in keep]
        return EulerIntegrand(pref, nuflat, polys, UF.xsyms, UF.kinsyms, UF.ring)
    end

    # ---- Feynman/Schwinger variant (eq 1.4, GL(1) fixed) --------------------
    g = gauge === nothing ? keep[end] : Int(gauge)
    g in keep || _graph_refuse(:BadInput,
        "msg" => "gauge must be a kept (non-pinched) 1-based edge index",
        "gauge" => g, "kept" => keep)
    length(keep) >= 2 || _graph_refuse(:BadInput,
        "msg" => "need >=2 propagators for a gauged Schwinger rep " *
                 "(one parameter is set to 1)")

    gpos = findfirst(==(g), keep)                 # position among kept edges
    keptnames = String[xn[i] for i in 1:length(keep) if i != gpos]
    rnames = vcat(keptnames, String[string(s) for s in UF.kinsyms])
    Rq, gq = polynomial_ring(QQ, rnames)
    nq = length(keptnames)
    imgs = Vector{QQMPolyRingElem}(undef, length(xn) + length(UF.kinsyms))
    j = 0
    for i in 1:length(xn)
        if i == gpos
            imgs[i] = one(Rq)                     # Cheng-Wu x_gauge -> 1
        else
            j += 1
            imgs[i] = gq[j]
        end
    end
    for i in 1:length(UF.kinsyms)
        imgs[length(xn)+i] = gq[nq+i]
    end
    U1 = evaluate(UF.U, imgs)
    F1 = evaluate(UF.F, imgs)

    eU = EpsExp(Rational{BigInt}(nutot) - (L + 1) * D2,
                Rational{BigInt}(L + 1))          # wl:2997
    eF = EpsExp(L * D2 - nutot, Rational{BigInt}(-L))   # wl:2996
    polys = Tuple{QQMPolyRingElem,EpsExp}[(U1, eU), (F1, eF)]
    gammas = vcat(Tuple{EpsExp,Int}[
            (EpsExp(Rational{BigInt}(nutot) - L * D2, Rational{BigInt}(L)), 1)],
        gammadenoms)                              # Γ(−eF)/∏Γ(ν_e), wl:3010
    pref = Prefactor(c, 0, gammas, gE)
    nuflat = EpsExp[EpsExp(Rational{BigInt}(nu[e] - 1)) for e in keep if e != g]
    vars = Symbol[Symbol(s) for s in keptnames]
    return EulerIntegrand(pref, nuflat, polys, vars, UF.kinsyms, Rq)
end

# ---------------------------------------------------------------------------
# JSON surface — schema "subtropica-euler-quad-v1" (fixtures/README.md, matches
# types.jl EulerIntegrand field-for-field; reference generator:
# scripts/make_fixtures.py). The poly-string emitter replicates
# make_fixtures.py mma_poly EXACTLY (explicit `*`, `^k` for k>=2,
# content-canonical rationals, terms sorted by exponent tuple DESCENDING
# over vars-then-kinvars) so that graph-front-end output is byte-comparable
# to the frozen Phase-A fixtures. The loader parses the same restricted
# positive-sum grammar (subset of CONTRACTS.md §d).
# ---------------------------------------------------------------------------
function _graph_qstr(q::Rational{BigInt})
    isone(denominator(q)) ? string(numerator(q)) :
        string(numerator(q)) * "/" * string(denominator(q))
end
_graph_qstr(q::Rational) = _graph_qstr(Rational{BigInt}(q))
function _graph_qstr(q::QQFieldElem)
    isone(denominator(q)) ? string(numerator(q)) :
        string(numerator(q)) * "/" * string(denominator(q))
end

_graph_epsjson(e::EpsExp) = Any[_graph_qstr(e.a), _graph_qstr(e.b)]

function _graph_qparse(s::AbstractString)::Rational{BigInt}
    t = strip(s)
    if occursin('/', t)
        pieces = split(t, '/')
        length(pieces) == 2 || _graph_refuse(:BadInput,
            "msg" => "bad rational string", "s" => String(t))
        return parse(BigInt, strip(pieces[1])) // parse(BigInt, strip(pieces[2]))
    end
    return parse(BigInt, t) // BigInt(1)
end

"""Lexicographic (elementwise) strict less-than on Int vectors."""
function _graph_lexless(a::Vector{Int}, b::Vector{Int})
    for i in 1:min(length(a), length(b))
        a[i] != b[i] && return a[i] < b[i]
    end
    return length(a) < length(b)
end

"""Poly -> canonical fixture string (make_fixtures.py mma_poly parity:
exponent-tuple-DESCENDING term order, require-positive coefficients)."""
function _graph_poly_string(p::QQMPolyRingElem)
    iszero(p) && return "0"
    names = String[string(s) for s in symbols(parent(p))]
    ts = Tuple{Vector{Int},QQFieldElem}[(exponent_vector(p, i), coeff(p, i))
                                        for i in 1:length(p)]
    sort!(ts; lt=(a, b) -> _graph_lexless(a[1], b[1]), rev=true)
    pieces = String[]
    for (ev, c) in ts
        c > 0 || _graph_refuse(:NonEuclidean,
            "msg" => "non-positive coefficient in emitted poly string " *
                     "(subtropica-euler-quad-v1 positivity contract)",
            "coeff" => string(c))
        parts = String[]
        for (vn, k) in zip(names, ev)
            if k == 1
                push!(parts, vn)
            elseif k >= 2
                push!(parts, vn * "^" * string(k))
            end
        end
        (isempty(parts) || !isone(c)) && pushfirst!(parts, _graph_qstr(c))
        push!(pieces, join(parts, "*"))
    end
    return join(pieces, " + ")
end

"""Parse the restricted positive-sum poly grammar into ring R."""
function _graph_parse_poly(s::AbstractString, R::QQMPolyRing)
    nameidx = Dict{String,Int}(string(v) => i for (i, v) in enumerate(symbols(R)))
    n = nvars(R)
    t = strip(s)
    t == "0" && return zero(R)
    isempty(t) && _graph_refuse(:BadInput, "msg" => "empty poly string")
    occursin('-', t) && _graph_refuse(:BadInput,
        "msg" => "poly grammar is positive-sum only (subtropica-euler-quad-v1)",
        "poly" => String(t))
    B = MPolyBuildCtx(R)
    for term in split(t, '+')
        c = Rational{BigInt}(1)
        ev = zeros(Int, n)
        facs = split(strip(term), '*')
        (isempty(facs) || any(f -> isempty(strip(f)), facs)) &&
            _graph_refuse(:BadInput, "msg" => "bad term", "term" => String(term))
        for f0 in facs
            f = strip(f0)
            if occursin(r"^\d+(/\d+)?$", f)
                c *= _graph_qparse(f)
            else
                m = match(r"^([A-Za-z][A-Za-z0-9_]*)(\^(\d+))?$", f)
                m === nothing && _graph_refuse(:BadInput,
                    "msg" => "unparseable factor", "factor" => String(f))
                haskey(nameidx, m.captures[1]) || _graph_refuse(:BadInput,
                    "msg" => "unknown variable in poly string",
                    "var" => String(m.captures[1]))
                ev[nameidx[m.captures[1]]] +=
                    m.captures[3] === nothing ? 1 : parse(Int, m.captures[3])
            end
        end
        push_term!(B, QQ(c), ev)
    end
    return finish(B)
end

"""
    quadruple_json(E::EulerIntegrand; id, title="", divergent=false,
                   notes="", provenance=Dict()) -> Dict{String,Any}

Emit the raw Euler quadruple as a schema subtropica-euler-quad-v1 dict
(fixtures/README.md). Values policy (DESIGN.md RT-iii): NO result values —
this emits the INTEGRAND only.
"""
function quadruple_json(E::EulerIntegrand; id::AbstractString,
                        title::AbstractString="", divergent::Bool=false,
                        notes::AbstractString="",
                        provenance::AbstractDict=Dict{String,Any}())
    pf = E.prefactor
    return Dict{String,Any}(
        "schema" => "subtropica-euler-quad-v1",
        "id" => String(id), "title" => String(title),
        "divergent" => divergent,
        "prefactor" => Dict{String,Any}(
            "c" => _graph_qstr(pf.c),
            "eps_power" => pf.eps_power,
            "gammas" => Any[Dict{String,Any}("arg" => _graph_epsjson(g),
                                             "power" => p)
                            for (g, p) in pf.gammas],
            "gammaE_eps" => _graph_qstr(pf.gammaE_eps)),
        "nu" => Any[_graph_epsjson(nui) for nui in E.nu],
        "polys" => Any[Dict{String,Any}("poly" => _graph_poly_string(P),
                                        "exp" => _graph_epsjson(e))
                       for (P, e) in E.polys],
        "vars" => String[string(v) for v in E.vars],
        "kinvars" => String[string(v) for v in E.kinvars],
        "notes" => String(notes),
        "provenance" => Dict{String,Any}(provenance))
end

"""Write a quadruple dict as JSON (indent 1, trailing newline — the
make_fixtures.py write_json format). Returns the path."""
function write_quadruple_json(path::AbstractString, d::AbstractDict)
    open(path, "w") do io
        write(io, JSON.json(d, 1))
        write(io, "\n")
    end
    return path
end

"""
    load_quadruple(path_or_dict) -> (E = EulerIntegrand, raw = Dict)

Load a subtropica-euler-quad-v1 fixture into an EulerIntegrand (ring =
QQ[vars..., kinvars...], the types.jl order; Nemo ring caching makes polys
comparable across independent loads of the same var list). The `divergent`
flag is surfaced in `raw` — loading a divergent fixture is legal (the
pipeline must refuse it LOUDLY downstream, Phase-A gate iv); comparing
values against it is verify.jl's refusal, not ours.
"""
function load_quadruple(d::AbstractDict)
    get(d, "schema", "") == "subtropica-euler-quad-v1" || _graph_refuse(:BadInput,
        "msg" => "not a subtropica-euler-quad-v1 object",
        "schema" => string(get(d, "schema", "")))
    vars = Symbol[Symbol(String(v)) for v in d["vars"]]
    kins = Symbol[Symbol(String(v)) for v in d["kinvars"]]
    R, _ = polynomial_ring(QQ, vcat(String[string(v) for v in vars],
                                    String[string(v) for v in kins]))
    nu = EpsExp[EpsExp(_graph_qparse(p[1]), _graph_qparse(p[2])) for p in d["nu"]]
    polys = Tuple{QQMPolyRingElem,EpsExp}[
        (_graph_parse_poly(pj["poly"], R),
         EpsExp(_graph_qparse(pj["exp"][1]), _graph_qparse(pj["exp"][2])))
        for pj in d["polys"]]
    pf = d["prefactor"]
    gammas = Tuple{EpsExp,Int}[
        (EpsExp(_graph_qparse(g["arg"][1]), _graph_qparse(g["arg"][2])),
         Int(g["power"])) for g in pf["gammas"]]
    pref = Prefactor(_graph_qparse(pf["c"]), Int(pf["eps_power"]), gammas,
                     _graph_qparse(pf["gammaE_eps"]))
    E = EulerIntegrand(pref, nu, polys, vars, kins, R)
    return (E=E, raw=Dict{String,Any}(d))
end
load_quadruple(path::AbstractString) = load_quadruple(JSON.parsefile(path))
