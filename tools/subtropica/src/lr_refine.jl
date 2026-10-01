# SPDX-License-Identifier: GPL-3.0-or-later
#
# Derived from HyperInt (Erik Panzer, "Algorithms for the symbolic integration of
# hyperlogarithms with applications to Feynman integrals", Comput. Phys. Commun. 188
# (2015) 148, arXiv:1403.3385; https://bitbucket.org/PanzerErik/hyperint),
# Copyright (C) 2014 Erik Panzer, licensed under the GNU General Public License
# v3.0 or later. This Julia transliteration of HyperInt's compatibility-graph
# reduction procedures (with F. Brown, arXiv:0910.0114, as the mathematical source)
# is Copyright (c) 2026 Anthropic, PBC; created by Matthew D. Schwartz, code written
# by Claude (Anthropic) under his supervision, and is distributed under the same
# GNU General Public License v3.0 or later. The rest of this package is MIT-licensed;
# see ../NOTICE.
#
# lr_refine.jl — compatibility-graph LR refinement.
# ---------------------------------------------------------------------------
# WHY THIS FILE EXISTS
#   hyperflint's find_lr_orders computes the SIMPLE ("Fubini_Lungo") polynomial
#   reduction: an UPPER bound on the letter set of an iterated hyperlogarithm
#   integration (SubTropica paper Sec 3.4.1 / App A.7; STFubiniLR at
#   SubTropica.wl:16524-16593 — leading coeff, constant term, discriminant,
#   ALL pairwise resultants, factored and intersected over elimination paths).
#   Because it is an upper bound, a "blocking letter" reported by that engine
#   may be SPURIOUS.  This file implements the sharper compatibility-graph
#   reduction so that a letter present in the simple bound but ABSENT from the
#   compatibility-graph bound is CERTIFIED spurious (both bounds are provable
#   upper bounds; the refined one is contained in the naive one).
#
# LAW / PROVENANCE (semantics transliterated, cited per rule below):
#   [B]  F. Brown, "On the periods of some Feynman integrals", arXiv:0910.0114,
#        Sec. 6:  eq (65) ([0,f]_x = const coeff, [inf,f]_x = leading coeff,
#        D_x(f) = lc^-1 [f,f']_x), Definition 68 (simple reduction),
#        Remark 71 (Fubini intersection S_[I] = ∩_e (S_[I\e])_e),
#        Definition 73 (compatibility-graph single reduction),
#        Definition 74 (Fubini + compatibility intersection, initial C =
#        complete graph), eq (71) (the exactly-three excluded spurious
#        resultant shapes = "four distinct grandparents").
#   [P]  E. Panzer, "Algorithms for the symbolic integration of
#        hyperlogarithms...", arXiv:1403.3385, Sec. 4.2 (cgReduction interface,
#        L[I] = [S_I, C_I]) and Example 4.2 (worked {(1+x)^2+y, y+z^2} case —
#        reproduced verbatim as a fixture in test_lr_refine.jl).
#   [H]  HyperInt.mpl (Panzer, bitbucket PanzerErik/hyperint, master):
#        cgSingleReductionFrancis (HyperInt.mpl:2913-2986; the DEFAULT,
#        HyperInt.mpl:79) and cgSingleReductionClique (HyperInt.mpl:2846-2911),
#        irreducibles (HyperInt.mpl:2614-2628, monomial/constant dropping),
#        cgReduction driver (HyperInt.mpl:2989-3086: parent gating by
#        max degree, lowest-degree-first, pairwise (S,C) intersection).
#   [ST] SubTropica.wl:15659-15691 (STpreparePolysAndPairsMaple): the initial
#        compatibility graph of a SUM of integrand terms is the union of the
#        per-term complete graphs (pairs = Subsets[polys,{2}] per term) — the
#        multi-group init used here; SubTropica.wl:16531-16556 (STFubiniLR
#        term1: lead + const + disc — matches [B] Def 68 singletons).
#
# SCOPE PIN (this pass): rational letters over Q(kinematics) ONLY.  Everything
#   is Nemo-native exact arithmetic in one QQMPolyRing whose variables are the
#   integration variables plus the kinematic symbols (EulerIntegrand ring
#   convention, src/types.jl).  Algebraic letters (types.jl AlgLetter /
#   HAlgLetter, sqrt/Root expressions) are REFUSED, typed
#   (LRRefineRefusal(:AlgebraicLetter)) — the FindRoots/Wm/Wp tier is the
#   engine cascade's job, not this file's.
#
# SOUNDNESS DIRECTION (what a verdict means):
#   * letter ABSENT from the refined bound  => CERTIFIED SPURIOUS (theorem-
#     grade: Brown Thm 77-class upper-bound property of the cg reduction).
#   * letter PRESENT in the refined bound   => "survives"; the bound is still
#     an over-approximation, so survival is necessary-but-not-sufficient
#     evidence that the letter/obstruction is real.  is_blocking_real == true
#     means exactly "not certifiable spurious by this method".
#
# OWNERSHIP: this module owns ONLY this file and test/test_lr_refine.jl
#   (DESIGN_B1.md Sec 2).  Include order: after src/types.jl + src/b1_types.jl.
#   Like types.jl this file is include()d inside `module SubTropica` after
#   `using Nemo, JSON` — no `using` statements here.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Typed refusal (SPEC_B1 §6 style: loud, typed, detail ALWAYS populated).
# Kinds: :AlgebraicLetter, :NotAPolynomial, :UnknownVariable, :TooManyVars,
#        :EmptyInput, :MixedRings.
# ---------------------------------------------------------------------------
struct LRRefineRefusal <: Exception
    kind   :: Symbol
    detail :: Dict{String,Any}
end
Base.showerror(io::IO, e::LRRefineRefusal) =
    print(io, "subtropica lr_refine: ", e.kind, " — ", e.detail, "; REFUSING")

# Hard cap on the subset DP (2^n nodes); production inputs here are n <= 7.
const _LRREFINE_MAX_VARS = 14

# ---------------------------------------------------------------------------
# Letter canonicalization.  Units in QQ[x...] are QQ^*, so "the same letter"
# means "equal up to a nonzero rational factor".  Canonical representative:
# monic w.r.t. Nemo's internal term order.  Constants and MONOMIALS are not
# letters and canonicalize to `nothing` (HyperInt irreducibles drops
# monomials, HyperInt.mpl:2618-2624; the 0/∞ coordinate hyperplanes are
# represented implicitly by the [0,f]/[inf,f] operations, [B] Lemma 67).
# ---------------------------------------------------------------------------
_is_monomial(f::QQMPolyRingElem) = length(f) == 1

function canon_letter(f::QQMPolyRingElem)
    iszero(f) && return nothing
    is_constant(f) && return nothing
    _is_monomial(f) && return nothing
    return divexact(f, leading_coefficient(f))
end

# Deterministic letter ordering for reproducible outputs / edge keys.
_letter_sortkey(p::QQMPolyRingElem) = (total_degree(p), length(p), string(p))
_sorted_letters(S) = sort!(collect(S); by = _letter_sortkey)

# ---------------------------------------------------------------------------
# irreducibles() — factor a list of polynomials over Q, canonicalize, drop
# constants and monomials; optionally intersect with `bound` (transliterates
# HyperInt.mpl:2614-2628 irreducibles(X, {}, bound); the bound path is
# SEMANTICS in the Francis rule — incompatible-pair resultants may only
# CONFIRM existing letters, never add new ones, HyperInt.mpl:2947-2955).
# ---------------------------------------------------------------------------
function _irreducibles(polys::Vector{QQMPolyRingElem};
                       bound::Union{Nothing,Set{QQMPolyRingElem}} = nothing)
    out = Set{QQMPolyRingElem}()
    for f in polys
        (iszero(f) || is_constant(f)) && continue
        if _is_monomial(f)
            continue   # monomial: every irreducible factor is a variable — drop
        end
        for (p, _e) in factor(f)
            c = canon_letter(p)
            c === nothing && continue
            if bound === nothing || (c in bound)
                push!(out, c)
            end
        end
    end
    return out
end

# ---------------------------------------------------------------------------
# Univariate-in-one-variable views: coefficient extraction, resultant
# (Sylvester determinant over the ambient mpoly ring), discriminant
# ([B] eq (65): D_x(f) = lc^-1 * Res_x(f, f')).
# ---------------------------------------------------------------------------
"""_coeff_in(f, vi, k) — coefficient of (ring var vi)^k in f, as an element of
the same ring (the variable vi does not appear in the result)."""
function _coeff_in(f::QQMPolyRingElem, vi::Int, k::Int)
    R = parent(f)
    ctx = MPolyBuildCtx(R)
    for (c, ev) in zip(coefficients(f), exponent_vectors(f))
        if ev[vi] == k
            ev2 = copy(ev); ev2[vi] = 0
            push_term!(ctx, c, ev2)
        end
    end
    return finish(ctx)
end

"""_resultant_in(f, g, vi) — Res_{x_vi}(f, g) via the Sylvester matrix over
the ambient ring ([B] Sec 6.2; degrees here are tiny, det is exact)."""
function _resultant_in(f::QQMPolyRingElem, g::QQMPolyRingElem, vi::Int)
    R = parent(f)
    n = degree(f, vi); m = degree(g, vi)
    n < 0 && return zero(R)          # f == 0
    m < 0 && return zero(R)
    fc = [_coeff_in(f, vi, k) for k in 0:n]   # fc[k+1] = coeff of vi^k
    gc = [_coeff_in(g, vi, k) for k in 0:m]
    if n == 0 && m == 0
        return one(R)
    elseif n == 0
        return fc[1]^m
    elseif m == 0
        return gc[1]^n
    end
    N = n + m
    Mat = zero_matrix(R, N, N)
    for i in 1:m, j in 0:n              # m rows of f coefficients
        Mat[i, i+j] = fc[n-j+1]
    end
    for i in 1:n, j in 0:m              # n rows of g coefficients
        Mat[m+i, i+j] = gc[m-j+1]
    end
    return det(Mat)
end

"""_discriminant_in(f, vi) — D_{x_vi}(f) = Res(f, ∂f)/lc ([B] eq (65); the
sign convention is irrelevant here since letters are tracked up to units)."""
function _discriminant_in(f::QQMPolyRingElem, vi::Int)
    R = parent(f)
    degree(f, vi) < 2 && return one(R)
    res = _resultant_in(f, derivative(f, vi), vi)
    lc  = _coeff_in(f, vi, degree(f, vi))
    ok, q = divides(res, lc)
    return ok ? q : res    # unit-level fallback; lc's factors are letters anyway
end

# ---------------------------------------------------------------------------
# Compatibility graph: node set S (canonical letters) + undirected edge set C.
# Pairs are stored sorted by _letter_sortkey so Set semantics are exact.
# ---------------------------------------------------------------------------
struct CGGraph
    S :: Set{QQMPolyRingElem}
    C :: Set{NTuple{2,QQMPolyRingElem}}
end
CGGraph() = CGGraph(Set{QQMPolyRingElem}(), Set{NTuple{2,QQMPolyRingElem}}())

function _edge(p::QQMPolyRingElem, q::QQMPolyRingElem)
    p == q && return nothing
    return _letter_sortkey(p) <= _letter_sortkey(q) ? (p, q) : (q, p)
end

"""_add_cartprod!(C, A, B) — HyperInt cartProd (HyperInt.mpl:2840-2843): all
unordered pairs {a,b}, a in A, b in B, a != b."""
function _add_cartprod!(C::Set{NTuple{2,QQMPolyRingElem}}, A, B)
    for a in A, b in B
        e = _edge(a, b)
        e === nothing || push!(C, e)
    end
    return C
end

function _add_clique!(C::Set{NTuple{2,QQMPolyRingElem}}, A)
    v = collect(A)
    for i in 1:length(v), j in i+1:length(v)
        e = _edge(v[i], v[j])
        e === nothing || push!(C, e)
    end
    return C
end

_graph_intersect(G::CGGraph, H::CGGraph) =
    CGGraph(intersect(G.S, H.S), intersect(G.C, H.C))

# ---------------------------------------------------------------------------
# cg_single_reduction — ONE projection (S,C) -> (S,C)_var.
# algorithm = :francis  (cgSingleReductionFrancis, HyperInt.mpl:2913-2986;
#                        the HyperInt DEFAULT, HyperInt.mpl:79) — realizes
#                        Brown Def 73 literally: new compatibilities from
#                        parent-sharing over ALL pairs {s1,s2},{s2,s3} with
#                        s_i in {0,∞,f_1..f_k}; incompatible-pair resultants
#                        are computed ONLY bounded by Snew (they cannot add
#                        letters, only edges).
# algorithm = :clique   (cgSingleReductionClique, HyperInt.mpl:2846-2911) —
#                        sharper: resultant-resultant edges only from
#                        TRIANGLES of C; 0/∞-shared edges only for C-pairs.
# Both prune exactly the Brown eq (71) shapes ("four distinct grandparents").
# Singleton descendants (lead/const/disc) are NEVER gated by C — a letter
# born as a leading coefficient survives every refinement (relevant for the
# recorded blocking-letter P verdict, see test_lr_refine.jl).
# var-free f: lcoeff = coeff0 = f, so f carries over as its own descendant
# under keys (:zero,f) and (:inf,f) (HyperInt.mpl:2859-2860 make no
# var-dependence exception) — edge inheritance then follows the 0/∞ rules.
# ---------------------------------------------------------------------------
function cg_single_reduction(G::CGGraph, vi::Int; algorithm::Symbol = :francis)
    algorithm in (:francis, :clique) ||
        throw(LRRefineRefusal(:NotAPolynomial,
            Dict{String,Any}("bad_algorithm" => algorithm)))
    S = _sorted_letters(G.S)
    C = G.C
    s = Dict{Any,Set{QQMPolyRingElem}}()
    nonlinear_seen = QQMPolyRingElem[]
    # --- singleton descendants (HyperInt.mpl:2855-2861 / 2921-2927) ---
    for f in S
        d = degree(f, vi)
        if d > 1
            push!(nonlinear_seen, f)
            s[(:disc, f)] = _irreducibles([_discriminant_in(f, vi)])
        end
        s[(:inf, f)]  = _irreducibles([_coeff_in(f, vi, max(d, 0))])
        s[(:zero, f)] = _irreducibles([_coeff_in(f, vi, 0)])
    end
    # --- compatible-pair resultants (HyperInt.mpl:2862-2868 / 2928-2934) ---
    respairs = NTuple{2,QQMPolyRingElem}[]
    for e in C
        f, g = e
        if degree(f, vi) > 0 && degree(g, vi) > 0
            s[(:res, f, g)] = _irreducibles([_resultant_in(f, g, vi)])
        else
            s[(:res, f, g)] = Set{QQMPolyRingElem}()
        end
        push!(respairs, e)
    end
    Snew = Set{QQMPolyRingElem}()
    for v in values(s); union!(Snew, v); end
    # --- Francis: ALL remaining pairwise resultants, factors bounded by Snew
    #     (edges only, no new letters; HyperInt.mpl:2947-2955) ---
    allpairs = respairs
    if algorithm === :francis
        allpairs = NTuple{2,QQMPolyRingElem}[]
        for i in 1:length(S), j in i+1:length(S)
            e = _edge(S[i], S[j]); e === nothing && continue
            push!(allpairs, e)
            haskey(s, (:res, e[1], e[2])) && continue
            if degree(e[1], vi) > 0 && degree(e[2], vi) > 0
                s[(:res, e[1], e[2])] = _irreducibles(
                    [_resultant_in(e[1], e[2], vi)]; bound = Snew)
            else
                s[(:res, e[1], e[2])] = Set{QQMPolyRingElem}()
            end
        end
    end
    # --- new compatibilities ---
    c = Set{NTuple{2,QQMPolyRingElem}}()
    # (1) within each factor set (HyperInt.mpl:2874 / 2958)
    for v in values(s); _add_clique!(c, v); end
    # (2) [0,f] x [inf,f] for every f (HyperInt.mpl:2876 / 2961)
    for f in S
        _add_cartprod!(c, s[(:zero, f)], s[(:inf, f)])
    end
    # (3)/(4)/(5): 0-0, inf-inf, res-x-own-coeffs; over C pairs (:clique,
    #     HyperInt.mpl:2877-2879) or ALL pairs (:francis, HyperInt.mpl:2962-2964)
    for e in allpairs
        f, g = e
        _add_cartprod!(c, s[(:zero, f)], s[(:zero, g)])
        _add_cartprod!(c, s[(:inf, f)],  s[(:inf, g)])
        if haskey(s, (:res, f, g))
            own = union(s[(:zero, f)], s[(:zero, g)], s[(:inf, f)], s[(:inf, g)])
            _add_cartprod!(c, s[(:res, f, g)], own)
        end
    end
    # (6) resultant-resultant edges from a SHARED parent — triangles
    #     {u1,u2,u3}: [u1,u2]x[u1,u3], [u2,u1]x[u2,u3], [u3,u1]x[u3,u2].
    #     :clique — only pairwise-compatible triples (HyperInt.mpl:2893-2895);
    #     :francis — all triples (HyperInt.mpl:2968-2969).  This is precisely
    #     what excludes Brown eq (71): [f1,f2]x[f3,f4] never share a parent.
    getres(a, b) = begin e = _edge(a, b); e === nothing ?
        Set{QQMPolyRingElem}() : get(s, (:res, e[1], e[2]), Set{QQMPolyRingElem}()) end
    for i in 1:length(S), j in i+1:length(S), k in j+1:length(S)
        u1, u2, u3 = S[i], S[j], S[k]
        if algorithm === :clique
            (_edge(u1,u2) in C && _edge(u1,u3) in C && _edge(u2,u3) in C) || continue
        end
        _add_cartprod!(c, getres(u1,u2), getres(u1,u3))
        _add_cartprod!(c, getres(u1,u2), getres(u2,u3))
        _add_cartprod!(c, getres(u1,u3), getres(u2,u3))
    end
    # (7) discriminant factors x var-free polynomials (HyperInt.mpl:2897-2903
    #     :clique needs {f,g} in C; HyperInt.mpl:2971-2976 :francis takes all)
    varfree = [f for f in S if degree(f, vi) == 0]
    for f in S
        degree(f, vi) > 1 || continue
        tgt = algorithm === :clique ?
            [g for g in varfree if (_edge(f, g) in C)] : varfree
        _add_cartprod!(c, s[(:disc, f)], tgt)
    end
    # (8) compatible-resultant factors x var-free polynomials
    #     (:clique — var-free g compatible with BOTH endpoints,
    #      HyperInt.mpl:2905-2908; :francis — union over C x all var-free,
    #      HyperInt.mpl:2978-2982)
    if algorithm === :francis
        resC = Set{QQMPolyRingElem}()
        for e in C; union!(resC, s[(:res, e[1], e[2])]); end
        _add_cartprod!(c, resC, varfree)
    else
        for e in C
            tgt = [g for g in varfree if (_edge(e[1], g) in C && _edge(e[2], g) in C)]
            _add_cartprod!(c, s[(:res, e[1], e[2])], tgt)
        end
    end
    # keep only edges between surviving nodes (bounded francis factor sets
    # are already subsets of Snew; cliques within s[e] too)
    return CGGraph(Snew, c), nonlinear_seen
end

# ---------------------------------------------------------------------------
# naive_single_reduction — Brown Def 68 / STFubiniLR (SubTropica.wl:16524-
# 16593): lead + const + disc + ALL pairwise resultants, NO compatibility
# tracking.  This reproduces the hyperflint find_lr_orders letter bound
# ("Fubini_Lungo") that the refinement is compared against.
# ---------------------------------------------------------------------------
function naive_single_reduction(S::Set{QQMPolyRingElem}, vi::Int)
    Sv = _sorted_letters(S)
    out = Set{QQMPolyRingElem}()
    for f in Sv
        d = degree(f, vi)
        if d > 1
            union!(out, _irreducibles([_discriminant_in(f, vi)]))
        end
        union!(out, _irreducibles([_coeff_in(f, vi, max(d, 0))]))
        union!(out, _irreducibles([_coeff_in(f, vi, 0)]))
    end
    for i in 1:length(Sv), j in i+1:length(Sv)
        if degree(Sv[i], vi) > 0 && degree(Sv[j], vi) > 0
            union!(out, _irreducibles([_resultant_in(Sv[i], Sv[j], vi)]))
        end
    end
    return out
end

# ---------------------------------------------------------------------------
# Subset DP drivers (transliterate cgReduction, HyperInt.mpl:2989-3086):
# for each subset I (ascending size), reduce every assigned parent I\{e}
# along e PROVIDED max_e deg(parent polys, e) <= max_degree (lowest degree
# first, HyperInt.mpl:3037-3065), then intersect pairwise ((S,C) ∩ rule of
# [B] Def 74).  Unreachable subsets stay unassigned (Example 4.2: L[{x}] is
# never computed).  max_degree = 1 is the HyperInt default and the
# find_lr_orders engine contract (letters deg <= 1 in each pivot).
# ---------------------------------------------------------------------------
function cg_reduction!(L::Dict{Set{Int},CGGraph}, elimvars::Vector{Int};
                       max_degree::Int = 1, algorithm::Symbol = :francis,
                       log::Union{Nothing,Vector{Dict{String,Any}}} = nothing)
    n = length(elimvars)
    for size in 1:n, I in _subsets_of_size(elimvars, size)
        haskey(L, I) && continue
        parents = Tuple{Int,Int}[]   # (var, maxdeg of parent set in var)
        for e in I
            P = setdiff(I, Set([e]))
            haskey(L, P) || continue
            H = maximum([degree(p, e) for p in L[P].S]; init = 0)
            H > max_degree && continue
            push!(parents, (e, H))
        end
        isempty(parents) && continue
        sort!(parents; by = last)     # lowest degree first (HyperInt.mpl:3046-3060)
        G = nothing
        used = Int[]
        for (e, _H) in parents
            P = setdiff(I, Set([e]))
            H, _nl = cg_single_reduction(L[P], e; algorithm = algorithm)
            G = G === nothing ? H : _graph_intersect(G, H)
            push!(used, e)
        end
        L[I] = G
        log === nothing || push!(log, Dict{String,Any}(
            "subset" => sort!(collect(I)), "parents_used" => used,
            "n_letters" => length(G.S), "n_edges" => length(G.C)))
    end
    return L
end

function naive_reduction!(L::Dict{Set{Int},Set{QQMPolyRingElem}},
                          elimvars::Vector{Int}; max_degree::Int = 1)
    n = length(elimvars)
    for size in 1:n, I in _subsets_of_size(elimvars, size)
        haskey(L, I) && continue
        G = nothing
        for e in I
            P = setdiff(I, Set([e]))
            haskey(L, P) || continue
            H = maximum([degree(p, e) for p in L[P]]; init = 0)
            H > max_degree && continue
            Hs = naive_single_reduction(L[P], e)
            G = G === nothing ? Hs : intersect(G, Hs)
        end
        G === nothing && continue
        L[I] = G
    end
    return L
end

function _subsets_of_size(v::Vector{Int}, k::Int)
    out = Set{Int}[]
    n = length(v)
    idx = collect(1:k)
    while true
        push!(out, Set(v[idx]))
        i = k
        while i >= 1 && idx[i] == n - k + i; i -= 1; end
        i == 0 && break
        idx[i] += 1
        for j in i+1:k; idx[j] = idx[j-1] + 1; end
    end
    return out
end

# ---------------------------------------------------------------------------
# Input plumbing: groups, rings, symbols.
# ---------------------------------------------------------------------------
function _check_ring(polys::Vector{QQMPolyRingElem})
    isempty(polys) && throw(LRRefineRefusal(:EmptyInput, Dict{String,Any}(
        "hint" => "no polynomials given")))
    R = parent(polys[1])
    for p in polys
        parent(p) === R || throw(LRRefineRefusal(:MixedRings, Dict{String,Any}(
            "hint" => "all polynomials must share one QQMPolyRing")))
    end
    return R
end

function _var_indices(R, vars::Vector{Symbol})
    syms = symbols(R)
    idx = Int[]
    for v in vars
        i = findfirst(==(v), syms)
        i === nothing && throw(LRRefineRefusal(:UnknownVariable,
            Dict{String,Any}("var" => v, "ring_symbols" => syms)))
        push!(idx, i)
    end
    length(unique(idx)) == length(idx) || throw(LRRefineRefusal(
        :UnknownVariable, Dict{String,Any}("hint" => "duplicate variables",
                                           "vars" => vars)))
    return idx
end

"""_initial_graph(groups) — S = irreducible factors of all inputs; C = union
of per-group complete graphs on that group's factors.  One flat group =
complete graph = [B] Def 74 init; several groups = SubTropica ct-sum init
(SubTropica.wl:15688, pairs = Subsets[polys,{2}] per locally-finite term)."""
function _initial_graph(groups::Vector{Vector{QQMPolyRingElem}})
    S = Set{QQMPolyRingElem}()
    C = Set{NTuple{2,QQMPolyRingElem}}()
    for g in groups
        gf = _irreducibles(g)
        union!(S, gf)
        _add_clique!(C, gf)
    end
    isempty(S) && throw(LRRefineRefusal(:EmptyInput, Dict{String,Any}(
        "hint" => "all inputs canonicalized away (constants/monomials only)")))
    return CGGraph(S, C)
end

_as_groups(polys::Vector{QQMPolyRingElem}) = [polys]
_as_groups(polys::Vector{Vector{QQMPolyRingElem}}) = polys

# ---------------------------------------------------------------------------
# parse_poly_qq — string -> QQMPolyRingElem, whitelisted grammar only
# (+ - * ^ / parentheses, integer literals, ring symbols, `**` accepted).
# sqrt/Root/fractional powers => :AlgebraicLetter (scope pin).
# ---------------------------------------------------------------------------
function parse_poly_qq(str::AbstractString, R::QQMPolyRing)
    s = replace(String(str), "**" => "^")
    ex = try
        Meta.parse(s)
    catch err
        throw(LRRefineRefusal(:NotAPolynomial,
            Dict{String,Any}("input" => String(str), "parse_error" => string(err))))
    end
    syms = symbols(R)
    gmap = Dict{Symbol,QQMPolyRingElem}(syms[i] => gen(R, i) for i in 1:nvars(R))
    return _expr_to_poly(ex, R, gmap, String(str))
end

function _expr_to_poly(ex, R::QQMPolyRing, gmap, orig::String)
    if ex isa Integer
        return R(QQ(ex))
    elseif ex isa Symbol
        haskey(gmap, ex) && return gmap[ex]
        throw(LRRefineRefusal(:UnknownVariable,
            Dict{String,Any}("var" => ex, "input" => orig)))
    elseif ex isa Expr && ex.head === :call
        op = ex.args[1]
        as = ex.args[2:end]
        if op === :+
            return sum(_expr_to_poly(a, R, gmap, orig) for a in as)
        elseif op === :- && length(as) == 1
            return -_expr_to_poly(as[1], R, gmap, orig)
        elseif op === :- && length(as) == 2
            return _expr_to_poly(as[1], R, gmap, orig) -
                   _expr_to_poly(as[2], R, gmap, orig)
        elseif op === :*
            return prod(_expr_to_poly(a, R, gmap, orig) for a in as)
        elseif op === :^ && length(as) == 2 && as[2] isa Integer && as[2] >= 0
            return _expr_to_poly(as[1], R, gmap, orig)^Int(as[2])
        elseif op === :/ && length(as) == 2
            den = _expr_to_poly(as[2], R, gmap, orig)
            is_constant(den) && !iszero(den) ||
                throw(LRRefineRefusal(:NotAPolynomial, Dict{String,Any}(
                    "input" => orig, "hint" => "non-constant denominator")))
            return _expr_to_poly(as[1], R, gmap, orig) *
                   R(inv(leading_coefficient(den)))
        elseif op in (:sqrt, :Sqrt, :Root, :root)
            throw(LRRefineRefusal(:AlgebraicLetter, Dict{String,Any}(
                "input" => orig,
                "hint" => "algebraic letters out of scope for this pass (FindRoots/Wm/Wp tier)")))
        end
        throw(LRRefineRefusal(:NotAPolynomial, Dict{String,Any}(
            "input" => orig, "bad_call" => string(op))))
    end
    throw(LRRefineRefusal(:NotAPolynomial, Dict{String,Any}(
        "input" => orig, "bad_node" => string(ex))))
end

# Typed refusal hooks for the algebraic-letter table objects (src/types.jl).
_letter_refusal(kind) = LRRefineRefusal(:AlgebraicLetter, Dict{String,Any}(
    "letter_kind" => kind,
    "hint" => "Wm/Wp/sqrt_disc letters are algebraic; rational-letter scope pin"))
is_blocking_real(polys, vars::Vector{Symbol}, letter::AlgLetter; kwargs...) =
    throw(_letter_refusal("AlgLetter"))
is_blocking_real(polys, vars::Vector{Symbol}, letter::HAlgLetter; kwargs...) =
    throw(_letter_refusal("HAlgLetter"))
blocking_certificate(polys, vars::Vector{Symbol}, letter::AlgLetter; kwargs...) =
    throw(_letter_refusal("AlgLetter"))
blocking_certificate(polys, vars::Vector{Symbol}, letter::HAlgLetter; kwargs...) =
    throw(_letter_refusal("HAlgLetter"))

# ---------------------------------------------------------------------------
# refine_letters(polys, vars, order) ->
#     (letters_certified, letters_spurious_removed, reduction_log)
#
#   polys : Vector{QQMPolyRingElem} (one group) or Vector{Vector{...}} (groups)
#   vars  : Vector{Symbol} — the integration variables (ring symbols); all
#           remaining ring symbols are kinematics.
#   order : Vector{Symbol} ⊆ vars — the elimination order to certify along.
#
#   letters_certified        : refined (compatibility-graph) letter bound —
#       the union over the prefix stages ∅, {o1}, {o1,o2}, ... of order of
#       the refined S_I ([P] Sec 4.1: the stage-k set bounds the alphabet of
#       the k-th partial integral).  Still an UPPER bound.
#   letters_spurious_removed : letters in the naive (find_lr_orders-class)
#       bound at those stages that the refinement REMOVED — each is
#       CERTIFIED spurious.
#   reduction_log            : Dict with per-stage detail (letters, removals,
#       admissible next variables, forced-step flags), algorithm, max_degree,
#       and citation strings.
# ---------------------------------------------------------------------------
function refine_letters(polys, vars::Vector{Symbol}, order::Vector{Symbol};
                        max_degree::Int = 1, algorithm::Symbol = :francis)
    groups = _as_groups(polys)
    R = _check_ring(reduce(vcat, groups))
    varidx = _var_indices(R, vars)
    isempty(order) && throw(LRRefineRefusal(:EmptyInput,
        Dict{String,Any}("hint" => "empty elimination order")))
    all(o -> o in vars, order) || throw(LRRefineRefusal(:UnknownVariable,
        Dict{String,Any}("hint" => "order must be a subset of vars",
                         "order" => order, "vars" => vars)))
    ordidx = _var_indices(R, order)
    length(ordidx) <= _LRREFINE_MAX_VARS || throw(LRRefineRefusal(:TooManyVars,
        Dict{String,Any}("n" => length(ordidx), "cap" => _LRREFINE_MAX_VARS)))

    G0 = _initial_graph(groups)
    Lcg = Dict{Set{Int},CGGraph}(Set{Int}() => G0)
    Lnv = Dict{Set{Int},Set{QQMPolyRingElem}}(Set{Int}() => copy(G0.S))
    dplog = Vector{Dict{String,Any}}()
    cg_reduction!(Lcg, ordidx; max_degree = max_degree, algorithm = algorithm,
                  log = dplog)
    naive_reduction!(Lnv, ordidx; max_degree = max_degree)

    syms = symbols(R)
    certified = Set{QQMPolyRingElem}()
    removed   = Set{QQMPolyRingElem}()
    stages = Vector{Dict{String,Any}}()
    walkable = true
    for k in 0:length(ordidx)
        I = Set(ordidx[1:k])
        stg = Dict{String,Any}("stage" => k,
            "eliminated" => [string(syms[i]) for i in ordidx[1:k]])
        if !haskey(Lcg, I)
            stg["computed"] = false
            stg["reason"] = "no admissible (deg <= $max_degree) projection reaches this stage"
            walkable = false
            push!(stages, stg)
            break
        end
        Scg = Lcg[I].S
        Snv = get(Lnv, I, Scg)   # naive missing can only happen if cg missing too
        union!(certified, Scg)
        rem_here = setdiff(Snv, Scg)
        union!(removed, rem_here)
        stg["computed"] = true
        stg["letters_refined"] = [string(p) for p in _sorted_letters(Scg)]
        stg["letters_naive"]   = [string(p) for p in _sorted_letters(Snv)]
        stg["letters_removed_here"] = [string(p) for p in _sorted_letters(rem_here)]
        if k < length(ordidx)
            nxt = ordidx[k+1]
            degs = Dict(string(p) => degree(p, nxt) for p in Scg)
            stg["next_var"] = string(syms[nxt])
            stg["next_var_maxdeg"] = maximum(values(degs); init = 0)
            stg["next_var_admissible"] = stg["next_var_maxdeg"] <= max_degree
            stg["blocking_letters_next"] =
                [string(p) for p in _sorted_letters(
                     [p for p in Scg if degree(p, nxt) > max_degree])]
            adm = [string(syms[e]) for e in varidx
                   if !(e in I) &&
                      maximum([degree(p, e) for p in Scg]; init = 0) <= max_degree]
            stg["admissible_next_vars"] = adm
        end
        push!(stages, stg)
    end

    reduction_log = Dict{String,Any}(
        "algorithm" => string(algorithm),
        "max_degree" => max_degree,
        "n_groups" => length(groups),
        "order" => [string(o) for o in order],
        "order_walkable" => walkable,
        "stages" => stages,
        "dp" => dplog,
        "semantics" => "refined and naive sets are BOTH upper bounds on the letter set; refined ⊆ naive; letters_spurious_removed are certified spurious",
        "citations" => [
            "Brown arXiv:0910.0114 Sec 6: Def 68 (simple), Def 73/74 (compatibility graphs + Fubini), eq (71) (excluded four-grandparent resultants)",
            "Panzer arXiv:1403.3385 Sec 4.2 (cgReduction), Example 4.2",
            "HyperInt.mpl:2846-3086 (cgSingleReductionClique/Francis + cgReduction driver; default Francis at HyperInt.mpl:79)",
            "SubTropica.wl:16524-16593 (STFubiniLR = naive baseline), SubTropica.wl:15659-15691 (multi-group initial pairs)"])
    return (_sorted_letters(certified), _sorted_letters(removed), reduction_log)
end

# ---------------------------------------------------------------------------
# blocking_certificate / is_blocking_real
#
# Adjudicates ONE candidate letter (the engine's "blocking letter") against
# the refined bound over ALL reachable elimination stages (every subset of
# `vars`, not just one order — [B] Def 74 is order-independent per subset).
#
#   status = :SURVIVES          letter in refined S_I for >= 1 reachable I
#                               => NOT certifiable spurious; blocking verdict
#                                  stands (upper-bound survival, not proof).
#   status = :SPURIOUS_REMOVED  letter in the naive S_I for some I, in NO
#                               refined S_I => CERTIFIED SPURIOUS.
#   status = :NOT_IN_BOUND      letter in neither bound at any reachable I.
#
# is_blocking_real returns (status == :SURVIVES); the certificate carries the
# witness stages either way.
# ---------------------------------------------------------------------------
function blocking_certificate(polys, vars::Vector{Symbol}, letter;
                              max_degree::Int = 1, algorithm::Symbol = :francis)
    groups = _as_groups(polys)
    R = _check_ring(reduce(vcat, groups))
    lp = letter isa AbstractString ? parse_poly_qq(letter, R) : letter
    lp isa QQMPolyRingElem || throw(LRRefineRefusal(:NotAPolynomial,
        Dict{String,Any}("letter_type" => string(typeof(letter)))))
    parent(lp) === R || throw(LRRefineRefusal(:MixedRings,
        Dict{String,Any}("hint" => "letter must live in the polys' ring")))
    lc = canon_letter(lp)
    lc === nothing && throw(LRRefineRefusal(:NotAPolynomial,
        Dict{String,Any}("hint" => "letter is a constant or monomial",
                         "letter" => string(lp))))
    varidx = _var_indices(R, vars)
    length(varidx) <= _LRREFINE_MAX_VARS || throw(LRRefineRefusal(:TooManyVars,
        Dict{String,Any}("n" => length(varidx), "cap" => _LRREFINE_MAX_VARS)))

    G0 = _initial_graph(groups)
    Lcg = Dict{Set{Int},CGGraph}(Set{Int}() => G0)
    Lnv = Dict{Set{Int},Set{QQMPolyRingElem}}(Set{Int}() => copy(G0.S))
    cg_reduction!(Lcg, varidx; max_degree = max_degree, algorithm = algorithm)
    naive_reduction!(Lnv, varidx; max_degree = max_degree)

    syms = symbols(R)
    stagename(I) = sort!([string(syms[i]) for i in I])
    _stagesort!(v) = sort!(v; by = x -> (length(x), join(x, ",")))
    in_refined = [stagename(I) for (I, G) in Lcg if lc in G.S]
    in_naive   = [stagename(I) for (I, S) in Lnv if lc in S]
    status = !isempty(in_refined) ? :SURVIVES :
             (!isempty(in_naive) ? :SPURIOUS_REMOVED : :NOT_IN_BOUND)
    return Dict{String,Any}(
        "letter" => string(lc),
        "status" => string(status),
        "real" => status === :SURVIVES,
        "stages_refined" => _stagesort!(in_refined),
        "stages_naive" => _stagesort!(in_naive),
        "stages_reachable_refined" => _stagesort!([stagename(I) for I in keys(Lcg)]),
        "algorithm" => string(algorithm),
        "max_degree" => max_degree,
        "meaning" => status === :SURVIVES ?
            "letter survives the compatibility-graph refinement: NOT certifiable spurious (upper-bound survival, not a positive proof of singularity)" :
            (status === :SPURIOUS_REMOVED ?
             "letter is in the naive Fubini bound but in NO refined stage: CERTIFIED SPURIOUS (Brown Def 73/74 upper-bound theorem)" :
             "letter is not produced by the reduction at any reachable stage"))
end

function is_blocking_real(polys, vars::Vector{Symbol}, letter;
                          max_degree::Int = 1, algorithm::Symbol = :francis)
    cert = blocking_certificate(polys, vars, letter;
                                max_degree = max_degree, algorithm = algorithm)
    return cert["real"]::Bool
end
