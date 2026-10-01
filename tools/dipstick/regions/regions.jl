#!/usr/bin/env julia
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision.
# regions.jl — dipstick member `regions`: region-completeness certifier for
# expansion-by-regions (formerly the standalone ebr-complete package).
# Identity: lower facets of Newt(G)⊂ℝ^{n+1} (last coord = small-param weight) tile the
# projection ⇒  Σ_R n!Vol(π(F_R)) = n!Vol(Newt(G))  and (Kouchnirenko)  Σ_R χ_reg(G_R)=χ_reg(G).
# Missing regions ⇒ LHS<RHS. Catches the single-region trap (a region-incomplete
# series converges, but to the wrong value).
#   dipstick regions <in.json> [more.json ...]
#   julia --project=<your Oscar project> regions.jl <in.json> [more.json ...]
# Several input files run in ONE process (one Oscar load); each report starts with
# its own "== dipstick regions: <file> ==" header.
# DIPSTICK_SKIP_CHI=1 (alias EBR_SKIP_CHI=1, accepted for one release): skip ALL
#   chi_reg Groebner work (vol-tiling certificate only; region lines print chi=NA).
#   For large families (ten or more edges) where the global chi Groebner basis
#   dominates the runtime AFTER volumes land. Consumers MUST treat skipped chi as
#   unverified (flag -> value-gate before use).
# Input JSON: {"G":str, "kin":[..], "x":[..], "kin_weight":{kin:int,..}, "regions":[[v..],..]?}
# Leviathan (chi_regulated engine): resolved from a sibling tools/leviathan, then
# upgrades/leviathan in this repo. Override with DIPSTICK_LEVIATHAN (the old name
# EBR_LEVIATHAN is still accepted as a fallback alias for one release).
let _lev = get(ENV, "DIPSTICK_LEVIATHAN", get(ENV, "EBR_LEVIATHAN", ""))
    if isempty(_lev)
        for _c in (joinpath(@__DIR__, "leviathan", "src", "Leviathan.jl"),
                   joinpath(@__DIR__, "..", "..", "leviathan", "src", "Leviathan.jl"),
                   joinpath(@__DIR__, "..", "..", "..", "upgrades", "leviathan", "src", "Leviathan.jl"))
            if isfile(_c)
                _lev = _c
                break
            end
        end
    end
    isempty(_lev) && error("Leviathan.jl not found — set DIPSTICK_LEVIATHAN to its path")
    include(_lev)
end
using .Leviathan, Random, JSON
import Oscar

veq(a, b) = length(a) == length(b) && all(Rational(a[i]) == Rational(b[i]) for i in 1:length(a))

function term_string(fam, c, e)  # rebuild a single G-monomial as a parseable string
    nk = length(fam.kin)
    parts = String[string(c)]
    for j in 1:nk            e[j]      > 0 && push!(parts, fam.kinnames[j] * "^" * string(e[j])) end
    for j in 1:length(fam.xs) e[nk+j]  > 0 && push!(parts, fam.xnames[j]   * "^" * string(e[nk+j])) end
    join(parts, "*")
end

const SKIP_CHI = get(ENV, "DIPSTICK_SKIP_CHI", get(ENV, "EBR_SKIP_CHI", "0")) == "1"

function chi_of_Gstring(Gstr, kin, x; rng)
    SKIP_CHI && return nothing
    fam = family(Gstr, kin, x)
    c = chi_regulated(fam; rng=rng)
    c === nothing ? 0 : c
end

chistr(c) = c === nothing ? "NA" : string(c)

function lower_facets(pts::Matrix{Int})
    P = Oscar.convex_hull(pts); n = size(pts, 2) - 1
    A, b = Oscar.halfspace_matrix_pair(Oscar.facets(P))
    out = Tuple{Vector{Rational{Int}},Vector{Int},Rational{Int}}[]
    for i in 1:Oscar.nrows(A)
        a = [Rational(A[i, j]) for j in 1:n+1]; bi = Rational(b[i])
        a[end] >= 0 && continue
        v = [a[j] // a[end] for j in 1:n]   # facet normal ∝ -(v,1): minimizer of w+v·e
        on = [k for k in 1:size(pts, 1) if sum(a[j] * pts[k, j] for j in 1:n+1) == bi]
        proj = Oscar.convex_hull(pts[on, 1:n])
        nv = Oscar.dim(proj) == n ? Rational(factorial(n) * Oscar.volume(proj)) : 0 // 1
        push!(out, (v, on, nv))
    end
    out
end

function main(path::AbstractString)
    inp = JSON.parsefile(path); rng = MersenneTwister(get(inp, "seed", 20260630))
    kin, x = String.(inp["kin"]), String.(inp["x"]); n = length(x); nk = length(kin)
    kw = Dict(k => Int(v) for (k, v) in get(inp, "kin_weight", Dict()))
    fam = family(inp["G"], kin, x)
    # monomial table: (coeff, full-exponent, x-exponent, t-weight)
    coefs, exps, xexp, tw = [], Vector{Int}[], Vector{Int}[], Int[]
    for (c, e) in zip(Oscar.coefficients(fam.G), Oscar.exponents(fam.G))
        push!(coefs, c); push!(exps, collect(Int, e))
        push!(xexp, [e[nk+j] for j in 1:n]); push!(tw, sum(get(kw, kin[j], 0) * e[j] for j in 1:nk))
    end
    M = length(coefs)
    pts = reduce(vcat, reshape(vcat(xexp[k], tw[k]), 1, n + 1) for k in 1:M)
    facets = lower_facets(pts)
    fullvol = Rational(factorial(n) * Oscar.volume(Oscar.convex_hull(pts[:, 1:n])))
    chi_full = chi_of_Gstring(inp["G"], kin, x; rng=rng)
    # report enumerated facets
    println("== dipstick regions: ", path, " ==")
    println("n=", n, "  monomials=", M, "  n!Vol(Newt(G))=", fullvol, "  chi_reg(G)=", chistr(chi_full))
    println("\n-- lower facets (independent enumeration) --")
    GRstr = String[]; chiR = Union{Int,Nothing}[]
    for (i, (v, on, nv)) in enumerate(facets)
        gs = join((term_string(fam, coefs[k], exps[k]) for k in on), " + ")
        push!(GRstr, gs)
        cr = chi_of_Gstring(gs, kin, x; rng=rng); push!(chiR, cr)
        scaleless = nv == 0 ? "  [scaleless]" : ""
        println("R", i, ": v=", v, "  n!Vol=", nv, "  chi=", chistr(cr), scaleless,
                "\n    G_R = ", length(gs) > 100 ? gs[1:100] * "..." : gs)
    end
    sv = sum(f[3] for f in facets; init=0//1)
    println("\nSigma n!Vol = ", sv, " vs ", fullvol, sv == fullvol ? "  [TILES OK]" : "  [TILE GAP!]")
    if SKIP_CHI
        println("Sigma chi   = SKIPPED (DIPSTICK_SKIP_CHI=1; vol certificate only)")
    else
        sc = sum(chiR; init=0)
        println("Sigma chi   = ", sc, " vs ", chi_full,
                sc == chi_full ? "  [CHI ADDITIVITY OK]" : "  [CHI MISMATCH (degenerate G_R?)]")
    end
    # user-supplied regions
    if haskey(inp, "regions")
        ureg = [Rational.(r) for r in inp["regions"]]
        println("\n-- user regions vs enumeration --")
        uvol, uchi, missing_f, extra = 0 // 1, 0, Int[], Int[]
        matched = falses(length(ureg))
        for (i, (v, on, nv)) in enumerate(facets)
            j = findfirst(u -> veq(u, v), ureg)
            if j === nothing
                push!(missing_f, i)
            else
                matched[j] = true; uvol += nv; uchi += something(chiR[i], 0)
            end
        end
        extra = [j for j in 1:length(ureg) if !matched[j]]
        for j in extra  # user region not a lower facet — compute its G_R anyway
            v = ureg[j]; d = denominator(sum(v) + 1)  # not used; just report
            on = let mn = minimum(sum(v[l] * xexp[k][l] for l in 1:n) + tw[k] for k in 1:M)
                [k for k in 1:M if sum(v[l] * xexp[k][l] for l in 1:n) + tw[k] == mn]
            end
            gs = join((term_string(fam, coefs[k], exps[k]) for k in on), " + ")
            cr = chi_of_Gstring(gs, kin, x; rng=rng); uvol_j = 0 // 1
            proj = Oscar.convex_hull(pts[on, 1:n])
            Oscar.dim(proj) == n && (uvol_j = Rational(factorial(n) * Oscar.volume(proj)))
            uvol += uvol_j; uchi += something(cr, 0)
            println("user R", j, " v=", v, " NOT a lower facet; n!Vol=", uvol_j, " chi=", chistr(cr))
        end
        println("user Sigma n!Vol = ", uvol, " / ", fullvol,
                "  user Sigma chi = ", SKIP_CHI ? "NA" : string(uchi), " / ", chistr(chi_full))
        ok = uvol == fullvol && isempty(missing_f)
        println("\nVERDICT: ", ok ? "PASS (region list complete)" :
                "FAIL (region list INCOMPLETE)")
        !isempty(missing_f) && println("  MISSING facets: ",
            join(("R$i v=$(facets[i][1]) n!Vol=$(facets[i][3])" for i in missing_f), "; "))
        !isempty(extra) && println("  user regions not matching any lower facet: ", extra)
        !ok && println("  vol deficit = ", fullvol - uvol, "  chi deficit = ",
                       SKIP_CHI ? "NA" : string(chi_full - uchi))
    end
end

if isempty(ARGS)
    println(stderr, "usage: dipstick regions <in.json> [more.json ...]   (see tools/dipstick/GUIDE.md, `regions` member)")
    exit(2)
end
for (i, path) in enumerate(ARGS)
    i > 1 && println()
    main(path)
end
