# MpolyFeed.jl — bivariate pre-cancel feed for connections stored as raw P/Q chunk lists.
#
# Problem this solves (measured on a 23x23 core connection): IBP/transport-derived on-curve
# kernels arrive as chunks eps^sh * P(eps,x)/Q(eps,x) with UNCANCELLED gcd(P,Q) and
# degrees up to (x:225, eps:43).  Feeding those through string parsing + Qex Horner makes
# Nemo's rational_function_field do bivariate gcds on the raw sizes: measured n=4 block
# build >40 min CPU, unusable.  Cure: build P,Q as FLINT multivariate polynomials
# (QQMPolyRingElem over QQ[eps,x]), gcd-cancel THERE (FLINT mpoly gcd is fast), and only
# then convert the small reduced numerator/denominator into the engine field Qex.
#
# API:
#   feed_matrix_from_chunks(ctx, n, chunks; shift_key="sh", pkey="P", qkey="Q") -> Qex matrix
#     chunks: iterable of Dict-like with keys i, j (0-based Ints), sh (eps_shift), and
#     P, Q = vectors over eps-power of ascending-x coefficient STRINGS ("a" or "a/b").
#   qex_from_mpoly_pair(ctx, N, D, R) -> Qex element of the reduced pair.
#
# Standalone feed layer: the core Counterweight modules do not depend on it.

module MpolyFeed

using Nemo
using ..RatFunc

export feed_matrix_from_chunks, qex_from_mpoly_pair, feed_cache_save, feed_cache_load

_qq(s::AbstractString) = begin
    if occursin("/", s)
        p = split(s, "/")
        QQFieldElem(parse(BigInt, p[1]), parse(BigInt, p[2]))
    else
        QQFieldElem(parse(BigInt, s))
    end
end

"Build one mpoly in R=QQ[eps,x] from rows: vector over eps-power of ascending-x coeff strings."
function _mpoly_from_rows(R, rows)
    b = MPolyBuildCtx(R)
    for (m, row) in enumerate(rows)          # m-1 = eps power
        for (k, cs) in enumerate(row)        # k-1 = x power
            scs = String(cs)
            scs == "0" && continue
            push_term!(b, _qq(scs), [m - 1, k - 1])
        end
    end
    return finish(b)
end

"Convert a reduced mpoly (in QQ[eps,x]) into ctx.Qex via precomputed generator powers."
function _qex_of_mpoly(ctx::Ctx, F, epspow::Vector, xpow::Vector)
    v = zero(ctx.Qex)
    for t in 1:length(F)
        ev = exponent_vector(F, t)
        c = coeff(F, t)
        # grow power caches on demand
        while length(epspow) <= ev[1]
            push!(epspow, epspow[end] * ctx.Qex(ctx.eps))
        end
        while length(xpow) <= ev[2]
            push!(xpow, xpow[end] * ctx.x)
        end
        v += ctx.Qex(c) * epspow[ev[1] + 1] * xpow[ev[2] + 1]
    end
    return v
end

"""
    _qex_pair_fast(ctx, N, D) -> Qex element N/D  (no fraction-field gcd)

PRECONDITION: (N,D) coprime in QQ[eps,x] — guaranteed by the FLINT mpoly gcd in
`feed_matrix_from_chunks`; by Gauss's lemma they are then coprime in Q(eps)[x].
AA's `//` would redo that gcd generically over Q(eps)[x]: measured 148.5s of a
150.6s n=4 feed.  Instead: build the
numerator/denominator as (Q(eps))[x] polys per-x-degree (coefficients = QQ[eps]
polys, denominator-free Qe elements) and construct the fraction directly.
Downstream behavior is IDENTICAL: `numerator`/`denominator` unit-normalize on
access (divide by lc(den), AA generic/Fraction.jl:27), which is exactly the
canonical form `//` would have stored for a coprime pair.
"""
function _qex_pair_fast(ctx::Ctx, N, D)
    Pe = parent(numerator(ctx.eps))       # QQ[eps]
    Px = parent(numerator(ctx.x))         # (Q(eps))[x]
    Qe = ctx.Qe
    function poly_of(F)
        store = Dict{Int,Dict{Int,QQFieldElem}}()
        maxx = 0
        for t in 1:length(F)
            evv = exponent_vector(F, t)
            d = get!(store, evv[2], Dict{Int,QQFieldElem}())
            d[evv[1]] = coeff(F, t)
            maxx = max(maxx, evv[2])
        end
        coeffs = Vector{elem_type(Qe)}(undef, maxx + 1)
        for k in 0:maxx
            if haskey(store, k)
                dd = store[k]
                me = maximum(keys(dd))
                coeffs[k+1] = Qe(Pe([get(dd, i, QQFieldElem(0)) for i in 0:me]))
            else
                coeffs[k+1] = zero(Qe)
            end
        end
        return Px(coeffs)
    end
    np = poly_of(N)
    dp = poly_of(D)
    iszero(dp) && error("zero denominator in reduced pair")
    fr = Nemo.AbstractAlgebra.Generic.FracFieldElem{typeof(np)}(np, dp)
    fr.parent = parent(getfield(ctx.x, :d))   # underlying FracField of Qex
    return ctx.Qex(fr)   # public wrap, RationalFunctionField.jl:577 — NO gcd
end

"Reduced pair (N,D) mpolys -> Qex element N/D.  (Now via the no-gcd fast path.)"
function qex_from_mpoly_pair(ctx::Ctx, N, D, R)
    return _qex_pair_fast(ctx, N, D)
end

const _CACHE_MAGIC = "MPOLYFEEDCACHE v1"

"""
    feed_cache_save(path, n, acc)

Serialize the REDUCED per-entry pairs (N,D) — the output of the expensive
chunk-sum + FLINT-gcd phase — to a sidecar text file.  Exact: every coefficient
is written as decimal numerator/denominator (BigInt round-trip, no string-float
path anywhere).  Format v1:
    MPOLYFEEDCACHE v1 / n / nentries / per entry: "i j #termsN #termsD" then one
    term per line: "num/den e_eps e_x" (N terms first, then D terms).
"""
function feed_cache_save(path::AbstractString, n::Int, acc)
    open(path, "w") do io
        println(io, _CACHE_MAGIC)
        println(io, n)
        println(io, length(acc))
        for ((i, j), (N, D)) in sort!(collect(acc); by=first)
            println(io, i, " ", j, " ", length(N), " ", length(D))
            for F in (N, D)
                for t in 1:length(F)
                    ev = exponent_vector(F, t)
                    c = coeff(F, t)
                    println(io, numerator(c), "/", denominator(c), " ",
                            ev[1], " ", ev[2])
                end
            end
        end
    end
    return path
end

"""
    feed_cache_load(path, R) -> (n, acc)

Reload reduced pairs into ring R = QQ[eps,x].  Exact BigInt/QQFieldElem
round-trip.  Callers MUST rerun their exact selftest on the rebuilt matrix
(the v4 driver does; the 2-point exact selftest compares against the raw
chunks, so it certifies the cache end-to-end).
"""
function feed_cache_load(path::AbstractString, R)
    acc = Dict{Tuple{Int,Int},Tuple{Any,Any}}()
    n = 0
    open(path, "r") do io
        magic = readline(io)
        magic == _CACHE_MAGIC || error("bad cache magic in $path: $magic")
        n = parse(Int, readline(io))
        k = parse(Int, readline(io))
        for _ in 1:k
            hdr = split(readline(io))
            i = parse(Int, hdr[1]); j = parse(Int, hdr[2])
            tn = parse(Int, hdr[3]); td = parse(Int, hdr[4])
            polys = Any[]
            for cnt in (tn, td)
                b = MPolyBuildCtx(R)
                for _ in 1:cnt
                    parts = split(readline(io))
                    nd = split(parts[1], "/")
                    c = QQFieldElem(parse(BigInt, nd[1]), parse(BigInt, nd[2]))
                    push_term!(b, c, [parse(Int, parts[2]), parse(Int, parts[3])])
                end
                push!(polys, finish(b))
            end
            acc[(i, j)] = (polys[1], polys[2])
        end
    end
    return n, acc
end

"""
    feed_matrix_from_chunks(ctx, n, chunks; cache=nothing) -> Qex n×n matrix

Entry (i,j) = sum over its chunks of eps^sh * P/Q, with all cancellation done in
QQ[eps,x] (FLINT mpoly gcd) BEFORE entering the Qex fraction field.  Negative eps_shift
goes to the denominator (mpolys have no negative exponents).

`cache`: optional sidecar path (serialize-once).  If the file exists the
expensive chunk-sum/gcd phase is SKIPPED and the reduced pairs are reloaded exactly;
otherwise the pairs are computed as before and saved there.  Backward compatible:
default `nothing` = original behavior, no code path touched.
"""
function feed_matrix_from_chunks(ctx::Ctx, n::Int, chunks;
                                 cache::Union{Nothing,AbstractString}=nothing)
    R, (ev, xv) = polynomial_ring(QQ, ["eps", "x"])
    local acc
    if cache !== nothing && isfile(cache)
        ncache, acc = feed_cache_load(cache, R)
        ncache == n || error("feed cache n=$ncache mismatches requested n=$n")
    else
    # accumulate per-entry reduced pairs
    acc = Dict{Tuple{Int,Int},Tuple{Any,Any}}()
    for ch in chunks
        i = Int(ch["i"]); j = Int(ch["j"]); sh = Int(get(ch, "sh", 0))
        P = _mpoly_from_rows(R, ch["P"])
        Q = _mpoly_from_rows(R, ch["Q"])
        iszero(Q) && error("chunk ($i,$j): zero denominator")
        if sh > 0
            P *= ev^sh
        elseif sh < 0
            Q *= ev^(-sh)
        end
        g = gcd(P, Q)
        if !isone(g)
            P = divexact(P, g); Q = divexact(Q, g)
        end
        if haskey(acc, (i, j))
            N0, D0 = acc[(i, j)]
            g2 = gcd(D0, Q)
            D0r = divexact(D0, g2); Qr = divexact(Q, g2)
            N = N0 * Qr + P * D0r
            D = D0 * Qr
            g3 = gcd(N, D)
            if !isone(g3)
                N = divexact(N, g3); D = divexact(D, g3)
            end
            acc[(i, j)] = (N, D)
        else
            acc[(i, j)] = (P, Q)
        end
    end
    if cache !== nothing
        feed_cache_save(cache, n, acc)
    end
    end
    A = zero_matrix(ctx.Qex, n, n)
    for ((i, j), (N, D)) in acc
        A[i + 1, j + 1] = _qex_pair_fast(ctx, N, D)
    end
    return A
end

end # module
