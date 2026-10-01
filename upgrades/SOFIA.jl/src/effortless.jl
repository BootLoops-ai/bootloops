# Odd symbol letters: port of Effortless.m (Antonela Matijašić, Julian
# Miczajka; MIT), bundled with SOFIA as EffortlessMarch2025/.
#
# An odd letter for a square root sqrt(Q) is W = (P - sqrt(Q))/(P + sqrt(Q))
# where P^2 - Q factors as c * (product of even letters), i.e. P^2 = Q + Prod
# is a perfect square for a candidate product Prod of even alphabet letters.
#
# Port decisions (documented divergences, both conservative):
# - Upstream tests perfect-squareness at 4 random integer points
#   (`AllowedProductQ`); we use exact `Nemo.is_square`.
# - Upstream pre-filters even letters via numeric zero-locus tests
#   (`AllowedEvenLetters`, degree < 5 only); we skip the pre-filter — it only
#   prunes candidates that the exact square test rejects anyway, so the
#   resulting letter set is identical or a superset.
# - Independence (`FindIndependentLetters`) is decided by exact linear algebra
#   on dlog-forms evaluated at deterministic rational points, replacing
#   NSolve at 50 random real configurations.

"""
    OddLetter

W = (P - sqrt(Q))/(P + sqrt(Q)): stores `P`, the squared root argument `Q`,
and the even-product witness `prod` with coefficient `c` (P^2 = Q + c*prod).
"""
struct OddLetter
    P::Any
    Q::Any
    c::Any
    prod::Any
end

"""
    find_odd_letters(Q, alphabet; coeffs=(4,-4)) -> Vector{OddLetter}

All odd letters for the squared root `Q`: candidate products of alphabet
letters whose total degree matches deg(Q), scaled by `c in coeffs`, such that
Q + c*prod is a perfect square (port of `FindOddLetters`/`ConstructProducts`/
`AllowedProductQ`).
"""
function find_odd_letters(Q, alphabet::AbstractVector; coeffs=(4, -4),
                          maxfactors::Int=4)
    R = parent(Q)
    q = Nemo.total_degree(Q)
    q > 0 || return OddLetter[]
    bydeg = Dict{Int,Vector{Any}}()
    for L in alphabet
        d = Nemo.total_degree(L)
        d > 0 && push!(get!(bydeg, d, Any[]), L)
    end
    degs = sort!(collect(keys(bydeg)))
    out = OddLetter[]
    seen = Set{Any}()
    for part in integer_partitions_capped(q, degs, maxfactors)
        for combo in product_combinations(part, bydeg)
            prod_ = one(R)
            for L in combo
                prod_ *= L
            end
            for c in coeffs
                cand = Q + c * prod_
                iszero(cand) && continue
                if Nemo.is_square(cand)
                    P = Nemo.sqrt(cand)
                    k = normal_form_key(P)
                    k in seen && continue
                    push!(seen, k)
                    push!(out, OddLetter(P, Q, c, prod_))
                end
            end
        end
    end
    return out
end

"""Partitions of `q` into parts drawn from `degs` (multiset, sorted,
at most `maxn` parts)."""
function integer_partitions_capped(q::Int, degs::Vector{Int}, maxn::Int)
    out = Vector{Vector{Int}}()
    function rec(rem, minidx, acc)
        if rem == 0
            push!(out, copy(acc))
            return
        end
        length(acc) >= maxn && return
        for i in minidx:length(degs)
            d = degs[i]
            d > rem && continue
            push!(acc, d)
            rec(rem - d, i, acc)
            pop!(acc)
        end
    end
    rec(q, 1, Int[])
    return out
end

"""All multisets of letters matching a degree partition (combinations with
repetition within each degree class)."""
function product_combinations(part::Vector{Int}, bydeg::Dict{Int,Vector{Any}})
    groups = Dict{Int,Int}()
    for d in part
        groups[d] = get(groups, d, 0) + 1
    end
    pools = Vector{Vector{Vector{Any}}}()
    for (d, k) in sort!(collect(groups))
        letters = bydeg[d]
        combos = Vector{Vector{Any}}()
        idx = ones(Int, k)
        n = length(letters)
        while true
            push!(combos, letters[idx])
            j = k
            while j >= 1 && idx[j] == n
                j -= 1
            end
            j == 0 && break
            idx[j] += 1
            for t in (j+1):k
                idx[t] = idx[j]
            end
        end
        push!(pools, combos)
    end
    isempty(pools) && return Vector{Vector{Any}}()
    out = Vector{Vector{Any}}()
    function rec(i, acc)
        if i > length(pools)
            push!(out, reduce(vcat, acc; init=Any[]))
            return
        end
        for c in pools[i]
            push!(acc, c)
            rec(i + 1, acc)
            pop!(acc)
        end
    end
    rec(1, Vector{Any}[])
    return out
end

# ---------- independence of letters via dlog forms ----------

"""
    dlog_odd(letter::OddLetter, vi) -> rational function

sqrt(Q) * d/dx_vi log((P - sqrt(Q))/(P + sqrt(Q))), which is rational:
  (Q dP/dx - P/2 dQ/dx) * 2 / (P^2 - Q)
(specialization of Effortless.m:302 `TakeDerivativeOfOddModSqrt` with p=1).
Returns (numerator, denominator) polynomials.
"""
function dlog_odd(l::OddLetter, vi::Int)
    dP = Nemo.derivative(l.P, vi)
    dQ = Nemo.derivative(l.Q, vi)
    num = 2 * (l.Q * dP) - l.P * dQ
    den = l.P^2 - l.Q
    return (num, den)
end

"""dlog of an even letter: (dL/dx, L)."""
dlog_even(L, vi::Int) = (Nemo.derivative(L, vi), L)

"""
    independent_letters(odd::Vector{OddLetter}, evens; npoints=6) -> Vector

Greedy independence filtering of odd letters sharing one root (port of
`FindIndependentLetters`): a letter is kept iff its dlog vector (sampled
exactly at deterministic rational points, all variables) is linearly
independent of the dlogs of the even alphabet and previously kept odd
letters. Exact rank computation over QQ.
"""
function independent_letters(odd::Vector{OddLetter}, evens::AbstractVector;
                             npoints::Int=8)
    isempty(odd) && return OddLetter[]
    R = parent(odd[1].P)
    nv = Nemo.nvars(R)
    # deterministic sample points (avoid zeros of denominators by retry)
    pts = Vector{Vector{Nemo.QQFieldElem}}()
    state = UInt64(0x9e3779b97f4a7c15)
    nextint() = begin
        state = state * 6364136223846793005 + 1442695040888963407
        Int(mod(state >> 33, 4000)) - 2000
    end
    alldens = vcat([l.P^2 - l.Q for l in odd], collect(evens))
    while length(pts) < npoints
        pt = [Nemo.QQ(nextint(), 1 + abs(nextint()) % 97) for _ in 1:nv]
        any(p -> iszero(Nemo.evaluate(p, pt)), alldens) && continue
        push!(pts, pt)
    end
    evalrf(nd, pt) = Nemo.evaluate(nd[1], pt) // Nemo.evaluate(nd[2], pt)
    # base span: even-letter dlogs
    function dlogvec(forms)
        v = Nemo.QQFieldElem[]
        for pt in pts, vi in 1:nv
            push!(v, evalrf(forms(vi), pt))
        end
        return v
    end
    rows = Vector{Vector{Nemo.QQFieldElem}}()
    for L in evens
        is_constant(L) && continue
        push!(rows, dlogvec(vi -> dlog_even(L, vi)))
    end
    kept = OddLetter[]
    for l in odd
        cand = dlogvec(vi -> dlog_odd(l, vi))
        if !in_rowspan(rows, cand)
            push!(rows, cand)
            push!(kept, l)
        end
    end
    return kept
end

"""Exact membership of `v` in the row span of `rows` over QQ."""
function in_rowspan(rows::Vector{Vector{Nemo.QQFieldElem}}, v::Vector{Nemo.QQFieldElem})
    isempty(rows) && return all(iszero, v)
    M = Nemo.matrix(Nemo.QQ, length(rows), length(rows[1]),
                    reduce(vcat, rows))
    Mv = Nemo.matrix(Nemo.QQ, length(rows) + 1, length(v),
                     vcat(reduce(vcat, rows), v))
    return Nemo.rank(M) == Nemo.rank(Mv)
end

"""
    effortless_odd_letters(singularities; doubleroots=true) -> Vector{OddLetter}

The full Effortless pipeline as SOFIA invokes it (`RunEffortlessIndependent`
with DoubleSquareRoots->True): roots are the non-square singularity
polynomials (and pairwise products), the even alphabet is the singularity
list itself; returns independent odd letters.
"""
function effortless_odd_letters(sing::AbstractVector; doubleroots::Bool=true,
                                maxfactors::Int=4)
    isempty(sing) && return OddLetter[]
    # roots: singularities that are not perfect squares (sqrt(sing) stays a
    # genuine root), plus pairwise products when doubleroots
    roots = Any[]
    for p in sing
        is_constant(p) && continue
        Nemo.is_square(p) || push!(roots, p)
    end
    if doubleroots
        n0 = length(roots)
        for i in 1:n0, j in (i+1):n0
            q = roots[i] * roots[j]
            Nemo.is_square(q) || push!(roots, q)
        end
    end
    out = OddLetter[]
    for Q in dedup_proportional(roots)
        odd = find_odd_letters(Q, sing; maxfactors)
        append!(out, independent_letters(odd, sing))
    end
    return out
end
