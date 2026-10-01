# staircase.jl — standard-monomial count from DRL leading exponents.
# Returns: 0 if unit ideal; `nothing` if not zero-dimensional (some variable has no
# pure-power leading monomial => Infinity/Indeterminate, the Type-1.1 detector);
# else the number of lattice points under the staircase.
# Counting is an exact interval-sliced recursion: slice along the last variable at the
# distinct generator exponents; between consecutive thresholds the (n-1)-variable slice
# ideal is constant, so each interval contributes (length) * (slice count). Cost scales
# with the staircase's threshold structure, not the box volume, so arbitrarily large
# staircase boxes count exactly. Counts wider than Int are returned as Int128.

"Drop leading exponents that are divisible by another (keep minimal generators)."
function _minimalize(les::Vector{Vector{Int}})
    keep = Vector{Vector{Int}}()
    for e in sort(les, by=sum)
        any(k -> all(e .>= k), keep) && continue
        push!(keep, e)
    end
    keep
end

function staircase_count(les::Vector{Vector{Int}}, nvars::Int)
    isempty(les) && return nothing
    any(e -> all(==(0), e), les) && return 0          # unit ideal
    les = _minimalize(les)
    # zero-dimensional <=> every variable has a pure-power leading monomial
    for j in 1:nvars
        any(e -> e[j] != 0 && all(k -> k == j || e[k] == 0, 1:nvars), les) ||
            return nothing                            # positive-dimensional
    end
    c = try
        _staircase_sliced(les, nvars)
    catch err
        err isa OverflowError && error("staircase count exceeds Int128")
        rethrow()
    end
    c <= Int128(typemax(Int)) ? Int(c) : c
end

"""
Exact standard-monomial count of a minimalized zero-dimensional staircase in variables
1..n. Recursion on the last variable: for exponent k the standard monomials with
x_n-degree k are those of the slice ideal {e[1:n-1] : e[n] <= k}, which only changes at
the distinct e[n] values, so sum (interval length) * (slice count) over the intervals
between consecutive thresholds. Slices stay zero-dimensional: the pure powers of the
other variables carry e[n] == 0 and enter every slice. Accumulates in checked Int128.
"""
function _staircase_sliced(les::Vector{Vector{Int}}, n::Int)::Int128
    any(e -> all(==(0), e), les) && return 0          # unit slice ideal
    n == 1 && return Int128(minimum(e[1] for e in les))
    # pure-power bound in variable n (standard monomials have e[n] < b)
    pp = [e[n] for e in les if all(k -> k == n || e[k] == 0, 1:n)]
    isempty(pp) && error("staircase slice lost zero-dimensionality (internal invariant)")
    b = minimum(pp)
    ts = sort!(unique!(Int[e[n] for e in les if e[n] < b]))
    # the other variables' pure powers sit at e[n] == 0, so intervals start at 0
    @assert !isempty(ts) && ts[1] == 0
    total = Int128(0)
    for (i, lo) in enumerate(ts)
        hi = i < length(ts) ? ts[i + 1] : b           # slice constant on [lo, hi-1]
        sub = _minimalize([e[1:n-1] for e in les if e[n] <= lo])
        total = Base.checked_add(total,
                    Base.checked_mul(Int128(hi - lo), _staircase_sliced(sub, n - 1)))
    end
    total
end

"a < b in degrevlex? (higher total degree wins; ties: last nonzero of a-b positive => a < b)"
function _drl_less(a::Vector{Int}, b::Vector{Int})
    sa, sb = sum(a), sum(b)
    sa != sb && return sa < sb
    for j in length(a):-1:1
        d = a[j] - b[j]
        d != 0 && return d > 0
    end
    false
end

"Leading exponents under DRL, computed directly (Oscar 1.7 kwarg API is inconsistent)."
function lead_exps(gb, R)
    out = Vector{Vector{Int}}()
    for g in collect(gb)
        iszero(g) && continue
        es = [collect(Int, e) for e in Nemo.exponent_vectors(g)]
        push!(out, reduce((a, b) -> _drl_less(a, b) ? b : a, es))
    end
    out
end

"chi of a (specialized, mod-p) ideal: f4 GB + staircase count."
function chi_of_ideal(I, R)
    gb = Oscar.groebner_basis_f4(I)
    staircase_count(lead_exps(gb, R), nvars(R))
end
