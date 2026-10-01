# Harmonic polylogarithms H(m⃗; z) (Remiddi–Vermaseren), via gpl.
#
# Convention matches the project's hand-rolled hpl_fast.py / hpl_basis.py:
#   df₀ = dz/z,  df₁ = dz/(1−z),  df₋₁ = dz/(1+z),
# i.e.  H(a⃗; z) = (−1)ᵖ G(a⃗; z) with p = #(aᵢ = +1) and aᵢ ∈ {0,±1}.
# Compressed m-vector: mᵢ > 0 ↦ 0^{mᵢ−1},1 ;  mᵢ < 0 ↦ 0^{|mᵢ|−1},−1 ;
# mᵢ = 0 is a literal 0-letter (so H(0,…,0;z) = logⁿz/n!).

"Expand a compressed m-vector (e.g. [2,1] ↦ [0,1,1], [-2] ↦ [0,-1])."
function hpl_expand(m::AbstractVector{<:Integer})
    a = Int[]
    for mi in m
        if mi == 0
            push!(a, 0)
        else
            append!(a, zeros(Int, abs(mi) - 1))
            push!(a, sign(mi))
        end
    end
    a
end

"""
    hpl(m, z; prec, cache=nothing) -> Acb

Harmonic polylog H(m⃗; z). `m` may be either the compressed weight vector
(e.g. `[2,1]`) or the expanded letter word over `{0,1,-1}` — disambiguated by
`expanded=true`.
"""
function hpl(m::AbstractVector{<:Integer}, z; prec::Int,
             cache::Union{Nothing,GPLCache} = nothing, expanded::Bool = false)
    a = expanded ? collect(Int, m) : hpl_expand(m)
    p = count(==(1), a)
    return (isodd(p) ? -1 : 1) * gpl(a, z; prec = prec, cache = cache)
end
