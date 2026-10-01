# Symbol-word and basis evaluation.
#
# A symbol word [s₁,…,sₙ] over a dlog alphabet {Lₛ} integrates to
#   ∫ dlog Lₛ₁ ⋯ dlog Lₛₙ  =  G(a₁,…,aₙ; z)
# when each letter is linear, Lₛ(z) ∝ (z − aₛ). The alphabet dictionary maps
# each letter symbol to a function `point -> aₛ::Number` (the GPL index, i.e.
# the root of the letter, possibly depending on the OTHER kinematic variables).
# `z` is the integration-variable value at `point`.
#
# This is the "symbol known → number" step. Batch evaluation reuses a single
# cache; for very large alphabets at a fixed point the LandauAlphabet
# `GPLContext` (suffix-cached series) is also provided as the fast path.

"""
    evaluate_symbol(word, alphabet, z, point; prec, cache) -> Acb

Evaluate the iterated dlog integral of `word` (a `Vector{Symbol}`) at the
upper endpoint `z`, with GPL indices `aᵢ = alphabet[word[i]](point)`.
"""
function evaluate_symbol(word::AbstractVector{Symbol},
                         alphabet::Dict{Symbol,<:Function},
                         z, point; prec::Int,
                         cache::Union{Nothing,GPLCache} = nothing)
    a = [alphabet[s](point) for s in word]
    gpl(a, z; prec = prec, cache = cache)
end

"""
    evaluate_basis(words, alphabet, z, point; prec) -> Vector{Acb}

Batch-evaluate a list of symbol words. Two paths:
  * default: a shared `GPLCache` (handles |z| ≥ min|aᵢ| via Hölder);
  * `series_only=true`: build one `LandauAlphabet.GPLContext` (suffix-cached
    LogSeries; fastest when |z| < min|aᵢ| and the basis is large).
"""
function evaluate_basis(words::AbstractVector,
                        alphabet::Dict{Symbol,<:Function},
                        z, point; prec::Int, series_only::Bool = false,
                        nterms::Int = 0)
    syms = sort!(collect(keys(alphabet)))
    idx  = Dict(s => i for (i, s) in enumerate(syms))
    avals = Acb[Acb(alphabet[s](point), prec = prec) for s in syms]
    zc = Acb(z, prec = prec)
    if series_only
        # fast LandauAlphabet path
        r = minimum(_lbf(abs(a)) for a in avals if !iszero(a); init = Inf)
        nt = nterms > 0 ? nterms :
             LandauAlphabet.gpl_terms_needed(_ubf(abs(zc)) / r,
                                              ceil(Int, prec * log10(2.0)))
        ctx = GPLContext(avals; nterms = nt, prec = prec)
        return Acb[LandauAlphabet.gpl(ctx, Int[idx[s] for s in w], zc) for w in words]
    else
        cache = GPLCache(prec)
        return Acb[evaluate_symbol(collect(Symbol, w), alphabet, z, point;
                                   prec = prec, cache = cache) for w in words]
    end
end
