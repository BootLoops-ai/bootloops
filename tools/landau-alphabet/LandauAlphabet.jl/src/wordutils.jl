# Words over a kernel alphabet: enumeration by weight, shuffle algebra, Lyndon words.
#
# A word is a Vector{Int} of letter ids (first letter = outermost integration,
# matching the standard G(a₁,…,aₙ;x) convention and the iterated-integral
# I(f₁,…,fₙ;q) convention of arXiv:1704.08895).

export Word, words_of_weight, shuffle_product, lyndon_words, shuffle_expand_in_lyndon

const Word = Vector{Int}

"""
    words_of_weight(letter_weights::Vector{Int}, w::Int) -> Vector{Word}

All words over letters `1:n` (with `letter_weights[i]` the weight of letter i)
whose total weight is exactly `w`. Letters of weight 0 are not allowed (use
weight ≥ 1; the modular letter `1` carries weight 1 in the τ-grading).
"""
function words_of_weight(letter_weights::Vector{Int}, w::Int)
    out = Word[]
    stack = [(Int[], 0)]
    while !isempty(stack)
        (pre, s) = pop!(stack)
        for (i, wi) in enumerate(letter_weights)
            wi <= 0 && error("letter weights must be positive")
            s2 = s + wi
            if s2 == w
                push!(out, vcat(pre, i))
            elseif s2 < w
                push!(stack, (vcat(pre, i), s2))
            end
        end
    end
    sort!(out)
    out
end

"""
    shuffle_product(u::Word, v::Word) -> Dict{Word,Int}

Shuffle product u ⧢ v as a multiset of words with integer multiplicities.
"""
function shuffle_product(u::Word, v::Word)
    out = Dict{Word,Int}()
    _shuffle!(out, Int[], u, 1, v, 1)
    out
end

function _shuffle!(out, pre, u, i, v, j)
    if i > length(u) && j > length(v)
        out[copy(pre)] = get(out, copy(pre), 0) + 1
        return
    end
    if i <= length(u)
        push!(pre, u[i]); _shuffle!(out, pre, u, i + 1, v, j); pop!(pre)
    end
    if j <= length(v)
        push!(pre, v[j]); _shuffle!(out, pre, u, i, v, j + 1); pop!(pre)
    end
end

"""
    lyndon_words(nletters, maxlen) -> Vector{Word}

All Lyndon words over the alphabet 1:nletters of length ≤ maxlen (Duval).
"""
function lyndon_words(nletters::Int, maxlen::Int)
    out = Word[]
    w = [1]
    while true
        length(w) <= maxlen && push!(out, copy(w))
        # extend periodically to maxlen, then increment
        x = copy(w)
        while length(x) < maxlen
            push!(x, x[length(x) % length(w) + 1])
        end
        while !isempty(x) && x[end] == nletters
            pop!(x)
        end
        isempty(x) && break
        x[end] += 1
        w = x
    end
    sort!(out, by = x -> (length(x), x))
    out
end

"""
    shuffle_expand_in_lyndon(word, lyndon_index) -> Dict{Vector{Word},Rational{BigInt}}

Not needed for the oracle suite; placeholder kept minimal. The constraint layer
imposes shuffle relations directly via `shuffle_product`.
"""
function shuffle_expand_in_lyndon end
