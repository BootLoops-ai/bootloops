# Integrable symbols over Q for multivariate dlog alphabets.
#
# An alphabet is given by irreducible polynomial FACTORS F₁…F_f in
# Q[x₁,…,x_v]; a LETTER is a multiplicative word over the factors (an integer
# exponent vector), so ratios like z̄/z are first-class letters.
#
# The symbol  Σ c_w L_{w₁} ⊗ ⋯ ⊗ L_{wₙ}  is integrable iff for every slot k
# and every environment (the other n−2 slots fixed) the adjacent wedge
# Σ c · dlog L_{w_k} ∧ dlog L_{w_{k+1}} vanishes as a rational two-form.
# All conditions are assembled exactly over Q.

export SymbolAlphabet, symbol_alphabet, integrable_symbols

struct SymbolAlphabet
    ring::QQMPolyRing
    factors::Vector{<:MPolyRingElem}
    letters::Vector{Vector{Int}}        # exponent vectors over factors
    letter_names::Vector{String}
end

"""
    symbol_alphabet(varnames, factor_strs, letters) -> SymbolAlphabet

`factor_strs` are polynomial expressions in the variables; `letters` is a
vector of (name, exponent-vector) pairs over those factors.
"""
function symbol_alphabet(varnames::Vector{String}, factor_strs::Vector{String},
                         letters::Vector{Tuple{String,Vector{Int}}})
    R, _ = polynomial_ring(QQ, varnames)
    facs = [eval_poly_expr(Meta.parse(p), R) for p in factor_strs]
    SymbolAlphabet(R, facs, [l[2] for l in letters], [l[1] for l in letters])
end

function eval_poly_expr(ex, R::QQMPolyRing)
    vars = Dict(string(s) => g for (s, g) in zip(symbols(R), gens(R)))
    _ev(e) = e isa Integer ? R(e) :
             e isa Symbol ? vars[string(e)] :
             e isa Expr && e.head == :call ?
                (e.args[1] == :+ ? reduce(+, map(_ev, e.args[2:end])) :
                 e.args[1] == :- ? (length(e.args) == 2 ? -_ev(e.args[2]) :
                                    _ev(e.args[2]) - _ev(e.args[3])) :
                 e.args[1] == :* ? reduce(*, map(_ev, e.args[2:end])) :
                 e.args[1] == :^ ? _ev(e.args[2])^e.args[3] :
                 error("op $(e.args[1])")) :
             error("bad expr $e")
    _ev(ex)
end

# Coordinates φ(i,j) ∈ Q^r of the factor wedges dlogFᵢ ∧ dlogFⱼ such that
# Q-relations among wedges ⟺ relations among coordinates.
function factor_wedge_coords(alpha::SymbolAlphabet)
    R = alpha.ring
    F = alpha.factors
    f = length(F)
    v = nvars(R)
    pairs = [(i, j) for i in 1:f for j in i+1:f]
    D = prod(F)^2
    # numerators of (dlogFᵢ ∧ dlogFⱼ)·D per dx_a∧dx_b component
    mons = Dict{Tuple{Int,Vector{Int}},Int}()
    polyvecs = Vector{Vector{elem_type(R)}}()
    for (i, j) in pairs
        comps = elem_type(R)[]
        ci = 0
        for a in 1:v, b in a+1:v
            ci += 1
            num = derivative(F[i], a) * derivative(F[j], b) -
                  derivative(F[i], b) * derivative(F[j], a)
            p = divexact(num * D, F[i] * F[j])
            push!(comps, p)
            for ex in exponent_vectors(p)
                haskey(mons, (ci, ex)) || (mons[(ci, ex)] = length(mons) + 1)
            end
        end
        push!(polyvecs, comps)
    end
    M = zero_matrix(QQ, length(pairs), max(length(mons), 1))
    for (r, pv) in enumerate(polyvecs)
        ci = 0
        for p in pv
            ci += 1
            for (cf, ex) in zip(coefficients(p), exponent_vectors(p))
                M[r, mons[(ci, ex)]] = cf
            end
        end
    end
    rk, RM = rref(M)
    # pivot columns of the rref → coordinate projection preserving row relations
    piv = Int[]
    for i in 1:rk
        j = findfirst(!iszero, [RM[i, c] for c in 1:ncols(RM)])
        push!(piv, j)
    end
    coords = Dict{Tuple{Int,Int},Vector{QQFieldElem}}()
    for (idx, (i, j)) in enumerate(pairs)
        coords[(i, j)] = [M[idx, c] for c in piv]
    end
    (coords, rk)
end

# wedge coordinates for letter pair (a,b)
function letter_wedge(alpha::SymbolAlphabet, coords, r::Int, a::Int, b::Int)
    ea, eb = alpha.letters[a], alpha.letters[b]
    out = [QQ(0) for _ in 1:r]
    f = length(alpha.factors)
    for i in 1:f, j in i+1:f
        kappa = ea[i] * eb[j] - ea[j] * eb[i]
        kappa == 0 && continue
        c = coords[(i, j)]
        for k in 1:r
            out[k] += kappa * c[k]
        end
    end
    out
end

"""
    integrable_symbols(alpha, w) -> (words, nullbasis::QQMatrix)

All weight-w integrable symbols over the alphabet. Columns of `nullbasis` are
coefficient vectors (indexed like `words`) spanning the integrable subspace.
"""
function integrable_symbols(alpha::SymbolAlphabet, w::Int)
    n = length(alpha.letters)
    # the per-environment splitting of the integrability condition (and the
    # very meaning of "basis of integrable symbols") requires the letters'
    # dlogs to be Q-linearly independent — i.e. the exponent matrix over the
    # (distinct irreducible) factors must have full row rank
    Emat = matrix(QQ, n, length(alpha.factors),
                  [QQ(alpha.letters[l][i]) for l in 1:n for i in 1:length(alpha.factors)])
    rank(Emat) == n || error("letters are multiplicatively dependent (exponent rank " *
                             "$(rank(Emat)) < $n); drop dependent letters or re-express the alphabet")
    words = words_of_weight(fill(1, n), w)
    w == 1 && return (words, identity_matrix(QQ, n))
    coords, r = factor_wedge_coords(alpha)
    W = Dict((a, b) => letter_wedge(alpha, coords, r, a, b) for a in 1:n for b in 1:n if a < b)
    widx = Dict(wd => i for (i, wd) in enumerate(words))
    envs = w == 2 ? [Int[]] : words_of_weight(fill(1, n), w - 2)
    nrows = (w - 1) * length(envs) * r
    A = zero_matrix(QQ, nrows, length(words))
    row = 0
    for k in 1:w-1, env in envs
        for comp in 1:r
            row += 1
            for a in 1:n, b in a+1:n
                cf = W[(a, b)][comp]
                iszero(cf) && continue
                wab = vcat(env[1:k-1], [a, b], env[k:end])
                wba = vcat(env[1:k-1], [b, a], env[k:end])
                A[row, widx[wab]] += cf
                A[row, widx[wba]] -= cf
            end
        end
    end
    nd, ns = nullspace(A)
    (words, ns)
end
