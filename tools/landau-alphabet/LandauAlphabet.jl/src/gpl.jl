# GPL / MPL kernel instance.
#
# G(a₁,…,aₙ; x) (Goncharov convention, a₁ outermost) evaluated through the
# unified iterated-integral evaluator with DlogKernel letters. Trailing zeros
# are automatically shuffle-regularized by the log-graded base-point
# prescription (G(0,…,0;x) = logⁿx/n!).
#
# Convergence: series at the base point requires |x| < min_{aᵢ≠0} |aᵢ|.
# Callers must choose sample points accordingly.

export GPLContext, gpl, gpl_terms_needed, const_value, MPLBasisFunction,
       basis_tower, basis_value, basis_label

struct GPLContext
    indices::Vector{Acb}            # distinct index values; letter i ↔ indices[i]
    ev::ItIntEvaluator
end

function GPLContext(indices::Vector{Acb}; nterms::Int, prec::Int)
    kernels = KernelForm[DlogKernel(a) for a in indices]
    GPLContext(indices, ItIntEvaluator(kernels; nterms = nterms, prec = prec))
end

"number of series terms for d digits at |x|/r ratio ρ"
gpl_terms_needed(rho::Float64, digits::Int) =
    rho >= 0.999 ? error("x too close to convergence boundary") :
    ceil(Int, digits * log(10) / log(1 / rho)) + 50

"""
    gpl(ctx, word_indices, x) -> Acb

Evaluate G(a⃗; x) where `word_indices` are positions into ctx.indices.
"""
gpl(ctx::GPLContext, w::Word, x::Acb) = eval_word(ctx.ev, w, x)

# ---------- transcendental-constant bookkeeping for the basis towers ----------

const CONST_WEIGHT = Dict(
    :one => 0, :pi => 1, :zeta2 => 2, :zeta3 => 3, :pi3 => 3, :pi2 => 2,
    :zeta4 => 4, :pi4 => 4, :zeta5 => 5, :zeta2zeta3 => 5, :pi5 => 5,
    :zeta6 => 6, :zeta3sq => 6)

function const_value(s::Symbol, prec::Int)
    P = Arb(π, prec = prec)
    s === :one && return Acb(1, prec = prec)
    s === :pi && return Acb(P, prec = prec)
    s === :pi2 && return Acb(P^2, prec = prec)
    s === :pi3 && return Acb(P^3, prec = prec)
    s === :pi4 && return Acb(P^4, prec = prec)
    s === :pi5 && return Acb(P^5, prec = prec)
    s === :zeta2 && return Acb(P^2 / 6, prec = prec)
    s === :zeta3 && return Acb(Arblib.zeta!(Arb(prec = prec), Arb(3, prec = prec)), prec = prec)
    s === :zeta4 && return Acb(P^4 / 90, prec = prec)
    s === :zeta5 && return Acb(Arblib.zeta!(Arb(prec = prec), Arb(5, prec = prec)), prec = prec)
    s === :zeta6 && return Acb(P^6 / 945, prec = prec)
    s === :zeta2zeta3 && return const_value(:zeta2, prec) * const_value(:zeta3, prec)
    s === :zeta3sq && return const_value(:zeta3, prec)^2
    error("unknown constant $s")
end

"""
A basis function  c × G(word; x)  with c a transcendental constant symbol.
`word` may be empty (pure constant).
"""
struct MPLBasisFunction
    csym::Symbol
    word::Word
end

basis_label(ctx::GPLContext, b::MPLBasisFunction) =
    (b.csym === :one ? "" : String(b.csym) * "*") *
    (isempty(b.word) ? "1" : "G(" * join([string(ctx.indices[i]) for i in b.word], ",") * ";x)")

basis_value(ctx::GPLContext, b::MPLBasisFunction, x::Acb, prec::Int) =
    const_value(b.csym, prec) * (isempty(b.word) ? Acb(1, prec = prec) : gpl(ctx, b.word, x))

"""
    basis_tower(nletters, weight; variant=:full) -> Vector{MPLBasisFunction}

Weight-graded MPL function basis for a one-variable alphabet with `nletters`
indices — each weight line pairs the length-w G-words with lower-length words
multiplied by the transcendental constants that restore the weight:
  * :unif — the uniform weight-`weight` line,
  * :full — union of :unif lines for weights 1..weight,
  * :sim  — like :full but dropping the odd-π terms.
Conventions: the weight-3 line uses π²·G (≡ ζ₂·G) and includes the pure
constant π³ (58 functions for 3 letters at weight ≤ 3, :full).
"""
function basis_tower(nletters::Int, weight::Int; variant::Symbol = :full)
    if variant === :full
        return reduce(vcat, [basis_tower(nletters, w; variant = :unif) for w in 1:weight])
    elseif variant === :sim
        return filter(b -> !(b.csym in (:pi, :pi3, :pi5)),
                      basis_tower(nletters, weight; variant = :full))
    end
    variant === :unif || error("variant must be :full, :sim or :unif")
    allwords(k) = k == 0 ? [Int[]] : vec([vcat(w, i) for w in allwords(k - 1), i in 1:nletters])
    out = MPLBasisFunction[]
    line = Dict(
        1 => [(:one, 1), (:pi, 0)],
        2 => [(:one, 2), (:pi, 1), (:zeta2, 0)],
        3 => [(:one, 3), (:pi, 2), (:pi2, 1), (:zeta3, 0), (:pi3, 0)],
        4 => [(:one, 4), (:pi, 3), (:zeta2, 2), (:pi3, 1), (:zeta3, 1), (:zeta4, 0)],
        5 => [(:one, 5), (:pi, 4), (:zeta2, 3), (:pi3, 2), (:zeta3, 2), (:zeta4, 1),
              (:zeta5, 0), (:zeta2zeta3, 0)],
        6 => [(:one, 6), (:pi, 5), (:zeta2, 4), (:pi3, 3), (:zeta3, 3), (:zeta4, 2),
              (:zeta5, 1), (:zeta2zeta3, 1), (:pi5, 1), (:zeta6, 0), (:zeta3sq, 0)],
    )
    haskey(line, weight) || error("weight $weight not tabulated")
    for (c, gw) in line[weight], w in allwords(gw)
        push!(out, MPLBasisFunction(c, w))
    end
    out
end
