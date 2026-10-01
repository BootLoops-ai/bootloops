# Unified numerical evaluator for regularized iterated integrals
#
#   I(ω_{w₁}, …, ω_{wₙ}; x) = ∫₀^x ω_{w₁}(λ) ∫₀^λ ω_{w₂}(λ') ⋯
#
# where each kernel one-form has at most a simple pole at the base point 0:
#   ω = ( pole/λ + Σ_{n≥0} taylor[n+1] λⁿ ) dλ .
# Divergences at the base point are regularized by the tangential-base-point /
# shuffle prescription: every primitive is taken with zero constant in the
# log-graded representation  F(λ) = Σ_k p_k(λ) Lᵏ, L = log λ.
# This reproduces simultaneously
#   * GPL trailing/leading-zero regularization: G(0;x) = log x, G(0,0;x) = log²x/2!,
#   * the modular convention I(1;q) = log q of arXiv:1704.08895.
#
# First letter of a word = OUTERMOST integration (G and I conventions).
#
# Accuracy: series truncated at `nterms`; an error ball is added from a
# geometric tail estimate (heuristic safety factor; certified evaluation is
# the job of the W4 kernel library Eichler.jl — documented engine limit).

export KernelForm, DlogKernel, SeriesKernel, ItIntEvaluator, eval_word, eval_words

abstract type KernelForm end

"dλ/(λ − a); a = 0 means the dlog-at-base letter dλ/λ."
struct DlogKernel <: KernelForm
    a::Acb
end

"(pole/λ + Σ taylor[n+1] λⁿ) dλ. For a modular form f, f(q)dq/q has pole = a₀(f)."
struct SeriesKernel <: KernelForm
    pole::Acb
    taylor::Vector{Acb}
end

# F(λ) = Σ_{k} blocks[k+1][n+1] λⁿ Lᵏ
struct LogSeries
    blocks::Vector{Vector{Acb}}
end

const ONE_LOGSERIES_KEY = Int[]

mutable struct ItIntEvaluator
    kernels::Vector{KernelForm}
    nterms::Int
    prec::Int
    cache::Dict{Word,LogSeries}   # suffix → primitive of that suffix word
end

function ItIntEvaluator(kernels::Vector{<:KernelForm}; nterms::Int, prec::Int)
    ItIntEvaluator(collect(KernelForm, kernels), nterms, prec,
                   Dict{Word,LogSeries}())
end

azero(prec) = Acb(0, prec = prec)

"upper bound of an Arb as a (point) Arb"
ub(x::Arb) = Arb(Arblib.ubound(Arf, x), prec = 64)
"Mag upper bound of an Arb"
magof(x::Arb) = Arblib.Mag(Arblib.ubound(Arf, x))

function one_logseries(N, prec)
    b = [azero(prec) for _ in 1:N]
    b[1] = Acb(1, prec = prec)
    LogSeries([b])
end

# --- kernel application: H(λ) = ω(λ)/dλ · F(λ), split as (pole·F/λ , analytic) ---

# returns (pole_coeff, analytic_blocks) where integrand = pole_coeff·F/λ + analytic
function apply_kernel(K::DlogKernel, F::LogSeries, N, prec)
    if iszero(K.a)
        return (Acb(1, prec = prec), nothing, F)
    end
    # 1/(λ−a)·F : per log-block recurrence  h_n = (h_{n−1} − f_n)/a
    blocks = Vector{Vector{Acb}}(undef, length(F.blocks))
    for (k, f) in enumerate(F.blocks)
        h = [azero(prec) for _ in 1:N]
        h[1] = -f[1] / K.a
        for n in 2:N
            h[n] = (h[n-1] - f[n]) / K.a
        end
        blocks[k] = h
    end
    return (azero(prec), LogSeries(blocks), F)
end

function apply_kernel(K::SeriesKernel, F::LogSeries, N, prec)
    blocks = Vector{Vector{Acb}}(undef, length(F.blocks))
    T = K.taylor
    for (k, f) in enumerate(F.blocks)
        h = [azero(prec) for _ in 1:N]
        for n in 0:N-1                      # h_n = Σ_{m≤n} T_m f_{n−m}
            s = azero(prec)
            for m in 0:min(n, length(T) - 1)
                Arblib.addmul!(s, T[m+1], f[n-m+1])
            end
            h[n+1] = s
        end
        blocks[k] = h
    end
    return (K.pole, LogSeries(blocks), F)
end

# --- integration ∫₀^λ with tangential-base-point regularization ---

# ∫₀^λ sⁿ⁻¹ logᵏ(s) ds = λⁿ Σ_{j=0}^{k} (−1)^{k−j} (k!/j!) n^{−(k−j+1)} logʲλ , n ≥ 1
# ∫₀^λ s⁻¹  logᵏ(s) ds = logᵏ⁺¹(λ)/(k+1)   (regularized: zero constant)
function integrate_logseries(pole_coeff::Acb, analytic::Union{LogSeries,Nothing},
                             F::LogSeries, N, prec)
    Kmax = length(F.blocks) + 1   # pole part can raise log degree by one
    out = [[azero(prec) for _ in 1:N] for _ in 1:Kmax]
    # pole part: pole_coeff · Σ_k Σ_n f_{k,n} λⁿ⁻¹ Lᵏ
    if !iszero(pole_coeff)
        for (kk, f) in enumerate(F.blocks)
            k = kk - 1
            out[kk+1][1] += pole_coeff * f[1] / (k + 1)   # n = 0 → Lᵏ⁺¹/(k+1)
            for n in 1:N-1
                c = pole_coeff * f[n+1]
                iszero(c) && continue
                for j in 0:k
                    coef = Acb((-1)^(k - j) * (factorial(big(k)) ÷ factorial(big(j))),
                               prec = prec) * inv(Acb(n, prec = prec))^(k - j + 1)
                    out[j+1][n+1] += c * coef
                end
            end
        end
    end
    # analytic part: Σ a_{k,n} λⁿ Lᵏ → ∫ uses the same formula with n → n+1
    if analytic !== nothing
        for (kk, a) in enumerate(analytic.blocks)
            k = kk - 1
            for n in 0:N-2
                c = a[n+1]
                iszero(c) && continue
                np = n + 1
                for j in 0:k
                    coef = Acb((-1)^(k - j) * factorial(big(k)) ÷ factorial(big(j)),
                               prec = prec) * inv(Acb(np, prec = prec))^(k - j + 1)
                    out[j+1][np+1] += c * coef
                end
            end
        end
    end
    # drop trailing all-zero log blocks
    kk = Kmax
    while kk > 1 && all(iszero, out[kk])
        kk -= 1
    end
    LogSeries(out[1:kk])
end

# --- word evaluation with suffix caching ---

function word_primitive(ev::ItIntEvaluator, w::Word)
    haskey(ev.cache, w) && return ev.cache[w]
    F = isempty(w[2:end]) ? one_logseries(ev.nterms, ev.prec) :
                            word_primitive(ev, w[2:end])
    pole, analytic, _ = apply_kernel(ev.kernels[w[1]], F, ev.nterms, ev.prec)
    G = integrate_logseries(pole, analytic, F, ev.nterms, ev.prec)
    ev.cache[w] = G
    G
end

"""
    eval_word(ev, word, x) -> Acb

Value of the regularized iterated integral of `word` at `x` (complex ball),
with a heuristic geometric tail bound added to the error radius.
"""
function eval_word(ev::ItIntEvaluator, w::Word, x::Acb)
    F = word_primitive(ev, w)
    L = log(x)
    N = ev.nterms
    total = azero(ev.prec)
    tail = Arb(0, prec = 64)
    ax = abs(x)
    for (kk, blk) in enumerate(F.blocks)
        # Horner in x
        s = azero(ev.prec)
        for n in N:-1:1
            Arblib.mul!(s, s, x)
            Arblib.add!(s, s, blk[n])
        end
        total += s * L^(kk - 1)
        # Tail estimate (heuristic, reviewed): coefficient magnitude from a
        # trailing WINDOW (≥13 wide — covers arithmetic-progression-supported
        # q-series like E₄(6τ)) with an additive floor so it is never exactly
        # zero, and geometric ratio from the EMPIRICAL coefficient growth
        # (≈ 1/min|aᵢ| for dlog letters), not from |x| alone.
        lo = max(1, N - 12)
        cN = Arb(0, prec = 64)
        for n in lo:N
            cN = max(cN, ub(abs(blk[n])))
        end
        cW = Arb(0, prec = 64)   # window lower in the series, for growth ratio
        for n in max(1, lo - 13):lo-1
            cW = max(cW, ub(abs(blk[n])))
        end
        growth = (cW > 0 && cN > 0) ? max(ub(cN / cW)^(Arb(1) / 13), Arb(1)) : Arb(1)
        rho = ub(ax) * growth
        if !(rho < 1)
            # outside (or too close to) the empirical convergence radius:
            # return an honest infinite-radius ball rather than a wrong one
            tail += Arb(Inf)
            continue
        end
        # tail ≈ Σ_{n>N} c_n xⁿ ≲ c_N·|x|^N · Σ_{m≥1} (growth·|x|)^m
        t = (cN + Arb(1e-300)) * ub(ax)^N * rho / (1 - rho) * 16 *
            ub(abs(L))^(kk - 1)
        tail += ub(t)
    end
    Arblib.add_error!(total, magof(tail))
    total
end

"""
    eval_words(ev, words, x) -> Vector{Acb}
"""
eval_words(ev::ItIntEvaluator, ws::Vector{Word}, x::Acb) =
    [eval_word(ev, w, x) for w in ws]
