# The equal-mass two-loop sunrise to all orders in ε through the Γ₁(6) modular
# representation (docs/conventions.md §5–6; AW's all-orders result). Units m² = μ² = 1, t = p²/m².

"""Boundary constants B^(k), k = 0..J-1 (ε-series, lead 0). Conventions §5."""
function sunrise_boundary_constants(prec::Int, J::Int)
    Jw = J + 3
    F3 = hyp2f1_at_r3(prec, Jw)                      # ε-poly of ₂F₁ at r₃
    F3c = [conj(c) for c in F3]                       # at r₃⁻¹ (real ε-Taylor in x)
    i_ = Acb(0, 1; prec)
    piA = Acb(π; prec)
    # e^{±iπε/3}
    eplus = eps_exp(EpsSeries(1, epoly_const(prec, Jw - 1, i_ * piA / 3)))
    eminus = eps_exp(EpsSeries(1, epoly_const(prec, Jw - 1, -i_ * piA / 3)))
    F3s = EpsSeries(0, F3); F3cs = EpsSeries(0, F3c)
    h = (eplus * F3s - eminus * F3cs) * eps_const(prec, Jw, 1 / i_)
    # S₁₁₁(2-2ε,0) = e^{2γε}Γ(1+2ε)(√3)^{-1-2ε}[3/(2ε²)Γ(1+ε)²/Γ(1+2ε)h - π/ε]
    γ = Arb(prec = prec); Arblib.const_euler!(γ)
    e2γ = eps_exp(EpsSeries(1, epoly_const(prec, Jw - 1, Acb(2 * γ; prec))))
    Γ1 = eps_gamma1p(prec, Jw, 1)
    Γ2 = eps_gamma1p(prec, Jw, 2)
    s3 = sqrt(Acb(3; prec))
    pw = eps_pow(prec, Jw, s3, -1, -2)               # (√3)^{-1-2ε}
    bracket = eps_shift(eps_const(prec, Jw, Acb(3//2; prec)) * Γ1 * Γ1 * eps_inv(Γ2) * h, -2) -
              eps_shift(eps_const(prec, Jw, piA), -1)
    S0 = e2γ * Γ2 * pw * bracket
    # The printed boundary equation gives the γ-stripped series Σ εʲ S^(j)(2,0)
    # (AW convention S = e^{-2γε} Σ εʲ S^(j)); the actual integral is e^{-2γε}·S0.
    # Σ ε^k B^(k) = (√3/2)·Γ(1+ε)^{-2}·e^{-2γε}·S0     [ψ₁(0)/π = 2/√3, L = 0]
    e2γm = eps_exp(EpsSeries(1, epoly_const(prec, Jw - 1, Acb(-2 * γ; prec))))
    Bser = eps_const(prec, Jw, s3 / 2) * eps_inv(Γ1 * Γ1) * e2γm * S0
    # the ε^{<0} coefficients must vanish — verify and strip
    for k in Bser.lead:-1
        Arblib.contains_zero(eps_coeff(Bser, k)) ||
            error("sunrise boundary: nonvanishing pole coefficient at ε^$k")
    end
    return [eps_coeff(Bser, k) for k in 0:J-1]
end

# ---- word machinery -----------------------------------------------------------

_kernel_series(G::Gamma16, s::Symbol) =
    s === :one ? qseries_const(G.ctx, 1) :
    s === :f2 ? G.f2 : s === :f3 ? G.f3 : s === :f4 ? G.f4 :
    error("unknown kernel $s")

"""Value (and value of the once-peeled word) of I(word...; q). Returns (I, dI) where
dI = q d/dq I = f₁(q)·I(rest; q)."""
function _word_value(G::Gamma16, word::Vector{Symbol}, q::Acb, logq::Acb,
                     cache::Dict{Vector{Symbol},Acb})
    isempty(word) && return (Acb(1; prec = G.ctx.prec), Acb(0; prec = G.ctx.prec))
    val = get!(cache, word) do
        iterated_integral_value([_kernel_series(G, s) for s in word], q, logq)
    end
    rest = word[2:end]
    restval = get!(cache, rest) do
        isempty(rest) ? Acb(1; prec = G.ctx.prec) :
            iterated_integral_value([_kernel_series(G, s) for s in rest], q, logq)
    end
    f1q = evaluate(_kernel_series(G, word[1]), q)
    return (val, f1q * restval)
end

_hword(j::Int) = reduce(vcat, [[:one, :f4] for _ in 1:j]; init = Symbol[])

"""All-orders sunrise S₁₁₁(2-2ε, t) (and its q-logarithmic derivative) at a kinematic
point t (ball, in units m²=1, Euclidean-region conventions). Returns a named tuple
(S, dSdlogq, q, logq, tauC, psi1) of EpsSeries/Acb with J stored ε-orders."""
function sunrise2(G::Gamma16, t::Acb; J::Int = 8,
                  B::Union{Nothing,Vector{Acb}} = nothing)
    prec = G.ctx.prec
    q, logq, tauC = nome_and_log_from_t(G, Acb(t; prec = prec))
    # independent self-check: certified Hauptmodul (eta-quotient route) at τ must contain t
    thm = hauptmodul_point(tauC)
    Arblib.overlaps(thm, Acb(t; prec = prec)) ||
        error("sunrise2: Hauptmodul(τ(t)) does not overlap t — branch failure")
    B === nothing && (B = sunrise_boundary_constants(prec, J))
    cache = Dict{Vector{Symbol},Acb}()
    # prefactor pieces
    If2, dIf2 = _word_value(G, [:f2], q, logq, cache)
    Γ1 = eps_gamma1p(prec, J, 1)
    expfac = eps_exp(EpsSeries(1, epoly_const(prec, J - 1, -If2)))   # e^{-ε I(f2;q)}
    P = Γ1 * Γ1 * expfac
    dP_dlq = EpsSeries(1, epoly_const(prec, J - 1, -evaluate(G.f2, q))) * P  # qd/dq P = -ε f2(q) P
    # F1 = Σ_j ε^{2j} I({1,f4}^j) - ½ ε^{2j+1} I({1,f4}^j,1)
    F1 = eps_zero(prec, J); dF1 = eps_zero(prec, J)
    Bs = eps_zero(prec, J)
    for k in 0:J-1
        Bs += eps_shift(eps_const(prec, J, B[k+1]), k)
    end
    for j in 0:div(J, 2)
        w = _hword(j)
        v, dv = _word_value(G, w, q, logq, cache)
        w1 = vcat(w, [:one])
        v1, dv1 = _word_value(G, w1, q, logq, cache)
        2j <= J - 1 && (F1 += eps_shift(eps_const(prec, J, v), 2j);
                        dF1 += eps_shift(eps_const(prec, J, dv), 2j))
        2j + 1 <= J - 1 && (F1 -= eps_shift(eps_const(prec, J, v1 / 2), 2j + 1);
                            dF1 -= eps_shift(eps_const(prec, J, dv1 / 2), 2j + 1))
    end
    # F3 = Σ_j ε^j Σ_{k≤j/2} I({1,f4}^k, 1, f3, {f2}^{j-2k})
    F3 = eps_zero(prec, J); dF3 = eps_zero(prec, J)
    for j in 0:J-1, k in 0:div(j, 2)
        w = vcat(_hword(k), [:one, :f3], fill(:f2, j - 2k))
        v, dv = _word_value(G, w, q, logq, cache)
        F3 += eps_shift(eps_const(prec, J, v), j)
        dF3 += eps_shift(eps_const(prec, J, dv), j)
    end
    R = F1 * Bs + F3
    dR = dF1 * Bs + dF3
    ψ = evaluate(G.psi1, q)
    dψ = evaluate(qdq(G.psi1), q)
    Pψ = eps_const(prec, J, ψ) * P
    S = Pψ * R
    dS = (eps_const(prec, J, dψ) * P + eps_const(prec, J, ψ) * dP_dlq) * R + Pψ * dR
    return (S = S, dSdlogq = dS, q = q, logq = logq, tauC = tauC, psi1 = ψ)
end

"""S₁₁₁(4-2ε, t) via the equal-mass dimension shift (conventions §6). Returns an
EpsSeries with lead = -2. AMFlow convention sunrise (k²-m²+i0 propagators) = -S₁₁₁."""
function sunrise4(G::Gamma16, t::Acb; J::Int = 10)
    prec = G.ctx.prec
    res2 = sunrise2(G, t; J = J)
    # dS/dt = (q dS/dq) / (q dt/dq)
    qdt = evaluate(qdq(G.t), res2.q)
    dSdt = res2.dSdlogq * eps_const(prec, J, inv(qdt))
    tA = Acb(t; prec = prec)
    A = (tA + 3) * (tA - 1) * (tA - 9)
    B0 = (tA - 1) * (tA - 9)
    B1 = tA^2 + 22 * tA - 87
    # T₁(4-2ε)² = [Γ(1+ε)/(ε(ε-1))]², lead -2
    Γ1 = eps_gamma1p(prec, J, 1)
    εm1 = eps_const(prec, J, -1) + eps_eps(prec, J)     # ε - 1
    T1 = eps_shift(Γ1 * eps_inv(εm1), -1)
    T1sq = T1 * T1
    tad = eps_const(prec, J, Acb(-6; prec)) + eps_shift(eps_const(prec, J, tA + 21), 1)
    one_ = eps_const(prec, J, 1)
    pref = eps_inv(eps_const(prec, J, 6) *
                   (one_ - 2 * eps_eps(prec, J)) *
                   (one_ - 3 * eps_eps(prec, J)) *
                   (eps_const(prec, J, 2) - 3 * eps_eps(prec, J)))
    main = eps_const(prec, J, A) * dSdt +
           (eps_const(prec, J, B0) + eps_shift(eps_const(prec, J, B1), 1)) * res2.S +
           3 * ((one_ - eps_eps(prec, J)) * (one_ - eps_eps(prec, J))) * tad * T1sq
    return pref * main
end
