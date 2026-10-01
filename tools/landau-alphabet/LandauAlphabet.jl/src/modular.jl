# The elliptic (Γ₁(6) sunrise/kite) kernel instance, per arXiv:1704.08895 and
# the G3 stratified ansatz strategy.
#
# Conventions (μ = m = 1, t in units of m²):
#   q₂ = exp(iπτ),  t(q₂) = −9 η(τ)⁴η(3τ/2)⁴η(6τ)⁴ / (η(τ/2)⁴η(2τ)⁴η(3τ)⁴)
#   ψ̂₁ := ψ₁/π = (2/√3)·u(q₂),  u rational with u(0)=1,
#   ψ̂₁² = 12 q t′(q) / (t(t−1)(t−9))      [from dt = ψ₁²/(iπW) dq/q]
#   q t′ = ψ̂₁² t(t−1)(t−9)/12.
#
# Curve-side kernel representations (all EXACT rational q₂-series):
#   g_{2,c} = q t′/(t−c),  c ∈ {0,1,9}          (dlogs of the Landau letters)
#   f₂ = g_{2,1} + g_{2,9} − ½ g_{2,0}
#   f₃ = −ψ̂₁³ t(t−1)(t−9)/24 = 3√3·f̂₃,  f̂₃ = −u³ t(t−1)(t−9)/81
#   f₁ = (t+3)ψ̂₁/(2√6),  f₄ = f₁⁴ = (t+3)⁴ψ̂₁⁴/576 = (t+3)⁴u⁴·(4/9)/576
#   g_{3,1} = (q t′)·ψ̂₁·t/(t−1) → rational ĝ₃₁ = u·qt′·t/(t−1)  [×(2/√3)]
#
# The maximal-cut/Γ₁(6) nome is q_C = −q₂ (maximal-cut convention); level-6 statements and
# the cusp form 6.4.a.a = η(τ_C)²η(2τ_C)²η(3τ_C)²η(6τ_C)² live in q_C.

export SunriseData, sunrise_data, qs_div_cancel, poly_in_t,
       solve_in_basis, level6_basis, sturm_bound_n, intersect_with_space,
       landau_weight2_candidates, landau_weight3_candidates, landau_weight4_candidates,
       ModularContext, modular_context, modular_value,
       cusp_constant, form_value, cusp_t_values, GAMMA16_CUSPS,
       qseries_value, q_of_t, psi1hat

# ---------- exact data pack ----------

struct SunriseData
    N::Int
    t::QS          # t(q₂)
    qtp::QS        # q t′(q₂)
    psi2::QS       # ψ̂₁² (constant term 4/3)
    u::QS          # u = (√3/2)ψ̂₁, u(0)=1, rational
    g20::QS; g21::QS; g29::QS; f2::QS
    f3hat::QS      # f₃/(3√3)
    f4::QS         # f₄ (rational)
    g31hat::QS     # g_{3,1}·(√3/2)  = u·qt′·t/(t−1)
    cusp4_qC::QS   # 6.4.a.a in q_C
end

"divide a/b allowing common leading zeros (exact cancellation)"
function qs_div_cancel(a::QS, b::QS)
    ka = findfirst(!iszero, a)
    kb = findfirst(!iszero, b)
    kb === nothing && error("division by zero series")
    ka === nothing && return qzero(length(a))
    k = kb - 1
    k == 0 && return qs_div(a, b)
    ka - 1 >= k || error("non-cancelling pole in series division")
    # truncate honestly: top k coefficients are unknown after the shift
    qs_div(a[k+1:end], b[k+1:end])
end

"evaluate a polynomial (coeff vector, constant first) on a series"
function poly_in_t(coeffs::Vector{QQFieldElem}, t::QS)
    N = length(t)
    acc = qzero(N)
    for c in reverse(coeffs)
        acc = qs_mul(acc, t)
        acc[1] += c
    end
    acc
end

function sunrise_data(N::Int)
    # t(q₂)
    lead, ser = eta_quotient_q2([(2, 4), (3, 4), (12, 4), (1, -4), (4, -4), (6, -4)], N)
    lead == 1 || error("t eta-quotient leading power ≠ 1")
    t = qs_scal(-9, qs_shift_power(ser, 1))
    qtp = qs_qdq(t)
    tm1 = copy(t); tm1[1] -= 1
    tm9 = copy(t); tm9[1] -= 9
    den = qs_mul(t, qs_mul(tm1, tm9))             # t(t−1)(t−9), leading q¹
    psi2 = qs_scal(12, qs_div_cancel(qtp, den))   # ψ̂₁²
    u = qs_sqrt_unit(qs_scal(QQ(3, 4), psi2))     # u² = (3/4)ψ̂₁²
    g20 = qs_div_cancel(qtp, t)
    g21 = qs_div(qtp, tm1)
    g29 = qs_div(qtp, tm9)
    f2 = qs_add(qs_add(g21, g29), qs_scal(QQ(-1, 2), g20))
    u3 = qs_mul(u, qs_mul(u, u))
    f3hat = qs_scal(QQ(-1, 81), qs_mul(u3, den))
    tp3 = copy(t); tp3[1] += 3
    tp34 = qs_mul(qs_mul(tp3, tp3), qs_mul(tp3, tp3))
    u4 = qs_mul(u3, u)
    f4 = qs_scal(QQ(16, 9 * 576), qs_mul(tp34, u4))  # ψ̂₁⁴ = (16/9)u⁴ → f₄ = (t+3)⁴u⁴/324
    g31hat = qs_mul(u, qs_div(qs_mul(qtp, t), tm1))
    lead4, ser4 = eta_quotient_q2([(1, 2), (2, 2), (3, 2), (6, 2)], N)
    lead4 == 1 || error("cusp form leading power ≠ 1")
    cusp4 = qs_shift_power(ser4, 1)
    # trim everything to the common valid length (divisions with cancellation
    # shorten series by the cancelled order)
    fields = (t, qtp, psi2, u, g20, g21, g29, f2, f3hat, f4, g31hat, cusp4)
    Nv = minimum(length.(fields))
    SunriseData(Nv, [f[1:Nv] for f in fields]...)
end

# ---------- exact linear algebra on q-series ----------

"""
    solve_in_basis(basis, target; nsolve, ncheck) -> Vector{QQFieldElem} or nothing

Solve target = Σ cᵢ basisᵢ using the first `nsolve` coefficients, then verify
exactly on `ncheck` coefficients. If the first `nsolve` rows underdetermine
the system (so the particular solution fails verification even though a valid
one exists), retries with twice the rows up to `ncheck`. Returns nothing only
when no solution verifies.
"""
function solve_in_basis(basis::Vector{QS}, target::QS;
                        nsolve::Int = 30, ncheck::Int = min(length(target), 200))
    m = length(basis)
    while true
        ns = min(nsolve, ncheck)
        A = zero_matrix(QQ, ns, m)
        b = zero_matrix(QQ, ns, 1)
        for i in 1:ns
            for j in 1:m
                A[i, j] = basis[j][i]
            end
            b[i, 1] = target[i]
        end
        fl, x = can_solve_with_solution(A, b, side = :right)
        if fl
            okall = true
            for i in 1:ncheck
                s = QQ(0)
                for j in 1:m
                    s += x[j, 1] * basis[j][i]
                end
                if s != target[i]
                    okall = false
                    break
                end
            end
            okall && return [x[j, 1] for j in 1:m]
        end
        ns >= ncheck && return nothing
        nsolve = 2 * ns
    end
end

"Sturm bound for M_k(Γ₁(6)) (index 24): k·24/12 = 2k, with margin."
sturm_bound_n(k::Int) = 2k + 12

"""
    level6_basis(sd, k) -> (names, Vector{QS})

Basis of M_k(Γ₁(6)) as exact q_C-series. Weight 2: B_{2,K}, K=2,3,6.
Weight 3 (character χ₋₃): E₃(Kτ;χ̄₁,χ̄₀), E₃(Kτ;χ̄₀,χ̄₁), K=1,2.
Weight 4: E₄(Kτ), K=1,2,3,6, plus the cusp form 6.4.a.a.
"""
function level6_basis(sd::SunriseData, k::Int)
    N = sd.N
    if k == 2
        return (["B22", "B23", "B26"], QS[B2K(2, N), B2K(3, N), B2K(6, N)])
    elseif k == 3
        return (["E3(1;chi3,1)", "E3(2;chi3,1)", "E3(1;1,chi3)", "E3(2;1,chi3)"],
                QS[eisenstein_E(3, -3, 1, 1, N), eisenstein_E(3, -3, 1, 2, N),
                   eisenstein_E(3, 1, -3, 1, N), eisenstein_E(3, 1, -3, 2, N)])
    elseif k == 4
        return (["E4(1)", "E4(2)", "E4(3)", "E4(6)", "cusp6.4.a.a"],
                QS[eisenstein_E(4, 1, 1, 1, N), eisenstein_E(4, 1, 1, 2, N),
                   eisenstein_E(4, 1, 1, 3, N), eisenstein_E(4, 1, 1, 6, N),
                   sd.cusp4_qC])
    end
    error("weight $k not implemented")
end

# ---------- Landau-data candidate spaces (the G3 stratified strategy) ----------

"""
    landau_weight2_candidates(sd) -> (names, q_C-series)

Weight-2 candidates from Landau data alone: dlog letters dt/(t−c) at the
finite Landau loci c ∈ {0,1,9}, transported to the modular side:
kernel = q t′/(t−c), mapped to the Γ₁(6) nome q_C = −q₂.
"""
function landau_weight2_candidates(sd::SunriseData)
    (["dlog(t)", "dlog(t-1)", "dlog(t-9)"],
     QS[qs_substitute_negq(sd.g20), qs_substitute_negq(sd.g21),
        qs_substitute_negq(sd.g29)])
end

"""
    landau_weight3_candidates(sd; maxdeg=3) -> (names, q_C-series)

Weight-3 candidates: ψ̂₁ × (q t′) × p(t)/(t(t−1)(t−9)) with deg p ≤ maxdeg,
rationally normalized as u·qt′·p(t)/(t(t−1)(t−9)). These are the one-forms a
single-elliptic-curve Landau analysis allows at this weight (simple poles at
the Landau loci, one period weight). Names record the t-monomial.
"""
function landau_weight3_candidates(sd::SunriseData; maxdeg::Int = 3)
    names = String[]
    out = QS[]
    tm1 = copy(sd.t); tm1[1] -= 1
    tm9 = copy(sd.t); tm9[1] -= 9
    den = qs_mul(sd.t, qs_mul(tm1, tm9))
    base = qs_mul(sd.u, sd.qtp)
    for d in 0:maxdeg
        coeffs = [QQ(i == d + 1 ? 1 : 0) for i in 1:d+1]
        num = poly_in_t(coeffs, sd.t)
        cand = qs_div_cancel(qs_mul(base, num), den)
        push!(names, "u*qt'*t^$d/(t(t-1)(t-9))")
        push!(out, qs_substitute_negq(cand))
    end
    (names, out)
end

"""
    landau_weight4_candidates(sd; maxdeg=4) -> (names, q_C-series)

Weight-4 candidates: ψ̂₁² × (q t′) × p(t)/(t(t−1)(t−9)), deg p ≤ maxdeg
(rationally normalized with u²).
"""
function landau_weight4_candidates(sd::SunriseData; maxdeg::Int = 4)
    names = String[]
    out = QS[]
    tm1 = copy(sd.t); tm1[1] -= 1
    tm9 = copy(sd.t); tm9[1] -= 9
    den = qs_mul(sd.t, qs_mul(tm1, tm9))
    u2 = qs_mul(sd.u, sd.u)
    base = qs_mul(u2, sd.qtp)
    for d in 0:maxdeg
        coeffs = [QQ(i == d + 1 ? 1 : 0) for i in 1:d+1]
        num = poly_in_t(coeffs, sd.t)
        cand = qs_div_cancel(qs_mul(base, num), den)
        push!(names, "u^2*qt'*t^$d/(t(t-1)(t-9))")
        push!(out, qs_substitute_negq(cand))
    end
    (names, out)
end

"""
    intersect_with_space(cands, basis; nrows, ncheck) -> Vector{Vector{QQFieldElem}}

All rational combinations c of `cands` with Σ cᵢ candᵢ ∈ span(basis), i.e. the
projection onto the candidate coefficients of ker[C | −B]. Solved on `nrows`
coefficients and exactly re-verified on `ncheck`.
"""
function intersect_with_space(cands::Vector{QS}, basis::Vector{QS};
                              nrows::Int = 40, ncheck::Int = 120)
    nc, nb = length(cands), length(basis)
    A = zero_matrix(QQ, nrows, nc + nb)
    for i in 1:nrows
        for j in 1:nc
            A[i, j] = cands[j][i]
        end
        for j in 1:nb
            A[i, nc+j] = -basis[j][i]
        end
    end
    nd, ns = nullspace(A)
    out = Vector{Vector{QQFieldElem}}()
    for k in 1:nd
        c = [ns[j, k] for j in 1:nc]
        all(iszero, c) && continue
        # verify: candidate combo solvable in basis to ncheck coefficients
        combo = qzero(min(ncheck, minimum(length.(cands))))
        for j in 1:nc
            iszero(c[j]) || (combo = qs_add(combo, qs_scal(c[j], cands[j][1:length(combo)])))
        end
        sol = solve_in_basis(basis, combo; nsolve = nrows, ncheck = length(combo))
        sol === nothing && error("intersection candidate failed exact verification")
        push!(out, c)
    end
    out
end

# ---------- numeric layer: forms, cusps, iterated integrals ----------

"QQFieldElem → Acb at given precision"
acbq(x::QQFieldElem, prec::Int) =
    Acb(Rational{BigInt}(BigInt(numerator(x)), BigInt(denominator(x))), prec = prec)

"evaluate an exact q-series at a complex ball q, with geometric tail estimate"
function qseries_value(series::QS, q::Acb; prec::Int = precision(q))
    s = Acb(0, prec = prec)
    for n in length(series):-1:1
        Arblib.mul!(s, s, q)
        s += acbq(series[n], prec)
    end
    aq = ub(abs(q))
    cmax = maximum(abs(Float64(series[n])) for n in max(1, length(series) - 10):length(series))
    tail = (cmax + 1) * aq^length(series) / max(1 - aq, Arb(1e-6)) * 64
    Arblib.add_error!(s, magof(tail))
    s
end

"""
    q_of_t(sd, tval; prec) -> Acb

Invert t(q₂) = tval by Newton iteration on the exact series (Euclidean t < 0
maps to small positive q₂).
"""
function q_of_t(sd::SunriseData, tval::Acb; prec::Int = precision(tval))
    qtp = sd.qtp
    # Newton on midpoints (point arithmetic — no ball inflation across
    # iterations), then a single residual-based error attachment at the end.
    mid(z) = Acb(Arb(Arblib.midref(real(z)), prec = prec),
                 Arb(Arblib.midref(imag(z)), prec = prec), prec = prec)
    q = mid(tval / Acb(-9, prec = prec))
    lastdq = Arb(1, prec = 64)
    for _ in 1:300
        f = qseries_value(sd.t, q; prec = prec) - tval
        df = qseries_value(qtp, q; prec = prec) / q
        dq = f / df
        q = mid(q - dq)
        nd = ub(abs(dq))
        # quadratic convergence stalls at the evaluation noise floor
        nd >= lastdq && Float64(lastdq) < 1e-10 && break
        lastdq = nd
    end
    res = qseries_value(sd.t, q; prec = prec) - tval
    df = qseries_value(qtp, q; prec = prec) / q
    Arblib.add_error!(q, magof(2 * ub(abs(res / df)) + Arb(BigFloat(2)^(-prec))))
    q
end

"ψ̂₁ = ψ₁/π = (2/√3)·u(q₂)"
psi1hat(sd::SunriseData, q::Acb; prec::Int = precision(q)) =
    2 / sqrt(Acb(3, prec = prec)) * qseries_value(sd.u, q; prec = prec)

"evaluate a q-series form numerically: f(τ) = Σ c_n q^n with q = e^{2πiwτ}"
function form_value(series::QS, tau::Acb, prec::Int; halfstep::Bool = false)
    # halfstep: series in q₂ = e^{iπτ}; else series in q_C = e^{2πiτ}
    q = exp((halfstep ? 1 : 2) * Acb(0, 1, prec = prec) * Arb(π, prec = prec) * tau)
    s = Acb(0, prec = prec)
    for n in length(series):-1:1
        Arblib.mul!(s, s, q)
        s += acbq(series[n], prec)
    end
    # geometric tail estimate (coefficients of holomorphic level-6 forms grow
    # polynomially; |q|<1 strictly in all our uses)
    aq = ub(abs(q))
    cmax = maximum(abs(Float64(series[n])) for n in max(1, length(series) - 10):length(series))
    tail = (cmax + 1) * aq^length(series) / max(1 - aq, Arb(1e-6)) * 64
    Arblib.add_error!(s, magof(tail))
    s
end

"Γ₁(6) cusps as (a, c, width)"
const GAMMA16_CUSPS = [(1, 0, 1), (0, 1, 6), (1, 2, 3), (1, 3, 2)]

"""
    cusp_constant(series, k, (a,c,width); y, M, prec, halfstep=false) -> Acb

Constant term of the weight-k form at the cusp a/c, computed by averaging
(cτ+d)^{−k} f(γτ) over M points on a horocycle of height y. The aliasing
error e^{−2πyM/width} is added to the error ball.
"""
function cusp_constant(series::QS, k::Int, cusp::Tuple{Int,Int,Int};
                       y::Float64 = 6.0, M::Int = 48, prec::Int = 512,
                       halfstep::Bool = false)
    a, c, h = cusp
    # γ = [a b; c d], ad − bc = 1
    g, d0, b0 = gcdx(a, -c)   # a*d0 − c*b0 = 1 when g == 1
    g == 1 || error("bad cusp")
    b, d = b0, d0
    s = Acb(0, prec = prec)
    for j in 0:M-1
        x = Arb(j * h, prec = prec) / M
        tau = Acb(x, Arb(y, prec = prec), prec = prec)
        gt = (a * tau + b) / (c * tau + d)
        s += form_value(series, gt, prec; halfstep = halfstep) / (c * tau + d)^k
    end
    s /= M
    err = exp(-2 * Arb(π, prec = 64) * y * M / h) * 1e6
    Arblib.add_error!(s, magof(err))
    s
end

"t-values at the Γ₁(6) cusps (identifies cusp ↔ Landau point); t series in q₂"
function cusp_t_values(sd::SunriseData; y::Float64 = 8.0, prec::Int = 256)
    vals = Acb[]
    tC = qs_substitute_negq(sd.t)   # t as q_C-series
    for (a, c, h) in GAMMA16_CUSPS
        g, d0, b0 = gcdx(a, -c)
        b, d = b0, d0
        tau = Acb(Arb(h, prec = prec) / 7, Arb(y, prec = prec), prec = prec)
        gt = (a * tau + b) / (c * tau + d)
        push!(vals, form_value(tC, gt, prec))
    end
    vals
end

# ---------- numeric modular iterated integrals ----------

struct ModularContext
    names::Vector{String}
    weights::Vector{Int}        # ε-grading degree of each letter
    ev::ItIntEvaluator
end

"""
    modular_context(sd, letters; prec, nterms) -> ModularContext

`letters` = vector of (name, qseries-in-q₂, scale::Acb, εdeg). Builds
SeriesKernels for forms f(q₂)·scale, as one-forms f dq₂/q₂.
"""
function modular_context(letters::Vector; prec::Int, nterms::Int)
    names = String[]
    weights = Int[]
    kernels = KernelForm[]
    # never silently zero-pad a kernel: shrink the evaluator to the shortest
    # supplied series so the tail heuristic operates at the true truncation
    minlen = minimum(length(ser) for (_, ser, _, _) in letters)
    if minlen < nterms
        @warn "modular_context: shortest letter has $minlen coefficients < nterms=$nterms; using nterms=$minlen"
        nterms = minlen
    end
    for (name, ser, scale, deg) in letters
        n = min(nterms, length(ser))
        pole = acbq(ser[1], prec) * scale
        taylor = [acbq(ser[i], prec) * scale for i in 2:n]
        push!(kernels, SeriesKernel(pole, taylor))
        push!(names, name)
        push!(weights, deg)
    end
    ModularContext(names, weights, ItIntEvaluator(kernels; nterms = nterms, prec = prec))
end

"I(word; q₂) — first letter outermost, per 1704.08895"
modular_value(ctx::ModularContext, w::Word, q2::Acb) = eval_word(ctx.ev, w, q2)
