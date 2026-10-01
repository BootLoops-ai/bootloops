# Fit layer: lattice-reduction analytic regression and integer-relation
# (PSLQ-style) fits, with acceptance gates: coefficient-height bound,
# residual at full precision, and cross-validation at held-out kinematic
# points.

export lattice_fit, integer_relation, FitResult, fit_residual

struct FitResult
    ok::Bool
    coeffs::Vector{QQFieldElem}   # f(x) = Σ cᵢ Bᵢ(x)
    denom::ZZRingElem             # common scale (u_{p+1})
    raw::Vector{ZZRingElem}
    height::ZZRingElem
    note::String
end

function arb_to_bigfloat(x::Arb, prec::Int)
    setprecision(BigFloat, prec + 64) do
        BigFloat(string(Arblib.midref(x)))   # decimal round-trip at high precision
    end
end

"round(10^s · x) for Arb x — ALL arithmetic at elevated BigFloat precision"
function scaled_int(x::Arb, s::Int, prec::Int)
    isfinite(x) || error("scaled_int: non-finite input ball")
    setprecision(BigFloat, prec + 64) do
        bf = BigFloat(string(Arblib.midref(x)))
        isfinite(bf) || error("scaled_int: midpoint not representable")
        ZZ(BigInt(round(bf * BigFloat(10)^s)))
    end
end

"""
    lattice_fit(fvals, Bvals; digits, height_bound=ZZ(10)^12) -> FitResult

Analytic regression of f(x⃗) = Σ cᵢ Bᵢ(x⃗) from values at p points
(`fvals::Vector{Arb}` length p, `Bvals::Matrix{Arb}` n×p) by LLL reduction:
the lattice adjoins an identity block for the integer coefficient vector to
the digit-scaled sampled values, so a short vector is a small-height candidate
relation. `digits` = trusted significant digits.

IMPORTANT: `FitResult.ok` is the height gate ONLY. A candidate fit must
additionally pass `fit_residual` at full precision on the fit points AND on
held-out cross-validation points before being accepted —
the height gate alone can pass noise when
n ≳ digits·p / log₁₀(height_bound).
"""
function lattice_fit(fvals::Vector{Arb}, Bvals::Matrix{Arb};
                     digits::Int, prec::Int = 4096,
                     height_bound::ZZRingElem = ZZ(10)^12)
    n, p = size(Bvals)
    length(fvals) == p || error("dimension mismatch")
    # scale: s = digits − Δmax (magnitudes from ball UPPER bounds, so wide or
    # zero-centered balls cannot produce NaN/underestimates)
    allv = vcat(fvals, vec(Bvals))
    dmax = maximum([Float64(log10(ub(abs(v)) + Arb(1e-300))) for v in allv])
    s = digits - ceil(Int, dmax)
    M = zero_matrix(ZZ, n + 1, p + n + 1)
    for j in 1:p
        M[1, j] = scaled_int(fvals[j], s, prec)
        for i in 1:n
            M[i+1, j] = scaled_int(Bvals[i, j], s, prec)
        end
    end
    for i in 1:n+1
        M[i, p+i] = 1
    end
    R = lll(M)
    # first row with nonzero f-coefficient
    for r in 1:n+1
        u = [R[r, p+i] for i in 1:n+1]
        if !iszero(u[1])
            coeffs = [QQ(-u[1+i], u[1]) for i in 1:n]
            h = maximum(vcat([abs(numerator(c)) for c in coeffs],
                             [abs(denominator(c)) for c in coeffs], [ZZ(1)]))
            ok = h <= height_bound
            return FitResult(ok, coeffs, ZZRingElem(u[1]), [ZZRingElem(x) for x in u], h,
                             ok ? "ok" : "height bound exceeded ($h)")
        end
    end
    FitResult(false, QQFieldElem[], ZZ(0), ZZRingElem[], ZZ(0), "no relation involving f found")
end

"""
    fit_residual(fr, fvals, Bvals) -> Arb

max_j |f(x_j) − Σ cᵢBᵢ(x_j)| evaluated in ball arithmetic — the full-precision
residual acceptance gate (use held-out points for cross-validation).
"""
function fit_residual(fr::FitResult, fvals::Vector{Arb}, Bvals::Matrix{Arb})
    n, p = size(Bvals)
    worst = Arb(0)
    for j in 1:p
        r = fvals[j]
        for i in 1:n
            r -= Arb(Rational{BigInt}(fr.coeffs[i]), prec = precision(fvals[j])) * Bvals[i, j]
        end
        a = abs(r)
        worst = max(worst, a)
    end
    worst
end

"""
    integer_relation(vals::Vector{Arb}; digits, height_bound) -> (ok, rel::Vector{ZZRingElem})

PSLQ-equivalent single-point integer relation Σ mᵢ valsᵢ ≈ 0 via LLL.

Acceptance combines three gates (any LLL row trivially satisfies the residual
bound by construction, so the residual alone is NOT discriminating):
height bound; residual consistent with zero in ball arithmetic; and the
exclusion bound n·log₁₀(height) ≤ 0.6·digits — spurious relations among
unrelated d-digit reals have height ~10^{d/n}, so a trustworthy relation must
sit far below that Minkowski floor.
"""
function integer_relation(vals::Vector{Arb}; digits::Int, prec::Int = 4096,
                          height_bound::ZZRingElem = ZZ(10)^10)
    n = length(vals)
    dmax = maximum([Float64(log10(ub(abs(v)) + Arb(1e-300))) for v in vals])
    s = digits - ceil(Int, dmax)
    M = zero_matrix(ZZ, n, n + 1)
    for i in 1:n
        M[i, 1] = scaled_int(vals[i], s, prec)
        M[i, 1+i] = 1
    end
    R = lll(M)
    rel = [R[1, 1+i] for i in 1:n]
    h = maximum(abs.(rel))
    # residual gate in ball arithmetic: |Σ mᵢvᵢ| must be consistent with zero at
    # ~the working digits (upper bound of the ball, not the midpoint)
    r = Arb(0, prec = precision(vals[1]))
    for i in 1:n
        r += Arb(BigInt(rel[i]), prec = precision(vals[1])) * vals[i]
    end
    resmag = BigFloat(Arblib.ubound(Arf, abs(r)))
    thresh = BigFloat(10)^(-(7 * digits) ÷ 10) * max(BigFloat(1), BigFloat(h))
    # exclusion bound (see docstring): noise relations have height ~10^{digits/n}
    confident = n * log10(Float64(BigInt(max(h, ZZ(2))))) <= 0.6 * digits
    ok = h <= height_bound && resmag < thresh && confident
    (ok, [ZZRingElem(x) for x in rel], r)
end
