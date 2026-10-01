# Certified Taylor transport ("Frobenius continuation") for linear ODEs
#   Σ_{i=0}^{r} p_i(x) y^{(i)}(x) = 0
# whose coefficients p_i are polynomials in x over the ring of truncated ε-power series
# with Acb ball coefficients (Jeps = 1 gives the plain Acb case).
#
# Rigor of the truncation: at an ordinary point x0 the Taylor coefficients c_n of any
# solution satisfy the (r+d)-term recurrence implied by the operator. We prove a
# geometric envelope ‖c_n‖ ≤ C·ρ⁻ⁿ for all n ≥ W_lo by interval verification:
#   * normalize the recurrence by p_{r}(x0)⁻¹ (ε-series inverse, ℓ₁ norms ‖·‖),
#   * bracket(m) = Σ_{(i,j)≠(r,0)} ‖p̃_ij‖ ρ^{r-i+j} γ_ij(m) with
#     γ_ij(m) = (m-j+i)!/(m-j)! · m!/(m+r)! ≤ (m+i)^i/(m+1)^r (i<r, decreasing in m),
#     γ_rj(m) ≤ 1 (j ≥ 1);
#   * if bracket(m_check) ≤ 1 for m_check = W_hi - r + 1 and C = max_{W_lo ≤ n ≤ W_hi}
#     ‖c_n‖ρⁿ with window width ≥ r + d, induction gives the envelope for all n > W_hi.
# The ℓ₁ norm on truncated ε-series is submultiplicative, and an envelope for the ℓ₁
# norm bounds every ε-coefficient, so the bound is valid order by order in ε.

# ---- ε-polynomial helpers (lead-0 truncated series, Vector{Acb} of length Jeps) ----

epoly_const(prec, J, v) = (c = [Acb(0; prec) for _ in 1:J]; c[1] = Acb(v; prec); c)
epoly_zero(prec, J) = [Acb(0; prec) for _ in 1:J]

function epoly_mul(a::Vector{Acb}, b::Vector{Acb})
    J = length(a); prec = Arblib.precision(a[1])
    c = [Acb(0; prec) for _ in 1:J]
    for i in 1:J, j in 1:(J - i + 1)
        c[i+j-1] += a[i] * b[j]
    end
    c
end
epoly_add(a, b) = [a[k] + b[k] for k in 1:length(a)]
epoly_sub(a, b) = [a[k] - b[k] for k in 1:length(a)]
epoly_scal(s::Acb, a) = [s * x for x in a]

function epoly_inv(a::Vector{Acb})
    J = length(a); prec = Arblib.precision(a[1])
    Arblib.contains_zero(a[1]) && error("epoly_inv: constant term contains zero")
    c = [Acb(0; prec) for _ in 1:J]
    c[1] = inv(a[1])
    for k in 2:J
        s = Acb(0; prec)
        for j in 2:k
            s += a[j] * c[k-j+1]
        end
        c[k] = -c[1] * s
    end
    c
end

"""ℓ₁ Mag norm of an ε-polynomial."""
function epoly_norm(a::Vector{Acb})
    s = Mag(0); t = Mag()
    for x in a
        Arblib.get!(t, x)
        Arblib.add!(s, s, t)
    end
    s
end

# ---- operator ----------------------------------------------------------------

"""Operator Σ p_i(x) y^{(i)}; coeffs[i+1] = vector of ε-polys, the x-coefficients of p_i."""
struct ODEOperator
    r::Int
    Jeps::Int
    prec::Int
    coeffs::Vector{Vector{Vector{Acb}}}
    sing::Vector{Acb}     # singular points in x (roots of leading coefficient + ∞ handled by caller)
end

"""Shift to local coordinates u = x - x0: p_i(x0 + u) via binomial re-expansion."""
function localize(op::ODEOperator, x0::Acb)
    out = Vector{Vector{Vector{Acb}}}(undef, op.r + 1)
    for i in 0:op.r
        p = op.coeffs[i+1]
        d = length(p) - 1
        q = [epoly_zero(op.prec, op.Jeps) for _ in 0:d]
        # p(x0+u) = Σ_j p_j (x0+u)^j ; expand each (x0+u)^j
        for j in 0:d
            # binomial: (x0+u)^j = Σ_k C(j,k) x0^{j-k} u^k
            for k in 0:j
                fac = Acb(binomial(big(j), big(k)); prec = op.prec) * x0^(j - k)
                q[k+1] = epoly_add(q[k+1], epoly_scal(fac, p[j+1]))
            end
        end
        out[i+1] = q
    end
    ODEOperator(op.r, op.Jeps, op.prec, out, op.sing)
end

"""Generate Taylor coefficients c_0..c_M at ordinary point (local operator `lop`),
from initial c_0..c_{r-1} (each an ε-poly)."""
function taylor_coeffs(lop::ODEOperator, c0::Vector{Vector{Acb}}, M::Int)
    r, prec, J = lop.r, lop.prec, lop.Jeps
    c = Vector{Vector{Acb}}(undef, M + 1)
    for i in 1:r
        c[i] = c0[i]
    end
    pr0inv = epoly_inv(lop.coeffs[r+1][1])
    for m in 0:(M - r)
        # solve coefficient of u^m: Σ_{i,j} p_ij ((m-j+i)!/(m-j)!) c_{m-j+i} = 0
        acc = epoly_zero(prec, J)
        for i in 0:r
            p = lop.coeffs[i+1]
            for j in 0:length(p)-1
                (i == r && j == 0) && continue
                m < j && continue
                idx = m - j + i
                idx > m + r - 1 && continue   # cannot happen except (r,0)
                ff = prod(big(m - j + 1):big(m - j + i); init = big(1))
                acc = epoly_add(acc, epoly_scal(Acb(ff; prec = prec), epoly_mul(p[j+1], c[idx+1])))
            end
        end
        ffr = prod(big(m + 1):big(m + r); init = big(1))
        cnew = epoly_scal(Acb(-1; prec) / Acb(ffr; prec = prec), epoly_mul(pr0inv, acc))
        c[m+r+1] = cnew
    end
    return c
end

"""Certified geometric envelope: returns (C, ok) with ‖c_n‖ ≤ C ρ⁻ⁿ for all n ≥ W_lo."""
function envelope(lop::ODEOperator, c::Vector{Vector{Acb}}, ρ::Mag)
    r, prec = lop.r, lop.prec
    d = maximum(length(p) - 1 for p in lop.coeffs)
    M = length(c) - 1
    W_hi = M
    W_lo = W_hi - (r + d)
    W_lo >= r || return (Mag(0), false)
    m_check = W_hi - r + 1
    pr0inv = epoly_inv(lop.coeffs[r+1][1])
    ρa = Arb(prec = prec); Arblib.set!(ρa, ρ)
    bracket = Arb(0; prec)
    for i in 0:r
        p = lop.coeffs[i+1]
        for j in 0:length(p)-1
            (i == r && j == 0) && continue
            nrm = Arb(prec = prec)
            Arblib.set!(nrm, epoly_norm(epoly_mul(pr0inv, p[j+1])))
            γ = if i == r
                Arb(1; prec)
            else
                Arb(big(m_check + i)^i; prec) / Arb(big(m_check + 1)^r; prec)
            end
            bracket += nrm * ρa^(r - i + j) * γ
        end
    end
    ok = bracket < 1
    C = Mag(0); t = Mag(); ρp = Mag()
    for n in W_lo:W_hi
        Arblib.pow!(ρp, ρ, UInt(n))
        Arblib.mul!(t, epoly_norm(c[n+1]), ρp)
        Arblib.cmp(t, C) > 0 && Arblib.set!(C, t)
    end
    return (C, ok)
end

"""One certified step: from (y, y', ..., y^{(r-1)})(x0) to the same at x0+h.
Auto-selects ρ from the certified distance to singular points; errors if |h| ≥ ρ/2
(caller subdivides). M = number of Taylor terms."""
function transport_step(op::ODEOperator, x0::Acb, h::Acb, Y::Vector{Vector{Acb}}; M::Int = 64)
    r, prec, J = op.r, op.prec, op.Jeps
    lop = localize(op, x0)
    # distance to singularities
    dist = Arb(prec = prec); first = true
    for s in op.sing
        ds = abs(s - x0)
        if first; dist = ds; first = false
        else
            dist = min(dist, ds)
        end
    end
    hm = Mag(); Arblib.get!(hm, h)
    # initial normalized coefficients c_i = y^{(i)}/i!
    c0 = [epoly_scal(Acb(1; prec) / Acb(factorial(big(i - 1)); prec), Y[i]) for i in 1:r]
    local c, C, ρ
    distf = Float64(Arblib.lbound(Arb, dist))
    hf = Float64(hm)
    hf < distf || error("transport_step: |h| not certifiably below singularity distance")
    found = false
    for θ in (0.75, 0.5, 0.35, 0.25, 0.18)
        ρf = distf * θ
        ρf > hf * 1.2 || continue
        ρ = Mag(ρf)
        # Taylor order: tail decays like (|h|/ρ)^M (envelope) but in practice like
        # (|h|/dist)^M; size M off the certified ratio with generous headroom.
        β = hf / ρf
        Mneed = ceil(Int, 1.3 * (prec + 60) * log(2.0) / log(1 / β)) + 2 * (op.r + 4)
        Muse = max(M, Mneed)
        c = taylor_coeffs(lop, c0, Muse)
        C, ok = envelope(lop, c, ρ)
        if ok
            found = true
            break
        end
    end
    found || error("transport_step: envelope verification failed at all trial radii")
    M = length(c) - 1
    # evaluate y^{(i)}(x0+h) = Σ_n n!/(n-i)! c_n h^{n-i} with certified tail
    out = Vector{Vector{Acb}}(undef, r)
    ρ_arb = Arb(prec = prec); Arblib.set!(ρ_arb, ρ)
    hma = Arb(prec = prec); Arblib.set!(hma, hm)
    for i in 0:r-1
        s = [Acb(0; prec) for _ in 1:J]
        for n in i:M
            ff = prod(big(n - i + 1):big(n); init = big(1))
            s = epoly_add(s, epoly_scal(Acb(ff; prec = prec) * h^(n - i), c[n+1]))
        end
        # tail: Σ_{n>M} n^i ‖c_n‖ |h|^{n-i} ≤ C/|h|^i Σ_{n>M} n^i (|h|/ρ)ⁿ
        θ = hma / ρ_arb
        Ca = Arb(prec = prec); Arblib.set!(Ca, C)
        # Σ_{n>M} n^i θⁿ ≤ (M+1)^i θ^{M+1}/(1 - θ((M+2)/(M+1))^i)  (require <1)
        den = 1 - θ * (Arb(M + 2; prec) / (M + 1))^i
        Arblib.is_positive(den) || error("transport_step: tail ratio ≥ 1")
        tail = Ca * Arb(big(M + 1)^i; prec) * θ^(M + 1) / den / hma^i
        tm = Mag(); Arblib.get!(tm, Arblib.ubound(tail))
        for k in 1:J
            Arblib.add_error!(s[k], tm)
        end
        out[i+1] = s
    end
    return out
end

"""Transport initial data Y along a piecewise-straight path (list of Acb anchor points).
Subdivides each leg adaptively so |h| < (certified distance to singularities)/2·θ."""
function transport(op::ODEOperator, path::Vector{Acb}, Y::Vector{Vector{Acb}}; M::Int = 64,
                   maxsteps::Int = 100000)
    prec = op.prec
    x = path[1]
    nsteps = 0
    for target in path[2:end]
        while !Arblib.overlaps(x, target) || !(x == target)
            (nsteps += 1) <= maxsteps || error("transport: step limit exceeded (wide input balls?)")
            # certified distance from x to singularities
            dist = minimum([abs(s - x) for s in op.sing])
            maxstep = Arblib.lbound(Arb, dist) / 8
            δ = target - x
            δm = abs(δ)
            if δm < maxstep
                Y = transport_step(op, x, δ, Y; M = M)
                x = target
                break
            else
                frac = maxstep / δm
                h = δ * Acb(Arblib.lbound(Arb, frac); prec = prec)
                Y = transport_step(op, x, h, Y; M = M)
                x = x + h
            end
        end
    end
    return Y, x
end
