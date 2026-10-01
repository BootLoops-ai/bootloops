# Taylor-coefficient extraction via Cauchy circles in ball arithmetic.
# Used for boundary-constant generation (e.g. the ε-expansion of the exact
# t=0 sunrise closed form of 1704.08895 eq. (boundary)).

export cauchy_coeffs, sunrise_boundary_constants

"""
    cauchy_coeffs(f, jmax; r=1//16, M=128, prec=1536) -> Vector{Acb}

Taylor coefficients c₀…c_jmax of an analytic f at 0 via the discrete Cauchy
integral on |ε| = r with M points. Aliasing error is estimated from the
decay of the computed top coefficients and added to the error balls.
"""
function cauchy_coeffs(f::Function, jmax::Int; r::Rational{Int} = 1 // 16,
                       M::Int = 128, prec::Int = 1536)
    rr = Arb(r, prec = prec)
    vals = Vector{Acb}(undef, M)
    for j in 0:M-1
        theta = 2 * Arb(π, prec = prec) * j / M
        eps = Acb(rr * cos(theta), rr * sin(theta), prec = prec)
        vals[j+1] = f(eps)
    end
    out = Vector{Acb}(undef, jmax + 1)
    vmax = maximum(ub(abs(v)) for v in vals)
    for k in 0:jmax
        s = Acb(0, prec = prec)
        for j in 0:M-1
            theta = -2 * Arb(π, prec = prec) * j * k / M
            s += vals[j+1] * Acb(cos(theta), sin(theta), prec = prec)
        end
        s /= (M * rr^k)
        # aliasing: |ĉ_k − c_k| ≤ Σ_{l≥1} |c_{k+lM}| r^{lM} ≲ vmax (r/ρ)^M r^{−k};
        # we bound with ρ ≥ 4r (radius of analyticity margin, checked by caller)
        err = vmax * (Arb(1, prec = 64) / 4)^M / rr^k * 16
        Arblib.add_error!(s, magof(err))
        out[k+1] = s
    end
    out
end

"""
    sunrise_boundary_constants(jmax; prec=1536) -> Vector{Acb}

S_111^{(j)}(2,0) for j = 0…jmax (μ = m = 1), from the exact closed form
(1704.08895): Σ εʲ S^{(j)} = e^{2γε} Γ(1+2ε) (√3)^{−1−2ε}
[ 3Γ(1+ε)²/(2ε²Γ(1+2ε)) · h(ε) − π/ε ],
h = (1/i)[(−r₃)^{−ε} ₂F₁(−2ε,−ε;1−ε;r₃) − (−r₃⁻¹)^{−ε} ₂F₁(−2ε,−ε;1−ε;r₃⁻¹)].
"""
function sunrise_boundary_constants(jmax::Int; prec::Int = 1536, M::Int = 128,
                                    r::Rational{Int} = 1 // 16)
    i_ = Acb(0, 1, prec = prec)
    r3 = exp(2 * Arb(π, prec = prec) / 3 * i_)
    r3i = inv(r3)
    sqrt3 = sqrt(Arb(3, prec = prec))
    gam = Arblib.const_euler!(Arb(prec = prec))
    function f(eps::Acb)
        F1 = Arblib.hypgeom_2f1!(Acb(prec = prec), -2eps, -eps, 1 - eps, r3,
                                 flags = 0, prec = prec)
        F2 = Arblib.hypgeom_2f1!(Acb(prec = prec), -2eps, -eps, 1 - eps, r3i,
                                 flags = 0, prec = prec)
        h = (exp(-eps * log(-r3)) * F1 - exp(-eps * log(-r3i)) * F2) / i_
        G1 = Arblib.gamma!(Acb(prec = prec), 1 + eps)
        G2 = Arblib.gamma!(Acb(prec = prec), 1 + 2eps)
        pref = exp(2 * gam * eps) * G2 * exp((-1 - 2eps) * log(Acb(sqrt3, prec = prec)))
        pref * (3 * G1^2 / (2 * eps^2 * G2) * h - Acb(Arb(π, prec = prec), prec = prec) / eps)
    end
    cauchy_coeffs(f, jmax; r = r, M = M, prec = prec)
end
