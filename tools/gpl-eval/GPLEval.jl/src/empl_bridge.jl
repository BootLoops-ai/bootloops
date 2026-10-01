# eMPL bridge — thin layer over Eichler.jl.
#
# For elliptic letters (Kronecker–Eisenstein g⁽ⁿ⁾(z,τ) kernels, or any kernel
# represented as a `QSeries` in the nome q), iterated integration is delegated
# to `Eichler.iterated_integral` (certified LogQSeries with rigorous tail).
# The kinematic point t is mapped to (q, log q, τ) via the `EllipticSector`'s
# τ-map (cusp-connected branch, t+i0 prescription).

"""
    empl(kernels::Vector{QSeries}, t, sec::EllipticSector) -> Acb

Certified value of the iterated Eichler integral I(k₁,…,kₙ; q(t)) at kinematic
point `t`, with the tangential-base-point regularisation. First kernel =
outermost integration (same convention as `gpl`).
"""
function empl(kernels::Vector{QSeries}, t, sec::EllipticSector)
    τ, q, L = tau_map(sec, t)
    V = iterated_integral(kernels)
    return evaluate(V, q, L)
end

"""
    empl(kernels::Vector{QSeries}, q::Acb, logq::Acb) -> Acb

Direct (q, log q) interface (caller supplies the nome and its log = 2πiτ).
"""
empl(kernels::Vector{QSeries}, q::Acb, logq::Acb) =
    evaluate(iterated_integral(kernels), q, logq)
