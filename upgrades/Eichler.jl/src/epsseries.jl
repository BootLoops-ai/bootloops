# Truncated Laurent series in ε with Acb ball coefficients. Purely formal: a value
# X = ε^{lead} Σ_{k=0}^{J-1} c_k ε^k, carried to a fixed number J of stored orders.
# All operations are exact order-by-order (no analytic tails needed).

struct EpsSeries
    lead::Int            # leading exponent
    c::Vector{Acb}       # length J, coefficients of ε^{lead}, ε^{lead+1}, ...
end

nord(x::EpsSeries) = length(x.c)
prec_of(x::EpsSeries) = Arblib.precision(x.c[1])

eps_zero(prec, J, lead = 0) = EpsSeries(lead, [Acb(0; prec) for _ in 1:J])
function eps_const(prec, J, v)
    c = [Acb(0; prec) for _ in 1:J]; c[1] = Acb(v; prec)
    EpsSeries(0, c)
end
function eps_eps(prec, J)  # the variable ε
    c = [Acb(0; prec) for _ in 1:J]; c[1] = Acb(1; prec)
    EpsSeries(1, c)
end

"""Trim/align two series to common leading exponent and length (min overlap)."""
function _align(x::EpsSeries, y::EpsSeries)
    prec = prec_of(x)
    lead = min(x.lead, y.lead)
    # number of valid orders: x valid for exponents [x.lead, x.lead+J-1]
    hi = min(x.lead + nord(x), y.lead + nord(y))
    J = hi - lead
    cx = [Acb(0; prec) for _ in 1:J]; cy = [Acb(0; prec) for _ in 1:J]
    for k in 1:nord(x)
        idx = x.lead + k - 1 - lead + 1
        1 <= idx <= J && (cx[idx] = x.c[k])
    end
    for k in 1:nord(y)
        idx = y.lead + k - 1 - lead + 1
        1 <= idx <= J && (cy[idx] = y.c[k])
    end
    return lead, cx, cy
end

function Base.:+(x::EpsSeries, y::EpsSeries)
    lead, cx, cy = _align(x, y)
    EpsSeries(lead, [cx[k] + cy[k] for k in 1:length(cx)])
end
Base.:-(x::EpsSeries) = EpsSeries(x.lead, [-ck for ck in x.c])
Base.:-(x::EpsSeries, y::EpsSeries) = x + (-y)

function Base.:*(x::EpsSeries, y::EpsSeries)
    prec = prec_of(x)
    J = min(nord(x), nord(y))
    c = [Acb(0; prec) for _ in 1:J]
    for i in 1:J, j in 1:(J - i + 1)
        c[i+j-1] += x.c[i] * y.c[j]
    end
    EpsSeries(x.lead + y.lead, c)
end
Base.:*(a::Acb, x::EpsSeries) = EpsSeries(x.lead, [a * ck for ck in x.c])
Base.:*(a::Union{Integer,Rational}, x::EpsSeries) = Acb(a; prec = prec_of(x)) * x

"""Multiplicative inverse (requires c[1] invertible); result has lead = -x.lead."""
function eps_inv(x::EpsSeries)
    prec = prec_of(x); J = nord(x)
    Arblib.contains_zero(x.c[1]) && error("eps_inv: leading coefficient contains zero")
    inv0 = inv(x.c[1])
    c = [Acb(0; prec) for _ in 1:J]
    c[1] = inv0
    for k in 2:J
        s = Acb(0; prec)
        for j in 2:k
            s += x.c[j] * c[k-j+1]
        end
        c[k] = -inv0 * s
    end
    EpsSeries(-x.lead, c)
end
Base.:/(x::EpsSeries, y::EpsSeries) = x * eps_inv(y)

"""exp of a series with lead ≥ 1 (no constant term). The input is known to absolute
order lead+nord-1; the output (lead 0) is exact to the same absolute order."""
function eps_exp(x::EpsSeries)
    prec = prec_of(x)
    x.lead >= 1 || error("eps_exp expects positive leading exponent")
    Jout = x.lead + nord(x)            # orders ε⁰ .. ε^{Jout-1}
    f = [Acb(0; prec) for _ in 1:Jout-1]   # f[i] = coeff of ε^i
    for k in 1:nord(x)
        i = x.lead + k - 1
        i <= Jout - 1 && (f[i] = x.c[k])
    end
    e = [Acb(0; prec) for _ in 1:Jout]; e[1] = Acb(1; prec)
    for k in 1:Jout-1                  # e_k = (1/k) Σ_{j=1}^{k} j f_j e_{k-j}
        s = Acb(0; prec)
        for j in 1:k
            s += j * f[j] * e[k-j+1]
        end
        e[k+1] = s / k
    end
    EpsSeries(0, e)
end

"""Γ(1 + a·ε) as an EpsSeries to J orders, via Arb's certified gamma_series."""
function eps_gamma1p(prec::Int, J::Int, a::Integer)
    P = AcbPoly(prec = prec)
    Arblib.set_coeff!(P, 0, Acb(1; prec))
    Arblib.set_coeff!(P, 1, Acb(a; prec))
    Q = AcbPoly(prec = prec)
    Arblib.gamma_series!(Q, P, J; prec = prec)
    EpsSeries(0, [n <= Arblib.degree(Q) ? Acb(Arblib.ref(Q, n); prec = prec) : Acb(0; prec) for n in 0:J-1])
end

"""x^(a + b·ε) ... helper: exp((a + bε)·log x) for positive real ball x."""
function eps_pow(prec::Int, J::Int, x::Acb, a::Int, b::Int)
    lx = log(x)
    c = [Acb(0; prec) for _ in 1:J-1]
    c[1] = b * lx
    e = eps_exp(EpsSeries(1, c))
    xa = x^a
    EpsSeries(e.lead, [xa * ck for ck in e.c])
end

"""Evaluate coefficient of ε^k (absolute exponent)."""
function eps_coeff(x::EpsSeries, k::Int)
    idx = k - x.lead + 1
    1 <= idx <= nord(x) ? x.c[idx] : Acb(0; prec = prec_of(x))
end

"""Shift: multiply by ε^m."""
eps_shift(x::EpsSeries, m::Int) = EpsSeries(x.lead + m, copy(x.c))
