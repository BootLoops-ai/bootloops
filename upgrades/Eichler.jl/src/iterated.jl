# Iterated integrals of modular forms, I(f1,...,fn; q), with the tangential base point
# regularization of AW (docs/conventions.md §4): constant terms integrate to powers of
# L = log q (for q = q_C this means L := 2πi·τ_C, NOT the principal branch).
#
# Representation: V(q) = Σ_{j=0}^{J} F_j(q) · L^j, F_j certified QSeries, plus a list of
# rigorous remainder records (A, j) encoding |R(q)| ≤ A·(|q|/x0)^{N+1}·(1+|L(q)|)^j along
# the integration ray q' = qσ, σ ∈ (0,1]; note |L(qσ)| = |L(q) + ln σ| ≤ |L(q)| + ln(1/σ).

struct LogQSeries
    levels::Vector{QSeries}            # levels[j+1] = F_j
    rem::Vector{Tuple{Mag,Int}}        # remainder records (A, j)
end

logqseries_one(ctx::SeriesContext) = LogQSeries([qseries_const(ctx, 1)], Tuple{Mag,Int}[])

"""S_{N,J} = Σ_{i=0}^{J} C(J,i)·i!/(N+1)^{i+1}  (path-integration constant, Mag upper)."""
function _SNJ(prec::Int, N::Int, J::Int)
    s = Arb(0; prec)
    for i in 0:J
        s += Arb(binomial(big(J), big(i)) * factorial(big(i)); prec) / Arb(big(N + 1)^(i + 1); prec)
    end
    m = Mag(); Arblib.get!(m, Arblib.ubound(s))
    return m
end

"""W(q) = ∫₀^q f(q') V(q') dq'/q'  with tangential-base-point regularization.
f must be a plain QSeries on the same context."""
function integrate_kernel(f::QSeries, V::LogQSeries)
    ctx = f.ctx
    prec, N = ctx.prec, ctx.N
    J = length(V.levels) - 1
    # U_j = f * F_j
    U = [f * Fj for Fj in V.levels]
    # output levels up to J+1
    out = [qseries_zero(ctx) for _ in 0:(J+1)]
    Ff = totalF(f)
    newrem = Tuple{Mag,Int}[]
    # carry forward old remainders: each (A, j) -> (A·F_f·S_{N,j}, j)
    for (A, j) in V.rem
        A2 = Mag(); Arblib.mul!(A2, A, Ff)
        Arblib.mul!(A2, A2, _SNJ(prec, N, j))
        push!(newrem, (A2, j))
    end
    for j in 0:J
        Uj = U[j+1]
        # constant term -> L^{j+1}/(j+1)
        u0 = qcoeff(Uj, 0)
        c = Acb(u0; prec); Arblib.div!(c, c, j + 1)
        pj1 = out[j+2]
        addc = qseries_const(ctx, c)
        out[j+2] = pj1 + addc
        # monomials n ≥ 1: ∫ qⁿ L^j dq/q = qⁿ Σ_{i=0}^{j} (-1)^i j!/(j-i)! n^{-(i+1)} L^{j-i}
        for i in 0:j
            p = AcbPoly(prec = prec)
            fac = (-1)^i * factorial(big(j)) ÷ factorial(big(j - i))
            for n in 1:Arblib.degree(Uj.coeffs)
                cn = Acb(Arblib.ref(Uj.coeffs, n); prec = prec)
                Arblib.mul!(cn, cn, Acb(fac; prec = prec))
                Arblib.div!(cn, cn, Acb(big(n)^(i + 1); prec = prec))
                Arblib.set_coeff!(p, n, cn)
            end
            # coefficient tail: |c_n·fac/n^{i+1}| ≤ |c_n|·j!  (since fac/n^{i+1} ≤ j!)
            tF = Mag(); Arblib.mul!(tF, Uj.tailF, Mag(factorial(big(j))))
            g = QSeries(ctx, p, tF)
            out[j-i+1] = out[j-i+1] + g
        end
        # remainder created by the coefficient tail of U_j under ∫·L^j:
        A = Mag(); Arblib.mul!(A, Uj.tailF, _SNJ(prec, N, j))
        push!(newrem, (A, j))
    end
    LogQSeries(out, newrem)
end

"""Evaluate V at the ball q, with L = logq supplied explicitly (= 2πi τ_C)."""
function evaluate(V::LogQSeries, q::Acb, logq::Acb)
    ctx = V.levels[1].ctx
    prec = ctx.prec
    val = Acb(0; prec)
    Lp = Acb(1; prec)
    for (j, Fj) in enumerate(V.levels)
        val += evaluate(Fj, q) * Lp
        Lp *= logq
    end
    # remainders
    qm = Mag(); Arblib.get!(qm, q)
    ratio = Mag(); Arblib.div!(ratio, qm, ctx.x0)
    rpow = Mag(); Arblib.pow!(rpow, ratio, UInt(ctx.N + 1))
    Lm = Mag(); Arblib.get!(Lm, logq)
    Lp1 = Mag(); Arblib.add!(Lp1, Lm, Mag(1))
    err = Mag(0)
    for (A, j) in V.rem
        e = Mag(); Arblib.pow!(e, Lp1, UInt(j))
        Arblib.mul!(e, e, A)
        Arblib.mul!(e, e, rpow)
        Arblib.add!(err, err, e)
    end
    Arblib.add_error!(val, err)
    return val
end

"""I(f₁,...,f_k; q) as a LogQSeries (fold from the innermost letter). Letters are
QSeries; pass `one` for the constant kernel 1."""
function iterated_integral(letters::Vector{QSeries})
    ctx = letters[1].ctx
    V = logqseries_one(ctx)
    for f in reverse(letters)
        V = integrate_kernel(f, V)
    end
    return V
end

"""Convenience: certified value of I(letters...; q) with logq = 2πiτ."""
function iterated_integral_value(letters::Vector{QSeries}, q::Acb, logq::Acb)
    evaluate(iterated_integral(letters), q, logq)
end
