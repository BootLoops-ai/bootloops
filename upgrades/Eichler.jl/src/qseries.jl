# Ball-certified truncated q-series.
#
# A QSeries represents a function f(q) = Σ_{n≥0} a_n qⁿ analytic on |q| < 1 through
#   * ball enclosures of a_0..a_N  (coeffs, an AcbPoly), and
#   * a rigorous tail bound: Σ_{n>N} |a_n| x0ⁿ ≤ tailF,
# where x0 is a fixed reference radius shared by all series in a SeriesContext.
# Every operation preserves both invariants, so evaluation at any |q| ≤ rmax < x0
# returns an honest enclosure: |Σ_{n>N} a_n qⁿ| ≤ tailF · (|q|/x0)^{N+1}.

struct SeriesContext
    prec::Int
    N::Int          # truncation order: indices 0..N stored
    x0::Mag         # reference radius (must exceed every evaluation radius)
    function SeriesContext(prec::Int, N::Int, x0::Mag)
        # x0 < 1 is load-bearing: qsubst_pow uses x0^{dn} ≤ x0ⁿ and evaluate uses
        # (r/x0)ⁿ ≤ (r/x0)^{N+1} for r < x0 (audit finding: unvalidated x0 > 1 breaks both)
        Arblib.cmp(x0, Mag(1)) < 0 || error("SeriesContext: reference radius x0 must be < 1")
        new(prec, N, x0)
    end
end

SeriesContext(prec::Int, N::Int; x0::Real = 0.25) = SeriesContext(prec, N, Mag(x0))

struct QSeries
    ctx::SeriesContext
    coeffs::AcbPoly  # length ≤ N+1
    tailF::Mag       # Σ_{n>N} |a_n| x0ⁿ ≤ tailF
end

Base.length(f::QSeries) = Arblib.degree(f.coeffs) + 1

# Stored coefficient a_n, n ≤ N only (beyond the stored degree the truncated polynomial
# coefficient is exactly 0; for n > N the true series coefficient is NOT known — error).
qcoeff(f::QSeries, n::Int) = (n <= f.ctx.N || error("qcoeff: n exceeds truncation order");
    n <= Arblib.degree(f.coeffs) ? Acb(Arblib.ref(f.coeffs, n)) : Acb(0; prec = f.ctx.prec))

# ---- norm helpers ------------------------------------------------------------

# upper bound for Σ_{n=lo}^{hi} |a_n| x0ⁿ
function _mass(p::AcbPoly, x0::Mag, lo::Int, hi::Int)
    s = Mag(0)
    xp = Mag(1)
    one_ = Mag(1)
    # build x0^lo
    pw = Mag()
    Arblib.pow!(pw, x0, UInt(lo))
    xp = pw
    t = Mag()
    for n in lo:min(hi, Arblib.degree(p))
        Arblib.get!(t, Arblib.ref(p, n))   # mag upper bound of |a_n|
        Arblib.mul!(t, t, xp)
        Arblib.add!(s, s, t)
        Arblib.mul!(xp, xp, x0)
    end
    return s
end

# F(f) ≥ Σ_{n≥0} |a_n| x0ⁿ  (total majorant mass)
totalF(f::QSeries) = (s = _mass(f.coeffs, f.ctx.x0, 0, f.ctx.N); Arblib.add!(s, s, f.tailF); s)

# H(f) ≥ Σ_{n>N/2} |a_n| x0ⁿ  (upper-half mass, used for product tail bounds)
halfF(f::QSeries) = (s = _mass(f.coeffs, f.ctx.x0, div(f.ctx.N, 2) + 1, f.ctx.N); Arblib.add!(s, s, f.tailF); s)

# ---- constructors ------------------------------------------------------------

function qseries_zero(ctx::SeriesContext)
    QSeries(ctx, AcbPoly(prec = ctx.prec), Mag(0))
end

function qseries_const(ctx::SeriesContext, c)
    p = AcbPoly(prec = ctx.prec)
    Arblib.set_coeff!(p, 0, Acb(c; prec = ctx.prec))
    QSeries(ctx, p, Mag(0))
end

"""QSeries from exact Nemo integer/rational coefficients c[1] = a_0, ..., c[N+1] = a_N
together with a rigorous tail bound tailF on Σ_{n>N}|a_n|x0ⁿ."""
function qseries_exact(ctx::SeriesContext, c::Vector{QQFieldElem}, tailF::Mag)
    p = AcbPoly(prec = ctx.prec)
    for n in 0:min(length(c) - 1, ctx.N)
        num = Arblib.set!(Arb(prec = ctx.prec), BigInt(numerator(c[n+1])))
        den = Arblib.set!(Arb(prec = ctx.prec), BigInt(denominator(c[n+1])))
        Arblib.div!(num, num, den)
        Arblib.set_coeff!(p, n, Acb(num; prec = ctx.prec))
    end
    QSeries(ctx, p, tailF)
end

# ---- arithmetic --------------------------------------------------------------

function Base.:+(f::QSeries, g::QSeries)
    @assert f.ctx === g.ctx
    p = AcbPoly(prec = f.ctx.prec)
    Arblib.add!(p, f.coeffs, g.coeffs)
    t = Mag(); Arblib.add!(t, f.tailF, g.tailF)
    QSeries(f.ctx, p, t)
end

function Base.:-(f::QSeries, g::QSeries)
    @assert f.ctx === g.ctx
    p = AcbPoly(prec = f.ctx.prec)
    Arblib.sub!(p, f.coeffs, g.coeffs)
    t = Mag(); Arblib.add!(t, f.tailF, g.tailF)
    QSeries(f.ctx, p, t)
end

function Base.:*(c::Acb, f::QSeries)
    p = AcbPoly(prec = f.ctx.prec)
    Arblib.mul!(p, f.coeffs, c)
    cm = Mag(); Arblib.get!(cm, c)
    t = Mag(); Arblib.mul!(t, f.tailF, cm)
    QSeries(f.ctx, p, t)
end
Base.:*(f::QSeries, c::Acb) = c * f
Base.:*(c::Union{Integer,Rational}, f::QSeries) = Acb(c; prec = f.ctx.prec) * f

function Base.:*(f::QSeries, g::QSeries)
    @assert f.ctx === g.ctx
    ctx = f.ctx
    p = AcbPoly(prec = ctx.prec)
    Arblib.mullow!(p, f.coeffs, g.coeffs, ctx.N + 1)
    # tail(fg) ≤ H(f)F(g) + F(f)H(g): every lost term a_j b_k (j+k>N) has j>N/2 or k>N/2,
    # and tail-of-f times anything is also covered by H(f)F(g).
    Ff, Fg, Hf, Hg = totalF(f), totalF(g), halfF(f), halfF(g)
    t1 = Mag(); Arblib.mul!(t1, Hf, Fg)
    t2 = Mag(); Arblib.mul!(t2, Ff, Hg)
    Arblib.add!(t1, t1, t2)
    QSeries(ctx, p, t1)
end

Base.:^(f::QSeries, k::Integer) = Base.power_by_squaring(f, k)
Base.one(f::QSeries) = qseries_const(f.ctx, 1)
Base.copy(f::QSeries) = QSeries(f.ctx, AcbPoly(f.coeffs), Mag(f.tailF))

"""Multiply by qᵐ (shift). Tail: Σ_{n+m>N}|a_n|x0^{n+m} = x0^m·Σ_{n>N-m}|a_n|x0ⁿ
≤ x0^m·(stored mass above N-m + tailF)."""
function qshift(f::QSeries, m::Int)
    m >= 0 || error("qshift: m must be nonnegative")
    ctx = f.ctx
    p = AcbPoly(prec = ctx.prec)
    for n in 0:min(Arblib.degree(f.coeffs), ctx.N - m)
        Arblib.set_coeff!(p, n + m, Acb(Arblib.ref(f.coeffs, n)))
    end
    xm = Mag(); Arblib.pow!(xm, ctx.x0, UInt(m))
    t = _mass(f.coeffs, ctx.x0, max(ctx.N - m + 1, 0), ctx.N)
    Arblib.add!(t, t, f.tailF)
    Arblib.mul!(t, t, xm)
    QSeries(ctx, p, t)
end

"""Substitute q -> q^d (for E_k(d·τ) style operators). d ≥ 1."""
function qsubst_pow(f::QSeries, d::Int)
    d == 1 && return f
    ctx = f.ctx
    p = AcbPoly(prec = ctx.prec)
    for n in 0:min(Arblib.degree(f.coeffs), div(ctx.N, d))
        Arblib.set_coeff!(p, d * n, Acb(Arblib.ref(f.coeffs, n)))
    end
    # coefficients a_n now live at dn; those with dn > N plus old tail (x0^{dn} ≤ x0ⁿ since x0<1)
    t = _mass(f.coeffs, ctx.x0, div(ctx.N, d) + 1, ctx.N)  # bound using x0ⁿ ≥ x0^{dn}
    Arblib.add!(t, t, f.tailF)
    QSeries(ctx, p, t)
end

"""q d/dq of the series, returned on the slightly smaller reference radius x0' = θ·x0
(θ = 9/10). Rigor: Σ_{n>N} n|a_n| x0'ⁿ = Σ_{n>N} (n θⁿ)|a_n| x0ⁿ ≤ (N+1)θ^{N+1}·tailF,
since n θⁿ is decreasing for n ≥ 1/ln(1/θ) ≈ 9.49, hence for all n ≥ N+1 when N ≥ 10
(asserted). Downstream evaluations must satisfy |q| < 0.9·x0."""
function qdq(f::QSeries)
    ctx = f.ctx
    @assert ctx.N >= 11
    θ = Mag(0.9)
    x0p = Mag(); Arblib.mul!(x0p, ctx.x0, θ)          # rounds UP
    # Audit fix: the stored radius x0p is up-rounded, so the effective ratio
    # θ_eff = x0p/x0 may exceed θ. Recompute an upper bound of the actual ratio
    # (Mag division rounds up) and use IT in the tail bound; certify monotonicity
    # via θ_ub ≤ 0.91 < exp(-1/(N+1)) for N ≥ 11.
    θub = Mag(); Arblib.div!(θub, x0p, ctx.x0)
    Arblib.cmp(θub, Mag(0.91)) < 0 || error("qdq: rounded radius ratio too large")
    ctx2 = SeriesContext(ctx.prec, ctx.N, x0p)
    p = AcbPoly(prec = ctx.prec)
    for n in 1:Arblib.degree(f.coeffs)
        c = Acb(Arblib.ref(f.coeffs, n); prec = ctx.prec)
        Arblib.mul!(c, c, n)
        Arblib.set_coeff!(p, n, c)
    end
    t = Mag(); Arblib.pow!(t, θub, UInt(ctx.N + 1))
    Arblib.mul!(t, t, Mag(ctx.N + 1))
    Arblib.mul!(t, t, f.tailF)
    QSeries(ctx2, p, t)
end

# ---- evaluation --------------------------------------------------------------

"""Evaluate with certified tail at q (an Acb ball with |q| rigorously < x0)."""
function evaluate(f::QSeries, q::Acb)
    ctx = f.ctx
    qm = Mag(); Arblib.get!(qm, q)
    # require |q| < x0
    Arblib.is_finite(q) || error("evaluate: q not finite")
    ratio = Mag(); Arblib.div!(ratio, qm, ctx.x0)
    Arblib.cmp(ratio, Mag(1)) < 0 || error("evaluate: |q| must be < x0")
    val = Acb(prec = ctx.prec)
    Arblib.evaluate!(val, f.coeffs, q)
    # tail: tailF · ratio^{N+1}
    t = Mag(); Arblib.pow!(t, ratio, UInt(ctx.N + 1))
    Arblib.mul!(t, t, f.tailF)
    Arblib.add_error!(val, t)
    return val
end

"""Evaluate Σ_{n≥1} a_n qⁿ/n (the regularized primitive ∫₀^q (f - a_0) dq'/q' ) — helper."""
function evaluate_primitive(f::QSeries, q::Acb)
    ctx = f.ctx
    p = AcbPoly(prec = ctx.prec)
    for n in 1:Arblib.degree(f.coeffs)
        c = Acb(Arblib.ref(f.coeffs, n); prec = ctx.prec)
        Arblib.div!(c, c, n)
        Arblib.set_coeff!(p, n, c)
    end
    g = QSeries(ctx, p, f.tailF)  # |a_n/n| ≤ |a_n| for n ≥ 1
    evaluate(g, q)
end
