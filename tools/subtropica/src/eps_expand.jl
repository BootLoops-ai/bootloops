# eps_expand.jl — ε-expansion (C8). Phase A = DEGENERATE
# form: single face, no subtraction — loggamma_series(Γ-prefactor) ×
# exp(εΣb·logP) multinomial. Coefficients are polynomials in Log[P_j] with
# rational-function coefficients — NO general Series machinery
# (spec: STfastEpsSeries, reference/SubTropica.wl:17380 — the short-name
# polynomial trick, wl:17404-17418: base polys stay OPAQUE `localPolyName[i]`
# tokens during Series; here the analogue is that Log[P_j] are opaque HAtom-
# style tokens and the ε-expansion is the exact exp-multinomial closed form).
# Footguns owned here: truncated-Laurent pole guard (R2: RATIO-gate that
# requested order + padding covers the deepest pole; known pitfall: a
# truncated Laurent series eats its top orders), sign-definiteness rule (RT-viii).
#
# PREFACTOR SEPARATION (matches upstream): the .wl keeps the Γ-prefactor
# OUTSIDE the per-order integrands and pads maxOrder by the prefactor pole
# depth (STIntegrateSubtractionNP, wl:11185-11187:
#   minOrder = |vars(last face)| − |vars(first face)|,
#   maxOrder = order + Max[0, −(min ε-order of prefactor)] );
# assembly (C12) convolves gamma_prefactor_series with the integrated
# per-order results. eps_expand therefore returns the POLY-PART LogIntegrand
# lists (prefactor excluded) plus the prefactor Laurent for the assembly layer (lr_dispatch).
#
# NOTE on hlogexpr's hlog_add/hlog_mul/hlog_scale: those own GENERAL HlogExpr
# arithmetic (rational-FUNCTION string coefs). Everything this file needs is
# the PURE-RATIONAL-coefficient subalgebra (γ_E/ζ/log-atoms with exact
# Rational{BigInt} coefs), computed in a private accumulator (_RExpr below)
# and converted to HlogExpr only at the API boundary — no reimplementation
# of hlogexpr's string-coef semantics.

# ---------------------------------------------------------------------------
# Typed refusals
# ---------------------------------------------------------------------------

"""
    EpsExpandRefusal(msg)

Typed refusal from the C8 ε-expansion layer (unsupported Γ-argument base,
non-integer constant exponent, Euclidean-positivity violation, pole-guard
trip, exact Γ-pole, ...). Never a silent skip (RT-ix discipline).
"""
struct EpsExpandRefusal <: Exception
    msg::String
end
Base.showerror(io::IO, e::EpsExpandRefusal) = print(io, "EpsExpandRefusal: ", e.msg)

"""
    SignIndefiniteError(factor, msg)

Typed refusal from the RT-viii sign-definiteness rule: the named polynomial
factor changes sign on the declared positive domain, so the (a·b)^(c·ε)
split is REFUSED — never silently proceed. `factor` is the offending factor
as a string.
"""
struct SignIndefiniteError <: Exception
    factor::String
    msg::String
end
Base.showerror(io::IO, e::SignIndefiniteError) =
    print(io, "SignIndefiniteError on factor `", e.factor, "`: ", e.msg)

# ---------------------------------------------------------------------------
# Private rational-coefficient symbolic-constant algebra (_RExpr)
# key = canonically sorted Vector{Pair{HAtom,Int}}; value = exact rational.
# ---------------------------------------------------------------------------

const _AtomKey = Vector{Pair{HAtom,Int}}
const _RExpr   = Dict{_AtomKey,Rational{BigInt}}

# Canonical atom ordering/merging is OWNED BY src/hlogexpr.jl (
# included before this file): _canon_atomvec (sort+merge+drop-zeros),
# _atomvec_lt, _atom_lt. Reused here — NOT reimplemented.

_akey(pairs::Pair...) = _canon_atomvec(Pair{HAtom,Int}[pairs...])
const _EMPTYKEY = _AtomKey()

function _radd!(e::_RExpr, key::_AtomKey, c::Rational{BigInt})
    iszero(c) && return e
    nc = get(e, key, Rational{BigInt}(0)) + c
    if iszero(nc)
        delete!(e, key)
    else
        e[key] = nc
    end
    return e
end

_rexpr(c::Rational{BigInt}) = iszero(c) ? _RExpr() : _RExpr(_EMPTYKEY => c)

function _rexpr_addto!(a::_RExpr, b::_RExpr, scale::Rational{BigInt}=Rational{BigInt}(1))
    for (k, c) in b
        _radd!(a, k, c * scale)
    end
    return a
end

function _rexpr_mul(a::_RExpr, b::_RExpr)
    out = _RExpr()
    for (ka, ca) in a, (kb, cb) in b
        _radd!(out, _canon_atomvec(vcat(ka, kb)), ca * cb)
    end
    return out
end

_rexpr_scale(a::_RExpr, q::Rational{BigInt}) =
    iszero(q) ? _RExpr() : _RExpr(k => c * q for (k, c) in a)

# ε-series of _RExpr: S[i] = coefficient of ε^(i-1) relative to a shift.
function _rser_mul(A::Vector{_RExpr}, B::Vector{_RExpr}, L::Int)
    C = [_RExpr() for _ in 1:L]
    for i in 1:min(L, length(A)), j in 1:min(L - i + 1, length(B))
        _rexpr_addto!(C[i + j - 1], _rexpr_mul(A[i], B[j]))
    end
    return C
end

"exp of an _RExpr series with EMPTY constant slot (Sd[1] == 0): term_m = term_{m-1}·Sd/m"
function _rser_exp_clean(Sd::Vector{_RExpr}, L::Int)
    @assert isempty(Sd[1]) "internal: _rser_exp needs zero constant term"
    E    = [_RExpr() for _ in 1:L]
    E[1] = _rexpr(Rational{BigInt}(1))
    term = [_RExpr() for _ in 1:L]
    term[1] = _rexpr(Rational{BigInt}(1))
    for m in 1:(L - 1)
        term = _rser_mul(term, Sd, L)
        term = [_rexpr_scale(t, Rational{BigInt}(1, m)) for t in term]
        for i in 1:L
            _rexpr_addto!(E[i], term[i])
        end
    end
    return E
end

# Pure-rational power series u[i] = coeff of ε^(i-1)
function _rsmul(u::Vector{Rational{BigInt}}, v::Vector{Rational{BigInt}}, L::Int)
    w = zeros(Rational{BigInt}, L)
    for i in 1:min(L, length(u)), j in 1:min(L - i + 1, length(v))
        w[i + j - 1] += u[i] * v[j]
    end
    return w
end

function _rsinv(u::Vector{Rational{BigInt}}, L::Int)
    @assert !iszero(u[1]) "internal: series inverse needs nonzero constant"
    v    = zeros(Rational{BigInt}, L)
    v[1] = 1 // u[1]
    for k in 2:L
        s = Rational{BigInt}(0)
        for j in 2:k
            j <= length(u) || break
            s += u[j] * v[k - j + 1]
        end
        v[k] = -s // u[1]
    end
    return v
end

function _rspow(u::Vector{Rational{BigInt}}, p::Int, L::Int)
    p == 0 && return [i == 1 ? Rational{BigInt}(1) : Rational{BigInt}(0) for i in 1:L]
    p < 0 && return _rspow(_rsinv(u, L), -p, L)
    w = u[1:min(L, length(u))]
    length(w) < L && append!(w, zeros(Rational{BigInt}, L - length(w)))
    out = [i == 1 ? Rational{BigInt}(1) : Rational{BigInt}(0) for i in 1:L]
    for _ in 1:p
        out = _rsmul(out, w, L)
    end
    return out
end

# _RExpr -> HlogExpr (API boundary). Coef strings are content-canonical
# rationals (CONTRACTS.md §d: gcd=1, q>0, sign on numerator — automatic for
# Rational{BigInt}); terms sorted by atom part (canonical kind order).
_ratstr(q::Rational{BigInt}) =
    denominator(q) == 1 ? string(numerator(q)) : string(numerator(q)) * "/" * string(denominator(q))

function _to_hlog(e::_RExpr; vars::Vector{Symbol}=Symbol[])
    keys_sorted = sort!(collect(keys(e)); lt = _atomvec_lt)
    terms = HTerm[]
    for k in keys_sorted
        c = e[k]
        iszero(c) && continue
        push!(terms, HTerm(_ratstr(c), k))
    end
    return HlogExpr(terms, vars)
end

# ---------------------------------------------------------------------------
# log of a positive rational as atoms: log(p/q) = Σ e_i log(prime_i) − ...
# prime 2 → HConst(:Log2) (CONTRACTS.md §d: Log2 token = Log[2]);
# other primes → HLogA("<prime>"). log Γ(1/2) = ½·HLogA("Pi") (our atom,
# never sent to HF — CONTRACTS.md §d last bullet).
# ---------------------------------------------------------------------------
function _addlog_rational!(e::_RExpr, r::Rational{BigInt}, scale::Rational{BigInt})
    r > 0 || throw(EpsExpandRefusal("log of non-positive rational $r"))
    (iszero(scale) || isone(r)) && return e
    for (v, sgn) in ((numerator(r), 1), (denominator(r), -1))
        v == 1 && continue
        for (p, ex) in factor(ZZ(v))
            atom = p == 2 ? HConst(:Log2) : HLogA(string(p))
            _radd!(e, _akey(atom => 1), scale * sgn * ex)
        end
    end
    return e
end

# log(r + b·ε) as an _RExpr series (length L), r > 0 rational:
#   log r + Σ_{k≥1} (−1)^{k+1} (b/r)^k ε^k / k
function _addlog_linear!(S::Vector{_RExpr}, r::Rational{BigInt}, b::Rational{BigInt}, L::Int)
    _addlog_rational!(S[1], r, Rational{BigInt}(1))
    iszero(b) && return S
    for k in 1:(L - 1)
        _radd!(S[k + 1], _EMPTYKEY, (-1)^(k + 1) * (b // r)^k // k)
    end
    return S
end

# ---------------------------------------------------------------------------
# loggamma core: ε-series of log Γ(a + b·ε), a > 0, 2a ∈ ℤ  (RT-vii).
# ψ-expansion around the base β ∈ {1, 1/2} plus the closed downward
# recursion log Γ(a+x) = log Γ(a−1+x) + log(a−1+x):
#   β = 1  : log Γ(1+x) = −γ_E x + Σ_{k≥2} (−1)^k ζ(k) x^k / k
#   β = 1/2: log Γ(1/2+x) = ½logπ + (−γ_E − 2log2) x
#            + Σ_{k≥2} (−1)^k (2^k − 1) ζ(k) x^k / k
# (ψ(1/2) = −γ_E − 2 log 2; ψ^{(n)}(1/2) = (−1)^{n+1} n! (2^{n+1}−1) ζ(n+1).)
# Other rational bases need Gauss-digamma atoms (π·cot, log sin) outside the
# HlogExpr atom set — typed refusal (Phase-A limitation, noted in STATUS).
# ---------------------------------------------------------------------------
function _loggamma_rseries(a::Rational{BigInt}, b::Rational{BigInt}, L::Int)
    a > 0 || throw(EpsExpandRefusal(
        "loggamma_series: argument constant part a=$a ≤ 0 — Γ-poles/negative shifts " *
        "are the caller's job (gamma_prefactor_series), not loggamma_series"))
    denominator(a) in (1, 2) || throw(EpsExpandRefusal(
        "loggamma_series: base a=$a unsupported — only integer/half-integer Γ-arguments " *
        "have ψ-values in the {γ_E, ζ(k), Log2, Log[Pi]} atom set (Phase A)"))
    β = denominator(a) == 1 ? Rational{BigInt}(1) : Rational{BigInt}(1, 2)
    S = [_RExpr() for _ in 1:L]
    aa = a
    while aa > β                       # closed downward recursion
        aa -= 1
        _addlog_linear!(S, aa, b, L)   # + log(aa + b·ε)
    end
    if β == 1
        L >= 2 && _radd!(S[2], _akey(HConst(:EulerGamma) => 1), -b)
        for k in 2:(L - 1)
            _radd!(S[k + 1], _akey(HZeta([k]) => 1), Rational{BigInt}((-1)^k) * b^k // k)
        end
    else
        _radd!(S[1], _akey(HLogA("Pi") => 1), Rational{BigInt}(1, 2))
        if L >= 2
            _radd!(S[2], _akey(HConst(:EulerGamma) => 1), -b)
            _radd!(S[2], _akey(HConst(:Log2) => 1), -2b)
        end
        for k in 2:(L - 1)
            _radd!(S[k + 1], _akey(HZeta([k]) => 1),
                   Rational{BigInt}((-1)^k) * (big(2)^k - 1) * b^k // k)
        end
    end
    return S
end

"""
    loggamma_series(arg::EpsExp, order::Int) -> Vector{HlogExpr}

Symbolic ε-series of log Γ(a+b·ε) through ε^order (coefficient of ε^k in
entry k+1), from the standard ψ-expansion at base β ∈ {1, 1/2} plus the
closed recursion log Γ(a+x) = log Γ(a−1+x) + log(a−1+x). Coefficients are
HlogExpr with HConst(:EulerGamma), HConst(:Log2), HZeta([k]), HLogA atoms
(RT-vii: EXPLICITLY C8's job — no Julia library gives Γ-series symbolically
in γ_E/ζ(k)). Requirements: a > 0 and 2a ∈ ℤ, else typed EpsExpandRefusal;
Γ-poles (a ≤ 0) surface via gamma_prefactor_series' pole bookkeeping, not
here. Unit-tested vs mpmath loggamma at ≥50 dps (test/test_expand.jl).
"""
function loggamma_series(arg::EpsExp, order::Int)
    order >= 0 || throw(EpsExpandRefusal("loggamma_series: order=$order < 0"))
    S = _loggamma_rseries(arg.a, arg.b, order + 1)
    return [_to_hlog(e) for e in S]
end

# ---------------------------------------------------------------------------
# Γ-prefactor Laurent (RT-vii + pole bookkeeping).
# Γ(a+bε) = M(ε) · Γ(β+bε) with M a finite product of (r+bε) factors,
#   a > β: M = ∏_{r=β..a−1} (r+bε)        (Γ(z+1) = zΓ(z) upward)
#   a < β: M = 1/∏_{r=a..β−1} (r+bε)      (r = 0 term ⇒ the 1/(bε) pole)
# M is handled as an exact rational ε-series with an explicit ε-power shift;
# only log Γ(β+bε) enters the transcendental exp. exp(S0) with
# S0 = Σ p_i·logΓ(β_i) = q_π·log π must have q_π ∈ ℤ (else the constant is
# π^{1/2}-irrational, unrepresentable as an HlogExpr term — typed refusal).
# ---------------------------------------------------------------------------

"""
    gamma_prefactor_series(pf::Prefactor, order::Int)
        -> (minorder::Int, coeffs::Vector{HlogExpr})

Full symbolic Laurent of a Prefactor c·ε^k·∏Γ(g_i)^{p_i}·e^{nγ_Eε}:
`coeffs[i]` = coefficient of ε^(minorder+i−1), i = 1..order+1. Explicit pole
bookkeeping: minorder = k + Σ_i p_i·(#{r=0 shift factors of Γ_i}·(−1))
(Γ-poles at non-positive-integer arguments), exact by construction (leading
coefficient provably nonzero). Built on the loggamma core (RT-vii).
Typed refusals: Γ at an exact pole (a ≤ 0 integer, b = 0), non-half-integer
bases, net π^{1/2} constants, c = 0.
"""
function gamma_prefactor_series(pf::Prefactor, order::Int)
    order >= 0 || throw(EpsExpandRefusal("gamma_prefactor_series: order=$order < 0"))
    iszero(pf.c) && throw(EpsExpandRefusal("gamma_prefactor_series: zero prefactor c=0"))
    L = order + 1

    # rational multiplicative part: c · ∏ M_i^{p_i}, with ε-power shift
    shift = pf.eps_power
    R     = [i == 1 ? pf.c : Rational{BigInt}(0) for i in 1:L]
    # transcendental additive part: Σ p_i logΓ(β_i+b_iε) + n·γ_E·ε
    S = [_RExpr() for _ in 1:L]
    L >= 2 && _radd!(S[2], _akey(HConst(:EulerGamma) => 1), pf.gammaE_eps)

    for (g, p) in pf.gammas
        p == 0 && continue
        a, b = g.a, g.b
        denominator(a) in (1, 2) || throw(EpsExpandRefusal(
            "gamma_prefactor_series: Γ($(a)+$(b)ε) base unsupported (need 2a ∈ ℤ, Phase A)"))
        β = denominator(a) == 1 ? Rational{BigInt}(1) : Rational{BigInt}(1, 2)
        # M as (mshift, mser)
        mshift = 0
        mser   = [i == 1 ? Rational{BigInt}(1) : Rational{BigInt}(0) for i in 1:L]
        if a > β
            r = β
            while r < a
                lin = zeros(Rational{BigInt}, L)
                lin[1] = r
                L >= 2 && (lin[2] = b)
                mser = _rsmul(mser, lin, L)
                r += 1
            end
        elseif a < β
            r = a
            while r < β
                if iszero(r)
                    iszero(b) && throw(EpsExpandRefusal(
                        "gamma_prefactor_series: Γ($(a)) at an exact pole (b=0)"))
                    mshift -= 1
                    mser = [x // b for x in mser]
                else
                    lin = zeros(Rational{BigInt}, L)
                    lin[1] = r
                    L >= 2 && (lin[2] = b)
                    mser = _rsmul(mser, _rsinv(lin, L), L)
                end
                r += 1
            end
        end
        shift += mshift * p
        R = _rsmul(R, _rspow(mser, p, L), L)
        # base logΓ(β + bε), weight p
        base = _loggamma_rseries(β, b, L)
        for i in 1:L
            _rexpr_addto!(S[i], base[i], Rational{BigInt}(p))
        end
    end

    # constant slot S0: only ½·log π entries possible; exp(S0) = π^{q_π}
    qpi = Rational{BigInt}(0)
    for (k, c) in S[1]
        (length(k) == 1 && k[1].first == HLogA("Pi") && k[1].second == 1) ||
            throw(EpsExpandRefusal("gamma_prefactor_series: unexpected constant atom in S0: $k"))
        qpi += c
    end
    denominator(qpi) == 1 || throw(EpsExpandRefusal(
        "gamma_prefactor_series: prefactor constant is π^($qpi) — irrational (odd net " *
        "half-integer Γ count), not representable as an HlogExpr term (Phase A refusal)"))
    E0key = iszero(qpi) ? _EMPTYKEY : _akey(HConst(:Pi) => Int(qpi))

    Sd = [i == 1 ? _RExpr() : S[i] for i in 1:L]
    T  = _rser_exp_clean(Sd, L)
    # apply E0 (π^{q_π}) and the rational series R
    F = [_RExpr() for _ in 1:L]
    for i in 1:L, j in 1:(L - i + 1)
        iszero(R[i]) && continue
        for (k, c) in T[j]
            _radd!(F[i + j - 1], _canon_atomvec(vcat(k, E0key)), c * R[i])
        end
    end
    return (shift, [_to_hlog(e) for e in F])
end

# ---------------------------------------------------------------------------
# RT-viii sign-definiteness rule (MANDATORY concrete form, DESIGN.md R1)
# ---------------------------------------------------------------------------

"""
    SignSplitResult

Result of `sign_split`: `P = sign · content · ∏ factors[i][1]^factors[i][2]`
with every stored factor SIGN-POSITIVE on the declared domain, `content > 0`
rational, `sign ∈ {+1,−1}`. `evidence[i] ∈ (:coefficient_signs,
:numeric_probe)` records which RT-viii tier certified factor i;
`provenance` carries the probe points/signs when tier 2 was used.
"""
struct SignSplitResult
    sign::Int
    content::Rational{BigInt}
    factors::Vector{Tuple{QQMPolyRingElem,Int}}
    evidence::Vector{Symbol}
    provenance::Dict{String,Any}
end

function _rand_domain_point(dom::Union{Nothing,Tuple{Rational{BigInt},Rational{BigInt}}})
    if dom === nothing                      # (0, ∞) positive default
        return Rational{BigInt}(rand(1:(2^20)), rand(1:(2^20)))
    else
        lo, hi = dom
        t = Rational{BigInt}(rand(1:(2^20 - 1)), 2^20)
        return lo + (hi - lo) * t           # strict interior
    end
end

function _normdomains(domains, vars, kinvars)
    out = Dict{Symbol,Tuple{Rational{BigInt},Rational{BigInt}}}()
    for (k, v) in domains
        sym = Symbol(k)
        lo  = Rational{BigInt}(v[1])
        hi  = Rational{BigInt}(v[2])
        (0 <= lo < hi) || throw(EpsExpandRefusal(
            "sign_split: domain $sym=($lo,$hi) invalid — need 0 ≤ lo < hi (positive-kinematics contract)"))
        out[sym] = (lo, hi)
    end
    return out
end

"""
    sign_split(P::QQMPolyRingElem, vars, kinvars;
               nprobe::Int=16, domains=Dict()) -> SignSplitResult

RT-viii sign-definiteness rule for (a·b)^(c·ε) splitting — MANDATORY
implementation (DESIGN.md R1; upstream analogue: the Factor-based
coeff/mon/pol split of STtoCoeffMonPols, reference/SubTropica.wl:10311).
After Nemo `factor`, for EACH polynomial factor:
 1. coefficient-sign test (`coefficients`/`exponent_vectors`): all
    coefficients one sign ⇒ sign-definite under positive kinematics
    (valid for any positive sub-domain), evidence `:coefficient_signs`;
 2. else numeric probe at `nprobe` (default 16, ≥8 per R1) random
    positive-rational points drawn FRESH each run from each variable's
    domain (`domains[var] = (lo,hi)` with 0 ≤ lo < hi; absent = (0,∞)):
    consistent sign ⇒ allowed, probed points+signs recorded in
    `provenance`, evidence `:numeric_probe`;
 3. ANY sign change ⇒ typed `SignIndefiniteError` NAMING the factor —
    never silently proceed. The ledger covers polynomial factors, not just
    numeric leading coefficients (an (x−y)²-type input REFUSES here).
Negative-definite factors are stored sign-flipped (positive) with the flip
folded into `sign`. Self-check: the returned decomposition reproduces P.
"""
function sign_split(P::QQMPolyRingElem, vars, kinvars; nprobe::Int=16, domains=Dict{Symbol,Any}())
    iszero(P) && throw(EpsExpandRefusal("sign_split: zero polynomial"))
    nprobe >= 8 || throw(EpsExpandRefusal("sign_split: nprobe=$nprobe < 8 (R1 floor)"))
    Rg    = parent(P)
    doms  = _normdomains(domains, vars, kinvars)
    ringsyms = [Symbol(string(s)) for s in symbols(Rg)]
    fac   = factor(P)
    u     = unit(fac)
    is_constant(u) || error("internal: non-constant unit from Nemo factor")
    ucoef = BigInt(numerator(leading_coefficient(u))) // BigInt(denominator(leading_coefficient(u)))
    sgn   = ucoef > 0 ? 1 : -1
    content = abs(ucoef)
    factors  = Tuple{QQMPolyRingElem,Int}[]
    evidence = Symbol[]
    prov     = Dict{String,Any}("nprobe" => nprobe, "probes" => Dict{String,Any}())

    for (f, m) in fac
        if is_constant(f)
            c = BigInt(numerator(leading_coefficient(f))) // BigInt(denominator(leading_coefficient(f)))
            content *= abs(c)^m
            c < 0 && isodd(m) && (sgn = -sgn)
            continue
        end
        # tier 1: coefficient-sign test (monomials positive on any positive domain)
        cs   = collect(coefficients(f))
        σ    = 0
        if all(c -> c > 0, cs)
            σ = 1
        elseif all(c -> c < 0, cs)
            σ = -1
        end
        if σ != 0
            push!(evidence, :coefficient_signs)
        else
            # tier 2: numeric probe on the declared domains, fresh randomness
            pts   = Vector{Vector{Rational{BigInt}}}()
            signs = Int[]
            for _ in 1:nprobe
                local val
                local pt
                ok = false
                for _ in 1:64
                    pt  = [_rand_domain_point(get(doms, s, nothing)) for s in ringsyms]
                    val = evaluate(f, [QQ(x) for x in pt])
                    if !iszero(val)
                        ok = true
                        break
                    end
                end
                ok || throw(EpsExpandRefusal("sign_split: factor `$f` vanished at 64 straight probe points"))
                push!(pts, pt)
                push!(signs, val > 0 ? 1 : -1)
            end
            if !all(==(signs[1]), signs)
                throw(SignIndefiniteError(string(f),
                    "sign change on declared domain across $nprobe probe points " *
                    "(RT-viii tier 3 — split refused; signs=$(signs))"))
            end
            σ = signs[1]
            prov["probes"][string(f)] = Dict("points" => pts, "signs" => signs)
            push!(evidence, :numeric_probe)
        end
        g = σ == 1 ? f : -f
        σ == -1 && isodd(m) && (sgn = -sgn)
        push!(factors, (g, m))
    end

    # self-check: exact reconstruction
    recon = Rg(QQ(sgn * content))
    for (g, m) in factors
        recon *= g^m
    end
    recon == P || error("internal: sign_split reconstruction mismatch")
    return SignSplitResult(sgn, content, factors, evidence, prov)
end

# ---------------------------------------------------------------------------
# Laurent pole guard (R2; known pitfall: truncated-Laurent pole padding)
# ---------------------------------------------------------------------------

"""
    laurent_pole_guard(prefactor_minorder::Int, requested_order::Int,
                       maxorder::Int) -> Nothing

R2 truncated-Laurent pole guard: verifies the poly-part expansion depth
`maxorder` covers requested_order + the prefactor pole padding
max(0, −prefactor_minorder) (wl:11186-11187 maxOrder rule). A too-shallow
maxorder means the deepest prefactor pole EATS the top requested orders
(a known pitfall) — typed refusal, never a silent truncation.
"""
function laurent_pole_guard(prefactor_minorder::Int, requested_order::Int, maxorder::Int)
    need = requested_order + max(0, -prefactor_minorder)
    maxorder >= need || throw(EpsExpandRefusal(
        "truncated-Laurent pole guard (R2): maxorder=$maxorder < requested $requested_order " *
        "+ prefactor pole padding $(max(0, -prefactor_minorder)) — deepest pole would eat " *
        "the top ε-orders"))
    return nothing
end

# ---------------------------------------------------------------------------
# Degenerate C8: the exp(εΣb·logP) multinomial expansion
# ---------------------------------------------------------------------------

"compositions of `total` into `parts` ordered nonnegative integers"
function _compositions(total::Int, parts::Int)
    parts == 0 && return total == 0 ? [Int[]] : Vector{Int}[]
    out = Vector{Int}[]
    for first in 0:total, rest in _compositions(total - first, parts - 1)
        push!(out, vcat([first], rest))
    end
    return out
end

function _int_exponent(q::Rational{BigInt}, what::String)
    denominator(q) == 1 || throw(EpsExpandRefusal(
        "eps_expand: $what has non-integer constant part $q — the ε⁰ layer would not be " *
        "a rational function (Phase A refusal)"))
    return Int(numerator(q))
end

"""
    eps_expand(quad::EulerIntegrand, order::Int;
               domains=Dict(), nprobe::Int=16)
        -> (minorder=0, integrands::Vector{Vector{LogIntegrand}},
            prefactor_minorder::Int, prefactor_coeffs::Vector{HlogExpr},
            maxorder::Int, requested_order::Int)

Degenerate-C8 Phase A expansion of the POLY PART ∏x^ν·∏P^(a+bε) into
per-ε-order LogIntegrand lists; `integrands[o+1]` = the ε^o coefficient
(minorder = 0 for the single degenerate face — wl:11185
minOrder = |vars(last face)| − |vars(first face)|). The Γ-prefactor is kept
SEPARATE exactly as upstream does (STIntegrateSubtractionNP, wl:11185-11187:
prefactor multiplies the integrated result; returned here as
`prefactor_minorder`/`prefactor_coeffs` for the C12 convolution in lr_dispatch), and
the poly part is padded to maxorder = order + max(0, −prefactor_minorder)
(wl:11186-11187), checked by `laurent_pole_guard` (R2).

Closed form (STfastEpsSeries, wl:17380, short-name trick wl:17404-17418 —
Log[P] stays an opaque token, no general Series):
  ∏x^ν ∏P^(a+bε) = A(x) · exp(ε·Σ_t c_t·Log[g_t]),
  ε^o coefficient = A · Σ_{|k|=o} (∏_t c_t^{k_t}/k_t!) ∏_t Log[g_t]^{k_t},
with A the rational function from integer exponent parts and {g_t} the
sign-positive letters from `sign_split` (RT-viii; ε-independent integer
powers bypass the split — integer powers need no positivity). ε-dependent
jacobian factors (1−u)^(1−TropI) enter the poly list like any P_j
(RT-extra-jac; MANDATORY fixture in test_expand.jl uses `domains` to declare
u ∈ (0,1)). Zero ε-orders are padded as EMPTY term lists (padFromSeries
analogue, wl:17385-17402).
"""
function eps_expand(quad::EulerIntegrand, order::Int;
                    domains=Dict{Symbol,Any}(), nprobe::Int=16)
    order >= 0 || throw(EpsExpandRefusal("eps_expand: order=$order < 0"))
    Rg = quad.ring
    length(quad.nu) == length(quad.vars) || throw(EpsExpandRefusal(
        "eps_expand: |nu|=$(length(quad.nu)) ≠ |vars|=$(length(quad.vars))"))
    nvars(Rg) == length(quad.vars) + length(quad.kinvars) || throw(EpsExpandRefusal(
        "eps_expand: ring has $(nvars(Rg)) vars ≠ |vars|+|kinvars|"))
    for (P, _) in quad.polys
        parent(P) === Rg || throw(EpsExpandRefusal("eps_expand: poly not in quad.ring"))
    end

    # prefactor Laurent + pole padding (wl:11186-11187)
    pmin0, _ = gamma_prefactor_series(quad.prefactor, order)
    pordr    = order + max(0, -pmin0)
    pmin, pcoeffs = gamma_prefactor_series(quad.prefactor, pordr - pmin0)
    pmin == pmin0 || error("internal: prefactor minorder not reproducible")
    minorder = 0                        # single degenerate face (wl:11185)
    maxo     = order + max(0, -pmin)
    laurent_pole_guard(pmin, order, maxo)

    num  = one(Rg)
    den  = one(Rg)
    coefrat = Rational{BigInt}(1)
    letters = Tuple{QQMPolyRingElem,Rational{BigInt}}[]
    function addletter(g::QQMPolyRingElem, c::Rational{BigInt})
        iszero(c) && return
        i = findfirst(t -> t[1] == g, letters)
        if i === nothing
            push!(letters, (g, c))
        else
            nc = letters[i][2] + c
            iszero(nc) ? deleteat!(letters, i) : (letters[i] = (g, nc))
        end
    end

    # monomial part ∏ x_i^(ν_a + ν_b ε)
    for (i, ν) in enumerate(quad.nu)
        na = _int_exponent(ν.a, "nu[$i]")
        na > 0 && (num *= gen(Rg, i)^na)
        na < 0 && (den *= gen(Rg, i)^(-na))
        iszero(ν.b) || addletter(gen(Rg, i), ν.b)
    end

    # poly part ∏ P_j^(a+bε)
    for (j, (P, e)) in enumerate(quad.polys)
        a = _int_exponent(e.a, "polys[$j] exponent")
        if iszero(e.b)
            # ε-independent INTEGER power: rational function regardless of sign,
            # no split, no positivity requirement.
            a > 0 && (num *= P^a)
            a < 0 && (den *= P^(-a))
            continue
        end
        if is_constant(P)
            c = BigInt(numerator(leading_coefficient(P))) // BigInt(denominator(leading_coefficient(P)))
            c > 0 || throw(EpsExpandRefusal(
                "eps_expand: constant poly factor $c ≤ 0 with ε-dependent exponent (Euclidean violation)"))
            coefrat *= c^a
            isone(c) || addletter(Rg(QQ(c)), e.b)
            continue
        end
        s = sign_split(P, quad.vars, quad.kinvars; nprobe = nprobe, domains = domains)
        s.sign == 1 || throw(EpsExpandRefusal(
            "eps_expand: polys[$j] = $P is NEGATIVE on the declared domain — violates the " *
            "declared-Euclidean positivity contract ((−F)^eF signs are the front-end's job)"))
        coefrat *= s.content^a
        isone(s.content) || addletter(Rg(QQ(s.content)), e.b)
        for (g, m) in s.factors
            ma = m * a
            ma > 0 && (num *= g^ma)
            ma < 0 && (den *= g^(-ma))
            addletter(g, m * e.b)
        end
    end

    # multinomial per ε-order (STfastEpsSeries closed form, wl:17380)
    T = length(letters)
    integrands = Vector{Vector{LogIntegrand}}()
    for o in 0:maxo
        terms = LogIntegrand[]
        for k in _compositions(o, T)
            q = coefrat
            logs = Pair{QQMPolyRingElem,Int}[]
            for t in 1:T
                k[t] == 0 && continue
                q *= letters[t][2]^k[t] // factorial(big(k[t]))
                push!(logs, letters[t][1] => k[t])
            end
            (T == 0 && o > 0) && continue
            push!(terms, LogIntegrand(num * QQ(q), den, logs,
                                      copy(quad.vars), copy(quad.kinvars)))
        end
        push!(integrands, terms)        # o with no terms => EMPTY list (explicit padding)
    end

    return (minorder = minorder, integrands = integrands,
            prefactor_minorder = pmin, prefactor_coeffs = pcoeffs,
            maxorder = maxo, requested_order = order)
end
