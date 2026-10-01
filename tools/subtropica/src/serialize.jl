# serialize.jl — Julia -> Mma-grammar strings for HF
# requests (CONTRACTS.md §d; input grammar of HF parse_expr,
# reference/hyperflint_cli_main.cpp:2593-2632).
#
# ROUND-TRIP LAW: every emitted string must round-trip through
# {"op":"parse_expr"} unchanged-modulo-canonicalization (standing unit test
# in test_bridge.jl; live examples: "-x^2" -> "-x^2",
# "x^-2" -> "1/x^2", Log[x+1] -> Hlog[x + 1,[0]] — canonicalization is
# allowed to rewrite, so the test gates on parse SUCCESS + semantic
# equality via the `eval` op at rational points, not on byte equality).
#
# Grammar rules implemented here (CONTRACTS.md §d):
#  - rationals CONTENT-CANONICAL: gcd(p,q)=1, q>0, sign on numerator
#    (Rat::parse does NOT normalize non-coprime content — handlers.cpp T3
#    caveat; Julia `Rational` and Nemo `QQFieldElem` are already canonical,
#    we only format them).
#  - explicit `*` always, `^` for powers (binds tighter than unary minus —
#    guaranteed by the v1.2.8 canary, CONTRACTS.md §a).
#  - identifiers [A-Za-z][A-Za-z0-9_]*, NO backtick contexts ever (upstream
#    strips them before emission, SubTropica.wl:12337-12345 — we validate
#    instead of stripping).
#  - Log[expr] only (rewritten internally to Hlog[arg,{0}]); Hlog[z,{...}];
#    never PolyLog.
#  - delta[var] is an OUTPUT-side token only: live probe shows
#    parse_expr REJECTS "delta[x4]" ("unexpected token '['"), so the
#    serializer refuses HDelta (it must never be sent toward HF).
#  - Our own symbolic atoms (EulerGamma, HZeta zeta tokens, HPeriod,
#    HAlgLetter) exist only INSIDE HlogExpr and are never emitted toward
#    HF (CONTRACTS.md §d) — mma_string(::HlogExpr) throws on them.

const _MMA_IDENT_RE = r"^[A-Za-z][A-Za-z0-9_]*$"

# SubTropicaSerializeError and the Log[rational] → prime-log canonical fold
# (`fold_log_atoms` + helpers) MOVED to hlogexpr.jl (canonical-atom law;
# canonical-atom law): `canonical_text` now folds at entry, so the fold
# must live where canonical_text lives — including in standalone contexts
# that include hlogexpr.jl WITHOUT serialize.jl (test_b1_parity.jl).
# This file still throws/uses both (same module; hlogexpr.jl included first).

"""
    mma_string(x) -> String

Serialize `x` to the HyperFLINT Mma input grammar (CONTRACTS.md §d).
Methods: `Integer`, `Rational`, `QQFieldElem`, `Symbol` (validated
identifier), `QQMPolyRingElem` (explicit `*`/`^`, deterministic Nemo term
order), `LogIntegrand` ((num/den)·∏Log[P]^k — flat-measure integrand of
CONTRACTS.md §b), `HlogExpr` (HF-parseable subset only: HConst(:Pi)/(:I),
HConst(:Log2) — emitted as `Log[2]`, the canonical Log2 atom decision —
HLogA, HHlog; throws a `SubTropicaSerializeError` on atoms HF must never
see: HZeta/HPeriod/HAlgLetter/HDelta/other HConst).

Round-trip law: every string this function emits must round-trip through
HF `parse_expr` (unit-tested in test_bridge.jl).
"""
function mma_string end

# ---- scalars ---------------------------------------------------------------

mma_string(n::Integer) = string(n)

# Julia Rational is content-canonical by construction (gcd=1, den>0).
function mma_string(q::Rational)
    denominator(q) == 1 && return string(numerator(q))
    return string(numerator(q)) * "/" * string(denominator(q))
end

# Nemo rational (canonical by construction; Nemo prints `p//q`, we need p/q).
function mma_string(q::QQFieldElem)
    d = denominator(q)
    d == 1 && return string(numerator(q))
    return string(numerator(q)) * "/" * string(d)
end

function mma_string(s::Symbol)
    st = string(s)
    occursin(_MMA_IDENT_RE, st) ||
        throw(SubTropicaSerializeError("symbol '$st' is not a valid Mma-grammar " *
            "identifier [A-Za-z][A-Za-z0-9_]* (backtick contexts are never " *
            "emitted — SubTropica.wl:12337-12345)"))
    return st
end

# ---- polynomials -----------------------------------------------------------

# One monomial "x^2*y" from an exponent vector (empty -> "").
function _mono_string(varnames::Vector{String}, expv::Vector{Int})
    parts = String[]
    for (v, k) in zip(varnames, expv)
        k == 0 && continue
        push!(parts, k == 1 ? v : v * "^" * string(k))
    end
    return join(parts, "*")
end

# Term as (negative?, magnitude-string) so sums join as "a - b" not "a+-b".
function _term_repr_poly(c::QQFieldElem, mono::String)
    neg = c < 0
    ca = neg ? -c : c
    if isempty(mono)
        return (neg, mma_string(ca))
    elseif ca == 1
        return (neg, mono)
    else
        return (neg, mma_string(ca) * "*" * mono)
    end
end

function _join_terms(reprs::Vector{Tuple{Bool,String}})
    isempty(reprs) && return "0"
    io = IOBuffer()
    for (i, (neg, s)) in enumerate(reprs)
        if i == 1
            neg && print(io, "-")
            print(io, s)
        else
            print(io, neg ? " - " : " + ")
            print(io, s)
        end
    end
    return String(take!(io))
end

"""
    mma_string(p::QQMPolyRingElem) -> String

Polynomial in the HF grammar: explicit `*`, `^`, content-canonical rational
coefficients, deterministic Nemo term order. Ring symbol names must match
the request's `vars` (caller's responsibility; names are validated as
Mma identifiers).
"""
function mma_string(p::QQMPolyRingElem)
    iszero(p) && return "0"
    varnames = [string(v) for v in symbols(parent(p))]
    for v in varnames
        occursin(_MMA_IDENT_RE, v) ||
            throw(SubTropicaSerializeError("ring variable '$v' is not a valid Mma identifier"))
    end
    reprs = Tuple{Bool,String}[]
    for (c, ev) in zip(coefficients(p), exponent_vectors(p))
        push!(reprs, _term_repr_poly(c, _mono_string(varnames, Vector{Int}(ev))))
    end
    return _join_terms(reprs)
end

# Fraction of two QQMPolyRingElem (Nemo frac field elements are already
# unit-normalized: monic denominator, sign on the numerator — verified:
# 2x/(4x^2) -> (1/2)/x, x/(-2y) -> (-1/2*x)/y).
# Emitted standalone-unambiguous: den==1 -> num; else "(num)/(den)".
function _frac_string(fr)
    num = numerator(fr)
    den = denominator(fr)
    isone(den) && return mma_string(num)
    return "(" * mma_string(num) * ")/(" * mma_string(den) * ")"
end

# True if `s` cannot be safely used as a factor in a `*` product without
# wrapping in parentheses: any top-level (bracket-depth-0) '+' or '-'
# beyond a leading unary minus. ("1/2", "-3*x" are safe: */ are
# left-associative at equal precedence; "s + t" is not.)
function _mul_unsafe(s::AbstractString)
    depth = 0
    for (i, c) in enumerate(s)
        if c in ('(', '[', '{')
            depth += 1
        elseif c in (')', ']', '}')
            depth -= 1
        elseif depth == 0 && i > 1 && (c == '+' || c == '-')
            return true
        end
    end
    return false
end

_paren_if_unsafe(s::AbstractString) = _mul_unsafe(s) ? "(" * s * ")" : String(s)

# ---- LogIntegrand ----------------------------------------------------------

"""
    mma_string(li::LogIntegrand) -> String

The locally-finite handoff integrand (CONTRACTS.md §b) as an HF `expr`
string: `(num)/(den)*Log[P_1]^k_1*...`. The measure is flat dx on
[0,∞)^|vars| and is NOT part of the string; integration order/vars are
passed separately in the request (`vars_int`). den==1 omits the division;
k==1 omits `^1`.
"""
function mma_string(li::LogIntegrand)
    parts = String[]
    num_s = mma_string(li.num)
    if isone(li.den)
        push!(parts, _paren_if_unsafe(num_s))
    else
        push!(parts, "(" * num_s * ")/(" * mma_string(li.den) * ")")
    end
    for (P, k) in li.logs
        k >= 1 || throw(SubTropicaSerializeError(
            "LogIntegrand log power must be ≥1 (types.jl contract), got $k"))
        s = "Log[" * mma_string(P) * "]"
        push!(parts, k == 1 ? s : s * "^" * string(k))
    end
    return join(parts, "*")
end

# ---- Log[rational] → prime-log basis fold: MOVED to hlogexpr.jl ----------
# (canonical-atom law: canonical_text now folds at
#  entry, hlogexpr.jl owns fold_log_atoms + helpers; mma_string below
#  still folds via the same function.)

# ---- HlogExpr (HF-parseable subset only) ------------------------------------

const _HF_EMITTABLE_CONSTS = (:Pi, :I)   # HF-native special symbols (CONTRACTS.md §d)

function _atom_mma_string(a::HAtom)
    if a isa HConst
        a.name in _HF_EMITTABLE_CONSTS && return string(a.name)
        # canonical Log2 atom : emit as Log[2] — the
        # grammar form parse_expr accepts (the bare "Log2" token is an
        # HF OUTPUT spelling, not guaranteed as parser input).
        a.name === :Log2 && return "Log[2]"
        throw(SubTropicaSerializeError("HConst(:$(a.name)) is an internal atom " *
            "and must never be emitted toward HF (CONTRACTS.md §d)"))
    elseif a isa HLogA
        return "Log[" * a.arg * "]"
    elseif a isa HHlog
        return "Hlog[" * a.z * ",{" * join(a.word, ",") * "}]"
    else
        throw(SubTropicaSerializeError("$(typeof(a)) atoms are internal " *
            "(HZeta/HPeriod/HAlgLetter) or output-only (HDelta — live probe " *
            "parse_expr rejects delta[..]); never emitted toward HF"))
    end
end

"""
    mma_string(e::HlogExpr) -> String

Emit the HF-parseable subset of an HlogExpr (coef · Pi/I/Log/Hlog powers).
Throws `SubTropicaSerializeError` on any atom HF must never see
(HZeta/HPeriod/HAlgLetter/HDelta/non-Pi-I HConst) — CONTRACTS.md §d.

Canonical-emission fold : `Log[rational]` atoms
are folded into the Log2/Log-prime basis first (`fold_log_atoms`), so
identically-zero combinations emit as the literal "0".
"""
function mma_string(e::HlogExpr)
    e = fold_log_atoms(e)
    isempty(e.terms) && return "0"
    reprs = Tuple{Bool,String}[]
    for t in e.terms
        c = t.coef
        neg = false
        if startswith(c, "-") && !_mul_unsafe(c[2:end])
            neg = true
            c = c[2:end]
        end
        factors = String[]
        (c == "1" && !isempty(t.atoms)) || push!(factors, _paren_if_unsafe(c))
        for (a, k) in t.atoms
            k == 0 && continue
            s = _atom_mma_string(a)
            push!(factors, k == 1 ? s : s * "^" * string(k))
        end
        push!(reprs, (neg, join(factors, "*")))
    end
    return _join_terms(reprs)
end
