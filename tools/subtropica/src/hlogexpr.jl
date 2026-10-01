# hlogexpr.jl — HF coef-string -> HlogExpr token parser +
# canonical text form + exact HlogExpr arithmetic + numeric verification
# hooks. Structs live in types.jl — this file is LOGIC ONLY.
#
# Contracts: CONTRACTS.md §a (token dictionary), §d (grammar, canonical
# ordering). Token semantics transliterated from the vendored
# reference/SubTropica.wl coef-string parser stHyperFlintCoefStringToMma
# (wl:12261-12321):
#   mzv_a_b_...  -> HZeta([a,b,...]), 'm' prefix = negative/alternating
#                   index (stMzvTokenToZeta, wl:12190-12202; longest token
#                   replaced first upstream, wl:12311-12316 — our tokenizer
#                   scans maximal identifiers so no prefix-eating issue)
#   Log2         -> HConst(:Log2)                  (wl:12318-12319; canonical
#                   Log2 ATOM DECISION: the ONE
#                   canonical atom for log 2 is HConst(:Log2) — the MZV table
#                   speaks "Log2", types.jl names it, loggamma_series emits it.
#                   HLogA("2") from `Log[2]` strings FOLDS to HConst(:Log2) at
#                   every canonicalization chokepoint (_norm_atom below), so
#                   the two spellings can never coexist in one term list.)
#   Wm_<i>/Wp_<i>-> HAlgLetter(:Wm/:Wp, i)         (wl:12287-12297; a token
#                   without a letter-table entry is a WIRING BUG upstream —
#                   STHyperFlint::algletter — and a typed error here)
#   WmOverWp_<i> -> Wm_i * Wp_i^-1                 (wl:12274-12285, replaced
#                   BEFORE Wm/Wp upstream to avoid prefix-eating)
#   sqrt_disc_<i>-> HAlgLetter(:sqrt_disc, i)      (wl:12299-12309; upstream
#                   substitutes Sqrt[<disc literal>] — we keep the atom and
#                   the AlgLetter table carries `disc`)
#   delta[var]   -> HDelta(var) (contour residues; OUTPUT-side only)
#   Pi, I        -> HConst; Log[q] -> HLogA; Hlog[z,{..}] -> HHlog
#
# Wm/Wp NORMALIZATION goes through HF's OWN ops (DESIGN.md RT-v — never a
# local Vieta reimplementation): see `normalize_wmwp` below. Only the
# token->struct PARSE lives here. Live-pinned op division of labor
# (v1.2.8 binary, poly u^2+3u+1):
#   simplify_with_vieta:    Wm_1*Wp_1 -> 1 (products/sums via Vieta;
#                           leaves Wm_1+Wp_1 alone but back_substitute
#                           maps it -> -3)
#   combine_wm_wp_ratios:   Wm_1/Wp_1 -> WmOverWp_1 (literal ratios ONLY —
#                           it does NOT recombine sqrt-form ratios, so it
#                           must run BEFORE back_substitute)
#   back_substitute:        Wm_1-Wp_1 -> -sqrt_disc_1; Wm_1+Wp_1 -> -3;
#                           Wm_1/Wp_1 -> (-1/2*sqrt_disc_1-3/2)/(1/2*..)
#                           (an atom-in-denominator RATIO our HlogExpr
#                           cannot represent — the reason for the order)
# => default chain: simplify_with_vieta, combine_wm_wp_ratios,
#    back_substitute.
#
# ZeroInfPeriod semantics (pinned by live probe + the T1-T3 smoke
# integrals): ZeroInfPeriod[{l1,...,ln}] = the shuffle-regularized constant
# term (log z -> 0) of G(l1,...,ln; z) as z -> +inf, same letter order,
# letter l = pole at t=l (dt/(t-l)). Verified: ["0","-1"]->mzv_2,
# ["0","-1","-1"]->mzv_3, ["-1"]->0, ["-2"]->-Log2,
# ["0","-2"]->1/2*Log2^2+mzv_2, ["0","-1","-2"]->-Log2*mzv_2+7/8*mzv_3 —
# all matched ≥60 dps against the independent ginac_gpl large-z
# extrapolation implemented in `zero_inf_period_ginac`.
# NOTE: CONTRACTS.md §a's example "['0','-1'] -> 1/2*mzv_2" does NOT match
# the live v1.2.8 binary (-> "mzv_2"); per reference/PROVENANCE.md the live
# probe wins. Recorded as a known upstream-example discrepancy.

# ---------------------------------------------------------------------------
# Typed errors (fail loudly: never return a default on failure)
# ---------------------------------------------------------------------------

struct SubTropicaParseError <: Exception
    msg::String
    input::String
end
Base.showerror(io::IO, e::SubTropicaParseError) =
    print(io, "SubTropicaParseError: ", e.msg, "\n  while parsing: ", e.input)

struct SubTropicaEvalError <: Exception
    msg::String
end
Base.showerror(io::IO, e::SubTropicaEvalError) =
    print(io, "SubTropicaEvalError: ", e.msg)

# Serialization-side typed error. Defined HERE (not serialize.jl) since
# `fold_log_atoms` (canonical-atom law, below) throws it and
# lives in this file so `canonical_text` and standalone contexts that
# include hlogexpr.jl WITHOUT serialize.jl (test_b1_parity.jl) carry the
# same canonical emission. serialize.jl (included after) keeps using it.
struct SubTropicaSerializeError <: Exception
    msg::String
end
Base.showerror(io::IO, e::SubTropicaSerializeError) =
    print(io, "SubTropicaSerializeError: ", e.msg)

# ---------------------------------------------------------------------------
# Canonical atom order (CONTRACTS.md §d: HConst by name, HZeta by index
# vector, HLogA by arg string, HHlog, HPeriod, HAlgLetter by (kind,idx),
# HDelta by var)
# ---------------------------------------------------------------------------

_atom_rank(::HConst) = 1
_atom_rank(::HZeta) = 2
_atom_rank(::HLogA) = 3
_atom_rank(::HHlog) = 4
_atom_rank(::HPeriod) = 5
_atom_rank(::HAlgLetter) = 6
_atom_rank(::HDelta) = 7

const _ALG_KIND_RANK = Dict{Symbol,Int}(:Wm => 1, :Wp => 2, :sqrt_disc => 3)

# three-way compare helpers
_cmp(a::T, b::T) where {T<:Union{String,Symbol,Int}} = a < b ? -1 : (a > b ? 1 : 0)
function _cmp_vec(a, b)
    for i in 1:min(length(a), length(b))
        c = _cmp(a[i], b[i])
        c != 0 && return c
    end
    return _cmp(length(a), length(b))
end

function _atom_cmp(a::HAtom, b::HAtom)
    c = _cmp(_atom_rank(a), _atom_rank(b))
    c != 0 && return c
    if a isa HConst
        return _cmp(string(a.name), string(b.name))
    elseif a isa HZeta
        return _cmp_vec(a.idx, b.idx)
    elseif a isa HLogA
        return _cmp(a.arg, b.arg)
    elseif a isa HHlog
        c = _cmp(a.z, b.z); c != 0 && return c
        return _cmp_vec(a.word, b.word)
    elseif a isa HPeriod
        return _cmp_vec(a.word, b.word)
    elseif a isa HAlgLetter
        c = _cmp(_ALG_KIND_RANK[a.kind], _ALG_KIND_RANK[b.kind]); c != 0 && return c
        return _cmp(a.idx, b.idx)
    elseif a isa HDelta
        return _cmp(string(a.var), string(b.var))
    end
    return 0
end
_atom_lt(a::HAtom, b::HAtom) = _atom_cmp(a, b) < 0

function _atomvec_cmp(a::Vector{Pair{HAtom,Int}}, b::Vector{Pair{HAtom,Int}})
    for i in 1:min(length(a), length(b))
        c = _atom_cmp(a[i].first, b[i].first)
        c != 0 && return c
        c = _cmp(a[i].second, b[i].second)
        c != 0 && return c
    end
    return _cmp(length(a), length(b))
end
_atomvec_lt(a, b) = _atomvec_cmp(a, b) < 0

# Canonical-representation fold : log 2 has
# exactly ONE atom spelling, HConst(:Log2). HLogA("2") — produced by `Log[2]`
# in HF coef strings before this fold existed — normalizes here; applied in
# _hv_atom (parser entry) and _canon_atomvec (arithmetic/merge chokepoint),
# so parse_coef output and hlog_add/hlog_mul results never carry HLogA("2").
_norm_atom(a::HAtom) = (a isa HLogA && a.arg == "2") ? HConst(:Log2) : a

# Canonicalize an atom-power vector: normalize reps, sort, merge equal
# atoms, drop power 0.
function _canon_atomvec(av::Vector{Pair{HAtom,Int}})
    isempty(av) && return Pair{HAtom,Int}[]
    v = sort(Pair{HAtom,Int}[_norm_atom(p.first) => p.second for p in av],
             lt=(p, q) -> _atom_lt(p.first, q.first))
    out = Pair{HAtom,Int}[]
    for p in v
        if !isempty(out) && _atom_cmp(out[end].first, p.first) == 0
            np = out[end].second + p.second
            pop!(out)
            np != 0 && push!(out, p.first => np)
        else
            p.second != 0 && push!(out, p)
        end
    end
    return out
end

# Canonical text of a single atom (CONTRACTS.md §d fixed forms).
function _atom_string(a::HAtom)
    if a isa HConst
        return string(a.name)
    elseif a isa HZeta
        return "mzv_" * join([i < 0 ? "m" * string(-i) : string(i) for i in a.idx], "_")
    elseif a isa HLogA
        return "Log[" * a.arg * "]"
    elseif a isa HHlog
        return "Hlog[" * a.z * ",{" * join(a.word, ",") * "}]"
    elseif a isa HPeriod
        return "ZeroInfPeriod[{" * join(a.word, ",") * "}]"
    elseif a isa HAlgLetter
        return (a.kind === :Wm ? "Wm_" : a.kind === :Wp ? "Wp_" : "sqrt_disc_") * string(a.idx)
    elseif a isa HDelta
        return "delta[" * string(a.var) * "]"
    end
    error("unreachable atom kind $(typeof(a))")
end

# ---------------------------------------------------------------------------
# Parse context: one Nemo QQ[vars...] fraction field for the coefficient
# rational functions; the algebraic-letter table for Wm/Wp validation.
# ---------------------------------------------------------------------------

struct _PCtx
    vars::Vector{Symbol}
    F::Any                      # fraction_field(QQ[vars...])
    varfr::Dict{Symbol,Any}     # var -> field element
    letters::Vector{AlgLetter}
    # Wm/Wp denominator router (gap closure): a function
    # String -> String that sends a Rat::parse-SAFE "1/(<poly in letters>)"
    # body through HF's OWN RT-v chain and returns the engine's normalized
    # string. `nothing` (default) = the old behavior: multi-term atom
    # denominators refuse typed. Armed only by `normalize_wmwp`.
    wmwp_router::Any
end

function _pctx(vars::Vector{Symbol}, letters::Vector{AlgLetter}=AlgLetter[];
               wmwp_router=nothing)
    length(unique(vars)) == length(vars) ||
        throw(SubTropicaParseError("duplicate variable in $(vars)", string(vars)))
    R, gens_ = polynomial_ring(QQ, String[string(v) for v in vars])
    F = fraction_field(R)
    varfr = Dict{Symbol,Any}(v => F(gens_[i]) for (i, v) in enumerate(vars))
    return _PCtx(vars, F, varfr, letters, wmwp_router)
end

# Routerless variant (parse of router OUTPUT must never re-route: no recursion).
_pctx_norouter(ctx::_PCtx) = _PCtx(ctx.vars, ctx.F, ctx.varfr, ctx.letters, nothing)

# Internal value: sparse sum over atom-power vectors with exact
# rational-function coefficients. Keys are canonical atom vectors.
const _HVal = Dict{Vector{Pair{HAtom,Int}},Any}

_hv_zero() = _HVal()
function _hv_const(ctx::_PCtx, q)
    v = _HVal()
    f = ctx.F(base_ring(ctx.F)(q))   # coerce rational constant via QQ[vars...]
    iszero(f) || (v[Pair{HAtom,Int}[]] = f)
    return v
end
_hv_atom(ctx::_PCtx, a::HAtom) = _HVal(Pair{HAtom,Int}[_norm_atom(a) => 1] => ctx.F(1))

function _hv_add(x::_HVal, y::_HVal)
    out = _HVal()
    for (k, v) in x
        out[k] = v
    end
    for (k, v) in y
        nv = haskey(out, k) ? out[k] + v : v
        if iszero(nv)
            delete!(out, k)
        else
            out[k] = nv
        end
    end
    return out
end

_hv_neg(x::_HVal) = _HVal(k => -v for (k, v) in x)

function _hv_mul(x::_HVal, y::_HVal)
    out = _HVal()
    for (ka, va) in x, (kb, vb) in y
        k = _canon_atomvec(vcat(ka, kb))
        v = va * vb
        nv = haskey(out, k) ? out[k] + v : v
        if iszero(nv)
            delete!(out, k)
        else
            out[k] = nv
        end
    end
    return out
end

# x^k for integer k. Negative k inverts: allowed for the atom-free single
# term and for a SINGLE term with atoms (atom powers negate — this is how
# HF puts algebraic letters in denominators, e.g. 1/Wp_1). Multi-term sums
# containing atoms (the back_substitute ratio shapes) are handled by, in
# order (Wm/Wp gap closure): (i) the table-driven sqrt_disc
# even-power fold (sqrt_disc_i^2 = disc — the engine NEVER substitutes
# disc, live probes T2-T5 in CONTRACTS.md §d); (ii) exact conjugate
# rationalization of a + b*sqrt_disc_i denominators; (iii) the armed
# wmwp_router (HF's OWN RT-v chain on a Rat::parse-safe "1/(D)" body —
# normalize_wmwp). With no router armed and (i)/(ii) insufficient the
# typed refusal below stands unchanged.
function _hv_pow(ctx::_PCtx, x::_HVal, k::Int, src::String)
    if k == 0
        return _hv_const(ctx, 1)
    elseif k > 0
        acc = x
        for _ in 2:k
            acc = _hv_mul(acc, x)
        end
        return acc
    else
        if length(x) == 1
            (av, c) = first(x)
            iszero(c) && throw(SubTropicaParseError("division by zero", src))
            inv1 = _HVal([a => -p for (a, p) in av] => inv(c))
            return _hv_pow(ctx, inv1, -k, src)
        end
        isempty(x) && throw(SubTropicaParseError("division by zero", src))
        inv1 = _hv_invert_multiterm(ctx, x, src)
        return _hv_pow(ctx, inv1, -k, src)
    end
end

# ---------------------------------------------------------------------------
# sqrt_disc canonicalization + multi-term denominator inversion
# (Wm/Wp algebraic-letter gap closure; CONTRACTS.md §d diagnosis)
# ---------------------------------------------------------------------------

# Algebraic-letter tokens anywhere in a string (word letters / coef strings).
const _ALGLETTER_TOKEN_RE = r"(?:^|[^A-Za-z0-9_])(?:Wm|Wp|WmOverWp|sqrt_disc)_\d+"

# disc of letter idx as a ctx fraction, or nothing when unknown/unusable
# (stub tables from _letters_of carry disc="" — folding is then skipped;
# skipping is exactness-preserving, the atom power is simply kept).
function _letter_disc_frac(ctx::_PCtx, idx::Int)
    for l in ctx.letters
        if l.idx == idx
            isempty(strip(l.disc)) && return nothing
            return try
                _parse_rat_frac(l.disc, ctx)
            catch e
                e isa SubTropicaParseError ? nothing : rethrow()
            end
        end
    end
    return nothing
end

"""
    _hv_fold_sqrtdisc(ctx, hv) -> _HVal

Table-driven canonical fold sqrt_disc_i^k -> disc^((k - mod(k,2))/2) *
sqrt_disc_i^mod(k,2) (floor-mod: k=-1 -> disc^-1 * sqrt_disc_i). Mirrors
upstream, where Sqrt[<disc literal>] is substituted into Mathematica and
auto-simplifies; our HlogExpr keeps the atom, so the defining relation
sqrt_disc^2 = disc (the response table's own `disc` field) is applied
here. Letters with unknown disc are left untouched. Exactness-preserving.
"""
function _hv_fold_sqrtdisc(ctx::_PCtx, hv::_HVal)
    out = _hv_zero()
    for (av, c) in hv
        nav = Pair{HAtom,Int}[]
        nc = c
        for (a, k) in av
            if a isa HAlgLetter && a.kind === :sqrt_disc && (k >= 2 || k <= -1)
                d = _letter_disc_frac(ctx, a.idx)
                if d === nothing || iszero(d)
                    push!(nav, a => k)
                else
                    r = mod(k, 2)                    # 0 or 1, also for k<0
                    nc *= d^div(k - r, 2)
                    r == 1 && push!(nav, a => 1)
                end
            else
                push!(nav, a => k)
            end
        end
        one_term = _HVal(_canon_atomvec(nav) => nc)
        out = _hv_add(out, one_term)
    end
    return out
end

# Invert a MULTI-TERM _HVal (denominator). Called from _hv_pow with k<0.
function _hv_invert_multiterm(ctx::_PCtx, x::_HVal, src::String)
    xf = _hv_fold_sqrtdisc(ctx, x)
    length(xf) == 1 && return _hv_pow(ctx, xf, -1, src)
    isempty(xf) && throw(SubTropicaParseError("division by zero", src))
    # (ii) conjugate rationalization: shape a + b*sqrt_disc_i (single idx,
    # power 1, known disc). 1/(a+b*s) = (a-b*s)/(a^2-b^2*disc) — exact,
    # atoms-free denominator by construction.
    if length(xf) == 2
        ks = collect(keys(xf))
        emptyk = findfirst(isempty, ks)
        if emptyk !== nothing
            other = ks[emptyk == 1 ? 2 : 1]
            if length(other) == 1 && other[1].first isa HAlgLetter &&
               other[1].first.kind === :sqrt_disc && other[1].second == 1
                d = _letter_disc_frac(ctx, other[1].first.idx)
                if d !== nothing && !iszero(d)
                    a = xf[Pair{HAtom,Int}[]]
                    b = xf[other]
                    den = a * a - b * b * d
                    iszero(den) && throw(SubTropicaParseError(
                        "conjugate rationalization hit a zero norm " *
                        "a^2 - b^2*disc — degenerate sqrt_disc denominator", src))
                    out = _HVal()
                    q1 = a / den
                    iszero(q1) || (out[Pair{HAtom,Int}[]] = q1)
                    q2 = -b / den
                    iszero(q2) || (out[other] = q2)
                    return out
                end
            end
        end
    end
    # (iii) route the denominator through HF's OWN RT-v chain (armed only
    # by normalize_wmwp). Transport law (live-refuted alternatives in
    # CONTRACTS.md §d): ONLY the single-top-level-division body
    # "1/(<'/'-free polynomial in the letters>)" is Rat::parse-safe —
    # the engine mis-parses trailing *factors after a /(...) group
    # (silently moves them INTO the denominator) and in-band-refuses its
    # own parenthesized "(N/D)" emission.
    if ctx.wmwp_router !== nothing && _hv_routable(xf)
        (S, scale) = _hv_ratsafe_string(ctx, xf)
        routed = String(ctx.wmwp_router("1/(" * S * ")"))
        rhv = try
            _hv_fold_sqrtdisc(ctx, _parse_hval(routed, _pctx_norouter(ctx)))
        catch e
            e isa SubTropicaParseError || rethrow()
            throw(SubTropicaParseError("Wm/Wp denominator was routed through " *
                "HF's simplify_with_vieta/combine_wm_wp_ratios/" *
                "back_substitute chain but the NORMALIZED form still fails " *
                "to parse: '$(routed)' ($(e.msg))", src))
        end
        # 1/D = scale * parse("1/(S)") with S = scale*D
        sc = _HVal()
        sc[Pair{HAtom,Int}[]] = scale
        return _hv_mul(sc, rhv)
    end
    throw(SubTropicaParseError("cannot invert a multi-term expression " *
        "containing transcendental/algebraic atoms — route it through " *
        "HF's simplify_with_vieta/combine_wm_wp_ratios/back_substitute " *
        "first (normalize_wmwp; DESIGN.md RT-v)", src))
end

# Routable = every atom is an HAlgLetter with POSITIVE power (the engine
# var pool knows Wm/Wp/WmOverWp/sqrt_disc only; mzv/Log/... in a
# denominator is contract drift and keeps the typed refusal).
function _hv_routable(hv::_HVal)
    for (av, _) in hv, (a, k) in av
        (a isa HAlgLetter && k > 0) || return false
    end
    return true
end

"""
    _hv_ratsafe_string(ctx, hv) -> (S::String, scale)

Emit `hv` (atoms: HAlgLetter only, positive powers) as a '/'-free
polynomial string S over the letter tokens with INTEGER-coefficient
kinematic polys, together with the exact `scale` (ctx.F element) such
that hv = (1/scale)*S. Guarantees "1/(S)" is a single-top-level-division
body — the only Rat::parse-safe division shape (CONTRACTS.md §d).
"""
function _hv_ratsafe_string(ctx::_PCtx, hv::_HVal)
    R = base_ring(ctx.F)
    # clear kinematic denominators
    Q = one(R)
    for (_, c) in hv
        Q = lcm(Q, denominator(c))
    end
    polys = Vector{Any}()
    keysv = Vector{Vector{Pair{HAtom,Int}}}()
    L = one(ZZ)
    for (av, c) in hv
        p = numerator(c) * divexact(Q, denominator(c))
        for q in coefficients(p)
            L = lcm(L, ZZ(denominator(q)))
        end
        push!(polys, p)
        push!(keysv, av)
    end
    reprs = Tuple{Bool,String}[]
    for (p, av) in zip(polys, keysv)
        ps = mma_string(p * QQ(L))
        neg = false
        if startswith(ps, "-") && !_mul_unsafe(ps[2:end])
            neg = true
            ps = ps[2:end]
        end
        mono = join([k == 1 ? _atom_string(a) : _atom_string(a) * "^" * string(k)
                     for (a, k) in av], "*")
        s = if isempty(mono)
            _paren_if_unsafe(ps)
        elseif ps == "1"
            mono
        else
            _paren_if_unsafe(ps) * "*" * mono
        end
        push!(reprs, (neg, s))
    end
    return (_join_terms(reprs), ctx.F(Q * R(QQ(L))))
end

_hv_div(ctx::_PCtx, x::_HVal, y::_HVal, src::String) =
    _hv_mul(x, _hv_pow(ctx, y, -1, src))

# Require atom-free; return the fraction-field element.
function _hv_to_frac(ctx::_PCtx, x::_HVal, what::String, src::String)
    isempty(x) && return ctx.F(0)
    (length(x) == 1 && isempty(first(x).first)) ||
        throw(SubTropicaParseError("$what must be free of transcendental atoms", src))
    return first(x).second
end

# Require an integer constant (for exponents).
function _hv_to_int(ctx::_PCtx, x::_HVal, src::String)
    fr = _hv_to_frac(ctx, x, "exponent", src)
    num = numerator(fr); den = denominator(fr)
    (is_constant(num) && isone(den)) ||
        throw(SubTropicaParseError("exponent must be an integer constant", src))
    q = constant_coefficient(num)
    isone(denominator(q)) ||
        throw(SubTropicaParseError("exponent must be an integer, got $q", src))
    return Int(numerator(q))
end

# ---------------------------------------------------------------------------
# Tokenizer (ASCII grammar: numbers, identifiers, + - * / ^ ( ) [ ] { } ,)
# ---------------------------------------------------------------------------

function _tokenize(s::AbstractString)
    toks = Tuple{Symbol,String}[]
    i = firstindex(s)
    n = lastindex(s)
    while i <= n
        c = s[i]
        if isspace(c)
            i = nextind(s, i)
        elseif isdigit(c)
            j = i
            while j <= n && isdigit(s[j])
                j = nextind(s, j)
            end
            push!(toks, (:num, String(s[i:prevind(s, j)])))
            i = j
        elseif 'A' <= c <= 'Z' || 'a' <= c <= 'z'
            j = i
            while j <= n && (isdigit(s[j]) || 'A' <= s[j] <= 'Z' || 'a' <= s[j] <= 'z' || s[j] == '_')
                j = nextind(s, j)
            end
            push!(toks, (:id, String(s[i:prevind(s, j)])))
            i = j
        elseif c in ('+', '-', '*', '/', '^', '(', ')', '[', ']', '{', '}', ',')
            push!(toks, (:op, string(c)))
            i = nextind(s, i)
        else
            throw(SubTropicaParseError("unexpected character '$c' at index $i", String(s)))
        end
    end
    push!(toks, (:eof, ""))
    return toks
end

mutable struct _PS
    toks::Vector{Tuple{Symbol,String}}
    i::Int
    ctx::_PCtx
    src::String
end
_peek(ps::_PS) = ps.toks[ps.i]
_advance!(ps::_PS) = (t = ps.toks[ps.i]; ps.i += 1; t)
function _expect_op!(ps::_PS, o::String)
    t = _advance!(ps)
    (t[1] === :op && t[2] == o) ||
        throw(SubTropicaParseError("expected '$o', got '$(t[2])' (token $(ps.i - 1))", ps.src))
    return t
end

# ---------------------------------------------------------------------------
# Recursive-descent parser, Mathematica precedence:
#   sum < product(*,/) < unary minus < power(^, right-assoc)
# so -x^2 == -(x^2) (the v1.2.8 canary pin, CONTRACTS.md §a) and
# 1/2*x == (1/2)*x (left-assoc equal precedence).
# ---------------------------------------------------------------------------

function _parse_sum(ps::_PS)
    t = _peek(ps)
    neg = false
    if t[1] === :op && (t[2] == "+" || t[2] == "-")
        _advance!(ps)
        neg = t[2] == "-"
    end
    acc = _parse_product(ps)
    neg && (acc = _hv_neg(acc))
    while true
        t = _peek(ps)
        (t[1] === :op && (t[2] == "+" || t[2] == "-")) || break
        _advance!(ps)
        rhs = _parse_product(ps)
        acc = _hv_add(acc, t[2] == "+" ? rhs : _hv_neg(rhs))
    end
    return acc
end

function _parse_product(ps::_PS)
    acc = _parse_factor(ps)
    while true
        t = _peek(ps)
        (t[1] === :op && (t[2] == "*" || t[2] == "/")) || break
        _advance!(ps)
        rhs = _parse_factor(ps)
        acc = t[2] == "*" ? _hv_mul(acc, rhs) : _hv_div(ps.ctx, acc, rhs, ps.src)
    end
    return acc
end

function _parse_factor(ps::_PS)
    t = _peek(ps)
    if t[1] === :op && (t[2] == "-" || t[2] == "+")
        _advance!(ps)
        v = _parse_factor(ps)
        return t[2] == "-" ? _hv_neg(v) : v
    end
    return _parse_power(ps)
end

function _parse_power(ps::_PS)
    base = _parse_primary(ps)
    t = _peek(ps)
    if t[1] === :op && t[2] == "^"
        _advance!(ps)
        expv = _parse_factor(ps)          # right-assoc; allows ^-2, ^(3)
        k = _hv_to_int(ps.ctx, expv, ps.src)
        return _hv_pow(ps.ctx, base, k, ps.src)
    end
    return base
end

const _MZV_TOKEN_RE = r"^mzv(?:_m?\d+)+$"
const _WMP_TOKEN_RE = r"^(Wm|Wp)_(\d+)$"
const _WMOVERWP_RE = r"^WmOverWp_(\d+)$"
const _SQRTDISC_RE = r"^sqrt_disc_(\d+)$"

function _letter_check(ps::_PS, idx::Int, tok::String)
    any(l -> l.idx == idx, ps.ctx.letters) ||
        throw(SubTropicaParseError("algebraic-letter token '$tok' has no entry in " *
            "the AlgLetter table — wiring bug (upstream STHyperFlint::algletter, " *
            "SubTropica.wl:12289-12297); pass the response's algebraic_letters " *
            "table via `letters`", ps.src))
end

# Parse a letter list into canonical strings (atom-free). Accepts BOTH
# our §d brace form {a,b,...} and HF's canonical-AST bracket form [a,b,...]
# (live parse_expr pin: "Hlog[x + 1,[0]]").
function _parse_brace_list(ps::_PS)
    t = _advance!(ps)
    (t[1] === :op && (t[2] == "{" || t[2] == "[")) ||
        throw(SubTropicaParseError("expected '{' or '[' letter list, got '$(t[2])'", ps.src))
    close = t[2] == "{" ? "}" : "]"
    out = String[]
    if _peek(ps)[1] === :op && _peek(ps)[2] == close
        _advance!(ps)
        return out
    end
    while true
        v = _parse_sum(ps)
        push!(out, _frac_string(_hv_to_frac(ps.ctx, v, "letter", ps.src)))
        t = _advance!(ps)
        t[1] === :op || throw(SubTropicaParseError("bad letter list", ps.src))
        t[2] == close && return out
        t[2] == "," || throw(SubTropicaParseError("expected ',' or '$close' in letter list", ps.src))
    end
end

function _parse_primary(ps::_PS)
    t = _advance!(ps)
    if t[1] === :num
        return _hv_const(ps.ctx, parse(BigInt, t[2]))
    elseif t[1] === :op && t[2] == "("
        v = _parse_sum(ps)
        _expect_op!(ps, ")")
        return v
    elseif t[1] === :id
        name = t[2]
        nxt = _peek(ps)
        has_bracket = nxt[1] === :op && nxt[2] == "["
        if name == "Log" && has_bracket
            _advance!(ps)
            argv = _parse_sum(ps)
            _expect_op!(ps, "]")
            fr = _hv_to_frac(ps.ctx, argv, "Log argument", ps.src)
            isone(fr) && return _hv_const(ps.ctx, 0)   # Log[1] = 0
            return _hv_atom(ps.ctx, HLogA(_frac_string(fr)))
        elseif name == "Hlog" && has_bracket
            _advance!(ps)
            zv = _parse_sum(ps)
            zs = _frac_string(_hv_to_frac(ps.ctx, zv, "Hlog argument", ps.src))
            _expect_op!(ps, ",")
            word = _parse_brace_list(ps)
            _expect_op!(ps, "]")
            return _hv_atom(ps.ctx, HHlog(zs, word))
        elseif name == "ZeroInfPeriod" && has_bracket
            _advance!(ps)
            word = _parse_brace_list(ps)
            _expect_op!(ps, "]")
            isempty(word) &&
                throw(SubTropicaParseError("empty ZeroInfPeriod word never stored (types.jl)", ps.src))
            return _hv_atom(ps.ctx, HPeriod(word))
        elseif name == "delta" && has_bracket
            _advance!(ps)
            vt = _advance!(ps)
            vt[1] === :id || throw(SubTropicaParseError("delta[..] takes a symbol", ps.src))
            _expect_op!(ps, "]")
            return _hv_atom(ps.ctx, HDelta(Symbol(vt[2])))
        elseif has_bracket
            throw(SubTropicaParseError("unknown function '$name[..]' in HF coef string", ps.src))
        elseif name == "Pi"
            return _hv_atom(ps.ctx, HConst(:Pi))
        elseif name == "I"
            return _hv_atom(ps.ctx, HConst(:I))
        elseif name == "EulerGamma"
            return _hv_atom(ps.ctx, HConst(:EulerGamma))
        elseif name == "Catalan"
            return _hv_atom(ps.ctx, HConst(:Catalan))
        elseif name == "Log2"                       # wl:12318-12319
            return _hv_atom(ps.ctx, HConst(:Log2))  # canonical Log2 atom decision
        elseif occursin(_MZV_TOKEN_RE, name)        # wl:12190-12202, 12311-12316
            idx = Int[]
            for part in split(name, '_')[2:end]
                if startswith(part, "m")
                    push!(idx, -parse(Int, part[2:end]))
                else
                    push!(idx, parse(Int, part))
                end
            end
            return _hv_atom(ps.ctx, HZeta(idx))
        elseif (m = match(_WMOVERWP_RE, name)) !== nothing   # wl:12274-12285
            i = parse(Int, m.captures[1])
            _letter_check(ps, i, name)
            return _HVal(_canon_atomvec(Pair{HAtom,Int}[HAlgLetter(:Wm, i) => 1,
                                        HAlgLetter(:Wp, i) => -1]) => ps.ctx.F(1))
        elseif (m = match(_WMP_TOKEN_RE, name)) !== nothing  # wl:12287-12297
            i = parse(Int, m.captures[2])
            _letter_check(ps, i, name)
            return _hv_atom(ps.ctx, HAlgLetter(m.captures[1] == "Wm" ? :Wm : :Wp, i))
        elseif (m = match(_SQRTDISC_RE, name)) !== nothing   # wl:12299-12309
            i = parse(Int, m.captures[1])
            _letter_check(ps, i, name)
            return _hv_atom(ps.ctx, HAlgLetter(:sqrt_disc, i))
        elseif haskey(ps.ctx.varfr, Symbol(name))
            v = _HVal()
            v[Pair{HAtom,Int}[]] = ps.ctx.varfr[Symbol(name)]
            return v
        else
            throw(SubTropicaParseError("unknown symbol '$name' — not a declared " *
                "kinematic variable and not a recognized HF token (declare it " *
                "in `vars`, or this is contract drift: print-loud refusal)", ps.src))
        end
    end
    throw(SubTropicaParseError("unexpected token '$(t[2])'", ps.src))
end

function _parse_hval(s::AbstractString, ctx::_PCtx)
    ps = _PS(_tokenize(s), 1, ctx, String(s))
    v = _parse_sum(ps)
    _peek(ps)[1] === :eof ||
        throw(SubTropicaParseError("trailing input at token $(ps.i): '$(_peek(ps)[2])'", ps.src))
    return v
end

# Parse a coef string that must be a pure rational function (no atoms).
function _parse_rat_frac(s::AbstractString, ctx::_PCtx)
    v = _parse_hval(s, ctx)
    return _hv_to_frac(ctx, v, "coefficient", String(s))
end

function _hv_to_hlogexpr(hv::_HVal, ctx::_PCtx)
    terms = HTerm[]
    for (av, fr) in hv
        iszero(fr) && continue
        push!(terms, HTerm(_frac_string(fr), copy(av)))
    end
    sort!(terms, lt=(a, b) -> _atomvec_lt(a.atoms, b.atoms))
    return HlogExpr(terms, copy(ctx.vars))
end

# ---------------------------------------------------------------------------
# Public parsing / canonical text API
# ---------------------------------------------------------------------------

"""
    parse_coef(s::AbstractString, vars::Vector{Symbol};
               letters::Vector{AlgLetter}=AlgLetter[]) -> HlogExpr

Parse one HF response coef string (`Rat::to_string` format: rational
function of kinematics + `Pi`, `I`, `Log[n]`/`Log2`, `delta[var]`,
`mzv_a_b_c` — 'm' = negative/alternating index — `Wm_i`/`Wp_i`/
`WmOverWp_i`/`sqrt_disc_i`, and `Hlog[z,{..}]`) into an `HlogExpr`.
Canonical-rep decision: both `Log2` and `Log[2]` produce
HConst(:Log2) — HLogA("2") never appears in parse output (`_norm_atom`).
Token dictionary per SubTropica.wl:12261-12321 (see file header).
`letters` = the response's `algebraic_letters` table; an algebraic-letter
token without a table entry is a typed error (upstream calls this a
wiring bug, wl:12289-12297).

Atoms may carry negative powers (single-term denominators like `1/Wp_1`).
A multi-term denominator containing atoms (the shape `back_substitute`
makes from ratios) is REFUSED with a pointer to `normalize_wmwp`
(DESIGN.md RT-v). Everything rational stays in the terms' canonical coef
strings (Nemo-normalized: content-canonical, monic denominator).
"""
function parse_coef(s::AbstractString, vars::Vector{Symbol};
                    letters::Vector{AlgLetter}=AlgLetter[])
    ctx = _pctx(vars, letters)
    # Exit fold: sqrt_disc_i^2 -> disc when the table knows
    # disc — ONE canonical spelling (the Log2-atom decision pattern); parse
    # output never carries sqrt_disc powers ≥2 or ≤-1 for disc-carrying
    # letters. Stub tables (disc="") are untouched.
    return _hv_to_hlogexpr(_hv_fold_sqrtdisc(ctx, _parse_hval(s, ctx)), ctx)
end

"""
    canonical_text(e::HlogExpr) -> String

Deterministic canonical text form (CONTRACTS.md §d): terms sorted by atom
part (fixed atom-kind order: HConst, HZeta, HLogA, HHlog, HPeriod,
HAlgLetter, HDelta), coefficients re-canonicalized through the Nemo
fraction field of `e.vars`, HF-style emission (`coef*atom^k`, coef 1
omitted, " + "/" - " joins). Byte-stable across runs; idempotent under
`parse_coef` (unit-tested).

Canonical-atom law (resolving the coordination REQUEST as option (a)): `Log[rational]` atoms are folded
into the Log2/Log-prime basis at entry (`fold_log_atoms`, below in this
file), so identically-zero combinations emit as the literal "0" from
BOTH emission chokepoints (`canonical_text` here, `mma_string` in
serialize.jl).
"""
function canonical_text(e::HlogExpr)
    e = fold_log_atoms(e)
    ctx = _pctx(e.vars, _letters_of(e))
    hv = _to_hval(e, ctx)
    ce = _hv_to_hlogexpr(hv, ctx)
    isempty(ce.terms) && return "0"
    reprs = Tuple{Bool,String}[]
    for t in ce.terms
        c = t.coef
        neg = false
        if startswith(c, "-") && !_mul_unsafe(c[2:end])
            neg = true
            c = c[2:end]
        end
        factors = String[]
        (c == "1" && !isempty(t.atoms)) || push!(factors, _paren_if_unsafe(c))
        for (a, k) in t.atoms
            s = _atom_string(a)
            push!(factors, k == 1 ? s : s * "^" * string(k))
        end
        push!(reprs, (neg, join(factors, "*")))
    end
    return _join_terms(reprs)
end

# Letters referenced by an HlogExpr (so canonical_text/_to_hval can re-parse
# its own emission without an external table).
function _letters_of(e::HlogExpr)
    idxs = Int[]
    for t in e.terms, (a, _) in t.atoms
        a isa HAlgLetter && push!(idxs, a.idx)
    end
    return [AlgLetter(i, "", :x, "", "", "", "") for i in sort(unique(idxs))]
end

# HlogExpr -> _HVal in a given context (coef strings re-parsed exactly).
function _to_hval(e::HlogExpr, ctx::_PCtx)
    hv = _hv_zero()
    for t in e.terms
        fr = _parse_rat_frac(t.coef, ctx)
        iszero(fr) && continue
        one_term = _HVal(_canon_atomvec(t.atoms) => fr)
        hv = _hv_add(hv, one_term)
    end
    return hv
end

# ---------------------------------------------------------------------------
# Exact HlogExpr arithmetic (one owner for the semantics — used by eps_expand's
# loggamma_series products and lr_dispatch's assembly sums)
# ---------------------------------------------------------------------------

function _union_vars(a::Vector{Symbol}, b::Vector{Symbol})
    out = copy(a)
    for v in b
        v in out || push!(out, v)
    end
    return out
end

function _union_letters(a::Vector{AlgLetter}, b::Vector{AlgLetter})
    out = copy(a)
    for l in b
        any(x -> x.idx == l.idx, out) || push!(out, l)
    end
    return out
end

"""
    hlog_add(a::HlogExpr, b::HlogExpr) -> HlogExpr

Exact sum. Kinematic contexts are unioned (a.vars then new vars of b);
coefficients are re-parsed into the union fraction field, merged on equal
atom parts, and re-emitted content-canonical.
"""
function hlog_add(a::HlogExpr, b::HlogExpr)
    ctx = _pctx(_union_vars(a.vars, b.vars), _union_letters(_letters_of(a), _letters_of(b)))
    return _hv_to_hlogexpr(_hv_add(_to_hval(a, ctx), _to_hval(b, ctx)), ctx)
end

"""
    hlog_mul(a::HlogExpr, b::HlogExpr) -> HlogExpr

Exact product (atom powers add; boundary ZeroInfPeriod shuffle products
collapse to ordinary products — SubTropica.wl:12926-12934). Contexts
unioned as in `hlog_add`.
"""
function hlog_mul(a::HlogExpr, b::HlogExpr)
    ctx = _pctx(_union_vars(a.vars, b.vars), _union_letters(_letters_of(a), _letters_of(b)))
    return _hv_to_hlogexpr(_hv_mul(_to_hval(a, ctx), _to_hval(b, ctx)), ctx)
end

"""
    hlog_scale(a::HlogExpr, c) -> HlogExpr

Exact scalar multiple. `c` may be an `Integer`, `Rational`, `QQFieldElem`,
or an Mma-grammar rational-function STRING over `a.vars` (parsed exactly;
atoms in `c` are refused).
"""
function hlog_scale(a::HlogExpr, c::Union{Integer,Rational,QQFieldElem})
    ctx = _pctx(a.vars, _letters_of(a))
    return _hv_to_hlogexpr(_hv_mul(_to_hval(a, ctx), _hv_const(ctx, QQ(c))), ctx)
end
function hlog_scale(a::HlogExpr, c::AbstractString)
    ctx = _pctx(a.vars, _letters_of(a))
    fr = _parse_rat_frac(c, ctx)
    hv = _HVal()
    iszero(fr) || (hv[Pair{HAtom,Int}[]] = fr)
    return _hv_to_hlogexpr(_hv_mul(_to_hval(a, ctx), hv), ctx)
end

# ---- Log[rational] → prime-log basis fold (canonical emission) ------------
#
# canonical emission folds every
# Log[rational] atom into the Log2/Log-prime basis so identically-zero
# combinations cancel TEXTUALLY (regression: the eq 4.3 ε⁻³ combination
# 3·Log2 − Log[1/2] + Log[1/4] − Log[4] must emit as the literal "0").
#
# Basis law (matches eps_expand.jl loggamma_series, the pre-existing
# prime-log convention): log 2 → HConst(:Log2) (the canonical Log2 atom
# decision, CONTRACTS.md §d); every other prime p → HLogA(string(p)).
# log(p/q) = Σ e_i·log(prime_i) with e_i ∈ Z from the exact factorization
# of numerator and denominator.
#
# Scope rules (all folds EXACT, never lossy):
#   * only HLogA atoms whose arg parses as a POSITIVE rational fold;
#     kinematic/symbolic args (Log[1+s], Log[Pi]) and non-positive
#     rationals pass through unchanged;
#   * Log[1]^k, k ≥ 1 ⇒ the whole term is exactly 0 (log 1 = 0) and is
#     dropped; Log[1]^k with k < 0 is a division by zero — LOUD refusal;
#   * single-prime args (Log[8] = 3·Log2, Log[1/3] = −Log[3]) fold for ANY
#     integer atom power via (e·Λ)^k = e^k·Λ^k;
#   * multi-prime args with NEGATIVE power (1/Log[6]) are unrepresentable
#     as HlogExpr terms — the atom is KEPT (never silently wrong);
#   * args with numerator/denominator above _FOLD_MAX_ABS are KEPT
#     (factorization cost guard; never expected from this pipeline).
#
# The fold is applied at BOTH emission chokepoints (—
# the coordination REQUEST resolved as option (a), block moved
# here from serialize.jl): `canonical_text` (entry fold, this file) and
# `mma_string(::HlogExpr)` (serialize.jl). Round-trip law holds: the
# emitted Log[2]/Log[p] forms parse back to HConst(:Log2)/HLogA("p")
# (parse_coef + _norm_atom), and the fold is idempotent on its own output.

const _FOLD_MAX_ABS = big(2)^96
const _FOLD_RAT_RE = r"^(-?[0-9]+)(?:/([0-9]+))?$"

# "p" / "p/q" (as emitted by _frac_string/_a3_ratstr/_b1_qstr) → Rational,
# else nothing. Whitespace-free canonical strings only — anything else is
# not a rational-constant atom arg.
function _fold_rat_arg(s::AbstractString)
    m = match(_FOLD_RAT_RE, s)
    m === nothing && return nothing
    num = parse(BigInt, m.captures[1])
    den = m.captures[2] === nothing ? big(1) : parse(BigInt, m.captures[2])
    iszero(den) && return nothing
    return num // den
end

# log q = Σ e·log p over primes, exact; q MUST be a positive rational.
# Returns Vector{Tuple{BigInt,Int}} sorted by prime (empty ⇔ q == 1).
function _log_prime_decomp(q::Rational{BigInt})
    q > 0 || throw(SubTropicaSerializeError(
        "_log_prime_decomp: argument $q not positive — caller must filter"))
    acc = Dict{BigInt,Int}()
    for (sgn, n) in ((1, numerator(q)), (-1, denominator(q)))
        n == 1 && continue
        for (p, e) in Nemo.factor(Nemo.ZZRingElem(n))
            pp = BigInt(p)
            acc[pp] = get(acc, pp, 0) + sgn * Int(e)
        end
    end
    out = Tuple{BigInt,Int}[(p, acc[p]) for p in sort!(collect(keys(acc)))]
    filter!(t -> t[2] != 0, out)
    return out
end

# The canonical prime-log atom (loggamma_series convention, eps_expand.jl):
# 2 → HConst(:Log2); other primes → HLogA("p").
_log_prime_atom(p::BigInt) = p == 2 ? HConst(:Log2) : HLogA(string(p))

# Foldable-atom probe: prime decomposition of a positive-rational Log arg,
# nothing when the atom must pass through unchanged.
function _foldable_log(a::HAtom)
    a isa HLogA || return nothing
    q = _fold_rat_arg(a.arg)
    q === nothing && return nothing
    q > 0 || return nothing
    (abs(numerator(q)) <= _FOLD_MAX_ABS && denominator(q) <= _FOLD_MAX_ABS) ||
        return nothing
    d = _log_prime_decomp(q)
    # already-canonical prime atom (Log[3], Log[5], ...) — passthrough, so
    # the fold is idempotent-by-identity on its own output. (HLogA("2") is
    # NOT canonical — its canonical spelling is HConst(:Log2) — and folds.)
    (length(d) == 1 && d[1][2] == 1 && a == _log_prime_atom(d[1][1])) &&
        return nothing
    return d
end

# Σ e_p·Λ_p as an HlogExpr (the fold of one first-power multi-prime atom).
_log_fold_expr(d::Vector{Tuple{BigInt,Int}}, vars::Vector{Symbol}) =
    HlogExpr(HTerm[HTerm(string(e), Pair{HAtom,Int}[_log_prime_atom(p) => 1])
                   for (p, e) in d], copy(vars))

"""
    fold_log_atoms(e::HlogExpr) -> HlogExpr

Fold every `Log[rational]` atom into the canonical Log2/Log-prime basis
(header block above): `Log[1/2] → −Log2`, `Log[4] → 2·Log2`,
`Log[6] → Log2 + Log[3]`, powers expanded multinomially through the exact
HlogExpr arithmetic (hlog_mul/hlog_add — the canonical-atom-law
chokepoints), so identically-zero combinations cancel. Returns `e`
UNCHANGED (same object) when nothing folds — byte-stable emission for
already-canonical input."""
function fold_log_atoms(e::HlogExpr)
    any(t -> any(pr -> _foldable_log(pr.first) !== nothing, t.atoms),
        e.terms) || return e
    acc = HlogExpr(HTerm[], copy(e.vars))
    for t in e.terms
        keep = Pair{HAtom,Int}[]
        factors = HlogExpr[]
        scale = Rational{BigInt}(1)
        zero_term = false
        for (a, k) in t.atoms
            d = _foldable_log(a)
            if d === nothing
                push!(keep, a => k)
            elseif isempty(d)                      # Log[1] = 0
                k > 0 && (zero_term = true; break)
                throw(SubTropicaSerializeError(
                    "fold_log_atoms: Log[1]^$k — log(1) = 0 in a denominator"))
            elseif length(d) == 1                  # (e·Λ)^k = e^k·Λ^k, any k
                (p, ex) = d[1]
                scale *= (Rational{BigInt}(ex))^k
                push!(keep, _log_prime_atom(p) => k)
            elseif k >= 1                          # multinomial expansion
                F = _log_fold_expr(d, e.vars)
                Fk = F
                for _ in 2:k
                    Fk = hlog_mul(Fk, F)
                end
                push!(factors, Fk)
            else                                   # 1/(Σ logs): unrepresentable — keep
                push!(keep, a => k)
            end
        end
        zero_term && continue
        piece = HlogExpr(HTerm[HTerm(t.coef, _canon_atomvec(keep))], copy(e.vars))
        isone(scale) || (piece = hlog_scale(piece, scale))
        for F in factors
            piece = hlog_mul(piece, F)
        end
        acc = hlog_add(acc, piece)
    end
    return acc
end


# ---------------------------------------------------------------------------
# Wm/Wp normalization through HF's OWN ops (DESIGN.md RT-v — no local Vieta)
# ---------------------------------------------------------------------------

"""
    _ratparse_safe(s) -> Bool

True iff `s` is in the transport class the engine's `Rat::parse` reads
with Mma-compatible semantics (live-probed, CONTRACTS.md §d):
either NO '/' at all (a polynomial body), or EXACTLY ONE '/' at bracket
depth 0 whose denominator group extends to end-of-string, with no
depth-0 '+'/'-' elsewhere (leading unary minus allowed). Everything else
is refused for whole-string transport: the engine in-band-errors on its
own parenthesized "(N/D)" coef emission and SILENTLY mis-parses trailing
factors after a "/(...)" group into the denominator
("N/(D)*mzv_2" -> N/(D*mzv_2) — live probe).
"""
function _ratparse_safe(s::AbstractString)
    depth = 0
    slash_at = 0
    nslash = 0
    plusminus = false
    cs = collect(s)
    for (i, c) in enumerate(cs)
        if c in ('(', '[', '{')
            depth += 1
        elseif c in (')', ']', '}')
            depth -= 1
        elseif c == '/'
            (depth == 0 && nslash == 0) || return false
            nslash = 1
            slash_at = i
        elseif depth == 0 && (c == '+' || c == '-') && i > 1
            plusminus = true
        end
    end
    nslash == 0 && return true       # polynomial body: sums are fine
    plusminus && return false        # sums + division: precedence not trusted
    # single depth-0 '/': everything after it must be ONE factor reaching
    # end-of-string (no depth-0 operators after the denominator group).
    depth = 0
    for i in (slash_at + 1):length(cs)
        c = cs[i]
        if c in ('(', '[', '{')
            depth += 1
        elseif c in (')', ']', '}')
            depth -= 1
        elseif depth == 0 && c in ('+', '-', '*', '/', '^')
            return false
        end
    end
    return true
end

# One pinned-order chain pass through HF's OWN ops on a Rat::parse-safe
# body. `allocs`/`vv` are the positional letter replay + PolyCtx vars.
function _wmwp_chain(s::AbstractString, allocs, vv, ops, cfg)
    r = String(s)
    for op in ops
        resp = hf_call(String(op), Dict{String,Any}("expr" => r,
            "allocations" => allocs, "vars" => vv); cfg=cfg)
        r = String(resp["result"])
    end
    return r
end

# Normalize the Wm/Wp SUB-MONOMIAL of every term through the chain (the
# numerator-side leg: Vieta products like Wm_1*Wp_1 and lone roots reduce
# engine-side; sqrt_disc powers fold locally afterwards). Chain results
# that fail to parse leave the original — exactness-preserving in both
# directions (same value either way); cached per distinct sub-monomial.
function _normalize_terms_via_chain(ctx::_PCtx, hv::_HVal, chain)
    cache = Dict{Vector{Pair{HAtom,Int}},Union{Nothing,_HVal}}()
    out = _hv_zero()
    rctx = _pctx_norouter(ctx)
    for (av, c) in hv
        wm = Pair{HAtom,Int}[p for p in av
             if p.first isa HAlgLetter && p.first.kind in (:Wm, :Wp)]
        if isempty(wm)
            out = _hv_add(out, _HVal(av => c))
            continue
        end
        rest = Pair{HAtom,Int}[p for p in av
               if !(p.first isa HAlgLetter && p.first.kind in (:Wm, :Wp))]
        nv = get!(cache, wm) do
            pos = [(a, k) for (a, k) in wm if k > 0]
            neg = [(a, -k) for (a, k) in wm if k < 0]
            f(t) = t[2] == 1 ? _atom_string(t[1]) :
                               _atom_string(t[1]) * "^" * string(t[2])
            S = isempty(pos) ? "1" : join(map(f, pos), "*")
            isempty(neg) || (S *= "/(" * join(map(f, neg), "*") * ")")
            routed = try
                chain(S)
            catch e
                e isa HFException || rethrow()
                return nothing            # engine refused: keep original
            end
            try
                _hv_fold_sqrtdisc(ctx, _parse_hval(String(routed), rctx))
            catch e
                e isa SubTropicaParseError || rethrow()
                nothing                   # unparseable normal form: keep original
            end
        end
        if nv === nothing
            out = _hv_add(out, _HVal(av => c))
        else
            term = _hv_mul(_HVal(_canon_atomvec(rest) => c), nv)
            out = _hv_add(out, term)
        end
    end
    return out
end

"""
    normalize_wmwp(coef::AbstractString, letters::Vector{AlgLetter};
                   vars::Vector{Symbol}=Symbol[], cfg::HFConfig=HFConfig(),
                   ops=("simplify_with_vieta","combine_wm_wp_ratios","back_substitute"))
        -> (e::HlogExpr, s::String)

Normalize a coef string over Wm_i/Wp_i algebraic letters by routing it
through HyperFLINT's OWN rewrite ops (DESIGN.md RT-v: never reimplement
the Vieta/Sqrt[disc] algebra of SubTropica.wl:12930-13005 locally). Each
op replays the AlgebraicLetterTable via the `allocations` field
(reference/hyperflint_cli_main.cpp:2634-2812; the table is PER-PROCESS and
cleared per invocation — CONTRACTS.md §a), so `letters` must be complete
with idx == position (1-based, validated here). The PolyCtx `vars` sent
engine-side always include each letter's OWN variable (root-cause fix: the allocation replay hard-fails — stderr + empty stdout,
an HFTransportError — when the letter var is missing; the eq 4.11 k0
reroute failure).

Default op order is live-pinned (see file header):
simplify_with_vieta (Vieta products/sums), THEN combine_wm_wp_ratios
(literal Wm_i/Wp_i ratios -> WmOverWp_i — must run before back_substitute,
which smears ratios into sqrt-disc-denominator shapes neither it nor
combine can recover), THEN back_substitute (residual Wm-Wp differences ->
sqrt_disc_i).

Transport (CONTRACTS.md §d diagnosis): the WHOLE
string goes through the chain only when `_ratparse_safe` — the engine
cannot read its own parenthesized "(N/D)" coef emission and silently
mis-parses "N/(D)*t". Otherwise the string is parsed locally with the
chain ARMED as the denominator router (`_PCtx.wmwp_router`: each
irreducible multi-term Wm/Wp denominator is routed as the single-division
body "1/(D)"), followed by the per-term Wm/Wp sub-monomial chain pass and
the sqrt_disc^2 -> disc table fold. Refusal is typed ONLY if a routed
NORMALIZED form still fails to parse. Returns the parsed HlogExpr and
the normalized string (engine string on the whole-string leg, canonical
text on the structural leg).
"""
function normalize_wmwp(coef::AbstractString, letters::Vector{AlgLetter};
                        vars::Vector{Symbol}=Symbol[], cfg::HFConfig=HFConfig(),
                        ops=("simplify_with_vieta", "combine_wm_wp_ratios",
                             "back_substitute"))
    ls = sort(letters, by=l -> l.idx)
    for (i, l) in enumerate(ls)
        l.idx == i || throw(SubTropicaEvalError("AlgLetter table gap: idx $(l.idx) " *
            "at position $i — the allocations replay is positional/1-based " *
            "(hyperflint_cli_main.cpp:2755-2760), Wm_<i> references would break"))
    end
    allocs = [Dict{String,Any}("polynomial" => l.poly, "var" => string(l.var)) for l in ls]
    vv = String[string(v) for v in vars]
    for l in ls                        # ROOT-CAUSE FIX (k0)
        string(l.var) in vv || push!(vv, string(l.var))
    end
    isempty(vv) && push!(vv, "x")      # handler defaults "x" (main.cpp:2647)
    chain = s -> _wmwp_chain(s, allocs, vv, ops, cfg)
    s0 = String(coef)
    # --- leg 1: whole-string transport (safe shapes only) -------------------
    if _ratparse_safe(s0)
        s1 = try
            chain(s0)
        catch e
            e isa HFException || rethrow()
            nothing                    # in-band engine refusal -> structural leg
        end
        if s1 !== nothing
            e1 = try
                ctx1 = _pctx(vars, ls)
                _hv_to_hlogexpr(_hv_fold_sqrtdisc(ctx1, _parse_hval(s1, ctx1)), ctx1)
            catch e
                e isa SubTropicaParseError || rethrow()
                nothing                # normal form not representable -> leg 2
            end
            e1 === nothing || return (e1, s1)
        end
    end
    # --- leg 2: structural routing (parser-armed chain) ---------------------
    ctx = _pctx(vars, ls; wmwp_router=chain)
    hv = _parse_hval(s0, ctx)
    hv = _normalize_terms_via_chain(ctx, hv, chain)
    hv = _hv_fold_sqrtdisc(ctx, hv)
    e = _hv_to_hlogexpr(hv, ctx)
    return (e, canonical_text(e))
end

# ---------------------------------------------------------------------------
# Numeric verification hooks (unit gate A(v) / RT-extra-period): ginac_gpl
# is the INDEPENDENT evaluator — this is our code with no upstream to diff.
# ---------------------------------------------------------------------------

const _GINAC_GPL_DEFAULT =
    joinpath(@__DIR__, "..", "..", "gpl-eval", "GPLEval.jl", "src", "ginac_gpl")  # sibling vendored package; override via SUBTROPICA_GINAC_GPL
_ginac_bin() = get(ENV, "SUBTROPICA_GINAC_GPL", _GINAC_GPL_DEFAULT)

const _RAT_LITERAL_RE = r"^([+-]?\d+)(?:/([+-]?\d+))?$"

function _parse_rat_literal(s::AbstractString)
    m = match(_RAT_LITERAL_RE, strip(s))
    m === nothing && throw(SubTropicaEvalError("letter '$s' is not a rational literal"))
    num = parse(BigInt, m.captures[1])
    den = m.captures[2] === nothing ? BigInt(1) : parse(BigInt, m.captures[2])
    iszero(den) && throw(SubTropicaEvalError("letter '$s' has zero denominator"))
    return num // den
end

_rat_arg(q::Rational{BigInt}) =
    denominator(q) == 1 ? string(numerator(q)) :
                          string(numerator(q)) * "/" * string(denominator(q))

"""
    ginac_gpl_G(letters::Vector{Rational{BigInt}}, z::AbstractString;
                digits::Int=70, timeout_s::Real=600) -> BigFloat

G(letters...; z) via the independent ginac_gpl CLI
(tools/gpl-eval/GPLEval.jl/src/ginac_gpl; override SUBTROPICA_GINAC_GPL).
`z` is passed verbatim (use exact integer/rational strings). Real-line
use only: throws if the imaginary part is not negligible. The result
BigFloat is parsed at ≥(digits+guard) working precision (known pitfall:
parsing at the ambient precision loses digits).
"""
function ginac_gpl_G(letters::Vector{Rational{BigInt}}, z::AbstractString;
                     digits::Int=70, timeout_s::Real=600)
    bin = _ginac_bin()
    isfile(bin) || throw(SubTropicaEvalError("ginac_gpl CLI not found at $bin " *
        "(set SUBTROPICA_GINAC_GPL)"))
    args = String[string(digits), String(z), "0"]
    for l in letters
        push!(args, _rat_arg(l), "0")
    end
    r = _spawn_capture(vcat(["bash", "-c",
            "ulimit -v 32505856; exec \"\$0\" \"\$@\"", bin], args), nothing;
            timeout_s=timeout_s)
    r.timed_out && throw(SubTropicaEvalError("ginac_gpl timed out after $(timeout_s)s (pid $(r.pid), killed by PID)"))
    lines = split(strip(r.stdout), '\n')
    length(lines) >= 2 || throw(SubTropicaEvalError("ginac_gpl bad output (rc=$(r.rc)): " *
        first(split(r.stderr * " ", '\n'))))
    prec = ceil(Int, digits * 3.33) + 64
    return setprecision(BigFloat, prec) do
        re = parse(BigFloat, replace(String(lines[1]), "E" => "e"))
        im = parse(BigFloat, replace(String(lines[2]), "E" => "e"))
        abs(im) <= max(abs(re), one(BigFloat)) * big(10.0)^(-digits + 8) ||
            throw(SubTropicaEvalError("ginac_gpl returned non-negligible imaginary part $im"))
        re
    end
end

"""
    zero_inf_period_ginac(word::Vector{String}; digits::Int=70) -> BigFloat

INDEPENDENT numeric value of ZeroInfPeriod[word] (rational letters only):
the shuffle-regularized constant term (log z -> 0) of G(word; z) as
z -> +∞ — semantics pinned by live probe (file header).
Method: evaluate G at n+1 nodes z_j = 10^(120+8j) via ginac_gpl and
Lagrange-extrapolate the degree-≤n polynomial in log z to log z = 0; the
O(1/z·log^n z) truncation is ≲10^-110 and the extrapolation loses ~n·log10(15)
digits, both covered by the +40-digit guard (prototype:
113-118 matched digits on 4 control words at digits=70). Used by the
test-corpus unit gate (every residual key cross-evaluated ≥60 dps,
DESIGN.md gate A(v)).
"""
function zero_inf_period_ginac(word::Vector{String}; digits::Int=70)
    isempty(word) && throw(SubTropicaEvalError("empty ZeroInfPeriod word"))
    letters = [_parse_rat_literal(w) for w in word]
    n = length(word)
    guard = 40
    prec = ceil(Int, (digits + guard) * 3.33) + 64
    return setprecision(BigFloat, prec) do
        E0, S = 120, 8
        Ls = BigFloat[]
        gs = BigFloat[]
        for j in 0:n
            push!(Ls, (E0 + S * j) * log(BigFloat(10)))
            push!(gs, ginac_gpl_G(letters, "1" * "0"^(E0 + S * j);
                                  digits=digits + guard))
        end
        c0 = zero(BigFloat)
        for j in 0:n
            w = one(BigFloat)
            for k in 0:n
                k == j && continue
                w *= (zero(BigFloat) - Ls[k+1]) / (Ls[j+1] - Ls[k+1])
            end
            c0 += w * gs[j+1]
        end
        c0
    end
end

# ---------------------------------------------------------------------------
# Word-level algebraic-letter boundary periods (closes
# the DEEPER layer of the Wm/Wp gap recorded in the polish-23 status note: the
# letters tier can mint Wm_i/Wp_i tokens INSIDE ZeroInfPeriod words, e.g.
# eq 4.11's ["-1","0","Wm_1"] over the Gaussian quadratic x2^2+1. These
# periods are algebraic-NUMBER-valued (the roots are constants), generally
# COMPLEX per word; conjugate Wm/Wp words cancel imaginary parts only in
# the assembled coefficient sum. There is no exact atom form for a complex-
# letter G here, so the SYMBOLIC result keeps the explicit HPeriod atom
# (resolve_periods strict=false path) and the EVALUATION leg is numeric:
# the same large-z Lagrange extrapolation as zero_inf_period_ginac at
# complex letters (the ginac_gpl CLI takes re/im argument pairs), with the
# roots resolved from the response's AlgLetter table. eval_symbolic
# (verify.jl) owns the numeric leg and the per-coefficient imaginary-part
# cancellation assertion. Promoted from the test-local evaluator in
# test_wmwp_integration.jl (which remains as the acceptance harness).
# Truncation floor: E0=120 ⇒ O(1/z·log^n z) ≈ 1e-110 — reliable to ~100
# digits, far above every ≥30d gate bar (same envelope as
# zero_inf_period_ginac's prototype pins).
# ---------------------------------------------------------------------------

# Complex-letter ginac_gpl call: letters as (re, im) DECIMAL-string pairs.
function _ginac_gpl_G_c(lets::Vector{Tuple{String,String}}, z::AbstractString;
                        digits::Int=70, timeout_s::Real=600)
    bin = _ginac_bin()
    isfile(bin) || throw(SubTropicaEvalError("ginac_gpl CLI not found at $bin " *
        "(set SUBTROPICA_GINAC_GPL)"))
    args = String[string(digits), String(z), "0"]
    for (re, im) in lets
        push!(args, re, im)
    end
    r = _spawn_capture(vcat(["bash", "-c",
            "ulimit -v 32505856; exec \"\$0\" \"\$@\"", bin], args), nothing;
            timeout_s=timeout_s)
    r.timed_out && throw(SubTropicaEvalError("ginac_gpl timed out after " *
        "$(timeout_s)s (pid $(r.pid), killed by PID)"))
    lines = split(strip(r.stdout), '\n')
    length(lines) >= 2 || throw(SubTropicaEvalError("ginac_gpl bad output " *
        "(rc=$(r.rc)): " * first(split(r.stderr * " ", '\n'))))
    prec = ceil(Int, digits * 3.33) + 64
    return setprecision(BigFloat, prec) do
        Complex{BigFloat}(parse(BigFloat, replace(String(lines[1]), "E" => "e")),
                          parse(BigFloat, replace(String(lines[2]), "E" => "e")))
    end
end

"""
    alg_letter_roots(l::AlgLetter; digits=70)
        -> (Wm, Wp, sqrt_disc)::NTuple{3,Complex{BigFloat}}

Numeric root values of one AlgLetter table entry under the MONIC law
(back_substitute pin `Wm − Wp → −sqrt_disc`; acceptance harness
test_wmwp_integration.jl): `Wm = (sum − sqrt(disc))/2`,
`Wp = (sum + sqrt(disc))/2`, both ROOT-CHECKED against
`x² − sum·x + product` to the working precision. Typed refusal on
non-monic letters (lc ≠ 1), stub tables (empty fields — pass the RESPONSE
table / LaurentSeries.letters, never `_letters_of` stubs), and
non-rational (kinematic) sum/product/disc.
"""
function alg_letter_roots(l::AlgLetter; digits::Int=70)
    for (fname, s) in (("sum", l.sum), ("product", l.product), ("disc", l.disc))
        isempty(strip(s)) && throw(SubTropicaEvalError("AlgLetter idx $(l.idx) " *
            "has an empty '$fname' field (stub table from _letters_of?) — " *
            "pass the response/LaurentSeries.letters table"))
    end
    strip(l.lc) == "1" || throw(SubTropicaEvalError("AlgLetter idx $(l.idx) is " *
        "non-monic (lc=$(l.lc)) — the monic root law Wm/Wp=(sum∓√disc)/2 " *
        "does not apply; extend the law first (never guess a normalization)"))
    s = _parse_rat_literal(l.sum)
    p = _parse_rat_literal(l.product)
    D = _parse_rat_literal(l.disc)
    prec = ceil(Int, (digits + 20) * 3.33) + 64
    return setprecision(BigFloat, prec) do
        bf = q -> BigFloat(numerator(q)) / BigFloat(denominator(q))
        sd = sqrt(Complex{BigFloat}(bf(D)))
        Wm = (bf(s) - sd) / 2
        Wp = (bf(s) + sd) / 2
        tol = big(10.0)^(-digits)
        for W in (Wm, Wp)
            abs(W^2 - bf(s) * W + bf(p)) <= tol ||
                throw(SubTropicaEvalError("AlgLetter idx $(l.idx) ROOT CHECK " *
                    "FAILED: table sum/product/disc mutually inconsistent " *
                    "(|W² − sum·W + product| > 1e-$(digits))"))
        end
        (Wm, Wp, sd)
    end
end

const _ZIP_ALG_CACHE = Dict{Tuple{Vector{Tuple{String,String}},Int},Complex{BigFloat}}()
const _WORD_ROOT_TOKEN_RE = r"^(Wm|Wp)_(\d+)$"

"""
    zero_inf_period_alg(word::Vector{String}, letters::Vector{AlgLetter};
                        digits::Int=70) -> Complex{BigFloat}

NUMERIC value of a ZeroInfPeriod word whose letters are rational literals
and/or LONE `Wm_i`/`Wp_i` root tokens (resolved via `alg_letter_roots`
from the AlgLetter table). Same shuffle-regularized large-z Lagrange
extrapolation as `zero_inf_period_ginac`, at complex letters. Generally
COMPLEX — the caller (eval_symbolic) owns the imaginary-part cancellation
assertion across the full coefficient. Typed refusal on any other word
letter shape (compound algebraic-letter expressions, `WmOverWp`/
`sqrt_disc` in a WORD = contract drift, kinematic letters). Cached per
(resolved letters, digits); conjugate symmetry (real base point ⇒
conjugate letters give the conjugate value) halves the engine calls.
"""
function zero_inf_period_alg(word::Vector{String}, letters::Vector{AlgLetter};
                             digits::Int=70)
    isempty(word) && throw(SubTropicaEvalError("empty ZeroInfPeriod word"))
    n = length(word)
    guard = 40
    prec = ceil(Int, (digits + guard) * 3.33) + 64
    lets = setprecision(BigFloat, prec) do
        out = Tuple{String,String}[]
        for w0 in word
            w = strip(w0)
            m = match(_WORD_ROOT_TOKEN_RE, w)
            if m === nothing
                occursin(_ALGLETTER_TOKEN_RE, w) && throw(SubTropicaEvalError(
                    "ZeroInfPeriod word letter '$w': only LONE Wm_i/Wp_i " *
                    "root tokens are supported in words (WmOverWp/sqrt_disc " *
                    "or compound letters in a WORD is contract drift)"))
                q = _parse_rat_literal(w)    # typed refusal on kinematics
                push!(out, (_rat_arg(q), "0"))
            else
                idx = parse(Int, m.captures[2])
                li = findfirst(l -> l.idx == idx, letters)
                li === nothing && throw(SubTropicaEvalError("ZeroInfPeriod " *
                    "word letter '$w': no AlgLetter idx $idx in the table " *
                    "— pass the response/LaurentSeries.letters table"))
                (Wm, Wp, _) = alg_letter_roots(letters[li]; digits=digits + guard)
                W = m.captures[1] == "Wm" ? Wm : Wp
                push!(out, (string(real(W)), string(imag(W))))
            end
        end
        out
    end
    all(l -> l == ("0", "0"), lets) && return Complex{BigFloat}(0)
        # all-zero word: (log z)^k/k! → 0 under the log z → 0 regularization
    key = (lets, digits)
    haskey(_ZIP_ALG_CACHE, key) && return _ZIP_ALG_CACHE[key]
    negim = Tuple{String,String}[(re, startswith(im, "-") ? String(im[2:end]) :
                 (im == "0" ? im : "-" * im)) for (re, im) in lets]
    haskey(_ZIP_ALG_CACHE, (negim, digits)) &&
        return (_ZIP_ALG_CACHE[key] = conj(_ZIP_ALG_CACHE[(negim, digits)]))
    val = setprecision(BigFloat, prec) do
        E0, S = 120, 8
        Ls = BigFloat[]
        gs = Complex{BigFloat}[]
        for j in 0:n
            push!(Ls, (E0 + S * j) * log(BigFloat(10)))
            push!(gs, _ginac_gpl_G_c(lets, "1" * "0"^(E0 + S * j);
                                     digits=digits + guard))
        end
        c0 = Complex{BigFloat}(0)
        for j in 0:n
            w = one(BigFloat)
            for k in 0:n
                k == j && continue
                w *= (zero(BigFloat) - Ls[k+1]) / (Ls[j+1] - Ls[k+1])
            end
            c0 += w * gs[j+1]
        end
        c0
    end
    _ZIP_ALG_CACHE[key] = val
    return val
end

"""
    mzv_bigfloat(idx::Vector{Int}; digits::Int=70) -> BigFloat

Numeric value of the HF `mzv_*` atom ζ(s_1,...,s_k) — HyperIntica/HyperInt
ASCENDING convention Σ_{0<n_1<...<n_k} ∏ σ_i^{n_i}/n_i^{|s_i|}, σ_i = -1
for 'm'-prefixed (negative) entries — via the independent ginac_gpl CLI:
ζ = (-1)^k · G(w; 1) with w built outermost-first
(w = ⨁_{t=k..1} [0^{|s_t|-1}, ∏_{j=t}^k σ_j]; standard G-representation of
multiple polylogs at ±1). The convention is not trusted on faith: it is
GATED numerically ≥60 dps against HF's own reductions in
test_hlogexpr.jl's period corpus. Divergent (s_k=+1 leading-letter-1)
input throws.
"""
function mzv_bigfloat(idx::Vector{Int}; digits::Int=70)
    isempty(idx) && throw(SubTropicaEvalError("empty mzv index"))
    k = length(idx)
    (abs(idx[end]) == 1 && idx[end] > 0) &&
        throw(SubTropicaEvalError("mzv_$(join(idx,'_')) diverges (outermost index +1)"))
    key = (copy(idx), digits)
    haskey(_MZV_CACHE, key) && return _MZV_CACHE[key]
    word = Rational{BigInt}[]
    for t in k:-1:1
        for _ in 1:(abs(idx[t]) - 1)
            push!(word, big(0) // 1)
        end
        sgn = prod(j -> idx[j] < 0 ? -1 : 1, t:k)
        push!(word, big(sgn) // 1)
    end
    v = ginac_gpl_G(word, "1"; digits=digits)
    v = isodd(k) ? -v : v
    _MZV_CACHE[key] = v
    return v
end
const _MZV_CACHE = Dict{Tuple{Vector{Int},Int},BigFloat}()

"""
    hlog_const_bigfloat(e::HlogExpr; digits::Int=70) -> BigFloat

Numeric value of a CONSTANT HlogExpr (real-valued): coefficients must be
rational constants; atoms HConst(:Pi/:EulerGamma/:Catalan/:Log2), HZeta
(via `mzv_bigfloat`/ginac), HLogA with rational argument. Typed refusal on
anything non-constant or non-real (kinematic coefs, HConst(:I), HHlog,
HPeriod, HAlgLetter, HDelta) — this is the verification-side evaluator for
the period unit gate, not a general evaluator.
"""
function hlog_const_bigfloat(e::HlogExpr; digits::Int=70)
    prec = ceil(Int, (digits + 20) * 3.33) + 64
    return setprecision(BigFloat, prec) do
        ctx = _pctx(e.vars, _letters_of(e))
        tot = zero(BigFloat)
        for t in e.terms
            fr = _parse_rat_frac(t.coef, ctx)
            num = numerator(fr); den = denominator(fr)
            (is_constant(num) && is_constant(den)) ||
                throw(SubTropicaEvalError("non-constant coefficient '$(t.coef)' — " *
                    "hlog_const_bigfloat evaluates constants only"))
            q = constant_coefficient(num) // constant_coefficient(den)
            v = BigFloat(numerator(q)) / BigFloat(denominator(q))
            for (a, p) in t.atoms
                v *= _atom_bigfloat(a; digits=digits)^p
            end
            tot += v
        end
        tot
    end
end

function _atom_bigfloat(a::HAtom; digits::Int=70)
    if a isa HConst
        a.name === :Pi && return BigFloat(pi)
        a.name === :EulerGamma && return BigFloat(Base.MathConstants.eulergamma)
        a.name === :Catalan && return BigFloat(Base.MathConstants.catalan)
        a.name === :Log2 && return log(BigFloat(2))
        throw(SubTropicaEvalError("HConst(:$(a.name)) has no real numeric value here " *
            "(I is refused: real-line evaluator)"))
    elseif a isa HZeta
        return mzv_bigfloat(a.idx; digits=digits)
    elseif a isa HLogA
        q = _parse_rat_literal(a.arg)
        q > 0 || throw(SubTropicaEvalError("Log[$(a.arg)] non-positive argument (real-line evaluator)"))
        return log(BigFloat(numerator(q)) / BigFloat(denominator(q)))
    end
    throw(SubTropicaEvalError("$(typeof(a)) atom is not a numeric constant " *
        "(resolve periods/letters first)"))
end

"""
    _matched_digits(a, b) -> Int

Floored count of matched decimal digits between two reals:
floor(-log10(|a-b| / max(|a|,|b|))), computed at (max input precision)+32
bits — NEVER at the ambient BigFloat precision, which would round both
operands together and report false equality (mplll ambient-dps footgun
class; caught live). If either side is exactly zero the
ABSOLUTE error is counted (relative digits are undefined against 0).
999 for exact equality. Internal helper for the period unit gate —
verify.jl owns the production-facing comparator.
"""
function _matched_digits(a::Real, b::Real)
    pa = a isa BigFloat ? precision(a) : 256
    pb = b isa BigFloat ? precision(b) : 256
    return setprecision(BigFloat, max(pa, pb) + 32) do
        x = BigFloat(a); y = BigFloat(b)
        x == y && return 999
        d = abs(x - y)
        if iszero(x) || iszero(y)
            return floor(Int, Float64(-log10(d)))   # absolute digits vs exact 0
        end
        return floor(Int, Float64(-log10(d / max(abs(x), abs(y)))))
    end
end
