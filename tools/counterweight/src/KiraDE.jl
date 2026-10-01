# KiraDE.jl — build the connection matrix A(eps,x) from Kira derive_dgl output.
#
# Kira's `derive_dgl` job writes, per master M_a and per kinematic variable v ∈ {s,t,...},
#     d/dv M_a = sum_b  c_{ab}(d,s,t)  M_b
# in the same text format as its reduction tables.  We parse those rows, substitute d = 4 - 2*eps,
# set the overall scale (pull out (-s) so the system is in the ratio x = t/s), and assemble the
# matrix A over the field tower Q(eps)(x).
#
# Two entry points:
#   load_kira_de(path, masters, var; subst) — parse a dgl file into Dict{(a,b)=>coeff string}.
#   build_A_from_rows(ctx, rows, masters; ...) — assemble the Qex matrix for the ratio variable.
#
# Because the engine's distinguished variable is the ratio x, the caller provides the scaling
# rule mapping (s,t) -> (overall, x).  For the standard 2-scale massless case s,t with x=t/s the
# s-derivative connection A_s and t-derivative A_t combine: in x the relevant single-variable
# connection is  A_x = s * A_t |_{t = x s}  (since dx = dt/s at fixed s), with the overall (-s)^a
# weight handled by a diagonal shift.  KiraDE exposes the raw coefficient strings and a sympy-
# free substitution helper; the heavy lifting (d->eps, t->x s, clearing s) is done with a small
# string->Nemo parser restricted to the rational functions Kira emits.

module KiraDE

using Nemo
using ..RatFunc

export load_kira_de, build_A_from_rows, parse_kira_coeff

"""
    load_kira_de(path; family="A1") -> Vector of (lhs::String, rhs::String, coeff::String)

Parse a Kira dgls file.  Format: blocks separated by blank lines; first line of a block is the
derivative-target integral, subsequent lines are `FAM[...] * (coeff)`.
"""
function load_kira_de(path::String; family::String="A1")
    txt = read(path, String)
    blocks = split(txt, r"\n\s*\n")
    rows = Tuple{String,String,String}[]
    for b in blocks
        isempty(strip(b)) && continue
        lines = [l for l in split(b, "\n") if !isempty(strip(l))]
        m = match(Regex("\\s*($(family)\\[[^\\]]*\\])"), lines[1])
        m === nothing && continue
        lhs = m.captures[1]
        for l in lines[2:end]
            mm = match(Regex("\\s*($(family)\\[[^\\]]*\\])\\s*\\*\\s*\\((.*)\\)\\s*\$"), l)
            mm === nothing && continue
            push!(rows, (lhs, mm.captures[1], mm.captures[2]))
        end
    end
    return rows
end

"""
    parse_kira_coeff(ctx, s_str; d_to_eps=true, svar, tvar) -> Qex element

Parse a Kira coefficient string (a rational function in d, s, t) into the engine field.
We use a minimal recursive-descent over Julia's own parser: Kira emits standard
arithmetic with `d`, `s`, `t`, integers, `+ - * / ^ ( )`.  We evaluate it symbolically by
substituting Nemo generators.  d -> 4 - 2*eps;  the scale is handled by the caller via
the substitution closure `kin` mapping the symbols to Qex values.
"""
function parse_kira_coeff(ctx::Ctx, s_str::AbstractString, kin::Dict{Symbol,<:Any})
    expr = Meta.parse(String(s_str))
    return _evalexpr(ctx, expr, kin)
end

function _evalexpr(ctx::Ctx, e, kin)
    if e isa Integer
        return ctx.Qex(e)
    elseif e isa Symbol
        haskey(kin, e) || error("unknown symbol in Kira coeff: $e")
        return kin[e]
    elseif e isa Expr
        if e.head == :macrocall
            # Big-integer literals: @int128_str / @big_str "NNN" — Julia parses
            # integers > Int64 as macrocall.  Evaluate the literal and lift to Qex.
            v = Core.eval(Main, e)
            return ctx.Qex(ZZ(v))
        elseif e.head == :call
            op = e.args[1]
            args = [_evalexpr(ctx, a, kin) for a in e.args[2:end]]
            if op == :+
                return reduce(+, args)
            elseif op == :-
                return length(args)==1 ? -args[1] : reduce(-, args)
            elseif op == :*
                return reduce(*, args)
            elseif op == :/ || op == ://
                # `//` appears when our own Nemo-printed matrices are fed back (round-trip
                # certify of an Anew produced by the epsform action)
                return args[1] // args[2]
            elseif op == :^
                base = args[1]
                p = e.args[3]
                return base ^ Int(p)
            else
                error("unsupported op $op")
            end
        else
            error("unsupported expr head $(e.head)")
        end
    elseif e isa AbstractFloat
        error("unexpected float in Kira coeff (expected exact rational)")
    else
        error("unsupported coeff token $(e)")
    end
end

"""
    build_A_from_rows(ctx, rows, masters; kin) -> Qex matrix

Assemble the n×n connection matrix for the distinguished variable, where `masters` is the
ordered Vector{String} of master labels (FAM[...] keys), `rows` the parsed (lhs,rhs,coeff)
list, and `kin` the symbol->Qex substitution (must define :d, :s, :t and whatever the ratio
scaling needs).  Entry A[a,b] = sum of parsed coeffs of rows with lhs=masters[a], rhs=masters[b].
"""
function build_A_from_rows(ctx::Ctx, rows, masters::Vector{String}, kin::Dict{Symbol,<:Any})
    n = length(masters)
    idx = Dict(m=>i for (i,m) in enumerate(masters))
    A = zero_matrix(ctx.Qex, n, n)
    for (lhs, rhs, c) in rows
        haskey(idx, lhs) || continue
        haskey(idx, rhs) || continue
        a = idx[lhs]; b = idx[rhs]
        A[a,b] += parse_kira_coeff(ctx, c, kin)
    end
    return A
end

end # module
