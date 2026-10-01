# Parser for the subset of Wolfram-Language InputForm that appears in SOFIA's
# inputs and outputs: PLD_database files, singularity lists, mass labels.
#
# Grammar covered: integers, rationals (a/b), symbols, indexed symbols m[1],
# function heads (Sqrt[...]), + - * / ^, parentheses, lists { ... },
# assignments `name = expr;` (returned as pairs).
#
# The AST is deliberately tiny: Int/Rational, Symbol, or WLCall(head, args).
# `tonemo` maps an AST to a Nemo multivariate polynomial given a variable map.

struct WLCall
    head::Symbol
    args::Vector{Any}
end

Base.:(==)(a::WLCall, b::WLCall) = a.head == b.head && a.args == b.args
Base.hash(a::WLCall, h::UInt) = hash(a.args, hash(a.head, hash(:WLCall, h)))

# ---------------- tokenizer ----------------

struct WLToken
    kind::Symbol   # :num, :sym, :op
    val::Any
end

function wltokenize(s::AbstractString)
    toks = WLToken[]
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
            push!(toks, WLToken(:num, parse(BigInt, s[i:prevind(s, j)])))
            i = j
        elseif isletter(c)
            j = i
            while j <= n && (isletter(s[j]) || isdigit(s[j]) || s[j] == '$')
                j = nextind(s, j)
            end
            push!(toks, WLToken(:sym, Symbol(s[i:prevind(s, j)])))
            i = j
        elseif c in ('+', '-', '*', '/', '^', '(', ')', '{', '}', '[', ']', ',', '=', ';')
            push!(toks, WLToken(:op, Symbol(c)))
            i = nextind(s, i)
        elseif c == '(' # unreachable; kept for clarity
            i = nextind(s, i)
        else
            error("wlparser: unexpected character $(repr(c)) at index $i")
        end
    end
    return toks
end

# ---------------- Pratt parser ----------------

mutable struct WLParser
    toks::Vector{WLToken}
    pos::Int
end

peek(p::WLParser) = p.pos <= length(p.toks) ? p.toks[p.pos] : nothing
advance!(p::WLParser) = (t = p.toks[p.pos]; p.pos += 1; t)

function expect!(p::WLParser, op::Symbol)
    t = advance!(p)
    (t.kind == :op && t.val == op) || error("wlparser: expected $op, got $(t.val)")
    return t
end

isop(t, op) = t !== nothing && t.kind == :op && t.val == op

# binding powers: + - (10), * / (20), unary - (25), ^ (30, right-assoc)
function parse_expr(p::WLParser, minbp::Int=0)
    t = advance!(p)
    local lhs
    if t.kind == :num
        lhs = t.val
    elseif t.kind == :sym
        lhs = t.val
        # indexed symbol or function call: sym[...]
        if isop(peek(p), Symbol("["))
            advance!(p)
            args = Any[]
            if !isop(peek(p), Symbol("]"))
                push!(args, parse_expr(p))
                while isop(peek(p), Symbol(","))
                    advance!(p)
                    push!(args, parse_expr(p))
                end
            end
            expect!(p, Symbol("]"))
            lhs = WLCall(t.val, args)
        end
    elseif t.kind == :op && t.val == Symbol("(")
        lhs = parse_expr(p)
        expect!(p, Symbol(")"))
    elseif t.kind == :op && t.val == Symbol("{")
        args = Any[]
        if !isop(peek(p), Symbol("}"))
            push!(args, parse_expr(p))
            while isop(peek(p), Symbol(","))
                advance!(p)
                push!(args, parse_expr(p))
            end
        end
        expect!(p, Symbol("}"))
        lhs = WLCall(:List, args)
    elseif t.kind == :op && t.val == :-
        arg = parse_expr(p, 25)
        lhs = arg isa Union{BigInt,Rational{BigInt}} ? -arg : WLCall(:Times, Any[BigInt(-1), arg])
    elseif t.kind == :op && t.val == :+
        lhs = parse_expr(p, 25)
    else
        error("wlparser: unexpected token $(t.val)")
    end

    while true
        t2 = peek(p)
        t2 === nothing && break
        t2.kind == :op || break
        op = t2.val
        if op == :+ || op == :-
            bp = 10
            bp <= minbp && break
            advance!(p)
            rhs = parse_expr(p, bp)
            rhs = op == :- ? WLCall(:Times, Any[BigInt(-1), rhs]) : rhs
            lhs = WLCall(:Plus, Any[lhs, rhs])
        elseif op == :* || op == :/
            bp = 20
            bp <= minbp && break
            advance!(p)
            rhs = parse_expr(p, bp)
            if op == :/
                if lhs isa Union{BigInt,Rational{BigInt}} && rhs isa Union{BigInt,Rational{BigInt}}
                    lhs = Rational{BigInt}(lhs) // Rational{BigInt}(rhs)
                else
                    lhs = WLCall(:Times, Any[lhs, WLCall(:Power, Any[rhs, BigInt(-1)])])
                end
            else
                lhs = WLCall(:Times, Any[lhs, rhs])
            end
        elseif op == :^
            bp = 30
            bp < minbp && break   # right-associative
            advance!(p)
            rhs = parse_expr(p, bp - 1)
            lhs = WLCall(:Power, Any[lhs, rhs])
        else
            break
        end
    end
    return lhs
end

"""Parse a single WL expression from a string."""
function wlparse(s::AbstractString)
    p = WLParser(wltokenize(s), 1)
    e = parse_expr(p)
    t = peek(p)
    (t === nothing || isop(t, Symbol(";"))) || error("wlparser: trailing tokens")
    return e
end

"""
Parse a file of `name = expr;` assignments (the PLD_database format).
Returns a Dict{Symbol,Any}. Comments `(* ... *)` are stripped first.
"""
function wlparsefile(path::AbstractString)
    src = read(path, String)
    src = replace(src, r"\(\*.*?\*\)"s => " ")
    toks = wltokenize(src)
    p = WLParser(toks, 1)
    out = Dict{Symbol,Any}()
    while peek(p) !== nothing
        t = advance!(p)
        t.kind == :sym || error("wlparser: expected assignment name, got $(t.val)")
        expect!(p, Symbol("="))
        out[t.val] = parse_expr(p)
        isop(peek(p), Symbol(";")) && advance!(p)
    end
    return out
end

# ---------------- AST -> Nemo polynomial ----------------

"""Evaluate a purely numeric AST to BigInt/Rational; return `nothing` if the
expression involves symbols."""
function wlnumeval(e)
    e isa Union{BigInt,Rational{BigInt}} && return e
    e isa WLCall || return nothing
    args = [wlnumeval(a) for a in e.args]
    any(isnothing, args) && return nothing
    if e.head == :Plus
        return sum(args)
    elseif e.head == :Times
        return prod(args)
    elseif e.head == :Power
        b, x = args
        x isa BigInt || return nothing
        v = x >= 0 ? Rational{BigInt}(b)^Int(x) : 1 // Rational{BigInt}(b)^Int(-x)
        return denominator(v) == 1 ? numerator(v) : v
    end
    return nothing
end

"""
    wlvariables(e) -> Set of variable keys

Collect symbols and indexed symbols (as `WLCall`) appearing in an AST,
excluding arithmetic heads.
"""
function wlvariables(e, acc=Set{Any}())
    if e isa Symbol
        push!(acc, e)
    elseif e isa WLCall
        if e.head in (:Plus, :Times, :Power, :List)
            foreach(a -> wlvariables(a, acc), e.args)
        else
            push!(acc, e)   # indexed symbol like m[1], treated as an atom
        end
    end
    return acc
end

"""
    varname(key) -> String

Canonical string name for a variable key (Symbol or indexed WLCall):
`m[1]` becomes "m1", plain symbols stay as-is.
"""
varname(k::Symbol) = string(k)
varname(k::WLCall) = string(k.head) * join(string.(k.args))

"""
    tonemo(e, vmap) -> polynomial

Convert an AST to an element of a Nemo ring, where `vmap` maps variable keys
(Symbols / indexed WLCalls) to ring generators. Powers must have non-negative
integer exponents.
"""
function tonemo(e, vmap::AbstractDict)
    if e isa BigInt
        R = parent(first(values(vmap)))
        return R(e)
    elseif e isa Rational{BigInt}
        R = parent(first(values(vmap)))
        return R(Nemo.QQ(e))
    elseif e isa Symbol
        haskey(vmap, e) || error("tonemo: unknown symbol $e")
        return vmap[e]
    elseif e isa WLCall
        if e.head == :Plus
            return sum(tonemo(a, vmap) for a in e.args)
        elseif e.head == :Times
            return prod(tonemo(a, vmap) for a in e.args)
        elseif e.head == :Power
            ex = wlnumeval(e.args[2])
            ex isa BigInt && ex >= 0 || error("tonemo: non-natural exponent $(e.args[2])")
            return tonemo(e.args[1], vmap)^Int(ex)
        elseif haskey(vmap, e)
            return vmap[e]
        else
            error("tonemo: unknown atom $(e.head)[...]")
        end
    else
        error("tonemo: cannot convert $(typeof(e))")
    end
end

"""
    loadsingularities(path) -> (polys, ring, varkeys)

Load a PLD_database file: builds a Nemo QQ multivariate ring over all
variables found and returns the singularity polynomials in it.
Variable keys are sorted by canonical name for reproducibility.
"""
function loadsingularities(path::AbstractString)
    d = wlparsefile(path)
    haskey(d, :singularities) || error("no `singularities` in $path")
    lst = d[:singularities]
    lst isa WLCall && lst.head == :List || error("`singularities` is not a list in $path")
    keys_ = sort!(collect(wlvariables(lst)); by=varname)
    R, gens_ = Nemo.polynomial_ring(Nemo.QQ, varname.(keys_))
    vmap = Dict{Any,Any}(zip(keys_, gens_))
    polys = [tonemo(a, vmap) for a in lst.args]
    return polys, R, keys_
end
