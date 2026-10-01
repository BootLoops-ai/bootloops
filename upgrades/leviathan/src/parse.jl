# parse.jl — tiny arithmetic-expression -> Oscar polynomial evaluator.
# Supports + - * ^ and / by numeric constants; integers; variables from a symbol table.

function poly_from_string(str::AbstractString, table::Dict{Symbol,T}, R) where {T}
    _pev(Meta.parse(String(str)), table, R)
end

function _pev(ex, tab, R)
    if ex isa Integer
        return R(Nemo.ZZ(ex))
    elseif ex isa Expr && ex.head == :macrocall && ex.args[end] isa String
        return R(Nemo.ZZ(ex.args[end]))   # @int128_str / @big_str literal (GlobalRef in args[1])
    elseif ex isa Symbol
        haskey(tab, ex) || error("unknown symbol $ex in G")
        return tab[ex]
    elseif ex isa Expr && ex.head == :call
        op = ex.args[1]
        as = ex.args[2:end]
        if op == :+
            return sum(_pev(a, tab, R) for a in as)
        elseif op == :-
            length(as) == 1 && return -_pev(as[1], tab, R)
            return _pev(as[1], tab, R) - sum(_pev(a, tab, R) for a in as[2:end])
        elseif op == :*
            return prod(_pev(a, tab, R) for a in as)
        elseif op == :^
            return _pev(as[1], tab, R)^Int(as[2])
        elseif op == :/
            den = as[2]
            den isa Integer || error("only numeric denominators supported, got $den")
            return _pev(as[1], tab, R) * inv(Nemo.QQ(den)) * one(R)
        elseif op == ://
            n, d = _pev(as[1], tab, R), _pev(as[2], tab, R)
            return R(Oscar.constant_coefficient(n) // Oscar.constant_coefficient(d))
        end
        error("unsupported operator $op")
    end
    error("unsupported expression $ex")
end
