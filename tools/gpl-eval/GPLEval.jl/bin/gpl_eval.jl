#!/usr/bin/env julia
# CLI:  julia --project bin/gpl_eval.jl --alphabet alpha.json --basis basis.json \
#             --point '{"s":-3,"t":-7,"z":0.5}' --prec 256 --out values.json
#
# alpha.json:  {"letters": {"l0": "0", "l1": "1", "lm1": "-1", "ls": "s+t"}, "zvar": "z"}
#   Each letter is a Julia expression in the kinematic variables (keys of --point);
#   it is the GPL INDEX aᵢ (root of the dlog letter). "zvar" names the integration
#   variable in --point.
# basis.json:  [["l0","l1"], ["l0","l1","lm1"], ...]   (list of symbol words)
# values.json: {"prec_bits":..., "values":[{"re":"...","im":"...","rad":"..."}...]}

using GPLEval, Arblib, JSON3

function _parse_args(args)
    d = Dict{String,String}()
    i = 1
    while i ≤ length(args)
        startswith(args[i], "--") || error("bad arg $(args[i])")
        d[args[i][3:end]] = args[i+1]; i += 2
    end
    d
end

# very small expression evaluator over a Dict of BigFloat values
function _eval_expr(s::AbstractString, env::Dict{String,BigFloat})
    ex = Meta.parse(s)
    ev(e) = e isa Number ? big(e) :
            e isa Symbol ? env[string(e)] :
            e isa Expr && e.head == :call ?
                (e.args[1] == :+ ? sum(ev.(e.args[2:end])) :
                 e.args[1] == :- ? (length(e.args) == 2 ? -ev(e.args[2]) :
                                    ev(e.args[2]) - ev(e.args[3])) :
                 e.args[1] == :* ? prod(ev.(e.args[2:end])) :
                 e.args[1] == :/ ? ev(e.args[2]) / ev(e.args[3]) :
                 e.args[1] == :^ ? ev(e.args[2])^ev(e.args[3]) :
                 e.args[1] == :sqrt ? sqrt(Complex(ev(e.args[2]))) :
                 error("op $(e.args[1])")) :
            error("bad expr $e")
    ev(ex)
end

function main(args)
    d = _parse_args(args)
    prec = parse(Int, get(d, "prec", "768"))
    digits = ceil(Int, prec * log10(2.0))
    alpha = JSON3.read(read(d["alphabet"], String))
    basis = JSON3.read(read(d["basis"], String))
    point = JSON3.read(d["point"])
    setprecision(BigFloat, prec + 64)
    env = Dict{String,BigFloat}(string(k) => big(Float64(v)) for (k, v) in pairs(point))
    # build alphabet :: Dict{Symbol,Function}
    A = Dict{Symbol,Function}()
    for (name, expr) in pairs(alpha.letters)
        let s = String(expr)
            A[Symbol(name)] = pt -> _eval_expr(s, pt)
        end
    end
    zvar = String(get(alpha, :zvar, "z"))
    z = env[zvar]
    words = [Symbol.(w) for w in basis]
    vals = evaluate_basis(words, A, z, env; prec = prec)
    out = Dict("prec_bits" => prec, "digits" => digits, "point" => point,
               "values" => [Dict("re" => string(real(v); digits = digits),
                                 "im" => string(imag(v); digits = digits),
                                 "rad" => string(Arblib.radref(real(v))))
                            for v in vals])
    haskey(d, "out") ? open(io -> JSON3.write(io, out), d["out"], "w") :
                       JSON3.write(stdout, out)
end

isinteractive() || main(ARGS)
