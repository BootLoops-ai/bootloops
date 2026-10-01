# amflow family: planar double box, ONE internal mass on the CENTRAL RUNG.
# 4 massless external legs. mt^2 = msq is the only mass.
#
# Routing mirrors the validated outer-mass dbox family (dbox_outermass_common.jl)
# but moves the mass: the 6 RING lines are MASSLESS, the rung (line 7,
# (k1-k2)^2) carries -msq. ISPs are lines 8,9 (massless).
#
#   D1(k1)=k1^2, D2=(k1+p1)^2, D3=(k1+p1+p2)^2,
#   D4(k2)=k2^2, D5=(k2+p1+p2)^2, D6=(k2+p1+p2+p3)^2,
#   D7=(k1-k2)^2 - msq   <-- the massive central rung
#   D8=(k1+p1+p2+p3)^2 (ISP), D9=(k2+p1)^2 (ISP)
# s=(p1+p2)^2, t=(p2+p3)^2.

using LandauBootstrap, Arblib, Nemo, JSON3
using LandauBootstrap: magof

const CR_FAMILY = AmflowFamily(
    "dbox_rungmass", ["k1", "k2"], ["p1", "p2", "p3", "p4"],
    Dict("p4" => "-p1 - p2 - p3"),
    Dict("p1^2" => "0", "p2^2" => "0", "p3^2" => "0", "p4^2" => "0",
         "(p1 + p2)^2" => "s", "(p2 + p3)^2" => "t"),
    ["k1^2", "(k1 + p1)^2", "(k1 + p1 + p2)^2",
     "k2^2", "(k2 + p1 + p2)^2", "(k2 + p1 + p2 + p3)^2",
     "(k1 - k2)^2 - msq",
     "(k1 + p1 + p2 + p3)^2", "(k2 + p1)^2"])
const CR_INDICES = [1, 1, 1, 1, 1, 1, 1, 0, 0]   # top sector, scalar

ratstr(x::Rational) = denominator(x) == 1 ? string(numerator(x)) :
                      string(numerator(x)) * "/" * string(denominator(x))
const CR_SAMPLES = "samples"

cr_tag(s::Rational, t::Rational) =
    "cr_s$(numerator(s))_$(denominator(s))_t$(numerator(t))_$(denominator(t))"

function cr_save(tag::String, res::Dict{Int,Acb}, meta::Dict)
    mkpath(CR_SAMPLES)
    d = Dict{String,Any}(string(k) => Dict("re" => string(real(v)),
                                           "im" => string(imag(v)))
                         for (k, v) in res)
    d["_meta"] = meta
    open(joinpath(CR_SAMPLES, tag * ".json"), "w") do io
        JSON3.write(io, d)
    end
end

function cr_load(tag::String, prec::Int)
    f = joinpath(CR_SAMPLES, tag * ".json")
    isfile(f) || return nothing
    d = JSON3.read(read(f, String))
    out = Dict{Int,Acb}()
    for (k, v) in d
        k == :_meta && continue
        out[parse(Int, string(k))] = Acb(parse_ball(string(v["re"]), prec),
                                         parse_ball(string(v["im"]), prec),
                                         prec = prec)
    end
    out
end
