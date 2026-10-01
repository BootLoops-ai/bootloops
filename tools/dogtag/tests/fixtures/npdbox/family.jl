# MISLABEL (kept as the specimen this tool exists to catch): this family is the PLANAR Smirnov double box (hep-ph/9905323), NOT the nonplanar crossed box its label claims. Proof: k2->-k2 maps this propagator set exactly onto Smirnov's; amflow gives Im=0 at s<0,t<0,u>0 (planar signature; a crossed box would have an open u-cut -> nonzero Im). It duplicates the planar double box family.
#
# amflow family: 2-loop massless nonplanar double box ("crossed box").
#
# Propagators (Tausk hep-ph/9909506 conventions, checked vs Gehrmann-Remiddi hep-ph/0101124):
#   D1 = k1^2
#   D2 = (k1 + p1)^2
#   D3 = (k1 + p1 + p2)^2
#   D4 = k2^2
#   D5 = (k2 + p3 + p4)^2  = (k2 - p1 - p2)^2
#   D6 = (k2 + p4)^2       = (k2 - p1 - p2 - p3)^2
#   D7 = (k1 + k2)^2       <- the CROSSED joining line (nonplanar)
# Two ISPs (for IBP completeness):
#   D8 = (k1 + p4)^2       = (k1 - p1 - p2 - p3)^2  (ISP)
#   D9 = (k2 + p1)^2                                  (ISP)
#
# Kinematics: p1^2=p2^2=p3^2=p4^2=0, p1+p2+p3+p4=0.
# s = (p1+p2)^2, t = (p2+p3)^2.
# Sampling: set s = -1 (Euclidean), t = -1/x_val  so x = s/t = x_val.
#
# Top-sector indices: [1,1,1,1,1,1,1,0,0] (all 7 propagators on, ISPs off).

using LandauBootstrap, Arblib, Nemo, JSON3
using LandauBootstrap: magof

const NP_FAMILY = AmflowFamily(
    "np2L_massless",
    ["k1", "k2"],
    ["p1", "p2", "p3", "p4"],
    Dict("p4" => "-p1 - p2 - p3"),
    Dict("p1^2" => "0", "p2^2" => "0", "p3^2" => "0", "p4^2" => "0",
         "(p1 + p2)^2" => "s", "(p2 + p3)^2" => "t"),
    ["k1^2",
     "(k1 + p1)^2",
     "(k1 + p1 + p2)^2",
     "k2^2",
     "(k2 + p3 + p4)^2",
     "(k2 + p4)^2",
     "(k1 + k2)^2",
     "(k1 + p4)^2",
     "(k2 + p1)^2"]
)
# Top-sector: all 7 physical propagators on, both ISPs off
const NP_INDICES = [1, 1, 1, 1, 1, 1, 1, 0, 0]

const SAMPLES_DIR = "samples"

ratstr(x::Rational) = denominator(x) == 1 ? string(numerator(x)) :
                      string(numerator(x)) * "/" * string(denominator(x))

# Tag: np_x<num>_<den>  (encoding x = s/t = num/den)
np_tag(x::Rational) = "np_x$(numerator(x))_$(denominator(x))"

function np_save(tag::String, res::Dict{Int,Acb}, meta::Dict)
    mkpath(SAMPLES_DIR)
    d = Dict{String,Any}(string(k) => Dict("re" => string(real(v)),
                                           "im" => string(imag(v)))
                         for (k, v) in res)
    d["_meta"] = meta
    open(joinpath(SAMPLES_DIR, tag * ".json"), "w") do io
        JSON3.write(io, d)
    end
end

function np_load(tag::String, prec::Int)
    f = joinpath(SAMPLES_DIR, tag * ".json")
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

# Massless families do NOT need MASSIVE_FAMILY_OPTS (no mass thresholds).
# NO_EXTRA_OPTS is the default (empty extra options).
const NO_EXTRA_OPTS = (top = Dict{String,Any}(), blackbox = Dict{String,Any}())
