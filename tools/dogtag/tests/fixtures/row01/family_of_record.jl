# Shared definitions for the OUTER-MASS planar double box oracle.
# Conventions:
#
#   1404.2922:  u = 4m^2/(-s), v = 4m^2/(-t);  bu = sqrt(1+u) etc.
#   dbox = GG{1,1,1,0}{1,0,1,1}{1}; in D=4 it is finite ("I5" of 2507.17815);
#   g10  = -(1/8) s^2 t bu buv I5   (Euclidean s,t < 0, m = 1)
#        =  (1/8) sqrt(s(s-4m^2)) sqrt(st(st-4m^2(s+t))) I5, pure weight 4,
#   boundary g10 -> 0 as u,v -> inf;  ANCHOR  g10(4,4) = 0.05303862978353642964
#   06023741052... (1404.2922 sec. Checks).
#   amflow J = I5 series in eps (measure d^Dk/(i pi^{D/2}) per loop); the
#   eps^0 coefficient IS I5 (integral finite, gammaE dressing irrelevant
#   at order eps^0).
#
# Alphabet: 12 letters of 2410.02424 Tab. tab:double_box_letters, expressed
# multiplicatively over the 12 (w,z) factors (rationalization
# u=(1-w^2)(1-z^2)/(w-z)^2, v=4wz/(w-z)^2, bu=(1-wz)/(w-z), bv=(w+z)/(w-z),
# buv=(1+wz)/(w-z); region u,v>0 <-> 0<z<w<1).

using LandauBootstrap, Arblib, Nemo, JSON3
using LandauBootstrap: magof, ub

const DBOM_ROOT = "LandauBootstrap"
const DBOM_SAMPLES = joinpath(DBOM_ROOT, "results", "samples")
const DBOM_BASIS_JSON = joinpath(DBOM_ROOT, "results", "dbox-outermass-basis.json")

# ----------------------------------------------------------------- alphabet --
const DBOM_FACTORS = ["w", "z", "1-w", "1+w", "1-z", "1+z", "w-z", "w+z",
                      "1-w*z", "1+w*z", "1-w+z+w*z", "1+w-z+w*z"]
# exponent vectors of L1..L12 over the factors (constants drop at dlog level)
const DBOM_LETTERS = [
    ("L1_u",    [0, 0, 1, 1, 1, 1, -2, 0, 0, 0, 0, 0]),
    ("L2_v",    [1, 1, 0, 0, 0, 0, -2, 0, 0, 0, 0, 0]),
    ("L3_1pu",  [0, 0, 0, 0, 0, 0, -2, 0, 2, 0, 0, 0]),
    ("L4_1pv",  [0, 0, 0, 0, 0, 0, -2, 2, 0, 0, 0, 0]),
    ("L5_upv",  [0, 0, 0, 0, 0, 0, -2, 0, 0, 0, 1, 1]),
    ("L6_bu",   [0, 0, 1, -1, -1, 1, 0, 0, 0, 0, 0, 0]),
    ("L7_bv",   [-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
    ("L8_buv",  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1]),
    ("L9_uvu",  [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
    ("L10_uvv", [0, 0, 1, -1, 1, -1, 0, 0, 0, 0, 0, 0]),
    ("L11_1puv",[0, 0, 0, 0, 0, 0, -2, 0, 0, 2, 0, 0]),
    ("L12_all", [-1, 1, 1, 1, -1, -1, 0, 0, 0, 0, 0, 0]),
]

# Galois parity (sigma_u, sigma_v, sigma_uv) of each letter; the function
# multiplying the prefactor 1/(bu*buv) must live in sector (1,0,1).
const DBOM_PARITY = [(0,0,0),(0,0,0),(0,0,0),(0,0,0),(0,0,0),
                     (1,0,0),(0,1,0),(0,0,1),(1,0,1),(0,1,1),(0,0,0),(1,1,1)]
const DBOM_TARGET_SECTOR = (1, 0, 1)

# first-entry combos (physical logarithmic branch cuts, 2410.02424):
# {L1/L3, L2/L4, L3 L8/(L5 L10), L4 L8/(L5 L9), L6, L7, L12}
const DBOM_FE_COMBOS = [
    ("C1", Dict(1 => 1, 3 => -1)),
    ("C2", Dict(2 => 1, 4 => -1)),
    ("C3", Dict(3 => 1, 8 => 1, 5 => -1, 10 => -1)),
    ("C4", Dict(4 => 1, 8 => 1, 5 => -1, 9 => -1)),
    ("C5", Dict(6 => 1)),
    ("C6", Dict(7 => 1)),
    ("C7", Dict(12 => 1)),
]

# published symbol (2410.02424 eq:symbol), S(I~_dbox), expanded to raw words
const DBOM_PUB_SYMBOL = [
    (-1, (6,1,6,9)), (1, (6,3,6,9)),
    (-1, (6,1,9,6)), (1, (6,3,9,6)),
    (1, (6,6,1,9)), (1, (6,6,2,9)), (-1, (6,6,3,9)), (-1, (6,6,5,9)),
    (1, (6,9,2,6)), (-1, (6,9,5,6)),
    (1, (6,6,8,6)),
    (1, (6,9,8,9)),
    (1, (7,10,2,6)), (-1, (7,10,5,6)),
    (1, (7,10,8,9)),
    (1, (7,7,1,9)), (-1, (7,7,5,9)),
    (1, (7,7,8,6)),
]

# ------------------------------------------------------------------- family --
# 1404.2922 two-loop family (def-2loopfamily), masses m^2 = msq on the six
# ring lines; k1: D1 D2 D3 (+ISP D4), k2: D1 D3 D4 (+ISP D2), rung massless.
# D1(k)=-k^2+m^2, D2(k)=-(k+p1)^2+m^2, D3(k)=-(k+p1+p2)^2+m^2, D4=-(k-p4)^2+m^2
# (overall signs of inverse propagators are irrelevant to amflow input, which
# takes (momentum^2 - mass^2)-type propagators; J differs by (-1)^7 = -1:
#   GG = - J  for indices (1,1,1,1,1,1,1,0,0)).
const DBOM_FAMILY = AmflowFamily(
    "dboxom", ["k1", "k2"], ["p1", "p2", "p3", "p4"],
    Dict("p4" => "-p1 - p2 - p3"),
    Dict("p1^2" => "0", "p2^2" => "0", "p3^2" => "0", "p4^2" => "0",
         "(p1 + p2)^2" => "s", "(p2 + p3)^2" => "t"),
    ["k1^2 - msq", "(k1 + p1)^2 - msq", "(k1 + p1 + p2)^2 - msq",
     "k2^2 - msq", "(k2 + p1 + p2)^2 - msq", "(k2 + p1 + p2 + p3)^2 - msq",
     "(k1 - k2)^2",
     "(k1 + p1 + p2 + p3)^2", "(k2 + p1)^2"])
const DBOM_INDICES = [1, 1, 1, 1, 1, 1, 1, 0, 0]

# Euclidean points (u,v), u=4msq/(-s), v=4msq/(-t); msq=1 -> s=-4/u, t=-4/v.
# Series evaluation of the basis needs min(u,v) > 1 (kernel convergence radius
# sqrt(min(u,v)) in mu); grid kept within [3/2, 9].
ratstr(x::Rational) = denominator(x) == 1 ? string(numerator(x)) :
                      string(numerator(x)) * "/" * string(denominator(x))
point_st(uv::Tuple{Rational{Int},Rational{Int}}) =
    Dict("s" => ratstr(-4 // uv[1]), "t" => ratstr(-4 // uv[2]), "msq" => "1")
dbom_tag(uv) = "dboxom_u$(numerator(uv[1]))_$(denominator(uv[1]))" *
               "_v$(numerator(uv[2]))_$(denominator(uv[2]))"

# g10 normalization from the amflow eps^0 value I5 at (u,v), m=1:
#   g10 = -(1/8) s^2 t bu buv I5,  s=-4/u, t=-4/v
function g10_from_I5(I5::Arb, uv; prec::Int)
    u = Arb(Rational{BigInt}(uv[1]), prec = prec)
    v = Arb(Rational{BigInt}(uv[2]), prec = prec)
    s = -4 / u
    t = -4 / v
    bu = sqrt(1 + u)
    buv = sqrt(1 + u + v)
    -(s^2 * t * bu * buv / 8) * I5
end

# ----------------------------------------------------------------- sample IO --
function dbom_save_sample(tag::String, res::Dict{Int,Acb}, meta::Dict)
    mkpath(DBOM_SAMPLES)
    d = Dict{String,Any}(string(k) => Dict("re" => string(real(v)),
                                           "im" => string(imag(v)))
                         for (k, v) in res)
    d["_meta"] = meta
    open(joinpath(DBOM_SAMPLES, tag * ".json"), "w") do io
        JSON3.write(io, d)
    end
end

function dbom_load_sample(tag::String, prec::Int)
    f = joinpath(DBOM_SAMPLES, tag * ".json")
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
