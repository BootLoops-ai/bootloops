# theta_g4_bridge.jl — ONE genus-4 acb_theta call through the Eichler wrapper
# src/siegel.jl (siegel_theta_all = FLINT acb_theta_all; used as-is, never
# edited).
# Bridge contract (abcount_s1.julia_theta_all_g4): EXACT dyadic
# mantissa/exponent integer pairs both ways; no decimal rounding anywhere.
#
# Input file lines:  prec P | g 4 | block re_man re_exp im_man im_exp rad_man rad_exp
# Output file lines: tau_echo k <8 ints> | theta idx <8 ints> | call_ns | load_ns | maxrss_bytes
# (the 8 ints per value: re_man re_exp im_man im_exp re_rad_man re_rad_exp im_rad_man im_rad_exp)

t_load0 = time_ns()
using Arblib
const SIEGEL_JL = get(ENV, "ABACUS_SIEGEL_JL", "")
isempty(SIEGEL_JL) && error("ABACUS_SIEGEL_JL unset — path to the Eichler src/siegel.jl wrapper is required")
include(SIEGEL_JL)
t_load1 = time_ns()

function arf_dump(rf; pbits::Int=4096)
    if iszero(rf)
        return (big(0), 0)
    end
    bf = BigFloat(rf; precision=pbits)
    e = exponent(bf)
    mi = BigInt(ldexp(bf, pbits - e))
    t = trailing_zeros(mi)
    mi >>= t
    return (mi, e - pbits + t)
end

function mag_dump(m; pbits::Int=256)
    rf = Arblib.Arf(prec=64)
    Arblib.set!(rf, m)
    return arf_dump(rf; pbits=pbits)
end

function exact_arb(man::BigInt, e::Int, prec::Int)
    x = Arb(man, prec=prec)
    Arblib.mul_2exp!(x, x, e)
    return x
end

inpath, outpath = ARGS[1], ARGS[2]
prec = 0
g = 0
blocks = Vector{NTuple{6,BigInt}}()
for line in eachline(inpath)
    w = split(line)
    isempty(w) && continue
    if w[1] == "prec"
        global prec = parse(Int, w[2])
    elseif w[1] == "g"
        global g = parse(Int, w[2])
    elseif w[1] == "block"
        push!(blocks, tuple([parse(BigInt, x) for x in w[2:7]]...))
    end
end
@assert g == 4 && length(blocks) == 4 && prec > 0

pb = prec + 64
tau = AcbMatrix(g, g; prec=prec)
for (k, b) in enumerate(blocks)
    rea = exact_arb(b[1], Int(b[2]), prec)
    ima = exact_arb(b[3], Int(b[4]), prec)
    z = Acb(rea, ima; prec=prec)
    if b[5] != 0
        r = exact_arb(b[5], Int(b[6]), prec)
        Arblib.add_error!(z, r)
    end
    tau[k, k] = z
end

z0 = AcbVector(g; prec=prec)
# ONE call, timed alone (pilot measurement)
t0 = time_ns()
th = siegel_theta_all(z0, tau, prec)
t1 = time_ns()

open(outpath, "w") do f
    for k in 1:g
        zz = tau[k, k]
        rm, re_ = arf_dump(Arblib.midref(Arblib.realref(zz)); pbits=pb)
        im_m, im_e = arf_dump(Arblib.midref(Arblib.imagref(zz)); pbits=pb)
        rrm, rre = mag_dump(Arblib.radref(Arblib.realref(zz)))
        irm, ire = mag_dump(Arblib.radref(Arblib.imagref(zz)))
        println(f, "tau_echo $k $rm $re_ $im_m $im_e $rrm $rre $irm $ire")
    end
    for (i, v) in enumerate(th)
        rm, re_ = arf_dump(Arblib.midref(Arblib.realref(v)); pbits=pb)
        im_m, im_e = arf_dump(Arblib.midref(Arblib.imagref(v)); pbits=pb)
        rrm, rre = mag_dump(Arblib.radref(Arblib.realref(v)))
        irm, ire = mag_dump(Arblib.radref(Arblib.imagref(v)))
        println(f, "theta $(i-1) $rm $re_ $im_m $im_e $rrm $rre $irm $ire")
    end
    println(f, "call_ns $(t1 - t0)")
    println(f, "load_ns $(t_load1 - t_load0)")
    println(f, "maxrss_bytes $(Sys.maxrss())")
end
