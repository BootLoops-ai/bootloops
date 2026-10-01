# period_g2_bridge.jl — S2 bridge: complex-chain certified genus-2 periods,
# full-matrix genus-g theta, Thomae/covariant tie data — through the Eichler
# wrapper src/siegel.jl (siegel_theta_all, thomae_sextic,
# igusa_clebsch, HyperellipticCurve; the Eichler files stay UNTOUCHED — this file
# assembles the wrapper's exposed-but-unassembled complex-chain primitives per
# battery/SEEDS.json S2_period_construction).
#
# Bridge contract (abcount_s2.bridge_call): EXACT dyadic mantissa/exponent
# integer pairs both ways; no decimal rounding anywhere (S1 bridge discipline).
#
# Input file lines (first line selects the mode):
#   mode periods | prec P | coeffs c0 .. c5 | chain i1 i2 i3 i4 i5
#   mode theta | prec P | g G | entry i j re_man re_exp im_man im_exp   (upper tri, exact)
#   mode thomae | prec P | entry i j ...                                 (2x2 tau)
#   mode covariants | prec P | coeffs c0 .. c5
# Output lines (8 ints per ball: re_man re_exp im_man im_exp re_rad_man
# re_rad_exp im_rad_man im_rad_exp):
#   periods:    root k <8> | pi k j <8> | pilot_integral_ns | sweep_ns | load_ns | maxrss_bytes
#   theta:      tau_echo i j <4 mid ints> | theta idx <8> | call_ns | maxrss_bytes
#   thomae:     chi5_nonzero true|false | troot k <8> | cov k <8>
#   covariants: cov k <8>

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

function dump_ball(io, tag, v; pbits::Int=4096)
    rm, re_ = arf_dump(Arblib.midref(Arblib.realref(v)); pbits=pbits)
    im_m, im_e = arf_dump(Arblib.midref(Arblib.imagref(v)); pbits=pbits)
    rrm, rre = mag_dump(Arblib.radref(Arblib.realref(v)))
    irm, ire = mag_dump(Arblib.radref(Arblib.imagref(v)))
    println(io, tag, " ", rm, " ", re_, " ", im_m, " ", im_e, " ",
            rrm, " ", rre, " ", irm, " ", ire)
end

# ---- input parse -------------------------------------------------------------

inpath, outpath = ARGS[1], ARGS[2]
mode = ""
prec = 0
g = 2
budget_ns = 0
coeffs = Int[]
chain = Int[]
entries = Tuple{Int,Int,BigInt,Int,BigInt,Int,BigInt,Int}[]
for line in eachline(inpath)
    w = split(line)
    isempty(w) && continue
    if w[1] == "mode"
        global mode = w[2]
    elseif w[1] == "prec"
        global prec = parse(Int, w[2])
    elseif w[1] == "g"
        global g = parse(Int, w[2])
    elseif w[1] == "budget_ns"
        global budget_ns = parse(Int, w[2])
    elseif w[1] == "coeffs"
        global coeffs = [parse(Int, x) for x in w[2:end]]
    elseif w[1] == "chain"
        global chain = [parse(Int, x) for x in w[2:end]]
    elseif w[1] == "entry"
        # 6 numbers = exact dyadic mid; 8 numbers = mid + radius (add_error)
        rm_, re_ = length(w) >= 9 ? (parse(BigInt, w[8]), parse(Int, w[9])) : (big(0), 0)
        push!(entries, (parse(Int, w[2]), parse(Int, w[3]),
                        parse(BigInt, w[4]), parse(Int, w[5]),
                        parse(BigInt, w[6]), parse(Int, w[7]), rm_, re_))
    end
end
@assert prec > 0 && mode != ""
pb = prec + 64

function tau_from_entries(gg::Int)
    tau = AcbMatrix(gg, gg; prec=prec)
    for (i, j, rm, re_, im_m, im_e, rad_m, rad_e) in entries
        z = Acb(exact_arb(rm, re_, prec), exact_arb(im_m, im_e, prec); prec=prec)
        if rad_m != 0
            Arblib.add_error!(z, exact_arb(rad_m, rad_e, prec))
        end
        tau[i, j] = z
        tau[j, i] = z
    end
    return tau
end

# ---- certified chord integral (SEEDS segment_integral_scheme) ----------------

function chord_period(lead, roots, ia, ib, k, prec)
    a, b = roots[ia], roots[ib]
    m = (a + b) / 2
    dl = (b - a) / 2
    others = [roots[l] for l in eachindex(roots) if l != ia && l != ib]
    us = Acb[]
    for r in others
        w1 = a - r; w2 = b - r
        u = w1 / abs(w1) + w2 / abs(w2)
        u = u / abs(u)
        (Arblib.is_positive(real(w1 / u)) && Arblib.is_positive(real(w2 / u))) ||
            error("FAIL-HALF-PLANE: certificate failed for chord ($ia,$ib)")
        push!(us, u)
    end
    K = sqrt(-lead * dl^2 * prod(us))
    integrand = let m = m, dl = dl, others = others, us = us, K = K, k = k, prec = prec
        (th; analytic::Bool = false) -> begin
            x = m + dl * sin(th)
            den = Acb(K; prec)
            bad = false
            for (r, u) in zip(others, us)
                z = (x - r) / u
                if analytic && !Arblib.is_positive(real(z))
                    bad = true
                    break
                end
                den *= sqrt(z)
            end
            if bad
                ind = Acb(prec = prec); Arblib.indeterminate!(ind); return ind
            end
            x^k * dl / den
        end
    end
    2 * Arblib.integrate(integrand, Acb(-Arb(π; prec) / 2; prec),
                         Acb(Arb(π; prec) / 2; prec);
                         check_analytic = true, prec = prec,
                         rtol = exp10(-0.28 * prec))
end

# ---- modes -------------------------------------------------------------------

open(outpath, "w") do f
    if mode == "periods"
        @assert length(coeffs) == 6 && length(chain) == 5
        C = HyperellipticCurve(coeffs, prec)
        ord = sortperm([(Float64(real(r)), Float64(imag(r))) for r in C.roots])
        rts = C.roots[ord]
        for (i, r) in enumerate(rts)
            dump_ball(f, "root $i", r; pbits=pb)
        end
        # ONE pilot integral timed alone; price the 15-integral sweep
        # BEFORE launching it (SEEDS S2_walls.pilot_integral_pricing)
        t0 = time_ns()
        P11 = chord_period(C.lead, rts, chain[1], chain[2], 0, prec)
        t1 = time_ns()
        if budget_ns > 0 && 15 * (t1 - t0) > budget_ns
            println(f, "pilot_integral_ns ", t1 - t0)
            println(f, "stop_wall projected_sweep_ns ", 15 * (t1 - t0),
                    " budget_ns ", budget_ns)
            println(f, "load_ns ", t_load1 - t_load0)
            println(f, "maxrss_bytes ", Sys.maxrss())
            exit(0)
        end
        P = Matrix{Acb}(undef, 2, 4)
        P[1, 1] = P11
        for j in 1:4, k in 0:1
            (j == 1 && k == 0) && continue
            P[k + 1, j] = chord_period(C.lead, rts, chain[j], chain[j + 1], k, prec)
        end
        t2 = time_ns()
        for k in 1:2, j in 1:4
            dump_ball(f, "pi $k $j", P[k, j]; pbits=pb)
        end
        println(f, "pilot_integral_ns ", t1 - t0)
        println(f, "sweep_ns ", t2 - t1)
        println(f, "load_ns ", t_load1 - t_load0)
        println(f, "maxrss_bytes ", Sys.maxrss())
    elseif mode == "theta"
        tau = tau_from_entries(g)
        for (i, j, rm, re_, im_m, im_e, _, _) in entries
            erm, ere = arf_dump(Arblib.midref(Arblib.realref(tau[i, j])); pbits=pb)
            eim, eie = arf_dump(Arblib.midref(Arblib.imagref(tau[i, j])); pbits=pb)
            println(f, "tau_echo $i $j $erm $ere $eim $eie")
        end
        z0 = AcbVector(g; prec=prec)
        t0 = time_ns()
        th = siegel_theta_all(z0, tau, prec)
        t1 = time_ns()
        for (i, v) in enumerate(th)
            dump_ball(f, "theta $(i-1)", v; pbits=pb)
        end
        println(f, "call_ns ", t1 - t0)
        println(f, "maxrss_bytes ", Sys.maxrss())
    elseif mode == "thomae"
        tau = tau_from_entries(2)
        fs, chi5 = thomae_sextic(tau, prec)
        println(f, "chi5_nonzero ", !Arblib.contains_zero(chi5) ? "true" : "false")
        cf(i) = (a = Acb(prec = prec); Arblib.get_coeff!(a, fs, i); a)
        scoeffs = [cf(i) for i in 0:6]
        Arblib.contains_zero(scoeffs[7]) &&
            error("FAIL-THOMAE-DEGENERATE: sextic leading coefficient contains 0")
        p = AcbPoly(prec = prec)
        for (i, c) in enumerate(scoeffs)
            Arblib.set_coeff!(p, i - 1, c)
        end
        roots = AcbVector(6; prec)
        n = ccall((:acb_poly_find_roots, Arblib.libflint), Clong,
            (Ptr{Arblib.acb_struct}, Ref{Arblib.acb_poly_struct}, Ptr{Cvoid}, Clong, Clong),
            roots, p, C_NULL, 0, prec)
        n == 6 || error("FAIL-THOMAE-ROOTS: could not certify 6 isolated roots")
        for i in 1:6
            dump_ball(f, "troot $i", roots[i]; pbits=pb)
        end
        ict = igusa_clebsch(scoeffs, prec)
        for (i, c) in enumerate(ict)
            dump_ball(f, "cov $i", c; pbits=pb)
        end
        println(f, "maxrss_bytes ", Sys.maxrss())
    elseif mode == "covariants"
        @assert length(coeffs) == 6
        C = HyperellipticCurve(coeffs, prec)
        icf = igusa_clebsch(C)
        for (i, c) in enumerate(icf)
            dump_ball(f, "cov $i", c; pbits=pb)
        end
        println(f, "maxrss_bytes ", Sys.maxrss())
    else
        error("FAIL-INTERNAL: unknown mode $mode")
    end
end
