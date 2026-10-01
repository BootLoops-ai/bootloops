# theta9_theta.jl — certified Siegel theta-constant engine core (theta9 member
# of the eichler family node): the certified layer on Eichler.jl's acb_theta
# wrapper.
#
# WHAT IT DOES: reads a plain-text request (key/value lines), evaluates ALL
# genus-g theta constants theta[ab](0, tau) SQUARED in certified ball
# arithmetic (FLINT acb_theta via Arblib.theta_all!, Kieffer 2023 — proved
# error bounds), forms the frame-covariant ratio observables of the Watson
# deep_theta class (ratio2 = th2[i]/th2[j]; pairprod = th2[a]th2[b]/th2[c]th2[d]
# over disjoint even-char pairs), and emits everything as JSON with (mid_re,
# mid_im, rad_re, rad_im) per value. Gates are SCRIPT-EMITTED booleans:
#   gate_imtau_posdef   — Im tau positive-definite, certified (Arb leading minors)
#   gate_tau_symmetric  — input symmetry (constructed symmetric; recorded)
#   gate_odd_vanish     — every odd characteristic's theta(0)^2 ball contains 0
#                         (a theorem at z=0: free internal consistency check)
# Input tau may carry an explicit ball radius (tau_rad) — mode U propagates the
# frame/transport uncertainty through to every observable rigorously; tau_rad 0
# is mode P (map evaluation at the exact midpoint).
#
# INPUT (line-based, '#' comments):
#   label <string>            g <int>            prec_bits <int>
#   tau_rad <decimal>         (added to re AND im of every tau entry)
#   tau_<i>_<j>_re <decimal>  tau_<i>_<j>_im <decimal>   for 1<=i<=j<=g
#   emit_ratio2 <0|1>         emit_pairprod <0|1>   (pairprod: g=2 only)
#   digits_out <int>          (midpoint print digits; default prec_bits*0.30-2)
# OUTPUT: JSON on the path given as ARGS[2].
# USAGE: julia --project=<Eichler.jl project> theta9_theta.jl in.txt out.json
#        (theta9_stub.py resolves the project from the EICHLER_PROJECT env var,
#         defaulting to the repo's upgrades/Eichler.jl)

using Arblib

function read_req(path::String)
    d = Dict{String,String}()
    for ln in eachline(path)
        s = strip(ln)
        (isempty(s) || startswith(s, "#")) && continue
        k, v = split(s, limit = 2)
        d[strip(k)] = strip(v)
    end
    return d
end

# characteristic parity, FLINT encoding (a1..ag b1..bg, MSB first)
function char_is_even(ab::Integer, g::Int)
    s = 0
    for j in 1:g
        s += ((ab >> (2g - j)) & 1) * ((ab >> (g - j)) & 1)
    end
    return iseven(s)
end

# zero-radius decimal print of an Arf midpoint, correctly rounded to `digits`
function mid_str(m::Arblib.ArfOrRef, prec::Int, digits::Int)
    z = Arb(prec = prec)
    Arblib.set!(z, m)
    s = string(z, digits = digits)
    # arb_get_str formats: "1.23", "[1.23 +/- 4e-99]", "[+/- 1e-5]"
    if startswith(s, "[")
        body = s[2:end-1]
        i = findfirst("+/-", body)
        core = i === nothing ? body : strip(body[1:first(i)-1])
        isempty(core) && (core = "0")
        return core
    end
    return s
end

rad_str(x::Arblib.ArbOrRef) = string(Float64(Arblib.radref(x)))

function ball_json(v::Acb, prec::Int, digits::Int)
    re = real(v); im = imag(v)
    string("{\"m_re\": \"", mid_str(Arblib.midref(re), prec, digits),
        "\", \"m_im\": \"", mid_str(Arblib.midref(im), prec, digits),
        "\", \"rad_re\": ", rad_str(re), ", \"rad_im\": ", rad_str(im), "}")
end

function main()
    length(ARGS) == 2 || error("usage: theta9_theta.jl <in.txt> <out.json>")
    req = read_req(ARGS[1])
    label = get(req, "label", "unlabeled")
    g = parse(Int, req["g"])
    prec = parse(Int, req["prec_bits"])
    digits = haskey(req, "digits_out") ? parse(Int, req["digits_out"]) :
        max(20, floor(Int, prec * 0.30) - 2)
    emit_r2 = get(req, "emit_ratio2", "0") == "1"
    emit_pp = get(req, "emit_pairprod", "0") == "1"
    tau_rad_s = get(req, "tau_rad", "0")

    # ---- build tau (symmetric by construction) with the declared input ball
    tau = AcbMatrix(g, g; prec)
    radd = Arb(tau_rad_s; prec)
    for i in 1:g, j in i:g
        re = Arb(req["tau_$(i)_$(j)_re"]; prec)
        im = Arb(req["tau_$(i)_$(j)_im"]; prec)
        if !Arblib.is_zero(radd)
            Arblib.add_error!(re, radd)
            Arblib.add_error!(im, radd)
        end
        v = Acb(re, im; prec)
        tau[i, j] = v
        tau[j, i] = v
    end

    # ---- gate: Im tau positive definite (Sylvester, certified Arb dets)
    Y = ArbMatrix(g, g; prec)
    for i in 1:g, j in 1:g
        Y[i, j] = imag(tau[i, j])
    end
    posdef = true
    for k in 1:g
        Mk = ArbMatrix(k, k; prec)
        for i in 1:k, j in 1:k
            Mk[i, j] = Y[i, j]
        end
        dk = Arb(prec = prec)
        Arblib.det!(dk, Mk)
        posdef &= Arblib.is_positive(dk)
    end

    # ---- all theta constants SQUARED, certified (FLINT acb_theta)
    z = AcbVector(g; prec)
    nch = 1 << (2g)
    th2v = AcbVector(nch; prec)
    Arblib.theta_all!(th2v, z, tau, 1, prec)   # sqr = 1
    th2 = [Acb(th2v[i]; prec) for i in 1:nch]

    evens = [ab for ab in 0:(nch - 1) if char_is_even(ab, g)]
    odds = [ab for ab in 0:(nch - 1) if !char_is_even(ab, g)]
    odd_vanish = all(Arblib.contains_zero(real(th2[ab + 1])) &&
                     Arblib.contains_zero(imag(th2[ab + 1])) for ab in odds)

    # ---- observables on the EVEN squared constants, position-indexed in the
    # even-list order (for g=2 this order equals the reference deep_theta
    # comparator's char order: FLINT evens [0,1,2,3,4,6,8,9,12,15] = positions 0..9)
    ne = length(evens)
    e2 = [th2[ab + 1] for ab in evens]

    io = IOBuffer()
    print(io, "{\n\"label\": \"", label, "\", \"g\": ", g,
        ", \"prec_bits\": ", prec, ", \"digits_out\": ", digits,
        ", \"tau_rad\": \"", tau_rad_s, "\",\n")
    print(io, "\"engine\": \"FLINT acb_theta via Arblib.theta_all! (Eichler siegel layer; certified balls)\",\n")
    print(io, "\"chars_even_flint\": [", join(evens, ","), "],\n")
    print(io, "\"gates\": {\"imtau_posdef\": ", posdef,
        ", \"tau_symmetric_by_construction\": true, \"odd_vanish_contain_zero\": ",
        odd_vanish, ", \"n_odd_checked\": ", length(odds), "},\n")
    print(io, "\"tau_echo\": [")
    firste = true
    for i in 1:g, j in i:g
        firste || print(io, ", ")
        firste = false
        print(io, "{\"ij\": [", i, ",", j, "], \"v\": ", ball_json(tau[i, j], prec, digits), "}")
    end
    print(io, "],\n\"theta_sq_even\": [")
    for (p, ab) in enumerate(evens)
        p > 1 && print(io, ", ")
        print(io, "{\"pos\": ", p - 1, ", \"char\": ", ab, ", \"v\": ",
            ball_json(e2[p], prec, digits), "}")
    end
    print(io, "]")
    if emit_r2
        print(io, ",\n\"ratio2\": [")
        first2 = true
        for i in 1:ne, j in 1:ne
            i == j && continue
            r = e2[i] / e2[j]
            first2 || print(io, ", ")
            first2 = false
            print(io, "{\"ij\": [", i - 1, ",", j - 1, "], \"v\": ", ball_json(r, prec, digits), "}")
        end
        print(io, "]")
    end
    if emit_pp
        g == 2 || error("emit_pairprod: g=2 fan only (explicit quad requests for higher g)")
        print(io, ",\n\"pairprod\": [")
        firstp = true
        for ia in 1:ne, ib in (ia + 1):ne, ic in 1:ne, id in (ic + 1):ne
            (ia == ic || ia == id || ib == ic || ib == id) && continue
            r = (e2[ia] * e2[ib]) / (e2[ic] * e2[id])
            firstp || print(io, ", ")
            firstp = false
            print(io, "{\"quad\": [", ia - 1, ",", ib - 1, ",", ic - 1, ",", id - 1,
                "], \"v\": ", ball_json(r, prec, digits), "}")
        end
        print(io, "]")
    end
    print(io, "\n}\n")
    open(ARGS[2], "w") do f
        write(f, String(take!(io)))
    end
    println("theta9_theta.jl DONE label=", label, " g=", g, " prec_bits=", prec,
        " posdef=", posdef, " odd_vanish=", odd_vanish)
end

main()
