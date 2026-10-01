# Production gauntlet: certified high-precision runs + cost-vs-digits curves.
# Writes gauntlet/results-production.json (flat text, parseable).

using Arblib, Eichler, Nemo
const E = Eichler

results = Dict{String,Any}()
out = open(joinpath(@__DIR__, "results-production.txt"), "w")
lk = ReentrantLock()
function rec(key, val)
    lock(lk) do
        println(out, key, " = ", val)
        flush(out)
    end
    println(key, " = ", val)
end

certified_digits(v::Arb) = floor(Int, max(0, Arblib.rel_accuracy_bits(v)) * 0.30103)
certified_digits(v::Acb) = floor(Int, max(0, Arblib.rel_accuracy_bits(v)) * 0.30103)

# how many digits we need, bits, q-series order at |q| = 0.08
specs = [(100, 512, 300), (200, 832, 550), (300, 1152, 800), (500, 1856, 1450)]

# ---- (1) Gamma16 build + sunrise cost curve ----------------------------------
for (D, prec, N) in specs
    t0 = time(); G = Gamma16(prec, N); t1 = time() - t0
    rec("gamma16_build_seconds_D$D", round(t1, digits = 2))
    t0 = time(); B = sunrise_boundary_constants(prec, 12); tB = time() - t0
    rec("boundary_constants_J12_seconds_D$D", round(tB, digits = 2))
    L2 = E.dirichlet_L_chim3(2, prec)
    target = Acb(3 * sqrt(Arb(3; prec = prec)) / 2 * L2; prec = prec)
    rec("gauntlet_a_overlap_D$D", Arblib.overlaps(B[1], target))
    rec("gauntlet_a_certified_digits_D$D", certified_digits(real(B[1])))
    t0 = time(); S4 = sunrise4(G, Acb(-1; prec = prec); J = 10, ); tS = time() - t0
    rec("sunrise4_J10_seconds_D$D", round(tS, digits = 2))
    rec("sunrise4_eps0_certified_digits_D$D", certified_digits(real(E.eps_coeff(S4, 0))))
    rec("sunrise4_eps7_certified_digits_D$D", certified_digits(real(E.eps_coeff(S4, 7))))
    if D == 500
        # record the 500-digit production values, eps^-2..7
        for k in -2:7
            v = real(-E.eps_coeff(S4, k))
            rec("sunrise4_amflowconv_eps$(k)_D500", string(v))
        end
        rec("B0_D500", string(real(B[1])))
        rec("L_chim3_2_D500", string(L2))
    end
end

# ---- (2) curve layer cost curve ----------------------------------------------
for (D, prec, N) in specs
    s = Nemo.QQ(7, 3)
    C = QuarticCurve([Nemo.QQ(0), -4 * s^2, s^2 - 4 * s, 2 * s, Nemo.QQ(1)], prec)
    ord = sortperm([Float64(real(r)) for r in C.roots])
    t0 = time(); Pq = period_quadrature(C, (ord[2], ord[3])); tq = time() - t0
    t0 = time(); ψ1, ψ2 = curve_periods(C); tk = time() - t0
    rec("period_quadrature_seconds_D$D", round(tq, digits = 2))
    rec("period_Kformula_seconds_D$D", round(tk, digits = 4))
    rec("period_quad_certified_digits_D$D", certified_digits(real(Pq)))
    rec("period_match_D$D", Arblib.overlaps(Pq, ψ1))
end

# ---- (3) PF transport cost ----------------------------------------------------
for (D, prec, N) in specs
    t0 = time()
    tr, dir = frobenius_transport_acm(Arb(2; prec = prec), Arb(5; prec = prec); M = 64)
    rec("pf_transport_seconds_D$D", round(time() - t0, digits = 2))
    rec("pf_transport_overlap_D$D", Arblib.overlaps(tr[1], dir))
    rec("pf_transport_certified_digits_D$D", certified_digits(real(tr[1])))
end

# ---- (4) gauntlet (c) residuals at 5 points, high precision -------------------
let prec = 1856
    for (i, sv) in enumerate([Nemo.QQ(1,2), Nemo.QQ(3), Nemo.QQ(25,7), Nemo.QQ(50), Nemo.QQ(7,3)])
        s = Arb(BigInt(numerator(sv)); prec = prec) / Arb(BigInt(denominator(sv)); prec = prec)
        r = picard_fuchs_residual_acm(s)
        rec("pf_residual_contains_zero_$i", Arblib.contains_zero(real(r)) && Arblib.contains_zero(imag(r)))
        rec("pf_residual_radius_$i", Float64(Arblib.radius(Arb, real(r))))
    end
end

# ---- (5) gauntlet (d) at 5 random rational points, 500 digits -----------------
let prec = 1856
    for (i, (num, den)) in enumerate([(7,3),(15,7),(101,13),(3,11),(89,2)])
        s = Nemo.QQ(num, den)
        C = QuarticCurve([Nemo.QQ(0), -4*s^2, s^2-4*s, 2*s, Nemo.QQ(1)], prec)
        ord = sortperm([Float64(real(r)) for r in C.roots])
        Pq = period_quadrature(C, (ord[2], ord[3]))
        ψ1, _ = curve_periods(C)
        r = [real(C.roots[j]) for j in ord]
        agm = Arblib.agm!(Acb(prec=prec), Acb(1;prec=prec), sqrt(Acb(1-((r[3]-r[2])*(r[4]-r[1]))/((r[4]-r[2])*(r[3]-r[1]));prec=prec)))
        ψagm = 4/sqrt((r[4]-r[2])*(r[3]-r[1])) * Acb(π;prec=prec)/(2*agm)
        rec("gauntlet_d_quad_vs_K_overlap_$i", Arblib.overlaps(Pq, ψ1))
        rec("gauntlet_d_K_vs_AGM_overlap_$i", Arblib.overlaps(ψ1, ψagm))
        rec("gauntlet_d_certified_digits_$i", certified_digits(real(Pq)))
    end
end

# ---- (6) puncture layer values (dictionary entries), 500 digits ---------------
let prec = 1856
    s = Acb(2; prec); t = Acb(-1//2; prec)
    t0 = time(); Z = E.abel_image_deformed(s, t; delta_exp = 1700); tZ = time() - t0
    rec("Z_xp_seconds_D500", round(tZ, digits = 2))
    rec("Z_xp_certified_digits_D500", certified_digits(Z))
    rec("Z_xp_re_D500", string(real(Z)))
    rec("Z_xp_im_D500", string(imag(Z)))
    w = E.varpi0(s, Acb(1; prec))
    rec("varpi0_s2_x1_D500", string(real(w)))
    rec("varpi0_certified_digits", certified_digits(real(w)))
end

# ---- (7) cusp dictionary at 500 digits ----------------------------------------
let prec = 1856
    t0 = time()
    dict = cusp_dictionary(prec; maxweight = 8)
    rec("cusp_dictionary_seconds_D500", round(time() - t0, digits = 2))
    E.write_dictionary(joinpath(@__DIR__, "cusp-dictionary-500.txt"), dict; digits = 540)
    rec("cusp_dictionary_entries", length(dict))
end

close(out)
println("PRODUCTION RUN COMPLETE")
