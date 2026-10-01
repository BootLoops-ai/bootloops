# Production part 2: true 500-digit deliverables + Zenodo threshold cross-checks.
using Arblib, Eichler, Nemo
const E = Eichler
out = open(joinpath(@__DIR__, "results-production2.txt"), "w")
rec(k, v) = (println(out, k, " = ", v); flush(out); println(k, " = ", v))
cd(v) = floor(Int, max(0, Arblib.rel_accuracy_bits(v isa Acb ? real(v) : v)) * 0.30103)

# 500-digit sunrise with full headroom
prec = 2496; N = 1450
t0 = time(); G16 = Gamma16(prec, N); rec("gamma16_build_seconds_D500hr", round(time()-t0, digits=2))
t0 = time(); B = sunrise_boundary_constants(prec, 12); rec("boundary_seconds_D500hr", round(time()-t0, digits=2))
rec("B0_certified_digits", cd(B[1]))
L2 = E.dirichlet_L_chim3(2, prec)
rec("gauntlet_a_overlap_500", Arblib.overlaps(B[1], Acb(3*sqrt(Arb(3;prec=prec))/2*L2; prec=prec)))
t0 = time(); S4 = sunrise4(G16, Acb(-1; prec); J = 10); rec("sunrise4_seconds_D500hr", round(time()-t0, digits=2))
for k in -2:7
    rec("sunrise4_eps$(k)_certified_digits", cd(E.eps_coeff(S4, k)))
end
rec("sunrise4_eps0_value_540dig", E.decimal_midpoint(real(-E.eps_coeff(S4,0)), 540))
rec("B0_value_540dig", E.decimal_midpoint(real(B[1]), 540))

# Zenodo threshold cross-checks (s = 4m², the G3-report interior dictionary point)
let p2 = 512
    w0 = E.varpi0(Acb(4;prec=p2), Acb(1;prec=p2))
    a0 = Arb("0.47250316546487902700221989863146743636324935057104918758306141468051"; prec=p2)
    rec("zenodo_a0_overlap", Arblib.overlaps(real(w0), a0))
    rec("zenodo_a0_absdiff", Float64(Arblib.ubound(Arb, abs(real(w0)-a0))))
    w0p = E.varpi0_prime(Acb(4;prec=p2), Acb(1;prec=p2))
    dsw0 = -w0/4 - w0p/4         # Zenodo: ds w0 = -w0/s - (m2/s) dm2 w0
    a1 = Arb("0.06513382235683734253462525034065980272226207906515132998094255579545"; prec=p2)
    rec("zenodo_a1_overlap_minus_dsw0", Arblib.overlaps(real(-dsw0), a1))
end

# G third-kind production values + interior dictionary at the threshold
let p2 = 1856
    t0 = time()
    Gv = third_kind_G(Acb(2; prec=p2), Acb(-1//2; prec=p2))
    rec("G_s2_tm05_seconds", round(time()-t0, digits=2))
    rec("G_s2_tm05_certified_digits", cd(Gv))
    rec("G_s2_tm05_value", E.decimal_midpoint(real(Gv), 500))
    dict = interior_dictionary(Arb(2; prec=p2), Arb(-1//2; prec=p2), p2)
    E.write_dictionary(joinpath(@__DIR__, "interior-dictionary-s2-500.txt"), dict; digits = 520)
    rec("interior_dictionary_entries", length(dict))
    # threshold point s = 4m² (curve non-degenerate): periods + derivative ring
    w0 = E.varpi0(Acb(4;prec=p2), Acb(1;prec=p2))
    w0p = E.varpi0_prime(Acb(4;prec=p2), Acb(1;prec=p2))
    rec("threshold_varpi0_500dig", E.decimal_midpoint(real(w0), 500))
    rec("threshold_dvarpi0_500dig", E.decimal_midpoint(real(w0p), 500))
end
close(out)
println("PRODUCTION2 COMPLETE")
