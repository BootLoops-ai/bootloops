# test_cuspvals.jl — acceptance tests A–G for src/cuspvals.jl (exact cusp
# values of the Γ₁(6) Eisenstein bases of src/modular.jl).
#
# Run:  julia --project=<pkgroot> test/test_cuspvals.jl
#
# NOTE: all numerics here are self-contained (Nemo AcbField horocycles,
# Arblib acb_modular_eta); we deliberately do not call modular.jl's
# form_value/cusp_constant numeric layer, so the check stays independent of it.

using LandauAlphabet
using Nemo
using Arblib

const PKGROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(PKGROOT, "src", "cuspvals.jl"))

const RESULTS = Pair{String,Bool}[]

function record!(name::String, ok::Bool, detail::String = "")
    push!(RESULTS, name => ok)
    println(rpad(ok ? "PASS" : "FAIL", 6), "[", name, "]  ", detail)
    ok
end

subcheck(ok::Bool, what::String) = (println("    ", ok ? "ok  " : "BAD ", what); ok)

function print_K_matrix(title::String, M, colnames::Vector{String})
    println("\n", title)
    rownames = ["(1,0)=inf", "(0,1)", "(1,2)", "(1,3)"]
    w = max(12, maximum(length.(colnames)) + 2)
    println(" "^12, join(rpad.(colnames, w)))
    for i in 1:size(M, 1)
        println(rpad(rownames[i], 12), join(rpad.(string.(M[i, :]), w)))
    end
end

# ---------------------------------------------------------------------------
# shared exact data
# ---------------------------------------------------------------------------

println("building exact data (sunrise_data, level-6 bases, cusp-value matrices) ...")
const sd = sunrise_data(140)
const names2, basis2 = level6_basis(sd, 2)
const names3, basis3 = level6_basis(sd, 3)
const names4, basis4 = level6_basis(sd, 4)
const M2 = cusp_values_level6(2)
const M3 = cusp_values_level6(3)
const M4 = cusp_values_level6(4)
const KK = cuspvals_field()[1]

print_K_matrix("cusp values, weight 2 (rows = cusps, widths [1,6,3,2]):", M2, names2)
print_K_matrix("cusp values, weight 3:", M3, names3)
print_K_matrix("cusp values, weight 4:", M4, names4)
println()

# ---------------------------------------------------------------------------
# numeric helpers (Nemo AcbField; exact-rational inputs, rigorous-enough
# truncation bookkeeping documented inline)
# ---------------------------------------------------------------------------

"embed x ∈ K = Q(√−3) into AcbField: s3 ↦ i√3"
function embed_K(x, CC::AcbField, RR::ArbField)
    CC(RR(coeff(x, 0))) + CC(RR(coeff(x, 1))) * onei(CC) * CC(sqrt(RR(3)))
end

"""
    horocycle_avg(series, k, cusp, width; y, Npts, prec) -> (Acb, nterms)

(1/N) Σ_j (cτ_j+d)^{−k} F(γτ_j), τ_j = j·width/Npts + iy — numeric n=0
Fourier coefficient of F|_kγ.  F is summed from its exact q-series; the
truncation length is chosen from the worst Im(γτ_j) on the horocycle:
Im(γτ) = y/|cτ+d|² ≥ y/(max_x (cx+d)² + c²y²) =: Im_min, so the q-decay rate
is rate = 2π·Im_min per term and we keep nterms ≈ 115/rate terms
(tail ≲ e^{−115}·poly ≪ 1e−40).  Aliasing error ~ e^{−2πyN/width} ≪ 1e−100.
"""
function horocycle_avg(series::Vector{QQFieldElem}, k::Int, cusp::Tuple{Int,Int},
                       width::Int; y::Rational{Int} = 5//2, Npts::Int = 96,
                       prec::Int = 256)
    a, c = cusp
    b, d = cusp_completion(cusp)
    CC = AcbField(prec)
    RR = ArbField(prec)
    dmax = max(Float64(d)^2, (Float64(c) * width + d)^2)       # max of (cx+d)² on [0,width]
    im_min = Float64(y) / (dmax + Float64(c)^2 * Float64(y)^2)  # worst Im(γτ_j)
    rate = 2 * pi * im_min
    nterms = clamp(ceil(Int, 115 / rate), 50, length(series))
    nterms < length(series) || error("series too short for horocycle: need ~$(ceil(Int,115/rate)), have $(length(series))")
    coeffs = [CC(series[n]) for n in 1:nterms]
    twopii = CC(2 * const_pi(RR)) * onei(CC)
    acc = zero(CC)
    for j in 0:Npts-1
        tau = CC(RR(QQ(j * width, Npts)), RR(QQ(y)))
        den = c * tau + d
        q = exp(twopii * (a * tau + b) / den)
        F = zero(CC)
        for n in nterms:-1:1
            F = F * q + coeffs[n]
        end
        acc += F / den^k
    end
    (acc / Npts, nterms)
end

"t(z) = 9 η(z)⁴η(6z)⁸/(η(2z)⁸η(3z)⁴) via Arblib's rigorous acb_modular_eta"
function t_via_eta(z::Acb, prec::Int)
    eta = w -> (r = Acb(0, prec = prec); Arblib.modular_eta!(r, w); r)
    9 * eta(z)^4 * eta(6 * z)^8 / (eta(2 * z)^8 * eta(3 * z)^4)
end

# ===========================================================================
# A. ∞-cusp column == exact a₀ coefficients from qseries.jl
# ===========================================================================
let ok = true
    for (Mk, basis, lbl) in ((M2, basis2, "wt2"), (M3, basis3, "wt3"), (M4, basis4, "wt4"))
        agree = all(Mk[1, j] == KK(basis[j][1]) for j in 1:length(basis))
        ok &= subcheck(agree, "$lbl: row(∞) == [series a₀] from level6_basis")
    end
    # dilation-independent literals
    ok &= subcheck(collect(M2[1, :]) == [KK(QQ(K - 1, 24)) for K in (2, 3, 6)],
                   "B_{2,K} at ∞ == (K−1)/24, K = 2,3,6")
    ok &= subcheck(collect(M3[1, :]) == [zero(KK), zero(KK), KK(QQ(-1, 9)), KK(QQ(-1, 9))],
                   "E₃(Kτ;χ₋₃,1) at ∞ == 0 and E₃(Kτ;1,χ₋₃) at ∞ == −B_{3,χ₋₃}/6 = −1/9")
    ok &= subcheck(KK(-LandauAlphabet.gen_bernoulli(3, -3) // QQ(6)) == M3[1, 3],
                   "−B_{3,χ₋₃}/6 from gen_bernoulli matches")
    ok &= subcheck(collect(M4[1, :]) == [KK(QQ(1, 240)), KK(QQ(1, 240)), KK(QQ(1, 240)), KK(QQ(1, 240)), zero(KK)],
                   "E₄(Kτ) at ∞ == 1/240 for all K; cusp form column 0")
    record!("A", ok, "∞-cusp column equals exact a₀ of the q-expansions")
end

# ===========================================================================
# B. weight-2 residue sum rule Σ_cusps width·a₀ = 0  (+ numeric horocycles)
# ===========================================================================
let ok = true
    for K in (2, 3, 6)
        s = sum(QQ(GAMMA16_WIDTHS[i]) * cusp_value_B2K(K, GAMMA16_CUSP_TUPLES[i]) for i in 1:4)
        ok &= subcheck(iszero(s), "Σ width·a₀(B_{2,$K}|γ) == 0 (exact)")
    end
    # independent numeric horocycle check of every weight-2 cusp value
    prec = 256
    CC = AcbField(prec); RR = ArbField(prec)
    tol = RR(QQ(1, 10)^20)
    worst = zero(RR)
    for (j, K) in enumerate((2, 3, 6))
        ser = B2K(K, 2000)
        for i in 1:4
            av, _ = horocycle_avg(ser, 2, GAMMA16_CUSP_TUPLES[i], GAMMA16_WIDTHS[i]; prec = prec)
            dif = abs(av - embed_K(M2[i, j], CC, RR))
            worst = max(worst, dif)
        end
    end
    ok &= subcheck(worst < tol, "numeric horocycle check of all 12 weight-2 cusp values, worst |diff| = $(Float64(worst))")
    record!("B", ok, "weight-2 residue sum rule (exact) + numeric cross-check")
end

# ===========================================================================
# C. Landau dlog letters: integer residue matrix, zero column sums,
#    exact divisor of t (eta-quotient orders) + Arblib-eta numeric t at cusps
# ===========================================================================
let ok = true
    expected = Dict("dlog(t)" => [-8, -4, 8], "dlog(t-1)" => [-9, -3, 3], "dlog(t-9)" => [-5, 5, -1])
    letters = ["dlog(t)", "dlog(t-1)", "dlog(t-9)"]
    series = Dict("dlog(t)" => sd.g20, "dlog(t-1)" => sd.g21, "dlog(t-9)" => sd.g29)
    coords = Dict{String,Vector{QQFieldElem}}()
    for nm in letters
        v = solve_in_basis(basis2, qs_substitute_negq(series[nm]))
        okv = v !== nothing && v == QQ.(expected[nm])
        ok &= subcheck(okv, "level-6 coordinates of $nm == $(expected[nm]) in [B22,B23,B26]")
        coords[nm] = v === nothing ? QQFieldElem[] : v
    end
    # weight-2 matrix entries are rational (no s3 part)
    ok &= subcheck(all(iszero(coeff(M2[i, j], 1)) for i in 1:4, j in 1:3), "M2 entries rational")
    M2qq = [coeff(M2[i, j], 0) for i in 1:4, j in 1:3]
    R = [QQ(GAMMA16_WIDTHS[i]) * sum(coords[letters[l]][j] * M2qq[i, j] for j in 1:3)
         for i in 1:4, l in 1:3]
    println("    R[cusp, letter] = width·a₀(letter|cusp):")
    for i in 1:4
        println("      ", rpad(string(GAMMA16_CUSP_TUPLES[i]), 7), [R[i, l] for l in 1:3])
    end
    ok &= subcheck(all(isone(denominator(x)) for x in R), "R is an integer matrix")
    ok &= subcheck(all(iszero(sum(R[:, l])) for l in 1:3), "columns of R sum to 0")

    # exact cross-check of the divisor pattern via eta-quotient orders:
    # t in q_C is the integer-argument eta quotient 9·η(z)⁴η(6z)⁸/(η(2z)⁸η(3z)⁴)
    # (this resolves the half-integer η(m·τ_C/2+..) arguments of the q₂-language
    # via η(z+1/2) = e^{iπ/24}η(2z)³/(η(z)η(4z)); we verify it EXACTLY in q-series)
    lead, ser = eta_quotient_q2([(1, 4), (6, 8), (2, -8), (3, -4)], sd.N)
    tC = qs_substitute_negq(sd.t)
    cand = qs_scal(9, qs_shift_power(ser, 1))
    ok &= subcheck(lead == 1 && tC[1:sd.N] == cand[1:sd.N],
                   "EXACT identity t(q_C) == 9·η(z)⁴η(6z)⁸/(η(2z)⁸η(3z)⁴)  ($(sd.N) coeffs)")
    tfac = [(1, 4), (2, -8), (3, -4), (6, 8)]
    ordt = [sum(QQ(r) * eta_order_at_cusp(m, GAMMA16_CUSP_TUPLES[i], GAMMA16_WIDTHS[i]) for (m, r) in tfac)
            for i in 1:4]
    ok &= subcheck(ordt == [QQ(1), QQ(0), QQ(-1), QQ(0)],
                   "exact divisor of t (eta orders) = (∞) − (1,2): t is a degree-1 hauptmodul")
    ok &= subcheck([R[i, 1] for i in 1:4] == [QQ(x) for x in (1, 0, -1, 0)],
                   "dlog(t) residue column == divisor of t")
    # since deg t = 1, the +1 rows of columns 2,3 PROVE t((0,1)) = 1, t((1,3)) = 9
    ok &= subcheck([R[i, 2] for i in 1:4] == [QQ(x) for x in (0, 1, -1, 0)] &&
                   [R[i, 3] for i in 1:4] == [QQ(x) for x in (0, 0, -1, 1)],
                   "residue pattern ⇒ t = 1 at cusp (0,1), t = 9 at cusp (1,3), t = ∞ at (1,2)")
    println("    fiber assignment (width w ⇒ Kodaira I_w): t=0 ↔ ∞-cusp (I₁), t=1 (pseudothreshold) ↔ (0,1) (I₆),")
    println("                                              t=∞ ↔ (1,2) (I₃), t=9 (threshold) ↔ (1,3) (I₂)")

    # numeric confirmation with Arblib's rigorous Dedekind eta (valid arbitrarily
    # close to the cusp: acb_modular_eta reduces to the fundamental domain)
    prec = 256
    y = 12.0
    tvals = Acb[]
    for i in 1:4
        a, c = GAMMA16_CUSP_TUPLES[i]
        b, d = cusp_completion((a, c))
        tau = Acb(Arb(0, prec = prec), Arb(y, prec = prec), prec = prec)
        z = (a * tau + b) / (c * tau + d)
        push!(tvals, t_via_eta(z, prec))
    end
    ok &= subcheck(abs(tvals[1]) < 1e-10, "numeric t near ∞-cusp ≈ 0      (|t| = $(Float64(abs(tvals[1]))))")
    ok &= subcheck(abs(tvals[2] - 1) < 1e-3, "numeric t near (0,1) ≈ 1      (|t−1| = $(Float64(abs(tvals[2]-1))))")
    ok &= subcheck(abs(tvals[3]) > 1e6, "numeric t near (1,2) → ∞      (|t| = $(Float64(abs(tvals[3]))))")
    ok &= subcheck(abs(tvals[4] - 9) < 1e-6, "numeric t near (1,3) ≈ 9      (|t−9| = $(Float64(abs(tvals[4]-9))))")
    record!("C", ok, "dlog letters: integer residues, zero column sums, t-divisor cross-checks")
end

# ===========================================================================
# D. weight-4 cross-check: lattice method == Hermite-decomposition method
# ===========================================================================
let ok = true
    for Kd in (1, 2, 3, 6), cu in GAMMA16_CUSP_TUPLES
        agree = cusp_value_E(4, 1, 1, Kd, cu) == KK(cusp_value_E4_hermite(Kd, cu))
        agree || subcheck(false, "mismatch at K=$Kd, cusp $cu")
        ok &= agree
    end
    ok && subcheck(true, "all 16 (dilation × cusp) E₄ values agree exactly")
    record!("D", ok, "E₄(Kτ): lattice sums == Hermite decomposition (exact, 4×4 cases)")
end

# ===========================================================================
# E. numeric horocycle averages for the weight-3 forms (25+ digits)
# ===========================================================================
let ok = true
    prec = 256
    CC = AcbField(prec); RR = ArbField(prec)
    tol = RR(QQ(1, 10)^20)
    specs = ((-3, 1, 1), (-3, 1, 2), (1, -3, 1), (1, -3, 2))   # (Dphi, Dpsi, Kdil), basis order
    worst = zero(RR)
    for (j, (Dphi, Dpsi, Kd)) in enumerate(specs)
        ser = eisenstein_E(3, Dphi, Dpsi, Kd, 2000)
        for i in 2:4                                            # cusps (0,1), (1,2), (1,3)
            cu = GAMMA16_CUSP_TUPLES[i]
            av, nt = horocycle_avg(ser, 3, cu, GAMMA16_WIDTHS[i]; y = 5//2, Npts = 96, prec = prec)
            dif = abs(av - embed_K(M3[i, j], CC, RR))
            agree = dif < tol
            agree || subcheck(false, "E₃($(Kd)τ;$Dphi,$Dpsi) at cusp $cu: |diff| = $(Float64(dif)) (nterms=$nt)")
            ok &= agree
            worst = max(worst, dif)
        end
    end
    ok &= subcheck(worst < tol, "all 12 weight-3 horocycle averages agree, worst |diff| = $(Float64(worst)) < 1e-20")
    record!("E", ok, "weight-3 numeric horocycle cross-check at 25+ digits (y=5/2, N=96)")
end

# ===========================================================================
# F. gate-(iii): f̂₃ coordinates (−1,8,0,0) recoverable from cusp data
# ===========================================================================
let ok = true
    v3 = solve_in_basis(basis3, qs_substitute_negq(sd.f3hat))
    cexp = [QQ(-1), QQ(8), QQ(0), QQ(0)]
    ok &= subcheck(v3 == cexp, "solve_in_basis(f̂₃) == (−1, 8, 0, 0)")
    M3n = matrix(KK, M3)
    ok &= subcheck(!iszero(det(M3n)), "weight-3 cusp-value matrix invertible over K (det = $(det(M3n)))")
    vF = M3n * matrix(KK, 4, 1, [KK(x) for x in cexp])          # cusp values of f̂₃
    println("    cusp values of f̂₃: ", [vF[i, 1] for i in 1:4])
    fl, sol = can_solve_with_solution(M3n, vF, side = :right)
    ok &= subcheck(fl && [sol[i, 1] for i in 1:4] == [KK(x) for x in cexp],
                   "M·c = v(f̂₃) solves exactly to (−1, 8, 0, 0)")
    record!("F", ok, "weight-3 cusp matrix invertible; f̂₃ recovered from cusp data")
end

# ===========================================================================
# G. weight-4 rank: rank 4, kernel = cusp-form direction (0,0,0,0,1)
# ===========================================================================
let ok = true
    M4n = matrix(KK, M4)
    ok &= subcheck(rank(M4n) == 4, "rank of the 4×5 weight-4 cusp-value matrix == 4")
    nd, ns = nullspace(M4n)
    ok &= subcheck(nd == 1, "right kernel is 1-dimensional")
    w = [ns[i, 1] for i in 1:5]
    ok &= subcheck(all(iszero, w[1:4]) && !iszero(w[5]),
                   "kernel = span{(0,0,0,0,1)} — the 6.4.a.a direction is cusp-data-blind")
    record!("G", ok, "weight-4 cusp-value matrix: rank 4, kernel = cusp-form direction")
end

# ---------------------------------------------------------------------------
println("\n================ SUMMARY ================")
allok = true
for (name, ok) in RESULTS
    println(rpad(name, 4), ok ? "PASS" : "FAIL")
    global allok &= ok
end
println("overall: ", allok ? "ALL PASS" : "FAILURES PRESENT")
exit(allok ? 0 : 1)
