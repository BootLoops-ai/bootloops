# Route B — INDEPENDENT oracle for c_alpha via Julia/Arb ball arithmetic, in
# z = -1/s, large-CW-circle monodromy at z=infty (= s=0 threshold).  GENERIC ORDER.
#
# Differs from Route A in: arithmetic (Arb balls), stepper (Eichler._theta_transport),
# coordinate (z), monodromy contour (large circle around ALL finite sings),
# seed (operator recursion, not multinomial).

# Eichler.jl project (ships in this repository under upgrades/Eichler.jl;
# EICHLER_PROJECT overrides — point it at your instantiated project).
import Pkg; Pkg.activate(get(ENV, "EICHLER_PROJECT",
    normpath(joinpath(@__DIR__, "..", "..", "upgrades", "Eichler.jl"))); io=devnull)
using Eichler; const E = Eichler
using Arblib

const HERE = @__DIR__
const DATAFILE = get(ENV, "PFDATA", "pf_cy3_data.jl")
include(joinpath(HERE, DATAFILE))
const PREC = parse(Int, get(ENV, "PREC", "512"))
const SDEC_STR = get(ENV, "SDEC", "1/4")
const SDEC_R = let s = split(SDEC_STR,'/'); length(s)==2 ? big(parse(Int,s[1]))//big(parse(Int,s[2])) : big(parse(Int,s[1]))//big(1); end
const NSER = parse(Int, get(ENV, "NSER", "2000"))
const NTHR = parse(Int, get(ENV, "NTHR", "200"))
const NLOOP = parse(Int, get(ENV, "NLOOP", "48"))
const OUT = get(ENV, "OUT", "ROUTEB.json")

println("Route B: Julia/Arb monodromy-projection in z=-1/s, prec=$PREC bits (~$(round(PREC*0.301)) d), order=$ORDER")
println("  data: $DATAFILE  msq=$MSQ  frac_exps=$FRAC_EXPS  sdec=$SDEC_R")

# ---- PFOperator from d/dz coeffs ----
coeffs = [[c//big(1) for c in DDZ_COEFFS[k+1]] for k in 0:ORDER]
L = E.PFOperator(PREC, coeffs)
println("  PFOperator: order=$(L.order), |sing|=$(length(L.sing))")
Rθ = E.theta_form(L)

# ---- INDEPENDENT seed: holomorphic ϖ_0 via OPERATOR RECURSION at z=0.
#      MUM indicial roots include 0 (mult>=order-?) and 1; resonance at n=1 -> a_1 free.
#      Pin a_0=1, a_1=-(Σ m_i^2) (BFKNS normalisation), recurse n>=2 with R0(n)!=0.
function frob_hol_series(Rθ, N::Int; a0=big(1)//1, a1=big(-SUM_MSQ)//1)
    S = length(Rθ) - 1
    Rpoly(s, x) = sum(Rθ[s+1][k+1] * x^k for k in 0:length(Rθ[s+1])-1)
    a = Vector{Rational{BigInt}}(undef, N+1)
    a[1] = a0
    R0 = n -> Rpoly(0, big(n))
    @assert R0(0) == 0 && R0(1) == 0  "indicial roots 0,1 expected"
    @assert Rpoly(1, big(0)) * a0 == 0  "n=1 resonance inconsistency"
    a[2] = a1
    for n in 2:N
        rhs = big(0)//1
        for s in 1:min(S, n)
            rhs -= Rpoly(s, big(n - s)) * a[n - s + 1]
        end
        r0n = R0(n)
        @assert r0n != 0  "unexpected resonance at n=$n"
        a[n+1] = rhs // r0n
    end
    return a
end
a = frob_hol_series(Rθ, NSER)
nck = min(length(BFKNS_A), NSER+1)
mismatch = count(i -> a[i] != BFKNS_A[i]//1, 1:nck)
println("  operator-recursion seed vs BFKNS multinomial (first $nck terms): mismatches=$mismatch")
@assert mismatch == 0

function theta_state(a, z::Acb, ord::Int)
    Y = [Acb(0; prec=PREC) for _ in 1:ord]
    zn = Acb(1; prec=PREC)
    for n in 0:length(a)-1
        an = Acb(a[n+1]; prec=PREC)
        t = an * zn
        nk = Acb(1; prec=PREC)
        for k in 1:ord
            Y[k] += nk * t
            nk *= n
        end
        zn *= z
    end
    Y
end

# Anchor: need |z_b| < min|z_sing| (convergence disk at MUM z=0).
# For CY3 (1,1,1,1,16) nearest sing |z|=1/64 -> s_b > 64; for K3 |z|=1/36 -> s_b > 36.
minsing = minimum(Float64(abs(Arblib.midpoint(Acb,s))) for s in L.sing if !Arblib.contains_zero(s))
sb_auto = parse(Int, get(ENV, "SB", string(ceil(Int, 1.6 / minsing))))
s_b = Acb(sb_auto; prec=PREC)
z_b = -1 / s_b
println("  min |z_sing| ~ $(round(minsing,sigdigits=4)) (=> s>$(round(1/minsing,digits=1))); anchor s_b=$sb_auto")
Yb = theta_state(a, z_b, ORDER)
println("  ϖ_0(-1/$sb_auto) rel_acc_bits=$(Arblib.rel_accuracy_bits(Yb[1]))")

# ---- transport z_b -> z_dec via UPPER-HALF z-plane ----
z_dec = Acb(-1//SDEC_R; prec=PREC)
sdec_arb = Arb(SDEC_R; prec=PREC)
lift = Acb(0, 1; prec=PREC) * Acb(1//2; prec=PREC)
path1 = Acb[z_b, z_b + lift, z_dec + lift, z_dec]
Yd, _ = E._theta_transport(L, path1, Yb)
println("  transported to z_dec=$(Float64(real(z_dec))): ϖ_0 rel_acc_bits=$(Arblib.rel_accuracy_bits(Yd[1]))")

# ---- monodromy at z=infty: CW large circle |z|=|z_dec|=1/sdec ----
twopi = 2 * Arb(pi; prec=PREC)
loop = Vector{Acb}(undef, NLOOP+1)
loop[1] = z_dec
for k in 1:NLOOP-1
    ph = Acb(0, -1; prec=PREC) * Acb(twopi * k / NLOOP; prec=PREC)
    loop[k+1] = z_dec * exp(ph)
end
loop[NLOOP+1] = z_dec

M = Matrix{Acb}(undef, ORDER, ORDER)
for col in 1:ORDER
    e = [Acb(k==col ? 1 : 0; prec=PREC) for k in 1:ORDER]
    Ye, _ = E._theta_transport(L, loop, e)
    for row in 1:ORDER; M[row, col] = Ye[row]; end
end
function mmul(A,B,n)
    C = [Acb(0; prec=PREC) for _ in 1:n, _ in 1:n]
    for i in 1:n, j in 1:n, k in 1:n; C[i,j] += A[i,k]*B[k,j]; end
    C
end
function mvmul(A,v,n)
    [sum(A[i,k]*v[k] for k in 1:n) for i in 1:n]
end
# det via LU-ish expansion (small n)
function det_acb(M,n)
    n==1 && return M[1,1]
    s = Acb(0; prec=PREC)
    for j in 1:n
        sub = [M[i,k] for i in 2:n, k in [c for c in 1:n if c!=j]]
        s += (-1)^(1+j) * M[1,j] * det_acb(sub, n-1)
    end
    s
end
detM = det_acb(M, ORDER)
println("  det(M_infty) = $detM")

# ---- spectral projectors ----
Iord = [Acb(i==j ? 1 : 0; prec=PREC) for i in 1:ORDER, j in 1:ORDER]
nilp = ORDER - length(FRAC_EXPS)
MmI = M .- Iord
MmIn = Iord
for _ in 1:nilp; global MmIn = mmul(MmIn, MmI, ORDER); end

eigs = [Acb(cospi(2*Arb(α;prec=PREC)), sinpi(2*Arb(α;prec=PREC)); prec=PREC) for α in FRAC_EXPS]

# ---- Phi_alpha theta_z-state on LOWER branch arg(t)=-pi at t=-sdec.
#      t^{alpha+n} on lower branch = exp(-i*pi*(alpha+n)) * sdec^{alpha+n}.
#      theta_t^k -> *(alpha+n)^k; theta_z^k = (-1)^k theta_t^k.
function phi_theta_state_lower(anum, aden, alpha::Rational, sdec::Arb, N::Int, ord::Int)
    Y = [Acb(0; prec=PREC) for _ in 1:ord]
    pi_a = Arb(pi; prec=PREC)
    for n in 0:min(N, length(anum)-1)
        an = Arb(anum[n+1]; prec=PREC) / Arb(aden[n+1]; prec=PREC)
        ex = Arb(alpha + n; prec=PREC)
        # t^ex on lower branch:
        tphase = Acb(cos(-pi_a*ex), sin(-pi_a*ex); prec=PREC)
        term = an * sdec^ex * tphase
        for k in 0:ord-1
            Y[k+1] += (k==0 ? Acb(1;prec=PREC) : Acb(ex;prec=PREC)^k) * term
        end
    end
    [Acb((-1)^k; prec=PREC) * Y[k+1] for k in 0:ord-1]
end

ANUM = [A0_NUM]; ADEN = [A0_DEN]
length(FRAC_EXPS) >= 2 && (push!(ANUM, A1_NUM); push!(ADEN, A1_DEN))

function digits_agree(a::Acb, b::Acb)
    d = a - b
    Arblib.contains_zero(d) && return Float64(PREC*0.301)
    return -log10(Float64(abs(Arblib.midpoint(Acb, d))) / max(Float64(abs(Arblib.midpoint(Acb, b))),1e-300))
end

results = Dict{String,Any}()
results["prec_bits"] = PREC
results["s_dec"] = string(SDEC_R)
results["seed_mismatches"] = mismatch
results["det_M"] = string(detM)
results["nilp"] = nilp
results["frac_exps"] = string.(FRAC_EXPS)
calpha = Dict{String,Any}()
for (ai, alpha) in enumerate(FRAC_EXPS)
    lam = eigs[ai]
    P = copy(MmIn)
    denom = (lam - 1)^nilp
    for (bi, beta) in enumerate(FRAC_EXPS)
        bi == ai && continue
        mu = eigs[bi]
        P = mmul(P, M .- mu .* Iord, ORDER)
        denom *= (lam - mu)
    end
    PY = [x / denom for x in mvmul(P, Yd, ORDER)]
    phi = phi_theta_state_lower(ANUM[ai], ADEN[ai], alpha, sdec_arb, NTHR, ORDER)
    crows = [PY[k] / phi[k] for k in 1:ORDER]
    rc = minimum(digits_agree(crows[1], crows[k]) for k in 2:ORDER)
    c = crows[1]
    println("\n  c_{$alpha} per row:")
    for k in 1:ORDER
        println("    row $k: $(crows[k])")
    end
    println("  row-consistency: $(round(rc,digits=1)) d")
    calpha[string(alpha)] = Dict(
        "c_re" => string(real(c)), "c_im" => string(imag(c)),
        "c_re_mid" => string(Arblib.midpoint(Arb,real(c))),
        "c_im_mid" => string(Arblib.midpoint(Arb,imag(c))),
        "c_acc_bits" => Arblib.rel_accuracy_bits(c),
        "row_consistency_d" => round(rc,digits=1),
    )
end
results["c_alpha"] = calpha

jesc(s) = replace(string(s), "\\"=>"\\\\", "\""=>"\\\"")
open(joinpath(HERE, OUT), "w") do f
    print(f, "{\n")
    pairs = collect(results)
    for (idx,(k,v)) in enumerate(pairs)
        sep = idx < length(pairs) ? "," : ""
        if k == "c_alpha"
            print(f, " \"c_alpha\": {\n")
            cap = collect(v)
            for (i2,(ka,va)) in enumerate(cap)
                print(f, "  \"$(jesc(ka))\": {")
                print(f, join(["\"$kk\": $(kk in ("c_acc_bits","row_consistency_d") ? vv : "\"$(jesc(vv))\"")" for (kk,vv) in va], ", "))
                print(f, "}$(i2<length(cap) ? "," : "")\n")
            end
            print(f, " }$sep\n")
        else
            print(f, " \"$k\": \"$(jesc(v))\"$sep\n")
        end
    end
    print(f, "}\n")
end
println("\n  wrote $OUT")
