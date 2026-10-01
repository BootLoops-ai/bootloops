# Shared definitions for the massless planar double-box oracle.
# Conventions are pinned to
# Smirnov hep-ph/9905323 eq. (2box):
#
#   ∫∫ d^dk d^dl / [(k+p1)^2 (k-p2)^2 k^2 (k-l)^2 (l+p1)^2 (l-p2)^2 (l-p1-p3)^2]
#     = (i π^{d/2} e^{-γ_E ε})^2 K(t/s, ε) / ((-s)^{2+2ε} (-t)),
#   s = (p1+p2)^2,  t = (p1+p3)^2,  x = t/s,  d = 4 - 2ε.
#
# amflow-cpp measure is d^dℓ/(iπ^{d/2}), so
#   J_amflow = e^{-2γ_E ε} K(x,ε) / ((-s)^{2+2ε} (-t)),
#   i.e.  K(x,ε) = e^{2γ_E ε} (-s)^{2+2ε} (-t) · J.
# Euclidean sampling at s = -1, t = -x (x>0):    K_j = x · [e^{2γε} J]_j.
# Physical  sampling at s = +1, t = -x (0<x<1), s+i0:
#   ln(-s-i0) = -iπ  ⇒  (-s)^{2+2ε} = e^{-2πiε},  K_j = x·[e^{(2γ-2πi)ε} J]_j.
#
# Known answer (Smirnov):  K = Σ_{j=-4}^{0} ε^j K^{(j)}(x) + O(ε),
# uniform weight 4+j at order ε^j, HPL alphabet {x, 1+x} (Henn 1304.1806).

using LandauBootstrap, Arblib, Nemo, JSON3
using LandauBootstrap: magof, ub

const DBOX_ROOT = "LandauBootstrap"
const DBOX_SAMPLES = joinpath(DBOX_ROOT, "results", "samples")

# ---------------------------------------------------------------- family ----
# Momentum routing and amf options taken verbatim from the amflow-cpp bench
# doublebox_sv (tools/bench/doublebox_sv_eps001_*.json), which PASSED the
# Tier-1 pristine audit in this environment.  The corner integral
# I_{1111111;00} is Smirnov's eq. (2box) integral under momentum shifts
# (k = -l1 - p1, l = -l2 - p1 there): same scalar double box, unit numerator,
# s = (p1+p2)^2, t = (p1+p3)^2.  Our own Smirnov routing hit the open
# audit-failure-#7 crash ("entry still has a pole at eta = 0 after the
# shift") with both default and MASSIVE_FAMILY_OPTS — recorded in PROGRESS.
const DBOX_FAMILY = AmflowFamily(
    "doublebox", ["l1", "l2"], ["p1", "p2", "p3", "p4"],
    Dict("p4" => "-p1 - p2 - p3"),
    Dict("p1^2" => "0", "p2^2" => "0", "p3^2" => "0", "p4^2" => "0",
         "(p1 + p2)^2" => "s", "(p1 + p3)^2" => "t"),
    ["l1^2", "(l1 + p1)^2", "(l1 + p1 + p2)^2", "(l2 + p1 + p2)^2",
     "(l2 + p1 + p2 + p4)^2", "l2^2", "(l1 - l2)^2",
     "(l1 + p3)^2", "(l2 + p1)^2"])
const DBOX_INDICES = [1, 1, 1, 1, 1, 1, 1, 0, 0]
const DBOX_AMF_OPTS = (
    top = Dict{String,Any}("ending_schemes" => ["Tradition", "SingleMass"],
                           "max_recursion_depth" => 16),
    blackbox = Dict{String,Any}("ibp_rank" => 0, "ibp_dot" => 1,
                                "integral_order" => 5))

# ---------------------------------------------------------------- points ----
# Euclidean: s=-1, t=-x, x>0 (K real there).  Last two = reference pair.
const XS_EUCLID = [1//7, 1//5, 1//4, 1//3, 1//2, 4//7, 3//5, 2//3, 5//6, 8//9,
                   2//5, 3//4]
const N_HOLDOUT = 2
# Physical s-channel: s=+1, t=-x, 0<x<1 (so u=-s-t<0); +i0 on s.
const XS_PHYS = [1//5, 1//3, 1//2, 2//3]

sample_tag(x::Rational, phys::Bool) =
    (phys ? "dbox_phys_x" : "dbox_eucl_x") * "$(numerator(x))_$(denominator(x))"

# ------------------------------------------------------- ε-series helpers ---
# series = Dict{Int,Acb} (order -> coefficient), truncated to orders ≤ jmax.
function series_mul(a::Dict{Int,Acb}, b::Dict{Int,Acb}, jmax::Int, prec::Int)
    out = Dict{Int,Acb}()
    for (ja, ca) in a, (jb, cb) in b
        j = ja + jb
        j > jmax && continue
        out[j] = get(out, j, Acb(0, prec = prec)) + ca * cb
    end
    out
end

"series of exp(c·ε) through ε^n"
function series_exp(c::Acb, n::Int, prec::Int)
    out = Dict{Int,Acb}()
    term = Acb(1, prec = prec)
    out[0] = term
    for k in 1:n
        term = term * c / k
        out[k] = term
    end
    out
end

"""
    K_from_J(J::Dict{Int,Acb}, x::Rational, phys::Bool, prec) -> Dict{Int,Acb}

K(x,ε) coefficients (orders -4..0) from an amflow J sample at
(s,t) = (-1,-x) Euclidean or (+1,-x) physical (+i0 on s).
"""
function K_from_J(J::Dict{Int,Acb}, x::Rational, phys::Bool, prec::Int)
    gamma = Arblib.const_euler!(Arb(prec = prec))
    c = Acb(2 * gamma, prec = prec)
    if phys
        # K multiplies by (-s)^{2+2ε}; at s = 1+i0: ln(-s-i0) = -iπ, so
        # (-s)^{2+2ε} = e^{-2πiε}
        c -= Acb(0, 2 * Arb(π, prec = prec), prec = prec)
    end
    pref = series_exp(c, 8, prec)
    KJ = series_mul(Dict(k => v for (k, v) in J), pref, 0, prec)
    Dict(j => Acb(x, prec = prec) * v for (j, v) in KJ if -4 <= j <= 0)
end

# ------------------------------------------------- Smirnov closed form ------
# Independent evaluation stack: Arb log/polylog (acb_polylog) + direct
# Nielsen series  S_{1,2}(z) = Σ_{n≥2} H_{n-1} z^n / n²,
#                 S_{2,2}(z) = Σ_{n≥2} H_{n-1} z^n / n³   (|z| < 1),
# with an attached geometric tail bound.  NO engine GPL code is used here.
function nielsen_series(p::Int, z::Acb, prec::Int)   # S_{1,p? } no: S_{q,2}, q=p
    # returns S_{p,2}(z) = Σ H_{n-1} z^n / n^{p+1}
    azu = Arblib.ubound(Arf, abs(z))
    az = Arb(azu, prec = 64)
    Float64(az) < 0.95 || error("nielsen_series: |z| too close to 1")
    N = ceil(Int, (prec * log(2) + 20) / log(1 / Float64(az))) + 30
    s = Acb(0, prec = prec)
    H = Arb(0, prec = prec)
    zp = Acb(1, prec = prec)
    for n in 1:N
        zp *= z                       # z^n
        n >= 2 && (s += H * zp / Acb(n, prec = prec)^(p + 1))
        H += Arb(1, prec = prec) / n  # H_n, used at n+1
    end
    # tail: Σ_{n>N} H_{n-1}|z|^n/n^{p+1} ≤ H-ish·|z|^{N+1}/(1-|z|); be generous
    tail = Arb(log(N + 1.0) + 2, prec = 64) * az^(N + 1) / (1 - az)
    Arblib.add_error!(s, magof(tail))
    s
end

polylog_acb(k::Int, z::Acb, prec::Int) =
    Arblib.polylog!(Acb(prec = prec), Acb(k, prec = prec), z)

"""
    smirnov_K(x::Acb, prec) -> Dict{Int,Acb}

The hep-ph/9905323 K^{(j)}(x), j = -4..0.  For physical points pass x as a
negative real Acb: principal branch log(x) = ln|x| + iπ realizes t/s + i0.
"""
function smirnov_K(x::Acb, prec::Int)
    P  = Arb(π, prec = prec)
    P2 = Acb(P^2, prec = prec)
    P4 = Acb(P^4, prec = prec)
    z3 = Acb(Arblib.zeta!(Arb(prec = prec), Arb(3, prec = prec)), prec = prec)
    L  = log(x)
    Lp = log(1 + x)
    li2 = polylog_acb(2, -x, prec)
    li3 = polylog_acb(3, -x, prec)
    li4 = polylog_acb(4, -x, prec)
    s12 = nielsen_series(1, -x, prec)   # S_{1,2}(-x)
    s22 = nielsen_series(2, -x, prec)   # S_{2,2}(-x)
    K = Dict{Int,Acb}()
    K[-4] = Acb(-4, prec = prec)
    K[-3] = 5 * L
    K[-2] = -2 * L^2 + Acb(5//2, prec = prec) * P2
    K[-1] = -Acb(2//3, prec = prec) * L^3 - Acb(11//2, prec = prec) * P2 * L +
            Acb(65//3, prec = prec) * z3 +
            (-4 * li3 + 4 * L * li2 + 2 * (L^2 + P2) * Lp)
    K[0]  = Acb(4//3, prec = prec) * L^4 + 6 * P2 * L^2 -
            Acb(88//3, prec = prec) * z3 * L + Acb(29//30, prec = prec) * P4 +
            (-4 * (s22 - L * s12) + 44 * li4 - 4 * (Lp + 6 * L) * li3 +
             2 * (L^2 + 2 * L * Lp + Acb(10//3, prec = prec) * P2) * li2 +
             (L^2 + P2) * Lp^2 -
             Acb(2//3, prec = prec) * (4 * L^3 + 5 * P2 * L - 6 * z3) * Lp)
    K
end

# ------------------------------------------------------------ fit basis -----
# Uniform-weight basis at weight k over GPL indices {0,-1} (positions 1,2 in
# the GPLContext) = the N_sim "reality" tower: ALL words, ζ-graded, odd-π
# constants dropped (K_j real for x > 0).
#
# RECORDED FINDING (provisional fits + closed-form witness): the
# first-entry condition does NOT cut the per-ε-order basis of K.  It holds
# for I (2-var words start with s or t only), but K = (-s)^{2+2ε}(-t)·I and
# the ε-dependent prefactor shuffles log(-s) letters across ε orders: the
# (s+t)-leading projections of a (1+x)-first x-word in K_j cancel against
# log^m(-s)⧢K_{j-m} cross terms, not within K_j.  Witness: expanding
# Smirnov's K^(-1) gives pure part -4G(0,0,0) + 4G(-1,0,0) — a (-1)-first
# pure word.  Fits with first-letter-x bases (pure sector or everywhere)
# FAIL at ε^-1/ε^0 with residual ~1e-6 while ε^-4..ε^-2 pass exactly.
#   B_k = ⋃_{ζ ∈ {1,ζ2,ζ3,ζ4}, wt(ζ)≤k} ζ · { G(a⃗; x) : a⃗ ∈ {0,-1}^{k-wt(ζ)} }
function dbox_cut_basis(k::Int)
    allw(j) = j == 0 ? [Int[]] :
              vec([collect(t) for t in Iterators.product(fill((1, 2), j)...)])
    out = MPLBasisFunction[]
    for (csym, cw) in ((:one, 0), (:zeta2, 2), (:zeta3, 3), (:zeta4, 4))
        cw > k && continue
        for w in allw(k - cw)
            push!(out, MPLBasisFunction(csym, w))
        end
    end
    out
end

# raw / sim tower sizes for the (N_raw -> M_residual) record
dbox_Nraw(k::Int) = k == 0 ? 1 : length(basis_tower(2, k; variant = :unif))
dbox_Nsim(k::Int) = k == 0 ? 1 :
    length(filter(b -> !(b.csym in (:pi, :pi3, :pi5)),
                  basis_tower(2, k; variant = :unif)))

weight_of_order(j::Int) = 4 + j   # ε^j coefficient of K has uniform weight 4+j

# ------------------------------------------------------------- sample IO ----
function save_sample(tag::String, res::Dict{Int,Acb}, meta::Dict)
    mkpath(DBOX_SAMPLES)
    d = Dict{String,Any}(string(k) => Dict("re" => string(real(v)),
                                           "im" => string(imag(v)))
                         for (k, v) in res)
    d["_meta"] = meta
    open(joinpath(DBOX_SAMPLES, tag * ".json"), "w") do io
        JSON3.write(io, d)
    end
end

function load_sample(tag::String, prec::Int)
    f = joinpath(DBOX_SAMPLES, tag * ".json")
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
