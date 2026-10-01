# types.jl — SHARED TYPES. Every module codes against these definitions as
# they are; a change here is a package-level design decision that touches
# every module, never a local edit.
# Semantics pinned in CONTRACTS.md; corrections RT-i..RT-x in DESIGN.md.
#
# This file is include()d inside `module SubTropica` after `using Nemo`,
# `using JSON` — do not add `using` here.

# ---------------------------------------------------------------------------
# ε-linear exponent: a + b·ε, exact rationals.
# Used for poly exponents in EulerIntegrand/LogIntegrand (incl. the
# ε-DEPENDENT jacobian exponents (1-u)^(1-TropI(ρ)) — RT-extra-jac) and for
# Γ-arguments in Prefactor (loggamma_series input — RT-vii).
# ---------------------------------------------------------------------------
struct EpsExp
    a::Rational{BigInt}   # constant part
    b::Rational{BigInt}   # coefficient of ε
end
EpsExp(a::Union{Integer,Rational}) = EpsExp(Rational{BigInt}(a), 0//1)
Base.:(==)(x::EpsExp, y::EpsExp) = x.a == y.a && x.b == y.b
Base.hash(x::EpsExp, h::UInt) = hash((x.a, x.b), hash(:EpsExp, h))
iseps(x::EpsExp) = !iszero(x.b)

# ---------------------------------------------------------------------------
# Prefactor: c · ε^k · ∏_i Γ(g_i)^{p_i} · e^{n·γ_E·ε}
# (rational constant, explicit overall ε power, Γ-factors with EpsExp
#  arguments and integer powers — negative = denominator, e.g. 1/∏Γ(ν_i) —
#  and the e^(L·γ_E·ε)-type normalization). Symbolic ε-expansion of the
#  Γ-factors is C8's job via loggamma_series (RT-vii).
# ---------------------------------------------------------------------------
struct Prefactor
    c::Rational{BigInt}                    # rational multiplier
    eps_power::Int                         # overall ε^k
    gammas::Vector{Tuple{EpsExp,Int}}      # (argument, power); power<0 = denominator
    gammaE_eps::Rational{BigInt}           # n in e^{n·γ_E·ε} normalization
end
Prefactor(c::Union{Integer,Rational}=1//1) =
    Prefactor(Rational{BigInt}(c), 0, Tuple{EpsExp,Int}[], 0//1)

# ---------------------------------------------------------------------------
# EulerIntegrand — THE Phase-A input: the raw Euler quadruple (RT-x).
#   integrand = prefactor · ∏_i x_i^(nu_i) · ∏_j P_j^(e_j),  e_j = a_j+b_j·ε
# over x_i ∈ [0,∞)^n, flat measure dx. NO graph front-end in Phase A;
# fixtures come from pySecDec LoopIntegralFromGraph glue (verify.jl).
# Polys live in ONE Nemo ring QQ[vars..., kinvars...] (variable order:
# integration vars first, in `vars` order, then kinvars in `kinvars` order).
# Region contract: all x_i and all kinematic coefficients POSITIVE
# (declared-Euclidean; region=:physical REFUSES in v1 — R4).
# ---------------------------------------------------------------------------
struct EulerIntegrand
    prefactor::Prefactor
    nu::Vector{EpsExp}                          # monomial exponents (may be ε-dependent)
                                                # MEASURE LAW:
                                                # `nu` is the FLAT-measure exponent (∏x^nu·dx — the Phase-A
                                                # quadruple/fixture boundary, frozen). The B1 tropical layer
                                                # (SPEC_B1 §0/§3.1: trop_on_ray/divergence_data/subtraction_terms)
                                                # reads nu in the DLOG chart, ν_dlog = ν_flat + 1; the bridge is
                                                # b1_driver.jl `_b1_dlog` — front-ends emit FLAT nu, never dlog.
    polys::Vector{Tuple{QQMPolyRingElem,EpsExp}} # (P_j, a_j + b_j·ε)
    vars::Vector{Symbol}                        # ORDERED integration variables
    kinvars::Vector{Symbol}                     # kinematic symbols (ring vars after `vars`)
    ring::QQMPolyRing                           # common parent of all P_j
end

# ---------------------------------------------------------------------------
# HlogExpr — closed-form symbolic values (HF coef strings parsed, our own
# γ_E/ζ series coefficients, assembled Laurent coefficients).
# Structure: sum of HTerm; each HTerm = exact rational-function coefficient
# (CONTENT-CANONICAL Mma-grammar string over kinvars — CONTRACTS.md §d)
# × ∏ atom^power. Deterministic canonical text form defined in §d.
# ---------------------------------------------------------------------------
abstract type HAtom end

struct HConst <: HAtom       # named constants: :Pi, :I, :EulerGamma, :Log2, :Catalan
    name::Symbol
end
struct HZeta <: HAtom        # mzv_a_b_c → ζ(a,b,c); NEGATIVE entry = alternating ('m' prefix)
    idx::Vector{Int}
end
struct HLogA <: HAtom        # Log[arg]; arg = canonical poly/rational Mma-grammar string
    arg::String
end
struct HHlog <: HAtom        # Hlog[z, {letters...}]; z and letters = canonical strings
    z::String
    word::Vector{String}
end
struct HPeriod <: HAtom      # residual ZeroInfPeriod[word] key (empty word never stored);
    word::Vector{String}     # evaluated via zero_inf_period op + ginac ≥60dps gate (A(v))
end
struct HAlgLetter <: HAtom   # Wm_i / Wp_i / sqrt_disc_i; kind ∈ (:Wm, :Wp, :sqrt_disc);
    kind::Symbol             # idx indexes the AlgLetter table carried alongside
    idx::Int                 # (LaurentSeries.letters after union; per-response before)
end
struct HDelta <: HAtom       # delta[var] contour residue (I·Pi·delta — Euclidean guard, R4)
    var::Symbol
end

for T in (:HConst, :HZeta, :HLogA, :HHlog, :HPeriod, :HAlgLetter, :HDelta)
    @eval Base.:(==)(x::$T, y::$T) =
        all(getfield(x, f) == getfield(y, f) for f in fieldnames($T))
    @eval Base.hash(x::$T, h::UInt) =
        hash(Tuple(getfield(x, f) for f in fieldnames($T)), hash($(QuoteNode(T)), h))
end

struct HTerm
    coef::String                      # canonical rational-function string (CONTRACTS.md §d)
    atoms::Vector{Pair{HAtom,Int}}    # atom => nonzero integer power, canonically ordered
end

struct HlogExpr
    terms::Vector{HTerm}
    vars::Vector{Symbol}              # kinematic symbol context of the coef strings
end
HlogExpr() = HlogExpr(HTerm[], Symbol[])   # zero

# ---------------------------------------------------------------------------
# Algebraic-letter table entry (HF response `algebraic_letters` field):
# quadratic factor bookkeeping for Wm_i/Wp_i/sqrt_disc_i.
# ---------------------------------------------------------------------------
struct AlgLetter
    idx::Int
    poly::String      # the quadratic polynomial (canonical string)
    var::Symbol       # its variable
    lc::String        # leading coefficient
    sum::String       # root sum   (Vieta)
    product::String   # root product (Vieta)
    disc::String      # discriminant
end

# ---------------------------------------------------------------------------
# LaurentSeries{T} — the pipeline result: Σ_{k≥0} coeffs[k+1]·ε^(minorder+k).
# T = HlogExpr for symbolic results; T = BigFloat/Complex for the numeric
# side of the verify comparator (same shape both sides).
# provenance: :analytic | :numerically_matched_sign (R4 quarantine — the
# latter REQUIRES evidence: matching points + floored digits, with the
# sign-matching points disjoint from ≥1 held-out verification point).
# ---------------------------------------------------------------------------
struct LaurentSeries{T}
    minorder::Int
    coeffs::Vector{T}
    letters::Vector{AlgLetter}        # unioned algebraic-letter table
    provenance::Symbol
    evidence::Dict{String,Any}        # points/digits/oracles; run tag; pins
end
LaurentSeries(minorder::Int, coeffs::Vector{T}) where {T} =
    LaurentSeries{T}(minorder, coeffs, AlgLetter[], :analytic, Dict{String,Any}())

# ---------------------------------------------------------------------------
# LogIntegrand — the locally-finite handoff object (CONTRACTS.md §b):
# per face, per ε-order, per counter-term:
#   (num/den)(face vars, kinvars) · ∏_j Log[P_j]^{k_j},
# flat measure dx on [0,∞)^{|vars|} (the /∏x_J division ALREADY applied,
# vars ∉ J ALREADY gauged to 1 — J per CONTRACTS.md WARNING-2).
# Individually possibly log-divergent; finite only in the face sum
# (⇒ check_divergences OFF for these inside the farm, ON everywhere else).
# ---------------------------------------------------------------------------
struct LogIntegrand
    num::QQMPolyRingElem
    den::QQMPolyRingElem
    logs::Vector{Pair{QQMPolyRingElem,Int}}  # Log[P]^k, k ≥ 1
    vars::Vector{Symbol}                     # remaining (face) integration vars, ordered
    kinvars::Vector{Symbol}
end

# ---------------------------------------------------------------------------
# HFConfig — HyperFLINT invocation configuration (CONTRACTS.md §a).
# Defaults = the READ-ONLY installed engine; override via env or kwargs.
# EVERY invocation runs under `ulimit -v 32505856` (the convention used
# throughout this package; the bridge enforces it in the spawned shell command).
# ---------------------------------------------------------------------------
Base.@kwdef struct HFConfig
    bin::String = get(ENV, "SUBTROPICA_HF_BIN",
        "hyperflint.sh")  # public fallback: PATH lookup — set SUBTROPICA_HF_BIN to your built engine wrapper
    mzv_data_path::String = get(ENV, "SUBTROPICA_MZV_DATA",
        "mzv_reductions.json")  # public fallback: cwd-relative — set SUBTROPICA_MZV_DATA explicitly
    omp_num_threads::Int = 4
    max_threads_per_call::Int = 1     # RSS lever: 1 ≈ 970 MB/call
    lr_time_budget_s::Int = 180
    narrow_ctx::Bool = false          # on narrow_ctx_insufficient: retry ONCE with false
    check_divergences::Union{Symbol,Bool} = :auto  # :auto = ON standalone, OFF in farm
    algebraic_letters::Bool = false   # FindRoots tier; escalation is the cascade's job
    timeout_s::Int = 3600
    # -- process transport  ----
    # :cli (DEFAULT + the validation oracle — process-isolated, immune to the
    # F2 bug class) or :cabi (in-process ccall on the F2-PATCHED static-
    # isolated lib in tools/subtropica/deps/ — farm throughput path, 24-4456x
    # per-op vs the CLI spawn floor). :cabi is OPT-IN only; ops without an
    # ABI entry point (parse_expr/eval/period/vieta ops — the C-ABI's gap G1,
    # see hf_bridge.jl)
    # always go through the CLI regardless. Timeouts are NOT enforceable
    # in-process: the farm kill unit is the WORKER PROCESS (by PID).
    # THREADS SHARING ONE HANDLE ARE FORBIDDEN (measured: FLINT aborts with a
    # corrupt-size allocation) — the bridge
    # serializes ccalls behind a lock; farm parallelism = processes.
    transport::Symbol = Symbol(get(ENV, "SUBTROPICA_HF_TRANSPORT", "cli"))
    cabi_lib::String = get(ENV, "SUBTROPICA_HF_CABI_LIB",
        joinpath(@__DIR__, "..", "deps", "libhf_cabi_static_f2.so"))  # rebuilt by deps/build_cabi.sh
end

# ---------------------------------------------------------------------------
# GateReport — verify_laurent output. Oracle set is PINNED (RT-iii):
# (:amflow, :hiprec_sectordecomp, :mpmath). :pysecdec only with
# sanity=true and NEVER in gates; fixture-value comparison FORBIDDEN in
# gates. digits are FLOORED-INTEGER matched digits per ε-order.
# ---------------------------------------------------------------------------
struct GateReport
    pass::Bool
    oracle::Symbol
    sanity_only::Bool                          # true ⇒ this report may NOT gate
    point::Dict{Symbol,Rational{BigInt}}       # exact rational kinematic point
    digits::Dict{Int,Int}                      # ε-order => floored matched digits
    details::Dict{String,Any}                  # controls, seeds, oracle provenance
end
