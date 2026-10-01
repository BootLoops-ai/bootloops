# verify.jl — verification harness (C13).
#
# ORACLE SET IS PINNED (DESIGN.md RT-iii): gate oracles = {:amflow
# (AMFlow-port, run under the memory-capped engine wrapper — the external
# verification harness runs it, NEVER this module), :hiprec_sectordecomp (the arbitrary-precision
# parametric evaluator, tools/longhand route A; the oracle key keeps its historical name),
# :mpmath (tanh-sinh direct quadrature, low-dim convergent — live here via
# scripts/make_fixtures.py --oracle)}. :pysecdec = sanity tier ONLY
# (~1e-3..1e-10), its reports can never gate. FORM5 is NOT an oracle.
#
# STRUCTURAL ENFORCEMENT (RT-iii): comparison targets come ONLY from
# oracle-run artifacts (schema "subtropica-oracle-run-v1", kind "oracle-run",
# provenance stamp required). `load_oracle_artifact` REFUSES fixture files
# and stamp-less files, so fixture-embedded values cannot reach the
# comparator by construction. Schemas: fixtures/README.md.
#
# Report discipline: FLOORED-integer matched digits per ε-order; the WORST
# order governs; >= 30 to PASS (memory: digit counts are honest floored
# integers, never decorated).
#
# Controls: the negative control is a ×(1+1e-25·π) coefficient mutation
# (NON-RATIONAL — known pitfall: rational perturbations alias into exact
# relations) applied to a SCRATCH COPY, never in place.
#
# Precision hygiene: all decimal strings are parsed inside
# setprecision(BigFloat, …) with a +35-digit guard (known pitfalls:
# parsing at import/ambient precision, dyadic-basepoint poisoning — points
# are EXACT rationals only).

const _FIXTURE_SCHEMA = "subtropica-euler-quad-v1"
const _ORACLE_SCHEMA  = "subtropica-oracle-run-v1"
const _GATE_ORACLES   = (:amflow, :hiprec_sectordecomp, :mpmath)
const _ALL_ORACLES    = (:amflow, :hiprec_sectordecomp, :mpmath, :pysecdec)

_fixtures_dir() = normpath(joinpath(@__DIR__, "..", "fixtures"))
_make_fixtures_py() = normpath(joinpath(@__DIR__, "..", "scripts", "make_fixtures.py"))

# ---------------------------------------------------------------------------
# small parsing helpers
# ---------------------------------------------------------------------------
"""
    parse_rational(s) -> Rational{BigInt}

Parse a content-canonical rational string `"p"` or `"p/q"` (fixtures/README.md).
EXACT rationals only — float/dyadic basepoints are refused (known pitfall:
dyadic basepoints poison PSLQ)."""
function parse_rational(s::AbstractString)
    occursin(r"^-?[0-9]+(/[0-9]+)?$", strip(s)) ||
        error("subtropica verify: non-rational token '$s' — exact rationals only " *
              "(dyadic-basepoint footgun)")
    parts = split(strip(s), '/')
    length(parts) == 1 ? Rational{BigInt}(parse(BigInt, parts[1])) :
        parse(BigInt, parts[1]) // parse(BigInt, parts[2])
end

_plain(x::AbstractDict) = Dict{String,Any}(String(k) => _plain(v) for (k, v) in x)
_plain(x::AbstractVector) = Any[_plain(v) for v in x]
_plain(x) = x

_prec_bits(dps::Int) = ceil(Int, 3.322 * (dps + 35)) + 64   # +35-digit guard

# ---------------------------------------------------------------------------
# fixture loading
# ---------------------------------------------------------------------------
"""
    load_fixture(path) -> Dict{String,Any}

Load and schema-check a `subtropica-euler-quad-v1` fixture (fixtures/README.md).
Returns the parsed Dict; carries NO result values by construction."""
function load_fixture(path::AbstractString)
    isfile(path) || error("subtropica verify: fixture not found: $path")
    fx = _plain(JSON.parsefile(path))
    get(fx, "schema", "") == _FIXTURE_SCHEMA ||
        error("subtropica verify: $path is not a $_FIXTURE_SCHEMA fixture")
    for k in ("id", "divergent", "prefactor", "nu", "polys", "vars", "kinvars",
              "provenance")
        haskey(fx, k) || error("subtropica verify: fixture $path missing key '$k'")
    end
    return fx
end

# ---------------------------------------------------------------------------
# oracle-run artifacts — the ONLY admissible comparison targets
# ---------------------------------------------------------------------------
"""
    OracleArtifact

Parsed oracle-run artifact (schema `subtropica-oracle-run-v1`). Fields:
`path`, `fixture_id`, `oracle::Symbol`, `sanity_only` (FORCED true for
:pysecdec), `dps`, `minorder`, `values::Dict{order=>BigFloat}` (real,
Phase-A Euclidean), `claimed_digits::Dict{order=>Int}` (honest per-order
accuracy cap), `point::Dict{Symbol,Rational{BigInt}}`, `provenance`."""
struct OracleArtifact
    path::String
    fixture_id::String
    oracle::Symbol
    sanity_only::Bool
    dps::Int
    minorder::Int
    values::Dict{Int,BigFloat}
    claimed_digits::Dict{Int,Int}
    point::Dict{Symbol,Rational{BigInt}}
    provenance::Dict{String,Any}
end

"""
    load_oracle_artifact(path) -> OracleArtifact

Load an oracle-run artifact — the ONLY source of comparison targets
(structural enforcement of DESIGN.md RT-iii). REFUSES, with a typed error:
fixture files (fixture-embedded values are FORBIDDEN comparison targets),
files without `schema=="subtropica-oracle-run-v1"` + `kind=="oracle-run"`,
files without a provenance stamp (generated_by + timestamp_utc + method),
and unknown oracles. `sanity_only` is FORCED true for :pysecdec."""
function load_oracle_artifact(path::AbstractString)
    isfile(path) || error("subtropica verify: oracle artifact not found: $path")
    j = _plain(JSON.parsefile(path))
    sch = get(j, "schema", "")
    if sch == _FIXTURE_SCHEMA || haskey(j, "polys")
        error("subtropica verify: REFUSED — $path is a FIXTURE, and fixture-embedded " *
              "values are FORBIDDEN as comparison targets (DESIGN.md RT-iii). " *
              "Oracle values must come from a $_ORACLE_SCHEMA oracle-run artifact.")
    end
    sch == _ORACLE_SCHEMA ||
        error("subtropica verify: REFUSED — $path has schema '$sch', not $_ORACLE_SCHEMA")
    get(j, "kind", "") == "oracle-run" ||
        error("subtropica verify: REFUSED — $path lacks kind=\"oracle-run\"")
    prov = get(j, "provenance", nothing)
    (prov isa AbstractDict &&
     all(haskey(prov, k) for k in ("generated_by", "timestamp_utc", "method"))) ||
        error("subtropica verify: REFUSED — $path lacks a provenance stamp " *
              "(generated_by/timestamp_utc/method): not an oracle-RUN artifact")
    oracle = Symbol(get(j, "oracle", ""))
    oracle in _ALL_ORACLES ||
        error("subtropica verify: REFUSED — unknown oracle '$(oracle)' in $path " *
              "(pinned set: $_ALL_ORACLES, RT-iii)")
    sanity = Bool(get(j, "sanity_only", false)) || oracle === :pysecdec  # forced
    dps = Int(j["dps"])
    lau = j["laurent"]
    minorder = Int(lau["minorder"])
    values = Dict{Int,BigFloat}()
    claimed = Dict{Int,Int}()
    setprecision(BigFloat, _prec_bits(dps)) do
        for c in lau["coeffs"]
            o = Int(c["order"])
            im = parse(BigFloat, get(c, "value_im", "0"))
            iszero(im) || error("subtropica verify: nonzero imaginary part at eps^$o " *
                                "in $path — Phase-A comparator is Euclidean/real only (R4)")
            values[o] = parse(BigFloat, c["value_re"])
            claimed[o] = Int(get(c, "claimed_digits", dps))
        end
    end
    point = Dict{Symbol,Rational{BigInt}}(
        Symbol(k) => parse_rational(v) for (k, v) in get(j, "point", Dict{String,Any}()))
    return OracleArtifact(String(path), String(j["fixture_id"]), oracle, sanity,
                          dps, minorder, values, claimed, point,
                          _plain(prov))
end

# ---------------------------------------------------------------------------
# rational verification-point picker
# ---------------------------------------------------------------------------
"""
    pick_point(kinvars; letters=QQMPolyRingElem[], signs=Dict{Symbol,Int}(),
               seed=20260708, height=97, maxtries=500)
        -> Dict{Symbol,Rational{BigInt}}

Deterministic rational Euclidean kinematic point away from letter
degeneracies (SubTropica.wl STVerify:6258 logic, transliterated).

THE PICKER, documented (task requirement):
 1. Candidates are EXACT positive rationals num/den with 2 ≤ num ≤ height+1
    and den drawn from a fixed odd-prime list, generated by a deterministic
    64-bit LCG seeded with `seed` (reproducible across runs and Julia
    versions; no global RNG).
 2. HARD EXCLUSIONS: value 1 (endpoint-representability trap), value 1/3
    (a known degenerate-slice pitfall), integers, and any value already used
    for another symbol (pairwise distinct).
 3. Every polynomial in `letters` (letter polynomials, discriminants — any
    degeneracy locus the caller knows) is evaluated EXACTLY over QQ at the
    candidate point; a zero anywhere rejects the point and the search
    continues. Exact Nemo arithmetic — no float leakage (dyadic-basepoint
    footgun).
 4. `signs[k] = -1` flips symbol k negative (for oracles whose Euclidean
    chart wants raw Mandelstams s = -sb < 0); default all-positive matches
    the types.jl declared-Euclidean contract.

Throws after `maxtries` candidates (letters over-constrained — raise
`height`)."""
function pick_point(kinvars::Vector{Symbol};
                    letters::Vector{QQMPolyRingElem}=QQMPolyRingElem[],
                    signs::Dict{Symbol,Int}=Dict{Symbol,Int}(),
                    seed::Integer=20260708, height::Int=97, maxtries::Int=500)
    state = UInt64(seed == 0 ? 88172645463325252 : seed)
    nextu = () -> (state = 0x5851f42d4c957f2d * state + 0x14057b7ef767814f;
                   (state >> 33))
    primes = BigInt[3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53]
    point = Dict{Symbol,Rational{BigInt}}()
    used = Set{Rational{BigInt}}()
    for k in kinvars
        placed = false
        for _ in 1:maxtries
            num = BigInt(2 + nextu() % UInt64(height))
            den = primes[1 + Int(nextu() % UInt64(length(primes)))]
            v = num // den
            (v == 1 || v == 1 // 3 || denominator(v) == 1 || v in used) && continue
            trial = copy(point)
            trial[k] = get(signs, k, 1) < 0 ? -v : v
            # letters whose ring vars are all assigned must be nonzero
            ok = true
            for p in letters
                names = Symbol.(String.(symbols(parent(p))))
                all(n -> haskey(trial, n), names) || continue
                vals = [QQ(trial[n]) for n in names]
                if iszero(evaluate(p, vals))
                    ok = false
                    break
                end
            end
            ok || continue
            point[k] = trial[k]
            push!(used, v)
            placed = true
            break
        end
        placed || error("subtropica verify: pick_point exhausted $maxtries candidates " *
                        "for $k — letter set over-constrained; raise height")
    end
    # final full-point letter sweep
    for p in letters
        names = Symbol.(String.(symbols(parent(p))))
        all(n -> haskey(point, n), names) ||
            error("subtropica verify: letter over vars $names but kinvars=$kinvars")
        iszero(evaluate(p, [QQ(point[n]) for n in names])) &&
            error("subtropica verify: pick_point internal error — degenerate final point")
    end
    return point
end

# ---------------------------------------------------------------------------
# floored-digit comparison
# ---------------------------------------------------------------------------
"""
    floored_matched_digits(a, b; cap=10^6) -> Int

FLOORED-integer matched decimal digits between two BigFloats:
`floor(-log10(|a-b| / max(|a|,|b|)))`, clamped to [0, cap]. Exact equality
(or both zero) returns `cap` — the caller supplies the honest cap
(the artifact's per-order `claimed_digits`)."""
function floored_matched_digits(a::BigFloat, b::BigFloat; cap::Int=10^6)
    scale = max(abs(a), abs(b))
    iszero(scale) && return cap
    rel = abs(a - b) / scale
    iszero(rel) && return cap
    return clamp(floor(Int, -log10(rel)), 0, cap)
end

"""
    compare_laurent(cand::LaurentSeries{BigFloat}, art::OracleArtifact;
                    min_digits=30) -> GateReport

The comparator. Per-ε-order floored matched digits, each capped by the
artifact's honest `claimed_digits`; the WORST order governs; PASS iff the
order sets agree exactly (truncated-Laurent pole-guard footgun: a missing
deepest pole must FAIL, never be skipped) AND worst ≥ `min_digits`.
`sanity_only` artifacts (pysecdec tier) produce reports that can never gate
(`can_gate` == false) regardless of digits."""
function compare_laurent(cand::LaurentSeries{BigFloat}, art::OracleArtifact;
                         min_digits::Int=30)
    cand_orders = collect(cand.minorder:(cand.minorder + length(cand.coeffs) - 1))
    art_orders = sort!(collect(keys(art.values)))
    details = Dict{String,Any}(
        "artifact_path" => art.path,
        "fixture_id" => art.fixture_id,
        "oracle_provenance" => art.provenance,
        "rule" => "floored digits per eps-order, capped at artifact claimed_digits; " *
                  "worst order governs; order sets must match exactly (pole guard)",
        "min_digits" => min_digits,
    )
    haskey(art.provenance, "mutated") &&
        (details["mutated_artifact"] = art.provenance["mutated"])
    digits = Dict{Int,Int}()
    if cand_orders != art_orders
        details["order_mismatch"] = "candidate $(cand_orders) vs oracle $(art_orders)"
        return GateReport(false, art.oracle, art.sanity_only, art.point, digits, details)
    end
    for (i, o) in enumerate(cand_orders)
        digits[o] = floored_matched_digits(cand.coeffs[i], art.values[o];
                                           cap=art.claimed_digits[o])
    end
    worst = minimum(values(digits))
    details["worst_digits"] = worst
    details["worst_order"] = first(sort([o for o in keys(digits) if digits[o] == worst]))
    pass = worst >= min_digits
    return GateReport(pass, art.oracle, art.sanity_only, art.point, digits, details)
end

"""
    can_gate(r::GateReport) -> Bool

True iff `r` may be used as a gate: PASSed AND not sanity-tier
(:pysecdec reports never gate — DESIGN.md RT-iii)."""
can_gate(r::GateReport) = r.pass && !r.sanity_only

"""
    verify_laurent(L::LaurentSeries{BigFloat}, artifact_path::AbstractString;
                   min_digits=30, expect_fixture_id=nothing) -> GateReport

Artifact-path front-end to `compare_laurent` — the ONLY way values enter the
comparator (structural RT-iii enforcement lives in `load_oracle_artifact`).
NOTE: the `verify_laurent(L, oracle::Symbol)` form in SubTropica.jl deliberately
remains a loud stub — a bare oracle NAME is not evidence; bring an
oracle-run artifact."""
function verify_laurent(L::LaurentSeries{BigFloat}, artifact_path::AbstractString;
                        min_digits::Int=30, expect_fixture_id=nothing)
    art = load_oracle_artifact(artifact_path)
    expect_fixture_id === nothing || String(expect_fixture_id) == art.fixture_id ||
        error("subtropica verify: artifact $(art.path) is for fixture " *
              "'$(art.fixture_id)', expected '$(expect_fixture_id)'")
    return compare_laurent(L, art; min_digits=min_digits)
end

"""
    verify_fixture(fixture_path, cand::LaurentSeries{BigFloat},
                   artifact_path; min_digits=30) -> GateReport

Fixture-aware wrapper: REFUSES loudly if the fixture is flagged
`divergent` (Phase-A gate iv — divergent inputs are fatal, never compared),
checks the artifact belongs to the fixture, then runs the comparator."""
function verify_fixture(fixture_path::AbstractString, cand::LaurentSeries{BigFloat},
                        artifact_path::AbstractString; min_digits::Int=30)
    fx = load_fixture(fixture_path)
    fx["divergent"] === true &&
        error("subtropica verify: REFUSING divergent fixture '$(fx["id"])' — " *
              "Phase-A gate (iv): the pipeline must refuse divergent inputs " *
              "LOUDLY; there is nothing to verify (any 'result' is a bug)")
    return verify_laurent(cand, artifact_path;
                          min_digits=min_digits, expect_fixture_id=fx["id"])
end

# ---------------------------------------------------------------------------
# oracle adapters
# ---------------------------------------------------------------------------
"""
    run_oracle(oracle::Symbol, fixture_path; dps=60, outdir, point=nothing)
        -> artifact path (String)

Oracle-run adapters (pinned set, RT-iii). Every produced value lands in an
oracle-run ARTIFACT, never in a fixture.

- `:mpmath` — LIVE: shells to `scripts/make_fixtures.py --oracle` (tanh-sinh
  quadrature; ≤2-var, convergent, ε-free fixtures) under
  `ulimit -v 32505856`; returns the artifact path (script's last stdout
  line). `point` = Dict{Symbol,Rational} for fixtures with kinvars.
- `:amflow` — **NOT-RUN stub**: live AMFlow-port execution is the
  verification harness's job (timed pilot first; memory-capped engine
  wrapper). This adapter always throws; consume the harness's amflow
  artifacts via `load_oracle_artifact` instead.
- `:hiprec_sectordecomp` — **NOT-RUN stub** in Phase A: the engine glue to
  tools/longhand (route A, per its PARAMETRIC.md) is a harness-side
  REQUEST; artifact consumption already works.
- `:pysecdec` — refused here: sanity tier only; its artifacts are accepted
  by `load_oracle_artifact` with `sanity_only` forced true."""
function run_oracle(oracle::Symbol, fixture_path::AbstractString;
                    dps::Int=60,
                    outdir::AbstractString=joinpath(dirname(abspath(fixture_path)), "oracle"),
                    point=nothing)
    if oracle === :amflow
        error("subtropica verify: AMFlow-port oracle NOT-RUN — live AMFlow calls are " *
              "the verification harness's job (timed pilot first; memory-capped " *
              "engine wrapper). Load its artifacts with load_oracle_artifact.")
    elseif oracle === :hiprec_sectordecomp
        error("subtropica verify: hiprec_sectordecomp oracle NOT-RUN in Phase A — " *
              "engine glue (tools/longhand, PARAMETRIC.md) is a harness-side " *
              "REQUEST; artifact consumption via load_oracle_artifact works today.")
    elseif oracle === :pysecdec
        error("subtropica verify: pysecdec is SANITY-ONLY (RT-iii) and has no live " *
              "adapter here; its artifacts load with sanity_only forced true.")
    elseif oracle === :mpmath
        script = _make_fixtures_py()
        isfile(script) || error("subtropica verify: $script missing")
        ptarg = ""
        if point !== nothing
            ptarg = " --point '" * join(["$k=$(numerator(v))/$(denominator(v))"
                                         for (k, v) in point], ",") * "'"
        end
        cmd = "ulimit -v 32505856; python3 '$script' --oracle '$fixture_path' " *
              "--dps $dps --out '$outdir'" * ptarg
        out = read(`bash -c $cmd`, String)
        path = String(strip(split(strip(out), '\n')[end]))
        isfile(path) || error("subtropica verify: mpmath oracle run produced no " *
                              "artifact (output: $out)")
        return path
    end
    error("subtropica verify: unknown oracle :$oracle — pinned set $_ALL_ORACLES (RT-iii)")
end

# ---------------------------------------------------------------------------
# mutation negative control (SCRATCH COPY — never in place)
# ---------------------------------------------------------------------------
"""
    mutate_coefficient(L::LaurentSeries{BigFloat};
                       order=nothing, times=nothing) -> LaurentSeries{BigFloat}

Negative-control mutation helper for the verification harness: returns a NEW
series (the input is never touched — required practice: mutation tests on
scratch copies only) with the ε^`order` coefficient multiplied by `times`
(default `1 + 1e-25·π`, NON-RATIONAL by construction — known pitfall:
rational perturbations alias into exact relations). `order` defaults to the
lowest (deepest-pole) order. A comparator that still PASSes the mutated
copy at min_digits=30 is broken (the mutation floors it at ~24 digits)."""
function mutate_coefficient(L::LaurentSeries{BigFloat};
                            order::Union{Nothing,Int}=nothing,
                            times::Union{Nothing,BigFloat}=nothing)
    o = order === nothing ? L.minorder : order
    idx = o - L.minorder + 1
    1 <= idx <= length(L.coeffs) ||
        error("subtropica verify: mutate_coefficient order $o outside " *
              "[$(L.minorder), $(L.minorder + length(L.coeffs) - 1)]")
    prec = max(precision(L.coeffs[idx]), 256)
    t = setprecision(BigFloat, prec) do
        times === nothing ? 1 + parse(BigFloat, "1e-25") * BigFloat(pi) : times
    end
    coeffs = copy(L.coeffs)
    coeffs[idx] = coeffs[idx] * t
    ev = copy(L.evidence)
    ev["mutated"] = "eps^$o coefficient x (1+1e-25*pi) NEGATIVE CONTROL — " *
                    "scratch copy, source series untouched"
    return LaurentSeries{BigFloat}(L.minorder, coeffs, copy(L.letters),
                                   L.provenance, ev)
end

"""
    mutate_coefficient(artifact_path::AbstractString, scratch_dir;
                       order=nothing, times=nothing) -> mutated path

File flavor: copies the oracle-run artifact into `scratch_dir` and mutates
the copy's ε^`order` value by `times` (default `1 + 1e-25·π`). REFUSES to
write into the source directory (never in place — a mutation test that
edits its source is a known pitfall). The mutated copy keeps its provenance stamp
plus a `mutated` marker so comparator reports show the control."""
function mutate_coefficient(artifact_path::AbstractString, scratch_dir::AbstractString;
                            order::Union{Nothing,Int}=nothing,
                            times::Union{Nothing,BigFloat}=nothing)
    art = load_oracle_artifact(artifact_path)   # validates schema + stamp
    src_dir = normpath(dirname(abspath(artifact_path)))
    dst_dir = normpath(abspath(scratch_dir))
    dst_dir == src_dir &&
        error("subtropica verify: mutate_coefficient REFUSES to mutate in place — " *
              "scratch_dir must differ from the artifact's directory")
    mkpath(dst_dir)
    j = _plain(JSON.parsefile(artifact_path))
    o = order === nothing ? art.minorder : order
    setprecision(BigFloat, _prec_bits(art.dps)) do
        t = times === nothing ? 1 + parse(BigFloat, "1e-25") * BigFloat(pi) : times
        found = false
        for c in j["laurent"]["coeffs"]
            if Int(c["order"]) == o
                c["value_re"] = string(parse(BigFloat, c["value_re"]) * t)
                found = true
            end
        end
        found || error("subtropica verify: no eps^$o coefficient in $artifact_path")
    end
    j["provenance"]["mutated"] = "eps^$o value_re x (1+1e-25*pi) NEGATIVE CONTROL; " *
                                 "scratch copy of $artifact_path (source untouched)"
    out = joinpath(dst_dir, replace(basename(artifact_path), ".json" => "") *
                            ".MUTATED.json")
    open(out, "w") do io
        write(io, JSON.json(j))
    end
    return out
end

# ---------------------------------------------------------------------------
# synthetic-truth controls (gate A(ii) building blocks)
# ---------------------------------------------------------------------------
"""
    synthetic_controls(; prec=320) -> Vector{NamedTuple}

The synthetic-truth control set. Each entry:
`(id, fixture, artifact, candidate)` — `fixture`/`artifact` are the
committed paths (fixtures/README.md inventory) and `candidate()` builds the
known-closed-form LaurentSeries{BigFloat} via Nemo/Arb at `prec` bits — an
implementation INDEPENDENT of the mpmath artifact chain (Arb/FLINT vs
mpmath+sympy), so a self-test PASS is a two-implementation agreement. The
mandatory NEGATIVE control is `mutate_coefficient` applied to a candidate
(or to a scratch copy of an artifact) — it must FAIL the comparator."""
function synthetic_controls(; prec::Int=320)
    fdir = _fixtures_dir()
    RR = ArbField(prec + 64)
    bf(x) = setprecision(() -> BigFloat(x), BigFloat, prec)
    return [
        (id = "synth_log2_1var",
         fixture = joinpath(fdir, "synth_log2_1var.json"),
         artifact = joinpath(fdir, "oracle", "synth_log2_1var.mpmath.json"),
         candidate = () -> LaurentSeries(0, BigFloat[bf(log(RR(2)))])),
        (id = "synth_zeta2_2var",
         fixture = joinpath(fdir, "synth_zeta2_2var.json"),
         artifact = joinpath(fdir, "oracle", "synth_zeta2_2var.mpmath.json"),
         candidate = () -> LaurentSeries(0, BigFloat[bf(zeta(RR(2)))])),
        (id = "synth_zeta3_beta_1var",
         fixture = joinpath(fdir, "synth_zeta3_beta_1var.json"),
         artifact = joinpath(fdir, "oracle", "synth_zeta3_beta_1var.mpmath.json"),
         candidate = () -> LaurentSeries(0, BigFloat[
             bf(one(RR)), bf(-3 * one(RR)),
             bf(9 - 2 * zeta(RR(2))),
             bf(-27 + 6 * zeta(RR(2)) + 6 * zeta(RR(3)))])),
    ]
end

# ---------------------------------------------------------------------------
# Phase-A stubs owned here (loud, typed — silence is never an option)
# ---------------------------------------------------------------------------
"""
    eval_symbolic(L::LaurentSeries{HlogExpr}, point; dps=60)
        -> LaurentSeries{BigFloat}

Evaluate the SYMBOLIC pipeline result numerically per ε-order at an EXACT
rational kinematic `point` (Dict{Symbol,Rational}; empty Dict for
constant results). Wired (stub-wiring rule
item 5) on the HlogExpr surface:

- coefficient strings re-parsed through the Nemo fraction field of
  `e.vars` (`_parse_rat_frac`) and evaluated EXACTLY over QQ at the point
  (no float leakage — dyadic-basepoint footgun); a vanishing denominator
  is a typed refusal (degenerate point — use `pick_point`).
- atoms: HConst(:Pi/:EulerGamma/:Catalan/:Log2), HZeta (INDEPENDENT
  ginac chain via `mzv_bigfloat`), HLogA (argument evaluated exactly at
  the point, must be > 0), HHlog (z and rational-function letters
  evaluated exactly at the point, then `ginac_gpl_G` — real-line guard
  inside). HPeriod atoms (resolve_periods strict=false explicit-atom
  path) evaluate NUMERICALLY : rational-letter words
  via `zero_inf_period_ginac`, algebraic-letter words via
  `zero_inf_period_alg` + the LaurentSeries.letters table; HAlgLetter
  atoms via `alg_letter_roots` (monic law, root-checked). Both are
  complex per TERM — the real-line law is enforced on each coefficient
  SUM (typed refusal if conjugate terms fail to cancel). HDelta /
  HConst(:I) still REFUSE typed (Euclidean guard R4).
- all decimal work at `_prec_bits(dps)` (= dps+35-digit guard against
  precision-at-parse pitfalls).

The output LaurentSeries{BigFloat} feeds `compare_laurent`/`verify_fixture`
directly (kinematic-fixture gating path)."""
function eval_symbolic(L::LaurentSeries{HlogExpr}, point; dps::Int=60)
    pt = Dict{Symbol,Rational{BigInt}}()
    for (k, v) in pairs(point)
        pt[Symbol(k)] = v isa Rational ? Rational{BigInt}(v) :
                        v isa Integer ? Rational{BigInt}(BigInt(v)) :
                        parse_rational(string(v))
    end
    coeffs = setprecision(BigFloat, _prec_bits(dps)) do
        BigFloat[_eval_hlogexpr_at(c, pt; dps=dps, letters=L.letters)
                 for c in L.coeffs]
    end
    ev = copy(L.evidence)
    ev["eval_symbolic_point"] = Dict(String(k) => "$(numerator(v))/$(denominator(v))"
                                     for (k, v) in pt)
    ev["eval_symbolic_dps"] = dps
    return LaurentSeries{BigFloat}(L.minorder, coeffs, copy(L.letters),
                                   L.provenance, ev)
end
eval_symbolic(L::LaurentSeries{BigFloat}, point; dps::Int=60) =
    error("subtropica verify: eval_symbolic takes the SYMBOLIC result " *
          "(LaurentSeries{HlogExpr}); numeric candidates go straight to " *
          "compare_laurent")

# exact QQ evaluation of a fraction-field element at the point
function _eval_frac_exact(fr, vars::Vector{Symbol},
                          pt::Dict{Symbol,Rational{BigInt}}, what::String)
    vals = [QQ(pt[v]) for v in vars]
    den = evaluate(denominator(fr), vals)
    iszero(den) && error("subtropica verify: eval_symbolic hit a vanishing " *
                         "denominator in $what at the point — degenerate " *
                         "kinematic point (use pick_point)")
    q = evaluate(numerator(fr), vals) // den
    return Rational{BigInt}(BigInt(numerator(q)), BigInt(denominator(q)))
end

function _eval_hlogexpr_at(e::HlogExpr, pt::Dict{Symbol,Rational{BigInt}};
                           dps::Int=60, letters::Vector{AlgLetter}=AlgLetter[])
    for v in e.vars
        haskey(pt, v) || error("subtropica verify: eval_symbolic point is " *
            "missing kinematic symbol '$v' (expr vars: $(e.vars))")
    end
    ctx = _pctx(e.vars, _letters_of(e))
    # COMPLEX accumulation (word-level Wm/Wp numeric leg):
    # HPeriod words over algebraic-letter roots and HAlgLetter atoms are
    # generally complex per TERM; the real-line law is enforced on the SUM —
    # residual imaginary parts must cancel across conjugate Wm/Wp terms, and
    # a non-negligible survivor is a typed refusal, never a silent real().
    tot = zero(Complex{BigFloat})
    for t in e.terms
        fr = _parse_rat_frac(t.coef, ctx)
        q = _eval_frac_exact(fr, e.vars, pt, "coefficient '$(t.coef)'")
        v = Complex{BigFloat}(BigFloat(numerator(q)) / BigFloat(denominator(q)))
        for (a, p) in t.atoms
            v *= _eval_atom_at(a, e.vars, pt, ctx; dps=dps, letters=letters)^p
        end
        tot += v
    end
    abs(imag(tot)) <= max(abs(real(tot)), one(BigFloat)) *
                      big(10.0)^(-(dps - 8)) ||
        throw(SubTropicaEvalError("eval_symbolic: non-negligible residual " *
            "imaginary part $(Float64(imag(tot))) in a coefficient — " *
            "conjugate Wm/Wp period terms did NOT cancel (real-line law; " *
            "check the letter table / word set)"))
    return real(tot)
end

function _eval_atom_at(a::HAtom, vars::Vector{Symbol},
                       pt::Dict{Symbol,Rational{BigInt}}, ctx; dps::Int=60,
                       letters::Vector{AlgLetter}=AlgLetter[])
    if a isa HPeriod
        # word-level numeric leg : residual ZeroInfPeriod
        # words kept as explicit atoms (resolve_periods strict=false) are
        # CONSTANTS — rational-letter words via the independent large-z
        # extrapolation, algebraic-letter words via the AlgLetter table
        # (complex; imag-cancellation asserted by the caller on the sum).
        if any(w -> occursin(_ALGLETTER_TOKEN_RE, w), a.word)
            return zero_inf_period_alg(a.word, letters; digits=dps)
        else
            return Complex{BigFloat}(zero_inf_period_ginac(a.word; digits=dps))
        end
    elseif a isa HAlgLetter
        li = findfirst(l -> l.idx == a.idx, letters)
        li === nothing && throw(SubTropicaEvalError("eval_symbolic: " *
            "HAlgLetter($(a.kind),$(a.idx)) has no entry in the result's " *
            "letter table — back-substitute (normalize_wmwp, RT-v) or pass " *
            "the LaurentSeries.letters table BEFORE numeric evaluation"))
        (Wm, Wp, sd) = alg_letter_roots(letters[li]; digits=dps)
        return a.kind === :Wm ? Wm : a.kind === :Wp ? Wp :
               a.kind === :sqrt_disc ? sd :
               throw(SubTropicaEvalError("eval_symbolic: unknown HAlgLetter " *
                   "kind $(a.kind)"))
    end
    if a isa HConst || a isa HZeta
        return _atom_bigfloat(a; digits=dps)   # typed refusal on :I inside
    elseif a isa HLogA
        arg = _eval_frac_exact(_parse_rat_frac(a.arg, ctx), vars, pt,
                               "Log argument '$(a.arg)'")
        arg > 0 || throw(SubTropicaEvalError(
            "eval_symbolic: Log[$(a.arg)] evaluates to $(arg) ≤ 0 at the " *
            "point — real-line evaluator (R4: no physical-region continuation)"))
        return log(BigFloat(numerator(arg)) / BigFloat(denominator(arg)))
    elseif a isa HHlog
        z = _eval_frac_exact(_parse_rat_frac(a.z, ctx), vars, pt,
                             "Hlog argument '$(a.z)'")
        letters = Rational{BigInt}[
            _eval_frac_exact(_parse_rat_frac(w, ctx), vars, pt,
                             "Hlog letter '$w'") for w in a.word]
        zs = denominator(z) == 1 ? string(numerator(z)) :
             string(numerator(z)) * "/" * string(denominator(z))
        return ginac_gpl_G(letters, zs; digits=dps + 10)
    end
    throw(SubTropicaEvalError(
        "eval_symbolic: $(typeof(a)) atom has no direct numeric value — " *
        "resolve ZeroInfPeriod keys (resolve_periods), back-substitute " *
        "algebraic letters (normalize_wmwp, RT-v), and keep delta residues " *
        "out of Euclidean results (R4) BEFORE numeric evaluation"))
end

"""
    mzv_provenance_gate(; nrules=25, dps=60,
                        table_path=HFConfig().mzv_data_path,
                        seed=20260901,
                        expect_sha256=get(ENV, "SUBTROPICA_MZV_SHA256", ""))
        -> GateReport

Gate A(iii)/R7: verify `nrules` random rules of the MZV reduction table
(`mzv_reductions.json`) LHS == RHS to ≥ `dps` matched digits via the
INDEPENDENT ginac chain (`mzv_bigfloat`/`hlog_const_bigfloat` — ginac_gpl
CLI, a different implementation from the HF engine that authored the
table). Both sides are evaluated at `dps + 30` working digits and compared
with `floored_matched_digits` capped at `dps + 20` (honest ceiling); the
WORST rule governs; PASS iff worst ≥ dps.

Provenance pin (R7): the table's sha256 is computed and recorded in
`details["table_sha256"]`; when `expect_sha256` is non-empty (default: the
`SUBTROPICA_MZV_SHA256` env var, template line in scripts/env.sh) a
mismatch is a typed refusal BEFORE any evaluation — any table change
re-triggers the gate against the new hash.

Sampling is DETERMINISTIC (seeded 64-bit LCG, same generator as
`pick_point`; no global RNG) and without replacement; the sampled lhs
tokens are recorded in `details["sampled_lhs"]`. Rules whose LHS is
divergent in the pinned ASCENDING convention (outermost index +1 — HF
tables carry shuffle-regularized entries; this ginac chain has no MZV
regularization) are EXCLUDED from the sample universe and counted loudly
in `details["excluded_divergent_lhs"]` — never silently skipped mid-run.

MANDATORY negative control (synthetic-truth law): the first sampled rule
is re-compared with its RHS value mutated by ×(1+1e-25·π) (NON-RATIONAL —
rational perturbations alias into exact relations); the mutated compare
must FAIL the dps bar (it floors at ~24 digits, hence the `dps ≥ 30`
requirement) or the gate throws — a comparator that passes the control is
broken. Control digits land in `details["negative_control_digits"]`.

Typed refusals: missing table (names `SUBTROPICA_MZV_DATA`), pin
mismatch, malformed pin, empty/absent `reductions`, dps < 30, missing
ginac_gpl CLI (names `SUBTROPICA_GINAC_GPL` — battery legs skip BY NAME
on this refusal), unparseable rule (a rule outside the reduction-output
grammar must be investigated, not skipped)."""
function mzv_provenance_gate(; nrules::Int=25, dps::Int=60,
                             table_path::AbstractString=HFConfig().mzv_data_path,
                             seed::Integer=20260901,
                             expect_sha256::AbstractString=get(ENV, "SUBTROPICA_MZV_SHA256", ""))
    nrules >= 1 || error("subtropica verify: mzv_provenance_gate nrules must be >= 1")
    dps >= 30 || error("subtropica verify: mzv_provenance_gate dps=$dps < 30 — " *
                       "below the house >=30d bar the (1+1e-25*pi) mutation " *
                       "control cannot bite; refusing")
    isfile(table_path) || error("subtropica verify: MZV reduction table not " *
        "found at $(table_path) — the engine's mzv_reductions.json is not " *
        "vendored; point SUBTROPICA_MZV_DATA (or table_path) at your engine " *
        "build's table (GUIDE.md ENGINE section)")
    table_sha = _sha256_hex(read(table_path))
    if !isempty(expect_sha256)
        occursin(r"^[0-9a-f]{64}$", expect_sha256) ||
            error("subtropica verify: expect_sha256 '$(expect_sha256)' is not " *
                  "64 lowercase hex chars — malformed pin (SUBTROPICA_MZV_SHA256)")
        expect_sha256 == table_sha ||
            error("subtropica verify: MZV table sha256 MISMATCH — " *
                  "$(table_path) hashes to $(table_sha), pinned " *
                  "$(expect_sha256) (SUBTROPICA_MZV_SHA256). The table " *
                  "changed: re-pin AND re-run this gate against the new table (R7)")
    end
    gbin = _ginac_bin()
    isfile(gbin) || error("subtropica verify: mzv_provenance_gate needs the " *
        "ginac_gpl CLI (independent evaluator) — not found at $(gbin); set " *
        "SUBTROPICA_GINAC_GPL")
    doc = _plain(JSON.parsefile(table_path))
    rules = get(doc, "reductions", Any[])
    isempty(rules) && error("subtropica verify: $(table_path) has no " *
        "\"reductions\" — not an MZV reduction table")
    # Sample universe: convergent-LHS rules only (see docstring); every
    # exclusion is recorded, never silent.
    universe = Int[]
    excluded = String[]
    for (i, r) in enumerate(rules)
        (r isa AbstractDict && haskey(r, "lhs") && haskey(r, "rhs")) ||
            error("subtropica verify: reduction rule $i lacks lhs/rhs — " *
                  "malformed table")
        idx = _token_zeta(String(r["lhs"]))
        if !isempty(idx) && abs(idx[end]) == 1 && idx[end] > 0
            push!(excluded, String(r["lhs"]))
        else
            push!(universe, i)
        end
    end
    isempty(universe) && error("subtropica verify: every rule in " *
        "$(table_path) has a divergent (shuffle-regularized) LHS — nothing " *
        "this ginac chain can evaluate; the gate cannot run")
    # Deterministic sample without replacement (pick_point's LCG; no global RNG).
    state = UInt64(seed == 0 ? 88172645463325252 : seed)
    nextu = () -> (state = 0x5851f42d4c957f2d * state + 0x14057b7ef767814f;
                   (state >> 33))
    k = min(nrules, length(universe))
    pool = copy(universe)
    picked = Int[]
    for _ in 1:k
        j = 1 + Int(nextu() % UInt64(length(pool)))
        push!(picked, pool[j])
        deleteat!(pool, j)
    end
    eval_digits = dps + 30
    cap = dps + 20
    digits = Dict{Int,Int}()
    sampled_lhs = String[]
    rule_digits = Dict{String,Int}()
    first_pair = nothing   # (lhs_val, rhs_val) of the first rule, for the control
    for (ord, i) in enumerate(picked)
        r = rules[i]
        lhs_tok = String(r["lhs"])
        push!(sampled_lhs, lhs_tok)
        lhs_val = mzv_bigfloat(_token_zeta(lhs_tok); digits=eval_digits)
        rhs_val = hlog_const_bigfloat(_parse_reduced(String(r["rhs"]));
                                      digits=eval_digits)
        d = setprecision(BigFloat, _prec_bits(eval_digits)) do
            floored_matched_digits(lhs_val, rhs_val; cap=cap)
        end
        digits[ord] = d
        rule_digits[lhs_tok] = d
        first_pair === nothing && (first_pair = (lhs_val, rhs_val))
    end
    # Mandatory negative control: mutated first rule must FAIL the bar.
    control_digits = setprecision(BigFloat, _prec_bits(eval_digits)) do
        t = 1 + parse(BigFloat, "1e-25") * BigFloat(pi)
        floored_matched_digits(first_pair[1], first_pair[2] * t; cap=cap)
    end
    control_digits < dps ||
        error("subtropica verify: mzv_provenance_gate NEGATIVE CONTROL " *
              "PASSED the comparator ($(control_digits) >= $(dps) digits on " *
              "a x(1+1e-25*pi) mutation) — the comparator is broken; " *
              "refusing to report a gate result")
    worst = minimum(values(digits))
    worst_lhs = first(sort([l for (l, d) in rule_digits if d == worst]))
    details = Dict{String,Any}(
        "gate" => "A(iii)/R7 MZV table provenance",
        "table_path" => String(table_path),
        "table_sha256" => table_sha,
        "table_sha256_pinned" => isempty(expect_sha256) ? "UNPINNED" :
                                 String(expect_sha256),
        "rule" => "per-rule floored matched digits LHS-vs-RHS via ginac, " *
                  "capped at dps+20; worst rule governs; PASS iff worst >= dps",
        "min_digits" => dps,
        "eval_digits" => eval_digits,
        "seed" => Int(seed),
        "nrules_requested" => nrules,
        "nrules_checked" => k,
        "n_rules_in_table" => length(rules),
        "excluded_divergent_lhs" => excluded,
        "sampled_lhs" => sampled_lhs,
        "rule_digits" => rule_digits,
        "worst_digits" => worst,
        "worst_lhs" => worst_lhs,
        "negative_control_digits" => control_digits,
        "evaluator" => "ginac_gpl CLI (mzv_bigfloat/hlog_const_bigfloat) — " *
                       "independent of the table-authoring engine",
    )
    return GateReport(worst >= dps, :ginac, false,
                      Dict{Symbol,Rational{BigInt}}(), digits, details)
end
