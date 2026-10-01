# assemble.jl — Result assembly (C12).
#
# Phase A = degenerate single face: per-ε-order HF results × prefactor
# Laurent (eps_expand.jl's loggamma_series output) → LaurentSeries{HlogExpr} with
# provenance metadata. The FOLD is written generally (multi-face pole
# padding alignment) because Phase B1's farm feeds the same code path; the
# upstream analog assembles per-face per-order sums into SeriesData with
# the e^(±L·γ_E·ε) normalization adjustment (SubTropica.wl:8860-8880) and
# orchestrates faces in STExpandIntegral (SubTropica.wl:10906, 11029);
# minOrder bookkeeping at SubTropica.wl:11188.
#
# POLICY (DESIGN.md R6, CONTRACTS.md §a check_divergences): Phase A has
# exactly ONE mode, :standalone, and check_divergences is HARD-ON there —
# assembling any result whose provenance does not carry
# `check_divergences == true` THROWS (typed DivergenceCheckViolation).
# This inverts the upstream silent-0-on-divergent behavior (confirmed on the
# binary, control fixture fixtures/divergent_xinv_1var.json;
# SubTropica.wl:12374-12384 policy, parity fix wl:12589).
# The face-sum mode (check_divergences OFF per counter-term) is Phase B1's
# farm and is REFUSED here, loudly.
#
# HlogExpr arithmetic: the canonical owner is hlogexpr.jl (hlog_add/
# hlog_mul — "implemented HERE so the semantics have one owner"). Those
# implementations are live, so the scaffold-era local fallback
# (_a3_hadd_local/_a3_hmul_local + probe-routing) was DELETED; _a3_hadd/
# _a3_hmul below delegate directly. The small _a3_* helpers that remain
# (tryrat/ratstr/norm_atoms/atom_key) are still used by the restricted
# mzv-reduction-output parser and the tests.

# ---------------------------------------------------------------------------
# Typed errors
# ---------------------------------------------------------------------------

"""
    DivergenceCheckViolation

Thrown when assembling results whose provenance does not certify
`check_divergences:true` (HARD-ON for Phase A standalone mode — DESIGN.md
R6 / CONTRACTS.md §a), or when the Phase-B-only `:face_sum` mode is
requested. Never a warning: the upstream package silently returns 0 on
divergent input with this flag off (confirmed on the binary).
"""
struct DivergenceCheckViolation <: Exception
    msg::String
end
Base.showerror(io::IO, e::DivergenceCheckViolation) =
    print(io, "DivergenceCheckViolation: ", e.msg)

"Missing/unusable provenance on assembly input (LR order, pins)."
struct AssembleProvenanceError <: Exception
    msg::String
end
Base.showerror(io::IO, e::AssembleProvenanceError) =
    print(io, "AssembleProvenanceError: ", e.msg)

"Inconsistent algebraic-letter tables at union (same (poly,var), different data)."
struct LetterTableMismatch <: Exception
    msg::String
end
Base.showerror(io::IO, e::LetterTableMismatch) =
    print(io, "LetterTableMismatch: ", e.msg)

"""
    LaurentTruncationError

The requested output order cannot be honestly formed from the available
input orders (truncated-Laurent pole-guard footgun class, DESIGN.md R2):
either an integrand-series or prefactor-series coefficient that the
convolution needs is missing. Callers must supply explicit zeros for
orders known to vanish — we never assume them.
"""
struct LaurentTruncationError <: Exception
    msg::String
end
Base.showerror(io::IO, e::LaurentTruncationError) =
    print(io, "LaurentTruncationError: ", e.msg)

"Atom (or coef string) with no pinned ginac-evaluable form."
struct GinshExportError <: Exception
    msg::String
end
Base.showerror(io::IO, e::GinshExportError) =
    print(io, "GinshExportError: ", e.msg)

"MZV reduction pass failure (op error, unknown symbol, unparseable output)."
struct MzvReduceError <: Exception
    msg::String
end
Base.showerror(io::IO, e::MzvReduceError) =
    print(io, "MzvReduceError: ", e.msg)

# ---------------------------------------------------------------------------
# FaceResult — one face's per-ε-order contributions (the general-fold input).
# Phase A produces exactly one; Phase B1's farm (C11) will produce one per
# face in Σ_div. NOTE: if Phase B wants this shared across
# modules, promote it to types.jl (design decision; flagged).
# ---------------------------------------------------------------------------

"""
    FaceResult(minorder, orders)

Per-face assembly input: `orders[k]` is the vector of contributions at
ε-order `minorder + k - 1`; each contribution is
`(expr::HlogExpr, letters::Vector{AlgLetter})` with `letters` the HF
response's own algebraic-letter table (per-response indexing — remapped
into the unioned table during assembly). An order with no contributions
is an explicit empty vector, never a skipped slot.
"""
struct FaceResult
    minorder::Int
    orders::Vector{Vector{Tuple{HlogExpr,Vector{AlgLetter}}}}
end

# ---------------------------------------------------------------------------
# Exact local HlogExpr helpers (+ probe-routing to hlogexpr.jl when available)
# ---------------------------------------------------------------------------

const _A3_RAT_RE = r"^[+-]?[0-9]+(?:/[0-9]+)?$"

"Parse a plain-rational coef string, else `nothing` (general rational
functions of kinematics stay opaque strings here)."
function _a3_tryrat(s::AbstractString)
    occursin(_A3_RAT_RE, s) || return nothing
    parts = split(s, '/')
    num = parse(BigInt, parts[1])
    den = length(parts) == 2 ? parse(BigInt, parts[2]) : big(1)
    den == 0 && return nothing
    return num // den
end

"CONTENT-CANONICAL rational string (CONTRACTS.md §d): gcd-reduced, q>0,
sign on the numerator, bare integer when q==1."
_a3_ratstr(r::Rational{BigInt}) =
    denominator(r) == 1 ? string(numerator(r)) :
                          string(numerator(r), "/", denominator(r))

# Deterministic atom ordering following the §d kind order (HConst, HZeta,
# HLogA, HHlog, HPeriod, HAlgLetter, HDelta). Final canonical text form is
# hlogexpr.jl's canonical_text; here we only need a stable sort.
_a3_atom_rank(::HConst) = 1
_a3_atom_rank(::HZeta) = 2
_a3_atom_rank(::HLogA) = 3
_a3_atom_rank(::HHlog) = 4
_a3_atom_rank(::HPeriod) = 5
_a3_atom_rank(::HAlgLetter) = 6
_a3_atom_rank(::HDelta) = 7
_a3_atom_str(a::HConst) = string(a.name)
_a3_atom_str(a::HZeta) = join(a.idx, ",")
_a3_atom_str(a::HLogA) = a.arg
_a3_atom_str(a::HHlog) = a.z * "|" * join(a.word, ",")
_a3_atom_str(a::HPeriod) = join(a.word, ",")
_a3_atom_str(a::HAlgLetter) = string(a.kind, "_", a.idx)
_a3_atom_str(a::HDelta) = string(a.var)
_a3_atom_key(a::HAtom) = (_a3_atom_rank(a), _a3_atom_str(a))

"Merge duplicate atoms (sum powers, drop zeros) and sort deterministically."
function _a3_norm_atoms(atoms::Vector{Pair{HAtom,Int}})
    acc = Dict{HAtom,Int}()
    order = HAtom[]
    for (a, p) in atoms
        if haskey(acc, a)
            acc[a] += p
        else
            acc[a] = p
            push!(order, a)
        end
    end
    out = Pair{HAtom,Int}[a => acc[a] for a in order if acc[a] != 0]
    sort!(out; by=pr -> _a3_atom_key(pr.first))
    return out
end

# --- HlogExpr arithmetic: DIRECT delegation to hlogexpr.jl's exact
# implementations (hlogexpr.jl). The scaffold-era local fallback
# (_a3_hadd_local/_a3_mulcoef/_a3_hmul_local, ~80 LOC) and its
# availability probe (historical private-arithmetic probes) were deleted
# — hlog_add/hlog_mul are live and
# suite-gated; two arithmetic implementations were a semantics-drift
# surface.
_a3_hadd(a::HlogExpr, b::HlogExpr) = hlog_add(a, b)
_a3_hmul(a::HlogExpr, b::HlogExpr) = hlog_mul(a, b)
function _a3_hscale(a::HlogExpr, c::Rational{BigInt})
    c == 0 && return HlogExpr(HTerm[], a.vars)
    return _a3_hmul(a, HlogExpr([HTerm(_a3_ratstr(c), Pair{HAtom,Int}[])], Symbol[]))
end
_a3_hscale(a::HlogExpr, c::Union{Integer,Rational}) = _a3_hscale(a, Rational{BigInt}(c))

_a3_hlog_iszero(e::HlogExpr) = isempty(e.terms)
_a3_hlog_zero() = HlogExpr()

# ---------------------------------------------------------------------------
# Letter-table union + atom remapping
# ---------------------------------------------------------------------------

"""
    union_letters(tables::Vector{Vector{AlgLetter}})
        -> (letters::Vector{AlgLetter}, remaps::Vector{Dict{Int,Int}})

Union algebraic-letter tables across HF responses. Entries are identified
by `(poly, var)`; a repeat must agree on lc/sum/product/disc — any
mismatch is an upstream bug and throws [`LetterTableMismatch`](@ref)
(Wm/Wp state is PER-PROCESS in HF, CONTRACTS.md §a — cross-response index
collisions are expected and are exactly what the remap fixes).
`remaps[i][old_idx] = unioned idx`; unioned entries are renumbered 1..n.
"""
function union_letters(tables::Vector{Vector{AlgLetter}})
    letters = AlgLetter[]
    seen = Dict{Tuple{String,Symbol},Int}()
    remaps = Dict{Int,Int}[]
    for tab in tables
        remap = Dict{Int,Int}()
        for L in tab
            key = (L.poly, L.var)
            j = get(seen, key, 0)
            if j == 0
                push!(letters, AlgLetter(length(letters) + 1, L.poly, L.var,
                                         L.lc, L.sum, L.product, L.disc))
                j = length(letters)
                seen[key] = j
            else
                M = letters[j]
                (L.lc == M.lc && L.sum == M.sum &&
                 L.product == M.product && L.disc == M.disc) ||
                    throw(LetterTableMismatch(
                        "same (poly,var)=(" * L.poly * "," * string(L.var) *
                        ") with different lc/sum/product/disc across " *
                        "responses — upstream bug, refusing"))
            end
            remap[L.idx] = j
        end
        push!(remaps, remap)
    end
    return (letters, remaps)
end

"""
    remap_letters(e::HlogExpr, remap::Dict{Int,Int}) -> HlogExpr

Rewrite per-response `HAlgLetter` indices into the unioned table's
indices. An index absent from the remap throws (the expr references a
letter its own response never declared).

HPeriod WORD tokens (`Wm_i`/`Wp_i`/`WmOverWp_i`/
`sqrt_disc_i` inside residual ZeroInfPeriod words — the letters tier
mints them live, eq 4.11) are remapped too, in ONE regex pass (no
sequential-substitution collisions). Before this fix the atom indices
were remapped while word tokens kept the per-response index — a
silent-wrong path the moment two responses union with an index shift
(the word-level numeric leg resolves word tokens against the UNIONED
table)."""
function remap_letters(e::HlogExpr, remap::Dict{Int,Int})
    isempty(remap) && return e
    terms = HTerm[]
    for t in e.terms
        atoms = Pair{HAtom,Int}[]
        for (a, p) in t.atoms
            if a isa HAlgLetter
                haskey(remap, a.idx) || throw(LetterTableMismatch(
                    "HAlgLetter idx $(a.idx) not present in its response's" *
                    " letter table — refusing"))
                push!(atoms, HAlgLetter(a.kind, remap[a.idx]) => p)
            elseif a isa HPeriod &&
                   any(w -> occursin(_ALGLETTER_TOKEN_RE, w), a.word)
                push!(atoms, HPeriod(String[_remap_word_tokens(w, remap)
                                            for w in a.word]) => p)
            else
                push!(atoms, a => p)
            end
        end
        push!(terms, HTerm(t.coef, atoms))
    end
    return HlogExpr(terms, e.vars)
end

const _WORD_ANY_TOKEN_RE = r"(Wm|Wp|WmOverWp|sqrt_disc)_(\d+)"

# Single-pass token remap inside one word letter (collision-safe: the
# whole string is rewritten in one `replace` sweep, so remap {1=>2,2=>1}
# cannot double-apply).
function _remap_word_tokens(w::String, remap::Dict{Int,Int})
    return replace(w, _WORD_ANY_TOKEN_RE => function (m)
        pm = match(_WORD_ANY_TOKEN_RE, m)
        i = parse(Int, pm.captures[2])
        haskey(remap, i) || throw(LetterTableMismatch(
            "ZeroInfPeriod word token '$m' references letter idx $i " *
            "absent from its response's letter table — refusing"))
        return pm.captures[1] * "_" * string(remap[i])
    end)
end

# ---------------------------------------------------------------------------
# Laurent padding/alignment (the general fold) + convolution
# ---------------------------------------------------------------------------

"""
    align_laurent(pieces::Vector{<:Tuple{Int,Vector{T}}}; zerofn=HlogExpr)
        -> (minorder::Int, padded::Vector{Vector{T}})

Pole-padding alignment across faces: bring every piece
`(minorder_i, coeffs_i)` to the common `minorder = min_i minorder_i` and a
common top order, padding with `zerofn()`. Phase A feeds a single piece;
Phase B1's face sum uses the same fold (faces reach different pole depths
— SubTropica.wl:11188 minOrder bookkeeping)."""
function align_laurent(pieces::Vector{<:Tuple{Int,<:Vector}}; zerofn=HlogExpr)
    isempty(pieces) && throw(ArgumentError("align_laurent: no pieces"))
    gmin = minimum(p[1] for p in pieces)
    gtop = maximum(p[1] + length(p[2]) - 1 for p in pieces)
    padded = [begin
                  (mo, cs) = p
                  front = [zerofn() for _ in 1:(mo - gmin)]
                  back  = [zerofn() for _ in 1:(gtop - (mo + length(cs) - 1))]
                  vcat(front, collect(cs), back)
              end
              for p in pieces]
    return (gmin, padded)
end

# ---------------------------------------------------------------------------
# Provenance pins
# ---------------------------------------------------------------------------

const _FILE_SHA_MEMO = Dict{String,String}()
function _file_sha256(path::AbstractString)
    get!(_FILE_SHA_MEMO, String(path)) do
        isfile(path) || return "unavailable:no-such-file"
        # GNU coreutils `sha256sum`, or stock-macOS `shasum` (Digest::SHA).
        # This stamps provenance AFTER the assembly has succeeded, so it
        # must never throw: any failure degrades to a loud unavailable:*
        # sentinel (same idiom as the missing-file case above).
        cmd = Sys.which("sha256sum") !== nothing ? `sha256sum $path` :
              Sys.which("shasum") !== nothing    ? `shasum -a 256 $path` :
              nothing
        cmd === nothing && return "unavailable:no-sha256-tool-on-path"
        out = try
            read(pipeline(cmd, stderr=devnull), String)
        catch err
            return "unavailable:sha256-shellout-failed:$(typeof(err))"
        end
        toks = split(out)
        (!isempty(toks) && occursin(r"^[0-9a-f]{64}$", toks[1])) ?
            String(toks[1]) : "unavailable:unexpected-sha256-output"
    end
end

# The wrapper is a .sh; the binary pin is the sibling ELF (env.sh comment).
function _hf_bin_pin(cfg::HFConfig)
    elf = endswith(cfg.bin, ".sh") ? joinpath(dirname(cfg.bin), "hyperflint") :
                                     cfg.bin
    isfile(elf) || (elf = cfg.bin)
    sha = _file_sha256(elf)
    startswith(sha, "unavailable") &&
        @warn "subtropica assemble: HF binary not found; provenance sha unavailable" path = elf
    return (elf, sha)
end

const _PINS_MEMO = Ref{Union{Nothing,Tuple{String,String}}}(nothing)
"Parse the two upstream SHAs out of reference/PROVENANCE.md (vendored, in-tree)."
function _upstream_pins()
    if _PINS_MEMO[] === nothing
        path = normpath(joinpath(@__DIR__, "..", "reference", "PROVENANCE.md"))
        isfile(path) || throw(AssembleProvenanceError(
            "reference/PROVENANCE.md missing at $(path) — cannot stamp " *
            "upstream pins; refusing to assemble unpinned results"))
        txt = read(path, String)
        mwl = match(r"source clone\s*\|\s*`([0-9a-f]{7,40})`", txt)
        mbin = match(r"prebuilt binary\s*\|\s*`([0-9a-f]{7,40})`", txt)
        (mwl === nothing || mbin === nothing) &&
            throw(AssembleProvenanceError(
                "could not parse upstream SHAs from $(path) — file format " *
                "changed; fix _upstream_pins()"))
        _PINS_MEMO[] = (String(mwl.captures[1]), String(mbin.captures[1]))
    end
    return _PINS_MEMO[]
end

# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------

"""
    assemble_faces(faces::Vector{FaceResult},
                   prefactor_series::Tuple{Int,Vector{HlogExpr}},
                   requested_order::Int;
                   evidence=Dict{String,Any}(), mode=:standalone,
                   cfg=HFConfig(), transport=_default_transport,
                   mzv=false)
        -> LaurentSeries{HlogExpr}

General assembly fold (upstream analog: STExpandIntegral orchestration,
SubTropica.wl:10906/11029 + SeriesData assembly wl:8860-8880):

1. POLICY gates: `mode=:standalone` is the ONLY Phase-A mode and requires
   `evidence["check_divergences"] === true` (throws
   [`DivergenceCheckViolation`](@ref) otherwise — the inputs were
   produced by HF calls and their flag state is part of their provenance)
   and `evidence["lr_order"]` present (the LR order used — R3 provenance;
   the 0-var vocabulary "no_integration_required" of SubTropica.wl:17794
   is accepted). `mode=:face_sum` (check_divergences OFF per ct) is the
   Phase-B1 farm and REFUSES here.
2. Letter-table union across every contribution (per-response Wm/Wp
   indices remapped into the unioned `LaurentSeries.letters`).
3. Per-face per-order contribution sums; pole-padding alignment across
   faces ([`align_laurent`](@ref)); face sum.
4. Laurent convolution with the symbolic prefactor series
   `(pmin, pcoeffs)` (Γ-prefactor expansion is eps_expand.jl's
   `loggamma_series`, RT-vii), truncated to `requested_order`, with a
   LOUD truncation guard (R2): any needed-but-missing input order throws
   [`LaurentTruncationError`](@ref) — vanishing orders must be explicit
   zeros.
5. Optional `mzv=true`: HF `apply_mzv_reductions` pass on every output
   coefficient ([`mzv_reduce`](@ref)); table provenance recorded.
6. Provenance stamped into `LaurentSeries.evidence`: HF binary sha256,
   both upstream SHAs (reference/PROVENANCE.md — note the adac2f72 /
   ead8c6e skew), LR order, check_divergences state, wall time,
   timestamps. Provenance tag is `:analytic` (the `:numerically_matched_
   sign` tag is exclusively the R4 continuation module's, Phase C+).
"""
function assemble_faces(faces::Vector{FaceResult},
                        prefactor_series::Tuple{Int,Vector{HlogExpr}},
                        requested_order::Int;
                        evidence::Dict=Dict{String,Any}(),
                        mode::Symbol=:standalone,
                        cfg::HFConfig=HFConfig(),
                        transport=_default_transport,
                        mzv::Bool=false)
    t0 = time()
    # --- 1. policy gates (LOUD; see docstring) ----------------------------
    if mode === :face_sum
        throw(DivergenceCheckViolation(
            "mode=:face_sum (check_divergences OFF per counter-term) is " *
            "the Phase-B1 farm (DESIGN.md C11) — not built; REFUSING"))
    end
    mode === :standalone || throw(ArgumentError(
        "assemble_faces: unknown mode $(mode) (Phase A: :standalone only)"))
    cd = get(evidence, "check_divergences", nothing)
    cd === true || throw(DivergenceCheckViolation(
        "standalone assembly requires evidence[\"check_divergences\"] === " *
        "true; got $(repr(cd)). check_divergences is HARD-ON in Phase A " *
        "(CONTRACTS.md §a policy; upstream silently returns 0 on divergent " *
        "input with it off — confirmed on the binary). Re-run the HF " *
        "calls with the flag armed; do NOT stamp the flag post hoc."))
    haskey(evidence, "lr_order") || throw(AssembleProvenanceError(
        "evidence[\"lr_order\"] missing — the LR order used is REQUIRED " *
        "provenance on every assembled series (DESIGN.md R3; use " *
        "\"no_integration_required\" for 0-var inputs, SubTropica.wl:17794)"))
    isempty(faces) && throw(ArgumentError("assemble_faces: no faces"))
    isempty(prefactor_series[2]) &&
        throw(ArgumentError("assemble_faces: empty prefactor series"))

    # --- 2. letter union ---------------------------------------------------
    tables = Vector{AlgLetter}[]
    for f in faces, blk in f.orders, (_, tab) in blk
        push!(tables, tab)
    end
    (letters, remaps) = union_letters(tables)

    # --- 3. per-face per-order sums, then aligned face fold ----------------
    ti = 0
    pieces = Tuple{Int,Vector{HlogExpr}}[]
    for f in faces
        coeffs = HlogExpr[]
        for blk in f.orders
            acc = _a3_hlog_zero()
            for (expr, _) in blk
                ti += 1
                acc = _a3_hadd(acc, remap_letters(expr, remaps[ti]))
            end
            push!(coeffs, acc)
        end
        push!(pieces, (f.minorder, coeffs))
    end
    (gmin, padded) = align_laurent(pieces)
    nI = length(padded[1])
    icoeffs = [reduce(_a3_hadd, (p[k] for p in padded); init=_a3_hlog_zero())
               for k in 1:nI]

    # --- 4. prefactor convolution + truncation guard (R2) ------------------
    (pmin, pcoeffs) = prefactor_series
    out_min = gmin + pmin
    requested_order >= out_min || throw(LaurentTruncationError(
        "requested_order=$(requested_order) is below the leading order " *
        "$(out_min) — nothing to return; check the pole offset"))
    imax = gmin + nI - 1
    pmax = pmin + length(pcoeffs) - 1
    need_i = requested_order - pmin
    need_p = requested_order - gmin
    need_i <= imax || throw(LaurentTruncationError(
        "convolution to ε^$(requested_order) needs integrand orders up to " *
        "ε^$(need_i) but inputs stop at ε^$(imax) — supply the missing " *
        "orders (explicit zeros if they vanish); never assume them " *
        "(truncated-Laurent pole guard, DESIGN.md R2)"))
    need_p <= pmax || throw(LaurentTruncationError(
        "convolution to ε^$(requested_order) needs prefactor orders up to " *
        "ε^$(need_p) but the prefactor series stops at ε^$(pmax) — extend " *
        "loggamma_series (eps_expand.jl) to the required order"))
    coeffs = HlogExpr[]
    for o in out_min:requested_order
        acc = _a3_hlog_zero()
        for k in gmin:imax
            j = o - k
            (pmin <= j <= pmax) || continue
            acc = _a3_hadd(acc, _a3_hmul(icoeffs[k - gmin + 1],
                                   pcoeffs[j - pmin + 1]))
        end
        push!(coeffs, acc)
    end

    # --- 5. optional MZV reduction pass ------------------------------------
    ev = Dict{String,Any}(pairs(evidence))
    if mzv
        coeffs = [mzv_reduce(c; cfg=cfg, transport=transport) for c in coeffs]
        ev["mzv_table_path"] = cfg.mzv_data_path
        ev["mzv_table_sha256"] = _file_sha256(cfg.mzv_data_path)
    end

    # --- 6. provenance stamp ------------------------------------------------
    (elf, binsha) = _hf_bin_pin(cfg)
    (wlsha, bingit) = _upstream_pins()
    ev["hf_bin_path"] = elf
    ev["hf_bin_sha256"] = binsha
    ev["upstream_wl_sha"] = wlsha          # vendored spec pin (adac2f72)
    ev["upstream_bin_git"] = bingit        # binary build pin (ead8c6e) — SKEW
    ev["mode"] = String(mode)
    ev["requested_order"] = requested_order
    ev["assembled_at"] = Libc.strftime("%Y-%m-%dT%H:%M:%S", time())
    ev["wall_s_assemble"] = time() - t0
    return LaurentSeries{HlogExpr}(out_min, coeffs, letters, :analytic, ev)
end

"""
    assemble(minorder::Int,
             order_results::Vector{Vector{Tuple{HlogExpr,Vector{AlgLetter}}}},
             prefactor_series::Tuple{Int,Vector{HlogExpr}},
             requested_order::Int; evidence=Dict{String,Any}(), kwargs...)
        -> LaurentSeries{HlogExpr}

Phase-A degenerate entry point: the SINGLE-face special case of
[`assemble_faces`](@ref) (the pinned C12 signature). `order_results[k]` =
contributions at ε-order `minorder + k - 1`. All policy gates, provenance
requirements and kwargs are those of `assemble_faces`.
"""
function assemble(minorder::Int, order_results, prefactor_series,
                  requested_order::Int; evidence::Dict=Dict{String,Any}(),
                  kwargs...)
    orders = [Tuple{HlogExpr,Vector{AlgLetter}}[
                  (expr, collect(AlgLetter, tab)) for (expr, tab) in blk]
              for blk in order_results]
    return assemble_faces([FaceResult(minorder, orders)],
                          prefactor_series, requested_order;
                          evidence=evidence, kwargs...)
end

# ---------------------------------------------------------------------------
# MZV reduction pass (HF apply_mzv_reductions; table provenance pinned)
# ---------------------------------------------------------------------------

_zeta_token(idx::Vector{Int}) =
    "mzv_" * join((i < 0 ? "m$(-i)" : string(i) for i in idx), "_")
function _token_zeta(tok::AbstractString)
    startswith(tok, "mzv_") || throw(MzvReduceError("not an mzv token: $tok"))
    idx = Int[]
    for p in split(tok[5:end], '_')
        if startswith(p, "m")
            push!(idx, -parse(Int, p[2:end]))
        else
            push!(idx, parse(Int, p))
        end
    end
    return idx
end

# Symbol pool of the reduction table (lhs tokens + every symbol appearing
# in a rhs + the basis). REQUIRED: apply_mzv_reductions builds its
# PolyCtx from the request's "vars" ONLY (no build_mzv_var_list —
# hyperflint_cli_main.cpp:2117ff), so a reduction whose OUTPUT symbol is
# missing from vars dies with an in-band Poly-parse error (live-verified
# observed: f=mzv_2_1 with vars=[mzv_2_1] → error "Poly: parse error:
# (-2*mzv_3)"; adding mzv_3 → result "-2*mzv_3").
const _MZV_POOL_MEMO = Dict{String,Vector{String}}()
function _mzv_symbol_pool(path::AbstractString)
    get!(_MZV_POOL_MEMO, String(path)) do
        isfile(path) || throw(MzvReduceError(
            "mzv table not found at $(path) (cfg.mzv_data_path)"))
        doc = JSON.parsefile(path)
        syms = Set{String}()
        for r in get(doc, "reductions", [])
            push!(syms, String(r["lhs"]))
            for m in eachmatch(r"[A-Za-z][A-Za-z0-9_]*", String(r["rhs"]))
                push!(syms, m.match)
            end
        end
        for b in get(doc, "basis", [])
            push!(syms, String(b))
        end
        sort!(collect(syms))
    end
end

const _MZV_RED_MEMO = Dict{Tuple{String,String},HlogExpr}()

function _mzv_monomial_reduce(f::String, tokens::Vector{String};
                              cfg::HFConfig, transport)
    key = (String(cfg.mzv_data_path), f)
    haskey(_MZV_RED_MEMO, key) && return _MZV_RED_MEMO[key]
    pool = _mzv_symbol_pool(cfg.mzv_data_path)
    vars = sort(unique(vcat(pool, tokens)))
    resp = transport("apply_mzv_reductions",
                     Dict{String,Any}("f" => f, "vars" => vars,
                                      "mzv_data_path" => cfg.mzv_data_path);
                     cfg=cfg)
    haskey(resp, "error") && throw(MzvReduceError(
        "apply_mzv_reductions in-band error on f=$(f): $(resp["error"])"))
    get(resp, "failed", false) === true &&
        throw(MzvReduceError("apply_mzv_reductions failed on f=$(f)"))
    haskey(resp, "result") || throw(MzvReduceError(
        "apply_mzv_reductions response missing \"result\" on f=$(f)"))
    red = _parse_reduced(String(resp["result"]))
    _MZV_RED_MEMO[key] = red
    return red
end

"""
    mzv_reduce(e::HlogExpr; cfg=HFConfig(), transport=_default_transport)
        -> HlogExpr

MZV reduction pass: rewrite every HZeta/Log2 monomial through HF's
`apply_mzv_reductions` (table =
`cfg.mzv_data_path`, default the pinned
the engine mzv_reductions.json data file (SUBTROPICA_MZV_DATA); basis =
Log2 + 9 irreducible MZVs). Monomials are reduced one at a time (the op
substitutes lhs tokens INSIDE monomials, live-verified: mzv_2_1^2 →
4*mzv_3^2) and memoized per (table, monomial). Non-MZV atoms (logs,
Hlogs, kinematics) ride along untouched. Table provenance
(path + sha256) must be recorded by the caller — `assemble_faces(mzv=true)`
does this automatically (R7)."""
function mzv_reduce(e::HlogExpr; cfg::HFConfig=HFConfig(),
                    transport=_default_transport)
    out = HlogExpr(HTerm[], e.vars)
    for t in e.terms
        mztoks = Pair{String,Int}[]
        rest = Pair{HAtom,Int}[]
        for (a, p) in t.atoms
            # Log2 arrives as HLogA("2") from hlogexpr.jl's parse_coef (§d:
            # Log2/Log[n] → HLogA) or as HConst(:Log2) (types.jl named
            # constant, e.g. from loggamma_series) — both are the table's
            # "Log2" symbol.
            if a isa HZeta || (a isa HConst && a.name === :Log2) ||
               (a isa HLogA && a.arg == "2")
                p > 0 || throw(MzvReduceError(
                    "negative/zero power $(p) on MZV atom — refusing"))
                tok = a isa HZeta ? _zeta_token(a.idx) : "Log2"
                push!(mztoks, tok => p)
            else
                push!(rest, a => p)
            end
        end
        base = HlogExpr([HTerm(t.coef, rest)], e.vars)
        if isempty(mztoks)
            out = _a3_hadd(out, base)
        else
            f = join((p == 1 ? tok : "$(tok)^$(p)" for (tok, p) in mztoks), "*")
            red = _mzv_monomial_reduce(f, String[tok for (tok, _) in mztoks];
                                       cfg=cfg, transport=transport)
            out = _a3_hadd(out, _a3_hmul(base, red))
        end
    end
    return out
end

# --- restricted parser for reduction output --------------------------------
# apply_mzv_reductions emits plain polynomials over the token ctx with
# rational coefficients (live-pinned formats: "-2*mzv_3",
# "1/2*Log2*mzv_2 - 1/4*mzv_3", "4*mzv_3^2", "mzv_2 - 2*mzv_3", "mzv_2^2";
# the table rhs additionally uses an outer "-(...)" wrap). Anything
# outside this grammar throws LOUDLY.
# STATUS: hlogexpr.jl's full parse_coef HAS landed (hlogexpr.jl `parse_coef`;
# it has no "not implemented" throw), so on this tree the probe below always
# answers true and _parse_reduced routes to parse_coef. _parse_mzv_poly is
# therefore NOT on the production path; it is retained as the pinned-grammar
# reference parser exercised directly by test/test_assemble.jl. The probe is
# kept as a runtime check of the engine, not as a live fallback switch.
const _A1_PARSE = Ref{Int}(-1)
function _hlog_parse_available()
    if _A1_PARSE[] < 0
        _A1_PARSE[] = try
            parse_coef("1", Symbol[])
            1
        catch err
            if err isa ErrorException && occursin("not implemented", err.msg)
                0
            else
                rethrow()
            end
        end
    end
    return _A1_PARSE[] == 1
end

function _parse_reduced(s::String)
    _hlog_parse_available() && return parse_coef(s, Symbol[])
    return _parse_mzv_poly(s)
end

_strip_outer_parens(s::AbstractString) = begin
    t = strip(s)
    while startswith(t, "(") && endswith(t, ")")
        depth = 0
        wraps = true
        for (i, c) in enumerate(t)
            c == '(' && (depth += 1)
            c == ')' && (depth -= 1)
            depth == 0 && i < lastindex(t) && (wraps = false; break)
        end
        wraps || break
        t = strip(t[2:end-1])
    end
    t
end

function _mzv_atom_of(tok::AbstractString)
    # Log2 → HConst(:Log2): the canonical Log2 atom (design decision,
    # CONTRACTS.md §d), matching hlogexpr.jl's parse_coef (_norm_atom fold) so both
    # routes of _parse_reduced agree.
    tok == "Log2" && return HConst(:Log2)
    tok == "Pi" && return HConst(:Pi)
    tok == "EulerGamma" && return HConst(:EulerGamma)
    tok == "Catalan" && return HConst(:Catalan)
    startswith(tok, "mzv_") && return HZeta(_token_zeta(tok))
    throw(MzvReduceError("unknown symbol \"$(tok)\" in reduction output — " *
                         "local restricted parser refuses (hlogexpr.jl's " *
                         "parse_coef handles the full grammar)"))
end

"Restricted local parse: rational-coefficient polynomial in
Log2/Pi/EulerGamma/Catalan/mzv_* tokens (see comment block above)."
function _parse_mzv_poly(s0::String)
    s = _strip_outer_parens(replace(s0, " " => ""))
    isempty(s) && throw(MzvReduceError("empty reduction output"))
    s == "0" && return HlogExpr()
    # handle a global "-(...)" wrap (table-rhs style)
    if startswith(s, "-(") && endswith(s, ")") &&
       _strip_outer_parens(s[2:end]) != s[2:end] # balanced wrap check
        inner = _parse_mzv_poly(String(s[2:end]))
        return _a3_hscale(inner, -1 // 1)
    end
    # split into signed terms at depth-0 +/- (not the leading sign)
    terms = HTerm[]
    depth = 0
    i = firstindex(s)
    seg_start = i
    seg_sign = 1 // 1
    if s[i] == '+'
        seg_start = nextind(s, i)
    elseif s[i] == '-'
        seg_sign = -1 // 1
        seg_start = nextind(s, i)
    end
    pos = seg_start
    function _push_seg(stop)
        seg = s[seg_start:stop]
        isempty(seg) && throw(MzvReduceError("dangling sign in \"$(s0)\""))
        coef = seg_sign * (1 // 1)
        atoms = Pair{HAtom,Int}[]
        for fac in split(_strip_outer_parens(seg), '*')
            r = _a3_tryrat(fac)
            if r !== nothing
                coef *= r
                continue
            end
            m = match(r"^([A-Za-z][A-Za-z0-9_]*)(?:\^([0-9]+))?$", fac)
            m === nothing && throw(MzvReduceError(
                "unparseable factor \"$(fac)\" in \"$(s0)\" — grammar " *
                "outside the pinned reduction-output forms; refusing"))
            p = m.captures[2] === nothing ? 1 : parse(Int, m.captures[2])
            push!(atoms, _mzv_atom_of(m.captures[1]) => p)
        end
        coef == 0 || push!(terms, HTerm(_a3_ratstr(Rational{BigInt}(coef)),
                                        _a3_norm_atoms(atoms)))
    end
    while pos <= lastindex(s)
        c = s[pos]
        if c == '('
            depth += 1
        elseif c == ')'
            depth -= 1
        elseif depth == 0 && (c == '+' || c == '-') && pos > seg_start
            prevc = s[prevind(s, pos)]
            if !(prevc in ('*', '/', '^', '(', '+', '-'))
                _push_seg(prevind(s, pos))
                seg_sign = (c == '-') ? -1 // 1 : 1 // 1
                seg_start = nextind(s, pos)
            end
        end
        pos = nextind(s, pos)
    end
    _push_seg(lastindex(s))
    return HlogExpr(terms, Symbol[])
end

# ---------------------------------------------------------------------------
# ginsh export (verification harness side-door)
# ---------------------------------------------------------------------------

# CONVENTION PINS (why the export is deliberately narrow):
#  * Depth-1 tokens are standard zeta: mzv_2 = ζ(2), mzv_3 = ζ(3) verified
#    to 60 digits vs independent mpmath quadrature (the T2/T3 smoke integrals).
#  * Depth-≥2 / alternating tokens are NOT naive ζ(a,b,...): live probe
#    the pinned binary gives mzv_2_1 → -2*mzv_3, whereas standard ζ(2,1) = +ζ(3)
#    (Euler) — i.e. HF's token value convention differs from GiNaC's
#    zeta({...}) by more than index order. Exporting them unverified would
#    fabricate constants. Run mzv_reduce first (kills every REDUCIBLE
#    token); irreducible basis tokens (mzv_1_m3, mzv_1_m5, mzv_1_1_m3,
#    mzv_1_1_1_m3, mzv_3_5) REFUSE until the zero_inf_period ↔ ginac
#    cross-evaluation (gate A(v), hlogexpr.jl/verify.jl) pins their values.
#  * Log[expr] ≡ Hlog[expr,{0}] upstream (CONTRACTS.md §d) matches
#    GiNaC's G(0;z) = log(z), so Hlog[z,{l...}] exports as G({l...},z).
#  * EulerGamma is `Euler` in ginsh; Log2 exports as log(2).

const _GINSH_CONST = Dict{Symbol,String}(
    :Pi => "Pi", :I => "I", :EulerGamma => "Euler",
    :Catalan => "Catalan", :Log2 => "log(2)")

function _ginsh_atom(a::HAtom)
    if a isa HConst
        haskey(_GINSH_CONST, a.name) && return _GINSH_CONST[a.name]
        throw(GinshExportError("named constant $(a.name) has no pinned " *
                               "ginsh form"))
    elseif a isa HZeta
        if length(a.idx) == 1 && a.idx[1] >= 2
            return "zeta($(a.idx[1]))"
        end
        throw(GinshExportError(
            "HZeta($(a.idx)): depth-≥2/alternating MZV tokens do NOT " *
            "follow the naive zeta({...}) convention (live pin: mzv_2_1 = " *
            "-2*mzv_3 vs standard ζ(2,1)=+ζ(3)). Run mzv_reduce first; " *
            "irreducible basis tokens await the zero_inf_period↔ginac " *
            "convention pin (gate A(v)). REFUSING rather than fabricating."))
    elseif a isa HLogA
        return "log(" * a.arg * ")"
    elseif a isa HHlog
        return "G({" * join(a.word, ",") * "}," * a.z * ")"
    elseif a isa HPeriod
        throw(GinshExportError("unresolved ZeroInfPeriod key " *
            "$(a.word) — resolve via resolve_periods (hlogexpr.jl) before export"))
    elseif a isa HAlgLetter
        throw(GinshExportError("unresolved algebraic letter " *
            "$(a.kind)_$(a.idx) — back-substitute via HF " *
            "simplify_with_vieta/back_substitute (RT-v) before export"))
    elseif a isa HDelta
        throw(GinshExportError("delta[$(a.var)] contour residue — " *
            "Euclidean-only export (R4); a delta atom here is a bug"))
    end
    throw(GinshExportError("unknown atom type $(typeof(a))"))
end

# Coef strings are Mma-grammar rational functions of kinematics
# (CONTRACTS.md §d) — that fragment is ginsh-compatible verbatim. Guard
# against un-parsed function tokens leaking through.
function _ginsh_coef_guard(coef::String)
    occursin(r"[A-Za-z][A-Za-z0-9_]*\[", coef) && throw(GinshExportError(
        "coef string \"$(coef)\" contains Mma bracket syntax — parse it " *
        "into atoms (hlogexpr.jl parse_coef) before ginsh export"))
    for bad in ("mzv_", "Wm_", "Wp_", "sqrt_disc_", "delta")
        occursin(bad, coef) && throw(GinshExportError(
            "coef string \"$(coef)\" contains unparsed token \"$(bad)\" — " *
            "refusing"))
    end
    return coef
end

"""
    to_ginsh(e::HlogExpr) -> String
    to_ginsh(L::LaurentSeries{HlogExpr}) -> Vector{Tuple{Int,String}}

Export to ginac-evaluable (ginsh) strings for the verification harness
(verify.jl / GPLEval-ginac oracle side). Convention pins and the deliberate
REFUSALS (depth-≥2 MZV tokens, unresolved periods/letters, deltas) are
documented in the source comment block above — this exporter never
fabricates a constant's value. Zero coefficients export as "0". The
series form returns one `(ε_order, string)` pair per coefficient."""
function to_ginsh(e::HlogExpr)
    _a3_hlog_iszero(e) && return "0"
    parts = String[]
    for t in e.terms
        fac = String["(" * _ginsh_coef_guard(t.coef) * ")"]
        for (a, p) in t.atoms
            s = _ginsh_atom(a)
            push!(fac, p == 1 ? s : (p > 0 ? "$(s)^$(p)" : "$(s)^($(p))"))
        end
        push!(parts, join(fac, "*"))
    end
    return join(parts, "+")
end

to_ginsh(L::LaurentSeries{HlogExpr}) =
    Tuple{Int,String}[(L.minorder + k - 1, to_ginsh(L.coeffs[k]))
                      for k in eachindex(L.coeffs)]

"""
    ginsh_series_string(L::LaurentSeries{HlogExpr}; eps="eps") -> String

Single ε-polynomial ginsh string `Σ (coeff_o)*eps^(o)` (convenience for
the harness; per-order comparison should use `to_ginsh(L)`)."""
ginsh_series_string(L::LaurentSeries{HlogExpr}; eps::String="eps") =
    join(("(" * s * ")*" * eps * "^(" * string(o) * ")"
          for (o, s) in to_ginsh(L)), "+")
