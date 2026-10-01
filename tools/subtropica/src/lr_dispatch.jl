# lr_dispatch.jl — Thin wrappers on HF LR-order ops (C10).
#
# PINS (CONTRACTS.md §a, live-verified against the v1.2.8 binary):
#  * find_lr_orders takes MULTI-GROUP `groups:[[..]]` per-addend semantics —
#    MANDATORY for counter-term sums (single-group `polys` on a ct-sum =>
#    NOLR by construction). This module NEVER emits the `polys` form.
#  * There is NO standalone verify_order op: order verification is the
#    OPTIONAL "verify_order":[...] FIELD of the find_lr_orders request
#    (handlers.cpp:741, VERIFY-ORDER mode). In verify mode the
#    envelope's best_order/score/nolr fields are INERT (handlers.cpp:858)
#    — the single binding bit is "order_is_lr", plus the
#    verify_blocking_* fields explaining a false.
#  * Search is time-budgeted branch-and-bound (HF_LR_TIME_BUDGET_S):
#    best_order may vary run-to-run, the ANSWER is order-independent by
#    theorem (DESIGN.md R3). Hence: persist best_order per integrand hash
#    and REPLAY it via the verify_order field on resume; any order change
#    is recorded so the ≥30d numeric gate re-verifies.
#  * Upstream persists bestOrder per face directory (SubTropica.wl:17821,
#    `Put[bestOrder, file2 <> "/bestOrder.m"]`); we mirror that with
#    best_order.json under a caller-supplied run dir (CONTRACTS.md §b
#    layout). Upstream's search loop aborts on NOLR (SubTropica.wl:17804
#    `If[bestOrder === NOLR, ... Abort[]]`) — but the INTEGRATION layer
#    has a confirmed silent-return-0 footgun on the no-contribution path
#    ($NoAlgebraicRootsContributions analog, SubTropica.wl:2553; confirmed
#    on the binary). Our port is the OPPOSITE: a missing LR
#    order is a LOUD, TYPED error (NoLROrderError), never a silent zero
#    and never a warning.
#
# Transport: every wrapper takes `transport(op::String, fields::Dict;
# cfg::HFConfig) -> Dict` (default = `hf_call` from hf_bridge.jl, CONTRACTS.md §a).
# Tests inject canned-JSON mocks so no long HF calls are needed.

# ---------------------------------------------------------------------------
# Typed errors — LOUD, never silent (see header).
# ---------------------------------------------------------------------------

"""
    NoLROrderError

Thrown when HF `find_lr_orders` reports `nolr:true`: NO linearly-reducible
integration order exists for the given groups. This error is deliberately
LOUD and typed — the upstream Mathematica integration layer silently
returns 0 in exactly this situation (SubTropica.wl:2553
\$NoAlgebraicRootsContributions analog; confirmed on the binary) — our
port REFUSES instead. If the input is a counter-term SUM, check that each
addend was passed as its OWN group (multi-group semantics): a ct-sum
flattened into one group is NOLR by construction (CONTRACTS.md §a).
"""
struct NoLROrderError <: Exception
    ngroups::Int
    nxvars::Int
    xvars::Vector{String}
    detail::String
end
function Base.showerror(io::IO, e::NoLROrderError)
    print(io, "NoLROrderError: NO linearly-reducible order exists ",
          "(ngroups=", e.ngroups, ", nxvars=", e.nxvars,
          ", xvars=", e.xvars, "). REFUSING — upstream Mathematica ",
          "silently returns 0 here; subtropica never does. ",
          "If this integrand is a counter-term SUM, pass one addend per ",
          "group (multi-group `groups` semantics); a flat single group on ",
          "a ct-sum is NOLR by construction. ", e.detail)
end

"""
    LRDispatchError

Malformed or in-band-error HF response on an LR-order op (CONTRACTS.md §a
taxonomy: `{"error":..}` / `{"failed":true}` / missing required fields —
never gate on exit code, rc is 0 even on handled failure).
"""
struct LRDispatchError <: Exception
    op::String
    msg::String
end
Base.showerror(io::IO, e::LRDispatchError) =
    print(io, "LRDispatchError[", e.op, "]: ", e.msg)

"""
    LRReplayError

Persisted best_order state is corrupt, belongs to a different integrand
hash, or fails verify replay under `on_stale=:error`.
"""
struct LRReplayError <: Exception
    path::String
    msg::String
end
Base.showerror(io::IO, e::LRReplayError) =
    print(io, "LRReplayError[", e.path, "]: ", e.msg)

# ---------------------------------------------------------------------------
# Transport default + envelope checks
# ---------------------------------------------------------------------------

# Default transport = the HF bridge (hf_bridge.jl). Tests inject mocks.
_default_transport(op::String, fields::Dict; cfg::HFConfig=HFConfig()) =
    hf_call(op, fields; cfg=cfg)

_strvec(v) = String[string(x) for x in v]

# In-band error taxonomy (print-loud-rc0-gates footgun class): gate on the
# presence of the needed fields AND absence of error/failed keys.
# (AbstractDict: JSON.jl 1.x parses responses into JSON.Object, not Dict.)
function _lr_check_envelope(resp::AbstractDict, needed::Vector{String})
    haskey(resp, "error") &&
        throw(LRDispatchError("find_lr_orders",
            "in-band error: " * string(resp["error"])))
    get(resp, "failed", false) === true &&
        throw(LRDispatchError("find_lr_orders", "in-band {\"failed\":true}"))
    for k in needed
        haskey(resp, k) || throw(LRDispatchError("find_lr_orders",
            "response missing required field \"" * k *
            "\" — malformed envelope (or pre-v1.2.8 binary; run hf_gate)"))
    end
    return nothing
end

function _lr_validate_groups(groups, xvars)
    isempty(groups) && throw(LRDispatchError("find_lr_orders",
        "empty `groups` — caller bug"))
    for (i, g) in enumerate(groups)
        isempty(g) && throw(LRDispatchError("find_lr_orders",
            "group $i is empty — caller bug"))
    end
    isempty(xvars) && throw(LRDispatchError("find_lr_orders",
        "empty `xvars` — nothing to order (the .wl vocabulary for the "
        * "0-var case is \"no_integration_required\", SubTropica.wl:17794; "
        * "handle it in the caller, do not dispatch)"))
    return nothing
end

function _lr_base_fields(groups, xvars, coeff_vars;
                         algebraic_letters::Bool, carry_discharge::Bool)
    fields = Dict{String,Any}(
        # ALWAYS the multi-group form; never "polys" (CONTRACTS.md §a).
        "groups"     => [_strvec(g) for g in groups],
        "xvars"      => _strvec(xvars),
        "coeff_vars" => _strvec(coeff_vars),
    )
    # Only arm optional flags when requested: the Strict envelope stays
    # byte-identical to the baseline (handlers.cpp regression gate).
    algebraic_letters && (fields["algebraic_letters"] = true)
    carry_discharge   && (fields["carry_discharge"] = true)
    return fields
end

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

"""
    find_lr_order(groups::Vector{<:AbstractVector}, xvars, coeff_vars;
                  cfg=HFConfig(), transport=_default_transport,
                  algebraic_letters=false, carry_discharge=false)
        -> (order::Vector{Symbol}, meta::Dict{String,Any})

`find_lr_orders` search wrapper (request/response pinned live;
schema_version 2). ALWAYS uses the MULTI-GROUP `groups` request form — a
counter-term sum must arrive as one addend per group; the single-group
`polys` form is never emitted (CONTRACTS.md §a). `nolr:true` throws the
LOUD typed [`NoLROrderError`](@ref) — never a silent skip (the upstream
package silently returns 0 in this spot; we are the opposite by design).

`meta` carries score/strategy/timing/hf_version/schema_version, the group
census (nXVars/nGroups/nPolys) and `root_polys` when
`algebraic_letters=true` (deg-2 polys HF will mint Wm/Wp letters for).
"""
function find_lr_order(groups, xvars, coeff_vars;
                       cfg::HFConfig=HFConfig(),
                       transport=_default_transport,
                       algebraic_letters::Bool=false,
                       carry_discharge::Bool=false)
    _lr_validate_groups(groups, xvars)
    fields = _lr_base_fields(groups, xvars, coeff_vars;
                             algebraic_letters=algebraic_letters,
                             carry_discharge=carry_discharge)
    resp = transport("find_lr_orders", fields; cfg=cfg)
    _lr_check_envelope(resp, ["nolr", "best_order"])
    if resp["nolr"] === true
        throw(NoLROrderError(length(groups), length(xvars), _strvec(xvars),
            "strategy=" * string(get(resp, "strategy", "?")) *
            ", timing_compute_s=" * string(get(resp, "timing_compute_s", "?"))))
    end
    order = Symbol.(_strvec(resp["best_order"]))
    isempty(order) && throw(LRDispatchError("find_lr_orders",
        "nolr:false but empty best_order — malformed response"))
    meta = Dict{String,Any}(
        "score"            => get(resp, "score", nothing),
        "strategy"         => get(resp, "strategy", nothing),
        "timing_compute_s" => get(resp, "timing_compute_s", nothing),
        "hf_version"       => get(resp, "hf_version", nothing),
        "schema_version"   => get(resp, "schema_version", nothing),
        "nXVars"           => get(resp, "nXVars", nothing),
        "nGroups"          => get(resp, "nGroups", nothing),
        "nPolys"           => get(resp, "nPolys", nothing),
    )
    haskey(resp, "root_polys") && (meta["root_polys"] = resp["root_polys"])
    haskey(resp, "carried_polys") && (meta["carried_polys"] = resp["carried_polys"])
    return (order, meta)
end

"""
    verify_lr_order(order, groups, xvars, coeff_vars;
                    cfg=HFConfig(), transport=_default_transport,
                    details=nothing) -> Bool

Replay a persisted best_order via the `"verify_order"` FIELD of
`find_lr_orders` (handlers.cpp:741; O(n) st_fubini_lr calls, no search —
this is NOT a standalone op, the source map's op list was wrong). Returns
the single binding bit `order_is_lr`. The envelope's best_order/score/nolr
are INERT in this mode and are deliberately ignored (handlers.cpp:858).
Pass a `details::Dict` to receive the `verify_blocking_step/degree/letter`
fields explaining a `false`. A response WITHOUT the `order_is_lr` key
throws (wrong-binary guard — run `hf_gate`).
"""
function verify_lr_order(order, groups, xvars, coeff_vars;
                         cfg::HFConfig=HFConfig(),
                         transport=_default_transport,
                         details::Union{Nothing,Dict}=nothing)
    _lr_validate_groups(groups, xvars)
    isempty(order) && throw(LRDispatchError("find_lr_orders",
        "empty `order` passed to verify_lr_order — caller bug"))
    fields = _lr_base_fields(groups, xvars, coeff_vars;
                             algebraic_letters=false, carry_discharge=false)
    fields["verify_order"] = _strvec(order)
    resp = transport("find_lr_orders", fields; cfg=cfg)
    _lr_check_envelope(resp, ["order_is_lr"])
    ok = resp["order_is_lr"] === true
    if details !== nothing
        for k in ("verify_malformed", "verify_blocking_step",
                  "verify_blocking_degree", "verify_forbidden_dep",
                  "verify_blocking_letter", "timing_compute_s")
            haskey(resp, k) && (details[k] = resp[k])
        end
    end
    if !ok
        @warn "subtropica verify_lr_order: order is NOT linearly reducible" order =
            _strvec(order) blocking_step = get(resp, "verify_blocking_step", nothing) blocking_letter =
            get(resp, "verify_blocking_letter", nothing)
    end
    return ok
end

# ---------------------------------------------------------------------------
# Persistence: best_order.json per integrand hash (R3)
# ---------------------------------------------------------------------------

# sha256 via a shellout (Project.toml is frozen, so no SHA stdlib dep;
# the shellout is exact and this is a cold-path helper). GNU coreutils
# ships `sha256sum`; stock macOS ships only `shasum` (Digest::SHA) — probe
# for either, and refuse with the typed error when neither is on PATH.
# Both tools emit identical digests, so persisted integrand keys are stable.
_sha256_cmd() = Sys.which("sha256sum") !== nothing ? `sha256sum` :
                Sys.which("shasum") !== nothing    ? `shasum -a 256` : nothing

function _sha256_hex(data::Vector{UInt8})
    cmd = _sha256_cmd()
    cmd === nothing && throw(LRDispatchError("sha256",
        "no sha256 tool on PATH (need `sha256sum` or `shasum`) — cannot " *
        "form integrand hash; refusing to persist under a weak key"))
    out = try
        read(pipeline(IOBuffer(data), cmd), String)
    catch err
        throw(LRDispatchError("sha256",
            "$(cmd.exec[1]) shellout failed ($(typeof(err))) — cannot form " *
            "integrand hash; refusing to persist under a weak key"))
    end
    toks = split(out)
    (!isempty(toks) && occursin(r"^[0-9a-f]{64}$", toks[1])) ||
        throw(LRDispatchError("sha256",
            "unexpected $(cmd.exec[1]) output: " * out))
    return String(toks[1])
end
_sha256_hex(s::AbstractString) = _sha256_hex(Vector{UInt8}(codeunits(s)))

"""
    integrand_hash(groups, xvars, coeff_vars) -> String

Deterministic content key (full sha256 hex) of the LR-search input triple.
The JSON payload is an ARRAY (not a Dict) so key order can never wobble.
Sensitive to group membership AND to xvar order (a different variable list
is a different search problem).
"""
integrand_hash(groups, xvars, coeff_vars) =
    _sha256_hex(JSON.json(Any["subtropica_lr_v1",
                              [_strvec(g) for g in groups],
                              _strvec(xvars), _strvec(coeff_vars)]))

"""
    best_order_path(run_dir, groups, xvars, coeff_vars) -> String

Canonical persistence location, mirroring the upstream per-face directory
layout (SubTropica.wl:17821 `bestOrder.m`; CONTRACTS.md §b):
`<run_dir>/integrands/<first-16-hex-of-hash>/best_order.json`.
"""
best_order_path(run_dir::AbstractString, groups, xvars, coeff_vars) =
    joinpath(run_dir, "integrands",
             first(integrand_hash(groups, xvars, coeff_vars), 16),
             "best_order.json")

"""
    persist_best_order(path, order, meta) -> path

Atomically write `best_order.json` (schema subtropica_best_order_v1): temp
file + `mv` so a killed writer never leaves a torn JSON for the resume
path to trip on.
"""
function persist_best_order(path::AbstractString, order, meta::Dict)
    mkpath(dirname(path))
    doc = Dict{String,Any}(
        "schema"     => "subtropica_best_order_v1",
        "best_order" => _strvec(order),
        "meta"       => meta,
        "saved_unix" => time(),
        "saved_at"   => Libc.strftime("%Y-%m-%dT%H:%M:%S", time()),
    )
    tmp = path * ".tmp.$(getpid())"
    open(tmp, "w") do io
        JSON.print(io, doc)
    end
    mv(tmp, path; force=true)
    return String(path)
end

"""
    load_best_order(path) -> (order::Vector{Symbol}, meta::Dict{String,Any})

Read a persisted best_order.json; schema mismatch or missing fields throw
[`LRReplayError`](@ref) (corrupt state must never be silently re-searched
away without the caller knowing — `ensure_lr_order` handles the policy).
"""
function load_best_order(path::AbstractString)
    isfile(path) || throw(LRReplayError(String(path), "no such file"))
    doc = try
        JSON.parsefile(path)
    catch err
        throw(LRReplayError(String(path), "unparseable JSON: $(err)"))
    end
    get(doc, "schema", "") == "subtropica_best_order_v1" ||
        throw(LRReplayError(String(path),
            "schema mismatch: " * string(get(doc, "schema", "MISSING"))))
    haskey(doc, "best_order") ||
        throw(LRReplayError(String(path), "missing best_order field"))
    order = Symbol.(_strvec(doc["best_order"]))
    meta = Dict{String,Any}(get(doc, "meta", Dict{String,Any}()))
    return (order, meta)
end

"""
    ensure_lr_order(groups, xvars, coeff_vars; run_dir,
                    cfg=HFConfig(), transport=_default_transport,
                    on_stale=:research,
                    algebraic_letters=false, carry_discharge=false)
        -> (order::Vector{Symbol}, meta::Dict{String,Any})

The resume driver (R3): if a best_order.json exists for this integrand
hash under `run_dir`, REPLAY it through the verify_order field (O(n), no
search) and return it (`meta["replayed"]=true`). If replay fails —
possible after an upstream pin change — behavior follows `on_stale`:
`:research` (default) runs a fresh search, re-persists, and records
`meta["order_changed_from"]` so the ≥30d numeric gate re-verifies;
`:error` throws [`LRReplayError`](@ref). No persisted state ⇒ fresh
search + persist.
"""
function ensure_lr_order(groups, xvars, coeff_vars;
                         run_dir::AbstractString,
                         cfg::HFConfig=HFConfig(),
                         transport=_default_transport,
                         on_stale::Symbol=:research,
                         algebraic_letters::Bool=false,
                         carry_discharge::Bool=false)
    on_stale in (:research, :error) ||
        throw(ArgumentError("ensure_lr_order: on_stale must be :research or :error"))
    h = integrand_hash(groups, xvars, coeff_vars)
    path = joinpath(run_dir, "integrands", first(h, 16), "best_order.json")
    stale_order = nothing
    if isfile(path)
        (order, meta) = load_best_order(path)
        stored = get(meta, "integrand_hash", "")
        stored == h || throw(LRReplayError(path,
            "stored integrand_hash $(stored) != computed $(h) — corrupt " *
            "or colliding state; refusing to replay"))
        if verify_lr_order(order, groups, xvars, coeff_vars;
                           cfg=cfg, transport=transport)
            meta = copy(meta)
            meta["replayed"] = true
            return (order, meta)
        end
        on_stale === :error && throw(LRReplayError(path,
            "persisted order $(order) failed verify_order replay " *
            "(NOT LR under the current binary) and on_stale=:error"))
        stale_order = order
    end
    (order, meta) = find_lr_order(groups, xvars, coeff_vars;
                                  cfg=cfg, transport=transport,
                                  algebraic_letters=algebraic_letters,
                                  carry_discharge=carry_discharge)
    meta["integrand_hash"] = h
    meta["replayed"] = false
    if stale_order !== nothing
        # R3: order changed => intermediate alphabet may change (answer
        # does not); the ≥30d gate MUST re-verify. Record loudly.
        meta["order_changed_from"] = _strvec(stale_order)
        @warn "subtropica ensure_lr_order: persisted order failed replay; re-searched" old =
            _strvec(stale_order) new = _strvec(order)
    end
    persist_best_order(path, order, meta)
    return (order, meta)
end

# ---------------------------------------------------------------------------
# Phase B wrappers: find_lr_orders_scan, factor_table, gauge derivation
# (CONTRACTS.md §a op contracts; request/response pinned against
# reference/hyperflint_bridge_handlers.cpp find_lr_orders_scan/factor_table).
# ---------------------------------------------------------------------------

# In-band error + required-field gate shared by the Phase B ops (never gate
# on rc — same taxonomy as _lr_check_envelope, parametrized by op name).
function _phaseb_check_envelope(op::String, resp::AbstractDict,
                                needed::Vector{String})
    haskey(resp, "error") &&
        throw(LRDispatchError(op, "in-band error: " * string(resp["error"])))
    get(resp, "failed", false) === true &&
        throw(LRDispatchError(op, "in-band {\"failed\":true}"))
    for k in needed
        haskey(resp, k) || throw(LRDispatchError(op,
            "response missing required field \"" * k *
            "\" — malformed envelope (or wrong-version binary; run hf_gate)"))
    end
    return nothing
end

"""
    find_lr_orders_scan(groups, xvars, coeff_vars, exps;
                        keep_rule="Strict", euler_filter=false,
                        max_orders=nothing, cfg=HFConfig(),
                        transport=_default_transport)
        -> (orders::Vector{Dict{String,Any}}, meta::Dict{String,Any})

`find_lr_orders_scan` wrapper (Cheng-Wu GAUGE SCAN, CONTRACTS.md §a): the
PRE-GAUGE engine over the RAW factor groups of the ungauged n-variable
system. `exps` = per-poly twist pairs `[[[a,b],...],...]`, shape matching
`groups` — REQUIRED, the projectivity gate (Σ a_i·d_i == −n, Σ b_i·d_i == 0)
is load-bearing. `keep_rule` is `"Strict"` (default) or `"FindRoots"`
(NECESSARY-only/speculative carried-sqrt tier).

Response semantics are DISTINGUISHED, never conflated (handler contract):
  * `projective:false` ⇒ the INPUT failed the projectivity gate — a typed
    [`LRDispatchError`](@ref) (caller bug in the twists; NOT a NOLR verdict).
  * `projective:true` + empty `orders` ⇒ no admissible (gauge, order) pair
    exists under the keep rule — the scan's NOLR analog, thrown as the LOUD
    typed [`NoLROrderError`](@ref) (upstream silently returns 0 in the
    integration layer; this port never does).
Returned `orders` entries carry `order::Vector{Symbol}`, `gauge::Symbol`,
plus the engine's score/carried_sqrts/kin_sqrts/terminal_quads fields,
score-ascending as emitted. `meta` carries projective/truncated/timing/
version/census fields; `truncated:true` means max_orders capped the list.
"""
function find_lr_orders_scan(groups, xvars, coeff_vars, exps;
                             keep_rule::AbstractString="Strict",
                             euler_filter::Bool=false,
                             max_orders::Union{Nothing,Int}=nothing,
                             cfg::HFConfig=HFConfig(),
                             transport=_default_transport)
    _lr_validate_groups(groups, xvars)
    keep_rule in ("Strict", "FindRoots") ||
        throw(LRDispatchError("find_lr_orders_scan",
            "unknown keep_rule \"$keep_rule\" (allowed: Strict, FindRoots)"))
    # exps shape must equal groups shape — mirror the engine's loud check
    # client-side so mock and live transports refuse identically.
    length(exps) == length(groups) ||
        throw(LRDispatchError("find_lr_orders_scan",
            "exps group count ($(length(exps))) != groups ($(length(groups)))"))
    for (g, (ge, gp)) in enumerate(zip(exps, groups))
        length(ge) == length(gp) ||
            throw(LRDispatchError("find_lr_orders_scan",
                "exps[$g] length $(length(ge)) != group size $(length(gp))"))
        for pr in ge
            (length(pr) == 2 && all(x -> x isa Integer, pr)) ||
                throw(LRDispatchError("find_lr_orders_scan",
                    "exps[$g] entry $(pr) is not an integer [a,b] twist pair"))
        end
    end
    fields = Dict{String,Any}(
        "groups"     => [_strvec(g) for g in groups],
        "xvars"      => _strvec(xvars),
        "coeff_vars" => _strvec(coeff_vars),
        "exps"       => [[Int[Int(pr[1]), Int(pr[2])] for pr in ge] for ge in exps],
    )
    # Only arm non-default knobs (baseline-envelope rule, as in _lr_base_fields).
    keep_rule != "Strict" && (fields["keep_rule"] = String(keep_rule))
    euler_filter && (fields["euler_filter"] = true)
    max_orders !== nothing && (fields["max_orders"] = max_orders)
    resp = transport("find_lr_orders_scan", fields; cfg=cfg)
    _phaseb_check_envelope("find_lr_orders_scan", resp,
                           ["projective", "truncated", "orders"])
    truncated = resp["truncated"] === true
    if resp["projective"] !== true
        throw(LRDispatchError("find_lr_orders_scan",
            "projectivity gate FAILED — the twist exponents do not satisfy " *
            "sum(a_i*d_i) == -n and sum(b_i*d_i) == 0; this is an INPUT " *
            "error (fix `exps`), NOT a NOLR verdict"))
    end
    raw = resp["orders"]
    isempty(raw) && throw(NoLROrderError(length(groups), length(xvars),
        _strvec(xvars),
        "gauge scan: projective input but NO admissible (gauge, order) pair " *
        "under keep_rule=$(keep_rule)" *
        (truncated ? " (truncated:true — raise max_orders)" : "") *
        ", timing_compute_s=" * string(get(resp, "timing_compute_s", "?"))))
    orders = Vector{Dict{String,Any}}()
    for (i, so) in enumerate(raw)
        (so isa AbstractDict && haskey(so, "order") && haskey(so, "gauge")) ||
            throw(LRDispatchError("find_lr_orders_scan",
                "orders[$i] lacks order/gauge — malformed response"))
        entry = Dict{String,Any}(
            "order" => Symbol.(_strvec(so["order"])),
            "gauge" => Symbol(String(so["gauge"])),
        )
        for k in ("score", "carried_sqrts", "kin_sqrts", "terminal_quads")
            haskey(so, k) && (entry[k] = so[k])
        end
        push!(orders, entry)
    end
    meta = Dict{String,Any}(
        "projective"       => true,
        "truncated"        => truncated,
        "keep_rule"        => String(keep_rule),
        "timing_compute_s" => get(resp, "timing_compute_s", nothing),
        "hf_version"       => get(resp, "hf_version", nothing),
        "schema_version"   => get(resp, "schema_version", nothing),
        "nXVars"           => get(resp, "nXVars", nothing),
        "nGroups"          => get(resp, "nGroups", nothing),
    )
    return (orders, meta)
end

"""
    factor_table(groups, xvars, coeff_vars, order;
                 algebraic_letters=false, max_pairs=nothing,
                 max_singletons=nothing, max_response_mb=nothing,
                 cfg=HFConfig(), transport=_default_transport)
        -> Dict{String,Any}

`factor_table` wrapper (CONTRACTS.md §a): single-chain replay along a
supplied LR `order` — factor-prediction table of monic pair differences and
per-letter coefficient/discriminant factor lists. ALWAYS the multi-group
`groups` request form (never `polys`); `order` must be a permutation of
`xvars` (checked client-side, mirroring the engine's loud refusal). Engine
guards are loud errors, never truncation — an in-band `error` (including
the max_pairs/max_singletons/max_response_mb caps) throws
[`LRDispatchError`](@ref). Returns the response table Dict: interned
`polys`, `stages`, `pairs`, `singletons`, `stats`, plus the order and
version fields.
"""
function factor_table(groups, xvars, coeff_vars, order;
                      algebraic_letters::Bool=false,
                      max_pairs::Union{Nothing,Int}=nothing,
                      max_singletons::Union{Nothing,Int}=nothing,
                      max_response_mb::Union{Nothing,Int}=nothing,
                      cfg::HFConfig=HFConfig(),
                      transport=_default_transport)
    _lr_validate_groups(groups, xvars)
    isempty(order) && throw(LRDispatchError("factor_table",
        "empty `order` — caller bug"))
    sort(_strvec(order)) == sort(_strvec(xvars)) ||
        throw(LRDispatchError("factor_table",
            "order $(_strvec(order)) is not a permutation of xvars " *
            "$(_strvec(xvars))"))
    fields = _lr_base_fields(groups, xvars, coeff_vars;
                             algebraic_letters=algebraic_letters,
                             carry_discharge=false)
    fields["order"] = _strvec(order)
    max_pairs !== nothing && (fields["max_pairs"] = max_pairs)
    max_singletons !== nothing && (fields["max_singletons"] = max_singletons)
    max_response_mb !== nothing && (fields["max_response_mb"] = max_response_mb)
    resp = transport("factor_table", fields; cfg=cfg)
    _phaseb_check_envelope("factor_table", resp,
                           ["polys", "stages", "pairs", "singletons", "stats"])
    table = Dict{String,Any}(
        "order"      => Symbol.(_strvec(get(resp, "order", _strvec(order)))),
        "polys"      => resp["polys"],
        "stages"     => resp["stages"],
        "pairs"      => resp["pairs"],
        "singletons" => resp["singletons"],
        "stats"      => resp["stats"],
        "hf_version"     => get(resp, "hf_version", nothing),
        "schema_version" => get(resp, "schema_version", nothing),
    )
    return table
end

"""
    derive_gauge_from_homogeneous_lr(letters, xvars;
                                     coeff_vars=String[], cfg=HFConfig(),
                                     transport=_default_transport)
        -> (gauge::Symbol, order::Vector{Symbol}, meta::Dict{String,Any})

Port of upstream `stDeriveGaugeFromHomogeneousLR` (SubTropica.wl:17097;
consumed by STIntegrate's `"GaugeStrategy" -> "Derive"`): for an eps-free
PROJECTIVE integrand, derive a Cheng-Wu gauge + LR integration order from
ONE rational LR-order search on the homogeneous singularity `letters` over
ALL `xvars` — the LAST variable of the full best order is the gauge (set
to 1), the remaining variables are the integration order. `letters` are
the distinct variable-bearing irreducible denominator factors (upstream's
`FactorList[Denominator[Together[...]]]` + DeleteCases step stays with the
caller — the graph front-end owns integrand factorization in this tree);
they are dispatched as ONE group, matching upstream's `{letters}` call.

Divergences from upstream, both deliberate: upstream returns `\$Failed` on
NOLR or on `DeriveTimeBudget` expiry — this port throws the LOUD typed
[`NoLROrderError`](@ref) instead (silence is never an option), and the
time budget rides the engine env (`HF_LR_TIME_BUDGET_S`), not an option.
Rational (FindRoots-off) tier only, as upstream. A best order that is not
a full permutation of `xvars` is a malformed response
([`LRDispatchError`](@ref)) — the gauge split needs the complete order.
"""
function derive_gauge_from_homogeneous_lr(letters, xvars;
                                          coeff_vars=String[],
                                          cfg::HFConfig=HFConfig(),
                                          transport=_default_transport)
    isempty(letters) && throw(LRDispatchError("find_lr_orders",
        "derive_gauge_from_homogeneous_lr: empty letter list — no " *
        "variable-bearing denominator factors (upstream returns \$Failed " *
        "here, wl:17102; this port refuses loudly instead)"))
    length(xvars) >= 2 || throw(LRDispatchError("find_lr_orders",
        "derive_gauge_from_homogeneous_lr needs >= 2 xvars (one gauge + at " *
        "least one integration variable); got $(length(xvars))"))
    (order, meta) = find_lr_order([_strvec(letters)], xvars, coeff_vars;
                                  cfg=cfg, transport=transport)
    sort(_strvec(order)) == sort(_strvec(xvars)) ||
        throw(LRDispatchError("find_lr_orders",
            "derive_gauge_from_homogeneous_lr: best_order $(_strvec(order)) " *
            "is not a full order over xvars $(_strvec(xvars)) — malformed " *
            "response (the gauge split needs the complete order)"))
    gauge = order[end]
    meta = copy(meta)
    meta["gauge_strategy"] = "derive_from_homogeneous_lr"
    meta["gauge"] = String(gauge)
    return (gauge, order[1:end-1], meta)
end
