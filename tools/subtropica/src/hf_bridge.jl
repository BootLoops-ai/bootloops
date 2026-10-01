# hf_bridge.jl — The HyperFLINT boundary (CONTRACTS.md §a).
# CLI-first transport: ONE flat JSON object on stdin of
# `hyperflint.sh eval-json`, parse stdout's LAST line only, errors are
# IN-BAND (never gate on rc — a handled failure still exits 0; control
# outcomes confirmed on the binary, see CONTRACTS.md §a).
# Every spawn under `ulimit -v 32505856` (the convention used throughout this
# package — built into the command line). Kill only by PID (we signal exactly
# our own child).
# narrow_ctx_insufficient => retry ONCE with HF_NARROW_CTX=0, no recursion
# (mirrors the engine's own guidance, handlers.cpp:2183-2201).
#
# Error taxonomy source: reference/hyperflint_bridge_handlers.cpp:2183-2231
# ({"narrow_ctx_insufficient":true} / {"divergent":true,"reason":..} /
# {"failed":true} / error_json_op) + main.cpp per-op {"error":..} paths.
# Version gate: DESIGN.md RT-vi (stamp + unary-minus canary; source builds
# mis-stamp 1.2.0 — stamp alone MISLEADS).

# ---------------------------------------------------------------------------
# Typed exceptions — in-band failures raise, NEVER return defaults.
# ---------------------------------------------------------------------------

abstract type HFException <: Exception end

"Transport-level failure: no parseable JSON on stdout's last line."
struct HFTransportError <: HFException
    op::String
    msg::String
    rc::Union{Int,Nothing}
    stderr_tail::String
end
Base.showerror(io::IO, e::HFTransportError) =
    print(io, "HFTransportError[$(e.op)]: ", e.msg, " (rc=", e.rc,
          "; stderr tail: ", e.stderr_tail, ")")

"In-band {\"error\":..} response (includes \"unknown op\") or a missing result field."
struct HFError <: HFException
    op::String
    msg::String
    response::Union{Dict,Nothing}
end
Base.showerror(io::IO, e::HFError) = print(io, "HFError[$(e.op)]: ", e.msg)

"In-band {\"failed\":true}: integration step failed (escalation is the CALLER's job — FindRoots cascade, C11)."
struct HFFailed <: HFException
    op::String
    request_summary::String
    response::Dict
end
Base.showerror(io::IO, e::HFFailed) =
    print(io, "HFFailed[$(e.op)]: {\"failed\":true} — integration step failed; ",
          "escalate (algebraic_letters tier / carry_discharge), request: ",
          e.request_summary)

"In-band {\"divergent\":true}: FATAL on standalone calls (CONTRACTS.md §a policy)."
struct HFDivergent <: HFException
    op::String
    reason::String
    response::Dict
end
Base.showerror(io::IO, e::HFDivergent) =
    print(io, "HFDivergent[$(e.op)]: ", e.reason)

"narrow_ctx_insufficient persisted after the single flag-off retry."
struct HFNarrowCtx <: HFException
    op::String
    response::Dict
end
Base.showerror(io::IO, e::HFNarrowCtx) =
    print(io, "HFNarrowCtx[$(e.op)]: narrow_ctx_insufficient with HF_NARROW_CTX=0 ",
          "— engine-side bug surface, refusing")

"Wall-clock timeout: the child was killed BY PID (SIGTERM, then SIGKILL)."
struct HFTimeout <: HFException
    op::String
    timeout_s::Float64
    pid::Int32
end
Base.showerror(io::IO, e::HFTimeout) =
    print(io, "HFTimeout[$(e.op)]: exceeded ", e.timeout_s, "s; child pid ",
          e.pid, " killed by PID")

"Version gate failure (stamp and/or behavioral canary — DESIGN.md RT-vi)."
struct HFGateError <: HFException
    msg::String
end
Base.showerror(io::IO, e::HFGateError) = print(io, "HFGateError: ", e.msg)

# ---------------------------------------------------------------------------
# Subprocess plumbing (shared with the ginac hooks in hlogexpr.jl)
# ---------------------------------------------------------------------------

"""
    _spawn_capture(argv, input; timeout_s=3600, env=nothing)
        -> (stdout, stderr, rc, timed_out, pid)

Run `argv` with optional stdin string, capturing stdout/stderr with async
readers (no pipe-full deadlock). On timeout: SIGTERM to exactly our child
PID, 5 s grace, then SIGKILL — kill only by PID. `env` is a
FULL replacement environment Dict (callers merge over `ENV`).
"""
function _spawn_capture(argv::Vector{String}, input::Union{Nothing,AbstractString};
                        timeout_s::Real=3600,
                        env::Union{Nothing,Dict{String,String}}=nothing)
    cmd = Cmd(argv)
    env === nothing || (cmd = setenv(cmd, env))
    out = Pipe(); err = Pipe()
    stdin_src = input === nothing ? devnull : IOBuffer(String(input))
    proc = run(pipeline(cmd, stdin=stdin_src, stdout=out, stderr=err), wait=false)
    pid = try Int32(getpid(proc)) catch; Int32(-1) end   # capture while alive
    close(out.in); close(err.in)
    out_task = @async read(out, String)
    err_task = @async read(err, String)
    t0 = time()
    while process_running(proc) && (time() - t0) < timeout_s
        sleep(0.02)
    end
    timed_out = process_running(proc)
    if timed_out
        kill(proc)                       # SIGTERM, our child's PID only
        t1 = time()
        while process_running(proc) && (time() - t1) < 5.0
            sleep(0.05)
        end
        process_running(proc) && kill(proc, 9)   # SIGKILL escalation
    end
    wait(proc)
    return (stdout=fetch(out_task), stderr=fetch(err_task),
            rc=proc.exitcode, timed_out=timed_out, pid=pid)
end

_bash_ulimit_argv(bin::String, args::Vector{String}) =
    vcat(["bash", "-c", "ulimit -v 32505856; exec \"\$0\" \"\$@\"", bin], args)

function _hf_env(cfg::HFConfig; narrow_ctx::Union{Nothing,Bool}=nothing)
    env = Dict{String,String}(k => v for (k, v) in ENV)
    env["OMP_NUM_THREADS"] = string(cfg.omp_num_threads)
    env["HF_MAX_THREADS_PER_CALL"] = string(cfg.max_threads_per_call)
    env["HF_LR_TIME_BUDGET_S"] = string(cfg.lr_time_budget_s)
    # HF_NARROW_CTX: only set explicitly when armed (or when force-disabling
    # on the retry path); leave the engine default otherwise (the environment
    # the control outcomes in CONTRACTS.md §a were established under).
    nc = narrow_ctx === nothing ? cfg.narrow_ctx : narrow_ctx
    if nc
        env["HF_NARROW_CTX"] = "1"
    elseif narrow_ctx === false
        env["HF_NARROW_CTX"] = "0"
    end
    return env
end

const _STDERR_TAIL = 400
_tail(s::String) = length(s) <= _STDERR_TAIL ? s : s[prevind(s, end, _STDERR_TAIL - 1):end]

# One raw exchange: CLI subprocess (default) or in-process C-ABI ccall
# (OPT-IN, transport=:cabi — see the C-ABI section below). Returns the
# parsed response object (CLI: LAST stdout line; ABI: the returned string).
function _hf_transport(req::AbstractDict, cfg::HFConfig;
                       timeout_s::Real, narrow_ctx::Union{Nothing,Bool}=nothing,
                       transport::Symbol=cfg.transport)
    op = String(get(req, "op", "?"))
    transport in (:cli, :cabi) || throw(HFError(op,
        "unknown transport $(repr(transport)) — :cli or :cabi", nothing))
    if transport === :cabi && haskey(_CABI_OPS, op)
        return _cabi_transport(req, cfg; narrow_ctx=narrow_ctx)
    end
    # :cli — or :cabi on an op with NO ABI entry point (parse_expr / eval /
    # period / vieta ops are CLI-only upstream, C-ABI gap G1): documented
    # CLI fallback (same handler bodies — the byte-identity class), never a
    # silent number.
    isfile(cfg.bin) || throw(HFTransportError(op,
        "HyperFLINT wrapper not found at $(cfg.bin) (set SUBTROPICA_HF_BIN)", nothing, ""))
    body = JSON.json(req)
    r = _spawn_capture(_bash_ulimit_argv(cfg.bin, ["eval-json"]), body;
                       timeout_s=timeout_s, env=_hf_env(cfg; narrow_ctx=narrow_ctx))
    r.timed_out && throw(HFTimeout(op, Float64(timeout_s), r.pid))
    lines = [l for l in split(r.stdout, '\n') if !isempty(strip(l))]
    isempty(lines) && throw(HFTransportError(op, "empty stdout", r.rc, _tail(r.stderr)))
    last_line = String(lines[end])          # CONTRACTS.md §a: LAST line only
    resp = try
        JSON.parse(last_line)
    catch e
        throw(HFTransportError(op, "unparseable last stdout line: " *
            _tail(last_line) * " ($(sprint(showerror, e)))", r.rc, _tail(r.stderr)))
    end
    resp isa AbstractDict ||   # JSON.jl v1 yields JSON.Object, an AbstractDict
        throw(HFTransportError(op, "last line is not a JSON object",
                               r.rc, _tail(r.stderr)))
    return resp
end

# ---------------------------------------------------------------------------
# C-ABI in-process transport (OPT-IN).
# References: deps/BUILD.md (the patched, statically isolated lib, its
# 8-symbol ABI surface and ownership rule, and the F2 patch: a per-call
# PolyCtx hygiene fix; we observed heap corruption when reusing the unpatched
# library in-process).
# LAWS (all measured, none negotiable):
#   * OPT-IN only: hf_call(...; transport=:cabi), HFConfig(transport=:cabi),
#     or SUBTROPICA_HF_TRANSPORT=cabi. The CLI stays the DEFAULT transport AND
#     the VALIDATION ORACLE — it is process-isolated and immune to the F2
#     bug class; production farms keep 5-10% of calls dual-run
#     (_cabi_call_raw vs _cli_call_raw, _cabi_normalize compare) as a
#     drift tripwire .
#   * ONLY the F2-PATCHED static-isolated lib is loadable (deps/BUILD.md):
#     the dlopen gate is strict version 1.2.8 + response envelope +
#     unary-minus canary + an F2 pair-repro probe. An unpatched library is
#     unsafe to reuse in-process; use the CLI transport (one process per call).
#   * ONE dlopen handle per PROCESS (RTLD_LOCAL — hf internals stay private
#     from Nemo/FLINT_jll). The farm unit is a WORKER PROCESS: setsid from
#     MAIN, ulimit -v per worker, kill by PID only.
#   * THREADS SHARING THE HANDLE ARE FORBIDDEN — in our tests 8 threads on
#     one handle abort inside FLINT, with or without the F2 patch (the cached
#     MZV-table rebuild is not mutex-protected). _CABI_LOCK serializes every ccall as a
#     belt-and-suspenders guard: accidental Threads use degrades to
#     safe-but-serial, never to a crash. Farm parallelism = PROCESSES.
#   * timeout_s is NOT enforced on this transport (an in-process call
#     cannot be killed safely); runaway calls die WITH the worker process.
#   * Ops with no ABI entry point (parse_expr/eval/period/vieta/... —
#     C-ABI gap G1) fall back to the CLI even under :cabi (documented in
#     _hf_transport; same upstream handler bodies).
# Measured why: 24-70x per op on small ops, 176-4497x on medium integrator
# ops vs the CLI's ~6-8 ms spawn+MZV-load floor (BENCH.md; re-measured
# post-F2-fix, speedup intact). NOT a wall-mover for genuine engine compute
# (a real LR-search wall is unchanged in-process — measured on the in-house
# divergent target).
# ---------------------------------------------------------------------------

# op name -> ABI symbol (the FULL public ABI surface, deps/BUILD.md).
const _CABI_OPS = Dict{String,Symbol}(
    "partial_fractions"   => :hf_partial_fractions,
    "linear_factors"      => :hf_linear_factors,
    "find_lr_orders"      => :hf_find_lr_orders,
    "find_lr_orders_scan" => :hf_find_lr_orders_scan,
    "factor_table"        => :hf_factor_table,
    "hyperflint"          => :hf_hyperflint_sym,
)

const _CABI_LOCK = ReentrantLock()               # serialize ALL ccalls
const _CABI_HANDLES = Dict{String,Ptr{Nothing}}() # lib path -> GATED handle

# One ccall exchange on an already-gated handle. Caller holds _CABI_LOCK.
# Ownership contract (deps/BUILD.md): response char* freed via hf_free_string
# ONLY; NULL = catastrophic alloc failure (fatal, typed).
function _cabi_exchange(h::Ptr{Nothing}, sym::Symbol, body::AbstractString)::String
    p = ccall(Libdl.dlsym(h, sym), Ptr{UInt8}, (Cstring,), String(body))
    p == C_NULL && throw(HFTransportError(String(sym),
        "C-ABI returned NULL — catastrophic allocation failure", nothing, ""))
    s = unsafe_string(p)
    ccall(Libdl.dlsym(h, :hf_free_string), Cvoid, (Ptr{UInt8},), p)
    return s
end

"""
    _cabi_handle(cfg::HFConfig) -> Ptr{Nothing}

dlopen-once-per-process handle for `cfg.cabi_lib`, GATED at first open
(memoized on success only): (1) `hf_version_string()` == "1.2.8" STRICT
(source builds mis-stamp 1.2.0 — DESIGN.md RT-vi footgun); (2) response
envelope `"hf_version":"1.2.8"` (double check, binding law);
(3) unary-minus canary through the IN-PROCESS parser (pfrac of -x^2);
(4) F2 pair-repro probe (three pfrac calls in this process, repeat must be
byte-identical and error-free — an unpatched lib fails/corrupts here).
"""
function _cabi_handle(cfg::HFConfig)::Ptr{Nothing}
    lock(_CABI_LOCK) do
        get!(_CABI_HANDLES, cfg.cabi_lib) do
            isfile(cfg.cabi_lib) || throw(HFGateError(
                "C-ABI lib not found at $(cfg.cabi_lib) — set SUBTROPICA_HF_CABI_LIB " *
                "or build it: tools/subtropica/deps/build_cabi.sh (deps/BUILD.md)"))
            h = Libdl.dlopen(cfg.cabi_lib)   # RTLD_LOCAL default
            try
                ver = unsafe_string(ccall(Libdl.dlsym(h, :hf_version_string),
                                          Cstring, ()))
                ver == "1.2.8" || throw(HFGateError(
                    "C-ABI hf_version_string() = $(repr(ver)), need exactly " *
                    "\"1.2.8\" (mis-stamp footgun: rebuild with -DHF_VERSION=1.2.8)"))
                can = _cabi_exchange(h, :hf_partial_fractions,
                    """{"op":"partial_fractions","f":"-x^2","vars":["x"],"var":"x"}""")
                occursin("\"hf_version\":\"1.2.8\"", can) || throw(HFGateError(
                    "C-ABI response envelope lacks \"hf_version\":\"1.2.8\" — " *
                    "stamp/envelope mismatch, refusing (canary resp: " * _tail(can) * ")"))
                (occursin("-x^2", can) && !occursin("\"error\"", can)) ||
                    throw(HFGateError("C-ABI unary-minus canary FAILED " *
                        "(pfrac(-x^2) resp: " * _tail(can) * ") — pre-precedence-fix " *
                        "build or broken parser; REFUSING"))
                p1 = """{"op":"partial_fractions","f":"1/(x*(1 + x))","vars":["x"],"var":"x"}"""
                p2 = """{"op":"partial_fractions","f":"1/(x*(1 + x)*(x - y))","vars":["x","y"],"var":"x"}"""
                r1 = _cabi_exchange(h, :hf_partial_fractions, p1)
                r2 = _cabi_exchange(h, :hf_partial_fractions, p2)
                r3 = _cabi_exchange(h, :hf_partial_fractions, p2)
                (!occursin("\"error\"", r1) && !occursin("\"error\"", r2) &&
                 r2 == r3) || throw(HFGateError(
                    "C-ABI F2 pair-repro probe FAILED (in-process reuse " *
                    "corrupts state) — this lib is NOT the F2-patched build; " *
                    "see the F2 fix report"))
            catch
                Libdl.dlclose(h)
                rethrow()
            end
            h
        end
    end
end

"""
    _cabi_call_raw(op::String, req_json; cfg=HFConfig(), narrow_ctx=nothing)
        -> String

RAW request/response string exchange over the in-process C ABI (the
byte-identity leg for tests and the production dual-run tripwire).
Env knobs (OMP_NUM_THREADS / HF_MAX_THREADS_PER_CALL / HF_LR_TIME_BUDGET_S
/ HF_NARROW_CTX) are set via `withenv` around the ccall — the engine reads
them PER REQUEST even in-process (measured: the budget knob tripped
in-process on the in-house divergent target). Serialized behind _CABI_LOCK.
"""
function _cabi_call_raw(op::String, req_json::AbstractString;
                        cfg::HFConfig=HFConfig(),
                        narrow_ctx::Union{Nothing,Bool}=nothing)
    sym = get(_CABI_OPS, op, nothing)
    sym === nothing && throw(HFError(op,
        "op '$op' has no C-ABI entry point (C-ABI gap G1: CLI-only " *
        "upstream) — use the CLI transport", nothing))
    h = _cabi_handle(cfg)
    envpairs = Pair{String,String}["OMP_NUM_THREADS" => string(cfg.omp_num_threads),
        "HF_MAX_THREADS_PER_CALL" => string(cfg.max_threads_per_call),
        "HF_LR_TIME_BUDGET_S" => string(cfg.lr_time_budget_s)]
    nc = narrow_ctx === nothing ? cfg.narrow_ctx : narrow_ctx
    if nc
        push!(envpairs, "HF_NARROW_CTX" => "1")
    elseif narrow_ctx === false
        push!(envpairs, "HF_NARROW_CTX" => "0")
    end
    lock(_CABI_LOCK) do
        withenv(envpairs...) do
            _cabi_exchange(h, sym, req_json)
        end
    end
end

# Parsed-object leg used by _hf_transport under transport=:cabi.
function _cabi_transport(req::AbstractDict, cfg::HFConfig;
                         narrow_ctx::Union{Nothing,Bool}=nothing)
    op = String(get(req, "op", "?"))
    s = _cabi_call_raw(op, JSON.json(req); cfg=cfg, narrow_ctx=narrow_ctx)
    resp = try
        JSON.parse(s)
    catch e
        throw(HFTransportError(op, "unparseable C-ABI response: " * _tail(s) *
            " ($(sprint(showerror, e)))", nothing, ""))
    end
    resp isa AbstractDict || throw(HFTransportError(op,
        "C-ABI response is not a JSON object", nothing, ""))
    hv = get(resp, "hf_version", nothing)
    hv == "1.2.8" || throw(HFGateError(
        "C-ABI response envelope hf_version = $(repr(hv)), need \"1.2.8\" " *
        "(per-response double check, binding law)"))
    return resp
end

# --- CLI-vs-ABI byte-identity helpers (tests + dual-run drift tripwire) ---
# Normalization = the repo's own snapshot-ctest strip: the ABI splices
# ,"schema_version":N,"hf_version":"..." after "op" (the CLI body lacks it
# for pfrac/lf/hyperflint), and timing_compute_s is wall-clock.
_strip_envelope(s::AbstractString) =
    replace(s, r",\"schema_version\":\d+,\"hf_version\":\"[^\"]*\"" => "")
_strip_timing(s::AbstractString) =
    replace(s, r",\"timing_compute_s\":[0-9.eE+-]+" => "")
_cabi_normalize(s::AbstractString) = _strip_timing(_strip_envelope(strip(s)))

"RAW CLI leg: last stdout line for a literal request string (oracle side)."
function _cli_call_raw(req_json::AbstractString; cfg::HFConfig=HFConfig(),
                       timeout_s::Real=cfg.timeout_s)
    isfile(cfg.bin) || throw(HFTransportError("raw",
        "HyperFLINT wrapper not found at $(cfg.bin)", nothing, ""))
    r = _spawn_capture(_bash_ulimit_argv(cfg.bin, ["eval-json"]), req_json;
                       timeout_s=timeout_s, env=_hf_env(cfg))
    r.timed_out && throw(HFTimeout("raw", Float64(timeout_s), r.pid))
    lines = [l for l in split(r.stdout, '\n') if !isempty(strip(l))]
    isempty(lines) && throw(HFTransportError("raw", "empty stdout", r.rc,
                                             _tail(r.stderr)))
    return String(lines[end])
end

# ---------------------------------------------------------------------------
# Request hygiene: the C++ JSON parser is REGEX-BASED AND FLAT
# (CONTRACTS.md §a) — refuse anything it cannot see BEFORE sending.
# ---------------------------------------------------------------------------

# Fields whose elements may be flat objects (allocations/regulator replay
# arrays — main.cpp:2634-2812, handlers.cpp evaluate_periods).
const _OBJECT_ARRAY_FIELDS = ("allocations", "regulator")

function _assert_flat_array(op::String, k::String, v::AbstractVector, depth::Int,
                            allow_obj::Bool)
    depth > 4 && throw(HFError(op, "field '$k': array nesting deeper than the " *
        "documented shapes — the HF regex parser cannot read it", nothing))
    for el in v
        if el isa Union{AbstractString,Number,Bool}
        elseif el isa AbstractVector
            _assert_flat_array(op, k, el, depth + 1, false)
        elseif el isa AbstractDict && allow_obj
            for (dk, dv) in el
                if dv isa Union{AbstractString,Number,Bool}
                elseif dv isa AbstractVector
                    _assert_flat_array(op, "$k.$dk", dv, depth + 2, false)
                else
                    throw(HFError(op, "field '$k' entry key '$dk': nested object — flat JSON only (CONTRACTS.md §a)", nothing))
                end
            end
        else
            throw(HFError(op, "field '$k': element of type $(typeof(el)) — flat JSON only (CONTRACTS.md §a)", nothing))
        end
    end
end

function _assert_flat(op::String, req::AbstractDict)
    for (k, v) in req
        ks = String(k)
        if v isa Union{AbstractString,Number,Bool}
        elseif v isa AbstractVector
            _assert_flat_array(op, ks, v, 1, ks in _OBJECT_ARRAY_FIELDS)
        else
            throw(HFError(op, "field '$ks' of type $(typeof(v)): the HF JSON " *
                "parser is regex-based and FLAT — never nest objects " *
                "(CONTRACTS.md §a)", nothing))
        end
    end
end

# Ops that consume the MZV reduction table: inject mzv_data_path explicitly
# (CONTRACTS.md §a "always pass it explicitly").
const _MZV_OPS = Set(["hyperflint", "evaluate_periods", "zero_inf_period",
                      "zero_one_period", "apply_mzv_reductions"])

# Caller-side misuse (missing required request fields) is a programmer
# error and raises a plain ErrorException BEFORE any spawn (the frozen
# scaffold test in runtests.jl relies on this); engine-side failures are
# the typed HFException taxonomy. Entries are any-of alternatives.
const _REQUIRED_FIELDS = Dict{String,Vector{Vector{String}}}(
    "parse_expr" => [["expr"]],
    "eval" => [["a"]],
    "hyperflint" => [["expr", "f", "wordlist"]],   # priority expr > f > wordlist
    "evaluate_periods" => [["regulator"]],
    "zero_inf_period" => [["word"]],
    "zero_one_period" => [["word"]],
    "simplify_with_vieta" => [["expr"]],
    "back_substitute" => [["expr"]],
    "combine_wm_wp_ratios" => [["expr"]],
    "find_lr_orders" => [["groups", "polys"]],
    "factor" => [["expr"]],
)

function _check_required(op::String, req::AbstractDict)
    alts = get(_REQUIRED_FIELDS, op, nothing)
    alts === nothing && return
    for group in alts
        any(k -> haskey(req, k), group) ||
            error("SubTropica.hf_call: op '$op' requires one of the fields " *
                  "$(group) (CONTRACTS.md §a) — refusing before spawn")
    end
end

# Known response payload key per op (checked present after the error
# taxonomy; ops not listed get no payload check beyond taxonomy).
const _RESULT_KEYS = Dict{String,String}(
    "eval" => "result", "hyperflint" => "result", "evaluate_periods" => "result",
    "zero_inf_period" => "result", "zero_one_period" => "result",
    "simplify_with_vieta" => "result", "back_substitute" => "result",
    "combine_wm_wp_ratios" => "result", "apply_mzv_reductions" => "result",
    "series_expansion" => "result", "partial_fractions" => "result",
    "parse_expr" => "canonical",                       # main.cpp:2593-2632
    "find_lr_orders" => "best_order",
    "find_lr_orders_scan" => "orders",
    "factor" => "factors",                             # main.cpp:149-169
    "algebraic_letters_allocate" => "idx",
)

# ---------------------------------------------------------------------------
# Version gate (RT-vi): stamp + behavioral canary, memoized per binary.
# ---------------------------------------------------------------------------

const _HF_GATE_OK = Dict{String,Bool}()

"""
    hf_gate(cfg::HFConfig=HFConfig()) -> Bool

Version gate = STAMP + BEHAVIORAL CANARY (DESIGN.md RT-vi). Stamp:
`hyperflint.sh --version` must print `HF_VERSION: v` with v ≥ 1.2.8
(live: `HF_VERSION: 1.2.8`). A stamp of 1.2.0/unknown means a
source build without `HF_VERSION` set at cmake — the error names that
footgun. Canary: `{"op":"eval","a":"-x^2","vars":["x"],"values":["3"]}`
must return `"-9"` (unary-minus precedence fix `-x^2 == -(x^2)`; pre-fix
builds silently return `"9"`). Either check failing throws `HFGateError`
— never a warning. Success is memoized per `cfg.bin`; failures are not
cached (a fixed environment may retry).
"""
function hf_gate(cfg::HFConfig=HFConfig())
    get(_HF_GATE_OK, cfg.bin, false) && return true
    isfile(cfg.bin) || throw(HFGateError("HyperFLINT wrapper not found at $(cfg.bin)"))
    # --- stamp ---
    r = _spawn_capture(_bash_ulimit_argv(cfg.bin, ["--version"]), nothing;
                       timeout_s=60, env=_hf_env(cfg))
    r.timed_out && throw(HFGateError("--version timed out"))
    m = match(r"HF_VERSION:\s*([0-9]+)\.([0-9]+)\.([0-9]+)", r.stdout * "\n" * r.stderr)
    m === nothing && throw(HFGateError("no HF_VERSION stamp in --version output " *
        "— pre-1.2.x binary or broken wrapper; stdout: " * _tail(r.stdout)))
    ver = VersionNumber(parse(Int, m.captures[1]), parse(Int, m.captures[2]),
                        parse(Int, m.captures[3]))
    ver >= v"1.2.8" || throw(HFGateError("HF_VERSION stamp $ver < 1.2.8. " *
        "KNOWN FOOTGUN (DESIGN.md RT-vi): source builds mis-stamp 1.2.0 unless " *
        "HF_VERSION is set at cmake — if this is a fresh lib/CLI build of new " *
        "code, fix the stamp at build time; a stamp-only pass/fail is never " *
        "trusted, but a FAILING stamp is fatal by design"))
    # --- behavioral canary ---
    resp = _hf_transport(Dict{String,Any}("op" => "eval", "a" => "-x^2",
                                          "vars" => ["x"], "values" => ["3"]),
                         cfg; timeout_s=60)
    res = get(resp, "result", nothing)
    if res == "9"
        throw(HFGateError("unary-minus canary FAILED: eval(-x^2 @ x=3) returned 9, " *
            "expected -9 — this binary predates the precedence fix that motivated " *
            "the v1.2.8 pin (version-stamp footgun: the stamp said $ver). REFUSING."))
    elseif res != "-9"
        throw(HFGateError("canary returned unexpected payload $(repr(res)) " *
            "(taxonomy keys: $(collect(keys(resp))))"))
    end
    _HF_GATE_OK[cfg.bin] = true
    return true
end

# ---------------------------------------------------------------------------
# hf_call — single invocation with full error taxonomy
# ---------------------------------------------------------------------------

"""
    hf_call(op::String, fields::Dict; cfg::HFConfig=HFConfig(),
            timeout_s::Real=cfg.timeout_s, gate::Bool=true,
            transport::Symbol=cfg.transport) -> Dict

Single HF invocation (CONTRACTS.md §a). Builds the flat request
(op + fields; `mzv_data_path` injected for MZV-consuming ops when absent),
sanity-checks flatness (the C++ JSON parser is regex-based), runs the
stamp+canary version gate once per binary (memoized), spawns
`bash -c 'ulimit -v 32505856; exec <bin> eval-json'` with the env contract
(OMP_NUM_THREADS, HF_MAX_THREADS_PER_CALL, HF_LR_TIME_BUDGET_S,
HF_NARROW_CTX when armed), parses stdout's LAST line, and classifies the
in-band taxonomy:

    {"error":..}                     -> HFError (fatal; incl. "unknown op")
    {"failed":true}                  -> HFFailed (cascade is the CALLER's job)
    {"divergent":true,..}            -> HFDivergent (FATAL standalone)
    {"narrow_ctx_insufficient":true} -> ONE retry with HF_NARROW_CTX=0,
                                        then HFNarrowCtx (no recursion)

Exit code is NEVER gated on (rc=0 even on handled failure — print-loud-rc0
footgun class); the known per-op result field must be present. Requests
combining HF_LAZY_SUM=1 with check_divergences are refused pre-flight
(engine refuses the combination — CONTRACTS.md §a).

`transport=:cabi` (OPT-IN; default `:cli` = the validation oracle) routes
ABI-covered ops (partial_fractions / linear_factors / find_lr_orders[_scan]
/ factor_table / hyperflint) through the in-process F2-patched C-ABI lib
(see the C-ABI section above for the LAWS: process-per-worker, ccalls
serialized, threads-on-one-handle forbidden, timeout_s not enforced
in-process); non-ABI ops fall back to the CLI (C-ABI gap G1). The
in-band error taxonomy is IDENTICAL on both transports (same upstream
handler bodies — CLI-vs-ABI byte-identity is test-gated in
test_cabi_transport.jl).
"""
function hf_call(op::String, fields::Dict; cfg::HFConfig=HFConfig(),
                 timeout_s::Real=cfg.timeout_s, gate::Bool=true,
                 transport::Symbol=cfg.transport)
    req = Dict{String,Any}("op" => op)
    for (k, v) in fields
        req[String(k)] = v
    end
    _check_required(op, req)     # programmer error -> ErrorException, pre-spawn
    gate && hf_gate(cfg)
    if op in _MZV_OPS && !haskey(req, "mzv_data_path")
        req["mzv_data_path"] = cfg.mzv_data_path
    end
    if get(ENV, "HF_LAZY_SUM", "0") == "1" && get(req, "check_divergences", false) === true
        throw(HFError(op, "HF_LAZY_SUM=1 combined with check_divergences is " *
            "REFUSED by the engine (CONTRACTS.md §a) — unset one before calling", nothing))
    end
    _assert_flat(op, req)
    resp = _hf_transport(req, cfg; timeout_s=timeout_s, transport=transport)
    # narrow-ctx: one retry with the flag explicitly off, never recursive.
    if get(resp, "narrow_ctx_insufficient", false) === true
        if cfg.narrow_ctx
            resp = _hf_transport(req, cfg; timeout_s=timeout_s, narrow_ctx=false,
                                 transport=transport)
            get(resp, "narrow_ctx_insufficient", false) === true &&
                throw(HFNarrowCtx(op, resp))
        else
            throw(HFNarrowCtx(op, resp))
        end
    end
    haskey(resp, "error") && throw(HFError(op, String(resp["error"]), resp))
    get(resp, "failed", false) === true &&
        throw(HFFailed(op, _tail(JSON.json(req)), resp))
    get(resp, "divergent", false) === true &&
        throw(HFDivergent(op, String(get(resp, "reason", "(no reason field)")), resp))
    rk = get(_RESULT_KEYS, op, nothing)
    if rk !== nothing && !haskey(resp, rk)
        throw(HFError(op, "response lacks the '$rk' payload field AND carries no " *
            "error key — contract drift, refusing to guess (keys: " *
            "$(collect(keys(resp))))", resp))
    end
    return resp
end

# ---------------------------------------------------------------------------
# hyperflint op wrapper + period resolution
# ---------------------------------------------------------------------------

function _parse_alg_letters(raw)::Vector{AlgLetter}
    out = AlgLetter[]
    raw isa AbstractVector || return out
    for e in raw
        e isa AbstractDict || continue
        push!(out, AlgLetter(Int(e["idx"]), String(get(e, "poly", "")),
                             Symbol(get(e, "var", "x")), String(get(e, "lc", "1")),
                             String(get(e, "sum", "0")), String(get(e, "product", "0")),
                             String(get(e, "disc", "0"))))
    end
    return out
end

"""
    hf_integrate(expr::String, vars_int::Vector{Symbol}, vars::Vector{Symbol};
                 cfg=HFConfig(), check_divergences::Bool=true,
                 algebraic_letters::Bool=false, carry_discharge::Bool=false,
                 timeout_s::Real=cfg.timeout_s)
        -> (terms::Vector{Tuple{HlogExpr,Vector{Vector{String}}}},
            letters::Vector{AlgLetter})

The `hyperflint` op (CONTRACTS.md §a; request per SubTropica.wl:12357):
`vars_int` is THE ordered integration order (each 0→∞, flat measure);
kinematic symbols = vars ∖ vars_int; `canonical_emission:true` always
(deterministic bytes). Every response coef string is parsed by
`parse_coef` against the response's `algebraic_letters` table; keys are
returned as raw letter-string words (each term = coef · ∏ ZeroInfPeriod;
empty key = pure constant — SubTropica.wl:12926-12934). Residual keys are
NOT evaluated here — see `resolve_periods`. The response's `"vars"` is the
AUGMENTED context (kin + MZV pool + Wm/Wp) and is never round-tripped.

check_divergences defaults hard-ON (standalone policy,
SubTropica.wl:12374-12384 + the silent-0-on-divergent control outcome,
CONTRACTS.md §a); `{"divergent":true}` raises `HFDivergent`. An EMPTY result
list with checks ON is a genuine zero (upstream: wl:12924).

Wm/Wp coef path (algebraic-letter gap closure, CONTRACTS.md §d):
a returned coef that `parse_coef` refuses AND that carries algebraic-letter
tokens is normalized through HF's OWN RT-v chain (`normalize_wmwp`,
pinned order simplify_with_vieta -> combine_wm_wp_ratios ->
back_substitute, the response's `algebraic_letters` table replayed) and
the NORMALIZED form is parsed; the refusal stays typed only if that form
still fails. Coefs without Wm/Wp tokens rethrow unchanged.
"""
function hf_integrate(expr::String, vars_int::Vector{Symbol},
                      vars::Vector{Symbol}; cfg::HFConfig=HFConfig(),
                      check_divergences::Bool=true,
                      algebraic_letters::Bool=false,
                      carry_discharge::Bool=false,
                      timeout_s::Real=cfg.timeout_s)
    issubset(vars_int, vars) ||
        throw(HFError("hyperflint", "vars_int ⊄ vars: $(setdiff(vars_int, vars))", nothing))
    fields = Dict{String,Any}(
        "expr" => expr,
        "vars_int" => [string(v) for v in vars_int],
        "vars" => [string(v) for v in vars],
        "check_divergences" => check_divergences,
        "algebraic_letters" => algebraic_letters,
        "carry_discharge" => carry_discharge,
        "canonical_emission" => true,
    )
    resp = hf_call("hyperflint", fields; cfg=cfg, timeout_s=timeout_s)
    letters = _parse_alg_letters(get(resp, "algebraic_letters", nothing))
    kin = Symbol[v for v in vars if !(v in vars_int)]
    terms = Tuple{HlogExpr,Vector{Vector{String}}}[]
    for t in resp["result"]
        cs = String(t["coef"])
        coefe = try
            parse_coef(cs, kin; letters=letters)
        catch e
            (e isa SubTropicaParseError && occursin(_ALGLETTER_TOKEN_RE, cs)) || rethrow()
            (ce, _) = normalize_wmwp(cs, letters; vars=kin, cfg=cfg)
            ce
        end
        key = Vector{Vector{String}}()
        for w in get(t, "key", Any[])
            push!(key, String[String(l) for l in w])
        end
        push!(terms, (coefe, key))
    end
    return (terms, letters)
end

"""
    zero_inf_period_value(word::Vector{String}; cfg=HFConfig(),
                          letters=AlgLetter[], vars::Vector{Symbol}=Symbol[],
                          strict::Bool=true) -> HlogExpr

Value of one residual ZeroInfPeriod word:

  * all-integer letters — HF's `zero_inf_period` op
    (main.cpp:2086-2099; live pin against the v1.2.8 binary: ["0","-1"]→mzv_2,
    ["0","-1","-1"]→mzv_3, ["-1"]→0, ["-2"]→-Log2; NOTE this supersedes
    CONTRACTS.md §a's "1/2*mzv_2" example, which does not match the live
    v1.2.8 binary — PROVENANCE.md: live probe wins). The op refuses
    non-integer letters in-band ("Phase 6c scope") → HFError.
  * single NON-integer letter ℓ — elementary closed form
    ZeroInfPeriod[{ℓ}] = -Log[-ℓ]: G(ℓ;z) = log((z-ℓ)/(-ℓ)) → constant
    term -log(-ℓ) (log z → 0), consistent with the live integer pins
    ["-2"]→-Log2, ["-1"]→0 and with the kinematic probe
    (coef -1/(a-1), key ["-a"]) ⇒ log(a)/(a-1). `vars` provides the
    kinematic context for parsing ℓ.
  * multi-letter non-integer word — REFUSED (strict=true, default) with a
    typed error: those are genuine polylogarithms of kinematics and belong
    to the Hlog/fibration pipeline, not the boundary-period path.
    strict=false keeps the word as an HPeriod atom instead (never silent).
"""
function zero_inf_period_value(word::Vector{String}; cfg::HFConfig=HFConfig(),
                               letters::Vector{AlgLetter}=AlgLetter[],
                               vars::Vector{Symbol}=Symbol[], strict::Bool=true)
    isempty(word) && throw(SubTropicaEvalError("empty ZeroInfPeriod word (empty key means constant term — never stored, types.jl)"))
    if any(w -> occursin(_ALGLETTER_TOKEN_RE, w), word)
        # Algebraic-letter boundary words (letters-tier keys carry Wm_i/Wp_i
        # — live, eq 4.11 k0/R'1): function-valued over the
        # letter roots; none of the rational closed-form branches below
        # apply (their _parse_rat_frac calls would die on the atom — a
        # WRONG error type, not a refusal). Typed refusal here, or an
        # explicit HPeriod atom under strict=false (never silent).
        strict && throw(SubTropicaEvalError("ZeroInfPeriod word $(word) carries " *
            "algebraic-letter tokens — boundary period over the Wm/Wp roots " *
            "with no exact atom form; pass strict=false to keep an explicit " *
            "HPeriod atom (eval_symbolic evaluates it numerically against " *
            "the LaurentSeries.letters table — zero_inf_period_alg)"))
        return HlogExpr([HTerm("1", Pair{HAtom,Int}[HPeriod(copy(word)) => 1])],
                        copy(vars))
    end
    all(w -> strip(w) == "0", word) && return HlogExpr(HTerm[], copy(vars))
        # all-zero word: (log z)^k/k! → 0 under the log z → 0 regularization
    if length(word) == 2 && count(w -> strip(w) == "0", word) == 1
        # ONE zero + one letter a = -c (rational or kinematic), weight 2:
        #   Z({0,a}) = reg const of G(0,a;z) = -Li2(z/a)
        #            = π²/6 + (1/2)·log²(c)      (Li2 inversion, z/a → -∞;
        #              consistent with the live pins {0,-1} → mzv_2 and the
        #              engine's {0,-2})
        #   Z({a,0}) = -Z({0,a})                 (shuffle: Z({a})·Z({0}) = 0)
        # Needed live (eq 4.3 face [2,4,10] check-OFF fallback →
        # residual {0,-1/2}); ginac-validated (see gate log). Non-positive
        # numeric c refuses loudly below.
        lz = strip(word[1]) == "0"                 # true: {0,a}; false: {a,0}
        a  = lz ? word[2] : word[1]
        ctx = _pctx(vars, letters)
        fr = _parse_rat_frac("-(" * a * ")", ctx)  # c = -a
        nu, de = numerator(fr), denominator(fr)
        if is_constant(nu) && is_constant(de)
            q = Rational{BigInt}(BigInt(numerator(leading_coefficient(nu))),
                                 BigInt(denominator(leading_coefficient(nu)))) //
                Rational{BigInt}(BigInt(numerator(leading_coefficient(de))),
                                 BigInt(denominator(leading_coefficient(de))))
            q > 0 || throw(SubTropicaEvalError("zero+letter ZeroInfPeriod " *
                "$(word): -a = $q ≤ 0 — non-positive rescale point; REFUSING"))
        end
        s = lz ? "1" : "-1"
        s2 = lz ? "1/2" : "-1/2"
        terms = HTerm[HTerm(s, Pair{HAtom,Int}[HZeta([2]) => 1])]
        isone(fr) ||
            push!(terms, HTerm(s2,
                Pair{HAtom,Int}[_norm_atom(HLogA(_frac_string(fr))) => 2]))
        return HlogExpr(terms, copy(vars))
    end
    nz = unique(w for w in word if strip(w) != "0")
    if length(word) == 3 && length(nz) == 1 && count(w -> strip(w) == "0", word) >= 1
        # Weight-3 words over {0, a} with ONE distinct nonzero letter a = -c.
        # Derivations (u = 1/(1+s) substitution + shuffle algebra; all
        # ginac-validated):
        #   Z({0,a,a}) =  ζ₃ − (1/6)log³c            (c=1 ⇒ mzv_3: live pin ✓)
        #   Z({0,0,a}) = −(1/6)log³c − ζ₂·log c      (c=1 ⇒ 0: live pin ✓)
        #   Z({a,0,a}) = −2ζ₃ − (1/6)log³c − ζ₂·log c        (shuffle a⧢{0,a})
        #   Z({a,a,0}) =  ζ₃ + (1/3)log³c + ζ₂·log c         (shuffle 0⧢{a,a})
        #   Z({0,a,0}) =  (1/3)log³c + 2ζ₂·log c             (shuffle 0⧢{0,a})
        #   Z({a,0,0}) = −(1/6)log³c − ζ₂·log c              (shuffle 0⧢{a,0})
        ctx = _pctx(vars, letters)
        fr = _parse_rat_frac("-(" * nz[1] * ")", ctx)   # c = -a
        nu, de = numerator(fr), denominator(fr)
        if is_constant(nu) && is_constant(de)
            q = Rational{BigInt}(BigInt(numerator(leading_coefficient(nu))),
                                 BigInt(denominator(leading_coefficient(nu)))) //
                Rational{BigInt}(BigInt(numerator(leading_coefficient(de))),
                                 BigInt(denominator(leading_coefficient(de))))
            q > 0 || throw(SubTropicaEvalError("weight-3 ZeroInfPeriod $(word): " *
                "-a = $q ≤ 0 — non-positive rescale point; REFUSING"))
        end
        pat = Tuple(strip(w) == "0" ? 0 : 1 for w in word)
        # coefficients (z3, l3, z2l) of ζ₃ + log³c + ζ₂·log c per pattern:
        (z3, l3, z2l) =
            pat == (0, 1, 1) ? (1//1, -1//6, 0//1)  :
            pat == (0, 0, 1) ? (0//1, -1//6, -1//1) :
            pat == (1, 0, 1) ? (-2//1, -1//6, -1//1) :
            pat == (1, 1, 0) ? (1//1, 1//3, 1//1)   :
            pat == (0, 1, 0) ? (0//1, 1//3, 2//1)   :
            pat == (1, 0, 0) ? (0//1, -1//6, -1//1) :
            throw(SubTropicaEvalError("weight-3 ZeroInfPeriod $(word): unhandled pattern $(pat)"))
        _rs(r::Rational) = denominator(r) == 1 ? string(numerator(r)) :
                           string(numerator(r)) * "/" * string(denominator(r))
        terms = HTerm[]
        iszero(z3) || push!(terms, HTerm(_rs(z3),
            _canon_atomvec(Pair{HAtom,Int}[HZeta([3]) => 1])))
        if !isone(fr)
            la = HLogA(_frac_string(fr))   # _canon_atomvec applies _norm_atom
            iszero(l3)  || push!(terms, HTerm(_rs(l3),
                _canon_atomvec(Pair{HAtom,Int}[la => 3])))
            iszero(z2l) || push!(terms, HTerm(_rs(z2l),
                _canon_atomvec(Pair{HAtom,Int}[HZeta([2]) => 1, la => 1])))
        end
        return HlogExpr(terms, copy(vars))
    end
    if all(w -> w == word[1], word)
        # Elementary closed form for EQUAL-letter words (integer or kinematic),
        # any multiplicity k: ZeroInfPeriod[{ℓ,…,ℓ}(k)] = (1/k!)·(G(ℓ;z)|_reg)^k
        # = (-log(-ℓ))^k / k!  — the shuffle power of the weight-1 period
        # (G(ℓ;z) = log((z-ℓ)/(-ℓ)); constant term -log(-ℓ) under the
        # log z → 0 regularization). Consistent with the live integer pins
        # ["-1"]→0, ["-2"]→-Log2. The HF op is NOT used here: the v1.2.8
        # engine only rescales single letters -1/-2 ("single-letter rescale
        # requires scale=2") — hit live on the eq 4.3 deep face
        # (letter -4, then word {-1/4,-1/4}; upstream v1.2.8 release notes). Non-positive
        # numeric -ℓ refuses LOUDLY (divergent/complex boundary period —
        # not this branch's business).
        k = length(word)
        ctx = _pctx(vars, letters)
        fr = _parse_rat_frac("-(" * word[1] * ")", ctx)   # -ℓ
        nu, de = numerator(fr), denominator(fr)
        if is_constant(nu) && is_constant(de)
            q = Rational{BigInt}(BigInt(numerator(leading_coefficient(nu))),
                                 BigInt(denominator(leading_coefficient(nu)))) //
                Rational{BigInt}(BigInt(numerator(leading_coefficient(de))),
                                 BigInt(denominator(leading_coefficient(de))))
            q > 0 || throw(SubTropicaEvalError("equal-letter ZeroInfPeriod " *
                "$(word): -ℓ = $q ≤ 0 — non-positive rescale point; REFUSING"))
        end
        isone(fr) && return HlogExpr(HTerm[], copy(vars)) # log(1) = 0
        c = Rational{BigInt}(big(-1)^k) // factorial(big(k))
        t = HTerm(denominator(c) == 1 ? string(numerator(c)) :
                  string(numerator(c)) * "/" * string(denominator(c)),
            Pair{HAtom,Int}[_norm_atom(HLogA(_frac_string(fr))) => k])
        return HlogExpr([t], copy(vars))
    elseif all(w -> occursin(r"^-?\d+$", strip(w)), word)
        resp = hf_call("zero_inf_period",
                       Dict{String,Any}("word" => word, "vars" => ["x"]); cfg=cfg)
        return parse_coef(String(resp["result"]), Symbol[]; letters=letters)
    elseif strict
        throw(SubTropicaEvalError("multi-letter non-integer ZeroInfPeriod word " *
            "$(word): HF's zero_inf_period is integer-letter only (live probe: " *
            "'non-integer letter (Phase 6c scope)') and kinematic boundary " *
            "words are function-valued — route through the Hlog pipeline; " *
            "pass strict=false to keep an explicit HPeriod atom"))
    else
        return HlogExpr([HTerm("1", Pair{HAtom,Int}[HPeriod(copy(word)) => 1])],
                        copy(vars))
    end
end

"""
    resolve_periods(terms; cfg=HFConfig(), letters=AlgLetter[],
                    strict::Bool=true, on_period=nothing) -> HlogExpr

Assemble `hf_integrate` output terms into a single HlogExpr: each term =
coef · ∏ⱼ ZeroInfPeriod[wordⱼ] (boundary periods commute; shuffle product
= ordinary product — SubTropica.wl:12926-12993). Word values come from
`zero_inf_period_value` (HF op / single-letter closed form / strict
refusal). Per-word results are cached within the call.

`on_period(word, value_expr)` is the verification hook (DESIGN.md
RT-extra-period / gate A(v)): verify.jl and the unit tests use it to
cross-evaluate EVERY resolved residual key against the independent
ginac_gpl evaluator (`zero_inf_period_ginac`) to ≥60 dps.

Wm/Wp-bearing coefs should be normalized through HF's own ops FIRST
(`normalize_wmwp`, DESIGN.md RT-v).
"""
function resolve_periods(terms; cfg::HFConfig=HFConfig(),
                         letters::Vector{AlgLetter}=AlgLetter[],
                         strict::Bool=true, on_period=nothing)
    acc = HlogExpr()
    cache = Dict{Vector{String},HlogExpr}()
    for (coefe, key) in terms
        val = coefe
        for word in key
            pv = get!(cache, word) do
                zero_inf_period_value(word; cfg=cfg, letters=letters,
                                      vars=coefe.vars, strict=strict)
            end
            on_period === nothing || on_period(word, pv)
            val = hlog_mul(val, pv)
        end
        acc = hlog_add(acc, val)
    end
    return acc
end
