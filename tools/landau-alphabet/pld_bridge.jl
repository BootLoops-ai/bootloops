#!/usr/bin/env julia
#=
pld_bridge.jl — thin bridge from the Landau Alphabet graph-spec JSON to PLD.jl
(Fevola–Mizera–Telen, arXiv:2311.14669) so the bounded-sympy engine and PLD.jl
can be CROSS-VALIDATED on the SAME input.

Usage:
    julia +1.10 pld_bridge.jl SPEC.json OUT.json [--method sym|num]
  or
    PLD_ENV=<your PLD julia env> julia pld_bridge.jl ...

Reads the same schema as landau_alphabet.py (edges = [[v1,v2,m2],...],
nodes = {vertex: [legs...]}).  Translates to PLD.jl's
    getPLD(edges, nodes; internal_masses=..., external_masses=..., method=:sym)
and writes the discriminant list under the same {"alphabet": [...], "faces": [...],
"honesty": {...}} schema (faces = PLD's per-weight records).

NOTE: PLD.jl assigns its OWN Mandelstam names (s,t for n=4; s_ij for n>=5); the
spec's "channels" map is therefore informational only — you may need a trivial
relabel when comparing to the sympy backend's output.

VERSION NOTE: PLD.jl + Oscar/GAP do NOT load on Julia 1.12
(GAP unsupported). This bridge auto-activates the working Julia-1.10 env
in the env named by ENV["PLD_ENV"] (known-good: Oscar 1.7.3 + HC 2.20). PLD.jl
may need an Oscar-1.x compat shim. If Oscar
fails to load, the bridge falls back to a HomotopyContinuation-only
discriminant (numerical) on the Lee-Pomeransky G=U+F polynomial:
per face (connected edge subset), the torus critical system
{G_face = 0, dG_face/dx_e = 0} is solved on a seeded generic rational LINE
in invariant space, and the surviving t-values are discriminant WITNESSES.
With --check CANDIDATES.json (an output of the sympy backend, or any JSON
with an "alphabet" list) each witness is matched against the candidate
letters: an UNMATCHED witness means no candidate covers that discriminant
surface (exit 2 — the cross-validation catch).  Witnesses are INTERIOR
pinches: second-type surfaces a DE-level alphabet prunes ARE witnessed
(box1l: s+t on the leading face), while endpoint-class coefficient letters
(box1l: s, t) produce no witness and stay the sympy chain's territory.
The fallback writes the same output schema with complete=false always: it
produces numerical witnesses, never a completeness certificate.  Flags:
--check FILE, --seed N, --tol X, --fallback (force the HC fallback even where
PLD loads).

FAST GROEBNER BACKEND: for the HEAVY faces that PLD.jl's numerical
homotopy errors on and the sympy resultant chain times out on (e.g. a 3-loop
K4 leading face, deg-267, and its 5-edge "3m+2gamma" subfaces), use the
dedicated msolve backend instead of this Julia bridge:
    python tools/landau-alphabet/pld_groebner_backend.py FACE_UF.json \
           --field s,t,m2 --backend msolve
It computes the per-face discriminant by Rabinowitsch-saturated elimination via
`msolve` (K4 leading: 32 s). It is also wired as an automatic escalation inside
landau_alphabet.py's discriminant_letters() (PLD_GROEBNER=auto|force|off).
See the Groebner-backend member section in this package's README.md.
=#

# --- environment auto-detect / activation -----------------------------------
import Pkg
const _PLD_ENV = get(ENV, "PLD_ENV", "")  # a julia env with PLD.jl deps (obtain PLD.jl upstream)
if VERSION < v"1.11" && isdir(_PLD_ENV)
    # Known-good env (Oscar 1.7.3 + GAP works on 1.10)
    Pkg.activate(_PLD_ENV; io=devnull)
elseif VERSION >= v"1.11"
    @warn ("PLD.jl/GAP is known-broken on Julia ≥1.11. Re-invoke this bridge "
           * "with `julia +1.10 ...` (juliaup channel 1.10 is installed). "
           * "Attempting anyway; will fall back to HC-only if Oscar fails.")
    isdir(_PLD_ENV) && Pkg.activate(_PLD_ENV; io=devnull)
end

haskey(ENV, "PLD_SRC") && push!(LOAD_PATH, ENV["PLD_SRC"])  # .../PLD/src checkout

_have_pld = true
try
    @eval using PLD
    @eval using Oscar
catch e
    @warn "PLD/Oscar failed to load ($(e)); falling back to HomotopyContinuation-only."
    global _have_pld = false
end
using JSON
using HomotopyContinuation

function load_spec(path)
    spec = JSON.parsefile(path)
    # --- edges -> Vector{Vector{Int}}, internal masses -------------------
    raw_edges = spec["edges"]
    edges = [[Int(e[1]), Int(e[2])] for e in raw_edges]
    # internal_masses: PLD wants :zero, :equal, :generic, or a vector of HC vars.
    # (PLD's getUF substitutes these into an HC expression, NOT an Oscar one.)
    msyms = unique([string(e[3]) for e in raw_edges if length(e) >= 3 && string(e[3]) != "0"])
    if isempty(msyms)
        internal_masses = :zero
    else
        mvars = Dict(s => HomotopyContinuation.Variable(Symbol(s)) for s in msyms)
        internal_masses = Any[]
        for e in raw_edges
            m = (length(e) >= 3) ? string(e[3]) : "0"
            push!(internal_masses, m == "0" ? 0 : mvars[m])
        end
    end
    # --- nodes -> Vector{Int} (FLAT list of vertices with external legs, in
    #     leg order). PLD's getUF expects nodes[i] = the VERTEX index that
    #     external leg i is attached to (a length-n flat Int vector), NOT the
    #     spec's {vertex: [legs...]} dict.  Invert the spec map.
    nd = spec["nodes"]
    leg2vtx = Dict{Int,Int}()
    for (v, legs) in nd
        for l in legs
            leg2vtx[Int(l)] = parse(Int, string(v))
        end
    end
    nodes = [leg2vtx[l] for l in sort(collect(keys(leg2vtx)))]
    # --- external masses -------------------------------------------------
    kin = get(spec, "kinematics", Dict())
    extm = get(kin, "external_masses", Dict())
    nlegs = length(nodes)
    if all(string(get(extm, string(i), "0")) == "0" for i in 1:nlegs)
        external_masses = :zero
    else
        # symbolic external masses M_i for the off-shell legs
        external_masses = Any[
            (string(get(extm, string(i), "0")) == "0" ? 0 :
             HomotopyContinuation.Variable(Symbol(string(get(extm, string(i), "M$i")))))
            for i in 1:nlegs ]
    end
    return spec, edges, nodes, internal_masses, external_masses
end

# ============================================================================
# HomotopyContinuation-only numerical fallback (the route promised in the
# VERSION NOTE above).  Everything below needs ONLY HomotopyContinuation +
# JSON + stdlib — no PLD.jl, no Oscar.
#
# Method: same graph-face lattice as the sympy backend (connected edge
# subsets, >=2 vertices), same Symanzik build driven by the spec's own
# channels map.  Per face, the Lee-Pomeransky polynomial G = U + F is
# restricted to a seeded generic rational line z(t) = z0 + t*z1 in invariant
# space and the SQUARE torus critical system
#     { G(x; z(t)) = 0,  dG/dx_e = 0  for every face parameter x_e }
# is solved by homotopy continuation in (x, t).  The t-values of the finite
# nonsingular torus solutions are numerical WITNESSES of the face's
# G-discriminant on the line.  With --check they are matched against a
# candidate letter set; an unmatched witness is a caught missing letter.
# ============================================================================
using Random

const _HCV = HomotopyContinuation.Variable

# HC's symbolic Expression/Variable types subtype Number, so `isa Number` does
# NOT mean "a plain constant" — use this everywhere a literal number must be
# told apart from a symbolic value.
_is_plain_num(x) = x isa Union{Integer,AbstractFloat,Rational,Complex}

_kin_ev(x::Number, vars) = x
function _kin_ev(x::Symbol, vars)
    haskey(vars, string(x)) || error("unknown symbol '$x' in a kinematics " *
        "expression (declare it in kinematics.invariants)")
    return vars[string(x)]
end
function _kin_ev(x::Expr, vars)
    x.head === :call || error("unsupported kinematics expression: $x")
    op = x.args[1]
    as = [_kin_ev(a, vars) for a in x.args[2:end]]
    op === :+ && return sum(as)
    op === :* && return prod(as)
    op === :- && return length(as) == 1 ? -as[1] : as[1] - as[2]
    op === :^ && return as[1]^as[2]
    op === :/ && return as[1] / as[2]
    error("unsupported operator '$op' in a kinematics expression")
end

"Parse a kinematics/letter string (sympy syntax; '**' for powers) into an
HC Expression over the invariant Variables in `vars`."
function _parse_kin_expr(s, vars::Dict{String,Any})
    return _kin_ev(Meta.parse(replace(string(s), "**" => "^")), vars)
end

"All k-subsets of v, lexicographic (k = 0 yields the one empty subset)."
function _ksubsets(v::Vector{Int}, k::Int)
    n = length(v)
    out = Vector{Vector{Int}}()
    k == 0 && return [Int[]]
    k > n && return out
    idx = collect(1:k)
    while true
        push!(out, v[idx])
        i = k
        while i >= 1 && idx[i] == n - k + i
            i -= 1
        end
        i == 0 && break
        idx[i] += 1
        for j in i+1:k
            idx[j] = idx[j-1] + 1
        end
    end
    return out
end

function _uf_root(par::Dict{Int,Int}, a::Int)
    while par[a] != a
        par[a] = par[par[a]]
        a = par[a]
    end
    return a
end

"Symanzik (U, F) of the face spanned by edge indices `Es` (others deleted).
Mirrors the sympy backend's spanning-tree / 2-forest construction, with the
spec's own channels map supplying every (sum of legs)^2 invariant."
function _face_symanzik(Es::Vector{Int}, edges, m2vals, chan, xv, nodes)
    vs = sort(unique(vcat([[edges[e][1], edges[e][2]] for e in Es]...)))
    length(vs) >= 2 || return nothing
    nV = length(vs)
    Lf = length(Es) - nV + 1
    Lf < 0 && return nothing
    # U: sum over spanning trees of the product of the DELETED edge params
    U = Expression(0)
    for T in _ksubsets(Es, nV - 1)
        par = Dict(v => v for v in vs)
        ok = true
        for e in T
            a, b = edges[e]
            ra, rb = _uf_root(par, a), _uf_root(par, b)
            if ra == rb
                ok = false
                break
            end
            par[ra] = rb
        end
        (ok && length(unique(_uf_root(par, v) for v in vs)) == 1) || continue
        term = Expression(1)
        for e in Es
            e in T || (term *= xv[e])
        end
        U += term
    end
    # F0: sum over spanning 2-forests of -(channel of one component)
    F0 = Expression(0)
    for fo in _ksubsets(Es, nV - 2)
        par = Dict(v => v for v in vs)
        ok = true
        for e in fo
            a, b = edges[e]
            ra, rb = _uf_root(par, a), _uf_root(par, b)
            if ra == rb
                ok = false
                break
            end
            par[ra] = rb
        end
        ok || continue
        comps = Dict{Int,Vector{Int}}()
        for v in vs
            push!(get!(comps, _uf_root(par, v), Int[]), v)
        end
        length(comps) == 2 || continue
        # channel of the cut: take either component's attached-leg set; if the
        # two disagree (a sub-face where a legged vertex is deleted), prefer
        # the NONZERO one — deterministic, and over-inclusive is safe for a
        # witness check
        c = 0
        for comp in values(comps)
            legs = Set{Int}()
            for v in comp, l in get(nodes, v, Int[])
                push!(legs, l)
            end
            cc = get(chan, legs, 0)     # unknown channel -> 0 (scaleless cut)
            if !(_is_plain_num(cc) && cc == 0)
                c = cc
                break
            end
        end
        (_is_plain_num(c) && c == 0) && continue
        term = Expression(1)
        for e in Es
            e in fo || (term *= xv[e])
        end
        F0 += term * (-c)
    end
    massterm = Expression(0)
    for e in Es
        m = m2vals[e]
        (_is_plain_num(m) && m == 0) && continue
        massterm += xv[e] * m
    end
    F = F0 + U * massterm
    return (U = U, F = F, L = Lf, G = U + F)
end

function run_hc_fallback(specpath::String, outpath::String;
                         check::String = "", seed::Int = 0, tol::Float64 = 1e-6)
    spec = JSON.parsefile(specpath)
    name = get(spec, "name", "graph")
    raw_edges = spec["edges"]
    E = length(raw_edges)
    edges = [(Int(e[1]), Int(e[2])) for e in raw_edges]
    kin = get(spec, "kinematics", Dict())
    inv_names = [string(s) for s in get(kin, "invariants", [])]
    isempty(inv_names) && error("the HC fallback needs kinematics.invariants in the spec")
    any(occursin(r"^x\d+$", n) for n in inv_names) &&
        error("invariant names of the form x<N> collide with the Feynman parameters")
    iv = [_HCV(Symbol(n)) for n in inv_names]
    kd = Dict{String,Any}(n => iv[i] for (i, n) in enumerate(inv_names))
    xv = [_HCV(Symbol("x", i)) for i in 1:E]
    tv = _HCV(Symbol("_tline_"))
    nI = length(iv)
    # --- edge masses and channels (spec-driven, mirroring the sympy backend) --
    m2vals = Any[length(e) >= 3 ? _parse_kin_expr(e[3], kd) : 0 for e in raw_edges]
    nodes = Dict{Int,Vector{Int}}()
    for (v, legs) in get(spec, "nodes", Dict())
        nodes[parse(Int, string(v))] = [Int(l) for l in legs]
    end
    all_legs = isempty(nodes) ? Int[] : sort(unique(vcat(values(nodes)...)))
    full = Set(all_legs)
    chan = Dict{Set{Int},Any}()
    for (k, v) in get(kin, "channels", Dict())
        legs = Set(parse(Int, strip(x)) for x in split(string(k), ","))
        ex = _parse_kin_expr(v, kd)
        chan[legs] = ex
        haskey(chan, setdiff(full, legs)) || (chan[setdiff(full, legs)] = ex)
    end
    extm = Dict{Int,Any}()
    for (l, m) in get(kin, "external_masses", Dict())
        extm[parse(Int, string(l))] = _parse_kin_expr(m, kd)
    end
    for l in all_legs
        haskey(chan, Set([l])) || (chan[Set([l])] = get(extm, l, 0))
        haskey(chan, setdiff(full, Set([l]))) || (chan[setdiff(full, Set([l]))] = get(extm, l, 0))
    end
    haskey(chan, Set{Int}()) || (chan[Set{Int}()] = 0)
    haskey(chan, full) || (chan[full] = 0)
    # --- seeded generic rational line in invariant space ----------------------
    rng = MersenneTwister(seed)
    z0 = [rand(rng, -9:9) for _ in 1:nI]
    z1 = [(rand(rng, Bool) ? 1 : -1) * rand(rng, 1:9) for _ in 1:nI]
    zline = [Expression(z0[i]) + tv * z1[i] for i in 1:nI]
    # --- candidate letters (--check): compile per-letter evaluators on the line
    cand_strs = String[]
    if !isempty(check)
        cj = JSON.parsefile(check)
        seen = Set{String}()
        addc = function (x)
            s = string(x)
            (s == "" || s in seen) && return
            push!(seen, s)
            push!(cand_strs, s)
        end
        for key in ("alphabet", "spurious", "first_entry", "last_entry")
            haskey(cj, key) && foreach(addc, cj[key])
        end
        if haskey(cj, "faces")
            for f in cj["faces"], key in ("letters_F", "letters_U")
                haskey(f, key) && foreach(addc, f[key])
            end
        end
        isempty(cand_strs) && error("--check file has no candidate letters (need an 'alphabet' list)")
    end
    probes = ComplexF64[0.31 + 0.77im, -1.23 + 0.41im, 2.17 - 0.89im]
    cands = Tuple{String,System,Float64}[]
    for s in cand_strs
        local ex
        try
            ex = _parse_kin_expr(s, kd)
        catch
            @warn "candidate letter '$s' is not a polynomial in the invariants; skipped"
            continue
        end
        _is_plain_num(ex) && continue
        exline = subs(ex, iv => zline)
        lsys = System([exline]; variables = [tv])
        scale = maximum(abs(lsys([p])[1]) for p in probes) + 1e-30
        push!(cands, (s, lsys, scale))
    end
    # --- face walk ------------------------------------------------------------
    println("[$name] HC-only numerical fallback: line z(t)=z0+t*z1, seed=$seed, " *
            "z0=$z0, z1=$z1 over $(inv_names)")
    flush(stdout)
    t0 = time()
    faces_out = Any[]
    alpha = String[]
    unmatched = Any[]
    n_res = 0
    n_fail = 0
    for k in E:-1:2
        for fidx in _ksubsets(collect(1:E), k)
            sg = _face_symanzik(fidx, edges, m2vals, chan, xv, nodes)
            sg === nothing && continue
            Gline = subs(sg.G, iv => zline)
            gvars = Set(variables(Gline))
            rec = Dict{String,Any}("face_edges" => ["x$(i)" for i in fidx],
                                   "n_edges" => k, "loops_L" => sg.L)
            if !(tv in gvars)
                rec["status"] = "no-kinematic-part (scaleless on the line)"
                rec["face_status"] = "scaleless"
                rec["t_witnesses"] = Any[]
                push!(faces_out, rec)
                continue
            end
            fx = [xv[i] for i in fidx if xv[i] in gvars]
            eqs = Expression[Gline]
            append!(eqs, [differentiate(Gline, x) for x in fx])
            svars = vcat(fx, [tv])
            sols = nothing
            hcseed = UInt32(0x5eed) + UInt32(mod(seed, 10000))
            try
                try
                    res = solve(System(eqs; variables = svars);
                                show_progress = false, seed = hcseed)
                    sols = solutions(res; only_nonsingular = true)
                catch err
                    # a face system with constant-in-x equations can have no
                    # fine mixed cells (polyhedral start impossible); retry
                    # with the total-degree start, which handles it
                    err isa OverflowError || rethrow()
                    res = solve(System(eqs; variables = svars);
                                show_progress = false, seed = hcseed,
                                start_system = :total_degree)
                    sols = solutions(res; only_nonsingular = true)
                end
            catch err
                rec["status"] = "hc-solve-failed: $(sprint(showerror, err))"
                rec["face_status"] = "over_budget"
                rec["t_witnesses"] = Any[]
                push!(faces_out, rec)
                n_fail += 1
                flush(stdout)
                continue
            end
            ts = ComplexF64[]
            for s in sols
                xs = s[1:end-1]
                tval = s[end]
                (all(abs.(xs) .> 1e-8) && all(abs.(xs) .< 1e10) && abs(tval) < 1e10) || continue
                any(abs(tval - u) < 1e-6 * (1 + abs(u)) for u in ts) && continue
                push!(ts, tval)
            end
            rec["t_witnesses"] = [[real(u), imag(u)] for u in ts]
            if !isempty(cands)
                matched_here = String[]
                for u in ts
                    best = ""
                    bestv = Inf
                    for (cs, lsys, scale) in cands
                        v = abs(lsys([u])[1]) / scale
                        if v < bestv
                            best = cs
                            bestv = v
                        end
                    end
                    if bestv < tol
                        push!(matched_here, best)
                        best in alpha || push!(alpha, best)
                    else
                        push!(unmatched, Dict("face_edges" => rec["face_edges"],
                                              "t" => [real(u), imag(u)],
                                              "closest_candidate" => best,
                                              "residual" => bestv))
                    end
                end
                rec["matched_letters"] = sort(unique(matched_here))
            end
            rec["status"] = "resolved-numerical"
            rec["face_status"] = "done"
            n_res += 1
            println("[face $(rec["face_edges"])] $(length(ts)) torus witness(es)" *
                    (isempty(cands) ? "" : ", matched: $(get(rec, "matched_letters", []))"))
            flush(stdout)
            push!(faces_out, rec)
        end
    end
    # --- honest output (same schema; complete stays FALSE by design) ----------
    out = Dict(
        "name" => name,
        "engine" => ("pld_bridge.jl -> HomotopyContinuation-only numerical fallback " *
                     "(Lee-Pomeransky G=U+F, per-face torus critical systems on a " *
                     "generic kinematic line)"),
        "method" => "hc-fallback-numerical",
        "graph" => spec,
        "complete" => false,
        "alphabet" => sort(alpha),
        "first_entry" => String[],
        "last_entry" => String[],
        "spurious" => String[],
        "faces" => faces_out,
        "honesty" => Dict(
            "complete" => false,
            "n_faces" => length(faces_out),
            "resolved" => n_res,
            "bounded_partial" => n_fail,
            "bounded_timeout" => 0,
            "unresolved_faces" => Any[],
            "global_deadline_hit" => false,
            "line" => Dict("seed" => seed, "invariants" => inv_names,
                           "z0" => z0, "z1" => z1, "tol" => tol),
            "unmatched_witnesses" => unmatched,
            "secs" => round(time() - t0; digits = 2),
            "note" => ("HC-only NUMERICAL fallback: per-face G-discriminant witnesses " *
                       "on ONE generic line. complete stays false by design -- " *
                       "witnesses confirm candidate letters (--check) but never " *
                       "certify completeness."),
        ),
    )
    open(outpath, "w") do io
        JSON.print(io, out, 2)
    end
    nw = sum(length(r["t_witnesses"]) for r in faces_out; init = 0)
    if n_fail > 0
        println("written -> $outpath  ($n_res faces resolved, $n_fail FAILED)")
        exit(2)
    elseif !isempty(cands) && !isempty(unmatched)
        println("written -> $outpath  CROSS-VALIDATION CATCH: $(length(unmatched)) " *
                "witness(es) match NO candidate letter (tol=$tol) -- the candidate " *
                "alphabet is missing letters")
        exit(2)
    elseif !isempty(cands)
        println("written -> $outpath  ($nw witnesses, all matched; letters witnessed: $(sort(alpha)))")
    else
        println("written -> $outpath  ($nw numerical witnesses; pass --check CANDIDATES.json to match letters)")
    end
    return nothing
end

function main()
    if length(ARGS) < 2
        println(stderr, "usage: julia pld_bridge.jl SPEC.json OUT.json [--method sym|num] " *
                        "[--check CANDIDATES.json] [--seed N] [--tol X] [--fallback]")
        exit(1)
    end
    specpath, outpath = ARGS[1], ARGS[2]
    method = :sym
    # Resume support: pass through getPLD's native resume kwargs.
    # --load_output OLDLOG  scrapes "codim: c, face: i/N, ..., discriminant: ..."
    # lines from a prior run's log; --codim_start/--face_start restart the face
    # walk there; --save_output NEWLOG appends per-face results for the NEXT resume.
    codim_start = -1
    face_start  = 1
    load_output = ""
    save_output = ""
    check = ""
    seed = 0
    tol = 1e-6
    force_fallback = false
    for i in 3:length(ARGS)
        if ARGS[i] == "--method" && i < length(ARGS)
            method = Symbol(ARGS[i+1])
        elseif ARGS[i] == "--codim_start" && i < length(ARGS)
            codim_start = parse(Int, ARGS[i+1])
        elseif ARGS[i] == "--face_start" && i < length(ARGS)
            face_start = parse(Int, ARGS[i+1])
        elseif ARGS[i] == "--load_output" && i < length(ARGS)
            load_output = ARGS[i+1]
        elseif ARGS[i] == "--save_output" && i < length(ARGS)
            save_output = ARGS[i+1]
        elseif ARGS[i] == "--check" && i < length(ARGS)
            check = ARGS[i+1]
        elseif ARGS[i] == "--seed" && i < length(ARGS)
            seed = parse(Int, ARGS[i+1])
        elseif ARGS[i] == "--tol" && i < length(ARGS)
            tol = parse(Float64, ARGS[i+1])
        elseif ARGS[i] == "--fallback"
            force_fallback = true
        end
    end

    if !_have_pld || force_fallback
        _have_pld || @warn ("PLD/Oscar unavailable on this Julia ($(VERSION)); running the "
               * "HomotopyContinuation-only numerical fallback on the Lee-Pomeransky "
               * "G=U+F polynomial. For the exhaustive symbolic PLD run, re-invoke "
               * "with `julia +1.10 $(PROGRAM_FILE) ...` and a working PLD env.")
        run_hc_fallback(specpath, outpath; check = check, seed = seed, tol = tol)
        return
    end

    spec, edges, nodes, im, em = load_spec(specpath)
    name = get(spec, "name", "graph")
    println("[$name] calling getPLD(edges=$(edges), nodes=$(nodes); method=$(method)) ...")
    t0 = time()

    # PLD returns (discriminants, pars, vars, U, F) in recent versions; older
    # versions return only the discriminant list.  Handle both.  If the call
    # ERRORS (Oscar/HC blowup, OOM), we MUST NOT vanish: write an honest
    # incomplete record so the run is never read as "covered everything".
    faces = Any[]
    alpha = String[]
    complete = true
    errmsg  = ""
    try
        res = getPLD(edges, nodes; internal_masses = im, external_masses = em,
                     method = method, verbose = true,
                     codim_start = codim_start, face_start = face_start,
                     load_output = load_output, save_output = save_output)
        disc = isa(res, Tuple) ? res[1] : res
        # Flatten + stringify the per-face/weight discriminants.
        for (i, d) in enumerate(disc)
            ls = String[]
            for f in (isa(d, Vector) ? d : [d])
                s = string(f)
                push!(ls, s)
                if !(s in alpha) && s != "1"
                    push!(alpha, s)
                end
            end
            push!(faces, Dict("face_id" => i, "letters" => ls,
                              "status" => "resolved", "face_status" => "done"))
        end
    catch e
        complete = false
        errmsg = string(e)
        @warn "getPLD failed; writing HONEST incomplete record" exception=(e, catch_backtrace())
    end

    out = Dict(
        "name" => name,
        "engine" => "pld_bridge.jl -> PLD.jl getPLD (Fevola-Mizera-Telen 2311.14669)",
        "method" => string(method),
        "graph" => spec,
        # TOP-LEVEL honesty flag: FALSE if getPLD did not return a full result.
        "complete" => complete,
        "alphabet" => alpha,
        "first_entry" => String[],   # PLD does not classify; left for the sympy backend
        "last_entry"  => String[],
        "spurious"    => String[],
        "faces" => faces,
        "honesty" => Dict(
            "complete" => complete,
            "n_faces" => length(faces),
            "resolved" => length(faces),
            "bounded_partial" => 0,
            "bounded_timeout" => 0,
            "unresolved_faces" => Any[],
            "global_deadline_hit" => false,
            "error" => errmsg,
            "secs" => round(time() - t0; digits = 2),
            "note" => complete ?
                "PLD.jl is exhaustive on the principal A-determinant; no per-face timeout." :
                "getPLD did NOT complete (see honesty.error); alphabet is PARTIAL — re-run / use the sympy backend per-face.",
        ),
    )
    open(outpath, "w") do io
        JSON.print(io, out, 2)
    end
    if complete
        println("written -> $outpath  ($(length(alpha)) letters)")
    else
        println("written -> $outpath  INCOMPLETE (getPLD failed: $errmsg)")
        exit(2)
    end
end

main()
