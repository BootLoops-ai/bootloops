#!/usr/bin/env julia
# counterweight_bridge.jl — JSON-in / JSON-out shim around the in-house Counterweight
# engine (Lee 1411.0911 over Nemo).  Called by canonical_form.wrap_epsfactor().
#
#   stdin / argv[1] : {"action": "certify"|"epsform",
#                      "n": N, "x": "x", "eps": "eps",
#                      "A": [["c11","c12",...], ...]}   — entries are rational
#                      functions in (x, eps) as plain strings.
#   stdout          : {"ok": bool, "fuchsian": bool, "epsform": bool, "dlog": bool,
#                      "singular_points": [...],
#                      "T": [[...]], "Atilde": [[...]], "Anew": [[...]],
#                      "report": {...}}
#
# This bridge is intentionally thin; the heavy lifting is in
# the Counterweight engine one directory up (src/{Fuchsia,Normalize}.jl), which has
# its own validation gate (test/validate_henn.jl, massless planar double box).

import JSON
include(joinpath(get(ENV, "COUNTERWEIGHT_ROOT",
                     normpath(joinpath(@__DIR__, ".."))),
                 "src", "Counterweight.jl"))
using .Counterweight
using .Counterweight.RatFunc, .Counterweight.Fuchsia, .Counterweight.Normalize
using .Counterweight.KiraDE: parse_kira_coeff
using Nemo

function _parse_mat(ctx, rows)
    n   = length(rows)
    kin = Dict{Symbol,Any}(:eps => ctx.Qex(ctx.eps), :x => ctx.x,
                           :d   => ctx.Qex(4 - 2*ctx.eps))
    A   = zero_matrix(ctx.Qex, n, n)
    for i in 1:n, j in 1:n
        s = replace(String(rows[i][j]), "**" => "^")
        A[i, j] = parse_kira_coeff(ctx, s, kin)
    end
    return A
end

function _str_mat(M)
    n, m = size(M)
    [[string(M[i, j]) for j in 1:m] for i in 1:n]
end

function main()
    path = length(ARGS) >= 1 ? ARGS[1] : nothing
    inp  = path === nothing ? JSON.parse(read(stdin, String)) :
                              JSON.parsefile(path)
    n    = Int(inp["n"])
    ctx  = make_ctx()                # Q(eps)(x); generator names fixed to eps,x
    A    = _parse_mat(ctx, inp["A"])
    act  = get(inp, "action", "epsform")

    cert = certify_epsform(ctx, A)
    out  = Dict{String,Any}(
        "fuchsian"        => cert["fuchsian"],
        "epsform"         => cert["epsform"],
        "dlog"            => get(cert, "dlog", false),
        "singular_points" => get(cert, "singular_points", String[]),
    )

    if act == "certify"
        out["ok"] = cert["epsform"]
        println(JSON.json(out)); return
    end

    if cert["epsform"]
        ok, B = is_epsform(ctx, A)
        out["ok"]     = true
        out["Atilde"] = _str_mat(B)
        out["Anew"]   = _str_mat(A)
        out["T"]      = _str_mat(identity_matrix(ctx.Qex, n))
        out["report"] = Dict("note" => "already eps-form")
        println(JSON.json(out)); return
    end

    r = try_epsfactor(ctx, A; max_balance_rounds=400)
    out["ok"]     = r.ok
    out["report"] = r.report
    if r.ok
        out["T"]      = _str_mat(r.T)
        out["Anew"]   = _str_mat(r.Anew)
        ok, B = is_epsform(ctx, r.Anew)
        out["Atilde"] = _str_mat(B)
        rec = certify_epsform(ctx, r.Anew)
        out["singular_points"] = get(rec, "singular_points", String[])
        out["dlog"]            = get(rec, "dlog", false)
    end
    println(JSON.json(out))
end

main()
