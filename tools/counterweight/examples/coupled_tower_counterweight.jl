# coupled_tower_counterweight.jl — eps-factorize a coupled-DE tower read from a Kira
# derive_dgl export.
#
# Input : a JSON export of Kira derive_dgl rows (d/dv M_i = sum_j A_v[i,j] M_j, coeff strings in d,s,t)
# Output: rotation_T.json  (rotation T raw->UT, eps-form Atilde, alphabet, certificate)
#
# Pipeline (Lee balance + Magnus, single ratio variable x = t/s at a probe s; gauge transported):
#   * Build A_t over Q(eps)(x) with d=4-2eps, s=sprobe, t=x  (A_x = sprobe*A_t at fixed s; the
#     overall scale weight enters only the diagonal e^{a eps log(-s)} which is the κ=1 strip and
#     is x-independent, so it does not affect the x-balances).
#   * Run try_epsfactor -> eps-form Atilde(x), gauge T(eps,x).
#   * Certify eps-form + dlog; report the alphabet (x-letters), confirm the engine's T reproduces
#     eps-form at a SECOND probe s (s-independence of the dlog letters after κ=1 strip).
#   * Export T and Atilde as coefficient strings in (eps, x) for the Python sample-applier.
#
# No Wolfram, no Sage. Pure Julia + Nemo (Flint).

include(joinpath(@__DIR__, "..", "src", "Counterweight.jl"))
using .Counterweight
using .Counterweight.Normalize, .Counterweight.Fuchsia, .Counterweight.RatFunc, .Counterweight.KiraDE
using Nemo
import JSON

# NOTE: no input DE ships with this example — it is a driver pattern.  Point
# TOWER_DE_JSON at your own Kira derive_dgl export to run it.
const DE  = get(ENV, "TOWER_DE_JSON", "de_tower.json")
const OUT = get(ENV, "TOWER_OUT_JSON", "rotation_T.json")

parse_de() = JSON.parsefile(DE)

# build the v-derivative connection over Qex with d=4-2eps, s=sprobe (rational), t=x
function build_Av(ctx::Ctx, de, which::String, sprobe::Int)
    n = de["n"]
    A = zero_matrix(ctx.Qex, n, n)
    kin = Dict{Symbol,Any}(:d => ctx.Qex(4) - 2*ctx.Qex(ctx.eps),
                           :s => ctx.Qex(sprobe),
                           :t => ctx.x)
    blk = de[which]
    for (key, cstr) in blk
        ij = split(key, ",")
        i = parse(Int, ij[1]) + 1; j = parse(Int, ij[2]) + 1
        A[i,j] += KiraDE.parse_kira_coeff(ctx, cstr, kin)
    end
    return A
end

# stringify a Qex matrix as nested coeff strings in (eps, x): entry -> Julia/Nemo string
function mat_to_strings(A)
    m, n = size(A)
    out = Vector{Vector{String}}()
    for i in 1:m
        row = String[]
        for j in 1:n
            push!(row, string(A[i,j]))
        end
        push!(out, row)
    end
    return out
end

function run_probe(de, sprobe::Int; rounds=600)
    ctx = make_ctx()
    At = build_Av(ctx, de, "A_t", sprobe)
    fu = is_fuchsian(ctx, At)
    pts, alg = singular_points_x(ctx, At)
    r = try_epsfactor(ctx, At; max_balance_rounds=rounds)
    return ctx, At, r, fu, sort(string.(pts)), string.(alg)
end

function main()
    de = parse_de()
    n = de["n"]
    println("tower n=$n; top-sector indices (0-based) = ", de["top_indices"])

    sprobe = -1
    ctx, At, r, fu, pts, alg = run_probe(de, sprobe)
    println("\n[1] s=$sprobe  Fuchsian(A_t in x)=$fu  finite x-singularities=$pts  alg-factors=$(length(alg))")
    println("    engine ok = ", r.ok)
    rep = r.report
    if !r.ok
        println("    STOP: ", get(rep,"stop","?"))
        println("    history (last): ", length(get(rep,"history",String[]))>0 ? last(rep["history"]) : "-")
        println("    balance_rounds: ", get(rep,"balance_rounds","-"))
    end

    result = Dict{String,Any}()
    result["n"] = n
    result["masters"] = de["masters"]
    result["top_indices"] = de["top_indices"]
    result["fuchsian_At"] = fu
    result["x_singularities"] = pts
    result["algebraic_factors"] = alg
    result["sprobe"] = sprobe

    if r.ok
        cert = certify_epsform(ctx, r.Anew)
        println("\n[2] RECOVERED eps-form certificate:")
        for k in ("fuchsian","singular_points","epsform","dlog","algebraic_factors")
            println("    $k = ", cert[k])
        end
        result["epsform"] = cert["epsform"]
        result["dlog"] = cert["dlog"]
        result["alphabet_x"] = cert["singular_points"]
        result["recovered_alg_factors"] = cert["algebraic_factors"]
        result["balance_rounds"] = get(rep, "balance_rounds", -1)
        result["history"] = get(rep, "history", String[])
        result["T"] = mat_to_strings(r.T)
        # Atilde = Anew/eps as eps-free dlog matrix
        ok, B = is_epsform(ctx, r.Anew)
        result["Atilde"] = mat_to_strings(B)

        # second-probe consistency: recovered alphabet must match at s=-2
        ctx2, At2, r2, _, pts2, _ = run_probe(de, -2)
        if r2.ok
            cert2 = certify_epsform(ctx2, r2.Anew)
            result["consistency_s2"] = Dict("ok"=>r2.ok, "alphabet_x"=>cert2["singular_points"],
                                            "match"=>Set(cert2["singular_points"])==Set(cert["singular_points"]))
            println("\n[3] second-probe s=-2 alphabet match: ", result["consistency_s2"]["match"])
        else
            result["consistency_s2"] = Dict("ok"=>false)
        end
    else
        result["epsform"] = false
        result["stop"] = get(rep, "stop", "?")
        result["balance_rounds"] = get(rep, "balance_rounds", -1)
        result["history"] = get(rep, "history", String[])
    end

    JSON.print(open(OUT,"w"), result)
    println("\nwrote ", OUT)
    return r
end

main()
