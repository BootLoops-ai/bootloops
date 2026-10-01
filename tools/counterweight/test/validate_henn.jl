# validate_henn.jl — HARD VALIDATION GATE (a): the massless planar double box, Henn 1304.1806.
#
# (1) Henn's canonical A certifies as eps-form + dlog, singular pts {0,-1}, alphabet {x,1+x}.
# (2) Scramble it with a nontrivial eps&x-dependent gauge -> generic non-UT Fuchsian A'.
# (3) Engine recovers an eps-form Atilde from A'.
# (4) Recovered residues are conjugate to Henn's a,b (same eigenvalues -> same canonical basis
#     up to the residual constant freedom).
#
# Run:  julia --project=<env-with-Nemo> test/validate_henn.jl
#       (or both files in one process: julia --project=. test/runtests.jl — the
#        registered battery; the include below is guarded so that reuses one load)

isdefined(Main, :Counterweight) ||
    include(joinpath(@__DIR__, "..", "src", "Counterweight.jl"))
using .Counterweight
using .Counterweight.Normalize, .Counterweight.Fuchsia, .Counterweight.RatFunc
using .Counterweight.HennDbox
using Nemo

# eigenvalue -> multiplicity dict for a QQ matrix (rational eigenvalues only; we only need to
# compare spectra, and Henn's residues have rational spectrum).
function _eig_mult(M)
    cp = charpoly(M)
    d = Dict{Any,Int}()
    for (q,e) in factor(cp)
        if degree(q) == 1
            r = -coeff(q,0)//coeff(q,1)
            d[r] = get(d, r, 0) + e
        else
            d[string(q)] = get(d, string(q), 0) + e  # irreducible block, compare as string
        end
    end
    return d
end

println("="^70)
println("VALIDATION GATE (a): massless planar double box — Henn 1304.1806")
println("="^70)

ctx = make_ctx()
A = HennDbox.henn_canonical(ctx)

println("\n[1] Henn canonical connection: certify eps-form")
cert = certify_epsform(ctx, A)
for k in ("fuchsian","singular_points","epsform","dlog","algebraic_factors")
    println("    $k = ", cert[k])
end
@assert cert["fuchsian"]  "Henn A not Fuchsian?!"
@assert cert["epsform"]   "Henn A not eps-form?!"
@assert cert["dlog"]      "Henn A not dlog?!"
@assert Set(cert["singular_points"]) == Set(["0","-1"]) "wrong singular set: $(cert["singular_points"])"
println("    => PASS: Henn basis is eps-form, dlog, singular {0,-1}, alphabet {x,1+x}.")

# record the eps^0-stripped residues (the symbol alphabet matrices a, b)
ok, B = is_epsform(ctx, A)
res0 = residue_matrix(ctx, B, ctx.Qe(0))
resm1 = residue_matrix(ctx, B, ctx.Qe(-1))
a0 = HennDbox.henn_a(); b0 = HennDbox.henn_b()
println("    residue at x=0 == a ? ", res0 == a0)
println("    residue at x=-1 == b ? ", resm1 == b0)

println("\n[2] Scramble to a generic non-UT Fuchsian basis A'")
Ascr, Bgauge = HennDbox.henn_scramble(ctx, A)
println("    A' Fuchsian? ", is_fuchsian(ctx, Ascr))
println("    A' eps-form? ", is_epsform(ctx, Ascr)[1], "  (expected false)")

println("\n[3] Run engine: reduce A' to eps-form")
t0 = time()
r = try_epsfactor(ctx, Ascr; max_balance_rounds=400)
dt = round(time()-t0, digits=2)
println("    ok = ", r.ok, "   (", dt, " s)")
if !r.ok
    println("    STOP reason: ", get(r.report,"stop","?"))
    println("    history: ", get(r.report,"history","?"))
    error("ENGINE FAILED to reduce the massless double box — VALIDATION GATE (a) FAIL")
end
println("    balance rounds = ", get(r.report,"balance_rounds","-"))
recert = certify_epsform(ctx, r.Anew)
for k in ("fuchsian","singular_points","epsform","dlog")
    println("    recovered $k = ", recert[k])
end
@assert recert["epsform"] && recert["dlog"]  "recovered system not eps-form/dlog"
@assert Set(recert["singular_points"]) == Set(["0","-1"]) "recovered singular set wrong"

println("\n[4] Gauge-equivalence to Henn: residue eigenvalues match")
ok2, Brec = is_epsform(ctx, r.Anew)
rec0  = residue_matrix(ctx, Brec, ctx.Qe(0))
recm1 = residue_matrix(ctx, Brec, ctx.Qe(-1))
ev_a0  = sort(string.(collect(keys(_eig_mult(a0)))))
ev_r0  = sort(string.(collect(keys(_eig_mult(rec0)))))
ev_b0  = sort(string.(collect(keys(_eig_mult(b0)))))
ev_rm1 = sort(string.(collect(keys(_eig_mult(recm1)))))
println("    eig(a)   = ", ev_a0)
println("    eig(res0)= ", ev_r0,  "   match? ", ev_a0 == ev_r0)
println("    eig(b)   = ", ev_b0)
println("    eig(resm1)=", ev_rm1, "   match? ", ev_b0 == ev_rm1)
@assert ev_a0 == ev_r0   "x=0 residue eigenvalues differ from Henn a"
@assert ev_b0 == ev_rm1  "x=-1 residue eigenvalues differ from Henn b"

println("\n" * "="^70)
println("VALIDATION GATE (a): PASS — engine reproduced the Henn massless double box")
println("eps-form (Atilde eps-independent, dlog, alphabet {x,1+x}) from a scrambled basis.")
println("="^70)
