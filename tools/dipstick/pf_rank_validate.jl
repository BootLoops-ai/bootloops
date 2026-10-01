# DIPSTICK member — `order` validation harness (timed pilot).
# pf_rank_validate.jl — validation harness (timed pilot).
include(joinpath(@__DIR__, "pf_rank.jl"))

# equal-mass sunrise (LP G, mm=1; sign of s matches test_sunrise.jl)
G_sun = "x1*x2+x1*x3+x2*x3 + s*x1*x2*x3 - (x1+x2+x3)*(x1*x2+x1*x3+x2*x3)"
# 3-loop banana equal-mass (K3 control): U=σ3, F=(Σx)σ3 - s·∏x
G_ban3 = "(x1*x2*x3+x1*x2*x4+x1*x3*x4+x2*x3*x4) + " *
         "(x1+x2+x3+x4)*(x1*x2*x3+x1*x2*x4+x1*x3*x4+x2*x3*x4) - s*x1*x2*x3*x4"
# 4-loop banana equal-mass (CY3): U=σ4, F=(Σx)σ4 - s·∏x
G_ban4 = "(x1*x2*x3*x4+x1*x2*x3*x5+x1*x2*x4*x5+x1*x3*x4*x5+x2*x3*x4*x5) + " *
         "(x1+x2+x3+x4+x5)*(x1*x2*x3*x4+x1*x2*x3*x5+x1*x2*x4*x5+x1*x3*x4*x5+x2*x3*x4*x5)" *
         " - s*x1*x2*x3*x4*x5"
# 3-loop ice-cone: G=U+F
U_icc = "x1*x2*x3+x1*x2*x4+x1*x2*x5+x1*x3*x4+x1*x3*x5+x2*x3*x4+x2*x3*x5"
F_icc = "-s*x4*x5*(x1*x2+x1*x3+x2*x3) + (x1+x2+x3+x4+x5)*(" * U_icc * ")"
G_icc = "(" * U_icc * ") + (" * F_icc * ")"

targets = [
    ("sunrise-eq (elliptic)",  G_sun,  ["s"], ["x1","x2","x3"],                "s", 3, 2),
    ("3L-banana-eq (K3)",      G_ban3, ["s"], ["x1","x2","x3","x4"],           "s", 4, 3),
    ("3L-ice-cone (2×sunrise)",G_icc,  ["s"], ["x1","x2","x3","x4","x5"],      "s", 5, 4),
    ("4L-banana-eq (CY3)",     G_ban4, ["s"], ["x1","x2","x3","x4","x5"],      "s", 5, 4),
]

results = []
for (nm,G,kn,xn,sv,Lm,exp_ord) in targets
    print(rpad(nm,26)); flush(stdout)
    t0=time()
    r = try pf_probe(G,kn,xn,sv;Lmax=Lm) catch e; (error=string(e),) end
    dt = round(time()-t0,digits=1)
    if haskey(r,:error)
        println("  ERROR: ",r.error); push!(results,(nm,:error,dt)); continue
    end
    ok = r.pf_order==exp_ord ? "PASS" : "FAIL"
    println("  hr=",r.holonomic_rank,"  pf_order=",r.pf_order," (exp ",exp_ord,") ",
            ok,"  2p=",r.pf_order_2p,"  red=",r.reducible,
            "  [",dt,"s; W=",r.nW_ncol,"]")
    push!(results,(nm,r,dt,exp_ord,ok))
    dt > 300 && (println("  >5min — stopping (time budget)"); break)
end
println("\nDONE")
