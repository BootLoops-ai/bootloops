# Gate 1 smoke: one-loop equal-mass bubble. Expect chi=3, factors {s, mm, s-4mm},
# additivity 3 = 1+1+1, all diagnostics PASS.
include(joinpath(@__DIR__, "src", "Leviathan.jl"))
using .Leviathan, Random

t0 = time()
rng = MersenneTwister(20260630)
fam = family("x1 + x2 + mm*(x1+x2)^2 - s*x1*x2", ["s", "mm"], ["x1", "x2"])

println("sectors: ", admissible_sectors(fam))
r = landau(fam; rng=rng)
print(report(r))
@assert r.chi_generic == 3
@assert r.chi_sectors[[1, 2]] == 1 && r.chi_sectors[[1]] == 1 && r.chi_sectors[[2]] == 1
expected = Set(["s", "mm", "s - 4*mm"])
got = Set(r.genuine)
println("expected genuine: ", expected)
@assert got == expected "genuine mismatch: got $got"

d = diagnose(fam; rng=rng, candidates=r.candidates)
for (k, v) in sort(collect(d), by=first)
    println(k, " => ", v[1] ? "PASS" : "FAIL", "   ", v[2])
end
@assert d["type1.1_sector_dim"][1]
@assert d["type2.1_chi_additivity"][1]
# type1.2 correctly FAILS on bubble: s=0 is nu=0-degenerate (second-type, regulated-lane-only).
@assert !d["type1.2_nodegeneracy"][1] &&
        all(occursin("s", g) for (_, g) in d["type1.2_nodegeneracy"][2].offenders)

println("SMOKE PASS in ", round(time() - t0, digits=1), " s (+ Oscar load)")
