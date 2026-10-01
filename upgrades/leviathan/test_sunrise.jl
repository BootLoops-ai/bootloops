# Gate 2 control: equal-mass sunrise (paper worked example).
# Expect chi_top(nu=0) = 4; Landau factors {s, mm, s-mm, s-9mm};
# nu=0 top-sector elimination: x0 cubic with c3 ~ 16*s, x1 quartic with
# c4 ~ 9*mm*(mm-s)^2*(9mm-s)  (factors {mm, s-mm, s-9mm}).
include(joinpath(@__DIR__, "src", "Leviathan.jl"))
using .Leviathan, Random, Oscar

t0 = time()
rng = MersenneTwister(20260630)
fam = family("x1*x2 + x1*x3 + x2*x3 + s*x1*x2*x3 - mm*(x1+x2+x3)*(x1*x2+x1*x3+x2*x3)",
             ["s", "mm"], ["x1", "x2", "x3"])
println("sectors: ", admissible_sectors(fam))

r = landau(fam; rng=rng)
print(report(r))
@assert r.chi_sectors[[1, 2, 3]] == 4 "chi_top != 4"
expected = Set(["s", "mm", "s - mm", "s - 9*mm"])
@assert Set(r.genuine) == expected "genuine mismatch: $(r.genuine)"
@assert isempty(r.spurious)

# paper's leading-coefficient CONTENT (nu=0, top sector). NOTE: degrees differ from the
# paper's (3,4): at nu=0 the x2<->x3 symmetry makes critical points share coordinates, so
# the exact generator of J_top ∩ K[x_i] is shorter; its lc carries the same factors
# (x0: s ~ 16s; x1: mm*(s-mm)^2*(s-9mm) ~ 9m^2(m^2-s)^2(9m^2-s)), which is the payload.
facs, info = sector_candidates(fam, [1, 2, 3]; regulated=false, rng=rng)
println("nu0 top-sector elim degrees: ", Dict(k => v.m for (k, v) in info))
println("  x0: m=", info["x0"].m, "  lc=", info["x0"].lc)
println("  x1: m=", info["x1"].m, "  lc=", info["x1"].lc)
lcfacs(v) = Set(string(f * inv(Oscar.leading_coefficient(f)))
                for (f, _) in Oscar.factor(Oscar.numerator(info[v].lc)) if !Oscar.is_constant(f))
@assert lcfacs("x0") == Set(["s"]) "x0 lc factors: $(lcfacs("x0"))"
@assert lcfacs("x1") == Set(["mm", "s - mm", "s - 9*mm"]) "x1 lc factors: $(lcfacs("x1"))"

d = diagnose(fam; rng=rng, candidates=r.candidates)
for (k, v) in sort(collect(d), by=first)
    println(k, " => ", v[1] ? "PASS" : "FAIL")
end
@assert all(v[1] for (_, v) in d)
println("SUNRISE PASS in ", round(time() - t0, digits=1), " s (+ Oscar load)")
