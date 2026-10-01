# builder.jl — graph -> Lee-Pomeransky G via PLD.jl's getUF (reused, not rewritten).
# Heavy deps (HomotopyContinuation) load only when this file is included.
# Usage:
#   include("builder.jl")
#   fam = family_from_graph(edges, nodes; internal_masses=..., external_masses=...)
# edges/nodes follow PLD conventions, e.g. sunrise: edges=[[1,2],[1,2],[1,2]], nodes=[1,2].

# Point PLD_PATH at src/PLD.jl from the PLD distribution
# (Fevola–Mizera–Telen, https://mathrepo.mis.mpg.de/PLD/ — obtain upstream).
const PLD_PATH = get(ENV, "PLD_PATH", "")
isempty(PLD_PATH) && error("builder.jl requires ENV[\"PLD_PATH\"] = /path/to/PLD/src/PLD.jl")
isdefined(Main, :PLD) || include(PLD_PATH)

"""
Build a Leviathan Family from a Feynman graph using PLD's Symanzik builder.
Kinematic parameter names are taken from PLD's output (masses + Mandelstams).
"""
function family_from_graph(edges, nodes; internal_masses=:zero, external_masses=:zero,
                           substitute=[], cut=Int[])
    U, F, pars, vars = Main.PLD.getUF(edges, nodes; internal_masses=internal_masses,
                                      external_masses=external_masses, substitute=substitute)
    kinnames = [string(p) for p in pars]
    xnames = [string(v) for v in vars]
    G = U + F
    # map the Oscar polynomial from PLD's ring into a fresh Leviathan family ring
    fam = family("0", kinnames, xnames; cut=cut)
    Rsrc = Oscar.parent(G)
    imgs = [Oscar.gens(fam.R)[findfirst(==(string(g)), vcat(kinnames, xnames))]
            for g in Oscar.gens(Rsrc)]
    h = Oscar.hom(Rsrc, fam.R, c -> Nemo.QQ(c), imgs)
    Family(fam.R, h(G), fam.kin, fam.xs, fam.kinnames, fam.xnames, fam.cut)
end
