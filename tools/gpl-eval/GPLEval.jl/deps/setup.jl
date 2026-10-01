#!/usr/bin/env julia
# One-time wiring for a fresh checkout:
#
#     julia GPLEval.jl/deps/setup.jl
#
# Develops the two in-repo dependencies into this project — LandauAlphabet
# (tools/landau-alphabet/LandauAlphabet.jl) and Eichler (upgrades/Eichler.jl) —
# and instantiates the registered deps (Arblib, JSON3).  Override the dependency
# locations with LANDAUALPHABET_JL_PATH / EICHLER_JL_PATH if you keep them
# outside this checkout.

using Pkg

const PKG  = normpath(joinpath(@__DIR__, ".."))            # GPLEval.jl/
const REPO = normpath(joinpath(PKG, "..", "..", ".."))     # repository root

const LB = get(ENV, "LANDAUALPHABET_JL_PATH",
               joinpath(REPO, "tools", "landau-alphabet", "LandauAlphabet.jl"))
const EI = get(ENV, "EICHLER_JL_PATH",
               joinpath(REPO, "upgrades", "Eichler.jl"))

for (name, envvar, path) in (("LandauAlphabet", "LANDAUALPHABET_JL_PATH", LB),
                             ("Eichler", "EICHLER_JL_PATH", EI))
    isfile(joinpath(path, "Project.toml")) ||
        error("$name not found at $path — set $envvar to its package directory")
end

Pkg.activate(PKG)
# Record the in-repo packages by RELATIVE path so the Manifest survives cloning
# the repository to any location; paths outside the repo stay as given.
relify(p) = startswith(normpath(p), REPO) ? relpath(normpath(p), PKG) : p
cd(PKG) do
    Pkg.develop([Pkg.PackageSpec(path = relify(LB)),
                 Pkg.PackageSpec(path = relify(EI))])
    Pkg.instantiate()
end
Pkg.status()
