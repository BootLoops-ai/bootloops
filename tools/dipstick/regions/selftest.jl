#!/usr/bin/env julia
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision.
# selftest.jl — the GATES control pair for the dipstick `regions` member
# (regions.jl, the expansion-by-regions completeness certifier), one command.
#
# Runs the certifier on three shipped inputs: the textbook large-Q bubble
# (tests/bubble_largeQ.json — three regions; must give "TILES OK", "CHI
# ADDITIVITY OK" and "VERDICT: PASS"), the positive control
# (tests/largemass_toy_complete.json — the full region list must give "TILES OK"
# and "VERDICT: PASS") and the negative control
# (tests/largemass_toy_incomplete.json — the region-I-only list must give
# "VERDICT: FAIL" and NAME the missing stratum in a "MISSING facets" line).
# All three run in ONE child julia process (one Oscar load, ~50 s); the reports
# are split on their "== dipstick regions: <file> ==" headers.
#
# Engine gate: regions.jl needs the Julia packages Oscar and JSON. This selftest
# looks for them first in the active environment, then in the project named by
# DIPSTICK_JULIA_PROJECT (if set), then in the committed project at
# ../../../upgrades/leviathan (the chi engine's own environment, which is
# exactly this stack). When none provides them it SKIPS by name — exit 0, with
# the exact instantiate command printed — the designed engine-absent behavior,
# not breakage. The same text runs from `python3 tools/dipstick/battery.py
# --legs regions` (which adds the no-julia named skip).
# Exit 0 = both controls behaved as documented, or a named skip;
# exit 1 = a genuine control failure, or a fail-closed refusal when the
# environment claims the packages but they fail to load.

const HERE = @__DIR__
const LEV_PROJECT = normpath(joinpath(HERE, "..", "..", "..", "upgrades", "leviathan"))
const ENV_PROJECT = get(ENV, "DIPSTICK_JULIA_PROJECT", "")
const NEEDED = ("Oscar", "JSON")
const HEADER = "== dipstick regions: "

# --- engine probe (cheap: package lookup only, nothing is loaded) -----------

function missing_in_active_env()
    for name in NEEDED
        id = Base.identify_package(name)
        (id === nothing || Base.locate_package(id) === nothing) && return name
    end
    return nothing
end

function project_provides(proj::String)
    isempty(proj) && return false
    isfile(joinpath(proj, "Project.toml")) || isfile(proj) || return false
    code = """
    for name in ("Oscar", "JSON")
        id = Base.identify_package(name)
        (id === nothing || Base.locate_package(id) === nothing) && exit(3)
    end
    """
    cmd = `$(Base.julia_cmd()) --startup-file=no --project=$proj -e $code`
    return success(run(pipeline(ignorestatus(cmd); stdout=devnull, stderr=devnull)))
end

function named_skip(missing_name)
    println("""
    SKIP (named engine gate): Julia package $missing_name is not available in this
    environment — the dipstick `regions` member (regions.jl) needs Oscar and JSON
    (neither the active Julia environment, nor DIPSTICK_JULIA_PROJECT, nor the
    committed project at ../../../upgrades/leviathan provides them). This skip is
    the designed engine-absent behavior for an engine-gated leg, not breakage.
    Nothing was verified.

    To set up from this clone (one time; installs Oscar + JSON from the General
    registry — Oscar is a large download):

        julia --project=$LEV_PROJECT -e 'import Pkg; Pkg.instantiate()'

    then rerun:  julia selftest.jl   (from this directory)
            or:  python3 tools/dipstick/battery.py --legs regions   (repo root)""")
    return 0
end

# --- controls ---------------------------------------------------------------

const LOAD_FAIL_SIGNS = ("not found in current path",
                         "does not seem to be installed",
                         "Failed to precompile",
                         "error while loading shared libraries")

function run_cases(projflag, json_names)
    args = String[]
    projflag !== nothing && push!(args, "--project=$projflag")
    push!(args, joinpath(HERE, "regions.jl"))
    for j in json_names
        push!(args, joinpath(HERE, "tests", j))
    end
    cmd = `$(Base.julia_cmd()) --startup-file=no $args`
    out = IOBuffer()
    p = run(pipeline(ignorestatus(cmd); stdout=out, stderr=out))
    return p.exitcode, String(take!(out))
end

# split the combined report into one section per input file (keyed by basename)
function sections(s::String)
    d = Dict{String,String}()
    cur = nothing; buf = IOBuffer()
    for rawline in split(s, '\n')
        line = strip(rawline, ['\r'])   # Oscar's banner suppression leaves a stray CR
        if startswith(line, HEADER)
            cur !== nothing && (d[cur] = String(take!(buf)))
            fname = strip(replace(line[length(HEADER)+1:end], r"==\s*$" => ""))
            cur = basename(String(fname)); buf = IOBuffer()
        end
        println(buf, line)
    end
    cur !== nothing && (d[cur] = String(take!(buf)))
    return d
end

function main()
    # pick the environment the certifier will run under; the child process does
    # not inherit a --project flag, so pass the chosen project explicitly
    projflag = nothing
    miss = missing_in_active_env()
    if miss === nothing
        projflag = dirname(Base.active_project())
    elseif project_provides(ENV_PROJECT)
        projflag = ENV_PROJECT
    elseif project_provides(LEV_PROJECT)
        projflag = LEV_PROJECT
    else
        return named_skip(miss)
    end

    cases = (("bubble_largeQ.json", :pass_chi),
             ("largemass_toy_complete.json", :pass),
             ("largemass_toy_incomplete.json", :fail))
    rc, s = run_cases(projflag, [c[1] for c in cases])
    if rc != 0 && any(sig -> occursin(sig, s), LOAD_FAIL_SIGNS)
        # fail-closed by name: the environment offered the packages but the
        # certifier could not LOAD them (this is not a control failure).
        println("ENGINE GATE (fail-closed): regions.jl could not load its Julia ",
                "dependencies (Oscar + JSON) even though the environment lists ",
                "them. Re-instantiate the environment, e.g.")
        println("    julia --project=$LEV_PROJECT -e 'import Pkg; Pkg.instantiate()'")
        println("engine said:")
        println(join(last(split(strip(s), '\n'), 2), '\n'))
        return 1
    end
    if rc != 0
        println("FAIL: regions.jl exited rc=", rc)
        println(join(last(split(strip(s), '\n'), 8), '\n'))
        return 1
    end
    sec = sections(s)
    ok = true
    for (json_name, want) in cases
        t = get(sec, json_name, "")
        if isempty(t)
            println("control FAIL: ", json_name, " -> no report section in the output")
            ok = false
            continue
        end
        if want == :pass_chi
            good = occursin("VERDICT: PASS (region list complete)", t) &&
                   occursin("[TILES OK]", t) && occursin("[CHI ADDITIVITY OK]", t)
            println(good ? "control ok: " : "control FAIL: ", json_name,
                    " -> ", good ? "TILES OK / CHI ADDITIVITY OK / VERDICT: PASS" :
                    "missing one of TILES OK / CHI ADDITIVITY OK / VERDICT: PASS")
            ok &= good
        elseif want == :pass
            good = occursin("VERDICT: PASS (region list complete)", t) &&
                   occursin("[TILES OK]", t)
            println(good ? "control ok: " : "control FAIL: ", json_name,
                    " -> ", good ? "TILES OK / VERDICT: PASS" : "no PASS verdict")
            ok &= good
        else
            good = occursin("VERDICT: FAIL (region list INCOMPLETE)", t) &&
                   occursin("MISSING facets", t)
            println(good ? "control ok: " : "control FAIL: ", json_name,
                    " -> ", good ? "VERDICT: FAIL with the missing facet named" :
                    "did not fail-and-name the missing facet")
            ok &= good
        end
    end
    println(ok ? "dipstick regions selftest PASS (bubble + positive + negative control)" :
                 "dipstick regions selftest FAIL")
    return ok ? 0 : 1
end

exit(main())
