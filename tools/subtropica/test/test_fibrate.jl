# test_fibrate.jl — (DESIGN_B1.md §2).  Tests for the
# scripts/fibrate/ SOUND fibration engine (gap-tool replacing HyperFLINT's
# unsound `fibration_basis` op, which is refused).
#
# Self-contained: shells out to the python CLI; needs NO subtropica src
# includes and NO Oscar (DESIGN_B1 §3 include-by-path pattern applies to
# files using the shared Julia types — this one doesn't).  Runs green
# standalone in the DEFAULT global Julia env:
#   ulimit -v 32505856; julia \
#     tools/subtropica/test/test_fibrate.jl
# and joins the shared suite automatically via test/runtests.jl.
#
# Every python invocation is prefixed `ulimit -v 32505856` (package convention).
#
# What is checked (hardening requirements):
#  (1) selfcheck: PSLQ control battery (positive + negative control) and
#      the adaptive-ZIP-series footgun guard (fixed-length dps is a
#      footgun: the series length must adapt to the word's minimum letter
#      magnitude) verified against exact log(26/9) truth.
#  (2) regression corpus: 10 words copied from the production pilot's
#      validated 370 (fixtures/regression_expected.json, source cited
#      inside) — exact canonical term equality, mandatory 2-point
#      validation digits, and the ':pslq_constants' provenance stamp.
#  (3) comparator negative control: a perturbed expected-term copy MUST
#      mismatch (print-loud-rc0-gates footgun class — the comparison can
#      fail, therefore its pass means something).
#  (4) typed refusals: on-contour letter, single validation point
#      (validation is mandatory, not tunable-away), malformed JSON.

using Test
using JSON

const _FB_DIR = joinpath(@__DIR__, "..", "scripts", "fibrate")
const _FB_CLI = joinpath(_FB_DIR, "fibrate.py")
const _FB_FIXDIR = joinpath(_FB_DIR, "fixtures")

"""Run `python3 fibrate.py <args>` under the house ulimit; return
(exitcode, stdout, stderr)."""
function _fb_run(args::Vector{String})
    # -B: no __pycache__ pollution inside the owned tool dir
    cmdstr = "ulimit -v 32505856; exec python3 -B " *
             join(vcat([_FB_CLI], args), " ")
    out = IOBuffer(); err = IOBuffer()
    p = run(pipeline(ignorestatus(`bash -c $cmdstr`),
                     stdout=out, stderr=err))
    return p.exitcode, String(take!(out)), String(take!(err))
end

# canonical term comparison (order-sensitive: the CLI and the fixture use
# the SAME canonical serializer, fibrate.serialize_rep)
_fb_terms_equal(a, b) = JSON.json(a) == JSON.json(b)

@testset "fibrate" begin

    @test isfile(_FB_CLI)
    @test isfile(joinpath(_FB_FIXDIR, "regression_words.json"))
    @test isfile(joinpath(_FB_FIXDIR, "regression_expected.json"))

    scratch = mktempdir()

    # ---- (1) selfcheck: control battery + footgun-guard regressions ------
    @testset "selfcheck battery" begin
        rc, out, err = _fb_run(["--selfcheck"])
        @test rc == 0
        rep = JSON.parse(out)
        @test rep["selfcheck"] == "PASS"
        @test rep["pslq_positive_control"] == "pass"
        @test rep["pslq_negative_control"] == "pass"
        # adaptive series length vs EXACT log(26/9): strictly-greater-
        # than-slack margin (>= 120 of 130 digits), not a soft print
        @test rep["adaptive_series_vs_log_truth_digits"] >= 120
        @test rep["geval_crosscheck_digits"] >= 50
    end

    # ---- (2) regression corpus vs pilot-validated fixtures ---------------
    outpath = joinpath(scratch, "regression_out.json")
    expected = JSON.parse(read(joinpath(_FB_FIXDIR, "regression_expected.json"),
                               String))
    doc = nothing
    @testset "regression corpus (10 pilot words)" begin
        rc, out, err = _fb_run([
            joinpath(_FB_FIXDIR, "regression_words.json"), "--out", outpath])
        @test rc == 0
        @test isfile(outpath)
        doc = JSON.parse(read(outpath, String))
        @test doc["tool"] == "subtropica/fibrate"
        @test length(doc["words"]) == 10
        # mandatory validation is declared and 2-point
        @test doc["validation"]["mandatory"] === true
        @test length(doc["validation"]["points"]) >= 2
        # PSLQ control battery report present with the full battery
        pc = doc["pslq_controls"]
        @test pc["positive_control"] == "pass"
        @test length(pc["dual_precision"]) == 2 &&
              pc["dual_precision"][1] != pc["dual_precision"][2]
        @test pc["height_cap"] == 10_000_000_000
        @test haskey(pc, "residual_gate")
        # per-word: exact canonical term equality vs the pilot fixture,
        # >= 85 validated digits at BOTH points, honest provenance stamp
        for (wo, we) in zip(doc["words"], expected["words"])
            @test wo["word"] == we["word"]
            @test _fb_terms_equal(wo["terms"], we["terms"])
            @test length(wo["validated_digits"]) >= 2
            @test minimum(collect(values(wo["validated_digits"]))) >= 85
            if we["own_rep_has_named_constants"]
                @test wo["used_pslq"] === true
                @test "pslq_constants" in wo["provenance"]
            end
            @test wo["used_pslq"] == ("pslq_constants" in wo["provenance"])
        end
        # the corpus contains PSLQ-constant words => top-level stamp
        @test "pslq_constants" in doc["provenance"]
    end

    # ---- (3) comparator NEGATIVE control (non-vacuous comparison) --------
    @testset "comparator negative control" begin
        @test doc !== nothing
        tampered = JSON.parse(read(outpath, String))   # fresh copy
        # perturb one coefficient of the first word that has terms
        idx = findfirst(w -> !isempty(w["terms"]), tampered["words"])
        @test idx !== nothing
        tampered["words"][idx]["terms"][1]["coef"] =
            tampered["words"][idx]["terms"][1]["coef"] * " + 1"
        @test !_fb_terms_equal(tampered["words"][idx]["terms"],
                           expected["words"][idx]["terms"])
    end

    # ---- (4) typed refusals ----------------------------------------------
    @testset "refusals are loud and nonzero-rc" begin
        # on-contour ZIP letter (the exact class where upstream
        # fibration_basis is unsound — must REFUSE, never guess)
        onc = joinpath(scratch, "oncontour.json")
        write(onc, "[[\"ze\"]]")
        rc, out, err = _fb_run([onc, "--out", "/dev/null"])
        @test rc == 2
        @test occursin("FATAL", err)
        @test occursin("on-contour", err)

        # validation is MANDATORY: a single validation point is refused
        okw = joinpath(scratch, "okword.json")
        write(okw, "[[\"-ze\"]]")
        rc2, out2, err2 = _fb_run([okw, "--validation-points", "2/5",
                                       "--out", "/dev/null"])
        @test rc2 == 2
        @test occursin("mandatory", err2)

        # malformed wordlist JSON
        bad = joinpath(scratch, "bad.json")
        write(bad, "not json")
        rc3, _, err3 = _fb_run([bad, "--out", "/dev/null"])
        @test rc3 == 2
        @test occursin("FATAL", err3)
    end
end
