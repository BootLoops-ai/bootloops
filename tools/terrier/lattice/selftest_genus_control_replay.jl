# Battery LAT-G (lattice/genus_enum.jl): replay the CONTROL genus — det 92,
# key 72d805a279f839ff (the positive-definite rank-5 genus with that
# discriminant form) — with the Aut-free Kneser engine + PARI POST, then gate
# on the reference numbers EXACTLY: 6 classes, Hecke mass 53/768, 0 root-free
# classes, verdict CONTROL-PASS.  genus_enum.jl's own control_check additionally
# enforces the qfisom bijection + per-class |Aut| against the reference receipt
# (any miss = throw = exit != 0).
# Needs: julia 1.10 with Oscar (TERRIER_JULIA_PROJECT), PARI/GP `gp` on PATH,
# and the reference receipt of the det-92 control genus, which is not included
# in the package: supply one with TERRIER_GENUS_CONTROL_RECEIPT (an
# out_<key8>.txt in the format genus_enum.jl writes: CLS| lines + HECKEMASS|).
# Run: nice julia +1.10 --project=$TERRIER_JULIA_PROJECT -t 1 selftest_genus_control_replay.jl
# (timeout generous: Oscar load dominates; control enum ~30 s post-load).
const CKEY = "72d805a279f839ff"
ref_rc = get(ENV, "TERRIER_GENUS_CONTROL_RECEIPT", "")  # NB: genus_enum.jl owns `REFRC`
isempty(ref_rc) && error("SKIP[genusenum]: TERRIER_GENUS_CONTROL_RECEIPT unset — the " *
    "reference receipt of the det-92 control genus is not included in the package; " *
    "point the variable at one to run")
isfile(ref_rc) || error("reference control receipt MISSING: $ref_rc")
ENV["TERRIER_GENUS_OUT"] = mktempdir()                 # receipts to scratch, never in-tree
empty!(ARGS); append!(ARGS, ["--control", "92", CKEY, ref_rc])
include(joinpath(@__DIR__, "genus_enum.jl")) # runs main(): enum + POST + control gate
import JSON
rec = JSON.parsefile(joinpath(ENV["TERRIER_GENUS_OUT"], "receipts", "CTRL_" * CKEY * ".json"))
ok = rec["verdict"] == "CONTROL-PASS" &&
     rec["n_classes"] == 6 &&                # reference: 6 classes
     rec["mass"] == "53//768" &&             # reference: exact Hecke mass
     rec["mass_check"] == true &&            # POST sum(1/|Aut|) == mass exactly
     rec["n_rootfree"] == 0                  # reference: all classes rooted
println("SELFTEST genus_control_replay: classes=", rec["n_classes"],
        " mass=", rec["mass"], " rootfree=", rec["n_rootfree"],
        " verdict=", rec["verdict"], ok ? " -- ALL PASS" : " -- FAIL")
exit(ok ? 0 : 1)
