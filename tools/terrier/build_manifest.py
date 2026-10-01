#!/usr/bin/env python3
"""Rebuild regression_manifest.json (sha-pins + battery registry). Run ONLY
after a full green selftest.py run; run-time receipts (test_receipts.json,
RUN_pins.json, SELFTEST_SWEEPCHASSIS.json, BATTERY_RECEIPT.json,
SMOKE_RECEIPT.json, out_selftest/) are regenerated and
checked by the batteries themselves and are deliberately NOT pinned.

The manifest pins the SHIPPED FIXTURES: the toy/known-answer data the
batteries check real results against. Batteries whose inputs are reference
data sets not included in the package stay registered — selftest.py prints
a loud NAMED SKIP naming the env var that supplies each one. The output is
deterministic: same tree bytes, same manifest."""
import glob, hashlib, json, os
ROOT = os.path.dirname(os.path.abspath(__file__))

sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()

art = sorted(glob.glob("periods/pipeline/cards/*.json", root_dir=ROOT)) + [
    "periods/pipeline/dkmm/operator_LS.json",
    "periods/pipeline/dkmm/series_curve.json"] + sorted(
    glob.glob("lattice/emit/bank/*", root_dir=ROOT)) + sorted(
    glob.glob("periods/pipeline/conipfv/out/*.json", root_dir=ROOT)) + [
    # ffp reference tables and pass-sets
    "periods/pipeline/ffp/nf_14_4_a_a.json",
    "periods/pipeline/ffp/nf_34_2_b_a.json",
    "periods/pipeline/ffp/scale/nf_180_4_a_e.json"] + sorted(
    glob.glob("periods/pipeline/ffp/scale/passsets/*.json",
              root_dir=ROOT)) + sorted(
    glob.glob("periods/pipeline/ffp/banks/*.jsonl", root_dir=ROOT)) + sorted(
    glob.glob("periods/pipeline/sunit/atlas/*.json", root_dir=ROOT)) + [
    "common/data/KNOWN_ORBIT.jsonl",
    "census/dedup/TEST_SEEDS.json"] + [
    f"periods/pipeline/cards/hv4-diag-L5_bank/{f}" for f in
    ("FRAME_RANK5.json", "SERIES_A.json", "SERIES_B.json",
     "BATTERY_RECEIPT_reference.json")]

B = lambda i, w, p, r="python", t=1800: dict(id=i, wing=w, path=p, runner=r,
                                             timeout_s=t)
bats = [B("cards", "periods", "periods/selftest_cards.py", t=300),
        B("geffseries", "periods", "periods/selftest_geffseries.py"),
        B("fluxcurves", "periods", "periods/selftest_fluxcurves.py"),
        B("opderive", "periods", "periods/selftest_opderive.py"),
        B("transport", "periods", "periods/selftest_transport.py"),
        B("verdicts", "periods", "periods/selftest_verdicts.py", t=700),
        B("envelope", "periods", "periods/selftest_envelope.py", t=600),
        B("receipts", "common", "common/selftest_receipts.py", t=300),
        B("controls", "common", "common/selftest_controls.py", t=300),
        B("frames", "common", "common/selftest_frames.py", t=300),
        B("certs", "common", "common/selftest_certs.py", t=300),
        B("verdict", "common", "common/selftest_verdict.py", t=300),
        B("lattice_L1", "lattice", "lattice/selftest_lattice.jl", "julia", 1200),
        B("nikulin_L1b", "lattice", "lattice/selftest_nikulin.jl", "julia", 1200),
        B("emit_L2", "lattice", "lattice/emit/selftest_emit.py", t=300),
        B("sweepchassis_L2b", "lattice", "lattice/emit/selftest_sweepchassis.py", t=600),
        B("pfaffian", "periods", "periods/selftest_pfaffian.py", t=600),
        B("conifold", "periods", "periods/selftest_conifold.py", t=300),
        B("ffp", "periods", "periods/selftest_ffp.py", t=600),
        B("sunit", "periods", "periods/selftest_sunit.py", t=900),
        B("dedup", "common", "common/selftest_dedup.py", t=300),
        # census wing batteries
        B("foldhnf", "census", "census/selftest_foldhnf.py", t=1800),
        B("dedup_census", "census", "census/selftest_dedup_census.py", t=900),
        B("distcert", "census", "census/selftest_distcert.py", t=900),
        B("controlsampler", "census", "census/selftest_controlsampler.py",
          t=900),
        B("census_smoke", "census", "census/selftest_census_smoke.py", t=300),
        # wing-map coverage wrappers (names track the wing_regressions
        # keys; ids stay the map's battery ids).
        B("summation", "periods",
          "periods/selftest_summation_reduced_u1s.py", t=600),
        B("slicecard", "periods",
          "periods/selftest_slice_operator_card.py", t=300),
        B("lpkill", "lattice",
          "lattice/selftest_lpkill_a1n_mustfail.py", t=1500),
        B("genusenum", "lattice",
          "lattice/selftest_genus_control_replay.jl", "julia", 3600),
        ]

wreg = {
 "dkmm_ball_reproduction": {"batteries": ["transport"], "gate":
  "route R 150 + vac 150 R re-run; W0/W0_vac ball midpoints string-EXACT "
  "(matched_digits 999) vs the reference data (TERRIER_KKLT_BANK); "
  "radius ratio in [0.9,1.12]"},
 "lattice_18_18": {"batteries": ["lattice_L1"], "gate":
  "validate.jl full battery 18/18 PASS in situ (writes its results report)"},
 "emit_glue_predicates": {"batteries": ["emit_L2"], "gate":
  "E1 44-row known-answer glue table exact (dual-integrality + trace + "
  "det-chain identities, exact arithmetic); E2 closed-form A2xA2 case; "
  "E3 tampered-glue + tampered-value mutations caught"},
 "reference_shard_replay": {"batteries": ["sweepchassis_L2b"], "gate":
  "SC3 reference-shard replay byte-identical (pinned clock); SC7 20/20 "
  "check_ids byte-compatible with the reference shard receipt "
  "(TERRIER_SWEEP_SHARD_RECEIPT)"},
 "planted_error_rc": {"batteries": ["sweepchassis_L2b"], "gate":
  "planted-error rc=8 (control + integrity legs); halt law rc=9 + clear+resume"},
 "mutation_control_trio": {"batteries": ["fluxcurves", "verdicts", "emit_L2"],
  "gate": "f3 xi-mutation MUST-FAIL (fluxcurves) + ball_gate radius-bar "
  "must-fail (verdicts) + tampered glue/value caught (emit E3)"},
 "pfaffian_ball_reproduction": {"batteries": ["pfaffian"], "gate":
  "ads-5-81 W0 ball (dps 60+90) byte-exact vs the reference receipts "
  "(TERRIER_PFAFFIAN_BANK); fresh-prime gate 0 diffs; fresh-prime mutation "
  "control caught"},
 "coni_notconi_control": {"batteries": ["conifold"], "gate":
  "CP0 routing: DKMM PFFV control returns NOT-CONI == reference; DKMM tower "
  "log-branch-flip + branch-swap mutations caught; transport fail-closed "
  "without coni_op"},
 "ffp_fingerprint_regression": {"batteries": ["ffp"], "gate":
  "int64 front door refuses p >= 2^31; brute q-battery exact; HV AESZ34 "
  "p53 + BCM AESZ4/AESZ11 p61 pass-sets REGENERATED == shipped reference; "
  "ground-truth p19 replay (b+40 control pinned 14/15); checker + BCM "
  "control batteries == shipped CHECKER_GATE/BCM_GATE; stripped-fire mutation "
  "caught"},
 "sunit_planted_point_gate": {"batteries": ["sunit"], "gate":
  "planted rational + quadratic points refound under a reduced documented "
  "net with HAND-DERIVED S_min EXACT (shipped synthetic control card); "
  "singular + out-of-box must-fail controls; known-orbit tagging law; "
  "truncation statement machine-readable (bounded net, NOT a finiteness "
  "proof); zero blowup alarms"},
 "census_fold_acceptance": {"batteries": ["foldhnf"], "gate":
  "G4 anchor-set EXACT: B=1 bin 6612 + window total 49544; 256 corner "
  "merges witnessed; 5 tampered-word must-fails; order invariance"},
 "census_twostack_red": {"batteries": ["dedup_census"], "gate":
  "8/8: T1-T8 synthetic (partition == proof-level truth both engines; "
  "RED-RECONCILED-MISS + RED-WITNESS-FAIL drills)"},
 "census_distcert_off": {"batteries": ["distcert"], "gate":
  "37/37: exact balls, certified sqrt/ln, z/WP charts, Fincke-Pohst "
  "exhaustive OFF certificates, planted near-miss balls certified OFF"},
 "census_sampler_exact": {"batteries": ["controlsampler"], "gate":
  "19/19: dyadic-coin AD rejection exact, CP brackets exact, "
  "p_U only-widens, no-exclusion tubes at closed-form rates"},
 "census_smoke_e2e": {"batteries": ["census_smoke"], "gate":
  "synthetic 2x2 window end-to-end: fold(48->4, witnessed+S/T-invariant) "
  "-> two-stack dedup GREEN-AGREE (2 transported merge back, plant "
  "own-class) -> dist_cert ON/OFF -> common receipt validated"},
 "summation_reduced_u1s": {"batteries": ["summation"], "gate":
  "reduced-size u1s: ratio-law 200; termwise 12 ops x 4^6 all zero; "
  "S6-equivariance; GF==lattice EXACT M=3 rho+theta; K4 entropy law; "
  "b7 Cauchy floors 4/4 (full-size gates need TERRIER_BOXOPS_BANK)"},
 "slice_operator_card": {"batteries": ["cards", "slicecard"], "gate":
  "hv4-diag-L5 bank battery PASS (fit==printed JKK op, exponents, MUM, "
  "pfaffian payload) + receipt fields == reference; cards.py n_slice >= 1"},
 "lpkill_a1n_mustfail": {"batteries": ["lpkill"], "gate":
  "A1^6 control line BYTE-matches the reference (strict h=2 KILLED, relaxed "
  "FEASIBLE n0=36) + d4 genus replay JSON-equal vs the reference receipt "
  "(TERRIER_R18_DIR)"},
 "genus_control_replay": {"batteries": ["genusenum"], "gate":
  "Aut-free Kneser control genus det 92: CONTROL-PASS (6 classes, mass "
  "53/768 exact, qfisom bijection + per-class |Aut| vs the reference "
  "receipt (TERRIER_GENUS_CONTROL_RECEIPT))"},
 "dedup_orbit_controls": {"batteries": ["dedup"], "gate":
  "packaged KNOWN_ORBIT 9 rows; 5 must-fire (incl. the X118 pullback "
  "of T3) + 3 must-not-fire orbit controls; newform backstop 4 fire + 1 "
  "unknown must-not-fire; explicit-receipt law (no ledger write without "
  "a path); CLI rc 10/0"}}

man = {"generated": "terrier " + open(os.path.join(ROOT, "VERSION")).read().strip(),
       "root": ".", "note":
       "shipped-fixture manifest: selftest.py verifies every pin FIRST and "
       "refuses battery replays on any miss (NAMED, no silent skips). "
       "Batteries whose inputs are reference data not included in the "
       "package print a NAMED SKIP naming the env var to set.",
       "artifacts": {p: sha(os.path.join(ROOT, p)) for p in art},
       "batteries": bats, "wing_regressions": wreg}
man["n_artifacts"] = len(man["artifacts"])
out = os.path.join(ROOT, "regression_manifest.json")
json.dump(man, open(out, "w"), indent=1, sort_keys=True)
print(f"wrote {out}: {man['n_artifacts']} pins, {len(bats)} batteries")
