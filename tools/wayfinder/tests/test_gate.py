#!/usr/bin/env python3
"""
Tests for tools/wayfinder/gate.py.

Synthetic vendored strings at 60 significant digits, all derived in-test
from mpmath closed forms (sqrt(2), pi*gammaE, zeta(3)) — no fabricated
digits.  Three mandated cases:

G1  exact match            -> PASS (matched >= 59, pair stable)
G2  wrong in last 5 digits -> FAIL (matched ~55 < 59; pair stable — the
                              vendored-string gate does the killing)
G3  two-precision unstable -> FAIL (pair agreement ~40 < 60 even though
                              one member matches the vendored string —
                              noise, not a result)

Plus: per-key report fields, overall verdict on the mixed dict, missing /
extra key handling, significant-digit counting edge cases, zero string,
and write_gate_report JSON round-trip.  Global mp.dps asserted untouched.

Degenerate-vendored negative controls (F2/F3 fixes):
Z1  nonzero computed vs vendored '0'      -> FAIL (was a vacuous PASS)
Z2  zero computed vs vendored '0', nonzero sibling scale -> PASS
Z3  zero-vendored with NO scale reference -> FAIL (un-gateable)
D1  1-digit nonzero vendored ('3')        -> FAIL even on agreement
P1  dps_pair (100,100)                    -> ValueError (one-precision gate)
P2  byte-identical pair strings           -> flagged in per_key

FEED-table controls (feed_gate_table — per-literal FEED gate
convention; pure formatting + charter-bar thresholds):
F1  good literal (stored 70, matched 69.8, pair 140)      -> PASS
F2  matched 55.2 < stored-1                               -> must-FAIL
F3  pair 60 < stored                                      -> must-FAIL
F4  stored 20 < charter bar 30                            -> must-FAIL
F5  unmeasured (matched None)                             -> must-FAIL
F6  all-good list -> overall PASS; empty list -> FAIL; malformed entry
    raises; table string carries header + verdict line; composes with
    gate_strings per_key output.
"""
import json
import os
import sys
import tempfile

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
sys.path.insert(0, PKG)

import gate  # noqa: E402

FAILS = []


def check(name, cond, msg=""):
    tag = "PASS" if cond else "FAIL"
    print("  [%s] %s%s" % (tag, name, (": " + msg) if msg else ""))
    if not cond:
        FAILS.append(name)


DPS0 = mp.mp.dps


def nstr_at(x, n):
    """n-significant-digit decimal string of an mp value (built at high
    precision so the string itself is correctly rounded)."""
    return mp.nstr(x, n, strip_zeros=False)


def corrupt_last(s, n=5):
    """Replace the last n digit characters of a decimal string (digit ->
    (digit+5) mod 10) — 'wrong in the last n digits'."""
    out = list(s)
    done = 0
    for i in range(len(out) - 1, -1, -1):
        if out[i].isdigit():
            out[i] = str((int(out[i]) + 5) % 10)
            done += 1
            if done == n:
                break
    return "".join(out)


def main():
    dps_pair = (80, 100)
    with mp.workdps(120):
        v_sqrt2 = mp.sqrt(2)
        v_pig = mp.pi * mp.euler
        v_z3 = mp.zeta(3)
        vendored = {
            "sqrt2": nstr_at(v_sqrt2, 60),
            "pig": nstr_at(v_pig, 60),
            "z3": nstr_at(v_z3, 60),
        }
        # G1 exact: both precisions correct
        c_sqrt2 = (nstr_at(v_sqrt2, 80), nstr_at(v_sqrt2, 100))
        # G2 wrong in last 5 digits of the 60d claim, but PAIR-STABLE
        bad = corrupt_last(vendored["pig"], 5)
        c_pig = (bad, bad)
        # G3 two-precision unstable: hi member drifts at ~40 digits
        drift = v_z3 * (1 + mp.mpf(10) ** (-40))
        c_z3 = (nstr_at(v_z3, 80), nstr_at(drift, 100))

    computed = {"sqrt2": c_sqrt2, "pig": c_pig, "z3": c_z3}
    rep = gate.gate_strings(computed, vendored, dps_pair)

    pk = rep["per_key"]
    check("G1-exact-PASS", pk["sqrt2"]["PASS"] is True,
          "matched=%s pair=%s" % (pk["sqrt2"]["matched_digits"],
                                  pk["sqrt2"]["pair_agree_digits"]))
    check("G1-stored-60", pk["sqrt2"]["stored_digits"] == 60)
    check("G1-matched-full", pk["sqrt2"]["matched_digits"] >= 59)

    check("G2-last5-FAIL", pk["pig"]["PASS"] is False,
          "matched=%s" % pk["pig"]["matched_digits"])
    check("G2-matched-window", 50 <= pk["pig"]["matched_digits"] <= 58,
          "matched=%s (expect ~55)" % pk["pig"]["matched_digits"])
    check("G2-pair-stable", pk["pig"]["pair_agree_digits"] >= 60,
          "pair=%s (identical strings)" % pk["pig"]["pair_agree_digits"])
    check("G2-reason", "matched" in (pk["pig"]["reason"] or ""))

    check("G3-unstable-FAIL", pk["z3"]["PASS"] is False)
    check("G3-pair-window", 35 <= pk["z3"]["pair_agree_digits"] <= 45,
          "pair=%s (expect ~40)" % pk["z3"]["pair_agree_digits"])
    check("G3-reason", "unstable" in (pk["z3"]["reason"] or ""))

    check("overall-FAIL", rep["overall"] == "FAIL" and rep["n_pass"] == 1)

    # all-pass dict -> overall PASS
    rep2 = gate.gate_strings({"sqrt2": c_sqrt2},
                             {"sqrt2": vendored["sqrt2"]}, dps_pair)
    check("overall-PASS", rep2["overall"] == "PASS")

    # missing + extra keys
    rep3 = gate.gate_strings({"sqrt2": c_sqrt2, "orphan": c_sqrt2},
                             {"sqrt2": vendored["sqrt2"],
                              "gone": vendored["z3"]}, dps_pair)
    check("missing-FAILs", rep3["overall"] == "FAIL"
          and rep3["missing_keys"] == ["gone"]
          and rep3["per_key"]["gone"]["PASS"] is False)
    check("extra-noted", rep3["extra_keys"] == ["orphan"])

    # significant-digit counting edges
    sd = gate._sig_digits
    check("sig-1.0", sd("1.0") == 2)
    check("sig-lead0", sd("0.000123") == 3)
    check("sig-signed-exp", sd("-3.1415e-07") == 5)
    check("sig-zero", sd("0.00") == 1)
    check("sig-60", sd(vendored["sqrt2"]) == 60)

    # ---- degenerate-vendored controls --------------------------------
    # Z1 nonzero computed vs vendored '0' must FAIL (was a vacuous PASS:
    # matched clamped to 0 >= n_stored-1 = 0). Sibling sqrt2 sets the scale.
    repz1 = gate.gate_strings({"nil": ("5.0", "5.0"), "sqrt2": c_sqrt2},
                              {"nil": "0", "sqrt2": vendored["sqrt2"]},
                              dps_pair)
    check("Z1-nonzero-vs-zero-FAIL",
          repz1["per_key"]["nil"]["PASS"] is False
          and repz1["overall"] == "FAIL",
          "reason=%s" % repz1["per_key"]["nil"]["reason"])
    # Z2 genuinely-zero computed vs vendored '0' with a nonzero sibling
    repz2 = gate.gate_strings({"nil": ("0.0", "0.0"), "sqrt2": c_sqrt2},
                              {"nil": "0.0", "sqrt2": vendored["sqrt2"]},
                              dps_pair)
    check("Z2-zero-vs-zero-PASS", repz2["overall"] == "PASS",
          "zero_bar=%s" % repz2["per_key"]["nil"].get("zero_bar"))
    # Z3 zero-vendored with NO scale reference (no nonzero sibling, no
    # explicit zero_scale) is un-gateable -> FAIL
    repz3 = gate.gate_strings({"nil": ("0.0", "0.0")}, {"nil": "0.0"},
                              dps_pair)
    check("Z3-no-scale-FAIL", repz3["overall"] == "FAIL"
          and "scale" in (repz3["per_key"]["nil"]["reason"] or ""))
    # Z4 explicit zero_scale/zero_digits make the lone zero gateable
    repz4 = gate.gate_strings({"nil": ("1e-70", "0.0")}, {"nil": "0"},
                              dps_pair, zero_scale="1.0", zero_digits=60)
    check("Z4-explicit-scale-PASS", repz4["overall"] == "PASS")
    repz5 = gate.gate_strings({"nil": ("1e-30", "1e-30")}, {"nil": "0"},
                              dps_pair, zero_scale="1.0", zero_digits=60)
    check("Z5-above-bar-FAIL", repz5["overall"] == "FAIL")
    # D1 nonzero 1-digit vendored string is un-gateable -> FAIL even when
    # the computed pair "agrees" (('7','7') vs '3' passed before the fix)
    repd1 = gate.gate_strings({"one": ("7", "7"), "sqrt2": c_sqrt2},
                              {"one": "3", "sqrt2": vendored["sqrt2"]},
                              dps_pair)
    check("D1-1digit-FAIL", repd1["per_key"]["one"]["PASS"] is False
          and "1 significant digit" in (repd1["per_key"]["one"]["reason"] or ""))
    repd2 = gate.gate_strings({"one": ("3", "3")}, {"one": "3"}, dps_pair)
    check("D1b-1digit-agree-FAIL", repd2["overall"] == "FAIL")
    # P1 equal precisions must raise (one-precision gate in disguise)
    raised = False
    try:
        gate.gate_strings({"sqrt2": c_sqrt2},
                          {"sqrt2": vendored["sqrt2"]}, (100, 100))
    except ValueError:
        raised = True
    check("P1-equal-dps-raises", raised)
    raised = False
    try:
        gate.gate_strings({"sqrt2": c_sqrt2},
                          {"sqrt2": vendored["sqrt2"]}, (100, 80))
    except ValueError:
        raised = True
    check("P1b-reversed-dps-raises", raised)
    # P2 byte-identical pair members flagged (G2 used (bad, bad))
    check("P2-identical-flagged",
          rep["per_key"]["pig"]["pair_identical_strings"] is True
          and rep["per_key"]["sqrt2"]["pair_identical_strings"] is False)

    # malformed computed entry -> per-key FAIL, not crash
    rep5 = gate.gate_strings({"sqrt2": nstr_at(v_sqrt2, 60)},
                             {"sqrt2": vendored["sqrt2"]}, dps_pair)
    check("not-a-pair-FAIL", rep5["per_key"]["sqrt2"]["PASS"] is False
          and "2-tuple" in rep5["per_key"]["sqrt2"]["reason"])

    # write_gate_report round-trip (tempfile only; gate_strings wrote nothing)
    fd, path = tempfile.mkstemp(suffix=".json", prefix="wayfinder_gate_")
    os.close(fd)
    try:
        gate.write_gate_report(rep, path)
        back = json.load(open(path))
        check("report-roundtrip", back["overall"] == "FAIL"
              and back["per_key"]["sqrt2"]["PASS"] is True
              and back["dps_pair"] == [80, 100])
    finally:
        os.unlink(path)

    # ---- FEED-table controls (feed_gate_table) --------------------------
    entries = [
        {"literal": "row8_k0_re", "stored_digits": 70,
         "matched_digits": 69.8, "pair_agree_digits": 140.0},           # F1
        {"literal": "bad_matched", "stored_digits": 70,
         "matched_digits": 55.2, "pair_agree_digits": 140.0},           # F2
        {"literal": "bad_pair", "stored_digits": 70,
         "matched_digits": 70.0, "pair_agree_digits": 60.0},            # F3
        {"literal": "below_bar", "stored_digits": 20,
         "matched_digits": 20.0, "pair_agree_digits": 25.0},            # F4
        {"literal": "unmeasured", "stored_digits": 70,
         "matched_digits": None, "pair_agree_digits": 140.0},           # F5
    ]
    frep = gate.feed_gate_table(entries)
    pl = {r["literal"]: r for r in frep["per_literal"]}
    check("F1-good-PASS", pl["row8_k0_re"]["PASS"] is True)
    check("F2-matched-FAIL", pl["bad_matched"]["PASS"] is False
          and "matched" in pl["bad_matched"]["reason"])
    check("F3-pair-FAIL", pl["bad_pair"]["PASS"] is False
          and "two-precision" in pl["bad_pair"]["reason"])
    check("F4-charter-bar-FAIL", pl["below_bar"]["PASS"] is False
          and "charter bar" in pl["below_bar"]["reason"])
    check("F5-unmeasured-FAIL", pl["unmeasured"]["PASS"] is False
          and "unmeasured" in pl["unmeasured"]["reason"])
    check("F-overall", frep["overall"] == "FAIL" and frep["n_pass"] == 1
          and frep["n_literals"] == 5)
    check("F-table-format", "FEED GATE: FAIL (1/5 literals)" in frep["table"]
          and "literal" in frep["table"].splitlines()[1]
          and frep["table"].count("FAIL") >= 5)
    # boundary case: matched == stored-1 and pair == stored PASS exactly
    fedge = gate.feed_gate_table([{"literal": "edge", "stored_digits": 30,
                                   "matched_digits": 29.0,
                                   "pair_agree_digits": 30.0}])
    check("F-edge-PASS", fedge["overall"] == "PASS")
    # all-good -> PASS; empty -> FAIL (an empty table gates nothing)
    fall = gate.feed_gate_table([entries[0]])
    check("F6-all-good-PASS", fall["overall"] == "PASS"
          and "FEED GATE: PASS (1/1 literals)" in fall["table"])
    check("F6-empty-FAIL", gate.feed_gate_table([])["overall"] == "FAIL")
    raised = False
    try:
        gate.feed_gate_table([{"literal": "x"}])
    except ValueError:
        raised = True
    check("F6-malformed-raises", raised)
    # composes with gate_strings per_key output (the intended wiring)
    fcomp = gate.feed_gate_table(
        [{"literal": k,
          "stored_digits": rep2["per_key"][k]["stored_digits"],
          "matched_digits": rep2["per_key"][k]["matched_digits"],
          "pair_agree_digits": rep2["per_key"][k]["pair_agree_digits"]}
         for k in rep2["per_key"]])
    check("F6-compose-gate_strings", fcomp["overall"] == "PASS",
          fcomp["table"].splitlines()[-1])


def test_gate_suite():
    """pytest entry point (the file is otherwise standalone-style): run the
    full check list; the FAILS accumulator must stay empty."""
    del FAILS[:]
    main()
    assert not FAILS, "gate self-checks failed: %s" % ",".join(FAILS)
    assert mp.mp.dps == DPS0, "global mp.dps mutated (import-dps footgun)"


if __name__ == "__main__":
    print("== test_gate ==")
    main()
    check("global-dps-untouched", mp.mp.dps == DPS0,
          "dps %d -> %d" % (DPS0, mp.mp.dps))
    print("== %s ==" % ("ALL PASS" if not FAILS else "FAILURES: " + ",".join(FAILS)))
    sys.exit(0 if not FAILS else 1)
