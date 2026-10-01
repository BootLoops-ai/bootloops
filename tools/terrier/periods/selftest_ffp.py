#!/usr/bin/env python3
"""selftest_ffp.py — battery for periods/ffp.py (finite-field fingerprint
engines; engine code in pipeline/ffp/).
F0 int64 front door: p >= 2^31 REFUSED; explicit big-prime flag passes
F1 brute-validation battery at tiny q (torus 7/11/13 deg1 + 5^2; twist12)
F2 HV AESZ34 pass-set REGENERATED at p=53 == stored (pilot2 law asserts in)
F3 BCM pass-sets REGENERATED at p=61 (AESZ4 h33 + AESZ11 h34) == stored
   on the exact fingerprint content (fires/alphas/counts/control);
   float_resid is machine-dependent in its low bits (libm/SIMD), so it is
   gated against the certified < 0.05 envelope, never compared byte-exact;
   unit-root congruence asserted at every z inside run()
F4 ground-truth p=19 replay (CDEvS table): Weil all-pass, b+40 control
   pinned at 14/15, fac-row splits recovered
F5 unit-root congruence gates h33+h34 at p=61 (importable-gate surface)
F6 checker regression battery == stored CHECKER_GATE.jsonl (per-case)
F7 BCM control battery PASS == stored BCM_GATE.jsonl verdict
F8 MUTATION control: control residue stripped from a pass-set -> caught
Packaged data: pipeline/ffp/{scale/passsets,banks}. Resource caps: ulimit
-v 32505856, nice >= 5."""
import copy, json, os, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ffp

FAIL = []


def gate(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name} {detail}", flush=True)
    if not ok:
        FAIL.append(name)


# Stripped from the byte-exact regen == stored comparisons: wall-clock, and
# the float residual certificate. The residual is a rounding-error bound
# (certified < 0.05 by an assert inside the engine); its trailing bits come
# from libm/SIMD dispatch and differ across machines, so it is gated against
# the envelope (run_F3), not byte-compared. Everything else in a pass-set is
# exact integer/string content and stays byte-compared.
VOLATILE = ("compute_seconds", "float_resid")


def nocs(d):
    return {k: v for k, v in d.items() if k not in VOLATILE}


def run_F0():
    try:
        ffp.hv_torus_counts(2147483659, 1)
        gate("F0a big prime refused", False, "no refusal raised")
    except ffp.BigPrimeRefusal as e:
        gate("F0a big prime refused", "REFUSED" in str(e)
             and "int64" in str(e), "p=2147483659 >= 2^31")
    ok = ffp.require_int64_prime(2147483659, allow_big_prime=True) \
        == 2147483659 and ffp.require_int64_prime(53) == 53
    gate("F0b flag + small prime pass front door", ok)


def run_F1():
    t0 = time.time()
    g = ffp.gate_brute_small_q()
    gate("F1 brute battery tiny q", g["pass"],
         f"{sorted(k for k in g if k != 'pass')} ({time.time()-t0:.0f}s)")


def run_F2():
    t0 = time.time()
    r = ffp.hv_passset_AESZ34(53)
    bank = json.load(open(os.path.join(ffp.PASSSET_DIR,
                                       "AESZ34.p53.json")))
    gate("F2 HV p53 pass-set regen == stored",
         json.dumps(nocs(r), sort_keys=True)
         == json.dumps(nocs(bank), sort_keys=True),
         f"fires {r.get('fire_count')} ({time.time()-t0:.0f}s)")


def run_F3():
    for tag, (op, hn, ctrl) in (("a", ("AESZ4", "h33", (-1, 5832))),
                                ("b", ("AESZ11", "h34", (-1, 432)))):
        t0 = time.time()
        r = ffp.bcm_passset(op, hn, ctrl, 61)
        bank = json.load(open(os.path.join(ffp.PASSSET_DIR,
                                           f"{op}.p61.json")))
        gate(f"F3{tag} BCM {op} p61 regen == stored",
             json.dumps(nocs(r), sort_keys=True)
             == json.dumps(nocs(bank), sort_keys=True)
             and r["float_resid"] < 0.05,
             f"ctrl fired {r['control']['fired']} "
             f"resid {r['float_resid']:.1e} (envelope 5e-2) "
             f"({time.time()-t0:.0f}s)")


def run_F4():
    g = ffp.gate_groundtruth_p19()
    gate("F4 ground-truth p19 replay", g["pass"]
         and g["n_weil_pass"] == g["n_smooth"] == 15
         and g["n_perturbed_pass"] == 14
         and not g["fac_rows_missed"],
         f"weil {g['n_weil_pass']}/15, b+40 control {g['n_perturbed_pass']}"
         " (pinned 14)")


def run_F5():
    for tag, hn in (("a", "h33"), ("b", "h34")):
        g = ffp.gate_unit_root_congruence(hn, 61)
        gate(f"F5{tag} unit-root congruence {hn} p61", g["pass"],
             f"{g['n_checked']} z checked, resid {g['float_resid']:.1e}")


def run_F6():
    t0 = time.time()
    b = ffp.checker_gate_battery()
    bank = [json.loads(l) for l in
            open(os.path.join(ffp.BANK_DIR, "CHECKER_GATE.jsonl"))]
    hdr, brow = bank[0], {r["name"]: r for r in bank[1:]}
    match = all(brow[r["name"]]["n_in"] == r["n_in"]
                and brow[r["name"]]["n_checkable"] == r["n_checkable"]
                and brow[r["name"]]["gate"] == r["gate"]
                for r in b["rows"])
    gate("F6 checker battery == stored CHECKER_GATE",
         b["verdict"] == hdr["verdict"] == "PASS" and match
         and b["n_primes"] == hdr["n_primes"],
         f"{len(b['rows'])} cases @ {b['n_primes']} primes "
         f"({time.time()-t0:.0f}s)")


def run_F7():
    t0 = time.time()
    b = ffp.bcm_gate_battery()
    hdr = json.loads(open(os.path.join(
        ffp.BANK_DIR, "BCM_GATE.jsonl")).readline())
    gate("F7 BCM control battery == stored BCM_GATE",
         b["verdict"] == hdr["verdict"] == "PASS",
         f"{len(b['rows'])} rows ({time.time()-t0:.0f}s)")


def run_F8():
    ps = ffp.load_passsets("AESZ4")
    mu = copy.deepcopy(ps)
    czr = -1 * pow(5832, 61 - 2, 61) % 61
    mu[61]["fires"].discard(czr)
    from fractions import Fraction
    r = ffp.check(dict(kind="rational", z=Fraction(-1, 5832)), mu)
    gate("F8 mutation control: stripped fire CAUGHT", not r["and_pass"],
         f"control residue {czr} removed at p=61 -> AND fails (must-fail)")


if __name__ == "__main__":
    T0 = time.time()
    run_F0()
    run_F1()
    run_F2()
    run_F3()
    run_F4()
    run_F5()
    run_F6()
    run_F7()
    run_F8()
    n = 12
    print(f"[selftest_ffp] {n - len(FAIL)}/{n} passed "
          f"({time.time()-T0:.0f}s wall)"
          + ("" if not FAIL else f"  FAILED: {FAIL}"))
    sys.exit(1 if FAIL else 0)
