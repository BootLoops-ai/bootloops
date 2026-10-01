#!/usr/bin/env python3
"""selftest_sunit.py — battery for periods/sunit.py (refined
S-integrality generator; pipeline/sunit generator code).

Data: the shipped atlas holds ONE synthetic control card
(pipeline/sunit/atlas/TOY7.json) whose singular divisor {z=1/7, z=-1}
makes every acceptance value derivable by hand — the expected S_min
sets below are computed in the card's own note, not copied from a run.

S0 census net defaults pinned (the frozen production caps)
S1 truncation statement is machine-readable + honest (bounded net,
   NOT a finiteness proof; blowup alarm threshold)
S2 planted-point gate: under a reduced documented net the generator
   REFINDS the planted rational point z=-1/7 with S_min EXACTLY
   {2,3,7} and the planted quadratic point 3+2*sqrt(2) (key 1,-6,1)
   with S_min EXACTLY {2} + unit flag — hand-derived, not copied from a run
S3 known-orbit dedup tagging: a planted known-orbit key is tagged
   KNOWN-ORBIT-REDISCOVERY; every other row stays untagged
S4 must-fail controls: out-of-box support REJECTED; singular point
   itself REJECTED (acceptance is exact, not a sieve)
S5 quadratic machinery: fund_unit(2)=(1,1,1) -> u^2 minpoly key
   1,-6,1; fund_unit(17)=(4,1,1) -> u^2 key 1,-66,1 (pure algebra)
S6 stats receipt: machine-readable, zero blowup alarms
Resource caps: ulimit -v 32505856, nice >= 5. Writes only to
pipeline/sunit/out_selftest/ (volatile, regenerated every run)."""
import json, os, shutil, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sunit

OUT = os.path.join(sunit.SUNIT_DIR, "out_selftest")
OP = "TOY7"
FAIL = []


def gate(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name} {detail}", flush=True)
    if not ok:
        FAIL.append(name)


def run_S0_S1():
    ts = sunit.truncation_statement()
    net = ts["net"]
    gate("S0 census net defaults pinned",
         net["A_MAX"] == 12 and net["H_RAT"] == 40.0
         and net["RAT_SUPP_MAX"] == 4 and net["H_POOL"] == 20.0
         and net["POOL_CAP"] == 600 and net["H_QUAD"] == 40.0
         and len(net["SQFREE_D"]) == 60
         and ts["box"] == [2, 3, 5, 7, 11, 13])
    gate("S1 truncation statement machine-readable + honest",
         "NOT A FINITENESS PROOF" in ts["completeness"]
         and ts["blowup_alarm_threshold"] == 10 ** 5
         and ts["law"].startswith("refined-conductor-boxS"))


def run_S4_S5():
    cards, _al = sunit.load_cards()
    fac = sunit.sing_factors(cards[OP])
    gate("S4a out-of-box support REJECTED (must-fail)",
         sunit.smin_rational(fac, 1, 17) is None
         and sunit.smin_rational(fac, 1, 2) == {2, 3, 5},
         "z=1/17 None; z=1/2 -> {2,3,5} in-box sanity")
    gate("S4b singular point itself REJECTED (must-fail)",
         sunit.smin_rational(fac, 1, 7) is None
         and sunit.smin_rational(fac, 1, 1, sgn=-1) is None,
         "z=1/7 and z=-1 are the card's singular values")
    import gen_sunit as GS
    u2a = GS.q_mul(*(sunit.fund_unit(2),) * 2, 2)
    u17 = sunit.fund_unit(17)
    u2b = GS.q_mul(u17, u17, 17)
    gate("S5a fundamental units + u^2 minpoly keys",
         sunit.fund_unit(2) == (1, 1, 1) and u2a == (3, 2, 1)
         and GS.q_minpoly_key(u2a, 2) == "1,-6,1"
         and u17 == (4, 1, 1) and u2b == (33, 8, 1)
         and GS.q_minpoly_key(u2b, 17) == "1,-66,1",
         "D=2: (1+sqrt2)^2 = 3+2*sqrt2; D=17: (4+sqrt17)^2 = 33+8*sqrt17")
    r = sunit.smin_quadratic((3, 2, 1), 2, fac)
    gate("S5b quadratic acceptance: S_min {2}, unit",
         r is not None and sorted(r[0]) == [2] and r[2] is True
         and r[3] == "1,-6,1",
         f"S_min={sorted(r[0]) if r else None} (hand-derived: both "
         "singular-factor norms = 8)")


def run_S2_S3_S6():
    if os.path.isdir(OUT):
        shutil.rmtree(OUT)
    sunit.configure_net(RAT_SUPP_MAX=2, H_POOL=6.0, POOL_CAP=100,
                        SQFREE_D=[2])
    t0 = time.time()
    # planted known-orbit key: exercises the tagging law on the toy op
    st, _cache = sunit.generate_op(OP, OUT,
                                   cache={"known": {OP: {"rat:-1/7"}}})
    rows = [json.loads(l) for l in open(os.path.join(OUT, f"{OP}.jsonl"))]
    byk = {r["match_key"]: r for r in rows}
    rat = byk.get("rat:-1/7")
    quad = byk.get("alg:1,-6,1")
    gate("S2 planted points refound, S_min EXACT (hand-derived)",
         rat is not None and rat["S_min"] == [2, 3, 7]
         and quad is not None and quad["S_min"] == [2]
         and quad["unit"] is True,
         f"reduced net, {len(rows)} rows ({time.time()-t0:.0f}s)")
    tagged = {r["match_key"] for r in rows
              if r.get("dedup") == "KNOWN-ORBIT-REDISCOVERY"}
    gate("S3 known-orbit tagging (planted key only)",
         tagged == {"rat:-1/7"}, sorted(tagged))
    gate("S6 stats receipt machine-readable, zero blowup alarms",
         st["rows"] == len(rows) and st["known_orbit_strips"] == 1
         and st["blowup_alarm"] is False
         and st["factors"] == [[-1, 7], [1, 1]],
         f"{st['rows']} rows, factors {st['factors']}")


if __name__ == "__main__":
    T0 = time.time()
    run_S0_S1()
    run_S4_S5()          # acceptance controls BEFORE the net is reduced
    run_S2_S3_S6()
    n = 9
    print(f"[selftest_sunit] {n - len(FAIL)}/{n} passed "
          f"({time.time()-T0:.0f}s wall)"
          + ("" if not FAIL else f"  FAILED: {FAIL}"))
    sys.exit(1 if FAIL else 0)
