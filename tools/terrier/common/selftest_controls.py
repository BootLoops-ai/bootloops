#!/usr/bin/env python3
"""battery: controls.py — hash placement, blind injection, mutation, planted."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import controls as C

n = [0]
def ok(name, cond):
    n[0] += 1
    print("[%s] %s" % ("PASS" if cond else "FAIL", name))
    assert cond, name

# hash placement: deterministic, in range, matches the pinned V1 definition
import hashlib
ref = int(hashlib.sha256(b"CTRL-POS-V1:S000001").hexdigest(), 16) % 8
ok("hpos == pinned CTRL-POS-V1", C.hpos("S000001", 7) == ref)
ok("hpos range", all(0 <= C.hpos("S%06d" % i, 5) <= 5 for i in range(50)))
ok("hid format", C.hid("S0", "CC_X").startswith("X") and len(C.hid("S0", "CC_X")) == 9)
ok("hid deterministic", C.hid("S0", "CC_X") == C.hid("S0", "CC_X"))
ctrls = [{"control_id": "CC_%d" % i, "kind": "positive",
          "req": {"n": i}, "expected": {"verdict": "ACCEPT"}} for i in range(6)]
o1 = [c["control_id"] for c in C.blind_order("S0", ctrls)]
o2 = [c["control_id"] for c in C.blind_order("S1", ctrls)]
ok("blind order deterministic", o1 == [c["control_id"] for c in C.blind_order("S0", ctrls)])
ok("blind order unit-dependent", o1 != o2 and sorted(o1) == sorted(o2))
# blind run: worker sees only anonymous check ids; one planted negative miss
seen = []
def runner(req, check_id):
    seen.append(check_id)
    return {"verdict": "ACCEPT" if req["n"] != 3 else "REJECT"}
res, missed = C.run_blind_controls(runner, "S0", ctrls)
ok("void rule: miss reported", missed == ["CC_3"])
ok("results carry check ids", all(r["check_id"].startswith("X") for r in res))
ok("no control_id leaked to worker", all(s.startswith("X") for s in seen))
# mutation harness: flipped sign MUST be detected; blind gate must NOT pass
state = {"xi": -270}
gate = lambda: (state["xi"] == -270) or (_ for _ in ()).throw(AssertionError("xi flip"))
m = C.mutation_must_fail(gate, lambda: state.update(xi=270),
                         lambda: state.update(xi=-270), name="xi-flip")
ok("mutation detected (raise)", m["detected"] and state["xi"] == -270)
m2 = C.mutation_must_fail(lambda: True, lambda: state.update(xi=270),
                          lambda: state.update(xi=-270), name="blind-gate")
ok("mutation-blind gate exposed", not m2["detected"])
rows, allok = C.mutation_battery([
    ("xi", gate, lambda: state.update(xi=270), lambda: state.update(xi=-270)),
    ("blind", lambda: True, lambda: None, None)])
ok("battery all_ok only if all detected", not allok and rows[0]["detected"])
# planted-error drill: detector passes clean, fires on planted survivor
det = lambda receipt: not receipt.get("survivors")
v = C.verify_planted(det, {"survivors": []}, {"survivors": [{"N_flux": 12}]})
ok("planted drill passes", v["pass"] and v["REHEARSAL_SYNTHETIC"] is True)
v2 = C.verify_planted(lambda r: True, {}, {"survivors": [1]})
ok("blind detector fails drill", not v2["pass"])
print("selftest_controls: %d/%d green" % (n[0], n[0]))
