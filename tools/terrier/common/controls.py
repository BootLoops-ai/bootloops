#!/usr/bin/env python3
"""controls.py — terrier common chassis: blind controls, mutation harness,
planted-error verification.

CTRL-POS-V1 / CTRL-ID-V1 hash placement, blind control injection, the
VOID_CONTROL_FAIL rule, the planted-survivor drill, and the mutation-testing
law: no discrimination claim without a demonstrated failure under a
deliberate flip. stdlib only.
Determinism: NO RNG here — placement and check ids are hash-determined.
"""
import hashlib

POS_TAG = "CTRL-POS-V1"
ID_TAG = "CTRL-ID-V1"


def hpos(unit_id, n, tag=POS_TAG):
    """blind control position: H(tag:unit) mod (n+1) — deterministic, unguessable
    by the worker, recorded in the receipt."""
    return int(hashlib.sha256(("%s:%s" % (tag, unit_id)).encode())
               .hexdigest(), 16) % (n + 1)


def hid(*parts, tag=ID_TAG):
    """hash-determined anonymous check_id — worker cannot tell control from
    witness."""
    s = tag + ":" + ":".join(map(str, parts))
    return "X" + hashlib.sha256(s.encode()).hexdigest()[:8]


def blind_order(unit_id, controls, id_key="control_id"):
    """deterministic hash-blinded injection order (sorted by hid(unit, id))."""
    return sorted(controls, key=lambda c: hid(unit_id, c[id_key]))


def check_expected(got, expected):
    """subset match: every expected key present and equal in got."""
    return all(got.get(k) == v for k, v in expected.items())


def run_blind_controls(runner, unit_id, controls, id_key="control_id"):
    """Run every control through runner(req, check_id) under anonymous ids.
    Returns (results, missed). ANY miss => enclosing unit is VOID_CONTROL_FAIL
    (void rule): caller must NOT record it as COMPLETE."""
    results, missed = [], []
    for c in blind_order(unit_id, controls, id_key):
        cid = hid(unit_id, c[id_key])
        got = runner(c["req"], cid)
        ok = check_expected(got, c["expected"])
        results.append({"control_id": c[id_key], "check_id": cid, "pass": ok,
                        "kind": c.get("kind"), "got": got})
        if not ok:
            missed.append(c[id_key])
    return results, missed


def mutation_must_fail(gate, mutate, restore=None, name=""):
    """Mutation-testing law: apply the deliberate flip; the gate MUST fail —
    raise AssertionError or return falsy. restore() always runs. Returns
    {"name", "detected", "mode", "detail"}; caller asserts detected. A gate
    that still passes under mutation discriminates NOTHING (e.g. a xi-sign
    flip or a wrong-sign trap must be caught)."""
    mutate()
    try:
        try:
            r = gate()
        except AssertionError as e:
            return {"name": name, "detected": True, "mode": "raised",
                    "detail": str(e)[:200]}
        return {"name": name, "detected": not bool(r), "mode": "returned",
                "detail": repr(r)[:200]}
    finally:
        if restore:
            restore()


def mutation_battery(cases):
    """cases: iterable of (name, gate, mutate, restore). Returns (rows, all_ok).
    all_ok iff EVERY deliberate flip was detected."""
    rows = [mutation_must_fail(g, m, r, name=n) for n, g, m, r in cases]
    return rows, all(row["detected"] for row in rows)


def verify_planted(detector, clean, planted):
    """Planted-error drill: detector(clean) must
    PASS and detector(planted) must FIRE (raise or return falsy). Both legs
    required. Result is loudly labeled synthetic — never leave a fake alarm
    looking real."""
    steps = []

    def leg(step, ok, detail):
        steps.append({"step": step, "pass": bool(ok), "detail": str(detail)[:160]})

    try:
        c = detector(clean)
        leg("clean_passes", bool(c), repr(c))
    except AssertionError as e:
        leg("clean_passes", False, e)
    try:
        p = detector(planted)
        leg("planted_fires", not bool(p), repr(p))
    except AssertionError as e:
        leg("planted_fires", True, e)
    return {"REHEARSAL_SYNTHETIC": True, "pass": all(s["pass"] for s in steps),
            "steps": steps}
