#!/usr/bin/env python3
"""selftest_sweepchassis.py — battery for sweepchassis.py (toy example).
Gates: SC1 shards COMPLETE + blind controls pass; SC2 resume-skip; SC3 completed-shard
replay BYTE-IDENTICAL (pinned clock, file sha equal); SC4 planted control error ->
VOID_CONTROL_FAIL rc=8; SC5 planted integrity error -> VOID_INTEGRITY rc=8; SC6 halt
law end-to-end (evidence freeze, quarantine, dispatch refusal rc=9, clear_halt,
resume to COMPLETE, candidate archived); SC7 CTRL-ID-V1 / CTRL-POS-V1 hash laws
BYTE-COMPATIBLE with a reference shard receipt produced by an independent run of
the same hash laws (not included in the package: point TERRIER_SWEEP_SHARD_RECEIPT
at a shard receipt JSON with fields shard_id, controls.tierC[].{control_id,
check_id} and controls.injection_pos_hash; read-only).
Writes SELFTEST_SWEEPCHASSIS.json next to this file; exit 0 iff all gates pass.
Run: ulimit -v 32505856; nice -n 5 python3 selftest_sweepchassis.py
"""
import hashlib, json, os, shutil, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sweepchassis as sc

WORK = os.path.join(HERE, "selftest_chassis_work")
REFERENCE = os.environ.get("TERRIER_SWEEP_SHARD_RECEIPT")
if not REFERENCE:
    sys.exit("REFUSE: env var TERRIER_SWEEP_SHARD_RECEIPT is unset — point it "
             "at a reference shard receipt JSON (not included in the package; read-only)")
PINNED_CLOCK = lambda: 1752192000.0   # byte-identical replay law (SC3)


def evaluate(req):
    """Toy predicate-chain evaluator: exact integer square oracle."""
    if req.get("op") == "square":
        return {"y": req["x"] * req["x"], "verdict": "OK"}
    return {"verdict": "UNKNOWN"}


def adjudicate(evaluate, sid, i, e):
    """Toy chain: ACCEPT iff n % 3 != 0; planted flags exercise void/halt paths."""
    v = "ACCEPT" if e["n"] % 3 else "REJECT"
    out = {"verdict": v, "record": {"i": i, "n": e["n"], "verdict": v}}
    if e.get("bad"):
        out["integrity"] = ["planted integrity mismatch at i=%d" % i]
    if e.get("survivor"):
        out["survivor"] = True
    return out


CONTROLS = [{"control_id": "CP_SQ_%d" % k, "kind": "positive",
             "req": {"op": "square", "x": k}, "expected": {"y": k * k}}
            for k in (2, 5, 9)]
ENTRIES = [{"n": n} for n in range(1, 41)]
SSHA = hashlib.sha256(json.dumps(ENTRIES, sort_keys=True).encode()).hexdigest()
STREAMS = [("TOY", ENTRIES, SSHA)]


def fresh(clock=None):
    shutil.rmtree(WORK, ignore_errors=True)
    os.makedirs(WORK)
    kw = {"clock": clock} if clock else {}
    return sc.Chassis(WORK, [os.path.join(HERE, "sweepchassis.py")],
                      label="SELFTEST-TOY", scope_line="toy battery only", **kw)


def main():
    gates = []

    def gate(name, ok, detail=""):
        gates.append({"gate": name, "pass": bool(ok), "detail": str(detail)})
        print("  %-34s %s  %s" % (name, "PASS" if ok else "FAIL", detail))

    # SC1: two shards run COMPLETE with blind controls; SC2: resume-skip
    ch = fresh(clock=PINNED_CLOCK)
    with ch.lock():
        rc, done, skip = ch.run(evaluate, adjudicate, STREAMS, CONTROLS, batch_size=25)
    r0 = json.load(open(ch.receipt_path("TOYB000")))
    gate("SC1.run_complete", rc == 0 and done == 2 and skip == 0,
         "rc=%d done=%d skip=%d" % (rc, done, skip))
    gate("SC1.controls_blind_pass", r0["controls"]["pass"]
         and len(r0["controls"]["tier"]) == 3 and r0["status"] == "COMPLETE",
         "3 controls hash-placed, status=%s" % r0["status"])
    rc2, done2, skip2 = ch.run(evaluate, adjudicate, STREAMS, CONTROLS, batch_size=25)
    gate("SC2.resume_skip", rc2 == 0 and done2 == 0 and skip2 == 2,
         "rc=%d done=%d skip=%d" % (rc2, done2, skip2))

    # SC3: replay a completed shard -> byte-identical receipt file (pinned clock)
    sha_before = sc.sha_file(ch.receipt_path("TOYB000"))
    ident = ch.replay(evaluate, adjudicate, "TOYB000", "TOY", SSHA, ENTRIES[:25],
                      CONTROLS)
    sha_after = sc.sha_file(ch.receipt_path("TOYB000"))
    gate("SC3.replay_byte_identical", ident and sha_before == sha_after,
         "sha %s == %s" % (sha_before[:12], sha_after[:12]))

    # SC4: planted control error -> VOID_CONTROL_FAIL, rc=8
    ch = fresh()
    planted = CONTROLS + [{"control_id": "CP_PLANT", "kind": "planted-error",
                           "req": {"op": "square", "x": 3}, "expected": {"y": 10}}]
    rc, done, _ = ch.run(evaluate, adjudicate, STREAMS, planted, batch_size=25)
    st = json.load(open(ch.receipt_path("TOYB000")))["status"]
    gate("SC4.planted_error_rc8", rc == 8 and st == "VOID_CONTROL_FAIL",
         "rc=%d status=%s" % (rc, st))

    # SC5: planted integrity mismatch -> VOID_INTEGRITY, rc=8
    ch = fresh()
    ents = [dict(e) for e in ENTRIES]
    ents[7]["bad"] = True
    rc, _, _ = ch.run(evaluate, adjudicate, [("TOY", ents, SSHA)], CONTROLS,
                      batch_size=25)
    st = json.load(open(ch.receipt_path("TOYB000")))["status"]
    gate("SC5.integrity_void_rc8", rc == 8 and st == "VOID_INTEGRITY",
         "rc=%d status=%s" % (rc, st))

    # SC6: halt law end-to-end (synthetic survivor drill, loudly rehearsal)
    ch = fresh()
    ents = [dict(e) for e in ENTRIES]
    ents[3]["survivor"] = True
    rc, _, _ = ch.run(evaluate, adjudicate, [("TOY", ents, SSHA)], CONTROLS,
                      batch_size=25, rehearsal=True)
    qp = os.path.join(WORK, "shards", "quarantine", "receipt_TOYB000.json")
    allowed, why = ch.dispatch_allowed()
    rc_ref, _, _ = ch.run(evaluate, adjudicate, STREAMS, CONTROLS, batch_size=25)
    att = os.path.join(WORK, "attention.jsonl")
    qq = json.load(open(qp)) if os.path.exists(qp) else {}
    gate("SC6.halt_fired_quarantined",
         rc == 9 and os.path.exists(ch.halt_file()) and os.path.exists(qp)
         and qq.get("status") == "QUARANTINED_SURVIVOR_REHEARSAL"
         and qq.get("REHEARSAL_SYNTHETIC") is True
         and not ch.receipt_ok("TOYB000", SSHA)
         and not allowed and rc_ref == 9 and os.path.exists(att),
         "rc=%d refused_rc=%d %s" % (rc, rc_ref, qq.get("status")))
    ch.clear_halt()
    ok2, _ = ch.dispatch_allowed()
    rc, done, _ = ch.run(evaluate, adjudicate, STREAMS, CONTROLS, batch_size=25)
    arch = os.listdir(os.path.join(WORK, "out", "adjudicated"))
    gate("SC6.clear_resume_complete", ok2 and rc == 0 and done == 2
         and len(arch) == 1 and arch[0].startswith("CANDIDATE_REHEARSAL_"),
         "re-ran %d shards; archived %d rehearsal candidate file(s)" % (done, len(arch)))

    # SC7: hash laws byte-compatible with the REFERENCE shard receipt (read-only)
    br = json.load(open(REFERENCE))
    sid, tc = br["shard_id"], br["controls"]["tierC"]
    ok_id = all(sc.hid(sid, c["control_id"]) == c["check_id"] for c in tc)
    order = [c["control_id"] for c in tc]
    ok_ord = order == sorted(order, key=lambda c: sc.hid(sid, c))
    ok_pos = sc.hpos(sid, len(tc)) == br["controls"]["injection_pos_hash"]
    gate("SC7.reference_checkids_byte_match", ok_id,
         "%d/%d CTRL-ID-V1 check_ids reproduced (shard %s)" % (
             sum(sc.hid(sid, c["control_id"]) == c["check_id"] for c in tc), len(tc), sid))
    gate("SC7.reference_injection_order", ok_ord, "hash-sorted order == reference order")
    gate("SC7.reference_pos_hash", ok_pos, "hpos(%s,%d)==%d" % (
        sid, len(tc), br["controls"]["injection_pos_hash"]))

    allpass = all(g["pass"] for g in gates)
    sc.atomic_write(os.path.join(HERE, "SELFTEST_SWEEPCHASSIS.json"),
                    {"schema": "sweepchassis-selftest-v1", "pass": allpass,
                     "gates": gates, "reference_receipt": os.path.basename(REFERENCE),
                     "code_sha256": sc.sha_file(os.path.join(HERE, "sweepchassis.py"))})
    shutil.rmtree(WORK, ignore_errors=True)
    print("selftest_sweepchassis: %s (%d/%d gates)"
          % ("PASS" if allpass else "FAIL", sum(g["pass"] for g in gates), len(gates)))
    return 0 if allpass else 1


if __name__ == "__main__":
    sys.exit(main())
