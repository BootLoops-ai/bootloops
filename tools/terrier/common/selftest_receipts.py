#!/usr/bin/env python3
"""battery: receipts.py — atomic write, sha-stream, code-sha, jsonl, resume law."""
import hashlib, json, os, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import receipts as R

n = [0]
def ok(name, cond):
    n[0] += 1
    print("[%s] %s" % ("PASS" if cond else "FAIL", name))
    assert cond, name

with tempfile.TemporaryDirectory() as d:
    # atomic json + sha-stream vs one-shot hash
    p = os.path.join(d, "bank.json")
    R.atomic_write_json(p, {"b": 2, "a": 1})
    ok("atomic_write round-trip", json.load(open(p)) == {"a": 1, "b": 2})
    ok("no tmp litter", not [f for f in os.listdir(d) if f.endswith(".tmp")])
    ok("sha_stream == sha256", R.sha_stream(p) ==
       hashlib.sha256(open(p, "rb").read()).hexdigest())
    # code_sha: deterministic, order-insensitive, content-sensitive
    q = os.path.join(d, "code2.py"); open(q, "w").write("y=2\n")
    p2 = os.path.join(d, "code1.py"); open(p2, "w").write("x=1\n")
    s1 = R.code_sha(["code1.py", "code2.py"], root=d)
    ok("code_sha order-insensitive", s1 == R.code_sha(["code2.py", "code1.py"], root=d))
    open(p2, "w").write("x=999\n")
    ok("code_sha content-sensitive", s1 != R.code_sha(["code1.py", "code2.py"], root=d))
    # jsonl append + keyed merge (update wins, sorted, atomic)
    jl = os.path.join(d, "bank.jsonl")
    R.append_jsonl(jl, {"id": "b", "v": 1}); R.append_jsonl(jl, {"id": "a", "v": 1})
    ok("append_jsonl", len(R.read_jsonl(jl)) == 2)
    R.merge_jsonl(jl, [{"id": "b", "v": 9}, {"id": "c", "v": 1}], key="id")
    got = R.read_jsonl(jl)
    ok("merge keyed+sorted", [e["id"] for e in got] == ["a", "b", "c"]
       and got[1]["v"] == 9)
    # receipt schema law
    cs = R.code_sha(["code1.py"], root=d)
    rec = R.make_receipt("terrier-test-v1", "COMPLETE",
                         {"code_sha256": cs}, payload=[1, 2])
    rp = os.path.join(d, "receipt_S0.json"); R.write_receipt(rp, rec)
    ok("complete + pin match", R.receipt_complete(rp, cs))
    ok("stale code-sha refused", not R.receipt_complete(rp, "deadbeef"))
    try:
        R.make_receipt("s", "VOID_CONTROL_FAIL", {"code_sha256": cs}); bad = False
    except AssertionError:
        bad = True
    ok("non-COMPLETE needs reason", bad)
    try:
        R.make_receipt("s", "COMPLETE", {}); bad = False
    except AssertionError:
        bad = True
    ok("pins need code_sha256", bad)
    v = R.make_receipt("s", "PENDING-DATA", {"code_sha256": cs}, reason="no card")
    R.write_receipt(os.path.join(d, "receipt_S1.json"), v)
    done, todo = R.resume_split(["S0", "S1", "S2"],
                                lambda u: os.path.join(d, "receipt_%s.json" % u), cs)
    ok("resume split", done == ["S0"] and todo == ["S1", "S2"])
print("selftest_receipts: %d/%d green" % (n[0], n[0]))
