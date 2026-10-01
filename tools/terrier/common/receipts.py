#!/usr/bin/env python3
"""receipts.py — terrier common chassis: receipt schema + atomic write/resume/sha.

One implementation of the receipt discipline; both wings import this instead
of re-growing it: atomic tmp+fsync+rename writes, receipt-as-checkpoint
resume, code-sha stamps, keyed jsonl merge, an fsync'd append-only log, and
per-unit receipts whose non-COMPLETE status carries an embedded reason
(PENDING-DATA / DROPPED-BUDGET). Battery: selftest_receipts.py. stdlib only.
"""
import hashlib, json, os, tempfile, time

# receipt law: the receipt IS the checkpoint. A unit is done iff its
# receipt parses, status is COMPLETE, and its code-sha pin matches running code.
STATUS_OK = ("COMPLETE",)
STATUS_BAD = ("VOID_CONTROL_FAIL", "QUARANTINED_SURVIVOR", "PENDING-DATA",
              "DROPPED-BUDGET")


def sha_stream(path, chunk=1 << 20):
    """sha256 of a file, streamed — never slurp multi-GB files."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for blk in iter(lambda: f.read(chunk), b""):
            h.update(blk)
    return h.hexdigest()


def code_sha(files, root=None):
    """sha256 over the sorted concatenation of the code files."""
    h = hashlib.sha256()
    for fn in sorted(files):
        p = os.path.join(root, fn) if root else fn
        with open(p, "rb") as f:
            h.update(f.read())
    return h.hexdigest()


def utc_now():
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def atomic_write_json(path, obj):
    """tmp + fsync + rename in the destination dir (crash-safe checkpoint)."""
    path = os.path.abspath(path)
    fd, tmp = tempfile.mkstemp(dir=os.path.dirname(path), suffix=".tmp")
    with os.fdopen(fd, "w") as f:
        json.dump(obj, f, indent=1, sort_keys=True, default=str)
        f.flush()
        os.fsync(f.fileno())
    os.rename(tmp, path)


def append_jsonl(path, obj):
    """fsync'd append-only jsonl — for logs, never checkpoints."""
    with open(path, "a") as f:
        f.write(json.dumps(obj, sort_keys=True, default=str) + "\n")
        f.flush()
        os.fsync(f.fileno())


def read_jsonl(path):
    if not os.path.exists(path):
        return []
    return [json.loads(l) for l in open(path) if l.strip()]


def merge_jsonl(path, entries, key):
    """merge entries into a keyed jsonl file, atomically."""
    byid = {e[key]: e for e in read_jsonl(path)}
    for e in entries:
        byid[e[key]] = e
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        for k in sorted(byid):
            f.write(json.dumps(byid[k], sort_keys=True, default=str) + "\n")
        f.flush()
        os.fsync(f.fileno())
    os.rename(tmp, path)
    return len(byid)


def make_receipt(schema, status, pins, **fields):
    """Receipt schema: schema + status + pins (MUST carry code_sha256) + payload.
    Non-COMPLETE terminal receipts must say why (reason= field)."""
    rec = dict(schema=schema, status=status, pins=dict(pins),
               generated_utc=utc_now(), **fields)
    validate_receipt(rec)
    return rec


def validate_receipt(rec):
    for k in ("schema", "status", "pins"):
        assert rec.get(k), "receipt missing %r" % k
    assert rec["pins"].get("code_sha256"), "receipt pins lack code_sha256 stamp"
    if rec["status"] not in STATUS_OK and not rec.get("reason"):
        raise AssertionError("non-COMPLETE receipt (%s) needs a reason= field"
                             % rec["status"])


def write_receipt(path, rec):
    validate_receipt(rec)
    atomic_write_json(path, rec)


def receipt_complete(path, expect_code_sha=None):
    """Resume law: done iff receipt parses, COMPLETE, and code-sha pin matches."""
    try:
        rec = json.load(open(path))
    except Exception:
        return False
    if rec.get("status") not in STATUS_OK:
        return False
    return (expect_code_sha is None or
            rec.get("pins", {}).get("code_sha256") == expect_code_sha)


def resume_split(unit_ids, receipt_path_fn, expect_code_sha=None):
    """(done, todo) partition by the receipt-as-checkpoint law."""
    done = [u for u in unit_ids
            if receipt_complete(receipt_path_fn(u), expect_code_sha)]
    todo = [u for u in unit_ids if u not in set(done)]
    return done, todo
