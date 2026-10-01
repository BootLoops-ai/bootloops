#!/usr/bin/env python3
"""emit_receipt.py — mechanical script-snapshot receipt emitter (vendored for
the deform_transport member): an emitter that writes a receipt should
snapshot the script it ran, not trust the file to sit still.

Every receipt JSON is emitted THROUGH this helper (or an equivalent that
does all three steps).  It (1) SNAPSHOTS the emitting script's bytes at
emission time (<member_dir>/code_snapshots/<name>.<sha16>.py —
content-addressed, idempotent), (2) writes the producer block pinning the
SNAPSHOT's sha (which can never drift, because the snapshot is immutable by
construction), and (3) refuses to emit if the producer block would be
incomplete.

Usage from an emitting script:
    from emit_receipt import emit
    emit(receipt_dict, out_path, component="MYCOMPONENT",
         script_path=__file__, member_dir="/path/to/member/dir")

The producer block written:
    producer.component / .script (the ORIGINAL path) / .script_sha256
    (of the bytes AS RUN) / .script_snapshot (relative snapshot path) /
    .stamp_utc
Verification for any future reader: sha256(snapshot bytes) == script_sha256,
always — the drift class is structurally extinct for receipts emitted here.
"""
import json, os, sys, hashlib, subprocess


def _sha256_file(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def emit(receipt, out_path, component, script_path, member_dir):
    if not isinstance(receipt, dict):
        raise SystemExit("emit_receipt: receipt must be a dict")
    sha = _sha256_file(script_path)
    outdir = os.path.dirname(os.path.abspath(out_path))
    if outdir:
        os.makedirs(outdir, exist_ok=True)
    snapdir = os.path.join(member_dir, "code_snapshots")
    os.makedirs(snapdir, exist_ok=True)
    base = os.path.basename(script_path)
    snap = os.path.join(snapdir, f"{base}.{sha[:16]}.py")
    if not os.path.exists(snap):
        with open(script_path, "rb") as f_in, open(snap, "wb") as f_out:
            f_out.write(f_in.read())
    # verify the snapshot before pinning it
    if _sha256_file(snap) != sha:
        raise SystemExit("emit_receipt: snapshot sha mismatch — refusing to emit")
    stamp = subprocess.run(["date", "-u", "+%Y-%m-%dT%H:%M:%SZ"],
                           capture_output=True, text=True).stdout.strip()
    receipt["producer"] = {
        "component": component,
        "script": script_path, "script_sha256": sha,
        "script_snapshot": os.path.relpath(snap, member_dir),
        "stamp_utc": stamp,
    }
    # inline completeness refusal (producer_lint --check semantics)
    for k in ("component", "script", "script_sha256", "stamp_utc"):
        if not receipt["producer"].get(k):
            raise SystemExit(f"emit_receipt: producer.{k} empty — refusing to emit")
    with open(out_path, "w") as f:
        json.dump(receipt, f, indent=1)
    return out_path


if __name__ == "__main__":
    print(__doc__)
