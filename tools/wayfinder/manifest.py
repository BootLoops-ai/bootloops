#!/usr/bin/env python3
"""
manifest.py — sha256 manifests for wayfinder input artifacts.

Why this exists:
  * Transport runs consume big on-disk DE files (A_of_y.json, monomial-json
    connections, kira_target.m chains). A manifest written when results are
    saved lets a cold resume PROVE the inputs are the ones the saved digits
    were computed from, instead of trusting a status note.
  * After any interruption or in-place mutation risk, a checksum diff of the
    on-disk data is the only reliable audit; check_manifest is that diff.

Contract (wayfinder shared API):
    sha256_file(path) -> str
    write_manifest(paths, out_json) -> dict   # {relpath: {sha256, bytes, mtime_iso}}
    check_manifest(manifest_json) -> list[str]  # mismatch strings, [] == OK

Conventions:
  * relpath keys are relative to the DIRECTORY OF THE MANIFEST FILE, so a
    manifest travels with its tree (runs/<problem>/<tag>/). Paths that cannot
    be expressed relative (different root) are stored absolute.
  * 'mtime_iso' is informational only — check_manifest compares bytes+sha256,
    never mtime (rsync/cp -p churn must not raise false alarms).
  * Pure stdlib; no mpmath here, nothing precision-dependent.
"""
import hashlib
import json
import os
from datetime import datetime, timezone

__all__ = ["sha256_file", "write_manifest", "check_manifest"]

_CHUNK = 1 << 20  # 1 MiB read blocks: A_of_y.json is 11 MB, kira_target.m up to ~7 MB


def sha256_file(path) -> str:
    """Hex sha256 of a file, streamed (never loads the file whole)."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(_CHUNK), b""):
            h.update(chunk)
    return h.hexdigest()


def _rel_key(abspath, base):
    """Manifest key for abspath: relative to base if possible, else absolute."""
    try:
        rel = os.path.relpath(abspath, base)
    except ValueError:  # e.g. different drive
        return abspath
    # refuse confusing keys that climb far out of the manifest tree yet are
    # longer than the absolute path — keep them absolute for readability
    if rel.startswith(".." + os.sep) and len(rel) >= len(abspath):
        return abspath
    return rel


def write_manifest(paths, out_json) -> dict:
    """Hash `paths`, write {relpath: {sha256, bytes, mtime_iso}} to out_json.

    Returns the manifest dict. Raises FileNotFoundError on a missing input and
    ValueError if two inputs collide onto the same relpath key.
    """
    out_json = os.path.abspath(out_json)
    base = os.path.dirname(out_json) or "."
    man = {}
    for p in paths:
        ap = os.path.abspath(p)
        st = os.stat(ap)  # raises FileNotFoundError — fail loudly, never skip
        key = _rel_key(ap, base)
        if key in man:
            raise ValueError(f"duplicate manifest key {key!r} (from {p!r})")
        man[key] = {
            "sha256": sha256_file(ap),
            "bytes": st.st_size,
            "mtime_iso": datetime.fromtimestamp(st.st_mtime, timezone.utc)
            .isoformat(timespec="seconds"),
        }
    tmp = out_json + ".tmp"
    with open(tmp, "w") as f:
        json.dump(man, f, indent=1, sort_keys=True)
        f.write("\n")
    os.replace(tmp, out_json)  # atomic: a crashed writer never leaves a torn manifest
    return man


def check_manifest(manifest_json) -> list:
    """Re-hash every entry of a manifest. Returns list of mismatch strings.

    Empty list == everything verifies. Mismatch kinds:
      'MISSING <key>'                       file gone
      'SIZE <key>: got G expected E'        byte count differs (cheap pre-check)
      'SHA256 <key>: got g.. expected e..'  content differs at equal size
    mtime is NOT checked (informational field only).
    """
    manifest_json = os.path.abspath(manifest_json)
    base = os.path.dirname(manifest_json) or "."
    with open(manifest_json) as f:
        man = json.load(f)
    mismatches = []
    for key, rec in sorted(man.items()):
        path = key if os.path.isabs(key) else os.path.join(base, key)
        if not os.path.isfile(path):
            mismatches.append(f"MISSING {key}")
            continue
        size = os.stat(path).st_size
        if size != rec["bytes"]:
            mismatches.append(f"SIZE {key}: got {size} expected {rec['bytes']}")
            continue
        digest = sha256_file(path)
        if digest != rec["sha256"]:
            mismatches.append(
                f"SHA256 {key}: got {digest[:16]}.. expected {rec['sha256'][:16]}.."
            )
    return mismatches
