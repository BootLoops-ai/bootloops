#!/usr/bin/env python3
"""
Tests for tools/wayfinder/manifest.py.

All mutation happens on scratch copies inside a tempdir (standing rule:
never mutation-test in place).

Run:  python3 tests/test_manifest.py        (or pytest tests/test_manifest.py)
"""
import hashlib
import json
import os
import shutil
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from manifest import sha256_file, write_manifest, check_manifest


def test_sha256_file_matches_hashlib():
    with tempfile.TemporaryDirectory() as td:
        p = os.path.join(td, "a.bin")
        payload = b"wayfinder manifest test\n" * 1000 + b"\x00\xff tail"
        with open(p, "wb") as f:
            f.write(payload)
        want = hashlib.sha256(payload).hexdigest()  # definitionally the reference
        got = sha256_file(p)
        assert got == want, (got, want)
    print("  PASS sha256_file == hashlib.sha256 on same bytes")


def test_roundtrip_clean():
    with tempfile.TemporaryDirectory() as td:
        a = os.path.join(td, "a.txt")
        sub = os.path.join(td, "sub")
        os.makedirs(sub)
        b = os.path.join(sub, "b.bin")
        with open(a, "w") as f:
            f.write("alpha\n")
        with open(b, "wb") as f:
            f.write(os.urandom(4096))
        out = os.path.join(td, "MANIFEST.json")
        man = write_manifest([a, b], out)
        # keys are relative to the manifest dir
        assert set(man) == {"a.txt", os.path.join("sub", "b.bin")}, set(man)
        on_disk = json.load(open(out))
        assert on_disk == man
        for rec in man.values():
            assert len(rec["sha256"]) == 64
            assert rec["bytes"] > 0
            assert "T" in rec["mtime_iso"]
        assert check_manifest(out) == []
    print("  PASS write/check round-trip clean")


def test_detects_mutation_size_and_content():
    with tempfile.TemporaryDirectory() as td:
        a = os.path.join(td, "a.txt")
        b = os.path.join(td, "b.bin")
        with open(a, "w") as f:
            f.write("alpha\n")
        with open(b, "wb") as f:
            f.write(b"0123456789abcdef")
        out = os.path.join(td, "MANIFEST.json")
        write_manifest([a, b], out)
        # size change
        with open(a, "a") as f:
            f.write("extra\n")
        bad = check_manifest(out)
        assert any(m.startswith("SIZE a.txt") for m in bad), bad
        # same-size content flip (the case mtime/size checks miss)
        with open(b, "r+b") as f:
            f.seek(3)
            f.write(b"X")
        bad = check_manifest(out)
        assert any(m.startswith("SHA256 b.bin") for m in bad), bad
        assert len(bad) == 2, bad
    print("  PASS SIZE + same-size SHA256 mutations detected")


def test_detects_missing_and_restores_clean():
    with tempfile.TemporaryDirectory() as td:
        a = os.path.join(td, "a.txt")
        with open(a, "w") as f:
            f.write("alpha\n")
        out = os.path.join(td, "MANIFEST.json")
        write_manifest([a], out)
        stash = os.path.join(td, "stash")
        shutil.move(a, stash)
        bad = check_manifest(out)
        assert bad == ["MISSING a.txt"], bad
        shutil.move(stash, a)
        assert check_manifest(out) == []
    print("  PASS MISSING detected; restore -> clean")


def test_duplicate_key_and_missing_input():
    with tempfile.TemporaryDirectory() as td:
        a = os.path.join(td, "a.txt")
        with open(a, "w") as f:
            f.write("x")
        out = os.path.join(td, "MANIFEST.json")
        try:
            write_manifest([a, a], out)
            raise AssertionError("expected ValueError on duplicate key")
        except ValueError:
            pass
        try:
            write_manifest([os.path.join(td, "nope.bin")], out)
            raise AssertionError("expected FileNotFoundError")
        except FileNotFoundError:
            pass
        assert not os.path.exists(out + ".tmp"), "torn tmp file left behind"
    print("  PASS duplicate-key + missing-input errors")


def test_real_artifact_smoke():
    """Manifest one real reference artifact (read-only!) if present."""
    # Reference fixture, not shipped (leg SKIPs when absent).
    cone = os.path.join(os.environ.get("WAYFINDER_REFERENCE_FIXTURES", ""),
                        "cone_As_t13.json")
    if not os.path.isfile(cone):
        print("  SKIP (cone_As_t13.json missing)")
        return
    with tempfile.TemporaryDirectory() as td:
        out = os.path.join(td, "MANIFEST.json")
        man = write_manifest([cone], out)
        assert list(man) == [cone], list(man)  # different root -> absolute key
        assert check_manifest(out) == []
    print("  PASS real-artifact smoke (absolute-key fallback)")


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        print(t.__name__)
        t()
    print(f"ALL {len(tests)} manifest test groups done")


if __name__ == "__main__":
    main()
