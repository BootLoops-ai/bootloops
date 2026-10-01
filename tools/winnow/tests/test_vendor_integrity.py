#!/usr/bin/env python3
"""Vendor integrity: the vendored RECEIPT files in ibplapper/receipt/ must
be BYTE-IDENTICAL to the canonical tools/trust/receipt copies in this
repository (the receipt member of the trust package)
(import-never-fork rule), and VENDOR_MANIFEST.json must match the same
canonical bytes, so drift on either side is loud.
"""
import hashlib
import os
import sys

from _common import check, finish, PKG

RECEIPT_SRC = os.path.normpath(os.path.join(PKG, os.pardir, "trust", "receipt"))
V = os.path.join(PKG, "ibplapper")


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


for f in ("core.py", "emitter.py", "detector.py", "WITNESS_FORMAT.md"):
    check(f"receipt_vendor_byte_identical::{f}",
          sha(os.path.join(RECEIPT_SRC, f))
          == sha(os.path.join(V, "receipt", f)))

# vendor manifest still matches the canonical sources
import json                                              # noqa: E402
man = json.load(open(os.path.join(V, "receipt", "VENDOR_MANIFEST.json")))
for f, h in man["files"].items():
    check(f"vendor_manifest_matches_canonical::{f}",
          sha(os.path.join(RECEIPT_SRC, f)) == h)

finish("test_vendor_integrity")
