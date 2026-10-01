#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision.
"""
Maintenance tool: re-pin the vendored fixtures after an intentional edit.

Every fixture under tests/fixtures/<row>/ is pinned by sha256 in three places,
and every battery refuses a fixture whose bytes drifted from its pins:

  1. tests/fixtures/<row>/PROVENANCE.json   files[].sha256 and files[].size
  2. tests/fixtures/CENSUS_BATTERY.json     "vendored_sha16" of every record whose
                                            "fixture" names the file
  3. the PINS / FIX dict literals in tests/test_*.py and topology_audit.py
     ("<row>/<file>": "<first 16 hex>" or "<row>/<file>": "<64 hex>")

This script recomputes the digest and size of every file PROVENANCE.json lists
and writes the new values into all three places.  record_sha256 / record_size
(the digest of the original corpus file the fixture was copied from) are never
touched, so an edited fixture stays marked as a modified copy.

Usage (from the package directory):
    python3 tests/repin_fixtures.py           # re-pin; prints every file whose pin changed
    python3 tests/repin_fixtures.py --check   # report drift only; exit 1 if any pin is stale

Run the batteries afterwards (python3 topology_audit.py --self-test; python3 -m
pytest tests): a re-pin makes the pins agree with the bytes on disk, it does not
make an edit correct.
"""
import hashlib
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
FIX = os.path.join(HERE, "fixtures")
MANIFEST = os.path.join(FIX, "CENSUS_BATTERY.json")
CODE_FILES = sorted(
    [os.path.join(HERE, f) for f in os.listdir(HERE) if f.startswith("test_") and f.endswith(".py")]
    + [os.path.join(HERE, "census_battery.py"), os.path.join(PKG, "topology_audit.py")])


def _sha256(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def _dump(path, obj):
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(obj, f, indent=1)
        f.write("\n")
    os.replace(tmp, path)


def scan():
    """{rel: {"sha256", "size", "old_sha256", "old_size", "prov"}} for every
    file a PROVENANCE.json lists; raises on a listed file that is absent."""
    out = {}
    for row in sorted(os.listdir(FIX)):
        prov_path = os.path.join(FIX, row, "PROVENANCE.json")
        if not os.path.isfile(prov_path):
            continue
        prov = json.load(open(prov_path))
        for ent in prov.get("files", []):
            rel = row + "/" + ent["path"]
            path = os.path.join(FIX, rel)
            if not os.path.isfile(path):
                raise SystemExit("listed in %s/PROVENANCE.json but absent: %s" % (row, rel))
            out[rel] = {"sha256": _sha256(path), "size": os.path.getsize(path),
                        "old_sha256": ent.get("sha256"), "old_size": ent.get("size"),
                        "prov": prov_path}
    return out


def repin_provenance(pins, write):
    changed = []
    by_prov = {}
    for rel, p in pins.items():
        by_prov.setdefault(p["prov"], []).append(rel)
    for prov_path, rels in sorted(by_prov.items()):
        prov = json.load(open(prov_path))
        row = os.path.basename(os.path.dirname(prov_path))
        dirty = False
        for ent in prov["files"]:
            rel = row + "/" + ent["path"]
            p = pins[rel]
            if ent.get("sha256") != p["sha256"] or ent.get("size") != p["size"]:
                changed.append((rel, ent.get("sha256"), p["sha256"]))
                ent["sha256"], ent["size"] = p["sha256"], p["size"]
                dirty = True
        if dirty and write:
            _dump(prov_path, prov)
    return changed


def repin_manifest(pins, write):
    if not os.path.isfile(MANIFEST):
        return []
    man = json.load(open(MANIFEST))
    changed = []

    def walk(node):
        if isinstance(node, dict):
            if "fixture" in node and "vendored_sha16" in node and node["fixture"] in pins:
                new = pins[node["fixture"]]["sha256"][:16]
                if node["vendored_sha16"] != new:
                    changed.append((node["fixture"], node["vendored_sha16"], new))
                    node["vendored_sha16"] = new
            for v in node.values():
                walk(v)
        elif isinstance(node, list):
            for v in node:
                walk(v)
    walk(man)
    if changed and write:
        _dump(MANIFEST, man)
    return changed


def repin_code(pins, write):
    changed = []
    for path in CODE_FILES:
        src = open(path).read()
        new_src = src
        for rel, p in pins.items():
            pat = re.compile(r'("' + re.escape(rel) + r'"\s*:\s*\n?\s*")([0-9a-f]{64}|[0-9a-f]{16})(")')

            def sub(m, p=p, rel=rel):
                old = m.group(2)
                new = p["sha256"] if len(old) == 64 else p["sha256"][:16]
                if old != new:
                    changed.append((os.path.relpath(path, PKG), rel, old, new))
                return m.group(1) + new + m.group(3)
            new_src = pat.sub(sub, new_src)
        if new_src != src and write:
            tmp = path + ".tmp"
            with open(tmp, "w") as f:
                f.write(new_src)
            os.replace(tmp, path)
    return changed


def main(argv):
    check = "--check" in argv
    write = not check
    pins = scan()
    c1 = repin_provenance(pins, write)
    c2 = repin_manifest(pins, write)
    c3 = repin_code(pins, write)
    verb = "stale" if check else "re-pinned"
    for rel, old, new in c1:
        print("PROVENANCE  %s  %s: %s -> %s" % (verb, rel, (old or "-")[:16], new[:16]))
    for rel, old, new in c2:
        print("MANIFEST    %s  %s: %s -> %s" % (verb, rel, old, new))
    for f, rel, old, new in c3:
        print("CODE        %s  %s [%s]: %s -> %s" % (verb, f, rel, old[:16], new[:16]))
    n = len(c1) + len(c2) + len(c3)
    print("%d fixture files scanned; %d pin(s) %s" % (len(pins), n, verb))
    return 1 if (check and n) else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
