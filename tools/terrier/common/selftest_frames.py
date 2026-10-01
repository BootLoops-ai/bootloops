#!/usr/bin/env python3
"""battery: frames.py — registry rows, evidence classes, toggles, freeze/load."""
import os, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from frames import Registry

n = [0]
def ok(name, cond):
    n[0] += 1
    print("[%s] %s" % ("PASS" if cond else "FAIL", name))
    assert cond, name

r = Registry("demo-project")
r.add("A1", "W = sqrt(2/pi) Int G3^Omega", "PINNED-V+N", group="normalization",
      evidence=[{"class": "SRC", "anchor": "1912 eqs.(1),(4)"},
                {"class": "PIL", "anchor": "T1 4 digits"}])
r.add("B9", "xi IN e^{-Kcs}", "PINNED-N", group="periods",
      evidence=[{"class": "REF", "anchor": "1D valley toggle"}])
r.add_toggle("B9", "xi-in", "2.0371e-8", "DKMM 2.037e-8")
r.add_toggle("B9", "xi-out", "2.0482e-8", "Broeckel/CMS")
r.add_toggle("B9", "xi-wrong-sign", "2.0595e-8", "NOTHING (trap)")
r.add("C5", "Q_O = chi_f/2 per geometry", "PIN-REQ", group="tadpole")
try:
    r.add("Z1", "bad", "FROZEN-ISH"); bad = False
except AssertionError:
    bad = True
ok("status vocabulary enforced", bad)
try:
    r.add("Z2", "bad ev", "OURS", evidence=[{"class": "WIKI", "anchor": "x"}]); bad = False
except AssertionError:
    bad = True
ok("evidence class enforced", bad)
ok("unpinned lists PIN-REQ", r.unpinned() == ["C5"])
r.require_pinned(["A1", "B9"])
try:
    r.require_pinned(["A1", "C5", "NOPE"]); bad = False
except AssertionError as e:
    bad = "C5" in str(e) and "NOPE" in str(e)
ok("require_pinned fail-closed", bad)
with tempfile.TemporaryDirectory() as d:
    p = os.path.join(d, "FREEZE.json")
    sha = r.freeze(p)
    try:
        r.add("Z3", "late", "OURS"); bad = False
    except AssertionError:
        bad = True
    ok("frozen registry rejects rows", bad)
    r2 = Registry.load(p)
    ok("load verifies sha", r2.frozen_sha == sha and len(r2.rows) == 3)
    ok("toggles survive round-trip", len(r2.rows["B9"]["toggles"]) == 3)
    raw = open(p).read().replace("2.0371e-8", "2.0372e-8")
    open(p, "w").write(raw)
    try:
        Registry.load(p); bad = False
    except AssertionError:
        bad = True
    ok("tampered freeze refused", bad)
md = r.to_markdown()
ok("markdown template", "| A1 |" in md and "PIN-REQ" in md and "xi-wrong-sign" in md)
print("selftest_frames: %d/%d green" % (n[0], n[0]))
