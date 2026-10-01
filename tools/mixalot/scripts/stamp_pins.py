"""Stamp _pins.py with the sha256 of each vendored file (build-time only;
run in the same write as any vendor refresh)."""
import hashlib
import os
import re

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VENDOR = os.path.join(HERE, "mixalot", "vendor")
PINS_PY = os.path.join(HERE, "mixalot", "_pins.py")


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


src = open(PINS_PY).read()
import importlib.util
spec = importlib.util.spec_from_file_location("_pins", PINS_PY)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
for rel, (source, _) in mod.PINS.items():
    h = sha(os.path.join(VENDOR, rel))
    pat = re.compile(r'("%s": \(\n?\s*"%s",\s*\n?\s*)(None|"[0-9a-f]{64}")'
                     % (re.escape(rel), re.escape(source)), re.S)
    src, n = pat.subn(lambda m: m.group(1) + '"%s"' % h, src)
    assert n == 1, rel
open(PINS_PY, "w").write(src)
print("stamped", len(mod.PINS), "pins")
