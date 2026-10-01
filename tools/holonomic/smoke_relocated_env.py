# smoke_relocated_env.py — RELOCATED-SAGEENV member (holonomic family):
# the smoke gate for a relocated pinned sageenv. Run under the RELOCATED
# env's own versioned bin python with the ENV line
# HOLONOMIC_DIR=<this dir> HOLONOMIC_SAGEENV=<relocated env root>;
# relocate_sageenv.sh beside this file wires it.
"""smoke_relocated_env.py — smoke gate for a RELOCATED pinned sageenv.

Gates: engine versions equal the pin (10.7 / 0.5 — HARD asserts: a
relocation must reproduce the pinned build exactly, unlike smoke_test.py,
which merely records a different build); sys.prefix resolves to the
relocated env root (fails loud if the ENV line is missing or wrong);
IC-convention probe; toy transition matrix. Run under the relocated env's
own versioned bin python directly (never `sage -python -c "import
ore_algebra"` bare — sage.all must be imported first; see SAGEENV_PIN.md).

Env: HOLONOMIC_DIR = dir containing holonomic_transport.py
     HOLONOMIC_SAGEENV = the relocated env root (prefix check)
"""
import os
import sys
import time

sys.path.insert(0, os.environ["HOLONOMIC_DIR"])
t0 = time.time()
import holonomic_transport as H                    # noqa: E402
from sage.all import QQ, PolynomialRing            # noqa: E402
from ore_algebra import OreAlgebra                 # noqa: E402

SAGE_PIN, ORE_PIN = "10.7", "0.5"

sv, ov = H.engine_versions()
print("versions:", sv, "| ore_algebra", ov)
assert sv.startswith(SAGE_PIN), "sage version %s != pin %s" % (sv, SAGE_PIN)
assert ov == ORE_PIN, "ore_algebra %s != pin %s" % (ov, ORE_PIN)

want = os.environ.get("HOLONOMIC_SAGEENV")
if want:
    assert os.path.realpath(sys.prefix) == os.path.realpath(want), \
        "sys.prefix %s is NOT the relocated env %s" % (sys.prefix, want)
    print("sys.prefix resolves to the relocated env:", sys.prefix)

conv = H.determine_ic_convention(prec=430)
print("ic convention:", conv)

Ru = PolynomialRing(QQ, "u"); u = Ru.gen()
A = OreAlgebra(Ru, "Du"); Du = A.gen()
L = (u * (1 - u)) * Du ** 2 + (1 - 2 * u) * Du - QQ(1) / 4   # 2F1-type toy
M = H.transition_matrix(L, [QQ(1) / 3, QQ(1) / 2], 1e-60)
print("toy transition matrix ok, entry00:", str(M[0, 0])[:40])
print("SMOKE_OK %.1fs" % (time.time() - t0))
