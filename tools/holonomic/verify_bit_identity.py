# verify_bit_identity.py — RELOCATED-SAGEENV member (holonomic family):
# the relocation-correctness bit-identity gate. Run on BOTH the origin and
# the relocated env (HOLONOMIC_DIR=<this dir>); the emitted VERIFY_SHA256
# lines must be EQUAL — deterministic 'naive' summation makes bit-identity
# the right gate grade.
"""verify_bit_identity.py — relocation-correctness gate, scaled-down control.

A fixed toy certified transport whose FULL-PRECISION serialization must be
bit-identical between the origin env and the relocated env. The
relocate_sageenv.sh helper runs this file on both envs and compares the
emitted VERIFY_SHA256 lines.

Toy: L = u(1-u) Du^2 + (1-2u) Du - 1/4 (2F1-type), transport
[1/3 -> 1/2], ics [1, 1/7] (exact rationals), eps 1e-60, prec 430,
algorithm='naive' (the reference build's only working path — see
SAGEENV_PIN.md Sec. 3).
"""
import hashlib
import os
import sys

sys.path.insert(0, os.environ["HOLONOMIC_DIR"])
import holonomic_transport as H                    # noqa: E402
from sage.all import QQ, PolynomialRing, RealBallField  # noqa: E402
from ore_algebra import OreAlgebra                 # noqa: E402

prec = 430
RBF = RealBallField(prec)
Ru = PolynomialRing(QQ, "u"); u = Ru.gen()
A = OreAlgebra(Ru, "Du"); Du = A.gen()
L = (u * (1 - u)) * Du ** 2 + (1 - 2 * u) * Du - QQ(1) / 4

res = H.certified_transport(L, QQ(1) / 3, QQ(1) / 2,
                            [RBF(1), RBF(QQ(1) / 7)],
                            eps=1e-60, prec=prec)
sv, ov = H.engine_versions()
lines = ["engine %s ore_algebra %s" % (sv, ov),
         "conv %s detoured %s" % (res["conv"], res["detoured"]),
         "enclosure_ok %s" % res["enclosure_ok"]]
def ser_ball(b):
    """Full-precision serialization of a real or complex ball: exact
    rational midpoint(s) + radius repr (dyadic mids are exact -> the
    serialization is bit-faithful)."""
    if hasattr(b, "real") and not hasattr(b.mid(), "exact_rational"):
        br, bi = b.real(), b.imag()
        return "(%s + %s*I) rad (%s, %s)" % (
            br.mid().exact_rational(), bi.mid().exact_rational(),
            br.rad(), bi.rad())
    return "%s rad %s" % (b.mid().exact_rational(), b.rad())


for i, b in enumerate(res["vector"]):
    bb = b if hasattr(b, "mid") else RBF(b)
    lines.append("comp%d %s" % (i, ser_ball(bb)))
blob = "\n".join(str(x) for x in lines)
print(blob)
print("VERIFY_SHA256 %s" % hashlib.sha256(blob.encode()).hexdigest())
