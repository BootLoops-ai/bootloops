"""smoke_test.py — tools/holonomic engine smoke (import + one gated example).

Run under Sage's Python (with ore_algebra installed into that Sage):

    sage -python smoke_test.py [receipt.json]

Checks:
  S1 engine import + versions (reference build: Sage 10.7 / ore_algebra 0.5 —
     recorded; a different build is recorded, not failed).
  S2 IC-convention probe returns a valid convention (recorded, never assumed).
  S3 gated example: L = (1+x)Dx^2 + Dx annihilates ln(1+x); certified
     transport [0, 1/2] at eps 1e-40 must CONTAIN ln(3/2) with radius <= 1e-38
     (reference run measured radius 2.9e-42 through the house layer; the raw
     transition-matrix entry certifies to 2.0e-47 — see reference/).
  S4 path-planning law: leading-coefficient root on the segment forces the
     complex-midpoint detour.
  S5 reference-build fact: the DEFAULT algorithm path raises TypeError
     (_small_primes_table signature) in the reference build — recorded; a
     build where it succeeds is reported as ENV-FACT-CHANGED, not a failure.
  S6 detoured VALUE gate: transport ln(1+x) from 0 to -3/2 — the leading
     root x=-1 sits on the segment, so the path detours through the upper
     half-plane and the continued value is genuinely complex, ln(1/2) + i*pi.
     The returned ball must be complex, CONTAIN that target (real AND
     imaginary parts), with radius <= 1e-38 — the imaginary part must
     survive, never be silently dropped.
"""
import json
import sys
import time

from sage.all import QQ, PolynomialRing, RealBallField, ComplexBallField
from ore_algebra import OreAlgebra

import holonomic_transport as H

OUT = sys.argv[1] if len(sys.argv) > 1 else "smoke_receipt.json"
rep = {"tool": "holonomic", "checks": {}}
t0 = time.time()
ok = True

# S1 versions
sv, ov = H.engine_versions()
rep["checks"]["S1_versions"] = {"sage": sv, "ore_algebra": ov,
                                "reference": {"sage": "10.7",
                                              "ore_algebra": "0.5"},
                                "match": sv.startswith("10.7") and ov == "0.5"}
print(f"[S1] sage {sv} / ore_algebra {ov}")

# S2 convention probe
conv = H.determine_ic_convention()
rep["checks"]["S2_ic_convention"] = conv
print(f"[S2] IC convention (determined, this session): {conv}")

# S3 gated example: ln(1+x)
Rx = PolynomialRing(QQ, "x")
x = Rx.gen()
A = OreAlgebra(Rx, "Dx")
Dx = A.gen()
RBF = RealBallField(430)
L = (1 + x) * Dx ** 2 + Dx
res = H.certified_transport(L, QQ(0), QQ(1) / 2,
                            [RBF(0), RBF(1)], eps=1e-40, conv=conv)
target = RBF(QQ(3) / 2).log()
diff = res["value"] - target
contains = bool(diff.contains_zero())
rad = float(res["value"].rad())
s3_ok = contains and rad <= 1e-38 and res["enclosure_ok"] \
    and not res["detoured"]
rep["checks"]["S3_ln32_gate"] = {
    "value": str(res["value"]), "target": "ln(3/2)",
    "contains_target": contains, "radius": rad, "radius_bar": 1e-38,
    "detoured": res["detoured"], "pass": s3_ok}
ok = ok and s3_ok
print(f"[S3] transport(ln(1+x), [0,1/2]) = {res['value']}  "
      f"contains ln(3/2)={contains}  rad={rad:.2e}  "
      f"{'PASS' if s3_ok else 'FAIL'}")

# S4 detour law: leading root at 1/2 inside [1/3, 1]
L2 = (2 * x - 1) * Dx ** 2 + Dx
path2, onpath2 = H.plan_path(L2, QQ(1) / 3, QQ(1))
s4_ok = len(onpath2) == 1 and len(path2) == 3
rep["checks"]["S4_detour_law"] = {"roots_on_path": [str(r) for r in onpath2],
                                  "path_len": len(path2), "pass": s4_ok}
ok = ok and s4_ok
print(f"[S4] root-on-segment detour: roots={onpath2} path_len={len(path2)}  "
      f"{'PASS' if s4_ok else 'FAIL'}")

# S5 default-algorithm build fact
try:
    L.numerical_transition_matrix([QQ(0), QQ(1) / 2], 1e-20)
    s5 = {"default_path": "SUCCEEDED",
          "note": "ENV-FACT-CHANGED: this build's compiled path works; all "
                  "recorded numbers in this package assume 'naive' — "
                  "re-measure before relying on the fast path"}
    print("[S5] default algorithm SUCCEEDED — ENV-FACT-CHANGED (see receipt)")
except TypeError as e:
    s5 = {"default_path": "TypeError (as in the reference build)",
          "error": str(e)[:160]}
    print(f"[S5] default algorithm raises TypeError as in the reference "
          f"build: {str(e)[:80]}")
except Exception as e:  # any other hard failure is still the broken path
    s5 = {"default_path": type(e).__name__, "error": str(e)[:160]}
    print(f"[S5] default algorithm raises {type(e).__name__}: {str(e)[:80]}")
rep["checks"]["S5_default_algorithm_probe"] = s5

# S6 detoured continuation keeps its imaginary part: ln(1+x) across x=-1
CBF = ComplexBallField(430)
res6 = H.certified_transport(L, QQ(0), QQ(-3) / 2,
                             [RBF(0), RBF(1)], eps=1e-40, conv=conv)
target6 = CBF(RBF(QQ(1) / 2).log(), RBF.pi())   # ln(1/2) + i*pi
val6 = res6["value"]
is_complex6 = hasattr(val6, "imag") and not val6.imag().is_zero()
contains6 = bool((CBF(val6) - target6).contains_zero())
rad6 = float(val6.rad())
s6_ok = res6["detoured"] and is_complex6 and contains6 \
    and rad6 <= 1e-38 and res6["enclosure_ok"]
rep["checks"]["S6_detoured_complex_value"] = {
    "value": str(val6), "target": "ln(1/2) + i*pi",
    "detoured": res6["detoured"], "value_is_complex": bool(is_complex6),
    "contains_target": contains6, "radius": rad6, "radius_bar": 1e-38,
    "pass": bool(s6_ok)}
ok = ok and s6_ok
print(f"[S6] detoured transport(ln(1+x), [0,-3/2]) = {val6}  "
      f"complex={is_complex6}  contains ln(1/2)+i*pi={contains6}  "
      f"rad={rad6:.2e}  {'PASS' if s6_ok else 'FAIL'}")

rep["wall_s"] = round(time.time() - t0, 2)
rep["stamp"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
rep["pass"] = bool(ok)
json.dump(rep, open(OUT, "w"), indent=1, default=str)
print(f"{'SMOKE PASS' if ok else 'SMOKE FAIL'}  ({rep['wall_s']}s)  "
      f"receipt -> {OUT}")
sys.exit(0 if ok else 1)
