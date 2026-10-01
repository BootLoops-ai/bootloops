#!/usr/bin/env python3
"""b3_entropy.py — B3 6-var entropy certificate over the first census box.
Radius law (README.md in this directory; EXACT Cauchy-Schwarz on multinom):
  sum_{|k|=m} c(k) |phi|^k <= s^{2m},  s = sum_I sqrt|phi_I|  => proven
  convergence for s < 1; per-order term ratio = s^2. K4 law:
  s^2 >= 0.9 anywhere on the box => stage-1 summation dead there.
Rows certified (ball-separated comparisons, flint prec 200; UNDECIDED=>FAIL):
  (a) 64 extreme tile centers (c0+-4)/1024 + (4/1024)i  [pilot instrument]
  (b) 65 Stage-A points: center + 64 real box corners (c0+-8)/1024
  (c) box supremum, square hull per real axis (Re c0+-8, Im +-8)/1024
      + polydisc hull row (|phi_I - c0_I| <= 8/1024)
Also: J_d(51) rho-jet majorant constants + closed-form tail table at M=50."""
import time
from fractions import Fraction as Fr
from flint import arb, ctx
import towers6 as tw
from towers6 import CQ

ctx.prec = 200
C0 = (11, 13, 15, 17, 19, 21)
K4 = Fr(9, 10)


def s2_of(sqmods):
    s = arb(0)
    for q in sqmods:
        s += (arb(q.numerator) / q.denominator).sqrt().sqrt()
    return s * s


def verdict(s2):
    """ball-separated K4 comparison; UNDECIDED => FAIL (fail-closed)."""
    thr = arb(K4.numerator) / K4.denominator
    if s2 < thr:
        return "PASS"
    if s2 > thr:
        return "FAIL"
    return "UNDECIDED->FAIL"


def main():
    t0 = time.time()
    rows = []
    # (a) 64 extreme tile centers (complex)
    for kbit in range(64):
        q = [Fr((C0[i] + (4 if (kbit >> i) & 1 else -4)) ** 2 + 16, 1024 ** 2)
             for i in range(6)]
        rows.append((f"tilectr_{kbit:02d}", s2_of(q)))
    # (b) 65 Stage-A points (real)
    rows.append(("stageA_center", s2_of([Fr(c ** 2, 1024 ** 2) for c in C0])))
    for kbit in range(64):
        q = [Fr((C0[i] + (8 if (kbit >> i) & 1 else -8)) ** 2, 1024 ** 2)
             for i in range(6)]
        rows.append((f"stageA_corner_{kbit:02d}", s2_of(q)))
    # (c) box suprema
    rows.append(("boxsup_square_hull",
                 s2_of([Fr((c + 8) ** 2 + 64, 1024 ** 2) for c in C0])))
    rows.append(("boxsup_polydisc",
                 s2_of([Fr((c + 8) ** 2, 1024 ** 2) for c in C0])))
    npass = nfail = 0
    worst = {}
    for tag, s2 in rows:
        v = verdict(s2)
        npass += v == "PASS"
        nfail += v != "PASS"
        fam = tag.rsplit("_", 1)[0] if tag[-1].isdigit() else tag
        if fam not in worst or float(s2) > float(worst[fam][1]):
            worst[fam] = (tag, s2, v)
        print(f"{tag:22s} s2={str(s2)[:28]:30s} K4(0.9): {v}", flush=True)
    print(f"\n[K4] rows={len(rows)} PASS={npass} FAIL={nfail}", flush=True)
    for fam, (tag, s2, v) in worst.items():
        print(f"[worst] {fam:20s} {tag}: s2={float(s2):.6f} margin_to_0.9="
              f"{0.9 - float(s2):+.4f} {v}", flush=True)
    # J_d(51) constants + tail table at M=50 for the three worst scales
    Ht = tw.harmonic_tables(200)
    J = tw.J_rho_exact(51, Ht)
    print("\n[J_d(51)] rho-jet majorant constants:",
          [f"{float(x):.4e}" for x in J], flush=True)
    for fam in ("tilectr", "stageA_corner", "boxsup_square_hull"):
        tag, s2, v = worst[fam]
        th, rh = tw.tails_at(s2, 50, Ht, arb)
        fmt = lambda L: [(f"{float(x):.2e}" if x is not None else "HYP-FAIL")
                         for x in L]
        print(f"[tails M=50] {tag} (s2={float(s2):.4f}): theta {fmt(th)} | "
              f"rho {fmt(rh)}", flush=True)
    print(f"\n[done] wall {time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
