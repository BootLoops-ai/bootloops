#!/usr/bin/env python3
# lockpick bilmine member — synthetic self-check: planted bilinear quadric recovered blind + committed negative NULL + ring-detector smoke. No external data, no env vars.
"""selftest.py — bilmine synthetic self-check (committed data only).

S1 planted quadric (C1-style, hand-provable): v0..v3 are committed
   sha512-derived complex numbers (~140 digits); v4 is constructed so that
   EXACTLY the primitive integer quadric

       3*(v0*v1) - 2*(v2*v3) + 1*(v4*v4) = 0

   holds (v4 = sqrt(2*v2*v3 - 3*v0*v1), principal branch). The blind mine
   over the 15 symmetric products must recover the canonical vector exactly
   on both precision legs (40/50 dps) and both engines (lattice + mp.pslq),
   with a clean audit and lattice rank 1.
S2 negative: the committed sha512-derived synthetic 5-vector (the C2
   construction), matched magnitudes — zero lattice hits at every lawful
   rung and pslq NULL, both legs.
S3 ring-detector smoke: the F1 stage's sealed absolute-floor Q(sqrt15)
   detector battery (5 committed cases) must pass 5/5.

Exit 0 PASS / 1 FAIL. Seconds-scale; safe from any cwd.
"""
import os, sys

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "f1"))
from run_control import mine_quadric, synthetic_pi, P5   # noqa: E402
from run_control_f1 import detector_smoke                # noqa: E402
from bilmine_lib import canonicalize                     # noqa: E402


def planted_vec5():
    """Committed 5-vector carrying exactly the quadric 3*v01 - 2*v23 + v44 = 0."""
    base = synthetic_pi([1, 1, 1, 1, 1])  # committed seeds BILMINE-NEG-1..5
    with mp.workdps(160):
        v = [mp.mpc(1, 0) + z for z in base[:4]]     # shift off the origin
        v4 = mp.sqrt(2 * v[2] * v[3] - 3 * v[0] * v[1])
        v.append(v4)
    return v


def want_vec():
    c = [0] * 15
    c[P5.index((0, 1))] = 3
    c[P5.index((2, 3))] = -2
    c[P5.index((4, 4))] = 1
    return list(canonicalize(c))


def lat_vecs(rec):
    return sorted(tuple(h["vec"]) for h in rec["lattice_hits"])


def main():
    ok_all = True

    # ---- S1: planted quadric, blind mine, both legs + both engines ----
    v = planted_vec5()
    want = tuple(want_vec())
    s1_ok = True
    for leg, d in (("lo", 40), ("hi", 50)):
        rec = mine_quadric(v, d, [10, 100], f"S1-{leg}")
        hits = lat_vecs(rec)
        p = rec["pslq_hit"]
        leg_ok = (rec["audit"]["audit_pass"]
                  and hits == [want]
                  and bool(p) and p["passes_both_rows"]
                  and tuple(p["vec"]) == want)
        print(f"S1-{leg}: audit={rec['audit']['audit_pass']} rank={len(hits)}"
              f" lattice_exact={hits == [want]}"
              f" pslq_exact={bool(p) and tuple(p['vec']) == want}"
              f" -> {'PASS' if leg_ok else 'FAIL'}")
        s1_ok = s1_ok and leg_ok
    ok_all = ok_all and s1_ok

    # ---- S2: committed negative must NULL everywhere ----
    with mp.workdps(160):
        mags = [max(abs(z), mp.mpf(1) / 4) for z in v]
    neg = synthetic_pi(mags)
    s2_ok = True
    for leg, d in (("lo", 40), ("hi", 50)):
        rec = mine_quadric(neg, d, [10, 100], f"S2-{leg}")
        p = rec["pslq_hit"]
        pslq_null = (not p) or (not p["passes_both_rows"])
        leg_ok = (len(rec["lattice_hits"]) == 0) and pslq_null
        print(f"S2-{leg}: lattice_hits={len(rec['lattice_hits'])}"
              f" pslq_null={pslq_null} -> {'PASS' if leg_ok else 'FAIL'}")
        s2_ok = s2_ok and leg_ok
    ok_all = ok_all and s2_ok

    # ---- S3: F1 ring-detector smoke (committed cases) ----
    det = detector_smoke()
    n_pass = sum(1 for t in det["tests"] if t["pass"])
    print(f"S3: detector smoke {n_pass}/{len(det['tests'])}"
          f" -> {'PASS' if det['pass'] else 'FAIL'}")
    ok_all = ok_all and det["pass"]

    print("bilmine selftest:", "PASS" if ok_all else "FAIL")
    return 0 if ok_all else 1


if __name__ == "__main__":
    sys.exit(main())
