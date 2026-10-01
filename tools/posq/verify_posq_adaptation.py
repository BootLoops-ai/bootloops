#!/usr/bin/env python3
"""verify_posq_adaptation.py -- run the posq engine (C kernel + certified
rule builder) as a consumer ADAPTATION through the shipped eras adversarial
battery (tools/eras/verify_eras_adversarial.py), per that
verifier's rerun law.

WHAT RUNS (see posq.py docstring for the adaptation):
  PHASE 0 -- posq domain gates (abort battery on any failure):
    G1 exact closed form at (709/2048, 10) == pinned RUNG1 surrogate
       Z_exact_str (sha256 pinned; EXACT Fraction identity, plus ln check).
    G2 packaged C kernel ball at (p0, lam0) CONTAINS the exact rational
       (relwidth kill-line 2^-100, the rung-1 law).
    G3 planted count corruption (0011: 6->5, 1110: 1->2 -- the pinned
       rung-1 surrogate plant) FIRES through the packaged kernel.
    G4 enclosure sanity: shell-1 c-bracket enclosure contains the exact
       value at the box center, the 4 exact corners, and 4 edge midpoints.
  PHASE 1 -- the eras battery: run_battery(PosqKernelAdaptation, shells,
    Ks=[0], plant_at=shell 1) -- ordering law, closed-box containment
    (corners/faces/band/2 fresh batches), vacuousness + point-precision
    gates, all plants (narrow x(1-1e-6), x0.5; shift width*1e-3, scale-
    adaptive abs; drop-dim for BOTH dims), precision probe (53-bit rule
    rebuild -> non-finite, caught), per-dim FD vs the exact closed-form
    derivative at gate 1.1e-12. result.ok is the pass signal.

SHELLS (vector (p, lambda), center (709/2048, 10); all inside the
monotone-bracket domain p < 1/2, lambda > 0):
    shell0  radii (1e-3, 5e-2)   measured relw ~0.42, width/spread ~6x
    shell1  radii (5e-3, 2.5e-1) measured relw ~2.7,  width/spread ~8x  [plants]
    shell2  radii (2e-2, 1.0)    measured relw ~47,   width/spread ~33x [FD box]
All informative under the 1e3x vacuous gate by construction (measured).

LAW (from the eras battery): a single-seed pass is NECESSARY, NOT
SUFFICIENT -- rerun with several --seed values. The run you make IS the
verification record; there is no stored log to consult.

Usage: python3 verify_posq_adaptation.py [--seed N] [--points N]
Exit: 0 PASS / 1 FAIL / 2 INDET (eras exit classes).
"""
import argparse
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ERAS = os.path.normpath(os.path.join(HERE, "..", "eras"))
sys.path.insert(0, HERE)
sys.path.insert(0, ERAS)

import posq                                    # noqa: E402
import verify_eras_adversarial as vb           # noqa: E402
from flint import arb, ctx                     # noqa: E402
from fractions import Fraction as Fr           # noqa: E402

SHELLS = [(("709/2048", "10"), ("1e-3", "5e-2")),
          (("709/2048", "10"), ("5e-3", "2.5e-1")),
          (("709/2048", "10"), ("2e-2", "1e0"))]
PLANT_AT = (1, 0)


def gates(adapt):
    """posq domain gates (PHASE 0). Returns ok bool; prints PASS/FAIL rows."""
    ok_all = True
    lines = []

    def chk(name, ok, got=""):
        nonlocal ok_all
        lines.append(f"{'PASS' if ok else 'FAIL'}  {name}  {got}")
        ok_all = ok_all and ok

    t0 = time.time()
    # G1: exact closed form == pinned rung-1 exact rational (sha pin)
    sha, Zx = posq.z_exact_sha256(adapt.exact)
    chk("G1 exact closed form sha256 == pinned RUNG1 Z_exact_str pin",
        sha == posq.Z_EXACT_STR_SHA256_PIN, sha[:16])
    old = ctx.prec
    ctx.prec = 320
    Zex = arb(Zx.numerator) / arb(Zx.denominator)
    lnZ = float(Zex.log().mid())
    ctx.prec = old
    chk("G1b ln Z_exact matches pinned -25.97566273247789",
        abs(lnZ - posq.LN_Z_EXACT_PIN) < 1e-9, f"{lnZ:.14f}")

    # G2: kernel containment at (p0, lam0)
    Zk = adapt._kernel_Z(posq.P0, posq.LAM0, 256)
    ctx.prec = 320
    relw2 = float((Zk.rad() / Zk.mid()).log() / arb(2).log())
    ctx.prec = old
    chk("G2 packaged kernel ball CONTAINS exact rational at (p0,lam0)",
        bool(Zk.contains(Zex)), f"relw2={relw2:.1f}")
    chk("G2b kernel relwidth beats the 2^-100 kill line", relw2 < -100)

    # G3: planted count corruption fires through the packaged kernel
    vec_bad = list(posq.SURR_VEC16)
    vec_bad[3] -= 1          # 0011: 6 -> 5
    vec_bad[14] += 1         # 1110: 1 -> 2
    import kernel_io as kio
    import bench_rung1_gauss as bg
    c = posq.c_of(posq.P0, posq.LAM0)
    rules = adapt._rules_at(c)
    base = os.path.join(adapt._tmp, "gate_plant")
    rp, jp, op = base + ".rules", base + ".job", base + ".out"
    ctx.prec = bg.PREC_BUILD + 64
    kio.write_rules(rp, {"u1": rules["u1"], "u2": rules["u2"],
                         "v3": rules["v3"]})
    kio.write_job(jp, "cat", 256, posq.P0.numerator, posq.P0.denominator,
                  0, 17, rp, op, [list(posq.SURR_VEC16), vec_bad])
    kio.run_kernel(jp, timeout=600)
    _, balls = kio.read_out(op)
    ctx.prec = 320
    cA = arb(c.numerator) / arb(c.denominator)
    Zh = cA ** 3 * balls[0] / 2
    Zb = cA ** 3 * balls[1] / 2
    dln = float((Zb.log() - Zh.log()).mid())
    ctx.prec = old
    for f in (rp, jp, op):
        os.remove(f)
    chk("G3 planted count corruption (0011:-1, 1110:+1) FIRES via kernel",
        abs(dln) > 1e-3, f"dln={dln:.4f}")

    # G4: enclosure contains exact values at center/corners/edge midpoints
    (p0s, rs) = SHELLS[1]
    enc = vb.norm_out(adapt.enclose(p0s, rs, 0))[0]
    pc, lc = posq.frs(p0s[0]), posq.frs(p0s[1])
    rp_, rl_ = posq.frs(rs[0]), posq.frs(rs[1])
    pts = [(pc, lc), (pc - rp_, lc - rl_), (pc - rp_, lc + rl_),
           (pc + rp_, lc - rl_), (pc + rp_, lc + rl_),
           (pc - rp_, lc), (pc + rp_, lc), (pc, lc - rl_), (pc, lc + rl_)]
    ok9 = True
    ctx.prec = 320
    for (pp, lam) in pts:
        Zp = adapt.exact.Z_fr(pp, lam)
        ok9 = ok9 and bool(enc.contains(arb(Zp.numerator) / arb(Zp.denominator)))
    ctx.prec = old
    chk("G4 shell-1 enclosure contains exact value at center+corners+edges",
        ok9, f"9/9={ok9}")

    out = "\n".join(lines + [f"POSQ GATES: {'PASS' if ok_all else 'FAIL'} "
                             f"({time.time()-t0:.1f}s)"])
    print(out, flush=True)
    return ok_all


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=987654321202607)
    ap.add_argument("--points", type=int, default=100)
    ap.add_argument("--fd-tol", default="1.1e-12")
    a = ap.parse_args()
    t00 = time.time()
    import kernel_io as kio
    try:
        kio.kernel_path()      # resolve (or locally build) the C kernel first
    except kio.KernelUnavailable as e:
        print(f"VERIFY INDET -- {e}; battery not run")
        return 2
    print(f"posq adaptation battery: seed {a.seed}, {a.points} points/shell/"
          f"batch, fd-tol {a.fd_tol} (eras battery: {ERAS})", flush=True)
    print("adaptation: posq C kernel point engine (rules (17,13,5) rebuilt + "
          "receipt-gated per point at prec 2048, sweeps at job prec), "
          "monotone c-bracket x ball-q closed-box enclosures, exact "
          "closed-form derivative jets; object = pinned rung-1 N=8 "
          "exact-rational surrogate; vector (p, lambda)", flush=True)
    adapt = posq.PosqKernelAdaptation()
    print("\n#### PHASE 0: posq domain gates ####", flush=True)
    if not gates(adapt):
        print("\nVERIFY FAIL -- posq domain gates failed; battery not run")
        return 1
    print("\n#### PHASE 1: eras adversarial battery ####", flush=True)
    cfg = vb.Cfg(seed=a.seed, points=a.points, fd_tol=a.fd_tol)
    res = vb.run_battery(adapt, SHELLS, [0], cfg, plant_at=PLANT_AT)
    print(f"\nNOTES ({len(vb.NOTES)}):")
    for n in vb.NOTES:
        print("  " + n)
    print(f"\nFATALS ({len(vb.FATALS)}):")
    for f in vb.FATALS:
        print("  " + f)
    import resource
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    if sys.platform == "darwin":
        rss /= 1024  # ru_maxrss is bytes on macOS, KiB on Linux
    ch = resource.getrusage(resource.RUSAGE_CHILDREN)
    cpu = (resource.getrusage(resource.RUSAGE_SELF).ru_utime
           + resource.getrusage(resource.RUSAGE_SELF).ru_stime
           + ch.ru_utime + ch.ru_stime)
    if vb.FATALS:
        verdict, code = "FAIL", 1
    elif not res.ok:
        verdict, code = "INDET (result.ok False without fatals)", 2
    else:
        verdict, code = "PASS", 0
    print(f"\nVERIFY {verdict} -- {vb.N_PASS[0]} checks passed, "
          f"{len(vb.NOTES)} notes, {len(vb.FATALS)} fatals; "
          f"coverage {res.cells_pass} informative PASS cells, "
          f"{res.cells_indet} no-information, {res.cells_declared} declared; "
          f"result.ok={res.ok} "
          f"({time.time()-t00:.1f}s wall, {cpu:.1f}s cpu incl. kernel, "
          f"maxrss {rss:.1f} MB)")
    return code


if __name__ == "__main__":
    sys.exit(main())
