#!/usr/bin/env python3
r"""VALIDATION GATE: reproduce the reference c_M = 1 (the arXiv:2601.16256
memory-sector constant) through the GENERALIZED membound path (generic
sector-weight engine, no hard-coded B1/B2).

Modes:
  --micro          time ONE representative quad (heaviest: I3M order-2, incl.
                   kappa2) and print the projected total (cost-rate probe —
                   run this FIRST to size a full run).
  (default)        full reduced-precision gate: all cores, assembly, checks.

Checks:
  j_1^(0) == 1/30 exactly; I1 leading pole == 1/(15(8pi)^4);
  I2 eps^-4, eps^-3 cancellation; c_M == 1; PSLQ [c_M, 1].
Digits demanded scale with --dps (an honest reduced-precision gate; a
high-precision run — dps 30, maxdegree 6 — has passed the same checks with
the pole at 28+ digits and c_M = 1 at 25+ digits, ~14 min wall farmed over
12 processes).
"""
import argparse, json, os, resource, sys, time
import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from membound.core import compute_core, make_integrand, KCache, domain_classes
from membound.spec import I1M_CORE, I3M_CORE, REGISTRY
from membound.farm import compute_cores_parallel
from membound import assemble as asm

T0 = time.time()
# the two registry cores of the validated P/memory slot; any user core farms
# through the same membound.farm path (python3 -m membound --spec core.json)
_SPECS = {"I1M": I1M_CORE, "I3M": I3M_CORE}


def rss_gb():
    self_r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    child_r = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    # ru_maxrss is kilobytes on Linux, bytes on macOS
    div = 1024 ** 3 if sys.platform == "darwin" else 1024 ** 2
    return max(self_r, child_r) / div


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dps", type=int, default=11)
    ap.add_argument("--maxdegree", type=int, default=4)
    ap.add_argument("--panels", default="0,0.5,2,8,inf")
    ap.add_argument("--micro", action="store_true")
    ap.add_argument("--jobs", type=int, default=None,
                    help="worker processes (default: min(12, cpu count))")
    ap.add_argument("--json-out", default=None)
    args = ap.parse_args()
    if args.jobs is None:
        args.jobs = max(1, min(12, os.cpu_count() or 1))

    mp.mp.dps = args.dps
    brk = [mp.inf if p == "inf" else mp.mpf(p) for p in args.panels.split(",")]
    log = print
    log(f"# membound gate_cM  dps={args.dps} maxdegree={args.maxdegree} panels={brk}")
    kc = KCache(args.dps)
    log(f"# kappa2 FD: workdps={kc.W} h=1e-{kc.W // 6} (dps-scaled)")

    if args.micro:
        f = make_integrand(I3M_CORE, ('same', 1, 1), 2, kc)
        t0 = time.time()
        val = mp.quad(f, brk, brk, maxdegree=args.maxdegree)
        dt = time.time() - t0
        # count unique quads: per core, per order, #classes
        nq = 0
        for spec in (I1M_CORE, I3M_CORE):
            for order in (0, 1, 2):
                nq += len(domain_classes(spec, order, kc))
        log(f"micro: I3M same(+,+) order-2 quad = {mp.nstr(val, args.dps)}  wall={dt:.1f}s")
        log(f"unique quads total = {nq}; projected quad wall ~= {nq * dt:.0f}s "
            f"(order-2 is the heaviest => upper bound); RSS={rss_gb():.2f}GB")
        return 0

    # --- full gate ---
    log("\n## omega-cores via generic engine")
    if args.jobs > 1:
        jj = compute_cores_parallel(_SPECS, args.dps, args.panels,
                                    args.maxdegree, args.jobs, log=log)
        j1, j3 = jj["I1M"], jj["I3M"]
    else:
        j1 = compute_core(I1M_CORE, args.dps, brk, args.maxdegree, log=log)
        j3 = compute_core(I3M_CORE, args.dps, brk, args.maxdegree, log=log)

    need_d = max(4, args.dps - 5)  # honest reduced-precision bar
    ok = True

    log("\n## eps-layer moments")
    for tag, j in (("j1", j1), ("j3", j3)):
        for k in (0, 1, 2):
            log(f"  {tag}^({k}) = {mp.nstr(j[k], args.dps)}")
    d_j10 = -mp.log10(abs(j1[0] - mp.mpf(1) / 30) / (mp.mpf(1) / 30))
    log(f"  CHECK j1^(0) vs 1/30: {mp.nstr(d_j10, 4)} digits (need > {need_d})")
    ok &= (d_j10 > need_d)

    I1 = asm.prefactor_I1() * asm.core_to_L(j1)
    I3 = asm.prefactor_I3() * asm.core_to_L(j3)
    I2 = asm.ibp_I2(I1, I3)
    log("\n## assembled Laurent series")
    log("  I1 =", I1)
    log("  I3 =", I3)
    log("  I2 =", I2)

    target = 1 / (15 * (8 * mp.pi) ** 4)
    d_I1 = -mp.log10(abs(I1.coeff(-1) / target - 1))
    log(f"\n  CHECK I1[eps^-1] vs 1/(15(8pi)^4): {mp.nstr(d_I1, 4)} digits (need > {need_d})")
    ok &= (d_I1 > need_d)

    for n in (-4, -3):
        rel = abs(I2.coeff(n) / I2.coeff(-2))
        log(f"  CHECK I2[eps^{n}] cancellation: rel = {mp.nstr(rel, 4)} "
            f"(need < 1e-{need_d})")
        ok &= (rel < mp.mpf(10) ** (-need_d))

    cM = -mp.mpf(6) * (8 * mp.pi) ** 4 / 5 * I2.coeff(-2)
    dev = abs(cM - 1)
    d_cM = -mp.log10(dev) if dev > 0 else mp.inf
    log(f"\n## RESULT  c_M = {mp.nstr(cM, args.dps)}")
    log(f"  |c_M - 1| = {mp.nstr(dev, 4)}  -> {mp.nstr(d_cM, 4)} digits vs the "
        f"reference c_M = 1  (need > {need_d})")
    ok &= (d_cM > need_d)

    try:
        # mpmath pslq needs >= 53-bit working prec; the DATA precision stays
        # capped by the quadrature dps (tol set from args.dps, not workdps)
        with mp.workdps(max(16, args.dps)):
            rel = mp.pslq([+cM, mp.mpf(1)], tol=mp.mpf(10) ** (-(args.dps - 3)),
                          maxcoeff=10 ** 4)
        log(f"  pslq([c_M, 1]) -> {rel}")
        ok &= (rel == [1, -1])
    except Exception as ex:
        log(f"  pslq failed: {ex}")
        ok = False

    wall = time.time() - T0
    log(f"\n## wall = {wall:.1f}s  RSS = {rss_gb():.2f}GB")
    log(f"## GATE: {'PASS' if ok else 'FAIL'} at dps={args.dps} (bar {need_d}d)")

    if args.json_out:
        ent = REGISTRY[("P", "memory", "t40-static")]
        bv = asm.boundary_vector(
            "P", "memory", "t40-static", ent["prescription"],
            layers=[{
                "k": ent["layer_k"],
                "slots": [
                    {"name": "I1M", "closure": "numeric(core)+gamma(prefactor)",
                     "eps_series": {str(n): mp.nstr(I1.coeff(n), args.dps)
                                    for n in range(-1, 2)},
                     "digits": float(d_I1), "notes": "pole gated vs 1/(15(8pi)^4)"},
                    {"name": "I3M", "closure": "numeric(core)+gamma(prefactor)",
                     "eps_series": {str(n): mp.nstr(I3.coeff(n), args.dps)
                                    for n in range(-1, 2)},
                     "digits": None, "notes": "unknown - measuring (no independent gate yet)"},
                    {"name": "I2M", "closure": "ibp(I1M,I3M)",
                     "eps_series": {str(n): mp.nstr(I2.coeff(n), args.dps)
                                    for n in range(-4, 0)},
                     "digits": float(d_cM),
                     "notes": "c_M slot; checked against the reference c_M=1"},
                ]}],
            meta={"dps": args.dps, "maxdegree": args.maxdegree,
                  "panels": [str(b) for b in brk], "wall_s": round(wall, 1),
                  "rss_gb": round(rss_gb(), 2), "gate": "PASS" if ok else "FAIL",
                  "bar_digits": need_d},
        )
        asm.emit(args.json_out, bv)
        log(f"## boundary vector -> {args.json_out}")

    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
