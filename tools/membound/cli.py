#!/usr/bin/env python3
r"""membound command line: eps-layer moments of a USER-SUPPLIED omega-core.

    python3 -m membound --spec core.json [--dps 11] [--maxdegree 4]
                        [--panels 0,0.5,2,8,inf] [--jobs N]
                        [--order-max 2] [--json-out out.json]

(run from the directory that contains the `membound/` package, e.g. tools/;
`python3 cli.py --spec ...` from inside the package directory is equivalent).

core.json (schema membound.core_spec.v1, see GUIDE.md INPUTS):
    {
      "name":   "mycore",
      "n_freq": 2,                       # 2 or 3 independent frequencies;
                                         # index n_freq+1 is the coupler
                                         # w3 = w1+w2  (or w4 = w1+w2+w3)
      "ret":    [[idx, p, q, flag], ...],# (0+ - i s w_idx)^(p - q eps),
                                         # flag in ret | adv | fey
      "klines": [idx, ...],              # |w_idx|^{-eps} K_eps(|w_idx|) lines
      "poly":   {"idx": n, ...},         # extra integer powers w_idx^n
      "proj":   "re" | "im"              # real-even or odd projection
    }

Output: the Laurent coefficients J^(0), J^(1), J^(2) of
    J(eps) = (2pi)^{-n_freq} Int d^n w (integrand expanded in eps)
           = J^(0) + J^(1) eps + J^(2) eps^2 + O(eps^3),
printed, and with --json-out written as a membound.core_moments.v1 record.
Size a run first with a low --dps / --maxdegree; cost grows steeply with both
(the c_M gate's two cores at dps 30 / maxdegree 6 took ~14 min over 12
processes).
"""
import argparse
import json
import os
import sys
import time

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from membound.core import OmegaCoreSpec, compute_core, KCache
    from membound.farm import compute_cores_parallel, parse_panels
else:
    from .core import OmegaCoreSpec, compute_core, KCache
    from .farm import compute_cores_parallel, parse_panels

import mpmath as mp

CONTRACT = "membound.core_moments.v1"


def load_spec(path):
    with open(path) as fh:
        d = json.load(fh)
    if d.get("contract") == "membound.core_spec.v1":
        d = {k: v for k, v in d.items() if k != "contract"}
    return OmegaCoreSpec.from_dict(d)


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="python3 -m membound",
        description="eps-layer moments J^(0..2) of a user-supplied "
                    "memory omega-core (see GUIDE.md INPUTS for the JSON "
                    "schema).")
    ap.add_argument("--spec", required=True,
                    help="core spec JSON file (membound.core_spec.v1)")
    ap.add_argument("--dps", type=int, default=11,
                    help="working decimal digits (default 11)")
    ap.add_argument("--maxdegree", type=int, default=4,
                    help="tanh-sinh maxdegree per axis (default 4)")
    ap.add_argument("--panels", default="0,0.5,2,8,inf",
                    help="magnitude-axis panel breakpoints (default "
                         "0,0.5,2,8,inf)")
    ap.add_argument("--jobs", type=int, default=None,
                    help="worker processes (default min(12, cpu count); "
                         "1 = serial library path)")
    ap.add_argument("--order-max", type=int, default=2, choices=(0, 1, 2),
                    help="highest eps order to compute (default 2)")
    ap.add_argument("--json-out", default=None,
                    help="write a membound.core_moments.v1 JSON record here")
    ap.add_argument("--quiet", action="store_true",
                    help="suppress per-quadrature progress lines")
    args = ap.parse_args(argv)
    if args.jobs is None:
        args.jobs = max(1, min(12, os.cpu_count() or 1))
    if args.dps < 4:
        ap.error("--dps must be >= 4")

    t0 = time.time()
    mp.mp.dps = args.dps
    try:
        spec = load_spec(args.spec)
        toks, brk = parse_panels(args.panels)
    except (OSError, ValueError, TypeError, KeyError) as ex:
        print(f"membound: error: {ex}", file=sys.stderr)
        return 2
    orders = tuple(range(args.order_max + 1))
    log = (lambda *a, **k: None) if args.quiet else print

    print(f"# membound  core={spec.name}  n_freq={spec.n_freq}  "
          f"layer_k={spec.layer_k}  proj={spec.proj}")
    print(f"# dps={args.dps} maxdegree={args.maxdegree} panels={toks} "
          f"jobs={args.jobs} orders={list(orders)}")
    kc = KCache(args.dps)
    log(f"# kappa2 FD: workdps={kc.W} h=1e-{kc.W // 6} (dps-scaled)")

    if args.jobs > 1:
        j = compute_cores_parallel({"core": spec}, args.dps, toks,
                                   args.maxdegree, args.jobs, log=log,
                                   orders=orders)["core"]
    else:
        j = compute_core(spec, args.dps, brk, args.maxdegree, log=log, kc=kc,
                         orders=orders)

    wall = time.time() - t0
    print("\n## eps-layer moments (Laurent coefficients of J(eps), "
          "measure (2pi)^-n included)")
    terms = []
    for k in orders:
        print(f"  J^({k}) = {mp.nstr(j[k], args.dps)}")
        terms.append(f"({mp.nstr(j[k], args.dps)})" + (f"*eps^{k}" if k else ""))
    print(f"  J(eps) = {' + '.join(terms)} + O(eps^{orders[-1] + 1})")
    print(f"## wall = {wall:.1f}s")

    if args.json_out:
        rec = {
            "contract": CONTRACT,
            "spec": spec.to_dict(),
            "layer_k": spec.layer_k,
            "layer_basis": "|omega|^{-k*eps}",
            "J": {str(k): mp.nstr(j[k], args.dps) for k in orders},
            "meta": {"dps": args.dps, "maxdegree": args.maxdegree,
                     "panels": toks, "jobs": args.jobs,
                     "wall_s": round(wall, 1),
                     "digits": None,
                     "notes": "unknown - measuring: no independent gate for "
                              "a user core; rerun at higher --dps/--maxdegree "
                              "and difference the two runs to measure digits"},
        }
        with open(args.json_out, "w") as fh:
            json.dump(rec, fh, indent=1)
        print(f"## core moments -> {args.json_out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
