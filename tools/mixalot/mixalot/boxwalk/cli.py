#!/usr/bin/env python3
"""mixalot.boxwalk CLI.  Usage (see README.md beside this file); from
tools/mixalot/ (or anywhere with tools/mixalot on sys.path):

  python3 -m mixalot.boxwalk selftest
  python3 -m mixalot.boxwalk plan    SPEC.json
  python3 -m mixalot.boxwalk produce SPEC.json OUTDIR [--nprimes N] [--batch K]
                                                      [--procs P]
  python3 -m mixalot.boxwalk fiber   SPEC.json [--ray CLASS] [--depth D]
  python3 -m mixalot.boxwalk verify  SPEC.json MANIFEST.json [--nfresh N]
                                                 [--replan] [--max-cells C]

`python3 path/to/mixalot/boxwalk/cli.py ...` works too (it puts tools/mixalot
on sys.path itself).  Farm launch: run `produce` under your machine's job
launcher (never raw nohup), e.g.
      python3 -m mixalot.boxwalk produce SPEC.json OUTDIR --nprimes 16
"""
import os, sys, json, argparse

if __package__ in (None, ''):
    # run as a plain script: tools/mixalot (three levels up) holds `mixalot`
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__)))))
from mixalot import boxwalk                                # noqa: E402
from mixalot.boxwalk import planner, driver, fiber         # noqa: E402


def main():
    ap = argparse.ArgumentParser(prog='mixalot.boxwalk')
    sub = ap.add_subparsers(dest='cmd', required=True)
    sub.add_parser('selftest')
    p1 = sub.add_parser('plan')
    p1.add_argument('spec')
    p1.add_argument('--slab-only', action='store_true')
    p2 = sub.add_parser('produce')
    p2.add_argument('spec')
    p2.add_argument('outdir')
    p2.add_argument('--nprimes', type=int, default=8)
    p2.add_argument('--batch', type=int, default=4)
    p2.add_argument('--procs', type=int, default=None)
    p2.add_argument('--slab-only', action='store_true')
    p3 = sub.add_parser('fiber')
    p3.add_argument('spec')
    p3.add_argument('--ray', default=None)
    p3.add_argument('--depth', type=int, default=80)
    p4 = sub.add_parser('verify')
    p4.add_argument('spec')
    p4.add_argument('manifest')
    p4.add_argument('--nfresh', type=int, default=2)
    p4.add_argument('--replan', action='store_true',
                    help='rebuild the plan and re-check program_sha')
    p4.add_argument('--max-cells', type=int, default=int(5e7))
    a = ap.parse_args()
    if a.cmd == 'selftest':
        from mixalot.boxwalk import selftest
        return selftest.main()
    spec = boxwalk.load_spec(a.spec)
    if a.cmd == 'plan':
        pl = planner.plan(spec, try_refill=not a.slab_only)
        print(json.dumps({'side': pl['side'], 'order': pl['order'],
                          'mem_gb_per_prime': pl['mem_gb_per_prime'],
                          'program_sha': pl['program_sha'],
                          'segments': [{'class': s['class'],
                                        'steps': s['steps'],
                                        'mech': s['mech']}
                                       for s in pl['segments']]}, indent=1))
        return 0
    if a.cmd == 'produce':
        man = driver.produce(spec, nprimes=a.nprimes, batch=a.batch,
                             procs=a.procs, outdir=a.outdir,
                             try_refill=not a.slab_only)
        return 0 if man['ok'] else 1
    if a.cmd == 'fiber':
        rec = fiber.emit_fiber_recurrence(spec, ray_class=a.ray,
                                          depth=a.depth)
        print(json.dumps(rec, indent=1))
        return 0 if rec['ok'] else 1
    if a.cmd == 'verify':
        from mixalot.boxwalk import verify
        pl = planner.plan(spec) if a.replan else None
        rep = verify.verify_manifest(spec, a.manifest, nfresh=a.nfresh,
                                     plan=pl, max_cells=a.max_cells)
        print(json.dumps(rep, indent=1))
        return 0 if rep['ok'] else 1


if __name__ == '__main__':
    sys.exit(main())
