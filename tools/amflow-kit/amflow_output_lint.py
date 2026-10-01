#!/usr/bin/env python3
"""amflow_output_lint.py — sanity validator for amflow_cli solve_integrals output.

Usage:
  amflow_output_lint.py [options] out_*.json ...

Checks per file:
  (a) nonzero coefficient count + fraction of integrals with any nonzero coef
      (--min-nonzero-frac, default 0.5; the all-important hard FAIL is frac==0)
  (b) Laurent depth: every NONZERO integral must carry >=2 nonzero eps-orders
      unless --allow-shallow (catches truncated/single-coef outputs)
  (c) --class vacuum|kinematic:
        vacuum    -> at least one integral must lead at eps-order < 0
                     (multi-loop vacuum families always have UV poles);
                     with --in JOB.json, every NONZERO dotted integral of a
                     multi-loop family must individually lead < 0
        kinematic -> no lead-sign constraint
  (d) all-zero output (every coefficient zero / no coefficients) = hard FAIL
  (e) every value parses (arb ball "[mid +/- rad]", "[+/- rad]", plain number,
      or [re,im] pair) and midpoints are finite
  (f) exit nonzero iff any FAIL

Calibrated on the shipped smoke fixture (tools/fixtures/amflow_smoke/
out_vac2Bprobe.json: 5 ints, leads -2..-1, one *genuinely* zero integral
[-1,2,2] — zero integrals are WARN, not FAIL) and on the output classes it
gates: HEALTHY = large tables with near-all-nonzero integrals and full
Laurent depth; BROKEN = all-zero tables and single-coefficient eps rows.
"""
import argparse
import json
import math
import re
import sys

BALL_RE = re.compile(r"^\[(?P<mid>[^\s\]]*)\s*\+/-\s*(?P<rad>[^\]]*)\]$")


def parse_arb(s):
    """Parse an arb-style string -> (midpoint float, ok). Accepts:
    "[mid +/- rad]", "[+/- rad]" (zero midpoint), plain decimal/rational."""
    s = str(s).strip()
    if s.startswith("["):
        m = BALL_RE.match(s)
        if not m:
            return float("nan"), False
        mid = m.group("mid").strip()
        if mid in ("", "+", "-"):  # "[+/- 1e-20]" — zero-midpoint ball
            return 0.0, True
        s = mid
    try:
        if "/" in s:
            num, den = s.split("/", 1)
            v = float(num) / float(den)
        else:
            v = float(s)
    except (ValueError, ZeroDivisionError):
        return float("nan"), False
    return v, math.isfinite(v)


def coef_value(c):
    """-> (re_mid, im_mid, parse_ok). value is {re,im} dict or [re,im] list."""
    v = c.get("value")
    if isinstance(v, dict):
        re_s, im_s = v.get("re", "0"), v.get("im", "0")
    elif isinstance(v, (list, tuple)) and len(v) >= 2:
        re_s, im_s = v[0], v[1]
    else:
        return float("nan"), float("nan"), False
    rm, rok = parse_arb(re_s)
    im, iok = parse_arb(im_s)
    return rm, im, (rok and iok)


def is_dotted(indices):
    return any(a >= 2 for a in indices if isinstance(a, (int, float)))


def lint_file(path, cls, allow_shallow, min_nonzero_frac, job=None):
    """-> (status_str, rows, messages). status in PASS/FAIL."""
    msgs, rows = [], []
    try:
        with open(path) as fh:
            data = json.load(fh)
    except (OSError, json.JSONDecodeError) as e:
        return "FAIL", rows, [f"unreadable/bad JSON: {e}"]

    result = data.get("result", data)
    if not isinstance(result, list) or not result:
        return "FAIL", rows, ["no 'result' integral list"]

    job_ints, n_loops = None, None
    if job is not None:
        job_ints = [i.get("indices", []) for i in job.get("integrals", [])]
        n_loops = len(job.get("family", {}).get("loops", []))
        if len(job_ints) != len(result):
            msgs.append(f"WARN --in integral count {len(job_ints)} != output {len(result)}")
            job_ints = None

    fail = False
    n_nonzero_int = 0
    tot_coef = tot_nonzero_coef = 0

    for k, integ in enumerate(result):
        coefs = integ.get("coefficients", [])
        nz_orders, bad_parse = [], 0
        for c in coefs:
            rm, im, ok = coef_value(c)
            if not ok:
                bad_parse += 1
                continue
            if rm != 0.0 or im != 0.0:
                nz_orders.append(c.get("order"))
        tot_coef += len(coefs)
        tot_nonzero_coef += len(nz_orders)
        lead = min(nz_orders) if nz_orders else None
        depth = len(nz_orders)
        flags = []
        if bad_parse:
            flags.append(f"FAIL:{bad_parse}-unparseable")        # (e)
            fail = True
        if not coefs:
            flags.append("FAIL:no-coefficients")
            fail = True
        if depth == 0:
            flags.append("ZERO")                                  # warn only
        else:
            n_nonzero_int += 1
            if depth < 2 and not allow_shallow:                   # (b)
                flags.append("FAIL:shallow-laurent")
                fail = True
            if (cls == "vacuum" and job_ints is not None and n_loops
                    and n_loops >= 2 and is_dotted(job_ints[k]) and lead >= 0):
                flags.append("FAIL:dotted-vacuum-lead>=0")        # (c) per-int
                fail = True
        rows.append((k, len(coefs), depth, lead, ",".join(flags) or "ok"))

    frac = n_nonzero_int / len(result)
    if n_nonzero_int == 0:                                        # (d)
        msgs.append("FAIL all-zero output: every coefficient of every integral is zero")
        fail = True
    elif frac < min_nonzero_frac:                                 # (a)
        msgs.append(f"FAIL nonzero-integral fraction {frac:.2f} < {min_nonzero_frac}")
        fail = True
    if cls == "vacuum" and n_nonzero_int > 0:                     # (c) file-level
        leads = [r[3] for r in rows if r[3] is not None]
        if not leads or min(leads) >= 0:
            msgs.append("FAIL class=vacuum: no integral leads at eps-order < 0 "
                        "(multi-loop vacuum must have UV poles)")
            fail = True
    msgs.append(f"integrals {len(result)} ({n_nonzero_int} nonzero, frac {frac:.2f}); "
                f"coefficients {tot_nonzero_coef}/{tot_coef} nonzero")
    return ("FAIL" if fail else "PASS"), rows, msgs


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("files", nargs="+", help="out_*.json from amflow_cli solve_integrals")
    ap.add_argument("--class", dest="cls", choices=["vacuum", "kinematic"],
                    default="kinematic", help="integral class rules (default kinematic)")
    ap.add_argument("--allow-shallow", action="store_true",
                    help="permit single-eps-order (depth-1) integrals")
    ap.add_argument("--min-nonzero-frac", type=float, default=0.5,
                    help="min fraction of integrals with any nonzero coef (default 0.5)")
    ap.add_argument("--in", dest="job", default=None,
                    help="optional matching in_*.json job spec (enables per-integral "
                         "dotted-vacuum lead check)")
    ap.add_argument("-q", "--quiet", action="store_true", help="summary lines only")
    args = ap.parse_args(argv)

    job = None
    if args.job:
        with open(args.job) as fh:
            job = json.load(fh)

    any_fail = False
    for path in args.files:
        status, rows, msgs = lint_file(path, args.cls, args.allow_shallow,
                                       args.min_nonzero_frac, job)
        any_fail |= (status == "FAIL")
        print(f"== {path}  [{status}]  class={args.cls}")
        for m in msgs:
            print(f"   {m}")
        if not args.quiet and rows:
            print(f"   {'int':>4} {'ncoef':>6} {'depth':>6} {'lead':>5}  flags")
            for k, nc, dp, ld, fl in rows:
                print(f"   {k:>4} {nc:>6} {dp:>6} {str(ld):>5}  {fl}")
    return 1 if any_fail else 0


if __name__ == "__main__":
    sys.exit(main())
