#!/usr/bin/env python3
"""tools/holonomic battery — engine-free reference check + Sage engine legs.

Default legs:
  reference  (no engine; runs everywhere): recompute ln(3/2) independently in
             pure Python integer arithmetic (2*atanh(1/5) with a rigorous
             tail bound) and verify that the recorded worked-example
             enclosures in reference/ln32_reference.json (a) CONTAIN it and
             (b) are internally consistent (certified radius <= printed
             radius, certified radius <= requested eps, recorded error bound
             <= twice the certified radius).
  smoke      (needs SageMath + ore_algebra): smoke_test.py under
             `sage -python` — live certified transport must contain ln(3/2)
             with radius <= 1e-38, plus the detour law, the IC-convention
             probe, and the detoured-VALUE gate (transport of ln(1+x) from
             0 to -3/2 must return a COMPLEX ball containing ln(1/2) + i*pi
             — the imaginary part must survive).

--full adds:
  ladder6    the rank-6 stage of ladder/rank_ladder.py: build the synthetic
             rank-6 hypergeometric operator, one ordinary transport at eps
             1e-40, then the validation stage (roundtrip ball must contain
             the identity; the transported known solution must contain the
             direct-series value).

Engine-absent legs SKIP by name and do not fail the battery: the smoke and
ladder legs need SageMath on PATH (install from https://www.sagemath.org or
conda-forge `sage`; HOLONOMIC_SAGE overrides the executable) with ore_algebra
installed into it (`sage -pip install ore_algebra`). A present engine with a
wrong result still FAILs. Exit 0 = no leg failed (skips do not fail).
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
SAGE = os.environ.get("HOLONOMIC_SAGE", "sage")


def ln32_independent(digits=70):
    """ln(3/2) = 2*atanh(1/5) as an exact Fraction with a rigorous error
    bound, in pure Python integer arithmetic (no float, no library):
    2*atanh(1/5) = 2 * sum_{k>=0} (1/5)^(2k+1) / (2k+1), summed in units of
    10^-digits with floor division. Each floor loses < 1 unit; the first
    dropped term is < 1 unit and the tail after it is geometric with ratio
    1/25, so the total error is < (terms + 2) units."""
    S = 10 ** digits
    total, k = 0, 0
    while True:
        t = (2 * S) // ((5 ** (2 * k + 1)) * (2 * k + 1))
        if t == 0:
            break
        total += t
        k += 1
    return Fraction(total, S), Fraction(k + 2, S)


def leg_reference():
    t0 = time.time()
    ref = json.load(open(os.path.join(HERE, "reference",
                                      "ln32_reference.json")))
    ln32, err = ln32_independent()
    out = {"ln32_digits_used": 70, "records": []}
    ok = True
    for rec in ref["records"]:
        mid = Fraction(rec["midpoint"])
        pr = Fraction(rec["printed_radius"])
        cr = Fraction(rec["certified_radius"])
        eps = Fraction(rec["eps"])
        contains = abs(mid - ln32) <= pr + err
        rad_consistent = cr <= pr
        eps_met = cr <= eps
        checks = {"contains_independent_ln32": contains,
                  "certified_radius_le_printed": rad_consistent,
                  "certified_radius_le_eps": eps_met}
        if "abs_err_vs_ln32_upper" in rec:
            eb = Fraction(rec["abs_err_vs_ln32_upper"])
            # containment => |ball - ln32| <= 2*rad (+ the 300-bit
            # comparison ball's own negligible radius)
            checks["recorded_err_bound_le_2rad"] = \
                eb <= 2 * cr + Fraction(1, 10 ** 60)
        rec_ok = all(checks.values())
        ok = ok and rec_ok
        out["records"].append({"source": rec["source"], "checks":
                               {k: bool(v) for k, v in checks.items()},
                               "pass": bool(rec_ok)})
    out["pass"] = bool(ok)
    out["wall_s"] = round(time.time() - t0, 2)
    return out


_SAGE_PROBE = []          # memo: probe sage once, not once per engine leg


def sage_missing():
    """Named-skip contract: the engine legs need SageMath + ore_algebra."""
    if _SAGE_PROBE:
        return _SAGE_PROBE[0]
    _SAGE_PROBE.append(_sage_missing_uncached())
    return _SAGE_PROBE[0]


def _sage_missing_uncached():
    if shutil.which(SAGE) is None:
        return (f"SageMath not on PATH ('{SAGE}'; install from "
                f"https://www.sagemath.org or conda-forge, then "
                f"`sage -pip install ore_algebra`; HOLONOMIC_SAGE overrides "
                f"the executable)")
    r = subprocess.run([SAGE, "-python", "-c",
                        "import sage.all; import ore_algebra"],
                       capture_output=True, text=True, timeout=600)
    if r.returncode != 0:
        tail = (r.stderr.strip().splitlines() or ["no stderr"])[-1]
        return (f"ore_algebra not importable under this sage "
                f"(`sage -pip install ore_algebra`; probe said: {tail[:120]})")
    return None


def leg_smoke(receipts):
    t0 = time.time()
    receipt = os.path.join(receipts, "holonomic_smoke_receipt.json")
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    r = subprocess.run([SAGE, "-python",
                        os.path.join(HERE, "smoke_test.py"), receipt],
                       capture_output=True, text=True, cwd=HERE, env=env,
                       timeout=1800)
    ok = r.returncode == 0 and "SMOKE PASS" in r.stdout
    return {"pass": bool(ok), "receipt": receipt,
            "tail": r.stdout.strip().splitlines()[-3:],
            "wall_s": round(time.time() - t0, 1)}


def leg_ladder6(receipts):
    t0 = time.time()
    work = tempfile.mkdtemp(prefix="holonomic_ladder_", dir=receipts)
    env = dict(os.environ, RL_RANK="6", RL_STAGES="ord40,validate",
               RL_OUT=work, PYTHONDONTWRITEBYTECODE="1")
    r = subprocess.run([SAGE, "-python",
                        os.path.join(HERE, "ladder", "rank_ladder.py")],
                       capture_output=True, text=True, cwd=work, env=env,
                       timeout=1800)
    rj = os.path.join(work, "ladder_R6.json")
    res = {"pass": False, "wall_s": None, "workdir": work}
    if r.returncode == 0 and os.path.exists(rj):
        d = json.load(open(rj))
        v = d.get("validate", {})
        res["pass"] = bool(v.get("roundtrip_contains_identity")
                           and v.get("known_solution_contains_zero")
                           and "ord40" in d.get("stages", {}))
        res["ic_convention"] = v.get("ic_convention")
        res["ord40_wall_s"] = d.get("stages", {}).get("ord40",
                                                      {}).get("wall_s")
    else:
        res["tail"] = (r.stdout + r.stderr).strip().splitlines()[-3:]
    res["wall_s"] = round(time.time() - t0, 1)
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true",
                    help="add the rank-6 ladder leg (build + transport + "
                         "validation under sage)")
    ap.add_argument("--legs", default=None,
                    help="comma list to override (reference,smoke,ladder6)")
    ap.add_argument("--receipts", default=os.getcwd(),
                    help="where the smoke receipt lands (default: run cwd)")
    args = ap.parse_args()
    legs = (args.legs.split(",") if args.legs else
            ["reference", "smoke"] + (["ladder6"] if args.full else []))
    res = {}
    for leg in legs:
        leg = leg.strip()
        if leg in ("smoke", "ladder6"):
            why = sage_missing()
            if why:
                res[leg] = {"pass": None, "skip": why, "wall_s": 0.0}
                print(f"{leg}: SKIP ({why})", flush=True)
                continue
        if leg == "reference":
            res[leg] = leg_reference()
        elif leg == "smoke":
            res[leg] = leg_smoke(args.receipts)
        elif leg == "ladder6":
            res[leg] = leg_ladder6(args.receipts)
        else:
            raise SystemExit(f"unknown leg {leg}")
        print(f"{leg}: {'PASS' if res[leg]['pass'] else 'FAIL'} "
              f"({res[leg]['wall_s']}s)", flush=True)
    outp = os.path.join(os.getcwd(), "holonomic_SUMMARY.json")
    json.dump(res, open(outp, "w"), indent=1)
    print("summary ->", outp)
    bad = [k for k, v in res.items() if "skip" not in v and not v["pass"]]
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
