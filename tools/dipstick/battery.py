#!/usr/bin/env python3
"""DIPSTICK battery — reproduction legs vs the reference fixtures.

Legs (each PASS/FAIL, receipts to --receipts):
  order_member   regress_order.jl against tools/dipstick/pf_rank.jl  -> 7/7 PASS +
                 the reference walls: sunrise order 2 / hr 7, 3L-banana order 3 / hr 15.
  order_shim     regress_order.jl against tools/pf_rank.jl (the forever-shim) ->
                 identical results (timings stripped) to order_member.
  count          `dipstick count 6` -> crit36.ms byte-identical to the reference
                 receipt; msolve (capped) -> parametrization degree 26.
  cycletype      `dipstick cycletype crit36.ms` -> per-prime .ms byte-identical
                 to the reference 6-prime files; cycle types equal to the types
                 independently re-derived (flint) from the REFERENCE .out eliminants; n=26.
  nonres         `dipstick nonres` -> stdout equal to the reference nonres.log
                 (43 exact-verified facets + the reference pairing verdict).
  audit          `dipstick audit` AND `dipstick audit --fast` on the vendored
                 toy fixture fixtures/audit_toy.ms (4 generic lines): both routes
                 must report crit = 3 (= the C(n-1,2) bounded-region cage) with
                 identical census/degR/arr_val.
  flatchi        `dipstick flatchi fixtures/audit_toy.ms` (affine route) -> chi =
                 q^2-4q+6, bounded 3 (equals the audit leg's crit count —
                 Varchenko), regions 11; plus `dipstick flatchi --quick`
                 (projective calibration cage) -> chi_top 42.
  regions        regions/selftest.jl: the expansion-by-regions completeness
                 certifier (regions/regions.jl) on its three shipped inputs in
                 ONE julia process — bubble_largeQ -> TILES OK / CHI ADDITIVITY
                 OK / VERDICT PASS; largemass_toy_complete -> TILES OK / VERDICT
                 PASS; largemass_toy_incomplete -> VERDICT FAIL with the missing
                 facet v=[1,1,1] NAMED.  julia + Oscar + JSON; without them the
                 leg SKIPs by name.

--tool-dir lets a MUTATION CONTROL point the same comparator at a mutated copy:
a mutated member must FAIL its leg.  Run: python3 battery.py [--legs a,b] ...

A leg whose external engine is absent (julia + Oscar for the order and regions
legs, msolve for count/cycletype, scipy for nonres) SKIPs with a named reason
and does not fail the battery; a wrong result with the engine present still
FAILs.
"""
import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
PRIMES = [1073741783, 1073741789, 1073741827, 1073741831, 1073741833, 1073741839]
_JP = os.environ.get("DIPSTICK_JULIA_PROJECT", "")
# DIPSTICK_JULIA overrides the julia invocation (e.g. "julia +1.10" to pin a
# juliaup channel); the default is whatever julia is on PATH.
JULIA = os.environ.get("DIPSTICK_JULIA", "julia").split() + \
    ([f"--project={_JP}"] if _JP else [])
MSOLVE = os.environ.get("DIPSTICK_MSOLVE",
                        os.path.join(HERE, "run_msolve_capped.sh"))


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def engine_missing(leg):
    """Cold-clone contract: a leg whose external engine is absent SKIPs with a
    named reason (by design, not breakage — see INSTALL.md); a leg whose engine
    is present but whose result is wrong still FAILs."""
    import importlib.util
    import shutil
    if leg in ("order_member", "order_shim", "regions"):
        if shutil.which(JULIA[0]) is None:
            return (f"julia not on PATH ('{JULIA[0]}'; the order and regions legs "
                    f"need julia + Oscar)")
    elif leg in ("count", "cycletype"):
        if not os.environ.get("DIPSTICK_MSOLVE") and \
                shutil.which(os.environ.get("MSOLVE_BIN", "msolve")) is None:
            return "msolve not on PATH (see toolkit/external/TOOLS.md; MSOLVE_BIN overrides)"
    elif leg == "nonres":
        if importlib.util.find_spec("scipy") is None:
            return "scipy not installed (nonres needs scipy.spatial.ConvexHull)"
    elif leg == "audit":
        for mod in ("sympy", "numpy"):
            if importlib.util.find_spec(mod) is None:
                return f"{mod} not installed (the audit members need sympy + numpy)"
    elif leg == "flatchi":
        if importlib.util.find_spec("sympy") is None:
            return "sympy not installed (the flatchi members need sympy)"
    return None


def cycle_types_from_out(outfile):
    """Factor the eliminant in a msolve .out at its prime -> sorted degree multiset.
    Same parse+factor route as the cycletype member (independent instance here,
    applied to the REFERENCE outputs)."""
    import flint
    s = open(outfile).read().replace("\n", " ")
    p = int(re.match(r"\[0, \[(\d+),", s).group(1))
    m = re.search(r"\[1, \[\[(\d+), \[([0-9, -]+)\]\]", s)
    coeffs = [int(x) % p for x in m.group(2).split(",")]
    fac = flint.nmod_poly(coeffs, p).factor()
    degs = sorted(int(g.degree()) for g, mult in fac[1] for _ in range(int(mult)))
    return p, int(m.group(1)), degs


def strip_timings(txt):
    # drop wall-clock tokens AND the "# tool under test: <path>" banner — the
    # target path differs BY DESIGN between the member and shim routes; every
    # result line must be identical.
    txt = re.sub(r"^# tool under test: .*$", "", txt, flags=re.M)
    return re.sub(r"\d+\.\d+s", "T", txt)


def leg_order(tool_dir, target, work, log):
    t0 = time.time()
    r = subprocess.run(JULIA + [f"{tool_dir}/regress_order.jl", target],
                       capture_output=True, text=True, cwd=work, timeout=3600)
    open(log, "w").write(r.stdout + r.stderr)
    if r.returncode != 0:
        for pat, why in (
                ("Package Oscar not found",
                 "Oscar not in the active Julia project (set DIPSTICK_JULIA_PROJECT)"),
                ("Leviathan.jl not found",
                 "Leviathan.jl not found (set DIPSTICK_LEVIATHAN)")):
            if pat in r.stdout + r.stderr:
                return {"pass": None, "skip": why,
                        "wall_s": round(time.time() - t0, 1)}
    ok = (r.returncode == 0 and "RESULT: 7/7 PASS" in r.stdout
          and re.search(r"c1 sunrise v1-probe: order=2 .*hr=7", r.stdout)
          and re.search(r"c7 3L-banana v1:\s+order=3 .*hr=15", r.stdout))
    return {"pass": bool(ok), "wall_s": round(time.time() - t0, 1),
            "normalized": strip_timings(r.stdout)}


def leg_count(tool_dir, fixtures, work, log):
    t0 = time.time()
    r = subprocess.run([sys.executable, f"{tool_dir}/dipstick", "count", "6"],
                       capture_output=True, text=True, cwd=work, timeout=600)
    ms_ok = os.path.exists(f"{work}/crit36.ms") and \
        sha(f"{work}/crit36.ms") == sha(f"{fixtures}/crit36.ms")
    t1 = time.time()
    env = dict(os.environ, MSOLVE_CAP_GB="30")
    r2 = subprocess.run([MSOLVE, "-f", "crit36.ms", "-o", "crit36.out", "-t", "8"],
                        capture_output=True, text=True, cwd=work, env=env, timeout=1800)
    msolve_wall = round(time.time() - t1, 2)
    deg = None
    if os.path.exists(f"{work}/crit36.out"):
        s = open(f"{work}/crit36.out").read().replace("\n", " ")
        m = re.search(r"\[1, \[\[(\d+), \[", s)
        deg = int(m.group(1)) if m else None
    out_sha_match = os.path.exists(f"{work}/crit36.out") and \
        sha(f"{work}/crit36.out") == sha(f"{fixtures}/crit36.out")
    open(log, "w").write(r.stdout + r.stderr + "\n--- msolve ---\n" + r2.stdout + r2.stderr)
    return {"pass": bool(ms_ok and deg == 26), "ms_byte_identical": bool(ms_ok),
            "degree": deg, "expect_degree": 26, "out_sha_match": bool(out_sha_match),
            "msolve_wall_s": msolve_wall, "wall_s": round(time.time() - t0, 1)}


def leg_cycletype(tool_dir, fixtures, work, log):
    t0 = time.time()
    import shutil
    shutil.copyfile(f"{fixtures}/crit36.ms", f"{work}/crit36.ms")
    r = subprocess.run([sys.executable, f"{tool_dir}/dipstick", "cycletype", "crit36.ms"],
                       capture_output=True, text=True, cwd=work, timeout=3600)
    open(log, "w").write(r.stdout + r.stderr)
    ms_ident = all(os.path.exists(f"{work}/gc_{p}.ms")
                   and sha(f"{work}/gc_{p}.ms") == sha(f"{fixtures}/crit36_{p}.ms")
                   for p in PRIMES)
    got = {int(m.group(1)): json.loads(m.group(2)) for m in
           re.finditer(r"p=(\d+): cycle type (\[[0-9, ]+\])", r.stdout)}
    ref, n_ref = {}, set()
    for p in PRIMES:
        pp, nn, degs = cycle_types_from_out(f"{fixtures}/crit36_{p}.out")
        assert pp == p
        ref[p] = degs
        n_ref.add(nn)
    types_ok = got == ref
    n_ok = n_ref == {26} and all(sum(v) == 26 for v in got.values())
    return {"pass": bool(ms_ident and types_ok and n_ok),
            "per_prime_ms_byte_identical": bool(ms_ident),
            "types_member": {str(k): v for k, v in sorted(got.items())},
            "types_from_reference_outs": {str(k): v for k, v in sorted(ref.items())},
            "n": sorted(n_ref), "wall_s": round(time.time() - t0, 1)}


def leg_nonres(tool_dir, fixtures, work, log):
    t0 = time.time()
    r = subprocess.run([sys.executable, f"{tool_dir}/dipstick", "nonres"],
                       capture_output=True, text=True, cwd=work, timeout=1800)
    open(log, "w").write(r.stdout + ("\n--- stderr ---\n" + r.stderr if r.stderr else ""))
    ref_log = open(f"{fixtures}/nonres.log").read()
    ok = r.stdout == ref_log
    return {"pass": bool(ok), "stdout_equals_reference_log": bool(ok),
            "facet_line_reproduced": "43 exact-verified distinct facets" in r.stdout,
            "wall_s": round(time.time() - t0, 1)}


def leg_audit(tool_dir, fixtures, work, log):
    """Both audit routes on the vendored toy fixture (4 generic lines in the
    plane): crit must equal 3 = C(n-1,2), the generic-arrangement cage, and the
    fast route must agree with the exact route field-for-field."""
    import ast
    import shutil
    t0 = time.time()
    shutil.copyfile(f"{fixtures}/audit_toy.ms", f"{work}/audit_toy.ms")
    outs = {}
    for name, argv in (("exact", ["audit", "audit_toy.ms"]),
                       ("fast", ["audit", "--fast", "audit_toy.ms"])):
        r = subprocess.run([sys.executable, f"{tool_dir}/dipstick"] + argv,
                           capture_output=True, text=True, cwd=work, timeout=600)
        outs[name] = {"rc": r.returncode, "stdout": r.stdout, "stderr": r.stderr}
        try:
            outs[name]["parsed"] = ast.literal_eval(r.stdout.strip().splitlines()[-1])
        except Exception:
            outs[name]["parsed"] = None
    open(log, "w").write(json.dumps({k: {kk: vv for kk, vv in v.items()
                                         if kk != "parsed"} for k, v in outs.items()},
                                    indent=1))
    ex, fa = outs["exact"]["parsed"], outs["fast"]["parsed"]
    crit_ok = bool(ex and fa and ex.get("crit") == 3 and fa.get("crit") == 3)
    agree = bool(ex and fa and all(ex.get(k) == fa.get(k)
                                   for k in ("lines", "census", "degR", "arr_val", "crit")))
    return {"pass": bool(crit_ok and agree and outs["exact"]["rc"] == 0
                         and outs["fast"]["rc"] == 0),
            "crit_exact": ex.get("crit") if ex else None,
            "crit_fast": fa.get("crit") if fa else None, "expect_crit": 3,
            "routes_agree": agree, "wall_s": round(time.time() - t0, 1)}


def leg_flatchi(tool_dir, fixtures, work, log):
    """Affine route on the toy fixture (chi = q^2-4q+6, bounded 3 = the audit
    leg's crit count, regions 11) + the quick projective calibration cage
    (chi_top 42)."""
    import shutil
    t0 = time.time()
    shutil.copyfile(f"{fixtures}/audit_toy.ms", f"{work}/audit_toy.ms")
    r1 = subprocess.run([sys.executable, f"{tool_dir}/dipstick", "flatchi",
                         "audit_toy.ms"],
                        capture_output=True, text=True, cwd=work, timeout=600)
    r2 = subprocess.run([sys.executable, f"{tool_dir}/dipstick", "flatchi",
                         "--quick"],
                        capture_output=True, text=True, cwd=work, timeout=600)
    open(log, "w").write(r1.stdout + r1.stderr + "\n--- quick cal ---\n"
                         + r2.stdout + r2.stderr)
    aff_chi = "chi_A(q) = q**2 - 4*q + 6" in r1.stdout
    aff_counts = "bounded (|chi(1)|): 3  regions: 11" in r1.stdout
    cal42 = "chi_top = 42 (expect 42)" in r2.stdout
    return {"pass": bool(r1.returncode == 0 and r2.returncode == 0 and aff_chi
                         and aff_counts and cal42),
            "affine_chi_ok": bool(aff_chi), "affine_counts_ok": bool(aff_counts),
            "projective_cal_42": bool(cal42), "wall_s": round(time.time() - t0, 1)}


def leg_regions(tool_dir, work, log):
    """The regions member's own control battery (regions/selftest.jl): three
    shipped inputs in one julia process.  selftest.jl itself SKIPs by name
    (exit 0) when neither the active environment, DIPSTICK_JULIA_PROJECT nor
    upgrades/leviathan provides Oscar + JSON; that is relayed as a SKIP here,
    not a PASS."""
    t0 = time.time()
    r = subprocess.run(JULIA + ["--startup-file=no",
                                f"{tool_dir}/regions/selftest.jl"],
                       capture_output=True, text=True, cwd=work, timeout=3600)
    out = r.stdout + r.stderr
    open(log, "w").write(out)
    if "SKIP (named engine gate)" in out:
        return {"pass": None,
                "skip": "Oscar + JSON not available to julia (regions/selftest.jl "
                        "printed the one-time instantiate command; see the log)",
                "wall_s": round(time.time() - t0, 1)}
    if "Leviathan.jl not found" in out:
        return {"pass": None, "skip": "Leviathan.jl not found (set DIPSTICK_LEVIATHAN)",
                "wall_s": round(time.time() - t0, 1)}
    controls = re.findall(r"^control ok: (\S+)", out, flags=re.M)
    ok = (r.returncode == 0 and "dipstick regions selftest PASS" in out
          and len(controls) == 3)
    return {"pass": bool(ok), "controls_ok": controls,
            "fail_closed": "ENGINE GATE (fail-closed)" in out,
            "wall_s": round(time.time() - t0, 1)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tool-dir", default=HERE)
    ap.add_argument("--fixtures", default=os.path.join(HERE, "fixtures"))
    ap.add_argument("--receipts", default=os.path.join(HERE, "fixtures"))
    ap.add_argument("--legs",
                    default="order_member,order_shim,count,cycletype,nonres,audit,flatchi,regions")
    ap.add_argument("--tag", default="")
    ap.add_argument("--summary-out", default=None,
                    help="summary json path (default: <tag>_SUMMARY.json in the run cwd; "
                         "--receipts is the read-only fixture root)")
    args = ap.parse_args()
    tag = args.tag or "battery"
    res = {}
    for leg in args.legs.split(","):
        why = engine_missing(leg)
        if why:
            res[leg] = {"pass": None, "skip": why, "wall_s": 0.0}
            print(f"{leg}: SKIP ({why})", flush=True)
            continue
        work = os.path.join(args.receipts, f"work_{tag}_{leg}")
        os.makedirs(work, exist_ok=True)
        log = os.path.join(args.receipts, f"{tag}_{leg}.log")
        if leg == "order_member":
            res[leg] = leg_order(args.tool_dir, f"{args.tool_dir}/pf_rank.jl", work, log)
        elif leg == "order_shim":
            # the forever-shim lives one level up from the tool dir (tools/)
            shim = os.path.join(os.path.dirname(os.path.abspath(args.tool_dir)),
                                "pf_rank.jl")
            res[leg] = leg_order(args.tool_dir, shim, work, log)
        elif leg == "count":
            res[leg] = leg_count(args.tool_dir, args.fixtures, work, log)
        elif leg == "cycletype":
            res[leg] = leg_cycletype(args.tool_dir, args.fixtures, work, log)
        elif leg == "nonres":
            res[leg] = leg_nonres(args.tool_dir, args.fixtures, work, log)
        elif leg == "audit":
            res[leg] = leg_audit(args.tool_dir, args.fixtures, work, log)
        elif leg == "flatchi":
            res[leg] = leg_flatchi(args.tool_dir, args.fixtures, work, log)
        elif leg == "regions":
            res[leg] = leg_regions(args.tool_dir, work, log)
        else:
            raise SystemExit(f"unknown leg {leg}")
        if res[leg].get("skip"):
            print(f"{leg}: SKIP ({res[leg]['skip']})", flush=True)
        else:
            print(f"{leg}: {'PASS' if res[leg]['pass'] else 'FAIL'} "
                  f"({res[leg]['wall_s']}s)", flush=True)
    if "normalized" in res.get("order_member", {}) and \
            "normalized" in res.get("order_shim", {}):
        ident = res["order_member"]["normalized"] == res["order_shim"]["normalized"]
        res["order_both_routes_identical"] = bool(ident)
        print(f"order both-routes identical (timings stripped): {ident}", flush=True)
    for k in ("order_member", "order_shim"):
        if k in res:
            res[k].pop("normalized", None)
    outp = args.summary_out or os.path.join(os.getcwd(), f"{tag}_SUMMARY.json")
    json.dump(res, open(outp, "w"), indent=1)
    print("summary ->", outp)
    bad = [k for k, v in res.items()
           if isinstance(v, dict) and "skip" not in v and not v.get("pass")]
    sys.exit(1 if bad or res.get("order_both_routes_identical") is False else 0)


if __name__ == "__main__":
    main()
