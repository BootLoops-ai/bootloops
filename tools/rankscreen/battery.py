#!/usr/bin/env python3
"""rankscreen synthetic battery — in-tree, self-contained, no external data.

Fixtures are built deterministically in code (fixtures.py) with ground
truths that follow from their construction, so the battery runs anywhere
the tool runs. Each leg prints PASS/FAIL (flushed); a report JSON is
written to the working directory; exactly one OVERALL line; exit 0 iff
every leg passes.

Legs:
  planted_control   known-truth staircase system: the screen must return
                    CERTIFIED-SCREEN with exactly the planted rank, pivot
                    sequence, empty inconsistent set and the closure flag,
                    identically on BOTH backends, and must agree with the
                    exact eliminator run on the same rows.
  zero_eq_one       planted inconsistent system (a literal 0 = 1 row plus
                    masked shifted combinations): the screen must flag every
                    planted row, on both backends and in agreement with the
                    exact eliminator, and must never report a closure
                    candidate.
  prime_disagree    a row engineered so different primes see different
                    ranks: the verdict must be ESCALATE-TO-EXACT — including
                    at k = 3 where two of the three primes agree (escalation
                    is required on ANY dissent; there is no majority vote).
  bad_denominator   a coefficient denominator divisible by a pool prime:
                    the prime must be detected, replaced and receipted in
                    primes_replaced (never a silent wrong verdict), and a
                    prime panel that is entirely unusable must be refused.

--tool-dir points the legs at another copy of the tool (default: the
directory of this file). --mutation-controls copies the tool to a scratch
directory under the working directory, applies ONE targeted sabotage per
leg, and requires exactly the sabotaged leg to FAIL; the copy is deleted
afterward and the installed tree is never modified.

Run from a scratch directory (the battery refuses to write inside the tool
tree):  python3 tools/rankscreen/battery.py [--legs a,b] [--summary-out F]
"""
import argparse
import importlib
import json
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
LEGS = ("planted_control", "zero_eq_one", "prime_disagree", "bad_denominator")


# ------------------------------------------------------------ tool import ---
def _load_tool(tool_dir):
    """Import screen/shim/fixtures from tool_dir (fresh, path-pinned)."""
    for m in ("screen", "shim", "fixtures", "modp_rref"):
        sys.modules.pop(m, None)
    sys.path.insert(0, tool_dir)
    screen = importlib.import_module("screen")
    shim = importlib.import_module("shim")
    fixtures = importlib.import_module("fixtures")
    modp = importlib.import_module("modp_rref")
    return screen, shim, fixtures, modp


def _exact_reference(shim, rows, labels):
    """The tool's exact eliminator (gmpy2), Fraction fallback if absent."""
    try:
        return shim.gauss_rank_mpq(rows, labels)
    except ImportError:
        piv, red, redr, incon = {}, [], [], []
        for (rdict0, rv), lab in zip(rows, labels):
            rdict = dict(rdict0)
            for col, ri in list(piv.items()):
                if rdict.get(col):
                    f = rdict[col] / red[ri][col]
                    for c2, v2 in red[ri].items():
                        nv = rdict.get(c2, 0) - f * v2
                        if nv:
                            rdict[c2] = nv
                        elif c2 in rdict:
                            del rdict[c2]
                    rv = rv - f * redr[ri]
            if rdict:
                piv[min(rdict)] = len(red)
                red.append(rdict)
                redr.append(rv)
            elif rv:
                incon.append(lab)
        return len(red), len(incon), piv, red, redr, incon


def _check(cond, msg, errs):
    if not cond:
        errs.append(msg)


# ------------------------------------------------------------------- legs ---
def leg_planted_control(screen, shim, fixtures, modp):
    rows, labels, ncols, nunk, truth = fixtures.planted_control()
    errs, out = [], {}
    for backend in ("dense", "sparse"):
        v = screen.screen_rows(rows, labels, ncols, k=2, backend=backend,
                               nunk=nunk, tag=f"battery-{backend}")
        out[backend] = v
        _check(v["verdict"] == "CERTIFIED-SCREEN",
               f"{backend}: verdict {v['verdict']}", errs)
        _check(v.get("rank") == truth["rank"],
               f"{backend}: rank {v.get('rank')} != {truth['rank']}", errs)
        _check(v.get("pivot_cols") == truth["pivot_cols"],
               f"{backend}: pivot cols differ from planted", errs)
        _check(v.get("pivot_rows") == truth["pivot_rows"],
               f"{backend}: pivot rows differ from planted", errs)
        _check(v.get("incon_idx") == [],
               f"{backend}: nonempty inconsistent set", errs)
        _check(v.get("closure_candidate") is True,
               f"{backend}: closure_candidate {v.get('closure_candidate')}",
               errs)
    rank, ninc, piv, _, _, incon = _exact_reference(shim, rows, labels)
    _check(rank == truth["rank"], f"exact rank {rank}", errs)
    _check(ninc == 0 and not incon, "exact found inconsistency", errs)
    seq = [c for c, _ in sorted(piv.items(), key=lambda cv: cv[1])]
    _check(seq == truth["pivot_cols"], "exact pivot sequence differs", errs)
    if errs:
        raise AssertionError("; ".join(errs))
    return (f"rank {truth['rank']}/{nunk} recovered on both backends + exact"
            f" (dense {out['dense']['wall_secs']}s,"
            f" sparse {out['sparse']['wall_secs']}s)")


def leg_zero_eq_one(screen, shim, fixtures, modp):
    rows, labels, ncols, nunk, truth = fixtures.planted_inconsistent()
    errs = []
    for backend in ("dense", "sparse"):
        v = screen.screen_rows(rows, labels, ncols, k=2, backend=backend,
                               nunk=nunk, tag=f"battery-{backend}")
        _check(v["verdict"] == "CERTIFIED-SCREEN",
               f"{backend}: verdict {v['verdict']}", errs)
        _check(v.get("n_inconsistent") == len(truth["incon_idx"]),
               f"{backend}: {v.get('n_inconsistent')} inconsistent, planted "
               f"{len(truth['incon_idx'])}", errs)
        _check(v.get("incon_idx") == truth["incon_idx"],
               f"{backend}: flagged rows differ from planted", errs)
        _check(v.get("incon_labels") ==
               [labels[i] for i in truth["incon_idx"]],
               f"{backend}: inconsistent labels differ", errs)
        _check(v.get("closure_candidate") is False,
               f"{backend}: inconsistent system passed as closure candidate",
               errs)
        _check(v.get("rank") == truth["rank"],
               f"{backend}: rank {v.get('rank')}", errs)
    rank, ninc, _, _, _, incon = _exact_reference(shim, rows, labels)
    _check(ninc == len(truth["incon_idx"]), f"exact ninc {ninc}", errs)
    _check(incon == [labels[i] for i in truth["incon_idx"]],
           "exact inconsistent labels differ", errs)
    if errs:
        raise AssertionError("; ".join(errs))
    return (f"{len(truth['incon_idx'])} planted 0=1 rows flagged on both "
            f"backends + exact; closure candidate refused")


def leg_prime_disagree(screen, shim, fixtures, modp):
    p1, p2, p3 = modp.PRIMES_DENSE[:3]
    rows, labels, ncols, nunk, truth = fixtures.planted_prime_disagreement(p1)
    errs = []
    v2 = screen.screen_rows(rows, labels, ncols, k=2, backend="dense",
                            primes=[p1, p2], tag="battery-k2")
    _check(v2["verdict"] == "ESCALATE-TO-EXACT",
           f"k=2 verdict {v2['verdict']}", errs)
    _check(v2["agree"] is False and v2.get("rank") is None,
           "k=2 disagreement not surfaced", errs)
    ranks2 = sorted(d["rank"] for d in v2.get("disagreement", []))
    _check(ranks2 == sorted([truth["rank_at_p"], truth["rank_generic"]]),
           f"k=2 per-prime ranks {ranks2}", errs)
    # two of three primes agree — escalation must STILL fire (no majority)
    v3 = screen.screen_rows(rows, labels, ncols, k=3, backend="dense",
                            primes=[p1, p2, p3], tag="battery-k3")
    _check(v3["verdict"] == "ESCALATE-TO-EXACT",
           f"k=3 (2-vs-1) verdict {v3['verdict']} — majority vote is "
           f"forbidden", errs)
    ranks3 = sorted(d["rank"] for d in v3.get("disagreement", []))
    _check(ranks3 == sorted([truth["rank_at_p"], truth["rank_generic"],
                             truth["rank_generic"]]),
           f"k=3 per-prime ranks {ranks3}", errs)
    if errs:
        raise AssertionError("; ".join(errs))
    return (f"rank {truth['rank_at_p']} vs {truth['rank_generic']} across "
            f"primes -> ESCALATE-TO-EXACT at k=2 and at k=3 (2-vs-1)")


def leg_bad_denominator(screen, shim, fixtures, modp):
    p1, p2 = modp.PRIMES_DENSE[:2]
    errs = []
    # a) pool prime divides a denominator: detect, replace, receipt
    rows, labels, ncols, nunk, truth = fixtures.planted_bad_denominator(p1)
    v = screen.screen_rows(rows, labels, ncols, k=2, backend="dense",
                           tag="battery-replace")
    rep = v.get("primes_replaced", [])
    _check(any(r["prime"] == p1 and "denominator" in r["why"] for r in rep),
           f"bad prime {p1} not receipted in primes_replaced: {rep}", errs)
    _check(p1 not in v["primes"],
           "bad prime used in the verdict panel", errs)
    _check(v["verdict"] == "CERTIFIED-SCREEN"
           and v.get("rank") == truth["rank"]
           and v.get("pivot_cols") == truth["pivot_cols"]
           and v.get("incon_idx") == [],
           f"verdict after replacement wrong: {v['verdict']} "
           f"rank {v.get('rank')}", errs)
    # b) every panel prime unusable: refuse, never emit a verdict
    rows2, labels2, ncols2, _, _ = fixtures.planted_bad_denominator(p1 * p2)
    refused = False
    try:
        screen.screen_rows(rows2, labels2, ncols2, k=2, backend="dense",
                           primes=[p1, p2], tag="battery-refuse")
    except modp.BadPrime as e:
        refused = "exhausted" in str(e)
    _check(refused, "fully-unusable prime panel did not refuse", errs)
    if errs:
        raise AssertionError("; ".join(errs))
    return (f"p|denominator detected: {p1} replaced + receipted, planted "
            f"truth recovered; all-unusable panel refused")


LEG_FN = {"planted_control": leg_planted_control,
          "zero_eq_one": leg_zero_eq_one,
          "prime_disagree": leg_prime_disagree,
          "bad_denominator": leg_bad_denominator}

# One targeted sabotage per leg: (file, exact original text, mutated text).
# Each disables precisely the mechanism its leg tests; applied only to a
# scratch copy by --mutation-controls.
MUTATIONS = {
    "planted_control": [
        ("modp_rref.py", "c0 = min(rdict)", "c0 = max(rdict)"),
        ("modp_rref.py", "c0 = int(nzc[0])", "c0 = int(nzc[-1])"),
    ],
    "zero_eq_one": [
        ("modp_rref.py", "elif rv:\n            incon.append(idx)",
         "elif False:\n            incon.append(idx)"),
        ("modp_rref.py", "elif buf[ncols] % p:\n            incon.append(idx)",
         "elif False:\n            incon.append(idx)"),
    ],
    "prime_disagree": [
        ("screen.py",
         "verdict = 'CERTIFIED-SCREEN' if agree else 'ESCALATE-TO-EXACT'",
         "verdict = 'CERTIFIED-SCREEN'"),
    ],
    "bad_denominator": [
        ("modp_rref.py", "if den % p == 0:", "if False and den % p == 0:"),
    ],
}


# ---------------------------------------------------------------- drivers ---
def run_leg_inproc(leg, tool_dir):
    screen, shim, fixtures, modp = _load_tool(tool_dir)
    return LEG_FN[leg](screen, shim, fixtures, modp)


def run_leg_subproc(leg, tool_dir):
    """One leg in a child process (isolates imports across tool copies)."""
    t0 = time.time()
    r = subprocess.run(
        [sys.executable, "-B", os.path.abspath(__file__), "--run-leg", leg,
         "--tool-dir", tool_dir],
        capture_output=True, text=True, timeout=600)
    detail = (r.stdout.strip().splitlines() or ["(no output)"])[-1]
    return {"leg": leg, "pass": r.returncode == 0,
            "wall_s": round(time.time() - t0, 2),
            "detail": detail if r.returncode == 0 else
            (r.stdout + r.stderr).strip()[-400:]}


def mutation_controls(tool_dir, legs, workdir):
    """Copy tool_dir to scratch, sabotage one mechanism per leg, require
    that leg to FAIL; delete the copy. The installed tree is untouched."""
    results = []
    for leg in legs:
        clone = os.path.join(workdir, f"mut_{leg}")
        if os.path.isdir(clone):
            shutil.rmtree(clone)
        os.makedirs(clone)
        for f in os.listdir(tool_dir):
            if f.endswith(".py") or f == "README.md":
                shutil.copy(os.path.join(tool_dir, f),
                            os.path.join(clone, f))
        for fname, old, new in MUTATIONS[leg]:
            path = os.path.join(clone, fname)
            src = open(path).read()
            assert src.count(old) == 1, \
                f"mutation anchor not unique in {fname}: {old!r}"
            open(path, "w").write(src.replace(old, new))
        res = run_leg_subproc(leg, clone)
        caught = not res["pass"]
        print(f"MUTATION {leg}: {'CAUGHT' if caught else 'MISSED'} "
              f"[{res['wall_s']}s]", flush=True)
        results.append({"leg": leg, "caught": caught,
                        "mutations": [f"{f}: {o!r} -> {n!r}"
                                      for f, o, n in MUTATIONS[leg]],
                        "wall_s": res["wall_s"]})
        shutil.rmtree(clone)
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tool-dir", default=HERE)
    ap.add_argument("--legs", default=",".join(LEGS))
    ap.add_argument("--run-leg", default=None, help="internal: one leg")
    ap.add_argument("--mutation-controls", action="store_true",
                    help="also sabotage each leg's mechanism on a scratch "
                         "copy and require that leg to FAIL")
    ap.add_argument("--summary-out", default=None,
                    help="report JSON path (default: "
                         "RANKSCREEN_BATTERY_SUMMARY.json in the cwd)")
    a = ap.parse_args()
    tool_dir = os.path.abspath(a.tool_dir)

    if a.run_leg:                                   # child mode
        try:
            detail = run_leg_inproc(a.run_leg, tool_dir)
            print(detail, flush=True)
            return 0
        except AssertionError as e:
            print(f"FAILED: {e}", flush=True)
            return 1

    cwd = os.path.abspath(os.getcwd())
    if cwd == tool_dir or cwd.startswith(tool_dir + os.sep):
        print("refusing to run inside the tool tree — run from a scratch "
              "directory", flush=True)
        return 2

    t0 = time.time()
    legs = [l for l in a.legs.split(",") if l]
    for l in legs:
        if l not in LEG_FN:
            raise SystemExit(f"unknown leg {l} (have {', '.join(LEGS)})")
    results, fails = [], []
    for leg in legs:
        res = run_leg_subproc(leg, tool_dir)
        results.append(res)
        if not res["pass"]:
            fails.append(leg)
        print(f"{'PASS' if res['pass'] else 'FAIL'}  {leg}  "
              f"[{res['wall_s']}s]  {res['detail'][:160]}", flush=True)
    report = {"tool": "rankscreen", "battery": "synthetic-v1",
              "tool_dir": tool_dir, "legs": results}
    if a.mutation_controls:
        mut = mutation_controls(tool_dir, legs, cwd)
        report["mutation_controls"] = mut
        missed = [m["leg"] for m in mut if not m["caught"]]
        fails += [f"mutation-missed:{l}" for l in missed]
    report["total_s"] = round(time.time() - t0, 1)
    report["overall"] = "PASS" if not fails else "FAIL"
    outp = a.summary_out or os.path.join(cwd,
                                         "RANKSCREEN_BATTERY_SUMMARY.json")
    json.dump(report, open(outp, "w"), indent=1)
    print(f"summary -> {outp}", flush=True)
    if fails:
        print(f"OVERALL FAIL: {fails} ({report['total_s']}s)", flush=True)
        return 1
    print(f"OVERALL PASS ({len(legs)} legs"
          + (f" + {len(legs)} mutation controls" if a.mutation_controls
             else "") + f", {report['total_s']}s)", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
