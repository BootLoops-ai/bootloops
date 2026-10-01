#!/usr/bin/env python3
"""Compare two amflow_cli solve_integrals output JSONs, integral by integral,
order by order.  Reports matched decimal digits (relative) per coefficient.

Usage: compare_amflow_json.py REF.json TEST.json [--min-digits 25] [--n-integrals N]
                              [--goal G | --goal GREF,GTEST | --goal-from-config REF_CFG[,TEST_CFG]
                               | --radii-only] [--strict-window]

Exit 0 iff every compared coefficient matches to >= min-digits
(zero-vs-zero counts as a pass; zero-vs-nonzero fails).

Two rules decide what is COMPARED (everything else is named and counted,
never failed):

NUMERICAL ZERO.  Every coefficient is read as an arb ball "[mid +/- rad]"
(a plain number or an empty ball has radius 0).  A pair whose two members
are BOTH below the zero threshold
    T = max(rad_ref, rad_test, 10^-(G-2)),   G = the lower of the two goals,
is a NUMERICAL ZERO: printed by name with both magnitudes and T, excluded
from the digit count and counted on the "numerical zeros: N (...)" line.
A pair with exactly one member above T is zero-vs-nonzero and FAILS (as
before, with its digit count printed).  A pair whose two midpoints are
exactly 0 is "both exactly 0 -> PASS" (as before).

THE GOAL has three sources, tried in this order, and the "goal:" line
names the one used:
  1. --goal G (one value for both files) or --goal GREF,GTEST; wins over
     every other source.
  2. --goal-from-config REF_CFG[,TEST_CFG]: the amflow_cli input configs
     (one path applies to both files); the goal is the config's top-level
     key "goal_digits" (a positive number; the amflow_cli input form), or
     "goal" when a config carries that key instead (the smoke fixtures'
     form); the key read is named.  A config path that does not exist
     exits 4 by name; a config without either key, or with a non-numeric
     one, exits 2 by name.
  3. with neither flag: a file carrying one of the keys in GOAL_KEYS at its
     top level or under "options" is read from there (the key named); the
     solve_integrals output form carries NONE of these keys (its "options"
     carry working_pre / chop_pre / rationalize_pre, working precisions,
     never read as the goal), so for such a file the config is looked up BY
     NAME beside it: out/<name>.json <-> configs/<name>.json (the config at
     ../configs/<name>.json relative to the out file's directory, read as
     in 2.); every location tried is named on the "goal:" line.
When no goal is derivable on either side the run REFUSES by name, exit 2:
  no goal: pass --goal G|GREF,GTEST, --goal-from-config REF_CFG[,TEST_CFG], or --radii-only
and compares nothing.  --radii-only is the explicit opt-in to the
threshold from the two radii alone (the previous behaviour without a
goal: a finite integral's pole-order zeros then read as "0.0 digits ->
FAIL"); it cannot be combined with --goal or --goal-from-config.

MISSING BY WINDOW.  An order present in one file only (the eps window a
solve returns depends on its goal) is "not compared (present in REF/TEST
only; outside the other member's window)", counted on the "orders not
compared: N (...)" line and never a FAIL.  --strict-window restores the
previous rule for such orders: "only one side, but zero -> PASS" when the
present member is exactly zero, otherwise "MISSING on one side, nonzero
-> FAIL".

OVERALL: PASS / FAIL over the compared components (a component = the re
or the im part of one order), printed with the counts beside it:
"OVERALL: PASS (10 compared, 4 numerical zeros excluded, 2 orders not
compared)".  Exit 0 iff every compared component is >= --min-digits (and,
under --strict-window, no order is missing with a nonzero member); 1
otherwise; 2 for a refused option, a malformed --goal, a config without a
usable goal_digits, or no derivable goal; 4 for a --goal-from-config path
that does not exist.
"""
import argparse
import json
import os
import sys
from mpmath import mp, mpf

mp.dps = 250

# keys read as the goal (decimal digits) when a file carries one: top level first, then under "options"
GOAL_KEYS = ("goal_digits", "goal")
# the keys read from an amflow_cli input config (--goal-from-config and the by-name lookup): top level only, first found
CONFIG_GOAL_KEYS = GOAL_KEYS
# the by-name lookup: out/<name>.json <-> configs/<name>.json
CONFIG_DIR_BY_NAME = "configs"
NO_GOAL_MSG = ("no goal: pass --goal G|GREF,GTEST, --goal-from-config REF_CFG[,TEST_CFG], or --radii-only")


def parse_arb(s):
    """Parse an arb string like '[1.23 +/- 4e-50]', '[+/- 4e-50]', '0' or '1.5'
    -> (mid, radius) as mpf; radius 0 for a plain number or an empty ball."""
    s = s.strip()
    if s.startswith("["):
        body = s[1:-1]
        if "+/-" in body:
            mid, rad = body.split("+/-", 1)
        else:
            mid, rad = body, "0"
        mid = mid.strip()
        rad = rad.strip()
        if mid in ("", "+", "-"):
            mid = "0"
        if rad == "":
            rad = "0"
        return mpf(mid), abs(mpf(rad))
    return mpf(s), mpf(0)


def coeffs_of(result_entry):
    """order -> ((re_mid, re_rad), (im_mid, im_rad))"""
    out = {}
    for c in result_entry["coefficients"]:
        v = c["value"]
        out[int(c["order"])] = (parse_arb(v["re"]), parse_arb(v["im"]))
    return out


def matched_digits(a, b):
    """Relative matched decimal digits between mpf a (ref) and b (test)."""
    if a == 0 and b == 0:
        return None  # exact zero match
    denom = max(abs(a), abs(b))
    diff = abs(a - b)
    if diff == 0:
        return 250
    return float(mp.log10(denom / diff))


def goal_of_file(doc):
    """(goal, key) from the first of GOAL_KEYS found at the top level, then under
    'options'; (None, None) when the file carries none (the solve_integrals output
    form carries none)."""
    for key in GOAL_KEYS:
        if isinstance(doc, dict) and isinstance(doc.get(key), (int, float)) and not isinstance(doc.get(key), bool) and doc[key] > 0:
            return int(doc[key]), key
    opts = doc.get("options") if isinstance(doc, dict) else None
    if isinstance(opts, dict):
        for key in GOAL_KEYS:
            if isinstance(opts.get(key), (int, float)) and not isinstance(opts.get(key), bool) and opts[key] > 0:
                return int(opts[key]), "options." + key
    return None, None


def goal_of_config(path):
    """(goal, key, error) from an amflow_cli input config: the first of CONFIG_GOAL_KEYS at
    the top level as a positive number -> (int, key, None); otherwise (None, None, the
    reason by name).  The caller has checked that the path exists."""
    try:
        doc = json.load(open(path))
    except (ValueError, UnicodeDecodeError) as e:
        return None, None, f"config {path} is not JSON ({e.__class__.__name__})"
    if not isinstance(doc, dict):
        return None, None, f"config {path} is not a JSON object"
    keys = [k for k in CONFIG_GOAL_KEYS if k in doc]
    if not keys:
        return None, None, (f"config {path} carries none of {', '.join(CONFIG_GOAL_KEYS)} "
                            f"at its top level")
    key = keys[0]
    g = doc[key]
    if isinstance(g, bool) or not isinstance(g, (int, float)) or g <= 0:
        return None, None, f"config {path} carries {key} = {g!r}, not a positive number"
    return int(g), key, None


def config_by_name(out_path):
    """the config looked up by name beside an out file: <dir of out>/../configs/<name>.json"""
    d = os.path.dirname(os.path.abspath(out_path))
    return os.path.join(os.path.dirname(d), CONFIG_DIR_BY_NAME, os.path.basename(out_path))


def parse_goal(text):
    """--goal G or GREF,GTEST as positive integers -> (GREF, GTEST); anything else is
    refused by argparse (rc 2) with the offending text named."""
    parts = text.split(",")
    if len(parts) not in (1, 2):
        raise argparse.ArgumentTypeError(
            f"expected G or GREF,GTEST as positive integers, got {text!r}")
    vals = []
    for p in parts:
        p = p.strip()
        if not p.isdigit() or int(p) <= 0:
            raise argparse.ArgumentTypeError(
                f"expected G or GREF,GTEST as positive integers, got {text!r}")
        vals.append(int(p))
    if len(vals) == 1:
        vals = vals * 2
    return vals[0], vals[1]


def parse_config_paths(text):
    """--goal-from-config REF_CFG or REF_CFG,TEST_CFG -> (REF_CFG, TEST_CFG); one path
    applies to both files; an empty member is refused by argparse (rc 2) by name."""
    parts = [p.strip() for p in text.split(",")]
    if len(parts) not in (1, 2) or any(p == "" for p in parts):
        raise argparse.ArgumentTypeError(
            f"expected REF_CFG or REF_CFG,TEST_CFG as file paths, got {text!r}")
    if len(parts) == 1:
        parts = parts * 2
    return parts[0], parts[1]


def zero_threshold(rad_a, rad_b, goal):
    """T = max(rad_ref, rad_test, 10^-(goal-2)); the goal term absent when no goal is known."""
    t = max(rad_a, rad_b)
    if goal is not None:
        t = max(t, mpf(10) ** (-(goal - 2)))
    return t


def fmt(x):
    return mp.nstr(x, 2, min_fixed=0, max_fixed=0)


def resolve_goals(args, ref_doc, test_doc):
    """(goal_ref, goal_test, goal_src, tried) by the three sources in order; goal_src is the
    text of the 'goal:' line's source clause, tried the by-name locations looked at and
    absent.  Exits 4 / 2 by name on a --goal-from-config path that does not exist / a
    config without a usable goal_digits (the by-name config likewise)."""
    prog = os.path.basename(sys.argv[0]) or "compare_amflow_json.py"
    if args.goal is not None:
        src = "--goal"
        if args.goal_from_config is not None:
            src += " (--goal-from-config not read: --goal wins)"
        return args.goal[0], args.goal[1], src, []
    if args.goal_from_config is not None:
        goals, toks = [], []
        for path, side in zip(args.goal_from_config, ("REF", "TEST")):
            if not os.path.isfile(path):
                print(f"{prog}: no config: {path} ({side}, --goal-from-config) does not exist or is not a file", file=sys.stderr)
                sys.exit(4)
            g, key, err = goal_of_config(path)
            if err:
                print(f"{prog}: no goal: {err} ({side}, --goal-from-config)", file=sys.stderr)
                sys.exit(2)
            goals.append(g)
            toks.append(f"--goal-from-config {path} {key} ({side})")
        return goals[0], goals[1], ", ".join(toks), []
    # neither flag: the file's own key, else the config by name beside it
    goals, toks, tried, all_file_keys = [], [], [], True
    for doc, path, side in ((ref_doc, args.ref, "REF"), (test_doc, args.test, "TEST")):
        g, key = goal_of_file(doc)
        if g is not None:
            goals.append(g)
            toks.append((f"{key} ({side})", f"file key {key} ({side})"))
            continue
        cfg = config_by_name(path)
        if os.path.isfile(cfg):
            g, key, err = goal_of_config(cfg)
            if err:
                print(f"goal: none ({side}: config by name {cfg} found; {err})")
                print(f"{prog}: no goal: {err} ({side}, config by name)", file=sys.stderr)
                sys.exit(2)
            goals.append(g)
            toks.append((None, f"config by name {cfg} {key} ({side})"))
            all_file_keys = False
        else:
            goals.append(None)
            tried.append(f"{cfg} ({side}, absent)")
    if all_file_keys:
        src = "file key " + ", ".join(t[0] for t in toks) if toks else "none"
    else:
        src = ", ".join(t[1] for t in toks)
    return goals[0], goals[1], src, tried


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ref")
    ap.add_argument("test")
    ap.add_argument("--min-digits", type=float, default=25.0)
    ap.add_argument("--n-integrals", type=int, default=0,
                    help="compare only the first N integrals (0 = all)")
    ap.add_argument("--goal", type=parse_goal, default=None, metavar="G|GREF,GTEST",
                    help="the goal (decimal digits) of REF and TEST: one value for both, or "
                         "GREF,GTEST; sets the zero threshold T = max(radii, 10^-(G-2)) with G "
                         "the lower of the two.  Wins over --goal-from-config and over a goal "
                         "read from the files (keys " + ", ".join(GOAL_KEYS) + " at the top "
                         "level or under 'options'; the solve_integrals output form carries none).")
    ap.add_argument("--goal-from-config", type=parse_config_paths, default=None,
                    metavar="REF_CFG[,TEST_CFG]",
                    help="read the goal from the amflow_cli input config(s): the top-level key "
                         + " or ".join(CONFIG_GOAL_KEYS) + " (a positive number; the key read is "
                         "named); one path applies to both files.  A path that does not exist "
                         "exits 4 by name; a config without a usable key exits 2 by name.  "
                         "Without this flag "
                         "and without --goal the config is looked up BY NAME beside each out "
                         "(out/<name>.json <-> " + CONFIG_DIR_BY_NAME + "/<name>.json) and every "
                         "location tried is named on the goal: line.")
    ap.add_argument("--radii-only", action="store_true",
                    help="the explicit opt-in to the zero threshold from the two radii alone "
                         "(the previous behaviour without a goal); without it a run with no "
                         "derivable goal refuses by name (exit 2) and compares nothing.  Cannot "
                         "be combined with --goal or --goal-from-config.")
    ap.add_argument("--strict-window", action="store_true",
                    help="an order present in one file only is 'MISSING on one side, nonzero "
                         "-> FAIL' (the previous rule) instead of 'not compared'")
    args = ap.parse_args()
    if args.radii_only and (args.goal is not None or args.goal_from_config is not None):
        ap.error("--radii-only cannot be combined with --goal or --goal-from-config")

    ref_doc = json.load(open(args.ref))
    test_doc = json.load(open(args.test))
    ref = ref_doc["result"]
    test = test_doc["result"]
    n = min(len(ref), len(test))
    if args.n_integrals:
        n = min(n, args.n_integrals)

    # the goal of each side: --goal, else --goal-from-config, else the file's key, else the config by name
    if args.radii_only:
        goal_ref, goal_test, goal_src, tried = None, None, "none", []
    else:
        goal_ref, goal_test, goal_src, tried = resolve_goals(args, ref_doc, test_doc)
    goals = [g for g in (goal_ref, goal_test) if g is not None]
    goal = min(goals) if goals else None
    if goal is None:
        if args.radii_only:
            print("goal: none (zero threshold = the two radii alone)")
        else:
            print("goal: none (tried by name: " + ", ".join(tried) + ")")
            print(f"{os.path.basename(sys.argv[0]) or 'compare_amflow_json.py'}: {NO_GOAL_MSG}",
                  file=sys.stderr)
            sys.exit(2)
    else:
        print(f"goal: REF {goal_ref if goal_ref is not None else 'none'}, "
              f"TEST {goal_test if goal_test is not None else 'none'} (from {goal_src}"
              + ("; tried by name: " + ", ".join(tried) if tried else "") + "); "
              f"zero threshold T = max(radii, 10^-({goal}-2) = {fmt(mpf(10) ** (-(goal - 2)))})")

    def side_label(side):
        g = goal_ref if side == "REF" else goal_test
        return f"present at goal {g} ({side}) only" if g is not None else f"present in {side} only"

    ok = True
    n_compared = 0
    num_zeros = []
    not_compared = []
    missing_strict = []
    for i in range(n):
        rc, tc = coeffs_of(ref[i]), coeffs_of(test[i])
        orders = sorted(set(rc) | set(tc))
        for o in orders:
            if o not in rc or o not in tc:
                name = f"int{i} eps^{o}"
                side = "REF" if o in rc else "TEST"
                if not args.strict_window:
                    print(f"{name}: not compared ({side_label(side)}; "
                          f"outside the other member's window)")
                    not_compared.append(name)
                    continue
                # --strict-window: the previous rule
                present = rc.get(o, tc.get(o))
                if present[0][0] == 0 and present[1][0] == 0:
                    print(f"{name}: only one side, but zero -> PASS")
                    continue
                print(f"{name}: MISSING on one side, nonzero -> FAIL")
                missing_strict.append(name)
                ok = False
                continue
            for part, ((a, ra), (b, rb)) in zip(("re", "im"), zip(rc[o], tc[o])):
                name = f"int{i} eps^{o} {part}"
                d = matched_digits(a, b)
                if d is None:
                    n_compared += 1
                    print(f"{name}: both exactly 0 -> PASS")
                    continue
                t = zero_threshold(ra, rb, goal)
                a_below, b_below = abs(a) < t, abs(b) < t
                if a_below and b_below:
                    print(f"{name}: NUMERICAL ZERO (|ref| {fmt(abs(a))}, |test| {fmt(abs(b))} "
                          f"below T {fmt(t)}) -> not compared")
                    num_zeros.append(name)
                    continue
                n_compared += 1
                if a_below or b_below:
                    lo, hi = ("ref", "test") if a_below else ("test", "ref")
                    print(f"{name}: {d:.1f} digits (ref={a}, test={b}) zero-vs-nonzero "
                          f"(|{lo}| below T {fmt(t)}, |{hi}| above) -> FAIL")
                    ok = False
                elif d >= args.min_digits:
                    print(f"{name}: {d:.1f} digits -> PASS")
                else:
                    print(f"{name}: {d:.1f} digits "
                          f"(ref={a}, test={b}) -> FAIL")
                    ok = False
    print(f"numerical zeros: {len(num_zeros)} ({', '.join(num_zeros)})")
    print(f"orders not compared: {len(not_compared)} ({', '.join(not_compared)})")
    if args.strict_window:
        print(f"orders missing on one side (--strict-window): {len(missing_strict)} "
              f"({', '.join(missing_strict)})")
    print(f"OVERALL: {'PASS' if ok else 'FAIL'} ({n_compared} compared, "
          f"{len(num_zeros)} numerical zeros excluded, {len(not_compared)} orders not compared)")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
