"""Shared paths + check plumbing for the ibplapper battery.

The archived reference artifacts (interface-table banks, the retrofit
registry, the wrong-table exhibit, the pinned kira Generate tree) are NOT
shipped with this tree; tests that replay them SKIP cleanly unless the
environment points at a copy:

    WINNOW_BANKED_ROOT  root of the archived reference runs
    WINNOW_T1_ROOT      pinned kira Generate artifact tree (T1 family)

Archived paths are read-only; scratch goes to the out dir passed on
argv."""
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
if PKG not in sys.path:
    sys.path.insert(0, PKG)

RL = os.environ.get("WINNOW_BANKED_ROOT")


def _rl(*parts):
    return os.path.join(RL, *parts) if RL else None


NOMAT = _rl("newmove_nomat_p0")
COPAIR = _rl("copair_day1")
D4 = _rl("strata_day4_smallrung", "staging")
D4OUT = _rl("strata_day4_smallrung", "out")
STAGEPACKS = _rl("lapper_bench_t1test", "stagepacks")
T1GEN = os.environ.get("WINNOW_T1_ROOT")


def skip_without(what, *roots):
    """Exit 0 with a SKIP notice when the archived roots are absent."""
    if any(r is None for r in roots):
        print(f"SKIP: {what} replays archived reference artifacts not"
              " shipped with this tree (set WINNOW_BANKED_ROOT /"
              " WINNOW_T1_ROOT to run it)")
        sys.exit(0)
REGISTRY = "e2d8acded38e6885"
PRIMES = [2147483647, 2147483629]
D0, E0 = 1234577, 87654321
SLICES4 = [(2147483647, 1234577, 87654321),
           (2147483647, 987654323, 1029384757),
           (2147483629, 1234577, 87654321),
           (2147483629, 987654323, 1029384757)]

FAILED = []
_T0 = time.time()


def out_dir(default):
    d = sys.argv[1] if len(sys.argv) > 1 else default
    os.makedirs(d, exist_ok=True)
    return d


def check(name, cond, extra=""):
    tag = "PASS" if cond else "FAIL"
    print(f"[{tag}] {name} {extra}", flush=True)
    if not cond:
        FAILED.append(name)
    return bool(cond)


def finish(label):
    wall = round(time.time() - _T0, 1)
    if FAILED:
        print(f"\n{label}: FAILURES {FAILED} ({wall}s)")
        sys.exit(1)
    print(f"\n{label}: ALL PASS ({wall}s)")
    sys.exit(0)


def jload(path):
    return json.load(open(path))
