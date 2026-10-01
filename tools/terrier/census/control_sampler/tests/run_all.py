"""control_sampler test battery (synthetic inputs only; touches NO census
data). Needs TERRIER_SEEDS_JSON (seeds file, not included in the package).
Run:  python3 tests/run_all.py   (single core, ~25 s)
"""
import importlib
import os
import sys
import time
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))

MODS = ["test_seeds", "test_cp_tables", "test_closedform", "test_no_exclusion"]


def main():
    n_pass, n_fail = 0, 0
    for mname in MODS:
        mod = importlib.import_module(mname)
        for name in sorted(dir(mod)):
            if not name.startswith("test_"):
                continue
            t0 = time.time()
            try:
                getattr(mod, name)()
                n_pass += 1
                print("PASS %-45s (%.2fs)" % (mname + "." + name,
                                              time.time() - t0))
            except Exception:
                n_fail += 1
                print("FAIL %s" % (mname + "." + name))
                traceback.print_exc()
    print("SUMMARY: %d pass / %d fail" % (n_pass, n_fail))
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
