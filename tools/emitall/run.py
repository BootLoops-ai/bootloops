#!/usr/bin/env python3
"""emitall CLI.

    python3 tools/emitall/run.py SPEC.json \
        [--out report.json] [--quiet]

Exit code: 0 when every quoted headline matches and every receipt is present
and fresh; 1 on ANY finding (mismatch, outward-band violation, eval error,
missing receipt, stale receipt). 2 on a malformed spec.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from emitall.engine import run_spec, SpecError  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description="emit-all headline recount")
    ap.add_argument("spec", help="claims spec (.json / .yaml)")
    ap.add_argument("--out", default=None, help="write JSON report here")
    ap.add_argument("--quiet", action="store_true", help="suppress the table")
    args = ap.parse_args(argv)
    try:
        _, code = run_spec(args.spec, out_path=args.out, quiet=args.quiet)
    except SpecError as e:
        print(f"SPEC ERROR: {e}", file=sys.stderr)
        return 2
    return code


if __name__ == "__main__":
    sys.exit(main())
