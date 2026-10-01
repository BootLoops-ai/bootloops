#!/usr/bin/env python3
"""clinch selftest — runs the acceptance battery from a scratch directory
and grades it against the cold-clone contract in GUIDE.md.

With CLINCH_REFERENCE_DIR set (the reference results, not shipped), the full
battery must pass. Without it, the self-contained legs (L1, L2, L5, L7, L8,
and L4's synthetic controls) must pass, while the receipt-reading legs (L3,
L6, and L4's receipt sub-leg) fail on the absent receipt files by design —
those are reported here as receipt-gated, not as breakage. Exit 0 = the
shipped contract holds; 1 = a self-contained leg broke.
"""
import json
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
BATTERY = os.path.join(HERE, "battery", "battery.py")
RECEIPT_GATED = ("L3", "L4", "L6")
RECEIPT_MARKS = ("L3_RESULTS.json", "CERT.json")
PUBLIC = ("L1", "L2", "L5", "L7", "L8")


def main():
    scratch = tempfile.mkdtemp(prefix="clinch_selftest_")
    r = subprocess.run([sys.executable, BATTERY], cwd=scratch,
                       capture_output=True, text=True, timeout=900)
    sys.stdout.write(r.stdout)
    if r.returncode == 0:
        print("clinch selftest PASS (full battery)")
        return 0
    if (os.environ.get("CLINCH_REFERENCE_DIR")
            or os.environ.get("CLINCH_RECEIPTS_ROOT")
            or os.environ.get("CLINCH_LANE_RECEIPTS")):
        sys.stderr.write(r.stderr)
        print("clinch selftest FAIL: CLINCH_REFERENCE_DIR is set, so the full "
              "battery must pass")
        return 1
    summary_path = os.path.join(scratch, "CLINCH_BATTERY_SUMMARY.json")
    if not os.path.isfile(summary_path):
        sys.stderr.write(r.stderr)
        print(f"clinch selftest FAIL: battery rc={r.returncode} and no "
              "summary written")
        return 1
    summary = json.load(open(summary_path))
    rows = {row["leg"].split()[0]: row for row in summary["results"]}
    bad = []
    for tag, row in rows.items():
        if row["status"] == "FAIL":
            gated = (tag in RECEIPT_GATED
                     and any(m in row["detail"] for m in RECEIPT_MARKS))
            if not gated:
                bad.append(f"{row['leg']}: {row['detail']}")
    for tag in PUBLIC:
        if rows.get(tag, {}).get("status") not in ("PASS", "REMAINDER"):
            bad.append(f"{tag} did not pass")
    if bad:
        print("clinch selftest FAIL:")
        for b in bad:
            print("  " + b)
        return 1
    print("clinch selftest PASS — self-contained legs green; receipt-gated "
          "legs (L3/L6 and L4's receipt sub-leg) skipped: "
          "CLINCH_REFERENCE_DIR unset (the reference results are not "
          "shipped; see GUIDE.md)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
