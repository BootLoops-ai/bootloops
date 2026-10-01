#!/usr/bin/env python3
"""Merge validation case logs + upstream timings into validation/REPORT.md."""
import os
import re
import sys

LOGDIR = sys.argv[1] if len(sys.argv) > 1 else "/tmp/sofia_val"

cases = {}
for line in open("validation/notebook_cases.tsv"):
    f = line.rstrip("\n").split("\t")
    i = len(cases)
    cases[i] = {"input": f[2][:100], "section": f[4] if len(f) > 4 else ""}

for line in open("validation/upstream_timings.tsv"):
    i, wall, sec = line.rstrip("\n").split("\t")
    cases[int(i)]["upstream_s"] = wall

rows = []
for i in sorted(cases):
    log = os.path.join(LOGDIR, f"case_{i}.log")
    c = cases[i]
    status, t, miss, extra, counts = "-", "", "", "", ""
    if os.path.exists(log):
        txt = open(log).read()
        m = re.search(r"-> (\w+)", txt)
        status = m.group(1) if m else ("timeout" if txt.strip() else "-")
        m = re.search(r"time ([0-9.]+)s", txt)
        t = m.group(1) if m else ""
        m = re.search(r"mine=(\d+) expected=(\d+) missing=(\d+) extra=(\d+)", txt)
        if m:
            counts = f"{m.group(1)}/{m.group(2)}"
            miss, extra = m.group(3), m.group(4)
    rows.append((i, c["section"][:40], status, counts, miss, extra,
                 t, c.get("upstream_s", "")))

with open("validation/REPORT.md", "w") as f:
    f.write("# Validation against SOFIA_examples.nb recorded outputs\n\n")
    f.write("Each case replays a worked example from the upstream notebook "
            "(same diagram, same options) and compares the candidate-singularity "
            "sets up to proportionality and variable naming. `superset` means "
            "every upstream singularity was reproduced plus extras (acceptable "
            "for candidate semantics); route-variance differences reflect the "
            "elimination-route dependence both implementations share (see the "
            "notebook's own automatic-vs-pinned LoopEdges demonstration).\n\n")
    f.write("| case | section | status | mine/exp | missing | extra | ours (s) | upstream (s) |\n")
    f.write("|---|---|---|---|---|---|---|---|\n")
    for r in rows:
        f.write("| " + " | ".join(str(x) for x in r) + " |\n")
    f.write("\nUpstream timings are from the notebook's EchoTiming cells "
            "(authors' machine); ours from this batch (single core per case).\n")
print("wrote validation/REPORT.md")
for r in rows:
    print(r)
