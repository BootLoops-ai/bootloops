#!/usr/bin/env python3
"""rankscreen battery entry — one command for the synthetic battery.

Runs battery.py (four legs + one mutation control per leg) from a scratch
directory, as the battery requires, and propagates its exit code. The battery
is self-contained: fixtures are built deterministically in code with planted
ground truths, so it runs anywhere the tool runs.
"""
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))

with tempfile.TemporaryDirectory(prefix="rankscreen_battery_") as scratch:
    r = subprocess.run(
        [sys.executable, "-B", os.path.join(HERE, "battery.py"),
         "--mutation-controls"], cwd=scratch)
sys.exit(r.returncode)
