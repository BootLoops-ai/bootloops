#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision.
"""dipstick `regions` member — standalone test entry: python3 selftest.py

Thin wrapper over selftest.jl (the GATES control pair for regions.jl, the
expansion-by-regions completeness certifier). The package battery runs the
same leg as `python3 tools/dipstick/battery.py --legs regions`; this file is
the member-local door. With no `julia` on PATH this is a loud NAMED skip
(exit 0) — the certifier is a Julia script and cannot run at all; that skip is
the designed engine-absent behavior for an engine-gated member, not breakage.
With julia present, selftest.jl takes over: it looks for the Oscar + JSON
packages in the active environment, then DIPSTICK_JULIA_PROJECT, then the
committed project at ../../../upgrades/leviathan, and itself SKIPS by name
(with the one-time instantiate command) when none provides them. Exit codes
follow selftest.jl: 0 = both controls green, or a named skip; 1 = a genuine
control failure or a named fail-closed refusal.

DIPSTICK_JULIA overrides the julia invocation (e.g. "julia +1.10" to pin a
juliaup channel), as for the rest of dipstick.
"""
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    julia_cmd = os.environ.get('DIPSTICK_JULIA', 'julia').split()
    if not julia_cmd or shutil.which(julia_cmd[0]) is None:
        print('dipstick regions selftest: SKIPPED (named engine gate) — no '
              '`julia` on PATH (%r). regions.jl is a Julia script (julia >= '
              '1.10, plus the Oscar + JSON packages; see tools/dipstick/GUIDE.md, '
              '`regions` member). Install julia (e.g. via juliaup), then rerun. '
              'This skip is the designed engine-absent behavior for an '
              'engine-gated member, not breakage. Nothing was verified.'
              % (julia_cmd[0] if julia_cmd else 'julia'))
        return 0
    return subprocess.run(julia_cmd + ['--startup-file=no',
                                       os.path.join(HERE, 'selftest.jl')],
                          cwd=HERE).returncode


if __name__ == '__main__':
    sys.exit(main())
