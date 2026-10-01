# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision.
"""python3 -m loomcheck [--selftest | --basis FILE ...] (cwd = tools/dogtag)."""
import sys

from .loomcheck import main

if __name__ == "__main__":
    sys.exit(main())
