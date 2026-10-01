#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision.
# Flat-path alias of tools/gatekeeper/probes/quad_probe.py: both the import form and
# the CLI form (python3 tools/quad_probe.py ...) at this path forward to
# gatekeeper.probes.quad_probe. The package is resolved from this file's own
# directory (tools/), whichever sys.path entry found the alias.
import os as _os
import sys as _sys
_d = _os.path.dirname(_os.path.abspath(__file__))
if _d not in _sys.path:
    _sys.path.insert(0, _d)
import gatekeeper.probes.quad_probe as _m
_sys.modules[__name__].__dict__.update({k: v for k, v in _m.__dict__.items() if not k.startswith('__')})
if __name__ == '__main__':
    _sys.exit(_m.main())
