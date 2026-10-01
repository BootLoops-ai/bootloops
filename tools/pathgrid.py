#!/usr/bin/env python3
"""Compatibility shim — canonical module: the tools/ratfit/ package.
coincidence_loci and design_grid live in tools/ratfit/__init__.py."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ratfit import coincidence_loci, design_grid  # noqa: F401
