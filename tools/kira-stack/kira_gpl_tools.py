#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision.
"""kira_gpl_tools — locate and import the two Kira-derived GPL utilities.

`pyred_weight.py` (codec for Kira's pyRed integral-weight layout) and
`ffsave_degree_census.py` (degree census of a Kira + FireFly `ff_save/` state)
were written by reading Kira's and FireFly's GPL sources, so they are
distributed under the GPL-3.0-or-later in the sibling repository
`kira`, directory `bootloops-tools/`, and not in this MIT tree.
This loader finds that directory and imports them on demand:

  1. the directory named by the environment variable BOOTLOOPS_KIRA_TOOLS,
     if set;
  2. otherwise <this repository's root>/../kira/bootloops-tools,
     which is where they sit when the two repositories are checked out side
     by side.

Nothing falls back silently: `load(name)` raises ImportError naming the
sibling repository and the variable when the file is not there, and
`available()` lets a battery skip its dependent legs by name.

    from kira_gpl_tools import load, available, script_path
    pw = load('pyred_weight')                  # module
    census = script_path('ffsave_degree_census')   # path, to run as a CLI
"""
import importlib.util
import os
import sys

ENV = 'BOOTLOOPS_KIRA_TOOLS'
NAMES = ('pyred_weight', 'ffsave_degree_census')
HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(HERE, os.pardir, os.pardir))
SIBLING = os.path.join(os.path.dirname(REPO_ROOT), 'kira',
                       'bootloops-tools')


def tools_dir():
    """The directory searched (env override first, else the sibling checkout)."""
    d = os.environ.get(ENV)
    return os.path.abspath(d) if d else SIBLING


def script_path(name, must_exist=True):
    if name not in NAMES:
        raise ValueError(f'kira_gpl_tools: unknown tool {name!r}; known: {NAMES}')
    p = os.path.join(tools_dir(), name + '.py')
    if must_exist and not os.path.isfile(p):
        raise ImportError(
            f'kira_gpl_tools: {name}.py not found in {tools_dir()}. It is a '
            f'GPL-3.0-or-later utility shipped in the sibling repository '
            f'kira (directory bootloops-tools/): check that '
            f'repository out beside this one (../kira) or set '
            f'{ENV} to the directory that holds it.')
    return p


def available(name=None):
    """True when the named tool (or every tool) can be loaded."""
    names = (name,) if name else NAMES
    return all(os.path.isfile(script_path(n, must_exist=False)) for n in names)


def load(name):
    """Import the tool as a module by path (the directory is also put on
    sys.path so the census can import the codec beside it)."""
    p = script_path(name)
    d = os.path.dirname(p)
    if d not in sys.path:
        sys.path.insert(0, d)
    key = f'_kira_gpl_{name}'
    if key in sys.modules:
        return sys.modules[key]
    spec = importlib.util.spec_from_file_location(key, p)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    sys.modules[key] = m
    return m


if __name__ == '__main__':
    d = tools_dir()
    print(f'kira_gpl_tools: directory {d} '
          f'({"env " + ENV if os.environ.get(ENV) else "sibling default"})')
    for n in NAMES:
        print(f'  {n:22s} {"present" if available(n) else "ABSENT"}')
    sys.exit(0 if available() else 1)
