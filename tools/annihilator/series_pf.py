#!/usr/bin/env python3
"""Import shim — the implementation lives in annihilator.py
(pf_from_series, rec_to_theta, ...)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from annihilator import *  # noqa: F401,F403
from annihilator import P31  # noqa: F401
if __name__ == '__main__':
    print("series_pf lives in tools/annihilator/annihilator.py — run that for the self-test.")
