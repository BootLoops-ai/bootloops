"""Import shim — the implementation lives in annihilator.py
(nullspace_mod_fast, find_recurrence_fast)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from annihilator import nullspace_mod_fast, find_recurrence_fast  # noqa: F401
