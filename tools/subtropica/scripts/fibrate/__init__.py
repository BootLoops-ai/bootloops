"""fibrate — sound fibration engine (subtropica gap-tool).

Replaces HyperFLINT's unsound `fibration_basis` op.
See README.md in this directory and fibrate.py's
module docstring.

Public API:
    from fibrate import Engine, run_wordlist, selfcheck, parse_letter
    from fibrate.zip_num import ZipEval
    from fibrate.gpl_num import GEval
"""
from .fibrate import (Engine, FibrateError, GEvalX, parse_letter,
                      run_wordlist, selfcheck, serialize_rep)
from .gpl_num import GEval
from .zip_num import ZipEval

__all__ = ["Engine", "FibrateError", "GEvalX", "GEval", "ZipEval",
           "parse_letter", "run_wordlist", "selfcheck", "serialize_rep"]
