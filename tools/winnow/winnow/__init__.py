"""Winnow — the exact elimination engine.

Alias package: the canonical import is ``ibplapper`` — receipts are named
``ibplapper.receipt`` and the shipped batteries pin that name; this module
re-exports it verbatim so ``import winnow`` works everywhere
``import ibplapper`` does.
"""
from ibplapper import *  # noqa: F401,F403
import ibplapper as _c
__version__ = _c.__version__
