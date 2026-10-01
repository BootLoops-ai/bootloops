"""python3 -m mixalot.boxwalk <selftest|plan|produce|fiber|verify> ...
(run with tools/mixalot on sys.path, e.g. from tools/mixalot/)."""
import sys

from .cli import main

sys.exit(main())
