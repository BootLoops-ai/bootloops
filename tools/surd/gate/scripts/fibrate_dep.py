"""Locate `zip_num` (shuffle-regularized ZIP_reg / G(w;1) series evaluator with adaptive series length), which lives in the sibling
package tools/subtropica (scripts/fibrate/zip_num.py, pure Python, mpmath only).  Search order: $SURD_FIBRATE, then
<repo>/tools/subtropica/scripts/fibrate resolved relative to this file.  Importers write `from fibrate_dep import zip_num`."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))          # <repo>/tools  (this file is tools/surd/gate/scripts/fibrate_dep.py)


def fibrate_dir():
    for c in (os.environ.get("SURD_FIBRATE"), os.path.join(TOOLS, "subtropica", "scripts", "fibrate")):
        if c and os.path.exists(os.path.join(c, "zip_num.py")):
            return os.path.abspath(c)
    raise ImportError("surd: zip_num.py not found (expected in the sibling package tools/subtropica/scripts/fibrate); "
                      "set SURD_FIBRATE to the directory that holds zip_num.py")


_P = fibrate_dir()
if _P not in sys.path:
    sys.path.insert(0, _P)
import zip_num          # noqa: E402
