"""baller.hygiene — precision-hygiene lints and ball-footgun checks.

Checks (baller.hygiene_checks):
  lint_fixed_wdps · lint_ambient_conv · lint_unkeyed_cache ·
  lint_ctx_restore · ftrim + RadiusWatch (radius discipline) ·
  bridge_mpf (re-round-safe mpf bridge).

MEMBER MODULES (beside this one; this wing is their only home):
  dps_lint     AST linter for the mpmath import-time-dps footgun class —
               module-level mp constants frozen at import-time dps=15,
               cross-module dps resets ([imported-sets-dps],
               [import-frozen-const], [pre-main-const]) and complex-step
               mp.diff on real-branched callables. Stdlib-only.
               CLI: python3 baller/dps_lint.py FILE_OR_DIR [...]   (exit 1 on
               findings, 2 on a bad path) · python3 -m baller.hygiene lint ...
               · --selftest (built-in battery).
  purity_scan  standalone finite-precision-literal scanner — AST+regex
               classification (MPF-OF-FLOAT / DEC-STRING / FLOAT-EXPR /
               auto-pure) over an arbitrary file list; report-only.
               CLI: python3 purity_scan.py <list-of-abs-paths.txt> <out.json>

`python3 -m baller.hygiene lint ARGS...` (tools/baller on sys.path) runs the
dps_lint CLI; `python3 -m baller.hygiene selftest` runs its built-in battery.
"""
from . import dps_lint
from . import purity_scan
from .hygiene_checks import (LINTS, lint_files, lint_fixed_wdps,
                             lint_ambient_conv, lint_unkeyed_cache,
                             lint_ctx_restore, ftrim, RadiusWatch,
                             RadiusBlowup, bridge_mpf, ctx_guard)

__all__ = ["dps_lint", "purity_scan", "LINTS", "lint_files",
           "lint_fixed_wdps", "lint_ambient_conv", "lint_unkeyed_cache",
           "lint_ctx_restore", "ftrim", "RadiusWatch", "RadiusBlowup",
           "bridge_mpf", "ctx_guard"]


def _main(argv):
    """`python3 -m baller.hygiene lint FILE_OR_DIR [...]` / `... selftest`."""
    import sys
    if argv and argv[0] == "lint":
        return dps_lint.main(argv[1:])
    if argv and argv[0] in ("selftest", "--selftest"):
        return dps_lint.selftest()
    print("usage: python3 -m baller.hygiene lint FILE_OR_DIR [FILE_OR_DIR ...]\n"
          "       python3 -m baller.hygiene selftest", file=sys.stderr)
    return 2


if __name__ == "__main__":
    import sys
    sys.exit(_main(sys.argv[1:]))
