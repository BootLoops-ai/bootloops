"""emitall — declarative emit-all harness for project headline claims.

A project writes ONE claims spec (JSON/YAML) mapping every quoted headline
to (receipt, field expression, comparison mode); the engine re-emits every
value from the stored receipts, compares at the quote's own register, and
exits nonzero on any finding. See README.md; example deployment in specs/.
"""
from .engine import (run_spec, load_spec, compare_printed,
                     compare_outward_band, fmt_like, SpecError)
from .expr import evaluate, resolve, resolve_pattern, Env, ExprError, ReceiptMissing
from .lint import lint_source, lint_file, lint_paths

__version__ = "1.0.0"

REGISTER = {
    "engine": "deterministic recount of stored receipts; no float participates "
              "in any outward-band accept/reject (exact Fraction); 'printed' "
              "mode is byte-compatible with the emitting-script comparison "
              "semantics it audits (fmt_like)",
    "lint": "template lint (emission-integrity rule): digit "
            "literals in emission prose are machine-caught findings "
            "(interpolation placeholders exempt; allow-literal pragma "
            "requires a reason), sampled/enumerated emission fields must be "
            "sorted()/seeded (guards against typed-not-emitted digits and "
            "PYTHONHASHSEED-order emissions)",
    "specs/toy_stars.json": "shipped example deployment — an INVENTED "
              "star-rating program over synthetic receipts in "
              "specs/toy_receipts/ (12 headlines: 11 MATCH + 1 emitted-only); "
              "exercised end-to-end by tests/test_all.py, clean and with "
              "planted corruptions",
    "battery.py": "paper claims battery (document-side verb): declarative "
              "FROZEN/NEVER/NEG/ABSTRACT rows over comment-stripped tex "
              "prose + wrap-proof-joined rendered text; fail-closed "
              "stale-render refusal (rc 2, never PASS against a stale "
              "render); --selftest is a gated 8-gate fixture battery "
              "(planted failure + clean twin per row class), ridden by "
              "tests/test_all.py under python3 and python3 -O",
    "paper_seams.py": "comment-swallow seam gate (document-side member beside "
              "battery.py): S1/S1b/S2/S3 rows over a master's \\input roster — "
              "a % comment block dropped inside a sentence strands the sentence "
              "head on the comment's last line; tex rows plus a pdftotext row "
              "confirmed in the de-commented tex; --allow names a legitimate "
              "lowercase sentence head; exit 1 on any tex-confirmed swallow; "
              "fixture battery tests/test_paper_seams.py over specs/seams/ "
              "(planted S1 -> 1, cured -> 0, head rule)",
}

__all__ = ["run_spec", "load_spec", "compare_printed", "compare_outward_band",
           "fmt_like", "evaluate", "resolve", "resolve_pattern", "Env",
           "ExprError", "ReceiptMissing", "SpecError", "REGISTER",
           "lint_source", "lint_file", "lint_paths", "__version__"]
