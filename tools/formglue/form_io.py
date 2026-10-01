#!/usr/bin/env python3
"""form_io — shared line-unwrap fix for FORM output (the #write mid-token
line-wrap footgun).

FORM wraps long output lines in two styles, both of which split TOKENS across
lines, so every consumer must join continuation lines BEFORE tokenizing:

  1. `#write <file> "%E", expr` (and expression dumps generally): the line is
     wrapped at a fixed width and the continuation line starts with
     INDENTATION, mid-token, with no continuation marker:
         FIX=deltalong4444^6+6*gammalong333*deltalong4444^5+15*
               gammalong333^2*deltalong4444^4+...
  2. stdout `Print` output (e.g. Format floatprecision floats): the wrapped
     line ends with a BACKSLASH and the continuation line is indented:
            FORACLE =
               1.20205690315959428539973816151144999076498629234049888179227\
               15553418382e+00;

Both styles pinned on real FORM 5.0.1 output
(tests/fixtures/routing_glue/{wrap_fixture.out,cap_float.stdout}).

`read_form_output` joins both continuation styles into clean logical lines:

  - a line ending in a backslash is ALWAYS continued by the next line
    (leading whitespace of the continuation stripped, backslash dropped);
  - a line STARTING with whitespace is a continuation of the previous
    logical line (leading whitespace stripped, joined with NO separator —
    FORM wraps mid-token, never between two space-separated word tokens),
    UNLESS the previous logical line is empty: blank lines reset wrapping,
    so records that legitimately begin indented after a blank line (FORM's
    own stdout print format: "   NAME =") start a fresh logical line.

Joining is done with no inserted space: in FORM expression output adjacent
tokens are always separated by operator characters, so a mid-token split is
the only case that matters and it must be joined seamlessly.

Consumers: tools/formglue/form_oracle.py (stdout float parse) and
tools/formglue/form_route.py (diagram-generator #write dumps). Regression test
with a GENUINE FORM-5-made wrapping fixture + mutation control:
tools/formglue/tests/test_form_io.py.
"""

from __future__ import annotations

import os

__all__ = ["read_form_output", "unwrap_form_text"]


def unwrap_form_text(text: str) -> str:
    """Join FORM's wrapped/indented continuation lines in `text` (see module
    docstring for the exact rules). Returns text of clean logical lines."""
    logical = []
    for raw in text.splitlines():
        line = raw.rstrip("\r")
        if logical and logical[-1].endswith("\\"):
            # backslash continuation: drop marker, strip continuation indent
            logical[-1] = logical[-1][:-1] + line.lstrip(" \t")
        elif logical and logical[-1] != "" and line[:1] in (" ", "\t") \
                and line.strip() != "":
            # indented wrap continuation of a non-blank logical line
            logical[-1] = logical[-1] + line.lstrip(" \t")
        else:
            logical.append(line)
    return "\n".join(logical)


def read_form_output(path_or_str) -> str:
    """Read FORM output from a file path or a raw string and return it with
    all wrapped continuation lines joined into clean logical lines.

    Input resolution (strict, documented): os.PathLike -> file; a str with no
    newline that names an existing file -> file; anything else -> treated as
    the output text itself.
    """
    if isinstance(path_or_str, os.PathLike):
        with open(os.fspath(path_or_str)) as f:
            return unwrap_form_text(f.read())
    if not isinstance(path_or_str, str):
        raise TypeError(f"path_or_str must be str or PathLike, "
                        f"got {type(path_or_str).__name__}")
    if "\n" not in path_or_str and os.path.isfile(path_or_str):
        with open(path_or_str) as f:
            return unwrap_form_text(f.read())
    return unwrap_form_text(path_or_str)
