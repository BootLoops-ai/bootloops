#!/usr/bin/env python3
"""make_regression — build the fibrate regression corpus fixtures.

SOURCE (cited, gated): the F^red1 pilot's differential-fibration memo
  engine_memo.pkl (env SUBTROPICA_PILOT_DIR points at it)
(370/370 words fibrated + validated to ~1e-90 at ze=2/5 and 3/8 against
the trusted step-1 ZIP semantics; the assembled formula passed the
38-42-digit held-out oracle gate).

This script COPIES the pilot's validated representations for 10 chosen
words (weights 1-5, with and without PSLQ constants) into
  regression_words.json     (CLI input)
  regression_expected.json  (expected canonical terms per word)
using fibrate.serialize_rep — the SAME canonical serializer the CLI uses,
so the Julia test (test/test_fibrate.jl) can compare exactly.

Run (only needed to regenerate; fixtures are checked in — requires
SUBTROPICA_PILOT_DIR pointing at the pilot run dir, which is not shipped):
  ulimit -v 32505856; SUBTROPICA_PILOT_DIR=... python3 make_regression.py
"""
import json
import os
import pickle
import sys

PILOT = os.environ.get("SUBTROPICA_PILOT_DIR")
if not PILOT:
    raise SystemExit("SKIP: SUBTROPICA_PILOT_DIR is not set (pilot run dir with "
                     "engine_memo.pkl; regeneration input, not shipped)")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))          # scripts/fibrate

from fibrate import parse_letter, serialize_rep    # noqa: E402

# 10 words from the pilot's validated 370 (unres_words.json), spanning
# weights 1-5; five with PSLQ-identified constants in their rep, five
# without.
WORDS = [
    ["-ze"],                                            # w1
    ["ze - 1"],                                         # w1
    ["1/(ze - 1)", "0"],                                # w2, consts
    ["-1", "-ze"],                                      # w2
    ["ze - 1", "0", "-1"],                              # w3, consts
    ["-1", "-ze", "-ze"],                               # w3, consts
    ["-1", "-1", "0", "-ze"],                           # w4, consts
    ["1/(ze - 1)", "-1", "1/(ze - 1)", "-1"],           # w4
    ["0", "-1", "-ze", "-ze", "-1"],                    # w5, consts
    ["-1", "-ze", "-1", "-1", "-1"],                    # w5, consts
]


def main():
    with open(os.path.join(PILOT, "engine_memo.pkl"), "rb") as f:
        memo = pickle.load(f)
    expected = []
    for wstr in WORDS:
        key = tuple(parse_letter(x) for x in wstr)
        if key not in memo:
            raise SystemExit(f"word {wstr} not in pilot memo — wrong corpus")
        rep = memo[key]
        used_pslq_own = any(K for (u, K), c in rep.items())
        expected.append({
            "word": wstr,
            "weight": len(wstr),
            "terms": serialize_rep(rep),
            # NOTE: pslq usage at the WORD level (own rep contains named
            # constants).  The CLI's used_pslq flag is TRANSITIVE (includes
            # rational constants fitted for subwords), so the test asserts
            # used_pslq >= this flag, and exact term equality regardless.
            "own_rep_has_named_constants": used_pslq_own,
        })
    src = {
        "source": "engine_memo.pkl (stored pilot memo)",
        "source_gates": "pilot gates (held-out oracle 38-42d; "
                        "per-word 2-point validation ~1e-90)",
        "generated_by": "make_regression.py",
    }
    with open(os.path.join(HERE, "regression_words.json"), "w") as f:
        json.dump(WORDS, f, indent=1)
        f.write("\n")
    with open(os.path.join(HERE, "regression_expected.json"), "w") as f:
        json.dump({**src, "words": expected}, f, indent=1)
        f.write("\n")
    print(f"wrote {len(WORDS)} regression words + expected terms")


if __name__ == "__main__":
    main()
