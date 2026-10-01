#!/usr/bin/env python3
"""m2-path gates for the kira adapter CoeffEvaluator (ordering-5 banks
symbolic in (d, m2) take an m2 path GUARDED BY CHARSET ONLY — expressions
without 'm' take the m2-free code path verbatim).

Legs:
  1. m2-free byte-identity: on the T1 Generate coefficient corpus (the
     battery's m2-free corpus of record; archived — leg skips without
     WINNOW_T1_ROOT) the m2-free literal-wrap regex and the m2-path regex
     produce CHAR-IDENTICAL wrapped strings, and evaluator values are
     IDENTICAL with and without an m2 value supplied.
  2. m2-bearing values vs an INDEPENDENT Fraction-arithmetic oracle at
     m2 = 17/5, d = 19999/5000 (mod p), full grammar (* + - / ^).
  3. A small system with rational-m2 coefficients (built FROM coefficient
     strings at m2 = 17/5 mod p) through System -> eliminate -> retrofit
     -> vendored RECEIPT verify; truth from Fraction arithmetic, never from
     the evaluator; one mutated witness MUST FAIL.
  4. Refusals: m2-bearing expr with no m2 value -> ValueError; alien chars
     -> ValueError (both paths).
"""
import _common  # noqa: F401  MUST be first: sys.path bootstrap for the package

import glob
import gzip
import json
import os
import re
from fractions import Fraction

import ibplapper as lap
from ibplapper.adapters.kira import CoeffEvaluator
from ibplapper.receipt.core import load_witness, verify_row
from _common import SLICES4, T1GEN, check, finish, out_dir

OUT = out_dir("/tmp/ibplapper_test_m2_coeffs")
P = 2147483647
D0 = (19999 * pow(5000, P - 2, P)) % P      # d  = 19999/5000 mod p
M20 = (17 * pow(5, P - 2, P)) % P           # m2 = 17/5 mod p
E0 = 87654321

# ---- leg 1: m2-free byte-identity on the T1 corpus ------------------------
OLD_WRAP = lambda e: re.sub(r"(\d+)", r"F(\1)", e.replace("^", "**"))  # noqa: E731
NEW_WRAP = lambda e: re.sub(r"(?<![a-zA-Z_])(\d+)", r"F(\1)",          # noqa: E731
                            e.replace("^", "**"))

if T1GEN is None:
    print("SKIP leg 1: the archived T1 corpus needs WINNOW_T1_ROOT")
else:
    corpus = set()
    for f in sorted(glob.glob(os.path.join(T1GEN, "tmp", "*",
                                           "SYSTEM_*.gz"))):
        lines = gzip.open(f, "rt").read().split("\n")
        i = 0
        while i < len(lines):
            if lines[i] != "Eq":
                i += 1
                continue
            nt = int(lines[i + 2])
            for j in range(nt):
                corpus.add(lines[i + 3 + j].split()[1])
            i += 3 + nt
    check("t1_corpus_harvested", len(corpus) > 10, f"n={len(corpus)}")
    check("t1_corpus_m2_free", all("m" not in e for e in corpus))
    check("wrap_char_identical_on_m2free_corpus",
          all(OLD_WRAP(e) == NEW_WRAP(e) for e in corpus))

    p1, d1, e1 = SLICES4[0]
    ev_plain = CoeffEvaluator(p1, d1, e1)
    ev_with = CoeffEvaluator(p1, d1, e1, m20=12345)
    check("values_identical_with_and_without_m2_env",
          all(ev_plain(e) == ev_with(e) for e in corpus))

# ---- leg 2: m2-bearing values vs independent Fraction oracle --------------
M2_F, D_F = Fraction(17, 5), Fraction(19999, 5000)
CASES = ["2*m2", "-2*m2", "2*m2+3", "m2*m2-4*m2+1", "3*m2/5-d/2",
         "m2^2", "-m2^3+2*d*m2-7", "(2*m2+3)/(5*m2)", "d*m2/(m2-1)",
         "1/m2", "m2/3", "17-m2"]
ev = CoeffEvaluator(P, D0, E0, m20=M20)


def frac_oracle(expr):
    v = eval(expr.replace("^", "**"),                       # noqa: S307
             {"m2": M2_F, "d": D_F, "__builtins__": {}})
    v = Fraction(v)
    return (v.numerator * pow(v.denominator, P - 2, P)) % P


for e in CASES:
    check(f"m2_value_vs_fraction_oracle::{e}", ev(e) == frac_oracle(e))

# compiled-cache path: second call must return the same value
check("m2_cache_stable", ev(CASES[0]) == ev(CASES[0]))

# ---- leg 3: small m2-bearing system end-to-end ----------------------------
A, B, C, E, F_ = "2*m2+3", "m2/3", "5-m2", "d*m2", "7*m2-2"
AF, BF, CF, EF, FF = (Fraction(2) * M2_F + 3, M2_F / 3, 5 - M2_F,
                      D_F * M2_F, 7 * M2_F - 2)
ROWS = [{203: 1, 101: (-ev("2*m2+3")) % P},                 # 203 = A*101
        {202: 1, 203: (-ev("m2/3")) % P, 102: (-ev("5-m2")) % P},
        {201: 1, 202: (-ev("d*m2")) % P, 101: (-ev("7*m2-2")) % P}]


def modp(fr):
    fr = Fraction(fr)
    return (fr.numerator * pow(fr.denominator, P - 2, P)) % P


TRUTH = {203: {101: modp(AF)},
         202: {101: modp(BF * AF), 102: modp(CF)},
         201: {101: modp(EF * BF * AF + FF), 102: modp(EF * CF)}}

sys_ = lap.System(ROWS, P, {201: 3, 202: 2, 203: 1}, forbid={101, 102})
res = lap.eliminate(sys_, lap.Schedule(policy="B2FT"), witnesses="retrofit",
                    witness_targets=[201, 202, 203],
                    witness_dir=os.path.join(OUT, "wits"))
for t in (201, 202, 203):
    rec = res.witnesses[t]
    check(f"m2_system_witness_certified::{t}",
          rec["status"] == "CERTIFIED" and rec["verify_pass"]
          and rec["c"] == TRUTH[t], rec.get("status"))

# mutated witness MUST FAIL the vendored verifier
raw = json.load(open(res.witnesses[202]["witness_path"]))
raw["lam"]["val"][0] = (raw["lam"]["val"][0] + 1) % P
bad = os.path.join(OUT, "w_mut.json")
json.dump(raw, open(bad, "w"))
ok, _ = verify_row(ROWS, load_witness(bad))
check("m2_system_mutated_witness_MUSTFAIL", not ok)

# ---- leg 4: refusals ------------------------------------------------------
try:
    CoeffEvaluator(P, D0, E0)("2*m2+3")
    check("m2_without_value_refused", False)
except ValueError:
    check("m2_without_value_refused", True)
for expr in ("2*m2+x", "m2%3"):
    try:
        ev(expr)
        check(f"alien_chars_refused::{expr}", False)
    except ValueError:
        check(f"alien_chars_refused::{expr}", True)
try:
    CoeffEvaluator(*SLICES4[0])("d+q")
    check("alien_chars_refused_m2free_lane", False)
except ValueError:
    check("alien_chars_refused_m2free_lane", True)

finish("test_m2_coeffs")
