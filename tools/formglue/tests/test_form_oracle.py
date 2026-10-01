#!/usr/bin/env python3
"""Selftest for tools/formglue/form_oracle.py — STRICT and rc-coded.

Knowns vs mpmath at TWO precisions (50d, 90d), an independent tornheim.py
cross-check, error-path checks (input validation, env-override run/parse
failures), and a MUTATION CONTROL (perturbed reference MUST make the
comparator fail — a comparator that cannot fail is void).

ENGINE GATING: the three legs that need a real float-capable FORM 5 probe the
engine ONCE per session (form_oracle.probe_float_support) and SKIP with a
named reason — which binary/version was found and what capability is missing —
when FORM is absent OR present but float-incapable (e.g. FORM 4.x, which has
no #StartFloat). Never a silent pass, never a spurious fail. The engine-free
legs (input validation, parse strictness, error paths against stub binaries)
always run.

Run modes:
  python3 tools/formglue/tests/test_form_oracle.py   -> prints per-test PASS/SKIP/FAIL,
                                               exits nonzero on ANY failure
  pytest tools/formglue/tests/test_form_oracle.py    -> same tests, pytest-discovered

No ambient-dps mutation anywhere (mp.workdps only). Scratch goes to a temp
dir (override with FORMGLUE_TEST_SCRATCH).
"""

import os
import stat
import sys
import traceback

import mpmath as mp

_HERE = os.path.dirname(os.path.abspath(__file__))
_PKG = os.path.dirname(_HERE)            # the formglue package dir
_TOOLS = os.path.dirname(_PKG)           # the repo tools/ dir
for _p in (_TOOLS, _PKG):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import form_oracle  # noqa: E402
from form_oracle import (FormParseError, FormRunError, euler, mzv,  # noqa: E402
                         zeta)

import tempfile
SCRATCH = os.environ.get("FORMGLUE_TEST_SCRATCH") or os.path.join(
    tempfile.gettempdir(), "formglue_oracle_selftest")
PRECISIONS = (50, 90)
REF_GUARD = 25  # extra dps for references / comparisons


def digits_agree(got_str: str, ref, dps: int) -> float:
    """Decimal digits of agreement between a returned string and a reference,
    computed at dps+REF_GUARD. STRICT: no clipping, no rounding-up."""
    with mp.workdps(dps + REF_GUARD):
        got = mp.mpf(got_str)
        ref = mp.mpf(ref) if not isinstance(ref, (mp.mpf, mp.ctx_mp_python.mpf)) else ref
        d = abs(got - ref)
        if d == 0:
            return float("inf")
        return float(-mp.log10(d / abs(ref)))


def _refs(dps: int):
    """Independent mpmath references at dps+REF_GUARD. All closed forms."""
    with mp.workdps(dps + REF_GUARD):
        return {
            "zeta2": (mp.zeta(2), "mpmath zeta(2)"),
            "zeta3": (mp.zeta(3), "mpmath zeta(3)"),
            "zeta4": (mp.zeta(4), "mpmath zeta(4)"),
            "mzv21": (mp.zeta(3), "Euler identity zeta(2,1)=zeta(3)"),
            "mzv31": (mp.pi ** 4 / 360, "zeta(3,1)=pi^4/360"),
            "mzv22": (mp.pi ** 4 / 120, "zeta(2,2)=pi^4/120"),
            "mzv211": (mp.zeta(4), "zeta(2,1,1)=zeta(4)"),
            "em1": (-mp.log(2), "euler_(-1)=-ln2"),
            "em2": (-mp.pi ** 2 / 12, "euler_(-2)=-pi^2/12"),
            "em21": (mp.zeta(3) / 8, "euler_(-2,1)=zeta(3)/8"),
        }


class _SkipTest(Exception):
    """Script-mode skip signal (pytest mode uses pytest.skip)."""


def _skip(reason: str):
    """Loud NAMED skip in either run mode."""
    if os.environ.get("PYTEST_CURRENT_TEST"):
        import pytest
        pytest.skip(reason)
    raise _SkipTest(reason)


def _require_float_engine():
    """Gate for legs that need a real float-capable FORM 5. Probes the engine
    once per session (cached in form_oracle by binary); on an absent or
    float-incapable FORM (e.g. FORM 4.x: no #StartFloat) this SKIPS with the
    probe's named reason instead of failing."""
    capable, detail = form_oracle.probe_float_support()
    if not capable:
        _skip(detail)


_CALLS = {
    "zeta2": lambda d: zeta(2, d, workdir=SCRATCH),
    "zeta3": lambda d: zeta(3, d, workdir=SCRATCH),
    "zeta4": lambda d: zeta(4, d, workdir=SCRATCH),
    "mzv21": lambda d: mzv((2, 1), d, workdir=SCRATCH),
    "mzv31": lambda d: mzv((3, 1), d, workdir=SCRATCH),
    "mzv22": lambda d: mzv((2, 2), d, workdir=SCRATCH),
    "mzv211": lambda d: mzv((2, 1, 1), d, workdir=SCRATCH),
    "em1": lambda d: euler((-1,), d, workdir=SCRATCH),
    "em2": lambda d: euler((-2,), d, workdir=SCRATCH),
    "em21": lambda d: euler((-2, 1), d, workdir=SCRATCH),
}


# ------------------------------------------------------------------- the tests

def test_knowns_two_precisions():
    """All knowns >= dps digits vs independent closed forms, at 50d AND 90d."""
    _require_float_engine()
    for dps in PRECISIONS:
        refs = _refs(dps)
        for name, call in _CALLS.items():
            got = call(dps)
            ref, prov = refs[name]
            d = digits_agree(got, ref, dps)
            assert d >= dps, (
                f"{name}@{dps}d: only {d:.1f} digits vs {prov} (bar {dps})")


def test_tornheim_independent_cross_check():
    """mzv(5,3) vs tornheim.mzv2(5,3) — independent Richardson-Hurwitz engine,
    dps-keyed cache. One precision (60d) to keep the ref cost ~1 s."""
    _require_float_engine()
    import tornheim
    dps = 60
    got = mzv((5, 3), dps, workdir=SCRATCH)
    with mp.workdps(dps + REF_GUARD):
        ref = tornheim.mzv2(5, 3)
    d = digits_agree(got, ref, dps)
    assert d >= dps, f"mzv(5,3)@{dps}d: only {d:.1f} digits vs tornheim (bar {dps})"


def test_mutation_control():
    """Perturb the zeta(3) reference by 1e-(dps-10): the digit comparator MUST
    fail at both precisions. If it passes, the comparator is broken -> FAIL."""
    _require_float_engine()
    for dps in PRECISIONS:
        got = zeta(3, dps, workdir=SCRATCH)
        with mp.workdps(dps + REF_GUARD):
            bad_ref = mp.zeta(3) + mp.mpf(10) ** (-(dps - 10))
        d = digits_agree(got, bad_ref, dps)
        assert d < dps, (
            f"MUTATION CONTROL FAILED at {dps}d: comparator still reports "
            f"{d:.1f} >= {dps} digits against a perturbed reference")
        # and the measured agreement must sit AT the planted defect size
        assert dps - 13 <= d <= dps - 7, (
            f"mutation@{dps}d: measured {d:.1f} digits, expected ~{dps - 10}")


def test_input_validation():
    for bad in [lambda: zeta(1, 50), lambda: mzv((1, 2), 50),
                lambda: euler((1, 2), 50), lambda: euler((0,), 50),
                lambda: mzv((), 50), lambda: mzv((2, 1), 0),
                lambda: mzv((2, 1), 50.0), lambda: mzv((12, 11), 50)]:
        try:
            bad()
        except ValueError:
            continue
        raise AssertionError(f"{bad} did not raise ValueError")


def test_env_override_and_run_error():
    """FORM5_BIN is read at call time; a broken binary must raise FormRunError."""
    old = os.environ.get("FORM5_BIN")
    os.environ["FORM5_BIN"] = "/bin/false"
    try:
        try:
            zeta(2, 30, workdir=SCRATCH)
        except FormRunError:
            pass
        else:
            raise AssertionError("FORM5_BIN=/bin/false did not raise FormRunError")
    finally:
        if old is None:
            os.environ.pop("FORM5_BIN", None)
        else:
            os.environ["FORM5_BIN"] = old


def test_parse_error_on_garbage_engine():
    """A fake engine emitting a malformed value must raise FormParseError."""
    os.makedirs(SCRATCH, exist_ok=True)
    fake = os.path.join(SCRATCH, "fake_form.sh")
    with open(fake, "w") as f:
        f.write("#!/bin/sh\nprintf '   FORACLE =\\n      garbage+1.0;\\n'\n")
    os.chmod(fake, os.stat(fake).st_mode | stat.S_IEXEC)
    old = os.environ.get("FORM5_BIN")
    os.environ["FORM5_BIN"] = fake
    try:
        try:
            zeta(2, 30, workdir=SCRATCH)
        except FormParseError:
            pass
        else:
            raise AssertionError("garbage engine output did not raise FormParseError")
    finally:
        if old is None:
            os.environ.pop("FORM5_BIN", None)
        else:
            os.environ["FORM5_BIN"] = old


def test_parse_strictness_units():
    """_parse_single unit checks: duplicates, missing, unevaluated leftovers."""
    ok = "   FORACLE =\n      1.20205690315959e+00;\n"
    assert form_oracle._parse_single(ok) == "1.20205690315959e+00"
    for bad in ["", ok + ok,
                "   FORACLE =\n      mzv_(5,3);\n",
                "   FORACLE =\n      1.2020e+00 + mzv_(2);\n"]:
        try:
            form_oracle._parse_single(bad)
        except FormParseError:
            continue
        raise AssertionError(f"parse accepted malformed output: {bad!r}")


TESTS = [test_knowns_two_precisions, test_tornheim_independent_cross_check,
         test_mutation_control, test_input_validation,
         test_env_override_and_run_error, test_parse_error_on_garbage_engine,
         test_parse_strictness_units]


def main() -> int:
    failures = skips = 0
    for t in TESTS:
        try:
            t()
            print(f"[PASS] {t.__name__}")
        except _SkipTest as e:
            skips += 1
            print(f"[SKIP] {t.__name__}: {e}")
        except Exception:
            failures += 1
            print(f"[FAIL] {t.__name__}")
            traceback.print_exc()
    print(f"TOTAL: {len(TESTS)} tests, {failures} failures, "
          f"{skips} named skips")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
