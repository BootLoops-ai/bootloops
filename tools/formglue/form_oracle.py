#!/usr/bin/env python3
"""form_oracle — thin strict wrapper around FORM 5's arbitrary-precision MZV engine.

Exposes zeta(n, dps), mzv(indices, dps), euler(indices, dps): each generates a
.frm program, runs the FORM 5 binary (float_support build), parses the printed
float STRICTLY, and returns an mpmath-ready string carrying >= dps correct
digits. An INDEPENDENT numeric oracle for zeta values / multiple zeta values /
alternating Euler sums — a different engine and algorithm class from mpmath
and tornheim.py, for cross-checking PSLQ boundary constants and ζ-ring fits.

Binary: FORM5_BIN env var, else `form` on PATH — must be FORM 5 built with
float support (FORM 5.0.1, float_support=true, flint=true; build from
github.com/form-dev/form with --enable-float). The env var is read AT CALL TIME
(no import-time freeze). Distro FORM 4.3 has NO float support and
cannot back this wrapper (it remains a good trace cross-check engine).

CONVENTIONS (pinned empirically by this tool's selftest;
validation receipts not shipped):

  mzv_(a1,a2,...,ak) = sum_{m1>m2>...>mk>=1}  m1^-a1 * m2^-a2 * ... * mk^-ak
      OUTER index FIRST. Pin: mzv_(2,1) == zeta(3) (Euler identity), measured
      59.9d/99.9d at 60d/100d; matches tools/formglue/tornheim.py mzv2(s,t) convention
      (mzv_(5,3) vs tornheim.mzv2(5,3): 60.2d/100.0d). Depth 3 pinned by
      mzv_(2,1,1) == zeta(4) (all 60 printed digits).
      Convergence requires a1 >= 2.

  euler_(s1,...,sk): a NEGATIVE index puts (-1)^m on that index's OWN
      summation variable (inner or outer alike).
      euler_(-s)   = sum_m (-1)^m / m^s = -altzeta(s);  euler_(-1) = -ln 2.
      euler_(-2,1) = sum_{m>n>=1} (-1)^m/(m^2 n) = zeta(3)/8   (validated vs
      CVZ-accelerated independent sums, 60.4d/99.7d; inner-sign euler_(2,-1)
      likewise 60.1d/100.6d). Convergence requires s1 != +1 (s1 = -1 is fine).

PRECISION MODEL: FORM's `Format floatprecision` prints P+1 significant digits
at `#StartFloat P d`; measured true agreement vs independent references is ~P
digits (print-width capped; internal precision is GMP limb-rounded above P).
This wrapper runs FORM at P = dps + GUARD (GUARD=10), so the returned string
carries >= dps correct digits. Consume with mp.workdps(dps): mp.mpf(s) —
never rely on ambient mp.mp.dps (dps-lint footgun class).

STRICTNESS: any FORM nonzero rc, missing/duplicated printed value, malformed
float, or unevaluated mzv_/euler_ leftover raises (FormRunError /
FormParseError). Weight (sum of |indices|) is capped at MAX_WEIGHT=22 =
the highest weight actually probed (at 100d); above that raise
rather than extrapolate trust.

NOT exposed (FORM 5.0.1 limitations): lin_/hpl_/mpl_ are
reserved names with NO numerical evaluation — do not lean on them for HPL
oracles; mzvhalf_ exists but is unexercised here.

Selftest: tests/test_form_oracle.py (strict, rc-coded, knowns vs mpmath
at two precisions + mutation control). Engine capability is probed once per
session via probe_float_support(): tests that need the float engine SKIP with
a named reason (binary + version found, capability missing) when the FORM on
hand is absent or float-incapable — never a silent pass, never a spurious fail.
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
import tempfile

# Package-relative import — never routed through the
# flat tools/ shims. A dirname(__file__) sys.path.insert is deliberately
# absent: inside the package it would alias formglue modules as top-level
# duplicates. Script mode (python3 tools/formglue/form_oracle.py) anchors
# the package PARENT (tools/) instead.
try:
    from .form_io import read_form_output  # noqa: E402  (shared #write/print unwrap fix)
except ImportError:  # direct script execution: no parent package
    _PKG_PARENT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if _PKG_PARENT not in sys.path:
        sys.path.insert(0, _PKG_PARENT)
    from formglue.form_io import read_form_output  # noqa: E402

FORM5_BIN_DEFAULT = "form"   # PATH lookup; override with env FORM5_BIN
GUARD = 10          # internal FORM digits above requested dps
MAX_WEIGHT = 22     # highest weight probed; raise only after probing higher
MAX_DPS = 1000      # sanity ceiling; raise deliberately if ever needed
_EXPR = "FORACLE"   # expression name used in generated .frm programs

# strict FORM float print: mantissa d.ddd... with exponent e[+-]dd
_FLOAT_RE = re.compile(r"^[-+]?\d+\.\d+e[+-]\d+$")


class FormOracleError(RuntimeError):
    """Base class for form_oracle failures."""


class FormRunError(FormOracleError):
    """FORM exited nonzero / could not be executed."""


class FormParseError(FormOracleError):
    """FORM output did not contain exactly one well-formed value."""


def _form_bin() -> str:
    # read at CALL time — overridable per-call via env, no import-time freeze
    return os.environ.get("FORM5_BIN", FORM5_BIN_DEFAULT)


def _check_dps(dps) -> int:
    if not isinstance(dps, int) or isinstance(dps, bool):
        raise ValueError(f"dps must be an int, got {dps!r}")
    if not (1 <= dps <= MAX_DPS):
        raise ValueError(f"dps={dps} outside [1, {MAX_DPS}]")
    return dps


def _check_indices(indices, kind: str):
    if isinstance(indices, int):
        indices = (indices,)
    indices = tuple(indices)
    if not indices:
        raise ValueError(f"{kind}: empty index tuple")
    for a in indices:
        if not isinstance(a, int) or isinstance(a, bool):
            raise ValueError(f"{kind}: non-integer index {a!r} in {indices}")
    if kind == "mzv":
        if any(a < 1 for a in indices):
            raise ValueError(f"mzv: indices must be >= 1, got {indices}")
        if indices[0] < 2:
            raise ValueError(f"mzv: leading index must be >= 2 (divergent), got {indices}")
    else:  # euler
        if any(a == 0 for a in indices):
            raise ValueError(f"euler: zero index in {indices}")
        if indices[0] == 1:
            raise ValueError(f"euler: leading index +1 diverges, got {indices}")
    weight = sum(abs(a) for a in indices)
    if weight > MAX_WEIGHT:
        raise ValueError(
            f"{kind}: weight {weight} > MAX_WEIGHT={MAX_WEIGHT} "
            f"(highest weight probed; raise deliberately after probing)")
    return indices, weight


def _run_form(frm_text: str, workdir=None, timeout: float = 300.0) -> str:
    """Write frm_text to a .frm in workdir (fresh tempdir if None), run FORM 5,
    return stdout. Raises FormRunError on any nonzero rc / exec failure."""
    binary = _form_bin()

    def _go(wd: str) -> str:
        frm = os.path.join(wd, "form_oracle_job.frm")
        with open(frm, "w") as f:
            f.write(frm_text)
        try:
            r = subprocess.run([binary, "form_oracle_job.frm"], cwd=wd,
                               capture_output=True, text=True, timeout=timeout)
        except FileNotFoundError:
            # `from None`: the chained stdlib traceback only obscures the
            # actual problem (no engine) — keep the report to one named line
            raise FormRunError(
                f"FORM binary not found ({binary}): install a float-capable "
                f"FORM 5 or point FORM5_BIN at one (see GUIDE.md)") from None
        except (OSError, subprocess.TimeoutExpired) as e:
            raise FormRunError(f"FORM launch failed ({binary}): {e}") from e
        if r.returncode != 0:
            raise FormRunError(
                f"FORM rc={r.returncode} ({binary})\n--- frm ---\n{frm_text}"
                f"\n--- stdout ---\n{r.stdout}\n--- stderr ---\n{r.stderr}")
        return r.stdout

    if workdir is not None:
        os.makedirs(workdir, exist_ok=True)
        return _go(workdir)
    with tempfile.TemporaryDirectory(prefix="form_oracle_") as wd:
        return _go(wd)


def _parse_single(stdout: str, name: str = _EXPR) -> str:
    """Extract exactly one printed float value for expression `name`. STRICT."""
    # join FORM's wrapped lines (backslash AND indented-continuation styles)
    # via the shared fix in tools/formglue/form_io.py
    txt = read_form_output(stdout)
    hits = re.findall(
        r"^\s+%s\s*=\s*\n?\s*([^;]+);" % re.escape(name), txt, flags=re.M)
    # duplicate detection must ALSO run on the RAW text: FORM's wrap format
    # cannot distinguish a duplicated record glued without blank-line
    # separators from an indented wrap continuation, so count the raw
    # record-start markers as well (strictness preserved after the unwrap)
    raw_marks = re.findall(
        r"^\s+%s\s*=" % re.escape(name), stdout, flags=re.M)
    if len(hits) != 1 or len(raw_marks) > 1:
        raise FormParseError(
            f"expected exactly 1 printed value for {name}, found "
            f"{max(len(hits), len(raw_marks))}\n--- stdout ---\n{stdout}")
    val = "".join(hits[0].split())
    if "mzv_" in val or "euler_" in val:
        raise FormParseError(
            f"{name} left (partially) UNevaluated: {val!r}\n--- stdout ---\n{stdout}")
    if not _FLOAT_RE.match(val):
        raise FormParseError(
            f"{name} value is not a well-formed FORM float: {val!r}"
            f"\n--- stdout ---\n{stdout}")
    return val


# ------------------------------------------------------------ capability probe

_PROBE_CACHE: dict = {}   # binary -> (capable: bool, detail: str), per session
_PROBE_REF = "1.64493406684822"   # zeta(2) to 15 digits — well inside the 20d probe


def probe_float_support(binary: str | None = None, timeout: float = 60.0):
    """Cheap one-shot capability probe: can this FORM evaluate floats?

    Runs a 1-line #StartFloat/Evaluate program (zeta(2) at 20 digits, ~ms) and
    checks the printed value. Returns (capable, detail); detail names the
    binary, the version banner found, and — when not capable — exactly what is
    missing. Cached per binary for the session. Never raises: absence and
    incapability are both reported as (False, named reason).

    This is a helper for callers (e.g. the selftest) to decide run-vs-skip;
    zeta()/mzv()/euler() themselves stay strict and raise on any engine
    problem.
    """
    binary = binary or _form_bin()
    hit = _PROBE_CACHE.get(binary)
    if hit is not None:
        return hit
    frm = ("#StartFloat 20d, MZV=2\n"
           "Format floatprecision;\n"
           f"L {_EXPR} = mzv_(2);\n"
           "Evaluate;\n"
           "Print +f;\n"
           ".end\n")
    verdict = None
    r = None
    try:
        with tempfile.TemporaryDirectory(prefix="form_probe_") as wd:
            with open(os.path.join(wd, "float_probe.frm"), "w") as f:
                f.write(frm)
            r = subprocess.run([binary, "float_probe.frm"], cwd=wd,
                               capture_output=True, text=True, timeout=timeout)
    except FileNotFoundError:
        verdict = (False,
                   f"no FORM binary found ({binary}): install a float-capable "
                   f"FORM 5 or point FORM5_BIN at one (see GUIDE.md)")
    except (OSError, subprocess.TimeoutExpired) as e:
        verdict = (False, f"FORM float probe could not run ({binary}): {e}")
    if verdict is None:
        m = re.search(r"FORM\s+[0-9][^\s]*", r.stdout)
        ver = m.group(0) if m else "FORM (version banner not recognized)"
        if os.sep not in binary:      # bare PATH name -> show what it resolved to
            import shutil
            binary = shutil.which(binary) or binary
        if r.returncode != 0:
            verdict = (False,
                       f"{ver} at {binary} has no float support "
                       f"(#StartFloat rejected, rc={r.returncode}); the float "
                       f"legs need FORM 5 built with --enable-float — install "
                       f"one or point FORM5_BIN at it (see GUIDE.md)")
        else:
            try:
                ok = _parse_single(r.stdout).startswith(_PROBE_REF)
            except FormParseError:
                ok = False
            if ok:
                verdict = (True, f"{ver} at {binary}: float-capable "
                                 f"(probe zeta(2)@20d verified)")
            else:
                verdict = (False,
                           f"{ver} at {binary}: float probe did not return "
                           f"zeta(2) (unevaluated or malformed output); the "
                           f"float legs need FORM 5 built with --enable-float "
                           f"(see GUIDE.md)")
    _PROBE_CACHE[binary] = verdict
    return verdict


def _eval_builtin(func: str, indices, dps: int, workdir=None) -> str:
    dps = _check_dps(dps)
    indices, weight = _check_indices(indices, func)
    P = dps + GUARD
    args = ",".join(str(a) for a in indices)
    frm = (f"#StartFloat {P}d, MZV={max(weight, 2)}\n"
           f"Format floatprecision;\n"
           f"L {_EXPR} = {func}_({args});\n"
           f"Evaluate;\n"
           f"Print +f;\n"
           f".end\n")
    return _parse_single(_run_form(frm, workdir=workdir))


# ------------------------------------------------------------------ public API

def mzv(indices, dps: int, workdir=None) -> str:
    """Multiple zeta value zeta(a1,...,ak), OUTER index first (see module
    docstring pin). Returns an mpmath-ready string with >= dps correct digits."""
    return _eval_builtin("mzv", indices, dps, workdir=workdir)


def zeta(n: int, dps: int, workdir=None) -> str:
    """Riemann zeta(n), integer n >= 2. Returns mpmath-ready string."""
    return mzv((n,), dps, workdir=workdir)


def euler(indices, dps: int, workdir=None) -> str:
    """Alternating Euler sum euler_(s1,...,sk): negative index = (-1)^m on that
    index's own summation variable (see module docstring pin). Returns an
    mpmath-ready string with >= dps correct digits."""
    return _eval_builtin("euler", indices, dps, workdir=workdir)


def to_mpf(s: str, dps: int):
    """Convenience: parse a returned string at an EXPLICIT dps (never ambient)."""
    import mpmath as mp
    with mp.workdps(_check_dps(dps)):
        return mp.mpf(s)


# ------------------------------------------------------------------------- CLI

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="FORM 5 zeta/MZV/euler-sum oracle (strict wrapper)",
        epilog="negative euler indices: put --dps BEFORE '--', e.g. "
               "form_oracle.py --dps 40 euler -- -2,1")
    ap.add_argument("func", choices=["zeta", "mzv", "euler"])
    ap.add_argument("indices", help="comma-separated integers, e.g. 5,3 or -2,1")
    ap.add_argument("--dps", type=int, required=True)
    args = ap.parse_args(argv)
    idx = tuple(int(t) for t in args.indices.split(","))
    if args.func == "zeta":
        if len(idx) != 1:
            ap.error("zeta takes exactly one index")
        print(zeta(idx[0], args.dps))
    elif args.func == "mzv":
        print(mzv(idx, args.dps))
    else:
        print(euler(idx, args.dps))
    return 0


if __name__ == "__main__":
    sys.exit(main())
