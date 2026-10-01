#!/usr/bin/env python3
"""Verify the FORM computation of the gg -> q qbar tree-level squared
amplitude against the textbook result, symbolically.

What it does:
  1. finds `form` on PATH (override with the FORM_BIN environment variable),
  2. runs gg2qq.frm and parses the resulting |M|^2 (summed over quark spins,
     physical gluon polarizations, and colors; d and N symbolic),
  3. averages over initial spins and colors (divide by 4 (N^2-1)^2), sets
     d = 4, N = 3, u = -s - t, and checks exact symbolic equality with the
     textbook reference (Ellis, Stirling & Webber, "QCD and Collider
     Physics"):
         |M|^2_avg / g^4 = (t^2+u^2)/(6 t u) - (3/8) (t^2+u^2)/s^2 ,
  4. mutation self-test: reruns FORM with -D MUT=1 (flipped three-gluon
     vertex sign) and -D MUT=2 (broken color Fierz identity) and requires
     the comparison to FAIL both times.  A gate that cannot fail is not
     a gate.

The comparison is computed from the FORM output with sympy each time this
script runs; the FORM result is not embedded here.  Exit codes:
  0 = all gates passed;  1 = a gate failed;  2 = form or sympy unavailable.
"""
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
FRM = HERE / "gg2qq.frm"


def bail(msg, rc):
    print(msg)
    sys.exit(rc)


def find_form():
    binname = os.environ.get("FORM_BIN", "form")
    if shutil.which(binname) is None:
        bail(
            "SKIP: `%s` is not on PATH, so the FORM run cannot be checked.\n"
            "Install FORM (https://github.com/form-dev/form) or point "
            "FORM_BIN at a FORM binary, then rerun." % binname,
            2,
        )
    return binname


def run_form(formbin, mut=None):
    """Run gg2qq.frm in a scratch dir; return the Msq expression string."""
    with tempfile.TemporaryDirectory() as tmp:
        shutil.copy(FRM, tmp)
        cmd = [formbin]
        if mut is not None:
            cmd += ["-D", "MUT=%d" % mut]
        cmd += ["gg2qq.frm"]
        proc = subprocess.run(
            cmd, cwd=tmp, capture_output=True, text=True, timeout=300
        )
        if proc.returncode != 0:
            bail(
                "FAIL: FORM exited with rc=%d\n%s\n%s"
                % (proc.returncode, proc.stdout, proc.stderr),
                1,
            )
        out = Path(tmp, "gg2qq.out").read_text()
    body = re.fullmatch(r"Msq=(.*);", re.sub(r"\s", "", out))
    if body is None:
        bail("FAIL: could not parse gg2qq.out", 1)
    return body.group(1)


def matches_textbook(body, sp):
    """True iff the FORM output equals the ESW reference exactly."""
    s, t, u, d, NF = sp.symbols("s t u d NF")
    expr = sp.sympify(
        body.replace("^", "**").replace("i_", "I"),
        locals=dict(s=s, t=t, u=u, d=d, NF=NF),
    )
    if expr.has(sp.I):
        return False
    # average over initial state: 2 helicities x (N^2-1) colors per gluon
    avg = expr / (4 * (NF**2 - 1) ** 2)
    got = avg.subs({d: 4, NF: 3, u: -s - t})
    ref = (
        (t**2 + u**2) / (6 * t * u) - sp.Rational(3, 8) * (t**2 + u**2) / s**2
    ).subs(u, -s - t)
    return sp.simplify(got - ref) == 0


def main():
    try:
        import sympy as sp
    except ImportError:
        bail("SKIP: sympy is required for the symbolic comparison "
             "(pip install sympy).", 2)
    formbin = find_form()

    body = run_form(formbin)
    if not matches_textbook(body, sp):
        bail("FAIL: FORM |M|^2 does NOT match the textbook result.", 1)
    print("PASS: FORM tree-level |M|^2 for gg -> q qbar matches "
          "(t^2+u^2)/(6tu) - (3/8)(t^2+u^2)/s^2 exactly (N=3, d=4, u=-s-t).")

    for mut, what in ((1, "three-gluon vertex sign"),
                      (2, "color Fierz identity")):
        if matches_textbook(run_form(formbin, mut), sp):
            bail("FAIL: mutation MUT=%d (%s) was NOT detected -- "
                 "the gate cannot fail." % (mut, what), 1)
        print("PASS: mutation MUT=%d (broken %s) detected -- comparison "
              "failed as required." % (mut, what))

    print("All gates passed.")


if __name__ == "__main__":
    main()
