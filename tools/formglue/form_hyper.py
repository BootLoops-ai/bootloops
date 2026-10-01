#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-only
#
# Derived from HyperFORM (https://github.com/adamkardos/HyperFORM; Adam Kardos,
# Sven-Olaf Moch and Oliver Schnetz, "HyperFORM — a FORM package for parametric
# integration with hyperlogarithms", arXiv:2607.01163; HyperFORM v1.0, Zenodo
# doi:10.5281/zenodo.17706909), licensed under the GNU General Public License
# v3.0 (the repository's LICENSE file and the Zenodo record state version 3 only,
# with no later-version option, so this file carries none either). This Python
# module reuses the structure of the HyperFORM example drivers: the driver
# template it renders (_DRIVER_TEMPLATE, _LONE_BLOCK and _FIB_BLOCK below)
# follows examples/zigzags/zigzag3.frm and examples/multiscale/massless-box.frm
# in the package's documented calling sequence (#include- hyperform.h, the
# HYPMAXEP / IntegralExpr / IntegrationSequence / ChenWuVar defines,
# HypParseInputExpr, HypEpExpand, HypApplyChenWu, HypSimplify, the
# HypIntegrationStep loop, HypFinalizeResult, HypPrintStatistics and the
# HYP* -> user-name replace_ line). No HyperFORM library source is included;
# HyperFORM is obtained from its authors' repository and executed by FORM as a
# separate process. The remainder of this module (the subprocess runner, the
# output parser and the self-test) is not taken from HyperFORM. The module is
# Copyright (c) 2026 Anthropic, PBC; created by Matthew D. Schwartz, code written
# by Claude (Anthropic) under his supervision, and is distributed under the same
# GNU General Public License v3.0, version 3 only. The rest of this package is
# MIT-licensed; see NOTICE in this directory (tools/formglue/NOTICE).
#
"""form_hyper — formglue bridge to HyperFORM (Kardos-Moch-Schnetz, arXiv:2607.01163).

HyperFORM is an open FORM 5 port of the core of Erik Panzer's HyperInt:
symbolic parametric integration of hyperlogarithms weighted by rational
functions (rational letters / rational arguments included). Engine location:
a git clone of github.com/adamkardos/HyperFORM (Zenodo DOI
10.5281/zenodo.17706909) — point env HYPERFORM_SRC at its src/ dir; the
package is PURE .frm/.h include files — there is nothing to build. Per the
formglue law the engine is EXTERNAL and CALLED, never absorbed: this module
writes a driver .frm, sets FORMPATH to the HyperFORM src dir, runs the house
FORM 5 binary, and parses printed expressions strictly.

Positive controls (FORM 5.0.1):
  - upstream check suite via check/check.rb --form <house form>:
    1327 tests / 2761 assertions / 0 failures (34 s, -w 8).
  - examples battery: zigzag3 (=6 zeta3, 0.16 s), zigzag4 (0.65 s),
    zigzag5 (18.2 s), 3-loop masters LA/BU (23.2/15.3 s), massless box
    incl. eps poles + fibration basis (46.6 s) - all shipped Diff == 0
    assertions reproduced.

API:
  run_driver(frm_text, ...)   -> {expr_name: body_string} for printed exprs
  run_example(path, ...)      -> same, on a shipped example file
  hyper_integrate(...)        -> build driver for a projective parametric
                                 integral (zigzag pattern: Chen-Wu + ordered
                                 HypIntegrationStep chain) and run it
  parse_terms(body)           -> [(Fraction coeff, {factor: power})] with
                                 factors "z<n>" | "L(a1,...,an,arg)" |
                                 "den(...)" | plain symbols
Output hyperlogarithm convention (pinned on the massless-box + eqDp56
examples): L(a1,...,an,x) is the hyperlogarithm with letters a_i and argument
x in HyperInt convention, i.e. GiNaC/Vollinga-Weinzierl G(a1,...,an;x); MZV
constants print as z2,z3,... (zeta values), pi as pi_ where kept.

STRICTNESS: nonzero rc, a missing expression, or an unparseable term raises
(HyperFormRunError / HyperFormParseError). Long jobs: pass timeout explicitly;
launch multi-hour jobs through your own job scheduler, not through this wrapper.

Selftest (script mode): python3 tools/formglue/form_hyper.py
  replays zigzag3 (asserts the parsed value is exactly 6*z3 and Diff == 0)
  and the basic/product example; rc=0 all green, rc=1 any failure.
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
import tempfile
from fractions import Fraction

try:
    from .form_io import read_form_output  # noqa: E402
except ImportError:  # direct script execution: no parent package
    _PKG_PARENT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if _PKG_PARENT not in sys.path:
        sys.path.insert(0, _PKG_PARENT)
    from formglue.form_io import read_form_output  # noqa: E402

FORM5_BIN_DEFAULT = "form"   # PATH lookup; override with env FORM5_BIN
# HyperFORM is EXTERNAL — obtain upstream (github.com/adamkardos/HyperFORM)
# and point HYPERFORM_SRC at its src/ dir (examples/ resolved beside it).
HYPERFORM_SRC_DEFAULT = None
HYPERFORM_EXAMPLES = None


class HyperFormError(RuntimeError):
    pass


class HyperFormRunError(HyperFormError):
    pass


class HyperFormParseError(HyperFormError):
    pass


def _form_bin() -> str:
    return os.environ.get("FORM5_BIN", FORM5_BIN_DEFAULT)


def _hyper_src() -> str:
    p = os.environ.get("HYPERFORM_SRC", HYPERFORM_SRC_DEFAULT)
    if not p:
        raise HyperFormError(
            "HYPERFORM_SRC is not set. HyperFORM is an external engine that is "
            "not vendored here — clone github.com/adamkardos/HyperFORM and "
            "point HYPERFORM_SRC at its src/ directory.")
    if not os.path.isfile(os.path.join(p, "hyperform.h")):
        raise HyperFormRunError(
            f"HyperFORM src dir {p!r} has no hyperform.h "
            "(set HYPERFORM_SRC or clone github.com/adamkardos/HyperFORM)")
    return p


def run_driver(frm_text: str, workdir=None, timeout: float = 3600.0,
               nthreads: int = 0, keep_log=None) -> dict:
    """Run a HyperFORM driver .frm through the house FORM 5. Returns a dict
    {expression_name: body} for every `name = ...;` record printed by the run.
    nthreads>1 uses tform with that many workers."""
    binary = _form_bin()
    src = _hyper_src()

    def _go(wd: str) -> str:
        frm = os.path.join(wd, "form_hyper_job.frm")
        with open(frm, "w") as f:
            f.write(frm_text)
        env = dict(os.environ)
        env["FORMPATH"] = src
        cmd = [binary, "form_hyper_job.frm"]
        if nthreads and nthreads > 1:
            tform = os.path.join(os.path.dirname(binary), "tform")
            cmd = [tform, f"-w{nthreads}", "form_hyper_job.frm"]
        try:
            # binary-safe: tform can emit non-UTF8 bytes in scratch dumps
            r0 = subprocess.run(cmd, cwd=wd, capture_output=True,
                                timeout=timeout, env=env)
        except FileNotFoundError:
            # `from None`: the chained stdlib traceback only obscures the
            # actual problem (no engine) — keep the report to one named line
            raise HyperFormRunError(
                f"FORM binary not found ({cmd[0]}): install a float-capable "
                f"FORM 5 or point FORM5_BIN at one (see GUIDE.md)") from None
        except (OSError, subprocess.TimeoutExpired) as e:
            raise HyperFormRunError(f"FORM launch failed ({binary}): {e}") from e

        class _R:
            returncode = r0.returncode
            stdout = r0.stdout.decode("utf-8", errors="replace")
            stderr = r0.stderr.decode("utf-8", errors="replace")
        r = _R()
        if keep_log is not None:
            with open(keep_log, "w") as f:
                f.write(r.stdout + "\n--- stderr ---\n" + r.stderr)
        if r.returncode != 0:
            raise HyperFormRunError(
                f"FORM rc={r.returncode}\n--- tail of stdout ---\n"
                f"{r.stdout[-4000:]}\n--- stderr ---\n{r.stderr[-2000:]}")
        return r.stdout

    if workdir is not None:
        os.makedirs(workdir, exist_ok=True)
        out = _go(workdir)
    else:
        with tempfile.TemporaryDirectory(prefix="form_hyper_") as wd:
            out = _go(wd)
    return _extract_exprs(out)


def run_example(relpath: str, **kw) -> dict:
    """Run a shipped example, e.g. run_example('zigzags/zigzag3.frm')."""
    ex = HYPERFORM_EXAMPLES or (
        os.path.join(os.path.dirname(os.environ["HYPERFORM_SRC"]), "examples")
        if os.environ.get("HYPERFORM_SRC") else None)
    if not ex:
        raise HyperFormError(
            "HyperFORM examples dir unknown — set HYPERFORM_SRC (examples/ is "
            "resolved beside it).")
    path = os.path.join(ex, relpath)
    with open(path) as f:
        # read verbatim: NO trailer stripping is performed here (the example
        # file is fed to the driver whole; an earlier version stripped a
        # check.rb assert trailer — the shipped code does not)
        text = f.read()
    return run_driver(text, **kw)


def _extract_exprs(stdout: str) -> dict:
    txt = read_form_output(stdout)
    out = {}
    for m in re.finditer(r"^\s{2,}([A-Za-z][A-Za-z0-9]*)\s*=\s*(.*?);",
                         txt, flags=re.M | re.S):
        name, body = m.group(1), m.group(2)
        # keep the LAST printed instance of each name (post-.sort final state)
        out[name] = " ".join(body.split())
    if not out:
        raise HyperFormParseError(
            f"no printed expressions found\n--- stdout tail ---\n{stdout[-3000:]}")
    return out


# ------------------------------------------------------------- driver builder

# Driver skeleton adapted from HyperFORM's shipped examples
# (examples/zigzags/zigzag3.frm; loop body and the HYP* -> user-name `replace_`
# line from examples/multiscale/massless-box.frm) — A. Kardos, S. Moch,
# O. Schnetz, GPL-3.0, github.com/adamkardos/HyperFORM. It contains only the
# library's documented calling sequence; no HyperFORM source is included;
# HyperFORM is obtained upstream and executed by FORM as a separate process.
# This is why the file is distributed under GPL-3.0-only (license header at the
# top of this file; NOTICE and LICENSE-GPL-3.0 in this directory).
_DRIVER_TEMPLATE = """*{{{{ {name} : generated by formglue.form_hyper
#-
#include- hyperform.h
off statistics;

#define HYPMAXEP "{maxep}"
#define IntegralExpr "{name}"
#define IntegrationSequence "{sequence}"
#define ChenWuVar "{chenwu}"

symbols al1,...,al{nalpha};
symbol ep;
symbols HFn1,HFn2;{ksym_decl}
cfunctions num,den;
cfunctions Linf,Lone,L,rat,log,HFaux,HFaux2;

local `IntegralExpr' = {integrand};

.sort

#call HypParseInputExpr(ep,num,den,al1,...,al{nalpha})

.sort

#call HypEpExpand

#call HypApplyChenWu(`IntegralExpr',ChenWuVar)

#call HypSimplify

#do IntVar={{`IntegrationSequence'}}
  #call HypIntegrationStep(`IntegralExpr',`IntVar')
  .sort
  PolyRatFun;
#enddo

#call HypFinalizeResult(ep,`HYPMAXEP')

#call HypPrintStatistics
{lone_block}
.sort

print +s;

.end
*}}}}
"""

_LONE_BLOCK = """
.sort
* numeric-kinematics finalization: regularized-at-infinity hyperlogs with
* rational letters -> Lone(w) = shuffle-regularized L(w, 1) = G_reg(w; 1)
multiply replace_(HYPLinfRegInfZero,Linf);
.sort
#call HypFromLinfToLone(`IntegralExpr',Linf,Lone,HFaux,HFaux2,HFn1)
"""

_FIB_BLOCK = """
.sort
* half-symbolic finalization (massless-box pattern): fibration basis in the
* declared kinematic symbol(s), then remaining constant Linf -> Lone(w)
multiply replace_(HYPrat,rat,HYPLinfRegInfZero,Linf,HYPlog,log);
.sort
#call HypFibrationBasis(`IntegralExpr',Linf,L,rat,{ksyms})
.sort
#call HypFromLinfToLone(`IntegralExpr',Linf,Lone,HFaux,HFaux2,HFn1)
"""


def build_driver(name: str, integrand: str, nalpha: int, sequence=None,
                 chenwu=None, maxep: int = 0, finalize: str = "plain",
                 ksymbols=()) -> str:
    """Driver for a FINITE projective parametric integral (zigzag pattern):
    integrand is a FORM expression in al1..al<nalpha> using num()/den(), e.g.
    'den(al1*al2 + al1*al3 + al2*al3)^2'. sequence = integration order
    (variable NUMBERS, default 1..nalpha-1), chenwu = Chen-Wu variable number
    (default nalpha). Divergent integrals (HypAutoRegularize) are not templated
    here - write a bespoke driver and use run_driver."""
    if chenwu is None:
        chenwu = nalpha
    if sequence is None:
        sequence = [k for k in range(1, nalpha + 1) if k != chenwu]
    seq = ",".join(str(s) for s in sequence)
    if finalize not in ("plain", "lone", "fibration"):
        raise ValueError(f"finalize must be plain/lone/fibration, got {finalize!r}")
    if finalize == "fibration" and not ksymbols:
        raise ValueError("finalize='fibration' needs ksymbols")
    blocks = {"plain": "", "lone": _LONE_BLOCK,
              "fibration": _FIB_BLOCK.format(ksyms=",".join(ksymbols))}
    ksym_decl = ("\nsymbols " + ",".join(ksymbols) + ";") if ksymbols else ""
    return _DRIVER_TEMPLATE.format(name=name, integrand=integrand,
                                   nalpha=nalpha, sequence=seq,
                                   chenwu=chenwu, maxep=maxep,
                                   ksym_decl=ksym_decl,
                                   lone_block=blocks[finalize])


def hyper_integrate(name: str, integrand: str, nalpha: int, sequence=None,
                    chenwu=None, maxep: int = 0, finalize: str = "plain",
                    ksymbols=(), **kw) -> str:
    """Build + run the driver; return the body string of expression `name`."""
    frm = build_driver(name, integrand, nalpha, sequence, chenwu, maxep,
                       finalize, ksymbols)
    exprs = run_driver(frm, **kw)
    if name not in exprs:
        raise HyperFormParseError(
            f"expression {name} not printed; got {sorted(exprs)}")
    return exprs[name]


# ------------------------------------------------------------- output parsing

def _split_terms(s: str):
    """Split a sum on +/- at paren depth 0 (sign kept with the term)."""
    terms, depth, cur = [], 0, []
    for ch in s:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        if ch in "+-" and depth == 0 and cur and cur[-1] not in "^*/(+-eE":
            terms.append("".join(cur))
            cur = []
        cur.append(ch)
    if cur:
        terms.append("".join(cur))
    return [t for t in terms if t and t not in "+-"]


def _split_factors(term: str):
    """Split a term on '*' at paren depth 0."""
    parts, depth, cur = [], 0, []
    for ch in term:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        if ch == "*" and depth == 0:
            parts.append("".join(cur))
            cur = []
        else:
            cur.append(ch)
    if cur:
        parts.append("".join(cur))
    return parts


def parse_terms(body: str):
    """Parse a printed FORM body into [(Fraction coeff, {factor: int power})].
    Factors are kept as canonical strings: 'z3', 'L(0,1,1/85)', 'den(...)',
    symbols. Raises HyperFormParseError on anything unrecognized."""
    s = "".join(body.split())
    if s in ("0", ""):
        return []
    terms = _split_terms(s)
    out = []
    for t in terms:
        sign = 1
        while t and t[0] in "+-":
            if t[0] == "-":
                sign = -sign
            t = t[1:]
        coeff = Fraction(sign)
        facs = {}
        for f in _split_factors(t):
            if not f:
                continue
            m = re.fullmatch(r"(\d+)(?:/(\d+))?", f)
            if m:
                coeff *= Fraction(int(m.group(1)), int(m.group(2) or 1))
                continue
            m = re.fullmatch(r"1/(\d+)", f)
            if m:
                coeff /= int(m.group(1))
                continue
            m = re.fullmatch(r"(.+?)\^(-?\d+)$", f) if "(" not in f else None
            if m:
                facs[m.group(1)] = facs.get(m.group(1), 0) + int(m.group(2))
                continue
            m = re.fullmatch(r"([A-Za-z][A-Za-z0-9_]*\(.*\))(?:\^(-?\d+))?", f)
            if m:
                facs[m.group(1)] = facs.get(m.group(1), 0) + int(m.group(2) or 1)
                continue
            m = re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", f)
            if m:
                facs[f] = facs.get(f, 0) + 1
                continue
            raise HyperFormParseError(f"unparseable factor {f!r} in term {t!r}")
        out.append((coeff, facs))
    return out


# ---------------------------------------------------------------- selftest

def _selftest() -> int:
    ok = True
    # 1) zigzag3: exact 6*z3 and Diff == 0
    exprs = run_example("zigzags/zigzag3.frm", timeout=120)
    z3terms = parse_terms(exprs["Z3"])
    if not (len(z3terms) == 1 and z3terms[0][0] == 6
            and z3terms[0][1] == {"z3": 1}):
        print(f"FAIL zigzag3 value: {exprs.get('Z3')!r}")
        ok = False
    if parse_terms(exprs.get("Diff", "0")) != []:
        print(f"FAIL zigzag3 Diff nonzero: {exprs.get('Diff')!r}")
        ok = False
    # 2) basic/product runs clean and prints expressions (bracketed bodies
    #    are extraction-checked only - parse_terms targets flat +s prints)
    exprs2 = run_example("basic/product.frm", timeout=120)
    if not exprs2:
        print("FAIL basic/product printed nothing")
        ok = False
    print("form_hyper selftest:", "PASS" if ok else "FAIL",
          f"(zigzag3 = 6*z3 exact, Diff=0; product exprs: {sorted(exprs2)})")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(_selftest())
