#!/usr/bin/env python3
"""mplll_lattice.py — integer-lattice LLL/BKZ backends for mplll.py.

Backends (in preference order for BKZ):
  1. fpylll (if importable)
  2. fplll CLI (`/usr/bin/fplll -a bkz -b β`; GMP big-int safe)
  3. Nemo.lll via a Julia subprocess (LLL only; flint big-int)
  4. pure-Python _int_lll (last resort)

All backends take/return list[list[int]] (rows = lattice vectors).

Per-reduction wall (the subprocess backends fplll-cli and Nemo.lll): env
MPLLL_FPLLL_TIMEOUT, read by _fplll_timeout() at every reduction — unset or 0
= UNBOUNDED (the default: a reduction ends by convergence or by its own
failure, never by this clock); a positive integer = seconds, an opt-in the
caller sets deliberately and states in its receipt; anything else is refused
by name (ValueError).  A reduction ended by the knob raises
MplllTimeoutExpired (a subprocess.TimeoutExpired, so existing handlers keep
catching it) whose text reads 'ended by MPLLL_FPLLL_TIMEOUT=<N> s, not by
convergence', and the same line goes to stderr, so the caller's receipt
carries it."""
import os, re, sys, json, math, shutil, tempfile, subprocess, time

try:
    from fpylll import IntegerMatrix, LLL as _fpylll_LLL
    _HAVE_FPYLLL = True
except Exception:
    _HAVE_FPYLLL = False
_FPLLL = shutil.which(os.environ.get("FPLLL", "fplll"))
_JULIA = os.environ.get("JULIA", "julia")
_JL_SCRIPT = r"""
using Nemo, JSON
D = JSON.parsefile(ARGS[1])
function red(d)
  n=length(d); m=length(d[1]); M=zero_matrix(ZZ,n,m)
  for i in 1:n, j in 1:m; M[i,j]=ZZ(d[i][j]); end
  L=lll(M); [[string(L[i,j]) for j in 1:m] for i in 1:n]
end
out = isa(D[1][1], AbstractArray) ? [red(d) for d in D] : red(D)
open(ARGS[2],"w") do f; JSON.print(f, out); end
"""

ENDED_BY = "ended by MPLLL_FPLLL_TIMEOUT=%s s, not by convergence"


class MplllTimeoutExpired(subprocess.TimeoutExpired):
    """Raised at the backend boundary when the MPLLL_FPLLL_TIMEOUT knob ends a
    reduction.  A subprocess.TimeoutExpired (callers that catch that class keep
    catching it); str() is the ended-by line so any receipt that prints the
    exception carries it."""

    def __init__(self, cmd, timeout, backend, what, output=None, stderr=None):
        super().__init__(cmd, timeout, output=output, stderr=stderr)
        self.backend = backend
        self.what = what
        self.ended_by = (ENDED_BY % timeout) + " [backend=%s, %s]" % (backend, what)

    def __str__(self):
        return self.ended_by


def _fplll_timeout():
    """The per-reduction wall from env MPLLL_FPLLL_TIMEOUT.
    unset or "0" -> None (unbounded, the default); a positive integer (ASCII
    digits) -> that many seconds (an opt-in the caller states in its receipt);
    anything else -> ValueError naming the knob."""
    raw = os.environ.get("MPLLL_FPLLL_TIMEOUT")
    if raw is None:
        return None
    s = raw.strip()
    if not (s.isascii() and s.isdigit()):
        raise ValueError("MPLLL_FPLLL_TIMEOUT=%r: expected unset, 0 (unbounded) "
                         "or a positive integer number of seconds" % raw)
    n = int(s)
    return n if n > 0 else None


def _ended_by(exc, timeout, backend, what):
    """Build the MplllTimeoutExpired for a subprocess.TimeoutExpired caught at a
    backend boundary and write its line to stderr; the caller raises it."""
    err = MplllTimeoutExpired(exc.cmd, timeout, backend, what,
                              output=exc.output, stderr=exc.stderr)
    sys.stderr.write("[mplll_lattice] %s\n" % err)
    return err


def _lll_fpylll(rows):
    A = IntegerMatrix(len(rows), len(rows[0]))
    for i, r in enumerate(rows):
        for j, x in enumerate(r):
            A[i, j] = int(x)
    _fpylll_LLL.reduction(A)
    return [[int(A[i, j]) for j in range(A.ncols)] for i in range(A.nrows)]


def _bkz_fpylll(rows, beta):
    from fpylll import BKZ
    A = IntegerMatrix(len(rows), len(rows[0]))
    for i, r in enumerate(rows):
        for j, x in enumerate(r):
            A[i, j] = int(x)
    BKZ.reduction(A, BKZ.Param(block_size=beta))
    return [[int(A[i, j]) for j in range(A.ncols)] for i in range(A.nrows)]


def _bkz_fplll_cli(rows, beta, timeout=None):
    """Shell to `fplll -a bkz -b β -f mpfr -p <bits>`. Big-int safe (GMP).
    timeout: seconds; None (default) = the MPLLL_FPLLL_TIMEOUT knob via
    _fplll_timeout() (unset/0 = unbounded)."""
    if timeout is None:
        timeout = _fplll_timeout()
    inp = "[" + "".join("[" + " ".join(str(int(x)) for x in r) + "]" for r in rows) + "]"
    bits = max(64, max(int(x).bit_length() for r in rows for x in r) + 64)
    args = [_FPLLL] + (["-a", "bkz", "-b", str(beta)] if beta > 0
                       else ["-a", "lll", "-m", "proved"]) \
        + ["-f", "mpfr", "-p", str(bits)]
    try:
        p = subprocess.run(args, input=inp, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired as e:
        raise _ended_by(e, timeout, "fplll-cli",
                        "lattice %dx%d" % (len(rows), len(rows[0]))) from e
    if p.returncode != 0:
        raise RuntimeError(f"fplll cli failed: {p.stderr[-400:]}")
    out, m = [], len(rows[0])
    for ln in p.stdout.splitlines():
        toks = re.findall(r'-?\d+', ln)
        if len(toks) == m:
            out.append([int(t) for t in toks])
    if len(out) != len(rows):
        raise RuntimeError(f"fplll cli parse: got {len(out)} rows, want {len(rows)}")
    return out


def _lll_nemo(rows, timeout=None, batch=False):
    """Shell to julia/Nemo for LLL. batch=True ⇒ rows is list-of-matrices.
    timeout: seconds; None (default) = the MPLLL_FPLLL_TIMEOUT knob via
    _fplll_timeout() (unset/0 = unbounded)."""
    if timeout is None:
        timeout = _fplll_timeout()
    with tempfile.TemporaryDirectory(prefix="mplll_") as td:
        inp, out, scr = (os.path.join(td, n) for n in ("in.json", "out.json", "run.jl"))
        enc = (lambda M: [[str(int(x)) for x in r] for r in M])
        json.dump([enc(m) for m in rows] if batch else enc(rows), open(inp, "w"))
        open(scr, "w").write(_JL_SCRIPT)
        try:
            r = subprocess.run([_JULIA, "--startup-file=no", scr, inp, out],
                               capture_output=True, text=True, timeout=timeout)
        except subprocess.TimeoutExpired as e:
            raise _ended_by(e, timeout, "nemo",
                            ("batch of %d" % len(rows)) if batch
                            else "lattice %dx%d" % (len(rows), len(rows[0]))) from e
        if r.returncode != 0:
            raise RuntimeError(f"Nemo.lll failed: {r.stderr[-800:]}")
        L = json.load(open(out))
    dec = (lambda M: [[int(x) for x in r] for r in M])
    return [dec(m) for m in L] if batch else dec(L)


def _lll_pure(rows):
    sys.path.insert(0, os.path.join(          # g2 retarget: gatekeeper lives in
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),  # tools/, one dir
        "gatekeeper"))                        # above this package (lockpick/)
    from heldout_cv import _int_lll
    return _int_lll([list(r) for r in rows])


def reduce_rows(rows, backend=None, bkz_beta=0):
    """LLL/BKZ-reduce integer row lattice.  Returns (reduced, backend_str, wall_s).
    BKZ via fpylll → fplll-cli → nemo-LLL fallback (noted in backend string)."""
    t0 = time.time()
    if bkz_beta and _HAVE_FPYLLL:
        R, b = _bkz_fpylll(rows, bkz_beta), f"fpylll-bkz{bkz_beta}"
    elif bkz_beta and _FPLLL and backend != "nemo":
        R, b = _bkz_fplll_cli(rows, bkz_beta), f"fplll-cli-bkz{bkz_beta}"
    elif backend == "fpylll" or (backend is None and _HAVE_FPYLLL):
        R, b = _lll_fpylll(rows), "fpylll"
    elif backend == "fplll-cli" or (backend is None and _FPLLL):
        R, b = _bkz_fplll_cli(rows, 0), "fplll-cli-lll"
    elif backend == "pure":
        R, b = _lll_pure(rows), "pure"
    else:
        try:
            R, b = _lll_nemo(rows), ("nemo-lll(no-bkz)" if bkz_beta else "nemo")
        except MplllTimeoutExpired:
            raise   # the knob's ending propagates; the pure fallback is not a stand-in for it
        except Exception as e:
            sys.stderr.write(f"[mplll_lattice] nemo failed ({e}); pure fallback\n")
            R, b = _lll_pure(rows), "pure"
    return R, b, time.time() - t0


def reduce_rows_batch(mats, backend=None, bkz_beta=0):
    """Reduce a list of row-matrices. One julia subprocess on nemo path; loops
    the fast CLI/fpylll backends otherwise. Returns (list[reduced], backend, wall_s)."""
    t0 = time.time()
    if bkz_beta and _HAVE_FPYLLL:
        Rs, b = [_bkz_fpylll(m, bkz_beta) for m in mats], f"fpylll-bkz{bkz_beta}"
    elif bkz_beta and _FPLLL and backend != "nemo":
        Rs, b = [_bkz_fplll_cli(m, bkz_beta) for m in mats], f"fplll-cli-bkz{bkz_beta}"
    elif backend == "fplll-cli" or (backend is None and _FPLLL):
        Rs, b = [_bkz_fplll_cli(m, 0) for m in mats], "fplll-cli-lll"
    elif backend == "fpylll" or (backend is None and _HAVE_FPYLLL):
        Rs, b = [_lll_fpylll(m) for m in mats], "fpylll"
    else:
        try:
            Rs, b = _lll_nemo(mats, batch=True), ("nemo-lll(no-bkz)" if bkz_beta else "nemo")
        except MplllTimeoutExpired:
            raise   # the knob's ending propagates; the pure fallback is not a stand-in for it
        except Exception as e:
            sys.stderr.write(f"[mplll_lattice] nemo batch failed ({e}); pure fallback\n")
            Rs, b = [_lll_pure(m) for m in mats], "pure"
    return Rs, b, time.time() - t0


lll_rows = reduce_rows   # backward-compat alias


def row_norm2(row):
    return sum(int(x) * int(x) for x in row)


def gcd_list(xs):
    g = 0
    for x in xs:
        g = math.gcd(g, abs(int(x)))
    return g
