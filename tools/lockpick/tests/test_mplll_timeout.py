"""Regression: the lockpick lattice backends' per-reduction wall is UNBOUNDED by
default; MPLLL_FPLLL_TIMEOUT is an opt-in, and a reduction the knob ends says
so by name.

Why: a finite default timeout would end a long proved-mode fplll reduction by
the clock, and running compute must never be killed by a clock it did not opt
into.  The default is 0
(unbounded), mapped to subprocess timeout=None (the unbounded sentinel of
subprocess.run; a huge integer is not used as a stand-in), the knob a
deliberate opt-in the caller states in its receipt, and a reduction ended by
the knob raises MplllTimeoutExpired (a subprocess.TimeoutExpired) whose text
reads 'ended by MPLLL_FPLLL_TIMEOUT=<N> s, not by convergence'.

Legs: the knob parser (unset / "0" -> None; "7" -> 7; "x" -> ValueError by
name); a PLANTED control with subprocess.run monkeypatched to raise
TimeoutExpired on both subprocess backends (the message text); a REAL timeout
against a sleeping stand-in fplll binary (the child is ended by the knob and
the message surfaces through reduce_rows); a real small LLL on the fplll-cli
backend under the default, with a spy asserting the subprocess.run call carried
timeout=None (skipped by name when fplll is not installed); a source-level
drift guard that every subprocess.run site in mplll_lattice.py passes the
knob's value.

Also: a nemo reduction the knob ends is not
absorbed by the nemo path's pure-LLL fallback in reduce_rows /
reduce_rows_batch (it propagates as MplllTimeoutExpired); the knob's value is
ASCII digits (a Unicode digit-like such as the superscript two, or fullwidth
digits, is refused by name).
"""
import inspect
import os
import stat
import subprocess
import time

import pytest

from lockpick import mplll_lattice as ML

KNOB = "MPLLL_FPLLL_TIMEOUT"
PHRASE = "not by convergence"
ROWS = [[1, 0, 0, 7], [0, 1, 0, 11], [0, 0, 1, 13]]


@pytest.mark.parametrize("value", [None, "0", " 0 ", "00"])
def test_unset_or_zero_is_unbounded(monkeypatch, value):
    if value is None:
        monkeypatch.delenv(KNOB, raising=False)
    else:
        monkeypatch.setenv(KNOB, value)
    assert ML._fplll_timeout() is None


@pytest.mark.parametrize("value,want", [("7", 7), ("1800", 1800), (" 3600 ", 3600)])
def test_positive_integer_is_seconds(monkeypatch, value, want):
    monkeypatch.setenv(KNOB, value)
    assert ML._fplll_timeout() == want


@pytest.mark.parametrize("value", ["x", "", "-5", "1.5", "1e9", "+7", "None", "\u00b2"])
def test_anything_else_is_refused_by_name(monkeypatch, value):
    monkeypatch.setenv(KNOB, value)
    with pytest.raises(ValueError) as ei:
        ML._fplll_timeout()
    assert KNOB in str(ei.value)


def _planted_run(*args, **kwargs):
    raise subprocess.TimeoutExpired(cmd=args[0], timeout=kwargs.get("timeout"))


def test_planted_timeout_fplll_cli_carries_the_ended_by_line(monkeypatch):
    monkeypatch.setenv(KNOB, "7")
    monkeypatch.setattr(ML.subprocess, "run", _planted_run)
    monkeypatch.setattr(ML, "_FPLLL", "/nonexistent/fplll")
    with pytest.raises(subprocess.TimeoutExpired) as ei:
        ML._bkz_fplll_cli(ROWS, 0)
    assert isinstance(ei.value, ML.MplllTimeoutExpired)
    assert str(ei.value).startswith("ended by MPLLL_FPLLL_TIMEOUT=7 s, " + PHRASE)
    assert "fplll-cli" in str(ei.value)


def test_planted_timeout_nemo_carries_the_ended_by_line(monkeypatch):
    monkeypatch.setenv(KNOB, "9")
    monkeypatch.setattr(ML.subprocess, "run", _planted_run)
    with pytest.raises(subprocess.TimeoutExpired) as ei:
        ML._lll_nemo(ROWS)
    assert str(ei.value).startswith("ended by MPLLL_FPLLL_TIMEOUT=9 s, " + PHRASE)
    assert "nemo" in str(ei.value)


def test_planted_timeout_stderr_line(monkeypatch, capsys):
    monkeypatch.setenv(KNOB, "7")
    monkeypatch.setattr(ML.subprocess, "run", _planted_run)
    with pytest.raises(subprocess.TimeoutExpired):
        ML._bkz_fplll_cli(ROWS, 0)
    err = capsys.readouterr().err
    assert "[mplll_lattice] ended by MPLLL_FPLLL_TIMEOUT=7 s, " + PHRASE in err


def test_real_timeout_against_a_sleeping_stand_in(monkeypatch, tmp_path):
    """The knob really ends the child: a stand-in 'fplll' that sleeps 30 s is
    ended after 1 s and the message surfaces through reduce_rows."""
    fake = tmp_path / "fplll"
    fake.write_text("#!/bin/sh\nsleep 30\n")
    fake.chmod(fake.stat().st_mode | stat.S_IXUSR)
    monkeypatch.setenv(KNOB, "1")
    monkeypatch.setattr(ML, "_FPLLL", str(fake))
    monkeypatch.setattr(ML, "_HAVE_FPYLLL", False)
    with pytest.raises(subprocess.TimeoutExpired) as ei:
        ML.reduce_rows(ROWS, backend="fplll-cli")
    assert str(ei.value).startswith("ended by MPLLL_FPLLL_TIMEOUT=1 s, " + PHRASE)


def test_nemo_ending_is_not_swallowed_by_the_pure_fallback(monkeypatch, tmp_path):
    """A Nemo.lll reduction the knob ends propagates out of reduce_rows
    and reduce_rows_batch as MplllTimeoutExpired; the nemo path's `except
    Exception` -> pure-LLL fallback does not absorb it (the failure mode would
    be a return with backend 'pure', the ended-by line on stderr only)."""
    fake = tmp_path / "julia"
    fake.write_text("#!/bin/sh\nexec sleep 30\n")
    fake.chmod(fake.stat().st_mode | stat.S_IXUSR)
    monkeypatch.setenv(KNOB, "1")
    monkeypatch.setattr(ML, "_JULIA", str(fake))
    monkeypatch.setattr(ML, "_HAVE_FPYLLL", False)
    for name, call in (("reduce_rows", lambda: ML.reduce_rows(ROWS, backend="nemo")),
                       ("reduce_rows_batch", lambda: ML.reduce_rows_batch([ROWS, ROWS], backend="nemo"))):
        t0 = time.time()
        try:
            _, backend, _ = call()
        except ML.MplllTimeoutExpired as e:
            assert str(e).startswith("ended by MPLLL_FPLLL_TIMEOUT=1 s, " + PHRASE)
            assert "nemo" in str(e)
            assert time.time() - t0 < 20, "ended by the knob, not by the stand-in's 30 s sleep"
            continue
        pytest.fail("%s swallowed the nemo ending and returned backend %r "
                    "instead of raising MplllTimeoutExpired" % (name, backend))


def test_real_small_lll_under_the_default_is_unbounded(monkeypatch):
    if not ML._FPLLL:
        pytest.skip("fplll not installed (shutil.which found no fplll binary)")
    monkeypatch.delenv(KNOB, raising=False)
    seen = []
    real_run = subprocess.run

    def spy(*args, **kwargs):
        seen.append(kwargs.get("timeout", "ABSENT"))
        return real_run(*args, **kwargs)

    monkeypatch.setattr(ML.subprocess, "run", spy)
    R, backend, wall = ML.reduce_rows(ROWS, backend="fplll-cli")
    assert backend == "fplll-cli-lll"
    assert len(R) == len(ROWS) and all(len(r) == len(ROWS[0]) for r in R)
    assert seen == [None], "the fplll-cli subprocess.run call must carry timeout=None under the default"


def test_every_subprocess_run_site_reads_the_knob():
    src = inspect.getsource(ML)
    assert src.count("subprocess.run(") == 2
    assert src.count("timeout=timeout") == 2
    assert src.count("timeout = _fplll_timeout()") == 2
    assert src.count("raise _ended_by(e, timeout,") == 2
    assert 'os.environ.get("MPLLL_FPLLL_TIMEOUT", "1800")' not in src
    sig = inspect.signature(ML._lll_nemo)
    assert sig.parameters["timeout"].default is None
    sig = inspect.signature(ML._bkz_fplll_cli)
    assert sig.parameters["timeout"].default is None
