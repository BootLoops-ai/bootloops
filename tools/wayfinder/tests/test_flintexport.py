#!/usr/bin/env python3
"""wayfinder.flintexport leg: the module's in-built exact-identity selftest
(5 legs: exact zero after expansion, gcd cancel, 3-var (d,s2,s12) pin identity,
canonical export shape + positive-lead law, non-integer-exponent refusal),
run under the package battery.  SKIPS BY NAME without python-flint or sympy —
flintexport is the one Wayfinder member that requires python-flint; the rest
of the package treats it as an optional accelerator.

Collected by `pytest tools/wayfinder`; also runnable as a script.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))


def _missing():
    try:
        import flint  # noqa: F401
    except ImportError:
        return "python-flint"
    try:
        import sympy  # noqa: F401
    except ImportError:
        return "sympy"
    return None


def test_flintexport_selftest():
    miss = _missing()
    if miss:
        import pytest
        pytest.skip(f"{miss} not installed: flintexport leg skipped")
    from wayfinder import flintexport            # lazy member
    fails = flintexport._selftest()
    assert fails == 0, f"flintexport selftest: {fails} identity legs failed"


def test_flintexport_lazy_not_loaded_by_package_import():
    """package import must stay python-flint/sympy-free: flintexport (and
    epslimit) load only on first touch.  Checked in a fresh interpreter so the
    battery's already-imported modules are left alone."""
    import subprocess
    tools_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
    code = ("import sys; sys.path.insert(0, %r); import wayfinder; "
            "bad = [m for m in ('wayfinder.flintexport', 'wayfinder.epslimit', 'sympy') "
            "if m in sys.modules]; print('eager:', bad); sys.exit(1 if bad else 0); "
            % tools_dir)
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert r.returncode == 0, f"package import pulled lazy members: {r.stdout}{r.stderr}"
    # and the lazy hook resolves both names
    code2 = ("import sys; sys.path.insert(0, %r); import wayfinder; "
             "from wayfinder import epslimit; assert hasattr(epslimit, 'extrapolate')" % tools_dir)
    r2 = subprocess.run([sys.executable, "-c", code2], capture_output=True, text=True)
    assert r2.returncode == 0, f"lazy epslimit hook failed: {r2.stderr}"


if __name__ == '__main__':
    miss = _missing()
    if miss:
        print(f"flintexport leg: SKIP ({miss} not installed)")
        sys.exit(0)
    test_flintexport_lazy_not_loaded_by_package_import()
    from wayfinder import flintexport
    sys.exit(flintexport._selftest())
